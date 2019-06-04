/*****************************************************************************
* ZrAux1.c - Aux tools to computes intersection curve of two surfaces.	     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Michael Barton and Gershon Elber		Ver 1.0, April 2008  *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"
#include "mvar_loc.h"

#define MVAR_ZERO_RECOVER_TRACE_ERR

#define MVAR_MVZR1D_EPS_WED_PROD	    1e-6
#define MVAR_MVZR1D_END_TOL_FACTOR	    2

static MvarPolylineStruct *MvarMVZR1DMergeTwoPoly(MvarPolylineStruct *Poly1, 
					          MvarPolylineStruct *Poly2, 
					          CagdRType Tol);
static int MvarMVZR1DPolyPtOnBound(MvarPolylineStruct *Poly,
				   CagdRType SubdivTol, 
				   CagdRType NumericTol, 
				   int BoundarySide, 
				   CagdRType BoundaryValue);
static int MvarMVZR1DPtCmpPoints(const MvarPtStruct *P1,
				 const MvarPtStruct *P2,
				 CagdRType SubdivTol, 
				 CagdRType NumericTol,
				 CagdRType *UsedTol);
static MvarPolylineStruct *MvarMVZR1DPrelimLink(MvarPolylineStruct **PolyList, 
					        MvarPolylineStruct *Poly,
					        CagdRType SubdivTol, 
					        CagdRType NumericTol,
					        CagdRType *UsedTol);
static MvarMVStruct **MvarMVZR1DCreateMVs(const CagdSrfStruct *Srf1, 
					  const CagdSrfStruct *Srf2);
static int MvarMVZR1DLinDep(MvarVecStruct **Vectors,
			    int Size, 
			    MvarMVZR1DAuxStruct *AS);
static MvarVecStruct *MvarMVZR1DWedgeProd(MvarVecStruct **Vectors,
					  int Size, 
					  MvarMVZR1DAuxStruct *AS);
static MvarPtStruct *MvarMVZR1DStepInTangDir(MvarMVStruct * const *MVs,  
					     MvarPtStruct *StartPoint,
					     MvarVecStruct *DirVec,
					     CagdRType Step,
					     MvarMVGradientStruct **MVGrads,
					     MvarMVZR1DAuxStruct *AS);
static MvarPtStruct *MvarMVZR1DCorrectionStep(MvarMVStruct * const *MVs,  
					      MvarPtStruct *StartPoint, 
					      MvarMVGradientStruct **MVGrads,
					      MvarMVZR1DAuxStruct *AS);
static int MvarMVZR1DCloseToIntersCrv(MvarMVStruct * const *MVs,
				      int Size, 
				      const MvarPtStruct *Point,
				      CagdRType Tol);
static CagdRType MvarMVZR1DEvalErr(MvarMVStruct * const *MVs,
				   int Size, 
				   const MvarPtStruct *Pt);

#ifdef DEBUG_MVAR_MVZR1D_LINK_NEIGH
static void MvarDbgMVZR1DPrintEndPtPlList(MvarPolylineStruct *PolyList);
static void MvarDbgMVZR1DExamineSegmentBndry(MvarPtStruct *PtList,
					     const MvarMVStruct *MV);
#endif /* DEBUG_MVAR_MVZR1D_LINK_NEIGH */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Merges two neighboring polylines together, in place.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly1:      1st polyline.						     *
*   Poly2:	2nd polyline.						     *
*   Tol:	The tolerance under which two points are proclaimed to be    *
*		identical.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *:  If the start/endpoint of Poly1 coincides with the *
*		        start/endpoint of Poly2, up to Tol, they are merged. *
*****************************************************************************/
static MvarPolylineStruct *MvarMVZR1DMergeTwoPoly(MvarPolylineStruct *Poly1, 
					          MvarPolylineStruct *Poly2, 
					          CagdRType Tol)
{
    MvarPtStruct *Pt1Last, *Pt2Last;

    if (!MvarPtCmpTwoPoints(Poly1 -> Pl, Poly2 -> Pl, Tol)) {
	Pt1Last = Poly1 -> Pl;
	Poly1 -> Pl = CagdListReverse(Poly1 -> Pl);
    }
    else if (!MvarPtCmpTwoPoints(Pt1Last = CagdListLast(Poly1 -> Pl),
				 Poly2 -> Pl, Tol)) {
    }
    else if (!MvarPtCmpTwoPoints(Poly1 -> Pl,
				 Pt2Last = CagdListLast(Poly2 -> Pl), Tol)) {
        IRIT_SWAP(MvarPolylineStruct *, Poly1, Poly2);
	Pt1Last = Pt2Last;
    } 
    else if (!MvarPtCmpTwoPoints(Pt1Last, Pt2Last, Tol)) { 
	Poly2 -> Pl = CagdListReverse(Poly2 -> Pl);
    }
    else {
#	ifdef DEBUG
	    fprintf(stderr, "MvarZeroMV1D: Polylines can not be merged.\n");
#	endif /* DEBUG */	

	return NULL;
    }

    if (MvarPtCmpTwoPoints(Pt1Last, Poly2 -> Pl, IRIT_EPS) == 0) {
	Poly1 -> Pl = CagdListAppend(Poly1 -> Pl, Poly2 -> Pl -> Pnext);
	Poly2 -> Pl -> Pnext = NULL;
    }
    else {
        Poly1 -> Pl = CagdListAppend(Poly1 -> Pl, Poly2 -> Pl);
	Poly2 -> Pl = NULL;
    }

    MvarPolylineFree(Poly2);

    return Poly1;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Classifies the mutual position of polyline and the boundary (of	     *
* its domain).								     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly:	    A polyline to classify.				     *
*   SubdivTol, NumericTol: The tolerances under which two points are         *
*                   proclaimed to be identical.				     *
*   BoundarySide:   The side of the domain.				     *
*   BoundaryValue:  The value at the boundary.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:	    0	neither startpoint nor endpoint of polyline lies     *
*			on the boundary.				     *
*		    1	startpoint lies on the boundary, endpoint not.	     *
*		    2	endpoint lies on the boundary, startpoint not.	     *
*		    3	startpoint and endpoint both lie on the boundary.    *
*****************************************************************************/
static int MvarMVZR1DPolyPtOnBound(MvarPolylineStruct *Poly,
				   CagdRType SubdivTol, 
				   CagdRType NumericTol, 
				   int BoundarySide,
				   CagdRType BoundaryValue)
{  
    MvarPtStruct
        *TempPt = CagdListLast(Poly -> Pl);
    int MidPt1 = AttrGetIntAttrib(Poly -> Pl -> Attr, "_MidMVZR1DPt") == TRUE,
        MidPt2 = AttrGetIntAttrib(TempPt -> Attr, "_MidMVZR1DPt") == TRUE;
    CagdRType
        Tol1 = MidPt1 ? SubdivTol : NumericTol,
        Tol2 = MidPt2 ? SubdivTol : NumericTol;

    return (IRIT_APX_EQ_EPS(Poly -> Pl -> Pt[BoundarySide],
			    BoundaryValue, Tol1) ? 1 : 0) +
           (IRIT_APX_EQ_EPS(TempPt -> Pt[BoundarySide],
			    BoundaryValue, Tol2) ? 2 : 0);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   A comparison function to examine if the given two points are the same.   *
* This special MVZR1D version is due to single points that we might merge at *
* subdivision tolerance if we stop at subdiv. tol.                           *
*									     *
* PARAMETERS:                                                                *
*   P1, P2:         Two multivariate points to compare.	                     *
*   SubdivTol, NumericTol:   The tolerance of the comparison.                *
*   UsedTol:        The actual tolerance used in the comparison (Single mid  *
*                   points use SubdivTol while normally NumericTol is used). *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:      0 if identical, -1 or +1 if first point is less than/greater   *
*	      than second point, in lexicographic order over dimensions.     *
*****************************************************************************/
static int MvarMVZR1DPtCmpPoints(const MvarPtStruct *P1,
				 const MvarPtStruct *P2,
				 CagdRType SubdivTol, 
				 CagdRType NumericTol,
				 CagdRType *UsedTol)
{
    int i,
        Dim = P1 -> Dim,
        MidPt1 = AttrGetIntAttrib(P1 -> Attr, "_MidMVZR1DPt") == TRUE,
        MidPt2 = AttrGetIntAttrib(P2 -> Attr, "_MidMVZR1DPt") == TRUE;
    CagdRType
        Tol = MidPt1 || MidPt2 ? SubdivTol : NumericTol;

    *UsedTol = Tol;

    if (Dim != P2 -> Dim)
	return FALSE;

    for (i = 0; i < Dim; i++) {
	if (!IRIT_APX_EQ_EPS(P1 -> Pt[i], P2 -> Pt[i], Tol))
	    return IRIT_SIGN(P1 -> Pt[i] - P2 -> Pt[i]);
    }

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Finds a polyline in the list, which is to be merged to a given polyline  *
* and remove it from the input list.					     *
*                                                                            *
* PARAMETERS:                                                                *
*   PolyList:	   List of polylines (found poly is removed from, in place). *
*   Poly:	   Polyline, to be merged with the one from the list.	     *
*   SubdivTol, NumericTol: The tolerances under which two points are         *
*		   proclaimed to be identical.				     *
*   UsedTol:       The actual tolerance used in the comparison (Single mid   *
*                  points use SubdivTol while normally NumericTol is used).  *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *:  Polyline which is suitable for Poly to be merged, *
*			NULL, if no such polyline exist.		     *
*****************************************************************************/
static MvarPolylineStruct *MvarMVZR1DPrelimLink(MvarPolylineStruct **PolyList, 
						MvarPolylineStruct *Poly,
						CagdRType SubdivTol, 
						CagdRType NumericTol,
						CagdRType *UsedTol)
{
    MvarPolylineStruct
        *PrevPoly = NULL,
        *TempPoly = *PolyList;
    MvarPtStruct
        *PtLast = CagdListLast(Poly -> Pl);

    while (TempPoly != NULL) {
        MvarPtStruct
	    *PtTempLast = CagdListLast(TempPoly -> Pl);

	if (!MvarMVZR1DPtCmpPoints(TempPoly -> Pl, Poly -> Pl,
				   SubdivTol, NumericTol, UsedTol) ||
	    !MvarMVZR1DPtCmpPoints(TempPoly -> Pl, PtLast,
				   SubdivTol, NumericTol, UsedTol) ||
	    !MvarMVZR1DPtCmpPoints(PtTempLast, Poly -> Pl,
				   SubdivTol, NumericTol, UsedTol) ||
	    !MvarMVZR1DPtCmpPoints(PtTempLast, PtLast,
				   SubdivTol, NumericTol, UsedTol)) {
	    if (PrevPoly == NULL)
	        *PolyList = (*PolyList) -> Pnext;
	    else
	        PrevPoly -> Pnext = TempPoly -> Pnext;

	    return TempPoly;
	}
	PrevPoly = TempPoly;
	TempPoly = TempPoly -> Pnext;
    }

    return NULL;
}

#ifdef DEBUG_MVAR_MVZR1D_LINK_NEIGH

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Print the end points of a list of polylines.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   PolyList:	    List of polylines.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void		                                                     *
*****************************************************************************/
static void MvarDbgMVZR1DPrintEndPtPlList(MvarPolylineStruct *PolyList)
{
    MvarPolylineStruct *Pl;

    for (Pl = PolyList; Pl != NULL; Pl = Pl -> Pnext) {
        int i;
        MvarPtStruct
	    *Pt1 = Pl -> Pl,
	    *Pt2 = CagdListLast(Pt1);

	printf("Plln from:");
	for (i = 0; i < Pt1 -> Dim; i++)
	    printf(" %12.10f", Pt1 -> Pt[i]);
	printf(" to:");
	for (i = 0; i < Pt2 -> Dim; i++)
	    printf(" %12.10f", Pt2 -> Pt[i]);
	printf("\n");
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Verify the given segment is closed or on MV's boundary.                  *
*                                                                            *
* PARAMETERS:                                                                *
*   PtList:	    One traced segment to examine...			     *
*   MV:             If PtList starting/ending on MV boundary or closed.      *
*                                                                            *
* RETURN VALUE:                                                              *
*   void		                                                     *
*****************************************************************************/
static void MvarDbgMVZR1DExamineSegmentBndry(MvarPtStruct *PtList,
					     const MvarMVStruct *MV)
{
    int i, j;
    CagdRType t1, t2;
    MvarPtStruct
        *LastPt = CagdListLast(PtList);

    for (i = 0; i < MV -> Dim; i++) {
        MvarMVDomain(MV, &t1, &t2, i);
	if (IRIT_APX_EQ(PtList -> Pt[i], t1) ||
	    IRIT_APX_EQ(PtList -> Pt[i], t2))
	    break; /* Starting point is on boundary. */
    }

    for (j = 0; j < MV -> Dim; j++) {
        MvarMVDomain(MV, &t1, &t2, j);
	if (IRIT_APX_EQ(LastPt -> Pt[j], t1) ||
	    IRIT_APX_EQ(LastPt -> Pt[j], t2))
	    break; /* End point is on boundary. */
    }

    if (i >= MV -> Dim || j >= MV -> Dim) {
        fprintf(stderr,
		"Incomplete segment created (not on boundary), %d points.\n", 
		CagdListLength(PtList));

	fprintf(stderr, "Plln from:");
	for (i = 0; i < PtList -> Dim; i++)
	    fprintf(stderr, " %6.3f", PtList -> Pt[i]);
	fprintf(stderr, " to:");
	for (i = 0; i < LastPt -> Dim; i++)
	    fprintf(stderr, " %6.3f", LastPt -> Pt[i]);

	fprintf(stderr, "\nMV Domain:");
	for (i = 0; i < MV -> Dim; i++) {
	    MvarMVDomain(MV, &t1, &t2, i);
	    fprintf(stderr, " [%6.3f, %6.3f]", t1, t2);
	}
	fprintf(stderr, "\n");
    }
}

#endif /* DEBUG_MVAR_MVZR1D_LINK_NEIGH */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Merges two neighboring lists of polylines that approximate a univariate  *
* intersection curve. ("inversion" operation to subdivision).		     *
*                                                                            *
* PARAMETERS:                                                                *
*   PolyList1:	    1st list of polylines.				     *
*   PolyList2:	    2nd list of polylines.				     *
*   SubdivTol, NumericTol: The tolerances under which two points are         *
*                   proclaimed to be identical.				     *
*   BoundarySide:   The side of the domain along which the lists are merged. *
*   BoundaryValue:  The value at the boundary.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *:   A list of polylines, merged PolyList1 and        *
*                           PolyList2.					     *
*****************************************************************************/
MvarPolylineStruct *MvarMVZR1DLinkNeighbours(MvarPolylineStruct *PolyList1, 
  					     MvarPolylineStruct *PolyList2, 
					     CagdRType SubdivTol, 
					     CagdRType NumericTol,
					     int BoundarySide, 
					     CagdRType BoundaryValue)
{
    MvarPolylineStruct *Pl, *Polys,
	*PlsTmp = NULL,
        *PolysOut = NULL;

#ifdef DEBUG_MVAR_MVZR1D_LINK_NEIGH
    if (PolyList1 != NULL && PolyList2 != NULL) {
        printf("\n\nFirst list (Side = %d, Param = %f):\n",
		BoundarySide, BoundaryValue);
	MvarDbgMVZR1DPrintEndPtPlList(PolyList1);
	printf("\nSecond list:\n");
	MvarDbgMVZR1DPrintEndPtPlList(PolyList2);
    }
#endif /* DEBUG_MVAR_MVZR1D_LINK_NEIGH */

    Polys = CagdListAppend(PolyList1, PolyList2);
    if (PolyList1 == NULL || PolyList2 == NULL)
        return Polys;

    /* Move to out-list all polys that are not on the share boundary. */
    while (Polys != NULL) {
        IRIT_LIST_POP(Pl, Polys);

	if (MvarMVZR1DPolyPtOnBound(Pl, SubdivTol, NumericTol,
				    BoundarySide, BoundaryValue) == 0) {
	    IRIT_LIST_PUSH(Pl, PolysOut);
	}
	else {
	    IRIT_LIST_PUSH(Pl, PlsTmp);
	}
    }
    Polys = PlsTmp;

    /* Start the merge process. */
    while (Polys != NULL) {
        CagdRType UsedTol;
	MvarPolylineStruct *Pl2;

        IRIT_LIST_POP(Pl, Polys);

	if ((Pl2 = MvarMVZR1DPrelimLink(&Polys, Pl, SubdivTol, NumericTol,
					&UsedTol)) == NULL) {
	    IRIT_LIST_PUSH(Pl, PolysOut);
	}
	else {
	    Pl = MvarMVZR1DMergeTwoPoly(Pl, Pl2, UsedTol);
	    IRIT_LIST_PUSH(Pl, Polys);
	}
    }

#ifdef DEBUG_MVAR_MVZR1D_LINK_NEIGH
    if (PolyList1 != NULL && PolyList2 != NULL) {
        printf("\nMerged list:\n");
	MvarDbgMVZR1DPrintEndPtPlList(PolysOut);
    }
#endif /* DEBUG_MVAR_MVZR1D_LINK_NEIGH */

    return PolysOut;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Allocates the memory required for a new MvarMVZR1DAuxStruct.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   NumOfZeroMVs:    Dim - 1 = number of zero constraints that describe the  *
*		 univariate-solution problem (4 equations with 3 unknowns,   *
*                for a regular SSI).					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVZR1DAuxStruct *:    An uninitialized auxiliary structure.          *
*****************************************************************************/
MvarMVZR1DAuxStruct *MvarMVZR1DAllocOnce(int NumOfZeroMVs)
{
    int i,
        ZNum1 = NumOfZeroMVs + 1;
    MvarMVZR1DAuxStruct
	*AuxStruct = (MvarMVZR1DAuxStruct *) IritMalloc(sizeof(
							 MvarMVZR1DAuxStruct));

    AuxStruct -> NumOfMVs = NumOfZeroMVs;

    AuxStruct -> OrthoBasis =
        (MvarVecStruct **) IritMalloc(sizeof(MvarVecStruct *) * ZNum1); 

    for (i = 0; i < ZNum1; i++)
        AuxStruct -> OrthoBasis[i] = MvarVecNew(ZNum1);

    AuxStruct -> TempVec = MvarVecNew(ZNum1);
    AuxStruct -> CorrVec = MvarVecNew(ZNum1);

    AuxStruct -> MinDmn = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));
    AuxStruct -> MaxDmn = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));

    AuxStruct -> MinDmn2 = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));
    AuxStruct -> MaxDmn2 = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));

    AuxStruct -> GradVecs =
	(MvarVecStruct **) IritMalloc(sizeof(MvarVecStruct *) * NumOfZeroMVs);

    for (i = 0; i < NumOfZeroMVs; i++)
        AuxStruct -> GradVecs[i] = MvarVecNew(ZNum1);

    AuxStruct -> SITTempVec = MvarVecNew(ZNum1);
    AuxStruct -> SITTanVec = MvarVecNew(ZNum1);

    AuxStruct -> A = (CagdRType *) IritMalloc(IRIT_SQR(ZNum1) * sizeof(CagdRType));
    AuxStruct -> x = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));
    AuxStruct -> b = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));
    AuxStruct -> bCopy = (CagdRType *) IritMalloc(ZNum1 * sizeof(CagdRType));

    AuxStruct -> TempList =
        (MvarVecStruct **) IritMalloc(sizeof(MvarVecStruct *) * ZNum1);

    for (i = 0; i < NumOfZeroMVs; i++)
        AuxStruct -> TempList[i] = MvarVecNew(ZNum1);	    	    

    AuxStruct -> Constraints = (MvarConstraintType *)
			 IritMalloc(sizeof(MvarConstraintType) * NumOfZeroMVs);

    for (i = 0; i < NumOfZeroMVs; i++)
        AuxStruct -> Constraints[i] = MVAR_CNSTRNT_ZERO;

    return AuxStruct;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Deallocates the memory required for a MvarMVZR1DAuxStruct.	    	     *
*                                                                            *
* PARAMETERS:                                                                *
*   AuxStruct:	Auxiliary structure to free.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void							             *
*****************************************************************************/
void MvarMVZR1DDeallocOnce(MvarMVZR1DAuxStruct *AuxStruct)
{
    int i;
    
    for (i = 0; i < AuxStruct -> NumOfMVs + 1; i++)
	MvarVecFree(AuxStruct -> OrthoBasis[i]);

    IritFree(AuxStruct -> OrthoBasis);

    MvarVecFree(AuxStruct -> TempVec);
    MvarVecFree(AuxStruct -> CorrVec);

    IritFree(AuxStruct -> MinDmn);
    IritFree(AuxStruct -> MaxDmn);
    IritFree(AuxStruct -> MinDmn2);
    IritFree(AuxStruct -> MaxDmn2);
 
    for (i = 0; i < AuxStruct -> NumOfMVs; i++)
	MvarVecFree(AuxStruct -> GradVecs[i]);	    	    

    IritFree(AuxStruct -> GradVecs);

    MvarVecFree(AuxStruct -> SITTempVec);
    MvarVecFree(AuxStruct -> SITTanVec);

    IritFree(AuxStruct -> A);
    IritFree(AuxStruct -> b);
    IritFree(AuxStruct -> bCopy);
    IritFree(AuxStruct -> x);

    for (i = 0; i < AuxStruct -> NumOfMVs; i++)
	MvarVecFree(AuxStruct -> TempList[i]);

    IritFree(AuxStruct -> TempList);

    IritFree(AuxStruct -> Constraints);

    IritFree(AuxStruct);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Eliminates similar points stored in the PtList, in place.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   PtList:	List of points, to use in place.			     *
*   Tol:	The tolerance under which two points are proclaimed to be    *
*		identical						     *
*									     *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:	New list with mutually different points.	     *
*****************************************************************************/
MvarPtStruct *MvarMVZR1DListOfDifferentPts(MvarPtStruct *PtList, 
					   CagdRType Tol)
{
    MvarPtStruct *p, *q, *LastQ, *TempPt;

    for (p = PtList; p != NULL && p -> Pnext != NULL; p = p -> Pnext) {
        for (LastQ = p, q = LastQ -> Pnext; q != NULL; ) {
	    if (MvarPtCmpTwoPoints(p, q,
				   Tol * MVAR_MVZR1D_END_TOL_FACTOR) == 0) {
	        TempPt = q;
		LastQ -> Pnext = q = TempPt -> Pnext;
		MvarPtFree(TempPt);
	    }
	    else {
		LastQ = q;
	        q = q -> Pnext;
	    }
	}
    }

    return PtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Checks the linearly dependency of vectors.	      			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Vectors:	Set of N vectors in N-dim space.			     *
*   Size:	The size N of the set.					     *
*   AS:	        Auxiliary MVZR1D structure that holds auxiliary data.	     *
*									     *
* RETURN VALUE:                                                              *
*   int:	TRUE if the vectors are linearly dependent, FALSE otherwise. *
*****************************************************************************/
static int MvarMVZR1DLinDep(MvarVecStruct **Vectors,
			    int Size, 
			    MvarMVZR1DAuxStruct *AS)
{
    int i, j; 
    CagdRType c;

    MVAR_VEC_COPY(AS -> OrthoBasis[0], Vectors[0]);

    for (i = 1; i < Size; i++) {
	MVAR_VEC_RESET(AS -> CorrVec);
	for (j = 0; j < i; j++) { 
	    MVAR_VEC_COPY(AS -> TempVec, AS -> OrthoBasis[j]);
	    c = MvarVecDotProd(AS -> OrthoBasis[j], Vectors[i]) /
		MvarVecDotProd(AS -> OrthoBasis[j], AS -> OrthoBasis[j]);
	    MvarVecScale(AS -> TempVec, -c);
	    MvarVecAdd(AS -> CorrVec, AS -> CorrVec, AS -> TempVec);
	}
	MvarVecAdd(AS -> OrthoBasis[i], Vectors[i], AS -> CorrVec);

	if (MvarVecLength(AS -> OrthoBasis[i]) < MVAR_MVZR1D_EPS_WED_PROD)
	    return TRUE;
    }
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes an orthogonal complement (wedge product) to the set of	     *
* Size vectors in higher dimensional linear space.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Vectors:   The set of Size vectors, each of dimension larger than Size.  *
*   Size:      The size of Vectors set.					     *
*   AS:	       Auxiliary MVZR1D structure that holds auxiliary data.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Wedge product.			                     *
*****************************************************************************/
static MvarVecStruct *MvarMVZR1DWedgeProd(MvarVecStruct **Vectors,
					  int Size, 
					  MvarMVZR1DAuxStruct *AS) 
{
    int i;
    MvarVecStruct *RandomVec,
        **OrthoBasis = AS -> OrthoBasis;

    assert(Size + 1 ==  Vectors[0] -> Dim);

    for (i = 0; i < Size; i++)
        MVAR_VEC_COPY(OrthoBasis[i], Vectors[i]);

    MvarVecWedgeProd(OrthoBasis, Size, &RandomVec, Size + 1, FALSE, NULL);

#ifdef DEBUG
    for (i = 0; i < Size; i++) {
        CagdRType
	    R = MvarVecDotProd(OrthoBasis[i], RandomVec);

	assert(fabs(R) < IRIT_EPS);
    }
#endif /* DEBUG */

    return RandomVec;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the wedge product gradients of multivariates MVs		     *
* at StartPoint and in this direction (TanVec) goes the distance Step.       *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	Multivariate structures to compute the wedge product of      *
*		its gradients at StartPoint				     *
*   StartPoint: Multivariate point to compute tangent line at		     *
*   DirVec:     Control direction that points the same half-space like 	     *
*		the tangent vector					     *
*   Step:       NewPoint = StartPoint + Step * (Unit)TanVec		     *
*   MVGrads:	Gradients of MVs					     *
*   AS:	        Auxiliary MVZR1D structure that holds auxiliary data.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:  NewPoint that lies at the distance Step from the 	     *
*                    the original point StartPoint in the direction of the   *
*	             wedge product vector of gradients of MVS at StartPoint. *
*****************************************************************************/
static MvarPtStruct *MvarMVZR1DStepInTangDir(MvarMVStruct * const *MVs,  
					     MvarPtStruct *StartPoint,
					     MvarVecStruct *DirVec,
					     CagdRType Step,
					     MvarMVGradientStruct **MVGrads,
					     MvarMVZR1DAuxStruct *AS)  
{
    int j,
	Dim = MVs[0] -> Dim;
    MvarPtStruct *NewPoint;
    MvarVecStruct *TanVec;

    NewPoint = MvarPtNew(Dim);

    for (j = 0; j < Dim - 1; j++) {
	CAGD_GEN_COPY(AS -> GradVecs[j] -> Vec, 
		      MvarMVEvalGradient(MVGrads[j], StartPoint -> Pt, 0),
		      sizeof(CagdRType) * Dim);
    }		
     
    TanVec = MvarMVZR1DWedgeProd(AS -> GradVecs, Dim - 1, AS);
    MvarVecNormalize(TanVec);

    /* TanVec has to point into the same half-space like the DirVec. */
    MvarVecScale(TanVec, MvarVecDotProd(TanVec, DirVec) < 0.0 ? -Step : Step);

    for (j = 0; j < NewPoint -> Dim; j++)
        NewPoint -> Pt[j] = StartPoint -> Pt[j] + TanVec -> Vec[j];

    MvarVecFree(TanVec);

#   ifdef DEBUG
    if (!MvarParamsInDomain(MVs[0], NewPoint -> Pt)) {
        fprintf(stderr, "MvarZeroMV1D: Step in the tangent direction is out of domain.\n");
    }
#   endif /* DEBUG */

    return NewPoint;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Correction step of multivariate Newton-Raphson. The intersection of	     *
* N = Dim hyperplanes is computed; Tangent hyperplanes to MVs[i]	     *
* i = 0..N-2 from the StartPoint; the last hyperplane is perpendicular       *
* to all gradient vectors of MVs[i] at StartPoint and passes through	     *
* this point.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	Multivariate structures					     *
*   StartPoint: Starting point from which we approach to the intersection    *
*		curve of MVs					  	     *
*   MVGrads:	Precomputed derivatives of multivariates MVs. In order to    *
*		speed up the computation, they are computed a priori and     * 
*		in the function are only evaluated			     *
*   AS:	        Auxiliary MVZR1D structure that holds auxiliary data.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:   NewPoint						     *
*****************************************************************************/
static MvarPtStruct *MvarMVZR1DCorrectionStep(MvarMVStruct * const *MVs,  
					      MvarPtStruct *StartPoint, 
					      MvarMVGradientStruct **MVGrads,
					      MvarMVZR1DAuxStruct *AS)
{
    int i, j,
	Dim = MVs[0] -> Dim;
    CagdRType *R;
    MvarVecStruct *TanVec;
    MvarPtStruct *NewPoint; 

    NewPoint = MvarPtNew(Dim);
 
    for (j = 0; j < Dim - 1; j++) {
	CAGD_GEN_COPY(AS -> GradVecs[j] -> Vec, 
		      MvarMVEvalGradient(MVGrads[j], StartPoint -> Pt, 0),
		      sizeof(CagdRType) * Dim);
    }	

    TanVec = MvarMVZR1DWedgeProd(AS -> GradVecs, Dim - 1, AS);
    MvarVecNormalize(TanVec);

    /* Creates a Dim * Dim linear system. */
    for (i = 0; i < Dim - 1; i++) {
	IRIT_GEN_COPY(&AS -> A[i * Dim], 
		      AS -> GradVecs[i] -> Vec, Dim * sizeof(CagdRType));
    }
    IRIT_GEN_COPY(&AS -> A[(Dim - 1) * Dim], TanVec -> Vec,
		  Dim * sizeof(CagdRType));

    IRIT_GEN_COPY(AS -> SITTempVec -> Vec, 
		  StartPoint -> Pt, Dim * sizeof(CagdRType));
    
    for (i = 0; i < Dim - 1; i++) {
	R = MvarMVEval(MVs[i], StartPoint -> Pt);
	AS -> b[i] = MvarVecDotProd(AS -> SITTempVec, AS -> GradVecs[i]) -
	             R[1];
    }
    AS -> b[Dim - 1] =  MvarVecDotProd(TanVec, AS -> SITTempVec);
  
    if (IritQRUnderdetermined(AS -> A, NULL, NULL, Dim, Dim)) {
	MvarPtFree(NewPoint);
	MvarVecFree(TanVec);
        return NULL;
    }
    IritQRUnderdetermined(NULL, NewPoint -> Pt, AS -> b, Dim, Dim);

    MvarVecFree(TanVec);

    return NewPoint;  
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test if the Point is "close enough" to the zero level of MVs.	     *
*   Close enough means that function value of MVs at Point is less than      *
*   given tolerance Tol for all MVs.					     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	Multivariate structures.				     *
*   Size:	Size of Mvars to test.					     *
*   Point:	Point, where the MVs are evaluated.			     *
*   Tol:	Numerical tolerance.				  	     *
*									     *
* RETURN VALUE:                                                              *
*   int:	TRUE if Point is "close enough".			     *
*****************************************************************************/
static int MvarMVZR1DCloseToIntersCrv(MvarMVStruct * const *MVs,
				      int Size, 
				      const MvarPtStruct *Point,
				      CagdRType Tol)
{
    int i;
    CagdRType *d;

    for (i = 0; i < Size; i++) {
	d = MvarMVEval(MVs[i], Point -> Pt); 
	if (IRIT_FABS(d[1]) > Tol)
	    return FALSE;
    }

    return TRUE;  
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Evaluate the error at the given position.                                *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	Multivariate structures.				     *
*   Size:	Size of Mvars to test.					     *
*   Pt:         Position where to evaluate the error.                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:   The error.                                                  *
*                                                                            *
* SEE ALSO:                                                                  *
*   MvarMVsZeros                                                             *
*                                                                            *
* KEYWORDS:                                                                  *
*   MvarMVZR1DEvalErr                                                        *
*****************************************************************************/
static CagdRType MvarMVZR1DEvalErr(MvarMVStruct * const *MVs,
				   int Size, 
				   const MvarPtStruct *Pt)
{
    int i;
    CagdRType
        Err = 0.0;

    for (i = 0; i < Size; i++) {
        CagdRType
	    *R = MvarMVEval(MVs[i], Pt -> Pt);

	Err += IRIT_FABS(R[1]);
    }

    return Err;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Traces a intersection curve of (Dim-1) multi-variables MVs in Dim space. *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	Multivariate structures.				     *
*   StartPoint:	Point, where the tracing starts.			     *
*   EndPoint:	Point, where the tracing is ended up.			     *
*   NumericTol:	Numerical tolerance.					     *
*   Step:	Distance in the direction of tangent vector.		     *
*   DirVec:     Control direction that points the same half-space like       *
*		the tangent vector					     *
*   MVGrads:	Gradients of MVs.					     *
*   AS:	        Auxiliary MVZR1D structure that holds auxiliary data.	     *
*   TraceError: Will be set to TRUE if the tracing failed.  Not a globably   *
*               fatal error as it is possible to subdivide & retry to trace. *
*									     *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:    List of MvarPoints that approximate the intersection  *
*		       curve.						     *
*****************************************************************************/
MvarPtStruct *MvarMVZR1DCurveTracing(MvarMVStruct * const *MVs,  
				     const MvarPtStruct *StartPoint,
				     const MvarPtStruct *EndPoint,
				     const MvarVecStruct *DirVec,
				     CagdRType Step, 
				     CagdRType NumericTol,
				     MvarMVGradientStruct **MVGrads,
				     MvarMVZR1DAuxStruct *AS,
				     int *TraceError)
{
    int i, Count,
	Dim = MVs[0] -> Dim;
    MvarPtStruct *TPoint, *TanDirPoint, *CorrPoint, *PtList, *LocStart;
    MvarVecStruct *ControlVec;
    CagdRType Error, NewError,
        Tol = Step * MVAR_MVZR1D_END_TOL_FACTOR;

    *TraceError = FALSE;

#   ifdef DEBUG_MVAR_MVZR1D_DMN_TRACING
    {
	CagdRType t, TMin, TMax;

	for (i = 0; i < MVs[0] -> Dim; i++) {
	    MvarMVDomain(MVs[0], &TMin, &TMax, i);
	    printf("Domain[%d] = [%12.10f  %12.10f]\n", i, TMin, TMax);
	}
	t = MvarMVEval(MVs[0], StartPoint -> Pt)[1];
	printf("\nStart Pt (%10.9lg):", t);
	for (i = 0; i < MVs[0] -> Dim; i++)
	    printf(" %10.9lf", StartPoint -> Pt[i]);
	t = MvarMVEval(MVs[0], EndPoint -> Pt)[1];
	printf("\nEnd   Pt (%10.9lg):", t);
	for (i = 0; i < MVs[0] -> Dim; i++)
	    printf(" %10.9lf", EndPoint -> Pt[i]);
	printf("\n");
    }
#   endif /* DEBUG_MVAR_MVZR1D_DMN_TRACING */

    PtList = MvarPtCopy(StartPoint);

    /*	If Step > Dist(StartPoint, EndPoint), return those 2 points. */
    if (MvarPtDistSqrTwoPoints(StartPoint, EndPoint) < IRIT_SQR(Step)) { 
	TPoint = MvarPtCopy(EndPoint);
	IRIT_LIST_PUSH(TPoint, PtList);
	TPoint = NULL;
	return PtList;	
    }

    LocStart = MvarPtCopy(StartPoint);
    ControlVec = MvarVecCopy(DirVec);
 
    do {
	/* Starting prediction stage - Step motion in tangent direction. */
	TanDirPoint = MvarMVZR1DStepInTangDir(MVs, LocStart, ControlVec, 
					      Step, MVGrads, AS);

	if (MvarParamsInDomain(MVs[0], TanDirPoint -> Pt))
	    TPoint = TanDirPoint;
	else {
	    if ((TPoint = MvarMVIntersPtOnBndry(MVs[0], LocStart,
		                                   TanDirPoint)) == NULL) {
	        MvarPtFree(TanDirPoint);
#	        ifdef DEBUG
	            fprintf(stderr, "Failed to find boundary intersection.\n");
#	        endif /* DEBUG */

		*TraceError = TRUE;
		break;
	    }
	    MvarPtFree(TanDirPoint);
	}

	if (!MvarPtCmpTwoPoints(LocStart, TPoint, NumericTol)) {
#	    ifdef DEBUG
	        fprintf(stderr, "Could move too little in tan dir!\n");
#	    endif /* DEBUG */

	    MvarPtFree(TPoint);
	    *TraceError = TRUE;
	    break;    
	}

	TanDirPoint = MvarPtCopy(TPoint);
	Error = MvarMVZR1DEvalErr(MVs, Dim - 1, TPoint);
#       ifdef DEBUG_MVZR1D_CORRECTION_STAGE
	    fprintf(stderr, "Initial Error is %10.9f\n", Error);
#       endif /* DEBUG_MVZR1D_CORRECTION_STAGE */

	/* Starting correction stage. */
	Count = 1;
	do {
	    if ((CorrPoint = MvarMVZR1DCorrectionStep(MVs, TPoint,
						      MVGrads, AS)) == NULL) {
	        *TraceError = TRUE;
	        break;
	    }

	    /* Test if the correction was successful. */
	    NewError = Error;
	    if (MvarParamsInDomain(MVs[0], CorrPoint -> Pt) &&
		(NewError = MvarMVZR1DEvalErr(MVs, Dim - 1,
					      CorrPoint)) < Error) {
#               ifdef DEBUG_MVZR1D_CORRECTION_STAGE
	            fprintf(stderr, "    Correction Error is %10.9f\n",
			    NewError);
#               endif /* DEBUG_MVZR1D_CORRECTION_STAGE */

		IRIT_GEN_COPY(TPoint -> Pt, CorrPoint -> Pt, 
			      Dim * sizeof(CagdRType));	    	 
	    }
	    else {
	        MvarPtFree(CorrPoint);
		CorrPoint = MvarPtInBetweenPoint(LocStart, TanDirPoint, 0.5);

#		ifdef DEBUG_MVZR1D_CORRECTION_STAGE
		    fprintf(stderr, "Correction step failed.\n");
#		endif /* DEBUG_MVZR1D_CORRECTION_STAGE */

		IRIT_GEN_COPY(TPoint -> Pt, CorrPoint -> Pt, 
			      Dim * sizeof(CagdRType));
		IRIT_GEN_COPY(TanDirPoint -> Pt, CorrPoint -> Pt, 
			      Dim * sizeof(CagdRType));
	    }
	    Error = NewError;

	    /* Cannot correct any more. */
	    if (!MvarPtCmpTwoPoints(LocStart, TPoint, NumericTol)) {
		IRIT_GEN_COPY(TPoint -> Pt, EndPoint -> Pt,
			      Dim * sizeof(CagdRType));

#		ifdef DEBUG
		    fprintf(stderr, "No progress in the tangent direction.\n");
#		endif /* DEBUG */
		
		*TraceError = TRUE;
		break; 
	    }		   
	    MvarPtFree(CorrPoint);
	}
	while (!*TraceError &&
	       Count++ < MVAR_NUMER_ZERO_NUM_STEPS &&
	       !MvarMVZR1DCloseToIntersCrv(MVs, Dim - 1, 
					   TPoint, NumericTol));

	if (Count >= MVAR_NUMER_ZERO_NUM_STEPS) {
#	    ifdef DEBUG
		fprintf(stderr, "Failed to correct after maximally allowed iterations.\n");
#	    endif /* DEBUG */
	    *TraceError = TRUE;
	}

	MvarPtFree(TanDirPoint);

	for (i = 0; i < Dim; i++) {
	    ControlVec -> Vec[i] = TPoint -> Pt[i] - LocStart -> Pt[i];
	}

	IRIT_GEN_COPY(LocStart -> Pt, TPoint -> Pt, Dim * sizeof(CagdRType));
	IRIT_LIST_PUSH(TPoint, PtList);
    }
    while (!*TraceError && MvarPtDistTwoPoints(TPoint, EndPoint) >= Tol);

    MvarPtFree(LocStart);
    MvarVecFree(ControlVec);

    IRIT_GEN_COPY(PtList -> Pt, EndPoint -> Pt, Dim * sizeof(CagdRType));

#   ifdef DEBUG_MVAR_MVZR1D_LINK_NEIGH4
    MvarDbgMVZR1DExamineSegmentBndry(PtList, MVs[0]);
#   endif /* DEBUG_MVAR_MVZR1D_LINK_NEIGH4 */

#   ifdef DEBUG_MVAR_MVZR1D_DMN_TRACING
    {
        printf("\nResult trace:\n");

        for (TPoint = PtList; TPoint != NULL; TPoint = TPoint -> Pnext) {
	    printf("[");
	    for (i = 0; i < MVs[0] -> Dim; i++)
	        printf(" %12.10f ", TPoint -> Pt[i]);
	    printf("]\n");
	}
    }
#   endif /* DEBUG_MVAR_MVZR1D_DMN_TRACING */

#   ifdef MVAR_ZERO_RECOVER_TRACE_ERR
    if (*TraceError && PtList != NULL) {
        /* Lets see if we can recover from this error. */
	if (!MvarPtCmpTwoPoints(CagdListLast(PtList), StartPoint, NumericTol) &&
	    !MvarPtCmpTwoPoints(PtList, EndPoint, NumericTol)) {
	    *TraceError = FALSE;
	}
    }
#   endif /* MVAR_ZERO_RECOVER_TRACE_ERR */

    return PtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   MVs is a system of N-1 constraints in N = Dim unknowns. In some middle   *
* (splitting) hyperplane of its domain, the well constrained system of	     *
*  (Dim-1) equations with (Dim-1) unknowns is solved.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	    (Dim - 1) Multivariate structures.			     *
*   AS:		    Auxiliary MVZR1D structure that holds auxiliary data.    *
*   SubdivTol:	    Subdivision tolerance.				     *
*   NumericTol:	    Numerical tolerance.				     *
*   BoundarySide:   Hyperplane direction.				     *
*   BoundaryValue:  Value, where MVs is evaluated at.			     *
*									     *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:	Intersection points; list of the 0-dimensional	     *
*			solutions of the system MVs in the middle	     *
*			hyperplane.					     *
*****************************************************************************/
MvarPtStruct *MvarMVZR1DMiddlePlaneCutPts(MvarMVStruct * const *MVs,
					  MvarMVZR1DAuxStruct *AS,
					  CagdRType SubdivTol,
					  CagdRType NumericTol,
					  int BoundarySide, 
					  CagdRType BoundaryValue)
{
    int j, k,
	Dim = MVs[0] -> Dim;
    MvarPtStruct *Aux,
	*MiddlePts = NULL,
	*ExtendedMiddlePts = NULL;
    MvarMVStruct
	**MiddleMVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) 
						                * (Dim - 1));

    for (j = 0; j < Dim - 1; j++)
	MiddleMVs[j] = MvarMVFromMV(MVs[j], BoundaryValue, BoundarySide);

    MiddlePts = MvarMVsZeros0D(MiddleMVs, AS -> Constraints,
			       Dim - 1, SubdivTol, -IRIT_FABS(NumericTol));

    /* The dimension of middle points is extended by 1 to Dim. */
    for (Aux = MiddlePts; Aux != NULL; Aux = Aux -> Pnext) {
        MvarPtStruct 
	    *NewPt = MvarPtNew(Dim);

	for (k = 0; k < BoundarySide; k++)
	    NewPt -> Pt[k] = Aux -> Pt[k];	    
	NewPt -> Pt[BoundarySide] = BoundaryValue;
	for (k = BoundarySide + 1; k < Dim; k++)
	    NewPt -> Pt[k] = Aux -> Pt[k-1];	    

	IRIT_LIST_PUSH(NewPt, ExtendedMiddlePts);
    }

    for (j = 0; j < Dim - 1; j++)
	MvarMVFree(MiddleMVs[j]);
    IritFree(MiddleMVs);

    MvarPtFreeList(MiddlePts);

    if (CagdListLength(ExtendedMiddlePts) > 1)
	return MvarMVZR1DListOfDifferentPts(ExtendedMiddlePts, NumericTol);
    else
	return ExtendedMiddlePts;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Boundary points are split into two lists. Splitting hyperplane is	     *
* perpendicular to BoundarySide axis and passes through BoundaryValue        *
* point.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   BoundaryPts:           Points to split.				     *
*   BoundarySide:          Hyperplane direction.			     *
*   BoundaryValue:         Value, where MVs is evaluated at.		     *
*   SplitPts0, SplitPts1:  The two new list.				     *
*									     *
* RETURN VALUE:                                                              *
*   void			                                             *
*****************************************************************************/
void MvarMVZR1DSplitBoundaryPts(const MvarPtStruct *BoundaryPts,
				int BoundarySide, 
				CagdRType BoundaryValue,
				MvarPtStruct **SplitPts0,
				MvarPtStruct **SplitPts1)
{
    *SplitPts0 = *SplitPts1 = NULL;

    if (BoundaryPts == NULL)
	return;
    else {
        MvarPtStruct *TempPt, *TempPt2,
	    *BoundaryPtsCp = MvarPtCopyList(BoundaryPts);

	while (BoundaryPtsCp != NULL) {
	    IRIT_LIST_POP(TempPt, BoundaryPtsCp);

	    if (IRIT_APX_EQ_EPS(TempPt -> Pt[BoundarySide], BoundaryValue,
				IRIT_UEPS)) {
	        TempPt-> Pt[BoundarySide] = BoundaryValue;
	        TempPt2 = MvarPtCopy(TempPt);
		IRIT_LIST_PUSH(TempPt, *SplitPts0);
		IRIT_LIST_PUSH(TempPt2, *SplitPts1);
	    }
	    else if (TempPt -> Pt[BoundarySide] < BoundaryValue){
		IRIT_LIST_PUSH(TempPt, *SplitPts0);
	    }
	    else {
		IRIT_LIST_PUSH(TempPt, *SplitPts1);
	    }	    
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes a vector, which points to the same half-space as the tangent    *
*   vector of the intersection curve (solution of MVs[i] = 0, i=0,..,N) at   *
*   the boundary point BoundaryPt.					     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	The system of N = Dim - 1 equations with N + 1 unknowns	     *
*		" MVs[i] = 0, i=0,..,N ".				     *
*   BoundaryPt: Point on the boundary of the domain of MVs                   *
*   AS:	        Auxiliary MVZR1D structure that holds auxiliary data.	     *
*   NumericTol: Numeric Tolerance.					     *
*									     *
* RETURN VALUE:                                                              *
*   MvarVecStruct *:	Vector that points into the domain of MVs.	     *
*****************************************************************************/
MvarVecStruct *MvarMVZR1DStartVec(MvarMVStruct * const *MVs, 
				  MvarPtStruct *BoundaryPt,
				  MvarMVZR1DAuxStruct *AS,
				  CagdRType NumericTol)
{   
    int i,
	Dim = MVs[0] -> Dim;
    MvarVecStruct
	*StartVec = MvarVecNew(Dim);

    MVAR_VEC_RESET(StartVec);

    MvarMVDomain(MVs[0], AS -> MinDmn, AS -> MaxDmn, -1);

    for (i = 0; i < Dim; i++) {
	if (BoundaryPt -> Pt[i] < AS -> MinDmn[i] + NumericTol)
	    StartVec -> Vec[i] = 1.0;
	else if (BoundaryPt -> Pt[i] > AS -> MaxDmn[i] - NumericTol)
	    StartVec -> Vec[i] = -1.0;
    }

    return StartVec;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Tests if the hyper-surfaces MVs are flat enough such that the	     *
*   intersection curve has no loop. Note, the intersection curve can	     *
*   consist of more than one segment.					     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	The system of N = Dim - 1 equations with N + 1 unknowns	     *
*		" MVs[i] = 0, i=0,..,N "				     *
*   AS:	        Auxiliary MVZR1D structure that holds auxiliary data.	     *
*									     *
* RETURN VALUE:                                                              *
*   int:	TRUE if the intersection curve has no loops.		     *
*****************************************************************************/
int MvarMVZR1DNoLoopTest(MvarMVStruct const * const *MVs, 
		         MvarMVZR1DAuxStruct *AS)
{
    int i, j, k, 
	Dim = MVs[0] -> Dim,
	PowerTwoDim = (int) pow(2, Dim - 1);
    MvarVecStruct *WedProd;
      
    for (i = 0; i < Dim - 1; i++) {
        MvarNormalConeStruct *Cone;

        if ((Cone = MVarMVNormalCone(MVs[i])) == NULL)	   
            return FALSE;

	CAGD_GEN_COPY(AS -> TempList[i] -> Vec, 
		      Cone -> ConeAxis -> Vec, sizeof(CagdRType) * Dim);

        /* Cone is valid: add plane orthogonal to cone axis to matrix A. */
        CAGD_GEN_COPY(&AS -> A[i * (Dim)], Cone -> ConeAxis -> Vec,
		      sizeof(CagdRType) * Dim);

        /* Take the anti-cone = cos(90 - angle) = sqrt(1-cos^2). */
        AS -> b[i] = sqrt(1.0 - IRIT_SQR(Cone -> ConeAngleCosine));

        MvarNormalConeFree(Cone);
    }

    /* Axes of tangent cones are linearly dependent. */
    if (MvarMVZR1DLinDep(AS -> TempList, Dim - 1, AS))
	 return FALSE;

    WedProd = MvarMVZR1DWedgeProd(AS -> TempList, Dim - 1, AS); 
    CAGD_GEN_COPY(&AS -> A[(Dim - 1) * Dim], WedProd -> Vec, 
		  sizeof(CagdRType) * Dim);
    MvarVecFree(WedProd);

    AS -> b[Dim - 1] = 0;

    IritQRUnderdetermined(AS -> A, NULL, NULL, Dim, Dim);
   
    for (i = 0; i < PowerTwoDim; i++) {
        CagdRType VecSqrLength;

        k = i;
        for (j = 0; j < Dim - 1; ++j) {
	    AS -> bCopy[j] = (k & 1) ? AS -> b[j] : -AS -> b[j];
            k >>= 1;
        }
	AS -> bCopy[Dim - 1] = 0;

        IritQRUnderdetermined(NULL, AS -> x, AS -> bCopy, Dim, Dim);

        for (VecSqrLength = 0.0, j = 0; j < Dim; ++j)
            VecSqrLength += IRIT_SQR(AS -> x[j]);

        if (VecSqrLength >= 1.0)
            return FALSE;
    }

    return TRUE;    
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   In order to plot domains where MvarMVZR1DCurveTracing was applied,       *
* polyline is provided by the domain. Used for debugging.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly:	Polyline, Poly -> PAux points to the domain.		     *
*   MV:	        Multivar structure with the domain.			     *
*									     *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *:	Polyline provided by domain.		     *
*****************************************************************************/
MvarPolylineStruct *MvarMVZR1DPolyWithDom(MvarPolylineStruct *Poly, 
  				          const MvarMVStruct *MV)
{
#ifdef DEBUG_MVAR_MVZR1D_SAVE_POLY_DOM
    int i,
        Dim = MV -> Dim;
    CagdRType *Domain;
   
    Domain = IritMalloc(2 * Dim * sizeof(CagdRType));

    MvarMVDomain(MVs[0], AS -> MinDmn2, AS -> MaxDmn2, -1);

    IRIT_GEN_COPY(&Domain[0], AS -> MinDmn2, Dim * sizeof(CagdRType));
    IRIT_GEN_COPY(&Domain[Dim], AS -> MaxDmn2, Dim * sizeof(CagdRType));

    Poly -> PAux = Domain;
#endif /* DEBUG_MVAR_MVZR1D_SAVE_POLY_DOM */

    return Poly;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates a polyn. system for SSI (surface-surface intersection) problem.  *
*   Two parametric surfaces MV1(u1,v1) and MV2(u2,v2) are given. The system  *
*   of 3 equations with 4 unknowns is created:		    		     *
*	x(u1,v1) - x(u1,v1) = 0						     *
*	y(u1,v1) - y(u1,v1) = 0						     *
*       z(u1,v1) - z(u1,v1) = 0.                                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Srf1, Srf2:	    Surfaces to be intersected.				     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct **:  The multivar system, its solution is desired curve.    *
*                     A static array of size 3 that returns 3 MVs,           *
*                     dynamically allocated.				     *
*****************************************************************************/
static MvarMVStruct **MvarMVZR1DCreateMVs(const CagdSrfStruct *Srf1, 
					  const CagdSrfStruct *Srf2)
{
    int i;
    CagdRType Min, Max;
    static MvarMVStruct *MVs[3];
    MvarMVStruct *TMV, *MVX1, *MVY1, *MVZ1, *MVW1, **TempMVs,
        *MVX2, *MVY2, *MVZ2, *MVW2, *MVTmp,
        *MV1 = MvarSrfToMV(Srf1),
        *MV2 = MvarSrfToMV(Srf2);
  
    if (MVAR_IS_BEZIER_MV(MV1)) {
	MVTmp = MvarCnvrtBzr2BspMV(MV1);
	MvarMVFree(MV1);
	MV1 = MVTmp;
    }

    if (MVAR_IS_BEZIER_MV(MV2)) {
	MVTmp = MvarCnvrtBzr2BspMV(MV2);
	MvarMVFree(MV2);
	MV2 = MVTmp;
    }

    TMV = MvarPromoteMVToMV2(MV1, 4, 0);
    MvarMVFree(MV1);
    MV1 = TMV;

    TMV = MvarPromoteMVToMV2(MV2, 4, 2);
    MvarMVFree(MV2);
    MV2 = TMV;
   
    for (i = 0; i < 2; i++) {
	MvarMVDomain(MV1, &Min, &Max, i);
	BspKnotAffineTrans2(MV2 -> KnotVectors[i],
			    MV2 -> Lengths[i] + MV2 -> Orders[i],
			    Min, Max);
    }

    for (i = 2; i < 4; i++) {
	MvarMVDomain(MV2, &Min, &Max, i);
	BspKnotAffineTrans2(MV1 -> KnotVectors[i],
			    MV1 -> Lengths[i] + MV1 -> Orders[i],
			    Min, Max);
    }

    TempMVs = MvarMVSplitScalar(MV1);
    MVW1 = TempMVs[0];
    MVX1 = TempMVs[1];
    MVY1 = TempMVs[2];
    MVZ1 = TempMVs[3];

    TempMVs = MvarMVSplitScalar(MV2);
    MVW2 = TempMVs[0];
    MVX2 = TempMVs[1];
    MVY2 = TempMVs[2];
    MVZ2 = TempMVs[3];

    MvarMVFree(MV1);
    MvarMVFree(MV2);

    if (MVW1 != NULL) {
	TMV = MvarMVMult(MVW1, MVX2); /* X2 = W1*X2. */
	MvarMVFree(MVX2);
	MVX2 = TMV;

	TMV = MvarMVMult(MVW1, MVY2); /* Y2 = W1*Y2. */
	MvarMVFree(MVY2);
	MVY2 = TMV;

	TMV = MvarMVMult(MVW1, MVZ2); /* Z2 = W1*Z2. */
	MvarMVFree(MVZ2);
	MVZ2 = TMV;

	MvarMVFree(MVW1);
    }

    if (MVW2 != NULL) {
	TMV = MvarMVMult(MVW2, MVX1); /* X1 = W2*X1. */
	MvarMVFree(MVX1);
	MVX1 = TMV;

	TMV = MvarMVMult(MVW2, MVY1); /* Y1 = W2*Y1. */
	MvarMVFree(MVY1);
	MVY1 = TMV;

	TMV = MvarMVMult(MVW2, MVZ1); /* Z1 = W2*Z1. */
	MvarMVFree(MVZ1);
	MVZ1 = TMV;

	MvarMVFree(MVW2);
    }

    MVs[0] = MvarMVSub(MVX1, MVX2);
    MVs[1] = MvarMVSub(MVY1, MVY2);
    MVs[2] = MvarMVSub(MVZ1, MVZ2);

    MvarMVFree(MVX1);
    MvarMVFree(MVY1);
    MvarMVFree(MVZ1);
    MvarMVFree(MVX2);
    MvarMVFree(MVY2);
    MvarMVFree(MVZ2);

    return MVs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes intersection curve of two surfaces.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf1, Srf2:	    Cagd surfaces to be intersected.			     M
*   Step:	    Step size for curve tracing.			     M
*   SubdivTol:	    The subdivision tolerance to use.			     M
*   NumericTol:	    The numerical tolerance to use.			     M
*									     *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: The list of polylines which approximate the curve. M
*			  Each polyline corresponds to the topologically     M
*			  isolated component of the curve and is in R^4, the M
*                         parametric spaces of both surfaces.   	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros1DMergeSingularPts, MvarMVsZeros1D                           M
*									     *
* KEYWORDS:                                                                  M
*   MvarSrfSrfInter, MvarSrfZeroSet	                                     M
*****************************************************************************/
MvarPolylineStruct *MvarSrfSrfInter(const CagdSrfStruct *Srf1, 
				    const CagdSrfStruct *Srf2,
				    CagdRType Step,
				    CagdRType SubdivTol,
				    CagdRType NumericTol)
{
    int i;
    CagdBBoxStruct BBox1, BBox2;
    MvarPolylineStruct *IntrCrv;
    MvarMVStruct **MVs;
    CagdSrfStruct
	*NewSrf1 = NULL,
	*NewSrf2 = NULL;

    if (CAGD_IS_PERIODIC_SRF(Srf1) ||
	CAGD_IS_PERIODIC_SRF(Srf2)) {
        MVAR_FATAL_ERROR(MVAR_ERR_RATIONAL_NO_SUPPORT);
	return NULL;
    }

    /* Examine bounding boxes for no-overlap for quick pruning. */
    CagdSrfBBox(Srf1, &BBox1);
    CagdSrfBBox(Srf2, &BBox2);
    if (BBox1.Max[0] < BBox2.Min[0] ||
	BBox1.Max[1] < BBox2.Min[1] ||
	BBox1.Max[2] < BBox2.Min[2] ||
	BBox2.Max[0] < BBox1.Min[0] ||
	BBox2.Max[1] < BBox1.Min[1] ||
	BBox2.Max[2] < BBox1.Min[2])
	return NULL;

    if (Srf1 -> PType != CAGD_PT_E3_TYPE && Srf1 -> PType != CAGD_PT_P3_TYPE)
        Srf1 = NewSrf1 = CagdCoerceSrfTo(Srf1,
					 CAGD_IS_RATIONAL_PT(Srf1 -> PType) ?
							  CAGD_PT_P3_TYPE : 
						          CAGD_PT_E3_TYPE,
				FALSE);

    if (Srf2 -> PType != CAGD_PT_E3_TYPE && Srf2 -> PType != CAGD_PT_P3_TYPE)
        Srf2 = NewSrf2 = CagdCoerceSrfTo(Srf2,
					 CAGD_IS_RATIONAL_PT(Srf2 -> PType) ?
							  CAGD_PT_P3_TYPE : 
						          CAGD_PT_E3_TYPE,
				FALSE);

    MVs = MvarMVZR1DCreateMVs(Srf1, Srf2);

    if (NewSrf1 != NULL)
        CagdSrfFree(NewSrf1);
    if (NewSrf2 != NULL)
        CagdSrfFree(NewSrf2);

#   ifdef DEBUG_MVAR_MVZR1D_MV_ZEROS
    {
        /* Test solution with the oD solver. */
        MvarPtStruct
	    *Pts = MvarMVsZeros0D(MVs, NULL, 3, SubdivTol, NumericTol);

	IntrCrv = MvarMatchPointListIntoPolylines(Pts, SubdivTol * 10);
    }
#   else
        IntrCrv = MvarMVsZeros1D(MVs, NULL, 3, Step, SubdivTol,
				 IRIT_ABS(NumericTol));
#   endif /* DEBUG_MVAR_MVZR1D_MV_ZEROS */

    for (i = 0; i < 3; i++)
        MvarMVFree(MVs[i]);	    	    

    return IntrCrv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the zeros of some surface.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Surface:	    Cagd surface to compute its zeros.			     M
*   Axis:	    Axis of Surface to consider its zeros.		     M
*   Step:	    Step size for curve tracing.			     M
*   SubdivTol:	    The subdivision tolerance to use.			     M
*   NumericTol:	    The numerical tolerance to use.			     M
*									     *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: The list of polylines which approximate the zeros. M
*			  Each polyline corresponds to the topologically     M
*			  isolated component of the curve and is in R^2, the M
*                         parametric spaces of the surface.   		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCrvZeroSet, MvarMVsZeros1D, MvarSrfSrfInter		             M
*									     *
* KEYWORDS:                                                                  M
*   MvarSrfZeroSet			                                     M
*****************************************************************************/
MvarPolylineStruct *MvarSrfZeroSet(const CagdSrfStruct *Surface,
			           int Axis,
			           CagdRType Step,
			           CagdRType SubdivTol,
			           CagdRType NumericTol)
{
    int i;
    MvarPolylineStruct *IntrCrv;
    MvarMVStruct *MVs[1];
    CagdSrfStruct *Srf[4];

    SymbSrfSplitScalar(Surface, &Srf[0], &Srf[1], &Srf[2], &Srf[3]);

    /* Converts to a multivariate. */
    MVs[0] = MvarSrfToMV(Srf[Axis]);

    for (i = 0; i <= 3; i++) {
        if (Srf[i] != NULL)
	    CagdSrfFree(Srf[i]);
    }

#   ifdef DEBUG_MVAR_MVZR1D_MV_ZEROS
    {
        /* Test solution with the 0D solver. */
        MvarPtStruct
	    *Pts = MvarMVsZeros0D(MVs, NULL, 1, SubdivTol, NumericTol);

	IntrCrv = MvarMatchPointListIntoPolylines(Pts, SubdivTol * 10);
    }
#   else
        IntrCrv = MvarMVsZeros1D(MVs, NULL, 1, Step, SubdivTol,
				 IRIT_ABS(NumericTol));
#   endif /* DEBUG_MVAR_MVZR1D_MV_ZEROS */

    MvarMVFree(MVs[0]);	    	    

    return IntrCrv;
}
