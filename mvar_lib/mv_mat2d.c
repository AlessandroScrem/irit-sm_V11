/******************************************************************************
* MV_MAT2D.c - computes the MAT of 2D closed curves.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, July 14.					      *
******************************************************************************/

#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "mvar_loc.h"

#define DEBUG_MAT_RETURN_TRITANS
#define DEBUG_MAT_RETURN_SPIRALS

#define MVAR_MAT_SPLIT_TOL	 1e-6
#define MVAR_MAT_GENERAL_TOL	 1e-10
#define MVAR_MAT_SOL_ALLOC_SIZE	 4

typedef struct MvarC2DMTriTanCircStruct {
    CagdPType Params;     /* The 3 parameters of the tri tangency location. */
    int Idx;            /* This curve is first, second, or third parameter. */
    CagdPType CircleCenter;
    CagdRType CircleRadius;
} MvarC2DMTriTanCircStruct;

typedef struct MvarC2DMTriTanCircArrayStruct {
    int NumTriTans;
    int MaxNumTriTans;
    MvarC2DMTriTanCircStruct *TriTans;              /* Array of parameters. */
} MvarC2DMTriTanCircArrayStruct;

static int MvarC2DMTriTans(const CagdCrvStruct *OCrv,
			   CagdCrvStruct *Crv1,
			   CagdCrvStruct *Crv2,
			   CagdCrvStruct *Crv3,
			   CagdRType SubdivTol,
			   CagdRType NumericTol,
			   CagdBType OneSideOrientation);
static CagdBType MvarC2DMSolutionVerified(const CagdCrvStruct *OCrv,
					  const CagdRType *Sol,
					  CagdPtStruct *Cntr,
					  CagdRType *Radius);
static int MvarC2DMUpdateTriTans(CagdCrvStruct *Crv,
				 const MvarPtStruct *MVPt,
				 int Idx,
				 const CagdPtStruct *Cntr,
				 CagdRType Radius);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the (inside) 2D MAT (medial axis transform) of a given closed   M
* planar self-intersection-free oriented B-spline curve.	             M
*                                                                            *
* PARAMETERS:                                                                M
*   OCrv:          B-spline curve to compute its 2D planar MAT.		     M
*   SubdivTol:     Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   NumericTol:    Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*   InvertOrientation:  Flip what is considered inside and outside.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   The MAT as a list of curves.	                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCrv2DMAT                                                             M
*****************************************************************************/
CagdCrvStruct *MvarCrv2DMAT(const CagdCrvStruct *OCrv,
			    CagdRType SubdivTol,
			    CagdRType NumericTol,
			    CagdBType InvertOrientation)
{
    int BndryCond;
    CagdCrvStruct *Crvs, *Crv, *NewCrvs, *Crv1,
        *OCrvRev = NULL;
#   if defined(DEBUG_MAT_RETURN_TRITANS) || defined(DEBUG_MAT_RETURN_SPIRALS)
    CagdCrvStruct *DebugCrv,
        *DebugRetCrvs = NULL;
#   endif /* DEBUG_MAT_RETURN_TRITANS || DEBUG_MAT_RETURN_SPIRALS */


    if (!CAGD_IS_BSPLINE_CRV(OCrv) || !CagdIsClosedCrv(OCrv)) {
        return NULL;
    }

    /* Split Crv at all C^1 discontinuities. */
    if (InvertOrientation) {
        OCrv = OCrvRev = CagdCrvReverse(OCrv);
        Crvs = CagdCrvSubdivAtAllC1Discont(OCrvRev, FALSE, 0.0);
    }
    else {
        Crvs = CagdCrvSubdivAtAllC1Discont(OCrv, FALSE, 0.0);
    }

    /* Split Crv into spirals. */
    NewCrvs = NULL;
    while (Crvs != NULL) {
        IRIT_LIST_POP(Crv, Crvs);

	if (Crv -> Order > 2) {
	    CagdPtStruct
	        *Pts = SymbCrvExtremCrvtrPts(Crv, MVAR_MAT_SPLIT_TOL);

	    if (Pts != NULL) {
		NewCrvs = CagdListAppend(NewCrvs,
				         CagdCrvSubdivAtParams(Crv, Pts,
							       IRIT_EPS,
							       &BndryCond));
		CagdPtFreeList(Pts);
		CagdCrvFree(Crv);
	    }
	    else {
	        NewCrvs = CagdListAppend(NewCrvs, Crv);
	    }
	}
	else {
	    NewCrvs = CagdListAppend(NewCrvs, Crv);
	}
    }
    Crvs = NewCrvs;

    /* Having only spiral segments now, iterate over all triplets and find  */
    /* all circles that are tangent to 3 different curves in the proper     */
    /* orientation.  Note a circle can be tangent to a spiral once only.    */
    for (Crv1 = Crvs; Crv1 != NULL; Crv1 = Crv1 -> Pnext) {
        CagdCrvStruct *Crv2, *Crv3;

        for (Crv2 = Crv1 -> Pnext; Crv2 != NULL; Crv2 = Crv2 -> Pnext) {
	    for (Crv3 = Crv2 -> Pnext; Crv3 != NULL; Crv3 = Crv3 -> Pnext) {
#		ifdef DEBUG_MAT_RETURN_TRITANS
		if ((DebugCrv = AttrGetRefPtrAttrib(Crv1 -> Attr, "_tricirc"))
		                                                    != NULL) {
		    IRIT_LIST_PUSH(DebugCrv, DebugRetCrvs);
		    AttrFreeOneAttribute(&Crv1 -> Attr, "_tricirc");
		}
#		endif /* DEBUG_MAT_RETURN_TRITANS */
	    }
	}
    }

    if (OCrvRev != NULL)
	CagdCrvFree(OCrvRev);

#   if defined(DEBUG_MAT_RETURN_TRITANS) || defined(DEBUG_MAT_RETURN_SPIRALS)
	Crvs = CagdListAppend(Crvs, DebugRetCrvs);
	for (Crv = Crvs; Crv != NULL; Crv = Crv -> Pnext) {
	    AttrSetRGBColor(&Crv -> Attr, (int) IritRandom(100, 255),
					  (int) IritRandom(100, 255),
					  (int) IritRandom(100, 255));
	}
	return Crvs;
#   else
	CagdCrvFreeList(Crvs);
	return NULL;
#   endif /* DEBUG_MAT_RETURN_TRITANS || DEBUG_MAT_RETURN_SPIRALS */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes all tri-tangent circles to a given set of three spiral curves,  *
* Crv?.                                                                      *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv1, Crv2, Crv3:  Three spiral curves to compute tri-tangent circles    *
*                      for.						     *
*   SubdivTol:     Tolerance of the subdivision process.  Tolerance is       *
*		   measured in the parametric space of the multivariates.    *
*   NumericTol:    Numeric tolerance of the numeric stage.  The numeric      *
*		   stage is employed only if NumericTol < SubdivTol.         *
*   OneSideOrientation:   TRUE to compute tri-tangencies on one side of the  *
*                  curves only.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  Number of valid solutions found.                                   *
*****************************************************************************/
static int MvarC2DMTriTans(const CagdCrvStruct *OCrv,
			   CagdCrvStruct *Crv1,
			   CagdCrvStruct *Crv2,
			   CagdCrvStruct *Crv3,
			   CagdRType SubdivTol,
			   CagdRType NumericTol,
			   CagdBType OneSideOrientation)
{
    int n = 0;
    MvarPtStruct *MVPt,
        *MVPts = MvarCircTanTo3Crvs(Crv1, Crv2, Crv3, SubdivTol, NumericTol,
				    OneSideOrientation);

    /* Spread the valid solutions into the different curves. */
    for (MVPt = MVPts; MVPt != NULL; MVPt = MVPt -> Pnext) {
	CagdRType Radius;
	CagdPtStruct Cntr;

	/* Only consider completely interior tri-tangent circles. */
	if (!MvarC2DMSolutionVerified(OCrv, MVPt -> Pt, &Cntr, &Radius))
	    continue;

	n++;
	MvarC2DMUpdateTriTans(Crv1, MVPt, 0, &Cntr, Radius);
	MvarC2DMUpdateTriTans(Crv2, MVPt, 1, &Cntr, Radius);
	MvarC2DMUpdateTriTans(Crv3, MVPt, 2, &Cntr, Radius);

#	ifdef DEBUG_MAT_RETURN_TRITANS
	{
	    AttrSetRefPtrAttrib(&Crv1 -> Attr, "_tricirc",
				BspCrvCreateCircle(&Cntr, Radius));
	}
#	endif /* DEBUG_MAT_RETURN_TRITANS */
    }

    MvarPtFreeList(MVPts);

    return n;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Verify if this tri-tangent circle solution is valid - that is, it does   *
* not contain any other portion of any of the three curves.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol:               Solution of a tri-tangency.                           *
*   Crv1, Crv2, Crv3:  The three circles to test if circle intersects with.  *
*   Cntr:              Will be updated with the circle center if verified.   *
*   Radius:            Will be updated with the circle radius if verified.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:      TRUE if solution valid, FALSE if not.                    *
*****************************************************************************/
static CagdBType MvarC2DMSolutionVerified(const CagdCrvStruct *OCrv,
					  const CagdRType *Sol,
					  CagdPtStruct *Cntr,
					  CagdRType *Radius)
{
    CagdRType *R, t, DistSqr;
    CagdPType Pt1, Pt2, Pt3, Pt;

    /* Update the position on three original primitive curves. */
    R = CagdCrvEval(OCrv, Sol[0]);
    CagdCoerceToE3(Pt1, &R, -1, OCrv -> PType);

    R = CagdCrvEval(OCrv, Sol[1]);
    CagdCoerceToE3(Pt2, &R, -1, OCrv -> PType);

    R = CagdCrvEval(OCrv, Sol[2]);
    CagdCoerceToE3(Pt3, &R, -1, OCrv -> PType);

    if (GMCircleFrom3Points(Cntr -> Pt, Radius, Pt1, Pt2, Pt3)) {
        DistSqr = IRIT_SQR(*Radius - MVAR_MAT_GENERAL_TOL);

        /* Check minimal distance of OCrv to circ center against Radius. */
        t = SymbDistCrvPoint(OCrv, Cntr -> Pt, TRUE, MVAR_MAT_GENERAL_TOL);
	R = CagdCrvEval(OCrv, t);
	CagdCoerceToE3(Pt, &R, -1, OCrv -> PType);
	return IRIT_PT2D_DIST_SQR(Pt, Cntr -> Pt) >= DistSqr;
    }
    else
        return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Updates the curve's local data structure with a new tri-tangency.        *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv:    Curve to update with a new tri-tangency information.	     *
*   MVPt:   Parameters of the tri-tangency.				     *
*   Idx:    Index of parameter of Crv: 0, 1, or 2.			     *
*   Cntr:   The circle center.						     *
*   Radius: The circle radius.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  Updated number of tri-tangencies in curve.                         *
*****************************************************************************/
static int MvarC2DMUpdateTriTans(CagdCrvStruct *Crv,
				 const MvarPtStruct *MVPt,
				 int Idx,
				 const CagdPtStruct *Cntr,
				 CagdRType Radius)
{
    MvarC2DMTriTanCircArrayStruct *Tans;

    /* Fetches the data structure and create new one if non found. */
    if ((Tans = (MvarC2DMTriTanCircArrayStruct *)
	            AttrGetRefPtrAttrib(Crv -> Attr, "_TriTans")) == NULL) {
         /* Create a new structure. */
        Tans = (MvarC2DMTriTanCircArrayStruct *)
	                    IritMalloc(sizeof(MvarC2DMTriTanCircArrayStruct));
	Tans -> TriTans = (MvarC2DMTriTanCircStruct *)
	                         IritMalloc(sizeof(MvarC2DMTriTanCircStruct) *
			                             MVAR_MAT_SOL_ALLOC_SIZE);
	Tans -> NumTriTans = 0;
	Tans -> MaxNumTriTans = MVAR_MAT_SOL_ALLOC_SIZE;
    }
    if (Tans -> NumTriTans == Tans -> MaxNumTriTans) {
        /* Reallocate more space. */
        Tans -> TriTans = (MvarC2DMTriTanCircStruct *)
	                   IritRealloc(Tans -> TriTans,
				       sizeof(MvarC2DMTriTanCircStruct) *
				              Tans -> MaxNumTriTans,
				       sizeof(MvarC2DMTriTanCircStruct) * 2 *
	                                      Tans -> MaxNumTriTans);
	Tans -> MaxNumTriTans *= 2;
    }
    IRIT_PT_COPY(Tans -> TriTans[Tans -> NumTriTans].Params, MVPt -> Pt);
    Tans -> TriTans[Tans -> NumTriTans].Idx = Idx;
    Tans -> TriTans[Tans -> NumTriTans].CircleRadius = Radius;
    IRIT_PT_COPY(Tans -> TriTans[Tans -> NumTriTans].CircleCenter, Cntr -> Pt);

    Tans -> NumTriTans++;

    AttrSetRefPtrAttrib(&Crv -> Attr, "_TriTans", Tans);

    return Tans -> NumTriTans;
}
