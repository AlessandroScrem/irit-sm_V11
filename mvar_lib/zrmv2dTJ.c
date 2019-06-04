/******************************************************************************
* ZrMV2DTJ.c - tools to handle T-Junction events, as part of the zero sets    *
*            computation of multivariates when the problem has MVs            *
*            representation and the expected solution  set is a two-manifold: *
*            surface(s).                                                      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, June 2014.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/geom_lib.h"

static CagdBType MvarZeroSearchTJInLoop(const MvarZeroTJunctionStruct* TJ,
					MvarPolylineStruct *Loop);
#if defined(DEBUG) && defined(DEBUG_TJ_LIST)
    static void MvarZeroSaveTJList(MvarZeroTJunctionStruct *TJList);
    static void MvarZeroPrintTJList(MvarZeroTJunctionStruct *TJList);
#endif /* DEBUG && DEBUG_TJ_LIST */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new multi-variate T-Junction.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   TJPrev:	The point before the junction.				     M
*   TJPt:	The T-Junction point.				             M
*   TJNext:	The point after the junction.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroTJunctionStruct *:    the T-Junction.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroTJFree, MvarZeroTJCopy                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroTJNew                                                            M
*****************************************************************************/
MvarZeroTJunctionStruct *MvarZeroTJNew(const MvarPtStruct *TJPrev,
				       const MvarPtStruct *TJPt,
				       const MvarPtStruct *TJNext)
{
    MvarZeroTJunctionStruct
	*TJ = (MvarZeroTJunctionStruct *) 
	                           IritMalloc(sizeof(MvarZeroTJunctionStruct));

    TJ -> Pnext = NULL;
    TJ -> IsHandled = FALSE;
    TJ -> TJuncPrev = MvarPtCopy(TJPrev);
    TJ -> TJunc = MvarPtCopy(TJPt);
    TJ -> TJuncNext = MvarPtCopy(TJNext);

    return TJ;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Copies a multi-variate T-Junction object.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   TJ:	The T-Junction object to copy.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroTJunctionStruct *:    the T-Junction copy.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroTJFree                                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroTJCopy                                                           M
*****************************************************************************/
MvarZeroTJunctionStruct *MvarZeroTJCopy(const MvarZeroTJunctionStruct *TJ)
{
    MvarZeroTJunctionStruct
	*TJCpy = (MvarZeroTJunctionStruct *) 
				  IritMalloc(sizeof(MvarZeroTJunctionStruct));

    TJCpy -> Pnext = NULL;
    TJCpy -> IsHandled = TJ -> IsHandled;
    TJCpy -> TJuncPrev = MvarPtCopy(TJ -> TJuncPrev);
    TJCpy -> TJunc = MvarPtCopy(TJ -> TJunc);
    TJCpy -> TJuncNext = MvarPtCopy(TJ -> TJuncNext);
    return TJCpy;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Copies a list of multi-variate T-Junction objects.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   TJList:	The T-Junction list to copy.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroTJunctionStruct *: The T-Junctions list copy.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroTJFree, MvarZeroTJCopy                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroTJCopyList                                                       M
*****************************************************************************/
MvarZeroTJunctionStruct *MvarZeroTJCopyList(
                                        const MvarZeroTJunctionStruct *TJList)
{
    MvarZeroTJunctionStruct *TJTemp, *NewTJList;

    if (TJList == NULL)
	return NULL;
    TJTemp = NewTJList = MvarZeroTJCopy(TJList);
    TJList = TJList -> Pnext;
    while (TJList) {
	TJTemp -> Pnext = MvarZeroTJCopy(TJList);
	TJTemp = TJTemp -> Pnext;
	TJList = TJList -> Pnext;
    }
    return NewTJList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Frees all slots of a multi-variate T-Junction structure.                 M
*                                                                            *
* PARAMETERS:                                                                M
*   TJ:      Multivariate T-Junction to free.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroTJFreeList                                                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroTJFree                                                           M
*****************************************************************************/
void MvarZeroTJFree(MvarZeroTJunctionStruct *TJ)
{
    if (TJ -> TJuncPrev != NULL)
	MvarPtFree(TJ -> TJuncPrev);
    if (TJ -> TJunc != NULL)
	MvarPtFree(TJ -> TJunc);
    if (TJ -> TJuncNext != NULL)
	MvarPtFree(TJ -> TJuncNext);
    IritFree(TJ);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deallocates and frees a list of  T-Junction structures.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   TJList:      Multivariate T-Junctions list to free.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroTJFree							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroTJFreeList                                                       M
*****************************************************************************/
void MvarZeroTJFreeList(MvarZeroTJunctionStruct *TJList)
{
    MvarZeroTJunctionStruct *Temp,
	*TJIter = TJList;

    while (TJIter != NULL) {
	Temp = TJIter -> Pnext;
	MvarZeroTJFree(TJIter);
	TJIter = Temp;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Handles the potential T-Junction events w.r.t a list of polylines.       *
* For each T-Junction candidate, traverses all polylines and finds the       *
* candidate's location. If needed, a copy of the candidate is inserted. If   *
* it is already there - it is tagged only.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Loops:	The solution lines.				             *
*   TJList:	The list of potential T-Junctions.			     *
*									     *
* RETURN VALUE:                                                              *
*   void	     							     *
*****************************************************************************/
void MvarZeroHandleTJunctions(MvarPolylineStruct *Loops, 
			      MvarZeroTJunctionStruct *TJList)
{
    int Dim = Loops -> Pl -> Dim;
    CagdRType 
	*MinDmn = (CagdRType *) IritMalloc(Dim * sizeof(CagdRType)),
	*MaxDmn = (CagdRType *) IritMalloc(Dim * sizeof(CagdRType));
    MvarZeroTJunctionStruct 
	*TJIter = TJList;
    MvarPolylineStruct *PolyIter;
    MvarZeroSubDmnInfoStruct *ProjInfo;

#if defined(DEBUG) && defined(DEBUG_TJ_LIST)
    MvarZeroSaveTJList(TJList);
    MvarZeroPrintTJList(TJList);
#endif /* DEBUG && DEBUG_TJ_LIST */

    while (TJIter != NULL) {
	PolyIter = Loops;
	while (PolyIter != NULL) {
	    ProjInfo = AttrGetRefPtrAttrib(PolyIter -> Attr, "_TrInfo");
	    MvarMVDomain(ProjInfo -> MVs[0], MinDmn, MaxDmn, -1);
	    if (MvarZeroIsPtOnDmnBndry(MinDmn, MaxDmn, TJIter -> TJunc)) {
		/* The potential T-Junction is indeed on the current         */
		/* domain's boundary. Search it and add if required:         */
		TJIter -> IsHandled = MvarZeroSearchTJInLoop(TJIter, PolyIter);
	    }
	    PolyIter = PolyIter -> Pnext;
	}
	TJIter = TJIter -> Pnext;
    }
    IritFree(MinDmn);
    IritFree(MaxDmn);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Searches a single polyline for a single T-Junction, assuming that if the *
* T-Junction is indeed there or needs to be added, then the previous and     *
* next points of the T-junction will be found first.                         *
*                                                                            *
* PARAMETERS:                                                                *
*   TJ:	    The T-Junction structure to search.                              *
*   Loop:   The polyline to search in.				             *
*									     *
* RETURN VALUE:                                                              *
*   CagdBType:	Success/Failure of the process.  	     		     *
*****************************************************************************/
static CagdBType MvarZeroSearchTJInLoop(const MvarZeroTJunctionStruct* TJ,
					MvarPolylineStruct *Loop)
{
    CagdBType IsLast, 
	TJHandled = FALSE,
	FoundNext = FALSE,
	FoundPrev = FALSE;
    MvarPtStruct
	*PtIter = Loop -> Pl;

    while (PtIter != NULL) {
	IsLast = PtIter -> Pnext == NULL;
	if (!(FoundPrev = MvarPtCmpTwoPoints(
	             TJ -> TJuncPrev, PtIter, MVAR_ZRMV2D_TJ_SEARCH_TOL) == 0))
	    FoundNext = MvarPtCmpTwoPoints(TJ -> TJuncNext, PtIter, 
	                                   MVAR_ZRMV2D_TJ_SEARCH_TOL) == 0;
	if (PtIter -> Pnext != NULL && (FoundPrev || FoundNext)) {
	    if (MvarPtCmpTwoPoints(TJ -> TJunc, PtIter -> Pnext, 
		                   MVAR_ZRMV2D_TJ_SEARCH_TOL) == 0) {
		/* Just tag it so it will not be purged, no adding needed. */
		AttrSetIntAttrib(&PtIter -> Pnext -> Attr, "TJunc", TRUE);
		TJHandled = TRUE;
		break;
	    }
	    else { /* Add the T-Junction as a new point. */
		MvarPtStruct
		    *NewPt = MvarPtCopy(TJ -> TJunc);

		assert(AttrGetIntAttrib(NewPt -> Attr, 
		                        "CrvSplitPt") == IP_ATTR_BAD_INT);
		AttrSetIntAttrib(&NewPt -> Attr, "TJunc", TRUE);
		if (!IsLast) {
		    NewPt -> Pnext = PtIter -> Pnext;
		    PtIter -> Pnext = NewPt;
		}
		else { /* Add as the last and as the first as well. */
		    MvarPtStruct
			*NewPtCpy = MvarPtCopy(NewPt);

		    PtIter -> Pnext = NewPt;
		    NewPt -> Pnext = NULL;
		    IRIT_LIST_PUSH(NewPtCpy, Loop -> Pl);
		    AttrSetIntAttrib(&NewPtCpy -> Attr, "TJunc", TRUE);
		}
		TJHandled = TRUE;
		break;
	    }
	}
	PtIter = PtIter -> Pnext;
    }
    return TJHandled;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test if a point belongs to the given domain's boundary, up to numeric    *
* tolerance.				                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   MinDmn:	    The min' values of the domain in each direction.	     *
*   MaxDmn:	    The max' values of the domain in each direction.	     *
*   Pt:             The  point to check.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType	    TRUE if point is on the boundary, FALSE otherwise.	     *
*****************************************************************************/
CagdBType MvarZeroIsPtOnDmnBndry(CagdRType *MinDmn,
				 CagdRType *MaxDmn,
				 const MvarPtStruct *Pt)
{
    CagdBType 
	FoundBoundaryDir = FALSE;
    int i, 
	Dim = Pt -> Dim;

    for (i = 0; i < Dim; i++) {
	if (IRIT_APX_EQ_EPS(Pt -> Pt[i], MinDmn[i],
	                    MVAR_ZRMV2D_PT_IN_FACE_TOL) || 
	    IRIT_APX_EQ_EPS(Pt -> Pt[i], MaxDmn[i],
			    MVAR_ZRMV2D_PT_IN_FACE_TOL)) {
	    FoundBoundaryDir = TRUE;
	    break;
	}
    }
    return FoundBoundaryDir;
}

#if defined(DEBUG) && defined (DEBUG_TJ_LIST)
/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates a file with the points of the global T-Junction list.            *
*                                                                            *
* PARAMETERS:                                                                *
*   void.				                                     *
*									     *
* RETURN VALUE:                                                              *
*   void.  	     		                                             *
*****************************************************************************/
static void MvarZeroSaveTJList(MvarZeroTJunctionStruct *TJList)
{
    int Handle = IPOpenDataFile("/temp/MVPts.itd", FALSE, TRUE);
    MvarZeroTJunctionStruct
	*TJIter = TJList;

    while (TJIter != NULL) {
	IPPutObjectToHandler(Handle, 
	                     IPGenPTObject(&(TJIter -> TJunc -> Pt[0]),
				           &(TJIter -> TJunc -> Pt[1]), 
				           &(TJIter -> TJunc -> Pt[2])));
	TJIter = TJIter -> Pnext;
    }
    IPCloseStream(Handle, TRUE);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Print a T-Junction list.                                                 *
*                                                                            *
* PARAMETERS:                                                                *
*   void.				                                     *
*									     *
* RETURN VALUE:                                                              *
*   void.  	     		                                             *
*****************************************************************************/
static void MvarZeroPrintTJList(MvarZeroTJunctionStruct *TJList)
{
    if (TJList == NULL)
	printf("List is empty");
    else {
	int i, 
	    TJCnt = 0,
	    Dim = TJList -> TJunc -> Dim;
	MvarZeroTJunctionStruct
	    *TJIter = TJList;

	while (TJIter != NULL) {
	    printf("TJ #%d: ", ++TJCnt);
	    for (i = 0; i < Dim; i++) {
		printf("%.6f ", TJIter -> TJunc -> Pt[i]);
	    }
	    printf("\n");
	    TJIter = TJIter -> Pnext;
	}
    }
}
#endif /* DEBUG && DEBUG_TJ_LIST */
