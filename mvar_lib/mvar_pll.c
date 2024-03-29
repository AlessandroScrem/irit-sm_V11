/******************************************************************************
* MVar_Pll.c - process polylines out of MvarPlnStruct.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Aug 2003.					      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "mvar_loc.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/allocate.h"

#define MVAR_MERGE_POLYLINES_REL_EPS 1e-6

#define TEST_CLOSEST_DIST(P1, P2, IsStart1, IsStart2) { \
		IrtRType DstSqr; \
		if ((DstSqr = MvarPtDistSqrTwoPoints(P1, P2)) < MinSqr) { \
	            MinSqr = DstSqr; \
		    *Start1 = IsStart1; \
		    *Start2 = IsStart2; \
	        } \
	    }

#define MV_TEST_CLOSEST_DIST(P1, P2) { \
		IrtRType DstSqr; \
		if ((DstSqr = MvarPtDistSqrTwoPoints(P1, P2)) < MinSqr) { \
	            MinSqr = DstSqr; \
	        } \
	    }

#define MV_SET_LAST_VERTEX(Pl) (Pl) -> PAux = MvarGetLastPt((Pl) -> Pl)
#define MV_SET_LAST_VERTEX2(Pl, PLast) (Pl) -> PAux = PLast
#define MV_GET_LAST_VERTEX(Pl) ((MvarPtStruct *) (Pl) -> PAux)

static void MVMergeGeomInitPolys(VoidPtr VPl);
static IrtRType MVMergeGeomDistSqr2Polys(VoidPtr VPl1, VoidPtr VPl2);
static IrtRType MVMergeGeomKeyPoly(VoidPtr VPl);
static int MVMergeGeomMerge2Polys(void **VPl1, void **VPl2);
static CagdRType MVDistanceBoundary(CagdRType *Pt,
				    int Dim,
				    IrtRType *ParamDomain);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Initialize the polylines to be merged.  Here, it amounts to caching the  *
* last vertex for fast access.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   VPl:    The polyline to initialize.                                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   void		                                                     *
*****************************************************************************/
static void MVMergeGeomInitPolys(VoidPtr VPl)
{
    MvarPolylineStruct
	*Pl = (MvarPolylineStruct *) VPl;

    MV_SET_LAST_VERTEX(Pl);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes minimal distance squared between end vertices of two polylines. *
*                                                                            *
* PARAMETERS:                                                                *
*   VPl1, VPl2:  To compute the minimal distance between the end points of.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtRType:   Minimal distance squared computed.			     *
*****************************************************************************/
static IrtRType MVMergeGeomDistSqr2Polys(VoidPtr VPl1, VoidPtr VPl2)
{
    IrtRType
	MinSqr = IRIT_INFNTY;
    MvarPolylineStruct
	*Pl1 = (MvarPolylineStruct *) VPl1,
	*Pl2 = (MvarPolylineStruct *) VPl2;
    MvarPtStruct
	*P1Start = Pl1 -> Pl,
	*P1End = MV_GET_LAST_VERTEX(Pl1),
        *P2Start = Pl2 -> Pl,
        *P2End = MV_GET_LAST_VERTEX(Pl2);

    MV_TEST_CLOSEST_DIST(P1Start, P2Start);
    MV_TEST_CLOSEST_DIST(P1End,   P2Start);
    MV_TEST_CLOSEST_DIST(P1Start, P2End);
    MV_TEST_CLOSEST_DIST(P1End,   P2End);

    return MinSqr;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returns a key to sort the polylines accordingly.			     *
*                                                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   VPl:  Polyline to return its key.                                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtRType:                                                                *
*****************************************************************************/
static IrtRType MVMergeGeomKeyPoly(VoidPtr VPl)
{
    MvarPolylineStruct
        *Pl = (MvarPolylineStruct *) VPl;

    /* Use the X coordinates of the first vertex as the key. */
    return Pl -> Pl -> Pt[0];
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Merges two polylines into one.  The address that references Pl2 is set   *
* to NULL after Pl2 is freed.  Pl1 will hold the merged poly.                *
*                                                                            *
* PARAMETERS:                                                                *
*   VPl1, VPl2:    The two polylines to merge into one.                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:     TRUE if successful, FALSE otherwise.                            *
*****************************************************************************/
static int MVMergeGeomMerge2Polys(void **VPl1, void **VPl2)
{
    MvarPolylineStruct
	**Pl1 = (MvarPolylineStruct **) VPl1,
	**Pl2 = (MvarPolylineStruct **) VPl2;
    MvarPtStruct
	*P1Start = (*Pl1) -> Pl,
	*P1End = MV_GET_LAST_VERTEX(*Pl1),
        *P2Start = (*Pl2) -> Pl,
        *P2End = MV_GET_LAST_VERTEX(*Pl2);
    IrtRType
	d00 = MvarPtDistSqrTwoPoints(P1Start, P2Start),
	d10 = MvarPtDistSqrTwoPoints(P1End,   P2Start),
	d01 = MvarPtDistSqrTwoPoints(P1Start, P2End),
	d11 = MvarPtDistSqrTwoPoints(P1End,   P2End);

    if (d00 < d10 && d00 < d01 && d00 < d11) {/* Merge Start Pl1 & Start Pl2.*/
        (*Pl1) -> Pl = P1Start = MvarPolyReverseList(P1End = P1Start);
    }
    else if (d10 < d01 && d10 < d11) {	       /* Merge End Pl1 & Start Pl2. */
    }
    else if (d01 < d11) {		       /* Merge Start Pl1 & End Pl2. */
        IRIT_SWAP(MvarPtStruct *, (*Pl1) -> Pl, (*Pl2) -> Pl);
	IRIT_SWAP(MvarPtStruct *, P1Start, P2Start);
	IRIT_SWAP(MvarPtStruct *, P1End, P2End);
    }
    else {				         /* Merge End Pl1 & End Pl2. */
        (*Pl2) -> Pl = P2Start = MvarPolyReverseList(P2End = P2Start);
    }

    /* Is merged point identical?  If so purge one of the two point. */
    if (!MvarPtCmpTwoPoints(P1End, P2Start, IRIT_UEPS)) {
	P1End -> Pnext = P2Start -> Pnext; 
	P2Start -> Pnext = NULL;
    }
    else {
	P1End -> Pnext = P2Start;
	(*Pl2) -> Pl = NULL;
    }
    MvarPolylineFree(*Pl2);
    *Pl2 = NULL;

    MV_SET_LAST_VERTEX2(*Pl1, P2End);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Merges separated multivariate polylines into longer ones, in place, as   M
* possible.								     M
*   Given a list of multivariate polylines, matches end points and merged as M
* possible multivariate polylines with common end points, in place.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Polys:      Multivariate polylines to merge, in place.                   M
*   Eps:	Epslion of similarity to merge multivariate polylines at.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *:  Merged as possible multivariate polylines.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMMergePolylines				                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarPolyMergePolylines, merge, polyline                                  M
*****************************************************************************/
MvarPolylineStruct *MvarPolyMergePolylines(MvarPolylineStruct *Polys,
					   IrtRType Eps)
{
    int i,
	NumOfPolys = CagdListLength(Polys);
    MvarPolylineStruct **PlsVec, *Pl;

    if (NumOfPolys < 2)
        return Polys;

    PlsVec = (MvarPolylineStruct **) IritMalloc(sizeof(MvarPolylineStruct *) *
					    NumOfPolys);
    for (i = 0, Pl = Polys; i < NumOfPolys; i++, Pl = Pl -> Pnext)
        PlsVec[i] = Pl;

    NumOfPolys = GMMergeGeometry((void **) PlsVec, NumOfPolys, Eps, IRIT_UEPS,
				 MVMergeGeomInitPolys,
				 MVMergeGeomDistSqr2Polys,
				 MVMergeGeomKeyPoly,
				 MVMergeGeomMerge2Polys);

    for (i = 1, Pl = Polys = PlsVec[0]; i < NumOfPolys; i++) {
        Pl -> Pnext = PlsVec[i];
	Pl = PlsVec[i];
    }
    Pl -> Pnext = NULL;

    IritFree(PlsVec);

    return Polys;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the last element in the list.                                    M
*                                                                            *
* PARAMETERS:                                                                M
*   Pts:    List of points to fetch last one.                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   Last element in list.                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarGetLastPt                                                            M
*****************************************************************************/
MvarPtStruct *MvarGetLastPt(MvarPtStruct *Pts)
{
    if (Pts == NULL)
	return NULL;

    while (Pts -> Pnext != NULL)
	Pts = Pts -> Pnext;

    return Pts;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A comparison function to examine if the given two points are the same.   M
*                                                                            *
* PARAMETERS:                                                                M
*   P1, P2:  Two multivariate points to compare.                             M
*   Eps:     The tolerance of the comparison.                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:      0 if identical, -1 or +1 if first point is less than/greater   M
*	      than second point, in lexicographic order over dimensions.     M
*             2 is returned if the dimensions are different.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarPtDistTwoPoints, MvarPtDistSqrTwoPoints                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarPtCmpTwoPoints                                                       M
*****************************************************************************/
int MvarPtCmpTwoPoints(const MvarPtStruct *P1,
		       const MvarPtStruct *P2,
		       CagdRType Eps)
{
    int i,
	Dim = P1 -> Dim;

    if (Dim != P2 -> Dim)
	return 2;

    for (i = 0; i < Dim; i++) {
	if (!IRIT_APX_EQ_EPS(P1 -> Pt[i], P2 -> Pt[i], Eps))
	    return IRIT_SIGN(P1 -> Pt[i] - P2 -> Pt[i]);
    }

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the Eucildean distance between two multivariate points.          M
*                                                                            *
* PARAMETERS:                                                                M
*   P1, P2:  Two points to compute the distance between.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:    Distance computed.                                         M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarPtCmpTwoPoints, MvarPtDistSqrTwoPoints                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarPtDistTwoPoints                                                      M
*****************************************************************************/
CagdRType MvarPtDistTwoPoints(const MvarPtStruct *P1, const MvarPtStruct *P2)
{
    return sqrt(MvarPtDistSqrTwoPoints(P1, P2));
}


/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes an in-between point (middle if t = 0.5) of given two points.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Pt1, Pt2:	Multivariate points.					     M
*   t:		Blending factor, 0 for Pt1, 0.5 for mid pt, 1 for Pt2.       M
*									     *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:  The in-between point of Pt1 and Pt2.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarPtCmpTwoPoints, MvarPtDistSqrTwoPoints                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarPtInBetweenPoint	                                             M
*****************************************************************************/
MvarPtStruct *MvarPtInBetweenPoint(const MvarPtStruct *Pt1,
				   const MvarPtStruct *Pt2,
				   CagdRType t)
{  
    int i,
	Dim = Pt1 -> Dim;
    CagdRType
        t1 = 1.0 - t;
    MvarPtStruct *InBetweenPt;

    if ((Pt1 == NULL || Pt2 == NULL) || (Pt1 -> Dim != Pt2 -> Dim))
	return NULL;

    InBetweenPt = MvarPtNew(Dim);

    for (i = 0; i < Dim; i++)
        InBetweenPt -> Pt[i] = Pt1 -> Pt[i] * t1 + Pt2 -> Pt[i] * t;	

    return InBetweenPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the Eucildean distance between two multivariate points.          M
*                                                                            *
* PARAMETERS:                                                                M
*   P1, P2:  Two points to compute the distance between.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:    Distance computed.                                         M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarPtCmpTwoPoints, MvarPtDistTwoPoints                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarPtDistSqrTwoPoints                                                   M
*****************************************************************************/
CagdRType MvarPtDistSqrTwoPoints(const MvarPtStruct *P1, const MvarPtStruct *P2)
{
    int i,
	Dim = P1 -> Dim;
    CagdRType
	Dist = 0.0;

    if (Dim != P2 -> Dim)
	return 0.0;

    for (i = 0; i < Dim; i++) {
	Dist += IRIT_SQR(P1 -> Pt[i] - P2 -> Pt[i]);
    }

    return Dist;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Connect the list of multivariate points into multivariate polylines by   M
* connecting the closest multivariate point pairs, until the distances       M
* between adjacent multivariate points/polylines is more than MaxTol.        M
*                                                                            *
* PARAMETERS:                                                                M
*   PtsList:      Point list to connect into multivariate polylines.         M
*   MaxTol:       Maximum distance allowed to connect multivariate points.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *:    Connected multivariate polylines, upto MaxTol   M
*			     tolerance.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMMatchPointListIntoPolylines2		                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMatchPointListIntoPolylines                                          M
*****************************************************************************/
MvarPolylineStruct *MvarMatchPointListIntoPolylines(const MvarPtStruct *PtsList,
						    IrtRType MaxTol)
{
    const MvarPtStruct *Pt;
    MvarPolylineStruct *Pl1, *PllList;

    /* Convert the poly list to polyline linked list with one point in   */
    /* each polyline so we can start and match-connect them.		 */
    for (PllList = NULL, Pt = PtsList; Pt != NULL; Pt = Pt -> Pnext) {
        Pl1 = MvarPolylineNew(MvarPtCopy(Pt));
	IRIT_LIST_PUSH(Pl1, PllList);
    }

    return MvarPolyMergePolylines(PllList, MaxTol);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of multivariate points into a list of cagd E3 points.    M
* Assumes all input points are E1, E2, or E3.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPts:     List of multivariate points to convert to list of =points.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:   A list of points.                                      M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPtsToPolys                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPtsToPts                                                      M
*****************************************************************************/
CagdPtStruct *MvarCnvrtMVPtsToPts(const MvarPtStruct *MVPts)
{
    const MvarPtStruct *MVPt;
    CagdPtStruct *Pt,
        *Pts = NULL;

    for (MVPt = MVPts; MVPt != NULL; MVPt = MVPt -> Pnext) {
        assert(MVPt -> Dim <= 3);
	Pt = CagdPtNew();
        IRIT_PT_COPY(Pt -> Pt, MVPt -> Pt);
	IRIT_LIST_PUSH(Pt, Pts);
    }

    return CagdListReverse(Pts);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of multivariate points into a list of control points.    M
* Assumes all points of same dimensions.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPts:     List of multivariate points to convert to list of ctlpts.     M
*   MergeTol:  If non negative, attempt to merge the data into polylines.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A list object of ctlpts.                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToPolys, MvarCnvrtMVPtsToPts                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPtsToPts                              M
*****************************************************************************/
IPObjectStruct *MvarCnvrtMVPtsToCtlPts(const MvarPtStruct *MVPts,
				       IrtRType MergeTol)
{
    int i, k;
    MvarPtStruct *MVPt;
    MvarPolylineStruct *MVPls, *MVPl;
    CagdPointType PType;
    IPObjectStruct *PObj,
	*ObjPlList = IPGenLISTObject(NULL);

    if (MVPts == NULL)
        return ObjPlList;

    PType = CAGD_MAKE_PT_TYPE(FALSE, IRIT_MIN(MVPts -> Dim,
			      CAGD_MAX_PT_COORD));

    if (MergeTol > 0)
        MVPls = MvarMatchPointListIntoPolylines(MVPts, MergeTol);
    else
	MVPls = MvarPolylineNew(MvarPtCopyList(MVPts));

    for (MVPl = MVPls, k = 0; MVPl != NULL; MVPl = MVPl -> Pnext) {
        IPObjectStruct
	    *ObjPtList = IPGenLISTObject(NULL);

	for (MVPt = MVPl -> Pl, i = 0; MVPt != NULL; MVPt = MVPt -> Pnext) {
	    IPListObjectInsert(ObjPtList, i++,
			       PObj = IPGenCTLPTObject(PType, &MVPt -> Pt[-1]));
	    IP_ATTR_SAFECOPY_ATTRS(PObj -> Attr, MVPt -> Attr);
	}
	IPListObjectInsert(ObjPtList, i, NULL);

	if (MVPl == MVPls && MVPls -> Pnext == NULL) {
	    /* Only one polyline - return it as a list of point obj. */
	    MvarPolylineFreeList(MVPls);
	    IPFreeObject(ObjPlList);
	    return ObjPtList;
	}
	else
	    IPListObjectInsert(ObjPlList, k++, ObjPtList);
    }
    IPListObjectInsert(ObjPlList, k, NULL);

    MvarPolylineFreeList(MVPls);

    return ObjPlList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of multivariate points into a list of polylines.         M
* Assumes all points of same dimensions.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPts:     List of multivariate points to convert to polylines.          M
*   MV:        A multivariate to evaluate through, if not NULL.		     M
*   MergeTol:  Tolerance to merge points into polylines.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A list object of polylines.                          M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPtsToPolys2                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPtsToPolys                                                    M
*****************************************************************************/
IPObjectStruct *MvarCnvrtMVPtsToPolys(const MvarPtStruct *MVPts,
				      const MvarMVStruct *MV,
				      IrtRType MergeTol)
{
    int i, k;
    MvarPtStruct *MVPt;
    MvarPolylineStruct *MVPls, *MVPl;
    CagdPointType PType;
    IPObjectStruct
	*RetList = IPGenLISTObject(NULL);

    if (MVPts == NULL)
	return RetList;

    if (MV != NULL)
        PType = (CagdPointType) MV -> PType;
    else
        PType = CAGD_MAKE_PT_TYPE(FALSE, IRIT_MIN(MVPts -> Dim,
				  CAGD_MAX_PT_COORD));

    MVPls = MvarMatchPointListIntoPolylines(MVPts, MergeTol);

    for (MVPl = MVPls, k = 0; MVPl != NULL; MVPl = MVPl -> Pnext) {
        IPObjectStruct *CtlPtObj,
	    *CtlPtList = IPGenLISTObject(NULL);

	for (MVPt = MVPl -> Pl, i = 0; MVPt != NULL; MVPt = MVPt -> Pnext) {
	    IrtRType *R;

	    if (MV != NULL)
	        R = MvarMVEval(MV, MVPt -> Pt);
	    else
	        R = &MVPt -> Pt[-1];

	    CtlPtObj = IPGenCTLPTObject(PType, R);
	    IPListObjectInsert(CtlPtList, i++, CtlPtObj);

	    IP_ATTR_SAFECOPY_ATTRS(CtlPtObj -> Attr, MVPt -> Attr);
	}
	IPListObjectInsert(CtlPtList, i, NULL);

	IPListObjectInsert(RetList, k++, CtlPtList);
    }
    IPListObjectInsert(RetList, k, NULL);

    MvarPolylineFreeList(MVPls);

    return RetList;
}

/*****************************************************************************
* DESCRIPTION:							 	     M
*   Connect a list of discrete points into a polylines.		 	     M
*									     M
* PARAMETERS:							 	     M
*   InPts:        A list of discrete points.				     M
*   FineNess:     Tolerance.				 		     M
*   Dim:          The dimension of discrete points, 1 to 3.	 	     M
*   ParamDomain:  The domain of the mvar points.                             M
*									     M
* RETURN VALUE:							 	     M
*   IPPolygonStruct *: Connected list of polylines.		 	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPtsToPolys                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPtsToPolys2                                                   M
*****************************************************************************/
IPPolygonStruct *MvarCnvrtMVPtsToPolys2(const MvarPtStruct *InPts,
					CagdRType FineNess,
					int Dim,
					IrtRType *ParamDomain)
{
    IPPolygonStruct *RetVal, RetVals,
	*PHead = NULL;
    IPVertexStruct 
	*VHead = NULL;
    MvarPtStruct *Pt, *Marker, *ConPts, *ConPt, *ConPt2,
	*Pts = MvarPtCopyList(InPts);

    RetVal = &RetVals;

    /* Connect the solution points into polylines. */
    for (Pt = Pts; Pt != NULL; ) {
	IrtRType 
	    Small = IRIT_INFNTY;
	MvarPtStruct 
	    *TargetPt = NULL;
	int Length = 0;
	    
	if (Pt -> Dim == -1) { 
	    /* We already visited this point. */
	    Pt = Pt -> Pnext;
	    continue;
	}

	/* We should find the point which is closest to the boundary and    */
	/* it becomes the first solution.			            */
	for (ConPt = Pt; ConPt != NULL; ConPt = ConPt -> Pnext) {
	    IrtRType 
		NearBoundary = 0.0;
			
	    if (ConPt -> Dim != -1) {
		NearBoundary = MVDistanceBoundary(ConPt -> Pt, Dim,
						  ParamDomain);

		if (Small > NearBoundary) {
		    Small = NearBoundary;
		    TargetPt = ConPt;
		}
	    }
	}

	if (TargetPt != NULL) {
	    ConPt = ConPts = MvarPtCopy(TargetPt);
	    TargetPt -> Dim = -1;
	    Length = 1;
	    Marker = Pt;
	}
	else {
#ifdef DEBUG_MVAR_OFFSET_NI_TRACING
	    fprintf(stderr, "This should not be happen\n");
#endif /* DEBUG_MVAR_OFFSET_NI_TRACING */
	    break;
	}
		
	while (TRUE) {
	    int i;
	    
	    Small = IRIT_INFNTY; 
	    TargetPt = NULL;
			
	    for (ConPt2 = Pt; ConPt2 != NULL; ConPt2 = ConPt2 -> Pnext) {
		if (ConPt2 -> Dim != -1) {
		    IrtRType
			Distance = 0.0;
		    
		    for (i = 0; i < Dim; i++) {
			Distance += (ConPt2 -> Pt[i] - ConPt -> Pt[i]) *
			            (ConPt2 -> Pt[i] - ConPt -> Pt[i]);
		    }

		    if (Distance < Small) {
			Small = Distance;
			TargetPt = ConPt2;
		    }
		}
	    }

	    /* If distance is larger than grid size, then disconnect this. */
	    if (Small > 3.0 * IRIT_SQR(FineNess) || TargetPt == NULL)
		break;
	    ConPt -> Pnext = MvarPtCopy(TargetPt);
	    Length++;

	    ConPt = ConPt -> Pnext;
	    ConPt -> Pnext = NULL;
	    TargetPt -> Dim = -1;
	}
	    
	if (Length > 100) {
	    VHead = NULL;
	    
	    /* Connect ConPts into a curve. */
	    for (ConPt = ConPts; ConPt != NULL; ConPt = ConPt -> Pnext) {
		VHead = IPAllocVertex2(VHead);
		IRIT_PT_SET(VHead -> Coord, ConPt -> Pt[0], ConPt -> Pt[1],
						       ConPt -> Pt[2]);
	    }
	    PHead = IPAllocPolygon(0, VHead, PHead);
	    VHead = NULL;
	}
	if (PHead != NULL) {
	    RetVal -> Pnext = PHead;
	    PHead = NULL;
	    RetVal = RetVal -> Pnext;
	}


	/* Free current stuffs. */
	MvarPtFreeList(ConPts);
	Pt = Marker;
    }
    RetVal -> Pnext = NULL;

    /* Free local stuffs. */
    MvarPtFreeList(Pts);
    
    return GMMergePolylines(RetVals.Pnext, FineNess * 2);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of MV polylines to a list of MV points.                  M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPlls:   MV polylines to convert to MV points.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:  MV points.                                              M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPolysToIritPolys                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPolysToMVPts                                                  M
*****************************************************************************/
MvarPtStruct *MvarCnvrtMVPolysToMVPts(const MvarPolylineStruct *MVPlls)
{
    const MvarPolylineStruct *MVPl;
    MvarPtStruct
        *MVPts = NULL;

    for (MVPl = MVPlls; MVPl != NULL; MVPl = MVPl -> Pnext) {
        MVPts = (MvarPtStruct *) CagdListAppend(MvarPtCopyList(MVPl -> Pl),
						MVPts);
    }

    return MVPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of MV polylines to a list of MV control points.          M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPlls:   MV polylines to convert to MV points.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   struct IPObjectStruct *:  MV control points.                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPolysToIritPolys                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPolysToCtlPts                                                 M
*****************************************************************************/
struct IPObjectStruct *MvarCnvrtMVPolysToCtlPts(const MvarPolylineStruct
								     *MVPlls)
{
    int i = 0;
    const MvarPolylineStruct *MVPl;
    IPObjectStruct
	*ObjPlList = IPGenLISTObject(NULL);

    for (MVPl = MVPlls; MVPl != NULL; MVPl = MVPl -> Pnext) {
        IPObjectStruct
 	    *ObjPl = MvarCnvrtMVPtsToCtlPts(MVPl -> Pl, -1.0);

	IPListObjectInsert(ObjPlList, i++, ObjPl);
    }
    IPListObjectInsert(ObjPlList, i, NULL);

    return ObjPlList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of mvar polylines into a list of irit polylines.         M
* Assumes all points of same dimensions.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPlls:     List of multivariate polylines to convert to Irit polylines. M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A list object of polylines (as lists of ctlpts).     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPtsToPolys2, MvarCnvrtMVPtsToPolys2   M
*   MvarCnvrtMVPolysToIritCrvs						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPolysToIritPolys                                              M
*****************************************************************************/
IPObjectStruct *MvarCnvrtMVPolysToIritPolys(const MvarPolylineStruct *MVPlls)
{
    const MvarPolylineStruct *MVPl;
    IPObjectStruct
        *AllObj = IPAllocObject("", IP_OBJ_LIST_OBJ, NULL);

    for (MVPl = MVPlls; MVPl != NULL; MVPl = MVPl -> Pnext) {
        MvarPtStruct *MVPt,
	    *MVPts = MVPl -> Pl;

	if (MVPts == NULL) {
	    /* An empty list. */
	    continue;
	}
	else if (MVPts -> Pnext == NULL) {
	    /* A single point in plls.  Return a single point object. */
	    IPListObjectAppend(AllObj,
		IPGenCTLPTObject(CAGD_MAKE_PT_TYPE(FALSE, MVPts -> Dim),
				 &MVPts -> Pt[-1]));
	}
	else {
	    int i = 0;
	    IPObjectStruct
		*CtlPtsObj = IPAllocObject("", IP_OBJ_LIST_OBJ, NULL);

	    for (MVPt = MVPts; MVPt != NULL; MVPt = MVPt -> Pnext) {
	        IPListObjectInsert(CtlPtsObj, i++,
		    IPGenCTLPTObject(CAGD_MAKE_PT_TYPE(FALSE, MVPt -> Dim),
				     &MVPt -> Pt[-1]));
	    }

	    IPListObjectInsert(CtlPtsObj, i, NULL);
	    IPListObjectAppend(AllObj, CtlPtsObj);
	}
    }

    IPListObjectAppend(AllObj, NULL);

    return AllObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of mvar polylines into a list of irit polylines.         M
* Assumes all points of same dimensions.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPlls:     List of multivariate polylines to convert to Irit polylines. M
*   IgnoreIndividualPts:  True to handle only polylines, ignoring points.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A polylines object, or if have individual points a   M
*                       list of two objects, as (PllnObjs, PntObjs).         M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPtsToPolys2, MvarCnvrtMVPtsToPolys    M
*   MvarCnvrtMVPolysToIritCrvs						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPolysToIritPolys2                                             M
*****************************************************************************/
IPObjectStruct *MvarCnvrtMVPolysToIritPolys2(const MvarPolylineStruct *MVPlls,
					     int IgnoreIndividualPts)
{
    const MvarPolylineStruct *MVPl;
    const CagdRType
        Z = 0.0;
    int Dim;
    IPObjectStruct *AllObjs,
        *PtsObj = NULL,
        *PllsObj = NULL;

    if (MVPlls == NULL)
	return NULL;
    Dim = MVPlls -> Pl -> Dim;
    for (MVPl = MVPlls; MVPl != NULL; MVPl = MVPl -> Pnext) {
        MvarPtStruct *MVPt,
	    *MVPts = MVPl -> Pl;

	if (MVPts == NULL) {
	    /* An empty list. */
	    continue;
	}
	else if (MVPts -> Pnext != NULL) {
	    IPPolygonStruct
		*Plln = IPAllocPolygon(0, NULL, NULL);

	    if (PllsObj == NULL) 
	        PllsObj = IPGenPOLYLINEObject(NULL);

	    Plln -> Pnext = PllsObj -> U.Pl;
	    PllsObj -> U.Pl = Plln;
	    for (MVPt = MVPts; MVPt != NULL; MVPt = MVPt -> Pnext) {
	        IPVertexStruct
		    *V = IPAllocVertex2(Plln -> PVertex);

		Plln -> PVertex = V;

		V -> Coord[0] = MVPt -> Pt[0];
		V -> Coord[1] = Dim > 1 ? MVPt -> Pt[1] : Z;
		V -> Coord[2] = Dim > 2 ? MVPt -> Pt[2] : Z;
	    }
	}
	else if (!IgnoreIndividualPts) {
	    /* A single point in plls.  Return a single point object. */
	    if (PtsObj == NULL) 
	        PtsObj = IPAllocObject("", IP_OBJ_LIST_OBJ, NULL);

	    IPListObjectAppend(PtsObj,
			       IPGenPTObject(&MVPts -> Pt[0],
					     Dim > 1 ? &MVPts -> Pt[1] : &Z,
					     Dim > 2 ? &MVPts -> Pt[2] : &Z));
	}
    }

    if (PtsObj != NULL && PllsObj != NULL) {
        AllObjs = IPAllocObject("", IP_OBJ_LIST_OBJ, NULL);
	IPListObjectAppend(AllObjs, PtsObj);
	IPListObjectAppend(AllObjs, PllsObj);
    }
    else
        AllObjs = PtsObj != NULL ? PtsObj : PllsObj;

    return AllObjs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of mvar polylines into a list of irit curves.            M
* Assumes all points of same dimensions.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPlls:     List of multivariate polylines to convert to Irit polylines. M
*   Order:      Order of constructed curves, typically 2 (linear).           M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   A list of curves.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVPtsToCtlPts, MvarCnvrtMVPolysToIritPolys2		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVPolysToIritCrvs                                               M
*****************************************************************************/
CagdCrvStruct *MvarCnvrtMVPolysToIritCrvs(const MvarPolylineStruct *MVPlls,
					  int Order)
{
    const MvarPolylineStruct *MVPl;
    CagdCrvStruct *Crv,
        *Crvs = NULL;
    
    for (MVPl = MVPlls; MVPl != NULL; MVPl = MVPl -> Pnext) {
        MvarPtStruct *MVPt,
	    *MVPts = MVPl -> Pl;

	if (MVPts == NULL || MVPts -> Pnext == NULL) {
	    /* An empty list. */
	    continue;
	}
	else {
	    int i, j,
		Dim = MVPts -> Dim,
		Len = CagdListLength(MVPts);

	    Crv = BspCrvNew(Len, Order > Len ? Len : Order,
			    CAGD_MAKE_PT_TYPE(FALSE, Dim));
	    BspKnotUniformOpen(Len, Order > Len ? Len : Order,
			       Crv -> KnotVector);

	    for (MVPt = MVPts, i = 0; MVPt != NULL; MVPt = MVPt -> Pnext, i++) {
		for (j = 0; j < Dim; j++)
		    Crv -> Points[j + 1][i] = MVPt -> Pt[j];
	    }
	    IRIT_LIST_PUSH(Crv, Crvs);
	}
    }

    Crvs = CagdListReverse(Crvs);

    return Crvs;
}

/*****************************************************************************
* DESCRIPTION:							 	     *
*   Computes distance from the point to the boundary.		 	     *
*									     *
* PARAMETERS:							 	     *
*   Pt:   Point to compute the distance.	    		 	     *
*   Dim:  The dimension of the point.					     *
*   ParamDomain: The domain of the mvar points.                              *
*									     *
* RETURN VALUE:							 	     *
*   Distance to the boundary.						     *
*****************************************************************************/
static CagdRType MVDistanceBoundary(CagdRType *Pt,
				    int Dim,
				    IrtRType *ParamDomain)
{
    IrtRType 
	L = IRIT_INFNTY;
    int i;

    for (i = 0; i < Dim; i++) {
	IrtRType
	    Li = IRIT_MIN(Pt[i] - ParamDomain[i * 2],
		     ParamDomain[i * 2 + 1] - Pt[i]);

	L = IRIT_MIN(L, Li);
    }
	
    return L;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of mvar triangles into a list of irit polygons. If the   M
* list of triangles consists of points of dimension higher than three,       M
* the polygons list is created as the projection on R^3, using the           M
* coordinates specified in Coords.				             M
* Assumes all points of same dimensions.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVTrs:     List of multivariate triangles to convert to Irit polygons.   M
*   Coords:    The required coordinates, or NULL if the dimension is three.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A list object of polygons.		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarIrit2DTrTo2DMVTrs, MvarCnvrtMVPolysToIritPolys			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCnvrtMVTrsToIritPolygons                                             M
*****************************************************************************/
IPObjectStruct *MvarCnvrtMVTrsToIritPolygons(const MvarTriangleStruct *MVTrs,
					     int *Coords)
{
    int i, HighDimIter, LowDimIter,
	Dim = MVTrs -> Dim,
        DoCirc = IPSetPolyListCirc(FALSE);
    const MvarTriangleStruct *MVTrIter;
    IPVertexStruct *NewVertex, *CurVertexList;
    IPPolygonStruct *NewPlgn;
    IPObjectStruct
	*AllObj = IPAllocObject("", IP_OBJ_POLY, NULL);

    IPSetPolyListCirc(DoCirc);

    for (MVTrIter = MVTrs; MVTrIter != NULL; MVTrIter = MVTrIter -> Pnext) {
#	ifdef MVAR_ZRMV2D_DBG_TR_CNVRT_PL
	    MvarZeroDbgTriangle(MVTrIter);
#	endif /* MVAR_ZRMV2D_DBG_TR_CNVRT_PL */

	/* Create the vertex list of the currently handled triangle. */
	CurVertexList = NULL;
	for (i = 0; i < 3; i++) {
	    NewVertex = IPAllocVertex(FALSE, NULL, NULL);
	    if (Coords == NULL) { /* MVTrs is already 3D, or simply first 3. */
		IRIT_GEN_COPY(NewVertex -> Coord, MVTrIter -> Vrtcs[i],
			      3 * sizeof(CagdRType));
		if (Dim == 3) {
		    IRIT_GEN_COPY(NewVertex -> Normal, MVTrIter -> Nrmls[i], 
		                  3 * sizeof(CagdRType));
		    if (!IRIT_PT_EQ_ZERO(NewVertex -> Normal))
			IP_SET_NORMAL_VRTX(NewVertex);
		}
	    }
	    else {   /* Need to project only the selected three coordinates. */
		HighDimIter = 0;
		for (LowDimIter = 0; LowDimIter < 3; LowDimIter++) {
		    while (HighDimIter < i)
			HighDimIter++;
		    /* Now the coordinate can be copied. */
		    NewVertex -> Coord[LowDimIter] = 
			MVTrIter -> Vrtcs[i][HighDimIter];
		}
	    }
	    
	    IRIT_LIST_PUSH(NewVertex, CurVertexList);
	}

	/* The vertex list should be cyclic: */
	if (DoCirc) {
	    IPGetLastVrtx(CurVertexList) -> Pnext = CurVertexList;
	}

	/* Each triangle corresponds to an IPPolygon: */
	NewPlgn = IPAllocPolygon(FALSE, CurVertexList, NULL);
	IRIT_LIST_PUSH(NewPlgn, AllObj -> U.Pl);
    }
    return AllObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of irit planar triangles (given as polygon objects)      M
* into a list of Mvar planar (2D) triangle structures.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   ObjTrs:    Irit polygons.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarTriangleStruct *:   A list triangles.		                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCnvrtMVTrsToIritPolygons, MvarCnvrtMVPolysToIritPolys		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarIrit2DTrTo2DMVTrs                                                    M
*****************************************************************************/
MvarTriangleStruct *MvarIrit2DTrTo2DMVTrs(IPObjectStruct *ObjTrs)
{
    int i;
    MvarTriangleStruct *CurTr,
	*MVTrs = NULL;
    IPVertexStruct *CurVertex;
    IPPolygonStruct
	*IritTrsIter = ObjTrs -> U.Pl;

    while (IritTrsIter != NULL) {
	CurTr = MvarTriangleNew(2);
	CurVertex = IritTrsIter -> PVertex;
	for (i = 0; i < 3; i++) {
	    IRIT_GEN_COPY(CurTr -> Vrtcs[i], CurVertex -> Coord, 
			  2 * sizeof(CagdRType));
	    IRIT_ZAP_MEM(CurTr -> Nrmls[i], 2 * sizeof(CagdRType));
	    CurVertex = CurVertex -> Pnext;
	}
	IRIT_LIST_PUSH(CurTr, MVTrs);
	IritTrsIter = IritTrsIter -> Pnext;
    }
    return MVTrs;
}


