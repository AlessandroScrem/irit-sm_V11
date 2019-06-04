/******************************************************************************
* Poly_Pts.c - polygonal data point and other filtering tools.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Feb 1997.					      *
******************************************************************************/

#include <math.h>
#include <stdio.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "geom_loc.h"
#include "inc_irit/bool_lib.h"
#include "inc_irit/ip_cnvrt.h"

#define GM_MERGE_POLYLINES_REL_EPS 1e-6

#define VERTEX_COPY(VDest, VSrc)  { IRIT_PT_COPY(VDest -> Coord, VSrc -> Coord); \
				    IRIT_PT_COPY(VDest -> Normal, VSrc -> Normal); \
				    VDest -> Tags = VSrc -> Tags; \
				    VDest -> Attr = \
				        IP_ATTR_COPY_ATTRS(VSrc -> Attr); }

#define COPY_RGB_ATTR(VDest, VSrc) { int RTmp, GTmp, BTmp; \
		if (AttrGetRGBColor(VSrc -> Attr, &RTmp, &GTmp, &BTmp)) \
		    AttrSetRGBColor(&VDest -> Attr, RTmp, GTmp, BTmp); }

#define COPY_UV_ATTR(VDest, VSrc) { \
    float *UV; \
    IRIT_STATIC_DATA AttribNumType \
        AttrUVValsID = ATTRIB_NAME_BAD_NUMBER; \
    IP_ATTR_INIT_UNIQUE_ID_NUM(AttrUVValsID, "uvvals"); \
    if ((UV = AttrGetUVAttrib2(VSrc -> Attr, AttrUVValsID)) != NULL) \
        AttrSetUVAttrib2(&VDest -> Attr, AttrUVValsID, UV[0], UV[1]); }

#define TEST_CLOSEST_DIST(V1, V2, IsStart1, IsStart2) { \
		IrtRType DstSqr; \
		if ((DstSqr = IRIT_PT_PT_DIST_SQR(V1 -> Coord, \
					     V2 -> Coord)) < MinSqr) { \
	            MinSqr = DstSqr; \
		    *Start1 = IsStart1; \
		    *Start2 = IsStart2; \
	        } \
	    }

typedef struct GMPlPlDistStruct {
    IPPolygonStruct *Pl;
    int Idx, MinIdx, MinStart1, MinStart2;
    IrtRType MinDistSqr;
} GMPlPlDistStruct;

static IrtRType GetPlPlDist(IPPolygonStruct *Pl1,
			    IPPolygonStruct *Pl2,
			    int *Start1,
			    int *Start2);
#if defined(ultrix) && defined(mips)
static int ComparePlPlDist(VoidPtr PlPlDist1, VoidPtr PlPlDist2);
#else
static int ComparePlPlDist(const VoidPtr PReal1, const VoidPtr PReal2);
#endif /* ultrix && mips (no const support) */
static void UpdateAdjPolys(IPVertexStruct *V1,
			   IPPolygonStruct *Pl1,
			   IPVertexStruct *V2,
			   IPPolygonStruct *Pl2);
static void InsertVertexToSplitList(IPVertexStruct *V,
				    IrtRType *Pos,
				    IrtRType t);
static IPVertexStruct *GetAdjVertex(IPVertexStruct *V,
				    IPPolygonStruct *Pl,
				    IPPolygonStruct *PAdj);
static IPVertexStruct *GMFindthirdVertexinTJun(const IPVertexStruct *V,
					       const IPPolygonStruct *Pl,
					       const IPVertexStruct *VAdj,
					       const IPPolygonStruct *PAdj,
					       IPPolygonStruct *Polys,
					       IPPolygonStruct **Pl2Ret,
					       IrtRType Eps);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new polygonal objects out of given one, that contains only     M
* polygons of upto n vertices.  Non convex polygons are split to convex one  M
* so the result will contain convex data only.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   PolyObj:   Polygonal object to split into up to n-gons.                  M
*   n:	       Maximal number of vertices.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A polygonal object containing polygons with upto n   M
*		vertices, representing the same model as PolyObj.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   ConvexPolyObjectN, GMConvertPolysToTriangles, GMConvertPolysToRectangles M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMConvertPolysToNGons                                                    M
*****************************************************************************/
IPObjectStruct *GMConvertPolysToNGons(IPObjectStruct *PolyObj, int n)
{
    IPPolygonStruct *Pl;
    int IsCirc = IPSetPolyListCirc(FALSE);

    IPSetPolyListCirc(IsCirc);	      /* Restore state, now that we know it. */

    n = IRIT_MAX(n, 3); /* No less than a triangle! */

    PolyObj = GMConvexPolyObjectN(PolyObj);

    for (Pl = PolyObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	IPVertexStruct *V,
	    *VHead = Pl -> PVertex;
	int j, m,
	    Len = IPVrtxListLen(VHead);

	for (j = 3, V = VHead; j < Len; j++, V = V -> Pnext) {
	    if (!GMCoplanar4Pts(V -> Coord, V -> Pnext -> Coord,
				V -> Pnext -> Pnext -> Coord,
				V -> Pnext -> Pnext -> Pnext -> Coord))
	    break;
	}
	if (j < Len)        /* Non planar data - split a triangle out of it. */
	    m = 3;
	else
	    m = n;

	if (Len > m) {
	    int i;
	    IPVertexStruct *VRest, *VPrev,
		*VLast = IPGetLastVrtx(VHead);
	    IPPolygonStruct *PlNew;

	    /* Find the splitting point. */
	    for (i = 1, VRest = VHead, VPrev = NULL;
		 i < m;
		 i++, VPrev = VRest, VRest = VRest -> Pnext);

	    /* Isolate the first polygon out of the original polygon. */
	    VPrev -> Pnext = IPAllocVertex2(IsCirc ? VHead : NULL);
	    VERTEX_COPY(VPrev -> Pnext, VRest);
	    IP_SET_INTERNAL_VRTX(VPrev -> Pnext);

	    /* Put the rest in a new polygon to be treated soon. */
	    PlNew = IPAllocPolygon(0, VRest, Pl -> Pnext);
	    IRIT_PLANE_COPY(PlNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlNew);
	    Pl -> Pnext = PlNew;
	    PlNew -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr);

	    VLast -> Pnext = IPAllocVertex2(IsCirc ? VRest : NULL);
	    VERTEX_COPY(VLast -> Pnext, VHead);
	    IP_SET_INTERNAL_VRTX(VLast -> Pnext);
	}
    }

    return PolyObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new polygonal objects out of given one, that contains only     M
* triangles.  Non convex polygons are split to convex one which, in turn,    M
* converted to triangles.  Collinear points are purged away.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   PolyObj:   Polygonal object to split into triangles.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A polygonal object containing only triangles         M
*		representing the same model as PolyObj.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   ConvexPolyObjectN, GMConvertPolysToNGons, GMLimitTrianglesEdgeLen        M
*   GMConvertPolysToTriangles2, GMConvertPolysToRectangles		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMConvertPolysToTriangles                                                M
*****************************************************************************/
IPObjectStruct *GMConvertPolysToTriangles(IPObjectStruct *PolyObj)
{
    IPPolygonStruct *Pl;
    int IsCirc = IPSetPolyListCirc(FALSE);

    IPSetPolyListCirc(IsCirc);	      /* Restore state, now that we know it. */

    PolyObj = GMConvexPolyObjectN(PolyObj);

    for (Pl = PolyObj -> U.Pl; Pl != NULL; ) {
	IPPolygonStruct
	    *PlNext = Pl -> Pnext;
	IPVertexStruct
	    *VHead = Pl -> PVertex;
	int Len = IPVrtxListLen(VHead),
	    FreeVHead = FALSE;

	if (Len > 3) {	                   /* Split into several triangles. */
	    int VLastTags,
                FirstIsDegenerated = FALSE;
	    IPVertexStruct *VLast,
		*VRest = VHead -> Pnext -> Pnext -> Pnext;
	    IPPolygonStruct
		*PlNew = NULL;

            /* Remove degenerate attribute. The new triangles will be marked */
            /* 'degenerated'eperately if needed.                             */
            if (IPGlblGenDegenPolys)
                AttrFreeOneAttribute(&Pl -> Attr, IP_ATTRIB_DEGEN_POLY);

	    /* Isolate the first triangle out of the original polygon. */
	    VHead -> Pnext -> Pnext -> Pnext = IsCirc ? VHead : NULL;
	    VLast = VHead -> Pnext -> Pnext;
	    VLastTags = VLast -> Tags;
	    IP_SET_INTERNAL_VRTX(VLast);
	    if (GMCollinear3Pts(VHead -> Coord,
				VHead -> Pnext -> Coord,
				VHead -> Pnext -> Pnext -> Coord)) {
                if (IPGlblGenDegenPolys) {
                    FirstIsDegenerated = TRUE;
                }
                else {
	            FreeVHead = TRUE;
		    Pl -> PVertex = NULL;
                }
	    }

	    /* Construct triangles out of the rest of the points. */
	    while (VRest != NULL && VRest != VHead) {
		IPVertexStruct
		    *VRestNext = VRest -> Pnext,
		    *V3 = IPAllocVertex2(NULL),
		    *V2 = IPAllocVertex2(V3),
		    *V1 = IPAllocVertex2(V2);
                int Collinear,
		    KeepDegenerated = FALSE;

		VERTEX_COPY(V1, VHead);
		VERTEX_COPY(V2, VLast);
		VERTEX_COPY(V3, VRest);

		if (IsCirc)
		    V3 -> Pnext = V1;

		IP_SET_INTERNAL_VRTX(V1);
		V2 -> Tags = VLastTags;
		if (VRest -> Pnext != NULL && VRest -> Pnext != VHead)
		    IP_SET_INTERNAL_VRTX(V3);
		else
		    V3 -> Tags = VRest -> Tags;

                Collinear = 
                    GMCollinear3Pts(V1 -> Coord, V2 -> Coord, V3 -> Coord);

                if (Collinear) {
                    if (IPGlblGenDegenPolys) {
                        KeepDegenerated = TRUE;
                    }
                    else {
		        IPFreeVertex(V1);
		        IPFreeVertex(V2);
                    }
		}
                if (!Collinear || KeepDegenerated) {
		    if (Pl -> PVertex == NULL) {
		        Pl -> PVertex = V1;
		    }
		    else {
		        PlNew = IPAllocPolygon(0, V1, PlNew);
			IRIT_PLANE_COPY(PlNew -> Plane, Pl -> Plane);
			IP_SET_PLANE_POLY(PlNew);
			PlNew -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr);
                        if (KeepDegenerated)
                            AttrSetIntAttrib(&PlNew -> Attr, 
					     IP_ATTRIB_DEGEN_POLY, TRUE);
		    }
		}

		VLast = V3;
		VLastTags = VRest -> Tags;

		IPFreeVertex(VRest);
		VRest = VRestNext;
	    }

	    if (PlNew != NULL) {
		Pl -> Pnext = PlNew;
		IPGetLastPoly(PlNew) -> Pnext = PlNext;
	    }

	    if (FreeVHead)
		IPFreeVertexList(VHead);

            if (FirstIsDegenerated)
                AttrSetIntAttrib(&Pl -> Attr, IP_ATTRIB_DEGEN_POLY, TRUE);
	}

	Pl = PlNext;
    }

    /* Purge empty polygons. */
    for (Pl = PolyObj -> U.Pl; Pl != NULL && Pl -> PVertex == NULL; ) {
        PolyObj -> U.Pl = Pl -> Pnext;
	IPFreePolygon(Pl);
	Pl = PolyObj -> U.Pl;
    }
    if ((Pl = PolyObj -> U.Pl) != NULL) {
        while (Pl -> Pnext != NULL) {
	    if (Pl -> Pnext -> PVertex == NULL) {
		IPPolygonStruct
		    *PlTmp = Pl -> Pnext;

		Pl -> Pnext = PlTmp -> Pnext;
		IPFreePolygon(PlTmp);
	    }
	    else
		Pl = Pl -> Pnext;
	}
    }

    /* Update the plane equation. */
    for (Pl = PolyObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        IrtVecType OldPlNrml;

	IRIT_VEC_COPY(OldPlNrml, Pl -> Plane);

	if (IPUpdatePolyPlane(Pl)) {
	    if (IRIT_DOT_PROD(OldPlNrml, Pl -> Plane) < 0.0)
	        IRIT_PLANE_SCALE(Pl -> Plane, -1.0);
	}
    }

    return PolyObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new polygonal objects out of given one, that contains only     M
* triangles.  Non convex polygons are split to convex one which, in turn,    M
* converted to triangles.  Collinear points are used and split at.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   PolyObj:   Polygonal object to split into triangles.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A polygonal object containing only triangles         M
*		representing the same model as PolyObj.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   ConvexPolyObjectN, GMConvertPolysToNGons, GMLimitTrianglesEdgeLen        M
*   GMConvertPolysToTriangles, GMConvertPolysToRectangles		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMConvertPolysToTriangles2                                               M
*****************************************************************************/
IPObjectStruct *GMConvertPolysToTriangles2(IPObjectStruct *PolyObj)
{
    IPPolygonStruct *Pl;
    int IsCirc = IPSetPolyListCirc(FALSE);

    IPSetPolyListCirc(IsCirc);	      /* Restore state, now that we know it. */

    PolyObj = GMConvexPolyObjectN(PolyObj);

    for (Pl = PolyObj -> U.Pl; Pl != NULL; ) {
        IPVertexStruct
	    *VHead = Pl -> PVertex;
	IPPolygonStruct
	    *PlNext = Pl -> Pnext;

	IPGetLastVrtx(VHead) -> Pnext = VHead;             /* Make circular. */

	while (IPVrtxListLen(VHead) > 3) {
	    IrtRType A,
	        ExtremeAngle = -2.0;
	    IrtVecType Vec1, Vec2;
	    IPVertexStruct *V, *VNext, *VNextNext, *V1, *V2, *V3,
	        *VBest = NULL;

	    /* Find a vertex with smallest turning angle to cut triangle at. */
	    V = VHead;
	    VNext = V -> Pnext;
	    VNextNext = VNext -> Pnext;
	    IRIT_VEC_SUB(Vec2, VNext -> Coord, V -> Coord);
	    IRIT_VEC_NORMALIZE(Vec2);

	    do {
		IRIT_VEC_COPY(Vec1, Vec2);
		IRIT_VEC_SUB(Vec2, VNextNext -> Coord, VNext -> Coord);
		IRIT_VEC_NORMALIZE(Vec2);

		A = IRIT_DOT_PROD(Vec1, Vec2);
		if (ExtremeAngle < A) {
		    ExtremeAngle = A;
		    VBest = V;
		}

		V = VNext;
		VNext = VNextNext;
		VNextNext = VNextNext -> Pnext;
	    }
	    while (V != VHead);

	    /* VBest and its two next vertices are the most collinear.      */
	    /* Cut starting from VBest -> Pnext.			    */
	    V = VHead = VBest -> Pnext;
	    assert(V != NULL);
	    VNext = V -> Pnext;
	    VNextNext = VNext -> Pnext;
	    V -> Pnext = VNext -> Pnext;
	    IP_SET_INTERNAL_VRTX(V);

	    /* Cut out a triangle here. */
	    V3 = IPAllocVertex2(NULL);
	    V2 = VNext;
	    V1 = IPAllocVertex2(V2);

	    /* Build the new triangle. */
	    V2 -> Pnext = V3;
	    if (IsCirc)
	        V3 -> Pnext = V1;
	    VERTEX_COPY(V1, V);
	    VERTEX_COPY(V3, VNextNext);
	    IP_SET_INTERNAL_VRTX(V3);
	    Pl -> Pnext = IPAllocPolygon(0, V1, Pl -> Pnext);
	    IRIT_PLANE_COPY(Pl -> Pnext -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(Pl -> Pnext);
	    Pl -> Pnext -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr);
	}

	/* Make original poly circular if so desired. */
	assert(VHead -> Pnext -> Pnext -> Pnext == VHead);
	if (!IsCirc)
	    VHead -> Pnext -> Pnext -> Pnext = NULL;
	Pl -> PVertex = VHead;

	Pl = PlNext;
    }

    /* Purge empty polygons. */
    for (Pl = PolyObj -> U.Pl; Pl != NULL && Pl -> PVertex == NULL; ) {
        PolyObj -> U.Pl = Pl -> Pnext;
	IPFreePolygon(Pl);
	Pl = PolyObj -> U.Pl;
    }
    if ((Pl = PolyObj -> U.Pl) != NULL) {
        while (Pl -> Pnext != NULL) {
	    if (Pl -> Pnext -> PVertex == NULL) {
		IPPolygonStruct
		    *PlTmp = Pl -> Pnext;

		Pl -> Pnext = PlTmp -> Pnext;
		IPFreePolygon(PlTmp);
	    }
	    else
		Pl = Pl -> Pnext;
	}
    }

    return PolyObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new polygonal objects out of given one, that contains only     M
* triangles.  Non convex polygons are split to convex one which, in turn,    M
* converted to triangles.  This version triangulates convex polygons by	     M
* connecting each vertex to the geometric centroid, hence co-linear pts are  M
* not purged, and do not create degenerate triangles. All other fields are   M
* not handled (normals, etc').						     M
*                                                                            *
* PARAMETERS:                                                                M
*   PolyObj:   Polygonal object to split into triangles.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   A polygonal object containing only triangles         M
*		representing the same model as PolyObj.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   ConvexPolyObjectN, GMConvertPolysToNGons, GMLimitTrianglesEdgeLen        M
*   GMConvertPolysToTriangles, GMConvertPolysToRectangles		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMConvertPolysToTrianglesIntrrPt                                         M
*****************************************************************************/
IPObjectStruct *GMConvertPolysToTrianglesIntrrPt(IPObjectStruct *PolyObj)
{
    IrtPtType TmpCentroid;
    IPPolygonStruct *PolygonIter, *NewTriangle,
	*TrianglePolygonList = NULL;
    int i, 
	IsCirc = IPSetPolyListCirc(FALSE);

    IPSetPolyListCirc(IsCirc);	      /* Restore state, now that we know it. */

    PolyObj = GMConvexPolyObjectN(PolyObj);
    PolygonIter = PolyObj -> U.Pl;
    while (PolygonIter != NULL) {
	int VNum = IPVrtxListLen(PolygonIter -> PVertex);
	/* The convention for the input polygons is that the start and end   */
	/* are the same (triangle has four points), and the output triangles */
	/* will have only the three distinct vertices.                       */
	if (VNum < 4) {
	    PolygonIter = PolygonIter -> Pnext;
	    continue;
	}
	if (VNum > 4) {
	    if (GMPolyCentroid(PolygonIter, TmpCentroid)) {
		IPVertexStruct *NewTrVertex, *CurPlVertex, *TrVertexList,
		    *VertexIter = PolygonIter -> PVertex;
		while (VertexIter -> Pnext != NULL) {
		    /* A valid next edge - defining a triangle together with */
		    /* the centroid, as a third vertex.                      */
		    TrVertexList = NULL;
		    for (i = 0; i < 2; i++) {
			switch (i) {
			    case 0:
			        CurPlVertex = VertexIter -> Pnext;
			        break;
			    case 1:
			        CurPlVertex = VertexIter;;
			        break;
			}
			NewTrVertex = IPAllocVertex2(NULL);
			IRIT_GEN_COPY(NewTrVertex -> Coord,
				      CurPlVertex -> Coord, 
				      3 * sizeof(CagdRType));
			IRIT_LIST_PUSH(NewTrVertex, TrVertexList);
		    }
		    NewTrVertex = IPAllocVertex2(NULL);
		    IRIT_GEN_COPY(NewTrVertex -> Coord, TmpCentroid, 
				  3 * sizeof(CagdRType));
		    IRIT_LIST_PUSH(NewTrVertex, TrVertexList);
		    NewTriangle = IPAllocPolygon(0, TrVertexList, NULL);
		    IRIT_LIST_PUSH(NewTriangle, TrianglePolygonList);
		    VertexIter = VertexIter -> Pnext;
		}
	    }
	    else /* Failed to compute the centroid for some reason. */
		assert(0);
	}
	else { 
	    /* This is already a triangle. Change to three vertices form and */
	    /* simply copy.                                                  */
	    IPVertexStruct *TmpV;

	    IRIT_LIST_POP(TmpV, PolygonIter -> PVertex);
	    IPFreeVertex(TmpV);
	    NewTriangle = IPCopyPolygon(PolygonIter);
	    IRIT_LIST_PUSH(NewTriangle, TrianglePolygonList);
	}
	PolygonIter = PolygonIter -> Pnext;
    }
    IPFreePolygonList(PolyObj -> U.Pl);
    PolyObj -> U.Pl = TrianglePolygonList;
    return PolyObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new polygonal objects out of given one, that contains only     M
* rectangles.  Non convex polygons with an empty kernel will generate self-  M
* intersecting results.							     M
*   By selecting a centroid location in the kernel of each polygon and       M
* connecting that centroid location to all the middle two adjacent edges,    M
* for all edges, n rectangular regions are defined for each n-gon, about     M
* half the (edge) size.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   PolyObj:    Polygonal object to split into rectangles, in place.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *: Return list of rectangulr polygons.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   ConvexPolyObjectN, GMConvertPolysToNGons, GMLimitTrianglesEdgeLen        M
*   GMConvertPolysToTriangles, GMConvertPolysToTriangles2		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMConvertPolysToRectangles                                               M
*****************************************************************************/
IPObjectStruct *GMConvertPolysToRectangles(IPObjectStruct *PolyObj)
{
    IPPolygonStruct *Pl, 
        *Pls = PolyObj -> U.Pl,
	*ResPl = NULL; 
    int IsCirc = IPSetPolyListCirc(FALSE); 

    IPSetPolyListCirc(IsCirc); 	     /* Restore state, now that we know it. */
    
    do {
      IPVertexStruct *V, *VN, *VNN, *VNNN, *V1, *VTmp; 
	int Len, 
	    Start = TRUE; 
	IrtPtType Mid; 

	IRIT_LIST_POP(Pl, Pls);
	V1 = Pl -> PVertex; 
	Len = IPVrtxListLen(V1);

	if (Len > 2) {
	    /* Compute a centroid point. */
	    if (!GMFindPtInsidePolyKernel(V1, Mid))
		GMComputeAverageVertex(V1, Mid, 1.0);

	    /*  Insert two vertices in the middle of every edge of Pl.       */
	    for (V = V1, VN = V1 -> Pnext;  
		 Start || (V != NULL && V != V1); 
		 V = VN, VN = VN -> Pnext) {
		Start = FALSE; 

		/* Insert two new vertices between V and V -> Pnext. */
		V -> Pnext = VTmp = IPAllocVertex(V -> Tags, NULL,
						  V -> Pnext); 
		IRIT_VEC_COPY(VTmp -> Normal, V -> Normal);
		IRIT_PT_BLEND(VTmp -> Coord, V -> Coord, 
			      VTmp -> Pnext -> Coord, 0.5); 

		V -> Pnext = VTmp = IPAllocVertex(V -> Tags, NULL,
						  V -> Pnext); 
		IRIT_VEC_COPY(VTmp -> Normal, VTmp -> Pnext -> Normal);
		IRIT_PT_COPY(VTmp -> Coord, VTmp -> Pnext -> Coord); 
	    }

	    /*   At this point every original edge has two identical         */
	    /* interior vertex.  We march on every 2nd such interior         */
	    /* vertex using V and extract rectangular domains.               */
	    /*   Add the centroid of the polygon with three boundary         */
	    /* vertices, forming a rectangle.				     */
	    for (V = V1 -> Pnext -> Pnext; 
		 V -> Pnext != V1; 
		 V = VNNN) {
		IPPolygonStruct *NewPl;

		VNN = V -> Pnext -> Pnext; 
		VNNN = VNN -> Pnext; 
		VNN -> Pnext = IPAllocVertex2(V);/* Make rectangle circular. */

		IRIT_PT_COPY(VNN -> Pnext -> Coord, Mid); 

		NewPl = IPAllocPolygon(Pl -> Tags, V, NULL); 
		IRIT_PLANE_COPY(NewPl -> Plane, Pl -> Plane); 
		IP_SET_PLANE_POLY(NewPl); 

		NewPl -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr); 
		IP_RST_BBOX_POLY(NewPl); 
		IRIT_LIST_PUSH(NewPl, ResPl);
	    }
	    V -> Pnext -> Pnext -> Pnext = IPAllocVertex2(V); 
	    IRIT_PT_COPY(V -> Pnext -> Pnext -> Pnext -> Coord, Mid); 
	    IRIT_LIST_PUSH(Pl, ResPl);
	}
	else
	    IPFreePolygon(Pl);
    }
    while (Pls != NULL);

    PolyObj -> U.Pl = ResPl;
    return PolyObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the distance between two polylines.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Pl1, Pl2:  To compute the minimal distance between the end points of.    *
*   Start1, Start2:  TRUE if start of polyline, FALSE if end of polyline.    *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtRType:   Minimal distance squared computed.			     *
*****************************************************************************/
static IrtRType GetPlPlDist(IPPolygonStruct *Pl1,
			    IPPolygonStruct *Pl2,
			    int *Start1,
			    int *Start2)
{
    IrtRType
	MinSqr = IRIT_INFNTY;
    IPVertexStruct
	*V1Start = Pl1 -> PVertex,
	*V1End = IPGetLastVrtx(V1Start),
        *V2Start = Pl2 -> PVertex,
        *V2End = IPGetLastVrtx(V2Start);

    if (V1Start == V1End) {
        if (V2Start == V2End) {
	    TEST_CLOSEST_DIST(V1Start, V2Start, 1, 1);
	}
	else {
	    TEST_CLOSEST_DIST(V1Start, V2Start, 1, 1);
	    TEST_CLOSEST_DIST(V1Start, V2End,   1, 0);
	}
    }
    else {
        if (V2Start == V2End) {
	    TEST_CLOSEST_DIST(V1Start, V2Start, 1, 1);
	    TEST_CLOSEST_DIST(V1End,   V2Start, 0, 1);
	}
	else {
	    TEST_CLOSEST_DIST(V1Start, V2Start, 1, 1);
	    TEST_CLOSEST_DIST(V1End,   V2Start, 0, 1);
	    TEST_CLOSEST_DIST(V1Start, V2End,   1, 0);
	    TEST_CLOSEST_DIST(V1End,   V2End,   0, 0);
	}
    }

    return MinSqr;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Merges separated polylines into longer ones, in place, as possible.      M
* Given a list of polylines, matches end points and merged as possible       M
* polylines with common end points, in place.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Polys:       Polylines to merge, in place.                               M
*   Eps:	 Epsilon of similarity to merge points at.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:  Merged as possible polylines.                        M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMMergeGeometry, MvarPolyMergePolylines, GMMergePolylines2               M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMMergePolylines, merge, polyline                                        M
*****************************************************************************/
IPPolygonStruct *GMMergePolylines(IPPolygonStruct *Polys, IrtRType Eps)
{
    int i,
	NumOfPolys = IPPolyListLen(Polys);
    IPPolygonStruct **PlsVec, *Pl;

    if (NumOfPolys < 2)
        return Polys;

    PlsVec = (IPPolygonStruct **) IritMalloc(sizeof(IPPolygonStruct *) *
					     NumOfPolys);
    for (i = 0, Pl = Polys; i < NumOfPolys; i++, Pl = Pl -> Pnext)
        PlsVec[i] = Pl;

    NumOfPolys = GMMergeGeometry((void **) PlsVec, NumOfPolys, Eps, IRIT_UEPS,
				 NULL, NULL, NULL, NULL);

    for (i = 1, Pl = Polys = PlsVec[0]; i < NumOfPolys; i++) {
        Pl -> Pnext = PlsVec[i];
	Pl = PlsVec[i];
    }
    Pl -> Pnext = NULL;

    IritFree(PlsVec);

    return Polys;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Routine to compare two GMPlPlDistStruct for sorting purposes.            *
*                                                                            *
* PARAMETERS:                                                                *
*   PlPlDist1, PlPlDist2:  Two pointers to PlPlDist structs.                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:   >0, 0, or <0 as the relation between the two distances (squared). *
*****************************************************************************/
#if defined(ultrix) && defined(mips)
static int ComparePlPlDist(VoidPtr PlPlDist1, VoidPtr PlPlDist2)
#else
static int ComparePlPlDist(const VoidPtr PlPlDist1, const VoidPtr PlPlDist2)
#endif /* ultrix && mips (no const support) */
{
    IrtRType
	Diff = ((GMPlPlDistStruct *) PlPlDist1) -> MinDistSqr -
               ((GMPlPlDistStruct *) PlPlDist2) -> MinDistSqr;

    return IRIT_SIGN(Diff);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Connect the list of points into polylines by connecting the closest      M
* point pairs, until the distances between adjacent points/polylines is more M
* than MaxTol.  Points are assumed to be in E3.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   PtsList:      Point list to connect into polylines.                      M
*   MaxTol:       Maximum distance allowed to connect to points.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   Connected polylines, upto MaxTol tolerance.         M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMMatchPointListIntoPolylines                                            M
*****************************************************************************/
IPPolygonStruct *GMMatchPointListIntoPolylines(IPObjectStruct *PtsList,
					       IrtRType MaxTol)
{
    int i, LastN, n;
    IrtRType
	MaxTolSqr = IRIT_SQR(MaxTol);
    IPPolygonStruct *PllList;
    IPObjectStruct *PObj;

    PtsList = IPCopyObject(NULL, PtsList, FALSE);
    IPCoercePtsListTo(PtsList, CAGD_PT_E3_TYPE);

    /* Convert the list object to polyline linked list with one vertex in   */
    /* each polyline so we can start and match-connect them.		    */
    for (PllList = NULL, i = 0;
	 (PObj = IPListObjectGet(PtsList, i++)) != NULL;
	 ) {
        IPVertexStruct *V;
	CagdRType
	    *Coords = PObj -> U.CtlPt.Coords;

	PllList = IPAllocPolygon(0, V = IPAllocVertex2(NULL), PllList);
	CagdCoerceToE3(V -> Coord, &Coords, -1, PObj -> U.CtlPt.PtType);
    }
    IPFreeObject(PtsList);

    n = IPPolyListLen(PllList);
    do {
	int j,
	    *InvIdxMap = (int *) IritMalloc(sizeof(int) * n);
	GMPlPlDistStruct
	    *Plls = (GMPlPlDistStruct *)
	        IritMalloc(sizeof(GMPlPlDistStruct) * n);

	LastN = n;

#       ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPtsToPlsSteps, FALSE)
	        IRIT_INFO_MSG_PRINTF("Num of Polylines = %5d\n", n);
	}
#       endif /* DEBUG */

	/* Make the list into an array. */
	for (i = 0; i < n; i++)
	    IRIT_LIST_POP(Plls[i].Pl, PllList);

	/* Compute minimal distance between all points/polylines. */
	for (i = 0; i < n; i++) {
	    IrtRType
		MinDistSqr = IRIT_INFNTY;
	    int MinStart1 = -1,
		MinStart2 = -1,
		MinIdx = -1;

	    for (j = i + 1; j < n; j++) {
	        int Start1 = 0,
		    Start2 = 0;
	        IrtRType
		    DistSqr = GetPlPlDist(Plls[i].Pl, Plls[j].Pl,
					  &Start1, &Start2);

		if (DistSqr < MinDistSqr) {
		    MinIdx = j;
		    MinStart1 = Start1;
		    MinStart2 = Start2;
		    MinDistSqr = DistSqr;  
		}
	    }

	    Plls[i].Idx = i;
	    Plls[i].MinIdx = MinIdx;
	    Plls[i].MinStart1 = MinStart1;
	    Plls[i].MinStart2 = MinStart2;
	    Plls[i].MinDistSqr = MinDistSqr;
	}

	/* Sort the array based on MinDistSqr slot and build an inverse map. */
	qsort(Plls, n - 1, sizeof(GMPlPlDistStruct), ComparePlPlDist);
	for (i = 0; i < n; i++)
	    InvIdxMap[Plls[i].Idx] = i;

	/* Merge all polyline we can (no merged polyline will be remerged). */
	for (i = 0, PllList = NULL;
	     i < n - 1 && Plls[i].MinDistSqr < MaxTolSqr;
	     i++) {
	    j = InvIdxMap[Plls[i].MinIdx];
	    if (Plls[i].Pl != NULL && Plls[j].Pl != NULL) {
		/* Merge plln index i with plln index j, j > i. */
	        if (Plls[i].MinStart1)
		    Plls[i].Pl -> PVertex =
		        IPReverseVrtxList2(Plls[i].Pl -> PVertex);
		if (!Plls[i].MinStart2)
		    Plls[j].Pl -> PVertex =
		        IPReverseVrtxList2(Plls[j].Pl -> PVertex);
		IPGetLastVrtx(Plls[i].Pl -> PVertex) -> Pnext =
		                                        Plls[j].Pl -> PVertex;
		Plls[j].Pl -> PVertex = NULL;
		IPFreePolygon(Plls[j].Pl);
		IRIT_LIST_PUSH(Plls[i].Pl, PllList);
		Plls[j].Pl = Plls[i].Pl = NULL;
	    }
	}

	/* Regroup into a linked list the rest of the polylines. */
	for (i = 0; i < n; i++)
	    if (Plls[i].Pl != NULL)
	        IRIT_LIST_PUSH(Plls[i].Pl, PllList);

	IritFree(Plls);
	IritFree(InvIdxMap);

	n = IPPolyListLen(PllList);
    }
    while (n < LastN);

    return PllList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Regularize a polygonal model by eliminating all T junction in the        M
* polygonal mesh.                                                            M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:            A polygonal object to regularize.                       M
*   SplitCollinear:  TRUE to also split polygons at collienar edges.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   Regularized object.                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMRegularizePolyModel                                                    M
*****************************************************************************/
IPObjectStruct *GMRegularizePolyModel(IPObjectStruct *PObj,
				      int SplitCollinear)
{
    IRIT_STATIC_DATA AttribNumType /* Persistent store to unique Attrs IDs. */
        AttVSListID = ATTRIB_NAME_BAD_NUMBER,
        AttrRegularID = ATTRIB_NAME_BAD_NUMBER;
    int TrisOnly = TRUE;
    IPPolygonStruct *Pl;

    IP_ATTR_INIT_UNIQUE_ID_NUM(AttVSListID, "_vslist");
    IP_ATTR_INIT_UNIQUE_ID_NUM(AttrRegularID, "regular");

    if (!IP_IS_POLY_OBJ(PObj) || !IP_IS_POLYGON_OBJ(PObj))
	return NULL;

    PObj = IPCopyObject(NULL, PObj, FALSE);
    IP_SET_POLYGON_OBJ(PObj);
    GMCleanUpPolygonList(&PObj -> U.Pl, IRIT_EPS);

    BoolGenAdjacencies(PObj);

    /* Find all edges in polygons to split according to adj. edges' sizes. */
    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	IPVertexStruct *VAdj,
	    *V = Pl -> PVertex;

	if (IPVrtxListLen(V) != 3)   /* Do we have triangles only in input? */
	    TrisOnly = FALSE;

	do {
	    if (V -> PAdj != NULL &&
		(VAdj = GetAdjVertex(V, Pl, V -> PAdj)) != NULL) {
		if ((IRIT_PT_APX_EQ(V -> Coord, VAdj -> Coord) &&
		     IRIT_PT_APX_EQ(V -> Pnext -> Coord,
				    VAdj -> Pnext -> Coord)) ||
		    (IRIT_PT_APX_EQ(V -> Coord, VAdj -> Pnext -> Coord) &&
		     IRIT_PT_APX_EQ(V -> Pnext -> Coord, VAdj -> Coord))) {
		    /* This edge is identical on both sides - regular! */
		    V -> PAdj = VAdj -> PAdj = NULL;     /* Clear adjacency. */
		}
		else {
		    /* Need to split this edge - update this adj. info. */
		    UpdateAdjPolys(V, Pl, VAdj, V -> PAdj);
		    V -> PAdj = VAdj -> PAdj = NULL;     /* Clear adjacency. */
		}
	    }
	    V = V -> Pnext;
	}
	while (V != NULL && V != Pl -> PVertex);
    }

    /* Insert the new vertices into the edges that require splitting. */
    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	IPVertexStruct *VSList, *VNext, *VMid,
	    *V = Pl -> PVertex;

	do {
	    VNext = V -> Pnext;

	    if ((VSList = AttrGetRefPtrAttrib2(V -> Attr,
					       AttVSListID)) != NULL) {
	        IPGetLastVrtx(VSList) -> Pnext = VNext;
		for (VMid = VSList; VMid != VNext; VMid = VMid -> Pnext) {
		    GMInterpVrtxNrmlBetweenTwo(VMid, V, V -> Pnext);
		    GMInterpVrtxRGBBetweenTwo(VMid, V, V -> Pnext);
		    GMInterpVrtxUVBetweenTwo(VMid, V, V -> Pnext);
		}
		V -> Pnext = VSList;
		AttrFreeOneAttribute(&V -> Attr, "_vslist");
	    }

	    V = VNext;
	}
	while (V != NULL && V != Pl -> PVertex);
    }

    if (SplitCollinear) {
        /* Split polygons that has collinear adj. edges due to vertex       */
        /* insertions.							    */
        Pl = GMSplitPolysAtCollinearVertices(PObj -> U.Pl);
	IPFreePolygonList(PObj -> U.Pl);
	PObj -> U.Pl = Pl;
    }

    if (TrisOnly) {
	IPObjectStruct
            *PRet = GMConvertPolysToTriangles(PObj);

	IPFreeObject(PObj);

	AttrSetObjectIntAttrib2(PRet, AttrRegularID, TRUE);

	return PRet;
    }
    else {
	AttrSetObjectIntAttrib2(PObj, AttrRegularID, TRUE);

	return PObj;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Update the two edges that are collinear with vertices to be split at     *
* as/if necessary.  Compute length along the edges' direction and insert     *
* into vertex split list if inside the edge.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   V1:       First vertex in edge (V1, V1 -> Pnext) of Pl1.                 *
*   Pl1:      First polygon holding V1.                                      *
*   V2:       First vertex in edge (V2, V2 -> Pnext) of Pl2.                 *
*   Pl2:      First polygon holding V2.                                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UpdateAdjPolys(IPVertexStruct *V1,
			   IPPolygonStruct *Pl1,
			   IPVertexStruct *V2,
			   IPPolygonStruct *Pl2)
{
    IrtRType Edge1Len, Edge2Len, t,
	*Pt1 = V1 -> Coord,
	*Pt1Next = V1 -> Pnext -> Coord,
	*Pt2 = V2 -> Coord,
	*Pt2Next = V2 -> Pnext -> Coord;
    IrtVecType Edge1Dir, Edge2Dir, Vec;

    /* Do the first, V2 inside V1, direction. */
    IRIT_VEC_SUB(Edge1Dir, Pt1Next, Pt1);
    Edge1Len = IRIT_VEC_LENGTH(Edge1Dir);

    IRIT_VEC_SUB(Vec, Pt2, Pt1);
    t = IRIT_VEC_LENGTH(Vec) / Edge1Len;
    if (IRIT_DOT_PROD(Vec, Edge1Dir) < 0.0)
        t = -t;
    if (0.0 < t && !IRIT_APX_EQ(0.0, t) && 1.0 > t && !IRIT_APX_EQ(1.0, t))
	InsertVertexToSplitList(V1, Pt2, t);

    IRIT_VEC_SUB(Vec, Pt2Next, Pt1);
    t = IRIT_VEC_LENGTH(Vec) / Edge1Len;
    if (IRIT_DOT_PROD(Vec, Edge1Dir) < 0.0)
        t = -t;
    if (0.0 < t && !IRIT_APX_EQ(0.0, t) && 1.0 > t && !IRIT_APX_EQ(1.0, t))
	InsertVertexToSplitList(V1, Pt2Next, t);

    /* Do the second, V1 inside V2, direction. */
    IRIT_VEC_SUB(Edge2Dir, Pt2Next, Pt2);
    Edge2Len = IRIT_VEC_LENGTH(Edge2Dir);

    IRIT_VEC_SUB(Vec, Pt1, Pt2);
    t = IRIT_VEC_LENGTH(Vec) / Edge2Len;
    if (IRIT_DOT_PROD(Vec, Edge2Dir) < 0.0)
        t = -t;
    if (0.0 < t && !IRIT_APX_EQ(0.0, t) && 1.0 > t && !IRIT_APX_EQ(1.0, t))
	InsertVertexToSplitList(V2, Pt1, t);

    IRIT_VEC_SUB(Vec, Pt1Next, Pt2);
    t = IRIT_VEC_LENGTH(Vec) / Edge2Len;
    if (IRIT_DOT_PROD(Vec, Edge2Dir) < 0.0)
        t = -t;
    if (0.0 < t && !IRIT_APX_EQ(0.0, t) && 1.0 > t && !IRIT_APX_EQ(1.0, t))
	InsertVertexToSplitList(V2, Pt1Next, t);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Insert a new vertex to split the edge (V, V -> Pnext) into two at, to    *
* the split vertex list.                                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   V:     Pointer on the edge (V, V -> Pnext) that Pos is a middle position *
*          to split at.							     *
*   Pos:   New location inside (V, V -> Pnext) to split at.                  *
*   t:     location along the edge where Pos is.			     *
*          (t = 0 is V, t = 1 is V -> Pnext).				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void InsertVertexToSplitList(IPVertexStruct *V,
				    IrtRType *Pos,
				    IrtRType t)
{
    IRIT_STATIC_DATA AttribNumType /* Persistent store to unique Attrs IDs. */
        AttVSListID = ATTRIB_NAME_BAD_NUMBER,
        AttrVSTID = ATTRIB_NAME_BAD_NUMBER;
    IPVertexStruct *VSList,
        *NewV = IPAllocVertex(V -> Tags, NULL, NULL);

    IP_ATTR_INIT_UNIQUE_ID_NUM(AttVSListID, "_vslist");
    IP_ATTR_INIT_UNIQUE_ID_NUM(AttrVSTID, "_vst");

    VSList = AttrGetRefPtrAttrib2(V -> Attr, AttVSListID);

    IRIT_PT_COPY(NewV -> Coord, Pos);
    IRIT_VEC_BLEND(NewV -> Normal, V -> Pnext -> Normal, V -> Normal, t);
    AttrSetRealAttrib2(&NewV -> Attr, AttrVSTID, t);

    /* Insert the new vertex into VSList at the proper order. */
    if (VSList == NULL) {
        VSList = NewV;
    }
    else {
        IrtRType Tl;
	IPVertexStruct *VSHead,
	    *VSPrev = NULL;

	/* Find the proper location to insert NewV into VSList. */
	VSHead = VSList;
        do {
	    Tl = AttrGetRealAttrib2(VSHead -> Attr, AttrVSTID);
	    if (IRIT_APX_EQ(Tl, t)) {
	        /* Drop this vertex - it is already in. */
	        IPFreeVertex(NewV);
	        return;
	    }

	    if (Tl > t)
	        break;

	    VSPrev = VSHead;
	    VSHead = VSHead -> Pnext;
	}
	while (VSHead != NULL);

	if (VSPrev == NULL) {
	    NewV -> Pnext = VSList;
	    VSList = NewV;
	}
	else {
	    NewV -> Pnext = VSPrev -> Pnext;
	    VSPrev -> Pnext = NewV;
	}
    }

    AttrSetRefPtrAttrib2(&V -> Attr, AttVSListID, VSList);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Splits the given polygons in vertices that connect two adjacent          M
* collinear edges. 							     M
*   Polygons are assumed convex other than this collinearity conditions.     M
*                                                                            *
* PARAMETERS:                                                                M
*   Pls:   List of polygons to split at collinear edges.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *: New list of polygons with no colinear adjacent edges. M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMSplitPolysAtCollinearVertices                                          M
*****************************************************************************/
IPPolygonStruct *GMSplitPolysAtCollinearVertices(IPPolygonStruct *Pls)
{
    IPPolygonStruct *Pl, *PlNext,
	*PlPrev = NULL;

    Pls = IPCopyPolygonList(Pls);

    for (Pl = Pls; Pl != NULL; ) {
        int WasSplit = FALSE;
	IPVertexStruct *VNext,
	    *V = Pl -> PVertex;

	do {
	    VNext = V -> Pnext;
	    if (GMCollinear3Pts(V -> Coord, VNext -> Coord,
				VNext -> Pnext -> Coord)) {
		/* Split polygon into two at VNext and chain second poly in. */
	        if ((PlNext = GMSplitPolyInPlaceAtVertex(Pl, VNext)) == NULL) {
		    /* Remove Pl from polygon list - a degenerated polygon. */
		    if (PlPrev == NULL) {
			Pls = Pls -> Pnext;
			IPFreePolygon(Pl);
			Pl = Pls;
		    }
		    else {
			PlPrev -> Pnext = Pl -> Pnext;
			IPFreePolygon(Pl);
			Pl = PlPrev -> Pnext;
		    }
		}
		else {
		    PlNext -> Pnext = Pl -> Pnext;
		    Pl -> Pnext = PlNext;
		}
		WasSplit = TRUE;
		break;
	    }
	    V = VNext;
	}
	while (V != NULL && V != Pl -> PVertex);

	if (!WasSplit) {
	    PlPrev = Pl;
	    Pl = Pl -> Pnext;
	}
    }

    return Pls;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Splits the given convex polygon, in place, into two, returning second    M
* half of the polygon while updating Pl to hold the first half.              M
*   Polygon is split so that VHead is on border between the two polygons.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Pl:      Convex polygon to split into two.                               M
*   VHead:   Vertex to split Pl at.	                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   The second half of the splitted polygon (first half M
*			 is returned, in place, in Pl).			     M
*			 This function returns a NULL if split failed due to M
*			 the fact that the polygon degenerated into a line.  M
*			 Pl is not affected if NULL is returned.             M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMSplitPolyInPlaceAt2Vertices				             M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMSplitPolyInPlaceAtVertex                                               M
*****************************************************************************/
IPPolygonStruct *GMSplitPolyInPlaceAtVertex(IPPolygonStruct *Pl,
					    IPVertexStruct *VHead)
{
    IPVertexStruct *VHead2, *V2, *VTmp,
	*VNext = VHead -> Pnext,
	*V = VNext -> Pnext;
    IPPolygonStruct *PlRet;

    do {
        if (V -> Pnext != VHead &&
	    !GMCollinear3Pts(VHead -> Coord, VNext -> Coord, V -> Coord)) {
	    /* The edge (VHead, V) is our splitting edge. */
	    VHead2 = IPAllocVertex(VHead -> Tags, NULL, VHead -> Pnext);
	    IRIT_PT_COPY(VHead2 -> Coord, VHead -> Coord);
	    IRIT_VEC_COPY(VHead2 -> Normal, VHead -> Normal);
	    VHead2 -> Attr = IP_ATTR_COPY_ATTRS(VHead -> Attr);

	    V2 = IPAllocVertex(V -> Tags, NULL, V -> Pnext);
	    IRIT_PT_COPY(V2 -> Coord, V -> Coord);
	    IRIT_VEC_COPY(V2 -> Normal, V -> Normal);
	    V2 -> Attr = IP_ATTR_COPY_ATTRS(V -> Attr);

	    for (VTmp = V2; VTmp -> Pnext != VHead; VTmp = VTmp -> Pnext);
	    VTmp -> Pnext = VHead2;

	    V -> Pnext = VHead;
	    IP_SET_INTERNAL_VRTX(V);
	    Pl -> PVertex = V;

	    VHead2 -> Pnext = V2;
	    IP_SET_INTERNAL_VRTX(VHead2);
	    PlRet = IPAllocPolygon(Pl -> Tags, V2, NULL);
	    IRIT_PLANE_COPY(PlRet -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlRet);
	    PlRet -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr);

	    IP_RST_BBOX_POLY(Pl);
	    IP_RST_BBOX_POLY(PlRet);

	    return PlRet;
	}

	V = V -> Pnext;
    }
    while (V != NULL && V -> Pnext != VHead);

    return NULL;	       /* This polygon has degenerated into a line. */
}


/*****************************************************************************
* DESCRIPTION:                                                               M
*   Splits the given convex polygon, in place, into two, returning second    M
* half of the polygon while updating Pl to hold the first half.              M
*                                                                            *
* PARAMETERS:                                                                M
*   Pl:      Convex polygon to split into two.                               M
*   V1:      First Vertex to split Pl at.	                             M
*   V2:	     Second Vertex to split Pl at.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   The second half of the splitted polygon (first half M
*			 is returned, in place, in Pl).			     M
*			 This function returns a NULL if split failed due to M
*			 the fact that the polygon degenerated into a line.  M
*			 Pl is not affected if NULL is returned.             M
*			 The second polygon is added as next to the first    M
*			 polygon.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMSplitPolyInPlaceAtVertex					             M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMSplitPolyInPlaceAt2Vertices                                            M
*****************************************************************************/
IPPolygonStruct *GMSplitPolyInPlaceAt2Vertices(IPPolygonStruct *Pl,
					       IPVertexStruct *V1,
					       IPVertexStruct *V2)
{
    int IPIsVertexListCirc = IPGetLastVrtx(Pl -> PVertex) -> Pnext != NULL;
    IPPolygonStruct *PlRet;
    IPVertexStruct *V11, *V21;

    /* Make list temporary circular, even if not. */
    if (!IPIsVertexListCirc)
        IPGetLastVrtx(Pl -> PVertex) -> Pnext = Pl -> PVertex;

    if (IRIT_PT_APX_EQ(V1 -> Coord, V2 -> Coord) || 
	IRIT_PT_APX_EQ(V1 -> Coord, V2 -> Pnext -> Coord) || 
	IRIT_PT_APX_EQ(V1 -> Pnext -> Coord, V2 -> Coord)) 
	return NULL;

    /* Duplicate the two vertices we are splitting at. */
    V11 = IPAllocVertex(V1 -> Tags, NULL, V1 -> Pnext);
    IRIT_PT_COPY(V11 -> Coord, V1 -> Coord);
    IRIT_VEC_COPY(V11 -> Normal, V1 -> Normal);
    V11 -> Attr = IP_ATTR_COPY_ATTRS(V1 -> Attr);

    V21 = IPAllocVertex(V2 -> Tags, NULL, V2 -> Pnext);
    IRIT_PT_COPY(V21 -> Coord, V2 -> Coord);
    IRIT_VEC_COPY(V21 -> Normal, V2 -> Normal);
    V21 -> Attr = IP_ATTR_COPY_ATTRS(V2 -> Attr);

    /* Make the two new links of the two pieces and create the second poly. */
    V1 -> Pnext = V21;
    V2 -> Pnext = V11;
	    
    IP_SET_INTERNAL_VRTX(V2);
    IP_SET_INTERNAL_VRTX(V1);
    PlRet = IPAllocPolygon(Pl -> Tags, V2, NULL);
    IRIT_PLANE_COPY(PlRet -> Plane, Pl -> Plane);
    IP_SET_PLANE_POLY(PlRet);
    PlRet -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr);

    IP_RST_BBOX_POLY(Pl);
    IP_RST_BBOX_POLY(PlRet);

    /* Make vertices' lists NULL terminated if input was NULL terminated. */
    if (!IPIsVertexListCirc) {
	IPGetLastVrtx(Pl -> PVertex) -> Pnext = NULL;
	IPGetLastVrtx(PlRet -> PVertex) -> Pnext = NULL;
    }

    PlRet -> Pnext = Pl -> Pnext;
    Pl -> Pnext = PlRet;

    return PlRet;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find out the adjacent edge to edge (V, V -> Pnext) in PAdj.              *
*                                                                            *
* PARAMETERS:                                                                *
*   V:         Of edge (V, V -> Pnext) to find adjacent edge on PAdj.        *
*   Pl:        Polygon containing vertex V and edge (V, V -> Pnext).	     *
*   PAdj:      The adjacent polygon to look in.                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPVertexStruct *:  The vertex the leads to the adjacent edge in PAdj.    *
*****************************************************************************/
static IPVertexStruct *GetAdjVertex(IPVertexStruct *V,
				    IPPolygonStruct *Pl,
				    IPPolygonStruct *PAdj)
{
    IPVertexStruct *VTmp;
    IrtVecType Dir, UDir;

    /* Lets see if we have the mutual pointer - a full adjacency. */
    VTmp = PAdj -> PVertex;
    do {
        if (VTmp -> PAdj == Pl)
	    return VTmp;

        VTmp = VTmp -> Pnext;
    }
    while (VTmp != NULL && VTmp != PAdj -> PVertex);

    /* Can still be a partial adjacency.  Check it now */
    IRIT_VEC_SUB(Dir, V -> Pnext -> Coord, V -> Coord);
    IRIT_VEC_COPY(UDir, Dir);
    IRIT_VEC_NORMALIZE(UDir);

    VTmp = PAdj -> PVertex;
    do {
	IrtPtType Pt1, Pt2;

        GMPointFromPointLine(VTmp -> Coord, V -> Coord, Dir, Pt1);
        GMPointFromPointLine(VTmp -> Pnext -> Coord, V -> Coord, Dir, Pt2);

	if (IRIT_PT_PT_DIST_SQR(VTmp -> Coord, Pt1) < IRIT_SQR(IRIT_EPS) &&
	    IRIT_PT_PT_DIST_SQR(VTmp -> Pnext -> Coord, Pt2) <
							 IRIT_SQR(IRIT_EPS)) {
	    IrtRType t1, t2;
	    IrtVecType Vec;

	    /* This edge is collinear with our edge - make sure these two   */
	    /* edge share a common ground.				    */
	    IRIT_VEC_SUB(Vec, Pt1, V -> Coord);
	    t1 = IRIT_DOT_PROD(Vec, UDir);
	    IRIT_VEC_SUB(Vec, Pt2, V -> Coord);
	    t2 = IRIT_DOT_PROD(Vec, UDir);
	    if (t1 > t2)
	        IRIT_SWAP(IrtRType, t1, t2);
	    if (t1 < 1.0 || t2 > 0.)
		return VTmp;
	}

        VTmp = VTmp -> Pnext;
    }
    while (VTmp != NULL && VTmp != PAdj -> PVertex);

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Splits all triangles that has edge length larger than MaxLen.  The       M
* output will have no edge in no triangle with length larger than Maxlen.    M
*                                                                            *
* PARAMETERS:                                                                M
*   OrigPls:   List of triangles.                                            M
*   MaxLen:    Maximum allowed length of an edge of a triangle.              M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *: Splitted polygons with edges smaller/equal to MaxLen. M
*                                                                            *
* SEE ALSO:                                                                  M
*   ConvexPolyObjectN, GMConvertPolysToNGons, GMConvertPolysToTriangles      M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMLimitTrianglesEdgeLen                                                  M
*****************************************************************************/
IPPolygonStruct *GMLimitTrianglesEdgeLen(const IPPolygonStruct *OrigPls,
					 IrtRType MaxLen)
{
    IrtRType
	MaxLenSqr = IRIT_SQR(MaxLen);
    IPPolygonStruct *Pl, *Pls;

    Pls = IPCopyPolygonList(OrigPls);

    for (Pl = Pls; Pl != NULL; ) {
	IrtRType D12, D23, D31;
	IPVertexStruct *V2, *V3, *Vl2,
	    *V1 = Pl -> PVertex;

	if (IPVrtxListLen(V1) != 3) {
	    GEOM_FATAL_ERROR(GEOM_ERR_TRIANGLES_ONLY);
	    return NULL;
	}

	V2 = V1 -> Pnext;
	V3 = V2 -> Pnext;

	D12 = IRIT_PT_PT_DIST_SQR(V1 -> Coord, V2 -> Coord);
	D23 = IRIT_PT_PT_DIST_SQR(V2 -> Coord, V3 -> Coord);
	D31 = IRIT_PT_PT_DIST_SQR(V3 -> Coord, V1 -> Coord);

	/* If we have an edge too long - duplicate the polygon as its next */
	/* and keep the two parts in Pl and pl -> Pnext.                   */
	if (D12 > MaxLenSqr || D23 > MaxLenSqr || D31 > MaxLenSqr) {
	    int DoNormals;

	    Pl -> Pnext = IPAllocPolygon(Pl -> Tags, IPCopyVertexList(V1),
					 Pl -> Pnext);
	    IRIT_PLANE_COPY(Pl -> Pnext -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(Pl -> Pnext);
	    Pl -> Pnext -> Attr = IP_ATTR_COPY_ATTRS(Pl -> Attr);
	    IP_RST_BBOX_POLY(Pl -> Pnext);

	    Vl2 = Pl -> Pnext -> PVertex;

	    if (IP_HAS_NORMAL_VRTX(V1) &&
		IP_HAS_NORMAL_VRTX(V2) &&
		IP_HAS_NORMAL_VRTX(V3)) {
		DoNormals = TRUE;
		IP_SET_NORMAL_VRTX(Vl2);
		IP_SET_NORMAL_VRTX(Vl2 -> Pnext);
		IP_SET_NORMAL_VRTX(Vl2 -> Pnext -> Pnext);
	    }
	    else {
		DoNormals = FALSE;
		IP_RST_NORMAL_VRTX(Vl2);
		IP_RST_NORMAL_VRTX(Vl2 -> Pnext);
		IP_RST_NORMAL_VRTX(Vl2 -> Pnext -> Pnext);
	    }

	    if (D12 >= D23 && D12 >= D31) {
		IRIT_PT_BLEND(Vl2 -> Coord, V1 -> Coord, V2 -> Coord, 0.5);

		GMInterpVrtxRGBBetweenTwo(Vl2, V1, V2);
		COPY_RGB_ATTR(V2, Vl2)
		GMInterpVrtxUVBetweenTwo(Vl2, V1, V2);
		COPY_UV_ATTR(V2, Vl2);

		if (DoNormals) {
		    GMInterpVrtxNrmlBetweenTwo(Vl2, V1, V2);
		    IRIT_VEC_COPY(V2 -> Normal, Vl2 -> Normal);
		}

	        IRIT_PT_COPY(V2 -> Coord, Vl2 -> Coord);

		IP_SET_INTERNAL_VRTX(V2);
		IP_SET_INTERNAL_VRTX(Vl2 -> Pnext -> Pnext);
	    }
	    else if (D23 >= D12 && D23 >= D31) {
		IRIT_PT_BLEND(Vl2 -> Pnext -> Coord, V2 -> Coord, V3 -> Coord, 0.5);

		GMInterpVrtxRGBBetweenTwo(Vl2 -> Pnext, V2, V3);
		COPY_RGB_ATTR(V3, Vl2 -> Pnext)
		GMInterpVrtxUVBetweenTwo(Vl2 -> Pnext, V2, V3);
		COPY_UV_ATTR(V3, Vl2 -> Pnext);

		if (DoNormals) {
		    GMInterpVrtxNrmlBetweenTwo(Vl2 -> Pnext, V2, V3);
		    IRIT_VEC_COPY(V3 -> Normal, Vl2 -> Pnext -> Normal);
		}

	        IRIT_PT_COPY(V3 -> Coord, Vl2 -> Pnext -> Coord);

		IP_SET_INTERNAL_VRTX(V3);
		IP_SET_INTERNAL_VRTX(Vl2);
	    }
	    else {
		IRIT_PT_BLEND(Vl2 -> Coord, V3 -> Coord, V1 -> Coord, 0.5);

		GMInterpVrtxRGBBetweenTwo(Vl2, V1, V3);
		COPY_RGB_ATTR(V3, Vl2)
		GMInterpVrtxUVBetweenTwo(Vl2, V1, V3);
		COPY_UV_ATTR(V3, Vl2)

		if (DoNormals) {
		    GMInterpVrtxNrmlBetweenTwo(Vl2, V1, V3);
		    IRIT_VEC_COPY(V3 -> Normal, Vl2 -> Normal);
		}

	        IRIT_PT_COPY(V3 -> Coord, Vl2 -> Coord);

		IP_SET_INTERNAL_VRTX(V2);
		IP_SET_INTERNAL_VRTX(Vl2);
	    }
	}
	else {
	    Pl = Pl -> Pnext;	 /* This polygon is small enough - skip it. */
	}
    }

    return Pls;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Affine transform the given UV coordinates in polygonal object PObj, in   M
* place.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:     A polygonal object to affine transform the UV vals.            M
*   Scale:    UV scale factors.						     M
*   Trans:    UV translational factors.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMGenUVValsForPolys                                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMAffineTransUVVals                                                      M
*****************************************************************************/
void GMAffineTransUVVals(IPObjectStruct *PObj,
			 const IrtRType Scale[2],
			 const IrtRType Trans[2])
{
    IRIT_STATIC_DATA AttribNumType /* Persistent store to unique Attrs IDs. */
        AttrUVValsID = ATTRIB_NAME_BAD_NUMBER;
    IPPolygonStruct *Pl;

    IP_ATTR_INIT_UNIQUE_ID_NUM(AttrUVValsID, "uvvals");

    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        IPVertexStruct
	    *V = Pl -> PVertex;

	do {
	    float
		*UV = AttrGetUVAttrib2(V -> Attr, AttrUVValsID);

	    UV[0] = (float) (UV[0] * Scale[0] + Trans[0]);
	    UV[1] = (float) (UV[1] * Scale[1] + Trans[1]);

	    V = V -> Pnext;
	}
	while (V != NULL && V != Pl -> PVertex);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Generates UV values for polygonal geometry, based on the XY(Z) 	     M
* coordinates of the geometry. Will NOT overwrite existing "uvvals", if any. M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:     A polygonal object to update UV vals for.                      M
*   UTextureRepeat, VTextureRepeat, WTextureRepeat:  Repeat texture factors. M
*   HasXYZScale:  If TRUE, WTextureRepeat is also valid - use XYZ coords.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMAffineTransUVVals                                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMGenUVValsForPolys                                                      M
*****************************************************************************/
void GMGenUVValsForPolys(IPObjectStruct *PObj,
			 IrtRType UTextureRepeat,
			 IrtRType VTextureRepeat,
			 IrtRType WTextureRepeat,
			 int HasXYZScale)
{
    IRIT_STATIC_DATA AttribNumType /* Persistent store to unique Attrs IDs. */
        AttrUVValsID = ATTRIB_NAME_BAD_NUMBER;
    int i, IgnoreAxis;
    float *OldUV;
    GMBBBboxStruct
	BBox = *GMBBComputeBboxObject(PObj);
    IPPolygonStruct *Pl;

    IP_ATTR_INIT_UNIQUE_ID_NUM(AttrUVValsID, "uvvals");

    if (HasXYZScale) {
	/* For each polygon, use the best out of the three XYZ coordinates. */
	for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	    int UseAxis[3];
	    IrtRType UV[2];
	    IrtVecType TextScl;
	    IPVertexStruct
		*V = Pl -> PVertex,
		*V2 = IPGetLastVrtx(V);
	    GMBBBboxStruct
		*PlBBox = GMBBComputeOnePolyBbox(Pl);

	    UseAxis[0] = UseAxis[1] = UseAxis[2] = 0;
	    do {
		int SameX = IRIT_APX_EQ(V -> Coord[0], V2 -> Coord[0]),
		    SameY = IRIT_APX_EQ(V -> Coord[1], V2 -> Coord[1]),
		    SameZ = IRIT_APX_EQ(V -> Coord[2], V2 -> Coord[2]);

		/* We must use the specified axes if other two are same. */
		UseAxis[0] += SameY && SameZ;
		UseAxis[1] += SameX && SameZ;
		UseAxis[2] += SameX && SameY;

		V2 = V;
		V = V -> Pnext;
	    }
	    while (V != NULL && V != Pl -> PVertex);

	    for (i = IgnoreAxis = 0; i < 3; i++) {
		if (UseAxis[i] == 0) {
		    IgnoreAxis = i;
		    break;
		}
	    }

	    /* Try to find the minimal bbox size that is also not used. */
	    if (IgnoreAxis >= 0) {
	        for (i = IgnoreAxis + 1; i < 3; i++) {
		    if (UseAxis[i] == 0 &&
			PlBBox -> Max[i] - PlBBox -> Min[i] <
			PlBBox -> Max[IgnoreAxis] - PlBBox -> Min[IgnoreAxis])
		        IgnoreAxis = i;
		}
	    }
	    else {
	        for (i = 0; i < 3; i++) {
		    if (PlBBox -> Max[i] - PlBBox -> Min[i] <
			PlBBox -> Max[IgnoreAxis] - PlBBox -> Min[IgnoreAxis])
		        IgnoreAxis = i;
		}
	    }

	    IRIT_VEC_RESET(TextScl);
	    switch (IgnoreAxis) {
	        case 0:
	            TextScl[1] = VTextureRepeat / (BBox.Max[1] - BBox.Min[1]);
		    TextScl[2] = WTextureRepeat / (BBox.Max[2] - BBox.Min[2]);
		    break;
		case 1:
		    TextScl[0] = UTextureRepeat / (BBox.Max[0] - BBox.Min[0]);
		    TextScl[2] = WTextureRepeat / (BBox.Max[2] - BBox.Min[2]);
		    break;
		case 2:
		    TextScl[0] = UTextureRepeat / (BBox.Max[0] - BBox.Min[0]);
		    TextScl[1] = VTextureRepeat / (BBox.Max[1] - BBox.Min[1]);
		    break;
	        default:
		    assert(0);
		    break;
	    }

	    /* And set all vertices of the polygon with UV. */
	    V = Pl -> PVertex;
	    do {
		switch (IgnoreAxis) {
		    case 0:
		        UV[0] = (V -> Coord[1] - BBox.Min[1]) * TextScl[1];
			UV[1] = (V -> Coord[2] - BBox.Min[2]) * TextScl[2];
			break;
		    case 1:
			UV[0] = (V -> Coord[0] - BBox.Min[0]) * TextScl[0];
			UV[1] = (V -> Coord[2] - BBox.Min[2]) * TextScl[2];
			break;
		    case 2:
			UV[0] = (V -> Coord[0] - BBox.Min[0]) * TextScl[0];
			UV[1] = (V -> Coord[1] - BBox.Min[1]) * TextScl[1];
			break;
	            default:
			assert(0);
			UV[0] = UV[1] = 1.0;
			break;
		}
		if ((OldUV = AttrGetUVAttrib2(V -> Attr,
					      AttrUVValsID)) != NULL) {
		    /* Use old values if has them but appl the scale. */
		    UV[0] = OldUV[0] * UTextureRepeat;
		    UV[1] = OldUV[1] * VTextureRepeat;
		}

		AttrSetUVAttrib2(&V -> Attr, AttrUVValsID, UV[0], UV[1]);

		V = V -> Pnext;
	    }
	    while (V != NULL && V != Pl -> PVertex);
	}
    }
    else {
	/* Globally find best two axes to employ to compute the UV values. */
	IgnoreAxis = 0;

	for (i = 1; i < 3; i++) {
	    if (BBox.Max[i] - BBox.Min[i] <
		BBox.Max[IgnoreAxis] - BBox.Min[IgnoreAxis])
		IgnoreAxis = i;
	}

	/* Keep the difference in Max slot, and scale with repeat vals. */
	for (i = 0; i < 3; i++)
	    BBox.Max[i] -= BBox.Min[i];

	switch (IgnoreAxis) {
	    case 0:
	        BBox.Max[1] = UTextureRepeat / BBox.Max[1];
		BBox.Max[2] = VTextureRepeat / BBox.Max[2];
		break;
	    case 1:
		BBox.Max[0] = UTextureRepeat / BBox.Max[0];
		BBox.Max[2] = VTextureRepeat / BBox.Max[2];
		break;
	    case 2:
		BBox.Max[0] = UTextureRepeat / BBox.Max[0];
		BBox.Max[1] = VTextureRepeat / BBox.Max[1];
		break;
	    default:
		assert(0);
		BBox.Max[0] = BBox.Max[1] = 1.0;
		break;
	}

	/* And set all vertices of all polygons with UV. */
	for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	    IPVertexStruct *V;
	    IrtRType UV[2];

	    V = Pl -> PVertex;
	    do {
	        switch (IgnoreAxis) {
		    case 0:
		        UV[0] = (V -> Coord[1] - BBox.Min[1]) * BBox.Max[1];
			UV[1] = (V -> Coord[2] - BBox.Min[2]) * BBox.Max[2];
			break;
		    case 1:
			UV[0] = (V -> Coord[0] - BBox.Min[0]) * BBox.Max[0];
			UV[1] = (V -> Coord[2] - BBox.Min[2]) * BBox.Max[2];
			break;
		    case 2:
			UV[0] = (V -> Coord[0] - BBox.Min[0]) * BBox.Max[0];
			UV[1] = (V -> Coord[1] - BBox.Min[1]) * BBox.Max[1];
			break;
	            default:
			assert(0);
			UV[0] = UV[1] = 1.0;
			break;
		}
		if ((OldUV = AttrGetUVAttrib2(V -> Attr,
					      AttrUVValsID)) != NULL) {
		    /* Use old values if has them but appl the scale. */
		    UV[0] = OldUV[0] * UTextureRepeat;
		    UV[1] = OldUV[1] * VTextureRepeat;
		}

		AttrSetUVAttrib2(&V -> Attr, AttrUVValsID, UV[0], UV[1]);

		V = V -> Pnext;
	    }
	    while (V != NULL && V != Pl -> PVertex);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find the third edge in the T junction VAdj on edge (V, V -> Pnext).      *
*   We assume here that no edge will have nore than one T junction on it.    *
*                                                                            *
* PARAMETERS:                                                                *
*   V, Pl:      Edge (V, V -> Next) and its polygon that has VAdj as an      *
*               T junction.						     *
*   Vadj, PAdj: The adjacent edge and polygon, with VAdj being a T junction. *
*   Polys:      All polygons in this polygonal mesh.                         *
*   Pl2Ret:     Will return polygon that holds third edge in T junction.     *
*   Eps:        Tolerance of comparisons.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPVertexStruct: Detected vertex V2 of the third edge (V2, V2 -> Pnext)   *
*                   in the T junction, or NULl if error.                     *
*****************************************************************************/
static IPVertexStruct *GMFindthirdVertexinTJun(const IPVertexStruct *V,
					       const IPPolygonStruct *Pl,
					       const IPVertexStruct *VAdj,
					       const IPPolygonStruct *PAdj,
					       IPPolygonStruct *Polys,
					       IPPolygonStruct **Pl2Ret,
					       IrtRType Eps)
{
    const IPVertexStruct *VThird1, *VThird2,
        *VNext = V -> Pnext != NULL ? V -> Pnext : Pl -> PVertex,
	*VAdjNext = VAdj -> Pnext != NULL ? VAdj -> Pnext : PAdj -> PVertex;
    IPPolygonStruct *Pl2;

    /* Find out relations between given two edges (which contains which).   */
    /* Reorder the data so (V, VNext) contains (VAdj, VAdjNext).            */
    if (IRIT_PT_PT_DIST_SQR(V -> Coord, V -> Pnext -> Coord) <
	IRIT_PT_PT_DIST_SQR(VAdj -> Coord, VAdj -> Pnext -> Coord)) {
        IRIT_SWAP(const IPVertexStruct *, V, VAdj);
        IRIT_SWAP(const IPVertexStruct *, VNext, VAdjNext);
        IRIT_SWAP(const IPPolygonStruct *, Pl, PAdj);
    }

    /* So now (V, VNext) contains (VAdj, VAdjNext) - find the two vertices  */
    /* of the missing third edge.					    */
    if (IRIT_PT_APX_EQ_EPS(V -> Coord, VAdj -> Coord, Eps)) {
        VThird1 = VNext;
	VThird2 = VAdjNext;
    }
    else if (IRIT_PT_APX_EQ_EPS(VNext -> Coord, VAdj -> Coord, Eps)) {
        VThird1 = V;
	VThird2 = VAdjNext;
    }
    else if (IRIT_PT_APX_EQ_EPS(V -> Coord, VAdjNext -> Coord, Eps)) {
        VThird1 = VNext;
	VThird2 = VAdj;
    }
    else if (IRIT_PT_APX_EQ_EPS(VNext -> Coord, VAdjNext -> Coord, Eps)) {
        VThird1 = V;
	VThird2 = VAdj;
    }
    else {
        assert(0);
	return NULL;
    }

    /* Search for this edge in all polygons. */
    for (Pl2 = Polys; Pl2 != NULL; Pl2 = Pl2 -> Pnext) {
        IPVertexStruct *V2Next,
	    *V2 = Pl2 -> PVertex;

	do {
	    V2Next = V2 -> Pnext != NULL ? V2 -> Pnext : Pl2 -> PVertex;

	    if ((IRIT_PT_APX_EQ_EPS(V2 -> Coord, VThird1 -> Coord, Eps) &&
		 IRIT_PT_APX_EQ_EPS(V2Next -> Coord, VThird2 -> Coord, Eps)) ||
		(IRIT_PT_APX_EQ_EPS(V2 -> Coord, VThird2 -> Coord, Eps) &&
		 IRIT_PT_APX_EQ_EPS(V2Next -> Coord, VThird1 -> Coord, Eps))) {
	        *Pl2Ret = Pl2;
	        return V2;
	    }

	    V2 = V2Next;
	}
	while (V2 != NULL && V2 != Pl2 -> PVertex);
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Function to search for T junctions in given polygonal mesh Polys.  The   M
* call back function TJunCB is invoked pn every such T junction.             M
*   This functions offers a naive and non-optimal (N^2) solution.            M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:     Polygonal mesh to search for T junctions in.                   M
*   TJuncCB:  Call back function to invoke with every detected T junction.   M
*             If E0 = (V0, V0 -> Pnext) is adjacent to both                  M
*             E1 = (V1, V1 -> Pnext) and E2 = (V2, V2 -> Pnext) so E2 is a T M
*             junction vertex on edge E0, the TJuncCB will be invoked as     M
*             TJuncCB(V0, V1, V2, Pl0, Pl1, Pl2).			     M
*   Eps:      Tolerance of testing.					     *
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    Number of detected T junctions, zero if none.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMIsTJunction                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMIdentifyTJunctions                                                     M
*****************************************************************************/
int GMIdentifyTJunctions(IPObjectStruct *PObj,
			 GMIdentifyTJunctionFuncType TJuncCB,
			 IrtRType Eps)
{
    int n = 0;
    IPPolygonStruct *Pl;

    if (!IP_IS_POLY_OBJ(PObj) || !IP_IS_POLYGON_OBJ(PObj))
	return 0;

    BoolGenAdjacencies(PObj);

    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	IPVertexStruct *VAdj,
	    *V = Pl -> PVertex;

	do {
	    if (V -> PAdj != NULL &&
		(VAdj = GetAdjVertex(V, Pl, V -> PAdj)) != NULL) {
		if ((IRIT_PT_APX_EQ(V -> Coord, VAdj -> Coord) &&
		     IRIT_PT_APX_EQ(V -> Pnext -> Coord, VAdj -> Pnext -> Coord)) ||
		    (IRIT_PT_APX_EQ(V -> Coord, VAdj -> Pnext -> Coord) &&
		     IRIT_PT_APX_EQ(V -> Pnext -> Coord, VAdj -> Coord))) {
		    /* This edge is identical on both sides. */
		    V -> PAdj = VAdj -> PAdj = NULL;   /* Clear adjacency. */

		}
		else if (IRIT_PT_PT_DIST_SQR(V -> Coord, V -> Pnext -> Coord) >
			 IRIT_PT_PT_DIST_SQR(VAdj -> Coord,
					     VAdj -> Pnext -> Coord)) {
		    /* A T junctions where (V, VNext) is the larger edge    */
		    /* and (Vadj, VAdjNext) is a subset ofg that - fetch    */
		    /* third edge and invoke TJuncCB.			    */
		    IPPolygonStruct *Pl2;
		    IPVertexStruct
		        *V2 = GMFindthirdVertexinTJun(V, Pl, VAdj, V -> PAdj,
						      PObj -> U.Pl, &Pl2, Eps);

		    if (V2 == NULL) {
		        GEOM_FATAL_ERROR(GEOM_ERR_CMPLX_T_JUNC);
		    }
		    else
			TJuncCB(V, VAdj, V2, Pl, V -> PAdj, Pl2);

		    n++;
		    V -> PAdj = VAdj -> PAdj = NULL;    /* Clear adjacency. */
		}
	    }
	    V = V -> Pnext;
	}
	while (V != NULL && V != Pl -> PVertex);
    }

    return n;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given a triangle Pl, divide and refine it along the edges as set by      M
* Refij.  The result can be between one (no refinement) and four (all edges  M
* are defined) triangles that are substituted in place in Pl and before      M
* Pl -> pnext.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   Pl:                  Triagle to refine if necessary.                     M
*   DeformVrtxFctrFunc:  Function to evaluate the deformation amount factor  M
*                        of a given vertex.			             M
*   Ref12, Ref23, Ref31: Booleans to set which edge must be refined.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:          TRUE if triangle underwent refinement, FALSE otherwise.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMRefineDeformedTriangle                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMRefineDeformedTriangle2                                                M
*****************************************************************************/
int GMRefineDeformedTriangle2(IPPolygonStruct *Pl,
			      GMPointDeformVrtxFctrFuncType DeformVrtxFctrFunc,
			      IrtBType Ref12,
			      IrtBType Ref23,
			      IrtBType Ref31)
{
    int i;
    IPVertexStruct *V12, *V23, *V31,
        *V1 = Pl -> PVertex,
        *V2 = V1 -> Pnext,
        *V3 = V2 -> Pnext;
    IrtBType
        VertexLoop = V3 -> Pnext == V1;
    IPPolygonStruct
        *PlsNew = NULL;

    if ((!Ref12 && !Ref23 && !Ref31) ||
	GMCollinear3Pts(V1 -> Coord, V2 -> Coord, V3 -> Coord))
        return FALSE;

    assert(V3 -> Pnext == NULL || V3 -> Pnext == V1);   /* Triangle please. */
    V3 -> Pnext = V1;                         /* Make vertex list circular. */

    /* Refine as necessary: */
    i = (Ref12 != FALSE) + (Ref23 != FALSE) + (Ref31 != FALSE);
    switch (i) {
        case 1:			     /* Refine one edge into two triangles. */
	    /* Make the refined edge V1V2:                                  */
	    if (Ref12) {
	    }
	    else if (Ref23) {
	        V1 = V2;
	    }
	    else { /* Ref31 */
	        V1 = V3;
	    }
	    Pl -> PVertex = V1;
	    V2 = V1 -> Pnext;
	    V3 = V2 -> Pnext;

	    /* Refine at V1V2 into two triangles.                           */
	    /*              * V3                                            */
	    /*             /|\                                              */
	    /*            / | \                                             */
	    /*           /  |  \                                            */
	    /*          /   |   \                                           */
	    /*      V1 *----*----* V2                                       */
	    /*             V12                                              */
	    V12 = IPAllocVertex2(V2);
	    GM_REF_UPDATE_INTRR_VRTX(V12, V1, V2);
	    V3 -> Pnext = VertexLoop ? V12 : NULL;
	    PlsNew = IPAllocPolygon(0, V12, Pl -> Pnext);      /* V12 V2 V3 */
	    IRIT_PLANE_COPY(PlsNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlsNew);

	    V12 = IPCopyVertex(V12);
	    V3 = IPCopyVertex(V3);
	    V1 -> Pnext = V12;
	    V12 -> Pnext = V3;
	    V3 -> Pnext = VertexLoop ? V1 : NULL;
	    assert(Pl -> PVertex == V1);                       /* V1 V12 V3 */
	    break;
        case 2:			   /* Refine two edges into tree triangles. */
	    /* Make the non-refined edge V3V1: */
	    if (!Ref31) {
	    }
	    else if (!Ref23) {
	        V1 = V3;
	    }
	    else { /* !Ref12 */
	        V1 = V2;
	    }
	    Pl -> PVertex = V1;
	    V2 = V1 -> Pnext;
	    V3 = V2 -> Pnext;

	    /* Refine at V1V2 and V2V3 into three triangles:                */
	    /*              * V3                                            */
	    /*             /|\                                              */
	    /*            / | \                                             */
	    /*           /  |  * V23                                        */
	    /*          /   | / \                                           */
	    /*         /    |/   \                                          */
	    /*     V1 *-----*-----* V2                                      */
	    /*             V12                                              */
	    V12 = IPAllocVertex2(V2);
	    GM_REF_UPDATE_INTRR_VRTX(V12, V1, V2);
	    V23 = IPAllocVertex2(VertexLoop ? V12 : NULL);
	    GM_REF_UPDATE_INTRR_VRTX(V23, V2, V3);
	    V2 -> Pnext = V23;
	    PlsNew = IPAllocPolygon(0, V12, Pl -> Pnext);     /* V12 V2 V23 */
	    IRIT_PLANE_COPY(PlsNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlsNew);

	    V23 = IPCopyVertex(V23);
	    V12 = IPCopyVertex(V12);
	    V23 -> Pnext = V3;
	    V12 -> Pnext = V23;
	    V3 -> Pnext = VertexLoop ? V12 : NULL;
	    PlsNew = IPAllocPolygon(0, V12, PlsNew);          /* V12 V23 V3 */
	    IRIT_PLANE_COPY(PlsNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlsNew);

	    V3 = IPCopyVertex(V3);
	    V12 = IPCopyVertex(V12);
	    V3 -> Pnext = VertexLoop ? V1 : NULL;
	    V12 -> Pnext = V3;
	    V1 -> Pnext = V12;
	    assert(Pl -> PVertex == V1);                 /* V1 V12 V3 in Pl */
	    break;
        case 3:
	    /* Refine all edges into four triangles:                        */
	    /*              * V3                                            */
	    /*             / \                                              */
	    /*            /   \                                             */
	    /*       V31 *-----* V23                                        */
	    /*          / \   / \                                           */
	    /*         /   \ /   \                                          */
	    /*     V1 *-----*-----* V2                                      */
	    /*             V12                                              */
	    V31 = IPAllocVertex2(NULL);
	    GM_REF_UPDATE_INTRR_VRTX(V31, V3, V1);
	    V23 = IPAllocVertex2(V31);
	    GM_REF_UPDATE_INTRR_VRTX(V23, V2, V3);
	    V12 = IPAllocVertex2(V23);
	    GM_REF_UPDATE_INTRR_VRTX(V12, V1, V2);
	    V31 -> Pnext = VertexLoop ? V12 : NULL;
	    PlsNew = IPAllocPolygon(0, V12, Pl -> Pnext);    /* V12 V23 V31 */
	    IRIT_PLANE_COPY(PlsNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlsNew);

	    V12 = IPCopyVertex(V12);
	    V23 = IPCopyVertex(V23);
	    V12 -> Pnext = V2;
	    V2 -> Pnext = V23;
	    V23 -> Pnext = VertexLoop ? V12 : NULL;
	    PlsNew = IPAllocPolygon(0, V12, PlsNew);          /* V12 V2 V23 */
	    IRIT_PLANE_COPY(PlsNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlsNew);

	    V31 = IPCopyVertex(V31);
	    V23 = IPCopyVertex(V23);
	    V31 -> Pnext = V23;
	    V23 -> Pnext = V3;
	    V3 -> Pnext = VertexLoop ? V31 : NULL;
	    PlsNew = IPAllocPolygon(0, V31, PlsNew);          /* V31 V23 V3 */
	    IRIT_PLANE_COPY(PlsNew -> Plane, Pl -> Plane);
	    IP_SET_PLANE_POLY(PlsNew);

	    V12 = IPCopyVertex(V12);
	    V31 = IPCopyVertex(V31);
	    V1 -> Pnext = V12;
	    V12 -> Pnext = V31;
	    V31 -> Pnext = VertexLoop ? V1 : NULL;
	    assert(Pl -> PVertex == V1);                      /* V12 V31 V1 */
	    break;
        default:
	    assert(0);
    }
    if (PlsNew != NULL)
        Pl -> Pnext = PlsNew;

#ifdef DEBUG
    {
        int j = i;

	PlsNew = Pl;
	while (j-- >= 0) {
	    V1 = PlsNew -> PVertex;
	    assert(V1 -> Pnext -> Pnext -> Pnext == NULL ||
		   V1 -> Pnext -> Pnext -> Pnext == V1);/* Triangle please. */
	    PlsNew = PlsNew -> Pnext;
	}
    }
#endif /* DEBUG */

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given a triangle in some deformation function and a function to evaluate M
* deformation amount per vertex, refine the triangle as necessary to make    M
* the deformation be accurate within tolerance Tol.                          M
*   Only the vertices and mid-edge points on triangle are examined.          M
*   Newly created (refined) triangles are appended in place between Pl and   M
* Pl -> Pnext.  Note Pl will also be modified in place if refinement occurs. M
*                                                                            *
* PARAMETERS:                                                                M
*   Pl:                  Triangle to refine if necessary.                    M
*   DeformVrtxFctrFunc:  Function to evaluate the deformation amount factor  M
*                        of a given vertex.			             M
*   DeformVrtxDirFunc:   Function to evaluate the deformation vector         M
*                        (direction and amount) of a given vertex.           M
*   DeviationTol:        Of deformation approximation.			     M
*   MaxEdgeLen:          to allow in a triangle.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:          TRUE if triangle underwent refinement, FALSE otherwise.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMRefineDeformedTriangle2                                                M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMRefineDeformedTriangle                                                 M
*****************************************************************************/
int GMRefineDeformedTriangle(IPPolygonStruct *Pl,
			     GMPointDeformVrtxFctrFuncType DeformVrtxFctrFunc,
			     GMPointDeformVrtxDirFuncType DeformVrtxDirFunc,
			     IrtRType DeviationTol,
			     IrtRType MaxEdgeLen)
{
    IrtBType Ref12, Ref23, Ref31;
    IrtRType *R, *R1, *R2, *R3;
    IrtPtType Pt1, Pt2,Pt3, Pt12, Pt23, Pt31;
    IrtVecType V;
    IPVertexStruct V12, V23, V31,
        *V1 = Pl -> PVertex,
        *V2 = V1 -> Pnext,
        *V3 = V2 -> Pnext;

    assert(V3 -> Pnext == NULL || V3 -> Pnext == V1);   /* Triangle please. */

    /* Compute expected deformation of the three vertices of the triangle. */
    if ((R1 = DeformVrtxDirFunc(V1)) != NULL)
        IRIT_PT_ADD(Pt1, V1 -> Coord, R1)
    else
        IRIT_PT_COPY(Pt1, V1 -> Coord);

    if ((R2 = DeformVrtxDirFunc(V2)) != NULL)
        IRIT_PT_ADD(Pt2, V2 -> Coord, R2)
    else
        IRIT_PT_COPY(Pt2, V2 -> Coord);

    if ((R3 = DeformVrtxDirFunc(V3)) != NULL)
        IRIT_PT_ADD(Pt3, V3 -> Coord, R3)
    else
        IRIT_PT_COPY(Pt3, V3 -> Coord);

    /* Compute expected deformation of the mid-point edges of the triangle. */
    IRIT_ZAP_MEM(&V12, sizeof(IPVertexStruct));
    GM_REF_UPDATE_INTRR_VRTX(&V12, V1, V2);

    IRIT_ZAP_MEM(&V23, sizeof(IPVertexStruct));
    GM_REF_UPDATE_INTRR_VRTX(&V23, V2, V3);

    IRIT_ZAP_MEM(&V31, sizeof(IPVertexStruct));
    GM_REF_UPDATE_INTRR_VRTX(&V31, V3, V1);

    if (R1 != NULL || R2 != NULL) {
        if ((R = DeformVrtxDirFunc(&V12)) != NULL)
	    IRIT_PT_ADD(Pt12, V12.Coord, R)
	else
	    IRIT_PT_COPY(Pt12, V12.Coord);
	IRIT_VEC_SUB(V, Pt2, Pt1);
	Ref12 = DeviationTol < GMDistPointLine(Pt12, Pt1, V);
    }
    else
        Ref12 = FALSE;

    if (R2 != NULL || R3 != NULL) {
        if ((R = DeformVrtxDirFunc(&V23)) != NULL)
	    IRIT_PT_ADD(Pt23, V23.Coord, R)
	else
	    IRIT_PT_COPY(Pt23, V23.Coord);
	IRIT_VEC_SUB(V, Pt3, Pt2);
	Ref23 = DeviationTol < GMDistPointLine(Pt23, Pt2, V);
    }
    else
        Ref23 = FALSE;

    if (R3 != NULL || R1 != NULL) {
        if ((R = DeformVrtxDirFunc(&V31)) != NULL)
	    IRIT_PT_ADD(Pt31, V31.Coord, R)
	else
	    IRIT_PT_COPY(Pt31, V31.Coord);
	IRIT_VEC_SUB(V, Pt1, Pt3);
	Ref31 = DeviationTol < GMDistPointLine(Pt31, Pt3, V);
    }
    else
        Ref31 = FALSE;

    /* Examine maximum edge lengths constraints: */
    if (IRIT_PT_PT_DIST(Pt1, Pt2) > MaxEdgeLen)
	Ref12 = TRUE;
    if (IRIT_PT_PT_DIST(Pt2, Pt3) > MaxEdgeLen)
	Ref23 = TRUE;
    if (IRIT_PT_PT_DIST(Pt3, Pt1) > MaxEdgeLen)
	Ref31 = TRUE;

    /* Compute the distances between the mid-points and the edges. */
    return GMRefineDeformedTriangle2(Pl, DeformVrtxFctrFunc,
				     Ref12, Ref23, Ref31);
}


