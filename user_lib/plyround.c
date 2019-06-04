/******************************************************************************
* PlyRound.c - rounding abilities to sharp edges in polygonal meshes.	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Jan. 15.					      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/bool_lib.h"

#define USER_MR_PERTURB_VAL	1.000301060

static int UserMRDistMeshVertices2Edge(IPPolygonStruct *Pl,
				       const IPPolygonStruct *EdgePlln);
static int UserMRFilterPolysByRadDist(IPPolygonStruct *Pl,
				     IrtRType RoundRadius,
				     IPPolygonStruct ***PlAtRad,
				     IPPolygonStruct ***PlBelowRad);
static IPPolygonStruct *UserMRExtractRulingAtRad(IPPolygonStruct **Pl1AtRad,
						 IrtRType RoundRadius,
						 IrtRType RulingExtent);
static int UserMRProjectVertices2RoundRadius(IPPolygonStruct **PlsVec,
					     const IPPolygonStruct *CenterPlln,
					     IrtRType RoundRadius,
					     IrtRType RoundShape);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Estimates the minimal distance of all given vertices in mesh Pl to       *
* polyline EdgePlln and store this minimal distance estimations as           *
* "_EdgeDist" attribute on the vertices.	                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Pl:       Polygonal mesh to project its vertices onto Crv and compute    *
*             their minimal distances to the curve.			     *
*   EdgePlln: Edge polyline to estimate the minimal distance of all vertices *
*             in Pl to.                                                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  TRUE if successful, FALSE otherwise.                               *
*****************************************************************************/
static int UserMRDistMeshVertices2Edge(IPPolygonStruct *Pl,
				       const IPPolygonStruct *EdgePlln)
{
    for ( ; Pl != NULL; Pl = Pl -> Pnext) {
        IPVertexStruct
	    *V = Pl -> PVertex;

	do {
	    CagdRType
	        EdgeDist = GMDistPolyPt2(EdgePlln, TRUE, V -> Coord,
					 NULL, FALSE);

	    assert(EdgeDist >= 0);

	    AttrSetRealAttrib(&V -> Attr, "_EdgeDist", EdgeDist);

	    V = V -> Pnext;
	}
	while (V != NULL && V != Pl -> PVertex);
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Searches input list Pl for polygons that are at distance RoundRadius or  *
* below it and return them in PlAtRad and PlBelowRad respectively.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Pls:         List of polygons to classify based on the "_EdgeDist"       *
*                attributes on their vertices if at RoundRadius or below it. *
*   RoundRadius: The rounding radius to seek and compare against.            *
*   PlsVecAtRad: A NULL terminated returned vector of all polygons that      *
*                are at distance RoundRadius.  Allocated here dynamically.   *
*   PlsVecBelowRad: A NULL terminated returned vector of all polygons that   *
*                closer than distance RoundRadius.  Allocated dynamically.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  TRUE if successful, FALSE otherwise.                               *
*****************************************************************************/
static int UserMRFilterPolysByRadDist(IPPolygonStruct *Pls,
				     IrtRType RoundRadius,
				     IPPolygonStruct ***PlsVecAtRad,
				     IPPolygonStruct ***PlsVecBelowRad)
{
    int AtRadIdx = 0,
        BelowRadIdx = 0,
        n = IPPolyListLen(Pls);
    IPPolygonStruct *Pl, **PlsAtRad, **PlsBelowRad;

    *PlsVecAtRad = PlsAtRad = IritMalloc(sizeof(IPPolygonStruct *) * (n + 1));
    *PlsVecBelowRad = PlsBelowRad =
                              IritMalloc(sizeof(IPPolygonStruct *) * (n + 1));

    for (Pl = Pls; Pl != NULL; Pl = Pl -> Pnext) {
        int HasMoreThanDist = FALSE,
	    HasLessThanDist = FALSE;
	IPVertexStruct
	    *V = Pl -> PVertex;

	do {
	    IrtRType
	        d = AttrGetRealAttrib(V -> Attr, "_EdgeDist");

	    assert(!IP_ATTR_IS_BAD_REAL(d));
	    if (d > RoundRadius)
	        HasMoreThanDist = TRUE;
	    if (d < RoundRadius)
	        HasLessThanDist = TRUE;

	    V = V -> Pnext;
	}
	while (V != Pl -> PVertex && V != NULL);

	if (HasLessThanDist) {
	    if (HasMoreThanDist) {
	        /* Has more & less than distance - place in AtRad vector. */
	        PlsAtRad[AtRadIdx++] = Pl;
	    }
	    else {
	        PlsBelowRad[BelowRadIdx++] = Pl;
	    }
	}
    }

    PlsAtRad[AtRadIdx] = NULL;
    PlsBelowRad[BelowRadIdx] = NULL;


#ifdef DEBUG_COLOR_POLYS
    {
        int i;

	for (Pl = Pls; Pl != NULL; Pl = Pl -> Pnext)
	    AttrSetRGBColor(&Pl -> Attr, 100, 255, 255);
	for (i = 0; i < AtRadIdx; i++)
	    AttrSetRGBColor(&PlsAtRad[i] -> Attr, 255, 255, 100);
	for (i = 0; i < BelowRadIdx; i++)
	    AttrSetRGBColor(&PlsBelowRad[i] -> Attr, 255, 100, 255);

	IPStdoutObject(IPGenPOLYObject(Pl), FALSE);
    }
#endif /* DEBUG_COLOR_POLYS */

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Extract a ruled surface from the polygons that are at distance           *
* RoundRadius.                                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   PlsVec:        NULL terminated vector of all the polygons that are at    *
*                  distance RoundRadius to compute a ruled surface from.     *
*   RoundRadius:   The rounding radius to seek and compare against.          *
*   RulingExtent:  The ruling extent to compute.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPPolygonStruct *:   The constructed ruled surface or NULL if error.     *
*****************************************************************************/
static IPPolygonStruct *UserMRExtractRulingAtRad(IPPolygonStruct **PlsVec,
						 IrtRType RoundRadius,
						 IrtRType RulingExtent)
{
    int i;
    IPPolygonStruct *Plln,*Pl,
        *RuledPls = NULL,
        *IsoLines = NULL;
    IPObjectStruct *PolyObj, *TriObj;

    for (i = 0; PlsVec[i] != NULL; i++) {
        IPVertexStruct *PrevV,
	    *NewLine = NULL,  /* Line crossing a triangle using 2 vertices. */
	    *V = PlsVec[i] -> PVertex,
	    *VLast = IPGetLastVrtx(V);
	IrtRType
	    PrevD = AttrGetRealAttrib(VLast -> Attr, "_EdgeDist");

	PrevV = VLast;
	if (VLast -> Pnext == NULL) {
	    /* Make vertex list cyclic temporarily */
	    VLast -> Pnext = V;
	}
	else
	    VLast = NULL;

	do {
	    IrtRType
	        d = AttrGetRealAttrib(V -> Attr, "_EdgeDist");

	    if ((PrevD - RoundRadius) * (d - RoundRadius) < 0.0) {
	        /* We found an edge crossing the RoundRadius distance. */
	        IrtRType
		    t = (RoundRadius - d) / (PrevD - d);

		NewLine = IPAllocVertex2(NewLine);
		IRIT_PT_BLEND(NewLine -> Coord, PrevV -> Coord, V -> Coord, t);
		if (IP_HAS_NORMAL_VRTX(V) && IP_HAS_NORMAL_VRTX(PrevV)) {
		    IRIT_VEC_BLEND(NewLine -> Normal,
				   PrevV -> Normal, V -> Normal, t);
		    IRIT_VEC_NORMALIZE(NewLine -> Normal);
		}
		else {
		    if (!IP_HAS_PLANE_POLY(PlsVec[i])) {
		        assert(0);                    /* Should not happen. */
			IPUpdatePolyPlane(PlsVec[i]);
		    }

		    IRIT_VEC_COPY(NewLine -> Normal, PlsVec[i] -> Plane);
		}
		IP_SET_NORMAL_VRTX(NewLine);
	    }

	    PrevD = d;
	    PrevV = V;
	    V = V -> Pnext;
	}
	while (V != PlsVec[i] -> PVertex);

	if (VLast != NULL) /* Recover non-circular linked list of vertices. */
	    VLast -> Pnext = NULL;

	if (NewLine != NULL &&		        /* Have exactly 2 vertices. */
	    NewLine -> Pnext != NULL &&
	    NewLine -> Pnext -> Pnext == NULL) {
	    /* Found a line segment at ROundRadius distance - add to pool. */
	    IsoLines = IPAllocPolygon(0, NewLine, IsoLines);
	}
	else if (NewLine != NULL)
	    IPFreeVertex(NewLine);
    }

    if (IsoLines == NULL)
        return NULL;

    IsoLines = GMMergePolylines(IsoLines, IRIT_EPS);

#   ifdef DEBUG_USER_MR_DUMP_ISO_LINES
	IPStdoutObject(IPGenPOLYLINEObject(IsoLines), FALSE);
#   endif /* DEBUG_USER_MR_DUMP_ISO_LINES */

    /* Building the rulings. */
    for (Plln = IsoLines; Plln != NULL; Plln = Plln -> Pnext) {
        IPVertexStruct *V;

	for (V = Plln -> PVertex;
	     V != NULL && V -> Pnext != NULL;
	     V = V -> Pnext) {
	    int VrtcsRvrsd;
	    IrtVecType V1, V2, V3, V4, N1, N2;
	
	    IRIT_VEC_SCALE2(N1, V -> Normal, RulingExtent);
	    IRIT_VEC_SCALE2(N2, V -> Pnext -> Normal, RulingExtent);

	    IRIT_VEC_ADD(V1, V -> Coord, N1);
	    IRIT_VEC_SUB(V2, V -> Coord, N1);
	    IRIT_VEC_SUB(V3, V -> Pnext -> Coord, N2);
	    IRIT_VEC_ADD(V4, V -> Pnext -> Coord, N2);

	    RuledPls = PrimGenPolygon4Vrtx(V1, V2, V3, V4, NULL, &VrtcsRvrsd,
					   RuledPls);
	}
    }

    IPFreePolygonList(IsoLines);

    /* Convert to triangles. */
    PolyObj = IPGenPOLYObject(RuledPls);
    TriObj = GMConvertPolysToTriangles(PolyObj);
    IPFreeObject(PolyObj);
    RuledPls = TriObj -> U.Pl;
    TriObj -> U.Pl = NULL;
    IPFreeObject(TriObj);

    /* Make sure all triangles has valid normal. */
    for (Pl = RuledPls; Pl != NULL; Pl = Pl -> Pnext) {
        IPUpdatePolyPlane(Pl);
    }

#   ifdef DEBUG_USER_MR_DUMP_RULINGS
	IPStdoutObject(IPGenPOLYObject(RuledPls), FALSE);
#   endif /* DEBUG_USER_MR_DUMP_ISO_RULINGS */

    return RuledPls;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Project vertices closer than RoundRadius to the edge curve towards the   *
* center curve so the vertices i approximately RoundRadius away from the     *
* center curve.                                                              *
*                                                                            *
* PARAMETERS:                                                                *
*   PlsVec:      A NULL terminated vector of polygons to update their        *
*                vertices.  						     *
*   CenterPlln:  Center polyline of the roundings to move vertices toward.   *
*   RoundRadius: The desired radius of the approximated blend.		     *
*                Holds only for meshes that meet at the shared edge          *
*                orthogonally.						     *
*   RoundShape:  Positive bias to affect the rounding Shape. 1.0 to have no  *
*                affect and values larger (smaller) than 1.0 to affect the   *
*                rounding shape.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:     TRUE if successful, FALSE otherwise.                            *
*****************************************************************************/
static int UserMRProjectVertices2RoundRadius(IPPolygonStruct **PlsVec,
					     const IPPolygonStruct *CenterPlln,
					     IrtRType RoundRadius,
					     IrtRType RoundShape)
{
    int i;

    for (i = 0; PlsVec[i] != NULL; i++) {
        IPVertexStruct
	    *V = PlsVec[i] -> PVertex;

	do {
	    IrtRType
	        EdgeDist = AttrGetRealAttrib(V -> Attr, "_EdgeDist");

	    /* If close to the (to be) rounded edge - compute the distance  */
	    /*  to the blend center curve and move vertex toward CenterCrv. */
	    if (EdgeDist < RoundRadius) {
		IrtPtType MinPt;
	        IrtRType
		    CenterDist = GMDistPolyPt2(CenterPlln, TRUE, V -> Coord,
					       MinPt, FALSE);

		if (CenterDist > RoundRadius) {
		    IrtRType
			t = RoundRadius / CenterDist;

		    if (RoundShape < IRIT_EPS)
			RoundShape = IRIT_EPS;
		    t = pow(t, RoundShape);

		    IRIT_PT_BLEND(V -> Coord, V -> Coord, MinPt, t);
		}
	    }

	    V = V -> Pnext;
	} while (V != NULL && V != PlsVec[i] -> PVertex);
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given two meshes, Pl1 and Pl2, sharing common boundary edge(s), Edge12,  M
* updated Pl1 and Pl2 in place and round them along Edge12.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Pl1, Pl2:    The two input meshes sharing edge Edge12 to round along.    M
*                in place.						     M
*   Edge12:      The common edge(s) to round along.  Can be a list of edges  M
*                to round around all of them.				     M
*   RoundRadius: The desired radius of the approximated blend.		     M
*   RoundShape:  Bias to affect the rounding size. 1.0 to have no affect,    M
*                and values larger (smaller) than 1.0 to enlarge (shrink)    M
*                the rounding size.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:  TRUE if successful, FALSE otherwise.                               M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMPolyMeshSmoothing                                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   User2PolyMeshRoundEdge                                                   M
*****************************************************************************/
int User2PolyMeshRoundEdge(IPPolygonStruct *Pl1,
			   IPPolygonStruct *Pl2,
			   const IPPolygonStruct *Edge12,
			   IrtRType RoundRadius,
			   IrtRType RoundShape)
{
    int OutputInterCurve;
    IrtRType RulingExtent;
    GMBBBboxStruct BBox, BBox1, BBox2;
    IPPolygonStruct **Pl1AtRad, **Pl1BelowRad, **Pl2AtRad, **Pl2BelowRad,
        *Pl1Ruling, *Pl2Ruling, *Edge12next;
    IPObjectStruct *PObj1, *PObj2, *PInter;

    RoundRadius *= USER_MR_PERTURB_VAL;          /* Use a general distance. */

    if (Edge12 -> Pnext != NULL) {
        /* While constant - change its Pnext temporarily as we must have a  */
        /* single edge.							    */
        Edge12next = Edge12 -> Pnext;
        ((IPPolygonStruct *) Edge12) -> Pnext = NULL;
    }
    else
        Edge12next = NULL;
    
    BBox1 = *GMBBComputePolyListBbox(Pl1);
    BBox2 = *GMBBComputePolyListBbox(Pl2);
    BBox = *GMBBMergeBbox(&BBox1, &BBox2);
    RulingExtent = IRIT_MAX(IRIT_MAX(BBox.Max[0] - BBox.Min[0],
				     BBox.Max[1] - BBox.Min[1]),
			    BBox.Max[2] - BBox.Min[2]);

    /* Traverse both meshes and compute the distance of the vertices from   */
    /* the given edge curve, by projecting the vertices on curve.	    */
    UserMRDistMeshVertices2Edge(Pl1, Edge12);
    UserMRDistMeshVertices2Edge(Pl2, Edge12);

    /* Find out all triangles that are at RoundRadius from EdgeCrv or less  */
    UserMRFilterPolysByRadDist(Pl1, RoundRadius, &Pl1AtRad, &Pl1BelowRad);
    UserMRFilterPolysByRadDist(Pl2, RoundRadius, &Pl2AtRad, &Pl2BelowRad);

    /* Extract a ruled surface orthogonal to the at-radius edge.            */
    Pl1Ruling = UserMRExtractRulingAtRad(Pl1AtRad, RoundRadius, RulingExtent);
    Pl2Ruling = UserMRExtractRulingAtRad(Pl2AtRad, RoundRadius, RulingExtent);

    /* Intersect the two ruling surfaces to find center edge of the blend.  */
    OutputInterCurve = BoolSetOutputInterCurve(TRUE);

    PObj1 = IPGenPOLYObject(Pl1Ruling);
    PObj2 = IPGenPOLYObject(Pl2Ruling);
    PInter = BooleanCUT(PObj1, PObj2);
    IPFreeObject(PObj1);
    IPFreeObject(PObj2);

    BoolSetOutputInterCurve(OutputInterCurve);

#   ifdef DEBUG_USER_MR_DUMP_CENTER_LINE
	IPStdoutObject(PInter, FALSE);
#   endif /* DEBUG_USER_MR_DUMP_CENTER_LINE */

    /* Project all vertices at RoundRadius or below it to distance of       */
    /* RoundRadius.	 Note we apply the change in place on input.        */
    UserMRProjectVertices2RoundRadius(Pl1AtRad, PInter -> U.Pl,
				      RoundRadius, RoundShape);
    UserMRProjectVertices2RoundRadius(Pl1BelowRad, PInter -> U.Pl,
				      RoundRadius, RoundShape);
    UserMRProjectVertices2RoundRadius(Pl2AtRad, PInter -> U.Pl,
				      RoundRadius, RoundShape);
    UserMRProjectVertices2RoundRadius(Pl2BelowRad, PInter -> U.Pl,
				      RoundRadius, RoundShape);

    IPFreeObject(PInter);

    IritFree(Pl1AtRad);
    IritFree(Pl1BelowRad);
    IritFree(Pl2AtRad);
    IritFree(Pl2BelowRad);

    if (Edge12next != NULL)
        ((IPPolygonStruct *) Edge12) -> Pnext = Edge12next;

    return TRUE;
}

#ifdef DEBUG_TEST_POLY_ROUND

void main(int argc, char **argv)
{
    int Handler;
    IrtRType
        Rad = 0.1;

    if (argc == 2 || argc == 3) {
        if (argc == 3) {
	    /* Fetch the radius. */
	    Rad = atof(argv[1]);
	    argc--;
	    argv++;
        }

	if ((Handler = IPOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IPObjectStruct *PObj;

	    IPSetFlattenObjects(FALSE);
	    PObj = IPGetObjects(Handler);

	    /* Done with file - close it. */
	    IPCloseStream(Handler, TRUE);

	    /* Process the geometry - compute the accumulated area. */
	    if (IP_IS_OLST_OBJ(PObj)) {
	        IPObjectStruct
		    *Mesh1 = IPListObjectGet(PObj, 0),
		    *Mesh2 = IPListObjectGet(PObj, 1),
		    *Edge12 = IPListObjectGet(PObj, 2);

		if (IP_IS_POLY_OBJ(Mesh1) && IP_IS_POLYGON_OBJ(Mesh1) &&
		    IP_IS_POLY_OBJ(Mesh2) && IP_IS_POLYGON_OBJ(Mesh2) &&
		    IP_IS_POLY_OBJ(Edge12) && IP_IS_POLYLINE_OBJ(Edge12)) {
		    fprintf(stderr, "Rounding the two meshes using Rad = %f\n",
			    Rad);

		    if (User2PolyMeshRoundEdge(Mesh1 -> U.Pl, Mesh2 -> U.Pl,
					       Edge12 -> U.Pl, Rad, 0)) {
		        IPStdoutObject(Mesh1, FALSE);
			IPStdoutObject(Mesh2, FALSE);
		    }
		}
		else
		    fprintf(stderr, "Read list object is not a polyGON/polygon/polyline.\n");
	    }
	    else
	        fprintf(stderr, "Read object is not a list object.\n");
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	}
    }
    else {
	fprintf(stderr, "Usage: PlyRound [Radius] geom.itd\n");
    }
}

#endif /* DEBUG_TEST_POLY_ROUND */
