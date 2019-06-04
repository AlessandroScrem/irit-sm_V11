/******************************************************************************
* geomsmth.c - functions to smooth poly data.				      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Daniel Ghosalker and Gershon Elber, July 2011.		      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/allocate.h"
#include "inc_irit/iritprsr.h"
#include "geom_loc.h"

static int GMIsInterLineLine2D(const IrtPtType A1, 
			       const IrtPtType A2, 
			       const IrtPtType B1, 
			       const IrtPtType B2,
			       IrtPtType InterPt);

static int GMPolymeshSmoothingMarkRadiusVertices(
					     const IPPolyVrtxIdxStruct *PVIdx,
					     unsigned char *IsMoveableVrtx,
					     IrtRType ExpansionRadius,
					     int Restricted);

/*****************************************************************************
* DESCRIPTION:                                                               M
*    Move designated Vertices (that are typically not on the boundary) of    M
* the polygons in PolyObj to new (averages of their 1-rings) positions,      M
* smoothing the shape of the mesh.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   PolyObj:     Polygonal object to smooth, in place.			     M
*   VerticesToRound: If not NULL, only vertices found in VerticesToRound are M
*                allowed to move and all other vertices are kept stationary. M
*                Overwrites AllowBndryRound status, if not NULL.	     M
*   AllowBndryRound:  If TRUE, allow the boundary to move. If FALSE,         M
*                boundary is constrained to be fixed.          		     M
*                   Can take affect only if VerticesToRound is NULL.         M
*   RoundingRadius:  If we have a restriction on the movable vertices        M
*                (either VerticesToRound is not NULL or AllowBndryRound is   M
*                FALSE) and RoundingRadius is positive, any other vertex     M
*                that is at a Euclidean distance of less than RoundingRadius M
*                from a restricted (AllowBndryRound = FALSE) or movable      M
*                (in VerticesToRound) vertex will also be tagged as such.    M
*   NumIters:    Number of times to perform this smoothing algorithm.	     M
*   BlendFactor: 1.0 to move the vertex all the way to average position.     M
*                Otherwise, (less than 1.0) to the proper fraction.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *: The smoothed out polygons, in place. Same as PolyObj.  M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMPolyMeshSmoothing							     M
*****************************************************************************/
IPObjectStruct *GMPolyMeshSmoothing(IPObjectStruct *PolyObj,
				    const IPPolygonStruct *VerticesToRound,
				    int AllowBndryRound,
				    IrtRType RoundingRadius,
				    int NumIters,
				    IrtRType BlendFactor)
{
    IPPolyVrtxIdxStruct
	*PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PolyObj, TRUE, 0);
    IPVertexStruct
	* const *Vertices = PVIdx -> Vertices;
    unsigned char
        *IsMovableVrtx = (unsigned char *) IritMalloc(sizeof(unsigned char) *
					                   PVIdx -> NumVrtcs);
    int i, j,
	**VPolys = PVIdx -> Polygons;
    IrtPtType
	*NewPoints = (IrtPtType *) IritMalloc(sizeof(IrtPtType) *
					                   PVIdx -> NumVrtcs);
    IPPolygonStruct *Pl;

    /* Update a vector of fixed/movable state for all vertices. */
    if (VerticesToRound != NULL) {
        const IPPolygonStruct *Pl;

        IRIT_ZAP_MEM(IsMovableVrtx, sizeof(unsigned char) *
		                                           PVIdx -> NumVrtcs);

	for (Pl = VerticesToRound; Pl != NULL; Pl = Pl -> Pnext) {
	    const IPVertexStruct
	        *V = Pl -> PVertex;

	    do {
	        for (i = 0; i < PVIdx -> NumVrtcs; i++) {
		    if (IRIT_PT_APX_EQ(Vertices[i] -> Coord, V -> Coord)) {
		        IsMovableVrtx[i] = TRUE;
			break;
		    }
		}

	        V = V -> Pnext;
	    }
	    while (V != NULL && V != Pl -> PVertex);
	}

	if (RoundingRadius > 0) {
	    /* Fix all vertices next to the boundary that should be kept    */
	    /* stationary as well.			       		    */
	    GMPolymeshSmoothingMarkRadiusVertices(PVIdx, IsMovableVrtx,
						  RoundingRadius, FALSE);
	}
    }
    else if (!AllowBndryRound) {
        /* Mark boundaries as not movable. */
	for (i = 0; i < PVIdx -> NumVrtcs; i++) {
	    IsMovableVrtx[i] = !IPCnvIsVertexBoundary(PVIdx, i);
	}

	if (RoundingRadius > 0) {
	    /* Fix all vertices next to the boundary that should be kept    */
	    /* stationary as well.			       		    */
	    GMPolymeshSmoothingMarkRadiusVertices(PVIdx, IsMovableVrtx,
						  RoundingRadius, TRUE);
	}
    }
    else {
	for (i = 0; i < PVIdx -> NumVrtcs; i++) {
	    /* All vertices are movable. */
	    IsMovableVrtx[i] = TRUE;
	}
    }

    /* Do the low pass filtering on the movable vertices. */
    for (j = 0; j < NumIters; j++) {
	for (i = 0; i < PVIdx -> NumVrtcs; i++) {
	    const int *NP;

	    if (IsMovableVrtx[i]) {
	        if ((NP = IPCnvPolyVrtxNeighbors(PVIdx, i, 1)) != NULL) {
		    IRIT_PT_COPY(NewPoints[i], PVIdx -> Vertices[i] -> Coord); 
		    GMComputeAverageVertex2(NP, PVIdx, NewPoints[i],
					    BlendFactor); 
		}
	    }
	}

	/* Copy back all smoothed data back into original poly mesh. */
	for (Pl = PolyObj -> U.Pl, i = 0; 
	     Pl != NULL; 
	     Pl = Pl -> Pnext, i++) {
	    IPVertexStruct *V, 
		*V1 = Pl -> PVertex; 
	    int j = 0; 

	    V = V1; 
	    do {
	        if (IsMovableVrtx[VPolys[i][j]]) {
		    IRIT_PT_COPY(V -> Coord, NewPoints[VPolys[i][j]]);
		}
		j++; 

		V = V -> Pnext; 
	    }
	    while (V != V1 && V != NULL); 
	}
    }

    /* Update polys' planes and clear vertices' normals. */
    for (Pl = PolyObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        IPVertexStruct *V, 
	    *V1 = Pl -> PVertex; 

	IPUpdatePolyPlane(Pl);
	V = V1; 
	do {
	    IP_RST_NORMAL_VRTX(V);

	    V = V -> Pnext; 
	}
	while (V != V1 && V != NULL); 
    }

    IritFree(IsMovableVrtx);
    IritFree(NewPoints); 
    IPPolyVrtxIdxFree(PVIdx); 

    return PolyObj;
}


/*****************************************************************************
* DESCRIPTION:                                                               *
*   Expand the restrictions on vertices by a ExpansionRadius in Euclidean    *
* space.   Can be either an expansion of fixed vertices (if restrict is      *
* TRUE) or expansion of allowed-to-move vertices (if Restrict is FALSE).     *
*                                                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   PVIdx:          The vertex indices structure we work with.               *
*   IsMovableVrtx:  Boolean vector with tags which vector is movable or not. *
*   ExpansionRadius: The expansion radius.                                   *
*   Restrict:       To expand further restricted (fixed) vertices if TRUE,   *
*                   or to expand further free-to-move vertices if FALSE.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:                                                                     *
*****************************************************************************/
static int GMPolymeshSmoothingMarkRadiusVertices(
					     const IPPolyVrtxIdxStruct *PVIdx,
					     unsigned char *IsMovableVrtx,
					     IrtRType ExpansionRadius,
					     int Restrict)
{
    int i, j;
    unsigned char
        *IsMovableVrtxExpanded = (unsigned char *) IritMalloc(sizeof(char) *
					                   PVIdx -> NumVrtcs);
    IrtRType
        ExpRadSqr = IRIT_SQR(ExpansionRadius);
    IPVertexStruct
	* const *Vertices = PVIdx -> Vertices;

    /* Keep original unchanged until the end so we do not propagate to all. */
    IRIT_GEN_COPY(IsMovableVrtxExpanded, IsMovableVrtx,
		  sizeof(char) * PVIdx -> NumVrtcs);

    if (Restrict) {
        /* Expand the restrictions to move to any vertex that is            */
        /* RoundRadius or less away from a restricted (fixed) vertex.       */
        for (i = 0; i < PVIdx -> NumVrtcs; i++) {
	    if (!IsMovableVrtx[i]) {
	        /* Find all vertices close to this vertex and fix as well.  */
	        for (j = 0; j < PVIdx -> NumVrtcs; j++) {
		    if (i == j ||
			!IsMovableVrtx[j] ||
			IRIT_PT_PT_DIST_SQR(Vertices[i] -> Coord,
					    Vertices[j] -> Coord) > ExpRadSqr)
		        continue;

		    /* A movable vertex next to a fixed vertex - fix it.    */
		    IsMovableVrtxExpanded[j] = FALSE;
		}
	    }
	}
    }
    else {
        /* Expand the allowed-to-move vertices to any vertex that is        */
        /* RoundRadius or less away from a movable vertex.		    */
        for (i = 0; i < PVIdx -> NumVrtcs; i++) {
	    if (IsMovableVrtx[i]) {
	        /* Find all vertices close to this vertex and make movable. */
	        for (j = 0; j < PVIdx -> NumVrtcs; j++) {
		    if (i == j ||
			IsMovableVrtx[j] ||
			IRIT_PT_PT_DIST_SQR(Vertices[i] -> Coord,
					    Vertices[j] -> Coord) > ExpRadSqr)
		        continue;

		    /* A fixed vertex next to movable vertex - make movable.*/
		    IsMovableVrtxExpanded[j] = TRUE;
		}
	    }
	}
    }

    /* Update the original movable vertices vector. */
    IRIT_GEN_COPY(IsMovableVrtx, IsMovableVrtxExpanded,
		  sizeof(char) * PVIdx -> NumVrtcs);
    IritFree(IsMovableVrtxExpanded);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the average location of the vertices in vertex indices NS and   M
* blend this average with CenterPoint with ratio BlendFactor (== CenterPoint M
* if 0).								     M
*                                                                            *
* PARAMETERS:                                                                M
*   NS:		  List of indices of vertices to average in PVIdx.           M
*   PVIdx:	  Vertex array data structure.				     M
*   CenterPoint:  Input center location into which average is to be blended. M
*   BlendFactor:  1.0 to move the vertex all the way to the average position M
*		  otherwise (less than 1.0) to the proper fraction.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if computed average is valid, FALSE otherwise.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMPolyCentroid, GMComputeAverageVertex                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMComputeAverageVertex2						     M
*****************************************************************************/
int GMComputeAverageVertex2(const int *NS, 
			    const IPPolyVrtxIdxStruct *PVIdx,
			    IrtPtType CenterPoint, 
			    IrtRType BlendFactor)
{
    int i,
	Count = 0; 
    IrtPtType AvgPt;

    IRIT_PT_RESET(AvgPt); 

    for (i = 0; NS[i] != -1; i++) {
	IPVertexStruct
	    *V = PVIdx -> Vertices[NS[i]];

        IRIT_PT_ADD(AvgPt, AvgPt, V -> Coord); 
	Count++; 
    }

    if (Count == 0)
        return FALSE;

    IRIT_PT_SCALE(AvgPt, 1.0 / Count);
    IRIT_PT_BLEND(CenterPoint, AvgPt, CenterPoint, BlendFactor);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the average location of the vertices in poly VS and blend this  M
* average with CenterPoint with ratio BlendFactor (== CenterPoint if 0).     M
*                                                                            *
* PARAMETERS:                                                                M
*   VS:           List of vertices of the Polygon.			     M
*   CenterPoint:  Input center location into which average is to be blended. M
*   BlendFactor:  1.0 to move the vertex all the way to the average position M
*                 otherwise (less than 1.0) to the proper fraction.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if computed average is valid, FALSE otherwise.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMPolyCentroid, GMComputeAverageVertex2                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMComputeAverageVertex						     M
*****************************************************************************/
int GMComputeAverageVertex(const IPVertexStruct *VS, 
			   IrtPtType CenterPoint, 
			   IrtRType BlendFactor)
{
    const IPVertexStruct
        *V = VS; 
    int Count = 0; 
    IrtPtType AvgPt;

    IRIT_PT_RESET(AvgPt); 

    do {
        IRIT_PT_ADD(AvgPt, AvgPt, V -> Coord); 
	Count++; 

	V = V -> Pnext; 
    }
    while (V != VS && V != NULL);  

    if (Count == 0)
        return FALSE;

    IRIT_PT_SCALE(AvgPt, 1.0 / Count);
    IRIT_PT_BLEND(CenterPoint, AvgPt, CenterPoint, BlendFactor);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Check if a line V1V2 and polygon VS interest, in 2D.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   VS:       Cyclic List of vertices of the Polygon.			     M
*   V1, V2:   The end points of the line.				     M
*   t:        The blending value from V1 to V2 where the intersection has    M
*             occurred.	 Can be NULL to ignore.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if line and polygon do interest, FALSE otherwise.	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMIsInterLinePolygon2D						     M
*****************************************************************************/
int GMIsInterLinePolygon2D(const IPVertexStruct *VS, 
			   const IrtPtType V1, 
			   const IrtPtType V2, 
			   IrtRType *t)
{
    const IPVertexStruct 
	*V = VS; 
 
    do {
	const IPVertexStruct
	    *VNext = V -> Pnext != NULL ? V -> Pnext : VS;

	if (GMIsInterLineLine2D(V -> Coord, VNext -> Coord, V1, V2, t)) {
	    assert(*t >= 0.0 && *t <= 1.0);
	    return TRUE;
	}

	V = V -> Pnext; 
    }
    while (V != NULL && V != VS); 

    return FALSE; 
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Check if the two lines A1A2 and B1B2 interest in the XY plane.	     M
*   End points intersections, up to IRIT_EPS, are ignored.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   A1, A2:   The two points of first line.				     M
*   B1, B2:   The two points of second line.				     M
*   t:        The blending value from V1 to V2 where the intersection has    M
*             occurred.	 Can be NULL to ignore.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if the two lines intersect, FALSE otherwise.		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMIsInterLineLine2D							     M
*****************************************************************************/
static int GMIsInterLineLine2D(const IrtPtType A1, 
			       const IrtPtType A2, 
			       const IrtPtType B1, 
			       const IrtPtType B2, 
			       IrtRType *t)
{
    IrtPtType Tmp1, Tmp2, A12, B12; 
    IrtRType T1, T2; 

    IRIT_PT2D_SUB(A12, A2, A1); 
    IRIT_PT2D_SUB(B12, B2, B1); 

    A12[2] = B12[2] = 0;

    if (GM2PointsFromLineLine(A1, A12, B1, B12, Tmp1, &T1, Tmp2, &T2)) {
	/* Check if the intersection is located in the actual segment. */
	if (t != NULL)
	    *t = T2;
	return T1 >= 0.0 && T1 <= 1.0 && T2 >= 0.0 && T2 <= 1.0;
    }

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes a point inside the kernel of the given polygon, if has any.     M
*                                                                            *
* PARAMETERS:                                                                M
*   VE:      Cyclic list of vertices of the Polygon.			     M
*   KrnlPt:  Computed interior kernel point.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if kernel is not empty, FALSE otherwise. 		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMFindPtInsidePolyKernel						     M
*****************************************************************************/
int GMFindPtInsidePolyKernel(const IPVertexStruct *VE, IrtPtType KrnlPt)
{
    IrtRType 
	TotalWeight = 0.0; 
    const IPVertexStruct
        *V = VE; 
    IrtVecType Nrml; 

    GMFindUnConvexPolygonNormal(VE, Nrml);

    IRIT_PT_RESET(KrnlPt); 

    do {
	IrtPtType Vec1;

	IRIT_PT_SUB(Vec1, V -> Pnext -> Coord, V -> Coord); 
	if (!IRIT_PT_APX_EQ_ZERO_EPS(Vec1, IRIT_UEPS)) {
	    const IPVertexStruct
	        *V2 = V -> Pnext; 
	    IrtRType
	        MinT = -IRIT_INFNTY, 
	        MaxT = IRIT_INFNTY; 

	    /* Clip edge (V, V -> Pnext) againt all other edges of VE. */
	    do {
		IrtPtType Vec2, Tmp1, Tmp2; 
		IrtRType T1, T2; 

		IRIT_PT_SUB(Vec2, V2 -> Pnext -> Coord, V2 -> Coord); 
		if (!IRIT_PT_APX_EQ_ZERO_EPS(Vec2, IRIT_UEPS)) {
		    if (GM2PointsFromLineLine(V -> Coord, Vec1, 
					      V2 -> Coord, Vec2, 
					      Tmp1, &T1, Tmp2, &T2)) {
 		        IrtRType DP;
			IrtPtType Cross;

			/* update the end of the edge of the Kernel. */
			GMVecCrossProd(Cross, Vec1, Vec2); 
			if ((DP = GMVecDotProd(Cross, Nrml)) > -IRIT_UEPS)
			    MaxT = IRIT_MIN(MaxT, T1); 
			else if (DP < IRIT_UEPS)
			    MinT = IRIT_MAX(MinT, T1); 
		    }
		    else if (!GMCollinear3Pts(V -> Coord, 
					      V2 -> Coord, 
					      V2 -> Pnext -> Coord)) {
		        IrtPtType Cross, Vec3; 

			IRIT_PT_SUB(Vec3, V2 -> Coord, V -> Coord); 
			GMVecCrossProd(Cross, Vec1, Vec3);

			/* This edge not in the Kernel */
			if (GMVecDotProd(Cross, Nrml) < IRIT_UEPS)
			    MaxT = MinT - 1.0; 
		    }
		    else if (V > V2)
			/* This edge appears twice so take it off. */
			MaxT = MinT - 1.0; 
		}
		V2 = V2 -> Pnext; 
	    }
	    while (V2 != V && MaxT - MinT > IRIT_UEPS);

	    /* if the edge not empty accumulate its middle point, weighted  */
	    /* according to the length of the edge.			    */
	    if (MaxT - MinT > IRIT_UEPS) {
		IrtRType
		    Len = IRIT_PT_LENGTH(Vec1),
		    Weight = (MaxT - MinT) * Len;
		IrtPtType MidPt;

		TotalWeight += Weight;

		/* Compute middle point on kernel edge. */
	        IRIT_PT_SCALE(Vec1, (MinT + MaxT) * 0.5); 
		IRIT_PT_ADD(MidPt, Vec1, V -> Coord); 

		/* Accumulate the middle point, weighted. */
		IRIT_PT_SCALE(MidPt, Weight); 
	        IRIT_PT_ADD(KrnlPt, KrnlPt, MidPt); 
	    }
	}
	V = V -> Pnext;  
    }
    while (V != VE); 

    if (TotalWeight == 0.0)			       /* Kernel is empty. */
        return FALSE;

    IRIT_PT_SCALE(KrnlPt, 1.0 / TotalWeight); 

    return TRUE; 
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Finds the normal of polygon that can be non-convex.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   VL:    List of vertices of the Polygon.				     M
*   Nrml:  Computed normal of the polygon. Not normalized.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void				                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   GMFindUnConvexPolygonNormal						     M
*****************************************************************************/
void GMFindUnConvexPolygonNormal(const IPVertexStruct *VL, IrtVecType Nrml)
{
    const IPVertexStruct
        *V = VL; 

    IRIT_PT_RESET(Nrml); 

    do {
        IrtPtType Tmp, C1, C2; 

	IRIT_PT_SUB(C1, V -> Coord, VL -> Coord); 
	IRIT_PT_SUB(C2, V -> Pnext -> Coord, VL -> Coord); 
	IRIT_CROSS_PROD(Tmp, C1, C2); 
	IRIT_PT_ADD(Nrml, Nrml, Tmp); 
	V = V -> Pnext; 
    }
    while (V != VL && V != NULL);  
}

