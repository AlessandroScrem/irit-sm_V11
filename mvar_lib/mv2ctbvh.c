/******************************************************************************
* mv2ctbvh.c - 2contact motion analysis for planar curves.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yong-Joon Kim, July 2014.					      *
******************************************************************************/

#include "inc_irit/misc_lib.h"
#include "mvar_loc.h"

typedef struct Mvar2CtMindistStruct {
    CagdRType U, V, Mindist;
    CagdPType Foot, P;
    int Id;
    CagdBType Bnegative;
} Mvar2CtMindistStruct;

static CagdRType Mvar2CtCrvRadius(CagdCrvStruct *Crv);
static int Mvar2CtComputeConvexity(Mvar2CtBVHStruct *Bvh, CagdRType t);
static void Mvar2CtLineLineDistAux(CagdPType X, 
				   CagdPType Y,    
				   CagdPType P, 
				   CagdPType A,
				   CagdPType Q, 
				   CagdPType B);
static CagdRType Mvar2CtLinePointDist2(CagdPType P,
				       CagdPType L[2],
				       CagdPType Foot);
static void Mvar2CtMakeSphere(CagdCrvStruct *Crv, Mvar2CtSphereStruct *SBV);
static void Mvar2CtMakeLSC(CagdCrvStruct *Crv, Mvar2CtLineStruct *LBV);
static void Mvar2CtMakeCone(CagdCrvStruct *NCrv, Mvar2CtConeStruct *TCone);
static void Mvar2CtSetNode(CagdCrvStruct *Crv,
			   CagdCrvStruct *NCrv,
			   CagdCrvStruct *Curvature, 
			   Mvar2CtBVNodeStruct **Node);
static CagdPtStruct* Mvar2CtRemoveDuplicate(CagdPtStruct *Pts, 
					    CagdRType Tol);
static void Mvar2CtBuildHierarchy(CagdCrvStruct *Crv,
				  CagdCrvStruct *NCrv,
				  CagdCrvStruct *Curvature,
				  CagdPtStruct *Pts,
				  int Convexity,
				  Mvar2CtBVNodeStruct **Node,
				  Mvar2CtBVHStruct *Bvh,
				  CagdRType SubdivTol,
				  CagdRType BvTol);
static void Mvar2CtBVNodeFree(Mvar2CtBVNodeStruct *Node);
static void Mvar2CtTransformBV(Mvar2CtBVNodeStruct *Node, 
			       Mvar2CtLineStruct *TBV, 
			       Mvar2CtSphereStruct *TSBV, 
			       CagdRType XTrans, 
			       CagdRType YTrans, 
			       CagdRType Rot);
static void Mvar2CtAuxPtCrvMindist(CagdPType P, 
				   Mvar2CtBVNodeStruct *Node,  
				   Mvar2CtMindistStruct *Info);
static void Mvar2CtPtCrvMindist(CagdPType P, 
				Mvar2CtBVHStruct **BvhBs,
				int BSize, 
				Mvar2CtMindistStruct *Info);
static void Mvar2CtAuxPenetrationDepth(Mvar2CtBVNodeStruct *Node, 
				       Mvar2CtBVHStruct **BvhBs,
				       int BSize,
				       CagdRType Xtrans,
				       CagdRType Ytrans,
				       CagdRType Rot,
				       Mvar2CtMindistStruct *Info);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute the maximum radius of a given curve.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv: Curve to compute radius.                                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: Maximum radius of a given curve.                              *
*****************************************************************************/
static CagdRType Mvar2CtCrvRadius(CagdCrvStruct *Crv)
{
    CagdRType Length,
        Result = 0;
    CagdPType P;
    int i;

    for (i = 0; i < Crv -> Length; ++i) {
	CagdCoerceToE2(P, Crv -> Points, i, Crv -> PType);

	Length = IRIT_PT2D_LENGTH(P);

	if (Length > Result)
	    Result = Length;
    }

    return Result;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Compute a convexity at given parameter t.                                 *
*                                                                            *
* PARAMETERS:                                                                *
*   Bvh:    BVH for a curve.                                                 *
*   t:  Parameter value to evaluate convexity.		                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: 1 for convex -1 for concave.                                        *
*****************************************************************************/
static int Mvar2CtComputeConvexity(Mvar2CtBVHStruct *Bvh, CagdRType t)
{
    CagdRType *R;
    CagdPType D, DD;

    R = CagdCrvEval(Bvh -> DCrv, t);
    CagdCoerceToE2(D, &R, -1, Bvh -> DCrv -> PType);

    R = CagdCrvEval(Bvh -> DDCrv, t);
    CagdCoerceToE2(DD, &R, -1, Bvh -> DDCrv -> PType);

    return IRIT_CROSS_PROD_2D(D, DD) > 0 ? 1 : -1;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Auxiliary routine to compute distance between two line segments.          *
*                                                                            *
* PARAMETERS:                                                                *
*   X, Y: Closest point of two line segments will be placed herein.          *
*   A: Origin of first line segment.                                         *
*   DiffA: Difference vector of first line segment.                          *
*   B: Origin of second line segment.                                        *
*   DiffB: Difference vector of second line segment.                         *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtLineLineDistAux(CagdPType X, 
			           CagdPType Y,    
			           CagdPType A, 
			           CagdPType DiffA,
			           CagdPType B, 
			           CagdPType DiffB) 
{
    CagdRType t, u, Denom;
    CagdRType T[3], AdotA, BdotB, AdotB, AdotT, BdotT;
    CagdRType Tmp[3];

    IRIT_PT2D_SUB(T,B,A);
    AdotA = IRIT_DOT_PROD_2D(DiffA, DiffA);
    BdotB = IRIT_DOT_PROD_2D(DiffB, DiffB);
    AdotB = IRIT_DOT_PROD_2D(DiffA, DiffB);
    AdotT = IRIT_DOT_PROD_2D(DiffA, T);
    BdotT = IRIT_DOT_PROD_2D(DiffB, T);
    /* t parameterizes ray A, DiffA.                                       */  
    /* u parameterizes ray B, DiffB.                                       */ 
    /* Compute t for the closest point on ray A, DiffA to B, DiffB.        */
    Denom = AdotA * BdotB - AdotB * AdotB;
    /* Denom = a^2 b^2 sin theta. */
    t = (AdotT * BdotB - BdotT * AdotB) / Denom;
    /* Clamp result so t is on the segment A, DiffA. */
    if (t < 0) 
	t = 0; 
    else if (t > 1) 
	t = 1; 
    else if (!(t >= 0 && t <= 1)) 
	t = 0;
    /* Find u for point on ray B, DiffB closest to point at t. */
    u = (t * AdotB - BdotT) / BdotB;
    /* If u is on segment B, DiffB, t and u correspond to closest points.  */ 
    /* Otherwise, clamp u, recompute and clamp t.                          */ 
    if (u <= 0) {
	IRIT_PT2D_COPY(Y, B);
	t = AdotT / AdotA;
	if (t <= 0) 
	    IRIT_PT2D_COPY(X, A);
	else if (t >= 1) {
	    IRIT_PT2D_ADD(X, A, DiffA);
	}
	else if (!(t > 0 && t < 1)) {
	    IRIT_PT2D_COPY(X, A);
	}
	else {
	    IRIT_PT2D_SCALE2(X, DiffA, t);
	    IRIT_PT2D_ADD(X, X, A);
	    Tmp[0] = Tmp[1] = 0;
	    Tmp[2] = IRIT_CROSS_PROD_2D(T, DiffA);
	}
    }
    else if (u >= 1) {
	IRIT_PT2D_ADD(Y, B, DiffB);
	t = (AdotB + AdotT) / AdotA;
	if (t <= 0) 
	    IRIT_PT2D_COPY(X, A);
	else if (t >= 1) {
	    IRIT_PT2D_ADD(X, A, DiffA);
	}
	else if (!(t > 0 && t < 1)) {
	    IRIT_PT2D_COPY(X, A);
	}
	else {
	    IRIT_PT2D_SCALE2(X, DiffA, t);
	    IRIT_PT2D_ADD(X, X, A);
	    IRIT_PT2D_SUB(T, Y, A);
	    Tmp[0] = Tmp[1] = 0;
	    Tmp[2] = IRIT_CROSS_PROD_2D(T, DiffA);
	}
    }
    else if (!(u > 0 && u < 1)) {
	IRIT_PT2D_COPY(Y, B);
	t = AdotT / AdotA;
	
	if (t <= 0) 
	    IRIT_PT2D_COPY(X, A);
	else if (t >= 1) {
	    IRIT_PT2D_ADD(X, A, DiffA);
	}
	else if (!(t > 0 && t < 1)) {
	    IRIT_PT2D_COPY(X, A);
	}
	else {
	    IRIT_PT2D_SCALE2(X, DiffA, t);
	    IRIT_PT2D_ADD(X, X, A);
	    Tmp[0] = Tmp[1] = 0;
	    Tmp[2] = IRIT_CROSS_PROD_2D(T, DiffA);
	}
    }
    else {
	IRIT_PT2D_SCALE2(Y, DiffB, u);
	IRIT_PT2D_ADD(Y, Y, B);

	if (t <= 0) {
	    IRIT_PT2D_COPY(X, A);
	    Tmp[0] = Tmp[1] = 0;
	    Tmp[2] = IRIT_CROSS_PROD_2D(T, DiffA);
	}
	else if (t >= 1) {
	    IRIT_PT2D_ADD(X, A, DiffA);
	    IRIT_PT2D_SUB(T, B, X);
	    Tmp[0] = Tmp[1] = 0;
	    Tmp[2] = IRIT_CROSS_PROD_2D(T, DiffB);
	}
	else if (!(t > 0 && t < 1)) { 
	    IRIT_PT2D_COPY(X, A);
	    Tmp[0] = Tmp[1] = 0;
	    Tmp[2] = IRIT_CROSS_PROD_2D(T, DiffB);
	}
	else {
	    IRIT_PT2D_SCALE2(X, DiffA, t);
	    IRIT_PT2D_ADD(X, X, A);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Compute distance between two line segments.                               M
*                                                                            *
* PARAMETERS:                                                                M
*   A, B: Two line segments to compute distance.                             M 
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType: Distance between A and B.                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtLineLineDist                                                      M
*****************************************************************************/
CagdRType Mvar2CtLineLineDist(Mvar2CtLineStruct *A, Mvar2CtLineStruct *B)
{
    CagdPType X, Y;
    CagdVType V1, V2;

    IRIT_PT2D_SUB(V1, A -> P[1], A -> P[0]);
    IRIT_PT2D_SUB(V2, B -> P[1], B -> P[0]);

    Mvar2CtLineLineDistAux(X, Y, A -> P[0], V1, B -> P[0], V2);

    return IRIT_PT2D_DIST(X, Y);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate LSC for a given curve.                                           *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv: Curve to compute LSC.                                               *
*   LBV : LSC for a given curve.                                             * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeLSC(CagdCrvStruct *Crv, Mvar2CtLineStruct *LBV)
{
    CagdRType Dist;
    int i;
    
    CagdCoerceToE3(LBV -> P[0], Crv -> Points, 0, Crv -> PType);
    CagdCoerceToE3(LBV -> P[1], Crv -> Points, 
	           Crv -> Length - 1, Crv -> PType);

    LBV -> Epsilon = 0;

    for (i = 1; i < Crv -> Length - 1; ++i){
	CagdPType P;
	CagdCoerceToE3(P, Crv -> Points, i, Crv -> PType);

	Dist = Mvar2CtLinePointDist(P, LBV -> P);

	if (Dist > LBV -> Epsilon)
	    LBV -> Epsilon = Dist;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Compute distance between point and line segment.                          M
*                                                                            *
* PARAMETERS:                                                                M
*   P: Point to compute distance.                                            M
*   L: Two end points of line segment.                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType: Minimum distance between a line segment and a point.          M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtLinePointDist                                                     M
*****************************************************************************/
CagdRType Mvar2CtLinePointDist(CagdPType P, CagdPType L[2])
{
    CagdPType D, YmP0; 
    CagdRType t, DdD;

    IRIT_PT2D_SUB(D, L[1], L[0]);
    IRIT_PT2D_SUB(YmP0, P, L[0]);
    
    t = IRIT_DOT_PROD_2D(D, YmP0);

    if (t <= 0) {
	return sqrt(IRIT_DOT_PROD_2D(YmP0, YmP0));
    }

    DdD = IRIT_DOT_PROD_2D(D, D);

    if (t >= DdD) {
	CagdPType YmP1;

	IRIT_PT2D_SUB(YmP1, P, L[1]);
	return sqrt(IRIT_DOT_PROD_2D(YmP1, YmP1));
    }
    return sqrt(IRIT_MAX(0, IRIT_DOT_PROD_2D(YmP0, YmP0) - t * t / DdD));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Compute distance between point and line segment.                          *
*                                                                            *
* PARAMETERS:                                                                *
*   P: Point to compute distance.                                            *
*   L: Line segment.                                                         *
*   Foot: Foot point of P onto L[2] will be placed herein.                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: Minimum distance between a line segment and a point.          *
*****************************************************************************/
static CagdRType Mvar2CtLinePointDist2(CagdPType P, 
				       CagdPType L[2], 
				       CagdPType Foot)
{
    CagdPType D, YmP0; 
    CagdRType t, DdD;

    IRIT_PT2D_SUB(D, L[1], L[0]);
    IRIT_PT2D_SUB(YmP0, P, L[0]);

    t = IRIT_DOT_PROD_2D(D, YmP0);

    if (t <= 0) {
	IRIT_PT2D_COPY(Foot, L[0]);
	return sqrt(IRIT_DOT_PROD_2D(YmP0, YmP0));
    }

    DdD = IRIT_DOT_PROD_2D(D, D);

    if (t >= DdD) {
	CagdPType YmP1; 
	IRIT_PT2D_SUB(YmP1, P, L[1]);
	IRIT_PT2D_COPY(Foot, L[1]);
	return sqrt(IRIT_DOT_PROD_2D(YmP1, YmP1));
    }

    IRIT_PT2D_SCALE2(Foot, L[0], 1 - t);
    IRIT_PT2D_SCALE2(D, L[1], t);
    IRIT_PT2D_ADD(Foot, Foot, D);

    return sqrt(IRIT_DOT_PROD_2D(YmP0, YmP0) - t * t / DdD);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate a bounding sphere for a given curve.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv: Curve to compute bounding sphere.                                   *
*   SBV: Bounding sphere will be placed herein.                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeSphere(CagdCrvStruct *Crv, Mvar2CtSphereStruct *SBV)
{
    CagdRType TMin, TMax, *R, Dist;
    int i;

    CagdCrvDomain(Crv, &TMin, &TMax);

    R = CagdCrvEval(Crv, (TMin + TMax) * 0.5);

    CagdCoerceToE3(SBV -> Center, &R, -1, Crv -> PType);

    SBV -> Radius = 0;

    for (i = 1; i < Crv -> Length - 1; ++i) {
	CagdPType P;

	CagdCoerceToE3(P, Crv -> Points, i, Crv -> PType);

	Dist = IRIT_PT2D_DIST(P, SBV -> Center);

	if (Dist > SBV -> Radius)
	    SBV -> Radius = Dist;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate normal cone for NCrv.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   NCrv: Normal curve.                                                      *
*   TCone : Generated normal cone will be placed herein.                     * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeCone(CagdCrvStruct *NCrv, Mvar2CtConeStruct *TCone)
{
    CagdVType V;
    CagdRType *R, Min, Max, Dot, Angle;
    int i;

    CagdCrvDomain(NCrv, &Min, &Max);
    R = CagdCrvEval(NCrv, (Min + Max) * 0.5);
    CagdCoerceToE2(TCone -> Axis, &R, -1, NCrv -> PType);

    TCone -> Axis[2] = 0;
    IRIT_PT2D_NORMALIZE(TCone -> Axis);
    TCone -> Angle = 0;

    for (i = 0; i < NCrv -> Length; ++i) {	
	CagdCoerceToE2(V, NCrv -> Points, i, NCrv -> PType);
	
	if (IRIT_PT2D_LENGTH(V) < IRIT_PT_NORMALIZE_ZERO)
	    continue;

	IRIT_PT2D_NORMALIZE(V);	

	Dot = IRIT_BOUND(IRIT_DOT_PROD_2D(V, TCone -> Axis), -1, 1);
	Angle = acos(Dot);

	if (Angle > TCone -> Angle)
	    TCone -> Angle = Angle;
    }   
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate bounding volume and set necessary information for the bounding   *
*  volume.                                                                   *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv: Crv to bound.                                                       *
*   NCrv: Normal curve for Crv.                                              *
*   Curvature: Curvature square curve for Crv.                               *
*   Node : Generated bounding volume will be placed herein.                  * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtSetNode(CagdCrvStruct *Crv, 
		           CagdCrvStruct *NCrv, 
		           CagdCrvStruct *Curvature, 
		           Mvar2CtBVNodeStruct **Node)
{
    Mvar2CtBVNodeStruct *PNode;
    CagdRType Min, Max;
    int i;

    PNode = (Mvar2CtBVNodeStruct *) IritMalloc(sizeof(Mvar2CtBVNodeStruct));
    PNode -> Left = NULL;
    PNode -> Right = NULL;

    CagdCrvDomain(Crv, &Min, &Max);

    PNode -> Min = Min;
    PNode -> Max = Max;

    Mvar2CtMakeSphere(Crv, &(PNode -> SBV));/* Generate bounding sphere. */
    Mvar2CtMakeCone(NCrv, &(PNode -> NCone));/* Generate normal cone. */
    Mvar2CtMakeLSC(Crv, &(PNode -> LBV));/* Generate LSC. */

    CagdCoerceToE2(PNode -> NMin, NCrv -> Points, 0, NCrv -> PType);
    PNode -> NMin[2] = 0;
      
    CagdCoerceToE2(PNode -> NMax, NCrv -> Points, 
	           NCrv -> Length - 1, NCrv -> PType);   
    PNode -> NMax[2] = 0;

    IRIT_PT2D_NORMALIZE(PNode -> NMin);/* End normal directions. */
    IRIT_PT2D_NORMALIZE(PNode -> NMax);

    PNode -> KMin = IRIT_INFNTY;
    PNode -> KMax = 0;

    PNode -> Length = IRIT_MAX(IRIT_PT2D_LENGTH(PNode -> LBV.P[0]), 
	                       IRIT_PT2D_LENGTH(PNode -> LBV.P[1])) 
			       + PNode -> LBV.Epsilon;

    for (i = 0; i < Curvature -> Length; ++i){ /* Curvature min max. */
        if (Curvature -> Points[1][i] / Curvature -> Points[0][i] <
	                                                       PNode -> KMin)
	    PNode -> KMin = 
	             Curvature -> Points[1][i] / Curvature -> Points[0][i];

	if (Curvature -> Points[1][i] / Curvature -> Points[0][i] >
	                                                       PNode -> KMax)
	    PNode -> KMax = 
	            Curvature -> Points[1][i] / Curvature -> Points[0][i];
    }
    if (PNode -> KMin < 0)
	PNode -> KMin = 0;

    PNode -> KMin = sqrt(PNode -> KMin);
    PNode -> KMax = sqrt(PNode -> KMax);
    (*Node) = PNode;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Set id of bounding volume recursively.                                    M
*                                                                            *
* PARAMETERS:                                                                M
*   Node:    Bounding volume to set id.                                      M
*   Id:  id value to set.		                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtSetNodeId                                                         M
*****************************************************************************/
void Mvar2CtSetNodeId(Mvar2CtBVNodeStruct *Node, int Id)
{
    Node -> Id = Id;

    if (Node -> Left != NULL) {
	Mvar2CtSetNodeId(Node -> Left, Id);
	Mvar2CtSetNodeId(Node -> Right, Id);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Find smallest bounding volume includes domain range Min~Max.              M
*                                                                            *
* PARAMETERS:                                                                M
*   Node:    Data structure for rotation.                                    M
*   Min, Max:  Domain range for testing.		                     M
*   Parent: Place for saving the result.                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtGetParentBVNode                                                   M
*****************************************************************************/
void Mvar2CtGetParentBVNode(Mvar2CtBVNodeStruct *Node,
	       		    CagdRType Min,
			    CagdRType Max,
			    Mvar2CtBVNodeStruct **Parent)
{
    if (Node -> Min <= Min && Node -> Max >= Max) {		
	(*Parent) = Node;

	if (Node -> Left != NULL) {
	    Mvar2CtGetParentBVNode(Node -> Left, Min, Max, Parent);
	    Mvar2CtGetParentBVNode(Node -> Right, Min, Max, Parent);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Remove duplicates from the point list.                                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Pts:  BVH for a curve.                                                   *
*   tol:  Tolerance for removal.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdPtStruct*: Point list without duplicates.                            *
*****************************************************************************/
static CagdPtStruct* Mvar2CtRemoveDuplicate(CagdPtStruct *Pts,  
				            CagdRType tol)
{
    CagdPtStruct
        *TPt = Pts,
        *NewList = NULL;
    CagdRType
        Prev = -1;

    Pts = CagdPtsSortAxis(Pts, 1);
    
    while (Pts != NULL) {
	IRIT_LIST_POP(TPt, Pts);
	
	if (IRIT_FABS(TPt -> Pt[0] - Prev) > tol) {

	    IRIT_LIST_PUSH(TPt, NewList);
	    Prev = TPt -> Pt[0];
	}
	else
	    CagdPtFree(TPt);
    }
   
    return (CagdPtStruct*) CagdListReverse(NewList);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Build a bounding volume hierarchy for a curve.                            *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv: Crv to build a BVH.                                                 *
*   NCrv: Normal curve for Crv.                                              *
*   Curvature: Curvature square curve for Crv.                               *
*   Pts: Parameter values for subdivision.                                   *
*   Convexity: convexity of curve. 1 for convex -1 for concave and 0 for non *
*               convex or concave curve.                                     *
*   Node : Bounding volume will be placed herein.                            * 
*   Bvh: BVH data structure.                                                 *
*   SubdivTol: Subdivision tolerance for the construction.                   *
*   BvTol: Bounding volume error tolerance for the construction.             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtBuildHierarchy(CagdCrvStruct *Crv,
				  CagdCrvStruct *NCrv,
				  CagdCrvStruct *Curvature,
				  CagdPtStruct *Pts,
				  int Convexity,
				  Mvar2CtBVNodeStruct **Node,
				  Mvar2CtBVHStruct *Bvh,
				  CagdRType SubdivTol,
				  CagdRType BvTol)
{ 
    Mvar2CtBVNodeStruct *PNode;
    CagdRType 
	Param = 0;
    CagdPtStruct *Prev, *TPt,
	*Next = NULL;
    CagdCrvStruct *SubCrv, *SubNCrv, *SubCurvature;
    int ListLength; 
  
    /* Set the bounding volume and necessary information. */
    Mvar2CtSetNode(Crv, NCrv, Curvature, Node);

    PNode = *Node;
    ListLength = CagdListLength(Pts);

    if (ListLength > 0) {
	/* If there is parameter value for subdivision                     */ 
	/* subdivide it at given parameter.                                */ 
	int Count =0;
	for (TPt = Pts, Prev = NULL; TPt != NULL; TPt = TPt -> Pnext) {
	    /* Find middle point. */
	    if (Count++ == ListLength / 2) {
		Param = TPt -> Pt[0];
		Next = TPt -> Pnext;

		if (Prev != NULL)
		    Prev -> Pnext = NULL;
		else
		    Pts = NULL;

		CagdPtFree(TPt);
		break;
	    }
	    Prev = TPt;
	}
    
        PNode -> Convexity = Convexity;
	/* Subdivide curves. */
	SubCrv = CagdCrvSubdivAtParam(Crv, Param);
	SubNCrv = CagdCrvSubdivAtParam(NCrv, Param);
	SubCurvature = CagdCrvSubdivAtParam(Curvature, Param);
	
	Mvar2CtBuildHierarchy(SubCrv, SubNCrv, SubCurvature, Pts, Convexity,
	                      &(PNode -> Left), Bvh, SubdivTol, BvTol);
	
	Mvar2CtBuildHierarchy(SubCrv -> Pnext, SubNCrv -> Pnext, 
	                      SubCurvature -> Pnext, Next, Convexity, 
			      &(PNode -> Right), Bvh, SubdivTol, BvTol);

	CagdCrvFreeList(SubCrv);
	CagdCrvFreeList(SubNCrv);
	CagdCrvFreeList(SubCurvature);
    }
    else {
	/* If there is no parameter, the curve segment                     */ 
	/* is either convex or concave.                                    */
	if (Convexity == 0)
	    /* Set convexity. */
	    Convexity = 
	        Mvar2CtComputeConvexity(Bvh,
					(PNode -> Min + PNode -> Max) * 0.5);
	
	PNode -> Convexity = Convexity;

	Param = (PNode -> Min + PNode -> Max) * 0.5;

	if (PNode -> LBV.Epsilon > BvTol && 
	    PNode -> Max - PNode -> Min > SubdivTol) {
	    SubCrv = CagdCrvSubdivAtParam(Crv, Param);
	    SubNCrv = CagdCrvSubdivAtParam(NCrv, Param);
	    SubCurvature = CagdCrvSubdivAtParam(Curvature, Param);

	    Mvar2CtBuildHierarchy(SubCrv, SubNCrv, SubCurvature,
				  NULL, Convexity,
				  &(PNode -> Left), Bvh, SubdivTol, BvTol);

	    Mvar2CtBuildHierarchy(SubCrv -> Pnext, SubNCrv -> Pnext, 
				  SubCurvature -> Pnext, NULL, Convexity,  
				  &(PNode -> Right), Bvh, SubdivTol, BvTol);

	    CagdCrvFreeList(SubCrv);
	    CagdCrvFreeList(SubNCrv);
	    CagdCrvFreeList(SubCurvature);	
	}
    }  
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Build a bounding volume hierarchy for a curve.                            M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv: Crv to build a BVH.                                                 M
*   SubdivTol: Subdivision tolerance for the construction.                   M
*   BvTol: Bounding volume error tolerance for the construction.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   Mvar2CtBVHStruct *: Bounding volume hierarchy for a curve                M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtBuildBVH                                                          M
*****************************************************************************/
Mvar2CtBVHStruct *Mvar2CtBuildBVH(CagdCrvStruct *Crv, 
			          CagdRType SubdivTol, 
			          CagdRType BvTol)
{
    Mvar2CtBVHStruct *Bvh = (Mvar2CtBVHStruct *) 
	                                 IritMalloc(sizeof(Mvar2CtBVHStruct));
    CagdPtStruct *Pts;
    int i;

    Bvh -> Crv = CagdCrvCopy(Crv);
    Bvh -> DCrv = CagdCrvDerive(Crv);
    Bvh -> DDCrv = CagdCrvDerive(Bvh -> DCrv);
    Bvh -> NCrv = SymbCrv2DUnnormNormal(Crv);
    Bvh -> Curvature = SymbCrv2DCurvatureSqr(Crv);
    Bvh -> CurvContacts=NULL;
    Bvh -> Tol = BvTol;
    
    Pts = SymbCrv2DInflectionPts(Crv, 1e-8); /* Divide it at inflections. */

    for (i = Bvh -> Crv -> Order; i < Bvh -> Crv -> Length; ++i) {
	/* Divide it at knot values. */
	CagdPtStruct *Pt = CagdPtNew();
	Pt -> Pt[0] = Crv -> KnotVector[i];

	IRIT_LIST_PUSH(Pt, Pts);
    }

    Pts = Mvar2CtRemoveDuplicate(Pts, SubdivTol);
    Bvh -> Radius = Mvar2CtCrvRadius(Crv);/* Compute radius of moving curve. */
    /* Build a BVH. */
    Mvar2CtBuildHierarchy(Bvh -> Crv, Bvh -> NCrv, Bvh -> Curvature, 
	                  Pts, 0, &(Bvh -> Root), Bvh, SubdivTol, BvTol);  

    return Bvh;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Free bounding volume recursively.                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:  Bounding volume to deallocate.                                    *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtBVNodeFree(Mvar2CtBVNodeStruct *Node)
{
    if (Node -> Left != NULL) {
	Mvar2CtBVNodeFree(Node -> Left);
	Mvar2CtBVNodeFree(Node -> Right);
    }
   
    IritFree(Node);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Free bounding volume hierarchy.                                           M
*                                                                            *
* PARAMETERS:                                                                M
*   Bvh:  BVH to deallocate.                                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtFreeBVH                                                           M
*****************************************************************************/
void Mvar2CtFreeBVH(Mvar2CtBVHStruct *Bvh)
{
    CagdCrvFree(Bvh -> Crv);
    CagdCrvFree(Bvh -> DCrv);
    CagdCrvFree(Bvh -> DDCrv);
    CagdCrvFree(Bvh -> NCrv);
    CagdCrvFree(Bvh -> Curvature);
    MvarPtFreeList(Bvh -> CurvContacts);

    Mvar2CtBVNodeFree(Bvh -> Root);

    IritFree(Bvh);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Transform bounding volumes for a given transformation.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:  Bounding volume to transform.                                     *
*   TBV:  Place for saving transformed LSC.                                  *
*   TSBV:  Place for saving transformed bounding sphere.                     *
*   XTrans:  XTrans x translation.                                           *
*   YTrans:  YTrans y translation.                                           *
*   Rot: rotation angle in radian.                                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtTransformBV(Mvar2CtBVNodeStruct *Node, 
			       Mvar2CtLineStruct *TBV, 
			       Mvar2CtSphereStruct *TSBV, 
			       CagdRType XTrans, 
			       CagdRType YTrans, 
			       CagdRType Rot)
{
    TBV -> P[0][0] = Node -> LBV.P[0][0] * cos(Rot) - 
	             Node -> LBV.P[0][1] * sin(Rot) + XTrans;
    TBV -> P[0][1] = Node -> LBV.P[0][0] * sin(Rot) + 
	             Node -> LBV.P[0][1] * cos(Rot) + YTrans;

    TBV -> P[1][0] = Node -> LBV.P[1][0] * cos(Rot) - 
	             Node -> LBV.P[1][1] * sin(Rot) + XTrans;
    TBV -> P[1][1] = Node -> LBV.P[1][0] * sin(Rot) + 
	             Node -> LBV.P[1][1] * cos(Rot) + YTrans;

    TSBV -> Center[0] = Node -> SBV.Center[0] * cos(Rot) - 
	                Node -> SBV.Center[1] * sin(Rot) + XTrans;
    TSBV -> Center[1] = Node -> SBV.Center[0] * sin(Rot) + 
	                Node -> SBV.Center[1] * cos(Rot) + YTrans;

    TBV -> Epsilon = Node -> LBV.Epsilon;
    TSBV -> Radius = Node -> SBV.Radius;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Auxiliary function for computing signed distance from P to obstacle       *
*  curves.                                                                   *
*                                                                            *
* PARAMETERS:                                                                *
*   P:   Point to compute distance.                                   	     *
*   Node:  Bounding volume for obstacle curves.                              *
*   Info: Place for saving distance information.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtAuxPtCrvMindist(CagdPType P, 
			           Mvar2CtBVNodeStruct *Node,  
			           Mvar2CtMindistStruct *Info)
{
    if (Node -> Left == NULL) {
        /* We reached leaf node. */
	CagdPType Foot, PQ;
	/* Compute distance between point and LSC. */
	CagdRType Dist = Mvar2CtLinePointDist2(P, Node -> LBV.P, Foot)
	                                             - Node -> LBV.Epsilon;
	
	/* If distance is smaller than current mindist update it. */
	if (Dist < Info -> Mindist) {
	    if (Dist < 0)
	        Info -> Mindist = 0;
	    else
	        Info -> Mindist = Dist;

	    IRIT_PT2D_COPY(Info -> Foot, Foot);
	    IRIT_PT2D_SUB(PQ, P, Foot);

	    if (IRIT_DOT_PROD_2D(Node -> NCone.Axis, PQ) < 0)
		/* Check direction for signed distance. */
		Info -> Bnegative = TRUE;
	    else
		Info -> Bnegative = FALSE;
	}
    }
    else {
	CagdRType Dist1, Dist2;

	/* Compute distance between point and LSC. */
	Dist1 = Mvar2CtLinePointDist(P, Node -> Left -> LBV.P)
	                                      - Node -> Left -> LBV.Epsilon;
	Dist2 = Mvar2CtLinePointDist(P, Node -> Right -> LBV.P)
	                                     - Node -> Right -> LBV.Epsilon;
	if (Dist1 < Dist2) {
	    /* Does this node have a chance to improve current mindist? */
	    if (Dist1 < Info -> Mindist)
		Mvar2CtAuxPtCrvMindist(P, Node -> Left, Info);
	    if (Dist2 < Info -> Mindist)
		Mvar2CtAuxPtCrvMindist(P, Node -> Right, Info);
	}
	else {
	     /* Does this node have a chance to improve current mindist? */
	    if (Dist2 < fabs(Info -> Mindist))
		Mvar2CtAuxPtCrvMindist(P, Node -> Right, Info);

	    if (Dist1 < fabs(Info -> Mindist))
		Mvar2CtAuxPtCrvMindist(P, Node -> Left, Info);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Compute a signed distance from P to obstacle curves.                      *
*                                                                            *
* PARAMETERS:                                                                *
*   P:   Point to compute distance.                                   	     *
*   BvhBs:  Bvh for obstacle curves.                                         *
*   BSize:  Number of obstacle curves.                                       *
*   Info: Place for saving distance information.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*                                                                            *
*****************************************************************************/
static void Mvar2CtPtCrvMindist(CagdPType P, 
			        Mvar2CtBVHStruct **BvhBs,
			        int BSize, 
			        Mvar2CtMindistStruct *Info)
{
    int i;

    Info -> Mindist = IRIT_INFNTY;

    for (i = 0; i < BSize; ++i)
	Mvar2CtAuxPtCrvMindist(P, BvhBs[i] -> Root, Info);

    if (Info -> Bnegative)
	Info -> Mindist = - Info -> Mindist;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Auxiliary function for computing a maximum penetration depth of moving    *
*  curve into obstacle curves.                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhBs:  Bvh for obstacle curves.                                         *
*   BSize:  Number of obstacle curves.                                       *
*   Xtrans: X translation.                                                   *
*   Ytrans: Y translation.                                                   *
*   Rot: Rotation angle in radian.                                           *
*   Info: Place for saving penetration information.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtAuxPenetrationDepth(Mvar2CtBVNodeStruct *Node, 
				       Mvar2CtBVHStruct **BvhBs,
				       int BSize,
				       CagdRType Xtrans,
				       CagdRType Ytrans,
				       CagdRType Rot,
				       Mvar2CtMindistStruct *Info)
{
    IRIT_STATIC_DATA Mvar2CtLineStruct TLBV1, TLBV2;
    IRIT_STATIC_DATA Mvar2CtSphereStruct TSBV1, TSBV2;
    IRIT_STATIC_DATA Mvar2CtMindistStruct PenethInfo1, PenethInfo2;
    Mvar2CtBVNodeStruct *Left, *Right;
    CagdRType Dist1, Dist2;
    CagdBType
        Result1 = TRUE, 
        Result2 = TRUE;

    if (Node -> Left == NULL) {
	/* Transform bounding volumes. */
	Mvar2CtTransformBV(Node, &TLBV1, &TSBV1, Xtrans, Ytrans, Rot);
	/* Compute penetration depth of center of bounding sphere. */
	Mvar2CtPtCrvMindist(TSBV1.Center, BvhBs, BSize, &PenethInfo1);

	if (PenethInfo1.Mindist < Info -> Mindist) {
	    IRIT_PT_COPY(Info -> P, TSBV1.Center);
	    IRIT_PT_COPY(Info -> Foot, PenethInfo1.Foot);
	}
    }

    else{	
	Left = Node -> Left;
	Right = Node -> Right;
	/* Transform bounding volumes. */ 
	Mvar2CtTransformBV(Left, &TLBV1, &TSBV1, Xtrans, Ytrans, Rot);
	Mvar2CtTransformBV(Right, &TLBV2, &TSBV2, Xtrans, Ytrans, Rot);
	/* Compute penetration depth of center of bounding sphere. */
	Mvar2CtPtCrvMindist(TSBV1.Center, BvhBs, BSize, &PenethInfo1);
	Mvar2CtPtCrvMindist(TSBV2.Center, BvhBs, BSize, &PenethInfo2);

	if (PenethInfo1.Mindist > 0) {
	    if (PenethInfo1.Mindist > TSBV1.Radius)
	     /* If penetration is positive and bigger than radius of */ 
             /* bounding sphere, this node does not have inter-penetration. */
		Result1 = FALSE;
	}

	if (Info -> Mindist > PenethInfo1.Mindist) {
	    /* If we find deeper penetration, update it. */
	    *Info = PenethInfo1;
	    IRIT_PT_COPY(Info -> P, TSBV1.Center);
	}

	Dist1 = IRIT_PT2D_DIST(PenethInfo1.Foot, TLBV1.P[0]) + TLBV1.Epsilon;
	Dist2 = IRIT_PT2D_DIST(PenethInfo1.Foot, TLBV1.P[1]) + TLBV1.Epsilon;

	if (IRIT_MAX(Dist1, Dist2) < IRIT_FABS(Info -> Mindist))
	    /* In this case this node cannot have deeper penetration */
	    /* than current penetration depth. */
	    Result1 = FALSE;

	if (PenethInfo2.Mindist > 0) {
	    if (PenethInfo2.Mindist > TSBV2.Radius)
		/* If penetration is positive and bigger than radius of */ 
		/* bounding sphere, this node does not have inter-penetration. */
		Result2 = FALSE;
	}

	if (Info -> Mindist > PenethInfo2.Mindist) {
	     /* If we find deeper penetration, update it. */
	    *Info = PenethInfo2;
	    IRIT_PT_COPY(Info -> P, TSBV2.Center);
	}

	Dist1 = IRIT_PT2D_DIST(PenethInfo2.Foot, TLBV2.P[0]) + TLBV2.Epsilon;
	Dist2 = IRIT_PT2D_DIST(PenethInfo2.Foot, TLBV2.P[1]) + TLBV2.Epsilon;

	if (IRIT_MAX(Dist1, Dist2) < IRIT_FABS(Info -> Mindist))
	    /* In this case this node cannot have deeper penetration */
	    /* than current penetration depth. */
	    Result2 = FALSE;

	/* Repeat the same procedure. */	
	if (Result1 && Result2) {
	    if (PenethInfo1.Mindist < PenethInfo2.Mindist) {
		Mvar2CtAuxPenetrationDepth(Left, BvhBs, BSize, 
		                           Xtrans, Ytrans, Rot, Info);
		Mvar2CtAuxPenetrationDepth(Right, BvhBs, BSize, 
		                           Xtrans, Ytrans, Rot, Info);
	    }
	    else {
		Mvar2CtAuxPenetrationDepth(Right, BvhBs, BSize, 
		                           Xtrans, Ytrans, Rot, Info);
		Mvar2CtAuxPenetrationDepth(Left, BvhBs, BSize, 
		                           Xtrans, Ytrans, Rot, Info);
	    }
	}
	else if (Result1 && !Result2) {
	    Mvar2CtAuxPenetrationDepth(Left, BvhBs, BSize, 
		                       Xtrans, Ytrans, Rot, Info);    
	}
	else if (!Result1 && Result2) {
	    Mvar2CtAuxPenetrationDepth(Right, BvhBs, BSize, 
		                       Xtrans, Ytrans, Rot, Info);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Compute a maximum penetration depth of moving curve into obstacle curves  M
*  using BVHs.                                                               M
*                                                                            *
* PARAMETERS:                                                                M
*   BvhA:   Bvh for moving curve.                                     	     M
*   BvhBs:  Bvh for obstacle curves.                                         M
*   BSize:  Number of obstacle curves.                                       M
*   Xtrans: X translation.                                                   M
*   Ytrans: Y translation.                                                   M
*   Rot: Rotation angle in radian.                                           M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType: Maximum penetration depth, negative if moving curve penetrate M
*              obstacle curves positive, otherwise.                          M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarStewartPlatformSolve                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtPenetrationDepth                                                  M
*****************************************************************************/
CagdRType Mvar2CtPenetrationDepth(Mvar2CtBVHStruct *BvhA, 
				  Mvar2CtBVHStruct **BvhBs,
				  int BSize,
				  CagdRType Xtrans,
				  CagdRType Ytrans,
				  CagdRType Rot)
{
    Mvar2CtMindistStruct Info;

    Info.Mindist = 0;

    Mvar2CtAuxPenetrationDepth(BvhA -> Root, BvhBs, BSize, 
			       Xtrans, Ytrans, Rot, &Info);
   
    return Info.Mindist;
}

