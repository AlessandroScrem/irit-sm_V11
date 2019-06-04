/******************************************************************************
* mv2ctaux.c - 2contact motion analysis for planar curves.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yong-Joon Kim, July 2014.					      *
******************************************************************************/

#include "inc_irit/misc_lib.h"
#include "mvar_loc.h"

#define MVAR2CT_PERIODIC_TOL 1e-10
#define MVAR2CT_CONNECT_ANGLE_PERIODIC

static int Mvar2CtSingleSubdivDir(Mvar2CtBVNodeStruct *ANode, 
				  Mvar2CtBVNodeStruct *BNode,
				  Mvar2CtCParamStruct *Cparam);
static int Mvar2CtDoubleSubdivDir(Mvar2CtBVNodeStruct *Nodes[4],
				  Mvar2CtCParamStruct *Cparam);
static int Mvar2CtTripleSubdivDir(Mvar2CtBVNodeStruct *Nodes[6],
				  Mvar2CtCParamStruct *Cparam);
static CagdBType Mvar2CtNormalOverlapDir(Mvar2CtBVNodeStruct *Nodes[4], 
					 CagdRType RMin, 
					 CagdRType RMax,
					 int MinMax,
					 int FixedDir);
static CagdBType Mvar2CtNodeOverlap(Mvar2CtBVNodeStruct *Nodes[4], 
				    CagdRType RMin,
				    CagdRType RMax);
static CagdBType Mvar2CtNodeOverlapDir(Mvar2CtBVNodeStruct *Nodes[4], 
				       CagdRType RMin, 
				       CagdRType RMax,
				       int MinMax,
				       int FixedDir);
static void Mvar2CtGetRequiredRot(Mvar2CtConeStruct *A, 
				  Mvar2CtConeStruct *B, 
				  CagdRType *RMin, 
				  CagdRType *RMax);
static CagdRType Mvar2CtMaxDomainLength(Mvar2CtBVNodeStruct *Nodes[4]);
static CagdRType Mvar2CtPeriodicDiff(CagdRType a, CagdRType b);
static CagdRType Mvar2CtPeriodicDist(MvarPtStruct *Pt1, MvarPtStruct *Pt2);
static CagdBType Mvar2CtIsSameBndry(MvarPtStruct *Pt1, MvarPtStruct *Pt2);
static int Mvar2CtIsConnectedTrace(MvarPolylineStruct *Poly1, 
				   MvarPolylineStruct *Poly2, 
				   CagdRType Numerictol);

static CagdBType Mvar2CtTraceCollideAux(MvarPtStruct *SPt1, 
					MvarPtStruct *EPt1, 
					MvarPtStruct *SPt2, 
					MvarPtStruct *EPt2, 
					int Type);

static void Mvar2CtMakeLSCRot(Mvar2CtBVNodeStruct *Node,
			      CagdRType RMin,
			      CagdRType RMax,
			      Mvar2CtLineStruct *LBV);
static void Mvar2CtMakeLSCRotPt(CagdPType P,
				CagdRType RMin,
				CagdRType RMax,
				Mvar2CtLineStruct *LBV);
static void Mvar2CtMakeLSCPair(Mvar2CtLineStruct *ALBV,
			       Mvar2CtLineStruct *BLBV,
			       Mvar2CtLineStruct *LBV);
static void Mvar2CtMakeLSCPairPt(Mvar2CtLineStruct *ALBV,
				 CagdPType P,
				 Mvar2CtLineStruct *LBV);


/*****************************************************************************
* DESCRIPTION:                                                               M
*   Check if two curves can have a same curvature by using BVH.              M
*                                                                            *
* PARAMETERS:                                                                M
*   ANode:   Bounding volume for a moving curve.                             M
*   BNode:   Bounding volume for a obstacle curve.                           M
*   Cparam:  Data structure for rotation.                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: FALSE if there is no curvature overlap TRUE otherwise.        M 
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtCurvatureOverlap                                                  M
*****************************************************************************/
CagdBType Mvar2CtCurvatureOverlap(Mvar2CtBVNodeStruct *ANode, 
				  Mvar2CtBVNodeStruct *BNode,
				  Mvar2CtCParamStruct *Cparam)
{
    int Dir;

    if (ANode -> KMin > BNode -> KMax)
	return FALSE;
    if (ANode -> KMax < BNode -> KMin)
	return FALSE;

    if (!Mvar2CtNormalOverlap(ANode, BNode, Cparam -> RMin, Cparam -> RMax))
	return FALSE;

    Dir = Mvar2CtSingleSubdivDir(ANode, BNode, Cparam);

    switch (Dir) {
        case 0:
	    return Mvar2CtCurvatureOverlap(ANode -> Left, BNode, Cparam) ||
	           Mvar2CtCurvatureOverlap(ANode -> Right, BNode, Cparam);
	    break; 
        case 1:
	    return Mvar2CtCurvatureOverlap(ANode, BNode -> Left, Cparam) ||
	           Mvar2CtCurvatureOverlap(ANode, BNode -> Right, Cparam);
	    break; 
        case 2:
	    return Mvar2CtCurvatureOverlap(ANode, BNode, Cparam -> Left) || 
	           Mvar2CtCurvatureOverlap(ANode, BNode, Cparam -> Right);
	    break; 
        case 3:
	    return TRUE;
    }

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Build a hierarchy for rotation data structure.                            M
*                                                                            *
* PARAMETERS:                                                                M
*   Circle: Unit circle.                                                     M
*   Node: New node will be allocated herein.                                 M
*   Min, Max: Domain for the node.                                           M
*   Tol: Tolerance for building hierarchy.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtBuildCParamHierarchy                                              M
*****************************************************************************/
void Mvar2CtBuildCParamHierarchy(CagdCrvStruct *Circle, 
  				 Mvar2CtCParamStruct **Node, 
				 CagdRType Min, 
				 CagdRType Max, 
				 CagdRType Tol)
{
    CagdRType *R;
    CagdPType P;

    (*Node) =  
	(Mvar2CtCParamStruct *) IritMalloc(sizeof(Mvar2CtCParamStruct));

    (*Node) -> TMin = Min;
    (*Node) -> TMax = Max;

    (*Node) -> Left = (*Node) -> Right = NULL;  

    if (Min == 0)
	(*Node) -> RMin = 0;
    else {
	R = CagdCrvEval(Circle, Min);
	CagdCoerceToE2(P, &R, -1, Circle -> PType);
	/* Compute rotation angle correspond to Min. */
	(*Node) -> RMin = Mvar2CtGetTheta(P[0], P[1]);
    }

    if (Max == 1)
	(*Node) -> RMax = M_PI_MUL_2;

    else { 
	R = CagdCrvEval(Circle, Max);
	CagdCoerceToE2(P, &R, -1, Circle -> PType);
	/* Compute rotation angle correspond to Max. */
	(*Node) -> RMax = Mvar2CtGetTheta(P[0], P[1]);
    }
   
    if (Max - Min > Tol) {
        Mvar2CtBuildCParamHierarchy(Circle, &((*Node) -> Left), 
	                            Min, (Min + Max) * 0.5, Tol);
	Mvar2CtBuildCParamHierarchy(Circle, &((*Node) -> Right), 
	                            (Min + Max) * 0.5, Max, Tol);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Find smallest bounding volume includes the domain Min-Max.                M
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
*   Mvar2CtGetParentCparam                                                   M
*****************************************************************************/
void Mvar2CtGetParentCparam(Mvar2CtCParamStruct *Node,
			    CagdRType Min,
			    CagdRType Max,
			    Mvar2CtCParamStruct **Parent)
{
    if (Node -> TMin <= Min && Node -> TMax >= Max) {		
        (*Parent) = Node;
	if (Node -> Left != NULL) {
	    Mvar2CtGetParentCparam(Node -> Left, Min, Max, Parent);
	    Mvar2CtGetParentCparam(Node -> Right, Min, Max, Parent);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Free Mvar2CtCParamStruct recursively.                                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Node: Mvar2CtCParamStruct node to free.                                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtFreeCparam                                                        M
*****************************************************************************/
void Mvar2CtFreeCparam(Mvar2CtCParamStruct *Node)
{
    if (Node -> Left != NULL) {
	Mvar2CtFreeCparam(Node -> Left);
	Mvar2CtFreeCparam(Node -> Right);
    }
    IritFree(Node);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Find a subdivision direction for pair of curves.                          *
*                                                                            *
* PARAMETERS:                                                                *
*   ANode: Bounding volume for moving curve.                                 *
*   BNode: Bounding volume for obstacle curve.                               *
*   Cparam: Data structure for rotation.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: direction for next subdivision.                                     *
*****************************************************************************/
static int Mvar2CtSingleSubdivDir(Mvar2CtBVNodeStruct *ANode,
			          Mvar2CtBVNodeStruct *BNode,
			          Mvar2CtCParamStruct *Cparam)
{
    CagdRType 
	Max = 0;
    int Dir = 3;

    if (ANode -> Max - ANode -> Min > Max && ANode -> Left != NULL) {
	Max = ANode -> Max - ANode -> Min;
	Dir = 0;
    }

    if (BNode -> Max - BNode -> Min > Max && BNode -> Left != NULL) {
	Max = BNode -> Max - BNode -> Min;
	Dir = 1;
    }

    if (Cparam -> TMax - Cparam -> TMin > Max && Cparam -> Left != NULL){
	Max = Cparam -> TMax - Cparam -> TMin;
	Dir = 2;
    }

    return Dir;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Find a subdivision direction for 2 contact computation.                   *
*                                                                            *
* PARAMETERS:                                                                *
*   Nodes: Bounding volumes for curves.                                      *
*   Cparam: Data structure for rotation.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: direction for next subdivision.                                     *
*****************************************************************************/
static int Mvar2CtDoubleSubdivDir(Mvar2CtBVNodeStruct *Nodes[4], 
				  Mvar2CtCParamStruct *Cparam)
{
    CagdRType 
	Max = 0;
    int i,
	Dir=5;

    for (i = 0; i < 4; ++i) {
	if (Nodes[i] -> Max - Nodes[i] -> Min > Max && 
	   Nodes[i] -> Left != NULL) {
	    Max = Nodes[i] -> Max - Nodes[i] -> Min;
	    Dir = i;
	}
    }

    if (Cparam -> TMax - Cparam -> TMin > Max && 
	Cparam -> Left != NULL) {
	Max = Cparam -> TMax - Cparam -> TMin;
	Dir = 4;
    }

    return Dir;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Find a subdivision direction for 3 contact computation.                   *
*                                                                            *
* PARAMETERS:                                                                *
*   Nodes: Bounding volumes for curves.                                      *
*   Cparam: Data structure for rotation.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: Direction for next subdivision.                                     *
*****************************************************************************/
static int Mvar2CtTripleSubdivDir(Mvar2CtBVNodeStruct *Nodes[6], 
				  Mvar2CtCParamStruct *Cparam)
{
    CagdRType 
	Max = 0;
    int i,
	Dir=7;

    for (i = 0; i < 6; ++i) {
	if (Nodes[i] -> Max - Nodes[i] -> Min > Max && 
	    Nodes[i] -> Left != NULL) {
	    Max = Nodes[i] -> Max - Nodes[i] -> Min;
	    Dir = i;
	}
    }

    if (Cparam -> TMax - Cparam -> TMin > Max && 
	Cparam -> Left != NULL) {
	Max = Cparam -> TMax - Cparam -> TMin;
	Dir=6;
    }

    return Dir;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Compute a rotation angle in radian correspond (x, y) vector.              M
*                                                                            *
* PARAMETERS:                                                                M
*   x, y: vector for computing rotation angle.                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType: Rotation angle in radian.                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtGetTheta                                                          M
*****************************************************************************/
CagdRType Mvar2CtGetTheta(CagdRType x, CagdRType y)
{
    CagdRType 
	Length = sqrt(x * x + y * y);

    x = x / Length;
    y = y / Length;

    x = IRIT_BOUND(x, -1, 1);
    y = IRIT_BOUND(y, -1, 1);

    if (y >= 0)
	return acos(x);
    
    else
	return M_PI_MUL_2 - acos(x);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Check the tangency condition for 2 contact points using BVH.              M
*                                                                            M
*                                                                            *
* PARAMETERS:                                                                M
*   ANode, BNode:  Bounding volumes to test.                                 M
*   RMin:      Minimum rotation angle in radian.		             M
*   RMax:      Maximum rotation angle in radian.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: FALSE if tangency condition does not satisfy TRUE             M
*              otherwise.                                                    M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtNormalOverlap                                                     M
*****************************************************************************/
CagdBType Mvar2CtNormalOverlap(Mvar2CtBVNodeStruct *ANode, 
			       Mvar2CtBVNodeStruct *BNode, 
			       CagdRType RMin, 
		               CagdRType RMax)
{
    CagdPType Axis;
    CagdRType Dot, Angle,
        Rot = (RMin + RMax) * 0.5,  
        AngleError = (RMax - RMin) * 0.5;

    Axis[0] = - (ANode -> NCone.Axis[0] * cos(Rot) - 
	                                ANode -> NCone.Axis[1] * sin(Rot));
    Axis[1] = - (ANode -> NCone.Axis[1] * cos(Rot) + 
	                                ANode -> NCone.Axis[0] * sin(Rot));

    Dot = IRIT_DOT_PROD_2D(Axis, BNode -> NCone.Axis);
    Dot = IRIT_BOUND(Dot, -1, 1);

    Angle = acos(Dot);

    if (Angle > ANode -> NCone.Angle + BNode -> NCone.Angle + AngleError)
	return FALSE;
    else 
	return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Check the tangency condition for 2 contact points at domain boundary      * 
*  using BVH.                                                                *
*                                                                            *
* PARAMETERS:                                                                *
*   Nodes:     Bounding volumes for curves.                                  *
*   RMin:      Minimum rotation angle in radian.		             *
*   RMax:      Maximum rotation angle in radian.                             *
*   MinMax:      0 if fixed parameter is minimum parameter of the bounding   *
                 volume, 1 if if fixed parameter is maximum parameter of the *
		 of the bounding volume.                                     *
*   FixedDir:   Fixed dimension.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: FALSE if tangency condition does not satisfy TRUE             *
*              otherwise .                                                   *
*****************************************************************************/
static CagdBType Mvar2CtNormalOverlapDir(Mvar2CtBVNodeStruct *Nodes[4], 
				         CagdRType RMin, 
				         CagdRType RMax,
				         int MinMax,
				         int FixedDir)
{
    CagdPType Axis1, Axis2;
    CagdRType Rot, Dot, Angle, AngleA, AngleB, AngleError, RotMin1,
        RotMax1, RotMin2, RotMax2;
    int i;

    if (FixedDir == 4) {
	/* If rotation fixed, RMin=RMax. */
	if (!MinMax)
	    RMax = RMin;
	else
	    RMin = RMax;    
    }    
    
    Rot = (RMin + RMax) * 0.5;
    /* Get required rotation for normal overlap. */
    Mvar2CtGetRequiredRot(&(Nodes[0] -> NCone), &(Nodes[1] -> NCone), 
	                  &RotMin1, &RotMax1);

    if (RotMin1 > RotMax1) {
	/* In this case, rotation includes (0 or 2Pi). */ 
	if (Rot < M_PI_MUL_2 - Rot)
	    /* Resolve it depending on the rot value. */ 
	    RotMin1 = RotMin1 - M_PI_MUL_2;
	else	   
	    RotMax1 = M_PI_MUL_2 + RotMax1;	
    }
     /* Get required rotation for normal overlap. */
    Mvar2CtGetRequiredRot(&(Nodes[2] -> NCone), &(Nodes[3] -> NCone), 
	                  &RotMin2, &RotMax2);

    if (RotMin2 > RotMax2) {
	/* In this case, rotation includes (0 or 2Pi). */ 
	if (Rot < M_PI_MUL_2 - Rot)	 
	    /* Resolve it using the rot value. */ 
	    RotMin2 = RotMin2 - M_PI_MUL_2;
	else
	    RotMax2 = M_PI_MUL_2 + RotMax2;  	
    }

    RMin = IRIT_MAX3(RotMin1, RotMin2, RMin);
    RMax = IRIT_MIN3(RotMax1, RotMax2, RMax);

    if (RMin > RMax)
	return FALSE;

    AngleError = (RMax - RMin) * 0.5;    
    
    Rot = (RMin + RMax) * 0.5;

    /* Check normal overlap. */
    for (i = 0; i < 2; ++i) {
	if (FixedDir == i * 2) {
	    if (!MinMax) {
		Axis1[0] = - (Nodes[i * 2] -> NMin[0] * cos(Rot) - 
		              Nodes[i * 2] -> NMin[1] * sin(Rot));
		Axis1[1] = - (Nodes[i * 2] -> NMin[1] * cos(Rot) + 
		              Nodes[i * 2] -> NMin[0] * sin(Rot));
	    }
	    else {
		Axis1[0] = - (Nodes[i * 2] -> NMax[0] * cos(Rot) - 
		              Nodes[i * 2] -> NMax[1] * sin(Rot));
		Axis1[1] = - (Nodes[i * 2] -> NMax[1] * cos(Rot) + 
		              Nodes[i * 2] -> NMax[0] * sin(Rot));
	    }

	    AngleA = 0;
	}
	else {    
	    Axis1[0] = - (Nodes[i * 2] -> NCone.Axis[0] * cos(Rot) -
		          Nodes[i * 2] -> NCone.Axis[1] * sin(Rot));
	    Axis1[1] = - (Nodes[i * 2] -> NCone.Axis[1] * cos(Rot) +
		          Nodes[i * 2] -> NCone.Axis[0] * sin(Rot));

	    AngleA = Nodes[i * 2] -> NCone.Angle;
	}

	if (FixedDir == i * 2 + 1) {
	    if (!MinMax)
		IRIT_PT2D_COPY(Axis2, Nodes[i * 2 + 1] -> NMin);	    
	    else
		IRIT_PT2D_COPY(Axis2, Nodes[i * 2 + 1] -> NMax);
	    AngleB = 0;
	}
	else {    
	    IRIT_PT2D_COPY(Axis2, Nodes[i * 2 + 1] -> NCone.Axis);
	    AngleB = Nodes[i * 2 + 1] -> NCone.Angle;
	}

	Dot = IRIT_BOUND(IRIT_DOT_PROD_2D(Axis1, Axis2), -1, 1);
	Angle = acos(Dot);

	if (Angle > AngleA + AngleB + AngleError)
	    return FALSE;
    }
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Check normal overlap for both direction.                                  M
*                                                                            *
* PARAMETERS:                                                                M
*   ANode, BNode:  Nodes for normal overlapping test.                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: TRUE if normal cones overlap FALSE otherwise.                 M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtNormalOverlapBoth                                                 M
*****************************************************************************/
CagdBType Mvar2CtNormalOverlapBoth(Mvar2CtBVNodeStruct *ANode, 
				   Mvar2CtBVNodeStruct *BNode)
{
    CagdRType Dot, Angle;

    Dot = IRIT_DOT_PROD_2D(ANode -> NCone.Axis, BNode -> NCone.Axis);
    Dot = IRIT_BOUND(Dot, -1, 1);
    Angle = acos(Dot);

    Angle = IRIT_MIN(Angle, M_PI - Angle); /* Check for both direction. */

    if (Angle > ANode -> NCone.Angle + BNode -> NCone.Angle)
	return FALSE;

    else 
	return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Check the intersection condition for 2 contact trace using BVH.           *
*                                                                            *
* PARAMETERS:                                                                *
*   Nodes:    Bounding volumes for curves.                                   *
*   RMin:      Minimum rotation angle in radian.		             *
*   RMax:      Maximum rotation angle in radian.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: FALSE if intersection condition does not satisfy TRUE         *
*              otherwise.                                                    *
*****************************************************************************/
static CagdBType Mvar2CtNodeOverlap(Mvar2CtBVNodeStruct *Nodes[4], 
			            CagdRType RMin,
			            CagdRType RMax)
{
    IRIT_STATIC_DATA Mvar2CtLineStruct ALBV1, ALBV2, Minkowski1, Minkowski2;
    CagdRType Dist;

     /* Generate LSCs(ALBV1, ALBV2) which bound moving curves. */  
    Mvar2CtMakeLSCRot(Nodes[0], RMin, RMax, &ALBV1);
    Mvar2CtMakeLSCRot(Nodes[2], RMin, RMax, &ALBV2);

    /* Generate LSCs(Minkowski1, Minkowski) which bounds Minkowski sum of  */ 
    /* moving curves and obstacle curves.                                  */
    Mvar2CtMakeLSCPair(&ALBV1, &(Nodes[1] -> LBV), &Minkowski1);
    Mvar2CtMakeLSCPair(&ALBV2, &(Nodes[3] -> LBV), &Minkowski2);
    
    Dist = Mvar2CtLineLineDist(&Minkowski1, &Minkowski2);

    if (Dist > Minkowski1.Epsilon + Minkowski2.Epsilon)
	return FALSE;
    else
	return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Check the intersection condition for 2 contact points at domain boundary  *
*  using  BVH.                                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   Nodes:    Bounding volumes for curves.                                   *
*   RMin:      Minimum rotation angle in radian.		             *
*   RMax:      Maximum rotation angle in radian.                             *
*   MinMax:      0 if fixed parameter is minimum parameter of the bounding   *
                 volume, 1 if if fixed parameter is maximum parameter of the *
		 of the bounding volume.                                     *
*   FixedDir:   Fixed dimension.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: FALSE if intersection condition does not satisfy TRUE         *
*              otherwise.                                                    *
*****************************************************************************/
static CagdBType Mvar2CtNodeOverlapDir(Mvar2CtBVNodeStruct *Nodes[4], 
				       CagdRType RMin, 
				       CagdRType RMax,
				       int MinMax,
				       int FixedDir)
{
    IRIT_STATIC_DATA Mvar2CtLineStruct ALBV1, ALBV2, Minkowski1, Minkowski2;
    CagdRType Dist, RotMin1, RotMin2, RotMax1, RotMax2;

    if (FixedDir == 4) {
	/* If rotation is fixed RMin=RMax. */
	if (!MinMax)
	    RMax = RMin;
	else
	    RMin = RMax;
    }
    else {
	/* Get required rotation for tangency condition. */
	Mvar2CtGetRequiredRot(&(Nodes[0] -> NCone), &(Nodes[1] -> NCone), 
	                      &RotMin1, &RotMax1);
	Mvar2CtGetRequiredRot(&(Nodes[2] -> NCone), &(Nodes[3] -> NCone), 
	                      &RotMin2, &RotMax2);

	RMin = IRIT_MAX3(RotMin1, RotMin2, RMin);
	RMax = IRIT_MIN3(RotMax1, RotMax2, RMax);
    }

    /* Generate LSCs(ALBV1, ALBV2) which bound moving curves. */  
    if (FixedDir == 0) {
	if (!MinMax)
	    Mvar2CtMakeLSCRotPt(Nodes[0] -> LBV.P[0], RMin, RMax, &ALBV1);
	else
	    Mvar2CtMakeLSCRotPt(Nodes[0] -> LBV.P[1], RMin, RMax, &ALBV1);
    }
    else
	 Mvar2CtMakeLSCRot(Nodes[0], RMin, RMax, &ALBV1);
   
    if (FixedDir == 2) {
	if (!MinMax)
	    Mvar2CtMakeLSCRotPt(Nodes[2] -> LBV.P[0], RMin, RMax, &ALBV2);
	else
	    Mvar2CtMakeLSCRotPt(Nodes[2] -> LBV.P[1], RMin, RMax, &ALBV2);
    }
    else
	Mvar2CtMakeLSCRot(Nodes[2], RMin, RMax, &ALBV2);
   
    /* Generate LSCs(Minkowski1, Minkowski) which bounds Minkowski sum of  */ 
    /* moving curves and obstacle curves.                                  */
    if (FixedDir == 1) {
	if (!MinMax)
	     Mvar2CtMakeLSCPairPt(&ALBV1, Nodes[1] -> LBV.P[0], &Minkowski1);

	else
	     Mvar2CtMakeLSCPairPt(&ALBV1, Nodes[1] -> LBV.P[1], &Minkowski1);
    }
    else
	Mvar2CtMakeLSCPair(&ALBV1, &(Nodes[1] -> LBV), &Minkowski1);

    if (FixedDir == 3) {
	if (!MinMax)
	    Mvar2CtMakeLSCPairPt(&ALBV2, Nodes[3] -> LBV.P[0], &Minkowski2);
	else
	    Mvar2CtMakeLSCPairPt(&ALBV2, Nodes[3] -> LBV.P[1], &Minkowski2);
    }
    else
	Mvar2CtMakeLSCPair(&ALBV2, &(Nodes[3] -> LBV), &Minkowski2);

    Dist = Mvar2CtLineLineDist(&Minkowski1, &Minkowski2);

    if (Dist > Minkowski1.Epsilon + Minkowski2.Epsilon)
	return FALSE;
    else
	return TRUE;
}


/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate LSC that bounds the sweeping region of P undergoes rotation      *
*  from RMin to RMax.                                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   P: Point to consider.                                                    *
*   RMin : Minimum rotation angle.                                           * 
*   RMax : Maximum rotation angle.                                           *
*   LBV : the sweeping region of P undergoes rotation                        *
*         from RMin to RMax.                                                 * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeLSCRotPt(CagdPType P,
           		        CagdRType RMin,
		                CagdRType RMax,
		                Mvar2CtLineStruct *LBV)
{
    CagdRType Length;

    LBV -> P[0][0] = - (P[0] * cos(RMin) - P[1] * sin(RMin));
    LBV -> P[0][1] = - (P[0] * sin(RMin) + P[1] * cos(RMin));

    LBV -> P[1][0] = - (P[0] * cos(RMax) - P[1] * sin(RMax));
    LBV -> P[1][1] = - (P[0] * sin(RMax) + P[1] * cos(RMax));

    Length = IRIT_PT2D_DIST(LBV -> P[0], LBV -> P[1]);

    LBV -> Epsilon = Length * (1 - cos((RMax - RMin) * 0.5));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate LSC that bounds the sweeping region of Node undergoes rotation   *
*  from RMin to RMax.                                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   Node: Bounding volume of the curve.                                      *
*   RMin : Minimum rotation angle.                                           * 
*   RMax : Maximum rotation angle.                                           *
*   LBV : LSC that bounds the sweeping region of Node undergoes rotation     *
*         from RMin to RMax will be placed herein.                           * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeLSCRot(Mvar2CtBVNodeStruct *Node,
            		      CagdRType RMin,
		              CagdRType RMax,
		              Mvar2CtLineStruct *LBV)
{
    IRIT_STATIC_DATA CagdPType P[4];
    CagdRType Dist1, Dist2;
    int i;

    for (i = 0; i < 2; ++i) {
	/* Compute Minkowski sum of two line segments. */
	P[i][0] = - (Node -> LBV.P[i][0] * cos(RMin) - 
	             Node -> LBV.P[i][1] * sin(RMin));
	P[i][1] = - (Node -> LBV.P[i][0] * sin(RMin) + 
	             Node -> LBV.P[i][1] * cos(RMin));

	P[i + 2][0] = - (Node -> LBV.P[i][0] * cos(RMax) - 
	                 Node -> LBV.P[i][1] * sin(RMax));
	P[i + 2][1] = - (Node -> LBV.P[i][0] * sin(RMax) + 
	                 Node -> LBV.P[i][1] * cos(RMax));
    }

    /* Set two end points. */
    if (IRIT_PT2D_DIST(P[0], P[3]) > IRIT_PT2D_DIST(P[1], P[2])) {
	IRIT_PT2D_COPY(LBV -> P[0], P[0]);
	IRIT_PT2D_COPY(LBV -> P[1], P[3]);

	LBV -> Epsilon = Node -> LBV.Epsilon;
	Dist1 = Mvar2CtLinePointDist(P[1], LBV -> P);
	Dist2 = Mvar2CtLinePointDist(P[2], LBV -> P);

	if (Dist1 > Dist2)
	    LBV -> Epsilon += Dist1;	
	else
	    LBV -> Epsilon += Dist2;

	LBV -> Epsilon += Node -> Length * (1 - cos((RMax - RMin) * 0.5));
    }
    else {
	IRIT_PT2D_COPY(LBV -> P[0], P[1]);
	IRIT_PT2D_COPY(LBV -> P[1], P[2]);

	LBV -> Epsilon = Node -> LBV.Epsilon;

	Dist1 = Mvar2CtLinePointDist(P[0], LBV -> P);
	Dist2 = Mvar2CtLinePointDist(P[3], LBV -> P);

	if (Dist1 > Dist2)
	    LBV -> Epsilon += Dist1;	
	else
	    LBV -> Epsilon += Dist2;

	LBV -> Epsilon += Node -> Length * (1 - cos((RMax - RMin) * 0.5));
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate LSC bounding Minkowski sum of Point and LSC.                     *
*                                                                            *
* PARAMETERS:                                                                *
*   ALBV: Line swept circle.                                                 *
*   P : Point.                                                               * 
*   LBV : LSC bounding Minkowski sum of ALBV and P will be placed herein.    * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeLSCPairPt(Mvar2CtLineStruct *ALBV,
		                 CagdPType P,
			         Mvar2CtLineStruct *LBV)
{
    int i;

    for (i = 0; i < 2; ++i)	
	IRIT_PT2D_ADD(LBV -> P[i], ALBV -> P[i], P);		

    LBV -> Epsilon = ALBV -> Epsilon;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Generate LSC bounding Minkowski sum of two LSCs.                          *
*                                                                            *
* PARAMETERS:                                                                *
*   ALBV, BLBV: Two LSC for Minkowski sum computation.                       *
*   LBV : LSC bounding Minkowski sum of ALBV and BLBV will be placed herein. * 
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtMakeLSCPair(Mvar2CtLineStruct *ALBV,
			       Mvar2CtLineStruct *BLBV,
			       Mvar2CtLineStruct *LBV)
{
    IRIT_STATIC_DATA CagdPType P[4];
    CagdRType Dist1, Dist2;
    int i, j;

    for (i = 0; i < 2; ++i) {
	for (j = 0; j < 2; ++j)
	    IRIT_PT2D_ADD(P[i * 2 + j], ALBV -> P[i], BLBV -> P[j]);
    }

    if (IRIT_PT2D_DIST(P[0], P[3]) > IRIT_PT2D_DIST(P[1], P[2])) {
	IRIT_PT2D_COPY(LBV -> P[0], P[0]);
	IRIT_PT2D_COPY(LBV -> P[1], P[3]);

	LBV -> Epsilon = ALBV -> Epsilon + BLBV -> Epsilon;
	Dist1 = Mvar2CtLinePointDist(P[1], LBV -> P);
	Dist2 = Mvar2CtLinePointDist(P[2], LBV -> P);

	if (Dist1 > Dist2)
	    LBV -> Epsilon += Dist1;	
	else
	    LBV -> Epsilon += Dist2;
    }
    else {
	IRIT_PT2D_COPY(LBV -> P[0], P[1]);
	IRIT_PT2D_COPY(LBV -> P[1], P[2]);

	LBV -> Epsilon = ALBV -> Epsilon + BLBV -> Epsilon;

	Dist1 = Mvar2CtLinePointDist(P[0], LBV -> P);

	Dist2 = Mvar2CtLinePointDist(P[3], LBV -> P);

	if (Dist1 > Dist2)
	    LBV -> Epsilon += Dist1;	
	else
	    LBV -> Epsilon += Dist2;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Compute required rotation range for normal overlapping.                   *
*                                                                            *
* PARAMETERS:                                                                *
*   A, B:   Normal cones to test.                                            * 
*   RMin:   Minimum required rotation angle in radian will be placed herein. * 
*   RMax:   Maximum required rotation angle in radian will be placed herein. *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtGetRequiredRot(Mvar2CtConeStruct *A, 
			          Mvar2CtConeStruct *B, 
			          CagdRType *RMin, 
			          CagdRType *RMax)
{
    CagdPType Axis;
    CagdRType Dot, Cross, Angle, AngleError;

    IRIT_PT2D_SCALE2(Axis, A -> Axis, -1.0);
    Dot = IRIT_DOT_PROD_2D(Axis, B -> Axis);
    Cross = IRIT_CROSS_PROD_2D(Axis, B -> Axis);

    Angle = acos(IRIT_BOUND(Dot, -1, 1));

    AngleError = A -> Angle + B -> Angle;

    if (Cross < 0) {
        *RMin = M_PI_MUL_2 - Angle - AngleError;
	*RMax = M_PI_MUL_2 - Angle + AngleError;
    }
    else {
	*RMin = Angle - AngleError;
	*RMax = Angle + AngleError;	
    }

    if (*RMin < 0)
	*RMin = M_PI_MUL_2 + (*RMin);
    if (*RMax > M_PI_MUL_2)
	*RMax = (*RMax) - M_PI_MUL_2;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Check if pair of curves satisfy curvature condition for 2 contact.        M
*                                                                            *
* PARAMETERS:                                                                M
*   Node1, Node2: Bounding volumes of curves to test.                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: TRUE if they don't satisfy the condition, FALSE otherwise.    M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtRejectbyCurvature                                                 M
*****************************************************************************/
CagdBType Mvar2CtRejectbyCurvature(Mvar2CtBVNodeStruct *Node1, 
				   Mvar2CtBVNodeStruct *Node2)
{
    if (Node1 -> Convexity < 0 && Node2 -> Convexity < 0)
	/* If both curves are concave we can purge it. */
        return TRUE;

    if (Node1 -> Convexity < 0 && Node2 -> Convexity > 0) {
	if (Node1 -> KMin > Node2 -> KMax)
	    /* Concave curve should have larger curvature. */
	    return TRUE;
    }

    else if (Node1 -> Convexity > 0 && Node2 -> Convexity <0 ) {
	if (Node1 -> KMax < Node2 -> KMin)
	    /* Concave curve should have larger curvature. */
	    return TRUE;
    }

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Check if two bounding volumes are connected.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   Node1, Node2:    Bounding volumes to test.                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: TRUE if connected, FALSE otherwise.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtIsConnectedNode                                                   M
*****************************************************************************/
CagdBType Mvar2CtIsConnectedNode(Mvar2CtBVNodeStruct *Node1, 
			         Mvar2CtBVNodeStruct *Node2)
{
    if (Node1 -> Max < Node2 -> Min && 
	!(Node1 -> Min == 0 && Node2 -> Max == 1) && 
        Node1 -> Id == Node2 -> Id)
	return FALSE;
    
    if (Node1 -> Min > Node2 -> Max && 
	!(Node1 -> Max == 1 && Node2 -> Min == 0) && 
        Node1 -> Id == Node2 -> Id)
	return FALSE;

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find a maximum domain length of 4 curves.                                *
*                                                                            *
* PARAMETERS:                                                                *
*   Nodes:    Bounding volumes for curves.                                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: Maximum domain length of 4 curve segments.                    *
*****************************************************************************/
static CagdRType Mvar2CtMaxDomainLength(Mvar2CtBVNodeStruct *Nodes[4])
{
    CagdRType Max=0;
    int i;

    for (i = 0; i < 4; ++i) {
	if (Nodes[i] -> Max - Nodes[i] -> Min > Max)
	    Max = Nodes[i] -> Max - Nodes[i] -> Min;
    }

    return Max;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Check if there exists 2 contact trace in the domain using BVH.           M
*                                                                            *
* PARAMETERS:                                                                M
*   Nodes:       Bounding volumes for curves.                                M
*   Cparam:      data structure for rotation.		                     M
*   Tol:         Subdivision tolerance of checking.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: FALSE if there is no solution in the domain, TRUE otherwise.  M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtCheck2CtTrace                                                     M
*****************************************************************************/
CagdBType Mvar2CtCheck2CtTrace(Mvar2CtBVNodeStruct *Nodes[4], 
			       Mvar2CtCParamStruct *Cparam,  
			       CagdRType Tol)
{    
    Mvar2CtBVNodeStruct *Left[4], *Right[4];
    int Dir, i;
    CagdBType 
	IsConnected = FALSE;
   
    if (Nodes[0] -> Min >= Nodes[2] -> Max && 
	Nodes[1] -> Id == Nodes[3] -> Id)
	return FALSE;

    if (Mvar2CtRejectbyCurvature(Nodes[0], Nodes[1]) ||
	Mvar2CtRejectbyCurvature(Nodes[2], Nodes[3]))
	/* Check curvature condition. */
	return FALSE;

    if (Mvar2CtIsConnectedNode(Nodes[0], Nodes[2]) && 
	Mvar2CtIsConnectedNode(Nodes[1], Nodes[3])) {
	if (Mvar2CtMaxDomainLength(Nodes) < Tol)
	    /* If curve pairs are connected and the domain size is         */
	    /* smaller than tolerance we assume there is no solution.      */
	    return FALSE;
	else
	    IsConnected = TRUE;
    }

    if (!IsConnected) {
	if (!Mvar2CtNormalOverlap(Nodes[0], Nodes[1], 
	    Cparam -> RMin, Cparam -> RMax) ||
	    !Mvar2CtNormalOverlap(Nodes[2], Nodes[3], 
	    Cparam -> RMin, Cparam -> RMax))
	    /* Check tangency condition. */
	    return FALSE;
	if (!Mvar2CtNodeOverlap(Nodes, Cparam -> RMin, Cparam -> RMax))
	    /* Check intersection condition. */
	    return FALSE;
    }

    Dir = Mvar2CtDoubleSubdivDir(Nodes, Cparam);

    if (Dir < 5) {
	if (Dir < 4) {
	    for (i = 0; i < 4; ++i) {
		if (i == Dir) {
		    Left[i] = Nodes[i] -> Left;
		    Right[i] = Nodes[i] -> Right;
		}
		else {
		    Left[i] = Nodes[i];
		    Right[i] = Nodes[i];
		}
	    }

	    return Mvar2CtCheck2CtTrace(Left, Cparam, Tol) || 
	           Mvar2CtCheck2CtTrace(Right, Cparam, Tol); 
	}
	else
	    return Mvar2CtCheck2CtTrace(Nodes, Cparam -> Left, Tol) || 
		   Mvar2CtCheck2CtTrace(Nodes, Cparam -> Right, Tol);	
    }
    else
        /* We reached leaf node, thus there might exist a solution in this */
        /* domain.							   */
        return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reduce the domain for solving 2 contact points at domain boundary using  M
*  BVH.                                                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   Nodes:    Bounding volumes for curves.                                   M
*   Cparam:      Data structure for rotation.		                     M
*   Min:      Current minimum domain values.                                 M
*   Max:      Current maximum domain values.   		                     M
*   MinMax:      0 if fixed parameter is minimum parameter of the bounding   M
                 volume, 1 if if fixed parameter is maximum parameter of the M
		 of the bounding volume.                                     M
*   FixedDir:   Fixed dimension.                                             M
*   Tol:         Subdivision tolerance of checking.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtReduce2CtDomain                                                   M
*****************************************************************************/
void Mvar2CtReduce2CtDomain(Mvar2CtBVNodeStruct *Nodes[4], 
			    Mvar2CtCParamStruct *Cparam, 
			    CagdRType Min[5], 
			    CagdRType Max[5],
			    int MinMax,
			    int FixedDir,
			    CagdRType Tol)
{
    int Dir, i;
    Mvar2CtBVNodeStruct *Left[4], *Right[4];
    CagdBType
	IsConnected = FALSE;

    if (Nodes[0] -> Min >= Nodes[2] -> Max && 
	Nodes[1] -> Id == Nodes[3] -> Id)
	return;

    if (Mvar2CtRejectbyCurvature(Nodes[0], Nodes[1]) ||
	Mvar2CtRejectbyCurvature(Nodes[2], Nodes[3]))
	/* Check curvature condition. */
	return;


    if (Mvar2CtIsConnectedNode(Nodes[0], Nodes[2]) && 
	Mvar2CtIsConnectedNode(Nodes[1], Nodes[3])) {
	    if (Mvar2CtMaxDomainLength(Nodes) < Tol)
		/* If curve pairs are connected and the domain size is */
		/* smaller than tolerance we assume there is no solution. */
		return;
	    else
		IsConnected = TRUE;
    }

    if (!IsConnected) {
	if (!Mvar2CtNormalOverlapDir(Nodes, Cparam -> RMin, 
	    Cparam -> RMax, MinMax, FixedDir))
	    /* Check tangency condition. */
	    return;
	if (!Mvar2CtNodeOverlapDir(Nodes, Cparam -> RMin, 
	    Cparam -> RMax, MinMax, FixedDir))
	    /* Check intersection condition. */
	    return;
    }

    Dir = Mvar2CtDoubleSubdivDir(Nodes, Cparam);

    if (Dir < 5) {
	if (Dir < 4) {
	    for (i = 0; i < 4; ++i) {
		if (i == Dir) {
		    Left[i] = Nodes[i] -> Left;
		    Right[i] = Nodes[i] -> Right;
		}
		else {
		    Left[i] = Nodes[i];
		    Right[i] = Nodes[i];
		}
	    }

	    if (Dir == FixedDir) {
		/*We only proceed to domain containing fixed parameter.*/
		if (!MinMax)
		    Mvar2CtReduce2CtDomain(Left, Cparam, Min, Max, 
		    MinMax, FixedDir, Tol);
		else
		    Mvar2CtReduce2CtDomain(Right, Cparam, Min, Max, 
		    MinMax, FixedDir, Tol);
	    }
	    else {		
		Mvar2CtReduce2CtDomain(Left, Cparam, Min, Max, 
		    MinMax, FixedDir, Tol); 
		Mvar2CtReduce2CtDomain(Right, Cparam, Min, Max, 
		    MinMax, FixedDir, Tol);
	    }
	}
	else {
	    if (FixedDir == 4) {
		if (!MinMax)
		    Mvar2CtReduce2CtDomain(Nodes, Cparam -> Left,
		    Min, Max,
		    MinMax, FixedDir, Tol); 
		else
		    Mvar2CtReduce2CtDomain(Nodes, Cparam -> Right,
		    Min, Max,
		    MinMax, FixedDir, Tol);
	    }
	    else {		
		Mvar2CtReduce2CtDomain(Nodes, Cparam -> Left, Min, Max, 
		    MinMax, FixedDir, Tol); 
		Mvar2CtReduce2CtDomain(Nodes, Cparam -> Right, Min, Max, 
		    MinMax, FixedDir, Tol);
	    }
	}
    }

    else {
	/* We reached leaf node, update domain information. */
	for (i = 0; i < 4; ++i) {
	    if (FixedDir == i) {
		if (!MinMax) {
		    Min[i] = Nodes[i] -> Min;
		    Max[i] = Nodes[i] -> Min;
		}
		else {
		    Min[i] = Nodes[i] -> Max;
		    Max[i] = Nodes[i] -> Max;
		}
	    }
	    else {		
		if (Nodes[i] -> Min < Min[i])
		    Min[i] = Nodes[i] -> Min;

		if (Nodes[i] -> Max > Max[i])
		    Max[i] = Nodes[i] -> Max;
	    }
	}

	if (FixedDir == 4) {
	    if (!MinMax) {
		Min[4] = Cparam -> TMin;
		Max[4] = Cparam -> TMin;
	    }
	    else {
		Min[4] = Cparam -> TMax;
		Max[4] = Cparam -> TMax;
	    }

	}
	else {
	    if (Min[4] > Cparam -> TMin)
		Min[4] = Cparam -> TMin;

	    if (Max[4] < Cparam -> TMax)
		Max[4] = Cparam -> TMax;
	}
    }

}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reduce the domain for solving 3 contact points using BVH.                M
*                                                                            *
* PARAMETERS:                                                                M
*   Nodes:    Bounding volumes for curves.                                   M
*   Cparam:      Data structure for rotation.		                     M
*   Min:         Current minimum domain values.                              M
*   Max:         Current maximum domain values.   		             M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtReduce3CtDomain                                                   M
*****************************************************************************/
void Mvar2CtReduce3CtDomain(Mvar2CtBVNodeStruct *Nodes[6], 
			    Mvar2CtCParamStruct *Cparam, 
			    CagdRType Min[7], 
			    CagdRType Max[7])
{
    int Dir, i;
    Mvar2CtBVNodeStruct *Left[6], *Right[6];

    /* Check tangency condition. */ 
    if (!Mvar2CtNormalOverlap(Nodes[0], Nodes[1], 
	                      Cparam -> RMin, Cparam -> RMax))
	return;

    if (!Mvar2CtNormalOverlap(Nodes[2], Nodes[3], 
	                      Cparam -> RMin, Cparam -> RMax))
	return;

    if (!Mvar2CtNormalOverlap(Nodes[4], Nodes[5], 
	                      Cparam -> RMin, Cparam -> RMax))
	return;

    /* Check intersection condition. */
    if (Mvar2CtRejectbyCurvature(Nodes[0], Nodes[1]) ||
	Mvar2CtRejectbyCurvature(Nodes[2], Nodes[3]) ||
	Mvar2CtRejectbyCurvature(Nodes[4], Nodes[5]))
	return;

    if (!Mvar2CtNodeOverlap(Nodes, Cparam -> RMin, Cparam -> RMax) ||
	!Mvar2CtNodeOverlap(Nodes + 2, Cparam -> RMin, Cparam -> RMax))
	return;

    Dir = Mvar2CtTripleSubdivDir(Nodes, Cparam);

    if (Dir < 7) {
	if (Dir < 6) {
	    for (i = 0; i < 6; ++i) {

		if (i == Dir) {
		    Left[i] = Nodes[i] -> Left;
		    Right[i] = Nodes[i] -> Right;
		}
		else {
		    Left[i] = Nodes[i];
		    Right[i] = Nodes[i];
		}
	    }

	    Mvar2CtReduce3CtDomain(Left, Cparam, Min, Max);
	    Mvar2CtReduce3CtDomain(Right, Cparam, Min, Max);
	}
	else {
	    Mvar2CtReduce3CtDomain(Nodes, Cparam -> Left, Min, Max);
	    Mvar2CtReduce3CtDomain(Nodes, Cparam -> Right, Min, Max);
	}
    }
    else {
	/* We reached leaf nodes. */ 
	/* Update current minimum and maximum value of the domain. */
	for (i = 0; i < 6; ++i) {
	    if (Nodes[i] -> Min < Min[i])
		Min[i] = Nodes[i] -> Min;

	    if (Nodes[i] -> Max > Max[i])
		Max[i] = Nodes[i] -> Max;
	}
	
	if (Min[6] > Cparam -> TMin)
	    Min[6] = Cparam -> TMin;

	if (Max[6] < Cparam -> TMax)
	    Max[6] = Cparam -> TMax;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reduce the domain for solving rotational extreme point for2 contact      M
* trace using BVH.                                                           M
*                                                                            *
* PARAMETERS:                                                                M
*   Nodes:    Bounding volumes for curves.                                   M
*   Cparam:      data structure for rotation.		                     M
*   Min:      Current minimum domain values.                                 M
*   Max:      Current maximum domain values.   		                     M
*   Tol:         Subdivision tolerance of checking.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtReduceRotExtremeDomain                                            M
*****************************************************************************/
void Mvar2CtReduceRotExtremeDomain(Mvar2CtBVNodeStruct *Nodes[4], 
				   Mvar2CtCParamStruct *Cparam, 
				   CagdRType Min[5], 
				   CagdRType Max[5], 
				   CagdRType Tol)
{
    Mvar2CtBVNodeStruct *Left[4], *Right[4];
    int Dir, i;
    CagdBType 
	IsConnected = FALSE;

    if (Nodes[0] -> Min >= Nodes[2] -> Max && Nodes[1] -> Id == Nodes[3] -> Id)
	return;
    if (Mvar2CtRejectbyCurvature(Nodes[0], Nodes[1]) ||
	Mvar2CtRejectbyCurvature(Nodes[2], Nodes[3]))
        /* Check curvature condition. */
 	return;

    if (Mvar2CtIsConnectedNode(Nodes[0], Nodes[2]) && 
       Mvar2CtIsConnectedNode(Nodes[1], Nodes[3])) {
	if (Mvar2CtMaxDomainLength(Nodes) < Tol)
	    /* If curve pairs are connected and the domain size is         */
	    /* smaller than tolerance we assume there is no solution.      */
	    return;
	else
	    IsConnected = TRUE;
    }

    if (!IsConnected) {
	if (!Mvar2CtNormalOverlap(Nodes[0], Nodes[1], 
	                          Cparam -> RMin, Cparam -> RMax) ||
	    !Mvar2CtNormalOverlap(Nodes[2], Nodes[3], 
	                          Cparam -> RMin, Cparam -> RMax) ||
	    !Mvar2CtNormalOverlapBoth(Nodes[1], Nodes[3]))
	    /* Check tangency condition. */
	    return;
	if (!Mvar2CtNodeOverlap(Nodes, Cparam -> RMin, Cparam -> RMax))
	    /* Check intersection condition. */
	    return;
    }

    Dir = Mvar2CtDoubleSubdivDir(Nodes, Cparam);

    if (Dir < 5) {
	if (Dir < 4) {
	    for (i = 0; i < 4; ++i){
		if (i == Dir) {
		    Left[i] = Nodes[i] -> Left;
		    Right[i] = Nodes[i] -> Right;
		}
		else {
		    Left[i] = Nodes[i];
		    Right[i] = Nodes[i];
		}
	    }

	    Mvar2CtReduceRotExtremeDomain(Left, Cparam, Min, Max, Tol); 
	    Mvar2CtReduceRotExtremeDomain(Right, Cparam, Min, Max, Tol); 
	}
	else {
	    Mvar2CtReduceRotExtremeDomain(Nodes, Cparam -> Left,
					  Min, Max, Tol); 
	    Mvar2CtReduceRotExtremeDomain(Nodes, Cparam -> Right,
					  Min, Max, Tol);
	}
    }

    else {
	 /* We reached leaf node, update domain information. */
	for (i = 0; i < 4; ++i) {
	    if (Nodes[i] -> Min < Min[i])
		Min[i] = Nodes[i] -> Min;

	    if (Nodes[i] -> Max > Max[i])
		Max[i] = Nodes[i] -> Max;
	}

	if (Min[4] > Cparam -> TMin)
	    Min[4] = Cparam -> TMin;

	if (Max[4] < Cparam -> TMax)
	    Max[4] = Cparam -> TMax;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Extract sub-region of multivariate correspond to a given sub domain.     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:	 Array of multivariates.				     M
*   MVNum:	 Number of constraints (may be updated).		     M
*   Min:      Minimum domain values.                                         M
*   Max:      Maximum domain values.		                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct **:  Extracted MVs.                                         M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtExtractMVRegion                                                   M
*****************************************************************************/
MvarMVStruct **Mvar2CtExtractMVRegion(MvarMVStruct **MVs, 
				      int MVNum, 
				      CagdRType *Min, 
				      CagdRType *Max)
{
    int i, j;
    MvarMVStruct **SubMVs1, **SubMVs2;

    SubMVs1 = (MvarMVStruct**) IritMalloc(MVNum * sizeof(MvarMVStruct*));
    SubMVs2 = (MvarMVStruct**) IritMalloc(MVNum * sizeof(MvarMVStruct*));

    for (j = 0; j < MVNum; ++j) {
	SubMVs1[j] = MvarMVCopy(MVs[j]);	  
    }

    for (i=0; i < MVs[0] -> Dim; ++i) {
	for (j = 0; j < MVNum; ++j) {
	    /* Extract sub MVs. */
	    SubMVs2[j] = MvarMVRegionFromMV(SubMVs1[j], Min[i], Max[i], i);

	    MvarMVFree(SubMVs1[j]);

	    SubMVs1[j] = SubMVs2[j];
	}
    }

    IritFree(SubMVs1);

    return SubMVs2;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Check if a point is located in a given domain.                            M
*                                                                            *
* PARAMETERS:                                                                M
*   Min: Minimum values of the domain.                                       M
*   Max: Maximum values of the domain.                                       M
*   MPt: Point to check.                                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: TRUE if MPt is in the domain, FALSE otherwise.                M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtInDomain                                                          M
*****************************************************************************/
CagdBType Mvar2CtInDomain(CagdRType *Min, 
		          CagdRType *Max, 
		          MvarPtStruct *MPt)
{
    int i;

    for (i = 0; i < MPt -> Dim; ++i) {
	if (MPt -> Pt[i] < Min[i] || MPt -> Pt[i] >= Max[i])
	    return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Check if a polyline passes a given domain.                                M
*                                                                            *
* PARAMETERS:                                                                M
*   Min: Minimum values of the domain.                                       M
*   Max: Maximum values of the domain.                                       M
*   MPoly: MvarPolyline to check.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: TRUE if MPoly passes the domain, FALSE otherwise.             M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtIsPassing                                                         M
*****************************************************************************/
CagdBType Mvar2CtIsPassing(CagdRType *Min, 
		           CagdRType *Max, 
		           MvarPolylineStruct *MPoly)
{
    MvarPtStruct *TPt;

    for (TPt = MPoly -> Pl; TPt != NULL; TPt = TPt -> Pnext) {	
	if (Mvar2CtInDomain(Min, Max, TPt))
	    return TRUE;
    }

    return FALSE;

}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Check whether the 2 contact points cause inter-penetration into obstacle   M
* curves. The 2 contact points having inter-penetration are removed from the M
* list.                                                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   MPts:  Linked list of curvature contacts.                                M
*   BvhA:     Bvh for moving curve.                                   	     M
*   BvhBs:    Bvh for obstacle curves.                                       M
*   BSize:   Number of obstacle curves.                                      M
*   Circle:   Unit circle.             		                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *: Linked list of curvature contact points having no         M
*                  penetration into obstacle curves.                         M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtValidate2Ct                                                       M
*****************************************************************************/
MvarPtStruct *Mvar2CtValidate2Ct(MvarPtStruct *MPts,
			       	 Mvar2CtBVHStruct *BvhA,
				 Mvar2CtBVHStruct **BvhBs,
				 int BSize,
				 CagdCrvStruct *Circle)
{
    MvarPtStruct *TPt, 
        *PtOut = NULL;
    int First;
    CagdRType x, y, Theta, *R, penetration;
    CagdPType P, Q, C;
   
    while (MPts != NULL) {
	IRIT_LIST_POP(TPt, MPts);

	/* Get the index of obstacle curve. */
	First = AttrGetIntAttrib(TPt -> Attr, "First");

	R = CagdCrvEval(BvhA -> Crv, TPt -> Pt[0]);
	CagdCoerceToE2(P, &R, -1, BvhA -> Crv -> PType);

	R = CagdCrvEval(BvhBs[First] -> Crv, TPt -> Pt[1]);
	CagdCoerceToE2(Q, &R, -1, BvhBs[First] -> Crv -> PType);

	R = CagdCrvEval(Circle, TPt -> Pt[4]);
	CagdCoerceToE2(C, &R, -1, Circle -> PType);

	x = Q[0] - P[0] * C[0] + P[1] * C[1]; /* Compute x, y translation. */
	y = Q[1] - P[1] * C[0] - P[0] * C[1];

	Theta = Mvar2CtGetTheta(C[0], C[1]); /* Compute rotation angle. */

	/* Compute inter-penetration. */
	penetration = 
	    Mvar2CtPenetrationDepth(BvhA, BvhBs, BSize, x, y, Theta);
	if (penetration < - 2.0 * BvhA -> Tol) {
	    /* If it causes inter-penetration remove it from the list. */
	     MvarPtFree(TPt);
	}	
	else {
	IRIT_LIST_PUSH(TPt, PtOut);
	}
    }

    return PtOut;

}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Check whether the curvature contact points cause inter-penetration into    M
* obstacle curves. The curvature contact points having inter-penetration are M
* removed from the list.                                                     M
*                                                                            *
* PARAMETERS:                                                                M
*   MPts:  Linked list of curvature contacts.                                M
*   BvhA:     Bvh for moving curve.                                   	     M
*   BvhBs:    Bvh for obstacle curves.                                       M
*   BSize:   Number of obstacle curves.                                      M
*   BIndex:  Index number of obstacle curve.             		     M
*   Circle:   Unit circle.             		                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *: Linked list of curvature contact points having no        M
*                  penetration into obstacle curves.                         M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtValidateCurvContact                                               M
*****************************************************************************/
MvarPtStruct *Mvar2CtValidateCurvContact(MvarPtStruct *MPts,
					 Mvar2CtBVHStruct *BvhA,
					 Mvar2CtBVHStruct **BvhBs,
					 int BSize,
					 int BIndex,
					 CagdCrvStruct *Circle)
{
    MvarPtStruct *TPt, 
	*PtOut = NULL;
    CagdRType x, y, Theta, *R, penetration;
    CagdPType P, Q, C;

    while (MPts != NULL) {
	IRIT_LIST_POP(TPt, MPts);

	R = CagdCrvEval(BvhA -> Crv, TPt -> Pt[0]);
	CagdCoerceToE2(P, &R, -1, BvhA -> Crv -> PType);

	R = CagdCrvEval(BvhBs[BIndex] -> Crv, TPt -> Pt[1]);
	CagdCoerceToE2(Q, &R, -1, BvhBs[BIndex] -> Crv -> PType);

	R = CagdCrvEval(Circle, TPt -> Pt[2]);
	CagdCoerceToE2(C, &R, -1, Circle -> PType);

	x = Q[0] - P[0] * C[0] + P[1] * C[1];/* Compute x, y traslation. */
	y = Q[1] - P[1] * C[0] - P[0] * C[1];

	Theta = Mvar2CtGetTheta(C[0], C[1]); /* Compute rotation angle. */

	/* Compute inter-penetration. */
	penetration = 
	    Mvar2CtPenetrationDepth(BvhA, BvhBs, BSize, x, y, Theta);

	if (penetration < - 2.0 * BvhA -> Tol) { 
	    /* If there is inter-penetration remove it from the list. */
	    MvarPtFree(TPt);
	}
	else { 
	   IRIT_LIST_PUSH(TPt, PtOut);
	}
    }

    return PtOut;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Check if 2contact traces inter-penetrate the obstacle. The trace having    M
* inter-penetration is removed from the list.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Polys:    2contact traces.	                                             M
*   BvhA:     Bvh for moving curve.                                   	     M
*   BvhBs:    Bvh for obstacle curves.                                       M
*   BSize:    Number of obstacle curves.             		             M
*   Circle:   Unit circle.                                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: 2contact motion curves.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtValidateTraces                                                    M
*****************************************************************************/
MvarPolylineStruct *Mvar2CtValidateTraces(MvarPolylineStruct *Polys,
				          Mvar2CtBVHStruct *BvhA,
				          Mvar2CtBVHStruct **BvhBs,
				          int BSize,
				          CagdCrvStruct *Circle)
{
    MvarPolylineStruct *TPoly, 
	*Polyout = NULL;
    MvarPtStruct *Pts[3]; /* We check three points start middle end. */
    CagdPType P, Q, C;
    CagdRType Theta, x, y, penetration;
    int i, First;

    while (Polys != NULL) {
	IRIT_LIST_POP(TPoly, Polys);

	Pts[0] = TPoly -> Pl;
	Pts[1] = Mvar2CtGetMiddlePt(TPoly -> Pl, CagdListLength(TPoly -> Pl));
	Pts[2] = (MvarPtStruct *) CagdListLast(TPoly -> Pl);

	First = AttrGetIntAttrib(TPoly->Attr, "First");

	for (i = 0; i < 3; ++i) {
	    CagdRType *R = CagdCrvEval(BvhA -> Crv, Pts[i] -> Pt[0]);
	    CagdCoerceToE2(P, &R, -1, BvhA -> Crv -> PType);
	    
	    R = CagdCrvEval(BvhBs[First] -> Crv, Pts[i] -> Pt[1]);
	    CagdCoerceToE2(Q, &R, -1, BvhBs[First] -> Crv -> PType);

	    R = CagdCrvEval(Circle, Pts[i] -> Pt[4]);
	    CagdCoerceToE2(C, &R, -1, Circle -> PType);

	    x = Q[0] - P[0] * C[0] + P[1] * C[1];
	    y = Q[1] - P[1] * C[0] - P[0] * C[1];

	    Theta = Mvar2CtGetTheta(C[0], C[1]);

	    /* Check whether this trace interpenetrates. */
	    penetration = 
	          Mvar2CtPenetrationDepth(BvhA, BvhBs, BSize, 
	                                  x, y, Theta);

	    if (penetration < - 2.0 * BvhA -> Tol) {
		MvarPolylineFree(TPoly);
		TPoly = NULL;
		break;
	    }
	}

	if (TPoly != NULL)
	    IRIT_LIST_PUSH(TPoly, Polyout);
    }

    return Polyout;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Test whether two traces have common contact point.                         *
*                                                                            *
* PARAMETERS:                                                                *
*   SPt1, EPt1: Start, end point of first trace.                             *
*   SPt2, EPt2: Start, end point of second trace.                            * 
*   Type: Type of common contact point                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: FALSE if there is no chance to exist common contact point btw *
*             two traces and TRUE otherwise.                                 *
*****************************************************************************/
static CagdBType Mvar2CtTraceCollideAux(MvarPtStruct *SPt1, 
			                MvarPtStruct *EPt1, 
			                MvarPtStruct *SPt2, 
			                MvarPtStruct *EPt2, 
			                int Type)
{
    int i, 
	Index1 = 0, 
	Index2 = 0;
    CagdRType Min1[3], Max1[3], Min2[3], Max2[3];
    MvarPtStruct *TPt;

    for (i = 0; i < 3; ++i) {
	Min1[i] = Min2[i] = IRIT_INFNTY;
	Max1[i] = Max2[i] = -IRIT_INFNTY;
    }

    switch (Type) {
        case 0:
	    Index1 = 0;
	    Index2 = 0;
	    break;
        case 1:
	    Index1 = 0;
	    Index2 = 2;
	    break;
        case 2:
	    Index1 = 2;
	    Index2 = 0;
	    break;
        case 3:
	    Index1 = 2;
	    Index2 = 2;
	    break;
    }

    for (TPt = SPt1; TPt != EPt1 -> Pnext; TPt = TPt -> Pnext) {
	/* Compute bounding box of trace1. */
	if (TPt -> Pt[Index1] < Min1[0])
	    Min1[0] = TPt -> Pt[Index1];
	if (TPt -> Pt[Index1 + 1] < Min1[1])
	    Min1[1] = TPt -> Pt[Index1 + 1];
	if (TPt -> Pt[4] < Min1[2])
	    Min1[2] = TPt -> Pt[4];

	if (TPt -> Pt[Index1] > Max1[0])
	    Max1[0] = TPt -> Pt[Index1];
	if (TPt -> Pt[Index1 + 1] > Max1[1])
	    Max1[1] = TPt -> Pt[Index1 + 1];
	if (TPt -> Pt[4] > Max1[2])
	    Max1[2] = TPt -> Pt[4];
    }


    for (TPt = SPt2; TPt != EPt2 -> Pnext; TPt = TPt -> Pnext) {
	/* Compute bouding box of trace2. */
	if (TPt -> Pt[Index2] < Min2[0])
	    Min2[0] = TPt -> Pt[Index2];
	if (TPt -> Pt[Index2 + 1] < Min2[1])
	    Min2[1] = TPt -> Pt[Index2 + 1];
	if (TPt -> Pt[4] < Min2[2])
	    Min2[2] = TPt -> Pt[4];

	if (TPt -> Pt[Index2] > Max2[0])
	    Max2[0] = TPt -> Pt[Index2];
	if (TPt -> Pt[Index2 + 1] > Max2[1])
	    Max2[1] = TPt -> Pt[Index2 + 1];
	if (TPt -> Pt[4] > Max2[2])
	    Max2[2] = TPt -> Pt[4];
    }

    for (i = 0; i < 3; ++i) {
	/* If there is no overlap between bounding boxes                   */
	/* we don't have common contact point.                             */
	if (Min1[i] > Max2[i])
	    return FALSE;
	if (Min2[i] > Max1[i])
	    return FALSE;
    }

    if (SPt1 -> Pnext == EPt1 && SPt2 -> Pnext == EPt2)
	/* If traces get short enough return true. */
	return TRUE;
    else {
	/* Divide longer trace and repet the same procedure. */
	MvarPtStruct *Temp = EPt1 -> Pnext;
	MvarPtStruct *MidPt;
	int Length1, Length2;

	EPt1 -> Pnext = NULL;

	Length1 = CagdListLength(SPt1);

	EPt1 -> Pnext = Temp;

	Temp = EPt2 -> Pnext;
	EPt2 -> Pnext = NULL;

	Length2 = CagdListLength(SPt2);

	EPt2 -> Pnext = Temp;	

	if (Length1 > Length2) {
	    MidPt = Mvar2CtGetMiddlePt(SPt1, Length1);

	    return Mvar2CtTraceCollideAux(SPt1, MidPt, SPt2, EPt2, Type) ||  
		   Mvar2CtTraceCollideAux(MidPt, EPt1, SPt2, EPt2, Type);
	}

	else{
	    MidPt = Mvar2CtGetMiddlePt(SPt2, Length2);

	    return Mvar2CtTraceCollideAux(SPt1, EPt1, SPt2, MidPt, Type) ||  
		   Mvar2CtTraceCollideAux(SPt1, EPt1, MidPt, EPt2, Type);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Find a middle point of Point list.                                         M
*                                                                            *
* PARAMETERS:                                                                M
*   PtList:   Linked list of the points.                                     M
*   Length:   Length of the list.                                            M 
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   Middle point of the list.                              M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtGetMiddlePt                                                       M
*****************************************************************************/
MvarPtStruct *Mvar2CtGetMiddlePt(MvarPtStruct *PtList, int Length)
{
    int i,
	Mid = Length / 2;
    MvarPtStruct 
	*MidPt = PtList;

    for (i = 0; i < Mid; ++i)
	MidPt = MidPt -> Pnext;

    return MidPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Test whether two traces Poly1, Poly2 have common contact point.            M
*                                                                            *
* PARAMETERS:                                                                M
*   Poly1, Poly2: Two traces to check.                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   int: 0 if there is no common contact 1~4 depending on the type of common M
*        contact point.                                                      M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtTraceCollide                                                      M
*****************************************************************************/
int Mvar2CtTraceCollide(MvarPolylineStruct *Poly1, 
			MvarPolylineStruct *Poly2)
{
    MvarPtStruct *Last1, *Last2;
    int First1, Second1, First2, Second2;

    First1 = AttrGetIntAttrib(Poly1 -> Attr, "First");  
    Second1 = AttrGetIntAttrib(Poly1 -> Attr, "Second");
    First2 = AttrGetIntAttrib(Poly2 -> Attr, "First");
    Second2 = AttrGetIntAttrib(Poly2 -> Attr, "Second"); 

    Last1 = (MvarPtStruct *) CagdListLast(Poly1 -> Pl);
    Last2 = (MvarPtStruct *) CagdListLast(Poly2 -> Pl);

    if (First1 == First2) {
	if (Mvar2CtTraceCollideAux(Poly1 -> Pl, Last1, Poly2 -> Pl, Last2, 0))
	    return 1;
    }

    if (First1 == Second2) {
	if (Mvar2CtTraceCollideAux(Poly1 -> Pl, Last1, Poly2 -> Pl, Last2, 1))
	    return 2;
    }

    if (Second1 == First2) {
	if (Mvar2CtTraceCollideAux(Poly1 -> Pl, Last1, Poly2 -> Pl, Last2, 2))
	    return 3;
    }

    if (Second1 == Second2) {	    
	if (Mvar2CtTraceCollideAux(Poly1 -> Pl, Last1, Poly2 -> Pl, Last2, 3))
	    return 4;
    }

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Swap the first and second parameter with third and fourth value .          M
*                                                                            *
* PARAMETERS:                                                                M
*   MPoly: Trace to swap.                                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtSwapTrace                                                         M
*****************************************************************************/
void Mvar2CtSwapTrace(MvarPolylineStruct *MPoly)
{
    MvarPtStruct *TPt;
    int First, Second;

    First = AttrGetIntAttrib(MPoly -> Attr, "First");
    Second = AttrGetIntAttrib(MPoly -> Attr, "Second");

    for (TPt = MPoly -> Pl; TPt != NULL; TPt = TPt -> Pnext) {
	/* Swap the first and second parameter with third and fourth value. */
	CagdRType Temp;
	
	Temp = TPt -> Pt[0];
	TPt -> Pt[0] = TPt -> Pt[2]; 
	TPt -> Pt[2] = Temp;

	Temp = TPt -> Pt[1];
	TPt -> Pt[1] = TPt -> Pt[3];
	TPt -> Pt[3] = Temp;
    }

    AttrSetIntAttrib(&(MPoly -> Attr), "First", Second);/* Swap id. */
    AttrSetIntAttrib(&(MPoly -> Attr), "Second", First);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Compute bounding box of a trace.                                           M
*                                                                            *
* PARAMETERS:                                                                M
*   Min, Max: Minimum and maximum values of Bounding box.                    M
*   SPt, EPt:  Start and end point of the trace.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtTraceBBox                                                         M
*****************************************************************************/
void Mvar2CtTraceBBox(CagdRType *Min, 
		      CagdRType *Max, 
		      MvarPtStruct *SPt, 
		      MvarPtStruct *EPt)
{
    int i;
    MvarPtStruct *TPt;

    for (i = 0; i < SPt -> Dim; ++i) {
	Min[i] = IRIT_INFNTY;
	Max[i] = -IRIT_INFNTY;
    }

    for (TPt = SPt; TPt != EPt -> Pnext; TPt = TPt -> Pnext) {
	for (i = 0; i < SPt -> Dim; ++i) {
	    if (TPt -> Pt[i] < Min[i])
		Min[i] = TPt -> Pt[i];

	    if (TPt -> Pt[i] > Max[i])
		Max[i] = TPt -> Pt[i];
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Compute periodic difference between two numbers a, b. 0 and 1 are assumed  *
* to be identical.                                                           *
*                                                                            *
* PARAMETERS:                                                                *
*   a, b: two numbers for computing difference.                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: periodic difference between a and b.                          *
*****************************************************************************/
static CagdRType Mvar2CtPeriodicDiff(CagdRType a, CagdRType b)
{
    if (b > a)
	return IRIT_MIN(IRIT_FABS(b - a), IRIT_FABS(1 - b + a));
    else
	return IRIT_MIN(IRIT_FABS(a - b), IRIT_FABS(1 - a + b));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Compute periodic distance between two points Pt1, Pt2. The parameters are  * 
* assumed to be periodic and compute distance in periodic space.             *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt1, Pt2: two points for distance computation.                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: periodic distance between Pt1 and Pt2.                        *
*****************************************************************************/
static CagdRType Mvar2CtPeriodicDist(MvarPtStruct *Pt1, MvarPtStruct *Pt2)
{
    CagdRType Diff, 
        Dist = 0;
    int i;

    for (i = 0; i < 4; ++i) {
	 /* Periodic distance. */
	Diff = Mvar2CtPeriodicDiff(Pt1 -> Pt[i], Pt2 -> Pt[i]);

	Dist += Diff*Diff;
    }
#ifdef MVAR2CT_CONNECT_ANGLE_PERIODIC
    Diff = Mvar2CtPeriodicDiff(Pt1 -> Pt[i], Pt2 -> Pt[i]);
#else
    Diff = IRIT_FABS(Pt1->Pt[4]-Pt2->Pt[4]);
#endif

    Dist += Diff*Diff;

    return sqrt(Dist);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Check if two points are on the same periodic boundary                      *
*  .                                                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt1, Pt2: two points for checking.                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if two points are on the same boundary FALSE otherwise.  *
*****************************************************************************/
static CagdBType Mvar2CtIsSameBndry(MvarPtStruct *Pt1, MvarPtStruct *Pt2)
{
    int i;

    for (i = 0; i < 5; ++i) {
	if (IRIT_FABS(IRIT_FABS(Pt1 -> Pt[i]-Pt2 -> Pt[i]) - 1) < 
	   MVAR2CT_PERIODIC_TOL)
	   return TRUE;
    }

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Connect 2Contact traces that are not connected due to the parametrization  *
* of the periodic curve.                                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly1, Poly2: two 2contact traces for checking the connectivity between  *
*   them.	                                                             *
*   Tol:   tolerance for the connection.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: 0 for disconnected traces and 1~8 for connected traces depending on *
*        how they connected.                                                 *
*****************************************************************************/
static int Mvar2CtIsConnectedTrace(MvarPolylineStruct *Poly1, 
				   MvarPolylineStruct *Poly2, 
				   CagdRType Tol)
{
    IRIT_STATIC_DATA MvarPtStruct 
	*SwapFirst = NULL, 
	*SwapLast = NULL;
    MvarPtStruct *Last1, *Last2;
    int RetVal, First1, First2, Second1, Second2;

    /* Index of first obstacle crv. */
    First1 = AttrGetIntAttrib(Poly1 -> Attr, "First"); 
    First2 = AttrGetIntAttrib(Poly2 -> Attr, "First");

    /* Index of second obstacle crv. */
    Second1 = AttrGetIntAttrib(Poly1 -> Attr, "Second"); 
    Second2 = AttrGetIntAttrib(Poly2 -> Attr, "Second");

    /* Check if two traces comes from same curve. */
    if (!(First1 == First2 && Second1 == Second2) &&
	!(First1 == Second2 && Second1 == First2)) 
	return 0;

    if (SwapLast == NULL) {
	SwapLast = MvarPtNew(5);
	SwapFirst = MvarPtNew(5);
    }

    Last1 = (MvarPtStruct *) CagdListLast(Poly1 -> Pl);
    Last2 = (MvarPtStruct *) CagdListLast(Poly2 -> Pl);

    /* Compute distance between endpoints. */
    if (Mvar2CtPeriodicDist(Poly1 -> Pl, Poly2 -> Pl) < Tol && 
	Mvar2CtIsSameBndry(Poly1 -> Pl, Poly2 -> Pl)) 
	return 1;
    if (Mvar2CtPeriodicDist(Last1, Last2) < Tol && 
	Mvar2CtIsSameBndry(Last1, Last2))
	return 2;
    if (Mvar2CtPeriodicDist(Last1, Poly2 -> Pl) < Tol &&
	Mvar2CtIsSameBndry(Last1, Poly2 -> Pl))
	return 3;
    if (Mvar2CtPeriodicDist(Poly1 -> Pl, Last2) < Tol &&
	Mvar2CtIsSameBndry(Poly1 -> Pl, Last2))
	return 4;

    /* We also need to check for the case of swapped trace. */
    SwapFirst -> Pt[0] = Poly2 -> Pl -> Pt[2];  
    SwapFirst -> Pt[1] = Poly2 -> Pl -> Pt[3]; 
    SwapFirst -> Pt[2] = Poly2 -> Pl -> Pt[0];
    SwapFirst -> Pt[3] = Poly2 -> Pl -> Pt[1];
    SwapFirst -> Pt[4] = Poly2 -> Pl -> Pt[4];

    SwapLast -> Pt[0] = Last2 -> Pt[2];
    SwapLast -> Pt[1] = Last2 -> Pt[3];
    SwapLast -> Pt[2] = Last2 -> Pt[0];
    SwapLast -> Pt[3] = Last2 -> Pt[1];
    SwapLast -> Pt[4] = Last2 -> Pt[4];

    if (Mvar2CtPeriodicDist(Poly1 -> Pl, SwapFirst) < Tol &&
	Mvar2CtIsSameBndry(Poly1 -> Pl, SwapFirst))
        RetVal = 5;
    else if (Mvar2CtPeriodicDist(Last1, SwapLast) < Tol &&
	     Mvar2CtIsSameBndry(Last1, SwapLast))
        RetVal = 6;
    else if (Mvar2CtPeriodicDist(Last1, SwapFirst) < Tol &&
	     Mvar2CtIsSameBndry(Last1, SwapFirst))
        RetVal = 7;
    else if (Mvar2CtPeriodicDist(Poly1 -> Pl, SwapLast) < Tol &&
	     Mvar2CtIsSameBndry(Poly1 -> Pl, SwapLast))
        RetVal = 8;
    else
        RetVal = 0;

    MvarPtFree(SwapFirst);
    MvarPtFree(SwapLast);

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Connect 2Contact traces that are not connected due to the parametrization  M
* of the periodic curve.                                                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Polys:    2contact traces.	                                             M
*   Tol:   tolerance for the connection.                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: connected 2contact motion curves.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CtConnectPeriodic                                                   M
*****************************************************************************/
MvarPolylineStruct *Mvar2CtConnectPeriodic(MvarPolylineStruct *Polys, 
				           CagdRType Tol)
{
    MvarPolylineStruct *TPoly, *TPoly2,
        *Polyout = NULL;
    MvarPtStruct *Reverse;

    while (Polys != NULL) {
	IRIT_LIST_POP(TPoly, Polys);

	for (TPoly2 = Polys; TPoly2 != NULL; TPoly2 = TPoly2 -> Pnext) {
	     /* Check if TPoly1 and TPoly2 are connected. */
	    int Type = Mvar2CtIsConnectedTrace(TPoly, TPoly2, Tol);

	    switch (Type) { /* Connect TPoly1 and TPoly2. */ 
	        case 1:
		    Reverse = (MvarPtStruct *) CagdListReverse(TPoly -> Pl);
		
		    TPoly2 -> Pl = 
		        (MvarPtStruct*) CagdListAppend(Reverse, TPoly2 -> Pl);
		
		    TPoly -> Pl=NULL;
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;
		    break;
	        case 2:
		    Reverse = (MvarPtStruct *) CagdListReverse(TPoly -> Pl);
		
		    TPoly2 -> Pl = 
		         (MvarPtStruct*) CagdListAppend(TPoly2 -> Pl, Reverse);
		
		    TPoly -> Pl = NULL;
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;		
		    break;
	        case 3:
		    TPoly2 -> Pl = 
		    (MvarPtStruct*) CagdListAppend(TPoly -> Pl, TPoly2 -> Pl);
		    TPoly -> Pl = NULL;
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;
		    break;
	        case 4:
		    TPoly2 -> Pl = 
		    (MvarPtStruct*) CagdListAppend(TPoly2 -> Pl, TPoly -> Pl);
		    TPoly -> Pl = NULL;	
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;
		    break;
	        case 5:
		    Mvar2CtSwapTrace(TPoly);
		    Reverse = (MvarPtStruct *) CagdListReverse(TPoly -> Pl);
		    TPoly2 -> Pl = 
		      (MvarPtStruct*) CagdListAppend(Reverse, TPoly2 -> Pl);
		    TPoly -> Pl = NULL;
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;
		    break;
	        case 6:
		    Mvar2CtSwapTrace(TPoly);
		    Reverse = (MvarPtStruct *) CagdListReverse(TPoly -> Pl);
		    TPoly2 -> Pl = 
		      (MvarPtStruct*) CagdListAppend(TPoly2 -> Pl, Reverse);
		    TPoly -> Pl = NULL;
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;		
		    break;
	        case 7:
		    Mvar2CtSwapTrace(TPoly);
		    TPoly2 -> Pl = 
		    (MvarPtStruct*) CagdListAppend(TPoly -> Pl, TPoly2 -> Pl);
	 	    TPoly -> Pl = NULL;
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;
		    break;
	        case 8:
		    Mvar2CtSwapTrace(TPoly);
		    TPoly2 -> Pl = 
		    (MvarPtStruct*) CagdListAppend(TPoly2 -> Pl, TPoly -> Pl);
		    TPoly -> Pl = NULL;		
		    MvarPolylineFree(TPoly);
		    TPoly = NULL;		
		    break;
	    }

	    if (TPoly == NULL)
		break;
	}

	if (TPoly != NULL)
	    IRIT_LIST_PUSH(TPoly, Polyout);
    }

    return Polyout;
}
