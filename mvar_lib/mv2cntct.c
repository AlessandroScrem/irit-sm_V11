/******************************************************************************
* mv2cntct.c - 2contact motion analysis for planar curves.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yong-Joon Kim, July 2014.					      *
******************************************************************************/

#include "inc_irit/misc_lib.h"
#include "mvar_loc.h"

#define MVAR2CT_BV_TOL		1e-6
#define MVAR2CT_BV_SUBTOL       1e-4
#define MVAR2CT_PAIR_TOL	1e-2

typedef struct Mvar2CtPairStruct {       /* BVNode Pair for Domain Purging. */
    Mvar2CtAABBStruct AABB;
    Mvar2CtBVNodeStruct *ANode, *BNode;
    MvarPtStruct *CurvContacts;
    Mvar2CtCParamStruct *Cparam;
    CagdBType Bsingle;    
    struct Mvar2CtPairStruct *Left, *Right;
} Mvar2CtPairStruct;

static int Mvar2CtPairSubdivDir(Mvar2CtPairStruct *Node1, CagdRType tol);
static int Mvar2CtPairSubdivDir2(Mvar2CtPairStruct *Node1,
				 Mvar2CtPairStruct *Node2);
static void Mvar2CtGetPairSubdivParam(Mvar2CtPairStruct *Node1, 
				      Mvar2CtPairStruct *Node2, 
				      int *Dir,
				      CagdRType *Param);
static CagdBType Mvar2CtAABBOverlap(Mvar2CtAABBStruct *A,
				    Mvar2CtAABBStruct *B);
static void Mvar2CtPairAABB(Mvar2CtBVNodeStruct *ANode, 
			    Mvar2CtBVNodeStruct *BNode, 
			    CagdRType RMin, 
			    CagdRType RMax, 
			    Mvar2CtAABBStruct *AABB);
static void Mvar2CtBuildPairHierarchy(Mvar2CtPairStruct **Node, 
				      Mvar2CtBVNodeStruct *ANode, 
				      Mvar2CtBVNodeStruct *BNode, 
				      Mvar2CtBVHStruct *BvhA, 
				      Mvar2CtBVHStruct **BvhBs,  
				      int BSize,
				      Mvar2CtCParamStruct *Cparam, 
				      CagdRType Tol);
static CagdBType Mvar2CtAddChildPair(Mvar2CtPairStruct *Node, 
				     CagdRType Tol);
static void Mvar2CtFreePair(Mvar2CtPairStruct *Node);
static CagdBType Mvar2CtTrivialReject(Mvar2CtPairStruct *Node1,
				      Mvar2CtPairStruct *Node2);
static CagdBType Mvar2CtIsConnectedPair(Mvar2CtPairStruct *Node1,
					Mvar2CtPairStruct *Node2);

static MvarPtStruct *Mvar2CtSolveRotExtreme(Mvar2CtPairStruct *Node1,
					    Mvar2CtPairStruct *Node2,
					    Mvar2CtBVHStruct *BvhA,
					    Mvar2CtBVHStruct **BvhBs,
					    CagdRType Subtol,
					    CagdRType Numerictol);
static void Mvar2CtRotationExtreme(Mvar2CtPairStruct *Node1,
				   Mvar2CtPairStruct *Node2,
				   Mvar2CtBVHStruct *BvhA,
				   Mvar2CtBVHStruct **BvhBs,
				   MvarPtStruct **MPts,
				   CagdRType Subtol,
				   CagdRType Numerictol);
static void Mvar2CtComputeCurvatureContact(Mvar2CtPairStruct *Node,
					   Mvar2CtBVHStruct *BvhA,
					   Mvar2CtBVHStruct *BvhB,
					   MvarPtStruct **CurvContacts,
					   CagdRType Subtol,
					   CagdRType Numerictol);
static MvarPolylineStruct * Mvar2CtSolveCurvatureTrace(Mvar2CtPairStruct *Node,
						       Mvar2CtBVHStruct *BvhA,
						       Mvar2CtBVHStruct *BvhB,
						       CagdRType Subtol,
						       CagdRType Numerictol);
static void Mvar2CtSetSingleRecursive(Mvar2CtPairStruct *Node);
static CagdBType Mvar2CtInPairDomain(Mvar2CtPairStruct *Node1, 
				     Mvar2CtPairStruct *Node2, 
				     MvarPtStruct *MPt);
static void Mvar2CtCheckSingleContact(Mvar2CtPairStruct *Node,
				      MvarPtStruct *CurvContact,
				      MvarPolylineStruct *CurvTrace,
				      CagdRType Tol);
static void Mvar2CtSetSingleContact(Mvar2CtPairStruct *Node,
				    Mvar2CtBVHStruct *BvhA,
				    Mvar2CtBVHStruct *BvhB,
				    MvarPtStruct *CurvContact,
				    CagdRType Subtol,
				    CagdRType Numerictol);
static MvarPtStruct *Mvar2CtBndrySolution(MvarMVStruct **MVs, 
					  MvarConstraintType *Constr, 
					  int MV_Num, 
					  Mvar2CtBVNodeStruct *Nodes[4], 
					  Mvar2CtCParamStruct *Cparam, 
					  CagdRType SubdivTol, 
					  CagdRType NumericTol);
static MvarPolylineStruct *Mvar2CtSolve2CtTrace(Mvar2CtPairStruct *Node1,
				                Mvar2CtPairStruct *Node2,
				                Mvar2CtBVHStruct *BvhA,
				                Mvar2CtBVHStruct **BvhBs,
				                CagdBType NoLoop,
						CagdRType Step,
				                CagdRType Subtol, 
				                CagdRType Numerictol);
static MvarPtStruct *Mvar2CtSolve3CtPoint(Mvar2CtBVHStruct *BvhA,
					  Mvar2CtBVHStruct **BvhBs,
					  int Ids[3],
					  CagdRType Min[7],
					  CagdRType Max[7],
					  CagdRType Subtol, 
					  CagdRType Numerictol);
static CagdBType Mvar2CtNoLoopTest(Mvar2CtPairStruct *Node1,
				   Mvar2CtPairStruct *Node2,
				   MvarPtStruct *RotExtremes);
static MvarPolylineStruct* Mvar2CtCompute2CtTrace(Mvar2CtPairStruct *Node1,
						  Mvar2CtPairStruct *Node2,
						  Mvar2CtBVHStruct *BvhA,
						  Mvar2CtBVHStruct **BvhBs,
						  MvarPtStruct *RotExtremes,
						  CagdRType Step,
						  CagdRType Subtol,
						  CagdRType Numerictol);

static void Mvar2CtCompute3CtPoint(Mvar2CtBVHStruct *BvhA, 
				   Mvar2CtBVHStruct **BvhBs, 
				   MvarPtStruct *SPts[3], 
				   MvarPtStruct *EPts[3],
				   Mvar2CtCParamStruct *CRoot,
				   int Ids[3], 
				   MvarPtStruct **Triples, 
				   CagdRType Subtol, 
				   CagdRType Numerictol);
static void Mvar2CtInsert3CtPoint(MvarPolylineStruct *Poly, MvarPtStruct *Pt);
static void Mvar2CtInsertCurvPoint(MvarPolylineStruct *Poly, MvarPtStruct *Pt);

static MvarPolylineStruct *Mvar2CtDivideAtCurvPointAux(MvarPolylineStruct *Polys);
static MvarPolylineStruct *Mvar2CtDivideAt3CtPointAux(MvarPolylineStruct
						                      *Polys);
static MvarPolylineStruct* Mvar2CtDivideAt3CtPoint(Mvar2CtBVHStruct *BvhA,
						   Mvar2CtBVHStruct **BvhBs,
						   Mvar2CtCParamStruct *CRoot,
						   MvarPolylineStruct
						                  *M2CtTraces,
						   CagdRType Subtol,
						   CagdRType Numerictol);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Perform AABB overlap test.                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   A, B:  AABBs for overlap test.                                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if overlap FALSE otherwise.                              *
*****************************************************************************/
static CagdBType Mvar2CtAABBOverlap(Mvar2CtAABBStruct *A, 
				    Mvar2CtAABBStruct *B)
{
    if (A -> Xmax < B -> Xmin)
	return FALSE;
    if (A -> Xmin > B -> Xmax)
	return FALSE;
    if (A -> Ymax < B -> Ymin)
	return FALSE;
    if (A -> Ymin > B -> Ymax)
	return FALSE;

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Generate a AABB for a given pair of curves.                              *
*                                                                            *
* PARAMETERS:                                                                *
*   ANode:  Bounding volume for moving curve.                                *
*   BNode:  Bounding volume for obstacle curves.                             *
*   RMin:  Minimum rotation.                                                 *
*   RMax:  Maximum rotation.                                                 *
*   AABB: Generated AABB will be placed herein.                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtPairAABB(Mvar2CtBVNodeStruct *ANode, 
		            Mvar2CtBVNodeStruct *BNode, 
		            CagdRType RMin, 
		            CagdRType RMax, 
		            Mvar2CtAABBStruct *AABB)
{
  
    CagdRType Rx1, Ry1, Rx2, Ry2, X, Y;
    /* Error of bounding volumes and rotation. */
    CagdRType
        Error = ANode -> LBV.Epsilon + BNode -> LBV.Epsilon +
	        ANode -> Length * (1 - cos((RMax - RMin) * 0.5));
    int i,j;

    AABB -> Xmin = AABB -> Ymin = IRIT_INFNTY;
    AABB -> Xmax = AABB -> Ymax = -IRIT_INFNTY;

    for (i = 0; i < 2; ++i) {
	/* Compute minkowski sum of LSCs. */ 
	Rx1 = cos(RMin) * ANode -> LBV.P[i][0] - 
	      sin(RMin) * ANode -> LBV.P[i][1];
	Ry1 = cos(RMin) * ANode -> LBV.P[i][1] + 
	      sin(RMin) * ANode -> LBV.P[i][0];

	Rx2 = cos(RMax) * ANode -> LBV.P[i][0] - 
	      sin(RMax) * ANode -> LBV.P[i][1];
	Ry2 = cos(RMax) * ANode -> LBV.P[i][1] + 
	      sin(RMax) * ANode -> LBV.P[i][0];

	for ( j = 0; j < 2; ++j) {
	    X = BNode -> LBV.P[j][0] - Rx1;
	    Y = BNode -> LBV.P[j][1] - Ry1;

	    if (X > AABB -> Xmax)
	        AABB -> Xmax = X;
	    if (X < AABB -> Xmin)
	        AABB -> Xmin = X;

	    if (Y > AABB -> Ymax)
	        AABB -> Ymax = Y;
	    if (Y < AABB -> Ymin)
	        AABB -> Ymin = Y;

	    X = BNode -> LBV.P[j][0] - Rx2;
	    Y = BNode -> LBV.P[j][1] - Ry2;

	    if (X > AABB -> Xmax)
	        AABB -> Xmax = X;
	    if (X < AABB -> Xmin)
	        AABB -> Xmin = X;

	    if (Y > AABB -> Ymax)
		AABB -> Ymax = Y;
	    if (Y < AABB -> Ymin)
		AABB -> Ymin = Y;
	}
    }
    /* Enlarge the AABB by error. */
    AABB -> Xmax += Error;
    AABB -> Xmin -= Error;
    AABB -> Ymax += Error;
    AABB -> Ymin -= Error;
    AABB -> Cx = (AABB -> Xmin + AABB -> Xmax) * 0.5;
    AABB -> Cy = (AABB -> Ymin + AABB -> Ymax) * 0.5;

    AABB -> Radius = sqrt((AABB -> Xmax - AABB -> Cx) * 
	                  (AABB -> Xmax - AABB -> Cx) +
	                  (AABB -> Ymax - AABB -> Cy) * 
			  (AABB -> Ymax - AABB -> Cy));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Find a subdivision direction for pair of curves.                          *
*                                                                            *
* PARAMETERS:                                                                *
*   Node: Pair of curves for first contact point.                            *
*   Tol:  Tolerance of subdivision.                                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: Direction for next subdivision.                                     *
*****************************************************************************/
static int Mvar2CtPairSubdivDir(Mvar2CtPairStruct *Node, CagdRType Tol)
{
    CagdRType Length1, Length2, Length3;

    Length1 = Node -> ANode -> Max - Node -> ANode -> Min;
    Length2 = Node -> BNode -> Max - Node -> BNode -> Min;
    Length3 = Node -> Cparam -> TMax - Node -> Cparam -> TMin;

    if (IRIT_MAX3(Length1, Length2, Length3) < Tol 
        && Node -> ANode -> Convexity * Node -> BNode -> Convexity != 0)
	/* If curves all small enough. */
	return 3;

    if (Length1 >= Length2 && Length1 >= Length3) 
	return 0;
    else if (Length2 >= Length3)
	return 1;
    else
	return 2;    
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Find a subdivision direction for two pairs of curves.                     *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1: Pair of curves for first contact point.                           *
*   Node2: Pair of curves for second contact point.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: direction for next subdivision.                                     *
*****************************************************************************/
static int Mvar2CtPairSubdivDir2(Mvar2CtPairStruct *Node1, 
				 Mvar2CtPairStruct *Node2)
{
    CagdRType Length1, Length2;
    CagdBType Bleaf1, Bleaf2;

    Bleaf1 = Node1 -> Left == NULL && Node1 -> Right == NULL;
    Bleaf2 = Node2 -> Left == NULL && Node2 -> Right == NULL;

    if (Bleaf1 && Bleaf2)
	/* If we reached leaf node. */
	return 2;
    else if (Bleaf1 && !Bleaf2)
	/* Node1 is already leaf. */
	return 1;
    else if (!Bleaf1 && Bleaf2)
	/* Node2 is already leaf. */
	return 0;
    else {
	Length1 = Node1 -> ANode -> Max - Node1 -> ANode -> Min 
	        + Node1 -> BNode -> Max - Node1 -> BNode -> Min 
		+ Node1 -> Cparam -> TMax - Node1 -> Cparam -> TMin; 

	Length2 = Node2 -> ANode -> Max - Node2 -> ANode -> Min 
		+ Node2 -> BNode -> Max - Node2 -> BNode -> Min 
		+ Node2 -> Cparam -> TMax - Node2 -> Cparam -> TMin; 

	if (Length1 >= Length2)
	    /* If Node1 has larger domain. */
	    return 0;
	else
	    return 1;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Find a subdivision direction and subdivision parameter for two pairs of   *
*  curves.                                                                   *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1: Pair of curves for first contact point.                           *
*   Node2: Pair of curves for second contact point.                          *
*   Dir: Direction of subdivision will be placed herein.                     *
*   Param: Parameter of subdivision will be placed herein.                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtGetPairSubdivParam(Mvar2CtPairStruct *Node1, 
			              Mvar2CtPairStruct *Node2, 
			              int *Dir,
			              CagdRType *Param)
{
    int Dir1, Dir2;

    /* Check which pair have larger domain. */
    Dir1 = Mvar2CtPairSubdivDir2(Node1, Node2);
    
    if (Dir1 == 0) {
	/* Check which curve (moving or obstacle) or rotation have         */
	/* larger domain.                                                  */
	Dir2 = Mvar2CtPairSubdivDir(Node1, 0);	
	switch (Dir2) {
	    case 0:
	        *Dir = 0;
	    
		if (Node1 -> Left != NULL)
		    *Param = Node1 -> Left -> ANode -> Max;
		else
		    *Param = Node1 -> Right -> ANode -> Min;	    
		return;
	    case 1:
	        *Dir = 1;
	    
		if (Node1 -> Left != NULL)
		    *Param = Node1 -> Left -> BNode -> Max;
		else
		    *Param = Node1 -> Right -> BNode -> Min;
		return;
	    case 2:
	        *Dir = 4;
	    
		if (Node1 -> Left != NULL)
		    *Param = Node1 -> Left -> Cparam -> TMax;	    
		else
		    *Param = Node1 -> Right -> Cparam -> TMin;
		return;
	}
    }
    else {
	/* Check which curve (moving or obstacle) or rotation have         */ 
	/* larger domain.                                                  */
	Dir2 = Mvar2CtPairSubdivDir(Node2, 0);

	switch(Dir2) {
	    case 0:
	        *Dir = 2;
	        if (Node2 -> Left != NULL)
		    *Param = Node2 -> Left -> ANode -> Max;
	        else
	 	    *Param = Node2 -> Right -> ANode -> Min;
	        return;
	    case 1:
	        *Dir = 3;
	        if (Node2 -> Left != NULL)
		    *Param = Node2 -> Left -> BNode -> Max;

	        else
		    *Param = Node2 -> Right -> BNode -> Min;
	        return;
	    case 2:
	        *Dir = 4;
	        if (Node2 -> Left != NULL)
		    *Param = Node2 -> Left -> Cparam -> TMax;
	        else
		    *Param = Node2 -> Right -> Cparam -> TMin;
	        return;
	}

    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Build a hierarchy for pair of curves recursively. For each pair we        *
*  compute maximum penetration depth at middle parameter and see if we can   *
*   purge corresponding domain.                                              *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:   New node will be allocated herein.                               *
*   ANode:  Bounding volume of moving curve.                                 *
*   BNode:  Bounding volume of obstacle curve.                               *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhBs:  Bvh for obstacle curves.                                         *
*   BSize:  Number of obstacle curves.                                       *
*   Cparam: Data structure for rotation.                                     *
*   Tol: Tolerance for allocation.                                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtBuildPairHierarchy(Mvar2CtPairStruct **Node, 
			              Mvar2CtBVNodeStruct *ANode, 
			              Mvar2CtBVNodeStruct *BNode, 
			              Mvar2CtBVHStruct *BvhA, 
			              Mvar2CtBVHStruct **BvhBs,  
			              int BSize,
			              Mvar2CtCParamStruct *Cparam, 
			              CagdRType Tol)
{  
    Mvar2CtAABBStruct AABB;
    CagdRType penetration;
    int Dir;
  
    if (!Mvar2CtNormalOverlap(ANode, BNode, Cparam -> RMin, Cparam -> RMax)) {
	/* Check tangency condition. */
	(*Node) = NULL;
	return;
    }
 
    /* Generate AABB for pair of curves. */
    Mvar2CtPairAABB(ANode, BNode, Cparam -> RMin, Cparam -> RMax, &AABB);

    /* Compute maximum penetration depth at middle transformation. */
    penetration =
        Mvar2CtPenetrationDepth(BvhA, BvhBs, BSize, AABB.Cx, AABB.Cy, 
				(Cparam -> RMin + Cparam -> RMax) * 0.5);
    /* See if we can purge domain. */ 
    if (penetration + AABB.Radius + BvhA -> Radius * 
	                     (Cparam -> RMax - Cparam -> RMin) * 0.5 > 0) {
	
	(*Node) = (Mvar2CtPairStruct *) IritMalloc(sizeof(Mvar2CtPairStruct));

	(*Node) -> ANode = ANode;
	(*Node) -> BNode = BNode;
	(*Node) -> AABB = AABB;
	(*Node) -> Cparam = Cparam;
	(*Node) -> Left = NULL;
	(*Node) -> Right = NULL;
	(*Node) -> Bsingle = FALSE;
	(*Node) -> CurvContacts = NULL;
	
	Dir = Mvar2CtPairSubdivDir((*Node), Tol);

	switch (Dir) {
	    case 0:
	        if (ANode -> Left != NULL) {	    
		    Mvar2CtBuildPairHierarchy(&((*Node) -> Left), 
			                      ANode -> Left, BNode, BvhA, 
					      BvhBs, BSize, Cparam, Tol);
		
		    Mvar2CtBuildPairHierarchy(&((*Node) -> Right), 
		                              ANode -> Right, BNode, BvhA, 
					      BvhBs, BSize, Cparam, Tol);
		
		    if ((*Node) -> Left == NULL && (*Node) -> Right == NULL) {
		        IritFree(*Node);
			(*Node) = NULL;
		    }
		}
	        break;
	    case 1:
	        if (BNode -> Left != NULL) {
		    Mvar2CtBuildPairHierarchy(&((*Node) -> Left), ANode, 
		                              BNode -> Left, BvhA, BvhBs, 
					      BSize, Cparam, Tol);
		    Mvar2CtBuildPairHierarchy(&((*Node) -> Right), ANode, 
		                              BNode -> Right, BvhA, BvhBs, 
					      BSize, Cparam, Tol);

		    if ((*Node) -> Left == NULL && (*Node) -> Right == NULL) {
		        IritFree(*Node);
			(*Node) = NULL;
		    }
		}
	        break;
	    case 2:
	        if (Cparam -> Left != NULL) {	   
		    Mvar2CtBuildPairHierarchy(&((*Node) -> Left), ANode, 
		                              BNode, BvhA, BvhBs, BSize, 
					      Cparam -> Left, Tol);
		    Mvar2CtBuildPairHierarchy(&((*Node) -> Right), ANode, 
		                              BNode, BvhA, BvhBs, BSize, 
					      Cparam -> Right, Tol);

		    if ((*Node) -> Left == NULL && (*Node) -> Right == NULL) {
		        IritFree(*Node);
			(*Node) = NULL;
		    }
		}
	        break;
	    case 3:
	    default:
	        return;
	}
    }
    else {
        (*Node) = NULL;
	return;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Allocate child node to pair struct.                                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Node: Mvar2CtPairStruct node to add child nodes.                         *
*   Tol: Tolerance for allocation.                                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if child nodes are allocated, FALSE otherwise.           *
*****************************************************************************/
static CagdBType Mvar2CtAddChildPair(Mvar2CtPairStruct *Node, 
			             CagdRType Tol)
{
    Mvar2CtPairStruct *Left, *Right;
    int Dir;

    /* Decide direction for subdivision. */
    Dir = Mvar2CtPairSubdivDir(Node, Tol);

    switch (Dir) {
        case 0:
	    Left = (Mvar2CtPairStruct*) IritMalloc(sizeof(Mvar2CtPairStruct));

	    Left -> ANode = Node -> ANode -> Left;
	    Left -> BNode = Node -> BNode;
	    Left -> Cparam = Node -> Cparam;
	    Mvar2CtPairAABB(Left -> ANode, Left -> BNode, 
	                    Left -> Cparam -> RMin, Left -> Cparam -> RMax,
			    &(Left -> AABB));
	    Left -> Left = NULL;
	    Left -> Right = NULL;
	    Left -> Bsingle = FALSE;
	    Left -> CurvContacts = NULL;

	    Right = 
                (Mvar2CtPairStruct*) IritMalloc(sizeof(Mvar2CtPairStruct));

	    Right -> ANode = Node -> ANode -> Right;
	    Right -> BNode = Node -> BNode;
	    Right -> Cparam = Node -> Cparam;
	
	    Mvar2CtPairAABB(Right -> ANode, Right -> BNode, 
	                    Right -> Cparam -> RMin, Right -> Cparam -> RMax,
	                    &(Right -> AABB));
	    Right -> Left = NULL;
	    Right -> Right = NULL;
	    Right -> Bsingle = FALSE;
	    Right -> CurvContacts = NULL;

	    Node -> Left = Left;
	    Node -> Right = Right;

	    return TRUE;
        case 1:
	    Left = 
	        (Mvar2CtPairStruct*) IritMalloc(sizeof(Mvar2CtPairStruct));
	    Left -> ANode = Node -> ANode;
	    Left -> BNode = Node -> BNode -> Left;
	    Left -> Cparam = Node -> Cparam;
	    Mvar2CtPairAABB(Left -> ANode, Left -> BNode, 
	                    Left -> Cparam -> RMin, Left -> Cparam -> RMax,
	                    &(Left -> AABB));
 	    Left -> Left = NULL;
	    Left -> Right = NULL;
	    Left -> Bsingle = FALSE;
	    Left -> CurvContacts = NULL;

	    Right = 
		(Mvar2CtPairStruct*) IritMalloc(sizeof(Mvar2CtPairStruct));

	    Right -> ANode = Node -> ANode;
	    Right -> BNode = Node -> BNode -> Right;
	    Right -> Cparam = Node -> Cparam;
	    Mvar2CtPairAABB(Right -> ANode, Right -> BNode, 
	                    Right -> Cparam -> RMin, Right -> Cparam -> RMax,
	                    &(Right -> AABB));
 	    Right -> Left = NULL;
	    Right -> Right = NULL;
	    Right -> Bsingle = FALSE;
	    Right -> CurvContacts = NULL;

	    Node -> Left = Left;
	    Node -> Right = Right;

	    return TRUE;	
        case 2:
	    Left = 
		(Mvar2CtPairStruct*) IritMalloc(sizeof(Mvar2CtPairStruct));
	    Left -> ANode = Node -> ANode;
	    Left -> BNode = Node -> BNode;
	    Left -> Cparam = Node -> Cparam -> Left;
	    
	    Mvar2CtPairAABB(Left -> ANode, Left -> BNode, 
	                    Left -> Cparam -> RMin, Left -> Cparam -> RMax,
	                    &(Left -> AABB));
	    Left -> Left = NULL;
	    Left -> Right = NULL;
	    Left -> Bsingle = FALSE;
	    Left -> CurvContacts = NULL;

	    Right = 
		(Mvar2CtPairStruct*) IritMalloc(sizeof(Mvar2CtPairStruct));

	    Right -> ANode = Node -> ANode;
	    Right -> BNode = Node -> BNode;
	    Right -> Cparam = Node -> Cparam -> Right;
	
	    Mvar2CtPairAABB(Right -> ANode, Right -> BNode, 
	                    Right -> Cparam -> RMin, Right -> Cparam -> RMax,
	                    &(Right -> AABB));
	    Right -> Left = NULL;
	    Right -> Right = NULL;
	    Right -> Bsingle = FALSE;
	    Right -> CurvContacts = NULL;

	    Node -> Left = Left;
	    Node -> Right = Right;

	    return TRUE;	
        case 3:
        default:
	    return FALSE;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Free Mvar2CtPairstruct recursively.                                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:    Mvar2CtPairstruct node to free.                                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtFreePair(Mvar2CtPairStruct *Node)
{
    if (Node == NULL)
	return;
	
    Mvar2CtFreePair(Node -> Left);
    Mvar2CtFreePair(Node -> Right);
    
    IritFree(Node);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Check if two pairs of curves can have 2 contact trace.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1, Node2: pairs of curves to test.                                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if there can be no solution, FALSE otherwise.            *
*****************************************************************************/
static CagdBType Mvar2CtTrivialReject(Mvar2CtPairStruct *Node1, 
				      Mvar2CtPairStruct *Node2)
{
    if (Node1 -> Cparam -> TMax <= Node2 -> Cparam -> TMin)
	/* If rotation of pairs doesn't overlap, we can purge it. */ 
	return TRUE;
    
    if (Node1 -> Cparam -> TMin >= Node2 -> Cparam -> TMax)
	/* If rotation of pairs doesn't overlap, we can purge it. */ 
	return TRUE;

    if (Node1 -> ANode -> Min >= Node2 -> ANode -> Max && 
	                  Node1 -> BNode -> Id == Node2 -> BNode -> Id)
	/* We assume that moving curve of Node1 has smaller parameter when */ 
	/* Node1 and Node2 are pairs from same obstacle curve.             */
	return TRUE;

    if (Mvar2CtRejectbyCurvature(Node1 -> ANode, Node1 -> BNode) || 
	Mvar2CtRejectbyCurvature(Node2 -> ANode, Node2 -> BNode))
	/* Check curvature condition. */
	return TRUE;

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Check if two pairs of curves are connected.                               *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1, Node2: Pair of curves to test.                                    *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if connected, FALSE otherwise.                           *
*****************************************************************************/
static CagdBType Mvar2CtIsConnectedPair(Mvar2CtPairStruct *Node1, 
					Mvar2CtPairStruct *Node2)
{
    if (Node1 -> BNode -> Id != Node2 -> BNode -> Id)
	return FALSE;

    if (Node1 -> ANode -> Max < Node2 -> ANode -> Min && 
       !(Node1 -> ANode -> Min == 0 && Node2 -> ANode -> Max == 1))
	return FALSE;

    if (Node1 -> ANode -> Min > Node2 -> ANode -> Max && 
       !(Node1 -> ANode -> Max == 1 && Node2 -> ANode -> Min == 0))
	return FALSE;

    if (Node1 -> BNode -> Max < Node2 -> BNode -> Min && 
	!(Node1 -> BNode -> Min == 0 && Node2 -> BNode -> Max == 1))
	return FALSE;
    if (Node1 -> BNode -> Min > Node2 -> BNode -> Max && 
	!(Node1 -> BNode -> Max == 1 && Node2 -> BNode -> Min == 0))
	return FALSE;

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes rotational extreme points for 2 Contact trace.                  *
*   The algebraic conditions for rotational extreme points for 2 contact     *
*   trace are following:                                                     *
*                                                                            *
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_x                                       * 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_x,    * 
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_y                                       * 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_y,    * 
*    det (Rot(theta) * CrvA'(u1), CrvB'(v1) = 0,                             *
*    det (Rot(theta) * CrvA'(u2), CrvB'(v2) = 0,                             *
*    det (CrvB'(v1), CrvB'(v2) = 0.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1:  Curve pair (moving curve, obstacle curve) for 1st contact point  *
*   Node2:  Curve pair (moving curve, obstacle curve) for 2st contact point  * 
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhBs:  Bvh for obstacle curves.                                         *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: Linked list of the solutions holding parameter            *
* of rotational extreme points for 2 contact traces.                         *
*****************************************************************************/
static MvarPtStruct *Mvar2CtSolveRotExtreme(Mvar2CtPairStruct *Node1,
				            Mvar2CtPairStruct *Node2,
				            Mvar2CtBVHStruct *BvhA,
				            Mvar2CtBVHStruct **BvhBs,
				            CagdRType Subtol,
				            CagdRType Numerictol)
{
    CagdCrvStruct *Circle, *Subcircle, *CrvA1, *CrvA2, *CrvB1, *CrvB2, 
	          *TCrv, *TCrv2;
    CagdRType Min[5], Max[5];
    CagdPType Translate;
    MvarMVStruct *MTemp1, *MTemp2, *MTemp3, *MTemp4, *MTemp5, *MTemp6, 
	*MCrvA1, *MCrvA2, *MCrvB1, *MCrvB2, *MCircle, *MDCrvA1, *MDCrvA2, 
	*MDCrvB1, *MDCrvB2, *MTCrv1, *MTCrv2, *MCrvA1Scalar[2], 
	*MCrvA2Scalar[2], *MDCrvA1Scalar[2], *MDCrvA2Scalar[2],  
	*MCrvB1Scalar[2], *MCrvB2Scalar[2], *MDCrvB1Scalar[2], 
	*MDCrvB2Scalar[2], *MCircleScalar[3], **MScalar, *MVs[9];
    MvarConstraintType Constr[9];
    MvarPtStruct
        *MPts = NULL;
    int i, 
	MV_Num = 8;
    Mvar2CtBVNodeStruct *Nodes[4];
    Mvar2CtCParamStruct *Cparam;

    if (Node1 -> Cparam -> TMin >= Node2 -> Cparam -> TMin &&
	Node1 -> Cparam -> TMax <= Node2 -> Cparam -> TMax)
	Cparam = Node1 -> Cparam;
    else
	Cparam = Node2 -> Cparam;

    for (i = 0; i < 5; ++i) {
	Min[i] = IRIT_INFNTY;
	Max[i] = -IRIT_INFNTY;
    }

    Nodes[0] = Node1 -> ANode;
    Nodes[1] = Node1 -> BNode;
    Nodes[2] = Node2 -> ANode;
    Nodes[3] = Node2 -> BNode;

    Circle = BspCrvCreateUnitCircle();
    BspKnotAffineTrans2(Circle -> KnotVector, 
	                Circle -> Length + Circle -> Order, 0, 1);

    Mvar2CtReduceRotExtremeDomain(Nodes, Cparam, Min, Max, Subtol);

    if (Min[0] == IRIT_INFNTY) {
	CagdCrvFree(Circle);
	return NULL;
    }

    if (Node1 -> BNode -> Id == Node2 -> BNode -> Id) {
	TCrv = BspCrvNew(2, 2, CAGD_PT_BASE); 
	TCrv2 = BspCrvNew(2, 2, CAGD_PT_BASE);

	BspKnotUniformOpen(2, 2, TCrv -> KnotVector);
	BspKnotUniformOpen(2, 2, TCrv2 -> KnotVector);

	TCrv -> Points[1][0] = Min[0];
	TCrv -> Points[1][1] = Max[0];

	TCrv2 -> Points[1][0] = Min[2];
	TCrv2 -> Points[1][1] = Max[2];

	MTemp1 = MvarCrvToMV(TCrv);
	MTCrv1 = MvarPromoteMVToMV2(MTemp1, 5, 0);
	MvarMVFree(MTemp1);

	MTemp1 = MvarCrvToMV(TCrv2);
	MTCrv2 = MvarPromoteMVToMV2(MTemp1, 5, 2);
	MvarMVFree(MTemp1);

	MvarMVSetAllDomains(MTCrv1, Min, Max, TRUE);
	MvarMVSetAllDomains(MTCrv2, Min, Max, TRUE);

	MV_Num = 8;
    }
    else {
	TCrv = NULL;
	TCrv2 = NULL;
	MTCrv1 = NULL;
	MTCrv2 = NULL;
	MV_Num = 7;
    }

    /* Extract curve region. */
    Subcircle = CagdCrvRegionFromCrv(Circle, Min[4], Max[4]);

    CrvA1 = CagdCrvRegionFromCrv(BvhA -> Crv, Min[0], Max[0]);
    CrvB1 = CagdCrvRegionFromCrv(BvhBs[Node1 -> BNode -> Id] -> Crv,
				 Min[1], Max[1]);
    CrvA2 = CagdCrvRegionFromCrv(BvhA -> Crv, Min[2], Max[2]);
    CrvB2 = CagdCrvRegionFromCrv(BvhBs[Node2 -> BNode -> Id] -> Crv,
				 Min[3], Max[3]);
    
    Constr[0] = Constr[1] = Constr[2] = Constr[3] = Constr[4] =
                                                            MVAR_CNSTRNT_ZERO;
    Constr[5] = Constr[6] = MVAR_CNSTRNT_NEGATIVE;
    Constr[7] = MVAR_CNSTRNT_POSITIVE;
    
    MTemp1 = MvarCrvToMV(CrvA1);
    MCrvA1 = MvarPromoteMVToMV2(MTemp1, 5, 0);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCrvA1, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(CrvB1);
    MCrvB1 = MvarPromoteMVToMV2(MTemp1, 5, 1);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCrvB1, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(CrvA2);
    MCrvA2 = MvarPromoteMVToMV2(MTemp1, 5, 2);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCrvA2, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(CrvB2);
    MCrvB2 = MvarPromoteMVToMV2(MTemp1, 5, 3);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCrvB2, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(Subcircle);
    MCircle = MvarPromoteMVToMV2(MTemp1, 5, 4);
    MvarMVFree(MTemp1);
    
    MvarMVSetAllDomains(MCircle, Min, Max, TRUE);

    MDCrvA1 = MvarMVDerive(MCrvA1, 0);
    MDCrvB1 = MvarMVDerive(MCrvB1, 1);
    MDCrvA2 = MvarMVDerive(MCrvA2, 2);
    MDCrvB2 = MvarMVDerive(MCrvB2, 3);

    MScalar = MvarMVSplitScalar(MCrvA1);

    MCrvA1Scalar[0] = MScalar[1];
    MCrvA1Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MCrvA2);

    MCrvA2Scalar[0] = MScalar[1];
    MCrvA2Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MCrvB1);

    MCrvB1Scalar[0] = MScalar[1];
    MCrvB1Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MCrvB2);

    MCrvB2Scalar[0] = MScalar[1];
    MCrvB2Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MDCrvA1);

    MDCrvA1Scalar[0] = MScalar[1];
    MDCrvA1Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MDCrvA2);

    MDCrvA2Scalar[0] = MScalar[1];
    MDCrvA2Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MDCrvB1);

    MDCrvB1Scalar[0] = MScalar[1];
    MDCrvB1Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MDCrvB2);

    MDCrvB2Scalar[0] = MScalar[1];
    MDCrvB2Scalar[1] = MScalar[2];

    MScalar = MvarMVSplitScalar(MCircle);

    MCircleScalar[0] = MScalar[0];
    MCircleScalar[1] = MScalar[1];
    MCircleScalar[2] = MScalar[2];
    if (MScalar[3] != NULL)
        MvarMVFree(MScalar[3]);

    /* Tangency constraint. */

    MTemp1 = MvarMVMult(MDCrvA1Scalar[0], MCircleScalar[1]);
    MTemp2 = MvarMVMult(MDCrvA1Scalar[1], MCircleScalar[2]);

    MTemp3 = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvA1Scalar[0], MCircleScalar[2]);
    MTemp2 = MvarMVMult(MDCrvA1Scalar[1], MCircleScalar[1]);

    MTemp4 = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB1Scalar[1], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB1Scalar[0], MTemp4);

    MVs[0] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB1Scalar[0], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB1Scalar[1], MTemp4);

    MVs[5] = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp1 = MvarMVMult(MDCrvA2Scalar[0], MCircleScalar[1]);
    MTemp2 = MvarMVMult(MDCrvA2Scalar[1], MCircleScalar[2]);

    MTemp3 = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvA2Scalar[0], MCircleScalar[2]);
    MTemp2 = MvarMVMult(MDCrvA2Scalar[1], MCircleScalar[1]);

    MTemp4 = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB2Scalar[1], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB2Scalar[0], MTemp4);

    MVs[1] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB2Scalar[0], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB2Scalar[1], MTemp4);

    MVs[6] = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    /* Intersection Constraints. */

    MTemp1 = MvarMVSub(MCrvA1Scalar[0], MCrvA2Scalar[0]);
    MTemp2 = MvarMVSub(MCrvA1Scalar[1], MCrvA2Scalar[1]);

    MTemp3 = MvarMVMult(MTemp1, MCircleScalar[1]);
    MTemp4 = MvarMVMult(MTemp2, MCircleScalar[2]);

    MTemp5 = MvarMVSub(MTemp3, MTemp4);

    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp3 = MvarMVMult(MTemp1, MCircleScalar[2]);
    MTemp4 = MvarMVMult(MTemp2, MCircleScalar[1]);

    MTemp6 = MvarMVAdd(MTemp3, MTemp4);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp1 = MvarMVSub(MCrvB1Scalar[0], MCrvB2Scalar[0]);

    MTemp3 = MvarMVMult(MTemp1, MCircleScalar[0]);

    MTemp2 = MvarMVSub(MCrvB1Scalar[1], MCrvB2Scalar[1]);

    MTemp4 = MvarMVMult(MTemp2, MCircleScalar[0]);

    MVs[2]= MvarMVSub(MTemp3, MTemp5);
    MVs[3]= MvarMVSub(MTemp4, MTemp6);
   
    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);
    MvarMVFree(MTemp5);
    MvarMVFree(MTemp6);

    MVs[4]= MvarMVCrossProd2D(MDCrvB1Scalar[0], MDCrvB1Scalar[1],
                              MDCrvB2Scalar[0], MDCrvB2Scalar[1]);

    if (MV_Num == 8) {
	Translate[0] = - Subtol;
	Translate[1] = 0;
	Translate[2] = 0;

	MVs[7] = MvarMVSub(MTCrv2, MTCrv1);

	MvarMVTransform(MVs[7], Translate, 1.0);
    }

    MPts = MvarMVsZeros0D(MVs, Constr, MV_Num, Subtol, Numerictol);

   for (i = 0; i < MV_Num; ++i)
       MvarMVFree(MVs[i]);
    
   if (MV_Num == 8) {
       CagdCrvFree(TCrv);
       CagdCrvFree(TCrv2);
       MvarMVFree(MTCrv1);
       MvarMVFree(MTCrv2);
   }

    CagdCrvFree(Circle);
    CagdCrvFree(Subcircle);
    CagdCrvFree(CrvA1);
    CagdCrvFree(CrvA2);
    CagdCrvFree(CrvB1);
    CagdCrvFree(CrvB2);

    MvarMVFree(MCrvA1);
    MvarMVFree(MCrvA2);
    MvarMVFree(MCrvB1);
    MvarMVFree(MCrvB2);

    MvarMVFree(MCircle);

    MvarMVFree(MDCrvA1);
    MvarMVFree(MDCrvA2);
    MvarMVFree(MDCrvB1);
    MvarMVFree(MDCrvB2);

    for (i = 0; i < 2; ++i) {
	MvarMVFree(MCrvA1Scalar[i]);
	MvarMVFree(MCrvA2Scalar[i]);
	MvarMVFree(MCrvB1Scalar[i]);
	MvarMVFree(MCrvB2Scalar[i]);

	MvarMVFree(MDCrvA1Scalar[i]);
	MvarMVFree(MDCrvA2Scalar[i]);
	MvarMVFree(MDCrvB1Scalar[i]);
	MvarMVFree(MDCrvB2Scalar[i]);

	MvarMVFree(MCircleScalar[i]);
    }

    MvarMVFree(MCircleScalar[2]);
   
    return MPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Compute rotational extreme points for 2 contact traces. We traverse curve  *
* pair struct (Mvar2ctPairStruct) and once it reaches leaf node we switch to *
* equation solver.                                                           *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1:  Curve pair (moving curve, obstacle curve) for 1st contact point. *
*   Node2:  Curve pair (moving curve, obstacle curve) for 2st contact point. *
*   BvhA:     Bvh for moving curve.                                   	     *
*   BvhBs:    Bvh for obstacle curves.                                       *
*   MPts:     Rotational extreme points will be pushed here.                 *
*   Subtol:    Subdivision tolerance for the computation.                    *
*   Numerictol:   Numerical tolerance for the computation.                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtRotationExtreme(Mvar2CtPairStruct *Node1,
			           Mvar2CtPairStruct *Node2,
			           Mvar2CtBVHStruct *BvhA,
			           Mvar2CtBVHStruct **BvhBs,
			           MvarPtStruct **MPts,
			           CagdRType Subtol,
			           CagdRType Numerictol)
{
    int Dir;
    MvarPtStruct *TPt,
        *RotExtreme = NULL;
  
    if (Node1 == Node2) {	
	if (Node1 -> ANode -> Convexity * Node1 -> BNode -> Convexity != 0)
	    /* Single convex/convex convex/concave pair cannot have        */
	    /* rotational extreme.                                         */
	    return;	

	if (Node1 -> Left != NULL)
	    Mvar2CtRotationExtreme(Node1 -> Left, Node2, BvhA, 
	                                   BvhBs, MPts, Subtol, Numerictol);
	  
	if (Node1 -> Right != NULL)
	    Mvar2CtRotationExtreme(Node1 -> Right, Node2, BvhA, 
	    BvhBs, MPts, Subtol, Numerictol);	   

	return;
    }

    if (!Mvar2CtNormalOverlapBoth(Node1 -> BNode, Node2 -> BNode) || 
	!Mvar2CtNormalOverlapBoth(Node1 -> ANode, Node2 -> ANode))
	/* Check tangency condition for rotational extreme points. */
	return;

    if (Mvar2CtTrivialReject(Node1, Node2))
	return;  

    if (Mvar2CtIsConnectedPair(Node1, Node2)) {
	if (Node1 -> ANode -> Convexity * Node2 -> ANode -> Convexity > 0 &&
	    Node1 -> BNode -> Convexity * Node2 -> BNode -> Convexity > 0)
	    /* Single convex/convex convex/concave pair                     */
	    /* cannot have rotational extreme.                              */
	    return;

	/* If moving curve segments are connected convex or concave curve   */
	/* we don't have rotational extreme.                                */
	if (Node1 -> ANode == Node2 -> ANode && 
	    Node1 -> ANode -> Convexity != 0)
	    return;
	/* If obstacle curve segments are connected convex or concave curve */
	/* we don't have rotational extreme.                                */
	if (Node1 -> BNode == Node2 -> BNode && 
	    Node1 -> BNode -> Convexity != 0)
	    return;

	if (Mvar2CtIsConnectedNode(Node1 -> ANode, Node2 -> ANode) &&  
	    Node1 -> ANode -> Convexity * Node2 -> ANode -> Convexity > 0)
	    return;
	if (Mvar2CtIsConnectedNode(Node1 -> BNode, Node2 -> BNode) &&  
	    Node1 -> BNode -> Convexity * Node2 -> BNode -> Convexity > 0)
	    return;

	Dir = Mvar2CtPairSubdivDir2(Node1, Node2);

	switch (Dir) {
	    case 0: 
	        if (Node1 -> Left != NULL)
		    Mvar2CtRotationExtreme(Node1 -> Left, Node2, BvhA, 
		                           BvhBs, MPts, Subtol, Numerictol);
	        if (Node1 -> Right != NULL)
		    Mvar2CtRotationExtreme(Node1 -> Right, Node2, BvhA, 
		                           BvhBs, MPts, Subtol, Numerictol);
	        break;
	    case 1:	    
	        if (Node2 -> Left!=NULL)
		    Mvar2CtRotationExtreme(Node1, Node2 -> Left, BvhA, 
		                           BvhBs, MPts, Subtol, Numerictol);
	        if (Node2 -> Right!=NULL)
		    Mvar2CtRotationExtreme(Node1, Node2 -> Right, BvhA, 
		                           BvhBs, MPts, Subtol, Numerictol);
	        break;	

	    case 2:
	        /* We reached leaf node, switch to multivariate solver. */
	        RotExtreme = Mvar2CtSolveRotExtreme(Node1, Node2, BvhA, BvhBs, 
		                                    Subtol, Numerictol);
	   
	        if (RotExtreme != NULL) {
		    for (TPt = RotExtreme; TPt != NULL; TPt = TPt -> Pnext) {
		        AttrSetIntAttrib(&(TPt -> Attr), 
	                "First", Node1 -> BNode -> Id);
		        AttrSetIntAttrib(&(TPt -> Attr), 
	                "Second", Node2 -> BNode -> Id);
		}
		
		(*MPts) = 
		    (MvarPtStruct *) CagdListAppend(RotExtreme, (*MPts));
	    }	   
	        return;
	}

    }
    else {
	if (Mvar2CtAABBOverlap(&(Node1 -> AABB), &(Node2 -> AABB))) {
	    if (Node1 -> ANode == Node2 -> ANode && 
		Node1 -> ANode -> Convexity != 0)
		/* If moving curve segments are connected convex or        */ 
		/* concave curve we don't have rotational extreme.         */
		return;
	    if (Node1 -> BNode == Node2 -> BNode && 
		Node1 -> BNode -> Convexity != 0)
		/* If obstacle curve segments are connected convex or      */ 
		/* concave curve we don't have rotational extreme.         */
		return;

	    if (Mvar2CtIsConnectedNode(Node1 -> ANode, Node2 -> ANode) &&  
		Node1 -> ANode -> Convexity * Node2 -> ANode -> Convexity > 0)
		return;
	    if (Mvar2CtIsConnectedNode(Node1 -> BNode, Node2 -> BNode) && 
		Node1 -> BNode -> Convexity * Node2 -> BNode -> Convexity > 0)
		return;

	    Dir = Mvar2CtPairSubdivDir2(Node1, Node2);

	    switch (Dir) {
	        case 0: 		
		    if (Node1 -> Left != NULL)
		        Mvar2CtRotationExtreme(Node1 -> Left, Node2, BvhA, 
		        BvhBs, MPts, Subtol, Numerictol);
		    if (Node1 -> Right != NULL)
		        Mvar2CtRotationExtreme(Node1 -> Right, Node2, BvhA, 
		        BvhBs, MPts, Subtol, Numerictol);
		    break;
	        case 1:		 
		    if (Node2 -> Left != NULL)
		        Mvar2CtRotationExtreme(Node1, Node2 -> Left, BvhA, 
		        BvhBs, MPts, Subtol, Numerictol);
		    if (Node2 -> Right != NULL)
		        Mvar2CtRotationExtreme(Node1, Node2 -> Right, BvhA, 
		        BvhBs, MPts, Subtol, Numerictol);		
		    break;
	        case 2:
		    /* We reached leaf node, switch to multivariate solver. */
		    RotExtreme = Mvar2CtSolveRotExtreme(Node1, Node2, BvhA, 
		                                        BvhBs, Subtol, 
							Numerictol);
		    if (RotExtreme != NULL) {
		        for (TPt = RotExtreme;
			     TPt != NULL; 
			     TPt = TPt -> Pnext) {
			    AttrSetIntAttrib(&(TPt -> Attr), 
			                     "First", Node1 -> BNode -> Id);
			    AttrSetIntAttrib(&(TPt -> Attr), 
			                     "Second", Node2 -> BNode -> Id);
			}

			(*MPts) = 
			    (MvarPtStruct *) CagdListAppend(RotExtreme,
							    (*MPts));
		    }
		    break;
	    }
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute curvature contact points for a given curve pair.                 *
*   The algebraic conditions for 2contact are following:                     *
*                                                                            *
*    k(u)^2 - k(v)^2 = 0,                                                    * 
*    k(u)*k'(u) - k(v)*k'(v)*dv/du = 0,                                      * 
*    det (Rot(theta) * CrvA'(u), CrvB'(v)= 0.			             *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:   Curve pair (moving curve, obstacle curve).                       *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhB:   Bvh for obstacle curve.                                          *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtstruct *: Linked list of curvature contact points.                 * 
*****************************************************************************/
static MvarPtStruct *Mvar2CtSolveCurvatureContact(Mvar2CtPairStruct *Node,
					           Mvar2CtBVHStruct *BvhA,
					           Mvar2CtBVHStruct *BvhB,
					           CagdRType Subtol,
					           CagdRType Numerictol)
{
    MvarMVStruct *MCurvatureA, *MCurvatureB, *MDCrvA, *MDCrvB, *MCircle, 
	         *MTemp1, *MTemp2, *MTemp3, *MTemp4, *MTemp5, *MTemp6,
                 *MVs[4], **MVScalar, *MDCrvAScalar[3], *MDCrvBScalar[3], 
		 *MCircleScalar[3], *MCurvatureAScalar[2], 
		 *MCurvatureBScalar[2]; 
    MvarPtStruct *MPts;
    MvarConstraintType Constr[4];
    CagdCrvStruct *Circle, *Subcircle, *DCrvA, *DCrvB, *CurvatureA, 
	          *CurvatureB;
    CagdRType Min[3], Max[3];
    int i;

    Min[0] = Node -> ANode -> Min;
    Min[1] = Node -> BNode -> Min;
    Min[2] = Node -> Cparam -> TMin;

    Max[0] = Node -> ANode -> Max;
    Max[1] = Node -> BNode -> Max;
    Max[2] = Node -> Cparam -> TMax;

    Circle = BspCrvCreateUnitCircle();
    BspKnotAffineTrans2(Circle -> KnotVector, Circle -> Order + 
	                Circle -> Length, 0, 1);

    /* Extract curve regions. */
    Subcircle = CagdCrvRegionFromCrv(Circle, Min[2], Max[2]);
    DCrvA = CagdCrvRegionFromCrv(BvhA -> DCrv, Min[0], Max[0]);
    DCrvB = CagdCrvRegionFromCrv(BvhB -> DCrv, Min[1], Max[1]);
    CurvatureA = CagdCrvRegionFromCrv(BvhA -> Curvature, Min[0], Max[0]);
    CurvatureB = CagdCrvRegionFromCrv(BvhB -> Curvature, Min[1], Max[1]);   

    Constr[0] = Constr[1] = Constr[2] = MVAR_CNSTRNT_ZERO;
    Constr[3] = MVAR_CNSTRNT_NEGATIVE;

    MTemp1 = MvarCrvToMV(DCrvA);
    MDCrvA = MvarPromoteMVToMV2(MTemp1, 3, 0);
    MvarMVFree(MTemp1);

    MTemp1 = MvarCrvToMV(DCrvB);
    MDCrvB = MvarPromoteMVToMV2(MTemp1, 3, 1);
    MvarMVFree(MTemp1);

    MTemp1 = MvarCrvToMV(CurvatureA);
    MCurvatureA = MvarPromoteMVToMV2(MTemp1, 3, 0);
    MvarMVFree(MTemp1);

    MTemp1 = MvarCrvToMV(CurvatureB);
    MCurvatureB = MvarPromoteMVToMV2(MTemp1, 3, 1);
    MvarMVFree(MTemp1);

    MTemp1 = MvarCrvToMV(Subcircle);
    MCircle = MvarPromoteMVToMV2(MTemp1, 3, 2);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MDCrvA, Min, Max, TRUE);
    MvarMVSetAllDomains(MDCrvB, Min, Max, TRUE);
    MvarMVSetAllDomains(MCurvatureA, Min, Max, TRUE);
    MvarMVSetAllDomains(MCurvatureB, Min, Max, TRUE);
    MvarMVSetAllDomains(MCircle, Min, Max, TRUE);

    /* Split the MVs into scalar components. */
    MVScalar = MvarMVSplitScalar(MDCrvA);

    MDCrvAScalar[1] = MVScalar[1];
    MDCrvAScalar[2] = MVScalar[2];
    if (MVScalar[3] != NULL)
        MvarMVFree(MVScalar[3]);

    MVScalar = MvarMVSplitScalar(MDCrvB);

    MDCrvBScalar[1] = MVScalar[1];
    MDCrvBScalar[2] = MVScalar[2];

    MVScalar = MvarMVSplitScalar(MCircle);

    MCircleScalar[0] = MVScalar[0];
    MCircleScalar[1] = MVScalar[1];
    MCircleScalar[2] = MVScalar[2];
    if (MVScalar[3] != NULL)
        MvarMVFree(MVScalar[3]);

    MVScalar = MvarMVSplitScalar(MCurvatureA);

    MCurvatureAScalar[0] = MVScalar[0];
    MCurvatureAScalar[1] = MVScalar[1];

    MVScalar = MvarMVSplitScalar(MCurvatureB);

    MCurvatureBScalar[0] = MVScalar[0];
    MCurvatureBScalar[1] = MVScalar[1];

    MTemp1 = MvarMVMult(MCurvatureAScalar[1], MCurvatureBScalar[0]);
    MTemp2 = MvarMVMult(MCurvatureAScalar[0], MCurvatureBScalar[1]);

    /* Equal curvature constraint. */
    MVs[0] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvAScalar[1], MCircleScalar[1]);
    MTemp2 = MvarMVMult(MDCrvAScalar[2], MCircleScalar[2]);

    MTemp3 = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvAScalar[1], MCircleScalar[2]);
    MTemp2 = MvarMVMult(MDCrvAScalar[2], MCircleScalar[1]);

    MTemp4 = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);


    MTemp1 = MvarMVMult(MTemp3, MDCrvBScalar[2]);
    MTemp2 = MvarMVMult(MTemp4, MDCrvBScalar[1]);

    /* Tangency constraint */
    MVs[1] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MTemp3, MDCrvBScalar[1]);
    MTemp2 = MvarMVMult(MTemp4, MDCrvBScalar[2]);

    MVs[3] = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);  
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp1 = MvarMVDerive(MVs[0], 0);
    MTemp2 = MvarMVDerive(MVs[0], 1);
    MTemp3 = MvarMVDerive(MVs[1], 0);
    MTemp4 = MvarMVDerive(MVs[1], 1);

    MTemp5 = MvarMVMult(MTemp1, MTemp4);
    MTemp6 = MvarMVMult(MTemp2, MTemp3);

    /* Equal curvature derivative constraint. */
    MVs[2] = MvarMVSub(MTemp5, MTemp6);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);  
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);
    MvarMVFree(MTemp5);
    MvarMVFree(MTemp6);

    MPts = MvarMVsZeros0D(MVs, Constr, 4, Subtol, -IRIT_FABS(Numerictol));

    for (i = 0; i < 4; ++i)
	MvarMVFree(MVs[i]);

    MvarMVFree(MDCrvA);
    MvarMVFree(MDCrvB);

    MvarMVFree(MCurvatureA);
    MvarMVFree(MCurvatureB);

    MvarMVFree(MCircle);

    MvarMVFree(MDCrvAScalar[1]);
    MvarMVFree(MDCrvAScalar[2]);
    MvarMVFree(MDCrvBScalar[1]);
    MvarMVFree(MDCrvBScalar[2]);
    MvarMVFree(MCircleScalar[0]);
    MvarMVFree(MCircleScalar[1]);
    MvarMVFree(MCircleScalar[2]);
    MvarMVFree(MCurvatureAScalar[0]);
    MvarMVFree(MCurvatureAScalar[1]);
    MvarMVFree(MCurvatureBScalar[0]);
    MvarMVFree(MCurvatureBScalar[1]);

    CagdCrvFree(DCrvA);
    CagdCrvFree(DCrvB);
    CagdCrvFree(Circle);
    CagdCrvFree(Subcircle);
    CagdCrvFree(CurvatureA);
    CagdCrvFree(CurvatureB);

    return MPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute equal curvature trace for a given curve pair.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:   Curve pair (moving curve, obstacle curve).                       *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhB:   Bvh for obstacle curve.                                          *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylinestruct *: Linked list of equal curvature trace.              * 
*****************************************************************************/
static MvarPolylineStruct *Mvar2CtSolveCurvatureTrace(Mvar2CtPairStruct *Node,
					              Mvar2CtBVHStruct *BvhA,
					              Mvar2CtBVHStruct *BvhB,
					              CagdRType Subtol,
					              CagdRType Numerictol)
{
    MvarMVStruct *MCurvatureA, *MCurvatureB, 
	         *MTemp1, *MTemp2, *MVs[1], **MVScalar,  
		 *MCurvatureAScalar[2], *MCurvatureBScalar[2]; 
    MvarConstraintType Constr[1];
    MvarPolylineStruct *MPolys;
    CagdCrvStruct *CurvatureA, *CurvatureB;
    CagdRType Min[2], Max[2];
    int i;

    Min[0] = Node -> ANode -> Min;
    Min[1] = Node -> BNode -> Min;

    Max[0] = Node -> ANode -> Max;
    Max[1] = Node -> BNode -> Max;
    
    /* Extract curve regions. */

    CurvatureA = CagdCrvRegionFromCrv(BvhA -> Curvature, Min[0], Max[0]);
    CurvatureB = CagdCrvRegionFromCrv(BvhB -> Curvature, Min[1], Max[1]);   

    Constr[0] = MVAR_CNSTRNT_ZERO;

    MTemp1 = MvarCrvToMV(CurvatureA);
    MCurvatureA = MvarPromoteMVToMV2(MTemp1, 2, 0);
    MvarMVFree(MTemp1);

    MTemp1 = MvarCrvToMV(CurvatureB);
    MCurvatureB = MvarPromoteMVToMV2(MTemp1, 2, 1);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCurvatureA, Min, Max, TRUE);
    MvarMVSetAllDomains(MCurvatureB, Min, Max, TRUE);

    MVScalar = MvarMVSplitScalar(MCurvatureA);

    MCurvatureAScalar[0] = MVScalar[0];
    MCurvatureAScalar[1] = MVScalar[1];

    MVScalar = MvarMVSplitScalar(MCurvatureB);

    MCurvatureBScalar[0] = MVScalar[0];
    MCurvatureBScalar[1] = MVScalar[1];

    /* Equal curvature condition. */
    MTemp1 = MvarMVMult(MCurvatureAScalar[1], MCurvatureBScalar[0]);
    MTemp2 = MvarMVMult(MCurvatureAScalar[0], MCurvatureBScalar[1]);

    MVs[0] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MPolys = MvarMVsZeros1D(MVs, Constr, 1, Subtol, Subtol, Numerictol);

    for (i = 0; i < 1; ++i)
	MvarMVFree(MVs[i]);

    MvarMVFree(MCurvatureA);
    MvarMVFree(MCurvatureB);

    MvarMVFree(MCurvatureAScalar[0]);
    MvarMVFree(MCurvatureAScalar[1]);
    MvarMVFree(MCurvatureBScalar[0]);
    MvarMVFree(MCurvatureBScalar[1]);

    CagdCrvFree(CurvatureA);
    CagdCrvFree(CurvatureB);

    return MPolys;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute curvature contact points. Traversing pair struct, we only        * 
*   compute the curvature contact of convex/concave or concave/convex pair.  *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:   Curve pair (moving curve, obstacle curve).                       *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhB:   Bvh for obstacle curve.                                          *
*   CurvContacts: Linked list of curvature contact.                          *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtComputeCurvatureContact(Mvar2CtPairStruct *Node,
				           Mvar2CtBVHStruct *BvhA,
				           Mvar2CtBVHStruct *BvhB,
				           MvarPtStruct **CurvContacts,
				           CagdRType Subtol,
			                   CagdRType Numerictol)
{
    if (Node -> ANode -> Convexity * Node -> BNode -> Convexity) {
	if (Node -> ANode -> Convexity * Node->BNode -> Convexity < 0 
	    && Mvar2CtCurvatureOverlap(Node -> ANode, Node -> BNode, 
	                                              Node -> Cparam)) {
	    /* we solve curvature contact for a convex/concave or          */ 
	    /* concave/convex pair.                                        */
	    MvarPtStruct
	        *Pts = Mvar2CtSolveCurvatureContact(Node, BvhA, BvhB, Subtol,
						    Numerictol);	   
	    
	    if (Pts != NULL)
	        (*CurvContacts) = (MvarPtStruct *) CagdListAppend(Pts, 
	                                                     (*CurvContacts));
	}
    }
    else {
	if (Node -> Left != NULL)
	    Mvar2CtComputeCurvatureContact(Node -> Left, BvhA, BvhB, 
					   CurvContacts, Subtol, Numerictol);

	if (Node -> Right != NULL)
	    Mvar2CtComputeCurvatureContact(Node -> Right, BvhA, BvhB, 
					   CurvContacts, Subtol, Numerictol);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Set single contactness recursively.                                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Node: curve pair to set the single contactness.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtSetSingleRecursive(Mvar2CtPairStruct *Node)
{
    Node -> Bsingle = TRUE;

    if (Node -> Left != NULL)
	Mvar2CtSetSingleRecursive(Node -> Left);
    if (Node -> Right!=NULL)
	Mvar2CtSetSingleRecursive(Node -> Right);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Check if a point is located in a domain of Node1 and Node2.               *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1: curve pair for a first contact point.                             *
*   Node2: curve pair for a second contact point.                            *
*   MPt: Point to check.                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if MPt is in the domain, FALSE otherwise.                *
*****************************************************************************/
static CagdBType Mvar2CtInPairDomain(Mvar2CtPairStruct *Node1, 
			             Mvar2CtPairStruct *Node2, 
			             MvarPtStruct *MPt)
{ 
    if (Node1 -> ANode -> Min > MPt -> Pt[0] ||
	Node1 -> ANode -> Max < MPt -> Pt[0])
	return FALSE;
    if (Node1 -> BNode -> Min > MPt -> Pt[1] ||
	Node1 -> BNode -> Max < MPt -> Pt[1])
	return FALSE;
    if (Node1 -> Cparam -> TMin > MPt -> Pt[4] ||
	Node1 -> Cparam -> TMax < MPt -> Pt[4])
	return FALSE;

    if (Node2 -> ANode -> Min > MPt -> Pt[2] ||
	Node2 -> ANode -> Max < MPt -> Pt[2])
	return FALSE;
    if (Node2 -> BNode -> Min > MPt -> Pt[3] ||
	Node2 -> BNode -> Max < MPt -> Pt[3])
	return FALSE;
    if (Node2 -> Cparam -> TMin > MPt -> Pt[4] ||
	Node2 -> Cparam -> TMax < MPt -> Pt[4])
	return FALSE;

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Set the single contact region of curve pairs. To have 2 contact          * 
*  configuration, curve pair should have same curvature at two different     * 
*  locations. We check if curve pair satisfy this condition and set the      *
*  single contact region.                                                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:   Curve pair (moving curve, obstacle curve).                       *
*   CurvContact: Linked list of curvature contact.                           *
*   CurvTrace: Linked list of equal curvature trace.                         *
*   Tol:   Numerical tolerance of computation.		                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtCheckSingleContact(Mvar2CtPairStruct *Node,
			              MvarPtStruct *CurvContact,
			              MvarPolylineStruct *CurvTrace,
			              CagdRType Tol)
{    
    CagdRType Min[3], Max[3];
    MvarPolylineStruct *TPoly;
    MvarPtStruct *TPt, 
	*CPt = NULL;
    int CurvContactCount = 0, 
	CurvTraceCount = 0;

    Min[0] = Node -> ANode -> Min;
    Min[1] = Node -> BNode -> Min;
    Min[2] = Node -> Cparam -> TMin;

    Max[0] = Node -> ANode -> Max;
    Max[1] = Node -> BNode -> Max;
    Max[2] = Node -> Cparam -> TMax;

    for (TPt = CurvContact; TPt != NULL; TPt = TPt -> Pnext) {
	/* Count how many curvature contact points in the domain. */
	if (Mvar2CtInDomain(Min, Max, TPt)) {
	    CPt = TPt;
	    CurvContactCount++;
	}

	if (CurvContactCount > 1)
	    break;
    }
  
    for (TPoly = CurvTrace; TPoly != NULL; TPoly = TPoly -> Pnext) {
	/* Count how many equal curvature trace passing this domain. */
	if (Mvar2CtIsPassing(Min, Max, TPoly))
	    CurvTraceCount++;

	if (CurvTraceCount > 1)
	    break;
    }
    
    if (CurvContactCount == 0 && CurvTraceCount < 2)
        /* If this curve pair doesn't satisfy necessary condition          */
	/* for 2contact.                                                   */
        Mvar2CtSetSingleRecursive(Node);

    else {
	if (Node -> Left == NULL && Node -> Right == NULL) {
	    if (CurvTraceCount == 1) {
		Node -> CurvContacts = MvarPtCopy(CPt);
	    }
	    else {
	        if (!Mvar2CtAddChildPair(Node, Tol))
		    Node -> Bsingle = TRUE;
	    }
	}
	
	if (Node -> Left != NULL)
	    Mvar2CtCheckSingleContact(Node -> Left, CurvContact, 
	                              CurvTrace, Tol);
	if (Node -> Right != NULL)
	    Mvar2CtCheckSingleContact(Node -> Right, CurvContact, 
	                              CurvTrace, Tol);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Set the single contact region of curve pairs.                            *
*                                                                            *
* PARAMETERS:                                                                *
*   Node:   Curve pair (moving curve, obstacle curve).                       *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhB:   Bvh for obstacle curve.                                          *
*   CurvContact: Linked list of curvature contact.                           *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtSetSingleContact(Mvar2CtPairStruct *Node,
			            Mvar2CtBVHStruct *BvhA,
			            Mvar2CtBVHStruct *BvhB,
			            MvarPtStruct *CurvContact,
			            CagdRType Subtol,
			            CagdRType Numerictol)
{
    MvarPolylineStruct
        *MPoly = NULL;
    MvarPtStruct *TPt;
    CagdRType Min[3], Max[3];

    if (Node -> ANode -> Convexity && Node -> BNode -> Convexity) {
	if (Node -> ANode -> Convexity * Node -> BNode -> Convexity > 0)
	    /* if both curve are convex there cannot exist 2 contact. */
	    Mvar2CtSetSingleRecursive(Node);	    

	if (Node -> ANode -> Convexity * Node -> BNode -> Convexity < 0) {
           /* If there exist a 2 contact configuration btw                 */
           /* convex/concave pair curve pair should share same curvatures  */
           /* at two different location.                                   */

	   if (Mvar2CtCurvatureOverlap(Node -> ANode, Node -> BNode, 
	                               Node -> Cparam)) {
                /* Find solutions where the curve pair shares              */ 
		/* the same curvature. */ 	     
 		MPoly = Mvar2CtSolveCurvatureTrace(Node, BvhA, BvhB, 
		                                   Subtol, Numerictol);
		Min[0] = Node -> ANode -> Min;
		Max[0] = Node -> ANode -> Max;
		Min[1] = Node -> BNode -> Min;
		Max[1] = Node -> BNode -> Max;
		Min[2] = Node -> Cparam -> TMin;
		Max[2] = Node -> Cparam -> TMax;

		for (TPt = CurvContact; TPt != NULL; TPt = TPt -> Pnext)
		{
		    if(Mvar2CtInDomain(Min, Max, TPt))
		    Mvar2CtInsertCurvPoint(MPoly, TPt);
		}
		/* Cut the curvature trace at curvature contact point. */ 
		MPoly = Mvar2CtDivideAtCurvPointAux(MPoly);

		/* Check for single contactness. */
		Mvar2CtCheckSingleContact(Node, CurvContact, MPoly, Subtol);
		MvarPolylineFreeList(MPoly);
	    }
	    else
		Mvar2CtSetSingleRecursive(Node);	    
	}
    }
    else {
	if (Node -> Left != NULL)
	    Mvar2CtSetSingleContact(Node -> Left, BvhA, BvhB, 
	                            CurvContact, Subtol, Numerictol);

	if (Node -> Right != NULL)
	    Mvar2CtSetSingleContact(Node -> Right, BvhA, BvhB, 
	                            CurvContact, Subtol, Numerictol);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute boundary solution for multivariate having 1 dimensional solution *
* For each boundary equation, we reduce the domain of the equation by using  *
* bounding volume hierarchy.                                                 *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:	 Array of multivariates.				     *
*   Constr:   Equality or inequality constraints.		             *
*   MV_Num:	 Number of constraints (may be updated).		     *
*   Nodes:    Bounding volumes for curves.                                   *
*   Cparam:      data structure for rotation.		                     *
*   Subtol:      Subdivision tolerance of Computation.                       *
*   Numerictol:  Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:  Boundary solution for multivariates.                    *
*****************************************************************************/
static MvarPtStruct *Mvar2CtBndrySolution(MvarMVStruct **MVs, 
				          MvarConstraintType *Constr, 
				          int MV_Num, 
				          Mvar2CtBVNodeStruct *Nodes[4], 
				          Mvar2CtCParamStruct *Cparam, 
				          CagdRType Subtol, 
				          CagdRType Numerictol)
{
    MvarMVStruct *BoundaryMVs[10], **SubMVs;
    CagdRType MinDmn[5], MaxDmn[5], Min[5], Max[5], SubMin[4], SubMax[4],
             *TempBound;
    /* Twice the same for TempBound = MinDmn, MaxDmn. */
    MvarPtStruct *Aux,
        *BoundaryIntersPts = NULL, 
        *StartEndPts = NULL;
    int l, i, j, k, s;

    MvarMVDomain(MVs[0], MinDmn, MaxDmn, -1);

    TempBound = MinDmn;

    for (l = 0; l < 2; l++) {	
	for (i = 0; i < 5; i++) {
	    for ( k = 0; k < 5; ++k) {
		Min[k] = IRIT_INFNTY;
		Max[k] = -IRIT_INFNTY;
	    }
	    /* Reduce the domain using BVH. */
	    Mvar2CtReduce2CtDomain(Nodes, Cparam, Min, Max, l, i, Subtol);

	    /* If there is valid domain stop processing. */
	    if (Min[0] == IRIT_INFNTY)
		continue;
 
	    for (k = 0, s = 0; k < 5; ++k) {
		if (k == i)
		    continue;
		SubMin[s] = Min[k];
		SubMax[s] = Max[k];		
		s++;
	    }

	    for (j = 0; j < MV_Num; j++)
		BoundaryMVs[j] = MvarMVFromMV(MVs[j], TempBound[i], i);	  
	   
	    /* Extract multivariate correspond to reduced domain. */
	    SubMVs = Mvar2CtExtractMVRegion(BoundaryMVs, MV_Num,
					    SubMin, SubMax);	
	   	
	    BoundaryIntersPts = MvarMVsZeros0D(SubMVs, Constr,
					       MV_Num, Subtol,
					       IRIT_FABS(Numerictol));

	    for (j = 0; j < MV_Num; j++)
		MvarMVFree(BoundaryMVs[j]);

	    for (j = 0; j < MV_Num; ++j)
		MvarMVFree(SubMVs[j]);

	    IritFree(SubMVs);

	    for (Aux = BoundaryIntersPts; Aux != NULL; Aux = Aux -> Pnext) {
		MvarPtStruct 
		    *NewPt = MvarPtNew(5);

		for (k = 0; k < i; k++)
		    NewPt -> Pt[k] = Aux -> Pt[k];	    
		    NewPt -> Pt[i] = TempBound[i];
		for (k = i + 1; k < 5; k++)
		    NewPt -> Pt[k] = Aux -> Pt[k - 1];	    

		IRIT_LIST_PUSH(NewPt, StartEndPts);
	    } 
	    MvarPtFreeList(BoundaryIntersPts);
	} 

	TempBound = MaxDmn;
    }

    return StartEndPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the 2 Contact trace between the moving curve and obstacle       *
*   curves.                                                                  *
*   The algebraic conditions for 2contact are following:                     *
*                                                                            *
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_x                                       * 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_x     * 
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_y                                       * 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_y     * 
*    det (Rot(theta) * CrvA'(u1), CrvB'(v1) = 0,                             *
*    det (Rot(theta) * CrvA'(u2), CrvB'(v2) = 0.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1:  Curve pair (moving curve, obstacle curve) for 1st contact point. *
*   Node2:  Curve pair (moving curve, obstacle curve) for 2st contact point. * 
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhBs:  Bvh for obstacle curves.                                         *
*   NoLoop: TRUE If there is no loop in the domain FALSE otherwise.          *
*   Step:         Step size to use in the numeric tracing.		     *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: Linked list of the solutions holding parameter     *
* of 2 contact traces and rotation in radian (u1, v1, u2, v2, theta).        *
*****************************************************************************/
static MvarPolylineStruct *Mvar2CtSolve2CtTrace(Mvar2CtPairStruct *Node1,
				                Mvar2CtPairStruct *Node2,
				                Mvar2CtBVHStruct *BvhA,
				                Mvar2CtBVHStruct **BvhBs,
				                CagdBType NoLoop,
						CagdRType Step, 
				                CagdRType Subtol, 
				                CagdRType Numerictol)
{
    CagdCrvStruct *Circle, *Subcircle, *CrvA1, *CrvB1, *CrvA2, *CrvB2, 
	          *TCrv, *TCrv2;
    CagdRType Min[5], Max[5], 
        Translate = -Subtol*0.1;
    Mvar2CtCParamStruct *Cparam;
    /* Get BVH for obstacle curve. */
    Mvar2CtBVHStruct
        *BvhB1 = BvhBs[Node1 -> BNode -> Id], 
        *BvhB2 = BvhBs[Node2 -> BNode -> Id];
    MvarMVStruct *MTemp1, *MTemp2, *MTemp3, *MTemp4, *MTemp5, *MTemp6, 
	         *MCrvA1, *MCrvA2, *MCrvB1, *MCrvB2, *MCircle, *MDCrvA1, 
		 *MDCrvA2, *MDCrvB1, *MDCrvB2, *MTCrv1, *MTCrv2,
		 *MCrvA1Scalar[3], *MCrvA2Scalar[3], *MDCrvA1Scalar[3], 
		 *MDCrvA2Scalar[3],  *MCrvB1Scalar[3], *MCrvB2Scalar[3], 
		 *MDCrvB1Scalar[3], *MDCrvB2Scalar[3], *MCircleScalar[3], 
		 **MScalar, *MVs[9];
    MvarConstraintType Constr[9];   
    Mvar2CtBVNodeStruct *Nodes[4];
    MvarPolylineStruct
        *MPoly = NULL;
    MvarPtStruct
        *Bndry = NULL;
    /* If obstacle curve are same curve we have one more constraint. */ 
    int i,
        MVNum = Node1 -> BNode -> Id == Node2 -> BNode -> Id ? 7 : 6;
   
    if (Node1 -> Cparam -> TMin >= Node2 -> Cparam -> TMin &&
	Node1 -> Cparam -> TMax <= Node2 -> Cparam -> TMax)
	Cparam = Node1 -> Cparam;
    else
	Cparam = Node2 -> Cparam;

    Min[0] = Node1 -> ANode -> Min;
    Min[1] = Node1 -> BNode -> Min;
    Min[2] = Node2 -> ANode -> Min;
    Min[3] = Node2 -> BNode -> Min;
    Min[4] = Cparam -> TMin;

    Max[0] = Node1 -> ANode -> Max;
    Max[1] = Node1 -> BNode -> Max;
    Max[2] = Node2 -> ANode -> Max;
    Max[3] = Node2 -> BNode -> Max;
    Max[4] = Cparam -> TMax;

    Circle = BspCrvCreateUnitCircle();
    BspKnotAffineTrans2(Circle -> KnotVector, 
	                Circle -> Order + Circle -> Length, 0, 1);

    /* Extract curve region for given domain. */
    CrvA1 = CagdCrvRegionFromCrv(BvhA -> Crv, Min[0], Max[0]);
    CrvB1 = CagdCrvRegionFromCrv(BvhB1 -> Crv, Min[1], Max[1]);
    CrvA2 = CagdCrvRegionFromCrv(BvhA -> Crv, Min[2], Max[2]);
    CrvB2 = CagdCrvRegionFromCrv(BvhB2 -> Crv, Min[3], Max[3]);
    Subcircle = CagdCrvRegionFromCrv(Circle, Min[4], Max[4]); 

    Constr[0] = Constr[1] = Constr[2] = Constr[3] = MVAR_CNSTRNT_ZERO;
    Constr[4] = Constr[5] = MVAR_CNSTRNT_NEGATIVE;
    Constr[6] = MVAR_CNSTRNT_POSITIVE;  

    if (Node1 -> BNode -> Id == Node2 -> BNode -> Id) { 
        /* If two contact points are located in same obstacle curve, we add */
        /* one more constraint so that we can exclude trivial solution.     */
	TCrv = BspCrvNew(2, 2, CAGD_PT_BASE); 
	TCrv2 = BspCrvNew(2, 2, CAGD_PT_BASE);

	BspKnotUniformOpen(2, 2, TCrv -> KnotVector);
	BspKnotUniformOpen(2, 2, TCrv2 -> KnotVector);

	TCrv -> Points[1][0] = Min[0];
	TCrv -> Points[1][1] = Max[0];

	TCrv2 -> Points[1][0] = Min[2];
	TCrv2 -> Points[1][1] = Max[2];

	MTemp1 = MvarCrvToMV(TCrv);
	MTCrv1 = MvarPromoteMVToMV2(MTemp1, 5, 0);
	MvarMVFree(MTemp1);

	MTemp1 = MvarCrvToMV(TCrv2);
	MTCrv2 = MvarPromoteMVToMV2(MTemp1, 5, 2);
	MvarMVFree(MTemp1);

	MvarMVSetAllDomains(MTCrv1, Min, Max, TRUE);
	MvarMVSetAllDomains(MTCrv2, Min, Max, TRUE);

	MVs[6] = MvarMVSub(MTCrv2, MTCrv1);

	MvarMVTransform(MVs[6], &Translate, 1.0);
    }
    else {
	TCrv = NULL;
	TCrv2 = NULL;
	MTCrv1 = NULL;
	MTCrv2 = NULL;
    }
      
    MTemp1 = MvarCrvToMV(CrvA1);
    MCrvA1 = MvarPromoteMVToMV2(MTemp1, 5, 0);
    MvarMVFree(MTemp1);
    
    MvarMVSetAllDomains(MCrvA1, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(CrvB1);
    MCrvB1 = MvarPromoteMVToMV2(MTemp1, 5, 1);
    MvarMVFree(MTemp1);
    
    MvarMVSetAllDomains(MCrvB1, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(CrvA2);
    MCrvA2 = MvarPromoteMVToMV2(MTemp1, 5, 2);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCrvA2, Min, Max, TRUE);

    MTemp1 = MvarCrvToMV(CrvB2);
    MCrvB2 = MvarPromoteMVToMV2(MTemp1, 5, 3);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCrvB2, Min, Max, TRUE);   

    MTemp1 = MvarCrvToMV(Subcircle);
    MCircle = MvarPromoteMVToMV2(MTemp1, 5, 4);
    MvarMVFree(MTemp1);
    
    MvarMVSetAllDomains(MCircle, Min, Max, TRUE);

    MDCrvA1 = MvarMVDerive(MCrvA1, 0);
    MDCrvB1 = MvarMVDerive(MCrvB1, 1);
    MDCrvA2 = MvarMVDerive(MCrvA2, 2);
    MDCrvB2 = MvarMVDerive(MCrvB2, 3);

    MScalar = MvarMVSplitScalar(MCrvA1);

    MCrvA1Scalar[1] = MScalar[1];
    MCrvA1Scalar[2] = MScalar[2];
    
    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MCrvA2);

    MCrvA2Scalar[1] = MScalar[1];
    MCrvA2Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MCrvB1);

    MCrvB1Scalar[1] = MScalar[1];
    MCrvB1Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MCrvB2);

    MCrvB2Scalar[1] = MScalar[1];
    MCrvB2Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MDCrvA1);

    MDCrvA1Scalar[1] = MScalar[1];
    MDCrvA1Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MDCrvA2);

    MDCrvA2Scalar[1] = MScalar[1];
    MDCrvA2Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MDCrvB1);

    MDCrvB1Scalar[1] = MScalar[1];
    MDCrvB1Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MDCrvB2);

    MDCrvB2Scalar[1] = MScalar[1];
    MDCrvB2Scalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    MScalar = MvarMVSplitScalar(MCircle);

    MCircleScalar[0] = MScalar[0];
    MCircleScalar[1] = MScalar[1];
    MCircleScalar[2] = MScalar[2];

    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    /* Tangency constraint. */

    MTemp1 = MvarMVMult(MDCrvA1Scalar[1], MCircleScalar[1]);
    MTemp2 = MvarMVMult(MDCrvA1Scalar[2], MCircleScalar[2]);

    MTemp3 = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvA1Scalar[1], MCircleScalar[2]);
    MTemp2 = MvarMVMult(MDCrvA1Scalar[2], MCircleScalar[1]);

    MTemp4 = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB1Scalar[2], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB1Scalar[1], MTemp4);

    MVs[0] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB1Scalar[1], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB1Scalar[2], MTemp4);

    MVs[4] = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp1 = MvarMVMult(MDCrvA2Scalar[1], MCircleScalar[1]);
    MTemp2 = MvarMVMult(MDCrvA2Scalar[2], MCircleScalar[2]);

    MTemp3 = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvA2Scalar[1], MCircleScalar[2]);
    MTemp2 = MvarMVMult(MDCrvA2Scalar[2], MCircleScalar[1]);

    MTemp4 = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB2Scalar[2], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB2Scalar[1], MTemp4);

    MVs[1] = MvarMVSub(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);

    MTemp1 = MvarMVMult(MDCrvB2Scalar[1], MTemp3);
    MTemp2 = MvarMVMult(MDCrvB2Scalar[2], MTemp4);

    MVs[5] = MvarMVAdd(MTemp1, MTemp2);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    /* Intersection Constraints. */

    MTemp1 = MvarMVSub(MCrvA1Scalar[1], MCrvA2Scalar[1]);
    MTemp2 = MvarMVSub(MCrvA1Scalar[2], MCrvA2Scalar[2]);

    MTemp3 = MvarMVMult(MTemp1, MCircleScalar[1]);
    MTemp4 = MvarMVMult(MTemp2, MCircleScalar[2]);

    MTemp5 = MvarMVSub(MTemp3, MTemp4);

    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp3 = MvarMVMult(MTemp1, MCircleScalar[2]);
    MTemp4 = MvarMVMult(MTemp2, MCircleScalar[1]);

    MTemp6 = MvarMVAdd(MTemp3, MTemp4);

    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);

    MTemp1 = MvarMVSub(MCrvB1Scalar[1], MCrvB2Scalar[1]);

    MTemp3 = MvarMVMult(MTemp1, MCircleScalar[0]);

    MTemp2 = MvarMVSub(MCrvB1Scalar[2], MCrvB2Scalar[2]);

    MTemp4 = MvarMVMult(MTemp2, MCircleScalar[0]);

    MVs[2]= MvarMVSub(MTemp3, MTemp5);
    MVs[3]= MvarMVSub(MTemp4, MTemp6);
   
    MvarMVFree(MTemp1);
    MvarMVFree(MTemp2);
    MvarMVFree(MTemp3);
    MvarMVFree(MTemp4);
    MvarMVFree(MTemp5);
    MvarMVFree(MTemp6);  
	
    Nodes[0] = Node1 -> ANode;
    Nodes[1] = Node1 -> BNode;
    Nodes[2] = Node2 -> ANode;
    Nodes[3] = Node2 -> BNode;
    
    if (NoLoop) {    
	/* If there is no loop we compute boundary solution first. */
	Bndry = Mvar2CtBndrySolution(MVs, Constr, MVNum, 
	                             Nodes, Cparam, Subtol, Numerictol);

	if (CagdListLength(Bndry) == 0)
	    /* If there is no loop and there exist no boundary solution,    */ 
	    /*  we don't have 2 contact trace.                              */
	    MPoly = NULL;
	else if (Node1 == Node2 && Node1 -> CurvContacts != NULL && 
	        CagdListLength(Bndry) == 1) {
	    /* If this domain contains curvature contact and                */ 
	    /* there is one boundary point we trace the solution            */
	    /* between curvature contact and boundary solution.             */
	    MvarPtStruct *TPt = MvarPtNew(5);
	    TPt -> Pt[0] = Node1 -> CurvContacts -> Pt[0];
	    TPt -> Pt[1] = Node1 -> CurvContacts -> Pt[1];
	    TPt -> Pt[2] = Node1 -> CurvContacts -> Pt[0];
	    TPt -> Pt[3] = Node1 -> CurvContacts -> Pt[1];
	    TPt -> Pt[4] = Node1 -> CurvContacts -> Pt[2];

	    Bndry -> Pnext = TPt;	    
	    MPoly = MvarMVsZeros1DOneTrace(MVs, Constr, MVNum, Bndry, 
		                           Step, Subtol, Numerictol);
	}
	else if (CagdListLength(Bndry)==1)
	    MPoly = NULL;

	else if (CagdListLength(Bndry) == 2)
            /* In this case we have only single component in the domain,   */ 
	    /* thus we directly trace the solution.                        */
	    MPoly = MvarMVsZeros1DOneTrace(MVs, Constr, MVNum, Bndry, 
	                                   Step, Subtol, Numerictol);
	else
	    MPoly = MvarMVsZeros1D(MVs, Constr, MVNum, 
	                           Step, Subtol, Numerictol);
    }

    else
	MPoly = MvarMVsZeros1D(MVs, Constr, MVNum, 
	                       Step, Subtol, Numerictol);

    MvarPtFreeList(Bndry);

    if (MVNum == 7) {
	MvarMVFree(MTCrv1);
	MvarMVFree(MTCrv2);

	CagdCrvFree(TCrv);
	CagdCrvFree(TCrv2);
    }

    for (i = 0; i < MVNum; ++i)
	MvarMVFree(MVs[i]);
    
    CagdCrvFree(Circle);
    CagdCrvFree(Subcircle);
    CagdCrvFree(CrvA1);
    CagdCrvFree(CrvA2);
    CagdCrvFree(CrvB1);
    CagdCrvFree(CrvB2);

    MvarMVFree(MCrvA1);
    MvarMVFree(MCrvA2);
    MvarMVFree(MCrvB1);
    MvarMVFree(MCrvB2);

    MvarMVFree(MCircle);

    MvarMVFree(MDCrvA1);
    MvarMVFree(MDCrvA2);
    MvarMVFree(MDCrvB1);
    MvarMVFree(MDCrvB2);

    for (i = 1; i < 3; ++i) {
	MvarMVFree(MCrvA1Scalar[i]);
	MvarMVFree(MCrvA2Scalar[i]);
	MvarMVFree(MCrvB1Scalar[i]);
	MvarMVFree(MCrvB2Scalar[i]);

	MvarMVFree(MDCrvA1Scalar[i]);
	MvarMVFree(MDCrvA2Scalar[i]);
	MvarMVFree(MDCrvB1Scalar[i]);
	MvarMVFree(MDCrvB2Scalar[i]);

	MvarMVFree(MCircleScalar[i]);
    }

    MvarMVFree(MCircleScalar[0]);

    return MPoly;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the 3 Contact points between the moving curve and obstacle      *
*   curves for a given domain.                                               *
*   The algebraic conditions for 2contact are following.                     *
*                                                                            *
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_x                                       * 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_x     * 
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_y                                       * 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_y     * 
*  [CrvB(v2) -Rot(theta) * CrvA(u2)]_x                                       * 
*                                 = [CrvB(v3) - Rot(theta) * CrvA(u3)]_x     * 
*  [CrvB(v2) -Rot(theta) * CrvA(u2)]_y                                       * 
*                                 = [CrvB(v3) - Rot(theta) * CrvA(u3)]_y     * 
*    det (Rot(theta) * CrvA'(u1), CrvB'(v1) = 0,                             *
*    det (Rot(theta) * CrvA'(u2), CrvB'(v2) = 0,                             *
*    det (Rot(theta) * CrvA'(u3), CrvB'(v3) = 0.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   BvhA:   Bvh for moving curve.                                     	     *
*   BvhBs:  Bvh for obstacle curves.                                         *
*   Ids: Ids of obstacle curves.                                   	     *
*   Min: Minimum domain value of each parameter.                             *
*   Max: Maximum domain value of each parameter.                             *
*   Subtol:    Subdivision tolerance of Computation.                         *
*   Numerictol:   Numerical tolerance of computation.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: Linked list of the solutions holding parameter value     *
* of 3 contact points and rotation in radian (u1, v1, u2, v2, u3,v3, theta). *
*****************************************************************************/
static MvarPtStruct *Mvar2CtSolve3CtPoint(Mvar2CtBVHStruct *BvhA,
				          Mvar2CtBVHStruct **BvhBs,
				          int Ids[3],
				          CagdRType Min[7],
		                          CagdRType Max[7],
				          CagdRType Subtol, 
				          CagdRType Numerictol)
{
    CagdCrvStruct *Circle, *Subcircle, *CrvA[3], *CrvB[3];
    int i;
    MvarMVStruct *MTemp1, *MTemp2, *MTemp3, *MTemp4, *MTemp5, *MTemp6, 
	*MCrvA[3], *MCrvB[3], *MCircle, *MDCrvA[3], *MDCrvB[3],
	*MCrvAScalar[3][3],*MCrvBScalar[3][3], *MDCrvAScalar[3][3], 
	*MDCrvBScalar[3][3],  *MCircleScalar[3], **MScalar, *MVs[10];
    MvarPtStruct *MPts;
    MvarConstraintType Constr[10];

    Constr[0] = Constr[1] = Constr[2] = Constr[3] = Constr[4] 
                                  = Constr[5] = Constr[6] = MVAR_CNSTRNT_ZERO;
    Constr[7] = Constr[8] = Constr[9] = MVAR_CNSTRNT_NEGATIVE;

    Circle = BspCrvCreateUnitCircle();
    BspKnotAffineTrans2(Circle -> KnotVector, Circle -> Order + 
	                Circle -> Length, 0, 1);

    for (i = 0; i < 3; ++i) {
	/* Extract curve regions correspond to a given domain. */ 
	CrvA[i] = CagdCrvRegionFromCrv(BvhA -> Crv, Min[i * 2], Max[i * 2]);
	CrvB[i] = CagdCrvRegionFromCrv(BvhBs[Ids[i]] -> Crv, 
	                               Min[i * 2 + 1], Max[i * 2 + 1]);
    }

    Subcircle = CagdCrvRegionFromCrv(Circle, Min[6], Max[6]);

    for (i = 0; i < 3; ++i) {
	MTemp1 = MvarCrvToMV(CrvA[i]);
	MCrvA[i] = MvarPromoteMVToMV2(MTemp1, 7, i * 2);
	MvarMVFree(MTemp1);

	MvarMVSetAllDomains(MCrvA[i], Min, Max, TRUE);

	MTemp1 = MvarCrvToMV(CrvB[i]);
	MCrvB[i] = MvarPromoteMVToMV2(MTemp1, 7, i * 2 + 1);
	MvarMVFree(MTemp1);

	MvarMVSetAllDomains(MCrvB[i], Min, Max, TRUE);
    }

    MTemp1 = MvarCrvToMV(Subcircle);
    MCircle = MvarPromoteMVToMV2(MTemp1, 7, 6);
    MvarMVFree(MTemp1);

    MvarMVSetAllDomains(MCircle, Min, Max, TRUE);

    for (i = 0; i < 3; ++i) {
	/* Split the MVs into scalar functions. */

	MDCrvA[i] = MvarMVDerive(MCrvA[i], i * 2);
	MDCrvB[i] = MvarMVDerive(MCrvB[i], i * 2 + 1);

	MScalar = MvarMVSplitScalar(MCrvA[i]);

	MCrvAScalar[i][1] = MScalar[1];
	MCrvAScalar[i][2] = MScalar[2];
	
	if (MScalar[3] != NULL)
	    MvarMVFree(MScalar[3]);

	MScalar = MvarMVSplitScalar(MCrvB[i]);

	MCrvBScalar[i][1] = MScalar[1];
	MCrvBScalar[i][2] = MScalar[2];

	if (MScalar[3] != NULL)
	    MvarMVFree(MScalar[3]);

	MScalar = MvarMVSplitScalar(MDCrvA[i]);

	MDCrvAScalar[i][1] = MScalar[1];
	MDCrvAScalar[i][2] = MScalar[2];

	if (MScalar[3] != NULL)
	    MvarMVFree(MScalar[3]);

	MScalar = MvarMVSplitScalar(MDCrvB[i]);

	MDCrvBScalar[i][1] = MScalar[1];
	MDCrvBScalar[i][2] = MScalar[2];

	if (MScalar[3] != NULL)
	    MvarMVFree(MScalar[3]);
    }

    MScalar = MvarMVSplitScalar(MCircle);

    MCircleScalar[0] = MScalar[0];
    MCircleScalar[1] = MScalar[1];
    MCircleScalar[2] = MScalar[2];
    
    if (MScalar[3] != NULL)
	MvarMVFree(MScalar[3]);

    /* Tangency constraint. */

    for (i = 0; i < 3; ++i) {
	MTemp1 = MvarMVMult(MDCrvAScalar[i][1], MCircleScalar[1]);
	MTemp2 = MvarMVMult(MDCrvAScalar[i][2], MCircleScalar[2]);

	MTemp3 = MvarMVSub(MTemp1, MTemp2);

	MvarMVFree(MTemp1);
	MvarMVFree(MTemp2);

	MTemp1 = MvarMVMult(MDCrvAScalar[i][1], MCircleScalar[2]);
	MTemp2 = MvarMVMult(MDCrvAScalar[i][2], MCircleScalar[1]);

	MTemp4 = MvarMVAdd(MTemp1, MTemp2);

	MvarMVFree(MTemp1);
	MvarMVFree(MTemp2);

	MTemp1 = MvarMVMult(MDCrvBScalar[i][2], MTemp3);
	MTemp2 = MvarMVMult(MDCrvBScalar[i][1], MTemp4);

	MVs[i] = MvarMVSub(MTemp1, MTemp2);

	MvarMVFree(MTemp1);
	MvarMVFree(MTemp2);

	MTemp1 = MvarMVMult(MDCrvBScalar[i][1], MTemp3);
	MTemp2 = MvarMVMult(MDCrvBScalar[i][2], MTemp4);

	MVs[i + 7] = MvarMVAdd(MTemp1, MTemp2);

	MvarMVFree(MTemp1);
	MvarMVFree(MTemp2);
	MvarMVFree(MTemp3);
	MvarMVFree(MTemp4);

    }

    /* Intersection Constraints. */

    for (i = 0; i < 2; ++i) {
	MTemp1 = MvarMVSub(MCrvAScalar[i][1], MCrvAScalar[i+1][1]);
	MTemp2 = MvarMVSub(MCrvAScalar[i][2], MCrvAScalar[i+1][2]);

	MTemp3 = MvarMVMult(MTemp1, MCircleScalar[1]);
	MTemp4 = MvarMVMult(MTemp2, MCircleScalar[2]);

	MTemp5 = MvarMVSub(MTemp3, MTemp4);

	MvarMVFree(MTemp3);
	MvarMVFree(MTemp4);

	MTemp3 = MvarMVMult(MTemp1, MCircleScalar[2]);
	MTemp4 = MvarMVMult(MTemp2, MCircleScalar[1]);

	MTemp6 = MvarMVAdd(MTemp3, MTemp4);

	MvarMVFree(MTemp1);
	MvarMVFree(MTemp2);
	MvarMVFree(MTemp3);
	MvarMVFree(MTemp4);

	MTemp1 = MvarMVSub(MCrvBScalar[i][1], MCrvBScalar[i + 1][1]);

	MTemp3 = MvarMVMult(MTemp1, MCircleScalar[0]);

	MTemp2 = MvarMVSub(MCrvBScalar[i][2], MCrvBScalar[i + 1][2]);

	MTemp4 = MvarMVMult(MTemp2, MCircleScalar[0]);

	MVs[3 + 2 * i]= MvarMVSub(MTemp3, MTemp5);
	MVs[3 + 2 * i + 1]= MvarMVSub(MTemp4, MTemp6);

	MvarMVFree(MTemp1);
	MvarMVFree(MTemp2);
	MvarMVFree(MTemp3);
	MvarMVFree(MTemp4);
	MvarMVFree(MTemp5);
	MvarMVFree(MTemp6);
    }

    MPts = MvarMVsZeros0D(MVs, Constr, 10, Subtol, Numerictol);

    CagdCrvFree(Circle);
    CagdCrvFree(Subcircle);

    for (i = 0; i < 3; ++i) {
	CagdCrvFree(CrvA[i]);
	CagdCrvFree(CrvB[i]);
	MvarMVFree(MCrvA[i]);
	MvarMVFree(MCrvB[i]);
	MvarMVFree(MDCrvA[i]);
	MvarMVFree(MDCrvB[i]);

	MvarMVFree(MCrvAScalar[i][1]);
	MvarMVFree(MCrvAScalar[i][2]);
	MvarMVFree(MCrvBScalar[i][1]);
	MvarMVFree(MCrvBScalar[i][2]);

	MvarMVFree(MDCrvAScalar[i][1]);
	MvarMVFree(MDCrvAScalar[i][2]);
	MvarMVFree(MDCrvBScalar[i][1]);
	MvarMVFree(MDCrvBScalar[i][2]);
	MvarMVFree(MCircleScalar[i]);
    }

    MvarMVFree(MCircle);

    for (i = 0; i < 10; ++i)
	MvarMVFree(MVs[i]);

    return MPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Check how many rotational extreme points in the domain. If there are less  *
*  than 2 extreme points, there exist no loop in the domain.                 *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1:  Curve pair (moving curve, obstacle curve) for 1st contact point. *
*   Node2:  Curve pair (moving curve, obstacle curve) for 2st contact point. *
*   RotExtremes: Rotational extreme points for 2 contact traces.             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: True if there is no loop exist in the domain, False otherwise.*
*****************************************************************************/
static CagdBType Mvar2CtNoLoopTest(Mvar2CtPairStruct *Node1,
			           Mvar2CtPairStruct *Node2,
			           MvarPtStruct *RotExtremes)
{
    MvarPtStruct *TPt;
    int Count = 0, 
	Id1 = Node1 -> BNode -> Id,
	Id2 = Node2 -> BNode -> Id;
    
    for (TPt = RotExtremes; TPt != NULL; TPt = TPt -> Pnext) {	
	if (Id1 != AttrGetIntAttrib(TPt -> Attr, "First") ||
	    Id2 != AttrGetIntAttrib(TPt -> Attr, "Second"))
	   continue;
        /* Count the number of rotational extreme in the domain. */
	if (Mvar2CtInPairDomain(Node1, Node1, TPt))
	    Count++;
    }

    if (Count < 2)
	return TRUE;
    else
	return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Compute 2 contact traces. We traverse curve pair struct (Mvar2ctPairStruct)*
* and once it reaches leaf node we switch to equation solver.                *
*                                                                            *
* PARAMETERS:                                                                *
*   Node1:  Curve pair (moving curve, obstacle curve) for 1st contact point. *
*   Node2:  Curve pair (moving curve, obstacle curve) for 2st contact point. *
*   BvhA:     Bvh for moving curve.                                   	     *
*   BvhBs:    Bvh for obstacle curves.                                       *
*   RotExtremes: Rotational extreme points for 2 contact traces.             *
*   Step:         Step size to use in the numeric tracing.                   *
*   Subtol:    Subdivision tolerance for 2 contact trace computation.        *
*   Numerictol:     Numerical tolerance for 2 contact trace computation.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: 2contact motion curves.                            *
*****************************************************************************/
static MvarPolylineStruct *Mvar2CtCompute2CtTrace(Mvar2CtPairStruct *Node1,
				                  Mvar2CtPairStruct *Node2,
				                  Mvar2CtBVHStruct *BvhA,
				                  Mvar2CtBVHStruct **BvhBs,
				                  MvarPtStruct *RotExtremes,
						  CagdRType Step,
				                  CagdRType Subtol,
				                  CagdRType Numerictol)
{   
    int Dir, c;
    CagdRType Param;
    MvarPolylineStruct *MPoly, 
        *LPoly = NULL, 
        *RPoly = NULL;
    Mvar2CtBVNodeStruct *Nodes[4];
    CagdBType NoLoop;

    if (Node1 == Node2) {
	if (Node1 -> Bsingle) 
	    /* If this pair has at most one contact */
	    /* there exist no 2 contact. */
	    return NULL;

	if (Node1 -> Left == NULL && Node1 -> Right == NULL) {
	 
	    /* No loop test. */
	    NoLoop = Mvar2CtNoLoopTest(Node1, Node2, RotExtremes);
	    return Mvar2CtSolve2CtTrace(Node1, Node2, BvhA, BvhBs, 
		                        NoLoop, Step, Subtol, Numerictol);
	}

	if (Node1 -> Left != NULL)
	    LPoly =  Mvar2CtCompute2CtTrace(Node1 -> Left, Node2, BvhA, 
	                        BvhBs, RotExtremes, Step, Subtol, Numerictol);
	if (Node1 -> Right!=NULL)
	    RPoly =  Mvar2CtCompute2CtTrace(Node1 -> Right, Node2, BvhA, 
	                        BvhBs, RotExtremes, Step, Subtol, Numerictol);   
	
	/* Get subdivision direction and parameter. */
	Mvar2CtGetPairSubdivParam(Node1, Node2, &c, &Param); 

	/* Link left and right. */
	return MvarMVZR1DLinkNeighbours(LPoly, RPoly, Subtol, 
	                                2.0 * Numerictol, c, Param);
    }

    if (Mvar2CtTrivialReject(Node1, Node2)) /* Reject trivial cases. */
	return NULL;  

    if (Mvar2CtIsConnectedPair(Node1, Node2)) {
	/* If connected pairs are both convex or both concave              */
	/* we don't have 2 contact.                                        */ 
	if (Node1 -> ANode -> Convexity * Node1 -> BNode -> Convexity > 0  &&
	    Node2 -> ANode -> Convexity * Node2 -> BNode -> Convexity > 0)
	    return NULL;

	Dir = Mvar2CtPairSubdivDir2(Node1, Node2);

	switch (Dir) {
	    case 0: 
	        if (Node1 -> Left != NULL)
		    LPoly =  Mvar2CtCompute2CtTrace(Node1 -> Left, Node2, 
		                                    BvhA, BvhBs, 
						    RotExtremes, Step, 
						    Subtol, Numerictol);
	        if (Node1 -> Right != NULL)
		    RPoly =  Mvar2CtCompute2CtTrace(Node1 -> Right, Node2, 
		                                    BvhA, BvhBs, 
						    RotExtremes, Step, 
						    Subtol, Numerictol);   	
	        /* Get subdivision direction and parameter. */
	        Mvar2CtGetPairSubdivParam(Node1, Node2, &c, &Param);

	        /* Link left and right. */
	        return MvarMVZR1DLinkNeighbours(LPoly, RPoly, Subtol, 
		                            2.0 * Numerictol, c, Param);
	    case 1:
	        if (Node2 -> Left != NULL)
		    LPoly =  Mvar2CtCompute2CtTrace(Node1, Node2 -> Left, 
		                                    BvhA, BvhBs, 
						    RotExtremes, Step, 
						    Subtol, Numerictol);
	        if (Node2 -> Right != NULL)
		    RPoly =  Mvar2CtCompute2CtTrace(Node1, Node2 -> Right, 
		                                    BvhA, BvhBs, 
						    RotExtremes, Step, 
						    Subtol, Numerictol);

	        /* Get subdivision direction and parameter. */
	        Mvar2CtGetPairSubdivParam(Node1, Node2, &c, &Param);
	        /* Link left and right. */
	        return MvarMVZR1DLinkNeighbours(LPoly, RPoly, Subtol, 
		                            2.0 * Numerictol, c, Param);
	    case 2:	   
	        /* We reached leaf node of Node1 and Node2. */
	        Nodes[0] = Node1 -> ANode;
	        Nodes[1] = Node1 -> BNode;
	        Nodes[2] = Node2 -> ANode;
	        Nodes[3] = Node2 -> BNode;

	        /* Check if there exist 2 contact trace using BVH. */ 
	        if (!Mvar2CtCheck2CtTrace(Nodes, Node1 -> Cparam, Subtol))
		    return NULL;
	        /* No loop test. */
	        NoLoop = Mvar2CtNoLoopTest(Node1, Node2, RotExtremes);

	        /* Solve constraint equation for 2 contact trace. */
	        MPoly = Mvar2CtSolve2CtTrace(Node1, Node2, BvhA, BvhBs, 
		                          NoLoop, Step, Subtol, Numerictol);
	        return MPoly;	   
	}
    }
    else {
	if (Mvar2CtAABBOverlap(&(Node1 -> AABB), &(Node2 -> AABB))) {
	    Dir = Mvar2CtPairSubdivDir2(Node1, Node2);
	    switch (Dir) {
	        case 0: 
		    if (Node1 -> Left != NULL)
		        LPoly =  Mvar2CtCompute2CtTrace(Node1 -> Left, Node2, 
			                                BvhA, BvhBs, 
							RotExtremes, Step, 
							Subtol, Numerictol);
		    if (Node1 -> Right != NULL)
		        RPoly =  Mvar2CtCompute2CtTrace(Node1 -> Right, Node2, 
			                                BvhA, BvhBs, 
							RotExtremes, Step, 
							Subtol, Numerictol);
	            /* Get subdivision direction and parameter. */
		    Mvar2CtGetPairSubdivParam(Node1, Node2, &c, &Param);
		    /* Link left and right. */
		    return MvarMVZR1DLinkNeighbours(LPoly, RPoly, Subtol, 
		                                    2.0 * Numerictol, c, 
						    Param);

	        case 1:
		    if (Node2 -> Left != NULL)
		        LPoly =  Mvar2CtCompute2CtTrace(Node1, Node2 -> Left, 
		                                        BvhA, BvhBs, 
						        RotExtremes, Step,
							Subtol, Numerictol);
		    if (Node2 -> Right != NULL)
		        RPoly =  Mvar2CtCompute2CtTrace(Node1, Node2 -> Right, 
		                                        BvhA, BvhBs, 
						        RotExtremes, Step, 
							Subtol, Numerictol);
		    /* Get subdivision direction and parameter. */
		    Mvar2CtGetPairSubdivParam(Node1, Node2, &c, &Param);

		    /* Link left and right. */
		    return MvarMVZR1DLinkNeighbours(LPoly, RPoly, Subtol, 
		                                    2.0 * Numerictol, c, 
						    Param);
	        case 2:
		    /* We reached leaf node of Node1 and Node2. */
		    Nodes[0] = Node1 -> ANode;
		    Nodes[1] = Node1 -> BNode;
		    Nodes[2] = Node2 -> ANode;
		    Nodes[3] = Node2 -> BNode;		

		    /* Check if there exist 2 contact trace using BVH. */ 
		    if (!Mvar2CtCheck2CtTrace(Nodes, Node1 -> Cparam, Subtol))
		        return NULL;
		    /* No loop test. */
		    NoLoop = Mvar2CtNoLoopTest(Node1, Node2, RotExtremes);
		    /* Solve constraint equation for 2 contact trace. */
		    MPoly = Mvar2CtSolve2CtTrace(Node1, Node2, BvhA, BvhBs, 
		                                 NoLoop, Step, Subtol, Numerictol);
		   
		    return MPoly;
	    }	    
	}
	return NULL;
    }
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Compute 3 contact points from triple of 2 contact traces. Direct attempt   *
*  to solve 3 contact points for entire domain is quite inefficient. Instead *
*  we find a common domain of 2 contact traces and solve 3 contact point     * 
*  from it.                                                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   BvhA:     Bvh for moving curve.                                   	     *
*   BvhBs:    Bvh for obstacle curves.                                       *
*   SPts:     start points of the 2 contact traces.                          *
*   EPts:     end point of the 2 contact traces.             		     *
*   CRoot:    Root of Mvar2CtCparam struct.            		             *
*   Ids:     Ids of obstacle curves.                                         *
*   Triples:  Linked list of resulting 3 contact points.                     *
*   Subtol:    Subdivision tolerance for 3Contact computation.               *
*   Numerictol:     Numerical tolerance for 3 contact computation.           *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtCompute3CtPoint(Mvar2CtBVHStruct *BvhA, 
		                   Mvar2CtBVHStruct **BvhBs, 
		                   MvarPtStruct *SPts[3], 
		                   MvarPtStruct *EPts[3],
			           Mvar2CtCParamStruct *CRoot,
		                   int Ids[3], 
		                   MvarPtStruct **Triples, 
		                   CagdRType Subtol, 
		                   CagdRType Numerictol)
{

    IRIT_STATIC_DATA CagdRType BBoxMin[3][5], BBoxMax[3][5], Min[7], Max[7],
	                       Min2[7], Max2[7];
    IRIT_STATIC_DATA int Length[3];
    MvarPtStruct *TPt, *MidPts[3], *MPt;
    Mvar2CtBVNodeStruct *Nodes[6];
    Mvar2CtCParamStruct *Cparam;
    CagdRType 
	MaxLength = 0;
    int i, j,
	MaxDir = 0;

    for (i = 0; i < 3; ++i)
	/* Compute bounding box of 2 contact traces. */
	Mvar2CtTraceBBox(BBoxMin[i], BBoxMax[i], SPts[i], EPts[i]);

    for (i = 0; i < 3; ++i) { 
	/* If there is no common domain, we don't have 3 contact. */
	if (BBoxMin[i][4] > BBoxMax[(i + 1) % 3][4])
	    return;
	if (BBoxMin[(i + 1) % 3][4] > BBoxMax[i][4])
	    return;
    }

    for (i = 0; i < 3; ++i) {
	for (j = 0; j < 2; ++j) {
	    /* If there is no common domain, we don't have 3 contact. */
	    if (BBoxMin[i][j + 2] > BBoxMax[(i + 1) % 3][j])
		return;
	    if (BBoxMin[(i + 1) % 3][j] > BBoxMax[i][j + 2])
		return;
	}
    }
    	
    for (i = 0; i < 3; ++i) {
	for (j = 0; j < 2; ++j) {
	    /* Extract domain. */
	    Min[i * 2 + j] = IRIT_MIN(BBoxMin[i][j], 
		                      BBoxMin[(i + 2) % 3][j + 2]);
	    Max[i * 2 + j] = IRIT_MAX(BBoxMax[i][j], 
		                      BBoxMax[(i + 2) % 3][j + 2]);
	}
	
	Length[i] = 0;
	
	for(TPt = SPts[i]; TPt != EPts[i]; TPt = TPt -> Pnext)
	    Length[i]++;
    }

    Min[6] = IRIT_MIN3(BBoxMin[0][4], BBoxMin[1][4], BBoxMin[2][4]);
    Max[6] = IRIT_MAX3(BBoxMax[0][4], BBoxMax[1][4], BBoxMax[2][4]);
    
    for (i = 0; i < 7; ++i) {	
	if (MaxLength < Max[i] - Min[i]) {
	    MaxLength = Max[i] - Min[i];
	    MaxDir = i;
	}
    }

    if (IRIT_MAX3(Length[0], Length[1], Length[2]) < 3) {	
	/* If the domain is small enough. */
	for (i = 0; i < 3; ++i) {
	    Mvar2CtGetParentBVNode(BvhA -> Root, Min[i * 2], 
		                   Max[i * 2], &(Nodes[2 * i]));
	    Mvar2CtGetParentBVNode(BvhBs[Ids[i]] -> Root, Min[i * 2 + 1], 
		                   Max[i * 2 + 1], &(Nodes[2 * i + 1]));
	}

	Mvar2CtGetParentCparam(CRoot, Min[6], Max[6],  &Cparam);

	for (i = 0; i < 7; ++i) {
	    Min2[i] = IRIT_INFNTY;
	    Max2[i] = -IRIT_INFNTY;
	}

	/* We further reduce the domain using BVH. */
	Mvar2CtReduce3CtDomain(Nodes, Cparam, Min2, Max2); 

	if (Min2[0] != IRIT_INFNTY) {
	    /* Solve 3 contact points */
	    for (i = 0; i < 7; ++i) {
		Min[i] = IRIT_MAX(Min[i], Min2[i]);
		Max[i] = IRIT_MIN(Max[i], Max2[i]);
		
		if (Min[i] >= Max[i])
		    return;
	    }
	    
	   MPt = Mvar2CtSolve3CtPoint(BvhA, BvhBs, Ids, Min, Max, 
	                              Subtol, Numerictol);
	    (*Triples) = (MvarPtStruct *) CagdListAppend(MPt, (*Triples));
	}
    }
    else {
	/* Subdivide 2 contact traces and repeat the procedure. */

	if (Length[0] > Length[1] && Length[0] > Length[2]) {
	    MidPts[0] = Mvar2CtGetMiddlePt(SPts[0], Length[0]);
	    MidPts[1] = EPts[1];
	    MidPts[2] = EPts[2];

	    Mvar2CtCompute3CtPoint(BvhA, BvhBs, SPts, MidPts, CRoot, 
		                   Ids, Triples, Subtol, Numerictol);

	    MidPts[1] = SPts[1];
	    MidPts[2] = SPts[2];

	    Mvar2CtCompute3CtPoint(BvhA, BvhBs, MidPts, EPts, CRoot, 
		                   Ids, Triples, Subtol, Numerictol);
	}

	else if (Length[1] > Length[2]) {
	    MidPts[1] = Mvar2CtGetMiddlePt(SPts[1], Length[1]);
	    MidPts[0] = EPts[0];
	    MidPts[2] = EPts[2];

	    Mvar2CtCompute3CtPoint(BvhA, BvhBs, SPts, MidPts, CRoot, 
		                   Ids, Triples, Subtol, Numerictol);

	    MidPts[0] = SPts[0];
	    MidPts[2] = SPts[2];

	    Mvar2CtCompute3CtPoint(BvhA, BvhBs, MidPts, EPts, CRoot, 
		                   Ids, Triples, Subtol, Numerictol);

	}
	else {
	    MidPts[2] = Mvar2CtGetMiddlePt(SPts[2], Length[2]);
	    MidPts[0] = EPts[0];
	    MidPts[1] = EPts[1];

	    Mvar2CtCompute3CtPoint(BvhA, BvhBs, SPts, MidPts, CRoot, 
		                   Ids, Triples, Subtol, Numerictol);

	    MidPts[0] = SPts[0];
	    MidPts[1] = SPts[1];
	    
	    Mvar2CtCompute3CtPoint(BvhA, BvhBs, MidPts, EPts, CRoot, 
		                   Ids, Triples, Subtol, Numerictol);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Insert a given 3 contact point into 2 contact trace.                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly: 2 contact traces to divide.                                        *
*   Pt:  3 contact point to insert.                                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtInsert3CtPoint(MvarPolylineStruct *Poly, MvarPtStruct *Pt)
{
    MvarVecStruct *Diff1, *Diff2;
    MvarPtStruct *TPt, *NewPt;
    int i;
    CagdRType Dist;

    Diff1 = MvarVecNew(5);
    Diff2 = MvarVecNew(5);

    for (TPt = Poly -> Pl; TPt -> Pnext != NULL; TPt = TPt -> Pnext) {
	for (i = 0; i < 5; ++i) {
	    Diff1 -> Vec[i] = Pt -> Pt[i] - TPt -> Pt[i];
	    Diff2 -> Vec[i] = Pt -> Pt[i] - TPt -> Pnext -> Pt[i];
	}

	Dist = MvarPtDistTwoPoints(TPt, TPt -> Pnext);

	/* Check if Pt is located between TPt and TPt->Pnext. */ 
	if (MvarVecLength(Diff1) < Dist && MvarVecLength(Diff2) < Dist 
	    && MvarVecDotProd(Diff1, Diff2) < 0) {

	    NewPt = MvarPtCopy(Pt);

	     /* Mark it as 3 contact. */
	    AttrSetIntAttrib(&(NewPt -> Attr), "TripleContact", 1);

	    NewPt -> Pnext = TPt -> Pnext;

	    TPt -> Pnext = NewPt;

	    MvarVecFree(Diff1);
	    MvarVecFree(Diff2);

	    return;
	}
    }

    MvarVecFree(Diff1);
    MvarVecFree(Diff2);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Insert a given curvature contact point into curvature trace.               *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly: curvature traces to divide.                                        *
*   Pt:  curvature contact point to insert.                                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Mvar2CtInsertCurvPoint(MvarPolylineStruct *Polys, MvarPtStruct *Pt)
{
    MvarVecStruct *Diff1, *Diff2;
    MvarPtStruct *TPt, *NewPt;
    MvarPolylineStruct *TPoly;
    int i;
    CagdRType Dist;

    Diff1 = MvarVecNew(2);
    Diff2 = MvarVecNew(2);


    for (TPoly = Polys; TPoly != NULL; TPoly = TPoly -> Pnext) {
	for (TPt = TPoly -> Pl; TPt -> Pnext != NULL; TPt = TPt -> Pnext) {
	    for (i = 0; i < 2; ++i) {
		Diff1 -> Vec[i] = Pt -> Pt[i] - TPt -> Pt[i];
		Diff2 -> Vec[i] = Pt -> Pt[i] - TPt -> Pnext -> Pt[i];
		}

		Dist = MvarPtDistTwoPoints(TPt, TPt -> Pnext);

		/* Check if Pt is located between TPt and TPt->Pnext. */ 
		if (MvarVecLength(Diff1) < Dist && MvarVecLength(Diff2) < Dist 
		    && MvarVecDotProd(Diff1, Diff2) < 0) {

		    NewPt = MvarPtNew(2);
		    NewPt->Pt[0] = Pt->Pt[0];
		    NewPt->Pt[1] = Pt->Pt[1];
	     /* Mark it as 3 contact. */
	    AttrSetIntAttrib(&(NewPt -> Attr), "CurvContact", 1);

		    NewPt -> Pnext = TPt -> Pnext;

		    TPt -> Pnext = NewPt;
	   
		    break;
		    }
	    }
    }
    MvarVecFree(Diff1);
    MvarVecFree(Diff2);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Divide Curvature traces at Curvature Contact points.                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Polys: curvature traces to divide.                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: curvature traces divided at curvature contacts.    *
*****************************************************************************/
static MvarPolylineStruct *Mvar2CtDivideAtCurvPointAux(MvarPolylineStruct *Polys)
{
    MvarPolylineStruct *TPoly, *NewPoly, 
	*Polyout = NULL;
    MvarPtStruct *TPt, *NewPt;

    while (Polys != NULL) {
	IRIT_LIST_POP(TPoly, Polys);

	for (TPt = TPoly -> Pl; TPt != NULL; TPt = TPt -> Pnext) {
            /* Is it 3 contact? */
	    if (AttrGetIntAttrib(TPt -> Attr, "CurvContact") == 1) { 
	        /* Divide it into two traces. */
		if (TPt -> Pnext != NULL && TPt != TPoly -> Pl) { 
		    NewPt = MvarPtCopy(TPt);
		    NewPt -> Pnext = TPt -> Pnext;
		    TPt -> Pnext = NULL;		   
		    NewPoly = MvarPolylineNew(NewPt);
		    IRIT_LIST_PUSH(NewPoly, Polys);

		    IRIT_LIST_PUSH(TPoly, Polyout);

		    TPoly = NULL;

		    break;
		}
	    }
	}

	if (TPoly != NULL)
	    IRIT_LIST_PUSH(TPoly, Polyout);
    }

    return Polyout;
}


/*****************************************************************************
* DESCRIPTION:                                                               *
* Divide 2Contact traces at 3 Contact points. Find a 3 contact inserted in   *
* 2 contact trace and divide it into two traces.                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Polys: 2 contact traces to divide.                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: 2contact motion curves divided at 3 contact points. *
*****************************************************************************/
static MvarPolylineStruct *Mvar2CtDivideAt3CtPointAux(MvarPolylineStruct *Polys)
{
    MvarPolylineStruct *TPoly, *NewPoly, 
	*Polyout = NULL;
    MvarPtStruct *TPt, *NewPt;
    int First, Second;

    while (Polys != NULL) {
	IRIT_LIST_POP(TPoly, Polys);

	for (TPt = TPoly -> Pl; TPt != NULL; TPt = TPt -> Pnext) {
            /* Is it 3 contact? */
	    if (AttrGetIntAttrib(TPt -> Attr, "TripleContact") == 1) { 
	        /* Divide it into two traces. */
		if (TPt -> Pnext != NULL && TPt != TPoly -> Pl) { 
		    NewPt = MvarPtCopy(TPt);
		    NewPt -> Pnext = TPt -> Pnext;

		    TPt -> Pnext = NULL;
		    First = AttrGetIntAttrib(TPoly -> Attr, "First");
		    Second = AttrGetIntAttrib(TPoly -> Attr, "Second");

		    NewPoly = MvarPolylineNew(NewPt);

		    AttrSetIntAttrib(&(NewPoly -> Attr), "First", First);
		    AttrSetIntAttrib(&(NewPoly -> Attr), "Second", Second);

		    IRIT_LIST_PUSH(NewPoly, Polys);

		    IRIT_LIST_PUSH(TPoly, Polyout);

		    TPoly = NULL;

		    break;
		}
	    }
	}

	if (TPoly != NULL)
	    IRIT_LIST_PUSH(TPoly, Polyout);
    }

    return Polyout;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Divide 2Contact traces at 3 Contact points. Test every triple of 2 contact *
* traces for computing 3 contact point. If there is a 3 contact, divide a    *
* traces at the 3 contact points.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   BvhA:     Bvh for moving curve.                                   	     *
*   BvhBs:    Bvh for obstacle curves.                                       *
*   CRoot:    Root of Mvar2CtCparam struct.            		             *
*   M2CtTraces: Linked list of 2 contact Traces.             		     *
*   Subtol:    Subdivision tolerance for 3Contact computation.               *
*   Numerictol:     Numerical tolerance for 3 contact computation.           *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: 2contact motion curves divided at 3 contact points. *
*****************************************************************************/
static MvarPolylineStruct *Mvar2CtDivideAt3CtPoint(Mvar2CtBVHStruct *BvhA,
			                           Mvar2CtBVHStruct **BvhBs,
			                           Mvar2CtCParamStruct *CRoot,
			                           MvarPolylineStruct *M2CtTraces,
			                           CagdRType Subtol,
			                           CagdRType Numerictol)
{   
     int Type, Type2, Type3, Ids[3];   		
     MvarPolylineStruct *TPoly1, *TPoly2, *TPoly3;
     MvarPtStruct *MVPt, *SPts[3], *EPts[3], 
	 *Triples = NULL;
    /* Check every triple (TPoly1, TPoly2, TPoly3). */
    for (TPoly1 = M2CtTraces; TPoly1 != NULL; TPoly1 = TPoly1 -> Pnext) { 
        for (TPoly2 = TPoly1 -> Pnext; 
	     TPoly2 != NULL; TPoly2 = TPoly2 -> Pnext) {
	     /* Check if two traces have common contact point. */  
	    Type = Mvar2CtTraceCollide(TPoly1, TPoly2); 

	    if (!Type)
		continue;

	    switch (Type) { 
	        /* We swap the traces so that common contact point becomes */  
	        /* the second contact point of TPoly1 and first contact    */ 
	        /* point of TPoly2.                                        */
	        case 1:
		    Mvar2CtSwapTrace(TPoly1); 
		    break;
	        case 2:
		    Mvar2CtSwapTrace(TPoly1);
		    Mvar2CtSwapTrace(TPoly2);
		    break;
	        case 3:
		    break;
	        case 4:
		    Mvar2CtSwapTrace(TPoly2);
		    break;
	    }
    
	    for (TPoly3 = TPoly2 -> Pnext;
		 TPoly3 != NULL; 
		 TPoly3 = TPoly3 -> Pnext) {
	        /* Check if two traces have common configuration. */ 
	        Type2 = Mvar2CtTraceCollide(TPoly1, TPoly3);  
		Type3 = Mvar2CtTraceCollide(TPoly2, TPoly3);
	      
		if (!Type2 || !Type3)
		    continue;	
		if (Type2 == 3 || Type2 == 4)
		    continue;
		if (Type3 == 1 || Type3 == 2)
		    continue;

		if (Type3 == 4)
		    Mvar2CtSwapTrace(TPoly3);

		/* Find ids of obstacle curves. */   
		Ids[0] = AttrGetIntAttrib(TPoly1 -> Attr, "First"); 
		Ids[1] = AttrGetIntAttrib(TPoly2 -> Attr, "First");
		Ids[2] = AttrGetIntAttrib(TPoly3 -> Attr, "First"); 
		
		SPts[0] = TPoly1 -> Pl;
		SPts[1] = TPoly2 -> Pl;
		SPts[2] = TPoly3 -> Pl;

		EPts[0] = (MvarPtStruct *) CagdListLast(SPts[0]);
		EPts[1] = (MvarPtStruct *) CagdListLast(SPts[1]);
		EPts[2] = (MvarPtStruct *) CagdListLast(SPts[2]);

		Triples = NULL;
		/* Compute 3 contact points. */
		Mvar2CtCompute3CtPoint(BvhA, BvhBs, SPts, EPts, CRoot, 
				       Ids, &Triples, Subtol, Numerictol);

		/* Remove duplicates. */
		Triples = MvarZeroFilterSolutionSet(Triples, NULL, NULL, 0,  
						    Numerictol * 10.0, FALSE,
						    TRUE, FALSE);

		if (Triples != NULL) {	   
		    for (MVPt = Triples; MVPt != NULL; MVPt = MVPt -> Pnext) {			
		        MvarPtStruct
			    *MPt = MvarPtNew(5);

			MPt -> Pt[0] = MVPt -> Pt[0];
			MPt -> Pt[1] = MVPt -> Pt[1];
			MPt -> Pt[2] = MVPt -> Pt[2];
			MPt -> Pt[3] = MVPt -> Pt[3];
			MPt -> Pt[4] = MVPt -> Pt[6];

			/* Insert the 3 contact point into trace. */
			Mvar2CtInsert3CtPoint(TPoly1, MPt); 

			MPt -> Pt[0] = MVPt -> Pt[2];
			MPt -> Pt[1] = MVPt -> Pt[3];
			MPt -> Pt[2] = MVPt -> Pt[4];
			MPt -> Pt[3] = MVPt -> Pt[5];
			MPt -> Pt[4] = MVPt -> Pt[6];
			
			/* Insert the 3 contact point into trace. */
			Mvar2CtInsert3CtPoint(TPoly2, MPt); 

			MPt -> Pt[0] = MVPt -> Pt[4];
			MPt -> Pt[1] = MVPt -> Pt[5];
			MPt -> Pt[2] = MVPt -> Pt[0];
			MPt -> Pt[3] = MVPt -> Pt[1];
			MPt -> Pt[4] = MVPt -> Pt[6];

			/* insert the 3 contact point into trace. */
			Mvar2CtInsert3CtPoint(TPoly3, MPt); 

			MvarPtFree(MPt);
		    }

		    MvarPtFreeList(Triples);
		}    
	    }
	}
    }

    return Mvar2CtDivideAt3CtPointAux(M2CtTraces);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the 2 Contact motion between two C^1 cont curves, CCrvA and     M 
*   CCrvB.                                                                   M
*   The curves are assumed to be C1 periodic curve with open end condition.  M
*   The algebraic conditions for 2contact are following (please refer the    M
*    paper,"Precise Continuous Contact Motion for Planar Freeform Geometric  M
*    Curves" for the details):                                               M
*                                                                            M
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_x                                       V 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_x     V 
*  [CrvB(v1) -Rot(theta) * CrvA(u1)]_y                                       V 
*                                 = [CrvB(v2) - Rot(theta) * CrvA(u2)]_y     V 
*    det (Rot(theta) * CrvA'(u1), CrvB'(v1) = 0,                             V
*    det (Rot(theta) * CrvA'(u2), CrvB'(v2) = 0.                             V
*                                                                            *
* PARAMETERS:                                                                M
*   CCrvA:     Moving curve.	                                             M
*   CCrvB:     Obstacle Curves.                                   	     M
*   Step:      Step size to use in the numeric tracing.                      M
*   Subtol:    Subdivision tolerance of Computation.                         M
*   Numerictol:   Numerical tolerance of computation.		             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: Linked list of the solutions holding parameter     M
*                         value of two contact points and rotation in        M
*                         radian (u1, v1, u2, v2, theta).		     M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Mvar2CntctCompute2CntctMotion, 2Contact Trace                            M
*****************************************************************************/
MvarPolylineStruct *Mvar2CntctCompute2CntctMotion(const CagdCrvStruct *CCrvA, 
					          const CagdCrvStruct *CCrvB, 
						  CagdRType Step,
					          CagdRType Subtol, 
					          CagdRType Numerictol)
{
    CagdCrvStruct *CrvA, **CrvBs, *Circle;
    const CagdCrvStruct *TCrv;
    Mvar2CtBVHStruct *BvhA, **BvhBs;
    Mvar2CtCParamStruct *CRoot;
    Mvar2CtPairStruct **PRoots;
    MvarPtStruct *TPt, 
        *RotExtremes = NULL;
    MvarPolylineStruct *MPoly, *TPoly, 
        *M2CtTraces = NULL;
    CagdRType TMinA, TMaxA, *TMinB, *TMaxB, *R;
    CagdPType C;
    int i, j, First, Second, 
	BSize = CagdListLength(CCrvB); 

    /* Make sure curves have open end condition and are compatible. */
    if (CAGD_IS_BEZIER_CRV(CCrvA))
	CrvA = CagdCnvrtBzr2BspCrv(CCrvA);
    else
	CrvA = CagdCnvrtBsp2OpenCrv(CCrvA);

    CagdCrvDomain(CrvA, &TMinA, &TMaxA);
    BspKnotAffineTrans2(CrvA -> KnotVector, 
	                CrvA -> Order + CrvA -> Length, 0, 1);

     /* Build a BVH for CrvA. */
    BvhA = Mvar2CtBuildBVH(CrvA, IRIT_MIN(MVAR2CT_BV_SUBTOL, Subtol * 0.1), 
	                   MVAR2CT_BV_TOL);
    Mvar2CtSetNodeId(BvhA -> Root, 0);

    CrvBs = (CagdCrvStruct **) IritMalloc(sizeof(CagdCrvStruct*) * BSize);
    TMinB = (CagdRType *) IritMalloc(sizeof(CagdRType) * BSize);
    TMaxB = (CagdRType *) IritMalloc(sizeof(CagdRType) * BSize);
    BvhBs = 
	(Mvar2CtBVHStruct **) IritMalloc(sizeof(Mvar2CtBVHStruct*) * BSize);
    PRoots = 
	(Mvar2CtPairStruct**) IritMalloc(sizeof(Mvar2CtPairStruct*) * BSize);

    for (TCrv = CCrvB, i = 0; TCrv != NULL; TCrv = TCrv -> Pnext, ++i) {
	/* Make sure curves have open end condition and are compatible. */
	if (CAGD_IS_BEZIER_CRV(CCrvB))
	    CrvBs[i] = CagdCnvrtBzr2BspCrv(TCrv);
	else
	    CrvBs[i] = CagdCnvrtBsp2OpenCrv(TCrv);

	CagdCrvDomain(CrvBs[i], &(TMinB[i]), &(TMaxB[i]));
	BspKnotAffineTrans2(CrvBs[i] -> KnotVector, 
	                    CrvBs[i] -> Order + CrvBs[i] -> Length, 0, 1);
	/* Build a BVH for CrvB. */
	BvhBs[i] = Mvar2CtBuildBVH(CrvBs[i], 
	                          IRIT_MIN(MVAR2CT_BV_SUBTOL, Subtol * 0.1), 
	                          MVAR2CT_BV_TOL); 
	Mvar2CtSetNodeId(BvhBs[i] -> Root, i);
    }

    Circle = BspCrvCreateUnitCircle();
    BspKnotAffineTrans2(Circle -> KnotVector, 
	                Circle -> Order + Circle -> Length, 0, 1);
    
    Mvar2CtBuildCParamHierarchy(Circle, &CRoot, 0, 1, Subtol * 0.1);
    BspKnotAffineTrans2(CrvA -> KnotVector, CrvA -> Order + 
	                CrvA -> Length, 0, 1);   

    for (i = 0; i < BSize; ++i) { 
       /* Build data structure for domain purging. */	
        Mvar2CtBuildPairHierarchy(&(PRoots[i]), BvhA -> Root, 
	                          BvhBs[i] -> Root, BvhA, 
				  BvhBs, BSize, CRoot, 
				  MVAR2CT_PAIR_TOL);
    } 

    for (i = 0; i < BSize; ++i) { /* Compute curvature contact points. */
	if (PRoots[i] == NULL)
	    continue;
	Mvar2CtComputeCurvatureContact(PRoots[i], BvhA, BvhBs[i],
	                               &(BvhBs[i] -> CurvContacts), 
				       Subtol, Numerictol);
        /* Delete redundant solutions. */
	BvhBs[i] -> CurvContacts = 
	    Mvar2CtValidateCurvContact(BvhBs[i] -> CurvContacts, 
	                               BvhA, BvhBs, BSize, i, Circle);

    } 



    for (i = 0; i < BSize; ++i) { /* Identify single contact domains. */	
	if (PRoots[i] == NULL)
	    continue;
	Mvar2CtSetSingleContact(PRoots[i], BvhA, BvhBs[i],
	                        BvhBs[i] -> CurvContacts, 
				Subtol, Numerictol);
    } 


    for (i = 0; i < BSize; ++i) {	
	if (PRoots[i] == NULL)
	    continue;	
	for (j = i; j < BSize; ++j) {/* Compute rotation extreme points. */
	    if (PRoots[j] == NULL)
		continue;
	    Mvar2CtRotationExtreme(PRoots[i], PRoots[j], BvhA, BvhBs, 
	                           &RotExtremes, Subtol, Numerictol);
	}
    }
    RotExtremes = Mvar2CtValidate2Ct(RotExtremes, BvhA, BvhBs, BSize, Circle);

    for (i = 0; i < BSize; ++i) {
	if (PRoots[i] == NULL)
	    continue;	
	for (j = i; j < BSize; ++j) { /* Solve 2Ct traces. */
	    if (PRoots[j] == NULL)
		continue;
	    MPoly = Mvar2CtCompute2CtTrace(PRoots[i], PRoots[j], BvhA, BvhBs, 
	                                   RotExtremes, Step, Subtol, Numerictol); 
	  
	    if (MPoly != NULL) {                
	        for (TPoly = MPoly; TPoly != NULL; TPoly = TPoly -> Pnext) {
		    AttrSetIntAttrib(&(TPoly -> Attr), "First", i);
		    AttrSetIntAttrib(&(TPoly -> Attr), "Second", j);
		}
	       
		M2CtTraces = 
		      (MvarPolylineStruct *) CagdListAppend(MPoly, M2CtTraces);
	    }
	}
    }

    if (CagdListLength(M2CtTraces) >= 3)
	M2CtTraces = Mvar2CtDivideAt3CtPoint(BvhA, BvhBs, CRoot, 
	M2CtTraces, Subtol, Numerictol); /* Compute 3Ct and divide it. */

    /* Delete inter-penetrating trace. */
    M2CtTraces = 
 	   Mvar2CtValidateTraces(M2CtTraces, BvhA, BvhBs,  BSize, Circle); 
 
    M2CtTraces = 
	   Mvar2CtConnectPeriodic(M2CtTraces, Subtol); 

    for (TPoly = M2CtTraces; TPoly != NULL; TPoly = TPoly -> Pnext) {	
	First = AttrGetIntAttrib(TPoly -> Attr, "First");
	Second = AttrGetIntAttrib(TPoly -> Attr, "Second");

	/* Recover original domain. */
	for (TPt = TPoly -> Pl; TPt != NULL; TPt = TPt -> Pnext) { 
	    TPt -> Pt[0] = 
		TMinA + (TMaxA - TMinA) * TPt -> Pt[0];
	    TPt -> Pt[1] = 
		TMinB[First] + 
			     (TMaxB[First] - TMinB[Second]) * TPt -> Pt[1];
            TPt -> Pt[2] = 
		TMinA + (TMaxA - TMinA) * TPt -> Pt[2];	
	    TPt -> Pt[3] = 
		TMinB[Second] + 
		             (TMaxB[Second] - TMinB[Second]) * TPt -> Pt[3];
	    /* Compute rotation angle. */
	    R = CagdCrvEval(Circle, TPt -> Pt[4]); 
	    CagdCoerceToE2(C, &R, -1, Circle -> PType);
	    TPt -> Pt[4] = Mvar2CtGetTheta(C[0], C[1]);
	}	
    }

    MvarPtFreeList(RotExtremes);
    CagdCrvFree(Circle);

    Mvar2CtFreeBVH(BvhA);

    for (i = 0; i < BSize; ++i) {
	Mvar2CtFreeBVH(BvhBs[i]);
	Mvar2CtFreePair(PRoots[i]);
	CagdCrvFree(CrvBs[i]);
    }

    Mvar2CtFreeCparam(CRoot);

    IritFree(CrvBs); 
    IritFree(TMinB); 
    IritFree(TMaxB); 
    IritFree(BvhBs); 
    IritFree(PRoots); 
    CagdCrvFree(CrvA);

    return M2CtTraces;
}
