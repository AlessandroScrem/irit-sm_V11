/******************************************************************************
* ZrMV2DTp.c - tools to compute zero sets of multivariates when the problem   *
*	     has MVs representation and the expected solution  set is a       *
*            two-manifold: surface(s). This file handles the topological      *
*	     guarantee, i.e. termination of subdivision, without the          *
*	     tesselation/triangulation.					      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, June2014.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/geom_lib.h"

#define MVAR_ZRMV2D_POLY_CHAIN_TOL 1e-14 
#define MVAR_ZRMV2D_POLY_POLY_CLOSED_TOL 1e-12
#define MVAR_ZRMV2D_CORNER_PT_TOL 1e-12
/* #define MVAR_DBG_2D_TOPO_CONCLUDE */

IRIT_GLOBAL_DATA int
    _MVGlblZeroOutputTypeLoops = FALSE;

static CagdBType MvarZeroSolverNoZeroTest2DMVs(MvarZeroPrblmStruct *Problem,
					       int Ind);
static void MvarZeroSolveBoundary2DMVs(MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverTopoTest2DMVs(MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep2DMVs(
					         MvarZeroPrblmStruct *Problem);
static MvarZeroPrblmStruct **MvarZeroSubdivProblem2DMVs(
					     MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol2DMVs(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions2DMVs(
						MvarZeroSolutionStruct *Sol1,
						MvarZeroSolutionStruct *Sol2, 
					        MvarMVDirType Dir,
					        CagdRType Param,
						MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverOne2OneProjMVs(MvarZeroPrblmStruct *Problem);
static MvarPtStruct * MvarZeroSolverEmbedPt(const MvarPtStruct *Pt,
			   		    MvarMVDirType NewCoord,
					    CagdRType NewCoordVal);
static MvarPolylineStruct* MvarZeroSolverEmbedPoly(
					        const MvarPolylineStruct *Poly,
					        MvarMVDirType Dir,
					        CagdRType Val);
static MvarNormalConeStruct *MvarZeroBasisNormalCone(int Dim,
						     MvarMVDirType Dir);
static MvarZeroTJunctionStruct *MvarZeroSplitBoundaryCrv(
					       MvarZeroPrblmStruct *Problem,
					       MvarMVDirType Dir,
					       CagdRType SplitVal,
					       const MvarPtStruct *SplitEndPts,
					       MvarPolylineStruct **SubPoly1,
					       MvarPolylineStruct **SubPoly2, 
					       CagdBType *Success);
static MvarPtStruct *MvarZeroSolverMeanPt(const MvarPtStruct *Pt1,
					  const MvarPtStruct *Pt2,
					  CagdRType MeanVal,
					  MvarMVDirType Dir);
static void MvarZeroTagCornerPts(MvarPolylineStruct *PolylineList, 
				 CagdRType *MinDmn, 
				 CagdRType *MaxDmn);
static MvarPtStruct *MvarZeroFirstLastPts(
                              const MvarPolylineStruct *PolyList);
static const MvarPtStruct *MvarZeroClosestPt(const MvarPtStruct *PtList,
				       const MvarPtStruct *Pt, 
				       CagdRType *MinSqDist);

IRIT_STATIC_DATA const MvarZeroSolverCallBackFcnStruct 
    ZrMv2DCB = {
        MvarZeroSolveBoundary2DMVs,
        MvarZeroSolverNoZeroTest2DMVs,
        MvarZeroSolverTopoTest2DMVs,
	MvarZeroSubdivProblem2DMVs,
	MvarZeroSolverNumericStep2DMVs,
	MvarZeroSolverSinglrSol2DMVs,
	MvarZeroSolverUniteSolutions2DMVs, 
	MvarZeroSolverOrganizeSol2DMVs,
	MvarZeroUpdateProblemDmnMVs,
	MvarZeroFirstSmoothUpdatesMVs,
	NULL
    };

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the callback functions for the MVs representation, 2D solution case.M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*		                                                             *
* RETURN VALUE: 							     M
*   void								     M
*		                                                             *
* KEYWORDS:								     M
*   MvarZeroSolverSetCallbackFcns2DMVs					     M
*****************************************************************************/
void MvarZeroSolverSetCallbackFcns2DMVs(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    IRIT_GEN_COPY(&PrblmIntrnl -> CallbackFunctions, &ZrMv2DCB,
		  sizeof(MvarZeroSolverCallBackFcnStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Testing if we can rule out the possibility of a solution for the problem *
* in its domain, by inspecting constraint number Ind. Callback function      *
* for the MVs case, 2D solution. 					     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*   Ind:	The index of the MV to be tested.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if domain can be purged (no zero can exist), FALSE      *
*		otherwise.						     *
*****************************************************************************/
static CagdBType MvarZeroSolverNoZeroTest2DMVs(MvarZeroPrblmStruct *Problem,
					       int Ind)
{
    return MvarZeroMVConstraintFail(Problem -> U.MVs[Ind],
				    Problem -> Constraints[Ind]);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Solving the MVs representation problem on the boundary of the domain. In *
* the 2D case, it means we are solving 2 * Dim under constrained problems,   *
* with one degree of freedom, namely a 1D solution curve, and chain them     *
* together for correct connectivity, up to the prescribed tolerance.         *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroSolveBoundary2DMVs(MvarZeroPrblmStruct *Problem)
{
    CagdBType 
	FreeSingleFaceSol = FALSE;
    int i, j, k, 
	Dim = MvarZeroSolverGetDmnDim(Problem),
	NumOfConstraints = Problem -> NumOfConstraints;
    CagdRType CurBoundaryVal,
	SubdivTol = Problem -> SubdivTol,
	NumericTol = Problem -> NumericTol,
	StepTol = Problem -> StepTol;
    MvarPolylineStruct *PolyIter;
    MvarMVStruct
	**MVs = Problem -> U.MVs,
	**BoundaryMVs = (MvarMVStruct **) 
			 IritMalloc(sizeof(MvarMVStruct *) * NumOfConstraints);
    MvarZeroPrblmStruct *SingleFaceProblem;
    MvarZeroSolutionStruct *SingleFaceSolution,
	*FullBoundarySolution = NULL;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    /* For each of the 2 * Dim faces of dimension Dim - 1, extract the       */
    /* problem with the univariate solution and solve it.		     */
    for (i = 0; i < 2 * Dim; i++) {
	CagdBType 
	    NewDir = !(i % 2);
	MvarMVDirType 
	    Dir = i / 2; /* Division of integers giving an integer. */

	/* Restrict all MVs to this hyper-plane, create new problem & solve. */
	CurBoundaryVal = NewDir ? PrblmIntrnl -> MVMinDmn[Dir] 
				: PrblmIntrnl -> MVMaxDmn[Dir];
	for (j = 0; j < NumOfConstraints; j++)
	    BoundaryMVs[j] = MvarMVFromMV(MVs[j], CurBoundaryVal, Dir);
	SingleFaceProblem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *)
								   BoundaryMVs,
						   NULL, NumOfConstraints,
						   Problem -> Constraints,
						   SubdivTol, NumericTol, 
						   StepTol, FALSE);
	SingleFaceSolution = MvarZeroSolver(SingleFaceProblem);
	MvarZeroSolverPrblmFree(SingleFaceProblem);

	/* Embed the univariate solution correctly in R^Dim, as it is now in */
	/* R^(Dim - 1). Add the required coordinate of the face just solved. */
	if (SingleFaceSolution != NULL) {
	    MvarPolylineStruct *NewPoly,
		*EmbeddedPolyList = NULL;

	    PolyIter = SingleFaceSolution -> U.Pl;
	    while (PolyIter != NULL) {
		NewPoly = MvarZeroSolverEmbedPoly(PolyIter, Dir, 
					          CurBoundaryVal);
		IRIT_LIST_PUSH(NewPoly, EmbeddedPolyList);
		PolyIter = PolyIter -> Pnext;
	    }
	    MvarPolylineFreeList(SingleFaceSolution -> U.Pl);
	    SingleFaceSolution -> U.Pl = EmbeddedPolyList;
	    if (SingleFaceSolution != NULL) {
		MvarZeroTagCornerPts(SingleFaceSolution -> U.Pl, 
		    PrblmIntrnl -> MVMinDmn, PrblmIntrnl -> MVMaxDmn);
	    }
	    if (FullBoundarySolution == NULL) { /* The first. */
		FullBoundarySolution = SingleFaceSolution;
		FreeSingleFaceSol = FALSE;
	    }
	    else { /* Not the first solution we concatenate. */
		CagdListAppend(FullBoundarySolution -> U.Pl, 
			       SingleFaceSolution -> U.Pl);
		/* The Solution struct, not the poly, will be deallocated. */
		FreeSingleFaceSol = TRUE; 
	    }
	}
	/* Free the allocated info of this face, and proceed to the next. */
	if (FreeSingleFaceSol)
	    /* The polylines were merged in places, free only the struct. */
	    MvarZeroSolverSolutionFree(SingleFaceSolution, FALSE);
	for (k = 0; k < NumOfConstraints; k++) {
	    MvarMVFree(BoundaryMVs[k]);
	}
    }

    /* Arrange the polylines with correct connectivity. This step must       */
    /* result in a set of closed loops for regular problems. If chaining     */
    /* fails - the Single Loop Boundary Test (SLBT) shall fail and further   */
    /* subdivision shall occur, even a single loop should've been detected.  */
    if (FullBoundarySolution != NULL) {
	FullBoundarySolution -> U.Pl = MvarPolyMergePolylines(
	                                          FullBoundarySolution -> U.Pl, 
	                                          MVAR_ZRMV2D_POLY_CHAIN_TOL);
	PrblmIntrnl -> NumOfBoundarySolutions = 
				  CagdListLength(FullBoundarySolution -> U.Pl);
    }
    else { /* Found no boundary solutions. */
	PrblmIntrnl -> NumOfBoundarySolutions = 0;
    }

    IritFree(BoundaryMVs);

    /* Update the Problem structure with the new information. */
    assert(PrblmIntrnl -> SolutionsOnBoundary == NULL);
    PrblmIntrnl -> SolutionsOnBoundary = FullBoundarySolution;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Topological guarantee in the MVs, 2D solution case.	Test for sufficient  *
* conditions that guarantee the solution set is homeomorphic to a closed,    *
* two dimensional disc, or an empty solution.	                 	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if the topology can be determined - either a closed     *
*               disc or an empty solution. FALSE if cannot guarantee.        *
*****************************************************************************/
static CagdBType MvarZeroSolverTopoTest2DMVs(MvarZeroPrblmStruct *Problem)
{
    CagdBType OneToOneProj;
    int PtsNum;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    switch (PrblmIntrnl -> NumOfBoundarySolutions) {
	case 0:
	    /* Test if the projection is one-to-one on some 2D plane and    */
	    /* update One2OneProjOnPlane2D to TRUE if so, FALSE otherwise.  */
	    if ((OneToOneProj = MvarZeroSolverOne2OneProjMVs(Problem))
		                                                  != FALSE) {
		/* Topology is guaranteed - empty. */
		PrblmIntrnl -> PurgeByNumericStep = TRUE;

#		if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
		    printf("\nThis domain is purged by Topo' test.\n");
#		endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */

		return TRUE;
	    }
	    else {
		/* Cannot be sure - there may be closed surfaces. */

#		if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
	            printf("\nThis domain failed 1-1 proj' test.\n");
#		endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */

		return FALSE;
	    }
	    break;
	case 1:
	    PtsNum = CagdListLength(
		            PrblmIntrnl -> SolutionsOnBoundary -> U.Pl -> Pl);
	    if (PtsNum < 4) { /* Suspected as a tangency point/location. */
	        /* Test if the projection is one-to-one on some 2D plane and */
	        /* set One2OneProjOnPlane2D to TRUE if so, FALSE otherwise.  */
		if ((OneToOneProj = MvarZeroSolverOne2OneProjMVs(Problem)) 
								!= FALSE) {
		    /* A solution cannot exist in the interior. */
		    PrblmIntrnl -> PurgeByNumericStep = TRUE;

#		    if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
		        printf("\nThis domain is purged by Topo' test.\n");
#		    endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */

		    return TRUE;
		}
		else {
		    /* Cannot be sure - there may be closed surfaces. */

#		    if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
		        printf("\nThis domain failed 1-1 proj' test.\n");
#		    endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */
		    return FALSE;
		}
	    }
	    else { /* A boundary polyline with at least 3 distinct pts. */
		/* Test if the projection is one-to-one on some 2D plane and */
	        /* set One2OneProjOnPlane2D to TRUE if so, FALSE otherwise.  */
		MvarPtStruct
		    *FirstPt = PrblmIntrnl -> SolutionsOnBoundary -> U.Pl -> Pl,
		    *LastPt = CagdListLast(FirstPt);
		if (MvarPtCmpTwoPoints(FirstPt, LastPt, 
		                       MVAR_ZRMV2D_POLY_POLY_CLOSED_TOL) != 0) {

#		    if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
		        printf("\nThis domain has a non loop boundary curve.\n");
#		    endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */

		    return FALSE;
		}
		if ((OneToOneProj = MvarZeroSolverOne2OneProjMVs(Problem))
								     != FALSE) {
		    /* Success: homeo' to a 2-D disc. */
#		    ifdef MVAR_DBG_2D_TOPO_CONCLUDE
#			ifdef DEBUG
			    printf("\nThis domain is topo' guaranteed.\n");
#			endif /* DEBUG */
		        (*Problem -> _Internal -> DBGCntGuaranteedDmns)++;
#		    endif /* MVAR_DBG_2D_TOPO_CONCLUDE */
		    return TRUE;
		}
		else {
#		    if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
		        printf("\nThis domain failed 1-1 proj' test.\n");
#		    endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */
		    return FALSE;
		}
	    }
	    break;
	default:
#	    if defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE)
		printf("\nThis domain has more than 1 boundary curve.\n");
#	    endif /* defined(DEBUG) && defined(MVAR_DBG_2D_TOPO_CONCLUDE) */
	    return FALSE;
    }  
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   The numeric improvement for the MVs representation problem, 2D case. We  *
* only make some preparations here such as saving the local MVs problem as   *
* an attribute, but no triangulation.					     *
*   The actual triangulation occurs after returning the union of all loops,  *
* and only then will the solution representation switch to triangles. See    *
* ZrMV2DTs.c, organize solution function.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The updated polyline solution.		     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep2DMVs(
					         MvarZeroPrblmStruct *Problem)
{
    int k;
    MvarPolylineStruct *Poly;
    MvarMVStruct **MVsCpy;
    MvarZeroSubDmnInfoStruct *TriangulationInfo;
    MvarZeroSolutionStruct *Sol;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    if (PrblmIntrnl -> PurgeByNumericStep)
	return NULL;
    if (PrblmIntrnl -> NumOfBoundarySolutions != 1) {
	/* Something went wrong. Shouldn't be here. */
	assert(0);
	return NULL;
    }

    Poly = PrblmIntrnl -> SolutionsOnBoundary -> U.Pl;
    if (Poly == NULL)
	return NULL;
    PrblmIntrnl -> SolutionsOnBoundary -> U.Pl = MvarPolylineCopy(Poly);

    /* If this is a loop intended to be triangulated - remember the local    */
    /* info' for later.							     */
    MVsCpy = (MvarMVStruct **) IritMalloc(
	                Problem -> NumOfConstraints * sizeof(MvarMVStruct *));
    for (k = 0; k < Problem -> NumOfConstraints; k++) {
	MVsCpy[k] = MvarMVCopy(Problem -> U.MVs[k]);
    }
    TriangulationInfo = MvarSubDmnInfoStructNew(
					    MVsCpy, 
					    PrblmIntrnl -> One2OneProjDirs[0],
					    PrblmIntrnl -> One2OneProjDirs[1]);
    AttrSetRefPtrAttrib(&Poly -> Attr, "_TrInfo", TriangulationInfo);
    Sol = MvarZeroSolverSolutionNew(NULL, Poly, NULL, 
				    MVAR_ZER_SLVR_SOLUTION_POLYLINE); 
    return Sol;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returning a solution in the case of reaching subdivision tolerance. In   *
* the 2D case we simply return nothing. 				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution.				     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol2DMVs(
					         MvarZeroPrblmStruct *Problem)
{
#ifdef MVAR_DBG_2D_TOPO_CONCLUDE
       (*Problem -> _Internal -> DBGCntSinglrDmns)++;
#ifdef DEBUG
       printf("\nThis domain reached subdiv tol!\n\n");
#endif /* DEBUG */
#endif /* MVAR_DBG_2D_TOPO_CONCLUDE */

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Unites two neighboring 2D solutions of a zero finding problem, MVs	     *
* representation, back to a new solution of the subdivided domain:           *
* Lists of closed loops are simply appended.                                 *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol1:	The first solution.					     *
*   Sol2:	The second solution.					     *
*   Dir:	The direction of uniting the solution. Not used.    	     *
*   Param:	The parameter value at which the subdivision occurred, along *
*		direction Dir. Not used.				     *
*   Problem:	The zero finding problems solved by the united solution.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The new solution, union of Sol1 and Sol2.     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions2DMVs(
						 MvarZeroSolutionStruct *Sol1,
						 MvarZeroSolutionStruct *Sol2, 
					         MvarMVDirType Dir,
					         CagdRType Param,
						 MvarZeroPrblmStruct *Problem)
{
    MvarPolylineStruct *NewPlList;
    MvarZeroTJunctionStruct *NewTJList;
    MvarZeroSolutionStruct *Solution;

#   ifdef MVAR_DBG_2D_TOPO_CONCLUDE
        if (Problem -> _Internal -> SubdivDepth == 0) {
	    printf("%d domains reached subdiv tol.\n", 
		   *Problem -> _Internal -> DBGCntSinglrDmns);
	    printf("%d domains are topo' guaranteed.\n", 
		   *Problem -> _Internal -> DBGCntGuaranteedDmns);
	}
#   endif /* MVAR_DBG_2D_TOPO_CONCLUDE */

    if (Sol1 == NULL) {
	if (Sol2 != NULL) {
	    /* Potential T-Junctions handed over from the problem to the solution: */
	    Sol2 -> TJList = CagdListAppend(Sol2 -> TJList,
		                            Problem -> _Internal -> TJList);
	    Problem -> _Internal -> TJList = NULL;
	}
	return Sol2;
    }
    if (Sol2 == NULL) {
	if (Sol1 != NULL) {
	    Sol1 -> TJList = CagdListAppend(Sol1 -> TJList,
		                            Problem -> _Internal -> TJList);
	    Problem -> _Internal -> TJList = NULL;
	}
	return Sol1;
    }
    NewPlList = CagdListAppend(Sol1 -> U.Pl, Sol2 -> U.Pl);
    NewTJList = CagdListAppend(Sol1 -> TJList, Sol2 -> TJList);
    NewTJList = CagdListAppend(NewTJList, Problem -> _Internal -> TJList);
    Problem -> _Internal -> TJList = NULL;

    Solution = MvarZeroSolverSolutionNew(NULL, NewPlList, NULL, 
				         MVAR_ZER_SLVR_SOLUTION_POLYLINE);
    /* Hand over the T-Junctions list, in place, to the united solution.     */
    /* Set the lists in Sol1, Sol2 to NULL to avoid the deallocation.        */
    Solution -> TJList = NewTJList;
    Sol1 -> TJList = NULL;
    Sol2 -> TJList = NULL;
    MvarZeroSolverSolutionFree(Sol1, FALSE);
    MvarZeroSolverSolutionFree(Sol2, FALSE);
    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Tests for a guarantee that the solution set can be projected into a      *
* a 2D plane in a one-to-one manner, by checking the orthogonal projection   *
* on each of the optional 2D coordinate planes, using the cone test. The     *
* Problem structure is updated with the result, TRUE if can guarantee, FALSE *
* otherwise.								     *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if a one-to-one projection can be guaranteed, FALSE      *
*	       otherwise.				                     *
*****************************************************************************/
static CagdBType MvarZeroSolverOne2OneProjMVs(MvarZeroPrblmStruct *Problem)
{
    int i, j, 
	Dim = MvarZeroSolverGetDmnDim(Problem);
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    MvarNormalConeStruct *ConesList, *NextCone, *ConeIter, 
	*BasisConeI, *BasisConeJ;

    /* Build all normal cones of all (Dim-2) MVs. */
    ConesList = MVarMVNormalCone(Problem -> U.MVs[0]);
    if (!ConesList) {
	return PrblmIntrnl -> One2OneProjOnPlane2D = FALSE;
    }
    ConeIter = ConesList;
    for (i = 1; i < Dim - 2; i++) {
	NextCone = MVarMVNormalCone(Problem -> U.MVs[i]);
	if (!NextCone) {
	    MvarNormalConeFreeList(ConesList);
	    return PrblmIntrnl -> One2OneProjOnPlane2D = FALSE;
	}
	ConeIter -> Pnext = NextCone;
	ConeIter = ConeIter -> Pnext;
    }

    /* The first Dim - 2 cones are ready. Complete the optional last two as */
    /* degenerated normal cones of opening angle of zero.                   */
    for (i = 0; i < Dim; i++) {
	BasisConeI = MvarZeroBasisNormalCone(Dim, i);
	ConeIter -> Pnext = BasisConeI;

	for (j = i + 1; j < Dim; j++) {
	    BasisConeJ = MvarZeroBasisNormalCone(Dim, j);
	    BasisConeI -> Pnext = BasisConeJ;
	    if (MvarConesOverlapAux(ConesList)) {
		MvarNormalConeFree(BasisConeJ);
		BasisConeI -> Pnext = NULL;
	    }
	    else { /* Success. Update and abort the process. */
		PrblmIntrnl -> One2OneProjOnPlane2D = TRUE;
		PrblmIntrnl -> One2OneProjDirs[0] = i;
		PrblmIntrnl -> One2OneProjDirs[1] = j;
		MvarNormalConeFree(BasisConeJ);
		BasisConeI -> Pnext = NULL;
		break;
	    }
	}

	MvarNormalConeFree(BasisConeI);
	ConeIter -> Pnext = NULL;

	if (PrblmIntrnl -> One2OneProjOnPlane2D == TRUE) {
	    break;
	}
    }
    MvarNormalConeFreeList(ConesList);
    return PrblmIntrnl -> One2OneProjOnPlane2D;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates a new point by adding a coordinate to the input point, at the    *
* prescribed value and direction.					     *
*									     *
* PARAMETERS:								     *
*   Pt:		    The point in R^(Dim - 1).				     *
*   NewCoord:	    The new coordinate.					     *
*   NewCoordVal:    The real value at the new, added coordinate.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: The new point.					     *
*****************************************************************************/
static MvarPtStruct *MvarZeroSolverEmbedPt(const MvarPtStruct *Pt,
			   		   MvarMVDirType NewCoord,
					   CagdRType NewCoordVal)
{
    int OldDim = Pt -> Dim,
	NewDim = OldDim + 1;
    MvarPtStruct 
	*NewPt = MvarPtNew(NewDim);

    assert(NewCoord >= 0 && NewCoord < NewDim);
    CAGD_GEN_COPY(NewPt -> Pt, Pt -> Pt, sizeof(CagdRType) * NewCoord);
    NewPt -> Pt[NewCoord] = NewCoordVal;
    CAGD_GEN_COPY(&NewPt -> Pt[NewCoord + 1], &Pt -> Pt[NewCoord],
	          sizeof(CagdRType) * (NewDim - NewCoord));
    NewPt -> Attr = IP_ATTR_COPY_ATTRS(Pt -> Attr);
    NewPt -> Pnext = Pt -> Pnext;

    return NewPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the output type of the bivariate solver: either a collection of     M
* loops (polylines) or a collection of triangles.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   IsPolyLines2DSolution:   New setting for the output type.                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting for the output type.	         	             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver, MvarZeroSolverOne2OneProjMVs			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZeros2DPolylines                                                  M
*****************************************************************************/
int MvarMVsZeros2DPolylines(int IsPolyLines2DSolution)
{
    int	OldVal = _MVGlblZeroOutputTypeLoops;

    _MVGlblZeroOutputTypeLoops = IsPolyLines2DSolution;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates a new polyline by adding a coordinate to all points of the       *
* input polyline, according to the prescribed direction and value.	     *
*   The order is reversed, hence if the order is important, the output is to *
* reversed again (not done here!).					     *
*									     *
* PARAMETERS:								     *
*   Poly:   The polyline in R^(Dim - 1).		                     *
*   Dir:    The new coordinate.					             *
*   Val:    The real value at the new, added coordinate.	             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: The new polyline.				     *
*****************************************************************************/
static MvarPolylineStruct *MvarZeroSolverEmbedPoly(
					        const MvarPolylineStruct *Poly,
					        MvarMVDirType Dir,
					        CagdRType Val)
{
    MvarPtStruct *NewPt,
	*TmpPtList = NULL,
	*PtIter = Poly -> Pl;
    MvarPolylineStruct *NewPoly;

    while (PtIter != NULL) {
	NewPt = MvarZeroSolverEmbedPt(PtIter, Dir, Val);
	IRIT_LIST_PUSH(NewPt, TmpPtList);
	PtIter = PtIter -> Pnext;
    }
    NewPoly = MvarPolylineNew(TmpPtList);
    return NewPoly;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Construct a normal cone structure of prescribed dimension, with the      *
* the standard basis element of direction Dir as the cone axis, and span     *
* angle zero.								     *
*									     *
* PARAMETERS:								     *
*   Dim: The dimension of the cone.					     *
*   Dir: The selected direction for the cone axis.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarNormalConeStruct *: The new cone.				     *
*****************************************************************************/
static MvarNormalConeStruct *MvarZeroBasisNormalCone(int Dim,
						     MvarMVDirType Dir)
{
    int i;
    MvarNormalConeStruct
	*Cone = MvarNormalConeNew(Dim);

    Cone -> ConeAngleCosine = 1.0;
    for (i = 0; i < Dim; i++) {
	Cone -> ConeAxis -> Vec[i] = i == Dir ? 1.0 : 0.0;
    }
    return Cone;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Subdivision of the zero finding problem.  Given a problem returns a      *
* vector of two subdivided problems, allocated dynamically.  The subdivided  *
* problems can be NULL.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroPrblmStruct **:  The array of two sub-problems.		     *
*****************************************************************************/
static MvarZeroPrblmStruct **MvarZeroSubdivProblem2DMVs(
                                                  MvarZeroPrblmStruct *Problem)
{
    CagdBType SubdivSuccess;
    int i, j,
	NumOfConstraints = Problem -> NumOfConstraints;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    MvarMVDirType 
	Dir = PrblmIntrnl -> MaxSideDir;
    CagdRType t, MaxSide, MaxSideMid, Perturb,
	SubdivTol = Problem -> SubdivTol,
	NumericTol = Problem -> NumericTol,
	Min = PrblmIntrnl -> MVMinDmn[Dir],
	Max = PrblmIntrnl -> MVMaxDmn[Dir];
    MvarPtStruct 
	*SplitPlanePolyStartEndPts = NULL;
    MvarPolylineStruct *SplitBoundaryPolys[2], *SplitPlanePoly,
        *SplitPlanePolyCpy;
    MvarMVStruct **MVs1, **MVs2, 
	**MVs = Problem -> U.MVs,
	**SplitPlaneMVs = (MvarMVStruct **) 
		         IritMalloc(sizeof(MvarMVStruct *) * NumOfConstraints);
    MvarZeroPrblmStruct *SplitPlaneProblem,
	**SubProblems = (MvarZeroPrblmStruct **)
		                 IritMalloc(2 * sizeof(MvarZeroPrblmStruct *));
    MvarZeroSolutionStruct *SplitPlaneSolution, *BoundarySol1, *BoundarySol2;
    			
    SubProblems[0] = SubProblems[1] = NULL;

    MaxSide = PrblmIntrnl -> MaxSide;
    MaxSideMid = PrblmIntrnl -> MVMinDmn[Dir] + 0.5 * MaxSide;

    /* Take the default value (small perturbation of the middle). */
    Perturb = MVAR_MVZR1D_SUBD_PRTRB * (1 + Dir * 0.1);
    assert(Perturb < 0.1);
    t = (0.5 + Perturb) * Max + (0.5 - Perturb) * Min;
    
    /* Find the solution on the splitting hyper-plane: */
    for (j = 0; j < NumOfConstraints; j++)
	SplitPlaneMVs[j] = MvarMVFromMV(MVs[j], t, Dir);

    SplitPlaneProblem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *)
								SplitPlaneMVs,
					       NULL, NumOfConstraints,
					       Problem -> Constraints, SubdivTol,
					       NumericTol, 
					       Problem -> StepTol, FALSE);
    SplitPlaneSolution = MvarZeroSolver(SplitPlaneProblem);

    SplitPlanePoly = SplitPlanePolyCpy = NULL;
    MvarZeroSolverPrblmFree(SplitPlaneProblem);

    /* Embed in one dim' higher domain, and make a copy. Each of these is to */
    /* be concatenated with the boundary solutions of the sub-domain.        */
    if (SplitPlaneSolution != NULL) {
	MvarPolylineStruct *NewPoly, 
	    *EmbeddedPolyList = NULL, 
	    *PolyIter = SplitPlaneSolution -> U.Pl;

#	ifdef MVAR_MVZR2D_DEBUG_SPLIT_POLYS
	    printf("\nThe solution on the split plane is:\n");
#	endif /* MVAR_MVZR2D_DEBUG_SPLIT_POLYS */

	while (PolyIter != NULL) {
	    NewPoly = MvarZeroSolverEmbedPoly(PolyIter, Dir, t);
	    IRIT_LIST_PUSH(NewPoly, EmbeddedPolyList);
	    PolyIter = PolyIter -> Pnext;
	}
	MvarPolylineFreeList(SplitPlaneSolution -> U.Pl);

	SplitPlanePoly = SplitPlaneSolution -> U.Pl = EmbeddedPolyList;
	MvarZeroTagCornerPts(SplitPlanePoly, PrblmIntrnl -> MVMinDmn, PrblmIntrnl -> MVMaxDmn);
	SplitPlanePolyCpy = MvarPolylineCopyList(SplitPlanePoly);
    }

    /* Split the boundary solution of the original problem. */
    SplitBoundaryPolys[0] = SplitBoundaryPolys[1] = NULL;
    if (SplitPlaneSolution != NULL) {
	SplitPlanePolyStartEndPts = 
	    MvarZeroFirstLastPts(SplitPlaneSolution -> U.Pl);
    }
    if (PrblmIntrnl -> SolutionsOnBoundary != NULL || 
	SplitPlaneSolution != NULL) {
	if (PrblmIntrnl -> SolutionsOnBoundary != NULL) {
	    SubdivSuccess = TRUE;
	    PrblmIntrnl -> TJList = MvarZeroSplitBoundaryCrv(
		                                     Problem, Dir, t, 
		                                     SplitPlanePolyStartEndPts,
				                     &SplitBoundaryPolys[0],
			                             &SplitBoundaryPolys[1], 
						     &SubdivSuccess);
	    if (!SubdivSuccess) {
		/* The split plane solution could not be resolved in a       */
		/* consistent manner w.r.t the boundary solutions split. We  */
		/* cannot recover from this, even by further subdivision.    */
		SubProblems[0] = SubProblems[1] = NULL;
		for (j = 0; j < NumOfConstraints; j++)
		    MvarMVFree(SplitPlaneMVs[j]);
		IritFree(SplitPlaneMVs);
		MvarPtFreeList(SplitPlanePolyStartEndPts);
		MvarZeroSolverSolutionFree(SplitPlaneSolution, FALSE);
		return SubProblems;
	    }
	}
	
	/* Finally - chain together with the split plane solution. */
	SplitBoundaryPolys[0] = CagdListAppend(SplitBoundaryPolys[0],
					       SplitPlanePoly);
	SplitBoundaryPolys[0] = MvarPolyMergePolylines(
	                                           SplitBoundaryPolys[0], 
	                                           MVAR_ZRMV2D_POLY_CHAIN_TOL);

	BoundarySol1 = MvarZeroSolverSolutionNew(
					    NULL, SplitBoundaryPolys[0], NULL,
					    MVAR_ZER_SLVR_SOLUTION_POLYLINE);
	SplitBoundaryPolys[1] = CagdListAppend(SplitBoundaryPolys[1], 
					       SplitPlanePolyCpy);
	SplitBoundaryPolys[1] = MvarPolyMergePolylines(
	                                             SplitBoundaryPolys[1], 
	                                             MVAR_ZRMV2D_POLY_CHAIN_TOL);

	BoundarySol2 = MvarZeroSolverSolutionNew(
					    NULL, SplitBoundaryPolys[1], NULL,
					    MVAR_ZER_SLVR_SOLUTION_POLYLINE);
    }
    else { /* No boundary solution and no split plane solution. */
	BoundarySol2 = BoundarySol1 = NULL;
    }    

    /* All is ready, subdivide and create the sub-problems. */
    MVs1 = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) *
					                    NumOfConstraints);
    MVs2 = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) *
					                    NumOfConstraints);
    for (i = 0; i < NumOfConstraints; i++) {
	MVs1[i] = MvarMVSubdivAtParam(MVs[i], t, Dir);
	MVs2[i] = MVs1[i] -> Pnext;
	MVs1[i] -> Pnext = NULL;
    }

    /* Remember the subdiv' info in the parent problem: */
    PrblmIntrnl -> SubdivDir = Dir;
    PrblmIntrnl -> ParamOfSubdiv = t;
    SubProblems[0] = MvarZeroSolverSubProblem(Problem, MVs1, NULL, 
					      BoundarySol1);
    SubProblems[1] = MvarZeroSolverSubProblem(Problem, MVs2, NULL, 
					      BoundarySol2);

    for (j = 0; j < NumOfConstraints; j++)
	MvarMVFree(SplitPlaneMVs[j]);
    IritFree(SplitPlaneMVs);
    MvarPtFreeList(SplitPlanePolyStartEndPts);
    /* Free solution wrapper, keeping the polyline (that are used above). */
    MvarZeroSolverSolutionFree(SplitPlaneSolution, FALSE);

    return SubProblems;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a direction and a split value along it, splits the polylines       *
* into two sublists, according to coordinate value at direction Dir (greater *
* or smaller than the split value). If a polyline crosses the split plane,   *
* the crossing segment is split at the plane, adding (if required) the       *
* corresponding end point and start point to the new, sub-polylines.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:        The problem, boundary solution curve of which we split.  *
*   Dir:	    The given direction.				     *
*   SplitVal:	    The value along direction Dir according to which the     *
*		    split is done.					     *
*   SplitEndPts:    The start and end points of the univariate solution      *
*		    curves contained in the split plane. Used for the        *
*		    refinement of the curve split location to the solution.  *
*   SubPoly1:	    The output sub-list to hold the polylines with points    *
*		    having a smaller value at coordinate Dir than the split  *
*		    value.						     *
*   SubPoly2:	    The output sub-list to hold the polylines with points    *
*		    having a greater value at coordinate Dir than the split  *
*		    value.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroTJunctionStruct *:	The list of new potential T-Junction events, *
*				due to the curve split operation.            *
*****************************************************************************/
static MvarZeroTJunctionStruct *MvarZeroSplitBoundaryCrv(
                                               MvarZeroPrblmStruct *Problem,
		                               MvarMVDirType Dir,
                                               CagdRType SplitVal,
				               const MvarPtStruct *SplitEndPts,
                                               MvarPolylineStruct **SubPoly1,
                                               MvarPolylineStruct **SubPoly2, 
					       CagdBType *Success)
{
    CagdBType CurLessThanVal, CurMoreThanVal, CurInPlane, IsNewPolyLessThanVal;
    int Dim = Problem -> U.MVs[0] -> Dim;
    CagdRType Val, SqrDist;
    MvarPtStruct *NewPt, *MeanPt, *MeanPtCpy,
	*CurPtList = NULL;
    const MvarPtStruct *ClosestSplitPlaneEndPt;
    MvarPolylineStruct *NewPoly, *PolyIter;
    MvarZeroTJunctionStruct *TJ, 
	*TJList = NULL;
    MvarZeroPrblmIntrnlStruct 
	*PrblmIntrnl = Problem -> _Internal;

    PolyIter = PrblmIntrnl -> SolutionsOnBoundary -> U.Pl;
    *Success = TRUE;
    
    while (PolyIter != NULL) {
	IsNewPolyLessThanVal = PolyIter -> Pl -> Pt[Dir] <= SplitVal;
	while (PolyIter -> Pl) {
	    IRIT_LIST_POP(NewPt, PolyIter -> Pl);
	    Val = NewPt -> Pt[Dir];

	    /* Update the status to the current point: */
	    CurLessThanVal = Val < SplitVal - MVAR_ZRMV2D_PT_IN_FACE_TOL;
	    CurMoreThanVal = Val > SplitVal + MVAR_ZRMV2D_PT_IN_FACE_TOL;
	    CurInPlane = !CurLessThanVal && !CurMoreThanVal;

	    if (CurInPlane) {/* In the split plane. Not considered a cross: */
		/* Assign the exact value to avoid numeric trouble. */
		NewPt -> Pt[Dir] = SplitVal;
		/* Add to the current list and don't change the list state. */
		IRIT_LIST_PUSH(NewPt, CurPtList);
	    }
	    else if (CurLessThanVal) {
		if (IsNewPolyLessThanVal) {
		    /* Did not cross the split plane. Add and continue. */
		    IRIT_LIST_PUSH(NewPt, CurPtList);
		}
		else { /* Crossing from "greater than" to "less than". */
		    MeanPt = MvarZeroSolverMeanPt(NewPt, CurPtList,
			                          SplitVal, Dir);
		    AttrSetIntAttrib(&MeanPt -> Attr, "CrvSplitPt", TRUE);

		    /* Updating the mean point to the solution: */
		    if (SplitEndPts == NULL) {
			#ifdef MVAR_DBG_SPLIT_PLANE
			printf("Warning: splitting boundary curve without split-plane solution.\n");
			#endif /* MVAR_DBG_SPLIT_PLANE */

			*SubPoly1 = *SubPoly2 = NULL;
			*Success = FALSE;
		    }
		    ClosestSplitPlaneEndPt = MvarZeroClosestPt(SplitEndPts, 
			                                       MeanPt, 
							       &SqrDist);
		    if (ClosestSplitPlaneEndPt != NULL) {
			assert(SqrDist < Problem -> SubdivTol);
			IRIT_GEN_COPY(MeanPt -> Pt,
				      ClosestSplitPlaneEndPt -> Pt, 
				      Dim * sizeof(CagdRType));
			/* Succeeded in updating it to the solution: */
			AttrFreeOneAttribute(&MeanPt -> Attr, "CrvSplitPt");
			AttrSetIntAttrib(&MeanPt -> Attr, "Corner", TRUE);
			/* A new T-Junction candidate, unless we're at the   */
			/* init domain.                                      */
			if (PrblmIntrnl -> SubdivDepth > 0 && 
			    !MvarZeroIsPtOnDmnBndry(PrblmIntrnl ->
						                 OrigMVMinDmn,
						    PrblmIntrnl ->
						                 OrigMVMaxDmn,
						    MeanPt)) {
			    TJ = MvarZeroTJNew(CurPtList, MeanPt, NewPt);
			    IRIT_LIST_PUSH(TJ, TJList);
			}
		    }
		    MeanPtCpy = MvarPtCopy(MeanPt);
		    
		    IRIT_LIST_PUSH(MeanPt, CurPtList);
		    NewPoly = MvarPolylineNew(CurPtList);
		    CurPtList = MeanPtCpy;
		    IRIT_LIST_PUSH(NewPt, CurPtList);
		    IRIT_LIST_PUSH(NewPoly, *SubPoly2);
		    IsNewPolyLessThanVal = TRUE;
		}
	    }
	    else if (CurMoreThanVal) {
		if (!IsNewPolyLessThanVal) {
		    /* Did not cross the split plane. */
		    IRIT_LIST_PUSH(NewPt, CurPtList);
		}
		else {
		    /* Crossing occurred from "less than" to "greater than". */
		    MeanPt = MvarZeroSolverMeanPt(CurPtList, NewPt,
						  SplitVal, Dir);
		    AttrSetIntAttrib(&MeanPt -> Attr, "CrvSplitPt", TRUE);

		    /* Updating the mean point to the solution: */
		    if (SplitEndPts == NULL) {
			#ifdef MVAR_DBG_SPLIT_PLANE
			printf("Warning: splitting boundary curve without split-plane solution.\n");
			#endif /* MVAR_DBG_SPLIT_PLANE */

			*SubPoly1 = *SubPoly2 = NULL;
			*Success = FALSE;
		    }
		    ClosestSplitPlaneEndPt = MvarZeroClosestPt(SplitEndPts, 
			                                       MeanPt, 
			                                       &SqrDist);
		    if (ClosestSplitPlaneEndPt != NULL) {
			assert(SqrDist < Problem -> SubdivTol);
			IRIT_GEN_COPY(MeanPt -> Pt,
				      ClosestSplitPlaneEndPt -> Pt, 
			              Dim * sizeof(CagdRType));
			/* Succeeded in updating it to the solution: */
			AttrFreeOneAttribute(&MeanPt -> Attr, "CrvSplitPt");
			AttrSetIntAttrib(&MeanPt -> Attr, "Corner", TRUE);
			/* A new T-Junction candidate, unless we're at the   */
			/* init domain.                                      */
			if (PrblmIntrnl -> SubdivDepth > 0 && 
			    !MvarZeroIsPtOnDmnBndry(PrblmIntrnl ->
						                 OrigMVMinDmn,
						    PrblmIntrnl ->
						                 OrigMVMaxDmn,
						    MeanPt)) {
			    TJ = MvarZeroTJNew(NewPt, MeanPt, CurPtList);
			    IRIT_LIST_PUSH(TJ, TJList);
			}
		    }
		    MeanPtCpy = MvarPtCopy(MeanPt);
		    
		    IRIT_LIST_PUSH(MeanPt, CurPtList);
		    NewPoly = MvarPolylineNew(CurPtList);
		    CurPtList = MeanPtCpy;
		    IRIT_LIST_PUSH(NewPt, CurPtList);
		    IRIT_LIST_PUSH(NewPoly, *SubPoly1);
		    IsNewPolyLessThanVal = FALSE;
		}
	    }
	}

	/* Poly is exhausted. Wrap up accordingly and proceed to the next. */
	NewPoly = MvarPolylineNew(CurPtList);
	CurPtList = NULL;
	if (IsNewPolyLessThanVal) { 
	    IRIT_LIST_PUSH(NewPoly, *SubPoly1);
	}
	else {
	    IRIT_LIST_PUSH(NewPoly, *SubPoly2);
	}
	PolyIter = PolyIter -> Pnext;
    }
    return TJList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create a point on the segment [Pt1,Pt2], when one of its coordinates is  *
* prescribed, assumed to be a mean value between the respective coordinates  *
* of Pt1 and Pt2.							     *
*									     *
* PARAMETERS:                                                                *
*   Pt1:	The point with the smaller value in direction Dir.	     *
*   Pt2:	The point with the greater value in direction Dir.	     *
*   MeanVal:	The required value for the new point in direction Dir.	     *
*   Dir:	The direction in which MeanVal is required to be assigned.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct:  The mean point, allocated dynamically.		     *
*****************************************************************************/
static MvarPtStruct *MvarZeroSolverMeanPt(const MvarPtStruct *Pt1,
					  const MvarPtStruct *Pt2,
					  CagdRType MeanVal,
					  MvarMVDirType Dir)
{
    int i,
	Dim = Pt1 -> Dim;
    CagdRType Weight1, Weight2;
    MvarPtStruct *MeanPt;

    assert((MeanVal >= Pt1 -> Pt[Dir] && MeanVal <= Pt2 -> Pt[Dir]) ||
	   (MeanVal <= Pt1 -> Pt[Dir] && MeanVal >= Pt2 -> Pt[Dir]));
    
    /* The closer point gets the larger weight, linearly. */
    Weight2 = (MeanVal - Pt1 -> Pt[Dir]) / (Pt2 -> Pt[Dir] - Pt1 -> Pt[Dir]);
    Weight1 = 1.0 - Weight2;

    MeanPt = MvarPtNew(Dim);
    for (i = 0; i < Dim; i++) {
	MeanPt -> Pt[i] = Weight1 * (Pt1 -> Pt[i]) + Weight2 * (Pt2 -> Pt[i]);
    }
    return MeanPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new multi-variate triangle.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Dim:      Number of dimensions of each vertex.		             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarTriangleStruct *:    An uninitialized triangle.                      M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarTriangleFree                                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarTriangleNew                                                          M
*****************************************************************************/
MvarTriangleStruct *MvarTriangleNew(int Dim)
{
    CagdRType *R;
    MvarTriangleStruct
	*Triangle = (MvarTriangleStruct *) 
                               IritMalloc(sizeof(MvarTriangleStruct) + 8 + 
					  sizeof(CagdRType) * Dim * 6);

    Triangle -> Pnext = NULL;
    Triangle -> Dim = Dim;

    R = (CagdRType *) ((((IritIntPtrSizeType) &Triangle[1]) + 7) & ~0x07);
    Triangle -> Vrtcs[0] = R;
    Triangle -> Vrtcs[1] = &R[Dim];
    Triangle -> Vrtcs[2] = &R[Dim * 2];

    Triangle -> Nrmls[0] = &R[Dim * 3];
    Triangle -> Nrmls[1] = &R[Dim * 4];
    Triangle -> Nrmls[2] = &R[Dim * 5];

    return Triangle;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deallocates and frees all slots of a multi-variate triangle structure.   M
*                                                                            *
* PARAMETERS:                                                                M
*   Tr:      Multivariate triangle to free.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarTriangleNew, MvarTriangleFreeList                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarTriangleFree                                                         M
*****************************************************************************/
void MvarTriangleFree(MvarTriangleStruct *Tr)
{
    IritFree(Tr);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deallocates and frees a list of triangle structures.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   TrList:      Multivariate triangles list to free.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarTriangleNew, MvarTriangleFree                                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarTriangleFreeList                                                     M
*****************************************************************************/
void MvarTriangleFreeList(MvarTriangleStruct *TrList)
{
    MvarTriangleStruct *Temp,
	*TrIter = TrList;

    while (TrIter != NULL) {
	Temp = TrIter -> Pnext;
	MvarTriangleFree(TrIter);
	TrIter = Temp;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Projecting a given list of polylines on any of its required coordinates, M
* as specified by the Coords vector of length projDim which is a subset of   M
* input dimension.							     M
*									     *
* PARAMETERS:                                                                M
*   PolyList:	 A list of polylines to project.			     M
*   Coords:      The required projection coordinates.  Ordered vector.       M
*   ProjDim:     The dimension of the required result.  Also the length of   M
*                Coords.					             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *:    The new list of polylines.			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverPolyProject                                                M
*****************************************************************************/
MvarPolylineStruct *MvarZeroSolverPolyProject(MvarPolylineStruct *PolyList,
					      int *Coords,
					      int ProjDim)
{
    int i, j;
    MvarPtStruct *PtIter,
	*NewPt = NULL,
	*PtList = NULL,
	*PrevProjPt = NULL;
    MvarPolylineStruct 
	*NewPoly = NULL,
	*ProjPolyList = NULL,
	*PrevProjPoly = NULL,
	*PolyIter = PolyList;

    while (PolyIter != NULL) {
	PtIter = PolyIter -> Pl;
	while (PtIter != NULL) {
	    NewPt = MvarPtNew(ProjDim);
	    j = 0;
	    for (i = 0; i < ProjDim; i++)
		NewPt -> Pt[i] = PtIter -> Pt[Coords[i]];
	    if (PtList == NULL) {
		PtList = NewPt;
	    }
	    else 
		PrevProjPt -> Pnext = NewPt;
	    PrevProjPt = NewPt;
	    PtIter = PtIter -> Pnext;
	}
	NewPt -> Pnext = NULL; /* Last point in current list. */
	NewPoly = MvarPolylineNew(PtList);
	PtList = NULL;
	if (ProjPolyList == NULL)
	    ProjPolyList = NewPoly;
	else
	    PrevProjPoly -> Pnext = NewPoly;
	PrevProjPoly = NewPoly;
	PolyIter = PolyIter -> Pnext;
    }
    NewPoly -> Pnext = NULL; /* Last projected polyline in list. */
    return ProjPolyList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a list of univariate solution curves that belong to a single       *
* hyper-plane (either a face of the domain, or a splitting hyper-plane),     *
* provide a "Corner" attribute to those start/end-points that lay on "edges" *
* i.e. faces of dimension n - 2. In other words: closed loops are untouched, *
* open curves have their entry and exit points tagged as corners.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   PolylineList:   The list of polyline solutions.			     *
*   MinDmn:	    The min' values of the domain in each direction.	     *
*   MaxDmn:	    The max' values of the domain in each direction.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroTagCornerPts(MvarPolylineStruct *PolylineList, 
				 CagdRType *MinDmn,
				 CagdRType *MaxDmn)
{
    MvarPtStruct *StartPt, *EndPt;
    MvarPolylineStruct
	*PolyIter = PolylineList;

    if (PolylineList == NULL)
	return;

    /* The start and end points are always corners later, unless they are    */
    /* equal, which means the polyline is a loop.                            */
    while (PolyIter != NULL) {
	StartPt = PolyIter -> Pl;
	EndPt = CagdListLast(PolyIter -> Pl);
	if (MvarPtCmpTwoPoints(StartPt, EndPt, MVAR_ZRMV2D_PTS_EQ_TOL)) {
	    /* Not a closed polyline, start and end must be corners. */
	    AttrSetIntAttrib(&StartPt -> Attr, "Corner", TRUE);
	    AttrSetIntAttrib(&EndPt -> Attr, "Corner", TRUE);
	}
	else {
	    assert(!MvarZeroIsPtOnEdge(StartPt, MinDmn, MaxDmn, 
		                       MVAR_ZRMV2D_CORNER_PT_TOL) && 
		   !MvarZeroIsPtOnEdge(EndPt, MinDmn, MaxDmn, 
		                       MVAR_ZRMV2D_CORNER_PT_TOL));
	}
	PolyIter = PolyIter -> Pnext;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Constructs a list of new point objects, consisting of copies of the      *
* start and end points of the given polylines.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   PolyList:   The list of polylines, to extract the start/end points from. *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: The list of start/end point copies.   		     *
*****************************************************************************/
static MvarPtStruct *MvarZeroFirstLastPts(const MvarPolylineStruct *PolyList)
{
    MvarPtStruct *NewPt,
	*FirstLastPts = NULL;
    const MvarPolylineStruct 
	*PolyIter = PolyList;

    while (PolyIter != NULL) {
	NewPt = MvarPtCopy(PolyIter -> Pl);
	IRIT_LIST_PUSH(NewPt, FirstLastPts);
	NewPt = MvarPtCopy(CagdListLast(PolyIter -> Pl));
	IRIT_LIST_PUSH(NewPt, FirstLastPts);
	PolyIter = PolyIter -> Pnext;
    }
    return CagdListReverse(FirstLastPts);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a fixed point and a list of points to check, returns the pointer   *
* to the point from the list, that is closest to the given fixed point.      *
*                                                                            *
* PARAMETERS:                                                                *
*   PtList:	The list of points, to extract the closest.	             *
*   Pt:		The constant point to check distance from.	             *
*   MinSqDist:	The minimal distance realized by the closest pt, squared.    *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: The closest point from the list.     		     *
*****************************************************************************/
static const MvarPtStruct *MvarZeroClosestPt(const MvarPtStruct *PtList,
					     const MvarPtStruct *Pt, 
					     CagdRType *MinSqDist)
{
    CagdRType TmpSqDist; 
    const MvarPtStruct 
	*PtIter = PtList,
	*ClosestPt = NULL;
	
    *MinSqDist = IRIT_INFNTY;
    while (PtIter != NULL) {
	if ((TmpSqDist = MvarPtDistSqrTwoPoints(Pt, PtIter)) < *MinSqDist) {
	    *MinSqDist = TmpSqDist;
	    ClosestPt = PtIter;
	}
	PtIter = PtIter -> Pnext;
    }
    return ClosestPt;
}
