/******************************************************************************
* ZrMV1D.c - tools to compute zero sets of multivariates when the problem has *
*	     MVs representation and the expected solution  set                *
*	     is a one-manifold: curve(s).				      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Michael Barton and Gershon Elber, Apr' 08, Yoni Mizrahi, Feb' 13.*
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"

#define MVAR_MVZR1D_SUBD_PRTRB	    1.301060e-2

IRIT_STATIC_DATA int
    GlblMvarInterMergeSingularPts = 1;

static CagdBType MvarZeroSolverNoZeroTest1DMVs(MvarZeroPrblmStruct *Problem,
					      int Ind);
static void MvarZeroSolveBoundary1DMVs(MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverTopoTest1DMVs(MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep1DMVs(
					       MvarZeroPrblmStruct *Problem);
static MvarZeroPrblmStruct **MvarZeroSubdivProblem1DMVs(MvarZeroPrblmStruct
							            *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol1DMVs(MvarZeroPrblmStruct
							             *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions1DMVs(
						MvarZeroSolutionStruct *Sol1,
						MvarZeroSolutionStruct *Sol2,
					        MvarMVDirType Dir,
					        CagdRType Param,
						MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol1DMVs
						  (MvarZeroSolutionStruct *Sol,
						   MvarZeroPrblmStruct 
					                             *Problem);


IRIT_STATIC_DATA const MvarZeroSolverCallBackFcnStruct 
    ZrMv1DCB = {
        MvarZeroSolveBoundary1DMVs,
        MvarZeroSolverNoZeroTest1DMVs,
        MvarZeroSolverTopoTest1DMVs,
	MvarZeroSubdivProblem1DMVs,
	MvarZeroSolverNumericStep1DMVs,
	MvarZeroSolverSinglrSol1DMVs,
	MvarZeroSolverUniteSolutions1DMVs,
	MvarZeroSolverOrganizeSol1DMVs,
	MvarZeroUpdateProblemDmnMVs,
	MvarZeroFirstSmoothUpdatesMVs,
	NULL
    };

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the callback functions for the MVs representation, 1D solution case.M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*		                                                             *
* RETURN VALUE: 							     M
*   void								     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroOrganizeMVs1DProblem	    	                             M
*		                                                             *
* KEYWORDS:								     M
*   MvarZeroSolverSetCallbackFcns1DMVs					     M
*****************************************************************************/
void MvarZeroSolverSetCallbackFcns1DMVs(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    IRIT_GEN_COPY(&PrblmIntrnl -> CallbackFunctions, &ZrMv1DCB,
		  sizeof(MvarZeroSolverCallBackFcnStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Some preliminary organization of the objects composing an MVs, 1D solution M
* zero finding problem. This routine is invoked by the problem structure     M
* construction routine, and should not be invoked when extracting a          M
* sub-problem from an existing one.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:		Array of multivariates.				     M
*   Constraints:	Equality or inequality constraints.		     M
*   NumOfMVs:		Number of constraints (may be updated).		     M
*									     *
* RETURN VALUE:                                                              M
*   MvarMVStruct **:  The updated array of MVs to be used.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverSetCallbackFcns1DMVs	    	                             M
*		                                                             *
* KEYWORDS:								     M
*   MvarZeroOrganizeMVs1DProblem					     M
*****************************************************************************/
MvarMVStruct **MvarZeroOrganizeMVs1DProblem(const MvarMVStruct * const *MVs,
					    MvarConstraintType *Constraints,
					    int *NumOfMVs)
{
    int i;
    MvarMVStruct *MV,
	**LclMVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *)
							         * *NumOfMVs);

    for (i = 0; i < *NumOfMVs; i++) {
	if (MVAR_IS_RATIONAL_MV(MVs[i])) {
	    /* Convert P1 point type to E1, in place. */
	    LclMVs[i] = MV = MvarMVCopy(MVs[i]);
	    MV -> PType = MVAR_PT_E1_TYPE;
#	    ifndef MVAR_MALLOC_STRUCT_ONCE
		IritFree(MV -> Points[0]);
#	    endif /* MVAR_MALLOC_STRUCT_ONCE */
	    MV -> Points[0] = NULL;
	}
	else
	    LclMVs[i] = MvarMVCopy(MVs[i]);

	if (MVAR_IS_BEZIER_MV(LclMVs[i])) {
	    MV = MvarCnvrtBzr2BspMV(LclMVs[i]);
	    MvarMVFree(LclMVs[i]);
	    LclMVs[i] = MV;
	}
    }

    return LclMVs;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test if we can rule out the possibility of a solution for the problem    *
* in its domain, by inspecting constraint number Ind. Callback function      *
* for the MVs case, 1D solution.					     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*   Ind:	The index of the MV to be tested.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if domain can be purged (no solution can exist), FALSE  *
*		otherwise.						     *
*****************************************************************************/
static CagdBType MvarZeroSolverNoZeroTest1DMVs(MvarZeroPrblmStruct *Problem,
					      int Ind)
{
    return MvarZeroMVConstraintFail(Problem -> U.MVs[Ind],
				    Problem -> Constraints[Ind]);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Solving the MVs representation problem on the boundary of the domain. In *
*   the 1D case, it means we are solving 2 * Dim fully constrained problems. *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*									     *
* RETURN VALUE:								     *
*   void								     *
*****************************************************************************/
static void MvarZeroSolveBoundary1DMVs(MvarZeroPrblmStruct *Problem)
{
    int i, j, k, l,
	Dim = MvarZeroSolverGetDmnDim(Problem);
    CagdRType *TempBound,
	SubdivTol = Problem -> SubdivTol,
	NumericTol = Problem -> NumericTol;
    MvarPtStruct *Aux,
	*StartEndPts = NULL,
	*BoundaryIntersPts = NULL;
    MvarMVStruct
	**MVs = Problem -> U.MVs,
	**BoundaryMVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) 
	* (Dim - 1));
    MvarMVZR1DAuxStruct
	*AS = Problem -> AS;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    MvarMVDomain(MVs[0], AS -> MinDmn, AS -> MaxDmn, -1);

    /* Twice the same for TempBound = MinDmn, MaxDmn. */
    TempBound = AS -> MinDmn;
    for (l = 0; l < 2; l++) {	
	for (i = 0; i < Dim; i++) {
	    for (j = 0; j < Dim - 1; j++)
		BoundaryMVs[j] = MvarMVFromMV(MVs[j], TempBound[i], i);

	    BoundaryIntersPts = MvarMVsZeros0D(BoundaryMVs, AS -> Constraints,
					       Dim - 1, SubdivTol,
					       IRIT_FABS(NumericTol));

	    /* The dimension of inters. points on the boundary is extended */
	    /* by 1 to Dim.						   */
	    for (Aux = BoundaryIntersPts; Aux != NULL; Aux = Aux -> Pnext) {
		MvarPtStruct 
		    *NewPt = MvarPtNew(Dim);

#		ifdef DEBUG_MVAR_MVZR1D_TEST_BNDRY_PTS
		{
		    CagdRType R, S, T;

		    R = MvarMVEval(BoundaryMVs[0], Aux -> Pt)[1];
		    S = MvarMVEval(BoundaryMVs[1], Aux -> Pt)[1];
		    T = MvarMVEval(BoundaryMVs[2], Aux -> Pt)[1];
		    fprintf(stderr, "Point %.8f %.8f %.8f [%.8f %.8f %.8f]\n",
			R, S, T, Aux -> Pt[0], Aux -> Pt[1], Aux -> Pt[2]);

		}
#		endif /* DEBUG_MVAR_MVZR1D_TEST_BNDRY_PTS */

		for (k = 0; k < i; k++)
		    NewPt -> Pt[k] = Aux -> Pt[k];	    
		NewPt -> Pt[i] = TempBound[i];
		for (k = i+1; k < Dim; k++)
		    NewPt -> Pt[k] = Aux -> Pt[k-1];	    

		IRIT_LIST_PUSH(NewPt, StartEndPts);
	    } 
	    for (j = 0; j < Dim - 1; j++)
		MvarMVFree(BoundaryMVs[j]);
	    MvarPtFreeList(BoundaryIntersPts);
	} 	
	TempBound = AS -> MaxDmn;
    }

    IritFree(BoundaryMVs);

    if (CagdListLength(StartEndPts) > 1)
	StartEndPts = MvarMVZR1DListOfDifferentPts(StartEndPts, NumericTol);

    PrblmIntrnl -> SolutionsOnBoundary = 
	MvarZeroSolverSolutionNew(NULL, NULL, StartEndPts, 
				  MVAR_ZER_SLVR_SOLUTION_PT_LIST);

    PrblmIntrnl -> NumOfBoundarySolutions = 
	PrblmIntrnl -> SolutionsOnBoundary ?
		CagdListLength(PrblmIntrnl -> SolutionsOnBoundary -> U.Pt) :
	        0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Topological guarantee in the MVs, 1D solution case. The tests involved   *
*   verify that no loops exist and that there are only an entry point and an *
*   exit point. Note that the tests assume that no situations of tangency    *
*   occur (between solution curve and domain boundary).			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if the solution curve is homeomorphic to an interval,   *
*		FALSE if cannot guarantee.				     *
*****************************************************************************/
static CagdBType MvarZeroSolverTopoTest1DMVs(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    /* The no loop property, if TRUE, is inherited. Test only if needed. */
    if (!PrblmIntrnl -> NoLoopTest1D)
	PrblmIntrnl -> NoLoopTest1D = 
	                  MvarMVZR1DNoLoopTest((MvarMVStruct const * const *)
					                    Problem -> U.MVs,
					        Problem -> AS);
    if (PrblmIntrnl -> NoLoopTest1D) {
	int n = PrblmIntrnl -> NumOfBoundarySolutions;

#ifdef DEBUG
	if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
	    printf("The current domain passed the No Loop Test.\n");
	    printf("There are %d solutions on the boundary.\n", n);
	}
#endif /* DEBUG */

	if (2 < n)
	    return FALSE;
	if (2 > n) 
	    /* If the No Loop Test is TRUE and less than two boundary    */
	    /* points, we can guarantee the topology: it is empty. The   */
	    /* actual purging shall occur in the numeric attempt, next.  */
	    return TRUE;
	else /* Exactly one entry point and one exit point. */
	    return TRUE;
    }
    return FALSE; /* Did not pass the No Loop Test. */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   The numeric improvement for the MVs representation problem, 1D case.     *
* Given the start and end points, reconstructs the curve, assuming it is     *
* guaranteed to be homeomorphic to an interval (single component, and no     *
* loops.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution, or NULL upon failure of the     *
*			       process.					     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep1DMVs(
					         MvarZeroPrblmStruct *Problem)
{
    int TraceError = FALSE;
    CagdRType 
	NumericTol = Problem -> NumericTol,
	Step = Problem -> StepTol;
    MvarVecStruct *StartVec;
    MvarPtStruct *SolutionPts,
	*BoundaryPts;
    MvarPolylineStruct *NewPoly;
    MvarMVStruct 
	**MVs = Problem -> U.MVs;
    MvarMVZR1DAuxStruct 
	*AS = Problem -> AS;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;   

    if (PrblmIntrnl -> NoLoopTest1D && 
	PrblmIntrnl -> NumOfBoundarySolutions < 2) {
#ifdef DEBUG
	    if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
		printf("This domain is purged: No Loop and less than 2 solutions on boundary.\n");
		printf("\n**********************************************************\n");
	    }
#endif /* DEBUG */
	PrblmIntrnl -> PurgeByNumericStep = TRUE;
	return NULL;
    }
    /* Make sure we are not here too early: */
    assert(PrblmIntrnl -> ConstructionComplete);

    BoundaryPts = PrblmIntrnl -> SolutionsOnBoundary -> U.Pt;
    StartVec = MvarMVZR1DStartVec(MVs, BoundaryPts,
				  AS, NumericTol);
    SolutionPts = MvarMVZR1DCurveTracing(MVs, BoundaryPts, 
					 BoundaryPts -> Pnext, 
					 StartVec, Step, NumericTol, 
					 PrblmIntrnl -> MVGradients, AS,
					 &TraceError);
    if (TraceError) {
#ifdef DEBUG
	if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
	    printf("Tracing error now.\n"); 
	}
#endif /* DEBUG */
	MvarPtFreeList(SolutionPts);
	MvarVecFree(StartVec);
	return NULL;
    }
#ifdef DEBUG
    if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
	if (SolutionPts != NULL) {
	    printf("Successful tracing now.\n");
	    printf("\n********************************************************\n");
	}
    }
#endif /* DEBUG */
    assert(SolutionPts != NULL);
    MvarVecFree(StartVec);
    /* Each polyline is provided by the current domain. */
    NewPoly = MvarPolylineNew(SolutionPts);
    NewPoly = MvarMVZR1DPolyWithDom(NewPoly, MVs[0]);
    return MvarZeroSolverSolutionNew(NULL, NewPoly, NULL, 
				     MVAR_ZER_SLVR_SOLUTION_POLYLINE);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Subdivision in the 1D, MVs case. Given a problem, returns a vector of    *
* two subdivided problems, allocated dynamically. The subdivided problems    *
* can be NULL.								     *
*									     *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroPrblmStruct **:  The array of two sub-problem.		     *
*****************************************************************************/
static MvarZeroPrblmStruct **MvarZeroSubdivProblem1DMVs(MvarZeroPrblmStruct
							            *Problem)
{
    int i, JLoc,
        NumOfMVs = Problem -> NumOfConstraints;
    CagdRType *MinDmn, *MaxDmn, t, MaxSide,
	NumericTol = Problem -> NumericTol,
	SubdivTol = Problem -> SubdivTol; 
    MvarPtStruct *MiddlePts1, *MiddlePts2, *SplitPts[2], *BoundaryPts1,
	*BoundaryPts2,
	*BoundaryPts = NULL;
    MvarMVStruct **MVs1, **MVs2,
	**MVs = Problem -> U.MVs;
    MvarZeroSolutionStruct *BoundarySol1, *BoundarySol2;
    MvarMVZR1DAuxStruct 
	*AS = Problem -> AS; 
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    MvarZeroPrblmStruct
	**SubProblems = (MvarZeroPrblmStruct **)
		                 IritMalloc(2 * sizeof(MvarZeroPrblmStruct *));
    			
    SubProblems[0] = SubProblems[1] = NULL;

    if (PrblmIntrnl -> SolutionsOnBoundary != NULL)
	BoundaryPts = PrblmIntrnl -> SolutionsOnBoundary -> U.Pt;
    MinDmn = PrblmIntrnl -> MVMinDmn;
    MaxDmn = PrblmIntrnl -> MVMaxDmn;
    JLoc = PrblmIntrnl -> MaxSideDir;
    MaxSide = PrblmIntrnl -> MaxSide;
    MVs1 = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * NumOfMVs);
    MVs2 = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * NumOfMVs);

    /* Subdivision at the middle of the domain, at    		     */
    /* t = (TMin + TMax) * 0.5, is more likely to create problems of */
    /* intersections along (subdivided) domains.  Subdivisions at a  */
    /* small perturbed domain is going to be more robust. 	     */
    t = MVAR_MVZR1D_SUBD_PRTRB * (1 + JLoc * 0.1);
    assert(t < 0.1);
    t = (0.5 + t) * MinDmn[JLoc] + (0.5 - t) * MaxDmn[JLoc];
    for (i = 0; i < NumOfMVs; i++) {
	MVs1[i] = MvarMVSubdivAtParam(MVs[i], t, JLoc);
	MVs2[i] = MVs1[i] -> Pnext;
	MVs1[i] -> Pnext = NULL;
    }

    /* Remember the subdiv info in the parent problem: */
    PrblmIntrnl -> SubdivDir = JLoc;
    PrblmIntrnl -> ParamOfSubdiv = t;

#ifdef DEBUG
    if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
	printf("\nSubdiv is at JLoc = %d, t = %.3f, into:\n", JLoc, t);
	printf("\nSubdomain 1:");
	for (i = 0; i < MVs1[0] -> Dim; i++) {
	    CagdRType t1, t2;

	    MvarMVDomain(MVs1[0], &t1, &t2, i);
	    printf("[%6.3f, %6.3f] ", t1, t2);
	}
	printf("\nSubdomain 2:");
	for (i = 0; i < MVs2[0] -> Dim; i++) {
	    CagdRType t1, t2;

	    MvarMVDomain(MVs2[0], &t1, &t2, i);
	    printf("[%6.3f, %6.3f] ", t1, t2);
	}
	printf("\n**********************************************************\n");
    }
#endif /* DEBUG */

    /* Find middle points and make sure they are not the same as the */
    /* input boundary points.				             */
    MiddlePts1 = MvarMVZR1DMiddlePlaneCutPts(MVs, AS, SubdivTol, 
				             NumericTol, JLoc, t);
    MiddlePts2 = MvarPtCopyList(MiddlePts1);
    MvarMVZR1DSplitBoundaryPts(BoundaryPts, JLoc, t, &SplitPts[0],
			       &SplitPts[1]);
    SplitPts[0] = CagdListAppend(MiddlePts1, SplitPts[0]);
    SplitPts[1] = CagdListAppend(MiddlePts2, SplitPts[1]);
    BoundaryPts1 = MvarMVZR1DListOfDifferentPts(SplitPts[0], NumericTol);
    BoundaryPts2 = MvarMVZR1DListOfDifferentPts(SplitPts[1], NumericTol);

#   ifdef DEBUG_DUMP_SINGULAR_CASES
    if (CagdListLength(BoundaryPts1) % 2 != 0 ||
	CagdListLength(BoundaryPts2) % 2 != 0) {
        MvarPtStruct *Pt;
        int Dim = MvarZeroSolverGetDmnDim(Problem);

	fprintf(stderr, "ERROR: Odd number of boundary points, JLoc = %d at t = %f\n",
		JLoc, t);

	for (Pt = BoundaryPts; Pt != NULL; Pt = Pt -> Pnext) {
	    fprintf(stderr, "\tInput Pts:");
	    for (i = 0; i < Dim; i++) {
	        fprintf(stderr, " %.5lf", Pt -> Pt[i]);
	    }

	    for (i = 0; i < Dim - 1; i++) {
	        CagdRType
		    *R = MvarMVEval(MVs[i], Pt -> Pt);

		if (IRIT_ABS(R[1]) < IRIT_UEPS)
		    fprintf(stderr, " [0]");		      
		else
		    fprintf(stderr, " [%.5g]", R[1]);
	    }

	    fprintf(stderr, "\n");
	}

	MiddlePts1 = MvarMVZR1DMiddlePlaneCutPts(MVs, AS, SubdivTol,
					         NumericTol, JLoc, t);
	for (Pt = MiddlePts1; Pt != NULL; Pt = Pt -> Pnext) {
	    fprintf(stderr, "\tMiddle Pts:");
	    for (i = 0; i < Dim; i++) {
	        fprintf(stderr, " %.5lf", Pt -> Pt[i]);
	    }

	    for (i = 0; i < Dim - 1; i++) {
	        CagdRType
		    *R = MvarMVEval(MVs[i], Pt -> Pt);

		if (IRIT_ABS(R[1]) < IRIT_UEPS)
		    fprintf(stderr, " [0]");		      
		else
		    fprintf(stderr, " [%.5g]", R[1]);
	    }

	    fprintf(stderr, "\n");
	}
	MvarPtFreeList(MiddlePts1);

	for (Pt = BoundaryPts1; Pt != NULL; Pt = Pt -> Pnext) {
	    fprintf(stderr, "\tPts1:");
	    for (i = 0; i < Dim; i++) {
	        fprintf(stderr, " %.8lf", Pt -> Pt[i]);
	    }
	    fprintf(stderr, "\n");
	}
	for (Pt = BoundaryPts2; Pt != NULL; Pt = Pt -> Pnext) {
	    fprintf(stderr, "\tPts2:");
	    for (i = 0; i < Dim; i++) {
	        fprintf(stderr, " %.8lf", Pt -> Pt[i]);
	    }
	    fprintf(stderr, "\n");
	}
    }
#   endif /* DEBUG_DUMP_SINGULAR_CASES */

    BoundarySol1 = MvarZeroSolverSolutionNew(NULL, NULL, BoundaryPts1, 
					    MVAR_ZER_SLVR_SOLUTION_PT_LIST);
    SubProblems[0] = MvarZeroSolverSubProblem(Problem, MVs1, NULL, 
					      BoundarySol1);
    BoundarySol2 = MvarZeroSolverSolutionNew(NULL, NULL, BoundaryPts2, 
					    MVAR_ZER_SLVR_SOLUTION_PT_LIST);
    SubProblems[1] = MvarZeroSolverSubProblem(Problem, MVs2, NULL, 
					      BoundarySol2);
    return SubProblems;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returning a solution in the case of reaching subdivision tolerance.      *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution.				     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol1DMVs(MvarZeroPrblmStruct
							             *Problem)
{
    MvarPolylineStruct *Segment;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    if (GlblMvarInterMergeSingularPts > 0) {
        /* If top level is smaller than subdivision tolerance: */
        if (PrblmIntrnl -> SubdivDepth == 0) {
            MVAR_ZERO_SLVR_APPLY(SolveBoundary)(Problem);
	}

        /* If we have exactly 2 boundary solutions - Attempt to trace. If    */
	/* failed, return the edge segment from start to end.                */
	if (PrblmIntrnl -> NumOfBoundarySolutions == 2)	{
	    MvarZeroSolutionStruct *Sol;

	    if (PrblmIntrnl -> ConstructionComplete == FALSE) {
		/* When under subdivision tolerance - it may happen that we  */
		/* do not have the gradients yet, etc.                       */
		if (!MVAR_ZERO_SLVR_APPLY(FirstSmoothUpdates)(Problem)) {
		    return NULL;
		}
	    }
	    Sol = MvarZeroSolverNumericStep1DMVs(Problem);

	    if (Sol != NULL)
		return Sol;
	    else {
		/* Sol is NULL, but we assume the two boundary points are   */
		/* indeed an entry and an exit. Return the line segment.    */
		Segment = MvarPolylineNew(MvarPtCopyList(
		                  PrblmIntrnl -> SolutionsOnBoundary -> U.Pt));
		return MvarZeroSolverSolutionNew(
		                              NULL, Segment, NULL, 
		                              MVAR_ZER_SLVR_SOLUTION_POLYLINE);
	    }
	}
	else
	    return NULL;
    }
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Unites two neighboring 1D solutions of a zero finding problem, MVs	     *
* representation, back to a new solution of the subdivided domain:           *
* the polylines are linked, making sure that in the common boundary we have  *
* identical points, so the linkage is consistent.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol1:	The first solution.					     *
*   Sol2:	The second solution.					     *
*   Dir:	The direction along which the subdivision occurred.	     *
*   Param:	The parameter value at which the subdivision occurred, along *
*		direction Dir.						     *
*   Problem:	The zero finding problem solved by the united solution.      *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The new solution, union of Sol1 and Sol2.     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions1DMVs(
						MvarZeroSolutionStruct *Sol1,
						MvarZeroSolutionStruct *Sol2,
					        MvarMVDirType Dir,
					        CagdRType Param,
						MvarZeroPrblmStruct *Problem)
{
    CagdRType
	SubdivTol = Problem -> SubdivTol,
	NumericTol = Problem -> NumericTol;
    MvarPolylineStruct *NewPoly, 
	*List1 = NULL,
	*List2 = NULL;
    MvarZeroSolutionStruct *Solution;

    if (Sol1 == NULL) {
	return Sol2;
    }
    if (Sol2 == NULL) {
	return Sol1;
    }

#ifdef DEBUG
    assert(Sol1 -> ActiveRepresentation == MVAR_ZER_SLVR_SOLUTION_POLYLINE &&
	   Sol2 -> ActiveRepresentation == MVAR_ZER_SLVR_SOLUTION_POLYLINE);
#endif /* DEBUG */

    /* Both are not NULL. Unite them in place, that is: use the structure of */
    /* Sol1, and only the U.Pt from Sol2. This explains the deallocation     */
    /* logic that is typical after returning from this routine.		     */
    List1 = Sol1 -> U.Pl;
    List2 = Sol2 -> U.Pl;
    NewPoly = MvarMVZR1DLinkNeighbours(List1, List2, SubdivTol, 
				       NumericTol, Dir, Param);
    Solution = MvarZeroSolverSolutionNew(NULL, NewPoly, NULL, 
				     MVAR_ZER_SLVR_SOLUTION_POLYLINE);
    MvarZeroSolverSolutionFree(Sol1, FALSE);
    MvarZeroSolverSolutionFree(Sol2, FALSE);
    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* In the 1D solution case, as for now we do nothing. Can be used for some    *
* organization like in the 0D case if needed in the future.		     *
*									     *
* PARAMETERS:                                                                *
*   Sol:	The solution to be organized.				     *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct:  The new, organized solution.		     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol1DMVs(
						 MvarZeroSolutionStruct *Sol,
						 MvarZeroPrblmStruct *Problem)
{
    if (Sol != NULL) {
	int NumOfMVs = Problem -> NumOfConstraints;
	CagdRType 
	    NumericTol = IRIT_ABS(Problem -> NumericTol);
	MvarPolylineStruct *MVPl, 
	    *ZeroSet = Sol -> U.Pl;
	MvarConstraintType
	    *Constraints = Problem -> Constraints;
	MvarMVStruct 
	    **LclMVs = Problem -> U.MVs;

	/* Inequality constraints verification: */
	for (MVPl = ZeroSet; MVPl != NULL; MVPl = MVPl -> Pnext) {
	    MVPl -> Pl = MvarZeroFilterSolutionSet(
					MVPl -> Pl,
					(const MvarMVStruct * const *) LclMVs,
					Constraints, NumOfMVs,
					IRIT_ABS(NumericTol), TRUE, FALSE, 
					TRUE);
	}

	/* Remove all empty polylines. */
	while (ZeroSet != NULL && ZeroSet -> Pl == NULL) {
	    MVPl = ZeroSet -> Pnext;
	    MvarPolylineFree(ZeroSet);
	    ZeroSet = MVPl;	    
	}
	Sol -> U.Pl = ZeroSet;
	if (ZeroSet != NULL) {
	    MvarPolylineStruct *MVPlTmp; 

	    for (MVPlTmp = ZeroSet; MVPlTmp -> Pnext != NULL; ) {
	        if (MVPlTmp -> Pnext -> Pl == NULL) {
		    MVPl = MVPlTmp -> Pnext -> Pnext;
		    MvarPolylineFree(MVPlTmp -> Pnext);
		    MVPlTmp -> Pnext = MVPl;
		}
		else {
		    MVPlTmp = MVPlTmp -> Pnext;
		}
	    }
	}
    }
    if (Sol != NULL && Sol -> U.Pl == NULL) {
	MvarZeroSolverSolutionFree(Sol, TRUE);
	return NULL;
    }
    return Sol;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the state of the singular points merger:                            M
*   If 0, singular locations are ignored (skipped).                          M
*   If 1, singular locations are merged using the subdivision tolerance      M
*         which improves the changes of a complete long merged curves.       M
*   If 2, singular locations are merged using the numeric tolerances         M
*         (like every other case) which means most likely they will be left  M
*         as isolated points.                                                M
*                                                                            *
* PARAMETERS:                                                                M
*   MergeSingularPts:   Set the desired state of singular points mergers.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    Old state.                                                       M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSrfSrfInter, MvarMVsZeros1D                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZeros1DMergeSingularPts                                           M
*****************************************************************************/
int MvarMVsZeros1DMergeSingularPts(int MergeSingularPts)
{
    int OldVal = MergeSingularPts;

    GlblMvarInterMergeSingularPts = MergeSingularPts;

    return OldVal;
}
