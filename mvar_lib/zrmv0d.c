/******************************************************************************
* ZrMV0D.c - tools to compute zero sets of multivariates when the problem has *
*	     MVs representation and the expected solution  set                *
*	     is a zero-manifold: a finite set of points.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, Feb. 13.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"

IRIT_GLOBAL_DATA MvarZeroSubdivTolActionType
    _MVGlblUponSubdivTol = MVAR_ZERO_NUMERIC_STEP_VERIFIED;

static CagdBType MvarZeroSolverNoZeroTest0DMVs(MvarZeroPrblmStruct *Problem,
					       int Ind);
static void MvarZeroSolveBoundary0DMVs(MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverTopoTest0DMVs(MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep0DMVs(
					         MvarZeroPrblmStruct *Problem);
static MvarZeroPrblmStruct **MvarZeroSubdivProblem0DMVs(MvarZeroPrblmStruct
							             *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol0DMVs(MvarZeroPrblmStruct
							             *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions0DMVs(
						 MvarZeroSolutionStruct *Sol1,
						 MvarZeroSolutionStruct *Sol2, 
					         MvarMVDirType Dir,
					         CagdRType Param,
						 MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol0DMVs(
						   MvarZeroSolutionStruct *Sol,
						   MvarZeroPrblmStruct 
					                             *Problem);

IRIT_STATIC_DATA const MvarZeroSolverCallBackFcnStruct 
    ZrMv0DCB = {
        MvarZeroSolveBoundary0DMVs,
        MvarZeroSolverNoZeroTest0DMVs,
        MvarZeroSolverTopoTest0DMVs,
	MvarZeroSubdivProblem0DMVs,
	MvarZeroSolverNumericStep0DMVs,
	MvarZeroSolverSinglrSol0DMVs,
	MvarZeroSolverUniteSolutions0DMVs,
	MvarZeroSolverOrganizeSol0DMVs,
	MvarZeroUpdateProblemDmnMVs,
	MvarZeroFirstSmoothUpdatesMVs,
	NULL
    };

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the callback functions for the MVs representation, 0D solution case.M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void								     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroOrganizeMVs0DProblem	    	                             M
*		                                                             *
* KEYWORDS:								     M
*   MvarZeroSolverSetCallbackFcns0DMVs					     M
*****************************************************************************/
void MvarZeroSolverSetCallbackFcns0DMVs(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    IRIT_GEN_COPY(&PrblmIntrnl -> CallbackFunctions, &ZrMv0DCB,
		  sizeof(MvarZeroSolverCallBackFcnStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Preliminary organization and scaling actions of the objects composing a  M
* zero finding problem.  Relevant to the MVs representation only. This       M
* routine is invoked by the problem structure construction routine, and      M
* should not be invoked when extracting a sub-problem from an existing one.  M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:		Array of multivariates.				     M
*   Constraints:	Equality or inequality constraints.		     M
*   NumOfMVs:		Number of constraints (may be updated).		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct **:  The updated array of MVs to be used.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverSetCallbackFcns0DMVs    	                             M
*		                                                             *
* KEYWORDS:								     M
*   MvarZeroOrganizeMVs0DProblem					     M
*****************************************************************************/
MvarMVStruct **MvarZeroOrganizeMVs0DProblem(const MvarMVStruct * const *MVs,
					    MvarConstraintType *Constraints,
					    int *NumOfMVs)
{
    CagdBType
	HasBezConst = FALSE,
	HasBspConst = FALSE;
    int i, l;
    MvarMVStruct **LclMVs;

#   ifdef MVAR_DEBUG_DUMP_INPUT
    for (l = 0; l < *NumOfMVs; l++) {
        printf("MVAR SOLVER INPUT *******************\nMV %d: %s type, ", l,
		MVAR_IS_BEZIER_MV(MVs[l]) ? "Bzr" : (MVAR_IS_BSPLINE_MV(MVs[l]) ? "Bsp" : "Other"));
	printf("orders: ");
	for (i = 0; i < MVs[l] -> Dim; i++)
	    printf("%d, ", MVs[l] -> Orders[i]);
	if (MVAR_IS_BSPLINE_MV(MVs[l])) {
	    printf("lengths: ");
	    for (i = 0; i < MVs[l] -> Dim; i++)
	        printf("%d, ", MVs[l] -> Lengths[i]);
	}
	printf("\n");
    }

    for (l = 0; l < *NumOfMVs; l++) {
        fprintf(stderr, "MV %d *******************************\n", l);
	MvarDbg(MVs[l]);
    }
#   endif /* MVAR_DEBUG_DUMP_INPUT */

    /* Make sure the parametric domain is identical at all multivariates. */
    for (l = 0; l < MVs[0] -> Dim; l++) {
	CagdRType Min, Max;

	MvarMVDomain(MVs[0], &Min, &Max, l);

	for (i = 1; i < *NumOfMVs; i++) {
	    CagdRType Min2, Max2;

	    MvarMVDomain(MVs[i], &Min2, &Max2, l);
	    if (!IRIT_APX_EQ(Min, Min2) || !IRIT_APX_EQ(Max, Max2)) {
		MVAR_FATAL_ERROR(MVAR_ERR_INCONS_DOMAIN);
		return NULL;
	    }
	}
    }

    /* Make sure we have domains (KVs) on all multivariates and all         */
    /* multivariates are polynomials,  either all Bezier or all B-spline.   */
    LclMVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct) * (*NumOfMVs));
    for (i = 0; i < *NumOfMVs; i++) {
	LclMVs[i] = MvarMVCopy(MVs[i]);
	if (MVAR_IS_BEZIER_MV(LclMVs[i])) {
	    MvarMVAuxDomainSlotReset(LclMVs[i]);
	    HasBezConst = TRUE;
	} 
	else if (MVAR_IS_BSPLINE_MV(MVs[i])) {
	    HasBspConst = TRUE;
	}
 
	if (MVAR_IS_RATIONAL_MV(LclMVs[i])) {
	    /* Convert P1 point type to E1, in local copy. */
#ifndef MVAR_MALLOC_STRUCT_ONCE
	    IritFree(LclMVs[i] -> Points[0]);
#endif /* MVAR_MALLOC_STRUCT_ONCE */
	    LclMVs[i] -> Points[0] = NULL;
	    LclMVs[i] -> PType = MVAR_PT_E1_TYPE;
	}
    }

    if (HasBspConst && HasBezConst) {
        MVAR_FATAL_ERROR(MVAR_ERR_CANNT_MIX_BSP_BEZ);
        return NULL;
    }

    /* Preconditionally scale all constraints to a 'reasonable' height. Note */
    /* such scale should not affect the zeros of the set of constraints.     */
    for (i = 0; i < *NumOfMVs; ) {
	int Length = MVAR_CTL_MESH_LENGTH(LclMVs[i]);
	CagdRType Scale, *R,
	    *Pts = LclMVs[i] -> Points[1],
	    Max = 0.0;

	for (R = &Pts[Length]; R-- != Pts; )
	    Max = IRIT_MAX(Max, IRIT_FABS(*R));

	/* For the unlikely case of a constraint that is identically zero. */
	if (Max == 0.0) {
	    /* Purge away this identically zero constraint. */
	    MvarMVFree(LclMVs[i]);
	    LclMVs[i] = LclMVs[*NumOfMVs - 1];
	    Constraints[i] = Constraints[*NumOfMVs - 1];
	    (*NumOfMVs)--;
	    continue;
	}
	else
	    i++;

	Scale = MVAR_ZERO_PRECOND_SCALE / Max;

	for (R = &Pts[Length]; R-- != Pts; )
	    *R *= Scale;
    }

    /* For the unlikely case of all constraints identically zero. */
    if (*NumOfMVs == 0) {
	IritFree(LclMVs);
	return NULL;
    }

    return LclMVs;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Testing if we can rule out the possibility of a solution for the problem   *
* in its domain, by inspecting constraint number Ind. Callback function      *
* for the MVs case.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*   Ind:	The index of the MV to be tested.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if domain can be purged (no solution can exist), FALSE  *
*		otherwise.						     *
*****************************************************************************/
static CagdBType MvarZeroSolverNoZeroTest0DMVs(MvarZeroPrblmStruct *Problem,
					       int Ind)
{
    /* The pairs of hyperplanes test may rule out the domain. If the test */
    /* is to be used (appropriate global is TRUE), apply it only for 0D   */
    /* problems and only at the first call (it looks at all constraints,  */
    /* not just number Ind). */
    if (Ind == 0 && MVAR_ZERO_SLVR_SOLUTION_SPACE_DIM(Problem) == 0 &&
	_MVGlblZeroApplyParallelHyperPlaneTest && 
	!MVarMVHyperPlanesTestForSol((MvarMVStruct const * const *) 
				                           Problem -> U.MVs, 
	                             Problem -> NumOfZeroConstraints))
	return TRUE;

    return MvarZeroMVConstraintFail(Problem -> U.MVs[Ind], 
				    Problem -> Constraints[Ind]);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Solving the problem on the boundary of the domain. De-facto we do        *
* nothing - the problem is already fully constrained. The solutions that     *
* that belong to the boundary should rise as part of the process.	     *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*		                                                             *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroSolveBoundary0DMVs(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    assert(PrblmIntrnl -> SolutionsOnBoundary == NULL);
    PrblmIntrnl -> SolutionsOnBoundary = NULL;
    PrblmIntrnl -> NumOfBoundarySolutions = 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Topological guarantee in the MVs, 0D solution case: here we test         *
* that at most a single solution can exist in the domain. This is done via   *
* the non-overlapping criteria of the tangent cones.                         *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if at most a single solution, FALSE otherwise.	     *
*****************************************************************************/
static CagdBType MvarZeroSolverTopoTest0DMVs(MvarZeroPrblmStruct *Problem)
{
    CagdBType
	ApplyNormalConeTest = _MVGlblZeroApplyNormalConeTest;    
    int NumOfZeroMVs = Problem -> NumOfZeroConstraints;
    MvarMVStruct 
	**MVs = Problem -> U.MVs;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    /* If we already guaranteed at most a single solution at a larger domain,*/
    /* then we are guaranteed here as well, no need to check:		     */
    if (PrblmIntrnl -> SingleSolutionTest0D)
	return TRUE;

    /* Check the normal cone overlapping criteria. */
    if (ApplyNormalConeTest && !MvarMVConesOverlap(MVs, NumOfZeroMVs)) {
#	ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintZeroPts, FALSE) {
	        int i,
		    Dim = MvarZeroSolverGetDmnDim(Problem);
		CagdRType Min, Max;

		IRIT_INFO_MSG("Zero by cones at");
		for (i = 0; i < Dim; i++) {
		    MvarGetSubdivParamDomains(MVs[0], &Min, &Max, i);
		    IRIT_INFO_MSG_PRINTF(" %f [%f %f]",
					(Min + Max) * 0.5, Min, Max);
		}
		IRIT_INFO_MSG("\n");
	    }
	}
#       endif /* DEBUG */
	/* Remember, for all future sub-domains, that this is guaranteed: */
	PrblmIntrnl -> SingleSolutionTest0D = TRUE;
	return TRUE;
    }
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   The numeric improvement for a single point in a domain of the zero       *
* finding problem, for the MVs, 0D solution case.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution: a single point or NULL.	     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep0DMVs(
					         MvarZeroPrblmStruct *Problem)
{
    int i,
	NumOfZeroMVs = Problem -> NumOfZeroConstraints;
    CagdRType
	NumericTol = Problem -> NumericTol;
    MvarPtStruct *ZeroPt;
    MvarMVStruct **LclMVs;
    MvarMVGradientStruct **Grads;
    MvarZeroSolutionStruct
	*Solution = NULL;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    /* Sanity check: we can arrive here only with guaranteed topology, or   */
    /* if subdivision tolerance has been reached, and this is the chosen    */
    /* action by the appropriate global flag:				    */
    if (PrblmIntrnl -> SingleSolutionTest0D == FALSE) {
	if (_MVGlblUponSubdivTol != MVAR_ZERO_NUMERIC_STEP_VERIFIED && 
	    PrblmIntrnl -> UnderSubdivTol) { 
	    MVAR_FATAL_ERROR(MVAR_ERR_ZER_NUMERIC_TOO_EARLY);
	    return NULL;
	}
    }
    /* For unlikely singular cases where the first smooth problem is not    */
    /* yet detected (should happen here only if we are under subdiv tol').  */
    if (!PrblmIntrnl -> ConstructionComplete) {
	MVAR_ZERO_SLVR_APPLY(FirstSmoothUpdates)(Problem);
    }
    LclMVs = PrblmIntrnl -> U.FirstSmoothMVs;

    Grads = Problem -> _Internal -> MVGradients;
    assert(Grads != NULL);  /* Otherwise we prep below and is a mem leak... */
    for (i = 0; i < NumOfZeroMVs; i++) {
	AttrSetRefPtrAttrib(&LclMVs[i] -> Attr, "Gradient",
			    Grads != NULL ? Grads[i]
			                  : MvarMVPrepGradient(LclMVs[i], TRUE));
    }
    ZeroPt = MvarZeroGenPtMidDmn(Problem, TRUE);

    /* Restrict domain to MVs domain. */
    {
	const MvarMVStruct
	    *MV = Problem -> U.MVs[0];
	CagdRType 
	    *MinDmn = (CagdRType *) IritMalloc(sizeof(CagdRType) * MV -> Dim),
	    *MaxDmn = (CagdRType *) IritMalloc(sizeof(CagdRType) * MV -> Dim);

	MvarMVDomain(MV, MinDmn, MaxDmn, -1);
	ZeroPt = MvarZero0DNumeric(ZeroPt, NULL,
				   (MvarMVStruct const * const *) LclMVs,
				   NumOfZeroMVs, NumericTol, MinDmn, MaxDmn);
	IritFree(MinDmn);
	IritFree(MaxDmn);
    }

    if (ZeroPt) {
        /* If domain extension was applied - make sure the point is in the  */
        /* domain of interest. Remember that if the domain was extended,    */
        /* the MVs that were extended are the first smooth ones, and they   */
        /* were the ones used in the numeric step. This might have returned */
        /* points outside the original problem's domain.		    */
	if (_MVGlblZeroDmnExtension > 0.0) {
	    CagdRType CurDimParam, CurMinDiff, CurMaxDiff,
		LclNumericTol = IRIT_ABS(NumericTol),
		*OrigMinDmn = PrblmIntrnl -> OrigMVMinDmn,
		*OrigMaxDmn = PrblmIntrnl -> OrigMVMaxDmn;

	    for (i = 0; i < ZeroPt -> Dim; i++) {
		/* Local convention: positive diff' means outside domain. */
		CurDimParam = ZeroPt -> Pt[i];
		CurMinDiff = OrigMinDmn[i] - CurDimParam;
		CurMaxDiff = CurDimParam - OrigMaxDmn[i];
		if (CurMinDiff <= 0.0 && CurMaxDiff <= 0.0)
		    /* So far - inside. Carry on. */
		    continue;
		else if (CurMinDiff > 0.0 && CurMinDiff < LclNumericTol) 
		    /* Outside but close enough. Trim to the min' boundary. */
		    ZeroPt -> Pt[i] = OrigMinDmn[i];
		else if (CurMaxDiff > 0.0 && CurMaxDiff < LclNumericTol)
		    /* Outside but close enough. Trim to the max' boundary. */
		    ZeroPt -> Pt[i] = OrigMaxDmn[i];
		else {		   /* Too far outside, at one of the sides. */
		    MvarPtFree(ZeroPt);
		    return NULL;
		}
	    }
	}

	Solution = MvarZeroSolverSolutionNew(NULL, NULL, ZeroPt, 
					     MVAR_ZER_SLVR_SOLUTION_PT_LIST);
    }

    return Solution;
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
static MvarZeroPrblmStruct **MvarZeroSubdivProblem0DMVs(MvarZeroPrblmStruct
							             *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    int i, j, l, HasInternalKnot,
	Dim = MvarZeroSolverGetDmnDim(Problem),
	NumOfMVs = Problem -> NumOfConstraints,
	NumOfZeroMVs = PrblmIntrnl -> NumOfZeroSubdivConstraints,
	Depth = PrblmIntrnl -> SubdivDepth;
    CagdRType Min, Max,
	SubdivTol = Problem -> SubdivTol;
    MvarMVStruct 
	**MVs = Problem -> U.MVs;
    MvarConstraintType
	*Constraints = Problem -> Constraints;
    MvarZeroPrblmStruct
	**SubProblems = (MvarZeroPrblmStruct **)
		                 IritMalloc(2 * sizeof(MvarZeroPrblmStruct *));
    			
    SubProblems[0] = SubProblems[1] = NULL;

#   ifdef MVAR_DEBUG_DEPTH
    if (Depth < MVAR_DEBUG_MAX_DEPTH)
        MvarSubdivLevel[Depth]++;
#   endif /* MVAR_DEBUG_DEPTH */

#   ifdef DEBUG_ONE_SOLUTION
    if (NumOfZeroMVs == 2 || NumOfZeroMVs == 3) {
        static int
	    Count = 1;
        CagdRType UMin, UMax, VMin, VMax,
	    WMin = 0.0,
	    WMax = 0.0,
	    u = 0.02223875187684,
	    v = 0.87209881542057,
	    w = 0.43073326427685;

        MvarGetSubdivParamDomains(MVs[0], &UMin, &UMax, 0);
        MvarGetSubdivParamDomains(MVs[0], &VMin, &VMax, 1);
	if (NumOfZeroMVs == 3)
	    MvarGetSubdivParamDomains(MVs[0], &WMin, &WMax, 2);

	if (UMin <= u && u <= UMax &&
	    VMin <= v && v <= VMax &&
	    (NumOfZeroMVs < 3 || (WMin <= w && w <= WMax)))
	    IRIT_INFO_MSG_PRINTF("In domain (%d) [%f %f %f] \t[%f %f %f]\n",
				  Count++, UMin, VMin, WMin, UMax, VMax, WMax);
    }
#   endif /* DEBUG_ONE_SOLUTION */

#   ifdef DEBUG_DUMP_DOMAINS
    {
        CagdRType UMin, UMax, VMin, VMax;

        MvarGetSubdivParamDomains(MVs[0], &UMin, &UMax, 0);
        MvarGetSubdivParamDomains(MVs[0], &VMin, &VMax, 1);

	fprintf(GlblDumpDomainsFile,
		"[OBJECT [RGB \"100,255,100\"] NONE    [POLYLINE 5\n\t[%f %f 0]\n\t[%f %f 0]\n\t[%f %f 0]\n\t[%f %f 0]\n\t[%f %f 0]\n    ]\n]\n",
		UMin, VMin, UMin, VMax, UMax, VMax, UMax, VMin, UMin, VMin);
    }
#   endif /* DEBUG_DUMP_DOMAINS */

    /*   If we got here then these patches may satisfy all constraints.      */
    /*   Subdivide them along their maximal parametric length dimension.     */
    /*   There is one exception to this rule and that is to subdivide first  */
    /* in internal knots so if the maximal parametric length dimension holds */
    /* no knots while another dir do has knots, the other dir will be used.  */
    l = -1;
    HasInternalKnot = FALSE;
    Max = Min = 0.0;
    for (i = 0; i < NumOfMVs && l < 0; i++) {
	for (j = 0; j < Dim && l < 0; j++) {
	    if (MVs[i] -> Lengths[j] != MVs[i] -> Orders[j]) {
		CagdRType TMin, TMax;

		/* Found a direction with internal knots. Save in l only if */
		/* this one is a larger domain than what was found so far.  */
		MvarGetSubdivParamDomains(MVs[0], &TMin, &TMax, j);
		if (TMax - TMin > Max - Min) {
		    Max = TMax;
		    Min = TMin;
		    l = j;
		    HasInternalKnot = TRUE;
		}
	    }
	}
    }

    for (i = 0; i < Dim; i++) {
	CagdRType MinCurr, MaxCurr;

	MvarGetSubdivParamDomains(MVs[0], &MinCurr, &MaxCurr, i);
        if (MaxCurr - MinCurr > Max - Min) {
	    int j;

	    /* If we got internal knots with a domain larger than SubdivTol */
	    /* make sure this direction is a direction with knots.          */
	    if (HasInternalKnot && Max - Min > SubdivTol) {
	        for (j = 0; j < NumOfMVs; j++) {
		    if (MVs[j] -> Lengths[i] != MVs[j] -> Orders[i])
		        break;
		}
	    }
	    else
	        j = 0;

	    if (j < NumOfMVs) {
	        l = i;
		Min = MinCurr;
		Max = MaxCurr;
	    }
        }
    }
    assert(l >= 0);

    if (Max - Min > SubdivTol) {
        CagdBType
	    WasReduction = FALSE;
	CagdRType t,
	    TMin = MVAR_IS_BEZIER_MV(MVs[0]) ? 0.0 : Min,
	    TMax = MVAR_IS_BEZIER_MV(MVs[0]) ? 1.0 : Max;
	MvarMVStruct **MVs1, **MVs2;

	/* MvarMVsReduceMvsDomains returns TMin/TMax in [0, 1] for Bezier   */
	/* and returns correct TMin/TMax for B-spline MVs.		    */
	/*   Also force subdivisions from time to time as domain reduction  */
	/* can fail to converge if applied by itself.			    */
	if (_MVGlblZeroApplyDomainReduction && (Depth & 0x03) != 0) {
	    CagdRType
	        OrigTMin = TMin,
	        OrigTMax = TMax;

	    if (!MvarMVsReduceMvsDomains(MVs, NumOfZeroMVs,
					 l, SubdivTol, &TMin, &TMax,
					 PrblmIntrnl -> ZeroMVsSameSpace))
		return SubProblems; /* The domain reduction ended up empty. */

	    WasReduction = !IRIT_APX_EQ(OrigTMin, TMin) ||
			   !IRIT_APX_EQ(OrigTMax, TMax);
	}

	if (WasReduction) {
	    MVs1 = (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));

	    for (i = 0; i < NumOfMVs; i++) {
	        MVs1[i] = MvarMVRegionFromMV(MVs[i], TMin, TMax, l);
		if (MvarZeroMVConstraintFail(MVs1[i], Constraints[i])) {
		    for ( ; i >= 0; i--)
		        MvarMVFree(MVs1[i]);

		    IritFree(MVs1);
		    return SubProblems; /* Both are NULL, will be purged. */
		}
	    }

	    SubProblems[0] = MvarZeroSolverSubProblem(Problem, MVs1, 
						      NULL, NULL);
	    return SubProblems;
	}
	else {
	    CagdBType
		CanPurgeDomain1 = FALSE,
		CanPurgeDomain2 = FALSE;

	    if (PrblmIntrnl -> ParamPerturb < (TMax - TMin) / 10.0)
	        t = (TMin + TMax) * 0.5 + PrblmIntrnl -> ParamPerturb;
	    else
	        t = (TMin + TMax) * 0.5;

	    /* Lets see if we have a B-spline multivariate with interior    */
	    /* knot in this direction.  If so pick up interior knot instead.*/
	    if (MVAR_IS_BSPLINE_MV(MVs[0])) {
	        for (i = 0; i < NumOfMVs; i++) {
		    if (MVs[i] -> Lengths[l] != MVs[i] -> Orders[l]) {
		        CagdRType
			    r = MVs[i] -> KnotVectors[l][
						   (MVs[i] -> Lengths[l] +
						    MVs[i] -> Orders[l]) >> 1];

			if (r > TMin + IRIT_EPS && r < TMax - IRIT_EPS) {
			    t = r;
			    break;
			}
		    }
		}
	    }

	    /* Ensure we have t within domain, so subdivision can be used. */
	    assert(t >= TMin && t <= TMax);

	    /* Remember the subdiv info in the father problem: */
	    PrblmIntrnl -> SubdivDir = l;
	    PrblmIntrnl -> ParamOfSubdiv = t;

	    MVs1 = (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));
	    MVs2 = (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));

	    for (i = 0; i < NumOfMVs; i++) {
	        if (CanPurgeDomain1) {
		    MVs1[i] = NULL;
		    MVs2[i] = MvarMVSubdivAtParamOneSide(MVs[i], t, l, FALSE);
		    CanPurgeDomain2 |= MvarZeroMVConstraintFail(MVs2[i],
							      Constraints[i]);
		}
		else if (CanPurgeDomain2) {
		    MVs1[i] = MvarMVSubdivAtParamOneSide(MVs[i], t, l, TRUE);
		    MVs2[i] = NULL;
		    CanPurgeDomain1 |= MvarZeroMVConstraintFail(MVs1[i],
							      Constraints[i]);
		}
		else {
		    MVs1[i] = MvarMVSubdivAtParam(MVs[i], t, l);
		    MVs2[i] = MVs1[i] -> Pnext;
		    MVs1[i] -> Pnext = NULL;
		    CanPurgeDomain1 |= MvarZeroMVConstraintFail(MVs1[i],
							      Constraints[i]);
		    CanPurgeDomain2 |= MvarZeroMVConstraintFail(MVs2[i],
							      Constraints[i]);
		}

		if (CanPurgeDomain1 && CanPurgeDomain2) {
		    /* Both domains are purged away - abort here. */
		    for ( ; i >= 0; i--) {
		        MvarMVFree(MVs1[i]);
			MvarMVFree(MVs2[i]);
		    }
		    IritFree(MVs1);
		    IritFree(MVs2);
		    return SubProblems; /* Both are NULL, will be purged. */
		}
	    }

	    if (CanPurgeDomain1) {
	        /* Free allocated MVs as we decided to purge this domain. */
		for (i = 0; i < NumOfMVs; i++) {
		    if (MVs1[i] == NULL)
		        break;
		    MvarMVFree(MVs1[i]);
		}
		IritFree(MVs1);
	    }
	    else {
		SubProblems[0] = MvarZeroSolverSubProblem(Problem, MVs1, 
							  NULL, NULL);
	    }

	    if (CanPurgeDomain2) {
	        /* Free allocated MVs as we decided to purge this domain. */
		for (i = 0; i < NumOfMVs; i++) {
		    if (MVs1[i] == NULL)
		        break;
		    MvarMVFree(MVs2[i]);
		}
		IritFree(MVs2);
	    }
	    else {
		SubProblems[1] = MvarZeroSolverSubProblem(Problem, MVs2, 
							  NULL, NULL);
	    }
	    return SubProblems;
	}
    }

#   ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintZeroPts, FALSE) {
	    CagdRType Min, Max;

	    IRIT_INFO_MSG("Zero by subdivision level at");
	    for (i = 0; i < Dim; i++) {
	        MvarGetSubdivParamDomains(MVs[0], &Min, &Max, i);
		IRIT_INFO_MSG_PRINTF(" %f [%f %f]",
				     (Min + Max) * 0.5, Min, Max);
	    }
	    IRIT_INFO_MSG("\n");
	}
    }
#   endif /* DEBUG */

    /*   If we got here then in all parametric directions, the domain is     */
    /* smaller than the SubdivTol.  Return the respective singular solution. */
    /* should not happen since the current function is to be invoked only    */
    /* above subdivision tolerance!					     */
    assert(0);
    return SubProblems;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returning a solution in the case of reaching subdivision tolerance.      *
*   In the MVs, 0D case the candidate is the center of the domain.	     *
*   Returning it or not is determined by the global flag.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution, in the center of the domain.    *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol0DMVs(MvarZeroPrblmStruct
							             *Problem)
{
    CagdRType Err;
    MvarPtStruct *Pt;
    MvarZeroSolutionStruct *Solution;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    switch (_MVGlblUponSubdivTol) {
	case MVAR_ZERO_ALWAYS_PURGE:
	    Solution = NULL;
	    break;
	case MVAR_ZERO_NEVER_PURGE:
	    Pt = MvarZeroGenPtMidDmn(Problem, TRUE);
	    Solution = MvarZeroSolverSolutionNew(NULL, 
				     NULL, Pt, MVAR_ZER_SLVR_SOLUTION_PT_LIST);
	    break;
	case MVAR_ZERO_RETURN_VERIFIED:
	    Pt = MvarZeroGenPtMidDmn(Problem, TRUE);
	    Err = MvarMVEvalErrorL1((MvarMVStruct const * const *)
				                            Problem -> U.MVs,
				    Pt -> Pt, 
	        		    PrblmIntrnl -> NumOfZeroSubdivConstraints);
	    if (Err < Problem -> NumericTol) 
		Solution = MvarZeroSolverSolutionNew(NULL, 
					       NULL, Pt, 
					       MVAR_ZER_SLVR_SOLUTION_PT_LIST);
	    else
		Solution = NULL;
	    break;
	case MVAR_ZERO_NUMERIC_STEP:
	    Solution = MvarZeroSolverNumericStep0DMVs(Problem);
	    break;
	case MVAR_ZERO_NUMERIC_STEP_VERIFIED:
	    if ((Solution = MvarZeroSolverNumericStep0DMVs(Problem)) != NULL) {
		Err = MvarMVEvalErrorL1((MvarMVStruct const * const *)
				            Problem -> _Internal -> U.FirstSmoothMVs,
					Solution -> U.Pt -> Pt, 
	        			PrblmIntrnl -> NumOfZeroSubdivConstraints);
		if (Err > IRIT_ABS(Problem -> NumericTol)) {
		    MvarZeroSolverSolutionFree(Solution, TRUE);
		    Solution = NULL;
		}
	    }
	    break;
	default:
	    Solution = NULL;
    }
    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Unites two 0D solutions of a zero finding problem, MVs representation,   *
* to a new solution, in place.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol1:	The first solution.					     *
*   Sol2:	The second solution.					     *
*   Dir:	Unused (in the 0D case).				     *
*   Param:	Unused (in the 0D case).				     *
*   Problem:	Unused (in the 0D case).				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct:  The new solution, union of Sol1 and Sol2.	     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions0DMVs(
						 MvarZeroSolutionStruct *Sol1,
						 MvarZeroSolutionStruct *Sol2,
					         MvarMVDirType Dir,
					         CagdRType Param,
						 MvarZeroPrblmStruct *Problem)
{
    if (Sol1 == NULL) 
	return Sol2;
    if (Sol2 == NULL)
	return Sol1;

#ifdef DEBUG
    assert(Sol1 -> ActiveRepresentation == MVAR_ZER_SLVR_SOLUTION_PT_LIST &&
	   Sol2 -> ActiveRepresentation == MVAR_ZER_SLVR_SOLUTION_PT_LIST);
#endif /* DEBUG */

    /* Both are not NULL. Unite them in place, that is: use the structure   */
    /* of Sol1, and only the U.Pt from Sol2. This explains the deallocation */
    /* logic that is typical after returning from this routine.		    */
    Sol1 -> U.Pt = (MvarPtStruct *) CagdListAppend(Sol1 -> U.Pt, Sol2 -> U.Pt);
    MvarZeroSolverSolutionFree(Sol2, FALSE);
    return Sol1;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a 0D solution of a zero finding MVs representation problem,        *
* filters the points that fail the inequality constraints or appear more     *
* than once (up to numeric tolerance), and sorts the point list.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol:	The solution to be organized.				     *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The new, organized solution.		     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol0DMVs(
						 MvarZeroSolutionStruct *Sol,
						 MvarZeroPrblmStruct *Problem)
{
    if (Sol != NULL) {
	int NumOfMVs = Problem -> NumOfConstraints;
	CagdRType 
	    NumericTol = Problem -> NumericTol;
	MvarPtStruct 
	    *ZeroSet = Sol -> U.Pt;
	MvarConstraintType
	    *Constraints = Problem -> Constraints;
	MvarMVStruct 
	    **LclMVs = Problem -> U.MVs;

#   ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugMvarFilteringZeroSet, FALSE) {
		MvarPtStruct *Pt;

		for (Pt = ZeroSet; Pt != NULL; Pt = Pt -> Pnext) {
		    int i;

		    IRIT_INFO_MSG_PRINTF("\n\tUnfiltered (Err = %.15f) = ",
			           AttrGetRealAttrib(Pt -> Attr, "RngError"));
		    for (i = 0; i < Pt -> Dim; i++)
			IRIT_INFO_MSG_PRINTF("%.14f  ", Pt -> Pt[i]);
		}
		IRIT_INFO_MSG("\nDone\n"); 
	    }
	}
#   endif /* DEBUG */

	/* Filter out points that fail on the inequality constraints or     */
	/* identical points in the zero set.			 	    */
	ZeroSet = MvarZeroFilterSolutionSet(ZeroSet,
					    (const MvarMVStruct * const *)
					                                LclMVs,
					    Constraints, NumOfMVs,
					    sqrt(fabs(NumericTol)) * 10.0,
					    FALSE, TRUE, FALSE);
#   ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugMvarFilteringZeroSet, FALSE) {
		MvarPtStruct *Pt;

		for (Pt = ZeroSet; Pt != NULL; Pt = Pt -> Pnext) {
		    int i;

		    IRIT_INFO_MSG("\n\tFiltered = ");
		    for (i = 0; i < Pt -> Dim; i++)
			IRIT_INFO_MSG_PRINTF("%.14f  ", Pt -> Pt[i]);
		}
	    }
	}
#   endif /* DEBUG */

	/* Sort the result based on first axis. */
	Sol -> U.Pt = ZeroSet;
	return Sol;
    }
    else
	return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the action taken by the 0D solver when reaching a sub-domain of     M
* size less than subdivision tolerance.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   SubdivTolAction:   New setting for the action type.                      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroSubdivTolActionType:       Old setting for the action type.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver, MvarZeroSolverSinglrSol0DMVs			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZerosSubdivTolAction                                                 M
*****************************************************************************/
MvarZeroSubdivTolActionType MvarZerosSubdivTolAction(
				 MvarZeroSubdivTolActionType SubdivTolAction)
{
    MvarZeroSubdivTolActionType	
	OldVal = _MVGlblUponSubdivTol;

    _MVGlblUponSubdivTol = SubdivTolAction;

    return OldVal;
}
