/******************************************************************************
* ZrET0D.c - tools to compute zero sets of multivariates expression trees     *
*            representation and the expected solution set is a                *
*	     zero-manifold: a finite set of points.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber and Yoni Mizrahi, Feb. 13.			      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"

IRIT_GLOBAL_DATA CagdBType 
    _MVGlblZeroETUseCommonExpr = TRUE,
    _MVGlblZeroETCnvrtBzrETs2MVs = TRUE;

static CagdBType MvarZeroSolverNoZeroTest0DExpTr(MvarZeroPrblmStruct *Problem,
						 int Ind);
static void MvarZeroSolveBoundary0DExpTr(MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverTopoTest0DExpTr(MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep0DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroPrblmStruct **MvarZeroSubdivProblem0DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol0DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions0DExpTr(
						 MvarZeroSolutionStruct *Sol1,
						 MvarZeroSolutionStruct *Sol2, 
					         MvarMVDirType Dir,
					         CagdRType Param,
						 MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol0DExpTr(
						   MvarZeroSolutionStruct *Sol,
						   MvarZeroPrblmStruct 
						                     *Problem);
IRIT_STATIC_DATA const MvarZeroSolverCallBackFcnStruct 
    ZrET0DCB = {
        MvarZeroSolveBoundary0DExpTr,
        MvarZeroSolverNoZeroTest0DExpTr,
        MvarZeroSolverTopoTest0DExpTr,
	MvarZeroSubdivProblem0DExpTr,
	MvarZeroSolverNumericStep0DExpTr,
	MvarZeroSolverSinglrSol0DExpTr,
	MvarZeroSolverUniteSolutions0DExpTr,
	MvarZeroSolverOrganizeSol0DExpTr,
	MvarZeroUpdateProblemDmnExpTr,
	MvarZeroFirstSmoothUpdatesExpTr
    };


static void MvarExprTreeEqnsUpdateCommonExprIdcs(
					     MvarExprTreeStruct *ET,
					     MvarExprTreeStruct **CommonExprs,
					     int NumOfCommonExpr);
static int MvarExprTreeToVector(MvarExprTreeStruct *ET,
				MvarExprTreeStruct **Vec,
				int Idx);
static void MvarExprTreeMarkUsedNodes(MvarExprTreeStruct *Expr);
static int MvarExprTreeEqnsSubdivAtParam(const MvarExprTreeEqnsStruct *Eqns,
					 CagdRType t,
					 MvarMVDirType Dir,
					 MvarExprTreeEqnsStruct **Eqns1,
					 MvarExprTreeEqnsStruct **Eqns2);
static MvarPtStruct *MvarExprTreeEqnsZeroFilter0DSolutionSet(
					   MvarPtStruct *ZeroSet,
					   const MvarExprTreeEqnsStruct *Eqns,
					   CagdRType Tol);

#ifdef DEBUG
static void MVDbgPrintfStdout(const char *Str);
#endif /* DEBUG */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the callback functions for the Expression Trees representation, 0D  M
* solution case.							     M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*		                                                             *
* RETURN VALUE: 							     M
*   void								     M
*									     *
* SEE ALSO:								     M
*									     M
*									     *
* KEYWORDS:                                                                  M
*   MvarZeroSolverSetCallbackFcns0DExpTr				     M
*****************************************************************************/
void MvarZeroSolverSetCallbackFcns0DExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    IRIT_GEN_COPY(&PrblmIntrnl -> CallbackFunctions, &ZrET0DCB,
		  sizeof(MvarZeroSolverCallBackFcnStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Copy the given expression trees and process and fetch common expressions M
* out to a separated common expressions' vector, all within the returned     M
* expression tree equations structure.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVETs:      Input mvar expression tree equations.                        M
*   NumOfMVETs: Number of input mvar expression tree equations, in MVETs.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarExprTreeEqnsStruct *:   Build set of equations with common exprs.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroOrganizeETs0DProblem		                             M
*****************************************************************************/
MvarExprTreeEqnsStruct *MvarZeroOrganizeETs0DProblem(
				      const MvarExprTreeStruct * const *MVETs,
				      int NumOfMVETs)
{
    int i, j, k, n, Res,
	MaxNumCommonExprs = _MVGlblZeroETUseCommonExpr ? NumOfMVETs : 0;
    MvarExprTreeStruct *Expr, **Exprs;
    MvarExprTreeEqnsStruct
	*Eqns = MvarExprTreeEqnsMalloc(NumOfMVETs, MaxNumCommonExprs);

#ifdef DEBUG_MVAR_ZERO_REPORT_ALL_ETS
    printf("\n\nEquations:\n");
    for (i = 0; i < NumOfMVETs; i++) {
	printf("\nEqn %d: ", i);
	MvarExprTreePrintInfo(MVETs[i], FALSE, TRUE, MVDbgPrintfStdout);
    }
    printf("\n");
#endif /* DEBUG_MVAR_ZERO_REPORT_ALL_ETS */

    /* Copy the constraints. */
    for (i = n = 0; i < NumOfMVETs; i++) {
        Eqns -> Eqns[i] = MvarExprTreeCopy(MVETs[i], FALSE, TRUE);
	n += MvarExprTreeSize(Eqns -> Eqns[i]);

	/* Make sure we only have Bsp leaves (and no Bzr leaves.) */
        Res = MvarExprTreeCnvrtBzr2BspMV(Eqns -> Eqns[i]);
        assert(Res);
    }

    if (!_MVGlblZeroETUseCommonExpr)
        return Eqns; 

    /* Build an auxiliary vector of all expressions in these expression     */
    /* trees, top level expressions first, and search in O(n^2) for all     */
    /* expressions that are the same, and mark them as such.	            */
    Exprs = (MvarExprTreeStruct **)
                                 IritMalloc(sizeof(MvarExprTreeStruct *) * n);
    for (i = k = 0; i < NumOfMVETs; i++)
        k = MvarExprTreeToVector(Eqns -> Eqns[i], Exprs, k);
    assert(k == n);

    /* Now search for similar expressions in O(n^2). */
    for (i = 0; i < n; i++) {
        Exprs[i] -> IAux = -1;
	Exprs[i] -> PAux = NULL;
    }

#   ifdef MVAR_NO_COMMON_EXPRESSION /* Define to test with no common exprs. */
    k = 0;
#   else
    for (i = k = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
	    if (Exprs[j] -> IAux > -1)
	        continue;

	    if (MvarExprTreesSame(Exprs[i], Exprs[j], IRIT_UEPS)) {
	        if (Exprs[i] -> IAux == -1)
		    Exprs[i] -> IAux = k++;

	        Exprs[j] -> IAux = Exprs[i] -> IAux;
	    }
	}
    }
#   endif /* MVAR_NO_COMMON_EXPRESSION */

    if (k > 0) {
        /* Having k common expressions - make sure we can hold them all. */
        MvarExprTreeEqnsReallocCommonExprs(Eqns, k);
	IRIT_ZAP_MEM(Eqns -> CommonExprs, k * sizeof(MvarExprTreeStruct *));

	/*  For every common expressions do:                                */
	/*  + Copy the common expression to the common expressions area.    */
	/*  + Remove expressions from original ETs and substitute with a    */
	/*    MVAR_ET_NODE_COMMON_EXPR NodeType that refers the common expr.*/
	for (i = 0; i < n; i++) {
	    int m;

	    /* skip if not a common expression or already processed. */
	    if (Exprs[i] -> IAux < 0 || Exprs[i] -> PAux != NULL)
	        continue;

	    /* If we are here, then this one is marked a common expression. */
	    j = Exprs[i] -> IAux;
	    assert(j < k);

	    if (Eqns -> CommonExprs[j] == NULL) {
	        /* Copy this expression to the CommonExprs vector - it is   */
	        /* the first time we encounter the expression in this loop. */
	        Eqns -> CommonExprs[j] = MvarExprTreeCopy(Exprs[i],
							  FALSE, TRUE);

		/* Cannot have a common expression as root node here. */
		assert(Eqns -> CommonExprs[j] -> NodeType !=
		                                    MVAR_ET_NODE_COMMON_EXPR);
	    }
	    else {
	        /* Other instances should have been removed by code below!  */
	        assert(Exprs[i] -> NodeType == MVAR_ET_NODE_COMMON_EXPR);
	    }

	    /* Remove all sub-expression of this one from vector, free      */
	    /* them, and make them a common expression reference node.      */
	    for (m = 0; m < n; m++) {
	        if (Exprs[m] -> IAux == j) {
		    Expr = Exprs[m];
		    MvarExprTreeFreeSlots(Expr, TRUE);
		    Expr -> Left = Eqns -> CommonExprs[j];
		    Expr -> Right = NULL;
		    Expr -> NodeType = MVAR_ET_NODE_COMMON_EXPR;
		    Expr -> PAux = Expr -> Left;      /* Mark as processed. */
		}
	    }
	}
	Eqns -> NumCommonExprs = k;

#	ifdef DEBUG
	    for (i = 0; i < k; i++)
	        assert(Eqns -> CommonExprs[i] != NULL);
#	endif /* DEBUG */
    }

    /* Time to free unused expressions, from the expressions vector. */
    for (i = 0; i < n ; i++)
	Exprs[i] -> IAux2 = 0;
    for (i = 0; i < NumOfMVETs; i++)
	MvarExprTreeMarkUsedNodes(Eqns -> Eqns[i]);
    for (i = 0; i < n ; i++)
	if (Exprs[i] -> IAux2 == 0)
	    MvarExprTreeFree(Exprs[i], TRUE);

    /* Note that not all k entries in CommonExprs will be valid and some can */
    /* be NULL as we could find a common sub expression in a common super    */
    /* expression, and we grab as common only the super expression.          */

#ifdef DEBUG_MVAR_REPORT_COMMON_EXPR
    for (i = 0; i < k; i++)
        assert(Eqns -> CommonExprs[i] != NULL);

    for (i = 0; i < n; i++) {
        if (Exprs[i] -> NodeType == MVAR_ET_NODE_COMMON_EXPR) {
	    assert(Exprs[i] -> Left != NULL && Exprs[i] -> Right == NULL);
	}
    }

    printf("\nFound %d common expression out of %d expressions in %d equations.\n",
	    k, n, NumOfMVETs);
    
    printf("Common Expressions:\n");
    for (i = 0; i < k; i++) {
        printf("\nCE %d (%08x): ", i, Eqns -> CommonExprs[i]);
	MvarExprTreePrintInfo(Eqns -> CommonExprs[i], TRUE, MVDbgPrintfStdout);
        printf("\nCE %d: ", i);
	MvarExprTreePrintInfo(Eqns -> CommonExprs[i], FALSE, MVDbgPrintfStdout);
    }
    printf("\n\nEquations:\n");
    for (i = 0; i < NumOfMVETs; i++) {
        printf("\nEqn %d: ", i);
	MvarExprTreePrintInfo(Eqns -> Eqns[i], TRUE, MVDbgPrintfStdout);
        printf("\nEqn %d: ", i);
	MvarExprTreePrintInfo(Eqns -> Eqns[i], FALSE, MVDbgPrintfStdout);
    }
    printf("\n");
#endif /* DEBUG_MVAR_REPORT_COMMON_EXPR */

    IritFree(Exprs);

    return Eqns;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocate a structure to hold NumEqns equations and at most               M
* MaxNumCommonExprs common expressions.                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   NumEqns:            Number of equations we have.                         M
*   MaxNumCommonExprs:  Maximum number of common expression we can initially M
*			hold.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarExprTreeEqnsStruct *:  Allocated structure.                          M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarExprTreeEqnsFree					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeEqnsMalloc                                                   M
*****************************************************************************/
MvarExprTreeEqnsStruct *MvarExprTreeEqnsMalloc(int NumEqns,
					       int MaxNumCommonExprs)
{
    MvarExprTreeEqnsStruct
        *Eqns = (MvarExprTreeEqnsStruct *)
				   IritMalloc(sizeof(MvarExprTreeEqnsStruct));

    /* Allocate the equations' slots. */
    Eqns -> Eqns = (MvarExprTreeStruct **)
		           IritMalloc(sizeof(MvarExprTreeStruct *) * NumEqns);
    IRIT_ZAP_MEM(Eqns -> Eqns, sizeof(MvarExprTreeStruct *) * NumEqns);
    Eqns -> NumEqns = NumEqns;
    Eqns -> NumZeroEqns = Eqns -> NumZeroSubdivEqns = -1;

    /* Allocate the common expressions' slots. */
    if (MaxNumCommonExprs > 0)
        Eqns -> CommonExprs = (MvarExprTreeStruct **)
		 IritMalloc(sizeof(MvarExprTreeStruct *) * MaxNumCommonExprs);
    else
        Eqns -> CommonExprs = NULL;
    Eqns -> NumCommonExprs = 0;
    Eqns -> MaxNumCommonExprs = MaxNumCommonExprs;

    /* Allocate the constraints' types. */
    Eqns -> ConstraintTypes = IritMalloc(sizeof(MvarConstraintType) * NumEqns);

    return Eqns;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Free all data allocated in the expression tree equation's structure.     M
*                                                                            *
* PARAMETERS:                                                                M
*   Eqns:   Data structure to free.                                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarExprTreeEqnsMalloc				     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeEqnsFree                                                     M
*****************************************************************************/
void MvarExprTreeEqnsFree(MvarExprTreeEqnsStruct *Eqns)
{
    int i;

    for (i = 0; i < Eqns -> NumCommonExprs; i++) {
        if (Eqns -> CommonExprs[i] != NULL)
	    MvarExprTreeFree(Eqns -> CommonExprs[i], FALSE);
    }

    if (Eqns -> CommonExprs != NULL)
        IritFree(Eqns -> CommonExprs);

    for (i = 0; i < Eqns -> NumEqns; i++)
        MvarExprTreeFree(Eqns -> Eqns[i], FALSE);
    IritFree(Eqns -> Eqns);

    IritFree(Eqns -> ConstraintTypes);

    IritFree(Eqns);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reallocate (increase) the number of common expressions give Eqns can     M
* hold, in place.						             M
*                                                                            *
* PARAMETERS:                                                                M
*   Eqns:      Set of equations to increase, in place, the number of coomon  M
*              expressions it can old.					     M
*   NewSize:   of vector of common expression, zero to double the size.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeEqnsReallocCommonExprs                                       M
*****************************************************************************/
void MvarExprTreeEqnsReallocCommonExprs(MvarExprTreeEqnsStruct *Eqns,
					int NewSize)
{
    if (NewSize == 0)
        NewSize = Eqns -> MaxNumCommonExprs * 2;
    else if (NewSize <= Eqns -> MaxNumCommonExprs)
        return;

    Eqns -> CommonExprs =
        IritRealloc(Eqns -> CommonExprs,
		    sizeof(MvarExprTreeStruct *) * Eqns -> MaxNumCommonExprs,
		    sizeof(MvarExprTreeStruct *) * NewSize);

    Eqns -> MaxNumCommonExprs = NewSize;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Traverse an expression tree and update every common expression index to  *
* a real reference to a (common) expression tree.                            *
*                                                                            *
* PARAMETERS:                                                                *
*   ET:         Expression tree to traverse and update common expresssions.  *
*   CommonExprs:      Vector of common expressions to reference from MVET.   *
*   NumOfCommonExpr:  Size of CommonExprs vector.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void MvarExprTreeEqnsUpdateCommonExprIdcs(
					     MvarExprTreeStruct *ET,
					     MvarExprTreeStruct **CommonExprs,
					     int NumOfCommonExpr)
{
    switch (ET -> NodeType) {
	case MVAR_ET_NODE_LEAF:
	    break;
	case MVAR_ET_NODE_ADD:
	case MVAR_ET_NODE_SUB:
	case MVAR_ET_NODE_MULT:
	case MVAR_ET_NODE_DOT_PROD:
	case MVAR_ET_NODE_CROSS_PROD:
	    MvarExprTreeEqnsUpdateCommonExprIdcs(ET -> Right, CommonExprs,
						 NumOfCommonExpr);
	case MVAR_ET_NODE_EXP:
	case MVAR_ET_NODE_LOG:
	case MVAR_ET_NODE_COS:
	case MVAR_ET_NODE_SQRT:
	case MVAR_ET_NODE_SQR:
	case MVAR_ET_NODE_NPOW:
	case MVAR_ET_NODE_RECIP:
	    MvarExprTreeEqnsUpdateCommonExprIdcs(ET -> Left, CommonExprs,
						 NumOfCommonExpr);
	    break;
	case MVAR_ET_NODE_COMMON_EXPR:
	    assert(ET -> IAux >= 0 && ET -> IAux < NumOfCommonExpr);
	    ET -> Left = CommonExprs[ET -> IAux];
	    break;
	default:
	    assert(0);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Place all nodes (top down) in one linear vector.		             *
*                                                                            *
* PARAMETERS:                                                                *
*   ET:            Expression tree to place all its node on Vec.	     *
*   Vec:           Where to place all nodes of ET, top down.		     *
*   Idx:           Index into vec where to place the next node.		     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:	   Last index of a placed node in vec.			     *
*****************************************************************************/
static int MvarExprTreeToVector(MvarExprTreeStruct *ET,
				MvarExprTreeStruct **Vec,
				int Idx)
{
    if (ET == NULL)
        return Idx;

    switch (ET -> NodeType) {
	case MVAR_ET_NODE_LEAF:
	    ET -> IAux2 = Idx;
	    Vec[Idx++] = ET;
	    return Idx;
	case MVAR_ET_NODE_ADD:
	case MVAR_ET_NODE_SUB:
	case MVAR_ET_NODE_MULT:
	case MVAR_ET_NODE_DOT_PROD:
	case MVAR_ET_NODE_CROSS_PROD:
	    Idx = MvarExprTreeToVector(ET -> Left, Vec, Idx);
	    Idx = MvarExprTreeToVector(ET -> Right, Vec, Idx);
	    ET -> IAux2 = Idx;
	    Vec[Idx++] = ET;  /* We want the smaller (sons) subtrees first. */
	    return Idx;
	case MVAR_ET_NODE_COMMON_EXPR:
	    /* Should not encounter COMMON EXPR nodes at this time. */
	default:
	    assert(0);
	    return 0;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Mark all nodes of this tree as 1 in the IAux2 slot.			     *
*                                                                            *
* PARAMETERS:                                                                *
*    Expr:     Root node of tree to mark all its nodes.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*    void                                                                    *
*****************************************************************************/
static void MvarExprTreeMarkUsedNodes(MvarExprTreeStruct *Expr)
{
    Expr -> IAux2 = 1;

    switch (Expr -> NodeType) {
	case MVAR_ET_NODE_LEAF:
	case MVAR_ET_NODE_COMMON_EXPR:
	    break;
	case MVAR_ET_NODE_ADD:
	case MVAR_ET_NODE_SUB:
	case MVAR_ET_NODE_MULT:
	case MVAR_ET_NODE_DOT_PROD:
	case MVAR_ET_NODE_CROSS_PROD:
	    MvarExprTreeMarkUsedNodes(Expr -> Left);
	    MvarExprTreeMarkUsedNodes(Expr -> Right);
	    break;
	default:
	    assert(0);
	    break;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Testing if we can rule out the possibility of a solution for the         *
* problem in current domain, by inspecting constraint number Ind. Callback   *
* function for the ETs case, 0D solution. 				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*   Ind:	The index of the ET to be tested.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if domain can be purged (no zero can exist), FALSE      *
*		otherwise.						     M
*****************************************************************************/
static CagdBType MvarZeroSolverNoZeroTest0DExpTr(MvarZeroPrblmStruct *Problem,
						 int Ind)
{
    const MvarBBoxStruct
        *BBox = MvarExprTreeBBox(Problem -> U.Eqns -> Eqns[Ind]);

    assert(BBox -> Dim == 1);   /* Expects scalar value to be returned. */

    switch (Problem -> Constraints[Ind]) {
        default:
	    assert(0);
        case MVAR_CNSTRNT_ZERO:
        case MVAR_CNSTRNT_ZERO_SUBDIV:
	    return BBox -> Min[0] > 0.0 || BBox -> Max[0] < 0.0;
        case MVAR_CNSTRNT_POSITIVE:
	    return BBox -> Max[0] < 0.0;
        case MVAR_CNSTRNT_NEGATIVE:
	    return BBox -> Min[0] > 0.0;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Solving the problem on the boundary of the domain. De-facto, here in 0D  *
* solution space, we do nothing (the problem is already fully constrained).  *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*                                                                            *
* KEYWORDS:								     *
*   MvarZeroSolveBoundary0DExpTr					     *
*****************************************************************************/
static void MvarZeroSolveBoundary0DExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    PrblmIntrnl -> SolutionsOnBoundary = NULL;
    PrblmIntrnl -> NumOfBoundarySolutions = 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Topological guarantee in the Expression Trees, 0D solution case.	     *
*									     *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if at most a single solution in current domain, FALSE   *
*               otherwise.						     *
*****************************************************************************/
static CagdBType MvarZeroSolverTopoTest0DExpTr(MvarZeroPrblmStruct *Problem)
{
    /* Normal cones for extression trees needs some revision as they are   */
    /* not complete (Dot/Cross products are missing for one.               */
    return FALSE;
    /* !MvarExprTreeConesOverlap(Problem -> U.Eqns); GERSHON */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   The recursive stage of the subdivision solver, invoked by the general    *
* solver, in the 0D, ETs case.   Given a problem, returns a vector of two    *
* subdivided problems, allocated dynamically.  The subdivided problems can   *
* be NULL.								     *
*									     *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroPrblmStruct **:  The array of two sub-problems, or NULL if can   *
*                            converts to MVs-only problem.		     *
*****************************************************************************/
static MvarZeroPrblmStruct **MvarZeroSubdivProblem0DExpTr(
					         MvarZeroPrblmStruct *Problem)
{
    int i, Dir;
    CagdRType Knot;
    MvarExprTreeEqnsStruct *Eqns1, *Eqns2,
        *Eqns = Problem -> U.Eqns;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    MvarZeroPrblmStruct
	**SubProblems = (MvarZeroPrblmStruct **)
		                 IritMalloc(2 * sizeof(MvarZeroPrblmStruct *));
    			
    for (i = 0; i < Eqns -> NumEqns; i++) {
        if ((Dir = MvarExprTreeInteriorKnots(Eqns -> Eqns[i], &Knot)) >= 0)
	    break;
    }

    if (i < Eqns -> NumEqns) {
        /* i'th Eqn has an interior knot in direction Dir - subdivide. */
        MvarExprTreeEqnsSubdivAtParam(Eqns, Knot, Dir, &Eqns1, &Eqns2);
	SubProblems[0] = MvarZeroSolverSubProblem(Problem, NULL, Eqns1, NULL);
	SubProblems[1] = MvarZeroSolverSubProblem(Problem, NULL, Eqns2, NULL);
    }
    else if (_MVGlblZeroETCnvrtBzrETs2MVs) {
        IritFree(SubProblems);    /* Caller should convert problem to MVs. */
        return NULL;
    }
    else {
	int i, j;
	CagdRType
	    Min = 0.0,
	    Max = 0.0;

        /* Find largest direction. */
        Dir = -1;
	for (i = 0; i < Eqns -> NumEqns; i++) {
	    CagdBType HasDmn;
	    CagdRType MinCurr, MaxCurr;

	    for (j = 0; j < Eqns -> NumEqns; j++)
	        if ((HasDmn = MvarETDomain(Eqns -> Eqns[j], &MinCurr, &MaxCurr,
					   i)) == TRUE)
		    break;
	    assert(HasDmn);

	    if (MaxCurr - MinCurr > Max - Min) {
	        Dir = i;
		Min = MinCurr;
		Max = MaxCurr;
	    }
	}

	if (Max - Min < Problem -> SubdivTol) {
	    /* Domain is small enough - should not get here in the 1st      */
	    /* place as the caller should have recognized domain too small. */
	    assert(0);
	    SubProblems[0] = SubProblems[1] = NULL;
	}
	else {
	    CagdRType t;

	    if (PrblmIntrnl -> ParamPerturb < (Max - Min) / 10.0)
	        t = (Min + Max) * 0.5 + PrblmIntrnl -> ParamPerturb;
	    else
	        t = (Min + Max) * 0.5;

	    MvarExprTreeEqnsSubdivAtParam(Eqns, t, Dir, &Eqns1, &Eqns2);
	    SubProblems[0] = MvarZeroSolverSubProblem(Problem, NULL, Eqns1,
						      NULL);
	    SubProblems[1] = MvarZeroSolverSubProblem(Problem, NULL, Eqns2,
						      NULL);
	}
    }

    /* EMPTY DEFINITION - DEFAULT ONLY */
    return SubProblems;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* The numeric improvement for a single point in a domain of the zero finding *
* problem, for the ETs, 0D solution case.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution - a single point or NULL.	     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep0DExpTr(
					        MvarZeroPrblmStruct *Problem)
{
    MvarPtStruct *ZeroPt;

    ZeroPt = MvarZeroGenPtMidDmn(Problem, TRUE);
    ZeroPt = MvarZero0DNumeric(ZeroPt, Problem -> U.Eqns, NULL,
			       Problem -> NumOfZeroConstraints,
			       Problem -> NumericTol, NULL, NULL);

    if (ZeroPt) {
	return MvarZeroSolverSolutionNew(NULL, NULL, ZeroPt, 
					 MVAR_ZER_SLVR_SOLUTION_PT_LIST);
    }
    else
        return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a set of equations over some domain, divided them all in Dir at t. *
*                                                                            *
* PARAMETERS:                                                                *
*   Eqns:    The equations to subdivide.                                     *
*   t:       Parameter to subdivide at.                                      *
*   Dir:     Direction of subdivision.                                       *
*   Eqns1:   First set of divided equations.				     *
*   Eqns2:   Second set of divided equations.                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    TRUE if successful, FALSE otherwise.                             *
*****************************************************************************/
static int MvarExprTreeEqnsSubdivAtParam(const MvarExprTreeEqnsStruct *Eqns,
					 CagdRType t,
					 MvarMVDirType Dir,
					 MvarExprTreeEqnsStruct **Eqns1,
					 MvarExprTreeEqnsStruct **Eqns2)
{
    int i;

    *Eqns1 = MvarExprTreeEqnsMalloc(Eqns -> NumEqns,
				    Eqns -> MaxNumCommonExprs);
    *Eqns2 = MvarExprTreeEqnsMalloc(Eqns -> NumEqns,
				    Eqns -> MaxNumCommonExprs);

    /* Subdivide the common expressions first, if any. */
    for (i = 0; i < Eqns -> NumCommonExprs; i++) {
        if (Eqns -> CommonExprs[i] == NULL) {
	    (*Eqns1) -> CommonExprs[i] = (*Eqns2) -> CommonExprs[i] = NULL;
	}
	else {
#ifdef DEBUG
	    int RetVal = 
#endif /* DEBUG */
	        MvarExprTreeSubdivAtParam(Eqns -> CommonExprs[i],
					  t, Dir,
					  &(*Eqns1) -> CommonExprs[i],
					  &(*Eqns2) -> CommonExprs[i]);

#ifdef DEBUG
	    assert(RetVal);
#endif /* DEBUG */

	    /* No need to go over the trees if no common eqns, recursively. */
	    if (Eqns -> NumCommonExprs > 0) {
	        MvarExprTreeEqnsUpdateCommonExprIdcs(
						 (*Eqns1) -> CommonExprs[i],
						 (*Eqns1) -> CommonExprs,
						 Eqns -> NumCommonExprs);
		MvarExprTreeEqnsUpdateCommonExprIdcs(
						 (*Eqns2) -> CommonExprs[i],
						 (*Eqns2) -> CommonExprs,
						 Eqns -> NumCommonExprs);
	    }
	}
    }

    /* Subdivide the equations themselves now. */
    for (i = 0; i < Eqns -> NumEqns; i++) {
#ifdef DEBUG
	    int RetVal = 
#endif /* DEBUG */
	        MvarExprTreeSubdivAtParam(Eqns -> Eqns[i], t, Dir,
					  &(*Eqns1) -> Eqns[i],
					  &(*Eqns2) -> Eqns[i]);

#ifdef DEBUG
	assert(RetVal);
#endif /* DEBUG */

        /* No need to go over the trees if no common equations.*/
        if (Eqns -> NumCommonExprs > 0) {
	    MvarExprTreeEqnsUpdateCommonExprIdcs((*Eqns1) -> Eqns[i],
						 (*Eqns1) -> CommonExprs,
						 Eqns -> NumCommonExprs);
	    MvarExprTreeEqnsUpdateCommonExprIdcs((*Eqns2) -> Eqns[i],
						 (*Eqns2) -> CommonExprs,
						 Eqns -> NumCommonExprs);
        }
    }

    (*Eqns1) -> NumCommonExprs =
        (*Eqns2) -> NumCommonExprs = Eqns -> NumCommonExprs;

    IRIT_GEN_COPY((*Eqns1) -> ConstraintTypes, Eqns -> ConstraintTypes,
	     sizeof(MvarConstraintType) * Eqns -> NumEqns);
    IRIT_GEN_COPY((*Eqns2) -> ConstraintTypes, Eqns -> ConstraintTypes,
	     sizeof(MvarConstraintType) * Eqns -> NumEqns);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Returning a solution in the case of reaching subdivision tolerance, in the *
* ETs, 0D case.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution.			             *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol0DExpTr(
					         MvarZeroPrblmStruct *Problem)
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
	    Solution = MvarZeroSolverNumericStep0DExpTr(Problem);
	    break;
	case MVAR_ZERO_NUMERIC_STEP_VERIFIED:
	    if ((Solution = MvarZeroSolverNumericStep0DExpTr(Problem)) != NULL) {
	        /* With the Assumption the MVs are defined... */
		Err = MvarMVEvalErrorL1((MvarMVStruct const * const *)
				                            Problem -> U.MVs,
					Solution -> U.Pt -> Pt, 
	        		        PrblmIntrnl -> NumOfZeroSubdivConstraints);
		if (Err > Problem -> NumericTol) {
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
*   Unites two 0D solutions of a zero finding ETs representation problem to  *
* a new solution, in place.						     *
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
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions0DExpTr(
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
*   Given a 0D solution of a zero finding ETs representation problem,        *
* filters close points and points that do not satisfy inequality             *
* constraints.								     *
*									     *
* PARAMETERS:                                                                *
*   Sol:	The solution to be organized.				     *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct:  The new, organized solution.		     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol0DExpTr(
						 MvarZeroSolutionStruct *Sol,
						 MvarZeroPrblmStruct *Problem)
{

    if (Sol != NULL) {
	CagdRType 
	    NumericTol = Problem -> NumericTol;
	MvarPtStruct 
	    *ZeroSet = Sol -> U.Pt;
	const MvarExprTreeEqnsStruct 
	    *Eqns = Problem -> U.Eqns;

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
	ZeroSet = MvarExprTreeEqnsZeroFilter0DSolutionSet(ZeroSet, Eqns,
					       sqrt(fabs(NumericTol)) * 10.0);
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
	ZeroSet = MvarPtSortListAxis(ZeroSet, 1);
	Sol -> U.Pt = ZeroSet;
	return Sol;
    }
    else
	return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Filters out identical points or points that fail the inequality          *
* constraints.                                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   ZeroSet:      The solution points to filter out.                         *
*   Eqns:         Multivariate ET constraints.                               *
*   Tol:          Tolerance to consider to points same in L^1 norm.          *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:   Filtered solution points.                              *
*****************************************************************************/
static MvarPtStruct *MvarExprTreeEqnsZeroFilter0DSolutionSet(
					   MvarPtStruct *ZeroSet,
					   const MvarExprTreeEqnsStruct *Eqns,
					   CagdRType Tol)
{
    int i, l;
    MvarPtStruct *Pt,
	*OutSet = NULL;
    MvarExprTreeStruct
        * const *ETs = Eqns -> Eqns;
    int NumEqns = Eqns -> NumEqns,
	NumZeroEqns = Eqns -> NumZeroEqns;

    /* Sort the result based on first axis. */
    ZeroSet = MvarPtSortListAxis(ZeroSet, 1);

    while (ZeroSet != NULL) {
        IRIT_LIST_POP(Pt, ZeroSet);

	if (AttrGetIntAttrib(Pt -> Attr, "Similar") == TRUE) {
	    MvarPtFree(Pt);
	}
	else {
	    MvarPtStruct *Pt2;
	    CagdBType
	        PurgePt = FALSE;

	    /* Lets see if we fail any inequality constraint. */
	    for (i = NumZeroEqns; i < NumEqns && !PurgePt; i++) {
                int NumOfCoord = ETs[i] -> PtSize;
		CagdRType
                    *R = MvarExprTreeEval(ETs[i], Pt -> Pt);

		switch (Eqns -> ConstraintTypes[i]) {
		    case MVAR_CNSTRNT_POSITIVE:
			if (NumOfCoord > 1) { /* Union implementation. */
			    PurgePt = TRUE;
			    for (l = 1; l <= NumOfCoord; l++) {
			        if (R[l] >= 0.0) {
				    PurgePt = FALSE;
				    break;
				}
			    }
			}
			else {
			    if (R[1] < 0.0)
			        PurgePt = TRUE;
			}
			break;
		    case MVAR_CNSTRNT_NEGATIVE:
		        if (NumOfCoord > 1) { /* Union implementation. */
			    PurgePt = TRUE;
			    for (l = 1; l <= NumOfCoord; l++) {
			        if (R[l] <= 0.0) {
				    PurgePt = FALSE;
				    break;
				}
			    }
			}
			else {
			    if (R[1] > 0.0)
			        PurgePt = TRUE;
			}
			break;
	            default:
		        break;
		}
	    }

	    if (PurgePt) {
	        MvarPtFree(Pt);
	    }
	    else {
	        for (Pt2 = ZeroSet; Pt2 != NULL; Pt2 = Pt2 -> Pnext) {
		    for (i = 0; i < Pt -> Dim; i++) {
		        if (!IRIT_APX_EQ_EPS(Pt -> Pt[i], Pt2 -> Pt[i], Tol))
			    break;
		    }

		    if (i >= Pt -> Dim) {
		        /* Pt and Pt2 are same - mark Pt2 as similar. */
		        AttrSetIntAttrib(&Pt2 -> Attr, "Similar", TRUE);
		    }
		}

		IRIT_LIST_PUSH(Pt, OutSet);
	    }
	}
    }

    return OutSet;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the exploitation of common expression extraction in expression      M
* trees.  If TRUE, the ETs are scanned for common expressions that are then  M
* processed once only, during the subdivision process.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   UseCommonExpr:   TRUE to use common expressions, FALSE otherwise.        M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting of common expressions' use.                       M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarExprTreesZeros, MvarExprTreeZerosCnvrtBezier2MVs		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeZerosUseCommonExpr                                           M
*****************************************************************************/
int MvarExprTreeZerosUseCommonExpr(int UseCommonExpr)
{
    int OldVal = _MVGlblZeroETUseCommonExpr;

    _MVGlblZeroETUseCommonExpr = UseCommonExpr;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the way expression trees of Bezier MVs are treated.  If TRUE, the   M
* ETs are converted into MVs and the regular MV zero solver is invoked.  If  M
* FALSE, the ETs are subdivided all the way to SUbdivTol.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Bezier2MVs:   TRUE to convert to MVs, FALSE to subdivide Bezier ETs.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting for Bezier conversion setting.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarExprTreesZeros, MvarExprTreeZerosUseCommonExpr			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeZerosCnvrtBezier2MVs                                         M
*****************************************************************************/
int MvarExprTreeZerosCnvrtBezier2MVs(int Bezier2MVs)
{
    int OldVal = _MVGlblZeroETCnvrtBzrETs2MVs;

    _MVGlblZeroETCnvrtBzrETs2MVs = Bezier2MVs;

    return OldVal;
}

#ifdef DEBUG

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Debug function to print a string to stdout.                              *
*                                                                            *
* PARAMETERS:                                                                *
*   const char *:  String to print.                                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void MVDbgPrintfStdout(const char *Str)
{
    printf(Str);
}

#endif /* DEBUG */
