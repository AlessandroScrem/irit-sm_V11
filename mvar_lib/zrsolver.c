/******************************************************************************
* ZrSolver.c - main entrance to tools and interface for the framework of the  *
*	       generic solver, for zero set finding problems of various       *
*              dimensions and representations.		                      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, Feb. 13.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/geom_lib.h"

#ifdef DEBUG
IRIT_GLOBAL_DATA int
    _DebugZeroSolverDmnNum = 0,
    _DebugZeroSolverBySolDim[3] = { FALSE, FALSE, FALSE },
    _DebugZeroZeroErr = FALSE;
#endif /* DEBUG */

IRIT_GLOBAL_DATA CagdRType
    _MVGlblZeroDmnExtension = 0.0;

static int MvarZeroVerifyConstraintTypes(const MvarConstraintType *Constraints,
					 int NumOfMVs,
					 int *NumOfZeroMVs,
					 int *NumOfZeroSubdivMVs,
					 int *NumOfNonZeroMVs);
static void MvarZeroVerifyBspConstraints(MvarMVStruct **MVs,
				         int NumOfMVs);
static MvarZeroPrblmStruct **MvarZeroC1DiscontSubdiv(
						 MvarZeroPrblmStruct *Problem);

static MvarZeroSolutionStruct *MvarZeroETs2MVs(const MvarExprTreeEqnsStruct
					                                *Eqns,
					       CagdRType SubdivTol,
					       CagdRType NumericTol);

#ifdef DEBUG
static void MvarZeroDBGDmn(MvarZeroPrblmStruct *Problem);
#endif /* DEBUG */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Interface function for a generic MV equation solver, 0D solutions:       M
* constructs the generic problem structure, calls the solver and extracts    M
* the solution point list from the generic solution structure.               M
*   The set of NumOfMVs constraints may consist of equality or inequality    M
* constraints, as prescribed by the constraints vector. All multivariates    M
* are assumed to be in the same parametric domain size and dimension.        M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*                 Can be NULL in which case all constraints are equality.    M
*   NumOfMVs:     Size of the MVs and Constraints vector.                    M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   NumericTol:   Numeric tolerance of the numeric stage.  Measured in the   M
*		  image space of the constraints function.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   List of points on the solution set.  Dimension of the  M
*		      points will be the same as the dimensions of all MVs.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZerosNormalConeTest, MvarMVsZerosDomainReduction,                 M
*   MvarMVsZerosVerifier, MvarETsZeros0D, MvarMVsZeros1D, MvarMVsZeros2D     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZeros0D                                                           M
*****************************************************************************/
MvarPtStruct *MvarMVsZeros0D(MvarMVStruct * const *MVs,
			     MvarConstraintType *Constraints,
			     int NumOfMVs,
			     CagdRType SubdivTol,
			     CagdRType NumericTol)
{
    CagdBType FreeConstraints;
    MvarZeroPrblmStruct *Problem;

    /* If we need to allocate the constraints (all equality) - do it now. */
    if (Constraints == NULL) {
	int i;
        Constraints = (MvarConstraintType *)
	                   IritMalloc(sizeof(MvarConstraintType) * NumOfMVs);
        for (i = 0; i < NumOfMVs; i++)
	    Constraints[i] = MVAR_CNSTRNT_ZERO;

        FreeConstraints = TRUE;
    }
    else
        FreeConstraints = FALSE;

    if ((Problem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *) MVs,
					  NULL, NumOfMVs, Constraints,
					  SubdivTol, NumericTol,
					  IRIT_INFNTY, FALSE)) != NULL) {
        MvarZeroPrblmIntrnlStruct
	    *PrblmIntrnl = Problem -> _Internal;
	MvarPtStruct
	    *ZeroSet = NULL;
	MvarZeroSolutionStruct *Solution;

	switch (PrblmIntrnl -> ExpectedSolutionDim) {
	    case 0:
		if ((Solution = MvarZeroSolver(Problem)) != NULL) {
		    assert(Solution -> ActiveRepresentation == 
					      MVAR_ZER_SLVR_SOLUTION_PT_LIST);
		    ZeroSet = Solution -> U.Pt;
		}
		MvarZeroSolverSolutionFree(Solution, FALSE);
		break;
	    default:
	        assert(0);
		break;
	}

	MvarZeroSolverPrblmFree(Problem);

	if (FreeConstraints)
	    IritFree(Constraints);

	return ZeroSet;
    }

    if (FreeConstraints)
        IritFree(Constraints);

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Interface function for a generic MV equation solver, 2D solutions, in    M
* the temporary case where the problem is solved using the 0D solver and     M
* then surface fitting to the solution points.                               M
*   Constructs the generic problem structure, calls the solver and extracts  M
* the solution point list from the generic solution structure. The surface   M
* fitting step is done by the specific problem algorithm outside the solver  M
* environment.                                                               M
*   The set of NumOfMVs constraints may consist of equality or inequality    M
* constraints, as prescribed by the constraints vector. All multivariates    M
* are assumed to be in the same parametric domain size and dimension.        M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*                 Can be NULL in which case all constraints are equality.    M
*   NumOfMVs:     Size of the MVs and Constraints vector.                    M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   NumericTol:   Numeric tolerance of the numeric stage.  Measured in the   M
*		  image space of the constraints function.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   List of points on the solution set.  Dimension of the  M
*		      points will be the same as the dimensions of all MVs.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCrvSrfBisectorApprox                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZeros2DBy0D                                                       M
*****************************************************************************/
MvarPtStruct *MvarMVsZeros2DBy0D(MvarMVStruct * const *MVs,
			         MvarConstraintType *Constraints,
			         int NumOfMVs,
			         CagdRType SubdivTol,
			         CagdRType NumericTol)
{
    CagdBType FreeConstraints;
    MvarZeroPrblmStruct *Problem;

    /* If we need to allocate the constraints (all equality) - do it now. */
    if (Constraints == NULL) {
	int i;

	Constraints = (MvarConstraintType *)
			  IritMalloc(sizeof(MvarConstraintType) * NumOfMVs);
	for (i = 0; i < NumOfMVs; i++)
	    Constraints[i] = MVAR_CNSTRNT_ZERO;

	FreeConstraints = TRUE;
    }
    else
	FreeConstraints = FALSE;

    if ((Problem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *) MVs,
					  NULL, NumOfMVs, Constraints,
					  SubdivTol, NumericTol,
					  IRIT_INFNTY, TRUE)) != NULL) {
        MvarZeroPrblmIntrnlStruct
	    *PrblmIntrnl = Problem -> _Internal;
	MvarPtStruct
	    *ZeroSet = NULL;
	MvarZeroSolutionStruct *Solution;

	switch (PrblmIntrnl -> ExpectedSolutionDim) {
	    case 0:
		if ((Solution = MvarZeroSolver(Problem)) != NULL) {
		    assert(Solution -> ActiveRepresentation == 
			                      MVAR_ZER_SLVR_SOLUTION_PT_LIST);
		    ZeroSet = Solution -> U.Pt;
		}
		MvarZeroSolverSolutionFree(Solution, FALSE);
		break;
	    default:
		assert(0);
		break;
	}

	MvarZeroSolverPrblmFree(Problem);

	if (FreeConstraints)
	    IritFree(Constraints);

	return ZeroSet;
    }

    if (FreeConstraints)
	IritFree(Constraints);

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Interface function for a generic MV equation solver, 1D solutions:       M
* constructs the generic problem structure, calls the solver and extracts    M
* the list of solution polylines from the generic solution structure.        M
*   The set of NumOfMVs constraints may consist of equality or inequality    M
* constraints, as prescribed by the constraints vector. All multivariates    M
* are assumed to be in the same parametric domain size and dimension.        M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*                 Can be NULL in which case all constraints are equality.    M
*   NumOfMVs:     Size of the MVs and Constraints vector.                    M
*   Step:         Step size to use in the numeric tracing.		     M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   NumericTol:   Numeric tolerance of the numeric stage.  Measured in the   M
*		  image space of the constraints function.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: The list of polylines which approximate the curve. M
*			  Each polyline corresponds to the topologically     M
*			  isolated component of the curve and is in R^k, the M
*                         unioned parametric spaces of all input MVs.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSrfSrfInter, MvarMVsZeros1DMergeSingularPts, MvarMVsZeros0D          M
*   MvarETsZeros0D, MvarMVsZeros2D, MvarMVsZeros1DOneTrace		     M
*									     *
* KEYWORDS:                                                                  M
*   MvarMVsZeros1D			                                     M
*****************************************************************************/
MvarPolylineStruct *MvarMVsZeros1D(MvarMVStruct * const *MVs,
			           MvarConstraintType *Constraints,
			           int NumOfMVs,
			           CagdRType Step,
			           CagdRType SubdivTol,
			           CagdRType NumericTol)
{
    CagdBType FreeConstraints;
    int i;
    MvarPolylineStruct 
	*IntrCrv = NULL;
    MvarZeroPrblmStruct *Problem;
    MvarZeroSolutionStruct *Solution;

    /* If we need to allocate the constraints (all equality) - do it now. */
    if (Constraints == NULL) {
        Constraints = (MvarConstraintType *)
	                   IritMalloc(sizeof(MvarConstraintType) * NumOfMVs);
        for (i = 0; i < NumOfMVs; i++)
	    Constraints[i] = MVAR_CNSTRNT_ZERO;

        FreeConstraints = TRUE;
    }
    else
        FreeConstraints = FALSE;

    Problem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *) MVs,
				     NULL, NumOfMVs, Constraints, 
				     SubdivTol, NumericTol, Step, FALSE);
    Solution = MvarZeroSolver(Problem);
    if (Solution != NULL) {
	assert(Solution -> ActiveRepresentation == 
					      MVAR_ZER_SLVR_SOLUTION_POLYLINE);
	IntrCrv = Solution -> U.Pl;
    }
    MvarZeroSolverSolutionFree(Solution, FALSE);
    MvarZeroSolverPrblmFree(Problem);

    if (FreeConstraints)
        IritFree(Constraints);

    return IntrCrv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Interface function similar to MvarMVsZeros1D but only for tracing a      M
* single univariate component between StartEndPts.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*                 Can be NULL in which case all constraints are equality.    M
*   NumOfMVs:     Size of the MVs and Constraints vector.                    M
*   StartEndPts:  Start/end points for the polyline solution to trace.       M
*   Step:         Step size to use in the numeric tracing.		     M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   NumericTol:   Numeric tolerance of the numeric stage.  Measured in the   M
*		  image space of the constraints function.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: The traced polyline which approximates solution,   M
*			  in R^k.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSrfSrfInter, MvarMVsZeros1DMergeSingularPts, MvarMVsZeros0D          M
*   MvarETsZeros0D, MvarMVsZeros2D,  MvarMVsZeros1D			     M
*									     *
* KEYWORDS:                                                                  M
*   MvarMVsZeros1DOneTrace		                                     M
*****************************************************************************/
MvarPolylineStruct *MvarMVsZeros1DOneTrace(MvarMVStruct * const *MVs,
					   MvarConstraintType *Constraints,
					   int NumOfMVs,
					   MvarPtStruct *StartEndPts,
					   CagdRType Step,
					   CagdRType SubdivTol,
					   CagdRType NumericTol)
{
    CagdBType FreeConstraints;
    int i;
    MvarPolylineStruct 
	*IntrCrv = NULL;
    MvarZeroPrblmStruct *Problem;
    MvarZeroSolutionStruct *Solution;
    MvarZeroPrblmIntrnlStruct *PrblmIntrnl;

    /* If we need to allocate the constraints (all equality) - do it now. */
    if (Constraints == NULL) {
        Constraints = (MvarConstraintType *)
	                   IritMalloc(sizeof(MvarConstraintType) * NumOfMVs);
        for (i = 0; i < NumOfMVs; i++)
	    Constraints[i] = MVAR_CNSTRNT_ZERO;

        FreeConstraints = TRUE;
    }
    else
        FreeConstraints = FALSE;

    Problem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *) MVs,
				     NULL, NumOfMVs, Constraints, 
				     SubdivTol, NumericTol, Step, FALSE);

    PrblmIntrnl = Problem -> _Internal;

    PrblmIntrnl -> MVGradients = (MvarMVGradientStruct **) 
                             IritMalloc(sizeof(MvarMVGradientStruct *) *
					Problem -> NumOfZeroConstraints);
    for (i = 0; i < Problem -> NumOfZeroConstraints; i++) {
	PrblmIntrnl -> MVGradients[i] = 
	    MvarMVPrepGradient(Problem -> U.MVs[i], FALSE);
    }

    PrblmIntrnl -> ConstructionComplete = TRUE;

    PrblmIntrnl -> SolutionsOnBoundary =
        MvarZeroSolverSolutionNew(NULL, NULL, MvarPtCopyList(StartEndPts), 
				  MVAR_ZER_SLVR_SOLUTION_PT_LIST);

    PrblmIntrnl -> NumOfBoundarySolutions = 2;

    Solution = MVAR_ZERO_SLVR_APPLY(NumericImprovement)(Problem);
    if (Solution != NULL) {
	assert(Solution -> ActiveRepresentation == 
					      MVAR_ZER_SLVR_SOLUTION_POLYLINE);
	IntrCrv = Solution -> U.Pl;
    }
    MvarZeroSolverSolutionFree(Solution, FALSE);
    MvarZeroSolverPrblmFree(Problem);

    if (FreeConstraints)
        IritFree(Constraints);

    return IntrCrv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Interface function for a generic MV equation solver, 2D solutions:       M
* constructs the generic problem structure, calls the solver and extracts    M
* the list of either solution triangles or solution polylines from the       M
* generic solution structure.                                                M
*   The set of NumOfMVs constraints must consist of equality constraints,    M
* prescribed by the constraints vector (inequalities not supported). All     M
* multivariates are assumed to be in the same parametric domain size and     M
* dimension.							             M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:	    Vector of multivariate constraints.                      M
*   Constraints:    Either an equality or an inequality type of constraint.  M
*		    Can be NULL in which case all constraints are equality.  M
*   NumOfMVs:	    Size of the MVs and Constraints vector.                  M
*   Step:	    Step size to use in the numeric tracing of boundary      M
*		    (univariate) solutions.		                     M
*   SubdivTol:	    Tolerance of the subdivision process.  Tolerance is      M
*		    measured in the parametric space of the multivariates.   M
*   NumericTol:	    Numeric tolerance of the numeric stage.  Measured in the M
*		    image space of the constraints function.                 M
*   MapPt2EuclidSp: A pointer to the function used for mapping the           M
*		    solution from the problem's parameter space to the 3D    M
*		    Euclidean space. If not provided, the solution points    M
*                   are returned in the dimension (and semantics) of the     M
*                   original parameter space.                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolyStruct *:   The list of either triangles or polylines, according M
*			to the required output types, which approximate the  M
*			surface.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSrfSrfInter, MvarMVsZeros0D, MvarETsZeros0D, MvarMVsZeros1D	     M
*									     *
* KEYWORDS:                                                                  M
*   MvarMVsZeros2D			                                     M
*****************************************************************************/
MvarPolyStruct *MvarMVsZeros2D(MvarMVStruct * const *MVs,
			       MvarConstraintType *Constraints,
			       int NumOfMVs,
			       CagdRType Step,
			       CagdRType SubdivTol,
			       CagdRType NumericTol, 
			       MvarMapPrm2EucCallBackFuncType MapPt2EuclidSp)
{
    CagdBType FreeConstraints;
    int i;
    MvarTriangleStruct
	*TrList = NULL;
    MvarPolylineStruct
	*PolylineList = NULL;
    MvarPolyStruct 
	*Poly = NULL;
    MvarZeroPrblmStruct *Problem;
    MvarZeroSolutionStruct *Solution;

    /* If we need to allocate the constraints (all equality) - do it now. */
    if (Constraints == NULL) {
	Constraints = (MvarConstraintType *)
	    IritMalloc(sizeof(MvarConstraintType) * NumOfMVs);
	for (i = 0; i < NumOfMVs; i++)
	    Constraints[i] = MVAR_CNSTRNT_ZERO;

	FreeConstraints = TRUE;
    }
    else
	FreeConstraints = FALSE;

    Problem = MvarZeroSolverPrblmNew((const MvarMVStruct * const *) MVs,
				     NULL, NumOfMVs, Constraints, 
				     SubdivTol, NumericTol, Step, FALSE);
    if (MapPt2EuclidSp != NULL) {
	Problem -> _Internal -> CallbackFunctions.MapPtParamSpace2EuclidSpace =
	                                                        MapPt2EuclidSp;
    }
    Solution = MvarZeroSolver(Problem);
    if (Solution != NULL) {
	if (_MVGlblZeroOutputTypeLoops)
	    PolylineList = Solution -> U.Pl;
	else
	    TrList = Solution -> U.Tr;
	Poly = MvarPolyNew(PolylineList, TrList);
    }
    MvarZeroSolverSolutionFree(Solution, FALSE);
    MvarZeroSolverPrblmFree(Problem);

    if (FreeConstraints)
	IritFree(Constraints);
   return Poly;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Interface function for a MV expression trees solver, 1D solutions.	     M
*   Computes the simultaneous solution of the given set of NumOfMVs          M
* expression trees constraints.  A constraint can be equality or ineqaulity  M
* as prescribed by the Constraints vector. 				     M
*   All multivariates are assumed to be in the same parametric domain size   M
* and dimension.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVETs:	  A vector of MV expression trees constraints.  Should be of M
*                 size Dim where Dim is the dimension of the domain of MVs.  M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*                 Can be NULL in which case all constraints are equality.    M
*   NumOfMVETs:   Total number of MV expression trees constraints.           M
*   SubdivTol:	  The subdivision tolerance to use.			     M
*   NumericTol:	  The numerical tolerance to use.			     M
*									     *
* RETURN VALUE:                                                              M
*   MvarPtStruct *: List of points on the solution set.  Dimension of the    M
*		    points will be the same as the dimensions of all MVETs.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros0D				                             M
*									     *
* KEYWORDS:                                                                  M
*   MvarETsZeros0D			                                     M
*****************************************************************************/
MvarPtStruct *MvarETsZeros0D(MvarExprTreeStruct * const *MVETs,
			     MvarConstraintType *Constraints,
			     int NumOfMVETs,
			     CagdRType SubdivTol,
			     CagdRType NumericTol)
{
    CagdBType FreeConstraints;
    MvarZeroPrblmStruct *Problem;

    /* If we need to allocate the constraints (all equality) - do it now. */
    if (Constraints == NULL) {
	int i;
        Constraints = (MvarConstraintType *)
	                 IritMalloc(sizeof(MvarConstraintType) * NumOfMVETs);
        for (i = 0; i < NumOfMVETs; i++)
	    Constraints[i] = MVAR_CNSTRNT_ZERO;

        FreeConstraints = TRUE;
    }
    else
        FreeConstraints = FALSE;

    if ((Problem = MvarZeroSolverPrblmNew(NULL,
					  (const MvarExprTreeStruct * const *)
									 MVETs,
					  NumOfMVETs, Constraints,
					  SubdivTol, NumericTol,
					  IRIT_INFNTY, FALSE)) != NULL) {
        MvarZeroPrblmIntrnlStruct
	    *PrblmIntrnl = Problem -> _Internal;
	MvarPtStruct
	    *ZeroSet = NULL;
	MvarZeroSolutionStruct *Solution;

	switch (PrblmIntrnl -> ExpectedSolutionDim) {
	    case 0:
		if ((Solution = MvarZeroSolver(Problem)) != NULL) {
		    assert(Solution -> ActiveRepresentation == 
					      MVAR_ZER_SLVR_SOLUTION_PT_LIST);
		    ZeroSet = Solution -> U.Pt;
		}
		MvarZeroSolverSolutionFree(Solution, FALSE);
		break;
	    default:
	        assert(0);
		break;
	}

	MvarZeroSolverPrblmFree(Problem);

	if (FreeConstraints)
	    IritFree(Constraints);

	return ZeroSet;
    }

    if (FreeConstraints)
        IritFree(Constraints);

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Parses the constraints vector and verify the types of constrains:        *
* ZERO constraints first, followed by ZERO_SUBDIV constraints, followed by   *
* non zero (POSITIVE or NEGATIVE) constraints.                               *
*                                                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   Constraints:        The vector of constraints.                           *
*   NumOfMvs:           Size of Constraints vector.                          *
*   NumOfZeroMVs:       Will be updated with the number of ZERO constraints. *
*   NumOfZeroSubdivMVs: Will be updated with number of ZERO and ZERO_SUBDIV  *
*                       cnstrnts.					     *
*   NumOfNonZeroMVs:    Will be updated with number of non ZERO cnstrnts.    *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:                                                                     *
*****************************************************************************/
static int MvarZeroVerifyConstraintTypes(const MvarConstraintType *Constraints,
					 int NumOfMVs,
					 int *NumOfZeroMVs,
					 int *NumOfZeroSubdivMVs,
					 int *NumOfNonZeroMVs)
{
    int i;

    *NumOfZeroMVs = *NumOfZeroSubdivMVs = *NumOfNonZeroMVs = 0;

    /* For the unlikely case of all constraints identically zero. */
    if (NumOfMVs == 0)
	return TRUE;

    /* Make sure all zero constraints are first, zero-subdiv are second,    */
    /* and count how many.  Then, non zero (i.e. positive) constraints.	    */
    for (i = 0; i < NumOfMVs; i++) {
        if (Constraints[i] != MVAR_CNSTRNT_ZERO)
	    break;
    }
    *NumOfZeroMVs = i;

    for ( ; i < NumOfMVs; i++) {
        if (Constraints[i] == MVAR_CNSTRNT_ZERO)
	    return FALSE;

        if (Constraints[i] != MVAR_CNSTRNT_ZERO_SUBDIV)
	    break;
    }
    *NumOfZeroSubdivMVs = i;
    *NumOfNonZeroMVs = NumOfMVs - i;

    for ( ; i < NumOfMVs; i++) {
	if (Constraints[i] == MVAR_CNSTRNT_ZERO ||
	    Constraints[i] == MVAR_CNSTRNT_ZERO_SUBDIV)
	    return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Construction of a new zero finding problem structure. Allocates the memory M
* and assigns slots that are already known upon construction.		     M
*   This routine should NOT be used when extracting a sub-problem from an    M
* existing problem (such as after domain reduction or subdivision): it	     M
* should only be called at very first construction (SubdivDepth == 0).       M
*   Scaling and organization takes place here that is not required for       M
* sub-problems.								     M
*									     *
* PARAMETERS:                                                                M
*   MVs:         Array of multivariates, NULL for ETs representations.       M
*   ETs:         Array of expression trees, NULL for MVs representations.    M
*   NumOfConstraints: Number of constraints in the problem (Total).          M
*   Constraints: Either an equality or an inequality type of constraint.     M
*   SubdivTol:   Tolerance of the subdivision process.  Tolerance is         M
*		 measured in the parametric space of the constraints.        M
*   NumericTol:  Tolerance of the numeric stage.			     M
*   StepTol:	 In 1D numeric stages (curve tracing)- the step size used.   M 
*   Solve2DBy0D: A temporary flag to indicate the use of the new 2D solver   M
*		 or not.					             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroPrblmStruct *:    The new problem structure or NULL if error.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverPrblmFree				                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverPrblmNew                                                   M
*****************************************************************************/
MvarZeroPrblmStruct *MvarZeroSolverPrblmNew(const MvarMVStruct * const *MVs, 
				            const MvarExprTreeStruct * const
									  *ETs,
				            int NumOfConstraints,
					    MvarConstraintType *Constraints,
					    CagdRType SubdivTol,
					    CagdRType NumericTol,
					    CagdRType StepTol,
					    CagdBType Solve2DBy0D)
{
    int i, l, NumOfZeroConstraints, NumOfZeroSubdivConstraints, NumOfNonZeroMVs;
    MvarMVStruct **LclMVs;
    MvarZeroPrblmStruct 
	*Problem = (MvarZeroPrblmStruct *)
		                     IritMalloc(sizeof(MvarZeroPrblmStruct));
    MvarZeroPrblmIntrnlStruct
        *PrblmIntrnl = (MvarZeroPrblmIntrnlStruct *)
			        IritMalloc(sizeof(MvarZeroPrblmIntrnlStruct));

    if (!MvarZeroVerifyConstraintTypes(Constraints, NumOfConstraints,
				       &NumOfZeroConstraints,
				       &NumOfZeroSubdivConstraints,
				       &NumOfNonZeroMVs)) {
        MVAR_FATAL_ERROR(MVAR_ERR_ZER_ORDER_CNSTRCT);
        return NULL;
    }

    /* Make sure we are reasonable. */
    SubdivTol = IRIT_MAX(SubdivTol, 1e-8);

    IRIT_ZAP_MEM(Problem, sizeof(MvarZeroPrblmStruct));
    IRIT_ZAP_MEM(PrblmIntrnl, sizeof(MvarZeroPrblmIntrnlStruct));
    Problem -> _Internal = PrblmIntrnl;

    Problem -> Constraints = Constraints;
    Problem -> NumOfConstraints = NumOfConstraints;
    Problem -> NumOfZeroConstraints = NumOfZeroConstraints;
    Problem -> StepTol = StepTol;
    Problem -> SubdivTol = SubdivTol;
    Problem -> NumericTol = NumericTol;
    PrblmIntrnl -> ParamOfSubdiv = IRIT_INFNTY;
    PrblmIntrnl -> SubdivDir = -1;
    PrblmIntrnl -> NumOfZeroSubdivConstraints = NumOfZeroSubdivConstraints;
    PrblmIntrnl -> DBGCntGuaranteedDmns = (int *) IritMalloc(sizeof(int));
    *PrblmIntrnl -> DBGCntGuaranteedDmns = 0;
    PrblmIntrnl -> DBGCntSinglrDmns = (int *) IritMalloc(sizeof(int));
    *PrblmIntrnl -> DBGCntSinglrDmns = 0;
    

    if (MVs != NULL && ETs == NULL) {
	Problem -> ActiveRepresentation = MVAR_ZER_SLVR_MVS;

	/* Make sure parametric domain is identical at all multivariates. */
	for (i = 1; i < NumOfConstraints; i++) {
	    if (MVs[0] -> Dim != MVs[i] -> Dim) {
	        MVAR_FATAL_ERROR(MVAR_ERR_INCONS_DOMAIN);
		IritFree(Problem);
		IritFree(PrblmIntrnl);
		return NULL;
	    }
	}
	for (l = 0; l < MVs[0] -> Dim; l++) {
	    CagdRType Min, Max;

	    MvarMVDomain(MVs[0], &Min, &Max, l);

	    for (i = 1; i < NumOfConstraints; i++) {
	        CagdRType Min2, Max2;

		MvarMVDomain(MVs[i], &Min2, &Max2, l);
		if (!IRIT_APX_EQ(Min, Min2) || !IRIT_APX_EQ(Max, Max2)) {
		    MVAR_FATAL_ERROR(MVAR_ERR_INCONS_DOMAIN);
		    IritFree(Problem);
		    IritFree(PrblmIntrnl);
		    return NULL;
		}
	    }
	}

	PrblmIntrnl -> ExpectedSolutionDim = 
	         IRIT_MAX(0, MVs[0] -> Dim - Problem -> NumOfZeroConstraints);
	if (Solve2DBy0D && PrblmIntrnl -> ExpectedSolutionDim == 2)
	    PrblmIntrnl -> ExpectedSolutionDim = 0;

	switch (PrblmIntrnl -> ExpectedSolutionDim) {
	    default:
	    case 0:
	        LclMVs = MvarZeroOrganizeMVs0DProblem(MVs, Constraints, 
						      &NumOfConstraints);
		break;
	    case 1:
	        LclMVs = MvarZeroOrganizeMVs1DProblem(MVs, Constraints,
						      &NumOfConstraints);
		break;
	    case 2:
		LclMVs = (MvarMVStruct **)
		    IritMalloc(sizeof(MvarMVStruct *) * NumOfConstraints);
		for (i = 0; i < NumOfConstraints; i++)
		    LclMVs[i] = MvarMVCopy(MVs[i]);
	}

        /* Scale MVs */ 
	for (i = 0; i < NumOfConstraints; i++) {
	    int Length = MVAR_CTL_MESH_LENGTH(LclMVs[i]);
	    CagdRType Scale, *R,
		*Pts = LclMVs[i] -> Points[1],
		Max = 0.0;

	    for (R = &Pts[Length]; R-- != Pts; )
		Max = IRIT_MAX(Max, IRIT_FABS(*R));

	    Scale = MVAR_ZERO_PRECOND_SCALE / Max;

	    for (R = &Pts[Length]; R-- != Pts; )
		*R *= Scale;
	}


	/* Make sure all MVs are B-spline type, and update in place if not: */
	MvarZeroVerifyBspConstraints(LclMVs, NumOfConstraints);

	/* Make sure we have a valid problem with non completely zero       */
	/* constraints. Can happen if all constraints are identically zero. */
	if (LclMVs == NULL) {
#	    ifdef DEBUG
		fprintf(stderr, "Warning: All constraints were purged, being singular (completely zero!?).\n");
#	    endif /* DEBUG */
	    IritFree(Problem);
	    IritFree(PrblmIntrnl);
	    return NULL;
	}
	else if (Problem -> NumOfConstraints != NumOfConstraints) {
#	    ifdef DEBUG
	        fprintf(stderr, "Warning: Some constraints were purged, being singular (completely zero!?).\n");
#	    endif /* DEBUG */
	    IritFree(Problem);
	    IritFree(PrblmIntrnl);
	    return NULL;
	}

	Problem -> U.MVs = LclMVs;
	PrblmIntrnl -> NumOfZeroSubdivConstraints = NumOfZeroSubdivConstraints;
	PrblmIntrnl -> HasC1Discont = MvarZeroHasC1Discont(
			    Problem -> U.MVs, Problem -> NumOfZeroConstraints, 
			    &(PrblmIntrnl -> C1DiscontInd), 
			    &(PrblmIntrnl -> C1DiscontParam));

	/* Since we are at depth zero, if the problem is C1 smooth, it is    */
	/* obviously the first time it is C1 smooth. This is important since */
	/* some actions are executed only on smooth problems, and only upon  */
	/* the first time along the subdivision tree the smoothness appears. */
	PrblmIntrnl -> IsFirstSmooth = !PrblmIntrnl -> HasC1Discont;
	PrblmIntrnl -> ParamPerturb = SubdivTol * MVAR_ZERO_PARAM_REL_PERTURB;
	PrblmIntrnl -> ZeroMVsSameSpace = MvarMVsZerosSameSpace(
							Problem -> U.MVs,
							NumOfZeroConstraints);

	/* If all domains are Bezier-like, and we use domain reduction with  */
	/* gradient preconditioning, make sure all functions share space.    */
	if (NumOfZeroConstraints == Problem -> U.MVs[0] -> Dim &&
	    _MVGlblZeroApplyDomainReduction &&
	    _MVGlblZeroApplyGradPreconditioning &&
	    !PrblmIntrnl -> ZeroMVsSameSpace) {/* Promote MVs to same space. */
	    /* First pass will leave MVs[0] as most common denominator. */
	    for (i = 1; i < NumOfZeroConstraints; i++)
	        MvarMakeMVsCompatible(&LclMVs[0], &LclMVs[i], TRUE, TRUE);

	    /* Second pass will elevate all MVs[i] to MVs[0] space. */
	    for (i = 1; i < NumOfZeroConstraints; i++)
	        MvarMakeMVsCompatible(&LclMVs[0], &LclMVs[i], TRUE, TRUE);

	    PrblmIntrnl -> ZeroMVsSameSpace = TRUE;
	}
    }
    else if ((MVs == NULL) && (ETs != NULL)) {
	Problem -> ActiveRepresentation = MVAR_ZER_SLVR_EXP_TREE;

	/* Make sure domain is identical in all MV expression trees. */
	for (i = 1; i < NumOfConstraints; i++) {
	    if (ETs[0] -> Dim != ETs[i] -> Dim) {
	        MVAR_FATAL_ERROR(MVAR_ERR_INCONS_DOMAIN);
		IritFree(Problem);
		IritFree(PrblmIntrnl);
		return NULL;
	    }
	}
	for (l = 0; l < ETs[0] -> Dim; l++) {
	    CagdRType Min, Max;

	    /* Find the domain in this l'th axis.  Note an ET can have no */
	    /* domain defined in some (dont care) axis.			  */
	    for (i = 0; i < NumOfConstraints; i++)
	        if (MvarETDomain(ETs[i], &Min, &Max, l))
		    break;

	    for (i = 1; i < NumOfConstraints; i++) {
	        CagdRType Min2, Max2;

		if (MvarETDomain(ETs[i], &Min2, &Max2, l)) {
		    if (!IRIT_APX_EQ(Min, Min2) || !IRIT_APX_EQ(Max, Max2)) {
		        MVAR_FATAL_ERROR(MVAR_ERR_INCONS_DOMAIN);
			IritFree(Problem);
			IritFree(PrblmIntrnl);
			return NULL;
		    }
		}
	    }
	}

	PrblmIntrnl -> ExpectedSolutionDim = 
	         IRIT_MAX(0, ETs[0] -> Dim - Problem -> NumOfZeroConstraints);

	switch (PrblmIntrnl -> ExpectedSolutionDim) {
	    default:
	    case 0:
	        Problem -> U.Eqns = MvarZeroOrganizeETs0DProblem(ETs,
							    NumOfConstraints);
		Problem -> U.Eqns -> NumEqns = NumOfConstraints;
		Problem -> U.Eqns -> NumZeroEqns = NumOfZeroConstraints;
		Problem -> U.Eqns -> NumZeroSubdivEqns =
						    NumOfZeroSubdivConstraints;
		IRIT_GEN_COPY(Problem -> U.Eqns -> ConstraintTypes,
			      Problem -> Constraints, 
			      sizeof(MvarConstraintType) * NumOfConstraints);
		break;

	}

	/* Can happen if all constraints are identically zero. */
	if (Problem -> U.Eqns == NULL) {
	    IritFree(Problem);
	    IritFree(PrblmIntrnl);
	    return NULL;
	}

	PrblmIntrnl -> NumOfZeroSubdivConstraints = NumOfZeroSubdivConstraints;
	PrblmIntrnl -> ParamPerturb = SubdivTol * MVAR_ZERO_PARAM_REL_PERTURB;
    }
    else {
	IritFree(Problem);
	IritFree(PrblmIntrnl);
	MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
	return NULL;
    }

    /* Some dimension specific initializations: */
    switch (PrblmIntrnl -> ExpectedSolutionDim) {
	case 0:
	    switch (Problem -> ActiveRepresentation) {
	        case MVAR_ZER_SLVR_MVS:
	            MvarZeroSolverSetCallbackFcns0DMVs(Problem);
		    break;
  	        case MVAR_ZER_SLVR_EXP_TREE:
		    MvarZeroSolverSetCallbackFcns0DExpTr(Problem);
		    break;
	        default:
		    assert(0);
		    MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
		    return NULL;
	    }
	    break;
	case 1:
	     switch (Problem -> ActiveRepresentation) {
	        case MVAR_ZER_SLVR_MVS:
		    MvarZeroSolverSetCallbackFcns1DMVs(Problem);
		    Problem -> AS = MvarMVZR1DAllocOnce(
				             Problem -> NumOfZeroConstraints);
		    break;
  	        case MVAR_ZER_SLVR_EXP_TREE:
		    MvarZeroSolverSetCallbackFcns1DExpTr(Problem);
		    break;
	        default:
		    assert(0);
		    MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
		    return NULL;
	    }
	    break;
	case 2:
	     switch (Problem -> ActiveRepresentation) {
	        case MVAR_ZER_SLVR_MVS:
		    MvarZeroSolverSetCallbackFcns2DMVs(Problem);
		    break;
	        case MVAR_ZER_SLVR_EXP_TREE:
		    MvarZeroSolverSetCallbackFcns2DExpTr(Problem);
		    break;
	        default:
		    assert(0);
		    MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
		    return NULL;
    	    }
	    break;
	default:
	    assert(0);
	    MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
	    return NULL;
    }

    MVAR_ZERO_SLVR_APPLY(UpdateDmn)(Problem);

    return Problem;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construction of a new zero finding sub-problem structure from an         M
* existing problem, for the given constraints defined on the sub-domain      M
* (the result of subdivision or domain reduction).			     M
*   Note that there are slots that are only allocated once, while some are   M
* are allocated at all depths.						     M
*									     *
* PARAMETERS:                                                                M
*   Problem:	The original problem, from which to extract the sub-problem. M
*   MVs:	Array of multivariates, NULL if not MVs representation.      M
*   Eqns:	Struct of expression trees, NULL if not ETs representation.  M
*   BoundarySol: The solution to the corresponding 2 * Dim problems of one   M
*		dimension less, defined on the boundary of the new, smaller  M
*		problem.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroPrblmStruct *:    The sub-problem.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverPrblmNew, MvarZeroOrganizeMVs0DProblem,	             M
*   MvarZeroOrganizeMVs1DProblem					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverSubProblem                                                 M
*****************************************************************************/
MvarZeroPrblmStruct *MvarZeroSolverSubProblem(
					   MvarZeroPrblmStruct const *Problem,
					   MvarMVStruct **MVs,
					   MvarExprTreeEqnsStruct *Eqns,
					   MvarZeroSolutionStruct *BoundarySol)
{
    MvarZeroPrblmStruct 
	*NewProblem = (MvarZeroPrblmStruct *) 
		                      IritMalloc(sizeof(MvarZeroPrblmStruct));
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal,
        *NewPrblmIntrnl = (MvarZeroPrblmIntrnlStruct *)
			        IritMalloc(sizeof(MvarZeroPrblmIntrnlStruct));

    IRIT_ZAP_MEM(NewProblem, sizeof(MvarZeroPrblmStruct));
    IRIT_ZAP_MEM(NewPrblmIntrnl, sizeof(MvarZeroPrblmIntrnlStruct));
    NewProblem -> _Internal = NewPrblmIntrnl;

    NewProblem -> ActiveRepresentation = Problem -> ActiveRepresentation;
    NewProblem -> AS = Problem -> AS; 
    NewProblem -> Constraints = Problem -> Constraints;
    NewProblem -> NumOfConstraints = Problem -> NumOfConstraints;
    NewProblem -> NumOfZeroConstraints = Problem -> NumOfZeroConstraints;
    NewProblem -> StepTol = Problem -> StepTol;
    NewProblem -> SubdivTol = Problem -> SubdivTol;
    NewProblem -> NumericTol = Problem -> NumericTol;

    NewPrblmIntrnl -> CallbackFunctions = PrblmIntrnl -> CallbackFunctions;
    NewPrblmIntrnl -> NumOfZeroSubdivConstraints = 
				    PrblmIntrnl -> NumOfZeroSubdivConstraints;
    NewPrblmIntrnl -> SubdivDepth = PrblmIntrnl -> SubdivDepth + 1;
    NewPrblmIntrnl -> DBGCntGuaranteedDmns = 
	PrblmIntrnl -> DBGCntGuaranteedDmns;
    NewPrblmIntrnl -> DBGCntSinglrDmns = PrblmIntrnl -> DBGCntSinglrDmns;


    /* Properties that are inherited to sub-domains: */
    NewPrblmIntrnl -> NoLoopTest1D = PrblmIntrnl -> NoLoopTest1D;
    NewPrblmIntrnl -> SingleSolutionTest0D = PrblmIntrnl -> SingleSolutionTest0D;
    NewPrblmIntrnl -> ExpectedSolutionDim = PrblmIntrnl -> ExpectedSolutionDim;
    NewPrblmIntrnl -> ConstructionComplete = PrblmIntrnl -> ConstructionComplete;
    NewPrblmIntrnl -> ParamPerturb = PrblmIntrnl -> ParamPerturb;
    NewPrblmIntrnl -> ZeroMVsSameSpace = PrblmIntrnl -> ZeroMVsSameSpace;
    NewPrblmIntrnl -> OrigMVMaxDmn = PrblmIntrnl -> OrigMVMaxDmn;
    NewPrblmIntrnl -> OrigMVMinDmn = PrblmIntrnl -> OrigMVMinDmn;

    NewPrblmIntrnl -> ParamOfSubdiv = IRIT_INFNTY;
    NewPrblmIntrnl -> SubdivDir = -1;

    switch (Problem -> ActiveRepresentation) {
        case MVAR_ZER_SLVR_MVS:
            NewProblem -> U.MVs = MVs;
	    MVAR_ZERO_SLVR_APPLY(UpdateDmn)(NewProblem);

	    if (PrblmIntrnl -> HasC1Discont) {       /* Still not C1 smooth. */
	        if (MvarZeroHasC1Discont(NewProblem -> U.MVs, 
					 NewProblem -> NumOfConstraints,
					 &(NewPrblmIntrnl -> C1DiscontInd),
					 &(NewPrblmIntrnl -> C1DiscontParam))) {
		    NewPrblmIntrnl -> HasC1Discont = TRUE;
		}
		else { /* C1 smooth, and for the first time. */
		    NewPrblmIntrnl -> IsFirstSmooth = TRUE;
		}
	    }

	    NewPrblmIntrnl -> U.FirstSmoothMVs =
	                                      PrblmIntrnl -> U.FirstSmoothMVs;
	    NewPrblmIntrnl -> MVGradients = PrblmIntrnl -> MVGradients;
	    NewPrblmIntrnl -> SolutionsOnBoundary = BoundarySol;
	    NewPrblmIntrnl -> NumOfBoundarySolutions =
	    NewPrblmIntrnl -> SolutionsOnBoundary ?
	        CagdListLength(NewPrblmIntrnl -> SolutionsOnBoundary -> U.Pt) :
	        0;

            break;
        case MVAR_ZER_SLVR_EXP_TREE:
	    NewProblem -> U.Eqns = Eqns;
	    MVAR_ZERO_SLVR_APPLY(UpdateDmn)(NewProblem);
	    break;
        default:
            assert(0);
    }

    return NewProblem;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deallocates and frees all slots of a problem structure of a zero         M
* finding problem. Constraints are never freed, since they are allocated     M
* externally. Note that some slots are freed only at depth zero, some are    M
* freed only at first C1 smooth instance.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:	    Problem structure to free.                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverPrblmNew, MvarZeroSolverSubProblem                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverPrblmFree                                                  M
*****************************************************************************/
void MvarZeroSolverPrblmFree(MvarZeroPrblmStruct *Problem)
{
    if (Problem != NULL) {
	MvarZeroPrblmIntrnlStruct
	    *PrblmIntrnl = Problem -> _Internal;
	int i,
	    n = Problem -> NumOfConstraints,
	    n1 = Problem -> NumOfZeroConstraints,
	    Depth = PrblmIntrnl -> SubdivDepth;

	if (Problem -> Attr != NULL)
	    _AttrFreeAttributeData(Problem -> Attr);

	switch (Problem -> ActiveRepresentation) {
	    case MVAR_ZER_SLVR_MVS:
	        for (i = 0; i < n; i++) {
		    MvarMVFree(Problem -> U.MVs[i]);
		}
		IritFree(Problem -> U.MVs);
		if (PrblmIntrnl -> IsFirstSmooth) {
		    if (PrblmIntrnl -> MVGradients != NULL) {
		        for (i = 0; i < n1; i++) {
			    MvarMVFreeGradient(PrblmIntrnl -> MVGradients[i]);
			}
			IritFree(PrblmIntrnl -> MVGradients);
		    }
		    if (PrblmIntrnl -> U.FirstSmoothMVs != NULL) {
		        for (i = 0; i < n; i++) 
			    MvarMVFree(PrblmIntrnl -> U.FirstSmoothMVs[i]);
			IritFree(PrblmIntrnl -> U.FirstSmoothMVs);
		    }
		}
		break;
	    case MVAR_ZER_SLVR_EXP_TREE:
		MvarExprTreeEqnsFree(Problem -> U.Eqns);
		break;
	    default:
		assert(0);
	}

	if (PrblmIntrnl -> SolutionsOnBoundary != NULL) 
	    MvarZeroSolverSolutionFree(PrblmIntrnl -> SolutionsOnBoundary,
				       TRUE);

	if (Problem -> AS != NULL && Depth == 0)
	    MvarMVZR1DDeallocOnce(Problem -> AS);

	if (PrblmIntrnl -> MVMinDmn)
	    IritFree(PrblmIntrnl -> MVMinDmn);
	if (PrblmIntrnl -> MVMaxDmn)
	    IritFree(PrblmIntrnl -> MVMaxDmn);
	if (PrblmIntrnl -> TJList != NULL)
	    MvarZeroTJFreeList(PrblmIntrnl -> TJList);
	if (Depth == 0) {
	    if (PrblmIntrnl -> OrigMVMinDmn)
		IritFree(PrblmIntrnl -> OrigMVMinDmn);
	    if (PrblmIntrnl -> OrigMVMaxDmn)
		IritFree(PrblmIntrnl -> OrigMVMaxDmn);
	    if (PrblmIntrnl -> DBGCntGuaranteedDmns != NULL)
		IritFree(PrblmIntrnl -> DBGCntGuaranteedDmns);
	    if (PrblmIntrnl -> DBGCntSinglrDmns != NULL)
		IritFree(PrblmIntrnl -> DBGCntSinglrDmns);
	}
	IritFree(PrblmIntrnl);
	IritFree(Problem);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new solution structure of a zero     M
* finding problem.							     M
*									     *
* PARAMETERS:                                                                M
*   Tr:	     A list of triangles, or NULL if the representation is not such. M
*   Pl:      A list of polylines, or NULL if the representation is not such. M
*   Pt:      A list of points, or NULL if the representation is not such.    M
*   Rep:     The active representation required for the new solution.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroSolutionStruct *:    The new solution or NULL if error.          M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverSolutionFree                                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverSolutionNew                                                M
*****************************************************************************/
MvarZeroSolutionStruct *MvarZeroSolverSolutionNew(MvarTriangleStruct *Tr,
						  MvarPolylineStruct *Pl,
						  MvarPtStruct *Pt,
						  Representation Rep)
{
    MvarZeroSolutionStruct 
	*Solution = NULL;

    if (!Pl && !Pt && !Tr)
	return NULL;
    Solution = (MvarZeroSolutionStruct *) 
				  IritMalloc(sizeof(MvarZeroSolutionStruct));
    Solution -> ActiveRepresentation = Rep;
    Solution -> Attr = NULL;
    Solution -> Pnext = NULL;
    Solution -> TJList = NULL;

    switch (Rep) {
	case MVAR_ZER_SLVR_SOLUTION_TR_LIST:
	    assert(Tr != NULL);
	    Solution -> U.Tr = Tr;
	    break;
        case MVAR_ZER_SLVR_SOLUTION_POLYLINE:
	    assert(Pl != NULL);
	    Solution -> U.Pl = Pl;
	    break;
        case MVAR_ZER_SLVR_SOLUTION_PT_LIST:
	    assert(Pt != NULL);
	    Solution -> U.Pt = Pt;
	    break;
        default:
	    MVAR_FATAL_ERROR(MVAR_ERR_ZER_SOL_CNSTRCT);
	    assert(0);
	    break;
    }

    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Deallocates and frees all slots of a solution structure of a zero          M
* finding problem. NOTE: the T-Junction list is not freed, as it is always   M
* freed as part of the problem deallocation.                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Solution:   Solution structure to free.                                  M
*   FreeUnion:	If TRUE, the points/polylines stored at in the solution are  M
*		freed, otherwise they are not. This is useful for the cases  M
*		of uniting solutions in place.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolverSolutionNew                                                M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverSolutionFree                                               M
*****************************************************************************/
void MvarZeroSolverSolutionFree(MvarZeroSolutionStruct *Solution,
				CagdBType FreeUnion)
{
    if (Solution != NULL) {
	if (Solution -> Attr != NULL)
	    _AttrFreeAttributeData(Solution -> Attr);

	if (FreeUnion) {
	    switch (Solution -> ActiveRepresentation) {
	        case MVAR_ZER_SLVR_SOLUTION_PT_LIST:
		    MvarPtFreeList(Solution -> U.Pt);
		    break;
	        case MVAR_ZER_SLVR_SOLUTION_POLYLINE:
		    MvarPolylineFreeList(Solution -> U.Pl);
		    break;
		case MVAR_ZER_SLVR_SOLUTION_TR_LIST:
		    MvarTriangleFreeList(Solution -> U.Tr);
		    break;
	        default:
		    assert(0);
	    }
	}
	if (Solution -> TJList != NULL)
	    MvarZeroTJFreeList(Solution -> TJList);
	IritFree(Solution);
    }
}

/*****************************************************************************
* DESCRIPTION:								     M
*   The general zero finding problem solver. Invokes various zero finding    M
* algorithms according to the function representations/solution set dim.     M
*   The solver follows as generally as possible the following paradigm:	     M
* 0. Recursively subdivide until the problem is C1 smooth.                   M
* 1. Try to rule out the possibility of a solution.			     M
* 2. If can't roule out, Check if the topology is guaranteed.		     M
* 3. If topology is guaranteed, improve/reconstruct numerically.	     M
* 4. If not, or if numeric step failed subdivide and solve recursively.      M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroSolutionStruct *: The point set (zero dimensional case) or       M
*			    piecewise linear approximation to the solution   M
*			    manifold (one/two dimensional  case).	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolver                                                           M
*****************************************************************************/
MvarZeroSolutionStruct *MvarZeroSolver(MvarZeroPrblmStruct *Problem)
{
    CagdBType 
	IsTopoGuaranteed = FALSE,
	CanPurgeDmn = FALSE;
    int i;
    MvarPtStruct *PtList;
    MvarZeroSolutionStruct
	*Sol1 = NULL,
	*Sol2 = NULL,
	*Solution = NULL;
    MvarZeroPrblmStruct **SubProblems;
    MvarZeroPrblmIntrnlStruct *PrblmIntrnl;

#ifdef DEBUG
    MvarZeroDBGDmn(Problem);
#endif /* DEBUG */

    if (Problem == NULL)
        return NULL;

    PrblmIntrnl = Problem -> _Internal;

    if (!PrblmIntrnl -> UnderSubdivTol) {
        /* Step 0: If C1 discontinuous, subdivide at the discont' and solve  */
        /* the sub-problems.						     */
        if (PrblmIntrnl -> HasC1Discont) {
	    SubProblems = MvarZeroC1DiscontSubdiv(Problem);
	    Sol1 = MvarZeroSolver(SubProblems[0]);
	    Sol2 = MvarZeroSolver(SubProblems[1]);
	    MvarZeroSolverPrblmFree(SubProblems[0]);
	    MvarZeroSolverPrblmFree(SubProblems[1]);
	    IritFree(SubProblems);
	    Solution = MVAR_ZERO_SLVR_APPLY(UniteSolutions)(
						 Sol1, Sol2,
				  		 PrblmIntrnl -> C1DiscontInd,
						 PrblmIntrnl -> C1DiscontParam,
						 Problem);
	    return PrblmIntrnl -> SubdivDepth == 0 ? 
	           MVAR_ZERO_SLVR_APPLY(OrganizeSolution)(Solution, Problem) :
	           Solution;
	}

	/* Step 1: Can we rule out the possibility of a solution?            */
	/* Note:							     */
	/* In fully constrained problems: This is checked in a full manner   */
	/* only at depth 0, since deeper down the subdiv tree we incorporate */
	/* it in the subdivision process (the "Early sign evaluation" idea). */
	if (PrblmIntrnl -> SubdivDepth == 0 || 
	    PrblmIntrnl -> ExpectedSolutionDim > 0 ||
	    Problem -> ActiveRepresentation != MVAR_ZER_SLVR_MVS) {
	    for (i = 0; i < Problem -> NumOfConstraints; i++) {
	        if ((CanPurgeDmn =
		    MVAR_ZERO_SLVR_APPLY(NoZeroCrossTest)(Problem, i))
		                                                   != FALSE) {
		    break;
		}
	    }
	}

        if (CanPurgeDmn) {
#	    ifdef DEBUG
	    if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
	        printf("\nDomain %d is purged by coeff signs of constraint %d ***********************\n",
		       _DebugZeroSolverDmnNum, i);
	    }
#	    endif /* DEBUG */
	    return NULL;
	}

	/* If cannot purge: Check if the problem is C1 smooth for the first  */
	/* time. If so, this is the right time for the effort of adding the  */
	/* gradients, adding the first smooth original constraints, and      */
	/* to solve the boundary.					     */
	if (PrblmIntrnl -> HasC1Discont == FALSE && 
	    PrblmIntrnl -> IsFirstSmooth == TRUE) {
	    if (!MVAR_ZERO_SLVR_APPLY(FirstSmoothUpdates)(Problem)) {
	        return NULL;
	    }
	}
	/* If we decide to use the Kantorovich test for 0D solutions,  we    */
	/* do not follow the general framework - we turn to it here:	     */
	if (_MVGlblZeroApplyKantorovichTest && 
				     PrblmIntrnl -> ExpectedSolutionDim == 0) {
	    PtList = MvarZeroGetRootsByKantorovich(
					       Problem -> U.MVs,
					       Problem -> Constraints,
					       Problem -> NumOfConstraints,
					       Problem -> NumOfZeroConstraints,
					       _MVGlblZeroApplyNormalConeTest,
					       Problem -> SubdivTol, 
					       PrblmIntrnl -> SubdivDepth,
					       PrblmIntrnl -> ZeroMVsSameSpace,
					       PrblmIntrnl -> ParamPerturb);
	    if (PtList != NULL) {
		Solution = MvarZeroSolverSolutionNew(
			   NULL, NULL, PtList, MVAR_ZER_SLVR_SOLUTION_PT_LIST);
		return PrblmIntrnl -> SubdivDepth == 0 ? 
		    MVAR_ZERO_SLVR_APPLY(OrganizeSolution)(Solution, Problem) :
		    Solution;
		
	    }
	}

	/* Do we have a call back test function?  If so call it. */
	if (_MVGlblZeroSubdivCallBackFunc != NULL &&
	    _MVGlblZeroSubdivCallBackFunc(Problem, PrblmIntrnl -> SubdivDepth)) {
	    return NULL;
	}

	/* Step 2: Is the topology guaranteed? */
	IsTopoGuaranteed = MVAR_ZERO_SLVR_APPLY(GuaranteedTopologyTest)(
								     Problem);

	/* Step 3: If topo' is guaranteed, improve/reconstruct numerically. */
	if (IsTopoGuaranteed) {
	    Solution = MVAR_ZERO_SLVR_APPLY(NumericImprovement)(Problem);
	    if (Solution != NULL || PrblmIntrnl -> PurgeByNumericStep) {
		return PrblmIntrnl -> SubdivDepth == 0 ? 
		    MVAR_ZERO_SLVR_APPLY(OrganizeSolution)(Solution, Problem) :
		    Solution;
	    }
	}

	/* Step 4: Subdivision. If topo' is not guaranteed or if the numeric */
	/* step failed, subdivide the domain and solve recursively.          */
	if ((SubProblems = MVAR_ZERO_SLVR_APPLY(SubdivProblem)(Problem))
	                                                            == NULL) {
	    switch (Problem -> ActiveRepresentation) {
	        case MVAR_ZER_SLVR_EXP_TREE:
		    /* This flag must be on for SubProblems to be NULL. In  */
		    /* such a case, we are to convert the ETs that holds    */
		    /* Bezier MVs only to Bezier MVs and solve MVs problem. */    
		    assert(_MVGlblZeroETUseCommonExpr);
		    return MvarZeroETs2MVs(Problem -> U.Eqns,
					   Problem -> SubdivTol,
					   Problem -> NumericTol);
	        default:
		    assert(0); /* Should not happen. */
		    return NULL;
	    }
	}

	Sol1 = MvarZeroSolver(SubProblems[0]);
	Sol2 = MvarZeroSolver(SubProblems[1]);
	MvarZeroSolverPrblmFree(SubProblems[0]);
	MvarZeroSolverPrblmFree(SubProblems[1]);
	IritFree(SubProblems);
	Solution = MVAR_ZERO_SLVR_APPLY(UniteSolutions)(
						 Sol1, Sol2,
				  		 PrblmIntrnl -> SubdivDir,
						 PrblmIntrnl -> ParamOfSubdiv,
						 Problem);
    }
    else { /* Subdivision tolerance is reached. */
        Solution = MVAR_ZERO_SLVR_APPLY(SingularSolution)(Problem);
    }

    if (PrblmIntrnl -> SubdivDepth == 0) 
        return MVAR_ZERO_SLVR_APPLY(OrganizeSolution)(Solution, Problem);
    else 
	return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Converts the given expression trees equations to multivariates and solve *
* using MVs zero solver.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   Eqns:         Equations/constraints with only MVs.			     *
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        *
*		  measured in the parametric space of the multivariates.     *
*   NumericTol:   Numeric tolerance of the numeric stage.  The numeric stage *
*		  is employed only if IRIT_FABS(NumericTol) < SubdivTol.     *
*		  If NumericTol is negative, points that fail to improve     *
*		  numerically are purged away.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:   The solution set - list of points.	     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroETs2MVs(const MvarExprTreeEqnsStruct
					                                *Eqns,
					       CagdRType SubdivTol,
					       CagdRType NumericTol)
{
    int i,
        NumOfMVs = Eqns -> NumEqns;
    MvarPtStruct *MVPts;
    MvarMVStruct
	**MVs = (MvarMVStruct **)
			       IritMalloc(sizeof(MvarMVStruct *) * NumOfMVs);

    /* Convert to a regular MV problem and solve. */
    for (i = 0; i < NumOfMVs; i++)
        MVs[i] = MvarExprTreeToMV(Eqns -> Eqns[i]);
    MvarMVUpdateConstDegDomains(MVs, NumOfMVs);

    MVPts = MvarMVsZeros0D(MVs, Eqns -> ConstraintTypes, NumOfMVs,
			   SubdivTol, NumericTol);

    for (i = 0; i < NumOfMVs; i++)
        MvarMVFree(MVs[i]);
    IritFree(MVs);

    return MvarZeroSolverSolutionNew(NULL, NULL, MVPts, 
				     MVAR_ZER_SLVR_SOLUTION_PT_LIST);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Extraction of the domain dimension of the MVs/ETs/other constraints in a   M
* zero finding problem.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:	The zero finding problem structure.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:  The dimension of the domain of the problem.	                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverGetDmnDim						     M
*****************************************************************************/
int MvarZeroSolverGetDmnDim(const MvarZeroPrblmStruct *Problem)
{
    switch (Problem -> ActiveRepresentation) {
        case MVAR_ZER_SLVR_MVS:
            return Problem -> U.MVs[0] -> Dim;
        case MVAR_ZER_SLVR_EXP_TREE:
	    return Problem -> U.Eqns -> Eqns[0] -> Dim;
        default:
	    assert(0);
	    MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
	    return -1;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* A preliminary recursive stage of the subdivision solver, invoked by the    M
* general solver - subdivides at the C1 discontinuity and solves recursively.M
*									     *
* PARAMETERS:                                                                M
*   Problem:	The zero finding problem structure.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroPrblmStruct **:  The array of two sub-problems.		     M
*									     *
* SEE ALSO:								     M
*   MvarZeroSolver							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroC1DiscontSubdiv						     M
*****************************************************************************/
static MvarZeroPrblmStruct **MvarZeroC1DiscontSubdiv(
						 MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    int i,
	JLoc = PrblmIntrnl -> C1DiscontInd,
	NumOfConstraints = Problem -> NumOfConstraints;
    CagdRType
	t = PrblmIntrnl -> C1DiscontParam;
    MvarMVStruct **MVs1, **MVs2,
	**MVs = Problem -> U.MVs;
    MvarZeroPrblmStruct
	**SubProblems = (MvarZeroPrblmStruct **)
		                IritMalloc(2 * sizeof(MvarZeroPrblmStruct *));
    			
    SubProblems[0] = SubProblems[1] = NULL;

    /* Subdivision at the discontinuity. */
    MVs1 = (MvarMVStruct **)
                        IritMalloc(sizeof(MvarMVStruct *) * NumOfConstraints);
    MVs2 = (MvarMVStruct **)
	                IritMalloc(sizeof(MvarMVStruct *) * NumOfConstraints);

    for (i = 0; i < NumOfConstraints; i++) {
	MVs1[i] = MvarMVSubdivAtParam(MVs[i], t, JLoc);
	MVs2[i] = MVs1[i] -> Pnext;
	MVs1[i] -> Pnext = NULL;
    }
#ifdef DEBUG
    if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]) {
	printf("\nC1 Discont!! Subdiv is at JLoc = %d, t = %.3f, into:\n", 
								     JLoc, t);
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
	printf("\n********************************************************\n");
    }
#endif /* DEBUG */

    SubProblems[0] = MvarZeroSolverSubProblem(Problem, MVs1, NULL, NULL);
    SubProblems[1] = MvarZeroSolverSubProblem(Problem, MVs2, NULL, NULL);
    return SubProblems;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Searches all MVs for parameter locations that are C1 discont.  If found  M
* update JLoc and t to this finding and returns TRUE. Same functionality as  M
* MvarSSIHasC1Discont, but for any number of MVs, not necessarily Dim - 1.   M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs: 	 The system of equations.				     M
*   NumOfMVs:    The number of equations.				     M
*   JLoc:        The direction in the domain with C1 discont, if has one.    M
*   t:           The parameter at the C1 discont. if has one.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:        TRUE if found a C1 discont., FALSE otherwise.          M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroHasC1Discont						     M
*****************************************************************************/
CagdBType MvarZeroHasC1Discont(MvarMVStruct * const *MVs,
			       int NumOfMVs,
			       int *JLoc,
			       CagdRType *t)
{
    int i, j,
	Dim = MVs[0] -> Dim;
    
    for (i = 0; i < NumOfMVs; i++) {
        MvarMVStruct const 
	    *MV = MVs[i];

	if (MVAR_IS_BEZIER_MV(MV))
	    continue;

        for (j = 0; j < Dim; j++) {
	    if (BspKnotC1Discont(MV -> KnotVectors[j],
				 MV -> Orders[j], MV -> Lengths[j], t)) {
	        /* Found a discontinuity - return it. */
	        *JLoc = j;
		return TRUE;
	    }
	}
    }
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Updtaes the problem with the domain and related data in the MVs case.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:	The zero finding problem structure.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroUpdateProblemDmnMVs                                              M
*****************************************************************************/
void MvarZeroUpdateProblemDmnMVs(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    int j, 
	JLoc = -1,
	Dim = Problem -> U.MVs[0] -> Dim;
    CagdRType *MinDmn, *MaxDmn,
	MaxSide = 0.0;

    PrblmIntrnl -> MVMinDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));
    PrblmIntrnl -> MVMaxDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));
    /* The original domain, only if we are here for the first time: */
    if (PrblmIntrnl -> SubdivDepth == 0) {
	PrblmIntrnl -> OrigMVMinDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));
	PrblmIntrnl -> OrigMVMaxDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));
	MvarMVDomain(Problem -> U.MVs[0],
		     PrblmIntrnl -> OrigMVMinDmn,
		     PrblmIntrnl -> OrigMVMaxDmn, -1);
    }

    MvarGetSubdivParamDomains(Problem -> U.MVs[0],
			      PrblmIntrnl -> MVMinDmn,
		              PrblmIntrnl -> MVMaxDmn, -1);
    MinDmn = PrblmIntrnl -> MVMinDmn;
    MaxDmn = PrblmIntrnl -> MVMaxDmn;

    for (j = 0; j < Dim; j++) {
        if (IRIT_ABS(MaxDmn[j] - MinDmn[j]) > MaxSide) {
	    MaxSide = IRIT_ABS(MaxDmn[j] - MinDmn[j]);
	    JLoc = j;
	}
    }
    PrblmIntrnl -> MaxSide = MaxSide;
    PrblmIntrnl -> MaxSideDir = JLoc;
    PrblmIntrnl -> UnderSubdivTol = MaxSide <= Problem -> SubdivTol;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Updates the problem with the domain and related data in the ETs case.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:	The zero finding problem structure.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroUpdateProblemDmnExpTr                                            M
*****************************************************************************/
void MvarZeroUpdateProblemDmnExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    int i, j,
	JLoc = -1,
	Dim = Problem -> U.Eqns -> Eqns[0] -> Dim;
    MvarExprTreeStruct
        **ETs = Problem -> U.Eqns -> Eqns;
    CagdRType *MinDmn, *MaxDmn,
	MaxSide = 0.0;

    PrblmIntrnl -> MVMinDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));
    PrblmIntrnl -> MVMaxDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));

    /* The original domain, only if we are here for the first time: */
    if (PrblmIntrnl -> SubdivDepth == 0) {
	PrblmIntrnl -> OrigMVMinDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));
	PrblmIntrnl -> OrigMVMaxDmn = (CagdRType *) 
					IritMalloc(Dim * sizeof(CagdRType));

	MinDmn = PrblmIntrnl -> OrigMVMinDmn;
	MaxDmn = PrblmIntrnl -> OrigMVMaxDmn;

	for (j = 0; j < ETs[0] -> Dim; j++) {
	    for (i = 0; i < Problem -> NumOfConstraints; i++)
	        if (MvarETDomain(ETs[i], &MinDmn[j], &MaxDmn[j], j))
		    break;
	    assert(i < Problem -> NumOfConstraints);/* 1 ET must have domain.*/
	}
    }

    MinDmn = PrblmIntrnl -> MVMinDmn;
    MaxDmn = PrblmIntrnl -> MVMaxDmn;

    for (j = 0; j < ETs[0] -> Dim; j++) {
        for (i = 0; i < Problem -> NumOfConstraints; i++)
	    if (MvarETDomain(ETs[i], &MinDmn[j], &MaxDmn[j], j))
	        break;
	assert(i < Problem -> NumOfConstraints);   /* 1 ET must have domain. */
    }

    for (j = 0; j < Dim; j++) {
        if (IRIT_ABS(MaxDmn[j] - MinDmn[j]) > MaxSide) {
	    MaxSide = IRIT_ABS(MaxDmn[j] - MinDmn[j]);
	    JLoc = j;
	}
    }
    PrblmIntrnl -> MaxSide = MaxSide;
    PrblmIntrnl -> MaxSideDir = JLoc;
    PrblmIntrnl -> UnderSubdivTol = MaxSide <= Problem -> SubdivTol;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Completes the construction of the MVs problem structure with the actions M
* that are performed only once: when attaining a C1 smooth problem for the   M
* first time, adds the gradients and the current MVs, later to be used in    M
* numeric step, and solves the problem on the boundary of the domain.        M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:	The zero finding problem structure.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:	TRUE if updates were successful, FALSE otherwise.	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroFirstSmoothUpdatesMVs                                            M
*****************************************************************************/
CagdBType MvarZeroFirstSmoothUpdatesMVs(MvarZeroPrblmStruct *Problem)
{
    CagdBType 
	CanExtendDmn = FALSE;
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;
    int i, l,
	NumOfZeroConstraints = Problem -> NumOfZeroConstraints;
    MvarMVStruct 
	**ExtendedDmnMVs = NULL;

    /* If any of the following is already updated, something is wrong. We    */
    /* are not visiting here for the first time.			     */
    if (PrblmIntrnl -> UnderSubdivTol == FALSE && 
	(PrblmIntrnl -> U.FirstSmoothMVs != NULL ||
	 PrblmIntrnl -> MVGradients != NULL ||
	 PrblmIntrnl -> SolutionsOnBoundary != NULL)) {
        MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
	return FALSE;
    }

    if (_MVGlblZeroDmnExtension > 0.0) {
	CanExtendDmn = MVAR_IS_BSPLINE_MV(Problem -> U.MVs[0]) && 
				       PrblmIntrnl -> ExpectedSolutionDim == 0;
	if (CanExtendDmn) {
	    CagdRType *Epsilons;

	    ExtendedDmnMVs = (MvarMVStruct **) IritMalloc(
			sizeof(MvarMVStruct *) * Problem -> NumOfConstraints);
	    Epsilons = IritMalloc(
			     sizeof(CagdRType) * Problem -> NumOfConstraints);
	    for (i = 0; i < Problem -> NumOfConstraints; i++) 
		Epsilons[i] = _MVGlblZeroDmnExtension * Problem -> SubdivTol;

	    for (i = 0; i < Problem -> NumOfConstraints; i++) {
		ExtendedDmnMVs[i] = MvarMVExtension(Problem -> U.MVs[i], 
						    NULL, NULL, Epsilons);
	    }
	    IritFree(Epsilons);
	}
    }

    /* If no domain extension is performed - simply copy the current MVs. */
    if (!CanExtendDmn) {
	/* Add the first smooth MVs: */
	PrblmIntrnl -> U.FirstSmoothMVs = (MvarMVStruct **) 
	      IritMalloc(sizeof(MvarMVStruct *) * Problem -> NumOfConstraints);

	for (i = 0; i < Problem -> NumOfConstraints; i++)
	    PrblmIntrnl -> U.FirstSmoothMVs[i] = 
					       MvarMVCopy(Problem -> U.MVs[i]);
    }
    else { /* The domain was extended. Use this as the first smooth MVs. */
	PrblmIntrnl -> U.FirstSmoothMVs = ExtendedDmnMVs;
    }

    switch (PrblmIntrnl -> ExpectedSolutionDim) {
	case 0:
	    /*   Several options for adding the gradients:                   */
	    /*   See if grads are to be stored along with the functions.     */
	    /*   A scalar function f of dim d will be converted to a 1+d     */
	    /*   vector function, holding original function f and the d      */
	    /*	 coefficients of the grad.				     */
	    /*   This 1+d storage will be used for normal cone tests, etc.   */
	    if (_MVGlblZeroApplyParallelHyperPlaneTest) {
		PrblmIntrnl -> MVGradients = (MvarMVGradientStruct **) 
		    IritMalloc(sizeof(MvarMVGradientStruct *) * Problem ->
			                                    NumOfConstraints);

		for (i = 0; i < NumOfZeroConstraints; i++) {
		    int Size = sizeof(CagdRType) * MVAR_CTL_MESH_LENGTH
			(Problem -> U.MVs[i]);
		    MvarMVGradientStruct
			*Grad = MvarMVPrepGradient(
				     PrblmIntrnl -> U.FirstSmoothMVs[i], TRUE);

		    assert(!MVAR_IS_RATIONAL_MV(Problem -> U.MVs[i]));

		    /*   Move the gradient vector to this MVs[i], in place.  */
		    /*   The MVs[i] was a scalar E1 function and now it will */
		    /*   become E(1+Dim).				     */
		    if (Problem -> U.MVs[i] -> Dim > MVAR_MAX_PT_COORD) {
			MVAR_FATAL_ERROR(MVAR_ERR_DIM_TOO_HIGH); 
			return FALSE;
		    }
		    for (l = 0; l < Problem -> U.MVs[i] -> Dim; l++) {
			Problem -> U.MVs[i] -> Points[l + 2] =
			                        (CagdRType *) IritMalloc(Size);
			IRIT_GEN_COPY(Problem -> U.MVs[i] -> Points[l + 2],
				      Grad -> MVGrad -> Points[l + 1], Size);

			Problem -> U.MVs[i] -> PType =
			    MVAR_MAKE_PT_TYPE(FALSE,
					      Problem -> U.MVs[i] -> Dim + 1);
		    }
		    /* Keep for the numeric marching step. */
		    PrblmIntrnl -> MVGradients[i] = Grad;
		}
	    }
	    else { /* Prepare the grads normally and store them. */
		PrblmIntrnl -> MVGradients = (MvarMVGradientStruct **) 
		    IritMalloc(sizeof(MvarMVGradientStruct *) * Problem ->
			                               NumOfZeroConstraints);
		for (i = 0; i < Problem -> NumOfZeroConstraints; i++) {
		    PrblmIntrnl -> MVGradients[i] = 
			MvarMVPrepGradient(PrblmIntrnl -> U.FirstSmoothMVs[i],
					   FALSE);
		}
	    }
	    break;
	case 1:	/* Same treatment as in the 2D case. */
	case 2:
	    PrblmIntrnl -> MVGradients = (MvarMVGradientStruct **) 
		IritMalloc(sizeof(MvarMVGradientStruct *) *
			                     Problem -> NumOfZeroConstraints);
	    for (i = 0; i < Problem -> NumOfZeroConstraints; i++) {
		PrblmIntrnl -> MVGradients[i] = 
		    MvarMVPrepGradient(Problem -> U.MVs[i], FALSE);
	    }
	    break;
    }
    MVAR_ZERO_SLVR_APPLY(SolveBoundary)(Problem);
    PrblmIntrnl -> ConstructionComplete = TRUE;
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Completes the construction of the ETs problem structure with the actions M
* that are performed only once: when attaining a C1 smooth problem for the   M
* first time. Adds the gradients and the current MVs, later to be used in    M
* numeric step, and solves the problem on the boundary of the domain.        M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:	The zero finding problem structure.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:	TRUE if updates were successful, false otherwise.	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroFirstSmoothUpdatesExpTr                                          M
*****************************************************************************/
CagdBType MvarZeroFirstSmoothUpdatesExpTr(MvarZeroPrblmStruct *Problem)
{
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the tolerance (or zero to disable) of the domain extension inside   M
* the multivariate subdivisions' zero set solver.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   DmnExt:   New setting for domain extension usage.	                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:       Old extensions tolerance for domain extension usage.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver, MvarMVsZerosDomainReduction,			     M
*   MvarMVsZerosGradPreconditioning					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosDmnExt	                                                     M
*****************************************************************************/
CagdRType MvarMVsZerosDmnExt(CagdRType DmnExt)
{
    CagdRType
        OldVal = _MVGlblZeroDmnExtension;

    _MVGlblZeroDmnExtension = DmnExt;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construct a point of the dimension as the given problem in the middle of M
* its parametric domain.                                                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:    To construct a point in the middle of its domain.            M
*   SingleSol:  If TRUE, this point is a single solution it MV domain.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:  The construct point in the middle of MV.                M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros                                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroGenPtMidDmn                                                      M
*****************************************************************************/
MvarPtStruct *MvarZeroGenPtMidDmn(const MvarZeroPrblmStruct *Problem,
				  int SingleSol)
{
    int l;
    MvarPtStruct
        *Pt = MvarPtNew(MvarZeroSolverGetDmnDim(Problem));
    const MvarMVStruct *MV;
    const MvarExprTreeStruct *ET;

    switch (Problem -> ActiveRepresentation) {
        case MVAR_ZER_SLVR_MVS:
	    MV = Problem -> U.MVs[0];
	    for (l = 0; l < MV -> Dim; l++) {
		CagdRType Min, Max;

		MvarGetSubdivParamDomains(MV, &Min, &Max, l);
		Pt -> Pt[l] = (Min + Max) * 0.5;
	    }
	    break;
        case MVAR_ZER_SLVR_EXP_TREE:
	    ET = Problem -> U.Eqns -> Eqns[0];
	    for (l = 0; l < ET -> Dim; l++) {
		CagdRType Min, Max;

		MvarETDomain(ET, &Min, &Max, l);
		Pt -> Pt[l] = (Min + Max) * 0.5;
	    }
	    break;
        default:
	    assert(0);
	    MVAR_FATAL_ERROR(MVAR_ERR_ZER_PRBLM_CNSTRCT);
	    return NULL;
    }

#   ifdef DEBUG_DUMP_DOMAINS
    fprintf(GlblDumpDomainsFile, "[OBJECT [RGB \"255,100,100\"] [Gray 0.5] [Width 0.02] NONE [CTLPT E%d", MV -> Dim);
    for (l = 0; l < MV -> Dim; l++)
	fprintf(GlblDumpDomainsFile, " %f", Pt -> Pt[l]);
    fprintf(GlblDumpDomainsFile, "]\n]\n");
#   endif /* DEBUG_DUMP_DOMAINS */

    AttrSetIntAttrib(&Pt -> Attr, "SingleSol", SingleSol);

    return Pt;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Apply a numerical improvement stage, as a first order minimization       M
* procedure of gradient computation and marching, in place.                  M
*                                                                            *
* PARAMETERS:                                                                M
*   ZeroPt:	  Approximated solution, derived from a subdivision process, M
*		  to improve in place.					     M
*   Eqns:         The constraints are given as Equations, if not NULL.       M
*   MVs:          Alternatively, the constraints are given as MVS.           M
*   NumMVs:       If MVs is not NULL, this specifies size of the MVs vector. M
*   NumericTol:   Tolerance of the numerical process.  Tolerance is measured M
*		  in the deviation of the scalar multivariates from their    M
*		  equality. Inequalities are ignored here.  Points that fail M
*		  to improve numerically are purged away.		     M
*   InputMinDmn, InputMaxDmn:  Optional domain restriction. Can be NULL in   M
*                 which, the MVS/Eqns domains are used.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   List of points on the solution set.  Dimension of the  M
*		      points will be the same as the dimensions of all MVs.  M
*		      Points that failed to improve are purged away.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarExprTreesZeros, MvarExprTreeEqnsZeros                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZero0DNumeric                                                        M
*****************************************************************************/
MvarPtStruct *MvarZero0DNumeric(MvarPtStruct *ZeroPt,
				const MvarExprTreeEqnsStruct *Eqns,
				MvarMVStruct const * const *MVs,
				int NumMVs,
				CagdRType NumericTol,
				const CagdRType *InputMinDmn,
				const CagdRType *InputMaxDmn)
{
    IRIT_STATIC_DATA int
        AllocDim = 0,
        AllocNumMVs = 0;
    IRIT_STATIC_DATA CagdRType
        *A = NULL,
        *x = NULL,
        *x2 = NULL,
        *b = NULL,
        *MinDmn = NULL,
        *MaxDmn = NULL;
    MvarPtStruct
	*ImprovedSet = NULL;
    int i, j, Dim,
        Count = 0,
        GoodMoves = 0;
    CagdRType NewRngErr, NewDmnErr,
        DmnErr = IRIT_INFNTY,
        RngErr = IRIT_INFNTY,
        RangeDiff = 0.0,
        StepSize = 1.0;

    NumericTol = IRIT_FABS(NumericTol);

    if (Eqns != NULL) {
        Dim = Eqns -> Eqns[0] -> Dim;
    }
    else {
        assert(MVs != NULL);
	Dim = MVs[0] -> Dim;
    }

    /* Make sure this is not an over-constrained system of equations. */
    if (Dim < NumMVs || NumMVs == 0)
	return ZeroPt;

    /* Make sure this is not an over-constrained system of equations. */
    if (AllocDim < Dim || AllocNumMVs < NumMVs) {
        if (A != NULL) {
	    IritFree(A);
	    IritFree(x);
	    IritFree(x2);
	    IritFree(b);
	    IritFree(MinDmn);
	    IritFree(MaxDmn);
	}
	AllocDim = 2 * Dim;
	AllocNumMVs = 2 * NumMVs;

	A = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim *
				                                 AllocNumMVs);
	x = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
	x2 = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
	b = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
	MinDmn = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
	MaxDmn = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
    }

    /* Get the domain of the constraints. */
    if (InputMinDmn != NULL && InputMaxDmn != NULL) {
        CAGD_GEN_COPY(MinDmn, InputMinDmn, sizeof(CagdRType) * Dim);
        CAGD_GEN_COPY(MaxDmn, InputMaxDmn, sizeof(CagdRType) * Dim);
    }
    else if (Eqns != NULL) {
        for (j = 0; j < Dim; j++) {
	    for (i = 0; i < NumMVs; i++) {
	        if (MvarETDomain(Eqns -> Eqns[i], &MinDmn[j], &MaxDmn[j], j))
		    break;
	    }
	    assert(i < NumMVs);
	}
    }
    else
        MvarMVDomain(MVs[0], MinDmn, MaxDmn, -1);

    do { /* The numeric iterations! */
#	ifdef DEBUG
        {
	    IRIT_IF_DEBUG_ON_PARAMETER(_DebugZeroZeroErr) {
	        if (DmnErr == IRIT_INFNTY && RngErr == IRIT_INFNTY)
		    IRIT_INFO_MSG("\n******************************\n");
		printf("\nRngErr=%.7g, DmnErr=%.7g, StpSz=%.7g at:\n\t",
		       RngErr, DmnErr, StepSize);
		for (i = 0; i < Dim; i++)
		    printf(" %.7g", ZeroPt -> Pt[i]);
		printf("\nEquations:\n");
	    }
	}
#	endif /* DEBUG */

	/* We now seek a point that is on each tangent hyperplanes.         */
	/* Hence, we have Dim linear equations and Dim dofs.                */
	for (i = 0; i < NumMVs; i++) {
	    /* Derive tangent hyperplane for multivariate at current pos. */
	    MvarPlaneStruct
	        *Pln = Eqns != NULL ?
	             MvarExprTreeEvalTanPlane(Eqns -> Eqns[i], ZeroPt -> Pt) :
		     MvarMVEvalTanPlane(MVs[i], ZeroPt -> Pt);

	    /* Copy the constraint into the matrix form. */
	    CAGD_GEN_COPY(&A[i * Dim], Pln -> Pln, sizeof(CagdRType) * Dim);
	    b[i] = -Pln -> Pln[Dim + 1];

#	    ifdef DEBUG
	    {
	        IRIT_IF_DEBUG_ON_PARAMETER(_DebugZeroZeroErr) {
		    for (j = 0; j < Dim; j++)
		        printf(" %.7g", Pln -> Pln[j]);
		    printf(" = %.7g\n", b[i]);
		}
	    }
#	    endif /* DEBUG */

	    /* Free the computed plane. */
	    MvarPlaneFree(Pln);
	}

	/* Solve, possibly under constrained, system of lin. equations. */
	if (IritQRUnderdetermined(A, NULL, NULL, NumMVs, Dim))
	    break;

	/* Add current position to b vector, seeking relative solution. */
	for (i = 0; i < NumMVs; i++) {
	    int j;
	    CagdRType
	        t = 0.0;

	    for (j = 0; j < Dim; j++)
	        t += A[i * Dim + j] * ZeroPt -> Pt[j];

	    b[i] -= t;
	}

	/* Solve for the new multivariate(s) coefficients. */
	IritQRUnderdetermined(NULL, x, b, NumMVs, Dim);

	/* And add the current location back, getting absolute result. */
	for (i = 0; i < Dim; i++)
	    x[i] = ZeroPt -> Pt[i] + x[i] * StepSize;

#	ifdef DEBUG
	{
	    IRIT_IF_DEBUG_ON_PARAMETER(_DebugZeroZeroErr) {
	        printf("\tNew solution: ");
		for (i = 0; i < Dim; i++)
		    printf(" %f", x[i]);
		printf("\n");
	    }
	}
#	endif /* DEBUG */

	for (i = 0; i < Dim; i++) {
	    /* Make sure we stay within our domain bounds. */
	    x[i] = IRIT_BOUND(x[i], MinDmn[i], MaxDmn[i]);
	}

	/* Now the million dollar question - is x a better solution!? */
	NewRngErr = MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumMVs, x);

#ifdef MVAR_ZERO_QUAD_SOLVER
	/* If convergence is slowed down, try a quadratic step. */
	if (NewRngErr != 0.0 && RngErr / NewRngErr < 2.0) {
	    int NumSols;
	    CagdRType NewRngErr2, NewRngErr3, A, B, C, Max, Sols[2];

	    for (i = 0; i < Dim; i++)
	        x2[i] = (x[i] + ZeroPt -> Pt[i]) * 0.5;
	    NewRngErr2 = MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumMVs, x2);

	    /* Build a quadratic fit through ZeroPt(t=0), x2(t=1/2), x(t=1) */
	    /* As "At^2 + Bt + C = 0".				            */
	    A = 2 * RngErr - 4 * NewRngErr2 + 2 * NewRngErr;
	    B = -3 * RngErr + 4 * NewRngErr2 - NewRngErr;
	    C = RngErr;
	    Max = IRIT_MAX(IRIT_MAX(IRIT_FABS(A), IRIT_FABS(B)),
			   IRIT_FABS(C));

	    if (IRIT_APX_EQ(A / Max, 0) && IRIT_APX_EQ(A / Max, 0)) {
	        /* Examine the constant equation case. */
	        if (IRIT_APX_EQ(B / Max, 0))
		    NumSols = -1;       /* Use the solution found before... */
		else {
		    /* Actually we have a linear function here... */
		    Sols[0] = -C / B;
		    NumSols = 1;
		}
	    }
	    else {
	        B /= A;
		C /= A;

		NumSols = GMSolveQuadraticEqn(B, C, Sols);
	    }

	    if (NumSols == 2) {
	        int j;
		CagdRType Vals[2];

		for (j = 0; j < 2; j++) {
		    for (i = 0; i < Dim; i++) {
		        x2[i] = ZeroPt -> Pt[i] +
			                   (x[i] - ZeroPt -> Pt[i]) * Sols[j];
			x2[i] = IRIT_BOUND(x2[i], MinDmn[i], MaxDmn[i]);
		    }

		    Vals[j] = MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumMVs, x2);
		}
		j = Vals[0] > Vals[1];
		if (Vals[j] < NewRngErr) {
		    /* We are improving with this quadratic step. */
		    for (i = 0; i < Dim; i++) {
		        x[i] = ZeroPt -> Pt[i] +
			                   (x[i] - ZeroPt -> Pt[i]) * Sols[j];
			x[i] = IRIT_BOUND(x[i], MinDmn[i], MaxDmn[i]);
		    }
		    NewRngErr = Vals[j];
		}
	    }
	    else {
	        if (NumSols == 0) {
		    /* Find minimum of parabola as the solution. */
		    Sols[0] = -B / 2;
		    NumSols = 1;
		}

		/* So now we typically have one solution: */
		if (NumSols == 1) {
		    for (i = 0; i < Dim; i++) {
		        x2[i] = ZeroPt -> Pt[i] +
			                   (x[i] - ZeroPt -> Pt[i]) * Sols[0];
			x2[i] = IRIT_BOUND(x2[i], MinDmn[i], MaxDmn[i]);
		    }
		}

		/* Compute the new location's error. If middle location */
		/* has smaller error, use that instead.		    */
		NewRngErr3 = MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumMVs, x2);
		if (NewRngErr > NewRngErr2 && NewRngErr > NewRngErr3) {
		    for (i = 0; i < Dim; i++)
		        x[i] = x2[i];
		    NewRngErr = NewRngErr2;
		}		        
	    }
	}
#endif /* MVAR_ZERO_QUAD_SOLVER */

	/* Estimate the domain's error as the difference between last       */
	/* parametric location and this one.			            */
	NewDmnErr = 0.0;
	for (i = 0; i < Dim; i++)
	    NewDmnErr += IRIT_FABS(ZeroPt -> Pt[i] - x[i]);

#	ifdef DEBUG
	{
	    IRIT_IF_DEBUG_ON_PARAMETER(_DebugZeroZeroErr) {
	        IRIT_INFO_MSG_PRINTF(
			    "Rng Err = %.7g (%.7g), Dmn Err = %.7g (%.7g)",
			    NewRngErr, RngErr, NewDmnErr, DmnErr);
		for (i = 0; i < Dim; i++)
		    printf(" %.7g", x[i]);
		printf("\n");
	    }
	}
#	endif /* DEBUG */

	if (NewRngErr < RngErr || NewRngErr == 0.0 || NewDmnErr == 0.0) {
	    RangeDiff = IRIT_FABS(RngErr - NewRngErr);
	    RngErr = NewRngErr;
	    DmnErr = NewDmnErr;
	    CAGD_GEN_COPY(ZeroPt -> Pt, x, sizeof(CagdRType) * Dim);

	    /* Increase step size if we have a sequence of good moves. */
	    if (GoodMoves++ > 5) {
	        StepSize = IRIT_MAX(1.0, StepSize * 2.0);
		GoodMoves = 0;
	    }
	}
	else {
	    if (StepSize < 0.5) {
	        /* Find a nearby point - perturb the solution pt. */
	        if (MvarEqnsPerturbZero(ZeroPt, RngErr, Eqns, MVs, NumMVs,
					MinDmn, MaxDmn, RngErr * StepSize))
		    StepSize *= 2.0;
		else
		    StepSize *= 0.5;
	    }
	    else {
	        StepSize *= 0.5;
		GoodMoves = 0;
	    }
	}

	/* When to terminate is a difficult question.  We clearly desire    */
	/* that the RngErr will be below the NumericTol but we also         */
	/* desire low DmnErr (in tangency intersections it might be much    */
	/* larger than the RngErr), so for the 1st few iterations ask       */
	/* for small DmnErr as well.				            */
    }
    while ((Count++ < 10 && DmnErr > NumericTol) ||
	   (RngErr > NumericTol &&
	    DmnErr > 0.0 &&
	    Count < MVAR_NUMER_ZERO_NUM_STEPS));

#   ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugZeroNumerZeroIters, FALSE) {
	    IRIT_INFO_MSG_PRINTF(
		    "\tNumeric step %ssuccessful after %d iters (%.8g %.8g)\n",
		    Count >= MVAR_NUMER_ZERO_NUM_STEPS ? "un" : "",
		    Count, DmnErr, RngErr);
	}
    }
#   endif /* DEBUG */

    if (RngErr > NumericTol) {
        MvarPtFree(ZeroPt);
	return NULL;
    }
    else {
        AttrSetRealAttrib(&ZeroPt -> Attr, "RngError", RngErr);
	IRIT_LIST_PUSH(ZeroPt, ImprovedSet);
	return ZeroPt;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the L1 error of the current position.                           *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:        The multivariates to evaluate error for.                     *
*   Params:     The location where the error is to be evaluated.             *
*   NumOfMVs:   Number of multivariates we have.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:     Computed error.                                           *
*****************************************************************************/
CagdRType MvarMVEvalErrorL1(MvarMVStruct const * const *MVs,
			    CagdRType *Params,
			    int NumOfMVs)
{
    int i;
    CagdRType
	Err = 0.0;

    for (i = 0; i < NumOfMVs; i++) {
	CagdRType
	    *R = MvarMVEval(MVs[i], Params);

	assert(!MVAR_IS_RATIONAL_MV(MVs[i]));

	Err += IRIT_FABS(R[1]);
    }

    return Err;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the L1 error of the current position.                           *
*                                                                            *
* PARAMETERS:                                                                *
*   Eqns:        The ETs to evaluate error for.		                     *
*   Params:     The location where the error is to be evaluated.             *
*   NumOfMVs:   Number of multivariates we have.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:     Computed error.                                           *
*****************************************************************************/
CagdRType MvarExprTreeEqnsEvalErrorL1(const MvarExprTreeEqnsStruct *Eqns,
				      CagdRType *Params)
{
    int i;
    CagdRType
	Err = 0.0;

    for (i = 0; i < Eqns -> NumEqns; i++) {
	CagdRType
	    *R = MvarExprTreeEval(Eqns -> Eqns[i], Params);

	Err += IRIT_FABS(R[1]);
    }

    return Err;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Look for a near by location that improves the solution.  This naive      *
* scheme is employed when the regular NR approach failed.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt:           To try and locally improve, in place.                      *
*   StartRngErr:  Current range error of parameteric location Pt.            *
*   Eqns:         The constraints are given as Equations, if not NULL.       *
*   MVs:          Alternatively, the constraints are given as MVS.           *
*   NumEqns:      Size of theEqns or  MVs vector.			     *
*   MinDmn, MaxDmn:  Domain of constraints.				     *
*   Step:         To take in perturbation.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:                                                                     *
*****************************************************************************/
int MvarEqnsPerturbZero(MvarPtStruct *Pt,
			CagdRType StartRngErr,
			const MvarExprTreeEqnsStruct *Eqns,
			MvarMVStruct const * const *MVs,
			int NumEqns,
			CagdRType *MinDmn,
			CagdRType *MaxDmn,
			CagdRType Step)
{
    int i, n,
        Dim = Eqns != NULL ? Dim = Eqns -> Eqns[0] -> Dim : MVs[0] -> Dim;
    CagdRType RngErr,
        x[MVAR_MAX_PT_SIZE];

    CAGD_GEN_COPY(x, Pt -> Pt, sizeof(CagdRType) * Dim);

    for (n = 0; n < 3; n++) {
        Step *= 0.5;

        /* Try +/- perturbations of Step size, in all Dim directions: */
	for (i = 0; i < Dim; i++) {
	    x[i] += Step;
	    if (x[i] > MaxDmn[i])
		x[i] = MaxDmn[i];
	    if ((RngErr = MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumEqns, x))
	                                                       > StartRngErr) {
	        x[i] -= Step * 2;
		if (x[i] < MinDmn[i])
		    x[i] = MinDmn[i];
		RngErr = MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumEqns, x);
	    }

	    if (RngErr < StartRngErr) {
		CAGD_GEN_COPY(Pt -> Pt, x, sizeof(CagdRType) * Dim);
		return TRUE;
	    }
	    else
	        x[i] = Pt -> Pt[i];
	}
    }

    return FALSE;
}

#ifdef DEBUG

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Enables to trace the subdivision tree by printing the domains and the    *
* coefficients of the MVs, but with control over which solution space dim    *
* the user wishes to debug. Supports MVs representation only.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	    The zero finding problem structure.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroDBGDmn(MvarZeroPrblmStruct *Problem)
{
    if (Problem == NULL)
	return;
    else {
	MvarZeroPrblmIntrnlStruct
	    *PrblmIntrnl = Problem -> _Internal;

	if (Problem -> ActiveRepresentation != MVAR_ZER_SLVR_MVS ||
	    !_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim])
	    return;
	else {    /* Domain debugging is required for the current problem. */
	    int k, i, 
		Dim = Problem -> U.MVs[0] -> Dim;
	    CagdRType t1, t2;

	    _DebugZeroSolverDmnNum++; 

	    printf("\nDomain num: %d is:\n", _DebugZeroSolverDmnNum);
	    for (i = 0; i < Problem -> U.MVs[0] -> Dim; i++) {
		t1 = PrblmIntrnl -> MVMinDmn[i];
		t2 = PrblmIntrnl -> MVMaxDmn[i];
		printf("[%6.3f, %6.3f] ", t1, t2);
	    }
	    if (_DebugZeroSolverBySolDim[PrblmIntrnl -> ExpectedSolutionDim]
		                                                       == 2) {
		for (i = 0; i < Problem -> NumOfConstraints; i++) {
		    printf("\nControl coeffs of MVs[%d]:\n", i);
		    for (k = 0; k < Problem -> U.MVs[i] -> SubSpaces[Dim]; k++) {
			printf("[%.4f] ", Problem -> U.MVs[i] -> Points[1][k]);
		    }
		    printf("\n");
		}
	    }
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prints a point list.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   PtList:   The point list to print.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
void MvarZR1DDbgPrintPtList(const MvarPtStruct * const PtList)
{
    int i, 
	PtCnt = 1,
	Dim = PtList -> Dim;
    const MvarPtStruct 
	*PtIter = PtList;

    printf("\n");
    while (PtIter != NULL) {
	printf("[%d] ", PtCnt++);
	for (i = 0; i < Dim; i++)
	    printf("%.4f  ", PtIter -> Pt[i]);
	printf("\n");
	PtIter = PtIter -> Pnext;
    }
    printf("\n");
}

#endif /* DEBUG */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the slots required for a sub-domain info' structure.           M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints, defined in the         M
*		  required sub-domain.					     M
*   ProjDir1:     First coordinate direction of the two directions w.r.t.    M
*		  which the IPT succeeded.				     M
*   ProjDir2:     Second coordinate direction of the two directions w.r.t.   M
*		  which the IPT succeeded.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroSubDmnInfoStruct *:   The new info' structure.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSubDmnInfoStructFree				                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarSubDmnInfoStructNew                                                  M
*****************************************************************************/
MvarZeroSubDmnInfoStruct *MvarSubDmnInfoStructNew(MvarMVStruct **MVs, 
						  MvarMVDirType ProjDir1,
						  MvarMVDirType ProjDir2)
{
    MvarZeroSubDmnInfoStruct 
	*InfoStruct = (MvarZeroSubDmnInfoStruct *) 
				  IritMalloc(sizeof(MvarZeroSubDmnInfoStruct));

    InfoStruct -> MVs = MVs;
    InfoStruct -> ProjDirs[0] = ProjDir1;
    InfoStruct -> ProjDirs[1] = ProjDir2;
    return InfoStruct;
}
/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deallocates a sub-domain info' structure.			             M
*                                                                            *
* PARAMETERS:                                                                M
*   InfoStruct:   The info' structure to be freed.	                     M
*   NumOfMVs:     Number of multi-variates held by the info' structure.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   void.					                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSubDmnInfoStructNew				                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarSubDmnInfoStructFree                                                 M
*****************************************************************************/
void MvarSubDmnInfoStructFree(MvarZeroSubDmnInfoStruct *InfoStruct,
			      int NumOfMVs)
{
    int i;

    for (i = 0; i < NumOfMVs; i++)
	MvarMVFree(InfoStruct -> MVs[i]);
    IritFree(InfoStruct -> MVs);
    IritFree(InfoStruct);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tests if all coefficients of a multi-variate are zero, up to prescribed  M
* tolerance.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:		The (scalar valued) multi-variate to be tested. 	     M
*   NumericTol: The tolerance under which a coefficient is considered zero.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType: TRUE if all zeros, FALSE otherwise.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroMVConstraintFail				                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverIsMVZero                                                   M
*****************************************************************************/
CagdBType MvarZeroSolverIsMVZero(MvarMVStruct *MV,
				 CagdRType NumericTol)
{
    CagdBType 
	IsZeroMV = TRUE;
    int i;

    int Length = MVAR_CTL_MESH_LENGTH(MV);
    for (i = 0; i < Length; i++) {
	IsZeroMV &= IRIT_APX_EQ_EPS(MV -> Points[1][i], 0.0, NumericTol);
	if (!IsZeroMV)
	    break;
    }
    return IsZeroMV;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Verifies that the given array of multi-variates are all B-spline, and    *
* converts to B-spline, in place, in case of Bezier type.                    *
*									     *
* PARAMETERS:                                                                *
*   MVs:	The array of multi-variates.				     *
*   NumOfMVs:	The number of multi-variates.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroVerifyBspConstraints(MvarMVStruct **MVs,
			                 int NumOfMVs)

{
    int i;
    MvarMVStruct *TmpMV;

    for (i = 0; i < NumOfMVs; i++) {
	if (MVs[i] -> GType == MVAR_BEZIER_TYPE) {
	    TmpMV = MvarCnvrtBzr2BspMV(MVs[i]);
	    MvarMVFree(MVs[i]);
	    MVs[i] = TmpMV;
	}
	else {
	    assert(MVs[i] -> GType == MVAR_BSPLINE_TYPE);
	}
    }
}
