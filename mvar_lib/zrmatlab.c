/******************************************************************************
* ZrMatlab.c - interface for zero set finding problems of various             * 
*              dimensions and representations via matlab environment.         *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, Jan. 14.					      *
******************************************************************************/

#include "mvar_loc.h"

static void MvarZeroBzr2GenDmnMapPt(MvarPtStruct *Pt, 
				    CagdRType *MinDmn,
				    CagdRType *MaxDmn);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Convert a zero finding prolem from matlab form to the solver's problem   M
* and solve.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   Eqns:       The constraints as recieved after the matlab parsing.        M
*   NumOfEqns:  The number of constraints.			             M
*   MaxVarsNum: The maximal number of unknowns appearing in the problem.     M
*   MinDmn:	The min end point of the domain in all directions. If NULL,  M
*		considered as all zeros.				     M
*   MaxDmn:	The max end point of the domain in all directions. If NULL,  M
*		considered as all ones.				             M
*   NumericTol: The required numeric tolerance of the solution.		     M
*   SubdivTol:  The subdivision tolerance.				     M
*   StepTol:    The step size for numeric tracing of curves.		     M
*   Constraints: A vector of constraints specifying equality/inequality.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarZeroSolutionStruct *:  The solution to the problem.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver					                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolveMatlabEqns                                                  M
*****************************************************************************/
MvarZeroSolutionStruct *MvarZeroSolveMatlabEqns(
					      MvarMatlabEqStruct **Eqns,
					      int NumOfEqns,
					      int MaxVarsNum,
					      CagdRType *MinDmn,
					      CagdRType *MaxDmn,
					      CagdRType NumericTol,
					      CagdRType SubdivTol,
					      CagdRType StepTol,
					      MvarConstraintType *Constraints)
{
    int i, j, k, CurMonomsNum, *CurPowersMat, *Lengths, *CurMaxPowers, 
	LinIndex, CurVarPower,
	Dim = MaxVarsNum;
    CagdBType
	DmnChange = FALSE;
    CagdRType *CurCoeffs;
    MvarPtStruct *PtIter;
    MvarPolylineStruct *PolyIter;
    MvarMVStruct *PwrMV, *TmpMV, **BzrMVs;
    MvarZeroPrblmStruct *Problem;
    MvarZeroSolutionStruct *Sol;

#ifdef DEBUG
    fprintf(stderr, "Now in MvarZeroSolveMatlabEqns");
#endif /* DEBUG */

    BzrMVs = (MvarMVStruct **) IritMalloc(NumOfEqns * sizeof(MvarMVStruct *));
    Lengths = (int *) IritMalloc(Dim * sizeof(int));

    /* Create the multivariates. First power basis, then Bezier: */
    for (i = 0; i < NumOfEqns; i++) {
	CurPowersMat = Eqns[i] -> PowersMat;
	CurCoeffs = Eqns[i] -> CoeffArr;
	CurMonomsNum = Eqns[i] -> NumOfMonoms;
	CurMaxPowers = Eqns[i] -> MaxPowers;

	for (j = 0; j < Dim; j++)
	    Lengths[j] = CurMaxPowers[j] + 1;

	PwrMV = MvarMVNew(Dim, MVAR_POWER_TYPE, MVAR_PT_E1_TYPE, Lengths);

	/* Remember that the monomials in the Matlab equation are not        */
	/* ordered according to their order in the power basis. To set       */
	/* the coefficients - obtain the current multi index and convert.    */
	for (k = 0; k < CurMonomsNum; k++) {
	    LinIndex = 0;
	    for (j = 0; j < Dim; j++) {
		CurVarPower = CurPowersMat[k * Dim + j];
		LinIndex += CurVarPower * PwrMV -> SubSpaces[j];
	    }
	    PwrMV -> Points[1][LinIndex] = CurCoeffs[k];
	}
	BzrMVs[i] = MvarCnvrtPwr2BzrMV(PwrMV);
	MvarMVFree(PwrMV);
    }

    /* If the required domain is not [0,1], edit the MVs as required: */
    for (j = 0; j < Dim; j++) {
	if (!IRIT_APX_EQ_EPS(MinDmn[j], 0.0, NumericTol) ||
			      (!IRIT_APX_EQ_EPS(MaxDmn[j], 1.0, NumericTol))) {
	    DmnChange = TRUE;
	    for (i = 0; i < NumOfEqns; i++) {
		TmpMV = MvarBzrMVRegionFromMV(BzrMVs[i],
					      MinDmn[j], MaxDmn[j], j);
		MvarMVFree(BzrMVs[i]);
		BzrMVs[i] = TmpMV;
	    }
	}
    }

    /* Create the problem and solve. */
#ifdef DEBUG
    printf("Creating problem and solving");
#endif /* DEBUG */

    Problem = MvarZeroSolverPrblmNew((const struct MvarMVStruct * const *) 
									BzrMVs,
				     NULL, NumOfEqns, Constraints, SubdivTol,
				     NumericTol, StepTol, FALSE);
    Sol = MvarZeroSolver(Problem);

    /* The solutions are in [0,1]. Map to the required domain if different. */
    if (Sol != NULL && DmnChange) {
	switch (Sol -> ActiveRepresentation) {
	    case MVAR_ZER_SLVR_SOLUTION_PT_LIST:
		PtIter = Sol -> U.Pt;
		while (PtIter != NULL) {
		    MvarZeroBzr2GenDmnMapPt(PtIter, MinDmn, MaxDmn);
		    PtIter = PtIter -> Pnext;
		}
		break;
	    case MVAR_ZER_SLVR_SOLUTION_POLYLINE:
		PolyIter = Sol -> U.Pl;
		while (PolyIter != NULL) {
		    PtIter = PolyIter -> Pl;
		    while (PtIter != NULL) {
			MvarZeroBzr2GenDmnMapPt(PtIter, MinDmn, MaxDmn);
			PtIter = PtIter -> Pnext;
		    }
		    PolyIter = PolyIter -> Pnext;
		}
		break;
	    default:
		assert(0);
		break;
	}
    }
    
    MvarZeroSolverPrblmFree(Problem);
    for (i = 0; i < NumOfEqns; i++) {
	MvarMVFree(BzrMVs[i]);
    }
    IritFree(BzrMVs);
    IritFree(Lengths);
    return Sol;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Maps a point from the domain [0,1]^Dim to the given domain, preserving   *
* the weights of the endpoints in each dimension (an affine map between the  *
* domains).								     *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt:	    The point in [0,1]^Dim					     *
*   MinDmn: The min endpoints of the new domain.			     *
*   MaxDmn: The max endpoints of the new domain.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroBzr2GenDmnMapPt(MvarPtStruct *Pt, 
				    CagdRType *MinDmn,
				    CagdRType *MaxDmn)
{
    int i,
	Dim = Pt -> Dim;

    for (i = 0; i < Dim; i++) 
	Pt -> Pt[i] = MinDmn[i] + Pt -> Pt[i] * (MaxDmn[i] - MinDmn[i]);
}
