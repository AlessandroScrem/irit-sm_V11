/******************************************************************************
* Control.c - analyze control functions algenraically.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Nov 13.					      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "mvar_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes a cycle of length CycLen of bouncing lines between the given    M
* scalar curve C(t) and the diagonal of the domain.			     M
*   A ray is bounced 'down' from the diagonal to C and from C, we bounce     M
* a horizontal ray to the diagonal.  Let Di, i == 1,k be the diagonal points M
* and Ci, i == 1,k the points on C.  Then, the cycle of length k = CycLen    M
* is:  (D1, C1, D2, C2, ... Dk, Ck).  The following algebraic constraints    M
* can be imposed:							     M
*    Ci(y) = Di+1(y), for all i, (horizontal move from curve to diagonal)    V
*    Di(x) = Ci(x),   for all i. (vertical move from diagonal to curve)      V
* Over all, we have 2n equations and 2n dofs.  We can reduce them to n       M
* equations and dofs (eliminating the diagonals):			     M
*    Ci(y) = Ci+1(x),   for all i.					     V
* Over all, we have now n equations and n dofs (n different parameters of C).M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:          A curve to compute a cycle of length CycLen for.           M
*   CycLen:       Length of sought cycle (also k above).	             M
*   SubdivTol:    Tolerance of the solution.  This tolerance is measured in  M
*		  the parametric space of the curves.			     M
*   NumerTol:     Numeric tolerance of a possible numeric improvment stage.  M
*		  The numeric stage is employed if NumericTol < SubdivTol.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:  Linked list of solutions, each holding the parameter    M
*		values of the different Cycles, if any.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCtrlComputeSrfNCycle						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCtrlComputeCrvNCycle                                                 M
*****************************************************************************/
MvarPtStruct *MvarCtrlComputeCrvNCycle(const CagdCrvStruct *Crv,
				       int CycLen,
				       CagdRType SubdivTol,
				       CagdRType NumerTol)
{
    int i;
    CagdCrvStruct *TCrv;
    MvarMVStruct **MVCnsts, **MVCrvs, *MVCrv;
    MvarPtStruct *Pts, *Pt,
        *FilteredPts = NULL;

    if (CAGD_IS_RATIONAL_CRV(Crv)) {
        MVAR_FATAL_ERROR(MVAR_ERR_RATIONAL_NO_SUPPORT);
	return NULL;
    }

    MVCnsts = (MvarMVStruct **)
                            IritMalloc(sizeof(MvarMVStruct *) * (CycLen + 1));
    MVCrvs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * CycLen);

    TCrv = CagdCoerceCrvTo(Crv, CAGD_PT_E2_TYPE, FALSE);
    MVCrv = MvarCrvToMV(TCrv);
    CagdCrvFree(TCrv);

    /* Create the CycLen different curves: */
    for (i = 0; i < CycLen; i++)
        MVCrvs[i] = MvarPromoteMVToMV2(MVCrv, CycLen, i);

    /* Now build the CycLen chained constraints of "Ci(y) = Ci+1(x)". */
    for (i = 0; i < CycLen; i++) {
        MvarMVStruct *MVCrv1X, *MVCrv1Y, *MVCrv2X, *MVCrv2Y, **MVSplit;
	
	MVSplit = MvarMVSplitScalar(MVCrvs[i]);
	MVCrv1X = MVSplit[1];
	MVCrv1Y = MVSplit[2];

	MVSplit = MvarMVSplitScalar(MVCrvs[(i + 1) % CycLen]);
	MVCrv2X = MVSplit[1];
	MVCrv2Y = MVSplit[2];

	/* Build a constraints of MVCrv1Y = MVCrv2X. */
        MVCnsts[i] = MvarMVSub(MVCrv1Y, MVCrv2X);

	MvarMVFree(MVCrv1X);
	MvarMVFree(MVCrv1Y);
	MvarMVFree(MVCrv2X);
	MvarMVFree(MVCrv2Y);
    }

    /* Solve the constraints: */
    Pts = MvarMVsZeros0D(MVCnsts, NULL, CycLen, SubdivTol, NumerTol);

    for (i = 0; i < CycLen; i++) {
        MvarMVFree(MVCrvs[i]);
        MvarMVFree(MVCnsts[i]);
    }
    IritFree(MVCnsts);
    IritFree(MVCrvs);

    /* Filter out cyclic solutions:  if (t1, t2, t3, ..., tn) is a solution  */
    /* so is (t2, t3, ..., tn, t1).  Filter out one of this cyclic sets  by  */
    /* only selecting the solution were t1 is the minimal value.	     */
    while (Pts != NULL) {
	IRIT_LIST_POP(Pt, Pts);

	for (i = 1; i < Pt -> Dim; i++) {
	    if (Pt -> Pt[0] > Pt -> Pt[i] &&
		!IRIT_APX_EQ_EPS(Pt -> Pt[0], Pt -> Pt[i], 10.0 * NumerTol)) {
		break;
	    }
	}

	if (i >= Pt -> Dim) {
	    IRIT_LIST_PUSH(Pt, FilteredPts);
	}
	else {
	    MvarPtFree(Pt);
	}
    }

#   ifdef DEBUG
    /* Verify the solutions. */
    for (Pt = FilteredPts; Pt != NULL; Pt = Pt -> Pnext) {
	CagdRType
	    *R = CagdCrvEval(Crv, Pt -> Pt[0]),
	    PrevVal = R[2];

        for (i = 1; i <= Pt -> Dim; i++) {
	    CagdRType
		*R = CagdCrvEval(Crv, Pt -> Pt[i % Pt -> Dim]);

	    assert(IRIT_APX_EQ(R[1], PrevVal));
	    PrevVal = R[2];
	}
    }
#   endif /* DEBUG */

    return FilteredPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes a cycle of length CycLen of bouncing lines between the given    M
* surface S(u, v) = (S1(u, v), S2(u, v)).				     M
*  The cycle of length k = CycLen is:					     M
*  S(u1, v1) -> (u2, v)							     V
*  S(ui, vi) -> (u(i+1), v(i+1))					     V
*  S(uk, vk) -> (u1, v1)						     V
* Over all, we have 2n equations and 2n dofs.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:          A surface to compute a cycle of length CycLen for.         M
*   CycLen:       Length of sought cycle (also k above).	             M
*   SubdivTol:    Tolerance of the solution.  This tolerance is measured in  M
*		  the parametric space of the surface.			     M
*   NumerTol:     Numeric tolerance of a possible numeric improvment stage.  M
*		  The numeric stage is employed if NumericTol < SubdivTol.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:  Linked list of solutions, each holding the parameter    M
*		values of the different Cycles, if any.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCtrlComputeCrvNCycle						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCtrlComputeSrfNCycle                                                 M
*****************************************************************************/
MvarPtStruct *MvarCtrlComputeSrfNCycle(const CagdSrfStruct *Srf,
				       int CycLen,
				       CagdRType SubdivTol,
				       CagdRType NumerTol)
{
    int i, j;
    CagdRType UMin, UMax, VMin, VMax;
    CagdSrfStruct *TSrf, *TSrf2;
    MvarMVStruct **MVCnsts, **MVSrfs, *MVSrf;
    MvarPtStruct *Pts, *Pt,
        *FilteredPts = NULL;

    if (CAGD_IS_RATIONAL_SRF(Srf)) {
        MVAR_FATAL_ERROR(MVAR_ERR_RATIONAL_NO_SUPPORT);
	return NULL;
    }

    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

    MVCnsts = (MvarMVStruct **)
                         IritMalloc(sizeof(MvarMVStruct *) * 2 * (CycLen + 1));
    MVSrfs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * CycLen);

    /* Place the parametrization itself in the last two dims, in TSrf. */
    TSrf = CagdCoerceSrfTo(Srf, CAGD_PT_E1_TYPE, FALSE);
    TSrf2 = CagdCoerceSrfTo(TSrf, CAGD_PT_E3_TYPE, TRUE);
    CagdSrfFree(TSrf);
    TSrf = CagdCoerceSrfTo(Srf, CAGD_PT_E4_TYPE, FALSE);
    CAGD_GEN_COPY(TSrf -> Points[3], TSrf2 -> Points[2],
		  sizeof(CagdRType) * TSrf -> ULength * TSrf -> VLength);
    CAGD_GEN_COPY(TSrf -> Points[4], TSrf2 -> Points[3],
		  sizeof(CagdRType) * TSrf -> ULength * TSrf -> VLength);
    CagdSrfFree(TSrf2);

    MVSrf = MvarSrfToMV(TSrf);
    CagdSrfFree(TSrf);

    /* Create the CycLen different surfaces, each with two parameters: */
    for (i = 0; i < CycLen; i++) {
        MVSrfs[i] = MvarPromoteMVToMV2(MVSrf, CycLen * 2, i * 2);

	/* Set all domains to the proper range. */
	for (j = 0; j < CycLen * 2; j += 2) {
	    MvarMVSetDomain(MVSrfs[i], UMin, UMax, j, TRUE);
	    MvarMVSetDomain(MVSrfs[i], VMin, VMax, j + 1, TRUE);
	}
    }
    MvarMVFree(MVSrf);

    /* Build CycLen chained constraints as "S(ui, vi) = (u(i+1), v(i+1))": */
    for (i = 0; i < CycLen; i++) {
        MvarMVStruct *MVSrf1X, *MVSrf1Y, *MVSrf1U, *MVSrf1V,
                     *MVSrf2X, *MVSrf2Y, *MVSrf2U, *MVSrf2V, **MVSplit;
	
	MVSplit = MvarMVSplitScalar(MVSrfs[i]);
	MVSrf1X = MVSplit[1];
	MVSrf1Y = MVSplit[2];
	MVSrf1U = MVSplit[3];
	MVSrf1V = MVSplit[4];

	MVSplit = MvarMVSplitScalar(MVSrfs[(i + 1) % CycLen]);
	MVSrf2X = MVSplit[1];
	MVSrf2Y = MVSplit[2];
	MVSrf2U = MVSplit[3];
	MVSrf2V = MVSplit[4];

	/* Build a constraints of MVSrf1X = MVSrf2U and MVSrf1Y = MVSrf2V. */
        MVCnsts[i * 2]     = MvarMVSub(MVSrf1X, MVSrf2U);
        MVCnsts[i * 2 + 1] = MvarMVSub(MVSrf1Y, MVSrf2V);

	MvarMVFree(MVSrf1X);
	MvarMVFree(MVSrf1Y);
	MvarMVFree(MVSrf1U);
	MvarMVFree(MVSrf1V);
	MvarMVFree(MVSrf2X);
	MvarMVFree(MVSrf2Y);
	MvarMVFree(MVSrf2U);
	MvarMVFree(MVSrf2V);
    }

    /* Solve the constraints: */
    Pts = MvarMVsZeros0D(MVCnsts, NULL, CycLen * 2, SubdivTol, NumerTol);

    for (i = 0; i < CycLen; i++) {
        MvarMVFree(MVSrfs[i]);
        MvarMVFree(MVCnsts[i * 2]);
        MvarMVFree(MVCnsts[i * 2 + 1]);
    }
    IritFree(MVCnsts);
    IritFree(MVSrfs);

    /* Filter out cyclic solutions:  if ((u1, v1), (u2, v2),...,(un, vn)) is */
    /* solution so is ((u2, v2),...,(un, vn), (u1, v1)).                     */
    /*   Filter out one of this cyclic sets  by only selecting the solution  */
    /* were u1 is the minimal value, for all u's.			     */
    while (Pts != NULL) {
	IRIT_LIST_POP(Pt, Pts);

	for (i = 1; i < Pt -> Dim / 2; i++) {
	    if (Pt -> Pt[0] > Pt -> Pt[i * 2] &&
		!IRIT_APX_EQ_EPS(Pt -> Pt[0], Pt -> Pt[i * 2],
				 10.0 * NumerTol)) {
		break;
	    }
	}

	if (i >= Pt -> Dim / 2) {
	    IRIT_LIST_PUSH(Pt, FilteredPts);
	}
	else {
	    MvarPtFree(Pt);
	}
    }

#   ifdef DEBUG
    /* Verify the solutions. */
    for (Pt = FilteredPts; Pt != NULL; Pt = Pt -> Pnext) {
        for (i = 0; i < Pt -> Dim - 1; i += 2) {
	    CagdRType
		*R = CagdSrfEval(Srf, Pt -> Pt[i], Pt -> Pt[i + 1]);

	    assert(IRIT_APX_EQ(R[1], Pt -> Pt[(i + 2) % Pt -> Dim]) &&
		   IRIT_APX_EQ(R[2], Pt -> Pt[(i + 3) % Pt -> Dim]));
	}
    }
#   endif /* DEBUG */

    return FilteredPts;
}
