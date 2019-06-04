/******************************************************************************
* FlankMil.c - flank milling analysis of freeform surfaces.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Jul. 14.					      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/symb_lib.h"
#include "mvar_loc.h"

// #define MVAR_EXAMINE_GLBL_SRF_FLANK_MIL

static CagdCrvStruct *MvarFlankMillEvalCrvToEuclidean(
						    const CagdCrvStruct *Crv,
						    const CagdSrfStruct *Srf,
						    int CrvSizeReduction,
						    CagdBType ReparamU);
static CagdCrvStruct *MvarFlankMillProcessOneSolution(
						  const CagdCrvStruct *Sol,
						  const CagdSrfStruct *Srf,
						  CagdCrvStruct **BndryUVCrvs);
static CagdCrvStruct *MvarFlankMillLineAnalyzeOnePath(
					      const CagdSrfStruct *Srf, 
					      const CagdCrvStruct *StripBndry,
					      const CagdCrvStruct *OrientField,
					      CagdRType Tolerance,
					      CagdBType DoNextBndry,
					      CagdRType SubdivTol,
					      CagdRType NumeircTol);
#ifdef MVAR_EXAMINE_GLBL_SRF_FLANK_MIL
static CagdCrvStruct *MvarFlankMillLineAnalyzeOnePath2(
					      const CagdSrfStruct *Srf, 
					      const CagdCrvStruct *StripBndry,
					      CagdRType Tolerance,
					      CagdRType SubdivTol,
					      CagdRType NumericTol);
#endif /* MVAR_EXAMINE_GLBL_SRF_FLANK_MIL */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes strips that can be gouging-free flanked milled so that Srf is   M
* completely covered within Tolerance.  That is, every location on Srf is    M
* within tolerance to at least one tool position.			     M
*   Computation is done assuming the tool is of radius zero.  I.e. a line.   M
*   This, as surface can be offseted by tool radius, reducing tool to line.  M
*   Also computation is perfromed from VMin boundary of Srf.                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:              To compute flank milling tool path for.                M
*   Tolerance:        Maximal allowed distance between Srf and tool.         M
*   StripBoundriesUV: List of UV curves in the parametric domain of Srf      M
*                     that delineates the boundaries of the strips.          M
*   CrvSizeReduction: A reduction in size of traced curve while ensuring     M
*                     the Tolerance conservatively.			     M
*   SubdivTol:        Tolerance of the subdivision process.  Tolerance is    M
*		      measured in the parametric space of the multivariates. M
*   NumericTol:       Numeric tolerance of the numeric stage.  The numeric   M
*		      stage is employed only if NumericTol < SubdivTol.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *: List of UV curves in the parametric domain of Srf to    M
*                    move the tool along and cover (machine) Srf to within   M
*                    Tolerance.                                              M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarFlankMillLineAnalyze                                                 M
*****************************************************************************/
CagdCrvStruct *MvarFlankMillLineAnalyze(const CagdSrfStruct *Srf,
					CagdRType Tolerance,
					CagdCrvStruct **StripBoundriesUV,
					int CrvSizeReduction,
					CagdRType SubdivTol,
					CagdRType NumericTol)
{
    CagdRType UMin, UMax, VMin, VMax;
    CagdUVType UV1, UV2;
    CagdCrvStruct *StripBndryE3,
        *OrientField = NULL,
        *ToolPathsUV = NULL;

    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
    UV1[0] = UMin;
    UV1[1] = VMin;
    UV2[0] = UMax;
    UV2[1] = VMin;
    *StripBoundriesUV = CagdMergeUvUv(UV1, UV2);    /* First VMin boundary. */
    StripBndryE3 = CagdCrvFromMesh(Srf, 0, CAGD_CONST_V_DIR);

    do {
        CagdCrvStruct *Sol, *StripToolPathUV, *NextStripBndryUV,
	    *NextStripBndryE3,
	    *Sol2 = NULL;

#	ifdef MVAR_EXAMINE_GLBL_SRF_FLANK_MIL
	if (OrientField == NULL) {
	    Sol2 = MvarFlankMillLineAnalyzeOnePath2(Srf, StripBndryE3,
						    Tolerance,
						    SubdivTol, NumericTol);

	    CagdDbg(Sol2);
	    CagdCrvFreeList(Sol2);
	    return NULL;
	}
#	endif /* MVAR_EXAMINE_GLBL_SRF_FLANK_MIL */

        Sol = MvarFlankMillLineAnalyzeOnePath(Srf, StripBndryE3,
					      OrientField, Tolerance,
					      TRUE, SubdivTol, NumericTol);

	if (Sol == NULL)
	    break;
	else if ((!IRIT_APX_EQ(Sol -> Points[1][0], UMin) &&
		  !IRIT_APX_EQ(Sol -> Points[1][0], UMax)) ||
		 (!IRIT_APX_EQ(Sol -> Points[1][Sol -> Length - 1], UMin) &&
		  !IRIT_APX_EQ(Sol -> Points[1][Sol -> Length - 1], UMax))) {
	    /* Solution does not span the entire U domain - try to solve    */
	    /* just for the middle of the strip tool path.		    */
	    Sol2 = MvarFlankMillLineAnalyzeOnePath(Srf, StripBndryE3,
						   OrientField, Tolerance,
						   FALSE,
						   SubdivTol, NumericTol);
	}

	/*   Split the solution into a tool path middle strip curve and the */
	/* next strip boundary curve.					    */
	StripToolPathUV = MvarFlankMillProcessOneSolution(Sol, Srf,
							  &NextStripBndryUV);
	if (Sol2 != NULL) {
	    CagdCrvStruct *StripToolPathUV2, *NextStripBndryUV2;

	    StripToolPathUV2 = MvarFlankMillProcessOneSolution(Sol2, Srf,
							   &NextStripBndryUV2);
	    assert(NextStripBndryUV2 == NULL);
	    CagdCrvFreeList(StripToolPathUV);
	    StripToolPathUV = StripToolPathUV2;
	    CagdCrvFreeList(Sol2);
	}

	CagdCrvFreeList(Sol);

	/* Always converted into a single curve, even if more than one: */
	NextStripBndryE3 = MvarFlankMillEvalCrvToEuclidean(NextStripBndryUV,
							   Srf,
							   CrvSizeReduction,
							   TRUE);
	if (OrientField != NULL)
	    CagdCrvFree(OrientField);
	if (NextStripBndryE3 == NULL) {
	    CagdCrvFree(StripBndryE3);
	    CagdCrvFreeList(ToolPathsUV);
	    CagdCrvFreeList(*StripBoundriesUV);
	    *StripBoundriesUV = NULL;
	    return NULL;
	}

	OrientField = SymbCrvSub(NextStripBndryE3, StripBndryE3);

	CagdCrvFree(StripBndryE3);
	StripBndryE3 = NextStripBndryE3;

	*StripBoundriesUV = CagdListAppend(*StripBoundriesUV,
					   NextStripBndryUV);
	ToolPathsUV = CagdListAppend(StripToolPathUV, ToolPathsUV);
    }
    while (TRUE);

    if (OrientField != NULL)
        CagdCrvFree(OrientField);

    ToolPathsUV = CagdListReverse(ToolPathsUV);
    *StripBoundriesUV = CagdListReverse(*StripBoundriesUV);

    return ToolPathsUV;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Evaluate the given parametric curve in Srf domain to Euclidean space.    *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv:   In Srf domain to evaluate into Euclidean space.                   *
*   Srf:   To use in the mapping to Euclidean space.                         *
*   CrvSizeReduction: A reduction in size of traced curve while ensuring     *
*                     the Tolerance conservatively.			     *
*   ReparamU:  TRUE to also reparameterize the curve to follow Srf U vals.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:  Euclidean space curve.                                 *
*****************************************************************************/
static CagdCrvStruct *MvarFlankMillEvalCrvToEuclidean(
						     const CagdCrvStruct *Crv,
						     const CagdSrfStruct *Srf,
						     int CrvSizeReduction,
						     CagdBType ReparamU)
{
    int i;
    CagdRType UMin, VMin, UMax, VMax;
    CagdCrvStruct *EucCrv,
        *ChainedCrv = NULL,
        *ReducedCrv = NULL;

    if (Crv == NULL || Crv -> Length < 2)
        return NULL;

    if (Crv != NULL && Crv -> Pnext != NULL) {
        CagdCrvStruct *NextCrv, *TCrv,
	    *CrvList = CagdCrvCopyList(Crv);

	/* Chain the curves into one, in case more than one curve is given. */
	while (TRUE) {
	    for (TCrv = CrvList, NextCrv = NULL;
		 TCrv != NULL;
		 TCrv = TCrv -> Pnext) {
	        if (AttrGetIntAttrib(TCrv -> Attr, "_used") != TRUE &&
		    (NextCrv == NULL ||
		     TCrv -> Points[1][0] < NextCrv -> Points[1][0]))
		    NextCrv = TCrv;
	    }

	    if (NextCrv == NULL)
	        break; /* Done - no more curves. */
	    else {
	        AttrSetIntAttrib(&NextCrv -> Attr, "_used", TRUE);
		
		if (ChainedCrv == NULL)
		    ChainedCrv = CagdCrvCopy(NextCrv);
		else {
		    TCrv = CagdMergeCrvCrv(ChainedCrv, NextCrv, TRUE);
		    CagdCrvFree(ChainedCrv);
		    ChainedCrv = TCrv;
		}
	    }
	}
	CagdCrvFreeList(CrvList);
	Crv = ChainedCrv;
    }

    if (Crv -> Length > CrvSizeReduction) {
        /* Reduce the curve to CrvSizeReduction. */
	int i;
        CagdRType Err;

        ReducedCrv = BspCrvFitLstSqr(Crv, 3, CrvSizeReduction,
				     FALSE, CAGD_UNIFORM_PARAM, TRUE,
				     TRUE, &Err);

	/* Shift control points toward UMin by Err to be conservative. */
	Err = IRIT_ABS(Err);
	for (i = 0; i < ReducedCrv -> Length; i++) {
	    ReducedCrv -> Points[2][i] -= Err;
	}
	Crv = ReducedCrv;
    }

    EucCrv =  BspCrvNew(Crv -> Length, Crv -> Order, CAGD_PT_E3_TYPE);
    BspKnotUniformOpen(Crv -> Length, Crv -> Order, EucCrv -> KnotVector);

    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

    for (i = 0; i < Crv -> Length; i++) {
	CagdRType *R;
	CagdPType Pt;

	if (Crv -> Points[1][i] < UMin || Crv -> Points[1][i] > UMax ||
	    Crv -> Points[2][i] < VMin || Crv -> Points[2][i] > VMax) {
	    /* Went outside the surface - probably least squares fit failed. */
	    if (ChainedCrv != NULL)
	        CagdCrvFree(ChainedCrv);
	    if (ReducedCrv != NULL)
	        CagdCrvFree(ReducedCrv);
	    CagdCrvFree(EucCrv);
	    return NULL;
	}
	else
	    R = CagdSrfEval(Srf, Crv -> Points[1][i], Crv -> Points[2][i]);

	CagdCoerceToE3(Pt, &R, -1, Srf -> PType);
	EucCrv -> Points[1][i] = Pt[0];
	EucCrv -> Points[2][i] = Pt[1];
	EucCrv -> Points[3][i] = Pt[2];
    }

    if (ReparamU) {
	int Len = Crv -> Length;
	CagdRType UMin, UMax, VMin, VMax;

	CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

        /* Fix the knots to follow that of the U parameter values. */
	if (Crv -> Points[1][0] > Crv -> Points[1][Len - 1]) {
	    EucCrv -> KnotVector[0] =
		EucCrv -> KnotVector[1] = Crv -> Points[1][Len - 1];
	    for (i = 1; i < Crv -> Length; i++) {
		EucCrv -> KnotVector[i + 1] =
		    IRIT_MAX(EucCrv -> KnotVector[i] +
						(i > 1 ? IRIT_EPS * 10 : 0.0),
			     Crv -> Points[1][Len - i]);
	    }
	    EucCrv -> KnotVector[i + 1] = 
		EucCrv -> KnotVector[i + 2] = 
		    IRIT_MAX(EucCrv -> KnotVector[i], Crv -> Points[1][1]);
	}
	else {
            EucCrv -> KnotVector[0] = 
		EucCrv -> KnotVector[1] = Crv -> Points[1][0];
	    for (i = 1; i < Crv -> Length; i++) {
	        EucCrv -> KnotVector[i + 1] =
		    IRIT_MAX(EucCrv -> KnotVector[i] +
					        (i > 1 ? IRIT_EPS * 10 : 0.0),
			     Crv -> Points[1][i - 1]);
	    }
	    EucCrv -> KnotVector[i + 1] = 
		EucCrv -> KnotVector[i + 2] =
		    IRIT_MAX(EucCrv -> KnotVector[i], Crv -> Points[1][i - 2]);
	}
	CagdCrvSetDomain(EucCrv, UMin, UMax);
	BspKnotVerifyKVValidity(EucCrv -> KnotVector, EucCrv -> Order,
				EucCrv -> Length, IRIT_EPS);
    }

    if (ChainedCrv != NULL)
        CagdCrvFree(ChainedCrv);
    if (ReducedCrv != NULL)
        CagdCrvFree(ReducedCrv);

    return EucCrv;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Splits the given solution curve in (u, v, t) space into two different    *
* (u, v) and (u, t) curves.                                                  *
*   Assumes a single solution than spans all of u domain (as a future        *
* possibility in case the surface is not convex, a lower envelope            *
* computation can be performed here over all input curves for the closest    *
* lower envelope to the previous strip boundary.                             *
*   Note this single solution can still consist of more than one curve.      *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol:         The solution curve in (u, v, t) space.                      *
*   Srf:         The surface for which all this computation take place.      *
*   BndryUVCrvs: The (u, t) split curve.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:   The (u, v) tool path strip middle curve.              *
*****************************************************************************/
static CagdCrvStruct *MvarFlankMillProcessOneSolution(
						  const CagdCrvStruct *SolList,
						  const CagdSrfStruct *Srf,
						  CagdCrvStruct **BndryUVCrvs)
{
    CagdCrvStruct
        *ToolPathCrvs = NULL;

    *BndryUVCrvs = NULL;

    for ( ; SolList != NULL; SolList = SolList -> Pnext) {
        const CagdCrvStruct
	    *Sol = SolList;
        CagdCrvStruct *BndryUV,
	    *SolRev = NULL,
	    *TP = BspCrvNew(Sol -> Length, Sol -> Order, CAGD_PT_E2_TYPE);

	if (Sol -> PType == CAGD_PT_E3_TYPE)
	    BndryUV = BspCrvNew(Sol -> Length, Sol -> Order, CAGD_PT_E2_TYPE);
	else
	    BndryUV = NULL;

	if (Sol -> Points[1][0] > Sol -> Points[1][Sol -> Length - 1]) {
	    /* Solution is from UMax to UMin - reverse it. */
	    Sol = SolRev = CagdCrvReverse(Sol);
	}

	BspKnotUniformOpen(TP -> Length, TP -> Order, TP -> KnotVector);
	IRIT_GEN_COPY(TP -> Points[1], Sol -> Points[1],
		      sizeof(CagdRType) * Sol -> Length);
	IRIT_GEN_COPY(TP -> Points[2], Sol -> Points[2],
		      sizeof(CagdRType) * Sol -> Length);
	IRIT_LIST_PUSH(TP, ToolPathCrvs);

	if (BndryUV != NULL) {
	    BspKnotUniformOpen(TP -> Length, TP -> Order,
			       BndryUV -> KnotVector);
	    IRIT_GEN_COPY(BndryUV -> Points[1], Sol -> Points[1],
			  sizeof(CagdRType) * Sol -> Length);
	    IRIT_GEN_COPY(BndryUV -> Points[2], Sol -> Points[3],
			  sizeof(CagdRType) * Sol -> Length);
	    IRIT_LIST_PUSH(BndryUV, *BndryUVCrvs);
	}

	if (SolRev != NULL)
	    CagdCrvFree(SolRev);
    }

    return ToolPathCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the next tool path, given current (previous) strip boundary by  *
* deriving the tool control location (tool path).                            *
*   Uses the fact that StripBndry and NextStripBoundary should be at a       *
* Tolerance distance from the line and in between these two curves, the tool *
* is in contact with Srf along some curve in S(u, v). Here, given StripBndry,*
*   This variant assumes the flank milling line is along the U parametric    *
* direction throughout.                                                      *
*                                                                            *
* PARAMETERS:                                                                *
*   Srf:          To compute flank milling tool path for.                    *
*   StripBndry:   End of previous strip.  A curve parameterized along the U  *
*                 direction of Srf.     				     *
*   OrientField:  If exists, prescribes where to look for solutions of the   *
*                 tool path (u, v) locations.				     *
*   Tolerance:    Maximal allowed distance between Srf and tool.             *
*   DoNextBndry:  TRUE to compute the next boundary of the strip as well.    *
*                 FALSE will only compute the middle of the strip tool path. *
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        *
*		  measured in the parametric space of the multivariates.     *
*   NumericTol:   Numeric tolerance of the numeric stage.  The numeric       *
*		  stage is employed only if NumericTol < SubdivTol.          *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct:   The next strip boundary(ies).			     *
*****************************************************************************/
static CagdCrvStruct *MvarFlankMillLineAnalyzeOnePath(
					      const CagdSrfStruct *Srf, 
					      const CagdCrvStruct *StripBndry,
					      const CagdCrvStruct *OrientField,
					      CagdRType Tolerance,
					      CagdBType DoNextBndry,
					      CagdRType SubdivTol,
					      CagdRType NumericTol)
{
    IRIT_STATIC_DATA MvarConstraintType
	Constraints[4] = { MVAR_CNSTRNT_ZERO,
			   MVAR_CNSTRNT_ZERO,
			   MVAR_CNSTRNT_POSITIVE,
			   MVAR_CNSTRNT_POSITIVE };
    int n,
        MultInterpFlag = BspMultComputationMethod(BSP_MULT_BEZ_DECOMP);
    CagdRType UMin, UMax, VMin, VMax;
    CagdSrfStruct
        *NSrf = SymbSrfNormalSrf(Srf);
    MvarMVStruct *MVTemp1, *MVTemp2, *MVTemp3, *MVTemp4, *MVSrf, *MVSrf2,
      *MVNSrf, *MVCrv, *MVOrient, *MVs[4];
    MvarPolylineStruct *TPath;
    CagdCrvStruct *TPathCrvs;

    IRIT_ZAP_MEM(MVs, sizeof(MvarMVStruct *) * 4);

    if (DoNextBndry) {
        /* Convert Srf(u, v) and Crv(u) to MVs(u, v, t).  Srf(u, t) is     */
        /* derived from srf(u, v) by reversing (swapping) the 2nd and 3rd  */
        /* parameter.						           */
        MVTemp1 = MvarSrfToMV(Srf);
	MVSrf = MvarPromoteMVToMV2(MVTemp1, 3, 0);
	MvarMVFree(MVTemp1);
	MVSrf2 = MvarMVReverse(MVSrf, 1, 2);

	MVTemp1 = MvarSrfToMV(NSrf);
	MVNSrf = MvarPromoteMVToMV2(MVTemp1, 3, 0);
	MvarMVFree(MVTemp1);
	CagdSrfFree(NSrf);

	MVTemp1 = MvarCrvToMV(StripBndry);
	MVCrv = MvarPromoteMVToMV2(MVTemp1, 3, 0);
	MvarMVFree(MVTemp1);

	/* Update mutual domains to be the same. */
	CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
	MvarMVSetDomain(MVCrv, VMin, VMax, 1, TRUE);
	MvarMVSetDomain(MVCrv, VMin, VMax, 2, TRUE);
	MvarMVSetDomain(MVSrf, VMin, VMax, 2, TRUE);
	MvarMVSetDomain(MVNSrf, VMin, VMax, 2, TRUE);
	MvarMVSetDomain(MVSrf2, VMin, VMax, 1, TRUE);

	if (OrientField != NULL) {
	    MVTemp1 = MvarCrvToMV(OrientField);
	    MVOrient = MvarPromoteMVToMV2(MVTemp1, 3, 0);
	    MvarMVFree(MVTemp1);
	    MvarMVSetDomain(MVOrient, VMin, VMax, 1, TRUE);
	    MvarMVSetDomain(MVOrient, VMin, VMax, 2, TRUE);
	}
	else
	    MVOrient = NULL;
    }
    else {
        /* Convert Srf(u, v) and Crv(u) to MVs(u, v). */
        MVSrf = MvarSrfToMV(Srf);
	MVNSrf = MvarSrfToMV(NSrf);
	CagdSrfFree(NSrf);

	MVTemp1 = MvarCrvToMV(StripBndry);
	MVCrv = MvarPromoteMVToMV2(MVTemp1, 2, 0);
	MvarMVFree(MVTemp1);

	/* Update mutual domains to be the same. */
	CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
	MvarMVSetDomain(MVCrv, VMin, VMax, 1, TRUE);

	if (OrientField != NULL) {
	    MVTemp1 = MvarCrvToMV(OrientField);
	    MVOrient = MvarPromoteMVToMV2(MVTemp1, 2, 0);
	    MvarMVFree(MVTemp1);
	    MvarMVSetDomain(MVOrient, VMin, VMax, 1, TRUE);
	}
	else
	    MVOrient = NULL;
    }

    /* Build the equation (parameters are ordered (u, v) or (u, v, t)):   */
    /*  < S(u, v) - C(u), n(u, v) >^2 = Tolerance^2 * || n(u, v) ||^2.    */
    /*  < S(u, t) - S(u, v), n(u, v) >^2 = Tolerance^2 * || n(u, v) ||^2. */
    /*  < S(u, t) - S(u, v),  S(u, v) - C(u) > > 0.                       */
    /* and if OrientField exists also add:                                */
    /*  < S(u, v) - C(u), OrientField > > 0.		                  */
    /* Note equations with t are only created if DoNextBndry is TRUE.     */
    MVTemp1 = MvarMVSub(MVSrf, MVCrv);
    MVTemp2 = MvarMVDotProd(MVTemp1, MVNSrf);
    MVTemp3 = MvarMVMult(MVTemp2, MVTemp2);
    MvarMVFree(MVTemp2);

    if (DoNextBndry) {
        MVTemp2 = MvarMVSub(MVSrf2, MVSrf);
	MVs[2] = MvarMVDotProd(MVTemp1, MVTemp2);
	MVs[3] = MVOrient != NULL ? MvarMVDotProd(MVTemp1, MVOrient) : NULL;
    }
    else
        MVs[1] = MVOrient != NULL ? MvarMVDotProd(MVTemp1, MVOrient) : NULL;

    MvarMVFree(MVTemp1);

    if (DoNextBndry) {
        MVTemp1 = MvarMVDotProd(MVTemp2, MVNSrf);
	MVTemp4 = MvarMVMult(MVTemp1, MVTemp1);
	MvarMVFree(MVTemp1);
	MvarMVFree(MVTemp2);
    }

    MVTemp1 = MvarMVDotProd(MVNSrf, MVNSrf);
    MVTemp2 = MvarMVScalarScale(MVTemp1, IRIT_SQR(Tolerance));
    MvarMVFree(MVTemp1);

    MVs[0] = MvarMVSub(MVTemp3, MVTemp2);
    if (DoNextBndry) {
        MVs[1] = MvarMVSub(MVTemp4, MVTemp2);
	MvarMVFree(MVTemp4);
	MvarMVFree(MVSrf2);
    }
    MvarMVFree(MVTemp2);
    MvarMVFree(MVTemp3);

    MvarMVFree(MVCrv);
    MvarMVFree(MVSrf);
    MvarMVFree(MVNSrf);
    if (MVOrient != NULL)
        MvarMVFree(MVOrient);

    /* Count how many valid constraints we have. */
    for (n = 0; n < 4 && MVs[n] != NULL; n++);

    /* Time to solve the equations. */
    TPath = MvarMVsZeros1D(MVs, DoNextBndry ? Constraints : &Constraints[1],
			   n, SubdivTol * 0.3, SubdivTol, NumericTol);
    while (--n >= 0)
        MvarMVFree(MVs[n]);


    TPathCrvs = MvarCnvrtMVPolysToIritCrvs(TPath, 3);
    MvarPolylineFreeList(TPath);

    BspMultComputationMethod(MultInterpFlag);

    return TPathCrvs;
}

#ifdef MVAR_EXAMINE_GLBL_SRF_FLANK_MIL

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the next tool path, given current strip boundary C(u) by        *
* deriving the surface location of the tool-surface tangential contacts.     *
*   Search all of S(u, v) for the locus of points, p(u) in S, such that for  *
 each p(u), the tangent plane of S at p(u) is Tolerance away from C(u).      *
*   This variant assumes the flank milling line is along the U parametric    *
* direction throughout.                                                      *
*   CURRENTLY NOT USED.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   Srf:          To compute flank milling tool path for.                    *
*   StripBndry:   End of previous strip.				     *
*   Tolerance:    Maximal allowed distance between Srf and tool.             *
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        *
*		  measured in the parametric space of the multivariates.     *
*   NumericTol:   Numeric tolerance of the numeric stage.  The numeric       *
*		  stage is employed only if NumericTol < SubdivTol.          *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct:   The next strip boundary(ies).			     *
*****************************************************************************/
static CagdCrvStruct *MvarFlankMillLineAnalyzeOnePath2(
					      const CagdSrfStruct *Srf, 
					      const CagdCrvStruct *StripBndry,
					      CagdRType Tolerance,
					      CagdRType SubdivTol,
					      CagdRType NumericTol)
{
    CagdRType UMin, UMax, VMin, VMax, TMin, TMax;
    CagdSrfStruct
        *NSrf = SymbSrfNormalSrf(Srf);
    MvarMVStruct *MVTemp1, *MVTemp2, *MVTemp3, *MVSrf, *MVNSrf, *MVCrv, *MVs[2];
    MvarPolylineStruct *TPath;
    CagdCrvStruct *TPathCrvs;

    /* Convert Srf(u, v) and Crv(w) to Trivars(u, v, w). */
    MVTemp1 = MvarSrfToMV(Srf);
    MVSrf = MvarPromoteMVToMV2(MVTemp1, 3, 0);
    MvarMVFree(MVTemp1);

    MVTemp1 = MvarSrfToMV(NSrf);
    MVNSrf = MvarPromoteMVToMV2(MVTemp1, 3, 0);
    MvarMVFree(MVTemp1);
    CagdSrfFree(NSrf);

    MVTemp1 = MvarCrvToMV(StripBndry);
    MVCrv = MvarPromoteMVToMV2(MVTemp1, 3, 2);
    MvarMVFree(MVTemp1);

    /* Update mutual domains to be the same. */
    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
    CagdCrvDomain(StripBndry, &TMin, &TMax);
    MvarMVSetDomain(MVCrv, UMin, UMax, 0, TRUE);
    MvarMVSetDomain(MVCrv, VMin, VMax, 1, TRUE);
    MvarMVSetDomain(MVSrf, TMin, TMax, 2, TRUE);
    MvarMVSetDomain(MVNSrf, TMin, TMax, 2, TRUE);

    /* Build the equations:					 	    */
    /*  f(u, v, w) = 							    */
    /*   < S(u, v) - C(w), n(u, v) >^2 - Tolerance^2 * || n(u, v) ||^2 = 0. */
    /*  df(u, v, w) / dw = 0.                                               */
    /* Note second eqn computes envelope of f(u, v, w) with respect to w.   */
    MVTemp1 = MvarMVSub(MVSrf, MVCrv);
    MVTemp2 = MvarMVDotProd(MVTemp1, MVNSrf);
    MVTemp3 = MvarMVMult(MVTemp2, MVTemp2);
    MvarMVFree(MVTemp1);
    MvarMVFree(MVTemp2);

    MVTemp1 = MvarMVDotProd(MVNSrf, MVNSrf);
    MVTemp2 = MvarMVScalarScale(MVTemp1, IRIT_SQR(Tolerance));
    MvarMVFree(MVTemp1);

    MVs[0] = MvarMVSub(MVTemp3, MVTemp2);
    MvarMVFree(MVTemp2);
    MvarMVFree(MVTemp3);

    MVs[1] = MvarMVDerive(MVs[0], 2);

    /* Time to solve the equations. */
    TPath = MvarMVsZeros1D(MVs, NULL, 2,
			   SubdivTol * 0.3, SubdivTol, NumericTol);
    MvarMVFree(MVs[0]);
    MvarMVFree(MVs[1]);

    TPathCrvs = MvarCnvrtMVPolysToIritCrvs(TPath, 3);
    MvarPolylineFreeList(TPath);

    return TPathCrvs;
}

#endif /* MVAR_EXAMINE_GLBL_SRF_FLANK_MIL */
