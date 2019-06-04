/******************************************************************************
* UnTrim.c - computes trimmed surfaces decomposition into tensor products.    *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Nov. 2013.					      *
******************************************************************************/

#include <assert.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "symb_loc.h"

#define TRIM_MATCH_MAX_SIZE	20
#define TRIM_MATCH_REDUCED_SIZE 3
#define TRIM_MATCH_SAMPLES_SIZE 20
#define TRIM_MATCH_JACOBIAN_ZERO 1e-10

static CagdCrvStruct *SymbCleanUnstrictEndMonotone(const CagdCrvStruct *Crvs,
						   int Dir);
static CagdRType SymbMatchRulingNorm(const CagdVType T1,
				     const CagdVType T2,
				     const CagdVType P1,
				     const CagdVType P2);
static CagdRType Symb2DSrfParamJacobian(const CagdSrfStruct *Srf);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Clean end points in given curves, if piecewise linear, that are constant *
* with respect to monotone direction.					     *
*                                                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   Crvs:        Curves to filter-out, in place, non strict monotone end     *
*                points.					             *
*   Dir:         1 for U or 2 for V.                                         *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:  Filtered curves.                                       *
*****************************************************************************/
static CagdCrvStruct *SymbCleanUnstrictEndMonotone(const CagdCrvStruct *Crvs,
						   int Dir)
{
    const CagdCrvStruct *Crv;
    CagdCrvStruct *TCrv,
        *NewCrvs = NULL;

    for (Crv = Crvs; Crv != NULL; Crv = Crv -> Pnext) {
        int i1, i2, j, n;
        CagdRType *C1Disconts, Pt0E2[2], PtE2[2];

	TCrv = CagdCrvCopy(Crv);
	C1Disconts = BspKnotAllC1Discont(TCrv -> KnotVector, TCrv -> Order,
					 TCrv -> Length, &n);

	if (C1Disconts == NULL) {
	    IRIT_LIST_PUSH(TCrv, NewCrvs);
	    continue;
	}

	/* Remove first identical value points from curve, if any. */
	for (i1 = 0; TCrv -> Length > TCrv -> Order && i1 < n; i1++) {
	    CagdCrvStruct
	        *Crv1 = CagdCrvSubdivAtParam(TCrv, C1Disconts[i1]);

	    CagdCoerceToE2(Pt0E2, Crv1 -> Points, 0, Crv1 -> PType);
	    for (j = 1; j < Crv1 -> Length; j++) {
	        /* Examine the control points in Dir if identical. */
	        CagdCoerceToE2(PtE2, Crv1 -> Points, j, Crv1 -> PType);
		if (!IRIT_APX_EQ(PtE2[Dir - 1], Pt0E2[Dir - 1]))
		    break;
	    }
	    if (j >= Crv1 -> Length) { /* Purge - pts are indeed identical. */
	        CagdCrvFree(TCrv);
		TCrv = Crv1 -> Pnext;
		CagdCrvFree(Crv1);
	    }
	    else {
	        CagdCrvFreeList(Crv1);
		break;
	    }
	}

	/* Remove last identical value points from curve. */
	for (i2 = n - 1; TCrv -> Length > TCrv -> Order && i2 > i1; i2--) {
	    CagdCrvStruct
	        *Crv1 = CagdCrvSubdivAtParam(TCrv, C1Disconts[i2]),
	        *Crv2 = Crv1 -> Pnext;

	    CagdCoerceToE2(Pt0E2, Crv2 -> Points, 0, Crv2 -> PType);
	    for (j = 1; j < Crv2 -> Length; j++) {
	        /* Examine the control points in Dir if identical. */
	        CagdCoerceToE2(PtE2, Crv2 -> Points, j, Crv2 -> PType);
		if (!IRIT_APX_EQ(PtE2[Dir - 1], Pt0E2[Dir - 1]))
		    break;
	    }
	    if (j >= Crv2 -> Length) { /* Purge - pts are indeed identical. */
	        CagdCrvFree(TCrv);
		TCrv = Crv1;
		TCrv -> Pnext = NULL;
		CagdCrvFree(Crv2);
	    }
	    else {
	        CagdCrvFreeList(Crv1);
		break;
	    }
	}

	IRIT_LIST_PUSH(TCrv, NewCrvs);
    };

    return CagdListReverse(NewCrvs);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the ruled norm for the matching, ensuring positivity of J:	     *
* < T1 x (P2 - P1), T2 x (P2 - P1) > must be positive and then the norm      *
* is 1.0 - MIN( || T1 x (P2 - P1) ||^2, || T2 x (P2 - P1) ||^2 )             *
*                                                         / || P2 - P1 ||^2. *
*                                                                            *
* PARAMETERS:                                                                *
*   T1: A pointer to unit tangent to the first curve at i-th point.          *
*   T2: A pointer to unit tangent to the second curve at j-th point.         *
*   P1:	A pointer to value of the first curve at i-th point.                 *
*   P2: A pointer to value of the second curve at j-th point.                *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: -1 for no matching or the cost of the matching for the point  *
*		between zero and one, the minimum (non negative) the better. *
*****************************************************************************/
static CagdRType SymbMatchRulingNorm(const CagdVType T1,
				     const CagdVType T2,
				     const CagdVType P1,
				     const CagdVType P2)
{
    CagdPType RuledVec;
    CagdRType CrossVec1, CrossVec2, d1, d2;

    assert(IRIT_APX_EQ(IRIT_PT2D_LENGTH(T1), 1.0) &&
	   IRIT_APX_EQ(IRIT_PT2D_LENGTH(T2), 1.0));

    IRIT_PT2D_SUB(RuledVec, P1, P2);
    IRIT_PT2D_SAFE_NORMALIZE(RuledVec);
    CrossVec1 = IRIT_CROSS_PROD_2D(T1, RuledVec);
    CrossVec2 = IRIT_CROSS_PROD_2D(T2, RuledVec);

    if (IRIT_PT2D_SQR_LENGTH(RuledVec) < IRIT_SQR(IRIT_UEPS))
	return 0.0;
    else if (CrossVec1 * CrossVec2 < 0.0) /* Force validity of match. */
	return -1.0;

    d1 = IRIT_SQR(CrossVec1);
    d2 = IRIT_SQR(CrossVec2);
    return 1.0 - IRIT_MIN(d1, d2);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Parameterize the area enclosed by given pair of curves that form a       M
* closed planar loop.				                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1, Crv2:  Two curves forming a loop 9and sharing starting and ending  M
*                locations (i.e. they are roughly monotone with respect to   M
*                each other.)						     M
*   Dir:         Direction where Crv1/2 were split (1 for U, 2 for V).	     M
*   ForceMatch:  Force match, even if the Jacobian might be negative.        M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *: 2D surface parameterizing region enclosed in Crv1/2,    M
*                    or NULL if failed.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrimDecompTrimSrf2Srfs, TrimCrvIsParameterizableDomain                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   Symb2DCrvParameterizing2Crvs                                             M
*****************************************************************************/
CagdSrfStruct *Symb2DCrvParameterizing2Crvs(const CagdCrvStruct *Crv1,
					    const CagdCrvStruct *Crv2,
					    int Dir,
					    CagdBType ForceMatch)
{
    CagdCrvStruct
        *Crv1f = Dir == 0 ? CagdCrvCopy(Crv1)
                          : SymbCleanUnstrictEndMonotone(Crv1, Dir),
        *Crv2f = Dir == 0 ? CagdCrvCopy(Crv2)
                          : SymbCleanUnstrictEndMonotone(Crv2, Dir);
    CagdSrfStruct *Srf;
    int BaseSize = IRIT_MIN(IRIT_MAX(Crv1f -> Length, Crv2f -> Length),
			    TRIM_MATCH_MAX_SIZE);
    CagdCrvStruct
        *Crv2fn = CagdMatchingTwoCurves(Crv1f, Crv2f,
					BaseSize * TRIM_MATCH_REDUCED_SIZE,
					BaseSize * TRIM_MATCH_SAMPLES_SIZE,
					3, FALSE, FALSE, FALSE,
					TRUE, SymbMatchRulingNorm);
    CagdRType
	Error = Crv2fn == NULL ? IRIT_INFNTY
			       : AttrGetRealAttrib(Crv2fn -> Attr, "_Error");

    if (Crv2fn != NULL) { /* Successful match. */
        CagdCrvFree(Crv2f);
	Crv2f = Crv2fn;
    }
    else if (!ForceMatch) {
        CagdCrvFree(Crv2f);
	CagdCrvFree(Crv1f);
	return NULL;
    }

    Srf = CagdRuledSrf(Crv1f, Crv2f, 2, 2);
    AttrSetRealAttrib(&Srf -> Attr, "_Error", Error);

    CagdCrvFree(Crv2f);
    CagdCrvFree(Crv1f);

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Fit a bivariate to the region enclosed by Crv by splitting at T1 and T2. M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:      The curve to fit a bivariate to the area enclosed in.          M
*   T1, T2:   Two parameters to split curve at and fit a match in between.   M
*   Dir:      Direction where Crv1/2 were split (0 None, 1 for U, 2 for V).  M
*   Error:    Error result in this fit.  Error can only be compared to the   M
*             error of similar invocations of this function.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   Fitted bivariate, or NULL if failed.                  M
*                                                                            *
* SEE ALSO:                                                                  M
*                                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   Symb2DCrvParamerize2Prms                                                 M
*****************************************************************************/
CagdSrfStruct *Symb2DCrvParamerize2Prms(const CagdCrvStruct *Crv,
					CagdRType T1,
					CagdRType T2,
					int Dir,
					CagdRType *Error)
{
    CagdRType TMin, TMax, RuledJacobian;
    CagdCrvStruct *Crv1, *Crv2, *TCrv;
    CagdSrfStruct *UVRuledSrf;

    CagdCrvDomain(Crv, &TMin, &TMax);

    if (T1 < TMin ||
	T1 > TMax ||
	T2 < TMin ||
	T1 > TMax ||
	(IRIT_APX_UEQ(T1, TMin) && IRIT_APX_UEQ(T2, TMax)) ||
	(IRIT_APX_UEQ(T2, TMin) && IRIT_APX_UEQ(T1, TMax)) ||
	IRIT_APX_EQ(T1, T2))
        return NULL;

    if (T1 > T2)
        IRIT_SWAP(CagdRType, T1, T2);

    if (IRIT_APX_EQ(T1, TMin)) {
        assert(!IRIT_APX_EQ(T2, TMax));

        /* Extreme domain in Dir starts at curve beginning. */
        Crv1 = CagdCrvSubdivAtParam(Crv,IRIT_MAX(T1, T2));
	Crv2 = Crv1 -> Pnext;
	Crv1 -> Pnext = NULL;
    }
    else if (IRIT_APX_EQ(T2, TMax)) {
        assert(!IRIT_APX_EQ(T1, TMin));

	/* Extreme domain in Dir ends at curve end. */
	Crv1 = CagdCrvSubdivAtParam(Crv, IRIT_MIN(T1, T2));
	Crv2 = Crv1 -> Pnext;
	Crv1 -> Pnext = NULL;
    }
    else {
        /* Extreme domain in Dir is interior to curve domain. */
        Crv1 = CagdCrvSubdivAtParam(Crv, IRIT_MIN(T1, T2));
	TCrv = Crv1 -> Pnext;
	Crv1 -> Pnext = NULL;
	Crv2 = CagdCrvSubdivAtParam(TCrv, IRIT_MAX(T1, T2));
	CagdCrvFree(TCrv);
	TCrv = CagdMergeCrvCrv(Crv2 -> Pnext, Crv1, FALSE);
	CagdCrvFree(Crv2 -> Pnext);
	Crv2 -> Pnext = NULL;
	CagdCrvFree(Crv1);
	Crv1 = TCrv;
    }

    /* Reverse one curve so they are parallel to each other. */
    TCrv = CagdCrvReverse(Crv1);
    CagdCrvFree(Crv1);
    Crv1 = TCrv;

    /* Do simple ruling. */
    UVRuledSrf = CagdRuledSrf(Crv1, Crv2, 2, 2);
    RuledJacobian = Symb2DSrfParamJacobian(UVRuledSrf);

    /* Try to forum a valid matching between Crv1 and Crv2 for ruling. */
    if (Crv1 != NULL && Crv2 != NULL) {
        CagdSrfStruct
            *UVSrf = Symb2DCrvParameterizing2Crvs(Crv1, Crv2, Dir, FALSE);
	CagdRType
	    Jacobian = UVSrf != NULL ? Symb2DSrfParamJacobian(UVSrf) : -1.0;

	CagdCrvFree(Crv1);
	CagdCrvFree(Crv2);
	if (Jacobian <= RuledJacobian + IRIT_UEPS) {
	    if (UVSrf != NULL)
		CagdSrfFree(UVSrf);
	    UVSrf = UVRuledSrf;
	}
	else
	    CagdSrfFree(UVRuledSrf);

	if (UVSrf != NULL) {
	    *Error = AttrGetRealAttrib(UVSrf -> Attr, "_Error");
	    return UVSrf;
	}
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Estimates the quality of the Jacobian of the given Srf as a mapping from *
* R^2 to R^2.                                                                *
*                                                                            *
* PARAMETERS:                                                                *
*   Srf:     A mapping from R^2 to R^2, to estimate its Jacobian's quality.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:   Negative if mapping self-intersects.  Otherwise, the ratio  *
*                between maximum and minimum assumed Jacobian values, or     *
*                zero if the Jacobian has locations with zero values.        *
*****************************************************************************/
static CagdRType Symb2DSrfParamJacobian(const CagdSrfStruct *Srf)
{
    CagdBBoxStruct BBox;
    BspMultComputationMethodType
	OldVal = BspMultComputationMethod(BSP_MULT_BEZ_DECOMP);
    CagdSrfStruct
        *J = Symb2DSrfJacobian(Srf);

    BspMultComputationMethod(OldVal);

    CagdSrfBBox(J, &BBox);

    CagdSrfFree(J);

    if (IRIT_FABS(BBox.Min[0]) < TRIM_MATCH_JACOBIAN_ZERO ||
	IRIT_FABS(BBox.Max[0]) < TRIM_MATCH_JACOBIAN_ZERO ||
	BBox.Min[0] * BBox.Max[0] > TRIM_MATCH_JACOBIAN_ZERO) {
        /* No negative Jacobian. */
	if (IRIT_FABS(BBox.Min[0]) < TRIM_MATCH_JACOBIAN_ZERO ||
	    IRIT_FABS(BBox.Max[0]) < TRIM_MATCH_JACOBIAN_ZERO)
	    return 0.0;
	else
	    return BBox.Max[0] < 0.0 ? BBox.Min[0] / BBox.Max[0] 
                                     : BBox.Max[0] / BBox.Min[0];
    }
    else
        return -1.0;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Examine curve Crv if its enclosed domain is bivariate parameterizable.   M
*                                                                            *
* PARAMETERS:                                                                M
*   UVCrv:        Curve to examine the ability to parameterize its interior  M
*                 domain as a 2D surface.			             M
*   Eps:          Tolerance of the decomposition.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   Planars surface parameterizing the domain enclosed.   M
*                 Typically one, possibly two surfaces, or NULL if failed.   M
*                                                                            *
* KEYWORDS:                                                                  M
*   Symb2DCrvParameterizeDomain                                              M
*****************************************************************************/
CagdSrfStruct *Symb2DCrvParameterizeDomain(const CagdCrvStruct *UVCrv,
					   CagdRType Eps)
{
    int Dir;
    CagdRType TMin, TMax,
        BestJ = -1.0,
	MinError = IRIT_INFNTY;
    CagdCrvStruct *Crv;
    CagdSrfStruct *TSrf, *TSrfs,
        *BestUVSrfs = NULL;

    Crv = CagdCoerceCrvTo(UVCrv, CAGD_IS_RATIONAL_CRV(UVCrv) ? CAGD_PT_P2_TYPE
	                                                     : CAGD_PT_E2_TYPE,
			  FALSE);
    CagdCrvDomain(Crv, &TMin, &TMax);

    /* Try the extreme values in the X and then Y direction: */
    for (Dir = 1; Dir <= 2; Dir++) {
        CagdRType *R, PtE2[2], Error,
	    Min = IRIT_INFNTY,
	    Max = -IRIT_INFNTY,
	    MinT = IRIT_INFNTY,
	    MaxT = IRIT_INFNTY;
	CagdSrfStruct
	    *UVSrf = NULL;
	CagdPtStruct *Ex,
	    *Extremes = SymbCrvExtremSet(Crv, Dir, Eps, FALSE);

	/* Add TMin and TMax for the tested values of extremes. */
	Ex = CagdPtNew();
	Ex -> Pt[0] = TMin;
	IRIT_LIST_PUSH(Ex, Extremes);

	Ex = CagdPtNew();
	Ex -> Pt[0] = TMax;
	IRIT_LIST_PUSH(Ex, Extremes);

	for (Ex = Extremes; Ex != NULL; Ex = Ex -> Pnext) {
	    R = CagdCrvEval(Crv, Ex -> Pt[0]);
	    CagdCoerceToE2(PtE2, &R, -1, Crv -> PType);
	    if (Min > PtE2[Dir - 1]) {
	        Min = PtE2[Dir - 1];
		MinT = Ex -> Pt[0];
	    }
	    if (Max < PtE2[Dir - 1]) {
	        Max = PtE2[Dir - 1];
		MaxT = Ex -> Pt[0];
	    }
	}
	if (Extremes != NULL)
	    CagdPtFreeList(Extremes);

	UVSrf = Symb2DCrvParamerize2Prms(Crv, MinT, MaxT, Dir, &Error);
	if (UVSrf != NULL && MinError > Error) {
	    MinError = Error;
	    if (BestUVSrfs != NULL)
	        CagdSrfFree(BestUVSrfs);
	    BestUVSrfs = UVSrf;
	}
	else
	    CagdSrfFree(UVSrf);
    }

    if (BestUVSrfs != NULL)
        BestJ = Symb2DSrfParamJacobian(BestUVSrfs);

    /* Try to fit 4-sided patches to the domain. */
    if ((TSrfs = CagdSrfFromNBndryCrvs(Crv, TRUE)) != NULL) {
        for (TSrf = TSrfs; TSrf != NULL; TSrf = TSrf -> Pnext) {
	    CagdRType
	        EstJ = Symb2DSrfParamJacobian(TSrf);

	    if (EstJ < 0.0 || (BestJ > 0 && EstJ > 0 && EstJ > BestJ))
	        break;
	}

	if (TSrf == NULL) {
	    CagdSrfFree(BestUVSrfs);
	    BestUVSrfs = TSrfs;
	}
	else
	    CagdSrfFreeList(TSrfs);	    
    }

    CagdCrvFree(Crv);

    return BestUVSrfs;
}
