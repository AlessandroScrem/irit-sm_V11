/******************************************************************************
* SymbZero.c - computes the zeros and extremes of a given object.	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, March 93.					      *
******************************************************************************/

#include "symb_loc.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/mvar_lib.h"

#define SYMB_ZERO_DEFAULT_IRIT_EPS  1e-6
#define SYMB_ZERO_IRIT_APX_EQ(x, y) (IRIT_FABS((x) - (y)) < GlblSetEpsilon * 10)
#define SYMB_SUBDIV_TOL		    1e-6

IRIT_STATIC_DATA CagdPtStruct
    *GlblPtList = NULL;
IRIT_STATIC_DATA CagdRType
    GlblSetEpsilon = SYMB_ZERO_DEFAULT_IRIT_EPS;

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes the zero set of a given curve, in given axis (0/1-3 for W/X-Z).   M
*   Returned is a list of the zero set points holding the parameter values   M
* at Pt[0] of each point.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:        To compute its zero set.                                     M
*   Axis:       The axis of Crv to compute zero set for, W = 0, X = 1, etc.  M
*   Epsilon:    Tolerance control.                                           M
*   NoSolsOnEndPts: If TRUE, solutions at the end of the domain are purged.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:   List of parameter values form which Crv is zero in     M
*                     axis Axis.                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvZeroSet, zero set, symbolic computation                           M
*****************************************************************************/
CagdPtStruct *SymbCrvZeroSet(const CagdCrvStruct *Crv,
			     int Axis,
			     CagdRType Epsilon,
			     CagdBType NoSolsOnEndPts)
{
    CagdRType TMin, TMax;
    CagdPtStruct *Pt, *Pts;

    if (Crv -> Order < 2)
        return NULL;

    assert(Axis >= 0 && Axis <= 3);

    GlblPtList = NULL;
    CagdCrvDomain(Crv, &TMin, &TMax);

    if (Crv -> Order == 2) {
        /* It is a polyline - solve for the zeros right here. */
	int i;
        CagdRType t,
	    *R = Crv -> Points[Axis];

	for (i = 1; i < Crv -> Length; i++) {
	    if (R[i] * R[i - 1] <= 0.0 &&
		(R[i] != 0.0 || R[i-1] != 0.0)) {
	        t = R[i - 1] / (R[i - 1] - R[i]);
		if (CAGD_IS_BSPLINE_CRV(Crv)) {
		    CagdRType
			*KV = &Crv -> KnotVector[1];

		    t = KV[i] * t + KV[i - 1] * (1.0 - t);
		}
		else {
		    /* A Bezier where Crv -> Length == 2 so t is it. */
		}

		if (NoSolsOnEndPts &&
		    (IRIT_APX_EQ_EPS(t, TMin, Epsilon) ||
		     IRIT_APX_EQ_EPS(t, TMax, Epsilon))) {
		    /* Ignore this solution on the boundary. */
		}
		else
		    SymbInsertNewParam(t);
	    }
	}
    }
    else {
        /* Use the multivariate solver for this univariate case... */
        if (IRIT_ABS(Epsilon) > SYMB_SUBDIV_TOL)
	    Epsilon = IRIT_SIGN(Epsilon) * SYMB_SUBDIV_TOL * 0.01;

        Pts = MvarCrvZeroSet(Crv, Axis, SYMB_SUBDIV_TOL,
			     Epsilon, TRUE);

	while (Pts != NULL) {
	    IRIT_LIST_POP(Pt, Pts);

	    if (NoSolsOnEndPts &&
		(IRIT_APX_EQ_EPS(Pt -> Pt[0], TMin, Epsilon) ||
		 IRIT_APX_EQ_EPS(Pt -> Pt[0], TMax, Epsilon))) {
	        /* Ignore this solution on the boundary. */
	    }
	    else
	        SymbInsertNewParam(Pt -> Pt[0]);

	    CagdPtFree(Pt);
	}
    }

    return GlblPtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the extrema set of a given curve, in given axis (0/1-3 for      M
* W/X-Z).								     M
*   Returned is a list of the extreme set points holding the parameter       M
* values at Pt[0] of each point.					     M
*   One could compute the derivative of the curve and find its zero set.     M
*   However, for rational curves, this will double the degree and slow down  M
* the computation considerably.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:      To compute its extrema set.                                    M
*   Axis:     The axis of Crv to compute extrema set for, W = 0, X = 1, etc. M
*   Epsilon:  Tolerance control.                                             M
*   NoSolsOnEndPts: If TRUE, solutions at the end of the domain are purged.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:   List of parameter values form which Crv has an         M
*                     extrema value in axis Axis.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvExtremSet, extrema set, symbolic computation                      M
*****************************************************************************/
CagdPtStruct *SymbCrvExtremSet(const CagdCrvStruct *Crv,
			       int Axis,
			       CagdRType Epsilon,
			       CagdBType NoSolsOnEndPts)
{
    CagdCrvStruct *CrvW, *CrvX, *CrvY, *CrvZ, *DCrv, *NewCrv,
	*ScalarCrv = NULL;
    CagdPtStruct *PtList;

    assert(Axis >= 0 && Axis <= 3);

    GlblSetEpsilon = Epsilon;

    if (Crv -> Order < 2)
        return NULL;
    else if (Crv -> Order == 2) {
        /* Find extrema by traversing control polygon - curve is a polyline. */
        int i,
	    Len = Crv -> Length;
        const CagdRType *KV, *R;

	GlblPtList = NULL;

	if (CAGD_IS_RATIONAL_CRV(Crv)) {
	    NewCrv = CagdCoerceCrvTo(Crv, CAGD_PT_E3_TYPE, FALSE);
	    R = NewCrv -> Points[Axis];
	    KV = &NewCrv -> KnotVector[1];
	}
	else {
	    NewCrv = NULL;
	    R = Crv -> Points[Axis];
	    KV = &Crv -> KnotVector[1];
	}

	for (i = 1; i < Len - 1; i++) {
	    if ((R[i] > R[i - 1] && R[i] >= R[i + 1]) ||
		(R[i] < R[i - 1] && R[i] <= R[i + 1]))
	        SymbInsertNewParam(KV[i]);
	}

	if (NewCrv != NULL)
	    CagdCrvFree(NewCrv);

	return GlblPtList;
    }

    SymbCrvSplitScalar(Crv, &CrvW, &CrvX, &CrvY, &CrvZ);

    switch (Axis) {
	case 0:
	    if (CrvW)
		ScalarCrv = CrvW;
	    else
		SYMB_FATAL_ERROR(SYMB_ERR_OUT_OF_RANGE);
	    break;
	case 1:
	    if (CrvX)
		ScalarCrv = CrvX;
	    else
		SYMB_FATAL_ERROR(SYMB_ERR_OUT_OF_RANGE);
	    break;
	case 2:
	    if (CrvY)
		ScalarCrv = CrvY;
	    else
		SYMB_FATAL_ERROR(SYMB_ERR_OUT_OF_RANGE);
	    break;
	case 3:
	    if (CrvZ)
		ScalarCrv = CrvZ;
	    else
		SYMB_FATAL_ERROR(SYMB_ERR_OUT_OF_RANGE);
	    break;
	default:
	    SYMB_FATAL_ERROR(SYMB_ERR_OUT_OF_RANGE);
    }

    if (Axis > 0) {
        NewCrv = SymbCrvMergeScalar(CrvW, ScalarCrv, NULL, NULL);
	if (CrvW)
	    CagdCrvFree(CrvW);
    }
    else
        NewCrv = ScalarCrv;

    if (CrvX)
	CagdCrvFree(CrvX);
    if (CrvY)
	CagdCrvFree(CrvY);
    if (CrvZ)
	CagdCrvFree(CrvZ);

    DCrv = CagdCrvDerive(NewCrv);

    PtList = SymbCrvZeroSet(DCrv, 1, Epsilon, NoSolsOnEndPts);

    CagdCrvFree(NewCrv);
    CagdCrvFree(DCrv);

    return PtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the constant set of a given curve, in the given axis (1-3 for   M
* X-Z).								 	     M
*   Returned is a list of the constant set points holding the parameter	     M
* values at Pt[0] of each point.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:      To compute its constant set.                                   M
*   Axis:     The axis of Crv to compute constant set for, X = 1, Y = 2, etc.M
*   Epsilon:  Tolerance control.                                             M
*   ConstVal:   The value at which to compute the constant set.		     M
*   NoSolsOnEndPts: If TRUE, solutions at the end of the domain are purged.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:   List of parameter values form which Crv has an         M
*                     value of ConstVal in axis Axis.                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvConstSet, constant set, zero set, symbolic computation            M
*****************************************************************************/
CagdPtStruct *SymbCrvConstSet(const CagdCrvStruct *Crv,
			      int Axis,
			      CagdRType Epsilon,
			      CagdRType ConstVal,
			      CagdBType NoSolsOnEndPts)
{
    CagdPType Trans;
    CagdCrvStruct
        *TCrv = CagdCrvCopy(Crv);
    CagdPtStruct *PtList;

    assert(Axis >= 1 && Axis <= 3);

    IRIT_PT_RESET(Trans);
    Trans[Axis - 1] = -ConstVal;
    CagdCrvTransform(TCrv, Trans, 1.0);

    PtList = SymbCrvZeroSet(TCrv, Axis, Epsilon, NoSolsOnEndPts);

    CagdCrvFree(TCrv);

    return PtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the zeros of low degree polynomial, analytically.               M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:   Low degree polynomial to derive its roots analytically.           M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:  list of zeros.                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbScalarCrvLowDegZeroSet                                               M
*****************************************************************************/
CagdPtStruct *SymbScalarCrvLowDegZeroSet(CagdCrvStruct *Crv)
{
    int j = 0;
    CagdRType TMin, TMax;
    CagdCrvStruct
        *BzrCrv = CAGD_IS_BEZIER_CRV(Crv) ? Crv
					  : CagdCnvrtBsp2BzrCrv(Crv),
        *PwrCrv = CagdCnvrtBzr2PwrCrv(BzrCrv);
    CagdRType Sols[4],
        *Pts = PwrCrv -> Points[1];

    GlblPtList = NULL;

    CagdCrvDomain(Crv, &TMin, &TMax);

    switch (Crv -> Order) {
        case 5:
	    if (IRIT_FABS(Pts[4]) > IRIT_EPS) {
	        j = GMSolveQuarticEqn(Pts[3] / Pts[4],
				      Pts[2] / Pts[4],
				      Pts[1] / Pts[4],
				      Pts[0] / Pts[4],
				      Sols);
		break;
	    }
        case 4:
	    if (IRIT_FABS(Pts[3]) > IRIT_EPS) {
	        j = GMSolveCubicEqn(Pts[2] / Pts[3],
				    Pts[1] / Pts[3],
				    Pts[0] / Pts[3],
				    Sols);
		break;
	    }
        case 3:
	    if (IRIT_FABS(Pts[2]) > IRIT_EPS) {
	        j = GMSolveQuadraticEqn(Pts[1] / Pts[2],
					Pts[0] / Pts[2],
					Sols);
		break;
	    }
        case 2:
	    if (IRIT_FABS(Pts[1]) > IRIT_EPS) {
	        j = 1;
		Sols[0] = -Pts[0] / Pts[1];
	    }
	    break;
    }
    CagdCrvFree(PwrCrv);

    while (--j >= 0) {
        if (Sols[j] > -GlblSetEpsilon && Sols[j] <= 1.0 + GlblSetEpsilon) {
	    Sols[j] = TMin + (TMax - TMin) * Sols[j];
	    SymbInsertNewParam(Sols[j]);
	}
    }

    if (BzrCrv != Crv)
        CagdCrvFree(BzrCrv);

    return GlblPtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Returns TRUE iff the Crv is not rational or rational with weights that are M
* entirely positive or entirely negative.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:       To examine for same sign weights, if any.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:   TRUE if no weights or all of same sign.                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvPosNegWeights, symbolic computation                               M
*****************************************************************************/
CagdBType SymbCrvPosNegWeights(const CagdCrvStruct *Crv)
{
    int i;
    CagdBType HasNeg, HasPos;
    CagdRType
	*Weights = Crv -> Points[0];

    if (Weights == NULL)
	return FALSE;				   /* Curve is not rational. */

    for (HasNeg = HasPos = FALSE, i = Crv -> Length - 1; i >= 0; i--) {
	HasNeg |= *Weights < 0.0;
	HasPos |= *Weights++ > 0.0;
    }

    return HasNeg && HasPos;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the poles of a rational surface, solving for the zeros of the   M
* surface's denominator.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:          Rational curves to extract its poles.                      M
*   Epsilon:	  The numerical tolerance to use.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:   The poles, as piecewise linear approximations.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdPointsHasPoles, MvarRationalCrvsPoles                                M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbSplitRationalCrvsPoles                                               M
*****************************************************************************/
CagdPtStruct *SymbSplitRationalCrvsPoles(const CagdCrvStruct *Crv,
					 CagdRType Epsilon)
{
    CagdPtStruct *Poles;

    if (!CAGD_IS_RATIONAL_SRF(Crv)) {
        return NULL;
    }

    Poles = SymbCrvZeroSet(Crv, 0, Epsilon, FALSE);

    return Poles;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Splits the given rational curves at all their poles.  Returned is a      M
* list of curves each of which has weights of the same (positive) sign.      M
*                                                                            *
* PARAMETERS:                                                                M
*   Crvs:      Rational curves to split at all its poles.		     M
*   Eps:       Tolerance of computation.                                     M
*   OutReach:  Clip end points of curves that goes to infinity at distance   M
*	       that is about OutReach from the origin.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   List of splitted curves.			             M
*                                                                            *
* SEE ALSO:                                                                  M
*    CagdPointsHasPoles, SymbCrvSplitPoleParams                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvsSplitPoleParams                                                  M
*****************************************************************************/
CagdCrvStruct *SymbCrvsSplitPoleParams(const CagdCrvStruct *Crvs,
				       CagdRType Eps,
				       CagdRType OutReach)
{
    const CagdCrvStruct *Crv;
    CagdCrvStruct
	*NoPolesCrvs = NULL;

    for (Crv = Crvs; Crv != NULL; Crv = Crv -> Pnext) {
        CagdCrvStruct
	    *TCrvs = SymbCrvSplitPoleParams(Crv, Eps, OutReach);

	NoPolesCrvs = CagdListAppend(TCrvs, NoPolesCrvs);
    }

    return NoPolesCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Splits the given rational curve at all its poles.  Returned is a list of M
* curves each of which has weights of the same (positive) sign.              M
*                                                                            *
* PARAMETERS:                                                                M
*   CCrv:      Rational curve to split at all its poles.		     M
*   Eps:       Tolerance of computation.                                     M
*   OutReach:  Clip end points of curves that goes to infinity at distance   M
*	       that is about OutReach from the origin.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   List of splitted curves.			             M
*                                                                            *
* SEE ALSO:                                                                  M
*    CagdPointsHasPoles, SymbCrvsSplitPoleParams                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvSplitPoleParams                                                   M
*****************************************************************************/
CagdCrvStruct *SymbCrvSplitPoleParams(const CagdCrvStruct *CCrv,
				      CagdRType Eps,
				      CagdRType OutReach)
{
    int i, k, Proximity,
	MaxCoord = CAGD_NUM_OF_PT_COORD(CCrv -> PType);
    CagdRType OrigTMin, OrigTMax, TMin, TMax, **Points, *WPts;
    CagdPtStruct *Pts;
    CagdCrvStruct *Crvs, *TCrv,
	*OutCrvs = NULL,
	*Crv = CagdCrvCopy(CCrv);

    if (!CAGD_IS_RATIONAL_CRV(Crv))
	return Crv;

    CagdCrvDomain(Crv, &OrigTMin, &OrigTMax);

    Points = Crv -> Points;
    WPts = Points[0];
    for (i = 0; i < Crv -> Length; i++) {
	if (IRIT_APX_EQ_EPS(WPts[i], 0.0, Eps)) {
	    for (k = 0; k <= MaxCoord; k++)
		Points[k][i] = 0.0;
	}
    }

    /* Possible to have poles - do we have both positive & negative weights? */
    if (!CagdPointsHasPoles(Crv -> Points, Crv -> Length) ||
	(Pts = SymbSplitRationalCrvsPoles(Crv, Eps)) == NULL)
        return Crv;

    Crvs = CagdCrvSubdivAtParams(Crv, Pts, Eps, &Proximity);
    CagdPtFreeList(Pts);
    CagdCrvFree(Crv);

    while (Crvs != NULL) {
        int i;
        CagdRType t;

	IRIT_LIST_POP(Crv, Crvs);
	WPts = Crv -> Points[0];

        /* Make sure all the weights are positive. */
        for (i = 0, t = 0.0; i < Crv -> Length; i++)
	    if (IRIT_FABS(t) < IRIT_FABS(WPts[i]))
	        t = WPts[i];

	if (t < 0.0) {
	    int j,
		MaxAxis = CAGD_NUM_OF_PT_COORD(Crv -> PType);
	    CagdRType
		**Points = Crv -> Points;

	    /* We are flipping all signs of all coefficients... */
	    for (i = 0; i < Crv -> Length; i++) {
		for (j = 1; j <= MaxAxis; j++)
		    Points[j][i] = -Points[j][i];
	    }
	}

	for (i = 0; i < Crv -> Length; i++)
	    WPts[i] = IRIT_FABS(WPts[i]);

	/* Clip end domain, following OutReach. */
	CagdCrvDomain(Crv, &TMin, &TMax);
	if (OutReach > 0.0 && TMax - TMin > 2 * OutReach) {
	    TCrv = CagdCrvRegionFromCrv(Crv,
			IRIT_APX_EQ(OrigTMin, TMin) ? TMin : TMin + OutReach,
			IRIT_APX_EQ(OrigTMax, TMax) ? TMax : TMax - OutReach);
	    IRIT_LIST_PUSH(TCrv, OutCrvs);
	}
	CagdCrvFree(Crv);
    }

    return OutCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Clear the global sorted list of points and return what was on that list    M
* before.  This sorted list is built via SymbInsertNewParam		     M
*                                                                            *
* PARAMETERS:                                                                M
*   None								     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:	The old point list.                                  M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbInsertNewParam                                                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbGetParamListAndReset		                                     M
*****************************************************************************/
CagdPtStruct *SymbGetParamListAndReset(void)
{
    CagdPtStruct
	*PtList = GlblPtList;

    GlblPtList = NULL;

    return PtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Insert a single t value into existing GlblPtList, provided no equal t      M
* value exists already in the list. List is ordered incrementally.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   t:         New value to insert to global GlblPtList list.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbGetParamListAndReset, SymbInsertNewParam2                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbInsertNewParam    		                                     M
*****************************************************************************/
void SymbInsertNewParam(CagdRType t)
{
    GlblPtList = SymbInsertNewParam2(GlblPtList, t);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Insert a single t value into given PtList, in place, provided no equal t   M
* value exists already in the list. List is ordered incrementally.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   PtList:    List to insert a new value t into.                            M
*   t:         New value to insert to PtList list.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:  Updated list, in place.                                 M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbGetParamListAndReset, SymbInsertNewParam                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbInsertNewParam2    		                                     M
*****************************************************************************/
CagdPtStruct *SymbInsertNewParam2(CagdPtStruct *PtList, CagdRType t)
{
    CagdPtStruct *PtTmp, *PtLast, *Pt;

    Pt = CagdPtNew();
    Pt -> Pt[0] = t;

    if (PtList) {
	for (PtTmp = PtList, PtLast = NULL;
	     PtTmp != NULL;
	     PtLast = PtTmp, PtTmp = PtTmp -> Pnext) {
	    if (SYMB_ZERO_IRIT_APX_EQ(PtTmp -> Pt[0], t)) {
	        IritFree(Pt);
		return PtList;
	    }
	    if (PtTmp -> Pt[0] > t)
	        break;
	}
	if (PtTmp) {
	    /* Insert the new point in the middle of the list. */
	    Pt -> Pnext = PtTmp;
	    if (PtLast)
		PtLast -> Pnext = Pt;
	    else
		PtList = Pt;
	}
	else {
	    /* Insert the new point as the last point in the list. */
	    PtLast -> Pnext = Pt;
	}
    }
    else
        PtList = Pt;

    return PtList;
}
