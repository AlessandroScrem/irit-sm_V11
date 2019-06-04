/******************************************************************************
* symb_cci.c - curve curve intersection related functions.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Fady Massarwi and Gershon Elber, July. 2014.		      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/symb_lib.h"

/* #define DEBUG_CALC_SUB_REGIONS_A_MAT */
/* #define DEBUG_PRINT_A_MATRIX */

static void SymbCalcCrvDerivativeAtPoint(const CagdCrvStruct *Crv1,
					 CagdRType Param,
					 CagdVType Res);
static int SymbCalcCrvsIntersectionDomains(const CagdCrvStruct *Crv1,
					   const CagdCrvStruct *Crv2,
					   CagdRType Eps,
					   CagdRType **InterDomains);
static void SymbGet2CrvsInterDAreaDCtlPtsAux(CagdRType *dAreadPts,
					    CagdCrvStruct *Crv1,
					    CagdCrvStruct *Crv2,
					    CagdRType Eps,
					    CagdCrvStruct *RefCrv);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Calculates the derivative (tangent vector) of a planar curve at          *
* a given point                                                              *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv1:       A planar curve (assume not NULL).                            *
*   Param:      Parameter value of the point.                                *
*   Res:        Output vector (Tu,Tv,0).                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void				  			             *
*****************************************************************************/
static void SymbCalcCrvDerivativeAtPoint(const CagdCrvStruct *Crv1,
					 CagdRType Param,
					 CagdVType Res)
{
    CagdCrvStruct
        *DerivCrv = CagdCrvDerive(Crv1);
    CagdRType
        *DrvVal = CagdCrvEval(DerivCrv, Param);

    Res[0] = DrvVal[1];
    Res[1] = DrvVal[2];
    Res[2] = 0.0;

    CagdCrvFree(DerivCrv);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Calculates the intersection domains of two curves, the curves must have  *
* the same orientation.                                                      *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv1:         First planar curve.                                        *
*   Crv2:         Second planar curve.                                       *
*   Eps:         Tolerance of computation.                                   *
*   InterDomains: Calculated intersecting domain(s) as a list of the         *
*                 following 4-tuple structure:				     *
*                              (u1_enter, u1_leave, u2_enter, u2_leave),     *
*                 sorted by u1 values.                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: Number of intersection domains.			             *
*****************************************************************************/
static int SymbCalcCrvsIntersectionDomains(const CagdCrvStruct *Crv1,
					   const CagdCrvStruct *Crv2,
					   CagdRType Eps,
					   CagdRType **InterDomains)
{
    int NumIntersectionDomains, i, IsShifted;
    CagdRType UStart1,UEnd1, UStart2, UEnd2;
    CagdVType TangentC1Vec, TangentC2Vec;
    CagdPtStruct *IntersectPts, *IntersectPtsHead;

    if (Crv1 == NULL || Crv2 == NULL) {
        *InterDomains = NULL;
        return 0;
    }

    IntersectPtsHead = IntersectPts = CagdCrvCrvInter(Crv1, Crv2, Eps);
    NumIntersectionDomains = CagdListLength(IntersectPts);

    /* Require always even number of intersection points. */
    if (NumIntersectionDomains == 0 || (NumIntersectionDomains % 2) != 0) {
	if (InterDomains) {
	    IritFree(InterDomains);
	}
	if (IntersectPtsHead)
	    CagdPtFreeList(IntersectPtsHead);
	return 0;
    }

    SymbCalcCrvDerivativeAtPoint(Crv1, IntersectPts -> Pt[0], TangentC1Vec);
    SymbCalcCrvDerivativeAtPoint(Crv2, IntersectPts -> Pt[1], TangentC2Vec);
    IsShifted = IRIT_CROSS_PROD_2D(TangentC1Vec, TangentC2Vec) > 0 ? 1 : 0;
    
    if (IsShifted) {
        /* Make the first intersection point last. */
	while (IntersectPts -> Pnext != NULL) {
	    IntersectPts = IntersectPts -> Pnext;
	}

	IntersectPts -> Pnext = IntersectPtsHead;
	IntersectPtsHead = IntersectPtsHead -> Pnext;
	IntersectPts -> Pnext -> Pnext = NULL;
    }

    /* Count intersection domains and not intersection points. */
    NumIntersectionDomains = NumIntersectionDomains / 2;

    /* (u1_start, u1_end, u2_start, u2_end) sorted according to u1. */    
    *InterDomains = IritMalloc(NumIntersectionDomains * 4 * sizeof(CagdRType));
    i = 0;
    IntersectPts = IntersectPtsHead;
    while (IntersectPts != NULL) {
        UStart1 = IntersectPts -> Pt[0];
        UStart2 = IntersectPts -> Pt[1];
        IntersectPts = IntersectPts -> Pnext;
        if (IntersectPts != NULL) {
            UEnd1 = IntersectPts -> Pt[0];
            UEnd2 = IntersectPts -> Pt[1];
            (*InterDomains)[i++] = UStart1;
            (*InterDomains)[i++] = UEnd1;
            (*InterDomains)[i++] = UStart2;
            (*InterDomains)[i++] = UEnd2;
            IntersectPts = IntersectPts -> Pnext;
        }
	else {
	    assert(0);     /* should have had even number of intersections. */
	}
    }
    
    if (IntersectPtsHead)
        CagdPtFreeList(IntersectPtsHead);

    return NumIntersectionDomains;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the intersection region(s) of given two curves if any.          M
* curves. The curves must have the same orientation.                         M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1, Crv2:  The two curves to intersect and compute intersecting        M
*                region(s).                                                  M
*   Eps:         Tolerance of computation.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:  List of closed intersecting region(s) or NULL if none. M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCrvCrvInter                                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbGet2CrvsIntersectionRegions                                          M
*****************************************************************************/
CagdCrvStruct *SymbGet2CrvsIntersectionRegions(const CagdCrvStruct *Crv1,
					       const CagdCrvStruct *Crv2,
					       CagdRType Eps)
{    
    CagdRType *InterDomains;    
    int NumIntersectionDomains, i;
    CagdCrvStruct
        *RetCrvs = NULL;

    NumIntersectionDomains = SymbCalcCrvsIntersectionDomains(Crv1, Crv2, Eps,
							     &InterDomains);

    if (NumIntersectionDomains == 0)
	return 0;

    for (i = 0; i < NumIntersectionDomains; ++i) {
        CagdRType *R, TMin, TMax;
	CagdPType Pt1End, Pt2Start, Pt2End;
        CagdCrvStruct *MergedCloseCrv,
	    *C1Part = CagdCrvRegionFromCrv(Crv1, 
					   InterDomains[4 * i + 0], 
					   InterDomains[4 * i + 1]),
	    *C2Part = CagdCrvRegionFromCrv(Crv2, 
					   InterDomains[4 * i + 2], 
					   InterDomains[4 * i + 3]);

	CagdCrvDomain(C1Part, &TMin, &TMax);
	R = CagdCrvEval(C1Part, TMax);
	CagdCoerceToE3(Pt1End, &R, -1, C1Part -> PType);

	CagdCrvDomain(C2Part, &TMin, &TMax);
	R = CagdCrvEval(C2Part, TMin);
	CagdCoerceToE3(Pt2Start, &R, -1, C1Part -> PType);
	R = CagdCrvEval(C2Part, TMax);
	CagdCoerceToE3(Pt2End, &R, -1, C1Part -> PType);

	if (IRIT_PT_PT_DIST_SQR(Pt1End, Pt2End) < 
	    IRIT_PT_PT_DIST_SQR(Pt1End, Pt2Start)) {
	    CagdCrvStruct
		*C2PartRev = CagdCrvReverse(C2Part);     /* Reverse C2Part. */

	    MergedCloseCrv = CagdMergeCrvCrv(C1Part, C2PartRev, FALSE);
	    CagdCrvFree(C2PartRev);
	}
	else
	    MergedCloseCrv = CagdMergeCrvCrv(C1Part, C2Part, FALSE);

	IRIT_LIST_PUSH(MergedCloseCrv, RetCrvs);

	CagdCrvFree(C1Part);
	CagdCrvFree(C2Part);
    }

    if (InterDomains)
        IritFree(InterDomains);

    return RetCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Calculates the total area of the intersection domains of two planar      M
* curves. The curves must have the same orientation.                         M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1:       First planar curve.                                          M
*   Crv2:       Second planar curve.                                         M
*   Eps:         Tolerance of computation.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:   Total area of intersection domains of the two curves.       M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbCrvEnclosedAreaEval, SymbGet2CrvsInterDAreaDCtlPts                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbGet2CrvsIntersectionAreas                                            M
*****************************************************************************/
CagdRType SymbGet2CrvsIntersectionAreas(const CagdCrvStruct *Crv1,
					const CagdCrvStruct *Crv2,
					CagdRType Eps)
{    
    CagdCrvStruct *Crv,
        *InterCrvs = SymbGet2CrvsIntersectionRegions(Crv1, Crv2, Eps);
    CagdRType
        TotalAreaSum = 0.0;

    for (Crv = InterCrvs; Crv != NULL; Crv = Crv -> Pnext) {
        TotalAreaSum += fabs(SymbCrvEnclosedAreaEval(Crv));
    }

    CagdCrvFreeList(InterCrvs);

    return TotalAreaSum;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   An auxiliary function for SymbGet2CrvsInterDAreaDCtlPts. Computes the    *
* change rate of the intersection regions of two curves depending on a       *
* change of a reference curve control points.                                *
*                                                                            *
* PARAMETERS:                                                                *
*   dAreadPts:  Output array of the change rate per each axis of the CP.     *
*   Crv1:       First planar curve.                                          *
*   Crv2:       Second planar curve.                                         *
*   RefCrv:     Reference curve.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void SymbGet2CrvsInterDAreaDCtlPtsAux(CagdRType *dAreadPts,
					     CagdCrvStruct *Crv1,
					     CagdCrvStruct *Crv2,
					     CagdRType Eps,
					     CagdCrvStruct *RefCrv)
{
    int i, j;

    for (i = 0; i < RefCrv -> Length; ++i) {
        for (j = 1; j <= 2; ++j) {
            RefCrv -> Points[j][i] += IRIT_EPS;

            dAreadPts[2 * i + j - 1] =
	                       SymbGet2CrvsIntersectionAreas(Crv1, Crv2, Eps);

            RefCrv -> Points[j][i] -= 2 * IRIT_EPS;

            dAreadPts[2 * i + j - 1] =
		(dAreadPts[2 * i + j - 1] - 
		 SymbGet2CrvsIntersectionAreas(Crv1, Crv2, Eps)) /
	                                                       (2 * IRIT_EPS);

            RefCrv -> Points[j][i] += IRIT_EPS;
        }
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Calculates the intersection domains of two curves, and the change rate   M
* of the total intersection area relative to control points change.          M
*   The curves must be planar and have the same orientation. 		     M
*   The intersections must be complete - that is, there are even number of   M
* intersections.					                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1:         First planar curve.                                        M
*   Crv2:         Second planar curve.                                       M
*   Eps:         Tolerance of computation.                                   M
*   InterDomains: Calculated domain as a list of the in the following        M
*                 4-tuple structure,					     M
*		      (u1_enter, u1_leave, u2_enter, u2_leave)               M
*                 sorted by u1 values.                                       M
*   dAreadPts:    Will be updated to contain the change rate of the          M
*                 intersection areas relative to changes in each control     M
*                 points of both curves.  The order in this array is as      M
*                 follows,					             M
*                 ( Crv1_ctpnt1_x, Crv1_ctpnt1_y,			     M
*                   Crv1_ctpnt2_x, Crv1_ctpnt2_y, ...                        M
*                   Crv2_ctpnt1_x, Crv2_ctpnt1_y, ...                        M
*                   Crv2_ctpnt2_x, Crv2_ctpnt2_y, ... ).                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int: Number of intersection domains of the two curves                    M
*        (the length of InterDomains is 4-times this return value).          M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbGet2CrvsIntersectionAreas 			                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbGet2CrvsInterDAreaDCtlPts                                            M
*****************************************************************************/
int SymbGet2CrvsInterDAreaDCtlPts(CagdCrvStruct *Crv1,
				  CagdCrvStruct *Crv2,
				  CagdRType Eps,
				  CagdRType **InterDomains,
				  CagdRType **dAreadPts)
{
    int DAreaListLen, NumIntersectionDomains;

    if (Crv1 == NULL ||
	Crv2 == NULL ||
	(NumIntersectionDomains =
	 SymbCalcCrvsIntersectionDomains(Crv1, Crv2, Eps,
					 InterDomains)) == 0) {
        *InterDomains = NULL;
        return 0;
    }

    DAreaListLen = (Crv1 -> Length + Crv2 -> Length) * 2;
    *dAreadPts = IritMalloc(DAreaListLen * sizeof(CagdRType));
    SymbGet2CrvsInterDAreaDCtlPtsAux(*dAreadPts, Crv1, Crv2, Eps, Crv1);
    SymbGet2CrvsInterDAreaDCtlPtsAux((*dAreadPts) + 2 * Crv1 -> Length,
				     Crv1, Crv2, Eps, Crv2);

#   ifdef DEBUG_CALC_SUB_REGIONS_A_MAT
    {
        int i;

	for (i = 0; i < NumIntersectionDomains; ++i) {
	    int Dim[2], j, k, l;
	    CagdRType
	        *A = SymbGetCrvSubRegionAlphaMatrix(Crv1,
						    (*InterDomains)[i * 4],
						    (*InterDomains)[i * 4 + 1],
						    Dim);

	    fprintf(stderr, "Extraction matrix (Subdomain [%8.5f :: %8.5f]):\n",
		    (*InterDomains)[i * 4], (*InterDomains)[i * 4 + 1]);
	    for (k = 0, j = 0; j < Dim[1]; j++) {
	        for (l = 0; l < Dim[0]; l++) {
		    fprintf(stderr, "%8.5f  ", A[k++]);
		}
		fprintf(stderr, "\n");
	    }
	}
    }
#   endif /* DEBUG_CALC_SUB_REGIONS_A_MAT */

    return NumIntersectionDomains;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the Alpha matrix that relates the control points of the curve   M
* region (t1, t2) to the control points of the input curve Crv of a larger   M
* domain.								     M
*   Compute the Alpha matrix by adding Order-1 knots at t1 and t2.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:      Curve to compute its region extraction Alpha matrix.           M
*             Assumed a non periodic curve.                                  M
*   t1, t2:   Of extracted domain.	                                     M
*   Dim:      The two dimensions of the returned matrix will be placed here. M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType *: Linearized 2D matrix, row by row, of size (CrvLen, RgnLen), M
*                where CrvLen is the number of control points in curve Crv,  M
*                and RgnLen is the size of the extracted curve region.       M
*                  Allocated dynamically.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspKnotEvalAlphaCoef				                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbGetCrvSubRegionAlphaMatrix                                           M
*****************************************************************************/
CagdRType *SymbGetCrvSubRegionAlphaMatrix(const CagdCrvStruct *Crv,
					  CagdRType t1,
					  CagdRType t2,
					  int *Dim)
{
    int i, j, k, Idx1, Idx2, Mult1, Mult2, RefTLen, MergedKVLen,
        Length = Crv -> Length,
        Order = Crv -> Order;
    CagdRType *MergedKV, *RefT, *RetMat,
        *KV = Crv -> KnotVector;
    BspKnotAlphaCoeffStruct *A;

    if (t1 > t2)
	IRIT_SWAP(CagdRType, t1, t2);

#   ifdef DEBUG
    {
        CagdRType TMin, TMax;

	BspCrvDomain(Crv, &TMin, &TMax);
	assert(TMin <= t1 && t1 < t2 && t2 <= TMax);
    }
#   endif /* DEBUG */

    /* Compute the current multiplicities of the knots at t1 and t2. */
    i = BspKnotLastIndexL(KV, Length + Order, t1);
    if (IRIT_APX_EQ(KV[i], t1))
	Idx1 = i;
    else if (IRIT_APX_EQ(KV[i + 1], t1))
	Idx1 = i + 1;
    else
        Idx1 = -1;
    Mult1 = Idx1 > 0 ? BspKnotMultiplicity(KV, Length + Order, Idx1) : 0;

    i = BspKnotLastIndexL(KV, Length + Order, t2);
    if (IRIT_APX_EQ(KV[i], t2))
	Idx2 = i;
    else if (IRIT_APX_EQ(KV[i + 1], t2))
	Idx2 = i + 1;
    else
        Idx2 = -1;
    Mult2 = Idx2 > 0 ? BspKnotMultiplicity(KV, Length + Order, Idx2) : 0;

    /* Build refinement knot vector and compute the alpha matrix for it. */
    assert(Mult1 <= Order - 1 && Mult2 <= Order - 1);
    RefTLen = (Order - 1) * 2 - Mult1 - Mult2;
    RefT = (CagdRType *) IritMalloc(sizeof(CagdRType) * RefTLen);
    for (i = 0; i < Order - 1 - Mult1; i++)
        RefT[i] = t1;
    for ( ; i < RefTLen; i++)
        RefT[i] = t2;

    MergedKV = BspKnotMergeTwo(KV, Length + Order,
			       RefT, RefTLen, 0, &MergedKVLen);
    A = BspKnotEvalAlphaCoef(Order, KV, Length, MergedKV,
			     MergedKVLen - Order, FALSE);

#   ifdef DEBUG_PRINT_A_MATRIX
    CagdDbgPrintAlphaMat(A);
#   endif /* DEBUG_PRINT_A_MATRIX */

    /* Sets control point indices of relevant extracted region and set Dmn. */
    Idx1 = BspKnotLastIndexL(MergedKV, MergedKVLen, t1);
    Idx2 = BspKnotLastIndexL(MergedKV, MergedKVLen, t2);
    Dim[0] = A -> Length;
    Dim[1] = Idx2 - Idx1 + 1;

    RetMat = (CagdRType *) IritMalloc(sizeof(CagdRType) * Dim[0] * Dim[1]);
    IRIT_ZAP_MEM(RetMat, sizeof(CagdRType) * Dim[0] * Dim[1]);

    for (k = 0, j = Idx1; j <= Idx2; j++) {
	for (i = 0; i < A -> Length; i++) {
	    RetMat[k++] = A -> Rows[i][j];
	}
    }

    IritFree(RefT);
    IritFree(MergedKV);
    BspKnotFreeAlphaCoef(A);

    return RetMat;
}
