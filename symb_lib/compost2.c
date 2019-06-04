/******************************************************************************
* Compost2.c - srf-srf composition computation (symbolic).		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Jul 2011.					      *
******************************************************************************/

#include <string.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/allocate.h"
#include "inc_irit/iritprsr.h"
#include "symb_loc.h"

#define SYMB_VERIFY_E2_CRV(Crv) \
    if (CAGD_NUM_OF_PT_COORD(Crv -> PType) < 2) { \
	CagdCrvStruct \
	  *CTmp = CagdCoerceCrvTo(Crv, CAGD_PT_E2_TYPE, FALSE);	\
	CagdCrvFree(Crv); \
	Crv = CTmp; \
    }

static CagdSrfStruct *SymbComposeSrfSrfAux(const CagdSrfStruct *Srf1,
					   const CagdSrfStruct *Srf2);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given surfaces Srf2 and Srf1, computes the composition Srf1(Srf2(u, v)).   M
*   Srf2 must be a two dimensional surface completely contained in the       M
* parametric domain of Srf1.                                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf1, Srf2:   The surfaces to compose. Srf1 must be a Bezier.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:    The resulting composition.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbDecomposeCrvCrv, SymbComposeSrfCrv, MvarComposeTVSrf		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbComposeSrfSrf, composition                                           M
*****************************************************************************/
CagdSrfStruct *SymbComposeSrfSrf(const CagdSrfStruct *Srf1,
				 const CagdSrfStruct *Srf2)
{
    switch (Srf1 -> GType) {
	case CAGD_SBEZIER_TYPE:
	    break;
	case CAGD_SBSPLINE_TYPE:
	    SYMB_FATAL_ERROR(SYMB_ERR_BZR_SRF_EXPECT);
	    break;
	case CAGD_SPOWER_TYPE:
	    SYMB_FATAL_ERROR(SYMB_ERR_POWER_NO_SUPPORT);
	    break;
	default:
	    SYMB_FATAL_ERROR(SYMB_ERR_UNDEF_SRF);
	    break;
    }

    switch (Srf2 -> GType) {
	case CAGD_SBEZIER_TYPE:
	    break;
	case CAGD_SBSPLINE_TYPE:
	    break;
	case CAGD_SPOWER_TYPE:
	    SYMB_FATAL_ERROR(SYMB_ERR_POWER_NO_SUPPORT);
	    break;
	default:
	    SYMB_FATAL_ERROR(SYMB_ERR_UNDEF_SRF);
	    break;
    }

    return SymbComposeSrfSrfAux(Srf1, Srf2);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Aux. function. Subdivides Srf2 until it is a Bezier surface, invokes the   *
* Bezier composition code on each, and merges them back to complete surface. *
*   At this point, the surface can be either Bezier or B-spline only.        *
*   Srf2 is assumed to have open end condition.	   		             *
*                                                                            *
* PARAMETERS:                                                                *
*   Srf1, Srf2:   The two surfaces to compose together.			     *
*                 Srf1 is assumed a Bezier surface.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdSrfStruct *:  The composed surface.                                  *
*****************************************************************************/
static CagdSrfStruct *SymbComposeSrfSrfAux(const CagdSrfStruct *Srf1,
					   const CagdSrfStruct *Srf2)
{
    CagdRType t;
    CagdSrfStruct *CmpsSrf2;

    if (Srf2 -> ULength > Srf2 -> UOrder) {
	/* Srf2 is not a Bezier segment in U. Subdivide, compute for each  */
	/* segment and merge back.				           */
	CagdSrfStruct *Srf2A, *Srf2B, *CmpsSrf2A, *CmpsSrf2B;

	t = Srf2 -> UKnotVector[(Srf2 -> UOrder + Srf2 -> ULength) >> 1];

	Srf2A = CagdSrfSubdivAtParam(Srf2, t, CAGD_CONST_U_DIR);
	Srf2B = Srf2A -> Pnext;
	Srf2A -> Pnext = NULL;

	CmpsSrf2A = SymbComposeSrfSrfAux(Srf1, Srf2A);
	CmpsSrf2B = SymbComposeSrfSrfAux(Srf1, Srf2B);
	CagdSrfFree(Srf2A);
	CagdSrfFree(Srf2B);

	CmpsSrf2 = CagdMergeSrfSrf(CmpsSrf2A, CmpsSrf2B,
				   CAGD_CONST_U_DIR, TRUE, FALSE);
	CagdSrfFree(CmpsSrf2A);
	CagdSrfFree(CmpsSrf2B);
    }
    else if (Srf2 -> VLength > Srf2 -> VOrder) {
	/* Srf2 is not a Bezier segment in V. Subdivide, compute for each  */
	/* segment and merge back.				           */
	CagdSrfStruct *Srf2A, *Srf2B, *CmpsSrf2A, *CmpsSrf2B;

	t = Srf2 -> VKnotVector[(Srf2 -> VOrder + Srf2 -> VLength) >> 1];

	Srf2A = CagdSrfSubdivAtParam(Srf2, t, CAGD_CONST_V_DIR);
	Srf2B = Srf2A -> Pnext;
	Srf2A -> Pnext = NULL;

	CmpsSrf2A = SymbComposeSrfSrfAux(Srf1, Srf2A);
	CmpsSrf2B = SymbComposeSrfSrfAux(Srf1, Srf2B);
	CagdSrfFree(Srf2A);
	CagdSrfFree(Srf2B);

	CmpsSrf2 = CagdMergeSrfSrf(CmpsSrf2A, CmpsSrf2B,
				   CAGD_CONST_V_DIR, TRUE, FALSE);
	CagdSrfFree(CmpsSrf2A);
	CagdSrfFree(CmpsSrf2B);
    }
    else {
    	/* Srf2 is a Bezier surface segment - compute its composition. */
	if (!CAGD_IS_BEZIER_SRF(Srf2)) {
	    CagdSrfStruct 
		*TSrf2 = CagdCnvrtBsp2BzrSrf(Srf2);

	    CmpsSrf2 = BzrComposeSrfSrf(Srf1, TSrf2);
	    CagdSrfFree(TSrf2);
	}
	else
    	    CmpsSrf2 = BzrComposeSrfSrf(Srf1, Srf2);
    }

    return CmpsSrf2;    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given surfaces Srf1 and Srf2, computes their composition     	     M
* Srf1(Srf2(a, b)).  Srf1 must be a Bezier.				     M
*   Srf2 must be a two dimensional surface completely contained in the       M
* parametric domain of Srf1.                                                 M
*   See: "Freeform surface analysis using a hybrid of symbolic and numeric   M
* computation" by Gershon Elber, PhD thesis, University of Utah, 1992.	     M
*   Compute the compositions by the products of:			     M
*									     M
* S(u, v) = S(u(a, b), v(a, b))						     V
*									     V
*            n   m							     V
*         = sum sum Pij Bi(u(a, b)) Bj(v(a, b))				     V
*           i=0 j=0							     V
*									     V
*            n   m	 n	     i		    n-i			     V
*         = sum sum Pij ( ) (u(a, b))  (1 - u(a, b))                         V
*           i=0 j=0	 i						     V
*                   	    m		j	       m-j		     V
*                          ( ) (v(a, b))  (1 - v(a, b))                      V
*                  	    j						     V
*									     M
* or if Srf2 is rational:						     M
*									     M
* S(u, v) = S(u(a, b), v(a, b))						     V
*									     V
*            n   m							     V
*         = sum sum Pij Bi(u(a, b)/w(a, b)) Bj(v(a, b)/w(a, b))		     V
*           i=0 j=0							     V
*									     V
*            n   m	 n	     i		          n-i		     V
*         = sum sum Pij ( ) (u(a, b))  (w(a, b) - u(a, b))                   V
*           i=0 j=0	 i						     V
*                   	    m		j	             m-j	     V
*                          ( ) (v(a, b))  (w(a, b) - v(a, b))    /           V
*                  	    j						     V
*                                 	     n+m			     V
*                                   (w(a, b)) 				     V
*                                                                            *
* PARAMETERS:                                                                M
*   Srf1, Srf2:     The two surfaces to compose.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:    The resulting composition.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeSrfSrf, MvarBzrComposeTVSrf				     M
*                                                                            *
* KEYWORDS:                                                                  M
*   BzrComposeSrfSrf, composition                                            M
*****************************************************************************/
CagdSrfStruct *BzrComposeSrfSrf(const CagdSrfStruct *Srf1,
				const CagdSrfStruct *Srf2)
{
    IRIT_STATIC_DATA CagdSrfStruct
	**SrfUFactors = NULL,
	**SrfVFactors = NULL;
    CagdBType
        Is1Rational = CAGD_IS_RATIONAL_SRF(Srf1),
	Is2Rational = CAGD_IS_RATIONAL_SRF(Srf2);
    int i, j, k, l,
        CmpsLength = 0,
	UOrder = Srf1 -> UOrder,
	VOrder = Srf1 -> VOrder,
	MaxCoord = CAGD_NUM_OF_PT_COORD(Srf1 -> PType);
    CagdRType * const
        *SPoints = Srf1 -> Points;
    CagdSrfStruct *SrfUV,
    	*CmpsSrf2 = NULL;

    if (CAGD_NUM_OF_PT_COORD(Srf2 -> PType) < 2) {
	SYMB_FATAL_ERROR(SYMB_ERR_WRONG_PT_TYPE);
	return NULL;
    }

    SrfUV = CagdCoerceSrfTo(Srf2,
			    Is2Rational ? CAGD_PT_P1_TYPE : CAGD_PT_E1_TYPE,
			    FALSE);

    if (SrfUFactors != NULL) {
        for (i = 0; SrfUFactors[i] != NULL; i++)
	    CagdSrfFree(SrfUFactors[i]);
	IritFree(SrfUFactors);
    }
    SrfUFactors = SymbComputeSurfacePowers(SrfUV, UOrder);

    CAGD_GEN_COPY(SrfUV -> Points[1], Srf2 -> Points[2],
		  sizeof(CagdRType) * Srf2 -> ULength * Srf2 -> VLength);

    if (SrfVFactors != NULL) {
        for (i = 0; SrfVFactors[i] != NULL; i++)
	    CagdSrfFree(SrfVFactors[i]);
	IritFree(SrfVFactors);
    }
    SrfVFactors = SymbComputeSurfacePowers(SrfUV, VOrder);

    CagdSrfFree(SrfUV);

    /* The main (double) loop (Compositions computation): */
    for (j = 0; j < VOrder; j++) {
        for (i = 0; i < UOrder; i++) {
	    CagdSrfStruct
	        *TSrf2 = SymbSrfMult(SrfUFactors[i], SrfVFactors[j]);
	    CagdRType
	        *TSrf2Points = TSrf2 -> Points[1];

	    for (k = !Is1Rational; k <= MaxCoord; k++) {
	        CagdRType *CmpsPoints,
		    SPt = SPoints[k][CAGD_MESH_UV(Srf1, i, j)];

		if (i == 0 && j == 0) {
		    if (k == !Is1Rational) {
		        /* Create the resulting surface. */
		        if (CAGD_IS_BEZIER_SRF(TSrf2))
			    CmpsSrf2 = BzrSrfNew(TSrf2 -> ULength,
						 TSrf2 -> VLength,
						 Srf1 -> PType);
			else {
			    CmpsSrf2 = BspSrfNew(TSrf2 -> ULength,
						 TSrf2 -> VLength,
						 TSrf2 -> UOrder,
						 TSrf2 -> VOrder,
						 Srf1 -> PType);
			    CAGD_GEN_COPY(CmpsSrf2 -> UKnotVector,
					  TSrf2 -> UKnotVector,
					  (TSrf2 -> ULength + TSrf2 -> UOrder) *
							     sizeof(CagdRType));
			    CAGD_GEN_COPY(CmpsSrf2 -> VKnotVector,
					  TSrf2 -> VKnotVector,
					  (TSrf2 -> VLength + TSrf2 -> VOrder) *
							     sizeof(CagdRType));
			}
			CmpsLength = CmpsSrf2 -> ULength * CmpsSrf2 -> VOrder;
		    }

		    CmpsPoints = CmpsSrf2 -> Points[k];
		    for (l = 0; l < CmpsLength; l++)
		        CmpsPoints[l] = TSrf2Points[l] * SPt;
		}
		else {
		    CmpsPoints = CmpsSrf2 -> Points[k];
		    for (l = 0; l < CmpsLength; l++)
		        CmpsPoints[l] += TSrf2Points[l] * SPt;
		}
	    }

	    CagdSrfFree(TSrf2);
	}
    }

    if (Is2Rational) {
	CagdSrfStruct *SrfX, *SrfY, *SrfZ, *SrfW, *TSrf,
	    *Srf2W = SymbSrfMult(SrfUFactors[UOrder], SrfVFactors[VOrder]);

	SymbSrfSplitScalar(CmpsSrf2, &SrfW, &SrfX, &SrfY, &SrfZ);
	CagdSrfFree(CmpsSrf2);
	if (SrfW != NULL) {
	    assert(Is1Rational);

	    TSrf = SymbSrfMult(Srf2W, SrfW);
	    CagdSrfFree(SrfW);
	    SrfW = TSrf;
	}
	else
	    SrfW = Srf2W;
	
	CmpsSrf2 = SymbSrfMergeScalar(SrfW, SrfX, SrfY, SrfZ);
	if (SrfX != NULL)
	    CagdSrfFree(SrfX);
	if (SrfY != NULL)
	    CagdSrfFree(SrfY);
	if (SrfZ != NULL)
	    CagdSrfFree(SrfZ);
	if (SrfW != NULL)
	    CagdSrfFree(SrfW);
    }

    for (i = 0; SrfUFactors[i] != NULL; i++)
        CagdSrfFree(SrfUFactors[i]);
    for (i = 0; SrfVFactors[i] != NULL; i++)
        CagdSrfFree(SrfVFactors[i]);

    IritFree(SrfUFactors);
    SrfUFactors = NULL;
    IritFree(SrfVFactors);
    SrfVFactors = NULL;

    return CmpsSrf2;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Computes the factors of the Bernstein polynomials where srf is a scalar    M
* surface, for i from 0 to n (degree):					     M
*									     M
*   n               n                 n-i            i			     V
*  B (srf(u, v)) = ( ) (1 - srf(u, v))    (srf(u, v))			     V
*   i               i							     V
*									     M
*   The surface srf(u, v) is a scalar, possibly a rational surface.	     M
*   If rational, the returned vector, index Order will contain               M
* wsrf(u, v)^n.		 						     M
* See: "Freeform surface analysis using a hybrid of symbolic and numeric     M
* computation" by Gershon Elber, PhD thesis, University of Utah, 1992.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:      Scalar surface to compute factors for.                         M
*   Order:    Order is n + 1.	                                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct **:   A vector of all possible factors for i equal to      M
*                       0 to n, allocated dynamically.                       M
*****************************************************************************/
CagdSrfStruct **SymbComputeSurfacePowers(const CagdSrfStruct *CSrf, int Order)
{
    IRIT_STATIC_DATA int
	LastOrderStaticAlloc = -1;
    IRIT_STATIC_DATA CagdSrfStruct
	**SrfFactors1_Srf = NULL,
	**SrfFactorsSrf = NULL;
    int i;
    CagdBType
	IsRational = CAGD_IS_RATIONAL_SRF(CSrf);
    CagdSrfStruct *Srf_1, *Srf,
	**SrfFactors = (CagdSrfStruct **)
			    IritMalloc((Order + 2) * sizeof(CagdSrfStruct *));
    CagdRType *Points,
	Translate = 0.0;

    if (LastOrderStaticAlloc < Order) {
	LastOrderStaticAlloc = Order * 2 + 1;

	if (SrfFactors1_Srf != NULL) {
	    IritFree(SrfFactors1_Srf);
	    IritFree(SrfFactorsSrf);
	}

	SrfFactors1_Srf = (CagdSrfStruct **)
		  IritMalloc(LastOrderStaticAlloc * sizeof(CagdSrfStruct *));
	SrfFactorsSrf = (CagdSrfStruct **)
		  IritMalloc(LastOrderStaticAlloc * sizeof(CagdSrfStruct *));
    }

    if (IsRational) {
	CagdSrfStruct *SrfY, *SrfZ, *SrfW, *CTmp;

	SymbSrfSplitScalar(CSrf, &SrfW, &Srf, &SrfY, &SrfZ);
	Srf_1 = SymbSrfSub(SrfW, Srf);

	/* Set SrfFactors[Order] to SrfW(u, v)^(Order-1) if rational. */
	SrfFactors[Order] = CagdSrfCopy(SrfW);
	for (i = 1; i < Order - 1; i++) {
	    CTmp = BzrSrfMult(SrfFactors[Order], SrfW);
	    CagdSrfFree(SrfFactors[Order]);
	    SrfFactors[Order] = CTmp;
	}
	SrfFactors[Order + 1] = NULL;

	CagdSrfFree(SrfW);
    }
    else {
	/* Prepare (1.0 - Srf). */
	Srf_1 = CagdSrfCopy(CSrf);
	Srf = CagdSrfCopy(CSrf);
	Points = Srf_1 -> Points[1];
    	for (i = 0; i < Srf -> ULength * Srf -> VLength; i++, Points++)
    	    *Points = 1.0 - *Points;

	SrfFactors[Order] = NULL;
    }

    for (i = 0; i < Order; i++) {
    	if (i == 0) {
    	    SrfFactors1_Srf[0] = NULL;
    	    SrfFactorsSrf[0] = NULL;
	}
    	else if (i == 1) {
    	    SrfFactors1_Srf[1] = Srf_1;
    	    SrfFactorsSrf[1] = CagdSrfCopy(Srf);
    	}
    	else {
    	    SrfFactors1_Srf[i] = SymbSrfMult(SrfFactors1_Srf[i - 1], Srf_1);
    	    SrfFactorsSrf[i] = SymbSrfMult(SrfFactorsSrf[i - 1], Srf);
    	}
    }

    for (i = 0; i < Order; i++) {
	if (i == 0) {
	    SrfFactors[i] = SrfFactors1_Srf[Order - 1];
	}
	else if (i == Order - 1) {
	    SrfFactors[i] = SrfFactorsSrf[Order - 1];
	}
	else {
	    SrfFactors[i] = SymbSrfMult(SrfFactors1_Srf[Order - 1 - i],
					SrfFactorsSrf[i]);
	}
    }

    for (i = 0; i < Order; i++) {
	CagdSrfTransform(SrfFactors[i],
			 &Translate, CagdIChooseK(i, Order - 1));
    }

    for (i = 1; i < Order - 1; i++) {
    	CagdSrfFree(SrfFactorsSrf[i]);
    	CagdSrfFree(SrfFactors1_Srf[i]);
    }

    CagdSrfFree(Srf);

    return SrfFactors;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Divides given surface Srf into two tensor product surfaces along a	     M
* general curve DivCrv that splits Srf into two.                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:      To devide into two new tensor product surfaces along DivCurv.  M
*             Must be a Bezier surface.                                      M
*   DivCrv:   2D curve in the UV space of Srf.  Must split the UV domain of  M
*             Srf into two and must start/end at opposite boundaries.        M
*                That is, if DivCrv starts at UMin, it must end at UMax.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   A list of two surfaces resulting from subdividing Srf M
*             along splitting curve DivCrv.                                  M
*                                                                            *
* SEE ALSO:    BzrComposeSrfSrf                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   BzrSrfSubdivAtCurve                                                      M
*****************************************************************************/
CagdSrfStruct *BzrSrfSubdivAtCurve(const CagdSrfStruct *Srf,
				   const CagdCrvStruct *DivCrv)
{
    CagdRType *R, TMin, TMax, UMin, UMax, VMin, VMax;
    CagdUVType UVStart, UVEnd, UV1, UV2;
    CagdCrvStruct *Crv1, *Crv2;
    CagdSrfStruct *UVSrf1, *UVSrf2, *Srf1, *Srf2;

    if (!CAGD_IS_BEZIER_SRF(Srf)) {
        SYMB_FATAL_ERROR(SYMB_ERR_BSP_SRF_EXPECT);
        return NULL;
    }
    CagdCrvDomain(DivCrv, &TMin, &TMax);

    R = CagdCrvEval(DivCrv, TMin);
    CagdCoerceToE2(UVStart, &R, -1, DivCrv -> PType);
    R = CagdCrvEval(DivCrv, TMax);
    CagdCoerceToE2(UVEnd, &R, -1, DivCrv -> PType);

    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
    if ((IRIT_APX_EQ(UVStart[0], UMin) && IRIT_APX_EQ(UVEnd[0], UMax)) ||
	(IRIT_APX_EQ(UVStart[0], UMax) && IRIT_APX_EQ(UVEnd[0], UMin))) {
        /* Division curve is from UMin to UMax (or vice versa). */
        if (UVStart[0] > UVEnd[0]) {
	    IRIT_PT2D_SWAP(UVStart, UVEnd);
	    Crv1 = CagdCrvReverse(DivCrv);
	}
	else
	    Crv1 = CagdCrvCopy(DivCrv);

	UV1[0] = UMin;
	UV1[1] = VMin;
	UV2[0] = UMax;
	UV2[1] = VMin;
	Crv2 = CagdMergeUvUv(UV1, UV2);
	UVSrf1 = CagdRuledSrf(Crv2, Crv1, 2, 2);

	UV1[0] = UMin;
	UV1[1] = VMax;
	UV2[0] = UMax;
	UV2[1] = VMax;
	Crv2 = CagdMergeUvUv(UV1, UV2);
	UVSrf2 = CagdRuledSrf(Crv1, Crv2, 2, 2);

	CagdCrvFree(Crv1);
	CagdCrvFree(Crv2);
    }
    else if ((IRIT_APX_EQ(UVStart[1], VMin) && IRIT_APX_EQ(UVEnd[1], VMax)) ||
	     (IRIT_APX_EQ(UVStart[1], VMax) && IRIT_APX_EQ(UVEnd[1], VMin))) {
        /* Division curve is from VMin to VMax (or vice versa). */
        if (UVStart[1] > UVEnd[1]) {
	    IRIT_PT2D_SWAP(UVStart, UVEnd);
	    Crv1 = CagdCrvReverse(DivCrv);
	}
	else
	    Crv1 = CagdCrvCopy(DivCrv);

	UV1[0] = UMin;
	UV1[1] = VMin;
	UV2[0] = UMin;
	UV2[1] = VMax;
	Crv2 = CagdMergeUvUv(UV1, UV2);
	UVSrf1 = CagdRuledSrf(Crv2, Crv1, 2, 2);

	UV1[0] = UMax;
	UV1[1] = VMin;
	UV2[0] = UMax;
	UV2[1] = VMax;
	Crv2 = CagdMergeUvUv(UV1, UV2);
	UVSrf2 = CagdRuledSrf(Crv1, Crv2, 2, 2);

	CagdCrvFree(Crv1);
	CagdCrvFree(Crv2);
    }
    else {
        SYMB_FATAL_ERROR(SYMB_ERR_INVALID_CRV);
	return NULL;
    }

    Srf1 = BzrComposeSrfSrf(Srf, UVSrf1);
    Srf2 = BzrComposeSrfSrf(Srf, UVSrf2);

    CagdSrfFree(UVSrf1);
    CagdSrfFree(UVSrf2);

    Srf1 -> Pnext = Srf2;

    return Srf1;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes a surface patch that is a general rectangle in the domain of    M
* given surface.                                                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:     In which the four UVij corners resides.                         M
*   UV00, UV01, UV10, UV11:  The four corners of the patch to extract.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   A patch inside Srf that is a rectangle through UVij   M
*                      in the domain of Srf.  Only the four boundaries will  M
*                      be precisely reconstructed and the resulting surface  M
8		       will be contained in Srf only if Srf is planar.       M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeSrfCrv, SymbComposeSrfSrf                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbComposeSrfPatch                                                      M
*****************************************************************************/
CagdSrfStruct *SymbComposeSrfPatch(const CagdSrfStruct *Srf,
				   const CagdUVType UV00,
				   const CagdUVType UV01,
				   const CagdUVType UV10,
				   const CagdUVType UV11)
{
    CagdCrvStruct *Crv1E3, *Crv2E3, *Crv3E3, *Crv4E3,
        *Crv1UV = CagdMergeUvUv(UV00, UV01),             /* Left. */
        *Crv2UV = CagdMergeUvUv(UV10, UV11),             /* Right. */
        *Crv3UV = CagdMergeUvUv(UV00, UV10),             /* Top reversed. */
	*Crv4UV = CagdMergeUvUv(UV01, UV11);             /* Bottom. */
    CagdSrfStruct *BSumSrf;

    /* If both V values are zero, will create an E1 curve. Remedy that: */
    SYMB_VERIFY_E2_CRV(Crv1UV);
    SYMB_VERIFY_E2_CRV(Crv2UV);
    SYMB_VERIFY_E2_CRV(Crv3UV);
    SYMB_VERIFY_E2_CRV(Crv4UV);

    Crv1E3 = SymbComposeSrfCrv(Srf, Crv1UV),
    Crv2E3 = SymbComposeSrfCrv(Srf, Crv2UV),
    Crv3E3 = SymbComposeSrfCrv(Srf, Crv3UV),
    Crv4E3 = SymbComposeSrfCrv(Srf, Crv4UV);

    CagdCrvFree(Crv1UV);
    CagdCrvFree(Crv2UV);
    CagdCrvFree(Crv3UV);
    CagdCrvFree(Crv4UV);

    BSumSrf = CagdBoolSumSrf(Crv1E3, Crv2E3, Crv3E3, Crv4E3);

    CagdCrvFree(Crv1E3);
    CagdCrvFree(Crv2E3);
    CagdCrvFree(Crv3E3);
    CagdCrvFree(Crv4E3);

    return BSumSrf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tile an input object, in place, (UTimes x VTimes) in the given bivariate M
* surface. Computation is made precise, using composition operations.        M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:       The object to map through the bivariate.		     M
*               Can be a (list of) curve(s) or surface(s) only.		     M
*                 If PObj is formed out of surface(s), DeformSrf must be a   M
*               Bezier surface.						     M
*   DeformSrf:  The mapping/deformation function from R2 to R2.		     M
*   UTimes, VTimes:  Number of times to tile the object in each axis.        M
*   FitObj:     TRUE to rescale PObj tile to precisely fit the domain        M
*		                              (UTimes x VTimes),	     M
*               FALSE to assume PObj is in [0,1]^2 when fitting domain.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  (UTimes x VTimes) mapped and deformed objects.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeSrfCrv, SymbComposeSrfSrf, TrivComposeTileObjectInTVBzr       M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbComposeTileObjectInSrf                                               M
*****************************************************************************/
IPObjectStruct *SymbComposeTileObjectInSrf(const IPObjectStruct *PObj,
					   const CagdSrfStruct *DeformSrf,
					   IrtRType UTimes,
					   IrtRType VTimes,
					   IrtBType FitObj)
{
    int u, v;
    CagdRType UMin, UMax, VMin, VMax;
    IPObjectStruct *PNewObj, *PTile, *PTmp,
        *PRetVal = IPGenLISTObject(NULL);
    GMBBBboxStruct
        *BBox = GMBBComputeBboxObject(PObj);
    IrtHmgnMatType Mat1, Mat2;

    /* We expect either curves, or surfaces with a Bezier deformation only. */
    if (IP_IS_CRV_OBJ(PObj)) {
    }
    else if (IP_IS_SRF_OBJ(PObj) && !CAGD_IS_BEZIER_SRF(DeformSrf)) {
        SYMB_FATAL_ERROR(SYMB_ERR_BZR_SRF_EXPECT);
	return NULL;
    }
    else {
        SYMB_FATAL_ERROR(SYMB_ERR_WRONG_INPUT);
	return NULL;
    }

    CagdSrfDomain(DeformSrf, &UMin, &UMax, &VMin, &VMax);

    if (FitObj) {
        /* Map the input object to the proper tile size, in DeformSrf. */
        MatGenMatTrans(-BBox -> Min[0], -BBox -> Min[1], 0.0, Mat1);
	MatGenMatScale((UMax - UMin) /
		           (UTimes * (BBox -> Max[0] - BBox -> Min[0])),
		       (VMax - VMin) /
		           (VTimes * (BBox -> Max[1] - BBox -> Min[1])),
		       1.0, Mat2);
	MatMultTwo4by4(Mat1, Mat1, Mat2);
    }
    else {
        /* Scale [0,1]^3 to fit DeformSrf domain as many times as needed. */
        MatGenMatScale((UMax - UMin) / UTimes, (VMax - VMin) / VTimes,
		       1.0, Mat1);
    }

    PTile = GMTransformObject(PObj, Mat1);

#ifdef DEBUG_TEST_BBOX
    BBox = GMBBComputeBboxObject(PTile);
    fprintf(stderr, "BBOX = [%f  %f  %f] :: [%f  %f  %f]\n",
	    BBox -> Min[0], BBox -> Min[1], BBox -> Min[2],
	    BBox -> Max[0], BBox -> Max[1], BBox -> Max[2]);
#endif /* DEBUG_TEST_BBOX */

    for (v = 0; v < VTimes; v++) {
        for (u = 0; u < UTimes; u++) {
	    CagdCrvStruct *Crv, *TCrv, *RetCrvs;
	    CagdSrfStruct *Srf, *TSrf, *RetSrfs;

	    MatGenMatTrans(UMin + u * (UMax - UMin) / UTimes,
			   VMin + v * (VMax - VMin) / VTimes,
			   0.0, Mat1);
	    PNewObj = GMTransformObject(PTile, Mat1);
	    switch (PTile -> ObjType) {
	        case IP_OBJ_CURVE:
		    RetCrvs = NULL;
		    for (Crv = PNewObj -> U.Crvs;
			 Crv != NULL;
			 Crv = Crv -> Pnext) {
		        TCrv = SymbComposeSrfCrv(DeformSrf, Crv);
			IRIT_LIST_PUSH(TCrv, RetCrvs);
		    }
		    PTmp = IPGenCRVObject(RetCrvs);
		    break;
	        case IP_OBJ_SURFACE:
		    RetSrfs = NULL;
		    for (Srf = PNewObj -> U.Srfs;
			 Srf != NULL;
			 Srf = Srf -> Pnext) {
		        TSrf = SymbComposeSrfSrf(DeformSrf, Srf);
			IRIT_LIST_PUSH(TSrf, RetSrfs);
		    }
		    PTmp = IPGenSRFObject(RetSrfs);
		    break;
	        default:
		    assert(0);
		    return NULL;
	    }
	    IPFreeObject(PNewObj);
	    IPListObjectAppend(PRetVal, PTmp);
	}
    }

    IPFreeObject(PTile);

    return PRetVal;
}
