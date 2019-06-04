/******************************************************************************
* TrivEXTR.c - Extrusion operator out of a given surface and a direction vec. *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Mar. 2000.					      *
******************************************************************************/

#include "triv_loc.h"
#include "inc_irit/mvar_lib.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs an extrusion trivariate in the Vector direction for the given   M
* surface. Input surface can be either a Bspline or a Bezier surface and     M
* the resulting output trivariate will be of the same type.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:     To exturde in direction specified by Vec.                       M
*   Vec:     Direction as well as magnitude of extursion.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:   An extrusion trivariate volume with Orders of the      M
*                      original Srf order and 2 in the extrusion direction.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdExtrudeSrf, TrivExtrudeTV2                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivExtrudeTV, trivariate constructors                                   M
*****************************************************************************/
TrivTVStruct *TrivExtrudeTV(const CagdSrfStruct *Srf, const CagdVecStruct *Vec)
{
    TrivTVStruct *TV;
    int i, j,
	MaxCoord = CAGD_NUM_OF_PT_COORD(Srf -> PType),
	Len = Srf -> ULength * Srf -> VLength;
    CagdPointType
	PType = Srf -> PType;
    CagdBType
	IsNotRational = !CAGD_IS_RATIONAL_SRF(Srf);
    CagdRType **TVPoints;
    CagdRType
	* const *SrfPoints = Srf -> Points;
    CagdRType
	const *Dir = Vec -> Vec;

    switch (PType) {
	case CAGD_PT_P2_TYPE:
	    PType = CAGD_PT_P3_TYPE;
	    break;
	case CAGD_PT_E2_TYPE:
	    PType = CAGD_PT_E3_TYPE;
	    break;
	case CAGD_PT_P3_TYPE:
	case CAGD_PT_E3_TYPE:
	    break;
	default:
	    TRIV_FATAL_ERROR(TRIV_ERR_UNSUPPORT_PT);
	    break;
    }

    switch (Srf -> GType) {
	case CAGD_SBEZIER_TYPE:
	    TV = TrivBzrTVNew(Srf -> ULength, Srf -> VLength, 2, PType);
	    break;
	case CAGD_SBSPLINE_TYPE:
	    TV = TrivBspTVNew(Srf -> ULength, Srf -> VLength, 2,
			      Srf -> UOrder, Srf -> VOrder, 2,
			      PType);
	    CAGD_GEN_COPY(TV -> UKnotVector, Srf -> UKnotVector,
			  sizeof(CagdRType) * (TV -> ULength + TV -> UOrder));
	    CAGD_GEN_COPY(TV -> VKnotVector, Srf -> VKnotVector,
			  sizeof(CagdRType) * (TV -> VLength + TV -> VOrder));
	    TV -> WKnotVector[0] = TV -> WKnotVector[1] = 0.0;
	    TV -> WKnotVector[2] = TV -> WKnotVector[3] = 1.0;
	    break;
	case CAGD_SPOWER_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_POWER_NO_SUPPORT);
	    return NULL;
	default:
	    TRIV_FATAL_ERROR(TRIV_ERR_UNDEF_SRF);
	    return NULL;
    }

    /* Copy the control mesh - first layer is exactly the same as the        */
    /* surface while second one is the same as first one translated by Vec.  */
    TVPoints = TV -> Points;

    for (i = IsNotRational; i <= MaxCoord; i++)	       /* First depth layer. */
	CAGD_GEN_COPY(TVPoints[i], SrfPoints[i], sizeof(CagdRType) * Len);

    /* Make a copy of the Second layer do we can "work" on it. */
    for (i = IsNotRational; i <= MaxCoord; i++)	      /* Second depth layer. */
	CAGD_GEN_COPY(&TVPoints[i][Len], SrfPoints[i],
		      sizeof(CagdRType) * Len);

    /* If the surface has lesser dimension (i.e. was 2D), Add zeros. */
    for (i = MaxCoord + 1; i <= 3; i++)
	for (j = 0; j < Len * 2; j++)
	    TVPoints[i][j] = 0.0;

    for (i = 1; i <= 3; i++)		      /* Translate the second layer. */
	for (j = Len; j < Len * 2; j++)
	    TVPoints[i][j] += IsNotRational ? Dir[i - 1] :
					       Dir[i - 1] * TVPoints[0][j];

    TRIV_SET_GEOM_TYPE(TV, TRIV_GEOM_EXTRUSION);

    return TV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a full circular twisted/rotated extrusion volume in the +Z      M
* direction for the given profile surface. Input surface can be either a     M
* Bspline or a Bezier surface.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:       To twist and extrude in the +Z direction.                     M
*   Rational:  TRUE to construct a rational (and precise) twist, FALSE to    M
*              approximate using polynomials.				     M
*   ZPitch:    The +Z amount for full 360 degrees.  If zero, the result will M
*              be a planar (degenerated) surface.  A negative value will     M
*              reverse the twist.			                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:   A twisted extrusion trivariate.                        M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivExtrudeSrf							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivZTwistExtrudeSrf, trivariate constructors                            M
*****************************************************************************/
TrivTVStruct *TrivZTwistExtrudeSrf(const CagdSrfStruct *Srf,
				   CagdBType Rational,
				   CagdRType ZPitch)
{
    int i, j, k, SrfLen;
    CagdRType **Pts;
    CagdRType * const *SrfPts, * const *CircPts;
    CagdPointType PType;
    CagdCrvStruct *Circ;
    TrivTVStruct *ExtTV;

    if (CAGD_IS_RATIONAL_SRF(Srf))
        Rational = TRUE;
    PType = Rational ? CAGD_PT_P3_TYPE : CAGD_PT_E3_TYPE;       
    Circ = Rational ? BspCrvCreateUnitCircle() : BspCrvCreateUnitPCircle();
    if (IRIT_SIGN(ZPitch) < 0.0) {
        IrtHmgnMatType Mat;
	CagdCrvStruct *TCrv;

	/* Reverse the twist direction. */
	MatGenMatScale(1.0, -1.0, 1.0, Mat);
	TCrv = CagdCrvMatTransform(Circ, Mat);
	CagdCrvFree(Circ);
	Circ = TCrv;
    }

    ExtTV = TrivBspTVNew(Srf -> ULength, Srf -> VLength, Circ -> Length,
			 Srf -> UOrder, Srf -> VOrder, Circ -> Order, PType);

    Pts = ExtTV -> Points;
    CircPts = Circ -> Points;
    SrfPts = Srf -> Points;
    SrfLen = Srf -> ULength * Srf -> VLength;
    for (i = 0; i < SrfLen; i++) {
        CagdRType PtLen, CosTheta, SinTheta;
        CagdPType Pt;

	CagdCoerceToE3(Pt, SrfPts, i, Srf -> PType);
	PtLen = IRIT_PT2D_LENGTH(Pt);
	CosTheta = PtLen == 0.0 ? 0.0 : Pt[0] / PtLen;
	SinTheta = PtLen == 0.0 ? 0.0 : Pt[1] / PtLen;

        for (j = 0, k = i; j < Circ -> Length; j++, k += SrfLen) {
	    CagdRType
	        WCirc = CAGD_IS_RATIONAL_CRV(Circ) ? CircPts[0][j] : 1.0;

	    Pts[1][k] = (CircPts[1][j] * CosTheta - CircPts[2][j] * SinTheta)
		                                              * PtLen / WCirc;
	    Pts[2][k] = (CircPts[1][j] * SinTheta + CircPts[2][j] * CosTheta)
		                                              * PtLen / WCirc;
	    Pts[3][k] = IRIT_ABS(ZPitch) * j / (Circ -> Length - 1);

	    if (Rational) {
	        int l;
	        CagdRType
		    Wgt = CircPts[0][j];

		if (CAGD_IS_RATIONAL_CRV(Srf))
		    Wgt *= SrfPts[0][i];

		for (l = 1; l <= 3; l++)
		    Pts[l][k] *= Wgt;

		Pts[0][k] = Wgt;
	    }
	}
    }

    BspKnotCopy(ExtTV -> UKnotVector, Srf -> UKnotVector,
		Srf -> UOrder + Srf -> ULength);
    BspKnotCopy(ExtTV -> VKnotVector, Srf -> VKnotVector,
		Srf -> VOrder + Srf -> VLength);
    BspKnotCopy(ExtTV -> WKnotVector, Circ -> KnotVector,
		Circ -> Order + Circ -> Length);

    CagdCrvFree(Circ);

    return ExtTV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs an extrusion trivariate along Crv, of the given surface:      M
* TV(u, v, t) = Srf(u, v) + Crv(t).					     M
*   Input curve/surface can be either a Bspline or a Bezier surface and      M
* the resulting output trivariate will be of the same type.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:     To exturde along the curve Crv.                                 M
*   Crv:     Curve along which to move Srf.  If Crv is  a line, reduces to   M
*            TrivExtrudeTV.				                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:   A trivariate volume with Orders of the original        M
*                     Srf/Crv orders.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivExtrudeTV2, CagdExtrudeSrf, SymbAlgebraicSumSrf                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivExtrudeTV2, trivariate constructors                                  M
*****************************************************************************/
TrivTVStruct *TrivExtrudeTV2(const CagdSrfStruct *Srf,
			     const CagdCrvStruct *Crv)
{
    TrivTVStruct *TV;
    MvarMVStruct *MV,
        *MVCrv1 = MvarCrvToMV(Crv),
        *MVCrv2 = MvarPromoteMVToMV2(MVCrv1, 3, 2),
        *MVSrf1 = MvarSrfToMV(Srf),
        *MVSrf2 = MvarPromoteMVToMV2(MVSrf1, 3, 0);

    MvarMVFree(MVCrv1);
    MvarMVFree(MVSrf1);

    if (!MvarMakeMVsCompatible(&MVCrv2, &MVSrf2, TRUE, TRUE) ||
	(MV = MvarMVAdd(MVCrv2, MVSrf2)) == NULL) {
        MvarMVFree(MVCrv2);
        MvarMVFree(MVSrf2);
	TRIV_FATAL_ERROR(TRIV_ERR_TVS_INCOMPATIBLE);
	return NULL;
    }

    MvarMVFree(MVCrv2);
    MvarMVFree(MVSrf2);

    TV = MvarMVToTV(MV);
    MvarMVFree(MV);

    return TV;
}
