/******************************************************************************
* TrivTRev.c - Trivariate of revolution out of a given profile.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Jul. 10.					      *
******************************************************************************/

#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "triv_loc.h"

#define MAX_SREV_ANGLE_ERR	(IRIT_UEPS / 4.0)
#define MAX_SREV_ANGLE_ITER	100

IRIT_STATIC_DATA CagdRType 
    CircRationalKnotVector[12] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4 },
    CircPolynomialKnotVector[17] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2,
							3, 3, 3, 4, 4, 4, 4 };

IRIT_STATIC_DATA CagdRType
    PolyApproxRotAngles[] = {
	0,
	28.911200818417, /* arctan(4 (sqrt(2) - 1) / 3) */
	90 - 28.911200818417,
	90
    };

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a trivariate of revolution around the Z axis of given surface.  M
* Resulting trivariate will be a B-spline trivariate, while input may be     M
* either a B-spline or a Bezier surface.                                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To create trivariate of revolution around Z with.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  Trivariate of revolution.                               M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivTVOfRev2, TrivTVOfRevAxis, TrivTVOfRevPolynomialApprox		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivTVOfRev, trivariate of revolution, trivariate constructors           M
*****************************************************************************/
TrivTVStruct *TrivTVOfRev(const CagdSrfStruct *Srf)
{
    int i, j, i9,
        ULen = Srf -> ULength,
        VLen = Srf -> VLength,
        Len = ULen * VLen;
    CagdPointType
	PType = Srf -> PType;
    CagdRType **TVPoints,
	sin45 = sin(M_PI / 4.0);
    CagdRType * const *SrfPoints;
    CagdSrfStruct *CpSrf;
    TrivTVStruct
        *TV = TrivBspTVNew(9, ULen, VLen, 3, Srf -> UOrder, Srf -> VOrder,
			   CAGD_PT_P3_TYPE);

    /* Make sure the surface resides in R^3. */
    if ((CAGD_NUM_OF_PT_COORD(PType) < 3)) {
        PType = CAGD_IS_RATIONAL_PT(PType) ? CAGD_PT_P3_TYPE
					   : CAGD_PT_E3_TYPE;
	Srf = CpSrf = CagdCoerceSrfTo(Srf, PType, FALSE);
    }
    else
        CpSrf = NULL;

    SrfPoints = Srf -> Points;

    IRIT_GEN_COPY(TV -> UKnotVector, CircRationalKnotVector,
		  sizeof(CagdRType) * 12);

    switch (Srf -> GType) {
	case CAGD_SBEZIER_TYPE:
	    BspKnotUniformOpen(ULen, Srf -> UOrder, TV -> VKnotVector);
	    BspKnotUniformOpen(VLen, Srf -> VOrder, TV -> WKnotVector);
	    break;
	case CAGD_SBSPLINE_TYPE:
	    CAGD_GEN_COPY(TV -> VKnotVector, Srf -> UKnotVector,
			  sizeof(CagdRType) * (CAGD_SRF_UPT_LST_LEN(Srf) +
					                       Srf -> UOrder));
	    CAGD_GEN_COPY(TV -> WKnotVector, Srf -> VKnotVector,
			  sizeof(CagdRType) * (CAGD_SRF_VPT_LST_LEN(Srf) +
					                       Srf -> VOrder));
	    break;
	case CAGD_SPOWER_TYPE:
	    CAGD_FATAL_ERROR(CAGD_ERR_POWER_NO_SUPPORT);
	    return NULL;
	default:
	    CAGD_FATAL_ERROR(CAGD_ERR_UNDEF_SRF);
	    return NULL;
    }

    TVPoints = TV -> Points;

    /* For each control point in original surface - generate 9 points that   */
    /* Form a circle perpendicular to the Z axis.			     */
    for (i = i9 = 0; i < Len; i++, i9 += 9) {
	TVPoints[0][i9] = 1.0;
	switch (PType) {
	    case CAGD_PT_P3_TYPE:
		TVPoints[0][i9] = SrfPoints[0][i];
	    case CAGD_PT_E3_TYPE:
		TVPoints[1][i9] = SrfPoints[1][i];
		TVPoints[2][i9] = SrfPoints[2][i];
		TVPoints[3][i9] = SrfPoints[3][i];
		break;
	    default:
		CAGD_FATAL_ERROR(CAGD_ERR_UNSUPPORT_PT);
		break;
	}

	/* Last point is exactly the same as first one in circle - copy it.  */
	for (j = 0; j <= 3; j++)
	    TVPoints[j][i9 + 8] = TVPoints[j][i9];

	/* The Z components are identical in all circle, while the XY        */
	/* components are rotated 45 degrees at a time:			     */
	for (j = 1; j < 8; j++) {
	    TVPoints[0][i9 + j] = TVPoints[0][i9];
	    TVPoints[1][i9 + j] = TVPoints[1][i9 + j - 1] * sin45 -
				  TVPoints[2][i9 + j - 1] * sin45;
	    TVPoints[2][i9 + j] = TVPoints[1][i9 + j - 1] * sin45 +
				  TVPoints[2][i9 + j - 1] * sin45;
	    TVPoints[3][i9 + j] = TVPoints[3][i9];
	}

	/* And finally we need to compensate for the W's on every second pt. */
	for (j = 1; j < 8; j += 2) {
	    TVPoints[0][i9 + j] *= sin45;
	    TVPoints[3][i9 + j] *= sin45;
	}
    }

    TRIV_SET_GEOM_TYPE(TV, TRIV_GEOM_TV_OF_REV);

    if (CpSrf != NULL)
        CagdSrfFree(CpSrf);

    return TV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a trivariate of revolution around vector Axis of given profile  M
* surface. Resulting trivariate will be a B-spline trivariate, while input   M
* may be either a B-spline or a Bezier surface.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To create trivariate of revolution around Axis.              M
*   AxisPoint:  Of axis of rotation of Srf.				     M
*   AxisVector: Of axis of rotation of Srf.				     M
*   PolyApprox: TRUE to construct a polynomial approximation volume of       M
*               revolution, FALSE to create precise rational volume.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  Trivariate of revolution.                               M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivTVOfRev, TrivTVOfRev2, TrivTVOfRevPolynomialApprox	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivTVOfRevAxis, trivariate of revolution, trivariate constructors       M
*****************************************************************************/
TrivTVStruct *TrivTVOfRevAxis(const CagdSrfStruct *Srf,
			      const TrivVType AxisPoint,
			      const TrivVType AxisVector,
			      CagdBType PolyApprox)
{
    IRIT_STATIC_DATA IrtPtType
        Origin = { 0.0, 0.0, 0.0 };
    IrtPtType ClosestPoint;
    IrtVecType UnitAxis;
    IrtHmgnMatType Mat, Mat2, InvMat;
    CagdSrfStruct *TSrf;
    TrivTVStruct *TV;

    /* Find closest point on Axis to origin (and use that as translation). */
    GMPointFromPointLine(Origin, AxisPoint, AxisVector, ClosestPoint);

    IRIT_VEC_COPY(UnitAxis, AxisVector);
    IRIT_VEC_NORMALIZE(UnitAxis);

    GMGenMatrixZ2Dir(Mat, UnitAxis);
    MatGenMatTrans(ClosestPoint[0], ClosestPoint[1], ClosestPoint[2],
		   Mat2);
    MatMultTwo4by4(Mat, Mat, Mat2);
    MatInverseMatrix(Mat, InvMat);		    /* Compute the inverse. */

    TSrf = CagdSrfMatTransform(Srf, InvMat); /* Bring Srf to Z axis of rot. */
    if (PolyApprox)
        TV = TrivTVOfRevPolynomialApprox(TSrf);
    else
        TV = TrivTVOfRev(TSrf);     /* Create the trivariate of revolution. */
    CagdSrfFree(TSrf);

    TrivTVMatTransform(TV, Mat);          /* Bring TV of rev. back to Axis. */

    CAGD_SET_GEOM_TYPE(TV, TRIV_GEOM_TV_OF_REV);

    return TV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a trivariate of revolution around the Z axis of the given       M
* profile surface from StartAngle to EndAngle. Resulting trivariate will be  M
* a B-spline surface, while input may be either a B-spline or a Bezier       M
* surface.							             M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To create trivariate of revolution around Z with.            M
*   PolyApprox: TRUE for a polynomial approximation, FALSE for a precise     M
*               rational construction.					     M
*   StartAngle: Starting Angle to consider rotating Srf from, in degrees.    M
*   EndAngle:   Terminating Angle to consider rotating Srf from, in degrees. M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  Trivariate of revolution.                               M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivTVOfRev, TrivTVOfRevAxis, TrivTVOfRevPolynomialApprox		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivTVOfRev2, surface of revolution, surface constructors                M
*****************************************************************************/
TrivTVStruct *TrivTVOfRev2(const CagdSrfStruct *Srf,
			   CagdBType PolyApprox,
			   CagdRType StartAngle,
			   CagdRType EndAngle)
{
    int j = 0;
    CagdRType StartParam, EndParam, *R, TMin, TMax, TMid, A;
    CagdCrvStruct
	*Circ = PolyApprox ? BspCrvCreateUnitPCircle()
                           : BspCrvCreateUnitCircle();
    TrivTVStruct *TVTmp, *TV;
    CagdMType ZRotMat;

    /* Find parameter values for starting and ending angles on a circle. */
    if (StartAngle > EndAngle)
	IRIT_SWAP(IrtRType, StartAngle, EndAngle)

    StartAngle = IRIT_DEG2RAD(StartAngle);
    MatGenMatRotZ1(StartAngle, ZRotMat);

    EndAngle = IRIT_DEG2RAD(EndAngle) - StartAngle;

    CagdCrvDomain(Circ, &TMin, &TMax);
    StartParam = TMin;
    do {
	/* Compute angle of mid of domain, between zero and 2Pi. */
	TMid = (TMin + TMax) * 0.5;
	R = CagdCrvEval(Circ, TMid);
	A = atan2(R[2], R[1]);
	if (A < 0.0)
	    A += M_PI_MUL_2;

	if (A > EndAngle)
	    TMax = TMid;
	else
	    TMin = TMid;
    }
    while (TMax - TMin > MAX_SREV_ANGLE_ERR && j++ < MAX_SREV_ANGLE_ITER);
    EndParam = (TMin + TMax) * 0.5;

    CagdCrvFree(Circ);

    TVTmp = PolyApprox ? TrivTVOfRevPolynomialApprox(Srf)
		       : TrivTVOfRev(Srf);
    TV = TrivTVRegionFromTV(TVTmp, StartParam, EndParam, CAGD_CONST_U_DIR);
    TrivTVFree(TVTmp);

    /* Rotate so it starts at StartAngle. */
    TrivTVMatTransform(TV, ZRotMat);

    TRIV_SET_GEOM_TYPE(TV, TRIV_GEOM_TV_OF_REV);

    return TV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a trivariate of revolution around the Z axis of given profile   M
* surface. Resulting trivariate will be a B-spline trivariate, while input   M
* may be either a B-spline or a Bezier surface.                              M
*   Resulting trivariate will be a polynomial B-spline trivariate,           M
* approximating a trivariate of revolution using a polynomial circle approx. M
* (See Faux & Pratt "Computational Geometry for Design and Manufacturing").  M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To approximate a trivariate of revolution around Z with.     M
*		Srf is assumed planar in a plane holding the Z axis.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  Trivariate of revolution approximation.                 M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivTVOfRev, TrivTVOfRev2, TrivTVOfRevAxis			             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivTVOfRevPolynomialApprox, trivariate of revolution, trivariate        M
*   constructors						             M
*****************************************************************************/
TrivTVStruct *TrivTVOfRevPolynomialApprox(const CagdSrfStruct *Srf)
{
    int i, j, i13,
        ULen = Srf -> ULength,
        VLen = Srf -> VLength,
        Len = ULen * VLen;
    CagdPointType
	PType = Srf -> PType;
    CagdRType **TVPoints;
    CagdRType * const *SrfPoints;
    CagdSrfStruct *CpSrf;
    TrivTVStruct
        *TV = TrivBspTVNew(13, ULen, VLen, 4, Srf -> UOrder, Srf -> VOrder,
			   CAGD_PT_E3_TYPE);

    if (CAGD_IS_RATIONAL_CRV(Srf)) {
	CAGD_FATAL_ERROR(CAGD_ERR_POLYNOMIAL_EXPECTED);
	return NULL;
    }

    /* Make sure the surface resides in R^3. */
    if (CAGD_NUM_OF_PT_COORD(PType) < 3)
	Srf = CpSrf = CagdCoerceSrfTo(Srf, PType = CAGD_PT_E3_TYPE, FALSE);
    else
        CpSrf = NULL;

    SrfPoints = Srf -> Points;

    IRIT_GEN_COPY(TV -> UKnotVector, CircPolynomialKnotVector,
		  sizeof(CagdRType) * 17);

    switch (Srf -> GType) {
	case CAGD_SBEZIER_TYPE:
	    BspKnotUniformOpen(ULen, Srf -> UOrder, TV -> VKnotVector);
	    BspKnotUniformOpen(VLen, Srf -> VOrder, TV -> WKnotVector);
	    break;
	case CAGD_SBSPLINE_TYPE:
	    CAGD_GEN_COPY(TV -> VKnotVector, Srf -> UKnotVector,
			  sizeof(CagdRType) * (CAGD_SRF_UPT_LST_LEN(Srf) +
					                       Srf -> UOrder));
	    CAGD_GEN_COPY(TV -> WKnotVector, Srf -> VKnotVector,
			  sizeof(CagdRType) * (CAGD_SRF_VPT_LST_LEN(Srf) +
					                       Srf -> VOrder));
	    break;
	case CAGD_SPOWER_TYPE:
	    CAGD_FATAL_ERROR(CAGD_ERR_POWER_NO_SUPPORT);
	    return NULL;
	default:
	    CAGD_FATAL_ERROR(CAGD_ERR_UNDEF_SRF);
	    return NULL;
    }

    TVPoints = TV -> Points;

    /* For each control point in original surface - generate 13 points that  */
    /* form a circle approximation perpendicular to the Z axis.		     */
    for (i = i13 = 0; i < Len; i++, i13 += 13) {
	int Quad;

	switch (PType) {
	    case CAGD_PT_E3_TYPE:
		TVPoints[1][i13] = SrfPoints[1][i];
		TVPoints[2][i13] = SrfPoints[2][i];
		TVPoints[3][i13] = SrfPoints[3][i];
		break;
	    default:
		CAGD_FATAL_ERROR(CAGD_ERR_UNSUPPORT_PT);
		break;
	}

	/* Last point is exactly the same as first one in circle - copy it.  */
	for (j = 1; j <= 3; j++)
	    TVPoints[j][i13 + 12] = TVPoints[j][i13];

	/* The Z components are identical in all circle, while the XY        */
	/* components are functions of PolyApproxRotAngles:		     */
	for (Quad = 0, j = 1; j < 12; j++) {
	    CagdRType Scl, Angle, CosAngle, SinAngle;

	    if (j % 3 == 0) {
		Quad++;
		Scl = 1.0;
	    }
	    else
	        Scl = 1.0 / cos(IRIT_DEG2RAD(PolyApproxRotAngles[1]));

	    Angle = Quad * 90 + PolyApproxRotAngles[j % 3];
	    CosAngle = cos(IRIT_DEG2RAD(Angle)) * Scl;
	    SinAngle = sin(IRIT_DEG2RAD(Angle)) * Scl;

	    TVPoints[1][i13 + j] = TVPoints[1][i13] * CosAngle -
	                           TVPoints[2][i13] * SinAngle;
	    TVPoints[2][i13 + j] = TVPoints[1][i13] * SinAngle +
	                           TVPoints[2][i13] * CosAngle;
	    TVPoints[3][i13 + j] = TVPoints[3][i13];
	}
    }

    CAGD_SET_GEOM_TYPE(TV, TRIV_GEOM_TV_OF_REV);

    if (CpSrf != NULL)
	CagdSrfFree(CpSrf);

    return TV;
}
