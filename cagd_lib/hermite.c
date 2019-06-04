/******************************************************************************
* Hermite.c - Hermite curve and surface constructors.            	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Apr. 95.					      *
******************************************************************************/

#include "cagd_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construct a cubic Bezier curve using the Hermite constraints - two       M
* positions and two tangents.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Pt1, Pt2:   Starting and end points of curve.                            M
*   Dir1, Dir2: Starting and end vectors of curve.                           M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   A cubic Bezier curve, satisfying the four constrants. M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdCubicHermiteCrv, Hermite                                             M
*****************************************************************************/
CagdCrvStruct *CagdCubicHermiteCrv(const CagdPType Pt1,
				   const CagdPType Pt2,
				   const CagdVType Dir1,
				   const CagdVType Dir2)
{
    int i;
    CagdCrvStruct
	*Crv = BzrCrvNew(4, CAGD_PT_E3_TYPE);
    CagdRType
	**Points = Crv -> Points;

    for (i = 0; i < 3; i++) {
	Points[i + 1][0] = Pt1[i];
	Points[i + 1][1] = Pt1[i] + Dir1[i] / 3.0;
	Points[i + 1][2] = Pt2[i] - Dir2[i] / 3.0;
	Points[i + 1][3] = Pt2[i];
    }

    return Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construct a cubic surface using the Hermite constraints - two            M
* positions and two tangents.  Other direction's degree depends on input.    M
*                                                                            *
* PARAMETERS:                                                                M
*   CPos1Crv, CPos2Crv:  Starting and end curves of surface.                 M
*   CDir1Crv, CDir2Crv:  Starting and end tangent fields surface.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   A cubic by something Bezier surface, satisfying the   M
*		       four constrants. The other something degree is the    M
*		       largest of the four given curves.		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdCubicHermiteSrf, Hermite                                             M
*****************************************************************************/
CagdSrfStruct *CagdCubicHermiteSrf(const CagdCrvStruct *CPos1Crv,
				   const CagdCrvStruct *CPos2Crv,
				   const CagdCrvStruct *CDir1Crv,
				   const CagdCrvStruct *CDir2Crv)
{
    int i, j, MaxAxis;
    CagdSrfStruct *Srf;
    CagdRType **Points;
    CagdCrvStruct *Pos1Crv, *Pos2Crv, *Dir1Crv, *Dir2Crv;

    Pos1Crv = CagdCrvCopy(CPos1Crv);
    Pos2Crv = CagdCrvCopy(CPos2Crv);
    Dir1Crv = CagdCrvCopy(CDir1Crv);
    Dir2Crv = CagdCrvCopy(CDir2Crv);

    if (!CagdMakeCrvsCompatible(&Pos1Crv, &Pos2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir1Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Pos2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir1Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir2Crv, TRUE, TRUE)) {
	CAGD_FATAL_ERROR(CAGD_ERR_CRV_FAIL_CMPT);
	CagdCrvFree(Pos1Crv);
	CagdCrvFree(Pos2Crv);
	CagdCrvFree(Dir1Crv);
	CagdCrvFree(Dir2Crv);
	return NULL;
    }

    if (CAGD_IS_BEZIER_CRV(Pos1Crv))
	Srf = BzrSrfNew(4, Pos1Crv -> Order, Pos1Crv -> PType);
    else {
	Srf = BspSrfNew(4, Pos1Crv -> Length,
			4, Pos1Crv -> Order, Pos1Crv -> PType);
	BspKnotUniformOpen(4, 4, Srf -> UKnotVector);
	CAGD_GEN_COPY(Srf -> VKnotVector, Pos1Crv -> KnotVector,
		      sizeof(CagdRType) *
		          (Pos1Crv -> Length + Pos1Crv -> Order +
			   (Pos1Crv -> Periodic ? Pos1Crv -> Order - 1 : 0)));
    }
    Points = Srf -> Points;

    MaxAxis = CAGD_NUM_OF_PT_COORD(Srf -> PType);

    for (j = 0; j < Pos1Crv -> Length; j++) {
	int Offset = 4 * j;

	for (i = !CAGD_IS_RATIONAL_SRF(Srf); i <= MaxAxis; i++) {
	    Points[i][Offset]     = Pos1Crv -> Points[i][j];
	    Points[i][Offset + 1] = Pos1Crv -> Points[i][j] +
					Dir1Crv -> Points[i][j] / 3.0;
	    Points[i][Offset + 2] = Pos2Crv -> Points[i][j] -
					Dir2Crv -> Points[i][j] / 3.0;
	    Points[i][Offset + 3] = Pos2Crv -> Points[i][j];
	}
    }

    CagdCrvFree(Pos1Crv);
    CagdCrvFree(Pos2Crv);
    CagdCrvFree(Dir1Crv);
    CagdCrvFree(Dir2Crv);

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construct a cubic surface using the Hermite constraints - two            M
* positions and two tangents.  Other direction's degree depends on input.    M
*                                                                            *
* PARAMETERS:                                                                M
*   CPos1Crv, CPos2Crv:    Starting and end curves of surface.               M
*   CDir1Crv, CDir2Crv:    Starting and end tangent fields surface.          M
*   C2Dir1Crv, C2Dir2Crv:  Starting and end 2nd derivative fields surface.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   A cubic by something Bezier surface, satisfying the   M
*		       four constrants. The other something degree is the    M
*		       largest of the four given curves.		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdQuinticHermiteSrf, Hermite                                           M
*****************************************************************************/
CagdSrfStruct *CagdQuinticHermiteSrf(const CagdCrvStruct *CPos1Crv,
				     const CagdCrvStruct *CPos2Crv,
				     const CagdCrvStruct *CDir1Crv,
				     const CagdCrvStruct *CDir2Crv,
				     const CagdCrvStruct *C2Dir1Crv,
				     const CagdCrvStruct *C2Dir2Crv)
{
    int i, j, MaxAxis;
    CagdSrfStruct *Srf;
    CagdRType **Points;
    CagdCrvStruct *Pos1Crv, *Pos2Crv, *Dir1Crv, *Dir2Crv,
        *Crvtr1Crv, *Crvtr2Crv;

    Pos1Crv = CagdCrvCopy(CPos1Crv);
    Pos2Crv = CagdCrvCopy(CPos2Crv);
    Dir1Crv = CagdCrvCopy(CDir1Crv);
    Dir2Crv = CagdCrvCopy(CDir2Crv);
    Crvtr1Crv = CagdCrvCopy(C2Dir1Crv);
    Crvtr2Crv = CagdCrvCopy(C2Dir2Crv);

    if (!CagdMakeCrvsCompatible(&Pos1Crv, &Pos2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir1Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Crvtr2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Crvtr2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Pos2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir1Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Dir2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Crvtr2Crv, TRUE, TRUE) ||
	!CagdMakeCrvsCompatible(&Pos1Crv, &Crvtr2Crv, TRUE, TRUE)) {
	CAGD_FATAL_ERROR(CAGD_ERR_CRV_FAIL_CMPT);
	CagdCrvFree(Pos1Crv);
	CagdCrvFree(Pos2Crv);
	CagdCrvFree(Dir1Crv);
	CagdCrvFree(Dir2Crv);
	CagdCrvFree(Crvtr1Crv);
	CagdCrvFree(Crvtr2Crv);
	return NULL;
    }

    if (CAGD_IS_BEZIER_CRV(Pos1Crv))
	Srf = BzrSrfNew(6, Pos1Crv -> Order, Pos1Crv -> PType);
    else {
	Srf = BspSrfNew(6, Pos1Crv -> Length,
			6, Pos1Crv -> Order, Pos1Crv -> PType);
	BspKnotUniformOpen(6, 6, Srf -> UKnotVector);
	CAGD_GEN_COPY(Srf -> VKnotVector, Pos1Crv -> KnotVector,
		      sizeof(CagdRType) *
		          (Pos1Crv -> Length + Pos1Crv -> Order +
			   (Pos1Crv -> Periodic ? Pos1Crv -> Order - 1 : 0)));
    }
    Points = Srf -> Points;

    MaxAxis = CAGD_NUM_OF_PT_COORD(Srf -> PType);

    for (j = 0; j < Pos1Crv -> Length; j++) {
	int Offset = 6 * j;

	for (i = !CAGD_IS_RATIONAL_SRF(Srf); i <= MaxAxis; i++) {
	    Points[i][Offset]     = Pos1Crv -> Points[i][j];
	    Points[i][Offset + 1] = Points[i][Offset] +
					Dir1Crv -> Points[i][j] / 5.0;
	    Points[i][Offset + 2] = -Points[i][Offset] +
	                             Points[i][Offset + 1] * 2.0 +
					Crvtr1Crv -> Points[i][j] / 20.0;

	    Points[i][Offset + 5] = Pos2Crv -> Points[i][j];
	    Points[i][Offset + 4] = Points[i][Offset + 5] -
					Dir2Crv -> Points[i][j] / 5.0;
	    Points[i][Offset + 3] = -Points[i][Offset + 5] +
	                             Points[i][Offset + 4] * 2.0 +
					Crvtr2Crv -> Points[i][j] / 20.0;
	}
    }

    CagdCrvFree(Pos1Crv);
    CagdCrvFree(Pos2Crv);
    CagdCrvFree(Dir1Crv);
    CagdCrvFree(Dir2Crv);
    CagdCrvFree(Crvtr1Crv);
    CagdCrvFree(Crvtr2Crv);

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Blends the two given surfaces along the u (first) parameter value.       M
*   Returned surface will start and Srf1(UMin) and terminate at Srf2(UMax).  M
*   Continuity is governed by the blending degree. 2 for C^0, 4 for C^1.     M
*   Odd degree will be rounded up to the next even degree. 		     M
*   Uses Hermite interpolation for the C^1 case.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf1:        First surface to blend.                                     M
*   Srf2:        Second surface to blend.                                    M
*   BlendDegree: Degree of the blending function in U. 2 for C^0, 4 for C^1, M
*                6 for C^2.						     M
*   TanScale:    If C^1, sets the tangency scale factor.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:    Resulting blended surface.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbBlendTwoSurfaces                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdBlendTwoSurfaces                                                     M
*****************************************************************************/
CagdSrfStruct *CagdBlendTwoSurfaces(const CagdSrfStruct *Srf1,
				    const CagdSrfStruct *Srf2,
				    int BlendDegree,
				    CagdRType TanScale)
{
    CagdVType TanScaleVec;
    CagdCrvStruct *Crv1, *Crv2, *DCrv1, *DCrv2, *D2Crv1, *D2Crv2;
    CagdSrfStruct *CpSrf1, *CpSrf2, *TSrf1, *TSrf2, *TSrf11, *TSrf22, *RetSrf;

    /* Make sure given surfaces span domain [0, 1]^2. */
    if (CAGD_IS_BEZIER_SRF(Srf1))
        CpSrf1 = CagdCnvrtBzr2BspSrf(Srf1);
    else
        CpSrf1 = CagdSrfCopy(Srf1);

    BspKnotAffineTrans2(CpSrf1 -> UKnotVector,
			CpSrf1 -> ULength + CpSrf1 -> UOrder, 0.0, 1.0);
    BspKnotAffineTrans2(CpSrf1 -> VKnotVector,
			CpSrf1 -> VLength + CpSrf1 -> VOrder, 0.0, 1.0);

    if (CAGD_IS_BEZIER_SRF(Srf2))
        CpSrf2 = CagdCnvrtBzr2BspSrf(Srf2);
    else
        CpSrf2 = CagdSrfCopy(Srf2);

    BspKnotAffineTrans2(CpSrf2 -> UKnotVector,
			CpSrf2 -> ULength + CpSrf2 -> UOrder, 0.0, 1.0);
    BspKnotAffineTrans2(CpSrf2 -> VKnotVector,
			CpSrf2 -> VLength + CpSrf2 -> VOrder, 0.0, 1.0);

    Crv1 = CagdCrvFromSrf(CpSrf1, 0.0, CAGD_CONST_V_DIR);
    Crv2 = CagdCrvFromSrf(CpSrf2, 1.0, CAGD_CONST_V_DIR);

    switch (BlendDegree) {
        case 1:
        case 2:
	    RetSrf = CagdRuledSrf(Crv1, Crv2, 2, 2);
	    break;
        case 3:
        case 4:
	    TSrf1 = CagdSrfDerive(CpSrf1, CAGD_CONST_V_DIR);
	    DCrv1 = CagdCrvFromSrf(TSrf1, 0.0, CAGD_CONST_V_DIR);
	    CagdSrfFree(TSrf1);

	    TSrf2 = CagdSrfDerive(CpSrf2, CAGD_CONST_V_DIR);
	    DCrv2 = CagdCrvFromSrf(TSrf2, 1.0, CAGD_CONST_V_DIR);
	    CagdSrfFree(TSrf2);

	    TanScaleVec[0] = TanScaleVec[1] = TanScaleVec[2] = TanScale;
	    CagdCrvScale(DCrv1, TanScaleVec);
	    CagdCrvScale(DCrv2, TanScaleVec);

	    RetSrf = CagdCubicHermiteSrf(Crv1, Crv2, DCrv1, DCrv2);
	    CagdCrvFree(DCrv1);
	    CagdCrvFree(DCrv2);
	    break;
        case 5:
        case 6:
	    TSrf1 = CagdSrfDerive(CpSrf1, CAGD_CONST_V_DIR);
	    DCrv1 = CagdCrvFromSrf(TSrf1, 0.0, CAGD_CONST_V_DIR);
	    TSrf11 = CagdSrfDerive(TSrf1, CAGD_CONST_V_DIR);
	    D2Crv1 = CagdCrvFromSrf(TSrf11, 0.0, CAGD_CONST_V_DIR);
	    CagdSrfFree(TSrf1);
	    CagdSrfFree(TSrf11);

	    TSrf2 = CagdSrfDerive(CpSrf2, CAGD_CONST_V_DIR);
	    DCrv2 = CagdCrvFromSrf(TSrf2, 1.0, CAGD_CONST_V_DIR);
	    TSrf22 = CagdSrfDerive(TSrf2, CAGD_CONST_V_DIR);
	    D2Crv2 = CagdCrvFromSrf(TSrf22, 1.0, CAGD_CONST_V_DIR);
	    CagdSrfFree(TSrf2);
	    CagdSrfFree(TSrf22);

	    TanScaleVec[0] = TanScaleVec[1] = TanScaleVec[2] = TanScale;
	    CagdCrvScale(DCrv1, TanScaleVec);
	    CagdCrvScale(DCrv2, TanScaleVec);

	    RetSrf = CagdQuinticHermiteSrf(Crv1, Crv2, DCrv1, DCrv2,
					   D2Crv1, D2Crv2);
	    CagdCrvFree(DCrv1);
	    CagdCrvFree(DCrv2);
	    CagdCrvFree(D2Crv1);
	    CagdCrvFree(D2Crv2);
	    break;
        default:
	    RetSrf = NULL;
    }

    CagdCrvFree(Crv1);
    CagdCrvFree(Crv2);

    CagdSrfFree(CpSrf1);
    CagdSrfFree(CpSrf2);

    return RetSrf;
}
