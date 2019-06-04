/******************************************************************************
* Compost3.c - crv-tv and srf-tv composition computation (symbolic).	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Sep 2013.					      *
******************************************************************************/

#include <string.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/symb_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "triv_loc.h"

static CagdSrfStruct *TrivComposeTVSrfAux(const TrivTVStruct *TV,
					  const CagdSrfStruct *Srf);
static CagdCrvStruct *TrivComposeTVCrvAux(const TrivTVStruct *TV,
					  const CagdCrvStruct *Crv);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tile an input object, in place, (UTimes x VTimes x WTimes) in the given  M
* trivariate. Computation is made precise, using composition operations.     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:       The object to map through the trivariate.                    M
*               Can be a (list of) curve(s) or surface(s) only.		     M
*   DeformTV:   The mapping/deformation function from R3 to R3.		     M
*   UTimes, VTimes, WTimes:  Number of times to tile the object in each      M
*               axis.						             M
*   FitObj:     TRUE to rescale PObj tile to precisely fit the domain        M
*		                              (UTimes x VTimes x WTimes),    M
*               FALSE to assume PObj is in [0,1]^3 when fitting domain.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  (UTimes x VTimes x WTimes) mapped and deformed        M
*               objects.					             M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivFFDTileObjectInTV, TrivFFDCtlMeshUsingTV, TrivFFDObjectTV            M
*   TrivComposeTileObjectInTVBzr					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivComposeTileObjectInTV                                                M
*****************************************************************************/
IPObjectStruct *TrivComposeTileObjectInTV(const IPObjectStruct *PObj,
					  const TrivTVStruct *DeformTV,
					  IrtRType UTimes,
					  IrtRType VTimes,
					  IrtRType WTimes,
					  IrtBType FitObj)
{
    CagdRType t;
    TrivTVStruct *DeformTVA, *DeformTVB;
    IPObjectStruct *CmpsDeformTV, *CmpsDeformTVA, *CmpsDeformTVB;

    if (DeformTV -> ULength > DeformTV -> UOrder) {
	/* DeformTV is not a Bezier in U. Subdivide, compute for each      */
	/* segment and merge back.				           */

	t = DeformTV -> UKnotVector[(DeformTV -> UOrder +
				                  DeformTV -> ULength) >> 1];

	DeformTVA = TrivTVSubdivAtParam(DeformTV, t, TRIV_CONST_U_DIR);
	DeformTVB = DeformTVA -> Pnext;
	DeformTVA -> Pnext = NULL;

	CmpsDeformTVA = TrivComposeTileObjectInTV(PObj, DeformTVA,
					     UTimes, VTimes, WTimes, FitObj);
	CmpsDeformTVB = TrivComposeTileObjectInTV(PObj, DeformTVB,
					     UTimes, VTimes, WTimes, FitObj);
	TrivTVFree(DeformTVA);
	TrivTVFree(DeformTVB);

	CmpsDeformTV = IPAppendListObjects(CmpsDeformTVA, CmpsDeformTVB);
	IPFreeObject(CmpsDeformTVA); /* Will not free sub element due to   */				     
	IPFreeObject(CmpsDeformTVB); /* ref. count.		           */
    }
    else if (DeformTV -> VLength > DeformTV -> VOrder) {
	/* DeformTV is not a Bezier in U. Subdivide, compute for each      */
	/* segment and merge back.				           */
	t = DeformTV -> VKnotVector[(DeformTV -> VOrder +
				                  DeformTV -> VLength) >> 1];

	DeformTVA = TrivTVSubdivAtParam(DeformTV, t, TRIV_CONST_V_DIR);
	DeformTVB = DeformTVA -> Pnext;
	DeformTVA -> Pnext = NULL;

	CmpsDeformTVA = TrivComposeTileObjectInTV(PObj, DeformTVA,
					    UTimes, VTimes, WTimes, FitObj);
	CmpsDeformTVB = TrivComposeTileObjectInTV(PObj, DeformTVB,
					    UTimes, VTimes, WTimes, FitObj);
	TrivTVFree(DeformTVA);
	TrivTVFree(DeformTVB);

	CmpsDeformTV = IPAppendListObjects(CmpsDeformTVA, CmpsDeformTVB);
	IPFreeObject(CmpsDeformTVA); /* Will not free sub element due to   */				     
	IPFreeObject(CmpsDeformTVB); /* ref. count.		           */
    }
    else if (DeformTV -> WLength > DeformTV -> WOrder) {
	/* DeformTV is not a Bezier in W. Subdivide, compute for each      */
	/* segment and merge back.				           */
	t = DeformTV -> WKnotVector[(DeformTV -> WOrder +
				                  DeformTV -> WLength) >> 1];

	DeformTVA = TrivTVSubdivAtParam(DeformTV, t, TRIV_CONST_W_DIR);
	DeformTVB = DeformTVA -> Pnext;
	DeformTVA -> Pnext = NULL;

	CmpsDeformTVA = TrivComposeTileObjectInTV(PObj, DeformTVA,
					    UTimes, VTimes, WTimes, FitObj);
	CmpsDeformTVB = TrivComposeTileObjectInTV(PObj, DeformTVB,
					    UTimes, VTimes, WTimes, FitObj);
	TrivTVFree(DeformTVA);
	TrivTVFree(DeformTVB);

	CmpsDeformTV = IPAppendListObjects(CmpsDeformTVA, CmpsDeformTVB);
	IPFreeObject(CmpsDeformTVA); /* Will not free sub element due to   */				     
	IPFreeObject(CmpsDeformTVB); /* ref. count.		           */
    }
    else {
	TrivTVStruct 
	    *DeformTVNew = NULL;

	if (!TRIV_IS_BEZIER_TV(DeformTV))
	    DeformTV = DeformTVNew = TrivCnvrtBsp2BzrTV(DeformTV);

	if (IP_IS_OLST_OBJ(PObj)) {
	    int i;
	    IPObjectStruct *PTmp, *PTmp2;

	    CmpsDeformTV = IPGenLISTObject(NULL);
	    for (i = 0; (PTmp = IPListObjectGet(PObj, i)) != NULL; i++) {
		PTmp2 = TrivComposeTileObjectInTVBzr(PTmp, DeformTV,
					     UTimes, VTimes, WTimes, FitObj);
		IPListObjectInsert(CmpsDeformTV, i, PTmp2);
	    }
	    IPListObjectInsert(CmpsDeformTV, i, NULL);
	}
	else {
            /* DeformTV is a Bezier surface segment - compute its composition. */
	    CmpsDeformTV = TrivComposeTileObjectInTVBzr(PObj, DeformTV,
					     UTimes, VTimes, WTimes, FitObj);
	}

	if (DeformTVNew)
	   TrivTVFree(DeformTVNew);
    }

    return CmpsDeformTV;    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tile an input object, in place, (UTimes x VTimes x WTimes) in the given  M
* Bezier trivariate. Computation is made precise, using composition          M
* operations.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:       The object to map through the trivariate.	             M
*               Can be a (list of) curve(s) or surface(s) only.		     M
*   DeformTV:   The mapping/deformation Bezier function from R3 to R3.	     M
*   UTimes, VTimes, WTimes:  Number of times to tile the object in each      M
*               axis.						             M
*   FitObj:     TRUE to rescale PObj tile to precisely fit the domain        M
*		                              (UTimes x VTimes x WTimes),    M
*               FALSE to assume PObj is in [0,1]^3 when fitting domain.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  (UTimes x VTimes x WTimes) mapped and deformed        M
*               objects.					             M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivFFDTileObjectInTV, TrivFFDCtlMeshUsingTV, TrivFFDObjectTV            M
*   TrivComposeTileObjectInTV						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivComposeTileObjectInTVBzr                                             M
*****************************************************************************/
IPObjectStruct *TrivComposeTileObjectInTVBzr(const IPObjectStruct *PObj,
					     const TrivTVStruct *DeformTV,
					     IrtRType UTimes,
					     IrtRType VTimes,
					     IrtRType WTimes,
					     IrtBType FitObj)
{
    int u, v, w;
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax;
    IPObjectStruct *PNewObj, *PTile, *PTmp,
        *PRetVal = IPGenLISTObject(NULL);
    GMBBBboxStruct
        *BBox = GMBBComputeBboxObject(PObj);
    IrtHmgnMatType Mat1, Mat2;

    if (!IP_IS_CRV_OBJ(PObj) &&
	!IP_IS_SRF_OBJ(PObj) &&
	!IP_IS_TRIMSRF_OBJ(PObj)) {
        TRIV_FATAL_ERROR(TRIV_ERR_CRV_OR_SRF_EXPECTED);
	return NULL;
    }

    TrivTVDomain(DeformTV, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);

    if (FitObj) {
        /* Map the input object to the proper tile size, in DeformTV. */
        MatGenMatTrans(-BBox -> Min[0],
		       -BBox -> Min[1],
		       -BBox -> Min[2], Mat1);
	MatGenMatScale((UMax - UMin) /
		           (UTimes * (BBox -> Max[0] - BBox -> Min[0])),
		       (VMax - VMin) /
		           (VTimes * (BBox -> Max[1] - BBox -> Min[1])),
		       (WMax - WMin) /
		           (WTimes * (BBox -> Max[2] - BBox -> Min[2])),
		       Mat2);
	MatMultTwo4by4(Mat1, Mat1, Mat2);
    }
    else {
        /* Scale [0,1]^3 to fit DeformTV domain as many times as needed. */
        MatGenMatScale((UMax - UMin) / UTimes,
		       (VMax - VMin) / VTimes,
		       (WMax - WMin) / WTimes,
		       Mat1);
    }

    PTile = GMTransformObject(PObj, Mat1);

#ifdef DEBUG_TEST_BBOX
    BBox = GMBBComputeBboxObject(PTile);
    fprintf(stderr, "BBOX = [%f  %f  %f] :: [%f  %f  %f]\n",
	    BBox -> Min[0], BBox -> Min[1], BBox -> Min[2],
	    BBox -> Max[0], BBox -> Max[1], BBox -> Max[2]);
#endif /* DEBUG_TEST_BBOX */

    for (w = 0; w < WTimes; w++) {
        for (v = 0; v < VTimes; v++) {
	    for (u = 0; u < UTimes; u++) {
	        CagdCrvStruct *Crv, *TCrv, *RetCrvs;
		CagdSrfStruct *Srf, *TSrf, *RetSrfs;
		TrimSrfStruct *TrimSrf, *TTrimSrf, *RetTrimSrfs;

		MatGenMatTrans(UMin + u * (UMax - UMin) / UTimes,
			       VMin + v * (VMax - VMin) / VTimes,
			       WMin + w * (WMax - WMin) / WTimes,
			       Mat1);
		PNewObj = GMTransformObject(PTile, Mat1);
		switch (PTile -> ObjType) {
  		    case IP_OBJ_CURVE:
		        RetCrvs = NULL;
		        for (Crv = PNewObj -> U.Crvs;
			     Crv != NULL;
			     Crv = Crv -> Pnext) {
			    TCrv = TrivComposeTVCrv(DeformTV, Crv);
			    IRIT_LIST_PUSH(TCrv, RetCrvs);
			}
			PTmp = IPGenCRVObject(RetCrvs);
		        break;
  		    case IP_OBJ_SURFACE:
		        RetSrfs = NULL;
		        for (Srf = PNewObj -> U.Srfs;
			     Srf != NULL;
			     Srf = Srf -> Pnext) {
			    TSrf = TrivComposeTVSrf(DeformTV, Srf);
			    IRIT_LIST_PUSH(TSrf, RetSrfs);
			}
			PTmp = IPGenSRFObject(RetSrfs);
		        break;
  		    case IP_OBJ_TRIMSRF:
		        RetTrimSrfs = NULL;
		        for (TrimSrf = PNewObj -> U.TrimSrfs;
			     TrimSrf != NULL;
			     TrimSrf = TrimSrf -> Pnext) {
			    TSrf = TrivComposeTVSrf(DeformTV, TrimSrf -> Srf);
			    TTrimSrf = TrimSrfNew(
				     TSrf,
				     TrimCrvCopyList(TrimSrf -> TrimCrvList),
				     TRUE);
			    IRIT_LIST_PUSH(TTrimSrf, RetTrimSrfs);
			}
			PTmp = IPGenTRIMSRFObject(RetTrimSrfs);
		        break;
  		    default:
		        assert(0);
			return NULL;
		}
		IPFreeObject(PNewObj);
		IPListObjectAppend(PRetVal, PTmp);
	    }
	}
    }

    IPFreeObject(PTile);

    return PRetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given curve Crv and trivar TV, computes the composition TV(Crv).           M
*   Crv must be a three dimensional curve completely contained in the        M
* parametric domain of TV.                                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   TV, Crv:   The trivar and curve to compose. TV must be a Bezier.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    The resulting composition.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeCrvCrv, SymbComposeSrfCrv, SymbComposeSrfSrf,                 M
*   TrivComposeTVSrf							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivComposeTVCrv, composition                                            M
*****************************************************************************/
CagdCrvStruct *TrivComposeTVCrv(const TrivTVStruct *TV,
				const CagdCrvStruct *Crv)
{
    switch (TV -> GType) {
	case TRIV_TVBEZIER_TYPE:
	    break;
	case TRIV_TVBSPLINE_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_BZR_TV_EXPECT);
	    break;
	case TRIV_TVPOWER_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_POWER_NO_SUPPORT);
	    break;
	default:
	    TRIV_FATAL_ERROR(TRIV_ERR_UNDEF_SRF);
	    break;
    }

    switch (Crv -> GType) {
	case CAGD_CBEZIER_TYPE:
	    break;
	case CAGD_CBSPLINE_TYPE:
	    break;
	case CAGD_CPOWER_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_POWER_NO_SUPPORT);
	    break;
	default:
	    TRIV_FATAL_ERROR(TRIV_ERR_UNDEF_CRV);
	    break;
    }

    return TrivComposeTVCrvAux(TV, Crv);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Aux. function. Subdivides Crv until it is a Bezier curve, invokes the      *
* Bezier composition code on each, and merges them back to complete curve.   *
*   At this point, the curve can be either Bezier or B-spline only.          *
*   Crv is assumed to have open end condition.	   		             *
*                                                                            *
* PARAMETERS:                                                                *
*   TV, Crv:   The trivar and curve to compose. TV must be a Bezier.         *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:  The composed curve.                                    *
*****************************************************************************/
static CagdCrvStruct *TrivComposeTVCrvAux(const TrivTVStruct *TV,
					  const CagdCrvStruct *Crv)
{
    CagdRType t;
    CagdCrvStruct *CmpsCrv;

    if (Crv -> Length > Crv -> Order) {
	/* Crv is not a Bezier segment. Subdivide, compute for each        */
	/* segment and merge back.				           */
	CagdCrvStruct *CrvA, *CrvB, *CmpsCrvA, *CmpsCrvB;

	t = Crv -> KnotVector[(Crv -> Order + Crv -> Length) >> 1];

	CrvA = CagdCrvSubdivAtParam(Crv, t);
	CrvB = CrvA -> Pnext;
	CrvA -> Pnext = NULL;

	CmpsCrvA = TrivComposeTVCrvAux(TV, CrvA);
	CmpsCrvB = TrivComposeTVCrvAux(TV, CrvB);
	CagdCrvFree(CrvA);
	CagdCrvFree(CrvB);

	CmpsCrv = CagdMergeCrvCrv(CmpsCrvA, CmpsCrvB, FALSE);
	CagdCrvFree(CmpsCrvA);
	CagdCrvFree(CmpsCrvB);
    }
    else {
    	/* Crv is a Bezier curve segment - compute its composition. */
	if (!CAGD_IS_BEZIER_CRV(Crv)) {
	    CagdCrvStruct 
		*TCrv = CagdCnvrtBsp2BzrCrv(Crv);

	    CmpsCrv = TrivBzrComposeTVCrv(TV, TCrv);
	    CagdCrvFree(TCrv);
	}
	else
    	    CmpsCrv = TrivBzrComposeTVCrv(TV, Crv);
    }

    return CmpsCrv;    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given curve Crv and Bezier trivariate TV, computes their composition     M
* TV(Crv(t)),    Crv(t) = (u(t), v(t), w(t)).				     M
*   Crv must be a three dimensional curve completely contained in the        M
* parametric domain of TV.                                                   M
*   Compute the compositions by the products of:			     M
*									     M
* TV(u, v, w) = TV(u(t), v(t), w(t))					     V
*									     V
*            n   m   l							     V
*         = sum sum sum Pijk Bi(u(t)) Bj(v(t)) Bk(w(t))			     V
*           i=0 j=0 k=0							     V
*				   					     V
*            n   m   l	      n	       i	   n-i			     V
*         = sum sum sum Pijk ( ) (u(t))  (1 - u(t))             	     V
*           i=0 j=0 k=0	      i						     V
*                   	         m	  j	      m-j		     V
*                               ( ) (v(t))  (1 - v(t))  	             V
*                  	         j					     V
*                   	            l	     k	         l-k		     V
*                                  ( ) (w(t))  (1 - w(t))       	     V
*                  	            k					     V
*                                                                            *
* PARAMETERS:                                                                M
*   TV, Crv:   The trivar and curve to compose. TV must be Bezier.           M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    The resulting composition.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeTVCrv, BzrComposeCrvCrv, TrivBzrComposeTVSrf		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivBzrComposeTVCrv, composition                                         M
*****************************************************************************/
CagdCrvStruct *TrivBzrComposeTVCrv(const TrivTVStruct *TV,
				   const CagdCrvStruct *Crv)
{
    IRIT_STATIC_DATA CagdCrvStruct
	**CrvUFactors = NULL,
        **CrvVFactors = NULL,
	**CrvWFactors = NULL;
    CagdBType
        IsRationalCrv = CAGD_IS_RATIONAL_CRV(Crv),
        IsRationalTV = TRIV_IS_RATIONAL_TV(TV),
        IsRational = IsRationalCrv | IsRationalTV;
    int i, j, k, l,
        CmpsLength = 0,
	UOrder = TV -> UOrder,
	VOrder = TV -> VOrder,
	WOrder = TV -> WOrder,
	MaxCoord = CAGD_NUM_OF_PT_COORD(TV -> PType);
    CagdRType * const
        *TVPoints = TV -> Points;
    CagdCrvStruct *CrvUVW,
    	*CmpsCrv = NULL;

    if (CAGD_NUM_OF_PT_COORD(Crv -> PType) < 3 || MaxCoord != 3) {
	TRIV_FATAL_ERROR(TRIV_ERR_UNSUPPORT_PT);
	return NULL;
    }

    CrvUVW = CagdCoerceCrvTo(Crv,
			     IsRationalCrv ? CAGD_PT_P1_TYPE : CAGD_PT_E1_TYPE,
			     FALSE);

    if (CrvUFactors != NULL) {
        for (i = 0; CrvUFactors[i] != NULL; i++)
	    CagdCrvFree(CrvUFactors[i]);
	IritFree(CrvUFactors);
    }
    CrvUFactors = SymbComputeCurvePowers(CrvUVW, UOrder);

    CAGD_GEN_COPY(CrvUVW -> Points[1], Crv -> Points[2],
		  sizeof(CagdRType) * Crv -> Length);

    if (CrvVFactors != NULL) {
        for (i = 0; CrvVFactors[i] != NULL; i++)
	    CagdCrvFree(CrvVFactors[i]);
	IritFree(CrvVFactors);
    }
    CrvVFactors = SymbComputeCurvePowers(CrvUVW, VOrder);

    CAGD_GEN_COPY(CrvUVW -> Points[1], Crv -> Points[3],
		  sizeof(CagdRType) * Crv -> Length);

    if (CrvWFactors != NULL) {
        for (i = 0; CrvWFactors[i] != NULL; i++)
	    CagdCrvFree(CrvWFactors[i]);
	IritFree(CrvWFactors);
    }
    CrvWFactors = SymbComputeCurvePowers(CrvUVW, WOrder);

    CagdCrvFree(CrvUVW);

    /* The main (triple) loop (Compositions computation): */
    for (k = 0; k < WOrder; k++) {
        for (j = 0; j < VOrder; j++) {
	    for (i = 0; i < UOrder; i++) {
	        int p,
		    Idx = TRIV_MESH_UVW(TV, i, j, k);
	        CagdCrvStruct
		    *TCrv1 = SymbCrvMult(CrvUFactors[i], CrvVFactors[j]),
		    *TCrv2 = SymbCrvMult(CrvWFactors[k], TCrv1);
		CagdRType
		    *TCrvPoints = TCrv2 -> Points[1];

		CagdCrvFree(TCrv1);

		for (p = !IsRationalTV; p <= MaxCoord; p++) {
		    CagdRType *CmpsPoints,
		        SPt = TVPoints[p][Idx];

		    if (i == 0 && j == 0 && k == 0) {
		        if (p == !IsRationalTV) {
			    /* Create the resulting surface. */
			    if (CAGD_IS_BEZIER_CRV(TCrv2))
			        CmpsCrv = BzrCrvNew(TCrv2 -> Length,
						    IsRational ?
						        CAGD_PT_P3_TYPE :
						        CAGD_PT_E3_TYPE);
			    else {
			        CmpsCrv = BspCrvNew(TCrv2 -> Length,
						    TCrv2 -> Order,
						    IsRational ?
						        CAGD_PT_P3_TYPE :
						        CAGD_PT_E3_TYPE);
				CAGD_GEN_COPY(CmpsCrv -> KnotVector,
					      TCrv2 -> KnotVector,
					      (TCrv2 -> Length +
					          TCrv2 -> Order) *
					          sizeof(CagdRType));
			    }

			    CmpsLength = CmpsCrv -> Length;

			    if (IsRational) {
				/* Reset the weights to 1.0. */
				CmpsPoints = CmpsCrv -> Points[0];
				for (l = 0; l < CmpsLength; l++)
				    CmpsPoints[l] = 1.0;
			    }
			}

			CmpsPoints = CmpsCrv -> Points[p];
			for (l = 0; l < CmpsLength; l++)
			    CmpsPoints[l] = TCrvPoints[l] * SPt;
		    }
		    else {
		        CmpsPoints = CmpsCrv -> Points[p];
			for (l = 0; l < CmpsLength; l++)
			    CmpsPoints[l] += TCrvPoints[l] * SPt;
		    }
		}

		CagdCrvFree(TCrv2);
	    }
	}
    }

    if (IsRationalCrv && !IsRationalTV) {
        CagdCrvStruct *CrvW, *CrvX, *CrvY, *CrvZ, *TCrv3,
	    *TCrv1 = SymbCrvMult(CrvUFactors[UOrder], CrvVFactors[VOrder]),
	    *TCrv2 = SymbCrvMult(TCrv1, CrvWFactors[WOrder]);

	CagdCrvFree(TCrv1);

	SymbCrvSplitScalar(CmpsCrv, &CrvW, &CrvX, &CrvY, &CrvZ);
	CagdCrvFree(CmpsCrv);

	if (CrvW != NULL) {
	    /* See if CrvW is identically one (in which case we ignore it). */
	    for (l = 0; l < CrvW -> Length; l++) {
		if (CrvW -> Points[1][l] != 1.0)
		    break;
	    }
	    if (l < CrvW -> Length) {
	        TCrv1 = SymbCrvMult(TCrv2, CrvW);
		CagdCrvFree(TCrv2);
		TCrv2 = TCrv1;
	    }

	    CagdCrvFree(CrvW);
	}

	TCrv3 = SymbCrvMergeScalar(TCrv2, CrvX, CrvY, CrvZ);
	CagdCrvFree(TCrv2);
	CagdCrvFree(CrvX);
	CagdCrvFree(CrvY);
	CagdCrvFree(CrvZ);
	CmpsCrv = TCrv3;
    }

    for (i = 0; CrvUFactors[i] != NULL; i++)
        CagdCrvFree(CrvUFactors[i]);
    for (i = 0; CrvVFactors[i] != NULL; i++)
        CagdCrvFree(CrvVFactors[i]);
    for (i = 0; CrvWFactors[i] != NULL; i++)
        CagdCrvFree(CrvWFactors[i]);

    IritFree(CrvUFactors);
    CrvUFactors = NULL;
    IritFree(CrvVFactors);
    CrvVFactors = NULL;
    IritFree(CrvWFactors);
    CrvWFactors = NULL;

    return CmpsCrv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given surface Srf and trivar TV, computes the composition TV(Srf).         M
*   Srf must be a three dimensional surface completely contained in the      M
* parametric domain of TV.                                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   TV, Srf:   The trivar and surface to compose. TV must be a Bezier.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:    The resulting composition.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeCrvCrv, SymbComposeSrfCrv, SymbComposeSrfSrf		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivComposeTVSrf, composition                                            M
*****************************************************************************/
CagdSrfStruct *TrivComposeTVSrf(const TrivTVStruct *TV,
				const CagdSrfStruct *Srf)
{
    switch (TV -> GType) {
	case TRIV_TVBEZIER_TYPE:
	    break;
	case TRIV_TVBSPLINE_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_BZR_TV_EXPECT);
	    break;
	case TRIV_TVPOWER_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_POWER_NO_SUPPORT);
	    break;
	default:
	    TRIV_FATAL_ERROR(TRIV_ERR_UNDEF_SRF);
	    break;
    }

    switch (Srf -> GType) {
	case CAGD_SBEZIER_TYPE:
	    break;
	case CAGD_SBSPLINE_TYPE:
	    break;
	case CAGD_SPOWER_TYPE:
	    TRIV_FATAL_ERROR(TRIV_ERR_POWER_NO_SUPPORT);
	    break;
	default:
	    TRIV_FATAL_ERROR(TRIV_ERR_UNDEF_SRF);
	    break;
    }

    return TrivComposeTVSrfAux(TV, Srf);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Aux. function. Subdivides Srf until it is a Bezier surface, invokes the    *
* Bezier composition code on each, and merges them back to complete surface. *
*   At this point, the surface can be either Bezier or B-spline only.        *
*   Srf is assumed to have open end condition.	   		             *
*                                                                            *
* PARAMETERS:                                                                *
*   TV, Srf:   The trivar and surface to compose. TV must be a Bezier.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdSrfStruct *:  The composed surface.                                  *
*****************************************************************************/
static CagdSrfStruct *TrivComposeTVSrfAux(const TrivTVStruct *TV,
					  const CagdSrfStruct *Srf)
{
    CagdRType t;
    CagdSrfStruct *CmpsSrf;

    if (Srf -> ULength > Srf -> UOrder) {
	/* Srf is not a Bezier segment in U. Subdivide, compute for each  */
	/* segment an`g41d merge back.				           */
	CagdSrfStruct *SrfA, *SrfB, *CmpsSrfA, *CmpsSrfB;

	t = Srf -> UKnotVector[(Srf -> UOrder + Srf -> ULength) >> 1];

	SrfA = CagdSrfSubdivAtParam(Srf, t, CAGD_CONST_U_DIR);
	SrfB = SrfA -> Pnext;
	SrfA -> Pnext = NULL;

	CmpsSrfA = TrivComposeTVSrfAux(TV, SrfA);
	CmpsSrfB = TrivComposeTVSrfAux(TV, SrfB);
	CagdSrfFree(SrfA);
	CagdSrfFree(SrfB);

	CmpsSrf = CagdMergeSrfSrf(CmpsSrfA, CmpsSrfB,
				  CAGD_CONST_U_DIR, TRUE, FALSE);
	CagdSrfFree(CmpsSrfA);
	CagdSrfFree(CmpsSrfB);
    }
    else if (Srf -> VLength > Srf -> VOrder) {
	/* Srf is not a Bezier segment in U. Subdivide, compute for each  */
	/* segment and merge back.				           */
	CagdSrfStruct *SrfA, *SrfB, *CmpsSrfA, *CmpsSrfB;

	t = Srf -> VKnotVector[(Srf -> VOrder + Srf -> VLength) >> 1];

	SrfA = CagdSrfSubdivAtParam(Srf, t, CAGD_CONST_V_DIR);
	SrfB = SrfA -> Pnext;
	SrfA -> Pnext = NULL;

	CmpsSrfA = TrivComposeTVSrfAux(TV, SrfA);
	CmpsSrfB = TrivComposeTVSrfAux(TV, SrfB);
	CagdSrfFree(SrfA);
	CagdSrfFree(SrfB);

	CmpsSrf = CagdMergeSrfSrf(CmpsSrfA, CmpsSrfB,
				  CAGD_CONST_V_DIR, TRUE, FALSE);
	CagdSrfFree(CmpsSrfA);
	CagdSrfFree(CmpsSrfB);
    }
    else {
    	/* Srf is a Bezier surface segment - compute its composition. */
	if (!CAGD_IS_BEZIER_SRF(Srf)) {
	    CagdSrfStruct 
		*TSrf = CagdCnvrtBsp2BzrSrf(Srf);

	    CmpsSrf = TrivBzrComposeTVSrf(TV, TSrf);
	    CagdSrfFree(TSrf);
	}
	else
    	    CmpsSrf = TrivBzrComposeTVSrf(TV, Srf);
    }

    return CmpsSrf;    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given surface Srf and Bezier trivariate TV, computes their composition   M
* TV(Srf(a, b)),  Srf(a, b) = (u(a, b), v(a, b), w(a, b)).		     M
*   Srf must be a three dimensional surface completely contained in the      M
* parametric domain of TV.                                                   M
*   Compute the compositions by the products of:			     M
*									     M
* TV(u, v, w) = TV(u(a, b), v(a, b), w(a, b))				     V
*									     V
*            n   m   l							     V
*         = sum sum sum Pijk Bi(u(a, b)) Bj(v(a, b)) Bk(w(a, b))	     V
*           i=0 j=0 k=0							     V
*				   					     V
*            n   m   l	      n	          i	         n-i		     V
*         = sum sum sum Pijk ( ) (u(a, b))  (1 - u(a, b))                    V
*           i=0 j=0 k=0	      i						     V
*                   	         m	     j	            m-j		     V
*                               ( ) (v(a, b))  (1 - v(a, b))                 V
*                  	         j					     V
*                   	            l	        k	       l-k	     V
*                                  ( ) (w(a, b))  (1 - w(a, b))              V
*                  	            k					     V
*                                                                            *
* PARAMETERS:                                                                M
*   TV, Srf:   The trivar and surface to compose. TV must be Bezier.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:    The resulting composition.                           M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbComposeTVSrf, BzrComposeSrfSrf, TrivBzrComposeTVCrv		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivBzrComposeTVSrf, composition                                         M
*****************************************************************************/
CagdSrfStruct *TrivBzrComposeTVSrf(const TrivTVStruct *TV,
				   const CagdSrfStruct *Srf)
{
    IRIT_STATIC_DATA CagdSrfStruct
	**SrfUFactors = NULL,
        **SrfVFactors = NULL,
	**SrfWFactors = NULL;
    CagdBType
        IsRationalSrf = CAGD_IS_RATIONAL_SRF(Srf),
	IsRationalTV = TRIV_IS_RATIONAL_TV(TV),
        IsRational = IsRationalSrf | IsRationalTV;
    int i, j, k, l,
        CmpsLength = 0,
	UOrder = TV -> UOrder,
	VOrder = TV -> VOrder,
	WOrder = TV -> WOrder,
	MaxCoord = CAGD_NUM_OF_PT_COORD(TV -> PType);
    CagdRType * const
        *TVPoints = TV -> Points;
    CagdSrfStruct *SrfUVW,
    	*CmpsSrf = NULL;

    if (CAGD_NUM_OF_PT_COORD(Srf -> PType) < 3 || MaxCoord != 3) {
	TRIV_FATAL_ERROR(TRIV_ERR_UNSUPPORT_PT);
	return NULL;
    }

    SrfUVW = CagdCoerceSrfTo(Srf,
			     IsRationalSrf ? CAGD_PT_P1_TYPE : CAGD_PT_E1_TYPE,
			     FALSE);

    if (SrfUFactors != NULL) {
        for (i = 0; SrfUFactors[i] != NULL; i++)
	    CagdSrfFree(SrfUFactors[i]);
	IritFree(SrfUFactors);
    }
    SrfUFactors = SymbComputeSurfacePowers(SrfUVW, UOrder);

    CAGD_GEN_COPY(SrfUVW -> Points[1], Srf -> Points[2],
		  sizeof(CagdRType) * Srf -> ULength * Srf -> VLength);

    if (SrfVFactors != NULL) {
        for (i = 0; SrfVFactors[i] != NULL; i++)
	    CagdSrfFree(SrfVFactors[i]);
	IritFree(SrfVFactors);
    }
    SrfVFactors = SymbComputeSurfacePowers(SrfUVW, VOrder);

    CAGD_GEN_COPY(SrfUVW -> Points[1], Srf -> Points[3],
		  sizeof(CagdRType) * Srf -> ULength * Srf -> VLength);

    if (SrfWFactors != NULL) {
        for (i = 0; SrfWFactors[i] != NULL; i++)
	    CagdSrfFree(SrfWFactors[i]);
	IritFree(SrfWFactors);
    }
    SrfWFactors = SymbComputeSurfacePowers(SrfUVW, WOrder);

    CagdSrfFree(SrfUVW);

    /* The main (triple) loop (Compositions computation): */
    for (k = 0; k < WOrder; k++) {
        for (j = 0; j < VOrder; j++) {
	    for (i = 0; i < UOrder; i++) {
	        int p,
		    Idx = TRIV_MESH_UVW(TV, i, j, k);
	        CagdSrfStruct
		    *TSrf1 = SymbSrfMult(SrfUFactors[i], SrfVFactors[j]),
		    *TSrf2 = SymbSrfMult(SrfWFactors[k], TSrf1);
		CagdRType
		    *TSrfPoints = TSrf2 -> Points[1];

		CagdSrfFree(TSrf1);

		for (p = !IsRationalTV; p <= MaxCoord; p++) {
		    CagdRType *CmpsPoints,
		        SPt = TVPoints[p][Idx];

		    if (i == 0 && j == 0 && k == 0) {
		        if (p == !IsRationalTV) {
			    /* Create the resulting surface. */
			    if (CAGD_IS_BEZIER_SRF(TSrf2))
			        CmpsSrf = BzrSrfNew(TSrf2 -> ULength,
						    TSrf2 -> VLength,
						    IsRational ?
						        CAGD_PT_P3_TYPE :
						        CAGD_PT_E3_TYPE);
			    else {
			        CmpsSrf = BspSrfNew(TSrf2 -> ULength,
						    TSrf2 -> VLength,
						    TSrf2 -> UOrder,
						    TSrf2 -> VOrder,
						    IsRational ?
						        CAGD_PT_P3_TYPE :
						        CAGD_PT_E3_TYPE);
				CAGD_GEN_COPY(CmpsSrf -> UKnotVector,
					      TSrf2 -> UKnotVector,
					      (TSrf2 -> ULength +
					          TSrf2 -> UOrder) *
					          sizeof(CagdRType));
				CAGD_GEN_COPY(CmpsSrf -> VKnotVector,
					      TSrf2 -> VKnotVector,
					      (TSrf2 -> VLength +
					          TSrf2 -> VOrder) *
					          sizeof(CagdRType));
			    }

			    CmpsLength = CmpsSrf -> ULength *
			                 CmpsSrf -> VLength;

			    if (IsRational) {
			        /* Reset the weights to 1.0. */
			        CmpsPoints = CmpsSrf -> Points[0];
				for (l = 0; l < CmpsLength; l++)
				    CmpsPoints[l] = 1.0;
			    }
			}

			CmpsPoints = CmpsSrf -> Points[p];
			for (l = 0; l < CmpsLength; l++)
			    CmpsPoints[l] = TSrfPoints[l] * SPt;
		    }
		    else {
		        CmpsPoints = CmpsSrf -> Points[p];
			for (l = 0; l < CmpsLength; l++)
			    CmpsPoints[l] += TSrfPoints[l] * SPt;
		    }
		}

		CagdSrfFree(TSrf2);
	    }
	}
    }

    if (IsRationalSrf && !IsRationalTV) {
        CagdSrfStruct *SrfW, *SrfX, *SrfY, *SrfZ, *TSrf3,
	    *TSrf1 = SymbSrfMult(SrfUFactors[UOrder], SrfVFactors[VOrder]),
	    *TSrf2 = SymbSrfMult(TSrf1, SrfWFactors[WOrder]);

	CagdSrfFree(TSrf1);

	SymbSrfSplitScalar(CmpsSrf, &SrfW, &SrfX, &SrfY, &SrfZ);
	CagdSrfFree(CmpsSrf);

	if (SrfW != NULL) {
	    /* See if SrfW is identically one (in which case we ignore it). */
	    for (l = 0; l < SrfW -> ULength * SrfW -> VLength; l++) {
		if (SrfW -> Points[1][l] != 1.0)
		    break;
	    }
	    if (l < SrfW -> ULength * SrfW -> VLength) {
	        TSrf1 = SymbSrfMult(TSrf2, SrfW);
		CagdSrfFree(TSrf2);
		TSrf2 = TSrf1;
	    }

	    CagdSrfFree(SrfW);
	}

	TSrf3 = SymbSrfMergeScalar(TSrf2, SrfX, SrfY, SrfZ);

	CagdSrfFree(TSrf2);
	CagdSrfFree(SrfX);
	CagdSrfFree(SrfY);
	CagdSrfFree(SrfZ);
	CmpsSrf = TSrf3;        
    }

    for (i = 0; SrfUFactors[i] != NULL; i++)
        CagdSrfFree(SrfUFactors[i]);
    for (i = 0; SrfVFactors[i] != NULL; i++)
        CagdSrfFree(SrfVFactors[i]);
    for (i = 0; SrfWFactors[i] != NULL; i++)
        CagdSrfFree(SrfWFactors[i]);

    IritFree(SrfUFactors);
    SrfUFactors = NULL;
    IritFree(SrfVFactors);
    SrfVFactors = NULL;
    IritFree(SrfWFactors);
    SrfWFactors = NULL;

    return CmpsSrf;
}
