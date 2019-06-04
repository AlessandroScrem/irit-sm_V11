/******************************************************************************
* triv_ffd.c - perform ffd (freeform deformations) using trivariates.	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber,		 Sep 2012.			      *
******************************************************************************/

#include <stdio.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "triv_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deform the given mesh in place, using the mapping that is defined by     M
* trivariate DeformTV.  Input points that are outside the domain of DeformTV M
* are coerced to the closest boundary.                                       M
*   Computation is approximated by mapping (control) points only.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Points:       The control mesh.                                          M
*   Length:       The length of the vector of Points.                        M
*   PType:        The point type of Points.                                  M
*   DeformTV:     The deformation mapping.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivFFDObjectTV                                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivFFDCtlMeshUsingTV                                                    M
*****************************************************************************/
void TrivFFDCtlMeshUsingTV(CagdRType **Points,
			   int Length,
			   CagdPointType PType,
			   const TrivTVStruct *DeformTV)
{
    int i;
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax;

    TrivTVDomain(DeformTV, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);

    for (i = 0; i < Length; i++) {
        CagdRType *R;
	CagdPType PtE3;

 	CagdCoerceToE3(PtE3, Points, i, PType);

	TRIV_COERCe_TO_DOMAIN(PtE3, UMin, UMax, VMin, VMax, WMin, WMax);
	R = TrivTVEval(DeformTV, PtE3[0], PtE3[1], PtE3[2]);
	CagdCoerceToE3(PtE3, &R, -1, DeformTV -> PType);

	switch (PType) {
	    case CAGD_PT_E9_TYPE:
	    case CAGD_PT_E8_TYPE:
	    case CAGD_PT_E7_TYPE:
	    case CAGD_PT_E6_TYPE:
	    case CAGD_PT_E5_TYPE:
	    case CAGD_PT_E4_TYPE:
	    case CAGD_PT_E3_TYPE:
	        Points[3][i] = PtE3[2];
	    case CAGD_PT_E2_TYPE:
	        Points[2][i] = PtE3[1];
	    case CAGD_PT_E1_TYPE:
	        Points[1][i] = PtE3[0];
		break;
	    case CAGD_PT_P9_TYPE:
	    case CAGD_PT_P8_TYPE:
	    case CAGD_PT_P7_TYPE:
	    case CAGD_PT_P6_TYPE:
	    case CAGD_PT_P5_TYPE:
	    case CAGD_PT_P4_TYPE:
	    case CAGD_PT_P3_TYPE:
	        Points[3][i] = PtE3[2] * Points[0][i];
	    case CAGD_PT_P2_TYPE:
	        Points[2][i] = PtE3[1] * Points[0][i];
	    case CAGD_PT_P1_TYPE:
	        Points[1][i] = PtE3[0] * Points[0][i];
		break;
	    default:
	        break;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deform an input object, in place, through the given trivariate.          M
*   Computation is approximated by mapping (control) points only.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:       The object to map through the trivariate, in place.          M
*   DeformTV:   The mapping/deformation function from R3 to R3.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   PObj, mapped/deformed object, in place.              M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivFFDCtlMeshUsingTV, TrivFFDTileObjectInTV                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivFFDObjectTV                                                          M
*****************************************************************************/
IPObjectStruct *TrivFFDObjectTV(IPObjectStruct *PObj,
				const TrivTVStruct *DeformTV)
{
    int i;
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax, *R, Pt[4];
    IPPolygonStruct *Pl;
    CagdCrvStruct *Crv;
    CagdSrfStruct *Srf;
    TrimSrfStruct *TSrf;
    TrivTVStruct *TV;
    TrngTriangSrfStruct *TriSrf;
    MdlModelStruct *Mdl;
    MvarMVStruct *MV;
    IPObjectStruct *PTmp;

    TrivTVDomain(DeformTV, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);

    switch (PObj -> ObjType) {
	case IP_OBJ_UNDEF:
	case IP_OBJ_MATRIX:
	case IP_OBJ_NUMERIC:
	case IP_OBJ_INSTANCE:
	case IP_OBJ_STRING:
	case IP_OBJ_PLANE:
	    break;
	case IP_OBJ_POLY:
	    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	        IPVertexStruct
		    *V = Pl -> PVertex;

		do {
		    TRIV_COERCe_TO_DOMAIN(V -> Coord,
					  UMin, UMax, VMin, VMax, WMin, WMax);
		    R = TrivTVEval(DeformTV,
				   V -> Coord[0],
				   V -> Coord[1],
				   V -> Coord[2]);
		    CagdCoerceToE3(V -> Coord, &R, -1, DeformTV -> PType);

		    V = V -> Pnext;
		}
		while (V != NULL && V != Pl -> PVertex);
	    }
	    break;
	case IP_OBJ_POINT:
	    TRIV_COERCe_TO_DOMAIN(PObj -> U.Pt, 
				  UMin, UMax, VMin, VMax, WMin, WMax);
	    R = TrivTVEval(DeformTV,
			   PObj -> U.Pt[0],
			   PObj -> U.Pt[1],
			   PObj -> U.Pt[2]);
	    CagdCoerceToE3(PObj -> U.Pt, &R, -1, DeformTV -> PType);
	    break;
	case IP_OBJ_VECTOR:
	    TRIV_COERCe_TO_DOMAIN(PObj -> U.Pt,
				  UMin, UMax, VMin, VMax, WMin, WMax);
	    R = TrivTVEval(DeformTV,
			   PObj -> U.Vec[0],
			   PObj -> U.Vec[1],
			   PObj -> U.Vec[2]);
	    CagdCoerceToE3(PObj -> U.Vec, &R, -1, DeformTV -> PType);
	    break;
	case IP_OBJ_CTLPT:
	    R = PObj -> U.CtlPt.Coords;
	    CagdCoercePointTo(&Pt[1], CAGD_PT_E3_TYPE,
			      (IrtRType * const *) &R,
			      -1, PObj -> U.CtlPt.PtType);
	    TRIV_COERCe_TO_DOMAIN(PObj -> U.Pt,
				  UMin, UMax, VMin, VMax, WMin, WMax);
	    R = TrivTVEval(DeformTV, Pt[0], Pt[1], Pt[2]);
	    CagdCoercePointTo(PObj -> U.CtlPt.Coords,
			      PObj -> U.CtlPt.PtType,
			      (IrtRType * const *) &R,
			      -1, CAGD_PT_E3_TYPE);
	    break;
	case IP_OBJ_LIST_OBJ:
	    for (i = 0; (PTmp = IPListObjectGet(PObj, i)) != NULL; i++)
	        TrivFFDObjectTV(PTmp, DeformTV);
	    break;
	case IP_OBJ_CURVE:
	    for (Crv = PObj -> U.Crvs; Crv != NULL; Crv = Crv -> Pnext)
	        TrivFFDCtlMeshUsingTV(Crv -> Points, Crv -> Length,
				      Crv -> PType, DeformTV);
	    break;
	case IP_OBJ_SURFACE:
	    for (Srf = PObj -> U.Srfs; Srf != NULL; Srf = Srf -> Pnext)
	        TrivFFDCtlMeshUsingTV(Srf -> Points,
				      Srf -> ULength * Srf -> VLength,
				      Srf -> PType, DeformTV);
	    break;
	case IP_OBJ_TRIMSRF:
	    for (TSrf = PObj -> U.TrimSrfs; TSrf != NULL; TSrf = TSrf -> Pnext) {
	        /* Note we ignore E3 trimmings curves - might be required. */
	        TrivFFDCtlMeshUsingTV(TSrf -> Srf -> Points,
				      TSrf -> Srf -> ULength *
				          TSrf -> Srf -> VLength,
				      TSrf -> Srf -> PType, DeformTV);
	    }
	    break;
	case IP_OBJ_TRIVAR:
	    for (TV = PObj -> U.Trivars; TV != NULL; TV = TV -> Pnext)
	        TrivFFDCtlMeshUsingTV(TV -> Points,
				      TV -> ULength *
				          TV -> VLength *
				              TV -> WLength,
				      TV -> PType, DeformTV);
	    break;
	case IP_OBJ_TRISRF:
	    for (TriSrf = PObj -> U.TriSrfs;
		 TriSrf != NULL;
		 TriSrf = TriSrf -> Pnext)
	        TrivFFDCtlMeshUsingTV(TriSrf -> Points,
				      TRNG_TRISRF_MESH_SIZE(TriSrf),
				      TriSrf -> PType, DeformTV);
	    break;
	case IP_OBJ_MODEL:
	    for (Mdl = PObj -> U.Mdls; Mdl != NULL; Mdl = Mdl -> Pnext) {
	        MdlTrimSrfStruct *TS;

		for (TS = Mdl -> TrimSrfList; TS != NULL; TS = TS -> Pnext) {
		    TrivFFDCtlMeshUsingTV(TS -> Srf -> Points,
					  TS -> Srf -> ULength *
					  TS -> Srf -> VLength,
					  TS -> Srf -> PType, DeformTV);
		}
	    }
	    break;
	    break;
	case IP_OBJ_MULTIVAR:
	    for (MV = PObj -> U.MultiVars; MV != NULL; MV = MV -> Pnext)
	        TrivFFDCtlMeshUsingTV(MV -> Points, MVAR_CTL_MESH_LENGTH(MV),
				      (CagdPointType) MV -> PType, DeformTV);
	    break;
	default:
	    IP_FATAL_ERROR(IP_ERR_UNDEF_OBJECT_FOUND);
    }

    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tile an input object, in place, (UTimes x VTimes x WTimes) in the given  M
* trivariate.  Computation is approximated using (control) points mapping.   M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:       The object to map through the trivariate, in place.          M
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
*   TrivFFDCtlMeshUsingTV, TrivFFDObjectTV, TrivComposeTileObjectInTV        M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivFFDTileObjectInTV                                                    M
*****************************************************************************/
IPObjectStruct *TrivFFDTileObjectInTV(const IPObjectStruct *PObj,
				      const TrivTVStruct *DeformTV,
				      IrtRType UTimes,
				      IrtRType VTimes,
				      IrtRType WTimes,
				      IrtBType FitObj)
{
    int u, v, w;
    CagdRType UMin, UMax, VMin, VMax, WMin, WMax;
    IPObjectStruct *PNewObj, *PTile,
        *PRetVal = IPGenLISTObject(NULL);
    GMBBBboxStruct
        *BBox = GMBBComputeBboxObject(PObj);
    IrtHmgnMatType Mat1, Mat2;

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
		MatGenMatTrans(UMin + u * (UMax - UMin) / UTimes,
			       VMin + v * (VMax - VMin) / VTimes,
			       WMin + w * (WMax - WMin) / WTimes,
			       Mat1);
		PNewObj = GMTransformObject(PTile, Mat1);
		TrivFFDObjectTV(PNewObj, DeformTV);
		IPListObjectAppend(PRetVal, PNewObj);
	    }
	}
    }

    IPFreeObject(PTile);

    return PRetVal;
}
