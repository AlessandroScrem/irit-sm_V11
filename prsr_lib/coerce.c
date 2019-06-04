/******************************************************************************
* coerce.c - coerce an object into a different type.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, January 1993.				      *
******************************************************************************/

#include <stdio.h>
#include "inc_irit/irit_sm.h"
#include "prsr_loc.h"
#include "inc_irit/bool_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/trim_lib.h"
#include "inc_irit/mdl_lib.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"

#define CAGD_PT_E_TYPES  case CAGD_PT_E1_TYPE: \
			 case CAGD_PT_E2_TYPE: \
			 case CAGD_PT_E3_TYPE: \
			 case CAGD_PT_E4_TYPE: \
			 case CAGD_PT_E5_TYPE: \
			 case CAGD_PT_E6_TYPE: \
			 case CAGD_PT_E7_TYPE: \
			 case CAGD_PT_E8_TYPE: \
			 case CAGD_PT_E9_TYPE: \
			 case CAGD_PT_E10_TYPE: \
			 case CAGD_PT_E11_TYPE: \
			 case CAGD_PT_E12_TYPE: \
			 case CAGD_PT_E13_TYPE: \
			 case CAGD_PT_E14_TYPE: \
			 case CAGD_PT_E15_TYPE: \
			 case CAGD_PT_E16_TYPE: \
			 case CAGD_PT_E17_TYPE: \
			 case CAGD_PT_E18_TYPE:

#define CAGD_PT_P_TYPES  case CAGD_PT_P1_TYPE: \
			 case CAGD_PT_P2_TYPE: \
			 case CAGD_PT_P3_TYPE: \
			 case CAGD_PT_P4_TYPE: \
			 case CAGD_PT_P5_TYPE: \
			 case CAGD_PT_P6_TYPE: \
			 case CAGD_PT_P7_TYPE: \
			 case CAGD_PT_P8_TYPE: \
			 case CAGD_PT_P9_TYPE: \
			 case CAGD_PT_P10_TYPE: \
			 case CAGD_PT_P11_TYPE: \
			 case CAGD_PT_P12_TYPE: \
			 case CAGD_PT_P13_TYPE: \
			 case CAGD_PT_P14_TYPE: \
			 case CAGD_PT_P15_TYPE: \
			 case CAGD_PT_P16_TYPE: \
			 case CAGD_PT_P17_TYPE: \
			 case CAGD_PT_P18_TYPE:

#define CAGD_PT_TYPES    \
			 CAGD_PT_E_TYPES \
			 CAGD_PT_P_TYPES

static IPObjectStruct *CurveReverse(const IPObjectStruct *CrvObj);
static IPObjectStruct *SurfaceReverse(const IPObjectStruct *SrfObj);
static IPObjectStruct *ModelReverse(const IPObjectStruct *SrfObj);

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Given a set of points, returns the list's common denominator that spans   M
* the space of all the points, taking into account type Type.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   PtObjList:  List of points.                                              M
*   Type:       Point type that we must span its space as well.              M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPointType:   Point type that spans the space of point type Type as   M
*                    well as all points in PtObjList.                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceCommonSpace, coercion                                            M
*****************************************************************************/
CagdPointType IPCoerceCommonSpace(IPObjectStruct *PtObjList,
				  CagdPointType Type)
{
    int i,
	Dim = CAGD_NUM_OF_PT_COORD(Type),
	Proj = CAGD_IS_RATIONAL_PT(Type);
    IPObjectStruct *PtObj;

    if (!IP_IS_OLST_OBJ(PtObjList)) {
	IP_FATAL_ERROR(IP_ERR_LIST_OBJ_EXPECTED);
	return CAGD_PT_NONE;
    }

    for (i = 0; (PtObj = IPListObjectGet(PtObjList, i++)) != NULL; ) {
	if (IP_IS_CTLPT_OBJ(PtObj)) {
	    Dim = IRIT_MAX(Dim, CAGD_NUM_OF_PT_COORD(PtObj -> U.CtlPt.PtType));
	    Proj |= CAGD_IS_RATIONAL_PT(PtObj -> U.CtlPt.PtType);
	}
	else if (IP_IS_POINT_OBJ(PtObj) || IP_IS_VEC_OBJ(PtObj)) {
	    Dim = IRIT_MAX(Dim, 3);
	}
	else {
	    IP_FATAL_ERROR(IP_ERR_PT_OBJ_EXPECTED);
	    return CAGD_PT_NONE;
	}
    }

    return CAGD_MAKE_PT_TYPE(Proj, Dim);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Converts a Gregory freeform into a Bezier freeform.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:      Gregory geometry to convert to Bezier geometry.               M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   Same geometry as PObj but in Gregory basis.          M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceGregoryToBezier                                                  M
*****************************************************************************/
IPObjectStruct *IPCoerceGregoryToBezier(const IPObjectStruct *PObj)
{
    if (PObj -> ObjType == IP_OBJ_TRISRF) {
	TrngTriangSrfStruct
	    *TriSrf = TrngCnvrtGregory2BzrTriSrf(PObj -> U.TriSrfs);

	return IPGenTRISRFObject(TriSrf);
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_TRI_SRF);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Converts a Bezier freeform into a power freeform.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:      Bezier geometry to convert to power geometry.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   Same geometry as PObj but in power basis.            M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceBezierToPower                                                    M
*****************************************************************************/
IPObjectStruct *IPCoerceBezierToPower(const IPObjectStruct *PObj)
{
    if (PObj -> ObjType == IP_OBJ_SURFACE) {
	CagdSrfStruct
	    *Srf = CagdCnvrtBzr2PwrSrf(PObj -> U.Srfs);

	return IPGenSRFObject(Srf);
    }
    else if (PObj -> ObjType == IP_OBJ_CURVE) {
	CagdCrvStruct
	    *Crv = CagdCnvrtBzr2PwrCrv(PObj -> U.Crvs);

	return IPGenCRVObject(Crv);
    }
    else if (PObj -> ObjType == IP_OBJ_MULTIVAR) {
	MvarMVStruct
	    *MV = MvarCnvrtBzr2PwrMV(PObj -> U.MultiVars);

	return IPGenMULTIVARObject(MV);
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_CRV_SRF_MV);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Converts a power freeform into a Bezier freeform.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:      Power geometry to convert to Bezier geometry.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   Same geometry as PObj but in Bezier form.            M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoercePowerToBezier                                                    M
*****************************************************************************/
IPObjectStruct *IPCoercePowerToBezier(const IPObjectStruct *PObj)
{
    if (PObj -> ObjType == IP_OBJ_SURFACE) {
	CagdSrfStruct
	    *Srf = CagdCnvrtPwr2BzrSrf(PObj -> U.Srfs);

	return IPGenSRFObject(Srf);
    }
    else if (PObj -> ObjType == IP_OBJ_CURVE) {
	CagdCrvStruct
	    *Crv = CagdCnvrtPwr2BzrCrv(PObj -> U.Crvs);

	return IPGenCRVObject(Crv);
    }
    else if (PObj -> ObjType == IP_OBJ_MULTIVAR) {
	MvarMVStruct
	    *MV = MvarCnvrtPwr2BzrMV(PObj -> U.MultiVars);

	return IPGenMULTIVARObject(MV);
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_TRI_SRF);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Converts a Bezier freeform into a Bspline freeform.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:      Bezier geometry to convert to Bspline geometry.               M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   Same geometry as PObj but as Bspline.                M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceBezierToBspline                                                  M
*****************************************************************************/
IPObjectStruct *IPCoerceBezierToBspline(const IPObjectStruct *PObj)
{
    if (PObj -> ObjType == IP_OBJ_SURFACE) {
	CagdSrfStruct
	    *Srf = CagdCnvrtBzr2BspSrf(PObj -> U.Srfs);

	return IPGenSRFObject(Srf);
    }
    else if (PObj -> ObjType == IP_OBJ_CURVE) {
	CagdCrvStruct
	    *Crv = CagdCnvrtBzr2BspCrv(PObj -> U.Crvs);

	return IPGenCRVObject(Crv);
    }
    else if (PObj -> ObjType == IP_OBJ_MULTIVAR) {
	MvarMVStruct
	    *MV = MvarCnvrtBzr2BspMV(PObj -> U.MultiVars);

	return IPGenMULTIVARObject(MV);
    }
    else if (PObj -> ObjType == IP_OBJ_TRIVAR) {
	TrivTVStruct
	    *TV = TrivCnvrtBzr2BspTV(PObj -> U.Trivars);

	return IPGenTRIVARObject(TV);
    }
    else if (PObj -> ObjType == IP_OBJ_TRISRF) {
	TrngTriangSrfStruct
	    *TriSrf = TrngCnvrtBzr2BspTriSrf(PObj -> U.TriSrfs);

	return IPGenTRISRFObject(TriSrf);
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_FREEFORM);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Convert a Bspline freeform into list of Bezier freeforms.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:     A Bspline geometry to convert to a Bezier geometry.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A Bezier geometry representing same geometry as PObj. M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceBsplineToBezier                                                  M
*****************************************************************************/
IPObjectStruct *IPCoerceBsplineToBezier(const IPObjectStruct *PObj)
{
    IPObjectStruct *PObjList, *PTmp;

    if (PObj -> ObjType == IP_OBJ_SURFACE) {
	CagdSrfStruct *TSrf,
	    *Srf = CagdCnvrtBsp2BzrSrf(PObj -> U.Srfs);

	if (Srf -> Pnext == NULL) {
	    PTmp = IPGenSRFObject(Srf);
	    IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);
	    return PTmp;
	}
	else {
	    int i;

	    PObjList = IPGenLISTObject(NULL);
	    PObjList -> Attr = IP_ATTR_COPY_ATTRS(PObj -> Attr);

	    for (i = 0; Srf != NULL; i++) {
	        IPListObjectInsert(PObjList, i, PTmp = IPGenSRFObject(Srf));
		PTmp -> Attr = IP_ATTR_COPY_ATTRS(Srf -> Attr);

		TSrf = Srf -> Pnext;
		Srf -> Pnext = NULL;
		Srf = TSrf;
	    }
	    IPListObjectInsert(PObjList, i, NULL);

	    return PObjList;
	}				   
    }
    else if (PObj -> ObjType == IP_OBJ_CURVE) {
	CagdCrvStruct *TCrv,
	    *Crv = CagdCnvrtBsp2BzrCrv(PObj -> U.Crvs);

	if (Crv -> Pnext == NULL) {
	    PTmp = IPGenCRVObject(Crv);
	    IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);
	    return PTmp;
	}
	else {
	    int i;

	    PObjList = IPGenLISTObject(NULL);

	    for (i = 0; Crv != NULL; i++) {
	        IPListObjectInsert(PObjList, i, PTmp = IPGenCRVObject(Crv));
		IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);

		TCrv = Crv -> Pnext;
		Crv -> Pnext = NULL;
		Crv = TCrv;
	    }
	    IPListObjectInsert(PObjList, i, NULL);

	    return PObjList;
	}				   
    }
    else if (PObj -> ObjType == IP_OBJ_MULTIVAR) {
	MvarMVStruct *TMV,
	    *MV = MvarCnvrtBsp2BzrMV(PObj -> U.MultiVars);

	if (MV -> Pnext == NULL) {
	    PTmp = IPGenMULTIVARObject(MV);
	    IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);
	    return PTmp;
	}
	else {
	    int i;

	    PObjList = IPGenLISTObject(NULL);

	    for (i = 0; MV != NULL; i++) {
	        IPListObjectInsert(PObjList, i,
				   PTmp = IPGenMULTIVARObject(MV));
		IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);

		TMV = MV -> Pnext;
		MV -> Pnext = NULL;
		MV = TMV;
	    }
	    IPListObjectInsert(PObjList, i, NULL);

	    return PObjList;
	}				   
    }
    else if (PObj -> ObjType == IP_OBJ_TRIVAR) {
	TrivTVStruct *TTV,
	    *TV = TrivCnvrtBsp2BzrTV(PObj -> U.Trivars);

	if (TV -> Pnext == NULL) {
	    PTmp = IPGenTRIVARObject(TV);
	    IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);
	    return PTmp;
	}
	else {
	    int i;

	    PObjList = IPGenLISTObject(NULL);

	    for (i = 0; TV != NULL; i++) {
	        IPListObjectInsert(PObjList, i, PTmp = IPGenTRIVARObject(TV));
		IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);

		TTV = TV -> Pnext;
		TV -> Pnext = NULL;
		TV = TTV;
	    }
	    IPListObjectInsert(PObjList, i, NULL);

	    return PObjList;
	}				   
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_FREEFORM);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Convert a Bspline freeform into list of Bezier freeforms.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:     A Bspline geometry to convert to a Bezier geometry.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A Bezier geometry representing same geometry as PObj. M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceTrimmedSrfToTrimmedBezier                                        M
*****************************************************************************/
IPObjectStruct *IPCoerceTrimmedSrfToTrimmedBezier(const IPObjectStruct *PObj)
{
    IPObjectStruct *PObjList, *PTmp;

    if (PObj -> ObjType == IP_OBJ_TRIMSRF) {
        TrimSrfStruct *TSrf,
	    *TrimBzrSrfs = TrimSrfCnvrt2BzrTrimSrf(PObj -> U.TrimSrfs);

	if (TrimBzrSrfs -> Pnext == NULL) {
	    PTmp = IPGenTRIMSRFObject(TrimBzrSrfs);
	    IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);
	    return PTmp;
	}
	else {
	    int i;

	    PObjList = IPGenLISTObject(NULL);
	    PObjList -> Attr = IP_ATTR_COPY_ATTRS(PObj -> Attr);

	    for (i = 0; TrimBzrSrfs != NULL; i++) {
	        IPListObjectInsert(PObjList, i,
				   PTmp = IPGenTRIMSRFObject(TrimBzrSrfs));
		PTmp -> Attr = IP_ATTR_COPY_ATTRS(TrimBzrSrfs -> Attr);

		TSrf = TrimBzrSrfs -> Pnext;
		TrimBzrSrfs -> Pnext = NULL;
		TrimBzrSrfs = TSrf;
	    }
	    IPListObjectInsert(PObjList, i, NULL);

	    return PObjList;
	}				   
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_TRIM_SRF);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Convert a trimmed B-spline freeform into untrimmed tensor product Bezier M
* freeforms.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:         A B-spline geometry to convert to a Bezier geometry.       M
*   ComposeE3:    TRUE to compose the tiles into TSrf, FALSE to return the   M
*                 surface tiles in the parametric domain of TSrf.            M
                                                                             *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A Bezier geometry representing same geometry as PObj. M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceTrimmedSrfToUnTrimmedBezier                                      M
*****************************************************************************/
IPObjectStruct *IPCoerceTrimmedSrfToUnTrimmedBezier(const IPObjectStruct *PObj,
						    int ComposeE3)
{
    IPObjectStruct *PObjList, *PTmp;
    int MinimizeSize = AttrGetObjectIntAttrib(PObj, "MinSize");

    if (PObj -> ObjType == IP_OBJ_TRIMSRF) {
        CagdSrfStruct *TSrf, *BzrSrfs;

	if (!IP_ATTR_IS_BAD_INT(MinimizeSize) && MinimizeSize)
	    AttrSetIntAttrib(&PObj -> U.TrimSrfs -> Attr, "MinSize", TRUE);

#	ifdef IP_COERCE_UNTRIM_OLD_FUNC
	    BzrSrfs = TrimSrfCnvrt2BzrRglrSrf(PObj -> U.TrimSrfs);
#	else
	    BzrSrfs = TrimSrfCnvrt2BzrRglrSrf2(PObj -> U.TrimSrfs, ComposeE3,
					       FALSE, TRUE, IRIT_EPS);
#	endif /* IP_COERCE_UNTRIM_OLD_FUNC */

	if (BzrSrfs == NULL)
	    return NULL;
	else if (BzrSrfs -> Pnext == NULL) {
	    PTmp = IPGenSRFObject(BzrSrfs);
	    IP_ATTR_SAFECOPY_ATTRS(PTmp -> Attr, PObj -> Attr);
	    return PTmp;
	}
	else {
	    int i;

	    PObjList = IPGenLISTObject(NULL);
	    PObjList -> Attr = IP_ATTR_COPY_ATTRS(PObj -> Attr);

	    for (i = 0; BzrSrfs != NULL; i++) {
	        IPListObjectInsert(PObjList, i,
				   PTmp = IPGenSRFObject(BzrSrfs));
		PTmp -> Attr = IP_ATTR_COPY_ATTRS(BzrSrfs -> Attr);

		TSrf = BzrSrfs -> Pnext;
		BzrSrfs -> Pnext = NULL;
		BzrSrfs = TSrf;
	    }
	    IPListObjectInsert(PObjList, i, NULL);

	    return PObjList;
	}				   
    }
    else {
        IP_FATAL_ERROR(IP_ERR_ONLY_TRIM_SRF);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Coerces a list of objects to Type.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   PtObjList:   Coerce points/vectors/control points in this list to Type.  M
*   Type:        A minimum space type to coerce to in PtObjList.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPointType:  The coercion type actually took place with in PtObjList. M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoercePtsListTo, coercion                                              M
*****************************************************************************/
CagdPointType IPCoercePtsListTo(IPObjectStruct *PtObjList, CagdPointType Type)
{
    int i;
    IPObjectStruct *TmpObj, *PtObj;
    CagdPointType
	PtType = IPCoerceCommonSpace(PtObjList, Type);

    if (PtType != CAGD_PT_NONE) {
	for (i = 0; (PtObj = IPListObjectGet(PtObjList, i++)) != NULL; ) {
	    if (IP_IS_CTLPT_OBJ(PtObj)) {
		TmpObj = IPCoerceObjectPtTypeTo(PtObj, PtType);
		PtObj -> U.CtlPt = TmpObj -> U.CtlPt;
		IPFreeObject(TmpObj);
	    }
	    else if (IP_IS_POINT_OBJ(PtObj) || IP_IS_VEC_OBJ(PtObj)) {
		TmpObj = IPCoerceObjectPtTypeTo(PtObj, PtType);
		IPReallocNewTypeObject(PtObj, IP_OBJ_CTLPT);
		PtObj -> U.CtlPt = TmpObj -> U.CtlPt;
		IPFreeObject(TmpObj);
	    }
	}
    }

    return PtType;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Coerces an object to a new object.  Mostly about point types coercions.    M
*   Points, vectors, control points and planes can always be coerced between M
* themselves using this routine by specifying the new object type desired    M
* such as IP_OBJ_PLANE or control point type like CAGD_PT_E4_TYPE. 	     M
*   Control points of curves and surfaces may be coerced to a new type by    M
* prescribing the needed point type as NewType, such as CAGD_PT_P2_TYPE.     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:      Object to coerce.                                             M
*   NewType:   New type which can be object type like IP_OBJ_VECTOR or point M
*              type like E2.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   Newly coerced object.                                M
*                                                                            *
* SEE ALSO:                                                                  M
*   IPCoerceObjectTo                                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceObjectPtTypeTo, coercion                                         M
*****************************************************************************/
IPObjectStruct *IPCoerceObjectPtTypeTo(const IPObjectStruct *PObj, int NewType)
{
    int i;
    IrtRType Pt[CAGD_MAX_PT_SIZE];
    const IrtRType *R;
    IPObjectStruct
	*NewObj = NULL;

    switch (PObj -> ObjType) {
	case IP_OBJ_POINT:
	case IP_OBJ_VECTOR:
	case IP_OBJ_PLANE:
	case IP_OBJ_CTLPT:
	    if (PObj -> ObjType == NewType)
	        return IPCopyObject(NULL, PObj, FALSE);

	    Pt[0] = 1.0;
	    for (i = 1; i < CAGD_MAX_PT_SIZE; i++)
	        Pt[i] = 0.0;
	    switch (PObj -> ObjType) {
		case IP_OBJ_POINT:
		    IRIT_PT_COPY(&Pt[1], PObj -> U.Pt);
		    break;
		case IP_OBJ_VECTOR:
		    IRIT_PT_COPY(&Pt[1], PObj -> U.Vec);
		    break;
		case IP_OBJ_PLANE:
		    IRIT_PLANE_COPY(&Pt[1], PObj -> U.Plane);
		    break;
		case IP_OBJ_CTLPT:
		    R = PObj -> U.CtlPt.Coords;
		    CagdCoercePointTo(&Pt[1], CAGD_MAX_E_POINT,
				      (IrtRType * const *) &R,
				      -1, PObj -> U.CtlPt.PtType);
		    break;
		default:
		    break;
	    }

	    switch (NewType) {
		case IP_OBJ_POINT:
		    NewObj = IPGenPTObject(&Pt[1], &Pt[2], &Pt[3]);
		    break;
		case IP_OBJ_VECTOR:
		    NewObj = IPGenVECObject(&Pt[1], &Pt[2], &Pt[3]);
		    break;
		case IP_OBJ_PLANE:
		    NewObj = IPGenPLANEObject(&Pt[1], &Pt[2], &Pt[3], &Pt[4]);
		    break;
		case IP_OBJ_CTLPT:
		    NewObj = IPGenCTLPTObject(CAGD_PT_E3_TYPE, Pt);
		    break;
		CAGD_PT_P_TYPES
		    if (PObj -> ObjType == IP_OBJ_CTLPT)
		        CagdCoercePointTo(Pt, CAGD_MAX_P_POINT,
					  (IrtRType * const *) &R,
					  -1, PObj -> U.CtlPt.PtType);
		CAGD_PT_E_TYPES
		    NewObj = IPGenCTLPTObject((CagdPointType) NewType, Pt);
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
		default:
		    break;		    
	    }
	    break;
	case IP_OBJ_CURVE:
	    switch (NewType) {
		CAGD_PT_TYPES
		    NewObj = IPGenCRVObject(
				   CagdCoerceCrvsTo(PObj -> U.Crvs,
						    (CagdPointType) NewType,
						    TRUE));
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
		default:
		    break;
	    }
	    break;
	case IP_OBJ_SURFACE:
	    switch (NewType) {
		CAGD_PT_TYPES
		    NewObj = IPGenSRFObject(
				  CagdCoerceSrfsTo(PObj -> U.Srfs,
						   (CagdPointType) NewType,
						   TRUE));
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
	        case IP_OBJ_MODEL:
		    NewObj = IPGenMODELObject(MdlCnvrtSrf2Mdl(PObj -> U.Srfs));
		    break;
		default:
		    break;
	    }
	    break;
	case IP_OBJ_TRIVAR:
	    switch (NewType) {
		CAGD_PT_TYPES
		    NewObj = IPGenTRIVARObject(
				   TrivCoerceTVsTo(PObj -> U.Trivars,
						   (CagdPointType) NewType));
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
		default:
		    break;
	    }
	    break;
	case IP_OBJ_TRIMSRF:
	    switch (NewType) {
		CAGD_PT_TYPES
		    NewObj = IPGenTRIMSRFObject(TrimSrfNew(
			CagdCoerceSrfTo(PObj -> U.TrimSrfs -> Srf,
					(CagdPointType) NewType, FALSE),
			TrimCrvCopyList(PObj -> U.TrimSrfs -> TrimCrvList),
			TRUE));
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
	        case IP_OBJ_MODEL:
		    NewObj = IPGenMODELObject(
				 MdlCnvrtTrimmedSrf2Mdl(PObj -> U.TrimSrfs));
		    break;
		default:
		    break;
	    }
	    break;
	case IP_OBJ_TRISRF:
	    switch (NewType) {
		CAGD_PT_TYPES
		    NewObj = IPGenTRISRFObject(
				  TrngCoerceTriSrfsTo(PObj -> U.TriSrfs,
						   (CagdPointType) NewType));
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
		default:
		    break;
	    }
	    break;
	case IP_OBJ_MULTIVAR:
	    switch (NewType) {
		CAGD_PT_TYPES
		    NewObj = IPGenMULTIVARObject(
				   MvarCoerceMVsTo(PObj -> U.MultiVars,
						   (MvarPointType) NewType));
		    break;
	        case CAGD_PT_MAX_SIZE_TYPE:		    
		    /* Should not happend - might happen if CagdPointType  */
		    /* is modified in cagd_lib.h.			   */
		    assert(0);
		default:
		    break;
	    }
	    break;
	default:
	    break;
    }

    return NewObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Coerce an object to a new object.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:        Object to coerce.                                           M
*   NewType:     New type for PObj.                                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   The newly coerced object.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   IPCoerceObjectPtTypeTo                                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPCoerceObjectTo                                                         M
*****************************************************************************/
IPObjectStruct *IPCoerceObjectTo(const IPObjectStruct *PObj, int NewType)
{
    CagdSrfStruct *Srf;
    TrivTVStruct *TV;
    IPObjectStruct
        *PObjCopy = NULL,
	*NewObj = NULL;

    if (PObj -> ObjType == IP_OBJ_LIST_OBJ) {
	int i, j;
	IPObjectStruct *PObjTmp;

	NewObj = IPGenLISTObject(NULL);
	for (i = j = 0; (PObjTmp = IPListObjectGet(PObj, i)) != NULL; i++) {
	    IPObjectStruct
	        *PObjTmpNew = IPCoerceObjectTo(PObjTmp, NewType);

	    if (PObjTmpNew != NULL)
	        IPListObjectInsert(NewObj, j++, PObjTmpNew);
	}
	IPListObjectInsert(NewObj, j++, NULL);

	return NewObj;
    }

    switch (NewType) {
	case CAGD_CPOWER_TYPE:
	case CAGD_SPOWER_TYPE:
	case MVAR_POWER_TYPE:
	case IP_COERCE_POWER_TYPE:
	    if ((PObj -> ObjType == IP_OBJ_CURVE &&
		 PObj -> U.Crvs -> GType == CAGD_CBSPLINE_TYPE) ||
		(PObj -> ObjType == IP_OBJ_SURFACE &&
		 PObj -> U.Srfs -> GType == CAGD_SBSPLINE_TYPE) ||
		(PObj -> ObjType == IP_OBJ_MULTIVAR &&
		 PObj -> U.MultiVars -> GType == MVAR_BSPLINE_TYPE)) {
		PObj = PObjCopy = IPCoerceBsplineToBezier(PObj);
	    }
	    if ((PObj -> ObjType == IP_OBJ_CURVE &&
		 PObj -> U.Crvs -> GType == CAGD_CBEZIER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_SURFACE &&
		 PObj -> U.Srfs -> GType == CAGD_SBEZIER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_MULTIVAR &&
		 PObj -> U.MultiVars -> GType == MVAR_BEZIER_TYPE))
		NewObj = IPCoerceBezierToPower(PObj);
	    break;
	case CAGD_CBEZIER_TYPE:
	case CAGD_SBEZIER_TYPE:
	case TRIV_TVBEZIER_TYPE:
	case MVAR_BEZIER_TYPE:
	case IP_COERCE_BEZIER_TYPE:
	    if ((PObj -> ObjType == IP_OBJ_CURVE &&
		 PObj -> U.Crvs -> GType == CAGD_CPOWER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_SURFACE &&
		 PObj -> U.Srfs -> GType == CAGD_SPOWER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_MULTIVAR &&
		 PObj -> U.MultiVars -> GType == MVAR_POWER_TYPE))
		NewObj = IPCoercePowerToBezier(PObj);
	    else if (PObj -> ObjType == IP_OBJ_TRISRF &&
		     PObj -> U.TriSrfs -> GType == TRNG_TRISRF_GREGORY_TYPE)
		NewObj = IPCoerceGregoryToBezier(PObj);
	    else if ((PObj -> ObjType == IP_OBJ_CURVE &&
		      PObj -> U.Crvs -> GType == CAGD_CBSPLINE_TYPE) ||
		     (PObj -> ObjType == IP_OBJ_SURFACE &&
		      PObj -> U.Srfs -> GType == CAGD_SBSPLINE_TYPE) ||
		     (PObj -> ObjType == IP_OBJ_TRIVAR &&
		      PObj -> U.Trivars -> GType == TRIV_TVBSPLINE_TYPE) ||
		     (PObj -> ObjType == IP_OBJ_MULTIVAR &&
		      PObj -> U.MultiVars -> GType == MVAR_BSPLINE_TYPE))
	        NewObj = IPCoerceBsplineToBezier(PObj);
	    else if (PObj -> ObjType == IP_OBJ_TRIMSRF)
	        NewObj = IPCoerceTrimmedSrfToTrimmedBezier(PObj);
	    break;
	case CAGD_CBSPLINE_TYPE:
	case CAGD_SBSPLINE_TYPE:
	case TRIV_TVBSPLINE_TYPE:
	case MVAR_BSPLINE_TYPE:
	case IP_COERCE_BSPLINE_TYPE:
	    if ((PObj -> ObjType == IP_OBJ_CURVE &&
		 PObj -> U.Crvs -> GType == CAGD_CPOWER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_SURFACE &&
		 PObj -> U.Srfs -> GType == CAGD_SPOWER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_MULTIVAR &&
		 PObj -> U.MultiVars -> GType == MVAR_POWER_TYPE)) {
	        PObj = PObjCopy = IPCoercePowerToBezier(PObj);
	    }
	    if ((PObj -> ObjType == IP_OBJ_CURVE &&
		 PObj -> U.Crvs -> GType == CAGD_CBEZIER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_SURFACE &&
		 PObj -> U.Srfs -> GType == CAGD_SBEZIER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_TRIVAR &&
		 PObj -> U.Trivars -> GType == TRIV_TVBEZIER_TYPE) ||
		(PObj -> ObjType == IP_OBJ_MULTIVAR &&
		 PObj -> U.MultiVars -> GType == MVAR_BEZIER_TYPE))
		NewObj = IPCoerceBezierToBspline(PObj);
	    break;
	case IP_COERCE_UNTRIMMED_TYPE:
	    if (PObj -> ObjType == IP_OBJ_TRIMSRF) {
		int ComposeE3 = AttrGetObjectIntAttrib(PObj, "Untrim2D") != 1;

	        NewObj = IPCoerceTrimmedSrfToUnTrimmedBezier(PObj, ComposeE3);
	    }
	    break;
	case IP_COERCE_UNIFORM_PERIODIC:
	    IP_FATAL_ERROR(IP_ERR_CNVRT_TO_PERIODIC);
	    return NULL;
	case IP_COERCE_UNIFORM_FLOAT:
	    if (IP_IS_CRV_OBJ(PObj)) {
		if (CAGD_IS_PERIODIC_CRV(PObj -> U.Crvs)) {
		    NewObj =
		        IPGenCRVObject(CagdCnvrtPeriodic2FloatCrv(PObj ->
								     U.Crvs));
		}
		else if (CAGD_IS_BSPLINE_CRV(PObj -> U.Crvs) &&
			 !BspCrvHasOpenEC(PObj -> U.Crvs)) {
		    NewObj = IPGenCRVObject(CagdCrvCopy(PObj -> U.Crvs));
		}
		else {
		    IP_FATAL_ERROR(IP_ERR_CNVRT_PER_TO_FLOAT);
		    return NULL;
		}
	    }
	    else if (IP_IS_SRF_OBJ(PObj)) {
		if (CAGD_IS_PERIODIC_SRF(PObj -> U.Srfs)) {
		    Srf = CagdCnvrtPeriodic2FloatSrf(PObj -> U.Srfs);
		    NewObj = IPGenSRFObject(Srf);
		}
		else if (CAGD_IS_BSPLINE_SRF(PObj -> U.Srfs) &&
			 !BspSrfHasOpenEC(PObj -> U.Srfs)) {
		    NewObj = IPGenSRFObject(CagdSrfCopy(PObj -> U.Srfs));
		}
		else {
		    IP_FATAL_ERROR(IP_ERR_CNVRT_PER_TO_FLOAT);
		    return NULL;
		}
	    }
	    else if (IP_IS_TRIVAR_OBJ(PObj)) {
		if (TRIV_IS_PERIODIC_TV(PObj -> U.Trivars)) {
		    TV = TrivCnvrtPeriodic2FloatTV(PObj -> U.Trivars);
		    NewObj = IPGenTRIVARObject(TV);
		}
		else if (TRIV_IS_BSPLINE_TV(PObj -> U.Trivars) &&
			 !TrivBspTVHasOpenEC(PObj -> U.Trivars)) {
		    NewObj = IPGenTRIVARObject(TrivTVCopy(PObj -> U.Trivars));
		}
		else {
		    IP_FATAL_ERROR(IP_ERR_CNVRT_PER_TO_FLOAT);
		    return NULL;
		}
	    }
	    else if (IP_IS_MVAR_OBJ(PObj)) {
	        if (MVAR_IS_BSPLINE_MV(PObj -> U.MultiVars)) {
		    NewObj = IPGenMULTIVARObject(
			      MvarCnvrtPeriodic2FloatMV(PObj -> U.MultiVars));
		}
		else {
		    IP_FATAL_ERROR(IP_ERR_CNVRT_BSP_TO_FLOAT);
		    return NULL;
		}
	    }
	    break;
	case IP_COERCE_UNIFORM_OPEN:
	    if (IP_IS_CRV_OBJ(PObj)) {
		if (CAGD_IS_BEZIER_CRV(PObj -> U.Crvs) ||
		    (CAGD_IS_BSPLINE_CRV(PObj -> U.Crvs) &&
		     BspCrvHasOpenEC(PObj -> U.Crvs))) {
		    NewObj = IPGenCRVObject(CagdCrvCopy(PObj -> U.Crvs));
		}
		else {
		    NewObj = IPGenCRVObject(BspCrvOpenEnd(PObj -> U.Crvs));
		}
	    }
	    else if (IP_IS_SRF_OBJ(PObj)) {
		if (CAGD_IS_BEZIER_SRF(PObj -> U.Srfs) ||
		    (CAGD_IS_BSPLINE_SRF(PObj -> U.Srfs) &&
		     BspSrfHasOpenEC(PObj -> U.Srfs))) {
		    NewObj = IPGenSRFObject(CagdSrfCopy(PObj -> U.Srfs));
		}
		else {
		    NewObj = IPGenSRFObject(BspSrfOpenEnd(PObj -> U.Srfs));
		}
	    }
	    else if (IP_IS_TRIVAR_OBJ(PObj)) {
		if (TRIV_IS_BEZIER_TV(PObj -> U.Trivars) ||
		    (TRIV_IS_BSPLINE_TV(PObj -> U.Trivars) &&
		     TrivBspTVHasOpenEC(PObj -> U.Trivars))) {
		    NewObj = IPGenTRIVARObject(TrivTVCopy(PObj -> U.Trivars));
		}
		else {
		    NewObj = IPGenTRIVARObject(TrivTVOpenEnd(PObj -> U.Trivars));
		}
	    }
	    else if (IP_IS_MVAR_OBJ(PObj)) {
	        if (MVAR_IS_BEZIER_MV(PObj -> U.MultiVars)) {
		    NewObj = IPGenMULTIVARObject(
					     MvarMVCopy(PObj -> U.MultiVars));
		}
		else {
		    NewObj = IPGenMULTIVARObject(
			          MvarCnvrtFloat2OpenMV(PObj -> U.MultiVars));
		}
	    }
	    break;
	case IP_OBJ_CURVE:
	    if (IP_IS_MVAR_OBJ(PObj)) {
		if (PObj -> U.MultiVars -> Dim == 1)
		    NewObj = IPGenCRVObject(MvarMVToCrv(PObj -> U.MultiVars));
		else
		    IP_FATAL_ERROR(IP_ERR_CNVRT_MV_NOT_UNIVAR);
	    }
	    break;
	case IP_OBJ_SURFACE:
	    if (IP_IS_MVAR_OBJ(PObj)) {
		if (PObj -> U.MultiVars -> Dim == 2)
		    NewObj = IPGenSRFObject(MvarMVToSrf(PObj -> U.MultiVars));
		else
		    IP_FATAL_ERROR(IP_ERR_CNVRT_MV_NOT_BIVAR);
	    }
	    break;
	case IP_OBJ_TRIVAR:
	    if (IP_IS_MVAR_OBJ(PObj)) {
		if (PObj -> U.MultiVars -> Dim == 3)
		    NewObj = IPGenTRIVARObject(MvarMVToTV(PObj -> U.MultiVars));
		else
		    IP_FATAL_ERROR(IP_ERR_CNVRT_MV_NOT_TRIVAR);
	    }
	    break;
	case IP_OBJ_MODEL:
	    {
		MdlModelStruct
		    *Mdl = NULL;

		if (IP_IS_SRF_OBJ(PObj)) {
		    Mdl = MdlCnvrtSrf2Mdl(PObj -> U.Srfs);
		}
		else if (IP_IS_TRIMSRF_OBJ(PObj)) {
		    Mdl = MdlCnvrtTrimmedSrf2Mdl(PObj -> U.TrimSrfs);
		}
		else {
		    IP_FATAL_ERROR(IP_ERR_CNVRT_TSRF_TO_MDL);
		}

		NewObj = Mdl == NULL ? NULL : IPGenMODELObject(Mdl);
	    }
	    break;
	case IP_OBJ_TRIMSRF:
	    {
		TrimSrfStruct
		    *TSrf = NULL;
		IPObjectStruct *TrimSrfObj;

		if (IP_IS_SRF_OBJ(PObj)) {
		    TSrf = TrimSrfNew(PObj -> U.Srfs, NULL, FALSE);
		    NewObj = IPGenTRIMSRFObject(TSrf);
		}
		else if (IP_IS_MODEL_OBJ(PObj)) {
		    if ((TSrf = MdlCnvrtMdl2TrimmedSrfs(PObj -> U.Mdls)) != NULL) {
			int i = 0;

			NewObj = IPGenLISTObject(NULL);
			while (TSrf != NULL) {
			    TrimSrfObj = IPGenTRIMSRFObject(TSrf);
			    TSrf = TSrf -> Pnext;
			    TrimSrfObj -> U.TrimSrfs -> Pnext = NULL;
			    IPListObjectInsert(NewObj, i++, TrimSrfObj);
			}
			IPListObjectInsert(NewObj, i, NULL);
		    }		    
		}
		else {
		    IP_FATAL_ERROR(IP_ERR_CNVRT_SRF_MDL_TO_TSRF);
		}
	    }
	    break;
	case IP_OBJ_MULTIVAR:
	    switch (PObj -> ObjType) {
	        case IP_OBJ_CURVE:
		    NewObj = IPGenMULTIVARObject(MvarCrvToMV(PObj -> U.Crvs));
		    break;
		case IP_OBJ_SURFACE:
		    NewObj = IPGenMULTIVARObject(MvarSrfToMV(PObj -> U.Srfs));
		    break;
		case IP_OBJ_TRIVAR:
		    NewObj = IPGenMULTIVARObject(MvarTVToMV(PObj -> U.Trivars));
		    break;
		default:
		    IP_FATAL_ERROR(IP_ERR_CNVRT_INVALID_GEOM_TO_MV);
		    break;
	    }
	    break;
	default:
	    NewObj = IPCoerceObjectPtTypeTo(PObj, NewType);
	    break;
    }

    if (PObjCopy != NULL)
        IPFreeObject(PObjCopy);

    if (NewObj == NULL)
	IP_FATAL_ERROR(IP_ERR_CNVRT_INVALID_COERCE);
    return NewObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Routine to reverse a curve.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   CrvObj:     Curve to reverse its parametrization.                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:    Reversed curve.                                     *
*                                                                            *
* KEYWORDS:                                                                  *
*   CurveReverse                                                             *
*****************************************************************************/
static IPObjectStruct *CurveReverse(const IPObjectStruct *CrvObj)
{
    const CagdCrvStruct *Crv;
    CagdCrvStruct *RCrv,
	*RevCrvs = NULL;

    for (Crv = CrvObj -> U.Crvs; Crv != NULL; Crv = Crv -> Pnext) {
        RCrv = CagdCrvReverse(Crv);
	IRIT_LIST_PUSH(RCrv, RevCrvs);
    }

    if (RevCrvs == NULL)
	return NULL;

    return IPGenCRVObject(CagdListReverse(RevCrvs));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Routine to reverse a surface by reversing the U direction.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   SrfObj:     Surface to reverse.                                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:    Reversed surface. The normal of the reversed        *
*                        surface is flipped with respected to the original   *
*                        surface, SrfObj, by 180 degrees.		     *
*****************************************************************************/
static IPObjectStruct *SurfaceReverse(const IPObjectStruct *SrfObj)
{
    IPObjectStruct
	*RetVal = NULL;

    if (IP_IS_SRF_OBJ(SrfObj)) {
        const CagdSrfStruct *Srf;
	CagdSrfStruct *RSrf,
	    *RevSrfs = NULL;

	for (Srf = SrfObj -> U.Srfs; Srf != NULL; Srf = Srf -> Pnext) {
	    RSrf = CagdSrfReverse(Srf);
	    IRIT_LIST_PUSH(RSrf, RevSrfs);
	}

	if (RevSrfs != NULL)
	    RetVal = IPGenSRFObject(CagdListReverse(RevSrfs));
    }
    else if (IP_IS_TRIMSRF_OBJ(SrfObj)) {
        const TrimSrfStruct *Srf;
        TrimSrfStruct*RSrf,
	    *RevSrfs = NULL;

	for (Srf = SrfObj -> U.TrimSrfs; Srf != NULL; Srf = Srf -> Pnext) {
	    RSrf = TrimSrfReverse(Srf);
	    IRIT_LIST_PUSH(RSrf, RevSrfs);
	}

	if (RevSrfs != NULL)
	    RetVal = IPGenTRIMSRFObject(CagdListReverse(RevSrfs));
    }

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Routine to reverse a model.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   MdlObj:     Model to reverse inside out.	                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:    Reversed model.                                     *
*                                                                            *
* KEYWORDS:                                                                  *
*   CurveReverse                                                             *
*****************************************************************************/
static IPObjectStruct *ModelReverse(const IPObjectStruct *MdlObj)
{
    return IPGenMODELObject(MdlModelNegate(MdlObj -> U.Mdls));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns a similar hierarchy to given one but with reversed semantics.    M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:    Input object to reverse. 		                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:     Reversed object.                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPReverseObject                                                          M
*****************************************************************************/
IPObjectStruct *IPReverseObject(IPObjectStruct *PObj)
{
    IPObjectStruct *RetVal;

    switch (PObj ->ObjType) {
        case IP_OBJ_NUMERIC:
	    RetVal = IPGenNUMValObject(-PObj -> U.R);
	    break;
        case IP_OBJ_POINT:
	    RetVal = IPCopyObject(NULL, PObj, FALSE);
	    IRIT_PT_SCALE(RetVal -> U.Pt, -1);
	    break;
        case IP_OBJ_VECTOR:
	    RetVal = IPCopyObject(NULL, PObj, FALSE);
	    IRIT_VEC_SCALE(RetVal -> U.Vec, -1);
	    break;
        case IP_OBJ_CTLPT:
	    {
	        int i;
		CagdCtlPtStruct *CtlPt;

		RetVal = IPCopyObject(NULL, PObj, FALSE);
		CtlPt = &RetVal -> U.CtlPt;
		for (i = 1; i <= CAGD_NUM_OF_PT_COORD(CtlPt -> PtType); i++)
		    CtlPt -> Coords[i] = -CtlPt -> Coords[i];
	    }
	    break;
        case IP_OBJ_PLANE:
	    RetVal = IPCopyObject(NULL, PObj, FALSE);
	    IRIT_VEC_SCALE(RetVal -> U.Plane, -1);
	    break;
        case IP_OBJ_STRING:
	    {
	        int i, j;
		char *Str;

		RetVal = IPCopyObject(NULL, PObj, FALSE);
		Str = RetVal -> U.Str;
		for (i = 0, j = (int) strlen(Str) - 1; i < j; i++, j--)
		    IRIT_SWAP(char, Str[i], Str[j]);
	    }
	    break;
        case IP_OBJ_MATRIX:
	    {
	        IrtRType
		    R = -1.0;

		RetVal = IPAllocObject("", IP_OBJ_MATRIX, NULL);
		MatScale4by4(*RetVal -> U.Mat, *PObj -> U.Mat, &R);
	    }
	    break;
        case IP_OBJ_POLY:
	    if (IP_IS_POLYGON_OBJ(PObj)) {
	        RetVal = BooleanNEG(PObj);
	    }
	    else {
	        IPPolygonStruct *Pl;

		RetVal = IPCopyObject(NULL, PObj, FALSE);
		for (Pl = RetVal -> U.Pl; Pl != NULL; Pl = Pl -> Pnext)
		    Pl -> PVertex = IPReverseVrtxList2(Pl -> PVertex);
	    }
	    break;
         case IP_OBJ_CURVE:
	     RetVal = CurveReverse(PObj);
	    break;
         case IP_OBJ_SURFACE:
         case IP_OBJ_TRIMSRF:
	     RetVal = SurfaceReverse(PObj);
	     break;
         case IP_OBJ_MODEL:
	     RetVal = ModelReverse(PObj);
	     break;
         case IP_OBJ_LIST_OBJ:
	     {
	         int i = 0;
		 IPObjectStruct *PTmp;

		 RetVal = IPGenListObject(IP_GET_OBJ_NAME(PObj), NULL, NULL);
		 while ((PTmp = IPListObjectGet(PObj, i)) != NULL)
		     IPListObjectInsert(RetVal, i++, IPReverseObject(PTmp));
		 IPListObjectInsert(RetVal, i, NULL);
	     }
	     break;
         default:
	     RetVal = IPCopyObject(NULL, PObj, FALSE);
	     break;
    }

    IPCopyObjectAuxInfo(RetVal, PObj);
    IP_SET_OBJ_NAME2(RetVal, IP_GET_OBJ_NAME(PObj));

    return RetVal;
}
