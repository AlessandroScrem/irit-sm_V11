/******************************************************************************
* Bsp-Gen.c - Bspline generic routines.					      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Aug. 90.					      *
******************************************************************************/

#include "cagd_loc.h"

#define CAGD_MESH_CONT_LENRATIO_EPS	0.05
#define CAGD_MESH_CONT_ANGULAR_EPS	0.999

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new Bspline surface.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   ULength:      Number of control points in the U direction.               M
*   VLength:      Number of control points in the V direction.               M
*   UOrder:       The order of the surface in the U direction.               M
*   VOrder:       The order of the surface in the V direction.               M
*   PType:        Type of control points (E2, P3, etc.).                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   An uninitialized freeform Bspline surface.            M
*                                                                            *
* SEE ALSO:                                                                  M
*   BzrSrfNew, BspPeriodicSrfNew, CagdSrfNew, CagdPeriodicSrfNew, TrimSrfNew M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfNew, allocation                                                    M
*****************************************************************************/
CagdSrfStruct *BspSrfNew(int ULength,
			 int VLength,
			 int UOrder,
			 int VOrder,
			 CagdPointType PType)
{
    CagdSrfStruct *Srf;

    if (ULength < UOrder || VLength < VOrder) {
	CAGD_FATAL_ERROR(CAGD_ERR_WRONG_ORDER);
	return NULL;
    }

    Srf = CagdSrfNew(CAGD_SBSPLINE_TYPE, PType, ULength, VLength);

    Srf -> UKnotVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
							   (UOrder + ULength));
    Srf -> VKnotVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
							   (VOrder + VLength));

    Srf -> UOrder = UOrder;
    Srf -> VOrder = VOrder;

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new, possibly periodic, Bspline      M
* surface.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   ULength:      Number of control points in the U direction.               M
*   VLength:      Number of control points in the V direction.               M
*   UOrder:       The order of the surface in the U direction.               M
*   VOrder:       The order of the surface in the V direction.               M
*   UPeriodic:    Is this surface periodic in the U direction?               M
*   VPeriodic:    Is this surface periodic in the V direction?               M
*   PType:        Type of control points (E2, P3, etc.).                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   An uninitialized freeform Bspline surface. If both    M
*                      UPeriodic and VPeriodic are FALSE, this function is   M
*                      identical to BspSrfNew.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspSrfNew, BzrSrfNew, CagdSrfNew, CagdPeriodicSrfNew, TrimSrfNew         M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspPeriodicSrfNew, allocation                                            M
*****************************************************************************/
CagdSrfStruct *BspPeriodicSrfNew(int ULength,
				 int VLength,
				 int UOrder,
				 int VOrder,
				 CagdBType UPeriodic,
				 CagdBType VPeriodic,
				 CagdPointType PType)
{
    CagdSrfStruct *Srf;

    if (ULength < UOrder || VLength < VOrder) {
	CAGD_FATAL_ERROR(CAGD_ERR_WRONG_ORDER);
	return NULL;
    }

    Srf = CagdPeriodicSrfNew(CAGD_SBSPLINE_TYPE, PType, ULength, VLength,
			     UPeriodic, VPeriodic);

    Srf -> UKnotVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
			   (UOrder + ULength + (UPeriodic ? UOrder - 1 : 0)));
    Srf -> VKnotVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
			   (VOrder + VLength + (VPeriodic ? VOrder - 1 : 0)));

    Srf -> UOrder = UOrder;
    Srf -> VOrder = VOrder;

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new Bspline curve.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Length:     Number of control points                                     M
*   Order:      The order of the curve                                       M
*   PType:      Type of control points (E2, P3, etc.).                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   An uninitialized freeform Bspline curve.              M
*                                                                            *
* SEE ALSO:                                                                  M
*   BzrCrvNew, BspPeriodicCrvNew, CagdCrvNew, CagdPeriodicCrvNew, TrimCrvNew M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvNew, allocation                                                    M
*****************************************************************************/
CagdCrvStruct *BspCrvNew(int Length, int Order, CagdPointType PType)
{
    CagdCrvStruct *Crv;

    if (Length < Order) {
	CAGD_FATAL_ERROR(CAGD_ERR_WRONG_ORDER);
	return NULL;
    }

    Crv = CagdCrvNew(CAGD_CBSPLINE_TYPE, PType, Length);

    Crv -> KnotVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
							     (Order + Length));

    Crv -> Order = Order;

    if (Crv -> Order == 2)
	CAGD_SET_GEOM_TYPE(Crv, CAGD_GEOM_LINEAR);
    else if (Crv -> Order == 1)
	CAGD_SET_GEOM_TYPE(Crv, CAGD_GEOM_CONST);

    return Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates the memory required for a new, possibly periodic, Bspline      M
* curve.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   Length:     Number of control points                                     M
*   Order:      The order of the curve                                       M
*   Periodic:   Is this curve periodic?                                      M
*   PType:      Type of control points (E2, P3, etc.).                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   An uninitialized freeform Bspline curve. If Periodic  M
*                      is FALSE, this function is identical to BspCrvNew.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   BzrCrvNew, BspCrvNew, CagdCrvNew, CagdPeriodicCrvNew, TrimCrvNew	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspPeriodicCrvNew, allocation                                            M
*****************************************************************************/
CagdCrvStruct *BspPeriodicCrvNew(int Length,
				 int Order,
				 CagdBType Periodic,
				 CagdPointType PType)
{
    CagdCrvStruct *Crv;

    if (Length < Order) {
	CAGD_FATAL_ERROR(CAGD_ERR_WRONG_ORDER);
	return NULL;
    }

    Crv = CagdPeriodicCrvNew(CAGD_CBSPLINE_TYPE, PType, Length, Periodic);

    Crv -> KnotVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
				(Order + Length + (Periodic ? Order - 1 : 0)));

    Crv -> Order = Order;

    return Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the parametric domain of a Bspline curve.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:       To get its parametric domain.                                 M
*   TMin:      Where to put the minimal domain's boundary.                   M
*   TMax:      Where to put the maximal domain's boundary.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCrvDomain                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvDomain, domain, parametric domain                                  M
*****************************************************************************/
void BspCrvDomain(const CagdCrvStruct *Crv, CagdRType *TMin, CagdRType *TMax)
{
    int k = Crv -> Order,
	Len = CAGD_CRV_PT_LST_LEN(Crv);

    *TMin = Crv -> KnotVector[k - 1];
    *TMax = Crv -> KnotVector[Len];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the parametric domain of a Bspline surface.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:       To get its parametric domain.                                 M
*   UMin:      Where to put the minimal U domain's boundary.                 M
*   UMax:      Where to put the maximal U domain's boundary.                 M
*   VMin:      Where to put the minimal V domain's boundary.                 M
*   VMax:      Where to put the maximal V domain's boundary.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdSrfDomain, TrimSrfDomain                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfDomain, domain, parametric domain                                  M
*****************************************************************************/
void BspSrfDomain(const CagdSrfStruct *Srf,
		  CagdRType *UMin,
		  CagdRType *UMax,
		  CagdRType *VMin,
		  CagdRType *VMax)
{
    int UOrder = Srf -> UOrder,
	VOrder = Srf -> VOrder,
	ULen = CAGD_SRF_UPT_LST_LEN(Srf),
	VLen = CAGD_SRF_VPT_LST_LEN(Srf);

    *UMin = Srf -> UKnotVector[UOrder - 1];
    *UMax = Srf -> UKnotVector[ULen];
    *VMin = Srf -> VKnotVector[VOrder - 1];
    *VMax = Srf -> VKnotVector[VLen];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns a curve with open end conditions, similar to given curve.	     M
*   Open end curve is computed by extracting a subregion from Crv that is    M
* the entire curve's parametric domain, by inserting multiple knots at the   M
* domain's boundary.                                                         M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:     To convert to a new curve with open end conditions.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *: Same curve as Crv but with open end conditions.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspSrfOpenEnd                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvOpenEnd, open end conditions                                       M
*****************************************************************************/
CagdCrvStruct *BspCrvOpenEnd(const CagdCrvStruct *Crv)
{
    CagdRType TMin, TMax;
    CagdCrvStruct *OpenCrv, *TCrv;

    CagdCrvDomain(Crv, &TMin, &TMax);

    if (CAGD_IS_PERIODIC_CRV(Crv)) {
        TCrv = CagdCnvrtPeriodic2FloatCrv(Crv);
	OpenCrv = CagdCrvRegionFromCrv(TCrv, TMin, TMax);
	CagdCrvFree(TCrv);
    }
    else
        OpenCrv = CagdCrvRegionFromCrv(Crv, TMin, TMax);

    CAGD_PROPAGATE_ATTR(OpenCrv, Crv);

    return OpenCrv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns a surface with open end conditions, similar to given surface.    M
*   Open end surface is computed by extracting a subregion from Srf that is  M
* the entire surface's parametric domain, by inserting multiple knots at the M
* domain's boundary.                                                         M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:     To convert to a new surface with open end conditions.  Input    M
*	     can also be periodic.				             M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *: Same surface as Srf but with open end conditions.       M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvOpenEnd                                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfOpenEnd, open end conditions                                       M
*****************************************************************************/
CagdSrfStruct *BspSrfOpenEnd(const CagdSrfStruct *Srf)
{
    CagdRType UMin, UMax, VMin, VMax;
    CagdSrfStruct *TSrf, *TSrf2, *OpenSrf;

    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

    if (CAGD_IS_PERIODIC_SRF(Srf)) {
        TSrf2 = CagdCnvrtPeriodic2FloatSrf(Srf);
	TSrf = CagdSrfRegionFromSrf(TSrf2, UMin, UMax, CAGD_CONST_U_DIR);
	CagdSrfFree(TSrf2);
    }
    else
        TSrf = CagdSrfRegionFromSrf(Srf, UMin, UMax, CAGD_CONST_U_DIR);

    OpenSrf = CagdSrfRegionFromSrf(TSrf, VMin, VMax, CAGD_CONST_V_DIR);

    CagdSrfFree(TSrf);

    CAGD_PROPAGATE_ATTR(OpenSrf, Srf);

    return OpenSrf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Scans the given knot vector of the given curve for a potential C0        M
* discontinuity.							     M
*   Looks for multiplicities in the knot sequence and then examine the mesh  M
* if indeed the mesh is discontinuous at that location.			     M
*   Assumes knot vectors has open end condition.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:       To examine its potential discontinuity.			     M
*   t:         Where to put the parameter value (knot) that can be C0        M
*              discontinuous.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:     TRUE if found a C0 discontinuity, FALSE otherwise.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvKnotC1Discont, BspCrvKnotC2Discont, BspKnotC0Discont,		     M
*   BspCrvMeshC1Continuous, BspCrvKnotC1Discont			             M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvKnotC0Discont, knot vectors, continuity, discontinuity             M
*****************************************************************************/
CagdBType BspCrvKnotC0Discont(const CagdCrvStruct *Crv, CagdRType *t)
{
    if (CAGD_IS_BEZIER_CRV(Crv) || CAGD_IS_POWER_CRV(Crv))
	return FALSE;
    else if (CAGD_IS_BSPLINE_CRV(Crv))
	return BspKnotC0Discont(Crv -> KnotVector, Crv -> Order, Crv -> Length, t);
    else
	return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Scans the given knot vector of the given curve for a potential C1        M
* discontinuity.							     M
*   Looks for multiplicities in the knot sequence and then examine the mesh  M
* if indeed the mesh is discontinuous at that location.			     M
*   Assumes knot vectors has open end condition.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:       To examine its potential discontinuity.			     M
*   t:         Where to put the parameter value (knot) that can be C1        M
*              discontinuous.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:     TRUE if found a C1 discontinuity, FALSE otherwise.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvKnotC0Discont, BspCrvKnotC2Discont, BspKnotC1Discont,              M
*   BspCrvMeshC1Continuous					             M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvKnotC1Discont, knot vectors, continuity, discontinuity             M
*****************************************************************************/
CagdBType BspCrvKnotC1Discont(const CagdCrvStruct *Crv, CagdRType *t)
{
    if (CAGD_IS_BEZIER_CRV(Crv) || CAGD_IS_POWER_CRV(Crv))
	return FALSE;
    else if (CAGD_IS_BSPLINE_CRV(Crv))
	return BspKnotC1Discont(Crv -> KnotVector, Crv -> Order, Crv -> Length, t);
    else
	return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Scans the given knot vector of the given curve for a potential C2        M
* discontinuity.							     M
*   Looks for multiplicities in the knot sequence and then examine the mesh  M
* if indeed the mesh is discontinuous at that location.			     M
*   Assumes knot vectors has open end condition.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:       To examine its potential discontinuity.			     M
*   t:         Where to put the parameter value (knot) that can be C1        M
*              discontinuous.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:     TRUE if found a C2 discontinuity, FALSE otherwise.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvKnotC0Discont, BspCrvKnotC1Discont, BspKnotC1Discont,              M
*   BspCrvMeshC1Continuous					             M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvKnotC2Discont, knot vectors, continuity, discontinuity             M
*****************************************************************************/
CagdBType BspCrvKnotC2Discont(const CagdCrvStruct *Crv, CagdRType *t)
{
    if (CAGD_IS_BEZIER_CRV(Crv) || CAGD_IS_POWER_CRV(Crv))
	return FALSE;
    else if (CAGD_IS_BSPLINE_CRV(Crv))
        return BspKnotC2Discont(Crv -> KnotVector, Crv -> Order, Crv -> Length, t);
    else
	return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Examine the control polygon of the given curve in index Idx for a real   M
* C1 discontinuity in the mesh.				                     M
*   This index will typically be for a knot multiplicity potential discont.  M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:       To examine its potential discontinuity.		             M
*   Idx:       Index where to examine the discontinuity.                     M
*   CosAngle:  If not NULL, updated with the cosine of the deviation angle.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:    TRUE if continuous there, FALSE otherwise.                 M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspKnotC1Discont                                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvMeshC1Continuous                                                   M
*****************************************************************************/
CagdBType BspCrvMeshC1Continuous(const CagdCrvStruct *Crv,
				 int Idx,
				 CagdRType *CosAngle)
{
    CagdRType * const
        *Points = Crv -> Points;
    CagdRType DotProd;
    CagdPType Pt0, Pt1, Pt2;
    CagdVType V1, V2;
    CagdPointType
	PType = Crv -> PType;

    CagdCoerceToE3(Pt0, Points, Idx - 1, PType);
    CagdCoerceToE3(Pt1, Points, Idx,     PType);
    CagdCoerceToE3(Pt2, Points, Idx + 1, PType);

    IRIT_PT_SUB(V1, Pt0, Pt1);
    IRIT_PT_SUB(V2, Pt1, Pt2);

    IRIT_VEC_SAFE_NORMALIZE(V1);
    IRIT_VEC_SAFE_NORMALIZE(V2);

    DotProd = IRIT_DOT_PROD(V1, V2);
    if (CosAngle != NULL)
	*CosAngle = DotProd;
    return DotProd > CAGD_MESH_CONT_ANGULAR_EPS;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Scans the given knot vector of the given surface for a potential C0        M
* discontinuity.							     M
*   Looks for multiplicities in the knot sequence and then examine the mesh  M
* if indeed the mesh is discontinuous at that location.			     M
*   Assumes knot vectors has open end condition.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:       To examine its potential discontinuity across Dir.            M
*   Dir:       Direction to examine the discontinuity across.                M
*   t:         Where to put the parameter value (knot) that can be C0        M
*              discontinuous.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:     TRUE if found a C0 discontinuity, FALSE otherwise.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspSrfKnotC1Discont, BspKnotC1Discont, BspSrfMeshC1Continuous            M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfKnotC0Discont, knot vectors, continuity, discontinuity             M
*****************************************************************************/
CagdBType BspSrfKnotC0Discont(const CagdSrfStruct *Srf,
			      CagdSrfDirType Dir,
			      CagdRType *t)
{
    int Order = Dir == CAGD_CONST_U_DIR ? Srf -> UOrder : Srf -> VOrder,
	Length = Dir == CAGD_CONST_U_DIR ? Srf -> ULength : Srf -> VLength;
    CagdRType
 	*KV = Dir == CAGD_CONST_U_DIR ? Srf -> UKnotVector
				      : Srf -> VKnotVector;

    return BspKnotC0Discont(KV, Order, Length, t);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Scans the given knot vector of the given surface for a potential C1        M
* discontinuity.							     M
*   Looks for multiplicities in the knot sequence and then examine the mesh  M
* if indeed the mesh is discontinuous at that location.			     M
*   Assumes knot vectors has open end condition.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:       To examine its potential discontinuity across Dir.            M
*   Dir:       Direction to examine the discontinuity across.                M
*   t:         Where to put the parameter value (knot) that can be C1        M
*              discontinuous.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:     TRUE if found a C1 discontinuity, FALSE otherwise.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspSrfKnotC0Discont, BspKnotC1Discont, BspSrfMeshC1Continuous            M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfKnotC1Discont, knot vectors, continuity, discontinuity             M
*****************************************************************************/
CagdBType BspSrfKnotC1Discont(const CagdSrfStruct *Srf,
			      CagdSrfDirType Dir,
			      CagdRType *t)
{
    int Order = Dir == CAGD_CONST_U_DIR ? Srf -> UOrder : Srf -> VOrder,
	Length = Dir == CAGD_CONST_U_DIR ? Srf -> ULength : Srf -> VLength;
    CagdRType
 	*KV = Dir == CAGD_CONST_U_DIR ? Srf -> UKnotVector
				      : Srf -> VKnotVector;

    return BspKnotC1Discont(KV, Order, Length, t);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Examine the mesh of the given surface across direction Dir in index      M
* of mesh Index for a real discontinuity in the mesh.                        M
*   This index will typically be for a knot multiplicity potential discont.  M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:       To examine its potential discontinuity across Dir.            M
*   Dir:       Direction to examine the discontiinuty across.                M
*   Idx:       Index where to examine the discontinuity.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:    TRUE if continuous there, FALSE otherwise.                 M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspKnotC1Discont, BspSrfIsC1DiscontAt                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfMeshC1Continuous                                                   M
*****************************************************************************/
CagdBType BspSrfMeshC1Continuous(const CagdSrfStruct *Srf,
				 CagdSrfDirType Dir,
				 int Idx)
{
    int i,
	ULength = Srf -> ULength,
	VLength = Srf -> VLength;
    CagdRType Len1, Len2, t,
	LenRatio = IRIT_INFNTY;
    CagdRType
	* const *Pts = Srf -> Points;
    CagdPType Pt0, Pt1, Pt2;
    CagdVType V1, V2;
    CagdPointType
	PType = Srf -> PType;

    switch (Dir) {
	case CAGD_CONST_U_DIR:
            for (i = 0; i < VLength; i++) {
	        CagdCoerceToE3(Pt0, Pts, CAGD_MESH_UV(Srf, Idx - 1, i), PType);
	        CagdCoerceToE3(Pt1, Pts, CAGD_MESH_UV(Srf, Idx, i),     PType);
	        CagdCoerceToE3(Pt2, Pts, CAGD_MESH_UV(Srf, Idx + 1, i), PType);

		IRIT_PT_SUB(V1, Pt0, Pt1);
		IRIT_PT_SUB(V2, Pt1, Pt2);

		Len1 = IRIT_PT_LENGTH(V1);
		Len2 = IRIT_PT_LENGTH(V2);

		if (Len1 < IRIT_EPS && Len2 < IRIT_EPS)
		    continue;

		if (LenRatio == IRIT_INFNTY && Len1 != 0 && Len2 != 0)
		    LenRatio = Len1 / Len2;
		else {
		    t = Len2 == 0 ? (Len1 == 0 ? LenRatio : IRIT_INFNTY)
				  : Len1 / Len2;
		    if (!IRIT_APX_EQ_EPS(LenRatio, t, CAGD_MESH_CONT_LENRATIO_EPS))
		        return FALSE;
		}

		if (Len1 > 0 && Len2 > 0) {
		    Len1 = 1.0 / Len1;
		    Len2 = 1.0 / Len2;
		    IRIT_PT_SCALE(V1, Len1);
		    IRIT_PT_SCALE(V2, Len2);

		    if (IRIT_DOT_PROD(V1, V2) < CAGD_MESH_CONT_ANGULAR_EPS)
		        return FALSE;
		}
	    }
	    break;
	case CAGD_CONST_V_DIR:
            for (i = 0; i < ULength; i++) {
	        CagdCoerceToE3(Pt0, Pts, CAGD_MESH_UV(Srf, i, Idx - 1), PType);
	        CagdCoerceToE3(Pt1, Pts, CAGD_MESH_UV(Srf, i, Idx),     PType);
	        CagdCoerceToE3(Pt2, Pts, CAGD_MESH_UV(Srf, i, Idx + 1), PType);

		IRIT_PT_SUB(V1, Pt0, Pt1);
		IRIT_PT_SUB(V2, Pt1, Pt2);

		Len1 = IRIT_PT_LENGTH(V1);
		Len2 = IRIT_PT_LENGTH(V2);

		if (Len1 < IRIT_EPS && Len2 < IRIT_EPS)
		    continue;

		if (LenRatio == IRIT_INFNTY && Len1 != 0 && Len2 != 0)
		    LenRatio = Len1 / Len2;
		else {
		    t = Len2 == 0 ? (Len1 == 0 ? LenRatio : IRIT_INFNTY)
				  : Len1 / Len2;
		    if (!IRIT_APX_EQ_EPS(LenRatio, t, CAGD_MESH_CONT_LENRATIO_EPS))
		        return FALSE;
		}

		if (Len1 > 0 && Len2 > 0) {
		    Len1 = 1.0 / Len1;
		    Len2 = 1.0 / Len2;
		    IRIT_PT_SCALE(V1, Len1);
		    IRIT_PT_SCALE(V2, Len2);

		    if (IRIT_DOT_PROD(V1, V2) < CAGD_MESH_CONT_ANGULAR_EPS)
		        return FALSE;
		}
	    }
	    break;
	default:
	    CAGD_FATAL_ERROR(CAGD_ERR_DIR_NOT_CONST_UV);
	    break;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Extends a B-spline curve. The domain of Crv is extended, such that the   M
* trace coincides with the input curve's trace over the original domain.     M
*   An interface function for the one-sided curve extension function.        M
*									     *
* PARAMETERS:                                                                M
*   OrigCrv:  The curve to be extended.					     M
*   ExtDirs:  A boolean array of size 2 to determine the required directions M
	      of extension  MinDmn, MaxDmn. A NULL here means (TRUE, TRUE).  M
*   Epsilon:  The length of the requested extension,  in the param. domain.  M
*   RemoveExtraKnots: If FALSE, the resulting curve will not have minimal    M
*	      multiplicity at the first internal knot on the extension side. M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    The new extended curve.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvExtraKnotRmv, BspSrfExtension, BspCrvExtensionOneSide              M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvExtension							     M
*****************************************************************************/
CagdCrvStruct *BspCrvExtension(const CagdCrvStruct *OrigCrv,
			       const CagdBType *ExtDirs,
			       CagdRType Epsilon,
			       CagdBType RemoveExtraKnots)
{
    CagdBType 
	BothDirs = ExtDirs == NULL;
    CagdCrvStruct *CTmp,
        *Crv = CagdCrvCopy(OrigCrv);

    if (BothDirs || ExtDirs[0]) {
	CTmp = BspCrvExtensionOneSide(Crv, TRUE,
				      Epsilon, RemoveExtraKnots);
	CagdCrvFree(Crv);
	Crv = CTmp;
    }
    if (BothDirs || ExtDirs[1]) {
	CTmp = BspCrvExtensionOneSide(Crv, FALSE,
				      Epsilon, RemoveExtraKnots);
	CagdCrvFree(Crv);
	Crv = CTmp;
    }
    return Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Extends a B-spline curve, at the min/max end. The domain of Crv is       M
* extended, such that the trace coincides with the original trace over the   M
* original domain. Assumes Crv has open end conditions.			     M
*									     *
* PARAMETERS:                                                                M
*   OrigCrv:  The curve to be extended.					     M
*   MinDmn:   TRUE for min domain extension, FALSE for max domain extension. M
*   Epsilon:  The length of the extension in the domain.		     M
*   RemoveExtraKnots: If FALSE, the resulting curve will not have minimal    M
*		multiplicity at the first internal knot on the extension     M
*		side.							     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    The new extended curve.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvExtraKnotRmv, BspSrfExtension	                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvExtensionOneSide						     M
*****************************************************************************/
CagdCrvStruct *BspCrvExtensionOneSide(const CagdCrvStruct *OrigCrv,
				      CagdBType MinDmn,
				      CagdRType Epsilon,
				      CagdBType RemoveExtraKnots)
{
    CagdBType
	IsRational = CAGD_IS_RATIONAL_CRV(OrigCrv);
    int i, j, InsertionsNum, MaxAxis, RmvInd,
	Mult = 1,
	/* Easier notation to match the indices in the algo. we use:        */
	k = OrigCrv -> Order - 1,	                /* k is the degree. */
	n = OrigCrv -> Length - 1,   /* n + 1 is number of control points.  */
	N = n + k + 1;		           /* N + 1 is the number of knots. */
    CagdRType Knot2Insert, TmpMaxDmnEnd, TmpMinDmnEnd, t, **ExtensionPoints;
    CagdCrvStruct *TmpBzr, *MinDmnCrv, *MaxDmnCrv, *Crv;

    if (! CAGD_IS_BSPLINE_CRV(OrigCrv)) {
        CAGD_FATAL_ERROR(CAGD_ERR_BZR_CRV_EXPECT);
	return NULL;
    }

    MaxAxis = CAGD_NUM_OF_PT_COORD(OrigCrv -> PType);
    /* Assume open end conditions. Find the multiplicity at the required    */
    /* knot we are to manipulate.                                           */
    if (MinDmn) {
	Knot2Insert = OrigCrv -> KnotVector[k + 1];
	Mult = BspKnotFindMult(OrigCrv -> KnotVector, OrigCrv -> Order, 
			       OrigCrv -> Length, Knot2Insert);
    }
    else {
	Knot2Insert = OrigCrv -> KnotVector[N - k - 1];
	i = N - k - 2;
	Mult = BspKnotFindMult(OrigCrv -> KnotVector, OrigCrv -> Order, 
			       OrigCrv ->Length, Knot2Insert);
    }
    /* Insert until the last segment is Bezier like: no internal knots and  */
    /* interpolating the end points: multiplicity k suffices at the         */
    /* required knot, one knot less than the actual end.                    */
    InsertionsNum = k - Mult;
    Crv = BspCrvKnotInsertNSame(OrigCrv, Knot2Insert, InsertionsNum);
    n = Crv -> Length - 1;	           /* Back to updated also indices. */
    N = n + k + 1;	

    /* Create the corresponding Bezier and perform Bezier extension.        */
    TmpBzr = BzrCrvNew(k + 1, Crv -> PType);
    if (MinDmn) {
	/* Control points are the first k + 1 */
	for (j = !IsRational; j <= MaxAxis; j++) {
	    for (i = 0; i < k + 1; i++) {
		TmpBzr -> Points[j][i] = Crv -> Points[j][i];
	    }
	}

	/* The domain of the (generalized) Bezier piece: */
	TmpMinDmnEnd = Crv -> KnotVector[0];
	t = TmpMinDmnEnd - Epsilon;

	/* Correct domain extension of the original curve: */
	for (i = 0; i < k + 1; i++) 
	    Crv -> KnotVector[i] = t;

	/* Remember that we use the generalized form (t - a) / (b - a): */
	t = (t - TmpMinDmnEnd) / (Knot2Insert - TmpMinDmnEnd);
    }
    else {
	/* Control points are the last k + 1 */
	for (j = !IsRational; j <= MaxAxis; j++) {
	    for (i = 0; i < k + 1; i++) {
		TmpBzr -> Points[j][k - i] = Crv -> Points[j][n - i];
	    }
	}
	TmpMaxDmnEnd = Crv -> KnotVector[N];
	t = TmpMaxDmnEnd + Epsilon;

	/* Correct domain extension of the original curve: */
	for (i = 0; i < k + 1; i++) 
	    Crv -> KnotVector[N - i] = t;

	/* Remember that we use the generalized form (t - a) / (b - a): */
	t = (t - Knot2Insert) / (TmpMaxDmnEnd - Knot2Insert);
    }
    MinDmnCrv = BzrCrvSubdivAtParam(TmpBzr, t);
    MaxDmnCrv = MinDmnCrv -> Pnext;

    /* We now update the control points of the original curve. */ 
    if (MinDmn) {
	ExtensionPoints = MaxDmnCrv -> Points;
	for (j = !IsRational; j <= MaxAxis; j++) {
	    CAGD_GEN_COPY(Crv -> Points[j], ExtensionPoints[j],
			  sizeof(CagdRType) * (k + 1));
	}
    }
    else {
	ExtensionPoints = MinDmnCrv -> Points;
	for (j = !IsRational; j <= MaxAxis; j++) {
	    CAGD_GEN_COPY(Crv -> Points[j] + n - k, ExtensionPoints[j], 
			  sizeof(CagdRType) * (k + 1));
	}
    }

    /* And now remove the extra knots: */
    if (RemoveExtraKnots) {
	CagdCrvStruct *TCrv;

	if (MinDmn) {
	    RmvInd = k + 1;
	    for (i = 0; i < InsertionsNum; i++) { 
		/* Removing from the min. end. */
		TCrv = BspCrvExtraKnotRmv(Crv, RmvInd);
		CagdCrvFree(Crv);
		Crv = TCrv;
	    }
	}
	else {
	    RmvInd = N - k - 1;
	    for (i = 0; i < InsertionsNum; i++) {
		/* Removing from the max. end, remember that the removal   */
		/* index is always N - k - 1, but N decreases!             */
		TCrv = BspCrvExtraKnotRmv(Crv, RmvInd--);
		CagdCrvFree(Crv);
		Crv = TCrv;
	    }
	}
    }
    CagdCrvFreeList(MinDmnCrv);
    CagdCrvFree(TmpBzr);

    return Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reverse operation of the knot insertion, assuming it is guaranteed that  M
*   the knot currently does not have its minimal multiplicity, that is:	     M
*   it is possible to remove it and maintain the exact trace.		     M
*   Indices and conventions follow Boehm's Knot Insertion algo, as in	     M
*   E. Cohen, R.F. Riesenfeld, G. Elber, "Geometric Modeling with Splines:   M
*   an Introduction", CH 07.						     M
*									     *
* PARAMETERS:                                                                M
*   Crv:       The curve to be updated.					     M
*   RmvIndex:  The index in the knot vector of the knot to be removed.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    The new updated curve.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvExtensionOneSide, BspSrfExtension	                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspCrvExtraKnotRmv							     M
*****************************************************************************/
CagdCrvStruct *BspCrvExtraKnotRmv(const CagdCrvStruct *Crv,
				  int RmvIndex)
{
    CagdBType
	IsRational = CAGD_IS_RATIONAL_CRV(Crv);
    int i, d, MaxAxis, 
	/* Easier notation to match the indices in the algo we use:         */
	k = Crv -> Order - 1,                   	/* k is the degree. */
	n = Crv -> Length - 1,/* n + 1 = number of points. (before removal).*/
	N = n + k + 1,	      /* N + 1 = number of knots. (before removal). */
	J = 0;
    CagdRType Enumerator, Denominator,
	Knot2Remove = Crv -> KnotVector[RmvIndex];
    CagdRType *RmvCoef;
    CagdCrvStruct *CrvRmvd; 

    /* Create output curve: Length decreased by one, same order, update KV. */
    CrvRmvd = BspCrvNew(n, k + 1, Crv -> PType);
    CAGD_GEN_COPY(CrvRmvd -> KnotVector, Crv -> KnotVector, 
		  sizeof(CagdRType) * RmvIndex);
    CAGD_GEN_COPY(CrvRmvd -> KnotVector + RmvIndex, 
		  Crv -> KnotVector + RmvIndex + 1,
		  sizeof(CagdRType) * (N - RmvIndex));

    /* Essentially we reverse Boehm's Knot Insertion. We need J, so we       */
    /* find it by pretending we are about to insert the knot we're removing. */
    /* (this is why we call now with N and not N + 1, etc.)*/
    J = BspKnotLastIndexLE(CrvRmvd -> KnotVector, N, Knot2Remove); 

    /* The removal coefficients are the ones appearing in Boehm's Knot       */
    /* Insertion, with the exception that we don't need the first one,       */
    /* under the assumption that exact removal is possible. It is an over    */
    /* constrained linear system, that we later solve. */
    RmvCoef = (CagdRType *) IritMalloc((k - 1) * sizeof(CagdRType));
    Denominator = CrvRmvd -> KnotVector[J + 1] - 
	          CrvRmvd -> KnotVector[J - k + 1];
    Enumerator = Knot2Remove - CrvRmvd -> KnotVector[J - k + 1];
    if ((IRIT_ABS(Denominator) < IRIT_EPS) || 
	(IRIT_ABS(1 - Enumerator / Denominator)) < IRIT_EPS) {
        CAGD_FATAL_ERROR(CAGD_ERR_NO_EXTRA_KNOTS);
	IritFree(RmvCoef);
	return NULL;
    }

    /* We now construct the coeff's of the removal matrix: a stripe matrix, */
    /* with 1 - RmvCoef[i] in its i'th diagonal entry, and RmvCoef[i] above */
    /* the diagonal. It is a k - 1 x k - 1 submatrix of the original matrix */
    /* in Boehm's algorithm, obtained by removing the first row, first      */
    /* column and last column.						    */
    for (i = 0; i < k - 1; i++) {
	Denominator = CrvRmvd -> KnotVector[i + J + 2] - 
	              CrvRmvd -> KnotVector[i + J - k + 2];
	Enumerator = Knot2Remove - CrvRmvd -> KnotVector[i + J - k + 2];
	if ((IRIT_ABS(Denominator) < IRIT_EPS) || 
	    (IRIT_ABS(1 - Enumerator / Denominator)) < IRIT_EPS) {
	    CAGD_FATAL_ERROR(CAGD_ERR_NO_EXTRA_KNOTS);
	    IritFree(RmvCoef);
	    return NULL;
	}
	else 
	    /* i in the array is i + J - k + 2 in Boehm's Algo!! */
	    RmvCoef[i] = Enumerator / Denominator;  
    }

    /* Control Points - Some index conventions:                             */
    /* Upon Insertion (Boehm's Algo):                                       */ 
    /* Points removed are indices: J - k + 1, ..., J - 1; (total: k - 1)    */
    /* New points are indices: J - k + 1,..., J; (total: k)                 */
    /* Therefore, upon removal, we pretend we have a double J index,        */
    /* first of which is removed, second one is not.                        */
    MaxAxis = CAGD_NUM_OF_PT_COORD(Crv -> PType);
    for (d = !IsRational; d <= MaxAxis; d++) {
	CAGD_GEN_COPY(CrvRmvd -> Points[d], Crv -> Points[d], 
		      sizeof(CagdRType) * (J - k + 1));

	/* So far we just copied, now we actually solve the linear system: */
	CrvRmvd -> Points[d][J - 1] = (1 / (1 - RmvCoef[k - 2])) * 
	    (Crv ->Points[d][J] - RmvCoef[k - 2] * Crv -> Points[d][J + 1]);
	/* The first one was special. We moved a known point to the RHS of */
	/* the linear system. */
	for (i = J - 2; i > J - k; i--) {
	    CrvRmvd -> Points[d][i] =
	        (1 / (1 - RmvCoef[i + 1 - (J - k + 2)])) * 
	                                           (Crv -> Points[d][i + 1] - 
		RmvCoef[i + 1 - (J - k + 2)] * CrvRmvd -> Points[d][i + 1]);
	}
	/* And now copying again. */
	CAGD_GEN_COPY(CrvRmvd -> Points[d] + J, Crv -> Points[d] + J + 1,
		      sizeof(CagdRType) * (n + 1 - J));
    }

    /* copy the rest of the fields. */
    CrvRmvd -> Attr = AttrCopyAttributes(Crv -> Attr);
    CrvRmvd -> GType = Crv -> GType;
    CrvRmvd -> Periodic = Crv -> Periodic;
    CrvRmvd -> Pnext = Crv -> Pnext;

    IritFree(RmvCoef);

    return CrvRmvd;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Extension of a B-spline surface, in any (or more than one) of the four   M
* optional directions of the 2D domain. The domain is extended, such that    M
* the trace coincides with the original trace over the original domain.      M
*   Assumes open end conditions (in both knot vectors u and v).		     M
*									     *
* PARAMETERS:                                                                M
*   OrigSrf:	The surface to be extended.				     M
*   ExtDirs:	A vector of four boolean values to set the extension         M
*               directions. The convention is	MinU, MinV, MaxU, MaxV.      M
*		if NULL, all four directions are extended.		     M
*   EpsilonU:   The length of the extension in the u direction.		     M
*   EpsilonV:   The length of the extension in the v direction.		     M
*   RemoveExtraKnots: If FALSE, the resulting surface will not have minimal  M
*		multiplicity at the first internal knot on the extension     M
*		side. This is boolean controls all extensions that were	     M
*		performed, one decision for all of them.		     M
*									     *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:    The new extended surface.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvExtensionOneSide, BspCrvExtraKnotRmv			             M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspSrfExtension							     M
*****************************************************************************/
CagdSrfStruct *BspSrfExtension(const CagdSrfStruct *OrigSrf,
			       const CagdBType *ExtDirs,
			       CagdRType EpsilonU,
			       CagdRType EpsilonV,
			       CagdBType RemoveExtraKnots)
{
    CagdBType
	IsRational = CAGD_IS_RATIONAL_SRF(OrigSrf),
	AllDirs = ExtDirs == NULL;
    int i, j, d,
	MaxAxis = CAGD_NUM_OF_PT_COORD(OrigSrf -> PType);
    CagdRType *CrvP, *SrfP;
    CagdCrvStruct *TmpIsoParamCrv, *CTmp;
    CagdSrfStruct *Srf;

    if (!CAGD_IS_BSPLINE_SRF(OrigSrf)) {
        CAGD_FATAL_ERROR(CAGD_ERR_BZR_SRF_EXPECT);
	return NULL;
    }

    Srf = CagdSrfCopy(OrigSrf);

    /* Perform the required extension in the u direction. For a double index */
    /* (i,j),  for each fixed j we have an isoparameter u curve extension.   */
    if (AllDirs || ExtDirs[0] || ExtDirs[2]) {
        for (j = 0; j < Srf -> VLength; j++) {
	    TmpIsoParamCrv = BspSrfCrvFromMesh(Srf, j, CAGD_CONST_V_DIR);

	    if (AllDirs || ExtDirs[0]) {
	        CTmp = BspCrvExtensionOneSide(TmpIsoParamCrv, TRUE,
					      EpsilonU, RemoveExtraKnots);
		CagdCrvFree(TmpIsoParamCrv);
		TmpIsoParamCrv = CTmp;
	    }
	    if (AllDirs || ExtDirs[2]) {
	        CTmp = BspCrvExtensionOneSide(TmpIsoParamCrv, FALSE,
					      EpsilonU, RemoveExtraKnots);
		CagdCrvFree(TmpIsoParamCrv);
		TmpIsoParamCrv = CTmp;
	    }

	    /* The Isocurve is ready. Update corresponding pts in surface:   */
	    for (d = !IsRational; d <= MaxAxis; d++) {
	        CrvP = TmpIsoParamCrv -> Points[d];
		SrfP = Srf -> Points[d] + j * CAGD_NEXT_V(Srf);
		for (i = 0; i < TmpIsoParamCrv -> Length; i++) {
		    *SrfP = *CrvP++;
		    SrfP += CAGD_NEXT_U(Srf);
		}
	    }

	    if (j == Srf -> VLength - 1 && RemoveExtraKnots)
	        CAGD_GEN_COPY(Srf -> UKnotVector, TmpIsoParamCrv -> KnotVector,
			 sizeof(CagdRType) * (Srf -> ULength + Srf -> UOrder));

	    CagdCrvFree(TmpIsoParamCrv);
	}
    }

    /* Perform the required extension in the u direction. For a double index */
    /* (i,j),  for each fixed i we have an isoparameter v curve extension.   */
    if (AllDirs || ExtDirs[1] || ExtDirs[3]){
        for (i = 0; i < Srf -> ULength; i++) {
	    TmpIsoParamCrv = BspSrfCrvFromMesh(Srf, i, CAGD_CONST_U_DIR);

	    if (AllDirs || ExtDirs[1]) {
	        CTmp = BspCrvExtensionOneSide(TmpIsoParamCrv, TRUE,
					      EpsilonV, RemoveExtraKnots);
		CagdCrvFree(TmpIsoParamCrv);
		TmpIsoParamCrv = CTmp;
	    }
	    if (AllDirs || ExtDirs[3]) {
	        CTmp = BspCrvExtensionOneSide(TmpIsoParamCrv, FALSE,
					      EpsilonV, RemoveExtraKnots);
		CagdCrvFree(TmpIsoParamCrv);
		TmpIsoParamCrv = CTmp;
	    }

	    /* Isocurve is ready. Update corresponding pts in the surface:   */
	    for (d = !IsRational; d <= MaxAxis; d++) {
		CrvP = TmpIsoParamCrv -> Points[d];
		SrfP = Srf -> Points[d] + i * CAGD_NEXT_U(Srf);
		for (j = 0; j < TmpIsoParamCrv -> Length; j++) {
		    *SrfP = *CrvP++;
		    SrfP += CAGD_NEXT_V(Srf);
		}
	    }

	    if (i == Srf -> ULength - 1 && RemoveExtraKnots)
		CAGD_GEN_COPY(Srf -> VKnotVector, TmpIsoParamCrv -> KnotVector,
			 sizeof(CagdRType) * (Srf -> VLength + Srf -> VOrder));

	    CagdCrvFree(TmpIsoParamCrv);
	}
    }

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a list of curves representing the B-spline basis functions of    M
* the given space (order and knot sequence).                                 M
*   All basis functions will be scaled to fit into the unit square [0, 1]^2. M
*   If the space is periodic, Length should reflect this in the input (i.e.  M
* length should be enlarged by order - 1).				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Order:    Of space to create basis functions for.                        M
*   Length:   Number of control points in this space.                        M
*   KV:       Knot sequence of this space, of length (Order + Length).       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   Length curves representing the basis functions.       M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspGenKnotsGeometryAsCurves                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspGenBasisFuncsAsCurves                                                 M
*****************************************************************************/
CagdCrvStruct *BspGenBasisFuncsAsCurves(int Order,
					int Length,
					const CagdRType *KV)
{
    int i, j,
	Order1 = Order - 1;
    CagdCrvStruct *BssCrv,
        *BssCrvs = NULL;
    CagdRType *R, *BssKV,
        Gap = 0.0;

    if (Order < 1 || Length < Order) {
	CAGD_FATAL_ERROR(CAGD_ERR_WRONG_ORDER);
	return NULL;
    }
    if (KV == NULL) {
        CAGD_FATAL_ERROR(CAGD_ERR_NO_KV_FOUND);
	return NULL;
    }

    BssCrv = BspCrvNew(Length + 2 * (Order - 1), Order, CAGD_PT_E2_TYPE);
    BssKV = BssCrv -> KnotVector;

    /* analyze the average gap between knots. */
    for (i = j = 0; i < Order + Length; i++) {
        if (KV[i + 1] - KV[i] > 0) {
	    j++;
	    Gap += KV[i + 1] - KV[i];
	}
    }
    Gap /= j;

    /* Copy the KV, add Order-1 new knots at each end, and map to [0, 1]. */
    CAGD_GEN_COPY(&BssKV[Order1], KV, sizeof(CagdRType) * (Order + Length));
    for (i = Order1; i > 0; i--)
        BssKV[i - 1] = BssKV[i] - Gap;
    for (i = BssCrv -> Length + 1; i < BssCrv -> Length + Order; i++)
        BssKV[i] = BssKV[i - 1] + Gap;

    BspKnotScale(BssKV, BssCrv -> Length + Order,
	         1.0 / (BssKV[BssCrv -> Length - Order1] - BssKV[Order1 * 2]));
    BspKnotTranslate(BssKV, BssCrv -> Length + Order, -BssKV[Order1 * 2]);

    /* Populate X axis with node values and zap Y axis of the curve. */
    R = BspKnotNodes(BssKV, BssCrv -> Length + Order, Order);
    CAGD_GEN_COPY(BssCrv -> Points[1], R, sizeof(CagdRType) * BssCrv -> Length);
    IritFree(R);
    IRIT_ZAP_MEM(BssCrv -> Points[2], sizeof(CagdRType) * BssCrv -> Length);

    for (j = 0; j < Length; j++) {
	CagdCrvStruct
	    *TCrv = CagdCrvCopy(BssCrv);

	TCrv -> Points[2][j + Order1] = 1.0;

	IRIT_LIST_PUSH(TCrv, BssCrvs);
    }

    CagdCrvFree(BssCrv);

    return BssCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a list of linear B-spline curves representing the B-spline knots M
* in the given space (order and knot sequence).                              M
*   All knots will be scaled to fit just below the unit square [0, 1]^2.     M
*   If the space is periodic, Length should reflect this in the input (i.e.  M
* length should be enlarged by order - 1).				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Order:    Of space to create the geometry of the knots for.              M
*   Length:   Number of control points in this space.                        M
*   KV:       Knot sequence of this space, of length (Order + Length).       M
*   SizeOfKnot:  The size of the plot knot.  Knots are created as triangles  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   Length curves representing the Length knots.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspGenBasisFuncsAsCurves                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   BspGenKnotsGeometryAsCurves                                              M
*****************************************************************************/
CagdCrvStruct *BspGenKnotsGeometryAsCurves(int Order,
					   int Length,
					   const CagdRType *KV,
					   CagdRType SizeOfKnot)
{
    /* Draw the knots. */
    int i,
	KnotNum = Length + Order;
    CagdCrvStruct
        *KnotCrvs = NULL,
        *KnotCrv = BspCrvNew(4, 2, CAGD_PT_E2_TYPE);
    CagdRType TMin, TMax, Dt,
        **Pts = KnotCrv -> Points,
        *MultiplicityKV = (CagdRType *) IritMalloc(sizeof(CagdRType) *
						                    KnotNum);

    BspKnotUniformOpen(4, 2, KnotCrv -> KnotVector);

    /* MultiplicityVec[k] = how many knots above knot k in display plot. */
    MultiplicityKV[0] = 0;
    for (i = 1; i < KnotNum; i++) {
	if (KV[i] == KV[i - 1])
	    MultiplicityKV[i] = MultiplicityKV[i - 1] + 1;
	else
	    MultiplicityKV[i] = 0;
    }

    /* Create a canonical knot as a triangle curve. */
    Pts[1][0] = 0.0;
    Pts[2][0] = 0.0;
    Pts[1][1] = SizeOfKnot * 0.25;
    Pts[2][1] = -SizeOfKnot;
    Pts[1][2] = SizeOfKnot * -0.25;
    Pts[2][2] = -SizeOfKnot;
    Pts[1][3] = 0.0;
    Pts[2][3] = 0.0;

    TMin = KV[Order - 1];
    TMax = KV[Length];
    Dt = TMax - TMin;

    for (i = 0; i < KnotNum; i++) {
	CagdPType Trans;
	CagdCrvStruct
	    *TCrv = CagdCrvCopy(KnotCrv);

	Trans[0] = (KV[i] - TMin) / Dt;
	Trans[1] = -MultiplicityKV[i] * SizeOfKnot;
	Trans[2] = 0.0;
	CagdCrvTransform(TCrv, Trans, 1.0);

	IRIT_LIST_PUSH(TCrv, KnotCrvs);
    }

    CagdCrvFree(KnotCrv);

    return KnotCrvs;
}
