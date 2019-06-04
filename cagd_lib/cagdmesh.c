/******************************************************************************
* CagdMesh.c - Extract surface control mesh/curve control polygon as polyline *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Aug. 90.					      *
******************************************************************************/

#include "cagd_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
* Extracts the control polygon of a curve as a polyline.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:        To extract a control polygon from.                           M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPolylineStruct *:  The control polygon of Crv.                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdCrv2CtrlPoly, control polygon                                        M
*****************************************************************************/
CagdPolylineStruct *CagdCrv2CtrlPoly(const CagdCrvStruct *Crv)
{
    int i,
	Length = Crv -> Length + (Crv -> Periodic != FALSE);
    CagdRType
	* const *CrvP = Crv -> Points;
    CagdPolylnStruct *NewPolyline;
    CagdPolylineStruct *P;

    P = CagdPolylineNew(Length);
    NewPolyline = P -> Polyline;

    for (i = 0; i < Length; i++) {
	CagdCoerceToE3(NewPolyline -> Pt, CrvP, i % Crv -> Length,
		       Crv -> PType);
	NewPolyline++;
    }
    return P;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Extracts the control mesh of a surface as a list of polylines.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To extract a control mesh from.		                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPolylineStruct *:  The control mesh of Srf.                          M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdSrf2KnotLines, CagdSrf2KnotCurves                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdSrf2CtrlMesh, control mesh                                           M
*****************************************************************************/
CagdPolylineStruct *CagdSrf2CtrlMesh(const CagdSrfStruct *Srf)
{
    int i, j,
	ULength = Srf -> ULength + (Srf -> UPeriodic != FALSE),
	VLength = Srf -> VLength + (Srf -> VPeriodic != FALSE);
    CagdRType
	* const *SrfP = Srf -> Points;
    CagdPolylnStruct *NewPolyline;
    CagdPolylineStruct *P,
	*PList = NULL;

    for (j = 0; j < VLength; j++) {	   /* Generate the rows of the mesh. */
	P = CagdPolylineNew(ULength);
	NewPolyline = P -> Polyline;

	for (i = 0; i < ULength; i++) {
	    CagdCoerceToE3(NewPolyline -> Pt, SrfP,
			   CAGD_MESH_UV(Srf, i % Srf -> ULength,
					     j % Srf -> VLength),
			   Srf -> PType);
            NewPolyline++;
	}
	IRIT_LIST_PUSH(P, PList);
    }

    for (i = 0; i < ULength; i++) {	   /* Generate the cols of the mesh. */
	P = CagdPolylineNew(VLength);
	NewPolyline = P -> Polyline;

	for (j = 0; j < VLength; j++) {
	    CagdCoerceToE3(NewPolyline -> Pt, SrfP,
			   CAGD_MESH_UV(Srf, i % Srf -> ULength,
					     j % Srf -> VLength),
			   Srf -> PType);
            NewPolyline++;
	}
	IRIT_LIST_PUSH(P, PList);
    }

    return PList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Extracts a polyline grid along all knot in U and V as a list of polylines. M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To extract a control mesh from.		                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPolylineStruct *:  The control mesh of Srf.                          M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdSrf2CtrlMesh, CagdSrf2KnotCurves                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdSrf2KnotLines, control mesh, knots, knot lines                       M
*****************************************************************************/
CagdPolylineStruct *CagdSrf2KnotLines(const CagdSrfStruct *Srf)
{
    int i, j, ULen, VLen, *Mult;
    CagdRType *R, *UVals, *VVals;
    CagdPolylnStruct *NewPolyline;
    CagdPolylineStruct *P,
	*PList = NULL;

    if (CAGD_IS_BEZIER_SRF(Srf)) {
        UVals = IritMalloc(sizeof(CagdRType) * 2);
        VVals = IritMalloc(sizeof(CagdRType) * 2);
	CagdSrfDomain(Srf, &UVals[0], &UVals[1], &VVals[0], &VVals[1]);
        ULen = 2;
        VLen = 2;
    }
    else if (CAGD_IS_BSPLINE_SRF(Srf)) {
	Mult = IritMalloc(sizeof(int) *
			  IRIT_MAX(Srf -> UOrder + Srf -> ULength,
				   Srf -> VOrder + Srf -> VLength));
        UVals = IritMalloc(sizeof(CagdRType) *
			   (Srf -> UOrder + Srf -> ULength));
	ULen = BspKnotsMultiplicityVector(Srf -> UKnotVector,
					  Srf -> UOrder + Srf -> ULength,
					  UVals, Mult);
        VVals = IritMalloc(sizeof(CagdRType) *
			   (Srf -> VOrder + Srf -> VLength));
	VLen = BspKnotsMultiplicityVector(Srf -> VKnotVector,
					  Srf -> VOrder + Srf -> VLength,
					  VVals, Mult);
	IritFree(Mult);
    }
    else {
	CAGD_FATAL_ERROR(CAGD_ERR_UNDEF_SRF);
	return NULL;
    }

    for (j = 0; j < VLen; j++) {	   /* Generate the rows of the mesh. */
	P = CagdPolylineNew(ULen);
	NewPolyline = P -> Polyline;

	for (i = 0; i < ULen; i++) {
	    R = CagdSrfEval(Srf, UVals[i], VVals[j]);
	    CagdCoerceToE3(NewPolyline -> Pt, &R, - 1, Srf -> PType);
            NewPolyline++;
	}
	IRIT_LIST_PUSH(P, PList);
    }

    for (i = 0; i < ULen; i++) {	   /* Generate the cols of the mesh. */
	P = CagdPolylineNew(VLen);
	NewPolyline = P -> Polyline;

	for (j = 0; j < VLen; j++) {
	    R = CagdSrfEval(Srf, UVals[i], VVals[j]);
	    CagdCoerceToE3(NewPolyline -> Pt, &R, - 1, Srf -> PType);
            NewPolyline++;
	}
	IRIT_LIST_PUSH(P, PList);
    }

    IritFree(UVals);
    IritFree(VVals);

    return PList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Extracts a list of isoparametric curves in a grid along all knots in U and M
* V.									     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To extract a control mesh from.		                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   Isoparameteric curves along all knot lines.           M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdSrf2CtrlMesh, CagdSrf2KnotLines                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdSrf2KnotCurves, control mesh, knots, knot lines, knot curves         M
*****************************************************************************/
CagdCrvStruct *CagdSrf2KnotCurves(const CagdSrfStruct *Srf)
{
    int i, j, ULen, VLen, *Mult;
    CagdRType *UVals, *VVals;
    CagdCrvStruct *Crv,
        *KnotCurves = NULL;

    if (CAGD_IS_BEZIER_SRF(Srf)) {
        UVals = IritMalloc(sizeof(CagdRType) * 2);
        VVals = IritMalloc(sizeof(CagdRType) * 2);
	CagdSrfDomain(Srf, &UVals[0], &UVals[1], &VVals[0], &VVals[1]);
        ULen = 2;
        VLen = 2;
    }
    else if (CAGD_IS_BSPLINE_SRF(Srf)) {
	Mult = IritMalloc(sizeof(int) *
			  IRIT_MAX(Srf -> UOrder + Srf -> ULength,
				   Srf -> VOrder + Srf -> VLength));
        UVals = IritMalloc(sizeof(CagdRType) *
			   (Srf -> UOrder + Srf -> ULength));
	ULen = BspKnotsMultiplicityVector(Srf -> UKnotVector,
					  Srf -> UOrder + Srf -> ULength,
					  UVals, Mult);
        VVals = IritMalloc(sizeof(CagdRType) *
			   (Srf -> VOrder + Srf -> VLength));
	VLen = BspKnotsMultiplicityVector(Srf -> VKnotVector,
					  Srf -> VOrder + Srf -> VLength,
					  VVals, Mult);
	IritFree(Mult);
    }
    else {
	CAGD_FATAL_ERROR(CAGD_ERR_UNDEF_SRF);
	return NULL;
    }

    for (j = 0; j < VLen; j++) {	   /* Generate the rows of the mesh. */
        Crv = CagdCrvFromSrf(Srf, VVals[j], CAGD_CONST_V_DIR);
	IRIT_LIST_PUSH(Crv, KnotCurves);
    }

    for (i = 0; i < ULen; i++) {	   /* Generate the cols of the mesh. */
        Crv = CagdCrvFromSrf(Srf, UVals[i], CAGD_CONST_U_DIR);
	IRIT_LIST_PUSH(Crv, KnotCurves);
    }

    IritFree(UVals);
    IritFree(VVals);

    return KnotCurves;
}
