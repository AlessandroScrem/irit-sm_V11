/******************************************************************************
* UsrCnvrt.c - general conversions.					      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, May 96.					      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/cagd_lib.h"
#include "user_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
* Returns a list of linear Bspline curves constructed from given polylines.  M
*                                                                            *
* PARAMETERS:                                                                M
*   Polys:      To convert to linear bspline curves.                         M
*   FilterDups: If TRUE, filters out duplicates points in polygon, in        M
*		place.							     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:  Linear Bspline curves representing Poly.               M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCnvrtPolyline2LinBspCrv, UserPolyline2LinBsplineCrv                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserPolylines2LinBsplineCrvs, linear curves, conversion                  M
*****************************************************************************/
CagdCrvStruct *UserPolylines2LinBsplineCrvs(const IPPolygonStruct *Polys,
					    CagdBType FilterDups)
{
    CagdCrvStruct
	*Crvs = NULL;
    const IPPolygonStruct *Poly;

    for (Poly = Polys; Poly != NULL; Poly = Poly -> Pnext) {
	CagdCrvStruct
	    *Crv = UserPolyline2LinBsplineCrv(Poly, FilterDups);

	if (Crv != NULL)
	    IRIT_LIST_PUSH(Crv, Crvs);
    }

    return Crvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Returns a linear Bspline curve constructed from given polyline.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Poly:       To convert to a linear bspline curve.                        M
*   FilterDups: If TRUE, filters out duplicates points, in polygon, in       M
*		place.							     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *: A linear Bspline curve representing Poly.               M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCnvrtPolyline2LinBspCrv, UserPolylines2LinBsplineCrvs                M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserPolyline2LinBsplineCrv, linear curves, conversion                    M
*****************************************************************************/
CagdCrvStruct *UserPolyline2LinBsplineCrv(const IPPolygonStruct *Poly,
					  CagdBType FilterDups)
{
    IPVertexStruct *V;
    int i, Length;
    CagdCrvStruct *Crv;
    CagdRType **Points;

    if (FilterDups && Poly -> PVertex != NULL) {
	for (V = Poly -> PVertex; V -> Pnext != NULL; ) {
	    if (IRIT_PT_APX_EQ(V -> Coord, V -> Pnext -> Coord)) {
		IPVertexStruct
		    *VNext = V -> Pnext -> Pnext;

		IPFreeVertex(V -> Pnext);
		V -> Pnext = VNext;
	    }
	    else
		V = V -> Pnext;

	    if (V == Poly -> PVertex)
		break;
	}
    }

    V = Poly -> PVertex;
    if ((Length = IPVrtxListLen(V)) >= 2) {
	Crv = BspCrvNew(Length, 2, CAGD_PT_E3_TYPE);
	Points = Crv -> Points;

	BspKnotUniformOpen(Length, 2, Crv -> KnotVector);
	BspKnotAffineTrans2(Crv -> KnotVector,
			    Crv -> Length + Crv -> Order, 0, 1);

	for (i = 0; i < Length; i++, V = V -> Pnext) {
	    Points[1][i] = V -> Coord[0];
	    Points[2][i] = V -> Coord[1];
	    Points[3][i] = V -> Coord[2];
	}

	return Crv;
    }
    else
	return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a polyline into an Irit polyline.                               M
*                                                                            *
* PARAMETERS:                                                                M
*   Poly:   Input polyline to convert into am Irit polyline.	             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   Converted polyline.                                 M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCnvrtLinBspCrv2Polyline, CagdCnvrtPolyline2LinBspCrv,		     M
*   CagdCnvrtPolyline2PtList, UserCagdPolylines2IritPolylines                M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserCagdPolyline2IritPolyline                                            M
*****************************************************************************/
IPPolygonStruct *UserCagdPolyline2IritPolyline(const CagdPolylineStruct *Poly)
{
    int i,
	PlLen = Poly -> Length;
    IPVertexStruct
        *VHead = NULL;

    for (i = 0; i < PlLen; i++) {
        VHead = IPAllocVertex2(VHead);
        IRIT_PT_COPY(VHead -> Coord, Poly -> Polyline[i].Pt);
    }

    if (VHead == NULL)
        return NULL;

    return IPAllocPolygon(0, IPReverseVrtxList2(VHead), NULL);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of polylines into an Irit polylines.                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Polys:   Input polylines to convert into Irit polylines.	             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   Converted polylines.                                M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCnvrtLinBspCrv2Polyline, CagdCnvrtPolyline2LinBspCrv,		     M
*   CagdCnvrtPolyline2PtList, UserCagdPolyline2IritPolyline                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserCagdPolylines2IritPolylines                                          M
*****************************************************************************/
IPPolygonStruct *UserCagdPolylines2IritPolylines(const CagdPolylineStruct
						                      *Polys)
{
    const CagdPolylineStruct *Poly;
    IPPolygonStruct *Pl,
        *RetPls = NULL;

    for (Poly = Polys; Poly != NULL; Poly = Poly -> Pnext) {
        Pl = UserCagdPolyline2IritPolyline(Poly);
	IRIT_LIST_PUSH(Pl, RetPls);
    }

    return IPReversePlList(RetPls);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a cagd linear curve into an Irit polyline.		             M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:   Input linear curve to convert into an Irit polyline.	             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   Converted polyline.                                 M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCnvrtLinBspCrv2Polyline, CagdCnvrtPolyline2LinBspCrv,		     M
*   CagdCnvrtPolyline2PtList, UserCagdPolyline2IritPolyline,                 M
*   UserCnvrtLinBspCrvs2IritPolylines					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserCnvrtLinBspCrv2IritPolyline                                          M
*****************************************************************************/
IPPolygonStruct *UserCnvrtLinBspCrv2IritPolyline(const CagdCrvStruct *Crv)
{
    CagdPolylineStruct
        *CagdPl = CagdCnvrtLinBspCrv2Polyline(Crv);
    IPPolygonStruct
        *Pl = UserCagdPolyline2IritPolyline(CagdPl);

    CagdPolylineFree(CagdPl);

    return Pl;    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list of cagd linear curves into Irit polylines.               M
*                                                                            *
* PARAMETERS:                                                                M
*   Crvs:   Input linear curves to convert into Irit polylines.	             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:   Converted polylines.                                M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdCnvrtLinBspCrv2Polyline, CagdCnvrtPolyline2LinBspCrv,		     M
*   CagdCnvrtPolyline2PtList, UserCagdPolyline2IritPolyline,                 M
*   UserCnvrtLinBspCrv2IritPolyline					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserCnvrtLinBspCrvs2IritPolylines                                        M
*****************************************************************************/
IPPolygonStruct *UserCnvrtLinBspCrvs2IritPolylines(const CagdCrvStruct *Crvs)
{
    const CagdCrvStruct *Crv;
    IPPolygonStruct *Pl,
        *RetPls = NULL;

    for (Crv = Crvs; Crv != NULL; Crv = Crv -> Pnext) {
        Pl = UserCnvrtLinBspCrv2IritPolyline(Crv);
	IRIT_LIST_PUSH(Pl, RetPls);
    }

    return IPReversePlList(RetPls);
}
