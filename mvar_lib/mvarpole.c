/******************************************************************************
* mvarpole.c - computes poles of rational forms.		              *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Dec. 2012.					      *
******************************************************************************/

#include <assert.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/geom_lib.h"
#include "mvar_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the poles of a rational surface, solving for the zeros of the   M
* surface's denominator.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:          Rational surface to extract its poles.                     M
*   SubdivTol:	  The subdivision tolerance to use.			     M
*   NumericTol:	  The numerical tolerance to use.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *:   The poles, as piecewise linear approximations.   M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdPointsHasPoles, SymbCrvsSplitPoleParams                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarRationalSrfsPoles                                                    M
*****************************************************************************/
MvarPolylineStruct *MvarRationalSrfsPoles(const CagdSrfStruct *Srf,
				          CagdRType SubdivTol,
				          CagdRType NumericTol)
{
    CagdSrfStruct *SrfW, *SrfX, *SrfY, *SrfZ;
    MvarMVStruct *MVs;
    MvarConstraintType Constr;
    MvarPolylineStruct *Poles;

    if (!CAGD_IS_RATIONAL_SRF(Srf)) {
        return NULL;
    }

    SymbSrfSplitScalar(Srf, &SrfW, &SrfX, &SrfY, &SrfZ);
    if (SrfX != NULL)
        CagdSrfFree(SrfX);
    if (SrfY != NULL)
        CagdSrfFree(SrfY);
    if (SrfZ!= NULL)
        CagdSrfFree(SrfZ);

    assert(SrfW != NULL);
    MVs = MvarSrfToMV(SrfW);
    CagdSrfFree(SrfW);

    Constr = MVAR_CNSTRNT_ZERO;

    Poles = MvarMVsZeros1D(&MVs, NULL, 1, SubdivTol, SubdivTol, NumericTol);

    MvarMVFree(MVs);

    return Poles;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Split the given rational surface at its poles, if any, by solving for    M
* the zeros of the surface's denominator.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:          Rational surface to split at its poles.                    M
*   SubdivTol:	  The subdivision tolerance to use.			     M
*   NumericTol:	  The numerical tolerance to use.			     M
*   OutReach:     Small offset to clip poles regions at, zero to disable.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrimSrfStruct *:   Trimmed surfaces, divided at all poles, or NULL if    M
*                      no poles.				             M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdPointsHasPoles, SymbCrvsSplitPoleParams                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarSrfSplitPoleParams                                                   M
*****************************************************************************/
TrimSrfStruct *MvarSrfSplitPoleParams(const CagdSrfStruct *Srf,
				      CagdRType SubdivTol,
				      CagdRType NumericTol,
				      CagdRType OutReach)
{
    MvarPolylineStruct
        *Poles = MvarRationalSrfsPoles(Srf, SubdivTol, NumericTol);
    IPObjectStruct *PObj;
    IPPolygonStruct *PlPoles;
    TrimSrfStruct *TSrfs;

    if (Poles == NULL)
        return NULL;

    PObj = MvarCnvrtMVPolysToIritPolys2(Poles, TRUE);
    MvarPolylineFreeList(Poles);
    PlPoles = PObj -> U.Pl;
    PObj -> U.Pl = NULL;
    IPFreeObject(PObj);

    if (OutReach > 0) {  /* Compute offsets instead of exact poles curves. */
        CagdRType UMin, UMax, VMin, VMax;
        IPPolygonStruct *Pl, *NewPl,
	    *NewPoles = NULL;

	CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

	for (Pl = PlPoles; Pl != NULL; Pl = Pl -> Pnext) {
	    int i, k, ClosedLoop;
	    IPVertexStruct
	        *VFirst = Pl -> PVertex,
	        *VLast = IPGetLastVrtx(VFirst);

	    if (VLast == VFirst)   /* Skip polylines with a single vertex. */
	        continue;
	    ClosedLoop = IRIT_PT_APX_EQ_E2(VFirst -> Coord, VLast -> Coord);

	    for (i = 0; i < 2; i++) {
		NewPl = GMPolyOffset(Pl, ClosedLoop,
				     i == 0 ? OutReach : -OutReach, NULL);

		VFirst = NewPl -> PVertex;
		VLast = IPGetLastVrtx(VFirst);

		/* Handle the end conditions. */
		if (ClosedLoop) {
		    for (k = 0 ; k < 2; k++) {
		        VFirst -> Coord[k] =
			    VLast -> Coord[k] = (VFirst -> Coord[k] +
						 VLast -> Coord[k]) * 0.5;
		    }
		}

		IRIT_LIST_PUSH(NewPl, NewPoles);
	    }
	}
	
	IPFreePolygonList(PlPoles);
	PlPoles = NewPoles;
    }

    TSrfs = TrimSrfsFromContours(Srf, PlPoles);
    IPFreePolygonList(PlPoles);

    return TSrfs;
}

