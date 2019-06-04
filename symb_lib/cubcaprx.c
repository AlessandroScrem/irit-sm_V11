/******************************************************************************
* cubcapprox.c - piecewise cubic approximation of freeform curves.  	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, July 13	  				      *
******************************************************************************/

#include "symb_loc.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes a piecewise cubic approximation to given freeform planar curve. M
*   The following steps are performed during this approximation process:     M
*    a. Fit the curve with a C^1 continuous cubic that is tangent to the     M
*       curves's end points.						     M
*    b. If fit is good enough, we stop.  Otherwise subdivide region into two M
*	and recursively invoke step 2 on both halfs.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:        2D Curve to approximate using piecewise cubics.              M
*   Tolerance:  Of approximation, in Hausdorff distance sense.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:  List of cubics approximating Crv to within Tolerance.  M
*                     List will be C^1.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbCrvBiArcApprox, SymbApproxCrvAsBzrCubics                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   SymbCrvCubicApprox                                                       M
*****************************************************************************/
CagdCrvStruct *SymbCrvCubicApprox(const CagdCrvStruct *Crv,
				  CagdRType Tolerance)
{
    CagdRType TMin, TMax, *R;
    CagdPType PtStart, PtEnd;
    CagdVecStruct TanStart, TanEnd;
    CagdBBoxStruct BBox;
    CagdCrvStruct *DistCrv, *DistCrvSqr, *CubicCrv, *SubdivCrvs;

    /* Computes curve's end positions and tangents. */
    CagdCrvDomain(Crv, &TMin, &TMax);
    R = CagdCrvEval(Crv, TMin);
    CagdCoerceToE3(PtStart, &R, -1, Crv -> PType);
    R = CagdCrvEval(Crv, TMax);
    CagdCoerceToE3(PtEnd, &R, -1, Crv -> PType);
    TanStart = *CagdCrvTangent(Crv, TMin, FALSE);
    IRIT_VEC_SCALE(TanStart.Vec, (TMax - TMin) / 1.0);
    TanEnd = *CagdCrvTangent(Crv, TMax, FALSE);
    IRIT_VEC_SCALE(TanEnd.Vec, (TMax - TMin) / 1.0);
      
    /* Create a cubic Hermite out of this data: */
    CubicCrv = CagdCubicHermiteCrv(PtStart, PtEnd, TanStart.Vec, TanEnd.Vec);
    CagdCrvSetDomain(CubicCrv, TMin, TMax);

    /* Compute an Hausdorff distance bound: */
    DistCrv = SymbCrvSub(CubicCrv, Crv);
    DistCrvSqr = SymbCrvDotProd(DistCrv, DistCrv);
    CagdCrvFree(DistCrv);
    CagdCrvBBox(DistCrvSqr, &BBox);
    CagdCrvFree(DistCrvSqr);

    if (IRIT_SQR(Tolerance) > BBox.Max[0]) {
	/* We are done - cubic fit is close enough.  Stop. */
	return CubicCrv;
    }
    else {                                          /* Divide and conquer. */
	CagdCrvFree(CubicCrv);
	SubdivCrvs = CagdCrvSubdivAtParam(Crv, (TMin + TMax) * 0.5);
	CubicCrv = CagdListAppend(SymbCrvCubicApprox(SubdivCrvs, Tolerance),
				  SymbCrvCubicApprox(SubdivCrvs -> Pnext,
						     Tolerance));
	CagdCrvFreeList(SubdivCrvs);
	return CubicCrv;
    }
}
