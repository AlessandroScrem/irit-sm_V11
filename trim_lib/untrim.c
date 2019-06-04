/******************************************************************************
* UnTrim.c - computes trimmed surfaces decomposition into tensor products.    *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Nov. 2013.					      *
******************************************************************************/

#include <assert.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/symb_lib.h"
#include "trim_loc.h"

static CagdSrfStruct *TrimOrientUVReparamSrf(const CagdSrfStruct *UVSrf);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a trimmed surface - subdivides it into trimmed Bezier surfaces (each M
* spanning domain [0, 1]^2).			                             M
*    Returns a list of trimmed Bezier surfaces, that together, are identical M
* geometrically to the input trimmed surface.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   TrimSrf:  A Bezier or a Bspline trimmed surface to convert to Bezier.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrimSrfStruct *:  The subdivided Bezier trimmed surfaces.                M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrimSrfSubdivAtParam                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrimSrfCnvrt2BzrTrimSrf, subdivision                                     M
*****************************************************************************/
TrimSrfStruct *TrimSrfCnvrt2BzrTrimSrf(TrimSrfStruct *TrimSrf)
{
    CagdBType OldPLCrvsVal;
    CagdSrfStruct
	*Srf = TrimSrf -> Srf;
    TrimSrfStruct *TrimSrf1, *TrimSrf2, *TrimSrf1Bzrs, *TrimSrf2Bzrs;

    if (CAGD_IS_BEZIER_SRF(Srf))
        return TrimSrfCopy(TrimSrf);
    else if (!CAGD_IS_BSPLINE_SRF(Srf)) {
	CAGD_FATAL_ERROR(CAGD_ERR_WRONG_SRF);
	return NULL;
    }

    if (Srf -> ULength > Srf -> UOrder) {
	CagdRType
	    t = Srf -> UKnotVector[(Srf -> ULength + Srf -> UOrder) >> 1];

	OldPLCrvsVal = TrimSrfSubdivAtParamForcePiecwiseLinearTrimCrvs(FALSE);
	TrimSrf1 = TrimSrfSubdivAtParam(TrimSrf, t, CAGD_CONST_U_DIR);
	TrimSrfSubdivAtParamForcePiecwiseLinearTrimCrvs(OldPLCrvsVal);

	TrimSrf2 = TrimSrf1 -> Pnext;
	TrimSrf1 -> Pnext = NULL;

	TrimSrf1Bzrs = TrimSrfCnvrt2BzrTrimSrf(TrimSrf1);
	TrimSrfFree(TrimSrf1);

	if (TrimSrf2 != NULL) {
	    TrimSrf2Bzrs = TrimSrfCnvrt2BzrTrimSrf(TrimSrf2);
	    TrimSrfFree(TrimSrf2);

	    return CagdListAppend(TrimSrf1Bzrs, TrimSrf2Bzrs);
	}
	else
	    return TrimSrf1Bzrs;
    }
    else if (Srf -> VLength > Srf -> VOrder) {
	CagdRType
	    t = Srf -> VKnotVector[(Srf -> VLength + Srf -> VOrder) >> 1];

	OldPLCrvsVal = TrimSrfSubdivAtParamForcePiecwiseLinearTrimCrvs(FALSE);
	TrimSrf1 = TrimSrfSubdivAtParam(TrimSrf, t, CAGD_CONST_V_DIR);
	TrimSrfSubdivAtParamForcePiecwiseLinearTrimCrvs(OldPLCrvsVal);

	TrimSrf2 = TrimSrf1 -> Pnext;
	TrimSrf1 -> Pnext = NULL;

	TrimSrf1Bzrs = TrimSrfCnvrt2BzrTrimSrf(TrimSrf1);
	TrimSrfFree(TrimSrf1);

	if (TrimSrf2 != NULL) {
	    TrimSrf2Bzrs = TrimSrfCnvrt2BzrTrimSrf(TrimSrf2);
	    TrimSrfFree(TrimSrf2);

	    return CagdListAppend(TrimSrf1Bzrs, TrimSrf2Bzrs);
	}
	else
	    return TrimSrf1Bzrs;
    }
    else {
        CagdRType UMin, UMax, VMin, VMax;

	CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);

	/* Map the domain of the trimming curves to be [0, 1]^2. */
	TrimSrf1 = TrimAffineTransTrimSrf(TrimSrf, 0.0, 1.0, 0.0, 1.0);
	Srf = TrimSrf1 -> Srf;

	/* Update the new surface to be a Bezier surface. */
	Srf -> GType = CAGD_SBEZIER_TYPE;
	IritFree(Srf -> UKnotVector);
	IritFree(Srf -> VKnotVector);
	Srf -> UKnotVector = NULL;
	Srf -> VKnotVector = NULL;

	/* Keep old Bspline domain, just in case. */
	AttrSetUVAttrib(&TrimSrf1 -> Attr, "BspDomainMin", UMin, VMin);
	AttrSetUVAttrib(&TrimSrf1 -> Attr, "BspDomainMax", UMax, VMax);

	return TrimSrf1;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Make sure the UV surface is properly oriented and correct.	             *
*                                                                            *
* PARAMETERS:                                                                *
*   UVSrf:   To examine its orientation and correct if needed.		     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdSrfStruct *:   Corrected, if necessary, UVSrf.                       *
*****************************************************************************/
static CagdSrfStruct *TrimOrientUVReparamSrf(const CagdSrfStruct *UVSrf)
{
    CagdRType UMin, UMax, VMin, VMax;
    CagdVecStruct *Nrml;

    CagdSrfDomain(UVSrf, &UMin, &UMax, &VMin, &VMax);
    Nrml = CagdSrfNormal(UVSrf, (UMin + UMax) * 0.5,
			        (VMin + VMax) * 0.5, FALSE);
    if (Nrml -> Vec[2] < 0.0)
        return CagdSrfCopy(UVSrf);
    else
	return CagdSrfReverse(UVSrf);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a trimmed surface - subdivides it into trimmed Bezier surfaces and   M
* convert the trimmed surface into regular Bezier/Bspline patches.           M
*    Returns a list of non-trimmed Bezier surfaces, that together, are       M
* identical geometrically to the input trimmed surface.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   TrimSrf:  A trimmed Bezier or a Bspline surface to convert to            M
*             non-trimmed Bezier surfaces.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:  The non trimmed regular Bezier surfaces.               M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrimSrfSubdivAtParam, TrimSrfCnvrt2BzrTrimSrf                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrimSrfCnvrt2BzrRglrSrf, subdivision                                     M
*****************************************************************************/
CagdSrfStruct *TrimSrfCnvrt2BzrRglrSrf(TrimSrfStruct *TrimSrf)
{
    TrimSrfStruct *TBSrf,
        *TrimBzrSrfs = TrimSrfCnvrt2BzrTrimSrf(TrimSrf);
    CagdSrfStruct
        *BezSrfs = NULL;
    int MinimizeSize = AttrGetIntAttrib(TrimSrf -> Attr, "MinSize");

    if (IP_ATTR_IS_BAD_INT(MinimizeSize))
        MinimizeSize = FALSE;

    for (TBSrf = TrimBzrSrfs; TBSrf != NULL; TBSrf = TBSrf -> Pnext) {
	TrimCrvSegStruct *TCSeg;
	CagdCrvStruct *UVCrv,
	     *UVCrvs = NULL;
	CagdSrfStruct *UVSrfs, *CompSrf;
	TrimCrvStruct *TCrvs;

#	ifdef DEBUG
	{
	    TCrvs = TrimChainTrimmingCurves2Loops(TBSrf -> TrimCrvList);

	    if (TCrvs -> Pnext != NULL) {
	        IRIT_INFO_MSG_PRINTF("Conversion of trimmed surfaces to Bezier assumes trimmed surfaces with single\ntrimming loop only.  Only first loop is considered.\n");
	    }
	    TrimCrvFreeList(TCrvs);
	}
#	endif /* DEBUG */

	/* Chain all the trimming curves into one linear list. */
	for (TCrvs = TBSrf -> TrimCrvList;
	     TCrvs != NULL;
	     TCrvs = TCrvs -> Pnext) {
	    for (TCSeg = TCrvs -> TrimCrvSegList;
		 TCSeg != NULL;
		 TCSeg = TCSeg -> Pnext) {
	        UVCrv = CagdCrvCopy(TCSeg -> UVCrv);
		IRIT_LIST_PUSH(UVCrv, UVCrvs);
	    }
	}

	/* Build bivariate surfaces from this closed loop.  Expects at most */
	/* two UV surfaces to be returned.				    */
	UVSrfs = CagdSrfFromNBndryCrvs(UVCrvs, MinimizeSize);
	if (UVSrfs != NULL) {
	    CagdSrfStruct *UVSrf;

	    /* Make sure the UVSrfs are properly oriented before composing. */
	    for (UVSrf = UVSrfs; UVSrf != NULL; UVSrf = UVSrf -> Pnext) {
	        CagdSrfStruct
		    *OrientedUVSrf = TrimOrientUVReparamSrf(UVSrf);

		CompSrf = SymbComposeSrfSrf(TBSrf -> Srf, OrientedUVSrf);
		CagdSrfFree(OrientedUVSrf);

		IRIT_LIST_PUSH(CompSrf, BezSrfs);
	    }

	    CagdSrfFreeList(UVSrfs);
	}
    }

    return CagdListReverse(BezSrfs);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Divides the given trimmed surface TSrf to a list of tensor products.     M
*   The result is tiling the original trimmed surface to within machine      M
* precision.								     M
*   The given TSrf surface is recursively divided until the trimming curves  M
* are simple, in which case the trimmed surface is converted into a regular  M
* tensor product surface.  A trimmed surface is considered simple if it has  M
* a single trimming loop that is double monotone with U or V (That is from   M
* the minimal point to the maximal point in U or V we have monotone progress M
* in two separated paths).						     M
*                                                                            *
* PARAMETERS:                                                                M
*   TSrf:         Trimmed surface to decompose into a list of tensor         M
*                 product (and no trimming) surfaces.			     M
*   ComposeE3:    TRUE to compose the tiles into TSrf, FALSE to return the   M
*                 surface tiles in the parametric domain of TSrf.            M
*   OnlyBzrSrfs:  TRUE to force only Bezier tensor products in result.       M
*   HigherOrderTrimmingCurves:  TRUE to keep trimming curves as higher order M
*                 resulting in precise higher order tensor product surfaces. M
*   Eps:          Tolerance of the decomposition.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   A list of regular tensor product surfaces that        M
*                 represents the same region as TSrf, to within machine      M
*                 precision.						     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrimCrvIsParameterizableDomain, Trim2DSrfFromDoubleMonotoneTrimCrv       M
*   TrimSrfCnvrt2BzrRglrSrf						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrimSrfCnvrt2BzrRglrSrf2                                                 M
*****************************************************************************/
CagdSrfStruct *TrimSrfCnvrt2BzrRglrSrf2(const TrimSrfStruct *TSrf,
					int ComposeE3,
					int OnlyBzrSrfs,
					int HigherOrderTrimmingCurves,
					CagdRType Eps)
{
    int OldHigherOrderTrimmingCurves = FALSE;
    CagdSrfStruct *Srf, *RetVal, *Srfs,
        *BzrSrfs = NULL;

    if (HigherOrderTrimmingCurves) {
        /* Do not convert to piecewise linear trimming curves. */
	OldHigherOrderTrimmingCurves =
	    TrimSrfSubdivAtParamForcePiecwiseLinearTrimCrvs(FALSE);
    }
    
    Srfs = TrimSrfCnvrt2TensorProdSrf(TSrf, ComposeE3, Eps);

    if (OnlyBzrSrfs) {
        /* Convert any B-spline surface to Beziers. */
        while (Srfs != NULL) {
	    IRIT_LIST_POP(Srf, Srfs);

	    if (CAGD_IS_BEZIER_SRF(Srf)) {
	        IRIT_LIST_PUSH(Srf, BzrSrfs);
	    }
	    else {
	        CagdSrfStruct
		    *BSrfs = CagdCnvrtBsp2BzrSrf(Srf);

		CagdSrfFree(Srf);
		BzrSrfs = CagdListAppend(BSrfs, BzrSrfs);
	    }
	}

	RetVal = BzrSrfs;
    }
    else
        RetVal = Srfs;

    if (HigherOrderTrimmingCurves) {
        TrimSrfSubdivAtParamForcePiecwiseLinearTrimCrvs(
					       OldHigherOrderTrimmingCurves);
    }

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Divides the given trimmed surface TSrf to a list of tensor product srfs. M
*   The result is tiling the original trimmed surface to within machine      M
* precision.								     M
*   The given TSrf surface is recursively divided until the trimming curves  M
* are simple, in which case the trimmed surface is converted into a regular  M
* tensor product surface.  A trimmed surface is considered simple if it has  M
* a single trimming loop that is double monotone with U or V (That is from   M
* the minimal point to the maximal point in U or V we have monotone progress M
* in two separated paths).						     M
*                                                                            *
* PARAMETERS:                                                                M
*   TSrf:         Trimmed surface to decompose into a list of tensor         M
*                 product (and no trimming) surfaces.			     M
*   ComposeE3:    TRUE to compose the tiles into TSrf, FALSE to return the   M
*                 surface tiles in the parametric domain of TSrf.            M
*   Eps:          Tolerance of the decomposition.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   A list of regular tensor product surfaces that        M
*                 represents the same region as TSrf, to within machine      M
*                 precision.						     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrimCrvIsParameterizableDomain, Trim2DSrfFromDoubleMonotoneTrimCrv       M
*   TrimSrfCnvrt2BzrRglrSrf, TrimSrfCnvrt2BzrRglrSrf2			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrimSrfCnvrt2TensorProdSrf                                               M
*****************************************************************************/
CagdSrfStruct *TrimSrfCnvrt2TensorProdSrf(const TrimSrfStruct *TSrf,
					  int ComposeE3,
					  CagdRType Eps)
{
    CagdRType SubdivVal;
    CagdSrfDirType SubdivDir;
    CagdBBoxStruct BBox1;
    CagdSrfStruct *Srfs1, *Srfs2,
        *Srf = TSrf -> Srf;
    TrimSrfStruct *TSrf1, *TSrf2;
    TrimCrvStruct
	*TrimCrvList = TrimChainTrimmingCurves2Loops(TSrf -> TrimCrvList);

#   ifdef DEBUG_TRIM_VIEW_DOMAIN
    {
	IRIT_STATIC_DATA int
	    PrgmIO = -1;
	CagdRType U1, U2, V1, V2;
	TrimSrfStruct *TSrfNext = TSrf -> Pnext;

	CagdSrfDomain(Srf, &U1, &U2, &V1, &V2);
	fprintf(stderr, "%f %f   %f %f\n", U1, U2, V1, V2);
	if (PrgmIO < 0) {
	    IPSocSrvrInit();    /* Initialize the listen socket for clients. */
#	    if defined(__WINNT__) && defined(_WIN64)
		PrgmIO = IPSocExecAndConnect(getenv("IRIT_DISPLAY64"),
#	    else
		PrgmIO = IPSocExecAndConnect(getenv("IRIT_DISPLAY"),
#	    endif /* __WINNT__ && _WIN64 */
					     getenv("IRIT_BIN_IPC") != NULL);
	    assert(PrgmIO >= 0);
	}
	((TrimSrfStruct *) TSrf) -> Pnext = NULL;
	IPSocWriteOneObject(PrgmIO, IPGenTRIMSRFObject((TrimSrfStruct *) TSrf));
	((TrimSrfStruct *) TSrf) -> Pnext = TSrfNext;
    }
#   endif /* DEBUG_TRIM_VIEW_DOMAIN */

    CagdCrvBBox(TrimCrvList -> TrimCrvSegList -> UVCrv, &BBox1);

    if (TrimCrvList == NULL) {
        /* This patch is not trimmed - just return it... */
        TrimCrvFreeList(TrimCrvList);
        return CagdSrfCopy(TSrf -> Srf);
    }
    else if (Srf -> UOrder != Srf -> ULength) {
        /* Cannot compose into B-spline surfaces - divide TSrf first. */
        SubdivVal = Srf -> UKnotVector[(Srf -> ULength + Srf -> UOrder) >> 1];
	SubdivDir = CAGD_CONST_U_DIR;
    }
    else if (Srf -> VOrder != Srf -> VLength) {
        /* Cannot compose into B-spline surfaces - divide TSrf first. */
        SubdivVal = Srf -> VKnotVector[(Srf -> VLength + Srf -> VOrder) >> 1];
	SubdivDir = CAGD_CONST_V_DIR;
    }
    else if (TrimCrvList -> Pnext != NULL) {
        CagdBBoxStruct BBox2;

        /* More than one trimming loop */
	CagdCrvBBox(TrimCrvList -> Pnext -> TrimCrvSegList -> UVCrv, &BBox2);

	if (BBox1.Max[0] < BBox2.Min[0]) {
	    /* Trimming curves do not overlap in X. */
	    SubdivVal = (BBox1.Max[0] + BBox2.Min[0]) * 0.5;
	    SubdivDir = CAGD_CONST_U_DIR;
	}
	else if (BBox2.Max[0] < BBox1.Min[0]) {
	    /* Trimming curves do not overlap in X. */
	    SubdivVal = (BBox2.Max[0] + BBox1.Min[0]) * 0.5;
	    SubdivDir = CAGD_CONST_U_DIR;
	}
	else if (BBox1.Max[1] < BBox2.Min[1]) {
	    /* Trimming curves do not overlap in Y. */
	    SubdivVal = (BBox1.Max[1] + BBox2.Min[1]) * 0.5;
	    SubdivDir = CAGD_CONST_V_DIR;
	}
	else if (BBox2.Max[1] < BBox1.Min[1]) {
	    /* Trimming curves do not overlap in Y. */
	    SubdivVal = (BBox2.Max[1] + BBox1.Min[1]) * 0.5;
	    SubdivDir = CAGD_CONST_V_DIR;
	}
	else {
	    /* Select middle of smaller trimming curve. */
	    if (BBox1.Max[0] - BBox1.Min[0] < BBox2.Max[0] - BBox2.Min[0])
	        SubdivVal = (BBox1.Max[0] + BBox1.Min[0]) * 0.5;
	    else
	        SubdivVal = (BBox2.Max[0] + BBox2.Min[0]) * 0.5;
	    SubdivDir = CAGD_CONST_U_DIR;
	}
    }
    else {
        /* Make sure we have domain [0, 1]^2 before converting to a Bezier. */
        CagdSrfStruct *UVSrfs, *UVSrf;
	TrimSrfStruct
	    *TSrfBzrDmn = TrimSrfSetDomain(TSrf, 0.0, 1.0, 0.0, 1.0);

	TrimCrvFreeList(TrimCrvList);
	TrimCrvList = TrimChainTrimmingCurves2Loops(TSrfBzrDmn -> TrimCrvList);

        if ((UVSrfs = Symb2DCrvParameterizeDomain(TrimCrvList -> TrimCrvSegList
						                     -> UVCrv,
						  Eps)) != NULL) {
	    CagdSrfStruct *TSrf, *OrientedUVSrf,
	        *RetSrfs = NULL,
	        *BzrSrf = CagdCnvrtBsp2BzrSrf(TSrfBzrDmn -> Srf);

	    assert(BzrSrf -> Pnext == NULL &&    /* Is a single Bezier srf. */
		   TrimCrvList != NULL &&        /* Has one trimming curve. */
		   TrimCrvList -> Pnext == NULL);

	    for (UVSrf = UVSrfs; UVSrf != NULL; UVSrf = UVSrf -> Pnext) {
		OrientedUVSrf = TrimOrientUVReparamSrf(UVSrf);
		if (ComposeE3) {
		    TSrf = SymbComposeSrfSrf(BzrSrf, OrientedUVSrf);
		    CagdSrfFree(OrientedUVSrf);
		    IRIT_LIST_PUSH(TSrf, RetSrfs);
		}
		else {
		    IRIT_LIST_PUSH(OrientedUVSrf, RetSrfs);
		}
	    }

	    TrimSrfFree(TSrfBzrDmn);
	    TrimCrvFreeList(TrimCrvList);
	    CagdSrfFree(BzrSrf);
	    CagdSrfFreeList(UVSrfs);

	    return RetSrfs;
	}
	else {  /* Too complex of a trimming curves - subdivided in middle. */
	    if (BBox1.Max[1] - BBox1.Min[1] > BBox1.Max[0] - BBox1.Min[0]) {
	        /* Subdiv. in V. */
	        SubdivVal = (BBox1.Min[1] + BBox1.Max[1]) * 0.5;
		SubdivDir = CAGD_CONST_V_DIR;
	    }
	    else {
	        /* Subdiv. in U. */
	        SubdivVal = (BBox1.Min[0] + BBox1.Max[0]) * 0.5;
		SubdivDir = CAGD_CONST_U_DIR;
	    }
	}

	TrimSrfFree(TSrfBzrDmn);
    }

    TSrf1 = TrimSrfSubdivAtParam((TrimSrfStruct *) TSrf,
				 SubdivVal, SubdivDir);
    TSrf2 = TSrf1 -> Pnext;

    Srfs1 = TrimSrfCnvrt2TensorProdSrf(TSrf1, ComposeE3, Eps);
    Srfs2 = TSrf2 != NULL ? TrimSrfCnvrt2TensorProdSrf(TSrf2, ComposeE3, Eps)
			  : NULL;

    TrimSrfFreeList(TSrf1);
    TrimCrvFreeList(TrimCrvList);

    return CagdListAppend(Srfs1, Srfs2);
}
