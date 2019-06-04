/******************************************************************************
* CagdBSum.c - Boolean sum surface constructor out of given four curves.      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Sep. 91.					      *
******************************************************************************/

#include "cagd_loc.h"

static CagdCrvStruct *CagdDivideAtC1DiscontsUptoFourBoundaries(
						   const CagdCrvStruct *Crvs);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a Boolean sum surface using the four provided boundary curves.  M
*   Curve's end points must meet at the four surface corners if surface      M
* boundary are to be identical to the four given curves.	             M
*                           bottom                                           V
*                      +----->-->-----+                                      V
*                      |              |                                      V
*                      ^              ^                                      V
*                      |              |                                      V
*                 left |              | right                                V
*                      |              |                                      V
*                      ^              ^                                      V
*                      |              |                                      V
*                      +----->-->-----+                                      V
*                            top                                             V
*                                                                            *
* PARAMETERS:                                                                M
*   CCrvLeft:   Left boundary curve of Boolean sum surface to be created.    M
*   CCrvRight:  Right boundary curve of Boolean sum surface to be created.   M
*   CCrvTop:    Top boundary curve of Boolean sum surface to be created.     M
*   CCrvBottom: Bottom boundary curve of Boolean sum surface to be created.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:  A Boolean sum surface constructed using given four     M
*                     curves.                                                M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdBoolSumSrf, Boolean sum, surface constructors                        M
*****************************************************************************/
CagdSrfStruct *CagdBoolSumSrf(const CagdCrvStruct *CCrvLeft,
			      const CagdCrvStruct *CCrvRight,
			      const CagdCrvStruct *CCrvTop,
			      const CagdCrvStruct *CCrvBottom)
{
    int i, j;
    CagdCrvStruct *Crv1, *Crv2, *CrvLeft, *CrvRight, *CrvTop, *CrvBottom;
    CagdSrfStruct *Ruled1, *Ruled2, *Ruled3, *Srf;
    CagdPtStruct Pt1, Pt2;
    CagdRType **SrfPoints, **Ruled1Pts, **Ruled2Pts, **Ruled3Pts;

    if (CAGD_IS_PERIODIC_CRV(CCrvLeft) ||
	CAGD_IS_PERIODIC_CRV(CCrvRight) ||
	CAGD_IS_PERIODIC_CRV(CCrvTop) ||
	CAGD_IS_PERIODIC_CRV(CCrvBottom)) {
	CAGD_FATAL_ERROR(CAGD_ERR_PERIODIC_NO_SUPPORT);
	return NULL;
    }

    if (CAGD_IS_BSPLINE_CRV(CCrvLeft) && !BspCrvHasOpenEC(CCrvLeft))
        CrvLeft = CagdCnvrtFloat2OpenCrv(CCrvLeft);
    else
        CrvLeft = CagdCrvCopy(CCrvLeft);
    if (CAGD_IS_BSPLINE_CRV(CCrvRight) && !BspCrvHasOpenEC(CCrvRight))
        CrvRight = CagdCnvrtFloat2OpenCrv(CCrvRight);
    else
        CrvRight = CagdCrvCopy(CCrvRight);
    if (CAGD_IS_BSPLINE_CRV(CCrvTop) && !BspCrvHasOpenEC(CCrvTop))
        CrvTop = CagdCnvrtFloat2OpenCrv(CCrvTop);
    else
        CrvTop = CagdCrvCopy(CCrvTop);
    if (CAGD_IS_BSPLINE_CRV(CCrvBottom) && !BspCrvHasOpenEC(CCrvBottom))
        CrvBottom = CagdCnvrtFloat2OpenCrv(CCrvBottom);
    else
        CrvBottom = CagdCrvCopy(CCrvBottom);

    if (CAGD_IS_BSPLINE_CRV(CrvLeft))
        BspKnotAffineTrans2(CrvLeft -> KnotVector,
			    CrvLeft -> Order + CrvLeft -> Length, 0.0, 1.0);
    if (CAGD_IS_BSPLINE_CRV(CrvRight))
        BspKnotAffineTrans2(CrvRight -> KnotVector,
			    CrvRight -> Order + CrvRight -> Length, 0.0, 1.0);
    if (CAGD_IS_BSPLINE_CRV(CrvTop))
        BspKnotAffineTrans2(CrvTop -> KnotVector,
			    CrvTop -> Order + CrvTop -> Length, 0.0, 1.0);
    if (CAGD_IS_BSPLINE_CRV(CrvBottom))
        BspKnotAffineTrans2(CrvBottom -> KnotVector,
			    CrvBottom -> Order + CrvBottom -> Length, 0.0, 1.0);

    /* The Left-Right and Top-Bottom curves should share same point/curve    */
    /* type as well as same order & knot vector (if B-spline representation).*/
    CagdMakeCrvsCompatible(&CrvLeft, &CrvRight, TRUE, TRUE);
    CagdMakeCrvsCompatible(&CrvTop, &CrvBottom, TRUE, TRUE);

    /* The Left-Right and Top-Bottom pairs must share same point/curve type. */
    CagdMakeCrvsCompatible(&CrvLeft, &CrvTop, FALSE, FALSE);
    CagdMakeCrvsCompatible(&CrvLeft, &CrvBottom, FALSE, FALSE);
    CagdMakeCrvsCompatible(&CrvRight, &CrvTop, FALSE, FALSE);
    CagdMakeCrvsCompatible(&CrvRight, &CrvBottom, FALSE, FALSE);

#ifdef DEBUG
    {
        /* Verify the corners. */
        CagdCoerceToE3(Pt1.Pt, CrvLeft -> Points, 0, CrvLeft -> PType);
        CagdCoerceToE3(Pt2.Pt, CrvTop -> Points, 0, CrvTop -> PType);
	assert(IRIT_PT_APX_EQ(Pt1.Pt, Pt2.Pt));

        CagdCoerceToE3(Pt1.Pt, CrvLeft -> Points, CrvLeft -> Length - 1,
		                                            CrvLeft -> PType);
        CagdCoerceToE3(Pt2.Pt, CrvBottom -> Points, 0, CrvBottom -> PType);
	assert(IRIT_PT_APX_EQ(Pt1.Pt, Pt2.Pt));

        CagdCoerceToE3(Pt1.Pt, CrvRight -> Points, 0, CrvRight -> PType);
        CagdCoerceToE3(Pt2.Pt, CrvTop -> Points, CrvTop -> Length - 1,
						          CrvBottom -> PType);
	assert(IRIT_PT_APX_EQ(Pt1.Pt, Pt2.Pt));

        CagdCoerceToE3(Pt1.Pt, CrvRight -> Points, CrvRight -> Length - 1,
		                                           CrvRight -> PType);
        CagdCoerceToE3(Pt2.Pt, CrvBottom -> Points, CrvBottom -> Length - 1,
							     CrvTop -> PType);
	assert(IRIT_PT_APX_EQ(Pt1.Pt, Pt2.Pt));
    }
#endif /* DEBUG */

    /* Now that the curves are in the right representation, form surface(s). */
    /* The two ruled surface between the respective curves, in right orders: */
    Ruled1 = CagdRuledSrf(CrvLeft, CrvRight, 2, 2);
    Ruled2 = CagdRuledSrf(CrvTop, CrvBottom, 2, 2);
    Srf = CagdSrfReverse2(Ruled2);
    CagdSrfFree(Ruled2);
    Ruled2 = Srf;
    CagdMakeSrfsCompatible(&Ruled1, &Ruled2, TRUE, TRUE, TRUE, TRUE);

    /* The ruled surface between the four corner points in the right orders. */
    CagdCoerceToE3(Pt1.Pt, CrvLeft -> Points, 0, CrvLeft -> PType);
    CagdCoerceToE3(Pt2.Pt, CrvLeft -> Points, CrvLeft -> Length - 1,
							CrvLeft -> PType);
    Crv1 = CagdMergePtPt(&Pt1, &Pt2);

    CagdCoerceToE3(Pt1.Pt, CrvRight -> Points, 0, CrvRight -> PType);
    CagdCoerceToE3(Pt2.Pt, CrvRight -> Points, CrvRight -> Length - 1,
							CrvRight -> PType);
    Crv2 = CagdMergePtPt(&Pt1, &Pt2);

    /* Should not change CrvLeft/Right only Crv1/2: */
    Ruled3 = CagdRuledSrf(Crv1, Crv2, 2, 2);
    if (CAGD_IS_BSPLINE_SRF(Ruled3)) {
        BspKnotAffineTrans2(Ruled3 -> UKnotVector,
			    Ruled3 -> UOrder + Ruled3 -> ULength, 0.0, 1.0);
	BspKnotAffineTrans2(Ruled3 -> VKnotVector,
			    Ruled3 -> VOrder + Ruled3 -> VLength, 0.0, 1.0);
    }
    CagdMakeSrfsCompatible(&Ruled1, &Ruled3, TRUE, TRUE, TRUE, TRUE);

    CagdCrvFree(Crv1);
    CagdCrvFree(Crv2);

    /* The boolean sum is equal to Ruled1 + Ruled2 - Ruled3. This boolean    */
    /* sum as computed is not exactly as defined in the literature for non   */
    /* uniform Bsplines since the ruled surfaces are computed with uniform   */
    /* distribution along the other axis even if it is non uniform.	     */
    if (CrvRight -> GType == CAGD_CBSPLINE_TYPE) {
	Srf = BspSrfNew(Ruled1 -> ULength, Ruled1 -> VLength,
			Ruled1 -> UOrder, Ruled1 -> VOrder, Ruled1 -> PType);
	BspKnotCopy(Srf -> UKnotVector, Ruled1 -> UKnotVector,
		    Ruled1 -> ULength + Ruled1 -> UOrder);
	BspKnotCopy(Srf -> VKnotVector, Ruled1 -> VKnotVector,
		    Ruled1 -> VLength + Ruled1 -> VOrder);
    }
    else if (CrvRight -> GType == CAGD_CBEZIER_TYPE)
	Srf = BzrSrfNew(Ruled1 -> ULength, Ruled1 -> VLength, Ruled1 -> PType);
    else
	CAGD_FATAL_ERROR(CAGD_ERR_POWER_NO_SUPPORT);
    SrfPoints = Srf -> Points;

    Ruled1Pts = Ruled1 -> Points;
    Ruled2Pts = Ruled2 -> Points;
    Ruled3Pts = Ruled3 -> Points;

    for (i = !CAGD_IS_RATIONAL_SRF(Srf);
	 i <= CAGD_NUM_OF_PT_COORD(Srf -> PType);
	 i++) {
	CagdRType
	    *Ruled1PtsPtr = Ruled1Pts[i],
	    *Ruled2PtsPtr = Ruled2Pts[i],
	    *Ruled3PtsPtr = Ruled3Pts[i],
	    *SrfPointsPtr = SrfPoints[i];

	 for (j = Srf -> ULength * Srf -> VLength; j > 0; j--)
	     *SrfPointsPtr++ =
		 *Ruled1PtsPtr++ + *Ruled2PtsPtr++ - *Ruled3PtsPtr++;
    }

    CagdSrfFree(Ruled1);
    CagdSrfFree(Ruled2);
    CagdSrfFree(Ruled3);

    CagdCrvFree(CrvTop);
    CagdCrvFree(CrvRight);
    CagdCrvFree(CrvBottom);
    CagdCrvFree(CrvLeft);

    CAGD_SET_GEOM_TYPE(Srf, CAGD_GEOM_BOOL_SUM);

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Constructs a Boolean sum surface using the single boundary curve.	     M
*   The curve is subdivided into four, equally spaced in parameter space,    M
* sub-regions which are used as the four curves to the Boolean sum           M
* constructor. See CagdBoolSumSrf.                                           M
*                                                                            *
* PARAMETERS:                                                                M
*   BndryCrv:   To be subdivided into four curves for a Boolean sum          M
*               construction.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:  A Boolean sum surface constructed using given curve    M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdOneBoolSumSrf, Boolean sum, surface constructors                     M
*****************************************************************************/
CagdSrfStruct *CagdOneBoolSumSrf(const CagdCrvStruct *BndryCrv)
{
    CagdRType TMin, TMax;
    CagdCrvStruct *Crvs, *TCrv, *CrvLeft, *CrvRight, *CrvTop, *CrvBottom;
    CagdSrfStruct *BoolSumSrf;

    CagdCrvDomain(BndryCrv, &TMin, &TMax);

    Crvs = CagdCrvSubdivAtParam(BndryCrv, TMin * 0.5 + TMax * 0.5);
    CrvLeft = CagdCrvSubdivAtParam(Crvs, TMin * 0.75 + TMax * 0.25);
    CrvRight = CagdCrvSubdivAtParam(Crvs -> Pnext, TMin * 0.25 + TMax * 0.75);

    CagdCrvFreeList(Crvs);

    CrvBottom = CrvLeft -> Pnext;
    CrvLeft -> Pnext = NULL;
    CrvTop = CrvRight -> Pnext;
    CrvRight -> Pnext = NULL;

    /* Reverse CrvBottom and CrvLeft: */
    TCrv = CagdCrvReverse(CrvTop);
    CagdCrvFree(CrvTop);
    CrvTop = TCrv;

    TCrv = CagdCrvReverse(CrvRight);
    CagdCrvFree(CrvRight);
    CrvRight = TCrv;

    BoolSumSrf = CagdBoolSumSrf(CrvLeft, CrvRight, CrvTop, CrvBottom);

    CagdCrvFree(CrvTop);
    CagdCrvFree(CrvRight);
    CagdCrvFree(CrvBottom);
    CagdCrvFree(CrvLeft);

    CAGD_SET_GEOM_TYPE(BoolSumSrf, CAGD_GEOM_BOOL_SUM);

    return BoolSumSrf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Properly reorder given curves of one closed loops, in place.  That is,    M
* input is assumed to define a complete loop where one curve ends where	     M
* another begins.  Compare end points in R^3 (XYZ).			     M
*    Input curves can be in arbitrary order and even partially reversed.     M
*    Orientation of loop will be following first curve in the list.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   UVCrvs:  Curves forming one loop to reorder so end point of one curve    M
*            is the beginning of the next curve.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    Curves properly reordered in loop.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdSrfFromNBndryCrvs					             M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdReorderCurvesInLoop					             M
*****************************************************************************/
CagdCrvStruct *CagdReorderCurvesInLoop(CagdCrvStruct *UVCrvs)
{
    CagdPType Pt1End;
    CagdCrvStruct *UVCrv, *UVCrnt,
	*UVRes = NULL;

    UVRes = UVCrnt = UVCrvs;
    UVCrvs = UVCrvs -> Pnext;
    UVRes -> Pnext = NULL;

    /* Evaluate the end point of the current curve. */
    CagdCoerceToE3(Pt1End, UVCrnt -> Points, UVCrnt -> Length - 1,
		   UVCrnt -> PType);

    while (UVCrvs != NULL) {
        CagdBType
	    NextRvrsrd = FALSE;
	CagdRType t,
	    Tol = IRIT_INFNTY;
	CagdPType Pt2Start, Pt2End;
	CagdCrvStruct
	    *NextCrv = NULL;

	for (UVCrv = UVCrvs; UVCrv != NULL; UVCrv = UVCrv -> Pnext) {
	    /* Evaluate the end points of the second curve. */
	    CagdCoerceToE3(Pt2Start, UVCrv -> Points, 0, UVCrv -> PType);
	    CagdCoerceToE3(Pt2End, UVCrv -> Points, UVCrv -> Length - 1,
			   UVCrv -> PType);

	    if ((t = IRIT_PT_PT_DIST_SQR(Pt1End, Pt2Start)) < Tol) {
	        NextCrv = UVCrv;
		NextRvrsrd = FALSE;
		Tol = t;
	    }
	    if ((t = IRIT_PT_PT_DIST_SQR(Pt1End, Pt2End)) < Tol) {
	        NextCrv = UVCrv;
		NextRvrsrd = TRUE;
		Tol = t;
	    }
	}
	assert(NextCrv != NULL);

	/* Remove NextCrv from list */
	if (NextCrv == UVCrvs) {
	    UVCrvs = UVCrvs -> Pnext;
	}
	else {
	    while (UVCrvs -> Pnext != NextCrv)
	        UVCrvs = UVCrvs -> Pnext;
	    UVCrvs -> Pnext = UVCrvs -> Pnext -> Pnext;
	}

	/* And chain in order. */
	if (NextRvrsrd) {
	    UVCrnt -> Pnext = CagdCrvReverse(NextCrv);
	    CagdCoerceToE3(Pt1End, NextCrv -> Points, 0, NextCrv -> PType);
	    CagdCrvFree(NextCrv);
	}
	else {
	    UVCrnt -> Pnext = NextCrv;
	    CagdCoerceToE3(Pt1End, NextCrv -> Points, NextCrv -> Length - 1,
			   NextCrv -> PType);
	}

	UVCrnt = UVCrnt -> Pnext;
    }

    return UVRes;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Break the given, less than four, curves into four (or a bit more) at C1  *
* discontinuities, if found any, seeking largest C1 discontinuities first.   *
*                                                                            *
* PARAMETERS:                                                                *
*   Crvs:   Less than four curves to break at C1 discontinuities into four.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:   List of four or a bit more curves, or NULL if failed. *
*****************************************************************************/
static CagdCrvStruct *CagdDivideAtC1DiscontsUptoFourBoundaries(
						    const CagdCrvStruct *Crvs)
{
    int n = CagdListLength(Crvs);
    CagdRType
        CosAngle = -0.9;

    /* Try to break curves at largest C1 discontinuities first. */
    while (TRUE) {
        int NewN = 0;
	const CagdCrvStruct *Crv;
        CagdCrvStruct *DivCrv,
	    *NewCrvs = NULL;

	for (Crv = Crvs; Crv != NULL; Crv = Crv -> Pnext) {
	    DivCrv = CagdCrvSubdivAtAllC1Discont(Crv, TRUE, CosAngle);
	    NewN += CagdListLength(DivCrv) - 1;       /* How many new segs? */

	    NewCrvs = CagdListAppend(NewCrvs, DivCrv);
	    if (n + NewN >= 4) {
	        /* Have enough - copy the rest and quit. */
	        NewCrvs = CagdListAppend(NewCrvs,
					 CagdCrvCopyList(Crv -> Pnext));
		return NewCrvs;
	    }
	}

	/* If we are here - we did not find enough C1 discontinuities -     */
	/* shrink the allowed angle.					    */
	CosAngle += 0.1;
	if (CosAngle >= 0.95) {		                 /* Last iteration. */
	    if (NewN > 0)
		return NewCrvs;
	    else {
		CagdCrvFreeList(NewCrvs);
		return NULL;
	    }
	}
	else
	    CagdCrvFreeList(NewCrvs);
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Builds tensor product surfaces that spans the given list of surface      M
* boundary curves in one closed loop.					     M
*   Can be 1 to 6 input boundary curves, and 1 to 2 surfaces are returned.   M
*                                                                            *
* PARAMETERS:                                                                M
*   Crvs:     To build 1/2 surfaces with these curves as boundaries.         M
*   MinimizeSize:   If true, minimize the size of the output, on expense of  M
*             accuracy.							     N
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   One or two planar surfaces spanning the curves.       M
*                                                                            *
* SEE ALSO:                                                                  M
*   CagdReorderCurvesInLoop					             M
*                                                                            *
* KEYWORDS:                                                                  M
*   CagdSrfFromNBndryCrvs	                                             M
*****************************************************************************/
CagdSrfStruct *CagdSrfFromNBndryCrvs(const CagdCrvStruct *Crvs,
				     CagdBType MinimizeSize)
{
    int n, l;
    CagdRType TMin, TMax;
    CagdPtStruct Pt1, Pt2;
    CagdCrvStruct *T4Crvs, *TCrv, *UVCrvs,
         *Crv1, *Crv2, *Crv3, *Crv4, *Crv23;
    CagdSrfStruct *UVSrf1,
        *UVSrf2 = NULL;

    UVCrvs = CagdReorderCurvesInLoop(CagdCrvCopyList(Crvs));

    n = CagdListLength(UVCrvs);

    /* Convert loop into 4 sided domain(s). */
    switch (n) {
        case 1:
	    /* Try to break at C^1Discontinuities and see if we get a      */
	    /* reasonable division.					   */
	    T4Crvs = CagdDivideAtC1DiscontsUptoFourBoundaries(UVCrvs);
	    l = T4Crvs ? CagdListLength(T4Crvs) : 0;
	    if (l >= 2 && l < 7) {
	        UVSrf1 = CagdSrfFromNBndryCrvs(T4Crvs, MinimizeSize);
		CagdCrvFreeList(T4Crvs);
	    }
	    else {
		UVSrf1 = CagdOneBoolSumSrf(UVCrvs);
	    }
	    break;
        case 2:
	    T4Crvs = CagdDivideAtC1DiscontsUptoFourBoundaries(UVCrvs);
	    l = T4Crvs ? CagdListLength(T4Crvs) : 0;
	    if (l >= 4 && l < 7) {
	        UVSrf1 = CagdSrfFromNBndryCrvs(T4Crvs, MinimizeSize);
		CagdCrvFreeList(T4Crvs);
	    }
	    else if (l > 2) {
	        CagdCrvFreeList(T4Crvs);

		UVSrf1 = CagdOneBoolSumSrf(UVCrvs);
	    }
	    else {
	        /* Create a ruled surface between the two curves. */
		TCrv = CagdCrvReverse(UVCrvs -> Pnext);

		UVSrf1 = CagdRuledSrf(UVCrvs, TCrv, 2, 2);

		CagdCrvFree(TCrv);
	    }
	    break;
        case 3:
	    T4Crvs = CagdDivideAtC1DiscontsUptoFourBoundaries(UVCrvs);
	    l = T4Crvs ? CagdListLength(T4Crvs) : 0;
	    if (l >= 4) {
	        UVSrf1 = CagdSrfFromNBndryCrvs(T4Crvs, MinimizeSize);
		CagdCrvFreeList(T4Crvs);
	    }
	    else {
	        /* Create a degenerated loop with one boundary singular. */
	        Crv1 = CagdCrvCopy(UVCrvs);
		Crv2 = CagdCrvCopy(UVCrvs -> Pnext);
		Crv3 = CagdCrvReverse(UVCrvs -> Pnext -> Pnext);
		CagdCoerceToE3(Pt1.Pt, Crv1 -> Points, 0, Crv1 -> PType);
		Crv4 = CagdMergePtPt(&Pt1, &Pt1);

		UVSrf1 = CagdBoolSumSrf(Crv1, Crv3, Crv4, Crv2);

		CagdCrvFree(Crv1);
		CagdCrvFree(Crv2);
		CagdCrvFree(Crv3);
		CagdCrvFree(Crv4);
	    }
	    break;
        case 4:
	    Crv1 = CagdCrvCopy(UVCrvs);
	    Crv2 = CagdCrvCopy(UVCrvs -> Pnext);
	    Crv3 = CagdCrvReverse(UVCrvs -> Pnext -> Pnext);
	    Crv4 = CagdCrvReverse(UVCrvs -> Pnext -> Pnext -> Pnext);

	    if (CagdEstimateCrvCollinearity(Crv1) < IRIT_EPS &&
		CagdEstimateCrvCollinearity(Crv3) < IRIT_EPS) {
	        /* Simply rule a surface between Crv2 and Crv4. */
	        UVSrf1 = CagdRuledSrf(Crv2, Crv4, 2, 2);
	    }
	    else if (CagdEstimateCrvCollinearity(Crv2) < IRIT_EPS &&
		     CagdEstimateCrvCollinearity(Crv4) < IRIT_EPS) {
	        /* Simply rule a surface between Crv1 and Crv3. */
	        UVSrf1 = CagdRuledSrf(Crv1, Crv3, 2, 2);
	    }
	    else {
	        if (MinimizeSize) {
		    /* Make the knot sequences the same, if possible. */
		    if (CAGD_IS_BSPLINE_CRV(Crv1) &&
			CAGD_IS_BSPLINE_CRV(Crv3) &&
			Crv1 -> Order == Crv3 -> Order &&
			Crv1 -> Length == Crv3 -> Length) {/* Force same KV.*/
		        CAGD_GEN_COPY(Crv1 -> KnotVector,
				      Crv3 -> KnotVector,
				      sizeof(CagdRType) * (Crv1 -> Order +
							   Crv1 -> Length));
		    }
		    /* Make the knot sequences the same, if possible. */
		    if (CAGD_IS_BSPLINE_CRV(Crv2) &&
			CAGD_IS_BSPLINE_CRV(Crv4) &&
			Crv2 -> Order == Crv4 -> Order &&
			Crv2 -> Length == Crv4 -> Length) {/* Force same KV.*/
		        CAGD_GEN_COPY(Crv2 -> KnotVector,
				      Crv4 -> KnotVector,
				      sizeof(CagdRType) * (Crv2 -> Order +
							   Crv2 -> Length));
		    }
		}

		UVSrf1 = CagdBoolSumSrf(Crv1, Crv3, Crv4, Crv2);
	    }

	    CagdCrvFree(Crv1);
	    CagdCrvFree(Crv2);
	    CagdCrvFree(Crv3);
	    CagdCrvFree(Crv4);
	    break;
        case 5:
	    /* Build first 4 sided loop: with crv1, crv2 & half of crv3. */
	    Crv1 = CagdCrvCopy(UVCrvs);
	    Crv2 = CagdCrvCopy(UVCrvs -> Pnext);
	    CagdCrvDomain(UVCrvs -> Pnext -> Pnext, &TMin, &TMax);
	    Crv23 = CagdCrvSubdivAtParam(UVCrvs -> Pnext -> Pnext,
					 (TMin + TMax) * 0.5);
	    Crv3 = CagdCrvReverse(Crv23);

	    CagdCoerceToE3(Pt1.Pt, Crv1 -> Points, 0, Crv1 -> PType);
	    CagdCoerceToE3(Pt2.Pt, Crv3 -> Points, 0, Crv3 -> PType);
	    Crv4 = CagdMergePtPt(&Pt1, &Pt2);

	    UVSrf1 = CagdBoolSumSrf(Crv1, Crv3, Crv4, Crv2);

	    CagdCrvFree(Crv1);
	    CagdCrvFree(Crv2);
	    CagdCrvFree(Crv3);
	    CagdCrvFree(Crv4);

	    /* Build second 4 sided loop: with half of crv3 and crv4/5. */ 
	    Crv1 = CagdCrvCopy(Crv23 -> Pnext);
	    Crv2 = CagdCrvCopy(UVCrvs -> Pnext -> Pnext -> Pnext);
	    Crv3 = CagdCrvReverse(UVCrvs -> Pnext -> Pnext -> Pnext -> Pnext);

	    CagdCoerceToE3(Pt1.Pt, Crv1 -> Points, 0, Crv1 -> PType);
	    CagdCoerceToE3(Pt2.Pt, Crv3 -> Points, 0, Crv3 -> PType);
	    Crv4 = CagdMergePtPt(&Pt1, &Pt2);

	    UVSrf2 = CagdBoolSumSrf(Crv1, Crv3, Crv4, Crv2);

	    CagdCrvFree(Crv1);
	    CagdCrvFree(Crv2);
	    CagdCrvFree(Crv3);
	    CagdCrvFree(Crv4);

	    CagdCrvFreeList(Crv23);
	    break;
        case 6:
	    /* Build first 4 sided loop: with crv1/2/3. */ 
	    Crv1 = CagdCrvCopy(UVCrvs);
	    Crv2 = CagdCrvCopy(UVCrvs -> Pnext);
	    Crv3 = CagdCrvReverse(UVCrvs -> Pnext -> Pnext);

	    CagdCoerceToE3(Pt1.Pt, Crv1 -> Points, 0, Crv1 -> PType);
	    CagdCoerceToE3(Pt2.Pt, Crv3 -> Points, 0, Crv3 -> PType);
	    Crv4 = CagdMergePtPt(&Pt1, &Pt2);

	    UVSrf1 = CagdBoolSumSrf(Crv1, Crv3, Crv4, Crv2);

	    CagdCrvFree(Crv1);
	    CagdCrvFree(Crv2);
	    CagdCrvFree(Crv3);
	    CagdCrvFree(Crv4);

	    /* Build second 4 sided loop: with crv4/5/6. */ 
	    Crv1 = CagdCrvCopy(TCrv = UVCrvs -> Pnext -> Pnext -> Pnext);
	    Crv2 = CagdCrvCopy(TCrv -> Pnext);
	    Crv3 = CagdCrvReverse(TCrv -> Pnext -> Pnext);

	    CagdCoerceToE3(Pt1.Pt, Crv1 -> Points, 0, Crv1 -> PType);
	    CagdCoerceToE3(Pt2.Pt, Crv3 -> Points, 0, Crv3 -> PType);
	    Crv4 = CagdMergePtPt(&Pt1, &Pt2);

	    UVSrf2 = CagdBoolSumSrf(Crv1, Crv3, Crv4, Crv2);

	    CagdCrvFree(Crv1);
	    CagdCrvFree(Crv2);
	    CagdCrvFree(Crv3);
	    CagdCrvFree(Crv4);
	    break;
       default:
	    IRIT_INFO_MSG_PRINTF("Loop too complex - contains more than 6 segments.\n");
	    UVSrf1 = UVSrf2 = NULL;
	    break;
    }

    CagdCrvFreeList(UVCrvs);

    if (UVSrf2 != NULL)
        ((CagdSrfStruct *) CagdListLast(UVSrf1)) -> Pnext = UVSrf2;
    return UVSrf1;
}
