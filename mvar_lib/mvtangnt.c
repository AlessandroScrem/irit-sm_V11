/******************************************************************************
* MvTangnt.c - Compute bi-tangents and tri-tangents of freeform surfaces.     *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, July 97.					      *
******************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "mvar_loc.h"
#include "inc_irit/user_lib.h"
#include "inc_irit/allocate.h"

#define SUBDIV_REL_DISTANCE_TOL		3
#define MVAR_CIRC_TAN_2CRVS_NUMER_TOL	1e-10

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes bi-tangents of freeform bivariate.                              M
* Let,									     M
*									     M
*	DMV = MV1(u, v) - MV2(r, s)					     V
*									     M
* then, computed the simultaneous solution of the following three equations: M
*									     M
*   d MV1   d MV1  d MV2	                 			     V
* < ----- x -----, ----- > = 0,				                     V
*     du      dv     dr				                             V
*									     M
*   d MV1   d MV1  d MV2	                 			     V
* < ----- x -----, ----- > = 0,				                     V
*     du      dv     ds				                             V
*									     M
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV > = 0.				                     V
*     du      dv					                     V
*									     M
*                                                                            *
* PARAMETERS:                                                                M
*   CMV1, CMV2:   The two multivariates to compute the bi-tangents for.      M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   NumericTol:   Numeric tolerance of the numeric stage.  The numeric stage M
*		  is employed only if NumericTol < SubdivTol.                M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *: Pllns on the bi-tangents of the two multivariates. M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbTangentToCrvAtTwoPts, MvarMVBiTangents2, MvarMVTriTangents           M
*   MvarMVTriTangentLine						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVBiTangents, bi-tangent                                             M
*****************************************************************************/
MvarPolylineStruct *MvarMVBiTangents(const MvarMVStruct *CMV1,
				     const MvarMVStruct *CMV2,
				     CagdRType SubdivTol,
				     CagdRType NumericTol)
{
    int i;
    MvarMVStruct *MV1, *MV2, *DMV, *MVTmp1,
        *DuMV1, *DvMV1, *DuMV2, *DvMV2, *MVs[3];
    MvarPolylineStruct *Solution;

    if (CMV1 == CMV2) {					  /* No self test. */
        MVAR_FATAL_ERROR(MVAR_ERR_TWO_SAME_MVS);
        return NULL;
    }

    if (MVAR_NUM_OF_MV_COORD(CMV1) != 3 || MVAR_NUM_OF_MV_COORD(CMV2) != 3) {
	MVAR_FATAL_ERROR(MVAR_ERR_PT_OR_LEN_MISMATCH);
	return NULL;
    }

    if (CMV1 -> GType != CMV2 -> GType) {
	MVAR_FATAL_ERROR(MVAR_ERR_SAME_GTYPE_EXPECTED);
	return NULL;
    }

    /* Bring both surfaces into a four-variate form. */
    if (CMV1 -> Dim == 2 && CMV2 -> Dim == 2) {
	MV1 = MvarPromoteMVToMV2(CMV1, 4, 0);   /* Four variate at axes 0,1. */
	MV2 = MvarPromoteMVToMV2(CMV2, 4, 2);   /* Four variate at axes 2,3. */

	/* Make sure domain are the same. */
	if (MV1 -> GType == MVAR_BSPLINE_TYPE) {
	    int i;
	    CagdRType Min, Max;

	    for (i = 0; i < 2; i++) {
		MvarMVDomain(MV1, &Min, &Max, i);
		BspKnotAffineTrans2(MV2 -> KnotVectors[i],
				    MV2 -> Lengths[i] + MV2 -> Orders[i],
				    Min, Max);
	    }
	    for (i = 2; i < 4; i++) {
		MvarMVDomain(MV2, &Min, &Max, i);
		BspKnotAffineTrans2(MV1 -> KnotVectors[i],
				    MV1 -> Lengths[i] + MV1 -> Orders[i],
				    Min, Max);
	    }
	}
    }
    else {
	MVAR_FATAL_ERROR(MVAR_ERR_GEOM_NO_SUPPORT);
	return NULL;
    }

    /* Compute the partial derivatives of the surfaces. */
    DuMV1 = MvarMVDerive(MV1, 0);
    DvMV1 = MvarMVDerive(MV1, 1);
    DuMV2 = MvarMVDerive(MV2, 2);
    DvMV2 = MvarMVDerive(MV2, 3);

    MVTmp1 = MvarMVCrossProd(DuMV1, DvMV1);

    MVs[0] = MvarMVDotProd(MVTmp1, DuMV2);
    MVs[1] = MvarMVDotProd(MVTmp1, DvMV2);
    DMV = MvarMVSub(MV1, MV2);
    MVs[2] = MvarMVDotProd(MVTmp1, DMV);

    MvarMVFree(MVTmp1);
    MvarMVFree(DuMV1);
    MvarMVFree(DvMV1);
    MvarMVFree(DuMV2);
    MvarMVFree(DvMV2);
    MvarMVFree(MV1);
    MvarMVFree(MV2);
    MvarMVFree(DMV);

    /* If this is a self bi-tangent test, make sure diagonal results are     */
    /* purged away by adding a fifth positive constraint for DMV distance.   */
    Solution = MvarMVsZeros1D(MVs, NULL, 3, SubdivTol / 10.0,
			      SubdivTol, NumericTol);

    for (i = 0; i < 3; i++)
        MvarMVFree(MVs[i]);

    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes tri-tangents of freeform bivariate. In other words, compute the M
* tangent plane at three points to the surface(s).                           M
* Let,									     M
*									     M
*	DMV12 = MV1(u, v) - MV2(r, s)					     V
*	DMV13 = MV1(u, v) - MV3(x, y)					     V
*	DMV23 = MV2(r, s) - MV3(x, y)					     V
*									     M
* then, compute the simultaneous solution of the following six equations:    M
*									     M
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV12 > = 0,				                     V
*     du      dv					                     V
*									     M
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     du      dv					                     V
*									     M
*   d MV2   d MV2		                 			     V
* < ----- x -----, DMV23 > = 0,				                     V
*     dr      ds					                     V
*                                                                            *
*   d MV2   d MV2		                 			     V
* < ----- x -----, DMV12 > = 0,				                     V
*     dr      ds					                     V
*                                                                            *
*   d MV3   d MV3		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     dx      dy					                     V
*                                                                            *
*   d MV3   d MV3		                 			     V
* < ----- x -----, DMV23 > = 0,				                     V
*     dx      dy					                     V
*                                                                            *
* PARAMETERS:                                                                M
*   CMV1, CMV2, CMV3:  The 3 multivariates to compute the tri-tangents for.  M
*		   If MV2 == MV3 == NULL, the self tri-tangents of MV1 are   M
*		   computed.						     M
*   Orientation:   0 for no effect, -1 or +1 for a request to get opposite   M
*		   or similar normal orientation bi tangencies only.         M
*   SubdivTol:     Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   NumericTol:    Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:    Points on the bi-tangents of the two multivariates.   M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbTangentToCrvAtTwoPts, MvarMVBiTangents, MvarMVBiTangents2            M
*   MvarMVTriTangentLine						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVTriTangents, tri-tangent                                           M
*****************************************************************************/
MvarPtStruct *MvarMVTriTangents(const MvarMVStruct *CMV1,
				const MvarMVStruct *CMV2,
				const MvarMVStruct *CMV3,
				int Orientation,
				CagdRType SubdivTol,
				CagdRType NumericTol)
{
    int i,
	Cnst = 0;
    CagdBType
	SelfTriTangent = FALSE;
    IRIT_STATIC_DATA CagdVType
	Translate = { 0.0, 0.0, 0.0 };
    MvarMVStruct *MV1, *MV2, *MV3, *DMV12, *DMV13, *DMV23,
	*MVTmp1, *MVTmp2, *MVs[11], *MV1Nrml, *MV2Nrml, *MV3Nrml;
    MvarConstraintType Constraints[11];
    MvarPtStruct *Solution;

    if (CMV2 == NULL && CMV3 == NULL) {			/* Do the self test. */
	CMV2 = CMV3 = CMV1;
	SelfTriTangent = TRUE;
    }
    if (CMV2 == NULL || CMV3 == NULL) {
	MVAR_FATAL_ERROR(MVAR_ERR_ONE_OR_THREE_EXPECTED);
	return NULL;
    }

    if (MVAR_NUM_OF_MV_COORD(CMV1) != 3 &&
	MVAR_NUM_OF_MV_COORD(CMV2) != 3 &&
	MVAR_NUM_OF_MV_COORD(CMV3) != 3) {
	MVAR_FATAL_ERROR(MVAR_ERR_PT_OR_LEN_MISMATCH);
	return NULL;
    }

    if (CMV1 -> GType != CMV2 -> GType || CMV1 -> GType != CMV3 -> GType) {
	MVAR_FATAL_ERROR(MVAR_ERR_SAME_GTYPE_EXPECTED);
	return NULL;
    }

    /* Bring all surfaces into a six-variate form. */
    if (CMV1 -> Dim == 2 && CMV2 -> Dim == 2 && CMV3 -> Dim == 2) {
	MV1 = MvarPromoteMVToMV2(CMV1, 6, 0);    /* Six variate at axes 0,1. */
	MV2 = MvarPromoteMVToMV2(CMV2, 6, 2);    /* Six variate at axes 2,3. */
	MV3 = MvarPromoteMVToMV2(CMV3, 6, 4);    /* Six variate at axes 4,5. */

	/* Make sure domain are the same. */
	if (MV1 -> GType == MVAR_BSPLINE_TYPE) {
	    int i;
	    CagdRType Min, Max;

	    for (i = 0; i < 2; i++) {
		MvarMVDomain(MV1, &Min, &Max, i);
		BspKnotAffineTrans2(MV2 -> KnotVectors[i],
				    MV2 -> Lengths[i] + MV2 -> Orders[i],
				    Min, Max);
		BspKnotAffineTrans2(MV3 -> KnotVectors[i],
				    MV3 -> Lengths[i] + MV3 -> Orders[i],
				    Min, Max);
	    }
	    for (i = 2; i < 4; i++) {
		MvarMVDomain(MV2, &Min, &Max, i);
		BspKnotAffineTrans2(MV1 -> KnotVectors[i],
				    MV1 -> Lengths[i] + MV1 -> Orders[i],
				    Min, Max);
		BspKnotAffineTrans2(MV3 -> KnotVectors[i],
				    MV3 -> Lengths[i] + MV3 -> Orders[i],
				    Min, Max);
	    }
	    for (i = 4; i < 6; i++) {
		MvarMVDomain(MV3, &Min, &Max, i);
		BspKnotAffineTrans2(MV1 -> KnotVectors[i],
				    MV1 -> Lengths[i] + MV1 -> Orders[i],
				    Min, Max);
		BspKnotAffineTrans2(MV2 -> KnotVectors[i],
				    MV2 -> Lengths[i] + MV2 -> Orders[i],
				    Min, Max);
	    }
	}
    }
    else {
	MVAR_FATAL_ERROR(MVAR_ERR_GEOM_NO_SUPPORT);
	return NULL;
    }

    DMV12 = MvarMVSub(MV1, MV2);
    DMV13 = MvarMVSub(MV1, MV3);
    DMV23 = MvarMVSub(MV2, MV3);

    /* Compute the partial derivatives of the first surface. */
    MVTmp1 = MvarMVDerive(MV1, 0);
    MVTmp2 = MvarMVDerive(MV1, 1);
    MV1Nrml = MvarMVCrossProd(MVTmp1, MVTmp2);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);
    MVs[Cnst] = MvarMVDotProd(MV1Nrml, DMV12);
    Constraints[Cnst++] = MVAR_CNSTRNT_ZERO;
    MVs[Cnst] = MvarMVDotProd(MV1Nrml, DMV13);
    Constraints[Cnst++] = MVAR_CNSTRNT_ZERO;

    /* Compute the partial derivatives of the second surface. */
    MVTmp1 = MvarMVDerive(MV2, 2);
    MVTmp2 = MvarMVDerive(MV2, 3);
    MV2Nrml = MvarMVCrossProd(MVTmp1, MVTmp2);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);
    MVs[Cnst] = MvarMVDotProd(MV2Nrml, DMV12);
    Constraints[Cnst++] = MVAR_CNSTRNT_ZERO;
    MVs[Cnst] = MvarMVDotProd(MV2Nrml, DMV23);
    Constraints[Cnst++] = MVAR_CNSTRNT_ZERO;

    /* Compute the partial derivatives of the third surface. */
    MVTmp1 = MvarMVDerive(MV3, 4);
    MVTmp2 = MvarMVDerive(MV3, 5);
    MV3Nrml = MvarMVCrossProd(MVTmp1, MVTmp2);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);
    MVs[Cnst] = MvarMVDotProd(MV3Nrml, DMV13);
    Constraints[Cnst++] = MVAR_CNSTRNT_ZERO;
    MVs[Cnst] = MvarMVDotProd(MV3Nrml, DMV23);
    Constraints[Cnst++] = MVAR_CNSTRNT_ZERO;

    if (Orientation) {
	MVs[Cnst] = MvarMVDotProd(MV1Nrml, MV2Nrml);
	Constraints[Cnst++] = Orientation > 0 ? MVAR_CNSTRNT_POSITIVE
					      : MVAR_CNSTRNT_NEGATIVE;
	MVs[Cnst] = MvarMVDotProd(MV1Nrml, MV3Nrml);
	Constraints[Cnst++] = Orientation > 0 ? MVAR_CNSTRNT_POSITIVE
					      : MVAR_CNSTRNT_NEGATIVE;
    }

    MvarMVFree(MV1Nrml);
    MvarMVFree(MV2Nrml);
    MvarMVFree(MV3Nrml);

    /* If this is a self bi-tangent test, make sure diagonal results are     */
    /* purged away by adding a fifth positive constraint for DMV distance.   */
    if (SelfTriTangent) {
	MVs[Cnst    ] = MvarMVDotProd(DMV12, DMV12);
	MVs[Cnst + 1] = MvarMVDotProd(DMV13, DMV13);
	MVs[Cnst + 2] = MvarMVDotProd(DMV23, DMV23);
	Translate[0] = -IRIT_SQR(SubdivTol * SUBDIV_REL_DISTANCE_TOL);
	MvarMVTransform(MVs[Cnst    ], Translate, 1.0);
	MvarMVTransform(MVs[Cnst + 1], Translate, 1.0);
	MvarMVTransform(MVs[Cnst + 2], Translate, 1.0);
	Constraints[Cnst++] = MVAR_CNSTRNT_POSITIVE;
	Constraints[Cnst++] = MVAR_CNSTRNT_POSITIVE;
	Constraints[Cnst++] = MVAR_CNSTRNT_POSITIVE;

	Solution = MvarMVsZeros0D(MVs, Constraints, Cnst, SubdivTol, NumericTol);
    }
    else {
	Solution = MvarMVsZeros0D(MVs, Constraints, Cnst, SubdivTol, NumericTol);
    }

    MvarMVFree(MV1);
    MvarMVFree(MV2);
    MvarMVFree(MV3);

    MvarMVFree(DMV12);
    MvarMVFree(DMV13);
    MvarMVFree(DMV23);

    for (i = 0; i < Cnst; i++)
	MvarMVFree(MVs[i]);

    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes all circles of prescribed radius that are tangent to given two  M
* curves.								     M
*   Solves for circles' centers P(x, y), using the following four equations  M
* in four unknowns (t, r, x, y), and R is the desired circle radius:	     M
*	||C1(t) - P||^2 = R^2,						     V
*	||C2(r) - P||^2 = R^2,						     V
*	< C1(t) - P, C1'(t) > = 0,					     V
*	< C2(t) - P, C2'(t) > = 0.					     V
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1, Crv2: The two curves to find the circles that is tangent to both.  M
*   Radius:     Of all the circle(s) that is tangent to Crv1/2.		     M
*   Tol:	Tolerance of approximation.  Subdiv Tol of MV Zeros.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:    List of the 4-tuples as (t, r, x, y).		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbCircTanTo2Crvs, MvarCircTanTo3Crvs	                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCircTanTo2Crvs, bi-tangent                                           M
*****************************************************************************/
MvarPtStruct *MvarCircTanTo2Crvs(const CagdCrvStruct *Crv1,
				 const CagdCrvStruct *Crv2,
				 CagdRType Radius,
				 CagdRType Tol)
{
    int i;
    CagdRType DmnMin[4], DmnMax[4], Trans;
    CagdBBoxStruct BBox1, BBox2;
    CagdSrfStruct *CircCenterSrf;
    MvarMVStruct *MTmp, *MTmp1, *MTmp2, *MCrv1, *MCrv2, *MDCrv1, *MDCrv2,
        *MCircCenter, *MVs[4];
    MvarPtStruct *Solution, *Pt;

    /* Compute the region in the XY plane where the circle-center can be    */
    /* as the expanded bbox pf the curves by (twice for safety) the radius. */
    CagdCrvBBox(Crv1, &BBox1);
    CagdCrvBBox(Crv2, &BBox2);

    for (i = 0; i < 2; i++) {
        if (BBox1.Min[i] > BBox2.Min[i])
	    BBox1.Min[i] = BBox2.Min[i] - Radius * 2.0;
	else
	    BBox1.Min[i] -= Radius * 2.0;

	if (BBox1.Max[i] < BBox2.Max[i])
	    BBox1.Max[i] = BBox2.Max[i] + Radius * 2.0;
	else
	    BBox1.Max[i] += Radius * 2.0;	
    }

    CircCenterSrf = CagdPrimPlaneSrf(BBox1.Min[0], BBox1.Min[1],
				     BBox1.Max[0], BBox1.Max[1], 0.0);

    MTmp = MvarCrvToMV(Crv1);  
    MCrv1 = MvarPromoteMVToMV2(MTmp, 4, 0);  
    MvarMVFree(MTmp);

    MTmp = MvarCrvToMV(Crv2);  
    MCrv2 = MvarPromoteMVToMV2(MTmp, 4, 1);  
    MvarMVFree(MTmp);

    MTmp = MvarSrfToMV(CircCenterSrf);  
    MCircCenter = MvarPromoteMVToMV2(MTmp, 4, 2);  
    MvarMVFree(MTmp);

    /* Make sure domains are consistent. */
    CagdCrvDomain(Crv1, &DmnMin[0], &DmnMax[0]);
    CagdCrvDomain(Crv2, &DmnMin[1], &DmnMax[1]);
    CagdSrfDomain(CircCenterSrf, &DmnMin[2], &DmnMax[2],
				&DmnMin[3], &DmnMax[3]);
    MCrv1 = MvarMVSetAllDomains(MCrv1, DmnMin, DmnMax, TRUE);
    MCrv2 = MvarMVSetAllDomains(MCrv2, DmnMin, DmnMax, TRUE);
    MCircCenter = MvarMVSetAllDomains(MCircCenter, DmnMin, DmnMax, TRUE);

    /* Build the constraints. */
    MDCrv1 = MvarMVDerive(MCrv1, 0);
    MDCrv2 = MvarMVDerive(MCrv2, 1);
    MTmp1 = MvarMVSub(MCrv1, MCircCenter);
    MTmp2 = MvarMVSub(MCrv2, MCircCenter);

    Trans = -IRIT_SQR(Radius);
    MVs[0] = MvarMVDotProd(MTmp1, MTmp1);
    MvarMVTransform(MVs[0], &Trans, 1.0);
    MVs[1] = MvarMVDotProd(MTmp2, MTmp2);
    MvarMVTransform(MVs[1], &Trans, 1.0);

    MVs[2] = MvarMVDotProd(MDCrv1, MTmp1);
    MVs[3] = MvarMVDotProd(MDCrv2, MTmp2);

    MvarMVFree(MDCrv1);
    MvarMVFree(MDCrv2);
    MvarMVFree(MTmp1);
    MvarMVFree(MTmp2);
    MvarMVFree(MCrv1);
    MvarMVFree(MCrv2);
    MvarMVFree(MCircCenter);

    Solution = MvarMVsZeros0D(MVs, NULL, 4, Tol, MVAR_CIRC_TAN_2CRVS_NUMER_TOL);

    for (i = 0; i < 4; i++)
	MvarMVFree(MVs[i]);

    /* Map the Circle centers to Eucliean space. */
    for (Pt = Solution; Pt != NULL; Pt = Pt -> Pnext) {
        CagdRType
	    *R = CagdSrfEval(CircCenterSrf, Pt -> Pt[2], Pt -> Pt[3]);

	Pt -> Pt[2] = R[1];
	Pt -> Pt[3] = R[2];
    }

    CagdSrfFree(CircCenterSrf);

    return Solution;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes all circles that are tangent to given three curves.	     M
*   Solves for circles' centers P(x, y), using the following process:        M
* Solve, symbolically for P in the following 2x2 system by Cremmer rule:     M
*	||C1(u) - P||^2 = ||C2(v) - P||^2,				     V
*	||C1(u) - P||^2 = ||C3(w) - P||^2,				     V
*    and substitute P into the following 3 euqations and solve 3 equations   M
* in (u, v, w):								     M
*	< C1(u) - P, C1'(u) > = 0,					     V
*	< C2(v) - P, C2'(v) > = 0,					     V
*	< C3(w) - P, C2'(w) > = 0.					     V
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1, Crv2, Crv3: The two curves to find the circles that is tangent to  M
*                  both.                                                     M
*   SubdivTol:     Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   NumericTol:    Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*   OneSideOrientation:   TRUE to compute tri-tangencies on one side of the  M
*                  curves only.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:    List of (u, v, w) solution points.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCircTanTo2Crvs		                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCircTanTo3Crvs, tri-tangent                                          M
*****************************************************************************/
MvarPtStruct *MvarCircTanTo3Crvs(const CagdCrvStruct *Crv1,
				 const CagdCrvStruct *Crv2,
				 const CagdCrvStruct *Crv3,
				 CagdRType SubdivTol,
				 CagdRType NumericTol,
				 CagdBType OneSideOrientation)
{
    IRIT_STATIC_DATA MvarConstraintType
        CTypes[7] = {
            MVAR_CNSTRNT_ZERO,
	    MVAR_CNSTRNT_ZERO,
	    MVAR_CNSTRNT_ZERO,
	    MVAR_CNSTRNT_POSITIVE,
	    MVAR_CNSTRNT_POSITIVE,
	    MVAR_CNSTRNT_POSITIVE,
	    MVAR_CNSTRNT_POSITIVE
        };
    CagdRType TMin, TMax,
        PSolNorm = -SubdivTol;
    MvarMVStruct *MVCrv1, *MVCrv2, *MVCrv3, *MVTan1, *MVTan2, *MVTan3,
	*MVVec[MVAR_MAX_PT_SIZE], *MVA1Split[MVAR_MAX_PT_SIZE],
        *MVA2Split[MVAR_MAX_PT_SIZE], *MVTmp1, *MVTmp2, *MVb1, *MVb2,
        *MVA1, *MVA2, *MVPDenom, *MVPNumer;
    MvarPtStruct *MVPts;
    CagdCrvStruct *DCrv;
    MvarBBoxStruct MVBBox;

    MVTmp1 = MvarCrvToMV(Crv1);
    MVCrv1 = MvarPromoteMVToMV2(MVTmp1, 3, 0);
    MvarMVFree(MVTmp1);

    MVTmp1 = MvarCrvToMV(Crv2);
    MVCrv2 = MvarPromoteMVToMV2(MVTmp1, 3, 1);
    MvarMVFree(MVTmp1);

    MVTmp1 = MvarCrvToMV(Crv3);
    MVCrv3 = MvarPromoteMVToMV2(MVTmp1, 3, 2);
    MvarMVFree(MVTmp1);

    /* Convert tangent curves. */
    DCrv = CagdCrvDerive(Crv1);
    MVTmp1 = MvarCrvToMV(DCrv);
    MVTan1 = MvarPromoteMVToMV2(MVTmp1, 3, 0);
    MvarMVFree(MVTmp1);
    CagdCrvFree(DCrv);

    DCrv = CagdCrvDerive(Crv2);
    MVTmp1 = MvarCrvToMV(DCrv);
    MVTan2 = MvarPromoteMVToMV2(MVTmp1, 3, 1);
    MvarMVFree(MVTmp1);
    CagdCrvFree(DCrv);

    DCrv = CagdCrvDerive(Crv3);
    MVTmp1 = MvarCrvToMV(DCrv);
    MVTan3 = MvarPromoteMVToMV2(MVTmp1, 3, 2);
    MvarMVFree(MVTmp1);
    CagdCrvFree(DCrv);

    /* Update mutual domains to be the same. */
    CagdCrvDomain(Crv1, &TMin, &TMax);
    MvarMVSetDomain(MVCrv2, TMin, TMax, 0, TRUE);
    MvarMVSetDomain(MVCrv3, TMin, TMax, 0, TRUE);
    MvarMVSetDomain(MVTan2, TMin, TMax, 0, TRUE);
    MvarMVSetDomain(MVTan3, TMin, TMax, 0, TRUE);
      
    CagdCrvDomain(Crv2, &TMin, &TMax);
    MvarMVSetDomain(MVCrv1, TMin, TMax, 1, TRUE);
    MvarMVSetDomain(MVCrv3, TMin, TMax, 1, TRUE);
    MvarMVSetDomain(MVTan1, TMin, TMax, 1, TRUE);
    MvarMVSetDomain(MVTan3, TMin, TMax, 1, TRUE);
      
    CagdCrvDomain(Crv3, &TMin, &TMax);
    MvarMVSetDomain(MVCrv1, TMin, TMax, 2, TRUE);
    MvarMVSetDomain(MVCrv2, TMin, TMax, 2, TRUE);
    MvarMVSetDomain(MVTan1, TMin, TMax, 2, TRUE);
    MvarMVSetDomain(MVTan2, TMin, TMax, 2, TRUE);

    /* Formulate distance constraints so that distance between Crv1 and Crv2 */
    /* is similar and the distance between Crv1 and Crv3 is similar:         */
    /*									     */
    /* < P - C1(u), P - C1(u) > = < P - C2(v), P - C2(v) >,		     */
    /* < P - C1(u), P - C1(u) > = < P - C3(w), P - C3(w) >,  or		     */
    /*									     */
    /* 2(C2(v) - C1(u)) P = C2(v)^2 - C1(u)^2,               		     */
    /* 2(C3(w) - C1(u)) P = C3(w)^2 - C1(u)^2, and solve for P = P(u, v, w). */
    MVA1 = MvarMVSub(MVCrv2, MVCrv1);
    MVTmp1 = MvarMVAdd(MVCrv2, MVCrv1);
    MVb1 = MvarMVDotProd(MVA1, MVTmp1);
    MvarMVTransform(MVA1, NULL, 2.0);
    MvarMVFree(MVTmp1);
    
    MVA2 = MvarMVSub(MVCrv3, MVCrv1);
    MVTmp1 = MvarMVAdd(MVCrv3, MVCrv1);
    MVb2 = MvarMVDotProd(MVA2, MVTmp1); 
    MvarMVTransform(MVA2, NULL, 2.0);
    MvarMVFree(MVTmp1);

    IRIT_GEN_COPY(MVA1Split, MvarMVSplitScalar(MVA1),
		  sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
    IRIT_GEN_COPY(MVA2Split, MvarMVSplitScalar(MVA2),
		  sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
    MvarMVFree(MVA1);
    MvarMVFree(MVA2);

    /* Solve for P = P(u, v, w). */
    IRIT_ZAP_MEM(MVVec, sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
    if (MVA1Split[0] != NULL || MVA2Split[0] != NULL ) {        /* Rational. */
        MvarMVStruct *MVRat[MVAR_MAX_PT_SIZE], *MVA1X, *MVA1Y, *MVA2X, *MVA2Y;

	IRIT_ZAP_MEM(MVRat, sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
	MVRat[0] = MVA1Split[0];
	MVRat[1] = MVA1Split[1];
	MVA1X = MvarMVMergeScalar(MVRat);
	MVRat[1] = MVA1Split[2];
	MVA1Y = MvarMVMergeScalar(MVRat);

	MVRat[0] = MVA2Split[0];
	MVRat[1] = MVA2Split[1];
	MVA2X = MvarMVMergeScalar(MVRat);
	MVRat[1] = MVA2Split[2];
	MVA2Y = MvarMVMergeScalar(MVRat);
	
        MVPDenom = MvarMVDeterminant2(MVA1X, MVA1Y, MVA2X, MVA2Y);
	MVVec[1] = MvarMVDeterminant2(MVb1, MVA1Y,
				      MVb2, MVA2Y);
	MVVec[2] = MvarMVDeterminant2(MVA1X, MVb1,
				      MVA1X, MVb2);
	MvarMVFree(MVA1X);
	MvarMVFree(MVA1Y);
	MvarMVFree(MVA2X);
	MvarMVFree(MVA2Y);
    }
    else {
        MVPDenom = MvarMVDeterminant2(MVA1Split[1], MVA1Split[2],
				      MVA2Split[1], MVA2Split[2]);
	MVVec[1] = MvarMVDeterminant2(MVb1, MVA1Split[2],
				      MVb2, MVA2Split[2]);
	MVVec[2] = MvarMVDeterminant2(MVA1Split[1], MVb1,
				      MVA2Split[1], MVb2);
    }

    if (MVA1Split[0] != NULL)
	MvarMVFree(MVA1Split[0]);
    MvarMVFree(MVA1Split[1]);
    MvarMVFree(MVA1Split[2]);
    if (MVA1Split[3] != NULL)
	MvarMVFree(MVA1Split[3]);

    if (MVA2Split[0] != NULL)
	MvarMVFree(MVA2Split[0]);
    MvarMVFree(MVA2Split[1]);
    MvarMVFree(MVA2Split[2]);
    if (MVA2Split[3] != NULL)
	MvarMVFree(MVA2Split[3]);

    MvarMVFree(MVb1);
    MvarMVFree(MVb2);

    if (MVAR_IS_RATIONAL_MV(MVVec[1])) {
	MvarMVStruct *MVRat[MVAR_MAX_PT_SIZE];

	assert(MVAR_IS_RATIONAL_MV(MVVec[2]));
	IRIT_ZAP_MEM(MVRat, sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);

	/* Decompose and re-assemble. */
	IRIT_GEN_COPY(MVA1Split, MvarMVSplitScalar(MVVec[1]),
		      sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
	IRIT_GEN_COPY(MVA2Split, MvarMVSplitScalar(MVVec[2]),
		      sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
	if (MvarMVsSame(MVA1Split[0], MVA2Split[0], NumericTol)) {
	    MVRat[0] = MVA1Split[0];
	    MVRat[1] = MVA1Split[1];
	    MVRat[2] = MVA2Split[1];
	    MVPNumer = MvarMVMergeScalar(MVRat);
	}
	else {
	    MVRat[0] = MvarMVMult(MVA1Split[0], MVA2Split[0]);
	    MVRat[1] = MvarMVMult(MVA1Split[1], MVA2Split[0]);
	    MVRat[2] = MvarMVMult(MVA1Split[0], MVA2Split[1]);
	    MVPNumer = MvarMVMergeScalar(MVRat);
	    MvarMVFree(MVRat[0]);
	    MvarMVFree(MVRat[1]);
	    MvarMVFree(MVRat[2]);
	}
	MvarMVFree(MVA1Split[0]);
	MvarMVFree(MVA1Split[1]);
	MvarMVFree(MVA2Split[0]);
	MvarMVFree(MVA2Split[1]);
    }
    else {
	MVPNumer = MvarMVMergeScalar(MVVec);
    }

#ifdef MVAR_DEBUG_DUMP_P_SOL
    /* Dump solution for P as a trivariate MV. A mem leak... */
    MVVec[0] = MVPDenom;
    MvarDbg(MvarMVMergeScalar(MVVec));
    MVVec[0] = NULL;
#endif /* MVAR_DEBUG_DUMP_P_SOL */

    MvarMVFree(MVVec[1]);
    MvarMVFree(MVVec[2]);

    /* Aim at positive weights as much as possible. */
    MvarMVBBox(MVPDenom, &MVBBox);
    if (IRIT_ABS(MVBBox.Min[0]) > IRIT_ABS(MVBBox.Max[0])) {
	CagdPType
	    Trans = { 0.0, 0.0, 0.0 };

	MvarMVTransform(MVPDenom, Trans, -1.0);
	MvarMVTransform(MVPNumer, Trans, -1.0);
    }

    /* Now formulate out the three following constraints with P's solution. */
    /*	< C1(u) - P, C1'(u) > = 0,					    */
    /*	< C2(v) - P, C2'(v) > = 0,					    */
    /*	< C3(w) - P, C3'(w) > = 0.					    */
    /* Note P can have poles, having both MVPDenom and MVPNumer vanish.     */
    /*   If orientation is to be observed, also add the three orientation   */
    /* positivity constraints for the three curves as:			    */
    /*  (C1(u) - P) x C1'(u) > 0,					    */
    /*  (C2(v) - P) x C2'(v) > 0,					    */
    /*  (C3(w) - P) x C3'(w) > 0.					    */

    MVTmp1 = MvarMVMultScalar(MVCrv1, MVPDenom);
    MVTmp2 = MvarMVSub(MVPNumer, MVTmp1);
    MVVec[0] = MvarMVDotProd(MVTmp2, MVTan1);
    if (OneSideOrientation)
        MVVec[4] = MvarMVCrossProdZ(MVTan1, MVTmp2);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);

    MVTmp1 = MvarMVMultScalar(MVCrv2, MVPDenom);
    MVTmp2 = MvarMVSub(MVPNumer, MVTmp1);
    MVVec[1] = MvarMVDotProd(MVTmp2, MVTan2);
    if (OneSideOrientation)
        MVVec[5] = MvarMVCrossProdZ(MVTan2, MVTmp2);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);

    MVTmp1 = MvarMVMultScalar(MVCrv3, MVPDenom);
    MVTmp2 = MvarMVSub(MVPNumer, MVTmp1);
    MVVec[2] = MvarMVDotProd(MVTmp2, MVTan3);
    if (OneSideOrientation)
        MVVec[6] = MvarMVCrossProdZ(MVTan3, MVTmp2);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);

    MvarMVFree(MVCrv1);
    MvarMVFree(MVCrv2);
    MvarMVFree(MVCrv3);
    MvarMVFree(MVTan1);
    MvarMVFree(MVTan2);
    MvarMVFree(MVTan3);

    /*   As P can have poles and hence introduces points at infinity (which */
    /* renders the solver slow), we add a heuristic requirement on the size */
    /* of the components of P to be at least subdivision level.		    */
    /*   Compute "PNumer^2 + PDenom^2 > PsolNorm":                          */
    MVTmp1 = MvarMVDotProd(MVPNumer, MVPNumer);
    MVTmp2 = MvarMVMult(MVPDenom, MVPDenom);
    MVVec[3] = MvarMVAdd(MVTmp1, MVTmp2);
    MvarMVTransform(MVVec[3], &PSolNorm, 1.0);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);

    /* Invoke the zero set solver. */    
    MVPts = MvarMVsZeros0D(MVVec, CTypes, OneSideOrientation ? 7 : 4,
			   SubdivTol, NumericTol);
    MvarMVFree(MVVec[3]);
    if (OneSideOrientation) {
        MvarMVFree(MVVec[4]);
        MvarMVFree(MVVec[5]);
        MvarMVFree(MVVec[6]);
    }

#ifdef DEBUG_DUMP_SOLUTIONS
    fprintf(stderr, "Num of sols = %d\n", CagdListLength(MVPts));
    {
        CagdRType *R, Radius, W;
	CagdPType Pt1, Pt2, Pt3;
	CagdPtStruct Cntr;
	MvarPtStruct *MVPt;

	for (MVPt = MVPts; MVPt != NULL; MVPt = MVPt -> Pnext) {
	    /* Update the position on three original primitive curves. */
	    R = CagdCrvEval(Crv1, MVPt -> Pt[0]);
	    CagdCoerceToE3(Pt1, &R, -1, Crv1 -> PType);

	    R = CagdCrvEval(Crv2, MVPt -> Pt[1]);
	    CagdCoerceToE3(Pt2, &R, -1, Crv2 -> PType);

	    R = CagdCrvEval(Crv3, MVPt -> Pt[2]);
	    CagdCoerceToE3(Pt3, &R, -1, Crv3 -> PType);

	    if (GMCircleFrom3Points(Cntr.Pt, &Radius, Pt1, Pt2, Pt3) || 1) {
		IPObjectStruct
		    *PCrv = IPGenCRVObject(BspCrvCreateCircle(&Cntr, Radius));

		AttrSetObjectRGBColor(PCrv, (int) IritRandom(100, 255),
				            (int) IritRandom(100, 255),
					    (int) IritRandom(100, 255));
		IPStderrObject(PCrv);
		IPFreeObject(PCrv);
	    }

	    R = MvarMVEval(MVPDenom, MVPt -> Pt);
	    W = R[1];
	    R = MvarMVEval(MVPNumer, MVPt -> Pt);
	    if (W == 0.0)
		fprintf(stderr, "Center at = %f  %f (%f)\n",
			R[1], R[2], W);
	    else
		fprintf(stderr, "Center at = %f  %f (%f)\n",
			R[1] / W, R[2] / W, W);
	}
    }
#endif /* DEBUG_DUMP_SOLUTIONS */

    MvarMVFree(MVVec[0]);
    MvarMVFree(MVVec[1]);
    MvarMVFree(MVVec[2]);

    MvarMVFree(MVPNumer);
    MvarMVFree(MVPDenom);

    return MVPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the constraints to solve for the tri-tangent line of three      M
* freeform bivariate. In other words, to compute the tangent line at three   M
* different points to the three surface(s).				     M
*   Let,								     M
*									     M
*	DMV12 = MV1(u, v) - MV2(s, t)					     V
*	DMV13 = MV1(u, v) - MV3(a, b)					     V
*	DMV23 = MV2(s, t) - MV3(a, b)					     V
*									     M
*   Then, compute the simultaneous solution of the following five equations: M
*									     M
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     du      dv					                     V
*									     *
*   d MV2   d MV2		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     dr      ds					                     V
*                                                                            *
*   d MV3   d MV3		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     da      db					                     V
*                                                                            *
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV12 > = 0,				                     V
*     du      dv					                     V
*                                                                            *
*   d MV3   d MV3		                 			     V
* < ----- x -----, DMV23 > = 0,				                     V
*     dx      dy					                     V
*                                                                            *
* PARAMETERS:                                                                M
*   CMV1, CMV2, CMV3: The 3 bivariates to compute the tri-tangent lines for. M
*   MVs:           To be populated with the computed MVS constraints.        M
*   Constraints:   To be populated with the constraints' types.              M
*                                                                            *
* RETURN VALUE:                                                              M
*   void								     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbTangentToCrvAtTwoPts, MvarMVBiTangents, MvarMVBiTangents2            M
*   MvarMVTriTangentLine						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVTriTangentLineCreateMVs, tri-tangent                               M
*****************************************************************************/
void MvarMVTriTangentLineCreateMVs(const MvarMVStruct *CMV1,
				   const MvarMVStruct *CMV2,
				   const MvarMVStruct *CMV3,
				   MvarMVStruct **MVs,
				   MvarConstraintType *Constraints)
{
    int i;
    CagdRType Min, Max;
    MvarMVStruct *MVSubUVST, *MVSubUVAB, *MVSubABST,
        *dUdMV1, *dVdMV1, *dSdMV2, *dTdMV2, *dAdMV3, *dBdMV3,
        *NormMV1, *NormMV2, *NormMV3,
	*MV1 = MvarPromoteMVToMV2(CMV1, 6, 0),
	*MV2 = MvarPromoteMVToMV2(CMV2, 6, 2),
	*MV3 = MvarPromoteMVToMV2(CMV3, 6, 4);

    /* Transform MV's domains to be common. */
    for (i = 0; i < 2; i++) {
	MvarMVDomain(MV1, &Min, &Max, i);
	BspKnotAffineTrans2(MV2 -> KnotVectors[i],
			    MV2 -> Lengths[i] + MV2 -> Orders[i],
			    Min, Max);
	BspKnotAffineTrans2(MV3 -> KnotVectors[i],
			    MV3 -> Lengths[i] + MV3 -> Orders[i],
			    Min, Max);
    }
    for (i = 2; i < 4; i++) {
	MvarMVDomain(MV2, &Min, &Max, i);
	BspKnotAffineTrans2(MV1 -> KnotVectors[i],
			    MV1 -> Lengths[i] + MV1 -> Orders[i],
			    Min, Max);
	BspKnotAffineTrans2(MV3 -> KnotVectors[i],
			    MV3 -> Lengths[i] + MV3 -> Orders[i],
			    Min, Max);
    }
    for (i = 4; i < 6; i++) {
	MvarMVDomain(MV3, &Min, &Max, i);
	BspKnotAffineTrans2(MV1 -> KnotVectors[i],
			    MV1 -> Lengths[i] + MV1 -> Orders[i],
			    Min, Max);
	BspKnotAffineTrans2(MV2 -> KnotVectors[i],
			    MV2 -> Lengths[i] + MV2 -> Orders[i],
			    Min, Max);
    }

    /* Compute some auxiliary MV's. */
    dUdMV1 = MvarMVDerive(MV1, 0);
    dVdMV1 = MvarMVDerive(MV1, 1);
    dSdMV2 = MvarMVDerive(MV2, 2);
    dTdMV2 = MvarMVDerive(MV2, 3);
    dAdMV3 = MvarMVDerive(MV3, 4);
    dBdMV3 = MvarMVDerive(MV3, 5);
    NormMV1 = MvarMVCrossProd(dUdMV1, dVdMV1);
    NormMV2 = MvarMVCrossProd(dSdMV2, dTdMV2);
    NormMV3 = MvarMVCrossProd(dAdMV3, dBdMV3);
    MVSubUVST = MvarMVSub(MV1, MV2);
    MVSubUVAB = MvarMVSub(MV1, MV3);
    MVSubABST = MvarMVSub(MV3, MV2);

    /* Build the constraints. */
    MVs[0] = MvarMVDotProd(MVSubUVST, NormMV1);
    MVs[1] = MvarMVDotProd(MVSubUVST, NormMV2);
    MVs[2] = MvarMVDotProd(MVSubUVAB, NormMV1);
    MVs[3] = MvarMVDotProd(MVSubUVAB, NormMV2);
    MVs[4] = MvarMVDotProd(MVSubABST, NormMV3);

    /* Add additional subdivision-only stage constraints. */
    MVs[5] = MvarMVDotProd(MVSubUVST, NormMV3);
    MVs[6] = MvarMVDotProd(MVSubUVAB, NormMV3);
    MVs[7] = MvarMVDotProd(MVSubABST, NormMV1);
    MVs[8] = MvarMVDotProd(MVSubABST, NormMV2);

    for (i = 0; i < 10; i++)
	Constraints[i] = i < 5 ? MVAR_CNSTRNT_ZERO
	                       : MVAR_CNSTRNT_ZERO_SUBDIV;

    MvarMVFree(MV1);
    MvarMVFree(MV2);
    MvarMVFree(MV3);
    MvarMVFree(dUdMV1);
    MvarMVFree(dVdMV1);
    MvarMVFree(dSdMV2);
    MvarMVFree(dTdMV2);
    MvarMVFree(dAdMV3);
    MvarMVFree(dBdMV3);
    MvarMVFree(NormMV1);
    MvarMVFree(NormMV2);
    MvarMVFree(NormMV3);
    MvarMVFree(MVSubUVST);
    MvarMVFree(MVSubUVAB);
    MvarMVFree(MVSubABST);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the constraints to solve for the tri-tangent line of three      M
* freeform bivariate. In other words, to compute the tangent line at three   M
* different points to the three surface(s).				     M
*   Let,								     M
*									     M
*	DMV12 = MV1(u, v) - MV2(s, t)					     V
*	DMV13 = MV1(u, v) - MV3(a, b)					     V
*	DMV23 = MV2(s, t) - MV3(a, b)					     V
*									     M
*   Then, compute the simultaneous solution of the following five equations: M
*									     M
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     du      dv					                     V
*									     *
*   d MV2   d MV2		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     dr      ds					                     V
*                                                                            *
*   d MV3   d MV3		                 			     V
* < ----- x -----, DMV13 > = 0,				                     V
*     da      db					                     V
*                                                                            *
*   d MV1   d MV1		                 			     V
* < ----- x -----, DMV12 > = 0,				                     V
*     du      dv					                     V
*                                                                            *
*   d MV3   d MV3		                 			     V
* < ----- x -----, DMV23 > = 0,				                     V
*     da      db					                     V
*                                                                            *
* PARAMETERS:                                                                M
*   CMV1, CMV2, CMV3: The 3 bivariates to compute the tri-tangent lines for. M
*   ETs:           To be populated with the computed ETs constraints.        M
*   Constraints:   To be populated with the constraints' types.              M
*                                                                            *
* RETURN VALUE:                                                              M
*   void								     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbTangentToCrvAtTwoPts, MvarMVBiTangents, MvarMVBiTangents2            M
*   MvarMVTriTangentLine, MvarMVTriTangentLineCreateMVs			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVTriTangentLineCreateETs, tri-tangent                               M
*****************************************************************************/
void MvarMVTriTangentLineCreateETs(const MvarMVStruct *CMV1,
				   const MvarMVStruct *CMV2,
				   const MvarMVStruct *CMV3,
				   MvarExprTreeStruct **ETs,
				   MvarConstraintType *Constraints)
{
    int i;
    CagdRType Min, Max;
    MvarMVStruct *dUdMV1, *dVdMV1, *dSdMV2, *dTdMV2, *dAdMV3, *dBdMV3,
	*MV1 = MvarPromoteMVToMV2(CMV1, 6, 0),
	*MV2 = MvarPromoteMVToMV2(CMV2, 6, 2),
	*MV3 = MvarPromoteMVToMV2(CMV3, 6, 4);
    MvarExprTreeStruct
        *NormET1, *NormET2, *NormET3, *MVSubUVST, *MVSubUVAB, *MVSubABST;

    for (i = 0; i < 2; i++) {
	MvarMVDomain(MV1, &Min, &Max, i);
	BspKnotAffineTrans2(MV2 -> KnotVectors[i],
			    MV2 -> Lengths[i] + MV2 -> Orders[i],
			    Min, Max);
	BspKnotAffineTrans2(MV3 -> KnotVectors[i],
			    MV3 -> Lengths[i] + MV3 -> Orders[i],
			    Min, Max);
    }
    for (i = 2; i < 4; i++) {
	MvarMVDomain(MV2, &Min, &Max, i);
	BspKnotAffineTrans2(MV1 -> KnotVectors[i],
			    MV1 -> Lengths[i] + MV1 -> Orders[i],
			    Min, Max);
	BspKnotAffineTrans2(MV3 -> KnotVectors[i],
			    MV3 -> Lengths[i] + MV3 -> Orders[i],
			    Min, Max);
    }
    for (i = 4; i < 6; i++) {
	MvarMVDomain(MV3, &Min, &Max, i);
	BspKnotAffineTrans2(MV1 -> KnotVectors[i],
			    MV1 -> Lengths[i] + MV1 -> Orders[i],
			    Min, Max);
	BspKnotAffineTrans2(MV2 -> KnotVectors[i],
			    MV2 -> Lengths[i] + MV2 -> Orders[i],
			    Min, Max);
    }

    /* Compute some auxiliary MV's. */
    dUdMV1 = MvarMVDerive(MV1, 0);
    dVdMV1 = MvarMVDerive(MV1, 1);
    dSdMV2 = MvarMVDerive(MV2, 2);
    dTdMV2 = MvarMVDerive(MV2, 3);
    dAdMV3 = MvarMVDerive(MV3, 4);
    dBdMV3 = MvarMVDerive(MV3, 5);

    NormET1 = MvarExprTreeSub(MvarExprTreeFromMV2(dUdMV1),
			      MvarExprTreeFromMV2(dVdMV1));
    NormET2 = MvarExprTreeSub(MvarExprTreeFromMV2(dSdMV2),
			      MvarExprTreeFromMV2(dTdMV2));
    NormET3 = MvarExprTreeSub(MvarExprTreeFromMV2(dAdMV3),
			      MvarExprTreeFromMV2(dBdMV3));
    MVSubUVST = MvarExprTreeSub(MvarExprTreeFromMV2(MV1),
				MvarExprTreeFromMV2(MV2));
    MVSubUVAB = MvarExprTreeSub(MvarExprTreeFromMV2(MV1),
				MvarExprTreeFromMV2(MV3));
    MVSubABST = MvarExprTreeSub(MvarExprTreeFromMV2(MV3),
				MvarExprTreeFromMV2(MV2));

    /* Build the constraints. */
    ETs[0] = MvarExprTreeDotProd(MVSubUVAB, NormET1);
    ETs[1] = MvarExprTreeDotProd(MVSubUVAB, NormET2);
    ETs[2] = MvarExprTreeDotProd(MVSubUVAB, NormET3);
    ETs[3] = MvarExprTreeDotProd(MVSubUVST, NormET1);
    ETs[4] = MvarExprTreeDotProd(MVSubUVST, NormET2);
    
    /* Add additional subdivision-only stage constraints. */
    ETs[5] = MvarExprTreeDotProd(MVSubUVST, NormET3);
    ETs[6] = MvarExprTreeDotProd(MVSubABST, NormET1);
    ETs[7] = MvarExprTreeDotProd(MVSubABST, NormET2);
    ETs[8] = MvarExprTreeDotProd(MVSubABST, NormET3);

    for (i = 0; i < 10; i++)
	Constraints[i] = i < 5 ? MVAR_CNSTRNT_ZERO
	                       : MVAR_CNSTRNT_ZERO_SUBDIV;

    MvarMVFree(MV1);
    MvarMVFree(MV2);
    MvarMVFree(MV3);
    MvarMVFree(dUdMV1);
    MvarMVFree(dVdMV1);
    MvarMVFree(dSdMV2);
    MvarMVFree(dTdMV2);
    MvarMVFree(dAdMV3);
    MvarMVFree(dBdMV3);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the bi-tangent line of two freeform curves.  If Crv1 == Crv2,   M
* the self bi-tangent is computed.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1, Crv2:    The 2 curves to compute their bi-tangent lines for.       M
*                  If Crv1 == Crv2, the self bi-tangent is computed.         M
*   SubdivTol:     Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   NumericTol:    Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   Pairs of parameter values on both curves, each pair    M
*                     defines a single bi-tangent.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbTangentToCrvAtTwoPts, MvarMVBiTangents, MvarMVBiTangents2            M
*   MvarMVTriTangentLineCreateMVs, MvarMVTriTangentLine			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVBiTangentLine, bi-tangent		                             M
*****************************************************************************/
MvarPtStruct *MvarMVBiTangentLine(const CagdCrvStruct *Crv1,
				  const CagdCrvStruct *Crv2,
				  CagdRType SubdivTol,
				  CagdRType NumericTol)
{

    IRIT_STATIC_DATA MvarConstraintType
        CTypes[3] = {
            MVAR_CNSTRNT_ZERO,
	    MVAR_CNSTRNT_ZERO,
	    MVAR_CNSTRNT_POSITIVE
        };
    CagdRType TMin1, TMax1, TMin2, TMax2;
    CagdCrvStruct *NewCrv1, *NewCrv2;
    MvarMVStruct *MVCrv1, *MVCrv2, *MVDCrv1, *MVDCrv2, *MVDiff,
        *MVTmp1, *MVTmp2, *MVTmp3, *MVs[3];
    MvarPtStruct *MVPt,
	*MVResult = NULL;

    /* Make sure the given curves are open end conditioned curve. */
    if (CAGD_IS_BEZIER_CRV(Crv1))
	NewCrv1 = CagdCnvrtBzr2BspCrv(Crv1);
    else if (CAGD_IS_BSPLINE_CRV(Crv1) && !BspCrvHasOpenEC(Crv1))
	NewCrv1 = BspCrvOpenEnd(Crv1);
    else
        NewCrv1 = CagdCrvCopy(Crv1);

    if (CAGD_IS_BEZIER_CRV(Crv2))
	NewCrv2 = CagdCnvrtBzr2BspCrv(Crv2);
    else if (CAGD_IS_BSPLINE_CRV(Crv2) && !BspCrvHasOpenEC(Crv2))
	NewCrv2 = BspCrvOpenEnd(Crv2);
    else
        NewCrv2 = CagdCrvCopy(Crv2);

    /* Make sure the domain is zero to one. */
    CagdCrvDomain(NewCrv1, &TMin1, &TMax1);
    BspKnotAffineTransOrder2(NewCrv1 -> KnotVector, NewCrv1 -> Order,
			     NewCrv1 -> Order + NewCrv1 -> Length,
			     0.0, 1.0);
    CagdCrvDomain(NewCrv2, &TMin2, &TMax2);
    BspKnotAffineTransOrder2(NewCrv2 -> KnotVector, NewCrv2 -> Order,
			     NewCrv2 -> Order + NewCrv2 -> Length,
			     0.0, 1.0);

    MVTmp1 = MvarCrvToMV(NewCrv1);
    MVCrv1 = MvarPromoteMVToMV2(MVTmp1, 2, 0);
    MvarMVFree(MVTmp1);
    CagdCrvFree(NewCrv1);

    MVTmp1 = MvarCrvToMV(NewCrv2);
    MVCrv2 = MvarPromoteMVToMV2(MVTmp1, 2, 1);
    MvarMVFree(MVTmp1);
    CagdCrvFree(NewCrv2);

    MVDCrv1 = MvarMVDerive(MVCrv1, 0);
    MVDCrv2 = MvarMVDerive(MVCrv2, 1);

    MVDiff = MvarMVSub(MVCrv1, MVCrv2);

    MVs[0] = MvarMVCrossProdZ(MVDiff, MVDCrv1);
    MVs[1] = MvarMVCrossProdZ(MVDiff, MVDCrv2);

    MvarMVFree(MVDCrv1);
    MvarMVFree(MVDCrv2);
    MvarMVFree(MVDiff);

    if (Crv1 == Crv2) {
        CagdRType
	    Trans = -SubdivTol * 5.0;

        /* Add a positivity constraint, purging the diagonal as a solution. */
        MVTmp1 = MvarBzrLinearInOneDir(2, 0, CAGD_PT_E1_TYPE);
        MVTmp2 = MvarBzrLinearInOneDir(2, 1, CAGD_PT_E1_TYPE);
	MVTmp3 = MvarMVSub(MVTmp1, MVTmp2);
	MVs[2] = MvarCnvrtBzr2BspMV(MVTmp3);
	MvarMVFree(MVTmp1);
	MvarMVFree(MVTmp2);
	MvarMVFree(MVTmp3);
	MvarMVTransform(MVs[2], &Trans, 1.0);

        MVResult = MvarMVsZeros0D(MVs, CTypes, 3, SubdivTol, NumericTol);

	MvarMVFree(MVs[2]);
    }
    else {
        MVResult = MvarMVsZeros0D(MVs, CTypes, 2, SubdivTol, NumericTol);
    }

    MvarMVFree(MVs[0]);
    MvarMVFree(MVs[1]);

    /* Map domain back to original curve's domain. */
    for (MVPt = MVResult; MVPt != NULL; MVPt = MVPt -> Pnext) {
        MVPt -> Pt[0] = MVPt -> Pt[0] * (TMax1 - TMin1) + TMin1;
        MVPt -> Pt[1] = MVPt -> Pt[1] * (TMax2 - TMin2) + TMin2;
    }

    return MVResult;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the tri-tangent line of three freeform bivariate. In other      M
* words, computes the tangent line at three different points to the three    M
* surface(s).  The result is a univariate (describing the line tri-tangent   M
* to the surfaces as it slides over the three surfaces).		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf1, Srf2, Srf3: The 3 bivariates to compute the tri-tangent lines for. M
*                  Assumed Bsplines surfaces with Open End Cond.             M
*   StepSize:      Tolerance of numeric tracing of univariate solution.      M
*   SubdivTol:     Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   NumericTol:    Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*   Euclidean:     True to return the result in Euclidean space, FALSE to    M
*                  return it in parameteric space.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:   Triplets of univariate solutions (three curves on the M
*                      three surfaces) of tri-tangent lines.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbTangentToCrvAtTwoPts, MvarMVBiTangents, MvarMVBiTangents2            M
*   MvarMVTriTangentLineCreateMVs					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVTriTangentLine, tri-tangent		                             M
*****************************************************************************/
CagdCrvStruct *MvarMVTriTangentLine(const CagdSrfStruct *Srf1,
				    const CagdSrfStruct *Srf2,
				    const CagdSrfStruct *Srf3,
				    CagdRType StepSize,
				    CagdRType SubdivTol,
				    CagdRType NumericTol,
				    int Euclidean)
{
    int i;
    CagdCrvStruct
        *AllCrvs = NULL;
    MvarMVStruct *MVs[10],
        *TmpMV1 = MvarSrfToMV(Srf1),
	*TmpMV2 = MvarSrfToMV(Srf2),
        *TmpMV3 = MvarSrfToMV(Srf3);
    MvarConstraintType Constraints[10];
    MvarPolylineStruct *MVResult, *MVIter;

    MvarMVTriTangentLineCreateMVs(TmpMV1, TmpMV2, TmpMV3, MVs, Constraints);
    MVResult = MvarMVsZeros1D(MVs, Constraints, 9,
			      StepSize, SubdivTol, NumericTol);

    MvarMVFree(TmpMV1);
    MvarMVFree(TmpMV2);
    MvarMVFree(TmpMV3);
    
    for (MVIter = MVResult; MVIter != NULL; MVIter = MVIter -> Pnext) {
	CagdRType *R;
	MvarPtStruct *MVPt,
	    *MVPts = MVIter -> Pl;
	IPVertexStruct *V1, *V2, *V3;
	IPPolygonStruct
	    *Pl1 = IPAllocPolygon(0, NULL, NULL),
	    *Pl2 = IPAllocPolygon(0, NULL, NULL),
	    *Pl3 = IPAllocPolygon(0, NULL, NULL);
	CagdCrvStruct *Crv1, *Crv2, *Crv3;

	for (MVPt = MVPts; MVPt != NULL; MVPt = MVPt -> Pnext) {
	    Pl1 -> PVertex = V1 = IPAllocVertex2(Pl1 -> PVertex);
	    Pl2 -> PVertex = V2 = IPAllocVertex2(Pl2 -> PVertex);
	    Pl3 -> PVertex = V3 = IPAllocVertex2(Pl3 -> PVertex);

	    R = CagdSrfEval(Srf1, MVPt -> Pt[0], MVPt -> Pt[1]);
	    CagdCoerceToE3(V1 -> Coord, &R, -1, Srf1 -> PType);
	    R = CagdSrfEval(Srf2, MVPt -> Pt[2], MVPt -> Pt[3]);
	    CagdCoerceToE3(V2 -> Coord, &R, -1, Srf2 -> PType);
	    R = CagdSrfEval(Srf3, MVPt -> Pt[4], MVPt -> Pt[5]);
	    CagdCoerceToE3(V3 -> Coord, &R, -1, Srf3 -> PType);
	}
	Pl1 -> PVertex = IPReverseVrtxList2(Pl1 -> PVertex);
	Pl2 -> PVertex = IPReverseVrtxList2(Pl2 -> PVertex);	
	Pl3 -> PVertex = IPReverseVrtxList2(Pl3 -> PVertex);

	Crv1 = UserPolyline2LinBsplineCrv(Pl1, 1);
	Crv2 = UserPolyline2LinBsplineCrv(Pl2, 1);
	Crv3 = UserPolyline2LinBsplineCrv(Pl3, 1);

	IPFreePolygon(Pl1);
	IPFreePolygon(Pl2);
	IPFreePolygon(Pl3);

	if (Crv1 != NULL && Crv2 != NULL && Crv3 != NULL) {
	    Crv1 -> Pnext = Crv2;
	    Crv2 -> Pnext = Crv3;
	    Crv3 -> Pnext = AllCrvs;
	    AllCrvs = Crv1;
	}
	else {
	    if (Crv1 != NULL)
		CagdCrvFree(Crv1);
	    if (Crv2 != NULL)
		CagdCrvFree(Crv2);
	    if (Crv3 != NULL)
		CagdCrvFree(Crv3);
	}
    }

    for (i = 0; i < 9; i++)
	MvarMVFree(MVs[i]);

    MvarPolylineFreeList(MVResult);

    return AllCrvs;
}
