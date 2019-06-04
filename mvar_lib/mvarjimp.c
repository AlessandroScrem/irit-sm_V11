/******************************************************************************
* mvarjimp.c - module to improve the jacobian of trivariate functions.	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber & Fady Massarwi, April 2013.			      *
******************************************************************************/

#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"
#include "inc_irit/mvar_lib.h"
#include "inc_irit/triv_iga.h"

#define MVAR_TV_JACOB_EPS	1e-10
#define DEBUG_PRINT_EXTREME_JACOBIAN_VALS

typedef enum {
    MVAR_EXTREME_POINT_INTERNAL = 0,
    MVAR_EXTREME_POINT_BOUNDARY
} MvarExtremePointAttrType;

static void MvarCalculateExtremePointsAux(
					const MvarMVStruct *MV,
					int PointDim,
					const CagdRType ParamsConstrValues[3],
					const CagdBType IsFixedDirection[3],
					MvarPtStruct **MVPtResultList);
static CagdRType MvarCalcJNumer(const MvarMVStruct *MV,
				const MvarPtStruct *ParamLoc);
static CagdRType MvarCalcDJDCoef(MvarMVStruct *MV, 
				 const MvarPtStruct *ParamLoc, 
				 int CtlPtIndex,
				 int Axis);
static int MvarPerformOneJImproveStep(TrivTVStruct *TV, IrtRType StepSize);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Auxiliary function of CalculateMVarExtremePoints.                        *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:			Input multivariate of interest.                      *
*   PointDim:		Dimension of the original multivariate		     *
                        (number of parameters).				     *
*   ParamsConstrValues: Fixed values of the three dimensions (u,v,w).	     *
*   IsFixedDirection:	Indicates at each axis if it is fixed or not.  	     *
*   MVPtResultList:	Accumulating extreme points list.   	             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarCalculateExtremePointsAux(
					const MvarMVStruct *MV,
					int PointDim,
					const CagdRType ParamsConstrValues[3],
					const CagdBType IsFixedDirection[3],
					MvarPtStruct **MVPtResultList)
{
    int i, DimAbsolouteIndex[3],
        NextParamIndex = 0,
        Dim = MV -> Dim;
    CagdRType BoundaryValues[3], DomainBoundaries[6];
    CagdBType LocalIsFixedDirection[3];
    MvarConstraintType Constraints[3];
    MvarPtStruct *TmpPt, *ExtremaPts;
    MvarMVStruct *DMV[3];

    assert(Dim <= 3);
    for (i = 0; i < Dim; ++i) {
	DMV[i] = MvarMVDerive(MV, i);
	Constraints[i] = MVAR_CNSTRNT_ZERO;
    }

    ExtremaPts = MvarMVsZeros0D(DMV, Constraints, Dim, 1e-4,
				MVAR_TV_JACOB_EPS);

    for (i = 0; i < Dim; ++i)
        MvarMVFree(DMV[i]);

    if (Dim == PointDim) {
        for (TmpPt = ExtremaPts; TmpPt != NULL; TmpPt = TmpPt -> Pnext) {
	    AttrSetIntAttrib(&TmpPt -> Attr, "extreme_type", 
			     MVAR_EXTREME_POINT_INTERNAL);
	}

	assert(*MVPtResultList == NULL);
	*MVPtResultList =
	        (MvarPtStruct *) CagdListAppend(*MVPtResultList, ExtremaPts);
    }
    else {
	MvarPtStruct
	    *ExtremaPtsIter = ExtremaPts;

	while (ExtremaPtsIter) {
	    MvarPtStruct
	        *NewPt = MvarPtNew(PointDim);
	    int i,
	        NextRelativeDim = 0;

	    for (i = 0; i < PointDim; ++i) {
		if (IsFixedDirection[i]) {
		    NewPt -> Pt[i] = ParamsConstrValues[i];
		}
		else {
		    NewPt -> Pt[i] = ExtremaPtsIter -> Pt[NextRelativeDim++];
		}
	    }

	    AttrSetIntAttrib(&NewPt -> Attr, "extreme_type", 
			     MVAR_EXTREME_POINT_BOUNDARY);

	    IRIT_LIST_PUSH(NewPt, *MVPtResultList);
	    ExtremaPtsIter = ExtremaPtsIter -> Pnext;
	}

	MvarPtFreeList(ExtremaPts);
    }

    for (i = 0; i < 3; ++i) {
	LocalIsFixedDirection[i] = IsFixedDirection[i];
	BoundaryValues[i] = ParamsConstrValues[i];
    }

    for (i = 0; i < Dim; ++i) {
	MvarMVDomain(MV, 
		     &DomainBoundaries[2 * i],
		     &DomainBoundaries[2 * i + 1], i);
    }

    for (i = 0; i < PointDim; ++i) {
	if (!IsFixedDirection[i]) {
	    DimAbsolouteIndex[NextParamIndex++] = i;
	}
    }

    switch (Dim) {
	case 3: {
	    TrivTVStruct
	        *JTV = MvarMVToTV(MV);
	    CagdSrfStruct 
	        **JBoundaries = TrivBndrySrfsFromTV(JTV);

	    TrivTVFree(JTV);

	    for (i = 0; i < 6; ++i) {
	        MvarMVStruct
		    *JBoundaryMVR = MvarSrfToMV(JBoundaries[i]);

		BoundaryValues[i / 2] = DomainBoundaries[i];
		LocalIsFixedDirection[i / 2] = TRUE;
		MvarCalculateExtremePointsAux(JBoundaryMVR, PointDim, 
					      BoundaryValues, 
					      LocalIsFixedDirection, 
					      MVPtResultList);
		CagdSrfFree(JBoundaries[i]);
		MvarMVFree(JBoundaryMVR);
		LocalIsFixedDirection[i / 2] = FALSE;
	    }
	    break;
	}
	case 2: {
	    CagdSrfStruct
	        *Srf = MvarMVToSrf(MV);
	    CagdCrvStruct
	        **BndryCrvs = CagdBndryCrvsFromSrf(Srf);

	    CagdSrfFree(Srf);

	    for (i = 0; i < 4; ++i) {
	        MvarMVStruct
		    *JBoundaryMVR = MvarCrvToMV(BndryCrvs[i]);

		BoundaryValues[DimAbsolouteIndex[i / 2]] = DomainBoundaries[i];
		LocalIsFixedDirection[DimAbsolouteIndex[i / 2]] = TRUE;
		MvarCalculateExtremePointsAux(JBoundaryMVR, PointDim, 
					      BoundaryValues, 
					      LocalIsFixedDirection, 
					      MVPtResultList);
		MvarMVFree(JBoundaryMVR);
		CagdCrvFree(BndryCrvs[i]);
		LocalIsFixedDirection[DimAbsolouteIndex[i / 2]] = FALSE;
	    }
	    break;
	}
	case 1: {
	    CagdRType TMin, TMax;
	    CagdCrvStruct
	        *Crv = MvarMVToCrv(MV);
	    MvarPtStruct *MinPt, *MaxPt;

	    CagdCrvDomain(Crv, &TMin, &TMax);

	    CagdCrvFree(Crv);

	    MinPt = MvarPtNew(PointDim);

	    for (i = 0; i < PointDim; ++i ) {
	        if (IsFixedDirection[i]) {
		    MinPt -> Pt[i] = ParamsConstrValues[i];
		}
	    }

	    MaxPt = MvarPtCopy(MinPt);

	    MinPt -> Pt[DimAbsolouteIndex[0]] = TMin;
	    MaxPt -> Pt[DimAbsolouteIndex[0]] = TMax;

	    AttrSetIntAttrib(&MaxPt -> Attr, "extreme_type", 
			     MVAR_EXTREME_POINT_BOUNDARY);
	    AttrSetIntAttrib(&MinPt -> Attr, "extreme_type", 
			     MVAR_EXTREME_POINT_BOUNDARY);

	    IRIT_LIST_PUSH(MinPt, *MVPtResultList);
	    IRIT_LIST_PUSH(MaxPt, *MVPtResultList);
	    break;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Removes duplicated points in a given multivariate points list, in place. M
*                                                                            *
* PARAMETERS:                                                                M
*   PtList: Input point list to make unique in place.                        M
*   Tol:    Equality tolerance on the different coefficient of the points.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void								     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMakeUniquePointsList, unique                                         M
*****************************************************************************/
void MvarMakeUniquePointsList(MvarPtStruct **PtList, CagdRType Tol)
{
    MvarPtStruct *Pt,
        *NewUniqueList = NULL;
    int IsFound, i;

    if (PtList == NULL || *PtList == NULL)
	return;
    
    while (*PtList != NULL) {
        /* if the current point is not found in NewUniqeList then add it. */
	MvarPtStruct 
	    *UniqueListIter = NewUniqueList;

	IRIT_LIST_POP(Pt, *PtList);
	IsFound = FALSE;
	while (UniqueListIter && !IsFound) {
	    for (i = 0; i < Pt -> Dim; ++i) {
		if (!IRIT_APX_EQ_EPS(UniqueListIter -> Pt[i], Pt -> Pt[i],
				     Tol))
		    break;
	    }

	    IsFound = i >= Pt -> Dim;
	    UniqueListIter = UniqueListIter -> Pnext;
	}

	if (IsFound) {
	    MvarPtFree(Pt);
	}
	else {
	    IRIT_LIST_PUSH(Pt, NewUniqueList);
	}
    }

    *PtList = CagdListReverse(NewUniqueList);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Calculates list of extreme points of a given scalar function.            M
*   Supports up to dimension 3 (i.e. trivariates).			     M
*   Note that the result, the returned list of parameteric locations, is a   M
* super set of the extreme values MV can assume as we decompose MV into low  M
* dimensional entities (i.e. boundary surfaces and curves) and an extrema in M
* a lower diemsnional entities does not mean it is an extreme in MV.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:		Input scalar function, represented as a multivariate, to     M
*               compute parametric locations of extreme values.              M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *: List of extreme multivariate points candidates	     M
*                   returns NULL in case of invalid input.   	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCalculateExtremePoints                                               M
*****************************************************************************/
MvarPtStruct *MvarCalculateExtremePoints(const MvarMVStruct *MV)
{
    if (!MV) {
	return NULL;
    }
    else {
	CagdBType IsFixedDirection[3];
        int i,
	    Dim = MV -> Dim;
	CagdRType ParamsConstrValues[3];
	MvarPtStruct
	    *PntList = NULL;

	for (i = 0; i < Dim; ++i) {
	    IsFixedDirection[i] = FALSE;
	}

	if (MVAR_NUM_OF_MV_COORD(MV) > 1) {
	    MvarMVStruct
		*MVTmp = MvarCoerceMVsTo(MV, MVAR_IS_RATIONAL_MV(MV) ?
							    MVAR_PT_P1_TYPE :
							    MVAR_PT_E1_TYPE);
							    
	    MvarCalculateExtremePointsAux(MVTmp, Dim, ParamsConstrValues, 
				          IsFixedDirection, &PntList);
	    MvarMVFree(MVTmp);
	}
	else
	    MvarCalculateExtremePointsAux(MV, Dim, ParamsConstrValues, 
				          IsFixedDirection, &PntList);

	MvarMakeUniquePointsList(&PntList, MVAR_TV_JACOB_EPS);

#	ifdef DEBUG_PRINT_EXTREME_PTS
	{
	    MvarPtStruct *TmpPt;

	    fprintf(stderr, "Extreme Pt list:\n");
	    for (TmpPt = PntList; TmpPt != NULL; TmpPt = TmpPt -> Pnext) {
		int i = AttrGetIntAttrib(TmpPt -> Attr, "extreme_type");

		fprintf(stderr, "Pt  (%d): %f %f %f\n", i,
			TmpPt -> Pt[0], TmpPt -> Pt[1], TmpPt -> Pt[2]);
	    }
	}
#	endif /* DEBUG_PRINT_EXTREME_PTS */

	return PntList;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Numerically evaluate the Jacobian of MV at a given parametric location.  *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:         The input multivariate.  Assumed non-rational.		     *
*   ParamLoc:   Parameteric location where to evaluate the Jacobian.         *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:                                                               *
*****************************************************************************/
static CagdRType MvarCalcJNumer(const MvarMVStruct *MV,
				const MvarPtStruct *ParamLoc)
{
    int i;
    IrtVecType Cross, Duvw[3];

    for (i = 0; i < 3; ++i) {
        int j,
	    i1 = i == 2 ? 0 : 2,
	    i2 = i == 0 ? 1 : 0;
	CagdRType *R;
	MvarMVStruct
	    *Srf = MvarMVFromMV(MV, ParamLoc -> Pt[i1], i1), /* UV, UV, VW. */
	    *Crv = MvarMVFromMV(Srf, ParamLoc -> Pt[i2], i2),   /* U, V, W. */
	    *DCrv = MvarMVDerive(Crv, 0);

	R = MvarMVEval(DCrv, &ParamLoc -> Pt[i]);
	for (j = 0; j < 3; ++j)	{
	    Duvw[i][j] = R[j + 1];
	}	

	MvarMVFree(Srf);
	MvarMVFree(Crv);
	MvarMVFree(DCrv);
    }

    GMVecCrossProd(Cross, Duvw[0], Duvw[1]);
    return GMVecDotProd(Cross, Duvw[2]);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Numerically estimate the derivative of the Jacobian with respect to      *
* one coefficient of one control point PointIndex in axis Axis.	             *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:         The input multivariate.  Assumed non-rational.		     *
*   ParamLoc:   Parameteric location where to evaluate the derivative.       *
*   CtlPtIndex: Index of control point to compute is Jacobian derivative.    *
*   Axis:       Axis in CtlPtIndex point to compute is Jacobian derivative.  *
*               1 for X, 2 for Y, 3 for Z.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:   The derivative, computed numerically.                       *
*****************************************************************************/
static CagdRType MvarCalcDJDCoef(MvarMVStruct *MV, 
				 const MvarPtStruct *ParamLoc, 
				 int CtlPtIndex,
				 int Axis)
{
    CagdRType J1, J2,
        OldCtlPtVal = MV -> Points[Axis][CtlPtIndex];

    assert(Axis <= CAGD_NUM_OF_PT_COORD(MV -> PType) &&
	   !MVAR_IS_RATIONAL_MV(MV));

    MV -> Points[Axis][CtlPtIndex] -= MVAR_TV_JACOB_EPS * 0.5;
    J1 = MvarCalcJNumer(MV, ParamLoc);
    MV -> Points[Axis][CtlPtIndex] += MVAR_TV_JACOB_EPS;
    J2 = MvarCalcJNumer(MV, ParamLoc);
    MV -> Points[Axis][CtlPtIndex] = OldCtlPtVal;

    return (J2 - J1) / MVAR_TV_JACOB_EPS;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Perform one numerical step of Jacobian improvement step.                 *
*                                                                            *
* PARAMETERS:                                                                *
*   TV:        To try and improve its Jacobian values, in place.             *
*   StepSize:  To use in numerical improvements.                             *
*   NumIters:  Number of numerical iterations to apply.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    TRUE if a step was effective, FALSE if failed.                   *
*****************************************************************************/
static int MvarPerformOneJImproveStep(TrivTVStruct *TV, IrtRType StepSize)
{
    CagdBType
        HasInteriorMax = FALSE;
    int a, b, c, k, Len, IsNotRational;
    MvarMVStruct 
	*MV = MvarTVToMV(TV),    
	*J = MvarCalculateTVJacobian(TV);
    MvarPtStruct *MinValUVW, *MaxValUVW,
        *JExtremePtsHead = MvarCalculateExtremePoints(J),
	*JExtremePts = JExtremePtsHead;
    CagdRType *PointValue,
        MinValS = IRIT_INFNTY,
        MaxValS = -IRIT_INFNTY;

#ifdef DEBUG_PRINT_EXTREME_PTS
    {
        MvarPtStruct *Pt;
	CagdRType Min[3], Max[3], Prm[3];

	MvarMVDomain(J, Min, Max, -1);

	for (Pt = JExtremePtsHead; Pt != NULL; Pt = Pt -> Pnext) {
	    int InternalExtr = AttrGetIntAttrib(Pt -> Attr, "extreme_type") 
					       == MVAR_EXTREME_POINT_INTERNAL;
	    IrtRType *R, V,
	        ExtrType = 0.0;

	    fprintf(stderr, "EPT = [%10.8lf  %10.8lf  %10.8lf] (%s)\n",
		    Pt -> Pt[0], Pt -> Pt[1], Pt -> Pt[2],
		    InternalExtr ? "Internal" : "Boundary");
	    R = MvarMVEval(J, Pt -> Pt);
	    V = R[1];

	    if (InternalExtr) {
	        for (a = 0; a < 100; a ++) { /* Examine local neighborhood. */
		    Prm[0] = IRIT_BOUND(Pt -> Pt[0] + IritRandom(-1e-5, 1e-5),
					Min[0], Max[0]);
		    Prm[1] = IRIT_BOUND(Pt -> Pt[1] + IritRandom(-1e-5, 1e-5),
					Min[1], Max[1]);
		    Prm[2] = IRIT_BOUND(Pt -> Pt[2] + IritRandom(-1e-5, 1e-5),
					Min[2], Max[2]);
		    R = MvarMVEval(J, Prm);
		    if (ExtrType == 0.0)
		        ExtrType = R[1] - V;
		    else
		        assert(ExtrType * (R[1] - V) > -IRIT_UEPS);
		}
	    }
	}
    }
#endif /* DEBUG_PRINT_EXTREME_PTS */

    if (JExtremePts == NULL) {
#       ifdef DEBUG_PRINT_EXTREME_JACOBIAN_VALS
            printf("No extreme Jacobian values detected\n");
#       endif /*  DEBUG_PRINT_EXTREME_JACOBIAN_VALS */
        MvarMVFree(MV);
	MvarMVFree(J);
        return FALSE;
    }

    /* Find internal minima and maxima among internal extremes candidates. */
    do {
	int AttrType = AttrGetIntAttrib(JExtremePts -> Attr, "extreme_type");

	/* Take into account only internal points. */
	if (AttrType != MVAR_EXTREME_POINT_INTERNAL) {
	    continue;
	}

	PointValue = MvarMVEval(J, JExtremePts -> Pt);
	if (MinValS > PointValue[1]) {
	    MinValS = PointValue[1];
	    MinValUVW = JExtremePts;
	}
	if (MaxValS < PointValue[1]) {
	    MaxValS = PointValue[1];
	    MaxValUVW = JExtremePts;
	}
    }
    while ((JExtremePts = JExtremePts -> Pnext) != NULL);
    HasInteriorMax = MinValS < MaxValS;

    MvarMVFree(J);

    if (MinValS == IRIT_INFNTY) {
#       ifdef DEBUG_PRINT_EXTREME_JACOBIAN_VALS
            printf("No internal extreme Jacobian values detected\n");
#       endif /*  DEBUG_PRINT_EXTREME_JACOBIAN_VALS */
        MvarMVFree(MV);
	MvarPtFreeList(JExtremePtsHead);
        return FALSE;
    }

#ifdef DEBUG_PRINT_EXTREME_JACOBIAN_VALS
    if (HasInteriorMax)
        printf("MinValue = %g (%g,%g,%g), MaxValue = %g (%g,%g,%g)\n",
	       MinValS,
	       MinValUVW -> Pt[0], MinValUVW -> Pt[1], MinValUVW -> Pt[2],
	       MaxValS,
	       MaxValUVW -> Pt[0], MaxValUVW -> Pt[1], MaxValUVW -> Pt[2]);
    else
        printf("MinValue = %g (%g,%g,%g), No Max Value\n",
	       MinValS,
	       MinValUVW -> Pt[0], MinValUVW -> Pt[1], MinValUVW -> Pt[2]);

#endif /* DEBUG_PRINT_EXTREME_JACOBIAN_VALS */

    Len = TRIV_TV_UPT_LST_LEN(TV) *
	  TRIV_TV_VPT_LST_LEN(TV) *
          TRIV_TV_WPT_LST_LEN(TV);
    IsNotRational = !CAGD_IS_RATIONAL_PT(TV -> PType),
    assert(CAGD_NUM_OF_PT_COORD(TV -> PType) == 3);

    /* Update only interior control points. */
    for (c = 1; c < TRIV_TV_WPT_LST_LEN(TV) - 1; c++) {
        for (b = 1; b < TRIV_TV_VPT_LST_LEN(TV) - 1; b++) {
	    for (a = 1; a < TRIV_TV_UPT_LST_LEN(TV) - 1; a++) {
	        int Idx = TRIV_MESH_UVW(TV, a, b, c);

	        /* Numerically compute derivative by adding some EPS to the */
	        /* Euclidean points at each direction.			    */
	        for (k = 1; k <= 3; k++) {
		    CagdRType 
		        DVal = MvarCalcDJDCoef(MV, MinValUVW, Idx, k) -
		               (HasInteriorMax ?
				   MvarCalcDJDCoef(MV, MaxValUVW, Idx, k) :
				   0.0);

		    TV -> Points[k][Idx] += StepSize * DVal;
		}
	    }
	}
    }

    MvarMVFree(MV);
    MvarPtFreeList(JExtremePtsHead);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Attempt the improve the Jacobian of the given trivariate by adjusting    M
* interior control-points of the TV.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   TV:        To try and improve (make more uniform) its Jacobian.          M
*   StepSize:  Numerical step size to move along the gradient (at the        M
*              Jacobian extreme values.)				     M
*   NumIters:  Number of numerical iterations to allow in the improvement    M
*              process.						             M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarCalculateTVJacobian			                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarTrivJacobianImprove                                                  M
*****************************************************************************/
void MvarTrivJacobianImprove(TrivTVStruct *TV,
			     CagdRType StepSize,
			     int NumIters)
{
    int i;

    for (i = 0; i < NumIters; i++) {
        if (!MvarPerformOneJImproveStep(TV, StepSize))
	    break;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Calculates the Jacobian of a given trivariate.		             M
*                                                                            *
* PARAMETERS:                                                                M
*   TV:   The input trivariate.                                              M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: The Jacobian in a multi variate representation.          M
*                                                                            *
* SEE ALSO:                                                                  M
*   Symb2DSrfJacobian, MvarTrivJacobianImprove, MvarCalculateExtremePoints   M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCalculateTVJacobian                                                  M
*****************************************************************************/
MvarMVStruct *MvarCalculateTVJacobian(const TrivTVStruct *TV)
{
    MvarMVStruct
        *MV = MvarTVToMV(TV),
	*DMVdU = MvarMVDerive(MV, 0),
	*DMVdV = MvarMVDerive(MV, 1),
	*DMVdW = MvarMVDerive(MV, 2),
        *DMVDuXDv = MvarMVCrossProd(DMVdU, DMVdV),
	*J = MvarMVDotProd(DMVDuXDv, DMVdW);

    MvarMVFree(MV);
    MvarMVFree(DMVdU);
    MvarMVFree(DMVdV);
    MvarMVFree(DMVdW);
    MvarMVFree(DMVDuXDv);
    
    return J;
}
