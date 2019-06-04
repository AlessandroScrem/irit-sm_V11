/******************************************************************************
* ZrAux.c - Auxiliary tools to compute zero sets of multivariates.	      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, July. 97.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/geom_lib.h"

typedef struct MvarMVParamDomainStruct {
    CagdRType Min, Max;					      /* The domain. */
} MvarMVParamDomainStruct;

IRIT_GLOBAL_DATA int
    _MVGlblZeroApplyDomainReduction = TRUE,
    _MVGlblZeroApplyGradPreconditioning = FALSE,
    _MVGlblZeroApplyParallelHyperPlaneTest = FALSE,
    _MVGlblZeroApplyNormalConeTest = TRUE;
IRIT_GLOBAL_DATA MvarMVsZerosSubdivCallBackFuncType
    _MVGlblZeroSubdivCallBackFunc = NULL;

#ifdef DEBUG_DUMP_DOMAINS
IRIT_STATIC_DATA FILE
    *GlblDumpDomainsFile = NULL;
#endif /* DEBUG_DUMP_DOMAINS */

static MvarMVStruct **MvarMVsOrthogonalizeGrads(MvarMVStruct **MVs,
						int NumOfMVs);
static void MVarMVHyperPlanesBound(const MvarMVStruct *MV,
				   CagdRType *Coeffs,
				   CagdRType *Coeff0Max,
				   CagdRType *Coeff0Min);

/* #define MVAR_DEBUG_DEPTH */
#ifdef MVAR_DEBUG_DEPTH
#define MVAR_DEBUG_MAX_DEPTH	200
IRIT_STATIC_DATA int
    MvarSubdivLevel[MVAR_DEBUG_MAX_DEPTH] = { 0 };
#endif /* MVAR_DEBUG_DEPTH */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the use (or not) of the normal cone tests inside the multivariate   M
* subdivisions' zero set solver.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   NormalConeTest:   New setting for normal cone testing usage.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting for normal cone testing usage.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarMVsZerosDomainReduction,				     M
*   MvarMVsZerosGradPreconditioning, MvarMVsZerosSetCallBackFunc	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosNormalConeTest                                               M
*****************************************************************************/
int MvarMVsZerosNormalConeTest(int NormalConeTest)
{
    int OldVal = _MVGlblZeroApplyNormalConeTest;

    _MVGlblZeroApplyNormalConeTest = NormalConeTest;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the use (or not) of gradient preconditiong option - application of  M
* and orthogonalization process over the gradients in multivariate	     M
* subdivisions' zero set solver.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   GradPreconditioning:   New setting for the gradient orthogonalization.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting of gradient orthogonalization.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarMVsZerosNormalConeTest, MvarMVsZerosDomainReduction,   M
*   MvarMVsZerosSetCallBackFunc						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosGradPreconditioning                                          M
*****************************************************************************/
int MvarMVsZerosGradPreconditioning(int GradPreconditioning)
{
    int OldVal = _MVGlblZeroApplyGradPreconditioning;

    _MVGlblZeroApplyGradPreconditioning = GradPreconditioning;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the use (or not) of the domain reduction option - Bezier (and       M
* (B-spline) clipping in multivariate subdivisions' zero set solver.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   DomainReduction:   New setting for the domain reduction option.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting for normal cone testing usage.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarMVsZerosNormalConeTest,				     M
*   MvarMVsZerosGradPreconditioning, MvarMVsZerosSetCallBackFunc             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosDomainReduction                                              M
*****************************************************************************/
int MvarMVsZerosDomainReduction(int DomainReduction)
{
    int OldVal = _MVGlblZeroApplyDomainReduction;

    _MVGlblZeroApplyDomainReduction = DomainReduction;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the use (or not) of the parallel plane termination criteria         M
* in multivariate subdivisions' zero set solver.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   ParallelHPlaneTest:   New setting for the domain reduction option.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting for normal cone testing usage.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarMVsZerosNormalConeTest,				     M
*   MvarMVsZerosGradPreconditioning, MvarMVsZerosDomainReduction,            M
*   MvarMVsZerosSetCallBackFunc						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosParallelHyperPlaneTest                                       M
*****************************************************************************/
int MvarMVsZerosParallelHyperPlaneTest(int ParallelHPlaneTest)
{
    int OldVal = _MVGlblZeroApplyParallelHyperPlaneTest;

    _MVGlblZeroApplyParallelHyperPlaneTest = ParallelHPlaneTest;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the use (or not) of a call back function that is invoked at every   M
* node of the subdivision tree process.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   SubdivCallBackFunc:   Call back function to use in the MV zeros'	     M
*			  subdivision stage.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVsZerosSubdivCallBackFuncType:       Old setting.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarMVsZerosDomainReduction,				     M
*   MvarMVsZerosGradPreconditioning, MvarMVsZerosNormalConeTest		     M
*   MvarMVsZerosParallelHyperPlaneTest		                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosSetCallBackFunc                                              M
*****************************************************************************/
MvarMVsZerosSubdivCallBackFuncType MvarMVsZerosSetCallBackFunc(
			MvarMVsZerosSubdivCallBackFuncType SubdivCallBackFunc)
{
    MvarMVsZerosSubdivCallBackFuncType
	OldVal = _MVGlblZeroSubdivCallBackFunc;

    _MVGlblZeroSubdivCallBackFunc = SubdivCallBackFunc;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Make sure all given MVs are in the same function space.                  M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   NumOfMVs:     Size of the MVs vector.		                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:    TRUE if in same function space, FALSE otherwise.           M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsSameSpace                                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosSameSpace                                                    M
*****************************************************************************/
CagdBType MvarMVsZerosSameSpace(MvarMVStruct **MVs, int NumOfMVs)
{
    int l;

    for (l = 1; l < NumOfMVs; l++) {
        if (!MvarMVsSameSpace(MVs[0], MVs[l], IRIT_EPS))
	    return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Extracts the subdivision parametric domain from a multivariate.          *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:         Multivariate to extract doamin from.			     *
*   Min, Max:   Subdivision parametric domain to extract.		     *
*   Dir:        Direction of sought domain.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
void MvarGetSubdivParamDomains(const MvarMVStruct *MV,
				      CagdRType *Min,
				      CagdRType *Max,
				      int Dir)
{
    if (MVAR_IS_BEZIER_MV(MV))
        MvarMVAuxDomainSlotGet(MV, Min, Max, Dir);
    else
	MvarMVDomain(MV, Min, Max, Dir);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Using the ideas of "Bezier-clipping" for both Bezier and/or B-spline,    *
* computes the possible reduced domain at both the minimum and the maximum   *
* of the domain, that is known to hold no zeros.			     *
*   Assumes all multivariates have the same domains and dimensions.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:        A vector of the multivariate constraints.                    *
*   NumOfMVs:   Size of the above multivariate constraints vector.           *
*   Dir:        Direction to compute domain reduction.			     *
*   SubdivTol:  Tolerance of the subdivision process.  Tolerance is          *
*	        measured in the parametric space of the multivariates.       *
*   TMin, TMax: Computed reduced domain.                                     *
*   SameSpace:  True if all MVs share the same function space.	             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:     TRUE if there might be some zeros, FALSE if no zero could exist *
*            in this case (i.e. all coefficients of the same sign, etc.)     *
*****************************************************************************/
int MvarMVsReduceMvsDomains(MvarMVStruct **MVs,
			    int NumOfMVs,
			    int Dir,
			    CagdRType SubdivTol,
			    CagdRType *TMin,
			    CagdRType *TMax,
			    CagdBType SameSpace)
{
    int i,
	Dim = MVs[0] -> Dim,
	*Indices = (int *) IritMalloc(sizeof(int) * Dim);
    CagdRType
        OrigTMin = *TMin,
        OrigTMax = *TMax;

    if (SameSpace && MVs[0] -> Dim == NumOfMVs) { 
        MvarMVStruct **NewMVs;

        /* Orthogonalize the (gradients of the) constraints at the center   */
	/* of the domain. 						    */
	if ((NewMVs = MvarMVsOrthogonalizeGrads(MVs, NumOfMVs)) != NULL) {
	    for (i = 0; i < NumOfMVs; i++) {
		MvarMVFree(MVs[i]);
		MVs[i] = NewMVs[i];
	    }
	}
    }

    for (i = 0; i < NumOfMVs; i++) {
	const MvarMVStruct
	    *MV = MVs[i];
	int NumOfSols, Index, j, k,
	    Length = MV -> Lengths[Dir],
	    Order = MV -> Orders[Dir],
	    Step = MVAR_NEXT_DIM(MV, Dir);
	CagdRType *Nodes, Sols[4], t;
	IPVertexStruct *V, *VEnv, *VNext,
	    *VUpperEnv = NULL,
	    *VLowerEnv = NULL;

	if (Length <= 1)
	    continue;

	if (MVAR_IS_BSPLINE_MV(MV)) {
	    Nodes = BspKnotNodes(MV -> KnotVectors[Dir],
				 Length + Order, Order);
	}
	else { /* MVAR_IS_BEZIER_MV(MV) */
	    Nodes = IritMalloc(sizeof(CagdRType) * Length);
	    for (j = 0; j < Length; j++)
	        Nodes[j] = j / ((CagdRType) (Length - 1));
	}

	IRIT_ZAP_MEM(Indices, sizeof(int) * Dim);

	/* Compute the convex hull of the control mesh.  For each node value */
	/* we need only take minimum and maximum so 2 coefficients per node. */
	for (k = j = 0; k < Length; k++) {
	    CagdRType
	        *Pts = MV -> Points[1],
	        Min = IRIT_INFNTY,
	        Max = -IRIT_INFNTY;

	    Indices[Dir] = k;
	    Index = Step * k; 
	    do {
	        if (Min > Pts[Index])
		    Min = Pts[Index];
	        if (Max < Pts[Index])
		    Max = Pts[Index];
	    }
	    while (MVAR_INC_SKIP_MESH_INDICES(MV, Indices, Dir, Index));

	    if (VUpperEnv != NULL &&
		IRIT_APX_EQ_EPS(VUpperEnv -> Coord[0], Nodes[k], IRIT_UEPS)) {
	        /* Duplicated node values can happen if Order identical     */
	        /* knots are presented in the interior of the doamin.       */
	        VUpperEnv -> Coord[1] = IRIT_MAX(VUpperEnv -> Coord[1], Max);
	        VLowerEnv -> Coord[1] = IRIT_MIN(VLowerEnv -> Coord[1], Min);
	    }
	    else {
	        VUpperEnv = IPAllocVertex2(VUpperEnv);
		VLowerEnv = IPAllocVertex2(VLowerEnv);
		VUpperEnv -> Coord[0] = VLowerEnv -> Coord[0] = Nodes[k];
		VUpperEnv -> Coord[1] = IRIT_FABS(Max) < IRIT_UEPS ? 0.0 : Max;
		VLowerEnv -> Coord[1] = IRIT_FABS(Min) < IRIT_UEPS ? 0.0 : Min;
		VUpperEnv -> Coord[2] = VLowerEnv -> Coord[2] = 0.0;
	    }
	}

	VUpperEnv = IPReverseVrtxList2(VUpperEnv);
	VLowerEnv = IPReverseVrtxList2(VLowerEnv);
	GMMonotonePolyConvex(VUpperEnv, TRUE);
	GMMonotonePolyConvex(VLowerEnv, FALSE);

	/* Chain the two envelopes into one circular list (convex hull). */
	VEnv = IPReverseVrtxList2(VUpperEnv);
	((IPVertexStruct *) IPGetLastVrtx(VEnv)) -> Pnext = VLowerEnv;
	((IPVertexStruct *) IPGetLastVrtx(VLowerEnv)) -> Pnext = VEnv;

#	ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugCHOutputPrint, FALSE) {
	        IPVertexDbg(VEnv);
	    }
	}
#	endif /* DEBUG */

	V = VEnv;
	NumOfSols = 0;
	do {
	    VNext = V -> Pnext;

	    if (!IRIT_PT_APX_EQ_E2_EPS(VNext -> Coord,
				       V -> Coord, IRIT_UEPS)) {
	        CagdRType
		    Sol = -IRIT_INFNTY;

		if (V -> Coord[1] == 0.0) {
		    Sol = V -> Coord[0];
		}
		else if (V -> Coord[1] * VNext -> Coord[1] < 0.0) {
		    t = VNext -> Coord[1] / (VNext -> Coord[1] -
					                    V -> Coord[1]);
		    Sol = V -> Coord[0] * t + VNext -> Coord[0] * (1.0 - t);
		}

		if (Sol > -IRIT_INFNTY) {
		    for (j = 0; j < NumOfSols; j++)
		        if (IRIT_APX_EQ_EPS(Sols[j],  Sol, IRIT_UEPS))
			    break;

		    if (j >= NumOfSols)
		        Sols[NumOfSols++] = Sol;
		}
	    }

	    V = VNext;
	}
	while (V != VEnv);

	IPFreeVertexList(VEnv);

	switch (NumOfSols) {
	    case 0:
	        /* CH is above or below the zero - no intersection. */
	        IritFree(Indices);
		IritFree(Nodes);

		*TMin = OrigTMin;
		*TMax = OrigTMax;

		return FALSE;
	    case 1:
	        /* Select a very narrow domain... */
	        IritFree(Indices);
		IritFree(Nodes);
	        *TMin = Sols[0] - IRIT_MAX(SubdivTol, 10.0 * IRIT_EPS);
		*TMax = Sols[0] + IRIT_MAX(SubdivTol, 10.0 * IRIT_EPS);
		*TMin = IRIT_MAX(OrigTMin, *TMin);
		*TMax = IRIT_MIN(OrigTMax, *TMax);

		return TRUE;
	    case 2:
	        if (Sols[0] > Sols[1]) {
		    IRIT_SWAP(CagdRType, Sols[0], Sols[1]);
		}

	        /* Make sure we miss no root due to numeric errors. */
		Sols[0] = IRIT_MAX(OrigTMin, Sols[0] - 10.0 * IRIT_UEPS);
		Sols[1] = IRIT_MIN(OrigTMax, Sols[1] + 10.0 * IRIT_UEPS);
		break;
	    default: /* NumOfSols > 2. */
	        /* Something is wrong in our solutions' computations. */
	        IritFree(Indices);
		IritFree(Nodes);
	        *TMin = OrigTMin;
		*TMax = OrigTMax;
#		ifdef DEBUG_NUM_OF_SOLS
		    MvarMVsReduceMvsDomains(MVs, NumOfMVs, Dir,
					    SubdivTol, TMin, TMax,
					    SameSpace);
		    MVAR_FATAL_ERROR(MVAR_ERR_INCONS_DOMAIN);
#		endif /* DEBUG_NUM_OF_SOLS */

		return TRUE;
	}

#	ifdef DEBUG
        {
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintReducedDomain, FALSE) {
	        int k;
	        MvarMVStruct *MVMin, *MVMax;
		CagdRType *R, *PtsMin, *PtsMax,
		    Min = Nodes == NULL ? 0.0 : Nodes[0],
		    Max = Nodes == NULL ? 1.0 : Nodes[Length - 1];

		if (Sols[0] > Min &&
		    !IRIT_APX_EQ_EPS(Sols[0], Min, IRIT_EPS)) {
		    MVMin = MvarMVRegionFromMV(MV, Min, Sols[0], Dir);
		    PtsMin = MVMin -> Points[1];

		    for (R = &PtsMin[1], k = MVAR_CTL_MESH_LENGTH(MVMin);
			 k-- > 1; ) {
		        if (*R++ * *PtsMin < -IRIT_EPS) {
			    printf("Error min in domain reduction computation\n");
			    MvarMVsReduceMvsDomains(MVs, NumOfMVs,
						    Dir, SubdivTol,
						    TMin, TMax,
						    SameSpace);
			}
		    }

		    MvarMVFree(MVMin);
		}

		if (Sols[1] < Max &&
		    !IRIT_APX_EQ_EPS(Sols[1], Max, IRIT_EPS)) {
		    MVMax = MvarMVRegionFromMV(MV, Sols[1], Max, Dir);
		    PtsMax = MVMax -> Points[1];
		
		    for (R = &PtsMax[1], k = MVAR_CTL_MESH_LENGTH(MVMax);
			 k-- > 1; ) {
		        if (*R++ * *PtsMax < -IRIT_EPS) {
			    printf("Error max in domain reduction computation\n");
			    MvarMVsReduceMvsDomains(MVs, NumOfMVs,
						    Dir, SubdivTol,
						    TMin, TMax,
						    SameSpace);
			}
		    }
		}
	    }
	}
#       endif /* DEBUG */

	IritFree(Nodes);

	/* The MV that gives the best (smallest) bound should be selected.   */
	*TMin = IRIT_MAX(*TMin, Sols[0]);
	*TMax = IRIT_MIN(*TMax, Sols[1]);
	if (*TMin > *TMax) {
	    /* The domains clipped from different constraints are disjoint.  */
	    /* No solution can exist here.				     */
	    IritFree(Indices);

	    return FALSE;
	}
    }

    IritFree(Indices);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given NumOfMVs multivariates, linearly blend them so that their	     *
* gradients are orthogonal at the center of their domain.		     *
*   This preconditioning of matrices is not only stabilizing the problem but *
* also speeds up the process by improving the domain-clipping and single     *
* solution test progress and chances of success.			     *
*   Assumes the number of multivariates equals to their dimension, and all   *
* MVs share the same function space (Orders, Lengths, and KVs.).	     *
*   Follows the paper "Subdivision methods for solving polynomial equations" *
* by Bernard Mourrain and Jean-Pascal Pavone.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:          Multivariate to manipulate, so their gradients are         *
*                 orthogonal at the center of the domain.		     *
*   NumOfZeroMVs: Number of multivariates.  Must be equal to the dimension   *
*		  of the MVs.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct **:  A locally allocated vector of NumOfMVs                 *
*		  multivariates holding NumOfMVs orthogonalized new          *
*		  constraints, following by NumOfMVs references to the old   *
*		  constraints.						     *
*****************************************************************************/
static MvarMVStruct **MvarMVsOrthogonalizeGrads(MvarMVStruct **MVs,
						int NumOfMVs)
{
    IRIT_STATIC_DATA int
	AllocDim = 0;
    IRIT_STATIC_DATA IrtRType
        *Params = NULL;
    IRIT_STATIC_DATA IrtGnrlMatType
	A = NULL,
	InvA = NULL;
    IRIT_STATIC_DATA MvarMVStruct
	**NewMVs = NULL;
    int i, j, k, l,
	Len = MVAR_CTL_MESH_LENGTH(MVs[0]),
	PtSize = MVAR_NUM_OF_MV_COORD(MVs[0]),
	Dim = MVs[0] -> Dim;
    CagdRType *R, Min, Max, *APtr, W, *Pts, *NewPts;

    assert(Dim == NumOfMVs);

    /* Verify the size of the cache we will use here. */
    if (AllocDim <= Dim) {
        if (Params != NULL) {
	    IritFree(A);
	    IritFree(InvA);
	    IritFree(Params);
	    IritFree(NewMVs);
	}

	AllocDim = Dim * 2;
	A = (IrtGnrlMatType) IritMalloc(sizeof(IrtRType) * IRIT_SQR(AllocDim));
	InvA = (IrtGnrlMatType) IritMalloc(sizeof(IrtRType) *
							   IRIT_SQR(AllocDim));
	Params = (CagdRType *) IritMalloc(sizeof(IrtRType) * AllocDim);
	NewMVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *)
							           * AllocDim);
    }

    /* Compute the center of the domain, from MVs[0]. */
    for (i = 0; i < Dim; i++) {
	MvarGetSubdivParamDomains(MVs[0], &Min, &Max, i);
	Params[i] = (Min + Max) * 0.5;
    }

    /* Place the gradients as rows into matrix A. */
    for (i = 0; i < Dim; i++) {
	R = MvarMVEvalGradient2(MVs[i], Params, NULL);
	CAGD_GEN_COPY(&A[i * Dim], R, sizeof(CagdRType) * Dim);
    }

    /* Invert A. If failed, we have collinear gradients - abort. */
    if (!MatGnrlInverseMatrix(A, InvA, Dim))
        return NULL;

#ifdef DEBUG_TEST_GRAD_MATS
    printf("Gradient Mat:\n");
    MatGnrlPrintMatrix(A, Dim, stdout);
    printf("Gradient Inverse Mat:\n");
    MatGnrlPrintMatrix(InvA, Dim, stdout);
#endif /* DEBUG_TEST_GRAD_MATS */

    /* Copy the Dim MVs to the returned vector's Dim slots.*/
    for (i = 0; i < Dim; i++)
        NewMVs[i] = MvarMVCopy(MVs[i]);

    /* Blend the MVs using the blending dictated by InvA. Note we can either */
    /* have scalar fields, or vector fields with points of dimension Dim+1,  */
    /* having the scalar field and the gradients merged together.	     */
    if (PtSize == 1) {
        for (i = 0; i < Dim; i++) {
	    APtr = &InvA[i * Dim];
	    W = APtr[0];
	    Pts = MVs[0] -> Points[1];

	    for (j = 0, NewPts = NewMVs[i] -> Points[1]; j < Len; j++)
	        *NewPts++ = W * *Pts++;

	    for (k = 1; k < Dim; k++) {
	        W = APtr[k];
		Pts = MVs[k] -> Points[1];

		for (j = 0, NewPts = NewMVs[i] -> Points[1]; j < Len; j++)
		    *NewPts++ += W * *Pts++;
	    }
	}
    }
    else {
        for (l = 1; l <= PtSize; l++) {
	    for (i = 0; i < Dim; i++) {
		APtr = &InvA[i * Dim];
		W = APtr[0];
		Pts = MVs[Dim] -> Points[l];

		for (j = 0, NewPts = NewMVs[i] -> Points[l]; j < Len; j++)
		    *NewPts++ = W * *Pts++;

		for (k = 1; k < Dim; k++) {
		    W = APtr[k];
		    Pts = MVs[Dim + k] -> Points[1];

		    for (j = 0, NewPts = MVs[i] -> Points[l]; j < Len; j++)
		        *NewPts++ += W * *Pts++;
		}
	    }
	}
    }

#ifdef DEBUG_TEST_GRADS_ORTHO
    for (i = 0; i < Dim; i++) {
	R = MvarMVEvalGradient2(NewMVs[i], Params, NULL);
	CAGD_GEN_COPY(&A[i * Dim], R, sizeof(CagdRType) * Dim);
	printf("Gradient %d :: [", i);
	for (j = 0; j < Dim; j++)
	    printf("%f  ", R[j]);
	printf("]\n");
    }
    for (i = 0; i < Dim; i++) {
        for (j = i + 1; j < Dim; j++) {
	    IrtRType
	        R = 0.0;

	    for (k = 0; k < Dim; k++)
	        R += A[i * Dim + k] * A[j * Dim + k];

	    if (!IRIT_APX_EQ(R, 0.0))
	        printf("None orthogonal gradients detected.\n");
	}
    }
#endif /* DEBUG_TEST_GRADS_ORTHO */

    return NewMVs;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes a pair of parallel bounding hyperplanes of a multivariate       *
* constraint.								     *
*   The algorithm:                                                           *
* 1. Evaluate the gradient N0 at the midpoint of the domain of MV (with an   * 
*    additional (d+1)'th coefficient equal to -1), and normalize it.         *
* 2. For each of the control (d+1)-dim points Pi of MV (i.e., MV             *
*    coefficients coerced to R^d+1):                                         *
*    a. Project Pi onto N0 (i.e., inner product with N0), computing	     *
*       min <Pi, N0> and max <Pi, N0> of that projection.                    *
*    b. The bounding (d+1)-dim hyperplanes are:				     *
*         <N0, P - Pmin> = 0 and <N0, P - Pmax> = 0.			     *
* 3. Therefore, the d-hyperplanes (after intersecting with plane x_{d+1}=0)  *
*    are:  <N0_clip, P> = min <Pi, N0> and <N0_clip, P> = max <Pi, N0>,      *
*          where N0_clip is N0 without the (d+1)-coordinate.                 *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:        Multivariate to compute two parallel bounding hyperplanes.    *
*   Coeffs:    The first d-1 coefficients of the d-dimensional bounding	     *
*	       hyperplanes.						     *
*   Coeff0Max: The unit coefficient of first bounding hyperplanes.	     *
*   Coeff0Min: The unit coefficient of second bounding hyperplane.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   None								     *
*****************************************************************************/
static void MVarMVHyperPlanesBound(const MvarMVStruct *MV,
				   CagdRType *Coeffs,
				   CagdRType *Coeff0Max,
				   CagdRType *Coeff0Min)
{
    int i, j, *IndicesVector,
	TotalLength = MVAR_CTL_MESH_LENGTH(MV),
	Dim = MV -> Dim,
	Index = 0;
    MvarVecStruct
        *UnitNormal = MvarVecNew(Dim + 1); 
    MvarVecStruct
        *D1Vec = MvarVecNew(Dim + 1); 
    MvarMVGradientStruct
	*Grad = NULL; 
    CagdRType *R, Min, Max, MinDotProd, MaxDotProd, DotProd,
	*Params = NULL,
	**NodeVectors = NULL;

    /* 1. Evaluate gradient at midpoint (& additional (d+1)-coordinate -1). */
    Params = (CagdRType *) IritMalloc(sizeof(IrtRType) * Dim);
    for (i = 0; i < Dim; i++) {
	MvarGetSubdivParamDomains(MV, &Min, &Max, i);
	Params[i] = (Min + Max) * 0.5;
    }

    if (MVAR_NUM_OF_MV_COORD(MV) == 1) {         /* No precomputed gradient. */
	Grad = MvarMVPrepGradient(MV, FALSE);

        /* Evaluate gradient at midpoint. */
	R = MvarMVEvalGradient(Grad, Params, 0);
        CAGD_GEN_COPY(UnitNormal -> Vec, R, sizeof(CagdRType) * Dim);
    }
    else if (MV -> Dim == MVAR_NUM_OF_MV_COORD(MV) - 1) {
        /* Gradient is embedded in Points[2] to Points[Dim + 1]. */
        
        /* Gradients are saved after the scalar value, in (1+Dim) vector. */
	R = MvarMVEval(MV, Params);
	CAGD_GEN_COPY(UnitNormal -> Vec, &R[2], sizeof(CagdRType) * Dim);
    }
    else {
        MVAR_FATAL_ERROR(MVAR_ERR_DIM_TOO_HIGH);   
	return;
    }
    UnitNormal -> Vec[Dim] = -1.0;
    MvarVecNormalize(UnitNormal);

    IritFree(Params);

    /* Project all control points of the mesh onto the normalized gradient.  */
    /*   The control (d+1)-points of the mesh are (i/k1, j/k2, ..., Pij) for */
    /* Bezier, scaled to the domain of MV.				     */
    /*   For BSpline, they are computed using the BspKnotNodes function.     */
    NodeVectors = (CagdRType **) IritMalloc(sizeof(CagdRType*) * MV -> Dim);

    /* Construct a nodes' vector for each dimension. */
    if (MVAR_IS_BEZIER_MV(MV)) {
        for (j = 0; j < MV -> Dim; j++) {
	    CagdRType
		*NodeVector = (CagdRType *) IritMalloc(sizeof(CagdRType) *
						       MV -> Orders[j]);

	    /* Domain needed for scaling. */
	    MvarGetSubdivParamDomains(MV, &Min, &Max, j);
            for (i = 0; i < MV -> Orders[j]; i++) {
	        NodeVector[i] = (i / (MV -> Orders[j] - 1.0))
							   * (Max - Min) + Min;
            }
            NodeVectors[j] = NodeVector;
        }
    }
    else {
        assert(MVAR_IS_BSPLINE_MV(MV));			    /* Sanity check. */
        for (j = 0; j < MV -> Dim; j++) {
            NodeVectors[j] = BspKnotNodes(MV -> KnotVectors[j],
					  MV -> Lengths[j] + MV -> Orders[j],
					  MV -> Orders[j]);
        }
    }

    /* 2a. Go over all control (d+1)-points, projecting and finding min and  */
    /*     max over N0 and find the extremum.				     */
    IndicesVector = (int *) IritMalloc(sizeof(int) * MV -> Dim);
    IRIT_ZAP_MEM(IndicesVector, sizeof(int) * MV -> Dim);

    /* IndicesVector holds the indices of the current iteration and is	     */
    /* incremented inside the loop.  For example if we're at control point   */
    /* P312 (of a trivariate), the x-coordinate will be NodeVectors[1][3]    */
    /* (for Bezier this will be 3/OrderX), the y-coordinate will be          */
    /* NodesVector[2][1] (for Bezier, 1/OrderY), and the z-coordinate will   */
    /* be NodesVector[3][2] (for Bezier, 2/OrderZ).			     */

    /* Initializing MinDotProd and MaxDotProd with the values corresponding  */
    /* to Index == 0.							     */
    for (i = 0; i < Dim; i++) {
	D1Vec -> Vec[i] = NodeVectors[i][IndicesVector[i]]; 
    }
    /* The d+1 coefficient is the scalar control point. */
    D1Vec -> Vec[Dim] = MV -> Points[1][Index];

    /* We have in D1Vec, the control (d+1)-point P00..0, project it on N0. */
    MinDotProd = MaxDotProd = MvarVecDotProd(D1Vec, UnitNormal);

    /* Loop iterating over all control points of multivariate. */
    for (Index = 1; Index < TotalLength; Index++) {
        /* Getting the IndicesVector. */
        if (!MvarMeshIndicesFromIndex(Index, MV, IndicesVector)) {
            MVAR_FATAL_ERROR(MVAR_ERR_UNDEFINE_ERR);   
	    return;
        }

        /* Coerce the scalar control point to a (d+1)-point.*/
        for (i = 0; i < Dim; i++) {
	    D1Vec -> Vec[i] = NodeVectors[i][IndicesVector[i]]; 
        }
        D1Vec -> Vec[Dim] = MV -> Points[1][Index];     /* The scalar ctlpt. */

        DotProd = MvarVecDotProd(D1Vec, UnitNormal);       /* project on N0. */

        if (DotProd < MinDotProd)
            MinDotProd = DotProd;
        if (DotProd > MaxDotProd)
            MaxDotProd = DotProd;
    }

    *Coeff0Min = MinDotProd; 
    *Coeff0Max = MaxDotProd; 

    /* Copy only first d coordinates of UnitNormal. */ 
    CAGD_GEN_COPY(Coeffs, UnitNormal -> Vec, Dim * sizeof(CagdRType));

    /* Free allocations. */
    if (Grad != NULL)
	MvarMVFreeGradient(Grad);

    for (j = 0; j < MV -> Dim; j++) {
        IritFree(NodeVectors[j]);
    }
    IritFree(NodeVectors);
    IritFree(IndicesVector);
    MvarVecFree(UnitNormal);
    MvarVecFree(D1Vec);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Checks whether a solution of a set of multivariate constraints, might    *
* exist inside the domain, using pairs of parallel hyperplanes.		     *
*   Algorithm:								     *
* 1. Construct bounding hyperplanes, for all MVs.			     *
* 2. Solve for intersections of 2^d vertices of polytope.		     *
* 3. Compare all intersections with every bounding hyperplane (halfspace)    *
*    on the boundary of the domain.				 	     *
* 4. Returns FALSE if all intersections are on one side of domain boundary.  *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:            Multivariates to check their solutions.		     *
*   NumOfZeroMVs:   Size of the vector MVs.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if solution is possible in the domain, FALSE if no      *
*               solution is possible in the domain.                          *
*****************************************************************************/
int MVarMVHyperPlanesTestForSol(MvarMVStruct const * const *MVs,
			        int NumOfZeroMVs)
{
    IRIT_STATIC_DATA int
	AllocDim = 0;
    IRIT_STATIC_DATA CagdRType *Solutions,
	*A = NULL,
	*x = NULL,
	*bMin = NULL,
	*bMax = NULL,
	*bCopy = NULL;
    int i = 0,
	j = 0,
	k = 0,
	Dim = MVs[0] -> Dim,
	PowerTwoDim = (int) pow(2, Dim);

    /* In order to save alloc/free run-time we use static variables. */
    if (AllocDim < Dim) {
        if (AllocDim > 0) {
            IritFree(A);
            IritFree(x);
            IritFree(bMin);
            IritFree(bMax);
            IritFree(bCopy);
            IritFree(Solutions);
        }

	AllocDim = Dim * 2;
        A = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim * AllocDim);
        x = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
        bMin = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
        bMax = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
        bCopy = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
        Solutions = (CagdRType *) IritMalloc(sizeof(CagdRType) * PowerTwoDim
					     * PowerTwoDim * AllocDim);
    }

    /* Construct the bounding hyperplanes of each MV, keeping them as rows   */
    /* A while the scalar coefficients of the pair are kept in bMin/bMax.    */
    for (i = 0; i < NumOfZeroMVs; i++) {
        MVarMVHyperPlanesBound(MVs[i], &A[i * Dim], &bMax[i], &bMin[i]);
    }

    /* Compute QR decomposition of matrix A. */
    if (IritQRUnderdetermined(A, NULL, NULL, Dim, Dim)) {
	return TRUE;    /* Something went wrong - return TRUE (cannot tell). */
    }

    /* Loop over 2^(d) combinations of b vector (000 -> ---, 111 -> +++).    */
    /*   Store the points to check in the end if they are all on one side of */
    /* the box.								     */

    for (i = 0; i < PowerTwoDim; i++) {
        /* Construct relevant copy of b (bMin/bMAx defined by binary rep.).  */
        k = i;
        for (j = 0; j < Dim; ++j) {
	    bCopy[j] = k & 1 ? bMax[j] : bMin[j];
            k >>= 1;
        }

        IritQRUnderdetermined(NULL, &Solutions[i * Dim], bCopy, Dim, Dim);
    }

    for (j = 0; j < Dim; ++j) {
        int CompToBox;
        CagdRType jMin, jMax;

        MvarGetSubdivParamDomains(MVs[0], &jMin, &jMax, j);
        CompToBox = (Solutions[j] > jMax)
				       ? 1 : ((Solutions[j] < jMin) ? -1 : 0);
        /* -1 for left of Min, 1 for right of Max, 0 inbetween Min and Max. */
        if (CompToBox != 0) {
            for (i = 0; i < PowerTwoDim; i++) {
                int CompToBoxTmp = (Solutions[i * Dim + j] > jMax)
			     ? 1 : ((Solutions[i * Dim + j] < jMin) ? -1 : 0);

                if (CompToBox != CompToBoxTmp)
                    break;                
            }

            if (i == PowerTwoDim) {
                /* All intersections are on one side of Box - no solution.  */
		/* is possible inside the domain.			    */
                return FALSE;
            }
        }
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A verification function to test the correctness of the solutions.        M
*   For mostly development/debugging purposes.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:           Input constraints.                                        M
*   NumOfZeroMVs:  Number of (zero only) constraints.                        M
*   Sols:          Linked lists of solutions found.                          M
*   NumerEps:      Numeric tolerance used in the solution.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver                                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosVerifier                                                     M
*****************************************************************************/
void MvarMVsZerosVerifier(MvarMVStruct * const  *MVs,
			  int NumOfZeroMVs,
			  MvarPtStruct *Sols,
			  CagdRType NumerEps)
{
    int i, j;
    MvarPtStruct *Sol;

    NumerEps = IRIT_FABS(NumerEps);

    for (Sol = Sols; Sol != NULL; Sol = Sol -> Pnext) {
	for (i = 0; i < NumOfZeroMVs; i++) {
	    CagdRType
	        *R = MvarMVEval(MVs[i], Sol -> Pt);

	    if (!IRIT_APX_EQ_EPS(R[1], 0.0, NumerEps)) {
	        printf("Invalid solution! [");
		for (j = 0; j < Sol -> Dim; j++)
		    printf("%f%s",
			   Sol -> Pt[j], j == Sol -> Dim - 1 ? "]" : ", ");
		printf("MV[%d] = %g\n", i, R[1]);
	    }
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Filters and purge solution points that are not satisfying some           M
* constraint in MVs.							     *
*                                                                            M
* PARAMETERS:                                                                M
*   MVPts:        Solution points to filter.				     M
*   MVs:          Constraints to test the solution points against.           M
*   Constraints:  Type of constraints.					     M
*   NumOfMVs:     Also number of constraints.				     M
*   Tol:          Tolerance the solution point must satisfy.		     M
*   CanHaveLoops: TRUE if point data can form loops in which case first and  M
*                 last point will be identical (and should not be purged).   M
*   SortSol:      TRUE to sort solution set first, based on the 1st axis.    M
*   InEqOnly:	  TRUE if inequality constraints should be checked only,     M
*                 FALSE if both equality and inequality should be checked.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   Filtered points that satisfy the given constraints.    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver                                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroFilterSolutionSet	                                             M
*****************************************************************************/
MvarPtStruct *MvarZeroFilterSolutionSet(MvarPtStruct *MVPts,
					const MvarMVStruct * const *MVs,
					const MvarConstraintType *Constraints,
					int NumOfMVs,
					IrtRType Tol,
					int CanHaveLoops,
					int SortSol, 
					CagdBType InEqOnly)
{
    int i, l,
        FirstPt = TRUE;
    MvarPtStruct *Pt,
	*OutSet = NULL;

    /* Sort the result based on first axis. */
    if (SortSol)
        MVPts = MvarPtSortListAxis(MVPts, 1);

#ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugMvarFilteringZeroSet, FALSE) {
	    fprintf(stderr, "********************** INITIAL SET OF SOLUTIONS ****************\n");

	    for (Pt = MVPts; Pt != NULL; Pt = Pt -> Pnext) {
	        int i;

		fprintf(stderr, "    [");
		for (i = 0; i < Pt -> Dim; i++)
		    fprintf(stderr, "%9.6f ", Pt -> Pt[i]);
		fprintf(stderr, "]\n");
	    }
	}
    }
#endif /* DEBUG */

    while (MVPts != NULL) {
        IRIT_LIST_POP(Pt, MVPts);

	if (AttrGetIntAttrib(Pt -> Attr, "Similar") == TRUE) {
	    MvarPtFree(Pt);
	}
	else {
	    MvarPtStruct *Pt2;
	    CagdBType
	        PurgePt = FALSE;

	    /* Lets see if we fail any constraint. */
	    for (i = 0; i < NumOfMVs && !PurgePt; i++) {
	        int NumOfCoord = MVAR_NUM_OF_MV_COORD(MVs[i]);
		CagdRType
                    *R = MvarMVEval(MVs[i], Pt -> Pt);

		switch (Constraints[i]) {
		    case MVAR_CNSTRNT_ZERO:
		         if (!InEqOnly && IRIT_ABS(R[1]) > Tol)
			    PurgePt = TRUE;
			 break;
		    case MVAR_CNSTRNT_POSITIVE:
			if (NumOfCoord > 1) { /* Union implementation. */
			    PurgePt = TRUE;
			    for (l = 1; l <= NumOfCoord; l++) {
			        if (R[l] >= 0.0) {
				    PurgePt = FALSE;
				    break;
				}
			    }
			}
			else {
			    if (R[1] < 0.0)
			        PurgePt = TRUE;
			}
			break;
		    case MVAR_CNSTRNT_NEGATIVE:
		        if (NumOfCoord > 1) { /* Union implementation. */
			    PurgePt = TRUE;
			    for (l = 1; l <= NumOfCoord; l++) {
			        if (R[l] <= 0.0) {
				    PurgePt = FALSE;
				    break;
				}
			    }
			}
			else {
			    if (R[1] > 0.0)
			        PurgePt = TRUE;
			}
			break;
	            default:
		        break;
		}
	    }

	    if (PurgePt) {
	        MvarPtFree(Pt);
	    }
	    else {
	        for (Pt2 = MVPts; Pt2 != NULL; Pt2 = Pt2 -> Pnext) {
		    for (i = 0; i < Pt -> Dim; i++) {
		        if (!IRIT_APX_EQ_EPS(Pt -> Pt[i], Pt2 -> Pt[i], Tol))
			    break;
		    }

		    /* Pt and Pt2 are same - mark Pt2 as similar. */
		    if (i >= Pt -> Dim) {
		        if (CanHaveLoops &&
			    FirstPt &&
			    Pt2 -> Pnext == NULL) {
			    /*   Make an exception if the first and last    */
			    /* points are the same and we can have loops.   */
			}
			else {
			    AttrSetIntAttrib(&Pt2 -> Attr, "Similar", TRUE);
			}
		    }
		}

		IRIT_LIST_PUSH(Pt, OutSet);
	    }
	}

	FirstPt = FALSE;
    }

#ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugMvarFilteringZeroSet, FALSE) {
	    fprintf(stderr, "********************** FINAL SET OF SOLUTIONS ****************\n");

	    for (Pt = OutSet; Pt != NULL; Pt = Pt -> Pnext) {
	        int i;

		fprintf(stderr, "    [");
		for (i = 0; i < Pt -> Dim; i++)
		    fprintf(stderr, "%9.6f ", Pt -> Pt[i]);
		fprintf(stderr, "]\n");
	    }
	}
    }
#endif /* DEBUG */

    return CagdListReverse(OutSet);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Test if the given multivariate may satisfy the constraint.               M
*   Examines the positivity/negativity of all coefficients in multivariate.  M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:         Multivariate to examine.                                     M
*   Constraint: Type of constraint - zero, pos., neg.                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:    TRUE if constraint cannot be satisfied, FALSE otherwise.   M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarZeroSolver                                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroMVConstraintFail                                                 M
*****************************************************************************/
CagdBType MvarZeroMVConstraintFail(const MvarMVStruct *MV,
				   MvarConstraintType Constraint)
{
    int l,
	Length = MVAR_CTL_MESH_LENGTH(MV),
	Pos = FALSE,
	NumOfZeros = 0,
	Neg = FALSE;
    CagdRType
	*R = MV -> Points[1];

    switch (Constraint) {
        case MVAR_CNSTRNT_ZERO:
        case MVAR_CNSTRNT_ZERO_SUBDIV:
	    while (--Length >= 0) {
	        Pos |= (*R > 0.0);
		Neg |= (*R < 0.0);
		NumOfZeros += (*R++ == 0.0);
		if (Neg && Pos) {
		    return FALSE;
		}
	    }
	    /* All the coefficients are identically zero.  While not the   */
	    /* mathematically correct answer, we purge this domain!        */
	    if (NumOfZeros == MVAR_CTL_MESH_LENGTH(MV))
	        return TRUE;

	    if (NumOfZeros > 0 && (Neg || Pos))
	        return FALSE;
	    break;
        case MVAR_CNSTRNT_POSITIVE:
	    while (--Length >= 0) {
	        if (*R++ >= 0.0)
		    return FALSE;
	    }
	    /* In inequality constraints, we do allow vector functions     */
	    /* with the semantics that one positive scalar function is     */
	    /* sufficient.  In other words, OR between scalar values.      */
	    if (Length < 0 && MVAR_NUM_OF_MV_COORD(MV) > 1) {
	        for (l = 2; l <= MVAR_NUM_OF_MV_COORD(MV); l++) {
		    R = MV -> Points[l];
		    Length = MVAR_CTL_MESH_LENGTH(MV);

		    while (--Length >= 0) {
		        if (*R++ >= 0.0)
			    return FALSE;
		    }
		}
	    }
	    break;
        case MVAR_CNSTRNT_NEGATIVE:
	    while (--Length >= 0) {
	        if (*R++ <= 0.0)
		    return FALSE;
	    }
	    /* In inequality constraints, we do allow vector functions   */
	    /* with the semantics that one negative scalar function is   */
	    /* sufficient.  In other words, OR between scalar values.    */
	    if (Length < 0 && MVAR_NUM_OF_MV_COORD(MV) > 1) {
	        for (l = 2; l <= MVAR_NUM_OF_MV_COORD(MV); l++) {
		    R = MV -> Points[l];
		    Length = MVAR_CTL_MESH_LENGTH(MV);

		    while (--Length >= 0) {
		        if (*R++ <= 0.0)
			    return FALSE;
		    }
		}
	    }
	    break;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the zeros of the given curve in the given axis.                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Curve:   To compute its zeros.                                           M
*   Axis:    of Crv to seek its zeros: 0 for W, 1 for X, 2 for Y, etc.       M
*   SubdivTol:  Tolerance of the solution.  This tolerance is measured in    M
*	        the parametric space of the curves.			     M
*   NumericTol: Numeric tolerance of a possible numeric improvement stage.   M
*	        The numeric stage is employed if NumericTol < SubdivTol.     M
*   FilterTangencies:  If TRUE, filter out tangencies at the zeros.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdPtStruct *:     List of zeros (parametric locations on CUrve).       M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarSrfZeroSet, MvarMVsZeros0D			                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarCrvZeroSet                                                           M
*****************************************************************************/
CagdPtStruct *MvarCrvZeroSet(const CagdCrvStruct *Curve,
			     int Axis,
			     CagdRType SubdivTol,
			     CagdRType NumericTol,
			     CagdBType FilterTangencies)
{
    int i;
    CagdRType R, TMin, TMax;
    CagdPType
	Trans = { 0.0, 0.0, 0.0 };
    CagdCrvStruct *Crv[4], *Crv0;
    MvarMVStruct *MVs[1];
    MvarPtStruct *MVPts;
    CagdPtStruct *Pts;
    CagdBBoxStruct BBox;

    SymbCrvSplitScalar(Curve, &Crv[0], &Crv[1], &Crv[2], &Crv[3]);
    Crv0 = Crv[Axis];
    Crv[Axis] = NULL;

    /* Converts to a multivariate, after some normalizations. */
    CagdCrvBBox(Crv0, &BBox);
    if ((R = BBox.Max[0] - BBox.Min[0]) == 0.0) {
	for (i = 0; i <= 3; i++) {
	    if (Crv[i] != NULL)
		CagdCrvFree(Crv[i]);
	}
	CagdCrvFree(Crv0);
	return NULL;					/* A constant curve. */
    }
    CagdCrvTransform(Crv0, Trans, 1.0 / R);        /* Scale to a ~unit size. */
    CagdCrvDomain(Crv0, &TMin, &TMax);
    if (CAGD_IS_BSPLINE_CRV(Crv0))
        BspKnotAffineTransOrder2(Crv0 -> KnotVector,  /* Make [0, 1] domain. */
				 Crv0 -> Order, Crv0 -> Length + Crv0 -> Order,
				 0.0, 1.0);

    MVs[0] = MvarCrvToMV(Crv0);

    for (i = 0; i <= 3; i++) {
        if (Crv[i] != NULL)
	    CagdCrvFree(Crv[i]);
    }

    /* Solve the constraints. */
    MVPts = MvarMVsZeros0D(MVs, NULL, 1, SubdivTol, NumericTol);
    Pts = MvarCnvrtMVPtsToPts(MVPts);

    MvarPtFreeList(MVPts);
    MvarMVFree(MVs[0]);

    /* Filter out tangency contact if so desired. */
    if (FilterTangencies) {
        CagdPtStruct *Pt,
	    *FilteredPts = NULL;
	CagdCrvStruct
	    *DCurve = CagdCrvDerive(Crv0),
	    *D2Curve = CagdCrvDerive(DCurve);

	while (Pts != NULL) {
	    CagdRType D, D2;

	    IRIT_LIST_POP(Pt, Pts);
	    D = CagdCrvEval(DCurve, Pt -> Pt[0])[0];
	    D2 = CagdCrvEval(D2Curve, Pt -> Pt[0])[0];

	    if (IRIT_ABS(D) < SubdivTol * 10.0 &&
		IRIT_ABS(D2) > SubdivTol * 10.0) {
		CagdPtFree(Pt);			      /* A tangency contact. */
	    }
	    else {
	        /* Convert [0, 1] to original curve's domain. */
	        Pt -> Pt[0] = Pt -> Pt[0] * (TMax - TMin) + TMin;
		IRIT_LIST_PUSH(Pt, FilteredPts);
	    }
	}

	CagdCrvFree(Crv0);
	CagdCrvFree(DCurve);
	CagdCrvFree(D2Curve);

        return CagdListReverse(FilteredPts);
    }
    else {
	CagdCrvFree(Crv0);
        return Pts;
    }
}
