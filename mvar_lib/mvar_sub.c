/******************************************************************************
* Mvar_Sub.c - Computes subdivision of multi-variates.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, May. 97.					      *
******************************************************************************/

#include "mvar_loc.h"

#ifdef DEBUG
static void MvarDbgPrintAlphaMat(BspKnotAlphaCoeffStruct *A);
#endif /* DEBUG */

static void MvarMVBzrSubdivCtlMesh(const MvarMVStruct *MV,
				   MvarMVStruct *LMV,
				   MvarMVStruct *RMV,
				   CagdRType t,
				   MvarMVDirType Dir);
static void MvarMVBzrSubdivCtlMeshOneSide(const MvarMVStruct *MV,
					  MvarMVStruct *LMV,
					  MvarMVStruct *RMV,
					  CagdRType t,
					  MvarMVDirType Dir,
					  IrtBType LeftSide);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a multi-variate, subdivides it at parameter value t in direction     M
* Dir.									     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:       Multi-Variate to subdivide.                                    M
*   t:        Parameter to subdivide at.                                     M
*   Dir:      Direction of subdivision.                                      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarBspMVSubdivAtParam, MvarBzrMVSubdivAtParam		             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVSubdivAtParam, multi-variates                                      M
*****************************************************************************/
MvarMVStruct *MvarMVSubdivAtParam(const MvarMVStruct *MV,
				  CagdRType t,
				  MvarMVDirType Dir)
{
    switch (MV -> GType) {
	case MVAR_BEZIER_TYPE:
	    return MvarBzrMVSubdivAtParam(MV, t, Dir);
	case MVAR_BSPLINE_TYPE:
	    return MvarBspMVSubdivAtParam(MV, t, Dir);
	default:
	    MVAR_FATAL_ERROR(MVAR_ERR_UNDEF_GEOM);
	    break;
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Given a multi-variate, Bezier subdivides it at parameter value t in        *
* direction Dir.  Note we could invoke this function also with a B-spline    *
* multivariate that has a Bezier knot sequence in Dir.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:       Bezier Multi-Variate to subdivide.                             *
*   t:        Parameter to subdivide at.                                     *
*   Dir:      Direction of subdivision.                                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. *
*****************************************************************************/
static void MvarMVBzrSubdivCtlMesh(const MvarMVStruct *MV,
				   MvarMVStruct *LMV,
				   MvarMVStruct *RMV,
				   CagdRType t,
				   MvarMVDirType Dir)
{
    CagdBType
	IsNotRational = !MVAR_IS_RATIONAL_MV(MV);
    int i, j, l, *RIndices, RIndex0,
	MaxCoord = MVAR_NUM_OF_MV_COORD(MV),
	Length = MV -> Lengths[Dir],
	Len = MVAR_CTL_MESH_LENGTH(MV),
        Step = MVAR_NEXT_DIM(MV, Dir),
	Dim = MV -> Dim,
        IsScalar = IsNotRational && MaxCoord == 1;
    CagdRType
	* const *Points = MV -> Points,
	**LPoints = LMV -> Points,
	**RPoints = RMV -> Points,
	t1 = 1.0 - t;

    /* Copy Points into R/LPoints, so we can apply the recursive algo. to. */
    for (j = IsNotRational; j <= MaxCoord; j++) {
        IRIT_GEN_COPY(RPoints[j], Points[j], Len * sizeof(CagdRType));
        IRIT_GEN_COPY(LPoints[j], Points[j], Len * sizeof(CagdRType));
    }

    /* Do the control mesh's subdivision. */
    RIndices = (int *) IritMalloc(sizeof(int) * Dim);
    IRIT_ZAP_MEM(RIndices, sizeof(int) * Dim);
    RIndex0 = 0;
    if (IsScalar) {
        CagdRType
	    *L1Points = LPoints[1],
	    *R1Points = RPoints[1];

        do {
	    int LIndex = RIndex0;

	    for (i = 1; i < Length; i++) {
	        CagdRType
		    *RPts = &R1Points[RIndex0];

		for (l = Length - i; l-- > 0; ) {
		    *RPts = *RPts * t1 + RPts[Step] * t;
		    RPts += Step;
		}

		/* Copy temporary result to LPoints: */
		L1Points[LIndex += Step] = R1Points[RIndex0];
	    }
	}
	while (MVAR_INC_SKIP_MESH_INDICES(RMV, RIndices, Dir, RIndex0));
    }
    else {
        do {
	    int RIndex, RIndex1,
	        LIndex = RIndex0;

	    for (i = 1; i < Length; i++) {
	        RIndex = RIndex0;

		for (l = Length - i; l-- > 0; ) {
		    RIndex1 = RIndex + Step;

		    for (j = IsNotRational; j <= MaxCoord; j++)
		        RPoints[j][RIndex] = RPoints[j][RIndex] * t1 +
		                             RPoints[j][RIndex1] * t;

		    RIndex = RIndex1;
		}

		/* Copy temporary result to LMV: */
		for (j = IsNotRational, LIndex += Step; j <= MaxCoord; j++)
		    LPoints[j][LIndex] = RPoints[j][RIndex0];
	    }
	}
	while (MVAR_INC_SKIP_MESH_INDICES(RMV, RIndices, Dir, RIndex0));
    }

    IritFree(RIndices);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a Bezier multi-variate, subdivides it at parameter value t in        M
* direction Dir.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:       Bezier Multi-Variate to subdivide.                             M
*   t:        Parameter to subdivide at.                                     M
*   Dir:      Direction of subdivision.                                      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarBspMVSubdivAtParam, MvarMVSubdivAtParam			             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarBzrMVSubdivAtParam, multi-variates                                   M
*****************************************************************************/
MvarMVStruct *MvarBzrMVSubdivAtParam(const MvarMVStruct *MV,
				     CagdRType t,
				     MvarMVDirType Dir)
{
    MvarMVStruct *RMV, *LMV;

    if (Dir < 0 || Dir >= MV -> Dim) {
	MVAR_FATAL_ERROR(MVAR_ERR_DIR_NOT_VALID);
	return NULL;
    }

    if (!MVAR_IS_BEZIER_MV(MV)) {
	MVAR_FATAL_ERROR(MVAR_ERR_UNDEF_MVAR);
	return NULL;
    }

    LMV = MvarBzrMVNew(MV -> Dim, MV -> Lengths, MV -> PType);
    RMV = MvarBzrMVNew(MV -> Dim, MV -> Lengths, MV -> PType);

    MvarMVBzrSubdivCtlMesh(MV, LMV, RMV, t, Dir);

    if (MvarMVAuxDomainSlotCopy(LMV, MV) &&
	MvarMVAuxDomainSlotCopy(RMV, MV)) {
        MvarMVAuxDomainSlotSetRel(LMV, 0.0, t, Dir);
        MvarMVAuxDomainSlotSetRel(RMV, t, 1.0, Dir);
    }

    LMV -> Pnext = RMV;

    CAGD_PROPAGATE_ATTR(RMV, MV);
    CAGD_PROPAGATE_ATTR(LMV, MV);

    return LMV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a Bspline multi-variate, subdivides it at parameter value t in       M
* direction Dir.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:       BsplineMulti-Variate to subdivide.                             M
*   t:        Parameter to subdivide at.                                     M
*   Dir:      Direction of subdivision.                                      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVSubdivAtParam, MvarBzrMVSubdivAtParam		                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarBspMVSubdivAtParam, multi-variates                                   M
*****************************************************************************/
MvarMVStruct *MvarBspMVSubdivAtParam(const MvarMVStruct *MV,
				     CagdRType t,
				     MvarMVDirType Dir)
{
    CagdBType
	IsNotRational = !MVAR_IS_RATIONAL_MV(MV);
    int i, j, KVLen, Index1, Index2, *Indices, Mult, RLength, LLength, Index,
	MaxCoord = MVAR_NUM_OF_MV_COORD(MV),
	Order = MV -> Orders[Dir],
	Length = MV -> Lengths[Dir];
    CagdRType *RefKV;
    MvarMVStruct *RMV, *LMV,
	*CpMV = NULL;
    BspKnotAlphaCoeffStruct *A;

    if (Dir < 0 || Dir >= MV -> Dim) {
	MVAR_FATAL_ERROR(MVAR_ERR_DIR_NOT_VALID);
	return NULL;
    }

    if (MVAR_IS_BSPLINE_MV(MV)) {
        if (MvarBspMVIsPeriodic(MV)) {
	    MV = CpMV = MvarCnvrtPeriodic2FloatMV(MV);
	    Length = MV -> Lengths[Dir];
	}

	RefKV = MV -> KnotVectors[Dir];
	KVLen = Order + Length;
	Index1 = BspKnotLastIndexL(RefKV, KVLen, t);
	if (Index1 + 1 < Order)
	    Index1 = Order - 1;
	Index2 = BspKnotFirstIndexG(RefKV, KVLen, t);
	if (Index2 > Length)
	    Index2 = Length;
	Mult = Order - 1 - (Index2 - Index1 - 1);

	MV -> Lengths[Dir] = Index1 + 1;
	LMV = MvarBspMVNew(MV -> Dim, MV -> Lengths, MV -> Orders, MV -> PType);
	MV -> Lengths[Dir] = Length - Index2 + Order;
	RMV = MvarBspMVNew(MV -> Dim, MV -> Lengths, MV -> Orders, MV -> PType);
	MV -> Lengths[Dir] = Length;

	/* Update the new knot vectors. */
	for (i = 0; i < MV -> Dim; i++) {
	    if (i == Dir) {
		CAGD_GEN_COPY(LMV -> KnotVectors[i],
			      MV -> KnotVectors[i],
			      sizeof(CagdRType) * (Index1 + 1));

		/* Close the knot vector with multiplicity Order: */
		for (j = Index1 + 1; j <= Index1 + Order; j++)
		    LMV -> KnotVectors[i][j] = t;

		CAGD_GEN_COPY(&RMV -> KnotVectors[i][Order],
			      &MV -> KnotVectors[i][Index2],
			      sizeof(CagdRType) *
				  (Length + Order - Index2));

		/* Make sure knot vector starts with multiplicity Order: */
		for (j = 0; j < Order; j++)
		    RMV -> KnotVectors[i][j] = t;
	    }
	    else {
		/* And copy the other direction(s)' knot vectors. */
		CAGD_GEN_COPY(LMV -> KnotVectors[i],
			      MV -> KnotVectors[i],
			      sizeof(CagdRType) * (MV -> Orders[i] +
						   MV -> Lengths[i]));
		CAGD_GEN_COPY(RMV -> KnotVectors[i],
			      MV -> KnotVectors[i],
			      sizeof(CagdRType) * (MV -> Orders[i] +
						   MV -> Lengths[i]));
	    }
	}
    }
    else {
	MVAR_FATAL_ERROR(MVAR_ERR_UNDEF_MVAR);
	return NULL;
    }

    if (BspKnotHasBezierKV(RefKV, Length, Order)) {
        CagdRType Min, Max;

        /* Do a Bezier control mesh subdivision. */
	MvarMVDomain(MV, &Min, &Max, Dir);

        MvarMVBzrSubdivCtlMesh(MV, LMV, RMV, (t - Min) / (Max - Min), Dir);
    }
    else {
        /* Do the B-spline control mesh subdivision. */
        if (Mult > 0) {
	    CagdRType Min, Max,
	        *NewKV = (CagdRType *) IritMalloc(sizeof(CagdRType) * Mult);

	    MvarMVDomain(MV, &Min, &Max, Dir);

	    CAGD_DOMAIN_T_VERIFY(t, Min, Max);
	    if (t == Max)
	        t -= CAGD_DOMAIN_IRIT_EPS;
	    for (i = 0; i < Mult; i++)
	        NewKV[i] = t;
	    A = BspKnotEvalAlphaCoefMerge(Order, RefKV, Length, NewKV, Mult,
					  FALSE);
	    IritFree(NewKV);
	}
	else
	    A = BspKnotEvalAlphaCoefMerge(Order, RefKV, Length, NULL, 0,
					  FALSE);

#ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintAlphaMat, FALSE)
		MvarDbgPrintAlphaMat(A);
	}
#endif /* DEBUG */

        /* Note that Mult can be negative in cases where original	    */
	/* multiplicity was order or more and we need to compensate here,   */
	/* since Alpha matrix will be just a unit matrix then.		    */
	Mult = Mult >= 0 ? 0 : -Mult;

	/* Update the control mesh. */
	Indices = (int *) IritMalloc(sizeof(int) * MV -> Dim);

	LLength = LMV -> Lengths[Dir];
	RLength = RMV -> Lengths[Dir];

	/* Do the left hand side. */
	IRIT_ZAP_MEM(Indices, sizeof(int) * MV -> Dim);
	Index = 0;
	do {
	    int LIndex = MvarGetPointsMeshIndices(LMV, Indices);

	    for (j = IsNotRational; j <= MaxCoord; j++) {
	        CagdRType
		    *LPts = &LMV -> Points[j][LIndex],
		    *Pts = &MV -> Points[j][Index];

		BspKnotAlphaLoopBlendStep(A, 0, LLength, Pts,
					  MVAR_NEXT_DIM(MV, Dir), -1, LPts,
					  MVAR_NEXT_DIM(LMV, Dir));
	    }
	}
	while (MVAR_INC_SKIP_MESH_INDICES(MV, Indices, Dir, Index));

	/* Do the right hand side. */
	IRIT_ZAP_MEM(Indices, sizeof(int) * MV -> Dim);
	Index = 0;
	do {
	    int RIndex = MvarGetPointsMeshIndices(RMV, Indices),
	        Offset = LLength - 1;

	    for (j = IsNotRational; j <= MaxCoord; j++) {
	        CagdRType
		    *RPts = &RMV -> Points[j][RIndex],
		    *Pts = &MV -> Points[j][Index];

		BspKnotAlphaLoopBlendStep(A, Offset + Mult,
					  RLength + Offset + Mult, Pts,
					  MVAR_NEXT_DIM(MV, Dir), -1, RPts,
					  MVAR_NEXT_DIM(RMV, Dir));
	    }
	}
	while (MVAR_INC_SKIP_MESH_INDICES(MV, Indices, Dir, Index));

	IritFree(Indices);

	BspKnotFreeAlphaCoef(A);
    }

    if (CpMV != NULL)
	MvarMVFree(CpMV);

    LMV -> Pnext = RMV;
    return LMV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a multi-variate, subdivides it at parameter value t in direction     M
* Dir.									     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:        Multi-Variate to subdivide.                                   M
*   t:         Parameter to subdivide at.                                    M
*   Dir:       Direction of subdivision.                                     M
*   LeftSide:  TRUE to only fetch left half, FALSE to fetch right half.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarBspMVSubdivAtParamOneSide, MvarBzrMVSubdivAtParamOneSide,	     M
*   MvarMVSubdivAtParam	                                   	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVSubdivAtParamOneSide, multi-variates                               M
*****************************************************************************/
MvarMVStruct *MvarMVSubdivAtParamOneSide(const MvarMVStruct *MV,
					 CagdRType t,
					 MvarMVDirType Dir,
					 IrtBType LeftSide)
{
    switch (MV -> GType) {
	case MVAR_BEZIER_TYPE:
	    return MvarBzrMVSubdivAtParamOneSide(MV, t, Dir, LeftSide);
	case MVAR_BSPLINE_TYPE:
	    return MvarBspMVSubdivAtParamOneSide(MV, t, Dir, LeftSide);
	default:
	    MVAR_FATAL_ERROR(MVAR_ERR_UNDEF_GEOM);
	    break;
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Given a multi-variate, Bezier subdivides it at parameter value t in        *
* direction Dir.  Note we could invoke this function also with a B-spline    *
* multivariate that has a Bezier knot sequence in Dir.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:        Bezier Multi-Variate to subdivide.                            *
*   t:         Parameter to subdivide at.                                    *
*   Dir:       Direction of subdivision.                                     *
*   LeftSide:  TRUE to only fetch left half, FALSE to fetch right half.      *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. *
*****************************************************************************/
static void MvarMVBzrSubdivCtlMeshOneSide(const MvarMVStruct *MV,
					  MvarMVStruct *LMV,
					  MvarMVStruct *RMV,
					  CagdRType t,
					  MvarMVDirType Dir,
					  IrtBType LeftSide)
{
    CagdBType
	IsNotRational = !MVAR_IS_RATIONAL_MV(MV);
    int i, l,
	MaxCoord = MVAR_NUM_OF_MV_COORD(MV),
	Length = MV -> Lengths[Dir],
	Len = MVAR_CTL_MESH_LENGTH(MV),
        Step = MVAR_NEXT_DIM(MV, Dir),
	Dim = MV -> Dim,
        IsScalar = IsNotRational && MaxCoord == 1;
    CagdRType
	* const *Points = MV -> Points,
	t1 = 1.0 - t;

    if (!IsScalar) {
        /* Do the old non-optimal way. */
        if (LeftSide) {
	    assert(RMV == NULL);
	    RMV = MvarMVCopy(MV);

	    MvarMVBzrSubdivCtlMesh(MV, LMV, RMV, t, Dir);

	    MvarMVFree(RMV);
	}
	else {
	    assert(LMV == NULL);
	    LMV = MvarMVCopy(MV);

	    MvarMVBzrSubdivCtlMesh(MV, LMV, RMV, t, Dir);

	    MvarMVFree(LMV);
	}
	return;
    }

    /* Do the control mesh's subdivision. */
    if (LeftSide) {
        int *LIndices = (int *) IritMalloc(sizeof(int) * Dim),
	    LIndex0 = 0;
	CagdRType
	    *L1Points = LMV -> Points[1];

	IRIT_ZAP_MEM(LIndices, sizeof(int) * Dim);
        IRIT_GEN_COPY(L1Points, Points[1], Len * sizeof(CagdRType));

	do {
	    CagdRType *LPts,
		*LPts0 = &L1Points[LIndex0];

	    for (i = 1; i < Length; i++) {
		LPts = LPts0 + Step * (Length - i - 1);
		for (l = Length - i; l-- > 0; ) {
		    LPts[Step] = *LPts * t1 + LPts[Step] * t;
		    LPts -= Step;
		}
		LPts0 += Step;
	    }
	}
	while (MVAR_INC_SKIP_MESH_INDICES(LMV, LIndices, Dir, LIndex0));

	IritFree(LIndices);
    }
    else {
        int *RIndices = (int *) IritMalloc(sizeof(int) * Dim),
	    RIndex0 = 0;
	CagdRType
	    *R1Points = RMV -> Points[1];

	IRIT_ZAP_MEM(RIndices, sizeof(int) * Dim);
        IRIT_GEN_COPY(R1Points, Points[1], Len * sizeof(CagdRType));

	do {
	    CagdRType *RPts,
		*RPts0 = &R1Points[RIndex0];

	    for (i = 1; i < Length; i++) {
		for (RPts = RPts0, l = Length - i; l-- > 0; ) {
		    *RPts = *RPts * t1 + RPts[Step] * t;
		    RPts += Step;
		}
	    }
	}
	while (MVAR_INC_SKIP_MESH_INDICES(RMV, RIndices, Dir, RIndex0));

	IritFree(RIndices);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a Bezier multi-variate, subdivides it at parameter value t in        M
* direction Dir.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:        Bezier Multi-Variate to subdivide.                            M
*   t:         Parameter to subdivide at.                                    M
*   Dir:       Direction of subdivision.                                     M
*   LeftSide:  TRUE to only fetch left half, FALSE to fetch right half.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarBspMVSubdivAtParamOneSide, MvarMVSubdivAtParam		             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarBzrMVSubdivAtParamOneSide, multi-variates                            M
*****************************************************************************/
MvarMVStruct *MvarBzrMVSubdivAtParamOneSide(const MvarMVStruct *MV,
					    CagdRType t,
					    MvarMVDirType Dir,
					    IrtBType LeftSide)
{
    MvarMVStruct *RMV, *LMV;

    if (Dir < 0 || Dir >= MV -> Dim) {
	MVAR_FATAL_ERROR(MVAR_ERR_DIR_NOT_VALID);
	return NULL;
    }

    if (!MVAR_IS_BEZIER_MV(MV)) {
	MVAR_FATAL_ERROR(MVAR_ERR_UNDEF_MVAR);
	return NULL;
    }

    if (LeftSide) {
        LMV = MvarBzrMVNew(MV -> Dim, MV -> Lengths, MV -> PType);
	RMV = NULL;

	if (MvarMVAuxDomainSlotCopy(LMV, MV))
	    MvarMVAuxDomainSlotSetRel(LMV, 0.0, t, Dir);
    }
    else {
        LMV = NULL;
        RMV = MvarBzrMVNew(MV -> Dim, MV -> Lengths, MV -> PType);

	if (MvarMVAuxDomainSlotCopy(RMV, MV))
	    MvarMVAuxDomainSlotSetRel(RMV, t, 1.0, Dir);
    }

    MvarMVBzrSubdivCtlMeshOneSide(MV, LMV, RMV, t, Dir, LeftSide);

    if (LeftSide) {
        CAGD_PROPAGATE_ATTR(LMV, MV);
	return LMV;
    }
    else {
        CAGD_PROPAGATE_ATTR(RMV, MV);
	return RMV;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Given a Bspline multi-variate, subdivides it at parameter value t in       M
* direction Dir.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:        BsplineMulti-Variate to subdivide.                            M
*   t:         Parameter to subdivide at.                                    M
*   Dir:       Direction of subdivision.                                     M
*   LeftSide:  TRUE to only fetch left half, FALSE to fetch right half.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A list of two multi-variates, result of the subdivision. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVSubdivAtParam, MvarBzrMVSubdivAtParam		                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarBspMVSubdivAtParamOneSide, multi-variates                            M
*****************************************************************************/
MvarMVStruct *MvarBspMVSubdivAtParamOneSide(const MvarMVStruct *MV,
					    CagdRType t,
					    MvarMVDirType Dir,
					    IrtBType LeftSide)
{
    CagdBType
	IsNotRational = !MVAR_IS_RATIONAL_MV(MV);
    int i, j, KVLen, Index1, Index2, *Indices, Mult, RLength, LLength, Index,
	MaxCoord = MVAR_NUM_OF_MV_COORD(MV),
	Order = MV -> Orders[Dir],
	Length = MV -> Lengths[Dir];
    CagdRType *RefKV;
    MvarMVStruct *RMV, *LMV,
	*CpMV = NULL;
    BspKnotAlphaCoeffStruct *A;

    if (Dir < 0 || Dir >= MV -> Dim) {
	MVAR_FATAL_ERROR(MVAR_ERR_DIR_NOT_VALID);
	return NULL;
    }

    if (MVAR_IS_BSPLINE_MV(MV)) {
        if (MvarBspMVIsPeriodic(MV)) {
	    MV = CpMV = MvarCnvrtPeriodic2FloatMV(MV);
	    Length = MV -> Lengths[Dir];
	}

	RefKV = MV -> KnotVectors[Dir];
	KVLen = Order + Length;
	Index1 = BspKnotLastIndexL(RefKV, KVLen, t);
	if (Index1 + 1 < Order)
	    Index1 = Order - 1;
	Index2 = BspKnotFirstIndexG(RefKV, KVLen, t);
	if (Index2 > Length)
	    Index2 = Length;
	Mult = Order - 1 - (Index2 - Index1 - 1);

	if (LeftSide) {
	    MV -> Lengths[Dir] = Index1 + 1;
	    LMV = MvarBspMVNew(MV -> Dim, MV -> Lengths,
			       MV -> Orders, MV -> PType);
	    RMV = NULL;
	}
	else {
	    MV -> Lengths[Dir] = Length - Index2 + Order;
	    LMV = NULL;
	    RMV = MvarBspMVNew(MV -> Dim, MV -> Lengths,
			       MV -> Orders, MV -> PType);
	}
	MV -> Lengths[Dir] = Length;

	/* Update the new knot vectors. */
	for (i = 0; i < MV -> Dim; i++) {
	    if (i == Dir) {
	        if (LeftSide) {
		    CAGD_GEN_COPY(LMV -> KnotVectors[i],
				  MV -> KnotVectors[i],
				  sizeof(CagdRType) * (Index1 + 1));

		    /* Close the knot vector with multiplicity Order: */
		    for (j = Index1 + 1; j <= Index1 + Order; j++)
		        LMV -> KnotVectors[i][j] = t;
		}
		else {
		    CAGD_GEN_COPY(&RMV -> KnotVectors[i][Order],
				  &MV -> KnotVectors[i][Index2],
				  sizeof(CagdRType) *
				            (Length + Order - Index2));

		    /* Make sure knot vector starts with multiplicity Order: */
		    for (j = 0; j < Order; j++)
		        RMV -> KnotVectors[i][j] = t;
		}
	    }
	    else {
	        if (LeftSide) {
		    /* And copy the other direction(s)' knot vectors. */
		    CAGD_GEN_COPY(LMV -> KnotVectors[i],
				  MV -> KnotVectors[i],
				  sizeof(CagdRType) * (MV -> Orders[i] +
						       MV -> Lengths[i]));
		}
		else {
		    CAGD_GEN_COPY(RMV -> KnotVectors[i],
				  MV -> KnotVectors[i],
				  sizeof(CagdRType) * (MV -> Orders[i] +
						       MV -> Lengths[i]));
		}
	    }
	}
    }
    else {
	MVAR_FATAL_ERROR(MVAR_ERR_UNDEF_MVAR);
	return NULL;
    }

    if (BspKnotHasBezierKV(RefKV, Length, Order)) {
        CagdRType Min, Max;

        /* Do a Bezier control mesh subdivision. */
	MvarMVDomain(MV, &Min, &Max, Dir);

        MvarMVBzrSubdivCtlMeshOneSide(MV, LMV, RMV, (t - Min) / (Max - Min),
				      Dir, LeftSide);

#	ifdef DEBUG_TEST_ON_SIDE_SUBDIV
	{
	    MvarMVStruct *LMV1;

	    /* Compare to regular subdivision. */
	    LMV1 = MvarBspMVSubdivAtParam(MV, t, Dir);
	    if (LeftSide) {
		LMV1 -> Pnext = NULL;
	        if (!MvarMVsSame(LMV, LMV1, IRIT_EPS)) {
		    MvarDbg(LMV);
		    MvarDbg(LMV1);
		}
	    }
	    else {
	        if (!MvarMVsSame(RMV, LMV1 -> Pnext, IRIT_EPS)) {
		    MvarDbg(RMV);
		    MvarDbg(LMV1 -> Pnext);
		}
	    }
	}
#	endif /* DEBUG_TEST_ON_SIDE_SUBDIV */
    }
    else {
        /* Do the B-spline control mesh subdivision. */
        if (Mult > 0) {
	    CagdRType Min, Max,
	        *NewKV = (CagdRType *) IritMalloc(sizeof(CagdRType) * Mult);

	    MvarMVDomain(MV, &Min, &Max, Dir);

	    CAGD_DOMAIN_T_VERIFY(t, Min, Max);
	    if (t == Max)
	        t -= CAGD_DOMAIN_IRIT_EPS;
	    for (i = 0; i < Mult; i++)
	        NewKV[i] = t;
	    A = BspKnotEvalAlphaCoefMerge(Order, RefKV, Length, NewKV, Mult,
					  FALSE);
	    IritFree(NewKV);
	}
	else
	    A = BspKnotEvalAlphaCoefMerge(Order, RefKV, Length, NULL, 0,
					  FALSE);

#ifdef DEBUG
	{
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintAlphaMat, FALSE)
		MvarDbgPrintAlphaMat(A);
	}
#endif /* DEBUG */

        /* Note that Mult can be negative in cases where original	    */
	/* multiplicity was order or more and we need to compensate here,   */
	/* since Alpha matrix will be just a unit matrix then.		    */
	Mult = Mult >= 0 ? 0 : -Mult;

	/* Update the control mesh. */
	Indices = (int *) IritMalloc(sizeof(int) * MV -> Dim);

	if (LeftSide) {
	    /* Do the left hand side. */
	    LLength = LMV -> Lengths[Dir];
	    IRIT_ZAP_MEM(Indices, sizeof(int) * MV -> Dim);
	    Index = 0;
	    do {
	        int LIndex = MvarGetPointsMeshIndices(LMV, Indices);

		for (j = IsNotRational; j <= MaxCoord; j++) {
		    CagdRType
		        *LPts = &LMV -> Points[j][LIndex],
		        *Pts = &MV -> Points[j][Index];

		    BspKnotAlphaLoopBlendStep(A, 0, LLength, Pts,
					      MVAR_NEXT_DIM(MV, Dir), -1, LPts,
					      MVAR_NEXT_DIM(LMV, Dir));
		}
	    }
	    while (MVAR_INC_SKIP_MESH_INDICES(MV, Indices, Dir, Index));
	}
	else {
	    /* Do the right hand side. */
	    RLength = RMV -> Lengths[Dir];
	    LLength = Index1 + 1;
	    IRIT_ZAP_MEM(Indices, sizeof(int) * MV -> Dim);
	    Index = 0;
	    do {
	        int RIndex = MvarGetPointsMeshIndices(RMV, Indices),
		    Offset = LLength - 1;

		for (j = IsNotRational; j <= MaxCoord; j++) {
		    CagdRType
		        *RPts = &RMV -> Points[j][RIndex],
		        *Pts = &MV -> Points[j][Index];

		    BspKnotAlphaLoopBlendStep(A, Offset + Mult,
					      RLength + Offset + Mult, Pts,
					      MVAR_NEXT_DIM(MV, Dir), -1, RPts,
					      MVAR_NEXT_DIM(RMV, Dir));
		}
	    }
	    while (MVAR_INC_SKIP_MESH_INDICES(MV, Indices, Dir, Index));
	}

	IritFree(Indices);

	BspKnotFreeAlphaCoef(A);
    }

    if (CpMV != NULL)
	MvarMVFree(CpMV);

    return LeftSide ? LMV : RMV;
}

#ifdef DEBUG

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prints the content of the alpha matrix.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   A:        Alpha matrix to print.                                         *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void MvarDbgPrintAlphaMat(BspKnotAlphaCoeffStruct *A)
{
    int i, j;

    IRIT_INFO_MSG_PRINTF("Order = %d, Length = %d\n",
			 A -> Order, A -> Length);
    for (i = 0; i < A -> Length; i++) {
	for (j = 0; j < A -> RefLength; j++)
	    IRIT_INFO_MSG_PRINTF(" %9.5g", A -> Rows[i][j]);
	IRIT_INFO_MSG("\n");
    }

    IRIT_INFO_MSG("    ");
    for (j = 0; j < A -> RefLength; j++)
	    IRIT_INFO_MSG_PRINTF(" %3d %3d |",
				 A -> ColIndex[j], A -> ColLength[j]);
	IRIT_INFO_MSG("\n");
}

#endif /* DEBUG */
