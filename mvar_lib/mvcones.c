/******************************************************************************
* MvCones.c - Tools to construct and intersect MV (anti-)cones & vectors.     *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Iddo Haniel and Gershon Elber, May 2005.			      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/extra_fn.h"

/* #define MV_CONES_TEST_OVERLAPS */

#ifdef MV_CONES_TEST_OVERLAPS
int MVConesNoOverlaps = 0,
    MVConesOverlapTests = 0;
#endif /* MV_CONES_TEST_OVERLAPS */

#define MVAR_2CONES_EXPAND_FACTOR	2
#define MVAR_2CONES_MAX_CONE_ANGLE	0.99 /* Max cone angle to do 2cones. */

static int MVMinSpanConeWithActiveVecs(MvarVecStruct *MVVecs,
				       int NumOfVecs,
				       MvarVecStruct **ActiveVecs,
				       int NumOfActiveVecs,
				       MvarNormalConeStruct *MVCone);
static CagdRType MVarMVMaximalDeviation(const MvarVecStruct *ConeAxis,
					CagdRType * const *GradPoints,
					int NumPoints,
					int *MaxDevIndex,
					CagdRType *MinLength,
					CagdRType *MaxLength);

/* Normal cone operations */
static MvarNormalConeStruct *MvarExprTreeNormalConeSum(
				       const MvarNormalConeStruct* ConeF,
				       const MvarNormalConeStruct* ConeG,
				       int Dim);
static MvarNormalConeStruct *MvarExprTreeNormalConeSub(
				       const MvarNormalConeStruct* ConeF,
				       const MvarNormalConeStruct* ConeG,
				       int Dim);
static MvarNormalConeStruct *MvarExprTreeNormalConeMul(
				       const MvarNormalConeStruct* ConeF,
				       const MvarNormalConeStruct* ConeG,
				       const MvarBBoxStruct* BBoxF,
				       const MvarBBoxStruct* BBoxG,
				       int Dim);
static MvarNormalConeStruct *MvarExprTreeNormalConeScale(
					const MvarNormalConeStruct* ConeF,
					const MvarBBoxStruct* BBoxGPrime,
					int Dim);
static MvarVecStruct *HyperplaneOrthoSystem(const MvarVecStruct* v);
static void MvarConesAssembleB(int CurrentPower,
			       int Dim,
			       const CagdRType *b,
			       CagdRType *bCopy);


#ifdef DEBUG
static void MvarMinSpanConeExhaustive(MvarVecStruct *MVVecs,
				      int NumOfVecs,
				      MvarVecStruct **MVActiveVecs,
				      int Level,
				      int First,
				      MvarNormalConeStruct *MVCone);
#endif /* DEBUG */

#define ADVANCE_VECTOR_COUNT(Vector, VecSizeStore, Zeroes) { \
    *VecSizeStore = MvarVecLength(Vector); \
    if (IRIT_UEPS < *VecSizeStore) \
        MvarVecScale(Vector++, 1.0 / *VecSizeStore++); \
    else \
        Zeroes++; \
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a multivariate normal cone structure.                         M
*                                                                            *
* PARAMETERS:                                                                M
*   Dim:  Dimension of the cone.                                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:     Constructed cone.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarNormalConeFree, MvarNormalConeFreeList, MvarNormalConeCopy           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarNormalConeNew                                                        M
*****************************************************************************/
MvarNormalConeStruct *MvarNormalConeNew(int Dim)
{
    MvarNormalConeStruct
	*NewCone = (MvarNormalConeStruct *)
			 	    IritMalloc(sizeof(MvarNormalConeStruct));

    NewCone -> ConeAxis = MvarVecNew(Dim);
    IRIT_ZAP_MEM((NewCone -> ConeAxis -> Vec), sizeof(CagdRType) * Dim);
    NewCone -> Attr = NULL;
    NewCone -> Pnext = NULL;
    NewCone -> AxisMinMax[0] = NewCone -> AxisMinMax[1] = 0.0;

    return NewCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Copy a multivariate normal cone structure.	                             M
*                                                                            *
* PARAMETERS:                                                                M
*   NormalCone:   Normal cone to copy.                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:     Copied normal cone.                          M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarNormalConeNew, MvarNormalConeCopyList, MvarNormalConeFree            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarNormalConeCopy                                                       M
*****************************************************************************/
MvarNormalConeStruct *MvarNormalConeCopy(const MvarNormalConeStruct 
					                         *NormalCone)
{
    int Dim = NormalCone -> ConeAxis -> Dim;
    MvarNormalConeStruct
	*NewCone = (MvarNormalConeStruct *)
			 	    IritMalloc(sizeof(MvarNormalConeStruct));

    NewCone -> ConeAxis = MvarVecNew(Dim);
    CAGD_GEN_COPY(NewCone -> ConeAxis -> Vec, NormalCone -> ConeAxis -> Vec,
		  sizeof(CagdRType) * Dim);
    NewCone -> ConeAngleCosine = NormalCone -> ConeAngleCosine;
    NewCone -> Attr = IP_ATTR_COPY_ATTRS(NormalCone -> Attr);
    NewCone -> Pnext = NULL;
    CAGD_GEN_COPY(NewCone -> AxisMinMax, NormalCone -> AxisMinMax,
		  sizeof(CagdRType) * 2);

    return NewCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Copy a list of multivariate normal cones' structures.                    M
*                                                                            *
* PARAMETERS:                                                                M
*   NormalCones:   List of multivariate normal cones to duplicate.           M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  Duplicated list of multivariate normal cones.   M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarNormalConeNew, MvarNormalConeCopy, MvarNormalConeFree	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarNormalConeCopyList, multi-variates                                   M
*****************************************************************************/
MvarNormalConeStruct *MvarNormalConeCopyList(const MvarNormalConeStruct
					                         *NormalCones)
{
    MvarNormalConeStruct *MVCTemp, *NewMVCList;

    if (NormalCones == NULL)
	return NULL;
    MVCTemp = NewMVCList = MvarNormalConeCopy(NormalCones);
    NormalCones = NormalCones -> Pnext;
    while (NormalCones != NULL) {
	MVCTemp -> Pnext = MvarNormalConeCopy(NormalCones);
	MVCTemp = MVCTemp -> Pnext;
	NormalCones = NormalCones -> Pnext;
    }
    return NewMVCList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Free a multivariate normal cone structure.                               M
*                                                                            *
* PARAMETERS:                                                                M
*   NormalCone:   Normal cone to free.                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarNormalConeNew, MvarNormalConeFreeList, MvarNormalConeCopy            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarNormalConeFree                                                       M
*****************************************************************************/
void MvarNormalConeFree(MvarNormalConeStruct *NormalCone)
{
    if (NormalCone == NULL)
	return;

    IP_ATTR_FREE_ATTRS(NormalCone -> Attr);
    MvarVecFree(NormalCone -> ConeAxis);
    IritFree(NormalCone);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Deallocates and frees a list of multi-variate cone structures.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   NormalConeList:    Multi-Variate cone list to free.                      M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarNormalConeFree, MvarNormalConeNew                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarNormalConeFreeList		                                     M
*****************************************************************************/
void MvarNormalConeFreeList(MvarNormalConeStruct *NormalConeList)
{
    MvarNormalConeStruct *NormalConeTemp;

    while (NormalConeList) {
	NormalConeTemp = NormalConeList -> Pnext;
	MvarNormalConeFree(NormalConeList);
	NormalConeList = NormalConeTemp;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Minimum spanning cone (MSC) computation of a set of vectors.	     M
* Find a central vector as the average of all given vectors and find the     M
* vector with maximal angular distance from it.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVVecs:          The set of vectors to compute their MSC.		     M
*   VecsNormalized:  TRUE if vectors are normalized, FALSE otherwise.        M
*   NumOfVecs:       Number of vectors in set MVVecs.			     M
*   MVCone:          Returns cone axis and cone cos angle of computed MSC.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:          TRUE if successful, FALSE otherwise.                       M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMMinSpanConeAvg		                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMinSpanConeAvg                                                       M
*****************************************************************************/
int MvarMinSpanConeAvg(MvarVecStruct *MVVecs,
		       int VecsNormalized,
		       int NumOfVecs,
		       MvarNormalConeStruct *MVCone)
{
    int i,
	Dim = MVVecs[0].Dim;
    IrtRType IProd,
	MinIProd = 1.0;
    MvarVecStruct *NrmlMVVecs,
	*ConeAxis = MVCone -> ConeAxis;

    if (NumOfVecs < 2) {
        MVAR_FATAL_ERROR(MVAR_ERR_MSC_TOO_FEW_PTS);
	return FALSE;
    }

    if (!VecsNormalized) {
        NrmlMVVecs = MvarVecArrayNew(NumOfVecs, Dim);
	IRIT_GEN_COPY(NrmlMVVecs, MVVecs, NumOfVecs * sizeof(MvarVecStruct));
	for (i = 0; i < NumOfVecs; i++) {
	    assert(MVVecs[i].Dim == Dim);     /* All dims must be the same. */
	    NrmlMVVecs[i].Attr = NULL;
	    NrmlMVVecs[i].Pnext = NULL;
	    MvarVecNormalize(&NrmlMVVecs[i]);
	}
    }
    else
        NrmlMVVecs = MVVecs;
    
    /* Compute the center of the cone. */
    MVAR_VEC_RESET(ConeAxis);
    for (i = 0; i < NumOfVecs; i++)
        MvarVecAdd(ConeAxis, ConeAxis, &NrmlMVVecs[i]);
    MvarVecNormalize(ConeAxis);

    /* Compute the aperture of the cone. */
    for (i = 0; i < NumOfVecs; i++) {
        IProd = MvarVecDotProd(ConeAxis, &NrmlMVVecs[i]);
	if (MinIProd > IProd)
	    MinIProd = IProd;
    }
    MVCone -> ConeAngleCosine = MinIProd;

    if (!VecsNormalized)
	MvarVecArrayFree(NrmlMVVecs, NumOfVecs);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Minimum spanning cone (MSC) computation of a set of vectors.	     M
* Algorithm is based on the Minimum Spanning Circle in Section 4.7 of        M
* "Computational Geometry, Algorithms and Applications" by M. de Berg et. al.M
*                                                                            *
* PARAMETERS:                                                                M
*   MVVecs:          The set of vectors to compute their MSC.  		     M
*   VecsNormalized:  TRUE if vectors are normalized, FALSE otherwise.        M
*   NumOfVecs:       Number of vectors in set MVVecs.			     M
*   MVCone:          Returns cone axis and cone cos angle of computed MSC.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:          TRUE if successful, FALSE otherwise.                       M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMMinSpanCone		                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMinSpanCone                                                          M
*****************************************************************************/
int MvarMinSpanCone(MvarVecStruct *MVVecs,
		    int VecsNormalized,
		    int NumOfVecs,
		    MvarNormalConeStruct *MVCone)
{
    int i, j, Found,
	Dim = MVVecs[0].Dim;
    IrtRType ConeCosAngle;
    MvarVecStruct *NrmlMVVecs, **ActiveVecs,
	*ConeAxis = MVCone -> ConeAxis;

    if (NumOfVecs < 2) {
        MVAR_VEC_COPY(ConeAxis, &MVVecs[0]);
	MVCone -> ConeAngleCosine = 1.0;
	return TRUE;
    }

    if (!VecsNormalized) {
        NrmlMVVecs = MvarVecArrayNew(NumOfVecs, Dim);
	IRIT_GEN_COPY(NrmlMVVecs, MVVecs, NumOfVecs * sizeof(MvarVecStruct));
	for (i = 0; i < NumOfVecs; i++) {
	    assert(MVVecs[i].Dim == Dim);     /* All dims must be the same. */
	    NrmlMVVecs[i].Attr = NULL;
	    NrmlMVVecs[i].Pnext = NULL;
	    MvarVecNormalize(&NrmlMVVecs[i]);
	}
    }
    else
        NrmlMVVecs = MVVecs;

    ActiveVecs = (MvarVecStruct **) IritMalloc(sizeof(MvarVecStruct *) * Dim);

    /* Compute an initial guess out of two non collinear vectors. */
    Found = FALSE;
    for (i = 0; !Found && i < NumOfVecs - 1; i++) {
        for (j = i + 1; !Found && j < NumOfVecs; j++) {
	    if (MvarVecDotProd(&NrmlMVVecs[i],
			       &NrmlMVVecs[j]) > -1.0 + IRIT_UEPS) {
	        MvarVecBlend(ConeAxis, &NrmlMVVecs[i], &NrmlMVVecs[j], 0.5);
		MvarVecNormalize(ConeAxis);
		Found = TRUE;

		/* Place i'th and j'th vectors in places 0 and 1. */
		if (i != 0) {
		    IRIT_SWAP(CagdRType *, NrmlMVVecs[i].Vec,
					   NrmlMVVecs[0].Vec);
		}
		if (j != 1) {
		    IRIT_SWAP(CagdRType *, NrmlMVVecs[j].Vec,
					   NrmlMVVecs[1].Vec);
		}
	    }
	}
    }
    if (!Found) {
        if (!VecsNormalized)
	    IritFree(NrmlMVVecs);
	IritFree(ActiveVecs);

	return FALSE;
    }

    ConeCosAngle = MvarVecDotProd(&NrmlMVVecs[0], ConeAxis);

    /* And examine the rest if inside. */
    for (i = 2; i < NumOfVecs; i++) {
	IrtRType
	    CosAng = MvarVecDotProd(&NrmlMVVecs[i], ConeAxis);

	if (CosAng < ConeCosAngle) {
	    ActiveVecs[0] = &NrmlMVVecs[i];
	    if (!MVMinSpanConeWithActiveVecs(NrmlMVVecs, i, ActiveVecs, 1,
					     MVCone)) {
	        if (!VecsNormalized)
		    MvarVecArrayFree(NrmlMVVecs, NumOfVecs);
		IritFree(ActiveVecs);

		return FALSE;
	    }
	    ConeCosAngle = MVCone -> ConeAngleCosine;
	}
    }
    MVCone -> ConeAngleCosine = ConeCosAngle;

#   ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugMVTestMSC, FALSE) {
	    for (i = j = 0; i < NumOfVecs; i++) {
	        int k;
		IrtRType
		    d = MvarVecDotProd(ConeAxis, &NrmlMVVecs[i]);

		if (IRIT_APX_EQ_EPS(d, ConeCosAngle, IRIT_UEPS))
		    j++;
	        else if (d < ConeCosAngle) {
		    printf("failed to fit a cone to vector %d = ", i);
		    for (k = 0; k < Dim; k++)
		        printf("%f ", NrmlMVVecs[i].Vec[k]);
		    printf("\n");
		}
	    }
	    if (j < 2 || j > Dim)       /* Should be supported by Dim vecs. */
	        printf("Cone is supported by %d vectors\n", j);
	}
    }
    {
	IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugMVTestMSCExhaustive, FALSE) {
	    MvarVecStruct
	        **MVActiveVecs = IritMalloc(Dim * sizeof(MvarVecStruct *));
	    MvarNormalConeStruct
		*MVConeExh = MvarNormalConeNew(Dim);

	    MVConeExh -> ConeAngleCosine = -1.0;          /* Initial guess. */
	    MvarMinSpanConeExhaustive(MVVecs, NumOfVecs, MVActiveVecs, 0, 0,
				      MVConeExh);

	    if (!IRIT_APX_EQ_EPS(MVCone -> ConeAngleCosine,
				 MVConeExh -> ConeAngleCosine, 1e-10) ||
		!IRIT_APX_EQ_EPS(MVCone -> ConeAxis -> Vec[0],
				 MVConeExh -> ConeAxis -> Vec[0], 1e-10) ||
		!IRIT_APX_EQ_EPS(MVCone -> ConeAxis -> Vec[1],
				 MVConeExh -> ConeAxis -> Vec[1], 1e-10) ||
		!IRIT_APX_EQ_EPS(MVCone -> ConeAxis -> Vec[2],
				 MVConeExh -> ConeAxis -> Vec[2], 1e-10))
	        printf("Exhaustive cone test returned a different answer.\n");

	    MvarNormalConeFree(MVConeExh);
	    IritFree(MVActiveVecs);
	}
    }

#   endif /* DEBUG */

    if (!VecsNormalized)
	MvarVecArrayFree(NrmlMVVecs, NumOfVecs);
    IritFree(ActiveVecs);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Auxiliary function of MvarMinSpanCone.                                   *
*                                                                            *
* PARAMETERS:                                                                *
*   MVVecs:          The set of vector to compute their MSC.                 *
*   NumOfVecs:       Number of vectors in set MVVecs.			     *
*   ActiveVecs:	     Active vectors that must be on the boundary of the MSC. *
*   NumOfActiveVecs: Number of vectors in ActiveVecs. Must be less than Dim. *
*   MVCone:          Holds ConeAxis and ConeCosAngle of computed MSC.        *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:	FALSE if failed, TRUE otherwise.                             *
*****************************************************************************/
static int MVMinSpanConeWithActiveVecs(MvarVecStruct *MVVecs,
				       int NumOfVecs,
				       MvarVecStruct **ActiveVecs,
				       int NumOfActiveVecs,
				       MvarNormalConeStruct *MVCone)
{
    int i, j,
	Dim = MVVecs[0].Dim;
    CagdRType
	*ConeCosAngle = &MVCone -> ConeAngleCosine;
    MvarVecStruct
	*ConeAxis = MVCone -> ConeAxis;

    /* Compute a trivial bound to first vector. */
    if (NumOfActiveVecs == 1) {
        if (MvarVecDotProd(&MVVecs[0], ActiveVecs[0]) < -1.0 + IRIT_UEPS)
	    return FALSE;
	MvarVecBlend(ConeAxis, &MVVecs[0], ActiveVecs[0], 0.5);
    }
    else if (NumOfActiveVecs == 2) {
        if (MvarVecDotProd(ActiveVecs[1], ActiveVecs[0]) < -1.0 + IRIT_UEPS)
	    return FALSE;
	MvarVecBlend(ConeAxis, ActiveVecs[1], ActiveVecs[0], 0.5);
    }
    else {
	MVHyperConeFromNPoints3(MVCone, ActiveVecs, NumOfActiveVecs);
    }

    MvarVecNormalize(ConeAxis);
    *ConeCosAngle = MvarVecDotProd(ActiveVecs[0], ConeAxis);

    /* And examine the rest if inside. */
    for (i = 0; i < NumOfVecs; i++) {
	IrtRType
	    CosAng = MvarVecDotProd(&MVVecs[i], ConeAxis);

	if (CosAng < *ConeCosAngle - IRIT_UEPS) {
	    ActiveVecs[NumOfActiveVecs] = &MVVecs[i];

	    if (NumOfActiveVecs == Dim - 1) {
	        CagdRType R1, R2;

	        /* Compute a hyper cone through Dim active vectors. */
		MVHyperConeFromNPoints(MVCone, ActiveVecs, Dim);

		/* Find a vector inside cone and reorient cone, if needed. */
		R1 = MvarVecDotProd(ConeAxis, ActiveVecs[0]);
		for (j = 0; j < i; j++) {
		    R2 = MvarVecDotProd(ConeAxis, &MVVecs[j]);
		    if (!IRIT_APX_EQ_EPS(R1, R2, IRIT_UEPS))
		        break;
		}
		if (j < i && R1 > R2) {
		    MvarVecScale(ConeAxis, -1);
		    *ConeCosAngle *= -1;
		}
	    }
	    else {
	        MVMinSpanConeWithActiveVecs(MVVecs, i, ActiveVecs,
					    NumOfActiveVecs + 1, MVCone);
	    }
	}
    }

    return TRUE;
}

#ifdef DEBUG

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Minimum spanning cone (MSC) computation of a set of vectors.	     *
* Algorithm is exhaustive from debugging purposes by examining all n-tuples. *
*                                                                            *
* PARAMETERS:                                                                *
*   MVVecs:          The set of unit vectors to compute their MSC. 	     *
*   NumOfVecs:       Number of vectors in set MVVecs.			     *
*   Level:	     Of recursion in this exhaustive search.		     *
*   First:	     Of first element in MVVecs to consider.		     *
*   MVCone:          Returns cone axis and cone cos angle of computed MSC.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   void		                                                     *
*****************************************************************************/
static void MvarMinSpanConeExhaustive(MvarVecStruct *MVVecs,
				      int NumOfVecs,
				      MvarVecStruct **MVActiveVecs,
				      int Level,
				      int First,
				      MvarNormalConeStruct *MVCone)
{
    int i,
	Dim = MVVecs[0].Dim;

    if (Level < Dim) {
        for (i = First; i < NumOfVecs; i++) {
	    MVActiveVecs[Level] = &MVVecs[i];
	    MvarMinSpanConeExhaustive(MVVecs, NumOfVecs, MVActiveVecs,
				      Level + 1, i + 1,	MVCone);
	}
    }

    if (Level > 1 && Level <= Dim) {
        int Pos = 0,
	    Neg = 0;
	CagdRType R1, R2;
        MvarNormalConeStruct
	    *MVConeTmp = MvarNormalConeNew(Dim);

	/* Construct the cone and examine its validity and optimality. */
	if (Level == Dim)
	    MVHyperConeFromNPoints(MVConeTmp, MVActiveVecs, Dim);
	else
	    MVHyperConeFromNPoints3(MVConeTmp, MVActiveVecs, Level);

	if ((R1 = MVConeTmp -> ConeAngleCosine) < MVCone -> ConeAngleCosine) {
	    MvarNormalConeFree(MVConeTmp);
	    return;			    /* We have a better cone by now. */
	}

	for (i = 0; i < NumOfVecs; i++) {
	    R2 = MvarVecDotProd(MVConeTmp -> ConeAxis, &MVVecs[i]);
	    if (!IRIT_APX_EQ_EPS(R1, R2, IRIT_UEPS)) {
	        if (R1 > R2)
		    Neg++;
		else
		    Pos++;
	    }
	    if (Pos && Neg) {
	        MvarNormalConeFree(MVConeTmp);
	        return;				   /* This cone is no good. */
	    }
	}
	if (Neg) { /* Reverse the cone direction and span of cone. */
	    MvarVecScale(MVConeTmp -> ConeAxis, -1);
	    MVConeTmp -> ConeAngleCosine = -MVConeTmp -> ConeAngleCosine;
	    if (MVConeTmp -> ConeAngleCosine < MVCone -> ConeAngleCosine) {
	        MvarNormalConeFree(MVConeTmp);
		return;			   /* We have a better cone by now. */
	    }
	}
	/* Keep this cone - best so far. */
	MVCone -> ConeAngleCosine = MVConeTmp -> ConeAngleCosine;
	MVAR_VEC_COPY(MVCone -> ConeAxis, MVConeTmp -> ConeAxis);

	MvarNormalConeFree(MVConeTmp);
    }
}

#endif /* DEBUG */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a hyper plane in R^n through n locations specified by Vecs.   M
*                                                                            *
* PARAMETERS:                                                                M
*   MVPlane:   The result is to be placed here.                              M
*   Vecs:      Input vectors, prescribing n locations in R^n.		     M
*   n: 	       Size of array Vecs.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVHyperConeFromNPoints                                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVHyperPlaneFromNPoints                                                  M
*****************************************************************************/
int MVHyperPlaneFromNPoints(MvarPlaneStruct *MVPlane,
			    MvarVecStruct * const *Vecs,
			    int n)
{
    int i,
	Dim = Vecs[0] -> Dim;
    CagdRType
        *A = (CagdRType *) IritMalloc(sizeof(CagdRType) * IRIT_SQR(Dim)),
        *x = MVPlane -> Pln,
        *b = (CagdRType *) IritMalloc(sizeof(CagdRType) * Dim);

    assert(Dim == n);

    for (i = 0; i < n; i++) {
	/* Copy the i'th vector to the matrix. */
        CAGD_GEN_COPY(&A[i * Dim], Vecs[i] -> Vec, sizeof(CagdRType) * Dim);

        b[i] = 1.0;    /* Assume free scalar coefficient can never be zero. */
    }

    /* Compute QR decomposition of matrix A. */
    if (IritQRUnderdetermined(A, NULL, NULL, Dim, Dim)) {
	return FALSE;				   /* Something went wrong. */
    }

    IritQRUnderdetermined(NULL, x, b, Dim, Dim);
    MvarPlaneNormalize(MVPlane);

    IritFree(A);
    IritFree(b);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Consturcts a hyper cone in R^n through n vectors specified by Vecs.      M
*                                                                            *
* PARAMETERS:                                                                M
*   MVCone:    The result is to be placed here.                              M
*   Vecs:      Input vectors, prescribing n locations in R^n.		     M
*   n: 	       Size of array Vecs.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVHyperPlaneFromNPoints, MVHyperConeFromNPoints2,			     M
*   MVHyperConeFromNPoints3						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVHyperConeFromNPoints                                                   M
*****************************************************************************/
int MVHyperConeFromNPoints(MvarNormalConeStruct *MVCone,
			   MvarVecStruct * const *Vecs,
			   int n)
{
    int Dim = Vecs[0] -> Dim;
    CagdRType R;
    MvarVecStruct
	*ConeAxis = MVCone -> ConeAxis;
    MvarPlaneStruct
	*MVPlane = MvarPlaneNew(Dim);

    if (!MVHyperPlaneFromNPoints(MVPlane, Vecs, Dim)) {
	MvarPlaneFree(MVPlane);
	return TRUE;
    }

    IRIT_GEN_COPY(ConeAxis -> Vec, MVPlane -> Pln, Dim * sizeof(CagdRType));
    R = MvarVecDotProd(ConeAxis, Vecs[0]);
    if (R < 0.0) {
        MvarVecScale(ConeAxis, -1);
	R = -R;
    }

    MVCone -> ConeAngleCosine = R;

    MvarPlaneFree(MVPlane);
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Consturcts a hyper cone in R^n through m (m < n) vectors specified by    M
* Vecs.									     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVCone:    The result is to be placed here.                              M
*   Vecs:      Input vectors, prescribing m locations in R^n.		     M
*   m: 	       Size of array Vecs, m < n.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVHyperPlaneFromNPoints, MVHyperConeFromNPoints, MVHyperConeFromNPoints3 M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVHyperConeFromNPoints2                                                  M
*****************************************************************************/
int MVHyperConeFromNPoints2(MvarNormalConeStruct *MVCone,
			    MvarVecStruct * const *Vecs,
			    int m)
{
    int i, j,
	n = Vecs[0] -> Dim;
    IrtGnrlMatType
	M = (IrtGnrlMatType) IritMalloc(IRIT_SQR(n + 1) * sizeof(IrtRType));
    CagdRType R,
	*b = &M[n * (n + 1)];      /* Use one allocation for both M and b. */

    assert(m < n && m > 0);

    /* Compute vectors that span the orthogonal space. */
    for (i = 0; i < m ; i++)
        CAGD_GEN_COPY(&M[i * n], Vecs[i] -> Vec, sizeof(CagdRType) * n);
    IRIT_ZAP_MEM(&M[m * n], sizeof(CagdRType) * n *(n - m));

    if (!MatGnrlOrthogonalSubspace(M, n)) {
        IritFree(M);
	return FALSE;
    }

    /* Update the last m-n rows so they will represent points on the hyper- */
    /* plane instead of vectors in plane as they do now, by adding Vec[0].  */
    for (i = m; i < n ; i++) {
	int ni = n * i;

        for (j = 0; j < n; j++) 
	    M[ni + j] += M[j];
    }

    /* Solve for the hyper-plane by computing QR decomposition of matrix M. */
    if (IritQRUnderdetermined(M, NULL, NULL, n, n)) {
	return FALSE;				   /* Something went wrong. */
    }
    for (j = 0; j < n; j++) 
        b[j] = 1.0;
    IritQRUnderdetermined(NULL, MVCone -> ConeAxis -> Vec, b, n, n);
    MvarVecNormalize(MVCone -> ConeAxis);
    R = MvarVecDotProd(MVCone -> ConeAxis, Vecs[0]);
    if (R < 0.0) {
        MvarVecScale(MVCone -> ConeAxis, -1);
	R = -R;
    }
    MVCone -> ConeAngleCosine = R;

#   ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugTestConeFromPts, FALSE) {
	    if (m == 2) {
	        MvarVecStruct
		    *VTemp = MvarVecNew(n);

		MvarVecBlend(VTemp, Vecs[1], Vecs[0], 0.5);
		MvarVecNormalize(VTemp);
		for (i = 0; i < n; i++) {
		    if (!IRIT_APX_EQ(VTemp -> Vec[i],
				MVCone -> ConeAxis -> Vec[i])) {
		        printf("Test cone case of m = 2 failed\n");
			break;
		    }
		}

		MvarVecFree(VTemp);
	    }
	}
    }
#   endif /* DEBUG */

    IritFree(M);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Consturcts a hyper cone in R^n through m (m < n) vectors specified by    M
* Vecs.	 Same functionality of MVHyperConeFromNPoints2 but more efficient,   M
* by solving for A A^T x = e, were e is [1, 1,..., 1], and having x being    M
* the linear combination of A's rows defining the cone axis.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVCone:    The result is to be placed here.                              M
*   Vecs:      Input vectors, prescribing m locations in R^n.		     M
*   m: 	       Size of array Vecs, m < n.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVHyperPlaneFromNPoints, MVHyperConeFromNPoints, MVHyperConeFromNPoints2 M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVHyperConeFromNPoints3                                                  M
*****************************************************************************/
int MVHyperConeFromNPoints3(MvarNormalConeStruct *MVCone,
			    MvarVecStruct * const *Vecs,
			    int m)
{
    int i, j,
	n = Vecs[0] -> Dim;
    IrtGnrlMatType
	M = (IrtGnrlMatType) IritMalloc(IRIT_SQR(m + 2) * sizeof(IrtRType));
    CagdRType R,
        *x = &M[m * (m + 2)],      /* Use one allocation for both M and b. */
	*e = &M[m * (m + 1)];

    assert(m < n && m > 0);

    /* Build M = A A^T, and e. */
    for (i = 0; i < m ; i++) {
	int mi = m * i;

        for (j = 0; j < m; j++) 
	    M[mi + j] = MvarVecDotProd(Vecs[i], Vecs[j]);

	e[i] = 1.0;
    }

    /* Solve for the hyper-plane by computing QR decomposition of matrix M. */
    if (IritQRUnderdetermined(M, NULL, NULL, m, m)) {
	return FALSE;				   /* Something went wrong. */
    }

    IritQRUnderdetermined(NULL, x, e, m, m);

    MVAR_VEC_RESET(MVCone -> ConeAxis);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
	    MVCone -> ConeAxis -> Vec[j] += Vecs[i] -> Vec[j] * x[i];
    }
    MvarVecNormalize(MVCone -> ConeAxis);
    R = MvarVecDotProd(MVCone -> ConeAxis, Vecs[0]);
    if (R < 0.0) {
        MvarVecScale(MVCone -> ConeAxis, -1);
	R = -R;
    }
    MVCone -> ConeAngleCosine = R;

#   ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugTestConeFromPts, FALSE) {
	    if (m == 2) {
	        MvarVecStruct
		    *VTemp = MvarVecNew(n);

		MvarVecBlend(VTemp, Vecs[1], Vecs[0], 0.5);
		MvarVecNormalize(VTemp);
		for (i = 0; i < n; i++) {
		    if (!IRIT_APX_EQ(VTemp -> Vec[i],
				MVCone -> ConeAxis -> Vec[i])) {
		        printf("Test cone case of m = 2 failed\n");
			break;
		    }
		}

		MvarVecFree(VTemp);
	    }
	}
    }
#   endif /* DEBUG */

#ifdef DEBUG_CMP_MVHYPERCONEFROMNPOINTS2
{
    MvarNormalConeStruct
        *MVC1 = MvarNormalConeNew(n);

    i = MVHyperConeFromNPoints2(MVC1, Vecs, m);
    if (i == 0 || !IRIT_APX_EQ(MVC1 -> ConeAngleCosine, MVCone -> ConeAngleCosine))
        printf("Error is cone computation 1\n");

    for (i = 0; i < n; i++)
        if (!IRIT_APX_EQ(MVC1 -> ConeAxis -> Vec[i],	MVCone -> ConeAxis -> Vec[i]))
	    printf("Error is cone computation 2\n");

    MvarNormalConeFree(MVC1);
}
#endif /* DEBUG_CMP_MVHYPERCONEFROMNPOINTS2 */

    IritFree(M);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   COmpute the maximal deviation angle between the vectors in GradPoints    *
* and the axis defined by ConeAxis.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   ConeAxis:     Axis to compute maximum angular deviation. Unit Length.    *
*   GradPoints:   The vectors to compare against ConeAxis.                   *
*   NumPoints:    Number of vectors in GradPoints.                           *
*   MaxDevIndex:  The index inGradPoints where maximal deviation occur will  *
*                 be kept here.						     *
*   MinLength, MaxLength:    Computed bounds on min/max lengths of vectors.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:    Maximal deviation, in cosine of the maximal angle, or 2.0  *
*                 if too large.						     *
*****************************************************************************/
static CagdRType MVarMVMaximalDeviation(const MvarVecStruct *ConeAxis,
					CagdRType * const *GradPoints,
					int NumPoints,
					int *MaxDevIndex,
					CagdRType *MinLength,
					CagdRType *MaxLength)
{
    int i, Index,
        Dim = ConeAxis -> Dim;
    CagdRType
        ConeAngleCosine = 2.0;
    MvarVecStruct
        *UnitVec = MvarVecNew(Dim);

    *MinLength = IRIT_INFNTY;
    *MaxLength = 0.0;
    *MaxDevIndex = -1;

    for (Index = 0; Index < NumPoints; Index++) {
        CagdRType InnerProd, VecLength;
	
	for (i = 0; i < Dim; ++i)
	    UnitVec -> Vec[i] = GradPoints[i][Index];
        VecLength = MvarVecLength(UnitVec);

	if (VecLength > IRIT_UEPS) {
	    InnerProd = MvarVecDotProd(ConeAxis, UnitVec) / VecLength;

            if (InnerProd < 0.0) {
                /* More than 90 degrees from axis of cone - abort. */
                MvarVecFree(UnitVec);
                return 2.0;
            }

	    if (ConeAngleCosine > InnerProd) {
	        *MaxDevIndex = i;
		ConeAngleCosine = InnerProd;
	    }

	    InnerProd *= VecLength;
	    if (InnerProd < *MinLength)
		*MinLength = InnerProd;
	    else if (*MaxLength < InnerProd)
		*MaxLength = InnerProd;
	}
	else
	    *MinLength = 0.0;
    }

    MvarVecFree(UnitVec);
    return ConeAngleCosine == 2.0 ? 2.0
			          : IRIT_BOUND(ConeAngleCosine, -1.0, 1.0);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the smallest principal direction of a set of normal vectors.    M
*                                                                            *
* PARAMETERS:                                                                M
*   SPDVec:   Vector to be updated with the smallest principla component.    M
*   ConeAxis: Axis of cone bounding all the normal vectors in GradPoints.    M
*             Assumed unit length.					     M
*   GradPoints:  The normal (Gradient) vectors to handle.		     M
*   TotalLength:  Number of normal vectors in GradPoints.		     M
*   Dim:      Dimension of Normal vectors.     				     N
*                                                                            *
* RETURN VALUE:                                                              M
*   void								     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarMVMaximalDeviation		                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVarSmallestPrincipalDirection                                           M
*****************************************************************************/
void MVarSmallestPrincipalDirection(MvarVecStruct *SPDVec,
				    MvarVecStruct *ConeAxis,
				    CagdRType * const *GradPoints,
				    int TotalLength,
				    int Dim)
{
    int i, j, k;
    IrtRType MinSV, Dp,
        *Sigma = (IrtRType *) IritMalloc(sizeof(IrtRType) * IRIT_SQR(Dim)),
        *V = (IrtRType *) IritMalloc(sizeof(IrtRType) * IRIT_SQR(Dim)),
        *U = (IrtRType *) IritMalloc(sizeof(IrtRType) * IRIT_SQR(Dim)),
        *S = (IrtRType *) IritMalloc(sizeof(IrtRType) * Dim);
    IrtGnrlMatType
        B = (IrtGnrlMatType) IritMalloc(sizeof(CagdRType) * Dim * TotalLength);
    MvarVecStruct
        *UnitVec = MvarVecNew(Dim),
        *MyConeAxis = MvarVecNew(Dim); 

    /* Build matrix B of size (Dim x TotalLength) of projected normals onto */
    /* the hyperplane orthogonal to ConeAxis.				    */
    for (i = 0; i< TotalLength; i++) {
	for (j = 0; j < Dim; j++)
	    UnitVec -> Vec[j] = GradPoints[j][i];
	MvarVecNormalize(UnitVec);

	IRIT_GEN_COPY(MyConeAxis -> Vec, ConeAxis -> Vec,
		      Dim * sizeof(IrtRType));

	/* Project UnitVec to hyperplane orthogonal to (My)ConeAxis: */
	MvarVecSub(UnitVec, UnitVec,
		   MvarVecScale(MyConeAxis,
				MvarVecDotProd(MyConeAxis, UnitVec)));

	for (j = 0; j < Dim; j++)
	    B[i + TotalLength * j] = UnitVec -> Vec[j];
    }

    /* Compute B^T * B, creating Sigma of size (Dim x Dim). */
    for (i = 0; i < Dim; i++) {
	for (j = 0; j < Dim; j++) {
	    IrtRType
	        Sum = 0.0;

	    for (k = 0; k < TotalLength; k++)
		Sum += B[k + TotalLength * i] * B[k + TotalLength * j];

	    Sigma[i + Dim * j] = Sum;
	}
    }

    /* Apply SVD to Sigma: */
    SvdMatrixNxN(Sigma, U, S, V, Dim);

    /* Find minimal singular value. */
    MinSV = S[0];
    j = 0;
    for (i = 1; i < Dim; i++) {
	assert(S[i] >= 0);

	if (S[i] < MinSV) {
	    for (k = 0; k < Dim; k++)
	        SPDVec -> Vec[k] = V[i + Dim * k];

	    /* One vector is expected to be the vector of ConeAxis itself  */
	    /* with a singular value of 0 as we projected the vectors to   */
	    /* the plane orthgoonal to ConeAxis.  Ignore that solution.    */
	    Dp = MvarVecDotProd(SPDVec, ConeAxis);   /* Both unit vectors. */

	    if (IRIT_ABS(Dp) < 0.5) {  /* Ignore if collinear to ConeAxis. */ 
		MinSV = S[i];
		j = i;
	    }
	}
    }

    for (i = 0; i < Dim; i++)
	SPDVec -> Vec[i] = V[i + Dim * j];

    MvarVecFree(UnitVec);
    MvarVecFree(MyConeAxis); 

    IritFree(Sigma);
    IritFree(V);
    IritFree(U);
    IritFree(S);
    IritFree(B);
}   

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a normal cone to a scalar multivariate MV.		     M
* Note the normal space of the trivariate is assumed of dimension one, and   M
* the gradient of the multivariate is assumed to span the normal space.	     M
*   If the input multivariate is not scalar it is assumed to be of point     M
* type E(1+Dim), where Dim is the dimension of the MV.  This scalar holds    M
* the gradient already in the Dim locations of the Points array in MV, in    M
* slots Points[2] to Points[Dim + 1], with Points[1] still holding the       M
* scalar field.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:   Multivariate to derive a cone bounding its normal space.           M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  A cone bounding the normal space of MV, or NULL M
*       if failed (i.e. angular span of normal space too large).             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarMVNormalCone2, MvarMVConesOverlap                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVarMVNormalCone                                                         M
*****************************************************************************/
MvarNormalConeStruct *MVarMVNormalCone(const MvarMVStruct *MV)
{
    CagdBType
	IsRational = MVAR_IS_RATIONAL_MV(MV);
    int TotalLength, MaxDevIndex;
    CagdRType * const *GradPoints;
    MvarMVGradientStruct *Grad;
    MvarNormalConeStruct *NormalCone;

    if (IsRational)			       /* No support for rationals. */
	return NULL;

    if (MVAR_NUM_OF_MV_COORD(MV) == 1) {        /* No precomputed gradient. */
        Grad = MvarMVBoundGradient(MV);
	GradPoints = &Grad -> MVGrad -> Points[1];
	TotalLength = MVAR_CTL_MESH_LENGTH(Grad -> MVGrad);
    }
    else if (MV -> Dim == MVAR_NUM_OF_MV_COORD(MV) - 1) {
        /* Gradient is embedded in Points[2] to Points[Dim + 1]. */
        Grad = NULL;
	GradPoints = &MV -> Points[2];
	TotalLength = MVAR_CTL_MESH_LENGTH(MV);
    }
    else {
        MVAR_FATAL_ERROR(MVAR_ERR_DIM_TOO_HIGH);   
	return NULL;
    }

    NormalCone = MVarMVNormalCone2(MV, GradPoints,
				   TotalLength, &MaxDevIndex);

    if (Grad != NULL)
	MvarMVFreeGradient(Grad);

    return NormalCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A second version of MVarMVNormalCone in which the control vectors are    M
* given directly.							     M
* Note the normal space of the trivariate is assumed of dimension one, and   M
* the gradient of the multivariate is assumed to span the normal space.	     M
*   If the input multivariate is not scalar it is assumed to be of point     M
* type E(1+Dim), where Dim is the dimension of the MV.  This scalar holds    M
* the gradient already in the Dim locations of the Points array in MV, in    M
* slots Points[2] to Points[Dim + 1], with Points[1] still holding the       M
* scalar field.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:           Multivariate to derive a cone bounding its normal space.   M
*   GradPoints:   Control vectors of gradient field.			     M
*   TotalLength:  Number of control vectors in gradient field.		     M
*   MaxDevIndex:  The index in GradPoints where maximal deviation occur will M
*                 be kept here.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  A cone bounding the normal space of MV, or NULL M
*       if failed (i.e. angular span of normal space too large).             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarMVNormalCone, MvarMVConesOverlap                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVarMVNormalCone2                                                        M
*****************************************************************************/
MvarNormalConeStruct *MVarMVNormalCone2(const MvarMVStruct *MV,
					CagdRType * const *GradPoints,
					int TotalLength,
					int *MaxDevIndex)
{
    int i,
	Dim = MV -> Dim,
	Index = 0;
    CagdRType
	MinLength = IRIT_INFNTY,
        MaxLength = 0.0;
    MvarVecStruct
        *UnitVec = MvarVecNew(Dim); 
    MvarNormalConeStruct
	*NormalCone = MvarNormalConeNew(Dim);

    /* Construct the central axis of the cone, ConeAxis, as an average of    */
    /* the control points of the gradient.                                   */
    for (Index = 0; Index < TotalLength; Index++) {
        for (i = 0; i < Dim; i++)
	    UnitVec -> Vec[i] = GradPoints[i][Index];

        MvarVecNormalize(UnitVec);

        MvarVecAdd(NormalCone -> ConeAxis, NormalCone -> ConeAxis, UnitVec);
    }

    MvarVecFree(UnitVec);

    if (!MvarVecNormalize(NormalCone -> ConeAxis)) { /* Zero len. avg. vec. */
	MvarNormalConeFree(NormalCone);
	return NULL;
    }

    /* Find the maximal angle, ConeAngleCosine, between ConeAxis and some   */
    /* vector in mesh.							    */
    if ((NormalCone -> ConeAngleCosine = MVarMVMaximalDeviation(
					    NormalCone -> ConeAxis,
					    GradPoints, TotalLength,
					    MaxDevIndex,
					    &MinLength, &MaxLength)) > 1.0) {
	MvarNormalConeFree(NormalCone);
	return NULL;
    }

    NormalCone -> AxisMinMax[0] = MinLength;
    NormalCone -> AxisMinMax[1] = MaxLength;

    return NormalCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a normal cone to a scalar multivariate MV.		     M
* Note the normal space of the trivariate is assumed of dimension one, and   M
* the gradient of the multivariate is assumed to span the normal space.	     M
*   If the input multivariate is not scalar it is assumed to be of point     M
* type E(1+Dim), where Dim is the dimension of the MV.  This scalar holds    M
* the gradient already in the Dim locations of the Points array in MV, in    M
* slots Points[2] to Points[Dim + 1], with Points[1] still holding the       M
* scalar field.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:         Multivariate to derive a cone bounding its normal space.     M
*   MainAxis:   Main axis (principal component) of the normal cone's         M
*               vectors distribution.  Valid only if success.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  A cone bounding the normal space of MV, or NULL M
*       if failed (i.e. angular span of normal space too large).             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarMVNormalConeMainAxis2, MVarMVNormalCone, MvarMVConesOverlap          M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVarMVNormalConeMainAxis                                                 M
*****************************************************************************/
MvarNormalConeStruct *MVarMVNormalConeMainAxis(const MvarMVStruct *MV,
					       MvarVecStruct **MainAxis)
{
    int TotalLength;
    CagdRType * const *GradPoints;
    MvarMVGradientStruct *Grad;
    MvarNormalConeStruct *NormalCone;

    if (MVAR_NUM_OF_MV_COORD(MV) == 1) {        /* No precomputed gradient. */
        Grad = MvarMVBoundGradient(MV);
	GradPoints = &Grad -> MVGrad -> Points[1];
	TotalLength = MVAR_CTL_MESH_LENGTH(Grad -> MVGrad);
    }
    else if (MV -> Dim == MVAR_NUM_OF_MV_COORD(MV) - 1) {
        /* Gradient is embedded in Points[2] to Points[Dim + 1]. */
        Grad = NULL;
	GradPoints = &MV -> Points[2];
	TotalLength = MVAR_CTL_MESH_LENGTH(MV);
    }
    else {
        MVAR_FATAL_ERROR(MVAR_ERR_DIM_TOO_HIGH);   
	return NULL;
    }

    NormalCone = MVarMVNormalConeMainAxis2(MV, GradPoints,
					   TotalLength, MainAxis);

    if (Grad != NULL)
	MvarMVFreeGradient(Grad);

    return NormalCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A second version of MVarMVNormalConeMainAxis in which the control points M
* are given directly.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:           Multivariate to derive a cone bounding its normal space.   M
*   GradPoints:   Control vectors of gradient field.			     M
*   TotalLength:  Number of control vectors in gradient field.		     M
*   MainAxis:     Main axis (principal component) of the normal cone's       M
*                 vectors distribution. Allocated and valid only if success. M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  A cone bounding the normal space of MV, or NULL M
*       if failed (i.e. angular span of normal space too large).             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarMVNormalConeMainAxis, MVarMVNormalCone, MvarMVConesOverlap           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVarMVNormalConeMainAxis2                                                M
*****************************************************************************/
MvarNormalConeStruct *MVarMVNormalConeMainAxis2(const MvarMVStruct *MV,
					        CagdRType * const *GradPoints,
						int TotalLength,
						MvarVecStruct **MainAxis)
{
    int i, MaxDevIndex,
        Dim = MV -> Dim;
    CagdRType Res;
    MvarNormalConeStruct *Cone;

    /* Compute regular normal cone: */
    if ((Cone = MVarMVNormalCone2(MV, GradPoints,
				  TotalLength, &MaxDevIndex)) == NULL) {
        return NULL;
    }

    *MainAxis = MvarVecNew(Dim);
    if (MaxDevIndex >=0) {		      /* Use it to derive MainAxis. */
	for (i = 0; i < Dim; i++)
	    (*MainAxis) -> Vec[i] = GradPoints[i][MaxDevIndex] -
	                                           Cone -> ConeAxis -> Vec[i];
	Res = 0.0;
    }
    else {  /* Fit a line through all the cloud of grad vectors and use it. */
        int Index;
        MvarVecStruct
	    *LinePos = MvarVecNew(Dim),
            *VecList = NULL;

	/* Convert grad vectors to a list so we can fit the line to them. */
	for (Index = 0; Index < TotalLength; Index++) {
	    MvarVecStruct
	        *UnitVec = MvarVecNew(Dim);

	    for (i = 0; i < Dim; i++)
	        UnitVec -> Vec[i] = GradPoints[i][Index];
	    MvarVecNormalize(UnitVec);
	    IRIT_LIST_PUSH(UnitVec, VecList);
	}

	Res = MvarLineFitToPts(VecList, *MainAxis, LinePos);
	MvarVecFreeList(VecList);
	MvarVecFree(LinePos);
    }
    MvarVecNormalize(*MainAxis);

    if (Res == IRIT_INFNTY) {
        MvarNormalConeFree(Cone);
	MvarVecFree(*MainAxis);
        return NULL;
    }
    else
        return Cone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes a 2cones bound to the normal field of multivariate MV.	     M
*  The 2cones bounds the normal field in the common intersection space.      M
*  The 2cones are computed using the regular normal cone by expanding in the M
* direction orthogonal to the cone axis and its main principal component.    M
*  The expansion is done an amount that is equal to regular cone radius      M
* times ExpandingFactor.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV:              To compute the normal 2cones for.	                     M
*   ExpandingFactor: Factor to expand placement of 2cones axes locations.    M
*   NumOfZeroMVs:    Number of zero type MVs in the problem we solve.	     M
*   Cone1, Cone2:    The two cones to compute or ConeAngle = M_PI if error.  M
*                      Can be NULL in which case no 2cones are computed -    M
*                    only the regular cone is computed.
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:   Regular normal cone if successful,	     M
*			      NULL otherwise.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbNormalConeForSrf, MvarNormalConeOverlap				     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVNormal2Cones, normals, normal bound		                     M
*****************************************************************************/
MvarNormalConeStruct *MvarMVNormal2Cones(const MvarMVStruct *MV,
					 CagdRType ExpandingFactor,
					 int NumOfZeroMVs,
					 MvarNormalConeStruct **Cone1,
					 MvarNormalConeStruct **Cone2)
{
    CagdBType
	IsRational = MVAR_IS_RATIONAL_MV(MV);
    int TotalLength, MaxDevIndex1, MaxDevIndex2;
    CagdRType MinLength, MaxLength, d, * const *GradPoints;
    MvarVecStruct *MainAxis, *Dir;
    MvarNormalConeStruct *Cone;
    MvarMVGradientStruct *Grad;

    if (IsRational) {
	MVAR_FATAL_ERROR(MVAR_ERR_RATIONAL_NO_SUPPORT);
	return NULL;
    }

    if (MVAR_NUM_OF_MV_COORD(MV) == 1) {       /* No precomputed gradient. */
        Grad = MvarMVBoundGradient(MV);
	GradPoints = &Grad -> MVGrad -> Points[1];
	TotalLength = MVAR_CTL_MESH_LENGTH(Grad -> MVGrad);
    }
    else if (MV -> Dim == MVAR_NUM_OF_MV_COORD(MV) - 1) {
        /* Gradient is embedded in Points[2] to Points[Dim + 1]. */
        Grad = NULL;
	GradPoints = &MV -> Points[2];
	TotalLength = MVAR_CTL_MESH_LENGTH(MV);
    }
    else {
        MVAR_FATAL_ERROR(MVAR_ERR_DIM_TOO_HIGH); 
	if (Cone1 != NULL && Cone2 != NULL)
	    *Cone1 = *Cone2 = NULL;
	return NULL;
    }

    if ((Cone = MVarMVNormalConeMainAxis2(MV, GradPoints, TotalLength,
					  &MainAxis)) == NULL) {
        /* Failed to compute the regular cone? */
        if (Grad != NULL)
	    MvarMVFreeGradient(Grad);
	if (Cone1 != NULL && Cone2 != NULL)
	    *Cone1 = *Cone2 = NULL;
        return NULL;
    }

    /* Cases were 2cones cannot be computed. */
    d = MvarVecDotProd(Cone -> ConeAxis, MainAxis);
    if (Cone1 == NULL ||
	Cone2 == NULL ||
	MV -> Dim < 2 ||
	MV -> Dim != NumOfZeroMVs ||
	Cone -> ConeAngleCosine > MVAR_2CONES_MAX_CONE_ANGLE ||
	IRIT_ABS(d) > 1.0 - IRIT_UEPS) {
        if (Grad != NULL)
	    MvarMVFreeGradient(Grad);
	MvarVecFree(MainAxis);
	if (Cone1 != NULL && Cone2 != NULL)
	    *Cone1 = *Cone2 = NULL;
	return Cone;
    }

    /* Find a direction orthogonal to both vectors. */
    Dir = MvarVecNew(MV -> Dim);

    MvarVecOrthogonal2(Dir, Cone -> ConeAxis, MainAxis);
    MvarVecFree(MainAxis);
    MvarVecNormalize(Dir);
    ExpandingFactor *= tan(acos(Cone -> ConeAngleCosine));
    MvarVecScale(Dir, ExpandingFactor);

    *Cone1 = MvarNormalConeNew(MV -> Dim);
    *Cone2 = MvarNormalConeNew(MV -> Dim);

    MvarVecAdd((*Cone1) -> ConeAxis, Cone -> ConeAxis, Dir);
    MvarVecSub((*Cone2) -> ConeAxis, Cone -> ConeAxis, Dir);
    MvarVecNormalize((*Cone1) -> ConeAxis);
    MvarVecNormalize((*Cone2) -> ConeAxis);

    MvarVecFree(Dir);

    /* Find the maximal angle, ConeAngleCosine, between ConeAxis and some   */
    /* vector in mesh.							    */
    (*Cone1) -> ConeAngleCosine = MVarMVMaximalDeviation((*Cone1) -> ConeAxis,
						      GradPoints, TotalLength,
						      &MaxDevIndex1,
						      &MinLength, &MaxLength);
    (*Cone2) -> ConeAngleCosine = MVarMVMaximalDeviation((*Cone2) -> ConeAxis,
						      GradPoints, TotalLength,
						      &MaxDevIndex2,
						      &MinLength, &MaxLength);

    if (Grad != NULL)
	MvarMVFreeGradient(Grad);

    return Cone;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Update, in place, bCopy with the actual values to solve for the next     *
* vertex to be examined if inside the unit sphere.                           *
*                                                                            *
* PARAMETERS:                                                                *
*   CurrentPower:  The current iteration of 2^n possibilities.               *
*   Dim:           Dimension of b and bCopy.				     *
*   b:             The current configuration of the planes to solve for.     *
*   bCopy:         Actual configuration to be updated in place, based on b   *
*                  and CurrentPower values.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void MvarConesAssembleB(int CurrentPower,
			       int Dim,
			       const CagdRType *b,
			       CagdRType *bCopy)
{
    int j,
        k = CurrentPower;

    for (j = 0; j < Dim; ++j) {
        bCopy[j] = k & 1 ? b[j] : -b[j];
	k >>= 1;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the tangency anti-cones of the set of normal cones,             M
* and returns whether they overlap or not.                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   ConesList: Cones in a list.                                              M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:  TRUE if overlap, FALSE if not.                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarConesOverlapAux, normals, normal bound		                     M
*****************************************************************************/
CagdBType MvarConesOverlapAux(const MvarNormalConeStruct *ConesList)
{
    IRIT_STATIC_DATA int
	AllocDim = 0;
    IRIT_STATIC_DATA CagdRType
	*A = NULL,
	*x = NULL,
	*b = NULL,
	*bCopy = NULL;
    int i = 0,
	Dim = ConesList -> ConeAxis -> Dim,
	PowerTwoDim = 1 << (Dim - 1);
    const MvarNormalConeStruct
        *Cone = ConesList;

    /* In order to save alloc/free run-time we use static variables. */
    if (AllocDim < Dim) {
        if (AllocDim > 0) {
            IritFree(A);
            IritFree(x);
            IritFree(b);
            IritFree(bCopy);
        }

	AllocDim = Dim * 2;
        A = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim * AllocDim);
        x = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
        b = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
        bCopy = (CagdRType *) IritMalloc(sizeof(CagdRType) * AllocDim);
    }

    for (i = 0; Cone; i++, Cone = Cone -> Pnext) {
	/* Add plane orthogonal to cone axis to matrix A. */
	CAGD_GEN_COPY(&A[i * Dim], Cone -> ConeAxis -> Vec,
		      sizeof(CagdRType) * Dim);

	/* Take the anti-cone = cos(90 - angle) = sqrt(1-cos^2). */
	b[i] = sqrt(1.0 - IRIT_SQR(Cone -> ConeAngleCosine));
    }

    /* We now have the matrix A containing planes orthogonal to cone's axes  */
    /* and the vector b of the expected solutions.                           */

    /* Compute QR decomposition of matrix A. */
    if (IritQRUnderdetermined(A, NULL, NULL, Dim, Dim)) {
	return TRUE; /* Something went wrong - return cones are overlapping. */
    }

    /* Loop over 2^(d-1) combinations of b vector (000 -> ---, 111 -> +++).  */
    /* If Qx=b returns a point that is out of unit hyper-sphere return TRUE  */
    /* meaning the cones overlap.					     */

    for (i = 0; i < PowerTwoDim; i++) {
        /* Construct relevant copy of b (+/- of b[j] defined by binary rep). */
	MvarConesAssembleB(i, Dim, b, bCopy);

        IritQRUnderdetermined(NULL, x, bCopy, Dim, Dim);

	if (MvarVecSqrLength2(x, Dim) >= 1.0) {
	    return TRUE;
	}
    }

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the tangency anti-cones of the set of multivariate constraints, M
* and returns whether they overlap or not.                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:            Multivariates to derive their tangency anti-cones.       M
*   NumOfZeroMVs:   Size of the vector MVs.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:  TRUE if overlap, FALSE if not.                               M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarMVNormalCone, MvarConesOverlapAux                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVConesOverlap                                                       M
*****************************************************************************/
CagdBType MvarMVConesOverlap(MvarMVStruct **MVs, int NumOfZeroMVs)
{
    int RetVal, i;
    MvarNormalConeStruct FirstCone,	          /* First cone is a dummy. */
        *Cones = &FirstCone,				
        *ConesStart = Cones;

    assert(NumOfZeroMVs > 0);

    for (i = 0; i < NumOfZeroMVs; i++) {
	if ((Cones -> Pnext = MvarMVNormal2Cones(MVs[i],
						 MVAR_2CONES_EXPAND_FACTOR,
						 NumOfZeroMVs, NULL, NULL))
	                                                           == NULL) {
	    /* Angular span is too large - return cones overlap. */
	    MvarNormalConeFreeList(ConesStart -> Pnext);
	    return TRUE;
	}

	Cones = Cones -> Pnext;
    }

    RetVal = MvarConesOverlapAux(ConesStart -> Pnext);

#ifdef MV_CONES_TEST_OVERLAPS
    MVConesNoOverlaps += RetVal ? 1 : 0;
    MVConesOverlapTests++;
#endif /* MV_CONES_TEST_OVERLAPS */

    MvarNormalConeFreeList(ConesStart -> Pnext);

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a normal cone of a sum.                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   ConeF: Normal cone of the the first summand.                             M
*   ConeG: Normal cone of the second summand.                                M
*   Dim:   Dimensions.                                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  The resulting normal cone.                      M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarExprTreeNormalCone, HyperplaneOrthoSystem                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeNormalConeSum                                                M
*****************************************************************************/
MvarNormalConeStruct *MvarExprTreeNormalConeSum(
					   const MvarNormalConeStruct *ConeF,
				  	   const MvarNormalConeStruct *ConeG,
					   int Dim)
{
    const int
        VecListLength = 16 * (Dim - 1) * (Dim - 1);
    int i, j, k, l,
        ZeroVecs = 0;
    MvarVecStruct
        *VecList = MvarVecArrayNew(VecListLength, Dim),
        *OrthoF = HyperplaneOrthoSystem(ConeF -> ConeAxis),
        *OrthoG = HyperplaneOrthoSystem(ConeG -> ConeAxis),
        *u = MvarVecNew(Dim),
        *v = MvarVecNew(Dim),
        *TmpVec = MvarVecNew(Dim),
	*CurVec = VecList;
    MvarNormalConeStruct
        *Result = MvarNormalConeNew(Dim);
    CagdRType T,
	*VecSizes = IritMalloc(sizeof(CagdRType) * VecListLength),
        *VecSize = VecSizes,
        MinAxis = IRIT_INFNTY,
        MaxAxis = 0.0;

    assert(VecList && OrthoF && OrthoG && u && v && Result && VecSizes);
    IRIT_ZAP_MEM(u -> Vec, sizeof(CagdRType) * u -> Dim);
    IRIT_ZAP_MEM(v -> Vec, sizeof(CagdRType) * v -> Dim);
    IRIT_ZAP_MEM(VecSizes, sizeof(CagdRType) * VecListLength);

    /* Align cones. */
    T = MvarVecDotProd(ConeF -> ConeAxis, ConeG -> ConeAxis);
    if (T < 0)
        MvarVecScale(ConeF -> ConeAxis, -1);
    T = 1.0 / ConeF -> ConeAngleCosine;
    T = sqrt((IRIT_SQR(T) - 1.0) * (Dim - 1));

    for (i = 1; i < Dim; i++)
	MvarVecScale(OrthoF + i, T);
    T = 1.0 / ConeG -> ConeAngleCosine;
    T = sqrt((IRIT_SQR(T) - 1.0) * (Dim - 1));

    for (i = 1; i < Dim; i++)
	MvarVecScale(OrthoG + i, T);

    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
	    for (k = 1; k < Dim; k++) {
		for (l = 1; l < Dim; l++) {
		    /* ++ */
		    MvarVecAdd(u, ConeF -> ConeAxis, OrthoF + k);
		    MvarVecAdd(v, ConeG -> ConeAxis, OrthoG + l);
		    MvarVecAdd(CurVec,
			       MvarVecScale(u, ConeF -> AxisMinMax[i]),
			       MvarVecScale(v, ConeG -> AxisMinMax[j]));
		    ADVANCE_VECTOR_COUNT(CurVec, VecSize, ZeroVecs);
		    IRIT_GEN_COPY(TmpVec -> Vec, u -> Vec,
				  sizeof(CagdRType) * Dim);

		    /* -+ */
		    MvarVecSub(u, ConeF -> ConeAxis, OrthoF + k);
		    MvarVecAdd(CurVec,
			MvarVecScale(u, ConeF -> AxisMinMax[i]), v);
		    ADVANCE_VECTOR_COUNT(CurVec, VecSize, ZeroVecs);

		    /* -- */
		    MvarVecSub(v, ConeG -> ConeAxis, OrthoG + l);
		    MvarVecAdd(CurVec,
			u, MvarVecScale(v, ConeG -> AxisMinMax[j]));
		    ADVANCE_VECTOR_COUNT(CurVec, VecSize, ZeroVecs);

		    /* +- */
		    MvarVecAdd(CurVec, TmpVec, v);
		    ADVANCE_VECTOR_COUNT(CurVec, VecSize, ZeroVecs);
		}
	    }
	}
    }

    if (!MvarMinSpanCone(VecList, TRUE, VecListLength - ZeroVecs, Result)) {
	MvarNormalConeFree(Result);
	Result = NULL;
    }
    else {
	if (ZeroVecs)
	    MinAxis = 0.0;
	for (i = 0; i < VecListLength-ZeroVecs; i++) {
	    T = MvarVecDotProd(VecList+i, Result -> ConeAxis) * VecSizes[i];
	    if (T < MinAxis)
		MinAxis = T;
	    else if (MaxAxis < T)
		MaxAxis = T;
	}

	Result -> AxisMinMax[0] = MinAxis;
	Result -> AxisMinMax[1] = MaxAxis;
    }

    MvarVecArrayFree(VecList, VecListLength);
    MvarVecArrayFree(OrthoF, Dim);
    MvarVecArrayFree(OrthoG, Dim);
    MvarVecFree(u);
    MvarVecFree(v);
    MvarVecFree(TmpVec);
    IritFree(VecSizes);

    return Result;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Constructs all spanning vectors of the hyperplane given by the normal.   *
*   Basically implements modified Graham Schmidt process.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   v: Normal of the hyperplane.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarVecStruct *:  Vectors of the spanning space.                         *
*                                                                            *
* SEE ALSO:                                                                  *
*   MvarExprTreeNormalConeSum, MatGnrlOrthogonalSubspace                     *
*****************************************************************************/
static MvarVecStruct *HyperplaneOrthoSystem(const MvarVecStruct *v)
{
    int i, j,
        Dim = v -> Dim,
        Saved = 0;
    MvarVecStruct
        *Basis = MvarVecArrayNew(Dim, Dim),
	*TmpVec = MvarVecNew(Dim + 1),
        *TmpVec1 = MvarVecNew(Dim);
    CagdRType VSize;

    assert(Basis);

    /* TmpVec is a bit larger, just to prevent memory overflows.           */
    /* The last element isn't really used. 				   */
    IRIT_ZAP_MEM(TmpVec -> Vec, sizeof(CagdRType) * (Dim + 1));
    TmpVec -> Vec[0] = 1.0;
    IRIT_GEN_COPY(Basis[0].Vec, v -> Vec, sizeof(CagdRType) * Dim);
    MvarVecNormalize(Basis);

    for (j = 1; j < Dim; j++) {
	IRIT_GEN_COPY(Basis[j].Vec, TmpVec -> Vec, sizeof(CagdRType) * Dim);
	for (i = 0;i < j; i++) {
	    IRIT_GEN_COPY(TmpVec1 -> Vec, Basis[i].Vec,
			  sizeof(CagdRType) * Dim);
	    MvarVecScale(TmpVec1, -MvarVecDotProd(Basis+i, Basis+j));
	    MvarVecAdd(Basis+j, Basis+j, TmpVec1);
	}
	VSize = MvarVecLength(Basis+j);

	TmpVec -> Vec[j -1 + Saved] = 0.0;
	TmpVec -> Vec[j + Saved] = 1.0;
	if (IRIT_FABS(VSize) < IRIT_UEPS) {
	    assert(!Saved);	/* In case something went really wrong. */
	    j--;
	    Saved = 1;
	}
	else
	    MvarVecScale(Basis+j, 1.0 / VSize);
    }

    MvarVecFree(TmpVec);
    MvarVecFree(TmpVec1);

    return Basis;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a normal cone of a difference.                                M
*                                                                            *
* PARAMETERS:                                                                M
*   ConeF: Normal cone of the minuend.                                       M
*   ConeG: Normal cone of the subtrahend.                                    M
*   Dim:   Dimensions.                                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  The resulting normal cone.                      M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarExprTreeNormalCone, MvarExprTreeNormalConeSum                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeNormalConeSub                                                M
*****************************************************************************/
MvarNormalConeStruct *MvarExprTreeNormalConeSub(
					const MvarNormalConeStruct *ConeF,
					const MvarNormalConeStruct *ConeG,
					int Dim)
{
    MvarNormalConeStruct *RetCone,
        *TmpCone = MvarNormalConeCopy(ConeG);
    int i;

    for (i = 0; i < Dim; i++)
	TmpCone -> ConeAxis -> Vec[i] = -TmpCone -> ConeAxis -> Vec[i];
    RetCone = MvarExprTreeNormalConeSum(ConeF, TmpCone, Dim);
    MvarNormalConeFree(TmpCone);
    return RetCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a normal cone of a multiplication.                            M
*                                                                            *
* PARAMETERS:                                                                M
*   ConeF: Normal cone of the first term.                                    M
*   ConeG: Normal cone of the second term.                                   M
*   BBoxF: Bounding box of the first term.                                   M
*   BBoxG: Bounding box of the second term.                                  M
*   Dim:   Dimensions.                                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  The resulting normal cone.                      M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarExprTreeNormalCone, MvarExprTreeNormalConeSum                        M
*   MvarExprTreeNormalConeScale                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeNormalConeMul                                                M
*****************************************************************************/
MvarNormalConeStruct *MvarExprTreeNormalConeMul(
					const MvarNormalConeStruct *ConeF,
					const MvarNormalConeStruct *ConeG,
					const MvarBBoxStruct *BBoxF,
					const MvarBBoxStruct *BBoxG,
					int Dim)
{
    MvarNormalConeStruct *RetCone,
        *Cone1 = NULL,
        *Cone2 = NULL;

    if (!(Cone1 = MvarExprTreeNormalConeScale(ConeF, BBoxG, Dim)))
	return NULL;

    if (!(Cone2 = MvarExprTreeNormalConeScale(ConeG, BBoxF, Dim))) {
	MvarNormalConeFree(Cone1);
	return NULL;
    }

    RetCone = MvarExprTreeNormalConeSum(Cone1, Cone2, Dim);

    MvarNormalConeFree(Cone1);
    MvarNormalConeFree(Cone2);

    return RetCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Scales normal cone (for composition and multiplication).                 M
*                                                                            *
* PARAMETERS:                                                                M
*   ConeF:       Normal cone of the function.                                M
*   BBoxGPrime:  Bounding box of the scale factor.                           M
*   Dim:         Dimensions.                                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  New scaled cone.		                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarExprTreeNormalCone, MvarExprTreeNormalConeMul                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeNormalConeScale                                              M
*****************************************************************************/
MvarNormalConeStruct *MvarExprTreeNormalConeScale(
					  const MvarNormalConeStruct *ConeF,
					  const MvarBBoxStruct *BBoxGPrime,
					  int Dim)
{
    MvarNormalConeStruct *ResCone;
    CagdRType
        Min = BBoxGPrime -> Min[0],
        Max = BBoxGPrime -> Max[0];

    assert(BBoxGPrime -> Dim == 1);

    if (Min * Max < -IRIT_EPS)
	return NULL;

    ResCone = MvarNormalConeCopy(ConeF);
    if (Min < 0.0) {
	ResCone -> AxisMinMax[0] *= -Max;
	ResCone -> AxisMinMax[1] *= -Min;
    }
    else {
	ResCone -> AxisMinMax[0] *= Min;
	ResCone -> AxisMinMax[1] *= Max;
    }

    return ResCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a normal cone to a scalar multivariate expression tree.	     M
*   If the input multivariate is not scalar it is assumed to be of point     M
* type E(1+Dim), where Dim is the dimension of the MV.  This scalar holds    M
* the gradient already in the Dim locations of the Points array in MV, in    M
* slots Points[2] to Points[Dim + 1], with Points[1] still holding the       M
* scalar field.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   Eqn:   Multivariate to derive a cone bounding its normal space.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarNormalConeStruct *:  The resulting normal cone.                      M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarExprTreeConesOverlap                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   MVarExprTreeNormalCone                                                   M
*****************************************************************************/
MvarNormalConeStruct *MVarExprTreeNormalCone(MvarExprTreeStruct *Eqn)
{
    int Dim = Eqn -> Dim;
    MvarNormalConeStruct
        *MVCone1 = NULL,
        *MVCone2 = NULL;
    const MvarBBoxStruct
        *BBox1 = NULL,
        *BBox2 = NULL;
    MvarBBoxStruct CompBBox;

    /* Prepare the cones */
    if (Eqn -> MVBCone) {
	MvarNormalConeFree(Eqn -> MVBCone);
	Eqn -> MVBCone = NULL;
    }
    switch (Eqn -> NodeType) {
	case MVAR_ET_NODE_LEAF:
	    return (Eqn -> MVBCone = MVarMVNormalCone(Eqn -> MV));
	case MVAR_ET_NODE_MULT:
	case MVAR_ET_NODE_SUB:
	case MVAR_ET_NODE_ADD:
	    if ((MVCone2 = MVarExprTreeNormalCone(Eqn -> Right)) == NULL)
		return NULL;
	    if (MVCone2 -> ConeAngleCosine > 1.0)
	        MVCone2 -> ConeAngleCosine = 1.0;
	case MVAR_ET_NODE_EXP:
	case MVAR_ET_NODE_LOG:
	case MVAR_ET_NODE_COS:
	case MVAR_ET_NODE_SQRT:
	case MVAR_ET_NODE_SQR:
	case MVAR_ET_NODE_NPOW:
	case MVAR_ET_NODE_RECIP:
	    if ((MVCone1 = MVarExprTreeNormalCone(Eqn -> Left)) == NULL)
		return NULL;
	    if (MVCone1 -> ConeAngleCosine > 1.0)
	        MVCone1 -> ConeAngleCosine = 1.0;
	    break;
	default:
	    /* DOT_PROD/CROSS_PROD are not supported here. */
	    return NULL;
    }

    /* Compute the resulting cone */
    switch(Eqn -> NodeType) {
	case MVAR_ET_NODE_ADD:
	    Eqn -> MVBCone = MvarExprTreeNormalConeSum(MVCone1, MVCone2, Dim);
	    break;
	case MVAR_ET_NODE_SUB:
	    Eqn -> MVBCone = MvarExprTreeNormalConeSub(MVCone1, MVCone2, Dim);
	    break;
	case MVAR_ET_NODE_MULT:
	    BBox1 = MvarExprTreeBBox(Eqn -> Left);
	    BBox2 = MvarExprTreeBBox(Eqn -> Right);
	    Eqn -> MVBCone = MvarExprTreeNormalConeMul(MVCone1, MVCone2, BBox1,
						       BBox2, Dim);
	    break;
	    /* Basically all scalar compositions follow and should be same */
	case MVAR_ET_NODE_EXP:
	case MVAR_ET_NODE_LOG:
	case MVAR_ET_NODE_COS:
	case MVAR_ET_NODE_SQRT:
	case MVAR_ET_NODE_SQR:
	case MVAR_ET_NODE_NPOW:
	case MVAR_ET_NODE_RECIP:
	    Eqn -> MVBCone = MvarExprTreeNormalConeScale(MVCone1,
		       MvarExprTreeCompositionDerivBBox(Eqn, &CompBBox), Dim);
	    break;
	default:
	    assert(0);
	    return FALSE;
    }

    if (Eqn -> MVBCone != NULL && Eqn -> MVBCone -> ConeAngleCosine > 1.0)
        Eqn -> MVBCone -> ConeAngleCosine = 1.0;

    return Eqn -> MVBCone;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the tangency anti-cones of the set of multivariate constraints, M
* and returns whether they overlap or not.                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   Eqns:         The MVETs constraints formated into equations with         *
*                 common expressions.					     *
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdBType:  TRUE if overlap, FALSE if not.                               M
*                                                                            *
* SEE ALSO:                                                                  M
*   MVarExprTreeNormalCone, MvarConesOverlapAux                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarExprTreeConesOverlap                                                 M
*****************************************************************************/
CagdBType MvarExprTreeConesOverlap(MvarExprTreeEqnsStruct *Eqns)
{
    CagdBType Overlap;
    int i = 0,
        NumEqns = Eqns -> NumEqns;
    MvarNormalConeStruct ConesStart,         /* Dummy to have initial setup. */
	*Cones = &ConesStart;

    /* Make a linked list from all the cones to need to check for overlap. */
    for (i = 0; i < NumEqns; i++) {
	if ((Cones -> Pnext = MVarExprTreeNormalCone(Eqns -> Eqns[i])) == NULL) {
	    /* Angular span is too large - return cones overlap. */
	    break;
	}
	Cones = Cones -> Pnext;
    }

    Overlap = i < NumEqns ? TRUE : MvarConesOverlapAux(ConesStart.Pnext);

    while ((Cones = ConesStart.Pnext) != NULL) {   /* break the linked list. */
	ConesStart.Pnext = Cones -> Pnext;
	Cones -> Pnext = NULL;
    }

    return Overlap;
}
