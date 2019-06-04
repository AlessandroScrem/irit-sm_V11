/******************************************************************************
* MvCones.c - Tools to construct and intersect MV (anti-)cones & vectors.     *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Iddo Haniel and Gershon Elber, May 2005.			      *
******************************************************************************/

#include "mvar_loc.h"

#define MVAR_ORTHO2_VALID_VAL		0.99
#define MVAR_VEC_NUMER_TOL		1e-10
#define MVAR_VEC_MAX_ORTHO_ITERS	10

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Add two multivariate vectors.                                            M
*                                                                            *
* PARAMETERS:                                                                M
*   VRes:    Result.  Can be one of V1 or V2.                                M
*   V1, V2:  Two input vectors.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecDotProd, MvarVecSqrLength, MvarVecLength, MvarVecScale            M
*   MvarVecNormalize, MvarVecBlend                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecAdd                                                               M
*****************************************************************************/
void MvarVecAdd(MvarVecStruct *VRes,
		const MvarVecStruct *V1,
		const MvarVecStruct *V2)
{
    int i;

    assert(V1 -> Dim == V2 -> Dim && VRes -> Dim == V2 -> Dim);

    for (i = 0; i < V1 -> Dim; i++)
	VRes -> Vec[i] = V1 -> Vec[i] + V2 -> Vec[i];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Add two multivariate vectors, second with scale: VRes = V1 + V2*Scale2.  M
*                                                                            *
* PARAMETERS:                                                                M
*   VRes:    Result.  Can be one of V1 or V2.                                M
*   V1, V2:  Two input vectors.						     M
*   Scale2:  Scaling factor of V2.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecDotProd, MvarVecSqrLength, MvarVecLength, MvarVecScale            M
*   MvarVecNormalize, MvarVecBlend                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecAddScale                                                          M
*****************************************************************************/
void MvarVecAddScale(MvarVecStruct *VRes,
		     const MvarVecStruct *V1,
		     const MvarVecStruct *V2,
		     CagdRType Scale2)
{
    int i;

    assert(V1 -> Dim == V2 -> Dim && VRes -> Dim == V2 -> Dim);

    for (i = 0; i < V1 -> Dim; i++)
	VRes -> Vec[i] = V1 -> Vec[i] + V2 -> Vec[i] * Scale2;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Subtract two multivariate vectors.                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   VRes:    Result.  Can be one of V1 or V2.                                M
*   V1, V2:  Two input vectors.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecDotProd, MvarVecSqrLength, MvarVecLength, MvarVecScale            M
*   MvarVecNormalize, MvarVecBlend, MvarVecAdd                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecSub                                                               M
*****************************************************************************/
void MvarVecSub(MvarVecStruct *VRes,
		const MvarVecStruct *V1,
		const MvarVecStruct *V2)
{
    int i;

    assert(V1 -> Dim == V2 -> Dim && VRes -> Dim == V2 -> Dim);

    for (i = 0; i < V1 -> Dim; i++)
	VRes -> Vec[i] = V1 -> Vec[i] - V2 -> Vec[i];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the dot product of two multivariate vectors.                     M
*                                                                            *
* PARAMETERS:                                                                M
*   V1, V2:  Two input vectors.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:  The dot product.                                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecSqrLength, MvarVecLength, MvarVecScale                M
*   MvarVecNormalize, MvarVecBlend                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecDotProd                                                           M
*****************************************************************************/
CagdRType MvarVecDotProd(const MvarVecStruct *V1, const MvarVecStruct *V2)
{
    int i;
    CagdRType
	DotProd = 0.0;

    assert(V1 -> Dim == V2 -> Dim);

    for (i = 0; i < V1 -> Dim; i++)
        DotProd += V1 -> Vec[i] * V2 -> Vec[i];

    return DotProd;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the length squared of a multivariate vector.                     M
*                                                                            *
* PARAMETERS:                                                                M
*   V:   Vector to compute its length.                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:   Computed length squared.                                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecDotProd, MvarVecLength, MvarVecScale                  M
*   MvarVecNormalize, MvarVecBlend, MvarVecSqrLength2                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecSqrLength                                                         M
*****************************************************************************/
CagdRType MvarVecSqrLength(const MvarVecStruct *V)
{
    return MvarVecDotProd(V, V);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the length squared of a multivariate vector.                     M
*                                                                            *
* PARAMETERS:                                                                M
*   v:    Vector to compute its length.                                      M
*   Dim:  Lentgth of vector v.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:   Computed length squared.                                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecDotProd, MvarVecLength, MvarVecScale                  M
*   MvarVecNormalize, MvarVecBlend, MvarVecSqrLength                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecSqrLength2                                                        M
*****************************************************************************/
CagdRType MvarVecSqrLength2(const CagdRType *v, int Dim)
{
    int j;
    CagdRType
        SumSqr = 0.0;

    for (j = 0; j < Dim; j++, v++)
        SumSqr += IRIT_SQR(*v);

    return SumSqr;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the length of a multivariate vector.      		             M
*                                                                            *
* PARAMETERS:                                                                M
*   V:   Vector to compute its length.                                       M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdRType:   Computed length. 	                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecDotProd, MvarVecSqrLength, MvarVecScale               M
*   MvarVecNormalize, MvarVecBlend                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecLength                                                            M
*****************************************************************************/
CagdRType MvarVecLength(const MvarVecStruct *V)
{
    return sqrt(MvarVecSqrLength(V));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Scale a given multivariate vector V by a scaling factor ScaleFactor.     M
*                                                                            *
* PARAMETERS:                                                                M
*   V:             Vector to scale, in place.				     M
*   ScaleFactor:   Scaling factor to use.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarVecStruct *:  The input vector, scaled.                              M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecDotProd, MvarVecSqrLength, MvarVecLength,             M
*   MvarVecNormalize, MvarVecBlend                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecScale                                                             M
*****************************************************************************/
MvarVecStruct *MvarVecScale(MvarVecStruct *V, CagdRType ScaleFactor)
{
    int i;

    for (i = 0; i < V -> Dim; i++)
        V -> Vec[i] *= ScaleFactor;

    return V;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the blend of the to given multivariate vectors as		     M
*	V1 * t + V2 * (1-t).				                     V
*                                                                            *
* PARAMETERS:                                                                M
*   VRes:    Result.  Can be one of V1 or V2.                                M
*   V1, V2:  Two input vectors to blend.				     M
*   t:	     Blending factor.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void			                                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecSqrLength, MvarVecLength, MvarVecScale                M
*   MvarVecNormalize, MvarVecDotProd                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecBlend                                                             M
*****************************************************************************/
void MvarVecBlend(MvarVecStruct *VRes,
		  const MvarVecStruct *V1,
		  const MvarVecStruct *V2,
		  CagdRType t)
{
    int i;
    CagdRType
	t1 = 1.0 - t;

    assert(V1 -> Dim == V2 -> Dim && VRes -> Dim == V2 -> Dim);

    for (i = 0; i < V1 -> Dim; i++)
        VRes -> Vec[i] = V1 -> Vec[i] * t + V2 -> Vec[i] * t1;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Normalize a given multivariate vector to a unit length, in place.        M
*                                                                            *
* PARAMETERS:                                                                M
*   V:   Vector to normalize.                                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:  TRUE if successful, FALSE if the input is the ZERO vector.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecAdd, MvarVecDotProd, MvarVecSqrLength, MvarVecLength,             M
*   MvarVecScale                                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecNormalize                                                         M
*****************************************************************************/
int MvarVecNormalize(MvarVecStruct *V)
{
    CagdRType
	VecLengthFactor = MvarVecLength(V);

    if (VecLengthFactor > IRIT_UEPS) {
	MvarVecScale(V, 1.0 / VecLengthFactor);
	return TRUE;
    }
    else
        return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Updates Dir to be the closest vector to Dir that is orthogonal to Vec.   M
*   In essence, apply a Graham Shmidt step.  Vectors need not be unit size.  M
*                                                                            *
* PARAMETERS:                                                                M
*   Dir:   Vector to update in place so it will be orthgoonal to Vec.        M
*   Vec:   Vector to make sure Dir is made orthogonal to.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if successful, FALSE otherwise.                              M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecOrthogonal2, MvarVecSetOrthogonalize, MvarVecWedgeProd            M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecOrthogonalize                                                     M
*****************************************************************************/
int MvarVecOrthogonalize(MvarVecStruct *Dir, const MvarVecStruct *Vec)
{
    int i,
        Dim = Dir -> Dim;
    CagdRType R,
        LenVec = MvarVecLength(Vec);

    assert(Vec -> Dim == Dim);

    if (LenVec == 0.0)
        return TRUE;

    R = MvarVecDotProd(Dir, Vec) / IRIT_SQR(LenVec);

    for (i = 0; i < Dim; i++)
        Dir -> Vec[i] -= Vec -> Vec[i] * R;

    assert(IRIT_APX_EQ(MvarVecDotProd(Dir, Vec), 0.0));

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Derives a unit vector Dir that is orthogonal to both Vec1, and Vec2.     M
*   Note that in R^2 there is no such vector, in R^3 only one such vector    M
* and in R^n, n > 3, there are infinitly many of which we find one.          M
*                                                                            *
* PARAMETERS:                                                                M
*   Dir:          Newly computed unit vector will be kept here.              M
*   Vec1, Vec2:   Two vectors we must be orthogonal to. Assumed unit length. M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if successful, FALSE otherwise.                              M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecOrthogonalize, MvarVecSetOrthogonalize, MvarVecWedgeProd          M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecOrthogonal2                                                       M
*****************************************************************************/
int MvarVecOrthogonal2(MvarVecStruct *Dir,
		       const MvarVecStruct *Vec1,
		       const MvarVecStruct *Vec2)
{
    int i, j,
	Dim = Dir -> Dim;
    MvarVecStruct *Vec2O;

    if (Dim < 3)
        return FALSE;

    IritRandomInit(301060); 	       /* Make the randomness consistent... */

    /* Make Vec2 and orthogonal vector to Vec1. */
    Vec2O = MvarVecCopy(Vec2);
    MvarVecOrthogonalize(Vec2O, Vec1);

    /* Create a random vector and make sure its projection on Vec1/2 is     */
    /* large enough (and if not randomize again.			    */
    for (i = 0; i < 1000; i++) {  /* Allow that many trials before failing. */
	CagdRType d1, d2;

        for (j = 0; j < Dim; j++)
	    Dir -> Vec[j] = IritRandom(-1.0, 1.0);

	MvarVecNormalize(Dir);
	d1 = MvarVecDotProd(Dir, Vec1);
	d2 = MvarVecDotProd(Dir, Vec2O);

	if (IRIT_ABS(d1) < MVAR_ORTHO2_VALID_VAL &&
	    IRIT_ABS(d2) < MVAR_ORTHO2_VALID_VAL &&
	    MvarVecOrthogonalize(Dir, Vec1) &&
	    MvarVecOrthogonalize(Dir, Vec2O)) {
	    assert(IRIT_APX_EQ(MvarVecDotProd(Dir, Vec1), 0.0) &&
		   IRIT_APX_EQ(MvarVecDotProd(Dir, Vec2), 0.0));

#	    ifdef DEBUG
	        if (i > 2)
	            fprintf(stderr, "Orthogonalize the vector in MvarVecOrthogonal2 after %d steps.\n", i);
#	    endif /* DEBUG */

	    MvarVecFree(Vec2O);

	    return TRUE;
	}
    }

#ifdef DEBUG
    fprintf(stderr, "Failed to orthogonalize the vector in MvarVecOrthogonal2.\n");
#endif /* DEBUG */

    MvarVecFree(Vec2O);

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Update the given set of vectors to be of unit size and orthogonal to     M
* each other.  Vectors are all assumed of the same dimension.		     M
*   In essence, apply a Graham Shmidt to all vectors.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Vecs:  Input vectors to make orthonormal.			             M
*   OrthoVecs:   Output vectors that span the same (sub) space as vec but    M
*                are orthogonal and unit length (orthognormal).  	     M
*                   Can be the same as Vecs.			             M
*   Size:  Number of vectors in Vecs and OrthoVecs vectors of vectors.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if successful, FALSE otherwise.                              M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecOrthogonal, MvarVecOrthogonal2, MvarVecWedgeProd                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecSetOrthogonalize                                                  M
*****************************************************************************/
int MvarVecSetOrthogonalize(const MvarVecStruct **Vecs,
			    MvarVecStruct **OrthoVecs,
			    int Size)
{
    int i, j, k, n;

    /* Do the Graham Schmidt process. */
    for (j = 0; j < Size; j++)
	MVAR_VEC_COPY(OrthoVecs[j], Vecs[j]);

    for (n = 0; n < MVAR_VEC_MAX_ORTHO_ITERS; n++) {
	CagdBType
	    Ortho = TRUE;

	for (j = 0; j < Size; j++) {
	    for (i = 0; i < j; i++) { 
		MvarVecAddScale(OrthoVecs[j], OrthoVecs[j], OrthoVecs[i],
				-MvarVecDotProd(OrthoVecs[j], OrthoVecs[i]));
	    }

	    if (!MvarVecNormalize(OrthoVecs[j])) {
		fprintf(stderr, "VecSetOrtho: Input vectors are linearly dependant.\n");
		for (k = 0; k < OrthoVecs[i] -> Dim; k++)/* Null the vector. */
		    OrthoVecs[i] -> Vec[k] = 0.0;
	    }
	}

	/* Evaluate the quality of the result. */
	for (i = 1; i < Size; i++) {
	    for (j = 0; j < i; j++) {
		CagdRType
		    R = MvarVecDotProd(OrthoVecs[i], OrthoVecs[j]);

		if (IRIT_ABS(R) > MVAR_VEC_NUMER_TOL)
		    Ortho = FALSE;
	    }
	}

	if (Ortho)
	    break;
    }

#ifdef DEBUG
    if (n > 0) {
	for (i = 1; i < Size; i++) {
	    for (j = 0; j < i; j++) {
		CagdRType
		    R = MvarVecDotProd(OrthoVecs[i], OrthoVecs[j]);

		if (R > IRIT_EPS) {
		    fprintf(stderr, "Inner product failed: %15.13f\n", R);
		}
	    }
	}
    }
#endif /* DEBUG */

    if (n >= MVAR_VEC_MAX_ORTHO_ITERS) {
	fprintf(stderr, "VecSetOrtho: Failed to orthogonolize the vectors.\n");
	return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes an orthogonal complement (wedge product) of dimension           M
* NewSize-Size to the given set of Size vectors.  			     M
*   Vectors is vector Vectors are in a higher dimensional linear space.	     M
* of at list dimension NewSize.                                              M
*                                                                            *
* PARAMETERS:                                                                M
*   Vectors:   The set of Size vectors, each of dimension larger than Size.  M
*              This set of vectors is modified in place - it is made         M
*              orthonormal.						     M
*   Size:      The size of Vectors set.					     M
*   NewVecs:   New vectors to allocate into here and update with orthogonal  M
*              complement subspace of dimension newSize-Size.		     M
*              Should be able to hold NewSize-Size pointers to vectors.      M
*   NewSize:   The new size of the set Vectors and NewVecs together.         M
*   CheckDet:  When Size + NewSize == Dim, this flag indicates if the        M
*              determinant of the Dim vectors is to be evaluated. Useful     M
*              for orientation issues.                                       M
*   DetVal:    Output value, the computed determinant if CheckDet is TRUE.   M
*									     *
* RETURN VALUE:                                                              M
*   CagdBType:	TRUE if successful, FALSE otherwise.	                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecOrthogonal, MvarVecOrthogonal2, MvarVecSetOrthogonalize           M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarVecWedgeProd                                                         M
*****************************************************************************/
CagdBType MvarVecWedgeProd(MvarVecStruct **Vectors,
			   int Size, 
			   MvarVecStruct **NewVecs,
			   int NewSize, 
			   CagdBType CheckDet,
			   CagdRType *DetVal)
{
    int i, n,
	OrigSize = Size,
	OrthoSize = NewSize - Size,
        Dim = Vectors[0] -> Dim;
    MvarVecStruct **AugmentedVectors;

    assert(Size < Dim + 1 && NewSize > Size && NewSize <= Dim);

    MvarVecSetOrthogonalize((const MvarVecStruct **) Vectors, Vectors, Size);

    AugmentedVectors = (MvarVecStruct **)
                               IritMalloc(sizeof(MvarVecStruct *) * NewSize);
    IRIT_GEN_COPY(AugmentedVectors, Vectors, sizeof(MvarVecStruct *) * Size);

    /* COmpute the NewSize-Size new vectors: */
    for (n = 0; n < OrthoSize; n++) {
        CagdBType
	    TryAnother = TRUE;
	MvarVecStruct
	    *RandomVec = MvarVecNew(Dim);

	while (TryAnother) {
	    for (i = 0; i < Dim; i++)
	        RandomVec -> Vec[i] = IritRandom(-1.0, 1.0);
	    TryAnother = FALSE;

	    for (i = 0; i < Size; i++) {
	        MvarVecAddScale(RandomVec, RandomVec, AugmentedVectors[i],
			     -MvarVecDotProd(RandomVec, AugmentedVectors[i]));
	    }

	    if (!MvarVecNormalize(RandomVec))
	        TryAnother = TRUE;

	    for (i = 0; i < Size; i++) {
	        CagdRType
		    R = MvarVecDotProd(AugmentedVectors[i], RandomVec);

		if (fabs(R) >= IRIT_EPS)
		    TryAnother = TRUE;
	    }
	}

#	ifdef DEBUG
	    for (i = 0; i < Size; i++) {
	        CagdRType
		    R = MvarVecDotProd(AugmentedVectors[i], RandomVec);

		assert(fabs(R) < IRIT_EPS);
	    }
#	endif /* DEBUG */

	AugmentedVectors[Size++] = NewVecs[n] = RandomVec;
    }

    if (CheckDet) {
	IrtGnrlMatType
	    RNBasisMat = (IrtGnrlMatType) IritMalloc(IRIT_SQR(Dim) *
	                                                     sizeof(IrtRType));
    
	for (i = 0; i < OrigSize; i++) {
	    IRIT_GEN_COPY(&RNBasisMat[Dim * i], Vectors[i] -> Vec, 
		          sizeof(CagdRType) * Dim);
	}

	for (i = OrigSize; i < NewSize; i++) {
	    IRIT_GEN_COPY(&RNBasisMat[Dim * i], NewVecs[i - OrigSize] -> Vec, 
		          sizeof(CagdRType) * Dim);
	}

	*DetVal = MatGnrlDetMatrix(RNBasisMat, Dim);

	IritFree(RNBasisMat);
    }
    IritFree(AugmentedVectors);
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Normalize a given multivariate plane's normal direction to a unit        M
* length, in place.						             M
*                                                                            *
* PARAMETERS:                                                                M
*   Pln:   Plane to normalize its normal direction.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:  TRUE if successful, FALSE if the input is the ZERO vector.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarVecNormalize                                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarPlaneNormalize                                                       M
*****************************************************************************/
int MvarPlaneNormalize(MvarPlaneStruct *Pln)
{
    int i,
	Dim = Pln -> Dim;
    CagdRType
	DProd = 0.0,
	*V = Pln -> Pln;

    for (i = 0; i < Dim; i++, V++)
	DProd += IRIT_SQR(*V);

    if (DProd != 0.0) {
        DProd = 1.0 / sqrt(DProd);

	for (i = 0, V = Pln -> Pln; i < Dim; i++)
	    *V++ *= DProd;

	return TRUE;
    }
    else
        return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Compute the intersection of a line and a hyperplane, in R^Dim.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   P, V:   Point and direction of line to intersect with hyperplane.        M
*	    Both P and V are of length Dim.				     M
*   Pln:    Hyperplane to intersect with line.                               M
*   Param:  Will be updated with the parameter along which the intersection  M
*           has occured, as Inter = P + V * t.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarVecStruct *:  The intersection point.                                M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarLinePlaneInter                                                       M
*****************************************************************************/
MvarVecStruct *MvarLinePlaneInter(const MvarVecStruct *P,
				  const MvarVecStruct *V,
				  const MvarPlaneStruct *Pln,
				  CagdRType *Param)
{
    int i,
	Dim = P -> Dim;
    CagdRType t,
	A = 0.0,       /* We will reduce problem to linear form: At + B = 0. */
        B = 0.0;
    MvarVecStruct *InterPt;

    /* Substitute the line as "P + V * t" into the plane and find t. */
    for (i = 0; i < Dim; i++) {
        A += V -> Vec[i] * Pln -> Pln[i];
	B += P -> Vec[i] * Pln -> Pln[i];
    }

    t = A == 0.0 ? IRIT_INFNTY : -(B + Pln -> Pln[Dim]) / A;

    /* Compute the intersecting location. */
    InterPt = MvarVecNew(Dim);
    for (i = 0; i < Dim; i++)
        InterPt -> Vec[i] = P -> Vec[i] + V -> Vec[i] * t;

    return InterPt;
}
