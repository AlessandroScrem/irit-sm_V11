/******************************************************************************
* ZrMvKant.c - Uniqueness of roots using Kantorovich theorem.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Carmi Grushko, Feb. 2011.					      *
******************************************************************************/

#include "mvar_loc.h"

IRIT_GLOBAL_DATA int
    _MVGlblZeroApplyKantorovichTest = FALSE,
     /* Should Kantorovich subdivided around the unique root's box ? */
    _MVGlblZeroApplyKantorovichTestAroundBox = FALSE;
IRIT_GLOBAL_DATA 
    IrtRType KantBallFractionThreshold = 0.1;


/* A square, symmetric matrix data structure, implemented using a            */
/* "triangular" array.							     */
typedef struct {
    CagdRType *Data;
    int N;
} TriangularArrayStruct;

static TriangularArrayStruct *TriangularArrayCreateAux(int N);
static void TriangularArrayDestroy(TriangularArrayStruct *Matrix);
static CagdRType *TriangularArrayAtAux(TriangularArrayStruct *Matrix,
				       int r,
				       int c);
static CagdRType MvarKantMaxSqrDerive(const MvarMVStruct *MV, 
				      MvarMVDirType DirectionB);
static void MvarKantCachePab(const MvarMVStruct *MV, 
    			     MvarMVGradientStruct* Gradient, 
		             TriangularArrayStruct *PCache);
static CagdRType MvarKantBoundSplineLipschitz(MvarMVDirType DirectionI,  
					      TriangularArrayStruct* PCache);
static CagdRType MvarKantBoundJacobian(MvarMVStruct * const *F, 
    				       MvarMVGradientStruct** Gradients, 
				       int Dim);
static CagdRType MvarKantMatrixInfNorm(IrtGnrlMatType Matrix, int Dim);
static IrtBType MvarKantComputeBetaEtha(MvarMVStruct * const * F, 
					MvarMVGradientStruct** Gradients, 
					int Dim, 
					CagdRType *Beta, 
					CagdRType *Etha, 
					CagdRType Pt[]);
static IrtBType MvarKantApplyKantorovichTest(MvarMVStruct * const * F, 
					     int n, 
					     CagdRType DivisionPoint[], 
					     CagdRType *DivisionRadius);
static void MvarZeroSubdivideAroundBoxAux2(MvarMVStruct *MVsSubdomain1[], 
					   MvarMVStruct *MVsSubdomain2[], 
					   int NumOfMVs, 
					   IrtRType t, 
					   int Dir, 
					   IrtBType KeepWorkingOnUpper);
static MvarMVStruct ***MvarKantSubdivideAroundBox(MvarMVStruct **MVs, 
						  int NumOfMVs, 
						  IrtRType DivisionPoint[], 
						  IrtRType DivisionRadius,
						  int* SubdomainMVsCount);
static MvarPtStruct *MvarZeroMVsSubdiv(MvarMVStruct **MVs,
				       MvarConstraintType *Constraints,
				       int NumOfMVs,
				       int NumOfZeroMVs,
				       int ApplyNormalConeTest,
				       CagdRType SubdivTol,
				       int Depth,
				       CagdBType SameSpace,
				       CagdRType ParamPerturb);
static MvarPtStruct *MvarZeroGenPtMidMvar(const MvarMVStruct *MV,
					  int SingleSol);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates an array of size that is enough for storage of an N by N         *
* square, symmetric matrix.                                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   N:   Dimension of the square, symmetric matrix.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType*:  A pointer to the create array                               *
*****************************************************************************/
static TriangularArrayStruct *TriangularArrayCreateAux(int N)
{
    TriangularArrayStruct
        *Matrix = (TriangularArrayStruct *)
                                    IritMalloc(sizeof(TriangularArrayStruct));

    Matrix -> Data = (CagdRType *) IritMalloc((N * (1 + N) / 2) *
					                   sizeof(CagdRType));
    Matrix -> N = N;
    return Matrix;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Destroys a square, symmetric matrix                                      *
*                                                                            *
* PARAMETERS:                                                                *
*   Matrix: The matrix to destroy.		                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   None                                                                     *
*****************************************************************************/
static void TriangularArrayDestroy(TriangularArrayStruct *Matrix)
{
    IritFree(Matrix -> Data);
    IritFree(Matrix);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Accessor for a square, symmetric matrix                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   Matrix:  The matrix to access.		                             *
*   r:  Zero-based row number.                                               *
*   c:  Zero-based column number.                                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType *:   A pointer to the (r, c) element.                          *
*****************************************************************************/
static CagdRType *TriangularArrayAtAux(TriangularArrayStruct *Matrix,
				       int r,
				       int c)
{
    int n = Matrix -> N;

    if (r > c)
	IRIT_SWAP(int, r, c);

    return &(Matrix -> Data[c - r + r * (n + n - r + 1) / 2]);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a direction, computes the derivative in this direction and finds   *
* the maximal absolute value of the control points. Assumes all control      *
* points are scalars.		                                             *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:          The multi-variate to operate on.                            *
*   Dir:         Direction to differentiate.                                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:   Maximal absolute value of the control points of the         *
*                derivative.                                                 *
*****************************************************************************/
static CagdRType MvarKantMaxSqrDerive(const MvarMVStruct *MV, 
				      MvarMVDirType Dir)
{
    CagdRType
        MaxP = 0;
    MvarMVStruct
        *dMV = MvarMVDerive(MV, Dir);
    int i,
	PCount = dMV -> SubSpaces[dMV -> Dim];

    for (i = 0; i < PCount; i++) {
	CagdRType
	    P = dMV -> Points[1][i];

	if (P < 0)
	    P = -P;
	if (P > MaxP)
	    MaxP = P;
    }

    MvarMVFree(dMV);

    return MaxP * MaxP;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Runs over all pairs of directions a,b such that a<b, and computes a      *
* bound on max{ d MV / da db } using the control points.                     *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:      The multi-variate to operate on.                                *
*   Pcache:  The square, symmetric matrix to write results into.             *
*                                                                            *
* RETURN VALUE:                                                              *
*   None                                                                     *
*****************************************************************************/
static void MvarKantCachePab(const MvarMVStruct *MV, 
    			     MvarMVGradientStruct* Gradient, 
		             TriangularArrayStruct *PCache)
{
    int a,
	n = MV -> Dim;

    for (a = 0; a < n; a++) {
	int b;
	MvarMVStruct
	    *dMVa = MvarMVDerive(MV, a);

	for (b = a; b < n; b++)
	    *TriangularArrayAtAux(PCache, a, b) = MvarKantMaxSqrDerive(dMVa,
								       b);
	MvarMVFree(dMVa);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes a bound for Lipschitz constant of d MV / da, based on the       *
* bounds on the second derivatives of MV in PCache.                          *
*                                                                            *
* PARAMETERS:                                                                *
*   Dir:         Direction of differentiation of MV.                         *
*   Pcache:      The square, symmetric matrix that holds bounds on           *
*	  	  max{ d MV / da db }.		                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:  A bound for the Lipschitz constant of \partial_{x_i} f_j.    *
*****************************************************************************/
static CagdRType MvarKantBoundSplineLipschitz(MvarMVDirType Dir,  
					      TriangularArrayStruct *PCache)
{
    int u,
	n = PCache -> N;
    CagdRType
        Sum = 0;

    for (u = 0; u < n; u++)
	Sum += *TriangularArrayAtAux(PCache, Dir, u);

    return sqrt(Sum);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes a bound for Lipschitz constant of the Jacobian of the system,   *
* based on Convex Hull property.                                             *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:    Set of functions defining the system.                            *
*   Dim:    Dimension of the problem (= number of function in the system).   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: A bound for the Lipschitz constant of the Jacobian of the     *
*              system.                                                       *
*****************************************************************************/
static CagdRType MvarKantBoundJacobian(MvarMVStruct * const *MVs, 
    				       MvarMVGradientStruct** Gradients, 
				       int Dim)
{    
    int i, j;
    CagdRType
        Max = 0;
    IrtGnrlMatType
        J = (IrtGnrlMatType) IritMalloc(Dim * Dim * sizeof(CagdRType));
    TriangularArrayStruct
        *Cache = TriangularArrayCreateAux(Dim);

    for (i = 0; i < Dim; i++) {
	MvarKantCachePab(MVs[i], Gradients[i], Cache);

	for (j = 0; j < Dim; j++) {
	    J[i * Dim + j] = MvarKantBoundSplineLipschitz(j, Cache);
	}
    }
	
    for (i = 0; i < Dim; i++) {		
	CagdRType
	    Res = 0;

	for (j = 0; j < Dim; j++)
	    Res += J[i * Dim + j];

	if (Res > Max)
	    Max = Res;
    }

    TriangularArrayDestroy(Cache);
    IritFree(J);

    return Max * sqrt((CagdRType) Dim);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the L-inifinity (max) norm of a (square) matrix                 *
*                                                                            *
* PARAMETERS:                                                                *
*   Matrix:   The matrix to operate on.                                      *
*   Dim:      Dimension of the matrix.                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:  L-inifinity norm of the matrix.                              *
*****************************************************************************/
static CagdRType MvarKantMatrixInfNorm(IrtGnrlMatType Matrix, int Dim)
{
    int i, j;
    CagdRType
        Max = 0;

    for (i = 0; i < Dim; i++) {
	CagdRType
	    Sum = 0;

	for (j = 0; j < Dim; j++)
	    Sum += IRIT_ABS(Matrix[i * Dim + j]);

	Max = IRIT_MAX(Max, Sum);
    }

    return Max;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes Beta and Etha of a system at a specified point. For more        *
* information, refer to the document about using the Kantorovich theorem.    *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:   Functions that define the system.                                 *
*   Dim:   Dimension of the problem.                                         *
*   Beta:  A pointer to receive the computed Beta.                           *
*   Etha:  A poitner to receieve the computed Etha.                          *
*   Pt:    The coordinates to evaluate Beta and Etha at.                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType:  TRUE if successfull, FALSE if the specified point cannot      *
*              be used (because the Jacobian is singular).                   *
*****************************************************************************/
static IrtBType MvarKantComputeBetaEtha(MvarMVStruct * const * MVs, 
					MvarMVGradientStruct **Gradients, 
					int Dim, 
					CagdRType *Beta, 
					CagdRType *Etha, 
					CagdRType Pt[])
{
    IrtBType Result;
    int i;
    IrtGnrlMatType 
	J = NULL, 
	MVsVal = NULL,
	invJ = NULL;

    /* Fill out a Jacobian matrix, and a MVs(Pt) matrix. */
    J = (IrtGnrlMatType) IritMalloc(sizeof(CagdRType) * Dim * Dim);
    MVsVal = (IrtGnrlMatType) IritMalloc(sizeof(CagdRType) * Dim);

    for (i = 0; i<Dim; i++) {
	CagdRType
	    *v = MvarMVEvalGradient(Gradients[i], Pt, 0);

	IRIT_GEN_COPY(J + i * Dim, v, Dim * sizeof(CagdRType));
	MVsVal[i] = v[Dim];
    }

    /* Compute inverse Jacobian. */
    if (Dim == 2) {
	IrtRType 
	    Det = J[0] * J[3] - J[1] * J[2];

	if (IRIT_ABS(Det) < IRIT_UEPS)
	    /* Jacobian is singular. */
	    Result = FALSE; 
	else {	    
	    Result = TRUE; 
	    IRIT_SWAP(IrtRType, J[0], J[3]);
	    J[1] = -J[1];
	    J[2] = -J[2];

	    J[0] /= Det;
	    J[1] /= Det;
	    J[2] /= Det;
	    J[3] /= Det;

	    invJ = J;
	}
    }
    else if (Dim == 3) {
	IrtRType 
	    A = J[4] * J[8] - J[5] * J[7],
	    B = J[5] * J[6] - J[3] * J[8],
	    C = J[3] * J[7] - J[4] * J[6],

	    D = J[2] * J[7] - J[1] * J[8],
	    E = J[0] * J[8] - J[2] * J[6],
	    F = J[1] * J[6] - J[0] * J[7],

	    G = J[1] * J[5] - J[2] * J[4],
	    H = J[2] * J[3] - J[0] * J[5],
	    K = J[0] * J[4] - J[1] * J[3],

	    Det = J[0] * A + J[1] * B + J[2] * C,

	    InvDet = 1/Det;

	if (IRIT_ABS(Det) < IRIT_UEPS)
	    /* Jacobian is singular */
	    Result = FALSE;
	else {
	    Result = TRUE; 
	    J[0] = A * InvDet;
	    J[1] = D * InvDet;
	    J[2] = G * InvDet;
	    J[3] = B * InvDet;
	    J[4] = E * InvDet;
	    J[5] = H * InvDet;
	    J[6] = C * InvDet;
	    J[7] = F * InvDet;
	    J[8] = K * InvDet;
	    invJ = J;
	}
    } 
    else {
        invJ = (IrtGnrlMatType) IritMalloc(sizeof(CagdRType) * Dim * Dim);
	Result = MatGnrlInverseMatrix(J, invJ, Dim);
    }

    if (Result) {
	IrtRType Sum;
	int r, c;
	IrtVecGnrlType 
	    VT = invJ;

	/* Compute Beta. */
	*Beta = MvarKantMatrixInfNorm(invJ, Dim);

	/* Compute Etha. */
	*Etha = 0;
	for (r=0; r<Dim; r++) {	    
	    Sum = 0;
	    for (c = 0; c < Dim; c++, VT++)
		Sum += *VT * MVsVal[c];

	    Sum = IRIT_ABS(Sum);
	    if (Sum > *Etha)
		*Etha = Sum;
	}
    }

    /* Free resources. */
    IritFree(MVsVal);
    if (invJ != J)
	IritFree(invJ);
    IritFree(J);

    return Result;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Apply the Kantorovich test to determine if there is an area in the       *
* currrent domain that has a unique root.                                    *
*                                                                            *
* PARAMETERS:                                                                *
*   F:            Functions that define the system.                          *
*   n:            Dimension of the problem (= number of functions).          *
*   DivisionPoint, DivisionRadius:  A pointer to return the Radius in which  *
*      there is a unique root, around the point specified in DivisionPoint.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType:  TRUE if there is unique root somewhere, FALSE otherwise.      *
*****************************************************************************/
static IrtBType MvarKantApplyKantorovichTest(MvarMVStruct * const * F, 
					     int n, 
					     CagdRType DivisionPoint[], 
					     CagdRType *DivisionRadius)
{
    int i,
        Attempts = 10;
    CagdRType Gamma,
	Beta = 1.0,
	Etha = 1.0,
        *Pt = (CagdRType *) IritMalloc(sizeof(CagdRType) * n),
	*MinArray = (CagdRType *) IritMalloc(sizeof(CagdRType) * n),
	*MaxArray = (CagdRType *) IritMalloc(sizeof(CagdRType) * n),
	MaxVolume = 0;
    MvarMVGradientStruct
        **Gradients = (MvarMVGradientStruct **)
                            IritMalloc(sizeof(MvarMVGradientStruct *) * n);

    for (i = 0; i < n; i++)
	 Gradients[i] = MvarMVPrepGradient(F[i], TRUE);

    /* Compute Gamma. */
    Gamma = MvarKantBoundJacobian(F, Gradients, n);

    /* Compute Beta, Etha at points sampled uniformly in the domain	    */
    /* of the first function (which should be the same for all functions.   */
    MvarMVDomain(F[0], MinArray, MaxArray, -1);

    *DivisionRadius = -1;

    for (; Attempts; Attempts--) {
	CagdRType Alpha;

	for (i = 0; i<n; i++)
	    Pt[i] = IritRandom(MinArray[i], MaxArray[i]);

	if (!MvarKantComputeBetaEtha(F, Gradients, n, &Beta, &Etha, Pt))
	    continue;                              /* Jacobian is singular. */

	/* Compute Alpha and the uniqueness radius, if exists. */
	Alpha = Beta * Etha * Gamma;

	if (Alpha <= 0.5) {
	    IrtBType 
		Reject = TRUE;
	    CagdRType
	        Radius = (1 + sqrt(1 - 2 * Alpha)) / Beta / Gamma,
	        volume = 1;

	    for (i = 0; i < n; i++) {
		volume *= IRIT_MIN(Pt[i] + Radius, MaxArray[i]) +
		          IRIT_MAX(Pt[i] - Radius, MinArray[i]);

		if ((MaxArray[i] - Pt[i]) > Radius ||
		    (Pt[i] - MinArray[i]) > Radius)
		    Reject = FALSE;
	    }

	    if (MaxVolume < volume) {
		volume = MaxVolume;
		*DivisionRadius = Radius;
		IRIT_GEN_COPY(DivisionPoint, Pt, sizeof(CagdRType) * n);
	    }

	    MaxVolume = IRIT_MAX(MaxVolume, volume);

	    if (Reject) {
		IritFree(MinArray);
		IritFree(MaxArray);
		IritFree(Pt);
		for (i = 0; i < n; i++)
		    MvarMVFreeGradient(Gradients[i]);
		IritFree(Gradients);
		return TRUE;
	    }
	}
    }

    IritFree(MinArray);
    IritFree(MaxArray);
    IritFree(Pt);
    for (i = 0; i < n; i++)
	MvarMVFreeGradient(Gradients[i]);
    IritFree(Gradients);
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Used in MvarZeroSubdivideAroundBoxAux. Subdivides MVsSubdomain1 along    *
*  direction Dir at coordinate t. Depending on KeppWorkingOnUpper,           *
* MVsSubdomain1 becomes the subdomain with a coordinate smaller than t in    *
* direction Dir, while MVsSubdomain2 becomes the subdomain with a coordinate *
* larger than t in that direction.		                             *
*                                                                            *
* PARAMETERS:                                                                *
*   MVsSubdomain1:      [in/out] Vector of multivariate constraints.         *
*   MVsSubdomain2:      [out]    Vector of multivariate constraints.         *
*   NumOfMVs:                    Size of the MVs and Constraints vector.     *
*   Dir:                         Subdivide along this direction.             *
*   t:                           Subdivide at this parametric coordinate.    *
*   KeepWorkingOnUpper:          If TRUE, MVsSubdomain1 becomes with smaller *
*				 than t in Dir, otherwise MVsSubdomain2      *
*                                                                            *
* RETURN VALUE:                                                              *
*   None                                                                     *
*****************************************************************************/
static void MvarZeroSubdivideAroundBoxAux2(MvarMVStruct *MVsSubdomain1[], 
					   MvarMVStruct *MVsSubdomain2[], 
					   int NumOfMVs, 
					   IrtRType t, 
					   int Dir, 
					   IrtBType KeepWorkingOnUpper)
{
    int i;

    if (KeepWorkingOnUpper) {
	for (i = 0; i < NumOfMVs; i++) {
	    MVsSubdomain2[i] = MvarMVSubdivAtParam(MVsSubdomain1[i], t, Dir);
	    MVsSubdomain1[i] = MVsSubdomain2[i] -> Pnext;
	    MVsSubdomain2[i] -> Pnext = NULL;
	}
    }
    else {
	for (i = 0; i < NumOfMVs; i++) {
	    MVsSubdomain1[i] = MvarMVSubdivAtParam(MVsSubdomain1[i], t, Dir);
	    MVsSubdomain2[i] = MVsSubdomain1[i] -> Pnext;
	    MVsSubdomain1[i] -> Pnext = NULL;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given an n-cube, subdivides the domain around this cube (the cube        *
* itself is "discarded"). MVs is destroyed in the process, while the         *
* resulting sub-domains are returned.                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   MVs:                     Multivariates to subdivide.                     *
*   NumOfMVs:                Size of the MVs and Constraints vector.         *
*   DivisionPoint[]:         Center of box to subdivide around.              *
*   DivisionRadius:          Radius (in L-inf norm) of the box to subdivide. *
*   SubdomainMVsCount: [out] number of resulting subdomains.                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct ***:        List of subdomains that are "around" the box    *
*****************************************************************************/
static MvarMVStruct ***MvarKantSubdivideAroundBox(MvarMVStruct **MVs, 
						  int NumOfMVs, 
						  IrtRType DivisionPoint[], 
						  IrtRType DivisionRadius,
						  int *SubdomainMVsCount)
{
    int Dir, mv, i,
	Dim = MVs[0] -> Dim,
        MVsArraySize = 0;
    MvarMVStruct
	/* Simple subdivision along some direction results in splitting MVs */
	/* into 2 new MVs; subdivision around an n-dimensional box results  */
	/* in up to 2*n new MVs.                                            */
        ***MVsArray = 
	    (MvarMVStruct ***) IritMalloc(2 * Dim * sizeof(MvarMVStruct **)),
        **InternalMVs = 
	    (MvarMVStruct **) IritMalloc(NumOfMVs * sizeof(MvarMVStruct *));
    IrtRType 
	*MinArray = (IrtRType *) IritMalloc(sizeof(IrtRType) * Dim),
	*MaxArray = (IrtRType *) IritMalloc(sizeof(IrtRType) * Dim);

    MvarMVDomain(MVs[0], MinArray, MaxArray, -1);

    for (i=0; i<NumOfMVs; i++)
	InternalMVs[i] = MvarMVCopy(MVs[i]);

    /* Create sub-domain MVs. */
    for (Dir=0; Dir<Dim; Dir++) {
	int i, sign;

	for (i=0, sign=1; i<2; i++, sign *= -1) {
	    IrtRType 
		t = DivisionPoint[Dir] + sign * DivisionRadius;

	    /* Subdivide in direction Dim at Pt+Radius and Pt-Radius. */
	    if (t < MaxArray[Dir]) {
		MVsArray[MVsArraySize] = 
		    (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));
		MvarZeroSubdivideAroundBoxAux2(InternalMVs, 
					       MVsArray[MVsArraySize], 
					       NumOfMVs, t, Dir, 
					       sign != 1);
		MVsArraySize++;
	    }
	}
    }

    for (mv = 0; mv < NumOfMVs; mv++)
	MvarMVFree(InternalMVs[mv]);

    IritFree(MinArray);
    IritFree(MaxArray);
    IritFree(InternalMVs);

    *SubdomainMVsCount = MVsArraySize;

    return MVsArray;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Handle all aspects of root-uniqueness using Kantorovich; return a        M
* candidate point where a unique root exists, and find the rest of the roots M
* recursively, using MvarZeroMVsSubdiv.                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*   NumOfMVs:     Size of the MVs and Constraints vector.                    M
*   NumOfZeroMVs: Number of zero or equality constraints.                    M
*   ApplyNormalConeTest:  TRUE to apply normal cones' single intersection    M
*		  tests.						     M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   Depth:        Of subdivision recursion.				     M
*   SameSpace:    True if all MVs share the same function space.	     M
*   ParamPerturb: The perturbation to applly to mid. subdivision location.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   List of points on the solution set.  Dimension of the  M
*		      points will be the same as the dimensions of all MVs.  M
*                     NULL is returned.					     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros                                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroGetRootsByKantorovich                                            M
*****************************************************************************/
MvarPtStruct *MvarZeroGetRootsByKantorovich(MvarMVStruct **MVs,
					    MvarConstraintType *Constraints,
					    int NumOfMVs,
					    int NumOfZeroMVs,
					    int ApplyNormalConeTest,
					    CagdRType SubdivTol,
					    int Depth,
					    CagdBType SameSpace,
					    CagdRType ParamPerturb)
{
    int Mv, i,
	SubdomainMVsCount = 0,
	Dim = MVs[0] -> Dim;
    MvarPtStruct
	*PtList = NULL;
    MvarMVStruct
	***SubdomainMVs = NULL;	

    if (!_MVGlblZeroApplyKantorovichTest)
        return NULL;

    /* Kantorovich Stopping Criterion. */
    if (Dim == NumOfZeroMVs) {
	IrtRType 
	    *DivisionPt = (IrtRType *) IritMalloc(sizeof(IrtRType) * Dim);
	IrtRType DivisionRadius;
	if (MvarKantApplyKantorovichTest(MVs, NumOfZeroMVs,
					    DivisionPt, &DivisionRadius)) {
	    IritFree( DivisionPt );
	    return MvarZeroGenPtMidMvar(MVs[0], TRUE);
	}

	if (_MVGlblZeroApplyKantorovichTestAroundBox &&
	    DivisionRadius > 0 &&
	    pow(DivisionRadius, Dim) / MvarMVVolumeOfDomain(MVs[0], Dim)  
		                            > KantBallFractionThreshold) {

	    SubdomainMVs = MvarKantSubdivideAroundBox(MVs, NumOfMVs,
					              DivisionPt, 
						      DivisionRadius,
					              &SubdomainMVsCount);

	    /* Set a candidate pt. in middle of the root-uniqueness box. */
	    PtList = MvarPtNew (Dim);
	    PtList -> Pt = DivisionPt;
	}
    }

    /* Recursively find roots in sub-domain MVs. */
    for (Mv = 0; Mv < SubdomainMVsCount; Mv++) {
	MvarPtStruct
	    *TempPtList = NULL;	

	TempPtList = MvarZeroMVsSubdiv(SubdomainMVs[Mv], Constraints,
				       NumOfMVs, NumOfZeroMVs,
				       ApplyNormalConeTest, SubdivTol,
				       Depth + 1, SameSpace, ParamPerturb);
        PtList = (MvarPtStruct *) CagdListAppend(PtList, TempPtList);
    }

    /* Clean up - free everything. */
    for (i = 0; i < SubdomainMVsCount; i++) {
	for (Mv = 0; Mv < NumOfMVs; Mv++)
	    MvarMVFree(SubdomainMVs[i][Mv]);

	IritFree(SubdomainMVs[i]);
    }
    IritFree(SubdomainMVs);        

    /* Return candidate points. */
    return PtList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the use (or not) of the Kantorovich test inside the multivariate    M
* subdivisions' zero set solver.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   KantorovichTest:   New setting for normal cone testing usage.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       Old setting for normal cone testing usage.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros, MvarMVsZerosDomainReduction,				     M
*   MvarMVsZerosGradPreconditioning, MvarMVsZerosSetCallBackFunc	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVsZerosKantorovichTest                                              M
*****************************************************************************/
int MvarMVsZerosKantorovichTest(int KantorovichTest)
{
    int OldVal = _MVGlblZeroApplyKantorovichTest;

    _MVGlblZeroApplyKantorovichTest = KantorovichTest;

    return OldVal;
}


/*****************************************************************************
* DESCRIPTION:                                                               M
*   Approximate a solution to the set of constraints, if any, using the      M
* subdivision of the parametric domains of the MVs.  Stops when the          M
* parametric domain is smaller than SubdivTol in all dimensions and returns  M
* a central point to that small multivariate patch.			     M
*   NOTE: The generic solver invokes this function only via the Kantorovich  M
*   option (global flag set accordingly), otherwise this function is unused. M
*                                                                            *
* PARAMETERS:                                                                M
*   MVs:          Vector of multivariate constraints.                        M
*   Constraints:  Either an equality or an inequality type of constraint.    M
*   NumOfMVs:     Size of the MVs and Constraints vector.                    M
*   NumOfZeroMVs: Number of zero or equality constraints.                    M
*   ApplyNormalConeTest:  TRUE to apply normal cones' single intersection    M
*		  tests.						     M
*   SubdivTol:    Tolerance of the subdivision process.  Tolerance is        M
*		  measured in the parametric space of the multivariates.     M
*   Depth:        Of subdivision recursion.				     M
*   SameSpace:    True if all MVs share the same function space.             M
*   ParamPerturb: The perturbation to applly to mid. subdivision location.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   List of points on the solution set.  Dimension of the  M
*		      points will be the same as the dimensions of all MVs.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVsZeros                                                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroMVsSubdiv                                                        M
*****************************************************************************/
static MvarPtStruct *MvarZeroMVsSubdiv(MvarMVStruct **MVs,
				       MvarConstraintType *Constraints,
				       int NumOfMVs,
				       int NumOfZeroMVs,
				       int ApplyNormalConeTest,
				       CagdRType SubdivTol,
				       int Depth,
				       CagdBType SameSpace,
				       CagdRType ParamPerturb)
{
    int i, j, l, HasInternalKnot,
	Dim = MVs[0] -> Dim;
    CagdRType Min, Max;

#   ifdef MVAR_DEBUG_DEPTH
    if (Depth < MVAR_DEBUG_MAX_DEPTH)
        MvarSubdivLevel[Depth]++;
#   endif /* MVAR_DEBUG_DEPTH */

#   ifdef DEBUG_ONE_SOLUTION
    if (NumOfZeroMVs == 2 || NumOfZeroMVs == 3) {
        static int
	    Count = 1;
        CagdRType UMin, UMax, VMin, VMax,
	    WMin = 0.0,
	    WMax = 0.0,
	    u = 0.02223875187684,
	    v = 0.87209881542057,
	    w = 0.43073326427685;

        MvarGetSubdivParamDomains(MVs[0], &UMin, &UMax, 0);
        MvarGetSubdivParamDomains(MVs[0], &VMin, &VMax, 1);
	if (NumOfZeroMVs == 3)
	    MvarGetSubdivParamDomains(MVs[0], &WMin, &WMax, 2);

	if (UMin <= u && u <= UMax &&
	    VMin <= v && v <= VMax &&
	    (NumOfZeroMVs < 3 || (WMin <= w && w <= WMax)))
	    IRIT_INFO_MSG_PRINTF("In domain (%d) [%f %f %f] \t[%f %f %f]\n",
				  Count++, UMin, VMin, WMin, UMax, VMax, WMax);
    }
#   endif /* DEBUG_ONE_SOLUTION */

#   ifdef DEBUG_DUMP_DOMAINS
    {
        CagdRType UMin, UMax, VMin, VMax;

        MvarGetSubdivParamDomains(MVs[0], &UMin, &UMax, 0);
        MvarGetSubdivParamDomains(MVs[0], &VMin, &VMax, 1);

	fprintf(GlblDumpDomainsFile,
		"[OBJECT [RGB \"100,255,100\"] NONE    [POLYLINE 5\n\t[%f %f 0]\n\t[%f %f 0]\n\t[%f %f 0]\n\t[%f %f 0]\n\t[%f %f 0]\n    ]\n]\n",
		UMin, VMin, UMin, VMax, UMax, VMax, UMax, VMin, UMin, VMin);
    }
#   endif /* DEBUG_DUMP_DOMAINS */

#   ifdef DEBUG
	/* Test if the multivariate may satisfy their constraints.  Examine */
	/* positivity/negativity of the set of values that multivariate has */
	/* and check against the prescribed constraint to the multivariate. */
	for (i = 0; i < NumOfMVs; i++) {
	    if (MvarZeroMVConstraintFail(MVs[i], Constraints[i])) {
 		/* Should not happen according to our invariance:  Input     */
		/* should pass purging test successfully!		     */
		assert(0);
		return NULL;			     /* No solutions exists. */
	    }
	}
#   endif /* DEBUG */

    /* If no solution is possible, by applying the pairs of hyperplanes      */
    /* bounding test, quit right here.					     */
    if (_MVGlblZeroApplyParallelHyperPlaneTest &&
	Dim == NumOfZeroMVs &&
	!MVarMVHyperPlanesTestForSol((MvarMVStruct const * const *) MVs,
				     NumOfZeroMVs))
        return NULL;

    /* Check the normal cone overlapping criteria.			     */
    if (ApplyNormalConeTest && !MvarMVConesOverlap(MVs, NumOfZeroMVs)) {
#	ifdef DEBUG
        {
	    IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintZeroPts, FALSE) {
	        CagdRType Min, Max;

		IRIT_INFO_MSG("Zero by cones at");
		for (i = 0; i < Dim; i++) {
		    MvarGetSubdivParamDomains(MVs[0], &Min, &Max, i);
		    IRIT_INFO_MSG_PRINTF(" %f [%f %f]",
					 (Min + Max) * 0.5, Min, Max);
		}
		IRIT_INFO_MSG("\n");
	    }
	}
#       endif /* DEBUG */

        return MvarZeroGenPtMidMvar(MVs[0], TRUE);
    }

    /* If we got here then these patches may satisfy all constraints.        */
    /* Subdivide them along their maximal parametric length dimension.	     */
    /*   There is one exception to this rule and that is to subdivide first  */
    /* in internal knots so if the maximal parametric length dimension holds */
    /* no knots while another dir do has knots, the other dir will be used.  */
    l = -1;
    HasInternalKnot = FALSE;
    Max = Min = 0.0;
    for (i = 0; i < NumOfMVs && l < 0; i++) {
        for (j = 0; j < Dim && l < 0; j++) {
	    if (MVs[i] -> Lengths[j] != MVs[i] -> Orders[j]) {
	        CagdRType TMin, TMax;

	        /* Found a direction with internal knots. Save in l only if */
		/* this one is a larger domain than what was found so far.  */
	        MvarGetSubdivParamDomains(MVs[0], &TMin, &TMax, j);
		if (TMax - TMin > Max - Min) {
		    Max = TMax;
		    Min = TMin;
		    l = j;
		    HasInternalKnot = TRUE;
		}
	    }
	}
    }

    for (i = 0; i < Dim; i++) {
	CagdRType MinCurr, MaxCurr;

	MvarGetSubdivParamDomains(MVs[0], &MinCurr, &MaxCurr, i);
        if (MaxCurr - MinCurr > Max - Min) {
	    int j;

	    /* If we got internal knots with a domain larger than SubdivTol */
	    /* make sure this direction is a direction with knots.          */
	    if (HasInternalKnot && Max - Min > SubdivTol) {
	        for (j = 0; j < NumOfMVs; j++) {
		    if (MVs[j] -> Lengths[i] != MVs[j] -> Orders[i])
		        break;
		}
	    }
	    else
	        j = 0;

	    if (j < NumOfMVs) {
	        l = i;
		Min = MinCurr;
		Max = MaxCurr;
	    }
        }
    }
    assert(l >= 0);

    if (Max - Min > SubdivTol) {
        CagdBType
	    WasReduction = FALSE;
	CagdRType t,
	    TMin = MVAR_IS_BEZIER_MV(MVs[0]) ? 0.0 : Min,
	    TMax = MVAR_IS_BEZIER_MV(MVs[0]) ? 1.0 : Max;
	MvarPtStruct *PtList1, *PtList2;
	MvarMVStruct **MVs1, **MVs2;

	/* Apply the kantorovich test if so desired. */
	PtList1 = MvarZeroGetRootsByKantorovich(MVs,Constraints, NumOfMVs, 
	                                        NumOfZeroMVs, 
						ApplyNormalConeTest, 
						SubdivTol, Depth, 
						SameSpace, ParamPerturb);
	if (PtList1 != NULL)
	    return PtList1;

	/* MvarMVsReduceMvsDomains returns TMin/TMax in [0, 1] for Bezier   */
	/* and returns correct TMin/TMax for B-spline MVs.		    */
	/*   Also force subdivisions from time to time as domain reduction  */
	/* can fail to converge if applied by itself.			    */
	if (_MVGlblZeroApplyDomainReduction && (Depth & 0x03) != 0) {
	    CagdRType
	        OrigTMin = TMin,
	        OrigTMax = TMax;

	    if (!MvarMVsReduceMvsDomains(MVs, NumOfZeroMVs,
					 l, SubdivTol, &TMin, &TMax,
					 SameSpace))
		return NULL;    /* If the domain reduction ended up empty. */

	    WasReduction = !IRIT_APX_EQ(OrigTMin, TMin) ||
			   !IRIT_APX_EQ(OrigTMax, TMax);
	}

	if (WasReduction) {
	    MVs1 = (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));

	    for (i = 0; i < NumOfMVs; i++) {
	        MVs1[i] = MvarMVRegionFromMV(MVs[i], TMin, TMax, l);
		if (MvarZeroMVConstraintFail(MVs1[i], Constraints[i])) {
		    for ( ; i >= 0; i--)
		        MvarMVFree(MVs1[i]);

		    IritFree(MVs1);
		    return NULL;
		}
	    }

	    PtList1 = MvarZeroMVsSubdiv(MVs1, Constraints, NumOfMVs,
					NumOfZeroMVs, ApplyNormalConeTest,
					SubdivTol, Depth + 1, SameSpace,
					ParamPerturb);

	    for (i = 0; i < NumOfMVs; i++)
	        MvarMVFree(MVs1[i]);

	    IritFree(MVs1);

	    return PtList1;
	}
	else {
	    CagdBType
		CanPurgeDomain1 = FALSE,
		CanPurgeDomain2 = FALSE;

	    if (ParamPerturb < (TMax - TMin) / 10.0)
	        t = (TMin + TMax) * 0.5 + ParamPerturb;
	    else
	        t = (TMin + TMax) * 0.5;

	    /* Lets see if we have a B-spline multivariate with interior    */
	    /* knot in this direction.  If so pick up interior knot instead.*/
	    if (MVAR_IS_BSPLINE_MV(MVs[0])) {
	        for (i = 0; i < NumOfMVs; i++) {
		    if (MVs[i] -> Lengths[l] != MVs[i] -> Orders[l]) {
		        CagdRType
			    r = MVs[i] -> KnotVectors[l][
						   (MVs[i] -> Lengths[l] +
						    MVs[i] -> Orders[l]) >> 1];

			if (r > TMin + IRIT_EPS && r < TMax - IRIT_EPS) {
			    t = r;
			    break;
			}
		    }
		}
	    }

	    /* Ensure we have t within domain, so subdivision can be used. */
	    assert(t >= TMin && t <= TMax);

	    MVs1 = (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));
	    MVs2 = (MvarMVStruct **) IritMalloc(NumOfMVs *
						      sizeof(MvarMVStruct *));

	    for (i = 0; i < NumOfMVs; i++) {
	        if (CanPurgeDomain1) {
		    MVs1[i] = NULL;
		    MVs2[i] = MvarMVSubdivAtParamOneSide(MVs[i], t, l, FALSE);
		    CanPurgeDomain2 |= MvarZeroMVConstraintFail(MVs2[i],
							      Constraints[i]);
		}
		else if (CanPurgeDomain2) {
		    MVs1[i] = MvarMVSubdivAtParamOneSide(MVs[i], t, l, TRUE);
		    MVs2[i] = NULL;
		    CanPurgeDomain1 |= MvarZeroMVConstraintFail(MVs1[i],
							      Constraints[i]);
		}
		else {
		    MVs1[i] = MvarMVSubdivAtParam(MVs[i], t, l);
		    MVs2[i] = MVs1[i] -> Pnext;
		    MVs1[i] -> Pnext = NULL;
		    CanPurgeDomain1 |= MvarZeroMVConstraintFail(MVs1[i],
							      Constraints[i]);
		    CanPurgeDomain2 |= MvarZeroMVConstraintFail(MVs2[i],
							      Constraints[i]);
		}

		if (CanPurgeDomain1 && CanPurgeDomain2) {
		    /* Both domains are purged away - abort here. */
		    for ( ; i >= 0; i--) {
		        MvarMVFree(MVs1[i]);
			MvarMVFree(MVs2[i]);
		    }
		    IritFree(MVs1);
		    IritFree(MVs2);

		    return NULL;
		}
	    }

	    if (!CanPurgeDomain1)
	        PtList1 = MvarZeroMVsSubdiv(MVs1, Constraints, NumOfMVs,
					    NumOfZeroMVs, ApplyNormalConeTest,
					    SubdivTol, Depth + 1, SameSpace, 
					    ParamPerturb);
	    else
	        PtList1 = NULL;

	    if (!CanPurgeDomain2)
	        PtList2 = MvarZeroMVsSubdiv(MVs2, Constraints, NumOfMVs,
					    NumOfZeroMVs, ApplyNormalConeTest,
					    SubdivTol, Depth + 1, SameSpace, 
					    ParamPerturb);
	    else
	        PtList2 = NULL;

	    for (i = 0; i < NumOfMVs; i++) {
	        MvarMVFree(MVs1[i]);
		MvarMVFree(MVs2[i]);
	    }
	    IritFree(MVs1);
	    IritFree(MVs2);

	    return (MvarPtStruct *) CagdListAppend(PtList1, PtList2);
	}
    }

#   ifdef DEBUG
    {
        IRIT_SET_IF_DEBUG_ON_PARAMETER(_DebugPrintZeroPts, FALSE) {
	    CagdRType Min, Max;

	    IRIT_INFO_MSG("Zero by subdivision level at");
	    for (i = 0; i < Dim; i++) {
	        MvarGetSubdivParamDomains(MVs[0], &Min, &Max, i);
		IRIT_INFO_MSG_PRINTF(" %f [%f %f]",
				     (Min + Max) * 0.5, Min, Max);
	    }
	    IRIT_INFO_MSG("\n");
	}
    }
#   endif /* DEBUG */

    /* If we got here then in all parametric directions, the domain is       */
    /* smaller than the SubdivTol.  Return central location of this patch.   */
    return MvarZeroGenPtMidMvar(MVs[0], FALSE);
}


/*****************************************************************************
* DESCRIPTION:                                                               *
*   Construct a point of the dimension as the given MV in the middle of its  *
* parametric domain.                                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   MV:         To construct a point in the middle of its domain.            *
*   SingleSol:  If TRUE, this point is a single solution it MV domain.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *:  The construct point in the middle of MV.                *
*****************************************************************************/
static MvarPtStruct *MvarZeroGenPtMidMvar(const MvarMVStruct *MV,
					  int SingleSol)
{
    int l;
    MvarPtStruct
        *Pt = MvarPtNew(MV -> Dim);

    for (l = 0; l < MV -> Dim; l++) {
	CagdRType Min, Max;

	MvarGetSubdivParamDomains(MV, &Min, &Max, l);
	Pt -> Pt[l] = (Min + Max) * 0.5;
    }

#   ifdef DEBUG_DUMP_DOMAINS
    fprintf(GlblDumpDomainsFile, "[OBJECT [RGB \"255,100,100\"] [Gray 0.5] [Width 0.02] NONE [CTLPT E%d", MV -> Dim);
    for (l = 0; l < MV -> Dim; l++)
	fprintf(GlblDumpDomainsFile, " %f", Pt -> Pt[l]);
    fprintf(GlblDumpDomainsFile, "]\n]\n");
#   endif /* DEBUG_DUMP_DOMAINS */

    AttrSetIntAttrib(&Pt -> Attr, "SingleSol", SingleSol);

    return Pt;
}
