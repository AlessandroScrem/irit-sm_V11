/******************************************************************************
* MvarStPl.c - Handling Stewart Platform constraints.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, May 14.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/geom_lib.h"

static MvarMVStruct *MvarStewartPlatPtFixedPtDistCnst(const MvarMVStruct *MVX,
						      const MvarMVStruct *MVY,
						      const MvarMVStruct *MVZ,
						      const CagdPType Pt2);
static MvarMVStruct *MvarStewartPlatformPtPtDistCnst(const MvarMVStruct *MV1X,
						     const MvarMVStruct *MV1Y,
						     const MvarMVStruct *MV1Z,
						     const MvarMVStruct *MV2X,
						     const MvarMVStruct *MV2Y,
						     const MvarMVStruct *MV2Z,
						     CagdRType EdgeLength);
static MvarMVStruct *MvarStewartPlatform2PtPtDistCnst(const MvarMVStruct *MV1,
						      const MvarMVStruct *MV2,
						      CagdRType EdgeLength);
static CagdCrvStruct *MvarST2BuildCircle(const CagdPType Pt1,
					 const CagdPType Pt2,
					 CagdRType D1,
					 CagdRType D2,
					 const CagdRType *OrientDir,
					 CagdBType Rational);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Build a quadratic constraint of a distance square between unknown point  *
* and a fixed point Pt2.   Distance is for now assumed zero.                 *
*                                                                            *
* PARAMETERS:                                                                *
*   MVX1, MVY1, MVZ1:   The X, Y, Z coefficients of the 1st unknown point.   *
*   Pt2:	     The second fixed location.                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:   The constructed constraint.                            *
*****************************************************************************/
static MvarMVStruct *MvarStewartPlatPtFixedPtDistCnst(const MvarMVStruct *MVX1,
						      const MvarMVStruct *MVY1,
						      const MvarMVStruct *MVZ1,
						      const CagdPType Pt2)
{
    CagdRType R;
    MvarMVStruct *MVTransX, *MVTransY, *MVTransZ,
        *MVTmp1, *MVTmp2, *MVTmp3, *MVTmp4, *MVTmp5;

    /* Translate the first point by -Pt2. */
    R = -Pt2[0];
    MVTransX = MvarMVCopy(MVX1);
    MvarMVTransform(MVTransX, &R, 1.0);

    R = -Pt2[1];
    MVTransY = MvarMVCopy(MVY1);
    MvarMVTransform(MVTransY, &R, 1.0);

    R = -Pt2[2];
    MVTransZ = MvarMVCopy(MVZ1);
    MvarMVTransform(MVTransZ, &R, 1.0);

    /* Compute L2 norm of the distance. */
    MVTmp1 = MvarMVMult(MVTransX, MVTransX);
    MVTmp2 = MvarMVMult(MVTransY, MVTransY);
    MVTmp3 = MvarMVMult(MVTransZ, MVTransZ);
    MVTmp4 = MvarMVAdd(MVTmp1, MVTmp2);
    MVTmp5 = MvarMVAdd(MVTmp3, MVTmp4);

    MvarMVFree(MVTransX);
    MvarMVFree(MVTransY);
    MvarMVFree(MVTransZ);
    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);
    MvarMVFree(MVTmp3);
    MvarMVFree(MVTmp4);

    return MVTmp5;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Build a quadratic constraint of a distance square between two unknown    *
* points.   Distance is given in EdgeLength.			             *
*                                                                            *
* PARAMETERS:                                                                *
*   MVX1, MVY1, MVZ1:   The X, Y, Z coefficients of the 1st unknown point.   *
*   MVX2, MVY2, MVZ2:   The X, Y, Z coefficients of the 2nd unknown point.   *
*   EdgeLength:         The desired distance between the two points.         *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:     The constructed constraint.                          *
*****************************************************************************/
static MvarMVStruct *MvarStewartPlatformPtPtDistCnst(const MvarMVStruct *MV1X,
						     const MvarMVStruct *MV1Y,
						     const MvarMVStruct *MV1Z,
						     const MvarMVStruct *MV2X,
						     const MvarMVStruct *MV2Y,
						     const MvarMVStruct *MV2Z,
						     CagdRType EdgeLength)
{
    CagdRType R;
    MvarMVStruct *MVTmp1, *MVTmp2, *MVTmp3, *MVTmp4,
                 *MVTmp5, *MVTmp6, *MVTmp7, *MVTmp8;
    
    /* Compute L2 norm of the distance. */
    MVTmp1 = MvarMVSub(MV1X, MV2X);
    MVTmp2 = MvarMVSub(MV1Y, MV2Y);
    MVTmp3 = MvarMVSub(MV1Z, MV2Z);
    MVTmp4 = MvarMVMult(MVTmp1, MVTmp1);
    MVTmp5 = MvarMVMult(MVTmp2, MVTmp2);
    MVTmp6 = MvarMVMult(MVTmp3, MVTmp3);
    MVTmp7 = MvarMVAdd(MVTmp4, MVTmp5);
    MVTmp8 = MvarMVAdd(MVTmp6, MVTmp7);

    /* Subtract the constrainted distance square. */
    R = -IRIT_SQR(EdgeLength);
    MvarMVTransform(MVTmp8, &R, 1.0);

    MvarMVFree(MVTmp1);
    MvarMVFree(MVTmp2);
    MvarMVFree(MVTmp3);
    MvarMVFree(MVTmp4);
    MvarMVFree(MVTmp5);
    MvarMVFree(MVTmp6);
    MvarMVFree(MVTmp7);

    return MVTmp8;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Build a quadratic constraint of a distance square between two unknown    *
* points.   Distance is given in EdgeLength.			             *
*                                                                            *
* PARAMETERS:                                                                *
*   MV1, MV2:   The XYZ functions of the two point positions as MVs.         *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:     The constructed constraint.                          *
*****************************************************************************/
static MvarMVStruct *MvarStewartPlatform2PtPtDistCnst(const MvarMVStruct *MV1,
						      const MvarMVStruct *MV2,
						      CagdRType EdgeLength)
{
    CagdRType R;
    MvarMVStruct *MVTmp1, *MVTmp2;

    /* Compute L2 norm of the distance. */
    MVTmp1 = MvarMVSub(MV1, MV2);
    MVTmp2 = MvarMVDotProd(MVTmp1, MVTmp1);

    /* Subtract the constrained distance square. */
    R = -IRIT_SQR(EdgeLength);
    MvarMVTransform(MVTmp2, &R, 1.0);

    MvarMVFree(MVTmp1);

    return MVTmp2;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Derive nine constraints for a stewart platform, so that given the        M
* lengths of the six rods connection the bottom and top bases, we can solve  M
* the problem.							             M
*   The bottom base is specified and so are the lengths of the edges of the  M
* top base.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   BottomBase:   Three vertices of triangle forming the base, denoted Pi.   M
*   TopBaseEdgeLengths:  Three edge lengths between three vertices Qi, of    M
*                 top base.						     M
*   WorkDomain:   The Euclidean space working domain min./max. range.        M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct **:  Dynamically allocated vector of nine constraints, in   M
*                 nine unknowns (bottom base with given vertices Pi, top     M
*                 vertices with unknown vertices Qi) as follows:	     M
*                 || Q1 - P1 || = L1,					     V
*                 || Q1 - P2 || = L2,					     V
*                 || Q2 - P2 || = L3,					     V
*                 || Q2 - P3 || = L4,					     V
*                 || Q3 - P3 || = L5,					     V
*                 || Q3 - P1 || = L6,					     V
*                 || Q1 - Q2 || = D1,					     V
*                 || Q2 - Q3 || = D2,					     V
*                 || Q3 - Q1 || = D3,					     V
*                 where Li are the six BaseConnectLengths and Dj are lengths M
*                 between adjacent edges of the top, moving, base.	     M
*                   Note the Li are not given and so in this building        M
*                 process the are assumed all zero.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarStewartPlatformSolve                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarStewartPlatformGenEqns                                               M
*****************************************************************************/
MvarMVStruct **MvarStewartPlatformGenEqns(const CagdPType BottomBase[3],
					  const CagdRType TopBaseEdgeLengths[6],
					  const CagdPType WorkDomain[2])
{
    MvarMVStruct *MVQ1X, *MVQ1Y, *MVQ1Z,
                 *MVQ2X, *MVQ2Y, *MVQ2Z,
                 *MVQ3X, *MVQ3Y, *MVQ3Z,
        **AllMVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * 9);

    /* The nine constraints: */
    MVQ1X = MvarBuildParamMV(9, 0, WorkDomain[0][0], WorkDomain[1][0]);
    MVQ1Y = MvarBuildParamMV(9, 1, WorkDomain[0][1], WorkDomain[1][1]);
    MVQ1Z = MvarBuildParamMV(9, 2, WorkDomain[0][2], WorkDomain[1][2]);
    MVQ2X = MvarBuildParamMV(9, 3, WorkDomain[0][0], WorkDomain[1][0]);
    MVQ2Y = MvarBuildParamMV(9, 4, WorkDomain[0][1], WorkDomain[1][1]);
    MVQ2Z = MvarBuildParamMV(9, 5, WorkDomain[0][2], WorkDomain[1][2]);
    MVQ3X = MvarBuildParamMV(9, 6, WorkDomain[0][0], WorkDomain[1][0]);
    MVQ3Y = MvarBuildParamMV(9, 7, WorkDomain[0][1], WorkDomain[1][1]);
    MVQ3Z = MvarBuildParamMV(9, 8, WorkDomain[0][2], WorkDomain[1][2]);

    AllMVs[0] = MvarStewartPlatPtFixedPtDistCnst(MVQ1X, MVQ1Y, MVQ1Z,
						 BottomBase[0]);
    AllMVs[1] = MvarStewartPlatPtFixedPtDistCnst(MVQ1X, MVQ1Y, MVQ1Z,
						 BottomBase[1]);
    AllMVs[2] = MvarStewartPlatPtFixedPtDistCnst(MVQ2X, MVQ2Y, MVQ2Z,
						 BottomBase[1]);
    AllMVs[3] = MvarStewartPlatPtFixedPtDistCnst(MVQ2X, MVQ2Y, MVQ2Z,
						 BottomBase[2]);
    AllMVs[4] = MvarStewartPlatPtFixedPtDistCnst(MVQ3X, MVQ3Y, MVQ3Z,
						 BottomBase[2]);
    AllMVs[5] = MvarStewartPlatPtFixedPtDistCnst(MVQ3X, MVQ3Y, MVQ3Z,
						 BottomBase[0]);

    AllMVs[6] = MvarStewartPlatformPtPtDistCnst(MVQ1X, MVQ1Y, MVQ1Z,
					        MVQ2X, MVQ2Y, MVQ2Z,
					        TopBaseEdgeLengths[0]);
    AllMVs[7] = MvarStewartPlatformPtPtDistCnst(MVQ2X, MVQ2Y, MVQ2Z,
					        MVQ3X, MVQ3Y, MVQ3Z,
					        TopBaseEdgeLengths[1]);
    AllMVs[8] = MvarStewartPlatformPtPtDistCnst(MVQ3X, MVQ3Y, MVQ3Z,
					        MVQ1X, MVQ1Y, MVQ1Z,
					        TopBaseEdgeLengths[2]);

    MvarMVFree(MVQ1X);
    MvarMVFree(MVQ1Y);
    MvarMVFree(MVQ1Z);
    MvarMVFree(MVQ2X);
    MvarMVFree(MVQ2Y);
    MvarMVFree(MVQ2Z);
    MvarMVFree(MVQ3X);
    MvarMVFree(MVQ3Y);
    MvarMVFree(MVQ3Z);

    return AllMVs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Derive nine constraints for a stewart platform, so that given the        M
* lengths of the six rods connection the bottom and top bases, we can solve  M
* the problem.							             M
*   The bottom base is specified and so are the lengths of the edges of the  M
* top base.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   AllCnstrnts:  Nine constraints without base connecting length set.       M
*   BaseConnectLengths:  Six specific edge-lengths of rods connection bases. M
*   WorkDomain:   The Euclidean space working domain min./max. range.        M
*   SubdivTol:    Tolerance of the solution.  This tolerance is measured in  M
*		  the parametric space of the surfaces.			     M
*   NumericTol:   Numeric tolerance of a possible numeric improvement stage. M
*		  The numeric stage is employed if NumericTol < SubdivTol.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   Solution(s) of the Qi locations, if any.		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarStewartPlatformGenEqns                                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarStewartPlatformSolve                                                 M
*****************************************************************************/
MvarPtStruct *MvarStewartPlatformSolve(const MvarMVStruct **AllCnstrnts,
				       const CagdRType BaseConnectLengths[6],
				       const CagdPType WorkDomain[2],
				       CagdRType SubdivTol,
				       CagdRType NumericTol)
{
    int i;
    MvarMVStruct *AllMVs[9];
    MvarPtStruct *Sols, *Sol;

    for (i = 0; i < 9; i++)
        AllMVs[i] = MvarMVCopy(AllCnstrnts[i]);

    /* Update the specific rod lengths given. */
    for (i = 0; i < 6; i++) {
        CagdRType
	    R = -BaseConnectLengths[i];

        MvarMVTransform(AllMVs[i], &R, 1.0);
    }

    /* And solve the constraints. */
    Sols = MvarMVsZeros0D(AllMVs, NULL, 9, SubdivTol, NumericTol);

    /* Map the solution locations to Euclidean space. */
    for (Sol = Sols; Sol != NULL; Sol = Sol -> Pnext) {
        for (i = 0; i < 9; i++) {
	    int j = i % 3;

	    Sol -> Pt[i] = Sol -> Pt[i] * (WorkDomain[1][j] - WorkDomain[0][j])
							    + WorkDomain[0][j];
	}
    }

    for (i = 0; i < 9; i++) {
        MvarMVFree(AllMVs[i]);
    }

    return Sols;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Constructs the circle (arc) in the plane orthogonal to line Pt1Pt2 so    *
* that the circle distance to Pt1 is D1 and to Pt2 is D2.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt1, Pt2:   Two points defining the line orthogonal to circle plane.     *
*   D1, D2:     Distances of constructed circle to Pt1 and Pt2, respectively.*
*   OrientDir:  If NULL full circles are derived, seeking all solutions      *
*                  including inverted and singular ones.		     *
*               If not NULL, assume direction to orient half-circle along,   *
*                  seeking only one consistent solution.		     *
*   Rational:   TRUE to represent the circles precisely as rational,         *
*               FALSE to approximate the circles as polynomials.             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:  Constructed circle or NULL if error.                   *
*****************************************************************************/
static CagdCrvStruct *MvarST2BuildCircle(const CagdPType Pt1,
					 const CagdPType Pt2,
					 CagdRType D1,
					 CagdRType D2,
					 const CagdRType *OrientDir,
					 CagdBType Rational)
{
    CagdRType Alpha, Height, Radius,
        D = IRIT_PT_PT_DIST(Pt1, Pt2);
    IrtVecType Dir;
    CagdPtStruct Center;
    IrtHmgnMatType Mat;
    CagdCrvStruct *Circ, *TCirc;

    /* Build the circle first in canonical position assuming Pt1 at the     */
    /* origin and Pt2 along the +Z axis.                                    */

    /* Use the cosine theorem to compute the angle in front of D2. */
    Alpha = (IRIT_SQR(D2) - IRIT_SQR(D1) - IRIT_SQR(D)) / (-2.0 * D1 * D);
    if (Alpha <= -1 || Alpha >= 1.0)
        return NULL;
    Alpha = acos(Alpha);
    Height = D1 * cos(Alpha);
    Radius = D1 * sin(Alpha);

    Center.Pt[0] = Center.Pt[1] = 0.0;
    Center.Pt[2] = Height;
    if (Rational)
        Circ = BspCrvCreateCircle(&Center, Radius);
    else
        Circ = BspCrvCreatePCircle(&Center, Radius);

    /* Rotate circle to the line Pt1Pt2 and translate the origin to Pt1. */
    IRIT_PT_SUB(Dir, Pt2, Pt1);
    if (OrientDir != NULL)
        GMGenTransMatrixZ2Dir2(Mat, Pt1, Dir, OrientDir, 1.0);
    else
        GMGenTransMatrixZ2Dir(Mat, Pt1, Dir, 1.0);

    TCirc = CagdCrvMatTransform(Circ, Mat);
    CagdCrvFree(Circ);

    return TCirc;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A second solver for Stewart platform that is optimized.  The original    M
* constraints are mapped to find a triangle (of TopEdgeLengths edges) whose  M
* vertices resides on three computed circles Ci, as follows,                 M
* 1. Compute circles Ci, i = 0, 1, 2, as the locus of points of distance     M
*    BotTopEdgeLengths[2i], and BotTopEdgeLengths[2i+1] from bottom base     M
*    points BottomBasePoints[i], BottomBasePoints[(i+1)%3], respectively.    M
* 2. Solve for the three constraints and three unknowns 		     M
*    || Ci - C[(i+1)%3] || = TopEdgeLengths[i].				     M
*                                                                            *
* PARAMETERS:                                                                M
*   BottomBasePoints:  The three points of the bottom base.		     M
*   BotTopEdgeLengths: The six edge lengths of rods connection bottom and    M
*                 top bases.						     M
*   TopEdgeLengths:  The three edge lengths of the top base.		     M
*   FullCircs:    TRUE to derive full circles, seeking inverted and singular M
*                 solutions.  FALSE to seek only one consistent solution.    M
*   Rational:     TRUE to represent the circles precisely as rational,       M
*                 FALSE to approximate the circles as polynomials.           M
*   SubdivTol:    Tolerance of the solution.  This tolerance is measured in  M
*		  the parametric space of the surfaces.			     M
*   NumericTol:   Numeric tolerance of a possible numeric improvement stage. M
*		  The numeric stage is employed if NumericTol < SubdivTol.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:  Computed solution locations.		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarStewartPlatformSolve, MvarStewartPlatformGenEqns                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarStewartPlatform2Solve                                                M
*****************************************************************************/
MvarPtStruct *MvarStewartPlatform2Solve(const CagdPType BottomBasePoints[3],
				        const CagdRType BotTopEdgeLengths[6],
				        const CagdRType TopEdgeLengths[3],
				        CagdBType FullCircs,
				        CagdBType Rational,
				        CagdRType SubdivTol,
				        CagdRType NumericTol)
{
    int i;
    CagdCrvStruct
        *Crv1 = MvarST2BuildCircle(BottomBasePoints[0], BottomBasePoints[1],
				   BotTopEdgeLengths[0], BotTopEdgeLengths[1],
				   NULL, Rational),
        *Crv2 = MvarST2BuildCircle(BottomBasePoints[1], BottomBasePoints[2],
				   BotTopEdgeLengths[2], BotTopEdgeLengths[3],
				   NULL, Rational),
        *Crv3 = MvarST2BuildCircle(BottomBasePoints[2], BottomBasePoints[0],
				   BotTopEdgeLengths[4], BotTopEdgeLengths[5],
				   NULL, Rational);
    MvarPtStruct *Sols, *Sol, *E3Sols, *E3Sol;
    MvarMVStruct *MV1, *MV2, *MV3, *MVTemp1, *MVs[3];

    if (Crv1 == NULL || Crv2 == NULL || Crv3 == NULL) {
        if (Crv1 != NULL)
	    CagdCrvFree(Crv1);
        if (Crv2 != NULL)
	    CagdCrvFree(Crv2);
        if (Crv3 != NULL)
	    CagdCrvFree(Crv3);

        return NULL;
    }
    CagdCrvSetDomain(Crv1, 0.0, 1.0);
    CagdCrvSetDomain(Crv2, 0.0, 1.0);
    CagdCrvSetDomain(Crv3, 0.0, 1.0);

    /* Now convert to constraints: || Crvi - Crvj || = TopEdgeLengths[i]. */
    MVTemp1 = MvarCrvToMV(Crv1);
    MV1 = MvarPromoteMVToMV2(MVTemp1, 3, 0);
    MvarMVFree(MVTemp1);
    MVTemp1 = MvarCrvToMV(Crv2);
    MV2 = MvarPromoteMVToMV2(MVTemp1, 3, 1);
    MvarMVFree(MVTemp1);
    MVTemp1 = MvarCrvToMV(Crv3);
    MV3 = MvarPromoteMVToMV2(MVTemp1, 3, 2);
    MvarMVFree(MVTemp1);

    MVs[0] = MvarStewartPlatform2PtPtDistCnst(MV1, MV2, TopEdgeLengths[0]);
    MVs[1] = MvarStewartPlatform2PtPtDistCnst(MV2, MV3, TopEdgeLengths[1]);
    MVs[2] = MvarStewartPlatform2PtPtDistCnst(MV3, MV1, TopEdgeLengths[2]);

    MvarMVFree(MV1);
    MvarMVFree(MV2);
    MvarMVFree(MV3);

    /* And solve the constraints. */
    Sols = MvarMVsZeros0D(MVs, NULL, 3, SubdivTol, NumericTol);

    /* Map the solution to Euclidean space. */
    for (Sol = Sols, E3Sols = NULL; Sol != NULL; Sol = Sol -> Pnext) {
        CagdRType *R;

	E3Sol = MvarPtNew(MVAR_PT_E9_TYPE);

	R = CagdCrvEval(Crv1, Sol -> Pt[0]);
	CagdCoerceToE3(&E3Sol -> Pt[0], &R, -1, Crv1 -> PType);
	R = CagdCrvEval(Crv2, Sol -> Pt[1]);
	CagdCoerceToE3(&E3Sol -> Pt[3], &R, -1, Crv2 -> PType);
	R = CagdCrvEval(Crv3, Sol -> Pt[2]);
	CagdCoerceToE3(&E3Sol -> Pt[6], &R, -1, Crv3 -> PType);

	IRIT_LIST_PUSH(E3Sol, E3Sols);
    }
    MvarPtFreeList(Sols);

    for (i = 0; i < 3; i++) {
        MvarMVFree(MVs[i]);
    }

    CagdCrvFree(Crv1);
    CagdCrvFree(Crv2);
    CagdCrvFree(Crv3);

    return E3Sols;
}
