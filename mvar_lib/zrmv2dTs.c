/******************************************************************************
* ZrMV2DTs.c - Tools to compute zero sets of multivariates when the problem   *
*	     has MVs representation and the expected solution  set is a       *
*            two-manifold: surface(s). This file handles the tessellation     *
*	     in the topologically guaranteed domains, from the known boundary *
*	     loops.							      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, June 2014.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/geom_lib.h"

#define MVAR_ZRMV2D_TR_EDGE_FACTOR 20.0   /* Triangle edge/subdiv tol ratio. */
#define MVAR_ZRMV2D_ORIENT_DET_TOL 1e-15 /* Smaller can't decide orientation.*/
#define MVAR_ZRMV2D_PURGE_CLOSE_PTS_TOL 1e-13      /* Minimal triangle edge. */
#define MVAR_ZRMV2D_DBG_TR_ELEV           /* Count back projection failures. */
#define MVAR_ZRMV2D_EXT_RATIO 2.0 /* For back projection MVs extension size. */   // yoni clean this
#define MVAR_ZRMV2D_NUMERIC_DIFF_STEP 1e-6

static void MvarZeroFilterLoopForTr(MvarPolylineStruct *Poly,
				    CagdRType MaxLen);
static void MvarZeroPurgeNonTrPts(MvarPolylineStruct *Poly);
static MvarPolylineStruct *MvarZeroPolyProjPlanarIn3D(
                                                const MvarPolylineStruct *Poly,
                                                const int *Coords);
static void MvarZeroValInDmn(const CagdRType *MinDmn,
			     const CagdRType *MaxDmn, 
			     CagdRType *Val, 
			     MvarMVDirType Dir, 
			     CagdRType Tol);
static CagdRType MvarZeroCheckTrOrientation(MvarTriangleStruct *Tr,
					    MvarZeroPrblmStruct *Problem);
static MvarPtStruct *MvarZeroEmbedElevatedPt(MvarZeroPrblmStruct *Problem,
				             MvarPtStruct *Pt,
					     CagdRType RestrictVal1, 
					     CagdRType RestrictVal2,
					     int *ProjCoords);
static MvarPtStruct *MvarZeroSearchElevPolyline(
					      MvarPolylineStruct *ElevPolyline,
					      CagdRType RestrictVal1, 
					      CagdRType RestrictVal2, 
					      int *ProjCoords);
static void MvarZeroBackProjMVs(MvarMVStruct **SubDmnMVs, 
				MvarMVStruct **NormalSpProblemMVs,
				int NumOfMVs,
				const CagdRType *MinDmn, 
				const CagdRType *MaxDmn, 
				CagdRType *RestrictVal1, 
				CagdRType *RestrictVal2, 
				int *ProjCoords, 
				CagdRType RestrictTol);
static void MvarZero2DSolOrientation(MvarTriangleStruct *ElevatedTr, 
				     MvarZeroPrblmStruct *Problem);
static void MvarZeroPreTrPolylineSol(MvarZeroSolutionStruct *Sol, 
				     MvarZeroPrblmStruct *Problem);
static CagdRType MvarZeroDistToNextTrPt(
		                  const MvarPtStruct ** const NextSegmentStart,
			          const MvarPtStruct *SegmentStart);
static void MvarZeroTagSegmentForTr(MvarPtStruct *SegmentStart,
				    MvarPtStruct *SegmentEnd, 
				    CagdRType StepLenInSegment);
static IPObjectStruct *MvarZeroTesselateLoopProj(
                                              MvarPolylineStruct *BoundaryPoly,
                                              CagdRType MaxEdgeLen, 
                                              const int *ProjCoords);
static CagdBType MvarZeroIsPtInDmn(const MvarPtStruct *Pt,
				   const CagdRType *MinDmn, 
				   const CagdRType *MaxDmn, 
				   const int *SkipDirs,
				   const CagdRType Tol);
static void MapPllnSol2EuclidSp(MvarZeroSolutionStruct *Sol, 
		                const MvarZeroPrblmStruct *Problem);
static MvarVecStruct *MvarZeroNormalVecAtPt(
                                           const MvarPtStruct *ParamSpPt, 
					   const MvarZeroPrblmStruct *Problem);
static MvarPolylineStruct *MapPl2EuclidSp(MvarPolylineStruct *Poly, 
					  const MvarZeroPrblmStruct *Problem);
static MvarVecStruct *MvarZeroParam2EuclidTangent(
                                           const MvarPtStruct *Pt,
                                           const MvarVecStruct *ParamTanVec, 
                                           const MvarZeroPrblmStruct *Problem);
static IrtGnrlMatType MvarZeroEvalNumericDiffMat(
                                             const MvarPtStruct *Pt, 
                                             MvarMapPrm2EucCallBackFuncType F, 
                                             CagdRType *MinDmn,
                                             CagdRType *MaxDmn, 
                                             int InDim, 
                                             int OutDim,
					     CagdBType *Success);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   In the 2D solution case, we convert the solution from a set of closed    *
* loops to a set of triangles. This means that the actual "numeric step"     *
* in fact occurs here, and not deeper down the subdivision tree.	     *
*									     *
* PARAMETERS:                                                                *
*   Sol:	The solution to be organized.				     *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct:  The new, organized solution.		     *
*****************************************************************************/
MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol2DMVs(
                                                  MvarZeroSolutionStruct *Sol,
                                                  MvarZeroPrblmStruct *Problem)
{
    CagdBType VertexElevFailed, Compute3DNormals, 
	MapParam2Euclid = Problem -> _Internal -> 
	                 CallbackFunctions.MapPtParamSpace2EuclidSpace != NULL;
    int i, j, *ProjCoords, VNum, PtInd, SolDim,
	ElevationFails = 0,
	NumOfMVs = Problem -> NumOfConstraints,
	Dim = Problem -> U.MVs[0] -> Dim;
    MvarConstraintType
	*Constraints = IritMalloc(NumOfMVs * sizeof(MvarConstraintType)); 
    CagdRType RestrictVal1, RestrictVal2, *Epsilons, 
	*MinDmn = (CagdRType *) IritMalloc(Dim * sizeof(CagdRType)),
	*MaxDmn = (CagdRType *) IritMalloc(Dim * sizeof(CagdRType)),
	NumericTol = Problem -> NumericTol,
	SubdivTol = Problem -> SubdivTol;
    MvarPtStruct
	*NewPt = NULL,
	**ElevPts = NULL, 
	**EuclidPts = NULL, 
	**SolPts = NULL;
    MvarVecStruct
	**VrtxNormalsR3 = NULL;
    MvarZeroSubDmnInfoStruct *ProjInfo;
    MvarPolylineStruct *PolyIter, *DummyPoly;
    MvarTriangleStruct *NewTr,
	*ElevatedTr = NULL;
    MvarMVStruct **SubDmnMVs, **NormalSpProblemMVs, **ExtendedDmnMVs;
    MvarZeroPrblmStruct *NormalSpProblem;
    MvarZeroSolutionStruct *NormalSpSolution, *OutTrSol;
    IPObjectStruct *IPObjPlaneTrs;

    Compute3DNormals = MapParam2Euclid || Dim == 3;

    if (Sol == NULL) {
	IritFree(Constraints);
	IritFree(MinDmn);
	IritFree(MaxDmn);
	return Sol;
    }
    else {
	MvarZeroHandleTJunctions(Sol -> U.Pl, Sol -> TJList);
	MvarZeroTJFreeList(Sol -> TJList);
	Sol -> TJList = NULL;
	MvarZeroPreTrPolylineSol(Sol, Problem);      /* Filtering pts only. */
    }

    if (_MVGlblZeroOutputTypeLoops) {             /* Return polyline loops. */
	IritFree(Constraints);
	IritFree(MinDmn);
	IritFree(MaxDmn);
	PolyIter = Sol -> U.Pl;
	while (PolyIter != NULL) {
	    ProjInfo = AttrGetRefPtrAttrib(PolyIter -> Attr, "_TrInfo");
	    MvarSubDmnInfoStructFree(ProjInfo, NumOfMVs);
	    PolyIter = PolyIter -> Pnext;
	}
	if (MapParam2Euclid)
	    MapPllnSol2EuclidSp(Sol, Problem);
	return Sol;
    }

    /* Triangles output required. */
    NormalSpProblemMVs = (MvarMVStruct **) 
	                         IritMalloc(NumOfMVs * sizeof(MvarMVStruct *));
    ExtendedDmnMVs = (MvarMVStruct **) 
	                         IritMalloc(sizeof(MvarMVStruct *) * NumOfMVs);
    IRIT_GEN_COPY(Constraints, Problem -> Constraints, 
	          NumOfMVs * sizeof(MvarConstraintType));

    Epsilons = IritMalloc(sizeof(CagdRType) * NumOfMVs);
    for (j = 0; j < NumOfMVs; j++) 
	Epsilons[j] = MVAR_ZRMV2D_EXT_RATIO * SubdivTol;

    PolyIter = Sol -> U.Pl;
    while (PolyIter != NULL) {
	IPVertexStruct *Cur2DVertex;
	IPPolyVrtxIdxStruct 
	    *UniqueVertexHelper = NULL;

	ProjInfo = AttrGetRefPtrAttrib(PolyIter -> Attr, "_TrInfo");
	SubDmnMVs = ProjInfo -> MVs;
	ProjCoords = ProjInfo -> ProjDirs;
	MvarMVDomain(SubDmnMVs[0], MinDmn, MaxDmn, -1);

	/* Project the polyline on the 2D plane, and triangulate. */
	IPObjPlaneTrs = MvarZeroTesselateLoopProj(PolyIter, 
	                                          MVAR_ZRMV2D_TR_EDGE_FACTOR *
	                                          Problem -> SubdivTol, 
	                                          ProjCoords);
	if (IPObjPlaneTrs -> U.Pl == NULL) {
	    /* Planar triangulation failed. Might occur in singular domains. */
	    MvarSubDmnInfoStructFree(ProjInfo, NumOfMVs);
	    DummyPoly = PolyIter;
	    PolyIter = PolyIter -> Pnext;
	    MvarPolylineFree(DummyPoly);
	    continue;
	}

	/* Convert to a unique vertex structure: */
	UniqueVertexHelper = IPCnvPolyToPolyVrtxIdxStruct(IPObjPlaneTrs, 0, 0);

	VNum = UniqueVertexHelper -> NumVrtcs;
	ElevPts = (MvarPtStruct **) IritMalloc(
	                                      (VNum) * sizeof(MvarPtStruct *));
	if (Compute3DNormals)
	    VrtxNormalsR3 = (MvarVecStruct **) IritMalloc(
	                                     (VNum) * sizeof(MvarVecStruct *));


	/* Elevate each vertex using the unique vertex array. */
	for (i = 0; i < VNum; i++) {
	    VertexElevFailed = FALSE;
	    Cur2DVertex = UniqueVertexHelper -> Vertices[i];
	    RestrictVal1 = Cur2DVertex -> Coord[0];
	    RestrictVal2 = Cur2DVertex -> Coord[1];

	    /* If this is a boundary point, no eq' solving needed, search    */
	    /* elevation point in the non-projected loop:                    */
	    NewPt = MvarZeroSearchElevPolyline(PolyIter, 
 		                               RestrictVal1, 
		                               RestrictVal2, 
		                               ProjCoords);

	    if (NewPt == NULL) { /* Back project by eq' solving. */
		MvarZeroBackProjMVs(SubDmnMVs, NormalSpProblemMVs,
		                    NumOfMVs, MinDmn, MaxDmn, 
		                    &RestrictVal1, &RestrictVal2, 
		                    ProjCoords, Problem -> SubdivTol);

		/* Domain extension, to increase chance of elevation       */
		/* success. Then create the problem and solve. 	           */
		for (j = 0; j < NumOfMVs; j++) {
		    ExtendedDmnMVs[j] = MvarMVExtension(NormalSpProblemMVs[j],
							NULL, NULL, Epsilons);
		}
		NormalSpProblem = MvarZeroSolverPrblmNew(
		                                 (const MvarMVStruct * const *)
		                                                ExtendedDmnMVs,
		                                  NULL, NumOfMVs, Constraints, 
		                                  NumericTol, SubdivTol * 1e-2,
		                                  IRIT_INFNTY, FALSE);
		NormalSpSolution = MvarZeroSolver(NormalSpProblem);
		if (NormalSpSolution == NULL) {
		    ElevationFails++;
		    VertexElevFailed = TRUE;
		}
		else { /* The back projection returned candidates. */
		    NewPt = MvarZeroEmbedElevatedPt(Problem,
			                            NormalSpSolution -> U.Pt,
			                            RestrictVal1, 
			                            RestrictVal2,
			                            ProjCoords);
		}
		MvarZeroSolverPrblmFree(NormalSpProblem);
		MvarZeroSolverSolutionFree(NormalSpSolution, TRUE);
		for (j = 0; j < NumOfMVs; j++) {
		    MvarMVFree(NormalSpProblemMVs[j]);
		    MvarMVFree(ExtendedDmnMVs[j]);
		}
	    }
	    ElevPts[i] = NewPt;
	}

	if (MapParam2Euclid) {
	    /* Each solution point is mapped to the 3D Euclidean space,      */
	    /* using the function provided by the problem:                   */
	    const CagdRType *Tmp;

	    EuclidPts = (MvarPtStruct **) IritMalloc(
		                              (VNum) * sizeof(MvarPtStruct *));
	    for (i = 0; i < VNum; i++) {
		EuclidPts[i] = NULL;
		if (ElevPts[i] != NULL) {
		    Tmp = MVAR_ZERO_SLVR_APPLY(MapPtParamSpace2EuclidSpace)
			                               (ElevPts[i] -> Pt, Dim);
		    EuclidPts[i] = MvarPtNew(3);
		    IRIT_GEN_COPY(EuclidPts[i] -> Pt, Tmp, 
			          3 * sizeof(CagdRType));
		}
		else
		    EuclidPts[i] = NULL;
	    } 
	}

	/* Compute the unit normal to the solution surface at the 3D         */
	/* Euclidean (or parameter) space points.                            */
	if (Compute3DNormals) {
	    for (i = 0; i < VNum; i++) {
		VrtxNormalsR3[i] = MvarZeroNormalVecAtPt(ElevPts[i], Problem);
	    }
	}

	/* Reconstruct the mesh structure from the vertex structure in 2D:   */
	/* use the index info in the 2D struct BUT OPERATE on the required   */
	/* arrays of Mvar points and normals, creating Mvar triangles.       */
	SolPts = MapParam2Euclid ? EuclidPts : ElevPts;
	SolDim = MapParam2Euclid ? 3 : Dim;

	for (i = 0; i < UniqueVertexHelper -> NumPlys; i++) {
	    NewTr = MvarTriangleNew(SolDim);
	    for (j = 0; j < 3; j++) {
		PtInd = UniqueVertexHelper -> Polygons[i][j];
		if (SolPts[PtInd] != NULL) {
		    IRIT_GEN_COPY(NewTr -> Vrtcs[j], 
			          SolPts[PtInd] -> Pt, 
			          SolDim * sizeof(CagdRType));
		}
		else {
		    MvarTriangleFree(NewTr);
		    NewTr = NULL;
		    break;
		}
		if (Compute3DNormals) {
		    /* A 3D normal vector for each point is available. */
		    if (VrtxNormalsR3[PtInd] != NULL) {
			IRIT_GEN_COPY(NewTr -> Nrmls[j], 
			              VrtxNormalsR3[PtInd] -> Vec, 
			              3 * sizeof(CagdRType));
		    }
		}
		else { /* Assign zero as the normal. */
		    IRIT_ZAP_MEM(NewTr -> Nrmls[j], Dim * sizeof(CagdRType));
		}
	    }
	    if (NewTr != NULL) {
		IRIT_LIST_PUSH(NewTr, ElevatedTr);
	    }
	}	
	MvarSubDmnInfoStructFree(ProjInfo, NumOfMVs);
	DummyPoly = PolyIter;
	PolyIter = PolyIter -> Pnext;
	MvarPolylineFree(DummyPoly);
	IPFreeObject(IPObjPlaneTrs);
	IPPolyVrtxIdxFree(UniqueVertexHelper);
	for (i = 0; i < VNum; i++) {
	    if (ElevPts[i] != NULL)
		MvarPtFree(ElevPts[i]);
	}
	IritFree(ElevPts);
	if (MapParam2Euclid) {
	    for (i = 0; i < VNum; i++) {
		if (EuclidPts[i] != NULL)
		    MvarPtFree(EuclidPts[i]);
	    } 
	  
	    IritFree(EuclidPts);
	}
	if (Compute3DNormals) {
	    for (i = 0; i < VNum; i++) {
		if (VrtxNormalsR3[i] != NULL)
		    MvarVecFree(VrtxNormalsR3[i]);
	    } 
	    IritFree(VrtxNormalsR3);
	}
    }
    IritFree(Epsilons);
    IritFree(Constraints);
    IritFree(NormalSpProblemMVs);
    IritFree(ExtendedDmnMVs);
    IritFree(MinDmn);
    IritFree(MaxDmn);
    MvarZeroSolverSolutionFree(Sol, FALSE); /* Polylines were already freed. */

    MvarZero2DSolOrientation(ElevatedTr, Problem);

#   ifdef MVAR_ZRMV2D_DBG_TR_ELEV 
    printf("Failed to elevate %d vertices.\n", ElevationFails);
#   endif /* MVAR_ZRMV2D_DBG_TR_ELEV */

    OutTrSol = MvarZeroSolverSolutionNew(ElevatedTr, NULL, NULL, 
				         MVAR_ZER_SLVR_SOLUTION_TR_LIST);
    return OutTrSol;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Tags point on a polyline loop with an attribute so they shall later      *
* be used as triangle vertices. Points are chosen according to their         *
* (cord) distance from the last chosen point, and in a manner that is        *
* orientation invariant.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly:		The polyline loop.			             *
*   MaxLen:	        The maximal edge length of the result (between       *
*			tags).						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.							 	     *
*****************************************************************************/
static void MvarZeroFilterLoopForTr(MvarPolylineStruct *Poly,
				    CagdRType MaxLen)
{
    CagdRType StepsInCurSegment, StepLenInSegment, CurSegmentLen;
    MvarPtStruct 
	*SegmentEnd = Poly -> Pl,
	*SegmentStart = Poly -> Pl;
    
    /* The current segment is the part of the polyline from the current      */
    /* "must-choose" pt to the next "must-choose" pt:                        */
    while (SegmentEnd -> Pnext != NULL) {
	CurSegmentLen = MvarZeroDistToNextTrPt((const MvarPtStruct ** const)
	                                                         &SegmentEnd,
					       SegmentStart);
	StepsInCurSegment = floor(CurSegmentLen / MaxLen) + 1.0;
	StepLenInSegment = CurSegmentLen / StepsInCurSegment;
	MvarZeroTagSegmentForTr(SegmentStart, SegmentEnd, 
	                        StepLenInSegment);
	SegmentStart = SegmentEnd;
    }
    /* Tag the last segment end (not tagged in the loop!) */
    AttrSetIntAttrib(&SegmentEnd -> Attr, "TRPt", TRUE);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Purging the point objects that are not tagged as triangulation points in *
* a polyline. If two triangulation points are too close, one of them is      *
* purged.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly:	    The polyline.					     *
*									     *
* RETURN VALUE:                                                              *
*   void	     							     *
*****************************************************************************/
static void MvarZeroPurgeNonTrPts(MvarPolylineStruct *Poly)
{
    CagdBType PurgeCurrent;
    MvarPtStruct 
	*PtIter = Poly -> Pl,
	*PrevPtIter = NULL;

    while (PtIter != NULL) { /* First pass - filter non-triangulation pts. */
	PurgeCurrent = AttrGetIntAttrib(PtIter -> Attr, "TRPt") != TRUE;
	if (PurgeCurrent) {
	    /* The first and the last should always be triangulation points. */
	    assert(PrevPtIter != NULL && PtIter -> Pnext != NULL);
	    PrevPtIter -> Pnext = PtIter -> Pnext;
	    MvarPtFree(PtIter);
	    PtIter = PrevPtIter -> Pnext;
	}
	else { /* No purge, just advance. */
	    PrevPtIter = PtIter;
	    PtIter = PtIter -> Pnext;
	}
    }

    PtIter = Poly -> Pl;
    PrevPtIter = NULL;
    while (PtIter != NULL) { /* Second pass - filter too close. */
	if (PrevPtIter != NULL && 
	    MvarPtDistTwoPoints(PrevPtIter, PtIter) < 
	    MVAR_ZRMV2D_PURGE_CLOSE_PTS_TOL) {
	    PrevPtIter -> Pnext = PtIter -> Pnext;
	    MvarPtFree(PtIter);
	    PtIter = PrevPtIter -> Pnext;
	}
	else {
	    PrevPtIter = PtIter;
	    PtIter = PtIter -> Pnext;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a polyline in R^Dim, projects it on the two required dimensions    *
* and assigns a third coordinate, that is always zero. The result is the     *
* planar projection, embedded in R^3, contained in the x,y plane.            *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly:   The polyline to project.	    				     *
*   Coords: The required two coordinates to project on.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPolylineStruct *: The new, planar polyline in R^3.                   *
*****************************************************************************/
static MvarPolylineStruct *MvarZeroPolyProjPlanarIn3D(
                                                const MvarPolylineStruct *Poly,
				                const int *Coords)
{
    CagdBType 
	IsFirst = TRUE;
    const MvarPtStruct *PtIter;
    MvarPtStruct *NewPt, *TmpPt, *PtList;

    if (Poly == NULL)
	return NULL;

    PtIter = Poly -> Pl;
    while (PtIter != NULL) {
	/* Create a 3D point with a zero 3rd coord', as expected by the      */
	/* triangulation routines.					     */
	NewPt = MvarPtNew(3);
	NewPt -> Pt[0] = PtIter -> Pt[Coords[0]];
	NewPt -> Pt[1] = PtIter -> Pt[Coords[1]];
	NewPt -> Pt[2] = 0.0;
	if (IsFirst) {
	    PtList = TmpPt = NewPt;
	    IsFirst = FALSE;
	}
	else {
	    TmpPt -> Pnext = NewPt;
	    TmpPt = NewPt;
	} 
	PtIter = PtIter -> Pnext;
    }
    return MvarPolylineNew(PtList);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test if a real value is in the given domain in a prescribed direction,   *
* and trim it to the domain boundary if it is not. The violation is assumed  *
* to be at most Tol.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   MinDmn:	    The min' values of the domain in each direction.	     *
*   MaxDmn:	    The max' values of the domain in each direction.	     *
*   Val:            The  value to check.				     *
*   Dir:	    The tested direction.			             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroValInDmn(const CagdRType *MinDmn,
			     const CagdRType *MaxDmn, 
			     CagdRType *Val, 
			     MvarMVDirType Dir, 
			     CagdRType Tol)
{
    if (*Val < MinDmn[Dir]) {
	if (*Val < MinDmn[Dir] - Tol) {
#	    ifdef MVAR_ZRMV2D_DBG_MV_RESTRICTION
		printf("Restriction value is too far out of domain.\n");
#	    endif /* MVAR_ZRMV2D_DBG_MV_RESTRICTION */

	    assert(0);
	}
	else
	    *Val = MinDmn[Dir];
    }
    if (*Val > MaxDmn[Dir]) {
	if (*Val > MaxDmn[Dir] + Tol) {
#	    ifdef MVAR_ZRMV2D_DBG_MV_RESTRICTION
		printf("Restriction value is too far out of domain.\n");
#	    endif /* MVAR_ZRMV2D_DBG_MV_RESTRICTION */

	    assert(0);
	}
	else
	    *Val = MaxDmn[Dir];
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test if orientation of a solution triangle is consistent with the	     *
* problem, by the convention that the triangle edges v1-v0, v2-v0 (in that   *
* order!) form a positive basis (determinant) for R^n together with the n-2  *
* gradients of the problem, when evaluated at v0.		             *
*                                                                            *
* PARAMETERS:                                                                *
*   Tr:		The triangle to test.					     *
*   Problem:	The original problem, used for gradient evaluation.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType:	The value of the respective determinant. If positive, the    *
*		the orientation is correct. If negative - should be flipped. *
*		If zero - either the gradients are linearly dependent or     *
*		the triangle is degenerate.				     *
*****************************************************************************/
static CagdRType MvarZeroCheckTrOrientation(MvarTriangleStruct *Tr,
				            MvarZeroPrblmStruct *Problem)
{
    int i, j,
	Dim = Tr -> Dim;         /* The problem's Dim may be different now. */
    CagdRType DetVal, *TmpGradVal;
    IrtGnrlMatType
        RNBasisMat = (IrtGnrlMatType) IritMalloc(IRIT_SQR(Dim) *
						            sizeof(IrtRType));

    /* First two rows of the orientation matrix are the triangle edges: */
    for (i = 0; i < 2; i++) {
	for (j = 0; j < Dim; j++) {
	    RNBasisMat[Dim * i + j] =
	                            Tr -> Vrtcs[i + 1][j] - Tr -> Vrtcs[0][j];
	}
    }

    /* The next rows are either the correctly oriented unit normal at       */
    /* vertex 0, or the Dim - 2 gradients evaluated at vertex 0, if Dim > 3.*/
    for (i = 2; i < Dim; i++) {
	if (Dim == 3) {
	    CAGD_GEN_COPY(&RNBasisMat[i * Dim], Tr -> Nrmls[0] , 
		          sizeof(CagdRType) * Dim);
	}
	else {
	    TmpGradVal = MvarMVEvalGradient(
				    Problem -> _Internal -> MVGradients[i - 2],
				    Tr -> Vrtcs[0], 0);
	    CAGD_GEN_COPY(&RNBasisMat[i * Dim], TmpGradVal, 
			  sizeof(CagdRType) * Dim);
	}
    }
    DetVal = MatGnrlDetMatrix(RNBasisMat, Dim);
    IritFree(RNBasisMat);
    return DetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test if there's an n-2 dimensional face of the n dimensional domain,     *
* that the point belongs to. For example - if n==3, this tests if the point  *
* belongs to an edge of the box  (in other words: are there at least two     *
* faces of dimension n-1 that the point contained in?)			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt:	    The point to test.						     *
*   MinDmn: The min' values of the domain in each direction.		     *
*   MaxDmn: The max' values of the domain in each direction.		     *
*   Tol:    The numeric tolerance to test against.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType: TRUE if there's such n-2 dim' face, FALSE otherwise.          *
*****************************************************************************/
CagdBType MvarZeroIsPtOnEdge(const MvarPtStruct *Pt,
			     const CagdRType *MinDmn,
			     const CagdRType *MaxDmn,
			     CagdRType Tol)
{
    int i,
	Dim = Pt -> Dim,
	FoundBoundaryCords = 0;
    CagdRType CurVal;

    for (i = 0; i < Dim; i++) {
	CurVal = Pt -> Pt[i];
	if (IRIT_APX_EQ_EPS(CurVal, MinDmn[i], Tol) || 
	    IRIT_APX_EQ_EPS(CurVal, MaxDmn[i], Tol)) {
	    FoundBoundaryCords++;
	}
	if (FoundBoundaryCords == 2)
	    return TRUE;
    }
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Inserting two additional coordinates at prescribed directions and values *
* to a given point.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	    The underlying equation solving problem.		     *
*   Pt:		    The point to be added two dimensions.		     *
*   RestrictVal1:   First new value.					     *
*   RestrictVal2:   Second new value.					     *
*   ProjCoords:	    The two new coordinate directions in terms of output pt. *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: The new point, two dimensions higher.		     *
*****************************************************************************/
static MvarPtStruct *MvarZeroEmbedElevatedPt(MvarZeroPrblmStruct *Problem,
				             MvarPtStruct *Pt,
					     CagdRType RestrictVal1, 
					     CagdRType RestrictVal2,
					     int *ProjCoords)

{
    int k, j,
	IndShift = 0,
	Dim = Pt -> Dim + 2;
    MvarPtStruct 
	*NewPt = MvarPtNew(Dim);

    for (k = 0; k < Dim; k++) {
	if (k == ProjCoords[0]) {
	    NewPt -> Pt[k] = RestrictVal1;
	    IndShift = 1;
	}
	else if (k == ProjCoords[1]) {
	    NewPt -> Pt[k] = RestrictVal2;
	    IndShift = 2;
	}
	else {
	    NewPt -> Pt[k] = Pt -> Pt[k - IndShift];
	}
    }

    /* If candidate is out of the original domain due extension, trim it.   */
    for (j = 0; j < Dim; j++) {
	MvarZeroValInDmn(Problem -> _Internal -> OrigMVMinDmn, 
		         Problem -> _Internal -> OrigMVMaxDmn, 
			 &(NewPt -> Pt[j]), j, Problem -> SubdivTol);
    }
    return NewPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Searching a polyline for a point that has two prescribed values at two   *
* given coordinate directions. The polyline is assumed to have at most one   *
* point with the above property, and if found - a copy is returned.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   ElevPolyline:	The polyline to search.				     *
*   RestrictVal1:	First required value.				     *
*   RestrictVal2:	Second required value.				     *
*   ProjCoords:		The coordinate directions of interest to search.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarPtStruct *: The point copy if found, or NULL otherwise. 	     *
*****************************************************************************/
static MvarPtStruct *MvarZeroSearchElevPolyline(
					      MvarPolylineStruct *ElevPolyline,
					      CagdRType RestrictVal1, 
					      CagdRType RestrictVal2, 
					      int *ProjCoords)
{
    MvarPtStruct
	*NewPt = NULL,
	*PtIter = ElevPolyline -> Pl;

    while (PtIter != NULL) {
	if (IRIT_APX_EQ_EPS(RestrictVal1, PtIter -> Pt[ProjCoords[0]], 
	                    MVAR_ZRMV2D_PTS_EQ_TOL) &&
	    IRIT_APX_EQ_EPS(RestrictVal2, PtIter -> Pt[ProjCoords[1]], 
	                    MVAR_ZRMV2D_PTS_EQ_TOL)) {
	    NewPt = MvarPtCopy(PtIter);
	    break;
	}
	PtIter = PtIter -> Pnext; 
    }
    return NewPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Restricting multivariates from Dim variables to Dim - 2 variables, by    *
* fixing two prescribed coordinates to given constant values.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   SubDmnMVs:		The input multivariates.			     *
*   NormalSpProblemMVs: The output multivariates (pointers preallocated).    *
*   NumOfMVs:		The number of input/output multivariates.	     *
*   MinDmn, MaxDmn:	The domain against which to validate the restriction *
*		        values, or trim them accordingly.		     *
*   RestrictVal1:	First required value, (may be trimmed to domain).    *
*   RestrictVal2:	Second required value, (may be trimmed to domain).   *
*   ProjCoords:		The coordinate directions to be set constant.        *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.							 	     *
*****************************************************************************/
static void MvarZeroBackProjMVs(MvarMVStruct **SubDmnMVs, 
				MvarMVStruct **NormalSpProblemMVs,
				int NumOfMVs,
			        const CagdRType *MinDmn, 
				const CagdRType *MaxDmn, 
				CagdRType *RestrictVal1, 
				CagdRType *RestrictVal2, 
				int *ProjCoords, 
				CagdRType RestrictTol)
{
    int j;
    MvarMVStruct *TmpMV;

    for (j = 0; j < NumOfMVs; j++) {
	/* Verify that small numeric changes in the planar triangulation did */
	/* not result in a location out of the current sub-domain.           */
	MvarZeroValInDmn(MinDmn, MaxDmn, RestrictVal1, ProjCoords[0], 
	                 RestrictTol);
	MvarZeroValInDmn(MinDmn, MaxDmn, RestrictVal2, ProjCoords[1], 
	                 RestrictTol);

	/* Restrict the original MV to an MV of two variables less. */
	TmpMV = MvarMVFromMV(SubDmnMVs[j], *RestrictVal1, ProjCoords[0]);
	NormalSpProblemMVs[j] = MvarMVFromMV(TmpMV, *RestrictVal2, 
					     ProjCoords[1] - 1);
	MvarMVFree(TmpMV);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Checking and correcting the orientation of a list of triangles, such     *
* that the triangle edges form a positively oriented basis with the grads    *
* of the problem (to which the triangles approximate the solution).	     *
*                                                                            *
* PARAMETERS:                                                                *
*   ElevatedTr:		The triangles list.			             *
*   Problem:	        The underlying equation solving problem.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.							 	     *
*****************************************************************************/
static void MvarZero2DSolOrientation(MvarTriangleStruct *ElevatedTr, 
				     MvarZeroPrblmStruct *Problem)
{
    int TrCnt = 0; 
#   ifdef MVAR_ZRMV2D_ORIENT_DET
        int SmallDetValCnt = 0;
#   endif /* MVAR_ZRMV2D_ORIENT_DET */
    CagdRType OrientationDetVal;
    MvarTriangleStruct 
	*TrIter = ElevatedTr;

    while (TrIter != NULL) {
	TrCnt++;
	OrientationDetVal = MvarZeroCheckTrOrientation(TrIter, Problem);
	if (OrientationDetVal < -MVAR_ZRMV2D_ORIENT_DET_TOL) { 
	    /* Flip the orientation: */
	    int l;

	    for (l = 0; l < TrIter -> Dim; l++) {
		IRIT_SWAP(CagdRType, TrIter -> Vrtcs[1][l],
			  TrIter -> Vrtcs[2][l]);
		IRIT_SWAP(CagdRType, TrIter -> Nrmls[1][l], 
			  TrIter -> Nrmls[2][l]);
	    }
	}
#	ifdef MVAR_ZRMV2D_ORIENT_DET
	    else if (OrientationDetVal < MVAR_ZRMV2D_ORIENT_DET_TOL) {
		SmallDetValCnt++;
	    }
#	endif /* MVAR_ZRMV2D_ORIENT_DET */

	TrIter = TrIter -> Pnext;
    }

#ifdef MVAR_ZRMV2D_ORIENT_DET
    printf("OrientationDet Val too small in %d out of %d triangles.\n\n", 
           SmallDetValCnt, TrCnt);
#endif /* MVAR_ZRMV2D_ORIENT_DET */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Preprocessing of a polyline loop, as a preparation step before it is     *
* triangulated later.				                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol:		The solution, still in polylines form.               *
*   Problem:	        The underlying equation solving problem.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.							 	     *
*****************************************************************************/
static void MvarZeroPreTrPolylineSol(MvarZeroSolutionStruct *Sol, 
				     MvarZeroPrblmStruct *Problem)
{
    MvarPolylineStruct
	*PolyIter = Sol -> U.Pl;

    while (PolyIter != NULL) { 
	MvarZeroFilterLoopForTr(PolyIter, 
		                MVAR_ZRMV2D_TR_EDGE_FACTOR * 
			        Problem -> SubdivTol);
	MvarZeroPurgeNonTrPts(PolyIter);
	PolyIter = PolyIter -> Pnext;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a polyline loop and a 2D plane into which it can be projected in a *
* one-to-one manner, returns a list of 2D triangles that are a tessellation  *
* of the 2D projected loop.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   BoundaryPoly:	The input, not yet projected, polyline.              *
*   MaxEdgeLen:	        An upper bound on the planar triangulation edges.    *
*   ProjCoords:		The two coordinate direction to project to.          *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *: The object containing the list of 2D triangles.        *
*****************************************************************************/
static IPObjectStruct *MvarZeroTesselateLoopProj(
                                              MvarPolylineStruct *BoundaryPoly,
                                              CagdRType MaxEdgeLen, 
                                              const int *ProjCoords)
{
    MvarPolylineStruct *ProjPoly;
    IPObjectStruct *TmpObjPoly, *ObjTrList;
    IPPolygonStruct *RefineObjTrList, *Plgns,
	*NewPlgns = NULL;

    /* Project the polyline on the 2D plane, and triangulate. */
    ProjPoly = MvarZeroPolyProjPlanarIn3D(BoundaryPoly, ProjCoords);
    TmpObjPoly = MvarCnvrtMVPolysToIritPolys2(ProjPoly, TRUE);

    /* Verify we only have three vertices or more. */
    for (Plgns = TmpObjPoly -> U.Pl; Plgns != NULL; ) {
	IPPolygonStruct *Plgn;

	IRIT_LIST_POP(Plgn, Plgns);
	if (IPVrtxListLen(Plgn -> PVertex) < 3) {
	    IPFreePolygon(Plgn);
	}
	else {
	    IRIT_LIST_PUSH(Plgn, NewPlgns);
	}
    }
    TmpObjPoly -> U.Pl = IPReversePlList(NewPlgns);

    ObjTrList = GMConvertPolysToTrianglesIntrrPt(TmpObjPoly);
    RefineObjTrList = GMLimitTrianglesEdgeLen(ObjTrList -> U.Pl, MaxEdgeLen);
    IPFreePolygonList(ObjTrList -> U.Pl);
    ObjTrList -> U.Pl = RefineObjTrList;

    IPFreeObject(TmpObjPoly);
    MvarPolylineFree(ProjPoly);

    return ObjTrList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find the next point along a polyline that must be chosen for             *
* triangulation, and measure the distance to that point along the polyline.  *
*                                                                            *
* PARAMETERS:                                                                *
*   NextSegmentStart:	The output pointer to the next must-choose point.    *
*   SegmentStart:	The current must-choose point to begin the search.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType: The total segment length.			             *
*****************************************************************************/
static CagdRType MvarZeroDistToNextTrPt(
			          const MvarPtStruct ** const NextSegmentStart,
				  const MvarPtStruct *SegmentStart)
{
    CagdBType MustChoose;
    CagdRType
	TotalDist = 0.0;
   
    if (SegmentStart -> Pnext == NULL) {
	/* Should not happen, but in case we are invoked with the last point */
	/* in the polyline:						     */
	return TotalDist;
    }

    else { /* Advance until encountering a "must-choose" point: */
	*NextSegmentStart = SegmentStart;
	MustChoose = FALSE;
	while (!MustChoose) {
	    /* Add the distance first: */
	    TotalDist += MvarPtDistTwoPoints(*NextSegmentStart, 
		                            (*NextSegmentStart) -> Pnext);
	    /* Advance and check if the next is a must-choose: */
	    *NextSegmentStart = (*NextSegmentStart) -> Pnext;
	    MustChoose = (*NextSegmentStart) -> Pnext == NULL || 
		         AttrGetIntAttrib((*NextSegmentStart) -> Attr, 
		                          "TJunc") == TRUE || 
			 AttrGetIntAttrib((*NextSegmentStart) -> Attr, 
					  "Corner") == TRUE;
	}
	return TotalDist;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Tags a subset of points from a polyline segment, as sparse as possible   *
* without exceeding the prescribed length between tagged points.             *
*                                                                            *
* PARAMETERS:                                                                *
*   SegmentStart:	The first point of the given segment.                *
*   SegmentEnd:	        The last point of the given segment.                 *
*   StepLenInSegment:	The maximal distance that should not be exceeded     *
*                       between two consecutive tagged points.		     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.                                                                    *
*****************************************************************************/
static void MvarZeroTagSegmentForTr(MvarPtStruct *SegmentStart,
				    MvarPtStruct *SegmentEnd, 
			            CagdRType StepLenInSegment)
{
    CagdRType DeltaDist, Ratio, CordLen,
	Dist = 0.0;
    MvarPtStruct
	*PtIter = SegmentStart;

    /* The first is tagged: */
    AttrSetIntAttrib(&PtIter -> Attr, "TRPt", TRUE); 

    while (PtIter != SegmentEnd) {
	/* Tag if traveled far enough. */
	Dist += CordLen = MvarPtDistTwoPoints(PtIter, PtIter -> Pnext);
	DeltaDist = Dist - StepLenInSegment;
	if (DeltaDist >= 0.0) {
	    /* Tagged pt should be the one closer to the exact distance: */
	    Ratio = DeltaDist / CordLen;
	    if (Ratio <= 0.5) {
		/* Tag the next pt, and update dist forward (positive). */
		AttrSetIntAttrib(&PtIter -> Pnext -> Attr, "TRPt", TRUE); 
		Dist = DeltaDist;
	    }
	    else { /* Ratio > 0.5 */
		/* Tag the current pt, and update dist back (negative). */
		AttrSetIntAttrib(&PtIter -> Attr, "TRPt", TRUE); 
		Dist = DeltaDist - CordLen;
	    }
	}
	PtIter = PtIter -> Pnext;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Post process of the solution in the polyline loop representation case,   *
* when a problem-specific mapping is required, from the parameter space      *
* where the equations are solved, to the Euclidean space where the equation  *
* originated from with the proper geometric semantics. The mapping is done   *
* in place.                                                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol:	        The solution in the parameter space, to be mapped.   *
*   Problem:	        The original equation solving problem.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.								     *
*****************************************************************************/
static void MapPllnSol2EuclidSp(MvarZeroSolutionStruct *Sol, 
			        const MvarZeroPrblmStruct *Problem)
{
    MvarPolylineStruct *PolyIter, *NewPoly, *ProjPolyList, *ProjPolyIter;

    assert(_MVGlblZeroOutputTypeLoops && 
	   Sol -> ActiveRepresentation == MVAR_ZER_SLVR_SOLUTION_POLYLINE);
    
    PolyIter = Sol -> U.Pl;
    ProjPolyIter = ProjPolyList = MapPl2EuclidSp(PolyIter, Problem);

    while (PolyIter != NULL) {
	NewPoly = MapPl2EuclidSp(PolyIter, Problem);
	ProjPolyIter -> Pnext = NewPoly;
	ProjPolyIter = ProjPolyIter -> Pnext;
	PolyIter = PolyIter -> Pnext;
    }
    ProjPolyIter -> Pnext = NULL;

    MvarPolylineFreeList(Sol -> U.Pl);
    Sol -> U.Pl = ProjPolyList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Mapping a single polyline from the parameter space of the problem to the *
* Euclidean space, using a problem-specific callback function, which maps    *
* a single point.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   Poly:	The solution polyline in the parameter space, to be mapped.  *
*   Problem:	The original equation solving problem.                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.								     *
*****************************************************************************/
static MvarPolylineStruct *MapPl2EuclidSp(MvarPolylineStruct *Poly, 
					  const MvarZeroPrblmStruct *Problem)
{
    int Dim = Problem -> U.MVs[0] -> Dim;
    const CagdRType *Tmp;
    MvarPtStruct *NewPt, *ProjPtIter,
	*PtIter = Poly -> Pl, 
	*NewPtList = MvarPtNew(3);

    Tmp = MVAR_ZERO_SLVR_APPLY(MapPtParamSpace2EuclidSpace)(PtIter -> Pt, Dim);
    IRIT_GEN_COPY(NewPtList -> Pt, Tmp, 3 * sizeof(CagdRType));
    ProjPtIter = NewPtList;
    PtIter = PtIter -> Pnext;

    while (PtIter != NULL) {
	NewPt = MvarPtNew(3);
	Tmp = MVAR_ZERO_SLVR_APPLY(MapPtParamSpace2EuclidSpace)
	                                                  (PtIter -> Pt, Dim);
	IRIT_GEN_COPY(NewPt -> Pt, Tmp, 3 * sizeof(CagdRType));
	ProjPtIter -> Pnext = NewPt;
	ProjPtIter = ProjPtIter -> Pnext;
	PtIter = PtIter -> Pnext;
    }
    ProjPtIter -> Pnext = NULL;

    return MvarPolylineNew(NewPtList);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a point on the solution manifold in the parameter space of the     *
* problem, obtain a 3D normal vector to the mapped solution surface in 3D    *
* Euclidean space, at the 3D image of the given point.                       *
*                                                                            *
* PARAMETERS:                                                                *
*   ParamSpPt:	    The location in the param' space on the solution.        *
*   Problem:	    The original equation solving problem.		     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarVecStruct *:	The computed 3D unit normal to the solution in  3D   *
*			Euclidean space, at the 3D image of Pt.              *
*****************************************************************************/
static MvarVecStruct *MvarZeroNormalVecAtPt(const MvarPtStruct *ParamSpPt, 
					    const MvarZeroPrblmStruct *Problem)
{
    CagdBType 
	MapParam2Euclid = Problem -> _Internal -> 
	                CallbackFunctions.MapPtParamSpace2EuclidSpace != NULL;
    int Dim = Problem -> U.MVs[0] -> Dim;
    CagdRType *TmpGradVal;
    MvarVecStruct *UnitNormal;
    
    if (ParamSpPt == NULL)
	return NULL;

    UnitNormal = MvarVecNew(3);

    if (MapParam2Euclid) {
	/* If we are here we cannot use the gradient(s) as the normal to the */
	/* to the solution surface at the given point. We must find two      */
	/* tangents, correctly oriented, map them to the tangents in the     */
	/* Euclidean space, and then take the cross product (this may be     */
	/* required with Dim == 3 sometimes as well).                        */

	int i;
	CagdRType Det;
	MvarVecStruct *EuclidTanVec1, *EuclidTanVec2, 
	    **TanVecs = (MvarVecStruct **) IritMalloc(
	                                          sizeof(MvarVecStruct *) * 2),
	    **GradVecs = (MvarVecStruct **) IritMalloc(
	                                  sizeof(MvarVecStruct *) * (Dim - 2));

	/* Step 1: compute Dim - 2 grads. */
	for (i = 0; i < Dim - 2; i++) {
	    TmpGradVal = MvarMVEvalGradient(
		                        Problem -> _Internal -> MVGradients[i],
		                        ParamSpPt -> Pt, 0);
	    GradVecs[i] = MvarVecNew(Dim);
	    CAGD_GEN_COPY(GradVecs[i] -> Vec, TmpGradVal, 
		          sizeof(CagdRType) * Dim); 
	}
	/* Step 2: obtain the grads ortho' complement. This pair spans the   */
	/*         tangent plane at the point, in the parameter space.       */
	if (MvarVecWedgeProd(GradVecs, Dim - 2, TanVecs, Dim, TRUE, &Det)) {
	    if (Det < 0.0) {
		IRIT_SWAP(MvarVecStruct *, TanVecs[0], TanVecs[1]);
	    }
	}
	else {
	    /* Failed to obtain tangents. */
	    for (i = 0; i < Dim - 2; i++)
		MvarVecFree(GradVecs[i]);
	    IritFree(GradVecs);
	    IritFree(TanVecs);
	    MvarVecFree(UnitNormal);
	    return NULL;
	}

	/* Step 3: Map the tangent pair to a 3D tangent pair. If F is the    */
	/*         mapping between the surfaces, then dF maps tangents at p  */
	/*         to tangents at F(p).                                      */
	EuclidTanVec1 = MvarZeroParam2EuclidTangent(ParamSpPt, TanVecs[0], 
	                                            Problem);
	EuclidTanVec2 = MvarZeroParam2EuclidTangent(ParamSpPt, TanVecs[1], 
	                                            Problem);

	/* Step 4: The Euclidean normal is the R^3 cross product of the      */
	/*         Euclidean tangents.                                       */
	IRIT_CROSS_PROD(UnitNormal -> Vec, EuclidTanVec1 -> Vec, 
	                EuclidTanVec2 -> Vec);
	if (!MvarVecNormalize(UnitNormal)) {
	    /* Close to zero cross product. Might be singular differential.  */
	    MvarVecFree(UnitNormal);
	    UnitNormal = NULL;
	}
	for (i = 0; i < 2; i++)
	    MvarVecFree(TanVecs[i]);
	IritFree(TanVecs);

	for (i = 0; i < Dim - 2; i++)
	    MvarVecFree(GradVecs[i]);
	IritFree(GradVecs);

	MvarVecFree(EuclidTanVec1);
	MvarVecFree(EuclidTanVec2);

	return UnitNormal;
    }
    else {
	assert(Dim == 3);
	/* The parameter space of the problem is where we are required to    */
	/* the solution at. In this case, the normal is the normalized grad. */
	/* This guarantees correct orientation of this normal vector. It     */
	/* does not guarantee the correct orientation of the triangle it     */
	/* it belongs to, which is handled elsewhere.                        */
	TmpGradVal = MvarMVEvalGradient(Problem -> _Internal -> MVGradients[0],
	                                ParamSpPt -> Pt, 0);
	CAGD_GEN_COPY(UnitNormal -> Vec, TmpGradVal, sizeof(CagdRType) * Dim);
	if (!MvarVecNormalize(UnitNormal)) {
	    /* Close to zero gradient. */
	    MvarVecFree(UnitNormal);
	    return NULL;
	}
    }
    return UnitNormal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given a tangent vector to the solution surface at a point in the         *
* parameter space of the problem, map it a tangent vector in the Euclidean   *
* space, at the Euclidean image of the point (on the Euclidean image of the  *
* solution surface). Using the mapping function F and the chain rule for     *
* differentiation, we simply implement the matrix multiplication of the      *
* differential matrix of F by the given tangent vector.                      *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt:		    The location in the param' space of the input tangent.   *
*   ParamTanVec:    The tangent vector in the parameter space.	             *
*   Problem:	    The original equation solving problem.		     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarVecStruct *:	The computed 3D tangent to the solution in  3D       *
*			Euclidean space, at the 3D image of Pt.              *
*****************************************************************************/
static MvarVecStruct *MvarZeroParam2EuclidTangent(
                                            const MvarPtStruct *Pt,
                                            const MvarVecStruct *ParamTanVec, 
					    const MvarZeroPrblmStruct *Problem)
{
    CagdBType DiffMatSuccess;
    int i, j, 
	Dim = Problem -> U.MVs[0] -> Dim;
    CagdRType 
	*MaxDmn = Problem -> _Internal -> OrigMVMaxDmn,
	*MinDmn = Problem -> _Internal -> OrigMVMinDmn;
    MvarMapPrm2EucCallBackFuncType
	F = Problem -> _Internal -> 
	                        CallbackFunctions.MapPtParamSpace2EuclidSpace;
    MvarVecStruct 
	*EuclidTanVec = MvarVecNew(3);
    IrtGnrlMatType 
	DFAtPt = MvarZeroEvalNumericDiffMat(Pt, F, MinDmn, MaxDmn, Dim, 3, 
	                                    &DiffMatSuccess);

    if (!DiffMatSuccess) {
	IritFree(DFAtPt);
	MvarVecFree(EuclidTanVec);
	return NULL;
    }
    /* Evaluate DFAtPt * ParamTanVec: */
    IRIT_ZAP_MEM(EuclidTanVec -> Vec, sizeof(CagdRType) * 3);
    for (i = 0; i < 3; i++) {
	for (j = 0; j < Dim; j++)
	    EuclidTanVec -> Vec[i] += 
	                          DFAtPt[i * Dim + j] * ParamTanVec -> Vec[j];
    }
    IritFree(DFAtPt);
    return EuclidTanVec;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Numeric approximation of the differential matrix of a given function at  *
* a given point. The central finite difference is used by default, while     *
* near the domain boundaries, the forward/backward difference is used.       *
*                                                                            *
* PARAMETERS:                                                                *
*   Pt:		The location to evaluate the differential at. Assumed in the *
*		domain.							     *
*   F:		The pointer to the function to be differentiated.	     *
*   MinDmn:	The min' edges of the domain of F.			     *
*   MaxDmn:	The max' edges of the domain of F.			     *
*   InDim:	The dimension of the domain of F.			     *
*   OutDim:	The dimension of the range of F.			     *
*   Success:	An output flag indicating successful or not computation.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtGnrlMatType:	The computed OutDim x InDim matrix (if successful).  *
*****************************************************************************/
static IrtGnrlMatType MvarZeroEvalNumericDiffMat(
                                             const MvarPtStruct *Pt, 
					     MvarMapPrm2EucCallBackFuncType F, 
					     CagdRType *MinDmn,
					     CagdRType *MaxDmn, 
					     int InDim, 
					     int OutDim, 
					     CagdBType *Success)
{
    CagdBType FWDNotInDmn, BWDNotInDmn;
    int i, j;
    CagdRType 
	*FWDPt = (CagdRType *) IritMalloc(sizeof(CagdRType) * InDim),
	*BWDPt = (CagdRType *) IritMalloc(sizeof(CagdRType) * InDim);
    const CagdRType
	*FAtFWDPt = NULL,
	*FAtBWDPt = NULL,
	*FAtPt = F(Pt -> Pt, InDim);
    MvarVecStruct
	*FValDiff = MvarVecNew(OutDim),
	*FValAtPt = MvarVecNew(OutDim),
	*FValAtFWDPt = MvarVecNew(OutDim),
	*FValAtBWDPt = MvarVecNew(OutDim);
    IrtGnrlMatType
	DFAtPt = (IrtGnrlMatType) IritMalloc(InDim * OutDim *sizeof(IrtRType));

    CAGD_GEN_COPY(FValAtPt -> Vec, FAtPt, OutDim * sizeof(CagdRType));
    *Success = TRUE;
    for (j = 0; j < InDim; j++) {
	/* Evaluating column j of the differential matrix at Pt. Find the    */
	/* numeric approximation of derivatives of F w.r.t. dimension j.     */
	CAGD_GEN_COPY(FWDPt, Pt -> Pt, InDim * sizeof(CagdRType));
	CAGD_GEN_COPY(BWDPt, Pt -> Pt, InDim * sizeof(CagdRType));
	FWDPt[j] += MVAR_ZRMV2D_NUMERIC_DIFF_STEP;
	BWDPt[j] -= MVAR_ZRMV2D_NUMERIC_DIFF_STEP;
	FWDNotInDmn = Pt -> Pt[j] + MVAR_ZRMV2D_NUMERIC_DIFF_STEP > MaxDmn[j];
	BWDNotInDmn = Pt -> Pt[j] - MVAR_ZRMV2D_NUMERIC_DIFF_STEP < MinDmn[j];

	if (!FWDNotInDmn && !BWDNotInDmn) {
	    /* Differentiation using central finite difference. */
	    FAtFWDPt = F(FWDPt, InDim);
	    CAGD_GEN_COPY(FValAtFWDPt -> Vec, FAtFWDPt, 
		          OutDim * sizeof(CagdRType));
	    FAtBWDPt = F(BWDPt, InDim);
	    CAGD_GEN_COPY(FValAtBWDPt -> Vec, FAtBWDPt, 
		          OutDim * sizeof(CagdRType));
	    MvarVecSub(FValDiff, FValAtFWDPt, FValAtBWDPt);
	    FValDiff = MvarVecScale(
		                  FValDiff, 
		                  1.0 / (2.0 * MVAR_ZRMV2D_NUMERIC_DIFF_STEP));
	}
	else if (BWDNotInDmn) {
	    /* Differentiation using forward finite difference. */
	    FAtFWDPt = F(FWDPt, InDim);
	    CAGD_GEN_COPY(FValAtFWDPt -> Vec, FAtFWDPt, 
		          OutDim * sizeof(CagdRType));
	    MvarVecSub(FValDiff, FValAtFWDPt, FValAtPt);
	    FValDiff = MvarVecScale(FValDiff, 
		                    1.0 / MVAR_ZRMV2D_NUMERIC_DIFF_STEP);
	}
	else if (FWDNotInDmn) {
	    /* Differentiation using backward finite difference. */
	    FAtBWDPt = F(BWDPt, InDim);
	    CAGD_GEN_COPY(FValAtBWDPt -> Vec, FAtBWDPt, 
		          OutDim * sizeof(CagdRType));
	    MvarVecSub(FValDiff, FValAtPt, FValAtBWDPt);
	    FValDiff = MvarVecScale(FValDiff, 
		                    1.0 / MVAR_ZRMV2D_NUMERIC_DIFF_STEP);
	}
	else {
	    /* Too small domain for this finite difference. */
	    MvarVecFree(FValDiff);
	    MvarVecFree(FValAtPt);
	    MvarVecFree(FValAtFWDPt);
	    MvarVecFree(FValAtBWDPt);
	    IritFree(DFAtPt);
	    IritFree(FWDPt);
	    IritFree(BWDPt);
	    *Success = FALSE;
	}
	for (i = 0; i < OutDim; i++)
	    DFAtPt[i * InDim + j] = FValDiff -> Vec[i];
    }
    MvarVecFree(FValDiff);
    MvarVecFree(FValAtPt);
    MvarVecFree(FValAtFWDPt);
    MvarVecFree(FValAtBWDPt);
    IritFree(FWDPt);
    IritFree(BWDPt);
    return DFAtPt;
}
