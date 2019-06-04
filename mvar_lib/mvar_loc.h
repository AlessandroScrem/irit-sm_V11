/******************************************************************************
* Mvar_loc.h - header file for the MVAR library.			      *
* This header is local to the library - it is NOT external.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, May. 97.					      *
******************************************************************************/

#ifndef MVAR_LOC_H
#define MVAR_LOC_H

#include <math.h>
#include <stdio.h>
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/mvar_lib.h"	    /* Include the external header as well. */

#define MVAR_BSCT_NUMER_TOL		1e-10

#define MVAR_NUMER_ZERO_NUM_STEPS	100
#define MVAR_ZERO_PARAM_REL_PERTURB	0.0000301060
#define MVAR_ZERO_PRECOND_SCALE		10
#define MVAR_MVZR1D_SUBD_PRTRB	        1.301060e-2
#define MVAR_ZRMV2D_PTS_EQ_TOL		1e-12
#define MVAR_ZRMV2D_PT_ON_EDGE_TOL	1e-12
#define MVAR_ZRMV2D_PT_IN_FACE_TOL	1e-12
#define MVAR_ZRMV2D_TJ_SEARCH_TOL       1e-12

/* #define DEBUG_MVAR_MVZR1D_DMN_TRACING */
/* #define DEBUG_MVAR_MVZR1D_LINK_NEIGH */

#define MVAR_ZERO_SLVR_APPLY(FunctionName) Problem -> \
    _Internal -> CallbackFunctions.FunctionName
#define MVAR_ZERO_SLVR_SOLUTION_SPACE_DIM(Problem) \
    IRIT_MAX(0, MvarZeroSolverGetDmnDim(Problem) - \
                     Problem -> _Internal -> NumOfZeroSubdivConstraints)
#define MVAR_EXPR_TREE_GET_COMMON_EXPR(ET) ( \
    assert(ET -> NodeType == MVAR_ET_NODE_COMMON_EXPR && \
	   ET -> Left != NULL),			 \
    ET -> Left)
#define MVAR_EVAL_NUMER_ERR_L1(Eqns, MVs, NumMVs, x) \
	(Eqns) != NULL ? MvarExprTreeEqnsEvalErrorL1((Eqns), (x)) : \
			 MvarMVEvalErrorL1((MVs), (x), NumMVs)

typedef enum {
    MVAR_BSCT_CV_CV = 1,
    MVAR_BSCT_CV_PT,
    MVAR_BSCT_NONE
} MvarBsctType;

typedef struct MvarVoronoiCrvStruct {
    struct MvarVoronoiCrvStruct *Pnext;
    MvarBsctType Type;
    CagdSrfStruct *F3;
    CagdCrvStruct *Crv1;
    CagdCrvStruct *Crv2; /* Crv2 will hold the second curve if Type is       */
                       /* CV_CV, and the rational bisector if Type is CV_PT. */
    CagdPType Pt;
} MvarVoronoiCrvStruct;

typedef struct MvarZeroTJunctionStruct {
    struct MvarZeroTJunctionStruct *Pnext;
    CagdBType IsHandled;             /* Already found a polyline and added. */
    MvarPtStruct *TJuncPrev;		 /* The point before the T-Junction */
    MvarPtStruct *TJunc;				  /* The T-Junction */
    MvarPtStruct *TJuncNext;		  /* The point after the T-Junction */
} MvarZeroTJunctionStruct;

/* Auxiliary structure for univariate solution spaces computation.  */
typedef struct MvarMVZR1DAuxStruct {           
    MvarVecStruct **OrthoBasis;
    MvarVecStruct *TempVec;
    MvarVecStruct *CorrVec;
    CagdRType *MinDmn;
    CagdRType *MaxDmn;
    CagdRType *MinDmn2;
    CagdRType *MaxDmn2;
    MvarVecStruct **GradVecs;
    MvarVecStruct *SITTempVec;/* SIT = StepInTangentVector+CorrStep related.*/
    MvarVecStruct *SITTanVec;
    CagdRType *A; 
    CagdRType *x;
    CagdRType *b;
    CagdRType *bCopy;
    MvarVecStruct **TempList;
    MvarConstraintType *Constraints;
    int NumOfMVs;
} MvarMVZR1DAuxStruct;

/* The zero set solver callback functions. These calls vary according to   */
/* the dimensionality of the solution space and the representation:        */
typedef struct MvarZeroSolverCallBackFcnStruct {
    void (*SolveBoundary)(MvarZeroPrblmStruct *Problem);    
    CagdBType (*NoZeroCrossTest)(MvarZeroPrblmStruct *Problem, 
				 int Ind);
    CagdBType (*GuaranteedTopologyTest)(MvarZeroPrblmStruct *Problem);
    MvarZeroPrblmStruct **(*SubdivProblem)(MvarZeroPrblmStruct *Problem);
    MvarZeroSolutionStruct *(*NumericImprovement)
				               (MvarZeroPrblmStruct *Problem);
    MvarZeroSolutionStruct *(*SingularSolution)
					       (MvarZeroPrblmStruct *Problem);
    MvarZeroSolutionStruct *(*UniteSolutions)(MvarZeroSolutionStruct *Sol1,
					      MvarZeroSolutionStruct *Sol2, 
					      MvarMVDirType Dir,
					      CagdRType Param,
					      MvarZeroPrblmStruct *Problem);
    MvarZeroSolutionStruct *(*OrganizeSolution)(MvarZeroSolutionStruct *Sol,
						MvarZeroPrblmStruct *Problem);
    void (*UpdateDmn)(MvarZeroPrblmStruct *Problem);
    CagdBType (*FirstSmoothUpdates)(MvarZeroPrblmStruct *Problem);
    MvarMapPrm2EucCallBackFuncType MapPtParamSpace2EuclidSpace;
} MvarZeroSolverCallBackFcnStruct;

/* The internal structure of the zero set solver: Private info for mvar_lib */
/* internal use.						            */
typedef struct MvarZeroPrblmIntrnlStruct {
    union {				
	MvarMVStruct **FirstSmoothMVs;	    /* Used in the 0d numeric step. */
	MvarExprTreeStruct **FirstSmoothETs;
    } U;
    MvarMVGradientStruct **MVGradients; 
    MvarZeroSolutionStruct *SolutionsOnBoundary;
    int NumOfBoundarySolutions;
    CagdBType UnderSubdivTol;
    CagdBType SingleSolutionTest0D;
    CagdBType NoLoopTest1D;
    CagdBType SingleComponentTest1D;
    CagdBType One2OneProjOnPlane2D;
    CagdBType HasC1Discont;
    CagdBType IsFirstSmooth;
    CagdBType PurgeByNumericStep;
    CagdBType ConstructionComplete;
    CagdBType ZeroMVsSameSpace;
    int C1DiscontInd;
    int ExpectedSolutionDim;	/* When the problem is regular. */
    int NumOfZeroSubdivConstraints;
    int SubdivDepth;
    int *DBGCntSinglrDmns;
    int *DBGCntGuaranteedDmns;
    MvarMVDirType One2OneProjDirs[2];
    CagdRType ParamPerturb;
    CagdRType *OrigMVMinDmn;	/* The domain of the original problem. */
    CagdRType *OrigMVMaxDmn;
    CagdRType *MVMinDmn;	/* The domain of the current problem. */
    CagdRType *MVMaxDmn;
    CagdRType MaxSide;
    MvarMVDirType MaxSideDir;
    MvarMVDirType SubdivDir;
    CagdRType ParamOfSubdiv;
    CagdRType C1DiscontParam;
    MvarZeroTJunctionStruct *TJList; /* Mesh T-Junctions, used in 2D solver. */
    MvarZeroSolverCallBackFcnStruct CallbackFunctions; 
} MvarZeroPrblmIntrnlStruct;

/* Additional info used in the 2D solution case. */
typedef struct MvarZeroSubDmnInfoStruct {
    MvarMVStruct **MVs; 		    /* Of the sub-domain's problem. */
    int ProjDirs[2];	/* Directions for which 1-1 projection is possible. */
} MvarZeroSubDmnInfoStruct;

/* Data structures for 2 contact motion analysis */

typedef struct Mvar2CtLineStruct {
    CagdPType P[2];
    CagdRType Epsilon;
} Mvar2CtLineStruct;

typedef struct Mvar2CtSphereStruct {
    CagdPType Center;
    CagdRType Radius;
} Mvar2CtSphereStruct;

typedef struct Mvar2CtAABBStruct {
    CagdRType Xmin, Xmax, Ymin, Ymax;
    CagdRType Cx, Cy, Radius;
} Mvar2CtAABBStruct;

typedef struct Mvar2CtConeStruct {
    CagdPType Axis;
    CagdRType Angle;
} Mvar2CtConeStruct;

typedef struct Mvar2CtCParamStruct {
    CagdRType TMin, TMax;
    CagdRType RMin, RMax;
    struct Mvar2CtCParamStruct *Left, *Right;
} Mvar2CtCParamStruct;

typedef struct Mvar2CtBVNodeStruct {
    int Convexity;
    Mvar2CtLineStruct LBV;	/* Bouding volume. */
    Mvar2CtSphereStruct SBV;	/* Bounding volume. */
    Mvar2CtConeStruct NCone;
    CagdRType KMin, KMax;	/* Curvature. */
    CagdRType Min, Max;		/* Domain. */
    CagdPType NMin, NMax;	/* Normal dir. */
    CagdRType Length;
    int Id;			/* Crv id. */
    struct Mvar2CtBVNodeStruct *Left, *Right;
} Mvar2CtBVNodeStruct;

typedef struct Mvar2CtBVHStruct {
    CagdCrvStruct *Crv;
    CagdCrvStruct *DCrv;
    CagdCrvStruct *DDCrv;
    CagdCrvStruct *NCrv;
    CagdCrvStruct *Curvature;

    MvarPtStruct *CurvContacts;
    Mvar2CtBVNodeStruct *Root;
    CagdRType Tol;
    CagdRType Radius;
} Mvar2CtBVHStruct;

IRIT_GLOBAL_DATA_HEADER CagdRType 
    MvarBsctSubdivTol, 
    MvarBsctNumerTol,
    MvarBsctUVTol;

/* Global vars defined in zero solver. */
IRIT_GLOBAL_DATA_HEADER CagdBType 
    _MVGlblZeroETUseCommonExpr,
    _MVGlblZeroETCnvrtBzrETs2MVs;

IRIT_GLOBAL_DATA_HEADER int
    _MVGlblZeroOutputTypeLoops,
    _MVGlblZeroApplyDomainReduction,
    _MVGlblZeroApplyGradPreconditioning,
    _MVGlblZeroApplyParallelHyperPlaneTest,
    _MVGlblZeroApplyNormalConeTest,
    _MVGlblZeroApplyKantorovichTest;

IRIT_GLOBAL_DATA_HEADER CagdRType
    _MVGlblZeroDmnExtension;

IRIT_GLOBAL_DATA_HEADER MvarZeroSubdivTolActionType
    _MVGlblUponSubdivTol;

IRIT_GLOBAL_DATA_HEADER MvarMVsZerosSubdivCallBackFuncType
    _MVGlblZeroSubdivCallBackFunc;

#ifdef DEBUG
IRIT_GLOBAL_DATA_HEADER int
    _DebugZeroSolverBySolDim[3],
    _DebugZeroSolverDmnNum,
    _DebugZeroZeroErr;
#endif /* DEBUG */

MvarZeroPrblmStruct *MvarZeroSolverPrblmNew(const MvarMVStruct * const *MVs, 
					    const MvarExprTreeStruct * const
					    *ETs,
					    int NumOfConstraints,
					    MvarConstraintType *Constraints,
					    CagdRType SubdivTol,
					    CagdRType NumericTol,
					    CagdRType StepTol,
					    CagdBType Solve2DBy0D);
MvarZeroPrblmStruct *MvarZeroSolverSubProblem
					 (MvarZeroPrblmStruct const *Problem,
					  MvarMVStruct **MVs,
					  MvarExprTreeEqnsStruct *Eqns,
					  MvarZeroSolutionStruct *BoundarySol);
void MvarZeroSolverPrblmFree(MvarZeroPrblmStruct *Problem);
MvarZeroSolutionStruct *MvarZeroSolverSolutionNew(MvarTriangleStruct *Tr,
						  MvarPolylineStruct *Pl,
						  MvarPtStruct *Pt,
						  Representation Rep);
MvarZeroSolutionStruct *MvarZeroSolverSolCpy(MvarZeroSolutionStruct const *Sol);
MvarZeroSolutionStruct *MvarZeroSolutionCpyList(MvarZeroSolutionStruct 
						const *SolutionList);
MvarPtStruct *MvarZeroGenPtMidDmn(const MvarZeroPrblmStruct *Problem,
				  int SingleSol);
MvarPtStruct *MvarZero0DNumeric(MvarPtStruct *ZeroPt,
				const MvarExprTreeEqnsStruct *Eqns,
				MvarMVStruct const * const *MVs,
				int NumMVs,
				CagdRType NumericTol,
				const CagdRType *InputMinDmn,
				const CagdRType *InputMaxDmn);
void MvarZeroSolverSolutionFree(MvarZeroSolutionStruct *Solution,
				CagdBType FreeUnion);
MvarZeroSolutionStruct *MvarZeroSolverOrganizeSol2DMVs
						(MvarZeroSolutionStruct *Sol,
						 MvarZeroPrblmStruct *Problem);
void MvarZeroHandleTJunctions(MvarPolylineStruct *Loops, 
			      MvarZeroTJunctionStruct *TJList);
MvarTriangleStruct *MvarTriangleNew(int Dim);
void MvarTriangleFree(MvarTriangleStruct *Tr);
void MvarTriangleFreeList(MvarTriangleStruct *TrList);

MvarZeroTJunctionStruct *MvarZeroTJNew(const MvarPtStruct *TJPrev,
				       const MvarPtStruct *TJPt,
				       const MvarPtStruct *TJNext);
void MvarZeroTJFreeList(MvarZeroTJunctionStruct *TJList);
void MvarZeroTJFree(MvarZeroTJunctionStruct *TJ);
MvarZeroTJunctionStruct *MvarZeroTJCopy(const MvarZeroTJunctionStruct *TJ);
MvarZeroTJunctionStruct *MvarZeroTJCopyList(
                                        const MvarZeroTJunctionStruct *TJList);
CagdBType MvarZeroIsPtOnEdge(const MvarPtStruct *Pt,
			     const CagdRType *MinDmn,
			     const CagdRType *MaxDmn,
			     CagdRType Tol);
CagdBType MvarZeroIsPtOnDmnBndry(CagdRType *MinDmn,
			         CagdRType *MaxDmn,
				 const MvarPtStruct *Pt);
CagdBType MvarZeroHasC1Discont(MvarMVStruct * const *MVs,
			       int NumOfMVs,
			       int *JLoc,
			       CagdRType *t);
void MvarZeroUpdateProblemDmnMVs(MvarZeroPrblmStruct *Problem);
void MvarZeroUpdateProblemDmnExpTr(MvarZeroPrblmStruct *Problem);
CagdBType MvarZeroFirstSmoothUpdatesMVs(MvarZeroPrblmStruct *Problem);
CagdBType MvarZeroFirstSmoothUpdatesExpTr(MvarZeroPrblmStruct *Problem);
MvarZeroSubDmnInfoStruct *MvarSubDmnInfoStructNew(MvarMVStruct **MVs, 
						  MvarMVDirType ProjDir1,
						  MvarMVDirType ProjDir2);
CagdBType MvarConesOverlapAux(const MvarNormalConeStruct *ConesList);
void MvarZR1DDbgPrintPtList(const MvarPtStruct * const PtList);
CagdBType MvarZeroSolverIsMVZero(MvarMVStruct *MV,
				 CagdRType NumericTol);
void MvarSubDmnInfoStructFree(MvarZeroSubDmnInfoStruct *InfoStruct,
			      int NumOfMVs);
CagdRType MvarMVEvalErrorL1(MvarMVStruct const * const *MVs,	
			    CagdRType *Params,
			    int NumOfMVs);
int MvarMVsReduceMvsDomains(MvarMVStruct **MVs,
			    int NumOfMVs,
			    int Dir,
			    CagdRType SubdivTol,
			    CagdRType *TMin,
			    CagdRType *TMax,
			    CagdBType SameSpace);
void MvarGetSubdivParamDomains(const MvarMVStruct *MV,		
			       CagdRType *Min,
			       CagdRType *Max,
			       int Dir);

/******************************************************************************
* This macro is called when the library has detected an unrecoverable error.  *
* Default action is to call MvarFatalError, but you may want to reroute this  *
* to invoke your handler and recover yourself (by long jump for example).     *
******************************************************************************/
#define MVAR_FATAL_ERROR(Msg)	MvarFatalError(Msg)

/******************************************************************************
* Voronoi cell computation.						      *
******************************************************************************/
int MvarBsctIsCurveLL(MvarVoronoiCrvStruct *Cv);
MvarVoronoiCrvStruct *MvarBsctApplyLL(MvarVoronoiCrvStruct *Cv);
int MvarBsctApplyCC(MvarVoronoiCrvStruct *Cv1, 
		    MvarVoronoiCrvStruct **CCFreeCrvs);
MvarVoronoiCrvStruct *MvarBsctPurgeAwayLLAndCCConstraints(MvarVoronoiCrvStruct 
							          *InputCrvs);
void MvarBsctTrimCurveBetween(MvarVoronoiCrvStruct *Cv, 
			      MvarPtStruct *Pt1, 
			      MvarPtStruct *Pt2, 
			      MvarVoronoiCrvStruct **TrimmedCurve);
void MvarBsctComputeLowerEnvelope(MvarVoronoiCrvStruct *inputCurves, 
				  MvarVoronoiCrvStruct **lowerEnvelope);
MvarPtStruct *MvarBsctImplicitCrvExtremeAliter(CagdSrfStruct *Srf,
					       CagdSrfDirType Dir,
					       CagdRType MvarBsctSubdivTol,
					       CagdRType MvarBsctNumerTol);
CagdSrfStruct *MvarBsctTrimSurfaceByUVBbox(CagdSrfStruct *Srf, 
					   CagdBBoxStruct UVBbox);
MvarPtStruct *MvarBsctNewFindZeroSetOfSrfAtParam(CagdSrfStruct *Srf, 
						 CagdRType Param,
						 CagdSrfDirType Dir, 
						 CagdRType MvarBsctSubdivTol,
						 CagdRType MvarBsctNumerTol,
						 CagdBType ShouldCheckEndPoints);
void MvarBsctSplitImplicitCrvToMonotonePieces(CagdSrfStruct *Srf,
					      CagdSrfStruct **OutLst,
					      CagdRType MvarBsctSubdivTol,
					      CagdRType MvarBsctNumerTol);
void MvarBsctComputeDenomOfP(CagdCrvStruct *Crv1Inp,	
			     CagdCrvStruct *Crv2Inp,	
			     CagdSrfStruct **DenomOut);
CagdRType *MvarComputeInterMidPoint(CagdCrvStruct *Crv1,
				    CagdRType t1,
				    CagdCrvStruct *Crv2,
				    CagdRType t2);
void MvarBsctComputeF3(CagdCrvStruct *Crv1Inp,                       
		       CagdCrvStruct *Crv2Inp,			           
		       CagdCrvStruct **Crv1Coerced,
		       CagdCrvStruct **Crv2Coerced,	
		       CagdSrfStruct **F3,
		       CagdSrfStruct **L1,			           
		       CagdSrfStruct **L2,			           
		       CagdSrfStruct **CC1,			           
		       CagdSrfStruct **CC2);
CagdRType *MvarBsctComputeXYFromBisTR(CagdCrvStruct *Crv1,
				      CagdRType t,
				      CagdCrvStruct *Crv2,
				      CagdRType r);
MvarPtStruct *MvarBsctSkel2DEqPts3Crvs(CagdCrvStruct *Crv1,
				       CagdCrvStruct *Crv2,
				       CagdCrvStruct *Crv3);
MvarVoronoiCrvStruct *MvarVoronoiCrvNew(void);
MvarVoronoiCrvStruct *MvarVoronoiCrvCopy(MvarVoronoiCrvStruct *Crv);
void MvarVoronoiCrvFree(MvarVoronoiCrvStruct *Crv);
void MvarVoronoiCrvFreeList(MvarVoronoiCrvStruct *CrvList);
MvarVoronoiCrvStruct *MvarVoronoiCrvReverse(MvarVoronoiCrvStruct *Crv);
int MvarBsctIsXSmaller(MvarPtStruct *P1, MvarPtStruct *P2);
void MvarBsctCurveLeft(MvarVoronoiCrvStruct *Cv, MvarPtStruct *Res);
void MvarBsctCurveRight(MvarVoronoiCrvStruct *Cv, MvarPtStruct *Res);
int MvarBsctCv1IsYSmallerAt(MvarVoronoiCrvStruct *Cv1, 
			    MvarVoronoiCrvStruct *Cv2, 
			    MvarPtStruct *MidPoint);
void MvarBsctGetAllIntersectionPoints(MvarVoronoiCrvStruct *Cv1, 
				      MvarVoronoiCrvStruct *Cv2, 
				      MvarPtStruct **Points);
void MvarBsctSplitCurve(MvarVoronoiCrvStruct *Cv, 
			MvarPtStruct *SplitPt, 
			MvarVoronoiCrvStruct **CvLeft,
			MvarVoronoiCrvStruct **CvRight);
CagdCrvStruct *MvarBsctTrimCrvPt(CagdCrvStruct *Crv, 
				 CagdRType *Pt, 
				 CagdRType Alpha,
				 CagdCrvStruct *BaseCrv);

/******************************************************************************
* Minimal spanning circle/sphere computation.				      *
******************************************************************************/
int MvarMSConstraints(CagdPtStruct *PointList,
		      struct IPObjectStruct **MVObjs,
		      int LenOfMVs,
		      CagdRType *Center,
		      CagdRType *Radius,
		      CagdRType SubdivTol,
		      CagdRType NumerTol);
MvarMVStruct **MvarMSConstraintsOfAPair(CagdPtStruct *PointList,
					struct IPObjectStruct **MVObjs,
					int LenOfMVs,
					int Dim);
MvarMVStruct **MvarMSConstraintsOfATriplet(CagdPtStruct *PointList,
					   struct IPObjectStruct **MVObjs,
					   int LenOfMVs,
					   MvarMVStruct **MVPDenom1, 
					   MvarMVStruct **MVPNumer1,
					   int Dim);
/******************************************************************************
* 2 contact motion computation.						      *
******************************************************************************/

CagdRType Mvar2CtLineLineDist(Mvar2CtLineStruct *A,
			      Mvar2CtLineStruct *B);
CagdRType Mvar2CtLinePointDist(CagdPType P, CagdPType L[2]);
Mvar2CtBVHStruct *Mvar2CtBuildBVH(CagdCrvStruct *Crv, 
				  CagdRType SubdivTol, 
				  CagdRType BvTol);
void Mvar2CtSetNodeId(Mvar2CtBVNodeStruct *Node, int Id);
void Mvar2CtFreeBVH(Mvar2CtBVHStruct *Bvh);
CagdRType Mvar2CtPenetrationDepth(Mvar2CtBVHStruct *BvhA, 
			          Mvar2CtBVHStruct **BvhBs,
			          int BSize,
				  CagdRType Xtrans,
				  CagdRType Ytrans,
				  CagdRType Rot);
CagdRType Mvar2CtGetTheta(CagdRType x, CagdRType y);
CagdBType Mvar2CtNormalOverlap(Mvar2CtBVNodeStruct *ANode, 
			       Mvar2CtBVNodeStruct *BNode, 
			       CagdRType RMin, 
			       CagdRType RMax);
CagdBType Mvar2CtNormalOverlapBoth(Mvar2CtBVNodeStruct *ANode, 
			           Mvar2CtBVNodeStruct *BNode);
void Mvar2CtBuildCParamHierarchy(CagdCrvStruct *Circle, 
				 Mvar2CtCParamStruct **Node, 
				 CagdRType Min, 
				 CagdRType Max, 
				 CagdRType Tol);
void Mvar2CtGetParentCparam(Mvar2CtCParamStruct *Node,
		            CagdRType Min,
		            CagdRType Max,
		            Mvar2CtCParamStruct **Parent);
void Mvar2CtGetParentBVNode(Mvar2CtBVNodeStruct *Node,
		            CagdRType Min,
			    CagdRType Max,
		            Mvar2CtBVNodeStruct **Parent);
void Mvar2CtFreeCparam(Mvar2CtCParamStruct *Node);
CagdBType Mvar2CtCheck2CtTrace(Mvar2CtBVNodeStruct *Nodes[4], 
			       Mvar2CtCParamStruct *Cparam,  
			       CagdRType Tol);
void Mvar2CtReduce2CtDomain(Mvar2CtBVNodeStruct *Nodes[4], 
			    Mvar2CtCParamStruct *Cparam, 
			    CagdRType Min[5], 
			    CagdRType Max[5],
			    int MinMax,
			    int FixedDir,
			    CagdRType Tol);
void Mvar2CtReduce3CtDomain(Mvar2CtBVNodeStruct *Nodes[6], 
			    Mvar2CtCParamStruct *Cparam, 
			    CagdRType Min[7], 
			    CagdRType Max[7]);
void Mvar2CtReduceRotExtremeDomain(Mvar2CtBVNodeStruct *Nodes[4], 
			           Mvar2CtCParamStruct *Cparam, 
				   CagdRType Min[5], 
			           CagdRType Max[5], 
				   CagdRType Tol);
MvarPolylineStruct * Mvar2CtConnectPeriodic(MvarPolylineStruct *Polys,
				            CagdRType Tol);
MvarMVStruct ** Mvar2CtExtractMVRegion(MvarMVStruct **MVs, 
				       int MV_Num, 
				       CagdRType *Min, 
				       CagdRType *Max);
CagdBType Mvar2CtIsConnectedNode(Mvar2CtBVNodeStruct *Node1,
				 Mvar2CtBVNodeStruct *Node2);
MvarPtStruct* Mvar2CtGetMiddlePt(MvarPtStruct *PtList, int Length);
int Mvar2CtTraceCollide(MvarPolylineStruct *Poly1,
			MvarPolylineStruct *Poly2);
void Mvar2CtSwapTrace(MvarPolylineStruct *MPoly);
MvarPolylineStruct * Mvar2CtValidateTraces(MvarPolylineStruct *Polys,
				           Mvar2CtBVHStruct *BvhA,
				           Mvar2CtBVHStruct **BvhBs,
				           int BSize,
				           CagdCrvStruct *Circle);
MvarPtStruct * Mvar2CtValidate2Ct(MvarPtStruct *MPts,
			          Mvar2CtBVHStruct *BvhA,
			          Mvar2CtBVHStruct **BvhBs,
			          int BSize,
			          CagdCrvStruct *Circle);
MvarPtStruct * Mvar2CtValidateCurvContact(MvarPtStruct *MPts,
				          Mvar2CtBVHStruct *BvhA,
				          Mvar2CtBVHStruct **BvhBs,
				          int BSize,
				          int BIndex,
					  CagdCrvStruct *Circle);
CagdBType Mvar2CtInDomain(CagdRType *Min,
		          CagdRType *Max,
		          MvarPtStruct *MPt);
CagdBType Mvar2CtIsPassing(CagdRType *Min,
		           CagdRType *Max,
		           MvarPolylineStruct *MPoly);
CagdBType Mvar2CtRejectbyCurvature(Mvar2CtBVNodeStruct *Node1, 
				   Mvar2CtBVNodeStruct *Node2);
CagdBType Mvar2CtCurvatureOverlap(Mvar2CtBVNodeStruct *ANode, 
			          Mvar2CtBVNodeStruct *BNode,
			          Mvar2CtCParamStruct *Cparam);
void Mvar2CtTraceBBox(CagdRType *Min, 
		      CagdRType *Max, 
		      MvarPtStruct *SPt, 
	              MvarPtStruct *EPt);

/******************************************************************************
* Zero set computation.							      *
******************************************************************************/
MvarPtStruct *MvarZeroFilterSolutionSet(MvarPtStruct *ZeroSet,
					const MvarMVStruct * const *MVs,
					const MvarConstraintType *Constraints,
					int NumOfMVs,
					CagdRType Tol,
					int CanHaveLoops,
					int SortSol, 
					CagdBType InEqOnly);
CagdBType MvarZeroMVConstraintFail(const MvarMVStruct *MV,
				   MvarConstraintType Constraint);
MvarPtStruct *MvarZeroGetRootsByKantorovich(MvarMVStruct **MVs,
                                            MvarConstraintType *Constraints,
					    int NumOfMVs,
					    int NumOfZeroMVs,
					    int ApplyNormalConeTest,
					    CagdRType SubdivTol,
					    int Depth,
					    CagdBType SameSpace,
					    CagdRType ParamPerturb);
void MvarZeroGetSubdivParamDomains(const MvarMVStruct *MV,
				   CagdRType *Min,
				   CagdRType *Max,
				   int Dir);

int MvarMVZR1DNoLoopTest(MvarMVStruct const * const *MVs, 
		         MvarMVZR1DAuxStruct *AS);
MvarMVZR1DAuxStruct *MvarMVZR1DAllocOnce(int Number);
int MvarZeroSolverGetDmnDim(MvarZeroPrblmStruct const *Problem);


MvarMVStruct **MvarZeroOrganizeMVs0DProblem(const MvarMVStruct * const *MVs,
					    MvarConstraintType *Constraints,
					    int *NumOfMVs);
MvarMVStruct **MvarZeroOrganizeMVs1DProblem(const MvarMVStruct * const *MVs,
					    MvarConstraintType *Constraints,
					    int *NumOfMVs);
MvarExprTreeEqnsStruct *MvarZeroOrganizeETs0DProblem(
				      const MvarExprTreeStruct * const *MVETs,
				      int NumOfMVETs);

void MvarZeroSolverSetCallbackFcns0DMVs(MvarZeroPrblmStruct *Problem);
void MvarZeroSolverSetCallbackFcns0DExpTr(MvarZeroPrblmStruct *Problem);
void MvarZeroSolverSetCallbackFcns1DMVs(MvarZeroPrblmStruct *Problem);
void MvarZeroSolverSetCallbackFcns1DExpTr(MvarZeroPrblmStruct *Problem);
void MvarZeroSolverSetCallbackFcns2DMVs(MvarZeroPrblmStruct *Problem);
void MvarZeroSolverSetCallbackFcns2DExpTr(MvarZeroPrblmStruct *Problem);

int MvarZeroSolverGetDmnDim(const MvarZeroPrblmStruct *Problem);

int MVarMVHyperPlanesTestForSol(MvarMVStruct const * const *MVs,
				int NumOfZeroMVs);
void MvarMVZR1DDeallocOnce(MvarMVZR1DAuxStruct *AS);
MvarPtStruct *MvarMVZR1DListOfDifferentPts(MvarPtStruct *PtList, 
     				           CagdRType Tol);
MvarPolylineStruct *MvarMVZR1DPolyWithDom(MvarPolylineStruct *Poly, 
				          const MvarMVStruct *MV);
MvarPtStruct *MvarMVZR1DCurveTracing(MvarMVStruct * const *MVs,  
				     const MvarPtStruct *StartPoint,
				     const MvarPtStruct *EndPoint,
				     const MvarVecStruct *DirVec,
				     CagdRType Step, 
				     CagdRType NumericTol,
				     MvarMVGradientStruct **MVGrads,
				     MvarMVZR1DAuxStruct *AS,
				     int *TraceError);
MvarVecStruct *MvarMVZR1DStartVec(MvarMVStruct * const *MVs, 
				  MvarPtStruct *BoundaryPt,
				  MvarMVZR1DAuxStruct *AS,
				  CagdRType NumericTol);
MvarPolylineStruct *MvarMVZR1DLinkNeighbours(MvarPolylineStruct *PolyList1, 
  					     MvarPolylineStruct *PolyList2, 
					     CagdRType SubdivTol, 
					     CagdRType NumericTol,
					     int BoundarySide, 
					     CagdRType BoundaryValue);
void MvarMVZR1DSplitBoundaryPts(const MvarPtStruct *BoundaryPts,
				int BoundarySide, 
				CagdRType BoundaryValue,
				MvarPtStruct **SplitPts0,
				MvarPtStruct **SplitPts1);
MvarPtStruct *MvarMVZR1DMiddlePlaneCutPts(MvarMVStruct * const *MVs,
					  MvarMVZR1DAuxStruct *AS,
					  CagdRType SubdivTol,
					  CagdRType NumericTol,
					  int BoundarySide, 
					  CagdRType BoundaryValue);

/******************************************************************************
* Zero set ET computation.						      *
******************************************************************************/
/* Expression trees' solving aux. routines. */
int MvarEqnsPerturbZero(MvarPtStruct *Pt,
			CagdRType StartRngErr,
			const MvarExprTreeEqnsStruct *Eqns,
			MvarMVStruct const * const *MVs,
			int NumEqns,
			CagdRType *MinDmn,
			CagdRType *MaxDmn,
			CagdRType Step);
CagdRType MvarExprTreeEqnsEvalErrorL1(const MvarExprTreeEqnsStruct *Eqns,
				      CagdRType *Params);

/* Expression trees' equations handling/building routine. */
MvarExprTreeEqnsStruct *MvarExprTreeEqnsMalloc(int NumEqns,
					       int MaxNumCommonExprs);
void MvarExprTreeEqnsFree(MvarExprTreeEqnsStruct *Eqns);
void MvarExprTreeEqnsReallocCommonExprs(MvarExprTreeEqnsStruct *Eqns,
					int NewSize);
MvarExprTreeEqnsStruct *MvarExprTreeEqnsBuild(MvarExprTreeStruct * const *MVETs,
					      const MvarConstraintType
					                   *ConstraintsTypes,
					      int NumOfMVETs);

#endif /* MVAR_LOC_H */
