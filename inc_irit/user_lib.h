/******************************************************************************
* User_lib.h - Header file for the User Interaction library.		      *
* This header is also the interface header to the world.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Mar. 90.					      *
******************************************************************************/

#ifndef USER_LIB_H
#define USER_LIB_H

#ifdef __WINNT__
#include <wchar.h>
#define USER_FONT_STR_CONST(Str)	L##Str
#define USER_FONT_STR_CPY(DStr, SStr)	wcscpy(DStr, SStr)
#define USER_FONT_STR_CAT(DStr, SStr)	wcscat(DStr, SStr)
#define USER_FONT_STR_DUP(Str)		_wcsdup(Str)
#define USER_FONT_STR_CHR(Str, Chr)	wcschr(Str, L##Chr)
#define USER_FONT_STR_LEN(Str)		wcslen(Str)
#define USER_FONT_IS_SPACE(c)		iswspace(c)
#define USER_FONT_TEXT2INT(Str)		_wtoi(Str)
#define USER_FONT_GET_WORD_ASCII(Str)	UserWChar2Ascii(Str)
#define USER_FONT_GET_WORD_UNICODE(Str) UserAscii2WChar(Str)
#else
#define USER_FONT_STR_CONST(Str)	Str
#define USER_FONT_STR_CPY(DStr, SStr)	strcpy(DStr, SStr)
#define USER_FONT_STR_CAT(DStr, SStr)	strcat(DStr, SStr)
#define USER_FONT_STR_DUP(Str)		strdup(Str)
#define USER_FONT_STR_CHR(Str, Chr)	strchr(Str, Chr)
#define USER_FONT_STR_LEN(Str)		strlen(Str)
#define USER_FONT_IS_SPACE(c)		isspace(c)
#define USER_FONT_TEXT2INT(Str)		atoi(Str)
#define USER_FONT_GET_WORD_ASCII(Str)	(Str)
#define USER_FONT_GET_WORD_UNICODE(Str) (Str)
#endif /* __WINNT__ */

#include "cagd_lib.h"
#include "geom_lib.h"
#include "iritprsr.h"

#define USER_HC_VEC_DRAW_SCALE	0.25

typedef enum {
    USER_ERR_WRONG_SRF,
    USER_ERR_MISSING_ATTRIB,
    USER_ERR_WRONG_ANGLE,
    USER_ERR_INVALID_PLANE,
    USER_ERR_RATIONAL_NO_SUPPORT,
    USER_ERR_NON_CRV_OBJ_IN_FONT,
    USER_ERR_NO_ADJ_INFO,
    USER_ERR_NO_NRML_INFO,
    USER_ERR_NO_CRVTR_INFO,
    USER_ERR_EXPCT_REG_TRIANG,
    USER_ERR_EXPCT_POLY_OBJ,
    USER_ERR_EXPCT_SRF_OBJ,
    USER_ERR_EXPCT_VRTX_NRMLS,
    USER_ERR_EXPCT_VRTX_UVS,
    USER_ERR_UNDEFINE_ERR,
    USER_ERR_WRONG_CTLPT_INDEX,
    USER_ERR_INVALID_SIZE,
    USER_ERR_INVALID_CURVE,
    USER_ERR_INVALID_SURFACE,
    USER_ERR_INVALID_TRIM_SRF,
    USER_ERR_INVALID_DIR,
    USER_ERR_INVALID_IMAGE_SIZE,
    USER_ERR_INVALID_KV_END_COND,
    USER_ERR_XRANGE_EMPTY,
    USER_ERR_YRANGE_EMPTY,
    USER_ERR_ZRANGE_EMPTY,

    USER_ERR_NC_INVALID_PARAM,
    USER_ERR_NC_INVALID_INTER,
    USER_ERR_NC_NO_POLYLINES,
    USER_ERR_NC_MIX_CRVS_PLLNS
} UserFatalErrorType;

typedef enum {					/* Type of surface marching. */
    USER_SRF_MARCH_ISO_PARAM,
    USER_SRF_MARCH_ISOCLINES,
    USER_SRF_MARCH_ORTHOCLINES,
    USER_SRF_MARCH_PRIN_CRVTR
} UserSrfMarchType;

typedef enum {
    USER_3D_SPREAD_RANDOM,
    USER_3D_SPREAD_DIAG_PLANE,
    USER_3D_SPREAD_DIAG_PLANE2,
    USER_3D_SPREAD_ANTI_DIAG_PLANE,
    USER_3D_SPREAD_ANTI_DIAG_PLANE2,
    USER_3D_SPREAD_ANTI_DIAG_PLANE3,
    USER_3D_SPREAD_DISCONT2PLANE,
    USER_3D_SPREAD_DISCONT4PLANE,
} User3DSpreadType;

typedef enum {
    USER_IMG_SHD_3D_BLOB_NO_COLOR,
    USER_IMG_SHD_3D_BLOB_GRAY_LEVEL,
    USER_IMG_SHD_3D_BLOB_COLOR,
} UserImgShd3dBlobColorType;

typedef enum {
    USER_CA_SPLIT_NONE = 0x0000,
    USER_CA_SPLIT_INFLECTION_PTS = 0x0001,
    USER_CA_SPLIT_MAX_CRVTR_PTS =  0x0002,
    USER_CA_SPLIT_C1DISCONT_PTS =  0x0004,
    USER_CA_SPLIT_REAL_C1DISCONT_PTS = 0x0008
} UserCASplitType;

typedef enum {
    USER_CA_INPUT_NONE = 0x0000,
    USER_CA_INPUT_POLYLINES =     0x0001,
    USER_CA_INPUT_CURVES =        0x0002,
    USER_CA_INPUT_TCRVS_IN_SRFS = 0x0004
} UserCAInputType;

typedef enum {
    USER_CA_UNDEF_TYPE,
    USER_CA_HANGING_TYPE,
    USER_CA_SIMPLE_TYPE,
    USER_CA_LOOP_TYPE,
    USER_CA_COMPLEX_TYPE,
    USER_CA_LEFTOVERS_TYPE
} UserCAObjType;

typedef enum {
    USER_CA_OPER_NONE,
    USER_CA_OPER_CREATE,
    USER_CA_OPER_COPY,
    USER_CA_OPER_FILTER_DUP,
    USER_CA_OPER_FILTER_TAN,
    USER_CA_OPER_SPLIT_CRV,
    USER_CA_OPER_BREAK_LIN,
    USER_CA_OPER_BREAK_INTER,
    USER_CA_OPER_BREAK_NEAR_PTS,
    USER_CA_OPER_UNION_CRV,
    USER_CA_OPER_LSTSQR_CRV,
    USER_CA_OPER_EVAL_CA,
    USER_CA_OPER_CLASSIFY,
    USER_CA_OPER_REPORT,
    USER_CA_OPER_OUTPUT,
    USER_CA_OPER_FREE,
} UserCAOpType;

/* Type of kinematic point constraints. */
typedef enum {
    USER_KNMTCS_PT_NONE = 0,
    USER_KNMTCS_PT_FIXED,
    USER_KNMTCS_PT_XY_PLANE,
    USER_KNMTCS_PT_X_DIRECTION,
    USER_KNMTCS_PT_Y_DIRECTION,
    USER_KNMTCS_PT_Z_DIRECTION,
    USER_KNMTCS_PT_MOVES_ALONG_CURVE,
    USER_KNMTCS_PT_MOVES_ALONG_SURFACE,
    USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE,
} UserKnmtcsMovabilityPointType;

/* Type of kinematic constraints. */
typedef enum {
    USER_KNMTCS_CONSTR_NONE = 0,
    USER_KNMTCS_CONSTR_DIST_PT_PT,
    USER_KNMTCS_CONSTR_DIST_PT_BAR,
    USER_KNMTCS_CONSTR_DIST_PT_CRV,
    USER_KNMTCS_CONSTR_DIST_PT_SRF,
    USER_KNMTCS_CONSTR_DIST_BAR_BAR,
    USER_KNMTCS_CONSTR_DIST_BAR_CRV,
    USER_KNMTCS_CONSTR_DIST_BAR_SRF,
    USER_KNMTCS_CONSTR_ANGLE,
    USER_KNMTCS_CONSTR_ORTHOGONALITY,
    USER_KNMTCS_CONSTR_TANGENCY,
    USER_KNMTCS_CONSTR_PARALLEL,
    USER_KNMTCS_CONSTR_ROT_CRV,

    USER_KNMTCS_CONSTR_MIN_DIST_PT_PT,/* In equality constraints from here. */
    USER_KNMTCS_CONSTR_MAX_DIST_PT_PT,
    USER_KNMTCS_CONSTR_XDIFF_POS,
    USER_KNMTCS_CONSTR_YDIFF_POS,
    USER_KNMTCS_CONSTR_ZDIFF_POS
} UserKnmtcsConstrType;

typedef struct UserFEKElementStruct {
    IrtRType k[2][2];			  /* (x, y) x (x, y) contributions. */
} UserFEKElementStruct;

typedef struct UserFECElementStruct {
    IrtRType c[2];			           /* (x, y) contributions. */
} UserFECElementStruct;

typedef struct UserFEInterIntervalStruct {
    struct UserFEInterIntervalStruct *Pnext;
    CagdRType T1Min, T1Max;		    /* Interval of overlap in Crv1. */
    CagdRType T2Min, T2Max;		    /* Interval of overlap in Crv2. */
    CagdRType Antipodal1, Antipodal2;  /* Locations of maximal penetration. */
    CagdVType ProjNrml;		    /* Direction to project penetration on. */
} UserFEInterIntervalStruct;

/* Curves arrangement holds a vector of curves, a vector of curves' end      */
/* points, and a vector of regions.					     */
/*   Tolerances and other aux. data are included as well.	             */
typedef struct UserCrvArngmntStruct {
    struct CagdCrvStruct *CagdCrvs;
    struct UserCAPtStruct *Pts;
    struct UserCACrvStruct *Crvs;
    struct UserCARegionStruct **Regions;
    IPObjectStruct *Output;                       /* CA output is kept here. */

    IrtRType EndPtEndPtTol;  /* Tolerance to consider crvs' end points same. */
    IrtRType InternalTol;    /* Internal tolerance for CCI, inflections etc. */
    IrtRType PlanarityTol;   /* Tolerance of curves to be considered planar. */
    IrtHmgnMatType XYZ2XYMat, XY2XYZMat;     /* General plane <--> XY plane. */
    IrtPlnType CrvsPlane;
    int ProjectOnPlane;/* TRUE to force crvs to be planar, on computed plane.*/
    int AllocSizePts;      /* Size of the allocated vectors of Pts and Crvs. */
    int AllocSizeCrvs;
    int NumOfPts;
    int NumOfOrigCrvs;                     /* Number of curves in the input. */
    int NumOfCrvs;			        /* Current number of curves. */
    int NumOfRegions;			       /* Current number of regions. */
    const char *Error; /* Last error string description will be placed here. */
} UserCrvArngmntStruct;

/* Structure which represent a kinematic point. */
typedef struct UserKnmtcsPtStruct {
    struct UserKnmtcsPtStruct *Pnext;
    UserKnmtcsMovabilityPointType Type;
    int Idx;
    CagdPType Pt;
    union {
        CagdCrvStruct *Crv;			    /* Pt moves along curve. */
        CagdSrfStruct *Srf;		 	  /* Pt moves along surface. */
    } U;
    CagdPtStruct Center;			          /* If Crv rotates. */
} UserKnmtcsPtStruct;

/* Structure which represent a kinematic bar. */
typedef struct UserKnmtcsBarStruct {
    struct UserKnmtcsBarStruct *Pnext;
    UserKnmtcsPtStruct *P1, *P2;      /* Starting & ending point of the bar. */
} UserKnmtcsBarStruct;

/* Structure which represent a kinematic constraint. */
typedef struct UserKnmtcsConstrStruct {
    struct UserKnmtcsConstrStruct *Pnext;
    UserKnmtcsConstrType Type;
    union{
        CagdRType distance;
        CagdRType angle;
    } V;
    union{
        struct {
            UserKnmtcsPtStruct *Pt1;			  /* Distance PT_PT. */
            UserKnmtcsPtStruct *Pt2;
        } DstPtPt;
        struct {
            UserKnmtcsPtStruct *Pt;			 /* Distance PT_BAR. */
            UserKnmtcsBarStruct *Bar;
        } DstPtBar;
        struct {
            UserKnmtcsPtStruct *Pt;			 /* Distance PT_CRV. */
	    UserKnmtcsPtStruct *CrvPt;			      /* Foot point. */
        } DstPtCrv;
	struct {
            UserKnmtcsPtStruct *Pt;			 /* Distance PT_SRF. */
	    UserKnmtcsPtStruct *SrfPt;			      /* Foot point. */
        } DstPtSrf;
        struct {
            UserKnmtcsBarStruct *Bar1;			/* Distance BAR_BAR. */
            UserKnmtcsBarStruct *Bar2;
        } DstBarBar;        
        struct {
            UserKnmtcsBarStruct *Bar;			/* Distance BAR_CRV. */
            CagdCrvStruct *Crv;
        } DstBarCrv;
        struct {
            UserKnmtcsBarStruct *Bar1;		    /* Angle, orthogonality. */
            UserKnmtcsBarStruct *Bar2;
        } Angle;
        struct {
            UserKnmtcsBarStruct *Bar;				/* Tangnecy. */
	    UserKnmtcsPtStruct *Pt;			   /* Contact point. */
        } Tan;
	struct {
            UserKnmtcsBarStruct *Bar;				/* Tangnecy. */
            CagdSrfStruct *Srf;
	    UserKnmtcsPtStruct *Pt;			   /* Contact point. */
        } TanSrf;
        struct {
            UserKnmtcsBarStruct *Bar1;			        /* Parallel. */
            UserKnmtcsBarStruct *Bar2;
        } Par;
    } C;
} UserKnmtcsConstrStruct;

/* Structure which represent all the kinematic simulations data. */
typedef struct UserKnmtcsStruct {
    CagdRType XMin;
    CagdRType XMax;
    CagdRType YMin;
    CagdRType YMax;
    CagdRType ZMin;
    CagdRType ZMax;
    int PtsNum;
    int BarsNum;
    int ConstraintsNum;
    struct UserKnmtcsPtStruct *Pts;		   /* Pointer to point list. */
    struct UserKnmtcsBarStruct *Bars;		     /* Pointer to bar list. */
    struct UserKnmtcsConstrStruct *Constraints;	     /* List of constraints. */
} UserKnmtcsStruct;

/* Font styles. */ 
typedef enum {
    USER_FONT_STYLE_REGULAR,
    USER_FONT_STYLE_ITALICS,
    USER_FONT_STYLE_BOLD,
    USER_FONT_STYLE_BOLD_ITALICS
} UserFontStyleType;

typedef enum {
    USER_FONT_3D_EDGE_NORMAL,
    USER_FONT_3D_EDGE_CHAMFER,
    USER_FONT_3D_EDGE_ROUND,
} UserFont3DEdgeType;

typedef enum {
    USER_FONT_ALIGN_LEFT,
    USER_FONT_ALIGN_CENTER,
    USER_FONT_ALIGN_RIGHT,
    USER_FONT_ALIGN_WIDE
} UserFontAlignmentType;

typedef enum {
    USER_FONT_OUTPUT_BEZIER_CRVS = 0,
    USER_FONT_OUTPUT_BSPLINE_CRVS,
    USER_FONT_OUTPUT_FILLED2D_POLYS,
    USER_FONT_OUTPUT_OUTLINE_FILLED2D_POLYS,
    USER_FONT_OUTPUT_SOLID3D_POLYS,
    USER_FONT_OUTPUT_FILLED2D_TSRFS,
    USER_FONT_OUTPUT_SOLID3D_TSRFS
} UserFontGeomOutputType;

typedef struct UserFontDimInfoStruct {
    IrtRType DescentLineHeight;     /* The four height lines of this font. */
    IrtRType BaseLineHeight;
    IrtRType MeanLineHeight;
    IrtRType AscentLineHeight;
    IrtRType SpaceWidth;              /* The estimated space width to use. */
    GMBBBboxStruct BBox;   /* Of a generic char 's' in this font/size etc. */
} UserFontDimInfoStruct;

#ifdef __WINNT__
typedef wchar_t *UserFontText;
typedef wchar_t UserFontChar;
#else
typedef char *UserFontText;
typedef char UserFontChar;
#endif /* __WINNT__ */

typedef char *UserFontName;

typedef struct UserFontWordLayoutStruct {
    struct UserFontWordLayoutStruct *Pnext;
    UserFontText Word;
    UserFontName FontName;
    UserFontStyleType FontStyle;
    IrtRType RelSize;                         /* Relative scale to the text. */
    UserFont3DEdgeType Font3DEdge;
    IrtPtType Font3DOptions;
    UserFontAlignmentType FontAlignment;
    IPObjectStruct *Geom;             /* The geometry representing the text. */
    GMBBBboxStruct BBox;       /* BBox of Geom, ignoring (X, Y) translation. */
    IrtRType X, Y;					   /* Word position. */
    IrtRType LeftOverSpace;  /* For last word in line only.  Otherwise zero. */
    IrtBType NewLine;			      /* A new line after this word. */
} UserFontWordLayoutStruct;

typedef void (*UserSetErrorFuncType)(UserFatalErrorType);
typedef int (*UserRegisterTestConverganceFuncType)(IrtRType CrntDist, int i);
typedef int (*UserCntrIsValidCntrPtFuncType)(const CagdSrfStruct *Srf,
					     CagdRType U,
					     CagdRType V);
typedef void (*UserHCEditDrawCtlPtFuncType)(int PtIndex,
					    int PtUniqueID,
					    IrtRType *Pos,
					    IrtRType *TanBack,
					    IrtRType *TanForward);

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/* Surface-Primitive Geometry (rays, points, etc.) interactions. */

VoidPtr IntrSrfHierarchyPreprocessSrf(const CagdSrfStruct *Srf,
				      IrtRType FineNess);
void IntrSrfHierarchyFreePreprocess(VoidPtr Handle);
CagdBType IntrSrfHierarchyTestRay(VoidPtr Handle,
				  CagdPType RayOrigin,
				  CagdVType RayDir,
				  CagdUVType InterUV);
CagdBType IntrSrfHierarchyTestPt(VoidPtr Handle,
				 CagdPType Pt,
				 CagdBType Nearest,
				 CagdUVType InterUV);

/* Surface-plane contouring. */

IPPolygonStruct *UserCntrSrfWithPlane(const CagdSrfStruct *Srf,
                                      const IrtPlnType Plane,
                                      IrtRType FineNess,
				      int UseSSI);
IPPolygonStruct *UserCntrEvalToE3(const CagdSrfStruct *Srf,
				  IPPolygonStruct *Cntrs,
				  UserCntrIsValidCntrPtFuncType
				                             ValidCntrPtFunc);

/* Linear Bsplines vs polylines convertion. */

CagdCrvStruct *UserPolyline2LinBsplineCrv(const IPPolygonStruct *Poly,
					  CagdBType FilterDups);
CagdCrvStruct *UserPolylines2LinBsplineCrvs(const IPPolygonStruct *Polys,
					    CagdBType FilterDups);
IPPolygonStruct *UserCagdPolyline2IritPolyline(const CagdPolylineStruct *Poly);
IPPolygonStruct *UserCagdPolylines2IritPolylines(const CagdPolylineStruct
						                      *Polys);

IPPolygonStruct *UserCnvrtLinBspCrv2IritPolyline(const CagdCrvStruct *Crv);
IPPolygonStruct *UserCnvrtLinBspCrvs2IritPolylines(const CagdCrvStruct *Crvs);

/* Surface cone decomposition. */

IPObjectStruct *UserSrfVisibConeDecomp(const CagdSrfStruct *Srf,
				       CagdRType Resolution,
				       CagdRType ConeSize);
TrimSrfStruct *UserVisibilityClassify(const IPObjectStruct *SclrSrf,
				      TrimSrfStruct *TrimmedSrfs);
IPObjectStruct *UserViewingConeSrfDomains(const CagdSrfStruct *Srf,
					  const CagdSrfStruct *NSrf,
					  const IPPolygonStruct *ConeDirs,
					  CagdRType SubdivTol,
					  CagdRType ConeSize,
					  CagdRType Euclidean);
IPPolygonStruct *UserSrfTopoAspectGraph(CagdSrfStruct *PSrf,
					CagdRType SubdivTol);

/* Surface marching. */

IPPolygonStruct *UserMarchOnSurface(UserSrfMarchType MarchType,
				    const CagdUVType UVOrig,
				    const CagdVType DirOrig,
				    const CagdSrfStruct *Srf,
				    const CagdSrfStruct *NSrf,
				    const CagdSrfStruct *DuSrf,
				    const CagdSrfStruct *DvSrf,
				    CagdRType Length,
				    CagdRType FineNess,
				    CagdBType ClosedInU,
				    CagdBType ClosedInV);
IPPolygonStruct *UserMarchOnPolygons(const IPObjectStruct *PObj,
				     UserSrfMarchType MarchType,
				     const IPPolygonStruct *PlHead,
				     IPVertexStruct *VHead,
				     CagdRType Length);

/* Curve/Surface visibility and accessibility. */

IPObjectStruct *UserCrvViewMap(const CagdCrvStruct *Crv,
			       const CagdCrvStruct *ViewCrv,
			       CagdRType SubTol,
			       CagdRType NumTol,
			       CagdBType TrimInvisible);
IPObjectStruct *UserCrvAngleMap(const CagdCrvStruct *Crv,
				CagdRType SubdivTol,
				CagdRType Angle);
IPObjectStruct *UserCrvOMDiagExtreme(const CagdCrvStruct *Crv,
				     const IPObjectStruct *OM,
				     int DiagExtRes);

CagdCrvStruct *UserCrvVisibleRegions(const CagdCrvStruct *Crv,
				     const CagdRType *View,
				     CagdRType Tolerance);

TrimSrfStruct *UserMoldReliefAngle2Srf(const CagdSrfStruct *Srf,
				       const CagdVType VDir,
				       CagdRType Theta,
				       int MoreThanTheta,
				       CagdRType SubdivTol);
CagdSrfStruct *UserMoldRuledRelief2Srf(const CagdSrfStruct *Srf,
				       const CagdVType VDir,
				       CagdRType Theta,
				       CagdRType SubdivTol);

/* Minimal distance to polylines/gons. */

IrtRType UserMinDistLineBBox(const IrtPtType LinePos,
			     const IrtVecType LineDir,
			     IrtBboxType BBox);
IrtRType UserMinDistLinePolygonList(const IrtPtType LinePos,
				    const IrtVecType LineDir,
				    IPPolygonStruct *Pls,
				    IPPolygonStruct **MinPl,
				    IrtPtType MinPt,
				    IrtRType *HitDepth,
				    IrtRType *IndexFrac);
IrtRType UserMinDistLinePolylineList(const IrtPtType LinePos,
				     const IrtVecType LineDir,
				     IPPolygonStruct *Pls,
				     int PolyClosed,
				     IPPolygonStruct **MinPl,
				     IrtPtType MinPt,
				     IrtRType *HitDepth,
				     IrtRType *IndexFrac);
IrtRType UserMinDistPointPolylineList(const IrtPtType Pt,
				      IPPolygonStruct *Pls,
				      IPPolygonStruct **MinPl,
				      IPVertexStruct **MinV,
				      int *Index);

/* Surface surface intersection. */

int UserSrfSrfInter(const CagdSrfStruct *Srf1,
		    const CagdSrfStruct *Srf2,
		    int Euclidean,
		    CagdRType Eps,
		    int AlignSrfs,
		    CagdCrvStruct **Crvs1,
		    CagdCrvStruct **Crvs2);

/* Jacobian of trivariates and zero set. */

IPObjectStruct *UserTVZeroJacobian(const TrivTVStruct *Tv,
				   CagdBType Euclidean,
				   int SkipRate,
				   const CagdRType Fineness[3]);
IPObjectStruct *UserTrivarZeros(const TrivTVStruct *Tv,
				const TrivTVStruct *TvEuclidean,
				int SkipRate,
				const CagdRType Fineness[3]);

/* Z direction collision. */

IrtRType UserTwoObjMaxZRelMotion(IPObjectStruct *PObj1,
				 IPObjectStruct *PObj2,
				 IrtRType FineNess,
				 int NumIters);

/* Create 3D geometric statues from a set of images. */

IPObjectStruct *UserMake3DStatueFrom2Images(const char *Image1Name,
					    const char *Image2Name,
					    int DoTexture,
					    const IPObjectStruct *Blob,
					    User3DSpreadType BlobSpreadMethod,
					    UserImgShd3dBlobColorType
						             BlobColorMethod,
					    int Resolution,
					    int Negative,
					    IrtRType Intensity,
					    IrtRType MinIntensity,
					    int MergePolys);
IPObjectStruct *UserMake3DStatueFrom3Images(const char *Image1Name,
					    const char *Image2Name,
					    const char *Image3Name,
					    int DoTexture,
					    const IPObjectStruct *Blob,
					    User3DSpreadType BlobSpreadMethod,
					    UserImgShd3dBlobColorType
						             BlobColorMethod,
					    int Resolution,
					    int Negative,
					    IrtRType Intensity,
					    IrtRType MinIntensity,
					    int MergePolys);
IPObjectStruct *User3DMicroBlobsFrom3Images(const char *Image1Name,
					    const char *Image2Name,
					    const char *Image3Name,
					    User3DSpreadType BlobSpreadMethod,
					    IrtRType Intensity,
					    const IrtVecType MicroBlobSpacing,
					    const IrtVecType RandomFactors,
					    int Resolution,
					    int Negative,
					    IrtRType CubeSize,
					    int MergePts);
IPPolygonStruct *User3DMicroBlobsTiling(IrtRType XZIntensity,
					IrtRType YZIntensity,
					IrtRType XYIntensity,
					const IrtVecType MicroBlobSpacing,
					const IrtVecType RandomFactors);
IPPolygonStruct *User3DMicroBlobsTiling2(IrtRType XZIntensity,
					 IrtRType YZIntensity,
					 IrtRType XYIntensity,
					 const IrtVecType MicroBlobSpacing,
					 const IrtVecType RandomFactors);
int *User3DMicroBlobsCreateRandomVector(int Size,
					User3DSpreadType BlobSpreadMethod,
					IrtBType FirstVec);
int **User3DMicroBlobsCreateRandomMatrix(int Size,
					 User3DSpreadType BlobSpreadMethod);

IPVertexStruct *User3DDitherSetXYTranslations(IPVertexStruct *Vrtcs);

IPObjectStruct *User3DDither2Images(const char *Image1Name,
				    const char *Image2Name,
				    int DitherSize,
				    int MatchWidth,
				    int Negate,
				    int AugmentContrast,
				    User3DSpreadType SpreadMethod,
				    IrtRType SphereRad,
				    IrtRType *AccumPenalty);
IPObjectStruct *User3DDither3Images(const char *Image1Name,
				    const char *Image2Name,
				    const char *Image3Name,
				    int DitherSize,
				    int MatchWidth,
				    int Negate,
				    int AugmentContrast,
				    User3DSpreadType SpreadMethod,
				    IrtRType SphereRad,
				    IrtRType *AccumPenalty);

/* Geometry registration. */

int UserRegisterTestConvergance(IrtRType Dist, int i);
IrtRType UserRegisterTwoPointSets(int n1,
				  IrtPtType *PtsSet1,
				  int n2,
				  IrtPtType *PtsSet2,
				  IrtRType AlphaConverge,
				  IrtRType Tolerance,
				  UserRegisterTestConverganceFuncType
				      RegisterTestConvergance,
				  IrtHmgnMatType RegMat);
IrtRType UserRegisterPointSetSrf(int n,
				 IrtPtType *PtsSet,
				 const CagdSrfStruct *Srf,
				 IrtRType AlphaConverge,
				 IrtRType Tolerance,
				 UserRegisterTestConverganceFuncType
				                    RegisterTestConvergance,
				 IrtHmgnMatType RegMat);

/* Bump mapping. */

IPObjectStruct *UserDDMPolysOverTrimmedSrf(const TrimSrfStruct *TSrf,
					   const IPObjectStruct *Texture,
					   IrtRType UDup,
					   IrtRType VDup,
					   int LclUV,
					   int Random);
IPObjectStruct *UserDDMPolysOverSrf(const CagdSrfStruct *Srf,
				    const IPObjectStruct *Texture,
				    IrtRType UDup,
				    IrtRType VDup,
				    int LclUV,
				    int Random);
IPObjectStruct *UserDDMPolysOverPolys(IPObjectStruct *PlSrf,
				      const IPObjectStruct *Texture,
				      IrtRType UDup,
				      IrtRType VDup,
				      int LclUV);

/* Freeform kernels. */

IPObjectStruct *UserSrfKernel(const CagdSrfStruct *Srf,
			      CagdRType SubdivTol,
			      int SkipRate);
IPObjectStruct *UserSrfParabolicLines(const CagdSrfStruct *Srf,
				      CagdRType Step,
				      CagdRType SubdivTol,
				      CagdRType NumericTol,
				      int Euclidean,
				      int DecompSrfs);
IPObjectStruct *UserSrfParabolicSheets(const CagdSrfStruct *Srf,
				       CagdRType Step,
				       CagdRType SubdivTol,
				       CagdRType NumericTol,
				       CagdRType SheetExtent);

/* Freeform umbilical and curvature analysis. */

MvarPtStruct *UserSrfUmbilicalPts(const CagdSrfStruct *Srf,
				  CagdRType SubTol,
				  CagdRType NumTol);
IPObjectStruct *UserSrfFixedCurvatureLines(const CagdSrfStruct *Srf,
					   CagdRType k1,
					   CagdRType Step,
					   CagdRType SubdivTol,
					   CagdRType NumericTol,
					   int Euclidean);
IPObjectStruct *UserCrvCrvtrByOneCtlPt(const CagdCrvStruct *Crv,
				       int CtlPtIdx,
				       CagdRType Min,
				       CagdRType Max,
				       CagdRType SubdivTol,
				       CagdRType NumerTol,
				       int Operation);

/* Polygonal mesh rounding. */

int User2PolyMeshRoundEdge(IPPolygonStruct *Pl1,
			   IPPolygonStruct *Pl2,
			   const IPPolygonStruct *Edge12,
			   IrtRType RoundRadius,
			   IrtRType RoundPower);

/* Image scaling by bivariate spline surface. */

IrtImgPixelStruct *IrtImgScaleImage(IrtImgPixelStruct *InImage,
				    int InMaxX,
				    int InMaxY,
				    int InAlpha,
				    int OutMaxX,
				    int OutMaxY,
				    int Order);

/* Warping of text using composition with surfaces. */

IPObjectStruct *UserWarpTextOnSurface(CagdSrfStruct *Srf,
				      const char *Txt,
				      IrtRType HSpace,
				      IrtRType VBase,
				      IrtRType VTop,
				      IrtRType Ligatures);

/* User inteface to construct piecewise planar cubic Hermite crvs. */

VoidPtr UserHCEditInit(CagdRType StartX, CagdRType StartY, CagdBType Periodic);
VoidPtr UserHCEditFromCurve(const CagdCrvStruct *Crv, CagdRType Tol);
int UserHCEditIsPeriodic(VoidPtr HC);
void UserHCEditSetPeriodic(VoidPtr HC, CagdBType Periodic);
CagdBType UserHCEditGetCtlPtCont(VoidPtr HC, int Index);
void UserHCEditSetCtlPtCont(VoidPtr HC, int Index, CagdBType Cont);
void UserHCEditSetDrawCtlptFunc(VoidPtr HC,
				UserHCEditDrawCtlPtFuncType CtlPtDrawFunc);
void UserHCEditDelete(VoidPtr HC);
VoidPtr UserHCEditCopy(VoidPtr HC);

int UserHCEditTranslate(VoidPtr HC, CagdRType Dx, CagdRType Dy);
int UserHCEditCreateAppendCtlpt(VoidPtr HC,
				CagdRType x,
				CagdRType y,
				int MouseMode);
int UserHCEditCreateDone(VoidPtr HC, CagdRType StartX, CagdRType StartY);
int UserHCEditInsertCtlpt(VoidPtr HC, CagdRType x, CagdRType y, CagdRType t);
int UserHCEditDeleteCtlpt(VoidPtr HC, CagdRType x, CagdRType y);
int UserHCEditUpdateCtl(VoidPtr HC,
			int CtlIndex,
			CagdBType IsPosition,
			CagdRType NewX,
			CagdRType NewY);
int UserHCEditMoveCtl(VoidPtr HC,
		      CagdRType OldX,
		      CagdRType OldY,
		      CagdRType NewX,
		      CagdRType NewY,
		      int MouseMode,
		      CagdRType *MinDist);
int UserHCEditMoveCtlPt(VoidPtr HC, 
			CagdRType OldX,
			CagdRType OldY,
			CagdRType NewX,
			CagdRType NewY,
			int MouseMode);
int UserHCEditMoveCtlTan(VoidPtr HC,
			 CagdRType OldX,
			 CagdRType OldY,
			 CagdRType NewX,
			 CagdRType NewY,
			 int MouseMode);
int UserHCEditIsNearCrv(VoidPtr HC,
			CagdRType x,
			CagdRType y,
			CagdRType *t,
			CagdRType Eps,
			int NormalizeZeroOne);
int UserHCEditIsNearCtlPt(VoidPtr HC,
			  CagdRType *x,
			  CagdRType *y,
			  int *Index,
			  int *UniqueID,
			  CagdRType Eps);
int UserHCEditIsNearCtlTan(VoidPtr HC,
			   CagdRType *x,
			   CagdRType *y,
			   int *Index,
			   int *UniqueID,
			   CagdBType *Forward,
			   CagdRType Eps);
CagdCrvStruct *UserHCEditGetCrvRepresentation(VoidPtr HC, int ArcLen);
int UserHCEditGetCtlPtTan(VoidPtr HC, int Index, CagdPType Pos, CagdPType Tan);
int UserHCEditGetNumCtlPt(VoidPtr HC);
int UserHCEditDrawCtlpts(VoidPtr HC, int DrawTans);
int UserHCEditMatTrans(VoidPtr HC, IrtHmgnMatType Mat);
int UserHCEditTransform(VoidPtr HC, CagdRType *Dir, CagdRType Scl);
int UserHCEditRelativeTranslate(VoidPtr HC, CagdRType *Dir);
int UserHCEditEvalDefTans(VoidPtr HC, int Index);

/* Functions to create machining NC tool path. */

IPObjectStruct *UserNCContourToolPath(const IPObjectStruct *PObj,
				      IrtRType Offset,
				      IrtRType ZBaseLevel,
				      IrtRType Tolerance,
				      IPNCGCodeUnitType Units);
IPObjectStruct *UserNCPocketToolPath(const IPObjectStruct *PObj,
				     IrtRType ToolRadius,
				     IrtRType RoughOffset,
				     IrtRType TPathSpace,
				     IrtRType TPathJoin,
				     IPNCGCodeUnitType Units,
				     int TrimSelfInters);

/* Functions related to finite elements' evaluations. */

UserFEKElementStruct *UserFEKBuildMat(CagdSrfStruct *Srf,
				      int IntegRes,
				      IrtRType E,
				      IrtRType Nu,
				      int *Size);
UserFEKElementStruct *UserFEKBuildMat2(CagdPType *Points,
				       int ULength,
				       int VLength,
				       int UOrder,
				       int VOrder,
				       CagdEndConditionType EndCond,
				       int IntegRes,
				       IrtRType E,
				       IrtRType Nu,
				       int *Size);
CagdBType UserFEPointInsideSrf(CagdSrfStruct *Srf, CagdPType Pt);
UserFEInterIntervalStruct *UserFEGetInterInterval(CagdCrvStruct *Crv1,
						  CagdSrfStruct *Srf1,
						  CagdCrvStruct *Crv2,
						  CagdSrfStruct *Srf2);
UserFECElementStruct *UserFEBuildC1Mat(CagdCrvStruct *Crv1,
				       CagdSrfStruct *Srf1,
				       CagdCrvStruct *Crv2,
				       CagdSrfStruct *Srf2,
				       int IntegRes);
UserFECElementStruct *UserFEBuildC1Mat2(CagdPType *Crv1Pts,
					int Crv1Length,
					int Crv1Order,
					CagdPType *Srf1Pts,
					int Srf1ULength,
					int Srf1VLength,
					int Srf1UOrder,
					int Srf1VOrder,
					CagdPType *Crv2Pts,
					int Crv2Length,
					int Crv2Order,
					CagdPType *Srf2Pts,
					int Srf2ULength,
					int Srf2VLength,
					int Srf2UOrder,
					int Srf2VOrder,
					CagdEndConditionType EndCond,
					int IntegRes);
UserFECElementStruct *UserFEBuildC2Mat(CagdCrvStruct *Crv1,
				       CagdSrfStruct *Srf1,
				       CagdCrvStruct *Crv2,
				       CagdSrfStruct *Srf2,
				       int IntegRes);
UserFECElementStruct *UserFEBuildC2Mat2(CagdPType *Crv1Pts,
					int Crv1Length,
					int Crv1Order,
					CagdPType *Srf1Pts,
					int Srf1ULength,
					int Srf1VLength,
					int Srf1UOrder,
					int Srf1VOrder,
					CagdPType *Crv2Pts,
					int Crv2Length,
					int Crv2Order,
					CagdPType *Srf2Pts,
					int Srf2ULength,
					int Srf2VLength,
					int Srf2UOrder,
					int Srf2VOrder,
					CagdEndConditionType EndCond,
					int IntegRes);
IrtRType UserFEEvalRHSC(UserFECElementStruct *C,
			CagdCrvStruct *Crv1,
			CagdCrvStruct *Crv2);

/* Curve arrangment. */

UserCrvArngmntStruct *UserCrvArngmnt(UserCAOpType Operation,
				     const UserCrvArngmntStruct *CA,
				     const void *Params[]);
UserCrvArngmntStruct *UserCrvArngmntCreate(const IPObjectStruct *PCrvs,
					   CagdRType EndPtEndPtTol,
					   CagdRType PlanarityTol,
					   int ProjectOnPlane,
					   int InputMaskType);
UserCrvArngmntStruct *UserCrvArngmntCopy(const UserCrvArngmntStruct *CA);
int UserCrvArngmntFree(UserCrvArngmntStruct *CA);
int UserCrvArngmntProcessEndPts(UserCrvArngmntStruct *CA);
int UserCrvArngmntClassifyConnectedRegions(UserCrvArngmntStruct *CA,
					   int IgnoreIteriorHangingCrvs);
int UserCrvArngmntIsContained(const UserCrvArngmntStruct *CA,
			      const CagdCrvStruct *InnerShape,
			      const CagdCrvStruct *OuterLoop);
int UserCrvArngmntIsContained2(const UserCrvArngmntStruct *CA,
			       const CagdPType Pt,
			       const CagdCrvStruct *Loop);
CagdCrvStruct *UserCrvArngmntGetCurves(UserCrvArngmntStruct *CA, int XYCurves);
int UserCrvArngmntRegions2Curves(const UserCrvArngmntStruct *CA,
				 int Merge,
				 int XYCurves,
				 IrtRType ZOffset);
int UserCrvArngmntRegionsTopology(const UserCrvArngmntStruct *CA,
				  int XYCurves,
				  IrtRType ZOffset);
void UserCrvArngmntReport(const UserCrvArngmntStruct *CA,
			  int DumpCurves,
			  int DumpPts,
			  int DumpRegions,
			  int DumpXYData);
int UserCrvArngmntOutput(const UserCrvArngmntStruct *CA,
			 int OutputStyle,
			 CagdRType Tolerance,
			 CagdRType ZOffset);

/* Belt curves creation around circles. */

IPObjectStruct *UserBeltCreate(IPVertexStruct *Circs,
			       IrtRType BeltThickness,
			       IrtRType BoundingArcs,
			       int ReturnCrvs,
			       int *Intersects,
			       const char **Error);

/* Ruled surface fitting. */

CagdSrfStruct *UserRuledSrfFit(const CagdSrfStruct *Srf,
			       CagdSrfDirType RulingDir,
			       CagdRType ExtndDmn,
			       int Samples,
			       CagdRType *Error,
			       CagdRType *MaxError);

/* Kinematics simulation using the MV solver. */

int UserKnmtcsSolveMotion(const UserKnmtcsStruct *System,
			  CagdRType NumTol,
			  CagdRType SubTol,
			  CagdRType Step,
			  int *SolDim,
			  CagdBType FilterSols);
int UserKnmtcsNumOfSolPts(int PolyIdx);
void UserKnmtcsEvalAtParams(int PolyIdx, int PtIdx);
CagdCrvStruct *UserKnmtcsEvalCrvTraces();
void UserKnmtcsFreeSol(void);
void UserKnmtcsSolveDone(void);

/* Font processing into curves and geometry. */

#ifdef __WINNT__
/* #define _FREETYPE_FONTS_                 On windows use native windows fonts. */
#define USER_FONT_DEFAULT_WIDTH	-1 

char *UserWChar2Ascii(const UserFontText Str);
UserFontText UserAscii2WChar(const char *Str);
IPObjectStruct *UserFontConvertFontToBezier(
			       const UserFontText Text,
			       const UserFontName FontName,
			       UserFontStyleType FontStyle,
			       IrtRType SpaceWidth,
			       int MergeToBsp,
			       const char *RootObjName);
#else
/* On other systems try to use the freetype library. */
#define _FREETYPE_FONTS_

IPObjectStruct *UserFontFTStringOutline2BezierCurves(
						 const UserFontText Text,
						 const UserFontName FontName,
						 IrtRType Spacing,
						 int MergeToBsp,
						 const char *RootObjName,
						 const char **ErrStr);
#endif /* __WINNT__ */

IPPolygonStruct *UserFontBspCrv2Poly(CagdCrvStruct *BspCrv,
				     IrtRType Tolerance);
CagdCrvStruct *UserFontBzrList2BspList(IPObjectStruct *BzrListObj,
				       IrtBType *HaveHoles);
IPObjectStruct *UserFontBspList2Plgns(IPObjectStruct *BspListObj,
				      IrtRType Tol,
				      const char *Name);
IPObjectStruct *UserFontBspList2TrimSrfs(IPObjectStruct *BspListObj,
					 IrtRType Tol,
					 const char *Name);
IPObjectStruct *UserFontBspList2Solids(IPObjectStruct *BspListObj,
				       UserFont3DEdgeType ExtStyle, 
				       IrtRType ExtOffset, 
				       IrtRType ExtHeight,
				       IrtRType Tol,
				       CagdBType GenTrimSrfs,
				       const char *Name);
int UserFontConvertTextToGeom(const UserFontText Text,
			      const UserFontName FontName,
			      UserFontStyleType FontStyle,
			      IrtRType FontSize,
			      IrtRType TextSpace,
			      UserFont3DEdgeType Text3DEdgeType,
			      const IrtRType Text3DSetup[3],
			      IrtRType Tolerance,
			      UserFontGeomOutputType OutputType,
			      IrtBType CompactOutput,
			      const char *PlacedTextBaseName,
			      IPObjectStruct **PlacedTextGeom,
			      char **ErrorSt);
int UserFontLayoutOverShape(const UserFontText Text,
			    const UserFontName FontName,
			    UserFontStyleType FontStyle,
			    IrtRType FontSize,
			    const IrtRType FontSpace[3],
			    IrtRType Tolerance,
			    UserFont3DEdgeType Text3DEdge,
			    const IrtRType Text3DSetup[2],
			    UserFontAlignmentType Alignment,
			    const IPPolygonStruct *BoundingPoly,
			    UserFontGeomOutputType OutputType,
			    IPObjectStruct **PlacedTextGeom,
			    char **ErrorStr);
int UserFontLayoutOverShape2(const UserFontText Text,
			     const UserFontName FontName,
			     UserFontStyleType FontStyle,
			     IrtRType FontSize,
			     const IrtRType FontSpace[3],
			     IrtRType Tolerance,
			     UserFont3DEdgeType Text3DEdge,
			     const IrtRType Text3DSetup[2],
			     UserFontAlignmentType Alignment,
			     const CagdCrvStruct *BoundingCrv,
			     UserFontGeomOutputType OutputType,
			     IPObjectStruct **PlacedTextGeom,
			     char **ErrorStr);
void UserFontLayoutOverShapeFree(struct UserFontWordLayoutStruct *Words);
struct UserFontWordLayoutStruct *UserFontLayoutOverShapeGenWords(
				    const UserFontText Text,
				    const UserFontName FontName,
				    UserFontStyleType FontStyle,
				    IrtRType FontSize,
				    const IrtRType FontSpace[3],
				    IrtRType Tolerance,
				    UserFont3DEdgeType Text3DEdge,
				    const IrtRType Text3DSetup[2],
				    UserFontAlignmentType Alignment,
				    const IPPolygonStruct *BoundingPoly,
				    UserFontGeomOutputType OutputType,
				    IrtBType CompactOutput,
				    const char *OutputBaseName,
				    UserFontDimInfoStruct *FontDims,
				    char **ErrorStr);
int UserFontLayoutOverShapePlaceWords(struct UserFontWordLayoutStruct *Words,
				      IrtRType FontSize,
				      const IrtRType FontSpace[3],
				      UserFontAlignmentType Alignment,
				      const UserFontDimInfoStruct *FontDims,
				      const IPPolygonStruct *BoundingPoly,
				      IPObjectStruct **PlacedTextGeom);

/* Error handling. */

UserSetErrorFuncType UserSetFatalErrorFunc(UserSetErrorFuncType ErrorFunc);
const char *UserDescribeError(UserFatalErrorType ErrorNum);
void UserFatalError(UserFatalErrorType ErrID);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /* USER_LIB_H */
