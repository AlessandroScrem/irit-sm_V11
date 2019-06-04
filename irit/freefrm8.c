/*****************************************************************************
*   "Irit" - the 3d (not only polygonal) solid modeller.		     *
*									     *
* Written by:  Gershon Elber				Ver 0.2, Mar. 1990   *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
*   Module to provide the required interfact for the cagd library for the    *
* free form surfaces and curves.					     *
*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "program.h"
#include "inc_irit/user_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/ext_lib.h"
#include "inc_irit/triv_iga.h"
#include "freeform.h"

static User3DSpreadType UserCnvrtInt2BlobSpreadMethod(int i);
static UserImgShd3dBlobColorType UserCnvrtInt2BlobColorMethod(int i);

typedef enum {
    IRIT_TRIV_IGA_CALL_ERROR = -1,

    IRIT_TRIV_IGA_CALL_NEW_ARNGMNT = 1,
    IRIT_TRIV_IGA_CALL_NEW_FIELD,
    IRIT_TRIV_IGA_CALL_SET_DEFAULT_SEEDING,
    IRIT_TRIV_IGA_CALL_SET_DEFAULT_DOMAIN,

    IRIT_TRIV_IGA_CALL_ADD_TV,
    IRIT_TRIV_IGA_CALL_EXTRUDE_TV,
    IRIT_TRIV_IGA_CALL_EXTRUDE_TV2,
    IRIT_TRIV_IGA_CALL_TV_OF_REV,
    IRIT_TRIV_IGA_CALL_TV_FROM_SRFS,
    IRIT_TRIV_IGA_CALL_TV_FROM_SRFS2,

    IRIT_TRIV_IGA_CALL_TV_REFINE,
    IRIT_TRIV_IGA_CALL_TV_DEG_RAISE,

    IRIT_TRIV_IGA_CALL_ARNGMNT_COMPLETE,

    IRIT_TRIV_IGA_CALL_PRINT_TV,
    IRIT_TRIV_IGA_CALL_GET_GLBL_MAX_IDS,
    IRIT_TRIV_IGA_CALL_GET_TV_CTLPT_ID_RANGE,
    IRIT_TRIV_IGA_CALL_GET_NUM_BZR_ELEMNTS,
    IRIT_TRIV_IGA_CALL_GET_BZR_ELEMNT_CTLPTS,
    IRIT_TRIV_IGA_CALL_GET_KNOT_INTERVAL,
    IRIT_TRIV_IGA_CALL_MOVE_CTLPT_POS,
    IRIT_TRIV_IGA_CALL_SET_CTLPT_POS,
    IRIT_TRIV_IGA_CALL_TV_EVAL,
    IRIT_TRIV_IGA_CALL_TV_EVAL_BASIS,

    IRIT_TRIV_IGA_CALL_ARNGMNT_FREE,

    IRIT_TRIV_IGA_CALL_GET_ID_BY_TV,
    IRIT_TRIV_IGA_CALL_GET_TV_BY_ID,
    IRIT_TRIV_IGA_CALL_GET_NEIGHBOR_TVS_IDS,
    IRIT_TRIV_IGA_CALL_GET_ALL_TVS_IDS,
    IRIT_TRIV_IGA_CALL_GET_TV,
    IRIT_TRIV_IGA_CALL_GET_TV_FACE_CTLPTS_INDICES,
    IRIT_TRIV_IGA_CALL_GET_TV_FACE_AS_SRF,
    IRIT_TRIV_IGA_CALL_GET_TV_CTLPTS_INDICES,
    IRIT_TRIV_IGA_CALL_GET_CTLPT,
    IRIT_TRIV_IGA_CALL_GET_MATERIAL,

    IRIT_TRIV_IGA_CALL_NEW_MATERIAL,
    IRIT_TRIV_IGA_CALL_EXPORT_TO_XML,
    IRIT_TRIV_IGA_CALL_LOAD_MATERIAL_FROM_XML,

    IRIT_TRIV_IGA_CALL_ADD_BNDRY_FACE,
    IRIT_TRIV_IGA_CALL_ADD_BNDRY_FACE2,
    IRIT_TRIV_IGA_CALL_ADD_BNDRY_FACE_BY_PT,
    IRIT_TRIV_IGA_CALL_GET_BNDRY_FACE_ID_BY_PT,
    IRIT_TRIV_IGA_CALL_GET_LAST_ERROR,
    IRIT_TRIV_IGA_CALL_DESCRIBE_ERROR
} IritTrivIGAFuncCallType;

static CagdCrvStruct *IritMapUVCrv2E3(const CagdCrvStruct *Crv,
				      const CagdSrfStruct *Srf);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Converts an integral integer value to a blob spread method.              *
*                                                                            *
* PARAMETERS:                                                                *
*   i:   Integral integer input.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   User3DSpreadType:    Converted data.                                     *
*****************************************************************************/
static User3DSpreadType UserCnvrtInt2BlobSpreadMethod(int i)
{
    switch (i) {
        default:
        case 1:
	    return USER_3D_SPREAD_RANDOM;
        case 2:
	    return USER_3D_SPREAD_DIAG_PLANE;
	case 3:
	    return USER_3D_SPREAD_DIAG_PLANE2;
        case 4:
	    return USER_3D_SPREAD_ANTI_DIAG_PLANE;
	case 5:
	    return USER_3D_SPREAD_ANTI_DIAG_PLANE2;
	case 6:
	    return USER_3D_SPREAD_ANTI_DIAG_PLANE3;
	case 7:
	    return USER_3D_SPREAD_DISCONT2PLANE;
        case 8:
	    return USER_3D_SPREAD_DISCONT4PLANE;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Converts an integral integer value to a blob color method.               *
*                                                                            *
* PARAMETERS:                                                                *
*   i:   Integral integer input.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   UserImgShd3dBlobColorType:    Converted data.                            *
*****************************************************************************/
static UserImgShd3dBlobColorType UserCnvrtInt2BlobColorMethod(int i)
{
    switch (i) {
        default:
        case 1:
	    return USER_IMG_SHD_3D_BLOB_NO_COLOR;
        case 2:
	    return USER_IMG_SHD_3D_BLOB_GRAY_LEVEL;    
        case 3:
	    return USER_IMG_SHD_3D_BLOB_COLOR;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates Resolution^2 blobs that looks like the 1st image (gray level)    M
* from the XZ plane and like the 2nd image from the YZ plane.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Image1Name:    Name of 1st image to load.                                M
*   Image2Name:    Name of 2nd image to load.                                M
*   RDoTexture:    TRUE to add 'uvvals' attributes, FALSE to add actual      M
*                  color to the vertices of the objels.                      M
*   Blob:          If specified used as the blob.  If NULL, a cross is used. M
*                  Blob must be of size one in each axis, probably centered  M
*                  around the origin.  Must be a list of 3 objects for blob  M
*                  coloring methods other that "No color".		     M
*   RBlobSpread:   Method of spreading the blobs.                            M
*   RBlobColor:    Method of coloring each blob.                             M
*   RResolution:   Resolution of created objects (n^2 objects are created).  M
*   RNegative:     Default (FALSE) is white blobs over dark background. If   M
*                  TRUE, assume dark blobs over white background.            M
*   Intensity:     The gray level affect on the blobs' scale.                M
*   MinIntensity:  Minimum intensity allowed.  Zero will collapse the blobl  M
*                  completely in one direction which will render it          M
*                  impossible to actually manufacture it.		     M
*   RMergePolys:   TRUE to merge all objects' polygons into one object.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A list (poly if MergePolys) object of Resolution^2    M
*		       spherical blobs.		                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserMake3DStatueFrom2Images                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   Model3DFrom2Images                                                       M
*****************************************************************************/
IPObjectStruct *Model3DFrom2Images(const char *Image1Name,
				   const char *Image2Name,
				   IrtRType *RDoTexture,
				   const IPObjectStruct *Blob,
				   IrtRType *RBlobSpread,
				   IrtRType *RBlobColor,
				   IrtRType *RResolution,
				   IrtRType *RNegative,
				   IrtRType *Intensity,
				   IrtRType *MinIntensity,
				   IrtRType *RMergePolys)
{
    User3DSpreadType
        BlobSpreadMethod = UserCnvrtInt2BlobSpreadMethod(
					 IRIT_REAL_PTR_TO_INT(RBlobSpread));
    UserImgShd3dBlobColorType
        BlobColorMethod = UserCnvrtInt2BlobColorMethod(
					 IRIT_REAL_PTR_TO_INT(RBlobColor));
	       
    if (!IP_IS_POLY_OBJ(Blob) && !IP_IS_OLST_OBJ(Blob))
        Blob = NULL;

    return UserMake3DStatueFrom2Images(Image1Name, Image2Name,
				       IRIT_REAL_PTR_TO_INT(RDoTexture), Blob,
				       BlobSpreadMethod, BlobColorMethod,
				       IRIT_REAL_PTR_TO_INT(RResolution),
				       IRIT_REAL_PTR_TO_INT(RNegative),
				       *Intensity, *MinIntensity,
				       IRIT_REAL_PTR_TO_INT(RMergePolys));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates Resolution^2 3D cross blobs that looks like the 1st image        M
* (gray level) from the XZ plane, like the 2nd image from the YZ plane, and  M
* like the 3rd image from the XY plane.					     M
*   A 3D cross is used as a blob.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Image1Name:    Name of 1st image to load.                                M
*   Image2Name:    Name of 2nd image to load.                                M
*   Image3Name:    Name of 3rd image to load.                                M
*   RDoTexture:    TRUE to add 'uvvals' attributes, FALSE to add actual      M
*                  color to the vertices of the objels.                      M
*   Blob:          If specified used as the blob. If NULL, a 3-cross is used.M
*                  Blob must be of size one in each axis, probably centered  M
*                  around the origin.  Should hold 3 polygons for "No color" M
*		   blob color method and should hold a list of 3 objects for M
*                  the other methods. 					     M
*   RBlobSpread:   Method of spreading the blobs.                            M
*   RBlobColor:    Method of coloring each blob.                             M
*   RResolution:   Resolution of created objects (n^2 objects are created).  M
*   RNegative:     Default (FALSE) is white blobs over dark background. If   M
*                  TRUE, assume dark blobs over white background.            M
*   Intensity:     The gray level affect on the blobs' scale.                M
*   MinIntensity:  Minimum intensity allowed.  Zero will collapse the blobl  M
*                  completely in one direction which will render it          M
*                  impossible to actually manufacture it.		     M
*   RMergePolys:   TRUE to merge all objects' polygons into one object.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A list (poly if MergePolys) object of Resolution^2    M
*		       spherical blobs.		                             M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserMake3DStatueFrom3Images                                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   Model3DFrom3Images                                                       M
*****************************************************************************/
IPObjectStruct *Model3DFrom3Images(const char *Image1Name,
				   const char *Image2Name,
				   const char *Image3Name,
				   IrtRType *RDoTexture,
				   const IPObjectStruct *Blob,
				   IrtRType *RBlobSpread,
				   IrtRType *RBlobColor,
				   IrtRType *RResolution,
				   IrtRType *RNegative,
				   IrtRType *Intensity,
				   IrtRType *MinIntensity,
				   IrtRType *RMergePolys)
{
    User3DSpreadType
        BlobSpreadMethod = UserCnvrtInt2BlobSpreadMethod(
					 IRIT_REAL_PTR_TO_INT(RBlobSpread));
    UserImgShd3dBlobColorType
        BlobColorMethod = UserCnvrtInt2BlobColorMethod(
					 IRIT_REAL_PTR_TO_INT(RBlobColor));

    if (!IP_IS_POLY_OBJ(Blob) && !IP_IS_OLST_OBJ(Blob))
        Blob = NULL;

    return UserMake3DStatueFrom3Images(Image1Name, Image2Name, Image3Name,
				       IRIT_REAL_PTR_TO_INT(RDoTexture),
				       Blob, BlobSpreadMethod, BlobColorMethod,
				       IRIT_REAL_PTR_TO_INT(RResolution),
				       IRIT_REAL_PTR_TO_INT(RNegative),
				       *Intensity, *MinIntensity,
				       IRIT_REAL_PTR_TO_INT(RMergePolys));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates Resolution^2 blobs that looks like the 1st image (gray level)    M
* from the XZ plane and like the 2nd image from the YZ plane.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Image1Name:  Name of 1st image to load.                                  M
*   Image2Name:  Name of 2nd image to load.                                  M
*   RDitherSize:  2, 3 or 4 for (2x2), (3x3) or (4x4) dithering.             M
*   RMatchWidth:  Width to allow matching in a row:			     M
*                           between pos[i] to pos[i +/- k],  k < MatchWidth. M
*   RNegate:     TRUE to negate the images before use.			     M
*   RAugmentContrast:  Number of iterations to add micro-pixels, to augment  M
*                the contrast, behind existing pixels.  Zero to disable.     M
*   SpreadMethod: If allowed (MatchWidth >= RowSize), selects initial random M
*                spread to use.						     M
*   SphereRad:   Radius of construct spherical blob, zero to return points.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  Micro blobs if SphereRad > 0, center points, if = 0.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserDither3D2Images                                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   MicroBlobsDither3DFrom2Images                                            M
*****************************************************************************/
IPObjectStruct *MicroBlobsDither3DFrom2Images(const char *Image1Name,
					      const char *Image2Name,
					      IrtRType *RDitherSize,
					      IrtRType *RMatchWidth,
					      IrtRType *RNegate,
					      IrtRType *RAugmentContrast,
					      IrtRType *RSpreadMethod,
					      IrtRType *SphereRad)
{
    IrtRType AccumPenalty;
    User3DSpreadType
        SpreadMethod = UserCnvrtInt2BlobSpreadMethod(
					 IRIT_REAL_PTR_TO_INT(RSpreadMethod));
    IPObjectStruct *PObj;

    PrimSetResolution(GetResolution(TRUE));

    PObj = User3DDither2Images(Image1Name, Image2Name,
			       IRIT_REAL_PTR_TO_INT(RDitherSize),
			       IRIT_REAL_PTR_TO_INT(RMatchWidth),
			       IRIT_REAL_PTR_TO_INT(RNegate),
			       IRIT_REAL_PTR_TO_INT(RAugmentContrast),
			       SpreadMethod, *SphereRad, &AccumPenalty);

    if (PObj == NULL)
	IRIT_NON_FATAL_ERROR("BFrom2Img: invalid input (images of different sizes!?)");
    
    IRIT_WNDW_FPRINTF2("BFrom2Img: total computed penalty = %f", AccumPenalty);

    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates Resolution^2 blobs that looks like the 1st image (gray level)    M
* from the XZ plane, like the 2nd image from the YZ plane, and like the 3nd  M
* image from the XY plane.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Image1Name:  Name of 1st image to load.                                  M
*   Image2Name:  Name of 2nd image to load.                                  M
*   Image3Name:  Name of 3rd image to load.                                  M
*   RDitherSize:  2, 3, or 4 for (2x2x2), (3x3x3), or (3x3x3) dithering.     M
*   RMatchWidth:  Width to allow matching in a row:			     M
*                           between pos[i] to pos[i +/- k],  k < MatchWidth. M
*   RNegate:     TRUE to negate the images before use.			     M
*   RAugmentContrast:  Number of iterations to add micro-pixels, to augment  M
*                the contrast, behind existing pixels.  Zero to disable.     M
*   SpreadMethod: If allowed (MatchWidth >= RowSize), selects initial random M
*                spread to use.						     M
*   SphereRad:   Radius of construct spherical blob, zero to return points.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  Micro blobs if SphereRad > 0, center points, if = 0.  M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserDither3D2Images                                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   MicroBlobsDither3DFrom2Images                                            M
*****************************************************************************/
IPObjectStruct *MicroBlobsDither3DFrom3Images(const char *Image1Name,
					      const char *Image2Name,
					      const char *Image3Name,
					      IrtRType *RDitherSize,
					      IrtRType *RMatchWidth,
					      IrtRType *RNegate,
					      IrtRType *RAugmentContrast,
					      IrtRType *RSpreadMethod,
					      IrtRType *SphereRad)
{
    IrtRType AccumPenalty;
    User3DSpreadType
        SpreadMethod = UserCnvrtInt2BlobSpreadMethod(
					 IRIT_REAL_PTR_TO_INT(RSpreadMethod));
    IPObjectStruct *PObj;

    PrimSetResolution(GetResolution(TRUE));

    PObj = User3DDither3Images(Image1Name, Image2Name, Image3Name,
			       IRIT_REAL_PTR_TO_INT(RDitherSize),
			       IRIT_REAL_PTR_TO_INT(RMatchWidth),
			       IRIT_REAL_PTR_TO_INT(RNegate),
			       IRIT_REAL_PTR_TO_INT(RAugmentContrast),
			       SpreadMethod, *SphereRad, &AccumPenalty);

    if (PObj == NULL)
	IRIT_NON_FATAL_ERROR("BFrom3Img: invalid input (images of different sizes!?)");
    else
        IRIT_WNDW_FPRINTF2("BFrom3Img: total computed penalty = %f",
			   AccumPenalty);

    AccumPenalty = 0;

    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Fit a ruled surface to the given general surface Srf.  The best fit is   M
* found using a dynmaic programming search over all possibly rulings, while  M
* each ruling line's distance is measured agains the surface.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   SrfObj:     To fit a ruled surface through.                              M
*   RDir:       Either the U or the V direction.  This is used only to       M
*               sample Srf and construct the possibly ruling lines.          M
*   ExtndDmn:   Amount to extended the selected sampled boundary curves.     M
*               Zero will not extend and match the ruling from the original  M
*               boundary to its maximum.				     M
*   RSamples:   Number of samples to compute the dynamic programming with.   M
*		Typically in the hundreds.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:    The fitted ruled surface, or NULL if error.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserRuledSrfFit                                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   RuledSrfFit                                                              M
*****************************************************************************/
IPObjectStruct *RuledSrfFit(IPObjectStruct *SrfObj,
			    IrtRType *RDir,
			    IrtRType *ExtndDmn,
			    IrtRType *RSamples)
{
    CagdRType Error, MaxError;
    CagdSrfStruct
        *RuledSrf = UserRuledSrfFit(SrfObj -> U.Srfs,
				    (CagdSrfDirType) IRIT_REAL_PTR_TO_INT(RDir),
				    *ExtndDmn,
				    IRIT_REAL_PTR_TO_INT(RSamples),
				    &Error,
				    &MaxError);

    if (RuledSrf == NULL)
        return NULL;
    else {
        IPObjectStruct
	    *PObj = IPGenSRFObject(RuledSrf);

	AttrSetObjectRealAttrib(PObj, "Error", Error);
	AttrSetObjectRealAttrib(PObj, "MaxError", MaxError);
	return PObj;
    }    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the orthogonal projection of the given curve on given surface.  M
*                                                                            *
* PARAMETERS:                                                                M
*   CrvObj:       The curve to project on SrfObj, orthogonally.		     M
*   SrfObj:       The surface to project CrvObj on.			     M
*   Tol:          Tolerance of the computation.			             M
*   Euclidean:    TRUE to return the curves in Euclidean space, FALSE in the M
*		  surface parametric domain.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectSTruct *:    Projected curves.  Can be several!                  M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVOrthoCrvProjOnSrf                                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   CrvOrthoProjOnSrf                                                        M
*****************************************************************************/
IPObjectStruct *CrvOrthoProjOnSrf(IPObjectStruct *CrvObj,
				  IPObjectStruct *SrfObj,
				  IrtRType *Tol,
				  IrtRType *Euclidean)
{
    IPPolygonStruct *Pl;
    const CagdSrfStruct
        *Srf = SrfObj -> U.Srfs;
    MvarPolylineStruct
        *Plls = MvarMVOrthoCrvProjOnSrf(CrvObj -> U.Crvs, Srf, *Tol);
    IPObjectStruct
        *PllsObj = MvarCnvrtMVPolysToIritPolys2(Plls, TRUE);
    CagdCrvStruct *Crvs;

    MvarPolylineFreeList(Plls);

    if (!IRIT_APX_EQ(*Euclidean, 0.0)) {         /* Evaluate to E3 points. */
        for (Pl = PllsObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	    IPVertexStruct *V;

	    for (V = Pl -> PVertex; V != NULL; V = V -> Pnext) {
	        CagdRType
		    *R = CagdSrfEval(Srf, V -> Coord[0], V -> Coord[1]);

		CagdCoerceToE3(V -> Coord, &R, -1, Srf -> PType);
	    }
	}
    }

    Crvs = UserPolylines2LinBsplineCrvs(PllsObj -> U.Pl, TRUE);
    IPFreeObject(PllsObj);

    return IPLnkListToListObject(Crvs, IP_OBJ_CURVE);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes tiling with rectangular regions to the domain bound by the      M
* given set of curves.                                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   CrvList:            List of curves that bound the domain.                M
*   RAngularDeviations: Maximal angular deviations, in degrees, on the       M
*                       boundary before forcing a split.		     M
*   RCurveOutputType:   1 for Polygonal rectangles,                          M
*                       2 for curves lists,				     M
*                       3 for surfaces.					     M
*   RSizeRectangle:     Approximated edge size of expected rectangles.       M
*   RNumSmoothingSteps: Low pass (smoothing) filtering to apply to the       M
*                       vertices, with respect to neighboring vertices.      M
*                       Number of smoothing iteration - zero to disable.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A list object that holds polygons, or curves, or      M
*                      surfaces, depending on RCurveOutputType.		     M
*                      NULL is returned in case of error (and ErrorMsg set). M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtExtC2SGeneral                                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   Crvs2RectRegionGeneral                                                   M
*****************************************************************************/
IPObjectStruct *Crvs2RectRegionGeneral(IPObjectStruct *CrvList,
				       IrtRType *RAngularDeviations,
				       IrtRType *RCurveOutputType,
				       IrtRType *RSizeRectangle,
				       IrtRType *RNumSmoothingSteps)
{
    int i = 0;
    const char *ErrMsg;
    CagdCrvStruct
        *Crvs = NULL;
    IPObjectStruct *PCrv, *PObj;

    while ((PCrv = IPListObjectGet(CrvList, i++)) != NULL) {
        if (IP_IS_CRV_OBJ(PCrv)) {
	    Crvs = CagdListAppend(Crvs, CagdCrvCopyList(PCrv -> U.Crvs));
        }
    }
    if (Crvs == NULL) {
        IRIT_NON_FATAL_ERROR("C2RectRgn: No curves were found in the input");
	return NULL;
    }

    PObj = IrtExtC2SGeneral(&Crvs, IRIT_REAL_PTR_TO_INT(RAngularDeviations),
			    *RSizeRectangle,
			    IRIT_REAL_PTR_TO_INT(RCurveOutputType),
			    *RSizeRectangle,
			    IRIT_REAL_PTR_TO_INT(RNumSmoothingSteps),
			    &ErrMsg, "C2RR");
    CagdCrvFreeList(Crvs);

    if (ErrMsg != NULL)
        IRIT_NON_FATAL_ERROR(ErrMsg);

    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes an extended freeform so the current domain is identical to the  M
* input freeform.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   FFObj:        The freeform to extend beyond its current domain.          M
*                 Either a curve or a surface.                               M
*   ExtDirs:      A list of two (for curves as Min, Max) or four (for        M
*                 as UMin, VMin, UMax, VMax) of boolean values to request    M
*                 extensions along these boundaries.                         M
*   ExtEps:       A list of one (for curve) or two (for surfaces for U and   M
*                 V) of the extension requested, in the parameteric domain.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectSTruct *:    Extended freeform.			             M
*                                                                            *
* SEE ALSO:                                                                  M
*   BspCrvExtension, BspSrfExtension                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   FreeformDomainExtend                                                     M
*****************************************************************************/
IPObjectStruct *FreeformDomainExtend(IPObjectStruct *FFObj,
				     IPObjectStruct *ExtDirs,
				     IPObjectStruct *ExtEps,
				     IrtRType *RemoveExtraKnots)
{
    IPObjectStruct *RetVal;

    if (IP_IS_CRV_OBJ(FFObj)) {
        CagdBType ExtDirsVec[2];
        CagdCrvStruct *ECrv;
	IPObjectStruct *ExtDir0, *ExtDir1, *ExtEps0;

	if (IPListObjectLength(ExtDirs) != 2 ||
	    IPListObjectLength(ExtEps) != 1) {
	    IRIT_NON_FATAL_ERROR("FFExtend: Expected two elements in dir extension and one in epsilons.");
	    return NULL;
	}

	ExtDir0 = IPListObjectGet(ExtDirs, 0);
	ExtDir1 = IPListObjectGet(ExtDirs, 1);
	if (!IP_IS_NUM_OBJ(ExtDir0) || !IP_IS_NUM_OBJ(ExtDir1)) {
	    IRIT_NON_FATAL_ERROR("FFExtend: Expected two numeric values in ExtDirs.");
	    return NULL;
	}
	ExtDirsVec[0] = !IRIT_APX_EQ(ExtDir0 -> U.R, 0.0);
	ExtDirsVec[1] = !IRIT_APX_EQ(ExtDir1 -> U.R, 0.0);

	ExtEps0 = IPListObjectGet(ExtEps, 0);
	if (!IP_IS_NUM_OBJ(ExtEps0)) {
	    IRIT_NON_FATAL_ERROR("FFExtend: Expected a numeric value in ExtEps.");
	    return NULL;
	}

	ECrv = BspCrvExtension(FFObj -> U.Crvs, ExtDirsVec,
			       ExtEps0 -> U.R,
			       IRIT_REAL_PTR_TO_INT(RemoveExtraKnots));

	RetVal = IPGenCRVObject(ECrv);
    }
    else if (IP_IS_SRF_OBJ(FFObj)) {
        CagdBType ExtDirsVec[4];
        CagdSrfStruct *ESrf;
	IPObjectStruct *ExtDir0, *ExtDir1, *ExtDir2, *ExtDir3, *ExtEps0, *ExtEps1;

	if (IPListObjectLength(ExtDirs) != 4 ||
	    IPListObjectLength(ExtEps) != 2) {
	    IRIT_NON_FATAL_ERROR("FFExtend: Expected four elements in dir extension and two in epsilons.");
	    return NULL;
	}

	ExtDir0 = IPListObjectGet(ExtDirs, 0);
	ExtDir1 = IPListObjectGet(ExtDirs, 1);
	ExtDir2 = IPListObjectGet(ExtDirs, 2);
	ExtDir3 = IPListObjectGet(ExtDirs, 3);
	if (!IP_IS_NUM_OBJ(ExtDir0) ||
	    !IP_IS_NUM_OBJ(ExtDir1) ||
	    !IP_IS_NUM_OBJ(ExtDir2) ||
	    !IP_IS_NUM_OBJ(ExtDir3)) {
	    IRIT_NON_FATAL_ERROR("FFExtend: Expected four numeric values in ExtDirs.");
	    return NULL;
	}
	ExtDirsVec[0] = !IRIT_APX_EQ(ExtDir0 -> U.R, 0.0);
	ExtDirsVec[1] = !IRIT_APX_EQ(ExtDir1 -> U.R, 0.0);
	ExtDirsVec[2] = !IRIT_APX_EQ(ExtDir2 -> U.R, 0.0);
	ExtDirsVec[3] = !IRIT_APX_EQ(ExtDir3 -> U.R, 0.0);

	ExtEps0 = IPListObjectGet(ExtEps, 0);
	ExtEps1 = IPListObjectGet(ExtEps, 1);
	if (!IP_IS_NUM_OBJ(ExtEps0) || !IP_IS_NUM_OBJ(ExtEps1)) {
	    IRIT_NON_FATAL_ERROR("FFExtend: Expected two numeric values in ExtEps.");
	    return NULL;
	}

	ESrf = BspSrfExtension(FFObj -> U.Srfs, ExtDirsVec,
			       ExtEps0 -> U.R, ExtEps1 -> U.R,
			       IRIT_REAL_PTR_TO_INT(RemoveExtraKnots));

	RetVal = IPGenSRFObject(ESrf);
    }
    else {
        assert(0); /* Should never get here. */
	return NULL;
    }

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A wrapper routine to read an image from a file, dither it and save into  M
* a file the dithered image.                                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   InputImage:     To dither.						     M
*   OutputImage:    Dithered image.                                          M
*   DitherSize:     Dithering matrices size: 2, 3, or 4.                     M
*   ErrorDiffusion: TRUE, to also diffuse the error in the image.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgDitherImage2                                                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritDitherImage                                                          M
*****************************************************************************/
void IritDitherImage(const char *InputImage,
		     const char *OututImage,
		     IrtRType *DitherSize,
		     IrtRType *ErrorDiffusion)
{
    if (!IrtImgDitherImage2(InputImage, OututImage,
			    IRIT_REAL_PTR_TO_INT(DitherSize),
			    IRIT_REAL_PTR_TO_INT(ErrorDiffusion)))
	IRIT_NON_FATAL_ERROR("DitherImage: Failed to dither input image.");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A wrapper routine to read an image from a file, dither it and save into  M
* a file the dithered image.                                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   PPolyObj:       A polygonal object to apply a subdivision scheme to.     M
*   RSubdivScheme:  The subdivision scheme to apply:			     M
*                   0 for Catmull Clark,				     M
*                   1 for Loop Scheme,					     M
*                   2 for Butterfly scheme.				     M
*   RNumIterations: Dithering matrices size: 2, 3, or 4.                     M
*   RSmoothNormals: Number of subdivision iterations to apply. At least one. M
*   RTrianglesOnly: TRUE for triangles only in the output.		     M
*   AdditionalParam:  As required by the different schemes:		     M
*                   For Butterfly scheme - the W coefficient.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  Refined polygonal model, using selected subdivision   M
*                      scheme.						     M
*                                                                            *
* SEE ALSO:                                                                  M
*   GMSubCatmullClark, GMSubLoop, GMSubButterfly                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritSubdivSurface                                                        M
*****************************************************************************/
IPObjectStruct *IritSubdivSurface(IPObjectStruct *PPolyObj,
				  IrtRType *RSubdivScheme,
				  IrtRType *RNumIterations,
				  IrtRType *RSmoothNormals,
				  IrtRType *RTrianglesOnly,
				  IrtRType *AdditionalParam)
{
    IPObjectStruct *TempObject,
        *ResultObject = NULL;
    int i,
	SubdivScheme = IRIT_REAL_PTR_TO_INT(RSubdivScheme),
        NumIterations = IRIT_REAL_PTR_TO_INT(RNumIterations),
        SmoothNormals = IRIT_REAL_PTR_TO_INT(RSmoothNormals),
        TrianglesOnly = IRIT_REAL_PTR_TO_INT(RTrianglesOnly);

    switch (SubdivScheme) {
	default:
	    assert(0);
	    IRIT_NON_FATAL_ERROR(
	        "PSubdiv: An invalid selected subdivision scheme ignored.");
	    return NULL;
	case 0: /* Catmull Clark. */
	    TempObject = PPolyObj;
	    for (i = 0; i < NumIterations; i++) {
	        ResultObject = GMSubCatmullClark(TempObject);
	        if (TempObject != PPolyObj)
		    IPFreeObject(TempObject);
		TempObject = ResultObject;
	    }
	    break;
	case 1: /* Loop. */
	    TempObject = PPolyObj;
	    for (i = 0; i < NumIterations; i++) {
	        ResultObject = GMSubLoop(TempObject);
	        if (TempObject != PPolyObj)
	    	    IPFreeObject(TempObject);
		TempObject = ResultObject;
	    }
	    break;
	case 2: /* Butterfly. */
	    TempObject = PPolyObj;
	    for (i = 0; i < NumIterations; i++) {
	        ResultObject = GMSubButterfly(TempObject, *AdditionalParam);
	        if (TempObject != PPolyObj)
	    	    IPFreeObject(TempObject);
		TempObject = ResultObject;
	    }
	    break;
    }

    if (TrianglesOnly) {
        IPObjectStruct
	    *PTmp = GMConvertPolysToTriangles(ResultObject);

	IPFreeObject(ResultObject);
	ResultObject = PTmp;
    }

    if (SmoothNormals)		     /* Compute smooth normals to vertices. */
        GMBlendNormalsToVertices(ResultObject -> U.Pl, 180);	
    else {	                                           /* Flat Normals. */
        IPPolygonStruct *Pl;
	IPVertexStruct *V;
	for (Pl = ResultObject -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	    V = Pl -> PVertex;
	    do {
	        IRIT_VEC_COPY(V -> Normal, Pl -> Plane);
		IP_SET_NORMAL_VRTX(V);
		V = V -> Pnext;
	    }
	    while (V != NULL && V != Pl -> PVertex);
	}
    }

    GMConvexPolyObject(ResultObject);

    return ResultObject;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Implements an interpreter interface to the Triv lib IGA code.            M
*                                                                            M
*                                                                            *
* PARAMETERS:                                                                M
*   CallNum:   Index of triv IGA function.                                   M
*   ParamList: The list of parameter of the TGA function.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  Return value of the IGA function.                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAInterface                                                         M
*****************************************************************************/
IPObjectStruct *TrivIGAInterface(IrtRType *CallNum, IPObjectStruct *ParamList)
{
    IRIT_STATIC_DATA int
        Talkative = FALSE;
    int i, *IDs,
        IGAFunc = IRIT_REAL_PTR_TO_INT(CallNum),
        ID = -1;
    const CagdCtlPtStruct *CtlPt;
    TrivTVStruct *TV;
    IPObjectStruct *Prm1, *Prm2, *Prm3, *Prm4, *Prm5, *Prm6, *Prm7, *Prm8,
        *Prm9, *TObj, *TObj2, *PTmp;

    if (Talkative) {
        fprintf(stderr, "TRIV IGA (%d) invoked\n", IGAFunc);
    }

    switch (IGAFunc) {
	case IRIT_TRIV_IGA_CALL_NEW_ARNGMNT:
	    if (IPListObjectLength(ParamList) == 1 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		IP_IS_NUM_OBJ(Prm1)) {
	        Talkative = IRIT_REAL_TO_INT(Prm1 -> U.R);
		ID = (int) TrivIGANewArrangement(&ID);
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (New Arrangement): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_NEW_FIELD:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_STR_OBJ(Prm2))
		ID = (int) TrivIGANewField(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					   Prm2 -> U.Str);
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (New Field): Invalid parameters.");
	    break;
        case IRIT_TRIV_IGA_CALL_SET_DEFAULT_SEEDING:
	    if (IPListObjectLength(ParamList) == 4 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4)) {
	        TrivTVDirType
		    Dir = TRIV_INT_TO_DIR(IRIT_REAL_TO_INT(Prm2 -> U.R));

	        ID = (int) TrivIGASetDefaultSeeding(IRIT_REAL_TO_INT(Prm1 -> U.R), 
						    Dir, Prm3 -> U.R,
						    IRIT_REAL_TO_INT(Prm4 -> U.R));
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Set Default Seeding): Invalid parameters.");
	    break;
        case IRIT_TRIV_IGA_CALL_SET_DEFAULT_DOMAIN:
	    if (IPListObjectLength(ParamList) == 4 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4)) {
	        TrivTVDirType
		    Dir = TRIV_INT_TO_DIR(IRIT_REAL_TO_INT(Prm2 -> U.R));

	        ID = (int) TrivIGASetDefaultDomain(IRIT_REAL_TO_INT(Prm1 -> U.R), 
						    Dir, Prm3 -> U.R,
						    Prm4 -> U.R);
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Set Default Seeding): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_ADD_TV:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL  &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_TRIVAR_OBJ(Prm2)) {
	        ID = AttrGetObjectIntAttrib(Prm2, "ID");
		if (IP_ATTR_IS_BAD_INT(ID))
		    ID = -1;
	        ID = TrivIGAAddTrivar(IRIT_REAL_TO_INT(Prm1 -> U.R),
				      Prm2 -> U.Trivars, ID);
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (New Trivar): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_EXTRUDE_TV:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_SRF_OBJ(Prm2) &&
		IP_IS_VEC_OBJ(Prm3)) {
	        ID = AttrGetObjectIntAttrib(Prm2, "ID");
		if (IP_ATTR_IS_BAD_INT(ID))
		    ID = -1;
	        ID = (int) TrivIGAExtrudeTV(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					    Prm2 -> U.Srfs, Prm3 -> U.Vec, ID);
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Extruded TV): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_EXTRUDE_TV2:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_SRF_OBJ(Prm2) &&
		IP_IS_CRV_OBJ(Prm3)) {
	        ID = AttrGetObjectIntAttrib(Prm2, "ID");
		if (IP_ATTR_IS_BAD_INT(ID))
		    ID = -1;
	        ID = (int) TrivIGAExtrudeTV2(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					     Prm2 -> U.Srfs, Prm3 -> U.Crvs,
					     ID);
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Extruded TV): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_TV_OF_REV:
	    if (IPListObjectLength(ParamList) == 7 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		(Prm6 = IPListObjectGet(ParamList, 5)) != NULL &&
		(Prm7 = IPListObjectGet(ParamList, 6)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_SRF_OBJ(Prm2) &&
		IP_IS_POINT_OBJ(Prm3) &&
		IP_IS_VEC_OBJ(Prm4) &&
		IP_IS_NUM_OBJ(Prm5) &&
		IP_IS_NUM_OBJ(Prm6) &&
		IP_IS_NUM_OBJ(Prm7)) {
	        ID = AttrGetObjectIntAttrib(Prm2, "ID");
		if (IP_ATTR_IS_BAD_INT(ID))
		    ID = -1;
	        ID = (int) TrivIGATVofRevol(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					    Prm2 -> U.Srfs,
					    Prm3 -> U.Pt,
					    Prm4 -> U.Vec,
					    Prm5 -> U.R,
					    Prm6 -> U.R,
					    Prm7 -> U.R != 0.0,
					    ID);
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (TV of Rev.): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_TV_FROM_SRFS:
	    if (IPListObjectLength(ParamList) == 4 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_OLST_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4)) {
	        CagdSrfStruct *Srf,
		    *Srfs = NULL;

		for (i = 0; (TObj = IPListObjectGet(Prm2, i++)) != NULL; ) {
		    if (IP_IS_SRF_OBJ(TObj)) {
		        Srf = CagdSrfCopy(TObj -> U.Srfs);
			IRIT_LIST_PUSH(Srf, Srfs);
		    }
		    else {
		        IRIT_NON_FATAL_ERROR("TRIVIGA (TV from Srfs.): Expected srfs in srf list.");
			CagdSrfFreeList(Srfs);
			return NULL;
		    }
		}

	        ID = AttrGetObjectIntAttrib(Prm2, "ID");
		if (IP_ATTR_IS_BAD_INT(ID))
		    ID = -1;
		ID = (int) TrivIGATVFromSurfaces(IRIT_REAL_TO_INT(Prm1 -> U.R),
						 Srfs, IRIT_REAL_TO_INT(Prm3 -> U.R),
						 Prm4 -> U.R != 0.0, ID);
		CagdSrfFreeList(Srfs);
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (TV from Srfs): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_TV_FROM_SRFS2:
	    if (IPListObjectLength(ParamList) == 5 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_SRF_OBJ(Prm2) &&
		IP_IS_OLST_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_NUM_OBJ(Prm5)) {
	        IrtHmgnMatType
		    *Mats = IritMalloc(IPListObjectLength(Prm3) *
				                     sizeof(IrtHmgnMatType));

		for (i = 0; (TObj = IPListObjectGet(Prm3, i)) != NULL; i++) {
		    if (IP_IS_MAT_OBJ(TObj)) {
		        IRIT_HMGN_MAT_COPY(Mats[i], TObj -> U.Mat);
		    }
		    else {
		        IRIT_NON_FATAL_ERROR("TRIVIGA (TV from Srfs 2.): Expected transforms in list.");
			IritFree(Mats);
			return NULL;
		    }
		}

	        ID = AttrGetObjectIntAttrib(Prm2, "ID");
		if (IP_ATTR_IS_BAD_INT(ID))
		    ID = -1;
		ID = (int) TrivIGATVFromSurfaces2(IRIT_REAL_TO_INT(Prm1 -> U.R), 
						  Prm2 -> U.Srfs, Mats,
						  IPListObjectLength(Prm3),
						  IRIT_REAL_TO_INT(Prm4 -> U.R),
						  IRIT_REAL_TO_INT(Prm5 -> U.R), ID);
		IritFree(Mats);
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA: Invalid parameters.");
	    break;

	case IRIT_TRIV_IGA_CALL_TV_REFINE:
	    if (IPListObjectLength(ParamList) == 4 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4)) {
	        TrivTVDirType
		    Dir = TRIV_INT_TO_DIR(IRIT_REAL_TO_INT(Prm3 -> U.R));

	        ID = TrivIGATVRefine(IRIT_REAL_TO_INT(Prm1 -> U.R), 
				     IRIT_REAL_TO_INT(Prm2 -> U.R),
			  	     Dir, Prm4 -> U.R) == NULL;
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Refine): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_TV_DEG_RAISE:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3)) {
	        TrivTVDirType
		    Dir = TRIV_INT_TO_DIR(IRIT_REAL_TO_INT(Prm3 -> U.R));

	        ID = TrivIGATDegreeRaise(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					 IRIT_REAL_TO_INT(Prm2 -> U.R),
					 Dir) == NULL;
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Deg. Raise): Invalid parameters.");
	    break;

	case IRIT_TRIV_IGA_CALL_ARNGMNT_COMPLETE:
	    if (IPListObjectLength(ParamList) == 1 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		IP_IS_NUM_OBJ(Prm1))
		ID = (int) TrivIGAArrangementComplete(IRIT_REAL_TO_INT(Prm1 -> U.R));
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Arrangement Complete): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_PRINT_TV:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
					        IRIT_REAL_TO_INT(Prm2 -> U.R));

		ID = (int) TrivIGAPrintTVContent(IRIT_REAL_TO_INT(Prm1 -> U.R),
						 TV);
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Print TV Content): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_GLBL_MAX_IDS:
	    if (IPListObjectLength(ParamList) == 1 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		IP_IS_NUM_OBJ(Prm1)) {
	        int *IDs = TrivIGAGetGlblMaxIDs(IRIT_REAL_TO_INT(Prm1 -> U.R));

		TObj = IPGenLISTObject(IPGenNUMValObject(IDs[0]));
		IPListObjectAppend(TObj, IPGenNUMValObject(IDs[1]));
		IPListObjectAppend(TObj, IPGenNUMValObject(IDs[2]));
		return TObj;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Max ID): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_TV_CTLPT_ID_RANGE:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
						IRIT_REAL_TO_INT(Prm2 -> U.R));
		int *IDs = TrivIGAGetCtlPtIDRange(IRIT_REAL_TO_INT(Prm1 -> U.R),
					          TV);

		TObj = IPGenLISTObject(IPGenNUMValObject(IDs[0]));
		IPListObjectAppend(TObj, IPGenNUMValObject(IDs[1]));
		return TObj;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Max ID): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_NUM_BZR_ELEMNTS:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        int NumU, NumV, NumW;
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
						IRIT_REAL_TO_INT(Prm2 -> U.R));
		
		ID = TrivIGAGetNumBzrElements(IRIT_REAL_TO_INT(Prm1 -> U.R),
					      TV, &NumU, &NumV, &NumW);
		TObj = IPGenLISTObject(IPGenNUMValObject(NumU));
		IPListObjectAppend(TObj, IPGenNUMValObject(NumV));
		IPListObjectAppend(TObj, IPGenNUMValObject(NumW));
		return TObj;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Number Bezier elements): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_BZR_ELEMNT_CTLPTS:
	    if (IPListObjectLength(ParamList) == 5 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_NUM_OBJ(Prm5)) {
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
						IRIT_REAL_TO_INT(Prm2 -> U.R));
		int Len = TV -> UOrder * TV -> VOrder * TV -> WOrder;  
		TrivIGACtrlPtStruct
		    *CtlPts = TrivIGAGetBzrElementCtrlPts(
						IRIT_REAL_TO_INT(Prm1 -> U.R),
						TV,
						IRIT_REAL_TO_INT(Prm3 -> U.R),
						IRIT_REAL_TO_INT(Prm4 -> U.R),
						IRIT_REAL_TO_INT(Prm5 -> U.R));

		if (CtlPts != NULL) {
		    TObj = IPGenLISTObject(NULL);
		    for (i = 0; i < Len; i++) {
			TObj2 = IPGenCTLPTObject(CtlPts[i].PtType,
						 CtlPts[i].Coord);
			IPListObjectAppend(TObj, TObj2);
		    }
		    IritFree(CtlPts);
		    return TObj;
		}
		else
		    IRIT_NON_FATAL_ERROR("TRIVIGA (Get Bezier Ctlpts): Failed to fetch ctlpts (wrong indices?).");
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Bezier Ctlpts): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_KNOT_INTERVAL:
	    if (IPListObjectLength(ParamList) == 4 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4)) {
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
						IRIT_REAL_TO_INT(Prm2 -> U.R));
		int Len = TRIV_TV_UPT_LST_LEN(TV) *
			  TRIV_TV_VPT_LST_LEN(TV) *
			  TRIV_TV_WPT_LST_LEN(TV);
	        TrivTVDirType
		    Dir = TRIV_INT_TO_DIR(IRIT_REAL_TO_INT(Prm3 -> U.R));
		const CagdRType
		    *Interval = NULL;

		if (Len > 0)
		    Interval = TrivIGAGetKnotInterval(
						IRIT_REAL_TO_INT(Prm1 -> U.R),
						TV, Dir,
						IRIT_REAL_TO_INT(Prm4 -> U.R));
		TObj = IPGenLISTObject(NULL);
		for (i = 0; i < Len; i++) {
		    IPListObjectAppend(TObj, IPGenNUMValObject(Interval[i]));
		}
		return TObj;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Knot Interval): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_MOVE_CTLPT_POS:
	case IRIT_TRIV_IGA_CALL_SET_CTLPT_POS:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_OLST_OBJ(Prm2)) {
	        int Len = IPListObjectLength(Prm2) / 2;
		TrivIGACtrlPtStruct
		    *Vals = (TrivIGACtrlPtStruct *)
				  IritMalloc(Len * sizeof(TrivIGACtrlPtStruct));

		for (i = 0; i < Len; i++) {
		    /* Expects pair's list: PtID1, CtlPt1, PtID2, CtlPt2,... */
		    PTmp = IPListObjectGet(Prm2, i * 2);
		    if (!IP_IS_NUM_OBJ(PTmp))
			break;
		    Vals[i].ID = (int) PTmp -> U.R;

		    PTmp = IPListObjectGet(Prm2, i * 2 + 1);
		    if (!IP_IS_CTLPT_OBJ(PTmp))
			break;

		    Vals[i].PtType  = PTmp -> U.CtlPt.PtType;
		    IRIT_GEN_COPY(Vals[i].Coord, PTmp -> U.CtlPt.Coords,
				  sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
		}

		if (i >= Len) {
		    if (IGAFunc == IRIT_TRIV_IGA_CALL_MOVE_CTLPT_POS)
		        ID = TrivIGAUpdateCtrlPtsPositions(
						IRIT_REAL_TO_INT(Prm1 -> U.R),
						Len, Vals);
		    else
		        ID = TrivIGASetCtrlPtsPositions(
						IRIT_REAL_TO_INT(Prm1 -> U.R),
						Len, Vals);
		    IritFree(Vals);
		}
		else {
		    IritFree(Vals);
		    IRIT_NON_FATAL_ERROR("TRIVIGA (Move/set Ctlpts): Non control point in list.");
		}
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Move/set Ctlpts): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_TV_EVAL:
	    if (IPListObjectLength(ParamList) == 9 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		(Prm6 = IPListObjectGet(ParamList, 5)) != NULL &&
		(Prm7 = IPListObjectGet(ParamList, 6)) != NULL &&
		(Prm8 = IPListObjectGet(ParamList, 7)) != NULL &&
		(Prm9 = IPListObjectGet(ParamList, 8)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_NUM_OBJ(Prm5) &&
		IP_IS_NUM_OBJ(Prm6) &&
		IP_IS_NUM_OBJ(Prm7) &&
		IP_IS_NUM_OBJ(Prm8) &&
		IP_IS_NUM_OBJ(Prm9)) {
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
						IRIT_REAL_TO_INT(Prm2 -> U.R));
	        const TrivIGACtrlPtStruct
		    *CtlPt = TrivIGATVEval(IRIT_REAL_TO_INT(Prm1 -> U.R),
					   TV,
					   IRIT_REAL_TO_INT(Prm3 -> U.R),
					   IRIT_REAL_TO_INT(Prm4 -> U.R),
					   IRIT_REAL_TO_INT(Prm5 -> U.R),
					   IRIT_REAL_TO_INT(Prm6 -> U.R),
					   Prm7 -> U.R,
					   Prm8 -> U.R,
					   Prm9 -> U.R);

		if (CtlPt != NULL) {
		    TObj = IPGenCTLPTObject(CtlPt -> PtType, CtlPt -> Coord);
		    return TObj;
		}
		else
		    IRIT_NON_FATAL_ERROR("TRIVIGA (TV Eval): Failed to evaluate.");
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (TV Eval): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_TV_EVAL_BASIS:
	    if (IPListObjectLength(ParamList) == 6 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		(Prm6 = IPListObjectGet(ParamList, 5)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_NUM_OBJ(Prm5) &&
		IP_IS_NUM_OBJ(Prm6)) {
	        int Len;
	        TrivTVDirType
		    Dir = TRIV_INT_TO_DIR(IRIT_REAL_TO_INT(Prm4 -> U.R));
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
						IRIT_REAL_TO_INT(Prm2 -> U.R));
	        const CagdRType *Basis;

		switch (Dir) {
		    case TRIV_CONST_U_DIR:
			Len = TV -> UOrder;
			break;
		    case TRIV_CONST_V_DIR:
			Len = TV -> VOrder;
			break;
		    case TRIV_CONST_W_DIR:
			Len = TV -> WOrder;
			break;
		    default:
		        Len = 0;
			IRIT_NON_FATAL_ERROR("TRIVIGA (Get Knot Interval): Invalid direction.");
			break;
		}

		if ((Basis = TrivIGATVEvalBasis(IRIT_REAL_TO_INT(Prm1 -> U.R),
						TV,
						IRIT_REAL_TO_INT(Prm3 -> U.R),
						Dir,
						IRIT_REAL_TO_INT(Prm5 -> U.R),
						Prm6 -> U.R)) != NULL) {
		    TObj = IPGenLISTObject(NULL);
		    for (i = 0; i < Len; i++) {
		        IPListObjectAppend(TObj, IPGenNUMValObject(Basis[i]));
		    }
		    return TObj;
		}
		else
		    IRIT_NON_FATAL_ERROR("TRIVIGA (TV Eval Basis): Failed to evaluate.");
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (TV Eval Basis): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_ARNGMNT_FREE:
	    if (IPListObjectLength(ParamList) == 1 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		IP_IS_NUM_OBJ(Prm1)) {
		ID = (int) TrivIGAFreeArrangement(IRIT_REAL_TO_INT(Prm1 -> U.R));
	        Talkative = FALSE;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Free Arrangement): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_ID_BY_TV:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        ID = TrivIGADataManagerGetTrivID((TrivTVStruct *)
							   ((long) Prm2 -> U.R));
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV ID): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_TV_BY_ID:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
					        IRIT_REAL_TO_INT(Prm2 -> U.R));

		if (TV != NULL) {
		    TObj = IPGenTRIVARObject(TrivTVCopy(TV));
		    AttrSetObjectIntAttrib(TObj, "OrigTVPtr", (long) TV);
		    return TObj;
		}
		else
		    IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV ID): TV not found.");
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV ID): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_NEIGHBOR_TVS_IDS:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm3)) {
	        int *IDs,
		    NeighborhoodType = IRIT_REAL_TO_INT(Prm3 -> U.R);
	        TrivIGAAdjacencyInfoStruct AdjInfo[6];
		TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
					       IRIT_REAL_TO_INT(Prm2 -> U.R));

		switch (NeighborhoodType) {
		    case 0:			    /* Facial neighborhood. */
			if (TrivIGAGetFaceNeighboringTVs(
						IRIT_REAL_TO_INT(Prm1 -> U.R), 
					        TV, AdjInfo)) {
			    TObj = IPGenLISTObject(NULL);
			    for (i = 0; i < 6; i++) {
				int Rvrs =
				         (AdjInfo[i].ReverseU ? 0x01 : 0x00) +
				         (AdjInfo[i].ReverseV ? 0x02 : 0x00) +
				         (AdjInfo[i].ReverseUwithV ? 0x04 : 0x00),
				    TVID = AdjInfo[i].AdjTV == NULL ?
					-1 :
					TrivIGADataManagerGetTrivID(
							    AdjInfo[i].AdjTV);

				TObj2 = IPGenLISTObject(
						     IPGenNUMValObject(TVID));
				IPListObjectAppend(TObj2,
						   IPGenNUMValObject(Rvrs));
				IPListObjectAppend(TObj2,
				     IPGenNUMValObject(AdjInfo[i].SameSpace));
				IPListObjectAppend(TObj, TObj2);
			    }
			    return TObj;
			}
			else {
			    IRIT_NON_FATAL_ERROR("TRIVIGA (Get Facial Neighbors): failed to get neighbors.");
			}
			break;
		    case 1:
		    case 2:
		        if (NeighborhoodType == 1)
			    IDs = TrivIGAGetEdgeNeighboringTVs(
						IRIT_REAL_TO_INT(Prm1 -> U.R), 
					        TV);
			else
			    IDs = TrivIGAGetVrtxNeighboringTVs(
						IRIT_REAL_TO_INT(Prm1 -> U.R), 
					        TV);
			if (IDs != NULL) {
			    TObj = IPGenLISTObject(NULL);
			    for (i = 0; IDs[i] != TRIV_IGA_INVALID_TV_ID; i++) {
			        IPListObjectAppend(TObj, IPGenNUMValObject(IDs[i]));
			    }
			    IritFree(IDs);
			    return TObj;
			}
			else {
			    IRIT_NON_FATAL_ERROR("TRIVIGA (Get edge/vrtx Neighbors): failed to get neighbors.");
			}
			break;
		    default:
			IRIT_NON_FATAL_ERROR("TRIVIGA (Get Neighbors): Invalid parameters:\n\texpected 0 (face), 1(edge) or 2(crnr) neighborhood, for third parameter\n");
			break;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get Neighbors): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_ALL_TVS_IDS:
	    if (IPListObjectLength(ParamList) == 1 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		IP_IS_NUM_OBJ(Prm1)) {
		if ((IDs = (int *) TrivIGAGetAllTVs(IRIT_REAL_TO_INT(Prm1 -> U.R))) != NULL) {
		    TObj = IPGenLISTObject(NULL);
		    for (i = 0; IDs[i] != TRIV_IGA_INVALID_TV_ID; i++) {
			IPListObjectAppend(TObj, IPGenNUMValObject(IDs[i]));
		    }
		    IritFree(IDs);
		    return TObj;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get All TVs Indices): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_TV:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        if ((TV = TrivIGAGetTV(
				    IRIT_REAL_TO_INT(Prm1 -> U.R), 
				    IRIT_REAL_TO_INT(Prm2 -> U.R))) != NULL) {
		    TObj = IPGenTRIVARObject(TrivTVCopy(TV));
		    return TObj;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_TV_FACE_CTLPTS_INDICES:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3)) {
	        if ((IDs = (int *) TrivIGAGetTVFaceCtlPtsIDs(
				    IRIT_REAL_TO_INT(Prm1 -> U.R),
				    IRIT_REAL_TO_INT(Prm2 -> U.R),
				    IRIT_REAL_TO_INT(Prm3 -> U.R))) != NULL) {
		    TObj = IPGenLISTObject(NULL);
		    for (i = 0; IDs[i] >= 0; i++) {
		        IPListObjectAppend(TObj, IPGenNUMValObject(IDs[i]));
		    }
		    IritFree(IDs);
		    return TObj;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV Face CtlPts indices): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_TV_FACE_AS_SRF:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3)) {
	        CagdSrfStruct *Srf;

	        if ((Srf = TrivIGAGetTVFaceAsSrf(
				    IRIT_REAL_TO_INT(Prm1 -> U.R),
				    IRIT_REAL_TO_INT(Prm2 -> U.R),
				    IRIT_REAL_TO_INT(Prm3 -> U.R))) != NULL) {
		    TObj = IPGenSRFObject(Srf);
		    return TObj;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV Face Srf): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_TV_CTLPTS_INDICES:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        if ((IDs = (int *) TrivIGAGetTVCtlPtsIndices(
				    IRIT_REAL_TO_INT(Prm1 -> U.R),
				    IRIT_REAL_TO_INT(Prm2 -> U.R))) != NULL) {
		    TObj = IPGenLISTObject(NULL);
		    for (i = 0; IDs[i] >= 0; i++) {
		        IPListObjectAppend(TObj, IPGenNUMValObject(IDs[i]));
		    }
		    IritFree(IDs);
		    return TObj;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV Ctlpts Indices): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_CTLPT:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
	        if ((CtlPt = TrivIGAGetCtlPt(
				    IRIT_REAL_TO_INT(Prm1 -> U.R), 
				    IRIT_REAL_TO_INT(Prm2 -> U.R))) != NULL) {
		    TObj = IPGenCTLPTObject(CtlPt -> PtType, CtlPt -> Coords);
		    return TObj;
		}
	    }
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get CtlPt): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_MATERIAL:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2))
	        ID = (int) TrivIGAGetMaterial(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					      IRIT_REAL_TO_INT(Prm2 -> U.R));
	    else
	        IRIT_NON_FATAL_ERROR("TRIVIGA (Get TV Material ID): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_NEW_MATERIAL:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_STR_OBJ(Prm2))
		ID = (int) TrivIGANewMaterial(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					      Prm2 -> U.Str);
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (New TV Material): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_EXPORT_TO_XML:
	    if (IPListObjectLength(ParamList) == 3 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_STR_OBJ(Prm2) &&
		IP_IS_STR_OBJ(Prm3))
		ID = (int) TrivIGAExportToXML(IRIT_REAL_TO_INT(Prm1 -> U.R),  
					      Prm2 -> U.Str,
					      Prm3 -> U.Str);
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Export to XML): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_LOAD_MATERIAL_FROM_XML:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_STR_OBJ(Prm2))
		ID = (int) TrivIGALoadMaterialFromXML(
					        IRIT_REAL_TO_INT(Prm1 -> U.R), 
						Prm2 -> U.Str);
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Load Material From XML): Invalid parameters.");
	    break;

	case IRIT_TRIV_IGA_CALL_ADD_BNDRY_FACE:
	    if (IPListObjectLength(ParamList) == 6 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		(Prm6 = IPListObjectGet(ParamList, 5)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_STR_OBJ(Prm5) &&
		IP_IS_NUM_OBJ(Prm6)) {
	        TrivTVBndryType Dir;
		TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
					       IRIT_REAL_TO_INT(Prm2 -> U.R));

		switch (IRIT_REAL_TO_INT(Prm3 -> U.R)) {
		    case 0:
			Dir = TRIV_U_MIN_BNDRY;
		        break;
		    case 1:
			Dir = TRIV_U_MAX_BNDRY;
		        break;
		    case 2:
			Dir = TRIV_V_MIN_BNDRY;
		        break;
		    case 3:
			Dir = TRIV_V_MAX_BNDRY;
		        break;
		    case 4:
			Dir = TRIV_W_MIN_BNDRY;
		        break;
		    case 5:
			Dir = TRIV_W_MAX_BNDRY;
		        break;
		    default:
		        Dir = TRIV_NO_BNDRY;
		        IRIT_NON_FATAL_ERROR("TRIVIGA (GAdd BNDRY): Invalid face.");
		        break;
		}

		ID = (int) TrivIGAAddBoundaryFace(
					        IRIT_REAL_TO_INT(Prm1 -> U.R), 
						TV, Dir,
					        IRIT_REAL_TO_INT(Prm4 -> U.R), 
						Prm5 -> U.Str,
						Prm6 -> U.R);
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Add Boundary Face): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_ADD_BNDRY_FACE2:
	    if (IPListObjectLength(ParamList) == 6 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		(Prm6 = IPListObjectGet(ParamList, 5)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_NUM_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_STR_OBJ(Prm5) &&
		IP_IS_NUM_OBJ(Prm6)) {
	        TrivTVBndryType Dir;
		TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
					       IRIT_REAL_TO_INT(Prm2 -> U.R));

		switch (IRIT_REAL_TO_INT(Prm3 -> U.R)) {
		    case 0:
			Dir = TRIV_U_MIN_BNDRY;
		        break;
		    case 1:
			Dir = TRIV_U_MAX_BNDRY;
		        break;
		    case 2:
			Dir = TRIV_V_MIN_BNDRY;
		        break;
		    case 3:
			Dir = TRIV_V_MAX_BNDRY;
		        break;
		    case 4:
			Dir = TRIV_W_MIN_BNDRY;
		        break;
		    case 5:
			Dir = TRIV_W_MAX_BNDRY;
		        break;
		    default:
		        Dir = TRIV_NO_BNDRY;
		        IRIT_NON_FATAL_ERROR("TRIVIGA (GAdd BNDRY): Invalid face.");
		        break;
		}

		ID = (int) TrivIGAAddBoundaryFace2(
					        IRIT_REAL_TO_INT(Prm1 -> U.R), 
					        TV, Dir, 
					        IRIT_REAL_TO_INT(Prm4 -> U.R), 
						Prm5 -> U.Str,
						Prm6 -> U.R);
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Add Boundary Face 2): Invalid parameters.");
	    break;
        case IRIT_TRIV_IGA_CALL_ADD_BNDRY_FACE_BY_PT:
	    if (IPListObjectLength(ParamList) == 6 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		(Prm3 = IPListObjectGet(ParamList, 2)) != NULL &&
		(Prm4 = IPListObjectGet(ParamList, 3)) != NULL &&
		(Prm5 = IPListObjectGet(ParamList, 4)) != NULL &&
		(Prm6 = IPListObjectGet(ParamList, 5)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2) &&
		IP_IS_POINT_OBJ(Prm3) &&
		IP_IS_NUM_OBJ(Prm4) &&
		IP_IS_STR_OBJ(Prm5) &&
		IP_IS_NUM_OBJ(Prm6)) {
		TrivTVStruct
		    *TV = TrivIGADataManagerGetTrivariate(
					       IRIT_REAL_TO_INT(Prm2 -> U.R));

		ID = (int) TrivIGAAddBoundaryFaceByPt(
					        IRIT_REAL_TO_INT(Prm1 -> U.R), 
						TV, Prm3 -> U.Pt,
					        IRIT_REAL_TO_INT(Prm4 -> U.R), 
						Prm5 -> U.Str,
						Prm6 -> U.R);
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Add Boundary Face): Invalid parameters.");
	    break;
        case IRIT_TRIV_IGA_CALL_GET_BNDRY_FACE_ID_BY_PT:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_POINT_OBJ(Prm2)) {
		IDs = TrivIGAGetBoundaryFaceByPt(IRIT_REAL_TO_INT(Prm1 -> U.R), 
						 NULL, Prm2 -> U.Pt);
		TObj = IPGenLISTObject(NULL);
		IPListObjectAppend(TObj, IPGenNUMValObject(IDs[0]));
		IPListObjectAppend(TObj, IPGenNUMValObject(IDs[1]));
		return TObj;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Boundary Face ID by Pt): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_GET_LAST_ERROR:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2))
		ID = (int) TrivIGAGetLastError(IRIT_REAL_TO_INT(Prm1 -> U.R), 
					       IRIT_REAL_TO_INT(Prm2 -> U.R));
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Get Last Error): Invalid parameters.");
	    break;
	case IRIT_TRIV_IGA_CALL_DESCRIBE_ERROR:
	    if (IPListObjectLength(ParamList) == 2 &&
		(Prm1 = IPListObjectGet(ParamList, 0)) != NULL &&
		(Prm2 = IPListObjectGet(ParamList, 1)) != NULL &&
		IP_IS_NUM_OBJ(Prm1) &&
		IP_IS_NUM_OBJ(Prm2)) {
		const char
		    *ErrStr = TrivIGADescribeError(IRIT_REAL_TO_INT(Prm2 -> U.R));

		TObj = IPGenSTRObject(ErrStr);
		return TObj;
	    }
	    else
		IRIT_NON_FATAL_ERROR("TRIVIGA (Describe Error): Invalid parameters.");
	    break;

        default:
        case IRIT_TRIV_IGA_CALL_ERROR:
	    IRIT_NON_FATAL_ERROR("TRIVIGA: Invalid function call.");
	    break;
    }

    return IPGenNUMValObject(ID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs geometry of the given text in teh desired font, font type,    M
* size and spacing.                                                          M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:            The text to layout.  A string.			     M
*   FontName:        The font name to use.  I.e. "Times new Roman".          M
*   FontStyle:       Font style to use.  I.e. italics.                       M
*   FontSpaceWidth:  Space added to individual chars.	                     M
*   Text3DEdgeType:  For 3D text geometry controls edges (i.e. chamferred.). M
*   Text3DSetup:     For 3D text, a vector of (Chamfer offset size, 3D       M
*                    extruded vertical distance).                            M
*   Tolerance:       Of approximating freeforms using polygons.              M
*   OutputType:      Selects the type of geometry to create:                 M
*                    0. Outline Bezier curves as in the font data.           M
*                    1. Outline B-spline curves (merging Bezier curves in 1).M
*                    2. Solid 2D polygons, for the outline geometry          M
*                       (polygons).				             M
*                    3. Both Solid 2D and B-spline outline (2+3 above).      M
*                    4. Full 3D text (polygons).			     M
*                    5. Solid 2D polygons, for the outline geometry          M
*                       ((trimmed) surfaces).			             M
*                    6. Full 3D text ((trimmed) surfaces).		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   The constructed geometry or NULl if failed.          M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontConvertFontToBezier, FTStringOutline2BezierCurves                M
*   LayoutTextOverShape							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   ConvertText2Geom                                                         M
*****************************************************************************/
IPObjectStruct *ConvertText2Geom(const char *Text,
				 const UserFontName FontName,
				 IrtRType *FontStyle,
				 IrtRType *FontSpaceWidth,
				 IrtRType *Text3DEdgeType,
				 IPObjectStruct *Text3DSetup,
				 IrtRType *Tolerance,
				 IrtRType *OutputType)
{
    char *ErrStr;
    IrtRType T3DSetup[2];
    UserFontStyleType FStyle;
    UserFontGeomOutputType OType;
    UserFont3DEdgeType EType;
    IPObjectStruct *TextGeom;

    switch (IRIT_REAL_PTR_TO_INT(FontStyle)) {
	case 0:
	default:
	    FStyle = USER_FONT_STYLE_REGULAR;
	    break;
	case 1:
	    FStyle = USER_FONT_STYLE_ITALICS;
	    break;
	case 2:
	    FStyle = USER_FONT_STYLE_BOLD;
	    break;
	case 3:
	    FStyle = USER_FONT_STYLE_BOLD_ITALICS;
	    break;
    }

    switch (IRIT_REAL_PTR_TO_INT(OutputType)) {
        default:
        case 0:
	    OType = USER_FONT_OUTPUT_BEZIER_CRVS;
	    break;
        case 1:
	    OType = USER_FONT_OUTPUT_BSPLINE_CRVS;
	    break;
        case 2:
  	    OType = USER_FONT_OUTPUT_FILLED2D_POLYS;
	    break;
        case 3:
  	    OType = USER_FONT_OUTPUT_OUTLINE_FILLED2D_POLYS;
	    break;
        case 4:
	    OType = USER_FONT_OUTPUT_SOLID3D_POLYS;
	    break;
	case 5:
	    OType = USER_FONT_OUTPUT_FILLED2D_TSRFS;
	    break;
	case 6:
	    OType = USER_FONT_OUTPUT_SOLID3D_TSRFS;
	    break;
    }

    switch (IRIT_REAL_PTR_TO_INT(Text3DEdgeType)) {
        default:
        case 0:
	    EType = USER_FONT_3D_EDGE_NORMAL;
	    break;
        case 1:
	    EType = USER_FONT_3D_EDGE_CHAMFER;
	    break;
        case 2:
  	    EType = USER_FONT_3D_EDGE_ROUND;
    }

    T3DSetup[0] = 0.01;
    T3DSetup[1] = 0.1;
    if (IP_IS_OLST_OBJ(Text3DSetup)) {
        const IPObjectStruct
	    *Num1 = IPListObjectGet(Text3DSetup, 0),
	    *Num2 = IPListObjectGet(Text3DSetup, 1);

	if (IP_IS_NUM_OBJ(Num1) && IP_IS_NUM_OBJ(Num2)) {
	    T3DSetup[0] = Num1 -> U.R;
	    T3DSetup[1] = Num2 -> U.R;
	}
	else
	    IRIT_NON_FATAL_ERROR("Txt2Geom: Expected a numeric item in list for 3D setup.\n");
    }
    else
        IRIT_NON_FATAL_ERROR("Txt2Geom: Expected a list object for 3D setup.\n");

    /* Convert the text to geometry: */

#   ifdef __WINNT__
        if (!UserFontConvertTextToGeom(UserAscii2WChar(Text), FontName,
				       FStyle, 1.0, *FontSpaceWidth, EType,
				       T3DSetup, *Tolerance, OType, FALSE,
				       NULL, &TextGeom, &ErrStr)) {


#   else
	if (!UserFontConvertTextToGeom((const UserFontText) Text, FontName,
				       FStyle, 1.0, *FontSpaceWidth, EType,
				       T3DSetup, *Tolerance, OType, FALSE,
				       NULL, &TextGeom, &ErrStr)) {
#   endif /* __WINNT__ */
	    fprintf(stderr, "Error in Freefont conversion: %s\n", ErrStr);
	    return NULL;
	}

    return TextGeom;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Layout the given text over the given bounding region.  Constructed Text  M
* will be controlled by FontName, FontStyle, FontSize, and FontSpaceWidth    M
* and will be aligned to follow TextAlignment.				     M
*   All construction is conducted over the XY plane, in 2D.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:            The text to layout.  A string.			     M
*   FontName:        The font name to use.  I.e. "Times new Roman".          M
*   FontStyle:       Font style to use.  I.e. italics.			     M
*   FontSize:        The size of the constructed text.                       M
*   FontSpace:       (WordWidth, SpaceWidth, LineHeight) spacing, specified  M
*   Tolerance:	     For 2D filled polygons and 3D solid text geometry.      M
*   Text3DEdgeType:  For 3D text geometry controls edges (i.e. chamferred.). M
*   Text3DSetup:     For 3D text, a vector of (Chamfer offset size, 3D       M
*                    extruded vertical distance).                            M
*                      In text font's point units.			     M
*   AlignmentType:   Text alignment, left, centered, etc.                    M
*   OutputType:      Output should be original Bezier curves, merged         M
*                    B-spline curves, filled planar polys, 3D polys, etc.    M
*   BoundingShape:   A closed region (a curve or a poly) to place text in.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectSTruct *:   Layed out geometry of text.		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontConvertFontToBezier, FTStringOutline2BezierCurves                M
*   UserFontLayoutOverShape, ConvertText2Bezier                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   LayoutTextOverShape                                                      M
*****************************************************************************/
IPObjectStruct *LayoutTextOverShape(const char *Text,
				    const UserFontName FontName,
				    IrtRType *FontStyle,
				    IrtRType *FontSize,
				    const IPObjectStruct *FontSpace,
				    IrtRType *Tolerance,
				    IrtRType *Text3DEdgeType,
				    const IPObjectStruct *Text3DSetup,
				    IrtRType *AlignmentType,
				    IrtRType *OutputType,
				    const IPObjectStruct *BoundingShape)
{
    char *ErrStr = "";
    int RetVal = FALSE;
    IrtRType T3DSetup[2], FSpace[3];
    IPObjectStruct
	*LayoutGeom = NULL;
    UserFontStyleType FStyle;
    UserFontAlignmentType AlignType;
    UserFontGeomOutputType GOType;
    UserFont3DEdgeType EType;

    switch (IRIT_REAL_PTR_TO_INT(FontStyle)) {
	default:
	case 0:
	    FStyle = USER_FONT_STYLE_REGULAR;
	    break;
	case 1:
	    FStyle = USER_FONT_STYLE_ITALICS;
	    break;
	case 2:
	    FStyle = USER_FONT_STYLE_ITALICS;
	    break;
	case 3:
	    FStyle = USER_FONT_STYLE_BOLD_ITALICS;
	    break;
    }

    switch (IRIT_REAL_PTR_TO_INT(AlignmentType)) {
	default:
	case 0:
	    AlignType = USER_FONT_ALIGN_LEFT;
	    break;
	case 1:
	    AlignType = USER_FONT_ALIGN_CENTER;
	    break;
	case 2:
	    AlignType = USER_FONT_ALIGN_RIGHT;
	    break;
	case 3:
	    AlignType = USER_FONT_ALIGN_WIDE;
	    break;
    }

    switch (IRIT_REAL_PTR_TO_INT(OutputType)) {
	default:
	case 0:
	    GOType = USER_FONT_OUTPUT_BEZIER_CRVS;
	    break;
	case 1:
	    GOType = USER_FONT_OUTPUT_BSPLINE_CRVS;
	    break;
	case 2:
	    GOType = USER_FONT_OUTPUT_FILLED2D_POLYS;
	    break;
	case 3:
	    GOType = USER_FONT_OUTPUT_SOLID3D_POLYS;
	    break;
	case 4:
	    GOType = USER_FONT_OUTPUT_FILLED2D_TSRFS;
	    break;
	case 5:
	    GOType = USER_FONT_OUTPUT_SOLID3D_TSRFS;
	    break;
    }

    switch (IRIT_REAL_PTR_TO_INT(Text3DEdgeType)) {
        default:
        case 0:
	    EType = USER_FONT_3D_EDGE_NORMAL;
	    break;
        case 1:
	    EType = USER_FONT_3D_EDGE_CHAMFER;
	    break;
        case 2:
  	    EType = USER_FONT_3D_EDGE_ROUND;
    }

    if (IP_IS_OLST_OBJ(FontSpace)) {
        IPObjectStruct
	    *NumerObj1 = IPListObjectGet(FontSpace, 0),
	    *NumerObj2 = IPListObjectGet(FontSpace, 1),
	    *NumerObj3 = IPListObjectGet(FontSpace, 2);

	if (IP_IS_NUM_OBJ(NumerObj1) &&
	    IP_IS_NUM_OBJ(NumerObj2) &&
	    IP_IS_NUM_OBJ(NumerObj3)) {
	    FSpace[0] = NumerObj1 -> U.R;
	    FSpace[1] = NumerObj2 -> U.R;
	    FSpace[2] = NumerObj3 -> U.R;
	}
	else {
	    fprintf(stderr, "Expected only a poly or a curve as a bounding shape.\n");
	    return NULL;
	}
    }

    if (IP_IS_OLST_OBJ(Text3DSetup)) {
        IPObjectStruct
	    *NumerObj1 = IPListObjectGet(Text3DSetup, 0),
	    *NumerObj2 = IPListObjectGet(Text3DSetup, 1);

	if (IP_IS_NUM_OBJ(NumerObj1) && IP_IS_NUM_OBJ(NumerObj2)) {
	    T3DSetup[0] = NumerObj1 -> U.R;
	    T3DSetup[1] = NumerObj2 -> U.R;
	}
	else {
	    fprintf(stderr, "TextLayShp: tExpected numeric values only in list.\n");
	    return NULL;
	}
    }

#   ifdef __WINNT__
        if (IP_IS_POLY_OBJ(BoundingShape))
	    RetVal = UserFontLayoutOverShape(
					UserAscii2WChar(Text), FontName,
					FStyle,	*FontSize, FSpace,
					*Tolerance, EType, T3DSetup, AlignType,
					BoundingShape -> U.Pl,
					GOType, &LayoutGeom, &ErrStr);
	else if (IP_IS_CRV_OBJ(BoundingShape))
	    RetVal = UserFontLayoutOverShape2(
					UserAscii2WChar(Text), FontName,
					FStyle,	*FontSize, FSpace,
					*Tolerance, EType, T3DSetup, AlignType,
					BoundingShape -> U.Crvs,
					GOType, &LayoutGeom, &ErrStr);
	else
	    fprintf(stderr, "TextLayShp: Expected only a poly or a curve as a bounding shape.\n");
#   else
        if (IP_IS_POLY_OBJ(BoundingShape))
	    RetVal = UserFontLayoutOverShape(
					(UserFontText) Text,
					(UserFontName) FontName,
					FStyle,	*FontSize, FSpace,
					*Tolerance, EType, T3DSetup, AlignType,
					BoundingShape -> U.Pl,
					GOType, &LayoutGeom, &ErrStr);
	else if (IP_IS_CRV_OBJ(BoundingShape))
	    RetVal = UserFontLayoutOverShape2(
					(UserFontText) Text,
					(UserFontName) FontName,
					FStyle,	*FontSize, FSpace,
					*Tolerance, EType, T3DSetup, AlignType,
					BoundingShape -> U.Crvs,
					GOType, &LayoutGeom, &ErrStr);
	else
	    fprintf(stderr, "TextLayOut: Expected only a poly or a curve as a bounding shape.\n");
#   endif /* __WINNT__ */

    if (!RetVal) {
        fprintf(stderr, "TextLayOut: Failed to place the desired text in the designated bounding region:\n%s\n", ErrStr);
    }

    return LayoutGeom;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes Srf(Crv), by evaluating Srf at all locations of Crv.            *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv:    In UV space of Srf to map to E3.                                 *
*   Srf:    To map Crv over it.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *:    E3 curve of Srf(Crv).                                *
*****************************************************************************/
static CagdCrvStruct *IritMapUVCrv2E3(const CagdCrvStruct *Crv,
				      const CagdSrfStruct *Srf)
{
    int i;
    const CagdRType * const 
	*Pts = (const CagdRType * const *) Crv -> Points;
    CagdCrvStruct *E3Crv;

    if (Crv -> Length == 2) {
        /* Assume it is the VMin initial boundary. */
        assert(Crv -> Points[2] == NULL ||
	       Crv -> Points[2][0] == Crv -> Points[2][1]); /* Const V crv. */
        return CagdCrvFromMesh(Srf, 0, CAGD_CONST_V_DIR);
    }

    E3Crv = BspCrvNew(Crv -> Length, Crv -> Order, CAGD_PT_E3_TYPE);
    BspKnotUniformOpen(E3Crv -> Length, E3Crv -> Order, E3Crv -> KnotVector);
		  
    for (i = 0; i < Crv -> Length; i++) {
        CagdRType
	    *R = CagdSrfEval(Srf, Pts[1][i], Pts[2] == NULL ? 0.0 : Pts[2][i]);
	CagdPType Pt;

	CagdCoerceToE3(Pt, &R, -1, Srf -> PType);
	E3Crv -> Points[1][i] = Pt[0];
	E3Crv -> Points[2][i] = Pt[1];
	E3Crv -> Points[3][i] = Pt[2];
    }

    return E3Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A wrapper routine to handle flank milling analysis of line surface       M
* contact tool path planning.                                                M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:           To plan flank milling tool path for.			     M
*   RTolerance:    Of ruled surfaces approximation of Srf using strips of    M
*                  flank milling.					     M
*   REuclidean:    TRUE to map the resulting curves to E3 space, over Srf.   M
*   CrvSizeReduction: A reduction in size of traced curve while ensuring     M
*                  the Tolerance conservatively.			     M
*   RSubdivTol:    Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   RNumericTol:   Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct: The planned tool path 9and strip boundaries.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarFlankMillLineAnalyze                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritFlankMillLineAnalyze                                                 M
*****************************************************************************/
IPObjectStruct *IritFlankMillLineAnalyze(const IPObjectStruct *Srf,
					 IrtRType *RTolerance,
					 IrtRType *REuclidean,
					 IrtRType *RCrvSizeReduction,
					 IrtRType *RSubdivTol,
					 IrtRType *RNumericTol)
{
    CagdCrvStruct *ToolPath, *StripBndries;

    if ((ToolPath = MvarFlankMillLineAnalyze(
				     Srf -> U.Srfs, *RTolerance, &StripBndries, 
				     IRIT_REAL_PTR_TO_INT(RCrvSizeReduction),
				     *RSubdivTol, *RNumericTol)) != NULL) {
        IPObjectStruct *PObj, *PCrv, *PList;
	CagdCrvStruct *Crv;

	/* Handle the tool path curves. */
	PList = IPGenLISTObject(NULL);
	for (Crv = ToolPath; Crv != NULL; Crv = Crv -> Pnext) {
	    if (IRIT_REAL_PTR_TO_INT(REuclidean))
	        PCrv = IPGenCRVObject(IritMapUVCrv2E3(Crv, Srf -> U.Srfs));
	    else
	        PCrv = IPGenCRVObject(CagdCrvCopy(Crv));
	    IPListObjectAppend(PList, PCrv);
	}
	PObj = IPGenLISTObject(PList);

	/* Handle the strip boundary curves curves. */
	PList = IPGenLISTObject(NULL);
	for (Crv = StripBndries; Crv != NULL; Crv = Crv -> Pnext) {
	    if (IRIT_REAL_PTR_TO_INT(REuclidean))
	        PCrv = IPGenCRVObject(IritMapUVCrv2E3(Crv, Srf -> U.Srfs));
	    else
	        PCrv = IPGenCRVObject(CagdCrvCopy(Crv));
	    IPListObjectAppend(PList, PCrv);
	}
	IPListObjectAppend(PObj, PList);

	CagdCrvFreeList(ToolPath);
	CagdCrvFreeList(StripBndries);

        return PObj;
    }
    else
	return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A wrapper routine to compute 2-contact motion analysis of two planar     M
* closed C^1 continuous curves.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Crv1, Crv2:    Two curves to compute 2-contact motion analysis for, when M
*                  one curve is rolling against the other curve.             M
*   RStepSize:     Tolerace of step size for numerical tracing.              M
*   RSubdivTol:    Tolerance of the subdivision process.  Tolerance is       M
*		   measured in the parametric space of the multivariates.    M
*   RNumericTol:   Numeric tolerance of the numeric stage.  The numeric      M
*		   stage is employed only if NumericTol < SubdivTol.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct: The 2-contact motion curves, of 5-tuples, prescribing    M
*                  the motion as (RotAngle, u1, v1, u2, v2), where u1/2 are  M
*                  the contact (parametric) locations on Crv1 and v1/2 are   M
*                  the contact (parametric) location of Crv2 and RotAngle is M
*                  the rotation 9in Z) needed of Crv2 to get to this         M
*                  2-contact state.				             M
*                                                                            *
* SEE ALSO:                                                                  M
*   Mvar2CntctCompute2CntctMotion                                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritCurves2ContactAnal                                                   M
*****************************************************************************/
IPObjectStruct *IritCurves2ContactAnal(const IPObjectStruct *Crv1,
				       const IPObjectStruct *Crv2,
				       IrtRType *RStpSize,
				       IrtRType *RSubdivTol,
				       IrtRType *RNumericTol)
{
    MvarPolylineStruct
        *MVPlls = Mvar2CntctCompute2CntctMotion(Crv1 -> U.Crvs, Crv2 -> U.Crvs,
						*RStpSize, 
						*RSubdivTol, *RNumericTol);
    IPObjectStruct *PObj;

    if (MVPlls == NULL)
        return NULL;

    PObj = MvarCnvrtMVPolysToIritPolys(MVPlls);
    MvarPolylineFreeList(MVPlls);
    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tile an input object, in place, (UTimes x VTimes (x WTimes)) in given    M
* bi/trivariate.  Computation is:					     M
* + approximated for trivariate DeformObj unless Precise is set and Tile is  M
*   a freeform curve or surface.					     M
* + precise for bivariate DeformObj, when Tile can only be curve(s) or       M
*   surface(s).								     M
*                                                                            *
* PARAMETERS:                                                                M
*   Tile:       The object to map through the trivariate, in place.          M
*   DeformObj:  The mapping/deformation function from R2/3 to R2/3.	     M
*   UTimes, VTimes, WTimes:  Number of times to tile the object in each      M
*               axis. WTimes is ignored for a surface DeformObj.             M
*   FitObj:     TRUE to rescale PObj tile to precisely fit the domain        M
*		                              (UTimes x VTimes (x WTimes)),  M
*               FALSE to assume PObj is in [0,1]^2/3 when fitting domain.    M
*   Precise:    IF TRUE and possible, compute the FFD precsiely, using       M
*               composition.  Ignored, for a bivariate DeformObj as input    M
*               Tile is either curves or surfaces only.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  (UTimes x VTimes x WTimes) mapped and deformed        M
*               objects.					             M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivFFDObjectTV, TrivFFDTileObjectInTV, SymbComposeTileObjectInSrf       M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritObjectTileFFD	                                                     M
*****************************************************************************/
IPObjectStruct *IritObjectTileFFD(const IPObjectStruct *Tile,
			          const IPObjectStruct *DeformObj,
			          IrtRType *UTimes,
			          IrtRType *VTimes,
			          IrtRType *WTimes,
			          IrtRType *FitObj,
			          IrtRType *Precise)
{
    IPObjectStruct
        *PRetVal = NULL;

    if (IP_IS_TRIVAR_OBJ(DeformObj)) {
        if (IRIT_REAL_PTR_TO_INT(Precise))
	    PRetVal = TrivComposeTileObjectInTV(Tile, DeformObj -> U.Trivars,
						*UTimes, *VTimes, *WTimes,
						IRIT_REAL_PTR_TO_INT(FitObj));
	else
	    PRetVal = TrivFFDTileObjectInTV(Tile, DeformObj -> U.Trivars,
					    *UTimes, *VTimes, *WTimes,
					    IRIT_REAL_PTR_TO_INT(FitObj));
    }
    else if (IP_IS_SRF_OBJ(DeformObj)) {
        if (!IP_IS_CRV_OBJ(Tile) && !IP_IS_SRF_OBJ(Tile)) {
	    fprintf(stderr, "Input tile must be either curve(s) or surface(s), for a surface defromation\n");
	    return NULL;
	}
        PRetVal = SymbComposeTileObjectInSrf(Tile, DeformObj -> U.Srfs,
					     *UTimes, *VTimes,
					     IRIT_REAL_PTR_TO_INT(FitObj));
    }

    return PRetVal;
}
/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given one or two meshes, PolyMeshes, and possibly (sharing boundary)     M
* edge(s), CommonEdges, construct an updated mesh that is rounded at the     M
* desired locations.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   RoundingMethod: 1 for global low pass filtering of one input mesh.       M
*                   2 for global low pass filtering of one input meshes.     M
*                   3 for rounding along shared CommonEdges between two      M
*                     input meshes.					     M
*   PolyMeshes:  A list object of 2 or 3 meshes sharing edges (and possibly  M
*                one vertex).						     M
*   CommonEdges: A list object of the shared edge(s):  ignored if            M
*                RoundingMethod = 1, a list of edges of any length if        M
*                RoundingMethod = 2 and one edge if two meshes are given in  M
*                PolyMeshes for RoundingMethod = 3.			     M
*   Params:      List of numeric parameters to control the rounding:         M
*                For RoundingMethod equal				     M
*                1) list(SmoothNormals, NumIterations, RoundPower,           M
*                        AllowBndryMove, RoundingRadius)                     M
*                   where SmoothNormals will compute smooth vertices normals M
*                   as a post-process, where NumIterations is a positive     M
*                   integer controlling how many iterations to apply,        M
*                   and RoundPower is a number between zero and one          M
*                   controlling smoothing affect (One for maximal effect).   M
*                     AllowBndryMove TRUE to allow the boundary vertices to  M
*                   move, FALSE to keep them fixed.			     M
*                     RoundingRadius affects the radius of influence from    M
*                   the designated (restricted/allowed to move) vertices.    M
*                3) list(SmoothNormals, RoundRadius, RoundPower)             M
*                   where SmoothNormals will compute smooth vertices normals M
*                   as a post-process, where RoundRadius is the desired      M
*                   radius of the approximated blend (precise only for       M
*                   orthogonal two meshes meeting at the common edge) and    M
*                   RoundPower is some Bias to affect the rounding size,     M
*                   with 1.0 to have no affect, and values larger (smaller)  M
*                   than 1.0 to enlarge (shrink) the rounding size.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectstruct *: Rounded meshes along shared edges(s) or NULL if error. M
*                                                                            *
* SEE ALSO:                                                                  M
*   User2PolyMeshRoundEdge, GMPolyMeshSmoothing, GMBlendNormalsToVertices    M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritPolyMeshRound                                                        M
*****************************************************************************/
IPObjectStruct *IritPolyMeshRound(IrtRType *RoundingMethods,
				  const IPObjectStruct *PolyMeshes,
			          const IPObjectStruct *CommonEdges,
			          const IPObjectStruct *Params)
{
    int SmoothNormals, OldCopyRef;
    IrtRType Param1, Param2, AllowBndryMove, RoundingRadius;
    IPObjectStruct *PMeshes, *Pl1, *Pl2;
    const IPObjectStruct *PTmp;

    PTmp = IPListObjectGet(Params, 0);
    if (PTmp != NULL && IP_IS_NUM_OBJ(PTmp))
        SmoothNormals = (int) PTmp -> U.R;
    else {
        fprintf(stderr, "Parameters' lists is wrong.\n");
	return NULL;
    }
    PTmp = IPListObjectGet(Params, 1);
    if (PTmp != NULL && IP_IS_NUM_OBJ(PTmp))
	Param1 = PTmp -> U.R;
    else {
        fprintf(stderr, "Parameters' lists is wrong.\n");
	return NULL;
    }
    PTmp = IPListObjectGet(Params, 2);
    if (PTmp != NULL && IP_IS_NUM_OBJ(PTmp))
	Param2 = PTmp -> U.R;
    else {
        fprintf(stderr, "Parameters' lists is wrong.\n");
	return NULL;
    }

    if (CommonEdges != NULL && IP_IS_OLST_OBJ(CommonEdges))
	CommonEdges = IPListObjectGet(CommonEdges, 0);

    switch (IRIT_REAL_PTR_TO_INT(RoundingMethods)) {
	default:
	    fprintf(stderr, "Round method type illegal.\n");
	    return NULL;
        case 1:
	    if (!IP_IS_POLY_OBJ(PolyMeshes)) {
	        fprintf(stderr, "Input does not contain a polygonal mesh.\n");
	        return NULL;
	    }

	    OldCopyRef = IPSetCopyObjectReferenceCount(FALSE);
	    Pl1 = IPCopyObject(NULL, PolyMeshes, FALSE);
	    IPSetCopyObjectReferenceCount(OldCopyRef);

	    PTmp = IPListObjectGet(Params, 3);
	    if (PTmp != NULL && IP_IS_NUM_OBJ(PTmp))
		AllowBndryMove = PTmp -> U.R;
	    else {
		fprintf(stderr, "Parameters' lists is wrong.\n");
		IPFreeObject(Pl1);
		return NULL;
	    }
	    PTmp = IPListObjectGet(Params, 4);
	    if (PTmp != NULL && IP_IS_NUM_OBJ(PTmp))
		RoundingRadius = PTmp -> U.R;
	    else {
	        fprintf(stderr, "Parameters' lists is wrong.\n");
	        IPFreeObject(Pl1);
		return NULL;
	    }

	    GMPolyMeshSmoothing(Pl1,
				CommonEdges == NULL ? NULL 
						    : CommonEdges -> U.Pl,
				(int) AllowBndryMove, RoundingRadius, 
				(int) Param1, Param2);
	    GMBlendNormalsToVertices(Pl1 -> U.Pl, 180);
	    return Pl1;
        case 2:
	    return NULL;
        case 3:
	    if (!IP_IS_OLST_OBJ(PolyMeshes)) {
	        fprintf(stderr, "Input does not contain 2 polygonal meshes in a list.\n");
	        return NULL;
	    }

	    OldCopyRef = IPSetCopyObjectReferenceCount(FALSE);
	    PMeshes = IPCopyObject(NULL, PolyMeshes, FALSE);
	    IPSetCopyObjectReferenceCount(OldCopyRef);
	    Pl1 = IPListObjectGet(PMeshes, 0);
	    Pl2 = IPListObjectGet(PMeshes, 1);

	    if (Pl1 == NULL ||
		!IP_IS_POLY_OBJ(Pl1) ||
		Pl2 == NULL ||
		!IP_IS_POLY_OBJ(Pl2)) {
	        fprintf(stderr, "Input does not contain 2 polygonal meshes.\n");
	        IPFreeObject(PMeshes);
	        return NULL;
	    }

	    if (CommonEdges == NULL ||
		!IP_IS_POLY_OBJ(CommonEdges) ||
		!IP_IS_POLYLINE_OBJ(CommonEdges)) {
	        fprintf(stderr, "Input does not contain a polyline edge.\n");
		IPFreeObject(PMeshes);
		return NULL;
	    }

	    if (!User2PolyMeshRoundEdge(Pl1 -> U.Pl, Pl2 -> U.Pl,
					CommonEdges -> U.Pl,
					Param1, Param2)) {
	        IPFreeObject(PMeshes);
		return NULL;	    
	    }
	    else {
	        GMBlendNormalsToVertices(Pl1 -> U.Pl, 180);
		GMBlendNormalsToVertices(Pl2 -> U.Pl, 180);
		return PMeshes;
	    }
    }
}
