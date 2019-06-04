/******************************************************************************
* Font3D.cpp - Converts outline font curves to 3D using solid extrusion.      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Roee Gavriel and Gershon Elber, May 2009.			      *
******************************************************************************/

#include <ctype.h>

#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/geom_lib.h"
#include "user_loc.h"

#define USER_FONT_3D_BBOX_EXPANSION	0.05

#define USER_FNT2BZR_IGNORE_V_TAG  0x10
#define USER_FNT2BZR_IS_IGNORED_VRTX(Vrtx) \
                                   (Vrtx -> Tags & USER_FNT2BZR_IGNORE_V_TAG)
#define USER_FNT2BZR_SET_IGNORED_VRTX(Vrtx) \
				  (Vrtx -> Tags |= USER_FNT2BZR_IGNORE_V_TAG)
#define USER_FNT2BZR_RST_IGNORED_VRTX(Vrtx) \
				 (Vrtx -> Tags &= ~USER_FNT2BZR_IGNORE_V_TAG)
#define USER_FNT2BZR_RST_POLY_IGNORED_VRTX(Pl) { \
    IPVertexStruct *V = Pl -> PVertex; \
    do  { \
        USER_FNT2BZR_RST_IGNORED_VRTX(V); \
	V = V -> Pnext; \
    } while (V != NULL && V != Pl -> PVertex); \
}

static void UserFontCrvOffset(CagdCrvStruct *Crv, CagdRType Offset);
static TrimSrfStruct *UserFontCreateTrimmedSrfFromCurves(CagdCrvStruct
							         *BspCrvList);

static IPObjectStruct *UserFontCnvrtPllns2Plgns(IPPolygonStruct *Pllns,
						int HaveHoles);
static IPPolygonStruct *UserCnvrtPllns2PlgnsAux(const TrimCrvStruct *Crv,
						const TrimCrvStruct
						                 *NestedCrvs);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   COnverts Text to geometry.  Note this function will behave differently   M
* on different platforms (i.e. Windows vs. Unix.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:           The string text to convert to geometry.                  M
*   FontName:       The font to use for conversion.                          M
*   FontStyle:      The font style (italic, bold, etc.)                      M
*   FontSize:       Created text scaling relative factor.                    M
*   TextSpace:      In-word relative spacing.                                M
*   Text3DEdgeType: For 3D text geometry controls edges (i.e. chamferred.).  M
*   Text3DSetup:    For 3D text, a vector of (Chamfer offset size, 3D        M
*                   extruded vertical distance).                             M
*   Tolerance:	    For 2D filled polygons and 3D solid text geometry.       M
*   OutputType:     Selects the type of geometry to create:                  M
*                   0. Outline Bezier curves as in the font data.            M
*                   1. Outline B-spline curves (merging Bezier curves in 1). M
*                   2. Solid 2D polygons, for the outline geometry           M
*                      (polygons).				             M
*                   3. Both Solid 2D and B-spline outline (2+3 above).       M
*                   4. Full 3D text (polygons).				     M
*                   5. Solid 2D polygons, for the outline geometry           M
*                      ((trimmed) surfaces).			             M
*                   6. Full 3D text ((trimmed) surfaces).		     M
*   CompactOutput:  If TRUE, merged trimming surfaces or polygons into       M
*                   single objects.					     M
*   PlacedTextBaseName:  Base name to use in placed text.  Can be NULL.      M
*   PlacedTextGeom: Output created geometry will be returned here.           M
*   ErrorStr:       Will be set with the error message if was one.           M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:       TRUE if successful, FALSE if conversion failed.               M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontConvertFontToBezier, UserFontFTStringOutline2BezierCurves        M
*   UserFontBspCrv2Poly,  UserFontBspList2Solids                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontConvertTextToGeom                                                M
*****************************************************************************/
int UserFontConvertTextToGeom(const UserFontText Text,
			      const UserFontName FontName,
			      UserFontStyleType FontStyle,
			      IrtRType FontSize,
			      IrtRType TextSpace,
			      UserFont3DEdgeType Text3DEdgeType,
			      const IrtRType Text3DSetup[2],
			      IrtRType Tolerance,
			      UserFontGeomOutputType OutputType,
			      IrtBType CompactOutput,
			      const char *PlacedTextBaseName,
			      IPObjectStruct **PlacedTextGeom,
			      char **ErrorStr)
{
    char Name[IRIT_LINE_LEN];
    int i, j;
    IPObjectStruct *TGeom, *PTmp1, *PTmp2, *PTmp3,
        *OutlineCurvesObjs;

    *ErrorStr = NULL;
    *PlacedTextGeom = NULL;

    if (PlacedTextBaseName != NULL)
	sprintf(Name, "%s_OL", PlacedTextBaseName);
    else
	strcpy(Name, "TxtOutline");

#   ifdef __WINNT__
        OutlineCurvesObjs = UserFontConvertFontToBezier(
				Text, FontName, FontStyle, TextSpace,
				OutputType != USER_FONT_OUTPUT_BEZIER_CRVS,
				Name);
	if (OutlineCurvesObjs == NULL) {
	    *ErrorStr = "Error in text conversion to geometry\n";
	    return FALSE;
	}
#   elif IRIT_HAVE_FREETYPE
        OutlineCurvesObjs = UserFontFTStringOutline2BezierCurves(
		                Text, FontName, TextSpace / 70.0,
				OutputType != USER_FONT_OUTPUT_BEZIER_CRVS,
				Name, (const char **) ErrorStr);
	if (OutlineCurvesObjs == NULL)
	    return FALSE;
#   else
	*ErrorStr = "No support for MS fonts or Freetype.\n";
	return FALSE;
#   endif /* __WINNT__ || IRIT_HAVE_FREETYPE */

    if (FontSize != 1.0) {
	IrtHmgnMatType Mat;

        if (IRIT_APX_EQ(FontSize, 0.0)) {
	    *ErrorStr = "Font size select is (almost) zero\n";
	    IPFreeObject(OutlineCurvesObjs);
	    return FALSE;	    
	}

	MatGenMatUnifScale(FontSize, Mat);
	TGeom = GMTransformObject(OutlineCurvesObjs, Mat);
	IPFreeObject(OutlineCurvesObjs);
	OutlineCurvesObjs = TGeom;	
    }

    if (CompactOutput) {
	IPObjectStruct *Letter;

        /* Go over the outline curves and merge all curves in a letter to */
        /* one curves object - was a list object.			  */
        for (i = 0;
	     (Letter = IPListObjectGet(OutlineCurvesObjs, i)) != NULL;
	     i++) {
	    if (IP_IS_CRV_OBJ(Letter)) {
	    }
	    else if (IP_IS_OLST_OBJ(Letter)) {
	        CagdCrvStruct *TCrv,
		    *Crvs = NULL;
		IPObjectStruct *CrvObj;

	        /* Convert to a curve object with several curves, in place. */
	        for (j = 0;
		     (CrvObj = IPListObjectGet(Letter, j)) != NULL;
		     j++) {
		    TCrv = CrvObj -> U.Crvs;
		    CrvObj -> U.Crvs = NULL;
		    IPFreeObject(CrvObj);
		    IRIT_LIST_PUSH(TCrv, Crvs);
		}
		/* Convert the list object to a curve object, in place. */
		IPListObjectInsert(Letter, 0, NULL);    /* Make list empty. */
		IPReallocNewTypeObject(Letter, IP_OBJ_CURVE);
		Letter -> U.Crvs = Crvs;
	    }
	    else
	        assert(0);
	}
    }

    switch (OutputType) {
        default:
        case USER_FONT_OUTPUT_BEZIER_CRVS:
        case USER_FONT_OUTPUT_BSPLINE_CRVS:
	    *PlacedTextGeom = OutlineCurvesObjs;
	    break;
        case USER_FONT_OUTPUT_FILLED2D_POLYS:
	    if (PlacedTextBaseName != NULL)
	        sprintf(Name, "%s_F2D", PlacedTextBaseName);
	    else
	        strcpy(Name, "TxtFille2D");
	    *PlacedTextGeom = UserFontBspList2Plgns(OutlineCurvesObjs,
						    Tolerance, Name);
	    IPFreeObject(OutlineCurvesObjs);
	    break;
        case USER_FONT_OUTPUT_OUTLINE_FILLED2D_POLYS:
	    if (PlacedTextBaseName != NULL)
	        sprintf(Name, "%s_F2D", PlacedTextBaseName);
	    else
	        strcpy(Name, "TxtPoly2D");
	    *PlacedTextGeom = IPGenListObject(PlacedTextBaseName, NULL, NULL);
	    TGeom = UserFontBspList2Plgns(OutlineCurvesObjs, Tolerance, Name);
	    for (i = 0;
		 (PTmp1 = IPListObjectGet(OutlineCurvesObjs, i)) != NULL &&
		 (PTmp2 = IPListObjectGet(TGeom, i)) != NULL;
		 i++) {
	        sprintf(Name, "%s_Mrg%d", PlacedTextBaseName, i);
		PTmp3 = IPGenListObject(Name, PTmp1, NULL);
		IPListObjectAppend(PTmp3, PTmp2);
		IPListObjectAppend(*PlacedTextGeom, PTmp3);
	    }
	    IPListObjectInsert(OutlineCurvesObjs, 0, NULL);
	    IPListObjectInsert(TGeom, 0, NULL);
	    IPFreeObject(TGeom);
	    IPFreeObject(OutlineCurvesObjs);
	    break;
        case USER_FONT_OUTPUT_SOLID3D_POLYS:
	    if (PlacedTextBaseName != NULL)
	        sprintf(Name, "%s_P3D", PlacedTextBaseName);
	    else
	        strcpy(Name, "TxtPoly3D");
	    *PlacedTextGeom = UserFontBspList2Solids(OutlineCurvesObjs,
						     Text3DEdgeType,
						     Text3DSetup[0],
						     Text3DSetup[1],
						     Tolerance, FALSE, Name);
	    IPFreeObject(OutlineCurvesObjs);
	    break;
	case USER_FONT_OUTPUT_FILLED2D_TSRFS:
	    if (PlacedTextBaseName != NULL)
	        sprintf(Name, "%s_T2D", PlacedTextBaseName);
	    else
	        strcpy(Name, "TxtTrim2D");
	    *PlacedTextGeom = UserFontBspList2TrimSrfs(OutlineCurvesObjs,
						       Tolerance, Name);
	    IPFreeObject(OutlineCurvesObjs);
	    break;
        case USER_FONT_OUTPUT_SOLID3D_TSRFS:
	    if (PlacedTextBaseName != NULL)
	        sprintf(Name, "%s_T3D", PlacedTextBaseName);
	    else
	        strcpy(Name, "TxtTrim3D");
	    *PlacedTextGeom = UserFontBspList2Solids(OutlineCurvesObjs,
						     Text3DEdgeType,
						     Text3DSetup[0],
						     Text3DSetup[1],
						     Tolerance, TRUE, Name);
	    IPFreeObject(OutlineCurvesObjs);
	    break;
    }

    return *PlacedTextGeom != NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Moving a B-spline Curve by Offset amount, in place, by moving all its    *
* control point in the curve normal direction, offset amount.		     *
*   Assuming CW for face out and CCW for face in.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   BspCrv:	B-spline curve to offset, in place.	                     *
*   Offset:	Offset amount.				                     *
*									     *
* RETURN VALUE:                                                              *
*   void					                             *
*****************************************************************************/
static void UserFontCrvOffset(CagdCrvStruct *Crv, CagdRType Offset)
{
    int i, Len;
    CagdRType
        Vec1[2],	/* Vec1 = Vector from current point to next point.  */
        Vec2[2],	/* Vec2 = Vector from prev. point to current point. */
        *Weights = CAGD_IS_RATIONAL_CRV(Crv) ? Crv -> Points[0] : NULL;

    if (Crv == NULL || Offset == 0)
	return;

    Len = Crv -> Length;
    Offset *= 0.5;

    /* Compute prev. unit normal. */
    Vec1[0] = Crv -> Points[1][Len - 2] - Crv -> Points[1][Len - 1];
    Vec1[1] = Crv -> Points[2][Len - 1] - Crv -> Points[2][Len - 2];
    IRIT_VEC2D_NORMALIZE(Vec1);

    for (i = 0; i < Len - 1; ++i) {
        IrtRType DProd, Scl;

        IRIT_VEC2D_COPY(Vec2, Vec1);

	/* Compute next. unit normal. */
	Vec1[0] = Crv -> Points[1][i] - Crv -> Points[1][i + 1];
	Vec1[1] = Crv -> Points[2][i + 1] - Crv -> Points[2][i];
	IRIT_VEC2D_NORMALIZE(Vec1);

	/*   Compute the (half the) angle between the two vectors but make  */
	/* sure we do no too sharp of a turn.				    */
	/*   The offset amount will be scaled by cos(acos(Angle) / 2).      */
	DProd = IRIT_DOT_PROD_2D(Vec1, Vec2);
	DProd = cos(acos(IRIT_BOUND(DProd, -1.0, 1.0)) * 0.5);

	Scl = (Weights ? Weights[i] : 1.0) * Offset /
	                                            IRIT_MAX(DProd, IRIT_EPS);
	Crv -> Points[1][i] += (Vec1[1] + Vec2[1]) * Scl;
	Crv -> Points[2][i] += (Vec1[0] + Vec2[0]) * Scl;
    }

    /* The last point is the same is the first one. */
    Crv -> Points[1][Len - 1] = Crv -> Points[1][0];
    Crv -> Points[2][Len - 1] = Crv -> Points[2][0];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converting a closed B-spline to a Polygon.	                             M
*                                                                            *
* PARAMETERS:                                                                M
*   BspCrv:	Pointer to the Curve.                                        M
*   Tolerance:	The tolerance used calculating the polygon.                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPPolygonStruct *:     Pointer to the new Polygon.                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontBspCrv2Poly                                                      M
*****************************************************************************/
IPPolygonStruct *UserFontBspCrv2Poly(CagdCrvStruct *BspCrv,
				     IrtRType Tolerance)
{
    IPPolygonStruct *Poly;
    IPVertexStruct *V;

    if (!BspCrv) 
	return NULL;

    Poly = IPCurve2Polylines(BspCrv, Tolerance, SYMB_CRV_APPROX_TOLERANCE);

    /* Because Object were built manually - need to add normals and plane. */
    Poly -> Plane[0] = Poly -> Plane[1] = Poly -> Plane[3] = 0.0;
    Poly -> Plane[2] = -1.0;
    IP_SET_PLANE_POLY(Poly);

    for (V = Poly -> PVertex; V != NULL; V = V -> Pnext) {
        IRIT_VEC_COPY(V -> Normal, Poly -> Plane);
	IP_SET_NORMAL_VRTX(V);
    };

    return Poly;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list object of Bezier curves into a list of Bspline curves    M
* by chaining adjacent Bezier curve segments into closed Bsplines loops.     M
*                                                                            *
* PARAMETERS:                                                                M
*   BzrListObj:   A list object of Bezier curves.                            M
*   HaveHoles:    Will be set to TRUE if holes are detected (examining the   M
*		  orientation of the loops).				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:    A list of Bspline loop curves.                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontBzrList2BspList                                                  M
*****************************************************************************/
CagdCrvStruct *UserFontBzrList2BspList(IPObjectStruct *BzrListObj,
				       IrtBType *HaveHoles)
{
    int j;
    IPObjectStruct *BzrCrvObj;
    CagdCrvStruct *BzrCrv, *TCrv,
        *BspCrv = NULL,
        *BspCrvList = NULL;

    *HaveHoles = FALSE;

    for (j = 0; (BzrCrvObj = IPListObjectGet(BzrListObj, j)) != NULL; j++) {
        BzrCrv = BzrCrvObj -> U.Crvs;

	if (BspCrv == NULL)
	    BspCrv = CagdCrvCopy(BzrCrv);
	else {
	    CagdPType Pt1, Pt2;

	    CagdCoerceToE2(Pt1, BspCrv -> Points, BspCrv -> Length - 1,
			   BspCrv -> PType);
	    CagdCoerceToE2(Pt2, BzrCrv -> Points, 0, BzrCrv -> PType);
	    if (IRIT_PT_APX_EQ_E2(Pt1, Pt2)) {
	        TCrv = CagdMergeCrvCrv(BspCrv, BzrCrv, FALSE);
		CagdCrvFree(BspCrv);
		BspCrv = TCrv;
	    }
	    else {
	        /* We getting here if the letter is build from more than   */
	        /* one loop.						   */
	        if (!TrimClassifyTrimCurveOrient(BspCrv)) {
		    /* The curve is CCW therefore it's a hole. */
		    *HaveHoles = TRUE;
		};
		IRIT_LIST_PUSH(BspCrv, BspCrvList);
		BspCrv = CagdCrvCopy(BzrCrv);
	    }
	}
    }

    if (BspCrv != NULL) {
        /* We getting here if the letter is build from more than one loop. */
        if (!TrimClassifyTrimCurveOrient(BspCrv))
	    *HaveHoles = TRUE;  /* The curve is CCW therefore it's a hole. */

	IRIT_LIST_PUSH(BspCrv, BspCrvList);
    }

    return BspCrvList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list object of Bspline merged closed curves (outline of the   M
* font) to polygons that fill the text.				    	     M
*                                                                            *
* PARAMETERS:                                                                M
*   BspListObj:  A list object of Bspline curves.                            M
*   Tol:	 Tolerance of approximation.				     M
*   Name:        Base name to use for constructed objects.  Can be NULL.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  The created polygons.		                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontBspList2TrimSrfs                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontBspList2Plgns                                                    M
*****************************************************************************/
IPObjectStruct *UserFontBspList2Plgns(IPObjectStruct *BspListObj,
				      IrtRType Tol,
				      const char *Name)
{
    char OName[IRIT_LINE_LEN_LONG];
    int i, j;
    IPPolygonStruct *Poly, *TPoly;
    IPObjectStruct *Letter, *CrvObj, *PlObj,
        *PolysObj = IPGenLISTObject(NULL);

    if (Name == NULL)
        Name = "Ltr";

    IP_SET_OBJ_NAME2(PolysObj, Name);

    for (i = 0;
	 (Letter = IPListObjectGet(BspListObj, i)) != NULL;
	 i++) {
        CagdCrvStruct *BspCrv, *TCrv,
	    *BspCrvList = NULL;

	/* Chain all curves forming the letter into one list. */
	if (IP_IS_CRV_OBJ(Letter)) {
	    BspCrvList = CagdCrvCopyList(Letter -> U.Crvs);
	}
	else if (IP_IS_OLST_OBJ(Letter)) {
	    for (j = 0;
		 (CrvObj = IPListObjectGet(Letter, j)) != NULL;
		 j++) {
	        TCrv = CagdCrvCopy(CrvObj -> U.Crvs);
		IRIT_LIST_PUSH(TCrv, BspCrvList);
	    }
	}
	else
	    assert(0);

	if (BspCrvList != NULL) {
	    for (BspCrv = BspCrvList, TPoly = NULL;
		 BspCrv != NULL;
		 BspCrv = BspCrv -> Pnext) {
	        Poly = UserFontBspCrv2Poly(BspCrv, Tol);
		IRIT_LIST_PUSH(Poly, TPoly);
	    }

	    /* Convert to polygons, using TPoly in place. */
	    PlObj = UserFontCnvrtPllns2Plgns(TPoly,
					     BspCrvList -> Pnext != NULL);

	    sprintf(OName, "%s_l%d", Name, i);
	    IP_SET_OBJ_NAME2(PlObj, OName);
	    IPListObjectAppend(PolysObj, PlObj);

	    CagdCrvFreeList(BspCrvList);
	}
    }

    return PolysObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list object of Bspline merged closed curves (outline of the   M
* font) to trimmed surfaces that fill the text.			    	     M
*                                                                            *
* PARAMETERS:                                                                M
*   BspListObj:  A list object of Bspline curves.                            M
*   Tol:	 Tolerance of approximation.				     M
*   Name:        Base name to use for constructed objects.  Can be NULL.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  The created trimmed surfaces.	                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontBspList2Plgns                                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontBspList2TrimSrfs                                                 M
*****************************************************************************/
IPObjectStruct *UserFontBspList2TrimSrfs(IPObjectStruct *BspListObj,
					 IrtRType Tol,
					 const char *Name)
{
    char OName[IRIT_LINE_LEN_LONG];
    int i, j;
    IPObjectStruct *Letter, *CrvObj, *TSrfObj,
        *TrimSrfsObj = IPGenLISTObject(NULL);

    if (Name == NULL)
        Name = "Ltr";

    IP_SET_OBJ_NAME2(TrimSrfsObj, Name);

    for (i = 0;
	 (Letter = IPListObjectGet(BspListObj, i)) != NULL;
	 i++) {
        CagdCrvStruct *TCrv,
	    *BspCrvList = NULL;

	/* Chain all curves forming the letter into one list. */
	if (IP_IS_CRV_OBJ(Letter)) {
	    BspCrvList = CagdCrvCopyList(Letter -> U.Crvs);
	}
	else if (IP_IS_OLST_OBJ(Letter)) {
	    for (j = 0;
		 (CrvObj = IPListObjectGet(Letter, j)) != NULL;
		 j++) {
	        TCrv = CagdCrvCopy(CrvObj -> U.Crvs);
		IRIT_LIST_PUSH(TCrv, BspCrvList);
	    }
	}
	else
	    assert(0);

	if (BspCrvList != NULL) {
	    TrimSrfStruct
	        *TSrf = UserFontCreateTrimmedSrfFromCurves(BspCrvList);

	    sprintf(OName, "%s_l%d", Name, i);
	    TSrfObj = IPGenTrimSrfObject(OName, TSrf, NULL);
	    IPListObjectAppend(TrimSrfsObj, TSrfObj);
	}
    }

    return TrimSrfsObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Constructs a trimemd surface representing one letter, out of a set of    *
* closed Bspline curves representing the outline of the letter.              *
*                                                                            *
* PARAMETERS:                                                                *
*   BspCrvList: List of closed B-spline curves, representing the outline.    *
*                                                                            *
* RETURN VALUE:                                                              *
*   TrimSrfStruct *:         A trimming surface filling the letter.          *
*****************************************************************************/
static TrimSrfStruct *UserFontCreateTrimmedSrfFromCurves(CagdCrvStruct
							          *BspCrvList)
{
    CagdRType Expand;
    CagdSrfStruct *Srf, *SrfAux;
    CagdBBoxStruct BBox;

    CagdCrvListBBox(BspCrvList, &BBox);
    Expand = IRIT_MAX(BBox.Max[0] - BBox.Min[0],
		      BBox.Max[1] - BBox.Min[1]) * 0.1;
    SrfAux = CagdPrimPlaneSrf(BBox.Min[0] - Expand,
			      BBox.Min[1] - Expand,
			      BBox.Max[0] + Expand,
			      BBox.Max[1] + Expand, 0.0);
    Srf = CagdCnvrtBzr2BspSrf(SrfAux);
    CagdSrfFree(SrfAux);

    BspKnotAffineTransOrder2(Srf -> UKnotVector,
			     Srf -> UOrder,
			     Srf -> ULength + Srf -> UOrder,
			     BBox.Min[0] - Expand,
			     BBox.Max[0] + Expand);
    BspKnotAffineTransOrder2(Srf -> VKnotVector,
			     Srf -> VOrder,
			     Srf -> VLength +  Srf -> VOrder,
			     BBox.Min[1] - Expand,
			     BBox.Max[1] + Expand);

    return TrimSrfNew2(Srf, BspCrvList, TRUE);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Convert closed polylines (outline of the font) to polygons that fill     *
* the text.                                                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   Pllns:      List of closed polylines to use in place and fill in.	     *
*               If HaveHoles, polylines are nested (with holes) and their    *
*               hierarchy must be classified. 				     *
*   HaveHoles:  True if polylines are nested in each other.                  *
*****************************************************************************/
static IPObjectStruct *UserFontCnvrtPllns2Plgns(IPPolygonStruct *Pllns,
						int HaveHoles)
{
    IPPolygonStruct *Pll;
    IPObjectStruct *RetVal, *PTmp;

    if (HaveHoles) {			    /* The letter have holes in it. */
        TrimCrvStruct *TCrv,
	    *LinCrvs = TrimPolylines2LinTrimCrvs(Pllns);

        IPFreePolygonList(Pllns);
	Pllns = NULL;

	TrimClassifyTrimmingLoops(&LinCrvs);

	for (TCrv = LinCrvs; TCrv != NULL; TCrv = TCrv -> Pnext) {
	    TrimCrvStruct
	        *NestedCrvs = (TrimCrvStruct *)
				AttrGetRefPtrAttrib(TCrv -> Attr, "_subTrims");

	    /* Convert nested polygons to simple polygons by connecting */
	    /* closest points on two polygons by a double bridge.       */
	    Pll = UserCnvrtPllns2PlgnsAux(TCrv, NestedCrvs);
	    Pllns = IPAppendPolyLists(Pll, Pllns);
	    TrimCrvFreeList(NestedCrvs);
	}
	TrimCrvFreeList(LinCrvs);
    }

    for (Pll = Pllns; Pll; Pll = Pll -> Pnext) {
	IPVertexStruct
	    *V = Pll -> PVertex;

	Pll -> Plane[0] = Pll -> Plane[1] = Pll -> Plane[3] = 0;
	Pll -> Plane[2] = -1;
	IP_SET_PLANE_POLY(Pll);

	do {
	    IRIT_VEC_COPY(V -> Normal, Pll -> Plane);
	    IP_SET_NORMAL_VRTX(V);

	    if (V -> Pnext == NULL)
		V -> Pnext = Pll -> PVertex; /* Verify we form closed loop. */
	    V = V -> Pnext;
	}
	while (V != Pll -> PVertex);
    }

    PTmp = IPGenPOLYObject(Pllns);
    RetVal = GMConvexPolyObjectN(PTmp);
    IPFreeObject(PTmp);

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Converted the possibly nested set of curve loops to a single simple      *
* polygon by created a two-sided bridge between the two closest points on    *
* the two polygons.                                                          *
*   Special care is made here not to reuse vertices that already participate *
* in previous bridges to prevent from crossings that yields wrong topology.  *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv:          Outer loop polygon.  Must exist.                           *
*   NestedCrvs:   Optional inner loop curves - can be NULL.                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPPolygonStruct *:    A single simple polygon with Crv as its outer loop *
*                         and NestedCrvs as its holes.                       *
*****************************************************************************/
static IPPolygonStruct *UserCnvrtPllns2PlgnsAux(const TrimCrvStruct *Crv,
						const TrimCrvStruct
						                  *NestedCrvs)
{
    IPVertexStruct *V1, *V1Next, *V2, *V2Next;
    CagdPolylineStruct
	*CagdPlln = CagdCnvrtLinBspCrv2Polyline(Crv -> TrimCrvSegList -> UVCrv);
    IPPolygonStruct *Pl2,
        *Pl = IPCagdPllns2IritPllns(CagdPlln);
    const TrimCrvStruct *NCrv;

    USER_FNT2BZR_RST_POLY_IGNORED_VRTX(Pl);

    for (NCrv = NestedCrvs; NCrv != NULL; NCrv = NCrv -> Pnext) {
	CagdPlln = CagdCnvrtLinBspCrv2Polyline(NCrv -> TrimCrvSegList -> UVCrv);
        Pl2 = IPCagdPllns2IritPllns(CagdPlln);
	USER_FNT2BZR_RST_POLY_IGNORED_VRTX(Pl2);

	GMDistPolyPoly(Pl, Pl2, &V1, &V2, USER_FNT2BZR_IGNORE_V_TAG);
	assert(V1 != NULL && V2 != NULL);
	USER_FNT2BZR_SET_IGNORED_VRTX(V1); /* Do not use these two vertices  */
	USER_FNT2BZR_SET_IGNORED_VRTX(V2); /* in further bridges.            */

	/* Make Pl2 closed. */
	IPGetLastVrtx(Pl2 -> PVertex) -> Pnext = Pl2 -> PVertex;

	/* Connect V1 to V2 and V2 to V1. */
	V1Next = IPCopyVertex(V1);
	V1Next -> Pnext = V1 -> Pnext;

	V2Next = IPCopyVertex(V2);
	V2Next -> Pnext = V2 -> Pnext;

	USER_FNT2BZR_SET_IGNORED_VRTX(V1Next); /* Do not use these two       */
	USER_FNT2BZR_SET_IGNORED_VRTX(V2Next);/* vertices in further bridges.*/

	V1 -> Pnext = V2Next;
	V2 -> Pnext = V1Next;
    }

    return Pl;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a list object of Bspline merged closed curves (outline of the   M
* font) to polygons that fill the text.				    	     M
*                                                                            *
* PARAMETERS:                                                                M
*   BspListObj:  A list object of Bspline curves.                            M
*   ExtStyle:    Type of 3D solid geometry: chamfered corners, rounded, etc. M
*   ExtOffset:   Size of the chamfer/rounding.				     M
*   ExtHeight:   Height of the 3D solid constructed text.		     M
*   Tol:	 Tolerance of approximation.				     M
*   GenTrimSrfs: TRUE to generate trimmed surfaces.  FALSE for polygons.     M
*   Name:        Base name to use for constructed objects.  Can be NULL.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  The created polygons.		                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontBspList2Solids                                                   M
*****************************************************************************/
IPObjectStruct *UserFontBspList2Solids(IPObjectStruct *BspListObj,
				       UserFont3DEdgeType ExtStyle, 
				       IrtRType ExtOffset, 
				       IrtRType ExtHeight,
				       IrtRType Tol,
				       CagdBType GenTrimSrfs,
				       const char *Name)
{
    int i, j;
    CagdVecStruct Vec;
    IPObjectStruct *Letter,
        *ResultObj = IPGenLISTObject(NULL);
    IrtVecType Translate;

    if (Name == NULL)
        Name = "Ltr";

    IP_SET_OBJ_NAME2(ResultObj, Name);

    IRIT_PT_RESET(Translate);
    IRIT_PT_RESET(Vec.Vec);

    for (i = 0;
	 (Letter = IPListObjectGet(BspListObj, i)) != NULL;
	 i++) {
        char OName[IRIT_LINE_LEN_LONG];
        CagdCrvStruct *BspCrv, *TCrv,
	    *BspCrvList = NULL;
	IPObjectStruct *CrvObj, *PlBaseObj, *PlBaseObj2,
	    *OneLetterGeom = IPGenLISTObject(NULL);

	sprintf(OName, "%s_l%d", Name, i);
	IP_SET_OBJ_NAME2(OneLetterGeom, OName);

	/* Chain all curves forming the letter into one list. */
	if (IP_IS_CRV_OBJ(Letter)) {
	    BspCrvList = CagdCrvCopyList(Letter -> U.Crvs);
	}
	else if (IP_IS_OLST_OBJ(Letter)) {
	    for (j = 0;
		 (CrvObj = IPListObjectGet(Letter, j)) != NULL;
		 j++) {
	        TCrv = CagdCrvCopy(CrvObj -> U.Crvs);
		IRIT_LIST_PUSH(TCrv, BspCrvList);
	    }
	}
	else
	    assert(0);

	if (BspCrvList != NULL) {
	    int SrfId = 0;
	    IPPolygonStruct *Poly, *Pl,
	        *Polys = NULL;

	    for (BspCrv = BspCrvList;
		 BspCrv != NULL;
		 BspCrv = BspCrv -> Pnext) {
	        int Pindex;
		CagdCrvStruct
		    *BspCrvForSrf1 = NULL,
		    *BspCrvForSrf2 = NULL;
		CagdSrfStruct *Srf;
		IPObjectStruct *SrfObj, *ExtObj, *PllnObj;

		Vec.Vec[2] = ExtHeight;

		if (ExtStyle == USER_FONT_3D_EDGE_NORMAL && !GenTrimSrfs) {
		    /* Construct the extruded sides as polygons. */
		    PllnObj = IPGenPOLYLINEObject(UserFontBspCrv2Poly(BspCrv,
								      Tol));

		    ExtObj = PrimGenEXTRUDEObject(PllnObj, Vec.Vec, FALSE);
		    sprintf(OName, "%s_l%dEX%d", Name, i, ++SrfId);
		    IP_SET_OBJ_NAME2(ExtObj, OName);

		    GMBlendNormalsToVertices(ExtObj -> U.Pl, 80);

		    IPFreeObject(PllnObj);
		}
		else {
		    /* Construct the extruded sides as surfaces. */

		    /* Get Srf orientation right. */
		    TCrv = CagdCrvReverse(BspCrv);
		    Srf = CagdExtrudeSrf(TCrv, &Vec);
		    CagdCrvFree(TCrv);
		    sprintf(OName, "%s_l%dEX%d", Name, i, ++SrfId);
		    ExtObj = IPGenSrfObject(OName, Srf, NULL);

		    /* Prepare two curves for bases, one possibly offset.*/
		    BspCrvForSrf1 = CagdCoerceCrvTo(BspCrv,
						    CAGD_PT_E3_TYPE, TRUE);

		    if (ExtStyle != USER_FONT_3D_EDGE_NORMAL)
		        UserFontCrvOffset(BspCrv, ExtOffset);
		    BspCrvForSrf2 = CagdCoerceCrvTo(BspCrv,
						    CAGD_PT_E3_TYPE, TRUE);
		    Translate[2] = -ExtOffset;
		    CagdCrvTransform(BspCrvForSrf2, Translate, 1.0);
		}
		IPListObjectAppend(OneLetterGeom, ExtObj);

		switch (ExtStyle) {
		    default:
		    case USER_FONT_3D_EDGE_NORMAL:
		        /* Do nothing if regular construction. */
			break;
		    case USER_FONT_3D_EDGE_CHAMFER:
		        /* Create the two chamfered (diagonal) surfaces. */
		        Srf = CagdRuledSrf(BspCrvForSrf1, BspCrvForSrf2, 2, 2);
			sprintf(OName, "%s_l%dCB%d", Name, i, ++SrfId);
			SrfObj = IPGenSrfObject(OName, Srf, NULL);
			IPListObjectAppend(OneLetterGeom, SrfObj);

			Translate[2] = ExtHeight;
			CagdCrvTransform(BspCrvForSrf1, Translate, 1.0);
			Translate[2] = ExtHeight + ExtOffset * 2;
			CagdCrvTransform(BspCrvForSrf2, Translate, 1.0);

			Srf = CagdRuledSrf(BspCrvForSrf2, BspCrvForSrf1, 2, 2);
			sprintf(OName, "%s_l%dCT%d", Name, i, ++SrfId);
			SrfObj = IPGenSrfObject(OName, Srf, NULL);
			IPListObjectAppend(OneLetterGeom, SrfObj);
			break;
		    case USER_FONT_3D_EDGE_ROUND:
		        /*   Create two ruled surface with three rows and   */
		        /* update  the middle row for the rounding.         */
		        Srf = CagdRuledSrf(BspCrvForSrf1, BspCrvForSrf2, 3, 3);
			for (Pindex = 0; Pindex < Srf -> ULength; ++Pindex) {
			    Srf -> Points[1][Pindex + Srf -> ULength] =
			                             Srf -> Points[1][Pindex];
			    Srf -> Points[2][Pindex + Srf -> ULength] =
			                             Srf -> Points[2][Pindex];
			    Srf -> Points[3][Pindex + Srf -> ULength] =
			       Srf -> Points[3][Pindex + Srf -> ULength * 2];
			}
			sprintf(OName, "%s_l%dRB%d", Name, i, ++SrfId);
			SrfObj = IPGenSrfObject(OName, Srf, NULL);
			IPListObjectAppend(OneLetterGeom, SrfObj);

			Translate[2] = ExtHeight;
			CagdCrvTransform(BspCrvForSrf1, Translate, 1.0);
			Translate[2] = ExtHeight + ExtOffset * 2;
			CagdCrvTransform(BspCrvForSrf2, Translate, 1.0);

			Srf = CagdRuledSrf(BspCrvForSrf2, BspCrvForSrf1, 3, 3);
			for (Pindex = 0; Pindex < Srf -> ULength; ++Pindex){
			    Srf -> Points[1][Pindex + Srf -> ULength] =
			        Srf -> Points[1][Pindex + Srf -> ULength * 2];
			    Srf -> Points[2][Pindex + Srf -> ULength] =
			        Srf -> Points[2][Pindex + Srf -> ULength * 2];
			    Srf -> Points[3][Pindex + Srf -> ULength] =
			                           Srf -> Points[3][Pindex];
			}
			sprintf(OName, "%s_l%dRT%d", Name, i, ++SrfId);
			SrfObj = IPGenSrfObject(OName, Srf, NULL);
			IPListObjectAppend(OneLetterGeom, SrfObj);
			break;
		}
		        
		/* Convert the curve to a polyline as a first stage to      */
		/* polygonalize the bases.				    */
		if (!GenTrimSrfs) {
		    Poly = UserFontBspCrv2Poly(BspCrv, Tol);
		    IRIT_LIST_PUSH(Poly, Polys);
		}

		CagdCrvFree(BspCrvForSrf1);
		CagdCrvFree(BspCrvForSrf2);
	    }

	    /* Convert all polylines of the base to polygons, using Polys,  */
	    /* in place.					            */
	    if (GenTrimSrfs) {
	        CagdRType Translate[3];

	        PlBaseObj = IPGenTRIMSRFObject(
				  UserFontCreateTrimmedSrfFromCurves(
						CagdCrvCopyList(BspCrvList)));
		PlBaseObj2 = IPReverseObject(PlBaseObj);

		/* Move the two bases to the right height. */
		Translate[0] = Translate[1] = 0.0;
		if (ExtStyle == USER_FONT_3D_EDGE_NORMAL) 
		    Translate[2] = ExtHeight;
		else
		    Translate[2] = ExtHeight + ExtOffset;
		CagdSrfTransform(PlBaseObj -> U.TrimSrfs -> Srf,
				 Translate, 1.0);

		if (ExtStyle != USER_FONT_3D_EDGE_NORMAL) {
		    Translate[2] = -ExtOffset;
		    CagdSrfTransform(PlBaseObj2 -> U.TrimSrfs -> Srf,
				     Translate, 1.0);
		}
	    }
	    else {
	        PlBaseObj = UserFontCnvrtPllns2Plgns(Polys,
						 BspCrvList -> Pnext != NULL);
		PlBaseObj2 = IPReverseObject(PlBaseObj);

		/* Update the two bases and move to the right height. */
		for (Pl = PlBaseObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
		    IPVertexStruct
		        *V = Pl -> PVertex;

		    do {
		        if (ExtStyle == USER_FONT_3D_EDGE_NORMAL) 
			    V -> Coord[2] = ExtHeight;
			else
			    V -> Coord[2] = ExtHeight + ExtOffset;

			V = V -> Pnext;
		    }
		    while (V != NULL && V != Pl -> PVertex);
		}
		for (Pl = PlBaseObj2 -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
		    IPVertexStruct
		        *V = Pl -> PVertex;

		    do {
		        V -> Normal[2] = 1.0;
			if (ExtStyle != USER_FONT_3D_EDGE_NORMAL) 
			    V -> Coord[2] = -ExtOffset;

			V = V -> Pnext;
		    }
		    while (V != NULL && V != Pl -> PVertex);
		}

		/* Update the planes of the polygons. */
		IPUpdatePolyPlane(PlBaseObj -> U.Pl);
		IPUpdatePolyPlane(PlBaseObj2 -> U.Pl);
	    }

	    /* Place the bases in the output list. */
	    sprintf(OName, "%s_l%dBT%d", Name, i, ++SrfId);
	    IP_SET_OBJ_NAME2(PlBaseObj, OName);
	    IPListObjectAppend(OneLetterGeom, PlBaseObj);

	    sprintf(OName, "%s_l%dBB%d", Name, i, ++SrfId);
	    IP_SET_OBJ_NAME2(PlBaseObj2, OName);
	    IPListObjectAppend(OneLetterGeom, PlBaseObj2);

	    CagdCrvFreeList(BspCrvList);
	}

	IPListObjectAppend(ResultObj, OneLetterGeom);
    }

    return ResultObj;
}
