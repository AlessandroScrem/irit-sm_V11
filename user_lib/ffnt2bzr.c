/******************************************************************************
* FFnt2Bzr.c - extracting system fonts as Bezier curves using the Freetype    *
* library, see www.freetype.org.					      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Nov 2010.					      *
******************************************************************************/

#ifdef IRIT_HAVE_FREETYPE

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <ft2build.h>
#include FT_FREETYPE_H
#include <freetype/ftoutln.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/miscattr.h"
#include "inc_irit/cagd_lib.h"
#include "user_loc.h"

#define USER_FONT_FT_IS_OFF_CRV_PT(Tag)	   (((Tag) & FT_CURVE_TAG_ON) == 0)
#define USER_FONT_FT_IS_CONIC_CRV_PT(Tag)  (((Tag) & FT_CURVE_TAG_CONIC) == 0)

typedef struct UserFontFTOutlineDecomposeDataStruct {
    CagdRType Dx;
    CagdRType Dy;
    CagdRType Scale;
} UserFontFTOutlineDecomposeDataStruct;

static FT_Vector GlblLastPos;
static CagdCrvStruct
    *GlblAllCrvs = NULL;

static int UserFontFTSetStartPoint(const FT_Vector *v1, void *User);
static int UserFontFTGenLinearBezier(const FT_Vector *v2, void *User);
static int UserFontFTGenQuadraticBezier(const FT_Vector *v2,
					const FT_Vector *v3,
					void *User);
static int UserFontFTGenCubicBezier(const FT_Vector *v2,
				    const FT_Vector *v3,
				    const FT_Vector *v4,
				    void *User);
static IPObjectStruct *UserFontFTOneCharOutline2BezierCurves(
						     FT_Outline *Outline,
						     CagdRType XShift,
						     CagdRType YShift,
						     CagdRType Scale,
						     int CharIndexInString,
						     int CharIndex,
						     const char *RootObjName);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Set the starting position.                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   v1:     Second control points of the linear Bezier from GlblLastPos.     *
*   User:   User data pointer.					             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    Returned error, 0 for success.                                   *
*****************************************************************************/
static int UserFontFTSetStartPoint(const FT_Vector *v1, void *User)
{
    GlblLastPos = *v1;

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create a linear Bezier.                                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   v2:     Second control points of the linear Bezier from GlblLastPos.     *
*   User:   User data pointer.					             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    Returned error, 0 for success.                                   *
*****************************************************************************/
static int UserFontFTGenLinearBezier(const FT_Vector *v2, void *User)
{
    UserFontFTOutlineDecomposeDataStruct
	*Data = (UserFontFTOutlineDecomposeDataStruct *) User;
    CagdCrvStruct
        *Crv = BzrCrvNew(2, CAGD_PT_E2_TYPE);

    Crv -> Points[1][0] = GlblLastPos.x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][0] = GlblLastPos.y * Data -> Scale + Data -> Dy;

    Crv -> Points[1][1] = v2 -> x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][1] = v2 -> y * Data -> Scale + Data -> Dy;

    GlblLastPos = *v2;

    IRIT_LIST_PUSH(Crv, GlblAllCrvs);

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create a quadratic Bezier.                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   v2, v3:   Ctl pts (with GlblLastPos as first one) of quadratic Bezier.   *
*   User:     User data pointer.				             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    Returned error, 0 for success.                                   *
*****************************************************************************/
static int UserFontFTGenQuadraticBezier(const FT_Vector *v2,
					const FT_Vector *v3,
					void *User)
{
    UserFontFTOutlineDecomposeDataStruct
	*Data = (UserFontFTOutlineDecomposeDataStruct *) User;
    CagdCrvStruct
        *Crv = BzrCrvNew(3, CAGD_PT_E2_TYPE);

    Crv -> Points[1][0] = GlblLastPos.x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][0] = GlblLastPos.y * Data -> Scale + Data -> Dy;

    Crv -> Points[1][1] = v2 -> x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][1] = v2 -> y * Data -> Scale + Data -> Dy;

    Crv -> Points[1][2] = v3 -> x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][2] = v3 -> y * Data -> Scale + Data -> Dy;

    GlblLastPos = *v3;

    IRIT_LIST_PUSH(Crv, GlblAllCrvs);

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create a cubic Bezier.                                                   *
*                                                                            *
* PARAMETERS:                                                                *
*   v2, v3, v4:   Ctl pts (with GlblLastPos as first one) of cubic Bezier.   *
*   User:         User data pointer.				             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    Returned error, 0 for success.                                   *
*****************************************************************************/
static int UserFontFTGenCubicBezier(const FT_Vector *v2,
				    const FT_Vector *v3,
				    const FT_Vector *v4,
				    void *User)
{
    UserFontFTOutlineDecomposeDataStruct
	*Data = (UserFontFTOutlineDecomposeDataStruct *) User;
    CagdCrvStruct
        *Crv = BzrCrvNew(4, CAGD_PT_E2_TYPE);

    Crv -> Points[1][0] = GlblLastPos.x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][0] = GlblLastPos.y * Data -> Scale + Data -> Dy;

    Crv -> Points[1][1] = v2 -> x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][1] = v2 -> y * Data -> Scale + Data -> Dy;

    Crv -> Points[1][2] = v3 -> x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][2] = v3 -> y * Data -> Scale + Data -> Dy;

    Crv -> Points[1][3] = v4 -> x * Data -> Scale + Data -> Dx;
    Crv -> Points[2][3] = v4 -> y * Data -> Scale + Data -> Dy;

    GlblLastPos = *v4;

    IRIT_LIST_PUSH(Crv, GlblAllCrvs);

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Given an outline of some character in Freetype outline style - convert   *
* to a set of curves in IRIT format.					     *
*                                                                            *
* PARAMETERS:                                                                *
*   Outline:      Outline to convert to IRIT format.                         *
*   XShift:       To apply to this character.				     *
*   YShift:       To apply to this character.				     *
*   Scale:        To apply to the character.				     *
*   CharIndexInString:  Index of this character in the converted string.     *
*   CharIndex:    Character code (I.e. Ascii).				     *
*   RootObjName:  Name of root object.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:  List of curves in IRIT format representing the char.  *
*****************************************************************************/
static IPObjectStruct *UserFontFTOneCharOutline2BezierCurves(
						     FT_Outline *Outline,
						     CagdRType XShift,
						     CagdRType YShift,
						     CagdRType Scale,
						     int CharIndexInString,
						     int CharIndex,
						     const char *RootObjName)
{
    static FT_Outline_Funcs
        OutlineFunc = {
            UserFontFTSetStartPoint,
	    UserFontFTGenLinearBezier,
	    UserFontFTGenQuadraticBezier,
	    UserFontFTGenCubicBezier,
	    0,
	    0
    };
    int i;
    char Name[IRIT_LINE_LEN];
    UserFontFTOutlineDecomposeDataStruct Data;
    IPObjectStruct *Letter, *PTmp;

    GlblAllCrvs = NULL;
    Data.Dx = XShift;
    Data.Dy = YShift;
    Data.Scale = Scale * 1e-4;

    if (FT_Outline_Decompose(Outline, &OutlineFunc, &Data))
	return NULL;

    Letter = IPLnkListToListObject(CagdListReverse(GlblAllCrvs),
				   IP_OBJ_CURVE);

    assert(IP_IS_OLST_OBJ(Letter));
    for (i = 0; (PTmp = IPListObjectGet(Letter, i++)) != NULL; ) {
        if (CharIndex > ' ' && CharIndex <= '~')
	    sprintf(Name, "%s_%c%d_%d",
		    RootObjName, CharIndex, CharIndexInString, i);
	else
	    sprintf(Name, "%s_U%d_%d_%d",
		    RootObjName, CharIndex, CharIndexInString, i);
	IP_SET_OBJ_NAME2(PTmp, Name);
    }

    if (CharIndex > ' ' && CharIndex <= '~')
        sprintf(Name, "%s_%c%d", RootObjName, CharIndex, CharIndexInString);
    else
        sprintf(Name, "%s_U%d_%d", RootObjName, CharIndex, CharIndexInString);
    IP_SET_OBJ_NAME2(Letter, Name);

    return Letter;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given a string and a font name - convert to a set of curves in IRIT      M
* format.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:          String to convert to Bezier outline.                      M
*   FontName:      Font to read the outline from.                            M
*   Spacing:       To apply between the characters.			     M
*   MergeToBsp:    TRUE to merge the Bezier curves into larger B-spline      M
*                  curves, FALSE to leave the original (font) Bezier curves. M
*   RootObjName:   Name of root object.					     M
*   ErrStr:        A string describing the error, if any.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *: List of curves in IRIT format representing the string. M
*                                                                            *
* SEE ALSO:                                                                  M
*   Freetype font library, http://www.freetype.org.                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontFTStringOutline2BezierCurves                                     M
*****************************************************************************/
IPObjectStruct *UserFontFTStringOutline2BezierCurves(
                                                const UserFontText Text,
						const UserFontName FontName,
						IrtRType Spacing,
						int MergeToBsp,
						const char *RootObjName,
						const char **ErrStr)
{
    int i, m, n;
    char ObjName[IRIT_LINE_LEN];
    FT_Error error;
    FT_Face face;
    FT_Library library;
    FT_GlyphSlot slot;
    FT_Matrix matrix;                              /* Transformation matrix. */
    FT_Vector pen;                                  /* Untransformed origin. */
    IPObjectStruct *PCrvs, *Letter;

    *ErrStr = NULL;

    error = FT_Init_FreeType(&library);               /* Initialize library. */
    if (error)
        return NULL;

    error = FT_New_Face(library, FontName, 0, &face );/* Create face object. */
    switch (error) {
        case 0:
	    break;
        case FT_Err_Unknown_File_Format:
	    *ErrStr = "Unknown font file format ignored";
	    return NULL;
        default:
	    *ErrStr = "Failed to open font file";
	    return NULL;
    }

    /* use 50pt at 100dpi */
    error = FT_Set_Char_Size(face, 50 * 64, 0,
			     100, 0 );                /* Set character size. */
    if (error) {
        *ErrStr = "Failed to scale font";
        return NULL;
    }

    slot = face -> glyph;

    /* Set up matrix. */
    matrix.xx = (FT_Fixed)(1 * 0x10000L);
    matrix.xy = (FT_Fixed)(0 * 0x10000L);
    matrix.yx = (FT_Fixed)(0 * 0x10000L);
    matrix.yy = (FT_Fixed)(1 * 0x10000L);

    /* The pen position in 26.6 cartesian space coordinates; */
    /* start at (300, 0) relative to the upper left corner   */
    pen.x = 0;
    pen.y = 0;

    /* Build a name for the object. */
    for (m = 0; m < 15; m++) {
        if (!isalnum(FontName[m]))
	    break;
	ObjName[m] = FontName[m];
    }
    ObjName[m++] = '_';
    for (i = 0, n = m; n < 30; ) {
        if (!isalnum(Text[i]))
	    break;
        ObjName[n++] = Text[i++];
    }
    ObjName[n++] = 0;

    PCrvs = IPGenListObject(ObjName, NULL, NULL);

    for (n = 0; n < (int) strlen(Text); n++) {
        /* Set transformation. */
        FT_Set_Transform(face, &matrix, &pen);

	/* Load glyph image into the slot (erase previous one). */
        error = FT_Load_Char(face, Text[n], FT_LOAD_DEFAULT);
	if (error) {
	    *ErrStr = "Failed to extract outline of character";
	    continue;                                      /* Ignore errors. */
	}

	Letter = UserFontFTOneCharOutline2BezierCurves(&slot -> outline,
					       n * Spacing, 0.0, 1.0,
					       n, Text[n], RootObjName);

	if (MergeToBsp) {
	    int j;
	    IrtBType HaveHoles;
	    IPObjectStruct *TmpObj, *Symbol;
	    CagdCrvStruct *Crv,
	        *Crvs = UserFontBzrList2BspList(Letter, &HaveHoles);

	    if (Crvs == NULL) {
		/* Ignore this one. */
	    }
	    else if (Crvs -> Pnext != NULL) {
		char Name[IRIT_LINE_LEN];

	        TmpObj = IPGenListObject(IP_GET_OBJ_NAME(Letter), NULL, NULL);
		IPFreeObject(Letter);
		Letter = TmpObj;

		j = 0;
		while (Crvs != NULL) {
		    IRIT_LIST_POP(Crv, Crvs);
		    TmpObj = IPGenCrvObject(IP_GET_OBJ_NAME(Letter),
					    Crvs, NULL);

		    if (Text[n] > ' ' && Text[n] <= '~')
		        sprintf(Name, "%s_%c%d_%d",
				RootObjName, Text[n], n, j);
		    else
		        sprintf(Name, "%s_U%d_%d_%d",
				RootObjName, Text[n], i, j);
		    Symbol = IPGenCrvObject(Name, Crv, NULL);
		    IPListObjectInsert(Letter, j++, Symbol);	   
		}
		IPListObjectInsert(Letter, j, NULL);
	    }
	    else {
	        TmpObj = IPGenCrvObject(IP_GET_OBJ_NAME(Letter), Crvs, NULL);
		IPFreeObject(Letter);
		Letter = TmpObj;
	    }
	}

	IPListObjectAppend(PCrvs, Letter);

	/* Increment pen position. */
	pen.x += slot -> advance.x;
	pen.y += slot -> advance.y;
    }

    FT_Done_Face(face);
    FT_Done_FreeType(library);

    return PCrvs;
}

#endif /* IRIT_HAVE_FREETYPE */
