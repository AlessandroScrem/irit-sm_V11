/******************************************************************************
* FontLOut.cpp - Performs font layout over a given polyline.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Dec 2013.					      *
******************************************************************************/

#include <ctype.h>

#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/geom_lib.h"
#include "user_loc.h"

#define USER_FONT_LOUT_MAX_INTERS	10

static int UserFontBoundaryLineInter(const IPPolygonStruct *BoundingPoly,
				     const IrtPtType LinePt,
				     IrtRType *InterXVals);

#if defined(ultrix) && defined(mips)
static int UserFontCmpReals(VoidPtr PReal1, VoidPtr PReal2);
#else
static int UserFontCmpReals(const VoidPtr PReal1, const VoidPtr PReal2);
#endif /* ultrix && mips (no const support) */

#ifdef __WINNT__

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts a wide (multi-short) string to an ascii (multi-byte) string.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Str:       Input wide (multi-short) string to convert to a               M
*              multi-byte ascii string.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   char *:  Converted string, allocated dynamically, or NULL if error.      M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserAscii2WChar                                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserWChar2Ascii                                                          M
*****************************************************************************/
char *UserWChar2Ascii(const UserFontText Str)
{
    int l,
        Len = (int) wcslen(Str);
    char
        *RetVal = (char *) malloc(sizeof(char) * Len + 1);

    if ((l = (int) wcstombs(RetVal, Str, Len)) > 0) {
        RetVal[l] = 0;
        return RetVal;
    }
    else {
        free(RetVal);
	return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts an ascii (multi-byte) string to a wide (multi-short) string.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Str:       Input multi byte ascii string to convert to a wide string.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   UserFontText:  Converted wide string, allocated dynamically.             M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserWChar2Ascii                                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserAscii2WChar                                                          M
*****************************************************************************/
UserFontText UserAscii2WChar(const char *Str)
{
    int Len = (int) strlen(Str);
    UserFontText
        RetVal = (UserFontText) malloc(sizeof(wchar_t) * Len + 2);

    mbstowcs(RetVal, Str, Len);
    RetVal[Len] = 0;

    return RetVal;
}

#endif /* __WINNT__ */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Layout the given text over the given bounding region.  Constructed Text  M
* will be controlled by FontName, FontStyle, FontSize, and FontSpace and     M
* will be aligned to follow TextAlignment.				     M
*   All construction is conducted over the XY plane, in 2D.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:            The text to layout.  A string.			     M
*   FontName:        The font name to use.  I.e. "Times new Roman".          M
*   FontStyle:       Font style to use.  I.e. italics.			     M
*   FontSize:        The size of the constructed text.                       M
*   FontSpace:       (WordWidth, SpaceWidth, LineHeight) spacing, specified  M
*                    in text font's point units.			     M
*   Tolerance:	     For 2D filled polygons and 3D solid text geometry.      M
*   Text3DEdge:      For 3D text geometry controls edges (i.e. chamferred.). M
*   Text3DSetup:     For 3D text, a vector of (Chamfer offset size, 3D       M
*                    extruded vertical distance).                            M
*   Alignment:       Text alignment, left, centered, etc.                    M
*   BoundingCrv:     A closed curve to place the text inside.                M
*   OutputType:      Output should be original Bezier curves, merged         M
*                    B-spline curves, filled planar polys, 3D polys, etc.    M
*   PlacedTextGeom:  The placed (geometry of the) text in.  Actually only    M
*                    the base line of the text is guaranteed to be in.       M
*		     NULL if error.					     M
*   ErrorStr:        Will be updated with an error description if was one,   M
*                    or with NULL if all went well.			     M
*									     M
* RETURN VALUE:                                                              M
*   int:     TRUE if all text fitted in the BoundingPoly, FALSE otherwise.   *
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontLayoutOverShape                                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontLayoutOverShape2                                                 M
*****************************************************************************/
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
			     char **ErrorStr)
{
    IPPolygonStruct
        *BoundingPoly = IPCurve2Polylines(BoundingCrv, Tolerance,
					  SYMB_CRV_APPROX_TOLERANCE);
    int RetVal = UserFontLayoutOverShape(Text, FontName, FontStyle, FontSize,
					 FontSpace, Tolerance,
					 Text3DEdge, Text3DSetup,
					 Alignment, BoundingPoly,
					 OutputType, PlacedTextGeom, ErrorStr);

    IPFreePolygon(BoundingPoly);

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Layout the given text over the given bounding region.  Constructed Text  M
* will be controlled by FontName, FontStyle, FontSize, and FontSpace and     M
* will be aligned to follow TextAlignment.				     M
*   All construction is conducted over the XY plane, in 2D.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:            The text to layout.  A string.			     M
*   FontName:        The font name to use.  I.e. "Times new Roman".          M
*   FontStyle:       Font style to use.  I.e. italics.			     M
*   FontSize:        The size of the constructed text.                       M
*   FontSpace:       (WordWidth, SpaceWidth, LineHeight) spacing, specified  M
*                    in text font's point units.			     M
*   Tolerance:	     For 2D filled polygons and 3D solid text geometry.      M
*   Text3DEdge:      For 3D text geometry controls edges (i.e. chamferred.). M
*   Text3DSetup:     For 3D text, a vector of (Chamfer offset size, 3D       M
*                    extruded vertical distance).                            M
*   Alignment:   Text alignment, left, centered, etc.                    M
*   BoundingPoly:    A closed region to place the text inside.               M
*   OutputType:      Selects the type of geometry to create:                 M
*                    1. Outline Bezier curves as in the font data.           M
*                    2. Outline B-spline curves (merging Bezier curves in 1).M
*                    3. Filling 2D polygons, for the outline geometry.       M
*                    4. Full 3D text.                                        M
*   PlacedTextGeom:  The result synthesized placed (geometry of the) text.   M
*                    Actually only the base line of the text is guaranteed   M
*                    to be in.						     M
*		     NULL if error.					     M
*   ErrorStr:        Will be updated with an error description if was one,   M
*                    or with NULL if all went well.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if all text fitted in the BoundingPoly, FALSE otherwise.   M
*            FALSE might also be returned if we failed to generate the text, M
*            in which case ErrorStr will hold some error description.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontLayoutOverShape2                                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontLayoutOverShape                                                  M
*****************************************************************************/
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
			    char **ErrorStr)
{
    int FittedInShape;
    UserFontDimInfoStruct FontDims;
    UserFontWordLayoutStruct
        *Words = UserFontLayoutOverShapeGenWords(
				    Text, FontName, FontStyle, 1.0,
				    FontSpace, Tolerance, Text3DEdge,
				    Text3DSetup, Alignment,
				    BoundingPoly, OutputType, FALSE,
				    "Text", &FontDims, ErrorStr);

    if (Words == NULL) {
        *PlacedTextGeom = NULL;
        return FALSE;
    }

    /* Time to place the text in the prescribed shape. */
    FittedInShape = UserFontLayoutOverShapePlaceWords(
					       Words, FontSize, FontSpace,
					       Alignment, &FontDims,
					       BoundingPoly, PlacedTextGeom);

    /* Free the data structures. */
    UserFontLayoutOverShapeFree(Words);

    return FittedInShape;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Free the allocated words and their geometry.                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Words:   Linked list of words to free.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontLayoutOverShape, UserFontLayoutOverShape2,                       M
*   UserFontLayoutOverShapeGenWords, UserFontLayoutOverShapePlaceWords.      M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontLayoutOverShapeFree                                              M
*****************************************************************************/
void UserFontLayoutOverShapeFree(struct UserFontWordLayoutStruct *Words)
{
    while (Words != NULL) {
        UserFontWordLayoutStruct *W;

        IRIT_LIST_POP(W, Words);
	IPFreeObject(W -> Geom);
	free(W -> Word);
	IritFree(W -> FontName);
	IritFree(W);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructed words as a preparation to layout the text over the given     M
* bounding region. 							     M
*   All construction is conducted over the XY plane, in 2D.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:            The text to layout.  A string.			     M
*   FontName:        The font name to use.  I.e. "Times new Roman".          M
*   FontStyle:       Font style to use.  I.e. italics.			     M
*   FontSize:        The size of the constructed text.                       M
*   FontSpace:       (WordWidth, SpaceWidth, LineHeight) spacing, specified  M
*                    in text font's point units.			     M
*   Tolerance:	     For 2D filled polygons and 3D solid text geometry.      M
*   Text3DEdge:      For 3D text geometry controls edges (i.e. chamferred.). M
*   Text3DSetup:     For 3D text, a vector of (Chamfer offset size, 3D       M
*                    extruded vertical distance).                            M
*   Alignment:       Text alignment, left, centered, etc.                    M
*   BoundingPoly:    A closed region to place the text inside.               M
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
*   CompactOutput:   If TRUE, merged trimming surfaces or polygons into      M
*                    single objects.					     M
*   OutputBaseName:  Base name to use in output geometry.		     M
*   FontDims:        Structure to be updated with font dimensions if all     M
*                    went well.						     M
*   ErrorStr:        Will be updated with an error description if was one,   M
*                    or with NULL if all went well.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   UserFontWordLayoutStruct *:     Generated words.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontLayoutOverShape, UserFontLayoutOverShape2,                       M
*   UserFontLayoutOverShapeFree, UserFontLayoutOverShapePlaceWords.          M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontLayoutOverShapeGenWords                                          M
*****************************************************************************/
UserFontWordLayoutStruct *UserFontLayoutOverShapeGenWords(
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
				    char **ErrorStr)
{
    int i, s, e, Len;
    UserFontText
	TextDup = USER_FONT_STR_DUP(Text);
    UserFontWordLayoutStruct *W,
        *Words = NULL;
    IPObjectStruct *Geom;

    *ErrorStr = NULL;

    /* Decompose string into words at white spaces: */
    s = 0;
    Len = (int) USER_FONT_STR_LEN(TextDup);
    do {
        /* Get to the first non-space character. */
	while (s < Len && USER_FONT_IS_SPACE(TextDup[s]))
	    s++;
	/* March until a space character is found.*/
	for (e = s; e < Len && !USER_FONT_IS_SPACE(TextDup[e]); e++);

	if (e > s) {
	    W = (UserFontWordLayoutStruct *)
	                         IritMalloc(sizeof(UserFontWordLayoutStruct));
	    IRIT_ZAP_MEM(W, sizeof(UserFontWordLayoutStruct));
	    W -> NewLine = TextDup[e] == '\n' || TextDup[e] == '\r';
	    TextDup[e] = 0;
	    W -> Word = USER_FONT_STR_DUP(&TextDup[s]);
	    W -> FontName = IritStrdup(FontName);
	    W -> FontStyle = FontStyle;
	    W -> RelSize = FontSize;
	    W -> Font3DEdge = Text3DEdge;
	    IRIT_PT2D_COPY(W -> Font3DOptions, Text3DSetup);
	    W -> FontAlignment = Alignment;

	    IRIT_LIST_PUSH(W, Words);
	}
	s = e + 1;
    }
    while (s < Len);

    Words = CagdListReverse(Words);
    free(TextDup);

#   ifdef DEBUG_USER_TEXT_LOUT_SHAPE
        for (W = Words; W != NULL; W = W -> Pnext) {
	    fprintf(stderr, "Word = \"%s\" [NL = %d]\n",
		    USER_FONT_GET_WORD_ASCII(W -> Word), W -> NewLine);
	}
#   endif /* DEBUG_USER_TEXT_LOUT_SHAPE */

    /* Place the word "shp" as first to measure the descent, base, mean,    */
    /* and ascent lines of this fonts (four different elevation).           */
    W = (UserFontWordLayoutStruct *)
	                         IritMalloc(sizeof(UserFontWordLayoutStruct));
    IRIT_ZAP_MEM(W, sizeof(UserFontWordLayoutStruct));
    W -> Word = USER_FONT_STR_DUP(USER_FONT_STR_CONST("shp"));
    IRIT_LIST_PUSH(W, Words);

    /* Convert the text to geometry and compute bounding boxes: */
    for (i = 0, W = Words; W != NULL; i++, W = W -> Pnext) {
        char BaseName[IRIT_LINE_LEN];

	/* Recall that i = 0 is for our "shp" pushed word: */
	sprintf(BaseName, "%s_%d", OutputBaseName, i);

	if (!UserFontConvertTextToGeom(W -> Word, FontName, FontStyle,
				       FontSize, FontSpace[0],
				       Text3DEdge, Text3DSetup,
				       Tolerance, OutputType, CompactOutput,
				       BaseName, &W -> Geom, ErrorStr)) {
	    return FALSE;
	}

	W -> BBox = *GMBBComputeBboxObject(W -> Geom);

#       ifdef DEBUG_USER_TEXT_LOUT_SHAPE
	    fprintf(stderr, "Word = \"%s\", BBox = [%.6f %.6f] : [%.6f %.6f]\n",
		    USER_FONT_GET_WORD_ASCII(W -> Word),
		    W -> BBox.Min[0], W -> BBox.Min[1],
		    W -> BBox.Max[0], W -> BBox.Max[1]);
#       endif /* DEBUG_USER_TEXT_LOUT_SHAPE */
    }

    /* Compute the descent, base, mean, and ascent lines of this fonts: */
    IRIT_LIST_POP(W, Words);
    FontDims -> DescentLineHeight = W -> BBox.Min[1];
    FontDims -> AscentLineHeight = W -> BBox.Max[1];
    Geom = IPListObjectGet(W -> Geom, 0);                 /* Fetch the 's'. */
    FontDims -> BBox = *GMBBComputeBboxObject(Geom);
    FontDims -> BaseLineHeight = FontDims -> BBox.Min[1];
    FontDims -> MeanLineHeight = FontDims -> BBox.Max[1];
    FontDims -> SpaceWidth = (FontDims -> BBox.Max[0] -
			      FontDims -> BBox.Min[0]) * 0.75;
    UserFontLayoutOverShapeFree(W);

    return Words;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Translate the given (geometry of) words so they fit into the given       M
* bounding shape BoundingPoly.                                               M
*                                                                            *
* PARAMETERS:                                                                M
*   Words:           The ordered list of words to place into bounding shape  M
*   FontSize:        The size of the constructed text.                       M
*   FontSpace:       (WordWidth, SpaceWidth, LineHeight) spacing, specified  M
*                    in text font's point units.			     M
*   Alignment:       Align each line of the placed words to the left,        M
*		     centered, etc.					     M
*   FontDims:        Dimensions information on this font.		     M
*   BoundingPoly:    The region to place the text in (actually only base     M
*                    line will be in).					     M
*   PlacedTextGeom:  Where to place the created geometry.		     M
*                       Organized as an object per character, grouped in     M
*                    list objects as words, all grouped in this list object. M
*                      word object will have line number, "TextLineCnt", and M
*                    word count in line, "TextWordCnt", attributes.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if all text fitted in the BoundingPoly, FALSE otherwise.   M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserFontLayoutOverShape, UserFontLayoutOverShape2,                       M
*   UserFontLayoutOverShapeFree, UserFontLayoutOverShapeGenWords.            M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontLayoutOverShapePlaceWords                                        M
*****************************************************************************/
int UserFontLayoutOverShapePlaceWords(UserFontWordLayoutStruct *Words,
				      IrtRType FontSize,
				      const IrtRType FontSpace[3],
				      UserFontAlignmentType Alignment,
				      const UserFontDimInfoStruct *FontDims,
				      const IPPolygonStruct *BoundingPoly,
				      IPObjectStruct **PlacedTextGeom)
{
    int i, n, wc;
    IrtRType y, AllInters[USER_FONT_LOUT_MAX_INTERS],
	FontLineHeight = (FontDims -> AscentLineHeight -
			  FontDims -> DescentLineHeight) * FontSize,
        LineHeight = FontLineHeight + FontSpace[2] / 70.0,
	SpaceWidth = FontDims -> SpaceWidth * FontSize + FontSpace[1] / 70.0;
    IrtPtType LinePt;
    GMBBBboxStruct BBox;
    UserFontWordLayoutStruct *W1, *W2, *W,
        *CrntWord = Words;
    IrtHmgnMatType ScaleMat;

    *PlacedTextGeom = NULL;
    BBox = *GMBBComputeOnePolyBbox(BoundingPoly);
    LinePt[0] = BBox.Min[0] - 1.0;
    LinePt[2] = 0.0;

    for (W = Words; W != NULL; W = W -> Pnext) {
        W -> LeftOverSpace = IRIT_INFNTY;      /* Reset the left over gaps. */
	W -> X = W -> Y = 0;
    }

    /* Generate the lines' level, top to bottom, and place the words. */
    for (y = BBox.Max[1] - FontLineHeight;
	 y >= BBox.Min[1] && CrntWord != NULL;
	 y -= LineHeight) {
        int i, n;

	LinePt[1] = y + FontLineHeight * 0.5;
	AllInters[0] = USER_FONT_LOUT_MAX_INTERS;
	if ((n = UserFontBoundaryLineInter(BoundingPoly, LinePt,
				                             AllInters)) > 0) {
	    assert(n % 2 == 0);		       /* Better be an even number. */

	    /* Sort all X intersection locations from min. to max. values. */
	    qsort(AllInters, n, sizeof(IrtRType), UserFontCmpReals);

	    /* Spread as many words as possible along the valid portions of */
	    /* the line - step over pairs of intersection locations.	    */
	    for (i = 0; i < n - 1; i += 2) {
	        IrtRType
		    X1 = AllInters[i],
		    X2 = AllInters[i + 1];
		UserFontWordLayoutStruct
		    *PrevWord = NULL;

		do {
		    IrtRType
		        WordWidth = (CrntWord -> BBox.Max[0] -
				     CrntWord -> BBox.Min[0]) * FontSize;

		    if (WordWidth < X2 - X1) {            /* Word fit here. */
		        CrntWord -> X = X1;
		        CrntWord -> Y = y;
			/* Advance X1 to end of word + words' space. */
			X1 = X1 + WordWidth + SpaceWidth;
			PrevWord = CrntWord;
			CrntWord = CrntWord -> Pnext;
		    }
		    else {
		        /* Cannot fit next word in (rest of) line. If there */
		        /* is a previous word (last in this line segment),  */
		        /* store the leftover gap in that word.             */
		        if (PrevWord != NULL)
			    PrevWord -> LeftOverSpace = X2 - (X1 - SpaceWidth);
		        break;
		    }

		    /* Force a new line here, if so desired: */
		    if (PrevWord != NULL && PrevWord -> NewLine)
			break;
		}
		while (X1 < X2 && CrntWord != NULL);

		if (CrntWord == NULL) {
		    PrevWord -> LeftOverSpace = X2 - (X1 - SpaceWidth);
		    break;                          /* We placed all words. */
		}
	    }	    
	}
    }

    /* Place the left over words, if any, below the last line. */
    for (W2 = NULL, W1 = CrntWord; W1 != NULL; W2 = W1, W1 = W1 -> Pnext) {
        IrtRType
	    WordWidth = W2 == NULL ? 0.0
				   : (W2 -> BBox.Max[0] - W2 -> BBox.Min[0])
								   * FontSize;

        W1 -> Y = BBox.Min[1] - LineHeight;
	W1 -> X = W2 == NULL ? BBox.Min[0] : W2 -> X + WordWidth + SpaceWidth;
    }

    /* The above placement is a left-alignment placement.  Need to          */
    /* compensate now if other alignment type is desired.		    */
    for (W1 = Words; W1 != NULL && W1 != CrntWord; ) {
        /* Isolate a chain of words in this line segment, between W1 and W2.*/
	for (W = W2 = W1;
	     W != NULL && W != CrntWord && W -> Y == W1 -> Y;
	     W2 = W, W = W -> Pnext) {
	    if (W -> LeftOverSpace != IRIT_INFNTY) {
		W2 = W;
		break;
	    }
	}

	if (W2 != NULL) {
	    if (W2 -> LeftOverSpace == IRIT_INFNTY) {
	        /* its is the last word in sequence. */
		W2 -> LeftOverSpace = 0;
	    }

	    switch (Alignment) {
	        default:
		    assert(0);
	        case USER_FONT_ALIGN_LEFT:
		    /* Nothing to do for this one. */
		    break;
	        case USER_FONT_ALIGN_CENTER:
	            /* Move right all words so gap on left & right is equal.*/
		    for (W = W1; W != W2 -> Pnext; W = W -> Pnext) {
		        W -> X += W2 -> LeftOverSpace * 0.5;
		    }
		    break;
	        case USER_FONT_ALIGN_RIGHT:
	            /* Move right all words so the gap on the right is zero.*/
		    for (W = W1; W != W2 -> Pnext; W = W -> Pnext) {
		        W -> X += W2 -> LeftOverSpace;
		    }
		    break;
	        case USER_FONT_ALIGN_WIDE:
	            /* Space words so no gap on left & right boundary. */
		    for (n = 1, W = W1; W != W2; n++, W = W -> Pnext);

		    if (n > 1) {
		        for (i = 0, W = W1;
			     W != W2 -> Pnext;
			     i++, W = W -> Pnext) {
			    W -> X += W2 -> LeftOverSpace * i / (n - 1.0);
			}
		    }
		    break;
	    }
	}

	W1 = W2 -> Pnext;
    }

    /* Generate the geometry at the newly placed positions. */
    MatGenMatUnifScale(FontSize, ScaleMat);
    *PlacedTextGeom = IPGenLISTObject(NULL);
    n = wc = 0; /* Line number and word count in line */
    for (i = 0, W = Words; W != NULL; W = W -> Pnext) {
        IrtHmgnMatType Mat;
	IPObjectStruct *PObj;

	MatGenMatTrans(W -> X, W -> Y, 0.0, Mat);
	MatMultTwo4by4(Mat, ScaleMat, Mat);

        IPListObjectInsert(*PlacedTextGeom, i++, 
			   PObj = GMTransformObject(W -> Geom, Mat));

	/* Update word count line line and line number attributes. */
	AttrSetObjectIntAttrib(PObj, "TextWordCnt", wc);
	AttrSetObjectIntAttrib(PObj, "TextLineCnt", n);
	if (W -> Pnext != NULL && !IRIT_APX_EQ(W -> Pnext -> Y, W -> Y)) {
	    n++;
	    wc = 0;
	}
	else
	    wc++;
    }
    IPListObjectInsert(*PlacedTextGeom, i++, NULL);      /* Close the list.*/

    return CrntWord == NULL; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes all the intersections of BoundingPoly polygons with a           *
* horizontal line through LinePt.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   BoundingPoly: Polygon to examine for intersections with horizontal line. *
*   LinePt:      A point to the left of BoundingPoly at the Y level of line. *
*   InterXVals:  X values of intersection locations.  A vector where the     *
*                first element must hold on entry its size.                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:                                                                     *
*****************************************************************************/
static int UserFontBoundaryLineInter(const IPPolygonStruct *BoundingPoly,
				     const IrtPtType LinePt,
				     IrtRType *InterXVals)
{
    int n;
    IrtRType FirstInterP;
    IPVertexStruct *FirstInterV;

    n = GMPolygonRayInter2(BoundingPoly, LinePt, 0,
			   &FirstInterV, &FirstInterP, InterXVals);

    return n;    
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Routine to compare two real numbers for sorting purposes.                *
*                                                                            *
* PARAMETERS:                                                                *
*   PReal1, PReal2:  Two pointers to real numbers.                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:   >0, 0, or <0 as the relation between the two reals.               *
*****************************************************************************/
#if defined(ultrix) && defined(mips)
static int UserFontCmpReals(VoidPtr PReal1, VoidPtr PReal2)
#else
static int UserFontCmpReals(const VoidPtr PReal1, const VoidPtr PReal2)
#endif /* ultrix && mips (no const support) */
{
    CagdRType
	Diff = (*((CagdRType *) PReal1)) - (*((CagdRType *) PReal2));

    return IRIT_SIGN(Diff);
}
