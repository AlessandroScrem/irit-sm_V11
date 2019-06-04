/******************************************************************************
* wfnt2bzr.c - extraction util. for windows outline font to bezier curves.    *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Masha Nikolski and Gershon Elber, January 2008.		      *
******************************************************************************/

#ifndef _FREETYPE_FONTS_

#include <ctype.h>

#define UNICODE          /* Support for wide characters (i.e. hebrew etc.). */

#ifdef __WINNT__
#include "windows.h"
#else
#error This file is only to be compiled under windows.
#endif /* __WINNT__ */

#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "user_loc.h"

#define USER_FONT_NUM_SYMBOLS	65536
#define USER_FONT_MAX_ORDER	4 
#define USER_FONT_SCALE		0.001

/* One curve of a symbol. */ 
typedef struct UserFontSymbolCurveStruct {
    IrtPtType Pt; 
    struct UserFontSymbolCurveStruct *Next;
} UserFontSymbolCurveStruct;

/* All curves of a symbol. */ 
typedef struct UserFontSymbolCurveListStruct {
    UserFontSymbolCurveStruct *Crv; 
    struct UserFontSymbolCurveListStruct *Next;
} UserFontSymbolCurveListStruct;

/* One symbol. */ 
typedef struct UserFontSymbolStruct {
#ifdef __WINNT__
    /* Windows measures. */
    GLYPHMETRICS GlyphMeasurements; 
    LPVOID Buffer;
    DWORD  BufferSize;
#endif

    /* Curves of the symbol. */
    UserFontSymbolCurveListStruct *CurvesList;
        
    /* Symbol width. */ 
    IrtRType Width;

    /* Symbol numeric identifier. */
    int Id;

    /* Do not extract the same symbol more than once, */
    /* if it appears numerous times in a string.      */ 
    CagdBType WasUpdated;

    /* Do not free the same symbol more than once, */ 
    /* if it appears numerous times in a string.   */ 
    CagdBType WasFreed;
} UserFontSymbolStruct;

/* Font structure for the current string. */
typedef struct UserFontStruct {
    /* Font transformation matrix. */
#ifdef __WINNT__
    MAT2 TransformMat;
#endif
    /* Text string (UNICODE). */
    wchar_t *TextString;
    /* Space width in the current font. */
    int SpaceWidth;
    /* Bounding box of the string. */
    IrtRType BBox[4];
    /* Amount of shift for each character in the string from the beginning. */
    IrtRType* Shifts;
    /* Information for each symbol in the current font (relevant for        */ 
    /* symbols that appear in the current string only).			    */
    UserFontSymbolStruct Symbols[USER_FONT_NUM_SYMBOLS];
} UserFontStruct;

static void UserFontSymbolToBzr(HDC *Dc,
				HFONT *Font,
				UserFontSymbolStruct *Symbol,
				MAT2 *TransformMat);
static IPObjectStruct *UserFontFont2Bzr(const UserFontText TextString,
					HDC *hdc,
					HFONT *fnt,
					int SpaceWidth,
					int MergeToBsp,
					const char *RootObjName);
static void UserFontExtractGeometry(UserFontSymbolStruct *Symbol); 
static HFONT UserFontCreateFontHandle(const UserFontName FontName, 
				      UserFontStyleType FontStyle,
				      int Height);
static double UserFontDoubleFromFixed(const FIXED f);
static FIXED  UserFontFixedFromDouble(const double d);
static UserFontSymbolCurveListStruct *UserFontAllocSymbolCurveList();
static UserFontSymbolCurveStruct *UserFontAllocSymbolCurve();
static void UserFontFreeSymbolCurveList(UserFontSymbolCurveListStruct *CrvList);
static void UserFontFreeSymbolCurve(UserFontSymbolCurveStruct *Crv);
static void UserFontFreeSymbol(UserFontSymbolStruct *Symb);
static void UserFontComputeBBox(UserFontStruct *Fnt);
static void UserFontGetShiftFromStartString(UserFontStruct *Fnt);
static int  UserFontGetRightDirectionLimits(UserFontStruct *Fnt,
					    int *VecBegins,
					    int *VecEnds);
static IPObjectStruct *UserFontGetTextCrvsAsIritCrvs(UserFontStruct *Fnt,
						     int MergeToBsp,
						     const char *RootObjName);
static CagdCrvStruct *UserFontGetSymbolCrvAsBezierCrv(UserFontSymbolCurveStruct
						                        *Symb,
						      IrtRType XTrans,
						      IrtRType YTrans,
						      IrtRType XScale,
						      IrtRType YScale);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Extracts Bezier curves from windows' fonts.      			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Text:         Text string.  					     M
*   FontName:     Font name.  						     M
*   FontStyle:    Font style.  						     M
*   SpaceWidth:   Space width.  Space, in points, added to individual chars. M
*   MergeToBsp:   TRUE to merge Bezier curves into larger B-spline curves,   M
*                 FALSE to leave the original (font) Bezier curves.	     M
*   RootObjName:  Name of root object.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *: Extracted text object or NULL in case of failure.      M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserFontConvertFontToBezier 			                     M
*****************************************************************************/
IPObjectStruct *UserFontConvertFontToBezier(const UserFontText Text,
					    const UserFontName FontName,
					    UserFontStyleType FontStyle,
					    IrtRType SpaceWidth,
					    int MergeToBsp,
					    const char *RootObjName)
{
    HDC hDC;
    HBITMAP memBM;
    HFONT hFont;
    HGDIOBJ oldFont;
    UINT bfsz;
    OUTLINETEXTMETRIC *MetrixBuffer;
    int lfHeight;
    IPObjectStruct *TextObj;

    hDC = CreateCompatibleDC(NULL);
    memBM = CreateCompatibleBitmap(hDC, 4, 4);
    SelectObject(hDC, memBM);    

    hFont = UserFontCreateFontHandle(FontName, FontStyle, 12);   
    oldFont = SelectObject(hDC, hFont);
    bfsz = GetOutlineTextMetrics(hDC, 0, NULL);

    MetrixBuffer = (OUTLINETEXTMETRIC *) IritMalloc(bfsz);
    GetOutlineTextMetrics(hDC, bfsz, MetrixBuffer);

    SelectObject(hDC, oldFont);
    DeleteObject(hFont);

    lfHeight = MetrixBuffer -> otmEMSquare;
    hFont = UserFontCreateFontHandle(FontName, FontStyle, lfHeight);

    if (hFont == NULL)         
	return NULL;   

    oldFont = SelectObject(hDC, hFont);

    if (&hFont == NULL)
	return NULL;

    TextObj = UserFontFont2Bzr(Text, &hDC, &hFont, (int) SpaceWidth,
			       MergeToBsp, RootObjName);

    SelectObject(hDC, oldFont);
    DeleteObject(hFont);
    DeleteObject(memBM);
    DeleteObject(hDC);

    IritFree(MetrixBuffer);

    return TextObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Extract the Bezier curves of a symbol.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Dc:     Handle to a device context.					     *
*   Font:   Windows font struct.				             *
*   Symbol: Pointer to symbol information struct.                            *
*   TransformMat: Font transformation matrix.                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void UserFontSymbolToBzr(HDC *Dc,
				HFONT *Font,
				UserFontSymbolStruct *Symbol,
				MAT2 *TransformMat)
{
    TEXTMETRIC tm;
    HFONT hOldFont;
    DWORD dwResult;

    if (Symbol -> WasUpdated)
	return;

    Symbol -> WasUpdated = TRUE;
    Symbol -> WasFreed = FALSE;

    /* A space character. */
    if (Symbol -> Id == 32) {	
	GetTextMetrics(*Dc, &tm);        
	Symbol -> Width = tm.tmAveCharWidth;
	return;
    }

    Symbol -> BufferSize = 0;
    Symbol -> Buffer = NULL;

    Symbol -> Width = 0;
    Symbol -> CurvesList = NULL;

    memset(&Symbol -> GlyphMeasurements, 0,
	   sizeof(Symbol -> GlyphMeasurements));

    hOldFont = (HFONT) SelectObject((*Dc), (*Font));

    /* Get size of the buffer needed for this glyph. */
    /* GetOutlineTextMetrics. */
    Symbol -> BufferSize = GetGlyphOutline(*Dc, Symbol -> Id,
					   GGO_NATIVE,
					   &Symbol -> GlyphMeasurements, 0,
					   NULL, TransformMat);
    
    if (Symbol -> BufferSize == GDI_ERROR) {
	SelectObject(*Dc, hOldFont);
	return;
    }

    /* Allocate the buffer. */
    Symbol -> Buffer = IritMalloc(Symbol -> BufferSize);
    if (Symbol -> Buffer  == NULL) {
	SelectObject(*Dc, hOldFont);
	return;
    }

    /* Get vector points for this glyph. */
    dwResult = GetGlyphOutline(*Dc, Symbol -> Id, 
			       GGO_NATIVE, &Symbol -> GlyphMeasurements,
			       Symbol -> BufferSize, Symbol -> Buffer, 
			       TransformMat);

    if (dwResult == GDI_ERROR) {
	SelectObject(*Dc, hOldFont);
	return;
    } 

    UserFontExtractGeometry(Symbol);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Extracts the geometry of a symbol.	        			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Symbol: Pointer to symbol information struct.                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserFontExtractGeometry                        		             *
*****************************************************************************/
static void UserFontExtractGeometry(UserFontSymbolStruct *Symbol)
{
    LPTTPOLYGONHEADER CurrPoly;      /* Pointer to current polygon in glyph. */
    LPTTPOLYCURVE CurrCrv;      /* Pointer to current poly curve in polygon. */
    DWORD
	/* Offset of current poly-curve header from start of polygon header. */
	CurrCrvOffset,
	/* Size of current poly-curve (depends on # of pts in poly curve).   */
	CrvSize, 
	/* Offset of current polygon header from start of buf. */
	CurrPolyOffset = 0;
    int i,
	ContourId = 0;
    IrtRType x, y, x1, y1, x2, y2;    
    IrtPtType PrevPt, FirstPt, p1, p2, p3, p4;        
    UserFontSymbolCurveListStruct *CrvList;

    if (Symbol -> Buffer == NULL)
	return;

    Symbol -> Width = Symbol -> GlyphMeasurements.gmCellIncX;

    /* Outer while loops for all polygon headers (can be more than one). */
    while (Symbol -> BufferSize >= CurrPolyOffset + sizeof(TTPOLYGONHEADER)) {
	/* Get pointer to start of the polygon. */
	CurrPoly = (LPTTPOLYGONHEADER)
	                       (((char *) Symbol -> Buffer) + CurrPolyOffset);
	ContourId++;

	x = UserFontDoubleFromFixed(CurrPoly -> pfxStart.x);
	y = UserFontDoubleFromFixed(CurrPoly -> pfxStart.y);
	PrevPt[0] = x;
	PrevPt[1] = y;

	IRIT_PT_COPY(FirstPt, PrevPt);

	/* Inner while loops for all poly-curves in one polygon.            */
	/* (A poly-curve is one or more polyline and/or QSpline records.)   */
	CurrCrvOffset = sizeof(TTPOLYGONHEADER);
	while (CurrPoly -> cb >= (CurrCrvOffset + sizeof(TTPOLYCURVE))) {
	    /* Get pointer to start of poly-curve. */
	    CurrCrv = (LPTTPOLYCURVE)
	                   (((char  *) Symbol->Buffer) + CurrPolyOffset
							     + CurrCrvOffset);

	    /* Test record type, draw polyline or series of Beziers         */
	    /* accordingly.						    */
	    switch (CurrCrv -> wType) {
		case TT_PRIM_LINE:
		    /* Draw polyline connecting pts. */
		    for (i = 0; i < CurrCrv -> cpfx; i++) {

			CrvList = UserFontAllocSymbolCurveList();
			CrvList -> Next = Symbol -> CurvesList;
			Symbol -> CurvesList = CrvList;
			
			IRIT_PT_COPY(CrvList -> Crv -> Pt, PrevPt);
			
			x = UserFontDoubleFromFixed(CurrCrv -> apfx[i].x);
			y = UserFontDoubleFromFixed(CurrCrv -> apfx[i].y);

			PrevPt[0] = x;
			PrevPt[1] = y;

			CrvList -> Crv -> Next = UserFontAllocSymbolCurve();
			IRIT_PT_COPY(CrvList -> Crv -> Next -> Pt, PrevPt);
		    }
		    break;
		case TT_PRIM_QSPLINE:
		{
		    /* Draw series of Beziers connecting points. But for    */
		    /* initial Bezier, grab last point on previous curve.   */

		    /* Three points defining spline. */		   
		    IRIT_PT_COPY(p3, PrevPt);
		    /* Draw QSplines. */
		    for (i = 0; i < CurrCrv -> cpfx - 1; i++) {
			/* p1 is 1st control -- last point in this or       */
			/* previous contour.				    */
			IRIT_PT_COPY(p1, p3);

			/* p2 is handle -- a point in the record. */
			x = UserFontDoubleFromFixed(CurrCrv -> apfx[i].x),
			y = UserFontDoubleFromFixed(CurrCrv -> apfx[i].y);

			p2[0] = x;
			p2[1] = y;
			
			/* p3 is 2nd control -- the midpoint of 2 spline     */
			/* points except for the last spline, when it is the */
			/* second to the last point in the spline record.    */
			 x1 = UserFontDoubleFromFixed(CurrCrv -> apfx[i+1].x),
			 y1 = UserFontDoubleFromFixed(CurrCrv -> apfx[i+1].y);

			if (i == (CurrCrv -> cpfx-2)) {
			    p3[0] = x1;
			    p3[1] = y1;
			}
			else {
			    p3[0] = 0.5 * (x1 + x);
			    p3[1] = 0.5 * (y1 + y);
			}
			CrvList = UserFontAllocSymbolCurveList();
			CrvList -> Next = Symbol -> CurvesList;
			Symbol -> CurvesList = CrvList;

			IRIT_PT_COPY(CrvList -> Crv->Pt, p1);
			CrvList -> Crv -> Next = UserFontAllocSymbolCurve();
			IRIT_PT_COPY(CrvList -> Crv->Next -> Pt, p2);
			CrvList -> Crv -> Next -> Next =
			                          UserFontAllocSymbolCurve();
			IRIT_PT_COPY(CrvList -> Crv -> Next -> Next -> Pt, p3);
			IRIT_PT_COPY(PrevPt, p3);
		    } 
		    break; 
		}
		default:
		{
		    IRIT_PT_COPY(p4, PrevPt);
		    /* Draw QSplines. */
		    for (i = 0; i < CurrCrv -> cpfx - 2; i++) {
			/* p1 is 1st control -- last point in this or       */
			/* previous contour.				    */
			IRIT_PT_COPY(p1, p4);
			/* p2 is handle -- a point in the record. */
			x = UserFontDoubleFromFixed(CurrCrv -> apfx[i].x),
			y = UserFontDoubleFromFixed(CurrCrv -> apfx[i].y);

			p2[0] = x;
			p2[1] = y;
			
			/* p3 is 2nd control -- the midpoint of 2 spline     */
			/* points except for the last spline, when it is the */
			/* second to the last point in the spline record.    */

			x1 = UserFontDoubleFromFixed(CurrCrv -> apfx[i+1].x),
			y1 = UserFontDoubleFromFixed(CurrCrv -> apfx[i+1].y);

			p3[0] = x1;			
			p3[1] = y1;
	
			x2 = UserFontDoubleFromFixed(CurrCrv -> apfx[i+2].x),
			y2 = UserFontDoubleFromFixed(CurrCrv -> apfx[i+2].y);

			p4[0] = x2;
			p4[1] = y2;

			CrvList = UserFontAllocSymbolCurveList();
			CrvList -> Next = Symbol -> CurvesList;
			Symbol -> CurvesList = CrvList;

			IRIT_PT_COPY(CrvList -> Crv -> Pt, p1);
			CrvList -> Crv -> Next = UserFontAllocSymbolCurve();
			IRIT_PT_COPY(CrvList -> Crv -> Next -> Pt, p2);
			CrvList -> Crv -> Next -> Next =
			                          UserFontAllocSymbolCurve();
			IRIT_PT_COPY(CrvList -> Crv -> Next -> Next -> Pt, p3);
			CrvList -> Crv -> Next -> Next -> Next =
			                          UserFontAllocSymbolCurve();
			IRIT_PT_COPY(CrvList -> Crv -> Next -> Next -> Next -> Pt,
				     p4);
			
			IRIT_PT_COPY(PrevPt, p4);
			i++;
		    } 
		    break; 
		}
	    } 

	    /* Increment curve offset so point to next poly-curve. */
	    CrvSize = sizeof(TTPOLYCURVE) + ((CurrCrv -> cpfx - 1)
		                                           * sizeof(POINTFX));
	    CurrCrvOffset += CrvSize;
	} 

	/* Ignore the degenerate case. */
	if (!(FirstPt[0] == PrevPt[0] && FirstPt[1] == PrevPt[1])) {
	    CrvList = UserFontAllocSymbolCurveList();
	    CrvList -> Next = Symbol -> CurvesList;
	    Symbol -> CurvesList = CrvList;

	    IRIT_PT_COPY(CrvList -> Crv -> Pt, PrevPt);
	    CrvList -> Crv -> Next = UserFontAllocSymbolCurve();
	    IRIT_PT_COPY(CrvList -> Crv -> Next -> Pt, FirstPt);
	}

	/* Increment polygon offset so point to next polygon. */
	CurrPolyOffset += CurrPoly -> cb;
    }	
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Convert from a fixed-point real number to double.                        *
*                                                                            *
* PARAMETERS:                                                                *
*   f:   A fixed-point real number struct.                                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   double: converted value.            				     *
*****************************************************************************/
static double UserFontDoubleFromFixed(const FIXED f)
{
    return (*((const long *) (&f))) / 65536.0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Convert from double to a fixed-point real number.                        *
*                                                                            *
* PARAMETERS:                                                                *
*   d: Double value.                      				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   FIXED: a fixed-point real number struct.                                 *
*****************************************************************************/
static FIXED UserFontFixedFromDouble(const double d)
{
    long
	l = (long) (d * 65536.0);

    return *(FIXED *) &l;
} 

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Allocation of a UserFontSymbolCurveListStruct structure.                 *
*                                                                            *
* PARAMETERS:                                                                *
*   None	                          				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   UserFontSymbolCurveListStruct *: Allocated memory.                       *
*****************************************************************************/
static UserFontSymbolCurveListStruct *UserFontAllocSymbolCurveList()
{
    UserFontSymbolCurveListStruct
        *NewList = (UserFontSymbolCurveListStruct *) 
                            IritMalloc(sizeof(UserFontSymbolCurveListStruct));
    NewList -> Crv = UserFontAllocSymbolCurve();
    NewList -> Next = NULL;

    return NewList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Allocation of a UserFontSymbolCurveStruct structure.                     *
*                                                                            *
* PARAMETERS:                                                                *
*   None	                          				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   UserFontSymbolCurveStruct *: Allocated memory.                           *
*****************************************************************************/
static UserFontSymbolCurveStruct *UserFontAllocSymbolCurve()
{
    UserFontSymbolCurveStruct
        *NewCurve = NULL;

    NewCurve = (UserFontSymbolCurveStruct *)
                                IritMalloc(sizeof(UserFontSymbolCurveStruct));
    NewCurve -> Next = NULL;

    return NewCurve;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Extract the Bezier curves of a font.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   TextString: String to perform the extraction for.			     *
*   hdc:        Handle to device context.				     *
*   fnt:        Windows font struct.					     *
*   SpaceWidth: Space width. If USER_FONT_DEFAULT_WIDTH then default space   *
*               is used.					             *
*   MergeToBsp: TRUE to merge the Bezier curves into larger B-spline curves, *
*               FALSE to leave the original (font) Bezier curves.	     *
*   RootObjName:  Name of root object.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct*: Text string as IPObjectStruct.                          *
*****************************************************************************/
static IPObjectStruct *UserFontFont2Bzr(const UserFontText TextString,
					HDC *hdc,
					HFONT *fnt,
					int SpaceWidth,
					int MergeToBsp,
					const char *RootObjName)
{
    TEXTMETRIC tm;
    UserFontStruct *NewFont;
    IPObjectStruct *TextObj;
    int NumberOfSymbols, CharVal, i;

    NewFont = (UserFontStruct *) IritMalloc(sizeof(UserFontStruct));

    for (i = 0; i < USER_FONT_NUM_SYMBOLS; i++) {
	NewFont -> Symbols[i].Id = i;
	NewFont -> Symbols[i].WasUpdated = FALSE;
    }

    GetTextMetrics(*hdc, &tm);        
    NewFont -> SpaceWidth = SpaceWidth;
    
    NumberOfSymbols = (int) wcslen(TextString) + 1;
    NewFont -> TextString = (wchar_t *)
                                IritMalloc(sizeof(wchar_t) * NumberOfSymbols);
    wcscpy(NewFont -> TextString, TextString);

    NewFont -> Shifts = NULL;
    NewFont -> TransformMat.eM11 =
        NewFont -> TransformMat.eM22 = UserFontFixedFromDouble(1.0);
    NewFont -> TransformMat.eM12 =
        NewFont -> TransformMat.eM21 = UserFontFixedFromDouble(0.0);
   
    for (i = 0; i < NumberOfSymbols - 1; ++i) {
	CharVal = (int) TextString[i];	
	if (!NewFont -> Symbols[CharVal].WasUpdated)
	    UserFontSymbolToBzr(hdc, fnt, &NewFont -> Symbols[CharVal],
			        &NewFont -> TransformMat);    	    
    }    

    UserFontComputeBBox(NewFont);
    TextObj = UserFontGetTextCrvsAsIritCrvs(NewFont, MergeToBsp, RootObjName);

    for (i = 0; i < NumberOfSymbols - 1; ++i) {
	CharVal = (int) TextString[i];
	UserFontFreeSymbol(&NewFont -> Symbols[CharVal]);	
    }    

    IritFree(NewFont -> TextString);
    IritFree(NewFont -> Shifts);
    IritFree(NewFont);

    return TextObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates a handle to extracted font.           			     *
*                                                                            *
* PARAMETERS:                                                                *
*   FontName:  Font name.  						     *
*   FontStyle: Font style.  						     *
*   Height:    Font height.  						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   HFONT:  A handle to extracted font.                                      *
*****************************************************************************/
static HFONT UserFontCreateFontHandle(const UserFontName FontName,
				      UserFontStyleType FontStyle,
				      int Height)
{
    HFONT hFont;
    WCHAR 
	*UniFontName = IritMalloc((strlen(FontName) + 1) * sizeof(WCHAR));

    MultiByteToWideChar(CP_OEMCP, 0, FontName, -1, UniFontName,
		        (int) (strlen(FontName) + 1));

    switch (FontStyle) {
	case USER_FONT_STYLE_REGULAR:
	    hFont = CreateFont(Height, 0, 10, 10, FW_REGULAR,
			       FALSE, FALSE, FALSE, 
			       DEFAULT_CHARSET, OUT_OUTLINE_PRECIS,
			       CLIP_DEFAULT_PRECIS,
			       DEFAULT_QUALITY, FF_ROMAN,
			       UniFontName);
	    break;
	case USER_FONT_STYLE_BOLD:
	    hFont = CreateFont(Height, 0, 10, 10, FW_BOLD,
			       FALSE, FALSE, FALSE, 
			       DEFAULT_CHARSET, OUT_OUTLINE_PRECIS,
			       CLIP_DEFAULT_PRECIS,
			       DEFAULT_QUALITY, FF_ROMAN,
			       UniFontName);
	    break;
	case USER_FONT_STYLE_ITALICS:
	    hFont = CreateFont(Height, 0, 10, 10, FW_REGULAR,
			       TRUE, FALSE, FALSE, 
			       DEFAULT_CHARSET, OUT_OUTLINE_PRECIS,
			       CLIP_DEFAULT_PRECIS,
			       DEFAULT_QUALITY, FF_ROMAN,
			       UniFontName);
	    break;
	case USER_FONT_STYLE_BOLD_ITALICS:
	    hFont = CreateFont(Height, 0, 10, 10, FW_BOLD,
			       TRUE, FALSE, FALSE, 
			       DEFAULT_CHARSET, OUT_OUTLINE_PRECIS,
			       CLIP_DEFAULT_PRECIS,
			       DEFAULT_QUALITY, FF_ROMAN,
			       UniFontName);
    }

    IritFree(UniFontName);

    return hFont;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the bbox of a text string. 				     *
*				 		                             *
* PARAMETERS:                                                                *
*   Fnt: Text string font structure.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void UserFontComputeBBox(UserFontStruct *Fnt)
{
    UserFontSymbolCurveListStruct *Crvs;
    UserFontSymbolCurveStruct *SymbolCrv;
    IrtRType Shift;
    int i, Val, Length;

    UserFontGetShiftFromStartString(Fnt);
    Length = (int) wcslen(Fnt -> TextString);
    Fnt -> BBox[0] = Fnt -> BBox[1] = IRIT_MAX_INT;
    Fnt -> BBox[2] = Fnt -> BBox[3] = -IRIT_MAX_INT;

    for (i = (int) Length - 1; i >= 0; i--) {
	Shift = Fnt -> Shifts[i];

	if ((Val = (int) (Fnt -> TextString[i])) == (int) ' ')
	    continue;
	
	Crvs = Fnt -> Symbols[Val].CurvesList;
	
	while (Crvs != NULL) {
	    SymbolCrv = Crvs -> Crv;
	    while (SymbolCrv != NULL) {
		if (Fnt -> BBox[0] > SymbolCrv -> Pt[0] + Shift)
		    Fnt -> BBox[0] = SymbolCrv -> Pt[0] + Shift;
		if (Fnt -> BBox[1] > SymbolCrv -> Pt[1])
		    Fnt -> BBox[1] = SymbolCrv -> Pt[1];
		if (Fnt -> BBox[2] < SymbolCrv -> Pt[0] + Shift)
		    Fnt -> BBox[2] = SymbolCrv -> Pt[0] + Shift;
		if (Fnt -> BBox[3] < SymbolCrv -> Pt[1])
		    Fnt -> BBox[3] = SymbolCrv -> Pt[1];
		SymbolCrv = SymbolCrv -> Next;
	    }
	    Crvs = Crvs -> Next;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute the shift of every symbol in m_strCurrString from the beginning  *
*   of the string.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   Fnt: Text string font structure.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void							 	     *
*****************************************************************************/
static void UserFontGetShiftFromStartString(UserFontStruct *Fnt) 
{
    int i, nCurrentLimit,
	CurrStringSize = (int) wcslen(Fnt -> TextString),
	*vecBegins = (int *) IritMalloc(sizeof(int) * CurrStringSize),
        *vecEnds = (int *) IritMalloc(sizeof(int) * CurrStringSize),
	nLimitsSize = UserFontGetRightDirectionLimits(Fnt, vecBegins, vecEnds);
    IrtRType 
	CurrShift = 0.0;

    Fnt -> Shifts = (IrtRType *) IritMalloc(sizeof(IrtRType) * CurrStringSize);

    if (nLimitsSize > 0)
	nCurrentLimit =	0;
    else {
	for (i = 0 ; i < CurrStringSize; i++) {
	    int Val = (int) (Fnt -> TextString[i]);

	    Fnt -> Shifts[i] = CurrShift;
	    
	    CurrShift += Fnt -> Symbols[Val].Width;		    
	}
	return;
    }

    for (i = 0 ; i < CurrStringSize; i++) {	
	if (nCurrentLimit < nLimitsSize	&& i ==	vecBegins[nCurrentLimit]) {
	    int j;

	    for (j = vecEnds[nCurrentLimit];
		j >= vecBegins[nCurrentLimit]; 
		j--) {
		    int Val2 = (int) (Fnt -> TextString[j]);

		    Fnt -> Shifts[j] = CurrShift;
		    CurrShift += Fnt -> Symbols[Val2].Width;		 
	    }			
	    i = vecEnds[nCurrentLimit++];
	}
	else {
	    int Val = (int) (Fnt -> TextString[i]);

	    Fnt -> Shifts[i] = CurrShift;
	    CurrShift += Fnt -> Symbols[Val].Width;	    
	}					
    }			

    IritFree(vecBegins);
    IritFree(vecEnds);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compute the shift of every word in string from the beginning	     *
*   and the end of the string.          				     *
*                                                                            *
* PARAMETERS:                                                                *
*   Fnt:       Text string font structure.				     *
*   vecBegins: Vector to fill in the shifts of beginnings of words.          *
*   vecEnds:   Vector to fill in the shifts of ends of words.	             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:   Number of words in a string.                	 	             *
*****************************************************************************/
static int UserFontGetRightDirectionLimits(UserFontStruct *Fnt,
					  int *vecBegins,
					  int *vecEnds)
{
    int i,
        nStringSize = (int) wcslen(Fnt -> TextString),
        nLimitsSize = 0,
	nLastIndexB = 0,
	nLastIndexE = 0;
    WORD
        dPrevType = C2_LEFTTORIGHT, 
	*dTypeArray = (WORD *) IritMalloc(sizeof(WORD) * nStringSize);

    GetStringTypeW(CT_CTYPE2, Fnt -> TextString, nStringSize, 
		   dTypeArray);
    
    for (i = 0; i < nStringSize; i++) {
	if (dTypeArray[i] == C2_RIGHTTOLEFT && dPrevType != dTypeArray[i]) {
	    vecBegins[nLastIndexB] = i;
	    nLastIndexB++;
	}
	if (dTypeArray[i] == C2_LEFTTORIGHT && dPrevType == C2_RIGHTTOLEFT) {
	    vecEnds[nLastIndexE] = i - 1;
	    nLastIndexE++;
	}
	if (dTypeArray[i] == C2_LEFTTORIGHT || dTypeArray[i] ==	C2_RIGHTTOLEFT)
	    dPrevType = dTypeArray[i];
    }

    if (dPrevType == C2_RIGHTTOLEFT)
	vecEnds[nLastIndexE] = nStringSize - 1;

    IritFree(dTypeArray);

    return (int) nLastIndexB;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Return an IPObjectStruct object of text string.             	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Fnt:	 Text string font structure.				     *
*   MergeToBsp:  TRUE to merge the Bezier curves into larger B-spline curves,*
*                FALSE to leave the original (font) Bezier curves.	     *
*   RootObjName: Name of root object.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *: Text string object.                         	     *
*****************************************************************************/
static IPObjectStruct *UserFontGetTextCrvsAsIritCrvs(UserFontStruct *Fnt,
						     int MergeToBsp,
						     const char *RootObjName)
{
    IPObjectStruct
        *Text = IPGenListObject(RootObjName, NULL, NULL);
    UserFontSymbolCurveListStruct *SymbCrvList;
    IPObjectStruct *Letter, *TmpObj, *Symbol;
    IrtRType
	*BBox = Fnt -> BBox;
    int Len = (int) wcslen(Fnt -> TextString),
	j = 0,
	i = 0,
	NumLetters = 0;

    for (i = 0; i < Len; i++) {
	if (Fnt -> TextString[i] != ' ') {
	    char Name[IRIT_LINE_LEN_LONG];

	    Letter = IPGenLISTObject(NULL);
	    if (iswalnum(Fnt -> TextString[i]))
		sprintf(Name, "%s_%c%d", RootObjName, Fnt -> TextString[i], i);
	    else
		sprintf(Name, "%s_%d_%d", RootObjName, Fnt -> TextString[i], i);
	    IP_SET_OBJ_NAME2(Letter, Name);

	    AttrSetObjectIntAttrib(Letter, "letter", Fnt -> TextString[i]);

	    SymbCrvList =
		Fnt -> Symbols[(int) Fnt -> TextString[i]].CurvesList;

	    j = 0;
	    while (SymbCrvList != NULL) {
  		Symbol = IPGenCRVObject(UserFontGetSymbolCrvAsBezierCrv(
		                                    SymbCrvList -> Crv,
						    Fnt -> Shifts[i] +
						        i * Fnt -> SpaceWidth,
						    0.0, USER_FONT_SCALE,
						    USER_FONT_SCALE));

		if (iswalnum(Fnt -> TextString[i]))
		    sprintf(Name, "%s_%c%d_%d",
			    RootObjName, Fnt -> TextString[i], i, j);
		else
		    sprintf(Name, "%s_%d_%d_%d",
			    RootObjName, Fnt -> TextString[i], i, j);
		IP_SET_OBJ_NAME2(Symbol, Name);
		IPListObjectInsert(Letter, j++, Symbol);

		SymbCrvList = SymbCrvList -> Next;
	    }
	    /* Reverse the list in place. */
	    IPListObjectInsert(Letter, j, NULL);
	    IPReverseListObj(Letter);     /* Reverse curves order in place. */

	    if (MergeToBsp) {
	        IrtBType HaveHoles;
	        CagdCrvStruct *Crv,
		    *Crvs = UserFontBzrList2BspList(Letter, &HaveHoles);

		if (Crvs -> Pnext != NULL) {
		    TmpObj = IPGenListObject(IP_GET_OBJ_NAME(Letter),
					     NULL, NULL);
		    IPFreeObject(Letter);
		    Letter = TmpObj;

		    j = 0;
		    while (Crvs != NULL) {
		        IRIT_LIST_POP(Crv, Crvs);
			TmpObj = IPGenCrvObject(IP_GET_OBJ_NAME(Letter), Crvs,
						NULL);

			if (iswalnum(Fnt -> TextString[i]))
			    sprintf(Name, "%s_%c%d_%d",
				    RootObjName, Fnt -> TextString[i], i, j);
			else
			    sprintf(Name, "%s_%d_%d_%d",
				    RootObjName, Fnt -> TextString[i], i, j);
			Symbol = IPGenCrvObject(Name, Crv, NULL);
			IPListObjectInsert(Letter, j++, Symbol);	   
		    }
		    IPListObjectInsert(Letter, j, NULL);
		}
		else {
		    TmpObj = IPGenCrvObject(IP_GET_OBJ_NAME(Letter),
					    Crvs, NULL);
		    IPFreeObject(Letter);
		    Letter = TmpObj;
		}
	    }

	    IPListObjectInsert(Text, NumLetters, Letter);
	    NumLetters++;
	}
    }

    IPListObjectInsert(Text, NumLetters, NULL);

    return Text;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Return an CagdCrvStruct * object of UserFontSymbolCurveStruct.     	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Symb:   A symbol.               					     *
*   XTrans: Translation amount in X.                            	     *
*   YTrans: Translation amount in Y.                            	     *
*   XScale: Scale amount in X.                                  	     *
*   YScale: Scale amount in Y.                                  	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdCrvStruct *: CagdCrvStruct curve.                         	     *
*****************************************************************************/
static CagdCrvStruct *UserFontGetSymbolCrvAsBezierCrv(UserFontSymbolCurveStruct
						                        *Symb,
						      IrtRType XTrans,
						      IrtRType YTrans,
						      IrtRType XScale,
						      IrtRType YScale)
{    
    int i = 0,
	Order = 0;
    IrtRType **Pts;
    CagdCrvStruct *Crv;
    UserFontSymbolCurveStruct
        *SymbCrv = Symb;

    while (SymbCrv != NULL) {
	SymbCrv = SymbCrv -> Next;
	Order++;
    }
    
    assert(Order >= 2 && Order <= USER_FONT_MAX_ORDER);

    Crv = BzrCrvNew(Order, CAGD_PT_E2_TYPE);
    Pts = Crv -> Points;
    
    SymbCrv = Symb;
    while (SymbCrv != NULL) {
	Pts[1][i] = (SymbCrv -> Pt[0] + XTrans) * XScale;
	Pts[2][i] = (SymbCrv -> Pt[1] + YTrans) * YScale;
	SymbCrv = SymbCrv -> Next;
	i++;
    }

    return Crv;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Frees UserFontSymbolCurveListStruct structure.                    	     *
*                                                                            *
* PARAMETERS:                                                                *
*   CrvList: Memory to free.           					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                        	     *
*****************************************************************************/
static void UserFontFreeSymbolCurveList(UserFontSymbolCurveListStruct *CrvList)
{
    UserFontSymbolCurveListStruct
	*Temp = CrvList;
    
    while (CrvList != NULL) {
	Temp = CrvList -> Next;
	UserFontFreeSymbolCurve(CrvList -> Crv);
	IritFree(CrvList);
	CrvList = Temp;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Frees UserFontSymbolCurve structure.                             	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Crv: Memory to free.           					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                        	     *
*****************************************************************************/
static void UserFontFreeSymbolCurve(UserFontSymbolCurveStruct *Crv)
{
    UserFontSymbolCurveStruct
	*Temp = Crv;

    while (Crv != NULL) {
	Temp = Crv -> Next;
	IritFree(Crv);
	Crv = Temp;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Frees UserFontSymbol structure.                               	     *
*                                                                            *
* PARAMETERS:                                                                *
*   Symb: Memory to free.           					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                        	     *
*****************************************************************************/
static void UserFontFreeSymbol(UserFontSymbolStruct *Symb)
{
    if (Symb -> WasUpdated && !Symb -> WasFreed && Symb -> Id!= 32) {
	UserFontFreeSymbolCurveList(Symb -> CurvesList);
	IritFree(Symb -> Buffer);
	Symb -> WasFreed = TRUE;
    }
}

#endif /* !_FREETYPE_FONTS_ */
