/*****************************************************************************
* Setting attributes for objects.					     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 0.2, Mar. 1990   *
*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "inc_irit/irit_sm.h"
#include "misc_loc.h"
#include "inc_irit/miscattr.h"

#define ATTR_ALLOCATE_NUM 1000

IRIT_STATIC_DATA IPAttributeStruct
    *AttrFreedList = NULL;

#define DEBUG_GET_ATTR_DUMP_FUNC_CALL
/*     fprintf(stderr, "Get Function %s, line %d called\n", \        */
/*             __FUNCTION__, __LINE__);				     */

#define DEBUG_SET_ATTR_DUMP_FUNC_CALL
/*    fprintf(stderr, "Set Function %s, line %d called\n", \         */
/*            __FUNCTION__, __LINE__);				     */

IRIT_STATIC_DATA const char
   **_AttrValidAttrList = NULL;

IRIT_GLOBAL_DATA int AttrIritColorTable[][3] = {
    {   0,    0,    0 },			/*  0. Black.        */
    {   0,    0,  255 },			/*  1. Blue.         */
    {   0,  255,    0 },			/*  2. Green.        */
    {   0,  255,  255 },			/*  3. Cyan.         */
    { 255,    0,    0 },			/*  4. Red.          */
    { 255,    0,  255 },			/*  5. Magenta.      */
    { 128,  128,    0 },			/*  6. Brown.        */
    { 128,  128,  128 },			/*  7. Lightgrey.    */
    {  64,   64,   64 },			/*  8. Darkgray.     */
    {  64,   64,  255 },			/*  9. Lightblue.    */
    {  64,  255,   64 },			/* 10. Lightgreen.   */
    {  64,  255,  255 },			/* 11. Lightcyan.    */
    { 255,   64,   64 },			/* 12. Lightred.     */
    { 255,   64,  255 },			/* 13. Lightmagenta. */
    { 255,  255,   64 },			/* 14. Yellow.       */
    { 255,  255,  255 }				/* 15. White.        */
};

IPAttributeStruct *_AttrMallocNumAttribute(AttribNumType AttribNum, 
					  IPAttributeType Type);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to fetch one of 16 basic colors based on its index.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Color:		Index of color bwteeen 0 and 15.	             M
*   Red, Green, Blue:   Component of RGB color.                              M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrGetRGBColor, AttrSetWidth, AttrGetWidth, M
*   AttrSetIntAttrib, AttrSetObjectRGBColor, AttrSetRGBColor		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetIndexColor, attributes, color, rgb                                M
*****************************************************************************/
void AttrGetIndexColor(int Color, int *Red, int *Green, int *Blue)
{
    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    Color = IRIT_BOUND(Color, 0, 15);
    *Red = AttrIritColorTable[Color][0];
    *Green = AttrIritColorTable[Color][1];
    *Blue = AttrIritColorTable[Color][2];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to set a color attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Where to place the color attribute.                           M
*   Color:     New color.		                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrGetColor, AttrSetRGBColor, AttrGetRGBColor, AttrSetWidth,            M
*   AttrGetWidth, AttrSetIntAttrib, AttrSetObjectColor			     M
*                                                                	     *
* KEYWORDS:                                                                  M
*   AttrSetColor, attributes, color	                                     M
*****************************************************************************/
void AttrSetColor(IPAttributeStruct **Attrs, int Color)
{
    IRIT_STATIC_DATA AttribNumType
        ColorAttribNum = ATTRIB_NAME_BAD_NUMBER;
    IPAttributeStruct *Attr;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL
  
    if (ColorAttribNum == ATTRIB_NAME_BAD_NUMBER)
        ColorAttribNum = AttrGetAttribNumber("color");

    Attr = AttrFindNumAttribute(*Attrs, ColorAttribNum);
    
    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_INT;
    }
    else {
	Attr = _AttrMallocNumAttribute(ColorAttribNum, IP_ATTR_INT);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.I = Color;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to return a color attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     For which we would like to know the color of.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:      Color or IP_ATTR_NO_COLOR if no color set.	             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrSetRGBColor, AttrGetRGBColor, AttrSetWidth,	     M
*   AttrGetWidth, AttrGetIntAttrib, AttrGetObjectColor			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetColor, attributes, color 	                                     M
*****************************************************************************/
int AttrGetColor(const IPAttributeStruct *Attrs)
{
    IRIT_STATIC_DATA AttribNumType
        ColorAttribNum = ATTRIB_NAME_BAD_NUMBER;
    IPAttributeStruct
        *Attr = NULL;

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (ColorAttribNum == ATTRIB_NAME_BAD_NUMBER) {
        if ((ColorAttribNum = AttrGetAttribNumber("color"))
						== ATTRIB_NAME_BAD_NUMBER)
	    return IP_ATTR_NO_COLOR;
    }

    Attr = AttrFindNumAttribute(Attrs, ColorAttribNum);

    if (Attr != NULL) {
	if (Attr -> Type == IP_ATTR_INT)
	    return Attr -> U.I;
	else if (Attr -> Type == IP_ATTR_STR)
	    return atoi(Attr -> U.Str);
	else
	    return IP_ATTR_NO_COLOR;
    }
    else
        return IP_ATTR_NO_COLOR;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to set an RGB color attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:		Where to place the TGB color attribute.              M
*   Red, Green, Blue:   Component of RGB color.                              M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrGetRGBColor, AttrSetWidth, AttrGetWidth, M
*   AttrSetIntAttrib, AttrSetObjectRGBColor				     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetRGBColor, attributes, color, rgb                                  M
*****************************************************************************/
void AttrSetRGBColor(IPAttributeStruct **Attrs, int Red, int Green, int Blue)
{
    IRIT_STATIC_DATA AttribNumType
        RGBAttribNum = ATTRIB_NAME_BAD_NUMBER;
    char SRGB[30];
    IPAttributeStruct *Attr;
    
    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (RGBAttribNum == ATTRIB_NAME_BAD_NUMBER)
        RGBAttribNum = AttrGetAttribNumber("rgb");

    sprintf(SRGB, "%d,%d,%d", Red, Green, Blue);

    Attr = AttrFindNumAttribute(*Attrs, RGBAttribNum);
    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> U.Str = IritStrdup(SRGB);
	Attr -> Type = IP_ATTR_STR;
    }
    else {
	Attr = _AttrMallocNumAttribute(RGBAttribNum, IP_ATTR_STR);
	Attr -> U.Str = IritStrdup(SRGB);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to return a RGB attribute.				     	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     		For which we would like to know the RGB of.          M
*   Red, Green, Blue:   Component of RGB color to initialize.                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if does have an RGB color attribute, FALSE otherwise.      M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrSetRGBColor, AttrSetWidth, AttrGetWidth, M
*   AttrSetIntAttrib, AttrGetObjectRGBColor,				     M
*   AttrSetRGBDoubleColor, AttrGetRGBDoubleColor			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRGBColor, attributes, color, rgb	                             M
*****************************************************************************/
int AttrGetRGBColor(const IPAttributeStruct *Attrs,
		    int *Red,
		    int *Green,
		    int *Blue)
{
    IRIT_STATIC_DATA AttribNumType 
        RGBAttribNum = ATTRIB_NAME_BAD_NUMBER;
    int i;
    IPAttributeStruct
        *p = NULL;

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (RGBAttribNum == ATTRIB_NAME_BAD_NUMBER) {
	if ((RGBAttribNum = AttrGetAttribNumber("rgb"))
						    == ATTRIB_NAME_BAD_NUMBER)
	    return FALSE;
    }
    
    p = AttrFindNumAttribute(Attrs, RGBAttribNum);

    if (NULL == p || p -> Type != IP_ATTR_STR)
	return FALSE;

    if (Red && Green && Blue)
        return sscanf(p -> U.Str, "%d,%d,%d", Red, Green, Blue) == 3;
    else
        return sscanf(p -> U.Str, "%d,%d,%d", &i, &i, &i) == 3;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to return a RGB attribute or COLOR attribute converted to RGB.     M
*   Beside, it can be used to get other color attributes: specular, etc.     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     		For which we would like to know the RGB of.          M
*   Name:     		Name of the attribute, if NULL default is taken.     M
*   Red, Green, Blue:   Component of RGB color to initialize.                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if does have any color attribute, FALSE otherwise.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrSetRGBColor, AttrSetWidth, AttrGetWidth, M
*   AttrSetIntAttrib, AttrGetObjectRGBColor, AttrGetRGBColor,		     M
*   AttrSetRGBDoubleColor, AttrGetRGBDoubleColor			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRGBColor2, attributes, color, rgb	                             M
*****************************************************************************/
int AttrGetRGBColor2(const IPAttributeStruct *Attrs, 
		     const char *Name,
		     int *Red, 
		     int *Green, 
		     int *Blue)
{
    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Name == NULL) {
	int Color;

	if (AttrGetRGBColor(Attrs, Red, Green, Blue))
	    return TRUE;
	Color = AttrGetColor(Attrs);
	if (IP_ATTR_NO_COLOR == Color)
	    return FALSE;
	if (Color < 0 || Color >= sizeof(AttrIritColorTable) /
					       sizeof(AttrIritColorTable[0]))
	    return FALSE;

	if (Red && Green && Blue) {
	    *Red = AttrIritColorTable[Color][0];
	    *Green = AttrIritColorTable[Color][1];
	    *Blue = AttrIritColorTable[Color][2];
	}
    }
    else {
	const char *p;

	if (NULL == (p = AttrGetStrAttrib(Attrs, Name)))
	    return FALSE;
	if (Red && Green && Blue)
	    return sscanf(p, "%d,%d,%d", Red, Green, Blue) == 3;
	return TRUE;
    }
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to return a RGB attribute.				     	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     		For which we would like to know the RGB of.          M
*   Red, Green, Blue:   Component of RGB color to initialize.                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if does have an RGB color attribute, FALSE otherwise.      M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrSetRGBColor, AttrSetWidth, AttrGetWidth, M
*   AttrSetIntAttrib, AttrGetObjectRGBColor, AttrSetRGBDoubleColor	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRGBDoubleColor, attributes, color, rgb                            M
*****************************************************************************/
int AttrGetRGBDoubleColor(const IPAttributeStruct *Attrs,
			  double *Red,
			  double *Green,
			  double *Blue)
{
    IRIT_STATIC_DATA AttribNumType 
        RGBAttribNum = ATTRIB_NAME_BAD_NUMBER;
    double d;
    IPAttributeStruct
        *p = NULL;

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (RGBAttribNum == ATTRIB_NAME_BAD_NUMBER) {
	if ((RGBAttribNum = AttrGetAttribNumber("rgbDouble"))
						    == ATTRIB_NAME_BAD_NUMBER)
	    return FALSE;
    }
    
    p = AttrFindNumAttribute(Attrs, RGBAttribNum);

    if (NULL == p || p -> Type != IP_ATTR_STR)
	return FALSE;

    if (Red && Green && Blue)
        return sscanf(p -> U.Str, "%lf,%lf,%lf", Red, Green, Blue) == 3;
    else
        return sscanf(p -> U.Str, "%lf,%lf,%lf", &d, &d, &d) == 3;
} 

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to set a Floating point 64 RGB color attribute.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:		Where to place the TGB color attribute.              M
*   Red, Green, Blue:   Component of RGB color.                              M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrGetRGBColor, AttrSetWidth, AttrGetWidth, M
*   AttrSetIntAttrib, AttrSetObjectRGBColor, AttrGetRGBDoubleColor	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetRGBDoubleColor, attributes, color, rgb                            M
*****************************************************************************/
void AttrSetRGBDoubleColor(IPAttributeStruct **Attrs,
			   double Red,
			   double Green,
			   double Blue)
{
    IRIT_STATIC_DATA AttribNumType
        RGBAttribNum = ATTRIB_NAME_BAD_NUMBER;
    char SRGB[70];
    IPAttributeStruct *Attr;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL
    
    if (RGBAttribNum == ATTRIB_NAME_BAD_NUMBER)
        RGBAttribNum = AttrGetAttribNumber("rgbDouble");

    sprintf(SRGB, "%.16lf,%.16lf,%.16lf", Red, Green, Blue);

    Attr = AttrFindNumAttribute(*Attrs, RGBAttribNum);
    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> U.Str = IritStrdup(SRGB);
	Attr -> Type = IP_ATTR_STR;
    }
    else {
	Attr = _AttrMallocNumAttribute(RGBAttribNum, IP_ATTR_STR);
	Attr -> U.Str = IritStrdup(SRGB);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
} 

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to set a width attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Where to place the width attribute.                           M
*   Width:     New width.		                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrSetRGBColor, AttrGetRGBColor,            M
*   AttrGetWidth, AttrSetRealAttrib, AttrSetObjectWidth			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetWidth, attributes, width	                                     M
*****************************************************************************/
void AttrSetWidth(IPAttributeStruct **Attrs, IrtRType Width)
{
    IRIT_STATIC_DATA AttribNumType 
        WidthAttrNum = ATTRIB_NAME_BAD_NUMBER;
    IPAttributeStruct *Attr;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (WidthAttrNum == ATTRIB_NAME_BAD_NUMBER)
        WidthAttrNum = AttrGetAttribNumber("width");

    Attr = AttrFindNumAttribute(*Attrs, WidthAttrNum);
    
    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_REAL;
    }
    else {
	Attr = _AttrMallocNumAttribute(WidthAttrNum, IP_ATTR_REAL);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.R = Width;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to return a width attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     For which we would like to know the width of.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtRType:  Width or IP_ATTR_NO_WIDTH if no width set.	             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetColor, AttrGetColor, AttrSetRGBColor, AttrGetRGBColor,            M
*   AttrSetWidth, AttrGetRealAttrib, AttrGetObjectWidth			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetWidth, attributes, width 	                                     M
*****************************************************************************/
IrtRType AttrGetWidth(const IPAttributeStruct *Attrs)
{
    IRIT_STATIC_DATA AttribNumType
        WidthAttrNum = ATTRIB_NAME_BAD_NUMBER;
    IPAttributeStruct
        *Attr = NULL;

    DEBUG_GET_ATTR_DUMP_FUNC_CALL
    
    if (WidthAttrNum == ATTRIB_NAME_BAD_NUMBER) {
	if ((WidthAttrNum = AttrGetAttribNumber("width"))
						== ATTRIB_NAME_BAD_NUMBER)
	    return IP_ATTR_NO_WIDTH;
    }
    
    Attr = AttrFindNumAttribute(Attrs, WidthAttrNum);

    if (Attr != NULL) {
	if (Attr -> Type == IP_ATTR_REAL)
	    return Attr -> U.R;
	else if (Attr -> Type == IP_ATTR_STR) {
	    IrtRType r;

	    if (sscanf("%lf", Attr -> U.Str, &r) == 1)
	        return r;
	    else
	        return IP_ATTR_NO_WIDTH;
	}
	else
	    return IP_ATTR_NO_WIDTH;
    }
    else
        return IP_ATTR_NO_WIDTH;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set an integer attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Attribute list where to place new attribute.                  M
*   Name:      Name of the newly introduced attribute.                       M
*   Data:      Ingeter attribute to save.                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib   M
*   AttrSetUVAttrib, AttrSetIntAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetIntAttrib, attributes                                             M
*****************************************************************************/
void AttrSetIntAttrib(IPAttributeStruct **Attrs, const char *Name, int Data)
{
    AttrSetIntAttrib2(Attrs, AttrGetAttribNumber(Name), Data);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set an integer attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:      Attribute list where to place new attribute.                 M
*   AttribNum:  Unique ID derived from name of requested attribute.          M
*   Data:       Ingeter attribute to save.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib   M
*   AttrSetUVAttrib, AttrSetIntAttrib					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetIntAttrib2, attributes                                            M
*****************************************************************************/
void AttrSetIntAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       int Data)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(*Attrs, AttribNum);

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_INT;
    }
    else {
	Attr = _AttrMallocNumAttribute(AttribNum, IP_ATTR_INT);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.I = Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get an integer attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   Name:     Name of requested attribute.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:      Found attribute, or IP_ATTR_BAD_INT if not found.              M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetIntAttrib2 M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetIntAttrib, attributes                                             M
*****************************************************************************/
int AttrGetIntAttrib(const IPAttributeStruct *Attrs, const char *Name)
{
    return AttrGetIntAttrib2(Attrs, AttrGetAttribNumber(Name));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get an integer attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   AttrNum:  Unique ID derived from name of requested attribute.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:      Found attribute, or IP_ATTR_BAD_INT if not found.              M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetIntAttrib  M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetIntAttrib2, attributes                                            M
*****************************************************************************/
int AttrGetIntAttrib2(const IPAttributeStruct *Attrs, AttribNumType AttrNum)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(Attrs, AttrNum);

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Attr != NULL) {
	if (Attr -> Type == IP_ATTR_INT)
	    return Attr -> U.I;
	else if (Attr -> Type == IP_ATTR_REAL)
	    return (int) Attr -> U.R;
	else if (Attr -> Type == IP_ATTR_STR)
	    return atoi(Attr -> U.Str);
    }

    return IP_ATTR_BAD_INT;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a pointer attribute.		 			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Attribute list where to place new attribute.                  M
*   Name:      Name of the newly introduced attribute.                       M
*   Data:      Pointer attribute to save.                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetRefPtrAttrib, AttrGetRefPtrAttrib,		     M
*   AttrSetPtrAttrib2							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetPtrAttrib, attributes                                             M
*****************************************************************************/
void AttrSetPtrAttrib(IPAttributeStruct **Attrs,
		      const char *Name,
		      VoidPtr Data)
{
    AttrSetPtrAttrib2(Attrs, AttrGetAttribNumber(Name), Data);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a pointer attribute.		 			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:      Attribute list where to place new attribute.                 M
*   AttribNum:  Unique ID derived from name of requested attribute.          M
*   Data:       Pointer attribute to save.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetRefPtrAttrib, AttrGetRefPtrAttrib,		     M
*   AttrSetPtrAttrib							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetPtrAttrib2, attributes                                            M
*****************************************************************************/
void AttrSetPtrAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       VoidPtr Data)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(*Attrs, AttribNum);

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_PTR;
    }
    else {
	Attr = _AttrMallocNumAttribute(AttribNum, IP_ATTR_PTR);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.Ptr = Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a pointer attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   Name:     Name of requested attribute.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   VoidPtr: Found attribute, or NULL if not found. 		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetRefPtrAttrib, AttrGetRefPtrAttrib,		     M
*   AttrGetPtrAttrib2							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetPtrAttrib, attributes                                             M
*****************************************************************************/
VoidPtr AttrGetPtrAttrib(const IPAttributeStruct *Attrs, const char *Name)
{
    return AttrGetPtrAttrib2(Attrs, AttrGetAttribNumber(Name));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a pointer attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   AttrNum:  Unique ID derived from name of requested attribute.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   VoidPtr: Found attribute, or NULL if not found. 		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetRefPtrAttrib, AttrGetRefPtrAttrib,		     M
*   AttrGetPtrAttrib2							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetPtrAttrib2, attributes                                            M
*****************************************************************************/
VoidPtr AttrGetPtrAttrib2(const IPAttributeStruct *Attrs,
			  AttribNumType AttrNum)
{
    IPAttributeStruct
	*Attr = AttrFindNumAttribute(Attrs, AttrNum);

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Attr != NULL && Attr -> Type == IP_ATTR_PTR)
        return Attr -> U.Ptr;
    else
        return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a pointer reference attribute.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Attribute list where to place new attribute.                  M
*   Name:      Name of the newly introduced attribute.                       M
*   Data:      Pointer attribute to save.                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,		     M
*   AttrSetRefPtrAttrib2						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetRefPtrAttrib, attributes                                          M
*****************************************************************************/
void AttrSetRefPtrAttrib(IPAttributeStruct **Attrs,
			 const char *Name,
			 VoidPtr Data)
{
    AttrSetRefPtrAttrib2(Attrs, AttrGetAttribNumber(Name), Data);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a pointer reference attribute.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:      Attribute list where to place new attribute.                 M
*   AttribNum:  Unique ID derived from name of requested attribute.          M
*   Data:       Pointer attribute to save.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrGetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,		     M
*   AttrSetRefPtrAttrib							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetRefPtrAttrib2, attributes                                         M
*****************************************************************************/
void AttrSetRefPtrAttrib2(IPAttributeStruct **Attrs,
			  AttribNumType AttribNum,
			  VoidPtr Data)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(*Attrs, AttribNum);

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_REFPTR;
    }
    else {
	Attr = _AttrMallocNumAttribute(AttribNum, IP_ATTR_REFPTR);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.RefPtr = Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a pointer reference attribute.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   Name:     Name of requested attribute.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   VoidPtr: Found attribute, or NULL if not found. 		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib			     M
*   AttrGetRefPtrAttrib2						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRefPtrAttrib, attributes                                          M
*****************************************************************************/
VoidPtr AttrGetRefPtrAttrib(const IPAttributeStruct *Attrs, const char *Name)
{
    return AttrGetRefPtrAttrib2(Attrs, AttrGetAttribNumber(Name));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a pointer reference attribute.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   AttrNum:  Unique ID derived from name of requested attribute.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   VoidPtr: Found attribute, or NULL if not found. 		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrSetRealAttrib, M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib			     M
*   AttrGetRefPtrAttrib							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRefPtrAttrib2, attributes                                         M
*****************************************************************************/
VoidPtr AttrGetRefPtrAttrib2(const IPAttributeStruct *Attrs,
			     AttribNumType AttrNum)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(Attrs, AttrNum);

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Attr != NULL && Attr -> Type == IP_ATTR_REFPTR)
        return Attr -> U.RefPtr;
    else
        return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a IrtRType attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Attribute list where to place new attribute.                  M
*   Name:      Name of the newly introduce dattribute.                       M
*   Data:      IrtRType attribute to save.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetRealAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetRealAttrib, attributes                                            M
*****************************************************************************/
void AttrSetRealAttrib(IPAttributeStruct **Attrs,
		       const char *Name,
		       IrtRType Data)
{
    AttrSetRealAttrib2(Attrs, AttrGetAttribNumber(Name), Data);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a IrtRType attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:      Attribute list where to place new attribute.                 M
*   AttribNum:  Unique ID derived from name of requested attribute.          M
*   Data:       IrtRType attribute to save.                                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrGetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrSetRealAttrib					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetRealAttrib2, attributes                                           M
*****************************************************************************/
void AttrSetRealAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       IrtRType Data)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(*Attrs, AttribNum);

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_REAL;
    }
    else {
	Attr = _AttrMallocNumAttribute(AttribNum, IP_ATTR_REAL);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.R = Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a IrtRType attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   Name:     Name of requested attribute.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtRType: Found attribute, or IP_ATTR_BAD_REAL if not found.             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrGetRealAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRealAttrib, attributes                                            M
*****************************************************************************/
IrtRType AttrGetRealAttrib(const IPAttributeStruct *Attrs, const char *Name)
{
    return AttrGetRealAttrib2(Attrs, AttrGetAttribNumber(Name));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a IrtRType attribute.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   AttrNum:  Unique ID derived from name of requested attribute.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtRType: Found attribute, or IP_ATTR_BAD_REAL if not found.             M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrGetRealAttrib					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetRealAttrib2, attributes                                           M
*****************************************************************************/
IrtRType AttrGetRealAttrib2(const IPAttributeStruct *Attrs,
			    AttribNumType AttrNum)
{
    IPAttributeStruct
	*Attr = AttrFindNumAttribute(Attrs, AttrNum);

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Attr != NULL) {
	if (Attr -> Type == IP_ATTR_REAL)
	    return Attr -> U.R;
	else if (Attr -> Type == IP_ATTR_INT)
	    return (IrtRType) Attr -> U.I;
	else if (Attr -> Type == IP_ATTR_STR) {
	    IrtRType r;

	    if (sscanf(Attr -> U.Str, "%lf", &r) == 1)
	        return r;
	}
    }

    return IP_ATTR_BAD_REAL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a UV attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Attribute list where to place new attribute.                  M
*   Name:      Name of the newly introduced attribute.                       M
*   U, V:      UV attribute to save.            	                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrGetRealAttrib, AttrSetStrAttrib,                  M
*   AttrGetStrAttrib, AttrSetUVAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetUVAttrib, attributes                                              M
*****************************************************************************/
void AttrSetUVAttrib(IPAttributeStruct **Attrs,
		     const char *Name,
		     IrtRType U,
		     IrtRType V)
{
  AttrSetUVAttrib2(Attrs, AttrGetAttribNumber(Name), U, V);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a UV attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:      Attribute list where to place new attribute.                 M
*   AttribNum:  Unique ID derived from name of requested attribute.          M
*   U, V:       UV attribute to save.            	                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrGetRealAttrib, AttrSetStrAttrib,                  M
*   AttrGetStrAttrib, AttrSetUVAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetUVAttrib2, attributes                                             M
*****************************************************************************/
void AttrSetUVAttrib2(IPAttributeStruct **Attrs,
		      AttribNumType AttribNum,
		      IrtRType U,
		      IrtRType V)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(*Attrs, AttribNum);

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> Type = IP_ATTR_UV;
    }
    else {
	Attr = _AttrMallocNumAttribute(AttribNum, IP_ATTR_UV);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
    Attr -> U.UV[0] = (float) U;
    Attr -> U.UV[1] = (float) V;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a UV attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   Name:     Name of requested attribute.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   float *:  Found attribute, or NULL if not found.    		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrGetUVAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetUVAttrib, attributes                                              M
*****************************************************************************/
float *AttrGetUVAttrib(const IPAttributeStruct *Attrs, const char *Name)
{
    return AttrGetUVAttrib2(Attrs, AttrGetAttribNumber(Name));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a UV attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   AttrNum:  Unique ID derived from name of requested attribute.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   float *:  Found attribute, or NULL if not found.    		     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrSetStrAttrib, AttrGetStrAttrib, AttrGetUVAttrib,  M
*   AttrSetUVAttrib, AttrGetUVAttrib					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetUVAttrib2, attributes                                             M
*****************************************************************************/
float *AttrGetUVAttrib2(const IPAttributeStruct *Attrs, AttribNumType AttrNum)
{
    IPAttributeStruct
	*Attr = AttrFindNumAttribute(Attrs, AttrNum);

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Attr != NULL) {
	if (Attr -> Type == IP_ATTR_UV)
	    return Attr -> U.UV;
	else if (Attr -> Type == IP_ATTR_STR) {
	    float Uv[2];

	    if (sscanf(Attr -> U.Str, "%f %f", &Uv[0], &Uv[1]) == 2) {
		_AttrFreeAttributeData(Attr);		/* Free the string. */

		Attr -> Type = IP_ATTR_UV;
	        Attr -> U.UV[0] = Uv[0];
	        Attr -> U.UV[1] = Uv[1];

		return Attr -> U.UV;
	    }
	}
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a string attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:     Attribute list where to place new attribute.                  M
*   Name:      Name of the newly introduced attribute.                       M
*   Data:      String attribute to save.                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrGetRealAttrib, AttrGetStrAttrib, AttrGetUVAttrib, M
*   AttrSetUVAttrib, AttrSetStrAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetStrAttrib, attributes                                             M
*****************************************************************************/
void AttrSetStrAttrib(IPAttributeStruct **Attrs,
		      const char *Name,
		      const char *Data)
{
    AttrSetStrAttrib2(Attrs, AttrGetAttribNumber(Name), Data);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to set a string attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:      Attribute list where to place new attribute.                 M
*   AttribNum:  Unique ID derived from name of requested attribute.          M
*   Data:       String attribute to save.                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrGetRealAttrib, AttrGetStrAttrib, AttrGetUVAttrib, M
*   AttrSetUVAttrib, AttrSetStrAttrib					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrSetStrAttrib2, attributes                                            M
*****************************************************************************/
void AttrSetStrAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       const char *Data)
{
    IPAttributeStruct
        *Attr = AttrFindNumAttribute(*Attrs, AttribNum);

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr) {
	_AttrFreeAttributeData(Attr);
	Attr -> U.Str = IritStrdup(Data);
	Attr -> Type = IP_ATTR_STR;
    }
    else {
	Attr = _AttrMallocNumAttribute(AttribNum, IP_ATTR_STR);
	Attr -> U.Str = IritStrdup(Data);
	Attr -> Pnext = *Attrs;
	*Attrs = Attr;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a string attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   Name:     Name of requested attribute.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *:   Found attribute, or NULL if not found.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrGetRealAttrib, AttrSetStrAttrib, AttrGetUVAttrib, M
*   AttrSetUVAttrib, AttrGetStrAttrib2					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetStrAttrib, attributes                                             M
*****************************************************************************/
const char *AttrGetStrAttrib(const IPAttributeStruct *Attrs, const char *Name)
{
    return AttrGetStrAttrib2(Attrs, AttrGetAttribNumber(Name));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to get a string attribute.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search for requested attribute.              M
*   AttrNum:  Unique ID derived from name of requested attribute.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *:   Found attribute, or NULL if not found.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   AttrSetIntAttrib, AttrGetIntAttrib, AttrSetPtrAttrib, AttrGetPtrAttrib,  M
*   AttrSetRealAttrib, AttrGetRealAttrib, AttrSetStrAttrib, AttrGetUVAttrib, M
*   AttrSetUVAttrib, AttrGetStrAttrib					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetStrAttrib2, attributes                                            M
*****************************************************************************/
const char *AttrGetStrAttrib2(const IPAttributeStruct *Attrs,
			      AttribNumType AttrNum)
{
    IPAttributeStruct
	*Attr = AttrFindNumAttribute(Attrs, AttrNum);

    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (Attr != NULL && Attr -> Type == IP_ATTR_STR)
        return Attr -> U.Str;
    else
        return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to aid in scanning a list of attributes.                           M
*   If TraceAttrs != NULL, a ptr to its attribute list is saved and the      M
* next attribute is returned every call until the end of the list is         M
* reached, in which NULL is returned.                                        M
*   FirstAttrs should be NULL in all but the first call in the sequence.     M
*   Attributes with names starting with an underscore '_' are assumed to be  M
* temporary or internal and are skipped.				     M 
*                                                                            *
* PARAMETERS:                                                                M
*   TraceAttrs:    If not NULL, contains the previously returned attribute.  M
*   FirstAttrs:    First attribute in list, usually NULL in all but the      M
*                  first invocation.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   const IPAttributeStruct *:  Next attribute in list.                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrTraceAttributes, attributes                                          M
*****************************************************************************/
const IPAttributeStruct *AttrTraceAttributes(
					 const IPAttributeStruct *TraceAttrs,
				         const IPAttributeStruct *FirstAttrs)
{
    DEBUG_GET_ATTR_DUMP_FUNC_CALL

    if (!_AttrNamesHashTbl)
        AttrInitHashTbl();

    if (FirstAttrs != NULL)
	TraceAttrs = FirstAttrs;
    else if (TraceAttrs != NULL)
	TraceAttrs = TraceAttrs -> Pnext;
    else
	return NULL;

    while (TraceAttrs) {
	const char 
	    *AttribName = AttrGetAttribName(TraceAttrs);

	if (AttribName[0] != '_') {
	    return TraceAttrs;
	}
	else
	    TraceAttrs = TraceAttrs -> Pnext;
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to convert an attribute to a string.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attr:      To convert to a string.                                       M
*   DataFileFormat:  If TRUE, the attribute is formated in the IRIT data     M
*		file format.  Otherwise, just the attribute value is	     M
*		returned as a string.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *:    A pointer to a static string representing/describing    M
*             the given attribute.					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Attr2String, attributes                                                  M
*****************************************************************************/
const char *Attr2String(const IPAttributeStruct *Attr, int DataFileFormat)
{
    IRIT_STATIC_DATA char Str[IRIT_LINE_LEN_XLONG];

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    Str[0] = 0;

    if (DataFileFormat) {
        const char  
	    *AttribName = AttrGetAttribName(Attr);

	switch (Attr -> Type) {
	    case IP_ATTR_INT:
	        sprintf(Str, "[%s %d]", AttribName, Attr -> U.I);
		break;
	    case IP_ATTR_REAL:
		sprintf(Str, "[%s %g]", AttribName, Attr -> U.R);
		break;
	    case IP_ATTR_UV:
		sprintf(Str, "[%s \"%g %g\"]", AttribName, 
			Attr -> U.UV[0], Attr -> U.UV[1]);

		break;
	    case IP_ATTR_STR:
		if (strchr(Attr -> U.Str, '"') != NULL ||
		    strchr(Attr -> U.Str, '\\') != NULL) {
		    /* Need to escape the quotes. */
		    int i, j,
		        Len = (int) strlen(Attr -> U.Str);

		    sprintf(Str, "[%s \"", AttribName);
		    j = (int) strlen(Str);

		    /* Need to escape the internal '"' or '\'. */
		    for (i = 0; i < Len; i++) {
		        if (Attr -> U.Str[i] == '"' ||
			    Attr -> U.Str[i] == '\\') {
			    Str[j++] = '\\';
			}
			Str[j++] = Attr -> U.Str[i];
		    }
		    Str[j] = 0;
		    strcat(Str, "\"]");
		}
		else if (!IRT_STR_ZERO_LEN(Attr -> U.Str))
		    sprintf(Str, "[%s \"%s\"]", AttribName, Attr -> U.Str);
		else
		    sprintf(Str, "[%s]", AttribName);
		break;
	    case IP_ATTR_OBJ:
		sprintf(Str, "[%s _OBJ_ATTR_NOT_CNVRTED_]", AttribName);
		break;
	    case IP_ATTR_PTR:
		sprintf(Str, "[%s _PTR_ATTR_NOT_CNVRTED_]", AttribName);
		break;
	    case IP_ATTR_REFPTR:
		sprintf(Str, "[%s _REFPTR_ATTR_NOT_CNVRTED_]", AttribName);
		break;
	    default:
		IRIT_FATAL_ERROR("Undefined attribute type 3");
		break;
	}
    }
    else {
	switch (Attr -> Type) {
	    case IP_ATTR_INT:
	        sprintf(Str, "%d", Attr -> U.I);
		break;
	    case IP_ATTR_REAL:
		sprintf(Str, "%g", Attr -> U.R);
		break;
	    case IP_ATTR_UV:
		sprintf(Str, "\"%g %g\"", Attr -> U.UV[0], Attr -> U.UV[1]);
		break;
	    case IP_ATTR_STR:
		if (strchr(Attr -> U.Str, '"') != NULL) {
		    /* Need to escape the quotes. */
		    int i, j,
		        Len = (int) strlen(Attr -> U.Str);

		    strcpy(Str, "\"");
		    j = (int) strlen(Str);

		    /* Need to quote the string or escape the internal " */
		    for (i = 0; i < Len; i++) {
			if (Attr -> U.Str[i] == '"')
			    Str[j++] = '\\';
			Str[j++] = Attr -> U.Str[i];
		    }
		    Str[j] = 0;
		    strcat(Str, "\"");
		}
		else if (!IRT_STR_ZERO_LEN(Attr -> U.Str))
		    sprintf(Str, "\"%s\"", Attr -> U.Str);
		else
		    sprintf(Str, "%s", Attr -> U.Str);
		if (IRT_STR_ZERO_LEN(Str))
		    strcpy(Str, "\"\"");
		break;
	    case IP_ATTR_OBJ:
		strcpy(Str, "_OBJ_ATTR_NOT_CNVRTED_");
		break;
	    case IP_ATTR_PTR:
		strcpy(Str, "_PTR_ATTR_NOT_CNVRTED_");
		break;
	    case IP_ATTR_REFPTR:
		strcpy(Str, "_REFPTR_ATTR_NOT_CNVRTED_");
		break;
	    default:
		IRIT_FATAL_ERROR("Undefined attribute type 4");
		break;
	}
    }

    return Str;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to reverse the given Attr list.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attr:   To reverse, in place.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPAttributeStruct *:  The reversed list, in place.			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrReverseAttributes, attributes                                        M
*****************************************************************************/
IPAttributeStruct *AttrReverseAttributes(IPAttributeStruct *Attr)
{
    IPAttributeStruct
	*NewAttrs = NULL;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    while (Attr) {
	IPAttributeStruct
	    *Pnext = Attr -> Pnext;

	Attr -> Pnext = NewAttrs;
	NewAttrs = Attr;

	Attr = Pnext;
    }

    return NewAttrs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to search for an attribute by Name.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search.                                      M
*   Name:     Attribute to search by this name.                              M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPAttributeStruct *:  Attribute if found, otherwise NULL.                M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrFindAttribute, attributes                                            M
*****************************************************************************/
IPAttributeStruct *AttrFindAttribute(const IPAttributeStruct *Attrs,
				     const char *Name)
{
    int AttribNum;

    if (Attrs == NULL)
	return NULL;

    if ((AttribNum = AttrGetAttribNumber(Name)) == ATTRIB_NAME_BAD_NUMBER)
        return NULL;
    else
        return AttrFindNumAttribute(Attrs, AttribNum);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Allocated a new attribute structure.                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   Name:   Name of newly created attribute.                                 M
*   Type:   Type of newly created attribute.                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPAttributeStruct *:  The newly created attribute.                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   _AttrMallocAttribute                                                     M
*****************************************************************************/
IPAttributeStruct *_AttrMallocAttribute(const char *Name, IPAttributeType Type)
{
    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    return _AttrMallocNumAttribute(AttrGetAttribNumber(Name), Type);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Allocated a new attribute structure.                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   AttribNum:  Number of newly created attribute.                           M
*   Type:       Type of newly created attribute.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPAttributeStruct *:  The newly created attribute.                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   _AttrMallocNumAttribute                                                  M
*****************************************************************************/
IPAttributeStruct *_AttrMallocNumAttribute(AttribNumType AttribNum, 
					   IPAttributeType Type)
{
    IPAttributeStruct *Attr;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (AttrFreedList != NULL) {
	Attr = AttrFreedList;
	AttrFreedList = AttrFreedList -> Pnext;
    }
    else {
	/* Allocate AllocateNumObj objects, returns first one as new   */
	/* and chain together the rest of them into the free list.     */
#ifdef DEBUG_ATTR_MALLOC
	Attr = (IPAttributeStruct *) IritMalloc(sizeof(IPAttributeStruct));
#else
	int i, AllocateNumObj;

	AllocateNumObj = getenv("IRIT_MALLOC") ? 1 : ATTR_ALLOCATE_NUM;

	if ((Attr = (IPAttributeStruct *) IritMalloc(sizeof(IPAttributeStruct)
					       * AllocateNumObj)) != NULL) {
	    for (i = 1; i < AllocateNumObj - 1; i++)
		Attr[i].Pnext = &Attr[i + 1];
	    Attr[AllocateNumObj - 1].Pnext = NULL;
	    if (AllocateNumObj > 1)
		AttrFreedList = &Attr[1];
	}
#endif /* DEBUG_ATTR_MALLOC */
    }

    Attr -> Type = Type;
    Attr -> _AttribNum = AttribNum;
    Attr -> Pnext = NULL;

    return Attr;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to remove and delete the attribute named Name from the given       M
* Attr list.                                                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:   To search for an attribute named Name and remove and delete it. M
*   Name:    Name of attribute to remove and delete.                         M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrFreeOneAttribute, attributes                                         M
*****************************************************************************/
void AttrFreeOneAttribute(IPAttributeStruct **Attrs, const char *Name)
{
    AttribNumType AttribNum;
    IPAttributeStruct *TmpAttr,
	*Attr = *Attrs;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if (Attr == NULL ||
	(AttribNum = AttrGetAttribNumber(Name)) == ATTRIB_NAME_BAD_NUMBER)
        return;

    if (Attr) {
	if (Attr -> _AttribNum == AttribNum) {
	    /* It is the first one - delete first in list. */
	    *Attrs = Attr -> Pnext;

	    _AttrFreeAttributeData(Attr);
#ifdef DEBUG_ATTR_MALLOC
	    IritFree(Attr);
#else
	    IRIT_LIST_PUSH(Attr, AttrFreedList);
#endif /* DEBUG_ATTR_MALLOC */
	}
	else {
	    /* Search the rest of the list and delete if found. */
	    while (Attr -> Pnext != NULL) {
		if (Attr -> Pnext -> _AttribNum == AttribNum) {
		    TmpAttr = Attr -> Pnext;
		    Attr -> Pnext = TmpAttr -> Pnext;

		    _AttrFreeAttributeData(TmpAttr);
#ifdef DEBUG_ATTR_MALLOC
		    IritFree(TmpAttr);
#else
		    IRIT_LIST_PUSH(TmpAttr, AttrFreedList);
#endif /* DEBUG_ATTR_MALLOC */
		}
		else
		    Attr = Attr -> Pnext;
	    }
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to remove anddelete all attributes of the given Attr list.         M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:   To remove and delete.                                           M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrFreeAttributes, attributes                                           M
*****************************************************************************/
void AttrFreeAttributes(IPAttributeStruct **Attrs)
{
    IPAttributeStruct *Attr, *Last;

    DEBUG_SET_ATTR_DUMP_FUNC_CALL

    if ((Last = *Attrs) == NULL)
	return;

    for (Attr = *Attrs; Attr != NULL; Attr = Attr -> Pnext) {
        if (Attr -> Pnext != NULL)
	    Last = Attr -> Pnext;

	_AttrFreeAttributeData(Attr);
    }

#ifdef DEBUG_ATTR_MALLOC
    while (*Attrs != NULL) {
        Attr = *Attrs;
	*Attrs = Attr -> Pnext;
	IritFree(Attr);
    }
#else
    /* Chain into attributes' free list. */
    Last -> Pnext = AttrFreedList;
    AttrFreedList = *Attrs;
#endif /* DEBUG_ATTR_MALLOC */

    *Attrs = NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the NULL terminated list of attribute's names to copy.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   AttrNames:      Vector of attribute names to use in attr list copy.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char **:   Old list of valid attributes to copy.                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrCopyValidAttrList, attributes                                        M
*****************************************************************************/
const char **AttrCopyValidAttrList(const char **AttrNames)
{
    const char
        **OldNames = _AttrValidAttrList;

    _AttrValidAttrList = AttrNames;

    return OldNames;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to copy an attribute list.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   Src:       Attribute list to duplicate.                                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPAttributeStruct *:   Duplicated attribute list.                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrCopyAttributes, attributes                                           M
*****************************************************************************/
IPAttributeStruct *AttrCopyAttributes(const IPAttributeStruct *Src)
{
    IPAttributeStruct
	*DestAttr = NULL,
	*Dest = NULL;

    for ( ; Src != NULL; Src = Src -> Pnext) {
	const char 
	    *AttribName = AttrGetAttribName(Src);

	if (_AttrValidAttrList != NULL) {
	    int i;

	    /* Have a list of valid attr name? Copy only if in the list. */
	    for (i = 0; _AttrValidAttrList[i] != NULL; i++)
	        if (stricmp(_AttrValidAttrList[i], AttribName) == 0)
		    break;
	    if (_AttrValidAttrList[i] == NULL)
	        continue;
	}
	else if (AttribName[0] == '_')  /* Do not copy internal attributes. */
	    continue;

	if (Dest == NULL) {	
	    Dest = DestAttr = AttrCopyOneAttribute(Src);
	}
	else {
	    DestAttr -> Pnext = AttrCopyOneAttribute(Src);

	    DestAttr = DestAttr -> Pnext;
	}
    }

    return Dest;
}
