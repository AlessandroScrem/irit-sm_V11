/*****************************************************************************
* Setting attributes for objects.					     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 0.2, Mar. 1990   *
*****************************************************************************/

#ifndef MISCATTR_H
#define MISCATTR_H

#include "irit_sm.h"

typedef enum {
    IP_ATTR_NONE,
    IP_ATTR_INT,
    IP_ATTR_REAL,
    IP_ATTR_UV,
    IP_ATTR_STR,
    IP_ATTR_OBJ,
    IP_ATTR_PTR,
    IP_ATTR_REFPTR
} IPAttributeType;

#define IP_ATTR_BAD_INT		-32767
#define IP_ATTR_BAD_REAL	1e30
#define IP_ATTR_NO_COLOR	999
#define IP_ATTR_NO_WIDTH	1e30

#define IP_ATTR_IS_BAD_INT(I)	((I) == IP_ATTR_BAD_INT)
#define IP_ATTR_IS_BAD_REAL(R)	((R) > IP_ATTR_BAD_REAL / 10.0)
#define IP_ATTR_IS_BAD_COLOR(C)	((C) == IP_ATTR_NO_COLOR)
#define IP_ATTR_IS_BAD_WIDTH(W)	((W) > IP_ATTR_NO_WIDTH / 10.0)

#define IP_ATTR_RESET_ATTRS(Attr)	{ (Attr) = NULL; }
#define IP_ATTR_FREE_ATTRS(Attr)	{ if ((Attr) != NULL) \
					      AttrFreeAttributes(&(Attr)); }
#define IP_ATTR_COPY_ATTRS(Attr) \
	(Attr) != NULL ? AttrCopyAttributes(Attr) : NULL;
#define IP_ATTR_COPY_ATTRS2(NewAttr, OldAttr) { \
	NewAttr = IP_ATTR_COPY_ATTRS(OldAttr); }
#define IP_ATTR_SAFECOPY_ATTRS(NewAttr, OldAttr) { \
	IP_ATTR_FREE_ATTRS(NewAttr); \
	NewAttr = IP_ATTR_COPY_ATTRS(OldAttr); }

#define IP_ATTR_ATTR_EXIST(Attrs, Name) (AttrFindAttribute(Attrs, Name) != NULL)

#define IP_ATTR_INIT_UNIQUE_ID_NUM(AttrID, StrID) \
    if (AttrID == ATTRIB_NAME_BAD_NUMBER) \
	AttrID = AttrGetAttribNumber(StrID);

#define IP_ATTR_IRIT_COLOR_TABLE_SIZE	16

typedef unsigned int AttribNumType;

#define ATTRIB_NAME_BAD_NUMBER    ((AttribNumType) - 1)

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

IRIT_GLOBAL_DATA_HEADER int AttrIritColorTable[][3];

IRIT_GLOBAL_DATA_HEADER struct IritHashTableStruct *_AttrNamesHashTbl;

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

/*****************************************************************************
* Attributes - an attribute has a name and can be one of the following:	     *
* an integer, real, string, or a pointer to an Object.			     *
*****************************************************************************/
typedef struct IPAttributeStruct {
    struct IPAttributeStruct *Pnext;
    IPAttributeType Type;
    union {
	char *Str;
   	int I;
	IrtRType R;
	float UV[2];
	struct IPObjectStruct *PObj;
	VoidPtr Ptr;
	VoidPtr RefPtr;
    } U;
    AttribNumType _AttribNum;
} IPAttributeStruct;

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

void AttrGetIndexColor(int Color, int *Red, int *Green, int *Blue);
void AttrSetColor(IPAttributeStruct **Attrs, int Color);
int AttrGetColor(const IPAttributeStruct *Attrs);
void AttrSetRGBColor(IPAttributeStruct **Attrs, int Red, int Green, int Blue);
int AttrGetRGBColor(const IPAttributeStruct *Attrs,
		    int *Red,
		    int *Green,
		    int *Blue);
int AttrGetRGBColor2(const IPAttributeStruct *Attrs, 
		     const char *Name,
		     int *Red, 
		     int *Green, 
		     int *Blue);
void AttrSetWidth(IPAttributeStruct **Attrs, IrtRType Width);
IrtRType AttrGetWidth(const IPAttributeStruct *Attrs);

void AttrSetIntAttrib(IPAttributeStruct **Attrs, const char *Name, int Data);
void AttrSetIntAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       int Data);
int AttrGetIntAttrib(const IPAttributeStruct *Attrs, const char *Name);
int AttrGetIntAttrib2(const IPAttributeStruct *Attrs, AttribNumType AttrNum);

void AttrSetRealAttrib(IPAttributeStruct **Attrs,
		       const char *Name,
		       IrtRType Data);
void AttrSetRealAttrib2(IPAttributeStruct **Attrs,
			AttribNumType AttribNum,
			IrtRType Data);
IrtRType AttrGetRealAttrib(const IPAttributeStruct *Attrs, const char *Name);
IrtRType AttrGetRealAttrib2(const IPAttributeStruct *Attrs,
			    AttribNumType AttrNum);

void AttrSetUVAttrib(IPAttributeStruct **Attrs,
		     const char *Name,
		     IrtRType U,
		     IrtRType V);
void AttrSetUVAttrib2(IPAttributeStruct **Attrs,
		      AttribNumType AttribNum,
		      IrtRType U,
		      IrtRType V);
float *AttrGetUVAttrib(const IPAttributeStruct *Attrs, const char *Name);
float *AttrGetUVAttrib2(const IPAttributeStruct *Attrs, AttribNumType AttrNum);

void AttrSetPtrAttrib(IPAttributeStruct **Attrs,
		      const char *Name,
		      VoidPtr Data);
void AttrSetPtrAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       VoidPtr Data);
VoidPtr AttrGetPtrAttrib(const IPAttributeStruct *Attrs, const char *Name);
VoidPtr AttrGetPtrAttrib2(const IPAttributeStruct *Attrs,
			  AttribNumType AttrNum);

void AttrSetRefPtrAttrib(IPAttributeStruct **Attrs,
			 const char *Name,
			 VoidPtr Data);
void AttrSetRefPtrAttrib2(IPAttributeStruct **Attrs,
			  AttribNumType AttribNum,
			  VoidPtr Data);
VoidPtr AttrGetRefPtrAttrib(const IPAttributeStruct *Attrs, const char *Name);
VoidPtr AttrGetRefPtrAttrib2(const IPAttributeStruct *Attrs,
			     AttribNumType AttrNum);

void AttrSetStrAttrib(IPAttributeStruct **Attrs,
		      const char *Name,
		      const char *Data);
void AttrSetStrAttrib2(IPAttributeStruct **Attrs,
		       AttribNumType AttribNum,
		       const char *Data);
const char *AttrGetStrAttrib(const IPAttributeStruct *Attrs, const char *Name);
const char *AttrGetStrAttrib2(const IPAttributeStruct *Attrs,
			      AttribNumType AttrNum);

const IPAttributeStruct *AttrTraceAttributes(
					  const IPAttributeStruct *TraceAttrs,
					  const IPAttributeStruct *FirstAttrs);
const char *Attr2String(const IPAttributeStruct *Attr, int DataFileFormat);

IPAttributeStruct *AttrReverseAttributes(IPAttributeStruct *Attr);

void AttrFreeOneAttribute(IPAttributeStruct **Attrs, const char *Name);
void AttrFreeAttributes(IPAttributeStruct **Attrs);

IPAttributeStruct *AttrFindAttribute(const IPAttributeStruct *Attrs,
				     const char *Name);

IPAttributeStruct *_AttrMallocAttribute(const char *Name,
					IPAttributeType Type);
IPAttributeStruct *_AttrMallocNumAttribute(AttribNumType AttribNum, 
					   IPAttributeType Type);
void _AttrFreeAttributeData(IPAttributeStruct *Attr);

const char **AttrCopyValidAttrList(const char **AttrNames);
IPAttributeStruct *AttrCopyAttributes(const IPAttributeStruct *Src);
IPAttributeStruct *AttrCopyOneAttribute(const IPAttributeStruct *Src);

/* Functions in miscatt1.c */
const char *AttrGetAttribName(const IPAttributeStruct *Attr);
AttribNumType AttrGetAttribNumber(const char *AttribName);
IPAttributeStruct *AttrFindNumAttribute(const IPAttributeStruct *Attrs, 
					AttribNumType AttrNum);
void AttrInitHashTbl(void);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /* MISCATTR_H */
