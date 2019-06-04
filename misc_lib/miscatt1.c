/*****************************************************************************
* MiscAtt1.c - handling hashing of string attributes as integers.	     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Ronen Lev and Gershon Elber		Ver 0.2, Mar. 1990   *
*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "inc_irit/irit_sm.h"
#include "misc_loc.h"
#include "inc_irit/miscattr.h"

#define ATTR_MAX_NAME_LEN 31

typedef struct _AttribNumInfoStruct {
    AttribNumType AttribNum;
    const char *AttrName;
} _AttribNumInfoStruct;

#define ATTRIB_NAME_HASH_SIZE	256
#define ATTRIB_NAME_HASH_SIZE1	255
#define REQUIRED_KEY_SHIFTS	24
#define REQUIRED_KEY_MASK	0x000000FF

IRIT_GLOBAL_DATA IritHashTableStruct
    *_AttrNamesHashTbl = NULL;

IRIT_STATIC_DATA int
    _AttrLastElemNumInRow[ATTRIB_NAME_HASH_SIZE];

static int AttrHashAttribName(const char *AttribName);
static int AttrHashCmpNameAux(VoidPtr Data1, VoidPtr Data2);
static int AttrHashCmpNumAux(VoidPtr Data1, VoidPtr Data2);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   This function returns the name of the attribute.                         M
*                                                                            *
* PARAMETERS:                                                                M
*   Attr:   The attribute structure.                                         M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *: The name of the attribute, "__undefined__" if don't exist. M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetAttribName                                                        M
*****************************************************************************/
const char *AttrGetAttribName(const IPAttributeStruct *Attr)
{
    int Key = ((Attr -> _AttribNum) >> REQUIRED_KEY_SHIFTS) & REQUIRED_KEY_MASK;
    _AttribNumInfoStruct CmpAttribNum, *ResAttribNum;

    if (!_AttrNamesHashTbl)
        return "__undefined__";

    CmpAttribNum.AttribNum = Attr -> _AttribNum;
    ResAttribNum = IritHashTableFind(_AttrNamesHashTbl, &CmpAttribNum,
				     AttrHashCmpNumAux, Key);
    if (ResAttribNum != NULL)
        return ResAttribNum -> AttrName;
    else
        return "__undefined__";
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   This function returns the matching attribute number to the specified     M
* attribute name. If the attribute name doesn't exist it is created.         M
*                                                                            *
* PARAMETERS:                                                                M
*   AttribName: The Attribute name to seek and create if not found.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   AttribNumType: The attribute number.                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrGetAttribNumber                                                      M
*****************************************************************************/
AttribNumType AttrGetAttribNumber(const char *AttribName)
{
    int Key = AttrHashAttribName(AttribName);
    _AttribNumInfoStruct CmpAttribName, *ResAttribName;

    if (_AttrNamesHashTbl == NULL)
        AttrInitHashTbl();

    /* Try to find this AttribName in the hash table. */
    CmpAttribName.AttrName = AttribName;
    ResAttribName = IritHashTableFind(_AttrNamesHashTbl, &CmpAttribName,
				      AttrHashCmpNameAux, Key);
    if (ResAttribName != NULL)
	return ResAttribName -> AttribNum;

    /* Insert a new copy of this AttribName into the hash table. */
    ResAttribName = 
        (_AttribNumInfoStruct *) IritMalloc(sizeof(_AttribNumInfoStruct));
    ResAttribName -> AttrName = IritStrdup(AttribName);
    if (!IritHashTableInsert(_AttrNamesHashTbl, ResAttribName,
			     AttrHashCmpNameAux, Key, FALSE)) {
        ResAttribName -> AttribNum = (Key << REQUIRED_KEY_SHIFTS) +
	                             (_AttrLastElemNumInRow[Key]++);
	return ResAttribName -> AttribNum;
    }
    else 
        IRIT_FATAL_ERROR("Error in the Attrib name hash table.");

    return ATTRIB_NAME_BAD_NUMBER;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Routine to search for an attribute by ID numeric index.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Attrs:    Attribute list to search.                                      M
*   AttrNum:  Attribute numeric index to search by this number.              M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPAttributeStruct *:  Attribute if found, otherwise NULL.                M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrFindNumAttribute                                                     M
*****************************************************************************/
IPAttributeStruct *AttrFindNumAttribute(const IPAttributeStruct *Attrs, 
					AttribNumType AttrNum)
{
    for ( ; Attrs != NULL; Attrs = Attrs -> Pnext)
        if (Attrs -> _AttribNum == AttrNum)
	    return (IPAttributeStruct *) Attrs;

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Calculate a hash value for an attribute name.                            *
* Currently returns the checksum of the attribute name.                      *
*                                                                            *
* PARAMETERS:                                                                *
*   const  char *: The attribute name to be hashed.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  The hash value.                                                    *
*****************************************************************************/
static int AttrHashAttribName(const char *AttribName)
{
    int Result = 0,
        CurCharInd = 0;
    const char
        *p = AttribName; 

    if (AttribName == NULL)
        return 0;

    while (*p && CurCharInd++ < ATTR_MAX_NAME_LEN)
	Result += (*p++ & ~0x20);    /* & with ~0x20 to make up/locase same. */

    return Result & ATTRIB_NAME_HASH_SIZE1;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Initialize the hash table for attributes names.                          M
*                                                                            *
* PARAMETERS:                                                                M
*   None                                                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   AttrInitHashTbl                                                          M
*****************************************************************************/
void AttrInitHashTbl(void)
{
    if (_AttrNamesHashTbl != NULL)
        return;

    _AttrNamesHashTbl = 
        IritHashTableCreate(0.0, ATTRIB_NAME_HASH_SIZE - 1, IRIT_EPS, 
			    ATTRIB_NAME_HASH_SIZE);

    IRIT_ZAP_MEM(_AttrLastElemNumInRow, sizeof(int) * ATTRIB_NAME_HASH_SIZE);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compare AttributeNumInfo by their string.                                *
*                                                                            *
* PARAMETERS:                                                                *
*   Data1: A pointer to a _AttribNumInfoStruct to be compared.               *
*   Data2: A pointer to a _AttribNumInfoStruct to be compared.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: The comparison result.                                              *
*****************************************************************************/
static int AttrHashCmpNameAux(VoidPtr Data1, VoidPtr Data2)
{
    return stricmp(((_AttribNumInfoStruct *) Data1) -> AttrName,
		   ((_AttribNumInfoStruct *) Data2) -> AttrName);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Compare AttributeNumInfo by their AttributeNumber.                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Data1: A pointer to a _AttribNumInfoStruct to be compared.               *
*   Data2: A pointer to a _AttribNumInfoStruct to be compared.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: The comparison result.                                              *
*****************************************************************************/
static int AttrHashCmpNumAux(VoidPtr Data1, VoidPtr Data2)
{
    return ((_AttribNumInfoStruct *) Data1) -> AttribNum - 
           ((_AttribNumInfoStruct *) Data2) -> AttribNum;
}
