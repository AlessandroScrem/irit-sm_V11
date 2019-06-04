/*****************************************************************************
*   "Irit" - the 3d (not only polygonal) solid modeller - DB support.	     *
*									     *
* Written by:  Gershon Elber				Ver 0.2, Sep. 2007   *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
*   Module to handle the data base of objects - fetch, insert, delete etc... *
*****************************************************************************/

#include <stdio.h>
#include "program.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "objects.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Print some useful information on the global list of objects.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   None			                                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritDBPrintAllObjs                                                       M
*****************************************************************************/
void IritDBPrintAllObjs(void)
{
    int i;

    for (i = 0; i < IritDBGetStackSize(); i++)
        PrintIritObjectList(IritDBGetStackObject(i));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A stab function IRIT.  Ignored.                                          M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:     Object to update its function creation parameters as a string. M
*   FuncID:      The ID of the function.				     M
*   ParamList:   Parameters' list of this specific instantiation of PObj.    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritDBUpdateParams                                                       M
*****************************************************************************/
void IritDBUpdateParams(IPObjectStruct *PObj,
			const char *FuncID,
			IPObjectStruct *ParamList)
{
    IRIT_WNDW_FPRINTF1("Error: Stab function - does nothing in IRIT mode and is ignored.\n");
}
