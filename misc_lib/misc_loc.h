/******************************************************************************
* Misc_loc.h - header file for the Misc library.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, Dec. 96.					      *
******************************************************************************/

#ifndef MISC_LOC_H
#define MISC_LOC_H

/******************************************************************************
* This macro is called when the library has detected an unrecoverable error.  *
* Default action is to call MiscFatalError, but you may want to reroute this  *
* to invoke your handler and recover yourself (by long jump for example).     *
******************************************************************************/
#define MISC_FATAL_ERROR(Msg, Desc)	MiscFatalError(Msg, Desc)

#include "inc_irit/misc_lib.h"
#include "inc_irit/miscattr.h"

IrtBType *_IrtImgVerifyAlignment(IrtBType *Data,
				 int *MaxX,
				 int *MaxY,
				 int Alpha);

#endif /* MISC_LOC_H */
