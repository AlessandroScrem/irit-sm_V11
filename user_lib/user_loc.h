/******************************************************************************
* User_loc.h - header file for the user interaction library.		      *
* This header is also the interface header to the world.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, May. 95.					      *
******************************************************************************/

#ifndef USER_LOC_H
#define USER_LOC_H

/******************************************************************************
* This macro is called when the library has detected an unrecoverable error.  *
* Default action is to call UserFatalError, but you may want to reroute this  *
* to invoke your handler and recover yourself (by long jump for example).     *
******************************************************************************/
#define USER_FATAL_ERROR(Msg)	UserFatalError(Msg)

#include "inc_irit/user_lib.h"
#include "inc_irit/user_gc.h"

IRIT_GLOBAL_DATA_HEADER jmp_buf UserGCStartOfProcess;
IRIT_GLOBAL_DATA_HEADER const int USER_GC_TIME_ZONE;

void UserGCSaveProcessedObjects(UserGCProblemDefinitionStruct *Problem);
void UserGCDeleteAllVisMaps(UserGCProblemDefinitionStruct *Problem);
void UserGCLoadProcessedObjects(UserGCProblemDefinitionStruct *Problem,
                                IPObjectStruct **GeoObj,
                                IPObjectStruct **Obstacles);
void UserGCDeleteProcessedObjects(UserGCProblemDefinitionStruct *Problem);
int UserGCGetOPsNum(UserGCProblemDefinitionStruct *Problem);
void UserGCPrepareScene(UserGCProblemDefinitionStruct *Problem);
void UserGCCreateVisMaps(UserGCProblemDefinitionStruct *Problem);
IrtImgPixelStruct *UserGCLoadVisMap2(UserGCProblemDefinitionStruct *Problem,
                                     int Index);
int UserGCInterpretOPGroupsSuggestion(UserGCProblemDefinitionStruct *Problem);

#endif /* USER_LOC_H */
