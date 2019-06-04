/******************************************************************************
* Filter to convert IRIT data files to a WebGL HTML file.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by:  Avi Kaplan				Ver 1.0, Feb 2012     *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/misc_lib.h"

#ifdef NO_CONCAT_STR 
IRIT_STATIC_DATA const char
    *VersionStr = "Irit2WGL		Version 11,	Gershon Elber,\n\
	 (C) Copyright 1989-2015 Gershon Elber, Non commercial use only.";
#else
IRIT_STATIC_DATA const char
    *VersionStr = "Irit2WGL		" IRIT_VERSION
	",	Gershon Elber,	" __DATE__ ",   " __TIME__ "\n"
	IRIT_COPYRIGHT ", Non commercial use only.";
#endif /* NO_CONCAT_STR */

IRIT_STATIC_DATA const char
    *CtrlStr = "irit2wgl l%- 4%- F%-PolyOpti|FineNess!d!F C%- w%-CanvasWidth!d h%-CanvasHeight!d b%-R|G|B!d!d!d W%- D%- P%- M%- d%-DrawMode!d T%- v%-ViewAngle!d p%-ProjectionMode!d a%-R|G|B!F!F!F o%-OutName!s z%- DFiles!*s";

IRIT_STATIC_DATA const char
    *OutFileName = "irit2wgl.html";

static void Irit2WglExit(int ExitCode);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Main module of Irit2Wgl - Read command line and do what is needed...	     M
*                                                                            *
* PARAMETERS:                                                                M
*   argc, argv:  Command line.                                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    Return code.                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   main                                                                     M
*****************************************************************************/
int main(int argc, char **argv)
{
    int Error,
	FineNessFlag = FALSE,
	HideCtrlBar = FALSE,
	CanvasWidthFlag = FALSE,
	CanvasWidth = 768,
	CanvasHeightFlag = FALSE,
	CanvasHeight = 512,
	BkRGBFlag = FALSE,
	BkRGB[3] = {0, 0, 0},
	ShowWorldAxes = FALSE,
	DisableDepthTest = FALSE,
	EnablePicking = FALSE,
	ShowModelAxes = FALSE,
	DrawModeFlag = FALSE,
	ModelTrans = FALSE,
	ViewAngleFlag = FALSE,
	ProjectionModeFlag = FALSE,
	AmbientRGBFlag = FALSE,
	VerFlag = FALSE,
	OutFileFlag = FALSE,
	NumFiles = 0;
    IrtRType
	AmbientRGB[3] = {0.2, 0.2, 0.2};
    Irit2WglDrawModeType
	DrawMode = IRIT2WGL_DRAW_MODE_WIREFRAME;
    Irit2WglViewAngleType
	ViewAngle = IRIT2WGL_VIEW_ANGLE_ORIGINAL;
    Irit2WglViewAngleType
	ProjectionMode = IRIT2WGL_PROJECTION_MODE_ORTHOGONAL;
    char **FileNames;
    IPObjectStruct
	*PObjects = NULL;

#ifdef DEBUG_IRIT_MALLOC
    IritInitTestDynMemory();
#endif /* DEBUG_IRIT_MALLOC */

    if ((Error = GAGetArgs(argc, argv, CtrlStr,
			   &IPFFCState.LinearOnePolyFlag,
			   &IPFFCState.FourPerFlat,
			   &FineNessFlag,
			   &IPFFCState.OptimalPolygons, &IPFFCState.FineNess,
			   &HideCtrlBar,
			   &CanvasWidthFlag, &CanvasWidth,
			   &CanvasHeightFlag, &CanvasHeight,
			   &BkRGBFlag,
			   &BkRGB[0], &BkRGB[1], &BkRGB[2],
			   &ShowWorldAxes, &DisableDepthTest, &EnablePicking,
			   &ShowModelAxes,
			   &DrawModeFlag, &DrawMode,
			   &ModelTrans,
			   &ViewAngleFlag, &ViewAngle,
			   &ProjectionModeFlag, &ProjectionMode,
			   &AmbientRGBFlag,
			   &AmbientRGB[0], &AmbientRGB[1], &AmbientRGB[2],
			   &OutFileFlag, &OutFileName,
			   &VerFlag, &NumFiles, &FileNames)) != 0) {
	GAPrintErrMsg(Error);
	GAPrintHowTo(CtrlStr);
	Irit2WglExit(1);
    }

    if (VerFlag) {
	IRIT_INFO_MSG_PRINTF("\n%s\n\n", VersionStr);
	GAPrintHowTo(CtrlStr);
	Irit2WglExit(0);
    }

    if (!NumFiles) {
	IRIT_WARNING_MSG("No data file names were given, exit.\n");
	GAPrintHowTo(CtrlStr);
	Irit2WglExit(1);
    }

    /* Wants UV coordinates for textures. */
    IPFFCState.ComputeUV = TRUE;

    /* Get the data files: */
    IPSetFlattenObjects(FALSE);
    if ((PObjects = IPGetDataFiles((const char **) FileNames,
				   NumFiles, TRUE, FALSE)) == NULL)
	Irit2WglExit(1);
    PObjects = IPResolveInstances(PObjects);

    IPWGLSaveFile(PObjects, OutFileFlag ? OutFileName : NULL, HideCtrlBar,
		  CanvasWidth, CanvasHeight, BkRGB, ShowWorldAxes,
		  DisableDepthTest, EnablePicking, ShowModelAxes, DrawMode,
		  ModelTrans, ViewAngle, ProjectionMode, AmbientRGB);

    Irit2WglExit(0);

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Irit2Wgl exit routine.                                                    *
*                                                                            *
* PARAMETERS:                                                                *
*   ExitCode:                                                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
void Irit2WglExit(int ExitCode)
{
    exit(ExitCode);
}

#ifdef DEBUG

/*****************************************************************************
* DESCRIPTION:                                                               *
*    Dummy function to link at debugging time.                               *
*                                                                            *
* PARAMETERS:                                                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*                                                                            *
* KEYWORDS:                                                                  *
*****************************************************************************/
void DummyLinkCagdDebug(void)
{
    IPDbg();
}

#endif /* DEBUG */
