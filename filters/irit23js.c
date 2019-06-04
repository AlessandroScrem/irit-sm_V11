/*****************************************************************************
* A filter to convert IRIT data files to a Three.js viewing package.         *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber/Jesse Oberstein		Ver 2.0, July 2014   *
*****************************************************************************/

#ifdef USE_VARARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif /* USE_VARARGS */
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/attribut.h"
#include "inc_irit/allocate.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/misc_lib.h"
#include "irit23js.h"

#ifdef NO_CONCAT_STR
IRIT_STATIC_DATA const char
    *VersionStr = "irit23JS 		Version 11,	Gershon Elber\
	/Jesse Oberstein,\n\
	(C) Copyright 1989-2015 Gershon Elber, Non commercial use only.";
#else
IRIT_STATIC_DATA const char
    *VersionStr = "irit23JS	" IRIT_VERSION ",	Gershon Elber\
	/Jesse Oberstein,	" __DATE__ ", " __TIME__ "\n" 
	IRIT_COPYRIGHT ", Non commercial use only.";
#endif /* NO_CONCAT_STR */

IRIT_STATIC_DATA const char
    *CtrlStr = "irit23JS l%- 4%- p%- F%-PolyOpti|FineNess!d!F i!-InFile!s\
		o!-OutName!s T%- t%-AnimTime!F z%-";

IRIT_STATIC_DATA int
    GlblTotalPolys = 0,
    GlblAccumVertices = 0;

IRIT_STATIC_DATA FILE 
    *GlblOutputJsFile, 
    *GlblFacesFile, 
    *GlblNormalsFile, 
    *GlblColorsFile, 
    *GlblMaterialsFile, 
    *GlblUVsFile;

IRIT_STATIC_DATA char 
    GlblFacesFileName[IRIT_LINE_LEN], 
    GlblNormalsFileName[IRIT_LINE_LEN], 
    GlblColorsFileName[IRIT_LINE_LEN], 
    GlblMaterialsFileName[IRIT_LINE_LEN], 
    GlblUVsFileName[IRIT_LINE_LEN],
    *GlblUVScales;

/* Functions that provide unique data to the output js file. */
static void DumpOneTraversedObject(IPObjectStruct *PObj, IrtHmgnMatType Mat);
static int PrintVertices(IPObjectStruct *PObj);
static void PrintFacesToFile(IPObjectStruct *PObj,
			     int PreviousShapeVerts,
			     const int *HasUVs);
static void PrintNormalsToFile(IPObjectStruct *PObj);
static void PrintColorsToFile(IPObjectStruct *PObj);
static void PrintMaterialsToFile(IPObjectStruct *PObj);
static char * PrintGlblUVScales(const char *Scales,
				const char *Scale1,
				const char *Scale2);
static void PrintUVsToFile(IPObjectStruct *PObj,
			   int PreviousShapeVerts,
			   int *HasUVs);
static void ReprintInformation(FILE *File,
			       const char *FileName,
			       const char *Header);

/* Prints out the ThreeJS viewer as an HTML file, and also prints the      */
/* necessary libraries for ThreeJS to function correctly:                  */
/* irit3js.js and iritOC.js.						   */
static void PrintLibraryFile(const char *OutputFileName,
			     const char *LibStr[],
			     int Lib3JS);
static void PrintOtherFiles(const char *InputFileName,
			    const char *OutputFileName,
			    int PerspectiveCamera);

/* Miscellaneous helper functions. */
static void ReplaceForwardSlashes(char *FileName);
static char *Real2Str(IrtRType R);
static void Irit23jsExit(int ExitCode);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Main module of skeletn1 - Read command line and do what is needed...	     M
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
	HasTime = FALSE,
	SrfFineNessFlag = FALSE,
	PerspectiveCameraFlag = FALSE,
	VerFlag = FALSE,
	OutFileFlag = FALSE,
	InFileFlag = FALSE;
    char
        *InFileName = NULL,
        *OutFileName = "irit23js.js";
    IrtRType CurrentTime;
    IPObjectStruct *PObjects;
    IrtHmgnMatType CrntViewMat;

    /* Here some useful parameters to play with in tessellating freeforms:   */
    IPFFCState.FineNess = 20;/* Resolution of tessellation, larger is finer. */
    IPFFCState.ComputeUV = FALSE;      /* Wants UV coordinates for textures. */
    IPFFCState.FourPerFlat = FALSE;  /* 4 polys per ~flat patch, 2 otherwise.*/
    IPFFCState.LinearOnePolyFlag = FALSE;  /* Linear srf generates one poly. */

#ifdef DEBUG_IRIT_MALLOC
    IritInitTestDynMemory();
#endif /* DEBUG_IRIT_MALLOC */
    if ((Error = GAGetArgs(argc, argv, CtrlStr,
			   &IPFFCState.LinearOnePolyFlag,
			   &IPFFCState.FourPerFlat, &PerspectiveCameraFlag,
			   &SrfFineNessFlag, &IPFFCState.OptimalPolygons,
			   &IPFFCState.FineNess,
			   &InFileFlag, &InFileName, &OutFileFlag,
			   &OutFileName,
			   &IPFFCState.Talkative,
			   &HasTime, &CurrentTime,
			   &VerFlag)) != 0) {

	GAPrintErrMsg(Error);
	GAPrintHowTo(CtrlStr);
	Irit23jsExit(1);
    }
    if (VerFlag) {
	IRIT_INFO_MSG_PRINTF("\n%s\n\n", VersionStr);
	GAPrintHowTo(CtrlStr);
	Irit23jsExit(0);
    }

    if (IPFFCState.Talkative)
	IRIT_INFO_MSG_PRINTF("%s triangles per flat will be created.\n",
		             IPFFCState.FourPerFlat ? "Four" : "Two");
    
    /* Append a .itd extension to the output file's name if it is missing. */
    if (strstr(InFileName, ".itd") == NULL) {
	char WithITDExtension[IRIT_LINE_LEN];

	strcpy(WithITDExtension, InFileName);
	InFileName = strcat(WithITDExtension, ".itd");
    }

    /* Append a .js extension to the output file's name if it is missing. */
    if (strstr(OutFileName, ".js") == NULL) {
	char WithJSExtension[IRIT_LINE_LEN];

	strcpy(WithJSExtension, OutFileName);
	OutFileName = strcat(WithJSExtension, ".js");
    }

    ReplaceForwardSlashes(InFileName);
    ReplaceForwardSlashes(OutFileName);

    /* Get the data files: */
    IPSetFlattenObjects(FALSE);
    if ((PObjects = IPGetDataFiles((const char **) &InFileName, 
				    1, TRUE, FALSE)) == NULL) {
	Irit23jsExit(1);
    }
    PObjects = IPResolveInstances(PObjects);

    if (IPWasPrspMat)
	MatMultTwo4by4(CrntViewMat, IPViewMat, IPPrspMat);
    else
	IRIT_GEN_COPY(CrntViewMat, IPViewMat, sizeof(IrtHmgnMatType));

    IPTraverseObjListHierarchy(PObjects, CrntViewMat, IPMapObjectInPlace);

    strcpy(GlblMaterialsFileName, OutFileName);
    GlblMaterialsFile = 
	fopen(strcat(GlblMaterialsFileName, "_materials.txt"), "w");
    fclose(GlblMaterialsFile);

    strcpy(GlblNormalsFileName, OutFileName);
    GlblNormalsFile = fopen(strcat(GlblNormalsFileName, "_normals.txt"), "w");
    fclose(GlblNormalsFile);

    strcpy(GlblColorsFileName, OutFileName);
    GlblColorsFile = fopen(strcat(GlblColorsFileName, "_colors.txt"), "w");
    fclose(GlblColorsFile);

    strcpy(GlblUVsFileName, OutFileName);
    GlblUVsFile = fopen(strcat(GlblUVsFileName, "_uvs.txt"), "w");
    fclose(GlblUVsFile);
    
    strcpy(GlblFacesFileName, OutFileName);
    if ((GlblFacesFile = fopen(strcat(GlblFacesFileName,
				      "_faces.txt"), "w")) == NULL) {
        fprintf(stderr,
		"irit23js: Failed to open file \"%s\" for writing.\n",
		GlblFacesFileName);
	return 1;
    }
    fclose(GlblFacesFile);
    
    /* Print out all the model's polygon data to the specified output file. */
    GlblOutputJsFile = fopen(OutFileName, "w");

    /* Print the header of the output js file. */
    fprintf(GlblOutputJsFile, "{\n\"metadata\": {}, \n\"scale\": 1.0, \n");
    IPTraverseObjListHierarchy(PObjects, CrntViewMat, DumpOneTraversedObject);
    /* Close the vertices array. */
    fprintf(GlblOutputJsFile, "\n],\n");
    ReprintInformation(GlblMaterialsFile, GlblMaterialsFileName, "materials");
    ReprintInformation(GlblNormalsFile, GlblNormalsFileName, "normals");
    ReprintInformation(GlblColorsFile, GlblColorsFileName, "colors");
    ReprintInformation(GlblUVsFile, GlblUVsFileName, "uvs");
    ReprintInformation(GlblFacesFile, GlblFacesFileName, "faces");

    /* Print the footer of the output js file. */
    fprintf(GlblOutputJsFile, "\"morphTargets\": [] \n} \n");
    fclose(GlblOutputJsFile);

    PrintOtherFiles(InFileName, OutFileName, PerspectiveCameraFlag);

    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function of IPTraverseObjListHierarchy. Called on every non    *
* list object found in hierarchy.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Non list object to handle.                                   *
*   Mat:        Transformation matrix to apply to this object.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpOneTraversedObject(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    int PreviousShapeVerts = 0;
    IPObjectStruct *PObjs;

    if (IP_IS_FFGEOM_OBJ(PObj))
        PObjs = IPConvertFreeForm(PObj, &IPFFCState);  /* Convert in place. */
    else {
	PObjs = PObj;
    }

    for (PObj = PObjs; PObj != NULL; PObj = PObj -> Pnext) {
	if (PObj -> ObjType == IP_OBJ_POLY) {
	    int *HasUVs;
	    
	    PreviousShapeVerts = GlblAccumVertices;
	    GlblAccumVertices = PrintVertices(PObj);
	    /* Allocate memory for checking if a vertex has UVs. */
	    HasUVs = (int *) malloc(sizeof(int) * GlblAccumVertices);
	    PrintUVsToFile(PObj, PreviousShapeVerts, HasUVs);
	    PrintFacesToFile(PObj, PreviousShapeVerts, HasUVs);
	    free(HasUVs);
	    PrintNormalsToFile(PObj);
	    PrintColorsToFile(PObj);
	    PrintMaterialsToFile(PObj);
	    GlblTotalPolys++;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prints the vertex data from given geometry object.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Object to print.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   The current number of accumulated vertices in the object.                *
*****************************************************************************/
static int PrintVertices(IPObjectStruct *PObj) 
{
    char VertLine[IRIT_LINE_LEN];
    IPVertexStruct **V;
    IPPolyVrtxIdxStruct *PVIdx;

    /* Dump vertices so that each identical vertex shows up once, Vrml style.*/
    PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 7);
    
    if (GlblAccumVertices == 0) {
	fprintf(GlblOutputJsFile, "\"vertices\": [ \n");
    }
    for (V = PVIdx -> Vertices; *V != NULL; V++) {
	sprintf(VertLine, "%s,%s,%s\n",
		Real2Str((*V) -> Coord[0]),
		Real2Str((*V) -> Coord[1]),
		Real2Str((*V) -> Coord[2]));

	GlblAccumVertices == 0 ?
	    fprintf(GlblOutputJsFile, "%s%s", "\t ", VertLine) :
	    fprintf(GlblOutputJsFile, "%s%s", "\t,", VertLine);

	GlblAccumVertices++;
    }    
    IPPolyVrtxIdxFree(PVIdx);
    return GlblAccumVertices;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prints the index data from given geometry object for all of its othe     *
* properties.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Object to print.                                             *
*   PreviousShapeVerts:     The amount of vertices in the previous polygon.  *
*   HasUVs:     An array holding values that determine if a vertex has UVs.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintFacesToFile(IPObjectStruct *PObj,
			     int PreviousShapeVerts,
			     const int *HasUVs)
{
    IPPolyVrtxIdxStruct *PVIdx;
    char Index[IRIT_LINE_LEN_SHORT], MaterIndx[IRIT_LINE_LEN_SHORT],
	VertIndxs[IRIT_LINE_LEN], ColorIndxs[IRIT_LINE_LEN],
	NormIndxs[IRIT_LINE_LEN], UVIndxs[IRIT_LINE_LEN];
    int **Pls, FaceVerts, Bitmask,
	ContainsUVs = FALSE;

    PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 7);
	
    for (Pls = PVIdx -> Polygons; *Pls != NULL; Pls++) {
        int *Pl = *Pls;
	int currentIndex;
	Bitmask = 0;
	FaceVerts = 0;

	IRIT_ZAP_MEM(VertIndxs, IRIT_LINE_LEN);
	IRIT_ZAP_MEM(NormIndxs, IRIT_LINE_LEN);
	IRIT_ZAP_MEM(ColorIndxs, IRIT_LINE_LEN);
	IRIT_ZAP_MEM(UVIndxs, IRIT_LINE_LEN);

	/* Assume at least one edge in polygon! */
	do {
	    *Pl += PreviousShapeVerts;
	    currentIndex = *Pl;
	    sprintf(Index, ",%d", *Pl++);
	    sprintf(VertIndxs, strcat(VertIndxs, Index));
	    sprintf(NormIndxs, strcat(NormIndxs, Index));
	    sprintf(ColorIndxs, strcat(ColorIndxs, Index));

	    if (HasUVs[currentIndex] == TRUE) {
		sprintf(UVIndxs, strcat(UVIndxs, Index));
		ContainsUVs = TRUE;
	    }
	    else {
		ContainsUVs = FALSE;
	    }
	    FaceVerts++;
	}
	while (*Pl >= 0);
		
	/* If the face has four vertices, add 1 to the bitmask.  If it has 
	only three, nothing is added to the bitmask. */
	if (FaceVerts == 4) {
	    Bitmask += 1;
	}

	/* Add 2 to the bitmask, and add the material index. */
	sprintf(MaterIndx, ",%d", GlblTotalPolys);
	strcat(VertIndxs, MaterIndx);
	Bitmask += 2;

	/* If the model has UVs, append the same vertex indicies time to 
	account for UV vertex indicies. */
	if (ContainsUVs == TRUE) {
	    strcat(VertIndxs, UVIndxs);
	    Bitmask += 8;
	}

	/* Append the same vertex indicies a second time to account for 
	normals indicies. */
	strcat(VertIndxs, NormIndxs);
	Bitmask += 32;

	/* Append the same vertex indicies a third time to account
	for color vertex indicies, and move to a new line. */
	strcat(VertIndxs, strcat(ColorIndxs, "\n"));
	Bitmask += 128;

	/* Print all the indices for the face, along with the bitmask prepended
	to the beginning of the line. */

	if ((GlblFacesFile = fopen(GlblFacesFileName, "a")) == NULL) {
	    fprintf(stderr,
		    "irit23js: Failed to open file \"%s\" for writing.\n",
		    GlblFacesFileName);
	    return;
	}

	fprintf(GlblFacesFile, " ,%s%s", Real2Str(Bitmask), VertIndxs);

	fclose(GlblFacesFile);
    }
    IPPolyVrtxIdxFree(PVIdx);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prints the normal data from given geometry object.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Object to print.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintNormalsToFile(IPObjectStruct *PObj) 
{
    IPVertexStruct **V;
    IPPolyVrtxIdxStruct *PVIdx;

    PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 7);

    if ((GlblNormalsFile = fopen(GlblNormalsFileName, "a")) == NULL) {
        fprintf(stderr, "irit23js: Failed to open file \"%s\" for writing.\n",
		GlblNormalsFileName);
        return;
    }

    for (V = PVIdx -> Vertices; *V != NULL; V++) {
	    fprintf(GlblNormalsFile, " ,%s,%s,%s\n",
		Real2Str((*V) -> Normal[0]),
		Real2Str((*V) -> Normal[1]),
		Real2Str((*V) -> Normal[2]));
    }
    fclose(GlblNormalsFile);
    IPPolyVrtxIdxFree(PVIdx);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prints the color data from given geometry object.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Object to print.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintColorsToFile(IPObjectStruct *PObj) 
{
    IPVertexStruct **V;
    IPPolyVrtxIdxStruct *PVIdx;
    int RGB[3];

    PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 7);

    if ((GlblColorsFile = fopen(GlblColorsFileName, "a")) == NULL) {
        fprintf(stderr, "irit23js: Failed to open file \"%s\" for writing.\n",
		GlblColorsFileName);
        return;
    }

    for (V = PVIdx -> Vertices; *V != NULL; V++) {
	if (AttrGetRGBColor2((*V) -> Attr, NULL,
			     &RGB[0], &RGB[1], &RGB[2]) == TRUE) {
		AttrGetRGBColor2((*V) -> Attr, NULL, 
				 &RGB[0], &RGB[1], &RGB[2]);
	}
	else if (AttrGetRGBColor2((PObj) -> Attr, NULL,
				  &RGB[0], &RGB[1], &RGB[2]) == TRUE) {
		AttrGetRGBColor2((PObj) -> Attr, NULL, 
				 &RGB[0], &RGB[1], &RGB[2]);
	}
	else {
		RGB[0] = 255;
		RGB[1] = 255;
		RGB[2] = 255;
	}
	fprintf(GlblColorsFile, " ,%d\n",
		((RGB[0] * 256 * 256) + (RGB[1] * 256) + RGB[2]));
    }
    fclose(GlblColorsFile);
    IPPolyVrtxIdxFree(PVIdx);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prints the material and texture data from given geometry object.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Object to print.                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintMaterialsToFile(IPObjectStruct *PObj) 
{
    int PTextureFlip = FALSE,
        NewImage = FALSE;
    const char
        *PTexture = NULL;
    char *p, PTextureFileName[IRIT_LINE_LEN_LONG];
    IrtRType Transp, PTextureSUV[3];
    IrtUVType SUV = { 1.0, 1.0 };

    Transp = AttrGetRealAttrib(PObj -> Attr, "transp");
   
    if ((GlblMaterialsFile = fopen(GlblMaterialsFileName, "a")) == NULL) {
        fprintf(stderr, "irit23js: Failed to open file \"%s\" for writing.\n",
		GlblMaterialsFileName);
        return;
    }

    PTexture = AttrGetStrAttrib(PObj -> Attr, "ptexture");
    fprintf(GlblMaterialsFile, " ,{ \n\
	    \"DbgColor\" : 15658734, \n\
	    \"DbgIndex\" : 0, \n\
	    \"DbgName\" : \"default\", \n");

    if (PTexture != '\0') {
	IrtImgParsePTextureString(PTexture, PTextureFileName, PTextureSUV,
				  &PTextureFlip, &NewImage);
	if ((p = strrchr(PTextureFileName, '.')) != NULL) {
	    *p = 0;
	}
	else {
	    p = PTextureFileName;
	}
	strcat(p, ".jpg");
	fprintf(GlblMaterialsFile, "\t\t\"mapDiffuse\": \"%s\", \n", 
	    PTextureFileName);
	SUV[0] = PTextureSUV[0];
	SUV[1] = PTextureSUV[1];
	GlblUVScales = PrintGlblUVScales(GlblUVScales,
					 Real2Str(SUV[0]), Real2Str(SUV[1]));
    }
    if (Transp != IP_ATTR_BAD_REAL) {
	fprintf(GlblMaterialsFile, "\t\t\"transparent\" : true,\n\
			      \"transparency\": %s, \n", Real2Str(1 - Transp));
    }
    fprintf(GlblMaterialsFile, "\t\t\"vertexColors\" : \"true\" \n\t} \n");
    fclose(GlblMaterialsFile);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prints the scale data for the UVs from a texture if it exists, from      *
* given geometry object.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   Scales:       The current string of all the scale data until now.        *
*   Scale1:       The horizontal UV scale value.			     *
*   Scale2:       The vertical UV scale value.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   The resulting string containing all the scale data so far, including the *
* newest added pair.							     *
*****************************************************************************/
static char *PrintGlblUVScales(const char *Scales,
			       const char *Scale1,
			       const char *Scale2)
{
    IRIT_STATIC_DATA char ScaleBuffer[IRIT_LINE_LEN_VLONG];

    sprintf(ScaleBuffer,
	    "%smaterial.materials[%d].map.repeat.set(%s, %s);\n\t\t\t\t", 
	    GlblTotalPolys == 0 ? "" : Scales,
	    GlblTotalPolys, Scale1, Scale2);

    return ScaleBuffer;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prints the UV data from given geometry object.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Object to print.                                             *
*   PreviousShapeVerts:     The amount of vertices in the previous polygon.  *
*   HasUVs:     An array holding values that determine if a vertex has UVs.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintUVsToFile(IPObjectStruct *PObj,
			   int PreviousShapeVerts,
			   int *HasUVs) 
{
    IPVertexStruct **V;
    IPPolyVrtxIdxStruct *PVIdx;
    int vertexCounter = PreviousShapeVerts;

    PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 7);

    for (V = PVIdx -> Vertices; vertexCounter < GlblAccumVertices; V++) {
	float *UVs;
		
	if ((GlblUVsFile = fopen(GlblUVsFileName, "a")) == NULL) {
	    fprintf(stderr,
		    "irit23js: Failed to open file \"%s\" for writing.\n",
		    GlblUVsFileName);
	    return;
	}

	if ((UVs = AttrGetUVAttrib((*V) -> Attr, "uvvals")) != NULL) {
	    fprintf(GlblUVsFile, " ,%s, %s\n",
		    Real2Str(UVs[0]), Real2Str(UVs[1]));
	    HasUVs[vertexCounter] = TRUE;
	}
	else {
	    fprintf(GlblUVsFile, " ,[]\n");
	    HasUVs[vertexCounter] = FALSE;
	}
	fclose(GlblUVsFile);
	vertexCounter++;
    }

    IPPolyVrtxIdxFree(PVIdx);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*    Retrieves the information that was printed to file for a property	     *
*  identified by the given header.                                           *
*									     *
* PARAMETERS:                                                                *
*   File:	  Pointer to the file the data is to be read from.           *
*   FileName:     Name of the file the data is to be read from.              *
*                 This file is removed as a side effect of this function.    *
*   Header:	  Name of the data being printed (Three.js-specific JSON     *
*		  header).						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void ReprintInformation(FILE *File,
			       const char *FileName,
			       const char *Header) 
{
    char InfoLine[IRIT_LINE_LEN_VLONG];
    int LineCounter = 0;

    fprintf(GlblOutputJsFile, "\"");
    fprintf(GlblOutputJsFile, Header);

    /* UV property is formatted with an extra bracket */
    if (strcmp(Header, "uvs") == 0)
	fprintf(GlblOutputJsFile, "\": [ [\n\t");
    else
	fprintf(GlblOutputJsFile, "\": [ \n\t");

    if ((File = fopen(FileName, "r")) == NULL) {
        fprintf(stderr, "irit23js: Failed to open file \"%s\" for reading.\n",
		FileName);
        return;
    }

    while (fgetc(File) != EOF) {
        fgets(InfoLine, IRIT_LINE_LEN_VLONG, File);
	if (LineCounter++ == 0)
	    InfoLine[0] = ' ';
	fprintf(GlblOutputJsFile, "%s\t", InfoLine);
    }

    fclose(File);

    /* UV property is formatted with an extra bracket */
    if (strcmp(Header, "uvs") == 0)
	fprintf(GlblOutputJsFile, "\n] ],\n");
    else
	fprintf(GlblOutputJsFile, "\n],\n");

    remove(FileName);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*    Prints the Three.js library file, as well as the OrbitControls.js       *
*  library file from strings embedded in separate C source files.            *
*									     *
* PARAMETERS:                                                                *
*   OutputFileName:     Name of the output library file.		     *
*   LibStr:		The library as strings, one per line.		     *
*   Lib3JS:		TRUE for printing 3js library, FALSE for OC library. *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintLibraryFile(const char *OutputFileName,
			     const char *LibStr[],
			     int Lib3JS)
{
    char *fn, FileName[IRIT_LINE_LEN_LONG];
    int i;
    FILE *File;

    strcpy(FileName, OutputFileName);
    if ((fn = strrchr(FileName, '\\')) != NULL) {
	*fn = 0;
	strcat(fn, Lib3JS ? "\\irit3js.js" : "\\iritOC.js");
    }
    else
        strcpy(FileName, Lib3JS ? "irit3js.js" : "iritOC.js");

    if ((File = fopen(FileName, "w")) == NULL) {
        fprintf(stderr, "irit23js: Failed to open file \"%s\" for writing.\n",
		FileName);
        return;
    }

    for (i = 0; LibStr[i] != NULL; i++) {
        fprintf(File, "%s\n", LibStr[i]);
    }
    fclose(File);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*    Prints irit3js.js, iritOC.js, and a viewer HTML file.		     *
*									     *
* PARAMETERS:                                                                *
*   InputFileName: Pointer to the file that the user gave as a command line  *
*		   argument with flag -i.				     *
*   OutputFileName: Pointer to the file that the user gave as a command line *
*		   argument with flag -o.				     *
*   PerspectiveCamera:  TRUE to create a persective camera, FALSE for an     *
*                  orthographic one.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrintOtherFiles(const char *InputFileName,
			    const char *OutputFileName,
			    int PerspectiveCamera)
{
    char HtmlFileName[IRIT_LINE_LEN_LONG], *h,
	FileOutputName[IRIT_LINE_LEN_LONG];
    const char *Camera, *MouseWheelEvent, *Zoom, *ResetCamera, *i, *o;
    FILE *HtmlFile;
    
    /* Print the Three.js library file */
    PrintLibraryFile(OutputFileName, (const char **) Irit23jsGlbl3JSString,
		     TRUE);

    /* Print the OrbitControls.js library file. */
    PrintLibraryFile(OutputFileName, (const char **) Irit23jsGlblOCString,
		     FALSE);

    /* Print the viewer HTML file and adjust file names accordinly. */
    strcpy(HtmlFileName, OutputFileName);
    strcpy(FileOutputName, OutputFileName);
    if ((h = strrchr(HtmlFileName, '.')) != NULL)
	*h = 0;

    if ((o = strrchr(OutputFileName, '\\')) != NULL)
	o = strrchr(OutputFileName, '\\') + 1;
    else
	o = FileOutputName;

    if ((i = strrchr(InputFileName, '\\')) != NULL)
	i = strrchr(InputFileName, '\\') + 1;
    else
	i = InputFileName;

    HtmlFile = fopen(strcat(HtmlFileName, "Viewer.html"), "w");

    if (PerspectiveCamera) {
        Camera = "camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 1, 1000);";
	MouseWheelEvent = "";
	Zoom = "";
	ResetCamera = "controls.reset();";
    }
    else {
        Camera = "camera = new THREE.OrthographicCamera(-2 * zoom, 2 * zoom, 2 * zoom * window.innerHeight / window.innerWidth, -2 * zoom * window.innerHeight / window.innerWidth, 1, 1000);";
	MouseWheelEvent = "renderer.domElement.addEventListener( 'mousewheel', mousewheel, false );\n\
			    renderer.domElement.addEventListener( 'DOMMouseScroll', mousewheel, false ); // firefox";
        Zoom = "// Uses the mousewheel event to control zoom for an orthographic camera.\n\
		function mousewheel( event ) {\n\
			event.preventDefault();\n\
			event.stopPropagation();\n\
			\n\
			var delta = 0;\n\
			\n\
			if ( event.wheelDelta ) { // WebKit / Opera / Explorer 9\n\
				delta = event.wheelDelta / 40;\n\
			} else if ( event.detail ) { // Firefox\n\
				delta = - event.detail / 3;\n\
			}\n\
			\n\
			var width = camera.right / zoom;\n\
			var height = camera.top / zoom;\n\
			\n\
			zoom -= delta * 0.05;\n\
			\n\
			camera.left = -zoom*width;\n\
			camera.right = zoom*width;\n\
			camera.top = zoom*height;\n\
			camera.bottom = -zoom*height;\n\
			\n\
			camera.updateProjectionMatrix();\n\
			\n\
			renderer.render( scene, camera );\n\
		}";
	ResetCamera = "zoom = 1;\n\
			camera = new THREE.OrthographicCamera(-2 * zoom, 2 * zoom, 2 * zoom * window.innerHeight / window.innerWidth, -2 * zoom * window.innerHeight / window.innerWidth, 1, 1000);\n\
			camera.position.x = 0;\n\
			camera.position.y = 0;\n\
			camera.position.z = 7;\n\
			\n\
			controls = new THREE.OrbitControls(camera, renderer.domElement);\n\
			controls.reset();";
    }

    fprintf(HtmlFile, Irit23jsGlblViewerString, 
	    PerspectiveCamera ? "Perspective" : "Orthographic", 
	    i, Camera, MouseWheelEvent, Zoom, ResetCamera, 
	    GlblUVScales == NULL ? "" : GlblUVScales, o);

    fclose(HtmlFile);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Replaces all the forward slashes in a filepath with backslashes.         *
*                                                                            *
* PARAMETERS:                                                                *
*   FileName:       The name of the file to examine for forward slashes.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void		                                                     *
*****************************************************************************/
static void ReplaceForwardSlashes(char *FileName) {
    size_t i;

    for (i = 0; i < strlen(FileName); i++) {
	if (FileName[i] == '/') {
	    FileName[i] = '\\';
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Converts a real number into a string.					     *
*   The routine maintains 6 different buffers simultanuously so up to 6      *
* calls can be issued from same printf...				     *
*                                                                            *
* PARAMETERS:                                                                *
*   R:          To convert to a string.                                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   char *:     A string representation for R.                               *
*****************************************************************************/
static char *Real2Str(IrtRType R)
{
    IRIT_STATIC_DATA int j, k,
	i = 0;
    IRIT_STATIC_DATA char Buffer[6][IRIT_LINE_LEN_SHORT];

    if (IRIT_FABS(R) < IRIT_EPS)
	R = 0.0;			    /* Round off very small numbers. */

    sprintf(Buffer[i], "%6g", R);

    for (k = 0; !isdigit(Buffer[i][k]) && k < IRIT_LINE_LEN; k++);
    if (k >= IRIT_LINE_LEN) {
	IRIT_WARNING_MSG_PRINTF("Conversion of real number (%f) failed.\n", R);
	Irit23jsExit(1);
    }

    for (j = (int) (strlen(Buffer[i]) - 1); Buffer[i][j] == ' ' && j > k; j--);
    if (strchr(Buffer[i], '.') != NULL)
	for (; Buffer[i][j] == '0' && j > k; j--);
    Buffer[i][j+1] = 0;

    j = i;
    i = (i + 1) % 6;
    return Buffer[j];
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Irit23jsExit exit routine.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   ExitCode:                                                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void Irit23jsExit(int ExitCode)
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

#ifdef IRIT23JS_STRING_FORMATTER

/*****************************************************************************
* DESCRIPTION:								     *
*   Used to convert a specified input file (in function as InputLib) to an   *
* output file consisting of an array of strings, terminating with NULL.      *
*									     *
*   Before converting a library, the following characters need to be	     *
* unescaped within the original InputLib.				     *
*     - % with %%							     *
*     - \ with \\							     *
*     - " with \"							     *
*									     *
*   NUM_LINES for OrbitControls.js is 614 (default is 38491 for Three.js)    *
*									     *
* PARAMETERS:                                                                *
*   NONE                                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
void PrintFormattedLibrary()
{
    FILE *OutputLib, *InputLib;
    static char Line[10000];
    int LineCounter,
	NUM_LINES = 38491,
	NSpecifier = _set_printf_count_output(1);

    /* Specify custom input/output file names here */
    InputLib = fopen("C:\\users\\jesse\\desktop\\three_copy.js", "r");
    OutputLib = fopen("C:\\users\\jesse\\desktop\\three_formatted.js", "w");
	for (LineCounter = 0; LineCounter < NUM_LINES; LineCounter++) {
		fgets(Line, sizeof(Line), InputLib);
		Line[strlen(Line) - 1] = ' ';
		fprintf(OutputLib, "\"");
		fprintf(OutputLib, strcat(strcat(strcat(Line, " \\"), "n"), "\","));
		fprintf(OutputLib, "\n");
    }
    fprintf(OutputLib, "NULL");
    fclose(InputLib);
    fclose(OutputLib);
    NSpecifier = _set_printf_count_output( 0 );
}

#endif /* IRIT23JS_STRING_FORMATTER */
