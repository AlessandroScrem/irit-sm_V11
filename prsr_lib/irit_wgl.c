/******************************************************************************
* Filter to convert IRIT data files to a WebGL HTML file.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by:  Avi Kaplan				Ver 1.0, Feb 2012     *
******************************************************************************/

/* There seems to be a bug with GL.LINE_LOOP use in drawElements. Until this */
/* bug is fixed, keep this define. The bug recreation can be found in        */
/* 'draw_elements_line_loop_bug.zip'. */
#define NASTY_PATCH_FOR_DRAW_ELEMENTS_LINE_LOOP_BUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/grap_lib.h"

#ifdef USE_VARARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif /* USE_VARARGS */

/* Maximum length of buffer. */
#define IRIT2WGL_MAX_BUFFER_LEN 1024

/* Maximum number of light sources. */
#define IRIT2WGL_MAX_LIGHT_SOURCE_NUM 8

/* Files used by this filter. */
#define IRIT2WGL_HTML_SKELETON_FILE_NAME "irit2wgl_skeleton.html"
#define IRIT2WGL_CSS_SKELETON_FILE_NAME "irit2wgl_skeleton.css"
#define IRIT2WGL_JS_SKELETON_FILE_NAME "irit2wgl_skeleton.js"

/* Object tag type. */
typedef enum {
    IRIT2WGL_MODEL_OBJECT_TAG_POLYLINE = 0,
    IRIT2WGL_MODEL_OBJECT_TAG_POLYGON
} Irit2WglModelObjectTagType;

/* Vertex struct. */
typedef struct Irit2WglVertexStruct {
    IrtPtType Pos;
    IrtNrmlType Normal;
    IrtUVType UV;
} Irit2WglVertexStruct;

/* Poly struct. */
typedef struct Irit2WglPolyStruct {
    struct Irit2WglPolyStruct *Pnext;

    int *Indices;
    int IndexCount;
} Irit2WglPolyStruct;

/* Texture struct. */
typedef struct Irit2WglTextureStruct {
    struct Irit2WglTextureStruct *Pnext;

    char *FileName;
    int Index;
} Irit2WglTextureStruct;

/* Model Object struct. */
typedef struct Irit2WglModelObjectStruct {
    struct Irit2WglModelObjectStruct *Pnext;

    char *Name;
    Irit2WglModelObjectTagType Tag;
    IrtRType RGBA[4];
    Irit2WglTextureStruct *Texture;
    IrtUVType SUV, MinUV, MaxUV;

    Irit2WglVertexStruct *Vertices;
    int VertexCount;

    Irit2WglPolyStruct *Polys;
    Irit2WglPolyStruct *LastPoly;
} Irit2WglModelObjectStruct;

/* Light source struct. */
typedef struct Irit2WglLightSourceStruct {
    IGLightType Pos;
    IrtVecType Diffuse;
    IrtRType DiffuseI;
    IrtVecType Specular;
    IrtRType SpecularI;
} Irit2WglLightSourceStruct;

/* Scene struct. */
typedef struct Irit2WglSceneStruct {
    Irit2WglModelObjectStruct *Objects;
    Irit2WglModelObjectStruct *LastObject;

    Irit2WglTextureStruct *Textures;
    Irit2WglTextureStruct *LastTexture;

    Irit2WglLightSourceStruct Lights[IRIT2WGL_MAX_LIGHT_SOURCE_NUM];
    int LightSourceSetCount;

    int NOP;
} Irit2WglSceneStruct;

IRIT_STATIC_DATA const char
    *GlblIritBinDir = NULL;

IRIT_STATIC_DATA int
    GlblHideCtrlBar,
    GlblCanvasWidth,
    GlblCanvasHeight,
    GlblBkRGB[3],
    GlblShowWorldAxes,
    GlblDisableDepthTest,
    GlblEnablePicking,
    GlblShowModelAxes,
    GlblDrawMode,
    GlblModelTrans,
    GlblViewAngle,
    GlblProjectionMode;

IRIT_STATIC_DATA IrtRType
    GlblAmbientRGB[3];

IRIT_STATIC_DATA IGLightType
    GlblLightSourceInitialPos[IRIT2WGL_MAX_LIGHT_SOURCE_NUM] = {
	{1.0, 1.0, 1.0, 1.0},
	{-1.0, -1.0, -1.0, 1.0},
	{-1.0, 1.0, 1.0, 1.0},
	{1.0, -1.0, -1.0, 1.0},
	{1.0, -1.0, 1.0, 1.0},
	{-1.0, 1.0, -1.0, 1.0},
	{1.0, 1.0, -1.0, 1.0},
	{-1.0, -1.0, 1.0, 1.0}
    };

IRIT_STATIC_DATA Irit2WglSceneStruct
    *GlblScene = NULL;

static void Irit2WglProcessOneTraversedObject(IPObjectStruct *PObj,
					      IrtHmgnMatType Mat);
static void Irit2WglProcessObject(IPObjectStruct *PObj);
static Irit2WglSceneStruct *Irit2WglNewScene();
static void Irit2WglFreeScene(Irit2WglSceneStruct *Scene);
static Irit2WglModelObjectStruct *Irit2WglNewModelObject(Irit2WglSceneStruct
							                *Scene,
							 const char *Name,
							 Irit2WglModelObjectTagType
							        ModelObjectTag,
							 IrtRType RGBA[4],
							 int MaxObjVertices,
							 const char
							     *PTextureFileName,
							 IrtUVType SUV);
static void Irit2WglAddVertex(Irit2WglSceneStruct *Scene,
			      Irit2WglModelObjectStruct *ModelObject,
			      IrtPtType Pos,
			      IrtNrmlType Normal,
			      IrtUVType UV);
static void Irit2WglAddPoly(Irit2WglSceneStruct *Scene,
			    Irit2WglModelObjectStruct *ModelObject,
			    int *PolyIndexArray);
static void Irit2WglSetLightSource(Irit2WglSceneStruct *Scene,
				   IGLightType LightPos,
				   IrtVecType LightColor);
static void Irit2WglDumpData(Irit2WglSceneStruct *Scene,
			     const char* OutputFileName);
static void Irit2WglDumpScripts(Irit2WglSceneStruct *Scene,
				FILE* HtmlOutputFile,
				const char* OutputFileName);
static void Irit2WglDumpJS(Irit2WglSceneStruct *Scene,
			   const char *OutputFileName);
static void Irit2WglDumpJSSetTextureData(Irit2WglSceneStruct *Scene,
					 FILE *JSOutputFile);
static void Irit2WglDumpJSSetModelData(Irit2WglSceneStruct *Scene,
				       FILE *JSOutputFile);
static void Irit2WglDumpJSInitCameraMatrices(Irit2WglSceneStruct *Scene,
					     FILE *JSOutputFile);
static void Irit2WglDumpJSSetLightSources(Irit2WglSceneStruct *Scene,
					  FILE *JSOutputFile);
static void Irit2WglDumpJSSetControlBarParams(Irit2WglSceneStruct *Scene,
					      FILE *JSOutputFile);
static void Irit2WglDumpJSSetMenuParams(Irit2WglSceneStruct *Scene,
					FILE *JSOutputFile);
static void Irit2WglDumpJSSetStatusLogParams(Irit2WglSceneStruct *Scene,
					     FILE *JSOutputFile);
static void Irit2WglDumpCSS(Irit2WglSceneStruct *Scene,
			    const char *OutputFileName);
static const char *Irit2WglViewAngleType2Str(Irit2WglViewAngleType ViewAngle);
static const char *Irit2WglProjectionModeType2Str(Irit2WglProjectionModeType
						              ProjectionMode);
static const char *Irit2WglModelObjectTagType2Str(Irit2WglModelObjectTagType
						              ModelObjectTag);
static void Irit2WglFatalError(const char *FatalErrorMsg);
#ifdef USE_VARARGS
static const char *Irit2WglFBuffer(char *va_alist, ...);
#else
static const char *Irit2WglFBuffer(char *Format, ...);
#endif /* USE_VARARGS */

/*****************************************************************************
* DESCRIPTION:                                                               *
* Dumps the data for html into WGLFileName (stdout in NULL).		     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:               To dump into file.                                   *
*   WGLFileName:        Where output should go to, "-" or NULL for stdout.   *
*   HideCtrlBar:        Hide the scene control bar?                          *
*   CanvasWidth:        Width of the canvas in pixels.                       * 
*   CanvasHeight:       Height of the canvas in pixels.                      *
*   BkRGB:              Background color.                                    *
*   ShowWorldAxes:      Show axes relative to world?                         *
*   DisableDepthTest:   Disable depth test?                                  *
*   EnablePicking:      Enable picking?                                      *
*   ShowModelAxes:      Show axes relative to model?                         *
*   DrawMode:           Draw mode.                                           *
*   ModelTrans:         Transform in model coordinates?                      *
*   ViewAngle:          View angle of the camera.                            *
*   ProjectionMode:     Projection mode.                                     *
*   AmbientRGB:         Global ambient intensity.                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:     TRUE if successful, FALSE otherwise.                            *
*****************************************************************************/
int IPWGLSaveFile(IPObjectStruct *PObj,
		  const char *WGLFileName,
		  int HideCtrlBar,
		  int CanvasWidth,
		  int CanvasHeight,
		  int BkRGB[3],
		  int ShowWorldAxes,
		  int DisableDepthTest,
		  int EnablePicking,
		  int ShowModelAxes,
		  Irit2WglDrawModeType DrawMode,
		  int ModelTrans,
		  Irit2WglViewAngleType ViewAngle,
		  Irit2WglProjectionModeType ProjectionMode,
		  IrtRType AmbientRGB[3])
{
    char *Path;
    IrtHmgnMatType CrntViewMat;

    /* Retrieve IRIT_PATH environment variable */
    if ((Path = getenv("IRIT_PATH")) != NULL) {
	/* Remove trailing '/' if has any. */
	Path = IritStrdup(Path);
	if (Path[strlen(Path) - 1] == '/' || Path[strlen(Path) - 1] == '\\')
	    Path[strlen(Path) - 1] = 0;
	GlblIritBinDir = Path;
    }

    /* Set globals. */
    GlblHideCtrlBar = HideCtrlBar;
    GlblCanvasWidth = CanvasWidth;
    GlblCanvasHeight = CanvasHeight;
    GlblBkRGB[0] = BkRGB[0];
    GlblBkRGB[1] = BkRGB[1];
    GlblBkRGB[2] = BkRGB[2];
    GlblShowWorldAxes = ShowWorldAxes;
    GlblDisableDepthTest = DisableDepthTest;
    GlblEnablePicking = EnablePicking;
    GlblShowModelAxes = ShowModelAxes;
    GlblDrawMode = DrawMode;
    GlblModelTrans = ModelTrans;
    GlblViewAngle = ViewAngle;
    GlblProjectionMode = ProjectionMode;
    GlblAmbientRGB[0] = AmbientRGB[0];
    GlblAmbientRGB[1] = AmbientRGB[1];
    GlblAmbientRGB[2] = AmbientRGB[2];

    if (IPWasPrspMat)
	MatMultTwo4by4(CrntViewMat, IPViewMat, IPPrspMat);
    else
	IRIT_GEN_COPY(CrntViewMat, IPViewMat, sizeof(IrtHmgnMatType));

    /* Create new scene. */
    GlblScene = Irit2WglNewScene();

    IPTraverseObjListHierarchy(PObj, CrntViewMat,
			       Irit2WglProcessOneTraversedObject);

    /* Dump data to create WGL file. */
    Irit2WglDumpData(GlblScene, WGLFileName);

    /* Free scene after use. */
    Irit2WglFreeScene(GlblScene);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function of IPTraverseObjListHierarchy. Called on every non    *
* list object found in hierarchy and processes it.                           *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:   Objects to process.                                              *
*   Mat:    Transformation matrix to apply to this object.                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
void Irit2WglProcessOneTraversedObject(IPObjectStruct *PObj,
				       IrtHmgnMatType Mat)
{
    IPObjectStruct *PObjs;

    if (IP_IS_FFGEOM_OBJ(PObj))
        PObjs = IPConvertFreeForm(PObj, &IPFFCState);  /* Convert in place. */
    else
	PObjs = PObj;

    for (PObj = PObjs; PObj != NULL; PObj = PObj -> Pnext) {
        Irit2WglProcessObject(PObj);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Process the given object and add it to the Web GL frame.                 * 
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:   Object to process.                                               *
*   Mat:    Transformation matrix to apply for the given object.             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
void Irit2WglProcessObject(IPObjectStruct *PObj)
{
    int PolyRedInt, PolyGreenInt, PolyBlueInt, Color,
	**Polygons = NULL,
        PTextureFlip = FALSE,
        NewImage = FALSE;
    char PTextureFileName[IRIT_LINE_LEN_LONG],
	*PTextureFormat = NULL;
    const char
	*PTexture = NULL;
    float
	*UVs = NULL;
    IrtRType Transp, RGBA[4], PTextureSUV[3];
    IPVertexStruct
	**Vertices = NULL;
    IPPolyVrtxIdxStruct
	*PVIdx = NULL;
    Irit2WglModelObjectStruct
	*ModelObject = NULL;
    Irit2WglModelObjectTagType ModelObjectTag;
    IrtPtType Pos;
    IrtNrmlType Normal;
    IrtUVType UV,
	SUV = { 1.0, 1.0 };

    if (PObj -> ObjType == IP_OBJ_POINT) {
	if (!IP_ATTR_IS_BAD_INT(AttrGetObjectIntAttrib(PObj,
						       "LIGHT_SOURCE"))) {
	    const char *p;
	    int i, Red, Green, Blue;
	    IGLightType LightPos;
	    IrtVecType LightColor;

	    for (i = 0; i < 3; i++)
		LightPos[i] = (float) PObj -> U.Pt[i];
	    LightPos[3] = (float)
		!(((p = AttrGetObjectStrAttrib(PObj, "TYPE")) != NULL) &&
		  (stricmp(p, "POINT_INFTY") == 0));
		    
	    if (AttrGetObjectRGBColor(PObj, &Red, &Green, &Blue)) {
		LightColor[0] = Red / 255.0;
		LightColor[1] = Green / 255.0;
		LightColor[2] = Blue/ 255.0;
	    }
	    else {
		LightColor[0] = LightColor[1] = LightColor[2] = 1.0;
	    }

	    Irit2WglSetLightSource(GlblScene, LightPos, LightColor);
	}

	return;
    }

    if (PObj -> ObjType != IP_OBJ_POLY) {
	IRIT_WARNING_MSG_PRINTF("Non polygonal object \"%s\" detected and ignored.\n",
	    IP_GET_OBJ_NAME(PObj));
	return;
    }

    PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 7);

    if (IP_IS_POLYLINE_OBJ(PObj)) {
	ModelObjectTag = IRIT2WGL_MODEL_OBJECT_TAG_POLYLINE;
    }
    if (IP_IS_POLYGON_OBJ(PObj)) {
	ModelObjectTag = IRIT2WGL_MODEL_OBJECT_TAG_POLYGON;
    }
    else if (IP_IS_POLYSTRIP_OBJ(PObj)) {
	IRIT_WARNING_MSG_PRINTF("POLYSTRIP type detected and ignored.\n");
	return;
    }
    else {
	IRIT_WARNING_MSG_PRINTF("Unsupported poly type detected and ignored.\n");
	return;
    }

    /* Set object's color and opacity (white and visible are the defaults). */
    RGBA[0] = RGBA[1] = RGBA[2] = RGBA[3] = 1.0;

    if (AttrGetRGBColor(PObj -> Attr,
			&PolyRedInt, &PolyGreenInt, &PolyBlueInt)) {
	RGBA[0] = PolyRedInt / 255.0;
	RGBA[1] = PolyGreenInt / 255.0;
	RGBA[2] = PolyBlueInt / 255.0;
    }
    else {
	Color = AttrGetColor(PObj -> Attr);
	if (Color != IP_ATTR_NO_COLOR) {
	    AttrGetIndexColor(Color, &PolyRedInt, &PolyGreenInt, &PolyBlueInt);
	    RGBA[0] = PolyRedInt / 255.0;
	    RGBA[1] = PolyGreenInt / 255.0;
	    RGBA[2] = PolyBlueInt / 255.0;
	}
    }

    Transp = AttrGetRealAttrib(PObj -> Attr, "transp");
    if (Transp != IP_ATTR_BAD_REAL)
	RGBA[3] = 1 - Transp;

    /* Set object's parametric texture (if exists). */
    PTexture = AttrGetStrAttrib(PObj -> Attr, "ptexture");
    if (PTexture != NULL) {
	IrtImgParsePTextureString(PTexture, PTextureFileName, PTextureSUV,
				  &PTextureFlip, &NewImage);
	PTextureFormat = strrchr(PTextureFileName, '.');
	if ((strcmp(PTextureFormat, ".bmp") != 0) &&
	    (strcmp(PTextureFormat, ".gif") != 0) &&
	    (strcmp(PTextureFormat, ".png") != 0) &&
	    (strcmp(PTextureFormat, ".jpg") != 0)) {
	    IRIT_WARNING_MSG_PRINTF("Unsupported ptexture image format (\"%s\") detected and replaced with \".png\".\n",
				    PTextureFormat);
	    strncpy(PTextureFormat, ".png", 4);
	}
	SUV[0] = PTextureSUV[0];
	SUV[1] = PTextureSUV[1];
    }

    /* New Model object. */
    ModelObject = Irit2WglNewModelObject(GlblScene, PObj -> ObjName,
					 ModelObjectTag, RGBA,
					 PVIdx -> NumVrtcs,
					 (PTexture != NULL) ? PTextureFileName : NULL,
					 SUV);

    /* Add vertices to model object. */
    for (Vertices = PVIdx -> Vertices; *Vertices != NULL; Vertices++) {
	IRIT_PT_COPY(Pos, (*Vertices) -> Coord);
	IRIT_VEC_COPY(Normal, (*Vertices) -> Normal);
	IRIT_UV_RESET(UV);
	if ((UVs = AttrGetUVAttrib((*Vertices) -> Attr, "uvvals")) != NULL) {
	    UV[0] = UVs[0];
	    UV[1] = UVs[1];
	}
	Irit2WglAddVertex(GlblScene, ModelObject, Pos, Normal, UV);
    }

    /* Add poly to Model object. */
    for (Polygons = PVIdx -> Polygons; *Polygons != NULL; Polygons++)
	Irit2WglAddPoly(GlblScene, ModelObject, *Polygons);

    IPPolyVrtxIdxFree(PVIdx);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Create new scene.                                                          M
*                                                                            *
* PARAMETERS:                                                                M
*   None								     *
*                                                                            *
* RETURN VALUE:                                                              M
*   Irit2WglSceneStruct *: New scene.                                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglNewScene, NewScene                                               M
*****************************************************************************/
Irit2WglSceneStruct *Irit2WglNewScene(void)
{
    Irit2WglSceneStruct
	*NewScene = NULL;
    int i, j;

    NewScene =
	(Irit2WglSceneStruct *) IritMalloc(sizeof(Irit2WglSceneStruct));

    NewScene -> Objects = NewScene -> LastObject = NULL;

    NewScene -> Textures = NewScene -> LastTexture = NULL;

    for (i = 0; i < IRIT2WGL_MAX_LIGHT_SOURCE_NUM; i++) {
	/* By default all light sources are point light sources. */
	for (j = 0; j < 4; j++)
	    NewScene -> Lights[i].Pos[j] = GlblLightSourceInitialPos[i][j];

	for (j = 0; j < 3; j++) {
	    NewScene -> Lights[i].Diffuse[j] = 1.0;
	    NewScene -> Lights[i].Specular[j] = 1.0;
	}

	/* By default 2 light sources are active. */
	if (i < 2) {
	    for (j = 0; j < 3; j++) {
		NewScene -> Lights[i].DiffuseI = 0.5;
		NewScene -> Lights[i].SpecularI = 0.2;
	    }
	} else {
	    for (j = 0; j < 3; j++) {
		NewScene -> Lights[i].DiffuseI = 0.0;
		NewScene -> Lights[i].SpecularI = 0.0;
	    }
	}
    }
    NewScene -> LightSourceSetCount = 0;

    NewScene -> NOP = 0;

    return NewScene;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Free the given scene.                                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglFreeScene, Free                                                  M
*****************************************************************************/
void Irit2WglFreeScene(Irit2WglSceneStruct *Scene)
{
    Irit2WglTextureStruct
	*TextureIter = NULL,
	*TempTexture = NULL;
    Irit2WglModelObjectStruct
	*ObjIter = NULL,
	*TempObj = NULL;
    Irit2WglPolyStruct
	*PolyIter = NULL,
	*TempPoly = NULL;

    /* Free textures. */
    TextureIter = Scene -> Textures;
    while (TextureIter != NULL) {
	TempTexture = TextureIter;
	TextureIter = TextureIter -> Pnext;

	IritFree(TempTexture -> FileName);

	IritFree(TempTexture);
    }

    /* Free objects. */
    ObjIter = Scene -> Objects;
    while (ObjIter != NULL) {
	TempObj = ObjIter;
	ObjIter = ObjIter -> Pnext;

	IritFree(TempObj -> Name);

	PolyIter = TempObj -> Polys;
	while (PolyIter != NULL) {
	    TempPoly = PolyIter;
	    PolyIter = PolyIter -> Pnext;
	    IritFree(TempPoly -> Indices);
	    IritFree(TempPoly);
	}

	IritFree(TempObj -> Vertices);

	IritFree(TempObj);
    }

    IritFree(Scene);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Create new model object.                                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   Name: Model object's name.                                               M
*   ModelObjectTag: Model object tag.                                        M
*   RGBA: RGBA values.                                                       M
*   MaxObjVertices: Maximum number of vertices in object.                    M
*   PTextureFileName: Texture file name.                                     M
*   SUV: UV texture scaling.                                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   Irit2WglModelObjectStruct *: New model object.                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglNewModelObject, NewModelObject                                   M
*****************************************************************************/
Irit2WglModelObjectStruct *Irit2WglNewModelObject(Irit2WglSceneStruct *Scene,
						  const char *Name,
						  Irit2WglModelObjectTagType
						               ModelObjectTag,
						  IrtRType RGBA[4],
						  int MaxObjVertices,
						  const char *PTextureFileName,
						  IrtUVType SUV)
{
    Irit2WglModelObjectStruct
	*NewModelObject = NULL;
    int i;

    NewModelObject = (Irit2WglModelObjectStruct *)
	IritMalloc(sizeof(Irit2WglModelObjectStruct));

    /* Object's tag. */
    NewModelObject -> Tag = ModelObjectTag;

    /* Object's name. */
    NewModelObject -> Name = NULL;
    if (Name != NULL)
	NewModelObject -> Name = IritStrdup(Name);

    /* Object's RGBA. */
    for (i = 0; i < 4; i++)
	NewModelObject -> RGBA[i] = RGBA[i];

    /* Object's maximum number of vertices. */
    NewModelObject -> Vertices = (Irit2WglVertexStruct *)
	IritMalloc(MaxObjVertices * sizeof(Irit2WglVertexStruct));

    /* Object's parametric texture data. */
    NewModelObject -> Texture = NULL;
    if (PTextureFileName != NULL) {
	NewModelObject -> Texture = Scene -> Textures;
	while (NewModelObject -> Texture != NULL) {
	    if (strcmp(PTextureFileName,
		       NewModelObject -> Texture -> FileName) == 0)
		break;
	    NewModelObject -> Texture = NewModelObject -> Texture -> Pnext;
	}
	if (NewModelObject -> Texture == NULL) {
	    NewModelObject -> Texture = (Irit2WglTextureStruct *)
		IritMalloc(sizeof(Irit2WglTextureStruct));
	    NewModelObject -> Texture -> Pnext = NULL;
	    NewModelObject -> Texture -> FileName =
		(char *) IritMalloc(strlen(PTextureFileName) + 1);
	    strcpy(NewModelObject -> Texture -> FileName, PTextureFileName);

	    if (Scene -> Textures == NULL) {
		NewModelObject -> Texture -> Index = 0;
		Scene -> LastTexture = Scene -> Textures
		    = NewModelObject -> Texture;
	    }
	    else {
		NewModelObject -> Texture -> Index
		    = Scene -> LastTexture -> Index + 1;
		Scene -> LastTexture -> Pnext = NewModelObject -> Texture;
		Scene -> LastTexture = NewModelObject -> Texture;
	    }
	}

	IRIT_UV_COPY(NewModelObject -> SUV, SUV);

	IRIT_UV_RESET(NewModelObject -> MinUV);
	IRIT_UV_RESET(NewModelObject -> MaxUV);
    }

    NewModelObject -> VertexCount = 0;
    NewModelObject -> Polys = NewModelObject -> LastPoly = NULL;

    NewModelObject -> Pnext = NULL;
    if (Scene -> Objects == NULL) {
	Scene -> Objects = Scene -> LastObject = NewModelObject;
    }
    else {
	Scene -> LastObject -> Pnext = NewModelObject;
	Scene -> LastObject = NewModelObject;
    }

    return NewModelObject;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Add new vertex to the given model object.                                  M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   ModelObject: Model object.                                               M
*   Pos: Vertex position.                                                    M
*   Normal: Vertex normal.                                                   M
*   UV: Vertex UV coordinates.                                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglAddVertex, Vertex, Pos, Normal, UV                               M
*****************************************************************************/
void Irit2WglAddVertex(Irit2WglSceneStruct *Scene,
		       Irit2WglModelObjectStruct *ModelObject,
		       IrtPtType Pos,
		       IrtNrmlType Normal,
		       IrtUVType UV)
{
    int i;

    IRIT_PT_COPY(ModelObject -> Vertices[ModelObject -> VertexCount].Pos, Pos);

    IRIT_VEC_COPY(ModelObject -> Vertices[ModelObject -> VertexCount].Normal,
		  Normal);

    IRIT_UV_COPY(ModelObject -> Vertices[ModelObject -> VertexCount].UV, UV);

    if (ModelObject -> VertexCount == 0) {
	IRIT_UV_COPY(ModelObject -> MinUV, UV);
	IRIT_UV_COPY(ModelObject -> MaxUV, UV);
    }

    for (i = 0; i < 2; i++) {
	if (ModelObject -> VertexCount != 0) {
	    if (UV[i] < ModelObject -> MinUV[i])
		ModelObject -> MinUV[i] = UV[i];
	    if (UV[i] > ModelObject -> MaxUV[i])
		ModelObject -> MaxUV[i] = UV[i];
	}
    }
    
    ModelObject -> VertexCount++;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Create new poly to the given model object.                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   ModelObject: Model object.                                               M
*   PolyIndexArray: Array of vertex indices defining the poly.               M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglAddPoly, Poly, PolyIndexArray, AddPoly                           M
*****************************************************************************/
void Irit2WglAddPoly(Irit2WglSceneStruct *Scene,
		     Irit2WglModelObjectStruct *ModelObject,
		     int *PolyIndexArray)
{
    Irit2WglPolyStruct
	*NewPoly = NULL;
    int i, *pi, IndexCount;
	
    NewPoly = (Irit2WglPolyStruct *) IritMalloc(sizeof(Irit2WglPolyStruct));

    IndexCount = 0;
    pi = PolyIndexArray;
    while (*(pi++) != -1)
	IndexCount++;

#ifdef NASTY_PATCH_FOR_DRAW_ELEMENTS_LINE_LOOP_BUG
    if (ModelObject -> Tag == IRIT2WGL_MODEL_OBJECT_TAG_POLYGON)
	IndexCount++;
#endif	/* NASTY_PATCH_FOR_DRAW_ELEMENTS_LINE_LOOP_BUG */ 

    NewPoly -> Indices = (int *)IritMalloc(IndexCount * sizeof(int));
    for (i = 0; i < IndexCount; i++)
	NewPoly -> Indices[i] = PolyIndexArray[i];

#ifdef NASTY_PATCH_FOR_DRAW_ELEMENTS_LINE_LOOP_BUG
    if (ModelObject -> Tag == IRIT2WGL_MODEL_OBJECT_TAG_POLYGON)
	NewPoly -> Indices[IndexCount - 1] = PolyIndexArray[0];
#endif	/* NASTY_PATCH_FOR_DRAW_ELEMENTS_LINE_LOOP_BUG */

    NewPoly -> IndexCount = IndexCount;

    NewPoly -> Pnext = NULL;
    if (ModelObject -> Polys == NULL) {
	ModelObject -> Polys = ModelObject -> LastPoly = NewPoly;
    }
    else {
	ModelObject -> LastPoly -> Pnext = NewPoly;
	ModelObject -> LastPoly = NewPoly;
    }

    Scene -> NOP++;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Set light source data for the given scene.                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   LightPos: Light source position.                                         M
*   LightColor: Light source diffuse color.                                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglSetLightSource, LightSource, LightPos, LightColor                M
*****************************************************************************/
void Irit2WglSetLightSource(Irit2WglSceneStruct *Scene,
			    IGLightType LightPos,
			    IrtVecType LightColor)
{
    int i;

    if (Scene -> LightSourceSetCount == IRIT2WGL_MAX_LIGHT_SOURCE_NUM) {
	IRIT_WARNING_MSG_PRINTF("%d light sources already set. Light source ignored.\n",
	    IRIT2WGL_MAX_LIGHT_SOURCE_NUM);
	return;
    }

    for (i = 0; i < 4; i++)
	Scene -> Lights[Scene -> LightSourceSetCount].Pos[i] = LightPos[i];

    for (i = 0; i < 3; i++) {
	Scene -> Lights[Scene -> LightSourceSetCount].Diffuse[i] = LightColor[i];
	Scene -> Lights[Scene -> LightSourceSetCount].Specular[i] = 1.0;
    }

    Scene -> Lights[Scene -> LightSourceSetCount].DiffuseI
	= IRIT_MAX(IRIT_MAX(LightColor[0], LightColor[1]), LightColor[2]);
    Scene -> Lights[Scene -> LightSourceSetCount].SpecularI = 0.2;

    (Scene -> LightSourceSetCount)++;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump WebGL HTML data.                                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   OutputFileName: Output file name.                                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpData, DumpData                                               M
*****************************************************************************/
void Irit2WglDumpData(Irit2WglSceneStruct *Scene, const char* OutputFileName)
{
    char Line[IRIT2WGL_MAX_BUFFER_LEN];
    FILE
	*HtmlOutputFile = stdout,
#	if defined(__WINNT__) || defined(OS2GCC)
	    *HtmlSkeletonFile = fopen(Irit2WglFBuffer("%s\\%s",
						      (GlblIritBinDir == NULL)
						      ? "." : GlblIritBinDir,
						      IRIT2WGL_HTML_SKELETON_FILE_NAME),
				      "r");
#		    else
	    *HtmlSkeletonFile = fopen(Irit2WglFBuffer("%s/%s",
						      (GlblIritBinDir == NULL)
						      ? "." : GlblIritBinDir,
						      IRIT2WGL_HTML_SKELETON_FILE_NAME),
				      "r");
#		    endif /* __WINNT__ || OS2GCC */
    if (HtmlSkeletonFile == NULL)
	Irit2WglFatalError("Couldn't open HTML skeleton file");

    if (OutputFileName != NULL)
	HtmlOutputFile = fopen(Irit2WglFBuffer("%s.html", OutputFileName), "w");

    while (fgets(Line, IRIT2WGL_MAX_BUFFER_LEN, HtmlSkeletonFile) != NULL) {
	if (strncmp(Line, "<!-- Dump: canvas -->", 15 + 6) == 0) {
	    fputs("        <canvas id=\"WebGL-canvas\" oncontextmenu=\"return false;\"",
		  HtmlOutputFile);
	    fputs(Irit2WglFBuffer("\n        width=\"%d\" height=\"%d\">\n",
		    GlblCanvasWidth, GlblCanvasHeight),
		  HtmlOutputFile);
	    fputs("        </canvas>\n", HtmlOutputFile);
	}
	else if (strncmp(Line, "<!-- Dump: scripts -->", 15 + 7) == 0) {
	    Irit2WglDumpScripts(GlblScene, HtmlOutputFile, OutputFileName);
	} else {
	    fputs(Line, HtmlOutputFile);
	}
    }
    fclose(HtmlOutputFile);

    fclose(HtmlSkeletonFile);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump WebGL HTML scripts related data.                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   HtmlOutputFile: HTML Output file.                                        M
*   OutputFileName: Output file name.                                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpScripts, Scripts, DumpSctipts                                M
*****************************************************************************/
void Irit2WglDumpScripts(Irit2WglSceneStruct *Scene,
			 FILE* HtmlOutputFile,
			 const char* OutputFileName)
{
    /* Dump JS. */
    fputs("\n", HtmlOutputFile);
    if (OutputFileName != NULL) {
	fputs(Irit2WglFBuffer(
	      "  <script type=\"text/javascript\" src=\"%s.js\"></script>\n",
	      OutputFileName),
	    HtmlOutputFile);
	Irit2WglDumpJS(Scene, OutputFileName);
    }
    else {
	fputs("  <script type=\"text/javascript\">\n", HtmlOutputFile);
	Irit2WglDumpJS(Scene, OutputFileName);
	fputs("  </script>\n", HtmlOutputFile);
    }
    fputs("  <noscript>Your browser doesn't support JavaScript.</noscript>\n",
	HtmlOutputFile);

    /* Dump CSS. */
    fputs("\n", HtmlOutputFile);
    if (OutputFileName != NULL) {
	fputs(Irit2WglFBuffer("  <link rel=\"stylesheet\" type=\"text/css\" href=\"%s.css\">\n",
	      OutputFileName),
	      HtmlOutputFile);
	Irit2WglDumpCSS(Scene, OutputFileName);	
    }
    else {
	fputs("  <style type=\"text/css\">\n", HtmlOutputFile);
	Irit2WglDumpCSS(Scene, OutputFileName);
	fputs("  </style>\n", HtmlOutputFile);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data.                                               M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   OutputFileName: Output file name.                                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJS, JS                                                       M
*****************************************************************************/
void Irit2WglDumpJS(Irit2WglSceneStruct *Scene, const char *OutputFileName)
{
    char Line[IRIT2WGL_MAX_BUFFER_LEN];
    FILE
	*JSOutputFile = stdout,
#	if defined(__WINNT__) || defined(OS2GCC)
	    *JSSkeletonFile = fopen(Irit2WglFBuffer("%s\\%s",
						    (GlblIritBinDir == NULL)
						    ? "." : GlblIritBinDir,
						    IRIT2WGL_JS_SKELETON_FILE_NAME),
				    "r");
#		    else
	    *JSSkeletonFile = fopen(Irit2WglFBuffer("%s/%s",
						    (GlblIritBinDir == NULL)
						    ? "." : GlblIritBinDir,
						    IRIT2WGL_JS_SKELETON_FILE_NAME),
				    "r");
#		    endif /* __WINNT__ || OS2GCC */
    if (JSSkeletonFile == NULL)
	Irit2WglFatalError("Couldn't open JS skeleton file");

    if (OutputFileName != NULL)
	JSOutputFile = fopen(Irit2WglFBuffer("%s.js", OutputFileName), "w");

    while (fgets(Line, IRIT2WGL_MAX_BUFFER_LEN, JSSkeletonFile) != NULL) {
	if (strncmp(Line, "// Dump: bkcol", 9 + 5) == 0) {
	    fputs("// Background color constant:\n", JSOutputFile);
	    fprintf(JSOutputFile, "var BK_RGB = [%d, %d, %d];\n",
		    GlblBkRGB[0], GlblBkRGB[1], GlblBkRGB[2]);
	}
	else if (strncmp(Line, "// Dump: texture", 9 + 7) == 0) {
	    fputs("// Set texture data.\n", JSOutputFile);
	    Irit2WglDumpJSSetTextureData(Scene, JSOutputFile);
	}
	else if (strncmp(Line, "// Dump: model", 9 + 5) == 0) {
	    fputs("// Set model data.\n", JSOutputFile);
	    Irit2WglDumpJSSetModelData(Scene, JSOutputFile);
	}
	else if (strncmp(Line, "// Dump: camera", 9 + 6) == 0) {
	    fputs("// Initialize camera matrices.\n", JSOutputFile);
	    Irit2WglDumpJSInitCameraMatrices(Scene, JSOutputFile);
	}
	else if (strncmp(Line, "// Dump: light", 9 + 5) == 0) {
	    fputs("// Set light sources.\n", JSOutputFile);
	    Irit2WglDumpJSSetLightSources(Scene, JSOutputFile);
	}
	else if (strncmp(Line, "// Dump: ctrl", 9 + 4) == 0) {
	    fputs("// Set control bar parameters.\n", JSOutputFile);
	    Irit2WglDumpJSSetControlBarParams(Scene, JSOutputFile);
	}
	else if (strncmp(Line, "// Dump: menu", 9 + 4) == 0) {
	    fputs("// Set menu parameters.\n", JSOutputFile);
	    Irit2WglDumpJSSetMenuParams(Scene, JSOutputFile);
	}
	else if (strncmp(Line, "// Dump: log", 9 + 3) == 0) {
	    fputs("// Set status log parameters.\n", JSOutputFile);
	    Irit2WglDumpJSSetStatusLogParams(Scene, JSOutputFile);
	} else {
	    fputs(Line, JSOutputFile);
	}
    }

    if (OutputFileName != NULL)
	fclose(JSOutputFile);

    fclose(JSSkeletonFile);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Texture data.                                M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSSetTextureData, TextureData                                M
*****************************************************************************/
void Irit2WglDumpJSSetTextureData(Irit2WglSceneStruct *Scene,
				  FILE *JSOutputFile)
{
    Irit2WglTextureStruct
	*TextureIter = NULL;

    fprintf(JSOutputFile, "Textures.prototype.setData = function() {\n");

    for (TextureIter = Scene -> Textures;
	TextureIter != NULL;
	TextureIter = TextureIter -> Pnext) {

	/* Add the next texture. */
	fprintf(JSOutputFile, "\tthis.add(\'%s\');\n",
	    TextureIter -> FileName);
    }

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Model data.                                  M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSSetModelData, ModelData                                    M
*****************************************************************************/
void Irit2WglDumpJSSetModelData(Irit2WglSceneStruct *Scene,
				FILE *JSOutputFile)
{
    int i, j;
    Irit2WglModelObjectStruct
	*ObjIter = NULL;
    Irit2WglPolyStruct
	*PolyIter = NULL;

    fprintf(JSOutputFile, "Model.prototype.setData = function() {\n");

    for (ObjIter = Scene -> Objects, i = 0;
	ObjIter != NULL;
	ObjIter = ObjIter -> Pnext, i++) {

	/* New object. */
	fprintf(JSOutputFile,
		"\n\tthis.objs[%d] = new ModelObject(%d);\n", i, i);

	/* Set name for object. */
	fprintf(JSOutputFile, "\tthis.objs[%d].name = \'%s\';\n",
		i, ObjIter -> Name);

	/* Set tag for object. */
	fprintf(JSOutputFile, "\tthis.objs[%d].tag = ModelObjectTag.%s;\n",
		i, Irit2WglModelObjectTagType2Str(ObjIter -> Tag));

	/* Set RBGA for object. */
	fprintf(JSOutputFile, "\tthis.objs[%d].RGBA = [%lf, %lf, %lf, %lf];\n",
		i,
		ObjIter -> RGBA[0],
		ObjIter -> RGBA[1],
		ObjIter -> RGBA[2],
		ObjIter -> RGBA[3]);

	/* Dump vertices. */
	fprintf(JSOutputFile, "\n\tthis.objs[%d].vertices = [\n", i);
	for (j = 0; j < ObjIter -> VertexCount; j++) {
	    fprintf(JSOutputFile, "\t\t%lf, %lf, %lf,%s\n",
		    ObjIter -> Vertices[j].Pos[0],
		    ObjIter -> Vertices[j].Pos[1],
		    ObjIter -> Vertices[j].Pos[2],
		    (j == ObjIter -> VertexCount - 1) ? "\n\t];" : "");
	}

	/* Dump normals. */
	fprintf(JSOutputFile, "\n\tthis.objs[%d].normals = [\n", i);
	for (j = 0; j < ObjIter -> VertexCount; j++) {
	    fprintf(JSOutputFile, "\t\t%lf, %lf, %lf,%s\n",
		    ObjIter -> Vertices[j].Normal[0],
		    ObjIter -> Vertices[j].Normal[1],
		    ObjIter -> Vertices[j].Normal[2],
		    (j == ObjIter -> VertexCount - 1) ? "\n\t];" : "");
	}

	/* Dump polys. */
	fprintf(JSOutputFile, "\n\tthis.objs[%d].polys = [\n", i);
	for (PolyIter = ObjIter -> Polys;
	    PolyIter != NULL;
	    PolyIter = PolyIter -> Pnext) {
	    fprintf(JSOutputFile, "\t\t[");
	    for (j = 0; j < PolyIter -> IndexCount; j++) {
		fprintf(JSOutputFile, "%d%s",
			PolyIter -> Indices[j],
			(j == PolyIter -> IndexCount - 1) ? "]" : ", "); 
	    }
	    fprintf(JSOutputFile, ",%s\n",
		    (PolyIter -> Pnext == NULL) ? "\n\t];" : "");
	}

	/* Dump texture data and UVs. */
	if (ObjIter -> Texture != NULL) {
	    IrtUVType
	        SUV = {
		    ObjIter -> SUV[0] / (ObjIter -> MaxUV[0] -
					 ObjIter -> MinUV[0]),
		    ObjIter -> SUV[1] / (ObjIter -> MaxUV[1] -
					 ObjIter -> MinUV[1])
	        };

	    fprintf(JSOutputFile,
		    "\n\tthis.objs[%d].texture.use = true;\n\tthis.objs[%d].texture.index = %d;\n",
		    i, i, ObjIter -> Texture -> Index);

	    fprintf(JSOutputFile, "\n\tthis.objs[%d].texture.UVs = [\n", i);
	    for (j = 0; j < ObjIter -> VertexCount; j++) {
		fprintf(JSOutputFile, "\t\t%lf, %lf,%s\n",
			ObjIter -> Vertices[j].UV[0] * SUV[0],
			ObjIter -> Vertices[j].UV[1] * SUV[1],
			(j == ObjIter -> VertexCount - 1) ? "\n\t];" : "");
	    }
	}
    }

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Camera matrices.                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSInitCameraMatrices, CameraMatrices                         M
*****************************************************************************/
void Irit2WglDumpJSInitCameraMatrices(Irit2WglSceneStruct *Scene,
				      FILE *JSOutputFile)
{
    fprintf(JSOutputFile, "Camera.prototype.initMatrices = function() {\n");

    /* Set view matrix. */
    fprintf(JSOutputFile,
	    "\tthis.viewMatrices[0] = [\n\t\t[%lf, %lf, %lf, %lf],\n\t\t[%lf, %lf, %lf, %lf],\n\t\t[%lf, %lf, %lf, %lf],\n\t\t[%lf, %lf, %lf, %lf]\n\t];\n",
	    IPViewMat[0][0], IPViewMat[0][1], IPViewMat[0][2], IPViewMat[0][3],
	    IPViewMat[1][0], IPViewMat[1][1], IPViewMat[1][2], IPViewMat[1][3],
	    IPViewMat[2][0], IPViewMat[2][1], IPViewMat[2][2], IPViewMat[2][3],
	    IPViewMat[2][0], IPViewMat[3][1], IPViewMat[3][2], IPViewMat[3][3]
    );

    /* Set perspective matrix. */
    fprintf(JSOutputFile,
	    "\tthis.prspMatrix = [\n\t\t[%lf, %lf, %lf, %lf],\n\t\t[%lf, %lf, %lf, %lf],\n\t\t[%lf, %lf, %lf, %lf],\n\t\t[%lf, %lf, %lf, %lf]\n\t];\n",
	    IPPrspMat[0][0], IPPrspMat[0][1], IPPrspMat[0][2], IPPrspMat[0][3],
	    IPPrspMat[1][0], IPPrspMat[1][1], IPPrspMat[1][2], IPPrspMat[1][3],
	    IPPrspMat[2][0], IPPrspMat[2][1], IPPrspMat[2][2], IPPrspMat[2][3],
	    IPPrspMat[2][0], IPPrspMat[3][1], IPPrspMat[3][2], IPPrspMat[3][3]
    );

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Light sources setting.                       M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSSetLightSources, LightSources                              M
*****************************************************************************/
void Irit2WglDumpJSSetLightSources(Irit2WglSceneStruct *Scene,
				   FILE *JSOutputFile)
{
    int i;

    fprintf(JSOutputFile, "Scene.prototype.setLightSources = function() {\n");

    /* Set light sources. */
    for (i = 0; i < IRIT2WGL_MAX_LIGHT_SOURCE_NUM; i++) {
	/* New light source. */
	fprintf(JSOutputFile, "\tthis.lights[%d] = new Light(%d);\n", i + 1,
		i + 1);

	/* Set light source type. */
	if (Scene -> Lights[i].Pos[3] > 0.0) {
	    fprintf(JSOutputFile,
		    "\tthis.lights[%d].type = LightType.POINT;\n",
		    i + 1);
	}
	else {
	    fprintf(JSOutputFile,
		    "\tthis.lights[%d].type = LightType.DIRECTIONAL;\n",
		    i + 1);
	}

	/* Set light source position. */
	fprintf(JSOutputFile, "\tthis.lights[%d].setPos([%lf, %lf, %lf]);\n",
		i + 1, Scene -> Lights[i].Pos[0], Scene -> Lights[i].Pos[1],
		Scene -> Lights[i].Pos[2]);

	/* Set light source direction. */
	fprintf(JSOutputFile, "\tthis.lights[%d].setDir([%lf, %lf, %lf]);\n",
		i + 1, Scene -> Lights[i].Pos[0], Scene -> Lights[i].Pos[1],
		Scene -> Lights[i].Pos[2]);

	/* Set light source attenuation. */
	fprintf(JSOutputFile,
		"\tthis.lights[%d].attK = LIGHT_ATTENUATION;\n",
		i + 1);

	/* Set light source diffuse. */
	fprintf(JSOutputFile, "\tthis.lights[%d].diffuse = [%lf, %lf, %lf];\n",
		i + 1, Scene -> Lights[i].Diffuse[0],
		Scene -> Lights[i].Diffuse[1], Scene -> Lights[i].Diffuse[2]);

	fprintf(JSOutputFile, "\tthis.lights[%d].diffuseI = %lf;\n",
		i + 1, Scene -> Lights[i].DiffuseI);

	/* Set light source specular. */
	fprintf(JSOutputFile, "\tthis.lights[%d].specular = [%lf, %lf, %lf];\n",
		i + 1, Scene -> Lights[i].Specular[0],
		Scene -> Lights[i].Specular[1], Scene -> Lights[i].Specular[2]);

	fprintf(JSOutputFile, "\tthis.lights[%d].specularI = %lf;\n",
		i + 1, Scene -> Lights[i].SpecularI);

	fprintf(JSOutputFile, "\n");
    }

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Control bar parameters.                      M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSSetControlBarParams, ControlBar                            M
*****************************************************************************/
void Irit2WglDumpJSSetControlBarParams(Irit2WglSceneStruct *Scene,
				       FILE *JSOutputFile)
{
    fprintf(JSOutputFile, "ControlBar.prototype.setParams = function() {\n");

    /* Show/Hide control bar. */
    fprintf(JSOutputFile, "\t// Show/Hide control bar.\n");

    fprintf(JSOutputFile, "\tthis.show = %s;\n",
	    (GlblHideCtrlBar ? "false" : "true"));

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Menu parameters.                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSSetMenuParams, Menu                                        M
*****************************************************************************/
void Irit2WglDumpJSSetMenuParams(Irit2WglSceneStruct *Scene,
				 FILE *JSOutputFile)
{
    fprintf(JSOutputFile, "Menu.prototype.setParams = function() {\n");

    /* Scene menu parameters. */
    fprintf(JSOutputFile, "\t// Scene menu parameters.\n");

    /* Show world axes? */
    fprintf(JSOutputFile, "\tthis.scene.showAxes = %s;\n",
	    (GlblShowWorldAxes ? "true" : "false"));

    /* Enable depth test? */
    fprintf(JSOutputFile, "\tthis.scene.enableDepthTest = %s;\n",
	    (GlblDisableDepthTest ? "false" : "true"));

    /* Enable picking? */
    fprintf(JSOutputFile, "\tthis.scene.enablePicking = %s;\n",
	    (GlblEnablePicking ? "true" : "false"));

    /* Model menu parameters. */
    fprintf(JSOutputFile, "\n\t// Model menu parameters.\n");

    /* Show model axes? */
    fprintf(JSOutputFile, "\tthis.model.showAxes = %s;\n",
	    (GlblShowModelAxes  ? "true" : "false"));

    /* Draw mode. */
    fprintf(JSOutputFile, "\tthis.model.drawMode = %#06x;\n", GlblDrawMode);

    /* World transformation? */
    fprintf(JSOutputFile, "\tthis.model.worldTrans = %s;\n",
	    (GlblModelTrans  ? "false" : "true"));

    /* Camera menu parameters. */
    fprintf(JSOutputFile, "\n\t// Camera menu parameters.\n");

    /* Camera view angle. */
    fprintf(JSOutputFile, "\tthis.camera.viewAngle = ViewAngle.%s;\n",
	    Irit2WglViewAngleType2Str(GlblViewAngle));

    /* Camera projection mode. */
    fprintf(JSOutputFile, "\tthis.camera.projMode = ProjectionMode.%s;\n",
	    Irit2WglProjectionModeType2Str(GlblProjectionMode));

    /* Light menu parameters. */
    fprintf(JSOutputFile, "\n\t// Light menu parameters.\n");

    /* Global ambient light. */
    fprintf(JSOutputFile, "\tthis.light.ambient = [%lf, %lf, %lf];\n",
	    GlblAmbientRGB[0], GlblAmbientRGB[1], GlblAmbientRGB[2]);

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump JS script related data - Status log parameters.                       M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   JSOutputFile: JS Output file.                                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpJSSetStatusLogParams, StatusLog                              M
*****************************************************************************/
void Irit2WglDumpJSSetStatusLogParams(Irit2WglSceneStruct *Scene,
				      FILE *JSOutputFile)
{
    fprintf(JSOutputFile, "StatusLog.prototype.setParams = function() {\n");

    /* Number of polygons. */
    fprintf(JSOutputFile, "\t// Number of polygons.\n");
    fprintf(JSOutputFile, "\tthis.NOP = %d;\n", Scene -> NOP);

    /* Resolution. */
    fprintf(JSOutputFile, "\n\t// Resolution.\n");
    fprintf(JSOutputFile, "\tthis.resolution = %lf;\n", IPFFCState.FineNess);

    fprintf(JSOutputFile, "};\n");
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Dump CSS script related data.                                              M
*                                                                            *
* PARAMETERS:                                                                M
*   Scene: Scene.                                                            M
*   OutputFileName: Output file name.                                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglDumpCSS, CSS                                                     M
*****************************************************************************/
void Irit2WglDumpCSS(Irit2WglSceneStruct *Scene,
		     const char *OutputFileName)
{
    char Line[IRIT2WGL_MAX_BUFFER_LEN];
    FILE
	*CSSOutputFile = stdout,
#	if defined(__WINNT__) || defined(OS2GCC)
	    *CSSSkeletonFile = fopen(Irit2WglFBuffer("%s\\%s",
						     (GlblIritBinDir == NULL)
						     ? "." : GlblIritBinDir,
						     IRIT2WGL_CSS_SKELETON_FILE_NAME),
				     "r");
#		    else
	    *CSSSkeletonFile = fopen(Irit2WglFBuffer("%s/%s",
						     (GlblIritBinDir == NULL)
						     ? "." : GlblIritBinDir,
						     IRIT2WGL_CSS_SKELETON_FILE_NAME),
				     "r");
#		    endif /* __WINNT__ || OS2GCC */
    if (CSSSkeletonFile == NULL)
	Irit2WglFatalError("Couldn't open CSS skeleton file");

    if (OutputFileName != NULL)
	CSSOutputFile = fopen(Irit2WglFBuffer("%s.css", OutputFileName), "w");

    while (fgets(Line, IRIT2WGL_MAX_BUFFER_LEN, CSSSkeletonFile) != NULL)
	fputs(Line, CSSOutputFile);

    if (OutputFileName != NULL)
	fclose(CSSOutputFile);

    fclose(CSSSkeletonFile);
}


/*****************************************************************************
* DESCRIPTION:                                                               M
* Convert view angle enum to string.                                         M
*                                                                            *
* PARAMETERS:                                                                M
*   ViewAngle: View angle.                                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *: View angle string.                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglViewAngleType2Str, ViewAngle                                     M
*****************************************************************************/
const char *Irit2WglViewAngleType2Str(Irit2WglViewAngleType ViewAngle)
{
    switch (ViewAngle) {
	case IRIT2WGL_VIEW_ANGLE_ORIGINAL:
	    return "ORIGINAL";
	case IRIT2WGL_VIEW_ANGLE_FRONT:
	    return "FRONT";
	case IRIT2WGL_VIEW_ANGLE_BACK:
	    return "BACK";
	case IRIT2WGL_VIEW_ANGLE_RIGHT:
	    return "RIGHT";
	case IRIT2WGL_VIEW_ANGLE_LEFT:
	    return "LEFT";
	case IRIT2WGL_VIEW_ANGLE_TOP:
	    return "TOP";
	case IRIT2WGL_VIEW_ANGLE_BOTTOM:
	    return "BOTTOM";
	default:
	    return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Convert projection mode enum to string.                                    M
*                                                                            *
* PARAMETERS:                                                                M
*   ProjectionMode: Projection mode.                                         M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *: Projection mode string.                                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglProjectionModeType2Str, ProjMode                                 M
*****************************************************************************/
const char *Irit2WglProjectionModeType2Str(Irit2WglProjectionModeType ProjectionMode)
{
    switch (ProjectionMode) {
	case IRIT2WGL_PROJECTION_MODE_ORTHOGONAL:
	    return "ORTHOGONAL";
	case IRIT2WGL_PROJECTION_MODE_PERSPECTIVE:
	    return "PERSPECTIVE";
	default:
	    return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Convert model object tag enum to string.                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   ModelObjectTag: Model object tag.                                        M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *: Object tag string.                                         M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglModelObjectTagType2Str, ObjTag                                   M
*****************************************************************************/
const char *Irit2WglModelObjectTagType2Str(Irit2WglModelObjectTagType ModelObjectTag)
{
    switch (ModelObjectTag) {
	case IRIT2WGL_MODEL_OBJECT_TAG_POLYLINE:
	    return "POLYLINE";
	case IRIT2WGL_MODEL_OBJECT_TAG_POLYGON:
	    return "POLYGON";
	default:
	    return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Trap Irit2Wgl errors right here. Gets an error description, print it and   M
* exit the program using exit.                                               M
*                                                                            *
* PARAMETERS:                                                                M
*   FatalErrorMsg:      Fatal error message.                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   Irit2WglFatalError                                                       M
*****************************************************************************/
void Irit2WglFatalError(const char *FatalErrorMsg)
{
    IRIT_WARNING_MSG_PRINTF("IRIT2WGL: %s\n", FatalErrorMsg);
    exit(-1);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Return formatted string according to the given argument list.              *
*                                                                            *
* PARAMETERS:                                                                *
*   va_alist:    Do "man stdarg".                                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
#ifdef USE_VARARGS
const char *Irit2WglFBuffer(char *va_alist, ...)
{
    static char Buffer[IRIT2WGL_MAX_BUFFER_LEN]
    char *Format;
    va_list ArgPtr;

    va_start(ArgPtr);
    Format = va_arg(ArgPtr, char *);
#else
const char *Irit2WglFBuffer(char *Format, ...)
{
    static char Buffer[IRIT2WGL_MAX_BUFFER_LEN];
    va_list ArgPtr;

    va_start(ArgPtr, Format);
#endif /* USE_VARARGS */

    vsprintf(Buffer, Format, ArgPtr);
    va_end(ArgPtr);

    return Buffer;
}
