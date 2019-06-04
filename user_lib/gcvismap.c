/*****************************************************************************
* GeomCovr.c - Calculating Geometric Covering.                               *
******************************************************************************
* Written by Nadav Shragai, February 2010.                                   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/rndr_lib.h"
#include "user_loc.h"

#include <setjmp.h>
#include <time.h>
#include <math.h>

typedef enum UserGCAxesType {
    USER_GC_X_AXIS = 0,
    USER_GC_Y_AXIS,
    USER_GC_Z_AXIS,
} UserGCAxesType;

IRIT_GLOBAL_DATA const int 
    USER_GC_TIME_ZONE = 2;
    
/* Parameters for the perspective matrix. */
IRIT_GLOBAL_DATA const IrtRType 
    USER_GC_NEAR = 0.01, 
    USER_GC_FAR = 100;
    
IRIT_GLOBAL_DATA const IrtPtType 
    USER_GC_INF_VEC = { IRIT_INFNTY, IRIT_INFNTY, IRIT_INFNTY };

/* For exiting inner function and returning without a solution. */
IRIT_GLOBAL_DATA jmp_buf UserGCStartOfProcess;

IRIT_STATIC_DATA const char 
    *USER_GC_UV_VALUES = "uvvals",
    *USER_GC_CLOSING_CIRCLE = "USER_GC_CLOSING_CIRCLE",
    *USER_GC_VIEW_MATRIX_NAME = "VIEW_MAT",
    *USER_GC_TERRAIN_SRF_NAME = "Terrain",
    *USER_GC_PROCESSED_GEO_OBJ_NAME = "ProcessedGeoObj",
    *USER_GC_PROCESSED_OBSTACLES_NAME = "ProcessedObstacles",
    *USER_GC_PERIMETER_CRV_FILE_FAIL = "Can't load the object file \"%s\". \n%s\n",
    *USER_GC_PERIMETER_CRV_EMPTY_FILE = "The file doesn't exist or doesn't contain any data.",
    *USER_GC_PERIMETER_CRV_BAD_CONTENT = "The file doesn't contain exactly the expected elements.",
    *USER_GC_PERIMETER_CRV_BAD_DIM = "The curve isn't a 2D curve.",
    *USER_GC_PERIMETER_CRV_NO_IN_DOMAIN = "The curve ins't in the surface parametric domain.",
    *USER_GC_NO_SOLVING_TITLE = "Can't start solving the problem",
    *USER_GC_NO_OPS = "No observations points could be created.",
    *USER_GC_RESOLUTION = "resolution",
    *USER_GC_U_RESOLUTION = "u_resolution",
    *USER_GC_V_RESOLUTION = "v_resolution",
    *USER_GC_ALREADY_TRIANGULATED_ATTRIB = "AlreadyTriangulated",
    *USER_GC_VISMAP_FILE_NAME_BASE = "VisMap";

IRIT_STATIC_DATA IrtVecType UserGCSixMainAxisDirection[] = {
    {  1,  0,  0 },
    { -1,  0,  0 },
    {  0,  1,  0 },
    {  0, -1,  0 },
    {  0,  0,  1 },
    {  0,  0, -1 }
};

IRIT_STATIC_DATA int UserGCVisMapLine, UserGCVisMapSize[2],
    UserGCSrf2PlgOptimal = TRUE;

IRIT_STATIC_DATA IrtRType 
    UserGCSrf2PlgFineness = 0.01;

IRIT_STATIC_DATA const IrtRType 
    USER_GC_ANGLE_OVARLAP = 5,
    /* If this value is changed, the number of maximum returned angles in    */
    /* UserGCDivideAngle might change and require change in parameters and   */
    /* description.                                                          */
    USER_GC_MAX_ANGLE = 120;

/* The general up direction for view matrix creation. */
IRIT_STATIC_DATA const IrtPtType 
    USER_GC_UP = { 0, 1, 0 };

/* Variables for the view map image creation (UserGCProblemDefinitionStruct  */
/* can't be given to the relevant functions).                                */
IRIT_STATIC_DATA IrtImgPixelStruct 
    *UserGCVisMap = NULL;

static void UserGCSaveObjectAndMatrix(char *FileName, 
				      IrtHmgnMatType Mat, 
				      IPObjectStruct *PObj);
static void UserGCCreatePrspMatrix(IrtRType ZAngle, 
                                   IrtRType XYAngle,
                                   IrtHmgnMatType PrspMat);
static void UserGCFindNewAxes(const UserGCObsPtSuggestionStruct *Op, 
                              const IrtVecType Up,
                              IrtPtType x,
                              IrtPtType y,
                              IrtPtType z);
static void UserGCCreateViewMatrix2(const UserGCObsPtSuggestionStruct *Op,
                                    IrtHmgnMatType ViewMat,
                                    const IrtVecType Up);
static void UserGCCreateViewMatrix(const UserGCObsPtSuggestionStruct *Op,
                                   IrtHmgnMatType ViewMat);
static void UserGCDivideAndCreateViewMatrices(UserGCObsPtSuggestionStruct *ObsPt,
                                              IrtRType OpeningInXY, 
                                              IrtRType OpeningInZ,
                                              int *ObsPtsNum,                              
                                              IrtRType *OpeningOutXY,
                                              IrtRType *OpeningOutZ,
                                              IrtHmgnMatType ViewMats[6]);
static int UserGCGetIntersection(IrtPtType Pt0,
				 IrtPtType Pt1,
				 UserGCAxesType Axis,
				 IrtRType Param,
				 IrtRType *Blend);
static IPPolygonStruct *UserGCCropPolygon(IPPolygonStruct *Polygon,
					  UserGCAxesType Axis,
					  int IsSmaller,
					  IrtRType Param); 
static IPPolygonStruct *UserGCCropTriangle(IPPolygonStruct *Triangle,
					   IrtRType AngleXY,
					   IrtRType AngleZ,
					   IrtRType Near,
					   IrtRType Far);
static void UserGCTransformAndCrop(IRndrPtrType Rend,
				   IPObjectStruct *Scene,
				   IrtHmgnMatType ViewMat,
				   IrtHmgnMatType PrspMat,
				   IrtRType AngleXY,
				   IrtRType AngleZ,
				   const char *PicSuffix,
				   int SaveObject);
static void UserGCScanObjects(IRndrPtrType Rend,
                              IPObjectStruct *Scene);
static IrtImgPixelStruct *UserGCCreateVisMap(
				       UserGCProblemDefinitionStruct *Problem,     
				       IPObjectStruct *Scene, 
                                       IrtHmgnMatType ViewMat,
				       IrtHmgnMatType PrspMat,
                                       IrtRType AngleXY,
                                       IrtRType AngleZ,
                                       const char *PicSuffix);
static void UserGCDivideAngle(IrtRType OpeningIn, 
                              int *AnglesNum, 
                              IrtRType Angles[3], 
                              IrtRType *OpeningOut);
static void UserGCDivideOP(UserGCObsPtSuggestionStruct *ObsPt,
                           IrtRType OpeningInXY, 
                           IrtRType OpeningInZ,
                           int *ObsPtsNum,
                           UserGCObsPtSuggestionStruct *ObsPts,
                           IrtRType *OpeningOutXY,
                           IrtRType *OpeningOutZ);
static IrtImgPixelStruct *UserGCCombineVisMaps(IrtImgPixelStruct *VisMap1, 
                                               IrtImgPixelStruct *VisMap2,
                                               int Width,
                                               int Height);
static IrtImgPixelStruct *UserGCCreateComplexVisMap(
				       UserGCProblemDefinitionStruct *Problem,     
				       IPObjectStruct *Scene, 
                                       UserGCObsPtGroupTypeStruct *OpGroup,
				       UserGCObsPtSuggestionStruct *Op,
				       int PicIndex);
static void UserGCReallocSuggestions(UserGCObsPtGroupTypeStruct *ObsPtGroup,
                                     int NewSuggestionsNumber);
static void UserGCSetSuggestionPointsFromArray(
				      UserGCObsPtGroupTypeStruct *ObsPtGroup,
				      IrtVecType *Array,
				      int Number);
static IPPolygonStruct *UserGCCreateTriangle(IPPolygonStruct *AddTo,
                                             IrtPtType PV1,
                                             IrtPtType PV2,
                                             IrtPtType PV3);
static void UserGCSetSuggestionPointsFromPolyhedron(
                                     UserGCProblemDefinitionStruct *Problem,
                                     UserGCObsPtGroupTypeStruct *ObsPtGroup);
static IPObjectStruct *UserGCFlatObj(IPObjectStruct *PObj, void *Param);
static int UserGCLoadSrfAndCrv(const char *FileName,
			       CagdCrvStruct **Crv,
			       CagdSrfStruct **Srf,
			       IPObjectStruct **Obstacles);
static int UserGCSetSuggestionPointsFromPerimeter(
                                      UserGCObsPtGroupTypeStruct *ObsPtGroup);
static IPPolygonStruct *UserGCSetColor(IPPolygonStruct *Pl, void *Param);
static void UserGCSetColorPlList(IPPolygonStruct *Pl, int r, int g, int b);
static IPObjectStruct *UserGCTriangulateObjAux(IPObjectStruct *PObj, 
                                               void *Param);
static IPObjectStruct *UserGCTriangulateObj(IPObjectStruct *PObj,
                                            IrtRType SrfPlgfine,
                                            int SrfPlgOptimal);
static int UserGCSetSuggestionPointsFromSurfaceGrid(
                                      UserGCProblemDefinitionStruct *Problem,
                                      UserGCObsPtGroupTypeStruct *ObsPtGroup);
static IPObjectStruct *UserGCPrepareObjAux(IPObjectStruct *PObj, 
                                           void *Param);
static IPObjectStruct *UserGCPrepareObj(IPObjectStruct *PObj, 
                                        int MapWidth, 
                                        int MapHeight,
                                        IPObjectStruct *PObj2);
static IrtImgImageType UserGCImgWriteSetType(const char *ImageType);
static int UserGcImgWriteOpenFile(const char **argv,
                                  const char *FName,
                                  int Alpha,
                                  int XSize,
                                  int YSize);
static void UserGCImgWritePutLine(IrtBType *Alpha, IrtImgPixelStruct *Pixels);
static void UserGCImgWriteCloseFile(void);
static FILE *UserGCOpenFile(const char *FileName, int ForWriting, char *ErrMsg);
static void UserGCSaveMatrix(int Index, IrtHmgnMatType Mat);
static void UserGCSaveVisMap(int Index, 
                             IrtImgPixelStruct *VisMap, 
                             int MapWidth, 
                             int MapHeightm,
                             UserGCSaveRgbImageToFileFuncType SaveImageFunc);
static void UserGCDeleteVisMap(int StartIndex, int EndIndex, char *Extension);
static int UserGCSaveObjectToFile(const char *FileName, 
                                  IPObjectStruct* PObj,
                                  int IsBinary,
                                  int AllList);
static int UserGCLoadObjectFromFile(const char *FileName,
                                    IPObjectStruct **PPObj,
                                    int IsBinary);
static IPObjectStruct *UserGCLoadProcessedObjects2(
                                     UserGCProblemDefinitionStruct *Problem);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Open the given file for reading or writing. Print warning using          *
* IRIT_WARNING_MSG_PRINTF if fails and then return NULL.                     *
*                                                                            *
* PARAMETERS:                                                                *
*   FileName:   The name of the file.                                        *
*   ForWriting: TRUE to open the file for writing. FALSE for reading.        *
*   ErrMsg:     The message to display if error occurs. It should be phrased *
*               so that this function could add the name of the file directly*
*               after it. e.g.: "Error opening file".                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   FILE*: The opened file or NULL if failed.                                *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCOpenFile                                                           *
*****************************************************************************/
static FILE *UserGCOpenFile(const char *FileName, int ForWriting, char *ErrMsg)
{
    FILE *File;

    if ((File = fopen(FileName, ForWriting ? "wb" : "rb")) == NULL) {
        IRIT_WARNING_MSG_PRINTF("%s \"%s\".", ErrMsg, FileName);
        return NULL;
    }
    else 
        return File;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Saving matrix with the given index to disk.                              *
*   If an error occur it longjmp to UserGCStartOfProcess.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Index:         The index of the matrix.                                  *
*   Mat:           The 
*   VisMap:        The visibility map to save.                               *
*   MapWidth:      The width of the image.                                   *
*   MapHeight:     The height of the image.                                  *
*   SaveImageFunc: Function for saving the image.                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCSaveMatrix(int Index, IrtHmgnMatType Mat)
{
    char FileName[IRIT_LINE_LEN];
    FILE *Out;
    IPObjectStruct *PObj;

    sprintf(FileName, "%s%d.imd", USER_GC_VISMAP_FILE_NAME_BASE, Index);
    Out = UserGCOpenFile(FileName, TRUE, "Can't open for writing file");
    if (!Out)
        longjmp(UserGCStartOfProcess, 1);
    PObj = IPGenMatObject("VIEW_MAT", Mat, NULL);
    IPPutObjectToFile(Out, PObj, FALSE);
    IPFreeObject(PObj);
    fclose(Out);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Saving visibilty map with the given index to disk.                       *
*   If an error occur it longjmp to UserGCStartOfProcess.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Index:         The index of the visibility map.                          *
*   VisMap:        The visibility map to save.                               *
*   MapWidth:      The width of the image.                                   *
*   MapHeight:     The height of the image.                                  *
*   SaveImageFunc: Function for saving the image.                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCSaveVisMap(int Index, 
                             IrtImgPixelStruct *VisMap, 
                             int MapWidth, 
                             int MapHeight,
                             UserGCSaveRgbImageToFileFuncType SaveImageFunc)
{
    char FileName[IRIT_LINE_LEN];

    sprintf(FileName, "%s%d", USER_GC_VISMAP_FILE_NAME_BASE, Index);

    if (!SaveImageFunc(FileName, VisMap, MapWidth, MapHeight))
        longjmp(UserGCStartOfProcess, 1);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Exposing UserGCLoadVisMap outside.                                       M
*   Loading visibilty map with the given index from disk where               M
* UserGCSaveVisMap saved it earlier.                                         M
*   If an error occur returns NULL.                                          M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:    The geometric covering problem to be solved.                 M
*   Index:      The index of the visibility map.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct *: The loaded visibility map (new allocated memory     M
*                        using IritMalloc) or NULL in case of an error.      M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserGCSaveVisMap, UserGCLoadVisMap2			                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCLoadVisMap					                     M
*****************************************************************************/
IrtImgPixelStruct *UserGCLoadVisMap(UserGCProblemDefinitionStruct *Problem, 
                                    int Index)
{
    if (setjmp(UserGCStartOfProcess)) {
        return NULL;
    }

    return UserGCLoadVisMap2(Problem, Index);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Loading visibilty map with the given index from disk where               M
* UserGCSaveVisMap saved it earlier.                                         M
*   If an error occur it longjmp to UserGCStartOfProcess.                    M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:    The geometric covering problem to be solved.                 M
*   Index:      The index of the visibility map.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct *: The loaded visibility map (new allocated memory     M
*                        using IritMalloc).                                  M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserGCSaveVisMap, UserGCDeleteVisMap		                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCLoadVisMap2					                     M
*****************************************************************************/
IrtImgPixelStruct *UserGCLoadVisMap2(UserGCProblemDefinitionStruct *Problem, 
                                     int Index)
{
    IrtImgPixelStruct *Res;
    int LocalWidth, LocalHeight;
    char FileName[IRIT_LINE_LEN];

    Problem -> DebugParams.Print("%s Loading image %d\n", 
              MiscISCGetTimeStamp(time(NULL), NULL, USER_GC_TIME_ZONE, FALSE), 
              Index);
    sprintf(FileName, "%s%d", USER_GC_VISMAP_FILE_NAME_BASE, Index);
    Res = Problem -> DebugParams.LoadImageFunc(FileName, &LocalWidth, 
                                               &LocalHeight);
    if (!Res)
        longjmp(UserGCStartOfProcess, 1);

    if ((LocalWidth != Problem -> SolvingParams.VisMapWidth) || 
        (LocalHeight != Problem -> SolvingParams.VisMapHeight)) {
        IRIT_WARNING_MSG_PRINTF(
            "Width(%d) or height(%d) of the image in the file \"%s\" are incompatible. \n", 
            LocalWidth, LocalHeight, FileName);
        IritFree(Res);
        longjmp(UserGCStartOfProcess, 1);
    }
    return Res;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Deleting all visibilty maps at the given index range from disk.          *
*   No error produced if any of the files doesn't exist.                     *
*                                                                            *
* PARAMETERS:                                                                *
*   StartIndex, EndIndex: The visibility maps indices range to delete.       *
*   Extension:            The extension of the visibility maps file names.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCDeleteVisMap(int StartIndex, int EndIndex, char *Extension)
{
    char FileName[IRIT_LINE_LEN];
    int i;

    for (i = StartIndex; i <= EndIndex; i++) {
        sprintf(FileName, "%s%d.%s", USER_GC_VISMAP_FILE_NAME_BASE, i,
                Extension);
        remove(FileName);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deleting all visibilty maps from disk.                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem: The geometric covering problem to be solved.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserGCDeleteVisMap					                     M
*									     *
* KEYWORDS:                                                                  M
*   UserGCDeleteAllVisMaps				                     M
*****************************************************************************/
void UserGCDeleteAllVisMaps(UserGCProblemDefinitionStruct *Problem) 
{
    UserGCDeleteVisMap(0, UserGCGetOPsNum(Problem) - 1,
		       Problem -> DebugParams.VisMapExtension);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Saving an object to file.                                                *
*                                                                            *
* PARAMETERS:                                                                *
*   FileName: The name of the file (without extension, the extension .ibd or *
*             .itd will added depneding on IsBinary).                        *
*   PObj:     The object to save. If NULL, an empty IP_OBJ_STRING object     *
*             file is created (an ugly patch to mark empty file).            *
*   IsBinary: Whether to save it as a binary file.                           *
*   AllList : To save all the list starting with PObj.                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: FALSE if failed saving (Error message is produced by the function). *
*****************************************************************************/
static int UserGCSaveObjectToFile(const char *FileName, 
                                  IPObjectStruct *PObj,
                                  int IsBinary,
                                  int AllList)
{
    FILE *Out;
    int FilterDegen;  
    char Name[IRIT_LINE_LEN];
    IPObjectStruct 
        *PObj2 = PObj;

    if (IsBinary) {
        sprintf(Name, "%s.ibd", FileName);
    }
    else {
        sprintf(Name, "%s.itd", FileName);
    }

    Out = UserGCOpenFile(Name, TRUE, "Can't open for writing file");
    if (!Out)
        return FALSE;

    FilterDegen = IPSetFilterDegen(FALSE);
    if (PObj == NULL)
        PObj2 = IPGenSTRObject("EmptyObject");

    if (AllList) {
        IPObjectStruct *PObj3;
        
        for (PObj3 = PObj; PObj3 != NULL; PObj3 = PObj3 -> Pnext)
            IPPutObjectToFile(Out, PObj3, IsBinary);
    }
    else
        IPPutObjectToFile(Out, PObj2, IsBinary);

    if (PObj == NULL)
        IPFreeObject(PObj2);

    IPSetFilterDegen(FilterDegen);
    fclose(Out);
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loading an object from file.                                             *
*                                                                            *
* PARAMETERS:                                                                *
*   FileName: IN, The name of the file.                                      *
*   PPObj:    OUT, The returned object, NULL if it the file is empty.        *
*   IsBinary: IN, Whether to load a binary file.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: FALSE if failed loading (Error message is produced by function).    *
*****************************************************************************/
static int UserGCLoadObjectFromFile(const char *FileName,
                                    IPObjectStruct **PPObj,
                                    int IsBinary)
{
    FILE *In;
    IPObjectStruct *PObj;
    int FlattenObjects, PropagateAttrs, Handler;
    IPProcessLeafObjType ProcessLeafFunc;
    char Name[IRIT_LINE_LEN];

    if (IsBinary) {
        sprintf(Name, "%s.ibd", FileName);
    }
    else {
        sprintf(Name, "%s.itd", FileName);
    }

    In = UserGCOpenFile(Name, FALSE, "Can't open for reading file");
    if (!In)
        return FALSE;
    
    /* Disable every thing. */
    FlattenObjects = IPSetFlattenObjects(FALSE);
    PropagateAttrs = IPSetPropagateAttrs(FALSE);
    ProcessLeafFunc = IPSetProcessLeafFunc(NULL);

    Handler = IPOpenStreamFromFile(In, TRUE, IsBinary, FALSE, FALSE);
    PObj = IPGetObjects(Handler);
    IPCloseStream(Handler, TRUE);
    fclose(In);
    IPSetProcessLeafFunc(ProcessLeafFunc);
    IPSetPropagateAttrs(PropagateAttrs);
    IPSetFlattenObjects(FlattenObjects);
    
    if (!PObj) {
        IRIT_WARNING_MSG_PRINTF("Can't read file \"%s\".", Name);
        return FALSE;
    }
    PObj -> Count = 1; /* Fix an error which may give Count value 2. */

    /* UserGCSaveObjectToFile marked this object as empty. */
    if (IP_IS_STR_OBJ(PObj)) { 
        IPFreeObject(PObj);
        PObj = NULL;
    }
    *PPObj = PObj;
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Saving both geometric object and obstacle object to disk so they can     M
* later be loaded using UserGCLoadProcessedObjects.                          M
*   If an error occur it longjmp to UserGCStartOfProcess.                    M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:   The geometric covering problem to be solved.                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserGCLoadProcessedObjects, UserGCDeleteProcessedObjects                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCSaveProcessedObjects                                               M
*****************************************************************************/
void UserGCSaveProcessedObjects(UserGCProblemDefinitionStruct *Problem)
{
    int IsBinary;

    if (Problem -> DebugParams.StoreObjectsBeforeOPs)
        IsBinary = FALSE;
    else
        IsBinary = TRUE;

    if (!UserGCSaveObjectToFile(USER_GC_PROCESSED_GEO_OBJ_NAME, 
                                Problem -> GeoObj, IsBinary, FALSE))
        longjmp(UserGCStartOfProcess, 1);

    if (!UserGCSaveObjectToFile(USER_GC_PROCESSED_OBSTACLES_NAME, 
				Problem -> Obstacles, IsBinary, FALSE))
        longjmp(UserGCStartOfProcess, 1);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Loading both geometric object and obstacle object from disk.             M
* (Those file were supposed to be saved earlier by UserGCSaveProcessedObject)M
*   If an error occur it longjmp to UserGCStartOfProcess.                    M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:   IN, The geometric covering problem to be solved.              M
*   GeoObj:    OUT, The loaded geometric object. If NULL it's ignored. The   M
*              geometric object can't be empty. This will be an error.       M
*   Obstacles: OUT, The loaded obstacles. If NULL, it's ignored. Obstacles   M
*              might be empty in which case *Obstacles will be NULL.         M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserGCSaveProcessedObjects, UserGCDeleteProcessedObjects                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCLoadProcessedObjects                                               M
*****************************************************************************/
void UserGCLoadProcessedObjects(UserGCProblemDefinitionStruct *Problem,
                                IPObjectStruct **GeoObj,
                                IPObjectStruct **Obstacles)
{
    int IsBinary;

    if (Problem -> DebugParams.StoreObjectsBeforeOPs)
        IsBinary = FALSE;
    else
        IsBinary = TRUE;

    if (GeoObj)
        if(!UserGCLoadObjectFromFile(
            USER_GC_PROCESSED_GEO_OBJ_NAME, GeoObj, IsBinary) || !(*GeoObj))
            longjmp(UserGCStartOfProcess, 1);

    if (Obstacles)
        if(!UserGCLoadObjectFromFile(USER_GC_PROCESSED_OBSTACLES_NAME, 
                                     Obstacles, IsBinary))
            longjmp(UserGCStartOfProcess, 1);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loading both geometric object and obstacle object from disk.             *
* (Those file were supposed to be saved earlier by UserGCSaveProcessedObject)*
*   If an error occur it longjmp to UserGCStartOfProcess.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:   The geometric covering problem to be solved.                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct*: The geometric object connected by Pnext to the obstacles*
*                    or NULL if no obstacles exist.                          *
*****************************************************************************/
static IPObjectStruct *UserGCLoadProcessedObjects2(
                                      UserGCProblemDefinitionStruct *Problem)
{
    IPObjectStruct *GeoObj, *Obstacles;

    UserGCLoadProcessedObjects(Problem, &GeoObj, &Obstacles);
    GeoObj -> Pnext = Obstacles;
    return GeoObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deleting the geometric object and obstacles files.                       M
* (Those file were supposed to be saved earlier by UserGCSaveProcessedObject)M
*   No error is reported in case of an error.                                M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem: The geometric covering problem to be solved.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   UserGCSaveProcessedObjects, UserGCLoadProcessedObjects                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCDeleteProcessedObjects                                             M
*****************************************************************************/
void UserGCDeleteProcessedObjects(UserGCProblemDefinitionStruct *Problem)
{
    char FileName[IRIT_LINE_LEN], *Extension;

    if (Problem -> DebugParams.StoreObjectsBeforeOPs)
        Extension = "itd";
    else
        Extension = "ibd";

    sprintf(FileName, "%s.%s", USER_GC_PROCESSED_GEO_OBJ_NAME, Extension);
    remove(FileName);
    sprintf(FileName, "%s.%s", USER_GC_PROCESSED_OBSTACLES_NAME, Extension);
    remove(FileName);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Saving the given matrix together with the object (Without the list of    *
*   pnext objects).                                                          *
*   Used for debug needs.                                                    *
*                                                                            *
* PARAMETERS:                                                                *
*   FileName: The name of the file to save.                                  *
*   Mat:      The matrix to save.                                            *
*   PObj:     The object to save.                                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCSaveObjectAndMatrix(char *FileName, 
				      IrtHmgnMatType Mat, 
				      IPObjectStruct *PObj) {

    FILE *Out;
    int i, j, FilterDegen;

    if ((Out = fopen(FileName, "w")) == NULL)
        assert("can't open file" == 0);

    if (Mat != NULL) {
        fprintf(Out,"[OBJECT VIEW_MAT\n\t[MATRIX");
        for(i = 0; i <= 3; i++) {
            fprintf(Out, "\n\t\t");
            for(j = 0; j <= 3; j++)
                fprintf(Out,"%-16.14lf ", Mat[i][j]);
        }
        fprintf(Out,"]\n]\n");
    }
    FilterDegen = IPSetFilterDegen(FALSE);
    IPPutObjectToFile(Out, PObj, FALSE);
    IPSetFilterDegen(FilterDegen);
    fclose(Out);
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Exposes UserGCCreatePrspMatrix.                                          M
*   Prepare prespective matrix for the given ZAngle and XYAngle.             M
*   The prespective matrix is used in order to consider how wide is the      M
*   opening of the observation points.                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   ZAngle:  IN, The opening in the Z axis of the viewer (which is y axis in M
*                our coordinate system).                                     M
*   XYAngle: IN, The opening in the xy plane of the viewer (which is zx      M
*                plane in our coordinate system).                            M
*   PrspMat: OUT, Prespective matrix to be used when evaluating the scene    M
*                 visibility map. No change if no perspective is required.   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposeCreatePrspMatrix					     M
*****************************************************************************/
void UserGCExposeCreatePrspMatrix(IrtRType ZAngle, 
                                  IrtRType XYAngle,
                                  IrtHmgnMatType PrspMat)
{
    UserGCCreatePrspMatrix(ZAngle, XYAngle, PrspMat);
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prepare prespective matrix for the given ZAngle and XYAngle.             *
*   The prespective matrix is used in order to consider how wide is the      *
*   opening of the observation points.                                       *
*                                                                            *
* PARAMETERS:                                                                *
*   ZAngle:  IN, The opening in the Z axis of the viewer (which is y axis in *
*                our coordinate system).                                     *
*   XYAngle: IN, The opening in the xy plane of the viewer (which is zx      *
*                plane in our coordinate system).                            *
*   PrspMat: OUT, Prespective matrix to be used when evaluating the scene    *
*                 visibility map. No change if no perspective is required.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCCreatePrspMatrix(IrtRType ZAngle, 
                                  IrtRType XYAngle,
                                  IrtHmgnMatType PrspMat)
{
    MatGenUnitMat(PrspMat);

    PrspMat[0][0] = 1 / tan((XYAngle / 2) * IRIT_DEG2RAD_CNVRT);
    PrspMat[1][1] = 1 / tan((ZAngle / 2) * IRIT_DEG2RAD_CNVRT);
    PrspMat[2][2] = (USER_GC_NEAR + USER_GC_FAR)/(USER_GC_NEAR - USER_GC_FAR);
    PrspMat[3][2] = 2 * (USER_GC_NEAR * USER_GC_FAR) /
                                                 (USER_GC_NEAR - USER_GC_FAR);
    PrspMat[2][3] = -1;
    PrspMat[3][3] = 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Get the axes of the coordinate system around the given observation point.*
*                                                                            *
* PARAMETERS:                                                                *
*   Op:    IN, The observation point.                                        *
*   Up:    IN, The general direction of Up or y.                             *
*   x,y,z: OUT, The returned axes.                                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCFindNewAxes(const UserGCObsPtSuggestionStruct *Op, 
                              const IrtVecType Up,
                              IrtPtType x,
                              IrtPtType y,
                              IrtPtType z)
{
    IrtPtType MinusZ;

    assert((IRIT_PT_LENGTH(Op -> Direction) > IRIT_PT_NORMALIZE_ZERO) && 
	   "The direction of the observation point can't be (0,0,0)");
    IRIT_PT_COPY(MinusZ, Op -> Direction);
    IRIT_PT_NORMALIZE(MinusZ);

    IRIT_CROSS_PROD(x, MinusZ, Up);
    if (IRIT_PT_LENGTH(x) < IRIT_PT_NORMALIZE_ZERO) {
        /* MinusZ is parallel to Up and there is no way to decide where is  */
        /* y. I therefore decide to set x to be arbitrary vector            */
        /* perpendicular to MinusZ.                                         */
        IrtPtType
	    Temp = { 0, 0, 1 };

        IRIT_CROSS_PROD(x, MinusZ, Temp);
        if (IRIT_PT_LENGTH(x) < IRIT_PT_NORMALIZE_ZERO) {
            /* Up seems to be parallel to {0,0,1} so we select another      */
            /* creator.                                                     */
            IRIT_PT_SET(Temp, 0, 1, 0);
            IRIT_CROSS_PROD(x, MinusZ, Temp);
        }
    }
    IRIT_PT_NORMALIZE(x);
    IRIT_CROSS_PROD(y, x, MinusZ);
    IRIT_PT_NORMALIZE(y);
    IRIT_PT_COPY(z, MinusZ);
    IRIT_PT_SCALE(z, -1);
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Exposing UserGCCreateViewMatrix2.                                        M
*   Prepare view matrix for the given observation point. Using this matrix   M
*   any scene can be modifed to look the way it will be seen from the given  M
*   observation point.                                                       M
*                                                                            *
* PARAMETERS:                                                                M
*   Op:         IN, The observation point from which the scene will be seen. M
*   ViewMat:    OUT, View matrix which set modification of the scene so it   M
*               will be seen from different location or/and direction.       M
*   Up:         IN, The general direction of the Up.                         M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposeCreateViewMatrix2			            	     M
*****************************************************************************/
void UserGCExposeCreateViewMatrix2(const UserGCObsPtSuggestionStruct *Op, 
                                   IrtHmgnMatType ViewMat,
                                   const IrtVecType Up)
{
    UserGCCreateViewMatrix2(Op, ViewMat, Up);
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prepare view matrix for the given observation point. Using this matrix   *
*   any scene can be modifed to look the way it will be seen from the given  *
*   observation point.                                                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Op:         IN, The observation point from which the scene will be seen. *
*   ViewMat:    OUT, View matrix which set modification of the scene so it   *
*               will be seen from different location or/and direction.       *
*   Up:         IN, The general direction of the Up.                         *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCCreateViewMatrix2(const UserGCObsPtSuggestionStruct *Op, 
                                    IrtHmgnMatType ViewMat,
                                    const IrtVecType Up)
{
    IrtPtType x, y, z;
    IrtHmgnMatType Mat;

    if (!IRIT_PT_APX_EQ(Op -> ObsPt, USER_GC_INF_VEC)) {
        MatGenMatTrans(-Op -> ObsPt[0], -Op -> ObsPt[1], -Op -> ObsPt[2], 
		       ViewMat);
    }
    else
        MatGenUnitMat(ViewMat);

    UserGCFindNewAxes(Op, Up, x, y, z);

    MatGenUnitMat(Mat);
    Mat[0][0] = x[0];
    Mat[0][1] = y[0];
    Mat[0][2] = z[0];
    Mat[1][0] = x[1];
    Mat[1][1] = y[1];
    Mat[1][2] = z[1];
    Mat[2][0] = x[2];
    Mat[2][1] = y[2];
    Mat[2][2] = z[2];

    MatMultTwo4by4(ViewMat, ViewMat, Mat);
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*  Exposing UserGCCreateViewMatrix.                                          M
*   Prepare view matrix for the given observation point. Using this matrix   M
* any scene can be modifed to look the way it will be seen from the given    M
* observation point.                                                         M
*   Call UserGCCreateViewMatrix2 with USER_GC_UP                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Op:         IN, The observation point from which the scene will be seen. M
*   ViewMat:    OUT, View matrix which set modification of the scene so it   M
*               will be seen from different location or/and direction.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposeCreateViewMatrix				             M
*****************************************************************************/
void UserGCExposeCreateViewMatrix(const UserGCObsPtSuggestionStruct *Op, 
                                  IrtHmgnMatType ViewMat)
{
    UserGCCreateViewMatrix(Op, ViewMat);
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prepare view matrix for the given observation point. Using this matrix   *
* any scene can be modifed to look the way it will be seen from the given    *
* observation point.                                                         *
*   Call UserGCCreateViewMatrix2 with USER_GC_UP                             *
*                                                                            *
* PARAMETERS:                                                                *
*   Op:         IN, The observation point from which the scene will be seen. *
*   ViewMat:    OUT, View matrix which set modification of the scene so it   *
*               will be seen from different location or/and direction.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCCreateViewMatrix(const UserGCObsPtSuggestionStruct *Op, 
                                   IrtHmgnMatType ViewMat)
{
    UserGCCreateViewMatrix2(Op, ViewMat, USER_GC_UP);
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Expose UserGCDivideAndCreateViewMatrices.                                M
*   Call UserGCDivideOP to divide the observation point and then return      M
* view matrices for each of them.                                            M
*                                                                            *
* PARAMETERS:                                                                m
*   ObsPt:        IN, The observation point to divide.                       M
*   OpeningInXY:  IN, The input angle opening in the xz plane.               M
*   OpeningInZ:   IN, The input angle opening in the y direction.            M
*   ObsPtsNum:    OUT, Number of observation points returned.                M
*   OpeningOutXY: OUT, The output angle opening in xz plane of all returned  M
*                 angles.                                                    M
*   OpeningOutZ:  OUT, The output angle opening in xz plane of all returned  M
*                 angles.                                                    M
*   ViewMats:     OUT, View matrices for the observation points. 6 at most.  M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposeDivideAndCreateViewMatrices			             M
*****************************************************************************/
void UserGCExposeDivideAndCreateViewMatrices(UserGCObsPtSuggestionStruct
					                               *ObsPt,
                                             IrtRType OpeningInXY, 
                                             IrtRType OpeningInZ,
                                             int *ObsPtsNum,
                                             IrtRType *OpeningOutXY,
                                             IrtRType *OpeningOutZ,
                                             IrtHmgnMatType ViewMats[6])
{
    UserGCDivideAndCreateViewMatrices(ObsPt, OpeningInXY, OpeningInZ, 
                                      ObsPtsNum, OpeningOutXY, OpeningOutZ,
                                      ViewMats);
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call UserGCDivideOP to divide the observation point and then return      *
* view matrices for each of them.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   ObsPt:        IN, The observation point to divide.                       *
*   OpeningInXY:  IN, The input angle opening in the xz plane.               *
*   OpeningInZ:   IN, The input angle opening in the y direction.            *
*   ObsPtsNum:    OUT, Number of observation points returned.                *
*   OpeningOutXY: OUT, The output angle opening in xz plane of all returned  *
*                 angles.                                                    *
*   OpeningOutZ:  OUT, The output angle opening in xz plane of all returned  *
*                 angles.                                                    *
*   ViewMat:      OUT, View matrices for the observation points. 6 at most.  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCDivideAndCreateViewMatrices(UserGCObsPtSuggestionStruct
					                               *ObsPt,
                                              IrtRType OpeningInXY, 
                                              IrtRType OpeningInZ,
                                              int *ObsPtsNum,
                                              IrtRType *OpeningOutXY,
                                              IrtRType *OpeningOutZ,
                                              IrtHmgnMatType ViewMats[6])
{
    IrtRType AnglesZ[3], Dummy2[3];
    IrtHmgnMatType InvMat;
    IrtPtType x, y, z;
    UserGCObsPtSuggestionStruct ObsPts[9];
    int i, ObsPtXYNum;

    /* Get all max 9 observation points. */
    UserGCDivideOP(ObsPt, OpeningInXY, OpeningInZ, ObsPtsNum, ObsPts, 
                   OpeningOutXY, OpeningOutZ);

    /* Get the observation points when rotating only on Z axis. That means  */
    /* only the division on z axis.                                         */
    /* I failed working both XY and Z together so I gave up on Z and for Z  */
    /* I use one angle up to 175 degrees.                                   */
    //UserGCDivideAngle(OpeningInZ, &Dummy, AnglesZ, Dummy2);
    assert(OpeningInZ <= 175);
    AnglesZ[0] = 0;

    /* Get the observation points when rotating only on XY plane. That      */
    /* means only the division on XY plane.                                 */
    UserGCDivideAngle(OpeningInXY, &ObsPtXYNum, Dummy2, Dummy2);

    /* We first calculate the division of vectors relative to the vector    */
    /* {0,0,-1}. Then we assume that {0,0,-1} is actually ObsPt.Direction   */
    /* and what we have is actually a view relative to ObsPt.Direction. So  */
    /* we need to rotate it so it will be relative to the scene. Usually we */
    /* use UserGCCreateViewMatrix to rotate a scene to be relative to       */
    /* ObsPt.direction. So the inverse of the created matrix will rotate    */
    /* ObsPt.direction back to be relative to the original scene.           */
    UserGCCreateViewMatrix(ObsPt, InvMat); 
    MatInverseMatrix(InvMat, InvMat); 
    
    for (i = 0; i < *ObsPtsNum; i++) {        
        if (i % ObsPtXYNum == 0) {
            IrtHmgnMatType MatZ;
            UserGCObsPtSuggestionStruct TempObsPt;
            IrtPtType
	        Pt = { 0, 0, -1 };

            MatGenMatRotX1(M_PI /180 * AnglesZ[i / ObsPtXYNum], MatZ);
            IRIT_GEN_COPY(&TempObsPt, ObsPt, 
                          sizeof(UserGCObsPtSuggestionStruct));
            MatMultVecby4by4(Pt, Pt, MatZ);/* Rotates Pt to the direction   */
                                           /* relatvie to ObsPt.Direction.  */
            /* Rotates Pt to be relative to the original scene. */
            MatMultVecby4by4(TempObsPt.Direction, Pt, InvMat);
            UserGCFindNewAxes(&TempObsPt, USER_GC_UP, x, y, z);
        }
        UserGCCreateViewMatrix2(&ObsPts[i], ViewMats[i], y);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find the intersection between a 3D line and a plane.                     *
*   If "axis" is x or y the input plane is perpendicular to the plane made   *
* of z and "axis". The the angle between the input plane and the negative z  *
* is "Param".                                                                *
*   If "axis" is z the input plane is prependicular to the z axis and the    *
* input plane intersect the z axis on coordinate "Param".                    *
*   The returned Blend is such that the intersection is Pt0 + Blend*(Pt1-Pt0)*
*                                                                            *
* PARAMETERS:                                                                *
*   Pt0, Pt1:  IN, Define the line.                                          *
*   Axis:      IN, if x or y then the input plane is perpendicular to the    *
*                  plane made of "Axis" and z.                               *
*                  If it is z then the input plane is perpendicular to the z *
*                  axis.                                                     *
*   Param:     IN, If "Axis" is x or y this is the angle of the input        *
*                  plane relative to negative z. If "Axis" is z this is the  *
*                  intersection of the input plane on the z axis.            *
*   Blend:     OUT, Recieved from the formula of the intersection point      *
*                  Pt0 + Blend*(Pt1-Pt0).                                    *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: FALSE if there is no intersection.                                  *
*****************************************************************************/
static int UserGCGetIntersection(IrtPtType Pt0,
				 IrtPtType Pt1,
				 UserGCAxesType Axis,
				 IrtRType Param,
				 IrtRType *Blend)
{
    IrtRType X0, Y0, X1, Y1, A0, B0, A1, B1, x, y;

    if (Axis == USER_GC_Z_AXIS) {
        if (IRIT_APX_EQ(Pt0[USER_GC_Z_AXIS], Pt1[USER_GC_Z_AXIS]))
            return FALSE;
        *Blend = (Param - Pt0[USER_GC_Z_AXIS])/
            (Pt1[USER_GC_Z_AXIS] - Pt0[USER_GC_Z_AXIS]);
        return TRUE;
    }

    /* x axis or y axis. */
    X0 = Pt0[Axis];
    Y0 = Pt0[USER_GC_Z_AXIS];
    X1 = Pt1[Axis];
    Y1 = Pt1[USER_GC_Z_AXIS];
    A1 = tan(IRIT_DEG2RAD(Param - 90));
    B1 = 0;
    if (!IRIT_APX_EQ(X1, X0)) {
        A0 = (Y1 - Y0)/(X1 - X0);
        B0 = Y0 - A0 * X0;
        if (IRIT_APX_EQ(A0, A1))
            return FALSE;    
        x = (B1 - B0)/(A0 - A1);
        *Blend = (x - X0)/(X1 - X0);
    }
    else {
        x = X1;
        y = A1 * x + B1;
        if (IRIT_APX_EQ(Y1, Y0))
            return FALSE;
        *Blend = (y - Y0)/(Y1 - Y0);
    }
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Crop the polygon Polygon by one of the planes of the trapeze of the      *
* perspective view.                                                          *
*   Polygon is assumed to be on one plane.                                   *
*   Polygon isn't touched at all.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   Polygon:   IN, The polygon to crop.                                      *
*   Axis:      IN, The axis which together with IsSmaller decides which side *
*              of PolyPlane should be kept.                                  *
*   IsSmaller: IN, If TRUE leaves the part of the polygon which is           *
*              below/left/farther. If FALSE leaves the part of the polygon   *
*              which is above/right/closer. It Depends on Axis.              *
*   Param:     For x and y axis it's the angle of the plane (Notice that for *
*              perspective view with opening angle A the angle expected here *
*              is either A/2 or -A/2). For z axis it's the z value of the    *
*              near/far planes.                                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPPolygonStruct*: The created polygon or or a the original               *
*                     polygon if no change was done or NULL if the polygon   *
*                     is entirely on the wrong side of PolyPlane.            *
*****************************************************************************/
static IPPolygonStruct *UserGCCropPolygon(IPPolygonStruct *Polygon,
                                          UserGCAxesType Axis,
                                          int IsSmaller,
                                          IrtRType Param) 
{
    IPPolygonStruct *NewPolygon;
    IPVertexStruct *V, 
        *PrevV = Polygon -> PVertex,
        *NewListV = IPAllocVertex2(NULL); /* Dummy vertex. */
    int Sign = IsSmaller ? -1 : 1,
        In = 0,
        EntirelyIn = TRUE;
    UserGCAxesType
        Coord0 = Axis;

    /* Closes a circle and looks at the last edge as well. */
    IPGetLastVrtx(Polygon -> PVertex) -> Pnext = 
                                            IPCopyVertex(Polygon -> PVertex); 
    AttrSetIntAttrib(&IPGetLastVrtx(Polygon -> PVertex) -> Attr, 
		     USER_GC_CLOSING_CIRCLE, TRUE);

    for (V = Polygon -> PVertex; V != NULL; V = V -> Pnext) {
        int NewIn;
        float *Uv[2], NewUv[2];
        IrtRType PlaneAxisValue;

        if (Axis == USER_GC_Z_AXIS) 
            PlaneAxisValue = Param;
        else 
            PlaneAxisValue = tan(IRIT_DEG2RAD(Param)) *
	                                       (-V -> Coord[USER_GC_Z_AXIS]);
        NewIn = V -> Coord[Coord0] * Sign >= PlaneAxisValue * Sign;
        /* We crossed the limit. Creates new vertex. */
        if ((V != Polygon -> PVertex) && (NewIn != In)) {
            IPVertexStruct 
                *NewV = IPAllocVertex2(NULL);
            IrtRType
		Factor = 0.0;

            /* We know in this place that there must be intersection. So no */
            /* need to check the result.                                    */
            UserGCGetIntersection(PrevV -> Coord, V -> Coord, Axis, Param, 
                                  &Factor);
            IRIT_PT_BLEND(NewV -> Coord, V -> Coord, PrevV -> Coord, Factor);
            Uv[0] = AttrGetUVAttrib(PrevV -> Attr, USER_GC_UV_VALUES);
            Uv[1] = AttrGetUVAttrib(V -> Attr, USER_GC_UV_VALUES);
            if ((Uv[0] != NULL) && (Uv[1] != NULL)) {
                NewUv[0] = Uv[0][0] + 
                    (float)Factor*(Uv[1][0] - Uv[0][0]);
                NewUv[1] = Uv[0][1] + 
                    (float)Factor*(Uv[1][1] - Uv[0][1]);
                AttrSetUVAttrib(&NewV -> Attr, USER_GC_UV_VALUES, 
                    NewUv[0], NewUv[1]);
            }
            IPGetLastVrtx(NewListV) -> Pnext = NewV;
        }
        if (NewIn) 
            IPGetLastVrtx(NewListV) -> Pnext = IPCopyVertex(V);
        else 
            EntirelyIn = FALSE;
        PrevV = V;
        In = NewIn;
    }

    /* Remove the last vertex previously added to Polygon. Break the circle.*/
    V = IPGetLastVrtx(Polygon -> PVertex);
    IPGetPrevVrtx(Polygon -> PVertex,V) -> Pnext = NULL;
    IPFreeVertex(V);

    /* The polygon is entirely on the wrong side. */
    if (NewListV -> Pnext == NULL) {
        IPFreeVertex(NewListV);
        return NULL;
    }

    /* Remove last vertex if it equals to the first. Break the circle. */
    V = IPGetLastVrtx(NewListV);

    if (AttrGetIntAttrib(V -> Attr,
			 USER_GC_CLOSING_CIRCLE) != IP_ATTR_BAD_INT) {
        IPGetPrevVrtx(NewListV,V) -> Pnext = NULL;
        IPFreeVertex(V);
    }

    /* The polygon is entirely in the good side. */
    if (EntirelyIn) {
        IPFreeVertexList(NewListV);
        return Polygon;
    }

    /* Copy all attributes and other stuff from the original polygon. */
    NewPolygon = IPCopyPolygon(Polygon);
    IPFreeVertexList(NewPolygon -> PVertex);
    NewPolygon -> PVertex = NewListV -> Pnext;
    IPFreeVertex(NewListV);
    return NewPolygon;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Crops a triangle if it is partially outside the trapeze of the           *
* perspective view defined by AngleXY, AngleZ, Near and Far.                 *
                                                                             *
* 3 possible outcomes:                                                       *
* - The triangle is entirely in the trapeze. The return value will the       *
*   original Triangle.                                                       *
* - The triangle is partially in the trapeze. The return value will be a list*
*   of triangles which are equal to the part of the original Triangle which  *
*   was inside the trapeze. The original Triangle isn't changed.             *
* - The triangle is entirely outside the box. The return value is NULL and   *
*   the original Triangle isn't changed.                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   Triangle:        The triangle to crop/purge away.                        *
*   AngleXY, AngleZ: The opening of the perspective view.                    *
*   Near, Far:       The near and far planes of the perspective view.        *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPPolygonStruct*: The created triangles or the original triangle if no   *
*                     change was done or NULL if the triangle is outside the *
*                     box.                                                   *
*****************************************************************************/
static IPPolygonStruct *UserGCCropTriangle( IPPolygonStruct *Triangle,
                                            IrtRType AngleXY,
                                            IrtRType AngleZ,
                                            IrtRType Near,
                                            IrtRType Far)
{
    IPPolygonStruct *Pl1, *Pl2, *TemplateTriangle;
    IPVertexStruct *V;
    int i;
    IrtRType Params[6] = {
        -AngleXY / 2,
        AngleXY / 2,
        -AngleZ / 2,
        AngleZ / 2,
        -Far,
        -Near,
    };

    Pl1 = Triangle;
    for (i = 0; i < 6; i++) { 
        Pl2 = UserGCCropPolygon(Pl1, (UserGCAxesType) (i / 2), 
				i % 2 == 1, Params[i]);
        if (Pl2 == NULL) {
            if (Pl1 != Triangle)
                IPFreePolygon(Pl1);
            return NULL;
        } 
        if ((Pl2 != Pl1) && (Pl1 != Triangle))
            IPFreePolygon(Pl1);
        Pl1 = Pl2;
    }

    /* Turn the polygon into triangles. */
    TemplateTriangle = IPCopyPolygon(Pl1);

    /* We create triangles starting from the second triangle which means the*/
    /* third and fourth vertices. The first trinagle is Pl1 itself  after   */
    /* leaving only its 3 first vertices.                                   */
    for (V = Pl1 -> PVertex -> Pnext -> Pnext;
	 V -> Pnext != NULL; 
	 V = V -> Pnext) {
        IPVertexStruct 
            *NewV = IPCopyVertex(Pl1 -> PVertex);
        IPPolygonStruct *Pl;

        NewV -> Pnext = IPCopyVertex(V);
        NewV -> Pnext -> Pnext = IPCopyVertex(V -> Pnext);
        Pl = IPCopyPolygon(TemplateTriangle);
        Pl -> PVertex = NewV;
        IPGetLastPoly(Pl1) -> Pnext = Pl;
    }

    /* Leaving only the first 3 vertices of NewTriangle. All the rest were  */
    /* added to the triangles above.                                        */
    IPFreeVertexList(Pl1 -> PVertex -> Pnext ->Pnext -> Pnext);
    Pl1 -> PVertex -> Pnext -> Pnext -> Pnext = NULL;
    IPFreePolygon(TemplateTriangle);

    return Pl1;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Transform Scene to its final stage (ViewMat, PrspMat and ScreenMat).     *
* Crop triangles for perspective view.                                       *
* This will modify Scene and it should be used again unless you know what you*
* are doing.                                                                 *
*                                                                            *
* PARAMETERS:                                                                *
*   Rend:        The render context.                                         *
*   Scene:       The scene to scan. Assuming it contain objects only in the  *
*                Scene->Pnext list and not the Scene->U.Lst list.            *
*                Should contain GeoObj followed on Pnext by obstacles or     *
*                NULL if no obstacle exist.                                  *
*   ViewMat:     The view matrix.                                            *
*   PrspMat:     Prespective matrix to be used when evaluating the scene     *
*                visibility map. NULL if no perspective is used (In which    *
*                case AngleXY and AngleZ are ignored as well. Also, no       *
*                cropping is done. In orthographic projection all the object *
*                is in the scene.)                                           *
*   AngleXY:     The Xy opening of the perspective view.                     *
*   AngleZ:      The Z opening of the perspective view.                      *
*   PicSuffix:   The name of this picture to be used as suffix in files saved*
*                for debug purposes. Should probably be made of the image    *
*                index.                                                      *
*   SaveObject:  Whether to save the object before rendering. For debug      *
*                purposes only.                                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCTransformAndCrop(IRndrPtrType Rend,
                                  IPObjectStruct *Scene,
                                  IrtHmgnMatType ViewMat,
			          IrtHmgnMatType PrspMat,
                                  IrtRType AngleXY,
                                  IrtRType AngleZ,
                                  const char *PicSuffix,
                                  int SaveObject)
{
    int i;
    IPPolygonStruct *Poly;
    IPObjectStruct *TempObstacles;
    char NameAfterViewMat[50], NameFinal[50];
    IrtHmgnMatType ScrnMat;

    /* We first set all matrices so that IRndrVisMapEnable's calculations   */
    /* will be correct. We need to do it before we starting our             */
    /* transformation because IRndrVisMapEnable look both at the Scene and  */
    /* at the matrices in Rend. We don't want double use of the matrices.   */
    IRndrSetViewPrsp(Rend, ViewMat, PrspMat, NULL);

    /* IRndrVisMapEnable calculates uv limits and xy limits. uv limits      */
    /* mustn't be calculated to obstacles. xy limits are important only for */
    /* mold which doesn't have any obstacles at all. Therefore we remove    */
    /* obstacles temporary.                                                 */
    TempObstacles = Scene -> Pnext; 
    Scene -> Pnext = NULL;
    IRndrVisMapEnable(Rend, Scene, 1, TRUE);
    Scene -> Pnext = TempObstacles;

    if (SaveObject) {
        sprintf(NameAfterViewMat, "GeoObjectAfterViewMat%s.itd", PicSuffix);
        sprintf(NameFinal, "GeoObjectFinal%s.itd", PicSuffix);
    }

    for (i = 0; i < 2; i++, Scene = Scene -> Pnext) {
        /* First doing only the ViewMat transformation. */
        for (Poly = Scene -> U.Pl; Poly != NULL; Poly = Poly -> Pnext) {
            IPVertexStruct *Vertex;
            
            for (Vertex = Poly -> PVertex;
		 Vertex != NULL;
		 Vertex = Vertex -> Pnext) {                
                MatMultPtby4by4(Vertex -> Coord, Vertex -> Coord, ViewMat);
            }            
        }
        if (PrspMat) {
            /* Cropping the triangle. */
            for (Poly = Scene -> U.Pl;
		 Poly != NULL;
		 Poly = Poly -> Pnext) {
                IPPolygonStruct *TempPoly;

                TempPoly = UserGCCropTriangle(Poly, AngleXY, AngleZ, 
                    USER_GC_NEAR, USER_GC_FAR);
                if (TempPoly == NULL) { /* Triangle is entirely outside. */
                    /* We mark the polygon as purged away so that in        */
                    /* IRndrPutTriangle it will be handled only by UV scan. */
                    AttrSetIntAttrib(&Poly -> Attr, IRNDR_ATTRIB_PURGED_POLY, 
                        TRUE);
                }
                else if (TempPoly != Poly) {/*Trianlge is partially outside.*/
                    IPPolygonStruct *Pl;

                    /* We mark the polygon as purged away so that in        */
                    /* IRndrPutTriangle it will be handled only by UV scan. */
                    /* It will be marked as not-rendered. This will put red */
                    /* color in all were it's rendered. However,            */
                    /* not-rendered is overwritten by any other value, so   */
                    /* the small triangles which come after it will         */
                    /* overwrite this values wherever they rendered, leaving*/
                    /* only out of the scop pixel to be not-rendered.       */
                    AttrSetIntAttrib(&Poly -> Attr, IRNDR_ATTRIB_PURGED_POLY, 
                        TRUE);
                    Pl = IPGetLastPoly(TempPoly);
                    Pl -> Pnext = Poly -> Pnext;
                    Poly -> Pnext = TempPoly;
                    Poly = Pl;
                }
            }
        }

        /* Debug issue. Saving the object. */
        if (SaveObject) 
            UserGCSaveObjectAndMatrix(NameAfterViewMat, NULL, Scene);

        /* For IRndrVertexTransform we need the matrices to be without the  */
        /* view matrix because we already applied the view matrices above.  */
        /* However, the current ScrnMat is the result of previous call to   */
        /* IRndrVisMapEnable and it should not be changed. So we first      */
        /* acquire it and then use it in later setting of the matrices.     */
        IRndrGetViewPrsp(Rend, NULL, NULL, ScrnMat);

        /* Doing the perspective transformation (if exist) and the screen   */
        /* transformation.                                                  */
        IRndrSetViewPrsp(Rend, NULL, PrspMat, ScrnMat);
        for (Poly = Scene -> U.Pl; Poly != NULL; Poly = Poly -> Pnext) {
            IPVertexStruct *Vertex;
            IrtRType Coord[4];

            for (Vertex = Poly -> PVertex;
		 Vertex != NULL;
		 Vertex = Vertex -> Pnext) {                

                /* In perspective mode we keep the coordinates before the   */
                /* perspective and screen transformation to be used in the  */
                /* UV rendering.                                            */
                if (PrspMat) {
                    AttrSetRealAttrib(&Vertex -> Attr,
                                      VIS_MAP_BEFORE_PRSP_X_ATTRIB, 
                                      Vertex -> Coord[0]);
                    AttrSetRealAttrib(&Vertex -> Attr, 
                                      VIS_MAP_BEFORE_PRSP_Y_ATTRIB, 
                                      Vertex -> Coord[1]);
                    AttrSetRealAttrib(&Vertex -> Attr, 
                                      VIS_MAP_BEFORE_PRSP_Z_ATTRIB, 
                                      Vertex -> Coord[2]);
                }

                IRndrVertexTransform(Rend, Vertex, Coord);
                IRIT_PT_COPY(Vertex -> Coord, Coord);
                if (PrspMat)
		    AttrSetRealAttrib(&Vertex ->  Attr, "_1/W", Coord[3]);
            }            
        }
        /* This will make sure rndr_lib won't do the transformation again. */
        AttrSetObjectIntAttrib(Scene, "_TRANSFORMED", TRUE);

        /* Now, we set all the matrices again to their right values for    */
        /* correct use in the rest of the algorithm such as for backface   */
        /* culling. We also use ScrnMat which contain the changes done by  */
        /* IRndrVisMapEnable.                                              */
        IRndrSetViewPrsp(Rend, ViewMat, PrspMat, ScrnMat);

        /* Debug issue. Saving the object. */
        if (SaveObject) 
            UserGCSaveObjectAndMatrix(NameFinal, NULL, Scene);

        if (Scene -> Pnext == NULL)
            break;

        /* Preparing for obstacle file which will be treated the next       */
        /* iteration.                                                       */
        if (SaveObject) {
            sprintf(NameAfterViewMat, "ObstaclesAfterViewMat%s.itd", 
                    PicSuffix);
            sprintf(NameFinal, "ObstaclesFinal%s.itd", PicSuffix);
        }
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Scane all the polygons in the given scene into the visibility map.       *
*   In the grand plan the objects suppose to already be after transformation *
* in this stage so no more transformation is done here. However, some        *
* changes may still happen to Scene. Don't use Scene again unless you know   *
* what you are doing (free it instead).                                      *
*                                                                            *
* PARAMETERS:                                                                *
*   Rend:        The render context.                                         *
*   Scene:       The scene to scan. Assuming it contain objects only in the  *
*                Scene->Pnext list and not the Scene->U.Lst list.            *
*                Should contain GeoObj followed on Pnext by obstacles or     *
*                NULL if no obstacle exist.                                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCScanObjects(IRndrPtrType Rend,
                              IPObjectStruct *Scene)
{
    int i;

    IRndrVisMapSetScanOnUV(Rend, TRUE);
    for (i = 0; i < 2; i++, Scene = Scene -> Pnext) {
        IPPolygonStruct *Poly;

        assert(IP_IS_POLY_OBJ(Scene) && IP_IS_POLYGON_OBJ(Scene));
        assert((Scene != NULL) && (Scene -> U.Pl != NULL));

        IRndrBeginObject(Rend, Scene, FALSE);
        for (Poly = Scene -> U.Pl;
	     Poly != NULL;
	     Poly = Poly -> Pnext) {
            IRndrPutTriangle(Rend, Poly);
        }
        IRndrEndObject(Rend);

        if (Scene -> Pnext == NULL)
            break;
        /* Preparing for obstacle file which will be treated the next       */
        /* iteration.                                                       */
        IRndrVisMapSetScanOnUV(Rend, FALSE);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create visibility map to the given scene with the given view and         *
* perspective matrices. In case PrspMat isn't NULL the vismap will see only  *
* objects between USER_GC_NEAR and USER_GC_FAR. The perspective projection   *
* itself is expected to be received in PrspMat.                              *
* Scene is freed during the process.                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:    The geometric covering problem to be solved.                 *
*   Scene:      The scene to scan. Freed during the process.                 *
*   ViewMat:    View matrix to be used when evaluating the scene visibility  *
*               map.                                                         *
*   PrspMat:    Prespective matrix to be used when evaluating the scene      *
*               visibility map. NULL if no perspective is used (In which case*
*               AngleXY and AngleZ are ignored as well).                     *
*   AngleXY:    The Xy opening of the perspective view.                      *
*   AngleZ:     The Z opening of the perspective view.                       *
*   PicSuffix:  The name of this picture to be used as suffix in files saved *
*               for debug purposes. Should probably be made of the image     *
*               index.                                                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtImgPixelStruct*: The visibility map of the given scene from the given *
*                       observation point.                                   *
*****************************************************************************/
static IrtImgPixelStruct *UserGCCreateVisMap(
				       UserGCProblemDefinitionStruct *Problem,
				       IPObjectStruct *Scene, 
                                       IrtHmgnMatType ViewMat,
				       IrtHmgnMatType PrspMat,
                                       IrtRType AngleXY,
                                       IrtRType AngleZ,
                                       const char *PicSuffix)
{
    IRndrPtrType Rend;
    /* Actually, since I'm using vismap the background color is overwritten */
    /* in IRndrInitialize and it's meaningless here.                        */
    IRndrColorType
        BackGround = { 255, 255, 255 };
    IrtImgPixelStruct *Result;
    /* In mold we use backface culling because this is a closed surface.    */
    /* Notice that this is applied only to the xyz scan. The uvz scan is    */
    /* always with backface culling.                                        */
    int BackFace = (PrspMat == NULL) ? TRUE : FALSE;

    /* Notice that I don't initialize background color because for          */
    /* visibility map it's initialized inside IRndrInitialized.             */
    Rend = IRndrInitialize(Problem -> SolvingParams.VisMapWidth,
                           Problem -> SolvingParams.VisMapHeight,
                           1, FALSE, TRUE, BackFace, BackGround , 0.0, TRUE);

    /* Set callback function that write the image to a buffer instead of    */
    /* a file.                                                              */
    IRndrSaveFileCB(Rend, UserGCImgWriteSetType, UserGcImgWriteOpenFile, 
		    UserGCImgWritePutLine, UserGCImgWriteCloseFile);
    IRndrSetShadeModel(Rend, IRNDR_SHADING_NONE);

    /* Must be before setting the matrices. */
    UserGCTransformAndCrop(Rend, Scene, ViewMat, PrspMat, AngleXY, AngleZ, 
                           PicSuffix, 
                           Problem -> DebugParams.StoreObjectsForOPs);

    /* Setting the TanAng, CriticAr and UvDilation. Must be after the       */
    /* initialization in IRndrVisMapEnable which is done inside             */
    /* UserGCTransformAndCrop.                                              */
    IRndrVisMapSetTanAngle(Rend,
		           cos((90 - Problem -> SolvingParams.VisMapTanAng)*
                           (M_PI / 180)));
    IRndrVisMapSetCriticAR(Rend, 
			   Problem -> SolvingParams.VisMapCriticAR);
    UserGCScanObjects(Rend, Scene);

    /* IRndrVisMapScan uses pointers which point into Scene. Therefore,     */
    /* Scene mustn't be freed before it.                                    */
    IRndrVisMapScan(Rend);
    IPFreeObjectList(Scene);
    IRndrSaveFileVisMap(Rend, ".", "dummy.ppm", "ppm");
    IRndrDestroy(Rend);
    Result = UserGCVisMap;
    UserGCVisMap = NULL;
    return Result;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Divide the angle OpeningIn into minimum number of angles covering the    *
* same OpeninigIn angle. Each angle is at most 120. This mean at most 3      *
* angles are returned.                                                       *
*   If OpeningIn is between 120 and 120 + USER_GC_ANGLE_OVARLAP one angle    *
* greater than 120 is still returned.                                        *
*   If OpeningIn is 360 (exactly 360) then there is overlap between the start*
* of the first opening and the end of the last opening.                      *
*                                                                            *
* PARAMETERS:                                                                *
*   OpeningIn:  IN, The input angle opening.                                 *
*   AnglesNum:  OUT, Number of angles returned.                              *
*   Angles:     OUT, The returned angles (3 at most).                        *
*   OpeningOut: OUT, The output angle opening of all returned angles.        *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCDivideAngle(IrtRType OpeningIn, 
                              int *AnglesNum, 
                              IrtRType Angles[3], 
                              IrtRType *OpeningOut) 
{
    IrtRType Step;

    assert(AnglesNum && Angles && OpeningOut && (OpeningIn <= 360) &&
	   (OpeningIn > 0));

    if (OpeningIn <= USER_GC_MAX_ANGLE + USER_GC_ANGLE_OVARLAP) {
        *AnglesNum = 1;
        Angles[0] = 0;
        *OpeningOut = OpeningIn;
        return;
    }

    *AnglesNum = (int) ceil(OpeningIn/USER_GC_MAX_ANGLE);

    if (OpeningIn == 360) /* Doing overlap between -0 to +360. */
        OpeningIn += USER_GC_ANGLE_OVARLAP;

    *OpeningOut = (OpeningIn + (*AnglesNum - 1)*USER_GC_ANGLE_OVARLAP) / 
                  (*AnglesNum);
    Angles[0] = (-OpeningIn / 2) + (*OpeningOut) / 2;
    Step = *OpeningOut - USER_GC_ANGLE_OVARLAP;
    Angles[1] = Angles[0] + Step;
    if (*AnglesNum > 2)
        Angles[2] = Angles[1] + Step; 
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Divide the observatiohn point into minimum number of observation points  M
* covering the same OpeningInXY and OpeningInZ. Each observation point cover M
* at most 120 degrees in XY. On Z there is no division and up to 175 degrees M
* are accepted. This mean at most 6 observation points are returned.         M
* If OpeningInXY is between 120 and 120+USER_GC_ANGLE_OVARLAP one            M
* observation point greater than 120 is still returned in that direction.    M
*                                                                            *
* PARAMETERS:                                                                M
*   ObsPt:        IN, The observation point to divide.                       M
*   OpeningInXY:  IN, The input angle opening in the xz plane.               M
*   OpeningInZ:   IN, The input angle opening in the y direction.            M
*   ObsPtsNum:    OUT, Number of observation points returned.                M
*   ObsPts:       OUT, The returned observation points (6 at most).          M
*   OpeningOutXY: OUT, The output angle opening in xz plane of all returned  M
*                 angles.                                                    M
*   OpeningOutZ:  OUT, The output angle opening in xz plane of all returned  M
*                 angles.                                                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposeDivideOP						     M
*****************************************************************************/
void UserGCExposeDivideOP(UserGCObsPtSuggestionStruct *ObsPt,
                          IrtRType OpeningInXY, 
                          IrtRType OpeningInZ,
                          int *ObsPtsNum,
                          UserGCObsPtSuggestionStruct *ObsPts,
                          IrtRType *OpeningOutXY,
                          IrtRType *OpeningOutZ)
{
    UserGCDivideOP(ObsPt, OpeningInXY, OpeningInZ, ObsPtsNum, ObsPts, 
                   OpeningOutXY, OpeningOutZ);
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Divide the observatiohn point into minimum number of observation points  *
* covering the same OpeningInXY and OpeningInZ. Each observation point cover *
* at most 120 degrees in XY. On Z there is no division and up to 175 degrees *
* are accepted. This mean at most 6 observation points are returned.         *
* If OpeningInXY is between 120 and 120+USER_GC_ANGLE_OVARLAP one            *
* observation point greater than 120 is still returned in that direction.    *
*                                                                            *
* PARAMETERS:                                                                *
*   ObsPt:        IN, The observation point to divide.                       *
*   OpeningInXY:  IN, The input angle opening in the xz plane.               *
*   OpeningInZ:   IN, The input angle opening in the y direction.            *
*   ObsPtsNum:    OUT, Number of observation points returned.                *
*   ObsPts:       OUT, The returned observation points (6 at most).          *
*   OpeningOutXY: OUT, The output angle opening in xz plane of all returned  *
*                 angles.                                                    *
*   OpeningOutZ:  OUT, The output angle opening in xz plane of all returned  *
*                 angles.                                                    *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCDivideOP(UserGCObsPtSuggestionStruct *ObsPt,
                           IrtRType OpeningInXY, 
                           IrtRType OpeningInZ,
                           int *ObsPtsNum,
                           UserGCObsPtSuggestionStruct *ObsPts,
                           IrtRType *OpeningOutXY,
                           IrtRType *OpeningOutZ)
{
    IrtRType AnglesXY[3], AnglesZ[3];
    int ObsPtZNum, ObsPtXYNum, z, XY, 
        Index = 0;
    IrtHmgnMatType InvMat;

    /* We first calculate the division of vectors relative to the vector    */
    /* {0,0,-1}. Then we assume that {0,0,-1} is actually ObsPt.Direction   */
    /* and what we have is actually a view relative to ObsPt.Direction. So  */
    /* we need to rotate it so it will be relative to the scene. Usually we */
    /* use UserGCCreateViewMatrix to rotate a scene to be relative to       */
    /* ObsPt.direction. So the inverse of the created matrix will rotate    */
    /* ObsPt.direction back to be relative to the original scene.           */
    UserGCCreateViewMatrix(ObsPt, InvMat); 
    MatInverseMatrix(InvMat, InvMat); 

    /* Divide the space to 120 degrees angles with 5 degrees of overlap.    */
    /* I failed working both XY and Z together so I gave up on Z and for Z  */
    /* I use one angle up to 175 degrees.                                   */
    //UserGCDivideAngle(OpeningInZ, &ObsPtZNum, AnglesZ, OpeningOutZ);
    assert(OpeningInZ <= 175);
    *OpeningOutZ = OpeningInZ;
    AnglesZ[0] = 0;
    ObsPtZNum = 1;

    UserGCDivideAngle(OpeningInXY, &ObsPtXYNum, AnglesXY, OpeningOutXY);
    *ObsPtsNum = ObsPtZNum*ObsPtXYNum;

    for (z = 0; z < ObsPtZNum; z++) {
        IrtHmgnMatType MatZ;

        MatGenMatRotX1(M_PI / 180 * AnglesZ[z], MatZ);
        for (XY = 0; XY < ObsPtXYNum; XY++) {
            IrtHmgnMatType MatXY, Mat;
            IrtPtType
	        Pt = { 0, 0, -1 };

            MatGenMatRotY1(M_PI / 180 * AnglesXY[XY], MatXY);
            MatMultTwo4by4(Mat, MatZ, MatXY);

            IRIT_GEN_COPY(&ObsPts[Index], ObsPt, 
			  sizeof(UserGCObsPtSuggestionStruct));
            MatMultVecby4by4(Pt, Pt, Mat); /* Rotates Pt to the direction   */
                                           /* relatvie to ObsPt.Direction.  */
            /* Rotates Pt to be relative to the original scene. */
            MatMultVecby4by4(ObsPts[Index].Direction, Pt, InvMat);
            Index++;
        }
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Combine the two visibility maps into one.                                *
*   If both visibility maps are NULL, return NULL. If only one is NULL       *
* returns the other unchanged.                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   VisMap1: IN OUT, The first visibility map. The visibility maps are       *
*                    combined into this visibility map.                      *
*   VisMap2: IN, The second visibility map.                                  *
*   Width:   IN, The width of the visibility map.                            *
*   Height:  IN, The height of the visibility map.                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtImgPixelStruct*: The combined visibility map which is VisMap1 (unless *
*                       VisMap1 is NULL and VisMap2 isn't in which case      *
*                       VisMap2 is returned).                                *
*****************************************************************************/
static IrtImgPixelStruct *UserGCCombineVisMaps(IrtImgPixelStruct *VisMap1, 
                                               IrtImgPixelStruct *VisMap2,
                                               int Width,
                                               int Height)
{
    int i, j,
        colorSize = sizeof(IrtImgPixelStruct);
    if (!VisMap1 && !VisMap2)
        return NULL;
    if (!VisMap1)
        return VisMap2;
    if (!VisMap2)
        return VisMap1;

    for (i = 0; i < Height; i++) {
        for (j = 0; j < Width; j++) {
            int Index = i * Width + j;
            /* If any is visible, the combination is visible. */
            if ((IRIT_GEN_CMP(&VisMap1[Index], &RNDR_VISMAP_VISIBLE_COLOR, 
			      colorSize) == 0) ||
                (IRIT_GEN_CMP(&VisMap2[Index], &RNDR_VISMAP_VISIBLE_COLOR, 
			      colorSize) == 0)) {
	        IRIT_GEN_COPY(&VisMap1[Index], &RNDR_VISMAP_VISIBLE_COLOR, 
			      colorSize);
            }
            /* If not visible, if any is mapped, the combination is mapped. */
            else if ((IRIT_GEN_CMP(&VisMap1[Index], &RNDR_VISMAP_MAPPED_COLOR, 
				   colorSize) == 0) ||
		     (IRIT_GEN_CMP(&VisMap2[Index], &RNDR_VISMAP_MAPPED_COLOR, 
				   colorSize) == 0)) {
	        IRIT_GEN_COPY(&VisMap1[Index], &RNDR_VISMAP_MAPPED_COLOR, 
			      colorSize);
            }

            /* If not visible and not mapped, if any is tangent, the        */
            /* combination is tangent.                                      */
            else if ((IRIT_GEN_CMP(&VisMap1[Index], &RNDR_VISMAP_TANGENT_COLOR,
				   colorSize) == 0) ||
		     (IRIT_GEN_CMP(&VisMap2[Index], &RNDR_VISMAP_TANGENT_COLOR,
				   colorSize) == 0)) {
	        IRIT_GEN_COPY(&VisMap1[Index], &RNDR_VISMAP_TANGENT_COLOR, 
			      colorSize);
            }

            /* If not visible, not mapped and not tangent, if any is poor   */
            /* angle, the combination is poor angle. Though tangent and     */
            /* poor angle are equal in important. It's arbitrary decision   */
            /* to put tangent before poor angle.                            */
            else if ((IRIT_GEN_CMP(&VisMap1[Index], &RNDR_VISMAP_POOR_AR_COLOR,
				   colorSize) == 0) ||
		     (IRIT_GEN_CMP(&VisMap2[Index], &RNDR_VISMAP_POOR_AR_COLOR,
				   colorSize) == 0)) {
	        IRIT_GEN_COPY(&VisMap1[Index], &RNDR_VISMAP_POOR_AR_COLOR, 
			      colorSize);
            }

            /* If no other then we take degenerated (which is what left). */
            else
                IRIT_GEN_COPY(&VisMap1[Index], &RNDR_VISMAP_DEGEN_COLOR, 
			      colorSize);
        }
    }
    return VisMap1;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create visibility map to the given scene from the given observation      *
* point. If the aperture is above 120 degrees in x or y, divide the creation *
* of the visibility maps to several creation and combine the result to one   *
* final visibility map.                                                      *
* Scene is freed during the process.                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:    The geometric covering problem to be solved.                 *
*   Scene:      The scene to scan. Freed during the process.                 *
*   OpGroup:    The observation point group which Op is one of its           *
*               suggestions.                                                 *
*   Op:         The observation point from which the scene will be seen.     *
*   PicIndex:   The index of the picture created. For debug purposes only.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtImgPixelStruct*: The visibility map of the given scene from the given *
*                       observation point.                                   *
*****************************************************************************/
static IrtImgPixelStruct *UserGCCreateComplexVisMap(
				       UserGCProblemDefinitionStruct *Problem,
				       IPObjectStruct *Scene, 
                                       UserGCObsPtGroupTypeStruct *OpGroup,
				       UserGCObsPtSuggestionStruct *Op,
				       int PicIndex)
{
    IrtImgPixelStruct 
        *VisMap = NULL;
    char PicSuffix[20];
    int PicsNum, i, PicSuffixLength;
    IrtRType OpeningXY, OpeningZ;
    IrtHmgnMatType PrspMat, ViewMats[9];

    sprintf(PicSuffix, "%d", PicIndex);
    PicSuffixLength = (int)strlen(PicSuffix);

    /* Orthographic mode. */
    if(IRIT_PT_APX_EQ(Op->ObsPt, USER_GC_INF_VEC)) {
        IrtHmgnMatType ViewMat;

        UserGCCreateViewMatrix(Op, ViewMat);
        return UserGCCreateVisMap(Problem, Scene, ViewMat, NULL, 0, 0, 
				  PicSuffix);
    }

    /* Perspective mode. */
    /* Divide the space to 120 degrees angles with 5 degrees of overlap.    */
    UserGCDivideAndCreateViewMatrices(Op, OpGroup -> ObsPtType.XYAngle, 
				      OpGroup -> ObsPtType.ZAngle,
				      &PicsNum, &OpeningXY, &OpeningZ, 
				      ViewMats);
    UserGCCreatePrspMatrix(OpeningZ, OpeningXY, PrspMat);

    for (i = 0; i < PicsNum; i++) {
        IrtImgPixelStruct *Res; 
        IPObjectStruct *SceneCopy; 

        if (Problem -> DebugParams.LoadObjectsFromDisk < 3)
            SceneCopy = IPCopyObjectList(Scene, FALSE); 
        else if (Scene != NULL) {
            SceneCopy = Scene;
            Scene = NULL;
        }
        else 
            SceneCopy = UserGCLoadProcessedObjects2(Problem);

        sprintf(PicSuffix + PicSuffixLength, "_%d", i);
        Res = UserGCCreateVisMap(Problem, 
            SceneCopy, ViewMats[i], PrspMat,
                OpeningXY, OpeningZ, PicSuffix);

#if 0 /* Keep temp visibility maps for each subimage of the visibility map. */
        {
            char name[100];
            sprintf(name, "Temp%s", PicSuffix);
            Problem -> DebugParams.SaveImageFunc(name, 
                Res, Problem -> SolvingParams.VisMapWidth,
                Problem -> SolvingParams.VisMapHeight);
        }
#endif
        if (VisMap == NULL) {
            VisMap = Res;
            Res = NULL;
        }         
        UserGCCombineVisMaps(VisMap, Res, 
                             Problem -> SolvingParams.VisMapWidth,
                             Problem -> SolvingParams.VisMapHeight);
        IritFree(Res);        
    }
    IPFreeObjectList(Scene);

    /* Remove ignored UV pixels given by the user. */
    if (Problem -> UVMap != NULL) {
        for (i = 0; i < Problem -> SolvingParams.VisMapHeight * 
            Problem -> SolvingParams.VisMapWidth; i++) {
            if ((Problem -> UVMap[i].r == 0) &&
                (Problem -> UVMap[i].b == 0) &&
                (Problem -> UVMap[i].b == 0)) {
                IRIT_GEN_COPY(&VisMap[i], &RNDR_VISMAP_DEGEN_COLOR, 
                    sizeof(RNDR_VISMAP_EMPTY_COLOR));
            }
        }
    }

    return VisMap;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Create visibility map of the scene for each observation point suggestion M
* and save it to disk. The scene is given in GeoObj and Obstacles in Problem.M
* In case of an error longjmp to UserGCStartOfProcess.                       M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:    The geometric covering problem to be solved (Contains the    M
*               required observation points, the geometic object, the        M
*               obstacles and other parameters).                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCCreateVisMaps							     M
*****************************************************************************/
void UserGCCreateVisMaps(UserGCProblemDefinitionStruct *Problem)
{
    /* If we don't load visibility map we need to create them. */
    if (!Problem -> DebugParams.LoadVisMaps) {
        IPObjectStruct *OriginScene;
        UserGCObsPtGroupTypeStruct **OpGroup;
        IPObjectStruct *Scene; 
        int Index = -1;

        /* Combining scene to be used in creating the visibility maps and   */
        /* deattach it from Problem.                                        */
        Problem -> GeoObj -> Pnext = Problem -> Obstacles;
        OriginScene = Problem -> GeoObj;
        if (Problem -> DebugParams.LoadObjectsFromDisk > 1) {
            Problem -> GeoObj = NULL;
            Problem -> Obstacles = NULL;
        }

        /* Creating all visibility maps, each visibility map is created,    */
        /* saved to disk and then freed. This is in order to reduce memory  */
        /* use during the creation of the rest of the visibility maps.      */
        for (OpGroup = Problem -> ObsPtsGroups; *OpGroup != NULL; OpGroup++) {
            UserGCObsPtSuggestionStruct **Op;
            IrtImgPixelStruct *VisMap;

            for (Op = (*OpGroup) -> Suggestions; *Op != NULL; Op++) {
                Index++;
                Problem -> DebugParams.Print("%s Creating image %d.\n", 
		                         MiscISCGetTimeStamp(time(NULL), NULL, 
                                         USER_GC_TIME_ZONE, FALSE), Index);
                if (Problem -> DebugParams.LoadObjectsFromDisk < 2)
                    Scene = IPCopyObjectList(OriginScene, FALSE); 
                else if (OriginScene != NULL) {
                    Scene = OriginScene;
                    OriginScene = NULL;
                }
                else if (Problem -> DebugParams.LoadObjectsFromDisk == 2) 
                    Scene = UserGCLoadProcessedObjects2(Problem);
                else 
                    Scene = NULL;
                VisMap = UserGCCreateComplexVisMap(Problem, Scene, *OpGroup,
						   *Op, Index);

                /* In case of longjmp we will have memory leak of VisMap. We */
                /* can live with that.                                       */
                UserGCSaveVisMap(Index, VisMap, 
				 Problem -> SolvingParams.VisMapWidth, 
				 Problem -> SolvingParams.VisMapHeight,
				 Problem -> DebugParams.SaveImageFunc);

                /* Store the view matrix if asked by user. */
                if (Problem -> DebugParams.StoreVisMapMatrix) {
                    IrtHmgnMatType ViewMat;

                    UserGCCreateViewMatrix(*Op, ViewMat);
                    UserGCSaveMatrix(Index, ViewMat);
                }
                IritFree(VisMap);
            }
        }
        if (Problem -> DebugParams.LoadObjectsFromDisk < 2)
            Problem -> GeoObj -> Pnext = NULL;
    }
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Expose UserGCPrepareScene.                                               M
*   Preparing the scene to be rendered.                                      M
*   The objects stored in Problem->GeoObj and Problem->Obstacles are         M
* destroyed and mustn't be accessed again. They are replaced by the results  M
* of this function.                                                          M
*   Both objects shoul should already contain only object of type polygon    M
* connected only by pnext (not by U.Lst).                                    M
*   Both objects are going through the processing mentioned in               M
* UserGCPrepareObj. If Problem -> DebugParams.LoadObjectsFromDisk is TRUE    M
* the new created objects are saved to disk using UserGCSaveProcessedObject  M
* and can be later accessed using UserGCLoadProcessedObject if.              M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem: The geometric covering problem to be solved.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   int: return FALSE if failed saving the object (when required).           M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposePrepareScene                                                 M
*****************************************************************************/
int UserGCExposePrepareScene(UserGCProblemDefinitionStruct *Problem) 
{
    if (setjmp(UserGCStartOfProcess)) {
        return FALSE;
    }
    UserGCPrepareScene(Problem);
    return TRUE;
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Preparing the scene to be rendered.                                      M
*   The objects stored in Problem->GeoObj and Problem->Obstacles are         M
* destroyed and mustn't be accessed again. They are replaced by the results  M
* of this function.                                                          M
*   Both objects shoul should already contain only object of type polygon    M
* connected only by pnext (not by U.Lst).                                    M
*   Both objects are going through the processing mentioned in               M
* UserGCPrepareObj. If Problem -> DebugParams.LoadObjectsFromDisk is TRUE    M
* the new created objects are saved to disk using UserGCSaveProcessedObject  M
* and can be later accessed using UserGCLoadProcessedObject if.              M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem: The geometric covering problem to be solved.                    M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCPrepareScene                                                       M
*****************************************************************************/
void UserGCPrepareScene(UserGCProblemDefinitionStruct *Problem)
{
    IPPolygonStruct *Pl;
    int i;

    Problem -> GeoObj = UserGCPrepareObj(Problem -> GeoObj, 
					 Problem -> SolvingParams.VisMapWidth, 
					 Problem -> SolvingParams.VisMapHeight,
					 Problem -> DebugParams.GeoObjOrig);

    if (Problem -> Obstacles != NULL)
        Problem -> Obstacles = UserGCPrepareObj(Problem -> Obstacles, 0, 0, 
						NULL);

    if (Problem -> DebugParams.LoadObjectsFromDisk > 0)
        UserGCSaveProcessedObjects(Problem);

    for (i = 0, Pl = Problem -> GeoObj -> U.Pl; Pl != NULL; 
        Pl = Pl -> Pnext, i++);
    Problem -> DebugParams.Print("The object contains %d polygons\n", i);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Free the suggestions members of the given UserGCObsPtGroupTypeStruct and *
* reallocate the given number of suggestions pointers.                       *
*                                                                            *
* PARAMETERS:                                                                *
*   ObsPtGroup:           The suggestion of this obsevation group will be    *
*                         freed.                                             *
*   NewSuggestionsNumber: The number of new suggestions to allocate.         *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCReallocSuggestions                                                 *
*****************************************************************************/
static void UserGCReallocSuggestions(UserGCObsPtGroupTypeStruct *ObsPtGroup,
                                     int NewSuggestionsNumber)
{
    UserGCObsPtSuggestionStruct **Suggestion;

    if (ObsPtGroup -> Suggestions != NULL) {
        for (Suggestion = ObsPtGroup -> Suggestions;
	     *Suggestion != NULL; 
	     Suggestion++) {
	    IritFree(*Suggestion);
        }
        IritFree(ObsPtGroup -> Suggestions);
    }
    ObsPtGroup -> Suggestions = (UserGCObsPtSuggestionStruct **)
        IritMalloc((NewSuggestionsNumber + 1) * 
                                      sizeof(UserGCObsPtSuggestionStruct *));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Empty the suggestions members of the given UserGCObsPtGroupTypeStruct    *
* and fill it with Suggestions from the Array which is array of vectors      *
* (directions from infinity).                                                *
*                                                                            *
* PARAMETERS:                                                                *
*   ObsPtGroup:  Will be emptied and filled with suggestions.                *
*   Array:       Array containing suggestions as directions from infinity.   *
*   Number:      Number of vectors in Array.                                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCSetSuggestionPointsFromArray                                       *
*****************************************************************************/
static void UserGCSetSuggestionPointsFromArray(
                                     UserGCObsPtGroupTypeStruct *ObsPtGroup,
				     IrtVecType *Array,
				     int Number)
{
    int i;

    UserGCReallocSuggestions(ObsPtGroup, Number + 1);

    for (i = 0; i <= Number - 1; i++) {
        UserGCObsPtSuggestionStruct
	    *Suggestion = (UserGCObsPtSuggestionStruct *)
	                      IritMalloc(sizeof(UserGCObsPtSuggestionStruct));

        IRIT_PT_COPY(Suggestion -> ObsPt, USER_GC_INF_VEC);
        IRIT_PT_SCALE2(Suggestion -> Direction, Array[i], -1);
        ObsPtGroup -> Suggestions[i] = Suggestion;
    }
    ObsPtGroup -> Suggestions[i] = NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create a polygon from the given vertices.                                *
*                                                                            *
* PARAMETERS:                                                                *
*   AddTo        : Add the new created polygon to the Pnext of this polygon. *
*                  (Assumes Pnext is NULL).                                  *
*   PV1, PV2, PV3: vertices of the polygon.                                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*                                                                            *
* SEE ALSO:                                                                  *
*   GCSolve                                                                  *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCCreateTriangle                                                     *
*****************************************************************************/
static IPPolygonStruct *UserGCCreateTriangle(IPPolygonStruct *AddTo,
                                             IrtPtType PV1,
                                             IrtPtType PV2,
                                             IrtPtType PV3)
{
    IPVertexStruct
        *PV = IPAllocVertex2(NULL);

    IRIT_PT_COPY(PV -> Coord, PV3);
    PV = IPAllocVertex2(PV);
    IRIT_PT_COPY(PV -> Coord, PV2);
    PV = IPAllocVertex2(PV);
    IRIT_PT_COPY(PV -> Coord, PV1);
    AddTo -> Pnext = IPAllocPolygon(0, PV, NULL);
    return AddTo -> Pnext;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Empty the suggestions members of the given UserGCObsPtGroupTypeStruct    *
* and fill it with Suggestions from its polyhedron.                          *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:    The geometric covering problem to be solved.                 *
*   ObsPtGroup: Will be emptied and filled with suggestions.                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   None                                                                     *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCSetSuggestionPointsFromPolyhedron                                  *
*****************************************************************************/
static void UserGCSetSuggestionPointsFromPolyhedron(
                                      UserGCProblemDefinitionStruct *Problem,
                                      UserGCObsPtGroupTypeStruct *ObsPtGroup)
{
    IPObjectStruct 
        *PObj = ObsPtGroup -> PredefinedSuggestions.Polyhedron.Poly;
    IPPolygonStruct *Old, *New;
    IPPolyVrtxIdxStruct *Vertices;
    IPVertexStruct **CurVertex;
    int i, Number;

    /* Start dividing each triangle to 4 triangles. It is assumed that the  */
    /* model contains only triangles.                                       */
    Old = PObj -> U.Pl;
    PObj -> U.Pl = NULL;
    New = IPAllocPolygon(0, NULL, NULL); /* Dummy to hold the list's start. */
    for (i=0; 
	 i<= ObsPtGroup -> PredefinedSuggestions.Polyhedron.Level - 1;
	 i++) {
        IPPolygonStruct *OldCur, *NewCur;

        NewCur = New;
        OldCur = Old;
        while (OldCur != NULL) {
            IrtPtType P2, P4, P6;
            IrtRType *P1, *P3, *P5;

            P1 = OldCur -> PVertex -> Coord;
            P3 = OldCur -> PVertex -> Pnext -> Coord;
            P5 = OldCur -> PVertex -> Pnext -> Pnext -> Coord;
            IRIT_PT_BLEND(P2, P1, P3, 0.5);
            IRIT_PT_BLEND(P4, P3, P5, 0.5);
            IRIT_PT_BLEND(P6, P5, P1, 0.5);
            IRIT_PT_NORMALIZE(P2);
            IRIT_PT_NORMALIZE(P4);
            IRIT_PT_NORMALIZE(P6);
            NewCur = UserGCCreateTriangle(NewCur, P1, P2, P6);
            NewCur = UserGCCreateTriangle(NewCur, P2, P3, P4);
            NewCur = UserGCCreateTriangle(NewCur, P4, P5, P6);
            NewCur = UserGCCreateTriangle(NewCur, P2, P4, P6);
            OldCur = OldCur -> Pnext;
        }
        IPFreePolygonList(Old);
        Old = New -> Pnext;
    }
    IPFreePolygon(New);

    PObj -> U.Pl = Old;
    Vertices =  IPCnvPolyToPolyVrtxIdxStruct(PObj, FALSE, 0);
    PObj -> U.Pl = NULL;
    IPFreeObject(PObj);

    for (CurVertex = Vertices -> Vertices, Number = 0;
	 *CurVertex != NULL;
	 CurVertex++, Number++);

    UserGCReallocSuggestions(ObsPtGroup, Number + 1);

    for (CurVertex = Vertices -> Vertices, i = 0; *CurVertex != NULL;
        CurVertex++, i++) {
        ObsPtGroup -> Suggestions[i] = (UserGCObsPtSuggestionStruct *)
            IritMalloc(sizeof(UserGCObsPtSuggestionStruct));
        IRIT_PT_COPY(ObsPtGroup -> Suggestions[i] -> ObsPt, USER_GC_INF_VEC);

        IRIT_PT_COPY(ObsPtGroup -> Suggestions[i] -> 
		     Direction, (**CurVertex).Coord);
    }
    ObsPtGroup -> Suggestions[i] = NULL;
    Problem -> DebugParams.Print("%d observation points were created for the SPHERE_POLYHEDRON\n", 
	                         i);
    IPPolyVrtxIdxFree(Vertices);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   This function accept IPObjectStruct. If the object contains lists        *
* using U.Lst, the function transform them into a list using Pnext           *
* (recursively) and returns the results. If the object is of any other type  *
* the function just returns the object as it is. If used in IPForEachObj2 it *
* will flatten an entire IPObjectStruct which uses U.Lst into a list using   *
* Pnext.                                                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:  The PObj to flatten.                                              *
*   Param: Ignored.                                                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct: The flat object.                                         *
*                                                                            *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCFlatObj                                                            *
*****************************************************************************/
static IPObjectStruct *UserGCFlatObj(IPObjectStruct *PObj, void *Param)
{
    IPObjectStruct **PPObj, *Res, 
        *FinalRes = NULL, 
        *LastObject = NULL;

    if (!IP_IS_OLST_OBJ(PObj))
        return PObj;

    for (PPObj = PObj -> U.Lst.PObjList; PPObj != NULL; PPObj++) {
        /* Copy the attributes of PObj to its children. */
        IPAttributeStruct *Attr1, *Attr2;

        Attr1 = AttrCopyAttributes(PObj -> Attr);
        for (Attr2 = Attr1; Attr2 -> Pnext != NULL; Attr2 = Attr2 -> Pnext);
        Attr2 -> Pnext = (*PPObj) -> Attr;
        (*PPObj) -> Attr = Attr1;

        Res = IPForEachObj2(*PPObj, UserGCFlatObj, NULL);
        if (FinalRes == NULL) {
            FinalRes = Res;
            LastObject = IPGetLastObj(Res);
        }
        else {
            LastObject -> Pnext = Res;
            LastObject = IPGetLastObj(LastObject);
        }
    }
    PObj -> U.Lst.PObjList[0] = NULL;
    IPFreeObject(PObj);
    return FinalRes;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Load a surface, a curve and optionally obstacles from the given file.    *
* The file must contain only a terrain surface and a 2D curve and obstacles  *
* if Obstacles isn't NULL. Every element must appear in seperate object.     *
* The terrain surface must be named as TERRAIN_SRF_NAME.                     *
* the 2D curve is defined over the parametric domain of the surface, creating*
* a 3D curve in the euclidean space.                                         *
* The obstacles are optional if Obstacles isn't NULL if Obstacles is NULL no *
* obstacles are allowed. The obstacles are either closed surfaces of closed  *
* polygon mesh.                                                              * 
* The file may also contain one or more matrices (only the last one is       *
* considered). The matrix isn't applied over the 2D domain as it is in the   *
* parametric domain and not the euclidean domain.                            *
*                                                                            *
*   throws GCMessageException if the content of the file doesn't fit the     *
* above describtion.                                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   FileName:  IN, The name of the file.                                     *
*   Crv:       OUT, The returned curve.                                      *
*   Srf:       OUT, The returned surface.                                    *
*   Obstacles: OUT, If not NULL, expecting to receive obstacles as well.     *
*              They must be closed surfaces or closed polygon mesh.          *
*                                                                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: return FALSE in case of an error after writing warning message.     *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCLoadSrfAndCrv                                                      *
*****************************************************************************/
static int UserGCLoadSrfAndCrv(const char *FileName,
                               CagdCrvStruct **Crv,
                               CagdSrfStruct **Srf,
                               IPObjectStruct **Obstacles) 
{

    IPObjectStruct *PObj, *Temp,
        *TempObstacles = NULL;
    CagdCrvStruct 
        *TempCrv = NULL;
    CagdSrfStruct 
        *TempSrf = NULL;
    IrtHmgnMatType Mat;
    CagdBBoxStruct BBox;
    IrtRType UMin, UMax, VMin, VMax;

    PObj = IPGetDataFiles(&FileName, 1, FALSE, FALSE);
    if (PObj == NULL) {
        IRIT_WARNING_MSG_PRINTF(USER_GC_PERIMETER_CRV_FILE_FAIL, FileName, 
                                USER_GC_PERIMETER_CRV_EMPTY_FILE);
        return FALSE;
    }

    /* Change U.lst to Pnext list. */
    PObj = IPForEachObj2(PObj, UserGCFlatObj, NULL);

    /* Transform any instance object to a full object (instance object just  */
    /* points to another object as a copy of the other object).              */
    PObj = IPResolveInstances(PObj);

    /* Check there is exactly one curve and one surface and get the last     */
    /* matrix object and use it on the curve and surface.                    */
    MatGenUnitMat(Mat);
    for (Temp = PObj; Temp != NULL; Temp = Temp -> Pnext) {
        if (IP_IS_MAT_OBJ(Temp) &&
	    (Temp -> ObjName != NULL) && 
            (strcmp(Temp -> ObjName, USER_GC_VIEW_MATRIX_NAME) == 0)) {
            IRIT_GEN_COPY(Mat, *(Temp -> U.Mat), sizeof(IrtHmgnMatType));
        }
        else if (IP_IS_CRV_OBJ(Temp) &&
		 (Temp -> U.Crvs != NULL) && 
		 (Temp -> U.Crvs -> Pnext == NULL) &&
		 (TempCrv == NULL)) {
            TempCrv = Temp -> U.Crvs;  
        }
        else if (IP_IS_SRF_OBJ(Temp) && 
            (stricmp( IP_GET_OBJ_NAME(Temp), USER_GC_TERRAIN_SRF_NAME) == 0) &&
            (Temp -> U.Srfs != NULL) && 
            (Temp -> U.Srfs -> Pnext == NULL)  &&
            (TempSrf == NULL)){
            TempSrf = Temp -> U.Srfs;  
        }
        else if ((Obstacles!= NULL) &&
		 ((IP_IS_SRF_OBJ(Temp) && 
		   (Temp -> U.Srfs != NULL) &&
		   (Temp -> U.Srfs -> Pnext == NULL)) ||
		  (IP_IS_POLY_OBJ(Temp) &&
		   IP_IS_POLYGON_OBJ(Temp) && 
		   (Temp -> U.Pl != NULL))))  {
            if (TempObstacles == NULL)
                TempObstacles = IPCopyObject(NULL, Temp, FALSE);
            else
                IPGetLastObj(TempObstacles) -> Pnext =
		                               IPCopyObject(NULL, Temp, FALSE);
        }
        else {
            IPFreeObjectList(PObj);
            IRIT_WARNING_MSG_PRINTF(USER_GC_PERIMETER_CRV_FILE_FAIL, FileName,
				    USER_GC_PERIMETER_CRV_BAD_CONTENT);
            return FALSE;
        }
    }    
    {
        /* Notice the 2D curve doesn't go through matrix transformations.    */
        /* However, there is some side effect in this function, probably the */
        /* change to bspline which is required. For now I just use this work */
        /* around.                                                           */
        IrtHmgnMatType Mat2;

        MatGenUnitMat(Mat2);
        TempCrv = CagdCrvMatTransform(TempCrv, Mat2);
    }
    TempSrf = CagdSrfMatTransform(TempSrf, Mat);
    if ((Obstacles != NULL) && (TempObstacles != NULL)) {
        IPObjectStruct
	    *Temp = TempObstacles;

        TempObstacles = GMTransformObjectList(TempObstacles, Mat);
        IPFreeObjectList(Temp);
    }
    IPFreeObjectList(PObj);

    /* Making sure the curve is 2D. */
    CagdCrvBBox(TempCrv, &BBox);
    if ((CAGD_NUM_OF_PT_COORD(TempCrv -> PType) > 3) ||
        ((CAGD_NUM_OF_PT_COORD(TempCrv -> PType) == 3) &&
        ((BBox.Min[2] != 0) || (BBox.Max[2] != 0)))) {
        CagdCrvFree(TempCrv);
        CagdSrfFree(TempSrf);
        IRIT_WARNING_MSG_PRINTF(USER_GC_PERIMETER_CRV_FILE_FAIL, FileName,
				USER_GC_PERIMETER_CRV_BAD_DIM);
        return FALSE;
    }

    /* Making sure the curve is in the surface parametric domain. */
    CagdSrfDomain(TempSrf, &UMin, &UMax, &VMin, &VMax);
    if ((BBox.Min[0] < UMin) || (BBox.Max[0] > UMax) ||
        (BBox.Min[1] < VMin) || (BBox.Max[1] > VMax)) {
        CagdCrvFree(TempCrv);
        CagdSrfFree(TempSrf);
        IPFreeObjectList(TempObstacles);
        IRIT_WARNING_MSG_PRINTF(USER_GC_PERIMETER_CRV_FILE_FAIL, FileName,
				USER_GC_PERIMETER_CRV_NO_IN_DOMAIN);
        return FALSE;
    }
    *Crv = TempCrv;
    *Srf = TempSrf;
    if (Obstacles != NULL)
        *Obstacles = TempObstacles;
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Empty the suggestions members of the given UserGCObsPtGroupTypeStruct and*
* fill it with Suggestions defined by the perimeter suggestion.              *
*   The curve of the perimeter is expected to be defined by a surface which  *
* its normals points inside (or down) and a trimming curve which its tangents*
* are pointing clockwise. The file                                           *
* ObsPtGroup -> PredefinedSuggestions.Perimeter.CrfFileName must answer the  *
* requirement written in UserGCLoadSrfAndCrv when Obstacles is NULL).        *
*                                                                            *
* PARAMETERS:                                                                *
*   ObsPtGroup:  Will be emptied and filled with suggestions.                *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: return FALSE if failed reading the curve file after writing warning *
*        message.                                                            *
*                                                                            *
* SEE ALSO:                                                                  *
*   GCSolve                                                                  *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCSetSuggestionPointsFromPerimeter                                   *
*****************************************************************************/
static int UserGCSetSuggestionPointsFromPerimeter(
                                       UserGCObsPtGroupTypeStruct *ObsPtGroup)
{
    CagdCrvStruct *Crv2D, *Crv, *TempCrv;
    CagdSrfStruct *Srf, *Obstacles;
    IrtRType TMin, TMax, *Lengths, TStep, *CurPt, PrevSample, SampleStep;
    int i, j, LengthsNum;
    IrtPtType PreviousPt;

    Obstacles = NULL;
    /* Loading the curve and the surface. */
    if (!UserGCLoadSrfAndCrv(
            ObsPtGroup -> PredefinedSuggestions.Perimeter.CrfFileName, &Crv2D, 
	    &Srf, NULL))
        return FALSE;

    /* Different functions transform periodic into floating and then perhaps*/
    /* to open end. All those changes cause unwanted changes in the         */
    /* parameterization. I solve it by initially do the transformation to   */
    /* open end by myself.                                                  */
    TempCrv = CagdCnvrtBsp2OpenCrv(Crv2D);
    CagdCrvFree(Crv2D);
    Crv2D = TempCrv;
    Crv = SymbComposeSrfCrv(Srf, Crv2D);

    /* Creating observation points on the curve. */

    /* Creating an array with LengthsNum samples (uniformly selected over t)*/
    /* on the curve. Instead of calculating the exact distance on the curve */
    /* between each sample, I calculate the direct distance between the two */
    /* samples. In each cell of the array I put the distance from start of  */
    /* that point (calculated as described above). It's equilant to turning */
    /* the curve into a polygon and calculate the distances on the polygon  */
    /* but now I also have those distances according to given t value.      */

    /* I use 10 samples for each required guard and not less than 10000.    */
    LengthsNum = IRIT_MAX(10000, 
             10 * ObsPtGroup -> PredefinedSuggestions.Perimeter.GuardsNumber);
    Lengths = (IrtRType *) IritMalloc((LengthsNum + 1) * sizeof(IrtRType));
    CagdCrvDomain(Crv, &TMin, &TMax);
    TStep = (TMax - TMin)/LengthsNum;
    Lengths[0] = 0;
    CurPt = CagdCrvEval(Crv, 0);
    if (CAGD_IS_RATIONAL_CRV(Crv))
        IRIT_PT_SCALE(CurPt+1, 1 / CurPt[0]);
    IRIT_PT_COPY(PreviousPt, CurPt + 1);

    for (i = 1; i < LengthsNum + 1; i++) {
        CurPt = CagdCrvEval(Crv, i * TStep);
        if (CAGD_IS_RATIONAL_CRV(Crv))
            IRIT_PT_SCALE(CurPt + 1, 1 / CurPt[0]);
        Lengths[i] = Lengths[i - 1] + IRIT_PT_PT_DIST(CurPt + 1, PreviousPt);
        IRIT_PT_COPY(PreviousPt, CurPt + 1);
    }

    /* Sampling in LengthsNum uniform locations over the created array. */
    UserGCReallocSuggestions(ObsPtGroup, LengthsNum + 1);
    SampleStep = Lengths[LengthsNum]/
        ObsPtGroup -> PredefinedSuggestions.Perimeter.GuardsNumber;
    PrevSample = -SampleStep; /* In order to make first sample at Length[0].*/
    j = 0;
    for (i = 0;
	 (i < LengthsNum + 1) && 
	   (PrevSample + SampleStep < Lengths[LengthsNum] - SampleStep / 2);
	 i++) {
        if (Lengths[i] >= PrevSample + SampleStep) {
            IrtRType t = i * TStep;
            IrtRType *Pt;
            IrtPtType Temp;
            UserGCObsPtSuggestionStruct
	        *Suggestion = (UserGCObsPtSuggestionStruct *)
	                     IritMalloc(sizeof(UserGCObsPtSuggestionStruct));
            CagdVecStruct *Normal, *Tangent;

            Pt = CagdCrvEval(Crv, t);
            IRIT_PT_COPY(Temp, Pt + 1);
            if (CAGD_IS_RATIONAL_CRV(Crv))
                IRIT_PT_SCALE(Temp, 1 / Pt[0]);

            /* Location of observation point. */
            IRIT_PT_COPY(Suggestion -> ObsPt, Temp);
            Suggestion -> ObsPt[1] += 
                ObsPtGroup -> PredefinedSuggestions.Perimeter.GuardsHeight;

            /* Direction of observation point. */
            Pt = CagdCrvEval(Crv2D, t);
            if (CAGD_IS_RATIONAL_CRV(Crv))
                IRIT_PT_SCALE(&Pt[1], 1 / Pt[0]);

            Normal = CagdSrfNormal(Srf, Pt[1], Pt[2], TRUE);
            Tangent = CagdCrvTangent(Crv, t, TRUE);
            IRIT_CROSS_PROD(Suggestion -> Direction, Tangent -> Vec,
			    Normal -> Vec);
            PrevSample = PrevSample + SampleStep;
            ObsPtGroup -> Suggestions[j] = Suggestion;
            j++;
        }
    }
    ObsPtGroup -> Suggestions[j] = NULL;
    CagdCrvFree(Crv);
    CagdCrvFree(Crv2D);
    CagdSrfFree(Srf);
    IritFree(Lengths);
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Set the the color in Param to the triangle Pl.                           *
*   Can be used with IPForEachPoly2 to set a certain color to all polygons.  *
*                                                                            *
* PARAMETERS:                                                                *
*   Pl:     The polygon to set its color.                                    *
*   Param:  IrtImgPixelStruct object with the color.                         *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:  The input polygon with the color attribute.           *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCSetColor                                                           *
*****************************************************************************/
static IPPolygonStruct *UserGCSetColor(IPPolygonStruct *Pl, void *Param)
{
    IrtImgPixelStruct
        *Pixel = (IrtImgPixelStruct *) Param;

    AttrSetRGBColor(&Pl -> Attr, Pixel -> r, Pixel -> g, Pixel -> b);
    return Pl;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Set the given color to all the polygon list.                             M
*                                                                            *
* PARAMETERS:                                                                M
*   Pl:    The polygons list to set their color.                             M
*   r,g,b: The color to set.                                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCSetColorPlList                                                     M
*****************************************************************************/
static void UserGCSetColorPlList(IPPolygonStruct *Pl, int r, int g, int b)
{
    IrtImgPixelStruct Pixel;

    Pixel.r = r;
    Pixel.g = g;
    Pixel.b = b;

    IPForEachPoly2(Pl, UserGCSetColor, &Pixel);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Accept an IPObjectStruct which is either surfaces, trimmed surfaces or   *
* polygons (Not polylines or points list. Polystrip is accepted). Any other  *
* type will cause an assert (Including IP_OBJ_LIST_OBJ). It will transform   *
* each of those objects to a polygon objects containing only triangles.      *
*   Objects connected with Pnext are ignored.                                *
*   The transformation is done in place.                                     *
*   Can be used with IPForEachObj2 to do it to a Pnext list.                 *
*   If PObj has color attribute, it's propgated to all the created polygons. *
*                                                                            *
*   The triangulation of surfaces and trimmed surfaces was copied from       *
* irender/parser.c and altered.                                              *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:  The object to be triangulated in place.                           *
*   Param: Ignored.                                                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:  The triangulated object.                              *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCTriangulateObjAux                                                  *
*****************************************************************************/
static IPObjectStruct *UserGCTriangulateObjAux(IPObjectStruct *PObj, 
                                               void *Param)
{
    IPPolygonStruct *Poly, 
        *FinalPoly = NULL;
    int GenTriOnly = IPSurface2PolygonsGenTriOnly(TRUE),
        GenDegenerated = IPSurface2PolygonsGenDegenPolys(TRUE);
    IrtRType RelativeFineNess;
    int r, g, b;

    /* Prevent unnecessarily returning this process. */ 
    if (AttrGetIntAttrib(PObj -> Attr, 
			 USER_GC_ALREADY_TRIANGULATED_ATTRIB) !=
	                                                    IP_ATTR_BAD_INT)
        return PObj;
    AttrSetIntAttrib(&PObj -> Attr, 
		     USER_GC_ALREADY_TRIANGULATED_ATTRIB, TRUE);

    RelativeFineNess = AttrGetObjectRealAttrib(PObj, USER_GC_RESOLUTION);
    if (IP_ATTR_IS_BAD_REAL(RelativeFineNess))
        RelativeFineNess = 1.0;
    
    if (IP_IS_SRF_OBJ(PObj)) {
        CagdSrfStruct *Surfaces, *Surface;

        Surfaces = PObj -> U.Srfs;
        for (Surface = Surfaces; Surface != NULL; Surface = Surface -> Pnext) {
            IrtRType t;

            t = AttrGetObjectRealAttrib(PObj, USER_GC_RESOLUTION);
            if (!IP_ATTR_IS_BAD_REAL(t))
                AttrSetRealAttrib(&Surface -> Attr, USER_GC_RESOLUTION, t);
            t = AttrGetObjectRealAttrib(PObj, USER_GC_V_RESOLUTION);
            if (!IP_ATTR_IS_BAD_REAL(t))
                AttrSetRealAttrib(&Surface -> Attr, USER_GC_RESOLUTION, t);

            Poly = IPSurface2Polygons(Surface,
                                        FALSE,
                                        RelativeFineNess *
                                            UserGCSrf2PlgFineness,
                                        TRUE,
                                        TRUE,
                                        UserGCSrf2PlgOptimal);
            if (FinalPoly == NULL)
                FinalPoly = Poly;
            else {
                IPGetLastPoly(Poly) -> Pnext = FinalPoly;
                FinalPoly = Poly;
            }
        }
        
        CagdSrfFreeList(Surfaces);
    }
    else if(IP_IS_TRIMSRF_OBJ(PObj)) {
        TrimSrfStruct *TrimSrfs, *TrimSrf;

        TrimSrfs = PObj -> U.TrimSrfs;
        for (TrimSrf = TrimSrfs; TrimSrf != NULL; TrimSrf = TrimSrf -> Pnext) {
            IrtRType t;

            t = AttrGetObjectRealAttrib(PObj, USER_GC_RESOLUTION);
            if (!IP_ATTR_IS_BAD_REAL(t)) {
                AttrSetRealAttrib(&TrimSrf -> Attr, USER_GC_U_RESOLUTION, t);
                AttrSetRealAttrib(&TrimSrf -> Srf -> Attr,
				  USER_GC_U_RESOLUTION, t);
            }
            t = AttrGetObjectRealAttrib(PObj, USER_GC_V_RESOLUTION);
            if (!IP_ATTR_IS_BAD_REAL(t)) {
                AttrSetRealAttrib(&TrimSrf -> Attr, USER_GC_V_RESOLUTION, t);
                AttrSetRealAttrib(&TrimSrf -> Srf -> Attr,
				  USER_GC_V_RESOLUTION, t);
            }

            Poly = IPTrimSrf2Polygons(TrimSrf,
                                        FALSE,
                                        RelativeFineNess *
                                            UserGCSrf2PlgFineness,
                                        TRUE,
                                        TRUE,
                                        UserGCSrf2PlgOptimal);
            if (FinalPoly == NULL)
                FinalPoly = Poly;
            else {
                IPGetLastPoly(Poly) -> Pnext = FinalPoly;
                FinalPoly = Poly;
            }
        }
        TrimSrfFreeList(TrimSrfs);
    }
    else if (IP_IS_POLY_OBJ(PObj) && 
             (IP_IS_POLYGON_OBJ(PObj) || IP_IS_POLYSTRIP_OBJ(PObj))) {
        int OldCirc;
        IPObjectStruct *PObj2;

        if (IP_IS_POLYSTRIP_OBJ(PObj)) {
            assert("Polystrip isn't yet supported." == 0);
        }

        /* Make sure all polygons are convex.                        */ 
        /* GMConvertPolysToTriangles already has call to             */
        /* GMConvexPolyObject. However, GMConvexPolyObject requires  */
        /* circular polygon list. So I need to surround              */
        /* GMConvertPolysToTriangles with IPOpenPolysToClosed and    */
        /* IPClosedPolysToOpen. However, I'm not sure about the      */
        /* influence of such act on other functions inside           */
        /* TriangulateInPlace. So my best action is to make all      */
        /* polygon marked as convex by GMConvexPolyObject here (in   */
        /* Tag field), and  that way GMConvexPolyObject inside       */
        /* GMConvertPolysToTriangles will recognize that they are    */
        /* convex and won't need the vertices to be in circular      */
        /* list. Notice, that in some cases if the original object   */
        /* was a free form, the polygons are already marked as       */
        /* convex in this point.                                     */
        OldCirc = IPSetPolyListCirc(TRUE);
        IPOpenPolysToClosed(PObj -> U.Pl);
        GMConvexPolyObject(PObj);
        IPSetPolyListCirc(OldCirc);
        IPClosedPolysToOpen(PObj -> U.Pl);

        /* Make sure all polygons are triangle. */
        PObj2 = GMConvertPolysToTriangles(PObj);
        IPFreePolygonList(PObj -> U.Pl);
        FinalPoly = PObj2 -> U.Pl;
        PObj2 -> U.Pl = NULL;
        IPFreeObject(PObj2);
    }
    IPSurface2PolygonsGenTriOnly(GenTriOnly);
    IPSurface2PolygonsGenDegenPolys(GenDegenerated);

    /* Propgate the color to all polygons. */
    if (AttrGetRGBColor(PObj -> Attr, &r, &g, &b))
        UserGCSetColorPlList(FinalPoly, r, g, b);
    PObj -> ObjType = IP_OBJ_POLY;
    IP_SET_POLYGON_OBJ(PObj);
    PObj -> U.Pl = FinalPoly;
    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Use IPForEachObj2 with GCTriangulateObjAux. See UserGCTriangulateObjAux  *
* for details.                                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:          The object to triangulate. This object can't be used      *
*                  anymore after this function returns. It mustn't be freed  *
*                  either. Use the returned object instead.                  *
*   SrfPlgfine:    The fineness in which the loaded object turned into       *
*                  polygons.                                                 *
*   SrfPlgOptimal: If true, the triangulation will be done using adaptive    *
*                  method.                                                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:  The triangulated object.                              *
*****************************************************************************/
static IPObjectStruct *UserGCTriangulateObj(IPObjectStruct *PObj,
                                            IrtRType SrfPlgfine,
                                            int SrfPlgOptimal)
{
    UserGCSrf2PlgFineness = SrfPlgfine;
    UserGCSrf2PlgOptimal = SrfPlgOptimal;
    return IPForEachObj2(PObj, UserGCTriangulateObjAux, NULL);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Empty the suggestions members of the given UserGCObsPtGroupTypeStruct    *
* and fill it with Suggestions defined by the surface grid suggestion.       *
*   The surface file expects to contain several objects. One is the surface  *
* which must be name "Terrain" (case insensitive). The other is a 2D curve   *
* over the domain of the surface - a trimming curve. Additional object will  *
* be treated as obstacles.                                                   *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:    The geometric covering problem to be solved.                 *
*   ObsPtGroup: Its suggestion list will be emptied and filled with          *
*               suggestions defined by the surface grid formal suggestion.   *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: return FALSE if failed reading the surface or no observation points *
*        were left. Also writing warning message.                            *
*                                                                            *
* KEYWORDS:                                                                  *
*   UserGCSetSuggestionPointsFromSurfaceGrid                                 *
*****************************************************************************/
static int UserGCSetSuggestionPointsFromSurfaceGrid(
                                      UserGCProblemDefinitionStruct *Problem,
                                      UserGCObsPtGroupTypeStruct *ObsPtGroup)
{
    IPPolygonStruct *PolyCrv;
    IPVertexStruct 
        *V = NULL;
    CagdCrvStruct *Crv, *TempCrv;
        CagdSrfStruct *Srf;
    IPObjectStruct *Obstacles;
    IrtRType TMin, TMax, TStep, *CurPt, UMin, UMax, VMin, VMax, UStep, VStep;
    int i, j, TotalGuards,
        UGuardsNum = 
            ObsPtGroup -> PredefinedSuggestions.SurfaceGrid.GuardsNumber[0],
        VGuardsNum = 
            ObsPtGroup -> PredefinedSuggestions.SurfaceGrid.GuardsNumber[1],
        LengthsNum = 1000;
    UserGCObsPtSuggestionStruct Suggestion;

    if (ObsPtGroup -> ObsPtType.XYAngle != 360) {
        IRIT_WARNING_MSG("SURFACE_GRID supports only XYAngle of 360 degrees.\n"
            "Therefore a 360 degrees will be used for the XYAngle.\n");
        ObsPtGroup -> ObsPtType.XYAngle = 360;
    }

    /* Loading the curve and the surface. */
    if (!UserGCLoadSrfAndCrv(
            ObsPtGroup -> PredefinedSuggestions.SurfaceGrid.SrfFileName,
	    &Crv, &Srf, &Obstacles))
        return FALSE;
    /* We assume without check that all surfaces are closed. */
    Obstacles = UserGCTriangulateObj(Obstacles, 0.0001, TRUE); 

    CagdCrvDomain(Crv, &TMin, &TMax);
    /* Different functions transform periodic into floating and then perhaps*/
    /* to open end. All those changes cause unwanted changes in the         */
    /* parameterization. I solve it by initially do the transformation to   */
    /* open end by myself.                                                  */
    TempCrv = CagdCnvrtBsp2OpenCrv(Crv);
    CagdCrvFree(Crv);
    Crv = TempCrv;

    /* Creating polygon curve out of the given curve. */ 
    TStep = (TMax - TMin)/LengthsNum;
    for (i = 0; i < LengthsNum; i++) {
        CurPt = CagdCrvEval(Crv, i * TStep);

        if (CAGD_IS_RATIONAL_CRV(Crv))
            IRIT_PT_SCALE(CurPt + 1, 1 / CurPt[0]);

        V = IPAllocVertex2(V);
        IRIT_PT_COPY(V -> Coord, CurPt + 1);    
    }
    PolyCrv = IPAllocPolygon(0, V, NULL);

    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
    UStep = (UMax - UMin)/(UGuardsNum - 1);
    VStep = (VMax - VMin)/(VGuardsNum - 1);

    /* Sampling the guard and taking those inside the surface domain. */
    UserGCReallocSuggestions(ObsPtGroup, VGuardsNum*UGuardsNum+1);
    TotalGuards = 0;
    for (i = 0; i < VGuardsNum; i++) {
        IrtRType
	    VValue = IRIT_MIN(VMin + i * VStep, VMax);
        IrtPtType Pt;

        for (j = 0; j < UGuardsNum; j++) {
            IrtRType 
                UValue = IRIT_MIN(UMin + j * UStep, UMax);
            IrtRType *Pt2;
            IPObjectStruct *PObj;

            IRIT_PT_SET(Pt, UValue, VValue, 0);

            Pt2 = CagdSrfEval(Srf, UValue, VValue);

            if (CAGD_IS_RATIONAL_SRF(Srf))
                IRIT_PT_SCALE(Pt2 + 1, 1 / Pt2[0]);

            Pt2[2] += 
                ObsPtGroup -> PredefinedSuggestions.SurfaceGrid.GuardsHeight;

            for (PObj = Obstacles; PObj != NULL; PObj = PObj -> Pnext)  
                if (GMPolygonPointInclusion3D(PObj -> U.Pl, Pt2 + 1)) 
                    break;

            if ((PObj == NULL) && GMPolygonPointInclusion(PolyCrv, Pt)) {
                IRIT_PT_COPY(Suggestion.ObsPt, Pt2 + 1);           
                ObsPtGroup -> Suggestions[TotalGuards] = 
                    (UserGCObsPtSuggestionStruct *)
                        IritMalloc(sizeof(UserGCObsPtSuggestionStruct));
                IRIT_GEN_COPY(ObsPtGroup -> Suggestions[TotalGuards], 
			      &Suggestion,
			      sizeof(UserGCObsPtSuggestionStruct));
                /* We only support 360 degrees in which case the direction */
                /* isn't important.                                        */
                IRIT_PT_SET(ObsPtGroup -> Suggestions[TotalGuards] -> 
			    Direction, 0, 0, -1);
                IRIT_PT_COPY(ObsPtGroup -> Suggestions[TotalGuards] ->
			     ObsPt, Pt2 + 1);
                TotalGuards++;
            }
        }
    }
    ObsPtGroup -> Suggestions[TotalGuards] = NULL;

    Problem -> DebugParams.Print(
                  "%d observation points were created for the SURFACE_GRID\n", 
	          TotalGuards);

    IPFreePolygon(PolyCrv);
    CagdCrvFree(Crv);
    CagdSrfFree(Srf);
    IPFreeObjectList(Obstacles);

    if (TotalGuards == 0) {
        IRIT_WARNING_MSG_PRINTF(USER_GC_NO_SOLVING_TITLE, USER_GC_NO_OPS);
        return FALSE;
    }

    return TRUE;
}

#ifdef USER_GC_EXPOSE_INNER_FUNCTIONALITIES

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Expose UserGCInterpretOPGroupsSuggestion.                                M
*   The observation points groups of this problem may contain implicitly     M
* defined observation points (such as 20 obsevation points around a sphere). M
*   This function will create explicitly defined observation points for all  M
* those observation points.                                                  M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:   The problem to process.				             M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if successful, FALSE otherwise.                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCExposeInterpretOPGroupsSuggestion                                  M
*****************************************************************************/
int UserGCExposeInterpretOPGroupsSuggestion(
                                       UserGCProblemDefinitionStruct *Problem)
{
    return UserGCInterpretOPGroupsSuggestion(Problem);
}

#endif /* USER_GC_EXPOSE_INNER_FUNCTIONALITIES */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   The observation points groups of this problem may contain implicitly     M
*   defined observation points (such as 20 obsevation points around a sphere)M
*   This function will create explicitly defined observation points for all  M
*   those observation points.                                                M
*                                                                            *
* PARAMETERS:                                                                M
*   Problem:   The problem to process.				             M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:   TRUE if successful, FALSE otherwise.                              M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCInterpretOPGroupsSuggestion                                        M
*****************************************************************************/
int UserGCInterpretOPGroupsSuggestion(UserGCProblemDefinitionStruct *Problem)
{
    UserGCObsPtGroupTypeStruct **ObsPtGroup;

    for (ObsPtGroup = Problem -> ObsPtsGroups;
	 *ObsPtGroup != NULL; 
	 ObsPtGroup++) {
        switch ((*ObsPtGroup) -> PredefinedSuggestions.PredefinedSuggestionType) {
            case USER_GC_CUSTOM_SUGGESTION:
                break;
            case USER_GC_SIX_MAIN_AXIS:
                UserGCSetSuggestionPointsFromArray(*ObsPtGroup,
						   UserGCSixMainAxisDirection,
						   6);
                break;
            case USER_GC_SPHERE_POLYHEDRON:
                UserGCSetSuggestionPointsFromPolyhedron(Problem, *ObsPtGroup);
                break;
            case USER_GC_SPHERE_COVER_4:
                UserGCSetSuggestionPointsFromArray(*ObsPtGroup,
						   GMSphereCoverVectors4,
						   4);
                break;
            case USER_GC_SPHERE_COVER_20:
                UserGCSetSuggestionPointsFromArray(*ObsPtGroup,
						   GMSphereCoverVectors20,
						   20);
                break;
            case USER_GC_SPHERE_COVER_50:
                UserGCSetSuggestionPointsFromArray(*ObsPtGroup,
						   GMSphereCoverVectors50,
						   50);
                break;
            case USER_GC_SPHERE_COVER_100:
                UserGCSetSuggestionPointsFromArray(*ObsPtGroup,
						   GMSphereCoverVectors100,
						   100);
                break;
            case USER_GC_SPHERE_COVER_130:
                UserGCSetSuggestionPointsFromArray(*ObsPtGroup,
						   GMSphereCoverVectors130,
						   130);
                break;
            case USER_GC_PERIMETER:
                if (!UserGCSetSuggestionPointsFromPerimeter(*ObsPtGroup))
                    return FALSE;
                break;
            case USER_GC_SURFACE_GRID:
                if (!UserGCSetSuggestionPointsFromSurfaceGrid(Problem, 
                                                              *ObsPtGroup))
                    return FALSE;
                break;                
            default:
                assert("Unknown observation point formal suggestion." == 0);
        }

        if ((*ObsPtGroup) -> PredefinedSuggestions.AddAntipodal) {
            int i, j, k, Len;
            UserGCObsPtSuggestionStruct **Suggestions;

            for (Suggestions = (*ObsPtGroup) -> Suggestions, Len = 0; 
		 *Suggestions != NULL;
		 Suggestions++, Len++);

            Suggestions = (UserGCObsPtSuggestionStruct **)
                          IritMalloc((2 * Len + 1) *
				     sizeof(UserGCObsPtSuggestionStruct *));
            IRIT_GEN_COPY(Suggestions, (*ObsPtGroup) -> Suggestions, 
			  Len * sizeof(UserGCObsPtSuggestionStruct*));
            IritFree((*ObsPtGroup) -> Suggestions);
            (*ObsPtGroup) -> Suggestions = Suggestions;
            
            i = Len;
            for (j = 0; j <= Len - 1; j++) {
                UserGCObsPtSuggestionStruct Suggestion;

                IRIT_PT_COPY(Suggestion.ObsPt, Suggestions[j] -> ObsPt);
                IRIT_PT_SCALE2(Suggestion.Direction,
			       Suggestions[j] -> Direction, -1);
                for (k = 0; k <= Len - 1; k++) {
                    if (IRIT_PT_APX_EQ(Suggestion.Direction, 
                            Suggestions[k] -> Direction) &&
                        IRIT_PT_APX_EQ(Suggestion.ObsPt,
                            Suggestions[k] -> ObsPt))
                        break;
                }
                if (k == Len) {
                    Suggestions[i] = (UserGCObsPtSuggestionStruct *)
                        IritMalloc(sizeof(UserGCObsPtSuggestionStruct));
                    IRIT_GEN_COPY(Suggestions[i], &Suggestion, 
				  sizeof(UserGCObsPtSuggestionStruct));
                    i++;
                }
            }
            Suggestions[i] = NULL;
        }
        (*ObsPtGroup) -> PredefinedSuggestions.PredefinedSuggestionType =
                                                    USER_GC_CUSTOM_SUGGESTION;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Auxilary function for UserGCPrepareObjAux. Used in IPForEachObj2 in       *
*  order to combine all polygons object to one polygon object.               *
*  If PObj isn't a polygon, free it and returns NULL.                        *
*  When using this function with IPForEachObj2 the returned IPObjectStruct   *
*  is empty IPObjectStruct.                                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:  The polygons to add to Param.                                     *
*   Param: A polygon IPObjectStruct to add PObj -> U.Pl to.                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *: Always return NULL (match pattern of call back func).  *
*****************************************************************************/
static IPObjectStruct *UserGCPrepareObjAux(IPObjectStruct *PObj, void *Param)
{
    IPPolygonStruct
        *Pl = IPGetLastPoly(((IPObjectStruct *) Param) -> U.Pl);

    if (Pl == NULL)
        ((IPObjectStruct *) Param) -> U.Pl = PObj -> U.Pl;
    else {
        Pl -> Pnext = PObj -> U.Pl;
    }

    PObj -> U.Pl = NULL;
    IPFreeObject(PObj);

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Prepares an object which means:                                          M
*     * Unites all polygon objects into one polygon object.                  M
*     * Change the uv values of all objects so they won't collide with each  M
*       other's uv spaces.                                                   M
*   The returned object may differ from PObj. PObj mustn't be accessed again M
*   after calling this function.                                             M
*                                                                            *
* PARAMETERS:                                                                M
*   PObj:                The object to prepare. The object should already    M
*                        contain only object of type polygon connected only  M
*                        by pnext (not by U.Lst). This object can't be used  M
*                        anymore after this function returns. Use the return M
*                        object instead.                                     M
*   MapWidth, MapHeight: The dimension of the visibility map. Used in order  M
*                        to arrange the UV domain of all PObj's objects. If  M
*                        any of them is 0, the function won't handle the UV  M
*                        domain of the objects.                              M
*   PObj2:               If it isn't NULL, this object is a list of objects  M
*                        which contains at least as much elements as PObj.   M
*                        Each element i in PObj2 will go through the same    M
*                        UV Transformations and scaling as element i in PObj.M
*                        Objects which aren't surfaces, trimmed surfaces or  M
*                        polygons are ignored.                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:   The prepared scene.                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCPrepareObj                                                         M
*****************************************************************************/
static IPObjectStruct *UserGCPrepareObj(IPObjectStruct *PObj, 
                                        int MapWidth, 
                                        int MapHeight,
                                        IPObjectStruct *PObj2)
{
    IPObjectStruct *PObjList;

    /* Arrange the UV values of all objects in PObj in to one non           */
    /* overlapping UV space. If PObj2 isn't NULL set its UV values with     */
    /* the same values as thouse of PObj.                                   */
    if ((MapWidth != 0) && (MapHeight != 0))
        IRndrVisMapPrepareUVValuesOfGeoObj(PObj, MapWidth, MapHeight, PObj2);

    /* Goes over all the objects in the list (using Pnext) and unites all   */
    /* polygon objects into one polygon object.                             */
    PObjList = IPGenPOLYObject(NULL);
    PObj = IPForEachObj2(PObj, UserGCPrepareObjAux, PObjList);
    IPFreeObjectList(PObj);
    return PObjList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Return number of observation points in the given problem.                M
*                                                                            *
* PARAMETERS:                                                                M
*    Problem:   The problem to count its observation points.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   int: The number of observation points in the problem.                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCGetOPsNum                                                          M
*****************************************************************************/
int UserGCGetOPsNum(UserGCProblemDefinitionStruct *Problem) 
{
    UserGCObsPtGroupTypeStruct **OpGroup;
    int Count = 0;

    for (OpGroup = Problem -> ObsPtsGroups; *OpGroup != NULL; OpGroup++) {
        UserGCObsPtSuggestionStruct **Op;

        for (Op = (*OpGroup) -> Suggestions; *Op != NULL; Op++) 
            Count++;
    }
    return Count;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Image type has no meaning since the image is written to a buffer.        *
*   Therefore this function just return an arbitrary image type.             *
*   It just retuns IRIT_IMAGE_PPM3_TYPE.                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   ImageType:  A string describing the image type.                          *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtImgImageType:  Returns IRIT_IMAGE_PPM3_TYPE as the detected type.     *
*****************************************************************************/
static IrtImgImageType UserGCImgWriteSetType(const char *ImageType)
{
    return IRIT_IMAGE_PPM3_TYPE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Prepare the variables required for keeping the image in a buffer.        *
*   The first three parameters aren't used. They are part of the required    *
*   signature for the use of this function.                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   argv:     Pointer to the name of this program. Not used.                 *
*   FName:    Filename to open. Not used.                                    *
*   Alpha:    Do we have aan alpha channel. Not used.                        *
*   XSize:    X dimension of the image.                                      *
*   YSize:    Y dimension of the image.                                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  TRUE if successful, FALSE otherwise. It never fails unless malloc  *
*         fails.                                                             *
*****************************************************************************/
static int UserGcImgWriteOpenFile(const char **argv,
                                  const char *FName,
                                  int Alpha,
                                  int XSize,
                                  int YSize)
{
    /* The semantic is that in order to use the created image (after it      */
    /* finished), one must either set UserGCVisMap to NULL or copy its       */
    /* content. Therefore, freeing it here, either frees a NULL pointer or   */
    /* frees a previously unused or copied image.                            */
    IritFree(UserGCVisMap);
    UserGCVisMap = (IrtImgPixelStruct *) IritMalloc(XSize * YSize * 
						    sizeof(IrtImgPixelStruct));
    if (!UserGCVisMap) {
        longjmp(UserGCStartOfProcess, 1);
    }

    UserGCVisMapSize[0] = XSize;
    UserGCVisMapSize[1] = YSize;
    UserGCVisMapLine = 0;
    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Store a line of pixels to the image buffer. Ignore the alpha channel.    *
*                                                                            *
* PARAMETERS:                                                                *
*   Alpha:  array of alpha values.                                           *
*   Pixels: array of color pixels.                                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void UserGCImgWritePutLine(IrtBType *Alpha, IrtImgPixelStruct *Pixels)
{
    int i,
        Start = UserGCVisMapLine * UserGCVisMapSize[0];
    IrtImgPixelStruct
        *VisMap = UserGCVisMap + Start;

    if (VisMap == NULL) {
        IRIT_WARNING_MSG_PRINTF("Error allocating memory in UserGCImgWritePutLine.\nThis message will be followed by exit of the program (this is the time to use debugger).");
	exit(-1);
    }

    for (i = 0; i <= UserGCVisMapSize[0] - 1; i++) {
        IRIT_GEN_COPY(&VisMap[i], &Pixels[i], sizeof(IrtImgPixelStruct));
    }
    UserGCVisMapLine++;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   No action is actually required to close the buffer image so this is just *
*   a dummy function.                                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   None                                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*                                                                            *
*****************************************************************************/
static void UserGCImgWriteCloseFile(void)
{
}

