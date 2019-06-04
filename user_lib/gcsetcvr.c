/*****************************************************************************
* gcsetcvr.c - Calculating Geometric Covering.                               *
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

typedef struct UserGCRgbToGrayStruct {
    IrtImgPixelStruct Rgb;
    MiscISCPixelType Gray;
} UserGCRgbToGrayStruct;

/* Set the transformation of RGB map to gray map. The first three nubmers   */
/* are the RGB original value and the forth is the gray destination value.  */
IRIT_STATIC_DATA const UserGCRgbToGrayStruct 
    USER_GC_COLOR_TO_GRAY_MAP[] = {
        { {0, 125, 0 }, 1 }, 
        { {255, 0, 0}, 0 }, 
        { {255, 255, 0}, 0 },
        { {255, 255, 255}, 0 }
    },
    USER_GC_VISMAP_TO_MAPPED_MAP[] = {
        { {0, 125, 0 }, 1 }, 
        { {255, 0, 0 }, 1 }, 
        { {255, 255, 0 }, 1 },
        { {255, 255, 255 }, 0 }
    };

IRIT_STATIC_DATA const int 
    USER_GC_COLOR_TO_GRAY_MAP_SIZE = sizeof(USER_GC_COLOR_TO_GRAY_MAP) / 
                                     sizeof(UserGCRgbToGrayStruct),
    USER_GC_VISMAP_TO_MAPPED_MAP_SIZE = sizeof(USER_GC_VISMAP_TO_MAPPED_MAP) / 
				        sizeof(UserGCRgbToGrayStruct);

static UserGCSolutionIndexStruct **UserGCSetSolutionSuggestion(
                                       UserGCProblemDefinitionStruct *Problem,
                                       UserGCObsPtGroupTypeStruct **OpGroups,
                                       int *SolutionOpsIndex,
                                       int SolutionOpsNum);
static MiscISCPixelType *UserGCFlatVisMap(
                                    int VisMapWidth,
                                    int VisMapHeight,
                                    const IrtImgPixelStruct *Pixels,
                                    MiscISCPixelType *VisMap,
                                    const UserGCRgbToGrayStruct *RgbToGrayMap,
                                    int SizeOfRgbToGraymap);
static int UserGCSolveSetcover(UserGCProblemDefinitionStruct *Problem, 
                               int **SolutionByIndex,
                               int *SolutionSize,
                               IrtRType *CoverPart);
static int UserGCSavePpmImageToFile(const char* FileName, 
                                    IrtImgPixelStruct* VisMap, 
                                    int Width, 
                                    int Height);
static IrtImgPixelStruct *UserGCLoadPpmImageFromFile(const char* FileName, 
                                                     int *Width, 
                                                     int *Height);
static int UserGCNullPrint(const char *Format, ...);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Return the indices of Suggestions in OpGroups wich appear in the         *
* solution and their visiblity maps.                                         *
*   The visibility map of each observation point points to one of VisMaps'   *
* visibility map and therefore shouldn't be freed.                           *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:           IN, The geometric covering problem to be solved.      *
*   OpGroups:          IN, The original suggestions used for solving the     *
*                          problem.                                          *
*   SolutionOpsIndex:  IN, Array with the indices of the suggestions which   *
*                          appear in the solution. The order of the          *
*                          suggestions (for which the index relay to) is     *
*                          starting from the suggestions in the first        *
*                          observation group of OpGroups and continue with   *
*                          the suggestions of the next observation group of  *
*                          OpGroups and so on.                               *
*   SolutionOpsNum:    IN, Number of observation points solving the problem  *
*                          (which is the number of indices in                *
*                          SolutionOpsIndex).                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   UserGCSolutionIndexStruct**: null terminated array of indices of         *
*       Suggestions in OpGroups which appears in the solution. Each elemnt   *
*       in the array contains the index of the observation point group, the  *
*       index of the suggestion in that observation point group and the      *
*       index given in SolutionOpsIndex.                                     *
*****************************************************************************/
static UserGCSolutionIndexStruct **UserGCSetSolutionSuggestion(
                                       UserGCProblemDefinitionStruct *Problem,
				       UserGCObsPtGroupTypeStruct **OpGroups,
				       int *SolutionOpsIndex,
				       int SolutionOpsNum)
{
    int i,
        SuggestionIndex = 0, 
        ObsPtGroupIndex = 0,
        SuggestionNum = UserGCGetOPsNum(Problem),
        *SuggestionIndices = (int *) IritMalloc(sizeof(int *) * SuggestionNum),
        *ObsPtGroupIndices = (int *) IritMalloc(sizeof(int *) * SuggestionNum);
    UserGCSolutionIndexStruct
        **Res = (UserGCSolutionIndexStruct **)
                        IritMalloc((SolutionOpsNum + 1) * 
                                   sizeof(UserGCSolutionIndexStruct *));

    Res[SolutionOpsNum] = NULL;

    for (i = 0; i <= SuggestionNum - 1; i++, SuggestionIndex++) {
        if (OpGroups[ObsPtGroupIndex] -> Suggestions[SuggestionIndex] == NULL){
            ObsPtGroupIndex++;
            SuggestionIndex = 0;
        }
        ObsPtGroupIndices[i] = ObsPtGroupIndex;
        SuggestionIndices[i] = SuggestionIndex;
    }

    for (i = 0; i <= SolutionOpsNum - 1; i++) {
        assert((SolutionOpsIndex[i] <= SuggestionNum) &&
	       "The solution contains solutions that doesn't appear in the original problem");

        Res[i] =  (UserGCSolutionIndexStruct*)
            IritMalloc(sizeof(UserGCSolutionIndexStruct));
        Res[i] -> ObsPtGroupIndex = ObsPtGroupIndices[SolutionOpsIndex[i]];
        Res[i] -> SuggestionIndex = SuggestionIndices[SolutionOpsIndex[i]];
        Res[i] -> Index = SolutionOpsIndex[i];
    }
    IritFree(SuggestionIndices);
    IritFree(ObsPtGroupIndices);

    return Res;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Turns an RGB visibilty maps into a gray map using RgbToGrayMap.          *
*   Print error message and exit if can't allocate memory (when VisMap is    *
* NULL).                                                                     *
*                                                                            *
* PARAMETERS:                                                                *
*   VisMapWidth:         IN, The width of the visibility map.                *
*   VisMapHeight:        IN, The height of the visibility map.               *
*   Pixels:              IN, The map.                                        *
*   VisMap:              OUT, Create the map in this buffer. If NULL,        *
                         allocates new buffer.                               *
*   RgbToGrayMap:        IN, Map from RGB color to gray color.               *
*   SizeOfRgbToGraymap:  IN, Size of RgbToGrayMap.                           *
*                                                                            *
* RETURN VALUE:                                                              *
*   MiscISCPixelType*: A gray visibility map (in a new allocated memory if   *
*                      VisMap is NULL).                                      *
*****************************************************************************/
static MiscISCPixelType *UserGCFlatVisMap(
				    int VisMapWidth,
				    int VisMapHeight,
				    const IrtImgPixelStruct *Pixels,
				    MiscISCPixelType *VisMap,
				    const UserGCRgbToGrayStruct *RgbToGrayMap,
				    int SizeOfRgbToGraymap)
{
    int i, j;

    if (VisMap == NULL) {
        VisMap = (MiscISCPixelType *) IritMalloc(VisMapWidth * VisMapHeight *
                                                 sizeof(MiscISCPixelType));    
        if (VisMap == NULL) {
            IRIT_WARNING_MSG_PRINTF("Error allocating memory in UserGCFlatVisMap.\nThis message will be followed by exit of the program (this is the time to use debugger).");
            exit(-1);
        }
    }

    for (i = 0; i <= VisMapWidth * VisMapHeight - 1; i++) {
        for (j = 0; j <= SizeOfRgbToGraymap - 1; j++) {
            if ((Pixels[i].r == RgbToGrayMap[j].Rgb.r) &&
                (Pixels[i].g == RgbToGrayMap[j].Rgb.g) &&
                (Pixels[i].b == RgbToGrayMap[j].Rgb.b)) {
                VisMap[i] = RgbToGrayMap[j].Gray;
                break;
            }
        }
    }
    return VisMap;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find the smallest set of pictures in VisMaps which cover the entire      *
* image. Pixels with empty/degenerated color are subtracted from the image so*
* that the image required to be cover doesn't include them.                  *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:         The geometric covering problem to be solved.            *
*   SolutionByIndex: OUT, The solution as indices of pictures by the order   *
*                    they were added to the solution.                        *
*   SolutionSize:    OUT, The size of the solution (size of SolutionByIndex).*
*   CoverPart:       OUT, The part of the image covered by the solution.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: FALSE if error occured.                                             *
*****************************************************************************/
static int UserGCSolveSetcover(UserGCProblemDefinitionStruct *Problem, 
                               int **SolutionByIndex,
                               int *SolutionSize,
                               IrtRType *CoverPart)
{
    MiscISCCalculatorPtrType Calc;
    int PicsNum, i,
        Res = FALSE;
    MiscISCPixelType 
        *GrayVisMap = NULL; 
    IrtImgPixelStruct
        *VisMap = NULL;

    PicsNum = UserGCGetOPsNum(Problem);
    /* Creating the calculator. */
    Calc = MiscISCNewCalculator(PicsNum, 
				Problem -> SolvingParams.VisMapWidth *
				Problem -> SolvingParams.VisMapHeight,
				MISC_ISC_BNW,
                                Problem -> DebugParams.Print);

    if (!Calc)
        return FALSE;

    /* Adding all the pictures to the calculator. */
    for (i = 0; i <= PicsNum - 1; i++) {
        IritFree(VisMap);
        VisMap = UserGCLoadVisMap2(Problem, i);
        /* Turn the visibility map to black and white in order to add it */
        /* to the calculator.                                            */
        GrayVisMap = UserGCFlatVisMap(Problem -> SolvingParams.VisMapWidth, 
				      Problem -> SolvingParams.VisMapHeight,
				      VisMap, GrayVisMap,
				      USER_GC_COLOR_TO_GRAY_MAP, 
				      USER_GC_COLOR_TO_GRAY_MAP_SIZE);
        MiscISCAddPicture(Calc, GrayVisMap);
    }

    /* Setting the image to be covered in the set cover problem.            */
    /* That image is a map with zeroes in all unmapped or degenerated       */
    /* locations.                                                           */
    GrayVisMap = UserGCFlatVisMap(Problem -> SolvingParams.VisMapWidth, 
				  Problem -> SolvingParams.VisMapHeight,
				  VisMap, GrayVisMap,
				  USER_GC_VISMAP_TO_MAPPED_MAP, 
				  USER_GC_VISMAP_TO_MAPPED_MAP_SIZE);
    MiscISCSetImageToCover(Calc, GrayVisMap);
    IritFree(GrayVisMap);
    IritFree(VisMap);

    /* Solving the set cover. */
    Problem -> DebugParams.Print("%s Starts solving the set cover.\n", 
	                         MiscISCGetTimeStamp(time(NULL), NULL, 
                                                     USER_GC_TIME_ZONE, 
                                                     FALSE));

    switch (Problem -> SolvingParams.SetCoverParams.Algorithm) {
        case USER_GC_GREEDY: {
            Res = MiscISCCalculateGreedy(Calc, SolutionByIndex, SolutionSize, 
				         CoverPart);
            break;
        }
        case USER_GC_EXHAUSTIVE: {
            Res = MiscISCCalculateExhaustive(Calc, 
                Problem -> SolvingParams.SetCoverParams.CoverLimit,
                Problem -> SolvingParams.SetCoverParams.SizeLimit,
                SolutionByIndex, SolutionSize, CoverPart);
            break;
        }
        case USER_GC_EXACT: {
            Res = MiscISCCalculateExact(Calc, 
                Problem -> SolvingParams.SetCoverParams.SizeLimit,
                SolutionByIndex, SolutionSize, CoverPart);
            break;
        }
        default: {
            assert("Unknown algorithm type." == 0);
        }
    }

    Problem -> DebugParams.Print("%s Finished solving the set cover.\n", 
	                         MiscISCGetTimeStamp(time(NULL), NULL, 
                                                     USER_GC_TIME_ZONE, 
                                                     FALSE));

    /* Freeing the calculator. */
    MiscISCFreeCalculator(Calc);

    return Res;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Saving image in a ppm 6 format.                                          M
*                                                                            *
* PARAMETERS:                                                                M
*   FileName: The name of the file to save (without an extension).           M
*   VisMap:   The visibility map to save.                                    M
*   Width:    Width of the map.                                              M
*   Height:   Height of the map.                                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   int: FALSE if failed saving (Error message is produced by the function). M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCSavePpmImageToFile, UserGCLoadPpmImageFromFile                     M
*****************************************************************************/
int UserGCSavePpmImageToFile(const char* FileName, 
                             IrtImgPixelStruct* VisMap, 
                             int Width, 
                             int Height)
{
    char FullFileName[IRIT_LINE_LEN_VLONG];
    int i;

    sprintf(FullFileName, "%s.ppm", FileName);

    IrtImgWriteSetType("ppm6");

    if (!IrtImgWriteOpenFile(NULL, FullFileName, FALSE, Width, Height)) {
        IRIT_WARNING_MSG_PRINTF("Can't open for writing file %s.\n", FileName);
        return FALSE;
    }

    for (i = 0; i <= Height - 1; i++)
        IrtImgWritePutLine(NULL, VisMap + i * Width);

    IrtImgWriteCloseFile();

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Loading image from a file into format of MiscISCPixelType.               M
*                                                                            M
* PARAMETERS:                                                                *
*   FileName: IN, The name of the file to load.                              M
*   Width:    OUT, The width of the image.                                   M
*   Height:   OUT, The height of the image.                                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct *: The loaded visibility map or NULL if failed (Error  M
*                        message is produced by the function).               M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCLoadPpmImageFromFile, UserGCSavePpmImageToFile                     M
*****************************************************************************/
IrtImgPixelStruct *UserGCLoadPpmImageFromFile(const char* FileName, 
                                              int *Width, 
                                              int *Height)
{
    int Alpha;
    IrtImgPixelStruct *Res;
    char FullFileName[IRIT_LINE_LEN_VLONG];

    sprintf(FullFileName, "%s.ppm", FileName);
    Alpha = FALSE;
    Res = IrtImgReadImage(FullFileName, Width, Height, &Alpha);
    (*Width)++;
    (*Height)++;
    if (Res == NULL)
        return NULL;

    return Res;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   A printf function which returns without doing anything.                  *
*                                                                            *
* PARAMETERS:                                                                *
*   None                                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: The number of printed characters. Always 0.                         *
*****************************************************************************/
static int UserGCNullPrint(const char *Format, ...)
{
    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Solves the geometric covering problem given by Problem and returns the   M
* solving observation points.                                                M
*   The geometric object and obastacles of Problem are going through further M
* processing and becoming compatible with the objects sent to the            M
* visibility map generator.                                                  M
*   The objects stored in Problem->GeoObj and Problem->Obstacles are         M
* destroyed and mustn't be accessed again. They are replaced by a processed  M
* version of themselves.                                                     M
*   Both objects shoul should already contain only object of type polygon    M
* connected only by pnext (not by U.Lst).                                    M
*   The function may fail, in which case errors will be anounced to the user M
* and it will return FALSE.                                                  M
*                                                                            M
* PARAMETERS:                                                                M
*   Problem:           IN OUT, The geometric covering problem to be solved.  M
*   SolutionOps:       OUT, Null terminated array of solving Suggestions.    M
*                      Each element in the array doesn't contain the         M
*                      suggestions itself but indices to the suggestion.     M
*                      The first two indices point to the location of the    M
*                      suggestion inside Problem -> ObsPtsGroups. Those are  M
*                      the index of the observation group and the index of   M
*                      the suggestion inside that obervation group.          M
*                      The third index is the index in the continues index   M
*                      of all suggestions together (The index in which the   M
*                      visibility map is saved to disk).                     M
*   CoverPart:         OUT, Will hold the part of the cover picture that we  M
*                      succeeded to cover.                                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   int: Return FALSE if encoutntered error.                                 M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserGCSolveGeoProblem                                                    M
*****************************************************************************/
int UserGCSolveGeoProblem(UserGCProblemDefinitionStruct *Problem,
                          UserGCSolutionIndexStruct *** SolutionOps,
                          IrtRType *CoverPart)
{
    int Res,
        PicsNum = -1, 
        *SolutionOpsIndex = NULL, 
        SolutionOpsNum = 0;

    if (setjmp(UserGCStartOfProcess)) {
        if (!Problem -> DebugParams.StoreVisMap)
            UserGCDeleteAllVisMaps(Problem);

        if ((!Problem -> DebugParams.StoreObjectsBeforeOPs) && 
            (Problem -> DebugParams.LoadObjectsFromDisk > 0))
            UserGCDeleteProcessedObjects(Problem);

        return FALSE;
    }

    if ((Problem -> DebugParams.LoadImageFunc == NULL) ||
        (Problem -> DebugParams.SaveImageFunc == NULL)) {
        strcpy(Problem -> DebugParams.VisMapExtension, "ppm");
        Problem -> DebugParams.LoadImageFunc = UserGCLoadPpmImageFromFile;
        Problem -> DebugParams.SaveImageFunc = UserGCSavePpmImageToFile;
    }

    if (Problem -> DebugParams.Print == NULL)
        Problem -> DebugParams.Print = UserGCNullPrint;

    /* Prepare the observation points. */
    if (!UserGCInterpretOPGroupsSuggestion(Problem))
        longjmp(UserGCStartOfProcess, 1);
    PicsNum = UserGCGetOPsNum(Problem); 

    /* Create the scene for render. */
    UserGCPrepareScene(Problem);

    /* Create all the visiblity maps. */
    UserGCCreateVisMaps(Problem);

    /* Solving the set cover. */
    if (!Problem -> DebugParams.DisableSetCover) {
        Res = UserGCSolveSetcover(Problem, &SolutionOpsIndex, 
            &SolutionOpsNum, CoverPart);
    }
    else {
        /* Don't solve, just load all pictures. */
        int i;

        SolutionOpsNum = PicsNum;
        SolutionOpsIndex = (int *)IritMalloc(sizeof(int) * SolutionOpsNum);
        for (i = 0; i <= SolutionOpsNum - 1; i++)
            SolutionOpsIndex[i] = i;
        Res = TRUE;
    }

    if (Res) {
        /* Creating the result from the solutions of the set cover. */
        *SolutionOps = UserGCSetSolutionSuggestion(Problem, 
            Problem -> ObsPtsGroups, SolutionOpsIndex, SolutionOpsNum);
        IritFree(SolutionOpsIndex);
    }
    else {
        *SolutionOps = (UserGCSolutionIndexStruct **) 
                             IritMalloc(sizeof(UserGCSolutionIndexStruct *));
        (*SolutionOps)[0] = NULL;
    }

    if (Problem -> DebugParams.LoadObjectsFromDisk > 1)
        UserGCLoadProcessedObjects(Problem, &Problem -> GeoObj, 
                                   &Problem -> Obstacles);

    if ((!Problem -> DebugParams.StoreObjectsBeforeOPs) && 
        (Problem -> DebugParams.LoadObjectsFromDisk > 0))
        UserGCDeleteProcessedObjects(Problem);

    if (!Problem -> DebugParams.StoreVisMap)
        UserGCDeleteAllVisMaps(Problem);

    return TRUE;
}
