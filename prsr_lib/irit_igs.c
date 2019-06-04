/*****************************************************************************
* Module to save IRIT data into IGES files.		        	     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber	 			 Ver 1.0, May 1998   *
*****************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <setjmp.h>
#include "inc_irit/irit_sm.h"
#include "prsr_loc.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/misc_lib.h"

#define IGES_MAX_LINE_LEN	1000
#define IGES_LINE_LEN		81
#define IGES_RECORD_TYPE_POS	72
#define IGES_RECORD_SEQNUM_POS	64
#define IGES_MAX_NAME_LEN	20
#define IGES_STATUS_BLANK	1000000
#define IGES_STATUS_PHY_DPND	  10000

#define IGES_BLANK_LINE(Line)	memset(Line, ' ', IGES_LINE_LEN)

#define IGES_BLANK_FILL(Line, c, Count) { \
				  int Len = (int) strlen(Line); \
				  memset(&Line[Len], ' ', \
					 IGES_LINE_LEN - Len); \
				  sprintf(&Line[IGES_RECORD_TYPE_POS], \
					  "%c%7d\n", c, Count++); \
			      }
#define IGES_PBLANK_FILL(Line, Entity, Count) { \
				  int Len = (int) strlen(Line); \
				  memset(&Line[Len], ' ', \
					 IGES_LINE_LEN - Len); \
				  sprintf(&Line[IGES_RECORD_SEQNUM_POS], \
					  "%8dP%7d\n", Entity, Count++); \
			      }

typedef struct IgesPSecStrRepStruct {
    struct IgesPSecStrRepStruct *Pnext;
    char *Str;
} IgesPSecStrRepStruct;

IRIT_STATIC_DATA char
    GlblLine[IGES_MAX_LINE_LEN + 1],
    GlblPSecLine[IGES_LINE_LEN];
IRIT_STATIC_DATA int
    GlblDSecCount = 1,
    GlblPSecCount = 1,
    GlblPSecStrRepCount = 0,
    GlblPSeqNum = 1,
    GlblSaveEucTrimCrvs = FALSE,
    GlblMoreFlag = FALSE;
IRIT_STATIC_DATA FILE
    *GlblIgesFile;
IRIT_STATIC_DATA IgesPSecStrRepStruct
    *GlblPSecStrRepHead = NULL,
    *GlblPSecStrRepTail = NULL;

static void IPMapModelInPlace(IPObjectStruct *PObj, IrtHmgnMatType Mat);
static void DumpDSection(IPObjectStruct *PObj, IrtHmgnMatType Mat);
static void DumpDOneObject(FILE *f, IPObjectStruct *PObj, int Status);
static IPObjectStruct *ProcessTrimSrf(FILE *f, TrimSrfStruct *TrimSrf);
static int GetObjColor(IPObjectStruct *PObj);
static void DumpPSection(IPObjectStruct *PObj, IrtHmgnMatType Mat);
static void DumpPOneObject(FILE *f, IPObjectStruct *PObj);
static int PrepPSectionObj(IPObjectStruct *PObj);
static void FreePSection(IPObjectStruct *PObj, IrtHmgnMatType Mat);
static int PrepPSectionPoly(IPPolygonStruct *Pl, int Polyline);
static void PrepPSecCtlPts(IrtRType **Points,
			   int Rational,
			   int Len,
			   char Delimiter);
static void PrepPSecVector(IrtRType *V, int Len, char Delimiter);
static void PrepPSecInteger(int i, char Delimiter);
static void PrepPSecReal(IrtRType r, char Delimiter);
static void PrepPSecFlush(void);
static void AppendPSecStr(char *Str);
static void FreePSecStr(IgesPSecStrRepStruct *p);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Dumps IRIT object as IGES 5.0 file.                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   PObject:       IritObject to dump as IGES 5.0 file.                      M
*   CrntViewMat:   The current viewing matrix to apply to the object.        M
*   IgesFileName:  Name of IGES file.					     M
*   Messages:      TRUE for warning messages.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:      TRUE if succesful, FALSE otherwise.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   IPIgesLoadFile, IPIgesSaveFileSetup                                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPIgesSaveFile                                                           M
*****************************************************************************/
int IPIgesSaveFile(const IPObjectStruct *PObject,
		   IrtHmgnMatType CrntViewMat,
		   const char *IgesFileName,
		   int Messages)
{
    int OldRefCountState,
	Count = 1;
    char Str[IRIT_LINE_LEN], IgesFName[IRIT_LINE_LEN];
    CagdBType
	OldBBoxTight = CagdTightBBox(TRUE);
    IrtHmgnMatType UnitMat;
    IPObjectStruct *PObj;

    GlblDSecCount = 1;
    GlblPSecCount = 1;
    GlblPSecStrRepCount = 0;
    GlblPSeqNum = 1;
    GlblMoreFlag = Messages;

    OldRefCountState = IPSetCopyObjectReferenceCount(FALSE);
    PObj = IPCopyObject(NULL, PObject, TRUE);
    PObj -> Pnext = NULL;
    IPSetCopyObjectReferenceCount(OldRefCountState);

    IPTraverseObjListHierarchy(PObj, CrntViewMat, IPMapModelInPlace);
    PObj = IPResolveInstances(PObj);

    if (IgesFileName != NULL && strncmp(IgesFileName, "-", 1) != 0) {
        strncpy(IgesFName, IgesFileName, IGES_MAX_NAME_LEN);
        IgesFName[IGES_MAX_NAME_LEN] = 0;
	if ((GlblIgesFile = fopen(IgesFileName, "w")) == NULL) {
	    if (GlblMoreFlag)
	        IRIT_INFO_MSG_PRINTF("Failed to open \"%s\".\n",
					IgesFileName);

	    CagdTightBBox(OldBBoxTight);

	    return FALSE;
	}
    }
    else {
	GlblIgesFile = stdout;
	IgesFName[0] = 0;
    }

    /* Print a comment header. */
    sprintf(GlblLine,
	    IRIT_EXP_STR("IGES file created via the IRIT solid modeling environment with irit2igs"));
    IGES_BLANK_FILL(GlblLine, 'S', Count);
    fprintf(GlblIgesFile, GlblLine);

    sprintf(GlblLine, "From %s", IgesFileName == NULL ? "-" : IgesFName);
    IGES_BLANK_FILL(GlblLine, 'S', Count);
    fprintf(GlblIgesFile, GlblLine);

    /* Print the global section. */
    Count = 1;
    sprintf(Str, "Irit %s", IRIT_VERSION);
    sprintf(GlblLine, "1H,,1H;,%dH%s,%dH%s,%dH%s,",
	    4, "Irit",
	    IgesFileName == NULL ? 6 : (int) strlen(IgesFName),
	    IgesFileName == NULL ? "stdout" : IgesFName,
	    (int) strlen(Str), Str);
    IGES_BLANK_FILL(GlblLine, 'G', Count);
    fprintf(GlblIgesFile, GlblLine);

    sprintf(GlblLine, "%dH%s,32,38,6,308,15,7HUnknown,1.0,2,2HMM,,,",
	    (int) strlen(Str), Str);
    IGES_BLANK_FILL(GlblLine, 'G', Count);
    fprintf(GlblIgesFile, GlblLine);
#if defined(sgi) || defined(SUN4) || defined(OS2GCC) || defined(__WINNT__) || defined(_INCLUDE_HPUX_SOURCE) || defined(OSF1DEC)
    {
	struct tm *Time;
	time_t t = time(0);
	Time = localtime(&t);

	sprintf(Str, "13H%02d%02d%02d.%02d%02d%02d",
		Time -> tm_year % 100, Time -> tm_mon, Time -> tm_mday,
		Time -> tm_hour, Time -> tm_min, Time -> tm_sec);
    }
#else
    strcpy(Str, "13HYYMMDD.HHMMSS");
#endif /* sgi || SUN4 || OS2GCC || __WINNT__ || _INCLUDE_HPUX_SOURCE || OSF1DEC */
    sprintf(GlblLine, "%s,0.000001,10.0,7HUnknown,7HUnknown,,,%s;",
	    Str, Str);
    IGES_BLANK_FILL(GlblLine, 'G', Count);
    fprintf(GlblIgesFile, GlblLine);

    MatGenUnitMat(UnitMat);

    /* First traversal - dump the iges 'D' section. */
    IPTraverseObjListHierarchy(PObj, UnitMat, DumpDSection);

    /* Second traversal - dump the iges 'P' section. */
    GlblPSecLine[0] = 0;
    GlblPSecCount = 1;
    IPTraverseObjListHierarchy(PObj, UnitMat, DumpPSection);

    /* Free auxiliary data. */
    IPTraverseObjListHierarchy(PObj, UnitMat, FreePSection);

    sprintf(GlblLine, "S%7dG%7dD%7dP%7d%40cT      1\n", 2, 3,
	    GlblDSecCount - 1, GlblPSecCount - 1, ' ');
    fprintf(GlblIgesFile, GlblLine);

    if (GlblIgesFile != stdout)
        fclose(GlblIgesFile);

    CagdTightBBox(OldBBoxTight);

    IPFreeObject(PObj);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Maps the object according to the given matrix, in place.                 *
*   If a model object, also map to a list of trimmed surfaces instead.       *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:      Object to map, in place, according to Mat.                    *
*   Mat:       Transformation matrix.                                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void IPMapModelInPlace(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    IPObjectStruct
	*Pnext = PObj -> Pnext,
	*MappedPObj = GMTransformObject(PObj, Mat);

    IPFreeObjectSlots(PObj);
    if (PObj -> ObjName)
        IritFree(PObj -> ObjName);

    IRIT_GEN_COPY(PObj, MappedPObj, sizeof(IPObjectStruct));
    PObj -> Pnext = Pnext;

    if (IP_IS_MODEL_OBJ(PObj)) {     /* Convert to a list of trimmed srfs. */
        TrimSrfStruct *TSrf,
	    *TrimSrfs = MdlCnvrtMdl2TrimmedSrfs(PObj -> U.Mdls);

	/* Chain trimming curves into closed loops. */
	for (TSrf = TrimSrfs; TSrf != NULL; TSrf = TSrf -> Pnext) {
	    TrimCrvStruct
	        *NewTCrvs = TrimChainTrimmingCurves2Loops(TSrf -> TrimCrvList);

	    TrimCrvFreeList(TSrf -> TrimCrvList);
	    TSrf -> TrimCrvList = NewTCrvs;
	}

	/* Convert PObj to list object, in place, & add trimmed surfaces. */
	IPReallocNewTypeObject(PObj, IP_OBJ_LIST_OBJ);

	while (TrimSrfs != NULL) {
	    IRIT_LIST_POP(TSrf, TrimSrfs);

	    IPListObjectAppend(PObj, IPGenTRIMSRFObject(TSrf));
	}
    }

    IPFreeObjectBase(MappedPObj);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the state of the IGES file save.                                    M
*                                                                            *
* PARAMETERS:                                                                M
*   SaveEucTrimCrvs:   TRUE to save Euclidean trimming curves as well.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IPIgesSaveFile                                                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   IPIgesSaveFileSetup                                                      M
*****************************************************************************/
void IPIgesSaveFileSetup(int SaveEucTrimCrvs)
{
    GlblSaveEucTrimCrvs = SaveEucTrimCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function of IPTraverseObjHierarchy. Called on every non        *
* list object found in hierarchy.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Non list object to handle.                                   *
*   Mat:        Transformation matrix to apply to this object.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpDSection(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    DumpDOneObject(GlblIgesFile, PObj, 0);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Dumps one object PObject to file f.                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   f:         File to dump object to.                                       *
*   PObj:      Object to dump to file f.                                     *
*   Status:    The status number of the object.                              *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpDOneObject(FILE *f, IPObjectStruct *PObj, int Status)
{
    int i,
	PSeqLen = 0;
    char
	ObjName[9];
    IPPolygonStruct *Pl;

    /* Prepare an object name. */
    strncpy(ObjName, IP_GET_OBJ_NAME(PObj), 8);
    for (i = (int) strlen(IP_GET_OBJ_NAME(PObj)); i < 8; i++)
	ObjName[i] = ' ';
    ObjName[8] = 0;

    switch (PObj -> ObjType) {
	case IP_OBJ_POLY:
	    PSeqLen = PrepPSectionObj(PObj);

	    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
		int PSeqLen1 = PrepPSectionPoly(Pl, IP_IS_POLYLINE_OBJ(PObj));

	        AttrSetIntAttrib(&Pl -> Attr, "_EntityNum", GlblDSecCount);
		fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
			IGS_ENTYPE_COPIOUS_DATA, GlblPSeqNum, Status, GlblDSecCount++);
		fprintf(f, "%8d       2%8d%8d      12       0       0%s       0D%7d\n",
			IGS_ENTYPE_COPIOUS_DATA, GetObjColor(PObj),
			PSeqLen1, ObjName, GlblDSecCount++);
		GlblPSeqNum += PSeqLen1;
	    }
	    break;
	case IP_OBJ_POINT:
	case IP_OBJ_VECTOR:
	case IP_OBJ_CTLPT:
	    PSeqLen = PrepPSectionObj(PObj);

	    AttrSetObjectIntAttrib(PObj, "_EntityNum", GlblDSecCount);
	    fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
		    IGS_ENTYPE_POINT, GlblPSeqNum, Status, GlblDSecCount++);
	    fprintf(f, "%8d       2%8d%8d       0       0       0%s       0D%7d\n",
		    IGS_ENTYPE_POINT, GetObjColor(PObj), PSeqLen,
		    ObjName, GlblDSecCount++);
	    break;
	case IP_OBJ_PLANE:
	    PSeqLen = PrepPSectionObj(PObj);

	    AttrSetObjectIntAttrib(PObj, "_EntityNum", GlblDSecCount);
	    fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
		    IGS_ENTYPE_PLANE, GlblPSeqNum, Status, GlblDSecCount++);
	    fprintf(f, "%8d       2%8d%8d       0       0       0%s       0D%7d\n",
		    IGS_ENTYPE_PLANE, GetObjColor(PObj), PSeqLen,
		    ObjName, GlblDSecCount++);
	    break;
	case IP_OBJ_CURVE:
	    PSeqLen = PrepPSectionObj(PObj);

	    AttrSetObjectIntAttrib(PObj, "_EntityNum", GlblDSecCount);
	    fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
		    IGS_ENTYPE_RATIONAL_BSPLINE_CURVE, GlblPSeqNum, Status, GlblDSecCount++);
	    fprintf(f, "%8d       2%8d%8d       0       0       0%s       0D%7d\n",
		    IGS_ENTYPE_RATIONAL_BSPLINE_CURVE,
		    GetObjColor(PObj), PSeqLen, ObjName, GlblDSecCount++);
	    break;
	case IP_OBJ_SURFACE:
	    PSeqLen = PrepPSectionObj(PObj);

	    AttrSetObjectIntAttrib(PObj, "_EntityNum", GlblDSecCount);
	    fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
		    IGS_ENTYPE_RATIONAL_BSPLINE_SURFACE, GlblPSeqNum, Status, GlblDSecCount++);
	    fprintf(f, "%8d       2%8d%8d       0       0       0%s       0D%7d\n",
		    IGS_ENTYPE_RATIONAL_BSPLINE_SURFACE,
		    GetObjColor(PObj), PSeqLen, ObjName, GlblDSecCount++);
	    break;
	case IP_OBJ_TRIMSRF:
	    if (PObj -> U.TrimSrfs == NULL) {
	        /* It is a trimming curve on a surface entity. */
	        PSeqLen = PrepPSectionObj(PObj);

		AttrSetObjectIntAttrib(PObj, "_EntityNum", GlblDSecCount);
	        fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
			IGS_ENTYPE_CURVE_ON_PARAM_SRF, GlblPSeqNum, Status, GlblDSecCount++);
	        fprintf(f, "%8d       2%8d%8d       0       0       0%s       0D%7d\n",
			IGS_ENTYPE_CURVE_ON_PARAM_SRF,
			GetObjColor(PObj), PSeqLen, ObjName, GlblDSecCount++);
	    }
	    else {
	        AttrSetObjectRefPtrAttrib(PObj, "_AuxGeom",
				          ProcessTrimSrf(f, PObj -> U.TrimSrfs));
		PSeqLen = PrepPSectionObj(PObj);

		AttrSetObjectIntAttrib(PObj, "_EntityNum", GlblDSecCount);
		fprintf(f, "%8d%8d               1       1       0       0        %08dD%7d\n",
			IGS_ENTYPE_TRIMMED_PARAM_SRF, GlblPSeqNum, Status, GlblDSecCount++);
		fprintf(f, "%8d       2%8d%8d       0       0       0%s       0D%7d\n",
			IGS_ENTYPE_TRIMMED_PARAM_SRF,
			GetObjColor(PObj), PSeqLen, ObjName, GlblDSecCount++);
	    }
	    break;
	case IP_OBJ_LIST_OBJ:
	    if (GlblMoreFlag)
	        IRIT_WARNING_MSG("Should not visit list objects at this point!\n");
	    break;
	case IP_OBJ_MATRIX:
	    break;
	case IP_OBJ_INSTANCE:
	    break;
	case IP_OBJ_STRING:
	    break;
	default:
	    break;
    }

    GlblPSeqNum += PSeqLen;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Process the sub entries of the trimmed surfaces.                         *
*                                                                            *
* PARAMETERS:                                                                *
*   f:           File to dump object to.                                     *
*   TrimSrf:     Trim surface to process its sub entries.                    *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:  List of auxiliary objects tha must be treated.        *
*****************************************************************************/
static IPObjectStruct *ProcessTrimSrf(FILE *f, TrimSrfStruct *TrimSrf)
{
    IrtRType MaxUWidth;
    TrimCrvStruct
	*TrimCrvs = TrimSrf -> TrimCrvList;
    IPObjectStruct *PObj, *PObjSrf, *PObjPrm, *PObjEuc, *PObjMaxUWidth,
	*PObjList = NULL;

    PObjSrf = IPGenSRFObject(CagdSrfCopy(TrimSrf -> Srf));
    DumpDOneObject(f, PObjSrf, IGES_STATUS_BLANK + IGES_STATUS_PHY_DPND);
    IRIT_LIST_PUSH(PObjSrf, PObjList);

    while (TrimCrvs != NULL) {
	CagdBBoxStruct BBox;

	if (TrimCrvs -> TrimCrvSegList -> Pnext != NULL) {
	    if (GlblMoreFlag)
	        IRIT_WARNING_MSG(
			"Only one trimming segment per curve is supported\n");
	    return NULL;
	}

	/* Convert the parametric space curve. */
	PObjPrm = IPGenCRVObject(CagdCrvCopy(TrimCrvs -> TrimCrvSegList -> UVCrv));
	DumpDOneObject(f, PObjPrm, IGES_STATUS_BLANK + IGES_STATUS_PHY_DPND);
	IRIT_LIST_PUSH(PObjPrm, PObjList);

	/* Convert the Euclidean space curve. */
	if (TrimCrvs -> TrimCrvSegList -> EucCrv != NULL)
	    PObjEuc = IPGenCRVObject(CagdCrvCopy(TrimCrvs -> TrimCrvSegList
						                    -> EucCrv));
	else
	    PObjEuc = IPGenCRVObject(TrimEvalTrimCrvToEuclid(TrimSrf,
					  TrimCrvs -> TrimCrvSegList -> UVCrv));
	if (GlblSaveEucTrimCrvs)
	    DumpDOneObject(f, PObjEuc, IGES_STATUS_BLANK + IGES_STATUS_PHY_DPND);
	IRIT_LIST_PUSH(PObjEuc, PObjList);

	/* Save a marker to the trimming curve on a surface entity. */
	PObj = IPGenTRIMSRFObject(NULL);
	AttrSetObjectRefPtrAttrib(PObj, "_TrimSrf", PObjSrf);
	AttrSetObjectRefPtrAttrib(PObj, "_PrmCrv", PObjPrm);
	AttrSetObjectRefPtrAttrib(PObj, "_EucCrv", PObjEuc);

	CagdCrvBBox(TrimCrvs -> TrimCrvSegList -> UVCrv, &BBox);
	AttrSetObjectRealAttrib(PObj, "_UWidth", BBox.Max[0] - BBox.Min[0]);

	DumpDOneObject(f, PObj, 0);
	IRIT_LIST_PUSH(PObj, PObjList);

	TrimCrvs = TrimCrvs -> Pnext;
    }

    /* Figure out the curve with the maximum U span and mark it. */
    PObjMaxUWidth = NULL;
    MaxUWidth = 0.0;
    for (PObj = PObjList; PObj != NULL; PObj = PObj -> Pnext) {
        if (IP_IS_TRIMSRF_OBJ(PObj)) {
	    IrtRType
	        UWidth = AttrGetObjectRealAttrib(PObj, "_UWidth");

	    if (UWidth > MaxUWidth) {
	        UWidth = MaxUWidth;
		PObjMaxUWidth = PObj;
	    }
	}
    }
    AttrSetObjectIntAttrib(PObjMaxUWidth, "_MaxUWidth", TRUE);

    PObjList = IPReverseObjList(PObjList);
    return PObjList;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Converts the color of the object to IGES color style.                    *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:     Object to extract its color.                                   *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:      IGES color index.                                              *
*****************************************************************************/
static int GetObjColor(IPObjectStruct *PObj)
{
    int R, G, B,
	Color = AttrGetObjectColor(PObj);

    if (AttrGetObjectRGBColor(PObj, &R, &G, &B)) {
	IRIT_STATIC_DATA short Colors[][3] = {
	    { 0,   0,   0   },  /* 0. default IGS color */
	    { 0,   0,   0   },  /* 1. IGS_COLOR_BLACK */
	    { 255, 0,   0   },  /* 2. IGS_COLOR_RED */
	    { 0,   255, 0   },  /* 3. IGS_COLOR_GREEN */
	    { 0,   0,   255 },  /* 4. IGS_COLOR_BLUE */
	    { 255, 255, 63  },  /* 5. IGS_COLOR_YELLOW */
	    { 255, 0,   255 },  /* 6. IGS_COLOR_MAGENTA */
	    { 0,   255, 255 },  /* 7. IGS_COLOR_CYAN */
	    { 255, 255, 255 }   /* 8. IGS_COLOR_WHITE */
	};
	int i,
	    MinDistSqr = 65536 * 3;

	for (i = 1; i <= 8; i++) {
	    if (IRIT_SQR(R - Colors[i][0]) +
		IRIT_SQR(G - Colors[i][1]) +
		IRIT_SQR(B - Colors[i][2]) < MinDistSqr) {
		MinDistSqr = IRIT_SQR(R - Colors[i][0]) +
			     IRIT_SQR(G - Colors[i][1]) +
			     IRIT_SQR(B - Colors[i][2]);
		Color = i;
	    }
	}

	return Color;
    }
    else
        return IGS_COLOR_WHITE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function of IPTraverseObjHierarchy. Called on every non        *
* list object found in hierarchy.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Non list object to handle.                                   *
*   Mat:        Transformation matrix to apply to this object.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpPSection(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    DumpPOneObject(GlblIgesFile, PObj);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Dumps one object PObject to file f.                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   f:         File to dump object to.                                       *
*   PObj:      Object to dump to file f.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpPOneObject(FILE *f, IPObjectStruct *PObj)
{
    int EntityNum = AttrGetObjectIntAttrib(PObj, "_EntityNum");
    char Line[IGES_LINE_LEN + 1];
    IgesPSecStrRepStruct
	*PSec = (IgesPSecStrRepStruct *) AttrGetObjectRefPtrAttrib(PObj,
								   "_PSec");

    if (IP_IS_TRIMSRF_OBJ(PObj)) {
        IPObjectStruct
	    *PTmp = (IPObjectStruct *) AttrGetObjectRefPtrAttrib(PObj,
								 "_AuxGeom");

	/* Dump the references to the trimming curves and trimmed srf first. */
	for ( ; PTmp != NULL; PTmp = PTmp -> Pnext)
	    DumpPOneObject(f, PTmp);
    }

    if (IP_IS_POLY_OBJ(PObj)) {
	IPPolygonStruct *Pl;

	for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	    IgesPSecStrRepStruct
		*P2Sec = (IgesPSecStrRepStruct *) AttrGetRefPtrAttrib(Pl -> Attr,
								      "_PSec");

	    EntityNum = AttrGetIntAttrib(Pl -> Attr, "_EntityNum");

	    while (P2Sec != NULL) {
		strcpy(Line, P2Sec -> Str);
		IGES_PBLANK_FILL(Line, EntityNum, GlblPSecCount);

		fprintf(GlblIgesFile, Line);

		P2Sec = P2Sec -> Pnext;
	    }
	}
    }
    else if (PSec != NULL) {
	while (PSec != NULL) {
	    strcpy(Line, PSec -> Str);
	    IGES_PBLANK_FILL(Line, EntityNum, GlblPSecCount);

	    fprintf(GlblIgesFile, Line);

	    PSec = PSec -> Pnext;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prepare the strings forming the P section of this object and store as      *
* attribute in object.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:      Object to dump to file f.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:       Number of lines in the P section.                             *
*****************************************************************************/
static int PrepPSectionObj(IPObjectStruct *PObj)
{
    int n;
    IrtRType UMin, UMax, VMin, VMax, TMin, TMax;
    IPObjectStruct *PTmp, *PTmp2;
    CagdCrvStruct *Crv;
    CagdSrfStruct *Srf;

    FreePSecStr(GlblPSecStrRepHead);
    GlblPSecStrRepHead = GlblPSecStrRepTail = NULL;
    GlblPSecStrRepCount = 0;

    switch (PObj -> ObjType) {
	case IP_OBJ_VECTOR:
	case IP_OBJ_CTLPT:
	    PTmp = IPCoerceObjectPtTypeTo(PObj, IP_OBJ_POINT);
	    PrepPSecInteger(IGS_ENTYPE_POINT, ',');
	    PrepPSecReal(PTmp -> U.Pt[0], ',');
	    PrepPSecReal(PTmp -> U.Pt[1], ',');
	    PrepPSecReal(PTmp -> U.Pt[2], ',');
	    PrepPSecInteger(0, ';');
	    PrepPSecFlush();
	    IPFreeObject(PTmp);
	    break;
	case IP_OBJ_POINT:
	    PrepPSecInteger(IGS_ENTYPE_POINT, ',');
	    PrepPSecReal(PObj -> U.Pt[0], ',');
	    PrepPSecReal(PObj -> U.Pt[1], ',');
	    PrepPSecReal(PObj -> U.Pt[2], ',');
	    PrepPSecInteger(0, ';');
	    PrepPSecFlush();
	    break;
	case IP_OBJ_PLANE:
	    PrepPSecInteger(IGS_ENTYPE_PLANE, ',');
	    PrepPSecReal(PObj -> U.Plane[0], ',');
	    PrepPSecReal(PObj -> U.Plane[1], ',');
	    PrepPSecReal(PObj -> U.Plane[2], ',');
	    PrepPSecReal(PObj -> U.Plane[3], ',');
	    PrepPSecInteger(0, ',');
	    PrepPSecInteger(0, ',');
	    PrepPSecInteger(0, ',');
	    PrepPSecInteger(0, ',');
	    PrepPSecInteger(0, ';');
	    PrepPSecFlush();
	    break;
	case IP_OBJ_CURVE:
	    Crv = CagdCoerceCrvTo(PObj -> U.Crvs,
			CAGD_IS_RATIONAL_CRV(PObj -> U.Crvs) ? CAGD_PT_P3_TYPE
							     : CAGD_PT_E3_TYPE,
				  FALSE);
	    if (CAGD_IS_BEZIER_CRV(Crv)) {
	        CagdCrvStruct
		    *TCrv = CagdCnvrtBzr2BspCrv(Crv);

		CagdCrvFree(Crv);
		Crv = TCrv;
	    }
	    if (CAGD_IS_PERIODIC_CRV(Crv)) {
	        CagdCrvStruct
		    *TCrv = CagdCnvrtPeriodic2FloatCrv(Crv);

		CagdCrvFree(Crv);
		Crv = TCrv;
	    }

	    PrepPSecInteger(IGS_ENTYPE_RATIONAL_BSPLINE_CURVE, ',');
	    PrepPSecInteger(Crv -> Length - 1, ',');
	    PrepPSecInteger(Crv -> Order - 1, ',');
	    PrepPSecInteger(0, ',');              /* Not necessarily planar. */
	    PrepPSecInteger(0, ',');              /* Open curve, bu default. */
	    PrepPSecInteger(CAGD_IS_RATIONAL_CRV(Crv) ? 0 : 1, ',');
	    PrepPSecInteger(Crv -> Periodic, ',');
	    PrepPSecVector(Crv -> KnotVector, Crv -> Length + Crv -> Order, ',');
	    PrepPSecCtlPts(Crv -> Points, CAGD_IS_RATIONAL_CRV(Crv),
			   Crv -> Length, ',');

	    CagdCrvDomain(Crv, &TMin, &TMax);
	    PrepPSecReal(TMin, ',');
	    PrepPSecReal(TMax, ';');

	    CagdCrvFree(Crv);
	    PrepPSecFlush();
	    break;
	case IP_OBJ_SURFACE:
	    Srf = CagdCoerceSrfTo(PObj -> U.Srfs,
			CAGD_IS_RATIONAL_SRF(PObj -> U.Srfs) ? CAGD_PT_P3_TYPE
							     : CAGD_PT_E3_TYPE,
					     FALSE);
	    if (CAGD_IS_BEZIER_SRF(Srf)) {
	        CagdSrfStruct
		    *TSrf = CagdCnvrtBzr2BspSrf(Srf);

		CagdSrfFree(Srf);
		Srf = TSrf;
	    }
	    if (CAGD_IS_PERIODIC_SRF(Srf)) {
	        CagdSrfStruct
		    *TSrf = CagdCnvrtPeriodic2FloatSrf(Srf);

		CagdSrfFree(Srf);
		Srf = TSrf;
	    }

	    PrepPSecInteger(IGS_ENTYPE_RATIONAL_BSPLINE_SURFACE, ',');
	    PrepPSecInteger(Srf -> ULength - 1, ',');
	    PrepPSecInteger(Srf -> VLength - 1, ',');
	    PrepPSecInteger(Srf -> UOrder - 1, ',');
	    PrepPSecInteger(Srf -> VOrder - 1, ',');
	    PrepPSecInteger(0, ','); /* Open curve, bu default in first dir. */
	    PrepPSecInteger(0, ',');/* Open curve, bu default in second dir. */
	    PrepPSecInteger(CAGD_IS_RATIONAL_SRF(Srf) ? 0 : 1, ',');
	    PrepPSecInteger(Srf -> UPeriodic ? 1 : 0, ',');
	    PrepPSecInteger(Srf -> VPeriodic ? 1 : 0, ',');

	    PrepPSecVector(Srf -> UKnotVector, Srf -> ULength + Srf -> UOrder,
			   ',');
	    PrepPSecVector(Srf -> VKnotVector, Srf -> VLength + Srf -> VOrder,
			   ',');
	    PrepPSecCtlPts(Srf -> Points, CAGD_IS_RATIONAL_SRF(Srf),
			   Srf -> ULength * Srf -> VLength, ',');

	    CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
	    PrepPSecReal(UMin, ',');
	    PrepPSecReal(UMax, ',');
	    PrepPSecReal(VMin, ',');
	    PrepPSecReal(VMax, ';');

	    CagdSrfFree(Srf);
	    PrepPSecFlush();
	    break;
	case IP_OBJ_TRIMSRF:
	    if (PObj -> U.TrimSrfs == NULL) {
		IPObjectStruct
		    *TrimSrfObj = (IPObjectStruct *)
				   AttrGetObjectRefPtrAttrib(PObj, "_TrimSrf"),
		    *PrmCrvObj = (IPObjectStruct *)
				    AttrGetObjectRefPtrAttrib(PObj, "_PrmCrv"),
		    *EucCrvObj = (IPObjectStruct *)
				    AttrGetObjectRefPtrAttrib(PObj, "_EucCrv");

	        /* Its a trimming curve on a surface entity. */
	        PrepPSecInteger(IGS_ENTYPE_CURVE_ON_PARAM_SRF, ',');

		/* Curve on a surface created in an unspecified way. */
	        PrepPSecInteger(0, ',');

	        PrepPSecInteger(AttrGetObjectIntAttrib(TrimSrfObj,
						       "_EntityNum"),
				',');
	        PrepPSecInteger(AttrGetObjectIntAttrib(PrmCrvObj,
						       "_EntityNum"),
				',');
		if (GlblSaveEucTrimCrvs)
		    PrepPSecInteger(AttrGetObjectIntAttrib(EucCrvObj,
							   "_EntityNum"),
				    ',');
		else
		    PrepPSecInteger(0, ',');

		/* Unspecified prefered representation. */
	        PrepPSecInteger(0, ';');
		PrepPSecFlush();
	    }
	    else {
	        /* Its a trimmed surface! */
	        PrepPSecInteger(IGS_ENTYPE_TRIMMED_PARAM_SRF, ',');
		if ((PTmp = (IPObjectStruct *)
		          AttrGetObjectRefPtrAttrib(PObj, "_AuxGeom")) == NULL)
		     break;

		/* Dump surface Data Entity index. */
		PrepPSecInteger(AttrGetObjectIntAttrib(PTmp, "_EntityNum"),
				',');

		/* Outer boundary is not necessarily surface boundary. */
	        PrepPSecInteger(1, ',');

		/* Number of trimming curves (Divide by 3 as we have Euc.    */
		/* crv entity, Param. crv entity, and trimming crv entity).  */
	        PrepPSecInteger(n = IPObjListLen(PTmp -> Pnext) / 3 - 1, ',');

		/* Find largest curve and make it the outer boundary curve. */
		for (PTmp2 = PTmp -> Pnext;
		     PTmp2 != NULL;
		     PTmp2 = PTmp2 -> Pnext) {
		    if (IP_IS_TRIMSRF_OBJ(PTmp2) &&
			AttrGetObjectIntAttrib(PTmp2, "_MaxUWidth") == TRUE) {
		        PrepPSecInteger(AttrGetObjectIntAttrib(PTmp2,
							       "_EntityNum"),
				       (char) (n == 0 ? ';' : ','));
			break;
		    }
		}

		/* Dump the references to the trimming curves. */
		for (PTmp2 = PTmp -> Pnext;
		     PTmp2 != NULL;
		     PTmp2 = PTmp2 -> Pnext) {
		    if (IP_IS_TRIMSRF_OBJ(PTmp2) &&
			AttrGetObjectIntAttrib(PTmp2, "_MaxUWidth") != TRUE)
		        PrepPSecInteger(AttrGetObjectIntAttrib(PTmp2,
							       "_EntityNum"),
				       (char) (--n == 0 ? ';' : ','));
		}
		PrepPSecFlush();
	    }
	    break;
	case IP_OBJ_LIST_OBJ:
	    if (GlblMoreFlag)
	        IRIT_WARNING_MSG("Should not visit list objects at this point!\n");
	    break;
	case IP_OBJ_MATRIX:
	    break;
	case IP_OBJ_INSTANCE:
	    break;
	case IP_OBJ_STRING:
	    break;
	default:
	    break;
    }

    AttrSetObjectRefPtrAttrib(PObj, "_PSec", GlblPSecStrRepHead);
    GlblPSecStrRepHead = GlblPSecStrRepTail = NULL;
    n = GlblPSecStrRepCount;
    GlblPSecStrRepCount = 0;

    return n;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function of IPTraverseObjHierarchy. Called on every non        *
* list object found in hierarchy.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Non list object to handle.                                   *
*   Mat:        Transformation matrix to apply to this object.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void FreePSection(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    IgesPSecStrRepStruct
	*PSec = (IgesPSecStrRepStruct *) AttrGetObjectRefPtrAttrib(PObj,
								   "_PSec");
    IPObjectStruct *PTmp,
        *PAuxGeom = (IPObjectStruct *) AttrGetObjectRefPtrAttrib(PObj,
								 "_AuxGeom");

    if (IP_IS_POLY_OBJ(PObj)) {
        IPPolygonStruct *Pl;

	for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	    IgesPSecStrRepStruct
		*P2Sec = (IgesPSecStrRepStruct *) AttrGetRefPtrAttrib(Pl -> Attr,
								      "_PSec");

	    FreePSecStr(P2Sec);
	}
    }
	
    FreePSecStr(PSec);

    if (PAuxGeom != NULL) {
        for (PTmp = PAuxGeom; PTmp != NULL; PTmp = PTmp -> Pnext)
	    FreePSection(PTmp, Mat);

        IPFreeObjectList(PAuxGeom);
	AttrFreeObjectAttribute(PObj, "_AuxGeom");
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Prepare the strings forming the P section of this polygon and store as     *
* attribute in polygon.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:      Object to dump to file f.                                     *
*   Polyline:  TRUE if polyline, FALSE for a polygon.                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:       Number of lines in the P section.                             *
*****************************************************************************/
static int PrepPSectionPoly(IPPolygonStruct *Pl, int Polyline)
{
    int n;
    IPVertexStruct
	*V = Pl -> PVertex;

    FreePSecStr(GlblPSecStrRepHead);
    GlblPSecStrRepHead = GlblPSecStrRepTail = NULL;
    GlblPSecStrRepCount = 0;

    PrepPSecInteger(IGS_ENTYPE_COPIOUS_DATA, ',');
    PrepPSecInteger(Polyline ? 12 : 2, ',');
    PrepPSecInteger(IPVrtxListLen(V), ',');
    do {
	PrepPSecReal(V -> Coord[0], ',');
	PrepPSecReal(V -> Coord[1], ',');
	PrepPSecReal(V -> Coord[2], (char) (V -> Pnext == NULL ? ';' : ','));

	V = V -> Pnext;
    }
    while (V != NULL && V != Pl -> PVertex);

    PrepPSecFlush();

    AttrSetRefPtrAttrib(&Pl -> Attr, "_PSec", GlblPSecStrRepHead);
    GlblPSecStrRepHead = GlblPSecStrRepTail = NULL;
    n = GlblPSecStrRepCount;
    GlblPSecStrRepCount = 0;

    return n;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Dumps a real value into the P section of entity.		             *
*                                                                            *
* PARAMETERS:                                                                *
*   Points:       Control points' array of curve/surface/etc.                *
*   Rational:	  If TRUE the freeform entity is rational.		     *
*   Len:          Length of each coordinate in the Points array.             *
*   Delimiter:    To use after the real.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrepPSecCtlPts(IrtRType **Points,
			   int Rational,
			   int Len,
			   char Delimiter)
{
    int i;
    IrtRType Weight;

    if (Rational)
	PrepPSecVector(Points[0], Len, Delimiter);
    else {
	for (i = 0; i < Len; i++)
	    PrepPSecReal(1.0, ',');
    }

    for (i = 0; i < Len; i++) {
	Weight = Rational ? Points[0][i] : 1.0;
	if (Weight < IRIT_UEPS)
	    Weight = IRIT_UEPS;

	PrepPSecReal(Points[1][i] / Weight, ',');
	PrepPSecReal(Points[2][i] / Weight, ',');
	PrepPSecReal(Points[3][i] / Weight,
		     (char) (i == Len - 1 ? Delimiter : ','));
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Dumps a vector of real values into the P section of entity.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   V:            Vector of reals to dump.                                   *
*   Len:          Length of vector of reals to dump.                         *
*   Delimiter:    To use after the real.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrepPSecVector(IrtRType *V, int Len, char Delimiter)
{
    int i;

    for (i = 0; i < Len; i++)
	PrepPSecReal(V[i], (char) (i == Len - 1 ? Delimiter : ','));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Dumps an integer value into the P section of entity.	             *
*                                                                            *
* PARAMETERS:                                                                *
*   i:            Integer value to dump                                      *
*   Delimiter:    To use after the integer.                                  *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrepPSecInteger(int i, char Delimiter)
{
    char Line[IGES_LINE_LEN];

    sprintf(Line, "%d%c", i, Delimiter);
    if (strlen(Line) + strlen(GlblPSecLine) >= IGES_RECORD_SEQNUM_POS - 1) {
	AppendPSecStr(GlblPSecLine);
	GlblPSecLine[0] = 0;
    }
    strcat(GlblPSecLine, Line);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Dumps a real value into the P section of entity.		             *
*                                                                            *
* PARAMETERS:                                                                *
*   r:            Real value to dump                                         *
*   Delimiter:    To use after the real.                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrepPSecReal(IrtRType r, char Delimiter)
{
    char Line[IGES_LINE_LEN];

    sprintf(Line, "%f%c", r, Delimiter);
    if (strlen(Line) + strlen(GlblPSecLine) >= IGES_RECORD_SEQNUM_POS - 1) {
	AppendPSecStr(GlblPSecLine);
	GlblPSecLine[0] = 0;
    }
    strcat(GlblPSecLine, Line);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Finish dumping P section of a given entity.                              *
*                                                                            *
* PARAMETERS:                                                                *
*   None	                                                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void PrepPSecFlush(void)
{
    if (!IRT_STR_ZERO_LEN(GlblPSecLine)) {
	AppendPSecStr(GlblPSecLine);
	GlblPSecLine[0] = 0;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Appends the given list to the currenly accumulated P section of this     *
* object.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   Str:        String to acculumate.                                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void AppendPSecStr(char *Str)
{
    if (GlblPSecStrRepHead == NULL) {
	GlblPSecStrRepTail =
	    GlblPSecStrRepHead = (IgesPSecStrRepStruct *)
	                             IritMalloc(sizeof(IgesPSecStrRepStruct));
    }
    else {
	GlblPSecStrRepTail -> Pnext =
	    (IgesPSecStrRepStruct *) IritMalloc(sizeof(IgesPSecStrRepStruct));
	GlblPSecStrRepTail = GlblPSecStrRepTail -> Pnext;
    }

    GlblPSecStrRepTail -> Pnext = NULL;
    GlblPSecStrRepTail -> Str = IritStrdup(Str);

    GlblPSecStrRepCount++;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Free the PSec strings link list.                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   p:   The PSec string linked list to free.                                *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void FreePSecStr(IgesPSecStrRepStruct *p)
{
    while (p != NULL) {
	IgesPSecStrRepStruct
	    *Pnext = p -> Pnext;

	IritFree(p -> Str);
	IritFree(p);

	p = Pnext;
    }
}
