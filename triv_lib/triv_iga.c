/*****************************************************************************
* Module to support IGA (isogeometric analysis) relayed queries of trivars.  *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber and Fady Massarwi		Ver 1.0, June 2012   *
*****************************************************************************/

#include "triv_loc.h"
#include "triv_loc_iga.h"

typedef struct TrivIGADataManager
{   
    TrivIGAArrangementStruct *Arrangements[TRIV_IGA_MAX_ARRANGEMENTS];
    TrivIGAArrangementID NextArrangementID;

    /* Materials is handled in the arrangement struct                       */
    /*    TrivIGAMaterialStruct *Materials[TRIV_IGA_MAX_MATERIALS];         */
    /* TrivIGAMaterialID NextMaterialID;				    */
    
    TrivIGATVStruct *IGATrivars[TRIV_IGA_MAX_TRIVARIATES];
    TrivIGATVID NextTrivariateID;
} TrivIGADataManager;

TrivIGADataManager
    *TrivIGAGlobalDataManager = NULL;

static void TrivIGAUpdateIGATVSlots(TrivIGAArrangementStruct *H,
				    TrivIGATVStruct *IGATV,
				    TrivTVStruct *TV,
				    CagdBType JustDerivatives);
static void TrivIGAResetTVIDS(TrivIGAArrangementStruct *H,
			      TrivIGATVStruct *IGATV);
static void TrivIGAFreeIGATVSlots(TrivIGATVStruct *IGATV, CagdBType FreeIDs);
static int TrivIGADeAllocArrangement(TrivIGAArrangementStruct *H);
static TrivIGAArrangementStruct *TrivIGAAllocateArrangement(void);
static TrivIGATVStruct *TrivIGATV2IGATV(TrivIGAArrangementStruct *H,
					const TrivTVStruct *TV);
static int TrivIGAPatchIdx2CtlPtIdx(TrivIGAArrangementStruct *H,
				    const CagdRType *KV,
				    int Len,
				    int PatchIdx);
static int TrivIGAUVW2RealUVW(TrivIGAArrangementStruct *H,
			      const TrivTVStruct *TV,
			      int IndexU,
			      int IndexV,
			      int IndexW,
			      CagdRType *U,
			      CagdRType *V,
			      CagdRType *W);
static const int *TrivIGAGetNeighboringCtlPtID(TrivIGAArrangementStruct *H,
					       const TrivTVStruct *TV,
					       const TrivTVStruct *SrcTV,
					       int SrcPtIndx);
static void TrivIGAAssignCtlPtsUniqueIDs(TrivIGAArrangementStruct *H,
					 TrivIGATVStruct *T);
static int TrivIGAUpdateOneCtlPt(TrivIGAArrangementStruct *H,
				 int ID,
				 CagdPointType PType,
				 const CagdRType *Coords,
				 int IsDelta);
static TrivIGACtrlPtStruct *TrivIGATVEvalBndryNormal(TrivIGAArrangementStruct *H,
						     const TrivTVStruct *TV,
						     CagdRType U,
						     CagdRType V,
						     CagdRType W,
						     TrivIGACtrlPtStruct *CtlPt);
static CagdBType TrivIGASameSrfs(const CagdSrfStruct *Srf,
				 const CagdSrfStruct *AdjSrf,
				 const TrivIGATVStruct *AdjTV,
				 int AdjBndryIdx,
				 TrivIGAAdjacencyInfoStruct *AdjInfo,
				 TrivIGAErrorType *Error);
static TrivIGAFieldType TrivIGAGetFieldValueType(const char *NamedType);
static char *TrivIGAExtractSubString(const char *Start, const char *End);
static int TrivIGASplitString(const char *Str, 
			      char DelimChar, 
			      char MapChar, 
			      char ***Keys, 
			      char ***Values,
			      int  *NumPairs);
static void TrivIGAAddBoundaryNode(TrivIGAArrangementStruct *H,
				   TrivIGATVStruct *IGATV,
				   int TVNodeIdx,
				   TrivIGANodeBoundaryType NodeBoundaryType,
				   const char *BoundaryAxisConditions,
				   CagdRType Value);
static void TrivIGAAddBoundaryNode(TrivIGAArrangementStruct *H,
				   TrivIGATVStruct *IGATV,
				   int TVNodeIdx,
				   TrivIGANodeBoundaryType NodeBoundaryType,
				   const char *BoundaryAxisConditions,
				   CagdRType Value);
static TrivIGAFieldType TrivIGAGetFieldValueType(const char *NamedType);
static void TrivIGAAddBoundaryNode(struct TrivIGAArrangementStruct *H,
				   TrivIGATVStruct *IGATV,
				   int TVNodeIdx,
				   TrivIGANodeBoundaryType NodeBoundaryType,
				   const char *BoundaryAxisConditions,
				   CagdRType Value);
static int TrivIGAGetBoundaryFaceByPtAux(const TrivTVStruct *TV,
					 const CagdPType Pt,
					 CagdRType *MinDistSqr);
static TrivIGADataManager *TrivIGADataManagerInit();

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returns substring without leading spaces.                                *
*                                                                            *
* PARAMETERS:                                                                *
*   Start:       Start position.			                     *
*   End:         End position.	  			                     *
*									     *
* RETURN VALUE:                                                              *
*   char *:     The allocated substring.                                     *
*****************************************************************************/
static char *TrivIGAExtractSubString(const char *Start, const char *End) 
{
    /* Don't include leading and ending spaces */
    char *Result;
    int ResultLen;

    while (' ' == (*Start)) {
	++Start;
    }

    while (' ' == (*End)) {
	--End;
    }

    ResultLen = (int) (End - Start + 1);
    if (ResultLen < 0)
	ResultLen = 0;
    Result = (char *) IritMalloc(sizeof(char) * (ResultLen + 1));
    strncpy(Result, Start, ResultLen);
    Result[ResultLen] = '\0';

    return Result;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   parses a string in the following format                                  *
*     key1 = value1, key2 = value 2, .....                                   *
*   and returns two arrays containing the keys and the values.               *
*                                                                            *
* PARAMETERS:                                                                *
*   Str:           String in the format described above.                     *
*   DelimChar:     Delimiter character (in the above example ',' ).          *
*   MapChar:       Mapping char (between key and value) (= in the example).  *
*   Keys:          Returned array of key strings.                            *
*   Values:        Returned Array of value strings.                          *
*   NumPairs:      Number of (key, value) pairs.                             *
*									     *
* RETURN VALUE:                                                              *
*   int:            TRUE on sucess, and FALSE otherwise.                     *
*****************************************************************************/
static int TrivIGASplitString(const char *Str, 
			      char DilimChar, 
			      char MapChar, 
			      char ***Keys, 
			      char ***Values,
			      int  *NumPairs)
{
    int i, StrLen,
        SearchState = 0,
        CurrentPair = 0;
    char CurrentChar, *TempStr, *LowerCaseTempStr;
    const char	
	*SubStrStart = Str,
	*StrIter = Str;


    if (Keys == NULL || Values == NULL || NumPairs == NULL)  {
	return FALSE;
    }

    /* Check how many dlimeters exists, and check validity. */
    *NumPairs = 0;
    StrLen = (int) strlen(Str);
    while ((CurrentChar = *StrIter++) != 0) {
	if (CurrentChar == DilimChar) {
	    if (SearchState == 0) {
	        /* Reached delimeter without facing map char. */
		return FALSE;
	    }
	    SearchState = 0;
	} 
	else if (CurrentChar == MapChar) {
	    if (SearchState == 0) {
		*NumPairs = *NumPairs + 1;
	    }
	    SearchState = 1;
	}
    }

    if (SearchState != 1) {
	return FALSE;
    }

    *Keys   = (char **) IritMalloc(sizeof(char *) * (*NumPairs));
    *Values = (char **) IritMalloc(sizeof(char *) * (*NumPairs));

    for (i = 0; i < *NumPairs; ++i) {
	(*Keys)[i] = NULL;
	(*Values)[i] = NULL;
    }

    SearchState = 0;
    StrIter = Str;
    while ((CurrentChar = *StrIter) != 0) {
	if (CurrentChar == DilimChar) {
	    (*Values)[CurrentPair++] = TrivIGAExtractSubString(SubStrStart,
							       StrIter - 1);
	    SubStrStart = StrIter + 1;
	    SearchState = 0;
	} 
	else if (CurrentChar == MapChar) {
	    if (SearchState == 0) {
		TempStr = TrivIGAExtractSubString(SubStrStart, StrIter - 1);
		LowerCaseTempStr = (char *)IritMalloc(sizeof(char) *
						       (strlen(TempStr) + 1));
		strcpy(LowerCaseTempStr, TempStr);
		IritStrLower(LowerCaseTempStr);

		(*Keys)[CurrentPair] = LowerCaseTempStr;
		if (TempStr)
		    IritFree(TempStr);
		SubStrStart = StrIter + 1;
	    }
	    SearchState = 1;
	}
	++StrIter;
    }

    (*Values)[CurrentPair] = TrivIGAExtractSubString(SubStrStart, StrIter-1);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new material and initializes it.			             M
*   Updates it otherwise.					             M
*                                                                            *
* PARAMETERS:                                                                M
*    MaterialStr: string representing the material in the following format   M
*    id = <id>, name = <name>, type = <type>, attrib1 = <attrib1>, ......    M
*    where id, name, type attributes are must.                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAMaterialStruct *: new created material 			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAParseMaterial                                                     M
*****************************************************************************/
TrivIGAMaterialStruct *TrivIGAParseMaterial(const char *MaterialStr)
{
    int NumPairs, i;
    char *LastDecimalPos, **Keys, **Values;
    CagdBType 
	ValidFormat = TRUE,
	NameExists  = FALSE,
	IdExists    = FALSE,
	TypeExists  = FALSE;
    TrivIGAMaterialStruct* 
	NewMaterial = (TrivIGAMaterialStruct*)
	                             IritMalloc(sizeof(TrivIGAMaterialStruct));

    NewMaterial -> NumProperties = 0;
    if (TrivIGASplitString(MaterialStr, ',', '=',
			   &Keys, &Values, &NumPairs)) {
	for (i = 0; i < NumPairs; ++i) {
	    if (!strcmp(Keys[i], "id")) {
		NewMaterial -> Id = (int) strtol(Values[i], 
		                                 &LastDecimalPos, 10);
		if (*LastDecimalPos != '\0') {
		    ValidFormat = FALSE;
		    break;
		}
		IdExists = TRUE;
	    } 
	    else if (!strcmp(Keys[i], "type")) {
		strcpy(NewMaterial -> Type, Values[i]);
		TypeExists = TRUE;
	    } 
	    else if (!strcmp(Keys[i], "name")) {
		strcpy(NewMaterial -> Name, Values[i]);
		NameExists = TRUE;
	    } 
	    else {
		strcpy(NewMaterial -> Properties[NewMaterial -> NumProperties].
							     Name, Keys[i]);
		strcpy(NewMaterial -> Properties[NewMaterial -> NumProperties].
							     Value, Values[i]);
		++NewMaterial -> NumProperties;
	    }
	}
    }

    if (!NameExists || !IdExists || !TypeExists) {
	ValidFormat = FALSE;
    }

    for (i = 0; i < NumPairs; ++i) {
	IritFree(Keys[i]);
	IritFree(Values[i]);
    }
    IritFree(Keys);
    IritFree(Values);

    if (!ValidFormat) {
	IritFree(NewMaterial);
	return NULL;
    }

    return NewMaterial;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Adds the given material to the given arrangement if not exists, and      M
* updates it otherwise.						             M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:   A handle on the IGA arrangement to process. 		     M
*   Material:	The material						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAMaterialID: Material ID or invalid if failed.                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAAddMaterial                                                       M
*****************************************************************************/
TrivIGAMaterialID TrivIGAAddMaterial(TrivIGAArrangementID ArgmntID,
				     TrivIGAMaterialStruct *Material)
{
    int i;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL || Material == NULL) {
	H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	return TRIV_IGA_INVALID_MATERIAL_ID;
    }

    for (i = 1; i < H -> NumMaterials; ++i) {
	if (H -> Materials[i] -> Id == Material -> Id) {
	    H -> LastError = TRIV_IGA_ERR_EXISTING_MATERIAL_ADDED; 
	    return TRIV_IGA_INVALID_MATERIAL_ID;
	}
    }
    H -> Materials[H -> NumMaterials] = Material;
    return H -> NumMaterials++;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Searches given arrangement and finds and returns the IGA TV of the given *
* TV, if any.								     *
*                                                                            *
* PARAMETERS:                                                                *
*   H:          A handle on the IGA arrangement to process. 		     *
*   TV:		The TV we seek its field.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   TrivIGATVStruct *:  Found field, or NULL if not found.                   *
*****************************************************************************/
static TrivIGATVStruct *TrivIGATV2IGATV(struct TrivIGAArrangementStruct *H,
					const TrivTVStruct *TV)
{
    const TrivIGAFieldStruct *F;

    if (H == NULL)
        return FALSE;
    if (H -> Fields == NULL || TV == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }

    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
        TrivIGATVStruct *T;

	for (T = F -> TVs; T != NULL; T = T -> Pnext) {
	    if (T -> TV == TV)
	        return T;
	}
    }

    return NULL;
}
			    
/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes the base (first) control point index of this PatchIdx in KV.    *
*   Indices of control points of PatchIdx in KV spans Base to Base+Order-1.  *
*                                                                            *
* PARAMETERS:                                                                *
*   H:        A handle on the IGA arrangement to process. 		     *
*   KV:       Knot vector to process and compute the first ctl pt index for. *
*   Len:      Length of KV.						     *
*   PatchIdx: The patch we seek its first ctl pt index, start from zero.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:   Base control point of this index, or -1 if error.                 *
*****************************************************************************/
static int TrivIGAPatchIdx2CtlPtIdx(struct TrivIGAArrangementStruct *H,
				    const CagdRType *KV,
				    int Len,
				    int PatchIdx)
{
    CagdRType *KnotVals;
    int i, Idx, *KnotMults;

    if (PatchIdx < 0) {
        H -> LastError = TRIV_IGA_ERR_INVALID_PATCH_INDEX;
        return -1;
    }

    KnotVals = (CagdRType *) IritMalloc(sizeof(CagdRType) * Len);
    KnotMults = (int *) IritMalloc(sizeof(int) * Len);

    Len = BspKnotsMultiplicityVector(KV, Len, KnotVals, KnotMults);
    if (PatchIdx > Len - 2) {
        H -> LastError = TRIV_IGA_ERR_INVALID_PATCH_INDEX;
        return -1;
    }

    for (i = 1, Idx = 0; i <= PatchIdx; i++)
        Idx += KnotMults[i];

    IritFree(KnotVals);
    IritFree(KnotMults);

    return Idx;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Maps the given parametric location, as an index of Bezier patch in TV    *
* in Dir and normalized [0, 1] parameter value inside that Bezier, into the  *
* real equivalent parameter of TV.                                           *
*                                                                            *
* PARAMETERS:                                                                *
*   H:         A handle on the IGA arrangement to process. 		     *
*   KV:        Knot vector to process and map domain.		             *
*   Len:       Length of KV.						     *
*   Order:     of the spline that uses this KV.				     *
*   Index:     The Bezier trivariate index in KV, starting from zero.        *
*   t:         Parameter values in [0, 1] in the designated Bezier interval. *
*              Updated in place to real parameter of (possibly B-spline) TV. *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:      TRUE if successful, FALSE otherwise.                           *
*****************************************************************************/
static int TrivIGAParam2RealParam(struct TrivIGAArrangementStruct *H,
				  const CagdRType *KV,
				  int Len,
				  int Order,
				  int Index,
				  CagdRType *t)
{
    int Base;
    CagdRType t1, t2;

    if (*t < 0.0 || *t > 1.0) {
        H -> LastError = TRIV_IGA_ERR_NOT_ZERO_ONE_PARAM;
        return FALSE;
    }

    Base = TrivIGAPatchIdx2CtlPtIdx(H, KV, Len, Index);
    if (Base == -1)
        return FALSE;

    /* Map the parameter values, now that we have the intervals. */
    t1 = KV[Base + Order - 1];
    t2 = KV[Base + Order];
    *t = IRIT_BLEND(t2, t1, *t);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Maps the given parametric location, as indices of Bezier patch in TV     *
* and normalized ([0, 1]^3) parameter values inside that Bezier, into the    *
* real equivalent parameters of TV.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   H:         A handle on the IGA arrangement to process. 		     *
*   TV:        Trivar to find the given parameters in its real space.        *
*   IndexU:    The Bezier trivariate index in U in TV, starting from zero.   *
*   IndexV:    The Bezier trivariate index in V in TV, starting from zero.   *
*   IndexW:    The Bezier trivariate index in W in TV, starting from zero.   *
*   U, V, W:   Parameter values to evaluate the designated Bezier trivar at. *
*              Given in [0, 1] of the relative Bezier domain and are updated *
*              in place to the real parameters of (possibly B-spline) TV.    *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:      TRUE if successful, FALSE otherwise.                           *
*****************************************************************************/
static int TrivIGAUVW2RealUVW(struct TrivIGAArrangementStruct *H,
			      const TrivTVStruct *TV,
			      int IndexU,
			      int IndexV,
			      int IndexW,
			      CagdRType *U,
			      CagdRType *V,
			      CagdRType *W)
{
    return TrivIGAParam2RealParam(H, TV -> UKnotVector,
				  TV -> ULength + TV -> UOrder,
				  TV -> UOrder, IndexU, U) &&
           TrivIGAParam2RealParam(H, TV -> VKnotVector,
				  TV -> VLength + TV -> VOrder,
				  TV -> VOrder, IndexV, V) &&
           TrivIGAParam2RealParam(H, TV -> WKnotVector,
				  TV -> WLength + TV -> WOrder,
				  TV -> WOrder, IndexW, W);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Creates a new structure to hold the IGA analysis arrangement.            M
*                                                                            *
* PARAMETERS:                                                                M
*   NewArgmntID:  If successful, will hold the ID of the newly allocated     M
*                 arrangement.						     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    0 if error or arrangement ID if valid.			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGANewArrangement                                                    M
*****************************************************************************/
int TrivIGANewArrangement(TrivIGAArrangementID *NewArgmntID)
{
    if (NewArgmntID == NULL) {
	return FALSE;
    }

    *NewArgmntID = TRIV_IGA_INVALID_ARRANGEMENT_ID;

    if (TrivIGAGlobalDataManager == NULL) {
	TrivIGAGlobalDataManager = TrivIGADataManagerInit();
	if (TrivIGAGlobalDataManager == NULL) {
	    return FALSE;
	}
    }
    
    *NewArgmntID = 
	TrivIGADataManagerAllocateArrangement(TrivIGAGlobalDataManager);
 
    return *NewArgmntID;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Allocate a new arrangement structure.                                    *
*                                                                            *
* PARAMETERS:                                                                *
*   None                                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   TrivIGAArrangementStruct *:  Allocated structure.                        *
*****************************************************************************/
static TrivIGAArrangementStruct *TrivIGAAllocateArrangement(void)
{
    int i;
    TrivIGAArrangementStruct
        *H = (TrivIGAArrangementStruct *)
                                IritMalloc(sizeof(TrivIGAArrangementStruct));

    H -> Fields = NULL;
    H -> UniqueGlblCtlPtIDMax = -1;
    H -> LastError = TRIV_IGA_ERR_NO_ERROR;
    H -> NumMaterials = 0;
    H -> BoundaryNodes = NULL;
    IRIT_ZAP_MEM(H -> Materials,
	         sizeof(TrivIGAMaterialStruct *) * TRIV_IGA_MAX_MATERIALS);

    /* For the 3 (U, V, W) direction, AlphaVal sets the (parametric)        */
    /* size-ratio between the last interval and the first.  NumOfIntervals  */
    /* specifies the number of intervals in each direction, with zero value */
    /* to deactivate.							    */
    for (i = 0; i < 3; i++) {
        H -> SeedingState.AlphaVal[i] = 1.0;
	H -> SeedingState.NumOfIntervals[i] = 0;
	H -> SeedingState.DmnMin[i] = H -> SeedingState.DmnMax[i] = 0.0;
    };

    return H;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returns the type of the field according to the field name                *
*                                                                            *
* PARAMETERS:                                                                *
*   NamedType:    The name field 		                             *
*                                                                            *
* RETURN VALUE:                                                              *
*  TrivIGAFieldType: TRIV_IGA_SCALAR_FIELD_TYPE for "pressure", "temp",      *
*                    and "temperature" and  TRIV_IGA_GEOMETRY_TYPE otherwise *
*****************************************************************************/
static TrivIGAFieldType TrivIGAGetFieldValueType(const char *NamedType)
{
    TrivIGAFieldType 
	Result = TRIV_IGA_GEOMETRY_TYPE;
    char 
	*LowerCaseNamedType = 
	        (char *) IritMalloc(sizeof(char) * (strlen(NamedType) + 1));

    strcpy(LowerCaseNamedType, NamedType);
    IritStrLower(LowerCaseNamedType);

    if (NamedType != NULL) {
	if (strstr(LowerCaseNamedType, "pressure") != NULL || 
	    strstr(LowerCaseNamedType, "temp") != NULL ||
	    strstr(LowerCaseNamedType, "temperature") != NULL) {
	    Result =  TRIV_IGA_SCALAR_FIELD_TYPE;
	}
    }

    if (LowerCaseNamedType != NULL)
	IritFree(LowerCaseNamedType);

    return Result;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Declares a new field.  Every field can consists of several (adjacent or  M
* not) trivariates and all subsequent call of TrivIGAAddTrivar2Field will be M
* placed in this new field.                                                  M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:          A handle on the IGA arrangement to process. 	     M
*   FieldAttributes:   Geometry or scalar/vector data fields.                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGANewField                                                          M
*****************************************************************************/
int TrivIGANewField(TrivIGAArrangementID ArgmntID,
		    const char *FieldAttributes)
{
    int NumPairs, i,
        MaterialID = -1;
    char *LastDecimalPos, **Keys, **Values;
    CagdBType 
	ValidFormat	= TRUE,
	MaterialExists  = FALSE,
	IdExists	= FALSE,
	TypeExists	= FALSE;
    TrivIGAFieldStruct *F;
    TrivIGAMaterialStruct 
	*FMaterial = NULL;
    struct TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL) {
        return FALSE;
    }

    if (!TrivIGASplitString(FieldAttributes, ',', '=',
			    &Keys, &Values, &NumPairs)) {
	H -> LastError = TRIV_IGA_ERR_INVALID_FIELD_ATTRIBUTES_FORMAT;
	return FALSE;
    }

    F = (TrivIGAFieldStruct *) IritMalloc(sizeof(TrivIGAFieldStruct));
    F -> NumProperties = 0;
    F -> TVs = NULL;
    F -> NamedType[0] = '\0';

    for (i = 0; i < NumPairs; ++i) {
	if (!strcmp(Keys[i], "id")) {
	    F -> ID = (int) strtol(Values[i], &LastDecimalPos, 10);
	    if (*LastDecimalPos != '\0') {
		ValidFormat = FALSE;
		break;
	    }
	    IdExists = TRUE;
	} 
	else if (!strcmp(Keys[i], "type")) {
	    strncpy(F -> NamedType , Values[i], TRIV_IGA_MAX_FIELD_TYPE_LEN - 1);
	    TypeExists = TRUE;
	} 
	else if (!strcmp(Keys[i], "material")) {
	    MaterialID = (int) strtol(Values[i], &LastDecimalPos, 10);
	    if (*LastDecimalPos != '\0') {
		ValidFormat = FALSE;
		break;
	    }
	    MaterialExists = TRUE;
	} 
	else {
	    assert(TRIV_IGA_MAX_XML_PROPERTY_LEN > F -> NumProperties);
	    strcpy(F -> Properties[F -> NumProperties].Name, Keys[i]);
	    strcpy(F -> Properties[F -> NumProperties].Value, Values[i]);
	    ++F -> NumProperties;
	}
    }

    if (!ValidFormat || !IdExists || !MaterialExists || !TypeExists) {
	IritFree(F);
	H -> LastError = TRIV_IGA_ERR_INVALID_FIELD_ATTRIBUTES_FORMAT;
	return FALSE;
    }

    for (i = 0; i < H -> NumMaterials; ++i) {
	if (H -> Materials[i] && H -> Materials[i] -> Id == MaterialID) {
	    FMaterial = H -> Materials[i];
	    break;
	}
    }

    if (!FMaterial) {
	IritFree(F);
	H -> LastError = TRIV_IGA_ERR_UNDEFINED_FIELD_MATERIAL;
	return FALSE;
    }
  
    F -> ValueType = TrivIGAGetFieldValueType(F -> NamedType);    
    F -> Material  = FMaterial;

    IRIT_LIST_PUSH(F, H -> Fields);

    for (i = 0; i < NumPairs; ++i) {
	IritFree(Keys[i]);
	IritFree(Values[i]);
    }
    IritFree(Keys);
    IritFree(Values);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Inserts a new trivariate into the IGA arrangement, in the current field, M
* in place.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.		     M
*   TV:          The Trivar to insert into the currently defined field, in   M
*                place.							     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVStruct *:   The allocated IGA TV, or NULL if error.	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGANewTV                                                             M
*****************************************************************************/
TrivIGATVStruct *TrivIGANewTV(TrivIGAArrangementID ArgmntID, TrivTVStruct *TV)
{
    TrivIGATVStruct *T;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL)
        return NULL;
    if (H -> Fields == NULL || TV == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }
    if (TRIV_IS_BEZIER_TV(TV)) {
	TrivTVStruct
	    *TVTmp = TrivCnvrtBzr2BspTV(TV);

	TrivTVFree(TV);
	TV = TVTmp;
    }

    T = (TrivIGATVStruct *) IritMalloc(sizeof(TrivIGATVStruct));
    TrivIGAUpdateIGATVSlots(H, T, TV, FALSE);
    T -> Pnext = NULL;

    /* Add new TV as last element: */
    if (H -> Fields -> TVs == NULL)
	H -> Fields -> TVs = T;
    else
        ((TrivIGATVStruct *) CagdListLast(H -> Fields -> TVs)) -> Pnext = T;

    return T;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Update all the slots of IGATV (which is assumed to hold no old data).    *
* using new trivar TV (which is used in place).                              *
*                                                                            *
* PARAMETERS:                                                                *
*   H:       The IGA arrangement to process.				     *
*   IGATV:   Structure to update all its slots.                              *
*   TV:      New TV to update IGATV with.                                    *
*   JustDerivatives:  TRUE to just recompute all partial derivative.         *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void TrivIGAUpdateIGATVSlots(TrivIGAArrangementStruct *H,
				    TrivIGATVStruct *IGATV,
				    TrivTVStruct *TV,
				    CagdBType JustDerivatives)
{
    int i,
	Len = TRIV_TV_UPT_LST_LEN(TV) *
	      TRIV_TV_VPT_LST_LEN(TV) *
	      TRIV_TV_WPT_LST_LEN(TV);

    IGATV -> TV = TV;

    /* Compute derivatives. */
    IGATV -> DuTV = TrivTVDeriveScalar(TV, TRIV_CONST_U_DIR);
    IGATV -> DvTV = TrivTVDeriveScalar(TV, TRIV_CONST_V_DIR);
    IGATV -> DwTV = TrivTVDeriveScalar(TV, TRIV_CONST_W_DIR);

    IGATV -> Du2TV = TrivTVDeriveScalar(IGATV -> DuTV, TRIV_CONST_U_DIR);
    IGATV -> Dv2TV = TrivTVDeriveScalar(IGATV -> DvTV, TRIV_CONST_V_DIR);
    IGATV -> Dw2TV = TrivTVDeriveScalar(IGATV -> DwTV, TRIV_CONST_W_DIR);

    IGATV -> DuDvTV = TrivTVDeriveScalar(IGATV -> DuTV, TRIV_CONST_V_DIR);
    IGATV -> DuDwTV = TrivTVDeriveScalar(IGATV -> DuTV, TRIV_CONST_W_DIR);
    IGATV -> DvDwTV = TrivTVDeriveScalar(IGATV -> DvTV, TRIV_CONST_W_DIR);

    if (JustDerivatives)
	return;

    IGATV -> Field = H -> Fields;     /* First field in H is current field. */
    IGATV -> CtlPtsIDs = (TrivIGACtlPtUniqueIDsStruct *)
                       IritMalloc(sizeof(TrivIGACtlPtUniqueIDsStruct) * Len);
    for (i = 0; i < 6; i++)
        IGATV -> Neighbors[i].AdjTV = NULL;

    TrivIGAResetTVIDS(H, IGATV);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Reset the IDs of all control points in the given IGA TV.                 *
*                                                                            *
* PARAMETERS:                                                                *
*   H:       The IGA arrangement to process.				     *
*   IGATV:   Structure to reset its control points IDs.                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void TrivIGAResetTVIDS(TrivIGAArrangementStruct *H,
			      TrivIGATVStruct *IGATV)
{
    int i,
	Len = TRIV_TV_UPT_LST_LEN(IGATV -> TV) *
	      TRIV_TV_VPT_LST_LEN(IGATV -> TV) *
	      TRIV_TV_WPT_LST_LEN(IGATV -> TV);

    /* Clear indexing setting. */
    for (i = 0; i < Len; i++)
        IGATV -> CtlPtsIDs[i].ID = -1;
    IGATV -> UniqueCtlPtIDMin = IGATV -> UniqueCtlPtIDMax = -1;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Inserts a new trivariate into the IGA arrangement, in the current field. M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.		     M
*   ExistingTV:  Current TV in the arrangement to update. Will be freed.     M
*   NewTV:       New TV to replace the existing TV.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVStruct *:  Reference to the updated IGA TV, or NULL if error.   M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAUpdateTV                                                          M
*****************************************************************************/
TrivIGATVStruct *TrivIGAUpdateTV(TrivIGAArrangementID ArgmntID,
				 TrivTVStruct *ExistingTV,
				 TrivTVStruct *NewTV)
{
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGATVStruct *IGATV;

    if (H == NULL)
        return NULL;

    IGATV = TrivIGATV2IGATV(H, ExistingTV);
    if (IGATV == NULL)
        return NULL;

    TrivIGAFreeIGATVSlots(IGATV, FALSE);
    TrivIGAUpdateIGATVSlots(H, IGATV, NewTV, TRUE);
    TrivTVFree(ExistingTV);

    return IGATV;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Search the given TV for control point Pt and return its index.           *
*                                                                            *
* PARAMETERS:                                                                *
*   H:         A handle on the IGA arrangement to process.		     *
*   TV:        Trivar to search for a control point similar to CtlPt.        *
*   SrcTV:     The current TV we process and scan its neighbors.             *
*                Can be the same as TV, in which case identical neighboring  *
*              ctlpts are scanned at lower indices only.		     *
*   SrcPtIndx: Index of the current ctlpt in the current TV we scan its      *
*              neighbors.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int *:    Pointer to the indices vector of this control point, or NULL   *
*             if not found.						     *
*****************************************************************************/
static const int *TrivIGAGetNeighboringCtlPtID(TrivIGAArrangementStruct *H,
					       const TrivTVStruct *TV,
					       const TrivTVStruct *SrcTV,
					       int SrcPtIndx)
{
    int i, j, k, l, Idx, IsNotRational, MaxCoord, ULen1, VLen1, WLen1;
    TrivIGATVStruct *T;

    if (TV == NULL || SrcTV == NULL)
	return NULL;

    if (TV -> PType != SrcTV -> PType) {
	return NULL;
    }

    IsNotRational = !CAGD_IS_RATIONAL_PT(TV -> PType);
    MaxCoord = CAGD_NUM_OF_PT_COORD(TV -> PType);
    ULen1 = TV -> ULength - 1;
    VLen1 = TV -> VLength - 1;
    WLen1 = TV -> WLength - 1;
  
    if (H == NULL)
        return NULL;

    T = TrivIGATV2IGATV(H, TV);
    if (T == NULL)
        return NULL;

    for (k = Idx = 0; k <= WLen1; k++) {
        for (j = 0; j <= VLen1; j++) {
	    for (i = 0; i <= ULen1; i++, Idx++) {
	        /* Scanning at the same TV with higher indices only. */
	        if (TV == SrcTV && Idx >= SrcPtIndx)
		    continue;

		for (l = IsNotRational; l <= MaxCoord; l++) {
		    if (!IRIT_APX_EQ(TV -> Points[l][Idx],
				     SrcTV -> Points[l][SrcPtIndx]))
		        break;
		}

		if (l > MaxCoord)
		    return &(T -> CtlPtsIDs[Idx].ID);
	    }
	}
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Assign unique IDs to the control points of given IGA trivar, possibly    *
* with neighbors (assigning same IDs if shared).			     *
*                                                                            *
* PARAMETERS:                                                                *
*   H:   A handle on the IGA arrangement to process.			     *
*   T:   Trivar to assign IDs uniquely for its ctl pts, taking into account  *
*        neighbors.						             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void TrivIGAAssignCtlPtsUniqueIDs(TrivIGAArrangementStruct *H,
					 TrivIGATVStruct *T)
{
    TrivTVStruct
        *TV = T -> TV;
    int PtIdx, i, j, k, l, Id,
        ULen1 = TV -> ULength - 1,
        VLen1 = TV -> VLength - 1,
        WLen1 = TV -> WLength - 1;
			
    T -> UniqueCtlPtIDMin = Id = H -> UniqueGlblCtlPtIDMax + 1;

    for (k = PtIdx = 0; k <= WLen1; k++) {
        for (j = 0; j <= VLen1; j++) {
	    for (i = 0; i <= ULen1; i++, PtIdx++) {
		CagdCtlPtStruct CtlPt;

		if (T -> CtlPtsIDs[PtIdx].ID >= 0) {
		    /* This control point was already assigned IDs from    */
		    /* its neighbor - skip it.				   */
		    continue;
		}
		
		CtlPt.PtType = TV -> PType;
		CagdCoercePointTo(CtlPt.Coords, CtlPt.PtType, TV -> Points,
				  TRIV_MESH_UVW(TV, i, j, k), TV -> PType);

		/* Go over all trivars and assign unique IDs to all control */
		/* points.						    */
		/*   If a boundary control point and has a neighboring      */
		/* control point in an adjacent TV that is identical, make  */
		/* them have the same ID.				    */
		for (l = 0; l < 6; l++) {
		    const int *NeighboringID;
		    const TrivTVStruct *AdjTV;

		    /* If there is not neighbor on this face, or neighbor   */
		    /* is in a different function space, skip it.	    */
		    if ((AdjTV = T -> Neighbors[l].AdjTV) == NULL ||
			!T -> Neighbors[l].SameSpace)
			continue;

		    NeighboringID = TrivIGAGetNeighboringCtlPtID(H, AdjTV,
								 TV, PtIdx);
		    if (NeighboringID != NULL && *NeighboringID >= 0) {
			T -> CtlPtsIDs[PtIdx].ID = *NeighboringID;
			break;
		    }
		}

		/* Assign new ID to this control point. */
		if (T -> CtlPtsIDs[PtIdx].ID < 0)
		    T -> CtlPtsIDs[PtIdx].ID = Id++;
	    }
	}
    }

    T -> UniqueCtlPtIDMax = H -> UniqueGlblCtlPtIDMax = Id - 1;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A function to signal the end of the initialization process - insertion   M
* of trivariates into (the fields of) the IGA arrangement.                   M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:   A handle on the IGA arrangement to process.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAArrangementComplete                                               M
*****************************************************************************/
int TrivIGAArrangementComplete(TrivIGAArrangementID ArgmntID)
{
    TrivIGAFieldStruct *F;
    TrivIGATVStruct *T;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL)
        return FALSE;

    /* Go over all trivars and compute neighborhoods. */
    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (T = F -> TVs; T != NULL; T = T -> Pnext) {
	    TrivIGAGetFaceNeighboringTVs(ArgmntID, T -> TV, T -> Neighbors);
	}
    }

    /* Go over all trivars and assign unique IDs to all control points.    */
    H -> UniqueGlblCtlPtIDMax = -1;  /* Reset the numbering to start over. */
    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (T = F -> TVs; T != NULL; T = T -> Pnext) {
	    TrivIGAResetTVIDS(H, T);
	}
    }
    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (T = F -> TVs; T != NULL; T = T -> Pnext) {
	    TrivIGAAssignCtlPtsUniqueIDs(H, T);
	}
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A debugging function to dump to stderr, the content of the given TV.     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.		     M
*   TV:          The Trivar to print to stderr.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAPrintTVContent                                                    M
*****************************************************************************/
int TrivIGAPrintTVContent(TrivIGAArrangementID ArgmntID,
			  const TrivTVStruct *TV)
{
    int i, j, MaxCoord;
    TrivIGATVStruct *ITV;
    TrivIGACtlPtUniqueIDsStruct *CtlPtsIDs;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL) {
	return FALSE;
    }

    if (TV == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return FALSE;
    }

    MaxCoord = CAGD_NUM_OF_PT_COORD(TV -> PType);
    ITV = TrivIGATV2IGATV(H, TV);
    if (ITV == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return FALSE;
    }
    CtlPtsIDs = ITV -> CtlPtsIDs;

    fprintf(stderr, "\n********** Dumping Trivar 0x%08ld: **********\n\n",
	    (long int) TV);

    if (TV -> GType != TRIV_TVBSPLINE_TYPE) {
        fprintf(stderr, "Trivar is not a Bspline trivar.  Ignored\n");
	return FALSE;
    }

    fprintf(stderr, "[TRIVAR BSPLINE %d %d %d %d %d %d %c%c\n",
	    TV -> ULength, TV -> VLength, TV -> WLength,
	    TV -> UOrder, TV -> VOrder, TV -> WOrder,
	    CAGD_IS_RATIONAL_PT(TV -> PType) ? 'P' : 'E',
	    MaxCoord + '0');

    /* Put out the knot vectors: */
    for (i = 0; i < 3; i++) {
        int Len, Periodic;
	CagdRType *KnotVector;

	switch (i) {
	    case 0:
	        KnotVector = TV -> UKnotVector;
		Periodic = TV -> UPeriodic;
		Len = TV -> ULength + TV -> UOrder
					+ (Periodic ? TV -> UOrder - 1 : 0);
		break;
	    case 1:
	        KnotVector = TV -> VKnotVector;
		Periodic = TV -> VPeriodic;
		Len = TV -> VLength + TV -> VOrder
					+ (Periodic ? TV -> VOrder - 1 : 0);
		break;
	    default:
	    case 2:
	        KnotVector = TV -> WKnotVector;
		Periodic = TV -> WPeriodic;
		Len = TV -> WLength + TV -> WOrder
					+ (Periodic ? TV -> WOrder - 1 : 0);
		break;
	}

	fprintf(stderr, Periodic ? "    [KVP" : "    [KV");
	for (j = 0; j < Len; j++) {
	     if (j && j % 5 == 0) {
	         fprintf(stderr, "\n\t");
	     }
	     fprintf(stderr, " %9.6f", KnotVector[j]);
	}
	fprintf(stderr, "]\n");
    }

    /* Put out the control mesh. */
    for (i = 0;
	 i < TV -> VLength * TV -> ULength * TV -> WLength;
	 i++) {
        if (i && i % TV -> ULength == 0)
	    fprintf(stderr, "\n");          /* Put empty line between raws. */
	if (i && i % TV -> UVPlane == 0)
	    fprintf(stderr, "\n");             /* Put 2 lns between planes. */

	fprintf(stderr, "    [");
	if (CAGD_IS_RATIONAL_PT(TV -> PType))
	    fprintf(stderr, "%9.6f ", TV -> Points[0][i]);
	for (j = 1; j <= MaxCoord; j++) {
	    fprintf(stderr, "%9.6f", TV -> Points[j][i]);
	    if (j < MaxCoord)
	        fprintf(stderr, " ");
	}

	fprintf(stderr, "]  \tIDs: {");

	for (j = 1; j <= MaxCoord; j++) {
	    fprintf(stderr, " %4d", CtlPtsIDs[i].ID);
	}
	fprintf(stderr, " }\n");
    }

    fprintf(stderr, "]\n");

    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the maximal IDS used in this arrangement.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:     Vector of maximal IDs,                                        M
*              (Max CtlPts ID, Max Trivars ID, Max Arrangement ID) if        M
*              successful, NULL otherwise.  Allocated statically.            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetGlblMaxIDs                                                     M
*****************************************************************************/
int *TrivIGAGetGlblMaxIDs(TrivIGAArrangementID ArgmntID)
{
    IRIT_STATIC_DATA int
        IDs[3];
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H != NULL && TrivIGAGlobalDataManager != NULL) {
        IDs[0] = H -> UniqueGlblCtlPtIDMax;
        IDs[1] = TrivIGAGlobalDataManager -> NextTrivariateID - 1;
        IDs[2] = TrivIGAGlobalDataManager -> NextArrangementID - 1;
	return IDs;
    }
    else {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the range of control point IDs used by this given TV.            M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.		     M
*   TV:           The Trivar to fetch its ID's range.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:     Vector of maximal IDs, (Max CtlPts ID, Max Trivars ID) if     M
*              successful, NULL otherwise.  Allocated statically.            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetCtlPtIDRange                                                   M
*****************************************************************************/
int *TrivIGAGetCtlPtIDRange(TrivIGAArrangementID ArgmntID,
			    const TrivTVStruct *TV)
{
    IRIT_STATIC_DATA int
        IDs[2];
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H != NULL) {
        const TrivIGATVStruct
	    *T = TrivIGATV2IGATV(H, TV);

	if (T != NULL) {
	    IDs[0] = T -> UniqueCtlPtIDMin;
	    IDs[1] = T -> UniqueCtlPtIDMax;
	    return IDs;
	}
	else {
	    H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	    return NULL;
	}
    }
    else {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the number of Bezier trivariates that exists along the three     M
* axes of (possibly B-spline) trivar TV in the IGA arrangement.              M
*   Assumes TV has open end conditions.                                      M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.		     M
*   TV:    The Trivar to fetch how many Bezier trivariates it has.	     M
*   NumU:  Number of Bezier trivariates in U.                                M
*   NumV:  Number of Bezier trivariates in V.                                M
*   NumW:  Number of Bezier trivariates in W.                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetNumBzrElements                                                 M
*****************************************************************************/
int TrivIGAGetNumBzrElements(TrivIGAArrangementID ArgmntID,
			     const TrivTVStruct *TV,
			     int *NumU,
			     int *NumV,
			     int *NumW)
{
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct
        *T = TrivIGATV2IGATV(H, TV);
    CagdRType *KnotVals;
    int Len, *KnotMults;

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return FALSE;
    }

    /* We need to count the number of unique interior knots in each axis. */
    Len = IRIT_MAX(IRIT_MAX(TV -> UOrder, TV -> VOrder), TV -> WOrder) +
          IRIT_MAX(IRIT_MAX(TV -> ULength, TV -> VLength), TV -> WLength);
    KnotVals = (CagdRType *) IritMalloc(sizeof(CagdRType) * Len);
    KnotMults = (int *) IritMalloc(sizeof(int) * Len);

    Len = BspKnotsMultiplicityVector(&TV -> UKnotVector[TV -> UOrder - 1],
				     TV -> ULength - TV -> UOrder + 2,
				     KnotVals, KnotMults);
    *NumU = Len - 1;               /* N unique values bound N-1 intervals. */

    Len = BspKnotsMultiplicityVector(&TV -> VKnotVector[TV -> VOrder - 1],
				     TV -> VLength - TV -> VOrder + 2,
				     KnotVals, KnotMults);
    *NumV = Len - 1;               /* N unique values bound N-1 intervals. */

    Len = BspKnotsMultiplicityVector(&TV -> WKnotVector[TV -> WOrder - 1],
				     TV -> WLength - TV -> WOrder + 2,
				     KnotVals, KnotMults);
    *NumW = Len - 1;               /* N unique values bound N-1 intervals. */

    IritFree(KnotVals);
    IritFree(KnotMults);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns a dynamically allocated vector of control points of the Bezier   M
* trivariate at the designated indices in TV.                                M
*   Assumes TV has open end conditions.                                      M
*   Control points will be returned in order, U changing first, W last.      M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.		     M
*   TV:     The Trivar to fetch the control points of one of its Beziers.    M
*   IndexU: The Bezier trivariate index in U in TV, starting from zero.      M
*   IndexV: The Bezier trivariate index in V in TV, starting from zero.      M
*   IndexW: The Bezier trivariate index in W in TV, starting from zero.      M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGACtrlPtStruct *:  A vector of control points of the Bezier trivar. M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetBzrElementCtrlPts                                              M
*****************************************************************************/
TrivIGACtrlPtStruct *TrivIGAGetBzrElementCtrlPts(
					 TrivIGAArrangementID ArgmntID,
					 const TrivTVStruct *TV,
					 int IndexU,
					 int IndexV,
					 int IndexW)
{
    int i, j, k, n, UBase, VBase, WBase, IsNotRational, NumOfCoords;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct
        *T = TrivIGATV2IGATV(H, TV);
    TrivIGACtrlPtStruct *CtlPts;

    if (TV == NULL || T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }

    CtlPts = (TrivIGACtrlPtStruct *) IritMalloc(sizeof(TrivIGACtrlPtStruct) *
				  TV -> UOrder * TV -> VOrder * TV -> WOrder);

    IsNotRational = !CAGD_IS_RATIONAL_PT(TV -> PType);
    NumOfCoords = CAGD_NUM_OF_PT_COORD(TV -> PType);

    /* Compute the Base (first) control point index in each direction. */
    UBase = TrivIGAPatchIdx2CtlPtIdx(H, TV -> UKnotVector,
				     TV -> ULength + TV -> UOrder, IndexU);
    VBase = TrivIGAPatchIdx2CtlPtIdx(H, TV -> VKnotVector,
				     TV -> VLength + TV -> VOrder, IndexV);
    WBase = TrivIGAPatchIdx2CtlPtIdx(H, TV -> WKnotVector,
				     TV -> WLength + TV -> WOrder, IndexW);

    if (UBase < 0 || VBase < 0 || WBase < 0) {
	H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	IritFree(CtlPts);
	return NULL;
    }

    /* Copy the control points. */
    for (n = 0, k = WBase; k < WBase + TV -> WOrder; k++) {
        for (j = VBase; j < VBase + TV -> VOrder; j++) {
	    for (i = UBase; i < UBase + TV -> UOrder; i++, n++) {
	        int l,
		    m = TRIV_MESH_UVW(TV, i, j, k);

	        CtlPts[n].PtType = TV -> PType;
		for (l = IsNotRational; l <= NumOfCoords; l++) {
		    CtlPts[n].Coord[l] = TV -> Points[l][m];
		}
		CtlPts[n].ID = T -> CtlPtsIDs[m].ID;
 	    }
	}
    }

    return CtlPts;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the knot sequence of TV that is used in the designated Bezier    M
* trivar index and direction.                                                M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.		     M
*   TV:     The Trivar to fetch its knot sequence.                           M
*   Dir:    The direction to fetch the knots: U, V, or W.		     M
*   BzrIntervalIndex:   Index of Bezier trivar, starting from zero.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   const CagdRType *:   The list of the (2*Order) knots, or NULL if error.  M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetKnotInterval                                                   M
*****************************************************************************/
const CagdRType *TrivIGAGetKnotInterval(TrivIGAArrangementID ArgmntID,
					const TrivTVStruct *TV,
					TrivTVDirType Dir,
					int BzrIntervalIndex)
{
    int Idx;
    const CagdRType *RetVal;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    switch (Dir) {
        case TRIV_CONST_U_DIR:
	    Idx = TrivIGAPatchIdx2CtlPtIdx(H, TV -> UKnotVector,
					   TV -> ULength + TV -> UOrder,
					   BzrIntervalIndex);
	    RetVal = Idx >= 0 ? &TV -> UKnotVector[Idx] : NULL;
	    break;
        case TRIV_CONST_V_DIR:
	    Idx = TrivIGAPatchIdx2CtlPtIdx(H, TV -> VKnotVector,
					   TV -> VLength + TV -> VOrder,
					   BzrIntervalIndex);
	    RetVal = Idx >= 0 ? &TV -> VKnotVector[Idx] : NULL;
	    break;
        case TRIV_CONST_W_DIR:
	    Idx = TrivIGAPatchIdx2CtlPtIdx(H, TV -> WKnotVector,
					   TV -> WLength + TV -> WOrder,
					   BzrIntervalIndex);
	    RetVal = Idx >= 0 ? &TV -> WKnotVector[Idx] : NULL;
	    break;
        default:
	    assert(0);
	    H -> LastError = TRIV_IGA_ERR_INVALID_DIR;
	    RetVal = NULL;
    }

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   A modifier function to update TV, which must be a geometry TV, with one  *
* delta value identified by ID.                                              *
*                                                                            *
* PARAMETERS:                                                                *
*   H:        A handle on the IGA arrangement to process.		     *
*   ID:       Unique ID of the value to update.                              *
*   PType:    Point type pf point to update.				     *
*   Coords:   Coefficients of point to update.                               *
*   IsDelta:  TRUE if Val is a delta value, FALSE a new actual value.        *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    TRUE if successful, FALSE otherwise.                             *
*****************************************************************************/
static int TrivIGAUpdateOneCtlPt(TrivIGAArrangementStruct *H,
				 int ID,
				 CagdPointType PType,
				 const CagdRType *Coords,
				 int IsDelta)
{
    int IsRational = CAGD_IS_RATIONAL_PT(PType),
        NumOfCoords = CAGD_NUM_OF_PT_COORD(PType);
    TrivIGATVStruct *IGATV;
    TrivIGAFieldStruct *F;

    if (H == NULL)
	return FALSE;

    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (IGATV = F -> TVs; IGATV != NULL; IGATV = IGATV -> Pnext) {
	    TrivTVStruct
		*TV = IGATV -> TV;
	    int i, k,
		Modified = FALSE,
		Len = TRIV_TV_UPT_LST_LEN(TV) *
		      TRIV_TV_VPT_LST_LEN(TV) *
	              TRIV_TV_WPT_LST_LEN(TV);
	    const TrivIGACtlPtUniqueIDsStruct
		*CtlPtsIDs = IGATV -> CtlPtsIDs;

	    for (i = 0; i < Len; i++) {
		if (CtlPtsIDs[i].ID == ID) {
		    /* Found this ID! */
		    if (TV -> PType != PType) {
		        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
		        return FALSE;
		    }
		  
		    for (k = !IsRational; k <= NumOfCoords; k++) {
			if (IsDelta)
			    TV -> Points[k][i] += Coords[k];
			else
			    TV -> Points[k][i] = Coords[k];
			Modified = TRUE;
		    }
		}
	    }

	    if (Modified) {
		/* Update all derivatives. */
		TrivIGAFreeIGATVSlots(IGATV, FALSE);
		TrivIGAUpdateIGATVSlots(H, IGATV, TV, TRUE);
	    }
	}
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A modifier function to update control points in an arrangement, with     M
* the delta values as specified by DeltaVals.  The control points are        M
* identified by their unique IDs and must all be in some TV(s).              M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.		     M
*   NumCtrlPts:   Length of vector DeltaVals.                 	             M
*   DeltaVals:    Vector of control points with delta values to set the      M
*                 relevant control points of TV with.			     M
*		    The control point type in Vals prescribes how many       M
*                 (Coord[i], ID[i]) pairs are in Vals which can be	     M
*                 arbitrary, as a regular int value.			     M
*		    The values are set by identified the ID in TV, for each  M
*                 (Coord[i], ID[i]) pair.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAUpdateCtrlPtsPositions                                            M
*****************************************************************************/
int TrivIGAUpdateCtrlPtsPositions(TrivIGAArrangementID ArgmntID,
				  int NumCtrlPts,
				  const TrivIGACtrlPtStruct *DeltaVals)
{
    int i;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL)
        return FALSE;

    for (i = 0; i < NumCtrlPts; i++) {
        const TrivIGACtrlPtStruct
	    *D = &DeltaVals[i];
    
	if (!TrivIGAUpdateOneCtlPt(H, D -> ID, D -> PtType, D -> Coord, TRUE))
	    return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A modifier function to replace control points' values in an arrangement, M
* with the values as specified by Vals.  The control points are              M
* identified by their unique IDs and must all be in some TV(s).              M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.		     M
*   NumCtrlPts:   Length of vector Vals.                 	             M
*   Vals:         Vector of control points to set the relevant control       M
*                 points of TV with.					     M
*		    The control point type in Vals prescribes how many       M
*                 (Coord[i], ID[i]) pairs are in Vals which can be	     M
*                 arbitrary, as a regular int value.			     M
*		    The values are set by identified the ID in TV, for each  M
*                 (Coord[i], ID[i]) pair.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGASetCtrlPtsPositions                                               M
*****************************************************************************/
int TrivIGASetCtrlPtsPositions(TrivIGAArrangementID ArgmntID,
			       int NumCtrlPts,
			       const TrivIGACtrlPtStruct *Vals)
{
    int i;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL)
        return FALSE;

    for (i = 0; i < NumCtrlPts; i++) {
        const TrivIGACtrlPtStruct
	    *D = &Vals[i];

	if (!TrivIGAUpdateOneCtlPt(H, D -> ID, D ->PtType, D -> Coord, FALSE))
	    return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Evaluate the given TV at the given location.  The evaluation location    M
* is identified by the relevant Bezier trivar (using IndexU/V/W) and the     M
* parameter values that are always normalized to be between zero and one.    M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:  A handle on the IGA arrangement to process.		     M
*   TV:        Trivar to evaluation a one location.                          M
*   EvalType:  What to evaluate:  position, derivatives, normals, etc.       M
*   IndexU:    The Bezier trivariate index in U in TV, starting from zero.   M
*   IndexV:    The Bezier trivariate index in V in TV, starting from zero.   M
*   IndexW:    The Bezier trivariate index in W in TV, starting from zero.   M
*   U, V, W:   Parameters values to evaluate the designated Bezier trivar at.M
*              In [0, 1] only.						     M
*                For many evaluations, U varying first and W last will yield M
*              the best performance.                                         M
*                                                                            *
* RETURN VALUE:                                                              M
*   const TrivIGACtrlPtStruct *:  A control point allocated statically with  M
*                                 the evaluation values, or NULL if error.   M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATVEval                                                            M
*****************************************************************************/
const TrivIGACtrlPtStruct *TrivIGATVEval(TrivIGAArrangementID ArgmntID,
					 const TrivTVStruct *TV,
					 TrivIGAEvalType EvalType,
					 int IndexU,
					 int IndexV,
					 int IndexW,
					 CagdRType U,
					 CagdRType V,
					 CagdRType W)
{
    IRIT_STATIC_DATA TrivIGACtrlPtStruct CtlPt;
    CagdRType *R;
    TrivTVStruct *TVEval,
        *TVEval2 = NULL,
        *TVEval3 = NULL,
        *TVEval4 = NULL;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct
        *T = TrivIGATV2IGATV(H, TV);

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }

    if (!TrivIGAUVW2RealUVW(H, TV, IndexU, IndexV, IndexW, &U, &V, &W))
        return NULL;

    switch (EvalType) {
        case TRIV_IGA_EVAL_VALUE:
	    TVEval = T -> TV;
	    break;
	case TRIV_IGA_EVAL_DU:
	    if (TRIV_IS_RATIONAL_TV(TV))
	        TVEval2 = T -> TV;
	    TVEval = T -> DuTV;
	    break;
	case TRIV_IGA_EVAL_DV:
	    if (TRIV_IS_RATIONAL_TV(TV))
	        TVEval2 = T -> TV;
	    TVEval = T -> DvTV;
	    break;
	case TRIV_IGA_EVAL_DW:
	    if (TRIV_IS_RATIONAL_TV(TV))
	        TVEval2 = T -> TV;
	    TVEval = T -> DwTV;
	    break;
	case TRIV_IGA_EVAL_D2U:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        TVEval2 = T -> DuTV;
	        TVEval3 = T -> DuTV;
	        TVEval4 = T -> TV;
	    }
	    TVEval = T -> Du2TV;
	    break;
	case TRIV_IGA_EVAL_D2V:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        TVEval2 = T -> DvTV;
	        TVEval3 = T -> DvTV;
	        TVEval4 = T -> TV;
	    }
	    TVEval = T -> Dv2TV;
	    break;
	case TRIV_IGA_EVAL_D2W:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        TVEval2 = T -> DwTV;
	        TVEval3 = T -> DwTV;
	        TVEval4 = T -> TV;
	    }
	    TVEval = T -> Dw2TV;
	    break;
	case TRIV_IGA_EVAL_DUDV:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        TVEval2 = T -> DuTV;
	        TVEval3 = T -> DvTV;
	        TVEval4 = T -> TV;
	    }
	    TVEval = T -> DuDvTV;
	    break;
	case TRIV_IGA_EVAL_DUDW:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        TVEval2 = T -> DuTV;
	        TVEval3 = T -> DwTV;
	        TVEval4 = T -> TV;
	    }
	    TVEval = T -> DuDwTV;
	    break;
	case TRIV_IGA_EVAL_DVDW:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        TVEval2 = T -> DvTV;
	        TVEval3 = T -> DwTV;
	        TVEval4 = T -> TV;
	    }
	    TVEval = T -> DvDwTV;
	    break;
	case TRIV_IGA_EVAL_BNRML:
	    if (TrivIGATVEvalBndryNormal(H, TV, U, V, W, &CtlPt))
	        return &CtlPt;
	    else
	        return NULL;
        default:
	    assert(0);
	    H -> LastError = TRIV_IGA_ERR_INVALID_EVAL_TYPE;
	    return NULL;
    }

    CtlPt.PtType = TVEval -> PType;
    R = TrivTVEval(TVEval, U, V, W);

    switch (EvalType) {
        case TRIV_IGA_EVAL_VALUE:
 	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        int l,
		    NumOfCoords = CAGD_NUM_OF_PT_COORD(TV -> PType);

		/* Coerce coordinate to Euclidean space. */
		for (l = 1; l <= NumOfCoords; l++) {
		    CtlPt.Coord[l] = R[l] / R[0];
		}
	    }    
	    else
	        CAGD_GEN_COPY(CtlPt.Coord, R,
			      sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
	    break;
	case TRIV_IGA_EVAL_DU:
	case TRIV_IGA_EVAL_DV:
	case TRIV_IGA_EVAL_DW:
 	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        int l,
		    NumOfCoords = CAGD_NUM_OF_PT_COORD(TV -> PType);
		CagdRType DR[CAGD_MAX_PT_SIZE];

		/* Keep 1st derivatives in DR and evaluate func values in R. */
		CAGD_GEN_COPY(DR, R, sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
		assert(TVEval2 != NULL);
		R = TrivTVEval(TVEval2, U, V, W);

		/* Evaluate the rational derivative using quotient rule. */
		for (l = 1; l <= NumOfCoords; l++) {
		    CtlPt.Coord[l] = (DR[l] * R[0] -
				      DR[0] * R[l]) / IRIT_SQR(R[0]) ;
		}
	    }    
	    else
	        CAGD_GEN_COPY(CtlPt.Coord, R,
			      sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
	    break;

	case TRIV_IGA_EVAL_D2U:
	case TRIV_IGA_EVAL_D2V:
	case TRIV_IGA_EVAL_D2W:
	case TRIV_IGA_EVAL_DUDV:
	case TRIV_IGA_EVAL_DUDW:
	case TRIV_IGA_EVAL_DVDW:
	    if (TRIV_IS_RATIONAL_TV(TV)) {
	        int l,
		    NumOfCoords = CAGD_NUM_OF_PT_COORD(TV -> PType);
		CagdRType DuR[CAGD_MAX_PT_SIZE], DvR[CAGD_MAX_PT_SIZE],
		    DuDvR[CAGD_MAX_PT_SIZE];

		/* Keep 2nd derivatives in DuDvR and evaluate 1st           */
		/*  derivatives and func values into R, DuR, and DvR.       */
		CAGD_GEN_COPY(DuDvR, R, sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
		assert(TVEval2 != NULL && TVEval3 != NULL && TVEval4 != NULL);
		R = TrivTVEval(TVEval2, U, V, W);
		CAGD_GEN_COPY(DuR, R, sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
		R = TrivTVEval(TVEval3, U, V, W);
		CAGD_GEN_COPY(DvR, R, sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
		R = TrivTVEval(TVEval4, U, V, W);

		/* Evaluate 2nd rational derivative using quotient rule. */
		for (l = 1; l <= NumOfCoords; l++) {
		    CtlPt.Coord[l] =
		        (DuDvR[l] * IRIT_SQR(R[0]) -
			 DuR[l] * DvR[0] * R[0] -
			 DvR[l] * DuR[0] * R[0] -
			 R[l] * (DuDvR[0] * R[0] - 2 * DuR[0] * DvR[0])) /
							IRIT_CUBE(R[0]);
		}
		break;
	    }
	    CAGD_GEN_COPY(CtlPt.Coord, R,
			  sizeof(CagdRType) * CAGD_MAX_PT_SIZE);
	    break;

        default:
	    assert(0);
	    H -> LastError = TRIV_IGA_ERR_INVALID_EVAL_TYPE;
	    return NULL;
    }

    return &CtlPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Evaluate the (non unit) normal of the TV at some boundary.               *
*                                                                            *
* PARAMETERS:                                                                *
*   H:         A handle on the IGA arrangement to process.		     *
*   TV:        Trivar to evaluation a one location.                          *
*   U, V, W:   Parameter values to evaluate normal at.  Must be on boundary. *
*   CtlPt:     Unnormalized normal vector will be kept here.                 *
*                                                                            *
* RETURN VALUE:                                                              *
*   TrivIGACtrlPtStruct *: CtlPt with normal if successful, FALSE otherwise. *
*****************************************************************************/
static TrivIGACtrlPtStruct *TrivIGATVEvalBndryNormal(
					 TrivIGAArrangementStruct *H,
					 const TrivTVStruct *TV,
					 CagdRType U,
					 CagdRType V,
					 CagdRType W,
					 TrivIGACtrlPtStruct *CtlPt)
{
    int IdxU, IdxV, IdxW;
    CagdRType *R;
    CagdVType V1, V2;
    const TrivIGATVStruct
        *T = TrivIGATV2IGATV(H, TV);
    TrivTVStruct *TV1, *TV2;

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }
  
    IdxU = BspKnotLastIndexLE(TV -> UKnotVector,
			      TV -> ULength + TV -> UOrder, U);
    IdxV = BspKnotLastIndexLE(TV -> VKnotVector,
			      TV -> VLength + TV -> VOrder, V);
    IdxW = BspKnotLastIndexLE(TV -> WKnotVector,
			      TV -> WLength + TV -> WOrder, W);
    if (IdxU < 0 || IdxV < 0 || IdxW < 0) {
        H -> LastError = TRIV_IGA_ERR_NOT_BNDRY_PARAMS;
	return FALSE;
    }

    if (IRIT_APX_EQ(TV -> UKnotVector[IdxU], U)) {
        TV1 = T -> DvTV;
        TV2 = T -> DwTV;
    }
    else if (IRIT_APX_EQ(TV -> VKnotVector[IdxV], V)) {
        TV1 = T -> DuTV;
        TV2 = T -> DwTV;
    }
    else if (IRIT_APX_EQ(TV -> WKnotVector[IdxW], W)) {
        TV1 = T -> DuTV;
        TV2 = T -> DvTV;
    }
    else {
        H -> LastError = TRIV_IGA_ERR_NOT_BNDRY_PARAMS;
	return NULL;
    }

    /* Evaluate the two partials and compute the normal as the cross prod. */
    R = TrivTVEval(TV1, U, V, W);
    CagdCoerceToE3(V1, &R, -1, TV1 -> PType);
    R = TrivTVEval(TV2, U, V, W);
    CagdCoerceToE3(V2, &R, -1, TV2 -> PType);

    IRIT_CROSS_PROD(&CtlPt -> Coord[1], V1, V2);
    CtlPt -> PtType = CAGD_PT_E3_TYPE;

    return CtlPt;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Evaluate the basis function values of given TV at the given location.    M
*   The evaluation location is identified by the relevant Bezier trivar      M
* (using IndexU/V/W) and the parameter values that are always normalized to  M
* be between zero and one.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:  A handle on the IGA arrangement to process.		     M
*   TV:        Trivar to evaluation a one location.                          M
*   EvalType:  What to evaluate:  position, derivatives, normals, etc.       M
*   Dir:       Direction to evaluate the basic function - U, V or W.         M
*   Index:     The Bezier trivariate index in Dir in TV, starting from zero. M
*   t:         Parameter values to evaluate the designated Bezier trivar at, M
*              in Dir, in [0, 1].					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   const CagdRType *:  Allocated statically vector of Order values holding  M
*                      the relevant basis function values, or NULL if error. M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATVEvalBasis                                                       M
*****************************************************************************/
const CagdRType *TrivIGATVEvalBasis(TrivIGAArrangementID ArgmntID,
				    const TrivTVStruct *TV,
				    TrivIGAEvalType EvalType,
				    TrivTVDirType Dir,
				    int Index,
				    CagdRType t)
{
    int Idx1st, Length, Order;
    CagdRType *KV, *RetVal;
    TrivTVStruct *TVEval;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct
        *T = TrivIGATV2IGATV(H, TV);

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }

    switch (EvalType) {
        case TRIV_IGA_EVAL_VALUE:
	    TVEval = T -> TV;
	    break;
	case TRIV_IGA_EVAL_DU:
	    TVEval = T -> DuTV;
	    break;
	case TRIV_IGA_EVAL_DV:
	    TVEval = T -> DvTV;
	    break;
	case TRIV_IGA_EVAL_DW:
	    TVEval = T -> DwTV;
	    break;
	case TRIV_IGA_EVAL_D2U:
	    TVEval = T -> Du2TV;
	    break;
	case TRIV_IGA_EVAL_D2V:
	    TVEval = T -> Dv2TV;
	    break;
	case TRIV_IGA_EVAL_D2W:
	    TVEval = T -> Dw2TV;
	    break;
	case TRIV_IGA_EVAL_DUDV:
	    TVEval = T -> DuDvTV;
	    break;
	case TRIV_IGA_EVAL_DUDW:
	    TVEval = T -> DuDwTV;
	    break;
	case TRIV_IGA_EVAL_DVDW:
	    TVEval = T -> DvDwTV;
	    break;
	case TRIV_IGA_EVAL_BNRML:
	    fprintf(stderr, "Normal boundary evaluation is not supported yet\n");
	    assert(0);
	    return NULL;
        default:
	    H -> LastError = TRIV_IGA_ERR_INVALID_EVAL_TYPE;
	    assert(0);
	    return NULL;
    }

    switch (Dir) {
        case TRIV_CONST_U_DIR:
	    KV = TVEval -> UKnotVector;
	    Length = TVEval -> ULength;
	    Order = TVEval -> UOrder;
	    break;
        case TRIV_CONST_V_DIR:
	    KV = TVEval -> VKnotVector;
	    Length = TVEval -> VLength;
	    Order = TVEval -> VOrder;
	    break;
        case TRIV_CONST_W_DIR:
	    KV = TVEval -> WKnotVector;
	    Length = TVEval -> WLength;
	    Order = TVEval -> WOrder;
	    break;
        default:
	    assert(0);
	    H -> LastError = TRIV_IGA_ERR_INVALID_DIR;
	    return NULL;
    }

    if (!TrivIGAParam2RealParam(H, KV, Length + Order, Order, Index, &t))
        return NULL;

    RetVal = BspCrvCoxDeBoorBasis(KV, Order, Length,
				  FALSE, t, &Idx1st);

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Examine if the two given surfaces are the same.  Allow the following:    *
* 1. U can be reversed in one surface, compared to the other.                *
* 2. V can be reversed in one surface, compared to the other.                *
* 3. U and V can be exchanged.						     *
*                                                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   Srf, AdjSrf:  The two surfaces to compare.                               *
*   AdjTV:        The other TV, considered as adjacent.			     *
*   AdjBndryIdx:  Index of adjacent boundary candidate.  Between 0 and 5.    *
*   AdjInfo:      Structure to update with the adjacency info if neighbors.  *
*   Error:        Will record errors, if any.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:    TRUE if same surfaces, FALSE otherwise.                    *
*****************************************************************************/
static CagdBType TrivIGASameSrfs(const CagdSrfStruct *Srf,
				 const CagdSrfStruct *AdjSrf,
				 const TrivIGATVStruct *AdjTV,
				 int AdjBndryIdx,
				 TrivIGAAdjacencyInfoStruct *AdjInfo,
				 TrivIGAErrorType *Error)
{
    CagdBType SameSpc,
        Found = FALSE;
    int Modified = 0;
    CagdSrfStruct *TSrf;

    *Error = TRIV_IGA_ERR_NO_ERROR;

    AdjInfo -> ReverseU =
        AdjInfo -> ReverseV =
            AdjInfo -> ReverseUwithV = FALSE;

    if (CagdSrfsSame(Srf, AdjSrf, TRIV_IGA_NEIGHBORHOOD_EPS)) {
        Found = SameSpc = TRUE;
	AdjInfo -> SameSpace = SameSpc;
    }

    if (!Found) {			            /* Try to reverse in U. */
        TSrf = CagdSrfReverseDir(AdjSrf, CAGD_CONST_U_DIR);

	if ((SameSpc = CagdSrfsSame(Srf, TSrf, TRIV_IGA_NEIGHBORHOOD_EPS)) ||
	    CagdSrfsSame2(Srf, TSrf, TRIV_IGA_NEIGHBORHOOD_EPS, &Modified)) {
	    Found = TRUE;
	    AdjInfo -> ReverseU = TRUE;
	    AdjInfo -> SameSpace = SameSpc;
	}

	CagdSrfFree(TSrf);
    }

    if (!Found) {			            /* Try to reverse in V. */
        TSrf = CagdSrfReverseDir(AdjSrf, CAGD_CONST_V_DIR);

	if ((SameSpc = CagdSrfsSame(Srf, TSrf, TRIV_IGA_NEIGHBORHOOD_EPS)) ||
	    CagdSrfsSame2(Srf, TSrf, TRIV_IGA_NEIGHBORHOOD_EPS, &Modified)) {
	    Found = TRUE;
	    AdjInfo -> ReverseV = TRUE;
	    AdjInfo -> SameSpace = SameSpc;
	}

	CagdSrfFree(TSrf);
    }

    if (!Found) {			        /* Try to reverse U with V. */
        TSrf = CagdSrfReverse2(AdjSrf);

	if ((SameSpc = CagdSrfsSame(Srf, TSrf, TRIV_IGA_NEIGHBORHOOD_EPS)) ||
	    CagdSrfsSame2(Srf, TSrf, TRIV_IGA_NEIGHBORHOOD_EPS, &Modified)) {
	    Found = TRUE;
	    AdjInfo -> ReverseUwithV = TRUE;
	    AdjInfo -> SameSpace = SameSpc;
	}

	CagdSrfFree(TSrf);
    }

    if (Found) {
        if (Modified == 3) {
	    /* Both surfaces same but one is not a sub-space of the other. */
	    *Error = TRIV_IGA_ERR_SRFS_NOT_SUBSPACES;
	    Found = FALSE;
	}

	AdjInfo -> AdjTV = AdjTV -> TV;
	AdjInfo -> AdjBndryIdx = AdjBndryIdx;
    }

     return Found;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given a TV in some field, returns the (up to) six facial neighboring     M
* trivariates of TV, if exists.  Notes:					     M
* 1. The adjacent TV can be the same as TV if TV is closed in some dir(s).   M
* 2. The returned adjacency structure will also state if the adjacent TV     M
*    needs to be reversed in U or in V or U and V should be reversed so      M
*    the boundary can be detected as adjacent.				     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID: A handle on the IGA arrangement to process.		     M
*   TV:       Trivar to seek its (up to six) neighbors.                      M
*   AdjInfo:  Already allocated vector of six entries to be updated herein.  M
*             Returns neighboring info for the six face-boundary surfaces of M
*             TV, if any.  Six neighbors are updated in the following order, M
*             Umin, Umax, Vmin, Vmax, Wmin, Wmax.		             M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAGetVrtxNeighboringTVs, TrivIGAGetEdgeNeighboringTVs               M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetFaceNeighboringTVs                                             M
*****************************************************************************/
int TrivIGAGetFaceNeighboringTVs(TrivIGAArrangementID ArgmntID,
				 const TrivTVStruct *TV,
				 TrivIGAAdjacencyInfoStruct *AdjInfo)
{
    int i, j;
    const TrivIGAFieldStruct *F;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct *TV2,
	*T = TrivIGATV2IGATV(H, TV);
    CagdSrfStruct *BndrySrfs[6], **BndrySrfs2;

    for (i = 0; i < 6; i++)
        IRIT_ZAP_MEM(&AdjInfo[i], sizeof(TrivIGAAdjacencyInfoStruct));

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return FALSE;
    }

    F = T -> Field;
    assert(F != NULL);

    /* Six faces in order: UMin, UMax, VMin, VMax, WMin, WMax. */
    IRIT_GEN_COPY(BndrySrfs, TrivBndrySrfsFromTV(TV),
		  sizeof(CagdSrfStruct *) * 6);

    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (TV2 = F -> TVs; TV2 != NULL; TV2 = TV2 -> Pnext) {	    
	    BndrySrfs2 = TrivBndrySrfsFromTV(TV2 -> TV);

	    for (i = 0; i < 6; i++) {     /* For the six faces of input TV. */
		if (AdjInfo[i].AdjTV != NULL)
		    continue;/* Already found adjacent TV for this boundary.*/

		for (j = 0; j < 6; j++) { /* Compare six boundaries of TV2. */
		    if (TV == TV2 -> TV && i == j)
			continue; /* Do not compare to the exact same face. */
		    if (TrivIGASameSrfs(BndrySrfs[i], BndrySrfs2[j], TV2,
					j, &AdjInfo[i],
					&H -> LastError)) { /* Also update! */
			break;
		    }
		}
	    }

	    for (i = 0; i < 6; i++)
		CagdSrfFree(BndrySrfs2[i]);
	}
    }
    
    for (i = 0; i < 6; i++)
        CagdSrfFree(BndrySrfs[i]);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given a TV in some field, returns the edge-neighboring trivariates of    M
* TV, if exists.							     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID: A handle on the IGA arrangement to process.		     M
*   TV:       Trivar to seek its edge neighbors.	                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:     A TRIV_IGA_INVALID_TV_ID terminated vector of edge	     M
*	       neighboring TV IDS or NULL if error.		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAGetFaceNeighboringTVs, TrivIGAGetVrtxNeighboringTVs               M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetEdgeNeighboringTVs                                             M
*****************************************************************************/
int *TrivIGAGetEdgeNeighboringTVs(TrivIGAArrangementID ArgmntID,
				  const TrivTVStruct *TV)
{
    int i, j, MaxNumNeighbors, *Neighbors, *RetNeighbors;
    const TrivIGAFieldStruct *F;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct *TV2,
	*T = TrivIGATV2IGATV(H, TV);
    CagdCrvStruct *BndryEdges, *BndryEdges2;

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return FALSE;
    }

    F = T -> Field;
    assert(F != NULL);

    MaxNumNeighbors = TrivIGAGlobalDataManager -> NextTrivariateID;
    Neighbors = (int *) IritMalloc(sizeof(int) * MaxNumNeighbors);
    IRIT_ZAP_MEM(Neighbors, sizeof(int) * MaxNumNeighbors);

    /* Fetch the twelve edges of TV. */
    BndryEdges = TrivBndryEdgesFromTV(TV);

    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (TV2 = F -> TVs; TV2 != NULL; TV2 = TV2 -> Pnext) {
	    int Found = FALSE;
	    CagdCrvStruct *BEdge, *BEdge2;

	    if (TV2 -> TV == TV) 
		continue;

	    /* Fetch the twelve edges of TV2. */
	    BndryEdges2 = TrivBndryEdgesFromTV(TV2 -> TV);

	    for (i = 0, BEdge = BndryEdges;
		 BEdge != NULL && !Found;
		 i++, BEdge = BEdge -> Pnext) {/* For 12 edges of input TV. */
		CagdCrvStruct
		    *RvrsedBEdge = CagdCrvReverse(BEdge);

		for (j = 0, BEdge2 = BndryEdges2;
		     BEdge2 != NULL && !Found;
		     j++, BEdge2 = BEdge2 -> Pnext) {/* Cmp 12 edges of TV2.*/
		    CagdCrvStruct
		        *PNext = BEdge -> Pnext,
		        *PNext2 = BEdge2 -> Pnext;

		    if (TV == TV2 -> TV && i == j)
			continue; /* Do not compare to the exact same edge. */

		    BEdge -> Pnext = BEdge2 -> Pnext = NULL;
		    if (CagdCrvsSame(BEdge, BEdge2,
				     TRIV_IGA_NEIGHBORHOOD_EPS) ||
			CagdCrvsSame(RvrsedBEdge,  BEdge2,
				     TRIV_IGA_NEIGHBORHOOD_EPS)) {
			Found = TRUE;
		    }
		    BEdge -> Pnext = PNext;
		    BEdge2 -> Pnext = PNext2;
		}

		CagdCrvFree(RvrsedBEdge);
	    }

	    if (Found) {
		int ID = (int) TrivIGADataManagerGetTrivID(TV2 -> TV);
    
		if (ID != TRIV_IGA_INVALID_TV_ID)
		    Neighbors[ID] = 1;
	    }

	    CagdCrvFreeList(BndryEdges2);
	}
    }

    CagdCrvFreeList(BndryEdges);

    RetNeighbors = (int *) IritMalloc(sizeof(int) * (MaxNumNeighbors + 1));
    for (i = j = 0; i < MaxNumNeighbors; i++) {
	if (Neighbors[i] != 0)
	    RetNeighbors[j++] = i;
    }
    RetNeighbors[j] = TRIV_IGA_INVALID_TV_ID;

    IritFree(Neighbors);

    return RetNeighbors;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Given a TV in some field, returns the vertex (corner) neighboring        M
* trivariates of TV, if exists.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID: A handle on the IGA arrangement to process.		     M
*   TV:       Trivar to seek its vertex (corner) neighbors.                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:     A TRIV_IGA_INVALID_TV_ID terminated vector of vertex          M
*              neighboring TV IDS or NULL if error.		             M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAGetFaceNeighboringTVs, TrivIGAGetEdgeNeighboringTVs               M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetVrtxNeighboringTVs                                             M
*****************************************************************************/
int *TrivIGAGetVrtxNeighboringTVs(TrivIGAArrangementID ArgmntID,
				  const TrivTVStruct *TV)
{
    int i, j, MaxNumNeighbors, *Neighbors, *RetNeighbors;
    const TrivIGAFieldStruct *F;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct *TV2,
	*T = TrivIGATV2IGATV(H, TV);
    CagdPtStruct *BndryCrnrs, *BndryCrnrs2;

    if (T == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return FALSE;
    }

    F = T -> Field;
    assert(F != NULL);

    MaxNumNeighbors = TrivIGAGlobalDataManager -> NextTrivariateID;
    Neighbors = (int *) IritMalloc(sizeof(int) * MaxNumNeighbors);
    IRIT_ZAP_MEM(Neighbors, sizeof(int) * MaxNumNeighbors);

    /* Fetch the eight corners of TV. */
    BndryCrnrs = TrivBndryCrnrsFromTV(TV);

    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (TV2 = F -> TVs; TV2 != NULL; TV2 = TV2 -> Pnext) {
	    int Found = FALSE;
	    CagdPtStruct *BCrnr, *BCrnr2;

	    if (TV2 -> TV == TV) 
		continue;

	    /* Fetch the eight corners of TV2. */
	    BndryCrnrs2 = TrivBndryCrnrsFromTV(TV2 -> TV);

	    for (i = 0, BCrnr = BndryCrnrs;
		 BCrnr != NULL && !Found;
		 i++, BCrnr = BCrnr -> Pnext) {/* For 12 edges of input TV. */
		for (j = 0, BCrnr2 = BndryCrnrs2;
		     BCrnr2 != NULL && !Found;
		     j++, BCrnr2 = BCrnr2 -> Pnext) {/* Cmp 12 edges of TV2.*/
		    if (TV == TV2 -> TV && i == j)
			continue; /* Do not compare to the exact same edge. */
		    if (IRIT_PT_APX_EQ(BCrnr -> Pt, BCrnr2 -> Pt)) {
			Found = TRUE;
		    }
		}
	    }

	    if (Found) {
		int ID = (int) TrivIGADataManagerGetTrivID(TV2 -> TV);
    
		if (ID != TRIV_IGA_INVALID_TV_ID)
		    Neighbors[ID] = 1;
	    }

	    CagdPtFreeList(BndryCrnrs2);
	}
    }

    CagdPtFreeList(BndryCrnrs);

    RetNeighbors = (int *) IritMalloc(sizeof(int) * MaxNumNeighbors);
    for (i = j = 0; i < MaxNumNeighbors; i++) {
	if (Neighbors[i] != 0)
	    RetNeighbors[j++] = i;
    }
    RetNeighbors[j] = TRIV_IGA_INVALID_TV_ID;

    IritFree(Neighbors);

    return RetNeighbors;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Free all auxiliary data structure allocated in H.                        M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:   A handle on the IGA arrangement to free.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.                            M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAFreeArrangement                                                   M
*****************************************************************************/
int TrivIGAFreeArrangement(TrivIGAArrangementID ArgmntID)
{
    int i,
	isEmpty = TRUE,
	Res = TrivIGADataManagerFreeArrangement(ArgmntID);

    if (TrivIGAGlobalDataManager != NULL) {
	for (i = 0; i < TRIV_IGA_MAX_ARRANGEMENTS; ++i) {
	    if (TrivIGAGlobalDataManager -> Arrangements[i] != NULL) {
		isEmpty = FALSE;
		break;
	    }
	}
	if (isEmpty == TRUE) {
	    IritFree(TrivIGAGlobalDataManager);
	    TrivIGAGlobalDataManager = NULL;
	}
    }

    return Res;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Free all slots of an IGA TV, but the original TV (and possibly the IDs). *
*                                                                            *
* PARAMETERS:                                                                *
*   IGATV:          IGA TV to free its slots,except the original TV.	     *
*   FreeIDs:        TRUE to also free the IDs array.		             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void TrivIGAFreeIGATVSlots(TrivIGATVStruct *IGATV, CagdBType FreeIDs)
{
    if (FreeIDs)
        IritFree(IGATV -> CtlPtsIDs);

    TrivTVFree(IGATV -> DuTV);
    TrivTVFree(IGATV -> DvTV);
    TrivTVFree(IGATV -> DwTV);
    TrivTVFree(IGATV -> Du2TV);
    TrivTVFree(IGATV -> Dv2TV);
    TrivTVFree(IGATV -> Dw2TV);
    TrivTVFree(IGATV -> DuDvTV);
    TrivTVFree(IGATV -> DuDwTV);
    TrivTVFree(IGATV -> DvDwTV);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Deallocate an arrangement.                                               *
*                                                                            *
* PARAMETERS:                                                                *
*   H:          A handle on the IGA arrangement to free.		     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:                                                                     *
*****************************************************************************/
static int TrivIGADeAllocArrangement(TrivIGAArrangementStruct *H)
{
    int i;
    TrivIGAFieldStruct *F, *Ftmp;
    TrivIGABoundaryNodeStruct *BNode, *BNodeTmp; 

    if (H == NULL)
        return FALSE;

    F = H -> Fields;
    while (F != NULL) {
	TrivIGATVStruct *Ttmp, *T;

        IRIT_LIST_POP(Ftmp, F);

	T = Ftmp -> TVs;
	while (T != NULL) {
	    IRIT_LIST_POP(Ttmp, T);

	    TrivIGAFreeIGATVSlots(Ttmp, TRUE);
	    TrivTVFree(Ttmp -> TV);

	    IritFree(Ttmp);
	}
	IritFree(Ftmp);
    }

    for (i = 0; i < H -> NumMaterials; i++)
	IritFree(H -> Materials[i]);

    BNode = H -> BoundaryNodes;
    while (BNode) {
	BNodeTmp = BNode -> Pnext;
	if (BNode -> BoundaryAxisConditions) {
	    IritFree(BNode -> BoundaryAxisConditions);
	}
	IritFree(BNode);
	BNode = BNodeTmp;
    }

    IritFree(H);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns anID describing the given error.                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:  A handle on the IGA arrangement to process.		     M
*   Reset:     TRUE to also reset the last error as a side effect.           M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAErrorType:     An ID with the error type,			     M
*                         with TRIV_IGA_ERR_NO_ERROR if no error.            M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGADescribeError                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetLastError, error handling                                      M
*****************************************************************************/
TrivIGAErrorType TrivIGAGetLastError(TrivIGAArrangementID ArgmntID,
				     int Reset)
{
    TrivIGAErrorType Err;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL)
        return TRIV_IGA_ERR_INVALID_ARRANGMENT;

    Err = H -> LastError;

    if (Reset)
        H -> LastError = TRIV_IGA_ERR_NO_ERROR;

    return Err;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns a sting describing the given error.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   ErrorNum:   Type of the error that was raised.                           M
*                                                                            *
* RETURN VALUE:                                                              M
*   const char *:     A string describing the error type.                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAGetLastError, TrivDescribeError                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADescribeError, error handling                                     M
*****************************************************************************/
const char *TrivIGADescribeError(TrivIGAErrorType ErrorNum)
{
    IRIT_STATIC_DATA TrivIGAErrorStruct ErrMsgs[] =
    {
	{ TRIV_IGA_ERR_NO_ERROR,		IRIT_EXP_STR("No error") },

	{ TRIV_IGA_ERR_INVALID_ARRANGMENT,   	IRIT_EXP_STR("Invalid arrangement") },
	{ TRIV_IGA_ERR_INVALID_INPUT,   	IRIT_EXP_STR("Invalid input given") },
	{ TRIV_IGA_ERR_RATIONAL_NO_SUPPORT,	IRIT_EXP_STR("Rational trivariates are not supported yet") },
	{ TRIV_IGA_ERR_BEZIER_NO_SUPPORT,	IRIT_EXP_STR("Insertion of Bezier trivariates is not supported") },
	{ TRIV_IGA_ERR_INVALID_PATCH_INDEX,	IRIT_EXP_STR("Invalid Bezier patch index") },
	{ TRIV_IGA_ERR_INVALID_DIR,		IRIT_EXP_STR("Invalid direction") },
	{ TRIV_IGA_ERR_ID_NOT_FOUND,		IRIT_EXP_STR("Unique ID not found") },
	{ TRIV_IGA_ERR_INVALID_EVAL_TYPE,	IRIT_EXP_STR("Invalid evaluation type requested") },
	{ TRIV_IGA_ERR_NOT_BNDRY_PARAMS,        IRIT_EXP_STR("Given parameters are not on patch boundary") },
	{ TRIV_IGA_ERR_NOT_ZERO_ONE_PARAM,	IRIT_EXP_STR("Given parameter is not in [0, 1] domain") },
	{ TRIV_IGA_ERR_EXISTING_MATERIAL_ADDED, IRIT_EXP_STR("Existing material has been added") },
	{ TRIV_IGA_ERR_INVALID_FIELD_ATTRIBUTES_FORMAT, IRIT_EXP_STR("Invalid fields attributes") },
	{ TRIV_IGA_ERR_UNDEFINED_FIELD_MATERIAL,IRIT_EXP_STR("Undefined field's mateiral") },
	{ TRIV_IGA_ERR_INCOMPATIBLE_PT_TYPE,	IRIT_EXP_STR("Incompatible point types given") },
	{ TRIV_IGA_ERR_ID_OUT_OF_RANGE,		IRIT_EXP_STR("ID is out of range") },
	{ TRIV_IGA_ERR_SRFS_NOT_COMPATIBLE,	IRIT_EXP_STR("Given surfaces are not compatible") },
	{ TRIV_IGA_ERR_SRFS_NOT_SUBSPACES,      IRIT_EXP_STR("Given two surfaces are not one sub-space of the other") },

	{ TRIV_IGA_ERR_UNDEFINE_ERR,		NULL },
    };
    int i = 0;

    for ( ; ErrMsgs[i].ErrorDesc != NULL; i++)
	if (ErrorNum == ErrMsgs[i].ErrorNum)
	    return ErrMsgs[i].ErrorDesc;

    return IRIT_EXP_STR("Unknown error occurred");
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Adds a boundary node to the arrangement.                                 *
*                                                                            *
* PARAMETERS:                                                                *
*   H:          A handle on the IGA arrangement to process.		     *
*   IGATV:      Trivariate of the face.                                      *
*   TVNodeIdx:  The boundary node Id.                                        *
*   NodeBoundaryType:   IGA boundary condition type.                         *
*   BoundaryAxisConditions:  String represents the relevant axes, eg. "xy".  *
*   Value:      Value of the boundary condition.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void TrivIGAAddBoundaryNode(TrivIGAArrangementStruct *H,
			           TrivIGATVStruct *IGATV,
			           int TVNodeIdx,
			           TrivIGANodeBoundaryType NodeBoundaryType,
			           const char *BoundaryAxisConditions,
			           CagdRType Value)
{
    CagdBType NodeFound;
    TrivIGABoundaryNodeStruct *BNodeTmp,
        *BNode = (TrivIGABoundaryNodeStruct *) 
                                IritMalloc(sizeof(TrivIGABoundaryNodeStruct));

    BNode -> NodeID = IGATV -> CtlPtsIDs[TVNodeIdx].ID;
    BNode -> BoundaryType = NodeBoundaryType;
    BNode -> Value = Value;

    if (BoundaryAxisConditions) {
	BNode -> BoundaryAxisConditions =  (char *) IritMalloc(sizeof(char) *
				        (1 + strlen(BoundaryAxisConditions)));

	strcpy(BNode -> BoundaryAxisConditions, BoundaryAxisConditions);
    }
    else {
	BNode -> BoundaryAxisConditions = NULL;
    }

    /* check if already exists */
    NodeFound = FALSE;
    BNodeTmp = H -> BoundaryNodes;
    while (BNodeTmp) {
	if (BNodeTmp -> NodeID == BNode -> NodeID) {
	    NodeFound = TRUE;
	    break;
	}
	BNodeTmp = BNodeTmp -> Pnext;
    }
    
    if (!NodeFound) {
	IRIT_LIST_PUSH(BNode, H -> BoundaryNodes);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Adds a boundary condition to a face in the arrangement.                  M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.		     M
*   TV:		   Trivariate of the face.                                   M
*   Boundary:      The boundary face on the trivariate (UMin, VMax, etc.).   M
*   NodeBoundary:  IGA boundary condition type.                              M
*   BoundaryAxisConditions: String represents the relevant axs, eg. "xy".    M
*   Value:         Value of the boundary condition.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    TRUE on success and	FALSE on failure.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAAddBoundaryFace2                                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAAddBoundaryFace                                                   M
*****************************************************************************/
int TrivIGAAddBoundaryFace(TrivIGAArrangementID ArgmntID,
			   const TrivTVStruct *TV,
			   TrivTVBndryType Boundary,
			   TrivIGANodeBoundaryType NodeBoundary,
			   const char *BoundaryAxisConditions,
			   CagdRType Value)
{
    TrivIGATVStruct *IGATV;
    int URange[2], VRange[2], WRange[2], i, j, k, Idx;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (H == NULL || TV == NULL)
	return FALSE;

    /* Check if TV exists in the arrangement trivariates */
    IGATV = TrivIGATV2IGATV(H, TV);
    if (IGATV == NULL) {
	H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	return FALSE;
    }

    URange[0] = VRange[0] = WRange[0] = 0;
    URange[1] = TRIV_TV_UPT_LST_LEN(TV);
    VRange[1] = TRIV_TV_VPT_LST_LEN(TV);
    WRange[1] = TRIV_TV_WPT_LST_LEN(TV);

    switch (Boundary) {
	case TRIV_U_MIN_BNDRY:
	    URange[1] = 1;
	    break;
	case TRIV_U_MAX_BNDRY:
	    URange[0] = URange[1] - 1;
	    break;
	case TRIV_V_MIN_BNDRY:
	    VRange[1] = 1;
	    break;
	case TRIV_V_MAX_BNDRY:
	    VRange[0] = VRange[1] - 1;
	    break;
	case TRIV_W_MIN_BNDRY:
	    WRange[1] = 1;
	    break;
	case TRIV_W_MAX_BNDRY:
	    WRange[0] = WRange[1] - 1;
	    break;

	default:
	    assert(0);
	    break;
    }

    for (k = WRange[0]; k < WRange[1]; k++) {
	for (j = VRange[0]; j < VRange[1]; j++) {
	    for (i = URange[0]; i < URange[1]; i++) {
		Idx = TRIV_MESH_UVW(TV, i, j, k);
		TrivIGAAddBoundaryNode(H, IGATV, Idx, NodeBoundary, 
		                       BoundaryAxisConditions, Value);
	    }
	}
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Adds a boundary condition to a face in the arrangement.                  M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.		     M
*   TV:	 	   Trivariate owning the face.  Can be NULL in which case we M
*                  search all TVs in arrangment for the closest face.        M
*   Pt:            A point in space to select the closest boundary.          M
*   NodeBoundary:  IGA boundary condition type.                              M
*   BoundaryAxisConditions: String represents the relevant axes, eg. "xy".   M
*   Value:         Value of the boundary condition.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    TRUE on success and	FALSE on failure.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAAddBoundaryFace2, TrivIGAGetBoundaryFaceByPt                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAAddBoundaryFaceByPt                                               M
*****************************************************************************/
int TrivIGAAddBoundaryFaceByPt(TrivIGAArrangementID ArgmntID,
			       const TrivTVStruct *TV,
			       const CagdPType Pt,
			       TrivIGANodeBoundaryType NodeBoundary,
			       const char *BoundaryAxisConditions,
			       CagdRType Value)
{
    TrivTVBndryType Boundary;
    int *IDs;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (TV == NULL)
	return FALSE;

    IDs = TrivIGAGetBoundaryFaceByPt(ArgmntID, TV, Pt);
    assert(IDs[0] == TrivIGADataManagerGetTrivID(TV));

    switch (IDs[1]) {
        case 0:
	    Boundary = TRIV_U_MIN_BNDRY;
	    break;
        case 1:
	    Boundary = TRIV_U_MAX_BNDRY;
	    break;
        case 2:
	    Boundary = TRIV_V_MIN_BNDRY;
	    break;
        case 3:
	    Boundary = TRIV_V_MAX_BNDRY;
	    break;
        case 4:
	    Boundary = TRIV_W_MIN_BNDRY;
	    break;
        case 5:
	    Boundary = TRIV_W_MAX_BNDRY;
	    break;
        default:
	    H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	    return FALSE;
    }

    return TrivIGAAddBoundaryFace(ArgmntID, TV, Boundary, NodeBoundary,
				  BoundaryAxisConditions, Value);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Find the face closest to Pt in a specific TV.		             *
*                                                                            *
* PARAMETERS:                                                                *
*   TV:		  Trivariate owning the face.  Can be NULL in which case we  *
*                 search all TVs in arrangement for the closest face.        *
*   Pt:           A point in space to select the closest boundary.           *
*   MinDistSqr:   Only faces closer than this so-far-min-distance-square are *
*                 updating here (and updating this distance square as well). *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:    Index (0 to 5) of closest face to Pt in TV.                      *
*           Returns -1 if no face is closer than input MinDistSqr.           *
*****************************************************************************/
static int TrivIGAGetBoundaryFaceByPtAux(const TrivTVStruct *TV,
					 const CagdPType Pt,
					 CagdRType *MinDistSqr)
{
    int i,
        MinDistIdx = -1;
    CagdSrfStruct
        **Srfs = TrivBndrySrfsFromTV(TV);

    /* Get the six boundary faces (surfaces) of the trivariate and seek the */
    /* face that Pt is distance-sqaure closer than MinDistSqr, if any.      */
     for (i = 0; i < 6; i++) {
        CagdRType *R, DistSqr,
	    *UV = MvarDistSrfPoint(Srfs[i], Pt, TRUE, 1e-3, 1e-10);
	CagdPType PtE3;

	if (UV == NULL) {
	    CagdSrfFree(Srfs[i]);
	    continue;
	}

	R = CagdSrfEval(Srfs[i], UV[0], UV[1]);
	CagdCoerceToE3(PtE3, &R, -1, Srfs[i] -> PType);
	DistSqr = IRIT_PT_PT_DIST_SQR(Pt, PtE3);
	if (*MinDistSqr > DistSqr) {
	    *MinDistSqr = DistSqr;
	    MinDistIdx = i;
	}

	CagdSrfFree(Srfs[i]);
    }

    return MinDistIdx;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Find the face closest to Pt in the arrangement or specific TV.           M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.		     M
*   TV:		   Trivariate owning the face.  Can be NULl in which case we M
*                  search all TVs in arrangement.                            M
*   Pt:            A point in space to select the closest boundary.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:    A list of two integers allocated statically, (TVID, FaceID).   M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAAddBoundaryFace2, TrivIGAAddBoundaryFaceByPt                      M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetBoundaryFaceByPt                                               M
*****************************************************************************/
int *TrivIGAGetBoundaryFaceByPt(TrivIGAArrangementID ArgmntID,
				const TrivTVStruct *TV,
				const CagdPType Pt)
{
    IRIT_STATIC_DATA int ReturnedIDs[2];
    int TVID, FaceID,
        MinFaceID = TRIV_IGA_INVALID_FACE_ID,
        MinTVID = TRIV_IGA_INVALID_TV_ID;
    CagdRType
        MinDistSqr = IRIT_INFNTY;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGAFieldStruct *Field;

    if (TV == NULL) {
        for (Field = H -> Fields; Field != NULL; Field = Field -> Pnext) {
	    TrivIGATVStruct *IGATV;

	    for (IGATV = Field -> TVs; IGATV != NULL; IGATV = IGATV -> Pnext) {
	        TVID = TrivIGADataManagerGetTrivID(IGATV -> TV);

		FaceID = TrivIGAGetBoundaryFaceByPtAux(IGATV -> TV, Pt,
						       &MinDistSqr);
		if (FaceID >= 0) {
		    MinFaceID = FaceID;
		    MinTVID = TVID;
		}
	    }
	}
    }
    else {
        FaceID = TrivIGAGetBoundaryFaceByPtAux(TV, Pt, &MinDistSqr);
	if (FaceID >= 0) {
	    MinFaceID = FaceID;
	    MinTVID = TrivIGADataManagerGetTrivID(TV);
	}
    }

    ReturnedIDs[0] = MinTVID;
    ReturnedIDs[1] = MinFaceID;

    return ReturnedIDs;    
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Allocates and initializes IGA arrangements and trivariates data manager. *
*                                                                            *
* PARAMETERS:                                                                *
*   None                                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   TrivIGADataManager: Pointer to the new allocated data manager in success *
*                       and NULL in failure.				     *
*****************************************************************************/
static TrivIGADataManager *TrivIGADataManagerInit()
{
    int i;
    TrivIGADataManager 
	*DM = (TrivIGADataManager *) IritMalloc(sizeof(TrivIGADataManager));
    
    DM -> NextTrivariateID = TRIV_IGA_INVALID_TV_ID + 1;
    DM -> NextArrangementID = TRIV_IGA_INVALID_ARRANGEMENT_ID + 1;
    
    for (i = 0; i < TRIV_IGA_MAX_ARRANGEMENTS; ++i) {
	DM -> Arrangements[i] = NULL;
    }

    for (i = 0; i < TRIV_IGA_MAX_TRIVARIATES; ++i) {
	DM -> IGATrivars[i] = NULL;
    }
    return DM;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates an IGA arrangement and associates it with new ID.              M
*                                                                            *
* PARAMETERS:                                                                M
*   DM: IGA Data manager	    				             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAArrangementID: ID of the new arrangement in success and           M
*                         INVALID_ARRANGEMENT_ID in failure	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerAllocateArrangement                                    M
*****************************************************************************/
TrivIGAArrangementID TrivIGADataManagerAllocateArrangement(
						       TrivIGADataManager *DM)
{
    TrivIGAArrangementID 
	NewArngmntID = TRIV_IGA_INVALID_ARRANGEMENT_ID;

    if (DM == NULL || DM -> NextArrangementID >= TRIV_IGA_MAX_ARRANGEMENTS) {
	return NewArngmntID;
    }
    
    NewArngmntID = DM -> NextArrangementID++;
    DM -> Arrangements[NewArngmntID] = TrivIGAAllocateArrangement();
    return NewArngmntID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Maps between arrangement ID and its actual structure's address.          M
*                                                                            *
* PARAMETERS:                                                                M
*   ArrngmntID:  ID of the arrangement.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAArrangementStruct *: Pointer to the arrangement structure if      M
*                               success and NULL in failure.	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerGetArrangement                                         M
*****************************************************************************/
TrivIGAArrangementStruct *TrivIGADataManagerGetArrangement(
					      TrivIGAArrangementID ArrngmntID)
{
    if (TrivIGAGlobalDataManager == NULL) {
	return NULL;
    }

    if (ArrngmntID < 0 ||
	ArrngmntID >= TrivIGAGlobalDataManager -> NextArrangementID) {
	return NULL;
    }

    return TrivIGAGlobalDataManager -> Arrangements[ArrngmntID];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Maps between arrangement structure address and its ID.                   M
*                                                                            *
* PARAMETERS:                                                                M
*   H:    Pointer to the arrangement structure.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAArrangementID: The arrangement id in success and NULL in failure  M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerGetArrangementID                                       M
*****************************************************************************/
TrivIGAArrangementID TrivIGADataManagerGetArrangementID(
					         TrivIGAArrangementStruct *H)
{
    int i;
    TrivIGAArrangementID 
	ArgmntID = TRIV_IGA_INVALID_ARRANGEMENT_ID;

    if (TrivIGAGlobalDataManager == NULL || H == NULL) {
	return ArgmntID;
    }

    for (i = 0; i < TrivIGAGlobalDataManager -> NextArrangementID; ++i) {
	if (TrivIGAGlobalDataManager -> Arrangements[i] == H) {
	    ArgmntID = i;
	    break;
	}
    }

    return ArgmntID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   DeAllocate given arrangement from the data manager.                      M
*                                                                            *
* PARAMETERS:                                                                M
*   ArrngmntID:  The arrangement ID to deallocate			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     TRUE if successful, FALSE otherwise.			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerFreeArrangement                                        M
*****************************************************************************/
int TrivIGADataManagerFreeArrangement(TrivIGAArrangementID ArrngmntID)
{
    TrivIGAArrangementStruct
	*Argmnt = TrivIGADataManagerGetArrangement(ArrngmntID);
    
    if (TrivIGAGlobalDataManager == NULL || Argmnt == NULL) {
	return FALSE;
    }

    if (TrivIGADeAllocArrangement(Argmnt) == FALSE) {
	return FALSE;
    }

    TrivIGAGlobalDataManager -> Arrangements[ArrngmntID] = NULL;
    return TRUE;
}


/*****************************************************************************
* DESCRIPTION:                                                               M
*   Allocates an ID to a given trivariate in the IGA data manager, if the    M
* trivariate already exists, its associated ID is returned.                  M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.		     M
*   TV:          Pointer to the trivariate.    		    	             M
*   ID:          if non negative, use this ID as the unique ID of TV.        M
*                Otherwise, assign a unique ID here.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:     ID of the trivariate TRUE if successful, and            M
*                    INVALID_IGA_TV_ID otherwise.			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerAddTrivariate                                          M
*****************************************************************************/
TrivIGATVID TrivIGADataManagerAddTrivariate(TrivIGAArrangementID ArgmntID, 
					    TrivTVStruct *TV,
					    int ID)
{
    TrivIGATVID 
	TVID = TRIV_IGA_INVALID_TV_ID;

    if (TrivIGAGlobalDataManager == NULL) {
	return TVID;
    }

    if ((TVID = TrivIGADataManagerGetTrivID(TV)) == TRIV_IGA_INVALID_TV_ID) {
	TrivIGATVStruct
	    *IGATV = TrivIGANewTV(ArgmntID, TV);

	if (IGATV != NULL) {
	    if (ID >= 0) {
		TVID = ID;
		TrivIGAGlobalDataManager -> NextTrivariateID = ID + 1;
	    }
	    else 
		TVID = TrivIGAGlobalDataManager -> NextTrivariateID++;

	    if (TRIV_IGA_MAX_TRIVARIATES <= TVID) {
		TrivIGAArrangementStruct 
		    *H = TrivIGADataManagerGetArrangement(ArgmntID);

	        H -> LastError = TRIV_IGA_ERR_ID_OUT_OF_RANGE;
		return TRIV_IGA_INVALID_MATERIAL_ID;
	    }
	    else
		TrivIGAGlobalDataManager -> IGATrivars[TVID] = IGATV;
	}
    }

    return TVID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns the ID of a given trivariate from IGA data manager.              M
*                                                                            *
* PARAMETERS:                                                                M
*   TV:   Pointer to the trivariate.	     		    	             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:     ID of the trivariate TRUE if successful, and            M
*                    TRIV_IGA_INVALID_TV_ID otherwise.			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerGetTrivID                                              M
*****************************************************************************/
TrivIGATVID TrivIGADataManagerGetTrivID(const TrivTVStruct *TV)
{
    int i;
    TrivIGATVID 
	TVID = TRIV_IGA_INVALID_TV_ID;

    if (TrivIGAGlobalDataManager == NULL) {
	return TVID;
    }

    for (i = 1; i < TrivIGAGlobalDataManager -> NextTrivariateID; ++i) {
	if (TrivIGAGlobalDataManager -> IGATrivars[i] != NULL &&
	    TrivIGAGlobalDataManager -> IGATrivars[i] -> TV == TV) {
	    TVID = i;
	    break;
	}
    }

    return TVID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Maps between trivariate ID and its structure's address in a given IGA    M
* data manager.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   TVID:   ID of the trivariate.     		    	                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:   Pointer to the trivariate structue if successfull and  M
*                     NULL otherwise.			    		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerGetTrivariate                                          M
*****************************************************************************/
TrivTVStruct *TrivIGADataManagerGetTrivariate(TrivIGATVID TVID)
{
    TrivIGATVStruct
        *IGATV = TrivIGADataManagerGetIGATrivariate(TVID);

    if (IGATV == NULL)
	return NULL;
    else
        return IGATV -> TV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Maps between trivariate ID and its structure's address in a given IGA    M
* data manager.								     M
*                                                                            *
* PARAMETERS:                                                                M
*   TVID:   ID of the trivariate.     		    	                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVStruct *: Pointer to IGA trivariate structue if successfull and M
*                      NULL otherwise.			    		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGADataManagerGetIGATrivariate                                       M
*****************************************************************************/
TrivIGATVStruct *TrivIGADataManagerGetIGATrivariate(TrivIGATVID TVID)
{
    if (TrivIGAGlobalDataManager == NULL) {
	return NULL;
    }

    if (TVID <= 0 || TVID >= TrivIGAGlobalDataManager -> NextTrivariateID) {
	return NULL;
    }

    return TrivIGAGlobalDataManager -> IGATrivars[TVID];
}

