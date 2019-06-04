/*****************************************************************************
* Module to support IGA (isogeometric analysis) related processing.          *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber and Fady Massarwi	 	 Ver 1.0, Sep 2013   *
*****************************************************************************/

#include "triv_loc_iga.h"

static TrivTVStruct *TrivIGAApplyDomainAndSeeding(TrivIGAArrangementID ArgmntID,
						  TrivTVStruct *TV);
/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the default seeding of all to-be-entered trivariates into the       M
* arrangement in the given direction (u, V, or W).			     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   Dir:	    U, V, W direction to set the default seeding to.         M
*   Alpha:          Ratio of last interval's length divided by first         M
*                   interval length, in set direction Dir.		     M
*   NumIntervals:   Number of interval to introduce.                         M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    The number of intervals that will be introduced in that axes.    M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGASetDefaultSeeding                                                 M
*****************************************************************************/
int TrivIGASetDefaultSeeding(TrivIGAArrangementID ArgmntID,
			     TrivTVDirType Dir,
			     CagdRType Alpha,
			     int NumIntervals)
{
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGASeedingStateStruct *SeedingState;

    if (H == NULL)
        return 0;

    SeedingState = &H -> SeedingState;

    if (NumIntervals < 0 || Alpha <= 0.0) {
	H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	return 0;
    }

    switch (Dir) {
	case TRIV_CONST_U_DIR:
	    SeedingState -> AlphaVal[0] = Alpha;
	    SeedingState -> NumOfIntervals[0] = NumIntervals;
	    break;
	case TRIV_CONST_V_DIR:
	    SeedingState -> AlphaVal[1] = Alpha;
	    SeedingState -> NumOfIntervals[1] = NumIntervals;
	    break;
	case TRIV_CONST_W_DIR:
	    SeedingState -> AlphaVal[2] = Alpha;
	    SeedingState -> NumOfIntervals[2] = NumIntervals;
	    break;
	default:
	    H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	    return 0;
    }

    return NumIntervals;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the default domain of all to-be-entered trivariates into the        M
* arrangement in the given direction (u, V, or W).			     M
*   If TMax <= TMin, then this option is disabled (default) and the doamin   M
* is kept as it is inserted into the arrangement.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.                M
*   Dir:	  U, V, W direction to set the default domain to.            M
*   Min, Max:     The domain to set.					     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    TRUE if successful, FALSE otherwise.                             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGASetDefaultDomain                                                  M
*****************************************************************************/
int TrivIGASetDefaultDomain(TrivIGAArrangementID ArgmntID,
			    TrivTVDirType Dir,
			    CagdRType Min,
			    CagdRType Max)
{
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGASeedingStateStruct *SeedingState;

    if (H == NULL)
        return FALSE;

    SeedingState = &H -> SeedingState;

    switch (Dir) {
	case TRIV_CONST_U_DIR:
	    SeedingState -> DmnMin[0] = Min > Max ? 0.0 : Min;
	    SeedingState -> DmnMax[0] = Min > Max ? 0.0 : Max;
	    break;
	case TRIV_CONST_V_DIR:
	    SeedingState -> DmnMin[1] = Min > Max ? 0.0 : Min;
	    SeedingState -> DmnMax[1] = Min > Max ? 0.0 : Max;
	    break;
	case TRIV_CONST_W_DIR:
	    SeedingState -> DmnMin[2] = Min > Max ? 0.0 : Min;
	    SeedingState -> DmnMax[2] = Min > Max ? 0.0 : Max;
	    break;
	default:
	    H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	    return FALSE;
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Apply, in place, seeding to the given trivariate, following the global   M
* arrangement seeding specifications.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   TV:	            Trivariate to apply seeding to.  Used in place.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:   Refined TV, after the seeding was applied.	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAApplyDomainAndSeeding                                             M
*****************************************************************************/
static TrivTVStruct *TrivIGAApplyDomainAndSeeding(TrivIGAArrangementID ArgmntID,
						  TrivTVStruct *TV)
{
    int i;
    CagdRType Min[3], Max[3];
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGASeedingStateStruct *SeedingState;
    TrivTVStruct *TmpTV,
        *RetTV = TV;

    if (H == NULL)
        return NULL;

    SeedingState = &H -> SeedingState;

    /* Apply domain change if so desired. */
    for (i = 0; i < 3; i++) {
	if (SeedingState -> DmnMin[i] < SeedingState -> DmnMax[i]) {
	    int Len, Order;
	    CagdRType *KV;

	    switch (i) {
		default:
		    assert(0);
		case 0:
		    KV = RetTV -> UKnotVector;
		    Order = RetTV -> UOrder;
		    Len = RetTV -> ULength;
		    break;
		case 1:
		    KV = RetTV -> VKnotVector;
		    Order = RetTV -> VOrder;
		    Len = RetTV -> VLength;
		    break;
		case 2:
		    KV = RetTV -> WKnotVector;
		    Order = RetTV -> WOrder;
		    Len = RetTV -> WLength;
		    break;
	    }
	    BspKnotAffineTransOrder2(KV, Order, Len + Order,
				     SeedingState -> DmnMin[i],
				     SeedingState -> DmnMax[i]);
	}
    }

    TrivTVDomain(TV, &Min[0], &Max[0], &Min[1], &Max[1], &Min[2], &Max[2]);

    for (i = 0; i < 3; i++) {
        int j, n1;
        CagdRType Step, Inc, Base, Sum, *RefKV;

        if (SeedingState -> NumOfIntervals[i] <= 1)
	    continue;

	n1 = SeedingState -> NumOfIntervals[i] - 1;

	/* Compute the intervals over the [0, 1] domain: */
	Step = Inc = pow(SeedingState -> AlphaVal[i], 1.0 / n1);
	Base = 1.0;

	RefKV = (CagdRType *) IritMalloc(sizeof(CagdRType) * n1);
	RefKV[0] = Base;
	for (j = 1; j < n1; j++) {
	    RefKV[j] = RefKV[j - 1] + Inc;
	    Inc *= Step;
	}
	Sum = RefKV[n1 - 1] + Inc;

	/* Normalize to exact domain. */
	for (j = 0; j < n1; j++)
	    RefKV[j] = Min[i] + (Max[i] - Min[i]) * RefKV[j] / Sum;

	TmpTV = TrivTVRefineAtParams(RetTV, TRIV_INT_TO_DIR(i),
				     FALSE, RefKV, n1);
	TrivTVFree(RetTV);
	RetTV = TmpTV;

	IritFree(RefKV);
    }

    return RetTV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Copy an existing new trivariate and add it to the give arrangement.      M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID: A handle on the IGA arrangement to process. 		     M
*   TV:       An existing trivariate to add (a copy thereof) to arrangement. M
*   ID:       ID to use or -1 to set a new ID.   			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:   INVALID_IGA_TV_ID if error, to TV ID if successful.       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAAddTrivar                                                         M
*****************************************************************************/
TrivIGATVID TrivIGAAddTrivar(TrivIGAArrangementID ArgmntID,
                             const TrivTVStruct *TV,
			     int ID)
{
    TrivTVStruct
        *TVCp = TrivTVCopy(TV);

    return TrivIGADataManagerAddTrivariate(ArgmntID,
				   TrivIGAApplyDomainAndSeeding(ArgmntID, TVCp),
				   ID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a new trivariate and add it to the give arrangement as an     M
* extrusion of the given surface.                                            M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:   A handle on the IGA arrangement to process. 		     M
*   Srf:        Surface to extrude.                                          M
*   Vec:        The extrusion vector.                                        M
*   ID:         ID to use or -1 to set a new ID.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:   INVALID_IGA_TV_ID if error, to TV ID if successful.       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAExtrudeTV                                                         M
*****************************************************************************/
TrivIGATVID TrivIGAExtrudeTV(TrivIGAArrangementID ArgmntID,
                             const CagdSrfStruct *Srf,
                             const IrtVecType Vec,
			     int ID)
{
    CagdVecStruct V;
    TrivTVStruct *TV;

    IRIT_VEC_COPY(V.Vec, Vec);
    TV = TrivExtrudeTV(Srf, &V);

    if (TV == NULL)
	return TRIV_IGA_INVALID_TV_ID;

    if (TRIV_IS_BEZIER_TV(TV)) {
	TrivTVStruct
	    *TV2 = TrivCnvrtBzr2BspTV(TV);

	TrivTVFree(TV);
	TV = TV2;
    }

    return TrivIGADataManagerAddTrivariate(ArgmntID,
				   TrivIGAApplyDomainAndSeeding(ArgmntID, TV),
				   ID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a new trivariate and add it to the give arrangement as an     M
* extrusion of the given surface along a curve.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:   A handle on the IGA arrangement to process. 		     M
*   Srf:        Surface to extrude.                                          M
*   Crv:        The extrusion curve.                                         M
*   ID:         ID to use or -1 to set a new ID.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:   INVALID_IGA_TV_ID if error, to TV ID if successful.       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAExtrudeTV2                                                        M
*****************************************************************************/
TrivIGATVID TrivIGAExtrudeTV2(TrivIGAArrangementID ArgmntID,
			      const CagdSrfStruct *Srf,
			      const CagdCrvStruct *Crv,
			      int ID)
{
    TrivTVStruct *TV;

    TV = TrivExtrudeTV2(Srf, Crv);

    if (TV == NULL)
	return TRIV_IGA_INVALID_TV_ID;

    if (TRIV_IS_BEZIER_TV(TV)) {
	TrivTVStruct
	    *TV2 = TrivCnvrtBzr2BspTV(TV);

	TrivTVFree(TV);
	TV = TV2;
    }

    return TrivIGADataManagerAddTrivariate(ArgmntID,
				   TrivIGAApplyDomainAndSeeding(ArgmntID, TV),
				   ID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a new trivariate and add it to the give arrangement as a      M
* volume of revolution of the given surface.                                 M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:              A handle on the IGA arrangement to process.       M
*   Srf:	           Surface to rotate.				     M
*   AxisPoint, AxisVector: Axis line of rotation's prescription.             M
*   StartAngle:            Starting angle in degs., 0 for a full circle.     M
*   EndAngle:              Terminating angle in degs., 360 for a full circle.M
*   IsRational:            TRUE to construct a rational precise volume of    M
*                          revolution, FALSE for polynomial approximation.   M
*   ID:		           ID to use or -1 to set a new ID.		     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:   INVALID_IGA_TV_ID if error, to TV ID if successful.       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATVofRevol                                                         M
*****************************************************************************/
TrivIGATVID TrivIGATVofRevol(TrivIGAArrangementID ArgmntID,
                             const CagdSrfStruct *Srf,
                             const IrtPtType AxisPoint,
                             const IrtVecType AxisVector,
                             CagdRType StartAngle,
                             CagdRType EndAngle,
                             CagdBType IsRational,
			     int ID)
{
    TrivTVStruct *TV;

    if (!IRIT_APX_EQ(StartAngle, 0.0) || !IRIT_APX_EQ(EndAngle, 360)) {
	if (!IRIT_APX_EQ(AxisPoint[0], 0.0) ||
	    !IRIT_APX_EQ(AxisPoint[1], 0.0) ||
	    !IRIT_APX_EQ(AxisVector[0], 0.0) ||
	    !IRIT_APX_EQ(AxisVector[1], 0.0)) {
	    TrivIGAArrangementStruct 
	        *H = TrivIGADataManagerGetArrangement(ArgmntID);

	    /* Do not support arbitrary angles and arbitrary position/axes. */
	    H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
	    return TRIV_IGA_INVALID_TV_ID;
	}

        TV = TrivTVOfRev2(Srf, !IsRational, StartAngle, EndAngle);
    }
    else
        TV = TrivTVOfRevAxis(Srf, AxisPoint, AxisVector, !IsRational);

    if (TV == NULL)
	return TRIV_IGA_INVALID_TV_ID;

    return TrivIGADataManagerAddTrivariate(ArgmntID,
				   TrivIGAApplyDomainAndSeeding(ArgmntID, TV),
				   ID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a new trivariate and add it to the give arrangement as a      M
* trivariate through surfaces of the given ordered surfaces.                 M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:          A handle on the IGA arrangement to process.           M
*   SrfList:           List of surfaces to construct a trivariate with.      M
*   OtherOrder:        Other, third, order of trivariate. 		     M
*   IsInterpolating:   TRUE for the trivariate to interpolate the surfaces,  M
*                      FALSE to approximate them.                            M
*   ID:		       ID to use or -1 to set a new ID.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:   INVALID_IGA_TV_ID if error, to TV ID if successful.       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATVFromSurfaces                                                    M
*****************************************************************************/
TrivIGATVID TrivIGATVFromSurfaces(TrivIGAArrangementID ArgmntID,
                                  const CagdSrfStruct *SrfList,
                                  int OtherOrder,
                                  CagdBType IsInterpolating,
				  int ID)
{
    TrivTVStruct *TV;

    if (IsInterpolating)
        TV = TrivTVInterpolateSrfs(SrfList, OtherOrder, CAGD_END_COND_OPEN,
				   CAGD_UNIFORM_PARAM, NULL);
    else
        TV = TrivTVFromSrfs(SrfList, OtherOrder, CAGD_END_COND_OPEN, NULL);

    if (TV == NULL)
	return TRIV_IGA_INVALID_TV_ID;

    return TrivIGADataManagerAddTrivariate(ArgmntID,
				   TrivIGAApplyDomainAndSeeding(ArgmntID, TV),
				   ID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Constructs a new trivariate and add it to the give arrangement as a      M
* trivariate through surfaces of the given ordered transformations of Srf.   M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.               M
*   Srf:           Base surface to transform in space & skin a volume thru.  M
*   Transforms:    Vector of transformation matrices to apply to Srf.        M
*   NumTransforms: Size of vector Transforms.                                M
*   OtherOrder:    Other order to use, in the skinning direction.            M
*   IsInterpolating: TRUE to interpolate surfaces, FALSE to approximate them.M
*   ID:		   ID to use or -1 to set a new ID.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID:   INVALID_IGA_TV_ID if error, to TV ID if successful.       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATVFromSurfaces2                                                   M
*****************************************************************************/
TrivIGATVID TrivIGATVFromSurfaces2(TrivIGAArrangementID ArgmntID,
                                   const CagdSrfStruct *Srf,
                                   IrtHmgnMatType Transforms[],
                                   int NumTransforms,
                                   unsigned int OtherOrder,
				   CagdBType IsInterpolating,
				   int ID)
{
    int i;
    TrivIGATVID RetVal;
    CagdSrfStruct *TSrf,
        *SrfList = NULL;

    for (i = NumTransforms - 1; i >= 0; i--) {
        TSrf = CagdSrfMatTransform(Srf, Transforms[i]);
	IRIT_LIST_PUSH(TSrf, SrfList);        
    }

    RetVal = TrivIGATVFromSurfaces(ArgmntID, SrfList,
				   OtherOrder, IsInterpolating, ID);

    CagdSrfFreeList(SrfList);

    return RetVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Refines a given TV (ID) in the arrangement.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.                M
*   TVID:         ID of the TV to refine in the arrangement.                 M
*   Dir:          Refinement direction - U, V or W.                          M
*   t:            Parameter at which to insert the knot.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVStruct *:   Reference to the IGA TV refined, or NULL if error.  M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATVRefine                                                          M
*****************************************************************************/
TrivIGATVStruct *TrivIGATVRefine(TrivIGAArrangementID ArgmntID,
				 TrivIGATVID TVID,
				 TrivTVDirType Dir,
				 CagdRType t)
{
    int i, Len;
    TrivTVStruct *TVRef,
        *TV = TrivIGADataManagerGetTrivariate(TVID);
    TrivIGATVStruct *IGATV;

    if (TV == NULL ||
	(TVRef = TrivTVRefineAtParams(TV, Dir, FALSE, &t, 1)) == NULL)
        return NULL;

    IGATV = TrivIGAUpdateTV(ArgmntID, TV, TVRef);

    /* Reset the control points IDs. */
    Len = TRIV_TV_UPT_LST_LEN(TVRef) *
	  TRIV_TV_VPT_LST_LEN(TVRef) *
	  TRIV_TV_WPT_LST_LEN(TVRef);
    IritFree(IGATV -> CtlPtsIDs);
    IGATV -> CtlPtsIDs = (TrivIGACtlPtUniqueIDsStruct *)
                       IritMalloc(sizeof(TrivIGACtlPtUniqueIDsStruct) * Len);
    for (i = 0; i < Len; i++)
        IGATV -> CtlPtsIDs[i].ID = -1;

    return IGATV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Degree raise a given TV (ID) in the arrangement.                         M
*                                                                            M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.                M
*   TVID:         ID of the TV to degree raise in the arrangement.           M
*   Dir:          Degree raising direction - U, V or W.                      M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVStruct *:   Reference to the IGA TV degree raised, or NULL if   M
*		         error.					             M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGATDegreeRaise                                                      M
*****************************************************************************/
TrivIGATVStruct *TrivIGATDegreeRaise(TrivIGAArrangementID ArgmntID,
				     TrivIGATVID TVID,
				     TrivTVDirType Dir)
{
    int i, Len;
    TrivTVStruct *TVDegRaise,
	*TV = TrivIGADataManagerGetTrivariate(TVID);
    TrivIGATVStruct *IGATV;

    if (TV == NULL || (TVDegRaise = TrivTVDegreeRaise(TV, Dir)) == NULL)
        return NULL;

    IGATV = TrivIGAUpdateTV(ArgmntID, TV, TVDegRaise);

    /* Reset the control points IDs. */
    Len = TRIV_TV_UPT_LST_LEN(TVDegRaise) *
	  TRIV_TV_VPT_LST_LEN(TVDegRaise) *
	  TRIV_TV_WPT_LST_LEN(TVDegRaise);
    IritFree(IGATV -> CtlPtsIDs);
    IGATV -> CtlPtsIDs = (TrivIGACtlPtUniqueIDsStruct *)
                       IritMalloc(sizeof(TrivIGACtlPtUniqueIDsStruct) * Len);
    for (i = 0; i < Len; i++)
        IGATV -> CtlPtsIDs[i].ID = -1;

    return IGATV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns a list of all trivariates in the given arrangement               M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.               M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGATVID *:   A dynamically allocated vector of IDs, terminated with  M
*                    INVALID_IGA_TV_ID.					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetAllTVs                                                         M
*****************************************************************************/
TrivIGATVID *TrivIGAGetAllTVs(TrivIGAArrangementID ArgmntID)
{
    int i, n;
    TrivIGATVID
	*IDs = NULL;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGAFieldStruct *Field;

    if (H == NULL)
        return NULL;
    if (H -> Fields == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }

    n = 10;
    IDs = (TrivIGATVID *) IritMalloc(sizeof(TrivIGATVID) * n);
    for (i = 0, Field = H -> Fields; Field != NULL; Field = Field -> Pnext) {
        TrivIGATVStruct *IGATV;

	for (IGATV = Field -> TVs; IGATV != NULL; IGATV = IGATV -> Pnext) {
	    IDs[i++] = TrivIGADataManagerGetTrivID(IGATV -> TV);
	    if (i >= n) {
		IDs = IritRealloc(IDs, sizeof(TrivIGATVID) * n,
				  sizeof(TrivIGATVID) * n * 2);
		n *= 2;
	    }
	}
    }
    IDs[i] = TRIV_IGA_INVALID_TV_ID;
    return IDs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Get a pointer to the TV, given its TV ID in the arrangement.             M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   TVID:	    The TV ID to search with.                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:    Found TV or NULL if not found.                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetTV                                                             M
*****************************************************************************/
TrivTVStruct *TrivIGAGetTV(TrivIGAArrangementID ArgmntID, TrivIGATVID TVID)
{
    return TrivIGADataManagerGetTrivariate(TVID);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Return indices of ALL control points in a face of designated trivariate. M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   TVID:	    The TV ID to seek its face.                              M
*   FaceID:         0 to 5 for (UMin, UMax, VMin, VMax, WMin, WMax)          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:         A vector of all IDs, terminated with -1 or NULL if error. M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetTVFaceCtlPtsIDs                                                M
*****************************************************************************/
int *TrivIGAGetTVFaceCtlPtsIDs(TrivIGAArrangementID ArgmntID,
			       TrivIGATVID TVID,
			       int FaceID)
{
    int i, j, k, n, m, Len;
    TrivIGATVID
	*IDs = NULL;
    TrivIGATVStruct
        *IGATV = TrivIGADataManagerGetIGATrivariate(TVID);
    TrivTVStruct *TV;
    TrivIGACtlPtUniqueIDsStruct *TVCtlPtsIDs;

    if (IGATV == NULL || (TV = IGATV -> TV) == NULL)
        return NULL;

    TVCtlPtsIDs = IGATV -> CtlPtsIDs;
    Len = TRIV_TV_UPT_LST_LEN(TV) *
	  TRIV_TV_VPT_LST_LEN(TV) *
	  TRIV_TV_WPT_LST_LEN(TV);

    IDs = (TrivIGATVID *) IritMalloc(sizeof(TrivIGATVID) * (Len + 1));

    for (k = n = m = 0; k < TV -> WLength; k++) {
        for (j = 0; j < TV -> VLength; j++) {
	    for (i = 0; i < TV -> ULength; i++, n++) {
	        switch (FaceID) {
		    case 0:
		        if (i == 0)
			    IDs[m++] = TVCtlPtsIDs[n].ID;
		        break;
		    case 1:
		        if (i == TV -> ULength - 1)
			    IDs[m++] = TVCtlPtsIDs[n].ID;
		        break;
		    case 2:
		        if (j == 0)
			    IDs[m++] = TVCtlPtsIDs[n].ID;
		        break;
		    case 3:
		        if (j == TV -> VLength - 1)
			    IDs[m++] = TVCtlPtsIDs[n].ID;
		        break;
		    case 4:
		        if (k == 0)
			    IDs[m++] = TVCtlPtsIDs[n].ID;
		        break;
		    case 5:
		        if (k == TV -> WLength - 1)
			    IDs[m++] = TVCtlPtsIDs[n].ID;
		        break;
		    default:
		        IritFree(IDs);
		        return NULL;
		}
	    }
	}
    }

    IDs[m] = TRIV_IGA_INVALID_TV_CTLPT_ID;

    return IDs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construct linear constraints hooking adjacent faces that do not share    M
* a common function space.                                                   M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   CallbackData:   Relevant data for the callback function		     M
*		    (such as XML information).				     M
*   NeighboringConstraintCallBack:  Call back function to call with the      M
*                   constraint as a string.                                  M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGenNeighboringConstraints                                         M
*****************************************************************************/
void TrivIGAGenNeighboringConstraints(TrivIGAArrangementID ArgmntID,
				      void *CallbackData,
				      TrivIGANeighboringConstraintCallBackType
				                NeighboringConstraintCallBack)
{
    int FaceID;
    const TrivIGAFieldStruct *F;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    const TrivIGATVStruct *TV;

    for (F = H -> Fields; F != NULL; F = F -> Pnext) {
	for (TV = F -> TVs; TV != NULL; TV = TV -> Pnext) {	    
	    TrivIGAAdjacencyInfoStruct AdjInfo[6];

	    if (TrivIGAGetFaceNeighboringTVs(ArgmntID, TV -> TV, AdjInfo)) {
	        for (FaceID = 0; FaceID < 6; FaceID++) {  /* 6 faces of TV. */
		    if (AdjInfo[FaceID].AdjTV != NULL &&
			AdjInfo[FaceID].SameSpace == FALSE) {
		        /* We have two common faces that do not share       */
			/* function space but still are the same - generate */
		        /* constraints for all control points not identical.*/
		        const TrivTVStruct
			    *AdjTV = AdjInfo[FaceID].AdjTV;
			int AdjFaceID = AdjInfo[FaceID].AdjBndryIdx;

			if (((IritIntPtrSizeType) AdjTV) >
			    ((IritIntPtrSizeType) TV -> TV) ||
			    (AdjTV == TV -> TV && AdjFaceID >= FaceID))
			    continue;           /* Do it once for the pair. */

			TrivIgaGenOneFaceNeighboringConstraints(
				      ArgmntID, NeighboringConstraintCallBack,
				      TV -> TV, FaceID, AdjTV, AdjFaceID,
				      CallbackData);
		    }
		}
	    }
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Construct linear constraints hooking one adjacent face that does not     M
* share a common function space with its neighbor while being the same.      M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   NeighboringConstraintCallBack:  Call back function to call with the      M
*                   constraint as a string.                                  M
*   TV1, FaceID1:   Description of first face.				     M
*   TV2, FaceID2:   Description of second face.				     M
*   CallbackData:   Relevant data for the callback function		     M
*		    (such as XML information).				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:  TRUE if successful, FALSE otherwise.                               M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIgaGenOneFaceNeighboringConstraints                                  M
*****************************************************************************/
int TrivIgaGenOneFaceNeighboringConstraints(
				     TrivIGAArrangementID ArgmntID,
				     TrivIGANeighboringConstraintCallBackType
					       NeighboringConstraintCallBack,
				     const TrivTVStruct *TV1,
				     int FaceID1,
				     const TrivTVStruct *TV2,
				     int FaceID2, 
				     void *CallbackData)
{
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    int i, i1, i2, n1, n2, IsRational, NumCoords,
        TV1ID = TrivIGADataManagerGetTrivID(TV1),
        TV2ID = TrivIGADataManagerGetTrivID(TV2),
        *CtlPtsIDs1 = TrivIGAGetTVFaceCtlPtsIDs(ArgmntID, TV1ID, FaceID1),
        *CtlPtsIDs2 = TrivIGAGetTVFaceCtlPtsIDs(ArgmntID, TV2ID, FaceID2);
    CagdSrfStruct *TSrf, *Srf1c, **Srf1All,
        *Srf1 = TrivIGAGetTVFaceAsSrf(ArgmntID, TV1ID, FaceID1),
        *Srf2 = TrivIGAGetTVFaceAsSrf(ArgmntID, TV2ID, FaceID2);

    if (CtlPtsIDs1 == NULL || CtlPtsIDs2 == NULL) {
        H -> LastError = TRIV_IGA_ERR_ID_NOT_FOUND;
        return FALSE;
    }

    /* Verify the lengths of the control meshes. */
#ifdef DEBUG
    for (n1 = 0;
	 CtlPtsIDs1[n1] != TRIV_IGA_INVALID_TV_CTLPT_ID &&
	     n1 <= Srf1 -> ULength * Srf1 -> VLength;
	 n1++);
    for (n2 = 0;
	 CtlPtsIDs2[n2] != TRIV_IGA_INVALID_TV_CTLPT_ID &&
	     n2 <= Srf2 -> ULength * Srf2 -> VLength;
	 n2++);
    assert(n1 == Srf1 -> ULength * Srf1 -> VLength &&
	   n2 == Srf2 -> ULength * Srf2 -> VLength);
#else
    n1 = Srf1 -> ULength * Srf1 -> VLength;
    n2 = Srf2 -> ULength * Srf2 -> VLength;
#endif /* DEBUG */


    if (n2 < n1) {
        /* Swap so Srf2 is the refined/degree raised version. */
        IRIT_SWAP(CagdSrfStruct *, Srf1, Srf2);
	IRIT_SWAP(int *, CtlPtsIDs1, CtlPtsIDs2);
	IRIT_SWAP(int, TV1ID, TV2ID);
	IRIT_SWAP(int, FaceID1, FaceID2);
	IRIT_SWAP(int, n1, n2);
    }

    /* Coerce both srfs to E1 as it is enough for what we are about to do. */
    IsRational = CAGD_IS_RATIONAL_SRF(Srf1);
    NumCoords = CAGD_NUM_OF_PT_COORD(Srf1 -> PType);

    TSrf = CagdCoerceSrfTo(Srf1, CAGD_PT_E1_TYPE, FALSE);
    CagdSrfFree(Srf1);
    Srf1 = TSrf;
    IRIT_ZAP_MEM(Srf1 -> Points[1], sizeof(CagdRType) * n1);
    TSrf = CagdCoerceSrfTo(Srf2, CAGD_PT_E1_TYPE, FALSE);
    CagdSrfFree(Srf2);
    Srf2 = TSrf;

    Srf1All = (CagdSrfStruct **) IritMalloc(sizeof(CagdSrfStruct *) * n1);
    for (i1 = 0; i1 < n1; i1++) {
	/* Elevate Srf1t to the space of Srf2 with only one ctlpt valid. */
	Srf1c = CagdSrfCopy(Srf1);
	Srf1c -> Points[1][i1] = 1.0;
	if (!CagdMakeSrfsCompatible(&Srf1c, &Srf2, TRUE, TRUE, TRUE, TRUE)) {
	    H -> LastError = TRIV_IGA_ERR_SRFS_NOT_COMPATIBLE;
	    return FALSE;
	}

	Srf1All[i1] = Srf1c;
    }

    /* Send the matrices as constraints: */
    for (i2 = 0; i2 < n2; i2++) {
        /* Do not generate constraints for identical points on both sides   */
        /* as they will be merged into one point anyway.                    */
        for (i = 0; i < n1; i++)
	    if (CtlPtsIDs2[i2] == CtlPtsIDs1[i])
	        break;
	if (CtlPtsIDs2[i2] == CtlPtsIDs1[i])
	    continue;

        NeighboringConstraintCallBack(H, IsRational, NumCoords,
				      CtlPtsIDs1, n1, CtlPtsIDs2[i2], i2,
				      (const CagdSrfStruct **) Srf1All,
				      CallbackData);
    }

    for (i1 = 0; i1 < n1; i1++)
        CagdSrfFree(Srf1All[i1]);
    IritFree(Srf1All);

    CagdSrfFree(Srf1);
    CagdSrfFree(Srf2);

    IritFree(CtlPtsIDs1);
    IritFree(CtlPtsIDs2);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Return indices of ALL control points in a face of designated trivariate. M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   TVID:	    The TV ID to seek its face.                              M
*   FaceID:         0 to 5 for (UMin, UMax, VMin, VMax, WMin, WMax)          M
*                                                                            *
* RETURN VALUE:                                                              M
*   CagdSrfStruct *:   Sought source of NULL if error.                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetTVFaceAsSrf                                                    M
*****************************************************************************/
CagdSrfStruct *TrivIGAGetTVFaceAsSrf(TrivIGAArrangementID ArgmntID,
				     TrivIGATVID TVID,
				     int FaceID)
{
    CagdSrfStruct *Srf;
    TrivIGATVStruct
        *IGATV = TrivIGADataManagerGetIGATrivariate(TVID);
    const TrivTVStruct *TV;

    if (IGATV == NULL || (TV = IGATV -> TV) == NULL)
        return NULL;

    switch (FaceID) {
        case 0:
	    Srf = TrivSrfFromMesh(TV, 0, TRIV_CONST_U_DIR);
	    break;
        case 1:
	    Srf = TrivSrfFromMesh(TV, TV -> ULength - 1, TRIV_CONST_U_DIR);
	    break;
        case 2:
	    Srf = TrivSrfFromMesh(TV, 0, TRIV_CONST_V_DIR);
	    break;
        case 3:
	    Srf = TrivSrfFromMesh(TV, TV -> VLength - 1, TRIV_CONST_V_DIR);
	    break;
        case 4:
	    Srf = TrivSrfFromMesh(TV, 0, TRIV_CONST_W_DIR);
	    break;
        case 5:
	    Srf = TrivSrfFromMesh(TV, TV -> WLength - 1, TRIV_CONST_W_DIR);
	    break;
        default:
	    return NULL;
    }

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Returns all the indices of all control points of the given TV.           M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.                 M
*   TVID:        ID of TV to get all its control points IDs.                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   int *:       A dynamically allocated vector of all IDS or NULL if error. M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetTVCtlPtsIndices                                                M
*****************************************************************************/
int *TrivIGAGetTVCtlPtsIndices(TrivIGAArrangementID ArgmntID, TrivIGATVID TVID)
{
    int i, Len;
    TrivIGATVID
	*IDs = NULL;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGATVStruct
        *IGATV = TrivIGADataManagerGetIGATrivariate(TVID);
    TrivTVStruct *TV;
    TrivIGACtlPtUniqueIDsStruct *TVCtlPtsIDs;

    if (H == NULL || IGATV == NULL)
        return NULL;

    TV = IGATV -> TV;
    TVCtlPtsIDs = IGATV -> CtlPtsIDs;
    Len = TRIV_TV_UPT_LST_LEN(TV) *
	  TRIV_TV_VPT_LST_LEN(TV) *
	  TRIV_TV_WPT_LST_LEN(TV);

    IDs = (TrivIGATVID *) IritMalloc(sizeof(TrivIGATVID) * (Len + 1));

    for (i = 0; i < Len; i++)
        IDs[i] = TVCtlPtsIDs[i].ID;
    IDs[Len] = TRIV_IGA_INVALID_TV_CTLPT_ID;

    return IDs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Fetches one control point from the arrangement. Not efficient!           M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.                 M
*   CtlPtID:     ID of control point to get.                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   const CagdCtlPtStruct *:   Fetched control point allocated statically.   M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetCtlPt                                                          M
*****************************************************************************/
const CagdCtlPtStruct *TrivIGAGetCtlPt(TrivIGAArrangementID ArgmntID,
				       int CtlPtID)
{
    IRIT_STATIC_DATA CagdCtlPtStruct
        CtlPt;
    int i, j, k;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGAFieldStruct *Field;
    TrivTVStruct *TV;
    TrivIGACtlPtUniqueIDsStruct *TVCtlPtsIDs;

    if (H == NULL)
        return NULL;
    if (H -> Fields == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return NULL;
    }

    for (i = 0, Field = H -> Fields; Field != NULL; Field = Field -> Pnext) {
        TrivIGATVStruct *IGATV;

	for (IGATV = Field -> TVs; IGATV != NULL; IGATV = IGATV -> Pnext) {
	    int Len;

	    TV = IGATV -> TV;
	    Len = TRIV_TV_UPT_LST_LEN(TV) *
	          TRIV_TV_VPT_LST_LEN(TV) *
	          TRIV_TV_WPT_LST_LEN(TV);
	    TVCtlPtsIDs = IGATV -> CtlPtsIDs;

	    for (j = 0; j < Len; j++) {
	        if (TVCtlPtsIDs[j].ID == CtlPtID) {
		    /* Fetches the i'th control point of this TV. */
		    CtlPt.PtType =  TV -> PType;
		    for (k = !TRIV_IS_RATIONAL_TV(TV);
			 k <= TRIV_NUM_OF_PT_COORD(TV );
			 k++)
		        CtlPt.Coords[k] = TV -> Points[k][j];

		    return &CtlPt;    
	        }
	    }
	}
    }

    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Fetches the material ID of a given TV ID.                                M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:    A handle on the IGA arrangement to process.                 M
*   TVID:        ID of TV to fetch its material ID.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAMaterialID:    The material fetched ID.                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAGetMaterial                                                       M
*****************************************************************************/
TrivIGAMaterialID TrivIGAGetMaterial(TrivIGAArrangementID ArgmntID,
				     TrivIGATVID TVID)
{
    int i;
    TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);
    TrivIGAFieldStruct *Field;
    TrivTVStruct
        *TV = TrivIGADataManagerGetTrivariate(TVID);

    if (H == NULL || TV == NULL)
        return TRIV_IGA_INVALID_MATERIAL_ID;
    if (H -> Fields == NULL) {
        H -> LastError = TRIV_IGA_ERR_INVALID_INPUT;
        return TRIV_IGA_INVALID_MATERIAL_ID;
    }

    for (i = 0, Field = H -> Fields; Field != NULL; Field = Field -> Pnext) {
        TrivIGATVStruct *IGATV;

	for (IGATV = Field -> TVs; IGATV != NULL; IGATV = IGATV -> Pnext) {
	    if (TV == IGATV -> TV) {
	        /* Found TV - return the field's material this TV is in. */
	        return Field -> Material -> Id;
	    }
	}
    }

    return TRIV_IGA_INVALID_MATERIAL_ID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Parses a new material description and create a new material.             M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:     A handle on the IGA arrangement to process.                M
*   MaterialStr:  material description to parse.                             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivIGAMaterialID:    New material ID or TRIV_IGA_INVALID_MATERIAL_ID if M
*                         failed.					     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGANewMaterial                                                       M
*****************************************************************************/
TrivIGAMaterialID TrivIGANewMaterial(TrivIGAArrangementID ArgmntID,
		                     const char *MaterialStr)
{
    TrivIGAMaterialID
	ID = TRIV_IGA_INVALID_MATERIAL_ID;
    TrivIGAMaterialStruct
	*Material = TrivIGAParseMaterial(MaterialStr);

    if (Material != NULL)
	ID = TrivIGAAddMaterial(ArgmntID, Material);

    return ID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Load the materials defined in the given XML file.                        M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:      A handle on the IGA arrangement to process.               M
*   FileName:      Name of material file to load.                            M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     Number of materials read.                                       M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGALoadMaterialXML                                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGALoadMaterialFromXML                                               M
*****************************************************************************/
int TrivIGALoadMaterialFromXML(TrivIGAArrangementID ArgmntID,
			       const char *FileName)
{
    int i, NumMaterials;
    TrivIGAMaterialStruct *Materials;

    if (!TrivIGALoadMaterialXML(FileName, &Materials, &NumMaterials))
	return 0;

    for (i = 0; i < NumMaterials; i++) {
	TrivIGAMaterialStruct
	    *Material = (TrivIGAMaterialStruct *) 
		                     IritMalloc(sizeof(TrivIGAMaterialStruct));

	*Material = Materials[i];
	TrivIGAAddMaterial(ArgmntID, Material);
    }

    IritFree(Materials);

    return NumMaterials;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Adds a boundary condition to a face in the arrangement.                  M
*                                                                            *
* PARAMETERS:                                                                M
   ArgmntID:       A handle on the IGA arrangement to process.		     M
*   TV:            Trivariate of the face.                                   M
*   BoundaryType:  Type of the boundary face on the trivariate.              M
*   NodeBoundaryType:   IGA boundary condition type.                         M
*   BoundaryAxisConditions: String represents the relevant axes, eg. "xy".   M
*   Value:         Value of the boundary condition.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:    TRUE on success and	FALSE on failure			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivIGAAddBoundaryFace                                                   M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAAddBoundaryFace2                                                  M
*****************************************************************************/
int TrivIGAAddBoundaryFace2(TrivIGAArrangementID ArgmntID,
			    const TrivTVStruct *TV,
			    TrivTVBndryType BoundaryType,
			    TrivIGANodeBoundaryType NodeBoundaryType,
			    const char *BoundaryAxisConditions,
			    CagdRType Value)
{
    fprintf(stderr, "TrivIGAAddBoundaryFace2 (%d) has been invoked - not implemented\n",
	    ArgmntID);
    return -1;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Function is not supported.                                               M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:       A handle on the IGA arrangement to process.              M
*   TV:             N.S.F.I.                                                 M
*   CtrlPointIndex: N.S.F.I.                                                 M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     N.S.F.I.                                                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAAddBoundaryNode                                                   M
*****************************************************************************/
int TrivIGAAddBoundaryNode(TrivIGAArrangementID ArgmntID,
                           TrivIGATVID TV,
                           int CtrlPointIndex)
{
    fprintf(stderr, "TrivIGAAddBoundaryNode (%d) has been invoked\n",
	    ArgmntID);
    return -1;
}

