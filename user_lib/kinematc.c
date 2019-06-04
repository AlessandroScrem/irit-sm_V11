/*****************************************************************************
* kinematc.c - Solves some kinematic motion problems using the MV solver.    *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Michael Barton and Gershon Elber		     November 2008   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/mvar_lib.h"
#include "inc_irit/grap_lib.h"
#include "user_loc.h"

#define USER_CNSTR_PT_PLANAR(P) \
    (P -> Type == USER_KNMTCS_PT_XY_PLANE || \
     ((P -> Type == USER_KNMTCS_PT_FIXED || \
       P -> Type == USER_KNMTCS_PT_X_DIRECTION || \
       P -> Type == USER_KNMTCS_PT_Y_DIRECTION || \
       P -> Type == USER_KNMTCS_PT_Z_DIRECTION) && \
       P -> Pt[2] == 0.0) || \
     (P -> Type == USER_KNMTCS_PT_MOVES_ALONG_CURVE && \
      CAGD_NUM_OF_PT_COORD(P -> U.Crv -> PType) < 3))

IRIT_STATIC_DATA MvarPolylineStruct
    *SolPolys = NULL;
IRIT_STATIC_DATA UserKnmtcsStruct GlblSystem;

CagdRType *UserKnmtcsDomainMin, *UserKnmtcsDomainMax;

static int UserKnmtcsCreateVarsIndex(UserKnmtcsStruct *System);
static int UserKnmtcsCreateEqnsIndex(UserKnmtcsStruct *System,
				     int *Inequalities);
static MvarMVStruct *UserKnmtcsProvideMVWithDomain(MvarMVStruct *MV, 
						   int NumOfVars,
						   CagdRType *DomainMin,
						   CagdRType *DomainMax);
static MvarMVStruct *UserKnmtcsCreateMVCAM(UserKnmtcsPtStruct *Pt);
static MvarMVStruct *UserKnmtcsCreateMVVecFromPt(UserKnmtcsPtStruct *Pt, 
						 int VarIdx, 
						 int NumOfUnknowns);
static MvarMVStruct *UserKnmtcsCreateMVVecFromBar(UserKnmtcsPtStruct *Pt1, 
						  UserKnmtcsPtStruct *Pt2,
						  int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVDstCAMPtPt(UserKnmtcsPtStruct *CAMPt, 
						  UserKnmtcsPtStruct *Pt,
						  CagdRType Distance,
						  int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVDstPtPt(UserKnmtcsPtStruct *Pt1, 
					       UserKnmtcsPtStruct *Pt2,
					       CagdRType Distance,
				     int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVAngleBarBar(UserKnmtcsBarStruct *Bar1, 
						   UserKnmtcsBarStruct *Bar2,
						   CagdRType Angle,
						   int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVOrthoBarBar(UserKnmtcsBarStruct *Bar1, 
						   UserKnmtcsBarStruct *Bar2,
						   int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVDstPtBar(UserKnmtcsPtStruct *Pt, 
						UserKnmtcsBarStruct *Bar,
						CagdRType Distance,
						int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVPar2Vecs(MvarMVStruct *MV1,
						MvarMVStruct *MV2);
static MvarMVStruct **UserKnmtcsCreateMVTanBarCrv(UserKnmtcsPtStruct *Pt, 
						  UserKnmtcsBarStruct *Bar,
						  int NumOfVars);
static MvarMVStruct **UserKnmtcsCreateMVDstPtCrv(UserKnmtcsPtStruct *CrvPt, 
						 UserKnmtcsPtStruct *Pt,
						 CagdRType Distance,
						 int NumOfVars);
static MvarMVStruct **UserKnmtcsCreateMVDstPtSrf(UserKnmtcsPtStruct *SrfPt, 
						 UserKnmtcsPtStruct *Pt,
						 CagdRType Distance,
						 int NumOfVars);
static MvarMVStruct *UserKnmtcsCreateMVParBarBar(UserKnmtcsBarStruct *Bar1, 
						 UserKnmtcsBarStruct *Bar2,
						 int NumOfVars);
static void UserKnmtcsCreateMotionDomains(UserKnmtcsStruct *System,
					  CagdRType *DomainMin,
					  CagdRType *DomainMax);
static void UserKnmtcsEvalOneSolPt(MvarPtStruct *SolPt,
				   UserKnmtcsPtStruct *MP);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Goes over points of the kinematic system and provides each unknown with  *
* a unique index. (Every point might be involved in more constraints and     *
* during the promoting process has to be shifted to the same position.).     *
*									     *
* PARAMETERS:                                                                *
*   System:	Kinematic mechanism.			 		     *
*									     *
* RETURN VALUE:                                                              *
*   int:    Number of unknowns.						     *
*****************************************************************************/
static int UserKnmtcsCreateVarsIndex(UserKnmtcsStruct *System)
{
    UserKnmtcsPtStruct *MP;
    int Index = 0;

    for (MP = System -> Pts; MP != NULL; MP = MP -> Pnext) {
	switch (MP -> Type) {
	    case USER_KNMTCS_PT_NONE:
		MP -> Idx = Index;
		Index = Index + 3;
		break;
	    case USER_KNMTCS_PT_XY_PLANE:
	    case USER_KNMTCS_PT_MOVES_ALONG_SURFACE:
		MP -> Idx = Index;
		Index = Index + 2;
		break;
	    case USER_KNMTCS_PT_Y_DIRECTION:
	    case USER_KNMTCS_PT_X_DIRECTION:
	    case USER_KNMTCS_PT_Z_DIRECTION:
	    case USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE:
	    case USER_KNMTCS_PT_MOVES_ALONG_CURVE:	    
		MP -> Idx = Index;
		Index = Index + 1;
		break;
	    case USER_KNMTCS_PT_FIXED:
		/* Point is fixed, none of its coordinates are unknowns. */
		break;
	    default:
		assert(0);	
		break;    
	}
    }
    return Index;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Goes over the kinematic system and counts the number of equations.	     *
*   Some constraints are expressed by more than one equation.		     *
*									     *
* PARAMETERS:                                                                *
*   System:        Kinematic mechanism.			 		     *
*   Inequalities:  Returns the number of inequality constraints found.	     *
*									     *
* RETURN VALUE:                                                              *
*   int:		Number of equations.				     *
*****************************************************************************/
static int UserKnmtcsCreateEqnsIndex(UserKnmtcsStruct *System,
				     int *Inequalities)
{
    UserKnmtcsConstrStruct *C;
    int Index = 0;

    *Inequalities = 0;

    for (C = System -> Constraints; C != NULL; C = C -> Pnext) {
	switch (C -> Type) {
	    case USER_KNMTCS_CONSTR_DIST_PT_PT:
	    case USER_KNMTCS_CONSTR_DIST_PT_BAR:
	    case USER_KNMTCS_CONSTR_ANGLE:
	    case USER_KNMTCS_CONSTR_ORTHOGONALITY:
	    case USER_KNMTCS_CONSTR_PARALLEL:
		Index = Index + 1;
		break;
	    case USER_KNMTCS_CONSTR_TANGENCY:
		Index = Index + 2;
		break;
	    case USER_KNMTCS_CONSTR_DIST_PT_CRV:
		Index = Index + 2;
		break;
	    case USER_KNMTCS_CONSTR_DIST_PT_SRF:
		Index = Index + 3;
		break;
	    case USER_KNMTCS_CONSTR_MIN_DIST_PT_PT:
	    case USER_KNMTCS_CONSTR_MAX_DIST_PT_PT:
	    case USER_KNMTCS_CONSTR_XDIFF_POS:
	    case USER_KNMTCS_CONSTR_YDIFF_POS:
	    case USER_KNMTCS_CONSTR_ZDIFF_POS:
	        (*Inequalities)++;
	        break;
	    default:
		assert(0);	
		break;    
	}
    }
    return Index;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Mvar is made sure to be in the B-spline representation and is            *
* transformed to the appropriate domain.				     *
*									     *
* PARAMETERS:                                                                *
*   MV:			  Multivariate system.				     *
*   NumOfVars:		  Number of variables.	 			     *
*   DomainMin, DomainMax: Desired domain.	 			     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	new allocated MV, as a Bspline in proper domain.     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsProvideMVWithDomain(MvarMVStruct *MV, 
						   int NumOfVars,
						   CagdRType *DomainMin,
						   CagdRType *DomainMax)
{
    int j;
    MvarMVStruct *NewMV, *MVTmp;

    NewMV = MvarMVCopy(MV);

    if (MVAR_IS_BEZIER_MV(MV)) {

	for (j = 0; j < NumOfVars; j++) {
	    MVTmp = MvarMVRegionFromMV(NewMV, DomainMin[j], DomainMax[j], j);
	    MvarMVFree(NewMV);
	    NewMV = MVTmp;
	}

	MVTmp = MvarCnvrtBzr2BspMV(NewMV);
	MvarMVFree(NewMV);
	NewMV = MVTmp;
    }

    for (j = 0; j < NumOfVars; j++) {
        BspKnotAffineTransOrder2(NewMV -> KnotVectors[j],
				 NewMV -> Orders[j],
				 NewMV -> Orders[j] + NewMV -> Lengths[j],
				 DomainMin[j], DomainMax[j]);
    }

    return NewMV;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates an MV that represents a y-coordiate of a CAM-moved point (a      *
* curve rotates, forcing the point move in the direction of y-axis.)	     *
*									     *
* PARAMETERS:                                                                *
*   Pt:			CAM moving point.		 		     *
*   NumOfVars:		Number of variables in the system.		     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate representing y-axis motion of a point.  *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVCAM(UserKnmtcsPtStruct *Pt)
{   
    int i;
    MvarMVStruct *MVTemp1, *MVTemp2, *MVY, **MVComp, *MV;
    CagdRType Translate[1];

    assert(Pt -> Type == USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE);

    /* Curve [x(t), y(t)] rotates along center [Cx, Cy]:                   */
    /* MVY = (x(t) - Cx)^2 + (y(t) - Cy)^2 + Cy.                           */

    MV = MvarCrvToMV(Pt -> U.Crv);
    MVComp = MvarMVSplitScalar(MV);
    MvarMVFree(MV);

    Translate[0] = -Pt -> Center.Pt[0];
    MvarMVTransform(MVComp[1], Translate, 1.0);

    Translate[0] = -Pt -> Center.Pt[1];
    MvarMVTransform(MVComp[2], Translate, 1.0);

    MVTemp1 = MvarMVMult(MVComp[1], MVComp[1]);
    MVTemp2 = MvarMVMult(MVComp[2], MVComp[2]);

    MVY = MvarMVAdd(MVTemp1, MVTemp2);

    MvarMVFree(MVTemp1);
    MvarMVFree(MVTemp2);

    for (i = 0; i < MVAR_MAX_PT_SIZE; i++) { 
	MvarMVFree(MVComp[i]);
    }	

    return MVY;  
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates a vector of Mvars which describes the kinematic point.	     *
*   For example: [x6, const] corresponds to a kinematic point, which moves   *
*   horizontally, and unknown x-coordinate is the 6th variable.(VarIdx=6)    *
*									     *
* PARAMETERS:                                                                *
*   Pt:			Kinematic point.		 		     *
*   VarIdx:		Variable index.		 			     *
*   NumOfUnknowns:	Number of variables in the system.		     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate representing Pt.			     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVVecFromPt(UserKnmtcsPtStruct *Pt, 
						 int VarIdx, 
						 int NumOfUnknowns)
{   
    int i,
	Length = 1,
	Lengths = 2;
    MvarMVStruct *MVVec[MVAR_MAX_PT_SIZE], *MVconst, *MVlin, *MVTmp, **MVTmps;
    CagdRType Translate[1];

    MVTmps = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * 4);

    assert(VarIdx < MVAR_MAX_PT_SIZE - 1);

    IRIT_ZAP_MEM(MVVec, sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);  

    /* Linear function f(x) = x. */
    MVlin = MvarBzrMVNew(1, &Lengths, MVAR_PT_E1_TYPE); 
    MVlin -> Points[1][0] = 0;
    MVlin -> Points[1][1] = 1;

    /* Constant function f(x) = 1. */
    MVconst = MvarBzrMVNew(0, &Length, MVAR_PT_E1_TYPE); 
    MVconst -> Points[1][0] = 1;

    Translate[0] = 0;

    switch (Pt -> Type) {
	case USER_KNMTCS_PT_FIXED:
	    MVTmp = MvarMVCopy(MVconst);
	    MvarMVTransform(MVTmp, Translate, Pt -> Pt[0]);
	    MVTmps[1] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, 0);
	    MvarMVFree(MVTmp);

	    MVTmp = MvarMVCopy(MVconst);
	    MvarMVTransform(MVTmp, Translate, Pt -> Pt[1]);
	    MVTmps[2] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, 0);
	    MvarMVFree(MVTmp);

	    MvarMVTransform(MVconst, Translate, Pt -> Pt[2]);
	    MVTmps[3] = MvarPromoteMVToMV2(MVconst, NumOfUnknowns, 0);
	    break;
	case USER_KNMTCS_PT_X_DIRECTION:	    
	    MVTmps[1] = MvarPromoteMVToMV2(MVlin, NumOfUnknowns, VarIdx);

	    MVTmp = MvarMVCopy(MVconst);
	    MvarMVTransform(MVTmp, Translate, Pt -> Pt[1]);
	    MVTmps[2] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, 0);
	    MvarMVFree(MVTmp);

	    MvarMVTransform(MVconst, Translate, Pt -> Pt[2]);
	    MVTmps[3] = MvarPromoteMVToMV2(MVconst, NumOfUnknowns, 0);
	    break;
	case USER_KNMTCS_PT_Y_DIRECTION:
	    MVTmp = MvarMVCopy(MVconst);
	    MvarMVTransform(MVTmp, Translate, Pt -> Pt[0]);
	    MVTmps[1] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, 0);
	    MvarMVFree(MVTmp);
	    
	    MVTmps[2] = MvarPromoteMVToMV2(MVlin, NumOfUnknowns, VarIdx);

	    MvarMVTransform(MVconst, Translate, Pt -> Pt[2]);
	    MVTmps[3] = MvarPromoteMVToMV2(MVconst, NumOfUnknowns, 0);
	    break;
	case USER_KNMTCS_PT_Z_DIRECTION:
	    MVTmp = MvarMVCopy(MVconst);
	    MvarMVTransform(MVTmp, Translate, Pt -> Pt[0]);
	    MVTmps[1] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, 0);
	    MvarMVFree(MVTmp);

	    MvarMVTransform(MVconst, Translate, Pt -> Pt[1]);
	    MVTmps[2]= MvarPromoteMVToMV2(MVconst, NumOfUnknowns, 0);

	    MVTmps[3] = MvarPromoteMVToMV2(MVlin, NumOfUnknowns, VarIdx);
	    break;
	case USER_KNMTCS_PT_XY_PLANE:
	    MVTmp = MvarMVCopy(MVlin);
	    MVTmps[1] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, VarIdx);
	    MvarMVFree(MVTmp);

	    MVTmps[2] = MvarPromoteMVToMV2(MVlin, NumOfUnknowns, VarIdx + 1);

	    MvarMVTransform(MVconst, Translate, 0);
	    MVTmps[3] = MvarPromoteMVToMV2(MVconst, NumOfUnknowns, 0);
	    break;
	case USER_KNMTCS_PT_NONE:
	    MVTmp = MvarMVCopy(MVlin);
	    MVTmps[1] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, VarIdx);
	    MvarMVFree(MVTmp);

	    MVTmp = MvarMVCopy(MVlin);
	    MVTmps[2] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, VarIdx + 1);
	    MvarMVFree(MVTmp);

	    MVTmps[3] = MvarPromoteMVToMV2(MVlin, NumOfUnknowns, VarIdx + 2);
	    break;
	case USER_KNMTCS_PT_MOVES_ALONG_CURVE:
	    {
	        MvarMVStruct *MV, **MVComp;

		if (CAGD_IS_BSPLINE_CRV(Pt -> U.Crv) &&
		    !BspCrvHasOpenEC(Pt -> U.Crv)) {
		    CagdCrvStruct
		        *TCrv = BspCrvOpenEnd(Pt -> U.Crv);

		    MV = MvarCrvToMV(TCrv);
		    CagdCrvFree(TCrv);
		}
		else
		    MV = MvarCrvToMV(Pt -> U.Crv);

		MVComp = MvarMVSplitScalar(MV);

		for (i = 1; i < 4; i++) { 
		    MVTmps[i] = MvarPromoteMVToMV2(MVComp[i], NumOfUnknowns,
						   VarIdx);
		}
		
		for (i = 0; i < MVAR_MAX_PT_SIZE; i++) { 
		    MvarMVFree(MVComp[i]);
		}
		MvarMVFree(MV);
	    }    
	    break;
	case USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE:
	    {
	        MvarMVStruct *MV, **MVComp;

		if (CAGD_IS_RATIONAL_CRV(Pt -> U.Crv))
		    USER_FATAL_ERROR(USER_ERR_RATIONAL_NO_SUPPORT);

		if (CAGD_IS_BSPLINE_CRV(Pt -> U.Crv) &&
		    !BspCrvHasOpenEC(Pt -> U.Crv))
		    USER_FATAL_ERROR(USER_ERR_INVALID_KV_END_COND);

		MV = MvarCrvToMV(Pt -> U.Crv);
		MVComp = MvarMVSplitScalar(MV);

		MVTmp = MvarMVCopy(MVconst);
		MvarMVTransform(MVTmp, Translate, Pt -> Center.Pt[0]);
		MVTmps[1] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, 0);
		MvarMVFree(MVTmp);    

		MvarMVTransform(MVconst, Translate, 0);
		MVTmps[3] = MvarPromoteMVToMV2(MVconst, NumOfUnknowns, 0);

		MVTmp = UserKnmtcsCreateMVCAM(Pt);
		MVTmps[2] = MvarPromoteMVToMV2(MVTmp, NumOfUnknowns, VarIdx);

		for (i = 0; i < MVAR_MAX_PT_SIZE; i++) { 
		    MvarMVFree(MVComp[i]);
		}
		MvarMVFree(MV);
	    }    
	    break;
	case USER_KNMTCS_PT_MOVES_ALONG_SURFACE:
	    {
	        MvarMVStruct *MV, **MVComp;

		if (CAGD_IS_RATIONAL_SRF(Pt -> U.Srf))
		    USER_FATAL_ERROR(USER_ERR_RATIONAL_NO_SUPPORT);

		if (CAGD_IS_BSPLINE_SRF(Pt -> U.Srf) &&
		    !BspSrfHasOpenEC(Pt -> U.Srf)) {
		    CagdSrfStruct
		        *TSrf = BspSrfOpenEnd(Pt -> U.Srf);

		    MV = MvarSrfToMV(TSrf);
		    CagdSrfFree(TSrf);
		}
		else
		    MV = MvarSrfToMV(Pt -> U.Srf);

    		MVComp = MvarMVSplitScalar(MV);

		for (i = 1; i < 4; i++) { 
		    MvarMVStruct
		        *Temp = MvarPromoteMVToMV2(MVComp[i], NumOfUnknowns,
						   VarIdx);

		    if (Temp != NULL)
		        MVTmps[i] = Temp;
		    else
		        MVTmps[i] = MvarMVCopy(MVComp[i]);
		}

		for (i = 0; i < MVAR_MAX_PT_SIZE; i++) { 
		    MvarMVFree(MVComp[i]);
		}
	    }	    
	    break;
	default:
	    /* Point is fixed, none of its coordinates are unknowns. */
	    break;    
    }

    MvarMVFree(MVlin);
    MvarMVFree(MVconst);

    for (i = 1; i < 4; i++) {    
	MVVec[i] = UserKnmtcsProvideMVWithDomain(MVTmps[i], NumOfUnknowns,
						 UserKnmtcsDomainMin,
						 UserKnmtcsDomainMax);
	MvarMVFree(MVTmps[i]);
    }

    IritFree(MVTmps);

    MVTmp = MvarMVMergeScalar(MVVec);

    for (i = 0; i < MVAR_MAX_PT_SIZE; i++)
	MvarMVFree(MVVec[i]);

    return MVTmp;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates an MV constraint which describes an axis distance for Pt1Pt2.    *
*									     *
* PARAMETERS:                                                                *
*   Pt1, Pt2:		Points representing a bar.	 		     *
*   Distance:		Axis distance (value) to preserve.		     *
*   NumOfVars:		Number of variables in the system.		     *
*   Axis:               Axis to consider.				     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate representing Bar.			     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateAxisDistFromBar(UserKnmtcsPtStruct *Pt1, 
						     UserKnmtcsPtStruct *Pt2,
						     CagdRType Distance,
						     int NumOfVars,
						     int Axis)
{
    int i;
    CagdRType Translate[1];
    MvarMVStruct *MVVec1, *MVVec2, *MVBar, **MVScalars;

    assert(Axis >= 1 && Axis <= 3);

    MVVec1 = UserKnmtcsCreateMVVecFromPt(Pt1, Pt1 -> Idx, NumOfVars);
    MVVec2 = UserKnmtcsCreateMVVecFromPt(Pt2, Pt2 -> Idx, NumOfVars);

    MVBar = MvarMVSub(MVVec1, MVVec2); 
    MvarMVFree(MVVec1);
    MvarMVFree(MVVec2);
    MVScalars = MvarMVSplitScalar(MVBar);
    MvarMVFree(MVBar);
    MVBar = MVScalars[Axis];

    for (i = 1; i < 4; i++) {
        if (i != Axis)
	    MvarMVFree(MVScalars[i]);
    }

    Translate[0] = -Distance;
    MvarMVTransform(MVBar, Translate, 1.0);

    return MVBar;  
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates an MV constraint which describes a kinematic bar Pt1Pt2.	     *
*									     *
* PARAMETERS:                                                                *
*   Pt1, Pt2:		Points representing a bar.	 		     *
*   NumOfVars:		Number of variables in the system.		     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate representing Bar.			     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVVecFromBar(UserKnmtcsPtStruct *Pt1, 
						  UserKnmtcsPtStruct *Pt2,
						  int NumOfVars)
{   
    MvarMVStruct *MVVec1, *MVVec2, *MVBar;

    MVVec1 = UserKnmtcsCreateMVVecFromPt(Pt1, Pt1 -> Idx, NumOfVars);
    MVVec2 = UserKnmtcsCreateMVVecFromPt(Pt2, Pt2 -> Idx, NumOfVars);

    MVBar = MvarMVSub(MVVec1, MVVec2); 

    MvarMVFree(MVVec1);
    MvarMVFree(MVVec2);

    return MVBar;  
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a point-point distance   *
* constraint, if one of the points is a CAM point - Special treatmnent is    *
* needed.								     *
*									     *
* PARAMETERS:                                                                *
*   CAMPt, Pt:		Two kinem. points.				     *
*   Distance:		Distance (value) to preserve. 			     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate vector.				     M
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVDstCAMPtPt(UserKnmtcsPtStruct *CAMPt, 
						  UserKnmtcsPtStruct *Pt,
						  CagdRType Distance,
						  int NumOfVars)
{
    MvarMVStruct *MV, *MVTemp1, *MVTemp2, *MVTemp3, *MVTemp4, *MVTemp5, *MVVec,
	*MVTemp6, *MVTemp7, *MVVecCAM, *MV1Split[MVAR_MAX_PT_SIZE], **MV2Split;
    CagdRType Translate[1], Cy;

    /* Center of the rotating curve [Cx; Cy].*/
    Cy = CAMPt -> Center.Pt[1];

    MVVecCAM = UserKnmtcsCreateMVVecFromPt(CAMPt, CAMPt -> Idx, NumOfVars);
    MVVec = UserKnmtcsCreateMVVecFromPt(Pt, Pt -> Idx, NumOfVars);

    IRIT_GEN_COPY(MV1Split, MvarMVSplitScalar(MVVecCAM),
		  sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
    MV2Split = MvarMVSplitScalar(MVVec);

    MVTemp1 = MvarMVSub(MV1Split[1], MV2Split[1]);
    MVTemp2 = MvarMVMult(MVTemp1, MVTemp1);

    Translate[0] = -Cy;
    MvarMVTransform(MV2Split[2], Translate, 1.0);

    MVTemp3 = MvarMVMult(MV2Split[2], MV2Split[2]);
    MVTemp4 = MvarMVAdd(MVTemp3, MV1Split[2]);	
    MVTemp5 = MvarMVAdd(MVTemp2, MVTemp4);

    Translate[0] = -IRIT_SQR(Distance);

    MvarMVTransform(MVTemp5, Translate, 1.0);

    MVTemp6 = MvarMVMult(MVTemp5, MVTemp5);
    MVTemp7 = MvarMVMult(MVTemp3, MV1Split[2]);

    MvarMVTransform(MVTemp7, NULL, 4.0);

    MV = MvarMVSub(MVTemp6, MVTemp7);

    MvarMVFree(MVTemp1);
    MvarMVFree(MVTemp2);
    MvarMVFree(MVTemp3);
    MvarMVFree(MVTemp4);
    MvarMVFree(MVTemp5);
    MvarMVFree(MVTemp6);
    MvarMVFree(MVTemp7);
    MvarMVFree(MVVecCAM);
    MvarMVFree(MVVec);

    return MV; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a point-point distance   *
* constraint.								     *
*									     *
* PARAMETERS:                                                                *
*   Pt1, 2:		Two kinem. points.				     *
*   Distance:		Distance (value) to preserve. 			     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate vector.				     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVDstPtPt(UserKnmtcsPtStruct *Pt1, 
					       UserKnmtcsPtStruct *Pt2,
					       CagdRType Distance,
					       int NumOfVars)
{
    MvarMVStruct *MVTemp, *MVBar, *MVTemp1;
    CagdRType Translate[1];

    assert(Pt1 -> Type != USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE ||
	   Pt2 -> Type != USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE);

    if (Pt1 -> Type != USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE &&
        Pt2 -> Type != USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE) {
	MVBar = UserKnmtcsCreateMVVecFromBar(Pt1, Pt2, NumOfVars);
     
	MVTemp = MvarMVDotProd(MVBar, MVBar);

	Translate[0] = -IRIT_SQR(Distance);
	MvarMVTransform(MVTemp, Translate, 1.0);	
    }
    /* The code is not complex here, the distance between CAM point  */
    /* and any other point types has to computed differently!        */
    else {
	if (Pt1 -> Type == USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE) {
	    MVTemp = UserKnmtcsCreateMVDstCAMPtPt(Pt1, Pt2, Distance,
						  NumOfVars);
	}
	else {
	    MVTemp = UserKnmtcsCreateMVDstCAMPtPt(Pt2, Pt1, Distance,
						  NumOfVars);
	}

	MVTemp1 = UserKnmtcsProvideMVWithDomain(MVTemp, NumOfVars,
						UserKnmtcsDomainMin,
						UserKnmtcsDomainMax);
	MvarMVFree(MVTemp);
	MVTemp = MVTemp1; 		
    }

    return MVTemp;   
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a bar-bar angle	     *
*   constraint.								     *
*									     *
* PARAMETERS:                                                                *
*   Bar1, Bar2:		Two kinem. bars.				     *
*   Angle:		Angle (value) to preserve in degrees. 		     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate vector.				     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVAngleBarBar(UserKnmtcsBarStruct *Bar1, 
						   UserKnmtcsBarStruct *Bar2,
						   CagdRType Angle,
						   int NumOfVars)
{
    MvarMVStruct *MVTemp1, *MVTemp2, 
	*MVBar1, *MVBar2, *MVDotProd, *MVDotProd2, *MVSqNorm1, *MVSqNorm2;
    CagdRType CosAngle;

    /* IF planar and the angular constraints is for 90 * n, n natural: */
    if (USER_CNSTR_PT_PLANAR(Bar1 -> P1) &&
	USER_CNSTR_PT_PLANAR(Bar1 -> P2) &&
	USER_CNSTR_PT_PLANAR(Bar2 -> P1) &&
	USER_CNSTR_PT_PLANAR(Bar2 -> P2)) {
	if (Angle == 90 || Angle == 270) {
	    return UserKnmtcsCreateMVOrthoBarBar(Bar1, Bar2, NumOfVars);
	}
	else if (Angle == 0 || Angle == 180) {
	    return UserKnmtcsCreateMVParBarBar(Bar1, Bar2, NumOfVars);
	}
    }

    MVBar1 = UserKnmtcsCreateMVVecFromBar(Bar1 -> P1, Bar1 -> P2, NumOfVars);
    MVBar2 = UserKnmtcsCreateMVVecFromBar(Bar2 -> P1, Bar2 -> P2, NumOfVars);

    CosAngle = cos(IRIT_DEG2RAD_CNVRT * Angle);

    MVDotProd = MvarMVDotProd(MVBar1, MVBar2);
    MVSqNorm1 = MvarMVDotProd(MVBar1, MVBar1);
    MVSqNorm2 = MvarMVDotProd(MVBar2, MVBar2);
    MVDotProd2 = MvarMVMult(MVDotProd, MVDotProd);

    MVTemp1 = MvarMVMult(MVSqNorm1, MVSqNorm2);

    MvarMVTransform(MVTemp1, NULL, IRIT_SQR(CosAngle));

    MVTemp2 = MvarMVSub(MVDotProd2, MVTemp1);

    MvarMVFree(MVTemp1);
    MvarMVFree(MVBar1);
    MvarMVFree(MVBar2);
    MvarMVFree(MVSqNorm1);
    MvarMVFree(MVSqNorm2);
    MvarMVFree(MVDotProd);
    MvarMVFree(MVDotProd2);

    return MVTemp2; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a bar-bar orthogonality  *
* constraint.								     *
*									     *
* PARAMETERS:                                                                *
*   Bar1, Bar2:		Two kinematic bars.				     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate vector.				     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVOrthoBarBar(UserKnmtcsBarStruct *Bar1, 
						   UserKnmtcsBarStruct *Bar2,
						   int NumOfVars)
{
    MvarMVStruct *MVTemp1, *MVBar1, *MVBar2;

    MVBar1 = UserKnmtcsCreateMVVecFromBar(Bar1 -> P1, Bar1 -> P2, NumOfVars);
    MVBar2 = UserKnmtcsCreateMVVecFromBar(Bar2 -> P1, Bar2 -> P2, NumOfVars);

    MVTemp1 = MvarMVDotProd(MVBar1, MVBar2);

    MvarMVFree(MVBar1);
    MvarMVFree(MVBar2);

    return MVTemp1; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a point-line distance    *
* constraint. Line is represented by Bar.				     *
*									     *
* PARAMETERS:                                                                *
*   Pt:		    Kinem. points.					     *
*   Bar:	    Kinem. bar.						     *
*   Distance:	    Distance (value) to preserve. 			     *
*   NumOfVars:	    Number of unknowns in the kinematic system.		     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Multivariate vector.				     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVDstPtBar(UserKnmtcsPtStruct *Pt, 
						UserKnmtcsBarStruct *Bar,
						CagdRType Distance,
						int NumOfVars)
{
    MvarMVStruct *MV, *MVv, *MVP, *MVA, *MVAP, *MVTemp1, *MVTemp2, *MVTemp3,
	*MVTemp4, *MVTemp5, *MVTemp6;

    MVv = UserKnmtcsCreateMVVecFromBar(Bar -> P1, Bar -> P2, NumOfVars);

    MVP = UserKnmtcsCreateMVVecFromPt(Pt, Pt -> Idx, NumOfVars);
    MVA = UserKnmtcsCreateMVVecFromPt(Bar -> P1, Bar -> P1 -> Idx, NumOfVars);

    /* Distance D of point P from line (A,v):   
	<A-P,A-P> * <v,v> - <v,A-P>^2 - D^2*<v,v> = 0. */

    MVAP = MvarMVSub(MVA, MVP);
    MVTemp1 = MvarMVDotProd(MVAP, MVAP);
    MVTemp2 = MvarMVDotProd(MVv, MVv);
    MVTemp3 = MvarMVMult(MVTemp1, MVTemp2);

    MVTemp4 = MvarMVDotProd(MVv, MVAP);
    MVTemp5 = MvarMVMult(MVTemp4, MVTemp4);

    MVTemp6 = MvarMVSub(MVTemp3, MVTemp5);

    MvarMVTransform(MVTemp2, NULL, IRIT_SQR(Distance));

    MV = MvarMVSub(MVTemp6, MVTemp2);

    MvarMVFree(MVTemp1);
    MvarMVFree(MVTemp2);
    MvarMVFree(MVTemp3);
    MvarMVFree(MVTemp4);
    MvarMVFree(MVTemp5);
    MvarMVFree(MVTemp6);
    MvarMVFree(MVv);
    MvarMVFree(MVP);
    MvarMVFree(MVA);
    MvarMVFree(MVAP);

    return MV; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a constraint that	     *
* multivariate vectors MV1 and MV2 are parallel.  			     *
*									     *
* PARAMETERS:                                                                *
*   MV1, MV2:		Vectors of multivariates.			     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct *:	Parallelism constraint multivariate.	    	     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVPar2Vecs(MvarMVStruct *MV1,
						MvarMVStruct *MV2)
{
    MvarMVStruct *MVTemp, *ScalarMVs1[MVAR_MAX_PT_SIZE], **ScalarMVs2;

    /* (u1, u2) * (-v2, v1) = 0. */ 
    if (MVAR_IS_RATIONAL_MV(MV1) || MVAR_IS_RATIONAL_MV(MV2)) {
        MVTemp = NULL;
	assert(0); 
    }
    else {
        IRIT_GEN_COPY(ScalarMVs1, MvarMVSplitScalar(MV1),
		      sizeof(MvarMVStruct *) * MVAR_MAX_PT_SIZE);
	ScalarMVs2 = MvarMVSplitScalar(MV2);

	MVTemp = MvarMVCrossProd2D(ScalarMVs1[1], ScalarMVs1[2],
				   ScalarMVs2[1], ScalarMVs2[2]);

	MVAR_FREE_SCALARS(ScalarMVs1);
	MVAR_FREE_SCALARS(ScalarMVs2);
    }

    return MVTemp; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates multivariates that corresponds to a bar-curve tangency           *
* constraint.	Curve is represented by Pt, a contact point between curve    *
* and bar.								     *
*									     *
* PARAMETERS:                                                                *
*   Pt:			Contact point between Bar and curve.		     *
*   Bar:		Kinematic bar.		 			     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct **:	Two constraints: contact and tangency.	    	     *
*****************************************************************************/
static MvarMVStruct **UserKnmtcsCreateMVTanBarCrv(UserKnmtcsPtStruct *Pt, 
						  UserKnmtcsBarStruct *Bar,
						  int NumOfVars)
{
    MvarMVStruct *MVBar, **MVTan, *MVC, *MVA, *MVDer, *MVAC;

    MVTan = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * 2);

    MVBar = UserKnmtcsCreateMVVecFromBar(Bar -> P1, Bar -> P2, NumOfVars);
    MVC = UserKnmtcsCreateMVVecFromPt(Pt, Pt -> Idx, NumOfVars);
    MVA = UserKnmtcsCreateMVVecFromPt(Bar -> P1, Bar -> P1 -> Idx, NumOfVars);
    MVDer = MvarMVDerive(MVC, Pt -> Idx);

    /* Tangency (direction) constraint C'(t)|| Bar. */
    MVTan[0] = UserKnmtcsCreateMVPar2Vecs(MVDer, MVBar);
   
    /* Colinearity constraint of points C(t), Bar -> P1, Bar -> P2. */
    MVAC = UserKnmtcsCreateMVVecFromBar(Pt, Bar -> P1, NumOfVars);
    MVTan[1] = UserKnmtcsCreateMVPar2Vecs(MVAC, MVBar);

    MvarMVFree(MVC);
    MvarMVFree(MVA);
    MvarMVFree(MVAC);
    MvarMVFree(MVDer);
    MvarMVFree(MVBar);

    return MVTan; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a point-curve distance   *
* constraint.	Curve is represented by a moving point CrvPt, a footpoint on *
* the curve.								     *
*									     *
* PARAMETERS:                                                                *
*   CrvPt:		Footpoint on the curve.				     *
*   Pt:			Point which remains constantly distant from a curve. *
*   Distance:		A distance between a point and a curve.		     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct **:	Two constraints: distance and orthogonality.   	     *
*****************************************************************************/
static MvarMVStruct **UserKnmtcsCreateMVDstPtCrv(UserKnmtcsPtStruct *CrvPt, 
						 UserKnmtcsPtStruct *Pt,
						 CagdRType Distance,
						 int NumOfVars)
{
    MvarMVStruct **MVs, *MVC, *MVDer, *MVPC;

    MVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * 2);

    /* Distance constraint between C(t) an P. */
    MVs[0] = UserKnmtcsCreateMVDstPtPt(CrvPt, Pt, Distance, NumOfVars);

    MVC = UserKnmtcsCreateMVVecFromPt(CrvPt, CrvPt -> Idx, NumOfVars);
    MVDer = MvarMVDerive(MVC, CrvPt -> Idx);
    MVPC = UserKnmtcsCreateMVVecFromBar(CrvPt, Pt, NumOfVars);

    /* Orthogonality constraint < C'(t), C - P > = 0. */
    MVs[1] = MvarMVDotProd(MVDer, MVPC);
	
    MvarMVFree(MVC);
    MvarMVFree(MVDer);
    MvarMVFree(MVPC);

    return MVs; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a point-surface distance *
* constraint.	Surface is represented by a moving point SrfPt, a footpoint  *
* the surface.								     *
*									     *
* PARAMETERS:                                                                *
*   SrfPt:	    Footpoint on the curve.				     *
*   Pt:		    Point which remains constantly distant from a surface.   *
*   Distance:	    A distance between a point and a surface.		     *
*   NumOfVars:	    Number of unknowns in the kinematic system.		     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct **:	Three constraints: distance and two orthogonalities. *
*****************************************************************************/
static MvarMVStruct **UserKnmtcsCreateMVDstPtSrf(UserKnmtcsPtStruct *SrfPt, 
						 UserKnmtcsPtStruct *Pt,
						 CagdRType Distance,
						 int NumOfVars)
{
    MvarMVStruct **MVs, *MVS, *MVDerU, *MVDerV, *MVPS;

    MVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) * 2);

    /* Distance constraint between S(u,v) and P. */
    MVs[0] = UserKnmtcsCreateMVDstPtPt(SrfPt, Pt, Distance, NumOfVars);

    MVS = UserKnmtcsCreateMVVecFromPt(SrfPt, SrfPt -> Idx, NumOfVars);
    MVDerU = MvarMVDerive(MVS, SrfPt -> Idx);
    MVDerV = MvarMVDerive(MVS, (SrfPt -> Idx)+1);
    MVPS = UserKnmtcsCreateMVVecFromBar(SrfPt, Pt, NumOfVars);

    /* Orthogonality constraints < dS/dU, S - P > = 0, < dS/dV, S - P > = 0. */
    MVs[1] = MvarMVDotProd(MVDerU, MVPS);
    MVs[2] = MvarMVDotProd(MVDerV, MVPS);
	
    MvarMVFree(MVS);
    MvarMVFree(MVDerU);
    MvarMVFree(MVDerV);
    MvarMVFree(MVPS);

    return MVs; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Creates scalar multivariate that corresponds to a bar-bar parallelism    *
* constraint.								     *
*									     *
* PARAMETERS:                                                                *
*   Bar1, Bar2:		Two kinematic bars.				     *
*   NumOfVars:		Number of unknowns in the kinematic system.	     *
*									     *
* RETURN VALUE:                                                              *
*   MvarMVStruct **:	Two constraints: contact and tangency.	    	     *
*****************************************************************************/
static MvarMVStruct *UserKnmtcsCreateMVParBarBar(UserKnmtcsBarStruct *Bar1, 
						 UserKnmtcsBarStruct *Bar2,
						 int NumOfVars)
{
    MvarMVStruct *MVBar1, *MVBar2, *MV;

    MVBar1 = UserKnmtcsCreateMVVecFromBar(Bar1 -> P1, Bar1 -> P2, NumOfVars);
    MVBar2 = UserKnmtcsCreateMVVecFromBar(Bar2 -> P1, Bar2 -> P2, NumOfVars);
    MV = UserKnmtcsCreateMVPar2Vecs(MVBar1, MVBar2);
    
    MvarMVFree(MVBar1);
    MvarMVFree(MVBar2);

    return MV; 
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Preliminary domain for the kinematic system is created.	 	     *
*									     *
* PARAMETERS:                                                                *
*   System:		    Kinematic system.		 		     *
*   DomainMin, DomainMax:   Domain to create.		 		     *
*									     *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void UserKnmtcsCreateMotionDomains(UserKnmtcsStruct *System,
					  CagdRType *DomainMin,
					  CagdRType *DomainMax)
{
    UserKnmtcsPtStruct *MP;

    /* XY should always be a valid domain. Z Can be empty if a 2D problem. */
    if (System -> XMax <= System -> XMin) {
        USER_FATAL_ERROR(USER_ERR_XRANGE_EMPTY);
    }
    if (System -> YMax <= System -> YMin) {
        USER_FATAL_ERROR(USER_ERR_YRANGE_EMPTY);
    }

    for (MP = System -> Pts; MP != NULL; MP = MP -> Pnext) {
	switch (MP -> Type) {
	    case USER_KNMTCS_PT_XY_PLANE:
		DomainMin[MP -> Idx] = System -> XMin;
		DomainMax[MP -> Idx] = System -> XMax;
		DomainMin[(MP -> Idx) + 1] = System -> YMin;
		DomainMax[(MP -> Idx) + 1] = System -> YMax;
		break;
	    case USER_KNMTCS_PT_NONE:
		DomainMin[MP -> Idx] = System -> XMin;
		DomainMax[MP -> Idx] = System -> XMax;
		DomainMin[(MP -> Idx) + 1] = System -> YMin;
		DomainMax[(MP -> Idx) + 1] = System -> YMax;
		DomainMin[(MP -> Idx) + 2] = System -> ZMin;
		DomainMax[(MP -> Idx) + 2] = System -> ZMax;
		if (System -> ZMax <= System -> ZMin) {
		    USER_FATAL_ERROR(USER_ERR_ZRANGE_EMPTY);
		}
		break;
	    case USER_KNMTCS_PT_X_DIRECTION:
		DomainMin[MP -> Idx] = System -> XMin;
		DomainMax[MP -> Idx] = System -> XMax;
		break;
	    case USER_KNMTCS_PT_Y_DIRECTION:
		DomainMin[MP -> Idx] = System -> YMin;
		DomainMax[MP -> Idx] = System -> YMax;
		break;
	    case USER_KNMTCS_PT_Z_DIRECTION:
		DomainMin[MP -> Idx] = System -> ZMin;
		DomainMax[MP -> Idx] = System -> ZMax;
		if (System -> ZMax <= System -> ZMin) {
		    USER_FATAL_ERROR(USER_ERR_ZRANGE_EMPTY);
		}
		break;
	    case USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE:
	    case USER_KNMTCS_PT_MOVES_ALONG_CURVE:
	        {
		    CagdRType TMin, TMax;

		    CagdCrvDomain(MP -> U.Crv, &TMin, &TMax);
		    DomainMin[MP -> Idx] = TMin;
		    DomainMax[MP -> Idx] = TMax;
		}
		break;
	    case USER_KNMTCS_PT_MOVES_ALONG_SURFACE:
	        {
		    CagdRType UMin, UMax, VMin, VMax;

		    CagdSrfDomain(MP -> U.Srf, &UMin, &UMax, &VMin, &VMax);
		    DomainMin[MP -> Idx] = UMin;
		    DomainMax[MP -> Idx] = UMax;
		    DomainMin[(MP -> Idx) + 1] = VMin;
		    DomainMax[(MP -> Idx) + 1] = VMax;
		}
		break;
	    case USER_KNMTCS_PT_FIXED:
		/* Point is fixed, none of its coordinates are unknowns. */
		break;
	    default:
		assert(0);	
		break;    
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Deallocates and frees the kinematic solution stored as a list of	     M
*   multi-variate polyline structures.					     M
*									     *
* PARAMETERS:                                                                M
*   None						 		     M
*									     *
* RETURN VALUE:                                                              M
*   void								     M
*									     *
* KEYWORDS:                                                                  M
*   UserKnmtcsFreeSol			                                     M
*****************************************************************************/
void UserKnmtcsFreeSol(void)
{ 
    MvarPolylineFreeList(SolPolys);
    SolPolys = NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reset all information given to the solver to initial state.		     M
*									     *
* PARAMETERS:                                                                M
*   None						 		     M
*									     *
* RETURN VALUE:                                                              M
*   void								     M
*									     *
* KEYWORDS:                                                                  M
*   UserKnmtcsSolveDone			                                     M
*****************************************************************************/
void UserKnmtcsSolveDone(void)
{
    IRIT_ZAP_MEM(&GlblSystem, sizeof(UserKnmtcsStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   For given planar kinematic mechanism, the system of constrains is 	     M
* created and solved. The list of solution points (in the parametric space)  M
* is returned.								     M
*									     *
* PARAMETERS:                                                                M
*   System:		Kinematic system.		 		     M
*   NumTol:		Numerical tolerance.		 		     M
*   SubTol:		Subdivision tolerance.		 		     M
*   Step:		"Animation" step.		 		     M
*   SolDim:             Dimension of the solution, even if invalid/error.    M
*   FilterSols:         True to filter and return only numerically           M
*                       converging solutions.				     M
*									     *
* RETURN VALUE:                                                              M
*   int:		Number of solution curves, if any. -1 if error.      M
*									     *
* KEYWORDS:                                                                  M
*   UserKnmtcsSolveMotion		                                     M
*****************************************************************************/
int UserKnmtcsSolveMotion(const UserKnmtcsStruct *System,
			  CagdRType NumTol,
		          CagdRType SubTol,
			  CagdRType Step,
			  int *SolDim,
			  CagdBType FilterSols)
{
    int i, NumOfVars, NumofInequalities, NumOfEqns, NumOfSols; 
    UserKnmtcsConstrStruct *CS;
    MvarMVStruct
	**MVs = NULL;
    MvarPtStruct *Pt,
	*IsolPts = NULL;
    MvarConstraintType *Constraints;
    MvarPolylineStruct *NewPoly;

    GlblSystem = *System;

    NumOfVars = UserKnmtcsCreateVarsIndex(&GlblSystem);
    NumOfEqns = UserKnmtcsCreateEqnsIndex(&GlblSystem, &NumofInequalities);
    *SolDim = NumOfVars - NumOfEqns;

    UserKnmtcsDomainMin =
        (CagdRType *) IritMalloc(NumOfVars * sizeof(CagdRType));
    UserKnmtcsDomainMax =
        (CagdRType *) IritMalloc(NumOfVars * sizeof(CagdRType));

    UserKnmtcsCreateMotionDomains(&GlblSystem,
				  UserKnmtcsDomainMin, UserKnmtcsDomainMax);

    MVs = (MvarMVStruct **) IritMalloc(sizeof(MvarMVStruct *) *
				       (NumOfEqns + NumofInequalities));
    Constraints = (MvarConstraintType *)
			  IritMalloc(sizeof(MvarConstraintType) *
				     (NumOfEqns + NumofInequalities));

    for (i = 0; i < NumOfEqns + NumofInequalities; i++)
        Constraints[i] = MVAR_CNSTRNT_ZERO;	       /* Default behaviour. */

    for (i = 0, CS = GlblSystem.Constraints;
	 CS != NULL;
	 CS = CS -> Pnext, i++) {
	switch (CS -> Type) {
	    case USER_KNMTCS_CONSTR_DIST_PT_PT:
		MVs[i] = UserKnmtcsCreateMVDstPtPt(CS -> C.DstPtPt.Pt1, 
						   CS -> C.DstPtPt.Pt2,
						   CS -> V.distance,
						   NumOfVars);
		break; 
	    case USER_KNMTCS_CONSTR_ANGLE:
		MVs[i] = UserKnmtcsCreateMVAngleBarBar(CS -> C.Angle.Bar1, 
						       CS -> C.Angle.Bar2,
						       CS -> V.angle,
						       NumOfVars);
		break;
	    case USER_KNMTCS_CONSTR_ORTHOGONALITY:
		MVs[i] = UserKnmtcsCreateMVOrthoBarBar(CS -> C.Angle.Bar1, 
						       CS -> C.Angle.Bar2,
						       NumOfVars);
		break;
	    case USER_KNMTCS_CONSTR_DIST_PT_BAR:
		MVs[i] = UserKnmtcsCreateMVDstPtBar(CS -> C.DstPtBar.Pt, 
						    CS -> C.DstPtBar.Bar,
						    CS -> V.distance,
						    NumOfVars);
		break;
	    case USER_KNMTCS_CONSTR_TANGENCY:
	        {
		    MvarMVStruct **MVTan;

		    if (CS -> C.Tan.Pt -> Type ==
			                    USER_KNMTCS_PT_MOVES_ALONG_CURVE) {
			MVTan = UserKnmtcsCreateMVTanBarCrv(CS -> C.Tan.Pt,
							    CS -> C.Tan.Bar,
							    NumOfVars);
			MVs[i] = MVTan[0];
			MVs[i + 1] = MVTan[1];

			IritFree(MVTan);

			i++;
		    }
		}
		break;
	    case USER_KNMTCS_CONSTR_PARALLEL:
		MVs[i] = UserKnmtcsCreateMVParBarBar(CS -> C.Par.Bar1,
						     CS -> C.Par.Bar2,
						     NumOfVars);
		break;
	    case USER_KNMTCS_CONSTR_DIST_PT_CRV:
	        {
		    MvarMVStruct
		        **MVDist = UserKnmtcsCreateMVDstPtCrv(
						      CS -> C.DstPtCrv.CrvPt, 
						      CS -> C.DstPtCrv.Pt,
						      CS -> V.distance,
						      NumOfVars);

		    MVs[i] = MVDist[0];
		    MVs[i + 1] = MVDist[1];

		    IritFree(MVDist);

		    i++;
		}
	        break;
	    case USER_KNMTCS_CONSTR_DIST_PT_SRF:
	        {
		    MvarMVStruct
		        **MVDist = UserKnmtcsCreateMVDstPtSrf(
						      CS -> C.DstPtSrf.SrfPt,
						      CS -> C.DstPtSrf.Pt,
						      CS -> V.distance,
						      NumOfVars);

		    MVs[i] = MVDist[0];
		    MVs[i + 1] = MVDist[1];
		    MVs[i + 2] = MVDist[2];

		    IritFree(MVDist);

		    i = i + 2;
		}
	        break;
	    case USER_KNMTCS_CONSTR_MIN_DIST_PT_PT:/* Inequality constraint.*/
		MVs[i] = UserKnmtcsCreateMVDstPtPt(CS -> C.DstPtPt.Pt1, 
						   CS -> C.DstPtPt.Pt2,
						   CS -> V.distance,
						   NumOfVars);
		Constraints[i] = MVAR_CNSTRNT_POSITIVE;
		break; 
	    case USER_KNMTCS_CONSTR_MAX_DIST_PT_PT:/* Inequality constraint.*/
		MVs[i] = UserKnmtcsCreateMVDstPtPt(CS -> C.DstPtPt.Pt1, 
						   CS -> C.DstPtPt.Pt2,
						   CS -> V.distance,
						   NumOfVars);
		Constraints[i] = MVAR_CNSTRNT_NEGATIVE;
		break; 
	    case USER_KNMTCS_CONSTR_XDIFF_POS:     /* Inequality constraint.*/
		MVs[i] = UserKnmtcsCreateAxisDistFromBar(CS -> C.DstPtPt.Pt1, 
							 CS -> C.DstPtPt.Pt2,
							 CS -> V.distance,
							 NumOfVars, 1);
		Constraints[i] = MVAR_CNSTRNT_POSITIVE;
		break; 
	    case USER_KNMTCS_CONSTR_YDIFF_POS:     /* Inequality constraint.*/
		MVs[i] = UserKnmtcsCreateAxisDistFromBar(CS -> C.DstPtPt.Pt1, 
							 CS -> C.DstPtPt.Pt2,
							 CS -> V.distance,
							 NumOfVars, 2);
		Constraints[i] = MVAR_CNSTRNT_POSITIVE;
		break; 
	    case USER_KNMTCS_CONSTR_ZDIFF_POS:     /* Inequality constraint.*/
		MVs[i] = UserKnmtcsCreateAxisDistFromBar(CS -> C.DstPtPt.Pt1, 
							 CS -> C.DstPtPt.Pt2,
							 CS -> V.distance,
							 NumOfVars, 3);
		Constraints[i] = MVAR_CNSTRNT_POSITIVE;
		break; 
	    default:
		break;
	}
    }

    /* Solution of previous mechanism is freed.*/
    if (SolPolys != NULL) {
	UserKnmtcsFreeSol();
	SolPolys = NULL;
    }

#ifdef DEBUG_KINEMATIC_CONSTRAINTS
    for (i = 0; i < NumOfEqns + NumofInequalities; i++) {
	fprintf(stderr, "************************* Eqn %d *************************\n", i);
	MvarDbg(MVs[i]);
    }
#endif /* DEBUG_KINEMATIC_CONSTRAINTS */

    /* Well constrained system -> MvarZeros. */
    if (NumOfEqns == NumOfVars) {
	IsolPts = MvarMVsZeros0D(MVs, Constraints,
				 NumOfEqns + NumofInequalities,
				 SubTol, NumTol);

	while (IsolPts != NULL) {
	    IRIT_LIST_POP(Pt, IsolPts);
	    NewPoly = MvarPolylineNew(Pt);
	    IRIT_LIST_PUSH(NewPoly, SolPolys);	    
	}
	NumOfSols = CagdListLength(SolPolys);

    }
    /* Univariate solution space -> MvarMVUnivarInter. */
    else if (NumOfEqns == NumOfVars - 1) {
	SolPolys = MvarMVsZeros1D(MVs, Constraints,
			          NumOfEqns + NumofInequalities,
				  Step, SubTol, NumTol);
	NumOfSols = CagdListLength(SolPolys);
    }
    else {
	/* The mechanism is defined incorrectly. */
	NumOfSols = -1;
    }

#ifdef DEBUG
    /* Verify the solution. */
    if (NumOfEqns == NumOfVars - 1) {
        int i, j;
	MvarPolylineStruct *SolPl;
	MvarPtStruct *SolPt;

	for (SolPl = SolPolys; SolPl != NULL; SolPl = SolPl -> Pnext) {
	    for (SolPt = SolPl -> Pl; SolPt != NULL; SolPt = SolPt -> Pnext) {
	        for (i = 0; i < NumOfEqns + NumofInequalities; i++) {
		    CagdRType
		        *R = MvarMVEval(MVs[i], SolPt -> Pt);

		    switch (Constraints[i]) {
		        case MVAR_CNSTRNT_ZERO:
			    if (!IRIT_APX_EQ_EPS(R[1], 0.0, NumTol)) {
			        printf("Invalid (zero) solution! [");
				for (j = 0; j < SolPt -> Dim; j++)
				    printf("%f%s",
					   SolPt -> Pt[j],
					   j == SolPt -> Dim - 1 ? "]" : ", ");
				printf("MV[%d] = %g\n", i, R[1]);
			    }
			    break;
		        case MVAR_CNSTRNT_POSITIVE:
			    if (R[1] < 0.0) {
			        printf("Invalid (positive) solution! [");
				for (j = 0; j < SolPt -> Dim; j++)
				    printf("%f%s",
					   SolPt -> Pt[j],
					   j == SolPt -> Dim - 1 ? "]" : ", ");
				printf("MV[%d] = %g\n", i, R[1]);
			    }
			    break;
		        case MVAR_CNSTRNT_NEGATIVE:
			    if (R[1] > 0.0) {
			        printf("Invalid (negative) solution! [");
				for (j = 0; j < SolPt -> Dim; j++)
				    printf("%f%s",
					   SolPt -> Pt[j],
					   j == SolPt -> Dim - 1 ? "]" : ", ");
				printf("MV[%d] = %g\n", i, R[1]);
			    }
			    break;
	                default:
			    break;
		    }
		}
	    }
	}
    }
#endif /* DEBUG */

    for (i = 0; i < NumOfEqns + NumofInequalities; i++)
	MvarMVFree(MVs[i]);

    IritFree(UserKnmtcsDomainMin);
    IritFree(UserKnmtcsDomainMax);
    IritFree(Constraints);

    IritFree(MVs);

    return NumOfSols;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   For a solution polyline at "PolyIdx"'s position, the number of solution  M
* points is returned.							     M
*									     *
* PARAMETERS:                                                                M
*   PolyIdx:		Polyline index.			 		     M
*									     *
* RETURN VALUE:                                                              M
*   int:		Number of points in the polyline.		     M
*									     *
* KEYWORDS:                                                                  M
*   UserKnmtcsNumOfSolPts		                                     M
*****************************************************************************/
int UserKnmtcsNumOfSolPts(int PolyIdx)
{
    MvarPolylineStruct *TempPoly;
    int i;

    assert(CagdListLength(SolPolys) > PolyIdx && PolyIdx >= 0);

    for (i = 0, TempPoly = SolPolys;
	 i < PolyIdx;
	 i++, TempPoly = TempPoly -> Pnext);

    return (int) CagdListLength(TempPoly -> Pl);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Evaluates the motions curves of all the points in the kinematic          M
* mechanism at the solution. 						     M
*									     *
* PARAMETERS:                                                                M
*   None				                                     M
*									     *
* RETURN VALUE:                                                              M
*   CagdCrvStruct *:  List of curves of all traces of all points in the      M
*		      mechanisms.					     M
*									     *
* KEYWORDS:                                                                  M
*   UserKnmtcsEvalCrvTraces		                                     M
*****************************************************************************/
CagdCrvStruct *UserKnmtcsEvalCrvTraces()
{
    MvarPtStruct *TempPt;
    MvarPolylineStruct *TempPoly;
    int i, k, n, m;
    UserKnmtcsPtStruct *MP;
    CagdCrvStruct *TraceCrvs, *Crv,
	*AllTracesCrvs = NULL;

    /* How many points (and hence traced curve) we have? */
    for (n = 0, MP = GlblSystem.Pts; MP != NULL; n++, MP = MP -> Pnext);

    for (k = 0, TempPoly = SolPolys;
	 TempPoly != NULL;
	 k++, TempPoly = TempPoly -> Pnext) {
	/* Compute the size of the curves. */
        for (m = 0, TempPt = TempPoly -> Pl;
	     TempPt != NULL;
	     m++, TempPt = TempPt -> Pnext);

        /* Allocate n curves of the proper length. */
        for (TraceCrvs = NULL, i = 0; i < n; i++) {
	    Crv = BspCrvNew(m, IRIT_MIN(m, 2), CAGD_PT_E3_TYPE);
	    BspKnotUniformOpen(m, IRIT_MIN(m, 2), Crv -> KnotVector);
	    AttrSetIntAttrib(&Crv -> Attr, "TraceIdx", k);
	    IRIT_LIST_PUSH(Crv, TraceCrvs);
	}
	TraceCrvs = CagdListReverse(TraceCrvs);

	/* Update the curves with the data. */
        for (i = 0, TempPt = TempPoly -> Pl;
	     TempPt != NULL;
	     i++, TempPt = TempPt -> Pnext) {
	    for (Crv = TraceCrvs, MP = GlblSystem.Pts;
		 MP != NULL;
		 Crv = Crv -> Pnext, MP = MP -> Pnext) {
	        UserKnmtcsEvalOneSolPt(TempPt, MP);
		Crv -> Points[1][i] = MP -> Pt[0];
		Crv -> Points[2][i] = MP -> Pt[1];
		Crv -> Points[3][i] = MP -> Pt[2];
	    }
	}

	AllTracesCrvs = CagdListAppend(AllTracesCrvs, TraceCrvs);
    }  

    return AllTracesCrvs;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Evaluates kinematic mechanism at given solution. Solution is stored as   M
* a linked list of polylines, each position is uniquely determined by        M
* indices of the polyline and the point inside this polyline		     M
*									     *
* PARAMETERS:                                                                M
*   PolyIdx:		Polyline index.			 		     M
*   PtIdx:		Point index.			 		     M
*									     *
* RETURN VALUE:                                                              M
*   void								     M
*									     *
* KEYWORDS:                                                                  M
*   UserKnmtcsEvalAtParams		                                     M
*****************************************************************************/
void UserKnmtcsEvalAtParams(int PolyIdx, int PtIdx)
{
    MvarPtStruct *TempPt;
    MvarPolylineStruct *TempPoly;
    int i, j;
    UserKnmtcsPtStruct *MP;

    for (i = 0, TempPoly = SolPolys;
	 i < PolyIdx;
	 i++, TempPoly = TempPoly -> Pnext);

    for (j = 0, TempPt = TempPoly -> Pl;
	 j < PtIdx;
	 j++, TempPt = TempPt -> Pnext);

    for (MP = GlblSystem.Pts; MP != NULL; MP = MP -> Pnext) {
        UserKnmtcsEvalOneSolPt(TempPt, MP);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Evaluates kinematic mechanism at given solution point in poly.           *
*									     *
* PARAMETERS:                                                                *
*   SolPt:		The given solution point.			     *
*   MP:                 The kinematic point to update.			     *
*									     *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void UserKnmtcsEvalOneSolPt(MvarPtStruct *SolPt,
				   UserKnmtcsPtStruct *MP)
{
    switch (MP -> Type) {
        case USER_KNMTCS_PT_XY_PLANE:
	    MP -> Pt[0] = SolPt -> Pt[MP -> Idx];
	    MP -> Pt[1] = SolPt -> Pt[MP -> Idx + 1];
	    MP -> Pt[2] = 0.0;
	    break;
        case USER_KNMTCS_PT_NONE:
	    MP -> Pt[0] = SolPt -> Pt[MP -> Idx];
	    MP -> Pt[1] = SolPt -> Pt[MP -> Idx + 1];
	    MP -> Pt[2] = SolPt -> Pt[MP -> Idx + 2];
	    break;
	case USER_KNMTCS_PT_Y_DIRECTION:
	    /* MP -> Pt[0] is given from user, y coordinate is       */
	    /* directly the solution parameter.			 */
	    MP -> Pt[1] = SolPt -> Pt[MP -> Idx];
	    break;
        case USER_KNMTCS_PT_X_DIRECTION:
	    /* MP -> Pt[1] is given from user, x coordinate is       */
	    /* directly the solution parameter.			 */
	    MP -> Pt[0] = SolPt -> Pt[MP -> Idx];
	    break;
        case USER_KNMTCS_PT_Z_DIRECTION:
	    /* MP -> Pt[2] is given from user, z coordinate is       */
	    /* directly the solution parameter.			 */
	    MP -> Pt[2] = SolPt -> Pt[MP -> Idx];
	    break;
        case USER_KNMTCS_PT_MOVES_ALONG_CURVE:
	{
	    CagdRType 
	        *Coord = CagdCrvEval(MP -> U.Crv,
				     SolPt -> Pt[MP -> Idx]);

	    if (CAGD_IS_RATIONAL_CRV(MP -> U.Crv)) {
	        MP -> Pt[0] = Coord[1] / Coord[0];
		MP -> Pt[1] = Coord[2] / Coord[0];
		MP -> Pt[2] = Coord[3] / Coord[0];
	    }
	    else {
	        MP -> Pt[0] = Coord[1];
		MP -> Pt[1] = Coord[2];
		MP -> Pt[2] = Coord[3];
	    } 
	    break;
	}
        case USER_KNMTCS_PT_MOVES_ALONG_ROT_CURVE:
	{
	    CagdRType 
	        *Coord = CagdCrvEval(MP -> U.Crv,
				     SolPt -> Pt[MP -> Idx]);
	    
	    if (CAGD_IS_RATIONAL_CRV(MP -> U.Crv)) {
	        MP -> Pt[0] = Coord[1] / Coord[0];
		MP -> Pt[1] = Coord[2] / Coord[0];
		MP -> Pt[2] = Coord[3] / Coord[0];
	    }
	    else {
	        MP -> Pt[0] = Coord[1];
		MP -> Pt[1] = Coord[2];
		MP -> Pt[2] = Coord[3];
	    }
	    break;
	}
        case USER_KNMTCS_PT_MOVES_ALONG_SURFACE:
	{
	    CagdRType 
	        *Coord = CagdSrfEval(MP -> U.Srf,
				     SolPt -> Pt[MP -> Idx], 
				     SolPt -> Pt[(MP -> Idx) + 1]);
	    
	    if (CAGD_IS_RATIONAL_SRF(MP -> U.Crv)) {
	        MP -> Pt[0] = Coord[1] / Coord[0];
		MP -> Pt[1] = Coord[2] / Coord[0];
		MP -> Pt[2] = Coord[3] / Coord[0];
	    }
	    else {
	        MP -> Pt[0] = Coord[1];
		MP -> Pt[1] = Coord[2];
		MP -> Pt[2] = Coord[3];
	    }
	    break;
	}
        case USER_KNMTCS_PT_FIXED:
	    /* Point is fixed, none of its coordinates are unknowns. */
	    break;
        default:
	    assert(0);	
	    break;    
    }    
}
