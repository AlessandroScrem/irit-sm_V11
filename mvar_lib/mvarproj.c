/******************************************************************************
* MvarProj.c - Compute the (orthogonal) projection of a curve on a surface.   *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber and Fady Massarwi, March 2011.		      *
******************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "mvar_loc.h"

#define MVAR_ORTH_PROJ_NUMER_TOL	1e-10

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the orthogonal projection of a curve C(t) on a surface S(u, v). M
*  That is, the projection is along the normal lines of the surface S.       M
*									     M
*  Computed as the univariate solution to the following two equations:	     M
*									     M
*                    dS		                 			     V
* < C(t) - S(u, v), ----- > = 0,				             V
*                    du				                             V
*									     M
*                    dS		                 			     V
* < C(t) - S(u, v), ----- > = 0.				             V
*                    dv				                             V
*                                                                            *
* PARAMETERS:                                                                M
*   Crv:          The curve to project on Srf, orthogonally.		     M
*   Srf:          The surface to project Crv on.			     M
*   Tol:          Tolerance of the computation.			             M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPolylineStruct *:    The projections in UV space of Srf.	     M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVOrthoIsoCrvProjOnSrf			 	   	             M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVOrthoCrvProjOnSrf, projection                                      M
*****************************************************************************/
MvarPolylineStruct *MvarMVOrthoCrvProjOnSrf(const CagdCrvStruct *Crv,
					    const CagdSrfStruct *Srf,
					    CagdRType Tol)
{
    CagdRType Min, Max;
    MvarMVStruct *MSrf, *MCrv, *MTmp, *DuMSrf, *DvMSrf, *MVs[2];
    MvarPolylineStruct *Solution;
    CagdSrfStruct *TSrf, *NewSrf;

    /* Convert the curve and surface to trivariates T(u, v, t). */

    if (CAGD_IS_BEZIER_SRF(Srf)) { /* Ensure building Bspline constraints...*/
        Srf = NewSrf = CagdCnvrtBzr2BspSrf(Srf);
    }
    else
        NewSrf = NULL;

    MTmp = MvarSrfToMV(Srf);
    MSrf = MvarPromoteMVToMV2(MTmp, 3, 0);  
    MvarMVFree(MTmp);

    TSrf = CagdSrfDerive(Srf, CAGD_CONST_U_DIR);
    MTmp = MvarSrfToMV(TSrf);  
    CagdSrfFree(TSrf);
    DuMSrf = MvarPromoteMVToMV2(MTmp, 3, 0);  
    MvarMVFree(MTmp);

    TSrf = CagdSrfDerive(Srf, CAGD_CONST_V_DIR);
    MTmp = MvarSrfToMV(TSrf);  
    CagdSrfFree(TSrf);
    DvMSrf = MvarPromoteMVToMV2(MTmp, 3, 0);  
    MvarMVFree(MTmp);

    /* Build the constraints. */
    if (CAGD_IS_BEZIER_CRV(Crv)) { /* Ensure building Bspline constraints...*/
        CagdCrvStruct
	    *TCrv = CagdCnvrtBzr2BspCrv(Crv);

	MTmp = MvarCrvToMV(TCrv);
	CagdCrvFree(TCrv);
    }
    else
        MTmp = MvarCrvToMV(Crv);
    MCrv = MvarPromoteMVToMV2(MTmp, 3, 2);
    MvarMVFree(MTmp);

    if (NewSrf != NULL)
	CagdSrfFree(NewSrf);

    /* Make domains the same. */
    MvarMVDomain(MSrf, &Min, &Max, 0);
    MTmp = MvarMVSetDomain(MCrv, Min, Max, 0, FALSE);
    MvarMVFree(MCrv);
    MCrv = MTmp;

    MvarMVDomain(MSrf, &Min, &Max, 1);
    MTmp = MvarMVSetDomain(MCrv, Min, Max, 1, FALSE);
    MvarMVFree(MCrv);
    MCrv = MTmp;

    MvarMVDomain(MCrv, &Min, &Max, 2);
    MTmp = MvarMVSetDomain(MSrf, Min, Max, 2, FALSE);
    MvarMVFree(MSrf);
    MSrf = MTmp;
    MTmp = MvarMVSetDomain(DuMSrf, Min, Max, 2, FALSE);
    MvarMVFree(DuMSrf);
    DuMSrf = MTmp;
    MTmp = MvarMVSetDomain(DvMSrf, Min, Max, 2, FALSE);
    MvarMVFree(DvMSrf);
    DvMSrf = MTmp;

    MTmp = MvarMVSub(MSrf, MCrv);
    MVs[0] = MvarMVDotProd(MTmp, DuMSrf);
    MVs[1] = MvarMVDotProd(MTmp, DvMSrf);
    MvarMVFree(MTmp);

    MvarMVFree(MCrv);
    MvarMVFree(MSrf);
    MvarMVFree(DuMSrf);
    MvarMVFree(DvMSrf);

    Solution = MvarMVsZeros1D(MVs, NULL, 2, Tol * 0.25, Tol,
			      MVAR_ORTH_PROJ_NUMER_TOL);

    MvarMVFree(MVs[0]);
    MvarMVFree(MVs[1]);

    return Solution;
}

/******************************************************************************
* DESCRIPTION:                                                                M
*   Computes the orthogonal projection of an isoparameteric curve of surface  M
* S1(r, t) at a fixed parameter value, RVal, into surface S2(u, v).           M
*   The projection is along the normal lines of S1, that contains the curve.  M
*									      M
*  Computed as the univariate solution to the following two equations:	      M
*									      M
*                         dS1                                                 V
* < S2(u, v) - S1(r, t), ----- > |          = 0,		              V
*                         dr     | (r=RVal)                                   V
*									      M
*                         dS1		                                      V
* < S2(u, v) - S1(r, t), ----- > |          = 0.		              V
*                         dt     | (r=RVal)                                   V
*                                                                             *
* PARAMETERS:                                                                 M
*   Srf1:         The surface to project from, orthogonaly.		      M
*   RVal:         The isoparametric value of the curve of Srf1 to project.    M
*   Dir:          Direction of isoparametric curve of S1.                     M
*   Srf2:         The surface to project to.                                  M
*   Tol:          Computation tolerance.                                      M
*                                                                             *
* RETURN VALUE:                                                               M
*   MvarPolylineStruct *:    The projections in UV space of Srf2.	      M
*                                                                             *
* SEE ALSO:                                                                   M
*   MvarMVOrthoCrvProjOnSrf				 		      M
*                                                                             *
* KEYWORDS:                                                                   M
*   MvarMVOrthoIsoCrvProjOnSrf, projection                                    M
******************************************************************************/
MvarPolylineStruct *MvarMVOrthoIsoCrvProjOnSrf(const CagdSrfStruct *Srf1,
					       const CagdRType RVal,
					       CagdSrfDirType Dir,
					       const CagdSrfStruct *Srf2,
					       CagdRType Tol)
{
    CagdRType Min, Max;
    MvarMVStruct *MSrf, *MCrv, *MTmp, *DuMSrf, *DvMSrf, *MVs[2];
    MvarPolylineStruct *Solution;
    CagdSrfStruct *TSrf, *NewSrf;
    CagdCrvStruct *Srf1KnotLine, *TCrv;
    
    /* Convert the curve and surface to trivariates T(u, v, t). */
    
    if (CAGD_IS_BEZIER_SRF(Srf2)) { /* Ensure building Bspline constraints...*/
        Srf2 = NewSrf = CagdCnvrtBzr2BspSrf(Srf2);
    }
    else
        NewSrf = NULL;
    
    MTmp = MvarSrfToMV(Srf2);
    MSrf = MvarPromoteMVToMV2(MTmp, 3, 0);
    MvarMVFree(MTmp);
    
    TSrf = CagdSrfDerive(Srf1, CAGD_CONST_U_DIR);
    TCrv =CagdCrvFromSrf(TSrf, RVal, Dir);
    CagdSrfFree(TSrf);
    MTmp = MvarCrvToMV(TCrv);
    CagdCrvFree(TCrv);
    DuMSrf = MvarPromoteMVToMV2(MTmp, 3, 2);
    MvarMVFree(MTmp);
    
    TSrf = CagdSrfDerive(Srf1, CAGD_CONST_V_DIR);
    TCrv =CagdCrvFromSrf(TSrf, RVal, Dir);
    CagdSrfFree(TSrf);
    MTmp = MvarCrvToMV(TCrv);
    CagdCrvFree(TCrv);
    DvMSrf = MvarPromoteMVToMV2(MTmp, 3, 2);
    MvarMVFree(MTmp);
    
    /* Build the constraints. */
    Srf1KnotLine = CagdCrvFromSrf(Srf1, RVal, Dir);
    if (CAGD_IS_BEZIER_CRV(Srf1KnotLine)) { /* Might not be relevant .... */
        CagdCrvStruct
        *TCrv = CagdCnvrtBzr2BspCrv(Srf1KnotLine);
        
        MTmp = MvarCrvToMV(TCrv);
        CagdCrvFree(TCrv);
    }
    else
        MTmp = MvarCrvToMV(Srf1KnotLine);
    
    CagdCrvFree(Srf1KnotLine);
    
    MCrv = MvarPromoteMVToMV2(MTmp, 3, 2);
    MvarMVFree(MTmp);
    
    if (NewSrf != NULL)
        CagdSrfFree(NewSrf);
    
    /* Make domains the same. Set u,v from Srf2 to knot line and derivatives  */
    /* of Srf1.                                                               */
    MvarMVDomain(MSrf, &Min, &Max, 0);
    MTmp = MvarMVSetDomain(MCrv, Min, Max, 0, FALSE);
    MvarMVFree(MCrv);
    MCrv = MTmp;
    
    MTmp = MvarMVSetDomain(DuMSrf, Min, Max, 0, FALSE);
    MvarMVFree(DuMSrf);
    DuMSrf = MTmp;
    
    MTmp = MvarMVSetDomain(DvMSrf, Min, Max, 0, FALSE);
    MvarMVFree(DvMSrf);
    DvMSrf = MTmp;
    
    MvarMVDomain(MSrf, &Min, &Max, 1);
    MTmp = MvarMVSetDomain(MCrv, Min, Max, 1, FALSE);
    MvarMVFree(MCrv);
    MCrv = MTmp;
    
    MTmp = MvarMVSetDomain(DuMSrf, Min, Max, 1, FALSE);
    MvarMVFree(DuMSrf);
    DuMSrf = MTmp;
    
    MTmp = MvarMVSetDomain(DvMSrf, Min, Max, 1, FALSE);
    MvarMVFree(DvMSrf);
    DvMSrf = MTmp;
    
    
    MvarMVDomain(MCrv, &Min, &Max, 2);
    MTmp = MvarMVSetDomain(MSrf, Min, Max, 2, FALSE);
    MvarMVFree(MSrf);
    MSrf = MTmp;
    
    MTmp = MvarMVSub(MSrf, MCrv);
    MVs[0] = MvarMVDotProd(MTmp, DuMSrf);
    MVs[1] = MvarMVDotProd(MTmp, DvMSrf);
    MvarMVFree(MTmp);
    
    MvarMVFree(MCrv);
    MvarMVFree(MSrf);
    MvarMVFree(DuMSrf);
    MvarMVFree(DvMSrf);
    
    Solution = MvarMVsZeros1D(MVs, NULL, 2, Tol * 0.25, Tol,
                              MVAR_ORTH_PROJ_NUMER_TOL);
    
    MvarMVFree(MVs[0]);
    MvarMVFree(MVs[1]);
    
    return Solution;
}
