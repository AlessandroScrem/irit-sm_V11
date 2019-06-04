/******************************************************************************
* SrfCrvtr.c - curvature computation of curves and surfaces.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, June 03.					      *
******************************************************************************/

#include "inc_irit/triv_lib.h"
#include "inc_irit/mvar_lib.h"
#include "user_loc.h"

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes the umbilical points of a given surface.                          M
*   The umbilicals are computed as the zeros of the function C = H^2 - K.    M
* C is never nagative and hence we actually solve for dC/du = dC/dv = 0 and  M
* test the values of C there.						     M
*   Hence (We only consider numerator of C which is sufficient for zeros),   M
*      C = H^2 - K = (2FM - EN - GL)^2 - 4(LN - M^2)(EG - F^2).		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:       Surface to compute its umbilical points, if any.              M
*   SubTol:    Subdivision tolerance of computation.			     M
*   NumTol:    Numerical tolerance of computation.			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarPtStruct *:   A list of UV parameters of umbilical points, or NULL   M
*		      if none.						     M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbSrfFff, SymbSrfSff, SymbSrfMeanCurvatureSqr, SymbSrfGaussCurvature,  M
*   MvarMVsZeros							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   UserSrfUmbilicalPts, curvature                                           M
*****************************************************************************/
MvarPtStruct *UserSrfUmbilicalPts(const CagdSrfStruct *Srf,
				  CagdRType SubTol,
				  CagdRType NumTol)
{
    CagdSrfStruct *DuSrf, *DvSrf, *SNormal, *STmp1, *STmp2, *STmp3, *STmp4,
	*HNumer, *FffDeterminant, *SffDeterminant, *CSrf,
	*FffE, *FffF, *FffG, *SffL, *SffM, *SffN;
    MvarPtStruct *MVPts, *MVPt, *MVPtTmp,
	*RetList = NULL;
    MvarMVStruct *Finals[2];

    SymbSrfFff(Srf, &DuSrf, &DvSrf, &FffE, &FffF, &FffG);
    SymbSrfSff(DuSrf, DvSrf, &SffL, &SffM, &SffN, &SNormal);
    CagdSrfFree(DuSrf);
    CagdSrfFree(DvSrf);
    CagdSrfFree(SNormal);

    FffDeterminant = SymbSrfDeterminant2(FffE, FffF, FffF, FffG);
    SffDeterminant = SymbSrfDeterminant2(SffL, SffM, SffM, SffN);

    STmp1 = SymbSrfMult(FffE, SffN);
    STmp2 = SymbSrfMult(FffG, SffL);
    STmp3 = SymbSrfMult(FffF, SffM);
    STmp4 = SymbSrfScalarScale(STmp3, 2.0);
    CagdSrfFree(STmp3);
    STmp3 = SymbSrfAdd(STmp1, STmp2);
    CagdSrfFree(STmp1);
    CagdSrfFree(STmp2);
    HNumer = SymbSrfSub(STmp3, STmp4);
    CagdSrfFree(STmp3);
    CagdSrfFree(STmp4);

    CagdSrfFree(FffE);
    CagdSrfFree(FffF);
    CagdSrfFree(FffG);
    CagdSrfFree(SffL);
    CagdSrfFree(SffM);
    CagdSrfFree(SffN);

    /* Combine the numerator of H and Fff and Sff together. */
    STmp1 = SymbSrfMult(HNumer, HNumer);
    CagdSrfFree(HNumer);

    STmp2 = SymbSrfMult(SffDeterminant, FffDeterminant);
    CagdSrfFree(FffDeterminant);
    CagdSrfFree(SffDeterminant);
    STmp3 = SymbSrfScalarScale(STmp2, 4.0);
    CagdSrfFree(STmp2);

    CSrf = SymbSrfSub(STmp1, STmp3);
    CagdSrfFree(STmp1);
    CagdSrfFree(STmp3);

    /* Construct the two derivatives of the C function, now that we have it. */
    STmp1 = CagdSrfDerive(CSrf, CAGD_CONST_U_DIR);
    STmp2 = CagdSrfDerive(CSrf, CAGD_CONST_V_DIR);

    Finals[0] = MvarSrfToMV(STmp1);
    Finals[1] = MvarSrfToMV(STmp2);
	
    MVPts = MvarMVsZeros0D(Finals, NULL, 2, SubTol, NumTol);

    MvarMVFree(Finals[0]);
    MvarMVFree(Finals[1]);

    /* Evaluate all solution points and pick those that has zero CSrf value. */
    while (MVPts != NULL) {
        CagdRType *R;

        MVPt = MVPts;
        MVPts = MVPts -> Pnext;
	MVPt -> Pnext = NULL;
        
	R = CagdSrfEval(CSrf, MVPt -> Pt[0], MVPt -> Pt[1]);
	if ((CAGD_IS_RATIONAL_SRF(CSrf) ? R[1] / R[0] : R[1]) < NumTol) {
	    for (MVPtTmp = RetList; MVPtTmp != NULL; MVPtTmp = MVPtTmp -> Pnext) {
		if (IRIT_PT_APX_EQ_E2_EPS(MVPt -> Pt, MVPtTmp -> Pt, NumTol))
		    break;
	    }
	    if (MVPtTmp == NULL) {
		IRIT_LIST_PUSH(MVPt, RetList);
	    }
	    else
	        MvarPtFree(MVPt);
	}
	else
	    MvarPtFree(MVPt);
    }

    return RetList;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes all locations (typically curves) on input surface Srf that have M
* one princple curvature value equal to a fixed value k1.		     M
*   Because   K1 k2 = K    and   ((k1 + k2)/2)^2 = H^2,			     V
* we have     ((k1 + K / k1)/2)^2 - H^2 = 0.				     V
*   Both K and H^2 are rational so let K = KN / KD (numerator / denominator) M
* and similarly let H^2 = HN / HD.					     M
*   Then, we need to solve for the zeros of				     M
* (k1 / 2 + KN / (KD * k1 * 2))^2 - HN / HD = 0, or the zeros of 	     M
* (K1 * KD * k1 + KN)^2 * HD - HN * (KD * k1 * 2)^2 = 0.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   Srf:        To compute lines of one principle curvature equal to k1.     M
*   k1:         The prescribed principle curvature.		             M
*   Step:	Step size for curve tracing.				     M
*   SubdivTol:	The subdivision tolerance to use.			     M
*   NumericTol:	The numerical tolerance to use.				     M
*   Euclidean:  TRUE to return curvature lines curves in Euclidean space,    M
*               FALSE in parameter space.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IPObjectStruct *:  A list of piecewise linear curves in Srf where one    M
*                      principle curvature is equal to k1.                   M
*                                                                            *
* SEE ALSO:                                                                  M
*   SymbSrfGaussCurvature, SymbSrfMeanCurvatureSqr			     M 
*                                                                            *
* KEYWORDS:                                                                  M
*   UserSrfFixedCurvatureLines, curvature	                             M
*****************************************************************************/
IPObjectStruct *UserSrfFixedCurvatureLines(const CagdSrfStruct *Srf,
					   CagdRType k1,
					   CagdRType Step,
					   CagdRType SubdivTol,
					   CagdRType NumericTol,
					   int Euclidean)
{
    CagdSrfStruct *KD, *KN, *HD, *HN, *STmp1, *STmp2, *STmp3, *STmp4,
        *KSrf = SymbSrfGaussCurvature(Srf, FALSE),
        *H2Srf = SymbSrfMeanCurvatureSqr(Srf);
    MvarPolylineStruct *MVPlls;
    MvarMVStruct *MV;
    IPPolygonStruct *Pl;
    IPObjectStruct
        *CrvtrLines = NULL;

    SymbSrfSplitScalar(KSrf, &KD, &KN, &STmp1, &STmp2);
    assert(STmp1 == NULL && STmp2 == NULL);
    SymbSrfSplitScalar(H2Srf, &HD, &HN, &STmp1, &STmp2);
    assert(STmp1 == NULL && STmp2 == NULL);
    CagdSrfFree(KSrf);
    CagdSrfFree(H2Srf);
    
    STmp1 = SymbSrfScalarScale(KD, IRIT_SQR(k1));
    STmp2 = SymbSrfAdd(STmp1, KN);
    STmp3 = SymbSrfMult(STmp2, STmp2);
    STmp4 = SymbSrfMult(STmp3, HD);
    CagdSrfFree(STmp1);
    CagdSrfFree(STmp2);
    CagdSrfFree(STmp3);

    STmp1 = SymbSrfScalarScale(KD, k1 * 2);
    STmp2 = SymbSrfMult(STmp1, STmp1);
    CagdSrfFree(STmp1);
    STmp1 = SymbSrfMult(STmp2, HN);
    CagdSrfFree(STmp2);

    STmp3 = SymbSrfSub(STmp4, STmp1);
    CagdSrfFree(STmp1);
    CagdSrfFree(STmp4);

    CagdSrfFree(KN);
    CagdSrfFree(KD);
    CagdSrfFree(HN);
    CagdSrfFree(HD);

    /* Compute the zeros of STmp3. */
    MV = MvarSrfToMV(STmp3);

    CagdSrfFree(STmp3);

    MVPlls = MvarMVsZeros1D(&MV, NULL, 1, Step, SubdivTol, NumericTol);
    MvarMVFree(MV);
    if (MVPlls != NULL) {
        CrvtrLines = MvarCnvrtMVPolysToIritPolys2(MVPlls, TRUE);
	MvarPolylineFreeList(MVPlls);

	if (CrvtrLines != NULL) {
	    CagdRType
		PurgeTol = sqrt(IRIT_ABS(NumericTol));

            SymbEvalSrfCurvPrep((CagdSrfStruct *) Srf, TRUE);

	    for (Pl = CrvtrLines -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
	        IPVertexStruct
		    *PrevV = NULL,
		    *V = Pl -> PVertex;

		while (V != NULL) {
		    CagdRType *R, Evalk1, Evalk2;
		    CagdVType EvalD1, EvalD2;

		    SymbEvalSrfCurvature(Srf, V -> Coord[0], V -> Coord[1],
					 TRUE, &Evalk1, &Evalk2,
					 EvalD1, EvalD2);

		    if (!IRIT_APX_EQ_EPS(k1, Evalk1, PurgeTol) &&
			!IRIT_APX_EQ_EPS(k1, Evalk2, PurgeTol)) {
			/* This point has no principle curvature at the     */
			/* designated value.  Can happen if we found -k1    */
			/* instead of k1.  Purge such solutions.            */
			if (PrevV == NULL) {
			    V = V -> Pnext;
			    IPFreeVertex(Pl -> PVertex);
			    Pl -> PVertex = V;
			}
			else {
			    PrevV -> Pnext = V -> Pnext;
			    IPFreeVertex(V);
			    V = PrevV -> Pnext;
			}
		    }
		    else {				 /* A valid vertex. */
		        if (Euclidean) {/* Eval and map to Euclidean space. */
			    R  = CagdSrfEval(Srf, V -> Coord[0], V -> Coord[1]);
			    CagdCoerceToE3(V -> Coord, &R, -1, Srf -> PType);
			}
			PrevV = V;
			V = V -> Pnext;
		    }
		}
	    }

	    SymbEvalSrfCurvPrep((CagdSrfStruct *) Srf, FALSE);

	    /* Clean up zero length polylines, in place. */
	    GMCleanUpPolylineList(&CrvtrLines -> U.Pl, IRIT_UEPS);
	}
    }

    return CrvtrLines;
}
