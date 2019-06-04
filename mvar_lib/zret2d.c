/******************************************************************************
* ZrET2D.c - tools to compute zero sets of multivariates when the problem has *
*	     Expression Trees representation and the expected solution  set   *
*	     is a two-manifold: surface(s).				      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, Feb. 13.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"

static CagdBType MvarZeroSolverNoZeroTest2DExpTr(MvarZeroPrblmStruct *Problem,
						 int Ind);
static void MvarZeroSolveBoundary2DExpTr(MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverTopoTest2DExpTr(MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep2DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol2DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions2DExpTr(
						MvarZeroSolutionStruct *Sol1,
						MvarZeroSolutionStruct *Sol2, 
					        MvarMVDirType Dir,
					        CagdRType Param,
						MvarZeroPrblmStruct *Problem);
static MvarZeroPrblmStruct **MvarZeroSubdivProblem2DExpTr(
						 MvarZeroPrblmStruct *Problem);

IRIT_STATIC_DATA const MvarZeroSolverCallBackFcnStruct 
    ZrET2DCB = {
        MvarZeroSolveBoundary2DExpTr,
        MvarZeroSolverNoZeroTest2DExpTr,
        MvarZeroSolverTopoTest2DExpTr,
	MvarZeroSubdivProblem2DExpTr,
	MvarZeroSolverNumericStep2DExpTr,
	MvarZeroSolverSinglrSol2DExpTr,
	MvarZeroSolverUniteSolutions2DExpTr,
	NULL,
	MvarZeroUpdateProblemDmnExpTr,
	MvarZeroFirstSmoothUpdatesExpTr
    };

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the callback functions for the ETs representation, 2D solution case.M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*		                                                             *
* RETURN VALUE:								     M
*   void								     M
*		                                                             *
* SEE ALSO:								     M
*   									     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverSetCallbackFcns2DExpTr				     M
*****************************************************************************/
void MvarZeroSolverSetCallbackFcns2DExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    IRIT_GEN_COPY(&PrblmIntrnl -> CallbackFunctions, &ZrET2DCB,
		  sizeof(MvarZeroSolverCallBackFcnStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Testing if we can rule out the possibility of a solution for the problem   *
* in its domain, by inspecting constraint number Ind. Callback function      *
* for the ETs case, 2D solution. 					     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*   Ind:	The index of the ET to be tested.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if domain can be purged (no zero can exist), FALSE      *
*		otherwise.						     *
*****************************************************************************/
static CagdBType MvarZeroSolverNoZeroTest2DExpTr(MvarZeroPrblmStruct *Problem,
						 int Ind)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Solving the ETs representation problem on the boundary of the domain. In *
*   the 2D case, it means we are solving 2 * Dim under constrained problems, *
*   with one degree of freedom, namely a 1D solution curve.		     *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*									     *
* KEYWORDS:								     *
*   MvarZeroSolveBoundary2DExpTr					     *
*****************************************************************************/
static void MvarZeroSolveBoundary2DExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    /* EMPTY DEFINITION - DEFAULT ONLY */
    PrblmIntrnl -> SolutionsOnBoundary = NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Topological guarantee in the ETs, 2D solution case.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if the solution curve is homeomorphic to a 2D disc,     *
*		FALSE if cannot guarantee.				     *
*****************************************************************************/
static CagdBType MvarZeroSolverTopoTest2DExpTr(MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return FALSE;    
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   The numeric improvement for the ETs representation problem, 2D case.     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution, or NULL upon failure of the     *
*			       process.					     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep2DExpTr(
					         MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Returning a solution in the case of reaching subdivision tolerance.        *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution.				     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol2DExpTr(
					         MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Unites two neighboring 2D solutions of a zero finding problem, ETs	     *
* representation, back to a new solution of the subdivided domain:           *
* the polylines are linked, making sure that in the common boundary we have  *
* identical curves, so the linkage is consistent.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol1:	The first solution.					     *
*   Sol2:	The second solution.					     *
*   Dir:	The direction along which the subdivision occured.	     *
*   Param:	The paramter value at which the subdivison occured, along    *
*		direction Dir.						     *
*   Problem:	The zero finding problems that are solved by Sol1	     *
*		and Sol2, before it was divided.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The new solution, union of Sol1 and Sol2.     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions2DExpTr(
						 MvarZeroSolutionStruct *Sol1,
						 MvarZeroSolutionStruct *Sol2, 
					         MvarMVDirType Dir,
					         CagdRType Param,
						 MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* The recursive stage of the subdivision solver, invoked by the general      *
* solver, in the 2D, ETs case.  Given a problem, returns a vector of two     *
* subdivided problems, allocated dynamically.  The subdivided problems can   *
* be NULL.								     *
*									     *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroPrblmStruct **:  The array of two sub-problems.		     *
*****************************************************************************/
static MvarZeroPrblmStruct **MvarZeroSubdivProblem2DExpTr(
						 MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY. */
    return NULL;
}
