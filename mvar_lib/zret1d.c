/******************************************************************************
* ZrET1D.c - tools to compute zero sets of multivariates when the problem has *
*	     Expression Trees representation and the expected solution  set   *
*	     is a one-manifold: curve(s).				      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Yoni Mizrahi, Feb. 13.					      *
******************************************************************************/

#include "mvar_loc.h"
#include "inc_irit/misc_lib.h"

static CagdBType MvarZeroSolverNoZeroTest1DExpTr(MvarZeroPrblmStruct *Problem,
						 int Ind);
static void MvarZeroSolveBoundary1DExpTr(MvarZeroPrblmStruct *Problem);
static CagdBType MvarZeroSolverTopoTest1DExpTr(MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep1DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol1DExpTr(
					        MvarZeroPrblmStruct *Problem);
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions1DExpTr(
						MvarZeroSolutionStruct *Sol1,
						MvarZeroSolutionStruct *Sol2, 
					        MvarMVDirType Dir,
					        CagdRType Param,
						MvarZeroPrblmStruct *Problem);
static MvarZeroPrblmStruct **MvarZeroSubdivProblem1DExpTr(
						 MvarZeroPrblmStruct *Problem);

IRIT_STATIC_DATA const MvarZeroSolverCallBackFcnStruct 
    ZrET1DCB = {
        MvarZeroSolveBoundary1DExpTr,
        MvarZeroSolverNoZeroTest1DExpTr,
        MvarZeroSolverTopoTest1DExpTr,
	MvarZeroSubdivProblem1DExpTr,
	MvarZeroSolverNumericStep1DExpTr,
	MvarZeroSolverSinglrSol1DExpTr,
	MvarZeroSolverUniteSolutions1DExpTr,
	NULL,
	MvarZeroUpdateProblemDmnExpTr,
	MvarZeroFirstSmoothUpdatesExpTr
    };

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets the callback functions for the ETs representation, 1D solution case.M
*									     *
* PARAMETERS:								     M
*   Problem:	The zero finding problem to be solved.			     M
*		                                                             *
* RETURN VALUE: 							     M
*   void								     M
*		                                                             *
* SEE ALSO:								     M
*   									     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarZeroSolverSetCallbackFcns1DExpTr				     M
*****************************************************************************/
void MvarZeroSolverSetCallbackFcns1DExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    IRIT_GEN_COPY(&PrblmIntrnl -> CallbackFunctions, &ZrET1DCB,
		  sizeof(MvarZeroSolverCallBackFcnStruct));
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Testing if we can rule out the possibility of a solution for the problem   *
* in its domain, by inspecting constraint number Ind. Callback function      *
* for the ETs case, 1D solution. 					     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*   Ind:	The index of the ET to be tested.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if domain can be purged (no zero can exist), FALSE      *
*		otherwise.						     *
*****************************************************************************/
static CagdBType MvarZeroSolverNoZeroTest1DExpTr(MvarZeroPrblmStruct *Problem,
						 int Ind)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Solving the ETs representation problem on the boundary of the domain. In *
*   the 1D case, it means we are solving 2 * Dim fully constrained problems. *
*									     *
* PARAMETERS:								     *
*   Problem:	The zero finding problem to be solved.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void MvarZeroSolveBoundary1DExpTr(MvarZeroPrblmStruct *Problem)
{
    MvarZeroPrblmIntrnlStruct
	*PrblmIntrnl = Problem -> _Internal;

    /* EMPTY DEFINITION - DEFAULT ONLY. */
    PrblmIntrnl -> SolutionsOnBoundary = NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Topological guarantee in the ETs, 1D solution case.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE if the solution curve is homeomorphic to an interval,   *
*		FALSE if cannot guarantee.				     *
*****************************************************************************/
static CagdBType MvarZeroSolverTopoTest1DExpTr(MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return FALSE;    
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   The numeric improvement for the ETs representation problem, 1D case.     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution, or NULL upon failure of the     *
*			       process.					     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverNumericStep1DExpTr(
					         MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Returning a solution in the case of reaching subdivision tolerance. In the *
* ETs, 1D case we do TBD.						     *
*                                                                            *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct *:  The solution.				     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverSinglrSol1DExpTr(
					         MvarZeroPrblmStruct *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY */
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Unites two neighboring 1D solutions of a zero finding problem, ETs	     *
* representation, back to a new solution of the subdivided domain:           *
* the polylines are linked, making sure that in the common boundary we have  *
* identical points, so the linkage is consistent.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Sol1:	The first solution.					     *
*   Sol2:	The second solution.					     *
*   Dir:	The direction along which the subdivision occured.	     *
*   Param:	The paramter value at which the subdivison occured, along    *
*		direction Dir.						     *
*   Problem:	The zero finding problems that are solved by Sol1 and Sol2,  *
*		before it was divided.					     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroSolutionStruct:  The new solution, union of Sol1 and Sol2.	     *
*****************************************************************************/
static MvarZeroSolutionStruct *MvarZeroSolverUniteSolutions1DExpTr(
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
* solver, in the 1D, ETs case.	  Given a problem, returns a vector of two   *
* subdivided problems, allocated dynamically.  The subdivided problems can   *
* be NULL.								     *
*									     *
* PARAMETERS:                                                                *
*   Problem:	The zero finding problem structure.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   MvarZeroPrblmStruct **:  The array of two sub-problems.		     *
*****************************************************************************/
static MvarZeroPrblmStruct **MvarZeroSubdivProblem1DExpTr(MvarZeroPrblmStruct
							            *Problem)
{
    /* EMPTY DEFINITION - DEFAULT ONLY. */
    return NULL;
}
