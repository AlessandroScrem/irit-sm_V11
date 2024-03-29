/*****************************************************************************
* Least squares fit a curve to random points.				     *
* We also see how to fork out a display device and communicate with it.      *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/misc_lib.h"

static char *CtrlStr =
    "Lst_Sqrs n%-#Pts!d d%-Degree!d f%-DOF!d p%-PrgmName!s h%-";

void main(int argc, char **argv)
{
    int i, Error, PrgmIO,
	NumOfPoints = 100,
	NumOfPtsFlag = FALSE,
	Degree = 3,
	DegreeFlag = FALSE,
	NumOfDOF = 10,
	NumOfDOFFlag = FALSE,
	PrgmFlag = FALSE,
	HelpFlag = FALSE;
    char *Err,
	*Program = getenv("IRIT_DISPLAY");

#ifdef __WINNT__
    if (Program == NULL)
	Program = "wntgdrvs -s-";
#endif /* __WINNT__ */
#ifdef __UNIX__
    if (Program == NULL)
	Program = "x11drvs -s-";
#endif /* __UNIX__ */

    if ((Error = GAGetArgs(argc, argv, CtrlStr,
			   &NumOfPtsFlag, &NumOfPoints,
			   &DegreeFlag, &Degree,
			   &NumOfDOFFlag, &NumOfDOF,
			   &PrgmFlag, &Program,
			   &HelpFlag)) != 0) {
	GAPrintErrMsg(Error);
	GAPrintHowTo(CtrlStr);
	exit(1);
    }

    if (HelpFlag) {
	GAPrintHowTo(CtrlStr);
	exit(0);
    }

    IPSocSrvrInit();            /* Initialize the listen socket for clients. */

    if ((PrgmIO = IPSocExecAndConnect(Program,
				      getenv("IRIT_BIN_IPC") != NULL)) >= 0) {
	char Line[IRIT_LINE_LEN];
	IPObjectStruct
	    *PClrObj = IPGenStrObject("command_", "clear", NULL);

	do {
	    CagdPtStruct
		*PtList = NULL;
	    IPPolygonStruct
		*PPoly = IPAllocPolygon(0, NULL, NULL);
	    CagdCrvStruct *Crv;
	    IPObjectStruct *PCrvObj, *PPolyObj;

	    for (i = 0; i < NumOfPoints; i++) {
		int j;
		IPVertexStruct *V;
		CagdPtStruct
		    *Pt = CagdPtNew();

		if (i == 0) {
		    for (j = 0; j < 3; j++)
			Pt -> Pt[j] = IritRandom(-1.0, 1.0);
		}
		else {
		    for (j = 0; j < 3; j++)
			Pt -> Pt[j] = PtList -> Pt[j] + IritRandom(-0.1, 0.1);
		}

		V = IPAllocVertex(0, NULL, PPoly -> PVertex);
		for (j = 0; j < 3; j++)
		    V -> Coord[j] = Pt -> Pt[j];
		PPoly -> PVertex = V;

		IRIT_LIST_PUSH(Pt, PtList);
	    }

	    Crv = BspCrvInterpPts(PtList, Degree + 1,
				  NumOfDOF, CAGD_UNIFORM_PARAM, FALSE);
	    CagdPtFreeList(PtList);

	    CagdCrvWriteToFile3(Crv, stdout, 0, "This is from LstSqrs", &Err);

	    /* Generate objects out of the geometry and set proper attrs. */
	    PCrvObj = IPGenCRVObject(Crv);
	    AttrSetObjectColor(PCrvObj, IG_IRIT_CYAN);

	    PPolyObj = IPGenPOLYObject(PPoly);
	    IP_SET_POLYLINE_OBJ(PPolyObj);
	    AttrSetObjectColor(PPolyObj, IG_IRIT_YELLOW);

	    /* Clear old data and display our curve and data. */
	    IPSocWriteOneObject(PrgmIO, PClrObj);
	    IPSocWriteOneObject(PrgmIO, PCrvObj);
	    IPSocWriteOneObject(PrgmIO, PPolyObj);

	    IPFreeObject(PCrvObj);
	    IPFreeObject(PPolyObj);

	    gets(Line);
	}
	while (Line[0] != 'q' && Line[0] != 'Q');

	IPSocDisConnectAndKill(TRUE, PrgmIO);
    }

    exit(0);
}
