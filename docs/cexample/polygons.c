/*****************************************************************************
* Read a surface and dump out a tesselated version of it.		     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/ip_cnvrt.h"
#include "inc_irit/allocate.h"

void main(int argc, char **argv)
{
    int Handler,
	FourPerFlat = TRUE, /* Settable parameters of IritSurface2Polygons. */
	FineNess = 20,
	ComputeUV = FALSE,
	ComputeNrml = FALSE,
	Optimal = FALSE;

    if ((Handler = IPOpenDataFile("-", TRUE, TRUE)) >= 0) {
	IPObjectStruct
	    *PObj = IPGetObjects(Handler);

	/* Done with file - close it. */
	IPCloseStream(Handler, TRUE);

	/* Process the surface into polygons. */
	if (IP_IS_SRF_OBJ(PObj)) {
	    IPPolygonStruct
	        *PPoly = IPSurface2Polygons(PObj -> U.Srfs, FourPerFlat,
					    FineNess, ComputeUV,
					    ComputeNrml, Optimal);
	    IPObjectStruct
	        *PObjPoly = IPGenPOLYObject(PPoly);

	    IPStdoutObject(PObjPoly, FALSE);

	    IPFreeObject(PObjPoly);
	}
	else
	    fprintf(stderr, "Read object is not a surface.\n");

	IPFreeObject(PObj);
    }
    else {
	fprintf(stderr, "Failed to read from stdin\n");
	exit(1);
    }

    exit(0);
}
