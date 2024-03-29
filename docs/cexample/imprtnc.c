/*****************************************************************************
*   Computes the importance of vertices and dump as vertex list + polygons   *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/allocate.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/ip_cnvrt.h"

IRIT_STATIC_DATA char
    *CtrlStr = "Imprtnc c%- h%- DFiles!*s";
IRIT_STATIC_DATA int GlblVrtxImportanceCount,
    GlblCrvtrInfoFlag = FALSE;
IRIT_STATIC_DATA IrtRType GlblVrtxImportanceVal;
IRIT_STATIC_DATA IPVertexStruct *GlblVrtxImportance;

static void DumpOneTraversedObject(IPObjectStruct *PObj, IrtHmgnMatType Mat);
static IPObjectStruct *SetCurvatureEstimates(IPObjectStruct *PObj);
static void DumpOneObjData(IPObjectStruct *PObj);
static void ProcessVertexImportance(IPVertexStruct *V1,
				    IPVertexStruct *V2,
				    IPPolygonStruct *Pl1,
				    IPPolygonStruct *Pl2);
static void GenPolyImportance(IPObjectStruct *PObj);

void main(int argc, char **argv)
{
    int NumFiles, Error,
	HelpFlag = FALSE;
    char **FileNames;
    IPObjectStruct *PObjects;
    IrtHmgnMatType CrntViewMat;

    if ((Error = GAGetArgs(argc, argv, CtrlStr, &GlblCrvtrInfoFlag,
			   &HelpFlag, &NumFiles, &FileNames)) != 0) {
	GAPrintErrMsg(Error);
	GAPrintHowTo(CtrlStr);
	exit(1);
    }

    if (HelpFlag) {
	fprintf(stderr, "This is Importance testing...\n");
	GAPrintHowTo(CtrlStr);
	exit(0);
    }

    /* Get the data files: */
    IPSetFlattenObjects(FALSE);
    if ((PObjects = IPGetDataFiles(FileNames, NumFiles, TRUE, FALSE)) == NULL)
	exit(1);
    PObjects = IPResolveInstances(PObjects);

    if (IPWasPrspMat)
	MatMultTwo4by4(CrntViewMat, IPViewMat, IPPrspMat);
    else
	IRIT_GEN_COPY(CrntViewMat, IPViewMat, sizeof(IrtHmgnMatType));

    /* Here some useful parameters to play with in tesselating freeforms: */
    IPFFCState.FineNess = 15; /* Resolution of tesselation, larger is finer. */
    IPFFCState.ComputeUV = TRUE;       /* Wants UV coordinates for textures. */
    IPFFCState.FourPerFlat = TRUE;/* 4 polygons per ~flat patch, 2 otherwise.*/
    IPFFCState.LinearOnePolyFlag = TRUE;   /* Linear srf generates one poly. */

    IPTraverseObjListHierarchy(PObjects, CrntViewMat, DumpOneTraversedObject);

    IPFreeObjectList(PObjects);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Update curvature property attribute to vertices of the given poly model. *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:   Poly model to udpate curvatrue info into.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:   A poly model, similar to PObj, but with curvature    *
*			information attached as attributes.                  *
*****************************************************************************/
static IPObjectStruct *SetCurvatureEstimates(IPObjectStruct *PObj)
{
    int OldCirc = IPSetPolyListCirc(TRUE);
    IPPolygonStruct *Pl;
    IPObjectStruct *PTmp1, *PTmp2;

    /* Convert to a regular polygonal model with triangles only. */
    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        if (IPVrtxListLen(Pl -> PVertex) != 3)
	    break;
    }
    PTmp2 = IPCopyObject(NULL, PObj, FALSE);
    GMVrtxListToCircOrLin(PTmp2 -> U.Pl, TRUE);
    if (Pl != NULL) {
        PTmp1 = GMConvertPolysToTriangles(PTmp2);
	IPFreeObject(PTmp2);
	PTmp2 = GMRegularizePolyModel(PTmp1, TRUE);
	IPFreeObject(PTmp1);
    }

    GMVrtxListToCircOrLin(PTmp2 -> U.Pl, FALSE);
    IPSetPolyListCirc(OldCirc);

    GMPlCrvtrSetCurvatureAttr(PTmp2 -> U.Pl, 1, TRUE);

    return PTmp2;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function of IPTraverseObjListHierarchy. Called on every non    *
* list object found in hierarchy.                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       Non list object to handle.                                   *
*   Mat:        Transformation matrix to apply to this object.               *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpOneTraversedObject(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    IPObjectStruct *PObjs;

    if (IP_IS_FFGEOM_OBJ(PObj))
        PObjs = IPConvertFreeForm(PObj, &IPFFCState);
    else
	PObjs = PObj;

    for (PObj = PObjs; PObj != NULL; PObj = PObj -> Pnext) {
	if (IP_IS_POLY_OBJ(PObj)) {
	    if (GlblCrvtrInfoFlag) {
	        IPObjectStruct
		    *PTri = SetCurvatureEstimates(PObj);

		DumpOneObjData(PTri);
		IPFreeObject(PTri);
	    }
	    else
	        DumpOneObjData(PObj);
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Convert the polygonal mesh into a vertex list with importance and the    *
* polygons as indices into this list.                                        *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:   A polygonal mesh to dump.                                        *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void DumpOneObjData(IPObjectStruct *PObj)
{
    IPPolyVrtxIdxStruct
	*PVIdx = IPCnvPolyToPolyVrtxIdxStruct(PObj, TRUE, 0);
    int	i,
	**Polygons = PVIdx -> Polygons;
    IPVertexStruct
	**Vertices = PVIdx -> Vertices;

    /* Compute importance for the vertices as "Imprt" attributes.  Note the */
    /* PVIdx data structure points on the original vertices so we are fine. */
    GenPolyImportance(PObj);

    /* Dump the vertices: */
    fprintf(stderr, "OBJECT \"%s\" - VERTICES:\n", PObj -> ObjName);
    for (i = 0; Vertices[i] != NULL; i++) {
        IrtRType R;
	float *uv;

	printf("%3d:  %6.3f %6.3f %6.3f", i,
	       Vertices[i] -> Coord[0], 
	       Vertices[i] -> Coord[1], 
	       Vertices[i] -> Coord[2]);
	if (IP_HAS_NORMAL_VRTX(Vertices[i]))
	    printf(" [%6.3f %6.3f %6.3f]",
		   Vertices[i] -> Normal[0], 
		   Vertices[i] -> Normal[1], 
		   Vertices[i] -> Normal[2]);
	if ((uv = AttrGetUVAttrib(Vertices[i] -> Attr, "uvvals")) != NULL)
	    printf(" {%6.3f %6.3f}", uv[0], uv[1]);

	R = AttrGetRealAttrib(Vertices[i] -> Attr, "Imprt");
	if (!IP_ATTR_IS_BAD_REAL(R))
	    printf(" (%6.3f)", R);

	if (GlblCrvtrInfoFlag) {
	    printf("\n\t<K1=%6.3f D1 = %s>\n\t<K2=%6.3f D2 = %s>",
		   AttrGetRealAttrib(Vertices[i] -> Attr, "K1Curv"),
		   AttrGetStrAttrib(Vertices[i] -> Attr, "D1"),
		   AttrGetRealAttrib(Vertices[i] -> Attr, "K2Curv"),
		   AttrGetStrAttrib(Vertices[i] -> Attr, "D2"));
	}

	printf("\n");
    }

    fprintf(stderr, "\nOBJECT \"%s\" - Vertices in polygons:\n",
	    PObj -> ObjName);
    for (i = 0; PVIdx -> PPolys[i] != NULL; i++) {
        IPPolyPtrStruct
	    *PPoly = PVIdx -> PPolys[i];

	printf("%3d: ", i);
	for ( ; PPoly != NULL; PPoly = PPoly -> Pnext) {
	    printf("%3d ", AttrGetIntAttrib(PPoly -> Poly -> Attr, "_PIdx"));
	}
	printf("\n");
    }

    /* Dump the polygons: */
    fprintf(stderr, "\nOBJECT \"%s\" - Polygons from vertices:\n",
	    PObj -> ObjName);
    for (i = 0; Polygons[i] != NULL; i++) {
	int *Pl = Polygons[i];

	printf("%3d: ", i);
	while (*Pl >= 0)
	    fprintf(stderr, " %5d", *Pl++);

	fprintf(stderr, "    -1\n");
    }

    IPPolyVrtxIdxFree(PVIdx);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Call back function for GMPolyAdjacncyVertex to process every edge of a   *
* given vertex.  The edge is provided as (V, V -> Pnext)		     *
*                                                                            *
* PARAMETERS:                                                                *
*   V1, V2:    Two vertices defining this edge.  Note the vertices are NOT   *
*	       necessarily chained together into a list.		     *
*   Pl1, Pl2:  The two polygons that share this edge.  The edge (V1, V2) is  *
*	       in both Pl1 and Pl2, with not necessarily the exact pointers  *
*	       IPVertexStruct of V1 and V2.				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void ProcessVertexImportance(IPVertexStruct *V1,
				    IPVertexStruct *V2,
				    IPPolygonStruct *Pl1,
				    IPPolygonStruct *Pl2)
{
    if (!IRIT_PT_APX_EQ_EPS(V1 -> Coord, GlblVrtxImportance -> Coord,
			    IRIT_EPS) &&
	!IRIT_PT_APX_EQ_EPS(V2 -> Coord, GlblVrtxImportance -> Coord,
			    IRIT_EPS))
        fprintf(stderr, "Edge does not match the given vertex - adj error!\n");

    if (Pl1 != NULL && Pl2 != NULL) {
        GlblVrtxImportanceCount++;
	GlblVrtxImportanceVal += acos(IRIT_DOT_PROD(Pl1 -> Plane,
						    Pl2 -> Plane));
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Computes importance for each vertex in the given polygonal mesh.         *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:         Polygonal model to process for importance.	             *
*                                                                            *
* RETURN VALUE:                                                              *
*   void			                                             *
*****************************************************************************/
static void GenPolyImportance(IPObjectStruct *PObj)
{
    VoidPtr
	PAdj = GMPolyAdjacncyGen(PObj, IRIT_EPS);
    IPPolygonStruct *Pl;

    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        int PlImportanceCount = 0;
        IrtRType
	    PlImportanceVal = 0.0;
        IPVertexStruct *V;

	for (V = Pl -> PVertex; V != NULL; V = V -> Pnext) {
	    GlblVrtxImportanceVal = 0.0;
	    GlblVrtxImportanceCount = 0;
	    GlblVrtxImportance = V;

	    GMPolyAdjacncyVertex(V, PAdj, ProcessVertexImportance);

	    if (GlblVrtxImportanceCount > 0) {
		GlblVrtxImportanceVal /= (GlblVrtxImportanceCount * M_PI);

		AttrSetRealAttrib(&V -> Attr, "Imprt", GlblVrtxImportanceVal);
	    }
	    else
	        fprintf(stderr, "Failed to compute Importance for vertex\n");
	}
    }

    GMPolyAdjacncyFree(PAdj);
}
