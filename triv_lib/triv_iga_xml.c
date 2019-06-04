/*****************************************************************************
* Module to support IGA (isogeometric analysis) XML parsing.                 *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber and Fady Massarwi		Ver 1.0, July 2013   *
*****************************************************************************/

#include "triv_loc.h"
#include "triv_loc_iga.h"

#ifdef IRIT_HAVE_XML_LIB

#include "roxml.h"

static void TrivIGANeighboringConstraintXML(struct TrivIGAArrangementStruct *H,
					    CagdBType IsRational,
					    int MaxCoord,
					    const int *FirstIDs,
					    int FirstIDsSize,
					    int SecondID,
					    int SecondIdx,
					    const CagdSrfStruct **Mat,
					    void *CallBackData);
static void TrivIGAAddMaterialsToXML(struct TrivIGAArrangementStruct *H,
				     node_t *Root);
static void TrivIGAAddGeometryToXML(struct TrivIGAArrangementStruct *H,
				    node_t *Root);
static void TrivIGAAddBoundaryToXML(struct TrivIGAArrangementStruct *H,
				    node_t *Root);
static CagdBType TrivIGAIsPointOfFieldType(struct TrivIGAArrangementStruct *H, 
					   int Id,
					   TrivIGAFieldType FieldType);
static CagdRType *TrivIGAGetUniquePoint(struct TrivIGAArrangementStruct *H, 
					int Id) ;
static int TrivIGAFieldCompFunc(const VoidPtr VPt1, const VoidPtr VPt2);


/*****************************************************************************
* DESCRIPTION:                                                               *
*   Syntehsizes one set of constraints for a control point on one srf face   *
* as a function of control points in the other srf face, for all axes.       *
*                                                                            *
*                                                                            *
* PARAMETERS:                                                                *
*   H:            A handle on the IGA arrangement.   		             *
*   IsRational:   TRUE if data is rational so apply for weights too.         *
*   MaxCoord:     Maximal coordinates to handle.  I.e. 3 for XYZ.	     *
*   FirstIDs:     IDs of control point in original srf.                      *
*   FirstIDsSize: Size of vector FirstIDs.				     *
*   SecondID:     ID of one control points of refined/deg. raised srf.       *
*   SecondIdx:    Index of one control points of refined/deg. raised srf.    *
*   Mat:          All values of coefficients of constraints in refined or    *
*                 deg. raised srfs.  Matrix is a vector of E1 surfaces.      *
*                   The length of this vector is the same as FirstIDsSize.   *
*   CallbackData: Relevant data for the callback function		     *
*		  (such as XML information).				     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void TrivIGANeighboringConstraintXML(struct TrivIGAArrangementStruct *H,
					    CagdBType IsRational,
					    int MaxCoord,
					    const int *FirstIDs,
					    int FirstIDsSize,
					    int SecondID,
					    int SecondIdx,
					    const CagdSrfStruct **Mat,
					    void *CallBackData)
{
    IRIT_STATIC_DATA char TmpBuff[TRIV_IGA_MAX_XML_ROW_LEN];
    int i, k;
    node_t
        *XmlParentNode = (node_t *)CallBackData,
        *CurrentNode, *CurrentParentNode, *ContactNode;

    for (k = !IsRational; k <= MaxCoord; k++) {
	char 
	    Axis = k == 0 ? 'w' : k + 'x' - 1;
#	ifdef DEBUG
	IrtRType Sum;
#	endif /* DEBUG */

	ContactNode = roxml_add_node(XmlParentNode, 0, ROXML_ELM_NODE, 
	                             "contact", NULL);
	roxml_add_node(ContactNode, 0, ROXML_ATTR_NODE, "type", 
	               "linear constraint");

	CurrentParentNode = roxml_add_node(ContactNode, 0, ROXML_ELM_NODE, 
	                                   "linear_constraint", NULL);

	CurrentNode = roxml_add_node(CurrentParentNode, 0, ROXML_ELM_NODE, 
	                             "node", "-1.0");
	sprintf(TmpBuff, "%d", SecondID + 1);         /* ID starts at one. */
	roxml_add_node(CurrentNode, 0, ROXML_ATTR_NODE, "id", TmpBuff);
	sprintf(TmpBuff, "%c", Axis);
	roxml_add_node(CurrentNode, 0, ROXML_ATTR_NODE, "bc", TmpBuff);
        
	for (i = 0; i < FirstIDsSize; i++) {
	    if (Mat[i] != NULL &&
		IRIT_ABS(Mat[i] -> Points[1][SecondIdx]) > IRIT_UEPS) {

		sprintf(TmpBuff, "%.13g", Mat[i] -> Points[1][SecondIdx]);
		CurrentNode = roxml_add_node(CurrentParentNode, 0, 
					     ROXML_ELM_NODE, "node", TmpBuff);
		sprintf(TmpBuff, "%d", FirstIDs[i] + 1);
		roxml_add_node(CurrentNode, 0, ROXML_ATTR_NODE, "id", 
		               TmpBuff);
		sprintf(TmpBuff, "%c", Axis);
		roxml_add_node(CurrentNode, 0, ROXML_ATTR_NODE, "bc", 
		               TmpBuff);
	    }
#	    ifdef DEBUG
	    {
	        IrtRType *R;

		/* Verify that the constraint sums to zero: */
		if (i == 0) {
		    R = TrivIGAGetUniquePoint(H, SecondID);
		    Sum = -R[k];
		}
		if (Mat[i] != NULL &&
		    IRIT_ABS(Mat[i] -> Points[1][SecondIdx]) > IRIT_UEPS) {
		    R = TrivIGAGetUniquePoint(H, FirstIDs[i]);
		    Sum += R[k] * Mat[i] -> Points[1][SecondIdx];
		}
		if (i == FirstIDsSize - 1) {
		    if (!IRIT_APX_EQ(Sum, 0.0))
		        assert(0);
		}
	    }
#	    endif /* DEBUG */
	}

	/* Add tolerance and penalty nodes, currently set to defaults as    */
	/* in the XML example file.					    */
	roxml_add_node(ContactNode, 0, ROXML_ELM_NODE, 
		       "tol", "1.0e-04");
	roxml_add_node(ContactNode, 0, ROXML_ELM_NODE, 
		       "penalty", "1.0e+03");
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Load the materials defined in the given XML file.                        M
*                                                                            *
* PARAMETERS:                                                                M
*   FileName:        XML file name.			                     M
*   Materials:       Array of loaded materials allocated dynamically.        M
*   NumMaterials:    The number of the loaded materials.                     M
*									     *
* RETURN VALUE:                                                              M
*   int:     TRUE on success, and FALSE in failure                           M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGALoadMaterialXML                                                   M
*****************************************************************************/
int TrivIGALoadMaterialXML(const char *FileName,
			   TrivIGAMaterialStruct **Materials,
			   int *NumMaterials)
{
    char TempStrBuffer[TRIV_IGA_MAX_MATERIAL_NAME_LEN];
    int i, j, TempStrBufferSize;
    node_t 
	*Root = NULL,
	*MaterialNode = NULL,
	*InternalMaterialNode = NULL,
	*InternalMaterialAttribNode = NULL,
	*InternalMaterialChildNode = NULL;

    if (Materials == NULL ||
	NumMaterials == NULL || 
	FileName == NULL) {
	return 0;
    }

    *Materials = NULL;
    *NumMaterials = 0;

    Root = roxml_load_doc((char *)FileName);

    if (Root == NULL)
	return FALSE;

    MaterialNode = roxml_get_chld(Root, "Material", 0);

    if (MaterialNode == NULL)
	return FALSE;

    *NumMaterials = roxml_get_chld_nb(MaterialNode);
    *Materials = (TrivIGAMaterialStruct *)
	          IritMalloc(sizeof(TrivIGAMaterialStruct) * (*NumMaterials));

    i = -1;
    while ((InternalMaterialNode = roxml_get_chld(MaterialNode, NULL, ++i))
								    != FALSE) {
	InternalMaterialAttribNode = roxml_get_attr(InternalMaterialNode,
	                                            "id", 0);
	if (InternalMaterialAttribNode == NULL) {
	    (*Materials)[i].Id = -1;
	    continue;
	}

	roxml_get_content(InternalMaterialAttribNode,
	                  TempStrBuffer,
	                  TRIV_IGA_MAX_MATERIAL_NAME_LEN,
	                  &TempStrBufferSize);

	(*Materials)[i].Id = atoi(TempStrBuffer);

	InternalMaterialAttribNode = roxml_get_attr(InternalMaterialNode,
	                                            "name", 0);

	if (InternalMaterialAttribNode != NULL) {
	    roxml_get_content(InternalMaterialAttribNode,
		              (*Materials)[i].Name,
		              TRIV_IGA_MAX_MATERIAL_NAME_LEN,
		              &TempStrBufferSize);
	}

	InternalMaterialAttribNode = roxml_get_attr(InternalMaterialNode,
	    "type", 0);
	if (InternalMaterialAttribNode != NULL) {
	    roxml_get_content(InternalMaterialAttribNode,
			      (*Materials)[i].Type,
			      TRIV_IGA_MAX_MATERIAL_NAME_LEN,
			      &TempStrBufferSize);
	}

	(*Materials)[i].NumProperties = roxml_get_chld_nb(InternalMaterialNode);
	j = -1;

	while ((InternalMaterialChildNode = 
	         roxml_get_chld(InternalMaterialNode, NULL, ++j)) != FALSE) {
	      roxml_get_name(InternalMaterialChildNode,
			     (*Materials)[i].Properties[j].Name,
			     TRIV_IGA_MAX_MATERIAL_NAME_LEN);
	      roxml_get_content(InternalMaterialChildNode,
				(*Materials)[i].Properties[j].Value,
				TRIV_IGA_MAX_MATERIAL_NAME_LEN,
				&TempStrBufferSize);
	}
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   comparison functor for fields, compares the ids of the fields.           *
*                                                                            *
* PARAMETERS:                                                                *
*   VPt1:   Pointer to the first field.                  		     *
*   VPt2:   Pointer to the second field.                                     *
*									     *
* RETURN VALUE:                                                              *
*   int:   0 if the same id, 1 if the first id is bigger, and -1 if the      *
*          second id is bigger.			                             *
*****************************************************************************/
static int TrivIGAFieldCompFunc(const VoidPtr VPt1, const VoidPtr VPt2)
{
    IrtRType
	Diff = (*(TrivIGAFieldStruct **) VPt1) -> ID - 
	       (*(TrivIGAFieldStruct **) VPt2) -> ID;

    return IRIT_SIGN(Diff);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Adds the materials defined in the given arrangement to the XML root.     *
*                                                                            *
* PARAMETERS:                                                                *
*   H:      A handle on the IGA arrangement.   		       		     *
*   Root:   XML root.		                                             *
*									     *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void TrivIGAAddMaterialsToXML(struct TrivIGAArrangementStruct *H,
				     node_t *Root)
{
    node_t *Materials;
    node_t *Material;
    int i, j, NumChildren;
    char MatID[TRIV_IGA_MAX_MATERIAL_NAME_LEN];
    const TrivIGAMaterialStruct *CurrentMaterial;

    /* Search for material node and update it. */
    Materials = roxml_get_chld(Root, "Material", 0);
    if (Materials != NULL) {
	NumChildren = roxml_get_chld_nb(Materials);
	for (i = NumChildren - 1; i >= 0; --i) {
	    roxml_del_node(roxml_get_chld(Materials, NULL, i));
	}
    } 
    else {
	Materials  = roxml_add_node(Root, 0, ROXML_ELM_NODE, 
	                            "Material", NULL);
    }

    for (i = 1; i < H -> NumMaterials; ++i) {
	CurrentMaterial = (H -> Materials[i]);
	sprintf(MatID, "%d", CurrentMaterial -> Id);
	Material  = roxml_add_node(Materials, 0, ROXML_ELM_NODE, 
	                           "material", NULL);
	roxml_add_node(Material, 0, ROXML_ATTR_NODE, "id"  , MatID);
	roxml_add_node(Material, 0, ROXML_ATTR_NODE, "type", 
	               (char *)CurrentMaterial -> Type);
	roxml_add_node(Material, 0, ROXML_ATTR_NODE, "name", 
	               (char *)CurrentMaterial -> Name);

	for (j = 0; j < CurrentMaterial -> NumProperties; ++j) {
	    roxml_add_node(Material,
		           0,
		           ROXML_ELM_NODE,
		           (char *)CurrentMaterial -> Properties[j].Name,
		           (char *)CurrentMaterial -> Properties[j].Value );
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Adds the geometry defined in the given arrangement to the given XML root.*
*                                                                            *
* PARAMETERS:                                                                *
*   H:      A handle on the IGA arrangement.             		     *
*   Root:   XML root.                             	                     *
*									     *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void TrivIGAAddGeometryToXML(struct TrivIGAArrangementStruct *H,
				    node_t *Root)
{
    IRIT_STATIC_DATA const CagdRType 
	ScalarCtrlPoint[4] = { 1.0f, 0.0f, 0.0f, 0.0f };
    IRIT_STATIC_DATA char TmpBuff[TRIV_IGA_MAX_XML_ROW_LEN],
	                  TmpBuffList[TRIV_IGA_MAX_XML_ROW_LEN];
    IRIT_STATIC_DATA node_t *Geometry, *Nodes, *Weights, *Elements, *Volumes,
	*Multivariate, *OneNode;
    int i, j, Len, FieldIndex, NumChildren, NumFields,
	WKnotsNum = 0,
	VKnotsNum = 0,
	UKnotsNum = 0,
	GenTVID = 1;
    CagdBType
	IsTheSameWeight = TRUE;
    const CagdRType *CtrlPnt;
    TrivIGATVStruct *FieldTV;
    struct TrivIGACtlPtUniqueIDsStruct *CtrlPntUinqIDStruct;
    CagdRType
	LastWeight = IRIT_INFNTY;
    TrivIGAFieldStruct *Field,
	**SortedFieldsArray = NULL;

    /* Remove Geometry node. */
    Geometry = roxml_get_chld(Root, "Geometry", 0);
    if (Geometry != NULL) {
	NumChildren = roxml_get_chld_nb(Geometry);
	for (i = NumChildren - 1; i >= 0; --i) {
	    roxml_del_node(roxml_get_chld(Geometry, NULL, i));
	}
    } 
    else {
	Geometry  = roxml_add_node(Root, 0, ROXML_ELM_NODE, "Geometry", NULL);
    }

    Nodes  = roxml_add_node(Geometry, 0, ROXML_ELM_NODE, "Nodes", NULL);

    /* Add Nodes. */
    for (i = 0; i <= H -> UniqueGlblCtlPtIDMax; ++i) {
	if (TrivIGAIsPointOfFieldType(H, i, TRIV_IGA_SCALAR_FIELD_TYPE)) {
	    CtrlPnt = ScalarCtrlPoint;
	} 
	else {
	    CtrlPnt = TrivIGAGetUniquePoint(H, i);
	}

	if (CtrlPnt) {
	    sprintf(TmpBuff, "%f,%f,%f", 
		    CtrlPnt[1] / CtrlPnt[0], 
		    CtrlPnt[2] / CtrlPnt[0], 
		    CtrlPnt[3] / CtrlPnt[0]);
	    OneNode = roxml_add_node(Nodes, 0, ROXML_ELM_NODE,
		                     "node", TmpBuff);
	    sprintf(TmpBuff, "%d", i + 1);
	    roxml_add_node(OneNode, 0, ROXML_ATTR_NODE, "id", TmpBuff);
	    if (IsTheSameWeight) {
		if (LastWeight == IRIT_INFNTY) {
		    LastWeight = CtrlPnt[0];
		} 
		else if (!IRIT_APX_EQ(LastWeight, CtrlPnt[0])) {
		    IsTheSameWeight = FALSE;
		}
	    }
	}
    }

    /* Add Weights. */
    if (!IsTheSameWeight) {
	Weights = roxml_add_node(Geometry, 0, ROXML_ELM_NODE, "Weights", NULL);
	for (i = 0; i <= H -> UniqueGlblCtlPtIDMax; ++i) {
	    CtrlPnt = TrivIGAGetUniquePoint(H, i);
	    if (CtrlPnt) {
		sprintf(TmpBuff, "%f", CtrlPnt[0]);
		OneNode = roxml_add_node(Weights, 0, ROXML_ELM_NODE, 
		                         "node", TmpBuff);
		sprintf(TmpBuff, "%d", i + 1);
		roxml_add_node(OneNode, 0, ROXML_ATTR_NODE, "id", TmpBuff);
	    }
	}
    }

    /* Add Elements. */
    Elements  = roxml_add_node(Geometry, 0, ROXML_ELM_NODE, "Elements", NULL);

    NumFields = CagdListLength(H -> Fields);
    if (NumFields == 0) {
	return;	
    }

    SortedFieldsArray = (TrivIGAFieldStruct **)
	IritMalloc(NumFields * sizeof(TrivIGAFieldStruct *));

    for (Field = H -> Fields, i = 0; 
	Field != NULL; 
	Field = Field -> Pnext, ++i) {
	    SortedFieldsArray[i] = Field;
    }

    /* sort the fields with increasing id order */
    qsort(SortedFieldsArray, NumFields, sizeof(TrivIGAFieldStruct *),
	  TrivIGAFieldCompFunc);

    for (FieldIndex = 0; FieldIndex < NumFields; ++FieldIndex) {
	Field = SortedFieldsArray[FieldIndex];
	sprintf(TmpBuff, "%d", Field -> ID);
	Volumes  = roxml_add_node(Elements, 0, ROXML_ELM_NODE, 
	                          "iga_volume", NULL);
	roxml_add_node(Volumes, 0, ROXML_ATTR_NODE, "id", TmpBuff);

	sprintf(TmpBuff, "%d", CagdListLength(Field -> TVs));
	roxml_add_node(Volumes, 0, ROXML_ATTR_NODE, "num_unknowns", TmpBuff);
	roxml_add_node(Volumes, 0, ROXML_ATTR_NODE, "num_multivars", TmpBuff);
	roxml_add_node(Volumes, 0, ROXML_ATTR_NODE, "type_formulation", 
	               Field -> NamedType);

	sprintf(TmpBuff, "%d", Field -> Material -> Id);
	roxml_add_node(Volumes, 0, ROXML_ATTR_NODE, "material", TmpBuff);

	for (j = 0; j < Field ->NumProperties; ++j) {
	    roxml_add_node(Volumes, 0, ROXML_ATTR_NODE, 
		           Field -> Properties[j].Name,
		           Field -> Properties[j].Value);
	}

	GenTVID = 1;		  /* The trivariates ids is local per field. */
	for (FieldTV = Field -> TVs; 
	     FieldTV != NULL; 
	     FieldTV = FieldTV -> Pnext) {
	    TrivTVStruct 
	        *CurrentTV = FieldTV -> TV;

	    Multivariate =  roxml_add_node(Volumes, 0, ROXML_ELM_NODE, 
					   "multivar", NULL);

	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, "type", 
			   TRIV_IS_RATIONAL_TV(CurrentTV) ? "NURBS" : 
			                                  "BSpline");

	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "name", Field -> NamedType);

	    sprintf(TmpBuff, "%d", GenTVID++);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "id", TmpBuff);

	    WKnotsNum = CurrentTV -> WLength + CurrentTV -> WOrder;
	    if ( CurrentTV -> WPeriodic )
	        --WKnotsNum;
	    sprintf(TmpBuff, "%d", WKnotsNum);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "n_knots_w", TmpBuff);

	    VKnotsNum = CurrentTV -> VLength + CurrentTV -> VOrder;
	    if (CurrentTV -> VPeriodic )
	        --VKnotsNum;
	    sprintf(TmpBuff, "%d", VKnotsNum);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "n_knots_v", TmpBuff);

	    UKnotsNum = CurrentTV -> ULength + CurrentTV -> UOrder;
	    if ( CurrentTV -> UPeriodic )
	        --UKnotsNum;
	    sprintf(TmpBuff, "%d", UKnotsNum);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "n_knots_u", TmpBuff);

	    sprintf(TmpBuff,"%d",CurrentTV -> WOrder - 1);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "degree_w", TmpBuff);

	    sprintf(TmpBuff, "%d", CurrentTV -> VOrder - 1);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "degree_v", TmpBuff);
	    
	    sprintf(TmpBuff, "%d", CurrentTV -> UOrder - 1);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "degree_u", TmpBuff);

	    sprintf(TmpBuff, "%d", CurrentTV -> WLength);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "n_nodes_w", TmpBuff);

	    sprintf(TmpBuff, "%d", CurrentTV -> VLength);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "n_nodes_v", TmpBuff);

	    sprintf(TmpBuff, "%d", CurrentTV -> ULength);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE, 
			   "n_nodes_u", TmpBuff);

	    sprintf(TmpBuff, "%d", 3);
	    roxml_add_node(Multivariate, 0, ROXML_ATTR_NODE,
			   "dim", TmpBuff);

	    /* U knots. */
	    TmpBuffList[0] = '\0';
	    for (i = 0; i < UKnotsNum - 1; ++i) {
	        sprintf(TmpBuff, "%f,", CurrentTV -> UKnotVector[i]);
		strcat(TmpBuffList, TmpBuff);
	    }	    
	    sprintf(TmpBuff, "%f", CurrentTV -> UKnotVector[UKnotsNum - 1]);
	    strcat(TmpBuffList, TmpBuff);
	    roxml_add_node(Multivariate, 0, ROXML_ELM_NODE, 
			   "knots_u", TmpBuffList);

	    /* V knots. */
	    TmpBuffList[0] = '\0';
	    for (i = 0; i < VKnotsNum - 1; ++i) {
	       sprintf(TmpBuff, "%f,", CurrentTV -> VKnotVector[i]);
	       strcat(TmpBuffList, TmpBuff);
	    }
	    sprintf(TmpBuff, "%f", CurrentTV -> VKnotVector[VKnotsNum - 1]);
	    strcat(TmpBuffList, TmpBuff);
	    roxml_add_node(Multivariate, 0, ROXML_ELM_NODE, 
			   "knots_v", TmpBuffList);

	    /* W knots. */
	    TmpBuffList[0] = '\0';
	    for (i = 0; i < WKnotsNum - 1; ++i) {
	        sprintf(TmpBuff, "%f,", CurrentTV -> WKnotVector[i]);
		strcat(TmpBuffList, TmpBuff);
	    }
	    sprintf(TmpBuff, "%f", CurrentTV -> WKnotVector[WKnotsNum - 1]);
	    strcat(TmpBuffList, TmpBuff);
	    roxml_add_node(Multivariate, 0, ROXML_ELM_NODE, "knots_w",
			   TmpBuffList);

	    /* Connectivity. */
	    TmpBuffList[0] = '\0';
	    Len = TRIV_TV_UPT_LST_LEN(CurrentTV) *
	          TRIV_TV_VPT_LST_LEN(CurrentTV) *
	          TRIV_TV_WPT_LST_LEN(CurrentTV);
	    CtrlPntUinqIDStruct = FieldTV -> CtlPtsIDs;
	    for (i = 0; i < Len - 1; ++i) {
	        sprintf(TmpBuff, "%d,", CtrlPntUinqIDStruct[i].ID + 1);
		strcat(TmpBuffList, TmpBuff);
	    }
	    sprintf(TmpBuff, "%d", CtrlPntUinqIDStruct[Len - 1].ID + 1);
	    strcat(TmpBuffList, TmpBuff);
	    roxml_add_node(Multivariate, 0, ROXML_ELM_NODE, 
			   "connectivity", TmpBuffList);
	}
    }

    IritFree(SortedFieldsArray);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Adds the boundaries defined in the given arrangement to the              *
* given XML root.							     *
*                                                                            *
* PARAMETERS:                                                                *
*   H:      A handle on the IGA arrangement.            		     *
*   Root:   XML root.                                                        *
*									     *
* RETURN VALUE:                                                              *
*   void								     *
*****************************************************************************/
static void TrivIGAAddBoundaryToXML(struct TrivIGAArrangementStruct *H,
				    node_t *Root)
{
    node_t 
	*XMLBoundaryNode = NULL, 
	*XMLFixedNode = NULL, 
	*XMLPrescribedNode = NULL,
	*XMLSelectedBoundaryNode = NULL,
	*CurrentNode = NULL;
    TrivIGABoundaryNodeStruct *BoundaryNodeTmp;
    char TmpBuff[4256], *NodeValueStr;

    XMLBoundaryNode  = 
	roxml_add_node(Root, 0, ROXML_ELM_NODE, "Boundary", NULL);

    TrivIGAGenNeighboringConstraints(TrivIGADataManagerGetArrangementID(H), 
				     XMLBoundaryNode, 
				     TrivIGANeighboringConstraintXML);

    BoundaryNodeTmp = H -> BoundaryNodes;
    while (BoundaryNodeTmp) {
	switch (BoundaryNodeTmp -> BoundaryType) {
	    case TRIV_IGA_BOUNDARY_FIXED:
	    {
	        if (XMLFixedNode == NULL) {
		    XMLFixedNode = roxml_add_node(XMLBoundaryNode, 0, 
						  ROXML_ELM_NODE, "fix", 
						  NULL);
		}
		XMLSelectedBoundaryNode = XMLFixedNode;
		NodeValueStr = NULL;
		break;
	    }
	    case TRIV_IGA_BOUNDARY_PRESCRIBED:
	    {
	        if (XMLPrescribedNode == NULL) {
		  XMLPrescribedNode = roxml_add_node(XMLBoundaryNode, 0, 
						     ROXML_ELM_NODE,
						     "prescribe", NULL);
		}
		XMLSelectedBoundaryNode = XMLPrescribedNode;
		sprintf(TmpBuff, "%f", BoundaryNodeTmp -> Value);
		NodeValueStr = TmpBuff;
		break;
	    }
	    default:
		assert(0);
		NodeValueStr = NULL;
		break;
	}

	CurrentNode = roxml_add_node(XMLSelectedBoundaryNode, 0, 
	                             ROXML_ELM_NODE, "node", NodeValueStr);
	sprintf(TmpBuff, "%d", BoundaryNodeTmp -> NodeID + 1);
	roxml_add_node(CurrentNode, 0, ROXML_ATTR_NODE, "id", TmpBuff);	
	roxml_add_node(CurrentNode, 0, ROXML_ATTR_NODE, "bc", 
	               BoundaryNodeTmp -> BoundaryAxisConditions);

	BoundaryNodeTmp = BoundaryNodeTmp -> Pnext;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Saves the input IGA arrangement as Febio XML file, loads template        M
* file and updates the material, geometry and boundary sections.             M
*                                                                            *
* PARAMETERS:                                                                M
*   ArgmntID:         A handle on the IGA arrangement to process.            M
*   FileName:         Output XML file name.				     M
*   TemplateFileName: Template XML file.                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:  TRUE if successful, FALSE otherwise.		                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivIGAExportToXML                                                       M
*****************************************************************************/
int TrivIGAExportToXML(TrivIGAArrangementID ArgmntID,
		       const char *FileName, 
		       const char *TemplateFileName)
{
    node_t *InternalNode, 
	*Root = roxml_load_doc((char *) TemplateFileName);
    struct TrivIGAArrangementStruct 
	*H = TrivIGADataManagerGetArrangement(ArgmntID);

    if (Root == NULL) {
	InternalNode = 
	    roxml_add_node(NULL, 0, ROXML_PI_NODE, "xml", 
	                   "version=\"1.0\" encoding=\"ISO-8859-1\"");

	Root = roxml_add_node(InternalNode, 0, ROXML_ELM_NODE, 
	                      "febio_spec", NULL);

	roxml_add_node(Root, 0, ROXML_ATTR_NODE, "version", "1.2");
	InternalNode  = roxml_add_node(Root, 0, ROXML_ELM_NODE, 
	                               "Module", NULL);

	roxml_add_node(InternalNode, 0, ROXML_ATTR_NODE, "type", "solid");
	InternalNode = roxml_add_node(Root, 0, ROXML_ELM_NODE, 
	                              "Control", NULL);
	roxml_add_node(InternalNode, 0, ROXML_ELM_NODE, "time_steps", "20");
	roxml_add_node(InternalNode, 0, ROXML_ELM_NODE, "step_size", "0.1");
    }

    TrivIGAAddMaterialsToXML(H, Root);
    TrivIGAAddGeometryToXML(H, Root);
    TrivIGAAddBoundaryToXML(H, Root);

    roxml_commit_changes(Root, (char *)FileName, NULL, 1);
    roxml_close(Root);

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Returns the coordinates of the given point id.                           *
*                                                                            *
* PARAMETERS:                                                                *
*   H:    A handle on the IGA arrangement.		       		     *
*   Id:   Node ID.                    		                             *
*									     *
* RETURN VALUE:                                                              *
*   CagdRType *: Pointer to node coordinates.                                *
*****************************************************************************/
static CagdRType* TrivIGAGetUniquePoint(struct TrivIGAArrangementStruct *H, 
					int Id) 
{
    IRIT_STATIC_DATA CagdRType ResPnt[CAGD_MAX_PT_SIZE];
    CagdBType IsNotRational;
    int i, j,  Len, MaxCoord;
    TrivIGAFieldStruct *Field;
    TrivIGATVStruct *FieldTV;
    struct TrivIGACtlPtUniqueIDsStruct *CtrlPntUinqIDStruct;

    ResPnt[0] = 1.0f;
    ResPnt[1] = 0.0f;
    ResPnt[2] = 0.0f;
    ResPnt[3] = 0.0f;
    for (Field = H -> Fields; Field != NULL; Field = Field -> Pnext) {
	for (FieldTV = Field -> TVs;
	     FieldTV != NULL;
	     FieldTV = FieldTV -> Pnext) {
	    IsNotRational = !TRIV_IS_RATIONAL_TV(FieldTV -> TV);
	    MaxCoord = CAGD_NUM_OF_PT_COORD(FieldTV -> TV -> PType);
	    Len = TRIV_TV_UPT_LST_LEN(FieldTV -> TV) *
		  TRIV_TV_VPT_LST_LEN(FieldTV -> TV) *
		  TRIV_TV_WPT_LST_LEN(FieldTV -> TV);

	    CtrlPntUinqIDStruct = FieldTV -> CtlPtsIDs;
	    for (i = 0; i < Len; ++i) {
		if (CtrlPntUinqIDStruct[i].ID == Id) {
		    for (j = IsNotRational; j <= MaxCoord; ++j) {
			ResPnt[j] = FieldTV -> TV -> Points[j][i];
		    }
		    return ResPnt;
		}
	    }		 
	}
    }
    return NULL;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Returns if the given point is participating in a field of the given type. *
*                                                                            *
* PARAMETERS:                                                                *
*   H:           A handle on the IGA arrangement.          		     *
*   Id:          Node ID.                                                    *
*   FieldType:   Field type.                                                 *
*									     *
* RETURN VALUE:                                                              *
*   CagdBType:  True if the point participates in the field of the given     *
*               type.                                                        *
*****************************************************************************/
static CagdBType TrivIGAIsPointOfFieldType(struct TrivIGAArrangementStruct *H, 
					   int Id,
					   TrivIGAFieldType FieldType)
{
    int i, Len;
    TrivIGAFieldStruct *Field;
    TrivIGATVStruct *FieldTV;
    CagdBType IsPointInField;
    struct TrivIGACtlPtUniqueIDsStruct *CtrlPntUinqIDStruct;
 
    for (Field = H -> Fields; Field != NULL; Field = Field -> Pnext) {
	IsPointInField = FALSE;
	for (FieldTV = Field -> TVs; 
	     FieldTV != NULL && !IsPointInField; 
	     FieldTV = FieldTV -> Pnext) {
	    Len = TRIV_TV_UPT_LST_LEN(FieldTV -> TV) *
		  TRIV_TV_VPT_LST_LEN(FieldTV -> TV) *
		  TRIV_TV_WPT_LST_LEN(FieldTV -> TV);

	    CtrlPntUinqIDStruct = FieldTV -> CtlPtsIDs;
	    for (i = 0; i < Len; ++i) {
		if ( CtrlPntUinqIDStruct[i].ID == Id ) {
		    if ( Field -> ValueType != FieldType ) {
			return FALSE;
		    }
		    IsPointInField = TRUE;
		    break;
		}
	    }		 
	}
    }

    return TRUE;
}

#else

/* Create stab functions that only issue a warning... */
int TrivIGALoadMaterialXML(const char *FileName,
			   TrivIGAMaterialStruct **Materials,
			   int *NumMaterials)
{
    fprintf(stderr, "Unsupported function TrivIGALoadMaterialXML called and ignored\n");
    return FALSE;
}

int TrivIGAExportToXML(TrivIGAArrangementID ArgmntID,
		       const char *FileName, 
		       const char *TemplateFileName)
{
    fprintf(stderr, "Unsupported function TrivIGAExportToXML called and ignored\n");
    return FALSE;
}

#endif /* IRIT_HAVE_XML_LIB */
