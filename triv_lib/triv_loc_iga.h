/*****************************************************************************
* Module to support IGA (isogeometric analysis) - Definitions file           *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber & Fady Massarwi		Ver 1.0, June 2012   *
*****************************************************************************/

#ifndef TRIV_LOC_IGA_H
#define TRIV_LOC_IGA_H

#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "triv_loc.h"
#include "inc_irit/triv_iga.h"
#include "inc_irit/mvar_lib.h"

#define TRIV_IGA_NEIGHBORHOOD_EPS	1e-4
#define TRIV_IGA_MAX_XML_ROW_LEN	164256

#define TRIV_IGA_MAX_ARRANGEMENTS 	1024
#define TRIV_IGA_MAX_TRIVARIATES 	(TRIV_IGA_MAX_ARRANGEMENTS * 128)

/* An IGA trivariate also consists of neighbors information, IDs, etc. */
typedef struct TrivIGACtlPtUniqueIDsStruct {
    int ID;
} TrivIGACtlPtUniqueIDsStruct;

typedef struct TrivIGATVStruct {
    struct TrivIGATVStruct *Pnext;
    TrivTVStruct *TV;
    TrivTVStruct *DuTV, *DvTV, *DwTV,              /* Derivative TVs to TV. */
	*Du2TV, *Dv2TV, *Dw2TV, *DuDvTV, *DuDwTV, *DvDwTV;
    TrivIGAAdjacencyInfoStruct Neighbors[6];/* The up to 6 neighbors of TV. */
    int UniqueCtlPtIDMin, UniqueCtlPtIDMax;
    struct TrivIGACtlPtUniqueIDsStruct *CtlPtsIDs;
    struct TrivIGAFieldStruct *Field;/* Back reference to the parent field. */
} TrivIGATVStruct;

/* An IGA field consists of several IGA trivariates. */
typedef struct TrivIGAFieldStruct {
    struct TrivIGAFieldStruct *Pnext;
    TrivIGATVStruct *TVs;                         /* A list of trivariates. */
    TrivIGAFieldType ValueType;
    TrivIGAMaterialStruct *Material;
    char NamedType[TRIV_IGA_MAX_FIELD_TYPE_LEN];
    TrivIGAXMLProperty Properties[TRIV_IGA_MAX_XML_PROPERTY_LEN];
    int NumProperties;
    int ID;
} TrivIGAFieldStruct;

typedef struct TrivIGABoundaryNodeStruct {
    struct TrivIGABoundaryNodeStruct *Pnext;
    int NodeID;
    TrivIGANodeBoundaryType BoundaryType;
    char *BoundaryAxisConditions;
    CagdRType Value;
} TrivIGABoundaryNodeStruct;

/* For the 3 (U, v, W) direction, Alpha sets the (parametric) size-ratio    */
/* between the last interval and the first.  NumOfIntervals specifies the   */
/* number of intervals in each direction, with zero value to deactivate.    */
typedef struct TrivIGASeedingStateStruct {
    CagdRType AlphaVal[3];
    int NumOfIntervals[3];
    CagdRType DmnMin[3], DmnMax[3];
} TrivIGASeedingStateStruct;

typedef struct TrivIGAArrangementStruct {
    TrivIGAFieldStruct *Fields;  /* List of IGA fields (a few trivariates). */
    TrivIGAErrorType LastError;
    int UniqueGlblCtlPtIDMax;
    TrivIGAMaterialStruct *Materials[TRIV_IGA_MAX_MATERIALS];
    int NumMaterials;
    TrivIGABoundaryNodeStruct *BoundaryNodes;
    TrivIGASeedingStateStruct SeedingState;
} TrivIGAArrangementStruct;

typedef struct TrivIGAErrorStruct {
    TrivFatalErrorType ErrorNum;
    const char *ErrorDesc;
} TrivIGAErrorStruct;

TrivIGAArrangementID TrivIGADataManagerAllocateArrangement();
TrivIGAArrangementStruct *TrivIGADataManagerGetArrangement( 
					     TrivIGAArrangementID ArrngmntID);
TrivIGAArrangementID TrivIGADataManagerGetArrangementID(
					     TrivIGAArrangementStruct *H);
int TrivIGADataManagerFreeArrangement(TrivIGAArrangementID ArrngmntID);
TrivIGATVID TrivIGADataManagerAddTrivariate(TrivIGAArrangementID ArgmntID, 
					    TrivTVStruct *TV,
					    int ID);

#endif /* TRIV_LOC_IGA_H */
