/*****************************************************************************
*   "Irit" - the 3d (not only polygonal) solid modeller.		     *
*									     *
* Written by:  Gershon Elber				Ver 0.2, Mar. 1990   *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
*   Module to evaluate the binary tree generated by the InptPrsr module.     *
*   All the objects are handled the same but the numerical one, which is     *
* moved as a IrtRType and not as an object (only internally within this	     *
* module) as it is frequently used and consumes much less memory this way.   *
*****************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "program.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "freeform.h"
#include "inptprsg.h"
#include "inptprsl.h"
#include "objects.h"
#include "overload.h"

IRIT_STATIC_DATA InptEvalUpdateAssignedObjFuncType
    GlblEvalUpdateAssignedObjCBFunc = NULL;

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets a call back function that is invoked on any object that just        M
* underwent an assignment operation.                                         M
*                                                                            *
* PARAMETERS:                                                                M
*   UpdateAssignObjFunc: New call back function or NULL to disable.          M
*                                                                            *
* RETURN VALUE:                                                              M
*   InptEvalUpdateAssignedObjFuncType:  Old call back function value.        M
*                                                                            *
* KEYWORDS:                                                                  M
*   InptPrsrSetUpdateAssignedObjectFunc                                      M
*****************************************************************************/
InptEvalUpdateAssignedObjFuncType InptPrsrSetUpdateAssignedObjectFunc(
		      InptEvalUpdateAssignedObjFuncType UpdateAssignObjFunc)
{
    InptEvalUpdateAssignedObjFuncType
	OldVal = GlblEvalUpdateAssignedObjCBFunc;

    GlblEvalUpdateAssignedObjCBFunc = UpdateAssignObjFunc;

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Do type checking to the given parsed tree - return type if found one       M
* or returns ERROR_EXPR if error in types was detected.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   Root:      Of parse tree to verify.                                      M
*   Level:     Of recursion. Used to identify top level (Level == 0).        M
*                                                                            *
* RETURN VALUE:                                                              M
*   IritExprType:  ERROR_EXPR if found error, or expression type otherwise.  M
*                                                                            *
* KEYWORDS:                                                                  M
*   InptPrsrTypeCheck                                                        M
*****************************************************************************/
IritExprType InptPrsrTypeCheck(ParseTree *Root, int Level)
{
    IritExprType Right, Left, Result;

    if (IP_IS_NUM_FUNCTION(Root -> NodeKind)) {/* Funcs returning Real Type: */
	if (IritEvalFuncParamMismatch(Root))
	    return ERROR_EXPR;
	return NUMERIC_EXPR;
    }

    if (IP_IS_GEN_FUNCTION(Root -> NodeKind)) {  /* Funcs returning nothing: */
	if (Level == 0) {
	    if (IritEvalFuncParamMismatch(Root))
		return ERROR_EXPR;
	    return NO_EXPR;
	}
	else {
	    IPGlblEvalError = IE_ERR_TYPE_MISMATCH;
	    UpdateCharError("Procedure ", Root -> NodeKind, Root);
	    return ERROR_EXPR;
	}
    }

    if (IP_IS_OBJ_FUNCTION(Root -> NodeKind)) {  /* Funcs returning objects: */
	if (IritEvalFuncParamMismatch(Root))
	    return ERROR_EXPR;
	return ObjFuncTable[Root -> NodeKind - IP_OBJ_FUNC_OFFSET].RetType;
    }

    switch (Root -> NodeKind) {
	case IP_TKN_PLUS:
	case IP_TKN_MINUS:
	case IP_TKN_MULT:
	case IP_TKN_DIV:
	case IP_TKN_POWER:
	    Right = InptPrsrTypeCheck(Root -> Right, Level + 1);
	    Left  = InptPrsrTypeCheck(Root -> Left,  Level + 1);
	    if (Right == ERROR_EXPR || Left == ERROR_EXPR)
		return ERROR_EXPR;
	    if (!OverLoadTypeCheck(Root -> NodeKind, Right, Left, &Result)) {
		IPGlblEvalError = IE_ERR_TYPE_MISMATCH;
                UpdateCharError("Operator ", Root -> NodeKind, Root);
		return ERROR_EXPR;
	    }
	    else
		return Result;
	case IP_TKN_UNARMINUS:
	    if ((Right = InptPrsrTypeCheck(Root -> Right, Level + 1))
							== ERROR_EXPR)
		return ERROR_EXPR;
	    else if (!OverLoadTypeCheck(Root -> NodeKind, Right, NO_EXPR,
								&Result)) {
		IPGlblEvalError = IE_ERR_TYPE_MISMATCH;
                UpdateCharError("Operator ", Root -> NodeKind, Root);
		return ERROR_EXPR;
	    }
	    else
		return Result;
	case IP_TKN_EQUAL:
	    if ((Right = InptPrsrTypeCheck(Root -> Right, Level + 1))
							== ERROR_EXPR)
		return ERROR_EXPR;
	    if (Root -> Left -> NodeKind != IP_TKN_PARAMETER) {
		IPGlblEvalError = IE_ERR_ASSIGN_LEFT_OP;
		InptPrsrPrintTree(Root -> Left, IPGlblCharData,
				  INPUT_LINE_LEN);
		return ERROR_EXPR;
	    }
	    else if (Root -> Left -> PObj -> ObjType == IP_OBJ_UNDEF)
		SET_TO_BE_ASSIGN_OBJ(Root -> Left -> PObj);
	    return Right;
	case IP_TKN_BOOL_AND:
	case IP_TKN_BOOL_OR:
	    if (InptPrsrTypeCheck(Root -> Right, Level + 1) != NUMERIC_EXPR ||
		InptPrsrTypeCheck(Root -> Left,  Level + 1) != NUMERIC_EXPR)
		return ERROR_EXPR;
	    return NUMERIC_EXPR;
	case IP_TKN_BOOL_NOT:
	    if (InptPrsrTypeCheck(Root -> Right, Level + 1) != NUMERIC_EXPR)
		return ERROR_EXPR;
	    return NUMERIC_EXPR;
	case IP_TKN_CMP_EQUAL:
	case IP_TKN_CMP_NOTEQUAL:
	case IP_TKN_CMP_LSEQUAL:
	case IP_TKN_CMP_GTEQUAL:
	case IP_TKN_CMP_LESS:
	case IP_TKN_CMP_GREAT:
	    Right = InptPrsrTypeCheck(Root -> Right, Level + 1);
	    Left  = InptPrsrTypeCheck(Root -> Left,  Level + 1);
	    if (Right == ERROR_EXPR || Left == ERROR_EXPR)
		return ERROR_EXPR;
	    return NUMERIC_EXPR;
	case IP_TKN_NUMBER:
	    return NUMERIC_EXPR;
	case IP_TKN_PARAMETER:
	    return InptPrsrObjType2Expr(Root, Root -> PObj);
	case IP_TKN_STRING:
	    return STRING_EXPR;
	case IP_USERFUNCDEF:
	    return ANY_EXPR;
	case IP_USERPROCDEF:
	    return NO_EXPR;
	case IP_USERINSTDEF:
	    return ANY_EXPR;
	case IP_TKN_COLON:
	    if (Root -> Left &&
		Root -> Left -> NodeKind == IP_TKN_EQUAL &&
		Root -> Left -> Right &&
		(Root -> Left -> Right -> NodeKind == IP_USERPROCDEF ||
		 Root -> Left -> Right -> NodeKind == IP_USERFUNCDEF)) {
		ParseTree *Body;

		/* A special form of function/procedure definition. */
		for (Body = Root;
		     Body -> NodeKind == IP_TKN_COLON;
		     Body = Body -> Right);
		if (Root -> Left -> Left -> NodeKind != IP_TKN_PARAMETER ||
		    Body == NULL ||
		    Body == Root)
		    return ERROR_EXPR;		    
	    }
	    else {
		Left  = InptPrsrTypeCheck(Root -> Left,  0);
		Right = InptPrsrTypeCheck(Root -> Right, 0);
		if (Right == ERROR_EXPR || Left == ERROR_EXPR)
		    return ERROR_EXPR;
	    }
	    return NO_EXPR;
	default:
	    /* Could only happen if InptPrsrTypeCheck is invoked from the    */
	    /* error printing function UpdateCharError. Do not recurse.      */
	    if (IPGlblEvalError == IPE_NO_ERR) {
		IPGlblEvalError = IE_ERR_FATAL_ERROR;
		UpdateCharError("Token ", Root -> NodeKind, Root);
	    }
	    return ERROR_EXPR;
    }
}

/* Disable the function with no prototype warning on Borland's compilers. */
#ifdef __BORLANDC__
#pragma warn -pro
#endif /* __BORLANDC__ */
/*****************************************************************************
* DESCRIPTION:                                                               M
*   Converts an object type to an expression type.                           M
*                                                                            *
* PARAMETERS:                                                                M
*   Root:  Of parse tree to convert.                                         M
*   PObj:  To convert its type to expression type.                           M
*                                                                            *
* RETURN VALUE:                                                              M
*   IritExprType:   The expression type.                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   InptPrsrObjType2Expr                                                     M
*****************************************************************************/
IritExprType InptPrsrObjType2Expr(ParseTree *Root, IPObjectStruct *PObj)
{
    switch (PObj -> ObjType) {
	case IP_OBJ_POLY:
	    return POLY_EXPR;
	case IP_OBJ_NUMERIC:
	    return NUMERIC_EXPR;
	case IP_OBJ_POINT:
	    return POINT_EXPR;
	case IP_OBJ_VECTOR:
	    return VECTOR_EXPR;
	case IP_OBJ_PLANE:
	    return PLANE_EXPR;
	case IP_OBJ_CTLPT:
	    return CTLPT_EXPR;
	case IP_OBJ_MATRIX:
	    return MATRIX_EXPR;
	case IP_OBJ_STRING:
	    return STRING_EXPR;
	case IP_OBJ_LIST_OBJ:
	    return OLST_EXPR;
	case IP_OBJ_CURVE:
	    return CURVE_EXPR;
	case IP_OBJ_SURFACE:
	    return SURFACE_EXPR;
	case IP_OBJ_TRIMSRF:
	    return TRIMSRF_EXPR;
	case IP_OBJ_MODEL:
	    return MODEL_EXPR;
	case IP_OBJ_TRIVAR:
	    return TRIVAR_EXPR;
	case IP_OBJ_MULTIVAR:
	    return MULTIVAR_EXPR;
	case IP_OBJ_INSTANCE:
	    if ((PObj = IritDBGetObjByName(PObj -> U.Instance -> Name))
		                                                   == NULL) {
		IPGlblEvalError = IE_ERR_UNDEF_INSTANCE;
		UpdateCharError("Token ", Root -> NodeKind, Root);
		return ERROR_EXPR;
	    }
	    return InptPrsrObjType2Expr(Root, PObj);
	case IP_OBJ_TRISRF:
	    return TRISRF_EXPR;
	default:
	    if (IS_TO_BE_ASSIGN_OBJ(PObj)) {
		/* A block of expressions separated by colons has  */
		/* an assignment followed by usage herein. No idea */
		/* of the type forces us to allow anything and     */
		/* hope for the best.				   */
		return ANY_EXPR;
	    }
	    else {
		IPGlblEvalError = IE_ERR_IP_OBJ_UNDEFINED;
		sprintf(IPGlblCharData,	IRIT_EXP_STR("Object = %s, Type %d"),
			IP_GET_OBJ_NAME(PObj), PObj -> ObjType);
		return ERROR_EXPR;
	    }
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Evaluates the given parsed and verified tree (via InptPrsrTypeCheck rtn).  M
*    The tree is modified, in place, during the evaluation process.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Root:      Parsed tree to evaluate.                                      M
*   Level:     Of recursion. Used to identify top level (Level == 0).        M
*                                                                            *
* RETURN VALUE:                                                              M
*   ParseTree *:   Result of evaluation.                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   InptPrsrEvalTree                                                         M
*****************************************************************************/
ParseTree *InptPrsrEvalTree(ParseTree *Root, int Level)
{
    int Index, NumOfParam,
	PrintIt = ((Level & IP_EVAL_BASE_MASK) == 0 &&
					    Root -> NodeKind != IP_TKN_EQUAL);
    const char *ErrorMsg;
    char *p, Name[IRIT_LINE_LEN_VLONG];
    ParseTree *TempL, *TempR, *Params[IP_MAX_PARAM],
	*RetVal = NULL;
    VoidPtr ParamPtrs[IP_MAX_PARAM];

    if (IP_IS_NUM_FUNCTION(Root -> NodeKind)) {/* Funcs which rets Real Type:*/
	Index = Root -> NodeKind - IP_NUM_FUNC_OFFSET;
	NumOfParam = NumFuncTable[Index].NumOfParam;

	switch (Root -> NodeKind) {
	    case IP_ASIN:  /* Real return functions with one real parameter. */
	    case IP_ACOS:
	    case IP_TAN:
	    case IP_COS:
	    case IP_EXP:
	    case IP_ABS:
	    case IP_LN:
	    case IP_LOG:
	    case IP_SIN:
	    case IP_SQRT:
	    case IP_ATAN:
	    case IP_FLOOR:
	    case IP_TIME:
		if (InptEvalFetchParameters(Root, NULL, 1, Level,
					    Params, ParamPtrs) != 1)
		    break;
		/* Use table entries to call the function directly. */
		Root -> PObj =
		    IPGenNUMValObject((NumFuncTable[Index].Func)
				        (Params[0] -> PObj -> U.R));
		RetVal = Root;
		break;

	    case IP_ATAN2:
	    case IP_FMOD:
	    case IP_POWER:
	    case IP_RANDOM:
		if (InptEvalFetchParameters(Root, NULL, 2, Level,
					    Params, ParamPtrs) != 2)
		    break;
		/* Use table entries to call the function directly. */
		Root -> PObj =
		    IPGenNUMValObject((NumFuncTable[Index].Func)
				        (Params[0] -> PObj -> U.R,
					 Params[1] -> PObj -> U.R));
		RetVal = Root;
		break;

	    default:
		if (InptEvalFetchParameters(Root,
			(FuncTableType *) &NumFuncTable[Index],
			NumOfParam, Level, Params, ParamPtrs) != NumOfParam)
		    break;

		/* Use table entries to call the function directly. */
		switch (NumFuncTable[Index].NumOfParam) {
		    case 0:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)());
			break;
		    case 1:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0]));
			break;
		    case 2:
			Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1]));
			break;
		    case 3:
			Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2]));
			break;
		    case 4:
			Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3]));
		        break;
		    case 5:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4]));
		        break;
		    case 6:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5]));
		        break;
		    case 7:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5],
						 ParamPtrs[6]));
		        break;
		    case 8:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5],
						 ParamPtrs[6], ParamPtrs[7]));
		        break;
		    case 9:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5],
						 ParamPtrs[6], ParamPtrs[7],
						 ParamPtrs[8]));
		        break;
		    case 10:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5],
						 ParamPtrs[6], ParamPtrs[7],
						 ParamPtrs[8], ParamPtrs[9]));
		        break;
		    case 11:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5],
						 ParamPtrs[6], ParamPtrs[7],
						 ParamPtrs[8], ParamPtrs[9],
						 ParamPtrs[10]));
		        break;
		    case 12:
		        Root -> PObj =
			    IPGenNUMValObject((NumFuncTable[Index].Func)
					        (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4], ParamPtrs[5],
						 ParamPtrs[6], ParamPtrs[7],
						 ParamPtrs[8], ParamPtrs[9],
						 ParamPtrs[10], ParamPtrs[11]));
		        break;
		    default:
			assert(0);
			break;
		}
		RetVal = Root;
		break;
	}
	if (RetVal && RetVal -> PObj)
	    RetVal -> PObj -> Count++;
    }
    else if (IP_IS_OBJ_FUNCTION(Root -> NodeKind)) {/* Funcs returning objs: */
	Index = Root -> NodeKind - IP_OBJ_FUNC_OFFSET;
	NumOfParam = ObjFuncTable[Index].NumOfParam;

	switch (Root -> NodeKind) {
	    case IP_COORD:
	        if (InptEvalFetchParameters(Root, NULL, 2, Level,
					    Params, ParamPtrs) != 2)
		    break;
	        /* Use table entries to call the function directly. */
		Root -> PObj =
		    (ObjFuncTable[Index].Func)
		        (Params[0] -> PObj, &Params[1] -> PObj -> U.R);
		if (Root -> PObj == NULL)
		    break;

		RetVal = Root;
		break;

	    case IP_NTH:
		if (InptEvalFetchParameters(Root, NULL, 2,
					    Level, Params, ParamPtrs) != 2)
		    break;
		/* Use table entries to call the function directly. */
		Root -> PObj =
		    (ObjFuncTable[Index].Func)
		        (Params[0] -> PObj, &Params[1] -> PObj -> U.R);
		if (Root -> PObj == NULL)
		    break;

		RetVal = Root;
		break;

	    case IP_LOAD:
		if (InptEvalFetchParameters(Root, NULL, 1,
					    Level, Params, ParamPtrs) != 1)
		    break;
		/* Use table entries to call the function directly. */
		Root -> PObj =
		    (ObjFuncTable[Index].Func)(Params[0] -> PObj -> U.Str, "");
		if (Root -> PObj == NULL) {
	            LoadSaveObjectParseError(&ErrorMsg);
		    IPGlblEvalError = IE_ERR_DATA_PRSR_ERROR;
		    strcpy(IPGlblCharData, ErrorMsg);
		    break;
		}

		RetVal = Root;
		break;

	    case IP_CONTOUR:
		if (InptEvalFetchParameters(Root, (FuncTableType *)
					    &ObjFuncTable[Index], 5,
					    Level, Params, ParamPtrs) == 5) {
		    Root -> PObj = (ObjFuncTable[Index].Func)
		                                (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 ParamPtrs[4]);
		}
		else if (InptEvalFetchParameters(Root, (FuncTableType *)
					    &ObjFuncTable[Index], 4,
					    Level, Params, ParamPtrs) == 4) {
		    Root -> PObj = (ObjFuncTable[Index].Func)
		                                (ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], ParamPtrs[3],
						 NULL);
		}
		else if (InptEvalFetchParameters(Root, (FuncTableType *)
						 &ObjFuncTable[Index], 3, Level,
						 Params, ParamPtrs) == 3) {
		    Root -> PObj = (ObjFuncTable[Index].Func)
						(ParamPtrs[0], ParamPtrs[1],
						 ParamPtrs[2], NULL, NULL);
		}
		else
		    Root -> PObj = NULL;

                if (Root -> PObj == NULL)
                    break;

                RetVal = Root;
		break;

	    case IP_CTLPT:
	    case IP_LIST:
		/* Use table entries to call the function directly. */
		Root -> PObj = (ObjFuncTable[Index].Func)(Root -> Right);
		if (Root -> PObj == NULL)
		    break;

		RetVal = Root;
		break;

	    default:
		if (InptEvalFetchParameters(Root,
			(FuncTableType *) &ObjFuncTable[Index],
			NumOfParam, Level, Params, ParamPtrs) != NumOfParam)
		    break;

		/* Use table entries to call the function directly. */
		switch (ObjFuncTable[Index].NumOfParam) {
		    case 0:
			Root -> PObj = (ObjFuncTable[Index].Func)();
			break;
		    case 1:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0]);
			break;
		    case 2:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1]);
			break;
		    case 3:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2]);
			break;
		    case 4:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3]);
			break;
		    case 5:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4]);
			break;
		    case 6:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5]);
			break;
		    case 7:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6]);
			break;
		    case 8:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7]);
			break;
		    case 9:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7], ParamPtrs[8]);
			break;
		    case 10:
			Root -> PObj = (ObjFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7], ParamPtrs[8],
				 ParamPtrs[9]);
			break;
		    case 11:
			Root -> PObj = (ObjFuncTable[Index].Func)
			    (ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
			    ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
			    ParamPtrs[6], ParamPtrs[7], ParamPtrs[8],
			    ParamPtrs[9], ParamPtrs[10]);
			break;
		    case 12:
			Root -> PObj = (ObjFuncTable[Index].Func)
			    (ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
			    ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
			    ParamPtrs[6], ParamPtrs[7], ParamPtrs[8],
			    ParamPtrs[9], ParamPtrs[10], ParamPtrs[11]);
			break;
		    default:
			assert(0);
			break;
		}
		if (Root -> PObj == NULL)
		    break;

		RetVal = Root;
		break;
	}
	if (RetVal && RetVal -> PObj)
	    RetVal -> PObj -> Count++;
    }
    else if (IP_IS_GEN_FUNCTION(Root -> NodeKind)) {  /* Funcs rets nothing: */
	Index = Root -> NodeKind - IP_GEN_FUNC_OFFSET;
	NumOfParam = GenFuncTable[Index].NumOfParam;

	switch (Root -> NodeKind) {
	    case IP_SAVE:
		if (InptEvalFetchParameters(Root, NULL, 2, Level,
					    Params, ParamPtrs) != 2)
		    break;

		/* Use table entries to call the function directly. */
		(GenFuncTable[Index].Func)(Params[0] -> PObj -> U.Str,
					   Params[1] -> PObj);

		/* Save the matrix. */
		strncpy(Name, Params[0] -> PObj -> U.Str, IRIT_LINE_LEN - 5);
		if ((p = strstr(Name, ".itd")) != NULL ||
		    (p = strstr(Name, ".ITD")) != NULL ||
		    (p = strstr(Name, ".icd")) != NULL ||
		    (p = strstr(Name, ".ICD")) != NULL) {
		    *p = 0;
		    strcat(Name, ".imd");
		    IritDispViewSaveMatrix(Name);
		}

	        if (LoadSaveObjectParseError(&ErrorMsg) != 0) {
		    IPGlblEvalError = IE_ERR_DATA_PRSR_ERROR;
		    strcpy(IPGlblCharData, ErrorMsg);
		    break;
	        }
		break;

	    case IP_FREE:
		if (InptEvalFetchParameters(Root, NULL, 1, Level, Params,
					    ParamPtrs) != 1 ||
	            (TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								   == NULL)
		    break;

	        if (!IP_VALID_OBJ_NAME(TempR -> PObj)) {
		    IPGlblEvalError = IE_ERR_FREE_SIMPLE;
                    UpdateCharError("Procedure ", IP_FREE, NULL);
		    RetVal = Root;
		    break;
	        }
	        /* Use table entries to call the function directly. */
	        (GenFuncTable[Index].Func)(TempR -> PObj);
	        TempR -> PObj = NULL;	    /* Make sure its disconnected... */
		break;

	    case IP_IF:
		switch (NumOfParam =
			          InptEvalCountNumParameters(Root -> Right)) {
		    case 2:
			InptEvalIfCondition(
			    InptEvalFetchParameter(Root -> Right, 0, 2),
			    InptEvalFetchParameter(Root -> Right, 1, 2),
			    NULL);
			break;
		    case 3:
			InptEvalIfCondition(
			    InptEvalFetchParameter(Root -> Right, 0, 3),
			    InptEvalFetchParameter(Root -> Right, 1, 3),
			    InptEvalFetchParameter(Root -> Right, 2, 3));
			break;
		    default:
			IPGlblEvalError = IE_ERR_NUM_PRM_MISMATCH;
			sprintf(IPGlblCharData,
				IRIT_EXP_STR("IF clause (2 or 3 expected, found %d)"),
				NumOfParam);
			break;
		}
		break;

	    case IP_FOR:
		InptEvalForLoop(InptEvalFetchParameter(Root -> Right, 0, 4),
				InptEvalFetchParameter(Root -> Right, 1, 4),
				InptEvalFetchParameter(Root -> Right, 2, 4),
				InptEvalFetchParameter(Root -> Right, 3, 4));
		break;

	    case IP_WHILE:
                InptEvalWhileLoop(InptEvalFetchParameter(Root -> Right, 0, 2),
                                  InptEvalFetchParameter(Root -> Right, 1, 2));
                break;

	    default:
		if (InptEvalFetchParameters(Root,
			(FuncTableType *) &GenFuncTable[Index],
			NumOfParam, Level, Params, ParamPtrs) != NumOfParam)
		    break;

		/* Use table entries to call the function directly. */
		switch (GenFuncTable[Index].NumOfParam) {
		    case 0:
			(GenFuncTable[Index].Func)();
			break;
		    case 1:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0]);
			break;
		    case 2:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1]);
			break;
		    case 3:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2]);
			break;
		    case 4:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3]);
			break;
		    case 5:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4]);
			break;
		    case 6:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5]);
			break;
		    case 7:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6]);
			break;
		    case 8:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7]);
			break;
		    case 9:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7], ParamPtrs[8]);
			break;
		    case 10:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7], ParamPtrs[8],
				 ParamPtrs[9]);
			break;
		    case 11:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7], ParamPtrs[8],
				 ParamPtrs[9], ParamPtrs[10]);
			break;
		    case 12:
			(GenFuncTable[Index].Func)
				(ParamPtrs[0], ParamPtrs[1], ParamPtrs[2],
				 ParamPtrs[3], ParamPtrs[4], ParamPtrs[5],
				 ParamPtrs[6], ParamPtrs[7], ParamPtrs[8],
				 ParamPtrs[9], ParamPtrs[10], ParamPtrs[11]);
			break;
		    default:
			assert(0);
			break;
		}
		break;
	}
	RetVal = Root;
	if (Root -> PObj) {
	    Root -> PObj -> ObjType = IP_OBJ_UNDEF;
	    Root -> PObj -> Count = 1;
	}
	else {
	    Root -> PObj = IPAllocObject("", IP_OBJ_UNDEF, NULL);
	    Root -> PObj -> Count++;
	}
    }
    else {
        switch (Root -> NodeKind) {		  /* The rest of the world. */
	    case IP_TKN_PLUS:
	    case IP_TKN_MINUS:
	    case IP_TKN_MULT:
	    case IP_TKN_DIV:
	    case IP_TKN_POWER:
	        if (((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								== NULL) ||
		    ((TempL = InptPrsrEvalTree(Root -> Left,  Level + 1))
								== NULL))
		    break;
		TempR = OverLoadEvalOper(Root, TempR, TempL,
					 &IPGlblEvalError, IPGlblCharData);
		RetVal = TempR;
		break;

	    case IP_TKN_UNARMINUS:
		if ((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								      == NULL)
		    break;
		TempR = OverLoadEvalOper(Root, TempR, NULL,
					 &IPGlblEvalError, IPGlblCharData);
		RetVal = TempR;
		break;

	    case IP_TKN_COLON:
		if (Root -> Left &&
		    Root -> Left -> NodeKind == IP_TKN_EQUAL &&
		    Root -> Left -> Right &&
		    (Root -> Left -> Right -> NodeKind == IP_USERPROCDEF ||
		     Root -> Left -> Right -> NodeKind == IP_USERFUNCDEF)) {
		    /* A special form of function/procedure definition. */
		    InptEvalDefineFunc(Root);
		}
		else {
		    InptPrsrEvalTree(Root -> Left, IP_EVAL_NEXT_LEVEL(Level));
		    InptPrsrEvalTree(Root -> Right, IP_EVAL_NEXT_LEVEL(Level));
		}
		break;

	    case IP_TKN_NUMBER:
		RetVal = Root;
		break;

	    case IP_TKN_PARAMETER:
		RetVal = Root;
		break;

	    case IP_TKN_STRING:
		RetVal = Root;
		break;

	    case IP_TKN_EQUAL:
		if ((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								== NULL)
		    break;
		TempL = Root -> Left;

		if (TempL -> PObj == TempR -> PObj) {
		    RetVal = TempR; /* A = A. */
		    break;
		}
		else if (TempL -> PObj == NULL) {
		    TempL -> PObj = IPAllocObject("", IP_OBJ_UNDEF, NULL);
		}
		else if (IPListObjectFind(TempR -> PObj, TempL -> PObj)) {
		    /* Object is in the list. To prevent a loop in the list */
		    /* structure, create a new object to hold the result.   */
		    strcpy(Name, IP_GET_OBJ_NAME(TempL -> PObj));
		    IritDBDeleteObject(TempL -> PObj, FALSE);
		    TempL -> PObj -> Count--;
		    IPFreeObject(TempL -> PObj);
		    
		    TempL -> PObj = IPAllocObject(Name, IP_OBJ_UNDEF, NULL);
		    TempL -> PObj -> Count++;
		    IritDBInsertObject(TempL -> PObj, FALSE);
		}
		if (GlblHandleDependencies) {
		    IPODObjectDpndncyStruct
			*Dpnds = TempL -> PObj -> Dpnds;

		    /* Saved dependencies before copy to recover original. */
		    TempL -> PObj -> Dpnds = NULL;
		    IPCopyObject(TempL -> PObj, TempR -> PObj, FALSE);
		    IPODFreeDependencies(TempL -> PObj -> Dpnds);
		    TempL -> PObj -> Dpnds = Dpnds;

		    InptEvalPropagateDependencies(TempL -> PObj, Root);
		}
		else
		    IPCopyObject(TempL -> PObj, TempR -> PObj, FALSE);

		/* Should we propagate names through the entire hierarchy? */
		if (GlblPropagateNames) {
		    IPPropagateObjectName(TempL -> PObj, NULL);
		    if (GlblMakeAllHierarchyVisible)
		        IritDBInsertHierarchy(TempL -> PObj, FALSE);
		}

		/* Have a call back function to update obj on assignment? */
		if (GlblEvalUpdateAssignedObjCBFunc != NULL)
		    GlblEvalUpdateAssignedObjCBFunc(TempL -> PObj);

		RetVal = TempL;
		break;

	    case IP_TKN_BOOL_AND:
	        if (((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								== NULL) ||
		    ((TempL = InptPrsrEvalTree(Root -> Left,  Level + 1))
								== NULL))
		    break;
		if (Root -> PObj)
		    Root -> PObj -> ObjType = IP_OBJ_NUMERIC;
		else {
		    Root -> PObj = IPAllocObject("", IP_OBJ_NUMERIC, NULL);
		    Root -> PObj -> Count++;
		}
		Root -> PObj -> U.R = (!IRIT_APX_EQ(TempR -> PObj -> U.R, 0.0) &&
				       !IRIT_APX_EQ(TempL -> PObj -> U.R, 0.0));
		RetVal = Root;
		break;

	    case IP_TKN_BOOL_OR:
	        if (((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								== NULL) ||
		    ((TempL = InptPrsrEvalTree(Root -> Left,  Level + 1))
								== NULL))
		    break;
		if (Root -> PObj)
		    Root -> PObj -> ObjType = IP_OBJ_NUMERIC;
		else {
		    Root -> PObj = IPAllocObject("", IP_OBJ_NUMERIC, NULL);
		    Root -> PObj -> Count++;
		}
		Root -> PObj -> U.R = (!IRIT_APX_EQ(TempR -> PObj -> U.R, 0.0) ||
				       !IRIT_APX_EQ(TempL -> PObj -> U.R, 0.0));
		RetVal = Root;
		break;

	    case IP_TKN_BOOL_NOT:
		if ((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								      == NULL)
		    break;
		if (Root -> PObj)
		    Root -> PObj -> ObjType = IP_OBJ_NUMERIC;
		else {
		    Root -> PObj = IPAllocObject("", IP_OBJ_NUMERIC, NULL);
		    Root -> PObj -> Count++;
		}
		Root -> PObj -> U.R = IRIT_APX_EQ(TempR -> PObj -> U.R, 0.0);
		RetVal = Root;
		break;

	    case IP_TKN_CMP_EQUAL:
	    case IP_TKN_CMP_NOTEQUAL:
	    case IP_TKN_CMP_LSEQUAL:
	    case IP_TKN_CMP_GTEQUAL:
	    case IP_TKN_CMP_LESS:
	    case IP_TKN_CMP_GREAT:
	        if (((TempR = InptPrsrEvalTree(Root -> Right, Level + 1))
								== NULL) ||
		    ((TempL = InptPrsrEvalTree(Root -> Left,  Level + 1))
								== NULL))
		    break;
		RetVal = InptEvalCompareObject(Root, TempL, TempR,
					   &IPGlblEvalError, IPGlblCharData);
		break;

	    case IP_USERINSTDEF:
		if (InptEvalCountNumParameters(Root -> Right) !=
		    Root -> UserFunc -> NumParams) {
		    IPGlblEvalError = IE_ERR_NUM_PRM_MISMATCH;
		    sprintf(IPGlblCharData,
			    IRIT_EXP_STR("%d expected in function \"%s\""),
			    Root -> UserFunc -> NumParams,
			    Root -> UserFunc -> FuncName);
		    break;
		}
		if (InptEvalFetchParameters(Root, NULL,
			Root -> UserFunc -> NumParams,
			Level, Params, ParamPtrs) !=
		    Root -> UserFunc -> NumParams)
		    break;

		RetVal = InptEvalUserFunc(Root, Params);
		break;
	}
    }

    if (PrintIt &&
	RetVal &&
	RetVal -> PObj &&
	RetVal -> PObj -> ObjType != IP_OBJ_UNDEF) {
	RetVal -> PObj -> Count--;	       /* Remove ref from this tree. */
	PrintIritObject(RetVal -> PObj);
	RetVal -> PObj -> Count++;		      /* Add reference back. */
    }

    return RetVal;
}

/* Restore the function with no prototype warning on Borland's compilers. */
#ifdef __BORLANDC__
#pragma warn .pro
#endif /* __BORLANDC__ */
