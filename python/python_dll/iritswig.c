/*****************************************************************************
*   Main module of "Irit" - the 3d (not only polygonal) solid modeller.      *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Usage:								     *
*   Irit [-t] [-g] [-i] [-z] [file[.irt]]				     *
*									     *
* Written by:  Gershon Elber				Ver 3.0, Apr. 1990   *
*****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "program.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/iritprsr.h"
#include "ctrl-brk.h"
#include "dosintr.h"
#include "inptprsg.h"
#include "inptprsl.h"
#include "objects.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/bool_lib.h"
#include "inc_irit/trim_lib.h"
#include "inc_irit/triv_lib.h"
#include "inc_irit/symb_lib.h"
#include "inc_irit/user_lib.h"

#ifdef USE_VARARGS
#include <varargs.h>
#else
#include <stdarg.h>
#endif /* USE_VARARGS */
#ifdef __WINNT__
#include <direct.h>
#endif /* __WINNT__ */

IRIT_GLOBAL_DATA IPObjectStruct
    *GlblObjList = NULL;		   /* All objects defined on system. */

IRIT_GLOBAL_DATA jmp_buf GlblLongJumpBuffer;	  /* Used in error recovery. */

IRIT_GLOBAL_DATA int
    GlblInterpProd	 = TRUE;

#define TRUE  1
#define FALSE 0
typedef int BOOL;

#define IRIT_TO_SWIG_LINE_LENGTH 500

static void Irit2SwigNumeric(FILE* h_file, FILE* c_file, FILE* i_file);
static void Irit2SwigIPObject(FILE* h_file, FILE* c_file, FILE* i_file);
static void Irit2SwigProcedures();
static void PrintNumericFuncDeclaration(FILE* file, int nIndex);
static void PrintNumericFuncDefinition(FILE* file, int nIndex);
static void PrintIPObjectFuncDeclaration(FILE* file, int nIndex, BOOL bToIFile);
static void PrintIPObjectFuncDefinition(FILE* file, int nIndex);
static void PrintProceduresFuncDeclaration(FILE* file, int nIndex, BOOL bToIFile);
static void PrintProceduresFuncDefinition(FILE* file, int nIndex);

static void Irit2SwigAddons(FILE *i_file, FILE *input_file);
static int SkipNumeric(char* FuncName);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Main module of IRIT - Read command line and do what is needed...	     M
*                                                                            *
* PARAMETERS:                                                                M
*   argc, argv:  Command line.                                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   main                                                                     M
*****************************************************************************/
int main(int argc, char **argv)
{
    FILE* h_file = fopen("iritpy_interface.h", "w");
    FILE* c_file = fopen("iritpy_interface.c", "w");
    FILE* i_file = fopen("irit.i", "w");
    FILE* inp_file1 = fopen("irit_to_swig1.in", "r");
    FILE* inp_file2 = fopen("python_link.h", "r");
    FILE* inp_file3 = fopen("irit_to_swig2.in", "r");

    Irit2SwigAddons(i_file, inp_file1);

    fprintf(c_file, "#include \"iritpy_interface.h\"\n\n");
    fprintf(c_file, "#include \"inc_irit/attribut.h\"\n");
    fprintf(c_file, "#include \"dosintr.h\"\n");
    fprintf(c_file, "#include \"program.h\"\n");
    fprintf(c_file, "#include \"inptprsg.h\"\n\n");
    fprintf(c_file, "#include \"inptprsl.h\"\n\n");

    fprintf(h_file, "#ifndef _IRITPY_INTERFACE_H_\n");
    fprintf(h_file, "#define _IRITPY_INTERFACE_H_\n\n");
    fprintf(h_file, "#include \"inc_irit/grap_lib.h\"\n");
    fprintf(h_file, "#include \"inc_irit/ext_lib.h\"\n");
    fprintf(h_file, "#include \"objects.h\"\n");
    fprintf(h_file, "#include \"freeform.h\"\n");
    fprintf(h_file, "#include \"bsc_geom.h\"\n");

    fprintf(h_file, "\n\n/* Numeric returning function: */\n");
    fprintf(i_file, "\n\n/* Numeric returning function: */\n");
    Irit2SwigNumeric(h_file, c_file, i_file);

    fprintf(h_file, "\n\n/* IPObject return functions: */\n");
    fprintf(i_file, "\n\n/* IPObject return functions: */\n");
    Irit2SwigIPObject(h_file, c_file, i_file);

    fprintf(h_file, "\n\n/* Procedures: */\n");
    fprintf(i_file, "\n\n/* Procedures: */\n");
    Irit2SwigProcedures(h_file, c_file, i_file);

    Irit2SwigAddons(i_file, inp_file2);
    Irit2SwigAddons(i_file, inp_file3);

    fprintf(h_file, "\n#endif\n");

    fclose(h_file);
    fclose(c_file);
    fclose(i_file);
    fclose(inp_file1);
    fclose(inp_file2);

    return 0;
}

void ExecOneLine(char *Line)
{
}

void IritExit0(void)
{
}

void IritIdleFunction(int MiliSeconds)
{
}

void IritPrintOneString(const char *Line)
{
    fprintf(stderr, "%s\n", Line);
}

void DefaultFPEHandler(int Type)
{
}


void IritExit(int ExitCode)
{
}

IPObjectStruct *SetIritState(char *Name, IPObjectStruct *Data)
{
    return NULL;
}

IPObjectStruct *IritPickCrsrClientEvent(IrtRType *RWaitTime)
{
    return NULL;
}

#ifndef __WINNT__
char *_strlwr(const char *Str)
{
    static char StrLwr[IRIT_LINE_LEN_LONG];

    strncpy(StrLwr, Str, IRIT_LINE_LEN_LONG);

    return IritStrLower(StrLwr);
}
#endif /* !__WINNT__ */

char *GetParamFromType(IritExprType Param)
{
    switch (Param) {
	case 0:
	    return "void";
	case NUMERIC_EXPR:
	    return "IrtRType";
	case POINT_EXPR | VECTOR_EXPR:
	    return "IrtPtType";
	case POINT_EXPR:
	    return "IrtPtType";
	case VECTOR_EXPR:
	    return "IrtVecType";
	case PLANE_EXPR:
	    return "IrtPlnType";
	case STRING_EXPR:
	    return "char *";
	case CTLPT_EXPR: /* Fall through. */
	case MATRIX_EXPR:
	case CURVE_EXPR: 
	case SURFACE_EXPR:         
	case OLST_EXPR: 
	case TRIMSRF_EXPR: 
	case TRIVAR_EXPR: 
	case INSTANCE_EXPR: 
	case TRISRF_EXPR: 
	case MULTIVAR_EXPR: 
	case POLY_CURVE_EXPR: 
	case OLST_GEOM_EXPR: 
	case ANY_EXPR:
	    return "IPObjectStruct*"; 
	default:
	    return "IPObjectStruct*"; /* all the |(or) cases */
    }
}


void Irit2SwigNumeric(FILE* file_h, FILE* file_c, FILE* file_i)
{
    int i;

    for (i = 0; i < NumFuncTableSize; i++) {
	if (SkipNumeric(NumFuncTable[i].FuncName)) {
	    continue;
	}

	/* Special case:
	Can't use sizeof, there id a python function by that name. */
	if (strcmp("SIZEOF", NumFuncTable[i].FuncName) == 0)  {
	    fprintf(file_i, "%%rename(%s) Py%s;\n", "SizeOf",
		    _strlwr(NumFuncTable[i].FuncName));
	}
	else {
	    fprintf(file_i, "%%rename(%s) Py%s;\n",
		    _strlwr(NumFuncTable[i].FuncName),
		    _strlwr(NumFuncTable[i].FuncName));
	}
		
	PrintNumericFuncDeclaration(file_i, i);
	fprintf(file_i, ";\n\n");

	PrintNumericFuncDeclaration(file_h, i);
	fprintf(file_h, ";\n");

	PrintNumericFuncDeclaration(file_c, i);
	fprintf(file_c, "\n");

	PrintNumericFuncDefinition(file_c, i);
	fprintf(file_c, "\n");

    }

}

void PrintNumericFuncDeclaration(FILE* file, int nIndex)
{
    int j;

    /* Special case:
    Discrepancy between parameter number 2 (should be IrtRType by the rules). */
    if (stricmp("meshsize", NumFuncTable[nIndex].FuncName) == 0) {
	fprintf(file, "IrtRType Pymeshsize(IPObjectStruct* Param1, IrtRType* Param2)");
	return;
    }
    /* Special case:
    Discrepancy between parameters number 2 and 3 (should be IrtRType by the rules). */
    else if (stricmp("hausdrpts", NumFuncTable[nIndex].FuncName) == 0) {
	fprintf(file, "IrtRType Pyhausdrpts(IPObjectStruct* Param1, IPObjectStruct* Param2, IrtRType* Param3, IrtRType* Param4)");
	return;
    }
    /* Special case:
    Discrepancy between parameters number 2 and 3 (should be IrtRType by the rules). */
    else if (stricmp("zcollide", NumFuncTable[nIndex].FuncName) == 0) {
	fprintf(file, "IrtRType Pyzcollide(IPObjectStruct* Param1, IPObjectStruct* Param2, IrtRType* Param3, IrtRType* Param4)");
	return;
    }
    /* Special case:
    Discrepancy between parameters number 2 (should be IrtRType by the rules). */
    else if (stricmp("tvolume", NumFuncTable[nIndex].FuncName) == 0) {
	fprintf(file, "IrtRType Pytvolume(IPObjectStruct* Param1, IrtRType* Param2)");
	return;
    }
    
    fprintf(file, "IrtRType ");
    fprintf(file, "Py");
    fprintf(file, _strlwr(NumFuncTable[nIndex].FuncName));
    fprintf(file, "(");
    for (j = 0; j < NumFuncTable[nIndex].NumOfParam - 1; j++) {
	fprintf(file, GetParamFromType((NumFuncTable[nIndex].ParamObjType)[j]));
	fprintf(file, " Param%d", j + 1);
	fprintf(file, ", ");
    }
    if ((NumFuncTable[nIndex].ParamObjType)[j] != 0) {
	fprintf(file, GetParamFromType((NumFuncTable[nIndex].ParamObjType)[j]));
	fprintf(file, " Param%d", j + 1);
    }

    fprintf(file, ")");    
}

void PrintNumericFuncDefinition(FILE* file, int nIndex)
{
    int j;

    fprintf(file, "{\n");
    fprintf(file, "\treturn ");
    fprintf(file, NumFuncTable[nIndex].CFuncName);
    fprintf(file, "(");
    for (j = 0; j < NumFuncTable[nIndex].NumOfParam - 1; j++) {
	fprintf(file, "Param%d", j + 1);
	fprintf(file, ", ");
    }
    fprintf(file, "Param%d", j + 1);
    fprintf(file, ");\n}\n");    
}



void Irit2SwigIPObject(FILE* file_h, FILE* file_c, FILE* file_i)
{
    int i;

    for (i = 0; i < ObjFuncTableSize; i++) {
	/* Skip those functions.
	We implement them in python. */
	if ((stricmp("ctlpt", ObjFuncTable[i].FuncName) == 0) ||
	    (stricmp("list", ObjFuncTable[i].FuncName) == 0)) {
	    continue;
	}
	if (stricmp("Contour", ObjFuncTable[i].FuncName) == 0) {
	    fprintf(file_i, "%%rename(Wrap%s) PyWrap%s;\n",
		    _strlwr(ObjFuncTable[i].FuncName),
		    _strlwr(ObjFuncTable[i].FuncName));
	}
	else {
	    fprintf(file_i, "%%rename(%s) Py%s;\n",
		    _strlwr(ObjFuncTable[i].FuncName),
		    _strlwr(ObjFuncTable[i].FuncName));
	}
	PrintIPObjectFuncDeclaration(file_i, i, TRUE);
	fprintf(file_i, ";\n\n");

	PrintIPObjectFuncDeclaration(file_h, i, FALSE);
	fprintf(file_h, ";\n");

	PrintIPObjectFuncDeclaration(file_c, i, FALSE);
	fprintf(file_c, "\n");

	PrintIPObjectFuncDefinition(file_c, i);
	fprintf(file_c, "\n");
    }
}

void PrintIPObjectFuncDeclaration(FILE* file, int nIndex, BOOL bToIFile)
{
    int j;

    /* Special case:
       contour is function which can get 3-5 parameters, so it should have 
       a wrapper Wrapcontour. */
    if (stricmp("contour", ObjFuncTable[nIndex].FuncName) == 0) {
	fprintf(file, "IPObjectStruct* PyWrapcontour(IPObjectStruct* Param1, IPObjectStruct* Param2, IPObjectStruct* Param3, IPObjectStruct* Param4, IPObjectStruct* Param5)");
	return;
    }

    fprintf(file, "IPObjectStruct* ");
    fprintf(file, "Py");
    fprintf(file, _strlwr(ObjFuncTable[nIndex].FuncName));
    fprintf(file, "(");
    for (j = 0; j < ObjFuncTable[nIndex].NumOfParam - 1; j++) {
	fprintf(file, GetParamFromType((ObjFuncTable[nIndex].ParamObjType)[j]));
	if ((ObjFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && bToIFile)
	    fprintf(file, "* INPUT");
	else if ((ObjFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && !bToIFile)
	    fprintf(file, "* Param%d", j + 1);
	else
	    fprintf(file, " Param%d", j + 1);
	fprintf(file, ", ");
    }
    if ((ObjFuncTable[nIndex].ParamObjType)[j] != 0) {
	fprintf(file, GetParamFromType((ObjFuncTable[nIndex].ParamObjType)[j]));
	if ((ObjFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && bToIFile)
	    fprintf(file, "* INPUT");
	else if ((ObjFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && !bToIFile)
	    fprintf(file, "* Param%d", j + 1);
	else
	    fprintf(file, " Param%d", j + 1);    	
    }
    fprintf(file, ")");    
}

void PrintIPObjectFuncDefinition(FILE* file, int nIndex)
{
    int j = 0;
    int NumOfLoops;

    /* Special case: contour is function which can get 3-5 parameters. */
    if (stricmp("contour", ObjFuncTable[nIndex].FuncName) == 0) 
	        fprintf(file, "{\n\treturn ContourFreeform(Param1, Param2 -> U.Plane, &Param3 -> U.R, Param4, Param5);\n}\n");
    else {
	  NumOfLoops = ObjFuncTable[nIndex].NumOfParam;

	  fprintf(file, "{\n");
	  fprintf(file, "\treturn ");
	  fprintf(file, ObjFuncTable[nIndex].CFuncName);
	  fprintf(file, "(");
	  for (j = 0; j < NumOfLoops - 1; j++) {
	      fprintf(file, "Param%d", j + 1);
	      fprintf(file, ", ");
	  }
	  if (NumOfLoops != 0) {
	      fprintf(file, "Param%d", j + 1);
	  }

	  fprintf(file, ");\n}\n");
      }
}

void Irit2SwigProcedures(FILE* file_h, FILE* file_c, FILE* file_i)
{
    int i;

    for (i = 0; i < GenFuncTableSize; i++) {
	/* Skip translating those functions, because we use python version. */
	if ((stricmp("IF", GenFuncTable[i].FuncName) == 0) ||
	    (stricmp("FOR", GenFuncTable[i].FuncName) == 0) ||
	    (stricmp("WHILE", GenFuncTable[i].FuncName) == 0) ||
	    (stricmp("EXEC", GenFuncTable[i].FuncName) == 0) ||
	    (stricmp("IDYNMEM", GenFuncTable[i].FuncName) == 0) ||
	    (stricmp("IQUERY", GenFuncTable[i].FuncName) == 0)) {
	    continue;
	}


	fprintf(file_i, "%%rename(%s) Py%s;\n",
		_strlwr(GenFuncTable[i].FuncName),
		_strlwr(GenFuncTable[i].FuncName));
	PrintProceduresFuncDeclaration(file_i, i, TRUE);
	fprintf(file_i, ";\n\n");

	PrintProceduresFuncDeclaration(file_h, i, FALSE);
	fprintf(file_h, ";\n");

	PrintProceduresFuncDeclaration(file_c, i, FALSE);
	fprintf(file_c, "\n");

	PrintProceduresFuncDefinition(file_c, i);
	fprintf(file_c, "\n");
    }      
}

void PrintProceduresFuncDeclaration(FILE* file, int nIndex, BOOL bToIFile)
{
    int j;

    fprintf(file, "void ");
    fprintf(file, "Py");		
    fprintf(file, _strlwr(GenFuncTable[nIndex].FuncName));
    fprintf(file, "(");
    for (j = 0; j < GenFuncTable[nIndex].NumOfParam - 1; j++) {
	fprintf(file, GetParamFromType((GenFuncTable[nIndex].ParamObjType)[j]));
	if ((GenFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && bToIFile)
	    fprintf(file, "* INPUT");
	else if ((GenFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && !bToIFile)
	    fprintf(file, "* Param%d", j + 1);
	else
	    fprintf(file, " Param%d", j + 1);
	fprintf(file, ", ");
    }
    if ((GenFuncTable[nIndex].ParamObjType)[j] != 0) {
	fprintf(file, GetParamFromType((GenFuncTable[nIndex].ParamObjType)[j]));
	if ((GenFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && bToIFile)
	    fprintf(file, "* INPUT");
	else if ((GenFuncTable[nIndex].ParamObjType)[j] == NUMERIC_EXPR && !bToIFile)
	    fprintf(file, "* Param%d", j + 1);
	else
	    fprintf(file, " Param%d", j + 1);
    }
    fprintf(file, ")");    
}

void PrintProceduresFuncDefinition(FILE* file, int nIndex)
{
    int j = 0;

    fprintf(file, "{\n");
    fprintf(file, "\t");
    fprintf(file, GenFuncTable[nIndex].CFuncName);
    fprintf(file, "(");
    for (j = 0; j < GenFuncTable[nIndex].NumOfParam - 1; j++) {
	fprintf(file, "Param%d", j + 1);
	fprintf(file, ", ");
    }
    if (GenFuncTable[nIndex].NumOfParam != 0) {
	fprintf(file, "Param%d", j + 1);
    }

    fprintf(file, ");\n}\n");    
}


void Irit2SwigAddons(FILE *i_file, FILE *input_file)
{
    char line[IRIT_TO_SWIG_LINE_LENGTH];

    while (fgets(line, IRIT_TO_SWIG_LINE_LENGTH, input_file) != NULL) {
	fprintf(i_file, "%s", line);
    }
}

/* Skip translating those functions, because we use python math library. */
int SkipNumeric(char* FuncName)
{    
    if ((stricmp("ACOS", FuncName) == 0) ||
	(stricmp("ASIN", FuncName) == 0) ||
	(stricmp("ATAN2", FuncName) == 0) ||
	(stricmp("ATAN", FuncName) == 0) ||
	(stricmp("COS", FuncName) == 0) ||
	(stricmp("EXP", FuncName) == 0) ||
	(stricmp("ABS", FuncName) == 0) ||
	(stricmp("FLOOR", FuncName) == 0) ||
	(stricmp("FMOD", FuncName) == 0) ||
	(stricmp("POWER", FuncName) == 0) ||
	(stricmp("LN", FuncName) == 0) ||
	(stricmp("LOG", FuncName) == 0) ||
	(stricmp("SIN", FuncName) == 0) ||
	(stricmp("SQRT", FuncName) == 0) ||
	(stricmp("TAN", FuncName) == 0)) {
	return 1;
    }
	
    return 0;
}
