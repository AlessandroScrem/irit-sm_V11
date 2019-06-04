/*****************************************************************************
*   "Irit" - the 3d (not only polygonal) solid modeller.		     *
*									     *
* Written by:  Gershon Elber				Ver 0.2, Mar. 1990   *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Main definition Header file for Irit - the 3d (not only polygonal) sm.     *
*****************************************************************************/

#ifndef	PROGRAM_H
#define	PROGRAM_H

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include <setjmp.h>	/* Used for the long jumps - to main iteration loop. */
#include "inc_irit/irit_sm.h"
#include "inc_irit/allocate.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/attribut.h"
#include "inc_irit/bool_lib.h"
#include "inc_irit/cagd_lib.h"	   /* We define curves and surfaces handles as well. */
#include "inc_irit/symb_lib.h"
#include "inc_irit/trim_lib.h"
#include "inc_irit/triv_lib.h"
#include "inc_irit/mvar_lib.h"

#define INPUT_LINE_LEN 10000   		/* Input parser maximum line length. */

#define DEFAULT_LOAD_COLOR 1  /* Default colors for object loaded using LOAD */
#define DEFAULT_PRIM_COLOR 4  /* primitives colors, respectively.	     */

#define DEFAULT_RESOLUTION	20	/* Used in Primitiv/Boolean modules. */
#define DEFAULT_DRAW_CTLPT	FALSE	/* If control mesh/poly to be drawn. */
#define DEFAULT_ECHOSRC		TRUE	            /* Echo sourced program? */
#define DEFAULT_DUMPLVL		1		        /* Min. information. */

#define IRIT_MACHINE_MSDOS		1
#define IRIT_MACHINE_SGI		2
#define IRIT_MACHINE_HP			3
#define IRIT_MACHINE_APOLLO		4
#define IRIT_MACHINE_SUN		5
#define IRIT_MACHINE_UNIX		6
#define IRIT_MACHINE_IBMOS2		7
#define IRIT_MACHINE_WINDOWS		8
#define IRIT_MACHINE_AMIGA		9
#define IRIT_MACHINE_CYGWIN		10
#define IRIT_MACHINE_MACOSX		11
#define IRIT_MACHINE_LINUX		12

#define KV_MIN_LEGAL		-9999
#define KV_UNIFORM_OPEN		-10000
#define KV_UNIFORM_FLOAT	-10001
#define KV_UNIFORM_PERIODIC	-10002
#define KV_UNIFORM_DISCONT_OPEN	-10003

#define IRIT_PROMPT	"Irit> "

#define GUI_PROMPT	"<<GUI>> "

#ifdef IRIT_QUIET_STRINGS
#define IRIT_WNDW_FPRINTF1(p1)
#define IRIT_WNDW_FPRINTF2(p1, p2)
#define IRIT_WNDW_FPRINTF3(p1, p2, p3)
#define IRIT_WNDW_FPRINTF4(p1, p2, p3, p4)
#define IRIT_WNDW_FPRINTF5(p1, p2, p3, p4, p5)
#define IRIT_WNDW_PUT_STR(Str)
#define IRIT_WNDW_PUT_STR2(Str)
#define IRIT_NON_FATAL_ERROR(Str)			IritQuietError()
#define IRIT_NON_FATAL_ERROR2(Str, p2)			IritQuietError()
#define IRIT_NON_FATAL_ERROR3(Str, p2, p3)		IritQuietError()
#define IRIT_NON_FATAL_ERROR4(Str, p2, p3, p4)		IritQuietError()
#define IRIT_NON_FATAL_ERROR5(Str, p2, p3, p4, p5)	IritQuietError()
#else
#define IRIT_WNDW_FPRINTF1		IritFprintf
#define IRIT_WNDW_FPRINTF2		IritFprintf
#define IRIT_WNDW_FPRINTF3		IritFprintf
#define IRIT_WNDW_FPRINTF4		IritFprintf
#define IRIT_WNDW_FPRINTF5		IritFprintf
#define IRIT_WNDW_PUT_STR(Str)		IritPutStr(Str)
#define IRIT_WNDW_PUT_STR2(Str)		IritPutStr2(Str)
#define IRIT_NON_FATAL_ERROR		IritNonFatalError
#define IRIT_NON_FATAL_ERROR2		IritNonFatalError
#define IRIT_NON_FATAL_ERROR3		IritNonFatalError
#define IRIT_NON_FATAL_ERROR4		IritNonFatalError
#define IRIT_NON_FATAL_ERROR5		IritNonFatalError
#endif /* IRIT_QUIET_STRINGS */

IRIT_GLOBAL_DATA_HEADER IPObjectStruct
    *GlblLastIritClientCursorEvent;

IRIT_GLOBAL_DATA_HEADER jmp_buf
    IritGlblLongJumpBuffer;		          /* Used in error recovery. */

IRIT_GLOBAL_DATA_HEADER FILE
    *GlblLogFile;		   /* If do log everything, it goes to here. */

IRIT_GLOBAL_DATA_HEADER BspMultComputationMethodType
    GlblBspProdMethod;

IRIT_GLOBAL_DATA_HEADER int
    GlblKeepClientsCursorEvents,
    GlblQuietMode,
    GlblCurrentDisplay,            /* Handlers to streams to display device. */
    GlblFlatLoadMode,          /* Hierarchy flattening in "load" operations. */
    GlblHandleDependencies,
    GlblPropagateNames, /* Update sons names of objs, when inserted into DB. */
    GlblMakeAllHierarchyVisible, /* Name all sons of new obj, when enter DB. */
    GlblStdinInteractive,
    GlblRunScriptAndQuit,
    GlblLoadColor,	      /* Default colors for object loaded using LOAD */
    GlblPrimColor,	       /* primitives colors, respectively.	     */
    GlblDoGraphics,		/* Control if running in graphics/text mode. */
    GlblGUIMode,				     /* Running under a GUI. */
    GlblFatalError,		  /* True if disaster in system - must quit! */
    GlblPrintLogFile,		     /* If TRUE everything goes to log file. */
    GlblPointLenAux;

IRIT_GLOBAL_DATA_HEADER char 
    *GlblHelpFileName,
    *GlblLineEditControl,
    *GlblStartFileName,			/* Name of startup file to executed. */
    *GlblLogFileName,					/* Name of log file. */
    *GlblFloatFormat;	      /* Controls the ways real numbers are printed. */

IRIT_GLOBAL_DATA_HEADER IrtRType
    GlblPointLen;			       /* Scaler for point if drawn. */

void IritIdleFunction(int MiliSeconds);
void IritPrintOneString(const char *Line);
void IritExit0(void);
void IritExit(int ExitCode);

/* Support functions. */

void IritExecOneLine(const char *Line);
IPObjectStruct *IritSetIritState(const char *Name, IPObjectStruct *Data);
void IritPrintInptPrsrError(void);
void IritInputWindowGetStr(char *Str, int Length);
void IritReset(void);

IritFatalMsgFuncType IritSetNonFatalErrorFunc(IritFatalMsgFuncType ErrorFunc);
#ifdef USE_VARARGS
void IritNonFatalError(const char *va_alist, ...);
#else
void IritNonFatalError(const char *Format, ...);
#endif /* USE_VARARGS */

void IritQuietError(void);
void IritDefaultBoolFatalError(BoolFatalErrorType ErrID);
void IritConvertData2Irt(IPObjectStruct *PObj);
void IritDefaultFPEHandler(int Type);
void IritAllFatalErrorTrap(void);
void IritDummyLinkDebug(void);

/* DB functions. */

int IritDBGetStackSize();
IPObjectStruct *IritDBGetStackObject(int Idx);
void IritDBFreeAll();
void IritDBFreeObject(IPObjectStruct *PObj);
void IritDBDeleteObject(IPObjectStruct *PObj, int Free);
void IritDBInsertObject(IPObjectStruct *PObj, int DelOld);
void IritDBInsertObjectLast(IPObjectStruct *PObj, int DelOld);
void IritDBInsertHierarchy(IPObjectStruct *PObj, int ByReference);
IPObjectStruct *IritDBGetObjByName(const char *Name);
IPObjectStruct *IritDBGetNextObj(IPObjectStruct *Crnt);
void IritDBValidateVariables(void);
void IritDBPrintAllObjs(void);
void IritDBPush(IPObjectStruct *NewDB);
void IritDBPop(int Free);
void IritDBUpdateParams(IPObjectStruct *PObj,
			const char *FuncID,
			IPObjectStruct *ParamList);

/* Prototypes of the I/O window functions: */

char *IritGetLineStdin(char *Line, const char *Prompt, int MaxLen);
void IritLogPrintSetup(IPObjectStruct *Set);
void IritPutStr(const char *Msg);
void IritPutStr2(const char *Msg);
#ifdef USE_VARARGS
void IritFprintf(const char *va_alist, ...);
#else
void IritFprintf(const char *Format, ...);
#endif /* USE_VARARGS */
IPObjectStruct *IritInputStdinObject(IrtRType *Type);
void IritDispViewSetDisplay(void);
void IritDispSetCrntDisplay(IrtRType *Ri);
void IritDispViewObject(IPObjectStruct *PObj);
void IritDispViewClearScreen(void);
void IritDispViewSaveMatrix(const char *FileName);
void IritDispViewDisconnect(void);
void IritDispViewExit(void);
double IritClientExecute(const char *PrgmName);
void IritClientClose(IrtRType *RHandler, IrtRType *RKillClient);
IPObjectStruct *IritClientRead(IrtRType *RHandler, IrtRType *RBlock);
IPObjectStruct *IritClientCursor(IrtRType *RWaitTime);
void IritClientWrite(IrtRType *RHandler, IPObjectStruct *PObj);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif	/* PROGRAM_H */
