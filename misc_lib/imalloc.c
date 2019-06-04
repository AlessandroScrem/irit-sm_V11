/******************************************************************************
* All mallocs from irit modules should be piped through this allocator.       *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
*					 Written by Gershon Elber, April 1993 *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(sgi) || defined(__WINNT__)
#include <malloc.h>
#endif /* sgi || __WINNT__ */
#ifdef __WINNT__
#include <crtdbg.h>
#ifdef IRIT_HAVE_SW_DBG
#include <StackWalkerLib.h>
#endif /* IRIT_HAVE_SW_DBG */
#endif /* __WINNT__ */
#include "inc_irit/irit_sm.h"
#include "misc_loc.h"

#define OVERWRITE_STAMP_START 0x12345678	/* An intgr32 */
#define OVERWRITE_STAMP_END1  0x92		/* An intgr32 */
#define OVERWRITE_STAMP_END4  0x301060bdL	/* An intgr32 */
#define OVERWRITE_STAMP_FREED 0xf0eef0eeL	/* An intgr32 */
#define IRIT_MALLOC_HASH_TABLE_SIZE 1048576
#define IRIT_MALLOC_HASH_ENTRY_SIZE 16

#define IRIT_DEBUG_MALLOC_NOT_MALLOCED		0x02
#define IRIT_DEBUG_MALLOC_FREE_TBL		0x04
#define IRIT_DEBUG_MALLOC_WNT_CRTDBG		0x08
#define IRIT_DEBUG_MALLOC_WNT_CRTDBG16		0x10
#define IRIT_DEBUG_MALLOC_STACK_TRACE		0x20

#ifdef IRIT_QUIET_STRINGS
#    define IRIT_ALLOC_ERROR(Msg, p)
#else
#    define IRIT_ALLOC_ERROR(Msg, p)  AllocError(Msg, p)
#endif /* IRIT_QUIET_STRINGS */

#define IRIT_MALLOC_DEBUG_HASH_KEY(p) \
	(((IritIntPtrSizeType) p) >> 4) & (IRIT_MALLOC_HASH_TABLE_SIZE - 1)

#ifdef DEBUG_IRIT_MALLOC
#define DEBUG_IRIT_MALLOC_PRINT
IRIT_STATIC_DATA int
    IritMallocInit = FALSE,
    IritDebugKeepStack = FALSE,
    IritDebugKeepStackStep = 1,
    IritMallocGlobalAllocID = 0,
    IritDebugMalloc = 0,
    IritDebugSearchID = 0;

#if defined(AMIGA) && defined(__SASC)
IRIT_STATIC_DATA VoidPtr __far
#else
IRIT_STATIC_DATA VoidPtr
#endif
    IritMallocDebugHashTable[IRIT_MALLOC_HASH_TABLE_SIZE]
			    [IRIT_MALLOC_HASH_ENTRY_SIZE];

static void AllocError(char *Msg, VoidPtr p);
static int Env2Int(char *EnvStr);
static int IritTestOneAllocPtr(VoidPtr p);
#endif /* DEBUG_IRIT_MALLOC */

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Tests the content of the dynamic memory allocated.                       M
*                                                                            *
* PARAMETERS:                                                                M
*   None								     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritMalloc, IritRealloc, IritInitTestDynMemory		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritTestAllDynMemory                                                     M
*****************************************************************************/
void IritTestAllDynMemory(void)
{
#ifdef DEBUG_IRIT_MALLOC
    int i, j;

    for (i = 0; i < IRIT_MALLOC_HASH_TABLE_SIZE; i++) {
        for (j = 0; j < IRIT_MALLOC_HASH_ENTRY_SIZE; j++) {
	    if (IritMallocDebugHashTable[i][j] != NULL) {
	        char
		    *q = (char *) IritMallocDebugHashTable[i][j];

		IritTestOneAllocPtr(q);
	    }
	}
    }
#endif /* DEBUG_IRIT_MALLOC */
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets/clear dynamic memory check marks and reports.  Initializes the      M
* dynamic memory testing routines, if not initialized until now.	     M
*                                                                            *
* PARAMETERS:                                                                M
*   Start:	1 to initial check mark, 2 to also request stack trace info, M
*		0 to terminate and dump result.				     M
*   KeepStackStep: If Start = 2, sets stack-keeping every n'th alloc.        M 
*   TrackAllocID:  if non-zero, also track the unique allocation ID.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritMalloc, IritRealloc, IritTestAllDynMemory		     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritCheckMarkDynMemory                                                   M
*****************************************************************************/
void IritCheckMarkDynMemory(IrtRType *Start,
			    IrtRType *KeepStackStep,
			    IrtRType *TrackAllocID)
{
#ifdef DEBUG_IRIT_MALLOC
    if (!IRIT_APX_EQ(*Start, 0.0)) {
	IRIT_ZAP_MEM(IritMallocDebugHashTable, sizeof(IritMallocDebugHashTable));
	IritDebugMalloc |= IRIT_DEBUG_MALLOC_FREE_TBL;
	IritDebugKeepStack = IRIT_APX_EQ(*Start, 2.0);
	IritDebugKeepStackStep = (int) *KeepStackStep;

	if (!IRIT_APX_EQ(*Start, 0.0)) {
	    IritDebugSearchID = (int) *TrackAllocID;


	    printf(IRIT_EXP_STR("Tracking allocation ID = %d"),
		   IritDebugSearchID);
	    if (IritDebugSearchID < 0) {
		printf(IRIT_EXP_STR(" - negative ID is invalid. Enter to continue:"));
		getchar();
	    }
	    else
		printf(IRIT_EXP_STR("\n"));
	}
    }
    else {
	int i, j,
	    n = 0;

	IritDebugKeepStack = FALSE;
	IritDebugMalloc = 0;

	for (i = 0; i < IRIT_MALLOC_HASH_TABLE_SIZE; i++) {
	    for (j = 0; j < IRIT_MALLOC_HASH_ENTRY_SIZE; j++) {
		if (IritMallocDebugHashTable[i][j] != NULL) {
		    AllocError("Unfreed object",
			       IritMallocDebugHashTable[i][j]);
		    n++;
		}
	    }
	}

	if (n > 0)
	    printf("\n%d blocks were detected as unfreed (out of %d allocations)\n",
		   n, IritMallocGlobalAllocID);

	if (IritDebugSearchID != 0)
	    printf("Tracking Alloacted ID %d\n", IritDebugSearchID);
    }
#endif /* DEBUG_IRIT_MALLOC */
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Initialize dynamic memory testing routines, using envvars.               M
*                                                                            *
* PARAMETERS:                                                                M
*   None                                                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritMalloc, IritRealloc, IritCheckMarkDynMemory		     M
*   IritInitTestDynMemory2						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritInitTestDynMemory                                                    M
*****************************************************************************/
void IritInitTestDynMemory(void)
{
#ifdef DEBUG_IRIT_MALLOC
    if (!IritMallocInit) {
        int DebugMalloc, DebugSearchAllocID;

        DebugMalloc = Env2Int(IRIT_EXP_STR("IRIT_MALLOC"));
        DebugSearchAllocID = Env2Int(IRIT_EXP_STR("IRIT_MALLOC_ID"));

	IritInitTestDynMemory2(DebugMalloc, DebugSearchAllocID);

	IritMallocInit = TRUE;
    }
#endif /* DEBUG_IRIT_MALLOC */
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Initialize dynamic memory testing routines.                              M
*                                                                            *
* PARAMETERS:                                                                M
*   DebugMalloc:         Bitwise control over what to test:		     M
*			 0x02 - List possibly freed yet not malloced ptrs.   M
*			 0x04 - Track unfreed memory blocks.		     M
*			 0x08 - Under windows use windows CrtDbg check.      M
*			 0x10 - Under windows use windows CrtDbg check every M
*			        16 mallocs.				     M
*			 0x20 - Under windows, enable stack info trace.      M
*                        See IRIT_DEBUG_MALLOC_* at the beginning of file.   M
*   DebugSearchAllocID:  Allocation ID to trace, 0 for no ID to trace.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritMalloc, IritRealloc, IritCheckMarkDynMemory		     M
*   IritInitTestDynMemory						     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritInitTestDynMemory2                                                   M
*****************************************************************************/
void IritInitTestDynMemory2(int DebugMalloc, int DebugSearchAllocID)
{
#ifdef DEBUG_IRIT_MALLOC
    IritDebugMalloc = DebugMalloc;
    IritDebugSearchID = DebugSearchAllocID;

#   ifdef __WINNT__
    IritDebugKeepStack = (IritDebugMalloc &
			  IRIT_DEBUG_MALLOC_STACK_TRACE) != 0;

    if (IritDebugMalloc & IRIT_DEBUG_MALLOC_WNT_CRTDBG ||
	IritDebugMalloc & IRIT_DEBUG_MALLOC_WNT_CRTDBG16) {
        int tmpDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);

	if (IritDebugMalloc & IRIT_DEBUG_MALLOC_WNT_CRTDBG)
	    tmpDbgFlag |= _CRTDBG_CHECK_ALWAYS_DF;
#   if _MSC_VER >= 1300                    /* VC++ 7.x or above. */
	else
	    tmpDbgFlag |= _CRTDBG_CHECK_EVERY_16_DF;
#   else
	else
	    tmpDbgFlag |= _CRTDBG_CHECK_ALWAYS_DF;
#   endif /* _MSC_VER */
	_CrtSetDbgFlag(tmpDbgFlag);

	printf(IRIT_EXP_STR("WINNT enabled CrtCheckMemory every malloc/free\n"));
    }
#   endif /* __WINNT__ */

    if (IritDebugSearchID != 0)
        printf(IRIT_EXP_STR("Tracking malloc ID %d\n"), IritDebugSearchID);

    if (IritDebugKeepStack)
        printf(IRIT_EXP_STR("Keeping stack info on every %d mallocs\n"),
	       IritDebugKeepStackStep);

    IRIT_ZAP_MEM(IritMallocDebugHashTable, sizeof(IritMallocDebugHashTable));
#endif /* DEBUG_IRIT_MALLOC */
}

#ifdef DEBUG_IRIT_MALLOC

/*****************************************************************************
* DESCRIPTION:                                                               *
* Routine to print an allocation error happened where pointer p is involved. *
*                                                                            *
* PARAMETERS:                                                                *
*   Msg:         A description of the error.                                 *
*   p:           Pointer that is involved in the error.                      *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void AllocError(char *Msg, VoidPtr p)
{
    unsigned char
	**q = (unsigned char **) p;
    IritIntPtrSizeType
        *r = (IritIntPtrSizeType *) &q[-4];		       /* 16 bytes. */

    printf("Memory alloc error (%d bytes), Ptr = 0x%lx (ID = %d): %s\n",
	   r[0], (IritIntPtrSizeType) p, r[1], Msg);

    /* Get the allocation location - file name and line number. */
    if (((IritIntPtrSizeType) q[-1]) == OVERWRITE_STAMP_START) {
        if (q[-5] == NULL)
	    printf("    \"%s\" at %s:%d\n",
		   q[-8], q[-7], (IritIntPtrSizeType) q[-6]);
	else {
	    printf("    \"%s\" at %s:%d\n    STACK:%s",
		   q[-8], q[-7], (IritIntPtrSizeType) q[-6], q[-5]);
	    if (IritDebugKeepStackStep > 1)
	        IritSleep(1000);
	}
    }
    else {
        printf("\nMemory around allocation is corrupted!\n");
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   COnverts a given environment string to integer.  The environment string  *
* could be in decimal or hex base in which it is prefixed with "0x".	     *
*                                                                            *
* PARAMETERS:                                                                *
*   EnvStr:    Environment string to convert to integer.                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:      Converted value.                                               *
*****************************************************************************/
static int Env2Int(char *EnvStr)
{
    int i;

    if (IRT_STR_NULL_ZERO_LEN(EnvStr) || (EnvStr = getenv(EnvStr)) == NULL)
	return 0;

    if (strncmp(EnvStr, "0x", 2) == 0) {
	if (sscanf(EnvStr, "%x", &i) == 1)
	    return i;
    }
    else {
	if (sscanf(EnvStr, "%d", &i) == 1)
	    return i;
    }

    return 0;
}
#endif /* DEBUG_IRIT_MALLOC */

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to reallocate dynamic memory for all IRIT program/tool/libraries.  M
*    All requests for dynamic memory should invoke this function.	     M
*    If the environment variable "IRIT_MALLOC" is set when an IRIT program   M
* is executed, the consistency of the dynamic memory is tested on every      M
* invokation of this routine. See IritTestAllDynMemory function for more.    M
*                                                                            *
* PARAMETERS:                                                                M
*   p:        Old pointer to reallocate.				     M
*   OldSize:  Size of old block pointed by p, zero if unknown.               M
*   NewSize:  Size of new block to allocate, in bytes.  Must be larger than  M
*	      OldSize.							     M
*                                                                            *
* RETURN VALUE:                                                              M
*   VoidPtr:  A pointer to the allocated block. A function calling this      M
*             may assume return value will never be NULL, since no more      M
*             memory cases are trapped locally.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritMalloc, IritCheckMarkDynMemory, IritInitTestDynMemory      M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritRealloc, allocation                                                  M
*****************************************************************************/
VoidPtr IritRealloc(VoidPtr p, unsigned OldSize, unsigned NewSize)
{
    VoidPtr
	NewP = IritMalloc(NewSize);

    assert(NewSize > OldSize);
    IRIT_GEN_COPY(NewP, p, OldSize ? OldSize : NewSize);
    
    IritFree(p);

    return NewP;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to allocate dynamic memory for all IRIT program/tool/libraries.    M
*    All requests for dynamic memory should invoke this function.	     M
*    If the environment variable "IRIT_MALLOC" is set when an IRIT program   M
* is executed, the consistency of the dynamic memory is tested on every      M
* invokation of this routine. See IritTestAllDynMemory function for more.    M
*                                                                            *
* PARAMETERS:                                                                M
*   Size:     Size of block to allocate, in bytes.                           M
*   ObjType:  This variable exists iff "#define DEBUG_IRIT_MALLOC".	     M
*	      This holds the object descriptions.			     M
*   FileName: This variable exists iff "#define DEBUG_IRIT_MALLOC".	     M
*	      This holds the file name where the call is invoked from.	     M
*   LineNum:  This variable exists iff "#define DEBUG_IRIT_MALLOC".	     M
*	      This holds the line number where the call is invoked from.     M
*                                                                            *
* RETURN VALUE:                                                              M
*   VoidPtr:  A pointer to the allocated block. A function calling this      M
*             may assume return value will never be NULL, since no more      M
*             memory cases are trapped locally.				     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritRealloc, IritCheckMarkDynMemory, IritInitTestDynMemory     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritMalloc, allocation                                                   M
*****************************************************************************/
#ifdef DEBUG_IRIT_MALLOC
#undef IritMalloc
VoidPtr IritMalloc(unsigned Size,
		   const char *ObjType,
		   const char *FileName,
		   int LineNum)
{
    VoidPtr
	p = malloc(Size + 40);
    char *r0;
    const char **q = (const char **) p;
    int *r;

    if (p == NULL)
	MISC_FATAL_ERROR(MISC_ERR_MALLOC_FAILED, "");

    /* First keep the allocation location - file name and line number. */
    q[0] = ObjType;
    q[1] = FileName;
    q[2] = (const char *) ((IritIntPtrSizeType) LineNum);
#   if defined(__WINNT__) && defined(IRIT_HAVE_SW_DBG)
	if (IritDebugKeepStack) {
	    IRIT_STATIC_DATA int
		KeepStackCount = 0;

	    if (KeepStackCount-- <= 0) {
	        q[3] = strdup(IritGetStackState());
		KeepStackCount = IritDebugKeepStackStep;
	    }
	    else
		q[3] = NULL;
	}
	else
	    q[3] = NULL;
#   else
	q[3] = NULL;
#   endif /* __WINNT__ && IRIT_HAVE_SW_DBG */
    r = (int *) &q[4];

    /* Now keep original requested size and put stamps at both ends. */
    r[0] = (int) Size;
    r[1] = ++IritMallocGlobalAllocID;
    r[2] = OVERWRITE_STAMP_START;
    r[3] = OVERWRITE_STAMP_START;
    p = &r[4]; /* Skip 32 bytes from malloc'ed pointer - align x8 ! */

    /* At the far end pad with bytes until we get to a x4 alignment. */
    r0 = (char *) (((IritIntPtrSizeType) p) + Size);
    r = (int *) (((((IritIntPtrSizeType) p) + Size + 3) >> 2) << 2);
    while (r0 != (char *) r)
        *r0++ = (char) OVERWRITE_STAMP_END1;
    r[0] = OVERWRITE_STAMP_END4;

    /* Lets see if we seek this specific pointer. */
    if (IritDebugSearchID && IritMallocGlobalAllocID == IritDebugSearchID) {
#	ifdef DEBUG_IRIT_MALLOC_PRINT
	    printf("Pointer 0x%08lx just allocated (ID = %d)\n",
		   (IritIntPtrSizeType) p, IritDebugSearchID);
#	    if defined(__WINNT__) && defined(IRIT_HAVE_SW_DBG)
		IritPrintStackState();
#	    endif /* __WINNT__ && IRIT_HAVE_SW_DBG */
	    abort();
#	endif /* DEBUG_IRIT_MALLOC_PRINT */
    }

    /* Place allocated pointer in global hash table, if so requested. */
    if (IritDebugMalloc & IRIT_DEBUG_MALLOC_FREE_TBL) {
        int i;
	IritIntPtrSizeType
	    Key = IRIT_MALLOC_DEBUG_HASH_KEY(p);
	VoidPtr
	    *v = IritMallocDebugHashTable[Key];

	for (i = 0; i < IRIT_MALLOC_HASH_ENTRY_SIZE; i++, v++) {
	    if (*v == NULL) {
	        *v = p;
		break;
	    }
	}
	if (i >= IRIT_MALLOC_HASH_ENTRY_SIZE)
	    printf( "Irit malloc debug hash table is full pointer 0x%08lx allocated\n",
		    (IritIntPtrSizeType) p);
    }

    *((char *) p) = 0; /* Dont leave OVERWRITE_STAMP_FREED intact. */

    return p;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Set the searched malloced ID. This function will take affect iff         M
* IritDebugMalloc is non zero as set via "IRIT_MALLOC" env.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   ID:		Allocation ID to seek and abort at.                          M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritFree, IritMalloc, IritRealloc, IritTestAllDynMemory,		     M
*   IritDebugMallocSearchPtr, IritDebugMallocAllocated      		     M
*   IritInitTestDynMemory, IritCheckMarkDynMemory			     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritDebugMallocSearchID                                                  M
*****************************************************************************/
void IritDebugMallocSearchID(int ID)
{
    IritDebugSearchID = ID;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Routine to free dynamic memory for all IRIT program/tool/libraries.        M
*    All requests to free dynamic memory should invoke this function.        M
*                                                                            *
* PARAMETERS:                                                                M
*   p:          Pointer to a block that needs to be freed.                   M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritMalloc, IritRealloc, IritCheckMarkDynMemory, IritInitTestDynMemory   M
*                                                                            *
* KEYWORDS:                                                                  M
*   IritFree, allocation                                                     M
*****************************************************************************/
#undef IritFree
void IritFree(VoidPtr p)
{
    unsigned char **r,
        *OldP = p;
    int ID;

#ifdef __WINNT__
    assert(((unsigned int) p) != 0xfeeefeee); /* On Win an invalid pointer. */
#endif /* __WINNT__ */

    if (p == NULL)
	return;

    ID = IritTestOneAllocPtr(p);

    r = (unsigned char **) p;
    if (r[-5] != NULL)
        free(r[-5]);                                         /* Stack info. */

    /* Lets see if we seek this specific pointer. */
    if (IritDebugSearchID && ID == IritDebugSearchID) {
#	ifdef DEBUG_IRIT_MALLOC_PRINT
	    printf("Pointer 0x%08lx just freed (ID = %d (current ID = %d))\n",
		   (IritIntPtrSizeType) p, ID, IritMallocGlobalAllocID);
#	    if defined(__WINNT__) && defined(IRIT_HAVE_SW_DBG)
		IritPrintStackState();
#	    endif /* __WINNT__ && IRIT_HAVE_SW_DBG */
	    abort();
#	endif /* DEBUG_IRIT_MALLOC_PRINT */
    }

    /* Mark as freed. */
    *((int *) p) = OVERWRITE_STAMP_FREED;

    free(&OldP[-32]);/* We use 32 bytes (8 ptrs) before returned malloc ptr. */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Test the validity of the given pointer.                                  *
*                                                                            *
* PARAMETERS:                                                                *
*   p:   Allocated pointer to test its validity.                             *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:  0 if invalid, alloc index otherwise.                               *
*****************************************************************************/
static int IritTestOneAllocPtr(VoidPtr p)
{
    unsigned char *r0,
        *OldP = p;
    IritIntPtrSizeType *r, Size, AllocIdx;

    if (p == NULL)
	return TRUE;

    /* Examine for existence in the allocated table. */
    if (IritDebugMalloc & IRIT_DEBUG_MALLOC_FREE_TBL) {
        int i;
	IritIntPtrSizeType
	    Key = IRIT_MALLOC_DEBUG_HASH_KEY(p);
	VoidPtr
	    *v = IritMallocDebugHashTable[Key];

	for (i = 0; i < IRIT_MALLOC_HASH_ENTRY_SIZE; i++, v++) {
	    if (*v == p)
	        break;
	}
	if (i >= IRIT_MALLOC_HASH_ENTRY_SIZE) {
	    if (IritDebugMalloc & IRIT_DEBUG_MALLOC_NOT_MALLOCED)
	        printf("A pointer that was never allocated, Ptr = 0x%lx (%lu)\n",
		       OldP, OldP);
	    return FALSE;
	}
	else
	    *v = NULL;    
    }

    /* Get the size and skip file name/line num/stack of allocation point. */
    r = (IritIntPtrSizeType *) &OldP[-16];
    Size = r[0];
    AllocIdx = r[1];

    /* Examine for an already freed data structure. */
    if (r[3] == OVERWRITE_STAMP_FREED) {
        IRIT_ALLOC_ERROR("Using an already freed object", OldP);
        return 0;
    }

    /* Examine for overwriting at the beginning. */
    if (r[2] != OVERWRITE_STAMP_START || r[3] != OVERWRITE_STAMP_START) {
        IRIT_ALLOC_ERROR("Dymanic memory overwrite before begining of allocated block", OldP);
        return 0;
    }

    /* Examine for overwriting at the end until we get to a x4 alignment. */
    r0 = ((unsigned char *) p) + Size;
    r = (IritIntPtrSizeType *)
                           (((((IritIntPtrSizeType) p) + Size + 3) >> 2) << 2);
    while (r0 != ((unsigned char *) r))
        if (*r0++ != OVERWRITE_STAMP_END1) {
	    IRIT_ALLOC_ERROR("Dymanic memory overwrite after end of allocated block", OldP);
	    return 0;
	}
    if (r[0] != OVERWRITE_STAMP_END4) {
	IRIT_ALLOC_ERROR("Dymanic memory overwrite after end of allocated block", OldP);
        return 0;
    }

    return AllocIdx;
}

#else

/* If we build under non 32 bit machine define them so the DLL prescribed by */
/* irit.def and/or iritD.def will link properly.			     */

#ifdef IritMalloc
#undef IritMalloc
#endif /* IritMalloc */
VoidPtr IritMalloc(unsigned Size,
		   const char *ObjType,
		   const char *FileName,
                   int LineNum)
{
    IRIT_FATAL_ERROR("No support for dynamic memory testing in 64 bits.");
    return NULL;
}
void IritDebugMallocSearchID(int ID)
{
    IRIT_FATAL_ERROR("No support for dynamic memory testing in 64 bits.");
}
#ifdef IritFree
#undef IritFree
#endif /* IritFree */
void IritFree(VoidPtr p)
{
    IRIT_FATAL_ERROR("No support for dynamic memory testing in 64 bits.");
}

#endif /* DEBUG_IRIT_MALLOC */
