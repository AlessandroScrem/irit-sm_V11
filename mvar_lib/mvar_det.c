/******************************************************************************
* Mvar_Det.c - Compute determinants of multi-variates.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber, May. 97.					      *
******************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "mvar_loc.h"

static void MvarGetSubMVMat(const MvarMVStruct * const *MVsMatrix, 
			    int MatSize,
			    const MvarMVStruct **SubMat,
			    int Row,
			    int Col,
			    CagdBType CreateCopy);

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes the expression of MV11 * MV22 - MV12 * MV21, which is a	     M
* determinant of a 2 by 2 matrix.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   MV11, MV12, MV21, MV22:  The four factors of the determinant.            M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A scalar field representing the determinant computation. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVDeterminant3, MvarSrfDeterminant2, MvarCrvDeterminant2	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVDeterminant2, determinant                                          M
*****************************************************************************/
MvarMVStruct *MvarMVDeterminant2(const MvarMVStruct *MV11,
				 const MvarMVStruct *MV12,
				 const MvarMVStruct *MV21,
				 const MvarMVStruct *MV22)
{
    MvarMVStruct
	*Prod1 = MvarMVMult(MV11, MV22),
	*Prod2 = MvarMVMult(MV21, MV12),
	*Add12 = MvarMVSub(Prod1, Prod2);

    MvarMVFree(Prod1);
    MvarMVFree(Prod2);
    return Add12;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes the expression of a 3 by 3 determinants.                          M
*                                                                            *
* PARAMETERS:                                                                M
*   MV11, MV12, MV13:  The nine factors of the determinant.                  M
*   MV21, MV22, MV23:			"				     M
*   MV31, MV32, MV33:			"				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A scalar field representing the determinant computation. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVDeterminant2, MvarMVDeterminant4, SymbSrfDeterminant3              M
*   SymbCrvDeterminant3							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVDeterminant3, determinant                                          M
*****************************************************************************/
MvarMVStruct *MvarMVDeterminant3(const MvarMVStruct *MV11,
				 const MvarMVStruct *MV12,
				 const MvarMVStruct *MV13,
				 const MvarMVStruct *MV21,
				 const MvarMVStruct *MV22,
				 const MvarMVStruct *MV23,
				 const MvarMVStruct *MV31,
				 const MvarMVStruct *MV32,
				 const MvarMVStruct *MV33)
{
    MvarMVStruct
	*Prod1 = MvarMVDeterminant2(MV22, MV23, MV32, MV33),
        *Prod1a = MvarMVMult(MV11, Prod1),
	*Prod2 = MvarMVDeterminant2(MV21, MV23, MV31, MV33),
        *Prod2a = MvarMVMult(MV12, Prod2),
	*Prod3 = MvarMVDeterminant2(MV21, MV22, MV31, MV32),
        *Prod3a = MvarMVMult(MV13, Prod3),
	*Sub12 = MvarMVSub(Prod1a, Prod2a),
	*Add123 = MvarMVAdd(Sub12, Prod3a);

    MvarMVFree(Prod1);
    MvarMVFree(Prod1a);
    MvarMVFree(Prod2);
    MvarMVFree(Prod2a);
    MvarMVFree(Prod3);
    MvarMVFree(Prod3a);
    MvarMVFree(Sub12);

    return Add123;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes the expression of a 4 by 4 determinants.                          M
*                                                                            *
* PARAMETERS:                                                                M
*   MV11, MV12, MV13, MV14:  The 16 factors of the determinant.              M
*   MV21, MV22, MV23, MV24:		    "				     M
*   MV31, MV32, MV33, MV34:		    "				     M
*   MV41, MV42, MV43, MV44:		    "				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A scalar field representing the determinant computation. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVDeterminant3, MvarMVDeterminant5				     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVDeterminant4, determinant                                          M
*****************************************************************************/
MvarMVStruct *MvarMVDeterminant4(const MvarMVStruct *MV11,
				 const MvarMVStruct *MV12,
				 const MvarMVStruct *MV13,
				 const MvarMVStruct *MV14,
				 const MvarMVStruct *MV21,
				 const MvarMVStruct *MV22,
				 const MvarMVStruct *MV23,
				 const MvarMVStruct *MV24,
				 const MvarMVStruct *MV31,
				 const MvarMVStruct *MV32,
				 const MvarMVStruct *MV33,
				 const MvarMVStruct *MV34,
				 const MvarMVStruct *MV41,
				 const MvarMVStruct *MV42,
				 const MvarMVStruct *MV43,
				 const MvarMVStruct *MV44)
{
    MvarMVStruct
	*Prod1 = MvarMVDeterminant3(MV22, MV23, MV24,
				    MV32, MV33, MV34,
				    MV42, MV43, MV44),
        *Prod1a = MvarMVMult(MV11, Prod1),
	*Prod2 = MvarMVDeterminant3(MV21, MV23, MV24,
				    MV31, MV33, MV34,
				    MV41, MV43, MV44),
        *Prod2a = MvarMVMult(MV12, Prod2),
	*Prod3 = MvarMVDeterminant3(MV21, MV22, MV24,
				    MV31, MV32, MV34,
				    MV41, MV42, MV44),
        *Prod3a = MvarMVMult(MV13, Prod3),
	*Prod4 = MvarMVDeterminant3(MV21, MV22, MV23,
				    MV31, MV32, MV33,
				    MV41, MV42, MV43),
        *Prod4a = MvarMVMult(MV14, Prod4),
	*Sub12 = MvarMVSub(Prod1a, Prod2a),
	*Sub34 = MvarMVSub(Prod3a, Prod4a),
	*Add1234 = MvarMVAdd(Sub12, Sub34);

    MvarMVFree(Prod1);
    MvarMVFree(Prod1a);
    MvarMVFree(Prod2);
    MvarMVFree(Prod2a);
    MvarMVFree(Prod3);
    MvarMVFree(Prod3a);
    MvarMVFree(Prod4);
    MvarMVFree(Prod4a);

    MvarMVFree(Sub12);
    MvarMVFree(Sub34);

    return Add1234;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Computes the expression of a 5 by 5 determinants.                          M
*                                                                            *
* PARAMETERS:                                                                M
*   MV11, MV12, MV13, MV14, MV15:  The 25 factors of the determinant.        M
*   MV21, MV22, MV23, MV24, MV25:		    "			     M
*   MV31, MV32, MV33, MV34, MV35:		    "			     M
*   MV41, MV42, MV43, MV44, MV45:		    "			     M
*   MV51, MV52, MV53, MV54, MV55:		    "			     M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A scalar field representing the determinant computation. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVDeterminant4							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVDeterminant5, determinant                                          M
*****************************************************************************/
MvarMVStruct *MvarMVDeterminant5(const MvarMVStruct *MV11,
				 const MvarMVStruct *MV12,
				 const MvarMVStruct *MV13,
				 const MvarMVStruct *MV14,
				 const MvarMVStruct *MV15,
				 const MvarMVStruct *MV21,
				 const MvarMVStruct *MV22,
				 const MvarMVStruct *MV23,
				 const MvarMVStruct *MV24,
				 const MvarMVStruct *MV25,
				 const MvarMVStruct *MV31,
				 const MvarMVStruct *MV32,
				 const MvarMVStruct *MV33,
				 const MvarMVStruct *MV34,
				 const MvarMVStruct *MV35,
				 const MvarMVStruct *MV41,
				 const MvarMVStruct *MV42,
				 const MvarMVStruct *MV43,
				 const MvarMVStruct *MV44,
				 const MvarMVStruct *MV45,
				 const MvarMVStruct *MV51,
				 const MvarMVStruct *MV52,
				 const MvarMVStruct *MV53,
				 const MvarMVStruct *MV54,
				 const MvarMVStruct *MV55)
{
    MvarMVStruct
	*Prod1 = MvarMVDeterminant4(MV22, MV23, MV24, MV25,
				    MV32, MV33, MV34, MV35,
				    MV42, MV43, MV44, MV45,
				    MV52, MV53, MV54, MV55),
        *Prod1a = MvarMVMult(MV11, Prod1),
	*Prod2 = MvarMVDeterminant4(MV21, MV23, MV24, MV25,
				    MV31, MV33, MV34, MV35,
				    MV41, MV43, MV44, MV45,
				    MV51, MV53, MV54, MV55),
        *Prod2a = MvarMVMult(MV12, Prod2),
	*Prod3 = MvarMVDeterminant4(MV21, MV22, MV24, MV25,
				    MV31, MV32, MV34, MV35,
				    MV41, MV42, MV44, MV45,
				    MV51, MV52, MV54, MV55),
        *Prod3a = MvarMVMult(MV13, Prod3),
	*Prod4 = MvarMVDeterminant4(MV21, MV22, MV23, MV25,
				    MV31, MV32, MV33, MV35,
				    MV41, MV42, MV43, MV45,
				    MV51, MV52, MV53, MV55),
        *Prod4a = MvarMVMult(MV14, Prod4),
	*Prod5 = MvarMVDeterminant4(MV21, MV22, MV23, MV24,
				    MV31, MV32, MV33, MV34,
				    MV41, MV42, MV43, MV44,
				    MV51, MV52, MV53, MV54),
        *Prod5a = MvarMVMult(MV15, Prod5),
	*Sub12 = MvarMVSub(Prod1a, Prod2a),
	*Sub34 = MvarMVSub(Prod3a, Prod4a),
	*Add1234 = MvarMVAdd(Sub12, Sub34),
	*Sub12345 = MvarMVAdd(Add1234, Prod5a);

    MvarMVFree(Prod1);
    MvarMVFree(Prod1a);
    MvarMVFree(Prod2);
    MvarMVFree(Prod2a);
    MvarMVFree(Prod3);
    MvarMVFree(Prod3a);
    MvarMVFree(Prod4);
    MvarMVFree(Prod4a);
    MvarMVFree(Prod5);
    MvarMVFree(Prod5a);

    MvarMVFree(Sub12);
    MvarMVFree(Sub34);
    MvarMVFree(Add1234);

    return Sub12345;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Computes the expression of a determinant of MVs, recursively.            M
* Expansion is always using the first row. Does not check for a better	     M
* option, such as a row with more zeros, etc. (Existing zeros means the zero M
* function as a multivariate).						     M
*                                                                            *
* PARAMETERS:                                                                M
*   MVsMatrix:  The matrix of MVs, as a linear array, with the following     M
*		double index to linear index convention -		     M
*	        [i][j] <----> [(i - 1) * MatSize + (j - 1)],		     M
*		where i, j = 1, 2, ..., MatSize (double index as in the      M
*		usual linear algebra convention).			     M
*   MatSize:    The size in each of the dimensions of the square MVs matrix. M
*                                                                            *
* RETURN VALUE:                                                              M
*   MvarMVStruct *: A scalar field representing the determinant computation. M
*                                                                            *
* SEE ALSO:                                                                  M
*   MvarMVDeterminant2							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   MvarMVDeterminant, determinant                                           M
*****************************************************************************/
MvarMVStruct *MvarMVDeterminant(const MvarMVStruct * const *MVsMatrix,
				int MatSize)
{
    assert(MatSize >= 1);

    if (MatSize == 1)
	return MvarMVCopy(MVsMatrix[0]);
    else if (MatSize == 2)
	return MvarMVDeterminant2(MVsMatrix[0], MVsMatrix[1],
				  MVsMatrix[2], MVsMatrix[3]);
    else {
	int i,
	    Sign = 1,
	    SubMatSize = MatSize - 1;
 	MvarMVStruct *TmpSubDetMV, *TmpMVToAdd, *TmpSumMV,
	    *SumMV = NULL;
	MvarMVStruct const
	    **TmpMVsSubMat = (MvarMVStruct const **) IritMalloc(
			     sizeof(MvarMVStruct *) * SubMatSize * SubMatSize);

	for (i = 0; i < MatSize; i++) {
	    MvarGetSubMVMat(MVsMatrix, MatSize, TmpMVsSubMat, 1, i + 1, FALSE);
	    TmpSubDetMV = MvarMVDeterminant((MvarMVStruct const * const *)
					            TmpMVsSubMat, SubMatSize);
	    TmpMVToAdd = MvarMVMult((const MvarMVStruct *) MVsMatrix[i],
				    TmpSubDetMV);

	    if (i == 0) {
		/* Refer to the computed MV and don't free it. Special case. */
		SumMV = TmpMVToAdd;
		MvarMVFree(TmpSubDetMV);
		Sign *= -1;
		continue;
	    }
	    else { /* i > 0 */
		/* Expand using the first row. Signs begin at '+' and        */
		/* alternate. The result is gradually constructed in SumMV.  */
		switch (Sign) {
		    case 1:
			TmpSumMV = MvarMVAdd(SumMV, TmpMVToAdd);
			MvarMVFree(SumMV);
			SumMV = TmpSumMV;
			break;
		    case -1:
			TmpSumMV = MvarMVSub(SumMV, TmpMVToAdd);
			MvarMVFree(SumMV);
			SumMV = TmpSumMV;
			break;
		    default:
			assert(0);
		}
		Sign *= -1;
		MvarMVFree(TmpMVToAdd);
		MvarMVFree(TmpSubDetMV);
	    }
	}

	IritFree((MvarMVStruct **) TmpMVsSubMat);        /* Pointers only. */
	return SumMV;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Extracts the (MatSize - 1) by (MatSize - 1) sub matrix of MVs, obtained  *
* by erasing one row and one column from the input MVs matrix of size        *
* MatSize by MatSize. Both are represented as a linear array of MVs.	     *
* The output can either refer to the original MVs, or to a new copy.	     *
*									     *
* PARAMETERS:                                                                *
*   MVsMatrix:	The original MVs matrix.			             *
*   SubMat:	The output, sub-matrix. Assumed to be allocated as required. *
*   MatSize:    The size of each of the dimensions of the square MVsMatrix.  *
*   Row:	The row index to erase.					     *
*   Col:	The col index to erase.					     *
*   CreateCopy:	Determines if the output refers to new copies of MVs objects *
*		or refers to the correct subset of the input objects.        *
*		(Set to FALSE when copies are not needed, for better         *
*		performance.						     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void.								     *
*****************************************************************************/
static void MvarGetSubMVMat(const MvarMVStruct * const *MVsMatrix, 
			    int MatSize,
			    const MvarMVStruct **SubMat,
			    int Row,
			    int Col,
			    CagdBType CreateCopy)
{
    int i, 
	CurRow, CurCol, 
	SubMatElems = 0;
    
    CurRow = 1;
    CurCol = 0;
    for (i = 0; i < MatSize * MatSize; i++) {
	/* Keep track of the double index format, and avoid the unneeded    */
	/* conversion at each iteration.				    */
	if (CurCol < MatSize)
	    CurCol++;
	else {			 /* Start a new row of the original Matrix. */
	    CurCol = 1;
	    CurRow++;
	}

	if (CurRow == Row || CurCol == Col)        /* An element for erase. */
	    continue;
	else {                   /* i is not an index of an erased element. */
	    SubMat[SubMatElems] = (MvarMVStruct *)
		       (CreateCopy ? MvarMVCopy(MVsMatrix[i]) : MVsMatrix[i]);
	    SubMatElems++;
	}
    }
    assert(SubMatElems == IRIT_SQR(MatSize - 1));
}
