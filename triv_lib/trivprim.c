/******************************************************************************
* trivprim.c - constructs non-singular 3D primitives.			      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Fady Massarwi and Gershon Elber,		 Dec 2014.	      *
******************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/cagd_lib.h"
#include "triv_loc.h"

static CagdBType TrivPrimCreateNonSingPlanarDiskSrfs(
					       const CagdVType Center,
					       CagdRType Radius,
					       CagdBType Rational,
					       CagdRType InternalQuadEdge,
					       CagdSrfStruct **InnerSrfs);

/*****************************************************************************
* DESCRIPTION:                                                               *
*   A surface constructor of a disk, centered at Center with radius Radius.  *
*   The disk is built from 5 non singular surfaces, an internal cube, and a  *
*   ruled surfaces between the disk outer circle patches and the cube edges. *
*                                                                            *
* PARAMETERS:                                                                *
*   Center:    Of the disk.                                                  *
*   Radius:    Of the disk.                                                  *
*   Rational:  TRUE for precise rational, FALSE for polynomial approximation.*
*   InternalQuadSize: size of the internal quad surface edge.                *
*   InnerSrfs: Output parameter, contains a vector of 5 surfaces, allocated  *
*              dynamocally. 					             *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdBType:  TRUE on success, FALSE otherwise.	                     *
*****************************************************************************/
static CagdBType TrivPrimCreateNonSingPlanarDiskSrfs(
					       const CagdVType Center,
					       CagdRType Radius,
 					       CagdBType Rational,
                                               CagdRType InternalQuadSize,
					       CagdSrfStruct **InnerSrfs)
{
    int i;
    CagdRType 
	InnerQuadHalfEdge = InternalQuadSize / 2.0;	
    CagdCrvStruct *BoxLines[4], *CircleArcs[4], *UnitCircle, *QuarterCrv;
    CagdPtStruct Pt1, Pt2;
    CagdSrfStruct *TmpSrf;  
    IrtHmgnMatType RotMat45, RotMat90;

    *InnerSrfs = (CagdSrfStruct *) IritMalloc(sizeof(CagdSrfStruct *) * 5);

    if (Rational) {
	UnitCircle = BspCrvCreateUnitCircle();
    } 
    else {
	UnitCircle = BspCrvCreateUnitPCircle();
    }
    
    QuarterCrv = CagdCrvRegionFromCrv(UnitCircle, 0.0, 1.0);
    CagdCrvFree(UnitCircle);
    MatGenMatRotZ1(IRIT_DEG2RAD(45), RotMat45);
    MatGenMatRotZ1(IRIT_DEG2RAD(90), RotMat90);
    CircleArcs[0] = CagdCrvMatTransform(QuarterCrv, RotMat45);
    CagdCrvTransform(CircleArcs[0], NULL, Radius);
    CagdCrvFree(QuarterCrv);
    
    for (i = 1; i < 4; ++i) {
        CircleArcs[i] = CagdCrvMatTransform(CircleArcs[i - 1], RotMat90);
    }

    for (i = 0; i < 4; ++i) {
	CagdCrvTransform(CircleArcs[i], Center, 1.0);
    }
    
    Pt1.Pt[0] = Center[0] + InnerQuadHalfEdge;
    Pt1.Pt[1] = Center[1] + InnerQuadHalfEdge;
    Pt1.Pt[2] = Center[2];
    
    Pt2.Pt[0] = Center[0] - InnerQuadHalfEdge;
    Pt2.Pt[1] = Center[1] + InnerQuadHalfEdge;
    Pt2.Pt[2] = Center[2];

    BoxLines[0] = CagdMergePtPt(&Pt2, &Pt1);
    
    Pt1.Pt[0] = Center[0] - InnerQuadHalfEdge;
    Pt1.Pt[1] = Center[1] - InnerQuadHalfEdge;
    Pt1.Pt[2] = Center[2];
    
    BoxLines[1] = CagdMergePtPt(&Pt2, &Pt1);

    Pt2.Pt[0] = Center[0] + InnerQuadHalfEdge;
    Pt2.Pt[1] = Center[1] - InnerQuadHalfEdge;
    Pt2.Pt[2] = Center[2];

    BoxLines[2] = CagdMergePtPt(&Pt1, &Pt2);

    Pt1.Pt[0] = Center[0] + InnerQuadHalfEdge;
    Pt1.Pt[1] = Center[1] + InnerQuadHalfEdge;
    Pt1.Pt[2] = Center[2];
    
    BoxLines[3] = CagdMergePtPt(&Pt1, &Pt2);
    
    TmpSrf = CagdBoolSumSrf(BoxLines[1], BoxLines[3], 
	                    BoxLines[0], BoxLines[2]);
    InnerSrfs[4] = CagdCoerceSrfsTo(TmpSrf, CAGD_PT_E3_TYPE, FALSE);
    CagdSrfFree(TmpSrf);
    
    /* before the ruling - invert curves to be the same direction with */
    /* the box lines.		       				       */
    QuarterCrv = CircleArcs[0];
    CircleArcs[0] = CagdCrvReverse(CircleArcs[0]);
    CagdCrvFree(QuarterCrv);

    QuarterCrv = CircleArcs[3];
    CircleArcs[3] = CagdCrvReverse(CircleArcs[3]);
    CagdCrvFree(QuarterCrv);
   
    for (i = 0; i < 4; ++i) {
        InnerSrfs[i] = CagdRuledSrf(BoxLines[i], CircleArcs[i], 2, 2);
    }

    for (i = 0; i < 5; ++i) {
	if (CAGD_IS_BEZIER_SRF(InnerSrfs[i])) {
	    CagdSrfStruct
	        *BSSrf = CagdCnvrtBzr2BspSrf(InnerSrfs[i]);

	    CagdSrfFree(InnerSrfs[i]);
	    InnerSrfs[i] = BSSrf;
	}
    }

    for (i = 0; i < 4; ++i) {
	CagdCrvFree(BoxLines[i]);
	CagdCrvFree(CircleArcs[i]);
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A trivariate constructor of a box, parallel to main axes.		     M
*                                                                            *
* PARAMETERS:                                                                M
*   MinX, MinY, MinZ:   Minimum range of box trivariate.                     M
*   MaxX, MaxY, MaxZ:   Maximum range of box trivariate.                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  A trivariate of a box model.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivNSPrimCone, TrivNSPrimCone2, TrivNSPrimCylinder, TrivNSPrimSphere,   M
*   TrivNSPrimTorus							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivNSPrimBox                                                            M
*****************************************************************************/
TrivTVStruct *TrivNSPrimBox(CagdRType MinX,
			    CagdRType MinY,
			    CagdRType MinZ,
			    CagdRType MaxX,
			    CagdRType MaxY,
			    CagdRType MaxZ)
{
    CagdVecStruct Vec;
    CagdSrfStruct
        *SrfXYMin = CagdPrimPlaneSrf(MaxX, MinY, MinX, MaxY, MinZ);
    TrivTVStruct *TV;

    Vec.Vec[0] = 0;
    Vec.Vec[1] = 0;
    Vec.Vec[2] = MaxZ - MinZ;
    
    TV = TrivExtrudeTV(SrfXYMin, &Vec);
        
    CagdSrfFree(SrfXYMin); 
   
    return TV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A trivariate constructor of a cylinder, centered at Center with radius   M
* Radius, and a Height height.						     M
*   The cylinder is built from 5 non singular ruled trivariates.             M
*                                                                            *
* PARAMETERS:                                                                M
*   Center:   Of the cylinder.                                               M
*   Radius:   Of the cylinder.                                               M
*   Height:   Of the cylinder.						     M
*   Rational: TRUE for precise rational, FALSE for polynomial approximation. M
*   InternalCubeSize: Size of the internal cube trivariate edge.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  A list of 5 trivariates, representing the cylinder.     M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivNSPrimBox, TrivNSPrimCone, TrivNSPrimCone2, TrivNSPrimSphere,        M
*   TrivNSPrimTorus							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivNSPrimCylinder                                                       M
*****************************************************************************/
TrivTVStruct *TrivNSPrimCylinder(const CagdVType Center,
				 CagdRType Radius,
				 CagdRType Height,
				 CagdBType Rational,
				 CagdRType InternalCubeSize)
{
    int i;    
    CagdSrfStruct *InnerSrfs[5];
    TrivTVStruct *CylParts[5];
    CagdVecStruct Vec;
    
    if (!TrivPrimCreateNonSingPlanarDiskSrfs(Center, Radius, Rational, 
					     InternalCubeSize,
					     InnerSrfs)) {
	return NULL;
    }

    Vec.Vec[0] = 0;
    Vec.Vec[1] = 0;
    Vec.Vec[2] = Height;

    for (i = 0; i < 5; ++i) {
        CylParts[i] = TrivExtrudeTV(InnerSrfs[i], &Vec);
	CagdSrfFree(InnerSrfs[i]);
    }
    
    for (i = 0; i < 4; ++i) {
        CylParts[i] -> Pnext = CylParts[i + 1];
    }
    CylParts[4] -> Pnext = NULL;
        
    return CylParts[0];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A trivariate constructor of a cone, centered at Center with radius       M
* Radius, and a Height height.						     M
*   The cone is built from 5 ruled trivariates. Note the cone will be        M
* (only) at the apex.						             M
*                                                                            *
* PARAMETERS:                                                                M
*   Center:   Of the cone.                                                   M
*   Radius:   Of the cone.                                                   M
*   Height:   Of the cone.						     M
*   Rational: TRUE for precise rational, FALSE for polynomial approximation. M
*   InternalCubeSize: Size of the internal cube trivariate edge.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  A list of 5 trivariates, representing the cone.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivNSPrimBox, TrivNSPriCmone2, TrivNSPrimCylinder, TrivNSPrimSphere,    M
*   TrivNSPrimTorus							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivNSPrimCone                                                           M
*****************************************************************************/
TrivTVStruct *TrivNSPrimCone(const CagdVType Center,
			     CagdRType Radius,
			     CagdRType Height,
			     CagdBType Rational,
			     CagdRType InternalCubeSize)
{
    return TrivNSPrimCone2(Center, Radius, IRIT_UEPS, Height, Rational,
			   IRIT_UEPS / 2.0);
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A trivariate constructor of a cone, centered at Center with major        M
* radius MajorRadius and minor radius MinorRadius, and a Height height.      M
*   The cone is built from 5 non singular ruled trivariates.                 M
*                                                                            *
* PARAMETERS:                                                                M
*   Center:       Of the cone.                                               M
*   MajorRadius:  Of the cone.                                               M
*   MinorRadius:  Of the cone.                                               M
*   Height:       Of the cone.						     M
*   Rational: TRUE for precise rational, FALSE for polynomial approximation. M
*   InternalCubeSize: Size of the internal cube trivariate edge.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  A list of 5 trivariates, representing the cone.         M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivNSPrimBox, TrivNSPrimCone, TrivNSPrimCylinder, TrivNSPrimSphere,     M
*   TrivNSPrimTorus							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivNSPrimCone2                                                          M
*****************************************************************************/
TrivTVStruct *TrivNSPrimCone2(const CagdVType Center,
			      CagdRType MajorRadius,
			      CagdRType MinorRadius,
			      CagdRType Height,
			      CagdBType Rational,
			      CagdRType InternalCubeSize)
{
    int i;    
    CagdSrfStruct *InnerSrfsMajor[5], *InnerSrfsMinor[5];
    TrivTVStruct *CylParts[5];
    CagdRType Trans[3];

    Trans[0] = 0.0;
    Trans[1] = 0.0;
    Trans[2] = Height;    

    if (!TrivPrimCreateNonSingPlanarDiskSrfs(Center, MajorRadius, Rational, 
					     InternalCubeSize,
					     InnerSrfsMajor)) {
	return NULL;
    }

    if (!TrivPrimCreateNonSingPlanarDiskSrfs(Center, MinorRadius, Rational, 
					     InternalCubeSize,
					     InnerSrfsMinor)) {
	return NULL;
    }

    for (i = 0; i < 5; ++i) {
	CagdSrfTransform(InnerSrfsMinor[i], Trans, 1.0);
	CylParts[i] = TrivRuledTV(InnerSrfsMajor[i], InnerSrfsMinor[i], 2 , 2);
	CagdSrfFree(InnerSrfsMajor[i]);
	CagdSrfFree(InnerSrfsMinor[i]);
    }    

    for (i = 0; i < 4; ++i) {
	CylParts[i] -> Pnext = CylParts[i + 1];
    }
    CylParts[4] -> Pnext = NULL;

    return CylParts[0];
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A trivariate constructor of a ball, centered at Center and radius        M
*   Radius, built from 7 non singular trivariates, an inner cube             M
*   and 6 ruled trivariate between patches of the ball's outer sphere        M
*   surface and the cube faces.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   Center:   Of the ball.	                                             M
*   Radius:   Of the ball.		                                     M
*   Rational: TRUE for precise rational, FALSE for polynomial approximation. M
*   InternalCubeSize: Size of the internal cube trivariate edge.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  A list of 7 trivariates, representing the sphere.       M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivNSPrimBox, TrivNSPrimCone, TrivNSPrimCone2, TrivNSPrimCylinder,      M
*   TrivNSPrimTorus							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivNSPrimSphere                                                         M
*****************************************************************************/
TrivTVStruct *TrivNSPrimSphere(const CagdVType Center,
			       CagdRType Radius,
			       CagdBType Rational,
			       CagdRType InternalCubeSize)
{
    /*    UMin, UMax, VMin, VMax, WMin, WMax   */
    /*      0     1     2     3     4     5    */
    int i,
        SphereCubeBndryMatching[6] = { 4, 3, 5, 2, 1, 0 };
    CagdRType NegCenter[3],
        InternalCubeHalfEdge = InternalCubeSize/2;
    /* x surfaces 90, 180, 270, 360 */
    CagdSrfStruct *InnerCubeBoundaries[6], *TmpSrf,
        *SphereSurfaces =  CagdPrimCubeSphereSrf(Center, Radius, Rational);
    IrtHmgnMatType RotMat;
    TrivTVStruct *RuledTVs[6],    
        *InternalCubeTV = TrivNSPrimBox(Center[0] - InternalCubeHalfEdge,
					Center[1] - InternalCubeHalfEdge,
					Center[2] - InternalCubeHalfEdge,
					Center[0] + InternalCubeHalfEdge,
					Center[1] + InternalCubeHalfEdge,
					Center[2] + InternalCubeHalfEdge);

    NegCenter[0] = -Center[0];
    NegCenter[1] = -Center[1];
    NegCenter[2] = -Center[2];

    IRIT_GEN_COPY(InnerCubeBoundaries, TrivBndrySrfsFromTV(InternalCubeTV),
                  sizeof(CagdSrfStruct *) * 6);
        
    for (i = 0; i < 6; ++i) {
        CagdSrfTransform(InnerCubeBoundaries[i], NegCenter, 1.0);
    }

    MatGenMatRotZ1(IRIT_DEG2RAD(180), RotMat);
    InnerCubeBoundaries[SphereCubeBndryMatching[0]] =
        CagdSrfMatTransform(InnerCubeBoundaries[SphereCubeBndryMatching[0]],
                            RotMat);

    MatGenMatRotY1(IRIT_DEG2RAD(180), RotMat);
    InnerCubeBoundaries[SphereCubeBndryMatching[1]] =
	CagdSrfMatTransform(InnerCubeBoundaries[SphereCubeBndryMatching[1]],
                            RotMat);

    MatGenMatRotX1(IRIT_DEG2RAD(90), RotMat);
    InnerCubeBoundaries[SphereCubeBndryMatching[4]] =
	CagdSrfMatTransform(InnerCubeBoundaries[SphereCubeBndryMatching[4]],
	                    RotMat);

    InnerCubeBoundaries[SphereCubeBndryMatching[5]] =
	CagdSrfMatTransform(InnerCubeBoundaries[SphereCubeBndryMatching[5]],
	                    RotMat);
    
    for (i = 0; i < 6; ++i) {
        CagdSrfTransform(InnerCubeBoundaries[i], Center, 1.0);
    }
    
    TmpSrf = SphereSurfaces;
    for (i = 0; i < 6; ++i) {
        RuledTVs[i] = 
	    TrivRuledTV(InnerCubeBoundaries[SphereCubeBndryMatching[i]],
                        TmpSrf, 1 , 2);
        TmpSrf = TmpSrf -> Pnext;
    }
    
    for (i = 0; i < 5; ++i) {
        RuledTVs[i] -> Pnext = RuledTVs[i+1];
    }
    RuledTVs[5] -> Pnext = NULL;
    
    InternalCubeTV -> Pnext = RuledTVs[0];
    return InternalCubeTV;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   A trivariate constructor of a torus, centered at Center with major       M
*   radius MajorRadius and minor radius MinorRadius. The torous is built     M
*   from 5 non singular trivariates of revolution.                           M
*                                                                            *
* PARAMETERS:                                                                M
*   Center:       Of the torus.                                              M
*   MajorRadius:  Of the torus.                                              M
*   MinorRadius:  Of the torus.                                              M
*   Rational: TRUE for precise rational, FALSE for polynomial approximation. M
*   InternalCubeSize: Size of the internal cube trivariate edge.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   TrivTVStruct *:  A list of 5 trivariates, representing the torus.        M
*                                                                            *
* SEE ALSO:                                                                  M
*   TrivNSPrimBox, TrivNSPrimCone, TrivNSPrimCone2, TrivNSPrimCylinder,      M
*   TrivNSPrimSphere							     M
*                                                                            *
* KEYWORDS:                                                                  M
*   TrivNSPrimTorus                                                          M
*****************************************************************************/
TrivTVStruct *TrivNSPrimTorus(const CagdVType Center,
			      CagdRType MajorRadius,
			      CagdRType MinorRadius,
			      CagdBType Rational,
			      CagdRType InternalCubeSize)
{
    int i;    
    CagdSrfStruct *InnerSrfs[5];
    TrivTVStruct *CylParts[5];    
    const CagdVType 
	Origin = { 0.0, 0.0, 0.0 };
    IrtHmgnMatType XRotMat, XTransMat, TranMat;
    
    MatGenMatRotX1(IRIT_DEG2RAD(90.0), XRotMat);
    MatGenMatTrans(MajorRadius, 0.0, 0.0, XTransMat);
    MatMultTwo4by4(TranMat, XTransMat, XRotMat);

    if (!TrivPrimCreateNonSingPlanarDiskSrfs(Origin, MinorRadius, Rational, 
					     InternalCubeSize, InnerSrfs)) {
	return NULL;
    }

    for (i = 0; i < 5; ++i) {
	CagdSrfStruct
	    *TSrf = CagdSrfMatTransform(InnerSrfs[i], TranMat);

	CagdSrfFree(InnerSrfs[i]);
	CylParts[i] = TrivTVOfRev(TSrf);
	TrivTVTransform(CylParts[i], (CagdRType *) Center, 1.0);
	CagdSrfFree(TSrf);
    }

    for (i = 0; i < 4; ++i) {
	CylParts[i] -> Pnext = CylParts[i + 1];
    }
    CylParts[4] -> Pnext = NULL;

    return CylParts[0];
}
