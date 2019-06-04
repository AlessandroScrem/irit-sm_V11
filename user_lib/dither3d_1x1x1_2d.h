/******************************************************************************
* dither3d_1x1x1_2d.h - 1x1x1 dithering matrix for 2 images.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Gershon Elber,					 July 2010.   *
******************************************************************************/

static const char Dither2Imgs3DySize1Penalty[2][2] = {
    {  0,  -1 },
    {  1,   0 }
};

static const unsigned char Dither2Imgs3DSize1Matrices[2][2][1] = {
    {
	{ 0 }, /* (0 x 0 (0 OnBits)) */
	{ 0 }, /* (0 x 1 (0 OnBits)) */
    },
    {
	{ 0 }, /* (1 x 0 (0 OnBits)) */
	{ 1 }  /* (1 x 1 (1 OnBits)) */
    }
};
