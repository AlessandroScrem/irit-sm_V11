/*****************************************************************************
* Reads one image in.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber			       Ver 1.0,	Dec. 1998    *
*****************************************************************************/

#include <math.h>
#include "inc_irit/irit_sm.h"
#include "misc_loc.h"

#ifdef IRIT_HAVE_URT_RLE
#define NO_DECLARE_MALLOC /* For rle.h */
#include <rle.h>
#include <rle_raw.h>
#endif /* IRIT_HAVE_URT_RLE */

#ifdef IRIT_HAVE_GIF_LIB
#include "gif_lib.h"
#endif /* IRIT_HAVE_GIF_LIB */

#ifdef IRIT_HAVE_PNG_LIB
#include "png.h"
#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif
#endif /* IRIT_HAVE_PNG_LIB */

#ifdef IRIT_HAVE_JPG_LIB
#include <jpeglib.h>
#include <jerror.h>
#endif/* IRIT_HAVE_JPG_LIB */

#ifdef __WINNT__
#include <io.h>
#include <fcntl.h>
#endif /* __WINNT__ */

typedef struct LoadedImagesStruct {
    struct LoadedImagesStruct *Pnext;
    char *FileName;
    int MaxX, MaxY, Alpha;
    IrtBType *Data;
} LoadedImagesStruct;

IRIT_STATIC_DATA unsigned int
    GlblImageXAlignment = 0xffffffff; /* No alignment by default. */

IRIT_STATIC_DATA LoadedImagesStruct
    *GlblLoadedImagesList = NULL;

static IrtBType *PPMReadImage(const char *File,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha);
static IrtBType *RLEReadImage(const char *File,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha);
static IrtBType *GIFReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha);
static IrtBType *PNGReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha);
static IrtBType *JPEGReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Sets an alignment of the width of the image. For example, OGL requires   M
* images to have width aligned on 4-bytes words (alignment 4).               M
*                                                                            *
* PARAMETERS:                                                                M
*   Alignment:   Word size alignment required.                               M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:     old alignment.                                                  M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgReadImage                                                          M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgReadImageXAlign                                                    M
*****************************************************************************/
int IrtImgReadImageXAlign(int Alignment)
{
    int OldVal = GlblImageXAlignment;

    GlblImageXAlignment = ~(Alignment - 1);

    return OldVal;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reads one image in from a file named ImageFileName.  The image is        M
* returned as a vector of RGBRGB... of size (MaxX+1) * (MaxY+1) * 3.         M
*                                                                            *
* PARAMETERS:                                                                M
*   ImageFileName:   Name of image to read.                                  M
*   MaxX:            Maximum X of read image is saved here.                  M
*   MaxY:            Maximum Y of read image is saved here.                  M
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  M
*		     and will return TRUE if successful in loading it.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct *:  A vector of RGBRGB... of size 		     M
*		(MaxX+1) * (MaxY+1) * 3 or NULL if failed.		     M
*		  If however, Alpha is requested and found RGBARGBA... is    M
*		returned as IrtImgRGBAPxlStruct.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgReadImage2, IrtImgReadImage3, IrtImgReadImageXAlign                M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgReadImage                                                          M
*****************************************************************************/
IrtImgPixelStruct *IrtImgReadImage(const char *ImageFileName,
				   int *MaxX,
				   int *MaxY,
				   int *Alpha)
{
    const char *Type;

    if (ImageFileName == NULL) {
        IRIT_FATAL_ERROR("Empty image file name to write to.");
        return NULL;
    }

    if ((Type = strrchr(ImageFileName, '.')) == NULL)
        Type = "";

    if (stricmp(Type, ".Z") == 0) {
        Type--;
	while (Type != ImageFileName && *Type != '.')
	    Type--;
    }

    if (stricmp(Type, ".ppm") == 0) {
        return (IrtImgPixelStruct *) PPMReadImage(ImageFileName, MaxX, MaxY,
						  Alpha);
    }
    else if (stricmp(Type, ".rle") == 0 || stricmp(Type, ".rle.Z") == 0) {
	return (IrtImgPixelStruct *) RLEReadImage(ImageFileName, MaxX, MaxY,
						  Alpha);
    }
    else if (stricmp(Type, ".gif") == 0) {
	return (IrtImgPixelStruct *) GIFReadImage(ImageFileName, MaxX, MaxY,
						  Alpha);
    }
    else if (stricmp(Type, ".png") == 0) {
	return (IrtImgPixelStruct *) PNGReadImage(ImageFileName, MaxX, MaxY,
						  Alpha);
    }
    else if ((stricmp(Type, ".jpg") == 0) || (stricmp(Type, ".jpeg") == 0)) {
	return (IrtImgPixelStruct *) JPEGReadImage(ImageFileName, MaxX, MaxY,
						  Alpha);
    }
    else {
	IRIT_WARNING_MSG_PRINTF(
	    IRIT_EXP_STR("Texture file \"%s\" must be image of type 'rle', 'ppm', 'gif', 'jpeg' or 'png'\n"),
	    ImageFileName);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Same as IrtImgReadImage2 but if a name of an image repeats itself, the   M
* image is read only once.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   ImageFileName:   Name of image to read.                                  M
*   MaxX:            Maximum X of read image is saved here.                  M
*   MaxY:            Maximum Y of read image is saved here.                  M
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  M
*		     and will return TRUE if successful in loading it.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct *:  A vector of RGBRGB... of size                      M
*		(MaxX+1) * (MaxY+1) * 3 or NULL if failed.		     M
*		  If however, Alpha is requested and found RGBARGBA... is    M
*		returned as IrtImgRGBAPxlStruct.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgReadImage, IrtImgReadImage3, IrtImgReadClrCache                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgReadImage2                                                         M
*****************************************************************************/
IrtImgPixelStruct *IrtImgReadImage2(const char *ImageFileName,
				    int *MaxX,
				    int *MaxY,
				    int *Alpha)
{
    LoadedImagesStruct *LoadedImage;
    IrtImgPixelStruct *Data;

    /* Search if we already loaded this image. */
    for (LoadedImage = GlblLoadedImagesList;
	 LoadedImage != NULL;
	 LoadedImage = LoadedImage -> Pnext) {
	if (strcmp(ImageFileName, LoadedImage -> FileName) == 0) {
	    *MaxX = LoadedImage -> MaxX;
	    *MaxY = LoadedImage -> MaxY;
	    *Alpha = LoadedImage -> Alpha;
	    return (IrtImgPixelStruct *) LoadedImage -> Data;
	}
    }

    if ((Data = IrtImgReadImage(ImageFileName, MaxX, MaxY, Alpha)) != NULL) {
	/* Add it to global list of loaded images. */
	LoadedImage = (LoadedImagesStruct *)
				    IritMalloc(sizeof(LoadedImagesStruct));
	LoadedImage -> FileName = IritStrdup(ImageFileName);
	LoadedImage -> MaxX = *MaxX;
	LoadedImage -> MaxY = *MaxY;
	LoadedImage -> Alpha = *Alpha;
	LoadedImage -> Data = (IrtBType *) Data;
	LoadedImage -> Pnext = GlblLoadedImagesList;
	GlblLoadedImagesList = LoadedImage;
    }

    return Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Same as IrtImgReadImage2 but if a name of an image repeats itself, the   M
* new image replaces the old one.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   ImageFileName:   Name of image to read.                                  M
*   MaxX:            Maximum X of read image is saved here.                  M
*   MaxY:            Maximum Y of read image is saved here.                  M
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  M
*		     and will return TRUE if successful in loading it.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct *:  A vector of RGBRGB... of size                      M
*		(MaxX+1) * (MaxY+1) * 3 or NULL if failed.		     M
*		  If however, Alpha is requested and found RGBARGBA... is    M
*		returned as IrtImgRGBAPxlStruct.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgReadImage, IrtImgReadImage2, IrtImgReadClrCache                    M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgReadImage3                                                         M
*****************************************************************************/
IrtImgPixelStruct *IrtImgReadImage3(const char *ImageFileName,
				    int *MaxX,
				    int *MaxY,
				    int *Alpha)
{
    LoadedImagesStruct *LoadedImage;
    IrtImgPixelStruct *Data;

    if ((Data = IrtImgReadImage(ImageFileName, MaxX, MaxY, Alpha)) != NULL) {
        /* Search if we already loaded this image. */
        for (LoadedImage = GlblLoadedImagesList;
	     LoadedImage != NULL;
	     LoadedImage = LoadedImage -> Pnext) {
	    if (strcmp(ImageFileName, LoadedImage -> FileName) == 0) {
	        IritFree(LoadedImage -> FileName);
		IritFree(LoadedImage -> Data);
		break;
	    }
	}

	if (LoadedImage == NULL) {                 /* Is a new image name. */
	    LoadedImage = (LoadedImagesStruct *)
				    IritMalloc(sizeof(LoadedImagesStruct));
	    LoadedImage -> Pnext = GlblLoadedImagesList;
	    GlblLoadedImagesList = LoadedImage;
	}

	LoadedImage -> FileName = IritStrdup(ImageFileName);
	LoadedImage -> MaxX = *MaxX;
	LoadedImage -> MaxY = *MaxY;
	LoadedImage -> Alpha = *Alpha;
	LoadedImage -> Data = (IrtBType *) Data;
    }

    return Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Clears the cache of read images.                                         M
*                                                                            *
* PARAMETERS:                                                                M
*   None                                                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgReadImage2, IrtImgReadClrOneImage                                  M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgReadClrCache                                                       M
*****************************************************************************/
void IrtImgReadClrCache(void)
{
    while (GlblLoadedImagesList != NULL) {
	LoadedImagesStruct
	    *LoadedImage = GlblLoadedImagesList;

	GlblLoadedImagesList = GlblLoadedImagesList -> Pnext;
	IritFree(LoadedImage -> FileName);
	IritFree(LoadedImage -> Data);
	IritFree(LoadedImage);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Removes one image, by name from the cache of images.                     M
*                                                                            *
* PARAMETERS:                                                                M
*   ImageName:  Name of image to remove.                                     M
*                                                                            *
* RETURN VALUE:                                                              M
*   void                                                                     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtImgReadImage2, IrtImgReadClrCache                                     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgReadClrOneImage                                                    M
*****************************************************************************/
void IrtImgReadClrOneImage(const char *ImageName)
{
    LoadedImagesStruct *LI;
 
    if (GlblLoadedImagesList == NULL)
	return;

    /* First image in list? */
    if (stricmp(ImageName, GlblLoadedImagesList -> FileName) == 0) {
	LI = GlblLoadedImagesList;
	GlblLoadedImagesList = GlblLoadedImagesList -> Pnext;
	IritFree(LI -> FileName);
	IritFree(LI -> Data);
	IritFree(LI);
	return;
    }
    
    for (LI = GlblLoadedImagesList; LI -> Pnext != NULL; LI = LI -> Pnext) {
	if (stricmp(ImageName, LI -> Pnext -> FileName) == 0) {
	    LoadedImagesStruct
		*LI2 = LI -> Pnext;

	    LI -> Pnext = LI2 -> Pnext;
	    IritFree(LI2 -> FileName);
	    IritFree(LI2 -> Data);
	    IritFree(LI2);
	    return;
	}
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reads the dimensions of one image in from a file named ImageFileName.    M
*                                                                            *
* PARAMETERS:                                                                M
*   ImageFileName:   Name of image to read.                                  M
*   Width:           The width of the image.				     M
*   Height:          The height of the image.				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:	     TRUE, if loaded successfully. Otherwise FALSE.	     M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgGetImageSize							     M
*****************************************************************************/
int IrtImgGetImageSize(const char *ImageFileName,
		       int *Width,
		       int *Height)
{
    int Alpha, MaxX, MaxY;
    IrtImgPixelStruct 
	*Image = IrtImgReadImage(ImageFileName, &MaxX, &MaxY, &Alpha);

    if (Image != NULL) {
        *Width = MaxX + 1;
	*Height = MaxY + 1;
	IritFree(Image);
	return TRUE;
    }
    
    IRIT_WARNING_MSG_PRINTF(
	IRIT_EXP_STR("IrtImgReadImageSize: Image file \"%s\" could not be fetched.\n"),
	ImageFileName);
    return FALSE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Verifies the alignment of the image.  Returned is either the input Data  M
* if aligned, or a new aligned copy (and Data is freed).                     M
*                                                                            *
* PARAMETERS:                                                                M
*   Data:    Image data to verify alignment, in place.                       M
*   MaxX:    Current dimension of image.  Might be changed after alignment.  M
*   MaxY:    Current dimension of image.				     M
*   Alpha:   Do we have Alpha in the image?				     M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtBType *: The verified image data.				     M
*                                                                            *
* KEYWORDS:                                                                  M
*   _IrtImgVerifyAlignment                                                   M
*****************************************************************************/
IrtBType *_IrtImgVerifyAlignment(IrtBType *Data,
				 int *MaxX,
				 int *MaxY,
				 int Alpha)
{
    int x, y,
	Width = 1 + *MaxX,
	Height = 1 + *MaxY;

    if ((Width & GlblImageXAlignment) != Width) {     /* We have alignments. */
	IrtBType *p, *q,
	    *AlignedData = IritMalloc((Alpha ? 4 : 3) * Width * Height);

	*MaxX = (Width & GlblImageXAlignment) - 1;              /* Align it. */

        for (y = 0, p = Data, q = AlignedData; y < Height; y++) {
	    for (x = 0; x < Width; x++) {
	        if (x <= *MaxX) {
		    /* Copy pixel. */
		    *q++ = *p++;
		    *q++ = *p++;
		    *q++ = *p++;
		    if (Alpha)
		        *q++ = *p++;
		}
		else
		    p += Alpha ? 4 : 3;                  /* Skip end pixels. */
	    }
	}

	IritFree(Data);
	Data = AlignedData;
    }

    return Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads image file in PPM format.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   ImageFileName:   Name of PPM image to read.                              *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType *:  Pointer to dynamicaly created image, NULL if non.           *
*****************************************************************************/
static IrtBType *PPMReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha)
{
    int x, y, Width, Height;
    char Line[IRIT_LINE_LEN_LONG], Line2[IRIT_LINE_LEN_LONG];
    IrtBType *Data;
    FILE *PPMLoad;

    *Alpha = FALSE;				        /* No alpha in PPM. */

#if defined(__WINNT__) || defined(__WINCE__)
    if ((PPMLoad = fopen(ImageFileName, "rb")) == NULL) {
#else
    if ((PPMLoad = fopen(ImageFileName, "r")) == NULL) {
#endif /* __WINNT__ || __WINCE__ */
        IRIT_WARNING_MSG_PRINTF("Failed to read PPM file \"%s\"\n",
				ImageFileName);
        return NULL;
    }

    fgets(Line2, IRIT_LINE_LEN_LONG - 1, PPMLoad);
    if (strncmp(Line2, "P3", 2) != 0 && strncmp(Line2, "P6", 2) != 0) {
        IRIT_WARNING_MSG_PRINTF("P3 or P6 expected, found \"%s\"\n", Line2);
        return NULL;
    }

    fgets(Line, IRIT_LINE_LEN_LONG - 1, PPMLoad);
    while (Line[0] == '#')
        fgets(Line, IRIT_LINE_LEN_LONG - 1, PPMLoad);
    sscanf(Line, "%d %d", &Width, &Height);
    if (Width < 0 || Width > 100000 || Height < 0 || Height > 100000) {
        IRIT_WARNING_MSG_PRINTF("Unrealistic image size %d by %d\n",
				Width, Height);
        return NULL;
    }
    /* Get the "255" line. */
    fgets(Line, IRIT_LINE_LEN_LONG - 1, PPMLoad);

    *MaxX = Width - 1;
    *MaxY = Height - 1;

    /* Allocate the image. */
    Data = IritMalloc(3 * Width * Height);

    if (strncmp(Line2, "P6", 2) == 0) {
	int LineSize = Width * 3;

	fread(Data, 3 * Width * Height, 1, PPMLoad);

	/* Swap the lines so YMin is YMax. */
	for (y = 0; y <= (Height >> 1); y++) {
	    IrtBType
	        *p1 = &Data[(*MaxY - y) * LineSize],
	        *p2 = &Data[y * LineSize];

	    for (x = LineSize; x > 0; x--, p1++, p2++) {
	        IRIT_SWAP(IrtBType, *p1, *p2);
	    }
	}
    }
    else { /* P3 */
	int LineSize = Width * 3;

	assert(strncmp(Line2, "P3", 2) == 0);

        for (y = 0; y < Height; y++) {
	    IrtBType
	        *p = &Data[(*MaxY - y) * LineSize];

	    for (x = 0; x < Width; x++) {
	        int r, g, b;

	        fscanf(PPMLoad, "%d %d %d", &r, &g, &b);
		*p++ = (IrtBType) r;
		*p++ = (IrtBType) g;
		*p++ = (IrtBType) b;
	    }
	}
    }

    fclose(PPMLoad);

    Data = _IrtImgVerifyAlignment(Data, MaxX, MaxY, *Alpha);

    return Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads image file in RLE format.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   ImageFileName:   Name of RLE image to read.                              *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType *:  Pointer to dynamicaly created image, NULL if non.           *
*****************************************************************************/
static IrtBType *RLEReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha)
{
#ifdef IRIT_HAVE_URT_RLE
    rle_hdr Header;
    rle_pixel **Rows;
    int Error, x, y;
    IrtBType *Data, *p;

    Header = *rle_hdr_init(NULL);
    Header.rle_file = rle_open_f_noexit("RleLoadImage",
					(char *) ImageFileName, "r");
    if (!Header.rle_file) {
        IRIT_WARNING_MSG_PRINTF("Failed to read RLE file \"%s\"\n",
				ImageFileName);
        return NULL;
    }

    if (Error = rle_get_setup(&Header)) {
        rle_get_error(Error, "RleLoadImage", (char *) ImageFileName);
        return NULL;
    }
    rle_row_alloc(&Header, &Rows);
    *MaxX = Header.xmax - Header.xmin;
    *MaxY = Header.ymax - Header.ymin;

    /* Get alpha only if requested to get it. */
    if (*Alpha)
        *Alpha = Header.alpha;

    if (*Alpha)
        Data = p = IritMalloc(4 * (*MaxX + 1) * (*MaxY + 1));
    else
        Data = p = IritMalloc(3 * (*MaxX + 1) * (*MaxY + 1));

    for (y = 0; y <= *MaxY; y++) {
        rle_getrow(&Header, Rows);
        for (x = 0; x <= *MaxX; x++) {
            *p++ = Rows[RLE_RED][x];
            *p++ = Rows[RLE_GREEN][x];
            *p++ = Rows[RLE_BLUE][x];
	    if (*Alpha)
	        *p++ = Rows[RLE_ALPHA][x];
        }
    }

    rle_close_f(Header.rle_file);

    Data = _IrtImgVerifyAlignment(Data, MaxX, MaxY, *Alpha);

    return Data;
#else
    IRIT_WARNING_MSG("Utah raster tool kit is not supported\n");
    return NULL;			   /* Silent the compiler's warning. */
#endif /* IRIT_HAVE_URT_RLE */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads image file in GIF format.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   ImageFileName:   Name of GIF image to read.                              *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType *:  Pointer to dynamicaly created image, NULL if non.           *
*****************************************************************************/
static IrtBType *GIFReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha)
{
#ifdef IRIT_HAVE_GIF_LIB
    int TranspColorIndex = -1;
    GifFileType *GifFileIn;
    IrtBType *Data, *p, *Line;
    GifRecordType RecordType;

    if ((GifFileIn = DGifOpenFileName(ImageFileName)) == NULL) {
        IRIT_WARNING_MSG_PRINTF("Failed to read GIF file \"%s\"\n",
				ImageFileName);
        return NULL;
    }

    /* Scan the content of the GIF file and load the image(s) in: */
    do {
	int i, j, l, ExtCode;
	GifByteType *Extension;
	ColorMapObject *ColorMap;
	GifColorType *ColorMapEntry;

	if (DGifGetRecordType(GifFileIn, &RecordType) == GIF_ERROR)
	    return NULL;

	switch (RecordType) {
	    case IMAGE_DESC_RECORD_TYPE:
		if (DGifGetImageDesc(GifFileIn) == GIF_ERROR)
		    return NULL;

		ColorMap =
		    (GifFileIn -> Image.ColorMap ? GifFileIn -> Image.ColorMap
						 : GifFileIn -> SColorMap);

		*MaxX = GifFileIn -> Image.Width - 1;
		*MaxY = GifFileIn -> Image.Height - 1;
		
		/* Allocate the image (padding it as well). */
		if (*Alpha)
		    *Alpha = TranspColorIndex >= 0;

		if (*Alpha)
		    Data = IritMalloc(4 * (GifFileIn -> Image.Width + 3) *
				          (GifFileIn -> Image.Height + 3));
		else
		    Data = IritMalloc(3 * (GifFileIn -> Image.Width + 3) *
				          (GifFileIn -> Image.Height + 3));
		Line = IritMalloc(GifFileIn -> Image.Width + 3);

		/* Read the image itself: */

		if (GifFileIn -> Image.Interlace) {
		    /* Interlaced images - offsets and jumps. */
		    static int
			InterlacedOffset[] = { 0, 4, 2, 1 },
			InterlacedJumps[] = { 8, 8, 4, 2 };

		    /* Need to perform 4 passes on the images: */
		    for (i = 0; i < 4; i++) {
			for (l = InterlacedOffset[i];
			     l < GifFileIn -> Image.Height;
			     l += InterlacedJumps[i]) {
			    if (DGifGetLine(GifFileIn, Line,
				      GifFileIn -> Image.Width) == GIF_ERROR) {
				IritFree(Data);
				return NULL;
			    }

			    p = &Data[(GifFileIn -> Image.Height - l - 1) *
				      GifFileIn -> Image.Width *
				      ((*Alpha) ? 4 : 3)];
			    for (j = 0; j < GifFileIn -> Image.Width; j++) {
			        ColorMapEntry = &ColorMap -> Colors[Line[j]];
				*p++ = ColorMapEntry -> Red;
				*p++ = ColorMapEntry -> Green;
				*p++ = ColorMapEntry -> Blue;
				if (*Alpha)
				    *p++ = Line[j] == TranspColorIndex ? 0 : 255;
			    }
			}
		    }
		}
		else {
		    for (i = 0; i < GifFileIn -> Image.Height; i++) {
		        if (DGifGetLine(GifFileIn, Line,
				      GifFileIn -> Image.Width) == GIF_ERROR) {
			    IritFree(Data);
			    return NULL;
			}

			p = &Data[(GifFileIn -> Image.Height - i - 1) *
				  GifFileIn -> Image.Width *
				  ((*Alpha) ? 4 : 3)];
			for (j = 0; j < GifFileIn -> Image.Width; j++) {
			    ColorMapEntry = &ColorMap -> Colors[Line[j]];
			    *p++ = ColorMapEntry -> Red;
			    *p++ = ColorMapEntry -> Green;
			    *p++ = ColorMapEntry -> Blue;
			    if (*Alpha)
			        *p++ = Line[j] == TranspColorIndex ? 0 : 255;
			}
		    }
		}
		IritFree(Line);
		break;
	    case EXTENSION_RECORD_TYPE:
		/* Skip extension blocks in file: */
		DGifGetExtension(GifFileIn, &ExtCode, &Extension);
		while (Extension != NULL) {
		    if (ExtCode == 249 &&
			Extension[0] == 4 &&
			(Extension[1] & 0x01) != 0) {
		        /* We have a transp. color specification - get it. */
		        TranspColorIndex = Extension[4];
		    }

		    DGifGetExtensionNext(GifFileIn, &Extension);
		}
		break;
	    case TERMINATE_RECORD_TYPE:
		break;
	    default:		    /* Should be traps by DGifGetRecordType. */
		break;
	}
    }
    while (RecordType != TERMINATE_RECORD_TYPE);

    DGifCloseFile(GifFileIn);

    Data = _IrtImgVerifyAlignment(Data, MaxX, MaxY, *Alpha);

    return Data;
#else
    IRIT_WARNING_MSG("GifLib tool kit is not supported\n");
    return NULL;			   /* Silent the compiler's warning. */
#endif /* IRIT_HAVE_GIF_LIB */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads image file in PNG format.  Based on the example from libpng.       *
*                                                                            *
* PARAMETERS:                                                                *
*   ImageFileName:   Name of PNG image to read.                              *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType *:  Pointer to dynamicaly created image, NULL if non.           *
*****************************************************************************/
static IrtBType *PNGReadImage(const char *ImageFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha)
{
#ifdef IRIT_HAVE_PNG_LIB
    int y, RowSize;
    IrtBType *Data;
    png_structp PngPtr;
    png_infop InfoPtr;
    unsigned int
	SigRead = 0;
    FILE *fp;

#if defined(__WINNT__) || defined(__WINCE__)
    if ((fp = fopen(ImageFileName, "rb")) == NULL) {
#else
    if ((fp = fopen(ImageFileName, "r")) == NULL) {
#endif /* __WINNT__ || __WINCE__ */
        IRIT_WARNING_MSG_PRINTF("Failed to read PNG file \"%s\"\n",
				ImageFileName);
        return NULL;
    }

    PngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING, png_voidp_NULL, 
				    png_error_ptr_NULL, png_error_ptr_NULL);
    if (PngPtr == NULL) {
        fclose(fp);
	return NULL;
    }

    /* Allocate/initialize the memory for image information.  REQUIRED. */
    InfoPtr = png_create_info_struct(PngPtr);
    if (InfoPtr == NULL) {
	fclose(fp);
	png_destroy_read_struct(&PngPtr, png_infopp_NULL, png_infopp_NULL);
	return NULL;
    }

    if (setjmp(png_jmpbuf(PngPtr))) {
        /* Free all of the memory associated with the PngPtr and InfoPtr. */
        png_destroy_read_struct(&PngPtr, &InfoPtr, png_infopp_NULL);
	fclose(fp);
	/* If we get here, we had a problem reading the file. */
	return NULL;
    }

    png_init_io(PngPtr, fp);

    /* If we have already read some of the signature. */
    png_set_sig_bytes(PngPtr, SigRead);

    /* Strip 16 bit/color files down to 8 bits/color. */
    png_set_strip_16(PngPtr);

    png_read_png(PngPtr, InfoPtr,
		 (*Alpha ? 0 : PNG_TRANSFORM_STRIP_ALPHA) |
		               PNG_TRANSFORM_EXPAND,
		 png_voidp_NULL);

    if (*Alpha)
        *Alpha = InfoPtr -> color_type == PNG_COLOR_TYPE_RGB_ALPHA;

    *MaxX = InfoPtr -> width - 1;
    *MaxY = InfoPtr -> height - 1;
    RowSize = InfoPtr -> rowbytes;
    Data = IritMalloc(RowSize * (*MaxY + 1));
    if (InfoPtr -> bit_depth == 8) {          /* 8 bits per color channel. */
        for (y = 0; y <= *MaxY; y++) {
	    IRIT_GEN_COPY(&Data[RowSize * y],
			  InfoPtr -> row_pointers[*MaxY - y],
			  RowSize);
	}
    }
    else {
        IRIT_WARNING_MSG("PNG Image is not 8 bits per color channel.");
	return NULL;
    }

    /* clean up after the read, and free any memory allocated - REQUIRED */
    png_destroy_read_struct(&PngPtr, &InfoPtr, png_infopp_NULL);

    /* close the file */
    fclose(fp);

    Data = _IrtImgVerifyAlignment(Data, MaxX, MaxY, *Alpha);

    /* that's it */
    return Data;
#else
    IRIT_WARNING_MSG("LibPng tool kit is not supported\n");
    return NULL;			   /* Silent the compiler's warning. */
#endif /* IRIT_HAVE_PNG_LIB */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads image file in JPEG format.                                         *
*                                                                            *
* PARAMETERS:                                                                *
*   ImageFileName:   Name of JPG image to read.                              *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType *:  Pointer to dynamicaly created image, NULL if non.           *
*****************************************************************************/
static IrtBType *JPEGReadImage(const char *ImageFileName,
			       int *MaxX,
			       int *MaxY,
			       int *Alpha)
{
#ifdef IRIT_HAVE_JPG_LIB
    FILE *JPEGLoad;
    struct jpeg_decompress_struct ContextInfo;
    struct jpeg_error_mgr JErrorManager;
    IrtBType *Line, *Data, *TmpData;

    ContextInfo.err = jpeg_std_error(&JErrorManager);
    jpeg_create_decompress(&ContextInfo);
	
#if defined(__WINNT__) || defined(__WINCE__)
    if ((JPEGLoad = fopen(ImageFileName, "rb")) == NULL) {
#else
    if ((JPEGLoad = fopen(ImageFileName, "r")) == NULL) {
#endif
        IRIT_WARNING_MSG_PRINTF("Failed to read Jpeg file \"%s\"\n",
				ImageFileName);
	return NULL;
    }
	
    jpeg_stdio_src(&ContextInfo, JPEGLoad);
    jpeg_read_header(&ContextInfo, TRUE);
	
    *MaxX = ContextInfo.image_width - 1;
    *MaxY = ContextInfo.image_height - 1;
    *Alpha = FALSE;

    jpeg_start_decompress(&ContextInfo);

    Data = (IrtBType *) IritMalloc(3 * ContextInfo.image_width *
				       ContextInfo.image_height);

    while (ContextInfo.output_scanline < ContextInfo.output_height) {
        Line = Data + 3 * ContextInfo.image_width *
	                  ContextInfo.output_scanline;
	jpeg_read_scanlines(&ContextInfo, &Line, 1);
    }

    Data = _IrtImgVerifyAlignment(Data, MaxX, MaxY, FALSE);
    TmpData = (IrtBType *)
        IrtImgFlipVerticallyImage((IrtImgPixelStruct *) Data, *MaxX, *MaxY,
				  FALSE);
    IritFree(Data);
    Data = TmpData;

    jpeg_finish_decompress(&ContextInfo);
    jpeg_destroy_decompress(&ContextInfo);
    return Data;
#else
    IRIT_WARNING_MSG("LibJpeg tool kit is not supported\n");
    return NULL;			   /* Silent the compiler's warning. */
#endif /* IRIT_HAVE_JPG_LIB */
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Parses the string of the "ptexture" attribute.			     M
* "ImageName {, S X Y {Z}} {, F} {, N}, {, H}, {, V}"      where	     M
* 1. X, Y, and possibly Z are scaling factors of how many times the image    M
*    should fit into the object,					     M
* 2. 'F' optionally request to flip X and Y axes of the image.		     M
* 3. 'N' optionally forces reload the image as a New image, even if an image M
*    by this exact same name was already loaded and cached.	             M
* 4. 'H' optionally request to filp the image pixels horizontally.	     M
* 5. 'V' optionally request to flip the image pixels vertically		     M
*                                                                            *
* PARAMETERS:                                                                M
*   PTexture:		The string of the "ptexture" attribute.              M
*   FName:		The texture file name will be placed here.	     M
*   Scale:		The scaling vector in XYZ or just XY if		     M
*			Z = IRIT_INFNTY.				     M
*   Flip:		If Image flipping was requested.		     M
*   NewImage:		if a new image of same name - replace in cache.      M
*   FlipHorizontally:	If horizontal Image flipping was requested.	     M
*   FlipVertically:	If vertical Image flipping was requested.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:         TRUE if parsed succesfully, FALSE otherwise.                M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgParsePTextureString2                                               M
*****************************************************************************/
int IrtImgParsePTextureString2(const char *PTexture,
			      char *FName,
			      IrtRType *Scale,
			      int *Flip,
			      int *NewImage,
			      int *FlipHorizontally,
			      int *FlipVertically)
{
    char *p;

    Scale[0] = Scale[1] = 1.0;
    Scale[2] = IRIT_INFNTY;
    *Flip = FALSE;
    *FlipHorizontally = FALSE;
    *FlipVertically = FALSE;
    *NewImage = FALSE;

    if (PTexture == NULL)
	return FALSE;

    strncpy(FName, PTexture, IRIT_LINE_LEN_LONG - 1);

    if ((p = strchr(FName, ',')) != NULL) {
	char *q;
	float Sx, Sy, Sz;

	*p++ = 0;		      /* Mark the end of the regular string. */

	if ((q = strchr(p, 'f')) != NULL || (q = strchr(p, 'F')) != NULL)
	    *Flip = TRUE;

	if ((q = strchr(p, 'h')) != NULL || (q = strchr(p, 'H')) != NULL)
	    *FlipHorizontally = TRUE;
	
	if ((q = strchr(p, 'v')) != NULL || (q = strchr(p, 'V')) != NULL)
	    *FlipVertically = TRUE;

	if ((q = strchr(p, 'n')) != NULL || (q = strchr(p, 'N')) != NULL)
	    *NewImage = TRUE;

	if (sscanf(p, "%f, %f, %f", &Sx, &Sy, &Sz) == 3 ||
	    (((q = strchr(p, 's')) != NULL ||
	      (q = strchr(p, 'S')) != NULL) &&
	     sscanf(q, "S %f %f %f", &Sx, &Sy, &Sz) == 3)) {
	    Scale[0] = Sx;
	    Scale[1] = Sy;
	    Scale[2] = Sz;
	}
	else if (sscanf(p, "%f, %f", &Sx, &Sy) == 2 ||
		 (((q = strchr(p, 's')) != NULL ||
		   (q = strchr(p, 'S')) != NULL) &&
		  sscanf(q, "S %f %f", &Sx, &Sy) == 2)) {
	    Scale[0] = Sx;
	    Scale[1] = Sy;
	}
    }

    return TRUE;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Parses the string of the "ptexture" attribute.			     M
* "ImageName {, S X Y {Z}} {, F} {, N}"      where			     M
* 1. X, Y, and possibly Z are scaling factors of how many times the image    M
*    should fit into the object,					     M
* 2. 'F' optionally request to flip X and Y axes of the image.		     M
* 3. 'N' optionally forces reload the image as a New image, even if an image M
*    by this exact same name was already loaded and cached.	             M
*                                                                            *
* PARAMETERS:                                                                M
*   PTexture:    The string of the "ptexture" attribute.                     M
*   FName:       The texture file name will be placed here.		     M
*   Scale:       The scaling vector in XYZ or just XY if Z = IRIT_INFNTY.    M
*   Flip:	 If Image flipping was requested.			     M
*   NewImage:    if a new image of same name - replace in cache.             M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:         TRUE if parsed succesfully, FALSE otherwise.                M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtImgParsePTextureString                                                M
*****************************************************************************/
int IrtImgParsePTextureString(const char *PTexture,
			      char *FName,
			      IrtRType *Scale,
			      int *Flip,
			      int *NewImage)
{
    int Flip1 = FALSE,
        Flip2 = FALSE;

    return IrtImgParsePTextureString2(PTexture, FName, Scale, Flip, NewImage, 
				      &Flip1, &Flip2);
}
