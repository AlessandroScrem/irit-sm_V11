/*****************************************************************************
* Reads one movie in.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber			       Ver 1.0,	Dec. 2012    *
*****************************************************************************/

#include <math.h>
#include "inc_irit/irit_sm.h"
#include "misc_loc.h"

#ifdef IRIT_HAVE_GIF_LIB
#include "gif_lib.h"
#endif /* IRIT_HAVE_GIF_LIB */

#ifdef __WINNT__
#include <io.h>
#include <fcntl.h>
#endif /* __WINNT__ */

#ifdef IRIT_HAVE_FFMPEG

#define inline __inline
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>

static int  FFMPEGFrameImageConvert(AVPicture *DestinationFrame,
				    enum PixelFormat DestinationFormat,
				    int DestinationWidth,
				    int DestinationHeight,
				    AVPicture *SourceFrame,
				    enum PixelFormat SourceFormat,
				    int SourceWidth,
				    int SourceHeight);
static void FFMPEGSaveFrame(AVFrame *DstFrame,
			    int Width,
			    int Height,
			    IrtImgPixelStruct **Frames);

#define FFMPEG_PRINT_ERROR_EXTRA(Msg, ResultCode, Extra) { \
    char ErrorStr[IRIT_LINE_LEN]; \
    IRIT_WARNING_MSG_PRINTF( \
	    IRIT_EXP_STR("FFMPEG ERROR: %s. " Msg ".\n"), \
	    av_make_error_string(ErrorStr, IRIT_LINE_LEN, ResultCode), \
	    Extra); }

#define FFMPEG_PRINT_ERROR(Msg, ResultCode) { \
    char ErrorStr[IRIT_LINE_LEN]; \
    IRIT_WARNING_MSG_PRINTF( \
	    IRIT_EXP_STR("FFMPEG ERROR: %s. " Msg ".\n"), \
	    av_make_error_string(ErrorStr, IRIT_LINE_LEN, ResultCode)); }

IRIT_STATIC_DATA int
    FFMPEGIsInitialized = FALSE;

#endif/* IRIT_HAVE_FFMPEG */

typedef struct LoadedMoviesStruct {
    struct LoadedMoviesStruct *Pnext;
    char *FileName;
    int MaxX, MaxY, Alpha;
    IrtBType **Data;
} LoadedMoviesStruct;

IRIT_STATIC_DATA LoadedMoviesStruct
    *GlblLoadedMoviesList = NULL;

static IrtBType **GIFReadMovie(const char *MovieFileName,
			      int *MaxX,
			      int *MaxY,
			       int *Alpha);
static IrtImgPixelStruct **FFMPEGReadMovie(const char *MovieFileName,
					   int *MaxX,
					   int *MaxY,
					   int *Alpha);
static int FFMPEGGetMovieProps(const char *MovieFileName,
			       int *Width,
			       int *Height);

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Reads a movie in from a file named MovieFileName.  The movies are        M
* returned as a vector of images: vector of RGBRGB... of size		     M
*                                                  (MaxX+1) * (MaxY+1) * 3.  M
*                                                                            *
* PARAMETERS:                                                                M
*   MovieFileName:   Name of movie to read.                                  M
*   MaxX:            Maximum X of read image is saved here.                  M
*   MaxY:            Maximum Y of read image is saved here.                  M
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  M
*		     and will return TRUE if successful in loading it.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct **:  A NULL terminated vector of images.  Each image   M
*               is itself a vector of RGBRGB... of size 		     M
*		(MaxX+1) * (MaxY+1) * 3 or NULL if failed.		     M
*		  If however, Alpha is requested and found RGBARGBA... is    M
*		returned as IrtImgRGBAPxlStruct.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtMovieReadMovie2, IrtMovieReadMovie3, IrtImgReadMovieXAlign            M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtMovieReadMovie                                                        M
*****************************************************************************/
IrtImgPixelStruct **IrtMovieReadMovie(const char *MovieFileName,
				      int *MaxX,
				      int *MaxY,
				      int *Alpha)
{
    const char *Type;

    if (MovieFileName == NULL) {
        IRIT_FATAL_ERROR("Empty movie file name to write to.");
        return NULL;
    }

    if ((Type = strrchr(MovieFileName, '.')) == NULL)
        Type = "";

    if (stricmp(Type, ".Z") == 0) {
        Type--;
	while (Type != MovieFileName && *Type != '.')
	    Type--;
    }

    if (stricmp(Type, ".gif") == 0) {
	return (IrtImgPixelStruct **) GIFReadMovie(MovieFileName, MaxX, MaxY,
						   Alpha);
    }
    else if (stricmp(Type, ".wmv") == 0 ||
	     stricmp(Type, ".avi") == 0 ||
	     stricmp(Type, ".mp4") == 0 ||
	     stricmp(Type, ".mkv") == 0 ||
	     stricmp(Type, ".mpg") == 0 ||
	     stricmp(Type, ".mpeg") == 0) {
	return FFMPEGReadMovie(MovieFileName, MaxX, MaxY, Alpha);
    }
    else {
	IRIT_WARNING_MSG_PRINTF(
	     IRIT_EXP_STR("Texture movie file \"%s\" must be of type 'gif', 'avi', 'wmv', 'mkv', 'mp4', 'mpg', or 'mpeg'\n"),
	     MovieFileName);
        return NULL;
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Same as IrtMovieReadMovie but if a name of a movie repeats itself, the   M
* movie is read only once.						     M
*                                                                            *
* PARAMETERS:                                                                M
*   MovieFileName:   Name of movie to read.                                  M
*   MaxX:            Maximum X of read image is saved here.                  M
*   MaxY:            Maximum Y of read image is saved here.                  M
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  M
*		     and will return TRUE if successful in loading it.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct **:  A NULL terminated vector of images.  Each image   M
*               is itself a vector of RGBRGB... of size 		     M
*		(MaxX+1) * (MaxY+1) * 3 or NULL if failed.		     M
*		  If however, Alpha is requested and found RGBARGBA... is    M
*		returned as IrtImgRGBAPxlStruct.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtMovieReadMovie, IrtImgReadMovie3, IrtMovieReadClrCache                M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtMovieReadMovie2                                                       M
*****************************************************************************/
IrtImgPixelStruct **IrtMovieReadMovie2(const char *MovieFileName,
				       int *MaxX,
				       int *MaxY,
				       int *Alpha)
{
    LoadedMoviesStruct *LoadedMovie;
    IrtImgPixelStruct **Data;

    /* Search if we already loaded this image. */
    for (LoadedMovie = GlblLoadedMoviesList;
	 LoadedMovie != NULL;
	 LoadedMovie = LoadedMovie -> Pnext) {
	if (strcmp(MovieFileName, LoadedMovie -> FileName) == 0) {
	    *MaxX = LoadedMovie -> MaxX;
	    *MaxY = LoadedMovie -> MaxY;
	    *Alpha = LoadedMovie -> Alpha;
	    return (IrtImgPixelStruct **) LoadedMovie -> Data;
	}
    }

    if ((Data = IrtMovieReadMovie(MovieFileName, MaxX, MaxY, Alpha)) != NULL) {
	/* Add it to global list of loaded images. */
	LoadedMovie = (LoadedMoviesStruct *)
				    IritMalloc(sizeof(LoadedMoviesStruct));
	LoadedMovie -> FileName = IritStrdup(MovieFileName);
	LoadedMovie -> MaxX = *MaxX;
	LoadedMovie -> MaxY = *MaxY;
	LoadedMovie -> Alpha = *Alpha;
	LoadedMovie -> Data = (IrtBType **) Data;
	LoadedMovie -> Pnext = GlblLoadedMoviesList;
	GlblLoadedMoviesList = LoadedMovie;
    }

    return Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Same as IrtMovieReadMovie2 but if a name of a movie repeats itself, the  M
* new movie replaces the old one.					     M
*                                                                            *
* PARAMETERS:                                                                M
*   MovieFileName:   Name of movie to read.                                  M
*   MaxX:            Maximum X of read image is saved here.                  M
*   MaxY:            Maximum Y of read image is saved here.                  M
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  M
*		     and will return TRUE if successful in loading it.       M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtImgPixelStruct **:  A NULL terminated vector of images.  Each image   M
*               is itself a vector of RGBRGB... of size 		     M
*		(MaxX+1) * (MaxY+1) * 3 or NULL if failed.		     M
*		  If however, Alpha is requested and found RGBARGBA... is    M
*		returned as IrtImgRGBAPxlStruct.			     M
*                                                                            *
* SEE ALSO:                                                                  M
*   IrtMovieReadMovie, IrtImgReadMovie2, IrtMovieReadClrCache                M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtMovieReadMovie3                                                       M
*****************************************************************************/
IrtImgPixelStruct **IrtMovieReadMovie3(const char *MovieFileName,
				       int *MaxX,
				       int *MaxY,
				       int *Alpha)
{
    LoadedMoviesStruct *LoadedMovie;
    IrtImgPixelStruct **Data;

    if ((Data = IrtMovieReadMovie(MovieFileName, MaxX, MaxY, Alpha)) != NULL) {
        /* Search if we already loaded this movie. */
        for (LoadedMovie = GlblLoadedMoviesList;
	     LoadedMovie != NULL;
	     LoadedMovie = LoadedMovie -> Pnext) {
	    if (strcmp(MovieFileName, LoadedMovie -> FileName) == 0) {
		int i;

	        IritFree(LoadedMovie -> FileName);
		for (i = 0; LoadedMovie -> Data[i] != NULL; i++)
		    IritFree(LoadedMovie -> Data[i]);
		IritFree(LoadedMovie -> Data);
		break;
	    }
	}

	if (LoadedMovie == NULL) {                 /* Is a new image name. */
	    LoadedMovie = (LoadedMoviesStruct *)
				    IritMalloc(sizeof(LoadedMoviesStruct));
	    LoadedMovie -> Pnext = GlblLoadedMoviesList;
	    GlblLoadedMoviesList = LoadedMovie;
	}

	LoadedMovie -> FileName = IritStrdup(MovieFileName);
	LoadedMovie -> MaxX = *MaxX;
	LoadedMovie -> MaxY = *MaxY;
	LoadedMovie -> Alpha = *Alpha;
	LoadedMovie -> Data = (IrtBType **) Data;
	LoadedMovie -> Pnext = GlblLoadedMoviesList;
	GlblLoadedMoviesList = LoadedMovie;
    }

    return Data;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads the properties of a movie file.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   MovieFileName:   Name of movie video file to read.                       *
*   Width:           The width of the video frames.			     *
*   Height:          The height of the video frames.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:	     TRUE, if loaded successfully. otherwise - FALSE.	     *
*****************************************************************************/
int IrtMovieGetMovieProps(const char *MovieFileName,
			  int *Width,
			  int *Height)
{
    const char *Type;
    int Alpha;

    if (MovieFileName == NULL) {
        IRIT_FATAL_ERROR("Empty movie file name to read from.");
	return FALSE;
    }

    if ((Type = strrchr(MovieFileName, '.')) == NULL)
        Type = "";

    if (stricmp(Type, ".Z") == 0) {
        Type--;
	while (Type != MovieFileName && *Type != '.')
	    Type--;
    }

    if (stricmp(Type, ".gif") == 0) {
	IrtImgPixelStruct 
	    **Video = (IrtImgPixelStruct **) GIFReadMovie(MovieFileName,
							  Width, Height,
							  &Alpha);

	if (Video == NULL) {
	    IRIT_WARNING_MSG_PRINTF(
	         IRIT_EXP_STR("Gif file \"%s\" could not be read.'\n"),
		 MovieFileName);
	    return FALSE;
	}

	(*Height)++;
	(*Width)++;
	IritFree(Video);
    }
    else if (stricmp(Type, ".wmv") == 0 ||
	     stricmp(Type, ".avi") == 0 ||
	     stricmp(Type, ".mp4") == 0 ||
	     stricmp(Type, ".mkv") == 0 ||
	     stricmp(Type, ".mpg") == 0 ||
	     stricmp(Type, ".mpeg") == 0) {
	return FFMPEGGetMovieProps(MovieFileName, Width, Height);
    }
    else {
	IRIT_WARNING_MSG_PRINTF(
	     IRIT_EXP_STR("Texture movie file \"%s\" must be of type 'gif', 'avi', 'wmv', 'mkv', 'mp4', 'mpg', or 'mpeg'\n"),
	     MovieFileName);
	return FALSE;
    }

    return TRUE;
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
*   IrtMovieReadMovie2                                                       M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtMovieReadClrCache                                                     M
*****************************************************************************/
void IrtMovieReadClrCache(void)
{
    while (GlblLoadedMoviesList != NULL) {
        int i;
	LoadedMoviesStruct
	    *LoadedMovie = GlblLoadedMoviesList;

	GlblLoadedMoviesList = GlblLoadedMoviesList -> Pnext;
	IritFree(LoadedMovie -> FileName);

	for (i = 0; LoadedMovie -> Data[i] != NULL; i++)
	    IritFree(LoadedMovie -> Data[i]);
	IritFree(LoadedMovie -> Data);

	IritFree(LoadedMovie);
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads image file in GIF format.                                          *
*                                                                            *
* PARAMETERS:                                                                *
*   MovieFileName:   Name of GIF image to read.                              *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtBType **:  Pointer to dynamicaly created movie - a NULL terminated    *
*                 vector of images, NULL if non.		             *
*****************************************************************************/
static IrtBType **GIFReadMovie(const char *MovieFileName,
			      int *MaxX,
			      int *MaxY,
			      int *Alpha)
{
#ifdef IRIT_HAVE_GIF_LIB
    int MaxNumImages = 3,
        NumImages = 0,
        TranspColorIndex = -1;
    GifFileType *GifFileIn;
    IrtBType *Data, *p, *Line, **Movie;
    GifRecordType RecordType;

    if ((GifFileIn = DGifOpenFileName(MovieFileName)) == NULL) {
        IRIT_WARNING_MSG_PRINTF("Failed to read GIF file \"%s\"\n",
				MovieFileName);
        return NULL;
    }

    Movie = (IrtBType **) IritMalloc(sizeof(IrtBType *) * MaxNumImages);

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

		Data = _IrtImgVerifyAlignment(Data, MaxX, MaxY, *Alpha);

		/* Push the next image on the movie's list of images. */
		if (NumImages > MaxNumImages - 2) {
		    Movie = (IrtBType **) IritRealloc(Movie,
				        sizeof(IrtBType *) * MaxNumImages,
				        sizeof(IrtBType *) * MaxNumImages * 2);
		    MaxNumImages *= 2;		    
		}
		Movie[NumImages++] = Data;
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

    Movie[NumImages] = NULL;
    return Movie;
#else
    IRIT_WARNING_MSG("GifLib tool kit is not supported\n");
    return NULL;			   /* Silent the compiler's warning. */
#endif /* IRIT_HAVE_GIF_LIB */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads video file in ANY format that FFMPEG library supports.             *
*   Namely: MPG, MPEG, AVI, WMV, MKV, MP4.				     *
*                                                                            *
* PARAMETERS:                                                                *
*   MovieFileName:   Name of video file to read.                             *
*   MaxX:            Maximum X of read image is saved here.                  *
*   MaxY:            Maximum Y of read image is saved here.                  *
*   Alpha:	     If TRUE, will attempt to load alpha channel if has any  *
*		     and will return TRUE if successful in loading it.	     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtImgPixelStruct **:  Pointer to dynamicaly created movie - a NULL      *
*                 terminated vector of images, NULL if non.		     *
*****************************************************************************/
static IrtImgPixelStruct **FFMPEGReadMovie(const char *MovieFileName,
					   int *MaxX,
					   int *MaxY,
					   int *Alpha)
{
#ifdef IRIT_HAVE_FFMPEG	
    AVFormatContext
	*FormatContext  = NULL;
    IrtImgPixelStruct **Movie;
    AVCodecContext				     /* Video identifiers: */
	*CodecContext   = NULL;
    AVCodec 
	*Codec = NULL;
    uint8_t				      /* Implementation variables: */
	*Buffer = NULL;
    AVFrame
	*SourceFrame, *DestinationFrame;
    IRIT_STATIC_DATA int IsFrameFinished;
    AVPacket Packet;
    int	Height, Width, 
	VideoStream = -1,
	VideoStreamIndex = 0,	
	NumFrames = 0, 
	MaxFrames = 4,
	ResultCode = 0;
    enum PixelFormat 
	DestinationFormat = PIX_FMT_RGB24;

    /* Initialize FFMPEG: */
    if (!FFMPEGIsInitialized) {
	/* Register all formats and codecs. */
	av_register_all();
	FFMPEGIsInitialized = TRUE;
    }

    /* Open video file: */
    if ((ResultCode = avformat_open_input(&FormatContext, MovieFileName,
					  NULL, NULL)) != 0) {
	FFMPEG_PRINT_ERROR_EXTRA("Couldn't open file %s", ResultCode,
				 MovieFileName);
	return NULL;
    }

    /* Retrieve stream information: */
    if ((ResultCode = avformat_find_stream_info(FormatContext, NULL)) < 0) {
	FFMPEG_PRINT_ERROR("Couldn't find stream information", ResultCode);
	return NULL;
    }

    /* Find the first video/audio stream. */
    for (VideoStreamIndex = 0;
	(unsigned int) VideoStreamIndex < FormatContext -> nb_streams &&
	VideoStream < 0;
	VideoStreamIndex++) {
	if (FormatContext -> streams[VideoStreamIndex] -> codec -> codec_type
	                                             == AVMEDIA_TYPE_VIDEO) {
	    VideoStream = VideoStreamIndex;
	}
    }
    if (VideoStream == -1) {
	IRIT_WARNING_MSG_PRINTF("FFMPEG ERROR: Couldn't find video stream.\n");
	return NULL;
    }

    /* Get a pointer to the codec context for the video and audio streams. */
    CodecContext = FormatContext -> streams[VideoStream] -> codec;

    /* Find the decoder for the video stream. */
    Codec = avcodec_find_decoder(CodecContext -> codec_id);
    if (!Codec) {
	IRIT_WARNING_MSG_PRINTF(
	    IRIT_EXP_STR("FFMPEG ERROR: Unsupported codec id: %d.\n"),
			 CodecContext->codec_id);
	return NULL;
    }

    /* Open codec */
    if ((ResultCode = avcodec_open2(CodecContext, Codec, NULL)) < 0) {
	FFMPEG_PRINT_ERROR_EXTRA("Couldn't open codec: %s", ResultCode,
				 CodecContext -> codec_name);
	return NULL;
    }

    /* Allocate video frame. */
    SourceFrame = av_frame_alloc();
    DestinationFrame = av_frame_alloc();

    Height = CodecContext -> height;
    Width = CodecContext -> width;

    /* Allocate the first movie's frames. */
    Movie = (IrtImgPixelStruct **) IritMalloc(sizeof(IrtBType *) * MaxFrames);

    Buffer = (uint8_t *) IritMalloc(sizeof(uint8_t) *
				    avpicture_get_size(DestinationFormat,
						       Width, Height));

    avpicture_fill((AVPicture *) DestinationFrame, Buffer,
		   DestinationFormat, Width, Height);

    while (av_read_frame(FormatContext, &Packet) >= 0) {
	/* Is this a Packet from the video stream? */
	if (Packet.stream_index == VideoStream) {
	    /* Decode video frame. */
	    avcodec_decode_video2(CodecContext, SourceFrame, &IsFrameFinished,
				  &Packet);
	    /* Did we get a video frame? */
	    if (IsFrameFinished) {
		/* Convert the image from its native format to RGB. */
		FFMPEGFrameImageConvert((AVPicture *) DestinationFrame,
					DestinationFormat, Width, Height,
					(AVPicture *) SourceFrame,
					CodecContext -> pix_fmt,
					Width, Height);
    		
		/* Push the next frame on the movie's list of frames. */
		if (NumFrames > MaxFrames - 2) {
		    Movie = (IrtImgPixelStruct **)
		            IritRealloc(Movie,
					sizeof(IrtBType *) * MaxFrames,
					sizeof(IrtBType *) * MaxFrames * 2);
		    MaxFrames *= 2;		    
		}

		FFMPEGSaveFrame(DestinationFrame, Width, Height,
				Movie + NumFrames++);
	    }
	}
        
	/* Free the Packet that was allocated by av_read_frame */
	av_free_packet(&Packet);
    }

    Movie[NumFrames] = NULL;

    IritFree(Buffer);

    /* Free the video frame. */
    av_free(SourceFrame);
    av_free(DestinationFrame);

    /* Close the codec. */
    avcodec_close(CodecContext);

    /* Close the video file. */
    avformat_close_input(&FormatContext);

    *MaxX = Width - 1;
    *MaxY = Height - 1;
    *Alpha = FALSE;

    return Movie;
#else
    IRIT_WARNING_MSG("FFMPEG library is not supported\n");
    return NULL;			   /* Silent the compiler's warning. */
#endif /* IRIT_HAVE_FFMPEG */
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Loads the properties of a movie video file in ANY format that FFMPEG     *
*   library supports. Namely: MPG, MPEG, AVI, WMV, MKV, MP4.		     *
*                                                                            *
* PARAMETERS:                                                                *
*   MovieFileName:   Name of video file to read.                             *
*   Width:           The width of the video frames.			     *
*   Height:          The height of the video frames.			     *
*                                                                            *
* RETURN VALUE:                                                              *
*   int:	     TRUE, if loaded successfully. otherwise - FALSE.	     *
*****************************************************************************/
static int FFMPEGGetMovieProps(const char *MovieFileName,
			       int *Width,
			       int *Height)
{
#ifdef IRIT_HAVE_FFMPEG
    AVFormatContext
	*FormatContext = NULL;
    AVCodecContext    /* Video identifiers: */
	*CodecContext = NULL;
    AVCodec 
	*Codec = NULL;
    int	VideoStream = -1,    /* Implementation variables: */
	VideoStreamIndex = 0,
	ResultCode = 0;

    /* Initialize FFMPEG: */
    if (!FFMPEGIsInitialized) {
	/* Register all formats and codecs. */
	av_register_all();
	FFMPEGIsInitialized = TRUE;
    }

    /* Open video file: */
    if ((ResultCode = 
	avformat_open_input(&FormatContext, MovieFileName,
			    NULL, NULL)) != 0) {
	FFMPEG_PRINT_ERROR_EXTRA("Couldn't open file %s", ResultCode,
				 MovieFileName);
	return FALSE;
    }

    /* Retrieve stream information: */
    if ((ResultCode = avformat_find_stream_info(FormatContext, NULL)) < 0) {
	FFMPEG_PRINT_ERROR("Couldn't find stream information", ResultCode);
	return FALSE;
    }

    /* Find the first video/audio stream. */
    for (VideoStreamIndex=0;
	 (unsigned int) VideoStreamIndex < FormatContext -> nb_streams &&
	     VideoStream < 0;
	 VideoStreamIndex++) {
	if (FormatContext -> streams[VideoStreamIndex] -> codec -> codec_type
	                                             == AVMEDIA_TYPE_VIDEO) {
	    VideoStream = VideoStreamIndex;
	}
    }
    if (VideoStream == -1) {
	IRIT_WARNING_MSG_PRINTF("FFMPEG ERROR: Couldn't find video stream.\n");
	return FALSE;
    }

    /* Get a pointer to the codec context for the video and audio streams. */
    CodecContext = FormatContext->streams[VideoStream] -> codec;

    /* Find the decoder for the video stream */
    Codec = avcodec_find_decoder(CodecContext -> codec_id);
    if (!Codec) {
	IRIT_WARNING_MSG_PRINTF(
	    IRIT_EXP_STR("FFMPEG ERROR: Unsupported codec id: %d.\n"),
			 CodecContext->codec_id);
	return FALSE;
    }

    /* Open codec. */
    if ((ResultCode = avcodec_open2(CodecContext, Codec, NULL))<0) {
	FFMPEG_PRINT_ERROR_EXTRA("Couldn't open codec: %s", ResultCode,
				 CodecContext->codec_name);
	return FALSE;
    }

    *Height = CodecContext -> height;
    *Width = CodecContext -> width;

    /* Close the codec. */
    avcodec_close(CodecContext);

    /* Close the video file. */
    avformat_close_input(&FormatContext);
    return TRUE;
#else
    IRIT_WARNING_MSG("FFMPEG library is not supported\n");
    return FALSE;
#endif /* IRIT_HAVE_FFMPEG */
}

/*****************************************************************************
* DESCRIPTION:                                                               M
*   Parses the string of the "pmovie" attribute.  Can be		     M
* "MovieName {, S X Y {Z}} {,F} {,R} {, T=tmin,tmax}	    where	     M
* 1. X, Y, and possibly Z are scaling factors of how many times the image    M
*    should fit into the object,					     M
* 2. F optionally request to flip X and Y axes of the image.		     M
* 3. R optionally request to restart the movie once it gets to the end.      M
* 4. T=tmin,tmax sets the time range of the movie.			     M
*                                                                            *
* PARAMETERS:                                                                M
*   PMovie:      The string of the "pmovie" attribute.                       M
*   FName:       The texture file name will be placed here.		     M
*   Scale:       The scaling vector in XYZ or just XY if Z = IRIT_INFNTY.    M
*   NewImage:    Force read an image even if in cache as the file is new?    M
*   Flip:	 If Image flipping was requested.			     M
*   Restart:     TRUE to restart the movie onces it ends.		     M
*   TimeSetup:   A vector of 3 slots to hold (tmin, tmax, dt) if specified   M
*                or all are zero if not specified.			     M
*   FlipHorizontally:	If horizontal Image flipping was requested.	     M
*   FlipVertically:	If vertical Image flipping was requested.	     M
*                                                                            *
* RETURN VALUE:                                                              M
*   int:         TRUE if parsed successfully, FALSE otherwise.               M
*                                                                            *
* KEYWORDS:                                                                  M
*   IrtMovieParsePMovieString                                                M
*****************************************************************************/
int IrtMovieParsePMovieString(const char *PMovie,
			      char *FName,
			      IrtRType *Scale,
			      int *NewImage,
			      int *Flip,
			      int *Restart,
			      IrtRType *TimeSetup,
			      int *FlipHorizontally,
			      int *FlipVertically)
{
    char *p;

    Scale[0] = Scale[1] = 1.0;
    Scale[2] = IRIT_INFNTY;
    *Flip = *Restart = *NewImage = FALSE;
    TimeSetup[0] = TimeSetup[1] = TimeSetup[2] = 0.0;

    if (PMovie == NULL)
	return FALSE;

    strncpy(FName, PMovie, IRIT_LINE_LEN_LONG - 1);

    if ((p = strchr(FName, ',')) != NULL) {
	char *q;
	float Sx, Sy, Sz;

	*p++ = 0;		      /* Mark the end of the regular string. */

	if ((q = strchr(p, 'f')) != NULL || (q = strchr(p, 'F')) != NULL)
	    *Flip = TRUE;

	if ((q = strchr(p, 'r')) != NULL || (q = strchr(p, 'R')) != NULL)
	    *Restart = TRUE;

	if ((q = strchr(p, 'n')) != NULL || (q = strchr(p, 'N')) != NULL)
	    *NewImage = TRUE;

	if ((q = strstr(p, "t=")) != NULL ||
	    (q = strstr(p, "T=")) != NULL) {
	    if (sscanf(&q[2], "%lf,%lf",
		       &TimeSetup[0], &TimeSetup[1]) != 2)
	        TimeSetup[0] = TimeSetup[1] = 0.0;
	}

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

#ifdef IRIT_HAVE_FFMPEG

/*****************************************************************************
* DESCRIPTION:                                                               *
*  Converts a video frame from the source format to the destination format.  *
*  This function is important for decoding frames and converting them to RGB.*
*                                                                            *
* PARAMETERS:                                                                *
*   DestinationFrame:	The output frame.                                    *
*   DestinationFormat:	The requested frame format.                          *
*   DestinationWidth:	The output frame's width.                            *
*   DestinationHeight:	The output frame's height.                           *
*   SourceFrame:	The input frame.                                     *
*   SourceFormat:	The input frame's format, to be changed.             *
*   SourceWidth:	The input frame's width.                             *
*   SourceHeight:	The input frame's height.                            *
*                                                                            *
* RETURN VALUE:                                                              *
*   int: TRUE - if the conversion is successful. Otherwise - FALSE.	     *
*****************************************************************************/
static int FFMPEGFrameImageConvert(AVPicture *DestinationFrame,
				   enum PixelFormat DestinationFormat,
				   int DestinationWidth,
				   int DestinationHeight,
				   AVPicture *SourceFrame,
				   enum PixelFormat SourceFormat,
				   int SourceWidth,
				   int SourceHeight)
{
    struct SwsContext *Context;

    /* create scaling context */
    Context = sws_getContext(SourceWidth, SourceHeight, SourceFormat,
                             DestinationWidth, DestinationHeight,
			     DestinationFormat, SWS_BILINEAR,
			     NULL, NULL, NULL);
    if (!Context) {
        IRIT_WARNING_MSG_PRINTF(IRIT_EXP_STR(
		"Impossible to create scale context for the conversion fmt:%s s:%dx%d -> fmt:%s s:%dx%d\n"),
                av_get_pix_fmt_name(SourceFormat), SourceWidth, SourceHeight,
                av_get_pix_fmt_name(DestinationFormat), DestinationWidth,
		DestinationHeight);
        return -1;
    }

    sws_scale(Context, (const uint8_t * const *) SourceFrame -> data,
	      SourceFrame -> linesize, 0, SourceHeight,
	      DestinationFrame -> data, DestinationFrame -> linesize);

    sws_freeContext(Context);
    return 0;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*    Adds the given AVFrame to the array of IrtImgPixelStruct frames.	     *
*                                                                            *
* PARAMETERS:                                                                *
*   DstFrame:	The frame to be saved.                                       *
*   Width:	The frame's width.                                           *
*   Height:	The frame's height.                                          *
*   Frames:	The array of IrtImgPixelStruct frames.                       *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void FFMPEGSaveFrame(AVFrame *DstFrame,
			    int Width,
			    int Height,
			    IrtImgPixelStruct **Frames)
{
    IrtImgPixelStruct
        *Frame = (IrtImgPixelStruct *) IritMalloc((Width + 3) * (Height + 3) *
						  sizeof(IrtImgPixelStruct));
    int i,j;

    /* Write pixel data. */
    for(i = 0; i < Height; i++) {
        for (j = 0; j < Width; j++) {
	    IrtImgPixelStruct
	        *IrtPixel = &Frame[(Height - 1 - i) * Width + j];
	    IrtBType
	        *FFPixel = DstFrame -> data[0] +
	                   i * DstFrame -> linesize[0] + j * 3;
		
	    IrtPixel->r = *(FFPixel++);
	    IrtPixel->g = *(FFPixel++);
	    IrtPixel->b = *(FFPixel++);
	}
    }
	
    Width--;
    Height--;
    Frame = (IrtImgPixelStruct *)_IrtImgVerifyAlignment((IrtBType *) Frame,
							&Width, &Height,
							FALSE);
    Frames[0] = Frame;
}

#endif /* IRIT_HAVE_FFMPEG */
