#pragma once

extern "C"
{
#include "jpeglib.h"
}
#include "stdio.h"
#include "string.h"
#include "geotiffio.h"
#include "xtiffio.h"

#ifndef __IMGUTIL__
#define __IMGUTIL__

typedef struct _IplImage
{
    int  width;             /* Image width in pixels.                           */
    int  height;            /* Image height in pixels.                          */
	int  nChannels;         /* Most functions support 1,2,3 or 4 channels */
    int  imageSize;         /* Image data size in bytes
                               (==image->height*image->widthStep
                               in case of interleaved data)*/
    unsigned char *imageData;        /* Pointer to aligned image data.         */
}
IplImage;

typedef struct _IuScalar
{
	double val[4];
}
IuScalar;

/*! Load jpg image from  file path */
IplImage* iu_LoadJpgImage(const char* filepath);
/*! Load tiff image from  file path */
IplImage* iu_LoadTifImage(const char* filepath);

/*! Read geotiff image with file path and  pointer to IplImage structure */
void iu_ReadTifImage(const char* filepath, IplImage* pImage);
/*! Save tiff image with file path and  pointer to IplImage structure */
void iu_SaveTifImage(const char* filepath, const IplImage* pImage);
/*! Save geotiff image with file path, pointer to IplImage structure, tiepoint[(0,0,0)->(x,y,z)], pixscale(xscale,yscale,zscale) and zone of utm */
void iu_SaveGTifImage(const char* filepath, const IplImage* pImage,
	const double tiepoints[6], const double pixscale[3], int utm_zone);

/*! Create a IplImage structure with width, height and nChanels */
IplImage* iu_CreateImage(int width, int height, int nChanels);

/*! Set a pixel data to image(row, column) */
void iu_Set2D(IplImage* pImage, int row, int column, IuScalar value);
/*! get a pixel data from image(row, column) */
IuScalar iu_Get2D(const IplImage* pImage, int row, int column);


IplImage* iu_LoadJpgImage(const char* filepath)
{
	IplImage* pImage = new IplImage;

	jpeg_decompress_struct cinfo;  
	jpeg_error_mgr jerr;  

	// STEP 1: StdError  
	cinfo.err = jpeg_std_error(&jerr);  

	// STEP 2: Create  
	jpeg_create_decompress(&cinfo);

	FILE* pf = fopen(filepath, "rb");

	if (pf != NULL)
	{
		// STEP 3: IO  
		jpeg_stdio_src(&cinfo, pf);

		// STEP 4: Header  
		jpeg_read_header(&cinfo, TRUE);

		pImage->width = cinfo.image_width;
		pImage->height = cinfo.image_height;
		pImage->nChannels = cinfo.num_components;
		pImage->imageSize = pImage->width * pImage->height * pImage->nChannels;
		pImage->imageData = new unsigned char [pImage->imageSize];

		if (pImage->imageData != NULL)
		{
			printf("OK.\nPrepare to decompress the image...\n");

			// STEP 5: Start
			jpeg_start_decompress(&cinfo);
			JSAMPROW row_pointer[1];

			// STEP 6: ReadScan
			while (cinfo.output_scanline < cinfo.output_height)  
			{  
				row_pointer[0] = &pImage->imageData[
					cinfo.output_scanline
					* cinfo.image_width
					* cinfo.num_components
				];
				jpeg_read_scanlines(&cinfo, row_pointer, 1);  
			}

			// STEP 7: Finish
			jpeg_finish_decompress(&cinfo);
		}

		fclose(pf);
	}

	return pImage;
}

IplImage* iu_LoadTifImage(const char* filepath)
{
	IplImage* pImage = new IplImage;

	// STEP 1: Open
	TIFF* tif = TIFFOpen(filepath, "r");

	if (tif)
	{
		uint16 ncn;
		// STEP 2: Get tiff key
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &pImage->width);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &pImage->height);
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &ncn);
		pImage->nChannels = ncn;
		pImage->imageSize = pImage->width * pImage->height * pImage->nChannels;

		// STEP 3: Get data
		pImage->imageData = new unsigned char[pImage->imageSize];
		for (int i=0; i<pImage->height; i++)
		{
			TIFFReadScanline(tif, &(pImage->imageData[i*pImage->width*pImage->nChannels]), i, 0);
		}


		// STEP 4: Close
		TIFFClose(tif);
	}

	return pImage;
}

void iu_ReadTifImage(const char* filepath, IplImage* pImage)
{
	// STEP 1: Open
	TIFF* tif = TIFFOpen(filepath, "r");

	if (tif)
	{
		uint16 ncn;
		// STEP 2: Get tiff key
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &pImage->width);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &pImage->height);
		TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &ncn);
		pImage->nChannels = ncn;
		pImage->imageSize = pImage->width * pImage->height * pImage->nChannels;

		// STEP 3: Get data
		pImage->imageData = new unsigned char[pImage->imageSize];
		for (int i=0; i<pImage->height; i++)
		{
			TIFFReadScanline(tif, &(pImage->imageData[i*pImage->width*pImage->nChannels]), i, 0);
		}


		// STEP 4: Close
		TIFFClose(tif);
	}
}

void iu_SaveTifImage(const char* filepath, const IplImage* pImage)
{
	// STEP 1: Creat and open
	TIFF* tif = TIFFOpen(filepath, "w");
	if (tif)
	{
		// STEP 2: Set title
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, pImage->width);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, pImage->height);
		TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, pImage->nChannels);

		// STEP 3: Set data
		unsigned char* bits = pImage->imageData;
		for (int i=0; i<pImage->height; i++)
		{
			TIFFWriteScanline(tif, &(bits[i*pImage->width*pImage->nChannels]), i, 0);
		}

		// STEP 4: Close and save
		TIFFClose(tif);
	}
}

void iu_SaveGTifImage(const char* filepath, const IplImage* pImage,
	const double tiepoints[6], const double pixscale[3], int utm_zone)
{
	// STEP 1: Creat and open
	TIFF* tif = XTIFFOpen(filepath, "w");
	if (!tif)
		return;
	GTIF* gtif = GTIFNew(tif);
	if (gtif)
	{
		// STEP 2: Set tiff key
		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, pImage->width);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, pImage->height);
		TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
		TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, pImage->nChannels);

		TIFFSetField(tif,TIFFTAG_GEOTIEPOINTS, 6, tiepoints);
		TIFFSetField(tif,TIFFTAG_GEOPIXELSCALE, 3, pixscale);
		
		// STEP 3: Set data
		unsigned char* bits = pImage->imageData;
		for (int i=0; i<pImage->height; i++)
		{
			TIFFWriteScanline(tif, &(bits[i*pImage->width*pImage->nChannels]), i, 0);
		}

		// STEP 4: Set geotiff key
		GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelProjected);
		GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
		GTIFKeySet(gtif, GTCitationGeoKey, TYPE_ASCII, 0, "My GeoTiff");
		
		GTIFKeySet(gtif, GeogCitationGeoKey, TYPE_ASCII, 0, "GCS_WGS_84 used.");
		GTIFKeySet(gtif, GeogLinearUnitsGeoKey, TYPE_SHORT, 1, Linear_Meter);
		GTIFKeySet(gtif, GeogAngularUnitsGeoKey, TYPE_SHORT, 1, Angular_Degree);
		GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);
		GTIFKeySet(gtif, GeogEllipsoidGeoKey, TYPE_SHORT, 1, Ellipse_WGS_84);
		GTIFKeySet(gtif, GeogPrimeMeridianGeoKey, TYPE_SHORT, 1, PM_Greenwich);

		GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, 32600 + utm_zone); //PCS_WGS84_UTM_zone__N = 32600 + utm_zone
		GTIFKeySet(gtif, ProjCoordTransGeoKey, TYPE_SHORT, 1, CT_TransverseMercator);
		GTIFKeySet(gtif, ProjFalseEastingGeoKey, TYPE_DOUBLE, 1, 500000.0);
		GTIFKeySet(gtif, ProjOriginLongGeoKey, TYPE_DOUBLE, 1, double(-177+(utm_zone-1)*6));
		GTIFKeySet(gtif, ProjScaleAtOriginGeoKey, TYPE_DOUBLE, 1, 0.9996);

		// STEP 4: Close and save
		GTIFWriteKeys(gtif);
		GTIFFree(gtif);
		XTIFFClose(tif);
	}
}

IplImage* iu_CreateImage(int width, int height, int nChanels)
{
	IplImage* pImage = new IplImage;
	pImage->width = width;
	pImage->height = height;
	pImage->nChannels = nChanels;
	pImage->imageSize = width * height * nChanels;
	pImage->imageData = new unsigned char[pImage->imageSize];
	memset(pImage->imageData, 0, sizeof(unsigned char) * pImage->imageSize);
	return pImage;
}

void iu_Set2D(IplImage* pImage, int row, int column, IuScalar value)
{
	int nChannels = pImage->nChannels;
	if (row >= pImage->height-1 || column >= pImage->width-1)
		return;
 	for (int n=0; n<nChannels; n++)
	{
		pImage->imageData[row*pImage->width*nChannels + nChannels*column + n] = (unsigned char)value.val[n];
	}
}

IuScalar iu_Get2D(const IplImage* pImage, int row, int column)
{
	int nChannels = pImage->nChannels;
	IuScalar value;
	for (int n=0; n<nChannels; n++)
	{
		value.val[n] = 0;
	}
	if (row >= pImage->height || column >= pImage->width)
		return value;
	for (int n=0; n<nChannels; n++)
	{
		value.val[n] = pImage->imageData[row*pImage->width*nChannels + nChannels*column + n];
	}
	return value;
}

#endif