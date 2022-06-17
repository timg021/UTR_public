//Header XA_tiff.h
//
//	CSIRO MIT
//
//	HEADER FILE TITLE:
//
//		TIFF file I/O facilities
//
//	COPYRIGHT:
//
//		XXX
//					All Rights Reserved
//
//
/*!
	\file		XA_tiff.h
	\brief		TIFF file I/O facilities
	\par		Description:
		This header contains functions that provide simple TIFF file I/O for XArray2D<T> objects.
		Only 8-bit, 16-bit and 32-bit uncompressed grey-scale TIFF I/O is supported.
*/
// The X-Y order here has been adjusted according to the XArray scheme

#if !defined XA_TIFF_H
#define XA_TIFF_H

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include <limits>
#include <omp.h>

#include "XArray2D.h"
#include "XArray3D.h"

//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
struct IXArray2D;

//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//

namespace // anonymous
{
	typedef unsigned char XARBYTE;
	typedef unsigned short XARWORD;
	typedef unsigned long XARDWORD;


	// TIFF datatypes
	const XARWORD TDT_BYTE		= 1; // byte (8-bit)
	const XARWORD TDT_ASCII	= 2; // zero terminated string (8-bit chars)
	const XARWORD TDT_SHORT	= 3; // short (16-bit)
	const XARWORD TDT_LONG		= 4; // long (32-bit)
	const XARWORD TDT_RATIONAL	= 5; // two longs representing the numerator and the denominator (64-bit)


	// TIFF tags (standard) sorted in ascending order
	const XARWORD TT_NSUBFILETYPE	= 0xFE;  // new subfile type
	const XARWORD TT_IMAGEWIDTH	= 0x100;
	const XARWORD TT_IMAGELENGTH	= 0x101;
	const XARWORD TT_BITSPERSAMPLE	= 0x102;
	const XARWORD TT_COMPRESSION	= 0x103;
	const XARWORD TT_PMETRICINTERP	= 0x106; // photometric interpretation
	const XARWORD TT_STRIPOFFSETS	= 0x111;
	const XARWORD TT_SAMPLESPERPIX	= 0x115;
	const XARWORD TT_ROWSPERSTRIP	= 0x116;
	const XARWORD TT_STRIPBYTECNTS	= 0x117; // strip byte counts
	const XARWORD TT_PLANARCONFIG	= 0x11C;
	const XARWORD TT_RESUNIT		= 0x128; // resolution unit
	const XARWORD TT_COLORMAP		= 0x140;


	// custom tags NOTE: these are temporary values
	const XARWORD TCT_XLO			= 32769;
	const XARWORD TCT_XHI			= 32770;
	const XARWORD TCT_YLO			= 32771;
	const XARWORD TCT_YHI			= 32772;
	const XARWORD TCT_ZLO			= 32773;
	const XARWORD TCT_ZHI			= 32774;

//---------------------------------------------------------------------------
//	MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
	typedef struct TIFFINFO_tag
	{
		XARWORD wBitsPerSample;
		XARWORD wCompression;
		XARDWORD dwImageLength;
		XARDWORD dwImageWidth;
		XARDWORD dwNewSubfileType;
		XARWORD wPhotoInterp; // photometric interpretation
		XARDWORD dwRowsPerStrip;
		XARWORD wSamplesPerPixel;
		XARDWORD *pdwStripOffsets;
		double dLoX;
		double dHiX;
		double dLoY;
		double dHiY;
		double dLoZ;
		double dHiZ;
	}TIFFINFO, *LPTIFFINFO;


	typedef struct TAGNAME_tag
	{
		char szTagName[50];
		XARWORD wTagId;
	}TAGNAME, *LPTAGNAME;


	const TAGNAME g_tnTagNames[] = 
	{ 
		{"NewSubfileType",				0xFE},
		{"SubfileType",					0xFF},
		{"ImageWidth",					0x100},
		{"ImageLength",					0x101},
		{"BitsPerSample",				0x102},
		{"Compression",					0x103},
		{"PhotometricInterpretation",	0x106},		
		{"Threshholding",				0x107},
		{"CellWidth",					0x108},
		{"CellLength",					0x109},
		{"FillOrder",					0x10A},
		{"DocumentName",				0x10D},
		{"ImageDescription",			0x10E},
		{"Make",						0x10F},
		{"Model",						0x110},		
		{"StripOffsets",				0x111},
		{"Orientation",					0x112},
		{"SamplesPerPixel",				0x115},	
		{"RowsPerStrip",				0x116},	
		{"StripByteCounts",				0x117},
		{"MinSampleValue",				0x118},
		{"MaxSampleValue",				0x119},
		{"XResolution",					0x11A},
		{"YResolution",					0x11B},
		{"PlanarConfiguration",			0x11C},
		{"PageName",					0x11D},
		{"XPosition",					0x11E},
		{"YPosition",					0x11F},
		{"FreeOffsets",					0x120},
		{"FreeByteCounts",				0x121},
		{"GrayResponseUnit",			0x122},
		{"GrayResponseCurve",			0x123},
		{"Group3Options",				0x124},
		{"Group4Options",				0x125},
		{"ResolutionUnit",				0x128},
		{"PageNumber",					0x129},
		{"ColorResponseCurves",			0x12D},
		{"Software",					0x131},
		{"DateTime",					0x132},	
		{"Artist",						0x13B},
		{"HostComputer",				0x13C},
		{"Predictor",					0x13D},
		{"WhitePoint",					0x13E},
		{"PrimaryChromaticities",		0x13F},
		{"ColorMap",					0x140},
		{"X-Low",						TCT_XLO},
		{"X-High",						TCT_XHI},
		{"Y-Low",						TCT_YLO},
		{"Y-High",						TCT_YHI},
		{"Z-Low",						TCT_ZLO},
		{"Z-High",						TCT_ZHI}
	};


	const XARWORD g_wTags = sizeof(g_tnTagNames)/sizeof(TAGNAME);


	typedef struct TIFFOPTIONS_tag
	{
		XARWORD wBitsPerSample; // number of bits per sample
		bool fWriteCustomTags; // TRUE to write custom tags
	}TIFFOPTIONS, *LPTIFFOPTIONS;

//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//

	//! Writes a tiff tag to a specified file
	//
	// PARAMS:
	// - pFile:		pointer to an open output file with appropriately positioned pointer
	// - wTag:		tag id
	// - wType:		tag type id
	// - dwByte:	length of the value in _bytes_
	// - dwOffset:  offset to where the value will be written
	// - pValue;	pointer to value
	// - pdwDeltaOffset; (pointer to) change of offset 
	inline void WriteTiffTag(FILE *pFile, XARWORD wTag, XARWORD wType, XARDWORD dwByteLen, XARDWORD dwOffset, void *pValue, XARDWORD* pdwDeltaOffset, const char* pchFilename)
	{
		XARDWORD dwLen; // length
		long lCurPos; // current file position
		char szEnding[] = {0, 0, 0};

		*pdwDeltaOffset = 0; // change of offset

		if(!fwrite(&wTag, sizeof(wTag), 1, pFile)) // write tag id
			throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");

		if(!fwrite(&wType, sizeof(wType), 1, pFile)) // write the data type
			throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");

		switch(wType)
		{
		case TDT_BYTE:
		case TDT_ASCII:
			dwLen = dwByteLen; // 8-bit
			break;

		case TDT_SHORT:
			dwLen = dwByteLen/2; // 16-bit
			break;

		case TDT_LONG: // 32-bit
			dwLen = dwByteLen/4;
			break;

		case TDT_RATIONAL: // 64-bit
			dwLen = dwByteLen/8;
			break;

		default:
			throw std::runtime_error("runtime_error in WriteTiffTag (unknown data type in wType)");
		}

		if(!fwrite(&dwLen, sizeof(dwLen), 1, pFile)) // write the length
			throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");

		if(dwByteLen <= 4) // we should write offset as the value
		{
			if(!fwrite(pValue, dwByteLen, 1, pFile)) // write the value
				throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");

			if(dwByteLen < 4)
				if(!fwrite(szEnding, 4 - dwByteLen, 1, pFile)) // write the value
					throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");

			*pdwDeltaOffset = 0;
			return;
		}

		if(!(dwOffset % 2)) // value must begin at a XARWORD boundary
		{
			dwOffset++;
			(*pdwDeltaOffset)++;
		}

		if(!fwrite(&dwOffset, sizeof(dwOffset), 1, pFile)) // write the offset
			throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");
		*pdwDeltaOffset += 4;
		
		lCurPos = ftell(pFile);
		fseek(pFile, dwOffset, SEEK_SET);
		if(!fwrite(pValue, dwByteLen, 1, pFile)) // write the value
			throw std::runtime_error("runtime_error in WriteTiffTag (error writing to output file)");
		fseek(pFile, lCurPos, SEEK_SET);

		return;
	}

	//! Reads a TIFF tag
	//
	// PARAMS:
	// - pFile:			pointer to the input file with with file pointer pointing to 
	//					the appropriate tag
	// - pwTagId:		pointer to rXAr2D XARWORD that will receive id of the tag
	// - pwType:		pointer to a XARWORD that will receive the tag type
	// - ppValue:		pointer to a pointer that will point to anewly allocated block of 
	//					memory containing the value if it doesn't fit into a XARDWORD
	inline void ReadTiffTag(FILE *pFile, XARWORD *pwTagId, XARWORD *pwType, XARDWORD *pdwLen, void **ppValue, const char* pchFilename)
	{
		XARDWORD dwLen;
		XARDWORD dwOffset;
		long lCurPos; // current position in the file
		
		if(!fread(pwTagId, 2, 1, pFile)) // read in the tag id
			throw std::runtime_error("runtime_error in ReadTiffTag (error reading input file)");

		if(!fread(pwType, 2, 1, pFile)) // read in the datatype
			throw std::runtime_error("runtime_error in ReadTiffTag (error reading input file)");

		if(!fread(&dwLen, 4, 1, pFile)) // read in the length
			throw std::runtime_error("runtime_error in ReadTiffTag (error reading input file)");
			
		*pdwLen = dwLen;

		if(!fread(&dwOffset, 4, 1, pFile)) // read in the offset
			throw std::runtime_error("runtime_error in ReadTiffTag (error reading input file)");

		// calculate the byte length
		switch(*pwType)
		{
		case TDT_BYTE:
		case TDT_ASCII:
			break;

		case TDT_SHORT:
			dwLen *= 2;
			break;

		case TDT_LONG:
			dwLen *= 4;
			break;

		case TDT_RATIONAL:
			dwLen *= 8;
			break;

		default:
			throw std::runtime_error("runtime_error in ReadTiffTag (unknown data type)");
		}

		// alloc memory for the value
		*ppValue = new XARBYTE[dwLen];
		if(!*ppValue)
			throw std::runtime_error("runtime_error in ReadTiffTag (out of memory)");

		if(dwLen <= 4) // if the value fits into XARDWORD, offset contains the actual value
		{
			memcpy(*ppValue, &dwOffset, dwLen);
			return;
		}
		else // the value is stored at some other location in the file
		{
			lCurPos = ftell(pFile); // remeber the current location
			fseek(pFile, dwOffset, SEEK_SET);
			if(!fread(*ppValue, dwLen, 1, pFile))
				throw std::runtime_error("runtime_error in ReadTiffTag (error reading input file)");
		}
		
		fseek(pFile, lCurPos, SEEK_SET); // go back to the saved location

		return;
	}

	//! Reads all of the tags that we care about
	//
	// PARAMS:
	// - pFile:			pointer to the input file with with file pointer pointing the start of IFD
	// - ptiInfo:		pointer to the TIFFINFO strucutre that will receive the tag values
	// - ptpParams		pointer to the T2GPARAMS structure containing params
	//
	inline void ReadTiffTags(FILE *pFile, TIFFINFO *ptiInfo, const char* pchFilename)
	{
		XARWORD wTags;
		XARWORD wCurTag;
		XARWORD wTagId;
		XARWORD wType;
		void *pValue0;
		XARBYTE* pValue;
		XARDWORD dwLen;
		XARDWORD dwNextIFD;

		if(!fread(&wTags, 2, 1, pFile))
			throw std::runtime_error("runtime_error in ReadTiffTags (error reading input file)");


		for(wCurTag = 0; wCurTag < wTags; wCurTag++)
		{
			ReadTiffTag(pFile, &wTagId, &wType, &dwLen, &pValue0, pchFilename);
			pValue = (XARBYTE*)pValue0;

			// process the documented tags
			switch(wTagId)
			{
			case TT_BITSPERSAMPLE: // must be short
				if(wType != TDT_SHORT)
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				ptiInfo->wBitsPerSample = *((XARWORD *)pValue);
				delete [] pValue;
				break;

			case TT_COMPRESSION: // must be short
				if(wType != TDT_SHORT)
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				ptiInfo->wCompression = *((XARWORD *)pValue);
				delete [] pValue;
				break;

			case TT_IMAGELENGTH: // must be short or long
				if(wType == TDT_SHORT)
					ptiInfo->dwImageLength = *((XARWORD *)pValue);	
				else if(wType == TDT_LONG)
					ptiInfo->dwImageLength = *((XARDWORD *)pValue);	
				else
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				delete [] pValue;
				break;

			case TT_IMAGEWIDTH: // must be short or long
				if(wType == TDT_SHORT)
					ptiInfo->dwImageWidth = *((XARWORD *)pValue);	
				else if(wType == TDT_LONG)
					ptiInfo->dwImageWidth = *((XARDWORD *)pValue);	
				else
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				delete [] pValue;
				break;

			case TT_NSUBFILETYPE: // must be long
				if(wType != TDT_LONG)
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				ptiInfo->dwNewSubfileType = *((XARDWORD *)pValue);
				delete [] pValue;
				break;
			
			case TT_PMETRICINTERP: // must be short
				if(wType != TDT_SHORT)
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				ptiInfo->wPhotoInterp = *((XARWORD *)pValue);
				delete [] pValue;
				break;

			case TT_ROWSPERSTRIP: // must be short or long
				if(wType == TDT_SHORT)
					ptiInfo->dwRowsPerStrip = *((XARWORD *)pValue);	
				else if(wType == TDT_LONG)
					ptiInfo->dwRowsPerStrip = *((XARDWORD *)pValue);	
				else
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				delete [] pValue;
				break;

			case TT_STRIPBYTECNTS:
				delete [] pValue;
				break;

			case TT_SAMPLESPERPIX: // must be short
				if(wType != TDT_SHORT)
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				ptiInfo->wSamplesPerPixel = *((XARWORD *)pValue);
				delete [] pValue;
				break;

			case TT_STRIPOFFSETS: // must be long
				if(wType != TDT_LONG)
					throw std::runtime_error("runtime_error in ReadTiffTags (invalid data type for tag)");
				ptiInfo->pdwStripOffsets = (XARDWORD *) pValue;	// we store the pointer but don't delete it
				break;

			case TT_RESUNIT: // must be short
				delete [] pValue;
				break;

			default:
				delete [] pValue;
				break;
			}
		}

		if(!fread(&dwNextIFD, 4, 1, pFile))// get the offset to the next IFD (there should be only one)
	//		throw std::runtime_error("runtime_error in ReadTiffTags (error reading input file)");
		//if(dwNextIFD) // another IFD exists, too bad, we'll ignore it
		//	MDWarning("ReadTIFFile(): TIFF file contains more than one image. Only the first image will be processed");
		1 == 1; //we'll ignore this
	}

	//! Verifies that the values are valid and are supported
	//
	// PARAMS:
	// - pFile:			pointer to the TIFFINFO strucutre tcontaining the tag values
	//
	inline void VerifyTiffTags(TIFFINFO *ptiInfo, const char* pchFilename)
	{
		if(ptiInfo->wBitsPerSample % 8 || ptiInfo->wBitsPerSample > 32) // only 8, 16 and 32 bits per sample TIFFs are supported
			throw std::runtime_error("runtime_error in VerifyTiffTags (only 8, 16 and 32 bits per sample TIFFs are supported)");

		if(ptiInfo->wCompression != 1) // we can't handle compression
			throw std::runtime_error("runtime_error in VerifyTiffTags (compressed TIFFs are not supported)");

		if(ptiInfo->dwImageLength <= 0 || ptiInfo->dwImageWidth <= 0) // image length and width must be non-zero
			throw std::runtime_error("runtime_error in VerifyTiffTags (image length and width must both be positive)");

		if(ptiInfo->dwNewSubfileType) // check the type of the current image. Exit() if it's a thumbnail or a transparency mask
			throw std::runtime_error("runtime_error in VerifyTiffTags (first subfile is not the full image)");

		//if(ptiInfo->wPhotoInterp > 1) // 0 and 1 are the values for grayscale images
		//	throw std::runtime_error("runtime_error in VerifyTiffTags (only grayscale images are supported)");

		if(ptiInfo->wSamplesPerPixel > 1) // grayscale images should have only one sample per pixel
			throw std::runtime_error("runtime_error in VerifyTiffTags (only grayscale images are supported)");

		if(!ptiInfo->pdwStripOffsets)
			throw std::runtime_error("runtime_error in VerifyTiffTags (the file is corrupt)");
	}

	//! Reads the tags of a 8, 16 or 32-bit gray-scale uncompressed TIFF file
	//
	// PARAMS:
	//pchFilename - full name of the file to be read
	//
	//NOTE: this program is based on the TIFF2GRD program by D.Ternovski
	inline TIFFINFO ReadTIFFInfo(const char* pchFilename)
	{

		XARDWORD dwTemp = 0;
		XARWORD wTemp = 0;
		TIFFINFO tiInfo;

		// init tiInfo with default values
		tiInfo.wBitsPerSample = 1;
		tiInfo.wCompression = 1;
		tiInfo.dwImageLength = 0;
		tiInfo.dwImageWidth = 0;
		tiInfo.dwNewSubfileType = 0;
		tiInfo.wPhotoInterp = 1;
		tiInfo.dwRowsPerStrip = 0xFFFFFFFF;
		tiInfo.wSamplesPerPixel = 1;
		tiInfo.pdwStripOffsets = NULL;

		xar::FilePtr pInputFile(pchFilename, "rb");
			throw std::runtime_error("runtime_error in ReadTIFFInfo (error opening input file)");

		if(!fread(&wTemp, 2, 1, pInputFile)) // read in the byte ordering XARWORD
			throw std::runtime_error("runtime_error in ReadTIFFInfo (error reading input file)");
		if(wTemp != 0x4949)
			throw std::runtime_error("runtime_error in ReadTIFFInfo (only Intel byte ordering is supported)");

		if(!fread(&wTemp, 2, 1, pInputFile)) // read in the tiff version
			throw std::runtime_error("runtime_error in ReadTIFFInfo (error reading input file)");
		if(wTemp != 42)
			throw std::runtime_error("runtime_error in ReadTIFFInfo (invalid TIFF version)");

		if(!fread(&dwTemp, 4, 1, pInputFile)) // read in the offset to the IFD
			throw std::runtime_error("runtime_error in ReadTIFFInfo (error reading input file)");

		fseek(pInputFile, dwTemp, SEEK_SET);

		ReadTiffTags(pInputFile, &tiInfo, pchFilename);
		VerifyTiffTags(&tiInfo, pchFilename);

		delete [] tiInfo.pdwStripOffsets; tiInfo.pdwStripOffsets = 0;

		return tiInfo;
	}

} // namespace anonymous closed


namespace xar
{
	//---------------------------------------------------------------------------
	//Function TIFFWriteFile
	//
	//	Writes XArray2D<T> object into an uncompressed grey-scale TIFF file
	//
	/*!
		\brief		Writes XArray2D<T> object into an uncompressed grey-scale TIFF file
		\param		rXAr2D	Reference to an XArray2D<T> object to be saved
		\param		pchFilename	Full name of the file to be written into
		\param		filetype eTIFF8, eTIFF16 or eTIFF32 for 8, 16 or 32 bit TIFF files respectively
		\param		bOptimize Determines if rXAr2D's data is stretched to the full range of the TIFF file (ignored for floating-point data)
		\exception  std::invalid_argument is thrown if XArray2D object is empty or contains complex values
		\exception  std::invalid_argument is thrown if 'filetype' parameter is not eTIFF8, eTIFF16 or eTIFF32
		\exception  std::runtime_error is thrown if any of the file write operations fail
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function writes a real 2D XArray object into an uncomressed 8, 16 or 32 bit
			TIFF file. The saved data can be 'optimized', i.e. the range of values contained
			in the XArray2D object can be linearly mapped onto the maximum data range available
			for the given TIFF format. The head of the object is not saved
	*/	
	// NOTE: this function is based on the GRD2TIFF program by D.Ternovski
	// NOTE: for floating point data (filetype == xar::eTIFF32), the bOptimize parameter is ignored
	template <class T> void TIFFWriteFile(const XArray2D<T>& rXAr2D, const char* pchFilename, _eFileType filetype, bool bOptimize = true)
	{
		if (rXAr2D.GetValuetype() == eXAFComplex || rXAr2D.GetValuetype() == eXADComplex)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in TIFFWriteFile (complex data)");

		if (rXAr2D.size() == 0)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in TIFFWriteFile (no data)");

		XARBYTE /*cTemp,*/ *pcTemp;
		XARWORD wTemp, *pwTemp;
		XARDWORD dwTemp/*, *pdwTemp*/;
		float* pfTemp;
		XARDWORD dwOffset; // current offset
		XARDWORD dwDeltaOffset = 0; //delta offset
		XARDWORD* pdwDeltaOffset = &dwDeltaOffset;
		index_t i,j;
		double umin, umax, tmax;
		double dPropConst(1.0); // proportionality constant - (2^BitsPerSample - 1 )/ (umax - umin)
		XARDWORD dwTags = 7; // number of tags
		short BytesPerSample; // = wBitsPerSample / 8;
		
		switch (filetype)
		{
		case eTIFF8: BytesPerSample = 1; break;
		case eTIFF16: BytesPerSample = 2; break;
		case eTIFF32: BytesPerSample = 4; break;
		default: throw std::invalid_argument("invalid_argument 'filetype' in TIFFWriteFile");
		}

		FilePtr pOutput(pchFilename, "wb");

		tmax = pow(2, BytesPerSample*8) - 1;
		if (bOptimize && BytesPerSample != 4)
		{
			umin = rXAr2D.Norm(eNormMin);
			umax = rXAr2D.Norm(eNormMax);
			dPropConst = (umax-umin) ?  tmax / (umax - umin) : 0;
		}
		
		if(!fwrite("II", 2, 1, pOutput)) // Intel bigendian format
			throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");

		wTemp = 42; // version number
		if(!fwrite(&wTemp, 2, 1, pOutput)) // write the version number
			throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");

		dwOffset = 8 + (XARDWORD)rXAr2D.size() * BytesPerSample; // IFD goes straight after the strip
		dwOffset += dwOffset % 2; // IFD must begin on a word boundary
		if(!fwrite(&dwOffset, 4, 1, pOutput)) // write the offset to IFD
			throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");

		// write the strip
		switch (BytesPerSample)
		{
		case 1:
			pcTemp = new XARBYTE [rXAr2D.GetDim2()];
			if (bOptimize)
				for(i = 0; i < rXAr2D.GetDim1(); i++)				
				{
					for(j = 0; j < rXAr2D.GetDim2(); j++)
					{
						pcTemp[j] = (XARBYTE)(dPropConst * (rXAr2D[i][j] - umin) + 0.5);
					}
					if (fwrite(pcTemp, BytesPerSample, rXAr2D.GetDim2(), pOutput) != rXAr2D.GetDim2())
					{
						delete pcTemp;
						throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
					}
				}
			else
				for(i = 0; i < rXAr2D.GetDim1(); i++)				
				{
					for(j = 0; j < rXAr2D.GetDim2(); j++)
					{
						if (rXAr2D[i][j] < 0) pcTemp[j] = 0;
						else if (rXAr2D[i][j] > 255) pcTemp[j] = 255;
						else pcTemp[j] = (XARBYTE)rXAr2D[i][j];
					}
					if (fwrite(pcTemp, BytesPerSample, rXAr2D.GetDim2(), pOutput) != rXAr2D.GetDim2())
					{
						delete pcTemp;
						throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
					}
				}
			delete pcTemp;
			break;
		case 2:
			pwTemp = new XARWORD [rXAr2D.GetDim2()];
			if (bOptimize)
				for(i = 0; i < rXAr2D.GetDim1(); i++)				
				{
					for(j = 0; j < rXAr2D.GetDim2(); j++)
					{
						pwTemp[j] = (XARWORD)(dPropConst * (rXAr2D[i][j] - umin) + 0.5);
					}
					if (fwrite(pwTemp, BytesPerSample, rXAr2D.GetDim2(), pOutput) != rXAr2D.GetDim2())
					{
						delete pwTemp;
						throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
					}
				}
			else
				for(i = 0; i < rXAr2D.GetDim1(); i++)				
				{
					for(j = 0; j < rXAr2D.GetDim2(); j++)
					{
						if (rXAr2D[i][j] < 0) pwTemp[j] = 0;
						else if (rXAr2D[i][j] > 65535) pwTemp[j] = 65535;
						else pwTemp[j] = (XARWORD)rXAr2D[i][j];
					}
					if (fwrite(pwTemp, BytesPerSample, rXAr2D.GetDim2(), pOutput) != rXAr2D.GetDim2())
					{
						delete pwTemp;
						throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
					}
				}
			delete pwTemp;
			break;
		case 4:
			pfTemp = new float [rXAr2D.GetDim2()];

			for(i = 0; i < rXAr2D.GetDim1(); i++)				
			{
				for(j = 0; j < rXAr2D.GetDim2(); j++)
				{
					pfTemp[j] = (float)(rXAr2D[i][j]);
				}
				if (fwrite(pfTemp, BytesPerSample, rXAr2D.GetDim2(), pOutput) != rXAr2D.GetDim2())
				{
					delete pfTemp;
					throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
				}
			}
			delete pfTemp;
			break;
		}
		
		// write the tags
		fseek(pOutput, dwOffset, SEEK_SET);
		dwOffset += 2 + dwTags*12; // offset to the start of values

		if(!fwrite(&dwTags, 2, 1, pOutput)) // write the number of tags
			throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
		
		// image width
		dwTemp = (XARDWORD)rXAr2D.GetDim2();
		WriteTiffTag(pOutput, TT_IMAGEWIDTH, TDT_LONG, 4, dwOffset, &dwTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		// image length
		dwTemp = (XARDWORD)rXAr2D.GetDim1();
		WriteTiffTag(pOutput, TT_IMAGELENGTH, TDT_LONG, 4, dwOffset, &dwTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		// Bits per sample
		wTemp = BytesPerSample * 8;
		WriteTiffTag(pOutput, TT_BITSPERSAMPLE, TDT_SHORT, 2, dwOffset, &wTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		// photometric interpretation
		wTemp = 1; // grayscale
		WriteTiffTag(pOutput, TT_PMETRICINTERP, TDT_SHORT, 2, dwOffset, &wTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		// strip offsets
		dwTemp = 8; // the strip goes straight after the image file header
		WriteTiffTag(pOutput, TT_STRIPOFFSETS, TDT_LONG, 4, dwOffset, &dwTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		// strip byte counts
		dwTemp = (XARDWORD)rXAr2D.size();
		WriteTiffTag(pOutput, TT_STRIPBYTECNTS, TDT_LONG, 4, dwOffset, &dwTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		// data type
		if (BytesPerSample == 4) wTemp = 3;  // floating point
		else wTemp = 1;  // unsigned integer data
		WriteTiffTag(pOutput, 0x153, TDT_SHORT, 2, dwOffset, &wTemp, &dwDeltaOffset, pchFilename);
		dwOffset += dwDeltaOffset;

		dwTemp = 0; // offset to the next IFD (there's only one)
		if(!fwrite(&dwTemp, 4, 1, pOutput))
			throw std::runtime_error("runtime_error in TIFFWriteFile (error writing to output file)");
	}


	//---------------------------------------------------------------------------
	//Function TIFFReadFile
	//
	//	Reads XArray2D<T> object from an uncompressed grey-scale TIFF file
	//
	/*!
		\brief		Reads XArray2D<T> object from an uncompressed grey-scale TIFF file
		\param		rXAr2D	Reference to an XArray2D<T> object to be read in
		\param		pchFilename	Full name of the file to be read from
		\exception  std::invalid_argument is thrown if XArray2D object is of complex type
		\exception  std::runtime_error is thrown if any of the file read operations fail
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function reads a real 2D XArray object from an uncomressed grey-scale 8, 16 or 32 bit
			TIFF file. The head of the object is not affected
	*/	
	// NOTE: this program is based on the TIFF2GRD program by D.Ternovski
	template <class T> void TIFFReadFile(XArray2D<T>& rXAr2D, const char* pchFilename)
	{
		if (rXAr2D.GetValuetype() == eXAFComplex || rXAr2D.GetValuetype() == eXADComplex)
			throw std::invalid_argument("invalid_argument 'rXAr2D' in TIFFReadFile (complex XArray)");

		XARBYTE *pcTemp;
		XARWORD wTemp = 0, *pwTemp;
		XARDWORD dwTemp = 0/*, *pdwTemp*/;
		float* pfTemp;
		TIFFINFO tiInfo;
		XARDWORD dwCurStrip;
		XARDWORD dwStrips;
		XARDWORD dwCurRow;
		XARDWORD dwCount, dwRow;

		// init tiInfo with default values
		tiInfo.wBitsPerSample = 1;
		tiInfo.wCompression = 1;
		tiInfo.dwImageLength = 0;
		tiInfo.dwImageWidth = 0;
		tiInfo.dwNewSubfileType = 0;
		tiInfo.wPhotoInterp = 1;
		tiInfo.dwRowsPerStrip = 0xFFFFFFFF;
		tiInfo.wSamplesPerPixel = 1;
		tiInfo.pdwStripOffsets = NULL;

		FilePtr pInputFile(pchFilename, "rb");

		if(!fread(&wTemp, 2, 1, pInputFile)) // read in the byte ordering XARWORD
			throw std::runtime_error("runtime_error in TIFFReadFile (error reading from file)");
		if(wTemp != 0x4949)
			throw std::runtime_error("runtime_error in TIFFReadFile (non-Intel byte ordering is not supported)");

		if(!fread(&wTemp, 2, 1, pInputFile)) // read in the tiff version
			throw std::runtime_error("runtime_error in TIFFReadFile (error reading from file)");
		if(wTemp != 42)
			throw std::runtime_error("runtime_error in TIFFReadFile (unsupported TIFF version)");

		if(!fread(&dwTemp, 4, 1, pInputFile)) // read in the offset to the IFD
			throw std::runtime_error("runtime_error in TIFFReadFile (error reading from file)");

		fseek(pInputFile, dwTemp, SEEK_SET);

		ReadTiffTags(pInputFile, &tiInfo, pchFilename);

		VerifyTiffTags(&tiInfo, pchFilename);
		if (tiInfo.wBitsPerSample > sizeof(T)*8) 
			throw std::runtime_error("runtime_error in TIFFReadFile (type T of XArray2D<T> is too small for the file data)");
		
		if(tiInfo.dwRowsPerStrip > tiInfo.dwImageLength) tiInfo.dwRowsPerStrip = tiInfo.dwImageLength;
		dwStrips = (tiInfo.dwImageLength + tiInfo.dwRowsPerStrip - 1)/tiInfo.dwRowsPerStrip;

		rXAr2D.Resize(tiInfo.dwImageLength, tiInfo.dwImageWidth, 0);		

		//DC shift may be necessary to fit an 'unsigned' type from TIFF file into a 'signed' type T
		//factor may be necessary if TIFF is a negative image
		short BytesPerSample = tiInfo.wBitsPerSample / 8;
		double dC, factor;

		if(tiInfo.wPhotoInterp) // 1: 0 is black
		{
			factor = 1;
			dC = 0.0;
			if (std::numeric_limits<T>::is_integer && BytesPerSample == sizeof(T)) 
				dC = std::numeric_limits<T>::min();
		}
		else // 0: 0 is white
		{
			factor = -1;
			dC = 0.0;
			if (std::numeric_limits<T>::is_integer && BytesPerSample == sizeof(T)) 
				dC = std::numeric_limits<T>::max();
		}


		// process all strips
		for(dwCurStrip = 0; dwCurStrip < dwStrips; dwCurStrip++)
		{
			fseek(pInputFile, tiInfo.pdwStripOffsets[dwCurStrip], SEEK_SET);
			 // and all rows in each strip
			switch (BytesPerSample)
			{
			case 1:
				pcTemp = new XARBYTE[tiInfo.dwImageWidth];
				for(dwCurRow = 0; dwCurRow < tiInfo.dwRowsPerStrip && ((dwCurStrip*tiInfo.dwRowsPerStrip + dwCurRow) < tiInfo.dwImageLength); dwCurRow++)
				{
					if(fread(pcTemp, 1, tiInfo.dwImageWidth, pInputFile) != tiInfo.dwImageWidth)
					{
						delete pcTemp;
						throw std::runtime_error("runtime_error in TIFFReadFile (error reading from file)");
					}
					dwRow = dwCurRow + tiInfo.dwRowsPerStrip * dwCurStrip;
					for(dwCount = 0; dwCount < tiInfo.dwImageWidth; dwCount++) // process all pixels in the row
						rXAr2D[dwRow][dwCount] = T(dC + factor * pcTemp[dwCount]); // convert to T
				}
				delete pcTemp;
				break;
			case 2:
				pwTemp = new XARWORD[tiInfo.dwImageWidth];
				for(dwCurRow = 0; dwCurRow < tiInfo.dwRowsPerStrip && ((dwCurStrip*tiInfo.dwRowsPerStrip + dwCurRow) < tiInfo.dwImageLength); dwCurRow++)
				{
					if(fread(pwTemp, 2, tiInfo.dwImageWidth, pInputFile) != tiInfo.dwImageWidth)
					{
						delete pwTemp;
						throw std::runtime_error("runtime_error in TIFFReadFile (error reading from file)");
					}
					dwRow = dwCurRow + tiInfo.dwRowsPerStrip * dwCurStrip;
					for(dwCount = 0; dwCount < tiInfo.dwImageWidth; dwCount++) // process all pixels in the row
						rXAr2D[dwRow][dwCount] = T(dC + factor * pwTemp[dwCount]); // convert to T
				}
				delete pwTemp;
				break;
			case 4:
				pfTemp = new float[tiInfo.dwImageWidth];
				for(dwCurRow = 0; dwCurRow < tiInfo.dwRowsPerStrip && ((dwCurStrip*tiInfo.dwRowsPerStrip + dwCurRow) < tiInfo.dwImageLength); dwCurRow++)
				{
					if(fread(pfTemp, 4, tiInfo.dwImageWidth, pInputFile) != tiInfo.dwImageWidth)
					{
						delete pfTemp;
						throw std::runtime_error("runtime_error in TIFFReadFile (error reading from file)");
					}
					dwRow = dwCurRow + tiInfo.dwRowsPerStrip * dwCurStrip;
					for(dwCount = 0; dwCount < tiInfo.dwImageWidth; dwCount++) // process all pixels in the row
						rXAr2D[dwRow][dwCount] = T(dC + factor * pfTemp[dwCount]); // convert to T
				}
				delete pfTemp;
				break;
			default:
				throw std::runtime_error("runtime_error in TIFFReadFile (unsupported number of bytes per sample)");
			}
		}
		
		delete [] tiInfo.pdwStripOffsets;
	}


	//! Reports the _eFileType of an 8, 16 or 32-bit grey-scale uncompressed TIFF file
	//
	// PARAMS:
	// pchFilename - full name of the file to be read
	//
	// NOTE: this program is based on the TIFF2GRD program by D.Ternovski
	inline _eFileType TIFFReadFileType(const char* pchFilename)
	{
		switch (ReadTIFFInfo(pchFilename).wBitsPerSample)
		{
		case 8: return xar::eTIFF8;
		case 16: return xar::eTIFF16;
		case 32: return xar::eTIFF32;
		default: return xar::eFTUnknown;
		}
	}

	//---------------------------------------------------------------------------
//Function TIFFWriteFileStack
//
//	Writes real XArray3D object into a stack of uncompressed grey-scale TIFF files
//
/*!
	\brief		Writes an XArray3D object into a stack of uncompressed grey-scale TIFF files
	\param		rXAr3D	Reference to an XArray3D<T> object to be saved
	\param		voutfilenames std::vector of std::strings with full names of the files to be written into
	\param		filetype eTIFF8, eTIFF16 or eTIFF32 for 8, 16 or 32 bit TIFF files respectively
	\param		bOptimize Determines if rXAr3D's data is stretched to the full range of the TIFF file (ignored for floating-point data)
	\exception  std::invalid_argument is thrown if XArray3D object contains complex values
	\exception  std::invalid_argument is thrown if 'filetype' parameter is not eTIFF8, eTIFF16 or eTIFF32
	\exception  std::runtime_error is thrown if any of the file write operations fail
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function writes a real 3D XArray object into a stack of uncompressed grey-scale TIFF files
*/
	template <class T> void TIFFWriteFileStack(XArray3D<T>& rXAr3D, vector<string> voutfilenames, _eFileType filetype, bool bOptimize = true)
	{
		if (rXAr3D.GetValuetype() == eXAFComplex || rXAr3D.GetValuetype() == eXADComplex)
			throw std::invalid_argument("invalid_argument 'rXAr3D' in TIFFWriteFileStack (complex data)");

		index_t nz = rXAr3D.GetDim1();
		if (nz != voutfilenames.size()) throw std::invalid_argument("invalid_argument 'rXAr3D' in XArData::TIFFWriteFilStacke (Dim1 of XArray3D object is different from the length of voutfilenames vector)");
		index_t ny = rXAr3D.GetDim2(); index_t nx = rXAr3D.GetDim3();

		#pragma omp parallel for shared (rXAr3D, voutfilenames, nz, ny, nx)
		for (int n = 0; n < nz; n++)
		{
			XArray2D<T> ipOut(ny, nx);
			for (index_t j = 0; j < ny; j++)
				for (index_t i = 0; i < nx; i++)
					ipOut[j][i] = rXAr3D[n][j][i];

			TIFFWriteFile(ipOut, voutfilenames[n].c_str(), filetype, bOptimize);
		}
	}



} //namespace xar closed


//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//


#endif	// XA_TIFF_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//