#pragma once
//Header XA_DICOM.h
//
//
//	HEADER FILE TITLE:
//
//		DICOM file I/O facilities for XArrays
//
/*!
	\file		XArFileIO_DICOM.h
	\brief		DICOM file I/O facilities for XArrays
	\par		Description:
		This header contains a class that provide DICOM file I/O for XArrayND<T> objects. 
*/
//---------------------------------------------------------------------------
//	INCLUDE FILES
//

//#include "XArFileIO.h"
//#include "XArray/XArray2D.h"
//#include "XArFileIO_Util.h"
//#include "XArFileIO_FilterProperties.h"
#include <limits>
#include <map>

#include "gdcmPixelFormat.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmRescaler.h"
#include "gdcmTag.h"
#include "gdcmTrace.h"
#include "gdcmImageRegionReader.h"
#include "gdcmImageHelper.h"
#include "gdcmWriter.h"
#include "gdcmTagKeywords.h"
#include "gdcmUIDGenerator.h"

#include <omp.h>
#include "XArray2D.h"
#include "XArray3D.h"

enum eImageFileType
{
	eImageFileType_DICOM_UINT8,
	eImageFileType_DICOM_UINT12,
	eImageFileType_DICOM_UINT16,
	eImageFileType_DICOM_UINT32,
	eImageFileType_DICOM_INT8,
	eImageFileType_DICOM_INT12,
	eImageFileType_DICOM_INT16,
	eImageFileType_DICOM_INT32,
	eImageFileType_DICOM_FLOAT16,
	eImageFileType_DICOM_FLOAT32,
	eImageFileType_DICOM_FLOAT64,
};

class XArFileIO_DICOM
{
public:
	static void ReadFile(xar::XArray2D<float>& rXAr2D, const std::string& strFileName);

	static gdcm::PixelFormat GetDICOMPixelFormatFromFile(const std::string& strFileName);

	template <class T> static void WriteFile(const xar::XArray2D<T>& rXAr2D, std::map<std::string, std::string> mTags, std::string& strFileName, eImageFileType FileType, bool bOptimize)
	{
		if (bOptimize && (FileType > 7) ) throw std::runtime_error("bOptimize == true cannot be used with non-integer FileType in XArFileIO_DICOM::WriteFile()");

		unsigned int width( static_cast<unsigned int>( rXAr2D.GetDim2() ) ); 
		unsigned int height( static_cast<unsigned int>( rXAr2D.GetDim1() ) );
		unsigned int dims[2] = { width, height };

		gdcm::ImageWriter writer;
		gdcm::Image &image( writer.GetImage() );
			
		// Always write files as "Grayscale", ie MONOCHROME2
		gdcm::PhotometricInterpretation pi( gdcm::PhotometricInterpretation::MONOCHROME2 );
		image.SetPhotometricInterpretation( pi );
		image.SetNumberOfDimensions( 2 );
		image.SetDimensions( dims );
			
		// Always write as "Little Endian", non-compressed
		image.SetTransferSyntax( gdcm::TransferSyntax::ExplicitVRLittleEndian );

		double dblSlope(1.0);
		double dblIntercept;
		if (bOptimize) dblIntercept = rXAr2D.Norm(xar::eNormMin);
		else dblIntercept = 0.0; // note that negative values in input data will be truncated away
		gdcm::PixelFormat pfTarget;
		std::vector<char> bufferImage;
		unsigned long len;

		switch( FileType )
		{
			case eImageFileType_DICOM_UINT8:
				pfTarget = gdcm::PixelFormat::UINT8;
				image.SetPixelFormat( pfTarget );
				len = image.GetBufferLength();
				bufferImage.resize(len);
				if (bOptimize)
					rXAr2D.Convert( (unsigned char *)( bufferImage.data() ), T(rXAr2D.Norm( xar::eNormMin )), T(rXAr2D.Norm( xar::eNormMax )), std::numeric_limits<unsigned char>::min(), std::numeric_limits<unsigned char>::max(), &dblSlope, &dblIntercept );
				else
					rXAr2D.Convert((unsigned char*)(bufferImage.data()), false, nullptr, nullptr);
				break;

			case eImageFileType_DICOM_UINT16:
				pfTarget = gdcm::PixelFormat::UINT16;
				image.SetPixelFormat( pfTarget );
				len = image.GetBufferLength();
				bufferImage.resize(len);
				if (bOptimize)
					rXAr2D.Convert( (unsigned short *)( bufferImage.data() ), T(rXAr2D.Norm( xar::eNormMin )), T(rXAr2D.Norm( xar::eNormMax )), std::numeric_limits<unsigned short>::min(), std::numeric_limits<unsigned short>::max(), &dblSlope, &dblIntercept );
				else
					rXAr2D.Convert((unsigned short*)(bufferImage.data()), false, nullptr, nullptr);
				break;

			//case eImageFileType_DICOM_FLOAT32:
			//	pfTarget = gdcm::PixelFormat::FLOAT32;
			//	image.SetPixelFormat( pfTarget );
			//	len = image.GetBufferLength();
			//	bufferImage.resize(len);
			//	rXAr2D.Convert( (float *)( bufferImage.data() ), T(rXAr2D.Norm( xar::eNormMin )), T(rXAr2D.Norm( xar::eNormMax )), float(rXAr2D.Norm(xar::eNormMin)), float(rXAr2D.Norm(xar::eNormMax)), &dblSlope, &dblIntercept );
			//	break;

			default:
				throw std::runtime_error( "Invalid DICOM format" );
		}

		image.SetSlope( 1.0 / dblSlope );
		image.SetIntercept( dblIntercept );

		gdcm::DataElement pixeldata( gdcm::Tag( 0x7fe0,0x0010 ) );
		pixeldata.SetByteValue( (char *)&bufferImage.front(),  image.GetBufferLength() );
		image.SetDataElement( pixeldata );

		writer.SetFileName( strFileName.c_str() );
		
		gdcm::File& file = writer.GetFile();
		gdcm::DataSet& ds = file.GetDataSet();

		gdcm::UIDGenerator uid;
		namespace kwd = gdcm::Keywords;
		kwd::FrameOfReferenceUID frameref;
		frameref.SetValue(uid.Generate());

		kwd::PatientName PatientName = { mTags["PatientName"].c_str() };
		ds.Insert(PatientName.GetAsDataElement());
		kwd::PatientID PatientID = { mTags["PatientID"].c_str() };
		ds.Insert(PatientID.GetAsDataElement());
		kwd::PatientBirthDate PatientBirthDate = { mTags["PatientBirthDate"].c_str() };
		ds.Insert(PatientBirthDate.GetAsDataElement());
		kwd::PatientSex PatientSex = { mTags["PatientSex"].c_str() };
		ds.Insert(PatientSex.GetAsDataElement());
		kwd::StudyID StudyID = { mTags["StudyID"].c_str() };
		ds.Insert(StudyID.GetAsDataElement());
		kwd::AccessionNumber AccessionNumber = { mTags["AccessionNumber"].c_str() };
		ds.Insert(AccessionNumber.GetAsDataElement());
		kwd::StudyInstanceUID StudyInstanceUID = { mTags["StudyInstanceUID"].c_str() };
		ds.Insert(StudyInstanceUID.GetAsDataElement());
		kwd::SeriesInstanceUID SeriesInstanceUID = { mTags["SeriesInstanceUID"].c_str() };
		ds.Insert(SeriesInstanceUID.GetAsDataElement());
		kwd::InstanceNumber InstanceNumber = { atoi(mTags["InstanceNumber"].c_str()) };
		ds.Insert(InstanceNumber.GetAsDataElement());
		kwd::NominalScannedPixelSpacing NominalScannedPixelSpacing = { atof(mTags["NominalScannedPixelSpacing"].c_str()) };
		ds.Insert(NominalScannedPixelSpacing.GetAsDataElement());
		// @@@@ the following lines do not compile with GDCM for some reason, so I had to resort to low-level code
		//kwd::WindowWidth WindowWidth = { atof(mTags["WindowWidth"].c_str()) };
		//ds.Insert(WindowWidth.GetAsDataElement());
		//kwd::WindowCenter WindowCenter = { atof(mTags["WindowCenter"].c_str()) };
		//ds.Insert(WindowCenter.GetAsDataElement());
		kwd::WindowWidth WindowWidth;
		gdcm::DataElement de(gdcm::Tag(0x0028, 0x1051));
		de.SetByteValue(mTags["WindowWidth"].c_str(), 8);
		de.SetVR(gdcm::Attribute<0x0028, 0x1051>::GetVR());
		ds.Insert(de);
		kwd::WindowCenter WindowCenter;
		gdcm::DataElement de1(gdcm::Tag(0x0028, 0x1050));
		de1.SetByteValue(mTags["WindowCenter"].c_str(), 8);
		de1.SetVR(gdcm::Attribute<0x0028, 0x1050>::GetVR());
		ds.Insert(de1);

		if( !writer.Write() )
			throw std::runtime_error( "Error writing DICOM file" );
	}
};

//Function WriteFileStackDICOM
//
//	Writes real XArray3D object into a stack of uncompressed grey-scale DICOM files
//
/*!
	\brief		Writes an XArray3D object into a stack of uncompressed grey-scale DICOM files
	\param		rXAr3D	Reference to an XArray3D<T> object to be saved
	\param		vTags vector of strings that will be written as DICOM tags  
	\param		voutfilenames std::vector of std::strings with full names of the files to be written into
	\param		filetype eDICOM8 or eDICOM16 for 8 or 16 bit DICOM files respectively
	\param		bOptimize Determines if rXAr3D's data is stretched to the full range of the DICOM file (ignored for floating-point data)
	\exception  std::invalid_argument is thrown if XArray3D object contains complex values
	\exception  std::invalid_argument is thrown if 'filetype' parameter is not eDICOM8 or eDICOM16
	\exception  std::runtime_error is thrown if any of the file write operations fail
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function writes a real 3D XArray object into a stack of uncompressed grey-scale DICOM files
*/
template <class T> void WriteFileStackDICOM(xar::XArray3D<T>& rXAr3D, std::map<std::string, std::string> mTags, std::vector<std::string> voutfilenames, xar::_eFileType filetype, bool bOptimize = true)
{
	if (rXAr3D.GetValuetype() == xar::eXAFComplex || rXAr3D.GetValuetype() == xar::eXADComplex)
		throw std::invalid_argument("invalid_argument 'rXAr3D' in WriteFileStackDICOM (complex data)");

	if (!(filetype == xar::eDICOM8 || filetype == xar::eDICOM16)) throw std::runtime_error("unimplemented filetype in WriteFileStackDICOM()");

	index_t nz = rXAr3D.GetDim1();
	if (nz != voutfilenames.size()) throw std::invalid_argument("invalid_argument 'rXAr3D' in XArData::DICOMWriteFilStacke (Dim1 of XArray3D object is different from the length of voutfilenames vector)");
	index_t ny = rXAr3D.GetDim2(); index_t nx = rXAr3D.GetDim3();

	eImageFileType FileType;
	switch (filetype)
	{
	case xar::eDICOM8: FileType = eImageFileType_DICOM_UINT8; break;
	case xar::eDICOM16: FileType = eImageFileType_DICOM_UINT16; break;
	//case xar::eDICOM32: FileType = eImageFileType_DICOM_FLOAT32; break; // this appears to be not implemented in GDCM ver.3.0.18 (Sep 2022)
	default: throw std::invalid_argument("invalid_argument 'filetype' in WriteFileStackDICOM");
	}

	bool bAbort(false);
	#pragma omp parallel for shared (rXAr3D, voutfilenames, nz, ny, nx)
	for (int n = 0; n < nz; n++)
	{
		if (bAbort) continue;
		std::map<std::string, std::string> mTags1(mTags);
		mTags1.insert({ std::string("ImageNumber"), std::to_string(n + 1) }); // 1-based (rather than 0-based) numbering
		mTags1.insert({ std::string("InstanceNumber"), std::to_string(n + 1) }); // 1-based (rather than 0-based) numbering
		try 
		{
			xar::XArray2D<T> ipOut(ny, nx);
			for (index_t j = 0; j < ny; j++)
				for (index_t i = 0; i < nx; i++)
					ipOut[j][i] = rXAr3D[n][j][i];

			XArFileIO_DICOM::WriteFile(ipOut, mTags1, voutfilenames[n], FileType, bOptimize);
		}
		catch (std::exception& E)
		{
			printf("\n\n!!!Exception: %s\n", E.what());
			bAbort = true;
		}
	}
	if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
}


