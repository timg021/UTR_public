#include "XA_DICOM.h"
#include "gdcmRescaler.h"

void XArFileIO_DICOM::ReadFile(	xar::XArray2D<float>& rXAr2D, const std::string& strFileName)
{
	std::vector<char> bufferImage;
	unsigned int nWidth, nHeight;
	std::string strPixelType;

	// Instantiate the reader:
	gdcm::ImageReader reader;
	reader.SetFileName( strFileName.c_str() );

	if( !reader.Read() )
		throw std::runtime_error("error reading DICOM file");

	// If we reach here, we know for sure only 1 thing:
	// It is a valid DICOM file (potentially an old ACR-NEMA 1.0/2.0 file)
	// (Maybe, it's NOT a Dicom image -could be a DICOMDIR, a RTSTRUCT, etc-)

	const gdcm::Image &image( reader.GetImage() );


	// Determine the file type from the DICOM pixel format
	switch(image.GetPixelFormat())
	{
		case gdcm::PixelFormat::UINT12:
		case gdcm::PixelFormat::INT12:
		case gdcm::PixelFormat::FLOAT16:
			throw std::runtime_error("Currently unsupported DICOM pixel format");
	}
		
	// Determine the file color space
	if(image.GetPhotometricInterpretation() != gdcm::PhotometricInterpretation::MONOCHROME1 && image.GetPhotometricInterpretation() != gdcm::PhotometricInterpretation::MONOCHROME2)
		throw std::runtime_error("Currently unsupported DICOM colorspace - only grayscale is supported");

	// Determine compression/Endiness from Transfer Syntax
	nHeight = image.GetRows();
	nWidth  = image.GetColumns();

	unsigned long len( image.GetBufferLength() );
	bufferImage.resize( len );

	// Get a pointer to the read image data
	if( !image.GetBuffer( &bufferImage.front() ) )
		throw std::runtime_error("error accessing DICOM raw data buffer");

	// prepare the array
	rXAr2D.Resize( nHeight, nWidth );	

	gdcm::Rescaler rescaler;
	gdcm::PixelFormat pfFloat( gdcm::PixelFormat::FLOAT32 );

	// Set slope, intercept and pixel format (always float)
	rescaler.SetIntercept( image.GetIntercept() );
	rescaler.SetSlope( image.GetSlope() );
	rescaler.SetPixelFormat( image.GetPixelFormat() );
	rescaler.SetTargetPixelType( pfFloat );
	rescaler.SetUseTargetPixelType( true );
		
	if( !rescaler.Rescale( (char *)( &rXAr2D.front() ), &bufferImage.front(), image.GetBufferLength() ) )
		throw std::runtime_error("error rescaling DICOM file");
}

gdcm::PixelFormat XArFileIO_DICOM::GetDICOMPixelFormatFromFile(const std::string& strFileName)
{
	// Instanciate the reader:
	gdcm::ImageReader reader;
	reader.SetFileName( strFileName.c_str() );
		
	if( !reader.Read() )
		throw std::runtime_error("error reading DICOM file");

	// If we reach here, we know for sure only 1 thing:
	// It is a valid DICOM file (potentially an old ACR-NEMA 1.0/2.0 file)
	// (Maybe, it's NOT a Dicom image -could be a DICOMDIR, a RTSTRUCT, etc-)

	const gdcm::Image &image( reader.GetImage() );

	// Determine the file type from the DICOM pixel format
	return image.GetPixelFormat();
}
