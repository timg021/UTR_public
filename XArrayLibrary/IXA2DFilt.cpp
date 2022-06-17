//Module IXA2DFilt.cpp
//
//
//	MODULE TITLE:
//
//		XArray2DFiltImpl class factory implementation
//
/*!
	\file		IXA2DFilt.cpp
	\brief		XArray2DFiltImpl class factory implementation
	\par		Description:
		This module contains class factory implementation for the class providing 2D 
		noise generation and filtering and related services for XArray2D objects. 
		This facility is required by client programs accessing XArray functionality
		via the interfaces
*/
//---------------------------------------------------------------------------
//	INCLUDE FILES
#include "XArray/XAr2DImp.h"
#include "XArray/XAr2DFiltImp.h" //for XArray2DFiltImpl

//---------------------------------------------------------------------------
//	LOCAL CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	LOCAL STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL IN-LINE FUNCTION DEFINITIONS
//
//---------------------------------------------------------------------------
//	GLOBAL FUNCTION PROTOTYPES
//
//---------------------------------------------------------------------------
//	LOCAL FUNCTION PROTOTYPES
//
//---------------------------------------------------------------------------
//	GLOBAL DATA DEFINITIONS
//
//---------------------------------------------------------------------------
//	STATIC DATA DEFINITIONS
//
/////////////////////////////////////////////////////////////////////////////
//

using namespace xar;

//! XArray2DFiltImpl<T> class factory; creates an XArray2DFiltImpl<T> object and returns an IXArray2DFiltPtr interface pointer to it
XA_API IXArray2DFiltPtr CreateXA2DFilt(IXArray2D* ixa2)
{
	const std::type_info& ti = typeid(*ixa2); 
	if (ti==typeid(xar::XArray2DImpl<char>)) 
		return IXArray2DFiltPtr(new XArray2DFiltImpl<char>(*reinterpret_cast<xar::XArray2D<char>*>(ixa2->GetXArray())));

	if (ti==typeid(xar::XArray2DImpl<short>))
		return IXArray2DFiltPtr(new XArray2DFiltImpl<short>(*reinterpret_cast<xar::XArray2D<short>*>(ixa2->GetXArray())));

	if (ti==typeid(xar::XArray2DImpl<long>))
		return IXArray2DFiltPtr(new XArray2DFiltImpl<long>(*reinterpret_cast<xar::XArray2D<long>*>(ixa2->GetXArray())));

	if (ti==typeid(xar::XArray2DImpl<float>))
		return IXArray2DFiltPtr(new XArray2DFiltImpl<float>(*reinterpret_cast<xar::XArray2D<float>*>(ixa2->GetXArray())));

	if (ti==typeid(xar::XArray2DImpl<double>))
		return IXArray2DFiltPtr(new XArray2DFiltImpl<double>(*reinterpret_cast<xar::XArray2D<double>*>(ixa2->GetXArray())));

	if (ti==typeid(xar::XArray2DImpl<fcomplex>))
		return IXArray2DFiltPtr(new XArray2DFiltImpl<fcomplex>(*reinterpret_cast<xar::XArray2D<fcomplex>*>(ixa2->GetXArray())));

	if (ti==typeid(xar::XArray2DImpl<dcomplex>))
		return IXArray2DFiltPtr(new XArray2DFiltImpl<dcomplex>(*reinterpret_cast<xar::XArray2D<dcomplex>*>(ixa2->GetXArray())));

	throw  std::invalid_argument("invalid_argument 'IXArray2D* ixa2' in CreateXA2DFilt()");
}
