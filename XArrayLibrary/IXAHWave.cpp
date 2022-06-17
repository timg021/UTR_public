//Module IXAHWave.cpp
//
//
//	MODULE TITLE:
//
//		WaveheadND class factory implementations
//
/*!
	\file		IXAHWave.cpp
	\brief		WaveheadND class factory implementations
	\par		Description:
		This module contains implementations for class factory functions for WaveheadND classes
*/
//---------------------------------------------------------------------------
//	INCLUDE FILES
#include "stdafx.h"

#include "XAHWave.h"
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

//! Creates a Wavehead1D object and returns an interface pointer to it
IXAHWave1D* CreateWavehead1D() { return new xar::Wavehead1D(); }

//! Creates a Wavehead2D object and returns an interface pointer to it
IXAHWave2D* CreateWavehead2D() { return new xar::Wavehead2D(); }

//! Creates a Wavehead3D object and returns an interface pointer to it
IXAHWave3D* CreateWavehead3D() { return new xar::Wavehead3D(); }

//! Creates a Wavehead object of the specified type and returns an interface pointer to it
IXAHead* CreateXAHead(const std::string& strHeadType)
{
	std::string temp = strHeadType.substr(0, 10);
	if (!temp.compare("Wavehead1D")) return CreateWavehead1D();
	if (!temp.compare("Wavehead2D")) return CreateWavehead2D();
	if (!temp.compare("Wavehead3D")) return CreateWavehead3D();
	throw std::invalid_argument("invalid_argument 'strHeadType' in CreateXAHead (unknown head type string)"); 
}
/////////////////////////////////////////////////////////////////////////////
//
//	End of Module
//
