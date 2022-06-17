//Module XAHWave.cpp
//

//
//	MODULE TITLE:
//
//		Definitions for the xar::WaveheadND classes
//
/*!
	\file		XAHWave.cpp
	\brief		Definitions for the xar::WaveheadND classes
	\par		Description:
		This module contains implementations of the xar::WaveheadND classes
*/

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "stdafx.h"
#include <sstream>
#include <iomanip>

#include "XAHWave.h"

using namespace xar;

//---------------------------------------------------------------------------
//	LOCAL CONSTANT DEFINITIONS
//
//! Percentage inaccuracy allowed in comparing Wavehead1Ds for equivalence
const double Wavehead1D::m_dblDelta = 1.e-4;
//! Percentage inaccuracy allowed in comparing Wavehead2Ds for equivalence
const double Wavehead2D::m_dblDelta = 1.e-4;
//! Percentage inaccuracy allowed in comparing Wavehead3Ds for equivalence
const double Wavehead3D::m_dblDelta = 1.e-4;

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

//---------------------------------------------------------------------------
//Class Wavehead1D
//

//! Copy constructor
Wavehead1D::Wavehead1D(const Wavehead1D& oh)
{
	SetData(oh.m_dblWl, oh.m_dblXlo, oh.m_dblXhi);
}

//! Main constructor
Wavehead1D::Wavehead1D(double dblWl, double dblXlo, double dblXhi)
{
	SetData(dblWl, dblXlo, dblXhi);
}

//! 'Copy' constructor from Wavehead1D object pointed by IXAHead*
Wavehead1D::Wavehead1D(const IXAHead* pHead)
{
	const Wavehead1D* pThat = dynamic_cast<const Wavehead1D*>(pHead);
	if (!pThat) throw std::invalid_argument("invalid_argument 'pHead' in Wavehead1D::Wavehead1D(const IXAHead*) (not a Wavehead1D*)"); 
	SetData(pThat->m_dblWl, pThat->m_dblXlo, pThat->m_dblXhi);
}

//! Sets head data using a generic IXAHead pointer
void Wavehead1D::SetHead(const IXAHead* pHead)
{ 
	const Wavehead1D* pThat = dynamic_cast<const Wavehead1D*>(pHead); 
	if (!pThat) throw std::invalid_argument("invalid_argument 'pHead' in Wavehead1D::SetHead (not a Wavehead1D*)"); 
	*this = *pThat; 
} 

//! Determines if two heads are 'equivalent' (not the same as memberwise equality)	
bool Wavehead1D::Equivalent(const IXAHead* pHead) const
{ 
	const Wavehead1D* pThat = dynamic_cast<const Wavehead1D*>(pHead); 
	if (!pThat) return false;
	return *this == *pThat; 
} 

//! Checks suitability of the data for Wavehead1D
void Wavehead1D::Validate(double dblWl, double dblXlo, double dblXhi) const
{
	if (dblWl <= 0) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead1D::Validate (non-positive dblWl)"); 
	if (dblXlo > dblXhi) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead1D::Validate (dblXlo > dblXhi)"); 
}

//! Returns head data as a formatted string
string Wavehead1D::GetString() const
{
	std::ostringstream ost;
	// NOTE: head type information saved in the string before ':' is used for dynamic creation by CreateXAHead()
	ost << std::scientific << std::setprecision(14) << "Wavehead1D: dblWl = " << m_dblWl <<
		"; dblXlo = " << m_dblXlo << "; dblXhi = " << m_dblXhi << ";";
	return ost.str();
}

//! Sets the head data from an appropriately formatted string
void Wavehead1D::SetString(const string& strData)
{
	index_t n0 = strData.find('=', 0);
	index_t n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead1D::SetString"); 
	double dblWl = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead1D::SetString"); 
	double dblXlo = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead1D::SetString"); 
	double dblXhi = strtod(strData.substr(n0+1, n1).c_str(), 0);
	SetData(dblWl, dblXlo, dblXhi);
}

//! Sets wavelength, lower and higher physical boundaries
void Wavehead1D::SetData(double dblWl, double dblXlo, double dblXhi)
{
	Validate(dblWl, dblXlo, dblXhi);
	m_dblWl = dblWl; 
	m_dblXlo = dblXlo; 
	m_dblXhi = dblXhi; 
}

//! Returns step in physical units
double Wavehead1D::GetStep(index_t NumPoints) const
{
	return (NumPoints > 1) ? (m_dblXhi - m_dblXlo) / NumPoints : 1;
}

//! Resizes the head in accordance with the owner XArray1D
void Wavehead1D::Resize(index_t OldNumPoints, index_t NewNumPoints)
{
	long DifNumPoints = long(NewNumPoints - OldNumPoints);
	long DifNumPoints2 = DifNumPoints / 2;
	double xst = GetStep(OldNumPoints);
	double xhi = m_dblXhi + DifNumPoints2 * xst;
	double xlo = m_dblXlo - (DifNumPoints - DifNumPoints2) * xst;
	if (xhi < xlo) 
	{
		if (xlo - xhi < m_dblDelta * max(fabs(xlo), fabs(xhi)))
		{
			xhi = xlo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'DifNumPoints' in Wavehead1D::Resize"); 
		}
	}
	m_dblXlo = xlo; 
	m_dblXhi = xhi; 
}

//! Trims the head in accordance with the owner XArray1D
void Wavehead1D::Trim(index_t OldNumPoints, index_t lngLeft, index_t lngRight)
{
	double xst = GetStep(OldNumPoints);
	double xlo = m_dblXlo + lngLeft * xst;
	double xhi = m_dblXhi - lngRight * xst;
	if (xhi < xlo) 
	{
		if (xlo - xhi < m_dblDelta * max(fabs(xlo), fabs(xhi)))
		{
			xhi = xlo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'lngLeft or lngRight' in Wavehead1D::Trim"); 
		}
	}
	m_dblXlo = xlo; 
	m_dblXhi = xhi; 
}

//! Pads the head in accordance with the owner XArray1D
void Wavehead1D::Pad(index_t OldNumPoints, index_t lngLeft, index_t lngRight)
{
	double xst = GetStep(OldNumPoints);
	m_dblXlo -= lngLeft * xst;
	m_dblXhi += lngRight * xst;
}

//! Moves the head (separately from the owner XArray1D)
void Wavehead1D::Move(index_t NumPoints, long lngMovePoints)
{
	double xst = GetStep(NumPoints);
	m_dblXlo += lngMovePoints * xst;
	m_dblXhi += lngMovePoints * xst;
}

//! Makes this Wavehead1D a (deep) copy of the rWHead
void Wavehead1D::operator=(const Wavehead1D& oh)
{
	if (&oh == this) return;
	SetData(oh.m_dblWl, oh.m_dblXlo, oh.m_dblXhi);
}

//! Compares two heads for equivalence
bool Wavehead1D::operator==(const Wavehead1D& oh) const
{
	if (fabs(m_dblWl - oh.m_dblWl) > m_dblDelta * m_dblWl) return false;
	if (fabs(m_dblXlo - oh.m_dblXlo) > m_dblDelta * fabs(m_dblXlo)) return false;
	if (fabs(m_dblXhi - oh.m_dblXhi) > m_dblDelta * fabs(m_dblXhi)) return false;
	return true;
}


//---------------------------------------------------------------------------
//Class Wavehead2D
//

//! Copy constructor
Wavehead2D::Wavehead2D(const Wavehead2D& oh)
{
	SetData(oh.m_dblWl, oh.m_dblYlo, oh.m_dblYhi, oh.m_dblXlo, oh.m_dblXhi);
}

//! Main constructor
Wavehead2D::Wavehead2D(double dblWl, double dblYlo, double dblYhi, double dblXlo, double dblXhi)
{
	SetData(dblWl, dblYlo, dblYhi, dblXlo, dblXhi);
}

//! 'Copy' constructor from a Wavehead2D object pointed to by pHead
Wavehead2D::Wavehead2D(const IXAHead* pHead)
{
	const Wavehead2D* pThat = dynamic_cast<const Wavehead2D*>(pHead);
	if (!pThat) throw std::invalid_argument("invalid_argument 'pHead' in Wavehead2D::Wavehead2D(const IXAHead*) (not a Wavehead2D*)"); 
	SetData(pThat->m_dblWl, pThat->m_dblYlo, pThat->m_dblYhi, pThat->m_dblXlo, pThat->m_dblXhi);
}

//! Sets head data using a generic IXAHead pointer
void Wavehead2D::SetHead(const IXAHead* pHead)
{ 
	const Wavehead2D* pThat = dynamic_cast<const Wavehead2D*>(pHead); 
	if (!pThat) throw std::invalid_argument("invalid_argument 'pHead' in Wavehead2D::SetHead (not a Wavehead2D*)"); 
	*this = *pThat; 
} 

//! Determines if two heads are 'equivalent' (not the same as memberwise equality)	
bool Wavehead2D::Equivalent(const IXAHead* pHead) const
{ 
	const Wavehead2D* pThat = dynamic_cast<const Wavehead2D*>(pHead); 
	if (!pThat) return false;
	return *this == *pThat; 
} 

//! Checks suitability of the data for Wavehead2D
void Wavehead2D::Validate(double dblWl, double dblYlo, double dblYhi, double dblXlo, double dblXhi) const
{
	if (dblWl <= 0) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead2D::Validate (non-positive dblWl)"); 
	if (dblYlo > dblYhi) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead2D::Validate (dblYlo > dblYhi)"); 
	if (dblXlo > dblXhi) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead2D::Validate (dblXlo > dblXhi)"); 
}

//! Returns head data as a formatted string
string Wavehead2D::GetString() const
{
	std::ostringstream ost;
	// NOTE: head type information saved in the string before ':' is used for dynamic creation by CreateXAHead()
	ost << std::scientific << std::setprecision(14) << "Wavehead2D: dblWl = " << m_dblWl << "; dblYlo = " << m_dblYlo << "; dblYhi = " << m_dblYhi
		<< "; dblXlo = " <<	m_dblXlo << "; dblXhi = " << m_dblXhi << ";";
	return ost.str();
}

//! Sets the head data from an appropriately formatted string
void Wavehead2D::SetString(const string& strData)
{
	index_t n0 = strData.find('=', 0);
	index_t n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead2D::SetString"); 
	double dblWl = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead2D::SetString"); 
	double dblYlo = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead2D::SetString"); 
	double dblYhi = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead2D::SetString"); 
	double dblXlo = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead2D::SetString"); 
	double dblXhi = strtod(strData.substr(n0+1, n1).c_str(), 0);
	SetData(dblWl, dblYlo, dblYhi, dblXlo, dblXhi);
}

//! Sets wavelength, lower and upper x- and y-boundaries
void Wavehead2D::SetData(double dblWl, double dblYlo, double dblYhi, double dblXlo, double dblXhi)
{
	Validate(dblWl, dblYlo, dblYhi, dblXlo, dblXhi);
	m_dblWl = dblWl; 
	m_dblYlo = dblYlo; 
	m_dblYhi = dblYhi; 
	m_dblXlo = dblXlo; 
	m_dblXhi = dblXhi; 
}

//! Returns y-step in physical units
double Wavehead2D::GetYStep(index_t NumYPoints) const
{
	return (NumYPoints > 1) ? (m_dblYhi - m_dblYlo) / NumYPoints : 1;
}

//! Returns x-step in physical units
double Wavehead2D::GetXStep(index_t NumXPoints) const
{
	return (NumXPoints > 1) ? (m_dblXhi - m_dblXlo) / NumXPoints : 1;
}

//! Resizes the head in accordance with the owner XArray2D
void Wavehead2D::Resize(index_t OldNumYPoints, index_t OldNumXPoints, index_t NewNumYPoints, index_t NewNumXPoints)
{
	long DifNumYPoints = long(NewNumYPoints - OldNumYPoints);
	long DifNumYPoints2 = DifNumYPoints / 2;
	double yst = GetYStep(OldNumYPoints);
	double yhi = m_dblYhi + DifNumYPoints2 * yst;
	double ylo = m_dblYlo - (DifNumYPoints - DifNumYPoints2) * yst;
	if (yhi < ylo) 
	{
		if (ylo - yhi < m_dblDelta * max(fabs(ylo), fabs(yhi)))
		{
			yhi = ylo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'DifNumYPoints' in Wavehead2D::Resize"); 
		}
	}
	m_dblYlo = ylo; 
	m_dblYhi = yhi; 

	long DifNumXPoints = long(NewNumXPoints - OldNumXPoints);
	long DifNumXPoints2 = DifNumXPoints / 2;
	double xst = GetXStep(OldNumXPoints);
	double xhi = m_dblXhi + DifNumXPoints2 * xst;
	double xlo = m_dblXlo - (DifNumXPoints - DifNumXPoints2) * xst;
	if (xhi < xlo) 
	{
		if (xlo - xhi < m_dblDelta * max(fabs(xlo), fabs(xhi)))
		{
			xhi = xlo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'DifNumXPoints' in Wavehead2D::Resize"); 
		}
	}
	m_dblXlo = xlo; 
	m_dblXhi = xhi; 
}

//! Transposes the head (swaps X and Y dimensions) in accordance with the owner XArray2D
void Wavehead2D::Transpose(void)
{
	std::swap(m_dblXlo, m_dblYlo);
	std::swap(m_dblXhi, m_dblYhi);
}

//! Trims the head in accordance with the owner XArray2D
void Wavehead2D::Trim(index_t OldNumYPoints, index_t OldNumXPoints, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight)
{
	double yst = GetYStep(OldNumYPoints);
	double ylo = m_dblYlo + lngYLeft * yst;
	double yhi = m_dblYhi - lngYRight * yst;
	if (yhi < ylo) 
	{
		if (ylo - yhi < m_dblDelta * max(fabs(ylo), fabs(yhi)))
		{
			yhi = ylo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'lngYLeft or lngYRight' in Wavehead2D::Trim"); 
		}
	}
	m_dblYlo = ylo; 
	m_dblYhi = yhi; 

	double xst = GetXStep(OldNumXPoints);
	double xlo = m_dblXlo + lngXLeft * xst;
	double xhi = m_dblXhi - lngXRight * xst;
	if (xhi < xlo) 
	{
		if (xlo - xhi < m_dblDelta * max(fabs(xlo), fabs(xhi)))
		{
			xhi = xlo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'lngXLeft or lngXRight' in Wavehead2D::Trim"); 
		}
	}
	m_dblXlo = xlo; 
	m_dblXhi = xhi; 
}

//! Pads the head in accordance with the owner XArray2D
void Wavehead2D::Pad(index_t OldNumYPoints, index_t OldNumXPoints, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight)
{
	double yst = GetYStep(OldNumYPoints);
	m_dblYlo -= lngYLeft * yst;
	m_dblYhi += lngYRight * yst;

	double xst = GetXStep(OldNumXPoints);
	m_dblXlo -= lngXLeft * xst;
	m_dblXhi += lngXRight * xst;
}

//! Moves the head (separately from the owner XArray2D)
void Wavehead2D::Move(index_t NumYPoints, index_t NumXPoints, long lngMoveYPoints, long lngMoveXPoints)
{
	double yst = GetYStep(NumYPoints);
	m_dblYlo += lngMoveYPoints * yst;
	m_dblYhi += lngMoveYPoints * yst;
	
	double xst = GetXStep(NumXPoints);
	m_dblXlo += lngMoveXPoints * xst;
	m_dblXhi += lngMoveXPoints * xst;
}

//! Makes this Wavehead2D a (deep) copy of the rWHead		
void Wavehead2D::operator=(const Wavehead2D& oh)
{
	if (&oh == this) return;
	SetData(oh.m_dblWl, oh.m_dblYlo, oh.m_dblYhi, oh.m_dblXlo, oh.m_dblXhi);
}

//! Compares two heads for equivalence
bool Wavehead2D::operator==(const Wavehead2D& oh) const
{
	if (fabs(m_dblWl - oh.m_dblWl) > m_dblDelta * m_dblWl) return false;
	if (fabs(m_dblYlo - oh.m_dblYlo) > m_dblDelta * fabs(m_dblYlo)) return false;
	if (fabs(m_dblYhi - oh.m_dblYhi) > m_dblDelta * fabs(m_dblYhi)) return false;
	if (fabs(m_dblXlo - oh.m_dblXlo) > m_dblDelta * fabs(m_dblXlo)) return false;
	if (fabs(m_dblXhi - oh.m_dblXhi) > m_dblDelta * fabs(m_dblXhi)) return false;
	return true;
}


//---------------------------------------------------------------------------
//Class Wavehead3D
//

//! Copy constructor
Wavehead3D::Wavehead3D(const Wavehead3D& oh)
{
	SetData(oh.m_dblWl, oh.m_dblZlo, oh.m_dblZhi, oh.m_dblYlo, oh.m_dblYhi, oh.m_dblXlo, oh.m_dblXhi);
}

//! Main constructor
Wavehead3D::Wavehead3D(double dblWl, double dblZlo, double dblZhi, double dblYlo, double dblYhi, double dblXlo, double dblXhi)
{
	SetData(dblWl, dblZlo, dblZhi, dblYlo, dblYhi, dblXlo, dblXhi);
}

//! 'Copy' constructor from a Wavehead3D object pointed to by pHead
Wavehead3D::Wavehead3D(const IXAHead* pHead)
{
	const Wavehead3D* pThat = dynamic_cast<const Wavehead3D*>(pHead);
	if (!pThat) throw std::invalid_argument("invalid_argument 'pHead' in Wavehead3D::Wavehead3D(const IXAHead*) (not a Wavehead3D*)"); 
	SetData(pThat->m_dblWl, pThat->m_dblZlo, pThat->m_dblZhi, pThat->m_dblYlo, pThat->m_dblYhi, pThat->m_dblXlo, pThat->m_dblXhi);
}

//! Sets head data using a generic IXAHead pointer
void Wavehead3D::SetHead(const IXAHead* pHead)
{ 
	const Wavehead3D* pThat = dynamic_cast<const Wavehead3D*>(pHead); 
	if (!pThat) throw std::invalid_argument("invalid_argument 'pHead' in Wavehead3D::SetHead (not a Wavehead3D*)"); 
	*this = *pThat; 
} 

//! Determines if two heads are 'equivalent' (not the same as memberwise equality)	
bool Wavehead3D::Equivalent(const IXAHead* pHead) const
{ 
	const Wavehead3D* pThat = dynamic_cast<const Wavehead3D*>(pHead); 
	if (!pThat) return false;
	return *this == *pThat; 
} 

//! Checks suitability of the data for Wavehead3D
void Wavehead3D::Validate(double dblWl, double dblZlo, double dblZhi, double dblYlo, double dblYhi, double dblXlo, double dblXhi) const
{
	if (dblWl <= 0) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead3D::Validate (non-positive dblWl)"); 
	if (dblZlo > dblZhi) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead3D::Validate (dblZlo > dblZhi)"); 
	if (dblYlo > dblYhi) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead3D::Validate (dblYlo > dblYhi)"); 
	if (dblXlo > dblXhi) 
		throw std::invalid_argument("invalid_argument '*this' in Wavehead3D::Validate (dblXlo > dblXhi)"); 
}

//! Returns head data as a formatted string
string Wavehead3D::GetString() const
{
	std::ostringstream ost;
	// NOTE: head type information saved in the string before ':' is used for dynamic creation by CreateXAHead()
	ost << std::scientific << std::setprecision(14) << "Wavehead3D: dblWl = " << m_dblWl << "; dblZlo = " << m_dblZlo << "; dblZhi = " << m_dblZhi
		<< "; dblYlo = " << m_dblYlo << "; dblYhi = " << m_dblYhi
		<< "; dblXlo = " <<	m_dblXlo << "; dblXhi = " << m_dblXhi << ";";
	return ost.str();
}

//! Sets the head data from an appropriately formatted string
void Wavehead3D::SetString(const string& strData)
{
	index_t n0 = strData.find('=', 0);
	index_t n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblWl = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblZlo = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblZhi = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblYlo = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblYhi = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblXlo = strtod(strData.substr(n0+1, n1).c_str(), 0);
	n0 = strData.find('=', n1);
	n1 = strData.find(';', n0);
	if (n0 == string::npos || n1 == string::npos) 
		throw std::invalid_argument("invalid_argument 'strData' in Wavehead3D::SetString"); 
	double dblXhi = strtod(strData.substr(n0+1, n1).c_str(), 0);
	SetData(dblWl, dblZlo, dblZhi, dblYlo, dblYhi, dblXlo, dblXhi);
}

//! Sets wavelength, lower and upper x-, y- and z-boundaries
void Wavehead3D::SetData(double dblWl, double dblZlo, double dblZhi, double dblYlo, double dblYhi, double dblXlo, double dblXhi)
{
	Validate(dblWl, dblZlo, dblZhi, dblYlo, dblYhi, dblXlo, dblXhi);
	m_dblWl = dblWl; 
	m_dblZlo = dblZlo; 
	m_dblZhi = dblZhi; 
	m_dblYlo = dblYlo; 
	m_dblYhi = dblYhi; 
	m_dblXlo = dblXlo; 
	m_dblXhi = dblXhi; 
}

//! Returns z-step in physical units
double Wavehead3D::GetZStep(index_t NumZPoints) const
{
	return (NumZPoints > 1) ? (m_dblZhi - m_dblZlo) / NumZPoints : 1;
}

//! Returns y-step in physical units
double Wavehead3D::GetYStep(index_t NumYPoints) const
{
	return (NumYPoints > 1) ? (m_dblYhi - m_dblYlo) / NumYPoints : 1;
}

//! Returns x-step in physical units
double Wavehead3D::GetXStep(index_t NumXPoints) const
{
	return (NumXPoints > 1) ? (m_dblXhi - m_dblXlo) / NumXPoints : 1;
}

//! Resizes the head in accordance with the owner XArray3D
void Wavehead3D::Resize(index_t OldNumZPoints, index_t OldNumYPoints, index_t OldNumXPoints, index_t NewNumZPoints, index_t NewNumYPoints, index_t NewNumXPoints)
{
	long DifNumZPoints = long(NewNumZPoints - OldNumZPoints);
	long DifNumZPoints2 = DifNumZPoints / 2;
	double zst = GetZStep(OldNumZPoints);
	double zhi = m_dblZhi + DifNumZPoints2 * zst;
	double zlo = m_dblZlo - (DifNumZPoints - DifNumZPoints2) * zst;
	if (zhi < zlo) 
	{
		if (zlo - zhi < m_dblDelta * max(fabs(zlo), fabs(zhi)))
		{
			zhi = zlo; // this is onlz to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'DifNumZPoints' in Wavehead3D::Resize"); 
		}
	}
	m_dblZlo = zlo; 
	m_dblZhi = zhi; 

	long DifNumYPoints = long(NewNumYPoints - OldNumYPoints);
	long DifNumYPoints2 = DifNumYPoints / 2;
	double yst = GetYStep(OldNumYPoints);
	double yhi = m_dblYhi + DifNumYPoints2 * yst;
	double ylo = m_dblYlo - (DifNumYPoints - DifNumYPoints2) * yst;
	if (yhi < ylo) 
	{
		if (ylo - yhi < m_dblDelta * max(fabs(ylo), fabs(yhi)))
		{
			yhi = ylo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'DifNumYPoints' in Wavehead3D::Resize"); 
		}
	}
	m_dblYlo = ylo; 
	m_dblYhi = yhi; 

	long DifNumXPoints = long(NewNumXPoints - OldNumXPoints);
	long DifNumXPoints2 = DifNumXPoints / 2;
	double xst = GetXStep(OldNumXPoints);
	double xhi = m_dblXhi + DifNumXPoints2 * xst;
	double xlo = m_dblXlo - (DifNumXPoints - DifNumXPoints2) * xst;
	if (xhi < xlo) 
	{
		if (xlo - xhi < m_dblDelta * max(fabs(xlo), fabs(xhi)))
		{
			xhi = xlo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'DifNumXPoints' in Wavehead3D::Resize"); 
		}
	}
	m_dblXlo = xlo; 
	m_dblXhi = xhi; 
}

//! Trims the head in accordance with the owner XArray3D
void Wavehead3D::Trim(index_t OldNumZPoints, index_t OldNumYPoints, index_t OldNumXPoints, index_t lngZLeft, index_t lngZRight, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight)
{
	double zst = GetZStep(OldNumZPoints);
	double zlo = m_dblZlo + lngZLeft * zst;
	double zhi = m_dblZhi - lngZRight * zst;
	if (zhi < zlo) 
	{
		if (zlo - zhi < m_dblDelta * max(fabs(zlo), fabs(zhi)))
		{
			zhi = zlo; // this is onlz to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'lngZLeft or lngZRight' in Wavehead3D::Trim"); 
		}
	}
	m_dblZlo = zlo; 
	m_dblZhi = zhi; 

	double yst = GetYStep(OldNumYPoints);
	double ylo = m_dblYlo + lngYLeft * yst;
	double yhi = m_dblYhi - lngYRight * yst;
	if (yhi < ylo) 
	{
		if (ylo - yhi < m_dblDelta * max(fabs(ylo), fabs(yhi)))
		{
			yhi = ylo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'lngYLeft or lngYRight' in Wavehead3D::Trim"); 
		}
	}
	m_dblYlo = ylo; 
	m_dblYhi = yhi; 

	double xst = GetXStep(OldNumXPoints);
	double xlo = m_dblXlo + lngXLeft * xst;
	double xhi = m_dblXhi - lngXRight * xst;
	if (xhi < xlo) 
	{
		if (xlo - xhi < m_dblDelta * max(fabs(xlo), fabs(xhi)))
		{
			xhi = xlo; // this is only to account for numerical precision
		}
		else
		{
			throw std::invalid_argument("invalid_argument 'lngXLeft or lngXRight' in Wavehead3D::Trim"); 
		}
	}
	m_dblXlo = xlo; 
	m_dblXhi = xhi; 
}

//! Pads the head in accordance with the owner XArray3D
void Wavehead3D::Pad(index_t OldNumZPoints, index_t OldNumYPoints, index_t OldNumXPoints, index_t lngZLeft, index_t lngZRight, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight)
{
	double zst = GetZStep(OldNumZPoints);
	m_dblZlo -= lngZLeft * zst;
	m_dblZhi += lngZRight * zst;

	double yst = GetYStep(OldNumYPoints);
	m_dblYlo -= lngYLeft * yst;
	m_dblYhi += lngYRight * yst;

	double xst = GetXStep(OldNumXPoints);
	m_dblXlo -= lngXLeft * xst;
	m_dblXhi += lngXRight * xst;
}

//! Moves the head (separately from the owner XArray3D)
void Wavehead3D::Move(index_t NumZPoints, index_t NumYPoints, index_t NumXPoints, long lngMoveZPoints, long lngMoveYPoints, long lngMoveXPoints)
{
	double zst = GetZStep(NumZPoints);
	m_dblZlo += lngMoveZPoints * zst;
	m_dblZhi += lngMoveZPoints * zst;

	double yst = GetYStep(NumYPoints);
	m_dblYlo += lngMoveYPoints * yst;
	m_dblYhi += lngMoveYPoints * yst;
	
	double xst = GetXStep(NumXPoints);
	m_dblXlo += lngMoveXPoints * xst;
	m_dblXhi += lngMoveXPoints * xst;
}

//! Makes this Wavehead3D a (deep) copy of the rWHead		
void Wavehead3D::operator=(const Wavehead3D& oh)
{
	if (&oh == this) return;
	SetData(oh.m_dblWl, oh.m_dblZlo, oh.m_dblZhi, oh.m_dblYlo, oh.m_dblYhi, oh.m_dblXlo, oh.m_dblXhi);
}

//! Compares two heads for equivalence
bool Wavehead3D::operator==(const Wavehead3D& oh) const
{
	if (fabs(m_dblWl - oh.m_dblWl) > m_dblDelta * m_dblWl) return false;
	if (fabs(m_dblZlo - oh.m_dblZlo) > m_dblDelta * fabs(m_dblZlo)) return false;
	if (fabs(m_dblZhi - oh.m_dblZhi) > m_dblDelta * fabs(m_dblZhi)) return false;
	if (fabs(m_dblYlo - oh.m_dblYlo) > m_dblDelta * fabs(m_dblYlo)) return false;
	if (fabs(m_dblYhi - oh.m_dblYhi) > m_dblDelta * fabs(m_dblYhi)) return false;
	if (fabs(m_dblXlo - oh.m_dblXlo) > m_dblDelta * fabs(m_dblXlo)) return false;
	if (fabs(m_dblXhi - oh.m_dblXhi) > m_dblDelta * fabs(m_dblXhi)) return false;
	return true;
}

/////////////////////////////////////////////////////////////////////////////
//
//	End of Module
//
