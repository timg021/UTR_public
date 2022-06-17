//Header XA_head1.h
//
//
//	HEADER FILE TITLE:
//
//		XArray1D-specific simple helper functions depending on IXAHWave1D
//
/*!
	\file		XA_head1.h
	\brief		XArray1D-specific simple helper functions depending on IXAHWave1D
	\par		Description:
		Contains simple helper functions requiring knowledge of both XArray1D<T>
		object and IXAHWave1D interface (to a Wavehead1D object)
*/
#ifndef XA_HEAD1_H
#define XA_HEAD1_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "IXAHWave.h"
#include "XArray1D.h"

namespace xar
{
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//
	//
	// The following functions are not strictly necessary, but they can make 
	// template-level programming with XArrays and Waveheads easier
	//

	//! Returns a pointer to the const IXAHWave1D interface or zero if it is unavailable
	template <class T> inline const IXAHWave1D* GetIXAHWave1D(const XArray1D<T>& rXAr1D)
	{
		return dynamic_cast<const IXAHWave1D*>(rXAr1D.GetHeadPtr());
	}

	//! Returns a pointer to the IXAHWave1D interface or zero if it is unavailable
	template <class T> inline IXAHWave1D* GetIXAHWave1D(XArray1D<T>& rXAr1D)
	{
		return dynamic_cast<IXAHWave1D*>(rXAr1D.GetHeadPtr());
	}

	//! Returns X-step or 1 if an IXAHWave1D interface is unavailable
	template <class T> double GetStep(const XArray1D<T>& rXAr1D)
	{
		const IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1D);
		if (!ph1) return 1;
		else return ph1->GetStep(rXAr1D.size());
	}

	//! Returns lower X-boundary or 0 if an IXAHWave1D interface is unavailable
	template <class T> inline double GetXlo(const XArray1D<T>& rXAr1D)
	{
		const IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1D);
		if (!ph1) return 0;
		else return ph1->GetXlo();
	}

	//! Resizes the 1D array and the head as well if it is present
	template <class T> void ResizeH(XArray1D<T>& rXAr1D, index_t iNewDim, T tVal = T())
	{
		IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1D);
		if (ph1) ph1->Resize(rXAr1D.size(), iNewDim);
		rXAr1D.Resize(iNewDim, tVal);
	}

	//! Extracts a sub-XArray1D into another XArray1D with an appropriate head
	template <class T> void GetSubarrayH(const XArray1D<T>& rXAr1DSrc, index_t iBegin, index_t iEnd, XArray1D<T>& rDestSubXArray)
	{
		rXAr1DSrc.GetSubarray(iBegin, iEnd, rDestSubXArray);
		const IXAHWave1D* ph1 = GetIXAHWave1D(rXAr1DSrc);
		if (ph1) 
		{
			IXAHWave1D* ph1new = CreateWavehead1D();
			double xlo = ph1->GetXlo() + GetStep(rXAr1DSrc) * iBegin;
			double xhi = xlo + GetStep(rXAr1DSrc) * (iEnd - iBegin - 1);
			ph1new->SetData(ph1->GetWl(), xlo, xhi);
			rDestSubXArray.SetHeadPtr(ph1new);
		}
	}

	
} // namespace xar closed
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XA_HEAD1_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
