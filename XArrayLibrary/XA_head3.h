//Header XA_head3.h
//
//
//	HEADER FILE TITLE:
//
//		XArray3D specific simple helper functions depending on IXAHWave3D
//
/*!
	\file		XA_head3.h
	\brief		XArray3D specific simple helper functions depending on IXAHWave3D
	\par		Description:
		Contains simple helper functions requiring knowledge of both XArray3D<T>
		object and IXAHWave3D interface (to a Wavehead3D object)
*/

#ifndef XA_HEAD3_H
#define XA_HEAD3_H

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "IXAHWave.h"
#include "XArray3D.h"

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

	//! Returns a pointer to the const IXAHWave3D interface or zero if it is unavailable
	template <class T> inline const IXAHWave3D* GetIXAHWave3D(const XArray3D<T>& rXAr3D)
	{
		return dynamic_cast<const IXAHWave3D*>(rXAr3D.GetHeadPtr());
	}

	//! Returns a pointer to the IXAHWave3D interface or zero if it is unavailable		
	template <class T> inline IXAHWave3D* GetIXAHWave3D(XArray3D<T>& rXAr3D)
	{
		return dynamic_cast<IXAHWave3D*>(rXAr3D.GetHeadPtr());
	}
	
	//! Returns Z-step or 1 if an IXAHWave3D interface is unavailable
	template <class T> inline double GetZStep(const XArray3D<T>& rXAr3D)
	{
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (!ph3) return 1;
		else return ph3->GetZStep(rXAr3D.GetDim1());
	}

	//! Returns Y-step or 1 if an IXAHWave3D interface is unavailable
	template <class T> inline double GetYStep(const XArray3D<T>& rXAr3D)
	{
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (!ph3) return 1;
		else return ph3->GetYStep(rXAr3D.GetDim2());
	}

	//! Returns X-step or 1 if an IXAHWave3D interface is unavailable
	template <class T> inline double GetXStep(const XArray3D<T>& rXAr3D)
	{
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (!ph3) return 1;
		else return ph3->GetXStep(rXAr3D.GetDim3());
	}

	//! Returns lower Z-boundary or 0 if an IXAHWave3D interface is unavailable
	template <class T> inline double GetZlo(const XArray3D<T>& rXAr3D)
	{
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (!ph3) return 0;
		else return ph3->GetZlo();
	}

	//! Returns lower Y-boundary or 0 if an IXAHWave3D interface is unavailable
	template <class T> inline double GetYlo(const XArray3D<T>& rXAr3D)
	{
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (!ph3) return 0;
		else return ph3->GetYlo();
	}

	//! Returns lower X-boundary or 0 if an IXAHWave3D interface is unavailable
	template <class T> inline double GetXlo(const XArray3D<T>& rXAr3D)
	{
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (!ph3) return 0;
		else return ph3->GetXlo();
	}

	//! Resizes the 3D array and the head as well if it is present		
	template <class T> inline void ResizeH(XArray3D<T>& rXAr3D, index_t iNewDim1, index_t iNewDim2, index_t iNewDim3, T tVal = T())
	{
		IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3D);
		if (ph3) ph3->Resize(rXAr3D.GetDim1(), rXAr3D.GetDim2(), rXAr3D.GetDim3(), iNewDim1, iNewDim2, iNewDim3);
		rXAr3D.Resize(iNewDim1, iNewDim2, iNewDim3, tVal);
	}

	//! Extracts a sub-XArray3D into another XArray3D with an appropriate head
	template <class T> void GetSubarrayH(const XArray3D<T>& rXAr3DSrc, index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, index_t iBeginDim3, index_t iEndDim3, XArray3D<T>& rDestSubXArray)
	{
		rXAr3DSrc.GetSubarray(iBeginDim1, iEndDim1, iBeginDim2, iEndDim2, iBeginDim3, iEndDim3, rDestSubXArray);
		const IXAHWave3D* ph3 = GetIXAHWave3D(rXAr3DSrc);
		if (ph3) 
		{
			IXAHWave3D* ph3new = CreateWavehead3D();
			double zlo = ph3->GetZlo() + GetZStep(rXAr3DSrc) * iBeginDim1;
			double zhi = zlo + GetZStep(rXAr3DSrc) * (iEndDim1 - iBeginDim1 - 1);
			double ylo = ph3->GetYlo() + GetYStep(rXAr3DSrc) * iBeginDim2;
			double yhi = ylo + GetYStep(rXAr3DSrc) * (iEndDim2 - iBeginDim2 - 1);
			double xlo = ph3->GetXlo() + GetXStep(rXAr3DSrc) * iBeginDim3;
			double xhi = xlo + GetXStep(rXAr3DSrc) * (iEndDim3 - iBeginDim3 - 1);
			ph3new->SetData(ph3->GetWl(), zlo, zhi, ylo, yhi, xlo, xhi);
			rDestSubXArray.SetHeadPtr(ph3new);
		}
	}

} // namespace xar closed

#endif // XA_HEAD3_H