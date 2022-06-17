//Header XA_fft2.h
//
//
//	HEADER FILE TITLE:
//
//		Two-dimensional complex FFT service class
//
//
/*!
	\file		XA_fft2.h
	\brief		Two-dimensional complex FFT service class
	\par		Description:
		This header contains a class that provides 2D complex Fast Fourier Transform (FFT)
		and related services for XArray2D<T> objects
*/
//
// NOTE!!!: in the functions below adapted from XDC, TOMO, XLI, etc., the order of x and y is 
// different from the one in XArray2D. This is handled here by swapping x and y related values at
// the beginning and at the end of the affected functions. Inside such functions the variable
// names, such as e.g. nx or ny, were retained in order not to disturb the well-tested code, 
// but the meaning have now changed. Such function are marked with the mark (old XY order!!!).
// This should be handled with care!

#pragma once

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XArray1D.h"
#include "XA_head2.h"
//#include "crtdbg.h"
#include "OouraFft.h"

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
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//Class XArray2DFFT<T>
//
//	Two-dimensional complex FFT service class
//
/*!
	\brief		Two-dimensional complex FFT service class
	\par		Description:
				This class template provides 2D complex Fast Fourier Transform (FFT) 
				and related services for XArray2D<T> objects	
	\remarks	An object of this class represents a thin 'wrapper' around a complex XArray2D
				object.	The wrapper	is an interface exposing several functions that provide 
				various 2D FFT services				
	\warning	This service class is special in that its template parameter T differs from the
				template parameter of the underlying XArray2D<complex<T> > object. This implies,
				in particular, that this class can only be instantiated for T=float or 
				T=double types.
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XArray2DFFT
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray2DFFT(XArray2D<std::complex<T> >& rXArray2D, bool bFFT = false) : m_rXArray2D(rXArray2D), m_bFFT(bFFT) 
		{ 
			GetValuetype(); 
			double temp = log2(m_rXArray2D.GetDim1());
			if (temp != int(temp)) throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T> constructor (m_dim1 is not an integer power of 2)");
			temp = log2(m_rXArray2D.GetDim2());
			if (temp != int(temp)) throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T> constructor (m_dim2 is not an integer power of 2)");
		}
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XArray2DFFT(const XArray2DFFT<T>& rCopy)
			: m_rXArray2D(rCopy.m_rXArray2D), m_bFFT(rCopy.m_bFFT) {}
	public:
		//! Destructor
		~XArray2DFFT() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XArray2DFFT<T>& rCopy);

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XArray2DFFT<T> objects for types T other than float or double
		//! Returns the xar::_eValueType corresponding to T
		static _eValueType GetValuetype(void);
		//! Returns a reference to the non-modifiable 'wrapped' XArray2D<complex<T> > object
		const XArray2D<std::complex<T> >& GetBaseObject() const { return m_rXArray2D; }
		//! Returns a reference to the 'wrapped' XArray2D<complex<T> > object
		XArray2D<std::complex<T> >& GetBaseObject() { return m_rXArray2D; }
		//! Sets the internal m_bFFT flag indicating if the wrapped m_rXArray2D is in the FFT state
		void SetFftFlag(bool bFFT) { m_bFFT = bFFT; }
		//! Gets the internal m_bFFT flag indicating if the wrapped m_rXArray2D is in the FFT state
		bool GetFftFlag() const { return m_bFFT; }

	// Operations
	public:
		//! Reshuffles the 'wrapped' complex XArray2D<T> object
		void Shuffle();
		//! Performs direct or inverse FFT of the 'wrapped' complex XArray2D object
		void FFT(bool bForward, bool bShuffle = true);
		//! Calculates 2D Kirchhoff integral
		void Kirchhoff(double dblDistance, bool bCheckValidity = true);
		//! Calculates 2D Fresnel integral
		void Fresnel(double dblDistance, bool bCheckValidity = true, double q2max = -1.0, double C3 = 0.0, double C5 = 0.0);
		//! Calculates 2D Fresnel integral with possible astigmatism
		void FresnelA(double dblDistance, double a1, double a2, bool bCheckValidity = true, double q2max = -1.0, double C3 = 0.0, double C5 = 0.0);
		void FresnelAold(double dblDistanceX, double dblDistanceY, bool bCheckValidity, double q2max, double C3, double C5);
		//! Calculates 2D Fresnel integral for long propagation distances
		void FresnelFar(double dblDistance, bool bCheckValidity = true);


	private:
	// Member variables	
		//! Reference to the 'wrapped' XArray2D<complex<T> > object that is being operated upon
		XArray2D<std::complex<T> >& m_rXArray2D;
		//! Determines if the wrapped array is already in the FFT form
		bool m_bFFT;
	};
} // end of namespace xar


//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//! Returns the xar::_eValueType corresponding to T=float
template<> xar::_eValueType xar::XArray2DFFT<float>::GetValuetype() { return eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> xar::_eValueType xar::XArray2DFFT<double>::GetValuetype() { return eXADouble; }


//! Assignment (declared protected to prohibit copying)
template <class T> void xar::XArray2DFFT<T>::operator=(const XArray2DFFT<T>& xaf2)
{ 
	if (this == &xaf2) return; 
	else
	{
		m_rXArray2D = xaf2.m_rXArray2D;
	}
}


//! Reshuffles the 'wrapped' XArray2D<T> object
// NOTE: twice shuffling = identity, hence inverse shuffle = shuffle.
// NOTE: don't mix up with (Un)Wrap which wraps a complex(!) Array2D into a std::real(!) 1D array
// NOTE!!!: changed from old to new XY order 24 April 2022!!!
template <class T> void xar::XArray2DFFT<T>::Shuffle()
{
	index_t ny = m_rXArray2D.GetDim1(), nyd2 = ny / 2, jj;
	index_t nx = m_rXArray2D.GetDim2(), nxd2 = nx / 2;

	for (index_t j = 0; j < nyd2; j++)
	{
		jj = j + nyd2;
		for (index_t i = 0; i < nxd2; i++)
			std::swap(m_rXArray2D[j][i], m_rXArray2D[jj][i + nxd2]);
	}
	for (index_t j = nyd2; j < ny; j++)
	{
		jj = j - nyd2;
		for (index_t i = 0; i < nxd2; i++)
			std::swap(m_rXArray2D[j][i], m_rXArray2D[jj][i + nxd2]);
	}
}


//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::FFT
//
//	Performs direct or inverse FFT of the 'wrapped' complex XArray2D object
//
/*!
	\brief		Performs direct or inverse FFT of the 'wrapped' complex XArray2D object
	\param		bForward	Determines if the direct FFT (true), or the inverse FFT (false) is performed
	\param		bShuffle Determines if the array is shuffled before and after the FFT transform
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function calculates FFT or inverse FFT of the 'wrapped' XArray2D<complex<T> > object
		using the Ooura FFT library
	\par		Example:
\verbatim
XArray2D<dcomplex> C(128, 64, dcomplex(1.0, 0.0)); // create a complex XArray2D object
XArray2DFFT<dcomplex> InterfaceFFT2D(C); // create a 2D FFT interface to the XArray2D object
InterfaceFFT2D.FFT(true); // perform direct FFT of the XArray2D object
\endverbatim
*/	
// NOTE!!!: new X-Y order !!!
template <class T> void xar::XArray2DFFT<T>::FFT(bool bForward, bool bShuffle)
{
	OouraFft<T> fft;

	double wl, ylo, yhi, xlo, xhi;
	IXAHWave2D* ph2 = GetIXAHWave2D(m_rXArray2D);
	if (ph2) //if the head implements IXAHWave2D
	{
		wl = ph2->GetWl();
		ylo = -0.5 / GetYStep(m_rXArray2D);
		yhi = 0.5 / GetYStep(m_rXArray2D);
		xlo = -0.5 / GetXStep(m_rXArray2D);
		xhi = 0.5 / GetXStep(m_rXArray2D);
	}

	if (bShuffle) Shuffle();
	fft.Complex2D((std::complex<T> *) &m_rXArray2D[0][0], 
		m_rXArray2D.GetDim1(), m_rXArray2D.GetDim2(), 
		bForward ? OouraFft<T>::eDirFwd : OouraFft<T>::eDirInv);
	if (bShuffle) Shuffle();
	else m_bFFT = !m_bFFT; // only change the internal m_bFFT flag if the wrapped array has been FFTed and not shuffled

	if (!bForward) 
	{
		T* u = reinterpret_cast<T*>(&(m_rXArray2D.front()));
		index_t nxy = m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2();
		T anxy = T(1.0) / nxy;
		for (index_t k=0; k<nxy*2; k++) u[k] *= anxy;
	}

	if (ph2 && bShuffle)
	{
		IXAHWave2D* ph2a = CreateWavehead2D();
		ph2a->SetData(wl, ylo, yhi, xlo, xhi);
		m_rXArray2D.SetHeadPtr(ph2a);
	}
}

//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::Kirchhoff
//
//	 Calculates 2D Kirchhoff integral
//
/*!
	\brief		Calculates 2D Kirchhoff integral
	\param		dblDistance	Propagation distance (in the same units as used in the Wavehead2D)
	\param		bCheckValidity	Determines the validity of the used implementation for given parameters
	\exception	std::invalid_argument is thrown if any of the two dimensions of the wrapped object
				is not an integer power of 2; or if the object does not have an associated Wavehead2D.
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This program calculates free-space propagation of the scalar complex amplitude
		defined by the 'wrapped' XArray2D object by evaluating the correponding 2D Kirchhoff
		integrals WITH EVANESCENT WAVES using 2D FFT and the spectral (Fourier space)
		representation of the Kirchhoff propagator.
		APPLICABILITY: It can be used when the z-distance between the object and image planes
		is NOT TOO LARGE : (1) sampling condition for the spectral space representation
		(dblDistance << a^2/lambda) & (dblDistance << b^2/lambda),
		or <=> (NFx >> 1) & (NFy >>1)) {NFx(y)=a(b)^2/lambda/dblDistance}.
		NOTE that here the image size in any image plane is always equal to the	initial image
		size in the object plane.
	\par		Example:
\verbatim
XArray2D<dcomplex> C0(128, 64, dcomplex(1.0, 0.0)); // create an incident plane wave
C0.SetHeadPtr(new Wavehead2D(0.0001, -100, 100, -50, 50)); // define the wavelength and physical boundaries
XArray2DFFT<dcomplex> InterfaceFFT2D(C0); // create a 2D FFT interface to the incident wave
InterfaceFFT2D.Kirchhoff(1.e+6); // calculate free space propagation (by 1 m, if units are microns)
\endverbatim
*/	
//  WARNING(old XY order!!!): internally, this function relates dim1 to X, and dim2 to Y, however, 
//	X and Y	are swapped at the entry point of this function
//
template <class T> void xar::XArray2DFFT<T>::Kirchhoff(double dblDistance, bool bCheckValidity)
{
	if (dblDistance==0) return;

	index_t nx = m_rXArray2D.GetDim1();
	index_t ny = m_rXArray2D.GetDim2();

	IXAHWave2D* ph2 = GetIXAHWave2D(m_rXArray2D);
	if (!ph2)
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::Kirchhoff (no Wavehead2D)");
	ph2->Validate();

	double wl = ph2->GetWl(); 
	double xlo = ph2->GetYlo();
	double xhi = ph2->GetYhi();
	double xst = GetYStep(m_rXArray2D); 
	double ylo = ph2->GetXlo();
	double yhi = ph2->GetXhi();
	double yst = GetXStep(m_rXArray2D); 

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap2 = (yhi - ylo) * (yhi - ylo);

	if (T(xap2) == T(0))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::Kirchhoff (xlo==xhi in Wavehead2D)");
	if (T(yap2) == T(0))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::Kirchhoff (ylo==yhi in Wavehead2D)");

//  srad=sqrt(1./(2.*xst)**2+1./(2.*yst)**2)
//  Nyquist frequency (spectral radius) of U0 = srad*wl
//	   if ((srad*wl).gt.1.0) then
//	       write (*,*) 'NOTE: Nyquist frequency > 1 !'
//	       write (*,*) 'Evanescent waves will be present !'
//	   end if

	double absdistance = fabs(dblDistance);
	double fnumx = xap2 / absdistance / wl;
	double fnumy = yap2 / absdistance / wl;
	if (bCheckValidity && (fnumx<10 || fnumy<10))
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::Kirchhoff (Fresnel number too low)");

//********* Fourier transforming initial amplitude u(i,j)

	OouraFft<T> fft;

	T* u = reinterpret_cast<T*>(&(m_rXArray2D.front()));
	if (!m_bFFT) fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(), m_rXArray2D.GetDim2(), OouraFft<T>::eDirFwd);

//********* Multiplying F[u] by the Kirchhoff_propagator

	index_t k, kj;	
	double csi2;
	double dcsi2 = 1.0 / xap2;
	double deta2 = 1.0 / yap2;
	double p2 = 1.0 / wl / wl;
	dcomplex ctemp;
    dcomplex cdfac = dcomplex(0.0, 1.0) * PI * 2.0 * dblDistance;

	for (long i = -long(nxd2); i < 0; i++)
	{
		csi2 = p2 - dcsi2 * i * i;
		kj = nxy2 + ny2 * i + ny2;
		for (long j = -long(nyd2); j < 0; j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(cdfac*std::sqrt(dcomplex(csi2-deta2*j*j)));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
		kj = nxy2 + ny2 * i;
	    for (long j = 0; j < long(nyd2); j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(cdfac*std::sqrt(dcomplex(csi2-deta2*j*j)));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
	}
	for (long i = 0; i < long(nxd2); i++)
	{
		csi2 = p2 - dcsi2 * i * i;
		kj = ny2 * i + ny2;
		for (long j = -long(nyd2); j < 0; j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(cdfac*std::sqrt(dcomplex(csi2-deta2*j*j)));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
		kj = ny2 * i;
		for (long j = 0; j < long(nyd2); j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(cdfac*std::sqrt(dcomplex(csi2-deta2*j*j)));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
	}


//********* inverse Fourier transforming F[u]*Kirchhoff_propagator

	fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(), 
		m_rXArray2D.GetDim2(), OouraFft<T>::eDirInv);
	T fact = T(1.0) / nxy;
	for (k = 0; k < nxy2; k++)	u[k] *= fact;
}


//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::Fresnel
//
//	 Calculates 2D Fresnel integral
//
/*!
	!!! CHANGED 04.06.2019 !!!	
	\brief		Calculates 2D Fresnel integral
	\param		dblDistance	Propagation distance (in the same units as used in the Wavehead2D)
	\param		bCheckValidity	Determines the validity of the used implementation for given parameters
	\param		q2max Defines the optional bandwidth limit (q2max <= 0.0 is interepreted as infinite aperture)
	\param		C3 optional 3rd-order spherical aberration (its dimensionality is the same as for dblDistance)
	\param		C5 optional 5th-order spherical aberration (its dimensionality is the same as for dblDistance)
	\exception	std::invalid_argument is thrown if any of the two dimensions of the wrapped object
				is not an integer power of 2; or if the object does not have an associated Wavehead2D.
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This program calculates paraxial free-space propagation of the scalar complex amplitude
		defined by the 'wrapped' XArray2D object by evaluating the correponding 2D Fresnel
		integrals using 2D FFT and the spectral (Fourier space) representation of the Fresnel propagator.
		APPLICABILITY: It can be used when the z-distance between the object and image planes
		is NOT TOO LARGE : (1) applicability condition of the spectral Fresnel approximation is
		(lambda << dx) & (lambda << dy), or (almost the same) Nyquist frequency of U0 << 1/lambda;
		(2) sampling condition for the spectral space representation (z << a^2/lambda) & (z << b^2/lambda),
		or <=> (NFx >> 1) & (NFy >>1)) {NFx(y)=a(b)^2/lambda/z}.
		NOTE that here the image size in any image plane is always equal to the	initial image
		size in the object plane.
	\par		Example:
\verbatim
XArray2D<dcomplex> C0(128, 64, dcomplex(1.0, 0.0)); // create an incident plane wave
C0.SetHeadPtr(new Wavehead2D(0.0001, -100, 100, -50, 50)); // define the wavelength and physical boundaries
XArray2DFFT<dcomplex> InterfaceFFT2D(C0); // create a 2D FFT interface to the incident wave
InterfaceFFT2D.Fresnel(1.e+6); // calculate free space propagation (by 1 m, if units are microns)
\endverbatim
*/	
//	WARNING(old XY order!!!): internally, this function relates dim1 to X, and dim2 to Y, however, 
//	X and Y	are swapped at the entry point of this function !!!! 4.6.2019 - I AM CHANGING IT NOW
//
template <class T> void xar::XArray2DFFT<T>::Fresnel(double dblDistance, bool bCheckValidity, double q2max, double C3, double C5)
{
	FresnelA(dblDistance, 0, 0, bCheckValidity, q2max, C3, C5);
}


//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::FresnelA
//
//	 Calculates 2D Fresnel integral
//
/*!
	\brief		Calculates 2D Fresnel integral with possible astigamtism
	\param		dblDistance propagation distance (in the same units as used in the Wavehead2D); in the presence of astigmatism, it is equal to (dblDistanceX + dblDistanceY) / 2
	\param		Z1mZ2d2 astigmatism parameter (dblDistanceX - dblDistanceY) / 2 (in the same units as used in the Wavehead2D)
	\param		phiA astigmatism angle (with the X axis, counterclockwise) (in radians)
	\param		bCheckValidity	Determines the validity of the used implementation for given parameters
	\param		q2max Defines the optional bandwidth limit (q2max <= 0.0 is interepreted as infinite aperture)
	\param		C3 optional 3rd-order spherical aberration (its dimensionality is the same as for dblDistanceX and Y)
	\param		C5 optional 5th-order spherical aberration (its dimensionality is the same as for dblDistanceX and Y)
	\exception	std::invalid_argument is thrown if any of the two dimensions of the wrapped object
				is not an integer power of 2; or if the object does not have an associated Wavehead2D.
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This program calculates paraxial free-space propagation of the scalar complex amplitude
		defined by the 'wrapped' XArray2D object by evaluating the correponding 2D Fresnel
		integrals using 2D FFT and the spectral (Fourier space) representation of the Fresnel propagator.
		APPLICABILITY: It can be used when the z-distance between the object and image planes
		is NOT TOO LARGE : (1) applicability condition of the spectral Fresnel approximation is
		(lambda << dx) & (lambda << dy), or (almost the same) Nyquist frequency of U0 << 1/lambda;
		(2) sampling condition for the spectral space representation (z << a^2/lambda) & (z << b^2/lambda),
		or <=> (NFx >> 1) & (NFy >>1)) {NFx(y)=a(b)^2/lambda/z}.
		NOTE that here the image size in any image plane is always equal to the	initial image
		size in the object plane.
	\par		Example:
\verbatim
XArray2D<dcomplex> C0(128, 64, dcomplex(1.0, 0.0)); // create an incident plane wave
C0.SetHeadPtr(new Wavehead2D(0.0001, -100, 100, -50, 50)); // define the wavelength and physical boundaries
XArray2DFFT<dcomplex> InterfaceFFT2D(C0); // create a 2D FFT interface to the incident wave
InterfaceFFT2D.Fresnel2(1.e+6, 1.1e+6); // calculate free space propagation (by 1 m for x and 1.1 m for y, if units are microns)
\endverbatim
*/
template <class T> void xar::XArray2DFFT<T>::FresnelA(double dblDistance, double Z1mZ2d2, double phiA, bool bCheckValidity, double q2max, double C3, double C5)
{
	if (dblDistance == 0) return;

	index_t nx = m_rXArray2D.GetDim2();
	index_t ny = m_rXArray2D.GetDim1();

	IXAHWave2D* ph2 = GetIXAHWave2D(m_rXArray2D);
	if (!ph2)
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelA (no Wavehead2D)");
	ph2->Validate();

	double wl = ph2->GetWl();
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = (yhi - ylo) / ny;
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = (xhi - xlo) / nx;

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;
	double xap = std::abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = std::abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::FresnelA (evanescent waves present)");

	//  Fresnel numbers NFx(y)=a(b)^2/wl/abs(dblDistance)=',fnumx,fnumy
	double absdistanceX = fabs(dblDistance + Z1mZ2d2);
	double absdistanceY = fabs(dblDistance - Z1mZ2d2);
	double fnumx = xap2 / absdistanceX / wl;
	double fnumy = yap2 / absdistanceY / wl;
	if (bCheckValidity && (fnumx < 10 || fnumy < 10))
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::FresnelA (Fresnel number too low)");

	//********* Fourier transforming initial amplitude u(i,j)

	OouraFft<T> fft;

	T* u = reinterpret_cast<T*>(&(m_rXArray2D.front()));
	if (!m_bFFT) fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(), m_rXArray2D.GetDim2(), OouraFft<T>::eDirFwd);


	//********* Multiplying F[u] by the Fresnel_propagator

	index_t k, kj;
	bool bInfAper(q2max <= 0); // infinite aperture
	bool bC35((C3 != 0) || (C5 != 0));
	double eta2, csi2, q2, q4, q6;
	double dcsi = 1.0 / xap;
	double dcsi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;
	dcomplex fac1 = 2.0 * dcomplex(0.0, 1.0) * PI * dblDistance / wl;
	double Z1 = dblDistance + Z1mZ2d2;
	double Z2 = dblDistance - Z1mZ2d2;
	double cosphiA = cos(phiA);
	double sinphiA = sin(phiA);
	double sin2phiA = 2.0 * sinphiA * cosphiA;
	dcomplex fac2x = -dcomplex(0.0, 1.0) * PI * wl * (Z1 * cosphiA * cosphiA + Z2 * sinphiA * sinphiA);
	dcomplex fac2y = -dcomplex(0.0, 1.0) * PI * wl * (Z2 * cosphiA * cosphiA + Z1 * sinphiA * sinphiA);
	dcomplex facxy = -dcomplex(0.0, 1.0) * PI * wl * Z1mZ2d2 * sin2phiA;
	facxy = facxy * deta * dcsi;
	dcomplex fac3 = -dcomplex(0.0, 1.0) * PI * pow(wl, 3) / 2.0 * C3;
	dcomplex fac5 = -dcomplex(0.0, 1.0) * PI * pow(wl, 5) / 3.0 * C5;
	dcomplex ctemp, etafacxy, eta2fac2y, fac2;
	T ttemp;

	for (long i = -long(nyd2); i < 0; i++)
	{
		kj = nxy2 + nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2 = fac1 + eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (bInfAper || q2 < q2max)
			{
				if (bC35)
				{
					q4 = q2 * q2;
					q6 = q4 * q2;
					ctemp = std::exp(fac2 + fac3 * q4 + fac5 * q6);
				}
				else ctemp = std::exp(fac2);
				ttemp = u[k] * T(ctemp.real()) - u[k + 1] * T(ctemp.imag());
				u[k + 1] = u[k] * T(ctemp.imag()) + u[k + 1] * T(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = T(0);
				u[k + 1] = T(0);
			}
		}
		kj = nxy2 + nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2 = fac1 + eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (bInfAper || q2 < q2max)
			{
				if (bC35)
				{
					q4 = q2 * q2;
					q6 = q4 * q2;
					ctemp = std::exp(fac2 + fac3 * q4 + fac5 * q6);
				}
				else ctemp = std::exp(fac2);
				ttemp = u[k] * T(ctemp.real()) - u[k + 1] * T(ctemp.imag());
				u[k + 1] = u[k] * T(ctemp.imag()) + u[k + 1] * T(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = T(0);
				u[k + 1] = T(0);
			}
		}
	}
	for (long i = 0; i < long(nyd2); i++)
	{
		kj = nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2 = fac1 + eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (bInfAper || q2 < q2max)
			{
				if (bC35)
				{
					q4 = q2 * q2;
					q6 = q4 * q2;
					ctemp = std::exp(fac2 + fac3 * q4 + fac5 * q6);
				}
				else ctemp = std::exp(fac2);
				ttemp = u[k] * T(ctemp.real()) - u[k + 1] * T(ctemp.imag());
				u[k + 1] = u[k] * T(ctemp.imag()) + u[k + 1] * T(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = T(0);
				u[k + 1] = T(0);
			}
		}
		kj = nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2 = fac1 + eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (bInfAper || q2 < q2max)
			{
				if (bC35)
				{
					q4 = q2 * q2;
					q6 = q4 * q2;
					ctemp = std::exp(fac2 + fac3 * q4 + fac5 * q6);
				}
				else ctemp = std::exp(fac2);
				ttemp = u[k] * T(ctemp.real()) - u[k + 1] * T(ctemp.imag());
				u[k + 1] = u[k] * T(ctemp.imag()) + u[k + 1] * T(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = T(0);
				u[k + 1] = T(0);
			}
		}
	}

	//********* inverse Fourier transforming F[u]*Kirchhoff_propagator

	fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(),
		m_rXArray2D.GetDim2(), OouraFft<T>::eDirInv);
	T fact = T(1.0 / double(nxy));
	for (k = 0; k < nxy2; k++)	u[k] *= fact;
}


//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::FresnelA
//
//	 Calculates 2D Fresnel integral
//
/*!
	\brief		Calculates 2D Fresnel integral with possible astigamtism
	\param		dblDistanceX Propagation distance corresponding to x (in the same units as used in the Wavehead2D)
	\param		dblDistanceY Propagation distance corresponding to y (in the same units as used in the Wavehead2D)
	\param		bCheckValidity	Determines the validity of the used implementation for given parameters
	\param		q2max Defines the optional bandwidth limit (q2max <= 0.0 is interepreted as infinite aperture)
	\param		C3 optional 3rd-order spherical aberration (its dimensionality is the same as for dblDistanceX and Y)
	\param		C5 optional 5th-order spherical aberration (its dimensionality is the same as for dblDistanceX and Y)
	\exception	std::invalid_argument is thrown if any of the two dimensions of the wrapped object
				is not an integer power of 2; or if the object does not have an associated Wavehead2D.
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This program calculates paraxial free-space propagation of the scalar complex amplitude
		defined by the 'wrapped' XArray2D object by evaluating the correponding 2D Fresnel
		integrals using 2D FFT and the spectral (Fourier space) representation of the Fresnel propagator.
		APPLICABILITY: It can be used when the z-distance between the object and image planes
		is NOT TOO LARGE : (1) applicability condition of the spectral Fresnel approximation is
		(lambda << dx) & (lambda << dy), or (almost the same) Nyquist frequency of U0 << 1/lambda;
		(2) sampling condition for the spectral space representation (z << a^2/lambda) & (z << b^2/lambda),
		or <=> (NFx >> 1) & (NFy >>1)) {NFx(y)=a(b)^2/lambda/z}.
		NOTE that here the image size in any image plane is always equal to the	initial image
		size in the object plane.
	\par		Example:
\verbatim
XArray2D<dcomplex> C0(128, 64, dcomplex(1.0, 0.0)); // create an incident plane wave
C0.SetHeadPtr(new Wavehead2D(0.0001, -100, 100, -50, 50)); // define the wavelength and physical boundaries
XArray2DFFT<dcomplex> InterfaceFFT2D(C0); // create a 2D FFT interface to the incident wave
InterfaceFFT2D.Fresnel2(1.e+6, 1.1e+6); // calculate free space propagation (by 1 m for x and 1.1 m for y, if units are microns)
\endverbatim
*/
template <class T> void xar::XArray2DFFT<T>::FresnelAold(double dblDistanceX, double dblDistanceY, bool bCheckValidity, double q2max, double C3, double C5)
{
	if (dblDistanceX == 0 && dblDistanceY == 0) return;

	index_t nx = m_rXArray2D.GetDim2();
	index_t ny = m_rXArray2D.GetDim1();

	if (log2(ny) != int(log2(ny)))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelA (m_dim1 is not an integer power of 2)");
	if (log2(nx) != int(log2(nx)))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelA (m_dim2 is not an integer power of 2)");

	IXAHWave2D* ph2 = GetIXAHWave2D(m_rXArray2D);
	if (!ph2)
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelA (no Wavehead2D)");
	ph2->Validate();

	double wl = ph2->GetWl();
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = (yhi - ylo) / ny;
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = (xhi - xlo) / nx;

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::FresnelA (evanescent waves present)");

	//  Fresnel numbers NFx(y)=a(b)^2/wl/abs(dblDistance)=',fnumx,fnumy
	double absdistanceX = fabs(dblDistanceX);
	double absdistanceY = fabs(dblDistanceY);
	double fnumx = xap2 / absdistanceX / wl;
	double fnumy = yap2 / absdistanceY / wl;
	if (bCheckValidity && (fnumx < 10 || fnumy < 10))
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::FresnelA (Fresnel number too low)");

	//********* Fourier transforming initial amplitude u(i,j)

	OouraFft<T> fft;

	T* u = reinterpret_cast<T*>(&(m_rXArray2D.front()));
	fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(),
		m_rXArray2D.GetDim2(), OouraFft<T>::eDirFwd);


	//********* Multiplying F[u] by the Fresnel_propagator

	index_t k, kj;
	bool bC35((C3 != 0) || (C5 != 0));
	double eta2, csi2, q2;
	double dcsi2 = 1.0 / xap2;
	double deta2 = 1.0 / yap2;
	//double dtemp = dblDistance / wl - floor(dblDistance / wl);
	//dcomplex fac1 = std::exp(dcomplex(0.0, 1.0) * PI * 2.0 * dtemp);
	dcomplex fac1 = dcomplex(0.0, 1.0) * PI * (dblDistanceX + dblDistanceY) / wl;
	dcomplex fac2x = -dcomplex(0.0, 1.0) * PI * wl * dblDistanceX;
	dcomplex fac2y = -dcomplex(0.0, 1.0) * PI * wl * dblDistanceY;
	dcomplex fac3 = -dcomplex(0.0, 1.0) * PI * pow(wl, 3) / 2.0 * C3;
	dcomplex fac5 = -dcomplex(0.0, 1.0) * PI * pow(wl, 5) / 3.0 * C5;
	dcomplex ctemp, eta2fac2x, fac2a;

	if (q2max > 0)
	{
		for (long i = -long(nyd2); i < 0; i++)
		{
			kj = nxy2 + nx2 * i + nx2;
			eta2 = deta2 * i * i;
			eta2fac2x = eta2 * fac2x;
			for (long j = -long(nxd2); j < 0; j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (q2 < q2max)
				{
					if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
					else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
					u[k] = T(std::real(ctemp));
					u[k + 1] = T(std::imag(ctemp));
				}
				else
				{
					u[k] = T(0);
					u[k + 1] = T(0);
				}
			}
			kj = nxy2 + nx2 * i;
			for (long j = 0; j < long(nxd2); j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (q2 < q2max)
				{
					if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
					else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
					u[k] = T(std::real(ctemp));
					u[k + 1] = T(std::imag(ctemp));
				}
				else
				{
					u[k] = T(0);
					u[k + 1] = T(0);
				}
			}
		}
		for (long i = 0; i < long(nyd2); i++)
		{
			kj = nx2 * i + nx2;
			eta2 = deta2 * i * i;
			eta2fac2x = eta2 * fac2x;
			for (long j = -long(nxd2); j < 0; j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (q2 < q2max)
				{
					if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
					else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
					u[k] = T(std::real(ctemp));
					u[k + 1] = T(std::imag(ctemp));
				}
				else
				{
					u[k] = T(0);
					u[k + 1] = T(0);
				}
			}
			kj = nx2 * i;
			for (long j = 0; j < long(nxd2); j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (q2 < q2max)
				{
					if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
					else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
					u[k] = T(std::real(ctemp));
					u[k + 1] = T(std::imag(ctemp));
				}
				else
				{
					u[k] = T(0);
					u[k + 1] = T(0);
				}
			}
		}
	}
	else
	{
		for (long i = -long(nyd2); i < 0; i++)
		{
			kj = nxy2 + nx2 * i + nx2;
			eta2 = deta2 * i * i;
			eta2fac2x = eta2 * fac2x;
			for (long j = -long(nxd2); j < 0; j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
				else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
				u[k] = T(std::real(ctemp));
				u[k + 1] = T(std::imag(ctemp));
			}
			kj = nxy2 + nx2 * i;
			for (long j = 0; j < long(nxd2); j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
				else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
				u[k] = T(std::real(ctemp));
				u[k + 1] = T(std::imag(ctemp));
			}
		}
		for (long i = 0; i < long(nyd2); i++)
		{
			kj = nx2 * i + nx2;
			eta2 = deta2 * i * i;
			eta2fac2x = eta2 * fac2x;
			for (long j = -long(nxd2); j < 0; j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
				else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
				u[k] = T(std::real(ctemp));
				u[k + 1] = T(std::imag(ctemp));
			}
			kj = nx2 * i;
			for (long j = 0; j < long(nxd2); j++)
			{
				k = kj + 2 * j;
				csi2 = dcsi2 * j * j;
				fac2a = eta2fac2x + csi2 * fac2y;
				q2 = csi2 + eta2;
				if (bC35) ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a + fac3 * q2 * q2 + fac5 * pow(q2, 3));
				else ctemp = dcomplex(u[k], u[k + 1]) * std::exp(fac1 + fac2a);
				u[k] = T(std::real(ctemp));
				u[k + 1] = T(std::imag(ctemp));
			}
		}
	}

	//********* inverse Fourier transforming F[u]*Kirchhoff_propagator

	fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(),
		m_rXArray2D.GetDim2(), OouraFft<T>::eDirInv);
	T fact = T(1.0 / double(nxy));
	for (k = 0; k < nxy2; k++)	u[k] *= fact;
}


//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::FresnelFar
//@@@@@@@@@@@@@ needs to be checked carefully
//
//	 Calculates 2D Fresnel integral for long propagation distances
//
/*!
	\brief		Calculates 2D Fresnel integral for long propagation distances
	\param		dblDistance	Propagation distance (in the same units as used in the Wavehead2D)
	\param		bCheckValidity	Determines the validity of the used implementation for given parameters
	\exception	std::invalid_argument is thrown if any of the two dimensions of the wrapped object
				is not an integer power of 2; or if the object does not have an associated Wavehead2D.
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This program calculates paraxial free-space propagation of the scalar complex amplitude
		defined by the 'wrapped' XArray2D object by evaluating the correponding 2D Fresnel
		integrals using 2D FFT and the conventional (direct space) representation of the Fresnel 
		propagator.
		APPLICABILITY: It can be used when the z-dblDistance between the object and image planes
		is NOT TOO SMALL (but it can be arbitrary large) :
		(1) applicability condition of the conventional Fresnel approximation is
	    (lambda << dx < a << z) & (lambda << dy < b << z),
	    a*b - aperture dimensions, dx&dy - grid step sizes;
		(2) sampling condition for the coordinate space representation
	    (z >> (dx)^2/lambda) & (z >> (dy)^2/lambda),
		or <=> (nx >> sqrt(NFx)) & (ny >> sqrt(NFy)) {NFx(y)=a(b)^2/lambda/z},
		nx=a/dx, ny=b/dy - number of grid points in a(b).
	    NOTE that unlike the Fresnel() program, where the aperture size
		in the image plane is always equal to the aperture size in the
		object plane, here the aperture size in the image plane is,
		generally, different (increasing with z) :
		image_aperture=(a*nx/NFx)*(b*ny/NFy),  steps= a/NFx & b/NFy.
		AUTOREPRODUCTION DISTANCE: NFx=nx & NFy=ny, i.e. z=a*dx/lambda=b*dy/lambda
	\par		Example:
\verbatim
XArray2D<dcomplex> C0(128, 64, dcomplex(1.0, 0.0)); // create an incident plane wave
C0.SetHeadPtr(new Wavehead2D(0.0001, -100, 100, -50, 50)); // define the wavelength and physical boundaries
XArray2DFFT<dcomplex> InterfaceFFT2D(C0); // create a 2D FFT interface to the incident wave
InterfaceFFT2D.FresnelFar(1.e+8); // calculate free space propagation (by 100 m, if units are microns)
\endverbatim
*/	
//	WARNING(old XY order!!!): internally, this function relates dim1 to X, and dim2 to Y, however, 
//	X and Y	are swapped at the entry point of this function
//
template <class T> void xar::XArray2DFFT<T>::FresnelFar(double dblDistance, bool bCheckValidity)
{
	if (m_bFFT) throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelFar (m_rXArray2D is in the FFT state)");

	if (dblDistance==0) return;

	IXAHWave2D* ph2 = GetIXAHWave2D(m_rXArray2D);
	if (!ph2)
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelFar (no Wavehead2D)");
	ph2->Validate();

	double wl = ph2->GetWl(); 
	double xlo = ph2->GetYlo();
	double xhi = ph2->GetYhi();
	double xst = GetYStep(m_rXArray2D);
	double ylo = ph2->GetXlo();
	double yhi = ph2->GetXhi();
	double yst = GetXStep(m_rXArray2D);

	if (T(xhi - xlo) == T(0))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelFar (xlo==xhi in Wavehead2D)");
	if (T(yhi - ylo) == T(0))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelFar (ylo==yhi in Wavehead2D)");

	Shuffle();
	// the head has been destroyed; ph2 is now invalid

	index_t nx = m_rXArray2D.GetDim1();
	index_t ny = m_rXArray2D.GetDim2();
	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;
	double xap = xhi - xlo;
	double yap = yhi - ylo;
	double xap2 = xap * xap;
	double yap2 = yap * yap;
	double xst2 = xst * xst;
	double yst2 = yst * yst;
	double wlz = wl * dblDistance;
	double swlz = sqrt(wlz);
	double sst = xst / swlz;
	double tst = yst / swlz;
	double sst2 = sst * sst;
	double tst2 = tst * tst;
	double smin = xlo / swlz;
	double smax = xhi / swlz;
	double tmin = ylo / swlz;
	double tmax = yhi / swlz;
	double pst = 1.0 / xap * swlz;
	double qst = 1.0 / yap * swlz;
	double pst2 = pst * pst;
	double qst2 = qst * qst;
	double pmin = -0.5 / sst;
	double pmax = -pmin - pst;
	double qmin = -0.5 / tst;
	double qmax = -qmin - qst;
	double xst1 = pst * swlz;
	double yst1 = pst * swlz;
	double xmin1 = pmin * swlz;
	double xmax1 = pmax * swlz;
	double ymin1 = qmin * swlz;
    double ymax1 = qmax * swlz;
	
	xlo = xmin1;
	xhi = xmax1;
	ylo = ymin1;
	yhi = ymax1;
	
	double xdis = xap * xst;
	double ydis = yap * yst;
	//if (xdis==ydis) printf("\nAutoreproduction dblDistance = %g", xdis/wl);

//  (lambda << dx < a << z) & (lambda << dy < b << z)
	double az = fabs(dblDistance);
	if (bCheckValidity && (0.1*xst<wl || 0.1*yst<wl))
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<T>::FresnelFar (grid step is too small)");
	if(bCheckValidity && (0.1*az<xap || 0.1*az<yap))
		throw std::invalid_argument("invalid_argument 'dblDistance' in XArray2DFFT<T>::FresnelFar (propagation distance is too small)");
//	Required: abs(z) > ~10*(xmax-xmin), abs(z) > ~10*(ymax-ymin).

//  (nx >> sqrt(NFx)) & (ny >> sqrt(NFy)) {NFx(y)=a(b)^2/lambda/z}
	double fnumx = xap2 / az / wl;
	double fnumy = yap2 / az / wl;
	double sfnx = sqrt(fnumx);
	double sfny = sqrt(fnumy);
	if (bCheckValidity && (sfnx>0.1*nx || sfny>0.1*ny))
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::FresnelFar (Fresnel number too high)");
//	Results will be unreliable!
//	Required: nx^2 > ~100*NFx, ny^2 > ~100*NFy.
//	Try NFRESNEL instead.

//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
    if (bCheckValidity && srad * wl > 1.0) 
		throw std::runtime_error("runtime_error in XArray2DFFT<T>::FresnelFar (evanescent waves present)");


//********* Multiplying u by the parabolic_term_1

	index_t k, kj;	
	double csi2;
	dcomplex ctemp;
	double dtemp = dblDistance / wl - floor(dblDistance / wl);
	dcomplex fac1 = -dcomplex(0.0, 1.0) * std::exp(dcomplex(0.0, 1.0) * PI * 2.0 * dtemp);
	dcomplex fac2 = dcomplex(0.0, 1.0) * PI;
	T* u = reinterpret_cast<T*>(&(m_rXArray2D.front()));

	for (long i = -long(nxd2); i<0; i++)
	{
		kj = nxy2 + ny2 * i + ny2;
		csi2 = sst2 * i * i;
		for (long j = -long(nyd2); j < 0; j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(fac2 * (tst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
		kj = nxy2 + ny2 * i;
	    for (long j = 0; j < long(nyd2); j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(fac2 * (tst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
	}
	for (long i = 0; i < long(nxd2); i++)
	{
		kj = ny2 * i + ny2;
		csi2 = sst2 * i * i;
		for (long j = -long(nyd2); j < 0; j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(fac2 * (tst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
		kj = ny2 * i;
		for (long j = 0; j < long(nyd2); j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * std::exp(fac2 * (tst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
	}

//********* Calculating inverse Fourier transform

	OouraFft<T> fft;

	fft.Complex2D((std::complex<T> *) u, m_rXArray2D.GetDim1(), 
		m_rXArray2D.GetDim2(), OouraFft<T>::eDirInv);

	double dsdt = sst * tst;
	for (index_t k = 0; k < nxy2; k++) u[k] *= T(dsdt);


//********* Multiplying F[u] by the parabolic_term_2


	for (long i = -long(nxd2); i < 0; i++)
	{
		kj = nxy2 + ny2 * i + ny2;
		csi2 = pst2 * i * i;
		for (long j = -long(nyd2); j < 0; j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * fac1 * std::exp(fac2 * (qst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
		kj = nxy2 + ny2 * i;
	    for (long j = 0; j < long(nyd2); j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * fac1 * std::exp(fac2 * (qst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
	}
	for (long i = 0; i < long(nxd2); i++)
	{
		kj = ny2 * i + ny2;
		csi2 = pst2 * i * i;
		for (long j = -long(nyd2); j < 0; j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * fac1 * std::exp(fac2 * (qst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
		kj = ny2 * i;
		for (long j = 0; j < long(nyd2); j++)
		{
			k = kj + 2 * j;
			ctemp = dcomplex(u[k], u[k+1]) * fac1 * std::exp(fac2 * (qst2 * j * j + csi2));
			u[k] = T(std::real(ctemp));
			u[k+1] = T(std::imag(ctemp));
		}
	}

	Shuffle();

	//create and attach a new head ( x and y reversed)
	ph2 = CreateWavehead2D();
	ph2->SetData(wl, xlo, xhi, ylo, yhi);
	m_rXArray2D.SetHeadPtr(ph2);
}

//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//

// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray2DFFT<float>;
	template class xar::XArray2DFFT<double>;
#endif

//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//

/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
