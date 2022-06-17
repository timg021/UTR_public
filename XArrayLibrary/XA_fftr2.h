//Header XA_fftr2.h
//
//
//	HEADER FILE TITLE:
//
//		Two-dimensional real FFT service class
//
/*!
	\file		XA_fftr2.h
	\brief		Two-dimensional real FFT service class
	\par		Description:
		This header contains a class that provides 2D real Fast Fourier Transform (FFT)
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

#include "OouraFft.h"
#include "XA_fft2.h"
// #include "XA_Optimisation.h"
#include "XA_ini.h"

namespace xar
{
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
	template <class T> class XArray2DMove;
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
//Class XArray2DFFTRe<T>
//
//	Two-dimensional real FFT service class
//
/*!
	\brief		Two-dimensional real FFT service class
	\par		Description:
				This class template provides 2D real Fast Fourier Transform (FFT) 
				and related services for XArray2D<T> objects	
	\remarks	An object of this class represents a thin 'wrapper' around a real XArray2D
				object.	The wrapper	is an interface exposing several functions that provide 
				various 2D FFT services				
	\warning	This service class can only be created for T=float or T=double.
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XArray2DFFTRe
	{
	// Typedefs
	public:
		typedef T type;

	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray2DFFTRe(XArray2D<T>& rXArray2D) : m_rXArray2D(rXArray2D) { GetValuetype(); }
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XArray2DFFTRe(const XArray2DFFTRe<T>& rCopy) : m_rXArray2D(rCopy.m_rXArray2D), m_vTemp(rCopy.m_vTemp) {}
	public:
		//! Destructor
		~XArray2DFFTRe() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XArray2DFFTRe<T>& rCopy);

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XArray2DFFTRe<T> objects for types T other than float or double
		//! Returns the xar::_eValueType corresponding to T
		static _eValueType GetValuetype(void);
		//#ifdef XA_INTEL_MKL
		//	static DFTI_CONFIG_VALUE GetValuetypeMKL(void);
		//#endif
		//! Returns a reference to the non-modifiable 'wrapped' XArray2D<T> object
		const XArray2D<T>& GetBaseObject() const { return m_rXArray2D; }
		//! Returns a reference to the 'wrapped' XArray2D<T> object
		XArray2D<T>& GetBaseObject() { return m_rXArray2D; }

	// Operations
	public:
		//! Reshuffles the 'wrapped' XArray2D<T> object
		void Shuffle();
		//! Performs direct or inverse FFT of the 'wrapped' real XArray2D object
		void FFTRe(XArray2D<std::complex<T> >& C, bool bForward, bool bCheck = true, bool bTruncate = true);
		//! Convolves the 'wrapped' real XArray2D object with another XArray2D object 
		void Convol(XArray2D<T>& rxa2DKer, bool bTrimBack = false);
		//! Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method
		void Deconvol(XArray2D<T>& rxa2DKer, double dblAlpha, bool bRescale, bool bTrimBack = false/*, void* ppMklDftiDescriptor = 0*/);
		//! Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method optimizing alpha
		//void DeconvolOpt(XArray2D<T>& rxa2DKer, double dblSigma, bool bRescale, double* pdblAlphaOpt);
		//! Iteratively deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method
		void DeconvolIter(XArray2D<T>& rxa2DKer, double dblAlpha, double dblSigma, bool bRescale, index_t& riNumIter);
		//! Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using Richardson-Lucy algorithm
		void RichLucy(XArray2D<T>& rxa2DKer, index_t& riNumIter, double dblSigma, bool bTrimBack = false);
		//! Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using Taylor series
		void DeconvolTaylor(XArray2D<T>& rxa2DKer, index_t iMaxDerOrder, double dblFiltSigmaY, double dblFiltSigmaX, double dblRst_dev_av, index_t& riNumIter);
		//! Phase-correlates the 'wrapped' real XArray2D object with another XArray2D object 
		T Correlate(XArray2D<T>& rxa2DOther, long& rlngShift1, long& rlngShift2, bool bTruncateThis = true, bool bTruncateRefArr = true);
		//! Correlates the 'wrapped' real XArray2D object with another XArray2D object 
		T CorrelateOld(XArray2D<T>& rxa2DOther, long& rlngShift1, long& rlngShift2, bool bTruncateThis, bool bTruncateRefArr);
		//! Filters the 'wrapped' real XArray2D object by convoluting it with a rectangular kernel
		void FilterRect(index_t iYWidth, index_t iXWidth);
		//!	Filters the 'wrapped' real XArray2D object by convoluting it with an Epanechnikov kernel
		void FilterEpanechnikov(double dblWidthY, double dblWidthX);
		//!	Filters the 'wrapped' real XArray2D object by convoluting it with an axially symmetrical Epanechnikov kernel
		void FilterSymEpanechnikov(double dblWidth);
		//! Filters the 'wrapped' real XArray2D object by convoluting it with a Gaussian kernel
		void FilterGauss(double dblSigmaY, double dblSigmaX);
		//! Filters the 'wrapped' real XArray2D object by convoluting it with a circularly symmetric Lorentz kernel
		void FilterLorentz(double dblSigma);
		//! Calculates 2D Fourier filter
		void FilterFourier(double SpectFraction);
		//! Filters the 'wrapped' real XArray2D object by applying rectangular window to its Fourier spectrum
		void FourFilterRect(index_t iYFreqCutoff, index_t iXFreqCutoff);
		//! Filters the 'wrapped' real XArray2D object by applying Gaussian window to its Fourier spectrum
		void FourFilterGauss(double dblFreqSigmaY, double dblFreqSigmaX);
		//! Filters the 'wrapped' real XArray2D object by applying Lorentzian window to its Fourier spectrum
		void FourFilterLorentz(double dblGamma, index_t iKy0, index_t iKx0);

		// Auxilliary function called by the optimization routine 'goldenC'
		T CalcOnce(T* ArrayPars, index_t What2Calc);
		
	// Implementation
	protected:

	private:
	// Member variables	
		//! Reference to the 'wrapped' XArray2D<T> object that is being operated upon
		XArray2D<T>& m_rXArray2D;
		//! Auxilliary vector used for passing parameters between some member functions
		std::vector<void*> m_vTemp;
	};

	//---------------------------------------------------------------------------
	//	TEMPLATE MEMBER DEFINITIONS
	//
	//! Returns the xar::_eValueType corresponding to T=float
	template<> inline _eValueType XArray2DFFTRe<float>::GetValuetype() { return eXAFloat; }
	//! Returns the xar::_eValueType corresponding to T=double
	template<> inline _eValueType XArray2DFFTRe<double>::GetValuetype() { return eXADouble; }
	//#ifdef XA_INTEL_MKL
		//! Returns the DFTI_CONFIG_VALUE corresponding to T=float
	//	inline DFTI_CONFIG_VALUE xar::XArray2DFFTRe<float>::GetValuetypeMKL() { return DFTI_SINGLE; }
		//! Returns the DFTI_CONFIG_VALUE corresponding to T=double
	//	inline DFTI_CONFIG_VALUE xar::XArray2DFFTRe<double>::GetValuetypeMKL() { return DFTI_DOUBLE; } 
	//#endif


	//! Assignment (declared protected to prohibit copying)
	template <class T> void XArray2DFFTRe<T>::operator=(const XArray2DFFTRe<T>& xaf2)
	{ 
		if (this == &xaf2) 
			return; 
		else
		{
			m_rXArray2D = xaf2.m_rXArray2D;
			m_vTemp = xaf2.m_vTemp;
		}
	}


	//! Reshuffles the 'wrapped' XArray2D<T> object
	// NOTE: twice shuffling = identity, hence inverse shuffle = shuffle.
	// NOTE: don't mix up with (Un)Wrap which wraps a complex(!) Array2D into a std::real(!) 1D array
	// NOTE!!!: old XY order!!!
	template <class T> void XArray2DFFTRe<T>::Shuffle()
	{
		index_t nx = m_rXArray2D.GetDim1(), nxd2 = nx / 2, i;
		index_t ny = m_rXArray2D.GetDim2(), nyd2 = ny / 2, j;

		for (i=0; i<nxd2; i++)
		{
			for (j=0; j<nyd2; j++) 
				std::swap(m_rXArray2D[i][j], m_rXArray2D[i+nxd2][j+nyd2]);

			for (j=nyd2; j<ny; j++)	
				std::swap(m_rXArray2D[i][j], m_rXArray2D[i+nxd2][j-nyd2]);
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FFTRe
	//
	//	Performs direct (inverse) FFT of (into) the 'wrapped' real XArray2D object
	//
	/*!
		\brief		Performs direct (inverse) FFT of (into) the 'wrapped' real XArray2D object
		\param		C	Target (source) complex XArray2D object with appropriate symmetries
		\param		bForward	Determines if the direct FFT (true), or the inverse FFT (false) is performed
		\param		bCheck	Determines if the array dimensions are checked to be integer powers of 2 (true)
		\param		bResizeToZero	Determines if the source object is resized to zero in this function
		\exception  std::invalid_argument is thrown if bCheck is true and the array dimensions are not integer powers of 2
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates FFT or inverse FFT of the 'wrapped' XArray2D<T> object
			using the Ooura FFT library. If bForward is true, then the forward FFT 
			is performed from the real 'wrapped' object into a complex XArray2D, C, with the 
			conjugate-symmetry property (if bResizeToZero is true, then  the real 'wrapped' object 
			is resized to zero). If bForward is false, then the inverse FFT is performed from the
			complex XArray2D, C, with the conjugate-symmetry property (not checked!) into the real
			'wrapped' object (if bResizeToZero is true, then C is resized to zero).
		\par		Example:
	\verbatim
	XArray2D<double> xa2D; // create a real XArray2D
	xa2D.Resize(4, 8); // resize it
	xa2D.Fill(0.0); // fill it with zeros
	xa2D[2][4] = 32; // set the centre value to 32
	XArray2DFFTRe<double> InterfaceFFTRe(xa2D); // create a XArray2DFFTRe interface to "xa2D"
	XArray2D<dcomplex> xc2D; // create a complex XArray2D (target of the real FFT)
	InterfaceFFTRe.FFTRe(xc2D, true); // calculate forward FFT of "xa2D" and save the result in "xc2D"
	\endverbatim
	*/	
	// NOTE!!!: old XY order!!!
	template <class T> void XArray2DFFTRe<T>::FFTRe(XArray2D<std::complex<T> >& C, bool bForward, bool bCheck, bool bResizeToZero)
	{
		index_t i, j, nx, ny, nxd2, nyd2;

		if (bForward) 
		{ 
			nx = m_rXArray2D.GetDim1(); 
			ny = m_rXArray2D.GetDim2(); 
		}
		else 
		{ 
			nx = C.GetDim1(); 
			ny = C.GetDim2(); 
		}

		if (bCheck)
		{
			i = 2;
		
			while (i<nx) 
				i *= 2;
		
			if (i != nx)
				if (bForward) 
					throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFTRe<T>::FFTRe (m_dim1 is not a power of 2)");
				else 
					throw std::invalid_argument("invalid_argument 'C' in XArray2DFFTRe<T>::FFTRe (m_dim1 is not a power of 2)"); 
		
			j = 2;
		
			while (j<ny) 
				j *= 2;
		
			if (j != ny)
				if (bForward) 
					throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFTRe<T>::FFTRe (m_dim2 is not a power of 2)");
				else 
					throw std::invalid_argument("invalid_argument 'C' in XArray2DFFTRe<T>::FFTRe (m_dim2 is not a power of 2)"); 
		}

		// transform the head and store the result
		double wl, xlo, xhi, ylo, yhi;
		IXAHead* pHead = 0;
		IXAHWave2D* ph2 = 0;
		
		if (bForward) 
			pHead = m_rXArray2D.GetHeadPtr();
		else 
			pHead = C.GetHeadPtr();
		
		if (pHead) //if there is a head
		{
			ph2 = dynamic_cast<IXAHWave2D*>(pHead);;
			if (ph2) // if the head implements IXAHWave2D
			{
				wl = ph2->GetWl();
				if (bForward)
				{
					xlo = -0.5 / GetYStep(m_rXArray2D);
					xhi = (0.5 - 1.0 / m_rXArray2D.GetDim1()) / GetYStep(m_rXArray2D);
					ylo = -0.5 / GetXStep(m_rXArray2D);
					yhi = (0.5 - 1.0 / m_rXArray2D.GetDim2()) / GetXStep(m_rXArray2D);
				}
				else
				{
					xlo = -0.5 / GetYStep(C);
					xhi = (0.5 - 1.0 / C.GetDim1()) / GetYStep(C);
					ylo = -0.5 / GetXStep(C);
					yhi = (0.5 - 1.0 / C.GetDim2()) / GetXStep(C);
				}
			}
		}

		nxd2 = nx / 2;
		nyd2 = ny / 2;

		if (!bForward) 
			m_rXArray2D.Resize(nx, ny);
		else 
			Shuffle();

		T* arr = &(m_rXArray2D.front());
		vector<T> vecSpeq(nx*2);
		T* speq = &(vecSpeq[0]);

		if (!bForward)
		{
			//wrap C into m_rXArray2D and speq
			for (i=0; i<nxd2; i++) 
			{ 
				speq[nx+2*i] = std::real(C[i][0]);  
				speq[nx+2*i+1] = std::imag(C[i][0]); 
			}

			for (i=nxd2; i<nx; i++) 
			{ 
				speq[2*i-nx] = std::real(C[i][0]); 
				speq[2*i-nx+1] = std::imag(C[i][0]); 
			}
		
			for (i=0; i<nxd2; i++)
				for (j=nyd2; j<ny; j++)
				{ 
						m_rXArray2D[i+nxd2][2*j-ny] = std::real(C[i][j]); 
						m_rXArray2D[i+nxd2][2*j-ny+1] = std::imag(C[i][j]); 
				}

			for (i=nxd2; i<nx; i++)
				for (j=nyd2; j<ny; j++)
				{ 
					m_rXArray2D[i-nxd2][2*j-ny] = std::real(C[i][j]); 
					m_rXArray2D[i-nxd2][2*j-ny+1] = std::imag(C[i][j]); 
				}
			if (bResizeToZero) 
				C.Truncate(); //C is no longer needed
		}

		OouraFft<T> fft;

		fft.Real2D(arr, speq, nx, ny, 
			bForward? OouraFft<T>::eDirFwd : OouraFft<T>::eDirInv);

		if (!bForward) 
		{
			T anxy = T(2.0 / nx / ny);
			
			for (index_t k=0; k<nx*ny; k++) 
				arr[k] *= anxy;
			
			Shuffle();
			
			if (ph2)
			{
				IXAHWave2D* ph2a = CreateWavehead2D();
				ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
				m_rXArray2D.SetHeadPtr(ph2a);
			}
		}

		if (bForward)
		{
			// wrap m_rXArray2D and speq into C
			C.Resize(nx, ny);

			for (i=0; i<nxd2; i++) C[i][0] = std::complex<T>(speq[nx+2*i], speq[nx+2*i+1]);
			for (i=nxd2; i<nx; i++) C[i][0] = std::complex<T>(speq[2*i-nx], speq[2*i-nx+1]);
		
			for (i=0; i<nxd2; i++)
				for (j=nyd2; j<ny; j++)
					C[i][j] = std::complex<T>(m_rXArray2D[i+nxd2][2*j-ny], m_rXArray2D[i+nxd2][2*j-ny+1]);

			for (i=nxd2; i<nx; i++)
				for (j=nyd2; j<ny; j++)
					C[i][j] = std::complex<T>(m_rXArray2D[i-nxd2][2*j-ny], m_rXArray2D[i-nxd2][2*j-ny+1]);
	
				for (j=1; j<nyd2; j++) C[0][j] = std::conj(C[0][ny-j]);
		
			for (i=1; i<nxd2; i++)
				for (j=1; j<nyd2; j++)
					C[i][j] = std::conj(C[nx-i][ny-j]);
		
			for (i=nxd2; i<nx; i++)
				for (j=1; j<nyd2; j++)
					C[i][j] = std::conj(C[nx-i][ny-j]);

			if (bResizeToZero) m_rXArray2D.Truncate();

			if (ph2)
			{
				IXAHWave2D* ph2a = CreateWavehead2D();
				ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
				C.SetHeadPtr(ph2a);
			}
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::Convol
	//
	//	Convolves the 'wrapped' real XArray2D object with another XArray2D object 
	//
	/*!
		\brief		Convolves the 'wrapped' real XArray2D object with another XArray2D object 
		\param		rxa2DKer	The 'kernel' XArray2D to convolve with
		\param		bTrimBack	Determines if the 'wrapped' real XArray2D object is trimmed back to its original size (it may be padded inside this function)
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the kernel
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function convolves (in-place) the 'wrapped' real XArray2D object with a 'kernel' represented
			by another real XArray2D object. On output the 'wrapped' XArray2D object is replaced 
			by the convolution. The result can be optionally trimmed back to the original
			dimensions of the 'wrapped' XArray2D object. The convolution is performed using the 
			FFT method via the Ooura FFT library. rxa2DKer is always resized to zero
			on output (it is spoiled anyway). Zero padding for preventing aliasing should be done if
			required outside this function. Only the padding to the nearest power of 2 is done if 
			necessary inside this function. If the 'wrapped' object and rxa2DKer have different sizes,
			then the convolved array may be slightly shifted in x or y, because of the natural uncertainty
			of the initial position of rxa2DKer relative to the 'wrapped' object. Note also that the step sizes
			are ignored in this	implementation of the convolution (dx and dy are set to 1, even if non-trivial
			values are available from the Wavehead2D).
	*/	
	// NOTE!!!: old XY order!!!
	template <class T> void XArray2DFFTRe<T>::Convol(XArray2D<T>& rxa2DKer, bool bTrimBack)
	{
		if (&GetBaseObject() == &rxa2DKer)
			throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Convol (rxa2DKer cannot be the same as *this)"); 

		double wl, xlo, xhi, ylo, yhi;
		IXAHead* pHead = m_rXArray2D.GetHeadPtr();
		IXAHWave2D* ph2 = 0;
		if (pHead) // if there is a head
		{
			ph2 = dynamic_cast<IXAHWave2D*>(pHead);
			
			if (ph2) // if the head implements IXAHWave2D
			{
				wl = ph2->GetWl();
				xlo = ph2->GetYlo();
				xhi = ph2->GetYhi();
				ylo = ph2->GetXlo();
				yhi = ph2->GetXhi();
			}
		}
		
		double xst = GetYStep(m_rXArray2D);
		double yst = GetXStep(m_rXArray2D);

		IXAHead* pHeadK = rxa2DKer.GetHeadPtr();
		IXAHWave2D* ph2K;
		if (pHeadK) // if there is a head
		{
			ph2K = dynamic_cast<IXAHWave2D*>(pHeadK);
			
			if (ph2K) // if the head implements IXAHWave2D
			{
				if (fabs(wl - ph2K->GetWl()) > 0.01 * wl)
					throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Convol (different wavelength)");
				
				if (fabs(xst - GetYStep(rxa2DKer)) > 0.01 * xst || fabs(yst - GetXStep(rxa2DKer)) > 0.01 * yst)
					throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Convol (different step)");
			}
		}

		index_t nxA = m_rXArray2D.GetDim1(), nxKer = rxa2DKer.GetDim1();
		index_t nyA = m_rXArray2D.GetDim2(), nyKer = rxa2DKer.GetDim2();
		index_t nx = (nxA >= nxKer) ? nxA : nxKer; //nx is max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t ny = (nyA >= nyKer) ? nyA : nyKer; //ny is max(m_rXArray2D.ny, rxa2DKer.ny)

		index_t i = 2;
		
		while (i < nx) 
			i *= 2;
		
		nx = i; // nx is the smallest power of 2 not less than max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t j = 2;
		
		while (j < ny) 
			j *= 2;
		
		ny = j; // ny is the smallest power of 2 not less than max(m_rXArray2D.ny, rxa2DKer.ny)

		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, 0);
		XArray2DMove<T> tmp1(rxa2DKer);
		tmp1.Pad((nx-nxKer)/2, nx-nxKer-(nx-nxKer)/2, (ny-nyKer)/2, ny-nyKer-(ny-nyKer)/2, 0);

		index_t nxd2 = nx / 2;
		index_t nyd2 = ny / 2;

		T* arr = &(m_rXArray2D.front());
		vector<T> vecSpeq(nx*2);
		T* speq = &(vecSpeq[0]);
		T* arrK = &(rxa2DKer.front());
		vector<T> vecSpeqK(nx*2);
		T* speqK = &(vecSpeqK[0]);
	
		// forward FFTs
		OouraFft<T> fft;

		fft.Real2D(arr, speq, nx, ny, OouraFft<T>::eDirFwd);
		fft.Real2D(arrK, speqK, nx, ny, OouraFft<T>::eDirFwd);

		// multiplication of the FFTs (implementation from NumRecipes, 2nd ed., p.531)
		double factor = 2. / (nx * ny), re, im;
		T* sp1 = arr;
		T* sp2 = arrK;
		
		for (i=0; i<nx*ny/2; i++)
		{
			re = sp1[0] * sp2[0] - sp1[1] * sp2[1];
			im = sp1[0] * sp2[1] + sp1[1] * sp2[0];
			sp1[0] = T(factor * re);
			sp1[1] = T(factor * im);
			sp1 += 2; sp2 += 2;
		}
		
		sp1 = speq;
		sp2 = speqK;
		
		for (i=0; i<nx; i++)
		{
			re = sp1[0] * sp2[0] - sp1[1] * sp2[1];
			im = sp1[0] * sp2[1] + sp1[1] * sp2[0];
			sp1[0] = T(factor * re);
			sp1[1] = T(factor * im);
			sp1 += 2; sp2 += 2;
		}
	
		rxa2DKer.Truncate();

		// inverse FFT
		fft.Real2D(arr, speq, nx, ny, OouraFft<T>::eDirInv);

		Shuffle();

		if (bTrimBack) 
		{
			XArray2DMove<T> tmp2(m_rXArray2D);
			tmp2.Trim((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2);
		}

	#if(0) //@@@@@@@ I do not see why the following bit was considered necessary
		if (ph2)
		{
			IXAHWave2D* ph2a = CreateWavehead2D();
			index_t nx = m_rXArray2D.GetDim1();
			index_t ny = m_rXArray2D.GetDim2();
			long nx0 = long((nx - nxA) / 2);
			long nx1 = long(nx - nxA - (nx - nxA) / 2);
			long ny0 = long((ny - nyA) / 2);
			long ny1 = long(ny - nyA - (ny - nyA) / 2);
			xlo -= xst * nx0;
			xhi += xst * nx1;
			ylo -= yst * ny0;
			yhi += yst * ny1;
			ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
			m_rXArray2D.SetHeadPtr(ph2a);
		}
	#endif
	}

	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::Deconvol
	//
	//	Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method
	//
	/*!
		\brief		Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method
		\param		rxa2DKer	The 'kernel' XArray2D to deconvolve with
		\param		dblAlpha	Regularizaton parameter 
		\param		bRescale	Determines if after the reguralized deconvolution the image is rescaled to the theoretical L1 norm
		\param		bTrimBack	Determines if the 'wrapped' real XArray2D object is trimmed back to its original size (it may be padded inside this function)
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the kernel
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function performs regularized (Wiener) deconvolution of the 'wrapped' real XArray2D object
			with a 'kernel' represented by another real XArray2D object. On output the 'wrapped' 
			XArray2D object is replaced by the deconvolution. The result can be optionally trimmed 
			back to the original dimensions of the 'wrapped' XArray2D object. The deconvolution is 
			performed using the FFT method via the Ooura FFT library and a simple 
			Wiener-type regularization. rxa2DKer is always resized to zero on output (it is spoiled
			anyway). Zero padding for preventing aliasing should be done if required outside this
			function. Only the padding to the nearest power of 2 is done if necessary inside this
			function. If the 'wrapped' object and rxa2DKer have different sizes, then the 
			deconvolved array may be slightly shifted in x or y, because of the natural uncertainty
			of the initial position of rxa2DKer relative to the 'wrapped' object.  Note also that the step sizes
			are ignored in this	implementation of the deconvolution (dx and dy are set to 1, even if non-trivial
			values are available from the Wavehead2D).
	*/	
	// NOTE!!!: old XY order!!!
	template <class T> void XArray2DFFTRe<T>::Deconvol(XArray2D<T>& rxa2DKer, double dblAlpha, bool bRescale, bool bTrimBack/*, void* ppMklDftiDescriptor*/)
	{
		if (dblAlpha<0) 
			throw std::invalid_argument("invalid_argument 'dblAlpha' in XArray2DFFTRe<T>::Deconvol (regularization parameter cannot be negative)"); 

		if (&GetBaseObject() == &rxa2DKer)
			throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Deconvol (rxa2DKer cannot be the same as *this)"); 

		double wl, xlo, xhi, ylo, yhi;
		IXAHead* pHead = m_rXArray2D.GetHeadPtr();
		IXAHWave2D* ph2 = 0;
		if (pHead) // if there is a head
		{
			ph2 = dynamic_cast<IXAHWave2D*>(pHead);
			if (ph2) // if the head implements IXAHWave2D
			{
				wl = ph2->GetWl();
				xlo = ph2->GetYlo();
				xhi = ph2->GetYhi();
				ylo = ph2->GetXlo();
				yhi = ph2->GetXhi();
			}
		}
		double xst = GetYStep(m_rXArray2D);
		double yst = GetXStep(m_rXArray2D);

		IXAHead* pHeadK = rxa2DKer.GetHeadPtr();
		IXAHWave2D* ph2K;
		if (pHeadK) // if there is a head
		{
			ph2K = dynamic_cast<IXAHWave2D*>(pHeadK);
			if (ph2K) // if the head implements IXAHWave2D
			{
				if (fabs(wl - ph2K->GetWl()) > 0.01 * wl)
					throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Deconvol (different wavelength)");
				if (fabs(xst - GetYStep(rxa2DKer)) > 0.01 * xst || fabs(yst - GetXStep(rxa2DKer)) > 0.01 * yst)
					throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Deconvol (different step)");
			}
		}
	
		double temp = rxa2DKer.Norm(eNormL2);
		dblAlpha *= (temp * temp); //!!! renormalization of dblAlpha

		// if I>=0 and Ker>=0, then ||I*Ker||_L1 = ||I||_L1 * ||Ker||_L1.
		// we want to preserve this despite the regularization
		double dblDeconvNormL1;
		if (bRescale)
		{
			double dblKerNormL1 = rxa2DKer.Norm(eNormL1);
			if (dblKerNormL1 == 0)
				throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::Deconvol (rxa2DKer is identical zero)"); 
			double dblINormL1 = m_rXArray2D.Norm(eNormL1);
			if (dblINormL1 == 0) return; // flat image
			dblDeconvNormL1 =  dblINormL1 / dblKerNormL1;
		}

		index_t nxA = m_rXArray2D.GetDim1(), nxKer = rxa2DKer.GetDim1();
		index_t nyA = m_rXArray2D.GetDim2(), nyKer = rxa2DKer.GetDim2();
		index_t nx = (nxA >= nxKer) ? nxA : nxKer; // nx is max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t ny = (nyA >= nyKer) ? nyA : nyKer; // ny is max(m_rXArray2D.ny, rxa2DKer.ny)

		index_t i = 2;
		while (i < nx) i *= 2;
		nx = i; // nx is the smallest power of 2 not less than max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t j = 2;
		while (j < ny) j *= 2;
		ny = j; // ny is the smallest power of 2 not less than max(m_rXArray2D.ny, rxa2DKer.ny)

		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, 0);
		XArray2DMove<T> tmp1(rxa2DKer);
		tmp1.Pad((nx-nxKer)/2, nx-nxKer-(nx-nxKer)/2, (ny-nyKer)/2, ny-nyKer-(ny-nyKer)/2, 0);

		index_t nxd2 = nx / 2;
		index_t nyd2 = ny / 2;

	/*
		#ifdef XA_INTEL_MKL
			// This implementation is very inefficient for multiple deconvolutions of fixed-size arrays
			//!!! It actually works slower than the NumRecipes implementation below

			// create complex copies (real DFT is not much more efficient, while more complicated)
			long status;
			static long ldim[2];
			XArray2D<std::complex<T> > xKer;
			MakeComplex(rxa2DKer, T(0), xKer, false);
			rxa2DKer.Truncate();
			XArray2D<std::complex<T> > xArr;
			MakeComplex(m_rXArray2D, T(0), xArr, false);
			void* pxArr(&(xArr.front())), *pxKer(&(xKer.front()));
		
			static DFTI_DESCRIPTOR *my_desc1_handle(0);
			if (!my_desc1_handle || ldim[0] != xArr.GetDim1() || ldim[1] != xArr.GetDim2())
			{
				if (my_desc1_handle) DftiFreeDescriptor(&my_desc1_handle);
				ldim[0] = xArr.GetDim1(); ldim[1] = xArr.GetDim2();
				// the following two functions could be called only once for multiple deconvolutions
				if (DftiCreateDescriptor(&my_desc1_handle, GetValuetypeMKL(), DFTI_COMPLEX, 2, ldim))
					throw std::runtime_error("runtime_error in XArray2DFFTRe<T>::Deconvol (Intel MKL error)"); 
				if (DftiCommitDescriptor(my_desc1_handle))
					throw std::runtime_error("runtime_error in XArray2DFFTRe<T>::Deconvol (Intel MKL error)"); 

			}

			// forward FFTs
			status = DftiComputeForward(my_desc1_handle, pxArr);
			status += DftiComputeForward(my_desc1_handle, pxKer);
			if(status)
			{
				DftiFreeDescriptor(&my_desc1_handle);
				throw std::runtime_error("runtime_error in XArray2DFFTRe<T>::Deconvol (Intel MKL error)"); 
			}

			// division of the FFTs (implementation following NumRecipes, 2nd ed., p.531)
			std::complex<T> ctemp, ctemp1;
			std::complex<T>* cpxArr = &(xArr.front());
			std::complex<T>* cpxKer = &(xKer.front());
			for (i=0; i<nx*ny; i++)
			{
				ctemp = *(cpxKer + i);
				ctemp1 = std::conj(ctemp);
				*(cpxArr + i) *= ctemp1 / (dblAlpha + ctemp * ctemp1);
			}

			// inverse FFT
			status = DftiComputeBackward(my_desc1_handle, pxArr);
			if(status)
			{
				DftiFreeDescriptor(&my_desc1_handle);
				throw std::runtime_error("runtime_error in XArray2DFFTRe<T>::Deconvol (Intel MKL error)"); 
			}
		
			Re(xArr, m_rXArray2D);

			// DFTI descriptor may be passed to the calling function for deallocation at the right moment
			if (ppMklDftiDescriptor) *(DFTI_DESCRIPTOR**)ppMklDftiDescriptor = my_desc1_handle;
			else if(DftiFreeDescriptor(&my_desc1_handle)) throw std::runtime_error("runtime_error in XArray2DFFTRe<T>::Deconvol (Intel MKL error)"); 

		#else
	*/
			T* arr = &(m_rXArray2D.front());
			vector<T> vecSpeq(nx*2);
			T* speq = &(vecSpeq[0]);
			T* arrK = &(rxa2DKer.front());
			vector<T> vecSpeqK(nx*2);
			T* speqK = &(vecSpeqK[0]);

		
			// forward FFTs
			OouraFft<T> fft;

			fft.Real2D(arr, speq, nx, ny, OouraFft<T>::eDirFwd);
			fft.Real2D(arrK, speqK, nx, ny, OouraFft<T>::eDirFwd);

			// division of the FFTs (implementation following NumRecipes, 2nd ed., p.531)
			double factor = 2. / (nx * ny), factor1, re, im;
			T* sp1 = arr;
			T* sp2 = arrK;
			for (i=0; i<nx*ny/2; i++)
			{
				re = sp1[0] * sp2[0] + sp1[1] * sp2[1];
				im = -sp1[0] * sp2[1] + sp1[1] * sp2[0];
				factor1 = factor / (dblAlpha + sp2[0] * sp2[0] + sp2[1] * sp2[1]);
				sp1[0] = T(factor1 * re);
				sp1[1] = T(factor1 * im);
				sp1 += 2; sp2 += 2;
			}
			sp1 = speq;
			sp2 = speqK;
			for (i=0; i<nx; i++)
			{
				re = sp1[0] * sp2[0] + sp1[1] * sp2[1];
				im = -sp1[0] * sp2[1] + sp1[1] * sp2[0];
				factor1 = factor / (dblAlpha + sp2[0] * sp2[0] + sp2[1] * sp2[1]);
				sp1[0] = T(factor1 * re);
				sp1[1] = T(factor1 * im);
				sp1 += 2; sp2 += 2;
			}
			rxa2DKer.Truncate();

			// inverse FFT
			fft.Real2D(arr, speq, nx, ny, OouraFft<T>::eDirInv);

	//	#endif

		Shuffle();

		// if I>=0 and Ker>=0, then ||I||_L1 = ||I*Ker||_L1 / ||Ker||_L1.
		// we want to preserve this despite the regularization
		if (bRescale)
		{
			double dblINormL1 = m_rXArray2D.Norm(eNormL1);
			if (dblINormL1 != 0) m_rXArray2D *= T(dblDeconvNormL1 / dblINormL1);
		}

		if (bTrimBack) 
		{
			XArray2DMove<T> tmp2(m_rXArray2D);
			tmp2.Trim((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2);
		}

	#if(0) //@@@@@@@ I do not see why the following bit was considered necessary
		if (ph2)
		{
			IXAHWave2D* ph2a = CreateWavehead2D();
			index_t nx = m_rXArray2D.GetDim1();
			index_t ny = m_rXArray2D.GetDim2();
			long nx0 = long((nx - nxA) / 2);
			long nx1 = long(nx - nxA - (nx - nxA) / 2);
			long ny0 = long((ny - nyA) / 2);
			long ny1 = long(ny - nyA - (ny - nyA) / 2);
			xlo -= xst * nx0;
			xhi += xst * nx1;
			ylo -= yst * ny0;
			yhi += yst * ny1;
			ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
			m_rXArray2D.SetHeadPtr(ph2a);
		}
	#endif
	}

#if(0)
	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::DeconvolOpt
	//
	//	Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method with optimum alpha
	//
	/*!
		\brief		Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method with optimum alpha
		\param		rxa2DKer	The 'kernel' XArray2D to deconvolve with
		\param		dblSigma	Relative standard deviation of the noise corresponding to the average intensity in the image
		\param		bRescale	Determines if after the reguralized deconvolution the image is rescaled to the theoretical L1 norm
		\param		pdblAlphaOpt Pointer to return the found optimal value of the regularization parameter
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the kernel
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function performs regularized (Wiener) deconvolution of the 'wrapped' real XArray2D object
			with a 'kernel' represented by another real XArray2D object. On output the 'wrapped' 
			XArray2D object is replaced by the deconvolution. The deconvolution is 
			performed using the FFT method via the Ooura FFT library and a simple 
			Wiener-type regularization with the regularization parameter 'alpha' optimized using
			the simple routine 'goldenC' adopted from Numerical Recipes. The optimization minimizes the
			difference between the number of points in the image and the chi-square difference between
			the original blurred image and the reconvolution of the deconvolution with the current value
			of the regularization parameter. rxa2DKer is always resized to zero on output (it is spoiled
			anyway). Zero padding for preventing aliasing should be done if required outside this
			function. Only the padding to the nearest power of 2 is done if necessary inside this
			function. If the 'wrapped' object and rxa2DKer have different sizes, then the 
			deconvolved array may be slightly shifted in x or y, because of the natural uncertainty
			of the initial position of rxa2DKer relative to the 'wrapped' object.  Note also that the step sizes
			are ignored in this	implementation of the deconvolution (dx and dy are set to 1, even if non-trivial
			values are available from the Wavehead2D).
	*/	
	template <class T> void XArray2DFFTRe<T>::DeconvolOpt(XArray2D<T>& rxa2DKer, double dblSigma, bool bRescale, double* pdblAlphaOpt)
	{
		if (&GetBaseObject() == &rxa2DKer)
			throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::DeconvolOpt (rxa2DKer cannot be the same as *this)"); 

		if (dblSigma <= 0)
			throw std::invalid_argument("invalid_argument 'dblSigma' in XArray2DFFTRe<T>::DeconvolOpt (must be positive)"); 

		XArray2D<T> CopyOf_m_rXArray2D(m_rXArray2D); // keep an unspoiled copy of the input image
		m_vTemp.resize(2, 0);
		m_vTemp[0] = &CopyOf_m_rXArray2D; // use m_pTemp[0] to pass the image copy to CalcOnce
		m_vTemp[1] = &rxa2DKer; // use m_vTemp[1] to pass the rxa2DKer argument to CalcOnce

		std::vector<T> vecArgs(4);	// vector of parameters to be passed to 'CalcOnce' via 'goldenC'
		T* mass = &(vecArgs.front());
		mass[0] = T(vecArgs.size() - 1); // number of other parameters in this array
		mass[1] = T(0);	// base10-log of the regularization parameter
		mass[2] = T(bRescale);	// rescale? the deconvolved image
		mass[3] = T(dblSigma * dblSigma * m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2()); // targeted reconstruction accuracy

		index_t num_varied = 1;	// logAlpha is the argument (no.1) to be optimized by 'goldenC'
		T tLogAlpha = mass[1];	// a median LogAlpha
		T tLogAlphaMin = T(-6); // the minimum LogAlpha
		T tLogAlphaMax = T(4);	// the maximum LogAlpha
		T tTolerance = T(0.2);  // optimization tolerance
		T tAlphaOpt = 0;		// storage for the returned optimum alpha

		T tRecAccuracy = XA_Optimisation::GoldenMinimumSearch< XArray2DFFTRe<T> >(
			tLogAlphaMin, 
			tLogAlpha, 
			tLogAlphaMax,
			*this, 
			mass, 
			num_varied, 
			tTolerance, 
			&tAlphaOpt, 
			0);

		// do the deconvolution with the found optimum alpha
		m_rXArray2D = CopyOf_m_rXArray2D; // restore the blurred image
		Deconvol(rxa2DKer, pow(10.0, double(tAlphaOpt)), bRescale, true); // deconvolve
		*pdblAlphaOpt = pow(10.0, double(tAlphaOpt)); // return the optimal alpha
	}
#endif


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::DeconvolIter
	//
	//	Iteratively deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method
	//
	/*!
		\brief		Iteratively deconvolves the 'wrapped' real XArray2D object with another XArray2D object using regularized (Wiener) method
		\param		rxa2DKer	The 'kernel' XArray2D to deconvolve with
		\param		dblAlpha	Regularizaton parameter 
		\param		riNumIter	Maximum number of iterations as input, and the actual number of iterations as output
		\param		dblSigma	Relative standard deviation of the noise corresponding to the average intensity in the image
		\param		bRescale	Determines if after the reguralized deconvolution the image is rescaled to the theoretical L1 norm
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the kernel
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function performs iterative regularized (Wiener) deconvolution of the 'wrapped' real XArray2D object
			with a 'kernel' represented by another real XArray2D object. On output the 'wrapped' 
			XArray2D object is replaced by the deconvolution. The deconvolution is 
			performed using the FFT method via the Ooura FFT library and a simple 
			Wiener-type regularization with the regularization parameter 'alpha'. rxa2DKer is always resized to zero on output (it is spoiled
			anyway). Zero padding for preventing aliasing should be done if required outside this
			function. Only the padding to the nearest power of 2 is done if necessary inside this
			function. If the 'wrapped' object and rxa2DKer have different sizes, then the 
			deconvolved array may be slightly shifted in x or y, because of the natural uncertainty
			of the initial position of rxa2DKer relative to the 'wrapped' object.  Note also that the step sizes
			are ignored in this	implementation of the deconvolution (dx and dy are set to 1, even if non-trivial
			values are available from the Wavehead2D).
	*/	
	template <class T> void XArray2DFFTRe<T>::DeconvolIter(XArray2D<T>& rxa2DKer, double dblAlpha, double dblSigma, bool bRescale, index_t& riNumIter)
	{
		if (&GetBaseObject() == &rxa2DKer)
			throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::DeconvolIter (rxa2DKer cannot be the same as *this)"); 
		if (dblAlpha < 0)
			throw std::invalid_argument("invalid_argument 'dblAlpha' in XArray2DFFTRe<T>::DeconvolIter (cannot be negative)"); 
		if (riNumIter <= 0)
			throw std::invalid_argument("invalid_argument 'riNumIter' in XArray2DFFTRe<T>::DeconvolIter (must be positive)"); 
		if (dblSigma < 0)
			throw std::invalid_argument("invalid_argument 'dblSigma' in XArray2DFFTRe<T>::DeconvolIter (cannot be negtaive)"); 
	
		XArray2D<T> I0(m_rXArray2D), Ker0(rxa2DKer), Deconv(m_rXArray2D);
		double dblSigma2 = dblSigma * dblSigma * m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2();
		double dblImageMin = m_rXArray2D.Norm(eNormMin);
		double dblTolerance = 0.01; // here we always expect the chi-2 to be decreasing
		double norma_resid;

		Deconv *= T(0);
		for (index_t i = 0; i < riNumIter; i++)
		{
			if (i==0) Deconvol(rxa2DKer, dblAlpha, bRescale, true);
			else Deconvol(rxa2DKer, dblAlpha, false, true);
			rxa2DKer = Ker0; 
			Deconv += m_rXArray2D;
			m_rXArray2D = Deconv;
			if (i == riNumIter - 1) break;
			Convol(rxa2DKer, true);
			rxa2DKer = Ker0; 
			if (dblSigma > 0.0) // check if the residuals are comparable with the noise
			{
				if (dblImageMin > 0) // no zeros or negative values in the image
					norma_resid = I0.Chi2(m_rXArray2D, 1.0, true) / dblSigma2 - 1.0;
				else // non-Poisson statistics (sigma is the same at all points)
					norma_resid = I0.Chi2(m_rXArray2D, 1.0, false) / dblSigma2 - 1.0;
				if (norma_resid < dblTolerance) 
				{
					riNumIter = i + 1;
					break;
				}
			}
			m_rXArray2D -= I0;
			m_rXArray2D *= T(-1.0);
		}
		m_rXArray2D = Deconv;
	}


	// Auxilliary function called by the optimization routine 'goldenC'
	// What2Calc determines which member function is called for the caller optimization routine
	// ArrayPars is the array of parameters to be passed to a member function
	// ArrayPars[1] = log(dblAlpha)	// base10-log of the regularization parameter
	// ArrayPars[2] = bRescale		// rescale? the deconvolved image
	// ArrayPars[3] = dblEpsilon	// dblSigma * sqrt(m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2());
	template <class T> T XArray2DFFTRe<T>::CalcOnce(T* ArrayPars, index_t What2Calc)
	{
		switch (What2Calc)
		{
		case 0:
			{
			if (ArrayPars[0] != 3)
				throw std::invalid_argument("invalid_argument 'ArrayPars' in XArray2DFFTRe<T>::CalcOnce (wrong number of parameters)"); 
			
			m_rXArray2D = *(reinterpret_cast<XArray2D<T>*>(m_vTemp[0])); // restore the blurred image
			XArray2D<T> xaKer = *(reinterpret_cast<XArray2D<T>*>(m_vTemp[1])); // restore the convolution kernel
			Deconvol(xaKer, pow(10.0, double(ArrayPars[1])), ArrayPars[2] != 0, true); // deconvolve
			xaKer = *(reinterpret_cast<XArray2D<T>*>(m_vTemp[1])); // restore the convolution kernel
			Convol(xaKer, true); // reconvolve
			xaKer = *(reinterpret_cast<XArray2D<T>*>(m_vTemp[0])); // use xaKer as a temp storage for the blurred image
			double dblRet;
			
			if (xaKer.Norm(eNormMin) > 0) // no zeros or negative values in the image
				dblRet = fabs(xaKer.Chi2(m_rXArray2D, 1.0, true) / double(ArrayPars[3]) - 1.0);
			else // non-Poisson statistics (sigma is the same at all points)
				dblRet = fabs(xaKer.Chi2(m_rXArray2D, 1.0, false) / double(ArrayPars[3]) - 1.0);
			
			return T(dblRet);
			}
			break;
		default:
			throw std::invalid_argument("invalid_argument 'What2Calc' in XArray2DFFTRe<T>::CalcOnce (only 0 allowed)"); 
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::RichLucy
	//
	//	Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using Richardson-Lucy algorithm
	//
	/*!
		\brief		Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using Richardson-Lucy algorithm
		\param		rxa2DKer	The 'kernel' XArray2D to deconvolve with
		\param		dblSigma	Stoppage criterion (average noise-to-signal ratio in the input image)
		\param		riNumIter		Maximum number of iterations
		\param		bTrimBack	Determines if the 'wrapped' real XArray2D object is trimmed back to its original size (it may be padded inside this function)
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the kernel
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function performs Richardson-Lucy iterative deconvolution of the 'wrapped' real XArray2D object
			with a 'kernel' represented by another real XArray2D object. On output the 'wrapped' 
			XArray2D object is replaced by the deconvolution. The result can be optionally trimmed 
			back to the original dimensions of the 'wrapped' XArray2D object. The deconvolution is 
			performed using the FFTR methods of the XArray2D. rxa2DKer remains unchanged since only its copy
			is used. Zero padding for preventing aliasing should be done if required outside this
			function. Only the padding to the nearest power of 2 is done if necessary inside this
			function. If the 'wrapped' object and rxa2DKer have different sizes, then the 
			deconvolved array may be slightly shifted in x or y, because of the natural uncertainty
			of the initial position of rxa2DKer relative to the 'wrapped' object
	*/	
	// NOTE!!!: old XY order!!!
	template <class T> void XArray2DFFTRe<T>::RichLucy(XArray2D<T>& rxa2DKer, index_t& riNumIter, double dblSigma, bool bTrimBack)
	{
		if (riNumIter <= 0) 
			throw std::invalid_argument("invalid_argument 'riNumIter' in XArray2DFFTRe<T>::RichLucy (number of iterations must be positive)"); 

		if (dblSigma < 0)
			throw std::invalid_argument("invalid_argument 'dblSigma' in XArray2DFFTRe<T>::RichLucy (must be non-negative)"); 

		if (&GetBaseObject() == &rxa2DKer)
			throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::RichLucy (rxa2DKer cannot be the same as *this)"); 

	//	if (m_rXArray2D.Norm(xar::eNormMin) <= 0) 
	//		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFTRe<T>::RichLucy (raw image must be positive)"); 

		double wl, xlo, xhi, ylo, yhi;
		IXAHead* pHead = m_rXArray2D.GetHeadPtr();
		IXAHWave2D* ph2 = 0;
		if (pHead) // if there is a head
		{
			ph2 = dynamic_cast<IXAHWave2D*>(pHead);
			if (ph2) // if the head implements IXAHWave2D
			{
				wl = ph2->GetWl();
				xlo = ph2->GetYlo();
				xhi = ph2->GetYhi();
				ylo = ph2->GetXlo();
				yhi = ph2->GetXhi();
			}
		}
		double xst = GetYStep(m_rXArray2D);
		double yst = GetXStep(m_rXArray2D);

		IXAHead* pHeadK = rxa2DKer.GetHeadPtr();
		IXAHWave2D* ph2K;
		if (pHeadK) // if there is a head
		{
			ph2K = dynamic_cast<IXAHWave2D*>(pHeadK);
			if (ph2K) // if the head implements IXAHWave2D
			{
				if (fabs(wl - ph2K->GetWl()) > 0.01 * wl)
					throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::RichLucy (different wavelength)");
				if (fabs(xst - GetYStep(rxa2DKer)) > 0.01 * xst || fabs(yst - GetXStep(rxa2DKer)) > 0.01 * yst)
					throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::RichLucy (different step)");
			}
		}
	
		double dblKerNormL1 = rxa2DKer.Norm(eNormL1);
		
		if (dblKerNormL1 == 0)
			throw std::invalid_argument("invalid_argument 'rxa2DKer' in XArray2DFFTRe<T>::RichLucy (rxa2DKer is identical zero)"); 
		
		if (dblKerNormL1 != 1) 
		{
			rxa2DKer /= T(dblKerNormL1);
			m_rXArray2D /= T(dblKerNormL1);
		}
		
		double dblINormL1 = m_rXArray2D.Norm(eNormL1);
		
		if (dblINormL1 == 0) 
			return; // flat image

		index_t nxA = m_rXArray2D.GetDim1(), nxKer = rxa2DKer.GetDim1();
		index_t nyA = m_rXArray2D.GetDim2(), nyKer = rxa2DKer.GetDim2();
		index_t nx = (nxA >= nxKer) ? nxA : nxKer; // nx is max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t ny = (nyA >= nyKer) ? nyA : nyKer; // ny is max(m_rXArray2D.ny, rxa2DKer.ny)

		index_t i = 2;
		while (i<nx) 
			i *= 2;

		nx = i; // nx is the smallest power of 2 not less than max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t j = 2;
		
		while (j<ny) 
			j *= 2;
		
		ny = j; // ny is the smallest power of 2 not less than max(m_rXArray2D.ny, rxa2DKer.ny)

		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, 0);
		XArray2DMove<T> tmp1(rxa2DKer);
		tmp1.Pad((nx-nxKer)/2, nx-nxKer-(nx-nxKer)/2, (ny-nyKer)/2, ny-nyKer-(ny-nyKer)/2, 0);

		XArray2D<T> image(m_rXArray2D);
		image *= 0;
		image += T(dblINormL1 / (nxA * nyA)); // uniform initial guess

		XArray2D<T> reconv(m_rXArray2D);
		XArray2D<T> Kernel_copy(rxa2DKer);
		XArray2D<T> data_copy(m_rXArray2D);
		XArray2DFFTRe<T> FFTR(reconv);
		XArray2DFFTRe<T> Correl(data_copy);
	
		long Shift1, Shift2;

		double dblSigma2 = dblSigma * dblSigma * m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2();
		double dblTolerance = 0.01; // allowed inaccuracy of chi-2 fitting (we expect chi-2 to be decreasing)
		double dblImageMin = m_rXArray2D.Norm(eNormMin);
		double norma_resid;
		for (index_t iter = 0; iter < riNumIter; iter++)
		{
			reconv = image;
			Kernel_copy = rxa2DKer;
			FFTR.Convol(Kernel_copy); // after that reconv contains reconvolved image		
			data_copy = m_rXArray2D;
			for (i = 0; i < data_copy.GetDim1(); i++)
			{
				for (j = 0; j < data_copy.GetDim2(); j++)
				{
					if (reconv[i][j] == 0) 
					{
						data_copy[i][j] = 0;
					}
					else
					{
						data_copy[i][j] /= reconv[i][j];
					}
				}
			}
			if (dblSigma > 0.0) // check if the residuals are comparable with the noise
			{
				if (dblImageMin > 0) // no zeros or negative values in the image
					norma_resid = m_rXArray2D.Chi2(reconv, 1.0, true) / dblSigma2 - 1.0;
				else // non-Poisson statistics (sigma is the same at all points)
					norma_resid = m_rXArray2D.Chi2(reconv, 1.0, false) / dblSigma2 - 1.0;
				if (norma_resid < dblTolerance)
				{
					riNumIter = iter + 1;
					break;
				}
			}
			Kernel_copy = rxa2DKer;
			double correl_max = Correl.CorrelateOld(Kernel_copy, Shift1, Shift2, false, false);
			image *= data_copy;
		}
		m_rXArray2D = image;

		if (bTrimBack) 
		{
			XArray2DMove<T> tmp2(m_rXArray2D);
			tmp2.Trim((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2);
		}

	#if(0) //@@@@@@@ I do not see why the following bit was considered necessary
		if (ph2)
		{
			IXAHWave2D* ph2a = CreateWavehead2D();
			index_t nx = m_rXArray2D.GetDim1();
			index_t ny = m_rXArray2D.GetDim2();
			long nx0 = long((nx - nxA) / 2);
			long nx1 = long(nx - nxA - (nx - nxA) / 2);
			long ny0 = long((ny - nyA) / 2);
			long ny1 = long(ny - nyA - (ny - nyA) / 2);
			xlo -= xst * nx0;
			xhi += xst * nx1;
			ylo -= yst * ny0;
			yhi += yst * ny1;
			ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
			m_rXArray2D.SetHeadPtr(ph2a);
		}
	#endif
	}

	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::DeconvolTaylor
	//
	//		Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using Taylor series
	//
	/*!
		\brief		Deconvolves the 'wrapped' real XArray2D object with another XArray2D object using Taylor series
		\param		rxa2DKer is the PSF to be deconvolved with
		\param		iMaxDerOrder is the maximum order of the Taylor series
		\param		dblSigmaY	is the half-width in pixels of the Gaussian filter along the first dimension
		\param		dblSigmaX	is the half-width in pixels of the Gaussian filter along the second dimension
		\param		dblRst_dev_av is the average noise-to-signal ratio in the input image
		\param		riNumIter contains the maximum number of iterations on input, and the actual number of iterations on output
		\exception	std::invalid_argument is thrown if the input parameters are inconsistent
		\return		\a None
		\par		Description:
			This function deconvolves the 'wrapped' image using its local 2D Taylor series decomposition and the integral
			moments of the PSF (convolution kernel). The maximum order of the Taylor series determines the 'degree' of 
			deconvolution as well as noise sensitivity. Here the 'zero order' deconvolution coincides with the original
			convoluted image, first order deconovolution depends on the 1st order partial derivatives of the image, etc.
			Gaussian (low-pass) filtering with defined parameters (iSigmaY and iSigmaX) is applied to the image prior to 
			calculation of partial derivatives. Deconvolution can be performed iteratively, applying the algorithm to the
			difference between the original input image and the convolution of the current solution with the PSF. If that
			difference is smaller than the average standard deviation of the noise (defined via the 'dblRst_dev_av' parameter),
			the iteration process is interrupted.
	*/	
	template<class T> void XArray2DFFTRe<T>::DeconvolTaylor(XArray2D<T>& rxa2DKer, index_t iMaxDerOrder, double dblSigmaY, double dblSigmaX,  double dblRst_dev_av, index_t& riNumIter)
	{
		if (iMaxDerOrder <= 0)
			throw std::invalid_argument("invalid arguments 'iMaxDerOrder' in XArray2DFFTRe<T>::DeconvolTaylor (must be positive)");
		
		if (dblSigmaY < 0 || dblSigmaX < 0)
			throw std::invalid_argument("invalid arguments 'dblSigmaX and dblSigmaY' in XArray2DFFTRe<T>::DeconvolTaylor (must be non-negative)");
		
		if (riNumIter <= 0)
			throw std::invalid_argument("invalid arguments 'riNumIter' in XArray2DFFTRe<T>::DeconvolTaylor (must be positive)");
		
		if (dblRst_dev_av < 0)
			throw std::invalid_argument("invalid arguments 'dblRst_dev_av' in XArray2DFFTRe<T>::DeconvolTaylor (must be non-negative)");
		
		if (fabs(GetXStep(m_rXArray2D) - GetXStep(rxa2DKer)) > 1.e-3 * GetXStep(m_rXArray2D) ||
			fabs(GetYStep(m_rXArray2D) - GetYStep(rxa2DKer)) > 1.e-3 * GetYStep(m_rXArray2D))
			throw std::invalid_argument("invalid arguments 'm_rXArray2D and rxa2DKer' in XArray2DFFTRe<T>::DeconvolTaylor (different steps)");

		index_t nxA = m_rXArray2D.GetDim1(), nxKer = rxa2DKer.GetDim1();
		index_t nyA = m_rXArray2D.GetDim2(), nyKer = rxa2DKer.GetDim2();
		index_t nx = (nxA>=nxKer) ? nxA : nxKer; // nx is max(m_rXArray2D.nx, rxa2DKer.nx)
		index_t ny = (nyA>=nyKer) ? nyA : nyKer; // ny is max(m_rXArray2D.ny, rxa2DKer.ny)

		{
			index_t i = 2;
			while (i<nx) i *= 2;
			nx = i; // nx is the smallest power of 2 not less than max(m_rXArray2D.nx, rxa2DKer.nx)
			index_t j = 2;
			while (j<ny) j *= 2;
			ny = j; // ny is the smallest power of 2 not less than max(m_rXArray2D.ny, rxa2DKer.ny)
		}

		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, 0);
		XArray2DMove<T> tmp1(rxa2DKer);
		tmp1.Pad((nx-nxKer)/2, nx-nxKer-(nx-nxKer)/2, (ny-nyKer)/2, ny-nyKer-(ny-nyKer)/2, 0);

		XArray2D<T> xaCoeff, xbCoeff;
		XArray2DMove<T> KerMove(rxa2DKer);
		KerMove.PSFMoments(iMaxDerOrder, xaCoeff);
		KerMove.InverseSeries(xaCoeff, xbCoeff);
		xbCoeff *= T(GetYStep(rxa2DKer) * GetXStep(rxa2DKer)); // required because XArray convolution ignores step sizes
		XArray2D<T> xaTemp, m_rXArray2D0(m_rXArray2D), xaOut(m_rXArray2D), rxa2DKer0(rxa2DKer);
		XArray2DMove<T> TempMove(xaTemp);
		XArray2DFFTRe TempFFTRe(xaTemp);
		
		double dblSigma2 = dblRst_dev_av * dblRst_dev_av * m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2();
		double dblErr, dblErr0 = (double)(m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2());
		double dblImageMin = m_rXArray2D0.Norm(eNormMin);
		double dblTolerance = 0.01; // allowed inaccuracy of chi-2 fitting (we assume that chi-2 is monotonically decreasing)

					//char aaa[128]; //@@@@@@@ for testing only

		xaOut *= T(0);	
		for (index_t k = 0; k < riNumIter; k++)
		{
			for (index_t i = 0; i <= iMaxDerOrder; i++)
			{
				for (index_t j = 0; j <= iMaxDerOrder - i; j++)
				{
					if (fabs(xbCoeff[i][j]) < 1.e-3 * fabs(xbCoeff[0][0])) continue; // save computer time
					xaTemp = m_rXArray2D;
					if (dblSigmaY + i > 0 || dblSigmaX + j > 0) TempFFTRe.FilterGauss(dblSigmaY + i, dblSigmaX + j);
					TempMove.Derivative(i, j);
					xaTemp *= xbCoeff[i][j];
					xaOut += xaTemp; // add to the accumulated solution
				}
			}
			if (k == riNumIter - 1) break;
			xaTemp = xaOut; // create a temporary to check the current accumulated solution
			TempFFTRe.Convol(rxa2DKer, true);
			rxa2DKer = rxa2DKer0; // restore rxa2DKer
			if (dblSigma2 > 0.0) // check if the residuals are comparable with the noise
			{
				if (dblImageMin > 0) // no zeros or negative values in the image
					dblErr = m_rXArray2D0.Chi2(xaTemp, 1.0, true) / dblSigma2 - 1.0;
				else // non-Poisson statistics (sigma is the same at all points)
					dblErr = m_rXArray2D0.Chi2(xaTemp, 1.0, false) / dblSigma2 - 1.0;
				// note that an additional stopping criterion, dblErr0 < dblErr0, sometimes helps
				if (dblErr < dblTolerance || fabs(dblErr) > fabs(dblErr0)) 
				{
					riNumIter = k + 1;
					break;
				}
				dblErr0 = dblErr;
			}
			xaTemp -= m_rXArray2D0; // (-1) * current_residual
			xaTemp *= T(-1);
			m_rXArray2D = xaTemp; // new_rhs
		}
	
		// replace the edge pixels spoiled by derivatives
		index_t iMask = index_t((float)iMaxDerOrder / 2.0f + 0.6);
		XArray2DMove<T> OutMove(xaOut);
		OutMove.Trim(iMask, iMask, iMask, iMask);
		OutMove.PadMirror(iMask, iMask, iMask, iMask);
		m_rXArray2D = xaOut;

		// trim back
		XArray2DMove<T> tmp2(m_rXArray2D);
		tmp2.Trim((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2);
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::Correlate
	//
	//	Correlates the 'wrapped' real XArray2D object with another XArray2D object 
	//
	/*!
		\brief		Phase-correlates the 'wrapped' real XArray2D object with another XArray2D object 
		\param		rxa2DOther	Reference XArray2D object to correlate with
		\param		rlngShift1	Detected shift along the first dimension
		\param		rlngShift2	Detected shift along the second dimension
		\param		bTruncateThis	Determines if the 'wrapped' real XArray2D object is resized to zero
		\param		bTruncateRefArr	Determines if the reference XArray2D object is resized to zero
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the reference
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function phase-correlates (in-place) the 'wrapped' real XArray2D object with a 'reference' represented
			by another real XArray2D object. Detected 'optimal' shift ('wrapped' with respect to the reference) is
			returned in (rlngShift1, rlngShift2). If bTruncateRefArr is true (default) rxa2DOther is resized to
			zero on output (it is spoiled anyway). If bTruncateThis is true (default),  the 'wrapped' real XArray2D
			object is also resized to zero on output, otherwise  the 'wrapped' real XArray2D object contains the
			phase-correlation matrix (possibly padded). The correlation is performed using the FFT method via the 
			Ooura FFT library.
	*/	
	// NOTE!!!: old XY order!!!
	template <class T> T XArray2DFFTRe<T>::Correlate(XArray2D<T>& rxa2DOther, long& rlngShift1, long& rlngShift2, bool bTruncateThis, bool bTruncateRefArr)
	{
		if (&GetBaseObject() == &rxa2DOther)
			throw std::invalid_argument("invalid_argument 'RefArr' in XArray2DFFTRe<T>::Correlate (cannot correlate array with itself)"); 

		double wl, xlo, xhi, ylo, yhi;
		IXAHead* pHead = m_rXArray2D.GetHeadPtr();
		IXAHWave2D* ph2 = 0;
		if (pHead) // if there is a head
		{
			ph2 = dynamic_cast<IXAHWave2D*>(pHead);
			if (ph2) // if the head implements IXAHWave2D
			{
				wl = ph2->GetWl();
				xlo = ph2->GetYlo();
				xhi = ph2->GetYhi();
				ylo = ph2->GetXlo();
				yhi = ph2->GetXhi();
			}
		}
		double xst = GetYStep(m_rXArray2D);
		double yst = GetXStep(m_rXArray2D);

		IXAHead* pHeadK = rxa2DOther.GetHeadPtr();
		IXAHWave2D* ph2K;
		if (pHeadK) // if there is a head
		{
			ph2K = dynamic_cast<IXAHWave2D*>(pHeadK);
			if (ph2K) // if the head implements IXAHWave2D
			{
				if (fabs(wl - ph2K->GetWl()) > 0.01 * wl)
					throw std::invalid_argument("invalid_argument 'rxa2DOther' in XArray2DFFTRe<T>::Correlate (different wavelength)");
				if (fabs(xst - GetYStep(rxa2DOther)) > 0.01 * xst || fabs(yst - GetXStep(rxa2DOther)) > 0.01 * yst)
					throw std::invalid_argument("invalid_argument 'rxa2DOther' in XArray2DFFTRe<T>::Correlate (different step)");
			}
		}

		index_t nxA = m_rXArray2D.GetDim1(), nyA = m_rXArray2D.GetDim2(), nx = nxA, ny = nyA;

		if (nxA != rxa2DOther.GetDim1() || nyA != rxa2DOther.GetDim2())
			throw std::invalid_argument("invalid_argument 'RefArr' in XArray2DFFTRe<T>::Correlate (array dimensions are different)"); 

		index_t i = 2;
		while (i<nx) i *= 2;
		nx = i; // nx is the smallest power of 2 not less than m_rXArray2D.nx
		index_t j = 2;
		while (j<ny) j *= 2;
		ny = j; // ny is the smallest power of 2 not less than m_rXArray2D.ny

		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, T(m_rXArray2D.Norm(eNormAver)));
		XArray2DMove<T> tmp1(rxa2DOther);
		tmp1.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, T(rxa2DOther.Norm(eNormAver)));

		double norm1 = m_rXArray2D.Norm(eNormL2);
		double norm2 = rxa2DOther.Norm(eNormL2);
		double norm = norm1 * norm2;

		index_t nxd2 = nx / 2;
		index_t nyd2 = ny / 2;

		T* arrA = &(m_rXArray2D.front());
		vector<T> vecSpeqA(nx*2);
		T* speqA = &(vecSpeqA[0]);
		T* arrB = &(rxa2DOther.front());
		vector<T> vecSpeqB(nx*2);
		T* speqB = &(vecSpeqB[0]);

		// forward FFTs
		OouraFft<T> fft;

		fft.Real2D(arrA, speqA, nx, ny, OouraFft<T>::eDirFwd);
		fft.Real2D(arrB, speqB, nx, ny, OouraFft<T>::eDirFwd);

		// multiplication of the FFTs (implementation following NumRecipes, 2nd ed., p.531)
		double factor = 2. / (nx * ny), re, im, abs12;
		T* sp1 = arrA;
		T* sp2 = arrB;
		for (i=0; i<nx*ny/2; i++)
		{
			re = sp1[0] * sp2[0] + sp1[1] * sp2[1];
			im = -sp1[0] * sp2[1] + sp1[1] * sp2[0];
			abs12 = sqrt(re * re + im * im);
			abs12 = (abs12 == 0.0) ? 0.0 : (1.0 / abs12);
			sp1[0] = T(factor * re * abs12);
			sp1[1] = T(factor * im * abs12);
			sp1 += 2; sp2 += 2;
		}
		sp1 = speqA;
		sp2 = speqB;
		for (i=0; i<nx; i++)
		{
			re = sp1[0] * sp2[0] + sp1[1] * sp2[1];
			im = -sp1[0] * sp2[1] + sp1[1] * sp2[0];
			abs12 = sqrt(re * re + im * im);
			abs12 = (abs12 == 0.0) ? 0.0 : (1.0 / abs12);
			sp1[0] = T(factor * re * abs12);
			sp1[1] = T(factor * im * abs12);
			sp1 += 2; sp2 += 2;
		}
		if (bTruncateRefArr) rxa2DOther.Truncate();

		// inverse FFT
		fft.Real2D(arrA, speqA, nx, ny, OouraFft<T>::eDirInv);

		Shuffle();

		// finding optimal shift
		T uamax = m_rXArray2D[0][0]; 
		rlngShift1 = 0; rlngShift2 = 0;
		for (i=0; i<nx; i++)
			for (j=0; j<ny; j++)
				if (m_rXArray2D[i][j]>uamax) { rlngShift1 = long(nxd2 - i); rlngShift2 = long(nyd2 - j); uamax = m_rXArray2D[i][j]; }
		//uamax /= T(norm);

		if (bTruncateThis) m_rXArray2D.Truncate();
	#if(0) // @@@@@@@@@@@ I don't see why the following bit was considered necessary
		else if (ph2)
		{
			IXAHWave2D* ph2a = CreateWavehead2D();
			index_t nx = m_rXArray2D.GetDim1();
			index_t ny = m_rXArray2D.GetDim2();
			long nx0 = long((nx - nxA) / 2);
			long nx1 = long(nx - nxA - (nx - nxA) / 2);
			long ny0 = long((ny - nyA) / 2);
			long ny1 = long(ny - nyA - (ny - nyA) / 2);
			xlo -= xst * nx0;
			xhi += xst * nx1;
			ylo -= yst * ny0;
			yhi += yst * ny1;
			ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
			m_rXArray2D.SetHeadPtr(ph2a);
		}
	#endif

		return uamax;
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::CorrelateOld
	//
	//	Correlates the 'wrapped' real XArray2D object with another XArray2D object 
	//
	/*!
		\brief		Correlates the 'wrapped' real XArray2D object with another XArray2D object 
		\param		rxa2DOther	Reference XArray2D object to correlate with
		\param		rlngShift1	Detected shift along the first dimension
		\param		rlngShift2	Detected shift along the second dimension
		\param		bTruncateThis	Determines if the 'wrapped' real XArray2D object is resized to zero
		\param		bTruncateRefArr	Determines if the reference XArray2D object is resized to zero
		\exception  std::invalid_argument is thrown if the 'wrapped' real XArray2D object and the reference
					object coincide, or if they have different wavelengths or step sizes
		\exception  std::runtime_error is thrown if there is not enough memory
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function correlates (in-place) the 'wrapped' real XArray2D object with a 'reference' represented
			by another real XArray2D object. Detected 'optimal' shift ('wrapped' with respect to the reference) is
			returned in (rlngShift1, rlngShift2). If bTruncateRefArr is true (default) rxa2DOther is resized to
			zero on output (it is spoiled anyway). If bTruncateThis is true (default),  the 'wrapped' real XArray2D
			object is also resized to zero on output, otherwise  the 'wrapped' real XArray2D object contains the
			phase-correlation matrix (possibly padded). The correlation is performed using the FFT method via the 
			Ooura FFT library.
	*/	
	// NOTE!!!: old XY order!!!
	template <class T> T XArray2DFFTRe<T>::CorrelateOld(XArray2D<T>& rxa2DOther, long& rlngShift1, long& rlngShift2, bool bTruncateThis, bool bTruncateRefArr)
	{
		if (&GetBaseObject() == &rxa2DOther)
			throw std::invalid_argument("invalid_argument 'RefArr' in XArray2DFFTRe<T>::CorrelateOld (cannot correlate array with itself)"); 

		double wl, xlo, xhi, ylo, yhi;
		IXAHead* pHead = m_rXArray2D.GetHeadPtr();
		IXAHWave2D* ph2 = 0;
		if (pHead) // if there is a head
		{
			ph2 = dynamic_cast<IXAHWave2D*>(pHead);
			if (ph2) // if the head implements IXAHWave2D
			{
				wl = ph2->GetWl();
				xlo = ph2->GetYlo();
				xhi = ph2->GetYhi();
				ylo = ph2->GetXlo();
				yhi = ph2->GetXhi();
			}
		}
		double xst = GetYStep(m_rXArray2D);
		double yst = GetXStep(m_rXArray2D);

		IXAHead* pHeadK = rxa2DOther.GetHeadPtr();
		IXAHWave2D* ph2K;
		if (pHeadK) // if there is a head
		{
			ph2K = dynamic_cast<IXAHWave2D*>(pHeadK);
			if (ph2K) // if the head implements IXAHWave2D
			{
				if (fabs(wl - ph2K->GetWl()) > 0.01 * wl)
					throw std::invalid_argument("invalid_argument 'rxa2DOther' in XArray2DFFTRe<T>::CorrelateOld (different wavelength)");
				if (fabs(xst - GetYStep(rxa2DOther)) > 0.01 * xst || fabs(yst - GetXStep(rxa2DOther)) > 0.01 * yst)
					throw std::invalid_argument("invalid_argument 'rxa2DOther' in XArray2DFFTRe<T>::CorrelateOld (different step)");
			}
		}

		index_t nxA = m_rXArray2D.GetDim1(), nyA = m_rXArray2D.GetDim2(), nx = nxA, ny = nyA;

		if (nxA != rxa2DOther.GetDim1() || nyA != rxa2DOther.GetDim2())
			throw std::invalid_argument("invalid_argument 'RefArr' in XArray2DFFTRe<T>::CorrelateOld (array dimensions are different)"); 

		index_t i = 2;
		while (i<nx) i *= 2;
		nx = i; // nx is the smallest power of 2 not less than m_rXArray2D.nx
		index_t j = 2;
		while (j<ny) j *= 2;
		ny = j; // ny is the smallest power of 2 not less than m_rXArray2D.ny

		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, T(m_rXArray2D.Norm(eNormAver)));
		XArray2DMove<T> tmp1(rxa2DOther);
		tmp1.Pad((nx-nxA)/2, nx-nxA-(nx-nxA)/2, (ny-nyA)/2, ny-nyA-(ny-nyA)/2, T(rxa2DOther.Norm(eNormAver)));

		double norm1 = m_rXArray2D.Norm(eNormL2);
		double norm2 = rxa2DOther.Norm(eNormL2);
		double norm = norm1 * norm2;

		index_t nxd2 = nx / 2;
		index_t nyd2 = ny / 2;

		T* arrA = &(m_rXArray2D.front());
		vector<T> vecSpeqA(nx*2);
		T* speqA = &(vecSpeqA[0]);
		T* arrB = &(rxa2DOther.front());
		vector<T> vecSpeqB(nx*2);
		T* speqB = &(vecSpeqB[0]);

		// forward FFTs
		OouraFft<T> fft;

		fft.Real2D(arrA, speqA, nx, ny, OouraFft<T>::eDirFwd);
		fft.Real2D(arrB, speqB, nx, ny, OouraFft<T>::eDirFwd);

		// multiplication of the FFTs (implementation following NumRecipes, 2nd ed., p.531)
		double factor = 2. / (nx * ny), re, im;
		T* sp1 = arrA;
		T* sp2 = arrB;
		
		for (i=0; i<nx*ny/2; i++)
		{
			re = sp1[0] * sp2[0] + sp1[1] * sp2[1];
			im = -sp1[0] * sp2[1] + sp1[1] * sp2[0];
			sp1[0] = T(factor * re);
			sp1[1] = T(factor * im);
			sp1 += 2; sp2 += 2;
		}
		sp1 = speqA;
		sp2 = speqB;
		
		for (i=0; i<nx; i++)
		{
			re = sp1[0] * sp2[0] + sp1[1] * sp2[1];
			im = -sp1[0] * sp2[1] + sp1[1] * sp2[0];
			sp1[0] = T(factor * re);
			sp1[1] = T(factor * im);
			sp1 += 2; sp2 += 2;
		}
		
		if (bTruncateRefArr) 
			rxa2DOther.Truncate();

		// inverse FFT
		fft.Real2D(arrA, speqA, nx, ny, OouraFft<T>::eDirInv);

		Shuffle();

		// finding optimal shift
		T uamax = m_rXArray2D[0][0]; 
		rlngShift1 = 0; rlngShift2 = 0;
		
		for (i=0; i<nx; i++)
			for (j=0; j<ny; j++)
				if (m_rXArray2D[i][j]>uamax) { rlngShift1 = long(nxd2 - i); rlngShift2 = long(nyd2 - j); uamax = m_rXArray2D[i][j]; }
		//uamax /= T(norm);

		if (bTruncateThis) m_rXArray2D.Truncate();
	#if(0) // @@@@@@@@@@@ I don't see why the following bit was considered necessary
		else if (ph2)
		{
			IXAHWave2D* ph2a = CreateWavehead2D();
			index_t nx = m_rXArray2D.GetDim1();
			index_t ny = m_rXArray2D.GetDim2();
			long nx0 = long((nx - nxA) / 2);
			long nx1 = long(nx - nxA - (nx - nxA) / 2);
			long ny0 = long((ny - nyA) / 2);
			long ny1 = long(ny - nyA - (ny - nyA) / 2);
			xlo -= xst * nx0;
			xhi += xst * nx1;
			ylo -= yst * ny0;
			yhi += yst * ny1;
			ph2a->SetData(wl, xlo, xhi, ylo, yhi); // x-y order reversed
			m_rXArray2D.SetHeadPtr(ph2a);
		}
	#endif

		return uamax;
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FilterRect
	//
	// Filters the 'wrapped' real XArray2D object by convoluting it with a rectangular kernel
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by convoluting it with a rectangular kernel
		\param		iYWidth	Width in pixels of the kernel along the first dimension
		\param		iXWidth	Width in pixels of the kernel along the second dimension
		\exception  std::invalid_argument is thrown if iYWidth or iXWidth are larger than the corresponding size of the 'wrapped' object
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements box-averaging filter in direct space. The filtering is performed
			using the FFT method. This implementation is simple and flexible, but it is suboptimal in
			terms of speed and RAM
	*/	
	template <class T> void XArray2DFFTRe<T>::FilterRect(index_t iYWidth, index_t iXWidth)
	{
		if (iYWidth > m_rXArray2D.GetDim1() || iXWidth > m_rXArray2D.GetDim2())
			throw std::invalid_argument("invalid_argument 'iYWidth or iXWidth' in XArray2DFFTRe<T>::FilterRect (filter width too large)");
		
#ifdef HAVE_UNIQUE_PTR
		std::unique_ptr<IXAHead> pHeadAuto;
#else
		std::unique_ptr<IXAHead> pHeadAuto;
#endif

		IXAHead* pHead( m_rXArray2D.GetHeadPtr() );

		if( pHead )
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif

		XArray2D<std::complex<T> > C, Cbox;
		XArray2D<T> box(GetBaseObject());
		box *= T(0);
		index_t ny2 = GetBaseObject().GetDim1() / 2;
		index_t nx2 = GetBaseObject().GetDim2() / 2;
		XArray2DMove<T> temp(box);
		temp.FillRect(ny2 - iYWidth / 2, ny2 - iYWidth / 2 + iYWidth, nx2 - iXWidth / 2, nx2 - iXWidth / 2 + iXWidth, T(1.0 / (iYWidth * iXWidth)));
		XArray2DFFTRe<T> Box(box);
		Box.FFTRe(Cbox, true);
		FFTRe(C, true);
		C *= Cbox;
		FFTRe(C, false);
		
		if (pHead) 
			m_rXArray2D.SetHeadPtr(pHeadAuto->Clone()); // restore the head
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FilterEpanechnikov
	//
	// Filters the 'wrapped' real XArray2D object by convoluting it with an Epanechnikov kernel
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by convoluting it with an Epanechnikov kernel
		\param		dblWidthY	width in pixels of the kernel along the first dimension
		\param		dblWidthX	width in pixels of the kernel along the second dimension
		\exception  std::invalid_argument is thrown if dblWidthY or dblWidthX are larger than the corresponding
					size of the 'wrapped' object or are non-positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements Epanechnikov filter in direct space. The filtering is performed
			using the FFT method. This implementation is simple and flexible, but it is suboptimal in
			terms of speed and RAM
	*/	
	template <class T> void XArray2DFFTRe<T>::FilterEpanechnikov(double dblWidthY, double dblWidthX)
	{
		if (dblWidthY > m_rXArray2D.GetDim1() || dblWidthX > m_rXArray2D.GetDim2())
			throw std::invalid_argument("invalid_argument 'dblWidthY or dblWidthX' in XArray2DFFTRe<T>::FilterEpanechnikov (filter width too large)");
		if (dblWidthY <= 0. || dblWidthX <= 0.)
			throw std::invalid_argument("invalid_argument 'dblWidthY or dblWidthX' in XArray2DFFTRe<T>::FilterEpanechnikov (filter width is non-positive)");
		
#ifdef HAVE_UNIQUE_PTR
		std::unique_ptr<IXAHead> pHeadAuto;
#else
		std::unique_ptr<IXAHead> pHeadAuto;
#endif

		IXAHead* pHead( m_rXArray2D.GetHeadPtr() );

		if( pHead )
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif
		XArray2D<T> Kernel(GetBaseObject());
		Kernel.Fill(1.);
		XArray2DMove<T> Move(Kernel);
		Move.EpanechnikovEnvelope(dblWidthY, dblWidthX);
		Convol(Kernel, true);
		
		if (pHead) 
			m_rXArray2D.SetHeadPtr(pHeadAuto->Clone()); // restore the head
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FilterSymEpanechnikov
	//
	// Filters the 'wrapped' real XArray2D object by convoluting it with an axially symmetrical Epanechnikov kernel
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by convoluting it with an axially symmetrical Epanechnikov kernel
		\param		dblWidth	diameter in pixels of the kernel
		\exception  std::invalid_argument is thrown if dblWidth is larger than the dimenasions
					of the 'wrapped' object or is non-positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements an axially symmetrical Epanechnikov filter in direct space. The filtering is performed
			using the FFT method. This implementation is simple and flexible, but it is suboptimal in
			terms of speed and RAM
	*/	
	template <class T> void XArray2DFFTRe<T>::FilterSymEpanechnikov( double dblWidth )
	{
		if( dblWidth > m_rXArray2D.GetDim1() || dblWidth > m_rXArray2D.GetDim2() )
			throw std::invalid_argument("invalid_argument 'dblWidth' in XArray2DFFTRe<T>::FilterSymEpanechnikov (filter width too large)");
		if( dblWidth <= 0. )
			throw std::invalid_argument("invalid_argument 'dblWidth' in XArray2DFFTRe<T>::FilterSymEpanechnikov (filter width is non-positive)");
		
#ifdef HAVE_UNIQUE_PTR
		std::unique_ptr<IXAHead> pHeadAuto;
#else
		std::unique_ptr<IXAHead> pHeadAuto;
#endif

		IXAHead* pHead( m_rXArray2D.GetHeadPtr() );

		if( pHead )
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif
		XArray2D<T> Kernel( GetBaseObject() );
		Kernel.Fill( 1. );
		XArray2DMove<T> Move( Kernel );
		Move.SymEpanechnikovEnvelope( dblWidth );
		Convol( Kernel, true );
		
		if( pHead ) 
			m_rXArray2D.SetHeadPtr( pHeadAuto->Clone() ); // restore the head
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FilterGauss
	//
	// Filters the 'wrapped' real XArray2D object by convoluting it with a Gaussian kernel
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by convoluting it with a Gaussian kernel
		\param		dblSigmaY	Half-width in pixels of the kernel along the first dimension
		\param		dblSigmaX	Half-width in pixels of the kernel along the second dimension
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements the conventional Gaussian filtering in direct space. The filtering
			is performed using the FFT method.
	*/	
	template <class T> void XArray2DFFTRe<T>::FilterGauss(double dblSigmaY, double dblSigmaX)
	{
		if (dblSigmaX < 0 || dblSigmaY < 0)
				throw std::invalid_argument("invalid_argument 'dblSigmaX or dblSigmaY' in XArray2DFFTRe<T>::FilterGauss (standard deviation must be non-negative)"); 
		if (dblSigmaX > 0 || dblSigmaY > 0)
		{
#ifdef HAVE_UNIQUE_PTR
			std::unique_ptr<IXAHead> pHeadAuto;
#else
			std::unique_ptr<IXAHead> pHeadAuto;
#endif

			IXAHead* pHead = GetBaseObject().GetHeadPtr();
			
			if (pHead) 
#ifdef HAVE_UNIQUE_PTR
				pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
				pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif
			else
			{
				IXAHWave2D* pHeadTemp = CreateWavehead2D();
				pHeadTemp->SetData(1, 0, (double)GetBaseObject().GetDim1(), 0, (double)GetBaseObject().GetDim2());
				GetBaseObject().SetHeadPtr(pHeadTemp);
			}

			double yst = GetYStep(GetBaseObject());
			double xst = GetXStep(GetBaseObject());
			XArray2D<std::complex<T> > C;
			FFTRe(C, true);
			XArray2DMove<std::complex<T> > temp(C);
			//temp.GaussEnvelope(1.0 / (tPI * dblSigmaY * yst), 1.0 / (tPI * dblSigmaX * xst));
			temp.GaussEnvelopeFourier(dblSigmaY * yst, dblSigmaX * xst);
			FFTRe(C, false);
			
			if (pHead) 
				GetBaseObject().SetHeadPtr(pHeadAuto->Clone()); // restore the head
			else 
				GetBaseObject().SetHeadPtr(0); // restore the absent head
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FilterLorentz
	//
	// Filters the 'wrapped' real XArray2D object by convoluting it with a circularly symmetric Lorentz kernel
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by convoluting it with a Lorentz kernel
		\param		dblSigma	Half-diameter in pixels of the kernel
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements the conventional Lorentz filtering in direct space. The filtering
			is performed using the FFT method. A circularly symmetric kernel is used: c/(2Pi sigma^2) * (1 + c r^2/sigma^2)^(-3/2), c = 2^(2/3) - 1.
	*/	
	template <class T> void XArray2DFFTRe<T>::FilterLorentz(double dblSigma)
	{
		if (dblSigma < 0)
			throw std::invalid_argument("invalid_argument 'dblSigma' in XArray2DFFTRe<T>::FilterLorentz (standard deviation must be non-negative)"); 
		if (dblSigma == 0)
			return; // do nothing

#ifdef HAVE_UNIQUE_PTR
		std::unique_ptr<IXAHead> pHeadAuto;
#else
		std::unique_ptr<IXAHead> pHeadAuto;
#endif

		IXAHead* pHead = GetBaseObject().GetHeadPtr();
			
		if (pHead) 
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif
		else
		{
			IXAHWave2D* pHeadTemp = CreateWavehead2D();
			pHeadTemp->SetData(1, 0, (double)GetBaseObject().GetDim1(), 0, (double)GetBaseObject().GetDim2());
			GetBaseObject().SetHeadPtr(pHeadTemp);
		}

		double yst = GetYStep(GetBaseObject());
		double xst = GetXStep(GetBaseObject());
		if (yst != xst)
			throw std::invalid_argument("invalid_argument 'this' in XArray2DFFTRe<T>::FilterLorentz (image must have square pixels)");

		//temporarily pad the wrapped array along both dimensions to the nearest larger-or-equal integer power of 2 (number of pixels)
		index_t nyA = GetBaseObject().GetDim1(), nxA = GetBaseObject().GetDim2();
		if (nyA < 2 || nxA < 2)
			throw std::invalid_argument("invalid_argument 'this' in XArray2DFFTRe<T>::FilterLorentz (image must have at least 2 pixels along each direction)");

		index_t ny = 2;
		while (ny < nyA) ny *= 2;
		index_t nx = 2;
		while (nx < nxA) nx *= 2;

		XArray2DMove<T> MoveA(GetBaseObject());
		MoveA.PadMirror((ny - nyA) / 2, (ny - nyA) - (ny - nyA) / 2, (nx - nxA) / 2, (nx - nxA) - (nx - nxA) / 2);

		XArray2D<std::complex<T> > C;
		FFTRe(C, true);
		XArray2DMove<std::complex<T> > temp(C);
		temp.ExponentialEnvelopeFourier(dblSigma * xst);
		FFTRe(C, false);

		//trim back to the original size (as before padding)
		MoveA.Trim((ny - nyA) / 2, (ny - nyA) - (ny - nyA) / 2, (nx - nxA) / 2, (nx - nxA) - (nx - nxA) / 2);
			
		if (pHead) 
			GetBaseObject().SetHeadPtr(pHeadAuto->Clone()); // restore the head
		else 
			GetBaseObject().SetHeadPtr(0); // restore the absent head
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FourFilterRect
	//
	// Filters the 'wrapped' real XArray2D object by applying rectangular window to its Fourier spectrum
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by applying rectangular window to its Fourier spectrum
		\param		iYFreqCutoff	Cut-off frequency corresponding to the first dimension
		\param		iXFreqCutoff	Cut-off frequency corresponding to the second dimension
		\exception  std::invalid_argument is thrown if iYFreqCutoff or iXFreqCutoff exceed the spectral width of the 'wrapped' object
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements low-pass Fourier filtering by preserving only the frequencies below the set cut-offs
	*/	
	template <class T> void XArray2DFFTRe<T>::FourFilterRect(index_t iYFreqCutoff, index_t iXFreqCutoff)
	{
		if (iYFreqCutoff > m_rXArray2D.GetDim1() / 2 || iXFreqCutoff > m_rXArray2D.GetDim2() / 2)
			throw std::invalid_argument("invalid_argument 'iYFreqCutoff or iXFreqCutoff' in XArray2DFFTRe<T>::FourFilterRect (cut-off too large)");
		
#ifdef HAVE_UNIQUE_PTR
			std::unique_ptr<IXAHead> pHeadAuto;
#else
			std::unique_ptr<IXAHead> pHeadAuto;
#endif

		IXAHead* pHead = GetBaseObject().GetHeadPtr();

		if (pHead) 
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif
		
		XArray2D<std::complex<T> > C;
		FFTRe(C, true);
		index_t ymask = C.GetDim1() / 2 - iYFreqCutoff;
		index_t xmask = C.GetDim2() / 2 - iXFreqCutoff;
		XArray2DMove<std::complex<T> > temp(C);
		temp.Mask(ymask, ymask, xmask, xmask, 0);
		FFTRe(C, false);
		
		if (pHead) 
			GetBaseObject().SetHeadPtr(pHeadAuto->Clone()); // restore the head
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FourFilterGauss
	//
	// Filters the 'wrapped' real XArray2D object by applying Gaussian window to its Fourier spectrum
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by applying Gaussian window to its Fourier spectrum
		\param		dblFreqSigmaY	Spectral half-width of the window corresponding to the first dimension
		\param		dblFreqSigmaX	Spectral half-width of the window corresponding to the second dimension
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements low-pass Fourier filtering by multiplying the frequencies by a Gaussian envelope
	*/	
	template <class T> void XArray2DFFTRe<T>::FourFilterGauss(double dblFreqSigmaY, double dblFreqSigmaX)
	{
#ifdef HAVE_UNIQUE_PTR
			std::unique_ptr<IXAHead> pHeadAuto;
#else
			std::unique_ptr<IXAHead> pHeadAuto;
#endif
		
		IXAHead* pHead = GetBaseObject().GetHeadPtr();

		if (pHead) 
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif

		XArray2D<std::complex<T> > C;
		FFTRe(C, true);
		XArray2DMove<std::complex<T> > temp(C);
		temp.GaussEnvelope(dblFreqSigmaY * GetYStep(C), dblFreqSigmaX * GetXStep(C));
		FFTRe(C, false);
		
		if (pHead) 
			GetBaseObject().SetHeadPtr(pHeadAuto->Clone()); // restore the head
	}


	//---------------------------------------------------------------------------
	//Function XArray2DFFTRe<T>::FourFilterLorentz
	//
	// Filters the 'wrapped' real XArray2D object by applying Lorentzian window to its Fourier spectrum
	//
	/*!
		\brief		Filters the 'wrapped' real XArray2D object by applying Lorentzian window to its Fourier spectrum
		\param		dblGamma	G-parameter of the Lorentzian distribution G / [(k-k0)^2 + G^2]
		\param		iKy0	ky-centre (in pixels) of the Lorentzian distribution in Fourier space
		\param		iKx0	kx-centre (in pixels) of the Lorentzian distribution in Fourier space
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function implements low-pass Fourier filtering by multiplying the frequencies by a Lorentzian envelope.
			Normalization is applied to preserve the L1-norm.
	*/	
	template <class T> void XArray2DFFTRe<T>::FourFilterLorentz(double dblGamma, index_t iKy0, index_t iKx0)
	{
#ifdef HAVE_UNIQUE_PTR
		std::unique_ptr<IXAHead> pHeadAuto;
#else
		std::unique_ptr<IXAHead> pHeadAuto;
#endif
		IXAHead* pHead = GetBaseObject().GetHeadPtr();
		
		if (pHead) 
#ifdef HAVE_UNIQUE_PTR
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#else
			pHeadAuto = std::unique_ptr<IXAHead>( pHead->Clone() ); // preserve the head if it exists
#endif
		
		double dblNorm1 = m_rXArray2D.Norm(eNormL1);
		XArray2D<std::complex<T> > C;
		FFTRe(C, true);
	
		// Lorentzian envelope G / [(k-k0)^2 + G^2]
		double fyst2 = xar::GetYStep(C) * xar::GetYStep(C);
		double fxst2 = xar::GetXStep(C) * xar::GetXStep(C);
		double dblCentreY = iKy0 + 0.5 * C.GetDim1();
		double dblCentreX = iKx0 + 0.5 * C.GetDim2();
		double yi, yi2, xj, xj2; 
		double dblGamma2 = dblGamma * dblGamma;
		
		for (index_t i = 0; i < C.GetDim1(); i++)
		{
			yi = i - dblCentreY;
			yi2 = yi * yi * fyst2 + dblGamma2;
			
			for (index_t j = 0; j < C.GetDim2(); j++)
			{
				xj = j - dblCentreX;
				xj2 = xj * xj * fxst2;
				C[i][j] *= T(dblGamma / (yi2 + xj2));
			}
		}

		FFTRe(C, false);
		m_rXArray2D *= dblNorm1 / m_rXArray2D.Norm(eNormL1);
		
		if (pHead) 
			GetBaseObject().SetHeadPtr(pHeadAuto->Clone()); // restore the head
	}




	//---------------------------------------------------------------------------
	//Function XArray2DFFT<T>::FourierFilt
	//@@@@@@@@@@@@@ needs to be checked carefully
	//
	// Fourier filtering
	//
	//
	//
	/*!
		\brief		Fourier filtering
		\param		SpectFraction defines the fraction of the Fourier spectrum to be used (the rest ofthe spectrum is set to zero)
		\exception	std::invalid_argument is thrown if any of the two dimensions of the wrapped object
					is not an integer power of 2; or if the object does not have an associated Wavehead2D.
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
	*/	
	//	WARNING(old XY order!!!): internally, this function relates dim1 to X, and dim2 to Y, however, 
	//	X and Y	are swapped at the entry point of this function
	//
	template <class T> void XArray2DFFTRe<T>::FilterFourier(double SpectFraction)
	{
		if (SpectFraction<=0 || SpectFraction>=1)
			throw std::invalid_argument("invalid_argument 'SpectFraction' in XArray2DFFT<T>::FourierFilt (must be >0 and <1)");

		if (m_rXArray2D.GetValuetype()!=xar::eXAFloat && m_rXArray2D.GetValuetype()!=xar::eXADouble)
			throw std::invalid_argument("invalid_argument 'this' (must be float or double 2D array)");

		long nx = (long)m_rXArray2D.GetDim1(), ny = (long)m_rXArray2D.GetDim2();
		long nx0 = nx, ny0 = ny;
		double energ = m_rXArray2D.Norm(xar::eNormL2);

		//pad to the nearest larger power of 2
		long i = 2;
		while (i<nx) i *= 2;
		nx = i;
		long j = 2;
		while (j<ny) j *= 2;
		ny = j;
		XArray2DMove<T> tmp(m_rXArray2D);
		tmp.Pad((nx-nx0)/2, nx-nx0-(nx-nx0)/2, (ny-ny0)/2, ny-ny0-(ny-ny0)/2, (T)m_rXArray2D.Norm(xar::eNormAver));

		// store x,y limits
		const IXAHWave2D *pW2Head = (dynamic_cast<const IXAHWave2D*>(m_rXArray2D.GetHeadPtr()));
		if (!pW2Head)
			throw std::invalid_argument("invalid_argument 'this' (Wave2Head is missing)");
		double Ylo=pW2Head->GetYlo();
		double Yhi=pW2Head->GetYhi();
		double Xlo=pW2Head->GetXlo();
		double Xhi=pW2Head->GetXhi();

		//do FFT
		XArray2D< std::complex<T> > C;
		MakeComplex(m_rXArray2D, T(0), C, false);
		XArray2DFFT<T> CFFTService(C);
		CFFTService.FFT(true);

		//filter
		long nxmask = long(double(nx) / 2.0 * (1.0 - SpectFraction) + 0.5);
		long nymask = long(double(ny) / 2.0 * (1.0 - SpectFraction) + 0.5);
		XArray2DMove<std::complex<T> > temp(C);
		//temp.Mask(nxmask-1, nx-nxmask-1, nymask-1, ny-nymask-1, 0);
		temp.Mask(nxmask, nxmask, nymask, nymask, 0);

		//do inverse FFT
		CFFTService.FFT(false);
		Re(C, m_rXArray2D);

		//trim back and restore energy
		XArray2DMove<T> tmp1(m_rXArray2D);
		tmp1.Trim((nx-nx0)/2, nx-nx0-(nx-nx0)/2, (ny-ny0)/2, ny-ny0-(ny-ny0)/2);

		// restore X,Y limits
		IXAHWave2D *newHead = CreateWavehead2D();
		pW2Head = (dynamic_cast<const IXAHWave2D*>(m_rXArray2D.GetHeadPtr()));
		*newHead = *pW2Head;
		newHead->SetData(pW2Head->GetWl(), Ylo, Yhi, Xlo, Xhi);
		m_rXArray2D.SetHeadPtr(newHead);

		double energ1 = m_rXArray2D.Norm(xar::eNormL2);
		if (energ1) m_rXArray2D *= (T)(energ / energ1);
	}

} // namespace xar closed
//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//

// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray2DFFTRe<float>;
	template class xar::XArray2DFFTRe<double>;
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
