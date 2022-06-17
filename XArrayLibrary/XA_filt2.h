//Header XA_filt2.h
//
//
//	HEADER FILE TITLE:
//
//		Two-dimensional direct-space filters, noise generators, etc
//
/*!
	\file		XA_filt2.h
	\brief		Two-dimensional direct-space filters, noise generators, etc
	\par		Description:
		This header contains a class that provides 2D filters, noise generators
		and related services for XArray2D<T> objects
*/
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
#include "XA_move2.h"
//#include "XA_nropt.h"
#ifdef HAVE_STD_UNIFORM_DIST
#include <random>
#endif
#ifdef HAVE_BOOST_UNIFORM_DIST
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#endif

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
	#define NEG 1
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//
	inline double NonZero(double x) { return (fabs(x) < 1e-5 ? 0 : 1); }
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//Class XArray2DFilt<T>
//
//	Two-dimensional noise generation and filtering service class
//
/*!
	\brief		Two-dimensional noise generation and filtering service class
	\par		Description:
				This class template provides 2D noise generators, filters and
				related services for XArray2D<T> objects
	\remarks	An object of this class represents a thin 'wrapper' around an XArray2D
				object.	The wrapper	is an interface exposing several functions that provide 
				various 2D noise generation and filtering services				
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XArray2DFilt
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray2DFilt(XArray2D<T>& rXArray2D) : m_rXArray2D(rXArray2D) { GetValuetype(); }
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XArray2DFilt(const XArray2DFilt<T>& rCopy) : m_rXArray2D(rCopy.m_rXArray2D), m_vTemp(rCopy.m_vTemp) {}
	public:
		//! Destructor
		~XArray2DFilt() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XArray2DFilt<T>& rCopy);

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XArray2DFilt<T> objects for unsupported types T
		//! Returns the xar::_eValueType corresponding to T
		static _eValueType GetValuetype(void);
		//! Returns a reference to the non-modifiable 'wrapped' XArray2D<T> object
		const XArray2D<T>& GetBaseObject() const { return m_rXArray2D; }
		//! Returns a reference to the 'wrapped' XArray2D<T> object
		XArray2D<T>& GetBaseObject() { return m_rXArray2D; }


	// Operations
	public:
		//! Replaces every value in the 'wrapped' object by the average of its iYPoins by iXPoints neighbours
		void AverageFilt(index_t iYPoints, index_t iXPoints);
		//! Replaces every value in the 'wrapped' object by the value of Laplacian at that point
		void MLaplacian(bool DivideByA = false);
		//! Replaces every value in the 'wrapped' object by the value of -div(rXArIntensity*grad) at that point
		void MDivIGrad(const XArray2D<T>& rXArIntensity);
		//! Adds Gaussian noise with given average relative standard deviation to the 'wrapped' object 
		void GaussNoise(double dblRelAverStndDeviation);
		//! Adds Poisson (photon statistics) noise with given average relative standard deviation to the 'wrapped' object 
		void PoissNoise(double dblRelAverStndDeviation);
		//! Denoises the 'wrapped' object using the SUSAN algorithm
		void SUSAN( double dblSigma, double dblThreshold );
		//! Denoises the 'wrapped' object using the SUSAN algorithm with optimization
		void SUSANopt( double& dblSigma, double& dblThreshold, double dblRelStDevAver, bool bUsePoissonStat, double dblFtol, index_t& riNumIter );
		//! Applies median filter with a rectangular mask
		void MedianFilt(index_t iMaskSizeY, index_t iMaskSizeX);
		//! Applies thresholded median filter with a rectangular mask
		void MedianFiltThresh(index_t iMaskSizeY, index_t iMaskSizeX, T LowThresh, T HighThresh);
		//! Applies thresholded 1D median filter
		void MedianFiltThresh1D(index_t iMaskSize, T LowThresh, T HighThresh, bool bRowWise = true);
		//! Applies masked median filter with a rectangular mask
		void MedianFiltMask(index_t iMaskSizeY, index_t iMaskSizeX, XArray2D<char>& rxaMask);
		//! Applies masked 1D median filter
		void MedianFiltMask1D(index_t iMaskSize, XArray2D<char>& rxaMask, bool bRowWise);
		//! Applies a ring filter of a specified size (must be odd)
		void RingFilt(index_t nFiltSize, bool fTransposedSino);
		//! Applies a zingers filter of a specified size (must be odd, by default 9) and threshold (by default 1.2)
		void ZingersFilt(index_t nFiltSize = 9, double dblThreshold = 1.2, bool bTransposed = false);
		//!	Applies a filter that removes vertical streaks
		void VertStreaksFilt(double dblThresh, index_t iWidth);
		//! Applies a TV-based denoising filter
		void TVFilt(
			double dblSTD,
			double dblBetaRel,
			double dblSupportRadius,
			T dblMuStep,
			size_t numMu, 
			T dblMuMin,
			size_t numIter,
			std::vector<double>& rvChi2,
			std::vector<double>& rvRegul);

#ifdef USE_SUSAN_OPT
		// This function is called by optimization routines (does not belong to user interface)
		double CalcOnce( double* x, index_t n_dim );
#endif

	// Implementation
	protected:
		// Median 3x3 filter (used in the SUSAN algorithm)
		T Median3x3(XArray2D<T>& rXArIn, index_t i, index_t j);
		// Optimisation cycle used in TVFilt
		void OptimisationCPU(
			xar::XArray2D<T>& rxaSlice,
			const xar::XArray2D<T>& rxaNoisySlice,
			double dblSTD,
			const std::vector<size_t>& rvIDX,
			const std::vector<T>& rvMu,
			double dblBeta,
			T dblMuMin);
		// Regularisation term
		double Regul(
			xar::XArray2D<T>& rxaSlice,
			unsigned int iNorm);
		// Goodness-of-fit term
		double Chi2a(
			const xar::XArray2D<T>& rxaArray1,
			const xar::XArray2D<T>& rxaArray2,
			double dblSTD);

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
//! Returns the _eValueType corresponding to T=char
template<> inline _eValueType XArray2DFilt<char>::GetValuetype() { return eXAChar; }
//! Returns the _eValueType corresponding to T=short
template<> inline _eValueType XArray2DFilt<short>::GetValuetype() { return eXAShort; }
//! Returns the _eValueType corresponding to T=long
template<> inline _eValueType XArray2DFilt<long>::GetValuetype() { return eXALong; }
//! Returns the _eValueType corresponding to T=float
template<> inline _eValueType XArray2DFilt<float>::GetValuetype() { return eXAFloat; }
//! Returns the _eValueType corresponding to T=double
template<> inline _eValueType XArray2DFilt<double>::GetValuetype() { return eXADouble; }
//! Returns the _eValueType corresponding to T=fcomplex
template<> inline _eValueType XArray2DFilt<xar::fcomplex>::GetValuetype() { return eXAFComplex; }
//! Returns the _eValueType corresponding to T=dcomplex
template<> inline _eValueType XArray2DFilt<xar::dcomplex>::GetValuetype() { return eXADComplex; }


//! Assignment (declared protected to prohibit copying)
template <class T> void XArray2DFilt<T>::operator=(const XArray2DFilt<T>& xaf2)
{ 
	if (this == &xaf2) 
		return; 
	else
	{
		m_rXArray2D = xaf2.m_rXArray2D;
		m_vTemp = xaf2.m_vTemp;
	}
}


//! Replaces every value in the 'wrapped' object by the average of its iYPoins by iXPoints neighbours
// Averaging for 2D arrays:
// m_rXArray2D - 2D array to be filtered
// iXPoints, iYPoints - total number of points involved in each averaging
//
// NOTE: this filter averages up to the very edges
//		using the reflection of the array over the boundary
// NOTE:!this filter uses a numerically very efficient algorithm of 
//      successive 1D averagings, first in columns, then in rows,
//		and with only one 1D auxilliary array
template <class T> void XArray2DFilt<T>::AverageFilt(index_t iYPoints, index_t iXPoints)
{
	index_t i, j, nx, ny, nxa2, nya2, nxa22, nya22;
	double an;
	T *temp;
	
	std::swap(iYPoints, iXPoints);
	nx = m_rXArray2D.GetDim1() - 1;
	ny = m_rXArray2D.GetDim2() - 1;
	(iXPoints==0) ? nxa2 = 0 : nxa2 = (iXPoints - 1) / 2;
	(iYPoints==0) ? nya2 = 0 : nya2 = (iYPoints - 1) / 2;
	
	if (nxa2==0 && nya2==0) 
		return;
	
	if (nxa2<0 || nx<nxa2*2 || nya2<0 || ny<nya2*2 || !nya2 && !nxa2) 
		throw std::invalid_argument("invalid_argument 'iXPoints or iYPoints' in AverageFilt (unsuitable parameters for averaging)");
	
	nxa22 = nxa2 * 2;
	nya22 = nya2 * 2;
	an = 1.0 / (nxa2 * 2 + 1) / (nya2 * 2 + 1);

	try 
	{ 
		temp = new T [xar::max(2*nx+nxa22+2, 2*ny+nya22+2)]; 
	}
	catch(...) 
	{ 
		throw std::runtime_error("cannot allocate memory in AverageFilt"); 
	}
	if (temp==0) 
		throw std::runtime_error("cannot allocate memory in AverageFilt");
	
	for (j=0; j<=ny; j++) // averaging in columns
	{
		for (i=0; i<=nx; i++) 
			temp[i+nxa2] = m_rXArray2D[i][j];

		for (i=0; i<=nxa2-1; i++) 
			temp[i] = temp[nxa22-1-i]; // top edge reflection

		for (i=nx+nxa2+1; i<=nx+nxa22; i++) 
			temp[i] = temp[2*nx+nxa22+1-i]; // bottom edge reflection

		m_rXArray2D[0][j] = T(0.0);
		
		for (i=0; i<=nxa22; i++) 
			m_rXArray2D[0][j] += temp[i]; // initial element

		for (i=1; i<=nx; i++) 
			m_rXArray2D[i][j] = m_rXArray2D[i-1][j] + temp[i+nxa22] - temp[i-1];
	}

	for (i=0; i<=nx; i++) // averaging in rows
	{
		for (j=0; j<=ny; j++) 
			temp[j+nya2] = m_rXArray2D[i][j];

		for (j=0; j<=nya2-1; j++) 
			temp[j] = temp[nya22-1-j]; //left edge reflection

		for (j=ny+nya2+1; j<=ny+nya22; j++) 
			temp[j]=temp[2*ny+nya22+1-j]; //right edge reflection

		m_rXArray2D[i][0] = T(0.0);

		for (j=0; j<=nya22; j++) 
			m_rXArray2D[i][0] += temp[j]; //initial element

		for (j=1; j<=ny; j++) 
			m_rXArray2D[i][j] = m_rXArray2D[i][j-1] + temp[j+nya22] - temp[j-1];
	}
		
	for (i=0; i<=nx; i++) // normalising
		for (j=0; j<=ny; j++) 
			m_rXArray2D[i][j] *= T(an);

	delete [] temp;
}

// Specialization for T=fcomplex
// Averaging for 2D arrays:
// m_rXArray2D - 2D array to be filtered
// iXPoints, iYPoints - total number of points involved in each averaging
//
// NOTE: this filter averages up to the very edges
//		using the reflection of the array over the boundary
// NOTE:!this filter uses a numerically very efficient algorithm of 
//      successive 1D averagings, first in columns, then in rows,
//		and with only one 1D auxilliary array
template<> inline void XArray2DFilt<fcomplex>::AverageFilt(index_t iYPoints, index_t iXPoints)
{
	index_t i, j, nx, ny, nxa2, nya2, nxa22, nya22;
	double an;
	fcomplex *temp;
	
	std::swap(iYPoints, iXPoints);
	nx = m_rXArray2D.GetDim1() - 1;
	ny = m_rXArray2D.GetDim2() - 1;
	(iXPoints==0) ? nxa2 = 0 : nxa2 = (iXPoints - 1) / 2;
	(iYPoints==0) ? nya2 = 0 : nya2 = (iYPoints - 1) / 2;
	
	if (nxa2==0 && nya2==0) 
		return;
	
	if (nxa2<0 || nx<nxa2*2 || nya2<0 || ny<nya2*2 || !nya2 && !nxa2) 
		throw std::invalid_argument("invalid_argument 'iXPoints or iYPoints' in AverageFilt (unsuitable parameters for averaging)");
	
	nxa22 = nxa2 * 2;
	nya22 = nya2 * 2;
	an = 1.0 / (nxa2 * 2 + 1) / (nya2 * 2 + 1);

	try 
	{ 
		temp = new fcomplex [xar::max(2*nx+nxa22+2, 2*ny+nya22+2)]; 
	}
	catch(...) 
	{ 
		throw std::runtime_error("cannot allocate memory in AverageFilt"); 
	}

	if (temp==0) 
		throw std::runtime_error("cannot allocate memory in AverageFilt");
	
	for (j=0; j<=ny; j++) // averaging in columns
	{
		for (i=0; i<=nx; i++) 
			temp[i+nxa2] = m_rXArray2D[i][j];

		for (i=0; i<=nxa2-1; i++) 
			temp[i] = temp[nxa22-1-i]; // top edge reflection

		for (i=nx+nxa2+1; i<=nx+nxa22; i++) 
			temp[i] = temp[2*nx+nxa22+1-i]; // bottom edge reflection

		m_rXArray2D[0][j] = 0.0f;

		for (i=0; i<=nxa22; i++) 
			m_rXArray2D[0][j] += temp[i]; // initial element
		
		for (i=1; i<=nx; i++) 
			m_rXArray2D[i][j] = m_rXArray2D[i-1][j] + temp[i+nxa22] - temp[i-1];
	}

	for (i=0; i<=nx; i++) // averaging in rows
	{
		for (j=0; j<=ny; j++) 
			temp[j+nya2] = m_rXArray2D[i][j];
		
		for (j=0; j<=nya2-1; j++) 
			temp[j] = temp[nya22-1-j]; //left edge reflection
		
		for (j=ny+nya2+1; j<=ny+nya22; j++) 
			temp[j]=temp[2*ny+nya22+1-j]; //right edge reflection

		m_rXArray2D[i][0] = 0.0f;

		for (j=0; j<=nya22; j++) 
			m_rXArray2D[i][0] += temp[j]; //initial element

		for (j=1; j<=ny; j++) 
			m_rXArray2D[i][j] = m_rXArray2D[i][j-1] + temp[j+nya22] - temp[j-1];
	}
		
	for (i=0; i<=nx; i++) // normalising
		for (j=0; j<=ny; j++) 
			m_rXArray2D[i][j] *= float(an);

	delete [] temp;
}


//! Replaces every value in the 'wrapped' object by the value of Laplacian at that point
// Replaces an Array2D object m_rXArray2D with -Laplacian(m_rXArray2D), continuing as const at the edges
// if DivideByA==true, also divides the result by m_rXArray2D
template <class T> void XArray2DFilt<T>::MLaplacian(bool DivideByA)
{
	index_t i, j;
	T a, b, c, d, e;
	
	index_t nx = m_rXArray2D.GetDim1(), nx1 = nx - 1;
	index_t ny = m_rXArray2D.GetDim2(), ny1 = ny - 1;

	double xst = GetYStep(m_rXArray2D);
	double yst = GetXStep(m_rXArray2D);

	T axst2 = T(-1.0 / (xst * xst));
	T ayst2 = T(-1.0 / (yst * yst));

	XArray2D<T> Temp(m_rXArray2D);

	for (i=0; i<nx; i++)
		for (j=0; j<ny; j++)
		{
			if (i==0) 
			{ 
				a = Temp[i+1][j]; 
				b = Temp[i][j]; 
			}
			else if (i==nx1) 
			{ 
				a = Temp[i][j]; 
				b = Temp[i-1][j]; 
			}
			else 
			{ 
				a = Temp[i+1][j]; 
				b = Temp[i-1][j]; 
			}

			if (j==0) 
			{ 
				c = Temp[i][j+1]; 
				d = Temp[i][j]; 
			}
			else if (j==ny1) 
			{ 
				c = Temp[i][j]; 
				d = Temp[i][j-1]; 
			}
			else 
			{ 
				c = Temp[i][j+1]; 
				d = Temp[i][j-1]; 
			}

			e = -T(2.0) * Temp[i][j];

			m_rXArray2D[i][j] = (a + b + e) * axst2 + (c + d + e) * ayst2;

			if (DivideByA)
				if (Temp[i][j]==T(0)) 
					throw std::invalid_argument("invalid_argument 'm_rXArray2D' in MLaplacian (zero value)");
				else 
					m_rXArray2D[i][j] /= Temp[i][j];
		}
}


// Specialization for T=fcomplex
// Replaces an Array2D object m_rXArray2D with -Laplacian(m_rXArray2D), continuing as const at the edges
// if DivideByA==true, also divides the result by m_rXArray2D
template<> inline void XArray2DFilt<fcomplex>::MLaplacian(bool DivideByA)
{
	index_t i, j;
	fcomplex a, b, c, d, e;
	
	index_t nx = m_rXArray2D.GetDim1(), nx1 = nx - 1;
	index_t ny = m_rXArray2D.GetDim2(), ny1 = ny - 1;

	double xst = GetYStep(m_rXArray2D);
	double yst = GetXStep(m_rXArray2D);

	float axst2 = float(-1.0 / (xst * xst));
	float ayst2 = float(-1.0 / (yst * yst));

	XArray2D<fcomplex> Temp(m_rXArray2D);

	for (i=0; i<nx; i++)
		for (j=0; j<ny; j++)
		{
			if (i==0) 
			{ 
				a = Temp[i+1][j]; 
				b = Temp[i][j]; 
			}
			else if (i==nx1) 
			{ 
				a = Temp[i][j]; 
				b = Temp[i-1][j]; 
			}
			else 
			{ 
				a = Temp[i+1][j]; 
				b = Temp[i-1][j]; 
			}

			if (j==0) 
			{ 
				c = Temp[i][j+1]; 
				d = Temp[i][j]; 
			}
			else if (j==ny1) 
			{ 
				c = Temp[i][j]; 
				d = Temp[i][j-1]; 
			}
			else 
			{ 
				c = Temp[i][j+1]; 
				d = Temp[i][j-1]; 
			}

			e = -2.0f * Temp[i][j];

			m_rXArray2D[i][j] = (a + b + e) * axst2 + (c + d + e) * ayst2;

			if (DivideByA)
				if (Temp[i][j]==0.0f) 
					throw std::invalid_argument("invalid_argument 'm_rXArray2D' in MLaplacian (zero value)");
				else 
					m_rXArray2D[i][j] /= Temp[i][j];
		}
}


//! Replaces every value in the 'wrapped' object by the value of -div(rXArIntensity*grad) at that point
// Replaces an Array2D object m_rXArray2D with -div(rXArIntensity(x,y)grad(m_rXArray2D)), continuing as const at the edges
template <class T> void XArray2DFilt<T>::MDivIGrad(const XArray2D<T>& rXArIntensity)
{
	index_t nx = m_rXArray2D.GetDim1(), nx1 = nx - 1;
	index_t ny = m_rXArray2D.GetDim2(), ny1 = ny - 1;
	double xst = GetYStep(m_rXArray2D);
	double yst = GetXStep(m_rXArray2D);
	
	if (nx != rXArIntensity.GetDim1() || ny != rXArIntensity.GetDim2())
		throw std::invalid_argument("invalid_argument 'rXArIntensity' in MDivIGrad (size differs from that of m_rXArray2D)");

	index_t i, j;
	T a, b, c, d, aI, bI, cI, dI, e;
	
	T axst2 = T(-1.0 / (xst * xst));
	T ayst2 = T(-1.0 / (yst * yst));
	T axyst = T(-0.25 / (xst * yst));

	XArray2D<T> Temp(m_rXArray2D);

	for (i=0; i<nx; i++)
		for (j=0; j<ny; j++)
		{
			if (i==0) 			
			{ 
				a=Temp[i+1][j]; 
				b=Temp[i][j]; 
				aI=rXArIntensity[i+1][j]; 
				bI=rXArIntensity[i][j]; 
			}
			else if (i==nx1) 			
			{ 
				a=Temp[i][j]; 
				b=Temp[i-1][j]; 
				aI=rXArIntensity[i][j]; 
				bI=rXArIntensity[i-1][j]; 
			}
			else 
			{ 
				a=Temp[i+1][j]; 
				b=Temp[i-1][j]; 
				aI=rXArIntensity[i+1][j]; 
				bI=rXArIntensity[i-1][j]; 
			}

			if (j==0) 
			{ 
				c=Temp[i][j+1]; 
				d=Temp[i][j]; 
				cI=rXArIntensity[i][j+1]; 
				dI=rXArIntensity[i][j]; 
			}
			else if (j==ny1) 
			{ 
				c=Temp[i][j]; 
				d=Temp[i][j-1]; 
				cI=rXArIntensity[i][j]; 
				dI=rXArIntensity[i][j-1]; 
			}
			else 
			{ 
				c=Temp[i][j+1]; 
				d=Temp[i][j-1]; 
				cI=rXArIntensity[i][j+1]; 
				dI=rXArIntensity[i][j-1]; 
			}

			e = -T(2.0) * Temp[i][j];

			m_rXArray2D[i][j] = axyst * ((aI-bI) * (a-b) + (cI-dI) * (c-d)) +
				rXArIntensity[i][j] * ((a + b + e) * axst2 + (c + d + e) * ayst2);				
		}
}


// Specialization for T=fcomplex
// Replaces an Array2D object m_rXArray2D with -div(rXArIntensity(x,y)grad(m_rXArray2D)), continuing as const at the edges
template<> inline void XArray2DFilt<fcomplex>::MDivIGrad(const XArray2D<fcomplex>& rXArIntensity)
{
	index_t nx = m_rXArray2D.GetDim1(), nx1 = nx - 1;
	index_t ny = m_rXArray2D.GetDim2(), ny1 = ny - 1;
	double xst = GetYStep(m_rXArray2D);
	double yst = GetXStep(m_rXArray2D);
	
	if (nx != rXArIntensity.GetDim1() || ny != rXArIntensity.GetDim2())
		throw std::invalid_argument("invalid_argument 'rXArIntensity' in MDivIGrad (size differs from that of m_rXArray2D)");

	index_t i, j;
	fcomplex a, b, c, d, aI, bI, cI, dI, e;
	
	float axst2 = float(-1.0 / (xst * xst));
	float ayst2 = float(-1.0 / (yst * yst));
	float axyst = float(-0.25 / (xst * yst));

	XArray2D<fcomplex> Temp(m_rXArray2D);

	for (i=0; i<nx; i++)
		for (j=0; j<ny; j++)
		{
			if (i==0) 
			{ 
				a=Temp[i+1][j]; 
				b=Temp[i][j]; 
				aI=rXArIntensity[i+1][j]; 
				bI=rXArIntensity[i][j]; 
			}
			else if (i==nx1) 
			{ 
				a=Temp[i][j]; 
				b=Temp[i-1][j]; 
				aI=rXArIntensity[i][j]; 
				bI=rXArIntensity[i-1][j]; 
			}
			else 
			{ 
				a=Temp[i+1][j]; 
				b=Temp[i-1][j]; 
				aI=rXArIntensity[i+1][j]; 
				bI=rXArIntensity[i-1][j]; 
			}

			if (j==0) 
			{ 
				c=Temp[i][j+1]; 
				d=Temp[i][j]; 
				cI=rXArIntensity[i][j+1]; 
				dI=rXArIntensity[i][j]; 
			}
			else if (j==ny1) 
			{ 
				c=Temp[i][j]; 
				d=Temp[i][j-1]; 
				cI=rXArIntensity[i][j]; 
				dI=rXArIntensity[i][j-1]; 
			}
			else 
			{ 
				c=Temp[i][j+1]; 
				d=Temp[i][j-1]; 
				cI=rXArIntensity[i][j+1]; 
				dI=rXArIntensity[i][j-1]; 
			}

			e = -2.0f * Temp[i][j];

			m_rXArray2D[i][j] = axyst * ((aI-bI) * (a-b) + (cI-dI) * (c-d)) +
				rXArIntensity[i][j] * ((a + b + e) * axst2 + (c + d + e) * ayst2);				
		}
}

#if(0)

//! Adds Gaussian noise with given average relative standard deviation to the 'wrapped' object 
// Addition of Gaussian photon noise to  intensity u-data of Array2D:
// NOTE: array elements are assumed to be intensity (not amplitude) !
// Relative standard deviation is defined for average intensity as 
// RelAverStndDeviation = 1 / sqrt(Nav), where Nav is the number of photons 
// registered at the point of average intensity.
// Our intensity data contains an unknown multiplicative scaling factor c,
// so that I = N / c. Then Iav = Nav / c = 1 / (c * RelAverStndDeviation^2),
// hence c = 1 / (Iav * RelAverStndDeviation^2). Now we can find the relative standard
// deviation of noise at an arbitrary pixel, 
// RelStndDeviation = 1 / sqrt(N) = 1 / sqrt(c * I) = RelAverStndDeviation * sqrt(Iav / I).
// I_noisy = I_clean + RelAverStndDeviation * sqrt(Iav_clean * I_clean) * Ran, where
// Ran is a Gaussian pseudo-random distribution with zero mean and unit variance.
template <class T> void XArray2DFilt<T>::GaussNoise( double dblRelAverStndDeviation )
{
	if( dblRelAverStndDeviation <= 0 )
		throw std::invalid_argument("invalid_argument 'dblRelAverStndDeviation' in GaussNoise (non-positive average relative standard deviation)");

	if( m_rXArray2D.Norm( xar::eNormMin ) < 0 )
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in GaussNoise (negative values in input intensity)");

	// Setup Mersenne twister engine & initialise from a random device
	std::random_device rd;
	std::mt19937 gen( rd() );
	std::normal_distribution<> distNormal( 0, 1 );

	double factor( dblRelAverStndDeviation * sqrt( m_rXArray2D.Norm( xar::eNormAver ) ) );

	for( T* p = m_rXArray2D.begin(); p != m_rXArray2D.end(); p++ )
		*p += static_cast<T>( factor * sqrt( static_cast<double>( *p ) ) * distNormal( gen ) );
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::GaussNoise(double dblRelAverStndDeviation)
{
	throw std::invalid_argument("invalid_argument '*this' in GaussNoise (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::GaussNoise(double dblRelAverStndDeviation)
{
	throw std::invalid_argument("invalid_argument '*this' in GaussNoise (not defined for dcomplex data)");
}


//! Adds Poisson (photon statistics) noise with given average relative standard deviation to the 'wrapped' object 
// Poisson (photon counting statistics) noise for u-data of Array2D:
// NOTE: array elements are assumed to be intensity (not amplitude) !
// Relative standard deviation is defined at the point of average intensity as 
// RelAverStndDeviation = 1 / sqrt(Nav), where Nav is the number of photons 
// registered at the point of average intensity.
// Our intensity data contains an unknown multiplicative scaling factor c,
// so that I = N / c. Then Iav = Nav / c = 1 / (c * RelAverStndDeviation^2),
// hence c = 1 / (Iav * RelAverStndDeviation^2). Now we can find the standard
// deviation of noise at an arbitrary pixel, 
// st_dev = sqrt(N) = sqrt(c * I) = sqrt(I / Iav) / RelAverStndDeviation
// I_noisy = Poiss(I / Iav / RelAverStndDeviation^2) / c , where
// Poiss(x) is a Poisson pseudo-random distribution with the mean and variance equal to x.
template <class T> void XArray2DFilt<T>::PoissNoise( double dblRelAverStndDeviation )
{
	if( dblRelAverStndDeviation <= 0 )
		throw std::invalid_argument("invalid_argument 'dblRelAverStndDeviation' in PoissNoise (non-positive relative standard deviation)");

	if( m_rXArray2D.Norm( xar::eNormMin ) < 0 )
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in PoissNoise (negative values in input intensity)");

	double dblAver( m_rXArray2D.Norm( xar::eNormAver ) );

	if( dblAver == 0 ) 
		return; // flat image (as there are no negative values)

	double dblC( 1.0 / ( dblAver * dblRelAverStndDeviation * dblRelAverStndDeviation ) );
	double dblAc( 1.0 / dblC );

	std::random_device rd;
	std::mt19937 gen( rd() );

	for( T* p = m_rXArray2D.begin(); p != m_rXArray2D.end(); p++ )
	{
		// Define Poisson distribution
		if( *p == T( 0 ) ) 
			continue;

		std::poisson_distribution<unsigned long> disPoisson( static_cast<double>( *p ) * dblC );
		*p = static_cast<T>( static_cast<double>( disPoisson( gen ) ) * dblAc );
	}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::PoissNoise(double dblRelAverStndDeviation)
{
	throw std::invalid_argument("invalid_argument '*this' in PoissNoise (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::PoissNoise(double dblRelAverStndDeviation)
{
	throw std::invalid_argument("invalid_argument '*this' in PoissNoise (not defined for dcomplex data)");
}


//---------------------------------------------------------------------------
//Function XArray2DFilt<T>::Median3x3
//
//	Service function performing 3x3 median filtering
//
//
//	\brief		Calculates smoothed image in point (i,j) of the input array using 3x3 median filter
//	\param		rXArIn input array which is to be smoothed 
//	\param		i defines the y-coordinate of the smoothed point
//	]param		j defines the x-coordinate of the smoothed point
//	\return		\a None
//	\par		Description:
//		This function smooths the wrapped array in the point (i,j) using 3x3 median filter
//		if the USAN area (denominator in the formula (36) in the original paper of S.M.Smith
//		and J.M.Brady) is equal to zero 
//
template<class T> T XArray2DFilt<T>::Median3x3(XArray2D<T>& rXArIn, index_t i, index_t j)
{
	long k, l;
	T p[8];
	T tmp;

	p[0] = rXArIn[i - 1][j - 1];
	p[1] = rXArIn[i - 1][j];
	p[2] = rXArIn[i - 1][j + 1];
	p[3] = rXArIn[i][j - 1];
	p[4] = rXArIn[i][j + 1];
	p[5] = rXArIn[i + 1][j - 1];
	p[6] = rXArIn[i + 1][j];
	p[7] = rXArIn[i + 1][j + 1];

	for(k = 0; k < 7; k++)
		for(l = 0; l < (7 - k); l++)
			if (p[l] > p[l + 1])
			{
				tmp = p[l];
				p[l] = p[l + 1];
				p[l + 1] = tmp;
			}

	return ( (p[3] + p[4]) / 2 );
}

// Specialization for T=fcomplex
template<> inline fcomplex XArray2DFilt<fcomplex>::Median3x3(XArray2D<fcomplex>& rXArIn, index_t i, index_t j)
{
	throw std::invalid_argument("invalid_argument '*this' in Median3x3 (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline dcomplex XArray2DFilt<dcomplex>::Median3x3(XArray2D<dcomplex>& rXArIn, index_t i, index_t j)
{
	throw std::invalid_argument("invalid_argument '*this' in Median3x3 (not defined for dcomplex data)");
}

//---------------------------------------------------------------------------
//Function XArray2DFilt<T>::SUSAN
//
//	Denoises the 'wrapped' object using the SUSAN algorithm
//
/*!
	\brief		Denoises the 'wrapped' object using the SUSAN algorithm
	\param		dblSigma controls the scale of the spatial filtering
	\param		dblThreshold is the brightness threshold (controls the scale of the "brightness smoothing")
	\return		\a None
	\par		Description:
		This function filters the wrapped array using the SUSAN algorithm. After Smith, S.M. and Brady, J.M., 
		"SUSAN - A New Approach to Low Level Image Processing", Int. Journal of Computer Vision, 23, pp.45-78 (1997).
		SUSAN = Smallest Univalue Segment Assimilating Nucleus.
		There are two thresholds which can be set at run-time. These are the brightness threshold (dblThreshold)
		and the distance threshold (dblSigma). In SUSAN smoothing dblSigma controls the size of the Gaussian mask.
		Increasing dblSigma gives more smoothing. The size of the mask is equal to 3 * dblSigma + 3.
		dblThreshold determines the maximum difference in greylevels between two pixels which allows them to be 
		considered part of the same "region" in the image. Reducing dblThreshold gives less smoothing, and vice versa.
		With SUSAN smoothing, more smoothing can also be obtained by iterating the algorithm several times. 
		This has a different effect	from varying dblSigma or dblThreshold.
*/	
template<class T> void XArray2DFilt<T>::SUSAN(double dblSigma, double dblThreshold)
{
	if (dblThreshold <= 0) 
		throw std::invalid_argument("invalid argument 'dblThreshold' in XArray2DFilt<T>::SUSAN (must be positive)");

	if (dblSigma <= 0) 
		throw std::invalid_argument("invalid argument 'dblSigma' in XArray2DFilt<T>::SUSAN (must be positive)");

	long i, j, k, l;
	long x_size = long(m_rXArray2D.GetDim2());
	long y_size = long(m_rXArray2D.GetDim1());
	T tThreshold = T(dblThreshold);

	long mask_size = (long)(1.5 * dblSigma + 0.5) + 1;

	if ( (2 * mask_size + 1 > x_size) || (2 * mask_size + 1 > y_size) )
		throw std::invalid_argument("invalid argument 'dblSigma' in XArray2DFilt<T>::SUSAN (mask size exceeds original matrix size)");

	XArray2D<T> Copy(m_rXArray2D);
	XArray2DMove<T> CopyMove(Copy);
	CopyMove.PadMirror(mask_size, mask_size, mask_size, mask_size);
	Copy /= tThreshold; // Now Copy contains relative intensities

	index_t n_max = (mask_size * 2) + 1;
	XArray2D<T> XArSpatial(n_max, n_max, T(0));

	// Filling in Spatial coefficientd matrix
	double temp = -0.5 / (dblSigma * dblSigma);
	for (i = -mask_size; i <= mask_size; i++)
		for (j = -mask_size; j <= mask_size; j++)
			XArSpatial[i + mask_size][j + mask_size] = (T)exp( (i * i + j * j) * temp );

	// main section
	T area, total, centre, temp1, res, zzz;
	for (i = mask_size; i < y_size + mask_size; i++)
	{
		for (j = mask_size; j < x_size + mask_size; j++)
		{
			area = 0;
			total = 0;
			centre = Copy[i][j];
			for (k = -mask_size; k <= mask_size; k++)
			{
				for(l = -mask_size; l <= mask_size; l++)
				{
					temp1 = Copy[i + k][j + l];
					res = temp1 - centre;
					zzz = XArSpatial[k + mask_size][l + mask_size] * T(exp(double(-res * res)));
					area += zzz;
					total += zzz * temp1;
				}
			}
			area -= 1;
			if (area == 0)
				m_rXArray2D[i - mask_size][j - mask_size] = tThreshold * Median3x3(Copy, i, j);
			else
				m_rXArray2D[i - mask_size][j - mask_size] = tThreshold * (total - centre) / area;
		}
	}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::SUSAN(double dblSigma, double dblThreshold)
{
	throw std::invalid_argument("invalid_argument '*this' in SUSAN (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::SUSAN(double dblSigma, double dblThreshold)
{
	throw std::invalid_argument("invalid_argument '*this' in SUSAN (not defined for dcomplex data)");
}

#ifdef USE_SUSAN_OPT
//---------------------------------------------------------------------------
//Function XArray2DFilt<T>::SUSANopt
//
// Denoises the 'wrapped' object using the SUSAN algorithm with optimization
//
/*!
	\brief		Denoises the 'wrapped' object using the SUSAN algorithm with optimization
	\param		dblSigma controls the scale of the spatial filtering, optimal value returned
	\param		dblThreshold controls the scale of the "brightness smoothing", optimal value returned
	\param		dblRelStDevAver relative standard deviation of the noise corresponding to the average intensity in the image
	\param		bUsePoissonStat if true, Poisson noise statistics is assumed, else additive Gaussian noise is assumed
	\param		dblFtol approximation tolerance (e.g. 0.1)
	\param		riNumIter	Maximum number of iterations as input, and the actual number of iterations as output
	\return		\a None
	\par		Description:
		This function filters the wrapped array using the SUSAN algorithm. After Smith, S.M. and Brady, J.M., 
		"SUSAN - A New Approach to Low Level Image Processing", Int. Journal of Computer Vision, 23, pp.45-78 (1997).
		SUSAN = Smallest Univalue Segment Assimilating Nucleus.
		There are two thresholds which can be set at run-time. These are the brightness threshold (dblThreshold)
		and the distance threshold (dblSigma). In SUSAN smoothing dblSigma controls the size of the Gaussian mask.
		Increasing dblSigma gives more smoothing. The size of the mask is equal to 3 * dblSigma + 3.
		dblThreshold determines the maximum difference in greylevels between two pixels which allows them to be 
		considered part of the same "region" in the image. Reducing dblThreshold gives less smoothing, and vice versa.
		The optimal values of dblSigma and dblThreshold are found inside this program as a function of noise level, 
		dblRelStDevAver, noise statistics, bUsePoissonStat, and optimization tolerance, dblFtol.
		On exit m_rXArray2D contains denoised image calculated by SUSAN with optimal parameters.
*/	
template<class T> void XArray2DFilt<T>::SUSANopt(double& dblSigma, double& dblThreshold, double dblRelStDevAver, bool bUsePoissonStat, double dblFtol, index_t& riNumIter)
{
	if (dblThreshold <= 0) 
		throw std::invalid_argument("invalid argument 'dblThreshold' in XArray2DFilt<T>::SUSANopt (cannot be zero)");

	if (dblSigma < 0) 
		throw std::invalid_argument("invalid argument 'dblSigma' in XArray2DFilt<T>::SUSANopt (cannot be negative)");

	double** p = new double*[3];
	for (int i = 0; i < 3; i++) p[i] = new double[2];
	double y[3];
	double x[2];

	index_t nmax = 100;
	
	index_t nfunk = 0;

	m_vTemp.resize(2, 0);
	m_vTemp[0] = &dblRelStDevAver; // use m_pTemp[0] to pass the dblRelStDevAver argument to CalcOnce
	m_vTemp[1] = &bUsePoissonStat; // use m_vTemp[1] to pass the bUsePoissonStat argument to CalcOnce

	x[0] = dblSigma;
	x[1] = dblThreshold;
	p[0][0] = x[0];
	p[0][1] = x[1];
	y[0] = CalcOnce(x, 2);

	double f_min = y[0];
	double dblSigmaOpt = x[0];
	double dblThresholdOpt = x[1];

	x[0] = 1.4 * dblSigma;
	x[1] = dblThreshold;
	p[1][0] = x[0];
	p[1][1] = x[1];
	y[1] = CalcOnce(x, 2);

	if (y[1] < f_min)
	{
		f_min = y[1];
		dblSigmaOpt = x[0];
		dblThresholdOpt = x[1];
	}

	x[0] = dblSigma;
	x[1] = 1.4 * dblThreshold;
	p[2][0] = x[0];
	p[2][1] = x[1];
	y[2] = CalcOnce(x, 2);

	if (y[2] < f_min)
	{
		f_min = y[2];
		dblSigmaOpt = x[0];
		dblThresholdOpt = x[1];
	}

	if (f_min > dblFtol)
	{
		amoebaC<xar::XArray2DFilt<T> >(p, y, 2, dblFtol, *this, &nfunk, nmax);
		dblSigma = fabs(p[0][0]);
		dblThreshold = fabs(p[0][1]);
//		f_min = y[0];
	}
	else
	{
		dblSigma = dblSigmaOpt;
		dblThreshold = dblThresholdOpt;
	}

	SUSAN(dblSigma, dblThreshold);

	for (int i = 0; i < 3; i++) delete[] p[i];
	delete[] p;

	riNumIter = nfunk;
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::SUSANopt(double& dblSigma, double& dblThreshold, double dblRelStDevAver, bool bUsePoissonStat, double dblFtol, index_t& riNumIter)
{
	throw std::invalid_argument("invalid_argument '*this' in SUSANopt (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::SUSANopt(double& dblSigma, double& dblThreshold, double dblRelStDevAver, bool bUsePoissonStat, double dblFtol, index_t& riNumIter)
{
	throw std::invalid_argument("invalid_argument '*this' in SUSANopt (not defined for dcomplex data)");
}


// Although n_dim-dimensional array is assumed, in this case n_dim is equal to 2 as SUSAN
// has two parameters. Therefore, n_dim is an excess parameter which shoukd be used only to 
// provide compatibility with 'amoeba' procedure in XA_nropt.h.
// To restrict sigma and threshold in SUSAN by positive values, here the absolute values are
// substituted.
template<class T> double XArray2DFilt<T>::CalcOnce(double* x, index_t n_dim)
{
	// Initially m_rXArray2D contains a noisy data
	XArray2D<T> Data(m_rXArray2D);	// Data contains copy of the initial data
	// x[0] is a dblSigma and x[1] is a dblThreshold for SUSAN filter
	double temp = x[1];

	if (fabs(temp) < 1e-5) 
		x[1] = 1e-5 * (temp < 0.0 ? -1.0 : 1.0);

	SUSAN(fabs(x[0]), fabs(x[1]));	// Now m_rXArray2D contains smoothed data
	std::swap(m_rXArray2D, Data);	// Now m_rXArray2D contains initial data and Data contains
									// smoothed data
	// m_vTemp[0] is used to deliver dblRelStDevAver value
	double dblRelStDevAver = *(reinterpret_cast<double*>(m_vTemp[0]));
	// m_vTemp[1] is used to deliver bUsePoissonStat value
	bool bUsePoissonStat = *(reinterpret_cast<bool*>(m_vTemp[1]));
	return fabs(m_rXArray2D.Chi2(Data, dblRelStDevAver, bUsePoissonStat) / (m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2()) - 1.0);
}

// Specialization for T=fcomplex
template<> inline double XArray2DFilt<fcomplex>::CalcOnce(double* x, index_t n_dim)
{
	throw std::invalid_argument("invalid_argument '*this' in CalcOnce (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline double XArray2DFilt<dcomplex>::CalcOnce(double* x, index_t n_dim)
{
	throw std::invalid_argument("invalid_argument '*this' in CalcOnce (not defined for dcomplex data)");
}
#endif

//! Applies median filter with a rectangular mask
template<class T> void XArray2DFilt<T>::MedianFilt(index_t iMaskSizeY, index_t iMaskSizeX)
{
	if (iMaskSizeY > m_rXArray2D.GetDim1() || iMaskSizeX > m_rXArray2D.GetDim2() || iMaskSizeY < 1 || iMaskSizeX < 1)
		throw std::invalid_argument("invalid argument 'iMaskSizeY or iMaskSizeX' in XArray2DFilt<T>::MedianFilt");

	index_t iOffsetY = iMaskSizeY / 2;
	index_t iOffsetY1 = iMaskSizeY - iMaskSizeY / 2 - 1;
	index_t iOffsetX = iMaskSizeX / 2;
	index_t iOffsetX1 = iMaskSizeX - iMaskSizeX / 2 - 1;
	XArray2D<T> xaSrc(m_rXArray2D); // create a temp copy
	XArray2DMove<T> xaSrcMove(xaSrc);
	xaSrcMove.PadMirror(iOffsetY, iOffsetY1, iOffsetX, iOffsetX1);

	XArray2D<T> xaTemp(iMaskSizeY, iMaskSizeX);

	for (index_t i = 0; i < m_rXArray2D.GetDim1(); i++)
		for (index_t j = 0; j < m_rXArray2D.GetDim2(); j++)
		{
			xaSrc.GetSubarray(i, i + iMaskSizeY, j, j + iMaskSizeX, xaTemp);
			//@@@@@@@ call to "std::sort" below should be replaced by smth more efficient (see std algorithms)
			std::sort(xaTemp.begin(), xaTemp.end());
			m_rXArray2D[i][j] = xaTemp[iOffsetY][iOffsetX];
		}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::MedianFilt(index_t iMaskSizeY, index_t iMaskSizeX)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFilt (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::MedianFilt(index_t iMaskSizeY, index_t iMaskSizeX)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFilt (not defined for dcomplex data)");
}

//! Applies median filter with a rectangular mask (only to pixels with the values below and equal to low threshold and above the high threshold)
template<class T> void XArray2DFilt<T>::MedianFiltThresh(index_t iMaskSizeY, index_t iMaskSizeX, T LowThresh, T HighThresh)
{
	if (iMaskSizeY > m_rXArray2D.GetDim1() || iMaskSizeX > m_rXArray2D.GetDim2() || iMaskSizeY < 1 || iMaskSizeX < 1)
		throw std::invalid_argument("invalid argument 'iMaskSizeY or iMaskSizeX' in XArray2DFilt<T>::MedianFiltThresh");

	index_t iOffsetY = iMaskSizeY / 2;
	index_t iOffsetY1 = iMaskSizeY - iMaskSizeY / 2 - 1;
	index_t iOffsetX = iMaskSizeX / 2;
	index_t iOffsetX1 = iMaskSizeX - iMaskSizeX / 2 - 1;
	XArray2D<T> xaSrc(m_rXArray2D); // create a temp copy
	XArray2DMove<T> xaSrcMove(xaSrc);
	xaSrcMove.PadMirror(iOffsetY, iOffsetY1, iOffsetX, iOffsetX1);

	XArray2D<T> xaTemp(iMaskSizeY, iMaskSizeX);

	for (index_t i = 0; i < m_rXArray2D.GetDim1(); i++)
		for (index_t j = 0; j < m_rXArray2D.GetDim2(); j++)
		{
			if (m_rXArray2D[i][j] > LowThresh && m_rXArray2D[i][j] <= HighThresh) continue;

			xaSrc.GetSubarray(i, i + iMaskSizeY, j, j + iMaskSizeX, xaTemp);
			std::sort(xaTemp.begin(), xaTemp.end());
			m_rXArray2D[i][j] = xaTemp[iOffsetY][iOffsetX];
		}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::MedianFiltThresh(index_t iMaskSizeY, index_t iMaskSizeX, fcomplex LowThresh, fcomplex HighThresh)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltThresh (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::MedianFiltThresh(index_t iMaskSizeY, index_t iMaskSizeX, dcomplex LowThresh, dcomplex HighThresh)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltThresh (not defined for dcomplex data)");
}

//! Applies 1D median filter (only to pixels with the values below and equal to low threshold and above the high threshold)
template<class T> void XArray2DFilt<T>::MedianFiltThresh1D(index_t iMaskSize, T LowThresh, T HighThresh, bool bRowWise)
{
	if (bRowWise && iMaskSize > m_rXArray2D.GetDim2() || !bRowWise && iMaskSize > m_rXArray2D.GetDim1() || iMaskSize < 3)
		throw std::invalid_argument("invalid argument 'iMaskSize' in XArray2DFilt<T>::MedianFiltThresh1D");

	index_t iOffset = iMaskSize / 2; // left/top offset
	index_t iOffset1 = iMaskSize - iMaskSize / 2 - 1; // right/bottom offset

	index_t dim1 = m_rXArray2D.GetDim1();
	index_t dim2 = m_rXArray2D.GetDim2();

	// process individual rows/columns 
	if (bRowWise) // this is equivalent to calling MedianFiltThresh(1, iMaskSize, LowThresh, HighThresh) but is more memory efficient
	{		
		// temporary row (mirror padded)
		xar::XArray1D<T> xaRow(dim2 + iMaskSize - 1);
		XArray1D<T> xaTemp(iMaskSize);

		for (index_t i = 0; i < dim1; i++)
		{
			for (index_t j = 0; j < iOffset; j++)
				xaRow[j] = m_rXArray2D[i][iOffset - j]; // mirror-padding at the left boundary
			for (index_t j = 0; j < dim2; j++)
				xaRow[j + iOffset] = m_rXArray2D[i][j];
			for (index_t j = 0; j < iOffset1; j++)
				xaRow[j + iOffset + dim2] = m_rXArray2D[i][dim2 - 2 - j]; // mirror-padding at the right boundary

			for (index_t j = 0; j < dim2; j++)
			{
				if (m_rXArray2D[i][j] > LowThresh && m_rXArray2D[i][j] <= HighThresh) continue;

				xaRow.GetSubarray(j, j + iMaskSize, xaTemp);
				std::sort(xaTemp.begin(), xaTemp.end());
				m_rXArray2D[i][j] = xaTemp[iOffset];
			}
		}
	}
	else // this is equivalent to calling MedianFiltThresh(iMaskSize, 1, LowThresh, HighThresh) but is more memory efficient
	{
		// temporary column (mirror padded)
		xar::XArray1D<T> xaCol(dim1 + iMaskSize - 1);
		XArray1D<T> xaTemp(iMaskSize);

		for (index_t j = 0; j < dim2; j++)
		{
			for (index_t i = 0; i < iOffset; i++)
				xaCol[i] = m_rXArray2D[iOffset - i][j]; // mirror-padding at the top boundary
			for (index_t i = 0; i < dim1; i++)
				xaCol[i + iOffset] = m_rXArray2D[i][j];
			for (index_t i = 0; i < iOffset1; i++)
				xaCol[i + iOffset + dim1] = m_rXArray2D[dim1 - 2 - i][j]; // mirror-padding at the bottom boundary

			for (index_t i = 0; i < dim1; i++)
			{
				if (m_rXArray2D[i][j] > LowThresh && m_rXArray2D[i][j] <= HighThresh) continue;

				xaCol.GetSubarray(i, i + iMaskSize, xaTemp);
				std::sort(xaTemp.begin(), xaTemp.end());
				m_rXArray2D[i][j] = xaTemp[iOffset];
			}
		}
	}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::MedianFiltThresh1D(index_t iMaskSize, fcomplex LowThresh, fcomplex HighThresh, bool bRowWise)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltThresh1D (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::MedianFiltThresh1D(index_t iMaskSize, dcomplex LowThresh, dcomplex HighThresh, bool bRowWise)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltThresh1D (not defined for dcomplex data)");
}

//! Applies median filter with a rectangular mask (only to pixels with non-zero mask values)
template<class T> void XArray2DFilt<T>::MedianFiltMask(index_t iMaskSizeY, index_t iMaskSizeX, XArray2D<char>& rxaMask)
{
	if (iMaskSizeY > m_rXArray2D.GetDim1() || iMaskSizeX > m_rXArray2D.GetDim2() || iMaskSizeY < 1 || iMaskSizeX < 1)
		throw std::invalid_argument("invalid argument 'iMaskSizeY or iMaskSizeX' in XArray2DFilt<T>::MedianFiltMask");

	index_t ny = m_rXArray2D.GetDim1();
	index_t nx = m_rXArray2D.GetDim2();

	if (rxaMask.GetDim1() != ny || rxaMask.GetDim2() != nx)
		throw std::invalid_argument("invalid argument 'rxaMask' in XArray2DFilt<T>::MedianFiltMask (the array and mask haved ifferent dimensions)");

	index_t iOffsetY = iMaskSizeY / 2;
	index_t iOffsetY1 = iMaskSizeY - iMaskSizeY / 2 - 1;
	index_t iOffsetX = iMaskSizeX / 2;
	index_t iOffsetX1 = iMaskSizeX - iMaskSizeX / 2 - 1;
	XArray2D<T> xaSrc(m_rXArray2D); // create a temp copy
	XArray2DMove<T> xaSrcMove(xaSrc);
	xaSrcMove.PadMirror(iOffsetY, iOffsetY1, iOffsetX, iOffsetX1);

	XArray2D<T> xaTemp(iMaskSizeY, iMaskSizeX);

	for (index_t i = 0; i < m_rXArray2D.GetDim1(); i++)
		for (index_t j = 0; j < m_rXArray2D.GetDim2(); j++)
		{
			if (rxaMask[i][j])
			{
				xaSrc.GetSubarray(i, i + iMaskSizeY, j, j + iMaskSizeX, xaTemp);
				std::sort(xaTemp.begin(), xaTemp.end());
				m_rXArray2D[i][j] = xaTemp[iOffsetY][iOffsetX];
			}
		}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::MedianFiltMask(index_t iMaskSizeY, index_t iMaskSizeX, XArray2D<char>& rxaMask)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltMask (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::MedianFiltMask(index_t iMaskSizeY, index_t iMaskSizeX, XArray2D<char>& rxaMask)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltMask (not defined for dcomplex data)");
}

//! Applies 1D median filter (only to pixels with non-zero mask values)
template<class T> void XArray2DFilt<T>::MedianFiltMask1D( index_t iMaskSize, XArray2D<char>& rxaMask, bool bRowWise )
{
	if( bRowWise && iMaskSize > m_rXArray2D.GetDim2() || !bRowWise && iMaskSize > m_rXArray2D.GetDim1() || iMaskSize < 3 )
		throw std::invalid_argument("invalid argument 'iMaskSize' in XArray2DFilt<T>::MedianFiltThresh1D");

	index_t ny( m_rXArray2D.GetDim1() );
	index_t nx( m_rXArray2D.GetDim2() );

	if( rxaMask.GetDim1() != ny || rxaMask.GetDim2() != nx )
		throw std::invalid_argument("invalid argument 'rxaMask' in XArray2DFilt<T>::MedianFiltMask1D (the array and mask haved ifferent dimensions)");

	index_t iOffset( iMaskSize / 2 ); // left/top offset
	index_t iOffset1( iMaskSize - iMaskSize / 2 - 1 ); // right/bottom offset

	index_t dim1( m_rXArray2D.GetDim1() );
	index_t dim2( m_rXArray2D.GetDim2() );

	// process individual rows/columns 
	if( bRowWise ) // this is equivalent to calling MedianFiltMask(1, iMaskSize, xaMask) but is more memory efficient
	{		
		// temporary row (mirror padded)
		xar::XArray1D<T> xaRow(dim2 + iMaskSize - 1);
		XArray1D<T> xaTemp( iMaskSize );

		for (index_t i = 0; i < dim1; i++)
		{
			for (index_t j = 0; j < iOffset; j++)
				xaRow[j] = m_rXArray2D[i][iOffset - j]; // mirror-padding at the left boundary
			for (index_t j = 0; j < dim2; j++)
				xaRow[j + iOffset] = m_rXArray2D[i][j];
			for (index_t j = 0; j < iOffset1; j++)
				xaRow[j + iOffset + dim2] = m_rXArray2D[i][dim2 - 2 - j]; // mirror-padding at the right boundary

			for (index_t j = 0; j < dim2; j++)
			{
				if( rxaMask[i][j] )
				{
					xaRow.GetSubarray( j, j + iMaskSize, xaTemp );
					std::sort( xaTemp.begin(), xaTemp.end() );
					m_rXArray2D[i][j] = xaTemp[iOffset];
				}
			}
		}
	}
	else // this is equivalent to calling MedianFiltMask(iMaskSize, 1, rxaMask) but is more memory efficient
	{
		// temporary column (mirror padded)
		xar::XArray1D<T> xaCol(dim1 + iMaskSize - 1);
		XArray1D<T> xaTemp(iMaskSize);

		for (index_t j = 0; j < dim2; j++)
		{
			for (index_t i = 0; i < iOffset; i++)
				xaCol[i] = m_rXArray2D[iOffset - i][j]; // mirror-padding at the top boundary
			for (index_t i = 0; i < dim1; i++)
				xaCol[i + iOffset] = m_rXArray2D[i][j];
			for (index_t i = 0; i < iOffset1; i++)
				xaCol[i + iOffset + dim1] = m_rXArray2D[dim1 - 2 - i][j]; // mirror-padding at the bottom boundary

			for (index_t i = 0; i < dim1; i++)
			{
				if( rxaMask[i][j] )
				{
					xaCol.GetSubarray( i, i + iMaskSize, xaTemp );
					std::sort( xaTemp.begin(), xaTemp.end() );
					m_rXArray2D[i][j] = xaTemp[iOffset];
				}
			}
		}
	}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::MedianFiltMask1D(index_t iMaskSize, XArray2D<char>& rxaMask, bool bRowWise)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltMask1D (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::MedianFiltMask1D(index_t iMaskSize, XArray2D<char>& rxaMask, bool bRowWise)
{
	throw std::invalid_argument("invalid_argument '*this' in MedianFiltMask1D (not defined for dcomplex data)");
}

//! Applies a ring filter of a specified size (must be odd)
template<class T> void XArray2DFilt<T>::RingFilt(index_t nFiltSize, bool fTransposedSino)
{
	XArray1D<T> arrAvg, arrFilter;
	long nAngles, nPixels;
	long nCurAngle, nCurPixel, i;
	long box;

	if(nFiltSize %2 !=0) nFiltSize++; //FiltSize must be odd (so it can be centred on a pixel)
	box = (long) (nFiltSize--)/2;	

	if(fTransposedSino)
	{
		nAngles = (long) m_rXArray2D.GetDim2();
		nPixels = (long) m_rXArray2D.GetDim1();
	}
	else
	{
		nAngles = (long) m_rXArray2D.GetDim1();
		nPixels = (long) m_rXArray2D.GetDim2();
	}

	// allocate the temp arrays
	arrAvg.Resize(nPixels, 0);
	arrFilter.Resize(nPixels, 0);

	//sum the array over angles and put result in arrAvg
	for(nCurAngle = 0; nCurAngle < nAngles; nCurAngle++)
	{
		for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
		{
			if(fTransposedSino)
				arrAvg[nCurPixel] += m_rXArray2D[nCurPixel][nCurAngle];
			else
				arrAvg[nCurPixel] += m_rXArray2D[nCurAngle][nCurPixel];
		}
	}
	// divide by the number of angles to get the average
	for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
		arrAvg[nCurPixel] = (T) (((double) arrAvg[nCurPixel]) / nAngles);

	//smooth results and copy into zeroth row of arrTemp
	for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
	{
		double sum = 0;
		int count = 0;
		
		for(i = -box; i <= box; i++)
		{
			int pixel = nCurPixel + i;
			if(pixel >= 0 && pixel < nPixels)
			{
				sum += arrAvg[pixel];
				count++;
			}
		}
		arrFilter[nCurPixel] = (T) (sum / (T) count);
	}

	//put diff between smooth and non-smooth average into arrFilter
	//This is the filter
	for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
		arrFilter[nCurPixel] = arrAvg[nCurPixel] - arrFilter[nCurPixel];
		
	//subtract arrFilter filter from the sinogram
	for(nCurAngle = 0; nCurAngle < nAngles; nCurAngle++)
	{
		for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
		{
			if(fTransposedSino)
				m_rXArray2D[nCurPixel][nCurAngle] -= arrFilter[nCurPixel];
			else
				m_rXArray2D[nCurAngle][nCurPixel] -= arrFilter[nCurPixel];
		}
	}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::RingFilt(index_t nFiltSize, bool fTransposedSino)
{
	throw std::invalid_argument("invalid_argument '*this' in RingFilt (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::RingFilt(index_t nFiltSize, bool fTransposedSino)
{
	throw std::invalid_argument("invalid_argument '*this' in RingFilt (not defined for dcomplex data)");
}

//! Applies a zingers filter of a specified size (must be odd, by default 9) and threshold (by default 1.2)
//	Processing is carried out for each individual detector row (if bTransposed == false) or column (if bTransposed == true)
//	To filter sinograms, apply exp(-) before calling this function, and -log() after
template<class T> void XArray2DFilt<T>::ZingersFilt(index_t nFiltSize, double dblThreshold, bool bTransposed)
{
	XArray1D<T> arrAvg;
	XArray1D<double> arrFilter;
	index_t nAngles, nPixels;
	index_t nCurAngle, nCurPixel;
	long box, i, pixel;

	if(nFiltSize%2 != 0) nFiltSize++; //FiltSize must be odd (so it can be centred on a pixel)
	box = (long) (nFiltSize--)/2;	
	if(bTransposed) { nAngles = m_rXArray2D.GetDim2(); nPixels = m_rXArray2D.GetDim1(); }
	else { nAngles = m_rXArray2D.GetDim1(); nPixels = m_rXArray2D.GetDim2(); }

	// allocate the temp arrays
	arrAvg.Resize(nPixels + 4);
	arrFilter.Resize(nPixels + 4);

	// sequentially process individual rows/columns in arrAvg
	for(nCurAngle = 0; nCurAngle < nAngles; nCurAngle++)
	{
		// select one row/column
		if(bTransposed)
			for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
				arrAvg[nCurPixel + 2] = m_rXArray2D[nCurPixel][nCurAngle];
		else
			for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
				arrAvg[nCurPixel + 2] = m_rXArray2D[nCurAngle][nCurPixel];

		arrAvg[0] = arrAvg[3]; arrAvg[1] = arrAvg[2];
		arrAvg[nPixels + 2] = arrAvg[nPixels + 1]; arrAvg[nPixels + 3] = arrAvg[nPixels];

		// smooth the row/column
		for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
		{
			double sum = 0;
			index_t count = 0;
			for(i = -box; i <= box; i++)
			{
				pixel = long(nCurPixel) + i;
				if(pixel >= 0 && pixel < nPixels) { sum += arrAvg[pixel + 2]; count++; }
			}
			arrFilter[nCurPixel + 2] = (sum / count);
		}

		// divide original row/column by the smoothed row/column and save in arrFilter
		for(nCurPixel = 2; nCurPixel < nPixels + 2; nCurPixel++)
			arrFilter[nCurPixel] = arrAvg[nCurPixel] / arrFilter[nCurPixel];

		arrFilter[0] = arrFilter[3]; arrFilter[1] = arrFilter[2];
		arrFilter[nPixels + 2] = arrFilter[nPixels + 1]; arrFilter[nPixels + 3] = arrFilter[nPixels];

		// threshold the arrFilter values
		for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
		{
			if (arrFilter[nCurPixel + 2] >= dblThreshold)
			{
//				arrAvg[nCurPixel + 2] = (arrAvg[nCurPixel] + arrAvg[nCurPixel + 4]) / T(2);
				double sum = 0;
				index_t count = 0;
				if (arrFilter[nCurPixel] < dblThreshold) { sum += arrAvg[nCurPixel]; count++; }
				if (arrFilter[nCurPixel + 1] < dblThreshold) { sum += arrAvg[nCurPixel + 1]; count++; }
				if (arrFilter[nCurPixel + 3] < dblThreshold) { sum += arrAvg[nCurPixel + 3]; count++; }
				if (arrFilter[nCurPixel + 4] < dblThreshold) { sum += arrAvg[nCurPixel + 4]; count++; }
				if (count > 0) arrAvg[nCurPixel + 2] = T(sum / count);
			}
		}

		// overwrite the original row/column with the filtered one
		if(bTransposed)
			for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
				m_rXArray2D[nCurPixel][nCurAngle] = arrAvg[nCurPixel + 2];
		else
			for(nCurPixel = 0; nCurPixel < nPixels; nCurPixel++)
				m_rXArray2D[nCurAngle][nCurPixel] = arrAvg[nCurPixel + 2];
	}
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::ZingersFilt(index_t nFiltSize, double dblThreshold, bool bTransposed)
{
	throw std::invalid_argument("invalid_argument '*this' in ZingersFilt (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::ZingersFilt(index_t nFiltSize, double dblThreshold, bool bTransposed)
{
	throw std::invalid_argument("invalid_argument '*this' in ZingersFilt (not defined for dcomplex data)");
}

//!	Applies a filter that removes vertical streaks
//	dblThresh is a threshold (the bigger the threshold the stronger is filtering; too big thrashold can result in artifacts)
//	As a rule of thumb, the threshold value should be chosen to be about the std of the streaks
//	iWidth is a width of the filter's window
template<class T> void XArray2DFilt<T>::VertStreaksFilt(double dblThresh, index_t iWidth)
{
	index_t ny = m_rXArray2D.GetDim1();
	index_t nx = m_rXArray2D.GetDim2();

	XArray2D<T> xaTemp;
	XArray2D<T> xaA(ny, iWidth);
	XArray2D<long> xaK(nx - iWidth, iWidth, 0);
	XArray2D<T> xaa(nx - iWidth, iWidth, 0);
	for (index_t i = 0; i < nx - iWidth; i++)
	{
		m_rXArray2D.GetSubarray(0, ny, i, i + iWidth, xaTemp);
		index_t counter = 0;
		for (index_t j = 0; j < ny; j++)
		{
			// calculate variance for each row
			double dblDisp = 0;
			double dblAver = 0;
			for (index_t k = 0; k < iWidth; k++)
			{
				double dblValue = double(xaTemp[j][k]);
				dblAver += dblValue;
				dblDisp += dblValue * dblValue;
			}
			dblAver /= iWidth;
			dblDisp /= iWidth;
			double dblVar = dblDisp - dblAver * dblAver;
			dblVar = sqrt(dblVar);

			if (dblVar < dblThresh)
			{
				for (index_t k = 0; k < iWidth; k++)
					xaA[counter][k] = xaTemp[j][k] - T(dblAver);
				counter++;	
			}
		}

		index_t iA = counter;
		if (iA > 0)
		{
			for (index_t k = 0; k < iWidth; k++)
			{
				xaK[i][k] = T(iA);
				double dblMean = 0;
				for (index_t j = 0; j < iA; j++) dblMean += xaA[j][k];
				xaa[i][k] = T(dblMean / iA);
			}
		}
	}

	std::vector<T> va(nx, 0);
	std::vector<long> vK(nx, 0);
#if(1)
	for (index_t n = 0; n < nx - iWidth; n++)
		for (index_t l = 0; l < iWidth; l++)
			if (xaK[n][l] > vK[n + l])
			{
				vK[n + l] = xaK[n][l];
				va[n + l] = xaa[n][l];
			}
#else	
	for (index_t n = 0; n < nx - iWidth; n++) {
		for (index_t l = 0; l < iWidth; l++) {
			if (xaK[n][l] > 0)
			{
				vK[n + l] += xaK[n][l];
				va[n + l] += xaa[n][l] * xaK[n][l];
			}
		}
	}
	for (index_t i = 0; i < nx; i++)
		if (vK[i] > 0) va[i] /= vK[i];
#endif

	for (index_t i = 0; i < ny; i++)
		for (index_t j = 0; j < nx; j++)
			m_rXArray2D[i][j] -= va[j];
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::VertStreaksFilt(double dblThresh, index_t iWidth)
{
	throw std::invalid_argument("invalid_argument '*this' in ZingersFilt (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::VertStreaksFilt(double dblThresh, index_t iWidth)
{
	throw std::invalid_argument("invalid_argument '*this' in ZingersFilt (not defined for dcomplex data)");
}

template<class T> void XArray2DFilt<T>::TVFilt(
	double dblSTD,
	double dblBetaRel,
	double dblSupportRadius,
	T dblMuStep,
	size_t numMu, 
	T dblMuMin,
	size_t numIter,
	std::vector<double>& rvChi2,
	std::vector<double>& rvRegul)
{
	xar::XArray2D<T> xaSlice(m_rXArray2D);

	///////////////////
	// Slice support //
	///////////////////
	index_t ny = xaSlice.GetDim1();
	index_t nx = xaSlice.GetDim2();
	double dblAxisJ = 0.5 * (nx - 1);//???
	xar::XArray2D<char> xaSupport(ny, nx, 0);
	double dblCentreX = dblAxisJ;
	double dblCentreY = 0.5 * (ny - 1);
//	double dblAxisJp = dblAxisJ - 0.5 * (nx - 1);
//	double x0dxstP = __min(0.5 * (nx - 1) - fabs(dblAxisJp), dblCentreY); // radius of the reconstructable circle (in pixels)
//	if (dblSupportRadius <= 0 || dblSupportRadius > x0dxstP) dblSupportRadius = x0dxstP;
	index_t numSupport;
	if (dblSupportRadius <= 0) {
		numSupport = nx * ny;
		xaSupport.Fill(1);
	}
	else {
		double dblRadius2 = dblSupportRadius * dblSupportRadius;
		numSupport = 0;
		for (index_t i = 0; i < ny; i++) {
			double y = i - dblCentreY;
			double y2 = y * y;
			for (index_t j = 0; j < nx; j++) {
				double x = j - dblCentreX;
				if (y2 + x * x < dblRadius2) {
					xaSupport[i][j] = 1;
					numSupport++;
				}
				//else xaSlice[i][j] = 0.f;
			}
		}
	}

	std::vector<size_t> GL_IDX(numSupport);
	size_t address = 0;
	for (index_t i = 0; i < ny; i++){
		for (index_t j = 0; j < nx; j++){
			if (xaSupport[i][j]) { GL_IDX[address] = i * nx + j; address++; }
		}
	}

#ifdef HAVE_STD_UNIFORM_DIST
	// Setup Mersenne twister engine & initialise from a random device
	std::random_device rd;
	std::mt19937 gen( rd() );
	std::uniform_real_distribution<> distUniform;
#endif
#ifdef HAVE_BOOST_UNIFORM_DIST
	boost::random::mt19937 gen;
	boost::random::uniform_real_distribution<> distUniform;
#endif

	for (index_t i = 0; i < GL_IDX.size() * 2; i++)
	{
		size_t i1 = size_t(distUniform(gen) * (GL_IDX.size() - 1) + 0.5);
		size_t i2 = size_t(distUniform(gen) * (GL_IDX.size() - 1) + 0.5);

		std::swap(GL_IDX[i1], GL_IDX[i2]);
	}

	///////////////////////////////////////////
	// Calculate chi^2/2 and total variation //
	///////////////////////////////////////////
	double chi2 = Chi2a(xaSlice, m_rXArray2D, dblSTD);
	double regul = Regul(xaSlice, 1);

	rvChi2.resize(numIter + 1);
	rvChi2[0] = chi2;
	rvRegul.resize(numIter + 1);
	rvRegul[0] = regul;

	double dblBeta = dblBetaRel / regul;

	/////////////////////////////////
	// Vector of deviations for mu //
	/////////////////////////////////
	std::vector<T> vMu(numMu);
	for (index_t i = 0; i < vMu.size(); i++) vMu[i] = (T)(-0.5 * dblMuStep * (numMu - 1) + dblMuStep * i);

	T MuMax = vMu[numMu - 1];
	T fltThresh = dblMuMin - MuMax;
	for (index_t i = 0; i < ny; i++)
		for (index_t j = 0; j < nx; j++)
			if (xaSlice[i][j] < fltThresh) xaSlice[i][j] = fltThresh;

	// Main loop
	for( index_t i = 1; i <= numIter; i++ )
	{
		OptimisationCPU(
			xaSlice,
			m_rXArray2D,
			dblSTD,
			GL_IDX,
			vMu,
			dblBeta,
			dblMuMin);

		rvChi2[i] = Chi2a(xaSlice, m_rXArray2D, dblSTD);
		rvRegul[i] = Regul(xaSlice, 1);
	}

	m_rXArray2D = xaSlice;
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::TVFilt(
	double dblSTD,
	double dblBetaRel,
	double dblSupportRadius,
	fcomplex dblMuStep,
	size_t numMu, 
	fcomplex dblMuMin,
	size_t numIter,
	std::vector<double>& rvChi2,
	std::vector<double>& rvRegul)
{
	throw std::invalid_argument("invalid_argument '*this' in TVFilt (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::TVFilt(
	double dblSTD,
	double dblBetaRel,
	double dblSupportRadius,
	dcomplex dblMuStep,
	size_t numMu, 
	dcomplex dblMuMin,
	size_t numIter,
	std::vector<double>& rvChi2,
	std::vector<double>& rvRegul)
{
	throw std::invalid_argument("invalid_argument '*this' in TVFilt (not defined for dcomplex data)");
}

template<class T> void XArray2DFilt<T>::OptimisationCPU(
	xar::XArray2D<T>& rxaSlice,
	const xar::XArray2D<T>& rxaNoisySlice,
	double dblSTD,
	const std::vector<size_t>& rvIDX,
	const std::vector<T>& rvMu,
	double dblBeta,
	T dblMuMin)
{
	index_t ny( rxaSlice.GetDim1() );
	index_t nx( rxaSlice.GetDim2() );
//	if (nx != ny)
//		throw std::runtime_error("runtime error in XArray2DFilt<T>::OptimisationCPU() (input slice MUST be square!!!)");
	if (nx == 0 || ny == 0)
		throw std::runtime_error("runtime error in XArray2DFilt<T>::OptimisationCPU() (input slice MUST have non-zero size!!!)");

	// Number of descrete values of mu
	index_t numMu = rvMu.size();
	// Number of pixels to be varied
	index_t numVar = rvIDX.size();

	// Main loop over the pixels
	//		int numProc = omp_get_num_procs();
	//#pragma omp parallel for 
	//num_threads(numProc)
	for( int iVar = 0; iVar < numVar; iVar++ )
	{
		size_t gl_idx = rvIDX[iVar];
		index_t nJ = gl_idx%nx;
		index_t nI = (gl_idx - nJ) / nx;
		if (nI == 0 || nI >= ny - 1 || nJ == 0 || nJ >= nx - 1) continue;

		// Initialise a vector of regularisation function increments
		std::vector<double> vSum(numMu, 0);

		T fltMuCtr = rxaSlice[nI][nJ];

		//////////////////////////////////////////////
		// CALCULATE THE REGULARISATION TERM CHANGE //
		//////////////////////////////////////////////
		if (dblBeta > 0)
		{
			// Contribution of the chosen pixel
			T dblBttm = rxaSlice[nI + 1][nJ];
			T dblRght = rxaSlice[nI][nJ + 1];
			T dblCntr = rxaSlice[nI][nJ];
			double dblTemp = sqrt(double(dblBttm - dblCntr) * double(dblBttm - dblCntr) + double(dblRght - dblCntr) * double(dblRght - dblCntr));
			for (index_t iMu = 0; iMu < numMu; iMu++)
			{
				T dblCntrVar = dblCntr + rvMu[iMu];
#if(NEG)
				if (dblCntrVar < dblMuMin) continue;
#endif
				double dblTest = sqrt(double(dblBttm - dblCntrVar) * double(dblBttm - dblCntrVar) + double(dblRght - dblCntrVar) * double(dblRght - dblCntrVar));
				vSum[iMu] = (dblTest - dblTemp);
			}
			// Contribution of the top pixel
			dblBttm = rxaSlice[nI][nJ];
			dblRght = rxaSlice[nI - 1][nJ + 1];
			dblCntr = rxaSlice[nI - 1][nJ];
			dblTemp = sqrt(double(dblBttm - dblCntr) * double(dblBttm - dblCntr) + double(dblRght - dblCntr) * double(dblRght - dblCntr));
			for (index_t iMu = 0; iMu < numMu; iMu++)
			{
				T dblBttmVar = dblBttm + rvMu[iMu];
#if(NEG)
				if (dblBttmVar < dblMuMin) continue;
#endif
				double dblTest = sqrt(double(dblBttmVar - dblCntr) * double(dblBttmVar - dblCntr) + double(dblRght - dblCntr) * double(dblRght - dblCntr));
				vSum[iMu] += (dblTest - dblTemp);
			}
			// Contribution of the left pixel
			dblBttm = rxaSlice[nI + 1][nJ - 1];
			dblRght = rxaSlice[nI][nJ];
			dblCntr = rxaSlice[nI][nJ - 1];
			dblTemp = sqrt(double(dblBttm - dblCntr) * double(dblBttm - dblCntr) + double(dblRght - dblCntr) * double(dblRght - dblCntr));
			for( index_t iMu = 0; iMu < numMu; iMu++ )
			{
				T dblRghtVar = dblRght + rvMu[iMu];
#if(NEG)
				if (dblRghtVar < dblMuMin) continue;
#endif
				double dblTest = sqrt(double(dblBttm - dblCntr) * double(dblBttm - dblCntr) + double(dblRghtVar - dblCntr) * double(dblRghtVar - dblCntr));
				vSum[iMu] += (dblTest - dblTemp);

				vSum[iMu] *= (nx * nx * dblBeta);
			}
		}

		///////////////////////////////////////////////
		// CALCULATE THE GOODNESS-OF-FIT TERM CHANGE //
		///////////////////////////////////////////////
		double norm = (1. / (dblSTD * dblSTD * nx * ny));
		T fltNoisyPixel = rxaNoisySlice[nI][nJ];
		double fltDelta = (double)(fltMuCtr - fltNoisyPixel);
		double fltDelta2 = fltDelta * fltDelta;	
		for (index_t iMu = 0; iMu < numMu; iMu++)
		{
			T fltTemp = fltMuCtr + rvMu[iMu];
#if(NEG)
			if (fltTemp < dblMuMin) { vSum[iMu] = 1e38; continue; }
#endif
			double fltdPtest = (double)(fltTemp - fltNoisyPixel);
			vSum[iMu] += 0.5 * (fltdPtest * fltdPtest - fltDelta2) * norm;
		}

		// FIND such Mu that results in minimum Sum
		index_t idx = 0;
		double minimum = 1e10;
		for (index_t iMu = 0; iMu < numMu; iMu++){
			if (vSum[iMu] < minimum)
			{
				minimum = vSum[iMu];
				idx = iMu;
			}
		}
		rxaSlice[nI][nJ] = fltMuCtr + rvMu[idx];
	} // iVar
}

// Specialization for T=fcomplex
template<> inline void XArray2DFilt<fcomplex>::OptimisationCPU(
	xar::XArray2D<fcomplex>& rxaSlice,
	const xar::XArray2D<fcomplex>& rxaNoisySlice,
	double dblSTD,
	const std::vector<size_t>& rvIDX,
	const std::vector<fcomplex>& rvMu,
	double dblBeta,
	fcomplex dblMuMin)
{
	throw std::invalid_argument("invalid_argument '*this' in OptimisationCPU (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline void XArray2DFilt<dcomplex>::OptimisationCPU(
	xar::XArray2D<dcomplex>& rxaSlice,
	const xar::XArray2D<dcomplex>& rxaNoisySlice,
	double dblSTD,
	const std::vector<size_t>& rvIDX,
	const std::vector<dcomplex>& rvMu,
	double dblBeta,
	dcomplex dblMuMin)
{
	throw std::invalid_argument("invalid_argument '*this' in OptimisationCPU (not defined for dcomplex data)");
}

template<class T> double XArray2DFilt<T>::Regul(
	xar::XArray2D<T>& rxaSlice,
	unsigned int iNorm)
{
	index_t ny = rxaSlice.GetDim1();
	index_t nx = rxaSlice.GetDim2();
//	if (nx != ny)
//		throw std::runtime_error("runtime error in XArray2DFilt<T>::Regul() (input slice MUST be square!!!)");
	if (nx == 0 || ny == 0)
		throw std::runtime_error("runtime error in XArray2DFilt<T>::Regul() (input slice MUST have non-zero size!!!)");

	// CALCULATE THE REGULARISATION TERM
	double 	dblReg = 0.;

	switch (iNorm)
	{
	case 0: // Penalized Weighted Least Squares
		for (index_t i = 1; i < ny - 1; i++){
			for (index_t j = 1; j < nx - 1; j++){
				double dblMuSlice = (double)rxaSlice[i][j];
				double dblTemp = (double)rxaSlice[i - 1][j - 1] - dblMuSlice;
				dblReg += 0.707 * dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i - 1][j] - dblMuSlice;
				dblReg += dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i - 1][j + 1] - dblMuSlice;
				dblReg += 0.707 * dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i][j - 1] - dblMuSlice;
				dblReg += dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i][j + 1] - dblMuSlice;
				dblReg += dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i + 1][j - 1] - dblMuSlice;
				dblReg += 0.707 * dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i + 1][j] - dblMuSlice;
				dblReg += dblTemp * dblTemp;
				dblTemp = (double)rxaSlice[i + 1][j + 1] - dblMuSlice;
				dblReg += 0.707 * dblTemp * dblTemp;
			}
		}
		dblReg *= (nx * ny);
		break;
	case 1: // Total Variation
		for (index_t i = 1; i < ny - 1; i++){
			for (index_t j = 1; j < nx - 1; j++){
				double dblCntr = (double)rxaSlice[i][j];
				double dblBttm = (double)rxaSlice[i + 1][j];
				double dblRght = (double)rxaSlice[i][j + 1];
				dblReg += sqrt((dblBttm - dblCntr) * (dblBttm - dblCntr) + (dblRght - dblCntr) * (dblRght - dblCntr));
			}
		}
		dblReg *= (nx * ny);
		break;
	case 2: // Total Boundary
		for (index_t i = 1; i < ny - 1; i++){
			for (index_t j = 1; j < nx - 1; j++){
				double dblMuSlice = (double)rxaSlice[i][j];
				// Top boundary
				dblReg += NonZero((double)rxaSlice[i - 1][j] - dblMuSlice);
				// Left boundary
				dblReg += NonZero((double)rxaSlice[i][j - 1] - dblMuSlice);
				// Right boundary
				dblReg += NonZero((double)rxaSlice[i][j + 1] - dblMuSlice);
				// Bottom boundary
				dblReg += NonZero((double)rxaSlice[i + 1][j] - dblMuSlice);
			}
		}
		break;
	default:
		throw std::invalid_argument("invalid_argument 'iNorm' in XArray2DFilt<T>::Regul() (unsupported norm)"); 
	}

	return dblReg;
}

// Specialization for T=fcomplex
template<> inline double XArray2DFilt<fcomplex>::Regul(
	xar::XArray2D<fcomplex>& rxaSlice,
	unsigned int iNorm)
{
	throw std::invalid_argument("invalid_argument '*this' in Regul (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline double XArray2DFilt<dcomplex>::Regul(
	xar::XArray2D<dcomplex>& rxaSlice,
	unsigned int iNorm)
{
	throw std::invalid_argument("invalid_argument '*this' in Regul (not defined for dcomplex data)");
}

template<class T> double XArray2DFilt<T>::Chi2a(
	const xar::XArray2D<T>& rxaArray1,
	const xar::XArray2D<T>& rxaArray2,
	double dblSTD)
{
	xar::index_t num1 = rxaArray1.GetDim1();
	xar::index_t num2 = rxaArray1.GetDim2();
	if (rxaArray2.GetDim1() != num1 || rxaArray2.GetDim2() != num2)
		throw std::invalid_argument("invalid arguments 'rxaArray1' and 'rxaArray2' in XArray2DFilt<T>::Chi2a() (different dimensions)");
	if (dblSTD <= 0)
		throw std::invalid_argument("invalid argument 'dblSTD' in XArray2DFilt<T>::Chi2a() (must be positive)");

	double chi2 = 0.;
	for (xar::index_t i = 0; i < num1; i++){
		for (xar::index_t j = 0; j < num2; j++){
			double temp = (double)(rxaArray1[i][j] - rxaArray2[i][j]);
			chi2 += temp * temp;
		}
	}

	return (0.5 * chi2 / (dblSTD * dblSTD * num1 * num2));
}

// Specialization for T=fcomplex
template<> inline double XArray2DFilt<fcomplex>::Chi2a(
	const xar::XArray2D<fcomplex>& rxaArray1,
	const xar::XArray2D<fcomplex>& rxaArray2,
	double dblSTD)
{
	throw std::invalid_argument("invalid_argument '*this' in Chi2a (not defined for fcomplex data)");
}

// Specialization for T=dcomplex
template<> inline double XArray2DFilt<dcomplex>::Chi2a(
	const xar::XArray2D<dcomplex>& rxaArray1,
	const xar::XArray2D<dcomplex>& rxaArray2,
	double dblSTD)
{
	throw std::invalid_argument("invalid_argument '*this' in Chi2a (not defined for dcomplex data)");
}

#endif

} // namespace xar closed



//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//


// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray2DFilt<char>;
	template class xar::XArray2DFilt<short>;
	template class xar::XArray2DFilt<long>;
	template class xar::XArray2DFilt<float>;
	template class xar::XArray2DFilt<double>;
	template class xar::XArray2DFilt<xar::fcomplex>;
	template class xar::XArray2DFilt<xar::dcomplex>;
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
