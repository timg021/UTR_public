//Header XA_TIE.h
//
//
//	HEADER FILE TITLE:
//
//		Phase retrieval algorithms based on the TIE approximations
//
//
//
/*!
	\file		XA_TIE.h
	\brief		Phase retrieval algorithms based on the TIE approximations
	\par		Description:
		This header contains a class that provides TIE-based phase retrieval
		services for XArray2D<T> objects
*/

#pragma once

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_move2.h"
#include "OouraFft.h"
//#include "unwrap_2d_ljmu.h"
//#include "XA_fft2.h"
//#include "XA_fftr2.h"
//#include "XA_lin2.h"
//#include "XA_filt2.h"
//#include "XA_Lapl.h"
//#include "XArFileIO_VOD.h"

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
//Class XA_2DTIE<T>
//
//	Phase retrieval algorithms based on the TIE approximations
//
/*!
	\brief		Phase retrieval algorithms based on the TIE approximations
	\par		Description:
				This class template provides TIE-based phase retrieval services
				for XArray2D<T> objects
	\remarks	An object of this class represents an interface exposing several functions
				that provide various TIE-based phase retrieval services				
	\remarks	This class can only be created for T=float or T=double
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XA_2DTIE
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XA_2DTIE() {}
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XA_2DTIE(const XA_2DTIE<T>& rCopy) { GetValuetype(); }
	public:
		//! Destructor
		~XA_2DTIE() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XA_2DTIE<T>& rCopy) {}

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XA_2DTIE<T> objects for unsupported types T
		//! Returns the xar::_eValueType corresponding to T
		static _eValueType GetValuetype(void);
		
	// Operations
	public:
		//! Retrieves phase from 1 image of a pure phase object using TIE approximation
		void TIE1(xar::XArray2D<T>& F, double R, bool UseFFT, long ncycle, double alpha);
		//! Retrieves phase&amplitude from 1 image of a 'single-material' object using TIE approximation
		void DP(xar::XArray2D<T>& F, double deltaoverbeta, double R);
		//! Replaces the modulus of a complex array by values from a new real-amplitude array
		void ReplaceModulus(xar::XArray2D< std::complex<T> >& F, xar::XArray2D<T>& A);
		//! Enforce the homogeneous object condition by changing the phase of the complex amplitude
		void Homogenise(xar::XArray2D< std::complex<T> >& F, double delta2beta);
		//! Enforce the homogeneous object condition by changing the modulus of the complex amplitude
		void Homogenise1(xar::XArray2D< std::complex<T> >& F, double delta2beta, double defocus, char* pinput_mask = 0);
		//! Enforce unit plane wave within a given vicinity of the boundary
		void EnforceSupport(xar::XArray2D< std::complex<T> >& F, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, std::complex<T> tMaskVal);
		//! Solves 2D phase unwrapping problem using unwrap_2d_ljmu.c module
		void CPhaseUnwrap(xar::XArray2D< std::complex<T> >& F, xar::XArray2D<T>& Fi, char* pinput_mask = 0);
		//! Retrieves phase&amplitude from 1 image of a 'single-material' object using TIE approximation and PSF deconvolution
		void DPDeconv(xar::XArray2D<T>& F, xar::XArray2D<T>& Ker, double deltaoverbeta, double R, double alpha);
		//! Retrieves phase&amplitude from 2 images at different defocus distances using TIE approximation
		void TIE2R(const xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, double R, bool UseFFT, long ncycle, double alpha);
		//! Retrieves phase&amplitude from 3 images at different defocus distances using TIE approximation
		void TIE3R(xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, const xar::XArray2D<T>& I2, double R, bool UseFFT, long ncycle, double alpha);
		//! Retrieves phase&amplitude from 2 images at different defocus distances using modified TIE approximation
		void ModTIE2R(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double R, double gamma2, long ncycle, double alpha);
		//! Retrieves phase&amplitude from 3 images at different defocus distances using modified TIE approximation
		void ModTIE3R(const xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, const xar::XArray2D<T>& I2, double R,  short ncycle, double averV);
		//! Retrieves phase&amplitude from 2 images at different wavelengths (energies) using TIE approximation
		void TIE2E(xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, double R, bool UseFFT, long ncycle, double alpha);
		//! Retrieves phase&amplitude from 3 images at different wavelengths (energies) using TIE approximation
		void TIE3E(xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double R, bool UseFFT, long ncycle, double alpha);
		//! Checks the validity of the TIE approximation for given image(s)
		bool TIE_Validity(const xar::XArray2D<T>& I0, const xar::XArray2D<T>& I1, const xar::XArray2D<T>& I2);
		//! Solves 2D phase unwrapping problem using TIE approximation
		void CPhaseTIE(const xar::XArray2D<std::complex<T> >& U, xar::XArray2D<T>& Fi, long ncycle, double alpha);
		//! Solves 2D phase unwrapping problem using TIE-Hom approximation
		void CPhaseTIE1(xar::XArray2D<std::complex<T> >& U, xar::XArray2D<T>& Fi, double defocus, double delta2beta);
		//! Qualitatively enhances phase-contrast images using unsharp masking
		void PhaseEnhancement(xar::XArray2D<T>& F, double deltaoverbeta=-1.0, double R=0.0, long blur_radius_in_pixels=0, double weighting_factor=0.3, double histogram_factor=0.99);
		//! Deblurrs an image of a single-component sample by means of defocus
		void Deblur(xar::XArray2D<T>& F, double deltaoverbeta, double R, double sigma);

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
	// Member functions

	};

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
#if(0)
//! Retrieves phase from 1 image of a pure phase object using TIE approximation
// TIE phase-amplitude retrieval from 1 image
// F - image; on exit it is replaced by the phase
// R - distance R-prime
// if UseFFT is true, then PoissonFFT is used for phase retrieval, else PoissonFMG is used
// alpha and ncycle - obvious meaning in the context of PoissonFFT or PoissonFMG
//
// NOTE!!!: the program assumes that F is background-corrected, i.e. DIVIDED by the 
//  corresponding INCIDENT intensity at wl (F = I/I_inc);
//  it solves equation -Laplacian(Fi)) = k*(F/<F0>-1)/R.
// NOTE!!!: the program also GENERALLY assumes that the object is non-absorbing, i.e. F.Norm(eNormAver)==1;
//  however F.Norm(eNormAver) is calculated in the program and used for normalization
//
// NOTE: PoissonFMG requires square arrays with xst=yst, PoissonFFT does not require that
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
 
template <class T> void XA_2DTIE<T>::TIE1(xar::XArray2D<T>& F, double R, bool UseFFT, long ncycle, double alpha)
{
	if (!GetIXAHWave2D(F))
		throw std::invalid_argument("invalid_argument 'F' in XA_2DTIE<T>::TIE1 (wavehead is missing)");

	 //enforce energy conservation
	double aver = F.Norm(xar::eNormAver);
	F -= T(aver);
	F *= T(xar::tPI / GetIXAHWave2D(F)->GetWl() / R / aver);
	XA_2DLapl<T> xaLapl(F);

	if (UseFFT) 
		xaLapl.PoissonFFT(alpha); 
	else 
		xaLapl.PoissonFMG(ncycle, alpha);
}
#endif

//! Retrieves phase&amplitude from 1 image of a 'single-material' object using TIE approximation
// Implements David Paganin's inversion algorithm applicable to in-line X-ray phase-contrast
// images of samples consisting of single component. Applicability condition is:
// beta(r)=const*delta(r), where n(r)=1-delta-i*beta is the complex refractive index.
// deltaoverbeta is delta/beta
// R is a propagation distance (in units used in F.head - default are microns)
// On exit F is replaced by the reconstructed INTENSITY (= exp(-4.0 * PI / wl * beta * thickness));
//
// NOTE!!!: the program assumes that F is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iinc, 1 must not be subtracted!!!)
// NOTE: this program uses the Ooura FFT library
// NOTE: this progam is a modification of PoissonFFT() - see invLa.h
//
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::DP( xar::XArray2D<T>& F, double deltaoverbeta, double R )
{
	// check input parameters
	if( !GetIXAHWave2D( F ) ) 
		throw std::invalid_argument( "invalid_argument 'F' in XA_2DTIE<T>::DP (wavehead is missing)" );

	GetIXAHWave2D( F )->Validate();
	
	if( deltaoverbeta <= 0 ) 
		throw std::invalid_argument( "invalid_argument 'deltaoverbeta' in XA_2DTIE<T>::DP (delta/beta must be positive)" );

	index_t nxF( F.GetDim1() );
	index_t nx( nxF );
	index_t nyF( F.GetDim2() );
	index_t ny( nyF );
	
	double wl( GetIXAHWave2D( F )->GetWl() ); 
	
	// pad to the nearest power of 2
	index_t i( 2 );
	
	while( i < nx ) 
		i *= 2;

	nx = i;
	index_t j( 2 );

	while( j < ny ) 
		j *= 2;

	ny = j;
	
	if( nx != nxF || ny != nyF ) 
	{
		xar::XArray2DMove<T> tmp(F);
		tmp.Pad( ( nx - nxF ) / 2, nx - nxF - ( nx - nxF ) / 2, ( ny - nyF ) / 2, ny - nyF - ( ny - nyF ) / 2, T( 1 ) );
	}

	T* arr( &F.front() );
	std::vector<T> vecSpeq( nx*2 );
	T* speq( &vecSpeq.front() );

	// forward FFT
	xar::OouraFft<T> fft;

	fft.Real2D( arr, speq, nx, ny, xar::OouraFft<T>::eDirFwd );

	// multiplication of the FFT by 1/(m^2+n^2)
	// the indexing follows that from the function FFTRe(...) - see vo_fftr2.h
	double factor( 2. / ( nx * ny ) );
	double a( GetIXAHWave2D( F )->GetYhi() - GetIXAHWave2D( F )->GetYlo() + GetYStep( F ) );			// X/Y dimensions reversed
	double b( GetIXAHWave2D( F )->GetXhi() - GetIXAHWave2D( F )->GetXlo() + GetXStep( F ) );			// X/Y dimensions reversed
	double a2( xar::PI * wl * R * deltaoverbeta / (a * a) / factor );
	double b2( xar::PI * wl * R * deltaoverbeta / (b * b) / factor );
	double alpha( 1.0 / factor );
	double dtemp;
	T temp;
	index_t nxd2( nx / 2 );
	index_t nyd2( ny / 2 );
	index_t i1;
	index_t j1;

	dtemp = alpha + b2 * nyd2 * nyd2;

	for( i = 0; i < nxd2; i++ ) 
	{ 
		i1 = nxd2 - i;
		temp = T( 1.0 / ( dtemp + a2 * i1 * i1 ) );
		speq[nx+2*i] *= temp;
		speq[nx+2*i+1] *= temp;
	}

	for( i = nxd2; i < nx; i++ ) 
	{ 
		i1 = i - nxd2;
		temp = T( 1.0 / ( dtemp + a2 * i1 * i1 ) );
		speq[2*i-nx] *= temp;
		speq[2*i-nx+1] *= temp;
	}
		
	for( i = 0; i < nxd2; i++ )
	{
		i1 = nxd2 - i;
		dtemp = alpha + a2 * i1 * i1; 
		
		for( j = nyd2; j < ny; j++ )
		{
			j1 = j - nyd2;
			temp = T( 1.0 / ( dtemp + b2 * j1 * j1 ) );
			F[i+nxd2][2*j-ny] *= temp; 
			F[i+nxd2][2*j-ny+1] *= temp;
		}
	}

	for( i = nxd2; i < nx; i++ )
	{
		i1 = i - nxd2;
		dtemp = alpha + a2 * i1 * i1; 
		
		for( j = nyd2; j < ny; j++ )
		{
			j1 = j - nyd2;
			temp = T( 1.0 / ( dtemp + b2 * j1 * j1 ) );
			F[i-nxd2][2*j-ny] *= temp; 
			F[i-nxd2][2*j-ny+1] *= temp; 
		}
	}

	// inverse FFT
	fft.Real2D( arr, speq, nx, ny, xar::OouraFft<T>::eDirInv );

	// trim back
	if( nx != nxF || ny != nyF )
	{
		xar::XArray2DMove<T> tmp( F );
		tmp.Trim( ( nx - nxF ) / 2, nx - nxF - ( nx - nxF ) / 2, ( ny - nyF ) / 2, ny - nyF - ( ny - nyF ) / 2 );
	}

	//thickness(r) = -(1/mu)*lnF
	//double amu = -1.0 / (4.0 * PI / wl * beta);
	//arr = F.GetArr();
	//try { for (i=0; i<F.GetNp(); i++) arr[i] = amu * log(arr[i]); }
	//catch (vo_math_exception& E)
	//{
	//	if (!strcmp(E.GetDescription(), "argument domain error")) throw vo_math_exception("negative absorption", "DP");
	//	else throw vo_math_exception(E.GetDescription(), "DP");
	//}
}

template <class T> void XA_2DTIE<T>::ReplaceModulus(xar::XArray2D< std::complex<T> >& F, xar::XArray2D<T>& A)
{
	// check input parameters
	if (!(F.GetDim1() == A.GetDim1() && F.GetDim2() == A.GetDim2()))
		throw std::invalid_argument("invalid_arguments 'F, A' in XA_2DTIE<T>::ReplaceModulus (arrays have different dimensions)");

	std::complex<T> *arrF = &F.front();
	T *arrA = &A.front();

	for (index_t i = 0; i < F.GetDim1() * F.GetDim2(); i++) arrF[i] *= (arrA[i] / abs(arrF[i]));
}

template <class T> void XA_2DTIE<T>::Homogenise(xar::XArray2D< std::complex<T> >& F, double delta2beta)
{
	T d2b = T(delta2beta), amp;
	std::complex<T>* arrF = &F.front();

	for (index_t i = 0; i < F.GetDim1() * F.GetDim2(); i++)
	{
		amp = abs(arrF[i]);
		arrF[i] = std::polar<T>(amp, d2b * log(amp));
	}
}

template <class T> void XA_2DTIE<T>::Homogenise1(xar::XArray2D< std::complex<T> >& F, double delta2beta, double defocus, char* pinput_mask)
{
	T b2d = T(1.0 / delta2beta);
	std::complex<T>* arrF = &F.front();

	xar::XArray2D<T> P(F.GetDim1(), F.GetDim2());
	//CArg(F, P);
	CPhaseUnwrap(F, P, pinput_mask); // this one works best in the context of GS-Hom
	//CPhaseTIE1(F, P, defocus, delta2beta);
	T* arrP = &P.front();

	for (index_t i = 0; i < F.GetDim1() * F.GetDim2(); i++)
		arrF[i] = std::polar<T>(exp(b2d * arrP[i]), arrP[i]);
}

template <class T> void XA_2DTIE<T>::EnforceSupport(xar::XArray2D< std::complex<T> >& F, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, std::complex<T> tMaskVal)
{
	XArray2DMove< std::complex<T> > xamove(F);
	xamove.Mask(iYLeft, iYRight, iXLeft, iXRight, tMaskVal);
}


#if (0)
// In order to use this function, the extern "C" module unwrap_2d_ljmu.c must be included into the project
// NOTE!!!: the description of input_mask in unwrap_2d_ljmu.c is incorrect: 0 corresponds to "good" points, 255 (unsigned char(-1)) corresponds to "bad" points
template<class T> void XA_2DTIE<T>::CPhaseUnwrap(xar::XArray2D< std::complex<T> >& F, xar::XArray2D<T>& Fi, char* pinput_mask)
{
	// prepare wrapped input phase array
	XArray2D<double> Fid(F.GetDim1(), F.GetDim2()), Fid1(F.GetDim1(), F.GetDim2());
	for (index_t i = 0; i < F.GetDim1(); i++)
		for (index_t j = 0; j < F.GetDim2(); j++)
			Fid[i][j] = (double)std::arg(F[i][j]);

	// prepare input mask
	unsigned char* pinput_mask1;
	XArray2D<char> input_mask;
	if (pinput_mask == 0)
	{
		input_mask.Resize(F.GetDim1(), F.GetDim2(), 0);
		pinput_mask1 = (unsigned char*)&input_mask.front();
	}
	else
		pinput_mask1 = (unsigned char*)pinput_mask;

	// do phase unwrapping
	unwrap2D(&Fid.front(), &Fid1.front(), pinput_mask1, (int)F.GetDim2(), (int)F.GetDim1(), 1, 1, 0, 0);

	// prepare output unwrapped phase array
	Fi.Resize(F.GetDim1(), F.GetDim2(), 0.0f);
	for (index_t i = 0; i < F.GetDim1(); i++)
		for (index_t j = 0; j < F.GetDim2(); j++)
			Fi[i][j] = (T)Fid1[i][j];
}



template <class T> void XA_2DTIE<T>::DPME(XArray2D<T>& I1, double d2b0, double d2b1, double R, short niter)
// does not give any interesting results!
{
	if (d2b0<=0 || d2b1<=0) throw std::invalid_argument("invalid_argument 'd2b0 or d2b1' in XA_2DTIE<T>::DPME (delta/beta must be positive)");
	if (niter<=0) throw std::invalid_argument("invalid_argument 'niter' in XA_2DTIE<T>::DPME (number of iterations must be positive)");

	vector<double> norms(niter+1);
	XArray2D<T> Temp;
	
	if (d2b1< d2b0) SWAP(d2b1, d2b0);
	double d2bst = (d2b1 - d2b0) / (niter - 1), d2b;

	for (long i=0; i<=niter; i++)
	{ 
		d2b = d2b0 + d2bst * i;
		Temp = I1;
		DP(Temp, d2b, R);
		norms[i] = Temp.Norm(7);
	}
	
	I1 = Temp;
}




//! Retrieves phase&amplitude from 1 image of a 'single-material' object using TIE approximation and PSF deconvolution
// Implements David Paganin's inversion algorithm applicable to in-line X-ray phase-contrast
// images of samples consisting of single component. Applicability condition is:
// beta(r)=const*delta(r), where n(r)=1-delta-i*beta is the complex refractive index.
// deltaoverbeta is delta/beta
// Simultaneously deconvolves with the image from Ker
// F - image; on exit it is replaced by the reconstructed INTENSITY (= exp(-4.0 * PI / wl * beta * thickness));
// Ker - convolution kernel
// R is a propagation distance (in units used in F.head - default are microns)
// alpha is the regularization parameter
//
// NOTE!!!: the program assumes that F is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iinc, 1 must not be subtracted!!!)
// NOTE: this program uses  the Ooura FFT library.
// NOTE: this progam is a modification of PoissonFFT() - see invLa.h
//
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
// @@@@@@@@@@@@@@@ NOTE!!!: consider renormalization compensating for regularization along the lines done in XArray2DFFT<T>::Deconvol
template <class T> void XA_2DTIE<T>::DPDeconv(xar::XArray2D<T>& F, xar::XArray2D<T>& Ker, double deltaoverbeta, double R, double alpha)
{
	if( &F==&Ker ) 
		throw std::invalid_argument("invalid_argument 'F and Ker' in XA_2DTIE<T>::DPDeconv (cannot deconvol XArray2D object with itself)");

	if( fabs(GetXStep(F) - GetXStep(Ker))>0.01*GetXStep(F) || fabs(GetYStep(F) - GetYStep(Ker))>0.01*GetYStep(F))
		throw std::invalid_argument("invalid_argument 'F and Ker' in XA_2DTIE<T>::DPDeconv (different steps)");

	if( deltaoverbeta <= 0 ) 
		throw std::invalid_argument("invalid_argument 'deltaoverbeta' in XA_2DTIE<T>::DPDeconv (delta/beta must be positive)");

	if( alpha < 0 ) 
		throw std::invalid_argument("invalid_argument 'alfa' in XA_2DTIE<T>::DPDeconv (alpha cannot be negative)");

	if( !GetIXAHWave2D( F ) ) 
		throw std::invalid_argument("invalid_argument 'F' in XA_2DTIE<T>::DPDeconv (wavehead is missing)");

	GetIXAHWave2D( F )->Validate();

	index_t nxF = F.GetDim1();
	index_t nyF = F.GetDim2();
	index_t nxKer = Ker.GetDim1();
	index_t nyKer = Ker.GetDim2();
	double wl( GetIXAHWave2D( F )->GetWl() ); 
	index_t nx = xar::max( nxF, nxKer );
	index_t ny = xar::max( nyF, nyKer );

	index_t i = 2;
	while (i<nx) i *= 2;
	nx = i; //nx is the smallest power of 2 not less than max(F.nx, Ker.nx)
	index_t j = 2;
	while (j<ny) j *= 2;
	ny = j; //ny is the smallest power of 2 not less than max(F.ny, Ker.ny)

	xar::XArray2DMove<T> tmp(F);
	tmp.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
	
	xar::XArray2DMove<T> tmp1(Ker);
	tmp1.Pad((nx-nxKer)/2, nx-nxKer-(nx-nxKer)/2, (ny-nyKer)/2, ny-nyKer-(ny-nyKer)/2, T(0));

	T* arr = &(F.front());
	std::vector<T> vecSpeq(nx*2);
	T* speq = &(vecSpeq[0]);
	T* arrK = &(Ker.front());
	std::vector<T> vecSpeqK(nx*2);
	T* speqK = &(vecSpeqK[0]);

	
	// forward FFTs
	xar::OouraFft<T> fft;

	fft.Real2D(arr, speq, nx, ny, xar::OouraFft<T>::eDirFwd);
	fft.Real2D(arrK, speqK, nx, ny, xar::OouraFft<T>::eDirFwd);

	//multiplication of the FFT by 1/((m^2+n^2)*Ker^(m,n))
	//the indexing follows that from the function FFTRe(...) - see vo_fftr2.h
	double factor = 2. / (nx * ny);
	double a( GetIXAHWave2D( F )->GetYhi() - GetIXAHWave2D( F )->GetYlo() + GetYStep( F ) );			// X/Y dimensions reversed
	double b( GetIXAHWave2D( F )->GetXhi() - GetIXAHWave2D( F )->GetXlo() + GetXStep( F ) );			// X/Y dimensions reversed
	double a2 = xar::PI * wl * R * deltaoverbeta / (a * a) / factor;
	double b2 = xar::PI * wl * R * deltaoverbeta / (b * b) / factor;
	double gamma = 1.0 / factor;
	double dtemp, hmnabs2;
	T temp0, temp1;
	index_t nxd2 = nx / 2, i1, nxi;
	index_t nyd2 = ny / 2, j1, nyj;

	//!!! renormalization of alpha
	double temp = Ker.Norm(xar::eNormL2);
	alpha *= (temp * temp) * gamma * (1.0 + xar::PI * wl * R * deltaoverbeta / (b * a));

	dtemp = gamma + b2 * nyd2 * nyd2;
	for (i=0; i<nxd2; i++) 
	{ 
		i1 = nxd2 - i;
		nxi = nx + 2 * i;
		hmnabs2 = speqK[nxi] * speqK[nxi] + speqK[nxi+1] * speqK[nxi+1];
		temp0 = T(1.0 / (alpha + hmnabs2 * (dtemp + a2 * i1 * i1)));
		temp1 = temp0 * (speq[nxi] * speqK[nxi] + speq[nxi+1] * speqK[nxi+1]);
		speq[nxi+1] = temp0 * (-speq[nxi] * speqK[nxi+1] + speq[nxi+1] * speqK[nxi]);
		speq[nxi] = temp1;
	}
	for (i=nxd2; i<nx; i++) 
	{ 
		i1 = i - nxd2;
		nxi = 2 * i - nx;
		hmnabs2 = speqK[nxi] * speqK[nxi] + speqK[nxi+1] * speqK[nxi+1];
		temp0 = T(1.0 / (alpha + hmnabs2 * (dtemp + a2 * i1 * i1)));
		temp1 = temp0 * (speq[nxi] * speqK[nxi] + speq[nxi+1] * speqK[nxi+1]);
		speq[nxi+1] = temp0 * (-speq[nxi] * speqK[nxi+1] + speq[nxi+1] * speqK[nxi]);
		speq[nxi] = temp1;
	}
		
	for (i=0; i<nxd2; i++)
	{
		i1 = nxd2 - i;
		dtemp = gamma + a2 * i1 * i1; 
		nxi = i + nxd2;
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2;
			nyj = 2 * j - ny;
			hmnabs2 = Ker[nxi][nyj] * Ker[nxi][nyj] + Ker[nxi][nyj+1] * Ker[nxi][nyj+1];
			temp0 = T(1.0 / (alpha + hmnabs2 * (dtemp + b2 * j1 * j1)));
			temp1 = temp0 * (F[nxi][nyj] * Ker[nxi][nyj] + F[nxi][nyj+1] * Ker[nxi][nyj+1]); 
			F[nxi][nyj+1] = temp0 * (-F[nxi][nyj] * Ker[nxi][nyj+1] + F[nxi][nyj+1] * Ker[nxi][nyj]); 
			F[nxi][nyj] = temp1;
		}
	}

	for (i=nxd2; i<nx; i++)
	{
		i1 = i - nxd2;
		dtemp = gamma + a2 * i1 * i1; 
		nxi = i - nxd2;
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2;
			nyj = 2 * j - ny;
			hmnabs2 = Ker[nxi][nyj] * Ker[nxi][nyj] + Ker[nxi][nyj+1] * Ker[nxi][nyj+1];
			temp0 = T(1.0 / (alpha + hmnabs2 * (dtemp + b2 * j1 * j1)));
			temp1 = temp0 * (F[nxi][nyj] * Ker[nxi][nyj] + F[nxi][nyj+1] * Ker[nxi][nyj+1]); 
			F[nxi][nyj+1] = temp0 * (-F[nxi][nyj] * Ker[nxi][nyj+1] + F[nxi][nyj+1] * Ker[nxi][nyj]); 
			F[nxi][nyj] = temp1;
		}
	}

	Ker.Truncate();

	//inverse FFT
	fft.Real2D(arr, speq, nx, ny, xar::OouraFft<T>::eDirInv);

	
	xar::XArray2DFFTRe<T> fftre(F);
	fftre.Shuffle();

	//trim back
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp(F);
		tmp.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
	}

	//thickness(r) = -(1/mu)*lnF
}



//! Retrieves phase&amplitude from 2 images at different defocus distances using TIE approximation
// TIE phase-amplitude retrieval from 2 images at different distances (in particular, I0 can be the object
//  plane intensity and I1 can be the in-line image at R2)
// I0 - image at a shorter distance R0; it is ASSUMED to be equal to the image(=object) plane intensity
// 	(note that for a non-absorbing object I0==incident_intensity)
// I1 - image at a longer distance R = R0+R2;  on exit it is replaced by the phase
// R - distance R-prime
// if UseFFT is true, then PoissonFFT is used for phase retrieval, else PoissonFMG is used
// alpha and ncycle - obvious meaning in the context of PoissonFFT or PoissonFMG
//
// NOTE!!!: energy conservation requires that I0.Norm(eNormAver)==I1.Norm(eNormAver);
// NOTE!!!: this is the 'exact TIE' solution algorithm for non-absorbing objects, but only an approximate
// one for objects with non-negligible absorption or non-uniform illumination;
//  it solves equation -Laplace(Fi)) = k*(I1/I0-1)/R'.
// NOTE!!!: it seems that more stable solution is obtained using -Laplace(Fi)) = k*(I1 -I0)/<I0>/R'
//
// NOTE: PoissonFMG requires square arrays with xst=yst, PoissonFFT does not require that
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::TIE2R(const xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, double R, bool UseFFT, long ncycle, double alpha)
{
	if (!GetIXAHWave2D(I0) || !GetIXAHWave2D(I1))
		throw std::invalid_argument("invalid_argument 'I0 or I1' in XA_2DTIE<T>::TIE2R (wavehead is missing)");
	
	double temp = I0.Norm(xar::eNormAver);
	double temp1 = I1.Norm(xar::eNormAver);
	if (fabs(temp - temp1) > 0.01 * temp)
		//throw std::invalid_argument("invalid_argument 'I0 and I1' in XA_2DTIE<T>::TIE2R (energy conservation does not hold)");
		I1 += T(temp - temp1);


	//I1 /= I0; //NOTE: experience seems to show that this division often leads to strong artefacts
	//I1 -= I1.Norm(eNormAver);
	I1 -= I0;
	I1 /= T(I0.Norm(xar::eNormAver));
	I1 *= T(xar::tPI / GetIXAHWave2D(I1)->GetWl() / R);

	XA_2DLapl<T> xaLapl(I1);
	if (UseFFT)
		xaLapl.PoissonFFT(alpha); 
	else 
		xaLapl.PoissonFMG(ncycle, alpha);
}



//! Retrieves phase&amplitude from 3 images at different defocus distances using TIE approximation
// TIE phase-amplitude retrieval from 3 images at different distances
// I0 - image at the middle distance R = R0;  it is the intenisty at R=R0
// I1 - image at the longest distance R = R0+R2;  on exit it is replaced by the phase at R=R0
// I2 - image at the shortest distance R0-R2; 
// R - distance R-prime (assumed to be the defocus distance between I0 and I1, and I2 and I0)
// if UseFFT is true, then PoissonFFT is used for phase retrieval, else PoissonFMG is used
// alpha and ncycle - obvious meaning in the context of PoissonFFT or PoissonFMG
//
// NOTE!!!: images I0, I1  and I2 HAVE TO BE BACKGROUND CORRECTED (i.e. DIVIDED by the 
//  corresponding INCIDENT intensity at wl (I/I_inc)), this algorithm is NOT self-normalizing, 
//  as it uses I0 and I2 - I1 separately.
// NOTE!!!: energy conservation requires that I0.Norm(eNormAver)==I1.Norm(eNormAver)==I2.Norm(eNormAver);
// NOTE!!!: this is the 'exact TIE' solution algorithm for absorbing objects even in the case of
//  non-uniform illumination; solves equation -div(I0*grad(Fi)) = k*(I1-I2)/R'.
//
// NOTE: PoissonFMG requires square arrays with xst=yst, PoissonFFT does not require that
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::TIE3R(xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, const xar::XArray2D<T>& I2, double R, bool UseFFT, long ncycle, double alpha)
{
	if (!GetIXAHWave2D(I0) || !GetIXAHWave2D(I1) || !GetIXAHWave2D(I2))
		throw std::invalid_argument("invalid_argument 'I0, I1 or I2' in XA_2DTIE<T>::TIE3R (wavehead is missing)");

	double temp = I0.Norm(xar::eNormAver);

	if (fabs(temp-I1.Norm(xar::eNormAver)) > 0.01 * temp) 
		throw std::invalid_argument("invalid_argument 'I0 and I1' in XA_2DTIE<T>::TIE3R (energy conservation does not hold)");

	if (fabs(temp-I2.Norm(xar::eNormAver)) > 0.01 * temp) 
		throw std::invalid_argument("invalid_argument 'I0 and I2' in XA_2DTIE<T>::TIE3R (energy conservation does not hold)");

	I1 -= I2; // !!! I2 is assumed to be the image at the shortest distance
	I1 *= T(xar::tPI / GetIXAHWave2D(I1)->GetWl() / (2.0 * R));
	xar::XArray2D<T> Itemp(I0); //stores I0, as it is spoiled by MDivIGradFFT()

	XA_2DLapl<T> xaLapl(I1);
	if (UseFFT) 
		xaLapl.MDivIGradFFT(I0, alpha);
	else 
		xaLapl.MDivIGradFMG(I0, ncycle, alpha);
	I0 = Itemp;
}



//! Retrieves phase&amplitude from 2 images at different defocus distances using modified TIE approximation
// Modified TIE phase-amplitude retrieval from 2 images at different distances
// I0 = object plane intensity (replaced by the "reconstructed" intensity" on exit)
// I1 = image plane intensity (replaced by the reconstructed phase on exit)
// R = object-to-image distance
// gamma2 = maximum estimated delta/beta ratio for the object
// ncycle = maximum allowed number of iterations
// sigma = L2 noise level in the input data (iterations are interrupted, when the current increment becomes smaller than sigma)
//
// NOTE: input intensities must be flat-field corrected
// !!!NOTE: the code of this function is DIFFERENT from that of the function ModTIE2R in vo_TIE.h
// This function (a) does not assume weak absorption in the object; (b) always uses gamma1 = 0 (for simplifying the equations)
// It solves the following equations:
// I1 / Iin = Q1 * Q2 - R * gamma2 / (2k) * div(Q1 * grad Q2),
// I0 / Iin = Q1 * Q2,
// by iterating as follows
// Q2(m+1) = [1 -  sigma2 * grad^2]^(-1){I1 / I0 * Q2(m)  - sigma2 * grad^2(Q2) + sigma2 / Q1(m) * div[Q1(m) * grad Q2(m)]},
// Q1(m+1) = I0  / [Iin * Q2(m)],
// where sigma2 = R * gamma2 / (2k).
// and the phase is retrieved at the end as follows:
// Phi0 = (gamma2 / 2) * log(Q2).
template <class T> void XA_2DTIE<T>::ModTIE2R(xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, double R, double gamma2, long ncycle, double sigma)
{
	#define MODTIE2R_TMP_CODE 0 // enables temporary code if the value is not zero

	if (gamma2 <= 0)
		throw std::invalid_argument("invalid_argument 'gamma2' in XA_2DTIE<T>::ModTIE2R (0 = gamma1 < gamma2 condition is violated)");

	if (ncycle <= 0)
		throw std::invalid_argument("invalid_argument 'ncycle' in XA_2DTIE<T>::ModTIE2R (max number of iterations is not positive)");

	const IXAHWave2D* ph = GetIXAHWave2D(I0);
	if (!ph)
		throw std::invalid_argument("invalid_argument 'I0' in XA_2DTIE<T>::ModTIE2R (wavehead is missing)");
	if (!(ph->Equivalent(GetIXAHWave2D(I1))))
		throw std::invalid_argument("invalid_argument 'I0 or I1' in XA_2DTIE<T>::ModTIE2R (different Waveheads)");

	//if (sigma <= 0)
		//throw std::invalid_argument("invalid_argument 'sigma' in XA_2DTIE<T>::ModTIE2R (noise level in input data is negative)");

	double dif0 = I0.Norm(xar::eNormL2);
	double dif1 = 0.99999 * dif0;
	if (fabs(dif0-I1.Norm(xar::eNormL2)) > 0.01 * dif0)
		throw std::invalid_argument("invalid_argument 'I0 or I1' in XA_2DTIE<T>::ModTIE2R (energy conservation does not hold)");

	T sigma2 = T(R * gamma2 * ph->GetWl() / (4.0 * xar::PI)); // R * gamma2 / (2k)

	if (I0.Norm(xar::eNormAver) <= 0)
		throw std::invalid_argument("invalid_argument 'I0' in XA_2DTIE<T>::ModTIE2R (average value of intput intensity is not positive)");
	if (I1.Norm(xar::eNormAver) <= 0)
		throw std::invalid_argument("invalid_argument 'I1' in XA_2DTIE<T>::ModTIE2R (average value of intput intensity is not positive)");

	double alpha = sigma * I0.Norm(xar::eNormL2); // relative noise level --> absolute noise level
	
	xar::XArray2D<T> Q1(I0), Q2, Q2a, Q2b, Q2prev;
	xar::XArray2DFilt<T> Q2aFilt(Q2a), Q2bFilt(Q2b);

#if(MODTIE2R_TMP_CODE)
	// temporary code
	xar::XArray2D<T> Q1a(Q1), Ph, Ph0;
	xarfileio::XArFileIO_VOD::ReadFile(Ph0, L"D:\\Method1_10m\\Fi0_aver0.vod");
	xar::XArray2DMove<T> Ph0move(Ph0);
	Ph0move.Transpose();
	double denom = Ph0.Norm(xar::eNormL2);
	xar::XArray1D<double> arrDif(ncycle+1), arrDif1(ncycle+1);
#endif

	// initialization
	long n(-1);
	//Q1.Fill(T(1));
	Q2 = I1; // in the future, different initialization can be tried
	Q2prev = I0;
	I1 /= I0;
	
	// iteration cycle
	while (++n < ncycle && dif1 < dif0 && dif1 > alpha)
	{
		// new iteration cycle
		if (n > 0)
		{
			dif0 = dif1;
			Q2a = Q2;
			Q2aFilt.MLaplacian();
			Q2a *= sigma2;
			Q2b = Q2;
			Q2bFilt.MDivIGrad(Q1);
			Q2b /= Q1;
			Q2b *= sigma2;
			Q2 *= I1;
			Q2 += Q2a;
			Q2 -= Q2b;
		}
		DP(Q2, gamma2, R); // new Q2 obtained
		Q1 = I0;
		Q1 /= Q2; // new Q1 obtained
		//check for convergence
		if (n > 0) { Q2prev -= Q2; dif1 = Q2prev.Norm(xar::eNormL2); }
		Q2prev = Q2; // store for the next cycle

#if(MODTIE2R_TMP_CODE)
		// temporary code
		// calculate and save the reconstruction error
		Q2a = Q2;
		Q2a.Log();
		Q2a *= T(0.5 * gamma2);
		Ph = Q2a; // output phase
		Ph -= T(Ph.Norm(xar::eNormAver));
		Ph -= Ph0;
		arrDif[n] = Ph.Norm(xar::eNormL2) / denom;
		arrDif1[n] = dif1;
#endif

	}

	Q2.Log();
	Q2 *= T(0.5 * gamma2);
	I1 = Q2; // output phase

#if(MODTIE2R_TMP_CODE)
	// temporary code
	xarfileio::XArFileIO_VOD::WriteFile(arrDif, L"D:\\Method1_10m\\AAArecerror.vod", xarfileio::WriteOptionsVOD_ASCII_NoOptimize);
	xarfileio::XArFileIO_VOD::WriteFile(arrDif1, L"D:\\Method1_10m\\AABrecerror.vod", xarfileio::WriteOptionsVOD_ASCII_NoOptimize);
#endif

}



//! Retrieves phase&amplitude from 3 images at different defocus distances using modified TIE approximation
// Modified TIE phase-amplitude retrieval from 3 images at different distances
// I0 - image at the middle distance R = R0;  it is the intenisty at R=R0
// I1 - image at the longest distance R = R0+R2;  on exit it is replaced by the phase at R=R0
// I2 - image at the shortest distance R0-R2; 
// R - distance R-prime (assumed to be the distance between I0 and I1, and I2 and I0)
// ncycle - obvious meaning in the context of SchrodFMG
// averV - regularization parameter expressed in percentage, see below
//
// NOTE!!!: images I0, I1  and I2 HAVE TO BE BACKGROUND CORRECTED (i.e. DIVIDED by the 
//  corresponding INCIDENT intensity at wl (I/I_inc)), this algorithm is NOT self-normalizing, 
//  as it uses I0 and I2 - I1 separately.
// NOTE!!!: energy conservation requires that I0.Norm(eNormAver)==I1.Norm(eNormAver)==I2.Norm(eNormAver);
// NOTE!!!: this is the 'exact ModTIE' solution algorithm for absorbing objects even in the case of
// non-uniform illumination; solves equation -Laplace(Psi) + D(I)Psi = 2*k*(I1^0.5-I2^0.5)/R'.
//
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::ModTIE3R(const xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, const xar::XArray2D<T>& I2, double R,  short ncycle, double averV)
{
	if (!GetIXAHWave2D(I0) || !GetIXAHWave2D(I1) || !GetIXAHWave2D(I2))
		throw std::invalid_argument("invalid_argument 'I0, I1 or I2' in XA_2DTIE<T>::ModTIE3R (wavehead is missing)");

	double temp = I0.Norm(xar::eNormAver);

	if (fabs(temp-I1.Norm(xar::eNormAver)) > 0.01 * temp) 
		throw std::invalid_argument("invalid_argument 'I0 and I1' in XA_2DTIE<T>::ModTIE3R (energy conservation does not hold)");

	if (fabs(temp-I2.Norm(xar::eNormAver)) > 0.01 * temp) 
		throw std::invalid_argument("invalid_argument 'I0 and I2' in XA_2DTIE<T>::ModTIE3R (energy conservation does not hold)");

	xar::XArray2D<T> Itemp(I2);
	I1 ^= 0.5;
	Itemp ^= 0.5;
	I1 -= Itemp;
	I1 *= T(xar::tPI / GetIXAHWave2D(I1)->GetWl() / R);
	if (averV < 100) // calculate true V=D2(I)
	{
		Itemp = I0;
		Itemp ^= 0.5;
		xar::XArray2DFilt<T> filt2D(Itemp);
		filt2D.MLaplacian(true);
		Itemp *= -1;
		if (averV > 0) // smooth true V for regularization
		{
			xar::XArray2DFFTRe<T> fftre(Itemp);
			fftre.FilterGauss((Itemp.GetDim1() - 1) * averV / 200.0, (Itemp.GetDim2() - 1) * averV / 200.0);
		}
	}
	else // replace V by a constant (averV-100) (in particular, averV==100 => geom.optics)
	{ 
		Itemp *= 0; 
		Itemp += T(averV - 100); 
	}
	XA_2DLapl<T> xaLapl(I1);
	xaLapl.SchrodFMG(Itemp, ncycle);
	Itemp = I0;
	Itemp ^= 0.5;
	I1 /= Itemp;
}



//! Retrieves phase&amplitude from 2 images at different wavelengths (energies) using TIE approximation
// TIE phase-amplitude retrieval from 2 images at different wavelengths
// (neglecting the term grad(fi)*grad(I))
// I0 and I1 - images (on exit they are replaced by the intensity and phase, respectively, at wl1)
// R = distance R-prime
// if UseFFT is true, then PoissonFFT is used for phase retrieval, else PoissonFMG is used
// alpha and ncycle - obvious meaning in the context of PoissonFFT or PoissonFMG
//
// NOTE!!!: the programs assumes that I0 and I1 are background-corrected, i.e. 
// DIVIDED by the corresponding INCIDENT intensities at wl[k] (Ik/Ik_inc, 1 must not be subtracted!!!)
//
// NOTE!!!: PoissonFMG requires square arrays with xst=yst, PoissonFFT does not require that
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::TIE2E(xar::XArray2D<T>& I0, xar::XArray2D<T>& I1, double R, bool UseFFT, long ncycle, double alpha)
{
	double wl[2];
	IXAHWave2D* ph0 = GetIXAHWave2D(I0);
	IXAHWave2D* ph1 = GetIXAHWave2D(I1);
	if (!ph0 || !ph1)
		throw std::invalid_argument("invalid_argument 'I0 or I1' in XA_2DTIE<T>::TIE2E (wavehead missing)");
	wl[0] = ph0->GetWl(); 
	wl[1] = ph1->GetWl();

	if (wl[0]==wl[1]) 
		throw std::invalid_argument("invalid_argument 'I0 and I1' in XA_2DTIE<T>::TIE2E (identical wavelengths)");

	double sigma = wl[1] / wl[0];
	double asigma2 = -1.0 / (sigma * sigma);

	xar::XArray2D<T> gtemp;

	I0.Log(); I1.Log();

	//reconstruct contact intensity
	gtemp = I1;
	gtemp *= T(asigma2);
	IXAHead* htemp1 = gtemp.GetHeadPtr()->Clone();
	gtemp.SetHeadPtr(GetIXAHWave2D(I0)->Clone());
	gtemp += I0;
	gtemp.SetHeadPtr(htemp1);
	gtemp *= -T(wl[0] / (wl[1] - wl[0]));
	gtemp.Exp();
	
	//reconstruct the phase
	I1 *= T(asigma2 / sigma);

	htemp1 = I1.GetHeadPtr()->Clone();
	I1.SetHeadPtr(GetIXAHWave2D(I0)->Clone());
	I1 += I0; 
	I1.SetHeadPtr(htemp1);
	I1 *= T(xar::tPI / R / (wl[1] - wl[0]) * sigma); //-Laplacian(phase[0])
	I1 -= T(I1.Norm(xar::eNormAver)); //!!! normalize so that Ph0 does not contain Z4
	
	XA_2DLapl<T> xaLapl(I1);
	if (UseFFT) 
		xaLapl.PoissonFFT(alpha); 
	else 
		xaLapl.PoissonFMG(ncycle, alpha);

	I0 = gtemp;
	I1.SetHeadPtr(GetIXAHWave2D(I0)->Clone()); //to change wl
}



//! Retrieves phase&amplitude from 3 images at different wavelengths (energies) using TIE approximation
// TIE phase-amplitude retrieval from 3 images at different wavelengths.
// I0, I1 and I2 - images (0=middle, 1=low_energy, 2=high_energy)
// (on exit I0 is replaced by the intensity and I1 by the phase at wl0)
// R = distance R-prime
// if UseFFT is true, then PoissonFFT is used for phase retrieval, else PoissonFMG is used
// alpha and ncycle - obvious meaning in the context of PoissonFFT or PoissonFMG
//
// NOTE!!!: the programs assumes that I0, I1 and I2 are all background-corrected, i.e.
// DIVIDED by the corresponding INCIDENT intensities at wl[k] (Ik/Ik_inc, 1 must not be subtracted!!!)
//
// NOTE!!!: PoissonFMG requires square arrays with xst=yst, PoissonFFT does not require that
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::TIE3E(
	xar::XArray2D<T>& I0, 
	xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& I2, 
	double R, 
	bool UseFFT, 
	long ncycle, 
	double alpha)
{
	// calculate inversion matrix A^(-1)
	IXAHWave2D* ph0 = GetIXAHWave2D(I0);
	IXAHWave2D* ph1 = GetIXAHWave2D(I1);
	IXAHWave2D* ph2 = GetIXAHWave2D(I2);
	
	if (!ph0 || !ph1 || !ph2)
		throw std::invalid_argument("invalid_argument 'I0, I1 or I2' in XA_2DTIE<T>::TIE3E (wavehead missing)");

	xar::XArray2D<double> A(3,3); //coefficients of the linear system of TIEs for wl0, wl1 and wl2
	double sig1, sig2; //sigma[i]=wl[i]/wl[0]
	double wl[3], g[3]; //g[i]=R*wl[i]/(2*PI)
	wl[0] = ph0->GetWl();
	wl[1] = ph1->GetWl();
	wl[2] = ph2->GetWl();
	
	if (wl[0]==wl[1] || wl[0]==wl[2] || wl[1]==wl[2])  
		throw std::invalid_argument("invalid_argument 'I0, I1 and I2' in XA_2DTIE<T>::TIE3E (some identical wavelengths)");
	
	sig1 = wl[1] / wl[0]; sig2 = wl[2] / wl[0];
	for (index_t i=0; i<3; i++) g[i] = R * wl[i] / xar::tPI;
	A[0][0] = -1.0; A[0][1] = g[0]; A[0][2] = g[0];
	A[1][0] = -pow(sig1,3); A[1][1] = g[1]*sig1; A[1][2] = g[1]*pow(sig1,4);
	A[2][0] = -pow(sig2,3); A[2][1] = g[2]*sig2; A[2][2] = g[2]*pow(sig2,4);	
	xar::XArray2DLin<double> lin2(A);
	lin2.Invert();

	// reconstruct contact intensity
	xar::XArray2D<T> gtemp, htemp;
	I0.Log(); I1.Log(); I2.Log();
	htemp = I0; htemp *= T(A[0][0]); gtemp = htemp; 
	IXAHead* htemp1 = gtemp.GetHeadPtr()->Clone();
	gtemp.SetHeadPtr(GetIXAHWave2D(I1)->Clone());
	htemp = I1; htemp *= T(A[0][1]); gtemp += htemp; 
	gtemp.SetHeadPtr(GetIXAHWave2D(I2)->Clone());
	htemp = I2; htemp *= T(A[0][2]); gtemp += htemp; 
	gtemp.SetHeadPtr(htemp1);
	gtemp *= -1;
	// sig1 = gtemp.Norm(-1); //maximum
	// if (sig1>log(1.0)) gtemp += (log(1.0)-sig1); //!!! normalize so that I0<=1.0
	gtemp.Exp(); // output intensity if activated, otherwise - attenuation

	// reconstruct the phase
	I0 *= T(A[1][0]); I1 *= T(A[1][1]); I2 *= T(A[1][2]);
	// we need Array2D<T> references to add XArray2D<T>s with different 'wl's
	htemp1 = I1.GetHeadPtr()->Clone();
	I1.SetHeadPtr(GetIXAHWave2D(I0)->Clone());
	I1 += I0; 
	I1.SetHeadPtr(GetIXAHWave2D(I2)->Clone());
	I1 += I2; // -Laplacian(phase[0])
	I1.SetHeadPtr(htemp1);
	// I1 -= I1.Norm(eNormAver); //!!! normalize so that Ph0 does not contain Z4 (optional)

	XA_2DLapl<T> xaLapl(I1);
	if (UseFFT) 
		xaLapl.PoissonFFT(alpha); 
	else 
		xaLapl.PoissonFMG(ncycle, alpha);
	
	I1.SetHeadPtr(GetIXAHWave2D(I0)->Clone()); // to change wl
	I0 = gtemp;
}
 


//! Checks the validity of the TIE approximation for given image(s)
// I0, I1, I2 = input images (I0 and I2 can be empty)
// returns true if the TIE validity condition is satisfied, false otherwise
//
// validity condition is: max_over_x,y |F(x,y) - 1| < 0.5
// F = (I1 - <I1>) / <I1> for TIE1 and DP; F = (I1 - I0) / I0 for TIE2R; F = (I1 - I2) / I0 for TIE3R
// F = (I0^sigma1^3 - I1) / I1 for TIE2E; F = (I2^sigma12^3 - I1) / I1 for TIE3E
template <class T> bool XA_2DTIE<T>::TIE_Validity(const xar::XArray2D<T>& I0, const xar::XArray2D<T>& I1, const xar::XArray2D<T>& I2)
{
	xar::XArray2D<T> F;

	if (I1.size()==0) 
		throw std::invalid_argument("invalid_argument 'I1' in XA_2DTIE<T>::TIE_Validity (first image not supplied)");

	if (I0.size()==0 && I2.size()==0) //TIE1 and DP case
	{
		F = I1; 
		F /= T(I1.Norm(xar::eNormAver));
		F -= T(1.0);
	}
	else
	{
		if (I0.size()==0 && I2.size()!=0) 
			throw std::invalid_argument("invalid_argument 'I0 and I2' in XA_2DTIE<T>::TIE_Validity (second image given without the zeroth)");

		if (I0.size()!=0 && I2.size()==0)
		{
			if (I0.SameHead(I1)) //TIE2R
			{
				F = I1; 
				F /= I0;
				F -= T(1.0);
			}
			else //TIE2E
			{
				const IXAHWave2D* ph1 = GetIXAHWave2D(I1);
				const IXAHWave2D* ph0 = GetIXAHWave2D(I0);

				if (!ph0 || !ph1)
					throw std::invalid_argument("invalid_argument 'I0 or I1' in XA_2DTIE<T>::TIE_Validity (wavehead missing)");

				double sigma1 = ph1->GetWl() / ph0->GetWl();
				
				if (sigma1 == 1) 
					throw std::invalid_argument("invalid_argument 'I0 and I1' in XA_2DTIE<T>::TIE_Validity (images have same wl, but different headers)");
				
				F = I0;
				F.SetHeadPtr(ph1->Clone());
				F ^= T(sigma1 * sigma1 * sigma1);
				F /= I1;
				F -= T(1.0);
			}
		}
		else
		{
			if (I0.SameHead(I1) && I0.SameHead(I2)) //TIE3R
			{
				F = I1;
				F -= I2;
				F /= I0;
			}
			else //TIE3E
			{
				const IXAHWave2D* ph1 = GetIXAHWave2D(I1);
				const IXAHWave2D* ph2 = GetIXAHWave2D(I2);
				
				if (!ph1 || !ph2)
					throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DTIE<T>::TIE_Validity (wavehead missing)");

				double sigma12 = ph1->GetWl() / ph2->GetWl();
				
				if (sigma12 == 1) 
					throw std::invalid_argument("invalid_argument 'I2 and I1' in XA_2DTIE<T>::TIE_Validity (images have same wl, but different headers)");
				
				F = I2;
				F.SetHeadPtr(ph1->Clone());
				F ^= T(sigma12 * sigma12 * sigma12);
				F /= I1;
				F -= T(1.0);
			}
		}
	}

	return F.Norm(xar::eNormC0) < 0.5;
}



//! Solves 2D phase unwrapping problem using TIE approximation
// deterministic extraction of a continuous phase from a 2D complex amplitude distribution
// by means of solution of equation -div[I*grad(Fi)] = -I * Im[-Laplacian(U)/U]
// U = input complex amplitude
// Fi = output phase distribution
// ncycle = number of V-cycles in the FMG solution of equation -div[I*grad(Fi)] = F
//
// NOTE!!!: in most cases it is sufficient to use the simpler, faster and more accurate function CArg
// instead of this function
//
// NOTE!!!: a better treatment of phase boundary conditions should be considered here by means of 
// 1D continuous phase extraction along the 1D boundary, supply of the boundary values to 
// a MDivIGrad function modified to use explicit Dirichlet boundary conditions. Note that a discontinuity
// (>2*PI) of phase values along a closed 1D boundary loop would indicate the presence of a vortex!
template <class T> void XA_2DTIE<T>::CPhaseTIE(const xar::XArray2D<std::complex<T> >& U, xar::XArray2D<T>& Fi, long ncycle, double alpha)
{
	xar::XArray2D<T> I;
	Abs2(U, I);
	xar::XArray2D<std::complex<T> > mLU(U);
	xar::XArray2DFilt<std::complex<T> > filt2D(mLU);
	filt2D.MLaplacian(true);
	Im(mLU, Fi);
	Fi *= I;
	XA_2DLapl<T> xaLapl(Fi);
	xaLapl.MDivIGradFMG(I, ncycle, alpha);
}
#endif

//! Solves 2D phase unwrapping problem using TIE-Hom approximation
// deterministic extraction of a continuous phase from a 2D complex amplitude distribution
//
// NOTE!!: this function changes the input complex array U
template <class T> void XA_2DTIE<T>::CPhaseTIE1(xar::XArray2D<std::complex<T> >& U, xar::XArray2D<T>& Fi, double defocus, double delta2beta)
{
	// make the modulus of U equal to 1
	for (index_t i = 0; i < U.GetDim1(); i++)
		for (index_t j = 0; j < U.GetDim2(); j++)
			U[i][j] /= abs(U[i][j]);

	// calculate forward propagated intensity
	XArray2DFFT<float> xafft2(U);
	xafft2.Fresnel(defocus*10);
	Abs2(U, Fi);

	// extract phase using TIE-Hom
	Fi -= T(1.0);
	DP(Fi, delta2beta*100, defocus*10);
	Fi *= T(delta2beta / 2.0);
}

#if (0)
//! Qualitatively enhances phase-contrast images using unsharp masking and optional 'refocusing'
// --------------------------------------------------------------------------------------------------------------------
// DESCRIPTION OF ALGORITHM :
// The algorithm implemented here falls into 5 consecutive steps :
// (0)  The object-plane intensity ('refocused' image) may be retrieved using the DP algorithm. The resultant image may be used for
//		(or instead of) the blurred image. By default, the 'refocusing' is not done, i.e. raw image is used
// (1)	Blur the 'refocused' image by convolution with a Gaussian kernel of sufficient width for all phase-contrast effects to 
//		be destroyed.
// (2)	Normalize the raw and blurred images so that they have the same range of greylevel values.  The blurred image is 
//		then multiplied by (weighting_factor - 1), where weighting_factor lies between 0 and 1 inclusive.  If weighting_factor=0 
//		then the function will return a phase-enhanced image which displays only edges and eliminates all plateau in the original image.  
//		If weighting_factor=1 then the function will return the raw image.  Intermediate values of weighting_factor interpolate between 
//		these extremes, allowing one to achieve a pleasing amalgam of the raw and edge-enhanced data.
// (3)	Calculate the sum of (i) the raw image and (ii) the reweighted blurred image to yield the phase-enhanced image.
// (4)	Now we need to "crunch" the image.  This is done by firstly calculating the histogram of greylevels in the image, and then
//		using this to calculate the "upper" and "lower" greyelevels, which are defined such that a certain fraction -- namely the 
//		"histogram_factor" parameter -- of greylevels lie between these extremes.  Greylevels which lie between these extremes
//		are left unchanged ; all other greylevels are altered to take the nearest extreme value.
// --------------------------------------------------------------------------------------------------------------------
// ADDITIONAL NOTES :
// (1)	If the image dimensions are not powers of 2, the image will be padded and then trimmed back .
// --------------------------------------------------------------------------------------------------------------------
template <class T> void XA_2DTIE<T>::PhaseEnhancement(
	xar::XArray2D<T>& F, 
	double deltaoverbeta, 
	double R, 
	long blur_radius_in_pixels, 
	double weighting_factor, 
	double histogram_factor)
{
// (0)	Check validity of input parameters
	// check input parameters
	IXAHWave2D* ph2 = GetIXAHWave2D(F);
	
	if (!ph2) 
		throw std::invalid_argument("invalid_argument 'F' in XA_2DTIE<T>::PhaseEnhancement (wavehead is missing)");
	ph2->Validate();
	
	if (weighting_factor < 0 || weighting_factor > 1)
		throw std::invalid_argument("invalid_argument 'weighting_factor' in XA_2DTIE<T>::PhaseEnhancement (weighting factor must be between 0 and 1 inclusive)");
	
	if (histogram_factor < 0 || histogram_factor > 1)
		throw std::invalid_argument("invalid_argument 'histogram_factor' in XA_2DTIE<T>::PhaseEnhancement (histogram_factor must be between 0 and 1 inclusive)");

// (0) Refocus
	xar::XArray2D<T> F_refocused(F); // 'refocused' image
	
	if (deltaoverbeta > 0) 
		DP(F_refocused, deltaoverbeta, R);

// (1) Blur 'refocused' image by convolution with Gaussian kernel of sufficient width for all phase-contrast effects to be destroyed.
//	   Gaussian kernel is xst*yst*1.0/[2*PI*sigmaX*sigmaY]*exp{-0.5*(x^2/sigmaX^2+y^2/sigmaY^2))}, where xst, yst, sigmaX, sigmaY 
//	   are respectively the four numbers stuffed in the "aparam" array below (description lifted from Tim Gureyev's routine "Gauss2D"). 
//	   NOTE!!!: the normalization (*xst*yst) is chosen to preserve the average after convolution.
	if (blur_radius_in_pixels < 0) // default blur radius = 5% of image size
		blur_radius_in_pixels = long((F.GetDim1() + F.GetDim2()) / 40.0);
	
	xar::XArray2DFFTRe<T> IFilter(F_refocused);
	IFilter.FilterGauss(blur_radius_in_pixels, blur_radius_in_pixels);

	// (2)	Normalize refocused and raw images so that the ratio of greylevel excursions is equal to the weighting factor - 1.
	F -= T(F.Norm(xar::eNormMin)); // remove the D.C from the image
	F_refocused -= T(F_refocused.Norm(xar::eNormMin)); // remove the D.C from the image
	double normalisation_factor = F_refocused.Norm(xar::eNormMax);
	
	if (normalisation_factor != 0) 
		normalisation_factor = F.Norm(xar::eNormMax)/normalisation_factor;

	F_refocused *= T(normalisation_factor); // raw and blurred images now have the same range of grey levels
	F_refocused *= T(weighting_factor - 1); 

// (3)	Calculate the phase-enhanced image.
	F += F_refocused;

//	(4)	Crunch the image histogram.  
	const index_t number_of_bins_for_histogram = 1000; // 1000 bins is plenty 
	std::vector<index_t> histogram(number_of_bins_for_histogram + 1, 0); // initialise histogram to zero
	index_t bin_number_for_intensity_of_this_pixel;
	F -= T(F.Norm(xar::eNormMin)); // make zero correspond to black in the image as it simplifies some of the expressions below
	T umax = T(F.Norm(xar::eNormMax));
	
	if (umax == 0) 
		return; // flat (featureless) image

	index_t i, j;
	double a = (double)number_of_bins_for_histogram / umax; // constant used in the loop below 
	
	for (i = 0; i < F.GetDim1(); i++)
		for (j = 0; j < F.GetDim2(); j++)
		{	
			bin_number_for_intensity_of_this_pixel = (index_t)( F[i][j] * a ); // work out which bin a given pixel belongs in ...
			histogram[bin_number_for_intensity_of_this_pixel]++; // ... and throw it into that bin	
		}
	
	// integrate the histogram
	for (i = 1; i <= number_of_bins_for_histogram; i++)
		histogram[i] += histogram[i - 1];
	
	// work out the truncation levels
	index_t lower_level = 0, upper_level = 0;
	index_t max_value_of_histogram = histogram[number_of_bins_for_histogram];
	
	for (i = 0; i <= number_of_bins_for_histogram; i++)
		if( histogram[i] >= (max_value_of_histogram * ((1.0 - histogram_factor) / 2.0)) )
		{
			lower_level = i;
			break;
		}
	
	for (i = number_of_bins_for_histogram; i > 0; i--)
		if( histogram[i] <= (max_value_of_histogram * ((1.0 + histogram_factor) / 2.0)) )
		{
			upper_level = i;
			break;
		}
	
	T max_greyscale_level_after_crunching = umax * ((T)upper_level) / ((T)number_of_bins_for_histogram);
	T min_greyscale_level_after_crunching = umax * ((T)lower_level) / ((T)number_of_bins_for_histogram);
	
	// here we actually do the image crunching
	for (i = 0; i < F.GetDim1(); i++)
		for (j = 0; j < F.GetDim2(); j++)
			if( F[i][j] < min_greyscale_level_after_crunching)
				F[i][j] = min_greyscale_level_after_crunching;
			else
				if( F[i][j] > max_greyscale_level_after_crunching)
					F[i][j] = max_greyscale_level_after_crunching;
}

//! Deblurrs an image of a single-component sample by means of defocus
//	Deblurrs the image by using the DP() phase retrieval and 2nd order Taylor deconvolution.
// deltaoverbeta is delta/beta
// R is a propagation distance (in units used in F.head - default are microns)
// sigma is a st.deviation of the Gaussian PSF
// On exit F is replaced by the deblurred object-plane INTENSITY
//
// NOTE!!!: the program assumes that F is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iinc, 1 must not be subtracted!!!)
// NOTE: this progam is a modification of PoissonFFT() - see invLa.h
//
// NOTE: if the image dimensions are not powers of 2, the image will be padded and then trimmed back.
template <class T> void XA_2DTIE<T>::Deblur(xar::XArray2D<T>& F, double deltaoverbeta, double R, double sigma)
{
	if (deltaoverbeta <= 0) 
		throw std::invalid_argument("invalid_argument 'deltaoverbeta' in XA_2DTIE<T>::Deblur (delta/beta must be positive)");

	if (R == 0) 
		throw std::invalid_argument("invalid_argument 'R' in XA_2DTIE<T>::Deblur (R cannot be zero)");

	if (sigma < 0) 
		throw std::invalid_argument("invalid_argument 'sigma' in XA_2DTIE<T>::Deblur (sigma cannot be negative)");

	IXAHWave2D* ph2 = GetIXAHWave2D(F);
	
	if (!ph2) 
		throw std::invalid_argument("invalid_argument 'F' in XA_2DTIE<T>::Deblur (wavehead is missing)");
	
	ph2->Validate();
	double wl = ph2->GetWl(); 
	double b2overalpha = -sigma * sigma * xar::tPI / deltaoverbeta / R / wl;

	xar::XArray2D<T> I0(F);
	DP(I0, deltaoverbeta, R);

	I0 *= T(1.0 + b2overalpha);
	F *= T(-b2overalpha);
	F += I0;
}
#endif

//! Returns the _eValueType corresponding to T=float
template<> inline xar::_eValueType XA_2DTIE<float>::GetValuetype() { return xar::eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> inline xar::_eValueType XA_2DTIE<double>::GetValuetype() { return xar::eXADouble; }

} // namespace xar closed

//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//

// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XA_2DTIE<float>;
	template class xar::XA_2DTIE<double>;
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
