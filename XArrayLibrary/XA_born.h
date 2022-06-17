//Header XA_Born.h
//
//
//	HEADER FILE TITLE:
//
//		Phase retrieval algorithms based on the 1st Born and Rytov approximations
//
/*!
	\file		XA_Born.h
	\brief		Phase retrieval algorithms based on the 1st Born and Rytov approximations
	\par		Description:
		This header contains a class that provides phase retrieval services 
		based on the 1st Born and Rytov approximations for XArray2D<T> objects
*/

#pragma once

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_move2.h"
#include "XA_fft2.h"
//#include "XA_Optimisation.h"

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
//Class XA_2DBorn<T>
//
//	Phase retrieval algorithms based on the 1st Born and Rytov approximations
//
/*!
	\brief		Phase retrieval algorithms based on the 1st Born and Rytov approximations
	\par		Description:
				This class template provides 1st Born and Rytov based phase retrieval
				services for XArray2D<T> objects
	\remarks	An object of this class represents an interface exposing several functions
				that provide various 1st Born and Rytov based phase retrieval services				
	\remarks	This class can only be created for T=float or T=double
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XA_2DBorn
	{
	// Typedefs
	public:
		typedef T type;
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XA_2DBorn() {}
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XA_2DBorn(const XA_2DBorn<T>& rCopy) : m_vTemp(rCopy.m_vTemp) { GetValuetype(); }
	public:
		//! Destructor
		~XA_2DBorn() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XA_2DBorn<T>& rCopy) {	if (this == &rCopy) return; else m_vTemp = rCopy.m_vTemp; }

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XA_2DBorn<T> objects for unsupported types T
		//! Returns the xar::_eValueType corresponding to T
		static xar::_eValueType GetValuetype(void);

	// Operations
	public:
		//! Retrieves phase from 1 image of a pure phase object using 1st Born approximation
		void Born1(xar::XArray2D<T>& I1, double R, double dblAlpha = 0);
		//! Retrieves phase from 1 image of a pure phase object using 1st Born approximation optimizing over rAlphaOpt
		//void Born1Opt(xar::XArray2D<T>& I1, double R, double sigma, double& rAlphaOpt);
		//! Retrieves phase from 1 image of a pure phase object using iterative 1st Born approximation
		void Born1Iter(
			const xar::XArray2D<T>& I1, 
			xar::XArray2D<T>& P0, 
			double R, 
			double dblAlpha, 
			double dblSigma, 
			index_t& riNumIter);

		//! Retrieves phase from 1 image of a pure phase object using Rytov approximation
		void Rytov1(xar::XArray2D<T>& I1, double R, double dblAlpha = 0);
		//! Retrieves phase&amplitude from 1 image of a 'single-material' object using 1st Born approximation
		void BornSC(xar::XArray2D<T>& I1, double R, double delta2beta, double dblAlpha, bool bFullFresnelHom);
		//! Retrieves phase&amplitude from 1 image of a 'single-material' object using 1st Born approximation optimizing over rAlphaOpt
		//void BornSCOpt(xar::XArray2D<T>& I1, double R, double delta2beta, double sigma, double& rAlphaOpt);
		//! Retrieves phase&amplitude from 1 image of a 'single-material' object using iterative 1st Born approximation
		void BornSCIter(
			const xar::XArray2D<T>& I1, 
			xar::XArray2D<T>& I0, 
			double R, 
			double delta2beta, 
			double dblAlpha, 
			double dblSigma, 
			index_t& riNumIter);
		//! Retrieves phase&amplitude from 1 image of a 'single-material' object using Rytov approximation
		void RytovSC(xar::XArray2D<T>& I1, double R, double delta2beta, double dblAlpha = 0);
		//! Retrieves phase&amplitude from 2 images at different defocus distances using 1st Born approximation
		void Born2R(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double RI1, double RI2, double dblAlpha = 0);
		//! Retrieves phase&amplitude from 2 images at different defocus distances using 1st Born approximation optimizing over rAlphaOpt
		//void Born2ROpt(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double RI1, double RI2, double sigma, double& rAlphaOpt);
		//! Retrieves phase&amplitude from 2 images at different defocus distances using iterative 1st Born approximation
		void Born2RIter(
			const xar::XArray2D<T>& I1, 
			const xar::XArray2D<T>& I2, 
			xar::XArray2D<T>& I0, 
			xar::XArray2D<T>& P0, 
			double R1, 
			double R2, 
			double dblAlpha, 
			double dblSigma, 
			index_t& riNumIter);

		//! Retrieves phase&amplitude from 2 images at different defocus distances using Rytov approximation
		void Rytov2R(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double RI1, double RI2, double dblAlpha = 0);
		//! Retrieves phase&amplitude from 2 images at different defocus distances using iterative algorithm by A.Pogany
		void AP2RIter(
			const xar::XArray2D<T>& I1, 
			const xar::XArray2D<T>& I2, 
			xar::XArray2D<std::complex<T> >& C0, 
			double R1, 
			double R2, 
			double dblAlpha, 
			short niter, 
			double beta);

		//! Retrieves phase&amplitude from 2 images at different wavelengths (energies) using 1st Born approximation
		void Born2E(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double R, double dblAlpha = 0);
		//! Retrieves phase&amplitude from 2 images at different wavelengths (energies) using Rytov approximation
		void Rytov2E(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double R, double dblAlpha = 0);

		// Auxilliary function called by the optimization routine 'goldenC'
		//T CalcOnce(T* ArrayPars, index_t What2Calc);

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
		//! Auxilliary vector used for passing parameters between some member functions
		std::vector<void*> m_vTemp;
	// Member functions

	};

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//! Returns the xar::_eValueType corresponding to T=float
template<> inline xar::_eValueType XA_2DBorn<float>::GetValuetype() { return xar::eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> inline xar::_eValueType XA_2DBorn<double>::GetValuetype() { return xar::eXADouble; }

//! Retrieves phase from 1 image of a pure phase object using 1st Born approximation
// Solves Fresnel integral equation in 1st Born approximation for a pure phase object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed phase (more precisely, by Im(Psi))
// R = propagation distance R'
// dblAlpha (DEFAULT=0) is the regularisation parameter
//
// NOTE!!!: the program assumes that I1 is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iin, 1 must not be subtracted!!!)
// NOTE: this program uses the Ooura FFT library
// NOTE: this progam is a modification of PoissonFFT() - see invLa.h
template <class T> void XA_2DBorn<T>::Born1(xar::XArray2D<T>& I1, double R, double dblAlpha)
{

	IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	if (dblAlpha < 0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Born1 (negative regularization parameter)");

	if (dblAlpha == 0) 
		dblAlpha = 1.e-10; // zero dblAlpha leads to artifacts

	double wl = ph2->GetWl(); 
	double xlo = ph2->GetYlo();
	double xhi = ph2->GetYhi();
	double xst = GetYStep(I1);
	double ylo = ph2->GetXlo();
	double yhi = ph2->GetXhi();
	double yst = GetXStep(I1);

	double Iout = I1.Norm(xar::eNormAver); // Iout - average transmission coefficient
	
	if (Iout <= 0)
		throw std::invalid_argument("invalid_argument 'I1' in XA_2DBorn<T>::Born1 (non-positive average intensity)");

	I1 /= T(Iout);
	I1 -= 1.0; // NOTE: (I1-1) must be << 1 for the Born1 approximation to be valid (check it outside)
	I1 *= 0.5;
	
	index_t nxF = I1.GetDim1(), nx = nxF;
	index_t nyF = I1.GetDim2(), ny = nyF;

	// pad to the nearest power of 2
	index_t i = 2;
	while (i<nx) i *= 2;
	nx = i;
	index_t j = 2;
	while (j<ny) j *= 2;
	ny = j;
	if (nx != nxF || ny != nyF) 
	{
		xar::XArray2DMove<T> tmp(I1);
		tmp.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
	}

	T* arr = &(I1.front());
	std::vector<T> vecSpeq(nx*2);
	T* speq = &(vecSpeq[0]);

	// forward FFT
	xar::OouraFft<T> fft;

	fft.Real2D(arr, speq, nx, ny, xar::OouraFft<T>::eDirFwd);

	// multiplication of the FFT by 1/sin(C(m^2+n^2))
	// the indexing follows that from the function FFTRe(...) - see XA_fftr2.h
	double factor = 2. / (nx * ny);
	double a = xhi - xlo + xst;
	double b = yhi - ylo + yst;
	double a2 = xar::PI * wl * R / (a * a);
	double b2 = xar::PI * wl * R / (b * b);
	double dtemp, dtemp1;
	T temp;
	index_t nxd2 = nx / 2, i1;
	index_t nyd2 = ny / 2, j1;

	dtemp = b2 * nyd2 * nyd2;
	
	for (i=0; i<nxd2; i++) 
	{ 
		i1 = nxd2 - i;
		dtemp1 = sin(dtemp + a2 * i1 * i1);
		temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
		speq[nx+2*i] *= temp;
		speq[nx+2*i+1] *= temp;
	}

	for (i=nxd2; i<nx; i++) 
	{ 
		i1 = i - nxd2;
		dtemp1 = sin(dtemp + a2 * i1 * i1);
		temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
		speq[2*i-nx] *= temp;
		speq[2*i-nx+1] *= temp;
	}

	for (i=0; i<nxd2; i++)
	{
		i1 = nxd2 - i;
		dtemp = a2 * i1 * i1; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2;
			dtemp1 = sin(dtemp + b2 * j1 * j1);
			temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
			I1[i+nxd2][2*j-ny] *= temp; 
			I1[i+nxd2][2*j-ny+1] *= temp;
		}
	}

	I1[0][0] = I1[0][1] = T(0); // integral of psi is equal to zero
	
	for (i=nxd2+1; i<nx; i++)  // cycles for i>nxd2 (i1>0)
	{
		i1 = i - nxd2;
		dtemp = a2 * i1 * i1; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2;
			dtemp1 = sin(dtemp + b2 * j1 * j1);
			temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
			I1[i-nxd2][2*j-ny] *= temp; 
			I1[i-nxd2][2*j-ny+1] *= temp; 
		}
	}

	for (j=nyd2+1; j<ny; j++) // cycle for i==nxd2 (i1==0)
	{
		j1 = j - nyd2;
		dtemp1 = sin(b2 * j1 * j1);
		temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
		I1[0][2*j-ny] *= temp; 
		I1[0][2*j-ny+1] *= temp; 
	}

	//inverse FFT
	fft.Real2D(arr, speq, nx, ny, xar::OouraFft<T>::eDirInv);
	
	//trim back
	if (nx != nxF || ny != nyF) 
	{
		xar::XArray2DMove<T> tmp1(I1);
		tmp1.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
	}

	// the average value of Im(Psi) must be zero
	I1 -= T(I1.Norm(xar::eNormAver));
}


#if(0)
//! Retrieves phase from 1 image of a pure phase object using 1st Born approximation optimizing over rAlphaOpt
// Solves Fresnel integral equation in 1st Born approximation for a pure phase object
// using the FFT method and optimizing the regularization parameter dblAlpha;
// I1 = image; on exit I1 is replaced by the reconstructed phase
// R = propagation distance R'
// sigma = relative standard deviation of the noise corresponding to the average intensity in the image
// rAlphaOpt = reference to return the found optimal value of the regulariztion parameter
//
// NOTE!!!: the program assumes that I1 is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iinc, 1 must not be subtracted!!!)
// NOTE: this program uses the Ooura FFT library
// NOTE: this program also uses modified routine golden() from NumRecipes - see xa_nropt.h
// NOTE: this progam is a modification of PoissonFFT() - see invLa.h
template <class T> void XA_2DBorn<T>::Born1Opt(
	xar::XArray2D<T>& I1, 
	double R, 
	double sigma, 
	double& rAlphaOpt)
{
	const IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	if (sigma <= 0)
		throw std::invalid_argument("invalid_argument 'sigma' in XA_2DBorn<T>::Born1Opt (must be positive)");

	if (rAlphaOpt <= 0)
		throw std::invalid_argument("invalid_argument 'rAlphaOpt' in XA_2DBorn<T>::Born1Opt (must be positive)");

	m_vTemp.resize(1, 0);
	m_vTemp[0] = &I1; // use m_pTemp[0] to pass the image to CalcOnce

	std::vector<T> vecArgs(6);	// vector of parameters to be passed to 'CalcOnce' via 'goldenC'
	T* mass = &(vecArgs.front());
	mass[0] = T(vecArgs.size() - 1); // number of other parameters in this array
	mass[1] = T(log10(rAlphaOpt));	// base10-log of the regularization parameter
	mass[2] = T(R);	// defocus distance
	mass[3] = T(sigma * sigma * I1.GetDim1() * I1.GetDim2()); // targeted reconstruction accuracy
	mass[4] = T(I1.Norm(xar::eNormMin));
	mass[5] = T(I1.Norm(xar::eNormAver));

	xar::index_t num_varied = 1;	// logAlpha is the argument (no.1) to be optimized by 'goldenC'
	T tLogAlpha = mass[1];	// a median LogAlpha
	T tLogAlphaMin = mass[1] - 2; // the minimum LogAlpha
	T tLogAlphaMax = mass[1] + 2;	// the maximum LogAlpha
	T tTolerance = T(0.2); // optimization tolerance
	T tAlphaOpt = 0;		// storage for the returned optimum dblAlpha

	T tRecAccuracy = xar::XA_Optimisation::GoldenMinimumSearch< XA::XA_2DBorn<T> >(
		tLogAlphaMin, 
		tLogAlpha, 
		tLogAlphaMax,
		*this, 
		mass, 
		num_varied, 
		tTolerance, 
		&tAlphaOpt, 
		0);

	// do the phase retrieval with the found optimum dblAlpha
	Born1(I1, R, pow(10.0, double(tAlphaOpt)));
	rAlphaOpt = pow(10.0, double(tAlphaOpt));
}



// Auxilliary function called by the optimization routine 'goldenC'
// What2Calc determines which member function is called for the caller optimization routine
// ArrayPars is the array of parameters to be passed to a member function
template <class T> T XA_2DBorn<T>::CalcOnce(T* ArrayPars, index_t What2Calc)
{
	switch (What2Calc)
	{
	case 0: // corresponds to Born1Opt
	{
		// ArrayPars[1] = log(dblAlpha)	// base10-log of the regularization parameter
		// ArrayPars[2] = R		// defocus distance
		// ArrayPars[3] = dblSigma2	// sigma * sigma * I1.GetDim1() * I1.GetDim2();
		// ArrayPars[4] = dblImageMin; // minimum of image values
		// ArrayPars[5] = dblImageAver; // average of image values
		if (ArrayPars[0] != 5)
			throw std::invalid_argument("invalid_argument 'ArrayPars' in XA_2DBorn<T><T>::CalcOnce (wrong number of parameters)"); 

		xar::XArray2D<T>* pI1 = reinterpret_cast<xar::XArray2D<T>*>(m_vTemp[0]);
		xar::XArray2D<T> I1 = *pI1; // restore the defocused image
		Born1(I1, double(ArrayPars[2]), pow(10.0, double(ArrayPars[1]))); // reconstruct object-plane phase
		xar::XArray2D< std::complex<T> > C0; // create complex object-plane amplitude
		MakeComplex(T(sqrt(ArrayPars[5])), I1, C0, true);
		xar::XArray2DFFT<T> CFFTService(C0);
		CFFTService.Kirchhoff(double(ArrayPars[2])); // calculate image
		Abs2(C0, I1);
		
		// the average value of the image intensity should stay constant
		assert(fabs(ArrayPars[5] - I1.Norm(xar::eNormAver)) < 0.01 * ArrayPars[5]);
		
		double ret;
		// we only compare the central quarter of the images
		xar::XArray2D<T> I1a = *pI1;
		xar::XArray2DMove<T> I1aMove(I1a);
		I1aMove.Trim(I1a.GetDim1()/4, I1a.GetDim1()/4, I1a.GetDim2()/4, I1a.GetDim2()/4);
		xar::XArray2DMove<T> I1Move(I1);
		I1Move.Trim(I1.GetDim1()/4, I1.GetDim1()/4, I1.GetDim2()/4, I1.GetDim2()/4);
		
		if (ArrayPars[4] > 0) // no zeros or negative values in the image
			//ret = fabs(pI1->Chi2(I1, 1.0, true) / double(ArrayPars[3]) - 1.0);
			ret = fabs(I1a.Chi2(I1, 1.0, true) / double(ArrayPars[3] / 4) - 1.0);
		else // non-Poisson statistics (sigma is the same at all points)
			//ret = fabs(pI1->Chi2(I1, 1.0, false) / double(ArrayPars[3]) - 1.0);
			ret = fabs(I1a.Chi2(I1, 1.0, false) / double(ArrayPars[3]/ 4) - 1.0);
		return T(ret);
	}
		break;
	case 1: // corresponds to BornSCOpt
		{
		// ArrayPars[1] = log(dblAlpha)	// base10-log of the regularization parameter
		// ArrayPars[2] = R		// defocus distance
		// ArrayPars[3] = delta2beta // delta/beta
		// ArrayPars[4] = dblSigma2	// sigma * sigma * I1.GetDim1() * I1.GetDim2();
		// ArrayPars[5] = dblImageMin; // minimum of image values
		// ArrayPars[6] = dblImageAver; // average of image values
		if (ArrayPars[0] != 6)
			throw std::invalid_argument("invalid_argument 'ArrayPars' in XA_2DBorn<T><T>::CalcOnce (wrong number of parameters)");

		xar::XArray2D<T>* pI1 = reinterpret_cast<xar::XArray2D<T>*>(m_vTemp[0]);
		xar::XArray2D<T> I1 = *pI1; // restore the defocused image

		BornSC(I1, double(ArrayPars[2]), double(ArrayPars[3]), pow(10.0, double(ArrayPars[1]))); // reconstruct object-plane phase

		xar::XArray2D<T> P0(I1);

		I1 ^= T(0.5); // intensity --> amplitude
		// ImPsi0 = delta2beta * (I0 / Iout - 1) / 2;
		P0.Log(); P0 *= T(0.5 * ArrayPars[3]); // intensity --> phase
		xar::XArray2D< std::complex<T> > C0; 
		MakeComplex(I1, P0, C0, true); // create complex object-plane amplitude
		xar::XArray2DFFT<T> CFFTService(C0);
		CFFTService.Kirchhoff(double(ArrayPars[2])); // calculate image
		Abs2(C0, I1);

		// the average value of the image intensity should stay constant
		assert(fabs(ArrayPars[6] - I1.Norm(xar::eNormAver)) < 0.01 * ArrayPars[6]);

		double ret;
		// we only compare the central quarter of the images
		P0 = *pI1;
		xar::XArray2DMove<T> P0Move(P0);
		P0Move.Trim(P0.GetDim1()/4, P0.GetDim1()/4, P0.GetDim2()/4, P0.GetDim2()/4);
		xar::XArray2DMove<T> I1Move(I1);
		I1Move.Trim(I1.GetDim1()/4, I1.GetDim1()/4, I1.GetDim2()/4, I1.GetDim2()/4);

		if (ArrayPars[5] > 0) // no zeros or negative values in the image
			//ret = fabs(pI1->Chi2(I1, 1.0, true) / double(ArrayPars[4]) - 1.0);
			ret = fabs(P0.Chi2(I1, 1.0, true) / double(ArrayPars[4] / 4) - 1.0);
		else // non-Poisson statistics (sigma is the same at all points)
			//ret = fabs(pI1->Chi2(I1, 1.0, false) / double(ArrayPars[4]) - 1.0);
			ret = fabs(P0.Chi2(I1, 1.0, false) / double(ArrayPars[4]/ 4) - 1.0);
		return T(ret);
		}
		break;
	case 2: // corresponds to Born2ROpt
		{
		// ArrayPars[1] = log(dblAlpha)	// base10-log of the regularization parameter
		// ArrayPars[2] = RI1		// first defocus distance
		// ArrayPars[3] = RI2		// second defocus distance
		// ArrayPars[4] = dblSigma2	// sigma * sigma * I1.GetDim1() * I1.GetDim2();
		// ArrayPars[5] = dblImageMin; // minimum of image values
		// ArrayPars[6] = dblImageAver; // average of image values
		if (ArrayPars[0] != 6)
			throw std::invalid_argument("invalid_argument 'ArrayPars' in XA_2DBorn<T><T>::CalcOnce (wrong number of parameters)"); 

		xar::XArray2D<T>* pI1 = reinterpret_cast<xar::XArray2D<T>*>(m_vTemp[0]);
		xar::XArray2D<T> I1 = *pI1; // restore the defocused image
		xar::XArray2D<T>* pI2 = reinterpret_cast<xar::XArray2D<T>*>(m_vTemp[1]);
		xar::XArray2D<T> I2 = *pI2; // restore the defocused image
		
		Born2R(I1, I2, double(ArrayPars[2]), double(ArrayPars[3]), pow(10.0, double(ArrayPars[1]))); // reconstruct object-plane phase

		if (ArrayPars[2] == 0) 
			I2 = *pI1; // restore the object-plane intensity

		if (ArrayPars[3] == 0) 
			I2 = *pI2; // restore the object-plane intensity

		I2 ^= T(0.5); // intensity --> amplitude
		xar::XArray2D<T> I1z, I2z;
		xar::XArray2D< std::complex<T> > C0; 
		MakeComplex(I2, I1, C0, true); // create complex object-plane amplitude
		xar::XArray2DFFT<T> CFFTService(C0);
		CFFTService.Kirchhoff(double(ArrayPars[2])); // calculate first image
		Abs2(C0, I1z);
		MakeComplex(I2, I1, C0, true); // restore complex object-plane amplitude
		CFFTService.Kirchhoff(double(ArrayPars[3])); // calculate second image
		Abs2(C0, I2z);

		// the average value of the image intensity should stay constant
		assert(fabs(ArrayPars[6] - I1z.Norm(xar::eNormAver)) < 0.01 * ArrayPars[6]);
		assert(fabs(ArrayPars[6] - I2z.Norm(xar::eNormAver)) < 0.01 * ArrayPars[6]);

		double ret = 0;
		// we only compare the central quarter of the images
		I1 = *pI1;
		xar::XArray2DMove<T> I1Move(I1);
		I1Move.Trim(I1.GetDim1()/4, I1.GetDim1()/4, I1.GetDim2()/4, I1.GetDim2()/4);
		xar::XArray2DMove<T> I1zMove(I1z);
		I1zMove.Trim(I1z.GetDim1()/4, I1z.GetDim1()/4, I1z.GetDim2()/4, I1z.GetDim2()/4);
		
		if (ArrayPars[5] > 0) // no zeros or negative values in the image
			//ret = fabs(pI1->Chi2(I1, 1.0, true) / double(ArrayPars[4]) - 1.0);
			ret = fabs(I1.Chi2(I1z, 1.0, true) / double(ArrayPars[4] / 4) - 1.0);
		else // non-Poisson statistics (sigma is the same at all points)
			//ret = fabs(pI1->Chi2(I1, 1.0, false) / double(ArrayPars[4]) - 1.0);
			ret = fabs(I1.Chi2(I1z, 1.0, false) / double(ArrayPars[4]/ 4) - 1.0);
		
		I2 = *pI2;
		xar::XArray2DMove<T> I2Move(I2);
		I2Move.Trim(I2.GetDim1()/4, I2.GetDim1()/4, I2.GetDim2()/4, I2.GetDim2()/4);
		xar::XArray2DMove<T> I2zMove(I2z);
		I2zMove.Trim(I2z.GetDim1()/4, I2z.GetDim1()/4, I2z.GetDim2()/4, I2z.GetDim2()/4);
		
		if (ArrayPars[5] > 0) // no zeros or negative values in the image
			//ret = fabs(pI2->Chi2(I2, 1.0, true) / double(ArrayPars[4]) - 1.0);
			ret = std::max(ret, fabs(I2.Chi2(I2z, 1.0, true) / double(ArrayPars[4] / 4) - 1.0));
		else // non-Poisson statistics (sigma is the same at all points)
			//ret = fabs(pI2->Chi2(I2, 1.0, false) / double(ArrayPars[4]) - 1.0);
			ret = std::max(ret, fabs(I2.Chi2(I2z, 1.0, false) / double(ArrayPars[4]/ 4) - 1.0));
		
		return T(ret);
		}
		break;
	default:
		throw std::invalid_argument("invalid_argument 'What2Calc' in XA_2DBorn<T>::CalcOnce (value not allowed)"); 
	}
}
#endif


//! Retrieves phase from 1 image of a pure phase object using iterative 1st Born approximation
// Attempts to iteratively improve initial phase guess using some initial approximation for a pure phase object
// and solving the equation I1 = I1(0) + 2Re{ I0(0) * Psi] }, where U0 = U0(0)(1 + Psi),
// using Born1 method
//
// I1 = image;
// P0 = initial approximation for the phase; on exit P0 is replaced by the corrected phase
// R = propagation distance R'
// dblAlpha = regularization parameter
// dblSigma = noise level
// riNumIter = reference to the number of iterations parameter
//
// NOTE!!!: the program assumes that I1 is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iinc, 1 must not be subtracted!!!)
template <class T> void XA_2DBorn<T>::Born1Iter(
	const xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& P0, 
	double R, 
	double dblAlpha, 
	double dblSigma, 
	index_t& riNumIter)
{
	const IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	if (dblSigma<0)
		throw std::invalid_argument("invalid_argument 'dblSigma' in XA_2DBorn<T>::Born1Iter (negative noise level)");

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Born1Iter (negative regularization parameter)");

	if (dblAlpha==0) 
		dblAlpha = 1.e-10; // zero dblAlpha leads to artifacts

	double dblNormdP0 = P0.Norm(xar::eNormL2), dblNormTemp;
	double dblSigma2 = dblSigma * dblSigma * I1.GetDim1() * I1.GetDim2();
	double dblImageMin = I1.Norm(xar::eNormMin);
	double dblTolerance = 0.01; // here we always expect the chi-2 to be decreasing
	double norma_resid, norma_resid0;
	double Iout = I1.Norm(xar::eNormAver);
	xar::XArray2D<T> I0(I1), dI, dP0;
	xar::XArray2D< std::complex<T> > C0;
	xar::XArray2DFFT<T> CFFTService(C0);
	I0.Fill(T(sqrt(Iout)));

	for (index_t i=0; i<riNumIter; i++)
	{
		MakeComplex(I0, P0, C0, true);
		CFFTService.Kirchhoff(R);
		Abs2(C0, dI);
		if (dblSigma > 0.0) // check if the residuals are comparable with the noise
		{
			if (dblImageMin > 0) // no zeros or negative values in the image
				norma_resid = I1.Chi2(dI, 1.0, true) / dblSigma2 - 1.0;
			else // non-Poisson statistics (sigma is the same at all points)
				norma_resid = I1.Chi2(dI, 1.0, false) / dblSigma2 - 1.0;
			
			if (norma_resid < dblTolerance) 
			{
				riNumIter = i;
				break;
			}
			
			if (i > 0 && norma_resid > norma_resid0) // last correction was bad 
			{
				P0 -= dP0; // new phase roll back
				riNumIter = i - 1;
				break;
			}
			norma_resid0 = norma_resid;
		}
		dI -= I1;
		dI *= -T(1.0);
		dI += T(1.0);
		Born1(dI, R, dblAlpha);
		dblNormTemp = dI.Norm(xar::eNormL2);
		
		if (dblNormTemp >= dblNormdP0) // convergence test
		{
			riNumIter = i;
			break;
		}
		
		dblNormdP0 = dblNormTemp;
		dP0 = dI; P0 += dI;
	}
}



//! Retrieves phase from 1 image of a pure phase object using Rytov approximation
// Solves Fresnel integral equation in Rytov approximation for a pure phase object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed phase
// R = propagation distance R'
// dblAlpha (DEFAULT=0) is the regularisation parameter
// NOTE: this program uses Ooura FFT library
template <class T> void XA_2DBorn<T>::Rytov1(xar::XArray2D<T>& I1, double R, double dblAlpha)
{
	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Rytov1 (negative regularization parameter)");
	
	if (I1.Norm(xar::eNormMin) <= 0)
		throw std::invalid_argument("invalid_argument 'I1' in XA_2DBorn<T>::Rytov1 (non-positive image intensity)");
	
	double Iout = I1.Norm(xar::eNormAver); // Iout - average transmission coefficient

	// transform I1 into J such that J / Iout - 1 = log(I1 / Iout), i.e. J = Iout * [log(I1 / Iout) + 1]
	I1 /= T(Iout);
	I1.Log();
	I1 += T(1);
	I1 *= T(Iout);
	
	Born1(I1, R, dblAlpha);
}



//! Retrieves phase&amplitude from 1 image of a 'single-material' object using 1st Born approximation
// Solves Fresnel integral equation in the first Born approximation for a single component object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed INTENSITY (more precisely, by Iout * (1 + 2 * Re(psi)))
// R = propagation distance R'
// delta2beta = delta / beta
// dblAlpha is the regularisation parameter
// bFullFresnel - if true will calculate full Fresnel-Hom retrieval, else will calculate 1stBorn-Hom retrieval
// NOTE!!!: the program assumes that I1 is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// NOTE!!!: if delta2beta==0, then this programs does 'pure absorption' retrieval
// NOTE: this program uses the Ooura FFT library
template <class T> void XA_2DBorn<T>::BornSC(
	xar::XArray2D<T>& I1, 
	double R, 
	double delta2beta, 
	double dblAlpha,
	bool bFullFresnelHom)
{
	if (delta2beta<0)
		throw std::invalid_argument("invalid_argument 'delta2beta' in XA_2DBorn<T>::BornSC (negative delta/beta parameter)");
	
	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::BornSC (negative regularization parameter)");
	
	if (dblAlpha==0) 
		dblAlpha = 1.e-10;
	
	double Iout = I1.Norm(xar::eNormAver); // Iout - average bulk transmission coefficient
	T tIout = T(Iout);
	
	double matangamma = -atan(delta2beta);

	IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	double wl = ph2->GetWl(); 
	double xlo = ph2->GetYlo();
	double xhi = ph2->GetYhi();
	double xst = GetYStep(I1); 
	double ylo = ph2->GetXlo();
	double yhi = ph2->GetXhi();
	double yst = GetXStep(I1);

	//I1 /= tIout;	
	//I1 -= T(1.0);
	//I1 *= T(0.5 / sqrt(1 + delta2beta * delta2beta));
	T* arrI1 = &I1.front();
	T tAbra = T(0.5 / sqrt(1 + delta2beta * delta2beta)); 
	if (bFullFresnelHom)
	{
		tAbra *= T(2);
		for (index_t i = 0; i < I1.size(); i++) arrI1[i] *= tAbra;
	}
	else
		for (index_t i = 0; i < I1.size(); i++) arrI1[i] = (arrI1[i] / tIout - T(1.0)) * tAbra;
	
	index_t nxF = I1.GetDim1(), nx = nxF;
	index_t nyF = I1.GetDim2(), ny = nyF;

	// pad to the nearest power of 2
	index_t i = 2;
	while (i<nx) i *= 2;
	nx = i;
	index_t j = 2;
	while (j<ny) j *= 2;
	ny = j;
	
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp(I1);
		tmp.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
	}

	T* arr = &(I1.front());
	std::vector<T> vecSpeq(nx*2);
	T* speq = &(vecSpeq[0]);

	// forward FFT
	xar::OouraFft<T> fft;

	fft.Real2D(arr, speq, nx, ny, xar::OouraFft<T>::eDirFwd);

	// multiplication of the FFT by c0/(c1*cos(pi*wl*z(m^2+n^2))-c2)
	// the indexing follows that from the function FFTRe(...) - see vo_fftr2.h
	double factor = 2. / (nx * ny);
	double a = xhi - xlo + xst;
	double b = yhi - ylo + yst;
	double a2 = xar::PI * wl * R / (a * a); // / factor;
	double b2 = xar::PI * wl * R / (b * b); // / factor;
	double dtemp, dtemp1;
	T temp;
	index_t nxd2 = nx / 2, i1;
	index_t nyd2 = ny / 2, j1;

	dtemp = b2 * nyd2 * nyd2 + matangamma;
	for (i=0; i<nxd2; i++) 
	{ 
		i1 = nxd2 - i;
		dtemp1 = cos(dtemp + a2 * i1 * i1);
		temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
		speq[nx+2*i] *= temp;
		speq[nx+2*i+1] *= temp;
	}

	for (i=nxd2; i<nx; i++) 
	{ 
		i1 = i - nxd2;
		dtemp1 = cos(dtemp + a2 * i1 * i1);
		temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
		speq[2*i-nx] *= temp;
		speq[2*i-nx+1] *= temp;
	}
		
	for (i=0; i<nxd2; i++)
	{
		i1 = nxd2 - i;
		dtemp = a2 * i1 * i1 + matangamma; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2;
			dtemp1 = cos(dtemp + b2 * j1 * j1);
			temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
			I1[i+nxd2][2*j-ny] *= temp; 
			I1[i+nxd2][2*j-ny+1] *= temp;
		}
	}

	//I1[0][0] = I1[0][1] = T(0); // integral of psi is equal to zero
	for (i=nxd2+1; i<nx; i++)  // cycles for i>nxd2 (i1>0)
	{
		i1 = i - nxd2;
		dtemp = a2 * i1 * i1 + matangamma; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2;
			dtemp1 = cos(dtemp + b2 * j1 * j1);
			temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
			I1[i-nxd2][2*j-ny] *= temp; 
			I1[i-nxd2][2*j-ny+1] *= temp; 
		}
	}

	for (j=nyd2; j<ny; j++) // cycle for i==nxd2 (i1==0)
	{
		j1 = j - nyd2;
		dtemp1 = cos(b2 * j1 * j1 + matangamma);
		temp = T(dtemp1 / (dblAlpha + dtemp1 * dtemp1) * factor);
		I1[0][2*j-ny] *= temp; 
		I1[0][2*j-ny+1] *= temp; 
	}

	//inverse FFT
	fft.Real2D(arr, speq, nx, ny, xar::OouraFft<T>::eDirInv);
	
	//trim back
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp1(I1);
		tmp1.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
	}

	// Output is I = A^2 * (1 + 2 * Re(psi))
	// the average value of Re(psi) must be zero
	//I1 -= I1.Norm(xar::eNormAver);
	//I1 *= T(2);
	//I1 += T(1);
	//I1 *= T(Iout);
	tAbra = T(I1.Norm(xar::eNormAver));
	arrI1 = &I1.front();
	if (!bFullFresnelHom)
		for (index_t i = 0; i < I1.size(); i++) arrI1[i] = ((arrI1[i] - tAbra) * T(2.0) + T(1.0)) * tIout;
}

#if(0)
//! Retrieves phase&amplitude from 1 image of a 'single-material' object using iterative 1st Born approximation optimizing over rAlphaOpt
// Solves Fresnel integral equation in 1st Born approximation for a single component object
// using the FFT method and optimizing the regularization parameter dblAlpha;
// I1 = image; on exit I1 is replaced by the reconstructed INTENSITY (!!!)
// R = propagation distance R'
// delta2beta = delta / beta
// sigma = relative standard deviation of the noise corresponding to the average intensity in the image
// rAlphaOpt = reference to return the found optimal value of the regulariztion parameter
//
// NOTE!!!: the program assumes that I1 is background-corrected, i.e. DIVIDED by the INCIDENT intensity
// (F=I1/Iinc, 1 must not be subtracted!!!)
// NOTE: this program uses the Ooura FFT library
// NOTE: this program also uses modified routine golden() from NumRecipes - see xa_nropt.h
template <class T> void XA_2DBorn<T>::BornSCOpt(
	xar::XArray2D<T>& I1, 
	double R, 
	double delta2beta, 
	double sigma, 
	double& rAlphaOpt)
{
	const IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	if (delta2beta<0)
		throw std::invalid_argument("invalid_argument 'delta2beta' in XA_2DBorn<T>::BornSCOpt (negative delta/beta parameter)");

	if (sigma <= 0)
		throw std::invalid_argument("invalid_argument 'sigma' in XA_2DBorn<T>::BornSCOpt (must be positive)");

	if (rAlphaOpt <= 0)
		throw std::invalid_argument("invalid_argument 'rAlphaOpt' in XA_2DBorn<T>::BornSCOpt (must be positive)");

	m_vTemp.resize(1, 0);
	m_vTemp[0] = &I1; // use m_pTemp[0] to pass the image to CalcOnce

	std::vector<T> vecArgs(7);	// vector of parameters to be passed to 'CalcOnce' via 'goldenC'
	T* mass = &(vecArgs.front());
	mass[0] = T(vecArgs.size() - 1); // number of other parameters in this array
	mass[1] = T(log10(rAlphaOpt));	// base10-log of the regularization parameter
	mass[2] = T(R);	// defocus distance
	mass[3] = T(delta2beta);	// delta/beta
	mass[4] = T(sigma * sigma * I1.GetDim1() * I1.GetDim2()); // targeted reconstruction accuracy
	mass[5] = T(I1.Norm(xar::eNormMin));
	mass[6] = T(I1.Norm(xar::eNormAver));

	index_t num_varied = 1;	// logAlpha is the argument (no.1) to be optimized by 'goldenC'
	T tLogAlpha = mass[1];	// a median LogAlpha
	T tLogAlphaMin = mass[1] - 2; // the minimum LogAlpha
	T tLogAlphaMax = mass[1] + 2;	// the maximum LogAlpha
	T tTolerance = T(0.2); // optimization tolerance
	T tAlphaOpt = 0;		// storage for the returned optimum dblAlpha

	T tRecAccuracy = xar::XA_Optimisation::GoldenMinimumSearch< XA::XA_2DBorn<T> >(
		tLogAlphaMin, 
		tLogAlpha, 
		tLogAlphaMax,
		*this, 
		mass, 
		num_varied, 
		tTolerance, 
		&tAlphaOpt, 
		1);

	// do the phase retrieval with the found optimum dblAlpha
	BornSC(I1, R, delta2beta, pow(10.0, double(tAlphaOpt)));
	rAlphaOpt = pow(10.0, double(tAlphaOpt));
}
#endif

//! Retrieves phase&amplitude from 1 image of a 'single-material' object using iterative 1st Born approximation
// Attempts to iteratively improve initial phase guess for a single component object
// and solving the equation I1 ~ I1(0) + 2Re{ I0(0) * Psi] }, where U0 = U0(0)(1 + Psi),
// using the BornSC method
//
// I1 = image;
// I0 = initial approximation for the intensity; on exit I0 is replaced by the corrected INTENSITY
// R = propagation distance R'
// delta2beta = delta / beta
// dblAlpha = regularization parameter
// dblSigma = noise level
// riNumIter = reference to the number of iterations parameter
template <class T> void XA_2DBorn<T>::BornSCIter(
	const xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& I0, 
	double R, 
	double delta2beta, 
	double dblAlpha, 
	double dblSigma, 
	index_t& riNumIter)
{
	if (delta2beta<0)
		throw std::invalid_argument("invalid_argument 'delta2beta' in XA_2DBorn<T>::BornSCIter (negative delta/beta parameter)");

	if (dblSigma<0)
		throw std::invalid_argument("invalid_argument 'dblSigma' in XA_2DBorn<T>::BornSCIter (negative noise level)");

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::BornSCIter (negative regularization parameter)");

	if (dblAlpha==0) 
		dblAlpha = 1.e-10;

	double Iout = I0.Norm(xar::eNormAver);
	
	if (fabs(Iout - I1.Norm(xar::eNormAver)) > 0.01 * Iout)
		throw std::invalid_argument("invalid_argument 'I0' in XA_2DBorn<T>::BornSCIter (different aver.intensity from I1)");

	xar::XArray2D<T> A0(I0), P0(I0), dI, dP0;
	A0 ^= T(0.5); // intensity --> amplitude
	P0.Log(); P0 *= T(0.5 * delta2beta); // intensity --> phase
	double dblNormdP0 = P0.Norm(xar::eNormL2), dblNormTemp;
	double dblSigma2 = dblSigma * dblSigma * I1.GetDim1() * I1.GetDim2();
	double dblImageMin = I1.Norm(xar::eNormMin);
	double dblTolerance = 0.01; // here we always expect the chi-2 to be decreasing
	double norma_resid, norma_resid0(0);
	xar::XArray2D< std::complex<T> > C0;
	xar::XArray2DFFT<T> CFFTService(C0);
	//XArData XArData; // @@@@@ for testing purposes only

	for (index_t i = 0; i < riNumIter; i++)
	{
		MakeComplex(A0, P0, C0, true);
		CFFTService.Kirchhoff(R);
		Abs2(C0, dI);

		if (dblSigma > 0.0) // check if the residuals are comparable with the noise
		{
			if (dblImageMin > 0) // no zeros or negative values in the image
				norma_resid = I1.Chi2(dI, 1.0, true) / dblSigma2 - 1.0;
			else // non-Poisson statistics (sigma is the same at all points)
				norma_resid = I1.Chi2(dI, 1.0, false) / dblSigma2 - 1.0;
			if (norma_resid < dblTolerance) 
			{
				riNumIter = i;
				break;
			}
			
			if (i > 0 && norma_resid > norma_resid0) // last correction was bad 
			{
				dI = dP0; P0 -= dI; // new phase roll back
				dI /= T(delta2beta); dI.Exp(); A0 /= dI; // new amplitude roll back
				I0 = A0; I0 ^= 2.0; // amplitude --> intensity	
				riNumIter = i - 1;
				break;
			}
			norma_resid0 = norma_resid;
		}

		dI -= I1;
		dI *= -T(1.0);
		dI += T(1.0);
		BornSC(dI, R, delta2beta, dblAlpha);
		dI -= T(1.0); dI *= T(0.5); dI /= I0; // dI becomes Re(psi0)
		dI *= T(delta2beta); 

		dblNormTemp = dI.Norm(xar::eNormL2);
		//if (dblNormTemp >= dblNormdP0) // convergence test
		//{
		//	riNumIter = i;
		//	break;
		//}

		dblNormdP0 = dblNormTemp; 
		dP0 = dI; 

		P0 += dI; // new phase
		dI /= T(delta2beta); 
		dI += T(1.0); 
		A0 *= dI; // new amplitude
		I0 = A0; 
		I0 ^= 2.0; // amplitude --> intensity	
	}
}



//! Retrieves phase&amplitude from 1 image of a 'single-material' object using Rytov approximation
// Solves Fresnel integral equation in Rytov approximation for a single component object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed INTENSITY (!!!)
// R = propagation distance R'
// delta2beta = delta / beta
// dblAlpha (DEFAULT=0) is the regularisation parameter
// NOTE!!!: if delta2beta==0, then this programs does 'pure absorption' retrieval
// NOTE: this program uses the Ooura FFT library
template <class T> void XA_2DBorn<T>::RytovSC(xar::XArray2D<T>& I1, double R, double delta2beta, double dblAlpha)
{
	if (delta2beta<0)
		throw std::invalid_argument("invalid_argument 'delta2beta' in XA_2DBorn<T>::RytovSC (negative delta/beta parameter)");

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::RytovSC (negative regularization parameter)");

	if (I1.Norm(xar::eNormMin) <= 0)
		throw std::invalid_argument("invalid_argument 'I1' in XA_2DBorn<T>::RytovSC (non-positive image intensity)");
	
	double Iout = I1.Norm(xar::eNormAver); // Iout - average transmission coefficient

	// transform I1 into J such that J / Iout - 1 = log(I1 / Iout), i.e. J = Iout * [log(I1 / Iout) + 1]
	I1 /= T(Iout);
	I1.Log();
	I1 += T(1);
	I1 *= T(Iout);

	BornSC(I1, R, delta2beta, dblAlpha);
}



//! Retrieves phase&amplitude from 2 images at different defocus distances using 1st Born approximation
// Solves Fresnel integral equation in the first Born approximation for an amplitude-phase object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed phase (more precisely, by Im(psi0))
// I2 = image; on exit I2 is replaced by the reconstructed intensity (more precisely, by Iout * (1 + 2 * Re(psi)))
// RI1 = propagation distance R' corresponding to image I1
// RI2 = propagation distance R' corresponding to image I2
// dblAlpha (DEFAULT=0) is the regularisation parameter
//
// NOTE!!!: the program assumes that I1 and I2 are background-corrected, i.e. DIVIDED by the INCIDENT intensity
// NOTE: this program uses the Ooura FFT library
template <class T> void XA_2DBorn<T>::Born2R(xar::XArray2D<T>& I1, xar::XArray2D<T>& I2, double RI1, double RI2, double dblAlpha)
{

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Born2R (negative regularization parameter)");
	if (dblAlpha==0) dblAlpha = 1.e-10;

	IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	double wl = ph2->GetWl(); 
	double xlo = ph2->GetYlo();
	double xhi = ph2->GetYhi();
	double xst = GetYStep(I1);
	double ylo = ph2->GetXlo();
	double yhi = ph2->GetXhi();
	double yst = GetXStep(I1);

	double Iout = I1.Norm(xar::eNormAver);

	if (fabs(Iout - I2.Norm(xar::eNormAver)) > 0.01 * fabs(Iout))
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Born2R (different averages)");

	I1 /= T(Iout);
	I1 -= 1.0;
	I1 *= 0.5;
	I2 /= T(Iout);
	I2 -= 1.0;
	I2 *= 0.5;
	
	index_t nxF = I1.GetDim1(), nx = nxF;
	index_t nyF = I1.GetDim2(), ny = nyF;

	// pad to the nearest power of 2
	index_t i = 2;
	while (i<nx) i *= 2;
	nx = i;
	index_t j = 2;
	while (j<ny) j *= 2;
	ny = j;
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp1(I1);
		tmp1.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
		xar::XArray2DMove<T> tmp2(I2);
		tmp2.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
	}

		
	T* arrI1 = &(I1.front());
	std::vector<T> vecSpeqI1(nx*2);
	T* speqI1 = &(vecSpeqI1[0]);
	T* arrI2 = &(I2.front());
	std::vector<T> vecSpeqI2(nx*2);
	T* speqI2 = &(vecSpeqI2[0]);

	
	// forward FFTs
	xar::OouraFft<T> fft;

	fft.Real2D(arrI1, speqI1, nx, ny, xar::OouraFft<T>::eDirFwd);
	fft.Real2D(arrI2, speqI2, nx, ny, xar::OouraFft<T>::eDirFwd);


	// multiplication of the FFT by M^(-1)
	// the indexing follows that from the function FFTRe(...) - see vo_fftr2.h
	T temp;
	double dtemp, dtemp1, dtemp2, dtempd, sa1, sa2, sd, ca1, ca2, arg1, arg2;
	double dR = RI1 - RI2;
	double factor = 2. / (nx * ny);
	double a = xhi - xlo + xst;
	double b = yhi - ylo + yst;
	dtemp = xar::PI * wl;
	dtemp1 = 1.0 / (a * a);
	dtemp2 = 1.0 / (b * b);
	double a1 = dtemp * RI1 * dtemp1;
	double b1 = dtemp * RI1 * dtemp2;
	double a2 = dtemp * RI2 * dtemp1;
	double b2 = dtemp * RI2 * dtemp2;
	double ad = dtemp * dR * dtemp1;
	double bd = dtemp * dR * dtemp2;

	index_t nxd2 = nx / 2, i1, ii, ii1;
	index_t nyd2 = ny / 2, j1, jj, jj1;

	dtemp = double(nyd2 * nyd2);
	dtemp1 = b1 * dtemp;
	dtemp2 = b2 * dtemp;
	dtempd = bd * dtemp;

	for (i=0; i<nxd2; i++) 
	{ 
		i1 = nxd2 - i; i1 *= i1;
		arg1 = dtemp1 + a1 * i1;
		sa1 = sin(arg1); ca1 = -cos(arg1);
		arg2 = dtemp2 + a2 * i1;
		sa2 = -sin(arg2); ca2 = cos(arg2);
		sd = sin(dtempd + ad * i1); sd *= factor / (dblAlpha + sd * sd);
		ii = nx + 2 * i; ii1 = ii + 1;
		speqI1[ii] *= T(sd); speqI1[ii1] *= T(sd);
		speqI2[ii] *= T(sd); speqI2[ii1] *= T(sd);
		temp = T(ca2 * speqI1[ii] + ca1 * speqI2[ii]);
		speqI2[ii] = T(sa2 * speqI1[ii] + sa1 * speqI2[ii]);
		speqI1[ii] = temp;
		temp = T(ca2 * speqI1[ii1] + ca1 * speqI2[ii1]);
		speqI2[ii1] = T(sa2 * speqI1[ii1] + sa1 * speqI2[ii1]);
		speqI1[ii1] = temp;
	}

	for (i=nxd2; i<nx; i++) 
	{ 
		i1 = i - nxd2; i1 *= i1;
		arg1 = dtemp1 + a1 * i1;
		sa1 = sin(arg1); ca1 = -cos(arg1);
		arg2 = dtemp2 + a2 * i1;
		sa2 = -sin(arg2); ca2 = cos(arg2);
		sd = sin(dtempd + ad * i1); sd *= factor / (dblAlpha + sd * sd);
		ii = 2 * i - nx; ii1 = ii + 1;
		speqI1[ii] *= T(sd); speqI1[ii1] *= T(sd);
		speqI2[ii] *= T(sd); speqI2[ii1] *= T(sd);
		temp = T(ca2 * speqI1[ii] + ca1 * speqI2[ii]);
		speqI2[ii] = T(sa2 * speqI1[ii] + sa1 * speqI2[ii]);
		speqI1[ii] = temp;
		temp = T(ca2 * speqI1[ii1] + ca1 * speqI2[ii1]);
		speqI2[ii1] = T(sa2 * speqI1[ii1] + sa1 * speqI2[ii1]);
		speqI1[ii1] = temp;
	}
		
	for (i=0; i<nxd2; i++)
	{
		i1 = nxd2 - i; i1 *= i1;
		ii = i + nxd2;
		dtemp1 = a1 * i1; 
		dtemp2 = a2 * i1; 
		dtempd = ad * i1; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2; j1 *= j1;
			arg1 = dtemp1 + b1 * j1;
			sa1 = sin(arg1); ca1 = -cos(arg1);
			arg2 = dtemp2 + b2 * j1;
			sa2 = -sin(arg2); ca2 = cos(arg2);
			sd = sin(dtempd + bd * j1); sd *= factor / (dblAlpha + sd * sd);
			jj = 2 * j - ny; jj1 = jj + 1;
			I1[ii][jj] *= T(sd); I1[ii][jj1] *= T(sd);
			I2[ii][jj] *= T(sd); I2[ii][jj1] *= T(sd);
			temp = T(ca2 * I1[ii][jj] + ca1 * I2[ii][jj]);
			I2[ii][jj] = T(sa2 * I1[ii][jj] + sa1 * I2[ii][jj]);
			I1[ii][jj] = temp;
			temp = T(ca2 * I1[ii][jj1] + ca1 * I2[ii][jj1]);
			I2[ii][jj1] = T(sa2 * I1[ii][jj1] + sa1 * I2[ii][jj1]);
			I1[ii][jj1] = temp;
		}
	}

	I1[0][0] = I1[0][1] = T(0); // integral of psi is equal to zero
	I2[0][0] = I2[0][1] = T(0); // integral of psi is equal to zero
	for (i=nxd2+1; i<nx; i++)  // cycles for i>nxd2 (i1>0)
	{
		i1 = i - nxd2; i1 *= i1;
		ii = i - nxd2;
		dtemp1 = a1 * i1; 
		dtemp2 = a2 * i1; 
		dtempd = ad * i1; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2; j1 *= j1;
			arg1 = dtemp1 + b1 * j1;
			sa1 = sin(arg1); ca1 = -cos(arg1);
			arg2 = dtemp2 + b2 * j1;
			sa2 = -sin(arg2); ca2 = cos(arg2);
			sd = sin(dtempd + bd * j1); sd *= factor / (dblAlpha + sd * sd);
			jj = 2 * j - ny; jj1 = jj + 1;
			I1[ii][jj] *= T(sd); I1[ii][jj1] *= T(sd);
			I2[ii][jj] *= T(sd); I2[ii][jj1] *= T(sd);
			temp = T(ca2 * I1[ii][jj] + ca1 * I2[ii][jj]);
			I2[ii][jj] = T(sa2 * I1[ii][jj] + sa1 * I2[ii][jj]);
			I1[ii][jj] = temp;
			temp = T(ca2 * I1[ii][jj1] + ca1 * I2[ii][jj1]);
			I2[ii][jj1] = T(sa2 * I1[ii][jj1] + sa1 * I2[ii][jj1]);
			I1[ii][jj1] = temp;
		}
	}

	for (j=nyd2+1; j<ny; j++) // cycle for i==nxd2 (i1==0)
	{
		j1 = j - nyd2; j1 *= j1;
		arg1 = b1 * j1;
		sa1 = sin(arg1); ca1 = -cos(arg1);
		arg2 = b2 * j1;
		sa2 = -sin(arg2); ca2 = cos(arg2);
		sd = sin(bd * j1); sd *= factor / (dblAlpha + sd * sd);
		jj = 2 * j - ny; jj1 = jj + 1;
		I1[0][jj] *= T(sd); I1[0][jj1] *= T(sd);
		I2[0][jj] *= T(sd); I2[0][jj1] *= T(sd);
		temp = T(ca2 * I1[0][jj] + ca1 * I2[0][jj]);
		I2[0][jj] = T(sa2 * I1[0][jj] + sa1 * I2[0][jj]);
		I1[0][jj] = temp;
		temp = T(ca2 * I1[0][jj1] + ca1 * I2[0][jj1]);
		I2[0][jj1] = T(sa2 * I1[0][jj1] + sa1 * I2[0][jj1]);
		I1[0][jj1] = temp;
	}

	// inverse FFT
	fft.Real2D(arrI1, speqI1, nx, ny, xar::OouraFft<T>::eDirInv);
	fft.Real2D(arrI2, speqI2, nx, ny, xar::OouraFft<T>::eDirInv);

	//trim back
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp1(I1);
		tmp1.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
		xar::XArray2DMove<T> tmp2(I2);
		tmp2.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
	}

	// the average value of Im(psi) must be zero
	I1 -= T(I1.Norm(xar::eNormAver));
	// the average value of Re(psi) must be zero
	I2 -= T(I2.Norm(xar::eNormAver));

	// I1 = Im(psi)
	// I2 = A^2 * (1 + 2 * Re(psi))
	I2 *= T(2); 
	I2 += T(1); 
	I2 *= T(Iout);
}

#if(0)
//! Retrieves phase&amplitude from 2 images at different distances from object using iterative 1st Born approximation optimizing over rAlphaOpt
// Solves Fresnel integral equation in 1st Born approximation for a general object
// using the FFT method and optimizing the regularization parameter dblAlpha;
// I1 = image; on exit I1 is replaced by the reconstructed phase (more precisely, by Im(psi0))
// I2 = image; on exit I2 is replaced by the reconstructed intensity (more precisely, by Iout * (1 + 2 * Re(psi)))
// RI1 = propagation distance R' corresponding to image I1
// RI2 = propagation distance R' corresponding to image I2
// sigma = relative standard deviation of the noise corresponding to the average intensity in the image
// rAlphaOpt = reference to return the found optimal value of the regulariztion parameter
//
// NOTE!!!: the program assumes that I1 and I2 are background-corrected, i.e. DIVIDED by the INCIDENT intensity
// NOTE: this program uses the Ooura FFT library
// NOTE: this program also uses modified routine golden() from NumRecipes - see xa_nropt.h
template <class T> void XA_2DBorn<T>::Born2ROpt(
	xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& I2, 
	double RI1, 
	double RI2, 
	double sigma, 
	double& rAlphaOpt)
{
	const IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();

	if (sigma <= 0)
		throw std::invalid_argument("invalid_argument 'sigma' in XA_2DBorn<T>::Born2ROpt (must be positive)");

	if (rAlphaOpt <= 0)
		throw std::invalid_argument("invalid_argument 'rAlphaOpt' in XA_2DBorn<T>::Born2ROpt (must be positive)");

	double Iout = I1.Norm(xar::eNormAver);

	if (fabs(Iout - I2.Norm(xar::eNormAver)) > 0.01 * fabs(Iout))
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Born2ROpt (different averages)");

	m_vTemp.resize(2, 0);
	m_vTemp[0] = &I1; // use m_pTemp[0] to pass the image to CalcOnce
	m_vTemp[1] = &I2; // use m_pTemp[1] to pass the image to CalcOnce

	std::vector<T> vecArgs(7);	// vector of parameters to be passed to 'CalcOnce' via 'goldenC'
	T* mass = &(vecArgs.front());
	mass[0] = T(vecArgs.size() - 1); // number of other parameters in this array
	mass[1] = T(log10(rAlphaOpt));	// base10-log of the regularization parameter
	mass[2] = T(RI1);	// first defocus distance
	mass[3] = T(RI2);	// second defocus distance
	mass[4] = T(sigma * sigma * I1.GetDim1() * I1.GetDim2()); // targeted reconstruction accuracy
	mass[5] = T(std::min(I1.Norm(xar::eNormMin),I2.Norm(xar::eNormMin)));
	mass[6] = T(Iout);

	index_t num_varied = 1;	// logAlpha is the argument (no.1) to be optimized by 'goldenC'
	T tLogAlpha = mass[1];	// a median LogAlpha
	T tLogAlphaMin = mass[1] - 2; // the minimum LogAlpha
	T tLogAlphaMax = mass[1] + 2;	// the maximum LogAlpha
	T tTolerance = T(0.2); // optimization tolerance
	T tAlphaOpt = 0;		// storage for the returned optimum dblAlpha

	T tRecAccuracy = xar::XA_Optimisation::GoldenMinimumSearch< XA::XA_2DBorn<T> >(
		tLogAlphaMin, 
		tLogAlpha, 
		tLogAlphaMax,
		*this, 
		mass, 
		num_varied, 
		tTolerance, 
		&tAlphaOpt, 
		2);

	// do the phase retrieval with the found optimum dblAlpha
	rAlphaOpt = pow(10.0, double(tAlphaOpt));
	Born2R(I1, I2, RI1, RI2, rAlphaOpt);
}
#endif

//! Retrieves phase&amplitude from 2 images at different defocus distances using iterative 1st Born approximation
// Attempts to iteratively improve initial phase guess for a phase-amplitude object
// and solving the equation I1 ~ I1(0) + 2Re{ I0(0) * Psi] }, where U0 = U0(0)(1 + Psi),
// using the Born2R method
//
// I1 = image at z=R1;
// I2 = image at z=R2;
// I0 = initial approximation for the object intensity; on exit I0 is replaced by the corrected intensity
// P0 = initial approximation for the object phase; on exit P0 is replaced by the corrected phase
// R1 = propagation distance R1
// R2 = propagation distance R2
// dblAlpha = regularization parameter
// dblSigma = noise level
// riNumIter = reference to the number of iterations parameter
template <class T> void XA_2DBorn<T>::Born2RIter(
	const xar::XArray2D<T>& I1, 
	const xar::XArray2D<T>& I2, 
	xar::XArray2D<T>& I0, 
	xar::XArray2D<T>& P0, 
	double R1, 
	double R2, 
	double dblAlpha, 
	double dblSigma, 
	index_t& riNumIter)
{
	if (R1 == R2)
		throw std::invalid_argument("invalid_argument 'R1 or R2' in XA_2DBorn<T>::Born2RIter (two propagation distances are equal)");

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Born2RIter (negative regularization parameter)");

	if (dblSigma<0)
		throw std::invalid_argument("invalid_argument 'dblSigma' in XA_2DBorn<T>::Born2RIter (negative noise level)");

	double Iout = I1.Norm(xar::eNormAver);

	if (fabs(Iout - I2.Norm(xar::eNormAver)) > 0.01 * Iout)
		throw std::invalid_argument("invalid_argument 'I2' in XA_2DBorn<T>::Born2RIter (different aver.intensity from I1)");

	if (dblAlpha==0) 
		dblAlpha = 1.e-10;

	bool bNoiseBreak(false);
	double dblNormdP0 = P0.Norm(xar::eNormL2), dblNormTemp;
	double dblSigma2 = dblSigma * dblSigma * I1.GetDim1() * I1.GetDim2();
	double dblImageMin = xar::min(I1.Norm(xar::eNormMin), I2.Norm(xar::eNormMin));
	double dblTolerance = 0.01; // here we always expect the chi-2 to be decreasing
	double norma_resid, norma_resid1, norma_resid2;
	xar::XArray2D<T> A0(I0), dI1, dI2, dP0, dA0;
	xar::XArray2D<std::complex<T> > C0;
	xar::XArray2DFFT<T> CFFTService(C0);

	A0 ^= 0.5; // intensity --> amplitude	

	for (index_t i=0; i<riNumIter; i++)
	{
		MakeComplex(A0, P0, C0, true);
		CFFTService.Kirchhoff(R1);
		Abs2(C0, dI1);
		if (dblSigma > 0.0) // check if the residuals are comparable with the noise
		{
			if (dblImageMin > 0) // no zeros or negative values in the image
				norma_resid = I1.Chi2(dI1, 1.0, true) / dblSigma2 - 1.0;
			else // non-Poisson statistics (sigma is the same at all points)
				norma_resid = I1.Chi2(dI1, 1.0, false) / dblSigma2 - 1.0;
			if (norma_resid < dblTolerance) bNoiseBreak = true;
			else bNoiseBreak = false;
			if (i > 0 && norma_resid > norma_resid1) // last correction was bad 
			{
				P0 -= dP0; // new phase roll back
				A0 /= dA0; // new amplitude roll back
				I0 = A0; I0 ^= 2.0; // amplitude --> intensity	
				riNumIter = i - 1;
				break;
			}
			norma_resid1 = norma_resid;
		}
		
		MakeComplex(A0, P0, C0, true);
		CFFTService.Kirchhoff(R2);
		Abs2(C0, dI2);
		
		if (dblSigma > 0.0) // check if the residuals are comparable with the noise
		{
			if (dblImageMin > 0) // no zeros or negative values in the image
				norma_resid = I2.Chi2(dI2, 1.0, true) / dblSigma2 - 1.0;
			else // non-Poisson statistics (sigma is the same at all points)
				norma_resid = I2.Chi2(dI2, 1.0, false) / dblSigma2 - 1.0;
			if (norma_resid < dblTolerance && bNoiseBreak) 
			{
				riNumIter = i;
				break;
			}
			if (i > 0 && norma_resid > norma_resid2) // last correction was bad 
			{
				P0 -= dP0; // new phase roll back
				A0 /= dA0; // new amplitude roll back
				I0 = A0; I0 ^= 2.0; // amplitude --> intensity	
				riNumIter = i - 1;
				break;
			}
			norma_resid2 = norma_resid;
		}

		dI1 -= I1;
		dI1 *= -T(1.0);
		dI1 += T(1.0);
		dI2 -= I2;
		dI2 *= -T(1.0);
		dI2 += T(1.0);
		Born2R(dI1, dI2, R1, R2, dblAlpha);
		dI1 /= I0; 
		dblNormTemp = dI1.Norm(xar::eNormL2);
		
		if (dblNormTemp >= dblNormdP0) // convergence test
		{
			riNumIter = i;
			break;
		}
		
		dblNormdP0 = dblNormTemp;
		dP0 = dI1; 
		P0 += dI1; // new phase
		dI2 -= T(1.0); 
		dI2 *= T(0.5); 
		dI2 /= I0; // Re[psi0]
		dI2 += T(1.0); 
		dA0 = dI2; 
		A0 *= dI2; //new amplitude
		I0 = A0; 
		I0 ^= 2.0; // amplitude --> intensity
	}
}



//! Retrieves phase&amplitude from 2 images at different defocus distances using Rytov approximation
// Solves Fresnel integral equation in Rytov approximation for an amplitude-phase object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed phase
// I2 = image; on exit I2 is replaced by the reconstructed intensity
// RI1 = propagation distance R' corresponding to image I1
// RI2 = propagation distance R' corresponding to image I2
// dblAlpha (DEFAULT=0) is the regularisation parameter
// NOTE: this program uses the Ooura FFT library
template <class T> void XA_2DBorn<T>::Rytov2R(
	xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& I2, 
	double RI1, 
	double RI2, 
	double dblAlpha)
{

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Rytov2R (negative regularization parameter)");

	if (I1.Norm(xar::eNormMin) <= 0 || I2.Norm(xar::eNormMin) <= 0)
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Rytov2R (non-positive image intensity)");

	double Iout = I1.Norm(xar::eNormAver); // Iout - average transmission coefficient
	
	if (fabs(Iout - I2.Norm(xar::eNormAver)) > 0.01 * fabs(Iout))
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Rytov2R (different averages)");

	// transform I1 into J1 such that J1 / Iout - 1 = log(I1 / Iout), i.e. J1 = Iout * [log(I1 / Iout) + 1]
	I1 /= T(Iout);
	I1.Log();
	I1 += T(1);
	I1 *= T(Iout);
	
	// transform I2 into J2 such that J2 / Iout - 1 = log(I2 / Iout), i.e. J2 = Iout * [log(I2 / Iout) + 1]
	I2 /= T(Iout);
	I2.Log();
	I2 += T(1);
	I2 *= T(Iout);

	Born2R(I1, I2, RI1, RI2, dblAlpha);
}



//! Retrieves phase&amplitude from 2 images at different defocus distances using iterative algorithm by A.Pogany
// Attempts to iteratively improve the initial phase and intensity guesses for a phase-amplitude object
// using the A.Pogany's iterative algorithm (based on the Guigay's formula for the Fresnel propagator),
// and using the Born2R method at each iteration
//
// I1 = image at z=R1;
// I2 = image at z=R2;
// C0 = initial approximation for the object c.amplitude; on exit C0 is replaced by the corrected c.amplitude
// R1 = propagation distance R1
// R2 = propagation distance R2
// dblAlpha = regularization parameter
// riNumIter = number of iterations
// beta = 'under-relaxation' parameter, 0 < beta <= 1.
template <class T> void XA_2DBorn<T>::AP2RIter(
	const xar::XArray2D<T>& I1, 
	const xar::XArray2D<T>& I2, 
	xar::XArray2D<std::complex<T> >& C0, 
	double R1, 
	double R2, 
	double dblAlpha, 
	short niter, 
	double beta)
{
	if (R1 == R2)
		throw std::invalid_argument("invalid_argument 'R1 or R2' in XA_2DBorn<T>::AP2RIter (two propagation distances are equal)");

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::AP2RIter (negative regularization parameter)");

	if (niter<0)
		throw std::invalid_argument("invalid_argument 'niter' in XA_2DBorn<T>::AP2RIter (negative number of iterations parameter)");

	if (beta<=0 || beta>1)
		throw std::invalid_argument("invalid_argument 'beta' in XA_2DBorn<T>::AP2RIter (must be >0 and <= 1)");

	if (dblAlpha==0) dblAlpha = 1.e-10;

	short i;
	xar::XArray2D<T>  I10, I20;
	xar::XArray2D<std::complex<T> > C1;
	xar::XArray2DFFT<T> CFFTService(C1);

	for (i=0; i<niter; i++)
	{
		C1 = C0;
		CFFTService.Kirchhoff(R1);
		Abs2(C1, I10);
		I10 -= I1;
		I10 *= T(-1);
		I10 += T(1);
		
		C1 = C0;
		CFFTService.Kirchhoff(R2);
		Abs2(C1, I20);
		I20 -= I2;
		I20 *= T(-1);
		I20 += T(1);

		Born2R(I10, I20, R1, R2, dblAlpha);
		I20.Log(); 
		I20 *= T(0.5); // intensity --> back to mu*T

		MakeComplex(I20, I10, C1, false);
		C1 *= T(beta);
		C0 += C1;
	}
}



//! Retrieves phase&amplitude from 2 images at different wavelengths (energies) using 1st Born approximation
// Solves Fresnel integral equation in the first Born approximation for an amplitude-phase object
// using the FFT method;
// I1 = image at wl0; on exit I1 is replaced by the reconstructed phase (more precisely, by Im(psi0))
// I2 = image at wl1; on exit I2 is replaced by the reconstructed intensity (more precisely, by Iout * (1 + 2 * Re(psi)))
// R = propagation distance R'
// dblAlpha (DEFAULT=0) is the regularisation parameter
// NOTE: this program uses the Ooura FFT library
//
// NOTE: this algorithms seems to require very close values of wl1 and wl2, and because of that is
// very sensitive to noise. This seems to be in a stark contrast with Born2R, for reasons unknown.
template <class T> void XA_2DBorn<T>::Born2E(
	xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& I2, 
	double R, 
	double dblAlpha)
{
	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Born2E (negative regularization parameter)");

	if (dblAlpha==0) dblAlpha = 1.e-10;

	IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();
	IXAHWave2D* ph22 = GetIXAHWave2D(I2);
	ph22->Validate();

	double wl = ph2->GetWl(); 
	double wl2 = ph22->GetWl(); 
	double xlo = ph2->GetYlo();
	double xhi = ph2->GetYhi();
	double xst = GetYStep(I1);
	double ylo = ph2->GetXlo();
	double yhi = ph2->GetXhi();
	double yst = GetXStep(I1);

	double sigma = wl2 / wl;
	double sigma3 = pow(sigma, 3);
	double Iout = I1.Norm(xar::eNormAver);
	double Iout2 = I2.Norm(xar::eNormAver);
	if (Iout <= 0 || Iout2 <= 0)
		throw std::invalid_argument("invalid_argument 'I1 or I1' in XA_2DBorn<T>::Born2E (non-positive average intensity)");

	if (fabs(log(Iout) - log(Iout2) / sigma3) > 0.1 * fabs(log(Iout)))
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Born2E (inconsistent averages)");

	I1 /= T(Iout);
	I1 -= 1.0;
	I1 *= 0.5;
	I2 /= T(Iout2);
	I2 -= 1.0;
	I2 *= 0.5;
	
	index_t nxF = I1.GetDim1(), nx = nxF;
	index_t nyF = I1.GetDim2(), ny = nyF;

	// pad to the nearest power of 2
	index_t i = 2;
	
	while (i<nx) 
		i *= 2;
	
	nx = i;
	index_t j = 2;
	
	while (j<ny) 
		j *= 2;
	ny = j;
	
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp1(I1);
		tmp1.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
		xar::XArray2DMove<T> tmp2(I2);
		tmp2.Pad((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2, T(0));
	}

	
	T* arrI1 = &(I1.front());
	std::vector<T> vecSpeqI1(nx*2);
	T* speqI1 = &(vecSpeqI1[0]);
	T* arrI2 = &(I2.front());
	std::vector<T> vecSpeqI2(nx*2);
	T* speqI2 = &(vecSpeqI2[0]);

	
	// forward FFTs
	xar::OouraFft<T> fft;

	fft.Real2D(arrI1, speqI1, nx, ny, xar::OouraFft<T>::eDirFwd);
	fft.Real2D(arrI2, speqI2, nx, ny, xar::OouraFft<T>::eDirFwd);

	// multiplication of the FFT by M^(-1)
	// the indexing follows that from the function FFTRe(...) - see vo_fftr2.h
	T temp;
	double dtemp, dtemp1, dtemp2, sa1, sa2, sd, ca1, ca2, arg1, arg2;
	double factor = 2. / (nx * ny);
	double a = xhi - xlo + xst;
	double b = yhi - ylo + yst;
	dtemp = xar::PI * wl;
	dtemp1 = 1.0 / (a * a);
	dtemp2 = 1.0 / (b * b);
	double a1 = dtemp * R * dtemp1;
	double b1 = dtemp * R * dtemp2;
	double a2 = sigma * a1;
	double b2 = sigma * b1;

	index_t nxd2 = nx / 2, i1, ii, ii1;
	index_t nyd2 = ny / 2, j1, jj, jj1;

	dtemp = double(nyd2 * nyd2);
	dtemp1 = b1 * dtemp;
	dtemp2 = b2 * dtemp;
	for (i=0; i<nxd2; i++) 
	{ 
		i1 = nxd2 - i; i1 *= i1;
		arg1 = dtemp1 + a1 * i1;
		sa1 = sin(arg1); ca1 = -cos(arg1);
		arg2 = dtemp2 + a2 * i1;
		sa2 = -sigma * sin(arg2); ca2 = sigma3 * cos(arg2);
		sd = -ca1 * sa2 + ca2 * sa1;
		sd *= factor / (dblAlpha + sd * sd);
		ii = nx + 2 * i; ii1 = ii + 1;
		speqI1[ii] *= T(sd); speqI1[ii1] *= T(sd);
		speqI2[ii] *= T(sd); speqI2[ii1] *= T(sd);
		temp = T(ca2 * speqI1[ii] + ca1 * speqI2[ii]);
		speqI2[ii] = T(sa2 * speqI1[ii] + sa1 * speqI2[ii]);
		speqI1[ii] = temp;
		temp = T(ca2 * speqI1[ii1] + ca1 * speqI2[ii1]);
		speqI2[ii1] = T(sa2 * speqI1[ii1] + sa1 * speqI2[ii1]);
		speqI1[ii1] = temp;
	}

	for (i=nxd2; i<nx; i++) 
	{ 
		i1 = i - nxd2; i1 *= i1;
		arg1 = dtemp1 + a1 * i1;
		sa1 = sin(arg1); ca1 = -cos(arg1);
		arg2 = dtemp2 + a2 * i1;
		sa2 = -sigma * sin(arg2); ca2 = sigma3 * cos(arg2);
		sd = -ca1 * sa2 + ca2 * sa1;
		sd *= factor / (dblAlpha + sd * sd);
		ii = 2 * i - nx; ii1 = ii + 1;
		speqI1[ii] *= T(sd); speqI1[ii1] *= T(sd);
		speqI2[ii] *= T(sd); speqI2[ii1] *= T(sd);
		temp = T(ca2 * speqI1[ii] + ca1 * speqI2[ii]);
		speqI2[ii] = T(sa2 * speqI1[ii] + sa1 * speqI2[ii]);
		speqI1[ii] = temp;
		temp = T(ca2 * speqI1[ii1] + ca1 * speqI2[ii1]);
		speqI2[ii1] = T(sa2 * speqI1[ii1] + sa1 * speqI2[ii1]);
		speqI1[ii1] = temp;
	}
		
	for (i=0; i<nxd2; i++)
	{
		i1 = nxd2 - i; i1 *= i1;
		ii = i + nxd2;
		dtemp1 = a1 * i1; 
		dtemp2 = a2 * i1; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2; j1 *= j1;
			arg1 = dtemp1 + b1 * j1;
			sa1 = sin(arg1); ca1 = -cos(arg1);
			arg2 = dtemp2 + b2 * j1;
			sa2 = -sigma * sin(arg2); ca2 = sigma3 * cos(arg2);
			sd = -ca1 * sa2 + ca2 * sa1; 
			sd *= factor / (dblAlpha + sd * sd);
			jj = 2 * j - ny; jj1 = jj + 1;
			I1[ii][jj] *= T(sd); I1[ii][jj1] *= T(sd);
			I2[ii][jj] *= T(sd); I2[ii][jj1] *= T(sd);
			temp = T(ca2 * I1[ii][jj] + ca1 * I2[ii][jj]);
			I2[ii][jj] = T(sa2 * I1[ii][jj] + sa1 * I2[ii][jj]);
			I1[ii][jj] = temp;
			temp = T(ca2 * I1[ii][jj1] + ca1 * I2[ii][jj1]);
			I2[ii][jj1] = T(sa2 * I1[ii][jj1] + sa1 * I2[ii][jj1]);
			I1[ii][jj1] = temp;
		}
	}

	I1[0][0] = I1[0][1] = T(0); // integral of psi is equal to zero
	I2[0][0] = I2[0][1] = T(0); // integral of psi is equal to zero
	for (i=nxd2+1; i<nx; i++)  //cycles for i>nxd2 (i1>0)
	{
		i1 = i - nxd2; i1 *= i1;
		ii = i - nxd2;
		dtemp1 = a1 * i1; 
		dtemp2 = a2 * i1; 
		for (j=nyd2; j<ny; j++)
		{
			j1 = j - nyd2; j1 *= j1;
			arg1 = dtemp1 + b1 * j1;
			sa1 = sin(arg1); ca1 = -cos(arg1);
			arg2 = dtemp2 + b2 * j1;
			sa2 = -sigma * sin(arg2); ca2 = sigma3 * cos(arg2);
			sd = -ca1 * sa2 + ca2 * sa1;
			sd *= factor / (dblAlpha + sd * sd);
			jj = 2 * j - ny; jj1 = jj + 1;
			I1[ii][jj] *= T(sd); I1[ii][jj1] *= T(sd);
			I2[ii][jj] *= T(sd); I2[ii][jj1] *= T(sd);
			temp = T(ca2 * I1[ii][jj] + ca1 * I2[ii][jj]);
			I2[ii][jj] = T(sa2 * I1[ii][jj] + sa1 * I2[ii][jj]);
			I1[ii][jj] = temp;
			temp = T(ca2 * I1[ii][jj1] + ca1 * I2[ii][jj1]);
			I2[ii][jj1] = T(sa2 * I1[ii][jj1] + sa1 * I2[ii][jj1]);
			I1[ii][jj1] = temp;
		}
	}

	for (j=nyd2+1; j<ny; j++) //cycle for i==nxd2 (i1==0)
	{
		j1 = j - nyd2; j1 *= j1;
		arg1 = b1 * j1;
		sa1 = sin(arg1); ca1 = -cos(arg1);
		arg2 = b2 * j1;
		sa2 = -sigma * sin(arg2); ca2 = sigma3 * cos(arg2);
		sd = -ca1 * sa2 + ca2 * sa1;
		sd *= factor / (dblAlpha + sd * sd);
		jj = 2 * j - ny; jj1 = jj + 1;
		I1[0][jj] *= T(sd); I1[0][jj1] *= T(sd);
		I2[0][jj] *= T(sd); I2[0][jj1] *= T(sd);
		temp = T(ca2 * I1[0][jj] + ca1 * I2[0][jj]);
		I2[0][jj] = T(sa2 * I1[0][jj] + sa1 * I2[0][jj]);
		I1[0][jj] = temp;
		temp = T(ca2 * I1[0][jj1] + ca1 * I2[0][jj1]);
		I2[0][jj1] = T(sa2 * I1[0][jj1] + sa1 * I2[0][jj1]);
		I1[0][jj1] = temp;
	}

	// inverse FFT
	fft.Real2D(arrI1, speqI1, nx, ny, xar::OouraFft<T>::eDirInv);
	fft.Real2D(arrI2, speqI2, nx, ny, xar::OouraFft<T>::eDirInv);

	// trim back
	if (nx != nxF || ny != nyF)
	{
		xar::XArray2DMove<T> tmp1(I1);
		tmp1.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
		xar::XArray2DMove<T> tmp2(I2);
		tmp2.Trim((nx-nxF)/2, nx-nxF-(nx-nxF)/2, (ny-nyF)/2, ny-nyF-(ny-nyF)/2);
	}

	// the average value of Im(psi) must be zero
	I1 -= T(I1.Norm(xar::eNormAver));
	// the average value of Re(psi) must be zero
	I2 -= T(I2.Norm(xar::eNormAver));

	// I1 = Im(psi)
	// I2 = A^2 * (1 + 2 * Re(psi))
	I2 *= T(2); 
	I2 += T(1); 
	I2 *= T(Iout);
	I2.SetHeadPtr(I1.GetHeadPtr()->Clone());
}



//! Retrieves phase&amplitude from 2 images at different wavelengths (energies) using Rytov approximation
// Solves Fresnel integral equation in Rytov approximation for an amplitude-phase object
// using the FFT method;
// I1 = image; on exit I1 is replaced by the reconstructed phase
// I2 = image; on exit I2 is replaced by the reconstructed intensity
// R = propagation distance R'
// dblAlpha (DEFAULT=0) is the regularisation parameter
// NOTE: this program uses the Ooura FFT library
//
// NOTE: this algorithms seems to require very close values of wl1 and wl2, and because of that is
// very sensitive to noise. This seems to be in a stark contrast with Rytov2R, for reasons unknown.
template <class T> void XA_2DBorn<T>::Rytov2E(
	xar::XArray2D<T>& I1, 
	xar::XArray2D<T>& I2, 
	double R, 
	double dblAlpha)
{

	if (dblAlpha<0)
		throw std::invalid_argument("invalid_argument 'dblAlpha' in XA_2DBorn<T>::Rytov2E (negative regularization parameter)");

	IXAHWave2D* ph2 = GetIXAHWave2D(I1);
	ph2->Validate();
	IXAHWave2D* ph22 = GetIXAHWave2D(I2);
	ph22->Validate();

	double wl = ph2->GetWl(); 
	double wl2 = ph22->GetWl(); 
	double sigma = wl2 / wl;
	double sigma3 = pow(sigma, 3);
	double Iout = I1.Norm(xar::eNormAver);
	double Iout2 = I2.Norm(xar::eNormAver);
	
	if (I1.Norm(xar::eNormMin) <= 0 || I2.Norm(xar::eNormMin) <= 0)
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Rytov2E (non-positive intensity)");
	
	if (fabs(log(Iout) - log(Iout2) / sigma3) > 0.1 * fabs(log(Iout)))
		throw std::invalid_argument("invalid_argument 'I1 or I2' in XA_2DBorn<T>::Rytov2E (inconsistent averages)");

	// transform I1 into J1 such that J1 / Iout - 1 = log(I1 / Iout), i.e. J1 = Iout * [log(I1 / Iout) + 1]
	I1 /= T(Iout);
	I1.Log();
	I1 += T(1);
	I1 *= T(Iout);
	// transform I2 into J2 such that J2 / Iout - 1 = log(I2 / Iout), i.e. J2 = Iout * [log(I2 / Iout) + 1]
	I2 /= T(Iout2);
	I2.Log();
	I2 += T(1);
	I2 *= T(Iout2);

	Born2E(I1, I2, R, dblAlpha);
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
	template class XA::XA_2DBorn<float>;
	template class XA::XA_2DBorn<double>;
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
