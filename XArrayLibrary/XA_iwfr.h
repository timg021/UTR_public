//Header XA_iwfr.h
//
//
//	HEADER FILE TITLE:
//
//		Phase retrieval algorithms based on the iwfr (Iterative Wave Function Reconstruction) method
//
/*!
	\file		XA_iwfr.h
	\brief		Phase retrieval algorithms based on the iwfr (Iterative Wave Function Reconstruction) method
	\par		Description:
		This header contains a class that provides phase retrieval services 
		based on the iwfr (Iterative Wave Function Reconstruction) method for XArray2D<T> objects
*/

#pragma once

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_spln2.h"
#include "XA_fft2.h"
#include "fftwd2c.h"
#include "fftwd2fc.h"
#include <omp.h>

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
//Class XA_IWFR<T>
//
//	Phase retrieval algorithms based on the 1st iwfr and Rytov approximations
//
/*!
	\brief		Phase retrieval algorithms based on the iwfr (Iterative Wave Function Reconstruction) method
	\par		Description:
				This class template provides iwfr (Iterative Wave Function Reconstruction) 
				services for XArray2D<T> objects
	\remarks	An object of this class represents an interface exposing several functions
				that provide iwfr (Iterative Wave Function Reconstruction) phase retrieval services				
	\remarks	This class can only be created for T=float or T=double
	\warning	Copying of objects of this class does not make sense and is prohibited
*/
	template <class T> class XA_IWFR
	{
	// Typedefs
	public:
		typedef T type;
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XA_IWFR() {}
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XA_IWFR(const XA_IWFR<T>& rCopy) { GetValuetype(); }
	public:
		//! Destructor
		~XA_IWFR() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XA_IWFR<T>& rCopy) { if (this == &rCopy) return; }

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XA_IWFR<T> objects for unsupported types T
		//! Returns the xar::_eValueType corresponding to T
		static xar::_eValueType GetValuetype(void);

	// Operations
	public:
		//! Retrieves complex wave function from a vector of defocused images
		void Iwfr(vector<XArray2D<T> >& vint0, XArray2D<std::complex<T> >& campOut, vector<double> vdefocusdist, double zmiddle, double k2maxo, double Cs3, double Cs5, int kmax, double epsilon, bool bVerboseOutput);
		//! Retrieves phase from several defocused images using L2 minimization of the CTF function
		XArray2D<T> CTFL2(vector<XArray2D<T> >& vint0, vector<double> vdefocusdist, double q2max, double Cs3, double Cs5, double alpha);
		//! Returns the conjugated phase retrieved from a single defocused image: -phase = -0.5 * NF * (Intensity / I0 - 1)
		XArray2D<T> Min05LogAmp(const XArray2D<T>& intIn, double NF = 1.0);
		//! Returns the conjugated phase retrieved from a single defocused image: -phase = -NF^(-1) + sqrt[NF^(-2) + (Intensity / I0) - 1]
		XArray2D<T> ConjPhaseGausBeam(const XArray2D<T>& intIn, double NF = 1.0);
		//! Inverts CTF for a pure phase or homogeneous object
		void InvertPhaseCTF(XArray2D<T>& int0, double defocusdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha);
		//! Calculates CTF-based free-space propagation for a pure phase object
		void ForwardPhaseCTF(XArray2D<T>& pha0, double defocusdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5);
		//! Inverts symmetrised CTF on the Ewald sphere for a pure phase object
		void InvertCTF_DT(const XArray2D<T>& LnIdIin, XArray2D<std::complex<T> >& cout, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC, int NumThreads);
		//! Inverts symmetrised CTF on the Ewald sphere for a pure phase object using FFTW library
		void InvertCTF_DT1(const XArray2D<T>& LnIdIin, XArray2D<std::complex<T> >& camp, Fftwd2c& fftw2D, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC, int NumThreads);
		//! Inverts symmetrised CTF on the Ewald sphere for a pure phase object using FFTWf (single-precision) library
		void InvertCTF_DT2(const XArray2D<T>& LnIdIin, XArray2D<std::complex<T> >& camp, Fftwd2fc& xafft, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC, int NumThreads);

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
//! Returns the xar::_eValueType corresponding to T=float
template<> inline xar::_eValueType XA_IWFR<float>::GetValuetype() { return xar::eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> inline xar::_eValueType XA_IWFR<double>::GetValuetype() { return xar::eXADouble; }

//! Retrieves phase from several defocused images using Iterative Wave Function Reconstruction method
// vint0 - input vector of defocused images
// campOut - output complex amplitude in the plane zmiddle
// vdefocusdist - vector of defocus distances
// zmiddle - z-coordinate of the defocus plane where the output is produced
// k2maxo - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// kmax - maximum number of iterations
// epsilon - minimum L2 error (interrupts the iterations)
// bVerboseOutput - if TRUE, then verbose diagnostic output is printed
template <class T> void XA_IWFR<T>::Iwfr(vector< XArray2D<T> >& vint0, XArray2D<std::complex<T> >& campOut, vector<double> vdefocusdist, double zmiddle, double k2maxo, double Cs3, double Cs5, int kmax, double epsilon, bool bVerboseOutput)
{
	int ndefocus = (int)vint0.size();
	if (ndefocus < 1) throw std::runtime_error("Error: input image vector empty in IWFR()");

	bool bAbort(false);
	double dtemp; // auxilliary variable
	double ssej(0.0), ssejm1(0.0); // current and previous average reconstruction errors
	vector<double> verr(ndefocus); // reconstruction errors in individual defocus planes
	vector<double> vint0_L1(ndefocus); // L1 norms of the initial defocused intensities
	vector<XArray2D<T> > vint(ndefocus); // iterated defocused intensities
	vector<XArray2D<std::complex<T> > > vcamp(ndefocus); // defocused complex amplitudes

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(vint0[0].GetHeadPtr()->Clone());
	IXAHWave2D* ph2 = GetIXAHWave2D(vint0[0]);
	ph2->Validate();
	index_t ny = vint0[0].GetDim1();
	index_t nx = vint0[0].GetDim2();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);

	vector<int> iSign(ndefocus);
	for (int n = 0; n < ndefocus; n++) iSign[n] = vdefocusdist[n] < zmiddle ? iSign[n] = -1 : iSign[n] = 1;

	for (int k = 0; k < kmax; k++)
	{
		// apply initial or newly reconstructed phases (the second case is equal to restoring the original moduli)
		// and propagate each defocused amplitude to the "middle" plane z = zmiddle
		#pragma omp parallel for
		for (int n = 0; n < ndefocus; n++)
		{
			try
			{
				if (k == 0) // create initial defocused complex amplitude
				{
					vint0_L1[n] = vint0[n].Norm(eNormL1);
					if (bVerboseOutput) printf("\nL1 norm of input defocused intensity no. %d = %g", n, vint0_L1[n]);
					if (vint0_L1[n] == 0) throw std::runtime_error("Error: input intensity file is empty in IWFR()");
					vint0[n] ^= 0.5; // intensity --> real amplitude
					MakeComplex(vint0[n], 0.0, vcamp[n], true); // apply initial zero phases
				}
				else // apply phases obtained on the previous iteration
					for (index_t j = 0; j < vcamp[n].GetDim1(); j++)
						for (index_t i = 0; i < vcamp[n].GetDim2(); i++)
						{
							dtemp = abs(vcamp[n][j][i]);
							if (dtemp) vcamp[n][j][i] *= vint0[n][j][i] / dtemp;
							else vcamp[n][j][i] = std::polar(vint0[n][j][i], 0.0);
						}
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(zmiddle - vdefocusdist[n], false, k2maxo, -iSign[n] * Cs3, -iSign[n] * Cs5); // propagate to z = 0
			}
			catch (std::runtime_error& E)
			{
				printf("\n\n!!!Exception: %s\n", E.what());
				bAbort = true;
			}
		}
		if (bAbort) throw std::runtime_error("at least one thread has thrown an exception in IWFR().");

		// average complex amplitudes in the middle plane
		campOut = vcamp[0];
		for (index_t n = 1; n < ndefocus; n++) campOut += vcamp[n];
		campOut /= double(ndefocus);

		// propagate the averaged complex amplitude from the middle plane to the individual defocus planes
		#pragma omp parallel for shared(campOut)
		for (int n = 0; n < ndefocus; n++)
		{
			try
			{
				vcamp[n] = campOut;
				xar::XArray2DFFT<double> xafft(vcamp[n]);
				xafft.Fresnel(vdefocusdist[n] - zmiddle, false, k2maxo, iSign[n] * Cs3, iSign[n] * Cs5); // propagate to z = z[n]
				Abs(vcamp[n], vint[n]);
				vint[n] -= vint0[n];
				verr[n] = pow(vint[n].Norm(eNormL2), 2.0) / vint0_L1[n];
			}
			catch (std::runtime_error& E)
			{
				printf("\n\n!!!Exception in IWFR(): %s\n", E.what());
				bAbort = true;
			}
		}
		if (bAbort) throw std::runtime_error("at least one thread has thrown an exception in IWFR().");

		// calculate the current reconstruction error and
		// if the difference with the previous error is smaller than the defined minimum
		// or, if the error started to increase, interrupt the iterations
		ssej = 0.0;
		for (index_t n = 0; n < ndefocus; n++) ssej += verr[n];
		ssej /= double(ndefocus);

		if (bVerboseOutput)
		{
			if (k == 0) printf("\nIteration number %d; SSE_aver error at 0th iteration = %g\n", k, ssej);
			else printf("\nIteration number %d; SSE_aver error difference with previous iteration = %g\n", k, ssejm1 - ssej);

			for (index_t n = 0; n < ndefocus; n++) printf("SSE(%zd) = %g ", n, verr[n]);
		}

		if (k > 0 && (ssejm1 - ssej) < epsilon) break;
		else ssejm1 = ssej;
	}
}

//! Retrieves phase from several defocused images using L2 minimization of the CTF function
// vint0 - input vector of defocused images
// vdefocusdist - vector of defocus distances
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// alpha - Tikhonov regularization parameter (as in the denominator of eq.(A6) of paper D. PAGANIN et al, J.Micros. 214 (2004) 51-61)
// returns output phase in the plane vdefocusdist[0]
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
// NOTE: it could be useful later to implement a modified version of this algorithm, where the phase is not simply retrieved at defocus[0], which is highly dependent on the image[0], 
// but instead it would be retrieved at each defocus[n], and then the complex amplitude would be averaged, e.g. at the middle plane, as in IWFR
template <class T> XArray2D<T> XA_IWFR<T>::CTFL2(vector<XArray2D<T> >& vint0, vector<double> vdefocusdist, double q2max, double Cs3, double Cs5, double alpha)
{
	bool bAper(q2max > 0);
	int ndefocus = (int)vint0.size();
	if (ndefocus < 2) throw std::runtime_error("Error: at least two defocused images are required in CTFL2()");

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(vint0[0].GetHeadPtr()->Clone());
	IXAHWave2D* ph2 = GetIXAHWave2D(vint0[0]);
	ph2->Validate();
	index_t ny = vint0[0].GetDim1();
	index_t nx = vint0[0].GetDim2();
	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	// convert input defocused images to contrast functions and complexify for FFT (should change the code to use real FFT at a later time)
	vector<XArray2D<std::complex<T> > > vcamp(ndefocus);
	for (int n = 0; n < ndefocus; n++)
	{
		vint0[n] -= vint0[n].Norm(eNormAver); // image --> contrast function * (-1)
		MakeComplex(vint0[n], T(0), vcamp[n], false);
	}

	index_t i = 2;
	while (i < ny) i *= 2;
	if (i != ny) throw std::invalid_argument("invalid_argument in CTFL2() (m_dim1 is not a power of 2)");
	index_t j = 2;
	while (j < nx) j *= 2;
	if (j != nx) throw std::invalid_argument("invalid_argument in CTFL2() (m_dim2 is not a power of 2)");

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::CTFL2 (evanescent waves present)");

	vector<int> iSign(ndefocus);
	vector<double> dblDistance(ndefocus);
	dblDistance[0] = 0.0;
	for (int n = 1; n < ndefocus; n++)
	{
		dblDistance[n] = vdefocusdist[n] - vdefocusdist[0];
		dblDistance[n] < 0 ? iSign[n] = -1 : iSign[n] = 1;
		if (dblDistance[n] == 0)
			throw std::runtime_error("runtime_error in XA_IWFR<T>::CTFL2 (duplicate propagation distance)");
	}

	//********* Fourier transforming initial amplitude u_n(i,j)
	vector<T*> u(ndefocus);
	OouraFft<T> fft;
	for (int n = 0; n < ndefocus; n++)
	{
		u[n] = reinterpret_cast<T*>(&(vcamp[n].front()));
		fft.Complex2D((std::complex<T> *) u[n], ny, nx, OouraFft<T>::eDirFwd);
	}

	// prepare output phase array
	XArray2D<T> fiOut(ny, nx, T(0));
	fiOut.SetHeadPtr(vint0[0].GetHeadPtr()->Clone());
	XArray2D<std::complex<T> > vcampOut{ MakeComplex(fiOut, T(0), false) };
	T* fiC;
	fiC = reinterpret_cast<T*>(&(vcampOut.front()));
	fft.Complex2D((std::complex<T> *) fiC, ny, nx, OouraFft<T>::eDirFwd);

	//********* Dividing F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi2 = 1.0 / xap2;
	double deta2 = 1.0 / yap2;

	vector<double> fac2(ndefocus);
	double fac2min = PI * wl * abs(dblDistance[1]);
	for (int n = 1; n < ndefocus; n++)
	{
		fac2[n] = PI * wl * dblDistance[n];
		if (abs(fac2[n]) < fac2min) fac2min = abs(fac2[n]);
	}

	if (alpha <= 0) alpha = 0.1 * pow(fac2min * std::min(dcsi2, deta2), 4);

	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	index_t k, kj;
	double eta2, q2;
	double dtemp, sintemp, sumsin2, costemp, sumuk, sumuk1, Cstemp(0), Cstemp0(0);

	for (long i = -long(nyd2); i < 0; i++)
	{
		kj = nxy2 + nx2 * i + nx2;
		eta2 = deta2 * i * i;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp0 = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					if (bC35) Cstemp = iSign[n] * Cstemp0;
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
		kj = nxy2 + nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp0 = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					if (bC35) Cstemp = iSign[n] * Cstemp0;
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
	}
	for (long i = 0; i < long(nyd2); i++)
	{
		kj = nx2 * i + nx2;
		eta2 = deta2 * i * i;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp0 = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					if (bC35) Cstemp = iSign[n] * Cstemp0;
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
		kj = nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			q2 = dcsi2 * j * j + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp0 = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				sumuk = sumuk1 = sumsin2 = 0;
				for (int n = 1; n < ndefocus; n++)
				{
					if (bC35) Cstemp = iSign[n] * Cstemp0;
					dtemp = fac2[n] * q2 + Cstemp;
					sintemp = sin(dtemp);
					costemp = cos(dtemp);
					sumsin2 += sintemp * sintemp;
					sumuk += 0.5 * (u[n][k] - u[0][k] * costemp) * sintemp;
					sumuk1 += 0.5 * (u[n][k + 1] - u[0][k + 1] * costemp) * sintemp;
				}
				dtemp = sumsin2 / (sumsin2 * sumsin2 + alpha);
				fiC[k] = T(sumuk * dtemp);
				fiC[k + 1] = T(sumuk1 * dtemp);
			}
			else
			{
				fiC[k] = T(0);
				fiC[k + 1] = T(0);
			}
		}
	}

	//********* inverse Fourier transforming
	fft.Complex2D((std::complex<T> *) fiC, ny, nx, OouraFft<T>::eDirInv);
	T fact = T(1.0) / nxy;
	for (k = 0; k < nxy2; k++)	fiC[k] *= fact;
	return Re(vcampOut); // relies on move assignment / constructor to avoid copying this object on return
}


//! Returns the conjugated phase retrieved from a single defocused image, i.e. returns -0.5 * NF * (Intensity / I0 - 1)
// intIn - input defocused image
// NF - Fresnel number NF = lambda * z / (2 * PI * sigma^2) 
// returns minus output phase (i.e. conjugated phase) in the same plane
// NOTE: here good results are obtained with "minus" log and positive d
// Typical value of dMult for light atoms (such as Oxygen) is 0.5, typical value for jheavy atoms (such as Au) is 0.1
template <class T> XArray2D<T> XA_IWFR<T>::Min05LogAmp(const XArray2D<T>& intIn, double NF)
{
	const double* pInt = &(intIn[0][0]);
	XArray2D<T> fiOut{ intIn };
	double* pPha = &(fiOut[0][0]);
	NF *= -0.5;
	for (index_t i = 0; i < fiOut.size(); i++)
		pPha[i] = NF * (pInt[i] - 1.0);

	return fiOut; // relies on move assignment / constructor to avoid copying this object on return
}


//! Returns the conjugated phase retrieved from a single defocused image, i.e. returns NF^(-1) - sqrt[NF^(-2) + (Intensity / I0) - 1]
// intIn - input defocused image
// NF - Fresnel number NF = lambda * z / (2 * PI * sigma^2) 
// returns minus output phase (i.e. conjugated phase) in the same plane
template <class T> XArray2D<T> XA_IWFR<T>::ConjPhaseGausBeam(const XArray2D<T>& intIn, double NF)
{
	XArray2D<T> fiOut{ intIn };
	if (NF == 0)
	{
		fiOut.Fill(T(0));
	}
	else
	{
		double invNF = 1.0 / NF;
		double invNF2 = invNF * invNF - 1.0;
		const double* pInt = &(intIn[0][0]);
		double* pPha = &(fiOut[0][0]);

		for (index_t i = 0; i < fiOut.size(); i++)
			pPha[i] = invNF - sqrt(abs(invNF2 + pInt[i]));
	}

	return fiOut; // relies on move assignment / constructor to avoid copying this object on return
}


//! Inverts CTF for a pure phase or homogeneous object
// int0 - input defocused image (it is "spoiled" inside this function)
// defocdist - defocus distance
// Z1mZ2d2 astigmatism parameter(dblDistanceX - dblDistanceY) / 2 (in the same units as used in the Wavehead2D)
// phiA astigmatism angle(with the X axis, counterclockwise) (in radians)
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// sigma - beta to delta ratio (1 / gamma), it is equal to zero for pure phase objects
// alpha - Tikhonov regularization parameter
// On return, the retrieved output phase in the plane -defocdist from the image plane replaces the input image
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
template <class T> void XA_IWFR<T>::InvertPhaseCTF(XArray2D<T>& int0, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha)
{
	//if (defocdist == 0) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertPhaseCTF() (propagation distance cannot be zero in this phase-contrast method)");

	bool bAper(q2max > 0);

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(int0.GetHeadPtr()->Clone());
	IXAHWave2D* ph2 = GetIXAHWave2D(int0);
	ph2->Validate();
	index_t ny = int0.GetDim1();
	index_t nx = int0.GetDim2();

	index_t i = 2;
	while (i < ny) i *= 2;
	if (i != ny) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertPhaseCTF() (m_dim1 is not a power of 2)");
	index_t j = 2;
	while (j < nx) j *= 2;
	if (j != nx) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertPhaseCTF() (m_dim2 is not a power of 2)");

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap = abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	double omega = atan(sigma);
	double cosomega = 1.0 / sqrt(1.0 + sigma * sigma);

	// I --> (1 - I / I_in) / 2, and then complexify for FFT (should change the code to use real FFT at a later time!)
	T aver = T(int0.Norm(eNormAver));
	if (aver <= 0) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertPhaseCTF() (average value in the input image is not positive)");
	XArray2D<std::complex<T> > camp;
	//int0 *= -0.5 / (aver * exp(2.0 * sigma)); // !!! this theoretically correct normalization seems incompatible with the weak object approximation, leading to poor results
	int0 *= -0.5 / aver;
	int0 += 0.5;
	MakeComplex(int0, T(0), camp, false);

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::InvertPhaseCTF() (evanescent waves present)");

	//********* Fourier transforming initial amplitude u_n(i,j)
	T* u(0);
	OouraFft<T> fft;
	u = reinterpret_cast<T*>(&(camp.front()));
	fft.Complex2D((std::complex<T> *) u, ny, nx, OouraFft<T>::eDirFwd);

	//********* Divide F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi = 1.0 / xap;
	double dcsi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;

	double fac2 = PI * wl * defocdist;
	if (alpha <= 0) alpha = 0.1 * pow(abs(fac2) * std::min(dcsi2, deta2), 2);

	double Z1 = defocdist + Z1mZ2d2;
	double Z2 = defocdist - Z1mZ2d2;
	double cosphiA = cos(phiA);
	double sinphiA = sin(phiA);
	double sin2phiA = 2.0 * sinphiA * cosphiA;
	double fac2x = PI * wl * (Z1 * cosphiA * cosphiA + Z2 * sinphiA * sinphiA);
	double fac2y = PI * wl * (Z2 * cosphiA * cosphiA + Z1 * sinphiA * sinphiA);
	double facxy = PI * wl * Z1mZ2d2 * sin2phiA;
	facxy = facxy * deta * dcsi;
	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	index_t k, kj;
	double eta2, csi2, q2;
	double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, Cstemp(0);
	T temp;

	for (long i = -long(nyd2); i < 0; i++)
	{
		kj = nxy2 + nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y + omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				temp = T(sintemp / (sin2 + alpha));
				u[k] *= temp;
				u[k + 1] *= temp;
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
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				temp = T(sintemp / (sin2 + alpha));
				u[k] *= temp;
				u[k + 1] *= temp;
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
		eta2fac2y = eta2 * fac2y + omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				temp = T(sintemp / (sin2 + alpha));
				u[k] *= temp;
				u[k + 1] *= temp;
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
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				temp = T(sintemp / (sin2 + alpha));
				u[k] *= temp;
				u[k + 1] *= temp;
			}
			else
			{
				u[k] = T(0);
				u[k + 1] = T(0);
			}
		}
	}

	//********* inverse Fourier transforming
	fft.Complex2D((std::complex<T> *) u, ny, nx, OouraFft<T>::eDirInv);
	T fact = T(cosomega / nxy);
	for (k = 0; k < nxy2; k++)	u[k] *= fact;
	int0 = Re(camp);
}

//! Calculates CTF-based free-space propagation for a pure phase object
// pha0 - input defocused image
// defocdist - defocus distance
// Z1mZ2d2 astigmatism parameter(dblDistanceX - dblDistanceY) / 2 (in the same units as used in the Wavehead2D)
// phiA astigmatism angle(with the X axis, counterclockwise) (in radians)
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// The propagated intensity in the plane defocdist from the object plane replaces the input phase
//
template <class T> void XA_IWFR<T>::ForwardPhaseCTF(XArray2D<T>& pha0, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5)
{
	bool bAper(q2max > 0);

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(pha0.GetHeadPtr()->Clone());
	IXAHWave2D* ph2 = GetIXAHWave2D(pha0);
	ph2->Validate();
	index_t ny = pha0.GetDim1();
	index_t nx = pha0.GetDim2();

	index_t i = 2;
	while (i < ny) i *= 2;
	if (i != ny) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::ForwardPhaseCTF() (m_dim1 is not a power of 2)");
	index_t j = 2;
	while (j < nx) j *= 2;
	if (j != nx) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::ForwardPhaseCTF() (m_dim2 is not a power of 2)");

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap = abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	// complexify input phase for FFT (should change the code to use real FFT at a later time!)
	XArray2D<std::complex<T> > camp;
	MakeComplex(pha0, T(0), camp, false);

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::ForwardPhaseCTF() (evanescent waves present)");

	//********* Fourier transforming initial amplitude u_n(i,j)
	T* u(0);
	OouraFft<T> fft;
	u = reinterpret_cast<T*>(&(camp.front()));
	fft.Complex2D((std::complex<T> *) u, ny, nx, OouraFft<T>::eDirFwd);

	//********* Multiply F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi = 1.0 / xap;
	double dcsi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;

	double fac2 = PI * wl * defocdist;
	double Z1 = defocdist + Z1mZ2d2;
	double Z2 = defocdist - Z1mZ2d2;
	double cosphiA = cos(phiA);
	double sinphiA = sin(phiA);
	double sin2phiA = 2.0 * sinphiA * cosphiA;
	double fac2x = PI * wl * (Z1 * cosphiA * cosphiA + Z2 * sinphiA * sinphiA);
	double fac2y = PI * wl * (Z2 * cosphiA * cosphiA + Z1 * sinphiA * sinphiA);
	double facxy = PI * wl * Z1mZ2d2 * sin2phiA;
	facxy = facxy * deta * dcsi;
	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	index_t k, kj;
	double eta2, csi2, q2;
	double etafacxy, eta2fac2y, fac2a, Cstemp(0);
	T temp;

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
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				temp = (T)sin(fac2a + Cstemp);
				u[k] *= temp;
				u[k + 1] *= temp;
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
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				temp = (T)sin(fac2a + Cstemp);
				u[k] *= temp;
				u[k + 1] *= temp;
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
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				temp = (T)sin(fac2a + Cstemp);
				u[k] *= temp;
				u[k + 1] *= temp;
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
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				temp = (T)sin(fac2a + Cstemp);
				u[k] *= temp;
				u[k + 1] *= temp;
			}
			else
			{
				u[k] = T(0);
				u[k + 1] = T(0);
			}
		}
	}

	//********* inverse Fourier transforming
	fft.Complex2D((std::complex<T> *) u, ny, nx, OouraFft<T>::eDirInv);
	T fact = T(1.0) / nxy;
	for (k = 0; k < nxy2; k++)	u[k] *= fact;
	pha0 = Re(camp);

	//********* finally, multiply by 2 and add 1 (assuming that I_in = 1)
	pha0 *= T(2.0);
	pha0 += T(1.0);
}

//! Inverts symmetrised CTF on the Ewald sphere for a monomoprphous object using Ooura FFT library
// LnIdIin - ln(I/I_in), where I is the defocused image
// cout - output solution in the form of 3D FFT of  delta * [4 * PI * sqrt(1 + sigma^2)] / wl, on the Ewald sphere in the Fourier space, in the "unshuffled form"
// defocdist - defocus distance
// Z1mZ2d2 astigmatism parameter(dblDistanceX - dblDistanceY) / 2 (in the same units as used in the Wavehead2D)
// phiA astigmatism angle(with the X axis, counterclockwise) (in radians)
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// sigma - beta to delta ratio (1 / gamma), it is equal to zero for pure phase objects
// alpha - Tikhonov regularization parameter
// bESCC - if true - Ewald sphere curvature correction is applied, if false - flat Ewald sphere is simulated
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
template <class T> void XA_IWFR<T>::InvertCTF_DT(const XArray2D<T>& LnIdIin, XArray2D<std::complex<T> >& cout, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC, int NumThreads)
{
	bool bAper(q2max > 0);

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(LnIdIin.GetHeadPtr()->Clone());
	const IXAHWave2D* ph2 = GetIXAHWave2D(LnIdIin);
	ph2->Validate();
	index_t ny = LnIdIin.GetDim1();
	index_t nx = LnIdIin.GetDim2();

	index_t i = 2;
	while (i < ny) i *= 2;
	if (i != ny) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT() (m_dim1 is not a power of 2)");
	index_t j = 2;
	while (j < nx) j *= 2;
	if (j != nx) throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT() (m_dim2 is not a power of 2)");

	if (NumThreads < 1)
		throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT() (number of threads must be >= 1)");
	//omp_set_num_threads(NumThreads); - !!! this should be called earlier in a sginle-thread region of the calling program

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap = abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	double omega = atan(sigma);

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::InvertCTF_DT() (evanescent waves present)");

	// complexify for FFT
	XArray2D<std::complex<T> > camp = MakeComplex(LnIdIin, T(0), false);
	XArray2DFFT<T> xafft(camp);
	xafft.Shuffle();

	T* u = reinterpret_cast<T*>(&(camp.front()));

	//********* Fourier transforming initial amplitude u_n(i,j)
	OouraFft<T> fft;
	fft.Complex2D((std::complex<T> *) u, ny, nx, OouraFft<T>::eDirInv); // NOTE we use inverse FFT in Ooura, because it is the one with the minus sign in the exponent

	//********* Divide F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi = 1.0 / xap;
	double dcsi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;

	double fac2 = PI * wl * defocdist;
	if (alpha <= 0) alpha = 0.1 * pow(abs(fac2) * std::min(dcsi2, deta2), 2);

	double Z1 = defocdist + Z1mZ2d2;
	double Z2 = defocdist - Z1mZ2d2;
	double cosphiA = cos(phiA);
	double sinphiA = sin(phiA);
	double sin2phiA = 2.0 * sinphiA * cosphiA;
	double fac2x = PI * wl * (Z1 * cosphiA * cosphiA + Z2 * sinphiA * sinphiA);
	double fac2y = PI * wl * (Z2 * cosphiA * cosphiA + Z1 * sinphiA * sinphiA);
	double facxy = PI * wl * Z1mZ2d2 * sin2phiA;
	facxy = facxy * deta * dcsi;
	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	cout.Resize(LnIdIin.GetDim1(), LnIdIin.GetDim2()); // filled with zeros
	T* v = reinterpret_cast<T*>(&(cout.front()));

	#pragma omp parallel for
	for (long i = -long(nyd2); i < 0; i++)
	{
		index_t k, kj;
		double eta2, csi2, q2;
		double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg, cosreg(0);

		kj = nxy2 + nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y - omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				v[k] = u[k] * sinreg + u[k + 1] * cosreg;
				v[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
			}
		}
		kj = nxy2 + nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				v[k] = u[k] * sinreg + u[k + 1] * cosreg;
				v[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
			}
		}
	}
	#pragma omp parallel for
	for (long i = 0; i < long(nyd2); i++)
	{
		index_t k, kj;
		double eta2, csi2, q2;
		double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg, cosreg(0);

		kj = nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y - omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				v[k] = u[k] * sinreg + u[k + 1] * cosreg;
				v[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
			}
		}
		kj = nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				v[k] = u[k] * sinreg + u[k + 1] * cosreg;
				v[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
			}
		}
	}

	XArray2DFFT<T> xafft1(cout);
	xafft1.Shuffle();
}

//! Inverts symmetrised CTF on the Ewald sphere for a monomorpous object using FFTW library
// LnIdIin - ln(I/I_in), where I is the defocused image
// camp - output solution in the form of 3D FFT of  delta * [4 * PI * sqrt(1 + sigma^2)] / wl, on the Ewald sphere in the Fourier space, in the "unshuffled form"
// xafft - reference to an externally created fftw2D object with correct nx and ny dimensions
// defocdist - defocus distance
// Z1mZ2d2 astigmatism parameter(dblDistanceX - dblDistanceY) / 2 (in the same units as used in the Wavehead2D)
// phiA astigmatism angle(with the X axis, counterclockwise) (in radians)
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// sigma - beta to delta ratio (1 / gamma), it is equal to zero for pure phase objects
// alpha - Tikhonov regularization parameter
// bESCC - if true - Ewald sphere curvature correction is applied, if false - flat Ewald sphere is simulated
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
// NOTE: differences with InvertCTF_DT: (a) InvertCTF_DT1 uses FFTW library instead of Ooura; 
//       (b) InvertCTF_DT1 allows any even input array dimensions, rather than only integer powers of 2; 
template <class T> void XA_IWFR<T>::InvertCTF_DT1(const XArray2D<T>& LnIdIin, XArray2D<std::complex<T> >& camp, Fftwd2c& xafft, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC, int NumThreads)
{
	if (LnIdIin.GetDim1() != xafft.GetDim1() || LnIdIin.GetDim2() != xafft.GetDim2())
		throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT1() (input 2D array and Fftwd2c object have different dimensions)");

	bool bAper(q2max > 0);

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(LnIdIin.GetHeadPtr()->Clone());
	const IXAHWave2D* ph2 = GetIXAHWave2D(LnIdIin);
	ph2->Validate();
	index_t ny = LnIdIin.GetDim1();
	index_t nx = LnIdIin.GetDim2();

	// check that the input array dimensions are even (FFTW works for odd dimensions as well, but my Shuffle() type routines require even dimensions)
	if (nx != 2 * int(nx / 2) || ny != 2 * int(ny / 2))
		throw std::invalid_argument("input array dimensions must be even in XA_IWFR<T>::InvertCTF_DT1()");

	if (NumThreads < 1)
		throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT1() (number of threads must be >= 1)");
	//omp_set_num_threads(NumThreads); - !!! this should be called earlier in a sginle-thread region of the calling program

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap = abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	double omega = atan(sigma);

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::InvertCTF_DT1() (evanescent waves present)");

	//********* Fourier transforming initial amplitude
	camp = MakeComplex(LnIdIin, T(0), false);
	camp.Shuffle();
	xafft.ForwardFFT(camp); // we re-use the same FFTW plan for multiple arrays of the same dimensions and perform in-place FFT on these "external" arrays
	T* u = reinterpret_cast<T*>(&(camp.front()));

	//********* Divide F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi = 1.0 / xap;
	double dcsi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;

	double fac2 = PI * wl * defocdist;
	if (alpha <= 0) alpha = 0.1 * pow(abs(fac2) * std::min(dcsi2, deta2), 2);

	double Z1 = defocdist + Z1mZ2d2;
	double Z2 = defocdist - Z1mZ2d2;
	double cosphiA = cos(phiA);
	double sinphiA = sin(phiA);
	double sin2phiA = 2.0 * sinphiA * cosphiA;
	double fac2x = PI * wl * (Z1 * cosphiA * cosphiA + Z2 * sinphiA * sinphiA);
	double fac2y = PI * wl * (Z2 * cosphiA * cosphiA + Z1 * sinphiA * sinphiA);
	double facxy = PI * wl * Z1mZ2d2 * sin2phiA;
	facxy = facxy * deta * dcsi;
	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	//#pragma omp parallel for
	for (long i = -long(nyd2); i < 0; i++)
	{
		index_t k, kj;
		double eta2, csi2, q2;
		double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg, cosreg(0);

		kj = nxy2 + nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y - omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp; 
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				dtemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = dtemp;
			}
		}
		kj = nxy2 + nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				dtemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = dtemp;
			}
		}
	}
	//#pragma omp parallel for
	for (long i = 0; i < long(nyd2); i++)
	{
		index_t k, kj;
		double eta2, csi2, q2;
		double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg, cosreg(0);

		kj = nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y - omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				dtemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = dtemp;
			}
		}
		kj = nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				dtemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = dtemp;
			}
		}
	}

	camp.Shuffle();
}

//! Inverts symmetrised CTF on the Ewald sphere for a monomorpous object using FFTWf (single-precision) library
// LnIdIin - ln(I/I_in), where I is the defocused image
// camp - output solution in the form of 3D FFT of  delta * [4 * PI * sqrt(1 + sigma^2)] / wl, on the Ewald sphere in the Fourier space, in the "unshuffled form"
// xafftf - reference to an externally created fftw2Df object with correct nx and ny dimensions
// defocdist - defocus distance
// Z1mZ2d2 astigmatism parameter(dblDistanceX - dblDistanceY) / 2 (in the same units as used in the Wavehead2D)
// phiA astigmatism angle(with the X axis, counterclockwise) (in radians)
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// sigma - beta to delta ratio (1 / gamma), it is equal to zero for pure phase objects
// alpha - Tikhonov regularization parameter
// bESCC - if true - Ewald sphere curvature correction is applied, if false - flat Ewald sphere is simulated
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
// NOTE: differences with InvertCTF_DT: (a) InvertCTF_DT1 uses FFTW library instead of Ooura; 
//       (b) InvertCTF_DT1 allows any even input array dimensions, rather than only integer powers of 2; 
template <class T> void XA_IWFR<T>::InvertCTF_DT2(const XArray2D<T>& LnIdIin, XArray2D<std::complex<T> >& camp, Fftwd2fc& xafftf, double defocdist, double Z1mZ2d2, double phiA, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC, int NumThreads)
{
	if (LnIdIin.GetDim1() != xafftf.GetDim1() || LnIdIin.GetDim2() != xafftf.GetDim2())
		throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT2() (input 2D array and Fftwd2c object have different dimensions)");

	bool bAper(q2max > 0);

	std::unique_ptr<IXAHead> pHead(nullptr);
	pHead.reset(LnIdIin.GetHeadPtr()->Clone());
	const IXAHWave2D* ph2 = GetIXAHWave2D(LnIdIin);
	ph2->Validate();
	index_t ny = LnIdIin.GetDim1();
	index_t nx = LnIdIin.GetDim2();

	// check that the input array dimensions are even (FFTW works for odd dimensions as well, but my Shuffle() type routines require even dimensions)
	if (nx != 2 * int(nx / 2) || ny != 2 * int(ny / 2))
		throw std::invalid_argument("input array dimensions must be even in XA_IWFR<T>::InvertCTF_DT2()");

	if (NumThreads < 1)
		throw std::invalid_argument("invalid_argument in XA_IWFR<T>::InvertCTF_DT2() (number of threads must be >= 1)");
	//omp_set_num_threads(NumThreads); - !!! this should be called earlier in a sginle-thread region of the calling program

	index_t nxd2 = nx / 2;
	index_t nyd2 = ny / 2;
	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	double wl = ph2->GetWl();
	double xlo = ph2->GetXlo();
	double xhi = ph2->GetXhi();
	double xst = ph2->GetXStep(nx);
	double ylo = ph2->GetYlo();
	double yhi = ph2->GetYhi();
	double yst = ph2->GetYStep(ny);
	double xap = abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;

	double omega = atan(sigma);

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in XA_IWFR<T>::InvertCTF_DT() (evanescent waves present)");

	//********* Fourier transforming initial amplitude
	camp = MakeComplex(LnIdIin, T(0), false);
	camp.Shuffle();
	xafftf.ForwardFFT(camp); // we re-use the same FFTW plan for multiple arrays of the same dimensions and perform in-place FFT on these "external" arrays
	T* u = reinterpret_cast<T*>(&(camp.front()));

	//********* Divide F[u] by the CTF
	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dcsi = 1.0 / xap;
	double dcsi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;

	double fac2 = PI * wl * defocdist;
	if (alpha <= 0) alpha = 0.1 * pow(abs(fac2) * std::min(dcsi2, deta2), 2);

	double Z1 = defocdist + Z1mZ2d2;
	double Z2 = defocdist - Z1mZ2d2;
	double cosphiA = cos(phiA);
	double sinphiA = sin(phiA);
	double sin2phiA = 2.0 * sinphiA * cosphiA;
	double fac2x = PI * wl * (Z1 * cosphiA * cosphiA + Z2 * sinphiA * sinphiA);
	double fac2y = PI * wl * (Z2 * cosphiA * cosphiA + Z1 * sinphiA * sinphiA);
	double facxy = PI * wl * Z1mZ2d2 * sin2phiA;
	facxy = facxy * deta * dcsi;
	double fac3 = PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = PI * pow(wl, 5) / 3.0 * Cs5;

	//#pragma omp parallel for
	for (long i = -long(nyd2); i < 0; i++)
	{
		index_t k, kj;
		double eta2, csi2, q2;
		double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg, cosreg(0), ttemp;

		kj = nxy2 + nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y - omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				ttemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = ttemp;
			}
		}
		kj = nxy2 + nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				ttemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = ttemp;
			}
		}
	}
	//#pragma omp parallel for
	for (long i = 0; i < long(nyd2); i++)
	{
		index_t k, kj;
		double eta2, csi2, q2;
		double etafacxy, eta2fac2y, fac2a, dtemp, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg, cosreg(0), ttemp;

		kj = nx2 * i + nx2;
		etafacxy = facxy * double(i);
		eta2 = deta2 * i * i;
		eta2fac2y = eta2 * fac2y - omega;
		for (long j = -long(nxd2); j < 0; j++)
		{
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				ttemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = ttemp;
			}
		}
		kj = nx2 * i;
		for (long j = 0; j < long(nxd2); j++)
		{
			k = kj + 2 * j;
			k = kj + 2 * j;
			csi2 = dcsi2 * j * j;
			fac2a = eta2fac2y + csi2 * fac2x + etafacxy * double(j);
			q2 = csi2 + eta2;
			if (!bAper || q2 < q2max)
			{
				if (bC35) Cstemp = fac3 * q2 * q2 + fac5 * pow(q2, 3);
				dtemp = fac2a + Cstemp;
				sintemp = sin(dtemp);
				sin2 = sintemp * sintemp;
				sinreg = T(sintemp / (sin2 + alpha));
				if (bESCC)
				{
					costemp = cos(dtemp);
					cos2 = costemp * costemp;
					cosreg = T(costemp / (cos2 + alpha));
				}
				ttemp = u[k] * sinreg + u[k + 1] * cosreg;
				u[k + 1] = u[k + 1] * sinreg - u[k] * cosreg;
				u[k] = ttemp;
			}
		}
	}

	camp.Shuffle();
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
	template class XA::XA_IWFR<float>;
	template class XA::XA_IWFR<double>;
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
