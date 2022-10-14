#include "XA_fft2.h"

//! Returns the xar::_eValueType corresponding to T=float
template<> xar::_eValueType xar::XArray2DFFT<float>::GetValuetype() { return eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> xar::_eValueType xar::XArray2DFFT<double>::GetValuetype() { return eXADouble; }


//---------------------------------------------------------------------------
//Function XArray2DFFT<T>::FresnelB
//
//	 Calculates 2D Fresnel integral (uses FFTWf)
//	!!!only works with T = float for now (because of the specificity of the FFTW library used here)
//
/*!
	\brief		Calculates 2D Fresnel integral with possible astigamtism (uses FFTW)
	\param		xafftf Reference to an externally created fftw2Df object with correct nx and ny dimensions
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
InterfaceFFT2D.FresnelB(...); // calculate free space propagation
\endverbatim
*/
template <> void xar::XArray2DFFT<float>::FresnelB(Fftwd2fc& xafftf, double dblDistance, double Z1mZ2d2, double phiA, bool bCheckValidity, double q2max, double C3, double C5)
{
	index_t nx = m_rXArray2D.GetDim2();
	index_t ny = m_rXArray2D.GetDim1();

	IXAHWave2D* ph2 = GetIXAHWave2D(m_rXArray2D);
	if (!ph2)
		throw std::invalid_argument("invalid_argument 'm_rXArray2D' in XArray2DFFT<float>::FresnelB (no Wavehead2D)");
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
		throw std::runtime_error("runtime_error in XArray2DFFT<float>::FresnelB (evanescent waves present)");

	if (dblDistance == 0) dblDistance = wl; // there may be still some work to be done even if dblDistance == 0, but we need to avoid division by zero

	//  Fresnel numbers NFx(y)=a(b)^2/wl/abs(dblDistance)=',fnumx,fnumy
	double absdistanceX = fabs(dblDistance + Z1mZ2d2);
	double absdistanceY = fabs(dblDistance - Z1mZ2d2);
	double fnumx = xap2 / absdistanceX / wl;
	double fnumy = yap2 / absdistanceY / wl;
	if (bCheckValidity && (fnumx < 10 || fnumy < 10))
		throw std::runtime_error("runtime_error in XArray2DFFT<float>::FresnelB (Fresnel number too low)");

	//********* Fourier transforming initial amplitude u(i,j)
	if (m_rXArray2D.GetValuetype() == eXAFComplex)
	{
		//Fftwd2fc xafftwdfc(reinterpret_cast<XArray2D<fcomplex>&>(m_rXArray2D), 1);
		xafftf.InverseFFT(m_rXArray2D); // we are keeping the code consistent with the use of OouraFFT in FresnelA(), taking into account that forward and inverse FFTs are swapped in the two libraries
	}
	else
		//if (m_rXArray2D.GetValuetype() == eXADComplex)
		//{
			//Fftwd2c xafftwdc(reinterpret_cast<XArray2D<dcomplex>&>(m_rXArray2D), 1);
			//xafftwdc.InverseFFT();  // we are keeping the code consistent with the use of OouraFFT in FresnelA(), taking into account that forward and inverse FFTs are swapped in the two libraries
		//}
		//else 
		throw std::runtime_error("runtime_error in XArray2DFFT<float>::FresnelB (unsupported data type)");

	float* u = reinterpret_cast<float*>(&(m_rXArray2D.front()));

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
	float ttemp;

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
				ttemp = u[k] * float(ctemp.real()) - u[k + 1] * float(ctemp.imag());
				u[k + 1] = u[k] * float(ctemp.imag()) + u[k + 1] * float(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = float(0);
				u[k + 1] = float(0);
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
				ttemp = u[k] * float(ctemp.real()) - u[k + 1] * float(ctemp.imag());
				u[k + 1] = u[k] * float(ctemp.imag()) + u[k + 1] * float(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = float(0);
				u[k + 1] = float(0);
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
				ttemp = u[k] * float(ctemp.real()) - u[k + 1] * float(ctemp.imag());
				u[k + 1] = u[k] * float(ctemp.imag()) + u[k + 1] * float(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = float(0);
				u[k + 1] = float(0);
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
				ttemp = u[k] * float(ctemp.real()) - u[k + 1] * float(ctemp.imag());
				u[k + 1] = u[k] * float(ctemp.imag()) + u[k + 1] * float(ctemp.real());
				u[k] = ttemp;
			}
			else
			{
				u[k] = float(0);
				u[k + 1] = float(0);
			}
		}
	}

	//********* inverse Fourier transforming F[u]*Kirchhoff_propagator

	if (m_rXArray2D.GetValuetype() == eXAFComplex)
	{
		//Fftwd2fc xafftwdfc(reinterpret_cast<XArray2D<fcomplex>&>(m_rXArray2D), 1);
		xafftf.ForwardFFT(m_rXArray2D); // we are keeping the code consistent with the use of OouraFFT in FresnelA(), taking into account that forward and inverse FFTs are swapped in the two libraries
	}
	//else
		//if (m_rXArray2D.GetValuetype() == eXADComplex)
		//{
			//Fftwd2c xafftwdc(reinterpret_cast<XArray2D<dcomplex>&>(m_rXArray2D), 1);
			//xafftwdc.ForwardFFT();  // we are keeping the code consistent with the use of OouraFFT in FresnelA(), taking into account that forward and inverse FFTs are swapped in the two libraries
		//}

	float fact = float(1.0 / double(nxy));
	for (k = 0; k < nxy2; k++)	u[k] *= fact;
}


template <> void xar::XArray2DFFT<double>::FresnelB(Fftwd2fc& xafftf, double dblDistance, double Z1mZ2d2, double phiA, bool bCheckValidity, double q2max, double C3, double C5)
{
	throw std::invalid_argument("invalid_argument in XArray2DFFT<double>::FresnelB (this function is not implemented for T = double)");
}

