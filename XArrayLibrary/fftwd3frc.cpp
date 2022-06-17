// Implementation of the Fftwd3frc class

#include "fftwd3frc.h"
#include "XA_head3.h"

Fftwd3frc::Fftwd3frc(int nz1, int ny1, int nx1, bool bRealSpace1, bool bInPlace1, bool bMeasure)
{
	if (nz1 <= 0 || ny1 <= 0 || nx1 <= 0) throw std::runtime_error("Error: non-positive array dimension passed to Fftwd3frc constructor.");
	nz = nz1; ny = ny1; nx = nx1;
	bInPlace = bInPlace1;
	bContainsData = false;
	bRealSpace = bRealSpace1;
	if (bInPlace)
	{
		pout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * GetNc());
		if (pout == nullptr) throw std::runtime_error("Error: failed to allocate complex 3D array in Fftwd3frc constructor.");
		pin = (float*)pout;
	}
	else
	{
		pin = (float*)fftwf_malloc(sizeof(float) * GetNr());
		if (pin == nullptr)
		{
			Cleanup();
			throw std::runtime_error("Error: failed to allocate real 3D array in Fftwd3frc constructor.");
		}
		pout = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * GetNc());
		if (pout == nullptr)
		{
			Cleanup();
			throw std::runtime_error("Error: failed to allocate complex 3D array in Fftwd3frc constructor.");
		}
	}
	uiflag = bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE;
	aplan = fftwf_plan_dft_r2c_3d(nz, ny, nx, pin, pout, uiflag);
	if (aplan == NULL)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to create forward plan in Fftwd3frc constructor.");
	}
	bplan = fftwf_plan_dft_c2r_3d(nz, ny, nx, pout, pin, uiflag);
	if (bplan == NULL)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to create inverse plan in Fftwd3frc constructor.");
	}
}


void Fftwd3frc::Cleanup()
{
	if (aplan != 0) { fftwf_destroy_plan(aplan); aplan = 0; }
	if (bplan != 0) { fftwf_destroy_plan(bplan); bplan = 0; }
	if (pin != 0) { if (!bInPlace) fftwf_free(pin); pin = 0; }
	if (pout != 0) { fftwf_free(pout); pout = 0; }
	nz = ny = nx = 0;
	bContainsData = false;
}


template <class T> void Fftwd3frc::GetRealXArray3D(xar::XArray3D<T>& aaa, bool bShuffle) const
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::GetRealXArray3D() cannot be called when Fftwd3frc object does not contain meaningful data");
	
	if (!bRealSpace)
		throw std::runtime_error("Fftwd3frc::GetRealXArray3D() cannot be called when Fftwd3frc object is in the Fourier space state");

	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the output real array are different from the internal one in Fftwd3frc::GetRealXArray3D()");

	int nx22 = GetNx2() * 2; // last dimension of the real array in the case bInPlace == true
	size_t m = 0;
	double fnorm = 1.0 / double(GetNr());
	if (bShuffle)
	{
		index_t nxd2 = aaa.GetDim3() / 2, nyd2 = aaa.GetDim2() / 2, nzd2 = aaa.GetDim1() / 2;
		index_t k1, j1, i1;
		for (index_t k = 0; k < nz; k++)
		{
			k < nzd2 ? k1 = k + nzd2 : k1 = k - nzd2;
			for (index_t j = 0; j < ny; j++)
			{
				j < nyd2 ? j1 = j + nyd2 : j1 = j - nyd2;
				for (index_t i = 0; i < nx; i++)
				{
					i < nxd2 ? i1 = i + nxd2 : i1 = i - nxd2;
					aaa[k1][j1][i1] = T(pin[m++] * fnorm); // we assume that the internal real array was produced by inverse FFT
				}
				if (bInPlace) for (int i = nx; i < nx22; i++) m++; // real array padding must be skipped
			}
		}
	}
	else
	{
		for (index_t k = 0; k < nz; k++)
			for (index_t j = 0; j < ny; j++)
			{
				for (index_t i = 0; i < nx; i++)
					aaa[k][j][i] = T(pin[m++] * fnorm); // we assume that the internal real array was produced by inverse FFT
				if (bInPlace) for (int i = nx; i < nx22; i++) m++; // real array padding must be skipped
			}
	}
}


template <class T> void Fftwd3frc::SetRealXArray3D(const xar::XArray3D<T>& aaa, bool bShuffle)
{
	if (!bRealSpace)
		throw std::runtime_error("Fftwd3frc::SetRealXArray3D() cannot be called when Fftwd3frc object is in the Fourier space state");

	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the input real array are different from the internal one in Fftwd3frc::SetRealXArray3D()");

	int nx22 = GetNx2() * 2;
	size_t m = 0;
	if (bShuffle)
	{
		index_t nxd2 = aaa.GetDim3() / 2, nyd2 = aaa.GetDim2() / 2, nzd2 = aaa.GetDim1() / 2;
		index_t k1, j1, i1;
		for (index_t k = 0; k < nz; k++)
		{
			k < nzd2 ? k1 = k + nzd2 : k1 = k - nzd2;
			for (index_t j = 0; j < ny; j++)
			{
				j < nyd2 ? j1 = j + nyd2 : j1 = j - nyd2;
				for (index_t i = 0; i < nx; i++)
				{
					i < nxd2 ? i1 = i + nxd2 : i1 = i - nxd2;
					pin[m++] = (float)aaa[k1][j1][i1];
				}
				if (bInPlace) for (int i = nx; i < nx22; i++) pin[m++] = 0.0f; // real array must padded
			}
		}
	}
	else
	{
		for (index_t k = 0; k < nz; k++)
			for (index_t j = 0; j < ny; j++)
			{
				for (index_t i = 0; i < nx; i++)
					pin[m++] = (float)aaa[k][j][i];
				if (bInPlace) for (int i = nx; i < nx22; i++) pin[m++] = 0.0f; // real array must padded
			}
	}

	bContainsData = true;
}


template <class T> void Fftwd3frc::GetComplexXArray3D(xar::XArray3D< std::complex<T> >& aaa, bool bShuffle) const
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::GetComplexXArray3D() cannot be called when Fftwd3frc object does not contain meaningful data");

	if (bRealSpace)
		throw std::runtime_error("Fftwd3frc::GetComplexXArray3D() cannot be called when Fftwd3frc object is in the real space state");

	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the output complex array disagree with the internal ones in Fftwd3frc::GetComplexXArray3D()");

	int nxc2 = GetNx2();
	size_t m = 0;
	if (bShuffle)
	{
		index_t nxd2 = aaa.GetDim3() / 2, nyd2 = aaa.GetDim2() / 2, nzd2 = aaa.GetDim1() / 2;
		index_t k1, j1, i1;
		for (index_t k = 0; k < nz; k++)
		{
			k < nzd2 ? k1 = k + nzd2 : k1 = k - nzd2;
			for (index_t j = 0; j < ny; j++)
			{
				j < nyd2 ? j1 = j + nyd2 : j1 = j - nyd2;
				for (index_t i = 0; i < nxd2; i++)
				{
					i1 = i + nxd2;
					aaa[k1][j1][i1] = std::complex<T>(pout[m][0], pout[m][1]);
					m++;
				}
				// and now the case of i = nx = 0
				aaa[k1][j1][0] = std::complex<T>(pout[m][0], pout[m][1]);
				m++;
			}
		}
		// and now the Hermitian conjugated "half"
		for (index_t k = 0; k < nz; k++)
		{
			k1 = (k == 0) ? 0 : nz - k;
			for (index_t j = 0; j < ny; j++)
			{
				j1 = (j == 0) ? 0 : ny - j;
				for (index_t i = 1; i < nxd2; i++)
					aaa[k1][j1][i] = std::conj(aaa[k][j][nx - i]);
			}
		}
	}
	else
	{
		for (index_t k = 0; k < nz; k++)
			for (index_t j = 0; j < ny; j++)
			{
				for (index_t i = 0; i < nxc2; i++)
				{
					aaa[k][j][i] = std::complex<T>(pout[m][0], pout[m][1]);
					m++;
				}
			}
		// and now the Hermitian conjugated "half"
		index_t k1, j1;
		for (index_t k = 0; k < nz; k++)
		{
			k1 = (k == 0) ? 0 : nz - k;
			for (index_t j = 0; j < ny; j++)
			{
				j1 = (j == 0) ? 0 : ny - j;
				for (index_t i = nxc2; i < nx; i++)
					aaa[k][j][i] = std::conj(aaa[k1][j1][nx - i]);
			}
		}
	}
}


template <class T> void Fftwd3frc::SetComplexXArray3D(const xar::XArray3D< std::complex<T> >& aaa, bool bShuffle)
{
	if (bRealSpace)
		throw std::runtime_error("Fftwd3frc::SetComplexXArray3D() cannot be called when Fftwd3frc object is in the real space state");

	if ((int)aaa.GetDim1() != nz || (int)aaa.GetDim2() != ny || (int)aaa.GetDim3() != nx)
		throw std::runtime_error("dimensions of the input complex array disagree with the internal ones in Fftwd3frc::SetComplexXArray3D()");

	int nxc2 = GetNx2();
	size_t m = 0;
	if (bShuffle)
	{
		index_t nxd2 = aaa.GetDim3() / 2, nyd2 = aaa.GetDim2() / 2, nzd2 = aaa.GetDim1() / 2;
		index_t k1, j1, i1;
		for (index_t k = 0; k < nz; k++)
		{
			k < nzd2 ? k1 = k + nzd2 : k1 = k - nzd2;
			for (index_t j = 0; j < ny; j++)
			{
				j < nyd2 ? j1 = j + nyd2 : j1 = j - nyd2;
				for (index_t i = 0; i < nxd2; i++)
				{
					i1 = i + nxd2;
					pout[m][0] = (float)aaa[k1][j1][i1].real();
					pout[m][1] = (float)aaa[k1][j1][i1].imag(); // arrays obtained with OouraFFT should be conjugated before this assignment
					m++;
				}
				// and finally the case of i = nx = 0
				pout[m][0] = (float)aaa[k1][j1][0].real();
				pout[m][1] = (float)aaa[k1][j1][0].imag();  // arrays obtained with OouraFFT should be conjugated before this assignment
				m++;
			}
		}
	}
	else
	{
		for (index_t k = 0; k < nz; k++)
			for (index_t j = 0; j < ny; j++)
				for (index_t i = 0; i < nxc2; i++)
				{
					pout[m][0] = (float)aaa[k][j][i].real();
					pout[m][1] = (float)aaa[k][j][i].imag();
					m++;
				}
	}

	bContainsData = true;
}



/*
void Fftwd3frc::PrintRealArray(const char* message)
{
	float* pin = GetReal();
	printf(message);
	size_t m = 0;
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx; i++)
				printf("\npin[%zd,%zd,%zd] = %g", k, j, i, pin[m++]);
}


void Fftwd3frc::PrintComplexArray(const char* message)
{
	printf(message);
	size_t m = 0;
	int nx2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
		for (index_t j = 0; j < ny; j++)
			for (index_t i = 0; i < nx2; i++)
			{
				printf("\npout[%zd,%zd,%zd] = (%g, %g)", k, j, i, pout[m][0], pout[m][1]); 
				m++; 
			}
}
*/

void Fftwd3frc::ForwardFFT()
// Forward FFT
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::ForwardFFT() cannot be called when Fftwd3frc object does not contain meaningful data");

	if (!bRealSpace)
		throw std::runtime_error("Fftwd3frc::ForwardFFT() cannot be called when Fftwd3frc object is in the Fourier space state");

	fftwf_execute(aplan); 
	bRealSpace = false;
}


void Fftwd3frc::InverseFFT()
// Inverse FFT
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::InverseFFT() cannot be called when Fftwd3frc object does not contain meaningful data");

	if (bRealSpace)
		throw std::runtime_error("Fftwd3frc::InverseFFT() cannot be called when Fftwd3frc object is in the real space state");

	fftwf_execute(bplan); 
	bRealSpace = true; 
}


template <class T> void Fftwd3frc::InverseMLaplacian(xar::XArray3D<T>& xa3, double alpha, double norm)
// Regularized inverse 3D (-Laplacian) via multiplication of the FFT of the 3D array by norm / [alpha + by 4 * PI^2 * (ksi^2 + eta^2 + dzeta^2)]
// alpha is the usual Tikhonov regularization parameter
// NOTE: input array is supposed to be represented in the real(!) space
// NOTE: compared to TIE-Hom, alpha = 4 * PI * sigma / (dz * wl), norm = alpha * sqrt(1 + sigma^2) / sigma = 4 * PI * sqrt(1 + sigma * sigma) / (dz * wl)
{
	if (!bRealSpace)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian() cannot be called when Fftwd3frc object is in the Fourier space state");

	if (!bInPlace)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian() should not be called when out-of-place mode is used");

	if (alpha <= 0)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian() cannot be called with alpha <= 0");

	if (nz != xa3.GetDim1() || ny != xa3.GetDim2() || nx != xa3.GetDim3())
		throw std::runtime_error("Error: input XArray3D in InverseLaplacian() has wrong dimensions.");

	IXAHWave3D* ph3 = GetIXAHWave3D(xa3);
	if (ph3 == 0) throw std::runtime_error("Error: input XArray3D in InverseLaplacian() does not have a Wave3D head.");
	double xlo = ph3->GetXlo();
	double xhi = ph3->GetXhi();
	double ylo = ph3->GetYlo();
	double yhi = ph3->GetYhi();
	double zlo = ph3->GetZlo();
	double zhi = ph3->GetZhi();

	// FFT of input array
	SetRealXArray3D(xa3, false);
	ForwardFFT();

	/// multiply FFT of input array by the FFT version of regularized inverse 3D Laplacian
	double fact = 4.0 * xar::PI * xar::PI;
	double dksi2 = fact / ((xhi - xlo) * (xhi - xlo));
	double deta2 = fact / ((yhi - ylo) * (yhi - ylo));
	double dzeta2 = fact / ((zhi - zlo) * (zhi - zlo));
	double dtemp, dk2, djk2;
	float ftemp;

	size_t m = 0;
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2, nxc2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;
		dk2 = k1 * k1 * dzeta2 + alpha;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nxc2; i++)
			{
				dtemp = i * i * dksi2 + djk2;
				ftemp = float(norm / dtemp);
				pout[m][0] *= ftemp;
				pout[m][1] *= ftemp;
				m++;
			}
		}
	}

	// inverse FFT of the product
	InverseFFT();
	GetRealXArray3D(xa3, false);
}


void Fftwd3frc::InverseMLaplacian1(double zaper, double yaper, double xaper, double alpha, double norm)
// Regularized inverse 3D (-Laplacian) via multiplication of the 3D array by norm / [alpha + by 4 * PI^2 * (ksi^2 + eta^2 + dzeta^2)]
// alpha is the usual Tikhonov regularization parameter
// 
// NOTE: input array is supposed to be represented in the Fourier(!) space
// NOTE: compared to TIE-Hom, alpha = 4 * PI * sigma / (dz * wl), norm = alpha * sqrt(1 + sigma^2) / sigma = 4 * PI * sqrt(1 + sigma * sigma) / (dz * wl)
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian1() cannot be called when Fftwd3frc object does not contain meaningful data");

	if (bRealSpace)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian1() cannot be called when Fftwd3frc object is in the real space state");

	if (alpha <= 0 || zaper <= 0 || yaper <= 0 || xaper <= 0)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian1() cannot be called with alpha, xaper, yaper or zaper <= 0");

	/// multiply FFT of input array by the FFT version of regularized inverse 3D Laplacian
	double fact = 4.0 * xar::PI * xar::PI;
	double dksi2 = fact / (xaper * xaper);
	double deta2 = fact / (yaper * yaper);
	double dzeta2 = fact / (zaper * zaper);
	double dtemp, dk2, djk2;
	float ftemp;

	size_t m = 0;
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2, nxc2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;
		dk2 = k1 * k1 * dzeta2 + alpha;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nxc2; i++)
			{
				dtemp = i * i * dksi2 + djk2;
				ftemp = float(norm / dtemp);
				pout[m][0] *= ftemp;
				pout[m][1] *= ftemp;
				m++;
			}
		}
	}
}


void Fftwd3frc::CTFcorrection1(xar::Wavehead3D head3D, double defocdist, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC)
// head3D - structure containing the wavelength and physical dimensions of the 3D array
// defocdist - defocus distance
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// sigma - beta to delta ratio (1 / gamma), it is equal to zero for pure phase objects
// alpha - Tikhonov regularization parameter
// bESCC - if true - Ewald sphere curvature correction is applied, if false - flat Ewald sphere is simulated
//
// NOTE: if alpha<=0 is given in the function call, we actually use alpha = 0.1 times the minimal non-zero CTF^4 - see code below.
// NOTE: this functions allows any even input array dimensions, rather than only integer powers of 2; 
// NOTE: it seems that when bESCC = true, this function may not perform correctly. A different function working on the whole(!) 3D array,
// (rather then on one-half of the array only, as this function does) seems to be performing much better - see CTFcorrection3D()
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::CTFcorrection1() cannot be called when Fftwd3frc object does not contain meaningful data");

	if (bRealSpace)
		throw std::runtime_error("Fftwd3frc::CTFcorrection1() cannot be called when Fftwd3frc object is in the real space state");

	bool bAper(q2max > 0);

	// check that the input array dimensions are even (FFTW works for odd dimensions as well, but my Shuffle() type routines require even dimensions)
	if (nx != 2 * int(nx / 2) || ny != 2 * int(ny / 2) || nz != 2 * int(nz / 2))
		throw std::invalid_argument("input array dimensions must be even in Fftwd3frc::CTFcorrection1()");

	index_t nx2 = nx * 2;
	index_t ny2 = ny * 2;
	index_t nxy = nx * ny;
	index_t nxy2 = nxy * 2;

	double wl = head3D.GetWl();
	double xlo = head3D.GetXlo();
	double xhi = head3D.GetXhi();
	double xst = head3D.GetXStep(nx);
	double ylo = head3D.GetYlo();
	double yhi = head3D.GetYhi();
	double yst = head3D.GetYStep(ny);
	double zlo = head3D.GetZlo();
	double zhi = head3D.GetZhi();
	double zst = head3D.GetZStep(ny);
	double xap = abs(xhi - xlo);
	double xap2 = (xhi - xlo) * (xhi - xlo);
	double yap = abs(yhi - ylo);
	double yap2 = (yhi - ylo) * (yhi - ylo);
	double zap = abs(zhi - zlo);
	double zap2 = (zhi - zlo) * (zhi - zlo);
	double xst2 = xst * xst;
	double yst2 = yst * yst;
	double zst2 = zst * zst;

	//  Nyquist frequency (spectral radius) of U0 = srad*wl
	double srad = sqrt(0.25 / xst2 + 0.25 / yst2 + 0.25 / zst2);
	if (srad * wl > 1.0)
		throw std::runtime_error("runtime_error in Fftwd3frc::CTFcorrection1() (evanescent waves present)");

	double omega = atan(sigma);

	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dksi = 1.0 / xap;
	double dksi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;
	double dzeta = 1.0 / zap;
	double dzeta2 = 1.0 / zap2;

	double fac2 = xar::PI * wl * defocdist;
	if (alpha <= 0) alpha = 0.1 * pow(abs(fac2) * std::min(dksi2, std::min(deta2, dzeta2)), 2);
	double fac3 = xar::PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = xar::PI * pow(wl, 5) / 3.0 * Cs5;

	double zeta2, eta2, ksi2, dtemp, zeta2fac2, eta2fac2, ksi2fac2, sintemp, sin2, sinreg, costemp, cos2, cosreg(0), Cstemp(0);
	float ftemp;

	size_t m = 0;
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2, nxc2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;

		zeta2 = dzeta2 * k1 * k1;
		zeta2fac2 = fac2 * zeta2 - omega;

		//dk2 = k1 * k1 * dzeta2 + alpha;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			
			eta2 = deta2 * j1 * j1;
			eta2fac2 = fac2 * eta2 + zeta2fac2;
			eta2 += zeta2;

			//djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nxc2; i++)
			{
				ksi2 = dksi2 * i * i;
				ksi2fac2 = fac2 * ksi2 + eta2fac2;
				ksi2 += eta2;
				
				if (!bAper || ksi2 < q2max)
				{
					if (bC35) Cstemp = fac3 * ksi2 * ksi2 + fac5 * pow(ksi2, 3);
					dtemp = ksi2fac2 + Cstemp;
					sintemp = sin(dtemp);
					sin2 = sintemp * sintemp;
					sinreg = sintemp / (sin2 + alpha);
					if (bESCC)
					{
						costemp = cos(dtemp);
						cos2 = costemp * costemp;
						cosreg = costemp / (cos2 + alpha);
					}
					ftemp = float(pout[m][0] * sinreg + pout[m][1] * cosreg);
					pout[m][1] = float(pout[m][1] * sinreg - pout[m][0] * cosreg);
					pout[m][0] = ftemp;
				}
				m++;
			}
		}
	}

}

template <class T> void Fftwd3frc::GaussFilter(xar::XArray3D<T>& xa3, double sigma)
// Convolution with a 3D Gaussian with the standard deviation equal to sigma
// !!! NOTE that this filter preserves the average pixel value
{
	if (!bRealSpace)
		throw std::runtime_error("Fftwd3frc::GausFilter() cannot be called when Fftwd3frc object is in the Fourier space state");

	if (!bInPlace)
		throw std::runtime_error("Fftwd3frc::GaussFilter() should not be called when out-of-place mode is used");

	if (nz != xa3.GetDim1() || ny != xa3.GetDim2() || nx != xa3.GetDim3())
		throw std::runtime_error("Error: input XArray3D in GaussFilter() has wrong dimensions.");

	IXAHWave3D* ph3 = GetIXAHWave3D(xa3);
	if (ph3 == 0) throw std::runtime_error("Error: input XArray3D in GaussFilter() does not have a Wave3D head.");
	double xlo = ph3->GetXlo();
	double xhi = ph3->GetXhi();
	double ylo = ph3->GetYlo();
	double yhi = ph3->GetYhi();
	double zlo = ph3->GetZlo();
	double zhi = ph3->GetZhi();

	double dnormA = xa3.Norm(xar::eNormAver);

	// FFT of input array
	SetRealXArray3D(xa3, false);
	ForwardFFT();

	/// multiply FFT of input array by the FFT transform of the 3D Gaussian
	double fact = 2.0 * xar::PI * xar::PI * sigma * sigma;
	double dksi2 = fact / ((xhi - xlo) * (xhi - xlo));
	double deta2 = fact / ((yhi - ylo) * (yhi - ylo));
	double dzeta2 = fact / ((zhi - zlo) * (zhi - zlo));
	double dk2, djk2;
	float ftemp;

	size_t m = 0;
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2, nc2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;
		dk2 = k1 * k1 * dzeta2;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nc2; i++)
			{
				ftemp = (float)exp(-(i * i * dksi2 + djk2));
				pout[m][0] *= ftemp;
				pout[m][1] *= ftemp;
				m++;
			}
		}
	}

	// inverse FFT of the product
	InverseFFT();
	GetRealXArray3D(xa3, false);

	// restore average value
	xa3 += float(dnormA - xa3.Norm(xar::eNormAver));
}


void Fftwd3frc::GaussFilter1(double zaper, double yaper, double xaper, double sigma)
// Convolution with a 3D Gaussian with the (real-space) standard deviation equal to sigma
{
	if (!bContainsData)
		throw std::runtime_error("Fftwd3frc::GaussFilter1() cannot be called when Fftwd3frc object does not contain meaningful data");

	if (bRealSpace)
		throw std::runtime_error("Fftwd3frc::GaussFilter1() cannot be called when Fftwd3frc object is in the real space state");

	if (sigma <= 0 || zaper <= 0 || yaper <= 0 || xaper <= 0)
		throw std::runtime_error("Fftwd3frc::InverseMLaplacian1() cannot be called with sigma, xaper, yaper or zaper <= 0");

	/// multiply FFT of input array by the FFT transform of the 3D Gaussian
	double fact = 2.0 * xar::PI * xar::PI * sigma * sigma;
	double dksi2 = fact / (xaper * xaper);
	double deta2 = fact / (yaper * yaper);
	double dzeta2 = fact / (zaper * zaper);
	double dk2, djk2;
	float ftemp;

	size_t m = 0;
	index_t k1, j1, nyd2 = ny / 2, nzd2 = nz / 2, nc2 = GetNx2();
	for (index_t k = 0; k < nz; k++)
	{
		k <= nzd2 ? k1 = k : k1 = nz - k;
		dk2 = k1 * k1 * dzeta2;
		for (index_t j = 0; j < ny; j++)
		{
			j <= nyd2 ? j1 = j : j1 = ny - j;
			djk2 = j1 * j1 * deta2 + dk2;
			for (index_t i = 0; i < nc2; i++)
			{
				ftemp = (float)exp(-(i * i * dksi2 + djk2));
				pout[m][0] *= ftemp;
				pout[m][1] *= ftemp;
				m++;
			}
		}
	}
}

// explicit instantiation of template functions (Visual C++ seems not clever enough to know when to do it automatically)
template void Fftwd3frc::GetRealXArray3D(xar::XArray3D<float>& aaa, bool bShuffle) const;
template void Fftwd3frc::GetRealXArray3D(xar::XArray3D<double>& aaa, bool bShuffle) const;
template void Fftwd3frc::SetRealXArray3D(const xar::XArray3D<float>& aaa, bool bShuffle);
template void Fftwd3frc::SetRealXArray3D(const xar::XArray3D<double>& aaa, bool bShuffle);
template void Fftwd3frc::GetComplexXArray3D(xar::XArray3D< std::complex<float> >& aaa, bool bShuffle) const;
template void Fftwd3frc::GetComplexXArray3D(xar::XArray3D< std::complex<double> >& aaa, bool bShuffle) const;
template void Fftwd3frc::SetComplexXArray3D(const xar::XArray3D< std::complex<float> >& aaa, bool bShuffle);
template void Fftwd3frc::SetComplexXArray3D(const xar::XArray3D< std::complex<double> >& aaa, bool bShuffle);
template void Fftwd3frc::GaussFilter(xar::XArray3D<float>& xa3, double sigma);
template void Fftwd3frc::GaussFilter(xar::XArray3D<double>& xa3, double sigma);
template void Fftwd3frc::InverseMLaplacian(xar::XArray3D<float>& xa3, double alpha, double norm);
template void Fftwd3frc::InverseMLaplacian(xar::XArray3D<double>& xa3, double alpha, double norm);