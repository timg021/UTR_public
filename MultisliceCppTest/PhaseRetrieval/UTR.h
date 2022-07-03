#pragma once
#include <complex.h>

#include "XA_ini.h"
#include "IXAHWave.h"
#include "XAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"

constexpr int natommax(500000); // maximum number of atomic positions that can be saved in an XYZ file

constexpr double dTargetSNR = 1.414213562373; // sqrt(2) targeted cut-off SNR in Fourier coefficients
constexpr float fTargetSNR = float(dTargetSNR); 

int FindPeaks(xar::XArray3D<float> K3out, double datomsizeXY, double datomsizeZ, int natommax, std::string XYZfilename); // this function is defined in FindPeaks.cpp

template <class T> void CTFcorrection3D(xar::XArray3D< std::complex<T> >& xar3D, xar::Wavehead3D head3D, double defocdist, double q2max, double Cs3, double Cs5, double sigma, double epsilon, bool bESCC)
// Performs 3D CTF correction with optional Ewald sphere correction
// 
// xar3D - input array in the Fourier(!) space, in the form of 3D FFT of  delta * [4 * PI * sqrt(1 + sigma^2)] / wl, on the Ewald sphere, in the "unshuffled form"
// head3D - structure containing the wavelength and physical dimensions of the 3D array
// defocdist - defocus distance
// q2max - maximum Fourier frequency (bandpass)
// Cs3 - third spherical aberration
// Cs5 - fifth spherical aberration
// sigma - beta to delta ratio (1 / gamma), it is equal to zero for pure phase objects
// epsilon - Tikhonov regularization parameter
// bESCC - if true - Ewald sphere curvature correction is applied, if false - flat Ewald sphere is simulated
//
// NOTE: input array is supposed to be represented in the Fourier(!) space in a non-wrapped state(!)
// NOTE: if epsilon<=0 is given in the function call, we actually use epsilon = 0.1 times the minimal non-zero CTF^4 - see code below.
// NOTE: this functions allows any even input array dimensions, rather than only integer powers of 2; 
{
	bool bAper(q2max > 0);

	index_t nz = xar3D.GetDim1(), ny = xar3D.GetDim2(), nx = xar3D.GetDim3();

	// check that the input array dimensions are even
	if (nx != 2 * int(nx / 2) || ny != 2 * int(ny / 2) || nz != 2 * int(nz / 2))
		throw std::invalid_argument("input array dimensions must be even in CTFcorrection3D()");

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
		throw std::runtime_error("runtime_error in CTFcorrection3D() (evanescent waves present)");

	double omega = atan(sigma);

	bool bC35((Cs3 != 0) || (Cs5 != 0));
	double dksi = 1.0 / xap;
	double dksi2 = 1.0 / xap2;
	double deta = 1.0 / yap;
	double deta2 = 1.0 / yap2;
	double dzeta = 1.0 / zap;
	double dzeta2 = 1.0 / zap2;

	double fac2 = xar::PI * wl * defocdist;
	if (epsilon <= 0) epsilon = 0.1 * pow(abs(fac2) * std::min(dksi2, std::min(deta2, dzeta2)), 2);
	double fac3 = xar::PI * pow(wl, 3) / 2.0 * Cs3;
	double fac5 = xar::PI * pow(wl, 5) / 3.0 * Cs5;

	int nxd2 = (int)(nx / 2), nyd2 = (int)(ny / 2), nzd2 = (int)(nz / 2);
	#pragma omp parallel for
	for (int k = 0; k < nz; k++)
	{
		int k1, j1, i1;
		double zeta2, eta2, ksi2, dtemp, zeta2fac2, eta2fac2, ksi2fac2, sintemp, sin2, costemp, cos2, Cstemp(0);
		T sinreg(0), cosreg(0);

		k1 = k - nzd2;

		zeta2 = dzeta2 * k1 * k1;
		zeta2fac2 = fac2 * zeta2 - omega;

		for (int j = 0; j < ny; j++)
		{
			j1 = j - nyd2;

			eta2 = deta2 * j1 * j1;
			eta2fac2 = fac2 * eta2 + zeta2fac2;
			eta2 += zeta2;

			for (int i = 0; i < nx; i++)
			{
				i1 = i - nxd2;

				ksi2 = dksi2 * i1 * i1;
				ksi2fac2 = fac2 * ksi2 + eta2fac2;
				ksi2 += eta2;

				if (!bAper || ksi2 < q2max)
				{
					if (bC35) Cstemp = fac3 * ksi2 * ksi2 + fac5 * pow(ksi2, 3);
					dtemp = ksi2fac2 + Cstemp;
					sintemp = sin(dtemp);
					sin2 = sintemp * sintemp;
					sinreg = T(sintemp / (sin2 + epsilon));
					if (bESCC)
					{
						costemp = cos(dtemp);
						cos2 = costemp * costemp;
						cosreg = T(costemp / (cos2 + epsilon));
					}
					xar3D[k][j][i] *= std::complex<T>(sinreg, -cosreg);
				}
			}
		}
	}

}

template <class T> void InverseMLaplacian3D(xar::XArray3D< std::complex<T> >& xar3D, double zaper, double yaper, double xaper, double alpha, double norm)
// Replaces the input 3D arrays with its regularized inverse (-Laplacian), 
// i.e. multiplies each point by norm / [alpha + 4 * PI^2 * (ksi^2 + eta^2 + dzeta^2))
// 
// xar3D - input array in the Fourier(!) space
// zpaper, yaper and xaper are physical apertures of the 3D array in the corresponding dimensions
// alpha is the usual Tikhonov regularization parameter
// norm is a normalization factor which can be chosen depending on the application
// 
// NOTE: input array is supposed to be represented in the Fourier(!) space in a non-wrapped state(!)
// NOTE: the array dimensions must be even
// NOTE: compared to TIE-Hom, alpha = 4 * PI / ((delta / beta) * dz * wl), norm = alpha
{
	if (alpha <= 0)
		throw std::runtime_error("InverseMLaplacian3D() cannot be called with alpha <= 0");

	if (xaper == 0 || yaper == 0 || zaper == 0)
		throw std::runtime_error("InverseMLaplacian3D() cannot be called with an aperture parameter = 0");

	index_t nz = xar3D.GetDim1(), ny = xar3D.GetDim2(), nx = xar3D.GetDim3();
	int nxd2 = int(nx / 2), nyd2 = int(ny / 2), nzd2 = int(nz / 2);
	if (nx != 2 * nxd2 || ny != 2 * nyd2 || nz != 2 * nzd2)
		throw std::invalid_argument("input array dimensions must be even in InverseMLaplacian3D()");

	double fact = 4.0 * xar::PI * xar::PI;
	double dksi2 = fact / (xaper * xaper);
	double deta2 = fact / (yaper * yaper);
	double dzeta2 = fact / (zaper * zaper);

	#pragma omp parallel for
	for (int k = 0; k < nz; k++)
	{
		int k1, j1, i1;
		double zeta2, eta2, ksi2;
		k1 = k - nzd2;
		zeta2 = k1 * k1 * dzeta2 + alpha;
		
		for (int j = 0; j < ny; j++)
		{
			j1 = j - nyd2;
			eta2 = j1 * j1 * deta2 + zeta2;
		
			for (int i = 0; i < nx; i++)
			{
				i1 = i - nxd2;
				ksi2 = dksi2 * i1 * i1 + eta2;
				xar3D[k][j][i] *= T(norm / ksi2);
			}
		}
	}
}


template <class T> void CT_3Dgridding(const xar::XArray2D<xar::dcomplex>& K2four, xar::XArray3D< std::complex<T> >& V3, xar::XArray3D<T>& Samp3, 
	double angleY, double angleZ, T fxlo, T fxst, T fylo, T fyst, T fzlo, T fzst, 
	std::complex<T> Fxc, std::complex<T> Fyc, std::complex<T> Fzc, double wl, bool bRotCentreShift, bool bESCC)
{
	int nx = (int)Samp3.GetDim3();
	int ny = (int)Samp3.GetDim2();
	int nz = (int)Samp3.GetDim1();
	int nx2 = int(nx) - 2, ny2 = int(ny) - 2;
	int nzd2 = int(nz / 2), nz2 = int(nz) - 2;
	index_t nxny = (nx * ny);
	double wl2 = 0.5 * wl;
	T sinangleY = (T)sin(angleY);
	T cosangleY = (T)cos(angleY);
	T sinangleZ = (T)sin(angleZ);
	T cosangleZ = (T)cos(angleZ);

	// calculate the coordinate illumination angle parameters
	double xxx;
	std::vector<T> xxx2(nx);
	std::vector<T> x_sinangleY(nx), x_cosangleY(nx);
	for (int i = 0; i < nx; i++)
	{
		xxx = fxlo + fxst * i;
		xxx2[i] = T(-wl2 * xxx * xxx);
		x_sinangleY[i] = T(xxx * sinangleY);
		x_cosangleY[i] = T(xxx * cosangleY);
	}

	double yyy;
	std::vector<T> yyy2(ny);
	std::vector<T> y_sinangleZ(ny), y_cosangleZ(ny);
	for (int j = 0; j < ny; j++)
	{
		yyy = fylo + fyst * j;
		yyy2[j] = T(-wl2 * yyy * yyy);
		y_sinangleZ[j] = T(yyy * sinangleZ);
		y_cosangleZ[j] = T(yyy * cosangleZ);
	}

	// qz rotated values are expressed via qx and qy - see below

	//#pragma omp parallel for shared (V3, Samp3, K2four, xxx2, yyy2, x_sinangleY, x_cosangleY, y_sinangleZ, y_cosangleZ)
	for (int j = 0; j < ny; j++)
	{
		int ii, jj, nn;
		xar::index_t nji;
		T xx, zz, xxx, yyy, zzz;
		T dx0, dx1, dy0, dy1, dz0, dz1, dxyz;
		std::complex<T> cd;
		T* pFilt3 = &(Samp3[0][0][0]);
		std::complex<T>* pV3 = &(V3[0][0][0]);

		for (int i = 0; i < nx; i++)
		{
			// inverse rotation around Y' axis
			if (bESCC) zz = xxx2[i] + yyy2[j];
			else zz = 0.0; // conventional CT
			xx = x_cosangleY[i] - zz * sinangleY; // qx coordinate after the rotation around Y'
			zzz = x_sinangleY[i] + zz * cosangleY; // qz coordinate after the rotation around Y'
			dz1 = (zzz - fzlo) / fzst;
			nn = (int)floor(dz1); if (nn < 0 || nn > nz2) continue;
			dz1 -= nn; dz0 = 1.0f - dz1;

			// inverse rotation around Z axis
			xxx = xx * cosangleZ + y_sinangleZ[j]; // qx coordinate after the rotation around Z
			yyy = -xx * sinangleZ + y_cosangleZ[j]; // qy coordinate after the rotation around Z

			dx1 = (xxx - fxlo) / fxst;
			ii = (int)floor(dx1); if (ii < 0 || ii > nx2) continue;
			dx1 -= ii; dx0 = 1.0f - dx1;
			dy1 = (yyy - fylo) / fyst;
			jj = (int)floor(dy1); if (jj < 0 || jj > ny2) continue;
			dy1 -= jj; dy0 = 1.0f - dy1;

			cd = K2four[j][i];
			// multiply by linear phase factors due to the shifted centre of rotation in the real space
			if (bRotCentreShift) cd *= exp(Fxc * xxx + Fyc * yyy + Fzc * zzz);

			dxyz = dx0 * dy0 * dz0;
			nji = nn * nxny + jj * nx + ii;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx1 * dy0 * dz0;
			nji += 1;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx0 * dy0 * dz1;
			nji += -1 + nxny;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx1 * dy0 * dz1;
			nji += 1;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx0 * dy1 * dz0;
			nji += -1 - nxny + nx;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx1 * dy1 * dz0;
			nji += 1;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx0 * dy1 * dz1;
			nji += -1 + nxny;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;

			dxyz = dx1 * dy1 * dz1;
			nji += 1;
			*(pFilt3 + nji) += dxyz;
			*(pV3 + nji) += dxyz * cd;
		}
	} // end of cycle updating the 3D potential or beta at the current rotational orientation "na"
}
