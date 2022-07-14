// Implementation of the Fftwd2fc class

#include "fftwd2fc.h"
//#include "omp.h"
//#include <string_view> // for std::swap


Fftwd2fc::Fftwd2fc(int ny1, int nx1, int nThreads, bool bMeasure)
//! Constructor allocating internal storage and creating "FFT plans"
{
	bExtStorage = false;
	if (ny1 <= 0 || nx1 <= 0) throw std::runtime_error("Error: non-positive array dimension passed to Fftwd2fc constructor.");
	if (nThreads <= 0) throw std::runtime_error("Error: non-positive number of threads in  Fftwd2fc constructor.");
	ny = ny1; nx = nx1;
	parr = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * ny * nx);
	if (parr == nullptr)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to allocate memory in Fftwd2fc constructor.");
	}
	uiflag = bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE;
	fftwf_plan_with_nthreads(nThreads);
	aplan = fftwf_plan_dft_2d(ny, nx, parr, parr, FFTW_FORWARD, uiflag);
	if (aplan == NULL)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to create forward plan in Fftwd2fc constructor.");
	}
	fftwf_plan_with_nthreads(nThreads);
	bplan = fftwf_plan_dft_2d(ny, nx, parr, parr, FFTW_BACKWARD, uiflag);
	if (bplan == NULL)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to create inverse plan in Fftwd2fc constructor.");
	}
}


Fftwd2fc::Fftwd2fc(xar::XArray2D<xar::fcomplex>& xarc2D, int nThreads, bool bMeasure)
//! Constructor using external storage and creating "FFT plans"
{
	if (nThreads <= 0) throw std::runtime_error("Error: non-positive number of threads in  Fftwd2fc constructor.");
	bExtStorage = true;
	ny = (int)xarc2D.GetDim1(); nx = (int)xarc2D.GetDim2();
	if (nx <= 0 || ny <= 0) throw std::runtime_error("Error: non-positive array dimension in  Fftwd2fc constructor.");
	parr = (fftwf_complex*)&xarc2D[0][0];
	uiflag = bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE;
	fftwf_plan_with_nthreads(nThreads);
	aplan = fftwf_plan_dft_2d(ny, nx, parr, parr, FFTW_FORWARD, uiflag);
	if (aplan == NULL)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to create forward plan in Fftwd2fc constructor.");
	}
	fftwf_plan_with_nthreads(nThreads);
	bplan = fftwf_plan_dft_2d(ny, nx, parr, parr, FFTW_BACKWARD, uiflag);
	if (bplan == NULL)
	{
		Cleanup();
		throw std::runtime_error("Error: failed to create inverse plan in Fftwd2fc constructor.");
	}
}


void Fftwd2fc::Cleanup()
{
	if (aplan != 0) { fftwf_destroy_plan(aplan); aplan = 0; }
	if (bplan != 0) { fftwf_destroy_plan(bplan); bplan = 0; }
	if (!bExtStorage && parr != 0) { fftwf_free(parr); parr = 0; }
	ny = nx = 0;
}


void Fftwd2fc::GetXArray2D(xar::XArray2D<xar::fcomplex>& aaa) const
{
	if (bExtStorage)
		throw std::runtime_error("Fftwd2fc::GetXArray2D() cannot be called in the external array storage mode");
	if ((int)aaa.GetDim1() != ny || (int)aaa.GetDim2() != nx)
		throw std::runtime_error("dimensions of the output complex array are different from the internal one in Fftwd2fc::GetXArray2D()");
	int m = 0;
	for (index_t j = 0; j < ny; j++)
		for (index_t i = 0; i < nx; i++)
		{ 
			aaa[j][i] = xar::fcomplex(parr[m][0], parr[m][1]); 
			m++; 
		}
}


void Fftwd2fc::SetXArray2D(const xar::XArray2D<xar::fcomplex>& aaa)
{
	if (bExtStorage)
		throw std::runtime_error("Fftwd2fc::SetXArray2D() cannot be called in the external array storage mode");
	if ((int)aaa.GetDim1() != ny || (int)aaa.GetDim2() != nx)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwd2fc::SetXArray2D()");
	int m = 0;
	for (index_t j = 0; j < ny; j++)
		for (index_t i = 0; i < nx; i++)
		{
			parr[m][0] = aaa[j][i].real();
			parr[m][1] = aaa[j][i].imag();
			m++;
		}
}

/*
void Fftwd2fc::PrintComplexArray(const char* message)
{
	printf(message);
	int m = 0;
	for (index_t j = 0; j < ny; j++)
		for (index_t i = 0; i < nx; i++)
		{
			printf("\npout[%zd,%zd] = (%g, %g)", j, i, parr[m][0], parr[m][1]); 
			m++; 
		}
}
*/

//! Reshuffles the internal array 
// NOTE: twice shuffling = identity, hence inverse shuffle = shuffle.
// NOTE: it has been checked that this shuffle is the same as xar::XArray2DFFT<double>::Shuffle()
void Fftwd2fc::Shuffle()
{
	int nyd2 = ny / 2, m, m1;
	int nxd2 = nx / 2;

	for (int j = 0; j < nyd2; j++)
	{
		m = nx * j;
		m1 = nx * (j + nyd2) + nxd2;
		for (int i = 0; i < nxd2; i++) 
			std::swap(parr[m++], parr[m1++]);
	}
	for (int j = nyd2; j < ny; j++)
	{
		m = nx * j;
		m1 = nx * (j - nyd2) + nxd2;
		for (int i = 0; i < nxd2; i++) 
			std::swap(parr[m++], parr[m1++]);
	}
}


// Forward transform using external arrays with an alrady existing plan (speeds up repeated transforms with the same dimensions)
void Fftwd2fc::ForwardFFT(xar::XArray2D<xar::fcomplex>& aaa)
{
	if ((int)aaa.GetDim1() != ny || (int)aaa.GetDim2() != nx)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwd2fc::ForwardFFT(xar::XArray2D<xar::fcomplex>& xarc2D)");

	fftwf_execute_dft(aplan, (fftwf_complex*)&aaa[0][0], (fftwf_complex*)&aaa[0][0]);
}


// Inverse transform using external arrays with an alrady existing plan (speeds up repeated transforms with the same dimensions)
void Fftwd2fc::InverseFFT(xar::XArray2D<xar::fcomplex>& aaa)
{
	if ((int)aaa.GetDim1() != ny || (int)aaa.GetDim2() != nx)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwd2fc::InverseFFT(xar::XArray2D<xar::fcomplex>& xarc2D)");

	fftwf_execute_dft(bplan, (fftwf_complex*)&aaa[0][0], (fftwf_complex*)&aaa[0][0]);
}



void Fftwd2fc::InverseMLaplacian1(xar::XArray2D<xar::fcomplex>& aaa, double alpha, double norm)
// Regularized inverse 2D (-Laplacian) via multiplication of the 2D array by norm / [alpha + 4 * PI^2 * (ksi^2 + eta^2)]
// alpha is the usual Tikhonov regularization parameter
// norm is a normalization factor which can be chosen depending on the application
// NOTE: compared to TIE-Hom, alpha = 4 * PI / ((delta / beta) * dz * wl), norm = alpha
// NOTE!!: this function does NOT do the final inverse FFT, leaving the array in the Fourier space in the unshuffled form
{
	if (alpha == 0)
		throw std::runtime_error("Fftwd2fc::InverseMLaplacian1() cannot be called with alpha == 0");

	if ((int)aaa.GetDim1() != ny || (int)aaa.GetDim2() != nx)
		throw std::runtime_error("dimensions of the input complex array are different from the internal one in Fftwd2fc::InverseMLaplacian1()");

	int nxd2 = int(nx / 2), nyd2 = int(ny / 2);
	if (nx != 2 * nxd2 || ny != 2 * nyd2)
		throw std::invalid_argument("input array dimensions must be even in Fftwd2fc::InverseMLaplacian1()");

	if (aaa.GetHeadPtr() == nullptr)
		throw std::runtime_error("input complex array does not have a header in Fftwd2fc::InverseMLaplacian1()");

	xar::Wavehead2D head2d = *dynamic_cast<const xar::Wavehead2D*>(aaa.GetHeadPtr());
	double xaper = head2d.GetXhi() - head2d.GetXlo();
	double yaper = head2d.GetYhi() - head2d.GetYlo();

	aaa.Shuffle();
	ForwardFFT(aaa);
	aaa.Shuffle();

	double fact = 4.0 * xar::PI * xar::PI;
	double dksi2 = fact / (xaper * xaper), ksi2;
	double deta2 = fact / (yaper * yaper), eta2;

	int i1, j1;
	for (int j = 0; j < ny; j++)
	{
		j1 = j - nyd2;
		eta2 = j1 * j1 * deta2 + alpha;

		for (int i = 0; i < nx; i++)
		{
			i1 = i - nxd2;
			ksi2 = dksi2 * i1 * i1 + eta2;
			aaa[j][i] *= float(norm / ksi2);
		}
	}

}
