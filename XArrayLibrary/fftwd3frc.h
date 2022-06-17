//Header fftwd3frc.h
//
//
//	HEADER FILE TITLE:
//
//	C++ wrapper for FFTW single-precision (float) FFT routines of real-valued 3D arrays (forward "real-to-complex" FFT and inverse "complex-to-real" FFT)
//
//	COPYRIGHT:
//
//		TEG 2019
//
//
/*!
	\file		fftwd3frc.h
	\brief		C++ wrapper for FFTW single-precision (float) FFT routines of real-valued 3D arrays
	\par		Description:
		This is C++ wrapper for FFTW single-precision (float) FFT routines of real-valued 3D arrays (forward "real-to-complex" FFT and inverse "complex-to-real" FFT).

	NOTE: this 3D real-to-complex FFT class always allocates internal storage and performs the operations on it, rather than on an externally provided array.
	This design has been selected, particularly because the structure of the internal real 3D array is a little complicated (it always contains some padding in
	the last dimension - see the FFTW library manual). It is therefore easier to use Set / Get functions to move the data in and out this class for processing.
	Because 3D FFT tends to be relatively slow, the movement of data in and out of this class should not affect the overall performance too much.
	As one consequence of this design, the "out-of-place" mode with separate internal arrays for "before" and "after" the transform rarely makes sense, 
	and so bInPlace should almost always be true.
	When bMeasure is equal to ESTIMATE, it takes a lot longer to create an object/plan, but does not seem to improve the FFT performance at all.
	
	NOTE: this FFT works in principle for any array dimensions, however: (a) it is the fastest when the dimensions are powers of 2; (b) when a dimension is odd,
	it can create problems for any "shuffle" like operation, which uses dimension/2 index - see e.g. Get/Set functions

	NOTE: in contrast to the fftwd3drc class, this class has several template functions, which gives it the ability to work not only on float, 
	but also on double 3D arrays. This can be useful when there is not enough RAM to create a copy of a double 3D array for FFTW processing
	
*/
#if !defined FFTWD3FRC_H
#define FFTWD3FRC_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//

#include "fftw3.h"
#include "XAHWave.h"
#include "XArray3D.h"

class Fftwd3frc
{
	// Enumerators
	// Structures
	// Constructors
public:
	//! Default constructor
	//Fftwd3frc() { if (!IsEmpty()) Cleanup(); uiflag = FFTW_ESTIMATE;  }
	//! Constructor allocating internal storage and creating "FFT plans"
	Fftwd3frc(int nz1, int ny1, int nx1, bool bRealSpace1 = true, bool bInPlace1 = true, bool bMeasure = false);

private:
	//! Copy constructor (declared private to prohibit copying)
	Fftwd3frc(const Fftwd3frc& other) {  }

public:
	//! Destructor, destroys the plans and deallocates internal storage
	~Fftwd3frc() { if (!IsEmpty()) Cleanup(); }

	// Operators
private:
	//! Assignment operator (declared private to prohibit copying)
	Fftwd3frc& operator=(const Fftwd3frc& other) { return *this; }

	// Attributes
public:

	// Operations
public:
	// Member functions
	// Get full size of the real array
	inline size_t GetNr() const { return size_t(nz) * size_t(ny) * size_t(nx); }
	// Get the last dimension of the complex array
	inline int GetNx2() const { return int(nx / 2 + 1); }
	// Get full size of the complex array
	inline size_t GetNc() const { return size_t(nz) * size_t(ny) * size_t(GetNx2()); }
	// Checks if the object is empty
	inline bool IsEmpty() const { return (nz == 0 && ny == 0 && nx == 0 && pin == 0 && pout == 0 && aplan == 0 && bplan == 0); }
	// Empties the object: may be called if one wants to release the memory before the destructor is called
	void Cleanup();
	// Get/Set arrays
	// Copies the internal real array into an external one of the same size (only works if bRealSpace)
	template <class T> void GetRealXArray3D(xar::XArray3D<T>& aaa, bool bShuffle) const;
	// Copies an external array of the same size into the internal real array (only works if bRealSpace)
	template <class T> void SetRealXArray3D(const xar::XArray3D<T>& aaa, bool bShuffle);
	// Copies the internal complex "half-sized" array into a full-size external complex Hermitian 3D array (only works if !bRealSpace)
	// NOTE: the output array is supposed to be "Nx2-sized", and this function fills it according to the Hermitian symmetry
	template <class T> void GetComplexXArray3D(xar::XArray3D< std::complex<T> >& aaa, bool bShuffle) const;
	// Sets the internal complex "half-sized" array from a full-size external complex Hermitian 3D array (only works if !bRealSpace)
	// NOTE: the input array is supposed to be "full size", and this function uses only the "first half" of the input Hermitian 3D array
	template <class T> void SetComplexXArray3D(const xar::XArray3D< std::complex<T> >& aaa, bool bShuffle);
	// Forward FFT (only works if bRealSpace)
	void ForwardFFT();
	// Inverse FFT (only works if !bRealSpace)
	void InverseFFT();
	// Print internal arrays
	//void PrintRealArray(const char* message);
	//void PrintComplexArray(const char* message);
	// Regularized inverse minus-Laplacian (only works if bRealSpace and bInPlace)
	template <class T> void InverseMLaplacian(xar::XArray3D<T>& xa3, double alpha, double norm);
	// Regularized inverse minus-Laplacian in Fourier space (only works if bContainsData && !bRealSpace)
	void InverseMLaplacian1(double zaper, double yaper, double xaper, double alpha, double norm);
	// Complex monomorphous CTF correction on the Ewald sphere (only works if bContainsData && !bRealSpace)
	void CTFcorrection1(xar::Wavehead3D head3D, double defocdist, double q2max, double Cs3, double Cs5, double sigma, double alpha, bool bESCC);
	// Convolution with a Gaussian (only works if bRealSpace and bInPlace)
	template <class T> void GaussFilter(xar::XArray3D<T>& xa3, double sigma);
	// Convolution with a Gaussian in Fourier space (only works if bContainsData && !bRealSpace, sigma is the REAL SPACE std)
	void GaussFilter1(double zaper, double yaper, double xaper, double sigma);

private:
	// Member variables
	int nz; // z dimension (= Dim1 for XArray3D, = n0 for fftw)
	int ny; // y dimension (= Dim2 for XArray3D, = n1 for fftw)
	int nx; // x dimension (= Dim3 for XArray3D, = n2 for fftw)
	unsigned int uiflag; // FFTW_MEASURE or FFTW_ESTIMATE

	float* pin; // pointer to the real array
	fftwf_complex* pout; // pointer to the complex array

	bool bInPlace; // indicates if in-place or out-of-place FFT is used (should be true almost always)
	bool bContainsData; // indicates if the object has been populated with meaningful data
	bool bRealSpace; // indicates if the current state corresponds to the real 3D array (real space) or complex 3D array (Fourier space)

	fftwf_plan aplan; // forward plan
	fftwf_plan bplan; // inverse plan
};


#endif	// FFTWD3FRC_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//