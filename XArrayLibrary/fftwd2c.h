//Header fftwd2c.h
//
//
//	HEADER FILE TITLE:
//
//	C++ wrapper for FFTW double-precision (double) in-place FFT routines of complex 2D arrays (forward "complex-to-complex" FFT and inverse "complex-to-complex" FFT)
//
//	COPYRIGHT:
//
//		TEG 2022
//
//
/*!
	\file		fftwd2c.h
	\brief		C++ wrapper for FFTW double-precision (double) in-place FFT routines of complex 2D arrays
	\par		Description:
		This is C++ wrapper for FFTW double-precision (double) in-place FFT routines of complex 2D arrays (forward "complex-to-complex" FFT and inverse "complex-to-complex" FFT).

		NOTE: this 2D FFT class can allocate internal storage and perform the operations on it, rather than on an externally provided array.
		In order to optimize the execution of FFTW routines, a "plan" should be "estimated" first, which takes time and messes up the storage volume. 
		When 2D FFT is called repeatedly, it is cheaper to reuse a single "optimized" plan and the same	internal stordage, 
		and shift external data in and out for processing, ratherthan use external array storage, renewing the "plans" every time.
*/
#if !defined FFTWD2C_H
#define FFTWD2C_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//

#include "fftw3.h"

#include "IXAHWave.h"
#include "XAHWave.h"
#include "XArray2D.h"

class Fftwd2c
{
	// Enumerators
	// Structures
	// Constructors
public:
	//! Default constructor
	//Fftwd2c() { if (!IsEmpty()) Cleanup(); uiflag = FFTW_ESTIMATE;  }
	//! Constructor allocating internal storage and creating "FFT plans"
	Fftwd2c(int ny1, int nx1, bool bMeasure = true);
	//! Constructor using external storage and creating "FFT plans"
	Fftwd2c(xar::XArray2D<xar::dcomplex>& xarc2D);
private:
	//! Copy constructor (declared private to prohibit copying)
	Fftwd2c(const Fftwd2c& other) {  }
public:
	//! Destructor, destroys the plans and deallocates internal storage
	~Fftwd2c() { if (!IsEmpty()) Cleanup(); }

	// Operators
private:
	//! Assignment operator (declared private to prohibit copying)
	Fftwd2c& operator=(const Fftwd2c& other) { return *this; }

	// Attributes
public:

	// Operations
public:
	// Member functions
	// Get the first dimension the complex array
	inline int GetDim1() const { return ny; }
	// Get the second dimension the complex array
	inline int GetDim2() const { return nx; }
	// Get full size of the complex array
	inline int GetNp() const { return ny * nx; }
	// Check if the object is empty
	inline bool IsEmpty() const { return (ny == 0 && nx == 0 && parr == 0 && aplan == 0 && bplan == 0); }
	// Empties the object: may be called if one wants to release the memory before the destructor is called
	void Cleanup();
	// Get/Set arrays
	// Copy the internal 2D array into an external one of the same size (only works if !bExtStorage)
	void GetXArray2D(xar::XArray2D<xar::dcomplex>& aaa) const;
	// Copy an external 2D array of the same size into the internal one (only works if !bExtStorage)
	void SetXArray2D(const xar::XArray2D<xar::dcomplex>& aaa);
	// Shuffle internal array (this is equivalent to xar::XArray2DFFT<double>::Shuffle())
	void Shuffle();
	// Forward FFT (NOTE: this corresponds to inverse FFT in OouraFFT and XArray2DFFT<double>::FFT)
	void ForwardFFT() { fftw_execute(aplan); }
	// Forward transform using external arrays with an alrady existing plan (speeds up repeated transforms with the same dimensions)
	void ForwardFFT(xar::XArray2D<xar::dcomplex>& xarc2D);
	// Inverse FFT (NOTE: this corresponds to forward FFT in OouraFFT and XArray2DFFT<double>::FFT)
	void InverseFFT() { fftw_execute(bplan); }
	// Inverse transform using external arrays with an already existing plan (speeds up repeated transforms with the same dimensions)
	void InverseFFT(xar::XArray2D<xar::dcomplex>& xarc2D);
	// Regularized inverse minus-Laplacian, leaving the input array in the unshuffled form in the Fourier space
	// NOTE: compared to TIE-Hom, alpha = 1 / [(delta/beta) * PI * dz * wl]
	void InverseMLaplacian1(xar::XArray2D<xar::dcomplex>& xarc2D, double alpha, double norm);
	// Print internal array
	void PrintComplexArray(const char* message);

private:
	// Member variables
	int ny; // y dimension (= Dim1 for XArray2D, = n0 for fftw)
	int nx; // x dimension (= Dim2 for XArray2D, = n1 for fftw)
	unsigned int uiflag; // FFTW_MEASURE or FFTW_ESTIMATE
	bool bExtStorage; // indicates if the array allocation is handled internally or externally

	fftw_complex* parr; // pointer to the complex array

	fftw_plan aplan; // forward plan
	fftw_plan bplan; // inverse plan
};


#endif	// FFTWD2C_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//