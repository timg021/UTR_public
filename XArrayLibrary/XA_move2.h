//Header XA_move2.h
//
//
//	HEADER FILE TITLE:
//
//		Two-dimensional geometrical transformations
//
/*!
	\file		XA_move2.h
	\brief		Two-dimensional geometrical transformations
	\par		Description:
		This class implements simple 'geometrical' transformations of
		XArray2D<T> objects
*/
#pragma once
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_head2.h"
#include "XA_fftr2.h"

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
//Class XArray2DMove
//
//	Two-dimensional geometrical transformations
//
/*!
	\brief		Two-dimensional geometrical transformations
	\par		Description:
				This class template defines a 'wrapper' around the XArray2D<T> object
				on which it operates; it contains functions implementing simple
				geometrical transformations of the 'wrapped' XArray2D<T> object
	\remarks    If a Wavehead2D is attached to the XArray2D<T> object, it is transformed too
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
*/
	template <class T> class XArray2DMove
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray2DMove(XArray2D<T>& rXAr2D) : m_rXArray2D(rXAr2D) {}
	protected:
		//! Copy constructor (declared protected to prohibit explicit copying)
		XArray2DMove(const XArray2DMove<T>& rCopy) : m_rXArray2D(rCopy.m_rXArray2D) {}
	public:
		//! Destructor 
		~XArray2DMove() {}
	
	// Operators
	protected:
		//! Assignment (declared protected to prohibit explicit copying)
		void operator=(const XArray2DMove<T>& rCopy);

	// Attributes
	public:
		//! Returns a reference to the non-modifiable 'wrapped' XArray2D<T> object
		const XArray2D<T>& GetBaseObject() const { return m_rXArray2D; }
		//! Returns a reference to the 'wrapped' XArray2D<T> object
		XArray2D<T>& GetBaseObject() { return m_rXArray2D; }

	// Operations
	public:
		//! Transposes the matrix
		void Transpose(void);
		//! Mirror-reflects the rows
		void FlipX(void);
		//! Trims elements from the edges of the array
		void Trim(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight);
		//! Adds new elements with a given value at the edges of the array
		void Pad(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tPadVal);
		//! Adds new elements at the edges of the array reflecting with respect to the boundaries
		void PadMirror(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight);
		//! Adds new elements with a given value at the edges of the array up to the nearest integer powers of 2
		void Pad2N(T tPadVal);
		//! Replaces elements at the edges of the array using a given value 
		void Mask(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tMaskVal);
		//! Replaces elements at the edges of the array using a given value, with a smooth transition
		void MaskSmooth(index_t iPix, double dSteepness, T tMaskVal);
		//! Moves elements of the array and fills the vacated positions with a given value
		void Move(long lngMoveYPoints, long lngMoveXPoints, T tFillVal);
		//! Fills a rectangular subarray with a given value
		void FillRect(index_t LowDim1, index_t HighDim1, index_t LowDim2, index_t HighDim2, T tFillVal);
		//! Imposes a Gaussian envelope over the existing elements in the array
		void GaussEnvelope(double dblSigmaY, double dblSigmaX);
		//! Imposes a Gaussian envelope (in the Fourier space) over the existing elements in the array
		void GaussEnvelopeFourier(double dblSigmaY, double dblSigmaX);
		//!	Imposes an Exponential envelope (in the Fourier space) over the existing elements in the array
		void ExponentialEnvelopeFourier(double dblSigma);
		//!	Imposes an Epanechnikov envelope over the existing elements in the array
		void EpanechnikovEnvelope(double dblWidthY, double dblWidthX);
		//!	Imposes an axially symmetrical Epanechnikov envelope over the existing elements in the array
		void SymEpanechnikovEnvelope(double dblWidth);
		//!	Filters the array using a mask
		void MaskFilter(const XArray2D<T>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX);
		//!	Filters the array using a mask
		void MaskFilterFourier(const XArray2D<T>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX);
		//!	Calculates partial derivatives of the array
		void Derivative(index_t iOrderY, index_t iOrderX);
		//! Calculates integral moments of the array
		T Moment(index_t iOrderY, index_t iOrderX, double dblCentreY, double dblCentreX) const;
		//! Divides the array by the average over the region
		void DoRegionNormalisation(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2);
		//! Divides each row of the array by the corresponding average value over its sub-row
		void DoRowWiseNormalisation(index_t iBeginDim2, index_t iEndDim2, int op = 0);
		//!	Converts the array from Cartesian to polar coordinates
		void Cart2Polar(bool& bOdd);
		//!	Converts the array from polar to Cartesian coordinates
		void Polar2Cart(bool bOdd);
		//!	Correlation coefficient between the 'wrapped' array and a reference array
		double CorrelationCoeff( const XArray2D<T>& rXArRef );
		//! Rebins the array using new pixel sizes and number of pixels
		void ReBin(double dblXlo, double dblXstep, index_t iNumX, double dblYlo, double dblYstep, index_t iNumY, bool bNormalize = true);
		//!	Applies modulus-gradient transform
		void ModGrad(void);
		//!	Applies modulus-gradient transform using right differences
		void ModGradRight(void);
		//! Perform a 4-connected erosion
		void Erode4(T tHigh = 1);
		//! Perform an 8-connected erosion
		void Erode8(T tHigh = 1);
		//! Perform a 4-connected dilation
		void Dilate4(T tHigh = 1);
		//! Perform an 8-connected dilation
		void Dilate8(T tHigh = 1);

		// The following functions are not designed for user interface
		// Calculates centred integral moments of the PSF (used in Taylor deconvolution)
		void PSFMoments(index_t iMaxOrder, XArray2D<T>& rMoments) const;
		// Calculates coefficients of the inverse Taylor series (used in Taylor deconvolution)
		static void InverseSeries(const XArray2D<T>& rAcoeff, XArray2D<T>& rBcoeff);

	private:
		//! Rebins a 1D array into a new one
		void ReBin1D(std::vector<double>& x, std::vector<T>& s, std::vector<double>& y, std::vector<T>& d, bool Initialize = true);

	// Overridables

	// Implementation
	private:
	// Member variables	
		//! Reference to the XArray2D<T> object that is being operated upon
		XArray2D<T>& m_rXArray2D;
	};


	//---------------------------------------------------------------------------
	//	TEMPLATE MEMBER DEFINITIONS
	//

	//! Assignment  (declared protected to prohibit explicit copying)
	template <class T> void XArray2DMove<T>::operator=(const XArray2DMove<T>& rCopy)
	{ 
		if (this == &rCopy) 
			return; 
		else 
			m_rXArray2D = rCopy.m_rXArray2D;
	}

	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Transpose
	//
	//	Transposes the matrix
	//
	/*!
		\brief		Transposes the matrix
		\param		None
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function transposes the 2D array (matrix), i.e. swaps the X- and Y-dimensions;
			it also calls the Transpose function on the IXAHWave2D head pointer if present
	*/	
	template <class T> void XArray2DMove<T>::Transpose(void)
	{
		XArray2D<T> xarTemp(m_rXArray2D.GetDim2(), m_rXArray2D.GetDim1());

		for(index_t j = 0; j < m_rXArray2D.GetDim2(); j++)
			for(index_t i = 0; i < m_rXArray2D.GetDim1(); i++)
				xarTemp[j][i] = m_rXArray2D[i][j];

		// set the same head
		xarTemp.SetHeadPtr(m_rXArray2D.GetHeadPtr() ? m_rXArray2D.GetHeadPtr()->Clone() : 0);
	
		// transform the head if it can be done
		if (GetIXAHWave2D(xarTemp)) 
			GetIXAHWave2D(xarTemp)->Transpose();

		m_rXArray2D.Swap(xarTemp);
	}

	//---------------------------------------------------------------------------
		//Function XArray2DMove<T>::FlipX
		//
		//	Flips the matrix along the second dimension
		//
		/*!
			\brief		Flips the matrix along the second dimension
			\param		None
			\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
						called from inside this function
			\return		\a None
			\par		Description:
				This function Flips the 2D array (matrix) along the second dimension, i.e. mirror-reflects its rows;
				it does not change the Wavehead2D
		*/
	template <class T> void XArray2DMove<T>::FlipX(void)
	{
		index_t nx1 = m_rXArray2D.GetDim2() - 1;
		index_t nxd2 = m_rXArray2D.GetDim2() / 2;

		for (index_t i = 0; i < m_rXArray2D.GetDim1(); i++)
			for (index_t j = 0; j < nxd2; j++)
				std::swap(m_rXArray2D[i][j], m_rXArray2D[i][nx1 - j]);
	}
	
	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Trim
	//
	//	Trims elements from the edges of the array
	//
	/*!
		\brief		Trims elements from the edges of the array
		\param		iYLeft	number of Y elements to be trimmed from the beginning
		\param		iYRight	number of Y elements to be trimmed from the end
		\param		iXLeft	number of X elements to be trimmed from the beginning
		\param		iXRight	number of X elements to be trimmed from the end
		\exception	std::invalid_argument is thrown if the number of elements to trim is either
					negative or exceeds the array size
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function 'trims' the defined number of elements from the left,	right
					top and bottom edges of the XArray2D<T> object; it also calls the Trim function
					on the IXAHWave2D head pointer if present
	*/	
	template <class T> void XArray2DMove<T>::Trim(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight)
	{
		if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0) 
			return;

		if (iXLeft + iXRight > m_rXArray2D.GetDim2()) 
			throw std::invalid_argument("invalid_argument 'iXLeft or iXRight' in XArray2DMove<T>::Trim"); 

		if (iYLeft + iYRight > m_rXArray2D.GetDim1()) 
			throw std::invalid_argument("invalid_argument 'iYLeft or iYRight' in XArray2DMove<T>::Trim"); 	

		index_t iDim2 = m_rXArray2D.GetDim2() - (iXLeft + iXRight);
		index_t iDim1 = m_rXArray2D.GetDim1() - (iYLeft + iYRight);
		XArray2D<T> xarTemp(iDim1, iDim2);
	
		for (index_t i = 0; i < iDim1; i++)
			for (index_t j = 0; j < iDim2; j++)
				xarTemp[i][j] = (m_rXArray2D)[i + iYLeft][j + iXLeft];

		// set the same head
		xarTemp.SetHeadPtr(m_rXArray2D.GetHeadPtr() ? m_rXArray2D.GetHeadPtr()->Clone() : 0);
	
		// transform the head if it can be done
		if (GetIXAHWave2D(xarTemp)) 
			GetIXAHWave2D(xarTemp)->Trim(m_rXArray2D.GetDim1(), m_rXArray2D.GetDim2(), iYLeft, iYRight, iXLeft, iXRight);

		m_rXArray2D.Swap(xarTemp);
	}

	//Function XArray2DMove<T>::Pad
	//
	//	Adds new elements with a given value at the edges of the array
	//
	/*!
		\brief		Adds new elements with a given value at the edges of the array
		\param		iYLeft	number of Y elements to be added at the beginning
		\param		iYRight	number of Y elements to be added at the end
		\param		iXLeft	number of X elements to be added at the beginning
		\param		iXRight	number of X elements to be added at the end
		\param		tPadVal	the value assigned to all new elements
		\exception	std::invalid_argument is thrown if the number of elements to pad is negative
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function 'pads' the XArray2D<T> object by adding the defined number of 
					elements at the left, right, top and bottom edges of the array;
					it also calls the Pad function on the IXAHWave2D head pointer if present
	*/	
	template <class T> void XArray2DMove<T>::Pad(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tPadVal)
	{
		if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0) 
			return;

		index_t iDim2 = m_rXArray2D.GetDim2() + (iXLeft + iXRight);
		index_t iDim1 = m_rXArray2D.GetDim1() + (iYLeft + iYRight);
		XArray2D<T> xarTemp(iDim1, iDim2, tPadVal);

		for (index_t i = 0; i < m_rXArray2D.GetDim1(); i++)
			for (index_t j = 0; j < m_rXArray2D.GetDim2(); j++)
				xarTemp[i + iYLeft][j + iXLeft] = m_rXArray2D[i][j];

		// set the same head
		xarTemp.SetHeadPtr(m_rXArray2D.GetHeadPtr() ? m_rXArray2D.GetHeadPtr()->Clone() : 0);
	
		// transform the head if it can be done
		if (GetIXAHWave2D(xarTemp)) 
			GetIXAHWave2D(xarTemp)->Pad(m_rXArray2D.GetDim1(), m_rXArray2D.GetDim2(), iYLeft, iYRight, iXLeft, iXRight);

		m_rXArray2D.Swap(xarTemp);
	}

	//Function XArray2DMove<T>::PadMirror
	//
	//	Adds new elements at the edges of the array reflecting with respect to the boundaries
	//
	/*!
		\brief		Adds new elements at the edges of the array reflecting with respect to the boundaries
		\param		iYLeft	number of Y elements to be added at the beginning
		\param		iYRight	number of Y elements to be added at the end
		\param		iXLeft	number of X elements to be added at the beginning
		\param		iXRight	number of X elements to be added at the end
		\exception	std::invalid_argument is thrown if the number of elements to pad is negative or exceeds the array size
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function 'pads' the XArray2D<T> object by adding the defined number of 
					elements at the left, right, top and bottom edges of the array mirror-reflecting existing elements
					with respect to the array boundaries; it also calls the Pad function on the IXAHWave2D head pointer 
					if present
	*/	
	template <class T> void XArray2DMove<T>::PadMirror(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight)
	{
		if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0) 
			return;

		if (iXLeft >= m_rXArray2D.GetDim2() || iXRight >= m_rXArray2D.GetDim2()) 
			throw std::invalid_argument("invalid_argument 'iXLeft or iXRight' in XArray2DMove<T>::PadMirror (too many X points to pad)"); 	

		if (iYLeft >= m_rXArray2D.GetDim1() || iYRight >= m_rXArray2D.GetDim1()) 
			throw std::invalid_argument("invalid_argument 'iYLeft or iYRight' in XArray2DMove<T>::PadMirror (too many Y points to pad)"); 	

		index_t iDim2 = m_rXArray2D.GetDim2() + (iXLeft + iXRight);
		index_t iDim1 = m_rXArray2D.GetDim1() + (iYLeft + iYRight);
		XArray2D<T> xarTemp(iDim1, iDim2);

		// fill the centre area
		index_t i, j;
		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
				xarTemp[i + iYLeft][j + iXLeft] = m_rXArray2D[i][j];

		// mirror 'reflection' into new points at left boundary
		index_t iTemp = 2 * iXLeft;
		for (j = 0; j < iXLeft; j++)
			for (i = 0; i < iDim1; i++)
				xarTemp[i][j] = xarTemp[i][iTemp - j];
		// mirror 'reflection' into new points at top boundary
		iTemp = 2 * iYLeft;
		for (i = 0; i < iYLeft; i++)
			for (j = 0; j < iDim2; j++)
				xarTemp[i][j] = xarTemp[iTemp - i][j];
		// mirror 'reflection' into new points at right boundary
		iTemp = 2 * (iDim2 - 1 - iXRight);
		for (j = iDim2 - iXRight; j < iDim2; j++)
			for (i = 0; i < iDim1; i++)
				xarTemp[i][j] = xarTemp[i][iTemp - j];
		// mirror 'reflection' into new points at bottom boundary
		iTemp = 2 * (iDim1 - 1 - iYRight);
		for (i = iDim1 - iYRight; i < iDim1; i++)
			for (j = 0; j < iDim2; j++)
				xarTemp[i][j] = xarTemp[iTemp - i][j];

		// set the same head
		xarTemp.SetHeadPtr(m_rXArray2D.GetHeadPtr() ? m_rXArray2D.GetHeadPtr()->Clone() : 0);

		// transform the head if it can be done
		if (GetIXAHWave2D(xarTemp)) 
			GetIXAHWave2D(xarTemp)->Pad(m_rXArray2D.GetDim1(), m_rXArray2D.GetDim2(), iYLeft, iYRight, iXLeft, iXRight);

		m_rXArray2D.Swap(xarTemp);
	}


	//Function XArray2DMove<T>::Pad2N
	//
	//	Adds new elements with a given value at the edges of the array up to the nearest integer powers of 2
	//
	/*!
		\brief		Adds new elements with a given value at the edges of the array up to the nearest integer powers of 2
		\param		tPadVal	the value assigned to all new elements
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function 'pads' the XArray2D<T> object by new elements with the given value
					at the left, right, top and bottom edges of the array increasing its dimensions to
					the nearest	larger powers of 2; it also calls the Pad function with the appropriate
					parameters on the IXAHWave2D head pointer if present
	*/	
	template <class T> void XArray2DMove<T>::Pad2N(T tPadVal)
	{
		index_t i = 2;
		while (i < m_rXArray2D.GetDim1()) i *= 2;
		index_t nypad = i - m_rXArray2D.GetDim1();
		index_t j = 2;
		while (j < m_rXArray2D.GetDim2()) j *= 2;
		index_t nxpad = j - m_rXArray2D.GetDim2();

		Pad(nypad / 2, nypad - nypad / 2, nxpad / 2, nxpad - nxpad / 2, tPadVal);
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Mask
	//
	//	Replaces elements at the edges of the array using a given value 
	//
	/*!
		\brief		Replaces elements at the edges of the array using a given value 
		\param		iYLeft	number of Y elements to be replaced at the beginning
		\param		iYRight	number of Y elements to be replaced at the end
		\param		iXLeft	number of X elements to be replaced at the beginning
		\param		iXRight	number of X elements to be replaced at the end
		\param		tMaskVal	the value to be assigned to all masked elements
		\exception	std::invalid_argument is thrown if the number of elements to mask is negative
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function 'masks' the defined number of elements at the left, right
					top and bottom edges of the XArray2D<T> object by assigning them a given
					value; the head is not affected
	*/	
	template <class T> void XArray2DMove<T>::Mask(index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tMaskVal)
	// !!! differs in behavior from the old VO-implementation
	{
		if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0) 
			return;

		if (iXLeft + iXRight > m_rXArray2D.GetDim2() || iYLeft + iYRight > m_rXArray2D.GetDim1())
			m_rXArray2D.Fill(tMaskVal);
		else
		{
			index_t i, j;
			for (i = 0; i < iYLeft; i++)
				for (j = 0; j < m_rXArray2D.GetDim2(); j++)
					m_rXArray2D[i][j] = tMaskVal;

			for (i = m_rXArray2D.GetDim1() - iYRight; i < m_rXArray2D.GetDim1(); i++)
				for (j = 0; j < m_rXArray2D.GetDim2(); j++)
					m_rXArray2D[i][j] = tMaskVal;

			for (i = iYLeft; i < m_rXArray2D.GetDim1()-iYRight; i++)
			{
				for (j = 0; j < iXLeft; j++) 
					m_rXArray2D[i][j] = tMaskVal;

				for (j = m_rXArray2D.GetDim2() - iXRight; j < m_rXArray2D.GetDim2(); j++) 
					m_rXArray2D[i][j] = tMaskVal;
			}
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>:: MaskSmooth(index_t iPix, T tMaskVal)
	//
	//	Replaces elements at the edges of the array using a given value with a smooth transition 
	//
	/*!
		\brief		Replaces elements at the edges of the array using a given value
		\param		iPix	number of elements to be modified at the edges
		\param		dSteepness	parameter controlling the steepness of transition to tMaskVal (normally between 1 and 100)
		\param		tMaskVal	the value to be assigned to all masked elements
		\return		\a None
		\par		Description:
					This function 'masks' the defined number of elements at the left, right
					top and bottom edges of the XArray2D<T> object by modifying them towards 
					a given value, with a smooth transition; the head is not affected
	*/
	template <class T> void XArray2DMove<T>::MaskSmooth(index_t iPix, double dSteepness, T tMaskVal) 
	{
		if (iPix <= 0) return;

		//int idist, jdist;
		int idist2, jdist2;
		index_t ny = m_rXArray2D.GetDim1();
		index_t nx = m_rXArray2D.GetDim2();
		int ny2 = int(ny / 2), nx2 = int(nx / 2);
		//int ny1 = int(ny - 1), nx1 = int(nx - 1);
		int iPix2 = int(0.5 * iPix + 0.5);
		int iRad = min(nx2 - iPix2, ny2 - iPix2); // radius length of the "transition" line
		double aPI = 1.0 / PI;
		double fact;

		if (iPix * 2 > m_rXArray2D.GetDim2() || iPix * 2 > m_rXArray2D.GetDim1())
			m_rXArray2D.Fill(tMaskVal);
		else
		{
			for (int j = 0; j < ny; j++)
			{
				//jdist = min(j, ny1 - j);
				jdist2 = (j - ny2) * (j - ny2);
				for (int i = 0; i < nx; i++)
				{
					//idist = min(jdist, min(i, nx1 - i));
					//if (idist > iPix) continue;
					//fact = 0.5 + aPI * atan(dSteepness * (idist - iPix2));
					idist2 = int(0.5 + sqrt(jdist2 + (i - nx2) * (i - nx2)));
					fact = 0.5 + aPI * atan(dSteepness * (idist2 - iRad));
					m_rXArray2D[j][i] = m_rXArray2D[j][i] * T(1.0 - fact) + tMaskVal * T(fact);
				}
			}
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Move
	//
	//	Moves elements of the array and fills the vacated positions with a given value
	//
	/*!
		\brief		Moves elements of the array and fills the vacated positions with a given value
		\param		lngMoveYPoints	the span of the Y translation; negative values correspond to translations upwards
		\param		lngMoveXPoints	the span of the X translation; negative values correspond to translations to the left
		\param		tFillVal	the value to be assigned to all vacated elements
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function 'moves' the the XArray2D<T> object by translating all the
					array elements by the defined number of positions in Y and X directions;
					it does not call any functions on the head. If the lngMoveYPoints is negative,
					the translation	is upwards, otherwise the array is translated downwards.
					If the lngMoveXPoints is negative, the translation	is to the left, otherwise
					the array is translated to the right. All vacated positions are filled with
					the defined value. The translation is performed	'in-place'
	*/	 
	template <class T> void XArray2DMove<T>::Move(long lngMoveYPoints, long lngMoveXPoints, T tFillVal)
	{
		if (lngMoveYPoints == 0 && lngMoveXPoints == 0) 
			return;

		// head does not move as the 'viewport' remains constant	

		if (fabs(lngMoveYPoints) >= m_rXArray2D.GetDim1() || fabs(lngMoveXPoints) >= m_rXArray2D.GetDim2())
		{
			m_rXArray2D.Fill(tFillVal);
			return;
		}

		index_t i, j;

		if (lngMoveYPoints>=0 && lngMoveXPoints>=0) 
		{
			// WARNING: cannot have i>=0 end-of-loop condition for the unsigned!
			for (i=index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i>0; i--) 
				for (j=index_t(m_rXArray2D.GetDim2()-lngMoveXPoints); j>0; j--)
					(m_rXArray2D)[i-1+lngMoveYPoints][j-1+lngMoveXPoints] = m_rXArray2D[i-1][j-1];
			for (i=0; i<index_t(lngMoveYPoints); i++) 
				for (j=0; j<m_rXArray2D.GetDim2(); j++) 
					m_rXArray2D[i][j] = tFillVal;
			for (i=index_t(lngMoveYPoints); i<m_rXArray2D.GetDim1(); i++) 
				for (j=0; j<index_t(lngMoveXPoints); j++)
					m_rXArray2D[i][j] = tFillVal;				
		}
		else if (lngMoveYPoints>=0 && lngMoveXPoints<0) 
		{
			lngMoveXPoints = -lngMoveXPoints;
			for (i=index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i>0; i--)
				for (j=0; j<index_t(m_rXArray2D.GetDim2()-lngMoveXPoints); j++)
					(m_rXArray2D)[i-1+index_t(lngMoveYPoints)][j] = (m_rXArray2D)[i-1][j+index_t(lngMoveXPoints)];
			for (i=0; i<index_t(lngMoveYPoints); i++) 
				for (j=0; j<m_rXArray2D.GetDim2(); j++) 
					m_rXArray2D[i][j] = tFillVal;
			for (i=lngMoveYPoints; i<m_rXArray2D.GetDim1(); i++) 
				for (j=index_t(m_rXArray2D.GetDim2()-lngMoveXPoints); j<m_rXArray2D.GetDim2(); j++) 
					m_rXArray2D[i][j] = tFillVal;
		}
		else if (lngMoveYPoints<0 && lngMoveXPoints>=0) 
		{
			lngMoveYPoints = -lngMoveYPoints;
			for (i=0; i<index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i++)
				for (j=index_t(m_rXArray2D.GetDim2()-lngMoveXPoints); j>0; j--)
					(m_rXArray2D)[i][j-1+index_t(lngMoveXPoints)] = (m_rXArray2D)[i+index_t(lngMoveYPoints)][j-1];
			for (i=index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i<m_rXArray2D.GetDim1(); i++) 
				for (j=0; j<m_rXArray2D.GetDim2(); j++) 
					m_rXArray2D[i][j] = tFillVal;
			for (i=0; i<index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i++) 
				for (j=0; j<index_t(lngMoveXPoints); j++) 
					m_rXArray2D[i][j] = tFillVal;
		}
		else if (lngMoveYPoints<0 && lngMoveXPoints<0) 
		{
			lngMoveYPoints = -lngMoveYPoints; lngMoveXPoints = -lngMoveXPoints;
			for (i=0; i<index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i++)
				for (j=0; j<index_t(m_rXArray2D.GetDim2()-lngMoveXPoints); j++)
					m_rXArray2D[i][j] = (m_rXArray2D)[i+index_t(lngMoveYPoints)][j+index_t(lngMoveXPoints)];
			for (i=index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i<m_rXArray2D.GetDim1(); i++) 
				for (j=0; j<m_rXArray2D.GetDim2(); j++)
					m_rXArray2D[i][j] = tFillVal;
			for (i=0; i<index_t(m_rXArray2D.GetDim1()-lngMoveYPoints); i++) 
				for (j=index_t(m_rXArray2D.GetDim2()-lngMoveXPoints); j<m_rXArray2D.GetDim2(); j++)
					m_rXArray2D[i][j] = tFillVal;
		}
	}

	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::FillRect
	//
	//	Fills a rectangular subarray with a given value
	//
	/*!
		\brief		Fills a rectangular subarray with a given value
		\param		LowDim1	lower Y boundary of the rectangle to be filled
		\param		HighDim1	upper Y boundary of the rectangle to be filled
		\param		LowDim2	lower X boundary of the rectangle to be filled
		\param		HighDim2	upper X boundary of the rectangle to be filled
		\param		tFillVal	the value to be assigned to all elements inside the rectangle
		\exception	std::invalid_argument is thrown if any rectangle boundary is negative or
					exceeds the	corresponding array boundary, or if lower rectangle boundary is
					larger or equal to the higher rectangle boundary
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function fills a defined rectangular region inside a 2D array with a
					given value. Lower boundaries are included in the filled rectangle, while the
					upper boundaries are excluded; the head is not affected
	*/	
	template <class T> void XArray2DMove<T>::FillRect(index_t LowDim1, index_t HighDim1, index_t LowDim2, index_t HighDim2, T tFillVal)
	{
		index_t i, j;

		if (HighDim1 > m_rXArray2D.GetDim1() || LowDim1 >= HighDim1)
			throw std::invalid_argument("invalid_argument 'LowDim1 or HighDim1' in XArray2DMove<T>::FillRect"); 

		if (HighDim2 > m_rXArray2D.GetDim2() || LowDim2 >= HighDim2)
			throw std::invalid_argument("invalid_argument 'LowDim2 or HighDim2' in XArray2DMove<T>::FillRect"); 

		for (i = LowDim1; i < HighDim1; i++)
			for (j = LowDim2; j < HighDim2; j++)
				m_rXArray2D[i][j] = tFillVal; 
	}

	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::GaussEnvelope
	//
	//	Imposes a Gaussian envelope over the existing elements in the array
	//
	/*!
		\brief		Imposes a Gaussian envelope over the existing elements in the array
		\param		dblSigmaY	standard deviation of the Gaussian in the Y direction
		\param		dblSigmaX	standard deviation of the Gaussian in the X direction
		\exception	std::invalid_argument is thrown if dblSigmaX or dblSigmaY are not positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function multiplies the existing element values in an XArray2D<T> object
					by a Gaussian distribution with the mean equal to the centre of the array and 
					the defined standard deviations in X and Y directions (i.e. multiplies array
					data by the function P(x,y) = exp{-0.5*(x^2/dblSigmaX^2+y^2/dblSigmaY^2))} ).
					The head is not affected
	*/	
	template <class T> void XArray2DMove<T>::GaussEnvelope(double dblSigmaY, double dblSigmaX)
	{
		if (dblSigmaX <= 0 || dblSigmaY <= 0)
				throw std::invalid_argument("invalid_argument 'dblSigmaX or dblSigmaY' in XArray2DMove<T>::GaussEnvelope (standard deviation must be positive)"); 

		index_t i, j;
		double yi, xj;
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double yst = xar::GetYStep(m_rXArray2D);
		double xst = xar::GetXStep(m_rXArray2D);
		double adblSigmaY = -0.5 * yst * yst / (dblSigmaY * dblSigmaY), y2adblSigmaY;
		double adblSigmaX = -0.5 * xst * xst / (dblSigmaX * dblSigmaX), x2adblSigmaX;
		double bsigma = 1.0; 
		// The normalization with bsigma = 1.0 preserves the average if applied in Fourier space.
		// The following normalization leads to preservation of the average after convolution
		// if applied in real space: bsigma = 1.0 / (tPI * dblSigmaX * dblSigmaY) * xst * yst.

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			yi = i - dblCentreY;
			y2adblSigmaY = yi * yi * adblSigmaY;
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				xj = j - dblCentreX;
				x2adblSigmaX = xj * xj * adblSigmaX;
				m_rXArray2D[i][j] *= T(bsigma * exp(x2adblSigmaX + y2adblSigmaY));
			}
		}
	}

	// This is a specialization of the function xar::XArray2DMove<T>::GaussEnvelope for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::GaussEnvelope(double dblSigmaY, double dblSigmaX)
	{
		index_t i, j;
		double yi, xj;
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double yst = xar::GetYStep(m_rXArray2D);
		double xst = xar::GetXStep(m_rXArray2D);
		double adblSigmaY = -0.5 * yst * yst / (dblSigmaY * dblSigmaY), y2adblSigmaY;
		double adblSigmaX = -0.5 * xst * xst / (dblSigmaX * dblSigmaX), x2adblSigmaX;
		double bsigma = 1.0; 
		// The normalization with bsigma = 1.0 preserves the average if applied in Fourier space.
		// The following normalization leads to preservation of the average after convolution
		// if applied in real space: bsigma = 1.0 / (tPI * dblSigmaX * dblSigmaY) * xst * yst.

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			yi = i - dblCentreY;
			y2adblSigmaY = yi * yi * adblSigmaY;
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				xj = j - dblCentreX;
				x2adblSigmaX = xj * xj * adblSigmaX;
				m_rXArray2D[i][j] *= float(bsigma * exp(x2adblSigmaX + y2adblSigmaY));
			}
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::GaussEnvelopeFourier
	//
	//	Imposes a Gaussian envelope (in the Fourier space) over the existing elements in the array
	//
	/*!
		\brief		Imposes a Gaussian envelope over the existing elements in the array
		\param		dblSigmaY	standard deviation of the Gaussian in the Y direction
		\param		dblSigmaX	standard deviation of the Gaussian in the X direction
		\exception	std::invalid_argument is thrown if dblSigmaX or dblSigmaY are not positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function multiplies the existing element values in an XArray2D<T> object
					by a Gaussian distribution with the mean equal to the centre of the array and 
					the defined standard deviations in X and Y directions (i.e. multiplies array
					data by the function P(u,v) = exp{-2*(u^2*dblSigmaX^2+v^2*dblSigmaY^2))} ).
					The head is not affected
	*/
	template <class T> void XArray2DMove<T>::GaussEnvelopeFourier(double dblSigmaY, double dblSigmaX)
	{
		if (dblSigmaX < 0 || dblSigmaY < 0)
				throw std::invalid_argument("invalid_argument 'dblSigmaX or dblSigmaY' in XArray2DMove<T>::GaussEnvelopeFourier (standard deviation must be non-negative)"); 

		index_t i, j;
		double u, v;
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double vst = xar::GetYStep(m_rXArray2D);
		double ust = xar::GetXStep(m_rXArray2D);
		double adblSigmaY = -2. * PI * PI * vst * vst * dblSigmaY * dblSigmaY, v2adblSigmaY;
		double adblSigmaX = -2. * PI * PI * ust * ust * dblSigmaX * dblSigmaX, u2adblSigmaX;

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			v = i - dblCentreY;
			v2adblSigmaY = v * v * adblSigmaY;
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				u = j - dblCentreX;
				u2adblSigmaX = u * u * adblSigmaX;
				m_rXArray2D[i][j] *= T(exp(u2adblSigmaX + v2adblSigmaY));
			}
		}
	}

	// This is a specialization of the function xar::XArray2DMove<T>::GaussEnvelopeFourier for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::GaussEnvelopeFourier(double dblSigmaY, double dblSigmaX)
	{
		index_t i, j;
		double u, v;
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double vst = xar::GetYStep(m_rXArray2D);
		double ust = xar::GetXStep(m_rXArray2D);
		double adblSigmaY = -2. * PI * PI * vst * vst * dblSigmaY * dblSigmaY, v2adblSigmaY;
		double adblSigmaX = -2. * PI * PI * ust * ust * dblSigmaX * dblSigmaX, u2adblSigmaX;

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			v = i - dblCentreY;
			v2adblSigmaY = v * v * adblSigmaY;
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				u = j - dblCentreX;
				u2adblSigmaX = u * u * adblSigmaX;
				m_rXArray2D[i][j] *= float(exp(u2adblSigmaX + v2adblSigmaY));
			}
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::ExponentialEnvelopeFourier
	//
	//	Imposes an Exponential envelope (in the Fourier space) over the existing elements in the array
	//
	/*!
		\brief		Imposes an Exponential envelope over the existing elements in the array
		\param		dblSigma	half-diameter (in pixels) of the symmetric Lorentz kernel (in the real space)
		\exception	std::invalid_argument is thrown if dblSigma is not positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function multiplies the existing element values in an XArray2D<T> object
					by an Exponential distribution with the mean equal to the centre of the array and 
					the defined half-diameter (i.e. multiplies array data by the function
					P(U) = exp[-2*Pi*dblSigma*U/sqrt(c)], U = sqrt(u^2 + v^2), c = 2^(2/3) - 1).
					The head is not affected
	*/
	template <class T> void XArray2DMove<T>::ExponentialEnvelopeFourier(double dblSigma)
	{
		if (dblSigma < 0)
			throw std::invalid_argument("invalid_argument 'dblSigma' in XArray2DMove<T>::ExponentEnvelopeFourier (half-diameter must be non-negative)"); 

		index_t i, j;
		double u, v;
		double vst = xar::GetYStep(m_rXArray2D);
		double ust = xar::GetXStep(m_rXArray2D);
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double adblSigma = -2. * PI * dblSigma * ust / sqrt(std::pow(2., 2./3.) - 1.), v2, UadblSigma;
		double dblRatio = vst / ust;

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			v = (i - dblCentreY) * dblRatio;
			v2 = v * v;
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				u = j - dblCentreX;
				UadblSigma = sqrt(v2 + u * u) * adblSigma;
				m_rXArray2D[i][j] *= T(exp(UadblSigma));
			}
		}
	}

	// This is a specialization of the function xar::XArray2DMove<T>::ExponentialEnvelopeFourier for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::ExponentialEnvelopeFourier(double dblSigma)
	{
		if (dblSigma < 0)
			throw std::invalid_argument("invalid_argument 'dblSigma' in XArray2DMove<T>::ExponentEnvelopeFourier (half-diameter must be non-negative)"); 

		index_t i, j;
		double u, v;
		double vst = xar::GetYStep(m_rXArray2D);
		double ust = xar::GetXStep(m_rXArray2D);
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double adblSigma = -2. * PI * dblSigma * ust / sqrt(std::pow(2., 2./3.) - 1.), v2, UadblSigma;
		double dblRatio = vst / ust;

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			v = (i - dblCentreY) * dblRatio;
			v2 = v * v;
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				u = j - dblCentreX;
				UadblSigma = sqrt(v2 + u * u) * adblSigma;
				m_rXArray2D[i][j] *= float(exp(UadblSigma));
			}
		}
	}

	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::EpanechnikovEnvelope
	//
	//	Imposes an Epanechnikov envelope over the existing elements in the array
	//
	/*!
		\brief		Imposes an Epanechnikov envelope over the existing elements in the array
		\param		dblWidthY	width (in pixels) of the Epanechnikov kernel in the Y direction
		\param		dblWidthX	width (in pixels) of the Epanechnikov kernel in the X direction
		\exception	std::invalid_argument is thrown if dblWidthX or dblWidthY are not positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function multiplies the existing element values in an XArray2D<T> object
					by an Epanechnikov distribution with the centre of mass equal to the centre of the array
					and the defined widths in X and Y directions (i.e. multiplies array
					data by the function P(x,y) = [1-(2x/dblWidthX)^2]+ * [1-(2y/dblWidthY)^2]+.
					The head is not affected
	*/	
	template <class T> void XArray2DMove<T>::EpanechnikovEnvelope(double dblWidthY, double dblWidthX)
	{
		if (dblWidthX <= 0 || dblWidthY <= 0)
			throw std::invalid_argument("invalid_argument 'dblWidthX or dblWidthY' in XArray2DMove<T>::EpanechnikovEnvelope (kernel width must be positive)"); 

		index_t i, j;
		double yi, xj;
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double dblRatioY = 2. / 3.  / (dblWidthY * dblWidthY), y2, Ky; //TEG: second moment of Epanechnikov equals sigma^2 / 3 (corresponding to 2*sigma^2 for Gaussians)
		double dblRatioX = 2. / 3.  /  (dblWidthX * dblWidthX), x2, Kx;

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			yi = i - dblCentreY;
			y2 = yi * yi * dblRatioY;
			Ky = (1. - y2);
			if (Ky <= 0.)
				for (j = 0; j < m_rXArray2D.GetDim2(); j++)
					m_rXArray2D[i][j] = T(0);
			else
				for (j = 0; j < m_rXArray2D.GetDim2(); j++)
				{
					xj = j - dblCentreX;
					x2 = xj * xj * dblRatioX;
					Kx = (1. - x2);
					if (Kx <= 0.) m_rXArray2D[i][j] = T(0);
					else m_rXArray2D[i][j] *= T(Ky * Kx);
				}
		}

		// The following normalization leads to preservation of the average after convolution if applied in real space
		double sum = m_rXArray2D.Norm(xar::eNormAver) * m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2();
		m_rXArray2D *= T(1./sum);
	}

	// This is a specialization of the function xar::XArray2DMove<T>::EpanechnikovEnvelope for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::EpanechnikovEnvelope(double dblWidthY, double dblWidthX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<fcomplex>::EpanechnikovEnvelope (not designed to work with complex data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::EpanechnikovEnvelope for T=dcomplex
	template<> inline void XArray2DMove<xar::dcomplex>::EpanechnikovEnvelope(double dblWidthY, double dblWidthX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<dcomplex>::EpanechnikovEnvelope (not designed to work with complex data)");
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::SymEpanechnikovEnvelope
	//
	//	Imposes an axially symmetrical Epanechnikov envelope over the existing elements in the array
	//
	/*!
		\brief		Imposes an axially symmetrical Epanechnikov envelope over the existing elements in the array
		\param		dblWidth	diameter (in pixels) of the Epanechnikov kernel
		\exception	std::invalid_argument is thrown if dblWidth is not positive
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
					This function multiplies the existing element values in an XArray2D<T> object
					by an Epanechnikov distribution with the centre of mass equal to the centre of the array
					and the defined widths in X and Y directions (i.e. multiplies array
					data by the function P(x,y) = [1-(2x/dblWidth)^2-(2y/dblWidth)^2]+.
					The head is not affected
	*/	
	template <class T> void XArray2DMove<T>::SymEpanechnikovEnvelope(double dblWidth)
	{
		if (dblWidth <= 0)
			throw std::invalid_argument("invalid_argument 'dblWidth' in XArray2DMove<T>::SymEpanechnikovEnvelope (kernel width must be positive)"); 

		index_t i, j;
		double yi, xj;
		double dblCentreY = 0.5 * m_rXArray2D.GetDim1(); // this centre position is v.important!
		double dblCentreX = 0.5 * m_rXArray2D.GetDim2(); // this centre position is v.important!
		double dblRatio = 2. / 3.  / (dblWidth * dblWidth), y2, Ky; //TEG: second moment of Epanechnikov equals sigma^2 / 3 (corresponding to 2*sigma^2 for Gaussians)
		double x2, K;

		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			yi = i - dblCentreY;
			y2 = yi * yi * dblRatio;
			Ky = 1. - y2;
			if (Ky <= 0.)
				for (j = 0; j < m_rXArray2D.GetDim2(); j++)
					m_rXArray2D[i][j] = T(0);
			else
				for (j = 0; j < m_rXArray2D.GetDim2(); j++)
				{
					xj = j - dblCentreX;
					x2 = xj * xj * dblRatio;
					K = Ky - x2;
					if (K <= 0.) m_rXArray2D[i][j] = T(0);
					else m_rXArray2D[i][j] *= T(K);
				}
		}

		// The following normalization leads to preservation of the average after convolution if applied in real space
		double sum = m_rXArray2D.Norm(xar::eNormAver) * m_rXArray2D.GetDim1() * m_rXArray2D.GetDim2();
		m_rXArray2D *= T(1./sum);
	}

	// This is a specialization of the function xar::XArray2DMove<T>::EpanechnikovEnvelope for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::SymEpanechnikovEnvelope(double dblWidth)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<fcomplex>::SymEpanechnikovEnvelope (not designed to work with complex data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::EpanechnikovEnvelope for T=dcomplex
	template<> inline void XArray2DMove<xar::dcomplex>::SymEpanechnikovEnvelope(double dblWidth)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<dcomplex>::SymEpanechnikovEnvelope (not designed to work with complex data)");
	}


	//---------------------------------------------------------------------------
	//Function XArray2DSavGol<T>::MaskFilter
	//
	//	Filters the array using a mask
	//
	/*!
		\brief		Filters the array using a mask
		\param		rXAr2DMask is the mask matrix
		\param		iMaskCentreY is the mask center Y coordinate
		\param		iMaskCentreX is the mask center X coordinate
		\exception	std::invalid_argument is thrown if input parameters are inconsistent
		\return		\a None
		\par		Description:
			This function applies a mask to filter the array or calculate its pointwise partial derivative. Uses
			mirror reflection at the boundaries. Can be used e.g. to average the array, calculate its Laplacian, etc.
	*/	
	template<class T> void XArray2DMove<T>::MaskFilter(const XArray2D<T>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		if (iMaskCentreY >= rXAr2DMask.GetDim1() || iMaskCentreX >= rXAr2DMask.GetDim2())
			throw std::invalid_argument("invalid argument iMaskCentreX or iMaskCentreY in xar::XArray2DMove<T>::MaskFilter (centre located outside the mask)");

		if (rXAr2DMask.GetDim1() > m_rXArray2D.GetDim1() || rXAr2DMask.GetDim2() > m_rXArray2D.GetDim2())
			throw std::invalid_argument("invalid argument rXAr2DMask in xar::XArray2DMove<T>::MaskFilter (mask too large for the image)");

		index_t iPointsLeft = iMaskCentreX;
		index_t iPointsRight = rXAr2DMask.GetDim2() - iMaskCentreX - 1;	
		index_t iPointsUp = iMaskCentreY;
		index_t iPointsDown = rXAr2DMask.GetDim1() - iMaskCentreY - 1;

		XArray2D<T> XArTemp(m_rXArray2D);
		XArray2DMove<T> XArTempMove(XArTemp);
		XArTempMove.PadMirror(iPointsUp, iPointsDown, iPointsLeft, iPointsRight);

		index_t i, j, k, l;
		for (i = 0; i < m_rXArray2D.GetDim1(); i++)
		{
			for (j = 0; j < m_rXArray2D.GetDim2(); j++)
			{
				T sum = 0;
				for (k = 0; k < rXAr2DMask.GetDim1(); k++)
				{
					for (l = 0; l < rXAr2DMask.GetDim2(); l++)
					{
						sum += XArTemp[i + k][j + l] * rXAr2DMask[k][l];
					}
				}
				m_rXArray2D[i][j] = sum;
			}
		}
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::MaskFilterFourier
	//
	//	Filters the array using a mask
	//
	/*!
		\brief		Filters the array using a mask
		\param		rXAr2DMask is the mask matrix
		\param		iMaskCentreY is the mask center Y coordinate
		\param		iMaskCentreX is the mask center X coordinate
		\exception	std::invalid_argument is thrown if input parameters are inconsistent
		\return		\a None
		\par		Description:
			This function applies a mask to filter the array or calculate its pointwise partial derivative. Uses
			mirror reflection at the boundaries. Can be used e.g. to average the array, calculate its Laplacian, etc.
	*/	
	template<class T> void XArray2DMove<T>::MaskFilterFourier(const XArray2D<T>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		if (iMaskCentreY >= rXAr2DMask.GetDim1() || iMaskCentreX >= rXAr2DMask.GetDim2())
			throw std::invalid_argument("invalid argument iMaskCentreX or iMaskCentreY in xar::XArray2DMove<T>::MaskFilter (centre located outside the mask)");

		if (rXAr2DMask.GetDim1() > m_rXArray2D.GetDim1() || rXAr2DMask.GetDim2() > m_rXArray2D.GetDim2())
			throw std::invalid_argument("invalid argument rXAr2DMask in xar::XArray2DMove<T>::MaskFilter (mask too large for the image)");

		index_t iPointsLeft = iMaskCentreX;
		index_t iPointsRight = rXAr2DMask.GetDim2() - iMaskCentreX - 1;	
		index_t iPointsUp = iMaskCentreY;
		index_t iPointsDown = rXAr2DMask.GetDim1() - iMaskCentreY - 1;
		index_t iSizeX = iPointsLeft + iPointsRight + 1;
		index_t iSizeY = iPointsUp + iPointsDown + 1;

		XArray2D<T> XArTemp(m_rXArray2D);
		XArray2DMove<T> XArTempMove(XArTemp);
		XArTempMove.PadMirror(iPointsUp, iPointsDown, iPointsLeft, iPointsRight);

		index_t ny = XArTemp.GetDim1();
		index_t nx = XArTemp.GetDim2();

		index_t Nx = 2;
		while (Nx<nx) Nx *= 2; // Nx is the smallest power of 2 not less than m_rXArray2D.nx
		index_t Ny = 2;
		while (Ny<ny) Ny *= 2; // Ny is the smallest power of 2 not less than m_rXArray2D.ny
		XArTempMove.Pad((Ny-ny)/2, Ny-ny-(Ny-ny)/2, (Nx-nx)/2, Nx-nx-(Nx-nx)/2, T(XArTemp.Norm(eNormAver)));

		XArray2D<T> XArKernel(XArTemp);
		XArKernel.Fill(T(0));
		index_t cy = (XArKernel.GetDim1() - 1) / 2;
		index_t cx = (XArKernel.GetDim2() - 1) / 2;

		for (index_t i = 0; i < iSizeY; i++)
			for (index_t j = 0; j < iSizeX; j++)
				XArKernel[cy + i - iPointsUp][cx + j - iPointsLeft] = rXAr2DMask[i][j];

		XArray2DFFTRe<T> FFTRe(XArTemp);
		long lngShift1, lngShift2;
		bool bTruncateThis = false;
		bool bTruncateRefArr = true;
		FFTRe.CorrelateOld(XArKernel, lngShift1, lngShift2, bTruncateThis, bTruncateRefArr);
		XArTempMove.Trim((Ny-ny)/2+iPointsUp, Ny-ny-(Ny-ny)/2+iPointsDown, (Nx-nx)/2+iPointsLeft, Nx-nx-(Nx-nx)/2+iPointsRight);

		m_rXArray2D = XArTemp;
	}

	// This is a specialization of the function xar::XArray2DMove<T>::MaskFilterFourier for T=char
	template<> inline void XArray2DMove<char>::MaskFilterFourier(const XArray2D<char>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<char>::MaskFilterFourier (does not work for integer data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::MaskFilterFourier for T=short
	template<> inline void XArray2DMove<short>::MaskFilterFourier(const XArray2D<short>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<short>::MaskFilterFourier (does not work for integer data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::MaskFilterFourier for T=long
	template<> inline void XArray2DMove<long>::MaskFilterFourier(const XArray2D<long>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<long>::MaskFilterFourier (does not work for integer data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::MaskFilterFourier for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::MaskFilterFourier(const XArray2D<xar::fcomplex>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<fcomplex>::MaskFilterFourier (does not work for complex data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::PSFMoments for T=dcomplex
	template<> inline void XArray2DMove<xar::dcomplex>::MaskFilterFourier(const XArray2D<xar::dcomplex>& rXAr2DMask, index_t iMaskCentreY, index_t iMaskCentreX)
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<dcomplex>::MaskFilterFourier (does not work for complex data)");
	}


	//---------------------------------------------------------------------------
	//Function XArray2DSavGol<T>::Derivative
	//
	//	Calculates partial derivatives of the array
	//
	/*!
		\brief		Calculates partial derivatives of the array
		\param		iOrderY is the y-order of the partial derivative (must not exceed 4)
		\param		iOrderX is the x-order of the partial derivative (must not exceed 4)
		\exception	std::invalid_argument is thrown if input parameters are inconsistent
		\return		\a None
		\par		Description:
			This function applies masks to calculate finite-difference partial derivatives of the array. It uses
			mirror reflection at the boundaries.
	*/	
	template<class T> void XArray2DMove<T>::Derivative(index_t iOrderY, index_t iOrderX)
	{
		if (iOrderY + iOrderX == 0) 
			return;

		index_t iMaskCentreY, iMaskCentreX;
		XArray2D<T> xaMask;

		// calculate partial y-derivative
		switch (iOrderY)
		{
		case 0:
			break;
		case 1:
			xaMask.Resize(3, 1, 0); 
			xaMask[0][0] = T(-0.5); 
			xaMask[2][0] = T(0.5);
			iMaskCentreY = 1; 
			iMaskCentreX = 0;
			break;
		case 2:
			xaMask.Resize(3, 1, 0); 
			xaMask[0][0] = T(1); 
			xaMask[1][0] = T(-2); 
			xaMask[2][0] = T(1);
			iMaskCentreY = 1; 
			iMaskCentreX = 0;
			break;
		case 3:
			xaMask.Resize(5, 1, 0); 
			xaMask[0][0] = T(-0.5); 
			xaMask[1][0] = T(1); 
			xaMask[3][0] = T(0.5); 
			xaMask[4][0] = T(-1);
			iMaskCentreY = 2; 
			iMaskCentreX = 0;
			break;
		case 4:
			xaMask.Resize(5, 1, 0); 
			xaMask[0][0] = T(1); 
			xaMask[1][0] = T(-4); 
			xaMask[2][0] = T(6); 
			xaMask[3][0] = T(-4); 
			xaMask[4][0] = T(1);
			iMaskCentreY = 2; 
			iMaskCentreX = 0;
			break;
		default:
			throw std::invalid_argument("invalid arguments 'iOrderY' in XArray2DMove<T>::Derivative (derivative orders must not exceed 4)");
		}
	
		if (iOrderY > 0) MaskFilter(xaMask, iMaskCentreY, iMaskCentreX);

		// calculate partial x-derivative
		xaMask *= T(0);
		switch (iOrderX)
		{
		case 0:
			break;
		case 1:
			xaMask.Resize(1, 3, 0); 
			xaMask[0][0] = T(-0.5); 
			xaMask[0][2] = T(0.5);
			iMaskCentreY = 0; 
			iMaskCentreX = 1;
			break;
		case 2:
			xaMask.Resize(1, 3, 0); 
			xaMask[0][0] = T(1); 
			xaMask[0][1] = T(-2); 
			xaMask[0][2] = T(1);
			iMaskCentreY = 0; 
			iMaskCentreX = 1;
			break;
		case 3:
			xaMask.Resize(1, 5, 0); 
			xaMask[0][0] = T(-0.5); 
			xaMask[0][1] = T(1); 
			xaMask[0][3] = T(0.5); 
			xaMask[0][4] = T(-1);
			iMaskCentreY = 0; 
			iMaskCentreX = 2;
			break;
		case 4:
			xaMask.Resize(1, 5, 0); 
			xaMask[0][0] = T(1); 
			xaMask[0][1] = T(-4); 
			xaMask[0][2] = T(6); 
			xaMask[0][3] = T(-4); 
			xaMask[0][4] = T(1);
			iMaskCentreY = 0; 
			iMaskCentreX = 2;
			break;
		default:
			throw std::invalid_argument("invalid arguments 'iOrderX' in XArray2DMove<T>::Derivative (derivative orders must not exceed 4)");
		}

		if (iOrderX > 0)
			MaskFilter(xaMask, iMaskCentreY, iMaskCentreX);

		double yst = GetYStep(m_rXArray2D);
		double xst = GetXStep(m_rXArray2D);
		m_rXArray2D *= T(1.0 / pow(yst, double(iOrderY)) / pow(xst, double(iOrderX)));
	}

	// specialization for T = fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::Derivative(index_t iOrderY, index_t iOrderX)
	{
		if (iOrderY + iOrderX == 0) return;

		index_t iMaskCentreY, iMaskCentreX;
		XArray2D<fcomplex> xaMask;

		// calculate partial y-derivative
		switch (iOrderY)
		{
		case 0:
			break;
		case 1:
			xaMask.Resize(3, 1, 0); 
			xaMask[0][0] = fcomplex(-0.5); 
			xaMask[2][0] = fcomplex(0.5);
			iMaskCentreY = 1; 
			iMaskCentreX = 0;
			break;
		case 2:
			xaMask.Resize(3, 1, 0);
			xaMask[0][0] = fcomplex(1); 
			xaMask[1][0] = fcomplex(-2); 
			xaMask[2][0] = fcomplex(1);
			iMaskCentreY = 1; 
			iMaskCentreX = 0;
			break;
		case 3:
			xaMask.Resize(5, 1, 0); 
			xaMask[0][0] = fcomplex(-0.5); 
			xaMask[1][0] = fcomplex(1); 
			xaMask[3][0] = fcomplex(0.5); 
			xaMask[4][0] = fcomplex(-1);
			iMaskCentreY = 2; 
			iMaskCentreX = 0;
			break;
		case 4:
			xaMask.Resize(5, 1, 0); 
			xaMask[0][0] = fcomplex(1); 
			xaMask[1][0] = fcomplex(-4); 
			xaMask[2][0] = fcomplex(6); 
			xaMask[3][0] = fcomplex(-4); 
			xaMask[4][0] = fcomplex(1);
			iMaskCentreY = 2; 
			iMaskCentreX = 0;
			break;
		default:
			throw std::invalid_argument("invalid arguments 'iOrderY' in XArray2DMove<fcomplex>::Derivative (derivative orders must not exceed 4)");
		}

		if (iOrderY > 0) 
			MaskFilter(xaMask, iMaskCentreY, iMaskCentreX);

		// calculate partial x-derivative
		xaMask *= fcomplex(0);
		switch (iOrderX)
		{
		case 0:
			break;
		case 1:
			xaMask.Resize(1, 3, 0); 
			xaMask[0][0] = fcomplex(-0.5); 
			xaMask[0][2] = fcomplex(0.5);
			iMaskCentreY = 0; 
			iMaskCentreX = 1;
			break;
		case 2:
			xaMask.Resize(1, 3, 0); 
			xaMask[0][0] = fcomplex(1); 
			xaMask[0][1] = fcomplex(-2); 
			xaMask[0][2] = fcomplex(1);
			iMaskCentreY = 0; 
			iMaskCentreX = 1;
			break;
		case 3:
			xaMask.Resize(1, 5, 0); 
			xaMask[0][0] = fcomplex(-0.5); 
			xaMask[0][1] = fcomplex(1); 
			xaMask[0][3] = fcomplex(0.5); 
			xaMask[0][4] = fcomplex(-1);
			iMaskCentreY = 0; 
			iMaskCentreX = 2;
			break;
		case 4:
			xaMask.Resize(1, 5, 0); 
			xaMask[0][0] = fcomplex(1); 
			xaMask[0][1] = fcomplex(-4); 
			xaMask[0][2] = fcomplex(6); 
			xaMask[0][3] = fcomplex(-4); 
			xaMask[0][4] = fcomplex(1);
			iMaskCentreY = 0; iMaskCentreX = 2;
			break;
		default:
			throw std::invalid_argument("invalid arguments 'iOrderX' in XArray2DMove<fcomplex>::Derivative (derivative orders must not exceed 4)");
		}

		if (iOrderX > 0)
			MaskFilter(xaMask, iMaskCentreY, iMaskCentreX);

		double yst = GetYStep(m_rXArray2D);
		double xst = GetXStep(m_rXArray2D);
		float aaa = (float)pow(yst, double(iOrderY));
		float bbb = (float)pow(xst, double(iOrderX));
		m_rXArray2D *= fcomplex(aaa * bbb, 0);
	}

	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Moment
	//
	//	Calculates integral moments of the array
	//
	/*!
		\brief		Calculates integral moments of the array
		\param		iOrderY is the Y order of the moment
		\param		iOrderX is the X order of the moment
		\param		dblCentreY is the Y coordinate of the centre of gravity
		\param		dblCentreX is the X coordinate of the centre of gravity
		\return		\a The moment
		\par		Description:
			This function returns the integral moment of the order (iOrderY, iOrderX) with respect
			to the 'centre of gravity' point defined by (dblCentreY, dblCentreX) using the trapezoidal formula
	*/	
	template<class T> T XArray2DMove<T>::Moment(index_t iOrderY, index_t iOrderX, double dblCentreY, double dblCentreX) const
	{
		index_t i, j;
		double yst = GetYStep(m_rXArray2D);
		double xst = GetXStep(m_rXArray2D);
		double yloc = GetYlo(m_rXArray2D) - dblCentreY;
		double xloc = GetXlo(m_rXArray2D) - dblCentreX;
		T sum, sumx;

		index_t iDim1 = m_rXArray2D.GetDim1();
		index_t iDim2 = m_rXArray2D.GetDim2();
		vector<T> yyc(iDim1), xxc(iDim2);
		for (i = 0; i < iDim1; i++) yyc[i] = (T)pow(yloc + yst * i, double(iOrderY));
		yyc[0] *= T(0.5); yyc[iDim1 - 1] *= T(0.5);
		for (j = 0; j < iDim2; j++) xxc[j] = (T)pow(xloc + xst * j, double(iOrderX));
		xxc[0] *= T(0.5); xxc[iDim2 - 1] *= T(0.5);

		sum = 0;
		for (i = 0; i < iDim1; i++)
		{
			sumx = 0;
			for (j = 0; j < iDim2; j++)
				sumx += xxc[j] * m_rXArray2D[i][j];

			sum += yyc[i] * sumx;
		}
		return sum * T(yst * xst);
	}

	// This is a specialization of the function xar::XArray2DMove<T>::Moment for T=fcomplex
	template<> inline fcomplex XArray2DMove<xar::fcomplex>::Moment(index_t iOrderY, index_t iOrderX, double dblCentreY, double dblCentreX) const
	{
		index_t i, j;
		float yst = (float)xar::GetYStep(m_rXArray2D);
		float xst = (float)xar::GetXStep(m_rXArray2D);
		float yloc = (float)(xar::GetYlo(m_rXArray2D) - dblCentreY);
		float xloc = (float)(xar::GetXlo(m_rXArray2D) - dblCentreX);
		xar::fcomplex sum, sumx;

		index_t iDim1 = m_rXArray2D.GetDim1();
		index_t iDim2 = m_rXArray2D.GetDim2();
		vector<float> yyc(iDim1), xxc(iDim2);
		for (i = 0; i < iDim1; i++) yyc[i] = (float)pow(yloc + yst * i, float(iOrderY));
		yyc[0] *= 0.5f; yyc[iDim1 - 1] *= 0.5f;
		for (j = 0; j < iDim2; j++) xxc[j] = (float)pow(xloc + xst * j, float(iOrderX));
		xxc[0] *= 0.5f; xxc[iDim2 - 1] *= 0.5f;

		sum = 0;
		for (i = 0; i < iDim1; i++)
		{
			sumx = 0;
			for (j = 0; j < iDim2; j++)
				sumx += xxc[j] * m_rXArray2D[i][j];

			sum += yyc[i] * sumx;
		}
		return sum * yst * xst;
	}


	// Calculates centred integral moments of the PSF (used in Taylor deconvolution)
	template<class T> void XArray2DMove<T>::PSFMoments(index_t iMaxOrder, XArray2D<T>& rMoments) const
	{
		// decomposition orders < 1 and > 10 are not really practical
		if( iMaxOrder < 1 || iMaxOrder > 10 )
			throw std::invalid_argument("invalid argument 'iMaxOrder' in XArray2DMove<T>::PSFMoments (maximum order too small or too large)");

		// if m_rXArray2D does not have a head, create a default one
		IXAHWave2D* ph2( GetIXAHWave2D( m_rXArray2D ) );
		bool bHeadChanged( false );

		if( ph2 == 0 )
		{
			ph2 = CreateWavehead2D();
			ph2->SetData(1, 0, static_cast<double>( m_rXArray2D.GetDim1() - 1 ), 0, static_cast<double>( m_rXArray2D.GetDim2() - 1 ) );
			m_rXArray2D.SetHeadPtr( ph2 );
			bHeadChanged = true;
		}

		rMoments.Resize(iMaxOrder + 1, iMaxOrder + 1);

		// calculate zero-order moment
		rMoments[0][0] = Moment(0, 0, 0, 0);
	
		if (rMoments[0][0] == 0)
			throw std::invalid_argument("invalid argument 'm_rXArray2D' in XArray2DMove<T>::PSFMoments (integral of the PSF is zero)");

		// calculate first-order moments
		rMoments[1][0] = Moment(1, 0, 0, 0);
		rMoments[0][1] = Moment(0, 1, 0, 0);

		// calculate higher-order centred moments
		vector<T> fact(iMaxOrder + 1);
		fact[0] = 1;
	
		for (index_t k = 1; k <= iMaxOrder; k++) 
			fact[k] = static_cast<T>(k) * fact[k-1];

		for (index_t m = 0; m <= iMaxOrder; m++)
			for (index_t n = 0; n <= iMaxOrder - m; n++)
			{
				if (m + n < 2) 
					continue;

				rMoments[m][n] = ((m + n) % 2 == 0 ? 1 : -1) / (fact[m] * fact[n]) * Moment(m, n, rMoments[1][0], rMoments[0][1]);
			}

		// assume the PSF to be centred at the 'centre of gravity' (this makes first moments equal to zero)
		rMoments[1][0] = rMoments[0][1] = T(0);

		// remove the temporary head if added
		if (bHeadChanged)
			m_rXArray2D.SetHeadPtr(0);
	}

	// This is a specialization of the function xar::XArray2DMove<T>::PSFMoments for T=fcomplex
	template<> inline void XArray2DMove<xar::fcomplex>::PSFMoments(index_t iMaxOrder, XArray2D<xar::fcomplex>& rMoments) const
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<fcomplex>::PSFMoments (does not work for complex data)");
	}

	// This is a specialization of the function xar::XArray2DMove<T>::PSFMoments for T=dcomplex
	template<> inline void XArray2DMove<xar::dcomplex>::PSFMoments(index_t iMaxOrder, XArray2D<xar::dcomplex>& rMoments) const
	{
		throw std::invalid_argument("invalid argument '*this' in XArray2DMove<dcomplex>::PSFMoments (does not work for complex data)");
	}


	//! Calculates coefficients of the inverse Taylor series (used in deconvolution)
	template<class T> void XArray2DMove<T>::InverseSeries(const XArray2D<T>& rAcoeff, XArray2D<T>& rBcoeff)
	{
		index_t iDim = rAcoeff.GetDim1();
		if (iDim != rAcoeff.GetDim2())
			throw std::invalid_argument("invalid argument 'rAcoeff' in XArray2DMove<T>::InverseSeries (rAcoeff is not square)");
		if (rAcoeff[0][0] == T(0)) 
			throw std::invalid_argument("invalid argument 'rAcoeff' in XArray2DMove<T>::InverseSeries (rAcoeff[0][0] is zero)");
		rBcoeff.Resize(iDim, iDim, T(0));
		rBcoeff[0][0] = T(1) / rAcoeff[0][0] ;

		index_t i, j, ii, jj;
		for (i = 0; i < iDim; i++)
			for (j = 0; j < iDim - i; j++)
			{
				if (i + j == 0) continue;
				for (ii = 0; ii < i; ii++)
					for (jj = 0; jj < j; jj++)
						rBcoeff[i][j] += rBcoeff[ii][jj] * rAcoeff[i-ii][j-jj];
				for (ii = 0; ii < i; ii++)
					rBcoeff[i][j] += rBcoeff[ii][j] * rAcoeff[i-ii][0];
				for (jj = 0; jj < j; jj++)
					rBcoeff[i][j] += rBcoeff[i][jj] * rAcoeff[0][j-jj];
				rBcoeff[i][j] *= -rBcoeff[0][0];
			}
	}


	//! Divides the array by the average over the region
	template<class T> void XArray2DMove<T>::DoRegionNormalisation(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2)
	{
		XArray2D<T> xaNormRegion;
		m_rXArray2D.GetSubarray( iBeginDim1, iEndDim1, iBeginDim2, iEndDim2, xaNormRegion );

		// Calculate the average value within the region
		double dblRegionAvg( xaNormRegion.Norm( xar::eNormAver ) );

		if( dblRegionAvg == 0 )
			throw std::invalid_argument("invalid arguments in XArray2DMove<T>::DoRegionNormalisation (bad normalization region - average value in the region is 0)");

		double dblInvRegionAvg = 1 / dblRegionAvg;

		// Divide the array by the region average
		m_rXArray2D *= static_cast<T>( dblInvRegionAvg );
	}


	//! Divides each row of the array by the corresponding average value over its sub-row
	//	iBeginDim2 and iEndDim2 is the starting and ending index of the region, respectively (both are inclusive!)
	template<class T> void XArray2DMove<T>::DoRowWiseNormalisation( index_t iBeginDim2, index_t iEndDim2, int op )
	{
		index_t dim1( m_rXArray2D.GetDim1() );
		index_t dim2( m_rXArray2D.GetDim2() );
		
		if( iBeginDim2 >= iEndDim2 || iEndDim2 > dim2 - 1 )
			throw std::invalid_argument("invalid arguments 'iBeginDim2' and 'iEndDim2' in XArray2DMove<T>::DoRowWiseNormalisation (the indexes are outside array's boundaries)");

		XArray2D<T> xaNormRegion;
		
		for( index_t i = 0; i < dim1; i++ )
		{
			// NOTE: GetSubarray() treats the end indeces as exclusive (hence + 1)
			m_rXArray2D.GetSubarray( i, i+1, iBeginDim2, iEndDim2 + 1, xaNormRegion );
			
			// Calculate the average value within the region
			double dblRegionAvg( xaNormRegion.Norm( xar::eNormAver ) );

			switch( op )
			{
			case 0: // division
				{
				if( dblRegionAvg == 0 )
					// If the row average is 0, set it to 1 to avoid division by zero condition
					dblRegionAvg = 1;
					//throw std::invalid_argument("invalid arguments in XArray2DMove<T>::DoRegionNormalisation (bad normalization row(s) - average value in the row(s) is 0)");

				double dblInvRegionAvg( 1 / dblRegionAvg );
			
				// Divide the row by the sub-row average
				for( index_t j = 0; j < dim2; j++ )
					m_rXArray2D[i][j] *= static_cast<T>( dblInvRegionAvg );
				}
				break;
			case 1: // subtraction
				{
				if( dblRegionAvg != 0 ) 
					// Subtract the sub-row average from the row
					for( index_t j = 0; j < dim2; j++ )
						m_rXArray2D[i][j] -= static_cast<T>( dblRegionAvg );
				}
				break;
			default:
				throw std::invalid_argument( "invalid argument 'op' in XArray2DMove<T>::DoRowWiseNormalisation() (non-supported operation)" );
			}
		}			
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Cart2Polar
	//
	//	In-place conversion of the array from Cartesian to polar coordinates
	//
	/*!
		\brief		Performs in-place conversion of the array from Cartesian to polar coordinates
		\param		bOdd (out) is an auxiliary parameter that specifies whether the size of the original array if odd (bOdd = true) or even (bOdd = false)
		\par		Description:
			This function performes in-place conversion of the array from the Cartesian to polar coordinates			
	*/	
	template<class T> void XArray2DMove<T>::Cart2Polar( bool& bOdd )
	{
		index_t nx( m_rXArray2D.GetDim2() );

		if( m_rXArray2D.GetDim1() != nx )
			throw std::invalid_argument( "Invalid argument 'm_rXArray2D' in XArray2DMove<T>::Cart2Polar() (MUST be square)" );

		IXAHWave2D* ph2( GetIXAHWave2D( m_rXArray2D ) );

		if( ph2 != 0 && ph2->GetXStep( nx ) != ph2->GetYStep( nx ) )
			throw std::invalid_argument( "Invalid argument 'm_rXArray2D' in XArray2DMove<T>::Cart2Polar() (MUST be square)" );

		bOdd = !(( nx % 2 ) == 0);

		index_t nr( bOdd ? ( nx - 1 ) / 2 + 1 : nx / 2 );
		index_t nphi( static_cast<index_t>( xar::tPI * nr + 0.5 ) );
		double dphi( xar::tPI / nphi );

		xar::XArray2D<T> XArPolar( nphi, nr );
		double num2( 0.5 * ( nx - 1 ) );

		for( index_t i = 0; i < nphi; i++ )
		{
			double phi( i * dphi - xar::PI2 ); //@@@

			for( index_t j = 0; j < nr; j++ )
			{
				double r( static_cast<double>( j ) );
				double x( r * cos( phi ) + num2 );
				double y( r * sin( phi ) + num2 );

				index_t ix( static_cast<index_t>( x ) );
				index_t iy( static_cast<index_t>( y ) );
				index_t ix1( ix + 1 );
				index_t iy1( iy + 1 );

				if( ix1 < nx && iy1 < nx )
				{
					double dx( x - static_cast<double>( ix ) );
					double dy( y - static_cast<double>( iy ) );
					double a00( ( 1. - dy ) * ( 1. - dx ) );
					double a01( ( 1. - dy ) * dx );
					double a10( dy * ( 1. - dx ) );
					double a11( dy * dx );
					XArPolar[i][j] = static_cast<T>( m_rXArray2D[iy][ix] * a00 + m_rXArray2D[iy][ix1] * a01 + m_rXArray2D[iy1][ix] * a10 + m_rXArray2D[iy1][ix1] * a11);
				}
			}
		}

		// if m_rXArray2D has a header, create a header
		if( ph2 != 0 )
		{
			IXAHWave2D* phNew( CreateWavehead2D() );
			phNew->SetData( ph2->GetWl(), 0, dphi * ( nphi - 1 ), 0, ph2->GetXStep( nx ) * ( nr - 1 ) );
			XArPolar.SetHeadPtr( phNew );
		}
		
		m_rXArray2D.Swap( XArPolar );
	}


	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::Polar2Cart
	//
	//	In-place conversion of the array from polar to Cartesian coordinates
	//
	/*!
		\brief		Performs in-place conversion of the array from polar to Cartesian coordinates
		\param		bOdd (in) is an auxiliary parameter that specifies whether the size of the array (in the Cartesian coordinates) if odd (bOdd = true) or even (bOdd = false)
		\par		Description:
			This function performes in-place conversion of the array from the polar to the Cartesian coordinates			
	*/	
	template<class T> void XArray2DMove<T>::Polar2Cart( bool bOdd )
	{
		index_t nphi( m_rXArray2D.GetDim1() );
		index_t nr( m_rXArray2D.GetDim2() );

		index_t nx( bOdd ? 2 * nr - 1 : 2 * nr );
		double invdphi( nphi / xar::tPI );

		xar::XArray2D<T> XArCart( nx, nx ); // a square array
		double num2( 0.5 * ( nx - 1 ) );

		for( index_t i = 0; i < nx; i++ )
		{
			double y( static_cast<double>( i ) - num2 );

			for( index_t j = 0; j < nx; j++ )
			{
				double x( static_cast<double>( j ) - num2 );
				double ro( sqrt( x * x + y * y ) );
				index_t iro( static_cast<index_t>( ro ) );//index_t(bOdd ? ro : ro - 0.5);
				index_t iro1( iro + 1 );

				if( iro1 < nr )
				{
					double phi;

					if( ro == 0 )
						phi = 0;
					else
					{
						phi = acos( x / ro );

						if( y < 0 )
							phi = xar::tPI - phi;

						phi += xar::PI2; //@@@

						if( phi >= xar::tPI )
							phi -= xar::tPI;

						//phi -= int(phi / xar::tPI) * xar::tPI;
						//if (phi < 0) phi += xar::tPI;
						phi *= invdphi;
					}

					double d_ro( ro - static_cast<double>( iro ) );
					index_t iphi( static_cast<index_t>( phi ) );
					index_t iphi1( iphi + 1 );

					if( iphi1 >= nphi )
						iphi1 -= nphi;

					double d_phi( phi - static_cast<double>( iphi ) );
					double a00( ( 1. - d_phi ) * ( 1. - d_ro ) );
					double a01( ( 1. - d_phi ) * d_ro );
					double a10( d_phi * ( 1. - d_ro ) );
					double a11( d_phi * d_ro );

					XArCart[i][j] = static_cast<T>( m_rXArray2D[iphi][iro] * a00 + m_rXArray2D[iphi][iro1] * a01 + m_rXArray2D[iphi1][iro] * a10 + m_rXArray2D[iphi1][iro1] * a11 );
				}
			}
		}

		IXAHWave2D* ph2( GetIXAHWave2D( m_rXArray2D ) );

		// if m_rXArray2D has a header, create a header
		if( ph2 != 0 )
		{
			IXAHWave2D* phNew( CreateWavehead2D() );
			phNew->SetData( ph2->GetWl(), 0, ph2->GetXStep( nr ) * ( nx - 1 ), 0, ph2->GetXStep( nr ) * (nx - 1) );
			XArCart.SetHeadPtr( phNew );
		}
		
		m_rXArray2D.Swap( XArCart );
	}

	//---------------------------------------------------------------------------
	//Function XArray2DMove<T>::CorrelationCoeff
	//
	//	Correlation coefficient of the array with a reference array
	//
	/*!
		\brief		Calculates the correlation coefficient between the array with a reference array
		\param		rXArRef is the reference array
		\par		Description:
			This function calculates the correlation coefficient between the array with a reference array
	*/	
	template<class T> double XArray2DMove<T>::CorrelationCoeff( const XArray2D<T>& rXArRef )
	{
		index_t ny( m_rXArray2D.GetDim1() );
		index_t nx( m_rXArray2D.GetDim2() );

		if (rXArRef.GetDim1() != ny || rXArRef.GetDim2() != nx)
			throw std::invalid_argument("Invalid arguments 'm_rXArray2D' and 'rXArRef' in XArray2DMove<T>::CorrelationCoeff() (different dimensions)");

		double averXAr = m_rXArray2D.Norm(xar::eNormAver);
		double averXArRef = rXArRef.Norm(xar::eNormAver);
		double stdXAr = m_rXArray2D.Norm(xar::eNormStdDev);
		double stdXArRef = rXArRef.Norm(xar::eNormStdDev);

		if (stdXAr == 0 || stdXArRef == 0)
			throw std::invalid_argument("Invalid argument 'm_rXArray2D' or/and 'rXArRef' in XArray2DMove<T>::CorrelationCoeff() (zero standard deviation)");

		xar::XArray2D<T> XArTemp1(m_rXArray2D);
		XArTemp1 -= averXAr;
		xar::XArray2D<T> XArTemp2(rXArRef);
		XArTemp2 -= averXArRef;
		IXAHWave2D* ph2( GetIXAHWave2D( XArTemp1 ) );
		if (ph2) XArTemp2.SetHeadPtr( ph2->Clone() );
		else
		{
			ph2 = GetIXAHWave2D( XArTemp2 );
			if (ph2) XArTemp1.SetHeadPtr( ph2->Clone() );
		}

		XArTemp1 *= XArTemp2;
		double cov = XArTemp1.Norm(xar::eNormAver);

		return cov / (stdXAr * stdXArRef);
	}

	//!	Rebins a 2D array into a new one (the original array is destroyed during the rebinning)
	//	dblXlo is the left boundary of the new array (microns)
	//	dblXstep is the horizontal step of the new array (microns)
	//	iNumX is the horizontal size of the new array (pixels)
	//	dblYlo is the bottom boundary of the new array (microns)
	//	dblYstep is the vertical step of the new array (microns)
	//	iNumY is the vertical size of the new array (pixels)
	template<class T> void xar::XArray2DMove<T>::ReBin(double dblXlo, double dblXstep, index_t iNumX, double dblYlo, double dblYstep, index_t iNumY, bool bNormalize)
	{
		IXAHWave2D* ph2( GetIXAHWave2D( m_rXArray2D ) );
		
		if( !ph2 )
			throw std::invalid_argument("invalid argument 'm_rXArray2D' in XArray2DMove<T>::ReBin (the array should have a header)");

		index_t numXold( m_rXArray2D.GetDim2() );
		index_t numYold( m_rXArray2D.GetDim1() );

		double xlo( ph2->GetXlo() );
		double ylo( ph2->GetYlo() );
		double Wl( ph2->GetWl() );

		double dx( GetXStep( m_rXArray2D ) );
		double dy( GetYStep( m_rXArray2D ) );

		std::vector<double> x, y;
		std::vector<T> s, d;

		// Rebining over the X-coordinate
		xar::XArray2D<T> temp( numYold, iNumX );
		x.resize( numXold + 1 );
		
		for (index_t i = 0; i <= numXold; i++) 
			x[i] = xlo + i * dx;

		y.resize(iNumX + 1);
		
		for( index_t i = 0; i <= iNumX; i++ ) 
			y[i] = dblXlo + i * dblXstep;

		s.resize(numXold);
		d.resize(iNumX);

		for (index_t i = 0; i < numYold; i++)
		{
			for (index_t j = 0; j < numXold; j++) s[j] = m_rXArray2D[i][j];
			ReBin1D(x, s, y, d, true);
			for (index_t j = 0; j < iNumX; j++) temp[i][j] = d[j];
		}
		//////////////////////////////////

		// Rebinning over the Y-coordinate
		m_rXArray2D.Resize( iNumY, iNumX );

		IXAHWave2D* phNew( CreateWavehead2D() );
		phNew->SetData( Wl, dblYlo, dblYlo + ( iNumY - 1 ) * dblYstep, dblXlo, dblXlo + ( iNumX - 1 ) * dblXstep );
		m_rXArray2D.SetHeadPtr( phNew );

		x.resize( numYold + 1) ;
		for( index_t i = 0; i <= numYold; i++ ) 
			x[i] = ylo + i * dy;
		
		y.resize( iNumY + 1 );

		for( index_t i = 0; i <= iNumY; i++ ) 
			y[i] = dblYlo + i * dblYstep;

		s.resize(numYold);
		d.resize(iNumY);

		for (index_t i = 0; i < iNumX; i++)
		{
			for (index_t j = 0; j < numYold; j++) s[j] = temp[j][i];
			ReBin1D(x, s, y, d, true);
			for (index_t j = 0; j < iNumY; j++) m_rXArray2D[j][i] = d[j];
		}
		//////////////////////////////////

		if (bNormalize)
			m_rXArray2D *= (1. / (dblXstep * dblYstep));
		else
			m_rXArray2D *= (1. / (dx * dy));
	}

	//!	Applies modulus-gradient transform
	template<class T> void xar::XArray2DMove<T>::ModGrad(void)
	{
		double invxst2 = std::pow(1 / xar::GetXStep(m_rXArray2D), 2);
		double invyst2 = std::pow(1 / xar::GetYStep(m_rXArray2D), 2);

		long ny = m_rXArray2D.GetDim1();
		long nx = m_rXArray2D.GetDim2();

		xar::XArray2D<T> xaModGrad(ny, nx);
		// set the same head
		xaModGrad.SetHeadPtr(m_rXArray2D.GetHeadPtr() ? m_rXArray2D.GetHeadPtr()->Clone() : 0);

		long zero(0);
		for( long i = 0; i < ny; i++ )
		{
			double cy( invyst2 * (i == 0 || i == ny-1 ? 1. : 0.25) );

			for( long j = 0; j < nx; j++ )
			{
				double cx( invxst2 * (j == 0 || j == nx-1 ? 1. : 0.25 ) );

				T xleft( m_rXArray2D[i][ std::max( zero, j - 1 ) ] );
				T xright( m_rXArray2D[i][ std::min( nx - 1, j + 1 ) ] );
				T yleft( m_rXArray2D[ std::max( zero, i - 1 ) ][j] );
				T yright( m_rXArray2D[ std::min( ny - 1, i + 1 ) ][j] );

				T dx = (xright > xleft ? xright - xleft : xleft - xright);
				T dy = (yright > yleft ? yright - yleft : yleft - yright);

				xaModGrad[i][j] = (T)sqrt( dx * (dx * cx) + dy * (dy * cy) );
			}
		}

		m_rXArray2D.Swap( xaModGrad );
	}

	//!	Applies modulus-gradient transform using right differences
	template<class T> void xar::XArray2DMove<T>::ModGradRight(void)
	{
		double invxst2 = std::pow(1 / xar::GetXStep(m_rXArray2D), 2);
		double invyst2 = std::pow(1 / xar::GetYStep(m_rXArray2D), 2);

		size_t ny = m_rXArray2D.GetDim1();
		size_t nx = m_rXArray2D.GetDim2();

		xar::XArray2D<T> xaModGrad(ny, nx, T(0));
		// set the same head
		xaModGrad.SetHeadPtr(m_rXArray2D.GetHeadPtr() ? m_rXArray2D.GetHeadPtr()->Clone() : 0);

		for( size_t i = 0; i < ny - 1; i++ )
		{
			for( size_t j = 0; j < nx - 1; j++ )
			{
				T cntr( m_rXArray2D[i][j] );
				T xright( m_rXArray2D[i][j + 1] );
				T yright( m_rXArray2D[i + 1][j] );

				T dx = (xright > cntr ? xright - cntr : cntr - xright);
				T dy = (yright > cntr ? yright - cntr : cntr - yright);

				xaModGrad[i][j] = (T)sqrt( dx * (dx * invxst2) + dy * (dy * invyst2) );
			}
		}

		m_rXArray2D.Swap( xaModGrad );
	}

	//!	Performs a 4-connected erosion
	template<class T> void xar::XArray2DMove<T>::Erode4(T tHigh)
	{
		xar::XArray2D<T> xaCopy(m_rXArray2D);
		size_t nx = xaCopy.GetDim2();
		size_t ny = xaCopy.GetDim1();

		// Mirror pad the image
		xar::XArray2DMove<T> Move(xaCopy);
		//Move.PadMirror(1, 1, 1, 1);
		Move.Pad(1, 1, 1, 1, tHigh);

		// Main loop
		for (size_t i = 1; i <= ny; i++) {
			for (size_t j = 1; j <= nx; j++) {
				if (xaCopy[i][j] != 0)
				{
					bool bFlag = (xaCopy[i-1][j] == 0 || xaCopy[i+1][j] == 0 || xaCopy[i][j-1] == 0 || xaCopy[i][j+1] == 0);
					m_rXArray2D[i-1][j-1] = (bFlag ? 0 : tHigh);
				}
			}
		}
	}


	//!	Performs an 8-connected erosion
	template<class T> void xar::XArray2DMove<T>::Erode8(T tHigh)
	{
		xar::XArray2D<T> xaCopy(m_rXArray2D);
		size_t nx = xaCopy.GetDim2();
		size_t ny = xaCopy.GetDim1();

		// Mirror pad the image
		xar::XArray2DMove<T> Move(xaCopy);
		Move.PadMirror(1, 1, 1, 1);

		// Main loop
		for (size_t i = 1; i <= ny; i++) {
			for (size_t j = 1; j <= nx; j++) {
				if (xaCopy[i][j] != 0)
				{
					bool bFlag = (xaCopy[i-1][j-1] == 0 || xaCopy[i-1][j] == 0 || xaCopy[i-1][j+1] == 0
						|| xaCopy[i][j-1] == 0 || xaCopy[i][j+1] == 0
						|| xaCopy[i+1][j-1] == 0 || xaCopy[i+1][j] == 0 || xaCopy[i+1][j+1] == 0);
					m_rXArray2D[i-1][j-1] = (bFlag ? 0 : tHigh);
				}
			}
		}
	}


	//!	Performs a 4-connected dilation
	template<class T> void xar::XArray2DMove<T>::Dilate4(T tHigh)
	{
		xar::XArray2D<T> xaCopy(m_rXArray2D);
		size_t nx = xaCopy.GetDim2();
		size_t ny = xaCopy.GetDim1();

		// Mirror pad the image
		xar::XArray2DMove<T> Move(xaCopy);
		//Move.PadMirror(1, 1, 1, 1);
		Move.Pad(1, 1, 1, 1, 0);

		// Main loop
		for (size_t i = 1; i <= ny; i++) {
			for (size_t j = 1; j <= nx; j++) {
				if (xaCopy[i][j] == 0)
				{
					bool bFlag = (xaCopy[i-1][j] != 0 || xaCopy[i+1][j] != 0 || xaCopy[i][j-1] != 0 || xaCopy[i][j+1] != 0);
					m_rXArray2D[i-1][j-1] = (bFlag ? tHigh : 0);
				}
				else
					m_rXArray2D[i-1][j-1] = tHigh;
			}
		}
	}


	//!	Performs an 8-connected dilation
	template<class T> void xar::XArray2DMove<T>::Dilate8(T tHigh)
	{
		xar::XArray2D<T> xaCopy(m_rXArray2D);
		size_t nx = xaCopy.GetDim2();
		size_t ny = xaCopy.GetDim1();

		// Mirror pad the image
		xar::XArray2DMove<T> Move(xaCopy);
		Move.PadMirror(1, 1, 1, 1);

		// Main loop
		for (size_t i = 1; i <= ny; i++) {
			for (size_t j = 1; j <= nx; j++) {
				if (xaCopy[i][j] == 0)
				{
					bool bFlag = (xaCopy[i-1][j-1] != 0 || xaCopy[i-1][j] != 0 || xaCopy[i-1][j+1] != 0
						|| xaCopy[i][j-1] != 0 || xaCopy[i][j+1] != 0
						|| xaCopy[i+1][j-1] != 0 || xaCopy[i+1][j] != 0 || xaCopy[i+1][j+1] != 0);
					m_rXArray2D[i-1][j-1] = (bFlag ? tHigh : 0);
				}
				else
					m_rXArray2D[i-1][j-1] = tHigh;
			}
		}
	}


	//! Rebins a 1D array into a new one using a distance-driven alrogithm
	//	s is the original data array
	//	x is the array containing the boundaries of the original data array
	//	d is the destination data array
	//	y is the array containing the boundaries of the destination array
	//	Initialize defines whether to initialize the destination array (true by default) or not
	template<class T> void xar::XArray2DMove<T>::ReBin1D( std::vector<double>& x, std::vector<T>& s, std::vector<double>& y, std::vector<T>& d, bool Initialize )
	{
		xar::index_t ns( s.size() );
		
		if( x.size() != ns + 1 )
			throw std::invalid_argument("invalid_argument 'x' in XArray2DMove<T>::ReBin1D() (wrong size)");
		index_t nd( d.size() );

		if( y.size() != nd + 1 )
			throw std::invalid_argument("invalid_argument 'y' in XArray2DMove<T>::ReBin1D() (wrong size)");

		if( Initialize )
			for( xar::index_t i = 0; i < nd; i++ ) 
				d[i] = T(0);

		index_t scount( 0 );
		index_t dcount( 0 );

		double* px = &(x[0]);
		double xm = *(px);
		px++;
		double xm1 = *(px);
		double* py = &(y[0]);
		double yn = *(py);
		py++;
		double yn1 = *(py);
		T* ps = &(s[0]);
		T* pd = &(d[0]);

		// Positioning to the first overlapping intervals [xm, xm1] and [yn, yn1]
		if (xm1 <= yn)
		{
			while (xm1 <= yn && scount < ns)
			{		
				xm = xm1;
				px++;
				xm1 = *(px);
				ps++;
				scount++;
			}
		}
		else if (yn1 <= xm)
		{
			while (yn1 <= xm && dcount < nd)
			{
				yn = yn1;
				py++;
				yn1 = *(py);
				pd++;
				dcount++;
			}
		}

		while( scount < ns && dcount < nd )
		{
			// Calculating the length of overlapping
			double l( std::min( xm1, yn1 ) - std::max( xm, yn ) );
			////

			// Contribution of the current source to the current detector
			*(pd) += l * *(ps);
			////

			if (xm1 < yn1)
			{
				xm = xm1;
				px++;
				xm1 = *(px);
				ps++;
				scount++;
			}
			else if (xm1 > yn1)
			{
				//			*(pd) /= (yn1 - yn);
				yn = yn1;
				py++;
				yn1 = *(py);
				pd++;
				dcount++;
			}
			else
			{
				//			*(pd) /= (yn1 - yn);
				xm = xm1;
				px++;
				xm1 = *(px);
				ps++;
				scount++;
				yn = yn1;
				py++;
				yn1 = *(py);
				pd++;
				dcount++;
			}
		}
		//	if (scount == ns) *(pd) /= (yn1 - yn);
	}
}  //namespace xar closed
//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings are generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.


// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray2DMove<char>;
	template class xar::XArray2DMove<short>;
	template class xar::XArray2DMove<long>;
	template class xar::XArray2DMove<float>;
	template class xar::XArray2DMove<double>;
	template class xar::XArray2DMove<xar::fcomplex>;
	template class xar::XArray2DMove<xar::dcomplex>;
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
