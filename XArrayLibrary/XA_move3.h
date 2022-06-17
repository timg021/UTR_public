//Header XA_move2.h
//
//
//	HEADER FILE TITLE:
//
//		Three dimensional geometrical transformations
//
/*!
	\file		XA_move3.h
	\brief		Three dimensional geometrical transformations
	\par		Description:
		This class implements simple 'geometrical' transformations of
		XArray3D<T> objects
*/
#pragma once
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_head3.h"

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
//Class XArray3DMove
//
//	Three dimensional geometrical transformations
//
/*!
	\brief		Three dimensional geometrical transformations
	\par		Description:
				This class template defines a 'wrapper' around the XArray3D<T> object
				on which it operates; it contains functions implementing simple
				geometrical transformations of the 'wrapped' XArray3D<T> object
	\remarks    If a Wavehead3D is attached to the XArray3D<T> object, it is transformed too
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
*/
	template <class T> class XArray3DMove
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray3DMove(XArray3D<T>& rXAr3D) : m_rXArray3D(rXAr3D) {}
	protected:
		//! Copy constructor (declared protected to prohibit explicit copying)
		XArray3DMove(const XArray3DMove<T>& rCopy) : m_rXArray3D(rCopy.m_rXArray3D) {}
	public:
		//! Destructor 
		~XArray3DMove() {}

	// Operators
	protected:
		//! Assignment (declared protected to prohibit explicit copying)
		void operator=(const XArray3DMove<T>& rCopy);

	// Attributes
	public:
		//! Returns a reference to the non-modifiable 'wrapped' XArray3D<T> object
		const XArray3D<T>& GetBaseObject() const { return m_rXArray3D; }
		//! Returns a reference to the 'wrapped' XArray3D<T> object
		XArray3D<T>& GetBaseObject() { return m_rXArray3D; }

	// Operations
	public:
		//! Trims elements from the edges of the array
		void Trim(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight);
		//! Adds new elements with a given value at the edges of the array
		void Pad(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tPadVal);
		//! Replaces values at the defined perifery of the array by the given value
		void Mask(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tMaskVal);
		//! Adds new elements with a given value at the edges of the array up to the nearest integer powers of 2
		void Pad2N(T tPadVal);
		//! Replaces elements at the edges of the array using a given value 
		// void Mask(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tMaskVal);
		//! Moves elements of the array and fills the vacated positions with a given value
		// void Move(long lngMoveZPoints, long lngMoveYPoints, long lngMoveXPoints, T tFillVal);
		//! Fills a rectangular subarray with a given value
		void FillRect(index_t LowDim1, index_t HighDim1, index_t LowDim2, index_t HighDim2, index_t LowDim3, index_t HighDim3, T tFillVal);
		//! Fills a vicinity of a given point with a given value, continuing periodically over the boundary where relevant
		void FillRectPeriodic(index_t kmax, index_t jmax, index_t imax, index_t karad, index_t jarad, index_t iarad, T tFillVal);
		//! Fills a vicinity of a given point with a given value, continuing periodically over the boundary where relevant
		void FillCylinderPeriodic(index_t k0, index_t j0, index_t i0, index_t karad, index_t jarad, index_t iarad, T tFillVal);
		// Fills all points OUTSIDE the defined 3D vicinity of a point with tFillVal; if the point is near the edge, the non-filled area may extend periodically across the boundaries
		void FillRectComplementPeriodic(index_t kmax, index_t jmax, index_t imax, index_t karad, index_t jarad, index_t iarad, T tFillVal);

		// Overridables

	// Implementation
	private:
	// Member variables	
		//! Reference to the XArray3D<T> object that is being operated upon
		XArray3D<T>& m_rXArray3D;
	};

}  //namespace xar closed

//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//

//! Assignment  (declared protected to prohibit explicit copying)
template <class T> void xar::XArray3DMove<T>::operator=(const XArray3DMove<T>& rCopy)
{ 
	if (this == &rCopy) 
		return; 
	else 
		m_rXArray3D = rCopy.m_rXArray3D;
}

//---------------------------------------------------------------------------
//Function XArray3DMove<T>::Trim
//
//	Trims elements from the edges of the array
//
/*!
	\brief		Trims elements from the edges of the array
	\param		iZLeft	number of Z elements to be trimmed from the beginning
	\param		iZRight	number of Z elements to be trimmed from the end
	\param		iYLeft	number of Y elements to be trimmed from the beginning
	\param		iYRight	number of Y elements to be trimmed from the end
	\param		iXLeft	number of X elements to be trimmed from the beginning
	\param		iXRight	number of X elements to be trimmed from the end
	\exception	std::invalid_argument is thrown if the number of elements to trim is either
				negative or exceeds the array size
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
				This function 'trims' the defined number of elements from the beginning
				and end of Z, Y and X dimensions of the XArray3D<T> object; 
				it also calls the Trim function	on the IXAHWave3D head pointer if present
*/	
template <class T> void xar::XArray3DMove<T>::Trim(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight)
{
	if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0 && iZLeft == 0 && iZRight == 0) 
		return;

	if ((iXLeft + iXRight) > m_rXArray3D.GetDim3()) 
		throw std::invalid_argument("invalid_argument 'iXLeft or iXRight' in XArray3DMove<T>::Trim"); 

	if ((iYLeft + iYRight) > m_rXArray3D.GetDim2()) 
		throw std::invalid_argument("invalid_argument 'iYLeft or iYRight' in XArray3DMove<T>::Trim"); 

	if ((iZLeft + iZRight) > m_rXArray3D.GetDim1()) 
		throw std::invalid_argument("invalid_argument 'iZLeft or iZRight' in XArray3DMove<T>::Trim"); 	

	index_t dim3 = m_rXArray3D.GetDim3() - (iXLeft + iXRight);
	index_t dim2 = m_rXArray3D.GetDim2() - (iYLeft + iYRight);
	index_t dim1 = m_rXArray3D.GetDim1() - (iZLeft + iZRight);
	XArray3D<T> xarTemp(dim1, dim2, dim3);

	for (index_t i = 0; i < dim1; i++)
		for (index_t j = 0; j < dim2; j++)
			for (index_t k = 0; k < dim3; k++)
				xarTemp[i][j][k] = m_rXArray3D[i + iZLeft][j + iYLeft][k + iXLeft];

	// set the same head
	xarTemp.SetHeadPtr(m_rXArray3D.GetHeadPtr() ? m_rXArray3D.GetHeadPtr()->Clone() : 0);

	// transform the head if it can be done
	if (GetIXAHWave3D(xarTemp)) 
		GetIXAHWave3D(xarTemp)->Trim(m_rXArray3D.GetDim1(), m_rXArray3D.GetDim2(), m_rXArray3D.GetDim3(), iZLeft, iZRight, iYLeft, iYRight, iXLeft, iXRight);

	m_rXArray3D.Swap(xarTemp);
}

//Function XArray3DMove<T>::Pad
//
//	Adds new elements with a given value at the edges of the array
//
/*!
	\brief		Adds new elements with a given value at the edges of the array
	\param		iZLeft	number of Z elements to be added at the beginning
	\param		iZRight	number of Z elements to be added at the end
	\param		iYLeft	number of Y elements to be added at the beginning
	\param		iYRight	number of Y elements to be added at the end
	\param		iXLeft	number of X elements to be added at the beginning
	\param		iXRight	number of X elements to be added at the end
	\param		tPadVal	the value assigned to all new elements
	\exception	std::invalid_argument is thrown if the number of elements to pad is negative
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
				This function 'pads' the XArray3D<T> object by adding the defined number of 
				elements at the beginning and end of Z, Y and X dimensions of the array;
				it also calls the Pad function on the IXAHWave3D head pointer if present
*/	
template <class T> void xar::XArray3DMove<T>::Pad(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tPadVal)
{
	if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0 && iZLeft == 0 && iZRight == 0) 
		return;

	index_t dim3 = m_rXArray3D.GetDim3() + (iXLeft + iXRight);
	index_t dim2 = m_rXArray3D.GetDim2() + (iYLeft + iYRight);
	index_t dim1 = m_rXArray3D.GetDim1() + (iZLeft + iZRight);
	XArray3D<T> xarTemp(dim1, dim2, dim3, tPadVal);

	for (index_t i = 0; i < m_rXArray3D.GetDim1(); i++)
		for (index_t j = 0; j < m_rXArray3D.GetDim2(); j++)
			for (index_t k = 0; k < m_rXArray3D.GetDim3(); k++)
				xarTemp[i + iZLeft][j + iYLeft][k + iXLeft] = m_rXArray3D[i][j][k];

	// set the same head
	xarTemp.SetHeadPtr(m_rXArray3D.GetHeadPtr() ? m_rXArray3D.GetHeadPtr()->Clone() : 0);

	// transform the head if it can be done
	if (GetIXAHWave3D(xarTemp)) 
		GetIXAHWave3D(xarTemp)->Pad(m_rXArray3D.GetDim1(), m_rXArray3D.GetDim2(), m_rXArray3D.GetDim3(), iZLeft, iZRight, iYLeft, iYRight, iXLeft, iXRight);

	m_rXArray3D.Swap(xarTemp);
}

//Function XArray3DMove<T>::Pad2N
//
//	Adds new elements with a given value at the edges of the array up to the nearest integer powers of 2
//
/*!
	\brief		Adds new elements with a given value at the edges of the array up to the nearest integer powers of 2
	\param		tPadVal	the value assigned to all new elements
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
				This function 'pads' the XArray3D<T> object by new elements with the given value
				at the beginning and end of Z, Y and X dimensions of the array increasing its dimensions
				to the nearest powers of 2; it also calls the Pad function with the appropriate
				parameters on the IXAHWave3D head pointer if present
*/	
template <class T> void xar::XArray3DMove<T>::Pad2N(T tPadVal)
{
	index_t i = 2;
	while (i < m_rXArray3D.GetDim1()) i *= 2;
	index_t nzpad = i - m_rXArray3D.GetDim1();
	index_t j = 2;
	while (j < m_rXArray3D.GetDim2()) j *= 2;
	index_t nypad = j - m_rXArray3D.GetDim2();
	index_t k = 2;
	while (k < m_rXArray3D.GetDim3()) k *= 2;
	index_t nxpad = k - m_rXArray3D.GetDim3();

	Pad(nzpad / 2, nzpad - nzpad / 2, nypad / 2, nypad - nypad / 2, nxpad / 2, nxpad - nxpad / 2, tPadVal);
}


//---------------------------------------------------------------------------
//Function XArray3DMove<T>::Mask
//
//	Replaces elements at the edges of the array using a given value 
//
/*!
	\brief		Replaces elements at the edges of the array using a given value 
	\param		iZLeft	number of Z elements to be replaced at the beginning
	\param		iZRight	number of Z elements to be replaced at the end
	\param		iYLeft	number of Y elements to be replaced at the beginning
	\param		iYRight	number of Y elements to be replaced at the end
	\param		iXLeft	number of X elements to be replaced at the beginning
	\param		iXRight	number of X elements to be replaced at the end
	\param		tMaskVal	the value to be assigned to all masked elements
	\exception	std::invalid_argument is thrown if the number of elements to mask is negative
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
				This function 'masks' the defined number of elements at the beginning and 
				end of Z, Y and X dimensions of the XArray3D<T> object by assigning
				them a given value; the head is not affected
*/	
template <class T> void xar::XArray3DMove<T>::Mask(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight, T tMaskVal)
{
	if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0 && iZLeft == 0 && iZRight == 0) 
		return;

	if (iXLeft + iXRight > m_rXArray3D.GetDim3() || iYLeft + iYRight > m_rXArray3D.GetDim2() || iZLeft + iZRight > m_rXArray3D.GetDim1())
	{
		m_rXArray3D.Fill(tMaskVal);
	}
	else
	{
		index_t i, j, k;
		for (k = 0; k < iZLeft; k++)
			for (j = 0; j < m_rXArray3D.GetDim2(); j++)
				for (i = 0; i < m_rXArray3D.GetDim3(); i++)
					m_rXArray3D[k][j][i] = tMaskVal;
		for (k = m_rXArray3D.GetDim1() - iZRight; k < m_rXArray3D.GetDim1(); k++)
			for (j = 0; j < m_rXArray3D.GetDim2(); j++)
				for (i = 0; i < m_rXArray3D.GetDim3(); i++)
					m_rXArray3D[k][j][i] = tMaskVal;
		for (k = 0; k < m_rXArray3D.GetDim1(); k++)
			for (j = 0; j < iYLeft; j++)
				for (i = 0; i < m_rXArray3D.GetDim3(); i++)
					m_rXArray3D[k][j][i] = tMaskVal;
		for (k = 0; k < m_rXArray3D.GetDim1(); k++)
			for (j = m_rXArray3D.GetDim2() - iYRight; j < m_rXArray3D.GetDim2(); j++)
				for (i = 0; i < m_rXArray3D.GetDim3(); i++)
					m_rXArray3D[k][j][i] = tMaskVal;
		for (k = 0; k < m_rXArray3D.GetDim1(); k++)
			for (j = 0; j < m_rXArray3D.GetDim2(); j++)
				for (i = 0; i < iXLeft; i++)
					m_rXArray3D[k][j][i] = tMaskVal;
		for (k = 0; k < m_rXArray3D.GetDim1(); k++)
			for (j = 0; j < m_rXArray3D.GetDim2(); j++)
				for (i = m_rXArray3D.GetDim3() - iXRight; i < m_rXArray3D.GetDim3(); i++)
					m_rXArray3D[k][j][i] = tMaskVal;
	}
}


#if(0)
//---------------------------------------------------------------------------
//Function XArray3DMove<T>::Move
//
//	Moves elements of the array and fills the vacated positions with a given value
//
/*!
	\brief		Moves elements of the array and fills the vacated positions with a given value
	\param		lngMoveZPoints	the span of the Z translation; negative values correspond to translations in the negative Z direction
	\param		lngMoveYPoints	the span of the Y translation; negative values correspond to translations in the negative Y direction
	\param		lngMoveXPoints	the span of the X translation; negative values correspond to translations in the negative X direction
	\param		tFillVal	the value to be assigned to all vacated elements
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
				This function 'moves' the the XArray3D<T> object by translating all the
				array elements by the defined number of positions towards the beginning
				or end of Z, Y and X dimension;	it does not call any functions on the head.
				All vacated positions are filled with the defined value. 
*/	 
template <class T> void xar::XArray3DMove<T>::Move(long lngMoveZPoints, long lngMoveYPoints, long lngMoveXPoints, T tFillVal)
{
	if (lngMoveZPoints == 0 && lngMoveYPoints == 0 && lngMoveXPoints == 0) 
		return;

	// head does not move as the 'viewport' remains constant
	
	if (fabs(lngMoveZPoints) >= m_rXArray3D.GetDim1() || fabs(lngMoveYPoints) >= m_rXArray3D.GetDim2() || fabs(lngMoveXPoints) >= m_rXArray3D.GetDim3())
	{
		Fill(Norm(eNormAver));
		return;
	}

	///@@@@@@@@@@@@@@@@
	index_t i, j;

	if (lngMoveYPoints>=0 && lngMoveXPoints>=0) 
	{
		for (i=m_rXArray3D.GetDim1()-lngMoveYPoints-1; i>=0; i--)
			for (j=m_rXArray3D.GetDim2()-lngMoveXPoints-1; j>=0; j--)
				m_rXArray3D[i+lngMoveYPoints][j+lngMoveXPoints] = m_rXArray3D[i][j];
		for (i=0; i<lngMoveYPoints; i++) 
			for (j=0; j<m_rXArray3D.GetDim2(); j++) 
				m_rXArray3D[i][j] = tFillVal;
		for (i=lngMoveYPoints; i<m_rXArray3D.GetDim1(); i++) 
			for (j=0; j<lngMoveXPoints; j++)
				m_rXArray3D[i][j] = tFillVal;				
	}
	else if (lngMoveYPoints>=0 && lngMoveXPoints<0) 
	{
		lngMoveXPoints = -lngMoveXPoints;
		for (i=m_rXArray3D.GetDim1()-lngMoveYPoints-1; i>=0; i--)
			for (j=0; j<m_rXArray3D.GetDim2()-lngMoveXPoints; j++)
				m_rXArray3D[i+lngMoveYPoints][j] = m_rXArray3D[i][j+lngMoveXPoints];
		for (i=0; i<lngMoveYPoints; i++) 
			for (j=0; j<m_rXArray3D.GetDim2(); j++) 
				m_rXArray3D[i][j] = tFillVal;
		for (i=lngMoveYPoints; i<m_rXArray3D.GetDim1(); i++) 
			for (j=m_rXArray3D.GetDim2()-lngMoveXPoints; j<m_rXArray3D.GetDim2(); j++) 
				m_rXArray3D[i][j] = tFillVal;
	}
	else if (lngMoveYPoints<0 && lngMoveXPoints>=0) 
	{
		lngMoveYPoints = -lngMoveYPoints;
		for (i=0; i<m_rXArray3D.GetDim1()-lngMoveYPoints; i++)
			for (j=m_rXArray3D.GetDim2()-lngMoveXPoints-1; j>=0; j--)
				m_rXArray3D[i][j+lngMoveXPoints] = m_rXArray3D[i+lngMoveYPoints][j];
		for (i=m_rXArray3D.GetDim1()-lngMoveYPoints; i<m_rXArray3D.GetDim1(); i++) 
			for (j=0; j<m_rXArray3D.GetDim2(); j++) 
				m_rXArray3D[i][j] = tFillVal;
		for (i=0; i<m_rXArray3D.GetDim1()-lngMoveYPoints; i++) 
			for (j=0; j<lngMoveXPoints; j++) 
				m_rXArray3D[i][j] = tFillVal;
	}
	else if (lngMoveYPoints<0 && lngMoveXPoints<0) 
	{
		lngMoveYPoints = -lngMoveYPoints; lngMoveXPoints = -lngMoveXPoints;
		for (i=0; i<m_rXArray3D.GetDim1()-lngMoveYPoints; i++)
			for (j=0; j<m_rXArray3D.GetDim2()-lngMoveXPoints; j++)
				m_rXArray3D[i][j] = m_rXArray3D[i+lngMoveYPoints][j+lngMoveXPoints];
		for (i=m_rXArray3D.GetDim1()-lngMoveYPoints; i<m_rXArray3D.GetDim1(); i++) 
			for (j=0; j<m_rXArray3D.GetDim2(); j++)
				m_rXArray3D[i][j] = tFillVal;
		for (i=0; i<m_rXArray3D.GetDim1()-lngMoveYPoints; i++) 
			for (j=m_rXArray3D.GetDim2()-lngMoveXPoints; j<m_rXArray3D.GetDim2(); j++)
				m_rXArray3D[i][j] = tFillVal;
	}
}
#endif


//---------------------------------------------------------------------------
//Function XArray3DMove<T>::FillRect
//
//	Fills a rectangular subarray with a given value
//
/*!
	\brief		Fills a rectangular subarray with a given value
	\param		LowDim1	lower Z boundary of the rectangle to be filled
	\param		HighDim1	upper Z boundary of the rectangle to be filled
	\param		LowDim2	lower Y boundary of the rectangle to be filled
	\param		HighDim2	upper Y boundary of the rectangle to be filled
	\param		LowDim3	lower X boundary of the rectangle to be filled
	\param		HighDim3	upper X boundary of the rectangle to be filled
	\param		tFillVal	the value to be assigned to all elements inside the rectangle
	\exception	std::invalid_argument is thrown if any rectangle boundary is negative or
				exceeds the	corresponding array boundary, or if lower rectangle boundary is
				larger or equal to the higher rectangle boundary
	\exception	std::exception and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
				This function fills a defined rectangular region inside a 3D array with a
				given value. Lower boundaries are included in the filled rectangle, while the
				upper boundaries are excluded; the head is not affected
*/	
template <class T> void xar::XArray3DMove<T>::FillRect(index_t LowDim1, index_t HighDim1, index_t LowDim2, index_t HighDim2, index_t LowDim3, index_t HighDim3, T tFillVal)
{
	index_t i, j, k, nn = 0;

	if ( HighDim1 > m_rXArray3D.GetDim1() || LowDim1 >= HighDim1)
		throw std::invalid_argument("invalid_argument 'LowDim1 or HighDim1' in XArray3DMove<T>::FillRect"); 
	if (HighDim2 > m_rXArray3D.GetDim2() || LowDim2 >= HighDim2)
		throw std::invalid_argument("invalid_argument 'LowDim2 or HighDim2' in XArray3DMove<T>::FillRect"); 
	if (HighDim3 > m_rXArray3D.GetDim3() || LowDim3 >= HighDim3)
		throw std::invalid_argument("invalid_argument 'LowDim3 or HighDim3' in XArray3DMove<T>::FillRect"); 

	for (i = LowDim1; i < HighDim1; i++)
		for (j = LowDim2; j < HighDim2; j++)
			for (k = LowDim3; k < HighDim3; k++)
				m_rXArray3D[i][j][k] = tFillVal; 
}


template <class T> void xar::XArray3DMove<T>::FillRectPeriodic(index_t k0, index_t j0, index_t i0, index_t karad, index_t jarad, index_t iarad, T tFillVal)
// Fills a 3D vicinity of a point with tFillVal; if the point is near the edge, the filling area may extend periodically across the boundaries
{
	index_t nz = m_rXArray3D.GetDim1(), ny = m_rXArray3D.GetDim2(), nx = m_rXArray3D.GetDim3();
	double dnz = double(nz), dny = double(ny), dnx = double(nx);
	
	int kk1, jj1, ii1;
	for (int kk = int(k0) - int(karad); kk <= int(k0 + karad); kk++)
	{
		kk1 = nmodm(kk, dnz);
		for (int jj = int(j0) - int(jarad); jj <= int(j0 + jarad); jj++)
		{
			jj1 = nmodm(jj, dny);
			for (int ii = int(i0) - int(iarad); ii <= int(i0 + iarad); ii++)
			{
				ii1 = nmodm(ii, dnx);
				m_rXArray3D[kk1][jj1][ii1] = tFillVal;
			}
		}
	}
}


template <class T> void xar::XArray3DMove<T>::FillCylinderPeriodic(index_t k0, index_t j0, index_t i0, index_t karad, index_t jarad, index_t iarad, T tFillVal)
// Fills a 3D vicinity of a point with tFillVal; if the point is near the edge, the filling area may extend periodically across the boundaries
{
	index_t nz = m_rXArray3D.GetDim1(), ny = m_rXArray3D.GetDim2(), nx = m_rXArray3D.GetDim3();
	double dnz = double(nz), dny = double(ny), dnx = double(nx);
	int ijrad2 = int((jarad * jarad + iarad * iarad) / 2.0 + 0.5);
	int kk1, jj1, ii1, jdist2, idist2;
	int kk0 = int(k0), kkarad = int(karad), jj0 = int(j0), jjarad = int(jarad), ii0 = int(i0), iiarad = int(iarad);

	for (int kk = kk0 - kkarad; kk <= kk0 + kkarad; kk++)
	{
		kk1 = nmodm(kk, dnz);
		for (int jj = jj0 - jjarad; jj <= jj0 + jjarad; jj++)
		{
			jj1 = nmodm(jj, dny);
			jdist2 = (jj - jj0) * (jj - jj0);
			for (int ii = ii0 - iiarad; ii <= ii0 + iiarad; ii++)
			{
				idist2 = (ii - ii0) * (ii - ii0);
				if (idist2 + jdist2 <= ijrad2)
				{
					ii1 = nmodm(ii, dnx);
					m_rXArray3D[kk1][jj1][ii1] = tFillVal;
				}
			}
		}
	}
}


template <class T> void xar::XArray3DMove<T>::FillRectComplementPeriodic(index_t k0, index_t j0, index_t i0, index_t karad, index_t jarad, index_t iarad, T tFillVal)
// Fills all points OUTSIDE the defined 3D vicinity of a point with tFillVal; if the point is near the edge, the non-filled area may extend periodically across the boundaries
{
	index_t nz = m_rXArray3D.GetDim1(), ny = m_rXArray3D.GetDim2(), nx = m_rXArray3D.GetDim3();
	double dnz = double(nz), dny = double(ny), dnx = double(nx);

	// create the set of indexes corresponding to the elements that need to be filled with tFillVal by knocking out the indexes of elements that should not be touched
	int i, j, k;
	vector<int> kindexes(nz), jindexes(ny), iindexes(nx);
	vector<int>::iterator kk0 = kindexes.begin(), jj0 = jindexes.begin(), ii0 = iindexes.begin(), kk, jj, ii;
	
	for (k = 0; k < nz; k++) kindexes[k] = k;
	for (j = 0; j < ny; j++) jindexes[j] = j;
	for (i = 0; i < nx; i++) iindexes[i] = i;
	
	for (k = int(k0) - int(karad); k <= int(k0 + karad); k++) kindexes[nmodm(k, dnz)] = -1;
	for (j = int(j0) - int(jarad); j <= int(j0 + jarad); j++) jindexes[nmodm(j, dny)] = -1;
	for (i = int(i0) - int(iarad); i <= int(i0 + iarad); i++) iindexes[nmodm(i, dnx)] = -1;
	
	// fill the appropriate elements with tFillVal
	for (kk = kindexes.begin(); kk != kindexes.end(); kk++)
	{
		if ((k = *kk) < 0) continue;
		for (j = 0; j < ny; j++)
			for (i = 0; i < nx; i++)
				m_rXArray3D[k][j][i] = tFillVal;
	}
	for (jj = jindexes.begin(); jj != jindexes.end(); jj++)
	{
		if ((j = *jj) < 0) continue;
		for (k = 0; k < nz; k++)
			for (i = 0; i < nx; i++)
				m_rXArray3D[k][j][i] = tFillVal;
	}
	for (ii = iindexes.begin(); ii != iindexes.end(); ii++)
	{
		if ((i = *ii) < 0) continue;
		for (k = 0; k < nz; k++)
			for (j = 0; j < ny; j++)
				m_rXArray3D[k][j][i] = tFillVal;
	}
}


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
	template class xar::XArray3DMove<char>;
	template class xar::XArray3DMove<short>;
	template class xar::XArray3DMove<long>;
	template class xar::XArray3DMove<float>;
	template class xar::XArray3DMove<double>;
	template class xar::XArray3DMove<xar::fcomplex>;
	template class xar::XArray3DMove<xar::dcomplex>;

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
