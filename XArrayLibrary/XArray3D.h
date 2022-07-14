//Header XArray3D.h
//
//
//	HEADER FILE TITLE:
//
//		Three dimensional XArray class
//
/*!
	\file		XArray3D.h
	\brief		Three dimensional XArray class
	\par		Description:
		This is a 3D version of the XArray class which implements the notion of a three-dimensional 
		resizable (dynamic) array with an optional head
*/
#if !defined XARRAY3D_H
#define XARRAY3D_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XArray.h"

//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
	namespace xar { template <class T> class XArray3D; }
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
namespace // anonymous
{
//---------------------------------------------------------------------------
//Class Index2D
//
//	Auxiliary class which is returned by the XArray3D<T>::operator[]
//
/*!
	\brief		Auxiliary class which is returned by the XArray3D<T>::operator[]
	\par		Description:
				This class is used for implementing the natural A[i][j][k] access to elements of XArray3D<T> by
				making A[i] an Index2D<T> object for which the operator[] is defined as in XArray2D<T> 
	\remarks	This class is intended for use only inside XArray3D<T> definition
*/
	template <class T> class Index2D
	{
	// Constructors
	public:
		//! Constructs modifiable Index2D<T> object corresponding to (*pXArray3D)[i]
		Index2D(xar::XArray3D<T>* pXArray3D, index_t i) : m_pXArray3D(pXArray3D), m_p_i_dim2_dim3(&(pXArray3D->front()) + i * pXArray3D->m_iDim2 * pXArray3D->m_iDim3) {}

	// Operators 
	public:
		//! Returns (*m_pXArray3D)[i][j], i.e. a modifiable pointer to the beginning of row (i, j)
		T* operator[] (index_t j) { return m_p_i_dim2_dim3 + j * m_pXArray3D->m_iDim3; }

	private:
	// Member variables	
		//! Pointer to the host XArray3D<T>
		xar::XArray3D<T>* m_pXArray3D;
		//! Pointer to the beginning of the plain (i), i.e. m_pXArray3D[i * pXArray3D->m_iDim2 * pXArray3D->m_iDim3]
		T* m_p_i_dim2_dim3; 
	};

//---------------------------------------------------------------------------
//Class constIndex2D
//
//	Auxiliary class which is returned by the XArray3D<T>::operator[] const
//
/*!
	\brief		Auxiliary class which is returned by the XArray3D<T>::operator[] const
	\par		Description:
				This class is used for implementing the natural, A[i][j][k], read-only access to elements
				of XArray3D<T> by making A[i] an constIndex2D<T> object for which the operator[] is defined
				as in XArray2D<T> 
	\remarks	This class is intended for use only inside XArray3D<T> definition
*/
	template <class T> class constIndex2D
	{
	// Constructors
	public:
		//! Constructs constIndex2D<T> object corresponding to the read-only pointer (*pXArray3D)[i] 
		constIndex2D(const xar::XArray3D<T>* const pXArray3D, index_t i) : m_pXArray3D(pXArray3D), m_p_i_dim2_dim3(&(pXArray3D->front()) + i * pXArray3D->m_iDim2 * pXArray3D->m_iDim3) {}

	// Operators 
	public:
		//! Returns (*m_pXArray3D)[i][j], i.e. a read-only pointer to the beginning of row (i, j)
		const T* operator[](index_t j) const { return m_p_i_dim2_dim3 + j * m_pXArray3D->m_iDim3; }

	private:
	// Member variables	
		//! Pointer to the host const XArray3D<T>
		const xar::XArray3D<T>* const m_pXArray3D;
		//! Pointer to the beginning of the plain (i) in XArray3D<T>, i.e. m_pXArray3D[i * pXArray3D->m_iDim2 * pXArray3D->m_iDim3]
		const T* m_p_i_dim2_dim3;
	};
} // end of anonymous namespace


namespace xar
{
//---------------------------------------------------------------------------
//Class XArray3D<T>
//
//	Three dimensional XArray class
//
/*!
	\brief		Three dimensional XArray class
	\par		Description:
				This class implements the notion of a three-dimensional resizable (dynamic) array with an optional head.
				Much of the class's functionality is implemented in XArray<T>.
	\remarks	This class can be instantiated (unlike its parent XArray<T>)
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
	\remarks    Functions that may require additional information about the head (in excess of what is available in the 
				IXAHead interface) are placed in separate classes (e.g. XArray3DMove<T>).
	\remarks	The C storage order is used for the array, i.e. the last index changes the fastest.
	\warning	The absence of appropriate specializations of the function GetValuetype() in the base class XArray<T>
				will prevent instantiation of XArray3D<T> object for new types T
*/
	template <class T> class XArray3D : public XArray<T>
	{
	// Friends
		friend class Index2D<T>;
		friend class constIndex2D<T>;

	// Enumerators
	// Structures
	// Constructors
	public:
		//! Default constructor
		XArray3D(void);
		//! Constructor with predefined sizes
		XArray3D(index_t iDim1, index_t iDim2, index_t iDim3, T tVal = T());
		//! Promotion from vector
		XArray3D(const vector<T>& rXArB, index_t iDim1, index_t iDim2, index_t iDim3);
		//! Move promotion from vector
		XArray3D(vector<T>&& rXArB, index_t iDim1, index_t iDim2, index_t iDim3);
		//! Copy constructor
		XArray3D(const XArray3D<T>& rXArray3D);
		//! Move constructor
		XArray3D(XArray3D<T>&& rXArray3D);
		//! Destructor (head is deleted in the base class)
		~XArray3D() {}

	// Operators 
	// Most are inherited from XArray<T>, others need to verify that the dimensions of the arguments are the same
	public:
		//! Provides read-only access to array's rows
		constIndex2D<T> operator[] (index_t i) const { return constIndex2D<T>(this, i); }
		//! Provides read and write access to array's rows
		Index2D<T> operator[] (index_t i) { return Index2D<T>(const_cast<XArray3D<T>*>(this), i); }
		//! Makes this array a (deep) copy of the rXArray3D
		void operator=(const XArray3D<T>& rXArray3D); 
		//! Move assignment operator from another XArray3D
		void operator=(XArray3D<T>&& rXArray3D) noexcept;
		//! Performs elementwise addition of the two arrays
		void operator+=(const XArray3D<T>& rXArray3D);
		//! Performs elementwise subtraction of the two arrays
		void operator-=(const XArray3D<T>& rXArray3D);
		//! Performs elementwise multiplication of the two arrays
		void operator*=(const XArray3D<T>& rXArray3D);
		//! Performs elementwise division of the two arrays (checks for division by zero)
		void operator/=(const XArray3D<T>& rXArray3D);
		//
		// The following statements are required because the above ones hide the corresponding 
		// XArray operators, and name resolution does not work across scopes 
		// (see Stroustrup, C++ Programming Language, 3rd ed., Sect.15.2.2)
		//
//#if _MSC_VER <= 1200 // VC++ 6.0 and lower
		//! Adds a given value to each array element
		void operator+=(T tVal) { XArray<T>::operator+=(tVal); }
		//! Subtracts a given value from each array element
		void operator-=(T tVal){ XArray<T>::operator-=(tVal); }
		//! Multiplies each array element by a given value
		void operator*=(T tVal){ XArray<T>::operator*=(tVal); }
		//! Divides each array element by a given value  (checks for division by zero)
		void operator/=(T tVal){ XArray<T>::operator/=(tVal); }
/*
#else
		//! Adds a given value to each array element
		using XArray<T>::operator+=; // brings in void XArray<T>::operator+=(tVal);
		//! Subtracts a given value from each array element
		using XArray<T>::operator-=; // brings in void XArray<T>::operator-=(tVal);
		//! Multiplies each array element by a given value
		using XArray<T>::operator*=; // brings in void XArray<T>::operator*=(tVal);
		//! Divides each array element by a given value  (checks for division by zero)
		using XArray<T>::operator/=; // brings in void XArray<T>::operator+=(tVal);
#endif
*/

	// Attributes
	public:
		//! Returns the dimensionality of the XArray2D<T> class (it is always equal to eDim2, i.e. two) 
		static _eArrayDim GetArraydim(void) { return eDim3; }
		//! Returns the first dimension (corresponds to Z)
		index_t GetDim1(void) const { return m_iDim1; }
		//! Returns the second dimension (corresponds to Y)
		index_t GetDim2(void) const { return m_iDim2; }
		//! Returns the third dimension (corresponds to X)
		index_t GetDim3(void) const { return m_iDim3; }
		//! Sets all three dimensions of the 3D array (reshapes the array)
		void SetDims(index_t iDim1, index_t iDim2, index_t iDim3); //this is a low-level function
		//! Returns true if and only if all three dimensions of the two arrays are the same
		bool IsSameShape(const XArray3D<T>& xa3) const 
			{ return m_iDim1==xa3.m_iDim1 && m_iDim2==xa3.m_iDim2 && m_iDim3==xa3.m_iDim3; }

	// Operations
	public:
		//! Provides slow (but checked and polymorphic) read access to an element
		double GetAt(index_t i, index_t j, index_t k) const { return XArray<T>::GetAt(i * m_iDim2 * m_iDim3 + j * m_iDim3 + k); }
		//! Provides slow (but checked and polymorphic) write access to an element
		void SetAt(index_t i, index_t j, index_t k, double dblVal) { XArray<T>::SetAt(i * m_iDim2 * m_iDim3 + j * m_iDim3 + k, dblVal); }
		//! Provides slow (but checked and polymorphic) read access to an element
		dcomplex GetCmplAt(index_t i, index_t j, index_t k) const { return XArray<T>::GetCmplAt(i * m_iDim2 * m_iDim3 + j * m_iDim3 + k); }
		//! Provides slow (but checked and polymorphic) write access to an element
		void SetCmplAt(index_t i, index_t j, index_t k, dcomplex cxdVal) { XArray<T>::SetCmplAt(i * m_iDim2 * m_iDim3 + j * m_iDim3 + k, cxdVal); }
		//! Accepts an external memory buffer with its contents and makes it the internal 3D array (head does not change)
		//void AcceptMemBuffer(T* ptBufBegin, index_t iDim1, index_t iDim2, index_t iDim3);
		//! Relinquishes the responsibility for the memory area occupied by the internal 3D array  (head is deleted)
		//void ReleaseMemBuffer();
		//! Changes the size of the array and fills the NEW elements with a given value (head is not affected)
		//  Call the non-member function ResizeH if the head should be resized as well
		void Resize(index_t iDim1, index_t iDim2, index_t iDim3, T tVal = T());
		//! Changes the size of the array (head is not affected)
		void ResizeA(const std::vector<index_t>& rvecNewSizes);
		//! Resizes the array to zero, FREES THE MEMORY and deletes m_pHead
		void Truncate();
		//! Swaps XArray3Ds and their heads
		void Swap(XArray3D<T>& rXArray3D);
		//! Extracts a sub-XArray3D into another XArray3D (head is not affected)
		//  Call the non-member function GetSubarrayH if the head should be returned as well
		void GetSubarray(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, index_t iBeginDim3, index_t iEndDim3, XArray3D<T>& rDestSubXArray) const; 
		//! Inserts a sub-XArray3D
		void SetSubarray(const XArray3D<T>& rSrcSubXArray, index_t iBeginDim1, index_t iBeginDim2, index_t iBeginDim3); 
		//! Finds the value and position of the minimum in a 3D array
		T Min3D(index_t& kmin, index_t& jmin, index_t& imin);
		//! Finds the value and position of the maximum in a 3D array
		T Max3D(index_t& kmax, index_t& jmax, index_t& imax);
		//! Trims elements from the edges of the array
		void Trim(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight);
		//! Returns rotationally averaged absolute values (1D radial section) of the 3D array
		vector<double> RotAbsAverage(int izc, int iyc, int ixc) const;

	// Overridables
	public:
	
	// Implementation
	protected:

		
	private:
	// Member variables	
		//! First dimension (corresponds to Z)		
		index_t m_iDim1;
		//! Second dimension (corresponds to Y)
		index_t m_iDim2;
		//! Third dimension (corresponds to X)
		index_t m_iDim3;
	};

} // namespace xar closed
//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
namespace xar
{
	//the following function definitions should ideally be placed into a .cpp file and 'export'ed

	//! Default constructor
	template <class T> XArray3D<T>::XArray3D()
	{
		m_iDim1 = m_iDim2 = m_iDim3 = 0;
	}

	//! Constructor with predefined sizes
	template <class T> XArray3D<T>::XArray3D(index_t iDim1, index_t iDim2, index_t iDim3, T val) 
		: XArray<T>(iDim1 * iDim2 * iDim3, val)
	{
		if (iDim1 < 0 || iDim2 < 0 || iDim3 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1, iDim2 or iDim3' in XArray3D<T>::XArray3D (negative dimension)"); 	
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
		m_iDim3 = iDim3;
	}

	//! Promotion from vector
	template <class T> XArray3D<T>::XArray3D(const vector<T>& va, index_t iDim1, index_t iDim2, index_t iDim3)
		: XArray<T>(va)
	{
		if (iDim1 < 0 || iDim2 < 0 || iDim3 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1, iDim2 or iDim3' in XArray3D<T>::XArray3D (negative dimension)"); 	
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
		m_iDim3 = iDim3;
	}

	//! Move promotion from vector
	template <class T> XArray3D<T>::XArray3D(vector<T>&& va, index_t iDim1, index_t iDim2, index_t iDim3)
		: XArray<T>(std::move(va))
	{
		if (iDim1 < 0 || iDim2 < 0 || iDim3 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1, iDim2 or iDim3' in XArray3D<T>::XArray3D (negative dimension)");
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
		m_iDim3 = iDim3;
	}

	//! Copy constructor
	template <class T> XArray3D<T>::XArray3D(const XArray3D<T>& xa3)
		: XArray<T>(xa3)
	{
		m_iDim1 = xa3.m_iDim1;
		m_iDim2 = xa3.m_iDim2;
		m_iDim3 = xa3.m_iDim3;
	}

	//! Move constructor
	template <class T> XArray3D<T>::XArray3D(XArray3D<T>&& xa3)
		: XArray<T>(std::move(xa3))
	{
		m_iDim1 = xa3.m_iDim1;
		m_iDim2 = xa3.m_iDim2;
		m_iDim3 = xa3.m_iDim3;
	}

	//! Changes the size of the array and fills the NEW elements with a given value (head is not affected)
	// Call non-member function ResizeH if the head should be resized as well
	template <class T> void XArray3D<T>::Resize(index_t iDim1, index_t iDim2, index_t iDim3, T val)
	{
		if (iDim1 < 0 || iDim2 < 0 || iDim3 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1, iDim2 or iDim3' in XArray3D<T>::Resize (negative dimension)"); 
		(*this).resize(iDim1 * iDim2 * iDim3, val);
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
		m_iDim3 = iDim3;
	}

	//! Changes the size of the array (head is not affected)
	// Call non-member function ResizeH if the head should be resized as well
	template <class T> void XArray3D<T>::ResizeA(const std::vector<index_t>& rvecNewSizes)
	{
		if (rvecNewSizes.size() != 3)
			throw std::invalid_argument("invalid_argument 'rvecNewSizes' in XArray3D<T>::ResizeA (wrong dimensionality)"); 
		Resize(rvecNewSizes[0], rvecNewSizes[1], rvecNewSizes[2]);
	}

	//! Resizes the array to zero, FREES THE MEMORY and deletes m_pHead	
	template <class T> void XArray3D<T>::Truncate()
	{
		XArray<T>::SetHeadPtr(0);
		(*this).truncate();
		m_iDim1 = 0;
		m_iDim2 = 0;
		m_iDim3 = 0;
	}

	//! Swaps XArray3Ds and their heads
	template <class T> void XArray3D<T>::Swap(XArray3D<T>& rXArray3D)
	{
		(*this).swap(rXArray3D);
		IXAHead* temp = XArray<T>::GetHeadPtr() ? XArray<T>::GetHeadPtr()->Clone() : 0;
		IXAHead* temp1 = rXArray3D.XArray<T>::GetHeadPtr() ? rXArray3D.XArray<T>::GetHeadPtr()->Clone() : 0;
		XArray<T>::SetHeadPtr(temp1); rXArray3D.XArray<T>::SetHeadPtr(temp);
		std::swap(m_iDim1, rXArray3D.m_iDim1);
		std::swap(m_iDim2, rXArray3D.m_iDim2);
		std::swap(m_iDim3, rXArray3D.m_iDim3);
	}

	//! Sets all three dimensions of the 3D array (reshapes the array)
	// This is a low-level function which should be used sparingly
	template <class T> void XArray3D<T>::SetDims(index_t iDim1, index_t iDim2, index_t iDim3)
	{
		if (iDim1 < 0 || iDim2 < 0 || iDim3 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1, iDim2 or iDim3' in XArray3D<T>::SetDims (negative dimension)"); 
		if (iDim1 * iDim2 * iDim3 != (*this).size())
			throw std::invalid_argument("invalid_argument 'iDim1, iDim2 or iDim3' in XArray3D<T>::SetDims (iDim1*iDim2 != (*this).size())");

		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
		m_iDim3 = iDim3;
	}

	//! Extracts a sub-XArray3D into another XArray3D (head is not affected)
	//  Call the non-member function GetSubarrayH if the head should be returned as well
	template <class T> void XArray3D<T>::GetSubarray(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, index_t iBeginDim3, index_t iEndDim3, XArray3D<T>& rDestSubXArray) const
	{
		if (iBeginDim1 >= iEndDim1 || iEndDim1 > GetDim1())
			throw std::invalid_argument("invalid argument 'iBeginDim1 or iEndDim1' in XArray3D<T>::GetSubarray");
		if (iBeginDim2 >= iEndDim2 || iEndDim2 > GetDim2())
			throw std::invalid_argument("invalid argument 'iBeginDim2 or iEndDim2' in XArray3D<T>::GetSubarray");
		if (iBeginDim3 >= iEndDim3 || iEndDim3 > GetDim3())
			throw std::invalid_argument("invalid argument 'iBeginDim3 or iEndDim3' in XArray3D<T>::GetSubarray");
		index_t iSize1 = iEndDim1 - iBeginDim1;
		index_t iSize2 = iEndDim2 - iBeginDim2;
		index_t iSize3 = iEndDim3 - iBeginDim3;
		if (rDestSubXArray.GetDim1() != iSize1 || rDestSubXArray.GetDim2() != iSize2 || rDestSubXArray.GetDim3() != iSize3) 
			rDestSubXArray.Resize(iSize1, iSize2, iSize3);
		for (index_t i = 0; i < iSize1; i++) 
		{
			const T* pTemp0 = &((*this).front()) + (iBeginDim1 + i) * GetDim2() * GetDim3();
			for (index_t j = 0; j < iSize2; j++)
			{
				const T* pTemp = pTemp0  + (iBeginDim2 + j) * GetDim3() + iBeginDim3;
				for (index_t k = 0; k < iSize3; k++)
				{
					rDestSubXArray[i][j][k] = *(pTemp + k);
				}
			}
		}
	}

	//! Trims elements from the edges of the array
	template <class T> void XArray3D<T>::Trim(index_t iZLeft, index_t iZRight, index_t iYLeft, index_t iYRight, index_t iXLeft, index_t iXRight)
	{
		if (iXLeft == 0 && iXRight == 0 && iYLeft == 0 && iYRight == 0 && iZLeft == 0 && iZRight == 0)
			return;

		if (iXLeft + iXRight > m_iDim3)
			throw std::invalid_argument("invalid_argument 'iXLeft or iXRight' in XArray3D<T>::Trim");

		if (iYLeft + iYRight > m_iDim2)
			throw std::invalid_argument("invalid_argument 'iYLeft or iYRight' in XArray3D<T>::Trim");

		if (iZLeft + iZRight > m_iDim1)
			throw std::invalid_argument("invalid_argument 'iXLeft or iXRight' in XArray3D<T>::Trim");


		index_t iDim3 = m_iDim3 - (iXLeft + iXRight);
		index_t iDim2 = m_iDim2 - (iYLeft + iYRight);
		index_t iDim1 = m_iDim1 - (iZLeft + iZRight);

		XArray3D<T> xarTemp(iDim1, iDim2, iDim3);

		for (index_t k = 0; k < iDim1; k++)
			for (index_t j = 0; j < iDim2; j++)
				for (index_t i = 0; i < iDim3; i++)
					xarTemp[k][j][i] = ((*this))[k + iZLeft][j + iYLeft][i + iXLeft];

		// set the same head
		xarTemp.SetHeadPtr((*this).GetHeadPtr() ? (*this).GetHeadPtr()->Clone() : 0);

		// transform the head if it can be done
		if (GetIXAHWave3D(xarTemp))
			GetIXAHWave3D(xarTemp)->Trim(m_iDim1, m_iDim2, m_iDim3, iZLeft, iZRight, iYLeft, iYRight, iXLeft, iXRight);

		(*this).Swap(xarTemp);
	}

	//! Returns rotationally averaged absolute values (1D radial section) of the 3D array
	template <class T> vector<double> XArray3D<T>::RotAbsAverage(int kc, int jc, int ic) const
	// (ixc, iyc, izc) - voxel indexes of the centre of rotation
	// returns 1D rotationally averaged radial section of the 3D array
	{
		if (ic < 0 || ic > m_iDim3 - 1)
			throw std::invalid_argument("invalid_argument 'ixc' in XArray3D<T>::RotAverage");

		if (jc < 0 || jc > m_iDim2 - 1)
			throw std::invalid_argument("invalid_argument 'iyc' in XArray3D<T>::RotAverage");

		if (kc < 0 || kc > m_iDim1 - 1)
			throw std::invalid_argument("invalid_argument 'izc' in XArray3D<T>::RotAverage");

		index_t nn = (index_t)(sqrt(m_iDim1 * m_iDim1 + m_iDim2 * m_iDim2 + m_iDim3 * m_iDim3) + 1);
		vector<double> vRadSection(nn, 0.0);
		vector<index_t> vWeights(nn, 0);

		long int k2, j2;
		int iR, iRmax(0);
		for (int k = 0; k < m_iDim1; k++)
		{
			k2 = (k - kc) * (k - kc);
			for (int j = 0; j < m_iDim2; j++)
			{
				j2 = k2 + (j - jc) * (j - jc);
				for (int i = 0; i < m_iDim3; i++)
				{
					iR = int(sqrt(j2 + (i - ic) * (i - ic)) + 0.5);
					if (iR > iRmax) iRmax = iR;
					vRadSection[iR] += std::abs(((*this))[k][j][i]);
					vWeights[iR]++;
				}
			}
		}

		for (index_t n = 0; n <= iRmax; n++) if (vWeights[n] != 0) vRadSection[n] /= (double)vWeights[n];
		vRadSection.resize(iRmax + 1);

		return vRadSection;
	}

	
	//! Inserts a sub-XArray3D
	template <class T> void XArray3D<T>::SetSubarray(const XArray3D<T>& rSrcSubXArray, index_t iBeginDim1, index_t iBeginDim2, index_t iBeginDim3)
	{
		if (iBeginDim1 + rSrcSubXArray.GetDim1() > GetDim1())
			throw std::invalid_argument("invalid argument 'iBeginDim1 or rSrcSubXArray' in XArray3D<T>::SetSubarray");
		if (iBeginDim2 + rSrcSubXArray.GetDim2() > GetDim2())
			throw std::invalid_argument("invalid argument 'iBeginDim2 or rSrcSubXArray' in XArray3D<T>::SetSubarray");
		if (iBeginDim3 + rSrcSubXArray.GetDim3() > GetDim3())
			throw std::invalid_argument("invalid argument 'iBeginDim3 or rSrcSubXArray' in XArray3D<T>::SetSubarray");
		for (index_t i = 0; i < rSrcSubXArray.GetDim1(); i++) 
		{
			T* pTemp0 = &((*this).front()) + (iBeginDim1 + i) * GetDim2() * GetDim3();
			for (index_t j = 0; j < rSrcSubXArray.GetDim2(); j++)
			{
				T* pTemp = pTemp0  + (iBeginDim2 + j) * GetDim3() + iBeginDim3;
				for (index_t k = 0; k < rSrcSubXArray.GetDim3(); k++)
				{
					*(pTemp + k) = rSrcSubXArray[i][j][k];
				}
			}
		}
		// heads are ignored
	}

	//! Finds the value and position of the minimum in a 3D array
	template <class T> T XArray3D<T>::Min3D(index_t& kmin, index_t& jmin, index_t& imin)
	{
		kmin = jmin = imin = 0;
		T amin = (*this)[kmin][jmin][imin];
		for (index_t kk = 0; kk < m_iDim1; kk++)
			for (index_t jj = 0; jj < m_iDim2; jj++)
				for (index_t ii = 0; ii < m_iDim3; ii++)
					if ((*this)[kk][jj][ii] < amin)
					{

						amin = (*this)[kk][jj][ii];
						kmin = kk; jmin = jj; imin = ii;
					}

		return amin;
	}

	//! Finds the value and position of the maximum in a 3D array
	template <class T> T XArray3D<T>::Max3D(index_t& kmax, index_t& jmax, index_t& imax)
	{
		kmax = jmax = imax = 0;
		T amax = (*this)[kmax][jmax][imax];
		for (index_t kk = 0; kk < m_iDim1; kk++)
			for (index_t jj = 0; jj < m_iDim2; jj++)
				for (index_t ii = 0; ii < m_iDim3; ii++)
					if ((*this)[kk][jj][ii] > amax)
					{

						amax = (*this)[kk][jj][ii];
						kmax = kk; jmax = jj; imax = ii;
					}

		return amax;
	}


	//***** XArray3D member operators

	//! Makes this array a (deep) copy of the argument
	template <class T> inline void XArray3D<T>::operator=(const XArray3D<T>& xa3)
	{
		if (this == &xa3) return;
		XArray<T>::operator=(xa3);
		m_iDim1 = xa3.m_iDim1;
		m_iDim2 = xa3.m_iDim2;
		m_iDim3 = xa3.m_iDim3;
	}

	//! Move assignment operator from another XArray3D
	template <class T> inline void XArray3D<T>::operator=(XArray3D<T>&& xa3) noexcept
	{
		if (this == &xa3) return;
		m_iDim1 = xa3.m_iDim1;
		m_iDim2 = xa3.m_iDim2;
		m_iDim3 = xa3.m_iDim3;
		XArray<T>::operator=(std::move(xa3));
	}

	//! Performs elementwise addition of the two arrays
	template <class T> inline void XArray3D<T>::operator+=(const XArray3D<T>& xa3)
	{
		if (!IsSameShape(xa3))
			throw std::invalid_argument("invalid_argument 'xa3' in XArray3D<T>::+= (different shape)"); 
		XArray<T>::operator+=(xa3);
	}

	//! Performs elementwise subtraction of the two arrays
	template <class T> inline void XArray3D<T>::operator-=(const XArray3D<T>& xa3)
	{
		if (!IsSameShape(xa3))
			throw std::invalid_argument("invalid_argument 'xa3' in XArray3D<T>::-= (different shape)"); 
		XArray<T>::operator-=(xa3);
	}

	//! Performs elementwise multiplication of the two arrays
	template <class T> inline void XArray3D<T>::operator*=(const XArray3D<T>& xa3)
	{
		if (!IsSameShape(xa3))
			throw std::invalid_argument("invalid_argument 'xa3' in XArray3D<T>::*= (different shape)"); 
		XArray<T>::operator*=(xa3);
	}

	//! Performs elementwise division of the two arrays (checks for division by zero)
	template <class T> inline void XArray3D<T>::operator/=(const XArray3D<T>& xa3)
	{
		if (!IsSameShape(xa3))
			throw std::invalid_argument("invalid_argument 'xa3' in XArray3D<T>::/= (different shape)"); 
		XArray<T>::operator/=(xa3);
	}

} // namespace xar closed

//---------------------------------------------------------------------------
//	RELATED IN-LINE NON-MEMBER DEFINITIONS
//
namespace xar
{
	//---------------------------------------------------------------------------
	//Function MakeComplex
	//
	//	Makes a complex XArray3D object C = A + ib or C = A * exp(ib) from a real XArray3D object A and a scalar b
	//
	/*!
		\brief		Makes a complex XArray3D object C = A + ib or C = A * exp(ib) from a real XArray3D object A and a scalar b
		\param		A	Real XArray3D object representing real part or modulus of a complex object to be constructed
		\param		b	Real scalar value representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray3D object constructed by this function
		\param		bMakePolar if \b true, then C = A * exp(ib) is constructed, else C = A + ib is constructed
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions called from 
					inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray3D object C = A + ib or C = A * exp(ib) from a real XArray3D object and a scalar.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter A
		\par		Example:
\verbatim
XArray3D<float> A(4, 5, 6, 1.0f);
XArray3D<fcomplex> C;
MakeComplex(A, 0.0f, C, false); // or
C = MakeComplex(A, 0.0f, false);
\endverbatim
	*/
	template <class T> void MakeComplex(const XArray3D<T>& A, T b, XArray3D< std::complex<T> >& C, bool bMakePolar)
	{
		MakeComplex((const XArray<T>&)A, b, (XArray< std::complex<T> >&)C, bMakePolar);
		C.SetDims(A.GetDim1(), A.GetDim2(), A.GetDim3());
	}

	template <class T> XArray3D< std::complex<T> > MakeComplex(const XArray3D<T>& A, T b, bool bMakePolar)
	{
		XArray3D< std::complex<T> > C;
		MakeComplex(A, b, C, bMakePolar);
		return C;  // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function MakeComplex
	//
	//	Makes a complex XArray3D object C = a + iB or C = a * exp(iB) from a real XArray3D object B and a scalar a
	//
	/*!
		\brief		Makes a complex XArray3D object C = a + iB or C = a * exp(iB) from a real XArray3D object B and a scalar a
		\param		a	Real scalar value representing real part or modulus of a complex object to be constructed
		\param		B	Real XArray3D object representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray3D object constructed by this function
		\param		bMakePolar if \b true, then C = a * exp(iB) is constructed, else C = a + iB is constructed
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions called from 
					inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray3D object C = a + iB or C = a * exp(iB) from a real XArray3D object and a scalar.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter B
		\par		Example:
\verbatim
XArray3D<float> B(4, 5, 6, 1.0f);
XArray3D<fcomplex> C;
MakeComplex(1.0f, B, C, false); // or
C = MakeComplex(1.0f, B, false);
\endverbatim
	*/
	template <class T> void MakeComplex(T a, const XArray3D<T>& B, XArray3D< std::complex<T> >& C, bool bMakePolar)
	{
		MakeComplex(a, (const XArray<T>&)B, (XArray< std::complex<T> >&)C, bMakePolar);
		C.SetDims(B.GetDim1(), B.GetDim2(), B.GetDim3());
	}
	
	template <class T> XArray3D< std::complex<T> > MakeComplex(T a, const XArray3D<T>& B, bool bMakePolar)
	{
		XArray3D< std::complex<T> > C;
		MakeComplex(a, B, C, bMakePolar);
		return C;  // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function MakeComplex
	//
	//	Makes a complex XArray3D object C = A + iB or C = A * exp(iB) from 2 real XArray3D objects, A and B
	//
	/*!
		\brief		Makes a complex XArray3D object C = A + iB or C = A * exp(iB) from 2 real XArray3D objects, A and B
		\param		A	Real XArray3D object representing real part or modulus of a complex object to be constructed
		\param		B	Real XArray3D object representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray3D object constructed by this function
		\param		bMakePolar if \b true, then C = A * exp(iB) is constructed, else C = A + iB is constructed
		\exception	std::invalid_argument is thrown if A and B have differented sizes or heads
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions called 
					from inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray3D object C = A + iB or C = A * exp(iB) from 2 real XArray3D objects.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter A
		\par		Example:
\verbatim
XArray3D<float> A(4, 5, 6, 1.0f), B(4, 5, 2.0f);
XArray3D<fcomplex> C;
MakeComplex(A, B, C, false); //or
C= MakeComplex(A, B, false);
\endverbatim
	*/
	template <class T> void MakeComplex(const XArray3D<T>& A, const XArray3D<T>& B, XArray3D< std::complex<T> >& C, bool bMakePolar)
	{
		MakeComplex((const XArray<T>&)A, (const XArray<T>&)B, (XArray< std::complex<T> >&)C, bMakePolar);
		C.SetDims(A.GetDim1(), A.GetDim2(), A.GetDim3());
	}

	template <class T> XArray3D< std::complex<T> > MakeComplex(const XArray3D<T>& A, const XArray3D<T>& B, bool bMakePolar)
	{
		XArray3D< std::complex<T> > C;
		MakeComplex(A, B, C, bMakePolar);
		return C;  // rely on move assignment / constructor to avoid copying this object on return
	}
	

	//---------------------------------------------------------------------------
	//Function Re
	//
	//	Calculates the real part of a complex XArray3D object
	//
	/*!
		\brief		Calculates the real part of a complex XArray3D object
		\param		C	Input complex XArray3D object
		\param		A	Output real XArray3D object to be made equal to the real part of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the real part of every element of a complex XArray3D object and 
			copies the result into a real XArray3D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray3D<dcomplex> C(4, 5, 6, dcomplex(1.0, 2.0));
XArray3D<double> A;
Re(C, A); // or
A = Re(C);
\endverbatim
	*/
	template <class T> void Re(const XArray3D< std::complex<T> >& C, XArray3D<T>& A)
	{ 
		Re((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2(), C.GetDim3());
	}
	
	template <class T> XArray3D<T> Re(const XArray3D< std::complex<T> >& C)
	{
		XArray3D<T> A;
		Re(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}

	
	//---------------------------------------------------------------------------
	//Function Im
	//
	//	Calculates the imaginary part of a complex XArray3D object
	//
	/*!
		\brief		Calculates the imaginary part of a complex XArray3D object
		\param		C	Input complex XArray3D object
		\param		A	Output real XArray3D object to be made equal to the imaginary part of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the imaginary part of every element of a complex XArray3D object
			and copies the result into a real XArray3D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray3D<dcomplex> C(4, 5, 6, dcomplex(1.0, 2.0));
XArray3D<double> A;
Im(C, A); // or
C = Im(A);
\endverbatim
	*/
	template <class T> void Im(const XArray3D< std::complex<T> >& C, XArray3D<T>& A)
	{ 
		Im((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2(), C.GetDim3());
	}
	
	template <class T> XArray3D<T> Im(const XArray3D< std::complex<T> >& C)
	{
		XArray3D<T> A;
		Im(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}

	//---------------------------------------------------------------------------
	//Function Abs
	//
	//	Calculates the modulus of a complex XArray3D object
	//
	/*!
		\brief		Calculates the modulus of a complex XArray3D object
		\param		C	Input complex XArray3D object
		\param		A	Output real XArray3D object to be made equal to the modulus of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the modulus of every element of a complex XArray3D object
			and copies the result into a real XArray3D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray3D<dcomplex> C(4, 5, 6, dcomplex(1.0, 2.0));
XArray3D<double> A; 
Abs(C, A); //or
A = Abs(C);
\endverbatim
	*/	
	template <class T> void Abs(const XArray3D< std::complex<T> >& C, XArray3D<T>& A)
	{ 
		Abs((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2(), C.GetDim3());
	}
	
	template <class T> XArray3D<T> Abs(const XArray3D< std::complex<T> >& C)
	{
		XArray3D<T> A;
		Abs(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function Arg
	//
	//	Calculates the argument of a complex XArray3D object
	//
	/*!
		\brief		Calculates the argument of a complex XArray3D object
		\param		C	Input complex XArray3D object
		\param		A	Output real XArray3D object to be made equal to the argument of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the argument of every element of a complex XArray3D object and copies
			the result into a real XArray3D object. The head for the output object A is copied from the input object C
		\par		Example:
\verbatim
XArray3D<dcomplex> C(4, 5, 6, dcomplex(1.0, 2.0));
XArray3D<double> A;
Arg(C, A); //or
A = Arg(C);
\endverbatim
	*/	
	template <class T> void Arg(const XArray3D< std::complex<T> >& C, XArray3D<T>& A)
	{ 
		Arg((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2(), C.GetDim3());
	}
	
	template <class T> XArray3D<T> Arg(const XArray3D< std::complex<T> >& C)
	{
		XArray3D<T> A;
		Arg(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function Abs2
	//
	//	Calculates the modulus of a complex XArray3D object
	//
	/*!
		\brief		Calculates the square modulus of a complex XArray3D object
		\param		C	Input complex XArray3D object
		\param		A	Output real XArray3D object to be made equal to the square modulus of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the square modulus of every element of a complex XArray3D object
			and copies the result into a real XArray3D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray3D<dcomplex> C(4, 5, 6, dcomplex(1.0, 2.0));
XArray3D<double> A;
Abs2(C, A); //or
A = Abs2(C);
\endverbatim
	*/	
	template <class T> void Abs2(const XArray3D< std::complex<T> >& C, XArray3D<T>& A)
	{ 
		Abs2((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2(), C.GetDim3());
	}
	
	template <class T> XArray3D<T> Abs2(const XArray3D< std::complex<T> >& C)
	{
		XArray3D<T> A;
		Abs2(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function CArg
	//
	//	Calculates the 1D-continuous phase of a complex XArray3D object
	//
	/*!
		\brief		Calculates the 1D-continuous phase of a complex XArray3D object
		\param		C	Input complex XArray3D object
		\param		A	Output real XArray3D object to be made equal to the phase of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the argument of every element of a complex XArray3D object using the 
			value of the argument of the preceding element to calculate the appropriate 2pi-multiple shift,
			and copies the result into a real XArray3D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray3D<dcomplex> C(4, 5, 6, dcomplex(1.0, 2.0));
XArray3D<double> A;
CArg(C, A); //or
A = CArg(C);
\endverbatim
	*/	
	template <class T> void CArg(const XArray3D< std::complex<T> >& C, XArray3D<T>& A) 
	{
		CArg((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2(), C.GetDim3());
	}

	template <class T> XArray3D<T> CArg(const XArray3D< std::complex<T> >& C)
	{
		XArray3D<T> A;
		CArg(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}


} //namespace xar closed

//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray3D<char>;
	template class xar::XArray3D<short>;
	template class xar::XArray3D<long>;
	template class xar::XArray3D<float>;
	template class xar::XArray3D<double>;
	template class xar::XArray3D<xar::fcomplex>;
	template class xar::XArray3D<xar::dcomplex>;
#endif

//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XARRAY3D_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
