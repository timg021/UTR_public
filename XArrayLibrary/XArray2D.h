//Header XArray2D.h
//
//
//	HEADER FILE TITLE:
//
//		Two-dimensional XArray class
//
/*!
	\file		XArray2D.h
	\brief		Two-dimensional XArray class
	\par		Description:
		This is a 2D version of the XArray class which implements the notion of a two-dimensional 
		resizable (dynamic) array with an optional head
*/
#if !defined XARRAY2D_H
#define XARRAY2D_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XArray.h"
#include "XArray1D.h"

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
//Class XArray2D<T>
//
//	Two-dimensional XArray class
//
/*!
	\brief		Two-dimensional XArray class
	\par		Description:
				This class implements the notion of a two-dimensional resizable (dynamic) array with an optional head.
				Much of the class's functionality is implemented in XArray<T>.
	\remarks	This class can be instantiated (unlike its parent XArray<T>)
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
	\remarks    Functions that may require additional information about the head (in excess of what is available in the 
				IXAHead interface) are placed in separate classes (e.g. XArray2DMove<T>).
	\remarks	The C storage order is used for the array, i.e. the last index changes the fastest.
	\warning	The absence of appropriate specializations of the function GetValuetype() in the base class XArray<T>
				will prevent instantiation of XArray2D<T> object for new types T
*/
	template <class T> class XArray2D : public XArray<T>
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Default constructor
		XArray2D(void);
		//! Constructor with predefined sizes
		XArray2D(index_t iDim1, index_t iDim2, T tVal = T());
		//! Promotion from vector
		XArray2D(const vector<T>& rXArB, index_t iDim1, index_t iDim2);
		//! Promotion from XArray
		XArray2D(const XArray<T>& rXArB, index_t iDim1, index_t iDim2);
		//! Move promotion from XArray
		XArray2D(XArray<T>&& rXArB, index_t iDim1, index_t iDim2);
		//! Copy constructor
		XArray2D(const XArray2D<T>& rXArray2D);
		//! Move constructor
		XArray2D(XArray2D<T>&& rXArray2D);
		//! Destructor (head is deleted in the base class)
		~XArray2D() {}

	// Operators 
	// Most are inherited from XArray<T>, others need to verify that the dimensions of the arguments are the same
	public:
		//! Provides read-only access to array's rows
		const T* operator[] (index_t i) const { return &((*this).front()) + i * m_iDim2; }
		//! Provides read and write access to array's rows
		T* operator[] (index_t i) { return &((*this).front()) + i * m_iDim2; }
		//! Makes this array a (deep) copy of the rXArray2D
		void operator=(const XArray2D<T>& rXArray2D); 
		//! Move assignment operator from another XArray2D
		void operator=(XArray2D<T>&& rXArray2D) noexcept;
		//! Performs elementwise addition of the two arrays
		void operator+=(const XArray2D<T>& rXArray2D);
		//! Performs elementwise subtraction of the two arrays
		void operator-=(const XArray2D<T>& rXArray2D);
		//! Performs elementwise multiplication of the two arrays
		void operator*=(const XArray2D<T>& rXArray2D);
		//! Performs elementwise division of the two arrays (checks for division by zero)
		void operator/=(const XArray2D<T>& rXArray2D);
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
		static _eArrayDim GetArraydim(void) { return eDim2; }
		//! Returns the first dimension (corresponds to Y)
		index_t GetDim1(void) const { return m_iDim1; }
		//! Returns the second dimension (corresponds to X)
		index_t GetDim2(void) const { return m_iDim2; }
		//! Sets both dimensions of the 2D array (reshapes the array)
		void SetDims(index_t iDim1, index_t iDim2);
		//! Returns true if and only if both dimensions of the two arrays are the same
		bool IsSameShape(const XArray2D<T>& xa2) const 
			{ return m_iDim1==xa2.m_iDim1 && m_iDim2==xa2.m_iDim2; }

	// Operations
	public:
		//! Provides slow (but checked and polymorphic) read access to an element
		double GetAt(index_t i, index_t j) const { return XArray<T>::GetAt(i * m_iDim2 + j); }
		//! Provides slow (but checked and polymorphic) write access to an element
		void SetAt(index_t i, index_t j, double dblVal) { XArray<T>::SetAt(i * m_iDim2 + j, dblVal); }
		//! Provides slow (but checked and polymorphic) read access to an element
		dcomplex GetCmplAt(index_t i, index_t j) const { return XArray<T>::GetCmplAt(i * m_iDim2 + j); }
		//! Provides slow (but checked and polymorphic) write access to an element
		void SetCmplAt(index_t i, index_t j, dcomplex cxdVal) { XArray<T>::SetCmplAt(i * m_iDim2 + j, cxdVal); }
		//! Accepts an external memory buffer with its contents and makes it the internal 2D array (head does not change)
		//void AcceptMemBuffer(T* ptBufBegin, index_t iDim1, index_t iDim2);
		//! Relinquishes the responsibility for the memory area occupied by the internal 2D array  (head is deleted)
		//void ReleaseMemBuffer();
		//! Changes the size of the array and fills the NEW elements with a given value (head is not affected)
		//  Call the non-member function ResizeH if the head should be resized as well
		void Resize(index_t iDim1, index_t iDim2, T tVal = T());
		//! Changes the size of the array (head is not affected)
		void ResizeA(const std::vector<index_t>& rvecNewSizes);
		//! Resizes the array to zero, FREES THE MEMORY and deletes m_pHead
		void Truncate();
		//! Swaps XArray2Ds and their heads
		void Swap(XArray2D<T>& rXArray2D);
		//! Returns the average value of elements of the array within the boundary layer with the width equal to iNumPoints
		T NormAverEdge(index_t iNumPoints) const;
		//! Extracts a sub-XArray2D into another XArray2D (head is not affected)
		//  Call the non-member function GetSubarrayH if the head should be returned as well
		void GetSubarray(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray2D<T>& rDestXArray) const; 
		//! Inserts a sub-XArray2D
		void SetSubarray(const XArray2D<T>& rSrcSubXArray, index_t iBeginDim1, index_t iBeginDim2); 
		//! Returns a cross-section along Y averaged over several columns
		//  Call the non-member function GetThickYSectionH if the head should be returned as well
		void GetThickYSection(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray1D<T>& rDestXArray1D) const; 
		//! Returns a cross-section along X averaged over several rows
		//  Call the non-member function GetThickXSectionH if the head should be returned as well
		void GetThickXSection(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray1D<T>& rDestXArray1D) const; 
		//! Shuffles XArray2D in the way that is required by most DFT routines
		// NOTE: Shuffle() is an involution, i.e. twice Shuffle() is equal to identity transform
		void Shuffle();

	// Overridables
	public:

	// Static functions
	public: 

	// Implementation
	protected:

	private:
	// Member variables	
		//! First dimension (corresponds to Y)
		index_t m_iDim1; 
		//! Second dimension (corresponds to X)
		index_t m_iDim2;

	// Member functions
	};

} // namespace xar closed
//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
namespace xar
{
	// The following function definitions should ideally be placed into a .cpp file and 'export'ed

	//! Default constructor
	template <class T> XArray2D<T>::XArray2D()
	{
		m_iDim1 = m_iDim2 = 0;
	}

	//! Constructor with predefined sizes
	template <class T> XArray2D<T>::XArray2D(index_t iDim1, index_t iDim2, T val) 
		: XArray<T>(iDim1 * iDim2, val)
	{
		if (iDim1 < 0 || iDim2 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::XArray2D (negative dimension)"); 	
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
	}

	//! Promotion from vector
	template <class T> XArray2D<T>::XArray2D(const vector<T>& va, index_t iDim1, index_t iDim2)
		: XArray<T>(va)
	{
		if (iDim1 < 0 || iDim2 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::XArray2D (negative dimension)"); 	
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
	}

	//! Promotion from XArray
	template <class T> XArray2D<T>::XArray2D(const XArray<T>& xa, index_t iDim1, index_t iDim2)
		: XArray<T>(xa)
	{
		if (iDim1 < 0 || iDim2 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::XArray2D (negative dimension)");
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
	}

	//! Move promotion from XArray
	template <class T> XArray2D<T>::XArray2D(XArray<T>&& xa, index_t iDim1, index_t iDim2)
		: XArray<T>(xa)
	{
		if (iDim1 < 0 || iDim2 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::XArray2D (negative dimension)");
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
	}

	//! Copy constructor
	template <class T> XArray2D<T>::XArray2D(const XArray2D<T>& xa2)
		: XArray<T>(xa2)
	{
		m_iDim1 = xa2.m_iDim1;
		m_iDim2 = xa2.m_iDim2;
	}

	//! Move constructor
	template <class T> XArray2D<T>::XArray2D(XArray2D<T>&& xa2)
		: XArray<T>(std::move(xa2))
	{
		m_iDim1 = xa2.m_iDim1;
		m_iDim2 = xa2.m_iDim2;
	}

	//! Accepts an external memory buffer with its contents and makes it the internal 2D array (head does not change)
	//template <class T> void XArray2D<T>::AcceptMemBuffer(T* ptBufBegin, index_t iDim1, index_t iDim2)
	//{
	//	if (iDim1 < 0 || iDim2 < 0)
	//		throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::AcceptMemBuffer (negative dimension)"); 
	//	(*this).acceptMemBuffer(ptBufBegin, iDim1 * iDim2);
	//	m_iDim1 = iDim1;
	//	m_iDim2 = iDim2;
	//}

	//! Relinquishes the responsibility for the memory area occupied by the internal 2D array  (head is deleted)
	//template <class T> void XArray2D<T>::ReleaseMemBuffer()
	//{
	//	(*this).releaseMemBuffer(); 
	//	XArray<T>::SetHeadPtr(0);
	//	m_iDim1 = 0;
	//	m_iDim2 = 0;
	//}

	//! Changes the size of the array and fills the NEW elements with a given value (head is not affected)
	// Call non-member function ResizeH if the head should be resized as well
	template <class T> void XArray2D<T>::Resize(index_t iDim1, index_t iDim2, T val)
	{
		if (iDim1 < 0 || iDim2 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::Resize (negative dimension)"); 
		(*this).resize(iDim1 * iDim2, val);
		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
	}

	//! Changes the size of the array (head is not affected)
	// Call non-member function ResizeH if the head should be resized as well
	template <class T> void XArray2D<T>::ResizeA(const std::vector<index_t>& rvecNewSizes)
	{
		if (rvecNewSizes.size() != 2)
			throw std::invalid_argument("invalid_argument 'rvecNewSizes' in XArray2D<T>::ResizeA (wrong dimensionality)"); 
		Resize(rvecNewSizes[0], rvecNewSizes[1]);
	}

	//! Resizes the array to zero, FREES THE MEMORY and deletes m_pHead
	template <class T> void XArray2D<T>::Truncate()
	{
		XArray<T>::SetHeadPtr(0);
		(*this).truncate();
		m_iDim1 = 0;
		m_iDim2 = 0;
	}

	//! Swaps XArray2Ds and their heads
	template <class T> void XArray2D<T>::Swap(XArray2D<T>& rXArray2D)
	{
		(*this).swap(rXArray2D);
		IXAHead* temp = XArray<T>::GetHeadPtr() ? XArray<T>::GetHeadPtr()->Clone() : 0;
		IXAHead* temp1 = rXArray2D.XArray<T>::GetHeadPtr() ? rXArray2D.XArray<T>::GetHeadPtr()->Clone() : 0;
		XArray<T>::SetHeadPtr(temp1); rXArray2D.XArray<T>::SetHeadPtr(temp);
		std::swap(m_iDim1, rXArray2D.m_iDim1);
		std::swap(m_iDim2, rXArray2D.m_iDim2);
	}

    //! Sets both dimensions of the 2D array (reshapes the array)
	// This is a low-level function which should be used sparingly
	template <class T> void XArray2D<T>::SetDims(index_t iDim1, index_t iDim2)
	{
		if (iDim1 < 0 || iDim2 < 0)
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::SetDims (negative dimension)"); 
		if (iDim1 * iDim2 != (*this).size())
			throw std::invalid_argument("invalid_argument 'iDim1 or iDim2' in XArray2D<T>::SetDims (iDim1*iDim2 != (*this).size())");

		m_iDim1 = iDim1;
		m_iDim2 = iDim2;
	}

	//! Returns the average value of elements of the array within the boundary layer with the width equal to iNumPoints
	template <class T> T XArray2D<T>::NormAverEdge(index_t iNumPoints) const
	{
		if (iNumPoints > GetDim1() || iNumPoints > GetDim2())
			throw std::invalid_argument("invalid argument 'iNumPoints' in XArray2D<T>::NormAverEdge");
		T tAver(0);
		
		for (index_t i = 0; i < GetDim1(); i++)
			for (index_t j = 0; j < iNumPoints; j++)
				tAver += (*this)[i][j];
		
		for (index_t i = 0; i < GetDim1(); i++)
			for (index_t j = GetDim2() - iNumPoints; j < GetDim2(); j++)
				tAver += (*this)[i][j];
		
		for (index_t i = 0; i < iNumPoints; i++)
			for (index_t j = 0; j < GetDim2(); j++)
				tAver += (*this)[i][j];
		
		for (index_t i = GetDim1() - iNumPoints; i < GetDim1(); i++)
			for (index_t j = 0; j < GetDim2(); j++)
				tAver += (*this)[i][j];

		tAver /= T(2 * GetDim1() * iNumPoints + 2 * GetDim2() * iNumPoints);

		return tAver;
	}

	//! Extracts a sub-XArray2D into another XArray2D (head is not affected)
	//  Call the non-member function GetSubarrayH if the head should be returned as well
	template <class T> void XArray2D<T>::GetSubarray(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray2D<T>& rDestXArray) const
	{
		if (iBeginDim1 >= iEndDim1 || iEndDim1 > GetDim1())
			throw std::invalid_argument("invalid argument 'iBeginDim1 or iEndDim1' in XArray2D<T>::GetSubarray");
		if (iBeginDim2 >= iEndDim2 || iEndDim2 > GetDim2())
			throw std::invalid_argument("invalid argument 'iBeginDim2 or iEndDim2' in XArray2D<T>::GetSubarray");
		index_t iSize1 = iEndDim1 - iBeginDim1;
		index_t iSize2 = iEndDim2 - iBeginDim2;
		if (rDestXArray.GetDim1() != iSize1 || rDestXArray.GetDim2() != iSize2) 
			rDestXArray.Resize(iSize1, iSize2);
		for (index_t i = 0; i < iSize1; i++) 
		{
			const T* pTemp = &((*this).front()) + (iBeginDim1 + i) * GetDim2() + iBeginDim2;
			for (index_t j = 0; j < iSize2; j++)
				rDestXArray[i][j] = *(pTemp + j);
		}
	}

	//! Inserts a sub-XArray2D
	template <class T> void XArray2D<T>::SetSubarray(const XArray2D<T>& rSrcSubXArray, index_t iBeginDim1, index_t iBeginDim2)
	{
		if (iBeginDim1 + rSrcSubXArray.GetDim1() > GetDim1())
			throw std::invalid_argument("invalid argument 'iBeginDim1 or rSrcSubXArray' in XArray2D<T>::SetSubarray");
		if (iBeginDim2 + rSrcSubXArray.GetDim2() > GetDim2())
			throw std::invalid_argument("invalid argument 'iBeginDim2 or rSrcSubXArray' in XArray2D<T>::SetSubarray");
		for (index_t i = 0; i < rSrcSubXArray.GetDim1(); i++) 
		{
			T* pTemp = &((*this).front()) + (iBeginDim1 + i) * GetDim2() + iBeginDim2;
			for (index_t j = 0; j < rSrcSubXArray.GetDim2(); j++)
				*(pTemp + j) = rSrcSubXArray[i][j];
		}
		// heads are ignored
	}

	//! Returns a cross-section along Y averaged over several columns
	//  Call the non-member function GetThickYSectionH if the head should be returned as well
	template <class T> void XArray2D<T>::GetThickYSection(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray1D<T>& rDestXArray1D) const
	{
		if (iBeginDim1 >= iEndDim1 || iEndDim1 > GetDim1())
			throw std::invalid_argument("invalid argument 'iBeginDim1 or iEndDim1' in XArray2D<T>::GetThickYSection");
		if (iBeginDim2 >= iEndDim2 || iEndDim2 > GetDim2())
			throw std::invalid_argument("invalid argument 'iBeginDim2 or iEndDim2' in XArray2D<T>::GetThickYSection");
		index_t iSize1 = iEndDim1 - iBeginDim1;
		index_t iSize2 = iEndDim2 - iBeginDim2;
		if (rDestXArray1D.size() != iSize1) 
			rDestXArray1D.Resize(iSize1);
		rDestXArray1D.Fill(T(0));
		for (index_t i = 0; i < iSize1; i++)
		{
			const T* pTemp = &((*this).front()) + (iBeginDim1 + i) * GetDim2() + iBeginDim2;
 			for (index_t j = 0; j < iSize2; j++)
					rDestXArray1D[i] += *pTemp++;
		}
	}

	//! Returns a cross-section along X averaged over several rows
	//  Call the non-member function GetThickXSectionH if the head should be returned as well
	template <class T> void XArray2D<T>::GetThickXSection(index_t iBeginDim1, index_t iEndDim1, index_t iBeginDim2, index_t iEndDim2, XArray1D<T>& rDestXArray1D) const
	{
		if (iBeginDim1 >= iEndDim1 || iEndDim1 > GetDim1())
			throw std::invalid_argument("invalid argument 'iBeginDim1 or iEndDim1' in XArray2D<T>::GetThickXSection");
		if (iBeginDim2 >= iEndDim2 || iEndDim2 > GetDim2())
			throw std::invalid_argument("invalid argument 'iBeginDim2 or iEndDim2' in XArray2D<T>::GetThickXSection");
		index_t iSize1 = iEndDim1 - iBeginDim1;
		index_t iSize2 = iEndDim2 - iBeginDim2;
		if (rDestXArray1D.size() != iSize2) 
			rDestXArray1D.Resize(iSize2);
		rDestXArray1D.Fill(T(0));
		index_t nx =  GetDim2();
		for (index_t j = 0; j < iSize2; j++)
		{
			const T* pTemp = &((*this).front()) + iBeginDim1 * nx + iBeginDim2 + j;
			for (index_t i = 0; i < iSize1; i++) 
					rDestXArray1D[j] += *(pTemp + i * nx);
		}
	}


	//***** XArray2D member operators

	//! Makes this array a (deep) copy of the argument
	template <class T> inline void XArray2D<T>::operator=(const XArray2D<T>& xa2)
	{
		if (this == &xa2) return;
		XArray<T>::operator=(xa2);
		m_iDim1 = xa2.m_iDim1;
		m_iDim2 = xa2.m_iDim2;
	}

	//! Move assignment operator from another XArray2D
	template <class T> inline void XArray2D<T>::operator=(XArray2D<T>&& xa2) noexcept
	{
		if (this == &xa2) return;
		m_iDim1 = xa2.m_iDim1;
		m_iDim2 = xa2.m_iDim2;
		XArray<T>::operator=(static_cast<XArray<T>&&>(xa2));
	}

	//! Performs elementwise addition of the two arrays
	template <class T> inline void XArray2D<T>::operator+=(const XArray2D<T>& xa2)
	{
		if (!IsSameShape(xa2))
			throw std::invalid_argument("invalid_argument 'xa2' in XArray2D<T>::+= (different shape)"); 
		XArray<T>::operator+=(xa2);
	}

	//! Performs elementwise subtraction of the two arrays
	template <class T> inline void XArray2D<T>::operator-=(const XArray2D<T>& xa2)
	{
		if (!IsSameShape(xa2))
			throw std::invalid_argument("invalid_argument 'xa2' in XArray2D<T>::-= (different shape)"); 
		XArray<T>::operator-=(xa2);
	}

	//! Performs elementwise multiplication of the two arrays
	template <class T> inline void XArray2D<T>::operator*=(const XArray2D<T>& xa2)
	{
		if (!IsSameShape(xa2))
			throw std::invalid_argument("invalid_argument 'xa2' in XArray2D<T>::*= (different shape)"); 
		XArray<T>::operator*=(xa2);
	}

	//! Performs elementwise division of the two arrays (checks for division by zero)
	template <class T> inline void XArray2D<T>::operator/=(const XArray2D<T>& xa2)
	{
		if (!IsSameShape(xa2))
			throw std::invalid_argument("invalid_argument 'xa2' in XArray2D<T>::/= (different shape)"); 
		XArray<T>::operator/=(xa2);
	}

	//! Shuffles XArray2D in the way that is required by most DFT routines
	template <class T> void XArray2D<T>::Shuffle()
	{
		if (m_iDim1 != 2 * int(m_iDim1 / 2) || m_iDim2 != 2 * int(m_iDim2 / 2))
			throw std::runtime_error("XArray2D<T>::Shuffle() should only be applied to arrays with even dimensions");

		index_t jj, nyd2 = m_iDim1 / 2, nxd2 = m_iDim2 / 2;
		for (index_t j = 0; j < nyd2; j++)
		{
			jj = j + nyd2;
			for (index_t i = 0; i < nxd2; i++)
				std::swap((*this)[j][i], (*this)[jj][i + nxd2]);
		}
		for (index_t j = nyd2; j < m_iDim1; j++)
		{
			jj = j - nyd2;
			for (index_t i = 0; i < nxd2; i++)
				std::swap((*this)[j][i], (*this)[jj][i + nxd2]);
		}
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
	//	Makes a complex XArray2D object C = A + ib or C = A * exp(ib) from a real XArray2D object A and a scalar b
	//
	/*!
		\brief		Makes a complex XArray2D object C = A + ib or C = A * exp(ib) from a real XArray2D object A and a scalar b
		\param		A	Real XArray2D object representing real part or modulus of a complex object to be constructed
		\param		b	Real scalar value representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray2D object constructed by this function
		\param		bMakePolar if \b true, then C = A * exp(ib) is constructed, else C = A + ib is constructed
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions called from 
					inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray2D object C = A + ib or C = A * exp(ib) from a real XArray2D object and a scalar.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter A
		\par		Example:
\verbatim
XArray2D<float> A(4, 5, 1.0f);
XArray2D<fcomplex> C; 
MakeComplex(A, 0.0f, C, false); //or
C = MakeComplex(A, 0.0f, C, false);
\endverbatim
	*/
	template <class T> void MakeComplex(const XArray2D<T>& A, T b, XArray2D< std::complex<T> >& C, bool bMakePolar)
	{
		MakeComplex((const XArray<T>&)A, b, (XArray< std::complex<T> >&)C, bMakePolar);
		C.SetDims(A.GetDim1(), A.GetDim2());
	}

	template <class T> XArray2D< std::complex<T> > MakeComplex(const XArray2D<T>& A, T b, bool bMakePolar)
	{
		XArray2D< std::complex<T> > C;
		MakeComplex(A, b, C, bMakePolar);
		return C;  // rely on move assignment / constructor to avoid copying this object on return
	}

	//---------------------------------------------------------------------------
	//Function MakeComplex
	//
	//	Makes a complex XArray2D object C = a + iB or C = a * exp(iB) from a real XArray2D object B and a scalar a
	//
	/*!
		\brief		Makes a complex XArray2D object C = a + iB or C = a * exp(iB) from a real XArray2D object B and a scalar a
		\param		a	Real scalar value representing real part or modulus of a complex object to be constructed
		\param		B	Real XArray2D object representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray2D object constructed by this function
		\param		bMakePolar if \b true, then C = a * exp(iB) is constructed, else C = a + iB is constructed
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions called from 
					inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray2D object C = a + iB or C = a * exp(iB) from a real XArray2D object and a scalar.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter B
		\par		Example:
\verbatim
XArray2D<float> B(4, 5, 1.0f);
XArray2D<fcomplex> C;
MakeComplex(1.0f, B, C, false); // or
C = MakeComplex(1.0f, B, false);
\endverbatim
	*/
	template <class T> void MakeComplex(T a, const XArray2D<T>& B, XArray2D< std::complex<T> >& C, bool bMakePolar)
	{
		MakeComplex(a, (const XArray<T>&)B, (XArray< std::complex<T> >&)C, bMakePolar);
		C.SetDims(B.GetDim1(), B.GetDim2());
	}
	
	template <class T> XArray2D< std::complex<T> > MakeComplex(T a, const XArray2D<T>& B, bool bMakePolar)
	{
		XArray2D< std::complex<T> > C;
		MakeComplex(a, B, C, bMakePolar);
		return C;  // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function MakeComplex
	//
	//	Makes a complex XArray2D object C = A + iB or C = A * exp(iB) from 2 real XArray2D objects, A and B
	//
	/*!
		\brief		Makes a complex XArray2D object C = A + iB or C = A * exp(iB) from 2 real XArray2D objects, A and B
		\param		A	Real XArray2D object representing real part or modulus of a complex object to be constructed
		\param		B	Real XArray2D object representing imaginary part or phase of a complex object to be constructed
		\param		C	Complex XArray2D object constructed by this function
		\param		bMakePolar if \b true, then C = A * exp(iB) is constructed, else C = A + iB is constructed
		\exception	std::invalid_argument is thrown if A and B have differented sizes or heads
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions called 
					from inside this function
		\return		\a None
		\par		Description:
			This function makes a complex XArray2D object C = A + iB or C = A * exp(iB) from 2 real XArray2D objects.
			The mode of construction (polar or normal) is determined by the last parameter. The head for the constructed object
			is copied from the parameter A
		\par		Example:
\verbatim
XArray2D<float> A(4, 5, 1.0f), B(4, 5, 2.0f);
XArray2D<fcomplex> C;
MakeComplex(A, B, C, false); // or
C = MakeComplex(A, B, false);
\endverbatim
	*/
	template <class T> void MakeComplex(const XArray2D<T>& A, const XArray2D<T>& B, XArray2D< std::complex<T> >& C, bool bMakePolar)
	{
		MakeComplex((const XArray<T>&)A, (const XArray<T>&)B, (XArray< std::complex<T> >&)C, bMakePolar);
		C.SetDims(A.GetDim1(), A.GetDim2()); 
	}

	template <class T> XArray2D< std::complex<T> > MakeComplex(const XArray2D<T>& A, const XArray2D<T>& B, bool bMakePolar)
	{
		XArray2D< std::complex<T> > C;
		MakeComplex(A, B, C, bMakePolar);
		return C;  // rely on move assignment / constructor to avoid copying this object on return
	}


	//---------------------------------------------------------------------------
	//Function Re
	//
	//	Calculates the real part of a complex XArray2D object
	//
	/*!
		\brief		Calculates the real part of a complex XArray2D object
		\param		C	Input complex XArray2D object
		\param		A	Output real XArray2D object to be made equal to the real part of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the real part of every element of a complex XArray2D object and 
			copies the result into a real XArray2D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray2D<dcomplex> C(4, 5, dcomplex(1.0, 2.0));
XArray2D<double> A;
Re(C, A); // or
A = Re(C);
\endverbatim
	*/
	template <class T> void Re(const XArray2D< std::complex<T> >& C, XArray2D<T>& A)
	{ 
		Re((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2());
	}

	template <class T> XArray2D<T> Re(const XArray2D< std::complex<T> >& C)
	{
		XArray2D<T> A;
		Re(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}

	//---------------------------------------------------------------------------
	//Function Im
	//
	//	Calculates the imaginary part of a complex XArray2D object
	//
	/*!
		\brief		Calculates the imaginary part of a complex XArray2D object
		\param		C	Input complex XArray2D object
		\param		A	Output real XArray2D object to be made equal to the imaginary part of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the imaginary part of every element of a complex XArray2D object
			and copies the result into a real XArray2D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray2D<dcomplex> C(4, 5, dcomplex(1.0, 2.0));
XArray2D<double> A;
Im(C, A); // or
A = Im(C);
\endverbatim
	*/
	template <class T> void Im(const XArray2D< std::complex<T> >& C, XArray2D<T>& A)
	{ 
		Im((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2());
	}

	template <class T> XArray2D<T> Im(const XArray2D< std::complex<T> >& C)
	{
		XArray2D<T> A;
		Im(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}

	//---------------------------------------------------------------------------
	//Function Abs
	//
	//	Calculates the modulus of a complex XArray2D object
	//
	/*!
		\brief		Calculates the modulus of a complex XArray2D object
		\param		C	Input complex XArray2D object
		\param		A	Output real XArray2D object to be made equal to the modulus of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the modulus of every element of a complex XArray2D object
			and copies the result into a real XArray2D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray2D<dcomplex> C(4, 5, dcomplex(1.0, 2.0));
XArray2D<double> A;
Abs(C, A); // or
A = Abs(C);
\endverbatim
	*/	
	template <class T> void Abs(const XArray2D< std::complex<T> >& C, XArray2D<T>& A)
	{ 
		Abs((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2());
	}
	
	template <class T> XArray2D<T> Abs(const XArray2D< std::complex<T> >& C)
	{
		XArray2D<T> A;
		Abs(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}

	
	//---------------------------------------------------------------------------
	//Function Arg
	//
	//	Calculates the argument of a complex XArray2D object
	//
	/*!
		\brief		Calculates the argument of a complex XArray2D object
		\param		C	Input complex XArray2D object
		\param		A	Output real XArray2D object to be made equal to the argument of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the argument of every element of a complex XArray2D object and copies
			the result into a real XArray2D object. The head for the output object A is copied from the input object C
		\par		Example:
\verbatim
XArray2D<dcomplex> C(4, 5, dcomplex(1.0, 2.0));
XArray2D<double> A;
Arg(C, A); // or
A = Arg(C);
\endverbatim
	*/	
	template <class T> void Arg(const XArray2D< std::complex<T> >& C, XArray2D<T>& A)
	{ 
		Arg((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2());
	}

	template <class T> XArray2D<T> Arg(const XArray2D< std::complex<T> >& C)
	{
		XArray2D<T> A;
		Arg(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}
	
	
	//---------------------------------------------------------------------------
	//Function Abs2
	//
	//	Calculates the modulus of a complex XArray2D object
	//
	/*!
		\brief		Calculates the square modulus of a complex XArray2D object
		\param		C	Input complex XArray2D object
		\param		A	Output real XArray2D object to be made equal to the square modulus of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the square modulus of every element of a complex XArray2D object
			and copies the result into a real XArray2D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray2D<dcomplex> C(4, 5, dcomplex(1.0, 2.0));
XArray2D<double> A;
Abs2(C, A); // or
A = Abs2(C);
\endverbatim
	*/	
	template <class T> void Abs2(const XArray2D< std::complex<T> >& C, XArray2D<T>& A)
	{ 
		Abs2((const XArray< std::complex<T> >&)C, (XArray<T>&)A);
		A.SetDims(C.GetDim1(), C.GetDim2());
	}
	
	template <class T> XArray2D<T> Abs2(const XArray2D< std::complex<T> >& C)
	{
		XArray2D<T> A;
		Abs2(C, A);
		return A; // rely on move assignment / constructor to avoid copying this object on return
	}

	
	//---------------------------------------------------------------------------
	//Function CArg
	//
	//	Calculates the 1D-continuous phase of a complex XArray2D object
	//
	/*!
		\brief		Calculates the 1D-continuous phase of a complex XArray2D object
		\param		C	Input complex XArray2D object
		\param		A	Output real XArray2D object to be made equal to the phase of C
		\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
					called from inside this function
		\return		\a None
		\par		Description:
			This function calculates the argument of every element of a complex XArray2D object using the 
			value of the argument of the preceding element to calculate the appropriate 2pi-multiple shift,
			and copies the result into a real XArray2D object. The head for the output object A is copied
			from the input object C
		\par		Example:
\verbatim
XArray2D<dcomplex> C(4, 5, dcomplex(1.0, 2.0));
XArray2D<double> A;
CArg(C, A); // or
A = CArg(C);
\endverbatim
	*/	
	template <class T> void CArg(const XArray2D< std::complex<T> >& C, XArray2D<T>& A) 
	{
		//CArg((const XArray< std::complex<T> >&)C, (XArray<T>&)A); OLD CODE
		index_t idim1 = C.GetDim1();
		index_t idim2 = C.GetDim2();
		A.resize(idim1 * idim2);
		A.SetDims(idim1, idim2);
		
		if (C.GetHeadPtr()) C.GetHeadPtr()->Validate();
		A.SetHeadPtr(C.GetHeadPtr() ? C.GetHeadPtr()->Clone() : 0);

		A[0][0] = carg(C[0][0], T(0));
		for (index_t j = 1; j < idim2; j++)
			A[0][j] = carg(C[0][j], A[0][j - 1]);
		for (index_t i = 1; i < idim1; i++)
		{
			if (i % 2) // odd i
			{
				A[i][idim2 - 1] = carg(C[i][idim2 - 1], A[i - 1][idim2 - 1]); // move down at the end
				for (index_t j = idim2 - 2; j >= 0 && j < idim2; j--)
					A[i][j] = carg(C[i][j], A[i][j + 1]); // move from right to left
			}
			else
			{
				A[i][0] = carg(C[i][0], A[i - 1][0]); // move down at the start
				for (index_t j = 1; j < idim2; j++)
					A[i][j] = carg(C[i][j], A[i][j - 1]); // move from left to right
			}
		}
	}

	template <class T> XArray2D<T> CArg(const XArray2D< std::complex<T> >& C)
	{
		XArray2D<T> A;
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
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray2D<char>;
	template class xar::XArray2D<short>;
	template class xar::XArray2D<long>;
	template class xar::XArray2D<float>;
	template class xar::XArray2D<double>;
	template class xar::XArray2D<xar::fcomplex>;
	template class xar::XArray2D<xar::dcomplex>;
#endif

//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XARRAY2D_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
