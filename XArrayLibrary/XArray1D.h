//Header XArray1D.h
//
//
//	HEADER FILE TITLE:
//
//		One-dimensional XArray class
//
/*!
	\file		XArray1D.h
	\brief		One-dimensional XArray class
	\par		Description:
		This is a 1D version of the XArray class which implements the notion of a one-dimensional 
		resizable (dynamic) array with an optional head
*/
#if !defined XARRAY1D_H
#define XARRAY1D_H
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XArray.h"

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
//Class XArray1D<T>
//
//	One-dimensional XArray class
//
/*!
	\brief		One-dimensional XArray class
	\par		Description:
				This class implements the notion of a one-dimensional resizable (dynamic) array with an optional head.
				Most of the class's functionality is simply inherited from XArray<T>.
	\remarks	This class can be instantiated (unlike its parent XArray<T>)
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
	\remarks    Functions that may require additional information about the head (in excess of what is available in the 
				IXAHead interface) are placed in separate classes (e.g. XArray1DMove<T>).
	\warning	The absence of appropriate specializations of the function GetValuetype() in the base class XArray<T>
				will prevent instantiation of XArray1D<T> object for new types T
*/
	template <class T> class XArray1D : public XArray<T>
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Default constructor
		XArray1D(void) {}
		//! Constructor with a predefined size
		explicit XArray1D(index_t NumPoints,  T tVal = T()) : XArray<T>(NumPoints, tVal) {} 
		//! Promotion from vector
		explicit XArray1D(const vector<T>& rvector) : XArray<T>(rvector) {}
		//! Construction from a raw memory buffer
		XArray1D(T* ptBufBegin, T* ptBufEnd) : XArray<T>(ptBufBegin, ptBufEnd) {} 
		//! Copy constructor
		XArray1D(const XArray1D<T>& rXArray1D) : XArray<T>(rXArray1D) {}
		//! Move constructor
		XArray1D(XArray1D<T>&& rXArray1D) : XArray<T>(std::move(rXArray1D)) {}
		//! Destructor (head is deleted in the base class)
		~XArray1D(void) {}

	// Operators (most member operators are simply inherited from XArray<T>)
	public:
		//! Makes this array a (deep) copy of the rXArray1D
		void operator=(const XArray1D<T>& rXArray1D) { XArray<T>::operator=(rXArray1D); }
		//! Move assignment operator from another XArray1D
		void operator=(XArray1D<T>&& rXArray1D) noexcept { XArray<T>::operator=(std::move(rXArray1D)); }
		//! Makes this array a (deep) copy of the rvector and sets the head pointer to 0
		void operator=(const vector<T>& rvector) { XArray<T>::operator=(rvector); }
		//! Move assignment operator from another vector<T>
		void operator=(vector<T>&& rvector) { XArray<T>::operator=(std::move(rvector)); }

	// Attributes
	public:
		//! Returns the dimensionality of the XArray1D<T> class (it is always equal to eDim1, i.e. one) 
		static _eArrayDim GetArraydim(void) { return eDim1; }

	// Operations
	public:
		//! Accepts an external memory buffer with its contents and makes it the internal 1D array (head does not change)
		//void AcceptMemBuffer(T* ptBufBegin, index_t BufSize);
		//! Relinquishes the responsibility for the memory area occupied by the internal 1D array  (head is deleted)
		//void ReleaseMemBuffer();
		//! Changes the size of the array and fills the NEW elements with a given value (head is not affected)
		//  Call the non-member function ResizeH if the head should be resized as well
		void Resize(index_t NumPoints, T tVal = T());
		//! Changes the size of the array (head is not affected)
		void ResizeA(const std::vector<index_t>& rvecNewSizes);
		//! Resizes the array to zero, FREES THE MEMORY and deletes m_pHead
		void Truncate(void);
		//! Swaps XArrays and their heads
		void Swap(XArray<T>& rXArray); 
		//! Extracts a sub-XArray1D into another XArray1D (head is not affected)
		//  Call the non-member function GetSubarrayH if the head should be returned as well
		void GetSubarray(index_t iBegin, index_t iEnd, XArray1D<T>& rDestSubXArray) const; 
		//! Inserts a sub-XArray1D
		void SetSubarray(const XArray1D<T>& rSrcSubXArray, index_t iBegin); 

	// Overridables
	public:
	
	// Implementation
	protected:

	private:
	// Member variables	
	// Member functions

	};

}
//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
namespace xar
{
	// the following function definitions should ideally be placed into a .cpp file and 'export'ed

	//***** XArray1D member operators

	//***** member functions

	//! Accepts an external memory buffer with its contents and makes it the internal 1D array (head does not change)
	//template <class T> void XArray1D<T>::AcceptMemBuffer(T* ptBufBegin, index_t BufSize)
	//{
	//	if (BufSize < 0)
	//		throw std::invalid_argument("invalid_argument 'BufSize' in XArray1D<T>::AcceptMemBuffer (negative dimension)"); 
	//	acceptMemBuffer(ptBufBegin, BufSize);
	//}

	//! Relinquishes the responsibility for the memory area occupied by the internal 1D array (head is deleted)
	//template <class T> void XArray1D<T>::ReleaseMemBuffer() 
	//{ 
	//	XArray<T>::SetHeadPtr(0);
	//	releaseMemBuffer(); 
	//}

	//! Changes the size of the array and fills the NEW elements with a given value (head is not affected).
	//  Call the non-member function ResizeH if the head should be resized as well
	template <class T> void XArray1D<T>::Resize(index_t newSize, T val)
	{ 
		if (newSize < 0)
			throw std::invalid_argument("invalid_argument 'newSize' in XArray1D<T>::Resize (negative dimension)"); 
		(*this).resize(newSize, val);
	}

	//! Changes the size of the array (head is not affected).
	//  Call the non-member function ResizeH if the head should be resized as well
	template <class T> void XArray1D<T>::ResizeA(const std::vector<index_t>& rvecNewSizes)
	{ 
		if (rvecNewSizes.size() != 1)
			throw std::invalid_argument("invalid_argument 'rvecNewSizes' in XArray1D<T>::ResizeA (wrong dimensionality)"); 
		Resize(rvecNewSizes[0]);
	}

	//! Resizes the array to zero, FREES THE MEMORY and deletes m_pHead
	template <class T> void XArray1D<T>::Truncate()
	{
		XArray<T>::SetHeadPtr(0);
		(*this).truncate();
	}

	//! Swaps XArrays and their heads
	template <class T> void XArray1D<T>::Swap(XArray<T>& rXArray)
	{
		swap(rXArray);
		IXAHead* temp =  XArray<T>::GetHeadPtr() ? XArray<T>::GetHeadPtr()->Clone() : 0;
		IXAHead* temp1 = rXArray.XArray<T>::GetHeadPtr() ? rXArray.XArray<T>::GetHeadPtr()->Clone() : 0;
		XArray<T>::SetHeadPtr(temp1); rXArray.XArray<T>::SetHeadPtr(temp);
	}

	//! Extracts a sub-XArray1D into another XArray1D (head is not affected)
	//  Call the non-member function GetSubarrayH if the head should be returned as well
	template <class T> void XArray1D<T>::GetSubarray(index_t iBegin, index_t iEnd, XArray1D<T>& rDestSubXArray) const
	{
		if (iBegin >= iEnd || iEnd > (*this).size())
			throw std::invalid_argument("invalid argument 'iBegin or iEnd' in XArray1D<T>::GetSubarray");
		index_t iSize = iEnd - iBegin;
		if (rDestSubXArray.size() != iSize) rDestSubXArray.resize(iSize);
		for (index_t i = 0; i < iSize; i++) rDestSubXArray[i] = (*this)[iBegin + i];
	}

	//! Inserts a sub-XArray1D
	template <class T> void XArray1D<T>::SetSubarray(const XArray1D<T>& rSrcSubXArray, index_t iBegin)
	{
		if (iBegin + rSrcSubXArray.size() > (*this).size())
			throw std::invalid_argument("invalid argument 'iBegin or rSrcSubXArray' in XArray1D<T>::SetSubarray");
		for (index_t i = 0; i < rSrcSubXArray.size(); i++) (*this)[iBegin + i] = rSrcSubXArray[i];
		// heads are ignored
	}


} // namespace xar closed
//---------------------------------------------------------------------------
//	RELATED IN-LINE NON-MEMBER DEFINITIONS
//
//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//
// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
#if(0)
	template class xar::XArray1D<char>;
	template class xar::XArray1D<short>;
	template class xar::XArray1D<long>;
	template class xar::XArray1D<float>;
	template class xar::XArray1D<double>;
	template class xar::XArray1D<xar::fcomplex>;
	template class xar::XArray1D<xar::dcomplex>;
#endif


//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
#endif	// XARRAY1D_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
