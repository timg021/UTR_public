//Header IXAHead.h
//
//
//	HEADER FILE TITLE:
//
//		Base interface to XArray head classes
//
/*!
	\file		IXAHead.h
	\brief		Base interface to XArray head classes
	\par		Description:
		This header defines base IXAHead interface, which contains functions that
		must be implemented by every XArray-compatible header class

*/
#if !defined IXAHEAD_H
#define IXAHEAD_H

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_ini.h"
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	TYPEDEFS
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
//	Note that although the following entity is declared as 'struct', it is in fact 
//	a pure abstract class (interface)
//---------------------------------------------------------------------------
//Struct IXAHead
//
//	IXAHead base interface class
//
/*!
	\brief		Base interface to XArray head classes
	\par		Description:
				This class is designed to provide common interface to all derived XA head classes
				and to allow manipulations of objects of derived header classes through generic IXAHead* pointers.
				The idea is to make the heads 'dumb', e.g. they will not have information about the 
				'owner' XArray<T> object (otherwise, the head would have to handle all possible template
				parameters, etc, of the XArray<T> classes). The heads are also only 'loosely' attached to
				the owner, i.e. the owner XArray<T> object has a member pointer of the generic IXAHead* type
				(by default this pointer is zero, which corresponds to XArray object with no head). Therefore,
				the owner knows about the attached head, but this knowledge is generally limited to the
				information accessible through the IXAHead interface ('enforced ignorance'). In particular,
				the XArray<T> object deletes the attached head in the ~XArray<T> destructor.
				This design solution implies that the attached head must be of a polymorphic type, 
				must have virtual functions, etc. The solution also implies certain trade-offs.
				On one hand, it allows XArray<T> to be ignorant about the implementation details of the
				head, and therefore to be used with any head derived from (implementing) IXAHead. On the
				other hand, it means that the static type safety is rather weak, e.g. there is nothing
				to stop one from attaching a Wavehead1D head to XArray2D<T> object. Alternatives could be
				to make XArray<T> classes aware of all possible types of heads that can be attached,
				or create separate set of XArray-derived classes for every new head type, which seems undesirable.
	\remarks	This class is an interface, i.e. it contains only pure virtual functions and no
				data members.
	\remarks	This interface must be implemented by every	XArray-compatible header class
	\remarks	This class, being an interface, is declared as 'struct' for convenience, compatibility
				with C-clients and following the common practice advocated e.g. by Microsoft
*/
struct XA_API IXAHead
{
// Constructors
	//! Virtual destructor (needed for proper clean-up)
	virtual ~IXAHead() {}

// Overridables
	//! Returns head type as an std::string
	virtual std::string GetType() const = 0; 
	//! Returns interface to a new copy of the head (this is used in copying of owner objects)
	virtual IXAHead* Clone() const = 0;
	//! Sets head data using a generic IXAHead pointer
	virtual void SetHead(const IXAHead* pHead) = 0;
	//! Determines if two heads are 'equivalent' (e.g. if two XArray<T> objects with these heads can be subtracted)
	virtual bool Equivalent(const IXAHead* pHead) const = 0;
	//! Checks internal consistency of the head (e.g. after reading the head from a file)
	virtual void Validate() const = 0; 
	//! Returns head data as a formatted string
	virtual std::string GetString() const = 0; 
	//! Sets the head data from an appropriately formatted string
	virtual void SetString(const std::string& strData) = 0; 
};
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
// Class Factory functions for base IXAHead-derived classes
// (currently implemented in IXAHWave.cpp)
//
//! Creates a Wavehead object of the specified type and returns an interface pointer to it
XA_API IXAHead* CreateXAHead(const std::string& strHeadType);


#endif //IXAHEAD_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
