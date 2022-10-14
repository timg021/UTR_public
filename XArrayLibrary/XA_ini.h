//Header XA_ini.h
//
//	HEADER FILE TITLE:
//
//		Common facilities for XArray and related classes
//
/*!
	\file		XA_ini.h
	\brief		Common facilities for XArray and related classes
	\par		Description:
		This header contains the basic auxiliary facilities required by most XArray and related classes.
*/
#if !defined XA_INI_H
#define XA_INI_H
//---------------------------------------------------------------------------
//	COMPILER DIRECTIVES
//
// The constant XA_DLL_EXPORTS can be declared in a project settings
#if defined XA_DLL_EXPORTS
	// This is used for producing a DLL exporting all of the XArray functionality
	#define XA_API __declspec(dllexport) 
#elif defined XAP_DLL_NET_EXPORTS
	// This is used for producing a DLL exporting all of the XArray functionality
	#define XA_API __declspec(dllexport) 
#elif defined XAP_DLL_IMPORT
	// This is used in projects importing all of the XArray functionality via interfaces from a DLL 
	#define XA_API __declspec(dllimport)
#else
	// This is used in projects that neither export nor import the complete XArray functionality to/from a DLL
	#define XA_API 
#endif

// The following #define is required to resolve the conflict between the min and max macros
// defined in <windef.h> and the corresponding functions defined in STL
#ifndef NOMINMAX
	#define NOMINMAX
#endif
// The above solution for min and max macros does not seem to work,
// so we still have to undefine the macros
#undef min
#undef max

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
// Headers universally used throughout the XAR class library
//
#include <memory>		// for std::unique_ptr, etc.
#include <string>		// for std::string
#include <complex>		// for std::complex
#include <vector>		// for std::vector
#include <stdexcept>	// for std exceptions
#include <assert.h>		// assert macro
//---------------------------------------------------------------------------
//	MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	TYPEDEFS
//
//! This is the only namespace of the XArray class library. Some types, e.g. interfaces, are defined outside the namespace.
namespace xar
{
	//! This type will usually be used in place of size_t
	//	We want the 'size type' in XArray-related classes to be in sync
	//	with the 'size type' of STL
	typedef std::string::size_type index_t;
	//! Type std::complex<float> is meant to be always available
	typedef std::complex<float> fcomplex;
	//! Type std::complex<double> is meant to be always available
	typedef std::complex<double> dcomplex;
//---------------------------------------------------------------------------
//	USING DECLARATIONS
//
	// Type std::string is meant to be always available in XArray classes.
	using std::string;
	// Type std::vector is meant to be always available in XArray classes.
	using std::vector; 
//---------------------------------------------------------------------------
//	FORWARD REFERENCES
//
//---------------------------------------------------------------------------
//	CONSTANT DEFINITIONS
//
	// Commonly used constants for XArray classes
	//
/*!
	\brief		Pi
	\par		Description:
		This is the usual Pi constant
*/
	constexpr double PI = 3.141592653589793;
	constexpr float PIf = float(PI);
/*!
	\brief		Half of Pi
	\par		Description:
		This is the usual Pi constant divided by two
*/
	constexpr double PI2 = PI / 2.0; 
/*!
	\brief		Two Pi
	\par		Description:
		This is the usual Pi constant times two
*/
	constexpr double tPI = 2.0 * PI;
/*!
	\brief		Pi divided by 180
	\par		Description:
		This constant is used to translate degrees into radians and back
*/
	constexpr double PI180 = PI / 180.0;
/*!
	\brief		Fixed distance between the corresponding real and complex valuetypes
	\par		Description:
		This constant is used to automatically find the corresponding complex data type, e.g. eXAFComplex,
		among the members of enum _eValueType, given the value of the corresponding real type, e.g. eXAFloat
*/
	constexpr short ValueTypeComplexShift = 2;
//---------------------------------------------------------------------------
//	ENUMERATED DATA TYPES
//
	//
	// Common enums
	//
	//! Dimensionality values allowed for XArray template classes
	enum _eArrayDim 
	{ 
		eADUnknown, eDim1, eDim2, eDim3, eADMax 
	};
	//! Value(element) types allowed for XArray template classes
	enum _eValueType
	{ 
		eVTUnknown, eXAChar, eXAShort, eXALong, eXAFloat, eXADouble, eXAFComplex, eXADComplex, eVTMax
	};
/*!
	\brief		Literal string array with names of value(element) types
	\par		Description:
		This constant array contains literal strings with names of value(element) types allowed for XArray template classes.
		It is designed to be used and has to be kept in sync with the enum _eValueType.
*/
	// This is, of course, not a enum, but a constant, but it depends on enum _eValueType, and eVTMax in particular
	const char* const chValueType[eVTMax] = 
	{ 
		"vtunknown", "xachar", "xashort", "xalong", "xafloat", "xadouble", "xafcomplex", "xadcomplex"
	};
	//! File types for which I/O operations have been implemented for XArray template classes
	enum _eFileType
	{ 
		eFTUnknown, eHDF5, eVODBIN, eVODASC, eGRDBIN, eGRDASC, eGRCBIN, eGRCASC, eTIFF8, eTIFF16, eTIFF32, eDICOM8, eDICOM16, 
		eJPEG, eJP2_8, eJP2_16, ePNG8, ePNG16, eBMP, eDAT, eFTMax
	};
	//! Various norms (metrics) of XArrays
	enum _eNormID
	{
		eNormTypeMin = -202, eNormTypeMax = -201, eNormIndexOfMin = -102, eNormIndexOfMax = -101,
		eNormYMin = -22, eNormYMax = -21, eNormXMin = -12, eNormXMax = -11, eNormMin = -2,
		eNormMax = -1, eNormC0 = 0, eNormL1 = 1, eNormL2 = 2, eNormAver = 3, eNormStdDev = 4,
		eNormL1N = 6, eNormL2N = 5, eNormEntropy = 7
	};
//---------------------------------------------------------------------------
//	STRUCTURE DEFINITIONS
//
	struct Pair2 { double v; int n; };

	int Pair2comp(const void* pP1, const void* pP2);

//---------------------------------------------------------------------------
//	IN-LINE FUNCTION DEFINITIONS
//
	//
	// Common helper functions
	//
	// NOTE that the following 2 functions have to be defined, because the MS STL implementation
	// does not seem to have the standard functions max() and min() in the <algorithm> header
	//! Calculates maximum of the two values
	template <class T> const T& max(const T& a, const T& b) { return (a > b) ? a : b; } 
	//! Calculates mimimum of the two values
	template <class T> const T& min(const T& a, const T& b) { return (a < b) ? a : b; }
	//! Calculates a modulo b for double values (non-cyclic for negative vs positive values)
	inline double amodb(double a, double b) { return a - b * (long)(a / b) ; } // note that amod(-3, 10) = -3, while nmodm(-3, 10) = 7
	//! Calculates a modulo b for float values (non-cyclic for negative vs positive values)
	inline float amodb(float a, float b) { return a - b * (long)(a / b) ; } // note that amod(-3, 10) = -3, while nmodm(-3, 10) = 7
	//! Calculates n modulo m for integer values (usually the second argument will be converted into double during the call; cyclic for all negative and positive values)
	inline int nmodm(int n, double m) { return (n - int(floor(double(n) / m) * m)); } // the return value is always inside [0, m) (the period)
	//! Calculates n! = n * (n - 1) * (n - 2) * ... * 3 * 2 * 1
	inline index_t factorial(index_t n) { return (n < 2) ? 1 : n * factorial(n - 1); }
	//! Returns the sign of a number
	template <class T> T sgn(const T& a) { return (a > 0) ? T(1) : ((a < 0) ? T(-1) : T(0)); }
	inline fcomplex sgn(const fcomplex& a) { throw std::invalid_argument("invalid argument 'a' in xar::sgn (fcomplex argument)"); }
	inline dcomplex sgn(const dcomplex& a) { throw std::invalid_argument("invalid argument 'a' in xar::sgn (dcomplex argument)"); }
	//
	// NOTE: it is suggested that std exception strings should be formatted as:
	// exception type name (e.g. invalid_argument) followed by the
	// argument or parameter name in quotes (e.g. 'double dblVal') followed by the
	// name of the function where the exception is being thrown from (e.g. Calculate())
	// optional comment in brackets (e.g. argument cannot be negative)
	// e.g. "invalid_argument 'double dblVal' in Calculate() (argument cannot be negative)"
	//
	// Complex-related typedefs and helper functions
	//
	//! Returns the modulus of the argument
	inline double fabs(fcomplex a) { return double(sqrt(a.real() * a.real() + a.imag() * a.imag())); }
	//! Returns the modulus of the argument
	inline double fabs(dcomplex a) { return sqrt(a.real() * a.real() + a.imag() * a.imag()); }
	//! Returns the modulus of the argument
	inline double fabs(char a) { return ::fabs(double(a)); }
	//! Returns the modulus of the argument
	inline double fabs(short a) { return ::fabs(double(a)); }
	//! Returns the modulus of the argument
	inline double fabs(long a) { return ::fabs(double(a)); }
	//! Returns the modulus of the argument
	inline double fabs(float a) { return ::fabs(a); }
	//! Returns the modulus of the argument
	inline double fabs(double a) { return ::fabs(a); }
	//! Returns the argument within [-inf,+inf] of a complex number calculating the 2*pi*n shift using the 'previous' phase value phaso.
	inline double carg(const dcomplex& z, const double phaso)
	{
		double argz = std::arg(z);
		double phaso2P = 6.283185307179586476925286766559 * floor(phaso * 0.15915494309189533576888376337251 + 0.5); //2Pi*floor((phase+Pi)/2Pi))
		double phasoM2P = phaso - phaso2P;
		//assert(phasoM2P >= -3.1415926535897932384626433832795 && phasoM2P <= 3.1415926535897932384626433832795);
		if (argz > phasoM2P + 3.1415926535897932384626433832795) argz -= 6.283185307179586476925286766559;
		else if (argz < phasoM2P - 3.1415926535897932384626433832795) argz += 6.283185307179586476925286766559;
		return argz + phaso2P;
	}
	//! Returns the argument within [-inf,+inf] of a complex number, calculating the 2*pi*n shift using the 'previous' phase value phaso.
	inline float carg(const fcomplex& z, const float phaso)
	{
		float argz = std::arg(z);
		double phaso2P = 6.283185307179586476925286766559 * floor(phaso * 0.15915494309189533576888376337251 + 0.5); //2Pi*floor((phase+Pi)/2Pi))
		double phasoM2P = phaso - phaso2P;
		//assert(phasoM2P >= -3.1415926535897932384626433832795f && phasoM2P <= 3.1415926535897932384626433832795f);
		if (argz > phasoM2P + 3.1415926535897932384626433832795f) argz -= 6.283185307179586476925286766559f;
		else if (argz < phasoM2P - 3.1415926535897932384626433832795f) argz += 6.283185307179586476925286766559f;
		return float(argz + phaso2P);
	}
//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//
//! Displays a simple message either in a console window or in a message box depending on the type of the project
XA_API void DMessage(const string message_string); // displays simple Message

// A variant of the DMessage function with a simple pre-formatting added 
//! Displays warnings formed from an exception's what() string
XA_API void DWhat(const string what_out);

//XA_API int _matherr(struct _exception* pE); //see xa_ini.cpp

double ExeTimer(bool initialize); //measures elapsed time
string ExeTimer1(); //measures elapsed time

} // namespace xar closed

using xar::index_t; // required for use outside the 'xar' namespace, e.g. in interfaces

#endif	// XA_INI_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
