//Header XAHWave.h
//
//
//	HEADER FILE TITLE:
//
//		Declarations for WaveheadND classes
//
/*!
	\file		XAHWave.h
	\brief		Declarations for WaveheadND classes
	\par		Description:
				This header contains  declarations for the Wavehead classes, concrete classes
				implementing IXHAWaveND interfaces
*/
#ifndef XAHWAVE_H
#define XAHWAVE_H

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "IXAHWave.h"

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
namespace xar
{
//---------------------------------------------------------------------------
//Class Wavehead1D
//
//	Concrete class implementing the IXAHWave1D interface
//
/*!
	\brief		Concrete class implementing the IXAHWave1D interface
	\par		Description:
				This class contains physical parameters corresponding to a 1D scalar wave
	\remarks	Default units for physical array bounds, etc., in Waveheads are MICRONS
	\remarks	This class is an implementation, it contains data members and can be instantiated
*/
	class Wavehead1D : public IXAHWave1D
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Default constructor
		Wavehead1D() { m_dblWl = 1.e-4; m_dblXlo = 0; m_dblXhi = 1; }
		//! Copy constructor	
		Wavehead1D(const Wavehead1D& rWHead);		
		//! Main constructor
		Wavehead1D(double dblWl, double dblXlo, double dblXhi);
		//! 'Copy' constructor from Wavehead1D object pointed by IXAHead*
		explicit Wavehead1D(const IXAHead* pHead);
		//! Virtual destructor
		virtual ~Wavehead1D() {}
	
	// Operators
	public:
		//! Makes this Wavehead1D a (deep) copy of the rWHead
		void operator=(const Wavehead1D& rWHead); 
		//! Compares two heads for equivalence (not the same as memberwise equality)
		bool operator==(const Wavehead1D& rWHead) const;
		//! Compares two heads for inequivalence (not the same as memberwise inequality)
		bool operator!=(const Wavehead1D& rWHead) const { return !(*this == rWHead); }

	// Attributes
	public:
		//! Returns head type as an std::string
		virtual std::string GetType() const { return std::string("Wavehead1D"); } 
		//! Returns wavelength
		virtual double GetWl() const { return m_dblWl; }
		//! Returns lower physical boundary
		virtual double GetXlo() const { return m_dblXlo; }
		//! Returns higher physical boundary
		virtual double GetXhi() const { return m_dblXhi; }
		//! Sets wavelength, lower and higher physical boundaries
		virtual void SetData(double dblWl, double dblXlo, double dblXhi);

	// Operations
	public:
		//
		// IXAHead implementation
		//
		//! Returns interface to a new copy of the head (this is used in copying of owner objects)
		virtual IXAHead* Clone() const { return new Wavehead1D(*this); }
		//! Sets head data using a generic IXAHead pointer
		virtual void SetHead(const IXAHead* pHead);
		//! Determines if two heads are 'equivalent' (not the same as memberwise equality)
		virtual bool Equivalent(const IXAHead* pHead) const;
		//! Checks internal consistency of the head (e.g. after reading the head from a file)
		virtual void Validate() const { Validate(m_dblWl, m_dblXlo, m_dblXhi); }
		//! Returns head data as a formatted string
		virtual string GetString() const;
		//! Sets the head data from an appropriately formatted string
		virtual void SetString(const string& strData);
		//
		// IXAHWave1D implementation (apart from attributes)
		//
		//! Returns step in physical units
		virtual double GetStep(index_t NumPoints) const;
		//! Resizes the head in accordance with the owner XArray1D
		virtual void Resize(index_t OldNumPoints, index_t NewNumPoints);
		//! Swaps the lower and upper physical boundaries
		virtual void Flip(void) {}
		//! Trims the head in accordance with the owner XArray1D
		virtual void Trim(index_t OldNumPoints, index_t lngLeft, index_t lngRight);
		//! Pads the head in accordance with the owner XArray1D
		virtual void Pad(index_t OldNumPoints, index_t lngLeft, index_t lngRight);
		//! Moves the head (separately from the owner XArray1D)
		virtual void Move(index_t NumPoints, long lngMovePoints);

	// Implementation
	protected:
		//! Checks suitability of the data for Wavehead1D
		void Validate(double dblWl, double dblXlo, double dblXhi) const;
		
	private:
	// Member variables	
		//! Wavelength (usually in microns)
		double m_dblWl; 
		//! Minimum of x-coordinates (usually in microns)
		double m_dblXlo;
		//! Maximum of x-coordinates (usually in microns)
		double m_dblXhi;
		//! Percentage inaccuracy allowed in comparing waveheads for equivalence
		static const double m_dblDelta;  
	
	// Member functions

	};

//---------------------------------------------------------------------------
//Class Wavehead2D
//
//	Concrete class implementing the IXAHWave2D interface
//
/*!
	\brief		Concrete class implementing the IXAHWave2D interface
	\par		Description:
				This class contains physical parameters corresponding to a 2D scalar wave
	\remarks	Default units for physical array bounds, etc., in Waveheads are MICRONS
	\remarks	This class is an implementation, it contains data members and can be instantiated
*/
	class Wavehead2D : public IXAHWave2D
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Default constructor
		Wavehead2D() { m_dblWl = 1.e-4; m_dblYlo = 0; m_dblYhi = 1; m_dblXlo = 0; m_dblXhi = 1; }
		//! Copy constructor
		Wavehead2D(const Wavehead2D& rWHead);			
		//! Main constructor
		Wavehead2D(double dblWl, double dblYlo, double dblYhi, double dblXlo, double dblXhi); 
		//! 'Copy' constructor from a Wavehead2D object pointed to by pHead
		explicit Wavehead2D(const IXAHead* pHead);
		//! Virtual destructor
		virtual ~Wavehead2D() {}

	// Operators
	public:
		//! Makes this Wavehead2D a (deep) copy of the rWHead
		void operator=(const Wavehead2D& rWHead);
		//! Compares two heads for equivalence (not the same as memberwise equality)
		bool operator==(const Wavehead2D& rWHead) const;
		//! Compares two heads for inequivalence (not the same as memberwise inequality)
		bool operator!=(const Wavehead2D& rWHead) const { return !(*this == rWHead); }

	// Attributes
	public:
		//! Returns head type as an std::string
		virtual std::string GetType() const { return std::string("Wavehead2D"); } 
		//! Returns wavelength
		virtual double GetWl() const { return m_dblWl; }
		//! Returns lower physical y-boundary
		virtual double GetYlo() const { return m_dblYlo; }
		//! Returns higher physical y-boundary
		virtual double GetYhi() const { return m_dblYhi; }
		//! Returns lower physical x-boundary
		virtual double GetXlo() const { return m_dblXlo; }
		//! Returns higher physical x-boundary
		virtual double GetXhi() const { return m_dblXhi; }
		//! Sets wavelength, lower and upper x- and y-boundaries
		virtual void SetData(double dblWl, double dblYlo, double dblYhi, double dblXlo, double dblXhi);

	// Operations
	public:
		//
		// IXAHead implementation
		//
		//! Returns interface to a new copy of the head (this is used in copying of owner objects)
		virtual IXAHead* Clone() const { return new Wavehead2D(*this); }
		//! Sets head data using a generic IXAHead pointer
		virtual void SetHead(const IXAHead* pHead);
		//! Determines if two heads are 'equivalent' (not the same as memberwise equality)
		virtual bool Equivalent(const IXAHead* pHead) const;
		//! Checks internal consistency of the head (e.g. after reading the head from a file)
		virtual void Validate() const { Validate(m_dblWl, m_dblYlo, m_dblYhi, m_dblXlo, m_dblXhi); }
		//! Returns head data as a formatted string
		virtual string GetString() const;
		//! Sets the head data from an appropriately formatted string
		virtual void SetString(const string& strData);
		//
		// IXAHWave2D implementation (apart from attributes)
		//
		//! Returns y-step in physical units
		virtual double GetYStep(index_t NumYPoints) const;
		//! Returns x-step in physical units
		virtual double GetXStep(index_t NumXPoints) const;
		//! Resizes the head in accordance with the owner XArray2D
		virtual void Resize(index_t OldNumYPoints, index_t OldNumXPoints, index_t NewNumYPoints, index_t NewNumXPoints);
		//! Transposes the head (swaps X and Y dimensions) in accordance with the owner XArray2D
		virtual void Transpose(void);
		//! Trims the head in accordance with the owner XArray2D
		virtual void Trim(index_t OldNumYPoints, index_t OldNumXPoints, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight);
		//! Pads the head in accordance with the owner XArray2D
		virtual void Pad(index_t OldNumYPoints, index_t OldNumXPoints, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight);
		//! Moves the head (separately from the owner XArray2D)
		virtual void Move(index_t NumYPoints, index_t NumXPoints, long lngMoveYPoints, long lngMoveXPoints);

	// Implementation
	protected:
		//! Checks suitability of the data for Wavehead2D
		void Validate(double dblWl, double dblYlo, double dblYhi, double dblXlo, double dblXhi) const;
	
	private:
	// Member variables	
		//! Wavelength (usually in microns)
		double m_dblWl;
		//! Minimum of y-coordinates (usually in microns)
		double m_dblYlo;
		//! Maximum of y-coordinates (usually in microns)
		double m_dblYhi;
		//! Minimum of x-coordinates (usually in microns)
		double m_dblXlo;
		//! Maximum of x-coordinates (usually in microns)
		double m_dblXhi;
		//! Percentage inaccuracy allowed in comparing waveheads for equality
		static const double m_dblDelta;
	
	// Member functions

	};

//---------------------------------------------------------------------------
//Class Wavehead3D
//
//	Concrete class implementing the IXAHWave3D interface
//
/*!
	\brief		Concrete class implementing the IXAHWave3D interface
	\par		Description:
				This class contains physical parameters corresponding to a 3D scalar wave
	\remarks	Default units for physical array bounds, etc., in Waveheads are MICRONS
	\remarks	This class is an implementation, it contains data members and can be instantiated
*/
	class Wavehead3D : public IXAHWave3D
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Default constructor
		Wavehead3D() { m_dblWl = 1.e-4; m_dblZlo = 0; m_dblZhi = 1; m_dblYlo = 0; m_dblYhi = 1; m_dblXlo = 0; m_dblXhi = 1; }
		//! Copy constructor
		Wavehead3D(const Wavehead3D& rWHead);
		//! Main constructor
		Wavehead3D(double dblWl, double dblZlo, double dblZhi, double dblYlo, double dblYhi, double dblXlo, double dblXhi);
		//! 'Copy' constructor from a Wavehead3D object pointed to by pHead
		explicit Wavehead3D(const IXAHead* pHead);
		//! Virtual destructor
		virtual ~Wavehead3D() {}

	// Operators
	public:
		//! Makes this Wavehead3D a (deep) copy of the rWHead
		void operator=(const Wavehead3D& rWHead);
		//! Compares two heads for equivalence (not the same as memberwise equality)
		bool operator==(const Wavehead3D& rWHead) const;
		//! Compares two heads for inequivalence (not the same as memberwise inequality)
		bool operator!=(const Wavehead3D& rWHead) const { return !(*this == rWHead); }

	// Attributes
	public:
		//! Returns head type as an std::string
		virtual std::string GetType() const { return std::string("Wavehead3D"); } 
		//! Returns wavelength
		virtual double GetWl() const { return m_dblWl; }
		//! Returns lower physical z-boundary
		virtual double GetZlo() const { return m_dblZlo; }
		//! Returns higher physical z-boundary
		virtual double GetZhi() const { return m_dblZhi; }
		//! Returns lower physical y-boundary
		virtual double GetYlo() const { return m_dblYlo; }
		//! Returns higher physical y-boundary
		virtual double GetYhi() const { return m_dblYhi; }
		//! Returns lower physical x-boundary
		virtual double GetXlo() const { return m_dblXlo; }
		//! Returns higher physical x-boundary
		virtual double GetXhi() const { return m_dblXhi; }
		//! Sets wavelength, lower and upper x-, y- and z-boundaries
		virtual void SetData(double dblWl, double dblZlo, double dblZhi, double dblYlo, double dblYhi, double dblXlo, double dblXhi);

		// Operations
	public:
		//
		// IXAHead implementation
		//
		//! Returns interface to a new copy of the head (this is used in copying of owner objects)
		virtual IXAHead* Clone() const { return new Wavehead3D(*this); }
		//! Sets head data using a generic IXAHead pointer
		virtual void SetHead(const IXAHead* pHead);
		//! Determines if two heads are 'equivalent' (not the same as memberwise equality)
		virtual bool Equivalent(const IXAHead* pHead) const;
		//! Checks internal consistency of the head (e.g. after reading the head from a file)
		virtual void Validate() const { Validate(m_dblWl, m_dblZlo, m_dblZhi, m_dblYlo, m_dblYhi, m_dblXlo, m_dblXhi); }
		//! Returns head data as a formatted string
		virtual string GetString() const;
		//! Sets the head data from an appropriately formatted string
		virtual void SetString(const string& strData);
		//
		// IXAHWave3D implementation (apart from attributes)
		//
		//! Returns z-step in physical units
		virtual double GetZStep(index_t NumYPoints) const;
		//! Returns y-step in physical units
		virtual double GetYStep(index_t NumYPoints) const;
		//! Returns x-step in physical units
		virtual double GetXStep(index_t NumXPoints) const;
		//! Resizes the head in accordance with the owner XArray3D
		virtual void Resize(index_t OldNumZPoints, index_t OldNumYPoints, index_t OldNumXPoints, index_t NewNumZPoints, index_t NewNumYPoints, index_t NewNumXPoints);
		//! Trims the head in accordance with the owner XArray3D
		virtual void Trim(index_t OldNumZPoints, index_t OldNumYPoints, index_t OldNumXPoints, index_t lngZLeft, index_t lngZRight, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight);
		//! Pads the head in accordance with the owner XArray3D
		virtual void Pad(index_t OldNumZPoints, index_t OldNumYPoints, index_t OldNumXPoints, index_t lngZLeft, index_t lngZRight, index_t lngYLeft, index_t lngYRight, index_t lngXLeft, index_t lngXRight);
		//! Moves the head (separately from the owner XArray3D)
		virtual void Move(index_t NumZPoints, index_t NumYPoints, index_t NumXPoints, long lngMoveZPoints, long lngMoveYPoints, long lngMoveXPoints);

	// Implementation
	protected:
		//! Checks suitability of the data for Wavehead3D
		void Validate(double dblWl, double dblZlo, double dblZhi, double dblYlo, double dblYhi, double dblXlo, double dblXhi) const;
	
	private:
	// Member variables	
		//! Wavelength (usually in microns)
		double m_dblWl;
		//! Minimum of z-coordinates (usually in microns)
		double m_dblZlo;
		//! Maximum of z-coordinates (usually in microns)
		double m_dblZhi;
		//! Minimum of y-coordinates (usually in microns)
		double m_dblYlo;
		//! Maximum of y-coordinates (usually in microns)
		double m_dblYhi;
		//! Minimum of x-coordinates (usually in microns)
		double m_dblXlo;
		//! Maximum of x-coordinates (usually in microns)
		double m_dblXhi;
		//! Percentage inaccuracy allowed in comparing waveheads for equality
		static const double m_dblDelta;  
	};

} // namespace xar closed
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//

#endif // XAHWAVE_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//
