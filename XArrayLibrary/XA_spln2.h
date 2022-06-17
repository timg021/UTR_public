//Header XA_spln2.h
//
//
//	HEADER FILE TITLE:
//
//		Two-dimensional spline interpolations
//
//
/*!
	\file		XA_spln2.h
	\brief		Two dimensional spline interpolations and related operations
	\par		Description:
		This class implements simple spline interpolations of XArray2D<T> objects and related operations.
*/
#pragma once
//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "XA_head2.h"


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
//Class XArray2DSpln
//
//	Two-dimensional spline interpolations
//
/*!
	\brief		Two-dimensional spline interpolations and related operations
	\par		Description:
				This class template defines a 'wrapper' around the XArray2D<T> object
				on which it operates; it contains functions implementing simple
				spline interpolations of the 'wrapped' XArray2D<T> object.
	\remarks	Unlike XArray2DMove<T>, most functions of this class place result in a separate 2D array.
	\remarks	This class is not designed for polymorphism and does not contain virtual functions
*/

	template <class T> class XArray2DSpln
	{
	// Enumerators
	// Structures
	// Constructors
	public:
		//! Constructor
		XArray2DSpln(const XArray2D<T>& rXAr2D);
	protected:
		//! Copy constructor (declared protected to prohibit copying)
		XArray2DSpln(const XArray2DSpln<T>& rCopy) : m_rXArray2D(rCopy.m_rXArray2D) { Initialize(); }
	public:
		//! Destructor 
		~XArray2DSpln() {}
	
	// Operators
	protected:
		//! Assignment (declared protected to prohibit copying)
		void operator=(const XArray2DSpln<T>& rCopy);

	// Attributes
	public:
		// NOTE: the absence of appropriate specializations of the following function
		// will prevent instantiation of XArray2DSpln<T> objects for types T other than float or double
		//! Returns the xar::_eValueType corresponding to T
		static _eValueType GetValuetype(void);
		//! Returns a reference to the non-modifiable 'wrapped' XArray2D<T> object
		const XArray2D<T>& GetBaseObject() const { return m_rXArray2D; }
		//! Returns a reference to the 'wrapped' XArray2D<T> object
		const XArray2D<T>& GetBaseObject() { return m_rXArray2D; }

	// Operations
	public:
		//! Produces an interpolated value at point (dblX, dblY)
		T Lin2_spln(double dblY, double dblX) const;
		//! Bilinearly interpolates 2D array in accordance with new step sizes
		void Interpolate(XArray2D<T>& xaResult, double dblYStep, double dblXStep) const;
		//! Bilinearly interpolates 2D array in accordance with new array dimensions
		void Interpolate(XArray2D<T>& xaResult, index_t iNewDim1, index_t iNewDim2) const;
		//! Rotates 2D array clockwise around a given point by a given angle 
		void Rotate(XArray2D<T>& xaResult, double dblTheta, double dblYCentre, double dblXCentre, T Backgr) const;
		//! Rotates 2D array clockwise around a given point by a given angle and returns the resultant rotated XArray2D<T> object
		XArray2D<T> Rotate(double dblTheta, double dblYCentre, double dblXCentre, T Backgr) const;
		//! Creates arrays of grid point positions(coordinates) for image rotation
		void RotPositions(XArray2D<float>& xaS, XArray2D<float>& xaT, XArray2D<short>& xaJ0, XArray2D<short>& xaI0, double dblTheta, double dblYCentre, double dblXCentre) const;
		//! Rotates 2D array using existing 2D arrays of rotated grid positions
		void Rotate0(XArray2D<T>& xaResult, const XArray2D<float>& xaS, const XArray2D<float>& xaT, const XArray2D<short>& xaJ0, const XArray2D<short>& xaI0, T Backgr) const;
		//! Creates index mapping for image rotation
		void RotIndexes(XArray2D<short>& xaJ0, XArray2D<short>& xaI0, double dblTheta, double dblYCentre, double dblXCentre) const;
		//! Rotates 2D array using pre-calculated index mapping
		void Rotate1(XArray2D<T>& xaResult, const XArray2D<short>& xaJ0, const XArray2D<short>& xaI0, T Backgr) const;
		//! (De)Magnifies 2D array by factor M with respect to the centre of dilation
		void Magnify(XArray2D<T>& xaResult, double dblM, double dblYCentre, double dblXCentre, T Backgr) const;
		//!	Tilts 2D array by the specified angle around the tilt axis
		void Tilt(XArray2D<T>& xaResult, double dblYAxis, double dblXAxis, double dblAngleYAxis, double dblTiltAngle, T Backgr) const;

	// Overridables

	// Implementation
	private:
	// Member variables	
		//! Reference to the 'wrapped' XArray2D<T> object that contains array of data at equidistant x-points 
		const XArray2D<T>& m_rXArray2D;
		//! Array dimensions
		index_t m_ny, m_nx;
		//! Wavehead2D data or its default equivalent
		double m_wl, m_ylo, m_yhi, m_xlo, m_xhi;
		//! Auxilliary values
		index_t m_ny1, m_nx1;  // array dimensions minus 1
		double m_xst, m_yst, m_ayst, m_axst; // array steps and their inverses
		static const double m_a05; // = std::numeric_limits<T>::is_integer ? 0.5 : 0.0;

	// Helper functions
		void Initialize();
	};

}  //namespace xar closed


//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//! Returns the xar::_eValueType corresponding to T=char
template<> inline xar::_eValueType xar::XArray2DSpln<char>::GetValuetype() { return eXAChar; }
//! Returns the xar::_eValueType corresponding to T=short
template<> inline xar::_eValueType xar::XArray2DSpln<short>::GetValuetype() { return eXAShort; }
//! Returns the xar::_eValueType corresponding to T=long
template<> inline xar::_eValueType xar::XArray2DSpln<long>::GetValuetype() { return eXALong; }
//! Returns the xar::_eValueType corresponding to T=float
template<> inline xar::_eValueType xar::XArray2DSpln<float>::GetValuetype() { return eXAFloat; }
//! Returns the xar::_eValueType corresponding to T=double
template<> inline xar::_eValueType xar::XArray2DSpln<double>::GetValuetype() { return eXADouble; }

// NOTE!!! Complex-valued arrays would require specialized templates for some functions of this class

template<> const double xar::XArray2DSpln<char>::m_a05 = 0.5;
template<> const double xar::XArray2DSpln<short>::m_a05 = 0.5;
template<> const double xar::XArray2DSpln<long>::m_a05 = 0.5;
template<> const double xar::XArray2DSpln<float>::m_a05 = 0.0;
template<> const double xar::XArray2DSpln<double>::m_a05 = 0.0;


//! Constructor
template <class T> xar::XArray2DSpln<T>::XArray2DSpln(const XArray2D<T>& rXAr2D) : m_rXArray2D(rXAr2D)
{
	GetValuetype();
	Initialize();
}


//! Initialize member variables (has to be called again after m_rXArray2D changes)
template <class T> void xar::XArray2DSpln<T>::Initialize()
{
	m_ny = m_rXArray2D.GetDim1();
	m_nx = m_rXArray2D.GetDim2();
	m_ny1 = m_ny - 1;
	m_nx1 = m_nx - 1;
	
	const IXAHWave2D* pHead = GetIXAHWave2D(m_rXArray2D);
	if (pHead)
	{
		m_wl = pHead->GetWl();
		m_ylo = pHead->GetYlo();
		m_yhi = pHead->GetYhi();
		m_xlo = pHead->GetXlo();
		m_xhi = pHead->GetXhi();
	}
	else
	{
		m_wl = 0.0001;
		m_ylo = 0;
		m_yhi = (double)m_ny;
		m_xlo = 0;
		m_xhi = (double)m_nx;
	}
	m_yst = GetYStep(m_rXArray2D);
	m_xst = GetXStep(m_rXArray2D);
	m_ayst = 1.0 / m_yst;
	m_axst = 1.0 / m_xst;
}


//! Assignment  (declared protected to prohibit copying)
template <class T> void xar::XArray2DSpln<T>::operator=(const XArray2DSpln<T>& rCopy)
{ 
	if (this == &rCopy) 
		return; 
	else
		throw std::logic_error("logic_error in XArray2DSpln<T>::operator= (array reference cannot be changed after construction)");
}

//---------------------------------------------------------------------------
//Function XArray2DSpln<T>::Lin2_spln
//
//	Produces an interpolated value at point (dblX, dblY)
//
/*!
	\brief		Produces an interpolated value at point (dblX, dblY)
	\param		dblY Y-coordinate of the interpolation point
	\param		dblX X-coordinate of the interpolation point
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a The interpolated value
	\par		Description:
		This function produces a value at the point (x,y), xlo<=x<=xhi, ylo<=y<=yhi, 
		by means of bilinear spline
		
*/	
template <class T> T xar::XArray2DSpln<T>::Lin2_spln(double dblY, double dblX) const
{
	double s = (dblY - m_ylo) * m_ayst;
	double t = (dblX - m_xlo) * m_axst;
	long j0 = long(s);
	long i0 = long(t);

	if (j0 >= long(m_ny1))
		if (dblY <= m_yhi) 
			j0--; //case dblY == yhi
		else 
			throw std::invalid_argument("invalid_argument 'dblY' in XArray2DSpln<T>::Lin2_spln ( dblY > yhi)");
	else 
		if (j0 < 0) 
			std::invalid_argument("invalid_argument 'dblY' in XArray2DSpln<T>::Lin2_spln ( dblY < ylo)");

	if (i0 >= long(m_nx1))
		if (dblX <= m_xhi) 
			i0--;  //case dblX == xhi
		else 
			throw std::invalid_argument("invalid_argument 'dblX' in XArray2DSpln<T>::Lin2_spln ( dblX > xhi)");
	else 
		if (i0 < 0) 
			throw std::invalid_argument("invalid_argument 'dblX' in XArray2DSpln<T>::Lin2_spln ( dblX < xlo)");

	s -= j0;
	t -= i0;

	return T((1.0 - t) * ((1.0 - s) * m_rXArray2D[j0][i0] + s * m_rXArray2D[j0+1][i0]) + t * ((1.0 - s) * m_rXArray2D[j0][i0+1] + s * m_rXArray2D[j0+1][i0+1]));
}


//---------------------------------------------------------------------------
//Function XArray2DSpln<T>::Interpolate
//
//	Bilinearly interpolates 2D array in accordance with new step sizes 
//
/*!
	\brief		Bilinearly interpolates 2D array in accordance with new step sizes
	\param		xaResult Resultant interpolated 2D array
	\param		dblY Y-coordinate of the interpolation point
	\param		dblX X-coordinate of the interpolation point
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function interpolates the 2D array according to new step sizes
		by means of bilinear spline
		
*/	
template <class T> void xar::XArray2DSpln<T>::Interpolate(XArray2D<T>& xaResult, double dblYStep, double dblXStep) const
{
	if (dblYStep <= 0 || dblXStep <= 0)
		throw std::invalid_argument("invalid_argument 'dblYStep or dblXStep' in XArray2DSpln<T>::Interpolate ( new step <= 0)");

	index_t ny1 = index_t((m_yhi - m_ylo) / dblYStep + 1);
	index_t nx1 = index_t((m_xhi - m_xlo) / dblXStep + 1);

	xaResult.Resize(ny1, nx1, 0);

	// create new head (note that the new head is created even if the array did not have a had before the interpolation)
	IXAHWave2D* pHead1 = CreateWavehead2D();
	pHead1->SetData(m_wl, m_ylo, m_ylo + dblYStep * (ny1 - 1), m_xlo, m_xlo + dblXStep * (nx1 - 1));
	xaResult.SetHeadPtr(pHead1);

	double dyst = dblYStep * m_ayst;
	double dxst = dblXStep * m_axst;

	long i0, j0;
	double t, t1, s, s1;

	for (index_t j = 0; j < ny1; j++)
	{
		s = j * dyst;
		j0 = long(s);
		if (j0 >= long(m_ny1)) j0 = long(m_ny1 - 1); // otherwise j0+1 may be out of bounds
		s -= j0;
		s1 = 1.0 - s;
		
		for (index_t i = 0; i < nx1; i++)
		{
			t = i * dxst;	
			i0 = long(t);
			if (i0 >= long(m_nx1)) i0 = long(m_nx1 - 1); // otherwise i0+1 may be out of bounds
			t -= i0;
			t1 = 1.0 - t;

			xaResult[j][i] = (T)(t1 * (s1 * m_rXArray2D[j0][i0] + s * m_rXArray2D[j0+1][i0]) + 
				t * (s1 * m_rXArray2D[j0][i0+1] + s * m_rXArray2D[j0+1][i0+1]) + m_a05);
		}
	}
}


//---------------------------------------------------------------------------
//Function XArray2DSpln<T>::Interpolate
//
//	Bilinearly interpolates 2D array in accordance with new array dimensions 
//
/*!
	\brief		Bilinearly interpolates 2D array in accordance with new array dimensions
	\param		xaResult Resultant interpolated 2D array
	\param		dblY Y-coordinate of the interpolation point
	\param		dblX X-coordinate of the interpolation point
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function interpolates the 2D array according to new array dimensions
		by means of bilinear spline
		
*/	
template <class T> void xar::XArray2DSpln<T>::Interpolate(XArray2D<T>& xaResult, index_t iNewDim1, index_t iNewDim2) const
{
	if (iNewDim1 <= 1 || iNewDim2 <= 1)
		throw std::invalid_argument("invalid_argument 'iNewDim1 or iNewDim2' in XArray2DSpln<T>::Interpolate ( new array dimension <= 1)");

	double dblYStep = (m_yhi - m_ylo) / (iNewDim1 - 1);
	double dblXStep = (m_xhi - m_xlo) / (iNewDim2 - 1);
	Interpolate(xaResult, dblYStep, dblXStep);
}

		
		
//---------------------------------------------------------------------------
//Function XArray2DSpln<T>::Rotate
//
//	Rotates 2D array clockwise around a given point by a given angle 
//
/*!
	\brief		Rotates 2D array clockwise around a given point by a given angle
	\param		xaResult Resultant rotated 2D array
	\param		dblTheta Rotation angle in degrees
	\param		dblYCentre Y-coordinate of the centre of rotation
	\param		dblXCentre X-coordinate of the centre of rotation
	\param		Backgr Value for filling "background" areas around the rotated array
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function rotates 2D array clockwise around a specified point by a specified angle
		using bilinear interpolation
		
*/	
template <class T> void xar::XArray2DSpln<T>::Rotate(XArray2D<T>& xaResult, double dblTheta, double dblYCentre, double dblXCentre, T Backgr) const
{
	long j0, i0;
	double y, x, ysint, ycost, s, t, s1, t1;

	dblTheta *= PI180;
	double st = sin(dblTheta);
	double ct = cos(dblTheta);

	double yloc = m_ylo - dblYCentre;
	double xloc = m_xlo - dblXCentre;
	
	if (xaResult.GetDim1() != m_ny || xaResult.GetDim2() != m_nx) xaResult.Resize(m_ny, m_nx);

	// create new head (note that the new head is created even if the array did not have a had before the rotation)
	IXAHWave2D* pHead1 = CreateWavehead2D();
	pHead1->SetData(m_wl, m_ylo, m_yhi, m_xlo, m_xhi);
	xaResult.SetHeadPtr(pHead1);

	for (index_t j = 0; j < m_ny; j++)
	{	
		y = yloc + m_yst * j;
		ysint = y * st - xloc;
		ycost = y * ct - yloc;

		for (index_t i = 0; i < m_nx; i++)
		{
			x = xloc + m_xst * i;
			s = (ycost - x * st) * m_ayst;
			t = (ysint + x * ct) * m_axst;

			j0 = long(floor(s));
			if (j0 < 0) { xaResult[j][i] = Backgr; continue; }
			else if (j0 >= long(m_ny1))
					if (j0 == m_ny1 && double(j0) == s) j0--;
					else { xaResult[j][i] = Backgr; continue; }

			s -= j0;
			s1 = 1.0 - s;

			i0 = long(floor(t));
			if (i0 < 0) { xaResult[j][i] = Backgr; continue; }
			else if (i0 >= long(m_nx1))
					if (i0 == m_nx1 && double(i0) == t) i0--;
					else { xaResult[j][i] = Backgr; continue; }

			t -= i0;
			t1 = 1.0 - t;

			xaResult[j][i] = (T)(t1 * (s1 * m_rXArray2D[j0][i0] + s * m_rXArray2D[j0+1][i0]) + 
				t * (s1 * m_rXArray2D[j0][i0+1] + s * m_rXArray2D[j0+1][i0+1]) + m_a05);
		}
	}
}

template <class T> xar::XArray2D<T> xar::XArray2DSpln<T>::Rotate(double dblTheta, double dblYCentre, double dblXCentre, T Backgr) const
{
	xar::XArray2D<T> xaResult;
	Rotate(xaResult, dblTheta, dblYCentre, dblXCentre, Backgr);
	return xaResult; // relying on move constructor or assignment to eliminate copying of xaResult upon the return
}


//! Creates arrays of grid point positions(coordinates) for image rotation
// The resultant arrays xaS, xaT, xaJ0, xaI0 are used in the function Rotate0
template <class T> void xar::XArray2DSpln<T>::RotPositions(XArray2D<float>& xaS, XArray2D<float>& xaT, XArray2D<short>& xaJ0, XArray2D<short>& xaI0, double dblTheta, double dblYCentre, double dblXCentre) const
{
	short i0, j0;
	float s, t;
	double y, x, ysint, ycost;

	dblTheta *= PI180;
	double st = sin(dblTheta);
	double ct = cos(dblTheta);

	double yloc = m_ylo - dblYCentre;
	double xloc = m_xlo - dblXCentre;
	
	if (xaS.GetDim1() != m_ny || xaS.GetDim2() != m_nx) xaS.Resize(m_ny, m_nx);
	if (xaT.GetDim1() != m_ny || xaT.GetDim2() != m_nx) xaT.Resize(m_ny, m_nx);
	if (xaJ0.GetDim1() != m_ny || xaJ0.GetDim2() != m_nx) xaJ0.Resize(m_ny, m_nx);
	if (xaI0.GetDim1() != m_ny || xaI0.GetDim2() != m_nx) xaI0.Resize(m_ny, m_nx);

	for (index_t j = 0; j < m_ny; j++)
	{	
		y = yloc + m_yst * j;
		ysint = y * st - xloc;
		ycost = y * ct - yloc;

		for (index_t i = 0; i < m_nx; i++)
		{
			x = xloc + m_xst * i;
			s = float((ycost - x * st) * m_ayst);
			t = float((ysint + x * ct) * m_axst);

			j0 = short(floor(s));
			if (j0 < 0 || j0 >= long(m_ny1)) { xaJ0[j][i] = -1; xaI0[j][i] = -1; continue; }
			else xaJ0[j][i] = j0;
			xaS[j][i] = s - j0;

			i0 = short(floor(t));
			if (i0 < 0 || i0 >= long(m_nx1)) { xaJ0[j][i] = -1; xaI0[j][i] = -1; continue; }
			else xaI0[j][i] = i0;
			xaT[j][i] = t - i0;
		}
	}
}


//! Rotates 2D array using existing 2D arrays of rotated grid positions
// This function works faster than Rotate when a large set of images needs to be rotated by the same angle
template <class T> void xar::XArray2DSpln<T>::Rotate0(XArray2D<T>& xaResult, const XArray2D<float>& xaS, const XArray2D<float>& xaT, const XArray2D<short>& xaJ0, const XArray2D<short>& xaI0, T Backgr) const
{
	if (xaResult.GetDim1() != m_ny || xaResult.GetDim2() != m_nx) xaResult.Resize(m_ny, m_nx);

	if (xaS.GetDim1() != m_ny || xaS.GetDim2() != m_nx || xaT.GetDim1() != m_ny || xaT.GetDim2() != m_nx ||
		xaJ0.GetDim1() != m_ny || xaJ0.GetDim2() != m_nx || xaI0.GetDim1() != m_ny || xaI0.GetDim2() != m_nx)
		throw std::invalid_argument("invalid_argument in XArray2DSpln<T>::Rotate0 (coordinate arrays have wrong dimensions)"); 

	short i0, j0;
	double s, t, s1, t1;

	for (index_t j = 0; j < m_ny; j++)
	{	
		for (index_t i = 0; i < m_nx; i++)
		{
			j0 = xaJ0[j][i];
			i0 = xaI0[j][i];

			if (j0 >= 0 && i0 >= 0)			
			{
				s = xaS[j][i];
				t = xaT[j][i];
				s1 = 1.0 - s;
				t1 = 1.0 - t;
				xaResult[j][i] = (T)(t1 * (s1 * m_rXArray2D[j0][i0] + s * m_rXArray2D[j0+1][i0]) + 
					t * (s1 * m_rXArray2D[j0][i0+1] + s * m_rXArray2D[j0+1][i0+1]) + m_a05);
			}
			else
				xaResult[j][i] = Backgr;
		}
	}
}


//! Creates index mapping for image rotation
// The resultant integer arrays xaJ0 and xaI0 are used in the function Rotate1
template <class T> void xar::XArray2DSpln<T>::RotIndexes(XArray2D<short>& xaJ0, XArray2D<short>& xaI0, double dblTheta, double dblYCentre, double dblXCentre) const
{
	if (long(m_ny) > std::numeric_limits<short>::max() || long(m_nx) > std::numeric_limits<short>::max())
		throw std::invalid_argument("invalid_argument in XArray2DSpln<T>::RotIndexes (input array too large)");

	short j0, i0;
	double y, x, ysint, ycost, s, t;

	dblTheta *= PI180;
	double st = sin(dblTheta);
	double ct = cos(dblTheta);

	double yloc = m_ylo - dblYCentre;
	double xloc = m_xlo - dblXCentre;
	
	if (xaJ0.GetDim1() != m_ny || xaJ0.GetDim2() != m_nx) xaJ0.Resize(m_ny, m_nx);
	if (xaI0.GetDim1() != m_ny || xaI0.GetDim2() != m_nx) xaI0.Resize(m_ny, m_nx);

	for (short j = 0; j < short(m_ny); j++)
	{	
		y = yloc + m_yst * j;
		ysint = y * st - xloc;
		ycost = y * ct - yloc;

		for (short i = 0; i < short(m_nx); i++)
		{
			x = xloc + m_xst * i;
			s = (ycost - x * st) * m_ayst;
			t = (ysint + x * ct) * m_axst;

			j0 = short(floor(s));
			if (j0 < 0 || j0 > short(m_ny1)) { xaJ0[j][i] = -1; xaI0[j][i] = -1; continue; }
			else xaJ0[j][i] = j0;

			i0 = short(floor(t));
			if (i0 < 0 || i0 > short(m_nx1)) { xaJ0[j][i] = -1; xaI0[j][i] = -1; continue; }
			else xaI0[j][i] = i0;
		}
	}
}


//! Rotation using existing indexes
// This function works much faster than Rotate when a large set of images needs to be rotated by the same angle,
// but it is slightly less accurate, as the bililnear interpolation is not performed
template <class T> void xar::XArray2DSpln<T>::Rotate1(XArray2D<T>& xaResult, const XArray2D<short>& xaJ0, const XArray2D<short>& xaI0, T Backgr) const
{
	if (xaResult.GetDim1() != m_ny || xaResult.GetDim2() != m_nx) 
		xaResult.Resize(m_ny, m_nx);

	if (xaJ0.GetDim1() != m_ny || xaJ0.GetDim2() != m_nx || xaI0.GetDim1() != m_ny || xaI0.GetDim2() != m_nx)
		throw std::invalid_argument("invalid_argument in XArray2DSpln<T>::Rotate1 (index arrays have wrong dimensions)"); 

	for (short j = 0; j < short(m_ny); j++)
	{	
		for (short i = 0; i < short(m_nx); i++)
		{
			if (xaJ0[j][i] >= 0 && xaI0[j][i] >= 0)
				xaResult[j][i] = m_rXArray2D[xaJ0[j][i]][xaI0[j][i]];
			else
				xaResult[j][i] = Backgr;
		}
	}
}


//---------------------------------------------------------------------------
//Function XArray2DSpln<T>::Magnify
//
//	(De)Magnifies 2D array by factor M with respect to the centre of dilation
//
/*!
	\brief		(De)Magnifies 2D array by factor M with respect to the centre of dilation
	\param		xaResult Resultant magnified 2D array
	\param		dblM Magnification factor
	\param		dblYCentre Y-coordinate of the centre of rotation
	\param		dblXCentre X-coordinate of the centre of rotation
	\param		Backgr Value for filling "background" areas around the (de)magnified array
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function (de)magnifies 2D array with respect to a specified point by a specified factor
		using bilinear interpolation. Only a portion of the (de)magnified image is calculated
		which is inside the original image boundaries
		
*/	
template <class T> void xar::XArray2DSpln<T>::Magnify(XArray2D<T>& xaResult, double dblM, double dblYCentre, double dblXCentre, T Backgr) const
{
	if (dblM <= 0)
			throw std::invalid_argument("invalid_argument 'dblM' in XArray2DSpln<T>::Magnify (non-positive magnification factor)");

	long j0, i0;
	double s, t, s1, t1;

	dblM = 1.0 / dblM;
	dblYCentre *= (1.0 - dblM);
	dblXCentre *= (1.0 - dblM);

	double ystM = m_yst * dblM;
	double xstM = m_xst * dblM;
	double yloM = m_ylo * dblM + dblYCentre - m_ylo;
	double xloM = m_xlo * dblM + dblXCentre - m_xlo;

	if (xaResult.GetDim1() != m_ny || xaResult.GetDim2() != m_nx) xaResult.Resize(m_ny, m_nx);

	// create new head (note that the new head is created even if the array did not have a had before the magnification)
	IXAHWave2D* pHead1 = CreateWavehead2D();
	pHead1->SetData(m_wl, m_ylo, m_yhi, m_xlo, m_xhi);
	xaResult.SetHeadPtr(pHead1);
	
	for (index_t j = 0; j < m_ny; j++)
	{
		s = (yloM + ystM * j) * m_ayst;

		j0 = long(floor(s));
		if (j0 < 0)
		{ 
			for (index_t i = 0; i < m_nx; i++) xaResult[j][i] = Backgr; 
			continue; 
		}
		else if (j0 >= long(m_ny1))
				if (j0 ==long( m_ny1) && double(j0) == s) j0--;
				else { for (index_t i = 0; i < m_nx; i++) xaResult[j][i] = Backgr; continue; }

		s -= j0;
		s1 = 1.0 - s;
		
		for (index_t i = 0; i < m_nx; i++)
		{
			t = (xloM + xstM * i) * m_axst;

			i0 = long(floor(t));
			if (i0 < 0) 
			{ 
				xaResult[j][i] = Backgr; 
				continue; 
			}
			else if (i0 >= long(m_nx1))
					if (i0 == m_nx1 && double(i0) == t) i0--;
					else { xaResult[j][i] = Backgr; continue; }

			t -= i0;
			t1 = 1.0 - t;

			xaResult[j][i] = (T)(t1 * (s1 * m_rXArray2D[j0][i0] + s * m_rXArray2D[j0+1][i0]) + 
				t * (s1 * m_rXArray2D[j0][i0+1] + s * m_rXArray2D[j0+1][i0+1]) + m_a05);
		}
	}
}


//---------------------------------------------------------------------------
//Function XArray2DSpln<T>::Tilt
//
//	Tilts 2D array by the specified angle around the tilt axis
//
/*!
	\brief		Tilts 2D array by the specified angle around the tilt axis
	\param		xaResult Resultant tilted 2D array
	\param		dblYAxis Y-coordinate of a point on a tilt axis
	\param		dblXAxis X-coordinate of a point on a tilt axis
	\param		dblAngleYAxis Angle between the tilt axis and Y-axis (degrees)
	\param		dblTiltAngle Tilt angle (degrees)
	\param		Backgr Value for filling "background" areas around the tilted array
	\exception	std::runtime_error and derived exceptions can be thrown indirectly by the functions
				called from inside this function
	\return		\a None
	\par		Description:
		This function tilts (linearly demagnifies) 2D array by the specified angle around the tilt axis
		which passes through a given point making a specified angle with Y-axis
		(demagnification factor is  equal to M=cos(dblTiltAngle)).
		The function uses bilinear interpolation. Only a portion of the tilted image is calculated
		which is inside the original image boundaries
		
*/	
template <class T> void xar::XArray2DSpln<T>::Tilt(XArray2D<T>& xaResult, double dblYAxis, double dblXAxis, double dblAngleYAxis, double dblTiltAngle, T Backgr) const
{
	dblTiltAngle *= PI180;
	double M = cos(dblTiltAngle);
	
	if (M == 0) 
		throw std::invalid_argument("invalid_argument 'dblTiltAngle' in XArray2DSpln<T>::Tilt (90*n degrees tilt is not allowed)");
	else
		M = 1.0 / M;
	
	double M1 = 1.0 - M;
	dblAngleYAxis *= PI180;
	double sg = sin(dblAngleYAxis);
	double cg = cos(dblAngleYAxis);
	double a2 = M1 * cg * cg;
	double b2 = M1 * sg * sg;
	double ab = M1 * sg * cg;
	double yshift = b2 * dblYAxis + ab * dblXAxis;
	double xshift = a2 * dblXAxis + ab * dblYAxis;
	double yfactor = M + a2;
	double xfactor = M + b2;

	if (xaResult.GetDim1() != m_ny || xaResult.GetDim2() != m_nx) xaResult.Resize(m_ny, m_nx);

	// create new head (note that the new head is created even if the array did not have a had before the tilt)
	IXAHWave2D* pHead1 = CreateWavehead2D();
	pHead1->SetData(m_wl, m_ylo, m_yhi, m_xlo, m_xhi);
	xaResult.SetHeadPtr(pHead1);

	double ystM = m_yst * yfactor;
	double ystab = -m_yst * ab;
	double yloM = m_ylo * yfactor - m_xlo * ab + yshift - m_ylo;
	double xstM = m_xst * xfactor;
	double xstab = -m_xst * ab;
	double xloM = m_xlo * xfactor - m_ylo * ab + xshift - m_xlo;

	long j0, i0;
	double y0, x0, s, t, s1, t1;
	
	for (index_t j = 0; j < m_ny; j++)
	{
		x0 = xloM + ystab * j;
		y0 = yloM + ystM * j;

		for (index_t i = 0; i < m_nx; i++)
		{
			s = (y0 + xstab * i) * m_ayst;
			t = (x0 + xstM * i) * m_axst;

			j0 = long(floor(s));
			if (j0 < 0) 
			{ 
				xaResult[j][i] = Backgr; 
				continue; 
			}
			else if (j0 >= long(m_ny1))
					if (j0 == m_ny1 && double(j0) == s) 
						j0--;
					else 
					{ 
						xaResult[j][i] = Backgr; 
						continue; 
					}


			s -= j0;
			s1 = 1.0 - s;

			i0 = long(floor(t));
			if (i0 < 0) 
			{ 
				xaResult[j][i] = Backgr; 
				continue; 
			}
			else if (i0 >= long(m_nx1))
					if (i0 == m_nx1 && double(i0) == t) 
						i0--;
					else 
					{ 
						xaResult[j][i] = Backgr; 
						continue; 
					}


			t -= i0;
			t1 = 1.0 - t;

			xaResult[j][i] = (T)(t1 * (s1 * m_rXArray2D[j0][i0] + s * m_rXArray2D[j0+1][i0]) + 
				t * (s1 * m_rXArray2D[j0][i0+1] + s * m_rXArray2D[j0+1][i0+1]) + m_a05);
		}
	}
}


//---------------------------------------------------------------------------
//	EXPLICIT TEMPLATE INSTANTIATION STATEMENTS
//

// The following causes explicit instantiation and catches compile errors that otherwise remain
// uncaught. Compiler warnings may be generated saying that classes have been already instantiated,
// but this is obviously not completely true, as plenty of errors remain uncaught when the
// following lines are commented out.
template class xar::XArray2DSpln<char>;
template class xar::XArray2DSpln<short>;
template class xar::XArray2DSpln<long>;
template class xar::XArray2DSpln<float>;
template class xar::XArray2DSpln<double>;

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