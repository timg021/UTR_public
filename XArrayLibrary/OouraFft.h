#pragma once

#include "XArray1D.h"
#include "XA_head2.h"


extern "C"
{
    void cdft2d_double(int, int, int, double **, double *, int *, double *);
    void rdft2d_double(int, int, int, double **, double *, int *, double *);
    void rdft2dsort_double(int, int, int, double **);
    void ddct2d_double(int, int, int, double **, double *, int *, double *);
    void ddst2d_double(int, int, int, double **, double *, int *, double *);
	
	void cdft2d_float(int, int, int, float **, float *, int *, float *);
    void rdft2d_float(int, int, int, float **, float *, int *, float *);
    void rdft2dsort_float(int, int, int, float **);
    void ddct2d_float(int, int, int, float **, float *, int *, float *);
    void ddst2d_float(int, int, int, float **, float *, int *, float *);

	void cdft_double(int, int, double *, int *, double *);
    void rdft_double(int, int, double *, int *, double *);
    void ddct_double(int, int, double *, int *, double *);
    void ddst_double(int, int, double *, int *, double *);
    void dfct_double(int, double *, double *, int *, double *);
    void dfst_double(int, double *, double *, int *, double *);
	
	void cdft_float(int, int, float *, int *, float *);
    void rdft_float(int, int, float *, int *, float *);
    void ddct_float(int, int, float *, int *, float *);
    void ddst_float(int, int, float *, int *, float *);
    void dfct_float(int, float *, float *, int *, float *);
    void dfst_float(int, float *, float *, int *, float *);
}


namespace xar
{

template <class T> class OouraFft
{
public:
	//! Transform direction
	enum _eDir
	{
		eDirFwd,	//!< Forward transform
		eDirInv,	//!< Inverse transform
	};
	//! Transform types
	enum _eType
	{
		eTypeComplexFft,	//!< Complex FFT
		eTypeRealFft,		//!< Real FFT
		eTypeDct,			//!< Discrete cosine transform
		eTypeDst,			//!< Discrete sine transform
		eTypeCos,			//!< Cosine transform
		eTypeSin,			//!< Sine transform
		eTypeMax
	};

public:
	OouraFft(void)
	{
	}


	// NOTE: the absence of appropriate specialisations of the following function
	// will prevent instantiation of OouraFft<T> objects for types T other than float or double
	//! Returns the xar::_eValueType corresponding to T
	static _eValueType GetValuetype(void);

	//! Complex 1-D forward FFT (in-place pointer form)
	void Complex1D(std::complex<T>* pcdData, size_t ulLength, _eDir eDir);

	//! Real 1-D forward FFT (in-place pointer form)
	void Real1D(T* pcdData, size_t ulLength, _eDir eDir);

	//! Complex 2-D forward FFT (in-place pointer form)
	void Complex2D(std::complex<T>* pcdData,
		size_t ulHeight, size_t ulWidth, _eDir eDir);

	//! Real 2-D forward FFT (in-place pointer form)
	void Real2D(T* pcdData, T *pcdSpeq, size_t ulHeight, 
		size_t ulWidth, _eDir eDir);


private:
	std::vector<int> m_vecBitRev;	//!< Bit reversal work area
	std::vector<T> m_vecTrigTbl;	//!< Sine/cosine lookup table
	std::vector<T> m_vecWork;		//!< Work area for 2-D transforms
	std::vector<T*> m_vecRow;		//!< Pointers to 2-D rows
	std::vector<T**> m_vecBlock;	//!< Pointers to 3-D blocks
	_eType m_eType;					//!< Latest transform type

	//! Prepare the 1-D work area and trigonometry table area
	void PrepareWorkArea1D(_eType eType, _eDir eDir, size_t ulDataCnt);

	//! Prepare a local vector of 2-D row pointers
	void PackRows(T* pdblData, size_t ulRowStep, size_t ulRowCnt);

	//! Prepare the 2-D work area and trigonometry table area
	void PrepareWorkArea2D(_eType eType, _eDir eDir, size_t ulHeight, 
		size_t ulWidth);
};
}


//	CComplex 2-D FFT (in-place pointer form)
//
/*!
	\param		pcdData		Transform data
	\param		ulHeight	Y-dimension of data
	\param		ulWidth		X-dimension of data
	\param		eDir		direction of the transform
	\remarks
		Data length of signal and transform arrays must be a power of two.
*/
template <class T> void xar::OouraFft<T>::Complex2D(std::complex<T>* pcdData,
	size_t ulHeight, size_t ulWidth, _eDir eDir)
{
	int nSign;

	if(eDir == eDirInv)
		nSign = -1;
	else
		nSign = 1;

	// Prepare the work area
	PrepareWorkArea2D(eTypeComplexFft, eDirFwd, ulHeight, ulWidth);
	
	// Prepare a local vector of row pointers
	PackRows(reinterpret_cast<T*>(pcdData), 2 * ulWidth, ulHeight);

	// Forward transform
	switch(GetValuetype())
	{
	case eXAFloat:
		::cdft2d_float(static_cast<int>(ulHeight),
			static_cast<int>(2 * ulWidth), nSign,
			(float **) &m_vecRow[0], (float *) &m_vecWork[0], &m_vecBitRev[0], 
			(float *) &m_vecTrigTbl[0]);
		break;

	case eXADouble:
		::cdft2d_double(static_cast<int>(ulHeight),
			static_cast<int>(2 * ulWidth), nSign,
			(double **) &m_vecRow[0], (double *) &m_vecWork[0], &m_vecBitRev[0], 
			(double *) &m_vecTrigTbl[0]);
		break;

	default:
		assert(false);
	}
}


//	Real 2-D FFT (in-place pointer form)
//
/*!
	\param		pcdData		Transform data
	\param		pcdSpeq		The Nyquist critical frequency values of the
							second frequency component (size = ulHeight + 2)
	\param		ulHeight	Y-dimension of data
	\param		ulWidth		X-dimension of data
	\param		eDir		direction of the transform
	\remarks
		Data length of signal and transform arrays must be a power of two.
*/
template <class T> void xar::OouraFft<T>::Real2D(T* pcdData, T *pcdSpeq,
		size_t ulHeight, size_t ulWidth, _eDir eDir)
{
	int nSign;
	size_t n1h, i;
	T x, y;

	if(eDir == eDirInv)
		nSign = -1;
	else
		nSign = 1;

	// Prepare the work area
	PrepareWorkArea2D(eTypeRealFft, eDirFwd, ulHeight, ulWidth);
	
	// Prepare a local vector of row pointers
	PackRows(reinterpret_cast<T*>(pcdData), ulWidth, ulHeight);


	if(eDir == eDirInv)	
	{
		// input ordering
		// take data from pcdData and pcdSpeq and shove them all into
		// pcdData in the weird Ooura ordering (there is no loss of information
		// because of symmetry
		// code is based on rdft2dsort except that we're using a Speq
		// array instead of an enlarged data array

		n1h = ulHeight >> 1;
		for (i = n1h + 1; i < ulHeight; i++) 
		{
            m_vecRow[i][0] = pcdSpeq[2*i + 1];
            m_vecRow[i][1] = pcdSpeq[2*i];
        }
        m_vecRow[0][1] = pcdSpeq[0];
        m_vecRow[n1h][1] = pcdSpeq[n1h*2];

	}

	// Forward transform
	switch(GetValuetype())
	{
	case eXAFloat:
		::rdft2d_float(static_cast<int>(ulHeight),
			static_cast<int>(ulWidth), nSign,
			(float **) &m_vecRow[0], (float *) &m_vecWork[0], &m_vecBitRev[0], 
			(float *) &m_vecTrigTbl[0]);
		break;

	case eXADouble:
		::rdft2d_double(static_cast<int>(ulHeight),
			static_cast<int>(ulWidth), nSign,
			(double **) &m_vecRow[0], (double *) &m_vecWork[0], &m_vecBitRev[0], 
			(double *) &m_vecTrigTbl[0]);
		break;

	default:
		assert(false);
	}

	if(eDir == eDirFwd)	
	{
		// output ordering
		// take data from pcdData that is arranged in the weird Ooura
		// ordering and re-arrange it in the normal, spilling over into
		// the Speq array (using the symmetry properties)
		// code is based on rdft2dsort except that we're using a Speq
		// array instead of an enlarged data array

		n1h = ulHeight >> 1;
		for (i = n1h + 1; i < ulHeight; i++) {
            y = m_vecRow[i][0];
            x = m_vecRow[i][1];
			pcdSpeq[2*i] = x;
			pcdSpeq[2*i + 1] = y;
			pcdSpeq[(ulHeight - i) * 2] = x;
			pcdSpeq[(ulHeight - i) * 2 + 1] = -y;
            m_vecRow[i][0] = m_vecRow[ulHeight - i][0];
            m_vecRow[i][1] = -m_vecRow[ulHeight - i][1];
        }
        pcdSpeq[0] = m_vecRow[0][1];
        pcdSpeq[1] = 0;
        m_vecRow[0][1] = 0;
        pcdSpeq[n1h * 2] = m_vecRow[n1h][1];
        pcdSpeq[n1h * 2 + 1] = 0;
        m_vecRow[n1h][1] = 0;

	}
}


//	Complex 1-D FFT (in-place pointer form)
//
/*!
	\param		pcdData		Transform data
	\param		ulLength	Number of elements in pdcData
	\param		eDir		direction of the transform
	\remarks
		Data length of signal and transform arrays must be a power of two.
*/
template <class T> void xar::OouraFft<T>::Complex1D(std::complex<T>* pcdData, 
					size_t ulLength, _eDir eDir)
{
	int nSign;

	if(eDir == eDirInv)
		nSign = -1;
	else
		nSign = 1;

	// Prepare the work area
	PrepareWorkArea1D(eTypeComplexFft, eDir, ulLength);
	
	// Forward transform
	switch(GetValuetype())
	{
	case eXAFloat:
		::cdft_float((int) (2 * ulLength), nSign, (float *) pcdData, 
			&m_vecBitRev[0], (float *) &m_vecTrigTbl[0]);
		break;

	case eXADouble:
		::cdft_double((int) (2 * ulLength), nSign, (double *) pcdData, 
			&m_vecBitRev[0], (double *) &m_vecTrigTbl[0]);
		break;

	default:
		assert(false);
	}
}

//	Real 1-D FFT (in-place pointer form)
//
/*!
	\param		pcdData		Transform data
	\param		ulLength	Number of elements in pdcData
	\param		eDir		direction of the transform
	\remarks
		Data length of signal and transform arrays must be a power of two.
*/
template <class T> void xar::OouraFft<T>::Real1D(T* pcdData, 
					size_t ulLength, _eDir eDir)
{
	int nSign;

	if(eDir == eDirInv)
		nSign = -1;
	else
		nSign = 1;

	// Prepare the work area
	PrepareWorkArea1D(eTypeRealFft, eDir, ulLength);
	
	// Forward transform
	switch(GetValuetype())
	{
	case eXAFloat:
		::rdft_float((int) ulLength, nSign, (float *) pcdData, 
			&m_vecBitRev[0], (float *) &m_vecTrigTbl[0]);
		break;

	case eXADouble:
		::rdft_double((int) ulLength, nSign, (double *) pcdData, 
			&m_vecBitRev[0], (double *) &m_vecTrigTbl[0]);
		break;

	default:
		assert(false);
	}
}



//! Prepare a local vector of 2-D row pointers
template <class T> void xar::OouraFft<T>::PackRows(T* pdblData, size_t ulRowStep,
						 size_t ulRowCnt)
{
	// Resize the 2-D row vector
	m_vecRow.resize(ulRowCnt);
	// Populate the vector with start-of-row pointers
	T** ppdblRow = &m_vecRow[0];
	size_t ulRowIdx;
	size_t ulRowOff = 0;
	for (ulRowIdx = 0; ulRowIdx < ulRowCnt; ulRowIdx++)
	{
		ppdblRow[ulRowIdx] = pdblData + ulRowOff;
		ulRowOff += ulRowStep;
	}
}


//---------------------------------------------------------------------------
//! Prepare the 1-D work area and trigonometry table area
template <class T> void xar::OouraFft<T>::PrepareWorkArea1D(_eType eType, _eDir eDir,
								  size_t ulDataCnt)
{
	size_t ulWorkDim = 0;
	size_t ulWorkLen = 0;
	size_t ulTrigDim = 0;
	size_t ulInnerDim = ulDataCnt;
	// Act on the transform type
	switch (eType)
	{
	case eTypeComplexFft:	// Complex FFT
		// ulTrigDim = n/2
		ulTrigDim = ulDataCnt / 2;
		// Calculate the required bit reversal work area length
		ulWorkLen = (size_t) ceil(2.0 + sqrt((double) ulDataCnt));
		break;

	case eTypeRealFft:		// Real FFT
		// ulTrigDim = n/2
		ulTrigDim = ulDataCnt / 2;
		// Calculate the required bit reversal work area length
		ulWorkLen = (size_t) ceil(2.0 + sqrt((double) ulDataCnt));
		break;

	case eTypeDct:			// Discrete cosine transform
		// ulWorkDim = n/2
		ulWorkDim = ulInnerDim / 2;
		// ulTrigDim = 5*n/4
		ulTrigDim = 5 * ulInnerDim / 4;
		// Calculate the required bit reversal work area length
		ulWorkLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		break;
	case eTypeDst:			// Discrete sine transform
		// ulWorkDim = n/2
		ulWorkDim = ulInnerDim / 2;
		// ulTrigDim = 5*n/4
		ulTrigDim = 5 * ulInnerDim / 4;
		// Calculate the required bit reversal work area length
		ulWorkLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		break;
	case eTypeCos:			// Cosine transform
		// Adjust the transform size 
		ulInnerDim -= 1;
		// ulWorkDim = n/4
		ulWorkDim = ulInnerDim / 4;
		// ulTrigDim = 5*n/8
		ulTrigDim = 5 * ulInnerDim / 8;
		// Calculate the required bit reversal work area length
		ulWorkLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		break;
	case eTypeSin:			// Sine transform
		// ulWorkDim = n/4
		ulWorkDim = ulInnerDim / 4;
		// ulTrigDim = 5*n/8
		ulTrigDim = 5 * ulInnerDim / 8;
		// Calculate the required bit reversal work area length
		ulWorkLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		break;
		break;
	default:
		break;
	}
		
	// Resize the work area and trigonometry table area
	if (ulWorkLen)
	{
		// Resize the bit reversal work area and force table initialisation
		if (m_vecBitRev.size() != ulWorkLen || m_eType != eType)
		{
			m_vecBitRev.resize(ulWorkLen);
			m_vecBitRev[0] = 0;
			m_eType = eType;
		}
		m_vecTrigTbl.resize(ulTrigDim);
	}
}


//---------------------------------------------------------------------------
//! Prepare the 2-D work area and trigonometry table area
template <class T> void xar::OouraFft<T>::PrepareWorkArea2D(_eType eType, _eDir eDir,
								  size_t ulHeight, size_t ulWidth)
{
	size_t ulWorkDim = 0;
	size_t ulWorkLen = 0;
	size_t ulBitRevLen = 0;
	size_t ulTrigDim = 0;
	size_t ulInnerDim = ulWidth;
	// Act on the transform type
	switch (eType)
	{
	case eTypeComplexFft:	// Complex FFT
		// ulWorkDim = max(n1, n2)
		ulWorkDim = max<size_t>(ulWidth, ulHeight);
		// ulTrigDim >= max(n1/2, n2/2)
		ulTrigDim = max<size_t>(ulWidth / 2, ulHeight / 2);
		// Calculate the required bit reversal work area length
		ulBitRevLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		// ulWorkLen = 8 * n1
		ulWorkLen = 8 * ulHeight;
		break;
	case eTypeRealFft:		// Real FFT
		// ulWorkDim = max(n1, n2/2)
		ulWorkDim = max<size_t>(ulHeight, ulWidth / 2);
		// ulTrigDim >= max(n1/2, n2/4) + n2/4
		ulTrigDim = max<size_t>(ulHeight / 2, ulWidth / 4) + ulWidth / 4;
		// Calculate the required bit reversal work area length
		ulBitRevLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		// ulWorkLen = 8 * n1
		ulWorkLen = 8 * ulHeight;
		break;
	case eTypeDct:			// Discrete cosine transform
		// Adjust the transform size 
		if (eDir == eDirInv)
		{
			ulInnerDim *= 2;
		}
		// ulWorkDim = max(n1/2, n2/2)
		ulWorkDim = max<size_t>(ulWidth / 2, ulInnerDim / 2);
		// ulTrigDim >= max(3*n1/2, 3*n2/2)
		ulTrigDim = max<size_t>(3 * ulWidth / 2, 3 * ulInnerDim / 4);
		// Calculate the required bit reversal work area length
		ulBitRevLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		// ulWorkLen = 4 * n1
		ulWorkLen = 4 * ulHeight;
		break;
	case eTypeDst:			// Discrete sine transform
		// Adjust the transform size 
		if (eDir == eDirInv)
		{
			ulInnerDim *= 2;
		}
		// ulWorkDim = max(n1/2, n2/2)
		ulWorkDim = max<size_t>(ulWidth / 2, ulInnerDim / 2);
		// ulTrigDim >= max(3*n1/2, 3*n2/2)
		ulTrigDim = max<size_t>(3 * ulWidth / 2, 3 * ulInnerDim / 4);
		// Calculate the required bit reversal work area length
		ulBitRevLen = static_cast<size_t>(ceil(2.0 + sqrt(static_cast<double>(ulWorkDim))));
		// ulWorkLen = 4 * n1
		ulWorkLen = 4 * ulHeight;
		break;
	default:
		break;
	}
	// Resize the work area and trigonometry table area
	if (ulBitRevLen)
	{
		// Resize the bit reversal work area and force table initialisation
		if (m_vecBitRev.size() != ulBitRevLen || m_eType != eType)
		{
			m_vecBitRev.resize(ulBitRevLen);
			m_vecBitRev[0] = 0;
			m_eType = eType;
		}
		m_vecWork.resize(ulWorkLen);
		m_vecTrigTbl.resize(ulTrigDim);
	}
}