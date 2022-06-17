//Header XA_nrrand.h
//
//
//	HEADER FILE TITLE:
//
//		Modified random generators and related routines from Numerical Recipes
//
/*!
	\file		XA_nrrand.h
	\brief		Modified random generators and related routines from Numerical Recipes
	\par		Description:
		This header contains various random number generators from the 2nd editions
		of "Numerical Recipes in C"	modified.
*/
#if !defined XA_NRRAND_H
#define XA_NRRAND_H

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
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

namespace // anonymous
{

	#define IA 16807
	#define IM 2147483647
	#define AM (1.0/IM)
	#define IQ 127773
	#define IR 2836
	#define NTAB 32
	#define NDIV (1+(IM-1)/NTAB)
	#define EPS 1.2e-7
	#define RNMX (1.0-EPS)

	//! Calculates random numbers with good statistical properties
	// Excellent random number generator, an improvement of 'ran0';
	// Passes all statistical tests until the number of calls >1E+8
	// NOTE!: Most recommended by NumRecipes.
	//
	// "Minimal" random number generator of Park and Miller with Bays-Durham
	// shuffle and added safeguards. Returns a uniform random deviate
	// 0.0<rand1<1.0. Call with 'idum' a NEGATIVE INTEGER to initialize;
	// thereafter do not alter 'idum' between successive deviates in a sequence.
	// RNMX should approximate the largest floating value that is less than 1.
	//
	// Numerical Recipes, 2nd edition, p.280.
	float ran1(long *idum)
	{
		int j;
		long k;
		static long iy=0;
		static long iv[NTAB];
		float temp;

		if (*idum <= 0 || !iy) {
			if (-(*idum) < 1) *idum=1;
			else *idum = -(*idum);
			for (j=NTAB+7;j>=0;j--) {
				k=(*idum)/IQ;
				*idum=IA*(*idum-k*IQ)-IR*k;
				if (*idum < 0) *idum += IM;
				if (j < NTAB) iv[j] = *idum;
			}
			iy=iv[0];
		}
		k=(*idum)/IQ;
		*idum=IA*(*idum-k*IQ)-IR*k;
		if (*idum < 0) *idum += IM;
		j=iy/NDIV;
		iy=iv[j];
		iv[j] = *idum;
		if ((temp=(float)(AM*iy)) > RNMX) return (float)RNMX;
		else return temp;
	}

	#undef IA
	#undef IM
	#undef AM
	#undef IQ
	#undef IR
	#undef NTAB
	#undef NDIV
	#undef EPS
	#undef RNMX


	//! Returns a normally distributed deviate with zero mean and unit variance
	// using ran1(idum) as the source of uniform deviates
	//
	// Numerical Recipes in C, 2nd ed., p.289
	static float gasdev(long *idum)
	{
		static int iset=0;
		static float gset;
		float fac,rsq,v1,v2;

		if  (iset == 0) {
			do {
				v1=(float)(2.0*ran1(idum)-1.0);
				v2=(float)(2.0*ran1(idum)-1.0);
				rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac=(float)sqrt(-2.0*log(rsq)/rsq);
			gset=v1*fac;
			iset=1;
			return v2*fac;
		} else {
			iset=0;
			return gset;
		}
	}


	//! Returns the value ln[Gamma(xx)] for xx>0.
	//
	// Numerical Recipes in C, 2nd ed., p.214
	float gammln(float xx)
	{
		double x,y,tmp,ser;
		static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
		int j;

		y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*log(tmp);
		ser=1.000000000190015;
		for (j=0;j<=5;j++) ser += cof[j]/++y;
		return (float)(-tmp+log(2.5066282746310005*ser/x));
	}


	//! Returns a random deviate from a Poisson distribution of a given mean
	// Returns as a floating-point number an integer value that is a random
	// deviate from a Poisson distribution of mean xm, using ran1(idum) as
	// a source of uniform random deviates.
	// Also uses the function gammln() to calculate natural log of Gamma-function
	//
	// NOTE!: for Poisson distribution mean = dispersion = xm
	//
	// Numerical Recipes in C, 2nd ed., p.294
	float poidev(float xm, long *idum)
	{
		static float sq,alxm,g,oldm=(-1.0);
		float em,t,y;

		if (xm < 12.0) {
			if (xm != oldm) {
				oldm=xm;
				g=(float)exp(-xm);
			}
			em = -1;
			t=1.0;
			do {
				++em;
				t *= ran1(idum);
			} while (t > g);
		} else {
			if (xm != oldm) {
				oldm=xm;
				sq=(float)sqrt(2.0*xm);
				alxm=(float)log(xm);
				g=xm*alxm-gammln(xm+(float)1.0);
			}
			do {
				do {
					y=(float)tan(3.14159265359*ran1(idum));
					em=sq*y+xm;
				} while (em < 0.0);
				em=(float)floor(em);
				t=(float)(0.9*(1.0+y*y)*exp(em*alxm-gammln(em+(float)1.0)-g));
			} while (ran1(idum) > t);
		}
		return em;
	}


} // namespace anonymous closed


//---------------------------------------------------------------------------
//	CLASS DECLARATIONS
//
//---------------------------------------------------------------------------
//	TEMPLATE MEMBER DEFINITIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL DATA DECLARATIONS
//
//---------------------------------------------------------------------------
//	EXTERNAL FUNCTION PROTOTYPES
//


#endif  // XA_NRRAND_H
/////////////////////////////////////////////////////////////////////////////
//
//	End of Header File
//