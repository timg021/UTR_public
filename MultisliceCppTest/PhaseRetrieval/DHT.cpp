// DHT.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <complex.h>
#include <chrono>
#include <omp.h>
#include <filesystem>

#include "IXAHWave.h"
#include "XAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "fftwd3drc.h"
#include "XA_spln2.h"
#include "XA_spln3.h"
#include "XA_iwfr.h"
#include "XA_tie.h"
#include "XA_tiff.h"
#include "XA_move2.h"

using namespace xar;

#define NATOMMAX 500000 // maximum number of atomic positions that can be saved in an XYZ file

// the following #define is a temporary fork of this program allowing one to correlate two sets of image frames (this fork is enabled when the macro CORRELATE_FRAMES is defined)
// The base name of the reference set of frames is given by the macro CORRELATE_FRAMES itself, as in the example below
// The output 2-column ASCI correlation file is always saved under the name "0Correlate.txt" in the same directory from which this executable is launched
//#define CORRELATE_FRAMES "Z:\\New_experimental_data_Sep2021\\ApoFramesFrom7A4M_assembly_DW\\Hand_flipped\\bpo.grd" 

int main(int argc, char* argv[])
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now(), start_timeNow;
	long int liFileReadTime(0), liFileWriteTime(0), li2DFFTtime(0), li3DFFTtime(0), liBackPropTime(0);

	try
	{
		printf("\nStarting Differential Holographic Tomography program ...");
		constexpr double dAtomWidth = 1.0; // "average width" of the potential distribution around a single atom in A, this may change in the future
		constexpr double dGamma = 0.1; // scaling parameter required in order to bring the maximum contrast for the realistic potentials to the theoretical Gaussian case 
		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;
		vector<Pair> v2shifts;
		vector<Pair> vastigm;

		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];

		string sInputParamFile("PhaseRetrieval.txt");
		if (argc > 1) sInputParamFile = argv[1];

		FILE* ff0 = fopen(sInputParamFile.c_str(), "rt");
		if (!ff0) throw std::runtime_error(string("Error: cannot open parameter file " + sInputParamFile + ".").c_str());
		else printf("\nReading input parameter file %s ...", sInputParamFile.c_str());

		// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
		while (true)
		{
			fgets(cline, 1024, ff0);
			if (!(cline[0] == '/' && cline[1] == '/')) break;
		}

		strtok(cline, "\n"); // 1. Verbose output during execution? Yes = 1, No = 0
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading verbose output parameter from input parameter file.");
		bool bVerboseOutput(true); // if this is TRUE, additional information is printed during execution
		(atoi(cparam) == 0 || atoi(cparam) == 1) ? bVerboseOutput = (bool)atoi(cparam) : throw std::runtime_error("Error: verbose output parameter must be 0 or 1 in input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2. Input file with rotation angles and defocus distances
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading file name with rotation angles and defocus distances from input parameter file.");
		string defocfile(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3.Optional file with indexes and weights for subseries selection or NONE
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading the name of a file with indexes and weights for subseries selection from input parameter file.");
		string subseriesfile(cparam);
		bool bSelectFrames(false);
		if (subseriesfile != string("NONE"))
		{
			if (GetFileExtension(subseriesfile) != string(".TXT"))
				throw std::runtime_error("Error: unrecognised filename extension in the file with indexes and weights for subseries selection in input parameter file.");
			bSelectFrames = true;
		}

		printf("\nReading defocus parameter file %s ...", defocfile.c_str());
		bool bRelion(false);
		double dTotalWeight(0);
		vector<index_t> vndefocusFull;
		vector<int> sfIndexes;
		vector<float> sfWeights;
		if (GetFileExtension(defocfile) == string(".TXT"))
			ReadDefocusParamsFile(defocfile, v2angles, vvdefocus, bVerboseOutput);
		else
			if (GetFileExtension(defocfile) == string(".RELIONNEW"))
			{
				ReadRelionDefocusParamsFile(defocfile, v2angles, vvdefocus, vastigm, v2shifts, bVerboseOutput);
				bRelion = true;
				if (bSelectFrames)
				{
					vndefocusFull.resize(v2angles.size()); // vector of numbers of defocus planes at different illumination angles
					for (index_t i = 0; i < v2angles.size(); i++) vndefocusFull[i] = vvdefocus[i].size();
					XArray1D<float> rXAr1D, rYAr1D;
					XArData::ReadFileDAT2(rXAr1D, rYAr1D, subseriesfile.c_str(), 0.02);
					if (rXAr1D.size() == 0) throw std::runtime_error("Error: the number of entries in the subseries selection file is equal to zero.");
					else if (rXAr1D.size() > v2angles.size()) throw std::runtime_error("Error: the number of entries in the subseries selection file is larger than the number of entries in the defocus file.");
					sfIndexes.resize(rXAr1D.size()); sfWeights.resize(rYAr1D.size());
					int sfMin = int(rXAr1D[0]), sfMax = int(rXAr1D[0]);
					for (int i = 0; i < rXAr1D.size(); i++)
					{
						sfIndexes[i] = int(rXAr1D[i]);
						if (sfIndexes[i] > sfMax) sfMax = sfIndexes[i];
						else if (sfIndexes[i] < sfMin) sfMin = sfIndexes[i];
					}
					if (sfMax > v2angles.size() - 1) throw std::runtime_error("Error: some frame indexes in the subseries selection file are larger than the number of entries - 1 in the defocus file.");
					if (sfMin < 0) throw std::runtime_error("Error: some frame indexes in the subseries selection file are negative.");
					for (int i = 0; i < rXAr1D.size(); i++)
					{
						sfWeights[i] = rYAr1D[i];
						dTotalWeight += rYAr1D[i];
					}

					vector<Pair> v2anglesSel, v2shiftsSel, vastigmSel;
					vector<vector <Pair> > vvdefocusSel;
					for (int i : sfIndexes)
					{
						v2anglesSel.push_back(v2angles[i]);
						v2shiftsSel.push_back(v2shifts[i]);
						vastigmSel.push_back(vastigm[i]);
						vvdefocusSel.push_back(vvdefocus[i]);
					}
					v2angles = v2anglesSel;
					v2shifts = v2shiftsSel;
					vastigm = vastigmSel;
					vvdefocus = vvdefocusSel;

					if (dTotalWeight == 0) throw std::runtime_error("Error: the sum of weights in the second column of the subseries selection file is equal to zero.");
					printf("\nThere are %zd entries in the subseries selection file.", v2angles.size());
				}
			}
			else throw std::runtime_error("Error: unrecognised filename extension of the defocus file in the input parameter file.");
		index_t nangles = v2angles.size(); // number of rotation steps 
		vector<index_t> vndefocus(nangles); // vector of numbers of defocus planes at different illumination angles
		index_t ndefocusmin, ndefocusmax;
		vndefocus[0] = vvdefocus[0].size();
		ndefocusmin = ndefocusmax = vndefocus[0];
		for (index_t i = 1; i < nangles; i++)
		{
			vndefocus[i] = vvdefocus[i].size();
			if (vndefocus[i] < ndefocusmin) ndefocusmin = vndefocus[i];
			if (vndefocus[i] > ndefocusmax) ndefocusmax = vndefocus[i];
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Input filename base of defocus series of the sample in TIFF, GRD, GRC or RAW format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		bool bRAWinput(false), bTIFFinput(false), bGRDinput(false), bGRCinput(false);
		string strTemp = GetFileExtension(filenamebaseIn);
		if (strTemp == string(".TIFF") || strTemp == string(".TIF")) bTIFFinput = true;
		else if (strTemp == string(".GRC")) bGRCinput = true;
		else if (strTemp == string(".RAW")) bRAWinput = true;
		else if (strTemp == string(".GRD")) bGRDinput = true; // this value is not used below at the moment, as GRD is considered to be a default when all other input formats are false
		else throw std::runtime_error("Error: input filename extension must be TIF, GRD, GRC or RAW.");
		if (bGRCinput) // check that there is only one input defocused complex amplitude file per each illumination angle
		{
			for (index_t i = 0; i < nangles; i++)
				if (vndefocus[i] != 1)
					throw std::runtime_error("Error: only one input defocused complex amplitude file per each illumination angle is allowed.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. RAW file parameters: HeaderLength(bytes) Endianness(0 = little, 1 = big) ElementLength(bytes) 
		// for RAW files, these parameter determine how the images are read
		bool bBigEndian;
		index_t nHeaderLength, nElementLength;
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam2, cparam3, cparam4) != 4)
			throw std::runtime_error("Error reading RAW file parameters from input parameter file.");
		nHeaderLength = atoi(cparam2);
		bBigEndian = (bool)atoi(cparam3);
		nElementLength = atoi(cparam4);
		if (nElementLength != sizeof(float) && nElementLength != sizeof(double))
			throw std::runtime_error("Unrecognised ElementLength in input parameter file (only sizeof(float) or sizeof(double) are allowed).");
		printf("\nRAW file parameters: HeaderLength = %zd, Endianness = %d, ElementLength = %zd", nHeaderLength, bBigEndian, nElementLength);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6. Width(pix) and height(pix) of input images
		// for RAW files, these parameter determine how the images are read; for TIFF, GRD and GRC files, these parameters may be used for trimming the frames
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::runtime_error("Error reading Width Height from input parameter file.");
		index_t nx = atoi(cparam);
		index_t ny = atoi(cparam1);
		printf("\nDimensions of input images: width = %zd, height = %zd (pixels)", nx, ny);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Pixel size in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
			throw std::runtime_error("Error reading pixel size from input parameter file.");
		double xst = atof(cparam);
		if (xst <= 0)
			throw std::runtime_error("Error: pixel size must be positive in input parameter file.");
		double yst = xst; // this is a provision for possible future extension to non-square pixels
		double xlo(0.0);
		double xhi = xst * nx;
		double ylo(0.0);
		double yhi = yst * ny;
		printf("\nPixel size = %g (Angstroms)", xst);
		printf("\nPhysical boundaries of input images: Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g (Angstroms)", xlo, xhi, ylo, yhi);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //8.Centre of rotation x, y and z shifts in Angstroms
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4)
			throw std::runtime_error("Error reading ccentre of rotation x, y and z shifts from input parameter file.");
		double dxc = atof(cparam);
		double dyc = atof(cparam1);
		double dzc = atof(cparam2);
		printf("\nCentre of rotation shifts in Angstroms: dXc = %g, dYc = %g, dZc = %g (Angstroms)", dxc, dyc, dzc);
		if (abs(dxc) > (xhi - xlo) || abs(dyc) > (yhi - ylo))
			printf("\nWARNING: X or Y shifts of the centre of rotation will take the centre outside the reconstruction volume.");
		bool bRotCentreShift = (dxc != 0 || dyc != 0 || dzc != 0) ? true : false;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9.Input data normalization factors f1 f2 (input -> (input / f1) + f2)
		double dNormFactor1(1.0), dNormFactor2(0.0);
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::runtime_error("Error reading input data normalization factors from input parameter file.");
		dNormFactor1 = atof(cparam);
		dNormFactor2 = atof(cparam1);
		printf("\nInput data normalization factors: f1 = %g, f2 = %g", dNormFactor1, dNormFactor2);
		if (dNormFactor1 == 0) throw std::runtime_error("The first normalization factor cannot be zero (leads to division by zero).");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10. Wavelength in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in Angstroms
		printf("\nWavelength = %g (A)", wl);
		if (wl < 0 || wl > 1)
			throw std::runtime_error("Error: wavelength value appears to be wrong.");
		double EE; // incident electron energy in volts (recalculated from the wavelength below)
		constexpr double hp = 6.62607004e-34; // Planck's constant (m2 kg / s)
		constexpr double cc = 299792458; // speed of light (m / s)
		constexpr double ee = 1.602176634e-19; // electron charge (coulomb)
		constexpr double m0 = 9.1093837015e-31; // electron rest mass (kg)
		constexpr long double mc2 = m0 * cc * cc; // mc^2
		long double chl2 = (long double)(cc * cc * hp * hp) / (long double)(wl * wl * 1.e-20);
		long double abra = sqrt(mc2 * mc2 + chl2);
		EE = double((-mc2 + abra) / (long double)(ee));
		printf("\nIncident electron energy E = %g (keV)", EE / 1000.0);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam);
		if (aobj != 0) printf("\nObjective aperture (half-angle) = %g (mrad)", aobj);
		else  printf("\nObjective aperture is infinite.");
		if (aobj < 0 || aobj > 1000)
			throw std::runtime_error("Error: objective aperture value appears to be wrong.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12. Spherical aberrations Cs3 and Cs5 in mm
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading spherical aberrations from input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> Angstroms
		Cs5 *= 1.e+7; // mm --> Angstroms

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13. Absorption coefficient (fraction of the real part)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading absorptioin coefficient from input parameter file.");
		double asigma = atof(cparam);
		printf("\nAbsorption coefficient (fraction of the real part) = %g", asigma);


		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14.Phase retrieval method: 1 = IWFR, 2 = CTFL2, 3 = Min05LogAmp, 4 = vCTF, 5 = symmetrized-CTF (sCTF)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading phase retrieval method from input parameter file.");
		int nPhaseRetrieval = atoi(cparam);
		printf("\nPhase retrieval method number selected = %d", nPhaseRetrieval);
		if (asigma != 0 && !(nPhaseRetrieval == 4 || nPhaseRetrieval == 5))
			throw std::runtime_error("Error: selected phase retrieval method is incompatibe non-zero absorption coefficient in the input parameter file.");
		if (nPhaseRetrieval == 1 || nPhaseRetrieval == 2)
		{
			if (ndefocusmin < 2)
				throw std::runtime_error("Error: selected phase retrieval method is incompatibe with a single defocus distance appearing in some lines of the defocus/orientations parameter file.");
		}
		else 
			if (nPhaseRetrieval == 3 || nPhaseRetrieval == 4 || nPhaseRetrieval == 5)
			{
				if (ndefocusmax != 1) // multiple images at different defocuses and the same orientation can only be used with Min05LogAmp and CTF-based methods if placed in different lines of the defocus parameter file and then treated independently
					throw std::runtime_error("Error: selected phase retrieval method is incompatibe with multiple defocus distances in some lines of the defocus/orientations parameter file.");
				if (nPhaseRetrieval != 3 && bGRCinput == true) // implemented CTF-based reconstruction methods currently work only with image intensities
					throw std::runtime_error("Error: selected reconstruction (phase retrieval) method has not been implemented for GRC input files.");
			}
			else throw std::runtime_error("Error: unknown phase retrieval method in the input parameter file.");
		if (nPhaseRetrieval != 5)
		{
			if (log2(nx) != int(log2(nx)) || log2(ny) != int(log2(ny)))
				throw std::runtime_error("Error: Width and Height parameters in input parameter file must be integer powers of 2.");
		}
		else
		{
			if (nx != 2 * int(nx / 2) || ny != 2 * int(ny / 2))
				throw std::runtime_error("Error: Width and Height parameters in input parameter file must be even.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15.Apply Ewald sphere curvature correction: 0=No or 1=Yes
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading Ewald sphere curvature correction switch from input parameter file.");
		switch (atoi(cparam))
		{
		case 0:
			printf("\nEwald sphere curvature correction won't be applied.");
			break;
		case 1:
			printf("\nEwald sphere curvature correction will be applied.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for Ewald sphere curvature correction switch in input parameter file.");
		}
		bool bESCC = (bool)atoi(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 16. Save or not phase phase-retrieved defocused complex amplitudes or filter matrix in files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading save_or_not phase-retrieved defocused complex amplitudes switch from input parameter file.");
		int nSaveDefocCAmpsOrFilterMatrix = atoi(cparam);
		if (nPhaseRetrieval == 4)
			printf("\nNo phase-retrieved complex amplitudes in the image planes will be created or saved in files");
		else
		{
			if (nSaveDefocCAmpsOrFilterMatrix == 1)
				if (nPhaseRetrieval == 5) printf("\n3D filter matrix will be saved in GRD files");
				else printf("\nPhase-retrieved complex amplitudes in the image planes will be saved in GRD files");
			else 
				if (nSaveDefocCAmpsOrFilterMatrix == 0)
					printf("\nPhase-retrieved complex amplitudes in the image planes or the filter matrix will not be saved in files");
				else
					throw std::runtime_error("Error: save_or_not phase-retrieved defocused complex amplitudes switch must be 0 or 1.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17. Maximal number of IWFR iterations or critical SNR
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading maximal number of iterations from input parameter file.");
		int itermax;
		double CriticalSNR;
		if (nPhaseRetrieval == 1)
		{
			itermax = atoi(cparam);
			if (itermax < 1)
				throw std::runtime_error("Error: the maximal number of iterations must be >= 1.");
			printf("\nMaximal number of iterations = %d", itermax);
		}
		else if (nPhaseRetrieval == 5)
		{
			CriticalSNR = atof(cparam);
			if (CriticalSNR < 0)
				throw std::runtime_error("Error: the critical SNR value must be non-negative.");
			printf("\nCritical SNR = %g", CriticalSNR);
		}


		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 18. Minimal phase reconstruction error(IWFR), Tikhonov regularization parameter alpha(CTFL2, vCTF, cCTF, sCTF) or multiplicative factor (Min05LogAmp)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading minimal phase reconstruction error/regulariztion parameter from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nMinimal phase reconstruciton error (IWFR), Tikhonov regularization parameter (CTFL2, vCTF, cCTF, sCTF) or multiplicative factor (Min05LogAmp) = %g", epsilon);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19. Output defocus distances min and max in Angstroms
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading output defocus distances from input parameter file.");
		double zlo = atof(cparam); // minimum output defocus in Angstroms - !!! will be corrected with dzextra below
		double zhi = atof(cparam1); // maximum output defocus in Angstroms - !!! will be corrected with dzextra below 
		double zst = xst; // this is a provision for possible future extensions to non-qubic voxels
		if (zlo > zhi) std::swap(zlo, zhi);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 20. Extra defocus for 3D reconstruction in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading extra defocus for 3D reconstruction parameter from input parameter file.");
		double dzextra = atof(cparam);
		printf("\nExtra defocus for 3D reconstruction = %g (Angstroms)", dzextra);
		if (dzextra == 0)
		{
			if (nPhaseRetrieval == 1 || nPhaseRetrieval == 2 || nPhaseRetrieval == 3) 
				throw std::runtime_error("Error: extra defocus = 0, which is unexpected for the selected reconstruction (phase retrieval) method.");
		}
		else
			if (nPhaseRetrieval == 4 || nPhaseRetrieval == 5) 
				throw std::runtime_error("ERROR: extra defocus != 0, which is unexpected for the phase retrieval method selected above.");
		double zlodz(zlo), zhidz(zhi);
		zlodz += dzextra; zhidz += dzextra;
		printf("\nOutput defocus distances: min = %g, max = %g, step = %g (Angstroms)", zlo, zhi, zst);
		if (abs(dzc) > (zhi - zlo))
				printf("\nWARNING: Z shift of the centre of rotation will take the centre out of the reconstruction volume.");

		index_t nz = index_t((zhidz - zlodz) / zst + 0.5); // number of defocus planes to propagate to
		if (nz <= 0)
			throw std::runtime_error("Error: number of output defocus planes must be positive.");
		vector<double> voutdefocus(nz); // vector of output defocus distances
		for (index_t n = 0; n < nz; n++) voutdefocus[n] = zlodz + zst * n;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //21. Enforce symmetry: not_apply(0), post-apply(1), or distribute input orientations(2)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading enforcing symmetry mode from input parameter file.");
		int imodeSymmetry = atoi(cparam);
		switch (imodeSymmetry)
		{
		case 0:
			printf("\nKnown symmetry won't be applied.");
			break;
		case 1:
			printf("\nInput orientations will be distributed according to known invariant rotational positions.");
			if (nPhaseRetrieval == 5) 
				throw std::runtime_error("Error: this symmetry enforcing mode has not been implemented for the sCTF phase retrieval method.");
			break;
		case 2:
			printf("\nKnown symmetry will be applied after initial reconstruction or import of a 3D potential.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for symmetry enforcing mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 22.Input file with rotation angles enforcing symmetry
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading file name with rotation angles enforcing symmetry from input parameter file.");
		index_t nanglesSym{ 1 }, nanglesIn1Sym{ nangles };
		vector<Pair> v2anglesSym;
		vector<vector <Pair> > vvdefocusSym;
		if (imodeSymmetry != 0)
			if (GetFileExtension(string(cparam)) == string(".TXT"))
			{
				printf("\nReading symmetry enforcing file %s ...", cparam);
				ReadDefocusParamsFile(cparam, v2anglesSym, vvdefocusSym, bVerboseOutput);
				nanglesSym = v2anglesSym.size(); // number of rotation steps in the symmetry file
				if (imodeSymmetry == 1)
				{
					nanglesIn1Sym = nangles / nanglesSym; // number of illumination angles in a set assigned to one symmetry angle
					if (nanglesIn1Sym < 1) throw std::runtime_error("Error: number of illumination angles is smaller than the number of symmetry angles, so the distribution into groups is impossible.");
					else printf("\nThere will be % zd illumination angles in each symmetry subgroup.", nanglesIn1Sym);
				}
			}
			else
				if (!string(cparam).empty())
					throw std::runtime_error("Error: filename extension for a file with rotation angles for enforcing symmetry must be .txt or the parameter should be empty.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 23. 3D Laplacian filter mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading 3D Laplacian filter mode from input parameter file.");
		int imodeInvLaplace = atoi(cparam);
		switch (imodeInvLaplace)
		{
		case 0:
			printf("\n3D Laplacian filter won't be applied.");
			break;
		case 1:
			printf("\n3D Laplacian filter will be applied.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for 3D Laplacian filter mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 24. Regularization parameter for inverse 3D Laplacian
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading regularization parameter for inverse 3D Laplacian from input parameter file.");
		double alpha = atof(cparam);
		printf("\nRegularization parameter for inverse 3D Laplacian = %g", alpha);
		if (imodeInvLaplace && alpha < 0)
			throw std::runtime_error("Error: regularization parameter for 3D Laplacian filter must be non-negative.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 25. Low-pass filter width in Angstroms, background subtraction value and lower threshold level in Volts
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::runtime_error("Error reading low-pass filter width, background subtraction value and lower threshold level from input parameter file.");
		double dlpfiltersize = atof(cparam);
		double dBackground = atof(cparam1);
		double dThreshold = atof(cparam2);
		printf("\nLow-pass filter width for 3D potential = %g (Angstroms)", dlpfiltersize);
		printf("\nBackground subtraction value for 3D potential = %g (Volts)", dBackground);
		printf("\nLower threshold level for 3D potential = %g (Volts)", dThreshold);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 26. Multislice reprojection of the 3D electrostatic potential mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading reprojectin of 3D potential mode from input parameter file.");
		int imodeReproj = atoi(cparam);
		switch (imodeReproj)
		{
		case 0:
			printf("\nMultislice reprojection of the 3D electrostatic potential won't be applied.");
			break;
		case 1:
			printf("\nMultislice reprojection of the 3D electrostatic potential will be applied.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for the multislice reprojection mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 27. Slice thickness for multislice reprojection in Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading slice thickness for multislice reprojection from input parameter file.");
		double sliceTh = atof(cparam);
		if (imodeReproj != 0)
		{
			printf("\nSlice thickness for multislice reprojection = %g (Angstroms)", sliceTh);
			if (imodeReproj && sliceTh < (4.0 * zst))
				throw std::runtime_error("Error: slice thickness for multislice reprojection must be 4 x z_step or larger.");
		}
#ifdef CORRELATE_FRAMES
		imodeReproj = 1;
#endif // CORRELATE_FRAMES

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 28. Peak localization mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading peak localization mode from input parameter file.");
		int imodePeaks = atoi(cparam);
		switch (imodePeaks)
		{
		case 0:
			printf("\nPeak localization in the 3D electrostatic potential won't be applied.");
			break;
		case 1:
			printf("\nPeak localization in the 3D electrostatic potential will be applied.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for the peak localization mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 29. Transverse and longitudinal side lengths for peak localization (in_Angstroms)
		if (sscanf(cline, "%s %s %s", ctitle, cparam ,cparam1) != 3) throw std::runtime_error("Error reading transverse and longitudinal side lengths for peak localization from input parameter file.");
		double datomsizeXY = atof(cparam), datomsizeZ = atof(cparam1);
		printf("\nTransverse and longitudinal side lengths for peak localization = %g %g (Angstroms)", datomsizeXY, datomsizeZ);
		if (imodePeaks && (int(datomsizeXY / zst + 0.5) < 2 || int(datomsizeZ / zst + 0.5) < 2))
			throw std::runtime_error("Error: transverse and longitudinal side lengths for peak localization must be 2 x z_step or larger.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 30. Output file name base in GRD or TIFF format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		printf("\nFile name base for 3D potential = %s", filenamebaseOut.c_str());
		strTemp = GetPathFromFilename(filenamebaseOut, false);
		if (!std::filesystem::exists(strTemp))
			throw std::runtime_error("Error: the specified file folder for 3D potential does not seem to exist.");
		bool bTIFFoutput;
		if (GetFileExtension(filenamebaseOut) == string(".TIFF") || GetFileExtension(filenamebaseOut) == string(".TIF")) bTIFFoutput = true;
		else if (GetFileExtension(filenamebaseOut) == string(".GRD")) bTIFFoutput = false;
		else throw std::runtime_error("Error: output filename extension must be TIF ot GRD.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 31. Import and reprocess existing 3D_potential files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading import and reprocess existing 3D potential files switch from input parameter file.");
		int imode3DPotential = atoi(cparam);
#ifdef CORRELATE_FRAMES
		imode3DPotential = 1;
#endif // CORRELATE_FRAMES
		switch (imode3DPotential)
		{
		case 0:
			printf("\nThis program will reconstruct 3D potential from defocused images and save it in output files.");
			break;
		case 1:
			printf("\nThis program will import existing 3D potential from files and reprocess it.");
			break;
		default:
			throw std::runtime_error("Error: unknown value of import and reprocess existing 3D potential files switch in input parameter file.");
		}

		if (imode3DPotential) // only read in pre-caculated 3D potential from "output" files and reprocess them
		{
			bGRDinput = bTIFFinput = bGRCinput = bRAWinput = false;
			string strTemp = GetFileExtension(filenamebaseOut);
			if (strTemp == string(".TIFF") || strTemp == string(".TIF")) bTIFFinput = true;
			else if (strTemp == string(".GRD")) bGRDinput = true;
			else if (strTemp == string(".RAW")) bRAWinput = true;
			else throw std::runtime_error("Error: input filename extension (in this mode - it is taken from the output filename template) must be TIF, GRD or RAW.");
		}
		else
			if (nPhaseRetrieval != 3 && epsilon < 0)
				throw std::runtime_error("Error: minimal phase error / regulalization parameter must be non-negative, unless -0.5LogAmp phase retrieval method is selected.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // // 32. Folder name for auxiliary files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading folder name for auxiliary files from input parameter file.");
		string folderAux = cparam;
		printf("\nFolder for auxiliary file output = %s", folderAux.c_str());
		if (folderAux.rfind('\\') != folderAux.size() - 1 && folderAux.rfind('/') != folderAux.size() - 1) // the last character is not '\' and not '/'
			folderAux.append("/");
		if (!std::filesystem::exists(folderAux))
			throw std::runtime_error("Error: the specified auxiliary file folder does not seem to exist.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 33. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
#ifdef CORRELATE_FRAMES
		nThreads = 1; // for some reason the program crashed in the multithreaded mode, but it is fast enough in the single-threade mode
#endif // CORRELATE_FRAMES
		printf("\nNumber of parallel threads = %d", nThreads);
		if (nThreads < 1)
			throw std::runtime_error("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		//************************************ end reading input parameters from file

		//************************************ create vectors of input and output file names

		vector<string> vinfilenamesTot, vinfilenames1Tot; // input filenames for the defocused intensities or complex amplitudes
		vector< vector<string> > vvinfilenames(nangles), vvinfilenames1(nangles); // same input filenames for the defocused series in the form of vector of vectors
		vector<string> vinfilenamesTotSel; // temporary storage for input filenames before selecting a subseries
		if (imode3DPotential) // only read in pre-caculated 3D potential from "output" files and reprocess them
		{
			printf("\nInput file name base for pre-existing 3D potential = %s", filenamebaseOut.c_str());
			FileNames(1, nz, filenamebaseOut, vinfilenamesTot); // create 1D array of input filenames to read the input 2D slices of a previously reconstructed 3D object
			vvinfilenames.resize(nz);
			for (index_t na = 0; na < nz; na++)
			{
				vvinfilenames[na].resize(1); // number of defocus planes at the current illumination angle
				vvinfilenames[na][0] = vinfilenamesTot[na];
			}
		}
		else // read in the input files as directed
		{
			printf("\nInput defocus series file name base = %s", filenamebaseIn.c_str());
			if (bSelectFrames)
			{
				FileNames2(vndefocusFull, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames
				for (int i : sfIndexes)	vinfilenamesTotSel.push_back(vinfilenamesTot[i]); // create an array of selected frames in the same order as the corresponding defocus parameters
				vinfilenamesTot = vinfilenamesTotSel;
			}
			else
				FileNames2(vndefocus, filenamebaseIn, vinfilenamesTot); // create "total 2D array" of input filenames

			index_t ndefcurrent = 0;
			for (index_t na = 0; na < nangles; na++)
			{
				vvinfilenames[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				for (index_t n = 0; n < vndefocus[na]; n++) vvinfilenames[na][n] = vinfilenamesTot[ndefcurrent++];
			}
		}

		string filenamebaseOutNew("BAD_STRING"), filenamebaseOutCAmp("BAD_STRING"), filenamebaseOutDefocCAmpOrFM("BAD_STRING"); // don't use it, unless it is redefined later
		vector<string> voutfilenamesTot; // output filenames for the reconstructed 3D potential
		if (imode3DPotential) // output the the renormalized 3D potential
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "R");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the renormalized 3D potential = %s", filenamebaseOutNew.c_str());
			FileNames(1, nz, filenamebaseOutNew, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the filtered and/or rescaled 3D object
		}
		else
		{
			FileNames(1, nz, filenamebaseOut, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
		}

		vector<string> voutfilenamesPeaksTot; // output filenames for the peak-localized reconstructed 3D potential
		if (imodePeaks) // output filenames for the peak-localized 3D potential
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "A");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the peak-localized 3D potential = %s", filenamebaseOutNew.c_str());
			FileNames(1, nz, filenamebaseOutNew, voutfilenamesPeaksTot); // create 1D array of output filenames to save 2D slices of the peak-localized 3D object
		}

		vector<string> voutfilenamesTotCAmp; // output filenames the reprojected defocused complex amplitudes
		vector< vector<string> > vvoutfilenamesCAmp(nangles);  // same output filenames for the reprojected defocused complex amplitudes in the form of vector of vectors
		if (imodeReproj) // create filenames for saving the reprojected intensities, based on the names of input defocused intensities
		{
			std::filesystem::path apath(filenamebaseIn);
			filenamebaseOutCAmp = apath.filename().string();
			filenamebaseOutCAmp.replace(filenamebaseOutCAmp.find_last_of("."), filenamebaseOutCAmp.length() - filenamebaseOutCAmp.find_last_of("."), "R.grc");
			filenamebaseOutCAmp = folderAux + filenamebaseOutCAmp;
			printf("\nOutput file name base for reprojected complex amplitudes = %s", filenamebaseOutCAmp.c_str());
			FileNames2(vndefocus, filenamebaseOutCAmp, voutfilenamesTotCAmp); // create "2D array" of output filenames for reprojected complex amplitudes

			// Also copy the filenames of the input defocused intensities for the purpose of comparing with the reprojected intensities
			// note that in the case imode3DPotential != 0 these files are different from input images (which come from the potential)
			printf("\nDefocus series file name base for comparing with re-projected images = %s", filenamebaseIn.c_str());
			if (bSelectFrames)
			{
				FileNames2(vndefocusFull, filenamebaseIn, vinfilenames1Tot); // create "total 2D array" of input filenames
				for (int i : sfIndexes)	vinfilenamesTotSel.push_back(vinfilenames1Tot[i]); // create an array of selected frames in the same order as the corresponding defocus parameters
				vinfilenames1Tot = vinfilenamesTotSel;
			}
			else
				FileNames2(vndefocus, filenamebaseIn, vinfilenames1Tot); // create "total 2D array" of input filenames

			// Create the vector of vectors of output filenames for reprojected complex amplitudes at different defocus distances at each illumination angle
			// and the vector of vectors of filenames of the input defocused images for the purpose of comparing with the reprojected intensities
			index_t ndefcurrent = 0;
			for (index_t na = 0; na < nangles; na++)
			{
				vvoutfilenamesCAmp[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				vvinfilenames1[na].resize(vndefocus[na]); // number of defocus planes at the current illumination angle
				for (index_t n = 0; n < vndefocus[na]; n++)
				{
					vvoutfilenamesCAmp[na][n] = voutfilenamesTotCAmp[ndefcurrent];
					vvinfilenames1[na][n] = vinfilenames1Tot[ndefcurrent];
					ndefcurrent++;
				}
			}
		}

		vector<string> voutfilenamesTotDefocCAmpOrFM; // output filenames the phase-retrieved defocused complex amplitudes or filter matrix
		if (nSaveDefocCAmpsOrFilterMatrix == 1) // create filenames for saving the phase-retrieved defocused complex amplitudes or filter matrix
		{
			if (nPhaseRetrieval == 5)
			{
				std::filesystem::path apath(filenamebaseOut);
				filenamebaseOutDefocCAmpOrFM = apath.filename().string();
				filenamebaseOutDefocCAmpOrFM.replace(filenamebaseOutDefocCAmpOrFM.find_last_of("."), filenamebaseOutDefocCAmpOrFM.length() - filenamebaseOutDefocCAmpOrFM.find_last_of("."), "FM.grd");
				filenamebaseOutDefocCAmpOrFM = folderAux + filenamebaseOutDefocCAmpOrFM;
				FileNames(nz, 1, filenamebaseOutDefocCAmpOrFM, voutfilenamesTotDefocCAmpOrFM);
			}
			else
			{
				std::filesystem::path apath(filenamebaseIn);
				filenamebaseOutDefocCAmpOrFM = apath.filename().string();
				filenamebaseOutDefocCAmpOrFM.replace(filenamebaseOutDefocCAmpOrFM.find_last_of("."), filenamebaseOutDefocCAmpOrFM.length() - filenamebaseOutDefocCAmpOrFM.find_last_of("."), "D.grc");
				filenamebaseOutDefocCAmpOrFM = folderAux + filenamebaseOutDefocCAmpOrFM;
				FileNames(nangles, 1, filenamebaseOutDefocCAmpOrFM, voutfilenamesTotDefocCAmpOrFM);
			}
			printf("\nOutput file name base for phase-retrieved defocused complex amplitudes or filter matrix = %s", filenamebaseOutDefocCAmpOrFM.c_str());
		}

		//************************************ end creating vectors of input and output file names


		//*********************************** start main calculations
		bool bAbort(false);

		std::unique_ptr<IXAHead> pHead(new Wavehead2D(wl, ylo, yhi, xlo, xhi)); // any header data in input image files will be ignored
		XArray3D<double> K3out; // big 3D reconstructed array (needs to fit into RAM alongside with with everything else)

		printf("\n\nPreparing FFTW library ...");
		// single Fftwd2c object is used for repeated 2D FFTs using FFTW later (the plan is created at this point and reused later)
		Fftwd2c fftw2D((int)ny, (int)nx, true);

#ifndef CORRELATE_FRAMES

		if (!imode3DPotential) // do phase retrieval and backpropagation prior to 3D filtering or reprojection and output
		{
			printf("\nPerforming phase retrieval and backpropagation into the 3D volume ...");

			index_t naSym{ 0 }, naSymCurrent{ 0 }; // indexes of illumination angles with respect to subdivision into symmetry subsets (only used when imodeSymmetry == 1)
			XArray3D<double> K3outsum; // accumulated K3out (only used when imodeSymmetry == 1)
			XArray3D<dcomplex> V3; // 3D potential in the Fourier space (only used in the case of nPhaseRetrieval == 5)
			XArray3D<double> Filt3; // universal 3D filter matrix in the Fourier space (only used in the case of nPhaseRetrieval == 5)	

			// allocate the large 3D output arrays
			if (nPhaseRetrieval == 5)
			{
				V3.Resize(nz, ny, nx);
				V3.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
				Filt3.Resize(nz, ny, nx);
				Filt3.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
			}
			else
			{
				K3out.Resize(nz, ny, nx);
				K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));

				if (imodeSymmetry == 1) // !!!@@@@@ the implementation of this mode must be redone, following the change in the omp mode
				{
					K3outsum.Resize(nz, ny, nx);
					K3outsum.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
				}
			}

			// start of cycle over illumination angles
			//#pragma omp parallel for //!!! @@@ for some reason, this parallelization does not work well. Try the MsctKirkland model instead.
			for (int na = 0; na < nangles; na++)
			{
				if (bAbort) continue;
				try
				{
					double angleZ = v2angles[na].a * PI180;
					double angleY = v2angles[na].b * PI180;
					double cosangleY = cos(angleY);
					double sinangleY = sin(angleY);
					double cosangleZ = cos(angleZ);
					double sinangleZ = sin(angleZ);
					if (bVerboseOutput) printf("\n\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
					else printf("\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);

					// the following three centre of rotation expressions (linear phase factors) are used only when nPhaseRetrieval == 5
					dcomplex Fxc, Fyc, Fzc;
					if (bRotCentreShift)
					{
						Fxc = dcomplex(0, -1) * tPI * (dxc * (1.0 - cosangleY * cosangleZ) - dyc * sinangleZ + dzc * sinangleY * cosangleZ);
						Fyc = dcomplex(0, -1) * tPI * (dxc * cosangleY * sinangleZ + dyc * (1.0 - cosangleZ) - dzc * sinangleY * sinangleZ);
						Fzc = dcomplex(0, -1) * tPI * (-dxc * sinangleY + dzc * (1.0 - cosangleY));
					}

					// start the execution timer per angle
					std::chrono::system_clock::time_point start_timeA = std::chrono::system_clock::now();

					index_t ndefocus = vndefocus[na]; // number of defocus planes at the current illumination angle
					vector<Pair> vdefocus = vvdefocus[na]; // vector of input Z" rotation angles and defocus positions at the current illumination angle
					vector<string> vinfilenames = vvinfilenames[na]; // input filenames of defocused images at the current illumination angle
					vector<string> voutfilenames(nz);

					double zout(0.0); // z plane at which campOut will be defined
					XArray2D<double> int0; // input defocused intensity image (only used if bCTF_CT == true)
					XArray2D<dcomplex> campOut; // complex amplitude in the "middle" defocus plane (used to store the input GRC data or the complex amplitude after phase retrieval)

					// do phase retrieval (if required)
					if (bGRCinput) // phase retrieval is not required - work with the provided complex amplitude
					{
						if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[0].a, vdefocus[0].b);
						zout = vdefocus[0].b;
						if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[0].c_str());
						XArData::ReadFileGRC(campOut, vinfilenames[0].c_str()); //	read the only one (!) input GRC file
						index_t ny0 = campOut.GetDim1();
						index_t nx0 = campOut.GetDim2();
						if (nx0 < nx || ny0 < ny)
							throw std::runtime_error("Error: the dimensions of the 2D array in the image file are smaller than the array dimensions in Width Height HeaderLength Endianness ElementLength line from input parameter file.");
						else if (nx0 > nx || ny0 > ny)
						{
							printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions (Width = %zd, Height = %zd) from input parameter file.", nx, ny);
							XArray2DMove<dcomplex> tmp2(campOut);
							tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
						}
						campOut.SetHeadPtr(pHead->Clone()); // we ignore the head read from file
						//double averInt = campOut.Norm(eNormL2N);
						//if (averInt == 0) throw std::runtime_error("Error: zero intensity values in the input file.");
						//if (averInt != 1.0) campOut /= averInt; // enforce unit average intensity

						if (vdefocus[0].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0)) // rotate around Z'', then shift along XY, if needed
						{
							XArray2D<double> ampRe{ Re(campOut) };
							XArray2D<double> ampIm{ Im(campOut) };
							double averRe = ampRe.NormAverEdge(5);
							double averIm = ampIm.NormAverEdge(5);
							// shift along X and/or Y back to the unshifted position
							if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
							{
								xar::XArray2DMove<double> xamoveRe(ampRe);
								xamoveRe.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), averRe);
								xar::XArray2DMove<double> xamoveIm(ampIm);
								xamoveIm.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), averIm);
							}
							// rotate input defocused complex amplitude around Z'' back to zero angle
							if (vdefocus[0].a != 0)
							{
								XArray2DSpln<double> xaSplnRe(ampRe), xaSplnIm(ampIm);
								ampRe = xaSplnRe.Rotate(-vdefocus[0].a, 0.5 * (yhi + ylo) + dyc, 0.5 * (xhi + xlo) + dxc, averRe); // expecting uniform background
								ampIm = xaSplnIm.Rotate(-vdefocus[0].a, 0.5 * (yhi + ylo) + dyc, 0.5 * (xhi + xlo) + dxc, averIm); // expecting uniform background
							}
							campOut = MakeComplex(ampRe, ampIm, false);
						}

						// optionally save rotated defocused and shifted complex amplitudes in GRC files
						if (nSaveDefocCAmpsOrFilterMatrix == 1)
							XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmpOrFM[na].c_str(), eGRCBIN);
					}
					else // do phase retrieval and find the complex wave amplitude in the plane z = zout (if required by the phase retrieval method)
					{
						// read defocused images from files
						vector<XArray2D<double>> vint0(ndefocus); // initial defocused intensities
						for (int n = 0; n < ndefocus; n++)
						{
							start_timeNow = std::chrono::system_clock::now();
							if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[n].a, vdefocus[n].b);
							if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[n].c_str());
							if (bRAWinput) XArData::ReadFileRAW(vint0[n], vinfilenames[n].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
							else if (bTIFFinput) TIFFReadFile(vint0[n], vinfilenames[n].c_str()); // read input TIFF files
							else XArData::ReadFileGRD(vint0[n], vinfilenames[n].c_str(), wl); //	read input GRD files
							liFileReadTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();

							index_t ny0 = vint0[n].GetDim1();
							index_t nx0 = vint0[n].GetDim2();
							if (nx0 < nx || ny0 < ny)
								throw std::runtime_error("Error: the dimensions of the 2D array in the image file are smaller than the array dimensions in Width Height HeaderLength Endianness ElementLength line from input parameter file.");
							else if (nx0 > nx || ny0 > ny)
							{
								printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions (Width = %zd, Height = %zd) from input parameter file", nx, ny);
								XArray2DMove<double> tmp2(vint0[n]);
								tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
							}
							vint0[n].SetHeadPtr(pHead->Clone());

							// renormalize input data
							if (dNormFactor1 != 1.0) vint0[n] /= dNormFactor1;
							if (dNormFactor2 != 0.0) vint0[n] += dNormFactor2;
							if (vint0[n].Norm(eNormMin) < 0)
							{
								printf("\nWARNING: negative values in the input intensity file.");
								vint0[n].ThresholdLow(0.0, 0.0);
							}
							//double averInt = vint0[n].Norm(eNormAver);
							//if (averInt == 0) throw std::runtime_error("Error: zero average value in the input file after normalization.");
							//if (averInt != 1.0) vint0[n] /= averInt; // enforce unit average intensity

							// rotate input defocused image around Z'' back to zero angle
							if (vdefocus[n].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
							{
								double aver = vint0[n].NormAverEdge(5);
								// shift along X and/or Y back to the unshifted position
								if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
								{
									XArray2DMove<double> xamove(vint0[n]);
									xamove.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), aver);
								}
								// rotate input defocused complex amplitude around Z'' back to zero angle
								if (vdefocus[n].a != 0)
								{
									XArray2D<double> vintTemp(vint0[n]);
									XArray2DSpln<double> xaSpln(vintTemp);
									xaSpln.Rotate(vint0[n], -vdefocus[n].a, 0.5 * (yhi + ylo) + dyc, 0.5 * (xhi + xlo) + dxc, aver); // expecting uniform background
								}
							}
						}

						// do phase retrieval
						switch (nPhaseRetrieval)
						{
						case 1: // IWFR phase retrieval in the central image plane (to be followed by back-propagation of the resultant complex amplitude)
						{
							vector<double> vdefocusdist(ndefocus);
							for (int n = 0; n < ndefocus; n++) vdefocusdist[n] = vdefocus[n].b;
							if (bVerboseOutput) printf("\nPerforming IWFR phase retrieval ...");
							zout = 0.0;
							for (size_t n = 0; n < ndefocus; n++) zout += vdefocus[n].b;
							zout /= double(ndefocus);
							XA_IWFR<double> xa_iwfr;
							xa_iwfr.Iwfr(vint0, campOut, vdefocusdist, zout, k2maxo, Cs3, Cs5, itermax, epsilon, bVerboseOutput);
							if (nSaveDefocCAmpsOrFilterMatrix == 1) XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmpOrFM[na].c_str(), eGRCBIN);
						}
						break;
						case 2: // CTFL2 phase retrieval in the central image plane (to be followed by back-propagation of the resultant complex amplitude)
						{
							vector<double> vdefocusdist(ndefocus);
							for (int n = 0; n < ndefocus; n++) vdefocusdist[n] = vdefocus[n].b;
							if (bVerboseOutput) printf("\nPerforming CTFL2 phase retrieval ...");
							zout = vdefocus[0].b; //CTFL2 always retrieves the phase in the first defocused plane
							XA_IWFR<double> xa_iwfr;
							XArray2D<double> fiOut, ampOut;
							ampOut = vint0[0]; // use this before the call to xa_iwfr.CTFL2 spoils all vint0 images
							ampOut ^= 0.5; // preserve the amplitude at the zeros distance for reuse after phase retrieval
							fiOut = xa_iwfr.CTFL2(vint0, vdefocusdist, k2maxo, Cs3, Cs5, epsilon);
							MakeComplex(ampOut, fiOut, campOut, true);
							if (nSaveDefocCAmpsOrFilterMatrix == 1) XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmpOrFM[na].c_str(), eGRCBIN);
						}
						break;
						case 3: // Min05LogAmp phase retrieval in the given image plane (to be followed by back-propagation of the resultant complex amplitude)
						{
							zout = vdefocus[0].b;
							double zoutSample = zout + 0.5 * (zhidz - zlodz); // we want to take into account the sample thickness here
							double sigma = (dAtomWidth / 2.355); // std of the Gaussian fit for the potential of one atom
							double NF = zoutSample != 0 ? tPI * sigma * sigma / (wl * zoutSample) : 1.0;
							XArray2D<double> fiOut, ampOut;
							ampOut = vint0[0]; ampOut ^= 0.5;
							XA_IWFR<double> xa_iwfr;
							if (bVerboseOutput)
							{
								//printf("\nPerforming Fact * {-0.5 * NF * (Intensity / I0 - 1)} phase retrieval ...");
								printf("\nPerforming Fact * {-0.5 * (Intensity / I0 - 1)} phase retrieval ...");
								printf("\nFresnel number = %g", NF);
							}
							//fiOut = xa_iwfr.Min05LogAmp(vint0[0], epsilon * NF);
							fiOut = xa_iwfr.Min05LogAmp(vint0[0], epsilon);
							MakeComplex(ampOut, fiOut, campOut, true);
							if (nSaveDefocCAmpsOrFilterMatrix == 1) XArData::WriteFileGRC(campOut, voutfilenamesTotDefocCAmpOrFM[na].c_str(), eGRCBIN);
						}
						break;
						case 4: // vCTF phase retrieval will be done directly into the object volume below
						case 5: // sCTF phase retrieval will be done directly onto the Ewald sphere below
						{
							int0 = vint0[0];
							zout = vdefocus[0].b;
						}
						break;
						}
					}

					// now start calculations of the output defocused images
					if (bVerboseOutput) printf("\nPropagating to output defocus planes ...");
				
					XArray2D<dcomplex> K2four; // symmetrized backpropagated contrast in the reciprocal space (only used in the case of nPhaseRetrieval == 5)
					XArray3D<double> K3defoc; // defocused contrast (this array is used in the case of nPhaseRetrieval != 5)
					if (nPhaseRetrieval != 5) K3defoc.Resize(nz, ny, nx); // !!! This large 3D array is created separately for each na, which has implications for multithreading

					if (nPhaseRetrieval == 1 || nPhaseRetrieval == 2 || nPhaseRetrieval == 3) // for these phase retrieval methods, do Fresnel backpropagation into the reconstructed volume
					{
						int iSign = -1; // we want to compensate for (i.e. revert) the known aberrations during the Fresnel backpropagation
						xar::XArray2DFFT<double> xafft(campOut);
						xafft.FFT(true, false); // forward FFT once, before propagating to different output planes (to avoid repeated FFTs and speed up calculations) 

						if (bESCC) // take Ewald sphere curvature into account and do backpropagation
						{
							#pragma omp parallel for shared(K3defoc, campOut)
							for (int n = 0; n < nz; n++) // cycle over backpropagation to different planes
							{
								if (bAbort) continue;
								XArray2D<dcomplex> campTemp(campOut); // we need to renew this at each n
								try
								{
									// back-propagate to the output defocused plane
									xar::XArray2DFFT<double> xafft1(campTemp, true);

									if (bRelion)
										xafft1.FresnelA(iSign * (zout - voutdefocus[n]), iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to z_out[n]
									else
										xafft1.Fresnel(iSign * (zout - voutdefocus[n]), true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to z_out[n]
								}
								catch (std::exception& E)
								{
									printf("\n\n!!!Exception: %s\n", E.what());
									bAbort = true;
								}
								for (int j = 0; j < ny; j++)
									for (int i = 0; i < nx; i++)
										K3defoc[n][j][i] = 1.0 - std::norm(campTemp[j][i]);
							}
							if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
						}
						else // assume flat Ewald sphere and do backprojection
						{
							// back-propagate to the central plane
							if (bRelion)
								xafft.FresnelA(iSign * (zout - voutdefocus[nz / 2]), iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to z_out[n]
							else
								xafft.Fresnel(iSign * (zout - voutdefocus[nz /2]), true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to z_out[n]

							// now do backprojection (copy the same contrast to all planes)
							#pragma omp parallel for shared(K3defoc, campOut)
							for (int n = 0; n < nz; n++)
								for (int j = 0; j < ny; j++)
									for (int i = 0; i < nx; i++)
											K3defoc[n][j][i] = 1.0 - std::norm(campOut[j][i]);
						}
					}
					else // for these phase retrieval methods, do variable or constant CTF correction
					{
						int iSign = 1; // all aberrations are compensated here by means of division by the relevant forward-propagation aberrated CTF
						if (nPhaseRetrieval == 4) // variable-distance CTF retrieval of the phase at different planes inside the reconstruction volume
						{
							if (bESCC) // take Ewald sphere curvature into account and do variable-distance CTF correction
							{
								#pragma omp parallel for shared(K3defoc, campOut)
								for (int n = 0; n < nz; n++)
								{
									if (bAbort) continue;
									XArray2D<double> int0temp;
									int0temp = int0; // we need to renew this for each n, as InversePhaseCTF spoils its first argument
									try
									{
										// propagate to the output defocused plane, i.e. do variable CTF correction (CTF backpropagation) of the projected phase
										XA_IWFR<double> xa_iwfr;
										if (bRelion)
											xa_iwfr.InvertPhaseCTF(int0temp, iSign * (voutdefocus[n] - zout), iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, k2maxo, iSign * Cs3, iSign * Cs5, asigma, epsilon);
										else
											xa_iwfr.InvertPhaseCTF(int0temp, iSign * (voutdefocus[n] - zout), 0, 0, k2maxo, iSign * Cs3, iSign * Cs5, asigma, epsilon);
									}
									catch (std::exception& E)
									{
										printf("\n\n!!!Exception: %s\n", E.what());
										bAbort = true;
									}
									for (int j = 0; j < ny; j++)
										for (int i = 0; i < nx; i++)
											K3defoc[n][j][i] = int0temp[j][i];
								}
								if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
							}
							else // assume flat Ewald sphere and do conventional CTF correction
							{
								XA_IWFR<double> xa_iwfr;
								if (bRelion)
									xa_iwfr.InvertPhaseCTF(int0, (voutdefocus[nz / 2] - zout), iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, k2maxo, iSign * Cs3, iSign * Cs5, asigma, epsilon);
								else
									xa_iwfr.InvertPhaseCTF(int0, (voutdefocus[nz / 2] - zout), 0, 0, k2maxo, iSign * Cs3, iSign * Cs5, asigma, epsilon);

								// now do backprojection (copy the same contrast to all planes)
								#pragma omp parallel for shared(K3defoc, int0)
								for (int j = 0; j < ny; j++)
									for (int i = 0; i < nx; i++)
										for (int n = 0; n < nz; n++)
											K3defoc[n][j][i] = int0[j][i];
							}
						}
						if (nPhaseRetrieval == 5) // symmetrized DT reconstruction on the Ewald sphere in the reciprocal space
						{
							int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0
							//int0 -= 1.0; // convert image intensity into -K
							start_timeNow = std::chrono::system_clock::now();
							XA_IWFR<double> xa_iwfr;
							if (bRelion)
								xa_iwfr.InvertCTF_DT1(int0, K2four, fftw2D, zout, iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, k2maxo, iSign * Cs3, iSign * Cs5, asigma, epsilon, bESCC, nThreads);
							else
								xa_iwfr.InvertCTF_DT1(int0, K2four, fftw2D,  zout, 0, 0, k2maxo, iSign * Cs3, iSign * Cs5, asigma, epsilon, bESCC, nThreads);
							li2DFFTtime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
						}
					}

					if (bSelectFrames) 
						if (nPhaseRetrieval == 5) K2four *= sfWeights[na];
						else K3defoc *= sfWeights[na];

					// update the 3D reconstructed object at the current illumination direction
					if (bVerboseOutput) printf("\nUpdating 3D reconstructed object ...");
					double xc = (xhi + xlo) / 2.0 + dxc; // x-coordinate of the centre of rotation
					double yc = (yhi + ylo) / 2.0 + dyc; // y-coordinate of the centre of rotation
					double zc = (zhidz + zlodz) / 2.0 + dzc; // z-coordinate of the centre of rotation

					start_timeNow = std::chrono::system_clock::now();
					if (nPhaseRetrieval != 5) // update the 3D potential in real space with the current back-propagated contrast dirstibution K3defoc
					{
						// calculate the coordinate illumination angle parameters
						double xxx;
						vector<double> x_sinangleY(nx), x_cosangleY(nx);
						for (int i = 0; i < nx; i++)
						{
							xxx = xlo + xst * i - xc;
							x_sinangleY[i] = xxx * sinangleY;
							x_cosangleY[i] = xxx * cosangleY;
						}

						double yyy;
						vector<double> y_sinangleZ(ny), y_cosangleZ(ny);
						for (int j = 0; j < ny; j++)
						{
							yyy = ylo + yst * j - yc;
							y_sinangleZ[j] = yyy * sinangleZ;
							y_cosangleZ[j] = yyy * cosangleZ;
						}
						double zzz;
						vector<double> z_sinangleY(nz), z_cosangleY(nz);
						for (int n = 0; n < nz; n++)
						{
							zzz = voutdefocus[n] - zc;
							z_sinangleY[n] = zzz * sinangleY;
							z_cosangleY[n] = zzz * cosangleY;
						}

						int ii, jj, nn;
						double xx;
						int nx2 = int(nx) - 2, ny2 = int(ny) - 2, nz2 = int(nz) - 2;
						double xx_sinangleZ, xx_cosangleZ, dx0, dx1, dy0, dy1, dz0, dz1, dz0K, dz1K;

						#pragma omp parallel for shared (K3out, K3defoc) private(xx, xxx, yyy, zzz, xx_sinangleZ, xx_cosangleZ, dx1, dx0, dy1, dy0, dz1, dz0, ii, jj, nn, dz0K, dz1K)
						for (int n = 0; n < nz; n++)
						{
							try
							{
								for (int i = 0; i < nx; i++)
								{
									if (bAbort) continue;
									// inverse rotation around Y' axis
									zzz = zc + z_cosangleY[n] + x_sinangleY[i]; // z coordinate after the inverse rotation around Y'
									dz1 = (zzz - zlodz) / zst;
									nn = (int)dz1; if (nn < 0 || nn > nz2) continue;
									dz1 -= nn; dz0 = 1.0 - dz1;

									xx = -z_sinangleY[n] + x_cosangleY[i]; // x - xc coordinate after the inverse rotation around Y'
									xx_sinangleZ = xx * sinangleZ;
									xx_cosangleZ = xx * cosangleZ;

									for (int j = 0; j < ny; j++)
									{
										// inverse rotation around Z axis
										xxx = xc + y_sinangleZ[j] + xx_cosangleZ; // x coordinate after the inverse rotation around Z
										dx1 = (xxx - xlo) / xst;
										ii = (int)dx1; if (ii < 0 || ii > nx2) continue;
										dx1 -= ii; dx0 = 1.0 - dx1;

										yyy = yc + y_cosangleZ[j] - xx_sinangleZ; // y coordinate after the inverse rotation around Z
										dy1 = (yyy - ylo) / yst;
										jj = (int)dy1; if (jj < 0 || jj > ny2) continue;
										dy1 -= jj; dy0 = 1.0 - dy1;

										dz0K = dz0 * K3defoc[n][j][i];
										dz1K = dz1 * K3defoc[n][j][i];
										K3out[nn][jj][ii] += dx0 * dy0 * dz0K;
										K3out[nn][jj][ii + 1] += dx1 * dy0 * dz0K;
										K3out[nn + 1][jj][ii] += dx0 * dy0 * dz1K;
										K3out[nn + 1][jj][ii + 1] += dx1 * dy0 * dz1K;
										K3out[nn][jj + 1][ii] += dx0 * dy1 * dz0K;
										K3out[nn][jj + 1][ii + 1] += dx1 * dy1 * dz0K;
										K3out[nn + 1][jj + 1][ii] += dx0 * dy1 * dz1K;
										K3out[nn + 1][jj + 1][ii + 1] += dx1 * dy1 * dz1K;
									}
								}
							}
							catch (std::exception& E)
							{
								printf("\n\n!!!Exception: %s\n", E.what());
								bAbort = true;
							}
						}
						if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

						// rotate the reconstructed volume if needed to distribute the illumination directions according to the symmetry group
						if (imodeSymmetry == 1)
						{
							naSym = int((na + 1) / nanglesIn1Sym); // integer division rounded to the closest lower integer; switches to the next one in steps of nanglesIn1Sym
							if (naSym > nanglesSym - 1) naSym -= nanglesSym * int(naSym / nanglesSym);
							if (na == nangles - 1) // last illumination angle, finalize this process
							{
								printf("\nFinished processing the remainder of all illumination angles, na = %d, symmetry_angle_no = %zd", na, naSymCurrent);
								printf("\nRotating around Z by %g (deg), around Y' by %g (deg) and around Z'' by %g (deg) ...", -v2anglesSym[naSymCurrent].a, -v2anglesSym[naSymCurrent].b, -vvdefocusSym[naSymCurrent][0].a);
								XArray3DSpln<double> xa3spln(K3out);
								XArray3D<double> K3outRot;
								xa3spln.Rotate3(K3outRot, -vvdefocusSym[naSymCurrent][0].a * PI180, -v2anglesSym[naSymCurrent].b * PI180, -v2anglesSym[naSymCurrent].a * PI180, nThreads); // !! rotate in the opposite direction
								if (K3outsum.size() == 0) K3out = K3outRot;
								else { K3outsum += K3outRot; K3out = K3outsum; } // put all the data back into K3out for further processing
							}
							else if (naSym != naSymCurrent) // new subset of illumination angles will be processed from the next iteration over na
							{
								printf("\nCompleted a group of illumination angles, na = %d, symmetry_angle_no = %zd", na, naSymCurrent);
								printf("\nRotating around Z by %g (deg), around Y' by %g (deg) and around Z'' by %g (deg) ...", -v2anglesSym[naSymCurrent].a, -v2anglesSym[naSymCurrent].b, -vvdefocusSym[naSymCurrent][0].a);
								XArray3DSpln<double> xa3spln(K3out);
								XArray3D<double> K3outRot;
								xa3spln.Rotate3(K3outRot, -vvdefocusSym[naSymCurrent][0].a * PI180, -v2anglesSym[naSymCurrent].b * PI180, -v2anglesSym[naSymCurrent].a * PI180, nThreads); // !! rotate in the opposite direction
								if (K3outsum.size() == 0) K3outsum = K3outRot; else K3outsum += K3outRot;
								K3out.Fill(0.0); // a separate accumulation file K3outsum is used here to avoid multiple rotations of the data in K3out during the accumulation over na
								naSymCurrent = naSym;
							}
						}
					}
					else // nPhaseRetrieval == 5 : update the reconstructed potential on the Ewald sphere with the current 2D CTF-corrected contrast in the reciprocal space
					{
						// calculate the coordinate illumination angle parameters
						int nxd2 = int(nx / 2);
						double fxst = 1.0 / (xhi - xlo), xxx, wl2 = 0.5 * wl;
						double fxlo = -fxst * nxd2;
						vector<double> xxx2(nx);
						vector<double> x_sinangleY(nx), x_cosangleY(nx);
						for (int i = 0; i < nx; i++)
						{
							xxx = fxlo + fxst * i;
							xxx2[i] = -wl2 * xxx * xxx;
							x_sinangleY[i] = xxx * sinangleY;
							x_cosangleY[i] = xxx * cosangleY;
						}

						int nyd2 = int(ny / 2);
						double fyst = 1.0 / (yhi - ylo), yyy;
						double fylo = -fyst * nyd2;
						vector<double> yyy2(ny);
						vector<double> y_sinangleZ(ny), y_cosangleZ(ny);
						for (int j = 0; j < ny; j++)
						{
							yyy = fylo + fyst * j;
							yyy2[j] = -wl2 * yyy * yyy;
							y_sinangleZ[j] = yyy * sinangleZ;
							y_cosangleZ[j] = yyy * cosangleZ;
						}

						int nzd2 = int(nz / 2);
						double fzst = 1.0 / (zhi - zlo);
						double fzlo = -fzst * nzd2;
						// qz rotated values are expressed via qx and qy - see below
				
						int nx2 = int(nx) - 2, ny2 = int(ny) - 2, nz2 = int(nz) - 2;

						index_t nxny = (nx * ny);
						//int ijnR2 = (int)pow(min(min(nxd2 - 1, nyd2 - 1), nzd2 - 1), 2); // square radius of the reconstruction sphere

						#pragma omp parallel for shared (V3, Filt3, K2four)
						for (int j = 0; j < ny; j++)
						{
							int ii, jj, nn;
							index_t nji;
							double dx0, dx1, dy0, dy1, dz0, dz1, dxyz, xx, xxx, yyy, zz, zzz;
							dcomplex cd;
							double* pFilt3 = &(Filt3[0][0][0]);
							dcomplex* pV3 = &(V3[0][0][0]);

							for (int i = 0; i < nx; i++)
							{
								// inverse rotation around Y' axis
								if (bESCC) zz = xxx2[i] + yyy2[j];
								else zz = 0.0; // conventional CT
								xx = x_cosangleY[i] - zz * sinangleY; // qx coordinate after the rotation around Y'
								zzz = x_sinangleY[i] + zz * cosangleY; // qz coordinate after the rotation around Y'
								dz1 = (zzz - fzlo) / fzst;
								nn = (int)floor(dz1); if (nn < 0 || nn > nz2) continue;
								dz1 -= nn; dz0 = 1.0 - dz1;

								// inverse rotation around Z axis
								xxx = xx * cosangleZ + y_sinangleZ[j]; // qx coordinate after the rotation around Z
								yyy = -xx * sinangleZ + y_cosangleZ[j]; // qy coordinate after the rotation around Z

								dx1 = (xxx - fxlo) / fxst;
								ii = (int)floor(dx1); if (ii < 0 || ii > nx2) continue;
								dx1 -= ii; dx0 = 1.0 - dx1;
								dy1 = (yyy - fylo) / fyst;
								jj = (int)floor(dy1); if (jj < 0 || jj > ny2) continue;
								dy1 -= jj; dy0 = 1.0 - dy1;

								// exclude points lying outside the reconstruction sphere
								//if ((ii - nxd2) * (ii - nxd2) + (jj - nyd2) * (jj - nyd2) + (nn - nzd2) * (nn - nzd2) >= ijnR2) continue;

								cd = K2four[j][i];
								// multiply by linear phase factors due to the shifted centre of rotation in the real space
								if (bRotCentreShift) cd *= exp(Fxc * xxx + Fyc * yyy + Fzc * zzz);
							
								dxyz = dx0 * dy0 * dz0;
								nji = nn * nxny + jj * nx + ii;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx1 * dy0 * dz0;
								nji += 1;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx0 * dy0 * dz1;
								nji += -1 + nxny;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx1 * dy0 * dz1;
								nji += 1;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx0 * dy1 * dz0;
								nji += -1 - nxny + nx;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx1 * dy1 * dz0;
								nji += 1;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx0 * dy1 * dz1;
								nji += -1 + nxny;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
								dxyz = dx1 * dy1 * dz1;
								nji += 1;
								*(pFilt3 + nji) += dxyz;
								*(pV3 + nji) += dxyz * cd;
							}
						}
						liBackPropTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
					} // end of cycle updating the 3D potential at the current rotational orientation "na"

					if (bVerboseOutput)
					{
						std::chrono::system_clock::time_point end_timeA = std::chrono::system_clock::now();
						printf("\nExecution time for this illumination angle = %zd ms.", std::chrono::duration_cast<std::chrono::milliseconds>(end_timeA - start_timeA).count());
					}

				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			} // end of cycle over illumination angles "na"
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

			// normalization of the reconstructed 3D distribution
			double dnorm;
			if (nPhaseRetrieval == 1 || nPhaseRetrieval == 2 || nPhaseRetrieval == 3)
			{
				dnorm = EE * dGamma / dAtomWidth;
				dnorm *= -wl / PI; // new CHR reconstruction formula without the inverse Laplacian

				if (bSelectFrames) dnorm /= dTotalWeight;
				else dnorm /= nangles; // 1/nangles plays the role of 1/(4*PI) in the theoretical formula 
			}
			if (nPhaseRetrieval == 4)
			{
				dnorm = wl * EE / PI / (dAtomWidth * 8.0); // conversion factor from Fi to V
				if (bSelectFrames) dnorm /= dTotalWeight;
				else dnorm /= nangles; // 1/nangles plays the role of 1/(4*PI) in the theoretical formula 
			}
			if (nPhaseRetrieval == 5)
			{
				printf("\n\nFiltering 3D distribution of the electrostatic potential in the Fourier space ...");
				double* pFilt3 = &(Filt3[0][0][0]);
				dcomplex* pV3 = &(V3[0][0][0]);
				double Filt3Min = Filt3.Norm(eNormMax), Filt3Max = Filt3.Norm(eNormMin);

				#pragma omp parallel for shared (V3, Filt3)
				for (int k = 0; k < nz; k++)
					for (index_t j = 0; j < ny; j++)
						for (index_t i = 0; i < nx; i++)
						{
							double x = Filt3[k][j][i];
							if (x == 0) continue;

							if (x < Filt3Min) Filt3Min = x;
							else if (x > Filt3Max) Filt3Max = x;

							x = (1 / x) * exp((x - CriticalSNR) / (x * x * x * x));

							// divide by the regulalized filter matrix
							Filt3[k][j][i] = x;
							V3[k][j][i] *= x;
							//V3[k][j][i] *= x / ( x * x + epsilon);
						}
				printf("\nFilt3-min = %f,  Filt3-max = %f", Filt3Min, Filt3Max);

				// optionally save the filter matrix in files
				if (nSaveDefocCAmpsOrFilterMatrix == 1)
				{
					printf("\nSaving the regularized inverse filter matrix in files ...");
					XArData::WriteFileStackGRD(Filt3, voutfilenamesTotDefocCAmpOrFM, eGRDBIN);
				}

				// delete Filt3 matrix to free memory prior to 3D FFT
				Filt3.Truncate();

				// enforce hermitian symmetry
				#pragma omp parallel for shared (V3)
				for (int k = 0; k < nz; k++)
				{
					index_t ka = (k == 0 ? 0 : nz - (index_t)k);
					for (index_t j = 0; j < ny; j++)
					{
						index_t ja = (j == 0 ? 0 : ny - j);
						for (index_t i = 0; i < nx / 2; i++)
						{
							index_t ia = (i == 0 ? 0 : nx - i);
							V3[k][j][ia] = 0.5 * (V3[k][j][ia] + std::conj(V3[ka][ja][i]));
						}
					}
				}

				// inverse FFT of the FFT(potential) on the Ewald sphere
				printf("\nAllocating Fftwd3drc object memory ...");
				// NOTE: it appears that multitheading actually slowing the FFTW 3D FFT down
				Fftwd3drc fft3((int)nz, (int)ny, (int)nx, 1, false); // create an Fftwd3drc object in the Fourier space state
				fft3.SetComplexXArray3D(V3, true);
				if (imodeInvLaplace)
				{
					printf("\nInverse Laplace filtering the reconstructed 3D object ...");
					fft3.InverseMLaplacian1(zhi - zlo, yhi - ylo, xhi - ylo, alpha, alpha);
				}
				V3.Truncate(); // free memory
				printf("\nInverse Fourier transforming the Fourier transform of the 3D distribution of the electrostatic potential ...");
				start_timeNow = std::chrono::system_clock::now();
				fft3.InverseFFT();
				li3DFFTtime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
				K3out.Resize(nz, ny, nx); // allocated memory
				K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
				fft3.GetRealXArray3D(K3out, true);
				fft3.Cleanup(); // free memory

				dnorm = wl * EE / (2.0 * PI * sqrt(1.0 + asigma * asigma)); // this is the norm for V = 2E delta
			}
			K3out *= dnorm;

		} // end of case if imode3DPotential == 0
		else // read in pre-existing 3D distribution of the electrostatic potential
		{
			printf("\n\nReading pre-existing 3D distribution of the electrostatic potential from files ...");

			K3out.Resize(nz, ny, nx);
			K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));

			// read input 2D slices
			XArray2D<double> ipIn;
			try
			{
				for (int n = 0; n < nz; n++)
				{
					if (bRAWinput) XArData::ReadFileRAW(ipIn, vinfilenamesTot[n].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
					else if (bTIFFinput) TIFFReadFile(ipIn, vinfilenamesTot[n].c_str()); // read input TIFF files
					else XArData::ReadFileGRD(ipIn, vinfilenamesTot[n].c_str(), wl); //	read input GRD files

					index_t ny0 = ipIn.GetDim1();
					index_t nx0 = ipIn.GetDim2();
					if (nx0 < nx || ny0 < ny)
						throw std::runtime_error("Error: the dimensions of the 2D array in the image file are smaller than the array dimensions in Width Height HeaderLength Endianness ElementLength line from input parameter file.");
					else if (nx0 > nx || ny0 > ny)
					{
						printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions Width = %zd Height = %zd from input parameter file", nx, ny);
						XArray2DMove<double> tmp2(ipIn);
						tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
					}
					ipIn.SetHeadPtr(pHead->Clone()); // we ignore the head read from file

					//if (dNormFactor1 != 1.0) ipIn /= dNormFactor1;
					//if (dNormFactor2 != 0.0) ipIn += dNormFactor2;

					for (index_t j = 0; j < ny; j++)
						for (index_t i = 0; i < nx; i++)
							K3out[n][j][i] = ipIn[j][i];
				}
			}
			catch (std::exception& E)
			{
				printf("\n\n!!!Exception: %s\n", E.what());
				bAbort = true;
			}
		}

		// enforce known symmetries
		if (imodeSymmetry == 2)
		{
			printf("\nEnforcing known symmetries of the reconstructed 3D potential ...");
			XArray3DSpln<double> xa3spln(K3out);
			XArray3D<double> K3outRot, K3outsum(K3out);
			for (int naSym = 1; naSym < nanglesSym; naSym++)
			{
				printf("\nRotating around Z by %g (deg), around Y' by %g (deg) and around Z'' by %g (deg) ...", -v2anglesSym[naSym].a, -v2anglesSym[naSym].b, -vvdefocusSym[naSym][0].a);
				xa3spln.Rotate3(K3outRot, -vvdefocusSym[naSym][0].a * PI180, -v2anglesSym[naSym].b * PI180, -v2anglesSym[naSym].a * PI180, nThreads); // !! rotate in the opposite direction
				K3outsum += K3outRot;
			}
			K3outsum /= double(nanglesSym);
			K3out = K3outsum;
		}

		// apply the regularized inverse 3D (-Laplacian)
		if (imodeInvLaplace && !(nPhaseRetrieval == 5 && !imode3DPotential)) // in the case nPhaseRetrieval == 5 and !imode3DPotential inverse Laplacian is calculated in the Fourier space
		{
			printf("\nInverse Laplace filtering the reconstructed 3D object ...");
			// NOTE: it appears that multitheading actually slowing the FFTW 3D FFT down
			Fftwd3drc fft3((int)nz, (int)ny, (int)nx, 1);
			fft3.InverseMLaplacian(K3out, alpha, alpha);
		}

		// high-pass filtering, background subtraction and thresholding of the reconstructed 3D distribution
		printf("\nHigh-pass filtering and thresholding the reconstructed 3D object ...");
		if (dlpfiltersize > 0 && dlpfiltersize < (xhi - xlo))
		{
			XArray3D<double> K3outTemp = K3out;
			// NOTE: it appears that multitheading actually slowing the FFTW 3D FFT down
			Fftwd3drc fft3((int)nz, (int)ny, (int)nx, 1);
			fft3.GaussFilter(K3outTemp, dlpfiltersize / 2.355);
			K3out -= K3outTemp;
		}
		K3out.ThresholdLow(dBackground, dThreshold);

		// output the 3D array
		printf("\n\nSaving the reconstructed 3D object into output files ...");
		start_timeNow = std::chrono::system_clock::now();
		if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesTot, eTIFF32);
		else XArData::WriteFileStackGRD(K3out, voutfilenamesTot, eGRDBIN);
		liFileWriteTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();

#endif // !CORRELATE_FRAMES

		// multislice reprojection
		if (imodeReproj)
		{
#ifndef CORRELATE_FRAMES
			printf("\nPerforming forward re-propagation through the 3D distribution of the electrostatic potential ...");
			double yc = (ylo + yhi) / 2.0 + dyc, xc = (xlo + xhi) / 2.0 + dxc;
			double dnorm3 = PI / (wl * EE) * zst;
			K3out *= dnorm3; // this converts the potential V into -2 * PI / lambda  * delta * zst, with zst required in the subsequent integration inside Multislice_eV
#endif
			xar::XArray3DSpln<double> spln3d(K3out);

			vector<vector<double> > vverr(nangles); // reconstruction errors in individual defocus planes at different illumination angles
			vector<double> verraver(nangles, 0.0); // reconstruction errors at different illumination angles averaged over the corresponding defocus distances
			vector<vector<double> > vvnorm0(nangles); // C0 norms of reprojected intensities-1 in individual defocus planes at different illumination angles
			vector<double> vnormaver0(nangles, 0.0); // C0 norms of reprojected intensities-1 different illumination angles averaged over the corresponding defocus distances
			vector<vector<double> > vvnorm1(nangles); // C0 norms of original intensities-1 in individual defocus planes at different illumination angles
			vector<double> vnormaver1(nangles, 0.0); // C0 norms of original intensities-1 different illumination angles averaged over the corresponding defocus distances

#ifdef CORRELATE_FRAMES
			string filenamebaseInRef(CORRELATE_FRAMES);
			printf("\nInput reference series file name base = %s", filenamebaseInRef.c_str());
			vector<string> vinfilenamesRef;
			FileNames2(vndefocus, filenamebaseInRef, vinfilenamesRef); // create "total 2D array" of input filenames
			XArray1D<double> vCorr(vndefocus.size());
			vCorr.SetHeadPtr(new Wavehead1D(wl, 0, (double)vndefocus.size()));
			printf("\n\n Correlating frames ...\n");
#endif //CORRELATE_FRAMES


#pragma omp parallel for shared (K3out, spln3d, vndefocus, vvdefocus, v2angles, vvoutfilenamesCAmp)
			for (int na = 0; na < nangles; na++)
			{
				if (bAbort) continue;
				try
				{
					vector<Pair> vdefocus;
					vector<string> voutfilenamesCAmp, vinfilenames1;
					double angleY(0);
					double angleZ(0);
#pragma omp critical
					{
						vdefocus = vvdefocus[na];
						vverr[na].resize(vndefocus[na]);
						vvnorm0[na].resize(vndefocus[na]);
						vvnorm1[na].resize(vndefocus[na]);
						voutfilenamesCAmp = vvoutfilenamesCAmp[na];
						vinfilenames1 = vvinfilenames1[na];
						angleZ = v2angles[na].a * PI180;
						angleY = v2angles[na].b * PI180;
						//if (!bVerboseOutput) printf("\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
					}
#ifndef CORRELATE_FRAMES
					// multislice propagate through the reconstructed 3D distribution of the electrostatic potential
					XArray2D<dcomplex> campOut(ny, nx, dcomplex(1.0, 0.0)); // "incident" complex amplitude
					campOut.SetHeadPtr(pHead->Clone());
					spln3d.Multislice_eV(campOut, angleZ, angleY, sliceTh, k2maxo);
#endif // !CORRELATE_FRAMES
					// propagate to defocuse planes, rotate around the illumination axis and save the defocused intensities
					XArray2D<double > vint0n, K0n, vint1n, K1n, ampRe, ampIm;
					for (int n = 0; n < vndefocus[na]; n++)
					{
#ifndef CORRELATE_FRAMES
						XArray2D<dcomplex> campOut1(campOut);
						xar::XArray2DFFT<double> xafft(campOut1);
						//int iSign = vdefocus[n].b < 0 ? iSign = -1 : iSign = 1;
						int iSign = +1; // now we want to reproduce(!) the known aberrations, rather than compensate them
						if (bRelion)
							xafft.FresnelA(vdefocus[n].b, iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to the current defocus distance
						else
							xafft.Fresnel(vdefocus[n].b, true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to the current defocus distance
						if (vdefocus[n].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0)) // rotate around Z'', then shift along XY, if needed
						{
							Re(campOut1, ampRe); Im(campOut1, ampIm);
							double averRe = ampRe.NormAverEdge(5);
							double averIm = ampIm.NormAverEdge(5);
							// rotate input defocused complex amplitude around Z''
							if (vdefocus[0].a != 0)
							{
								XArray2DSpln<double> xaSplnRe(ampRe), xaSplnIm(ampIm);
								ampRe = xaSplnRe.Rotate(vdefocus[n].a, yc, xc, averRe); // expecting uniform background
								ampIm = xaSplnIm.Rotate(vdefocus[n].a, yc, xc, averIm); // expecting uniform background
							}
							// shift along X and/or Y
							if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
							{
								xar::XArray2DMove<double> xamoveRe(ampRe), xamoveIm(ampIm);
								xamoveRe.Move((long)floor(v2shifts[na].b / yst + 0.5), (long)floor(v2shifts[na].a / xst + 0.5), averRe);
								xamoveIm.Move((long)floor(v2shifts[na].b / yst + 0.5), (long)floor(v2shifts[na].a / xst + 0.5), averIm);
							}
							MakeComplex(ampRe, ampIm, campOut1, false);
						}
						Abs(campOut1, vint1n);
						Abs2(campOut1, K1n); K1n -= 1.0;
#endif // !CORRELATE_FRAMES
#ifdef CORRELATE_FRAMES
						// read in the reference image
						XArData::ReadFileGRD(vint1n, vinfilenamesRef[na].c_str(), wl);
						double dmeanRef = vint1n.Norm(eNormAver);
						vint1n -= dmeanRef;
						double stdRef = vint1n.Norm(eNormStdDev);
						if (stdRef != 0) vint1n /= stdRef;
						double dL2_1 = vint1n.Norm(eNormL2);
						pHead.reset(vint1n.GetHeadPtr()->Clone());
#endif //CORRELATE_FRAMES
						// read in the original intensity
						if (GetFileExtension(vinfilenames1[n]) == string(".GRC"))
						{
							XArray2D<dcomplex> campOut2;
							XArData::ReadFileGRC(campOut2, vinfilenames1[n].c_str());
							index_t ny0 = campOut2.GetDim1();
							index_t nx0 = campOut2.GetDim2();
							if (nx0 < nx || ny0 < ny)
								throw std::runtime_error("Error: the dimensions of the 2D array in the image file are smaller than the array dimensions in Width Height HeaderLength Endianness ElementLength line from input parameter file.");
							else if (nx0 > nx || ny0 > ny)
							{
								printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions Width = %zd Height = %zd from input parameter file", nx, ny);
								XArray2DMove<dcomplex> tmp2(campOut2);
								tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
							}
							campOut2.SetHeadPtr(pHead->Clone()); // we ignore the head read from file
							//campOut2 /= campOut2.Norm(eNormL2N); // enforce unit average intensity
							Abs(campOut2, vint0n);
							Abs2(campOut2, K0n); K0n -= 1.0;
						}
						else
						{
							if (GetFileExtension(vinfilenames1[n]) == string(".RAW"))
								XArData::ReadFileRAW(vint0n, vinfilenames1[n].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
							else if (GetFileExtension(vinfilenames1[n]) == string(".TIF") || GetFileExtension(vinfilenames1[n]) == string(".TIFF"))
								TIFFReadFile(vint0n, vinfilenames1[n].c_str());
							else
								XArData::ReadFileGRD(vint0n, vinfilenames1[n].c_str(), wl);
							index_t ny0 = vint0n.GetDim1();
							index_t nx0 = vint0n.GetDim2();
							if (nx0 < nx || ny0 < ny)
								throw std::runtime_error("Error: the dimensions of the 2D array in the image file are smaller than the array dimensions in Width Height HeaderLength Endianness ElementLength line from input parameter file.");
							else if (nx0 > nx || ny0 > ny)
							{
								printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions Width = %zd Height = %zd from input parameter file", nx, ny);
								XArray2DMove<double> tmp2(vint0n);
								tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
							}
							vint0n.SetHeadPtr(pHead->Clone()); // we ignore the head read from file

							if (dNormFactor1 != 1.0) vint0n /= dNormFactor1;
							if (dNormFactor2 != 0.0) vint0n += dNormFactor2;
							if (vint0n.Norm(eNormMin) < 0)
							{
								//printf("\nWARNING: negative values in the input intensity file.");
								vint0n.ThresholdLow(0.0, 0.0);
							}
							//double averInt = vint0n.Norm(eNormAver);
							//if (averInt == 0) throw std::runtime_error("Error: zero average value in the input file after normalization.");
							//if (averInt != 1.0) vint0n /= averInt; // enforce unit average intensity
#ifndef CORRELATE_FRAMES
							K0n = vint0n; K0n -= 1.0;
							vint0n ^= 0.5;
#endif // !CORRELATE_FRAMES
#ifdef CORRELATE_FRAMES
							double dmeanFr = vint0n.Norm(eNormAver);
							vint0n -= dmeanFr;
							double stdFr = vint0n.Norm(eNormStdDev);
							if (stdFr != 0) vint0n /= stdFr;
							double dL2_0 = vint0n.Norm(eNormL2);
							vint0n *= vint1n;
							vCorr[na] = vint0n.Norm(eNormAver) * (double)vint0n.size() / (dL2_0 * dL2_1);
							printf("\rFrame no. %d, correlation coefficient = %g", na, vCorr[na]);
							fflush(stdout);
#endif // CORRELATE_FRAMES
						}
						// compare the reprojected and the original defocused intensities
#ifndef CORRELATE_FRAMES
						vvnorm0[na][n] = K0n.Norm(eNormC0);
						vnormaver0[na] += vvnorm0[na][n];
						vvnorm1[na][n] = K1n.Norm(eNormC0);
						vnormaver1[na] += vvnorm1[na][n];
						vint1n -= vint0n;
						vverr[na][n] = pow(vint1n.Norm(eNormL2), 2.0) / pow(vint0n.Norm(eNormL2), 2.0); // this error norm is modelled after the one in IWFR
						verraver[na] += vverr[na][n];

#pragma omp critical
						{
							if (bVerboseOutput)
							{
								printf("\n\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
								printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (Angstroms)", vdefocus[n].a, vdefocus[n].b);
								printf("\nWriting output file %s ...", voutfilenamesCAmp[n].c_str());
							}
						}
						XArData::WriteFileGRC(campOut1, voutfilenamesCAmp[n].c_str(), eGRCBIN); //	write output GRC files
#endif // !CORRELATE_FRAMES
							}
					verraver[na] /= double(vndefocus[na]);
					vnormaver0[na] /= double(vndefocus[na]);
					vnormaver1[na] /= double(vndefocus[na]);
						}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
					}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

#ifdef CORRELATE_FRAMES
			XArData::WriteFileDAT(vCorr, "0Correlate.txt");
#endif // CORRELATE_FRAMES
#ifndef CORRELATE_FRAMES
			// calculate the current reconstruction error and
			// if the difference with the previous error is smaller than the defined minimum
			// or, if the error started to increase, interrupt the iterations
			double ssej = 0.0;
			for (index_t na = 0; na < nangles; na++) ssej += verraver[na];
			ssej /= double(nangles);

			printf("\n\nAverage relative error in reprojected defocused intensities = %g", ssej);
			for (index_t na = 0; na < nangles; na++) printf("\nIllum_angle_no. = %zd, C0(orig_contrast) = %g, C0(reproj_contrast) = %g, Rel.error = %g ", na, vnormaver0[na], vnormaver1[na], verraver[na]);

			// restore the 3D potential V in case we need to work with it further
			dnorm3 = 1.0 / dnorm3;
			K3out *= dnorm3;
#endif // !CORRELATE_FRAMES
				}

#ifndef CORRELATE_FRAMES
		// finding atomic positions
		if (imodePeaks)
		{
			printf("\nSearching for peak positions in the 3D distribution of the electrostatic potential ...");

			// exclude the points located outside the reconstruction volume from the subsequent peak search
			double xR = (xhi - xlo) / 2.0; // x-radius
			double yR = (yhi - ylo) / 2.0; // y-radius
			double zR = (zhidz - zlodz) / 2.0; // z-radius
			double R2 = (xR * xR + yR * yR + zR * zR) / 3.0;
			double K3min = K3out.Norm(eNormMin);

#pragma omp parallel for shared(K3out, K3min, xR, yR, zR, R2, xst, yst, zst)
			for (int k = 0; k < nz; k++)
			{
				double zzz = -zR + zst * k;
				zzz *= zzz;
				for (int j = 0; j < ny; j++)
				{
					double yyy = -yR + yst * j;
					yyy *= yyy; yyy += zzz;
					for (int i = 0; i < nx; i++)
					{
						double xxx = -xR + xst * i;
						xxx *= xxx; xxx += yyy;
						if (xxx > R2)
							K3out[k][j][i] = K3min; // mark points outside the reconstruction volume for exclusion from the peak search
					}
				}
			}

			int katom = int(datomsizeZ / zst + 0.5), jatom = int(datomsizeXY / yst + 0.5), iatom = int(datomsizeXY / xst + 0.5);
			int natom(0);
			vector<int> vimax, vjmax, vkmax;
			vector<float> xa, ya, za;

			// search for peaks
			// NOTE that we exclude one-atomsize-wide vicinity of the outer boundary from the search, as we expect artefacts there
#pragma omp parallel for shared(K3out, natom, katom, jatom, iatom, vimax, vjmax, vkmax)
			for (int k = katom; k < nz - katom * 2; k += katom)
			{
				for (int j = jatom; j < ny - jatom * 2; j += jatom)
				{
					for (int i = iatom; i < nx - iatom * 2; i += iatom)
					{
						double K3max(0);
						int kmax(0), jmax(0), imax(0);
						for (int kk = k; kk < k + katom; kk++)
							for (int jj = j; jj < j + jatom; jj++)
								for (int ii = i; ii < i + iatom; ii++)
									if (K3out[kk][jj][ii] > K3max)
									{
										K3max = K3out[kk][jj][ii];
										K3out[kmax][jmax][imax] = 0.0;
										kmax = kk; jmax = jj; imax = ii;
									}
									else K3out[kk][jj][ii] = 0.0;
#pragma omp critical
						if (K3max > 0)
						{
							natom++;
							vimax.push_back(imax);
							vjmax.push_back(jmax);
							vkmax.push_back(kmax);
						}
					}
				}
			}

			if (natom == 0) printf("\n\n!!WARNING: no peaks have been found!");
			else
			{

				// create vectors of peak location coordinates
				xa.resize(natom); ya.resize(natom); za.resize(natom);
				vector<Pair2> K3maxPair(natom);
				for (int n = 0; n < natom; n++)
				{
					xa[n] = float(xlo + xst * vimax[n]);
					ya[n] = float(ylo + yst * vjmax[n]);
					za[n] = float(zst * vkmax[n]);
					K3maxPair[n].v = K3out[vkmax[n]][vjmax[n]][vimax[n]];
					K3maxPair[n].n = n;
				}

				// sort the located peaks according to their height, in order to improve the subsequent procedure of elimination of closely located peaks
				printf("\nSorting the located peaks in descending order according to their heights ...");
				std::qsort(&K3maxPair[0], (size_t)natom, sizeof(Pair2), Pair2comp);

				// exclude the smaller one from each pair of peak positions that are located closer than datomsize to each other (e.g. in adjacent corners of neigbouring cubes)
				double datomsize2 = datomsizeXY * datomsizeXY;
				vector<int> vimax1, vjmax1, vkmax1, Znum1;
				vector<float> xa1, ya1, za1, occ1, wobble1;
				printf("\nEliminating adjacent peaks in the reconstructed 3D distribution ...");
				int n, m;
				for (int nn = natom - 1; nn >= 0; nn--)
				{
					if (K3maxPair[nn].v != 0.0)
					{
						// we can already count this maximum in, as it is guaranteed to be larger than all subsequent ones
						n = K3maxPair[nn].n;
						Znum1.push_back(6);
						vimax1.push_back(vimax[n]);
						vjmax1.push_back(vjmax[n]);
						vkmax1.push_back(vkmax[n]);
						xa1.push_back(xa[n]);
						ya1.push_back(ya[n]);
						za1.push_back(za[n]);
						occ1.push_back((float)K3maxPair[nn].v);
						wobble1.push_back(0);
#pragma omp parallel for private(m)
						for (int mm = nn - 1; mm >= 0; mm--)
						{
							if (K3maxPair[mm].v != 0.0)
							{
								m = K3maxPair[mm].n;
								if ((xa[n] - xa[m]) * (xa[n] - xa[m]) + (ya[n] - ya[m]) * (ya[n] - ya[m]) + (za[n] - za[m]) * (za[n] - za[m]) < datomsize2)
									K3maxPair[mm].v = 0.0;
							}
						}
					}
				}
				natom = (int)Znum1.size();
				printf("\n%d peak positions have been found in total.", natom);

				// copy the located isolated peaks back into the 3D potential distribution
				K3out.Fill(0.0);
				for (int n = 0; n < natom; n++)
					K3out[vkmax1[n]][vjmax1[n]][vimax1[n]] = occ1[n];

				// exclude the peaks located outside the minimal ball
				//double xac(50.0), yac(50.0), zac(50.0), xarad(41.0);
				//double xarad2 = xarad * xarad;
				//for (int n = 0; n < natom; n++)
				//		if ((xa[n] - xac) * (xa[n] - xac) + (ya[n] - yac) * (ya[n] - yac) + (za[n] - zac) * (za[n] - zac) > xarad2) 
				//			K3out[vkmax[n]][vjmax[n]][vimax[n]] = 0.0f;

				printf("\nSaving the reconstructed peak-localized 3D object into output files ...");
				if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesPeaksTot, eTIFF32);
				else XArData::WriteFileStackGRD(K3out, voutfilenamesPeaksTot, eGRDBIN);

				if (natom > NATOMMAX)
					printf("\nNot saving the reconstructed peak-localized 3D object into an XYZ file, since too many peaks (%d) have been found.", natom);
				else
				{
					printf("\nSaving the reconstructed peak-localized 3D object into an XYZ file ...");
					size_t dotpos = filenamebaseOutNew.find_last_of(".");
					string fileoutXYZ = filenamebaseOutNew.replace(dotpos + 1, filenamebaseOutNew.size() - dotpos - 1, "xyz");
					SaveXYZfile(fileoutXYZ, float(xhi - xlo), float(zhi - zlo), (int)natom, &Znum1[0], &xa1[0], &ya1[0], &za1[0], &occ1[0], &wobble1[0]);
				}
			}
		} // end of peak localization module
#endif // !CORRELATE_FRAMES
		if (bVerboseOutput)
		{
			printf("\nTotal file read time = %ld", liFileReadTime);
			printf("\nTotal file write time = %ld", liFileWriteTime);
			printf("\nTotal file 2D FFT and CTF correction time = %ld", li2DFFTtime);
			printf("\nTotal file 3D FFT time = %ld", li3DFFTtime);
			printf("\nTotal back propagation time = %ld", liBackPropTime);
		}

		// global cleanup of FFTW space (is must not be called from destructors of individual FFTW objects)
	}
	catch (std::runtime_error& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}
	catch (std::exception& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}

	std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
	printf("\n\nMain program finished. Execution time = %zd s.", std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	printf("\nPress any key to exit..."); getchar();
	return 0;
}