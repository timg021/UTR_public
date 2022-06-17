// UHR.cpp : This file contains the 'main' function which implements Unified Tomographic Reconstruction method

#include <chrono>
#include <omp.h>
#include <filesystem>

#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_spln2.h"
#include "XA_spln3.h"
#include "XA_iwfr.h"
#include "XA_tie.h"
#include "XA_tiff.h"
#include "XA_move2.h"
#include "fftwd3frc.h"

#include "UTR.h" // various include headers, constants and templates

using namespace xar;

int main(int argc, char* argv[])
{
	// start the execution timer
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now(), start_timeNow;
	long int liFileReadTime(0), liFileWriteTime(0), li2DFFTtime(0), li3DFFTtime(0), liBackPropTime(0), liPeakLocalizeTime(0); // time counters

	try
	{
		printf("\nStarting Unified Tomographic Reconstruction program ...");

		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;
		vector<Pair> v2shifts;
		vector<Pair> vastigm;

		//************************************ read input parameters from file
		// read input parameter file
		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024], cparam3[1024], cparam4[1024];

		string sInputParamFile("UTR.txt");
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
		if (ndefocusmax > 1)
			throw std::runtime_error("Error: this program is currently implemented only for defocus parameter files with one defocus distance per orientation.");

		// check if conditions for 3D CTF correction are satisfied
		bool bSingleDefocusDistance(true); // indicates if all defocus distances are the same
		int iAstigmatism(0); // 0 = no astigmatism, 1 = same astigmatism at all orientations, 2 = variable astigmatism
		double zoutAver(vvdefocus[0][0].b); // average defocus
		for (index_t i = 1; i < nangles; i++) // check if all first defocus distances are the same
		{
			if (bSingleDefocusDistance && (vvdefocus[i][0].b != vvdefocus[0][0].b)) bSingleDefocusDistance = false;
			zoutAver += vvdefocus[i][0].b;
		}
		zoutAver /= nangles;
		if (bRelion) // check for astigmatism
		{
			for (index_t i = 0; i < nangles; i++)
			{
				if (iAstigmatism == 0 && (vastigm[i].a != 0 || vastigm[i].b != 0)) iAstigmatism = 1;
				if (iAstigmatism != 2 && (vastigm[i].a != vastigm[0].a || vastigm[i].b != vastigm[0].b)) iAstigmatism = 2;
			}
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4. Input filename base of defocus series of the sample in TIFF, GRD or RAW format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		bool bRAWinput(false), bTIFFinput(false), bGRDinput(false);
		string strTemp = GetFileExtension(filenamebaseIn);
		if (strTemp == string(".TIFF") || strTemp == string(".TIF")) bTIFFinput = true;
		else if (strTemp == string(".RAW")) bRAWinput = true;
		else if (strTemp == string(".GRD")) bGRDinput = true; // this value is not used below at the moment, as GRD is considered to be a default when all other input formats are false
		else throw std::runtime_error("Error: input filename extension must be TIF, GRD or RAW.");

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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6. Input data normalization factors f1 f2 (input -> (input / f1) + f2)
		double dNormFactor1(1.0), dNormFactor2(0.0);
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::runtime_error("Error reading input data normalization factors from input parameter file.");
		dNormFactor1 = atof(cparam);
		dNormFactor2 = atof(cparam1);
		printf("\nInput data normalization factors: f1 = %g, f2 = %g", dNormFactor1, dNormFactor2);
		if (dNormFactor1 == 0) throw std::runtime_error("The first normalization factor cannot be zero (leads to division by zero).");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Width(pix) and height(pix) of input images
		// for RAW files, these parameter determine how the images are read; for TIFF and GRD files, these parameters may be used for trimming the frames
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::runtime_error("Error reading Width Height from input parameter file.");
		index_t nx = atoi(cparam), nx0;
		index_t ny = atoi(cparam1), ny0;
		printf("\nDimensions of input images (possibly, after trimming or padding): width = %zd, height = %zd (pixels)", nx, ny);
		int nxd2 = int(nx / 2), nyd2 = int(ny / 2);
		if (nx != 2 * nxd2 || ny != 2 * nyd2)
			throw std::runtime_error("Error: width and height parameters in input parameter file must be even.");
		int nx2 = int(nx) - 2, ny2 = int(ny) - 2;
		index_t nxny = (nx * ny);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Modality and units of length(UL): TEM and Angstroms(0) or X-ray CT and microns(1)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
			throw std::runtime_error("Error reading units of length from input parameter file.");
		int iModality = atoi(cparam);
		if (iModality == 0) // electrons and Angstroms (TEM)
			printf("\nUnits of length (UL) are set to Angstroms (TEM imaging mode will be used)");
		else if (iModality == 1) // X-rays and microns (X-ray CT)
				printf("\nUnits of length (UL) are set to microns (X-ray CT imaging mode will be used)");
			else
				throw std::runtime_error("Error: modality and units of length parameter in input parameter file must be 0 for TEM/Angstroms or 1 for X-ray CT/microns.");
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9. Pixel size in UL
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
		printf("\nPixel size = %g (UL)", xst);
		printf("\nPhysical boundaries of input images: Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g (UL)", xlo, xhi, ylo, yhi);
		// NOTE that the input 2D images will be trimmed or padded to the (ny, nx) size, if they don't have this size already
		double fxst = 1.0 / (xhi - xlo), fyst = 1.0 / (yhi - ylo);
		double fxlo = -fxst * nxd2, fylo = -fyst * nyd2;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10. Output defocus distances min and max in UL
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading output defocus distances from input parameter file.");
		double zlo = atof(cparam); // minimum output defocus in UL - !!! will be corrected with dzextra below
		double zhi = atof(cparam1); // maximum output defocus in UL - !!! will be corrected with dzextra below 
		double zst = xst; // this is a provision for possible future extensions to non-qubic voxels
		if (zlo > zhi) std::swap(zlo, zhi);
		index_t nz = index_t((zhi - zlo) / zst + 0.5); // number of defocus planes to propagate to
		if (nz <= 0)
			throw std::runtime_error("Error: number of z steps is not positive.");
		int nzd2 = int(nz / 2), nz2 = int(nz) - 2;
		if (nz != 2 * nzd2)
			throw std::runtime_error("Error: number of z points must be even.");
		vector<double> voutdefocus(nz); // vector of output defocus distances
		for (index_t n = 0; n < nz; n++) voutdefocus[n] = zlo + zst * n;
		double fzst = 1.0 / (zhi - zlo);
		double fzlo = -fzst * nzd2;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //11. Centre of rotation x, y and z shifts in UL
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4)
			throw std::runtime_error("Error reading ccentre of rotation x, y and z shifts from input parameter file.");
		double dxc = atof(cparam);
		double dyc = atof(cparam1);
		double dzc = atof(cparam2);
		printf("\nCentre of rotation shifts in UL: dXc = %g, dYc = %g, dZc = %g (UL)", dxc, dyc, dzc);
		if (abs(dxc) > (xhi - xlo) || abs(dyc) > (yhi - ylo))
			printf("\nWARNING: X or Y shifts of the centre of rotation will take the centre outside the reconstruction volume.");
		bool bRotCentreShift = (dxc != 0 || dyc != 0 || dzc != 0) ? true : false;
		double xc = (xhi + xlo) / 2.0 + dxc; // x-coordinate of the centre of rotation
		double yc = (yhi + ylo) / 2.0 + dyc; // y-coordinate of the centre of rotation
		double zc = (zhi + zlo) / 2.0 + dzc; // z-coordinate of the centre of rotation

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12. Wavelength in UL
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading wavelength from input parameter file.");
		double wl = atof(cparam); // wavelength in UL
		printf("\nWavelength = %g (A)", wl);
		if (wl <= 0 || wl > 1)
			throw std::runtime_error("Error: wavelength value appears to be wrong.");
		double EE; // incident electron energy in volts (recalculated from the wavelength below)
		if (iModality == 0)
		{
			constexpr double hp = 6.62607004e-34; // Planck's constant (m2 kg / s)
			constexpr double cc = 299792458; // speed of light (m / s)
			constexpr double ee = 1.602176634e-19; // electron charge (coulomb)
			constexpr double m0 = 9.1093837015e-31; // electron rest mass (kg)
			constexpr long double mc2 = m0 * cc * cc; // mc^2
			long double chl2 = (long double)(cc * cc * hp * hp) / (long double)(wl * wl * 1.e-20);
			long double abra = sqrt(mc2 * mc2 + chl2);
			EE = double((-mc2 + abra) / (long double)(ee));
			printf("\nIncident electron energy E = %g (keV)", EE / 1000.0);
		}
		else
		{
			EE = 0.001239841662513396 / wl;
			printf("\nIncident X-ray energy E = %g (keV)", EE);
		}
		double DOF = xst * xst / (2.0 * wl); // depth of field
		printf("\nDepth of field = %g (UL)", DOF);
		if (DOF >= (zhi - zlo)) printf("\nDepth of field is larger than the thickness of the reconstruction volume - Ewald sphere curvature will be unimportant.");
		else printf("\nDepth of field is smaller than the thickness of the reconstruction volume - Ewald sphere curvature may be important.");
		double wl2 = 0.5 * wl;

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13. Objective aperture in mrad
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading objective aperture from input parameter file.");
		double aobj = atof(cparam);
		if (aobj != 0) printf("\nObjective aperture (half-angle) = %g (mrad)", aobj);
		else  printf("\nObjective aperture is infinite.");
		if (aobj < 0 || aobj > 1000)
			throw std::runtime_error("Error: objective aperture value appears to be wrong.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14. Spherical aberrations Cs3 and Cs5 in mm
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading spherical aberrations from input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> UL
		Cs5 *= 1.e+7; // mm --> UL
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15.Phase_retrieval_method:_CTF-2D(0),_CTF-3D(1),_TIE-Hom-3D(2)_or_none(3): 2
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading CTF correction phase retrieval method parameter from input parameter file.");
		int iCTFcorrectionMode = atoi(cparam);
		switch (iCTFcorrectionMode)
		{
		case 0:
			printf("\nNo CTF correction will be applied.");
			break;
		case 1:
			printf("\nTIE-Hom correction will be applied in 2D to the input images.");
			if (bRelion) printf("\nWARNING: spherical aberrations and astigmatism have not been implemented yet in TIE-Hom phase retrieval and will be ignored.");
			break;
		case 2:
			printf("\nTIE-Hom correction will be applied in 3D to the CT-reconstructed volume.");
			if (bRelion) printf("\nWARNING: spherical aberrations and astigmatism have not been implemented yet in TIE-Hom phase retrieval and will be ignored.");
			break;
		case 3:
			printf("\nCTF correction will be applied in 2D to the input images.");
			break;
		case 4:
			printf("\nCTF correction will be applied in 3D to the CT-reconstructed volume.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for phase retrieval parameter in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 16.Apply Ewald sphere curvature correction: 0=No or 1=Yes
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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17. Absorption fraction beta / delta
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading absorptioin coefficient from input parameter file.");
		double asigma = atof(cparam);
		printf("\nAbsorption fraction beta / delta = %g", asigma);
		if (asigma < 0) throw std::runtime_error("Error: absorption fraction in input parameter file cannot be negative.");
		if (asigma == 0 && (iModality == 1 || iCTFcorrectionMode == 1 || iCTFcorrectionMode == 2)) 
			throw std::runtime_error("Error: absorption fraction in input parameter file must be positive in X-ray CT mode or when using TIE-Hom.");
		//if (iCTFcorrectionMode == 1 || iCTFcorrectionMode == 2) printf("\nEffective regularization parameter alpha for inverse 3D Laplacian = %g", asigma / (PI * zoutAver * wl));

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 18. Average noise to signal ratio (1/SNR) in input images (1/sqrt(Nphot))
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading average noise-to-signal ratio in input images from input parameter file.");
		double dInputNSR = atof(cparam); // average NSR in the input images
		if (dInputNSR < 0)
			throw std::runtime_error("Error: NSR in input images must be non-negative.");
		printf("\nAverage NSR in input images = %g", dInputNSR);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19. Tikhonov regularization parameter
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading minimal phase reconstruction error/regulariztion parameter from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nTikhonov regularization parameter = %g", epsilon);
		if (epsilon < 0)
			throw std::runtime_error("Error: the regulalization parameter must be non-negative.");
		if (epsilon == 0 && asigma == 0)
			throw std::runtime_error("Error: when beta / delta = 0, epsilon must be positive.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //20. Enforce symmetry: not_apply(0), post-apply(1), or distribute input orientations(2)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading enforcing symmetry mode from input parameter file.");
		int imodeSymmetry = atoi(cparam);
		switch (imodeSymmetry)
		{
		case 0:
			printf("\nKnown symmetry won't be applied.");
			break;
		case 1:
			//printf("\nInput orientations will be distributed according to known invariant rotational positions.");
				throw std::runtime_error("Error: symmetry enforcing mode has not been implemented in this program yet.");
			break;
		case 2:
			printf("\nKnown symmetry will be applied after initial reconstruction or import of a 3D potential or beta.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for symmetry enforcing mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 21.Input file with rotation angles enforcing symmetry
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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 22. 3D Laplacian filter mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading 3D Laplacian filter mode from input parameter file.");
		int imodeInvLaplace = atoi(cparam);
		switch (imodeInvLaplace)
		{
		case 0:
			printf("\n3D Laplacian filter won't be applied.");
			break;
		case 1:
			if (iCTFcorrectionMode == 1 || iCTFcorrectionMode == 2) printf("\nWARNING: 3D Laplacian filter won't be applied, because TIE-Hom phase retrieval has been already selected.");
			else printf("\n3D Laplacian filter will be applied.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for 3D Laplacian filter mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 23. Low-pass filter width in UL, background subtraction value and lower threshold level in Volts
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::runtime_error("Error reading low-pass filter width, background subtraction value and lower threshold level from input parameter file.");
		double dlpfiltersize = atof(cparam);
		double dBackground = atof(cparam1);
		double dThreshold = atof(cparam2);
		printf("\nLow-pass filter width for 3D potential or beta = %g (UL)", dlpfiltersize);
		printf("\nBackground subtraction value for 3D potential or beta = %g (Volts)", dBackground);
		printf("\nLower threshold level for 3D potential or beta = %g (Volts)", dThreshold);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 24. Peak localization mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading peak localization mode from input parameter file.");
		int imodePeaks = atoi(cparam);
		switch (imodePeaks)
		{
		case 0:
			printf("\nPeak localization in the 3D electrostatic potential or beta won't be applied.");
			break;
		case 1:
			printf("\nPeak localization in the 3D electrostatic potential or beta will be applied.");
			break;
		default:
			throw std::runtime_error("Error: unknown value for the peak localization mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 25. Transverse and longitudinal side lengths for peak localization (in_UL)
		if (sscanf(cline, "%s %s %s", ctitle, cparam ,cparam1) != 3) throw std::runtime_error("Error reading transverse and longitudinal side lengths for peak localization from input parameter file.");
		double datomsizeXY = atof(cparam), datomsizeZ = atof(cparam1);
		printf("\nTransverse and longitudinal side lengths for peak localization = %g %g (UL)", datomsizeXY, datomsizeZ);
		if (imodePeaks && (int(datomsizeXY / zst + 0.5) < 2 || int(datomsizeZ / zst + 0.5) < 2))
			throw std::runtime_error("Error: transverse and longitudinal side lengths for peak localization must be 2 x z_step or larger.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 26. Output file name base in GRD or TIFF format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		printf("\nFile name base for 3D potential or beta = %s", filenamebaseOut.c_str());
		strTemp = GetPathFromFilename(filenamebaseOut, false);
		if (!std::filesystem::exists(strTemp))
			throw std::runtime_error("Error: the specified file folder for 3D potential or beta does not seem to exist.");
		bool bTIFFoutput;
		if (GetFileExtension(filenamebaseOut) == string(".TIFF") || GetFileExtension(filenamebaseOut) == string(".TIF")) bTIFFoutput = true;
		else if (GetFileExtension(filenamebaseOut) == string(".GRD")) bTIFFoutput = false;
		else throw std::runtime_error("Error: output filename extension must be TIF ot GRD.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 27. Save or not the sampling matrix in files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading save_or_not inverse sampling matrix switch from input parameter file.");
		int nSaveDefocCAmpsOrSamplingMatrix = atoi(cparam);
		if (nSaveDefocCAmpsOrSamplingMatrix == 1)
			printf("\n3D sampling matrix will be saved in GRD files");
		else
			printf("\nSampling matrix will not be saved in files");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 28. Import and reprocess existing 3D potential or beta files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading import and reprocess existing 3D potential or beta files switch from input parameter file.");
		int imode3DPotential = atoi(cparam);
		switch (imode3DPotential)
		{
		case 0:
			printf("\nThis program will reconstruct 3D potential or beta from defocused images and save it in output files.");
			break;
		case 1:
			printf("\nThis program will import existing 3D potential or beta from files and reprocess it.");
			break;
		default:
			throw std::runtime_error("Error: unknown value of import and reprocess existing 3D potential or beta files switch in input parameter file.");
		}

		if (imode3DPotential == 1) // only read in pre-caculated 3D potential or beta from "output" files and reprocess them
		{
			bGRDinput = bTIFFinput = bRAWinput = false;
			string strTemp = GetFileExtension(filenamebaseOut);
			if (strTemp == string(".TIFF") || strTemp == string(".TIF")) bTIFFinput = true;
			else if (strTemp == string(".GRD")) bGRDinput = true;
			else if (strTemp == string(".RAW")) bRAWinput = true;
			else throw std::runtime_error("Error: input filename extension (in this mode - it is taken from the output filename template) must be TIF, GRD or RAW.");
		}

		if (imode3DPotential != 1 && (iCTFcorrectionMode == 2 || iCTFcorrectionMode == 4) && (!bSingleDefocusDistance || iAstigmatism != 0))
		{
			printf("\n\nWARNING: different defocus distances and/or non-zero astigmatism have been detected in the input file %s.", defocfile.c_str());
			printf("\n3D CTF correction have only been implemented for constant defocus distance and no astigmatism, possibly leading to inaccurate results in this case.");
			printf("\nDo you still want to continue (y = yes, any other symbol = no): ");
			int c = getchar();
			if (c != 'y' && c != 'Y') exit(1);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // // 29. Folder name for auxiliary files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading folder name for auxiliary files from input parameter file.");
		string folderAux = cparam;
		printf("\nFolder for auxiliary file output = %s", folderAux.c_str());
		if (folderAux.rfind('\\') != folderAux.size() - 1 && folderAux.rfind('/') != folderAux.size() - 1) // the last character is not '\' and not '/'
			folderAux.append("/");
		if (!std::filesystem::exists(folderAux))
			throw std::runtime_error("Error: the specified auxiliary file folder does not seem to exist.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 30. Number of parallel threads
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading number of parallel threads from input parameter file.");
		int nThreads = atoi(cparam);
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
		if (imode3DPotential) // only read in pre-caculated 3D potential or beta from "output" files and reprocess them
		{
			printf("\nInput file name base for pre-existing 3D potential or beta = %s", filenamebaseOut.c_str());
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
		vector<string> voutfilenamesTot; // output filenames for the reconstructed 3D potential or beta
		if (imode3DPotential) // output the the renormalized 3D potential or beta
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "R");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the renormalized 3D potential or beta = %s", filenamebaseOutNew.c_str());
			FileNames(1, nz, filenamebaseOutNew, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the filtered and/or rescaled 3D object
		}
		else
		{
			FileNames(1, nz, filenamebaseOut, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
		}

		vector<string> voutfilenamesPeaksTot; // output filenames for the peak-localized reconstructed 3D potential or beta
		if (imodePeaks) // output filenames for the peak-localized 3D potential or beta
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "A");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the peak-localized 3D potential or beta = %s", filenamebaseOutNew.c_str());
			FileNames(1, nz, filenamebaseOutNew, voutfilenamesPeaksTot); // create 1D array of output filenames to save 2D slices of the peak-localized 3D object
		}

		vector<string> voutfilenamesTotDefocCAmpOrFM; // output filenames for sampling matrix
		if (nSaveDefocCAmpsOrSamplingMatrix == 1) // create filenames for saving the sampling matrix
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutDefocCAmpOrFM = apath.filename().string();
			filenamebaseOutDefocCAmpOrFM.replace(filenamebaseOutDefocCAmpOrFM.find_last_of("."), filenamebaseOutDefocCAmpOrFM.length() - filenamebaseOutDefocCAmpOrFM.find_last_of("."), "FM.grd");
			filenamebaseOutDefocCAmpOrFM = folderAux + filenamebaseOutDefocCAmpOrFM;
			FileNames(nz, 1, filenamebaseOutDefocCAmpOrFM, voutfilenamesTotDefocCAmpOrFM);

			printf("\nOutput file name base for sampling matrix = %s", filenamebaseOutDefocCAmpOrFM.c_str());
		}

		//************************************ end creating vectors of input and output file names


		//*********************************** start main calculations
		bool bAbort(false);

		std::unique_ptr<IXAHead> pHead(new Wavehead2D(wl, ylo, yhi, xlo, xhi)); // any header data in input image files will be ignored

		printf("\n\nPreparing FFTW library ...");
		// single Fftwd2c object is used for repeated 2D FFTs using FFTW later (the plan is created at this point and reused later)
		Fftwd2c fftw2D((int)ny, (int)nx, true);

		XArray3D<float> K3out; // big 3D reconstructed array (needs to fit into RAM alongside with with everything else)

		if (!imode3DPotential) // do phase retrieval and backpropagation prior to 3D filtering and output
		{
			index_t naSym{ 0 }, naSymCurrent{ 0 }; // indexes of illumination angles with respect to subdivision into symmetry subsets (only used when imodeSymmetry == 1)
			XArray3D<fcomplex> V3; // 3D potential or beta in the Fourier space
			XArray3D<float> Samp3; // universal 3D sampling matrix in the Fourier space 
			XArray3D<float> K3outsum; // accumulated K3out (only used when imodeSymmetry == 1)

			// allocate the large 3D output arrays
			printf("\nAllocating 3D arrays ...");
			V3.Resize(nz, ny, nx);
			V3.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
			Samp3.Resize(nz, ny, nx);
			Samp3.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));

			if (imodeSymmetry == 1) // !!!@@@@@ the implementation of this mode must be rechecked and redone
			{
				K3outsum.Resize(nz, ny, nx);
				K3outsum.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
			}

			// start of cycle over illumination angles
			printf("\nPerforming phase retrieval and backpropagation into the 3D volume ...");
			#pragma omp parallel for 
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
					if (bVerboseOutput) printf("\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);
					else printf("\n*** Illumination angle[%d] = (%g, %g) (degrees)", na, angleZ / PI180, angleY / PI180);

					// Centre of rotation expressions (linear phase factors)
					dcomplex Fxc, Fyc, Fzc;
					if (bRotCentreShift)
					{
						Fxc = dcomplex(0, -1) * tPI * (dxc * (1.0 - cosangleY * cosangleZ) - dyc * sinangleZ + dzc * sinangleY * cosangleZ);
						Fyc = dcomplex(0, -1) * tPI * (dxc * cosangleY * sinangleZ + dyc * (1.0 - cosangleZ) - dzc * sinangleY * sinangleZ);
						Fzc = dcomplex(0, -1) * tPI * (-dxc * sinangleY + dzc * (1.0 - cosangleY));
					}

					// start the execution timer per angle
					std::chrono::system_clock::time_point start_timeA;
					if (bVerboseOutput) start_timeA  = std::chrono::system_clock::now();

					//index_t ndefocus = vndefocus[na]; // number of defocus planes at the current illumination angle
					// NOTE: in the current version of the program, ndefocus is always 1
					vector<Pair> vdefocus = vvdefocus[na]; // vector of input Z" rotation angles and defocus positions at the current illumination angle
					vector<string> vinfilenames = vvinfilenames[na]; // input filenames of defocused images at the current illumination angle
					vector<string> voutfilenames(nz);

					XArray2D<double> int0; // input defocused intensity image
					double zout = vdefocus[0].b; // defocus distance

					// read defocused images from files
					if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of file read time
					//if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (UL)", vdefocus[0].a, vdefocus[0].b);
					//if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[0].c_str());
					if (bRAWinput) XArData::ReadFileRAW(int0, vinfilenames[0].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
					else if (bTIFFinput) TIFFReadFile(int0, vinfilenames[0].c_str()); // read input TIFF files
					else XArData::ReadFileGRD(int0, vinfilenames[0].c_str(), wl); //	read input GRD files
					if (bVerboseOutput) liFileReadTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();

					int ny0 = int0.GetDim1();
					int nx0 = int0.GetDim2();
					if (nx0 != nx || ny0 != ny)
					{
						if (nx0 < nx && ny0 <= ny || ny0 < ny && nx0 <= nx)
						{
							//printf("\nWARNING: 2D array from the image file will be padded with 1.0 to the array dimensions (Width = %zd, Height = %zd).", nx, ny);
							XArray2DMove<double> tmp2(int0);
							tmp2.Pad((ny - ny0) / 2, ny - ny0 - (ny - ny0) / 2, (nx - nx0) / 2, nx - nx0 - (nx - nx0) / 2, 1.0);
							// NOTE that the reconstructed 3D volume will NOT be trimmed back before saving, as it gets too complicated and error prone
						}
						else if (nx0 > nx && ny0 >= ny || ny0 > ny && nx0 >= nx)
						{
							//printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions (Width = %zd, Height = %zd).", nx, ny);
							XArray2DMove<double> tmp2(int0);
							tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
							// NOTE that the reconstructed 3D volume will NOT be padded back before saving, as it does not seem to make sense
						}
						else
							throw std::runtime_error("Error: the dimensions of the 2D array in the image file are inconsistent with the array dimensions in input parameter file (both dimensions should be equal, smaller or larger simultaneously).");
					}
					int0.SetHeadPtr(pHead->Clone());
					// NOTE that after the above trimming or padding, nx, ny, xlo, xhi, ylo, yhi parameters will correspond to those defined after parameter 7 above

					// renormalize input data
					if (dNormFactor1 != 1.0) int0 /= dNormFactor1;
					if (dNormFactor2 != 0.0) int0 += dNormFactor2;
					if (int0.Norm(eNormMin) < 0)
					{
						//printf("\nWARNING: negative values in the input intensity file.");
						int0.ThresholdLow(0.0, 0.0);
					}
					//double averInt = int0.Norm(eNormAver);
					//if (averInt == 0) throw std::runtime_error("Error: zero average value in the input file after normalization.");
					//if (averInt != 1.0) int0 /= averInt; // enforce unit average intensity

					// rotate input defocused image around Z'' back to zero angle
					if (vdefocus[0].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
					{
						double aver = int0.NormAverEdge(5);
						// shift along X and/or Y back to the unshifted position
						if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
						{
							XArray2DMove<double> xamove(int0);
							xamove.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), aver);
						}
						// rotate input defocused complex amplitude around Z'' back to zero angle
						if (vdefocus[0].a != 0)
						{
							XArray2D<double> vintTemp(int0);
							XArray2DSpln<double> xaSpln(vintTemp);
							xaSpln.Rotate(int0, -vdefocus[0].a, 0.5 * (yhi + ylo) + dyc, 0.5 * (xhi + xlo) + dxc, aver); // expecting uniform background
						}
					}

					// now do 2D CTF correction
					XArray2D<dcomplex> K2four; // symmetrized backpropagated contrast in the reciprocal space

					if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now();  // start of 2D FFT time

					XA_IWFR<double> xa_iwfr;
					switch (iCTFcorrectionMode)
					{
					case 1: // 2D TIE-Hom correction
						if (zout != 0) // if zout == 0, no TIE-Hom correction is applied (corresponding to the conventional=absorption CT case)
						{
							//if (bVerboseOutput)	printf("\nApplying 2D TIE-Hom correction to the input image ...");
							//XA_2DTIE<double> xa_tie;
							//xa_tie.DP(int0, 1.0 / asigma, zout);
							int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0 (after the flat field correction)
							K2four = MakeComplex(int0, 0.0, false);
							double alphatemp = 4.0 * PI * asigma / (zout * wl);
							double normtemp = alphatemp * sqrt(1.0 + asigma * asigma) / asigma;
							fftw2D.InverseMLaplacian1(K2four, alphatemp, normtemp);
						}
						break;
					case 0: // no CTF correction
					case 2: // 3D TIE-Hom correction later
					case 4: // 3D CTF correction later
						//if (bVerboseOutput && iCTFcorrectionMode != 1) printf("\nCalculating 2D FFT of the input image ...");
						int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0 (after the flat field correction)
						K2four = MakeComplex(int0, 0.0, false);
						K2four.Shuffle();
						fftw2D.ForwardFFT(K2four);
						K2four.Shuffle();
						break;
					case 3: // 2D CTF correction
						//if (bVerboseOutput) printf("\nApplying 2D CTF correction to the input image ...");
						int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0 (after the flat field correction)
						if (bRelion)
							xa_iwfr.InvertCTF_DT1(int0, K2four, fftw2D, zout, vastigm[na].a, vastigm[na].b * PI180, k2maxo, Cs3, Cs5, asigma, epsilon, bESCC, 1);
						else
							xa_iwfr.InvertCTF_DT1(int0, K2four, fftw2D, zout, 0, 0, k2maxo, Cs3, Cs5, asigma, epsilon, bESCC, 1);
						// the output (K2four) is normalized as delta * 4 * PI * sqrt(1 + sigma^2) / wl
						break;
					}
					if (bSelectFrames) K2four *= sfWeights[na];

					if (bVerboseOutput) li2DFFTtime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();


					// Update 3D object at the current illumination direction on the Ewald sphere in 3D Fourier space
					if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of back-propagation time
					//if (bVerboseOutput) printf("\nUpdating 3D reconstructed object on the Ewald sphere ...");
					CT_3Dgridding(K2four, V3, Samp3, angleY, angleZ, fxlo, fxst, fylo, fyst, fzlo, fzst, Fxc, Fyc, Fzc, wl, bRotCentreShift, bESCC);
					if (bVerboseOutput)
					{
						liBackPropTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
						std::chrono::system_clock::time_point end_timeA = std::chrono::system_clock::now();
						//printf("\nExecution time for this illumination angle = %zd ms.", std::chrono::duration_cast<std::chrono::milliseconds>(end_timeA - start_timeA).count());
					}
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			} // end of cycle over illumination angles "na"
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");

			// optionally save the sampling matrix in files
			if (nSaveDefocCAmpsOrSamplingMatrix == 1)
			{
				printf("\n\nSaving the sampling matrix in files %s, etc. ...", voutfilenamesTotDefocCAmpOrFM[0].c_str());
				XArData::WriteFileStackGRD(Samp3, voutfilenamesTotDefocCAmpOrFM, eGRDBIN);
			}

			// renormalization for the detected sampling and noise filtering
			printf("\n\nNormalizing for the detected sampling and applying noise filter ...");
			dInputNSR *= sqrt(nx * ny); // change of the STD of noise after the 2D DFT; later it will be multiplied by sqrt(F)
#if 0
			//@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  start of temporary test code
			K3out = Abs(V3);
			XArData::WriteFileStackGRD(K3out, voutfilenamesTot, eGRDBIN);
			Samp3.SetHeadPtr(new Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
			//CTFcorrectedNoise.SetHeadPtr(new Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
			Samp3 ^= 0.5f;
			//Samp3 *= CTFcorrectedNoise;
			Samp3 *= float(dInputNSR);
			XArData::WriteFileStackGRD(Samp3, voutfilenamesTotDefocCAmpOrFM, eGRDBIN);
			exit(-22);
			//@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end of temporary test code
#endif
			
			// we want abs(V3) >  sqrt(2) * dInputNSR * sqrt(nx * ny) * sqrt(F), i.e. (signal + noise)^2 > 2 * variance(noise)
			#pragma omp parallel for shared (V3, Samp3)
			for (int k = 0; k < nz; k++)
			{
				float x, y;
				for (int j = 0; j < ny; j++)
					for (int i = 0; i < nx; i++)
					{
						x = Samp3[k][j][i];
						if (x == 0) { V3[k][j][i] = 0.0f; continue; }
						if (dInputNSR == 0)
						{
							V3[k][j][i] *= (1 / x);
						}
						else
						{
							y = abs(V3[k][j][i]) / float(dInputNSR * sqrt(x)); // SNR in this Fourier coefficient
							V3[k][j][i] *= (1 / x) * (0.5f + atan(100.0f * (y - fTargetSNR)) / PIf); // suppress V3 elements with SNR < TargetSNR
							//if (y < fTargetSNR) V3[k][j][i] = 0.0f; // "hard" cut-off filter
						}
					}
			}

			// delete Samp3 matrix to free memory prior to 3D FFT
			Samp3.Truncate();

			// do 3D CTF correction if necessary
			if (zoutAver != 0)
			{
				if (iCTFcorrectionMode == 2)
				{
					printf("\nPerforming 3D TIE-Hom correction ...");
					double alphatemp = 4.0 * PI * asigma / (zoutAver * wl);
					double normtemp = alphatemp * sqrt(1.0 + asigma * asigma) / asigma;
					InverseMLaplacian3D(V3, zhi - zlo, yhi - ylo, xhi - ylo, alphatemp, normtemp);
					// the output(V3) is normalized as delta * 4 * PI * sqrt(1 + sigma ^ 2) / wl
				}
				else if (iCTFcorrectionMode == 4)
				{
					printf("\nPerforming 3D CTF correction ...");
					CTFcorrection3D(V3, Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi), zoutAver, k2maxo, Cs3, Cs5, asigma, epsilon, bESCC);
				}
			}

			// inverse Laplace filtering
			if (imodeInvLaplace && (iCTFcorrectionMode != 1 && iCTFcorrectionMode != 2))
			{
				printf("\nInverse Laplace filtering the 3D distribution of the electrostatic potential or beta ...");
				InverseMLaplacian3D(V3, zhi - zlo, yhi - ylo, xhi - ylo, epsilon, epsilon);
			}

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
						V3[k][j][ia] = 0.5f * (V3[k][j][ia] + std::conj(V3[ka][ja][i]));
					}
				}
			}

			// inverse FFT of the FFT[potential or beta] on the Ewald sphere
			printf("\nAllocating Fftwd3frc object ...");
			Fftwd3frc fft3f((int)nz, (int)ny, (int)nx, false, true, false); // create an Fftwd3frc object in the Fourier space state
			fft3f.SetComplexXArray3D(V3, true);
			V3.Truncate(); // free memory
			printf("\nInverse Fourier transforming the Fourier transform of the 3D distribution of the electrostatic potential or beta ...");
			if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of inverse 3D FFT time
			fft3f.InverseFFT();
			if (bVerboseOutput) li3DFFTtime = (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
			K3out.Resize(nz, ny, nx); // allocated memory
			K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
			fft3f.GetRealXArray3D(K3out, true);
			fft3f.Cleanup(); // free memory

			float fnorm = (float)(wl * EE / (2.0 * PI * sqrt(1.0 + asigma * asigma))); // this is the conversion factor for V = 2E delta
			if (iModality == 1)
			{
				fnorm *= (float)(asigma / (2.0 * EE)); // this is the conversion factor for beta = sigma delta = sigma * V / (2E)
				if (iCTFcorrectionMode == 2 || iCTFcorrectionMode == 3) fnorm *= -1;
			}
			K3out *= fnorm;

		} // end of case if imode3DPotential == 0
		else // read in pre-existing 3D distribution of the electrostatic potential or beta
		{
			printf("\n\nReading pre-existing 3D distribution of the electrostatic potential or beta from files ...");

			K3out.Resize(nz, ny, nx);
			K3out.SetHeadPtr(new xar::Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));

			// read input 2D slices
			if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of file read time
			#pragma omp parallel for shared (K3out)
			for (int n = 0; n < nz; n++)
			{
				if (bAbort) continue;
				try
				{
					XArray2D<float> ipIn;
					if (bRAWinput) XArData::ReadFileRAW(ipIn, vinfilenamesTot[n].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
					else if (bTIFFinput) TIFFReadFile(ipIn, vinfilenamesTot[n].c_str()); // read input TIFF files
					else XArData::ReadFileGRD(ipIn, vinfilenamesTot[n].c_str(), wl); //	read input GRD files

					index_t ny0 = ipIn.GetDim1();
					index_t nx0 = ipIn.GetDim2();
					if (nx0 < nx || ny0 < ny) // there seems to be no need to allow padding when re-processing previously reconstructed 3D distributions
						throw std::runtime_error("Error: the dimensions of the 2D array in the image file are smaller than the array dimensions in Width Height HeaderLength Endianness ElementLength line from input parameter file.");
					else if (nx0 > nx || ny0 > ny)
					{
						printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions Width = %zd Height = %zd from input parameter file", nx, ny);
						XArray2DMove<float> tmp2(ipIn);
						tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
					}

					//if (dNormFactor1 != 1.0) ipIn /= dNormFactor1;
					//if (dNormFactor2 != 0.0) ipIn += dNormFactor2;

					for (index_t j = 0; j < ny; j++)
						for (index_t i = 0; i < nx; i++)
							K3out[n][j][i] = ipIn[j][i];
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
			if (bVerboseOutput) liFileReadTime = (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
		}

		// enforce known symmetries
		if (imodeSymmetry == 2)
		{
			printf("\nEnforcing known symmetries of the reconstructed 3D potential or beta ...");
			XArray3DSpln<float> xa3spln(K3out);
			XArray3D<float> K3outRot, K3outsum(K3out);
			for (int naSym = 1; naSym < nanglesSym; naSym++)
			{
				printf("\nRotating around Z by %g (deg), around Y' by %g (deg) and around Z'' by %g (deg) ...", -v2anglesSym[naSym].a, -v2anglesSym[naSym].b, -vvdefocusSym[naSym][0].a);
				xa3spln.Rotate3(K3outRot, -vvdefocusSym[naSym][0].a * PI180, -v2anglesSym[naSym].b * PI180, -v2anglesSym[naSym].a * PI180, nThreads); // !! rotate in the opposite direction
				K3outsum += K3outRot;
			}
			K3outsum /= float(nanglesSym);
			K3out = K3outsum;
		}

		// apply the regularized inverse 3D (-Laplacian) and high-pass filter
		if (imodeInvLaplace && imode3DPotential || (dlpfiltersize > 0 && dlpfiltersize < (xhi - xlo)))
		{
			if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of inverse 3D FFT time

			Fftwd3frc fft3f((int)nz, (int)ny, (int)nx);

			if (imodeInvLaplace && imode3DPotential) // in the case !imode3DPotential inverse Laplacian is calculated in the Fourier space
			{
				printf("\nInverse Laplace filtering the reconstructed 3D object ...");
				fft3f.InverseMLaplacian(K3out, epsilon, epsilon);
			}

			if (dlpfiltersize > 0 && dlpfiltersize < (xhi - xlo))
			{
				printf("\nHigh-pass filtering the reconstructed 3D object ...");
				XArray3D<float> K3outTemp(K3out);
				fft3f.GaussFilter(K3outTemp, dlpfiltersize / 2.355);
				K3out -= K3outTemp;
			}

			if (bVerboseOutput) li3DFFTtime = (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
		}

		if (dBackground != 0 || K3out.Norm(eNormMin) < dThreshold)
		{
			printf("\nThresholding the reconstructed 3D object ...");
			K3out.ThresholdLow((float)dBackground, (float)dThreshold);
		}

		// output the 3D array
		printf("\n\nSaving the reconstructed 3D object into output files %s, etc. ...", voutfilenamesTot[0].c_str());
		if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of file write time
		if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesTot, eTIFF32);
		else XArData::WriteFileStackGRD(K3out, voutfilenamesTot, eGRDBIN);
		if (bVerboseOutput) liFileWriteTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();

		// find peaks in the reconstructed 3D distribution
		if (imodePeaks)
		{
			printf("\nSearching for peak positions in the 3D distribution of the electrostatic potential or beta ...");
			if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of peak localization time

			int natom = FindPeaks(K3out, datomsizeXY, datomsizeZ, natommax, filenamebaseOutNew);

			printf("\nSaving the reconstructed peak-localized 3D object into output files %s, etc. ...", voutfilenamesPeaksTot[0].c_str());
			if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesPeaksTot, eTIFF32);
			else XArData::WriteFileStackGRD(K3out, voutfilenamesPeaksTot, eGRDBIN);

			if (bVerboseOutput) liPeakLocalizeTime = (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
		} // end of peak localization module

		if (bVerboseOutput)
		{
			printf("\nTotal file read time = %ld (ms)", liFileReadTime);
			printf("\nTotal file write time = %ld  (ms)", liFileWriteTime);
			printf("\nTotal 2D FFT and CTF correction time = %ld (ms)", li2DFFTtime);
			printf("\nTotal backpropagation time = %ld (ms)", liBackPropTime);
			printf("\nTotal 3D FFT time = %ld (ms)", li3DFFTtime);
			printf("\nTotal peak localization time = %ld (ms)", liPeakLocalizeTime);
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

