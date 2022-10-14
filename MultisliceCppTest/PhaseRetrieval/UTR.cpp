// UHR.cpp : This file contains the 'main' function which implements Unified Tomographic Reconstruction method

#include <chrono>
#include <omp.h>
#include <filesystem>
#include <map>

#include "XArray2D.h"
#include "XA_data.h"
#include "XA_file.h"
#include "XA_fft2.h"
#include "XA_spln2.h"
#include "XA_spln3.h"
#include "XA_iwfr.h"
#include "XA_tie.h"
#include "XA_tiff.h"
#include "XA_DICOM.h"
#include "XA_move2.h"
#include "fftwd3frc.h"

#include "FindPeaks.h" // a function for finding peak values in a 3D distribution
#include "UTR.h" // various include headers, constants and templates

using namespace xar;


string ReadString(FILE* ff0)
{
	char cline[1024], ctitle[1024], cparam[1024];
	fgets(cline, 1024, ff0); strtok(cline, "\n");
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading input text file.");
	return string(cparam);
}


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

		bool bXrayPar(false);
		index_t iDotPosition = sInputParamFile.find_last_of(".");
		if (iDotPosition == string::npos) throw std::runtime_error("Error: filename extension not found in the command-line argument input parameter file");
		string sInputParamFileExt = sInputParamFile.substr(iDotPosition, 4);
		if (sInputParamFileExt == string(".xri") || sInputParamFileExt == string(".XRI")) bXrayPar = true;
		else if (sInputParamFileExt == string(".txt") || sInputParamFileExt == string(".TXT")) bXrayPar = false;
		else throw std::runtime_error("Error: unknown filename extension in the command-line argument input parameter file");

		FILE* ff0 = fopen(sInputParamFile.c_str(), "rt");
		if (!ff0) throw std::runtime_error(string("Error: cannot open parameter file " + sInputParamFile + ".").c_str());

		// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
		printf("\nReading input parameter file %s ...", sInputParamFile.c_str());
		while (true)
		{
			fgets(cline, 1024, ff0);
			if (!(cline[0] == '/' && cline[1] == '/')) break;
		}

		int iModality(1);
		if (!bXrayPar)
		{
			strtok(cline, "\n"); // 1. Modality and units of length(UL): TEM and Angstroms(0) or X-ray CT and microns(1)
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
				throw std::runtime_error("Error reading units of length from input parameter file.");
			iModality = atoi(cparam);
		}
		if (iModality == 0) // electrons and Angstroms (TEM)
			printf("\nUnits of length (UL) are set to Angstroms (TEM imaging mode will be used)");
		else if (iModality == 1) // X-rays and microns (X-ray CT)
			printf("\nUnits of length (UL) are set to microns (X-ray CT imaging mode will be used)");
		else
			throw std::runtime_error("Error: modality and units of length parameter in input parameter file must be 0 for TEM/Angstroms or 1 for X-ray CT/microns.");


		if (!bXrayPar) fgets(cline, 1024, ff0);
		strtok(cline, "\n"); // 2. Verbose output during execution? Yes = 1, No = 0
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading verbose output parameter from input parameter file.");
		bool bVerboseOutput(true); // if this is TRUE, additional information is printed during execution
		(atoi(cparam) == 0 || atoi(cparam) == 1) ? bVerboseOutput = (bool)atoi(cparam) : throw std::runtime_error("Error: verbose output parameter must be 0 or 1 in input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3. Input file with rotation angles and defocus distances
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading file name with rotation angles and defocus distances from input parameter file.");
		string defocfile(cparam);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4.Optional file with indexes and weights for subseries selection or NONE
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
		{
			ReadDefocusParamsFile(defocfile, v2angles, vvdefocus, bVerboseOutput);
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
					vvdefocusSel.push_back(vvdefocus[i]);
				}
				v2angles = v2anglesSel;
				vvdefocus = vvdefocusSel;
				if (dTotalWeight == 0) throw std::runtime_error("Error: the sum of weights in the second column of the subseries selection file is equal to zero.");
				printf("\nThere are %zd entries in the subseries selection file.", v2angles.size());
			}
		}
		else if (GetFileExtension(defocfile) == string(".RELIONNEW"))
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
		bool bRotZonly(true), bRotYonly(true), bRotZ2only(true);
		vndefocus[0] = vvdefocus[0].size();
		ndefocusmin = ndefocusmax = vndefocus[0];
		for (index_t i = 1; i < nangles; i++)
		{
			vndefocus[i] = vvdefocus[i].size();
			if (vndefocus[i] < ndefocusmin) ndefocusmin = vndefocus[i];
			if (vndefocus[i] > ndefocusmax) ndefocusmax = vndefocus[i];
			if (v2angles[i].a != 0) { bRotYonly = false; bRotZ2only = false; }
			if (v2angles[i].b != 0) { bRotZonly = false; ; bRotZ2only = false; }
			for (index_t j = 0; j < vndefocus[i]; j++) if (vvdefocus[i][j].a != 0) { bRotYonly = false; bRotZonly = false; }
		}
		if (ndefocusmax > 1)
			throw std::runtime_error("Error: this program is currently implemented only for defocus parameter files with one defocus distance per orientation.");
		if (bRotZonly) printf("\nDefocus parameter file contains orientations with rotations around the Z axis only.");
		if (bRotYonly) printf("\nDefocus parameter file contains orientations with rotations around the Y axis only.");
		if (bRotZ2only) printf("\nDefocus parameter file contains orientations with rotations around the Z'' axis only.");

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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5. Input filename base of defocus series of the sample in TIFF, GRD or RAW format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading defocus series file name base from input parameter file.");
		string filenamebaseIn = cparam;
		bool bRAWinput(false), bTIFFinput(false), bGRDinput(false);
		string strTemp = GetFileExtension(filenamebaseIn);
		if (strTemp == string(".TIFF") || strTemp == string(".TIF")) bTIFFinput = true;
		else if (strTemp == string(".RAW")) bRAWinput = true;
		else if (strTemp == string(".GRD")) bGRDinput = true; // this value is not used below at the moment, as GRD is considered to be a default when all other input formats are false
		else throw std::runtime_error("Error: input filename extension must be TIF, GRD or RAW.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6. RAW file parameters: HeaderLength(bytes) Endianness(0 = little, 1 = big) ElementLength(bytes) 
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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7. Input data normalization factors f1 f2 (input -> (input / f1) + f2)
		double dNormFactor1(1.0), dNormFactor2(0.0);
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::runtime_error("Error reading input data normalization factors from input parameter file.");
		dNormFactor1 = atof(cparam);
		dNormFactor2 = atof(cparam1);
		printf("\nInput data normalization factors: f1 = %g, f2 = %g", dNormFactor1, dNormFactor2);
		if (dNormFactor1 == 0) throw std::runtime_error("The first normalization factor cannot be zero (leads to division by zero).");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8. Width(x, pix) and height(y, pix) of input images
		// for RAW files, these parameter determine how the images are read; for TIFF and GRD files, these parameters may be used for trimming the frames
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
			throw std::runtime_error("Error reading the width and height of input images in pixels from input parameter file.");
		index_t nx = atoi(cparam);
		index_t ny = atoi(cparam1);
		printf("\nDimensions of input images (possibly, after trimming or padding): width = %zd, height = %zd (pixels)", nx, ny);
		int nx2 = int(nx) - 2, ny2 = int(ny) - 2;
		index_t nxny = (nx * ny);
		int nxd2 = int(nx / 2), nyd2 = int(ny / 2);
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9. Pixel size in UL
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
			throw std::runtime_error("Error reading pixel size from input parameter file.");
		double xst = atof(cparam);
		if (xst <= 0)
			throw std::runtime_error("Error: pixel size must be positive in input parameter file.");
		double yst = xst; // this is a provision for possible future extension to non-square pixels
		double xlo(-0.5 * xst * nx);
		double xhi = xlo + xst * nx;
		double ylo(-0.5 * yst * ny);
		double yhi = ylo + yst * ny;
		printf("\nPixel size = %g (UL)", xst);
		printf("\nPhysical boundaries of input images: Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g (UL)", xlo, xhi, ylo, yhi);
		// NOTE that the input 2D images will be trimmed or padded to the (ny, nx) size, if they don't have this size already
		double fxst = 1.0 / (xhi - xlo), fyst = 1.0 / (yhi - ylo);
		double fxlo = -fxst * nxd2, fylo = -fyst * nyd2;

		double zlo(xlo); // minimum output defocus in UL (relative to the centre of the reconstruction volume)
		double zhi(xhi); // maximum output defocus in UL (relative to the centre of the reconstruction volume)
		double zst(xst); // this is a provision for possible future extensions to non-qubic voxels
		double dzextra(0); // defocus offset parameter
		if (!bXrayPar)
		{
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10. Output defocus distances min and max in UL
			if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::runtime_error("Error reading output defocus distances from input parameter file.");
			zlo = atof(cparam); // minimum output defocus in UL
			zhi = atof(cparam1); // maximum output defocus in UL
			dzextra = atof(cparam2);
			if (zlo > zhi) std::swap(zlo, zhi);
		}
		index_t nz = index_t((zhi - zlo) / zst + 0.5); // number of defocus planes to propagate to
		printf("\nDimensions of the output volume: width = %zd, height = %zd, thickness = %zd (pixels)", nx, ny, nz);
		printf("\nPhysical boundaries of the output volume: Xmin = %g, Xmax = %g, Ymin = %g, Ymax = %g, Zmin = %g, Zmax = %g (UL)", xlo, xhi, ylo, yhi, zlo, zhi);
		printf("\nDefocus offset = %g (UL)", dzextra);
		if (nz <= 0)
			throw std::runtime_error("Error: number of z steps is not positive.");
		int nzd2 = int(nz / 2), nz2 = int(nz) - 2;
		double fzst = 1.0 / (zhi - zlo);
		double fzlo = -fzst * nzd2;
		zoutAver -= dzextra;

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
		index_t ic = index_t((xc - xlo) / xst + 0.5), jc = index_t((yc - ylo) / yst + 0.5), kc = index_t((zc - zlo) / zst + 0.5);

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

		double aobj(0);
		if (!bXrayPar)
		{
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13. Objective aperture in mrad
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading objective aperture from input parameter file.");
			aobj = atof(cparam);
			if (aobj != 0) printf("\nObjective aperture (half-angle) = %g (mrad)", aobj);
			else  printf("\nObjective aperture is infinite.");
			if (aobj < 0 || aobj > 1000)
				throw std::runtime_error("Error: objective aperture value appears to be wrong.");
		}
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth

		double Cs3(0), Cs5(0);
		if (!bXrayPar)
		{
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14. Spherical aberrations Cs3 and Cs5 in mm
			if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading spherical aberrations from input parameter file.");
			Cs3 = atof(cparam);
			Cs5 = atof(cparam1);
			printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
			Cs3 *= 1.e+7; // mm --> UL
			Cs5 *= 1.e+7; // mm --> UL
		}
		
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

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17. Absorption fraction beta / delta (or -delta / beta in the case of .xri parameter file)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading absorption coefficient from input parameter file.");
		double asigma = atof(cparam);
		if (bXrayPar)
		{
			printf("\ndelta / beta = %g (n = 1 - delta + i * beta)", asigma);
			if (asigma <= 0) throw std::runtime_error("Error: delta / beta must be positive in the case of .xri parameter file (X-ray imaging mode).");
			asigma = -1.0 / asigma;
		}
		else printf("\nAbsorption fraction beta/delta = %g", asigma);
		if (iModality == 0 && asigma < 0) throw std::runtime_error("Error: beta/delta must be non-negative in TEM mode.");
		if (iModality == 1 && asigma >= 0) throw std::runtime_error("Error: beta/delta must be negative in hard X-ray mode (unless .xri parameter file is used).");
		if (asigma >= 0 && (iCTFcorrectionMode == 1 || iCTFcorrectionMode == 2)) 
			throw std::runtime_error("Error: beta/delta must be negative in TIE-Hom phase retrieval.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 18. Average noise to signal ratio (1/SNR) in input images (1/sqrt(Nphot))
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading average noise-to-signal ratio in input images from input parameter file.");
		double dInputNSR = atof(cparam); // average NSR in the input images
		if (dInputNSR < 0)
			throw std::runtime_error("Error: NSR in input images must be non-negative.");
		printf("\nAverage NSR in input images = %g", dInputNSR);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19. Tikhonov regularization parameter
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading minimal phase reconstruction error/regulariztion parameter from input parameter file.");
		double epsilon = atof(cparam);
		printf("\nTikhonov regularization epsilon parameter = %g", epsilon);
		if (epsilon < 0)
			throw std::runtime_error("Error: the regulalization parameter must be non-negative.");
		if (epsilon == 0 && asigma == 0)
			throw std::runtime_error("Error: when beta/delta == 0, epsilon must be positive.");

		int imodeSymmetry(0);
		index_t nanglesSym{ 1 }, nanglesIn1Sym{ nangles };
		vector<Pair> v2anglesSym;
		vector<vector <Pair> > vvdefocusSym;
		if (!bXrayPar)
		{
			fgets(cline, 1024, ff0); strtok(cline, "\n"); //20. Enforce symmetry: not_apply(0), post-apply(1), or distribute input orientations(2)
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading enforcing symmetry mode from input parameter file.");
			imodeSymmetry = atoi(cparam);
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
		}

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

		int imodePeaks(0);
		double datomsizeXY, datomsizeZ;
		if (!bXrayPar)
		{
			fgets(cline, 1024, ff0); strtok(cline, "\n"); // 24. Peak localization mode
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading peak localization mode from input parameter file.");
			imodePeaks = atoi(cparam);
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
			if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading transverse and longitudinal side lengths for peak localization from input parameter file.");
			datomsizeXY = atof(cparam);
			datomsizeZ = atof(cparam1);
			printf("\nTransverse and longitudinal side lengths for peak localization = %g %g (UL)", datomsizeXY, datomsizeZ);
			if (imodePeaks && (int(datomsizeXY / xst + 0.5) < 2 || int(datomsizeZ / zst + 0.5) < 2))
				throw std::runtime_error("Error: transverse and longitudinal side lengths for peak localization must be 2 x z_step or larger.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 26. Output file name base in GRD, TIFF or DICOM format
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading output file name base from input parameter file.");
		string filenamebaseOut = cparam;
		printf("\nFile name base for 3D potential or beta = %s", filenamebaseOut.c_str());
		strTemp = GetPathFromFilename(filenamebaseOut, false);
		if (!std::filesystem::exists(strTemp))
			throw std::runtime_error("Error: the specified file folder for 3D potential or beta does not seem to exist.");
		bool bTIFFoutput(false), bDICOMoutput(false);
		std::map<string, string> mTags; // DICOM tags
		if (GetFileExtension(filenamebaseOut) == string(".TIFF") || GetFileExtension(filenamebaseOut) == string(".TIF")) bTIFFoutput = true;
		else if (GetFileExtension(filenamebaseOut) == string(".DCM"))
		{
			bDICOMoutput = true;
			FILE* ff0 = fopen("UTR_DICOM.txt", "rt"); // default text file containing pre-defined DICOM tags
			if (!ff0)
			{
				printf("\nWARNING: unable to open UTR_DICOM.txt file; using default values for DICOM tags!");
				mTags.insert({ "PatientName", string("UTR") });
				mTags.insert({ "PatientID", string("UTR") });
				mTags.insert({ "PatientBirthDate", string("19800101") });
				mTags.insert({ "PatientSex", string("Female") });
				mTags.insert({ "StudyID", string("UTR") });
				mTags.insert({ "AccessionNumber", string("UTR_coronal") });
				mTags.insert({ "StudyInstanceUID", mTags["StudyID"] });
				mTags.insert({ "SeriesInstanceUID", mTags["AccessionNumber"] });
				mTags.insert({ "NominalScannedPixelSpacing", "0.1" });
				mTags.insert({ "WindowWidth", "4095" });
				mTags.insert({ "WindowCenter", "2047" });
			}
			else
			{
				printf("\nReading UTR_DICOM.txt file ...");
				// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
				printf("\nReading input parameter file %s ...", sInputParamFile.c_str());
				while (true)
				{
					fgets(cline, 1024, ff0);
					if (!(cline[0] == '/' && cline[1] == '/')) break;
				}
				strtok(cline, "\n");
				if (sscanf(cline, "%s %s", ctitle, cparam) != 2) 
					throw std::runtime_error("Error reading input text file.");
				mTags.insert({ "PatientName", string(cparam)});
				mTags.insert({ "PatientID", ReadString(ff0) });
				mTags.insert({ "PatientBirthDate", ReadString(ff0) });
				mTags.insert({ "PatientSex", ReadString(ff0) });
				mTags.insert({ "StudyID", ReadString(ff0) });
				mTags.insert({ "AccessionNumber", ReadString(ff0) });
				mTags.insert({ "StudyInstanceUID", ReadString(ff0) });
				mTags.insert({ "SeriesInstanceUID", ReadString(ff0) });
				mTags.insert({ "NominalScannedPixelSpacing", ReadString(ff0) });
				mTags.insert({ "WindowWidth", ReadString(ff0) });
				mTags.insert({ "WindowCenter", ReadString(ff0) });
				for (auto it = mTags.begin(); it != mTags.end(); ++it)
					printf("\n %s %s.", (it->first).c_str(), (it->second).c_str());
			}
			fclose(ff0);
		}
		else if (GetFileExtension(filenamebaseOut) != string(".GRD"))
			throw std::runtime_error("Error: output filename extension must be TIF/TIFF, DCM or GRD.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 27. Minimum-in, maximum-in and maximum-out for 16-bit output (all zeros trigger 32-bit_output)
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::runtime_error("Error reading minimum-in, maximum-in and maximum-out for 16-bit output from input parameter file.");
		float fbit16_minin = (float)atof(cparam);
		float fbit16_maxin = (float)atof(cparam1);
		long lbit16_maxout = atoi(cparam2);
		bool b16bitout = (fbit16_minin == 0 && fbit16_maxin == 0 && lbit16_maxout == 0) ? false : true;
		if (bTIFFoutput || bDICOMoutput)
		{
			if (b16bitout)
			{
				printf("\nMin-in = %g, max-in = %g and max_out = %ld (for 16-bit file output).", fbit16_minin, fbit16_maxin, lbit16_maxout);
				if (fbit16_minin >= fbit16_maxin || lbit16_maxout <= 0 || lbit16_maxout > 65535)
					throw std::runtime_error("Unsuitable values for minimum-in, maximum-in or maximum-out for 16-bit output from input parameter file.");
			}
			else 
				if (bTIFFoutput) printf("\n32-bit (floating point) TIFF output file format will be used.");
				else throw std::runtime_error("only 16-bit data can be saved in DICOM files.");
		}
		else // GRD output
			if (b16bitout) throw std::runtime_error("16-bit output cannot be saved in GRD files.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 28. Euler rotation angles (Z, Y',Z") in degrees for 3D output
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::runtime_error("Error reading Euler rotation angles for 3D output from input parameter file.");
		double angleZout = -atof(cparam); // we multiply by (-1) here to make it consistent with the behaviour of ReadDefocusParameterFile(), ReadRelionParameterFile() and pdb.exe
		double angleY1out = -atof(cparam1); // we multiply by (-1) here to make it consistent with the behaviour of ReadDefocusParameterFile(), ReadRelionParameterFile() and pdb.exe
		double angleZ2out = -atof(cparam2); // we multiply by (-1) here to make it consistent with the behaviour of ReadDefocusParameterFile(), ReadRelionParameterFile() and pdb.exe
		bool b3Drotout = (angleZout == 0 && angleY1out == 0 && angleZ2out == 0) ? false : true;
		if (b3Drotout)
			printf("\nReconstructed 3D potential or beta will be rotated by (%g, %g, %g) degrees for output.", angleZout, angleY1out, angleZ2out);
		else
			printf("\nReconstructed 3D potential or beta will NOT be rotated for output.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 29.Slice output 3D volume in XY(0), XZ(1), YZ(2), YX(3), ZX(4) or ZY(5) 2D-planes
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading reslicing parameter from input parameter file.");
		int iResliceMode = atoi(cparam);
		if (iResliceMode == 0) printf("\nOutput 3D volume will be sliced in XY planes");
		else if (iResliceMode == 1) printf("\nOutput 3D volume will be sliced in XZ planes");
		else if (iResliceMode == 2) printf("\nOutput 3D volume will be sliced in YZ planes");
		else if (iResliceMode == 3) printf("\nOutput 3D volume will be sliced in YX planes");
		else if (iResliceMode == 4) printf("\nOutput 3D volume will be sliced in ZX planes");
		else if (iResliceMode == 5) printf("\nOutput 3D volume will be sliced in ZY planes");
		else throw std::runtime_error("Error: unknown value for the output 3D volume reslicing parameter.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 30.Flip direction of X, Y, and / or Z axes in output 3D volume
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4) throw std::runtime_error("Error reading axes flipping parameters for 3D output from input parameter file.");
		bool bXFlip = atoi(cparam) == 0 ? false : (atoi(cparam) == 1 ? true : throw std::runtime_error("Unknown value of X axis flipping parameter in input parameter file."));
		bool bYFlip = atoi(cparam1) == 0 ? false : (atoi(cparam1) == 1 ? true : throw std::runtime_error("Unknown value of Y axis flipping parameter in input parameter file."));
		bool bZFlip = atoi(cparam2) == 0 ? false : (atoi(cparam2) == 1 ? true : throw std::runtime_error("Unknown value of Z axis flipping parameter in input parameter file."));

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 31. Make thick output slices: 0=no, 1=MIP, 2=averaging; slice thickness in UL
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading thick output slices info from input parameter file.");
		int iThickSlicesOut = atoi(cparam);
		int kThickStep(0);
		double dOutSliceThickness = atof(cparam1);
		if (iThickSlicesOut)
		{
			kThickStep = (int)(dOutSliceThickness / zst + 0.5);
			if (kThickStep <= 1) throw std::runtime_error("Slice thickness is smaller than or equal to the voxel thickness.");
		}
		switch (iThickSlicesOut)
		{
		case 0: printf("\nOutput slices will have single-voxel thickness."); 
			break;
		case 1: 
			printf("\nThick output slices will be obtained using MIP over %d voxels.", kThickStep); 
			break;
		case 2: 
			printf("\nThick output slices will be obtained by averaging over %d voxels.", kThickStep); 
			break;
		default: 
			throw std::runtime_error("Error: unknown value for thick output slices mode in input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 32. Save or not the sampling matrix in files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading save_or_not inverse sampling matrix switch from input parameter file.");
		int nSaveDefocCAmpsOrSamplingMatrix = atoi(cparam);
		if (nSaveDefocCAmpsOrSamplingMatrix == 1)
			printf("\n3D sampling matrix will be saved in GRD files");
		else
			printf("\nSampling matrix will not be saved in files");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 33. Import and reprocess existing 3D potential or beta files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading import and reprocess existing 3D potential or beta files switch from input parameter file.");
		int imodeReprocess = atoi(cparam);
		switch (imodeReprocess)
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

		if (imodeReprocess == 1) // only read in pre-caculated 3D potential or beta from "output" files and reprocess them
		{
			bGRDinput = bTIFFinput = bRAWinput = false;
			string strTemp = GetFileExtension(filenamebaseOut);
			if (strTemp == string(".TIFF") || strTemp == string(".TIF")) bTIFFinput = true;
			else if (strTemp == string(".GRD")) bGRDinput = true;
			else if (strTemp == string(".RAW")) bRAWinput = true;
			else
				if (strTemp == string(".DCM"))
				{
					bDICOMoutput = true;
					printf("\nWARNING: input filename extension (in this mode - it is taken from the output filename template) must be TIF, GRD or RAW.");
					printf("\nPress 't' if you want to substitute TIFF extension for input files, 'g' for GRD, 'r' for RAW or any other symbol to exit: ");
					int c = getchar(); 
					if (c == 't' || c == 'T') { bTIFFinput = true; filenamebaseOut.replace(filenamebaseOut.find_last_of("."), filenamebaseOut.length() - filenamebaseOut.find_last_of("."), ".tif"); }
					else if (c == 'g' || c == 'G') { bGRDinput = true; filenamebaseOut.replace(filenamebaseOut.find_last_of("."), filenamebaseOut.length() - filenamebaseOut.find_last_of("."), ".grd"); }
					else if (c == 'r' || c == 'R') { bRAWinput = true; filenamebaseOut.replace(filenamebaseOut.find_last_of("."), filenamebaseOut.length() - filenamebaseOut.find_last_of("."), ".raw"); }
					else exit(1);
				}
			else throw std::runtime_error("Error: input filename extension (in this mode - it is taken from the output filename template) must be TIF, GRD, DCM or RAW.");
		}
		else
		{
			if (nx != 2 * nxd2 || ny != 2 * nyd2)
				throw std::runtime_error("Error: width and height parameters in input parameter file must be even.");
			if (nz != 2 * nzd2)
				throw std::runtime_error("Error: number of z points must be even.");
		}

		if (imodeReprocess != 1 && (iCTFcorrectionMode == 2 || iCTFcorrectionMode == 4) && (!bSingleDefocusDistance || iAstigmatism != 0))
		{
			printf("\n\nWARNING: different defocus distances and/or non-zero astigmatism have been detected in the input file %s.", defocfile.c_str());
			printf("\n3D CTF correction have only been implemented for constant defocus distance and no astigmatism, possibly leading to inaccurate results in this case.");
			printf("\nDo you still want to continue (y = yes, any other symbol = no): ");
			int c = getchar();
			if (c != 'y' && c != 'Y') exit(1);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // // 34. Folder name for auxiliary files
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading folder name for auxiliary files from input parameter file.");
		string folderAux = cparam;
		printf("\nFolder for auxiliary file output = %s", folderAux.c_str());
		if (folderAux.rfind('\\') != folderAux.size() - 1 && folderAux.rfind('/') != folderAux.size() - 1) // the last character is not '\' and not '/'
			folderAux.append("/");
		if (!std::filesystem::exists(folderAux))
			throw std::runtime_error("Error: the specified auxiliary file folder does not seem to exist.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 35. Number of parallel threads
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
		if (imodeReprocess) // only read in pre-caculated 3D potential or beta from "output" files and reprocess them
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

		string filenamebaseOutSMatrix("BAD_STRING"); // don't use it, unless it is redefined later
		vector<string> voutfilenamesTotDefocCAmpOrFM; // output filenames for sampling matrix
		if (nSaveDefocCAmpsOrSamplingMatrix == 1) // create filenames for saving the sampling matrix
		{
			std::filesystem::path apath(filenamebaseOut);
			filenamebaseOutSMatrix = apath.filename().string();
			filenamebaseOutSMatrix.replace(filenamebaseOutSMatrix.find_last_of("."), filenamebaseOutSMatrix.length() - filenamebaseOutSMatrix.find_last_of("."), "FM.grd");
			filenamebaseOutSMatrix = folderAux + filenamebaseOutSMatrix;
			FileNames(nz, 1, filenamebaseOutSMatrix, voutfilenamesTotDefocCAmpOrFM);

			printf("\nOutput file name base for sampling matrix = %s", filenamebaseOutSMatrix.c_str());
		}

		//************************************ end creating vectors of input and output file names


		//*********************************** start main calculations
		bool bAbort(false);
		bool bPadded(false); // this indicates if input images have been padded, and so the output images need to be trimmed before saving to files
		index_t nx0(0), ny0(0); // these numbers will store the actual dimensions of input images, in case they need to be restored later

		std::unique_ptr<IXAHead> pHead(new Wavehead2D(wl, ylo, yhi, xlo, xhi)); // any header data in input image files will be ignored

		printf("\n\nPreparing FFTW library ...");
		if (!fftw_init_threads())
			std::runtime_error("Error: multi-threaded initialization of FFTW library failed.");
		// single Fftwd2c object is used for repeated 2D FFTs using FFTW later (the plan is created at this point and reused later)
		Fftwd2fc fftw2Df((int)ny, (int)nx, 1, true); // we don't use mutlithreaded 2D FFT here, as the code calling these functions below is already mutlithreaded

		XArray3D<float> K3out; // big 3D reconstructed array (will be allocated after Samp3 matrix is truncated)

		if (!imodeReprocess) // do phase retrieval and backpropagation prior to 3D filtering and output
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
					fcomplex Fxc, Fyc, Fzc;
					if (bRotCentreShift)
					{
						Fxc = fcomplex(0, -1) * float(tPI * (dxc * (1.0 - cosangleY * cosangleZ) - dyc * sinangleZ + dzc * sinangleY * cosangleZ));
						Fyc = fcomplex(0, -1) * float(tPI * (dxc * cosangleY * sinangleZ + dyc * (1.0 - cosangleZ) - dzc * sinangleY * sinangleZ));
						Fzc = fcomplex(0, -1) * float(tPI * (-dxc * sinangleY + dzc * (1.0 - cosangleY)));
					}

					// start the execution timer per angle
					std::chrono::system_clock::time_point start_timeA;
					if (bVerboseOutput) start_timeA  = std::chrono::system_clock::now();

					//index_t ndefocus = vndefocus[na]; // number of defocus planes at the current illumination angle
					// NOTE: in the current version of the program, ndefocus is always 1
					vector<Pair> vdefocus = vvdefocus[na]; // vector of input Z" rotation angles and defocus positions at the current illumination angle
					vector<string> vinfilenames = vvinfilenames[na]; // input filenames of defocused images at the current illumination angle
					vector<string> voutfilenames(nz);

					XArray2D<float> int0; // input defocused intensity image
					double zout = vdefocus[0].b - dzextra; // defocus distance

					// read defocused images from files
					if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of file read time
					//if (bVerboseOutput) printf("\nz'' rotation angle = %g (degrees), defocus distance = %g (UL)", vdefocus[0].a, vdefocus[0].b);
					//if (bVerboseOutput) printf("\nReading input file %s ...", vinfilenames[0].c_str());
					if (bRAWinput) XArData::ReadFileRAW(int0, vinfilenames[0].c_str(), ny, nx, nHeaderLength, nElementLength, bBigEndian);
					else if (bTIFFinput) TIFFReadFile(int0, vinfilenames[0].c_str()); // read input TIFF files
					else XArData::ReadFileGRD(int0, vinfilenames[0].c_str(), wl); //	read input GRD files
					if (bVerboseOutput) liFileReadTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();

					if (ny0 == 0) ny0 = int0.GetDim1();
					else if (ny0 != int0.GetDim1()) throw std::runtime_error("Error: y-dimension in the input file is different from that in previous input files.");
					if (nx0 == 0) nx0 = int0.GetDim2();
					else if (nx0 != int0.GetDim2()) throw std::runtime_error("Error: x-dimension in the input file is different from that in previous input files.");

					if (nx0 != nx || ny0 != ny)
					{
						if (nx0 < nx && ny0 <= ny || ny0 < ny && nx0 <= nx)
						{
							//printf("\nWARNING: 2D array from the image file will be padded with 1.0 to the array dimensions (Width = %zd, Height = %zd).", nx, ny);
							XArray2DMove<float> tmp2(int0);
							tmp2.Pad((ny - ny0) / 2, ny - ny0 - (ny - ny0) / 2, (nx - nx0) / 2, nx - nx0 - (nx - nx0) / 2, 1.0f);
							bPadded = true;
							// NOTE that the reconstructed 3D volume WILL be trimmed back before saving
						}
						else if (nx0 > nx && ny0 >= ny || ny0 > ny && nx0 >= nx)
						{
							//printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions (Width = %zd, Height = %zd).", nx, ny);
							XArray2DMove<float> tmp2(int0);
							tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
							// NOTE that the reconstructed 3D volume will NOT be padded back before saving, as it does not seem to make sense
						}
						else
							throw std::runtime_error("Error: the dimensions of the 2D array in the image file are inconsistent with the array dimensions in input parameter file (both dimensions should be equal, smaller or larger simultaneously).");
					}
					int0.SetHeadPtr(pHead->Clone());
					// NOTE that after the above trimming or padding, nx, ny, xlo, xhi, ylo, yhi parameters will correspond to those defined after parameter 7 above

					// renormalize input data
					if (dNormFactor1 != 1.0) int0 /= float(dNormFactor1);
					if (dNormFactor2 != 0.0) int0 += float(dNormFactor2);
					if (int0.Norm(eNormMin) < 0)
					{
						//printf("\nWARNING: negative values in the input intensity file.");
						int0.ThresholdLow(0.0f, 0.0f);
					}

					// rotate input defocused image around Z'' back to zero angle
					if (vdefocus[0].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
					{
						float aver = int0.NormAverEdge(5);
						// shift along X and/or Y back to the unshifted position
						if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
						{
							XArray2DMove<float> xamove(int0);
							xamove.Move((long)floor(-v2shifts[na].b / yst + 0.5), (long)floor(-v2shifts[na].a / xst + 0.5), aver);
						}
						// rotate input defocused complex amplitude around Z'' back to zero angle
						if (vdefocus[0].a != 0)
						{
							XArray2D<float> vintTemp(int0);
							XArray2DSpln<float> xaSpln(vintTemp);
							xaSpln.Rotate(int0, -vdefocus[0].a, yc, xc, aver); // expecting uniform background
						}
					}

					// now do 2D CTF correction
					XArray2D<fcomplex> K2four; // symmetrized backpropagated contrast in the reciprocal space

					if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now();  // start of 2D FFT time

					XA_IWFR<float> xa_iwfr;
					switch (iCTFcorrectionMode)
					{
					case 0: // no CTF correction
					case 2: // 3D TIE-Hom correction later
					case 4: // 3D CTF correction later
						//if (bVerboseOutput && iCTFcorrectionMode != 1) printf("\nCalculating 2D FFT of the input image ...");
						int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0 (after the flat field correction)
						K2four = MakeComplex(int0, 0.0f, false);
						K2four.Shuffle();
						fftw2Df.ForwardFFT(K2four);
						K2four.Shuffle();
						break;
					case 1: // 2D TIE-Hom correction
						if (zout != 0) // if zout == 0, no TIE-Hom correction is applied (corresponding to the conventional=absorption CT case)
						{
							//if (bVerboseOutput)	printf("\nApplying 2D TIE-Hom correction to the input image ...");
							int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0 (after the flat field correction)
							K2four = MakeComplex(int0, 0.0f, false);
							//normalization of the inverse 2D minus-Laplacian here corresponds to the normalization of the InvertCTF_DT1() function below
							double alphatemp = -4.0 * PI * asigma / (zoutAver * wl); // asigma should be negative, so alpha is positive
							double normtemp = -alphatemp * sqrt(1.0 + asigma * asigma) / asigma;
							fftw2Df.InverseMLaplacian1(K2four, alphatemp, normtemp);
						}
						break;
					case 3: // 2D CTF correction
						//if (bVerboseOutput) printf("\nApplying 2D CTF correction to the input image ...");
						int0.Log0(); // convert image intensity into ln(I / I_in), assuming I_in = 1.0 (after the flat field correction)
						if (bRelion)
							xa_iwfr.InvertCTF_DT2(int0, K2four, fftw2Df, zout, vastigm[na].a, vastigm[na].b * PI180, k2maxo, Cs3, Cs5, asigma, epsilon, bESCC, 1);
						else
							xa_iwfr.InvertCTF_DT2(int0, K2four, fftw2Df, zout, 0, 0, k2maxo, Cs3, Cs5, asigma, epsilon, bESCC, 1);
						// the output (K2four) is normalized as delta * 4 * PI * sqrt(1 + sigma^2) / wl / (xst * yst), the (xst * yst) term here is due to 2D DFT
						break;
					}
					if (bSelectFrames) K2four *= sfWeights[na];

					if (bVerboseOutput) li2DFFTtime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();


					// Update 3D object at the current illumination direction on the Ewald sphere in 3D Fourier space
					if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of back-propagation time
					//if (bVerboseOutput) printf("\nUpdating 3D reconstructed object on the Ewald sphere ...");
					CT_3Dgridding<float>(K2four, V3, Samp3, angleY, angleZ, float(fxlo), float(fxst), float(fylo), float(fyst), float(fzlo), float(fzst), Fxc, Fyc, Fzc, wl, bRotCentreShift, bESCC);
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
			dInputNSR *= sqrt(nx * ny); // change of the STD of noise after the 2D DFT; later it will be multiplied by sqrt(Samp3)
			if (dInputNSR == 0) // the case of noise-free input images
			{
				#pragma omp parallel for shared (V3, Samp3)
				for (int k = 0; k < nz; k++)
				{
					float x;
					for (int j = 0; j < ny; j++)
						for (int i = 0; i < nx; i++)
						{
							x = Samp3[k][j][i];
							if (x == 0) V3[k][j][i] = 0.0f;
							else V3[k][j][i] *= (1 / x);
						}
				}
			}
			else // the case of noisy input images
			{
				switch (iCTFcorrectionMode)
				{
				case 0: // noise filtering in the case of no CTF correction
				case 2: // noise filtering in the case of 3D TIE-Hom correction later
				case 3: // noise filtering in the case of 2D CTF correction - somehow this version of noise filtering seems to be working the best for this case
				case 4: // noise filtering in the case of 3D CTF correction later
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
								y = abs(V3[k][j][i]) / float(dInputNSR * sqrt(x)); // SNR in this Fourier coefficient
								V3[k][j][i] *= (1 / x) * (0.5f + atan(100.0f * (y - fTargetSNR)) / PIf); // suppress V3 elements with SNR < TargetSNR
								//if (y < fTargetSNR) V3[k][j][i] = 0.0f; // "hard" cut-off filter
							}
					}
					break;
				// tests appear to indicate that noise filtering using the simple algorithm above works very similar to the following algortithm for case 1
				// with very small epsilon. If epsilon is increased below, the filtering actually becomes less severe (i.e. the output images become noisier)
				case 1: // noise filtering in the case of 2D TIE-Hom correction
					{
						double norm = sqrt(1.0 + asigma * asigma) * dInputNSR;
						double fac2 = xar::PI * wl * zoutAver;
						double dksi2 = fac2 / ((xhi - xlo) * (xhi - xlo));
						double deta2 = fac2 / ((yhi - ylo) * (yhi - ylo));
						double dzeta2 = fac2 / ((zhi - zlo) * (zhi - zlo));

						#pragma omp parallel for shared (V3, Samp3)
						for (int k = 0; k < nz; k++)
						{
							int k1, j1, i1;
							double zeta2, eta2;
							k1 = k - nzd2;
							zeta2 = k1 * k1 * dzeta2 - asigma;
							float x, y, fInputNSR;
							for (int j = 0; j < ny; j++)
							{
								j1 = j - nyd2;
								eta2 = j1 * j1 * deta2 + zeta2;
								for (int i = 0; i < nx; i++)
								{
									x = Samp3[k][j][i];
									if (x == 0) { V3[k][j][i] = 0.0f; continue; }

									i1 = i - nxd2;
									fInputNSR = float(norm * sqrt(x) / (dksi2 * i1 * i1 + eta2 + epsilon)); // dInputNSR is contained in the "norm" term

									y = abs(V3[k][j][i]) / fInputNSR; // SNR in this Fourier coefficient
									V3[k][j][i] *= (1 / x) * (0.5f + atan(100.0f * (y - fTargetSNR)) / PIf); // suppress V3 elements with SNR < TargetSNR
									//if (y < fTargetSNR) V3[k][j][i] = 0.0f; // "hard" cut-off filter
								}
							}
						}
					}
					break;
				/*
				// tests appear to indicate that noise filtering using the simple algorithm above (case 0:) works better than the following more accurate algorithm
				case 3: // noise filtering in the case of 2D CTF correction
					{
						double omega = atan(asigma);
						double fac2 = xar::PI * wl * zoutAver;
						double dksi2 = fac2 / ((xhi - xlo) * (xhi - xlo));
						double deta2 = fac2 / ((yhi - ylo) * (yhi - ylo));
						double dzeta2 = fac2 / ((zhi - zlo) * (zhi - zlo));

						#pragma omp parallel for shared (V3, Samp3)
						for (int k = 0; k < nz; k++)
						{
							int j1, i1;
							int k1 = k - nzd2;
							float x, y, fInputNSR;
							double eta2, sintemp, sin2;
							double zeta2 = k1 * k1 * dzeta2 - omega;
							for (int j = 0; j < ny; j++)
							{
								j1 = j - nyd2;
								eta2 = j1 * j1 * deta2 + zeta2;
								for (int i = 0; i < nx; i++)
								{
									x = Samp3[k][j][i];
									if (x == 0) { V3[k][j][i] = 0.0f; continue; }

									i1 = i - nxd2;
									sintemp = sin(dksi2 * i1 * i1 + eta2);
									sin2 = sintemp * sintemp;
									fInputNSR = float(dInputNSR * sqrt(x) * sintemp / (sin2 + epsilon));

									y = abs(V3[k][j][i]) / fInputNSR; // SNR in this Fourier coefficient
									V3[k][j][i] *= (1 / x) * (0.5f + atan(100.0f * (y - fTargetSNR)) / PIf); // suppress V3 elements with SNR < TargetSNR
									//if (y < fTargetSNR) V3[k][j][i] = 0.0f; // "hard" cut-off filter
								}
							}
						}
					}
					break;
				*/
				}
			} // end of noisy input images 3D filtration case

			// delete Samp3 matrix to free memory prior to 3D FFT
			Samp3.Truncate();

			// do 3D CTF correction if necessary
			if (zoutAver != 0)
			{
				if (iCTFcorrectionMode == 2)
				{
					printf("\nPerforming 3D TIE-Hom correction ...");
					//normalization of the inverse 2D minus-Laplacian here corresponds to the normalization of the CTFcorrection3D() function below
					double alphatemp = -4.0 * PI * asigma / (zoutAver * wl); // asigma should be negative, so alpha is positive
					double normtemp = -alphatemp * sqrt(1.0 + asigma * asigma) / asigma;
					InverseMLaplacian3D(V3, zhi - zlo, yhi - ylo, xhi - ylo, alphatemp, normtemp);
				}
				else if (iCTFcorrectionMode == 4)
				{
					printf("\nPerforming 3D CTF correction ...");
					CTFcorrection3D(V3, Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi), zoutAver, k2maxo, Cs3, Cs5, asigma, epsilon, bESCC);
					// the output(V3) is normalized here as delta * 4 * PI * sqrt(1 + sigma ^ 2) / wl
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
			if (!fftwf_init_threads())
				std::runtime_error("Error: multi-threaded initialization of FFTW library failed.");
			Fftwd3frc fft3f((int)nz, (int)ny, (int)nx, nThreads, false, true, false); // create an Fftwd3frc object in the Fourier space state
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

			// suppression of points outside the reconstruction sphere and compensation for the "reverse cupping" artefact
			index_t nxR2 = (ic <= nxd2) ? ic * ic : (nx - ic) * (nx - ic);
			index_t nyR2 = (jc <= nyd2) ? jc * jc : (ny - jc) * (ny - jc);
			index_t nzR2 = (kc <= nzd2) ? kc * kc : (nz - kc) * (nz - kc);
			index_t nR2;
			if (bRotYonly) nR2 = std::min(nxR2, nzR2);
			else if (bRotZonly || bRotZ2only) nR2 = std::min(nxR2, nyR2);
			else nR2 = std::min(min(nxR2, nyR2), nzR2);
			double dtemp = 2.0 / double(kc * kc + jc * jc + ic * ic);
			#pragma omp parallel for
			for (int k = 0; k < nz; k++)
			{
				index_t kk2(0), jj2, ii2;
				if (!bRotZonly && !bRotZ2only) kk2 = (k - kc) * (k - kc); // if the rotation is around Z or Z" axis only, the Z extent of the reconstructed volume should not be limited
				for (index_t j = 0; j < ny; j++)
				{
					if (!bRotYonly) jj2 = (j - jc) * (j - jc) + kk2; // if the rotation is around Y axis only, the Y extent of the reconstructed volume should not be limited
					else jj2 = kk2;
					for (index_t i = 0; i < nx; i++)
					{
						ii2 = (i - ic) * (i - ic) + jj2;
						if (ii2 > nR2)
							K3out[k][j][i] = 0.0f;
						else
							K3out[k][j][i] *= float(exp(dtemp * ii2));
					}
				}
			}

			// final normalization of the result
			// V3 is normalized at this point as delta * 4 * PI * sqrt(1 + sigma ^ 2) * zst / wl - see the explanation that follows:
			// V3 in the Fourier space was obtained from K2four by 3D gridding (the gridding was self-normalized later by dividing by Samp3).
			// K2four was normalized as delta * 4 * PI * sqrt(1 + sigma^2) / wl / (xst * yst), the (xst * yst) term was due to 2D DFT.
			// Indeed, 2D forward DFT produces an extra factor of 1 / (xst * yst), compared to the true 2D Fourier transform in the integral form.
			// V3 was then 3D inverse DFT transformed. Note that 3D inverse DFT produces an extra factor of 1 / (dksi * deta * dzeta), 
			// compared to the true inverse 3D Fourier transform in the integral form. Overall, we then have
			// 1 / (xst * dksi * yst * deta * dzeta) = (xaper / xst) * (yaper / yst) * zaper = nx * ny * zaper. 
			// Finally, GetRealXArray3D() contains the normalization factor 1 / (nx * ny * nz), and we end up with the extra factor zaper / nz = zst.
			float fnorm;
			if (iModality == 0)
				fnorm = (float)(EE * wl / (2.0 * PI * sqrt(1.0 + asigma * asigma) * zst)); // this is the conversion factor for V = 2E delta
			if (iModality == 1)
				fnorm = (float)(asigma * wl / (4.0 * PI * sqrt(1.0 + asigma * asigma) * zst)); // this is the conversion factor for beta = sigma * delta
				//fnorm = (float)(wl / (4.0 * PI * sqrt(1.0 + asigma * asigma))); // this is the conversion factor for delta
			K3out *= fnorm;

		} // end of case if imodeReprocess == 0
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
					if (nx0 != nx || ny0 != ny)
					{
						if (nx0 < nx && ny0 <= ny || ny0 < ny && nx0 <= nx)
						{
							//printf("\nWARNING: 2D array from the image file will be padded with 0.0 to the array dimensions (Width = %zd, Height = %zd).", nx, ny);
							XArray2DMove<float> tmp2(ipIn);
							tmp2.Pad((ny - ny0) / 2, ny - ny0 - (ny - ny0) / 2, (nx - nx0) / 2, nx - nx0 - (nx - nx0) / 2, 0.0);
						}
						else if (nx0 > nx && ny0 >= ny || ny0 > ny && nx0 >= nx)
						{
							//printf("\nWARNING: 2D array from the image file will be trimmed to the array dimensions (Width = %zd, Height = %zd).", nx, ny);
							XArray2DMove<float> tmp2(ipIn);
							tmp2.Trim((ny0 - ny) / 2, ny0 - ny - (ny0 - ny) / 2, (nx0 - nx) / 2, nx0 - nx - (nx0 - nx) / 2);
							// NOTE that the reconstructed 3D volume will NOT be padded back before saving, as it does not seem to make sense
						}
						else
							throw std::runtime_error("Error: the dimensions of the 2D array in the image file are inconsistent with the array dimensions in input parameter file (both dimensions should be equal, smaller or larger simultaneously).");
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
		if (imodeInvLaplace && imodeReprocess || (dlpfiltersize > 0 && dlpfiltersize < (xhi - xlo)))
		{
			if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of inverse 3D FFT time

			Fftwd3frc fft3f((int)nz, (int)ny, (int)nx, nThreads, true, true, false); // create an Fftwd3frc object in the real space state

			if (imodeInvLaplace && imodeReprocess) // in the case !imodeReprocess inverse Laplacian is calculated in the Fourier space
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

		// rotate the 3D array for output if required
		if (b3Drotout)
		{
			XArray3D<float> K3outRot(K3out);
			XArray3DSpln<float> xa3spln(K3outRot);
			printf("\nRotating around Z by %g (deg), around Y' by %g (deg) and around Z'' by %g (deg) ...", angleZout, angleY1out, angleZ2out);
			// this corresponds to the same resultant rotation as in pdb_run.cpp
			xa3spln.Rotate3(K3out, angleZout * PI180, angleY1out * PI180, angleZ2out * PI180, nThreads);
		}

		// trim back the reconstructed volume if input images were padded
		if (bPadded)
		{
			K3out.Trim((nx - nx0) / 2, nx - nx0 - (nx - nx0) / 2, (ny - ny0) / 2, ny - ny0 - (ny - ny0) / 2, (nx - nx0) / 2, nx - nx0 - (nx - nx0) / 2);
			nz = nx0;
			ny = ny0;
			nx = nx0;
			nxny = nx * ny;
			xlo = -0.5 * xst * nx;
			xhi = xlo + xst * nx;
			ylo = -0.5 * yst * ny;
			yhi = ylo + yst * ny;
			if (iModality) { zlo = xlo; zhi = xhi; } // X-ray case
			else // TEM case
			{ 
				zlo += double(index_t((nx - nx0) / 2)) * zst;  
				zlo -= double(index_t(nx - nx0 - (nx - nx0) / 2)) * zst; 
			}
			K3out.SetHeadPtr(new Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
		}

		// in the TEM case, place the reconstructed volume into hte positive octant in 3D Cartesian space
		if (iModality == 0)
		{
			xhi -= xlo; xlo = 0.0;
			yhi -= ylo; ylo = 0.0;
			zhi -= zlo; zlo = 0.0;
			K3out.SetHeadPtr(new Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
		}

		// find peaks in the reconstructed 3D distribution
		if (imodePeaks)
		{
			printf("\nSearching for peak positions in the 3D distribution of the electrostatic potential or beta ...");
			if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of peak localization time

			int natom = FindPeaks(K3out, datomsizeXY, datomsizeZ, natommax, filenamebaseOut);

			vector<string> voutfilenamesPeaksTot; // output filenames for the peak-localized reconstructed 3D potential or beta
			std::filesystem::path apath(filenamebaseOut);
			string filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "A");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			FileNames(1, nz, filenamebaseOutNew, voutfilenamesPeaksTot); // create 1D array of output filenames to save 2D slices of the peak-localized 3D object
			printf("\nSaving the reconstructed peak-localized 3D object into output files %s, etc. ...", voutfilenamesPeaksTot[0].c_str());

			if (bTIFFoutput) TIFFWriteFileStack(K3out, voutfilenamesPeaksTot, eTIFF32);
			else XArData::WriteFileStackGRD(K3out, voutfilenamesPeaksTot, eGRDBIN);

			if (bVerboseOutput) liPeakLocalizeTime = (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
		} // end of peak localization module

		// reslice and flip axes in the output 3D array
		if (iResliceMode || bXFlip || bYFlip || bZFlip)
		{
			printf("\nReslicing the output 3D volume ...");
			K3out = K3out.Reslice(iResliceMode, bXFlip, bYFlip, bZFlip);   
			nz = K3out.GetDim1(); ny = K3out.GetDim2(); nx = K3out.GetDim3();
		}

		// create thick slices if required
		if (iThickSlicesOut > 0)
		{
			printf("\nCreating thick output slices ...");
			if (kThickStep > nz) throw std::runtime_error("Slice thickness is larger than the sample thickness.");
			index_t nzThick = int(nz / kThickStep);

			XArray3D<float> K3outThick(nzThick, ny, nx);
			K3outThick.SetHeadPtr(K3out.GetHeadPtr()->Clone());

			#pragma omp parallel for
			for (int j = 0; j < ny; j++)
			{
				index_t k;
				float ftemp;
				for (index_t i = 0; i < nx; i++)
					for (index_t k1 = 0; k1 < nzThick; k1++)
					{
						k = k1 * kThickStep;
						ftemp = K3out[k][j][i];
						if (iThickSlicesOut == 1) // MIP
						{
							for (index_t n = 1; n < kThickStep; n++) if (K3out[k + n][j][i] > ftemp) ftemp = K3out[k + n][j][i];
						}
						else // averaging over the slice thickness
						{
							for (index_t n = 1; n < kThickStep; n++)  ftemp += K3out[k + n][j][i];
							ftemp /= float(kThickStep);
						}
						K3outThick[k1][j][i] = ftemp;
					}
			}
			K3out = K3outThick;
			nz = nzThick; // thin slices can disappear completely at this stage
		}

		// creat 16-bit output if required
		XArray3D<long> K3outInt; // we declare this "long", even though we use it for 16-bit output, because XArray<T> is not defined for T = unsigned short
		if (b16bitout) // 16-bit output
		{
			printf("\nCreating 16-bit output slices ...");
			K3outInt.Resize(nz, ny, nx);
			K3out.Convert(&K3outInt[0][0][0], fbit16_minin, fbit16_maxin, long(0), lbit16_maxout, nullptr, nullptr);
		}

		vector<string> voutfilenamesTot; // output filenames for the reconstructed 3D potential or beta
		if (imodeReprocess) // output the the renormalized 3D potential or beta
		{
			std::filesystem::path apath(filenamebaseOut);
			string filenamebaseOutNew = apath.filename().string();
			filenamebaseOutNew.insert(filenamebaseOutNew.find_last_of("."), "R");
			if (bDICOMoutput) filenamebaseOutNew.replace(filenamebaseOutNew.find_last_of("."), filenamebaseOutNew.length() - filenamebaseOutNew.find_last_of("."), ".dcm");
			filenamebaseOutNew = folderAux + filenamebaseOutNew;
			printf("\nOutput file name base for the renormalized 3D potential or beta = %s", filenamebaseOutNew.c_str());
			if (bDICOMoutput)
				FileNames(1, nz, filenamebaseOutNew, voutfilenamesTot, true); // create 1D array of output filenames to save 2D slices of the filtered and/or rescaled 3D object
			else
				FileNames(1, nz, filenamebaseOutNew, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the filtered and/or rescaled 3D object
		}
		else
		{
			if (bDICOMoutput)
				FileNames(1, nz, filenamebaseOut, voutfilenamesTot, true); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
			else
				FileNames(1, nz, filenamebaseOut, voutfilenamesTot); // create 1D array of output filenames to save 2D slices of the reconstructed 3D object
		}

		printf("\n\nSaving the reconstructed 3D object into output files %s, etc. ...", voutfilenamesTot[0].c_str());
		if (bVerboseOutput) start_timeNow = std::chrono::system_clock::now(); // start of file write time

		if (bTIFFoutput)
		{
			if (b16bitout) // 16-bit TIFF output
				TIFFWriteFileStack(K3outInt, voutfilenamesTot, eTIFF16, false); // 16-bit thin or thick slices
			else // 32-bit TIFF output
				TIFFWriteFileStack(K3out, voutfilenamesTot, eTIFF32); // 32-bit thin or thick slices
		}
		else if (bDICOMoutput)
		{
			if (b16bitout) // 16-bit DICOM output
				WriteFileStackDICOM(K3outInt, mTags, voutfilenamesTot, eDICOM16, false); // 16-bit thin or thick slices
			else
				// 32-bit DICOM output (this appears to be not implemented in GDCM ver.3.0.18 (Sep 2022))
				// we should not be able to get here (there are checks at the start to prevent this), but in case we do, we throw an exceptioin here
				throw std::runtime_error("Error: floating-point data cannot be saved in DICOM files.");
		}
		else // GRD output
			if (b16bitout)
				// we should not be able to get here (there are checks at the start to prevent this), but in case we do, we throw an exceptioin here
				throw std::runtime_error("Error: 16-bit data data cannot be saved in GRD files.");
			else
				XArData::WriteFileStackGRD(K3out, voutfilenamesTot, eGRDBIN); // 32-bit thin or thick slices

		if (bVerboseOutput) liFileWriteTime += (long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_timeNow).count();
		

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
		//fftw_cleanup_threads();
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

