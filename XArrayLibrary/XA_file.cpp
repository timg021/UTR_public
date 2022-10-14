//Module xa_file.cpp
//
//
//	MODULE TITLE:
//
//		Common file facilities for XArray and related classes
//
/*!
	\file		xa_file.cpp
	\brief		Common file facilities for XArray and related classes
	\par		Description:
		This module contains implementations of the functions declared in xa_file.h
*/

//---------------------------------------------------------------------------
//	INCLUDE FILES
//
#include "stdafx.h"

#include <algorithm>
#include <locale>

#include "XA_file.h"

//---------------------------------------------------------------------------
//	LOCAL CONSTANT DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL MACRO DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL ENUMERATED DATA TYPES
//
//---------------------------------------------------------------------------
//	LOCAL STRUCTURE DEFINITIONS
//
//---------------------------------------------------------------------------
//	LOCAL IN-LINE FUNCTION DEFINITIONS
//
//---------------------------------------------------------------------------
//	GLOBAL FUNCTION PROTOTYPES
//
//---------------------------------------------------------------------------
//	LOCAL FUNCTION PROTOTYPES
//
//---------------------------------------------------------------------------
//	GLOBAL DATA DEFINITIONS
//
//---------------------------------------------------------------------------
//	STATIC DATA DEFINITIONS
//
/////////////////////////////////////////////////////////////////////////////

xar::string str_toupper(xar::string s) 
{
	std::transform(s.begin(), s.end(), s.begin(),
		[](unsigned char c) { return std::toupper(c); }
	);
	return s;
}


xar::string xar::GetFileExtension(const string& filename, bool Convert2Upper)
// returns ".EXT"
{
	string strTemp = GetFilenameFromPath(filename);
	index_t i = filename.rfind('.');
	if (i != string::npos) strTemp = filename.substr(i);
	else return string("");
	if (Convert2Upper) strTemp = str_toupper(strTemp);
	return strTemp;
}


xar::string xar::GetFilenameFromPath(const string& strPath, bool Convert2Upper)
{
	index_t i = strPath.rfind('\\');
	if (i == string::npos)
	{
		i = strPath.rfind('/');
		if (i == string::npos) return strPath; // no path present
	}
	if (i == string::npos) return strPath; // no path present
	if (i == strPath.size() - 1) return string(""); // no filename present
	string strTemp = strPath.substr(i + 1);
	if (Convert2Upper) strTemp = str_toupper(strTemp);
	return strTemp;
}


xar::string xar::GetPathFromFilename(const string& strFilename, bool Convert2Upper)
// returns "path\"
{
	index_t i = strFilename.rfind('\\');
	if (i == string::npos)
	{
		i = strFilename.rfind('/');
		if (i == string::npos) return string(""); // no path present
	}
	if (i == strFilename.size() - 1) return strFilename; // no filename present
	string strTemp = strFilename.substr(0, i + 1);
	if (Convert2Upper) strTemp = str_toupper(strTemp);
	return strTemp;
}


bool xar::DoesFileExist(const char* filename)
{
	return std::filesystem::exists(filename);
	//int fh = _open(filename, _O_RDONLY);
	//if (fh==-1) return false;
	//else { _close(fh); return true; }
}


xar::string xar::SortFiles(const string& infiles)
// sorts long strings of comma-separated filenames (the first one containing the full path)
// infiles - input filenames
// returns alphabetically sorted infiles
{
	index_t i, i0, i1;
	string infilesnames(GetFilenameFromPath(infiles, false));	
	list<string> listInputFiles;

	if (infilesnames.size() <= 0) return infiles;

	i = i0 = i1 = 0;
	while (i != string::npos)
	{
		i = infilesnames.find(',', i0);
		if (i == string::npos) i1 = infilesnames.size();
		else i1 = i;
		listInputFiles.push_back(infilesnames.substr(i0, i1 - i0));
		i0 = i1 + 1;
	}

	listInputFiles.sort();

	string result = GetPathFromFilename(infiles, false);
	list<string>::iterator iter = listInputFiles.begin();
	while (iter != listInputFiles.end())
	{
		result += *iter++;
		result += ',';
	}
	result.erase(result.size() - 1);

	return result;
}


index_t xar::NumberOfFiles(const string& infiles)
// analyses long strings of comma-separated filenames (the first one containing the full path)
// returns the number of input files 
{
	index_t i = 0, nfiles = 0;

	if (infiles.size() <= 0) return 0;
	while (i != string::npos)
	{
		nfiles++;
		i = infiles.find_first_of(',', i+1);
	}

	return nfiles;
}


index_t xar::IOFileLists(const string& infiles, const string& outfile, list<string>& listInputFiles, list<string>& listOutputFiles, int* piNoProblem)
// analyses long strings of comma-separated filenames (the first one containing the full path)
// and fills the lists listInputFiles and listOutputFiles with input and output filenames
// infiles - input filenames
// outfile - output filename used to extract the output path, prefix and suffix to be added
//				to input filenames when used as output filenames
// piNoProblem - pointer to an integer variable to return the "error code":
//				0 = no problems
//				1 = duplicate output filenames generated
//				2 = some output files already exist
//				3 = duplicate output filenames generated AND some output files already exist
// returns the number of input files (== number of output files)
{
	index_t i, i0, i1, j, nfiles;
	int iProblem1(0), iProblem2(0);
	string infilesnames(GetFilenameFromPath(infiles, false));	
	string csTemp, csOutPath, csOutSuffix, csOutPrefix;

	listInputFiles.clear();
	listOutputFiles.clear();

	if (infilesnames.size() <= 0) return 0;

	csOutPath = GetPathFromFilename(outfile, false);
	if (csOutPath.empty())
		csOutPath = GetPathFromFilename(infiles, false);
	csTemp = GetFilenameFromPath(outfile, false);
	j = csTemp.rfind('.');
	if (j == string::npos) j = csTemp.size();
	csOutPrefix = csTemp.substr(0, j); // "NAME"
	csOutSuffix = csTemp.substr(j); // ".EXT"
	
	i = i0 = i1 = 0;
	nfiles = 0;
	while (i != string::npos)
	{
		nfiles++;
		i = infilesnames.find(',', i0);
		if (i == string::npos) i1 = infilesnames.size();
		else i1 = i;
		csTemp = infilesnames.substr(i0, i1 - i0);
		listInputFiles.push_back(csTemp);
		csTemp.insert(0, csOutPrefix); // insert NAME on the left
		if (csOutSuffix.size() > 0) // insert or replace ".EXT" on the right
		{
			j = csTemp.rfind('.');
			if (j == string::npos) // no extension in the input name
				csTemp.insert(csTemp.size(), csOutSuffix);
			else 
				csTemp.replace(j, csTemp.size() - j, csOutSuffix);
		}
		if (piNoProblem) 
			if (std::find(listOutputFiles.begin(), listOutputFiles.end(), csTemp) != listOutputFiles.end()) 
				iProblem1 = 1;
		listOutputFiles.push_back(csTemp);
		csTemp = csOutPath + csTemp;
		if (piNoProblem && DoesFileExist(csTemp.c_str())) iProblem2 = 2;
		i0 = i1 + 1;
	}

	if (piNoProblem) *piNoProblem = iProblem1 + iProblem2;
	return nfiles;
}


index_t xar::IOFileSequence(const string& infiles, const string& outbase, index_t nOutFiles, list<string>& listInputFiles, list<string>& listOutputFiles, bool bIndexFormatInOut, int* piNoProblem)
// analyses long strings of comma-separated filenames (the first one containing the full path)
// creates a sequence of output filenames containing consecutive numbers
// fills the lists listInputFiles and listOutputFiles with input and output filenames
// infiles - input filenames
// outbase - output base filename used to extract the output path, prefix and suffix for the output filenames
// nOutFiles - number of output files
// bIndexFormatInOut - if true, the number of digits in the output file indexes is determined from the number of digits in the index of the last input file
//			if false, the number of digits in the output file indexes is determined from the total number of output files
// piNoProblem - pointer to an integer variable to return the "error code":
//				0 = no problems
//				2 = some output files already exist
// returns the number of input files (!= number of output files !!!)
{
	vector<index_t> vIndexes(nOutFiles);
	for (index_t i = 0; i < nOutFiles; i++) vIndexes[i] = i;
	return IOFileSet(infiles, outbase, vIndexes, listInputFiles, listOutputFiles, bIndexFormatInOut, piNoProblem);
}


index_t xar::IOFileSet(const string& infiles, const string& outbase, vector<index_t> vIndexes, list<string>& listInputFiles, list<string>& listOutputFiles, bool bIndexFormatInOut, int* piNoProblem)
// analyses long strings of comma-separated filenames (the first one containing the full path)
// creates a sequence of output filenames with given indexes
// fills the lists listInputFiles and listOutputFiles with input and output filenames
// infiles - input filenames
// outbase - output base filename used to extract the output path, prefix and suffix for the output filenames
// vIndexes - vector of indexes for the output files
// bIndexFormatInOut - if true, the number of digits in the output file indexes is determined from the number of digits in the index of the last input file
//			if false, the number of digits in the output file indexes is determined from the total number of output files
// piNoProblem - pointer to an integer variable to return the "error code":
//				0 = no problems
//				2 = some output files already exist
// returns the number of input files (!= number of output files !!!)
{
	index_t i, i0, i1, j, nfiles;
	int iProblem2(0);
	string infilesnames(GetFilenameFromPath(infiles, false));	
	string csTemp, csOutPath, csOutSuffix, csOutPrefix;

	listInputFiles.clear();
	listOutputFiles.clear();

    if (infilesnames.size() <= 0) return 0;

	csOutPath = GetPathFromFilename(outbase, false);
	if (csOutPath.empty())
		csOutPath = GetPathFromFilename(infiles, false);
	csTemp = GetFilenameFromPath(outbase, false);
	j = csTemp.rfind('.');
	csOutPrefix = csTemp.substr(0, j); // "NAME"
	if (j == string::npos) 
		csOutSuffix = ".VOD"; // default extension = ".VOD"
	else
		csOutSuffix = csTemp.substr(j); // ".EXT"
	
	// create list of input files
	i = i0 = i1 = 0;
	nfiles = 0;
	while (i != string::npos)
	{
		nfiles++;
		i = infilesnames.find(',', i0);
		if (i == string::npos) i1 = infilesnames.size();
		else i1 = i;
		csTemp = infilesnames.substr(i0, i1 - i0);
		listInputFiles.push_back(csTemp);
		i0 = i1 + 1;
	}

	// find the required number of digits in the output slice filenames
	int ndigits(0);
	if (bIndexFormatInOut)
	{
		csTemp = listInputFiles.back();
		i1 = csTemp.find_last_of('.') - 1;
		i0 = csTemp.find_last_not_of("0123456789", i1);
		ndigits = int(i1) - int(i0);
	}

	//create list of output files
	char outformat[64];
	// sort vIndexes in the ascending order
	std::sort(vIndexes.begin(), vIndexes.end()); 
	 // remove duplicate elements from vIndexes
	index_t nOutputFiles = std::unique(vIndexes.begin(), vIndexes.end()) - vIndexes.begin();
	// determine the number of digits in the output file indexes from the total number of output files
	if (!ndigits) ndigits=(int)((vIndexes[nOutputFiles-1] > 0) ? log10((float)vIndexes[nOutputFiles-1]) + 1 : 1);
	sprintf(outformat,"%%s%%.%dd%s",ndigits,csOutSuffix.c_str());
	vector<char> vTemp(csOutPrefix.size() + ndigits + csOutSuffix.size() + 1, '\0');

	for (index_t n = 0; n < nOutputFiles; n++)
	{
		sprintf(&vTemp[0], outformat, csOutPrefix.c_str(), vIndexes[n]);
		csTemp = &vTemp[0];
		listOutputFiles.push_back(csTemp);

		csTemp = csOutPath + csTemp;

		if (piNoProblem && DoesFileExist(csTemp.c_str())) iProblem2 = 2;
	}

	if (piNoProblem) *piNoProblem = iProblem2;

	return nfiles;
}


index_t xar::IOFileString2List(const string& infiles, list<string>& listInputFiles, int* piNoProblem)
// analyses long strings of comma-separated filenames (the first one containing the full path)
// and fills the list listInputFiles with input filenames
// infiles - input filenames
// piNoProblem - pointer to an integer variable to return the "error code":
//				0 = no problems
//				1 = duplicate output filenames generated
//				2 = some output files already exist
//				3 = duplicate output filenames generated AND some output files already exist
// returns the number of input files (== number of output files)
{
	index_t i, i0, i1, nfiles;
	int iProblem1(0), iProblem2(0);
	string infilesnames(GetFilenameFromPath(infiles, false));	
	string csInPath(GetPathFromFilename(infiles, false));
	string csTemp;

	if (infilesnames.size() <= 0) return 0;

	listInputFiles.clear();

	i = i0 = i1 = 0;
	nfiles = 0;
	while (i != string::npos)
	{
		nfiles++;
		i = infilesnames.find(',', i0);
		if (i == string::npos) i1 = infilesnames.size();
		else i1 = i;
		csTemp = infilesnames.substr(i0, i1 - i0);
		if (piNoProblem && !iProblem1) 
			if (std::find(listInputFiles.begin(), listInputFiles.end(), csTemp) != listInputFiles.end()) 
				iProblem1 = 1;
		listInputFiles.push_back(csTemp);
		csTemp = csInPath + csTemp;
		if (piNoProblem && !iProblem2 && DoesFileExist(csTemp.c_str())) iProblem2 = 2;
		i0 = i1 + 1;
	}

	if (piNoProblem) *piNoProblem = iProblem1 + iProblem2;
	return nfiles;
}


index_t xar::FindRefractiveIndices(const string strMaterialLegendFilename, const vector<float> vEnergies, vector<vector<float> >& vOutDelta, vector<vector<float> >& vOutBeta)
// strMaterialLegendFilename - full name of a 2-column file containing the association between material indices and the corresponding files with material data;
//						NOTE!!!: integer indices of material components must start from 0, be consecutive and in ascending order (i.e. 0, 1, 2, ...)
//						The material files are assumed to be in Adisp.exe (A.Stevenson) 4-column format containing
//						(1) energies in keV with 0.1 keV steps, (2) Mu in cm^-1, (3) T for 10% trans. in cm, (4) Phi in rad/mm;
//						NOTE that the two first rows are occupied by column titles.
// vEnergies - vector of energies present in the spectrum (sorted in ascending order)
// vOutDelta - output array of delta values indexed by vMaterialIndices in the 1st index, and by energies in the 2nd index
// vOutBeta - output array of beta values indexed by vMaterialIndices in the 1st index, and by energies in the 2nd index
// returns the number of material indices in strMaterialLegend file
{
	// read Material Legend file which relates material indices to material file names
	char buf[2048], buf1[129], buf2[1024], buf3[129], buf4[129];
	FilePtr fp(strMaterialLegendFilename.c_str(), "rt");

	index_t nMaterials = 0;
	for (;;) 	// first, count the number of rows in the file
	{ 
		if (fgets(buf, 2048, fp) == 0)
			break;
		else
		{
			nMaterials++;
		}
	}
	if (nMaterials == 0) throw std::runtime_error("runtime error in XArCbTomo<T>::FindRefractiveIndices() (bad material legend file)");
	rewind (fp);
	vector<string> vMaterialFiles(nMaterials);
	for (index_t i = 0; i < nMaterials; i++) 	// now, read the file data
	{
		fgets(buf, 2048, fp);
		sscanf(buf, "%s %s", buf1, buf2);
		if (index_t(strtod(buf1, 0)) != i) throw std::runtime_error("runtime error in XArCbTomo<T>::FindRefractiveIndices() (bad material legend file: material indices must be consecutive)");
		vMaterialFiles[i] = string(buf2);
	}

	// resize arrays
	vOutDelta.resize(nMaterials);
	vOutBeta.resize(nMaterials);
	index_t nEnergies = vEnergies.size();
	for (index_t i = 0; i < nMaterials; i++) 
	{ 
		vOutDelta[i].resize(nEnergies, 0.0);
		vOutBeta[i].resize(nEnergies, 0.0);
	}

	// read material file 
	double lambda, dlambda;
	vector<double> Energy, beta, delta;
	for (index_t i = 0; i < nMaterials; i++) 
	{ 
		FilePtr fp1(vMaterialFiles[i].c_str(), "rt"); // can throw an exception
		// first, count the number of rows in the file
		index_t iSize = 0;
		for (;;) 
		{ 
			if (fgets(buf, 2048, fp1) == 0)
				break;
			else
			{
				iSize++;
			}
		}
		if (iSize == 0) throw std::runtime_error("runtime error in XArCbTomo<T>::FindRefractiveIndices() (bad material file)");
		// now, read the file data
		rewind (fp1);
		Energy.resize(iSize - 2); beta.resize(iSize - 2); delta.resize(iSize - 2);
		for (index_t j = 0; j < 2; j++) // skip column titles
			fgets(buf, 2048, fp1);
		for (index_t j = 0; j < iSize - 2; j++)
		{
			fgets(buf, 2048, fp1);
			sscanf(buf, "%s %s %s %s", buf1, buf2, buf3, buf4);
			Energy[j] = strtod(buf1, 0); // energy in keV
			beta[j] = strtod(buf2, 0); // actually this is mu in cm^-1
			delta[j] = strtod(buf4, 0); // actually this is phi in rad/mm
			lambda = 12.396e-4 / Energy[j]; // wavelength in microns
			beta[j] *= lambda / (4.0 * PI) * 1.e-4;
			delta[j] *= -lambda / tPI * 1.e-3; 
		}

		// fill in the output arrays
		index_t k = 0, j;
		for (j = 0; j < iSize - 2; j++)
		{
			if(Energy[j] >= vEnergies[k])
			{
				dlambda = Energy[j] / vEnergies[k];
				vOutBeta[i][k] = float(pow(dlambda, 4) * beta[j]);
				vOutDelta[i][k] = float(pow(dlambda, 2) * delta[j]);
				k++;
			}
			if (k >= nEnergies) break;
		}
		if (j == iSize - 2) throw std::runtime_error("runtime error in XArCbTomo<T>::FindRefractiveIndices() (X-ray energy outside the range present in material file)");
	}

	return nMaterials;
}


index_t xar::ReadSpectrumFile(const string strSpectrumFilename, vector<float>& vOutEnergies, vector<float>& vOutCounts)
// strSpectrumFilename - full name of a 2-column spectrum file containing photon energies in keV in the first column and the corresponding number of "particles" in the second column;
//						NOTE: the conversion efficiency of the detector may be implicitly included by multiplying the numbers in the second column
// vOutEnergies - output vector of energies present in the spectrum (sorted in ascending order)
// vOutCounts - output vector containing either photon fluencies or detetor counts depending on the desired interpretation
// returns the number of different energies in the spectrum file
{
	char buf[2048], buf1[1024], buf2[1024];
	FilePtr fp(strSpectrumFilename.c_str(), "rt");

	// first, count the number of rows in the file
	index_t nEnergies = 0;
	for (;;)
	{ 
		if (fgets(buf, 2048, fp) == 0)
			break;
		else
		{
			nEnergies++;
		}
	}
	if (nEnergies == 0) throw std::runtime_error("runtime error in XArCbTomo<T>::ReadSpectrumFile() (bad spectrum file)");

	// now, read the file data
	rewind (fp);
	vector<float> vEnergies(nEnergies), vCounts(nEnergies); // unsorted vectors
	for (index_t i = 0; i < nEnergies; i++)
	{
		fgets(buf, 2048, fp);
		sscanf(buf, "%s %s", buf1, buf2);
		vEnergies[i] = (float)strtod(buf1, 0); // energy in keV;
		vCounts[i] = (float)strtod(buf2, 0); // counts
		if (vCounts[i] <= 0 || vEnergies[i] <= 0)
			throw std::runtime_error("bad data in 'strSpectrumFile' in xar::ReadSpectrumFile() (all energies and count numbers must be positive)"); 
	}

	// sort
	vOutEnergies.resize(nEnergies);
	vOutCounts.resize(nEnergies);
	vOutEnergies = vEnergies;
	std::sort(vOutEnergies.begin(), vOutEnergies.end()); 
	for (index_t i = 0; i < nEnergies; i++)
		for (index_t j = 0; j < nEnergies; j++)
			if (vEnergies[j] == vOutEnergies[i]) 
			{ 
				vOutCounts[i] = vCounts[j]; 
				break; 
			}

	return nEnergies;
}


void xar::ReadDefocusParamsFile(string difile, vector<Pair>& v2angles, vector<vector<Pair> >& vvdefocus, bool bVerboseOutput)
// Reads defocus parameter data from a text file with each line containing two illumination direction angles (around Z and Y' axes) in degrees, followed by 
//			alternating values of molecule rotation angle around Z" axis and defocus distances at this angle in Angstroms,
//			with all values separated by a single white space and with the new line symbol at the end of each line
// difile - input file name
// v2angles - output vector of pairs (z, y') each corresponding to illumination direction defined by the rotation angles around Z and Y' axes
// vvdefocus - output vector of vectors of pairs (z", d) each containing a rotation angle around Z" axis (illumination axis) and a defocus distance at the angle triplet (z, y', z")
// any number of comment lines (starting with //) can be present inside the file, and will be ignored by this function
// !!!NOTE that this program changes the signs of all angles read from the input file in order to account for the fact that the input file is expected to contain
// !!!NOTE that the defocus distances are assumed to be measured from the centre of a molecule or nanoparticle, which should generally coincide with the centre of rotation
// the rotation angles for the coordinate system (around the immobile sample/structure), while we will be rotating the sample/structure instead of the coordinate system
{
	char cline[2049];
	char buffer[1024];
	double dtemp;
	vector<Pair> vdefocus(0); // vector of defocus distances at a given rotation angle
	vector<size_t> vwhite(0); // vector of white spaces separating different defocus distances (there should be exactly one white space before each defocus distance and no spaces at the end)

	FILE* ff0 = fopen(difile.c_str(), "rt");
	if (!ff0) throw std::runtime_error((string("Error: cannot open input file ") + difile + string(".")).c_str());
	//fgets(cline, 1024, ff0); // 1st line - comment

	Pair pair;
	v2angles.resize(0);
	vvdefocus.resize(0);
	index_t nline(0);
	std::locale locl;
	while (fgets(cline, 2048, ff0) != NULL) // read lines, each line consisting of three rotation angles followed by one or more defocus distance
	{
		nline++;
		if (cline[0] == '/' && cline[1] == '/') continue; // skip any comment lines
		if (!isdigit(cline[0]) && cline[0] != '-' && cline[0] != '+' && cline[0] != '.') // non-commented lines should contain only numbers
		{
			throw std::runtime_error((string("Bad first character in a line of input file ") + difile + string(".")).c_str());
		}
		strtok(cline, "\n");
		strtok(cline, "\r");
		vwhite.resize(1); vwhite[0] = 0;
		for (size_t i = 1; i < strlen(cline); i++) if (std::isspace(cline[i], locl)) vwhite.push_back(i); // count the number of different pairs of z" angles and defocus distances in the input parameter file
		if (vwhite.size() % 2)
		{
			printf("\nline no. = %zd, number of white spaces + 1 = %zd", nline, vwhite.size());
			throw std::runtime_error((string("Error reading input file (odd number of entries in a row) ") + difile + string(".")).c_str());
		}
		vwhite.push_back(strlen(cline)); // add one more entry corresponding to the end of the parameter string
		int ndefocus = int(vwhite.size() - 3); 
		if (ndefocus < 1) break; // assume that this is the end of file
		if ((ndefocus % 2) != 0 ) // this number should be >0 and even, as it is supposed to contain pairs of values (z", d)
		{
			fclose(ff0); // close input file
			throw std::runtime_error((string("Incorrect number of entries in a line of input file ") + difile + string(".")).c_str());
		}
		else ndefocus /= 2; // number of detected defocus distances
		vdefocus.resize(ndefocus); // vector of pairs (z", d) of z" angles and defocus distances (double precision numbers)
		if (bVerboseOutput) printf("\nRotations angles: ");
		for (size_t j = 0; j < 2; j++)
		{
			for (size_t i = vwhite[j]; i < vwhite[j + 1]; i++)
				buffer[i - vwhite[j]] = cline[i];
			buffer[vwhite[j + 1] - vwhite[j]] = '\0'; // string terminator
			dtemp = atof(buffer);
			j ? pair.b = -dtemp : pair.a = -dtemp; // !!! we change angle signs here, as we will use these to rotate the sample rather than the coordinate system
		}
		if (bVerboseOutput) printf("Z angle = %g, Y' angle = %g (degrees)", -pair.a, -pair.b);
		v2angles.push_back(pair);
		if (bVerboseOutput) printf("\nDefocus distances from the centre of rotation (%d in total): ", ndefocus);
		for (size_t j = 0; j < ndefocus * 2; j++)
		{
			for (size_t i = vwhite[j + 2]; i < vwhite[j + 3]; i++)
				buffer[i - vwhite[j + 2]] = cline[i];
			buffer[vwhite[j + 3] - vwhite[j + 2]] = '\0'; // string terminator
			dtemp = atof(buffer);
			(j % 2) ? (vdefocus[j / 2]).b = dtemp : (vdefocus[j / 2]).a = -dtemp; // odd entries are MINUS Z" angles, even entries are defocus distances
			if (j %2 && bVerboseOutput) printf("Z'' angle = %g (degrees), defocus = %g (A); ", -(vdefocus[j / 2]).a, (vdefocus[j / 2]).b);
		}
		vvdefocus.push_back(vdefocus);
	}
	if (bVerboseOutput) printf("\n%zd rotational positions in total in the input file.", vvdefocus.size());

	if (vvdefocus.size() == 0) throw std::runtime_error((string("Error: defocus parameter file ") + difile + string(" has zero valid entries.")).c_str());

	fclose(ff0); // close input file
}


void xar::ReadRelionDefocusParamsFile(string difile, vector<Pair>& v2angles, vector<vector<Pair> >& vvdefocus, vector<Pair>& vastigm, vector<Pair>& v2shifts, bool bVerboseOutput)
// Reads RELION defocus parameter data from a text file with each line containing 
// |particle# | filepath | shiftx | shifty | df1 | df2  | astigm | rot | tilt | psi |
// separated by a single white space and with the new line symbol at the end of each line
// difile - input file name
// v2angles - output vector of pairs (z, y') each corresponding to illumination direction defined by the rotation angles around Z and Y' axes in degrees
// v2shifts - X and Y shift of the image in Angstroms
// vvdefocus - output vector of vectors of pairs (z", d) each containing a rotation angle around Z" axis (illumination axis) and a defocus distance at the angle triplet (z, y', z")
// vastigm - astigmatism angle of the image in degrees - see J. Zivanov et al, IUCrJ 7, 2020 (doi.org/10.1107/S2052252520000081)
//			we assume that the input parameters are: df1 = Z1, df2 = Z2, astigm = phi_A, in the notation of Appendix A of the above paper,
//			the distances Z1 and Z2 are assumed to be in Angstroms, while the phi_A angle is assumed to be expressed in degrees,
//			and we transform them on output into: vvdefocus[n][0].b = (Z1 + Z2) / 2, vastigm[n].a = (Z1 - Z2) / 2, vastigm[n].b = phi_A
// any number of comment lines (starting with //) can be present inside the file, and will be ignored by this function
// !!!NOTE that this program changes the signs of all angles and XY shifts read from the input file in order to account for the fact that the input file is expected to contain
// !!!NOTE that the defocus distances are assumed to be measured from the centre of a molecule or nanoparticle, which should generally coincide with the centre of rotation
// the rotation angles for the coordinate system (around the immobile sample/structure), while we will be rotating the sample/structure instead of the coordinate system
{
	char cline[2049];
	char buffer[1024];
	double df1, df2, phiA, shiftx, shifty, rot, tilt, psi;
	vector<Pair> vdefocus(1); // vector of defocus distances at a given rotation angle

	FILE* ff0 = fopen(difile.c_str(), "rt");
	if (!ff0) throw std::runtime_error((string("Error: cannot open input file ") + difile + string(".")).c_str());
	//fgets(cline, 1024, ff0); // 1st line - comment

	Pair pair0, pair, pair1;
	v2angles.resize(0);
	vvdefocus.resize(0);
	v2shifts.resize(0);
	vastigm.resize(0);
	int nlineTemp;
	while (fgets(cline, 2048, ff0) != NULL) // read lines, each line consisting of three rotation angles followed by one or more defocus distance
	{
		if (cline[0] == '/' && cline[1] == '/') continue; // skip any comment lines
		if (!isdigit(cline[0]) && cline[0] != '-' && cline[0] != '+' && cline[0] != '.' && cline[0] != '\n') // non-commented lines should contain only numbers
			throw std::runtime_error((string("Error reading input file ") + difile + string(".")).c_str());
		strtok(cline, "\n");
		if (sscanf(cline, "%d %s %lg %lg %lg %lg %lg %lg %lg %lg", &nlineTemp, buffer, &shiftx, &shifty, &df1, &df2, &phiA, &rot, &tilt, &psi) != 10)
			break; // assume that this is the end of file

		pair.a = -rot; pair.b = -tilt;
		v2angles.push_back(pair);
		vdefocus[0].a = -psi; vdefocus[0].b = -0.5 * (df1 + df2);
		vvdefocus.push_back(vdefocus);
		pair0.a = -0.5 * (df1 - df2);
		pair0.b = phiA;
		vastigm.push_back(pair0);
		pair1.a = -shiftx; pair1.b = -shifty;
		v2shifts.push_back(pair1);

		if (bVerboseOutput) printf("\nZ_angle = %g Y'_angle = %g ", -pair.a, -pair.b);
		if (bVerboseOutput) printf("Z''_angle = %g (degrees), (DefocX+DefocY)/2 = %g (A), ", -vdefocus[0].a, -vdefocus[0].b);
		if (bVerboseOutput) printf("(DefocX-DefocY)/2 = %g (A), phiA = %g (deg), X-shift = %g (A), Y-shift = %g (A)", -pair0.a, pair0.b, -pair1.a, -pair1.b);
	}
	
	if (vvdefocus.size() == 0)
		throw std::runtime_error((string("Error reading input file ") + difile + string(" (number of entries in a line is not equal to 10?).")).c_str());
	else
		printf("\n%zd rotational positions read in total from the input file.", vvdefocus.size());

	fclose(ff0); // close input file
}


void xar::FileNames(index_t nangles, index_t ndefocus, string filenamebase, vector<string>& output, bool OneBased)
// Creates a sequence of file names properly indexed by rotation angles and defocus distances (using the same algorithm as in MultisliceK.cpp)
{
	if (ndefocus < 1 || nangles < 1)
		throw std::runtime_error("bad number of angles and/or defocus distances in FileNames()");

	char buffer[128];
	string strAngle, outfilename_i, outfilename_j;

	output.resize(ndefocus * nangles); // vector of full output filenames

	// create formatting string to add properly formatted indexes at the end of the output file names
	index_t i_dot = filenamebase.rfind('.'), nfieldA_length, nfieldB_length;
	char ndig[8];
	string myformat("");
	if (ndefocus > 1)
	{
		nfieldA_length = 1 + index_t(log10(double(ndefocus - 1))); //maximum number of digits corresponding to defocuses in the output file name
		sprintf(ndig, "%zd", nfieldA_length); //convert the calculated maximum number of digits corresponding to defocuses into a string, e.g. 5 into "5"
		myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded defocus indexes into file names
	}
	if (nangles > 1)
	{
		nfieldB_length = 1 + index_t(log10(double(nangles - 1))); //maximum number of digits corresponding to angles in the output file name
		sprintf(ndig, "%zd", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
		myformat += "_%0" + string(ndig) + "d"; //construct format string for inserting two 0-padded angle indexes into file names - see usage below
	}

	for (index_t i = 0; i < nangles; i++)
	{
		for (index_t j = 0; j < ndefocus; j++)
		{
			outfilename_j = filenamebase;
			if (OneBased) // numbering starts from 1
			{ 
				if (ndefocus == 1 && nangles > 1) sprintf(buffer, myformat.data(), i + 1);
				else if (ndefocus > 1 && nangles == 1) sprintf(buffer, myformat.data(), j + 1);
				else sprintf(buffer, myformat.data(), j + 1, i + 1);
			}
			else // numbering starts from 0
			{
				if (ndefocus == 1 && nangles > 1) sprintf(buffer, myformat.data(), i);
				else if (ndefocus > 1 && nangles == 1) sprintf(buffer, myformat.data(), j);
				else sprintf(buffer, myformat.data(), j, i);
			}
			outfilename_j.insert(i_dot, buffer);
			output[i * ndefocus + j] = outfilename_j;
		}
	}
}



void xar::FileNames2(vector<index_t> vndefocus, string filenamebase, vector<string>& output)
// Creates a sequence of file names properly indexed by rotation angles and defocus distances
// Similar to FileNames, but the number of defocus distances is allowed to be different at different rotation angles
// This function can also handle different zero padding and different starting indexes
{
	char buffer[128];
	string strAngle, outfilename_i, outfilename_j;
	index_t nangles = vndefocus.size(), nstart = 0, nstart1 = 0;
	if (nangles < 1) throw std::runtime_error("bad number of rotation angles in FileNames2()");

	bool bZeroPadded(true); // filename indexes are considered zero padded by default
	string stemp = filenamebase, stemp1, filenamebase0 = filenamebase;
	int idot = (int)stemp.rfind('.'); // position of the dot
	int ipos1, ipos = idot; // to become the position of the first digit in the index before the dot
	while (isdigit(stemp[--ipos])); ipos++; // position of the first digit
	int nfieldA_length = 0, nfieldB_length = int(idot - ipos); // length of the index for angles
	if (nfieldB_length > 0)
	{
		if (nfieldB_length == 1) bZeroPadded = false; // one digit in the filename index => no zero padding
		else if (stemp[ipos] != '0') bZeroPadded = false; // first digit in a multi-digit index is not zero => no zero padding
		stemp1 = stemp.substr(ipos, nfieldB_length);
		nstart = atoi(stemp1.c_str()); // starting index
		if (stemp[ipos - 1] == '_') // possible presence of the defocus index
		{
			ipos1 = ipos - 1; // to become the position of the first digit in the defocus index
			while (isdigit(stemp[--ipos1])); ipos1++; // position of the first digit
			nfieldA_length = int(ipos - 1 - ipos1); // length of the index for defocuses
			if (nfieldA_length > 0)
			{
				stemp1 = stemp.substr(ipos1, nfieldA_length);
				nstart1 = atoi(stemp1.c_str()); // starting index
			}
		}
	}
	if (nfieldA_length > 0)
		filenamebase0 = filenamebase.erase(ipos1, nfieldA_length + 1 + nfieldB_length); // filenamebase with both indexes removed
	else if (nfieldB_length > 0)
		filenamebase0 = filenamebase.erase(ipos, nfieldB_length); // filenamebase with removed index

	// calculate the total number of defocus distances, i.e. the total number of output files
	index_t ndefocustot(0), ndefocusmax(0);
	for (index_t i = 0; i < nangles; i++)
	{
		if (vndefocus[i] < 1) throw std::runtime_error("bad number of defocus distances in FileNames2()");
		ndefocustot += vndefocus[i];
		if (vndefocus[i] > ndefocusmax) ndefocusmax = vndefocus[i];
	}
	output.resize(ndefocustot); // vector of full output filenames

	// create formatting string to add properly formatted indexes at the end of the output file names
	char ndig[8];
	string myformat("");
	if (bZeroPadded)
	{
		if (ndefocusmax > 1)
		{
			int nfieldA_length1 = 1 + int(log10(double(ndefocusmax + nstart - 1))); //maximum number of digits corresponding to defocuses in the output file name
			if (nfieldA_length1 > nfieldA_length) nfieldA_length = nfieldA_length1; // overwrite the length of zero padding if it was absent or too short in filenamebase
			sprintf(ndig, "%d", nfieldA_length); //convert the calculated maximum number of digits corresponding to defocuses into a string, e.g. 5 into "5"
			myformat = "%0" + string(ndig) + "d"; //construct format string for inserting 0-padded defocus indexes into file names
		}
		if (nangles > 1)
		{
			int nfieldB_length1 = 1 + int(log10(double(nangles + nstart - 1))); //maximum number of digits corresponding to angles in the output file name
			if (nfieldB_length1 > nfieldB_length) nfieldB_length = nfieldB_length1; // overwrite the length of zero padding if it was absent or too short in filenamebase
			sprintf(ndig, "%d", nfieldB_length); //convert the calculated maximum number of digits corresponding to angles into a string, e.g. 3 into "3"
			if (ndefocusmax > 1) myformat += "_%0" + string(ndig) + "zd"; //construct format string for inserting two 0-padded angle indexes into file names - see usage below
			else myformat += "%0" + string(ndig) + "zd";
		}
	}
	else
	{
		if (ndefocusmax > 1)
		{
			myformat = "%zd";
		}
		if (nangles > 1)
		{
			if (ndefocusmax > 1) myformat += "_%zd";
			else myformat = "%zd";
		}
	}

	index_t ncurrent(0);
	idot = (int)filenamebase0.rfind('.');
	for (index_t i = nstart; i < nangles + nstart; i++)
	{
		for (index_t j = nstart1; j < vndefocus[i - nstart] + nstart1; j++)
		{
			outfilename_j = filenamebase0;
			if (ndefocusmax == 1 && nangles > 1) sprintf(buffer, myformat.data(), i);
			else if (ndefocusmax > 1 && nangles == 1) sprintf(buffer, myformat.data(), j);
			else sprintf(buffer, myformat.data(), j, i);
			outfilename_j.insert(idot, buffer);
			output[ncurrent++] = outfilename_j;
		}
	}
}
