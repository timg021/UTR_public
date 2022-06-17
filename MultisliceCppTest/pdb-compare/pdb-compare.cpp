#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <xstring>
#include <vector>
#include "Hungarian.h"
#include "pdb.h"

using namespace std;

int main(int argc, char* argv[])
{
	char pdbfile[1024], pdbfile1[1024];
	char outfile[1024];
	char cfileinfo[1024];
	char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024]; // auxiliary storage

	// Read input parameter file pdb.txt
	printf("\nStarting pdb-compare program ...\n");

	std::string sInputParamFile("pdb-compare.txt");
	if (argc > 1) sInputParamFile = argv[1];
	FILE* ffpar = fopen(sInputParamFile.c_str(), "rt");
	if (!ffpar)
	{
		printf("\nError: cannot open parameter file %s!\n", sInputParamFile.c_str());
		return -1;
	}
	else printf("\nReading input parameter file %s ...\n", sInputParamFile.c_str());

	// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
	while (true)
	{
		fgets(cline, 1024, ffpar);
		if (!(cline[0] == '/' && cline[1] == '/')) break;
	}

	// line 1
	strtok(cline, "\n"); // 1nd parameter: input test file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile) != 2)
	{
		printf("\n!!!Error reading input test file name from input parameter file.");
		return -1;
	}

	// line 2
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 2nd parameter: input test file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input test file type from input parameter file.");
		return -1;
	}
	int nfiletype = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype = atoi(cparam);
	if (nfiletype != 0 && nfiletype != 1 && nfiletype != 2)
	{
		printf("\n!!!Unknown input test file type in input parameter file.");
		return -1;
	}

	// line 3
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 3rd parameter: Euler rotation angles in degrees
	if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4)
	{
		printf("\n!!!Error reading Euler rotation angles from input parameter file.");
		return -1;
	}
	double anglez = -atof(cparam);
	double angley = -atof(cparam1);
	double anglez2 = -atof(cparam2);
	if (nfiletype != 2 && (anglez || angley || anglez2))
		printf("\n!!!WARNING: Euler rotation angles will be ignored, as they can only be applied to input files in Kirkland XYZ format.");
	else
		if (nfiletype == 2 && (anglez || angley || anglez2))
		{
			printf("\nThe structure will be rotated by %g degrees around the z axis.", anglez);
			printf("\nThe structure will be rotated by %g degrees around the y' axis.", angley);
			printf("\nThe structure will be rotated by %g degrees around the z'' axis.", anglez2);
		}
	anglez *= 3.141592653589793 / 180.0;
	angley *= 3.141592653589793 / 180.0;
	anglez2 *= 3.141592653589793 / 180.0;

	// line 4
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 4th parameter: input reference file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile1) != 2)
	{
		printf("\n!!!Error reading input reference file name from input parameter file.");
		return -1;
	}

	// line 5
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 5th parameter: input reference file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input reference file type from input parameter file.");
		return -1;
	}
	int nfiletype1 = 0; // input file type 0 - for PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file.
	nfiletype1 = atoi(cparam);
	if (nfiletype1 != 0 && nfiletype1 != 1 && nfiletype1 != 2)
	{
		printf("\n!!!Unknown input reference file type in input parameter file.");
		return -1;
	}

	// line 6
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //6th parameter: maximum distance in Angstroms below which a found atom match is considered acceptable
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading maximum match distance from input parameter file.");
		return -1;
	}
	double DistMax = 0.0;
	DistMax = atof(cparam);
	if (DistMax <= 0)
	{
		printf("\n!!!Error: maximum match distance must be positive.");
		return -1;
	}
	double DistMax2 = DistMax * DistMax;

	// line 7
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 7th parameter: output data sort
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading output data sort type from input parameter file.");
		return -1;
	}
	int noutsort = 0; // output file sort 0 - no sort, 1 - sort by ascending order of z coordinate, 2- sort by distances between the matched test and reference atoms in the ascending order
	noutsort = atoi(cparam);
	switch (noutsort)
	{
	case 0: printf("\nThe atoms in the input structure will be sorted by the atomic numbers in the descending order."); break;
	case 1: printf("\nThe atoms in the input structure will be sorted by the z coordinate in the ascending order."); break;
	case 2: printf("\nThe atoms in the input structure will be sorted by the distances between the matched test and reference atoms in the ascending order."); break;
	default: printf("\n!!!Unknown output data sort type in the input parameter file."); return -1;
	}

	// line 8
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 8th parameter: Hungarian algorithm switch
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading Hungarian algorithm switch from input parameter file.");
		return -1;
	}
	int nHungarian = 0; // 0 = don't use Hungarian algorithm, 1 = use it
	nHungarian = atoi(cparam);
	switch (nHungarian)
	{
	case 0: printf("\nHungarian algorithm won't be used for atom matching."); break;
	case 1: printf("\nHungarian algorithm will be used for atom matching."); break;
	default: printf("\n!!!Unknown Hungarian algorithm switch value in the input parameter file."); return -1;
	}

	// line 9
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 8th parameter: output file name
	if (sscanf(cline, "%s %s", ctitle, outfile) != 2)
	{
		printf("\n!!!Error reading output file name from input parameter file.");
		return -1;
	}

	// line 10
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 9th parameter: free form first line for the output file
	if (sscanf(cline, "%s %s", ctitle, cfileinfo) != 2)
	{
		printf("\n!!!Error reading free-form line for the output file from input parameter file.");
		return -1;
	}

	fclose(ffpar);


	// read input file0 (test structure)
	printf("\n");
	float xctblength, yctblength, zctblength;
	pdbdata pd;
	pdbdata_init(&pd);
	if (nfiletype == 0)
	{
		printf("\nReading input file1 %s in PDB format ...", pdbfile);
		if (read_pdb(pdbfile, &pd) == -1) return -1; // read PDB file
	}
	else if (nfiletype == 1)
	{
		printf("\nReading input file1 %s in Vesta XYZ export format ...", pdbfile);
		if (data_from_VestaXYZfile(pdbfile, &pd) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype == 2)
	{
		printf("\nReading input file1 %s in Kirkland XYZ format ...", pdbfile);
		if (data_from_KirklandXYZfile(pdbfile, &pd, &xctblength, &yctblength, &zctblength) == -1) return -1; // read Kirkland XYZ file

		if ((anglez || angley || anglez2) && (xctblength <= 0 || yctblength <= 0 || zctblength <= 0))
		{
			printf("\n!!!ERROR: CT box sizes are not all positive in the input file %s.", pdbfile);
			return -1;
		}
	}
	if (pd.natoms < 2)
	{
		printf("\n!!!ERROR: fewer than 2 atoms in the input file %s.", pdbfile);
		return -1;
	}


	// read input file1 (reference structure)
	pdbdata pd1;
	pdbdata_init(&pd1);
	if (nfiletype1 == 1)
	{
		printf("\nReading input file2 %s in Vesta XYZ export format ...", pdbfile1);
		if (data_from_VestaXYZfile(pdbfile1, &pd1) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype1 == 0)
	{
		printf("\nReading input file2 %s in PDB format ...", pdbfile1);
		if (read_pdb(pdbfile1, &pd1) == -1) return -1; // read PDB file
	}
	else if (nfiletype1 == 2)
	{
		printf("\nReading input file2 %s in Kirkland XYZ format ...", pdbfile1);
		if (data_from_KirklandXYZfile(pdbfile1, &pd1) == -1) return -1; // read Kirkland XYZ file
	}
	if (pd1.natoms < 2)
	{
		printf("\n!!!ERROR: fewer than 2 atoms in the input file %s.", pdbfile1);
		return -1;
	}

	//translate element symbols into atomic numbers
	int* ia = (int*)malloc(pd.natoms * sizeof(int));
	if (nfiletype != 2)
	{
		if (pdb_atomnumbers(&pd, ia))
		{
			printf("\n!!!Error encountered while finding atomic number for a given element name!!!\n");
			return -1;
		}
	}
	else
		for (int i = 0; i < pd.natoms; i++) ia[i] = pd.adata[i].serial;

	int* ia1 = (int*)malloc(pd1.natoms * sizeof(int));
	if (nfiletype1 != 2)
	{
		if (pdb_atomnumbers(&pd1, ia1))
		{
			printf("\n!!!Error encountered while finding atomic number for a given element name!!!\n");
			return -1;
		}
	}
	else
		for (int i = 0; i < pd1.natoms; i++) ia1[i] = pd1.adata[i].serial;

	// check that the reference data does not contain hydrogen atoms (matching of hydrogen atoms is considered counter-productive since the test positions will be "wasted" on them)
	for (int i = 0; i < pd1.natoms; i++)
		if (ia1[i] == 1)
		{
			printf("\n!!!Error: hydrogen atom found in the reference data - matching of hydrogen atoms in not supported!!!\n");
			return -1;
		}

	// sort entries by atom number in descending order
	pdb_bubbleSort1(&pd, ia);
	pdb_bubbleSort1(&pd1, ia1);

	// open the output file early in order to start writing some results there as they are obtained
	FILE* ff = fopen(outfile, "wt");
	if (ff == NULL)
	{
		printf("\nERROR: cannot open %s file for writing!!!", outfile);
		return -2;
	}
	fprintf(ff, "For each atom from %s file, locating the closest atom in %s file.\n", pdbfile, pdbfile1);
	fprintf(ff, "Each entry contains the following data:\n");
	fprintf(ff, "1st column contains atomic numbers (reference data)\n");
	fprintf(ff, "2nd, 3rd, and 4th columns contain x, y and z atomic coordinates, respectively (reference data)\n");
	fprintf(ff, "5th column contains the L2 distance between atoms in the test and the matched reference data.\n");
	fprintf(ff, "6th column contains atomic numbers (matched test data)\n");
	fprintf(ff, "7th, 8th, and 9th columns contain x, y and z atomic coordinates, respectively (matched test data)\n");
	fprintf(ff, "10th column contains original data for the test file\n");
	fprintf(ff, "%s", cfileinfo); // free-form file info line
	fprintf(ff, "\n");

	// clean the entries for the future matched indexes
	for (int j = 0; j < pd.natoms; j++) pd.adata[j].tempFactor = -1; // this parameter will contain the index of the found best match
	for (int i = 0; i < pd1.natoms; i++) pd1.adata[i].tempFactor = -1; // this parameter will contain the index of the found best match

	// rotate the molecule
	if (anglez || angley || anglez2)
	{
		printf("\nRotating the test structure ...");
		double cosanglez = cos(anglez);
		double sinanglez = sin(anglez);
		double cosangley = cos(angley);
		double sinangley = sin(angley);
		double cosanglez2 = cos(anglez2);
		double sinanglez2 = sin(anglez2);

		double xc = 0.5 * xctblength;
		double yc = 0.5 * yctblength;
		double zc = 0.5 * zctblength;

		double xk, yk, zk;

		for (int i = 0; i < pd.natoms; i++)
		{
			xk = pd.adata[i].x;
			yk = pd.adata[i].y;
			zk = pd.adata[i].z;

			//rotation around Z axis
			pd.adata[i].x = xc + (xk - xc) * cosanglez + (-yk + yc) * sinanglez;
			pd.adata[i].y = yc + (xk - xc) * sinanglez + (yk - yc) * cosanglez;

			//rotation around Y' axis
			xk = pd.adata[i].x;
			pd.adata[i].x = xc + (xk - xc) * cosangley + (zk - zc) * sinangley;
			pd.adata[i].z = zc + (-xk + xc) * sinangley + (zk - zc) * cosangley;

			//rotation around Z" axis
			xk = pd.adata[i].x;
			yk = pd.adata[i].y;
			pd.adata[i].x = xc + (xk - xc) * cosanglez2 + (-yk + yc) * sinanglez2;
			pd.adata[i].y = yc + (xk - xc) * sinanglez2 + (yk - yc) * cosanglez2;
		}
	}

	printf("\nCalculating distances between atoms in the %s test file ...", pdbfile);
	double disttot(0), rmintot = AtomDist(pd.adata[0], pd.adata[1]);
	for (int i = 0; i < pd.natoms - 2; i++)
	{
		double r, rmin = AtomDist(pd.adata[i], pd.adata[i + 1]);
		for (int j = i + 2; j < pd.natoms; j++)
		{
			r = AtomDist(pd.adata[i], pd.adata[j]);
			if (r < rmin) rmin = r;
		}
		if (rmin < rmintot) rmintot = rmin;
		disttot += rmin;
	}
	disttot /= double(pd.natoms - 1);
	printf("\nAverage distance between closest pairs of atoms in the test structure from %s = %g.", pdbfile, disttot);
	printf("\nMinimum distance between closest pairs of atoms in the test structure from %s = %g.", pdbfile, rmintot);

	printf("\nCalculating distances between atoms in the %s reference file ...", pdbfile1);
	disttot = 0.0; rmintot = AtomDist(pd1.adata[0], pd1.adata[1]);
	for (int i = 0; i < pd1.natoms - 2; i++)
	{
		double r, rmin = AtomDist(pd1.adata[i], pd1.adata[i + 1]);
		for (int j = i + 2; j < pd1.natoms; j++)
		{
			r = AtomDist(pd1.adata[i], pd1.adata[j]);
			if (r < rmin) rmin = r;
		}
		if (rmin < rmintot) rmintot = rmin;
		disttot += rmin;
	}
	disttot /= double(pd1.natoms - 1);
	printf("\nAverage distance between closest pairs of atoms in the reference structure from %s = %g.", pdbfile1, disttot);
	printf("\nMinimum distance between closest pairs of atoms in the reference structure from %s = %g.", pdbfile1, rmintot);


	printf("\n\nCalculating all pair-wise distances between atoms in %s test file and those in %s reference file ...", pdbfile, pdbfile1);
	vector< vector<double> > mR2(pd1.natoms);
	for (int i = 0; i < pd1.natoms; i++) mR2[i].resize(pd.natoms);
	for (int i = 0; i < pd1.natoms; i++)
	{
		for (int j = 0; j < pd.natoms; j++)
		{
			mR2[i][j] = AtomDist2(pd1.adata[i], pd.adata[j]);
		}
	}

	// for each atom in the reference structure pd1 find the closest atom in the test structure pd, and save the results in the new "double" structure pd2 
	if (nHungarian)
	{
		//THIS IS VERY SLOW - may be OK for 100 atoms or so, but definitely not for 10,000
		printf("\nUsing Hungarian algorithm to find the optimal bipartate matching ...");
		HungarianAlgorithm HungAlgo;
		vector<int> assignment;
		double cost = HungAlgo.Solve(mR2, assignment);

		for (int i = 0; i < pd1.natoms; i++)
		{
			pd1.adata[i].tempFactor = assignment[i]; // index of the closest atom from the test strcture
			pd1.adata[i].occupancy = sqrt(mR2[i][assignment[i]]); // distance to the closest atom
		}
	}
	else
	{
		// this is a simple and fast, but imperfect, alternative to optimal bipartate matching
		printf("\nRunning a simple sequential search for closest atom pairs ...");
		for (int i = 0; i < pd1.natoms; i++)
		{
			int jmin = -1;
			double r2min = std::numeric_limits<double>::max();
			for (int j = 0; j < pd.natoms; j++)
			{
				if (pd.adata[j].tempFactor != -1) continue; // if this test atom has been already matched, skip it
				if (mR2[i][j] < r2min)
				{
					r2min = mR2[i][j];
					jmin = j;
				}
			}
			if (r2min < DistMax2) // if the found minimal distance is smaller than the defined upper limit, then mark the match
			{
				pd1.adata[i].tempFactor = jmin; // index of the closest atom from the test strcture
				pd1.adata[i].occupancy = sqrt(mR2[i][jmin]); // distance to the closest atom
				pd.adata[jmin].tempFactor = i; // mark this test atom as already matched
			}
		}
	}

	// calculate average distance and std
	int nodup = 0;
	double adist = 0.0, maxdist = 0.0;
	for (int i = 0; i < pd1.natoms; i++)
	{
		if (pd1.adata[i].tempFactor == -1) continue;
		nodup++;
		adist += pd1.adata[i].occupancy;
		if (pd1.adata[i].occupancy > maxdist) maxdist = pd1.adata[i].occupancy;
	}
	if (nodup == 0) adist = 0.0; else adist /= (double)nodup;
	double stddist = 0.0;
	for (int i = 0; i < pd1.natoms; i++)
	{
		if (pd1.adata[i].tempFactor == -1) continue;
		stddist += (pd1.adata[i].occupancy - adist) * (pd1.adata[i].occupancy - adist);
	}
	if (nodup <= 1) stddist = 0.0; else stddist = sqrt(stddist / ((double)(nodup - 1)));
	printf("\n%d out of total %d atoms in the reference structure (i.e. %g percent) have been uniquely matched with atoms of the test structure.", nodup, pd1.natoms, (double)nodup / (double)pd1.natoms * 100.0);
	printf("\nThe result contains %d false negatives (%g percent) and %d false positives (%g percent).", pd1.natoms - nodup, 100.0 * (double)(pd1.natoms - nodup) / (double)pd1.natoms, pd.natoms - nodup, 100.0 * (double)(pd.natoms - nodup) / (double)pd.natoms);
	printf("\nAverage distance between the matched test and template atoms = %g.", adist);
	printf("\nMaximum distance between the matched test and template atoms = %g.", maxdist);
	printf("\nStandard deviation of the distance between the matched test and template atoms = %g.", stddist);
	printf("\n");
	fprintf(ff, "\n%d atoms (out of total %d atoms) from the reference structure (i.e. %g percent) have been uniquely matched with atoms of the test structure.", nodup, pd1.natoms, (double)nodup / (double)pd1.natoms * 100.0);
	fprintf(ff, "\nThe result contains %d false negatives (%g percent) and %d false positives (%g percent).", pd1.natoms - nodup, 100.0 * (double)(pd1.natoms - nodup) / (double)pd1.natoms, pd.natoms - nodup, 100.0 * (double)(pd.natoms - nodup) / (double)pd.natoms);
	fprintf(ff, "\nAverage distance between the matched test and template atoms = %g.", adist);
	fprintf(ff, "\nMaximum distance between the matched test and template atoms = %g.", maxdist);
	fprintf(ff, "\nStandard deviation of the distance between the matched test and template atoms = %g.", stddist);
	fprintf(ff, "\n");

	// list all atoms from the reference file which have not been matched by any of the atoms from the test file
	for (int i = 0; i < pd1.natoms; i++)
		if (pd1.adata[i].tempFactor == -1)
		{
			//printf("\n@@@ Reference atom no. %d, atom type = %d, atom positions = (%g, %g, %g) has not been matched.", i, ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z);
			fprintf(ff, "\n@@@ Reference atom no. %d, atom type = %d, atom positions = (%g, %g, %g) has not been matched.", i, ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z);
		}

	// sort output records: 0 - no sort, 1 - sort by ascending order of z coordinate, 2- sort by distances between the matched test and reference atoms in the ascending order
	// NOTE that we don't need to worry about sorting the test structure, because after the reference data is sorted, the matched atoms of the test structure will still be linked through pd1.adata[i].tempFactor
	if (noutsort == 1)
	{
		printf("\nSorting the output in ascending order of z coordinate in the reference structure ...");
		pdb_bubbleSort2(&pd1, 0, pd.natoms);
	}
	else if (noutsort == 2)
	{
		printf("\nSorting the output by distances between the matched reference and test atoms in ascending order ...");
		pdb_bubbleSort3a(&pd1, 0, pd.natoms);
	}

	// output to the target file
	printf("\nWriting output comparison file %s ...", outfile);
	fprintf(ff, "\n");
	int jmin;
	for (int i = 0; i < pd1.natoms; i++)
	{
		if (pd1.adata[i].tempFactor == -1) continue;
		jmin = (int)pd1.adata[i].tempFactor;
		fprintf(ff, "\n%d %f %f %f %f", ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z, pd1.adata[i].occupancy); // reference atom first
		fprintf(ff, " %d %f %f %f %f ", ia[jmin], pd.adata[jmin].x, pd.adata[jmin].y, pd.adata[jmin].z, pd.adata[jmin].occupancy); // matched nearest test atom second
	}
	fclose(ff);
	free(ia1);
	free(ia);

	printf("\nFinished!\n");

	return 0;
}