#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <ctime>
#include "pdb.h"
#include "AddIce.h"
#include "Seb.h"
#include "XA_file.h"
#include "Hungarian.h"

// This program reads a PDB, Vesta XYZ or Kirkland XYZ file, 
// optionally centers the position of the "structure" within a given "enclosing cube" [0, ctblength] x [0, ctblength] x [0, ctblength],
// optionally rotates the "structure" by a given angle around the z and y' axes, 
// optionally removes "duplicate" atomes, and
// outputs the data in the form of a muSTEM, Vesta XYZ or Kirkland XYZ file, while it also
// sorts all atoms in descending order with respect to atom numbers, and sorts all atoms with the same number in ascending order with respect to z coordinate
//
// This programs uses the "miniball" code written by Martin Kutz <kutz@math.fu-berlin.de> and Kaspar Fischer <kf@iaeth.ch> for finding the minimal ball enclosing the structure
// See https://github.com/hbf/miniball for authors, theory and copyright information

using namespace xar;

void CreatePdbFromOrientations(pdbdata& pd, vector<Pair> v2angles);

int main(int argc, char* argv[])
{
	printf("\nStarting pdb program ...");
	int i, j, k;
	char pdbfile[1024]; // input file name
	char outfile[1024]; // output file name
	char cfileinfo[1024]; // free-form 1st line in the output file

	char cline[1024], ctitle[1024], cparam[1024], cparam1[1024]; // auxiliary storage

	// Read input parameter file pdb.txt
	std::string sInputParamFile("pdb.txt");
	if (argc > 1) sInputParamFile = argv[1];
	FILE* ffpar = fopen(sInputParamFile.c_str(), "rt");
	if (!ffpar)
	{
		printf("\nError: cannot open parameter file %s!\n", sInputParamFile.c_str());
		return -1;
	}
	else printf("\nReading input parameter file %s ...", sInputParamFile.c_str());


	// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
	while (true)
	{
		fgets(cline, 1024, ffpar);
		if (!(cline[0] == '/' && cline[1] == '/')) break;
	}

	// line 1
	strtok(cline, "\n"); // 1nd parameter: input file name
	if (sscanf(cline, "%s %s", ctitle, pdbfile) != 2)
	{
		printf("\n!!!Error reading input file name from input parameter file.");
		return -1;
	}
	printf("\nInput structure file = %s", pdbfile);

	// line 2
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 2nd parameter: input file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading input file type from input parameter file.");
		return -1;
	}
	int nfiletype = 0; // input file type 0 - PDB input file, 1 - Vesta XYZ input file, 2 - Kirkland XYZ file, 3 - text file with orientations/defocuses.
	nfiletype = atoi(cparam);
	switch (nfiletype)
	{
	case 0: printf("\nInput structure file format = PDB."); break;
	case 1: printf("\nInput structure file format = Vesta_XYZ."); break;
	case 2: printf("\nInput structure file format = Kirkland_XYZ."); break;
	case 3: printf("\nInput structure file format = text file with orientations/defocuses."); break;
	default: printf("\n!!!Unknown input file type in input parameter file."); return -1;
	}

	// line 3
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //3rd parameter: optional carbon support layer thickness and transversal extent
	if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3)
	{
		printf("\n!!!Error reading carbon support layer thickness and width from input parameter file.");
		return -1;
	}
	double cThick(-1), cWidth(-1);
	cThick = atof(cparam); cWidth = atof(cparam1);
	if (cThick > 0)
	{
		printf("\nCarbon support layer thickness = %g, width and length = %g (Angstroms).", cThick, cWidth);
		if (cWidth <= 0)
		{
			printf("\n!!!carbon support width %g is not positive in pdb.txt!!!", cWidth);
			return -1;
		}
	}
	else printf("\nNo carbon support layer will be created.");

	// line 4
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //4th parameter: CT cube side length
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading CT cube side length from input parameter file.");
		return -1;
	}
	double ctblength = 0; // CT box length
	ctblength = atof(cparam);
	printf("\nEnclosing cube side length = %g (Angstroms).", ctblength);

	// line 5
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //5th parameter: rotation angle around z
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading rotation angle around z axis from input parameter file.");
		return -1;
	}
	double anglez = 0; // rotation angle 
	anglez = -atof(cparam);
	printf("\nThe structure, together with the optional carbon support, will be rotated by %g degrees around the z axis.", anglez);
	anglez *= 3.141592653589793 / 180.0;
	double cosanglez = cos(anglez);
	double sinanglez = sin(anglez);

	// line 6
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //6th parameter: rotation angle around y' (i.e. around y axis after the initial rotation around z axis)
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading rotation angle around y' axis from input parameter file.");
		return -1;
	}
	double angley = 0; // rotation angle 
	angley = -atof(cparam);
	printf("\nThe structure, together with the optional carbon support, will be rotated by %g degrees around the y' axis.", angley);
	angley *= 3.141592653589793 / 180.0;
	double cosangley = cos(angley);
	double sinangley = sin(angley);

	// line 7 
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //7th parameter: rotation angle around z" (this rotation can also be done later within the 2D projected images)
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading rotation angle around z'' axis from input parameter file.");
		return -1;
	}
	double anglez2 = 0; // rotation angle 
	anglez2 = -atof(cparam);
	printf("\nThe structure, together with the optional carbon support, will be rotated by %g degrees around the z'' axis.", anglez2);
	anglez2 *= 3.141592653589793 / 180.0;
	double cosanglez2 = cos(anglez2);
	double sinanglez2 = sin(anglez2);

	// line 8
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 8th parameter: apply the miniball centre and diameter
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading the centre in the miniball option from input parameter file.");
		return -1;
	}
	int nApplyMiniball = 0; // 0 - do not shift the structure, 1 - shift the structure to centre it in the miniball
	nApplyMiniball = atoi(cparam);
	switch (nApplyMiniball)
	{
	case 0: printf("\nThe structure will not be shifted."); break;
	case 1: printf("\nThe structure will be shifted to centre it in the determined minimal containing ball."); break;
	default: printf("\n!!!Unknown value for the centring in the miniball option in the input parameter file."); return -1;
	}

	// line 9
	fgets(cline, 1024, ffpar); strtok(cline, "\n");  //9th parameter: maximum distance to remove duplicates
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading maximum distance to remove duplicates from input parameter file.");
		return -1;
	}
	double dupdist = 0.0;
	dupdist = atof(cparam);
	if (dupdist < 0) dupdist = 0.0;
	if (dupdist > 0) printf("\nAny atom located closer than %g Angstroms from its neighbours will be removed from the structure.", dupdist);
	else printf("\nNo checks of the inter-atomic distance will be applied.");
	double dupdist2 = dupdist * dupdist;

	// line 10
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 10th parameter: output data sort
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading output data sort type from input parameter file.");
		return -1;
	}
	int noutsort = 0; // output file sort 0 - no sort, 1 - sort by ascending order of z coordinate, 2- sort by descending order of the occupancy values
	noutsort = atoi(cparam);
	switch (noutsort)
	{
	case 0: printf("\nThe atoms in the input structure will be sorted by the atomic numbers in the descending order."); break;
	case 1: printf("\nThe atoms in the input structure will be sorted by the z coordinate in the ascending order."); break;
	case 2: printf("\nThe atoms in the input structure will be sorted by the occupancy parameter in the descending order."); break;
	default: printf("\n!!!Unknown output data sort type in the input parameter file."); return -1;
	}

	// line 11
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 11th parameter: output file name
	if (sscanf(cline, "%s %s", ctitle, outfile) != 2)
	{
		printf("\n!!!Error reading output file name from input parameter file.");
		return -1;
	}
	printf("\nThe output file name = %s.", outfile);

	// line 12
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 12th parameter: output file type
	if (sscanf(cline, "%s %s", ctitle, cparam) != 2)
	{
		printf("\n!!!Error reading output file type from input parameter file.");
		return -1;
	}
	int noutfiletype = 0; // output file type 0 -  for Kirkland-format XYZ file, 1 - for muSTEM-foram XTL file
	noutfiletype = atoi(cparam);
	switch (noutfiletype)
	{
	case 0: printf("\nOutput structure file format = muSTEM."); break;
	case 1: printf("\nOutput structure file format = Vesta_XYZ."); break;
	case 2: printf("\nOutput structure file format = Kirkland_XYZ."); break;
	default: printf("\n!!!Unknown output file type in input parameter file."); return -1;
	}

	// line 13
	fgets(cline, 1024, ffpar); strtok(cline, "\n"); // 13th parameter: free form first line for the output file
	if (sscanf(cline, "%s %s", ctitle, cfileinfo) != 2)
	{
		printf("\n!!!Error reading free-form line for the output file from input parameter file.");
		return -1;
	}
	printf("\nThe free-form information line in the output file will be %s.", cfileinfo);
	printf("\n");

	fclose(ffpar); // close input parameter file


	// read input PDB, Vesta or Kirkland XYZ file
	pdbdata pd;
	pdbdata_init(&pd);
	if (nfiletype == 0)
	{
		printf("\nReading input file %s in PDB format ...", pdbfile);
		if (read_pdb(pdbfile, &pd) == -1) return -1; // read PDB file
	}
	else if (nfiletype == 1)
	{
		printf("\nReading input file %s in Vesta export XYZ format ...", pdbfile);
		if (data_from_VestaXYZfile(pdbfile, &pd) == -1) return -1; // read Vesta export XYZ file
	}
	else if (nfiletype == 2)
	{
		printf("\nReading input file %s in Kirkland XYZ format ...\n", pdbfile);
		if (data_from_KirklandXYZfile(pdbfile, &pd) == -1) return -1; // read Kirkland XYZ file
	}
	else if (nfiletype == 3)
	{
		printf("\nReading input file %s in text format ...\n", pdbfile);
		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;
		vector<Pair> v2shifts;
		vector<Pair> vastigm;
		try
		{
			if (GetFileExtension(string(pdbfile)) == string(".TXT"))
				ReadDefocusParamsFile(pdbfile, v2angles, vvdefocus, false);
			else
				if (GetFileExtension(string(pdbfile)) == string(".RELIONNEW"))
					ReadRelionDefocusParamsFile(pdbfile, v2angles, vvdefocus, vastigm, v2shifts, false);
				else
				{
					printf("\n!!!Error: orientation text file can only have .txt or .RELLIONNEW filename extension.");
					return -1;
				}
		}
		catch (std::runtime_error& E)
		{
			printf("\n\n!!!Exception: %s\n", E.what());
			return -1;
		}

		index_t nangles = v2angles.size(); // number of rotation steps 
		CreatePdbFromOrientations(pd, v2angles);
	}
	//print_pdb_data(&pd);

	// find the spatial bounds of the structure
	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin = xmax = pd.adata[0].x;
	ymin = ymax = pd.adata[0].y;
	zmin = zmax = pd.adata[0].z;

	for (i = 1; i < pd.natoms; i++)
	{
		if (pd.adata[i].x < xmin) xmin = pd.adata[i].x;
		else if (pd.adata[i].x > xmax) xmax = pd.adata[i].x;
		if (pd.adata[i].y < ymin) ymin = pd.adata[i].y;
		else if (pd.adata[i].y > ymax) ymax = pd.adata[i].y;
		if (pd.adata[i].z < zmin) zmin = pd.adata[i].z;
		else if (pd.adata[i].z > zmax) zmax = pd.adata[i].z;
	}
	double xsize = xmax - xmin;
	double ysize = ymax - ymin;
	double zsize = zmax - zmin;
	printf("\nThere were %d atoms in the input structure file.", pd.natoms);
	printf("\nThe x, y and z extent of the structure in the input file were:");
	printf("\n  [%g, %g], [%g, %g], and [%g, %g], respectively.", xmin, xmax, ymin, ymax, zmin, zmax);

	// translate element symbols into atomic numbers
	int* ia = (int*)malloc(pd.natoms * sizeof(int));
	if (nfiletype == 2 || nfiletype == 3)
	{
		for (i = 0; i < pd.natoms; i++) ia[i] = pd.adata[i].serial;
		pdb_symbols(&pd, ia);
	}
	else if (pdb_atomnumbers(&pd, ia))
	{
		printf("\n!!!Error encountered while finding atomic number for a given element name!!!\n");
		return -1; // the function pdb_atomnumbers will print its error messages before exiting
	}

	// optionally add carbon support layer
	if (cThick > 0)
	{
		unsigned long iseed = std::abs((long)time(NULL));
		int natoms0 = pd.natoms;
		int ntotal = AddCarbon(cThick, cWidth, pd, &iseed, xmin, xmax, ymin, ymax, zmin, zmax);
		printf("\nThis program has added a carbon support layer with %d atoms.", ntotal - natoms0);
		printf("\nThe combined structure with the included carbon support now contains %d atoms.", ntotal);
		// update atom numbers array
		free(ia); ia = (int*)malloc(pd.natoms * sizeof(int)); pdb_atomnumbers(&pd, ia);
		xsize = xmax - xmin;
		ysize = ymax - ymin;
		zsize = zmax - zmin;
		printf("\nThe x, y and z extent of the structure is now:");
		printf("\n  [%g, %g], [%g, %g], and [%g, %g], respectively.", xmin, xmax, ymin, ymax, zmin, zmax);
	}

	// Miniball code - find the minimum 3D ball containing the molecule (surprisingly non-trivial!)
	// See https://github.com/hbf/miniball for authors, theory and copyright information
	typedef double FT;
	typedef Seb::Point<FT> Point;
	typedef std::vector<Point> PointVector;
	typedef Seb::Smallest_enclosing_ball<FT> Miniball;
	PointVector S;
	std::vector<double> coords(3);
	for (int i = 0; i < pd.natoms; ++i)
	{
		coords[0] = pd.adata[i].x;
		coords[1] = pd.adata[i].y;
		coords[2] = pd.adata[i].z;
		S.push_back(Point(3, coords.begin()));
	}
	Miniball mb(3, S);
	double diam = floor(2.0 * mb.radius()) + 1.0; // we want to round the diameter up and ensure that atoms won't stick out after rotation due to numerical precision
	Miniball::Coordinate_iterator center_it = mb.center_begin();
	double xc0 = center_it[0];
	double yc0 = center_it[1];
	double zc0 = center_it[2];
	#ifdef _DEBUG
	printf("\n");
	mb.verify();
	#endif
	printf("\nThe centre of the minimal ball enclosing the structure is located at the point x = %g, y = %g, z = %g.", xc0, yc0, zc0);
	printf("\nThe diameter of this minimal ball is %g.", diam);


	// TEG: we make all coordinates non-negative, as it seems that Kirkland's software expects this,
	// and we also centre each of XYZ coordinates within the interval [0, ctblength], so that the sample can remain within the cube during rotation
	// and we rotate the structure around the z axis by the anglez and around the y axis by the angley set in the input parameter file

	if (ctblength < diam)
	{
		if (nApplyMiniball)
		{

			ctblength = diam;
			printf("\n!!! Enclosing_cube_side_length parameter in the input parameter file was too small and has been replaced with %g.", diam);
		}
		else
		{
			printf("\n!!!Error: enclosing_cube_side_length parameter in the input parameter file is smaller than the diameter of the minimum ball enclosing the structure.");
			return -1;
		}
	}

	xsize = ysize = zsize = ctblength;

	double xc = 0.5 * ctblength;
	double yc = 0.5 * ctblength;
	double zc = 0.5 * ctblength;

	double xshift(0.0), yshift(0.0), zshift(0.0);
	if (nApplyMiniball)
	{
		xshift = xc - xc0;
		yshift = yc - yc0;
		zshift = zc - zc0;
		printf("\nThe structure will be shifted into the enclosing cube with side length %g and positive coordinates.", ctblength);
		printf("\nThe new centre of the minimal ball containing the structure will be located at the point x = %g, y = %g, z = %g.\n", xc, yc, zc);
	}
	else
	{
		printf("\nThe structure will NOT be shifted.");
		printf("\nThe enclosing cube has the side length %g.", ctblength);
		printf("\nThe centre of rotation is x = %g, y = %g, z = %g.\n", xc, yc, zc);
	}

	// center and rotate the molecule
	double xk, yk, zk;
	for (i = 0; i < pd.natoms; i++)
	{
		// centering
		xk = pd.adata[i].x + xshift;
		yk = pd.adata[i].y + yshift;
		zk = pd.adata[i].z + zshift;

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

	// sort entries by atom number in descending order
	pdb_bubbleSort1(&pd, ia);

	// find the number of different atom types (this assumes that all entries have been sorted in descending order)
	int natypes = 1;
	for (i = 0; i < pd.natoms - 1; i++)
		if (ia[i + 1] != ia[i]) natypes++;

	// find the number of atoms for each type
	int* vtypes = (int*) malloc(natypes * sizeof(int));
	j = 0;
	vtypes[0] = 1;
	for (i = 0; i < pd.natoms - 1; i++)
	{
		if (ia[i + 1] != ia[i]) { j++; vtypes[j] = 1; }
		else vtypes[j]++; // count in the next(!) atom
	}

	// print out the numbers of atoms for each atom type present in the input file
	printf("\nThe structure contains %d different types of atoms", natypes);
	k = 0;
	for (j = 0; j < natypes; j++)
		for (i = 0; i < vtypes[j]; i++)
		{
			if (i == 0) printf("\nThere are %d %s atoms (with atomic number = %d) in the structure.", vtypes[j], pd.adata[k].element, ia[k]);
			k++;
		}
	printf("\n");

	// sort atoms of each type by z coordinate in ascending order or occupancy values in descending order
	if (noutsort == 1)
	{
		k = 0;
		for (j = 0; j < natypes; j++)
		{
			if (pdb_bubbleSort2(&pd, k, k + vtypes[j]) != 0) return -1;
			k += vtypes[j];
		}
	}
	else if (noutsort == 2)
	{
		k = 0;
		for (j = 0; j < natypes; j++)
		{
			if (pdb_bubbleSort3(&pd, k, k + vtypes[j]) != 0) return -1;
			k += vtypes[j];
		}
	}

    // find and remove duplicates
	double r2;
	for (i = 0; i < pd.natoms; i++)
	{
		if (pd.adata[i].tempFactor == -777) continue; // if this atom has already been marked as duplicate, skip it
		for (j = i + 1; j < pd.natoms; j++)
		{
			if (pd.adata[j].tempFactor == -777) continue; // don't compare atom with those already classified as duplicates
			r2 = (pd.adata[i].x - pd.adata[j].x) * (pd.adata[i].x - pd.adata[j].x) + (pd.adata[i].y - pd.adata[j].y) * (pd.adata[i].y - pd.adata[j].y) + (pd.adata[i].z - pd.adata[j].z) * (pd.adata[i].z - pd.adata[j].z);
			if (r2 < dupdist2) // atom j is too close to atom i
			{
				if (pd.adata[i].occupancy < pd.adata[j].occupancy) 
				{
					pd.adata[i].tempFactor = -777; // mark this atom as duplicate
					printf("\n!!! Found duplicate atom entries no. %d and %d (0 based), distance = %g Angstroms.", i, j, sqrt(r2));
					printf("\n!!!   Atom no. %d is removed", i);
					break;
				}  
				else
				{
					pd.adata[j].tempFactor = -777; // mark this atom as duplicate
					printf("\n!!! Found duplicate atom entries no. %d and %d (0 based), distance = %g Angstroms.", i, j, sqrt(r2));
					printf("\n!!!   Atom no. %d is removed", j);
				}
			}
		}
	}

	// count the number of non-duplicate atoms
	int natoms1 = 0;
	for (i = 0; i < pd.natoms; i++) if (pd.adata[i].tempFactor != -777) natoms1++;
	int nduplicates = pd.natoms - natoms1;
	printf("\nIn total, %d duplicate atoms found in the structure.", nduplicates);

	// create structure containing the output structure with removed duplicates
	pdbdata pd1;
	pdbdata_init(&pd1);
	pd1.natoms = natoms1;
	pdballocate(&pd1);
	int* ia1 = (int*) malloc(pd1.natoms * sizeof(int));
	j = 0;
	for (i = 0; i < pd.natoms; i++) 
		if (pd.adata[i].tempFactor != -777)
		{
			ia1[j] = ia[i];
			pd1.adata[j] = pd.adata[i];
			j++;
		}

	// find the number of different atom types (this assumes that all entries have been sorted in descending order)
	int natypes1 = 1;
	for (i = 0; i < pd1.natoms - 1; i++)
		if (ia1[i + 1] != ia1[i]) natypes1++;

	// find the number of atoms for each type
	int* vtypes1 = (int*) malloc(natypes1 * sizeof(int));
	j = 0;
	vtypes1[0] = 1;
	for (i = 0; i < pd1.natoms - 1; i++)
	{
		if (ia1[i + 1] != ia1[i]) { j++; vtypes1[j] = 1; }
		else vtypes1[j]++; // count in the next(!) atom
	}

	// print out the numbers of atoms for each atom type present after the duplicated have been exluded
	if (nduplicates > 0)
	{
		k = 0;
		for (j = 0; j < natypes1; j++)
			for (i = 0; i < vtypes1[j]; i++)
			{
				if (i == 0) printf("\nThere are now %d %s atoms (with atomic number = %d) in the structure.", vtypes1[j], pd1.adata[k].element, ia1[k]);
				k++;
			}
	}
	printf("\n");


	// output to the target file
	FILE* ff = fopen(outfile, "wt");
	if (noutfiletype == 0) // muSTEM format output
	{
		printf("\nWriting output file %s in muSTEM XTL format ...", outfile);

		fprintf(ff, "%s\n", cfileinfo); // free-form file info line
		fprintf(ff, "%f %f %f 90.0 90.0 90.0\n", xsize, ysize, zsize); // if npositive == 1, ctblength is written as the unit cell dimension
		fprintf(ff, "%f \n", 200.0); // probe energy
		fprintf(ff, "%d \n", natypes1); // number of atoms types

		k = 0;
		for (j = 0; j < natypes1; j++)
		{
			fprintf(ff, "%s\n", pd1.adata[k].element); // next atom type info
			fprintf(ff, "%d %f 1.0 0.0\n", vtypes1[j], (double)ia[k]); // number of atoms of this type, atom number for this type
			for (i = 0; i < vtypes1[j]; i++)
			{
				fprintf(ff, "%f %f %f\n", pd1.adata[k].x / ctblength, pd1.adata[k].y / ctblength, pd1.adata[k].z / ctblength); // fractional coordinates of an atom
				k++;
			}
		}
		
		free(vtypes1);
	}
	else if (noutfiletype == 1) // Vesta XYZ output
	{
		printf("\nWriting output file %s in Vesta XYZ format ...", outfile);
		fprintf(ff, "%d\n", pd1.natoms);
		fprintf(ff, "%s\n", cfileinfo); // free-form file info line
		for (i = 0; i < pd1.natoms; i++)
		{
			if (pd1.adata[i].occupancy <= 0 || pd1.adata[i].occupancy > 1) pd1.adata[i].occupancy = 1.0;
			fprintf(ff, "%s %f %f %f %f\n", pd1.adata[i].element, pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z, pd1.adata[i].occupancy);
		}
	}
	else if (noutfiletype == 2) // Kirkland XYZ output
	{
		printf("\nWriting output file %s in Kirkland's XYZ format ...", outfile);
		fprintf(ff, "%s\n", cfileinfo); // free-form file info line
		fprintf(ff, "%f %f %f\n", xsize, ysize, zsize); // if npositive == 1, ctblength is written as the unit cell dimension
		for (i = 0; i < pd1.natoms; i++)
		{
			//if (pd1.adata[i].occupancy <= 1.e-7 || pd1.adata[i].occupancy > 1) pd1.adata[i].occupancy = 1.0;
			fprintf(ff, "%d %f %f %f %f\n", ia1[i], pd1.adata[i].x, pd1.adata[i].y, pd1.adata[i].z, pd1.adata[i].occupancy);
		}
		fprintf(ff, "%i\n\n\n", -1);
	}

	fclose(ff);
	free(ia);

	printf("\nFinished!\n");

	return 0;
}


void CreatePdbFromOrientations(pdbdata& pd, vector<Pair> v2angles)
{
	pd.natoms = (int)v2angles.size();
	pdballocate(&pd);

	vector<float> x(v2angles.size());
	vector<float> y(v2angles.size());
	vector<float> z(v2angles.size());

	float xc(0.0f), yc(0.0f), zc(0.0f);
	float xxx, yyy, zzz;

	for (index_t na = 0; na < v2angles.size(); na++)
	{
		Pair angle = v2angles[na];
		angle.a *= PI180; angle.b *= PI180;
		float sinZ = (float)sin(angle.a), cosZ = (float)cos(angle.a);
		float sinY = (float)sin(angle.b), cosY = (float)cos(angle.b);

		x[na] = y[na] = 0.0f; z[na] = 1.0f; // initial position of the illumination vector prior to rotation

		// rotation around Y axis (by the Y' rotation angle)
		xxx = xc + (x[na] - xc) * cosY + (-z[na] + zc) * sinY;
		zzz = zc + (x[na] - xc) * sinY + (z[na] - zc) * cosY;
		x[na] = xxx; z[na] = zzz;

		// rotation around Z axis (by the Z rotation angle)
		xxx = xc + (x[na] - xc) * cosZ + (y[na] - yc) * sinZ;
		yyy = yc + (-x[na] + xc) * sinZ + (y[na] - yc) * cosZ;
		x[na] = xxx; y[na] = yyy;

		pd.adata[na].x = x[na];
		pd.adata[na].y = y[na];
		pd.adata[na].z = z[na];
		pd.adata[na].serial = 1;
		pd.adata[na].atomName[0] = 'H'; pd.adata[na].atomName[1] = '\0';
		pd.adata[na].occupancy = 1.0;
		pd.adata[na].element[0] = 'H'; pd.adata[na].element[1] = '\0';
	}
}