#include <stdexcept>
#include "AddIce.h"
#include "slicelib.hpp"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// The code in this file is the result of some minor modifications of the code borrowed from the module pdb2xyz.cpp
// from the project temsim of the computem solution created by E.J. Kirkland (https://sourceforge.net/projects/computem/),
// see the relevant copyright / usage conditions at https://sourceforge.net/projects/computem/
// No changes to copyright, usage or ownership are claimed here. In particular, no warranties are given whatsoever and any
// use of this code will be entirely at the risk of the user

const int NVAL = 5;
enum { OX = 0, OY = 1, OZ = 2, Oocc = 3, Oznum = 4 };
const int NSMAX = 1000;   // max number of slices
const int NCMAX = 1024;   // max characters in file names
const int NZMAX = 103;    // max atomic number Z
const double RMIN = 1.4;  // min separation distance of atoms in ice layer (in Angstroms)

void hydrocoord(std::vector<double>& ctest, std::vector<double>& coord, double ax, double by, double cz, unsigned long* piseed);
bool testcoord(std::vector<double>& ctest, std::vector< std::vector<double> >& coord,
	const int ncoord, int* pos, const double rmin);
bool testInside(std::vector<double>& ctest,
	std::vector<double>& pdbmin, std::vector<double>& pdbmax);
void insertcoord(std::vector<double>& ctest, std::vector<std::vector<double> >& coord,
	const int ncoord, const int pos);
void newcoord(std::vector<double>& c,
	const double xmax, const double ymax, const double zmax, unsigned long* piseed);
void sortByZd(float* x, float* y, float* z, float* occ,
	int* Znum, int natom);


/*------------------------ fillIce() ------------------------*/
/*
	fill a vol. with amorphous ice as random coord.
	in the range (0,0,0) to (ax,by,cz)

	add ice mol. everywhere avoiding occupied positions
	until density/number outside PDB bounding box is correct
	(i.e. can't calc. ice density around rough edges of PDB structure
	but assume its the same as in vac. if every position is equally
	added to)

	start with existing PDB coord. and work around it

  coord[][] = start with existing list of PDB atoms sorted by z = coord[][OZ]
				will merge in a set of ICE coord to fill rest of vol.
				dimensions ncoord x NVAL
  nMolIce = number of ice mol. coordinates outside of PDB vol to generate
  np = number of start coord from PDB atoms
  rmin = minimum separation distance

  ax*by*cz = total vol
  (xmin,ymin,zmin) -> (xmax,ymax,zmax) = pdbmin -> pdbmax
		= vol of PDB to exclude for ice density calculation

  its hard to calculate vol of PDB with random edges so exclude
  a bounding box enclosing the PDB

  assumed globals
	 NVAL

*/
int fillIce(std::vector< std::vector<double> >& coord, const int nMolIce, const int np,
	double ax, double by, double cz,
	std::vector<double>& pdbmin, std::vector<double>& pdbmax, double rmin, unsigned long* piseed)
{
	int ic, i, ncoord;
	std::vector<double> ctest(NVAL), hcoord(NVAL);

	// initial number of coord. = number in PDB structure
	ncoord = np;    //  total number of atoms PDB+ice
	ic = 1;         // number of ice mol. added
	//@@@@@@@@@@@@@@ TEMP
	//printf("\n@@@@@@ nMolIce = %d", nMolIce);
	do {
		// get new coordinate of ice mol (actually only the O atom)
		newcoord(ctest, ax, by, cz, piseed);
		ctest[Oocc] = 1.0; //  occupancy = 1.0
		ctest[Oznum] = 8;   // Z=8 for oxygen
		//cout << "new coord= " << ctest[0] << ", " << 
		//           ctest[1] << ", " <<ctest[2] << endl;   // ??? diagnostic

		//  add it to the list if its OK - this also sorts
		if (testcoord(ctest, coord, ncoord, &i, RMIN) == true) {
			insertcoord(ctest, coord, ncoord, i);
			ncoord++;

			// count mol. outside PDB bounding box
			if (0 == testInside(ctest, pdbmin, pdbmax)) {
				ic++;
				//if (ic % 1000 == 0) cout << "\r nmol= " << ic << flush;
				//@@@@@@@@@@@@@@ TEMP
				//printf("\n@@@ ic = %d", ic);
			}

			//  add 2 x H atoms at random positions near O atom
			//  should constrain bond angle but ignore for simplicity
			//   insert is slow - might be faster to insert 3 at one time (more code needed)
			hcoord[Oocc] = 1.0;    //  occupancy = 1.0
			hcoord[Oznum] = 1;      //  Z=1 for hydrogen
			hydrocoord(hcoord, ctest, ax, by, cz, piseed);
			insertcoord(hcoord, coord, ncoord, i);
			ncoord++;
			hydrocoord(hcoord, ctest, ax, by, cz, piseed);
			insertcoord(hcoord, coord, ncoord, i);
			ncoord++;
		}
	} while (ic < nMolIce);

	return ncoord;

}  // end fillIce()


/*------------------------ fillSolid() ------------------------*/
/*
	fill a vol. with an amorphous solid with random coord.
	in the range (0,0,0) to (ax,by,cz)

  coord[][] = to get list of existing sorted by z = coord[][OZ]
				dimensions ncoord x NVAL
  ncoord = number of coordinates to generate
  rmin = minimum separation distance

  assumed globals
	 NVAL

*/
int fillSolid(std::vector< std::vector<double> >& coord, const int ncoord,
	double ax, double by, double cz, double rmin, unsigned long* piseed)
{
	int ic, i;
	std::vector<double> ctest(NVAL);

	newcoord(coord[0], ax, by, cz, piseed);   // start coord. list
	ic = 1;
	do {
		// get new coordinate 
		newcoord(ctest, ax, by, cz, piseed);
		//cout << "new coord= ", << ctest[0] << ", " << 
		//           ctest[1] << ", " <<ctest[2] << endl;

		//  add it to the list if its OK - this also sorts
		if (testcoord(ctest, coord, ic, &i, RMIN) == true) {
			insertcoord(ctest, coord, ic, i);
			ic++;
			//if (ic % 1000 == 0) cout << "\r natom= " << ic << flush;
		}
	} while (ic < ncoord);

	return ic;

}  // end fillSolid()

/*------------------------ hydrocoord() ------------------------*/
/*
	generate coordinates for Hydrogen around Oxygen atom

  NOTE: this ignores possible overlap of H atoms
		and the bonding angle of H atoms

  H bond length is < min. sep. distance so overlap unlikely

  ctest[] = to get H coord.
  coord[] = position of O atom

 */
void hydrocoord(std::vector<double>& ctest, std::vector<double>& coord, double ax, double by, double cz, unsigned long* piseed)
{
	long i;
	double dx, d = 0.0;
	double ax1 = ax - 0.1, by1 = by - 0.1, cz1 = cz - 0.1;

	/* Hydrogen-Oxygen bond length in Angstroms */
	static const double bondlength = 0.958;

	//  get random x,y,z offset at spec. bond length
	for (i = 0; i < 3; i++) {
		ctest[i] = dx = ranflat(piseed) - 0.5;
		d += dx * dx;
	}
	d = bondlength / sqrt(d);

	for (i = 0; i < 3; i++) {
		ctest[i] = coord[i] + ctest[i] * d;
		if (ctest[i] < 0.1) ctest[i] = 0.1; //@@@@ TEG code: we don't allow negative coordinates
	}
	if (ctest[0] > ax1) ctest[0] = ax1; //@@@@ TEG code: we don't allow H atoms to stick out of the box
	if (ctest[1] > by1) ctest[1] = by1; //@@@@ TEG code: we don't allow H atoms to stick out of the box
	if (ctest[2] > cz1) ctest[2] = cz1; //@@@@ TEG code: we don't allow H atoms to stick out of the box

	//  probably should generate a 2nd H coord at fixed bond angle here (???)

	return;
} /* end hydrocoord() */

/*------------------------ insertcoord() ------------------------*/
/*
	insert a new coordinate in the list

  ctest[] = new coord to insert
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = position to insert coord. at

  insert new coord at index pos and move all the rest down one

  assumed globals
	 NVAL

	 might use std::vector<> insert here sometime - needs an iterator though
*/
void insertcoord(std::vector<double>& ctest, std::vector<std::vector<double> >& coord,
	const int ncoord, const int pos)
{
	int i, j, n;

	//  check that there is enough memory left
	n = (int)coord.size();
	if (ncoord >= n) {
		throw std::runtime_error("Error: too big a number of new coordinates in insertcoord().");
		//cout << "ncoord = " << ncoord << " too big in insertcoord(), exit...." << endl;
		//exit(0);
	}

	for (i = ncoord; i > pos; i--) {
		for (j = 0; j < NVAL; j++)
			coord[i][j] = coord[i - 1][j];
	}

	for (j = 0; j < NVAL; j++)
		coord[pos][j] = ctest[j];

} // end insertcoord() 


/*------------------------ newcoord() ------------------------*/
/*
	generate a new random coordinate inside the required volume

	xmax, ymax, zmax = volume size

*/
void newcoord(std::vector<double>& c,
	const double xmax, const double ymax, const double zmax, unsigned long* piseed)
{
	c[OX] = xmax * ranflat(piseed);
	c[OY] = ymax * ranflat(piseed);
	c[OZ] = zmax * ranflat(piseed);
	return;

} // end newcoord() 

/*----------------- sortByZd() ------------------------------

from slicelib.cpp but converted from float to double

	improved Shell sort modeled after prog. 6.5 (pg. 274) of
	R. Sedgewick, "Algorithms in C", 3rd edit. Addison-Wesley 1998

	x[], y[], z[]   = atom coordinates
	occ[]           = occupancy of each atom
	Znum[]          = atomic number of each atom
	natom           = number of atoms
*/
void sortByZd(float* x, float* y, float* z, float* occ,
	int* Znum, int natom)
{
	int i, j, h, Znum2;
	//float x2, y2, z2, occ2;
	float x2, y2, z2, occ2;

	for (h = 1; h <= (natom - 1) / 9; h = 3 * h + 1);

	for (; h > 0; h /= 3)
		for (i = h; i < natom; i++) {
			j = i;
			x2 = x[i];
			y2 = y[i];
			z2 = z[i];
			occ2 = occ[i];
			Znum2 = Znum[i];
			while ((j >= h) && (z2 < z[j - h])) {
				x[j] = x[j - h];
				y[j] = y[j - h];
				z[j] = z[j - h];
				occ[j] = occ[j - h];
				Znum[j] = Znum[j - h];
				j -= h;
			}
			x[j] = x2;
			y[j] = y2;
			z[j] = z2;
			occ[j] = occ2;
			Znum[j] = Znum2;
		}

	/*  Test sort routine -- DELETE this after awhile  */
	for (i = 1; i < natom; i++) {
		if (z[i - 1] > z[i]) throw std::runtime_error("Error: bad sort in sortByZd().");
	}

}  /* end sortByZ() */


/*------------------------ testcoord() ------------------------*/
/*
	test the new coordinate to see if its too close to an existing
	coordinate (closer than RMIN)

	A straight search of the whole list is very slow (proportional to N^2)
	and was not pratical for a large set of coordinates (i.e. took
	an absurd amount of CPU time).

	This routines assumes that the list of existing coordinates is
	sorted wrt one coord (z in this case).  Once the position of the new
	coordinate is located in the list (using a binary search), then
	only the small range of coord. near this point need to be
	tested.  This makes the test dramatically faster.

  ctest[] = new coord to test
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = returned position to insert coord.

  assumed globals
	 OX, OY, OZ, RMIN

*/
bool testcoord(std::vector<double>& ctest, std::vector< std::vector<double> >& coord,
	const int ncoord, int* pos, const double rmin)
{
	/* this is the more -sophisticated version sorted by Z */
	long i, j, k;
	double d, dx, dy, dz, dz2, rmin2, z, range;

	/* ---- find postion of ctest[] in coord[][] using a binary search  ---
		 i should get the position to insert (between i and i+1)
		 this assumes coord[][] is sorted by coord OZ
	*/
	//  printf("testcoord() top, ncoord= %d\n", ncoord );

	z = ctest[OZ];
	if (z <= coord[0][OZ]) i = j = 0;
	else if (z >= coord[ncoord - 1][OZ]) i = j = ncoord;
	else {
		i = 0;
		j = ncoord - 1;
		do {
			k = (i + j) / 2;
			if (z < coord[k][OZ])  j = k;
			else if (z >= coord[k][OZ]) i = k;
		} while ((j - i) > 1);
	}

	/* now that we have the position of the new point
	   we only have to explore within RMIN of this point */

	rmin2 = rmin * rmin;
	range = 2.0 * rmin2;   /* add a little safety margin */
	k = j;
	while (k >= 0) {
		dx = ctest[OX] - coord[k][OX];
		dy = ctest[OY] - coord[k][OY];
		dz = ctest[OZ] - coord[k][OZ];
		dz2 = dz * dz;
		d = dx * dx + dy * dy + dz2;
		if (d <= rmin2) return(false);
		if (dz2 > range) break;
		k--;
	}
	k = j;
	while (k < ncoord) {
		dx = ctest[OX] - coord[k][OX];
		dy = ctest[OY] - coord[k][OY];
		dz = ctest[OZ] - coord[k][OZ];
		dz2 = dz * dz;
		d = dx * dx + dy * dy + dz2;
		if (d <= rmin2) return(false);
		if (dz2 > range) break;
		k++;
	}

	/* if it gets to here then this is a good point */

	*pos = j;

	return(true);

} // end testcoord() 

/*------------------------ testInside() ------------------------*/
/*
  test if position ctest=(x,y,z) is inside
 (xmin,ymin,zmin) -> (xmax,ymax,zmax) = pdbmin -> pdbmax
		= vol of PDB to exclude for ice density calculation

*/
bool testInside(std::vector<double>& ctest,
	std::vector<double>& pdbmin, std::vector<double>& pdbmax)
{
	if ((ctest[OX] > pdbmin[OX]) && (ctest[OX] < pdbmax[OX])
		&& (ctest[OY] > pdbmin[OY]) && (ctest[OY] < pdbmax[OY])
		&& (ctest[OZ] > pdbmin[OZ]) && (ctest[OZ] < pdbmax[OZ])) return (true);
	else return(false);
}


//-------------------------------------------------------------------------------------------------
int AddIce(float iceThick, float ctblength, int natom, int** pZnum, float** px, float** py, float** pz, float** pocc, float** pwobble, float wobbleaver, unsigned long* piseed)
// this function assumes that the input structure comes enclosed in the cube ctblength^3
{
	// recenter the molecule in the thicker z-slab
	float zoff = float(0.5 * iceThick - 0.5 * ctblength);
	int* Znum(*pZnum);
	float* x(*px), * y(*py), * z(*pz), * occ(*pocc), * wobble(*pwobble);
	for (int i = 0; i < natom; i++) z[i] += zoff;

	sortByZd(x, y, z, occ, Znum, natom);

	const double ICE_WEIGHT = 18.0152;     // molecular weight for water (in gm/mole) 
	const double ICE_DENSITY = 0.9167;    //  density in gm/cm^3
	const double ICE_FILL_FRACTION = 1.0;  // filling fraction for this density
	const double RMIN = 1.4;          // min separation distance of C (in Angstroms)
	const double NAV = 6.0225e23;     // Avagadro's number (#/mole)

	double density = ICE_FILL_FRACTION * NAV * ICE_DENSITY * (1.0e-24) / ICE_WEIGHT; //  # atoms/Angs^3 
	float ax = ctblength;
	float by = ctblength;
	float cz = iceThick;

	//  allocate memory for much more than needed - enough for all ice plus PDB
	//    not efficient but easier to program (sorry?)
	int nMolIce = (int)(ax * by * cz * density);  //  total number of ice molecules without PDB
	int nAtomsIce = 3 * nMolIce;  //  total number of ice atoms without PDB
	int ntotal = 2 * nAtomsIce + natom;  //  should be much more than needed  
	std::vector<double> tmp(NVAL);
	std::vector< std::vector<double> > coord(ntotal, tmp); // bigger than needed 

	//  copy PDB coord. into common format to use non-overlap subroutines 
	float xmin(x[0]), xmax(x[0]);
	float ymin(y[0]), ymax(y[0]);
	float zmin(z[0]), zmax(z[0]);
	for (int i = 0; i < natom; i++)
	{
		coord[i][OX] = x[i];
		coord[i][OY] = y[i];
		coord[i][OZ] = z[i];
		coord[i][Oocc] = occ[i];
		coord[i][Oznum] = Znum[i];
		if (x[i] < xmin) xmin = x[i];
		if (x[i] > xmax) xmax = x[i];
		if (y[i] < ymin) ymin = y[i];
		if (y[i] > ymax) ymax = y[i];
		if (z[i] < zmin) zmin = z[i];
		if (z[i] > zmax) zmax = z[i];
	}
	if (zmax > cz) throw std::runtime_error("Error: z-extent of the molecule is larger than the thickness of the ice layer.");

	//  fill the whole vol. with ice excluding the PDB aotm positions
	//    by adding H2O to existing PDB coord.
	//  ice may go into the holes in the PDB so only the density
	//  outside the PDB vol can be calculated.

	double vol = ax * by * cz - (zmax - zmin) * (ymax - ymin) * (xmax - xmin);  // vol excluding PDB bounding box
	nMolIce = (int)(vol * density);  //  total number of ice mol outside PDB bound. box

	std::vector<double> pdbmin(NVAL), pdbmax(NVAL);
	pdbmin[OX] = xmin;
	pdbmin[OY] = ymin;
	pdbmin[OZ] = zmin;
	pdbmax[OX] = xmax;
	pdbmax[OY] = ymax;
	pdbmax[OZ] = zmax;

	//  now do all the work 
	natom = fillIce(coord, nMolIce, natom, ax, by, cz, pdbmin, pdbmax, RMIN, piseed);

	//  copy PDB coord. back into the original simple format
	free(x); free(y); free(z); free(occ); free(wobble); free(Znum);
	x = (float*)malloc1D(natom, sizeof(float), "x"); *px = x;
	y = (float*)malloc1D(natom, sizeof(float), "y"); *py = y;
	z = (float*)malloc1D(natom, sizeof(float), "z"); *pz = z;
	occ = (float*)malloc1D(natom, sizeof(float), "occ"); *pocc = occ;
	wobble = (float*)malloc1D(natom, sizeof(float), "wobble"); *pwobble = wobble;
	Znum = (int*)malloc1D(natom, sizeof(int), "Znum"); *pZnum = Znum;
	for (int i = 0; i < natom; i++) {
		x[i] = (float)coord[i][OX];
		y[i] = (float)coord[i][OY];
		z[i] = (float)coord[i][OZ];
		occ[i] = (float)coord[i][Oocc];
		Znum[i] = (int)coord[i][Oznum];
		wobble[i] = wobbleaver;
	}

return natom;
}

//--------------------------------------------------------------
int AddCarbon(double cThick, double cWidth, pdbdata& pd, unsigned long* piseed, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax)
// Adds a 'support' parallelepiped with (x, y, z) dimensions ("cWidth", "cWidth", "cThick") filled with randomly positioned atoms of carbon,
// then updates the molecular structure with the added carbon atoms.
// The support parallelepiped is placed immediately below the molecule and centred with respect to the molecule's transverse position
// "natom" is the original number of atoms in the molecule
// returns "ntotal" which is equal to the sum of the original "natom" and the number of added carbon atoms
{
	const double NAV = 6.0225e23;     // Avagadro's number (#/mole)
	//  CFILL_FRACTION=1 makes flat support with almost no structure
	//      make <1 to get a little surface roughness
	const double CWEIGHT = 12.01115;       // molec. weight of carbon
	const double CDENSITY = 2.0;           // density in gm/cm^3 approx. for amorphous C
	const double CFILL_FRACTION = 0.9;     // filling fraction for this density
	const double RMIN = 1.4;			   // min separation distance of C (in Angstroms)

	// fill the substrate parallelepiped with carbon atoms
	double density = CFILL_FRACTION * NAV * CDENSITY * (1.0e-24) / CWEIGHT; //  # atoms/Angs^3 
	double ax = cWidth;
	double by = cWidth;
	int ncarbon = (int)(ax * by * cThick * density);
	std::vector<double> tmp(NVAL);
	std::vector< std::vector<double> > coord(ncarbon, tmp);
	ncarbon = fillSolid(coord, ncarbon, ax, by, cThick, RMIN, piseed);

	// shift the support to the bottom centre of the molecule
	double xc = 0.5 * (xmin + xmax);
	double yc = 0.5 * (ymin + ymax);
	double xshift = xc - 0.5 * cWidth; // the transverse centre of the substrate should coincide with the transverse centre of the molecule
	double yshift = yc - 0.5 * cWidth; // the transverse centre of the substrate should coincide with the transverse centre of the molecule
	//double zshift = zmin - cThick - RMIN; // the highest atom of the substrate should sit at (zmin - RMIN)
	double zshift = zmin - cThick - 1.0; // the highest atom of the substrate should sit at (zmin - 1.0)
	for (int i = 0; i < ncarbon; i++)
	{
		coord[i][OX] += xshift;
		coord[i][OY] += yshift;
		coord[i][OZ] += zshift;
	}

	//  copy PDB coord. of the molecule's atoms into common format
	int natom = pd.natoms;
	std::vector<double> tmp0(NVAL);
	std::vector< std::vector<double> > coord0(natom, tmp0);
	std::vector<std::string> elem(natom);
	for (int i = 0; i < natom; i++)
	{
		coord0[i][OX] = pd.adata[i].x;
		coord0[i][OY] = pd.adata[i].y;
		coord0[i][OZ] = pd.adata[i].z;
		coord0[i][Oocc] = pd.adata[i].occupancy;
		coord0[i][Oznum] = pd.adata[i].serial;
		elem[i] = pd.adata[i].element;
	}

	//  copy PDB coord. back into the original simple format
	int ntotal = natom + ncarbon;
	pd.natoms = ntotal;
	pdballocate(&pd);
	for (int i = 0; i < natom; i++) {
		pd.adata[i].x = coord0[i][OX];
		pd.adata[i].y = coord0[i][OY];
		pd.adata[i].z = coord0[i][OZ];
		pd.adata[i].occupancy = coord0[i][Oocc];
		pd.adata[i].serial = (int)coord0[i][Oznum];
		strcpy(pd.adata[i].element, elem[i].c_str());
	}
	for (int i = 0; i < ncarbon; i++) {
		int i1 = i + natom;
		pd.adata[i1].x = coord[i][OX];
		pd.adata[i1].y = coord[i][OY];
		pd.adata[i1].z = coord[i][OZ];
		pd.adata[i1].occupancy = 1.0;
		pd.adata[i1].serial = 6;
		strcpy(pd.adata[i1].element, "C\0");
	}

	// update the bounds of the new joint structure
	if (xc - 0.5 * cWidth < xmin) xmin = xc - 0.5 * cWidth;
	if (xc + 0.5 * cWidth > xmax) xmax = xc + 0.5 * cWidth;
	if (yc - 0.5 * cWidth < ymin) ymin = yc - 0.5 * cWidth;
	if (yc + 0.5 * cWidth > ymax) ymax = yc + 0.5 * cWidth;
	zmin -= cThick;

	return ntotal;
}
