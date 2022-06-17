
/*
 * pdb.h      : reading molecule info from a pdb file
 * A.V. Martin Jan 2013
 *
 * http://www.wwpdb.org/documentation/format33/v3.3.html
 *
 * Contains :
 *            definition of pdbdata structure
 *
 *
 */

/*
  ATOM Record Format

  COLUMNS        DATA  TYPE    FIELD        DEFINITION
  -------------------------------------------------------------------------------------
  1 -  6        Record name   "ATOM  "
  7 - 11        Integer       serial       Atom  serial number.
  13 - 16        Atom          name         Atom name.
  17             Character     altLoc       Alternate location indicator.
  18 - 20        Residue name  resName      Residue name.
  22             Character     chainID      Chain identifier.
  23 - 26        Integer       resSeq       Residue sequence number.
  27             AChar         iCode        Code for insertion of residues.
  31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
  39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
  47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
  55 - 60        Real(6.2)     occupancy    Occupancy.
  61 - 66        Real(6.2)     tempFactor   Temperature  factor.
  77 - 78        LString(2)    element      Element symbol, right-justified.
  79 - 80        LString(2)    charge       Charge  on the atom.
*/





#ifndef PDB_H
#define PDB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
  
  char rname[7];
  int  serial;
  char atomName[5];
  char altLoc;
  char resName[4];
  char chainID;
  int  resSeq;
  char iCode;
  double x;
  double y;
  double z;
  double occupancy;
  double tempFactor;
  char element[3];
  char charge[3];

} atomdata;

typedef struct{
  
  int natoms;

  atomdata * adata;

} pdbdata;

/*
 *  Initilize values of a pdbdata structure
 *  Values that require pdbcheck are set to 0
 *  Arrays that require pdbcheck are set to NULL
 */
void pdbdata_init( pdbdata * p);

/*
 *  Read all the atom information into a pdbdata structure
 */
int read_pdb( char * pdbname, pdbdata * pd );


/*
 *   Put info for each atom line int pdbdata
 */
int pdb_atomdata_from_file(char * pdbfname, pdbdata * pd);
int data_from_VestaXYZfile(char* pdbfname, pdbdata* pd);
int data_from_KirklandXYZfile(char* pdbfname, pdbdata* pd, float* pxlen = nullptr, float* pylen = nullptr, float* pzlen = nullptr);


/*
 * Counts the number of ATOM lines in the file
 */
int pdb_check(char * pdbfname, int * natoms);

/*
 * Allocate the memory for pdbdata.adata
 */
void pdballocate(pdbdata * pd);

/*
 * Free the memory for pdbdata.adata
 */
void pdbfree(pdbdata * pd);

void initialise_atomdata( atomdata * ad );

/*
 *  Returns 1 if the first word of the line is ATOM, and 0 otherwise
 */
int isATOMline(char * line);

/*
 *  Reads a line by tokenizing using strtok
 *  If it is an ATOM line, it prints each token
 */
void readline( char * line);


/*
 *  Reads an ATOM line and stores relevant information
 */
void readATOMline( char * line, atomdata * ad);

/*
 * print all the information contained in pdfdata.adata
 */
void print_pdb_data( pdbdata * pd );

/*
 *  Count the number of atoms of different types
 */
void pdb_chemical_composition( pdbdata * pd );

/*
 *  Sort ATOM lines in the ascending order according to element names
 */
void pdb_bubbleSort(pdbdata* pd);

/*
 *  Sort ATOM lines in the descending order according to atomic numbers
 */
void pdb_bubbleSort1(pdbdata* pd, int* ia);

/*
 *  Sort ATOM lines in the ascending order according to z coordinate
 */
int pdb_bubbleSort2(pdbdata* pd, int iFirst, int iLast);

/*
 *  Sort ATOM lines in the descending order according to occupancies
 */
int pdb_bubbleSort3(pdbdata* pd, int iFirst, int iLast);

/*
 *  Sort ATOM lines in the ascending order according to occupancies
 */
int pdb_bubbleSort3a(pdbdata* pd, int iFirst, int iLast);

/*
 *  Compare two element names
 */
int pdb_greater(char a[3], char b[3]);

/*
 *  Translate element names into atomic numbers
 */
int pdb_atomnumbers(pdbdata* pd, int ia[]);

/*
*  Translate atomic number into the element symbol
*/
int pdb_symbols(pdbdata* pd, int ia[]);

/*
*  Calculate L2 distance between two atomdata points
*/
inline double AtomDist(atomdata adata, atomdata bdata)
{
	return sqrt((adata.x - bdata.x) * (adata.x - bdata.x) + (adata.y - bdata.y) * (adata.y - bdata.y) + (adata.z - bdata.z) * (adata.z - bdata.z));
}

/*
*  Calculate square L2 distance between two atomdata points
*/
inline double AtomDist2(atomdata adata, atomdata bdata)
{
	return (adata.x - bdata.x) * (adata.x - bdata.x) + (adata.y - bdata.y) * (adata.y - bdata.y) + (adata.z - bdata.z) * (adata.z - bdata.z);
}

#endif
