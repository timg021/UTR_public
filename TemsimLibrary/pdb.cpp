
/*
 *  pdb.h   : read a pdb file and store data in a structure
 *  A. Martin January 2013
 *
 */

#include "pdb.h"
#include <ctype.h>


/*
 * Routines to write:
 * read the pdb
 * check pdb for number of atoms (I have the feeling I wanted to check something else in this routine too
 * Initialize the pdbdata struct
 * Initilize pdbdata arrays (after checking for the number of atoms)
 */


/*
 *  Initilize values of a pdbdata structure
 *  Values that require pdbcheck are set to 0
 *  Arrays that require pdbcheck are set to NULL
 */
void pdbdata_init( pdbdata * pd)
{

  pd->natoms = 0;
  pd->adata = NULL;
}



/*
 *  Read all the atom information into a pdbdata structure
 */
int read_pdb( char * pdbname, pdbdata * pd )
{

  int natoms;

  /*
   *  do a pdb check to get the number of atoms
   */
  if (pdb_check( pdbname, &natoms) == -1)
	  return -1;
  
  /*
   *  Allocate memory
   */
  printf("\nNumber of atoms %d\n", natoms);
 
  pd->natoms = natoms;

  pdballocate( pd );


  /*
   * Read atom data into the pdbdata
   */
  pdb_atomdata_from_file(pdbname, pd);

  /*
   *  Print the data in pd
   */
  //print_pdb_data( pd );
  return 0;

}



/*
 *   Put info for each atom line int pdbdata
 */
int pdb_atomdata_from_file(char * pdbfname, pdbdata * pd)
{

  FILE * file = fopen ( pdbfname, "r" );
  char line[256];
  char line_temp[256];
  int i=0;
  line[0] = '\0';

  int count =0;

  if (file != NULL){
    while( fgets(line, sizeof line, file) != NULL) {
     
      //fputs ( line, stdout ); /* write the line */
      strcpy( line_temp, line );

      if (isATOMline(line_temp) == 1) {
	readATOMline( line, &(pd->adata[count]) );
	count ++;
      }
      i++;
      //printf("Line %d ; ATOM count %i\n", i, count);  
      
    }
   
    fclose(file);
  } 
  else {
    printf("pdb file %s was not found (for pdb_check).\n", pdbfname);
	return -1;
  }

  return 0;
}


/*
 *   //@@@@ added by TEG for reading Vesta-export XYZ files; 11 June 2019
 */
int data_from_VestaXYZfile(char* pdbfname, pdbdata* pd)
{

	FILE* file = fopen(pdbfname, "r");
	char line[256];
	float x, y, z, occup;
	char elem[16];

	int count = 0, natoms = 0;

	if (file != NULL) {
		
		fgets(line, sizeof line, file); // first line must contain the number of atoms
		sscanf(line, "%d", &natoms);
		printf("Number of atoms %d\n", natoms);
		pd->natoms = natoms;
		pdballocate(pd);

		fgets(line, sizeof line, file); // second line contains arbitrary description text

		while (fgets(line, sizeof line, file) != NULL) 
		{
			sscanf(line, "%s %f %f %f %f", &elem[0], &x, &y, &z, &occup);
			strncpy(pd->adata[count].element, elem, 2); pd->adata[count].element[2] = '\0';
			pd->adata[count].x = x;
			pd->adata[count].y = y;
			pd->adata[count].z = z;
			pd->adata[count].occupancy = occup;
			count++;
		}

		if (count != natoms)
		{
			printf("!!!wrong number of atom lines (in data_from_VestaXYZfile).\n");
			return -1;
		}

		fclose(file);
	}
	else {
		printf("!!!file %s was not found (in data_from_VestaXYZfile).\n", pdbfname);
		return -1;
	}

	return 0;
}


/*
 *   //@@@@ added by TEG for reading Kirkland XYZ files; 13 October 2019
 */
int data_from_KirklandXYZfile(char* pdbfname, pdbdata* pd, float* pxlen, float* pylen, float* pzlen)
{

	FILE* file = fopen(pdbfname, "r");
	char line[256];
	float x, y, z, occup;

	int count = 0, natoms = 0, anum, neof = 0;

	if (file != NULL) {

		while (fgets(line, sizeof line, file) != NULL)
		{
			sscanf(line, "%d", &neof);
			if (neof == -1) break;
			count++;
		}
		fclose(file);

		natoms = count - 2; // the first 2 lines in Kirkland XYZ file do not contain atom entries

		if (natoms < 1 || neof != -1)
		{
			printf("!!!problem reading file %s (in data_from_KirklandXYZfile).\n", pdbfname);
			return -1;
		}

		pd->natoms = natoms;
		pdballocate(pd);

		file = fopen(pdbfname, "r");

		fgets(line, sizeof line, file); // first line contains arbitrary info
		fgets(line, sizeof line, file); // second line contains XYZ sizes of the molecule
		if (pxlen == nullptr || pylen == nullptr || pzlen == nullptr)
			sscanf(line, "%f %f %f", &x, &y, &z);
		else
			sscanf(line, "%f %f %f", pxlen, pylen, pzlen);
        pd->adata[0].tempFactor = 0; //TEG!!!: I have no recollections why I introduced this line earlier. I am keeping it for now, since it seems harmless

		for (count = 0; count < natoms; count++)
		{
			fgets(line, sizeof line, file);
			sscanf(line, "%d %f %f %f %f", &anum, &x, &y, &z, &occup); // Kirkland XYZ file may also contain a column with Debye-Waller factor, but we ignore it
			pd->adata[count].serial = anum;
			pd->adata[count].x = x;
			pd->adata[count].y = y;
			pd->adata[count].z = z;
			pd->adata[count].occupancy = occup;
		}

		fclose(file);
	}
	else {
		printf("!!!file %s was not found (in data_from_KirklandXYZfile).\n", pdbfname);
		return -1;
	}

	return 0;
}


/*
 * Counts the number of ATOM lines in the file
 */

int pdb_check(char * pdbfname, int * natoms){

  FILE * file = fopen ( pdbfname, "r" );
  char line[256];
  int i=0;
  line[0] = '\0';

  int count =0;

  if (file != NULL){
    while( fgets(line, sizeof line, file) != NULL) {
     
      //fputs ( line, stdout ); /* write the line */
      
      if (isATOMline(line) == 1) {
	count ++;
      }
      i++;
      //printf("Line %d ; ATOM count %i\n", i, count);
           
    }
    
    *natoms = count;

    fclose(file);
  } 
  else {
    printf("pdb file %s was not found (for pdb_check).\n", pdbfname);
	return -1;
  }

  return 0;
}


/*
 *  allocate an array of atomdata structures in a pdbdata structure
 */
void pdballocate(pdbdata * pd)
{

  
  if (pd->adata != NULL){
    free(pd->adata);
    pd->adata = NULL;
    
    } 
 
  pd->adata = (atomdata * ) malloc( pd->natoms * sizeof(atomdata) );
  int i;
  for (i=0;i<pd->natoms;i++){
    initialise_atomdata( &(pd->adata[i]) );
  }

}

/*
 *  free data for pdbdata.adata
 */
void pdbfree(pdbdata * pd)
{

  free(pd->adata);
  pd->adata = NULL;
 
}

void initialise_atomdata( atomdata * ad ){

  ad->rname[0] = '\0';
  ad->atomName[0] = '\0';
  ad->resName[0] = '\0';
  ad->element[0] = '\0';
  ad->charge[0] = '\0';
  ad->serial = 0;
  ad->altLoc = '\0';
  ad->chainID = '\0';
  ad->x = 0;
  ad->y = 0;
  ad->z = 0;
  ad->occupancy = 0.0;
  ad->tempFactor = 0.0;
}


/*
 *  Returns 1 if the first word of the line is ATOM, and 0 otherwise
 */

int isATOMline(char * line)
{

  char * token;
  //token = "\0";

  if (strncmp(line, "HETATM", 6) == 0) return 1;

  token = strtok(line," =:\n"); 

  if (strcmp(token,"ATOM") == 0) {
    return 1;
  }
  else {
    return 0;
  }
 
}


/*
 *  Reads a line by tokenizing using strtok
 *  If it is an ATOM line, it prints each token
 */

void readline( char * line){

  char * token;
  //token = "\0";

  token = strtok(line," =:\n"); 

  if (strcmp(token,"ATOM") == 0) {
  
    while ( token != NULL){
      
      fputs ( token, stdout ); /* write the next part of the line */
      fputs ( " ", stdout ); /* space between tokens */
      token = strtok(NULL," =:\n"); 
    }
    fputs ( "\n", stdout ); /* new line */

  }
 

}




/*
 *  Reads an ATOM line and stores relevant information
 */
void readATOMline( char * line, atomdata * ad)
{

  int i;
  int pos, len;

  /*
   *  Copy the record name
   */
  pos = 0;
  len = 6;
  char rname[7];
  for (i=0;i<len;i++){
    rname[i] = line[pos+i];
  }
  pos += len;
  rname[len] = '\0';
  strcpy( ad->rname, rname );

  /*
   *  Copy the atom serial number
   */
  len = 5;
  char serial[6];
  for (i=0;i<len;i++){
    serial[i] = line[pos+i];
  }
  pos += len;
  serial[len] = '\0';
  ad->serial = atoi(serial);


  /*
   *  12th position is blank
   */
  pos ++;


  /*
   *  Copy the atom name
   */
  len = 4;
  char aname[5];
  for (i=0;i<len;i++){
    aname[i] = line[pos+i];
  }
  pos += len;
  aname[len] = '\0';
  strcpy( ad->atomName, aname );

  /*
   *  Copy the alternative location indicator
   */
  ad->altLoc = line[pos];
  pos ++;


  /*
   *  Copy the residual name
   */
  len = 3;
  char resName[4];
  for (i=0;i<len;i++){
    resName[i] = line[pos+i];
  }
  pos += len;
  resName[len] = '\0';
  strcpy( ad->resName, resName );

  /*
   *  21st position is blank
   */
  pos ++;
  
  /*
   *  Copy the chain ID
   */
  ad->chainID = line[pos];
  pos ++;
  

  /*
   *  Copy the residue sequence number
   */
  len = 4;
  char resSeq[5];
  for (i=0;i<len;i++){
    resSeq[i] = line[pos+i];
  }
  pos += len;
  resSeq[len] = '\0';
  ad->resSeq = atoi( resSeq );


  /*
   *  Code for the insertion of residues
   */
  ad->iCode = line[pos];
  pos ++;

  /*
   *  28-30 is blank
   */
  pos += 3;



  /*
   *  the x co-ordinate of the atom
   */
  len = 8;
  char x[9];
  for (i=0;i<len;i++){
    x[i] = line[pos+i];
  }
  pos += len;
  x[len] = '\0';
  ad->x = atof( x );


  /*
   *  the y co-ordinate of the atom
   */
  len = 8;
  char y[8];
  for (i=0;i<len;i++){
    y[i] = line[pos+i];
  }
  pos += len;
  ad->y = atof( y );

  
  /*
   *  the z co-ordinate of the atom
   */
  len = 8;
  char z[8];
  for (i=0;i<len;i++){
    z[i] = line[pos+i];
  }
  pos += len;
  ad->z = atof( z );


  /*
   *  occupancy
   */
  len = 6;
  char occ[6];
  for (i=0;i<len;i++){
    occ[i] = line[pos+i];
  }
  pos += len;
  ad->occupancy = atof( occ );


  /*
   *  temperature factor
   */
  len = 6;
  char tempfact[6];
  for (i=0;i<len;i++){
    tempfact[i] = line[pos+i];
  }
  pos += len;
  ad->tempFactor = atof( tempfact );


  /*
   *  67-76 is blank
   */
  pos += 10;


  /*
   *  element
   */
  len = 2;
  char element[3];
  for (i=0;i<len;i++){
    element[i] = line[pos+i];
  }
  pos += len;
  element[len] = '\0';
  strcpy( ad->element, element );



  /*
   *  charge
   */
  len = 2;
  char charge[3];
  for (i=0;i<len;i++){
    charge[i] = line[pos+i];
  }
  pos += len;
  charge[len] = '\0';
  strcpy( ad->charge, charge );


}


/*
 * print all the information contained in pdfdata.adata
 */
void print_pdb_data( pdbdata * pd )
{

  int i;

  for (i=0;i<pd->natoms;i++){
    
    printf("just resName %s\n", pd->adata[i].resName);
    
    printf("atom %d : %s | %d | %s | %c | %s | %c | %d | %c | %f | %f | %f | %f | %f | %s | %s\n",
	   i, pd->adata[i].rname, pd->adata[i].serial, pd->adata[i].atomName, 
	   pd->adata[i].altLoc, pd->adata[i].resName, pd->adata[i].chainID, 
	   pd->adata[i].resSeq, pd->adata[i].iCode, 
	   pd->adata[i].x, pd->adata[i].y, pd->adata[i].z, 
	   pd->adata[i].occupancy, pd->adata[i].tempFactor, 
	   pd->adata[i].element, pd->adata[i].charge);

  }

}


void pdb_chemical_composition( pdbdata * pd ){

  int i;
  int n[6];
  for(i=0;i<6;i++){
    n[i] = 0;
  }

  for (i=0;i<pd->natoms;i++){
    //printf("%s\n", pd->adata[i].element);
    if(pd->adata[i].element[1] == 'H'){
      n[0]++;
    }
    if(pd->adata[i].element[1] == 'C'){
      n[1]++;
    }

    if(pd->adata[i].element[1] == 'N'){
      n[2]++;
    }
    if(pd->adata[i].element[1] == 'O'){
      n[3]++;
    }
    if(pd->adata[i].element[1] == 'P'){
      n[4]++;
    }
    if(pd->adata[i].element[1] == 'S'){
      n[5]++;
    }
  }

  for(i=0;i<6;i++){
    printf("%d ", n[i]);
  }
  printf("\n");
}


void pdb_swap(atomdata* xp, atomdata* yp)
{
	atomdata temp = *xp;
	*xp = *yp;
	*yp = temp;
}


int pdb_greater(char a[3], char b[3])
{
	if (a[0] > b[0]) return 1;
	else if (a[0] < b[0]) return 0;
	else if (a[1] > b[1]) return 1; // if we are here, than the first symbols are the same
	else return 0;
}


void pdb_bubbleSort(pdbdata* pd)
{
	int n = pd->natoms;
	int i, j;
	for (i = 0; i < n - 1; i++)
		// Last i elements are already in place    
		for (j = 0; j < n - i - 1; j++)
			if (pdb_greater(pd->adata[j].element, pd->adata[j+1].element))
				pdb_swap(pd->adata + j, pd->adata + j + 1);
}


void pdb_bubbleSort1(pdbdata* pd, int* ia)
{
	int n = pd->natoms;
	int i, j, iatemp;
	for (i = 0; i < n - 1; i++)
		// Last i elements are already in place    
		for (j = 0; j < n - i - 1; j++)
			if (ia[j] < ia[j + 1])
			{
				pdb_swap(pd->adata + j, pd->adata + j + 1);
				iatemp = ia[j]; ia[j] = ia[j + 1]; ia[j + 1] = iatemp;
			}
}


int pdb_bubbleSort2(pdbdata* pd, int iFirst, int iLast)
{
	int n = iLast - iFirst;
	if (n <= 0 || iFirst < 0 || iFirst > pd->natoms || iLast < 0 || iLast > pd->natoms)
	{
		printf("\n!!!Bad input parameters in pdb_bubbleSort2!!!"); return -1;
	}
	int i, j;
	for (i = iFirst; i < iLast - 1; i++)
		// Last i elements are already in place    
		for (j = iFirst; j < iLast - (i - iFirst) - 1; j++)
			if (pd->adata[j].z > pd->adata[j + 1].z)
				pdb_swap(pd->adata + j, pd->adata + j + 1);
	return 0;
}


int pdb_bubbleSort3(pdbdata* pd, int iFirst, int iLast)
{
	int n = iLast - iFirst;
	if (n <= 0 || iFirst < 0 || iFirst > pd->natoms || iLast < 0 || iLast > pd->natoms)
	{
		printf("\n!!!Bad input parameters in pdb_bubbleSort3!!!"); return -1;
	}
	int i, j;
	for (i = iFirst; i < iLast - 1; i++)
		// Last i elements are already in place    
		for (j = iFirst; j < iLast - (i - iFirst) - 1; j++)
			if (pd->adata[j].occupancy < pd->adata[j + 1].occupancy)
				pdb_swap(pd->adata + j, pd->adata + j + 1);
	return 0;
}


int pdb_bubbleSort3a(pdbdata* pd, int iFirst, int iLast)
{
	int n = iLast - iFirst;
	if (n <= 0 || iFirst < 0 || iFirst > pd->natoms || iLast < 0 || iLast > pd->natoms)
	{
		printf("\n!!!Bad input parameters in pdb_bubbleSort3a!!!"); return -1;
	}
	int i, j;
	for (i = iFirst; i < iLast - 1; i++)
		// Last i elements are already in place    
		for (j = iFirst; j < iLast - (i - iFirst) - 1; j++)
			if (pd->adata[j].occupancy > pd->adata[j + 1].occupancy)
				pdb_swap(pd->adata + j, pd->adata + j + 1);
	return 0;
}



int pdb_atomnumbers(pdbdata* pd, int ia[])
{
	int i, j;
	for (i = 0; i < pd->natoms; i++)
	{
		// translate the element symbol into atomic number
		ia[i] = 0; j = 0;
		while (isspace(pd->adata[i].element[j])) j++; // skip leading white spaces
		if (j > 2)
		{
			printf("\n!!!Too many white spaces %s!!!", pd->adata[i].element);
			return -1;
		}
		if (pd->adata[i].element[j] == 'H') ia[i] = 1; // H = hydrogen
		else if (pd->adata[i].element[j] == 'C')
		{
			ia[i] = 6; // C = carbon
			if (toupper(pd->adata[i].element[j + 1]) == 'A') ia[i] = 20; // CA = calcium
			else if (toupper(pd->adata[i].element[j + 1]) == 'O') ia[i] = 27; // CO = cobalt
		}
		else if (pd->adata[i].element[j] == 'N')
		{
			ia[i] = 7; // N = nitrogen
			if (toupper(pd->adata[i].element[j + 1]) == 'A') ia[i] = 11; // NA = sodium
			else if (toupper(pd->adata[i].element[j + 1]) == 'I') ia[i] = 28; // NI = nickel
		}
		else if (pd->adata[i].element[j] == 'O') ia[i] = 8; // O = oxygen
		else if (pd->adata[i].element[j] == 'P')
		{
			ia[i] = 15; // P = posphorus
			if (toupper(pd->adata[i].element[j + 1]) == 'T') ia[i] = 78; // PT = platinum
			if (toupper(pd->adata[i].element[j + 1]) == 'B') ia[i] = 82; // PB = led
		}
		else if (pd->adata[i].element[j] == 'S') ia[i] = 16; // S = sulphur
		else if (pd->adata[i].element[j] == 'F') ia[i] = 26; // assuming FE = ferrum
		else if (pd->adata[i].element[j] == 'M')
		{
			if (toupper(pd->adata[i].element[j + 1]) == 'G') ia[i] = 12; // MG = magnesium
			else if (toupper(pd->adata[i].element[j + 1]) == 'N') ia[i] = 25; // MN = manganese
			else { printf("\n!!!Unknown atom type %s!!!", pd->adata[i].element); return -1; }
		}
		else if (pd->adata[i].element[j] == 'Z') ia[i] = 30; // assuming ZN = zinc
		else if (pd->adata[i].element[j] == 'A')
		{
			if (toupper(pd->adata[i].element[j + 1]) == 'U') ia[i] = 79; // AU = gold
			else if (toupper(pd->adata[i].element[j + 1]) == 'G') ia[i] = 47; // AG = silver
			else if (toupper(pd->adata[i].element[j + 1]) == 'L') ia[i] = 13; // AL = aluminum
			else { printf("\n!!!Unknown atom type %s!!!", pd->adata[i].element); return -1; }
		}
		else if (ia[i] == 0)
		{
			printf("\n!!!Unknown atom type %s!!!", pd->adata[i].element);
			return -1;
		}
	}
	return 0;
}


int pdb_symbols(pdbdata* pd, int ia[])
{
	int i;
	for (i = 0; i < pd->natoms; i++)
	{
		// translate atomic number into the element symbol
		if (ia[i] == 1) { pd->adata[i].element[0] = 'H'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 6) { pd->adata[i].element[0] = 'C'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 7) { pd->adata[i].element[0] = 'N'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 8) { pd->adata[i].element[0] = 'O'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 11) { pd->adata[i].element[0] = 'N'; pd->adata[i].element[1] = 'A'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 12) { pd->adata[i].element[0] = 'M'; pd->adata[i].element[1] = 'G'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 13) { pd->adata[i].element[0] = 'A'; pd->adata[i].element[1] = 'L'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 15) { pd->adata[i].element[0] = 'P'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 16) { pd->adata[i].element[0] = 'S'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 20) { pd->adata[i].element[0] = 'C'; pd->adata[i].element[1] = 'A'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 25) { pd->adata[i].element[0] = 'M'; pd->adata[i].element[1] = 'N'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 26) { pd->adata[i].element[0] = 'F'; pd->adata[i].element[1] = 'E'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 27) { pd->adata[i].element[0] = 'C'; pd->adata[i].element[1] = 'O'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 28) { pd->adata[i].element[0] = 'N'; pd->adata[i].element[1] = 'I'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 30) { pd->adata[i].element[0] = 'Z'; pd->adata[i].element[1] = '\0'; }
		else if (ia[i] == 47) { pd->adata[i].element[0] = 'A'; pd->adata[i].element[1] = 'G'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 78) { pd->adata[i].element[0] = 'P'; pd->adata[i].element[1] = 'T'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 79) { pd->adata[i].element[0] = 'A'; pd->adata[i].element[1] = 'U'; pd->adata[i].element[2] = '\0'; }
		else if (ia[i] == 82) { pd->adata[i].element[0] = 'P'; pd->adata[i].element[1] = 'B'; pd->adata[i].element[2] = '\0'; }
		else
		{
			printf("\n!!!Unknown atom number %d!!!", ia[i]);
			return -1;
		}
	}
	return 0;
}
