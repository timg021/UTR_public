#ifndef ADDICE_H   // only include this file if its not already
#define ADDICE_H   // remember that this has been included

#include <string>
#include <vector>
#include "pdb.h"

int AddIce(float iceThick, float ctblength, int natom, int** pZnum, float** px, float** py, float** pz, float** pocc, float** pwobble, float wobbleaver, unsigned long* piseed);
int AddCarbon(double cThick, double cWidth, pdbdata& pd, unsigned long* piseed, double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax);

#endif //ADDICE_H
