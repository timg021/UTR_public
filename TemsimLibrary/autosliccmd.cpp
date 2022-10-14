/*
      *** autosliceCMD.cpp ***

Front end for autoslic.cpp containing a command line user interface.

------------------------------------------------------------------------
Copyright 1998-2012 Earl J. Kirkland

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 
------------------------------------------------------------------------

  ANSI C and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  gcc -O -fopenmp -o autoslic autoslic.c slicelib.o
                       tiffsubs.o  cfpix.o -lfftw3f

  Transmit an electron wave through a specimen using the
  multislce method with automatic slicing.  Read in the (x,y,z)
  coordinates of the whole specimen and break into slices
  on-the-fly.

  started 24-july-1996 E. Kirkland
  working 19feb-1997 ejk
  last revised 19-feb-1997 ejk
  added look-up-table vzatomLUT() for 3X-4X increase 
        in speed 23-may-1997 ejk
  put bandwith limit inside trlayer() 1-oct-1997 ejk
  added Gaussian thermal displacements 1-oct-1997 ejk
  removed /sqrt(3) from Thermal rms displacements 
    to be consistent with Int'l X-ray tables 22-dec-1997 ejk
  corrected zmin/max error with thermal displac. 24-dec-1997 ejk
  fixed small aliasing problem 5-jan-1998 ejk
  added unit cell replication option and moved ReadXYZcoord()
    into slicelib.c  11-jan-1998 ejk
  added astigmatism and modify to use different set of
    random offsets on each illum. angle with partial coherence
         5-feb-1998 ejk
  fix typo in z range message with partial coherence and
    thermal vibrations 9-july-1998 ejk
  update memory allocation routines 19-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  fixed bug in zmin/zmax calculation in coherent mode
     (move to after sortByZ() - it was before ) 8-jan-2002 ejk
  add cross section option (in non-partial coherence mode only)
        27-may-2005 ejk
  convet to faster sortByZ() 8-feb-2006 ejk
  move sortbyz() to slicelib.c 5-sep-2006 ejk
  add echo on y position in pixels for xz mode 4-may-2007 ejk
  update data type of nxl,nyl to be consistent with new tiffsubs
     17-jul-2007 ejk
  move xz depthpix save to be after transmit+propagate to get a
     full slice and proper anti-aliasing and also be consisten
     with what you get doing it by hand  and increase possible
     slices output (nz was off by one) 24-jan-2008 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better (consistent with autostem) 23-mar-2008 ejk
  take small things out of loop in trlayer() 14-may-2008 ejk
  parameterize vzatomLUT() vs r^2 instead of r to avoid a lot of sqrt()
      calls (a little faster)  6-jun-2008 ejk
  move vzatomLUT() to slicelib.c  11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  add Cs5 (and Cs->Cs3) 15-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4
     and convert to 4 char TAB size formatting 21-feb-2010 ejk
  add parallel computing of a few parts 21-feb-2010 ejk
  start conversion to faster FFTW 24-feb-2010 ejk
  move some things into slicelibW.c to share 6-mar-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add option to average over many frozen phonon
      configurations 3-aug-2010 ejk
  add multipole aberrations to probe 12-may-2011 ejk
  start conversion to floatTIFF.cpp and C++ 28-may-2012 ejk
  working 3-jun-2012 ejk
  convert to cfpix/fftw class from raw fftw 13-nov-2012 to 21-nov-2012 ejk
  move calculation into a class with separate command line front end
      29-may-2013 ejk
  change RNG seed argument to referenece so it get updated for 
      successive calls 21-sep-2013 ejk
  fix %ld to %d format in printout of aslice.nillum  12-oct-2013 ejk
  add param[] for mode 15-oct-2013 ejk
  convert to streams and strings 22-mar-2014 ejk
  fix minor formatting issue 2-jun-2014 ejk

  ax,by,cz  = unit cell size in x,y
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 4 characters 
  
*/

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <istream>
#include <math.h>
#include <cmath>
#include <ctime>
#include <chrono>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output

#include "XArray2D.h"
#include "XA_data.h"
#include "XA_fft2.h"
#include "XA_spln2.h"
#include "XA_nrrand.h"
#include "XA_tiff.h"
#include "XA_move2.h"
#include "XA_data.h"
#include "AddIce.h"

//!!!! temporary (to check the effect of WPOA in free-space propagation)
//#include "XA_iwfr.h"

using namespace std;

#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // misc. routines for multislice
#include "floatTIFF.hpp"   // file I/O routines in TIFF format
#include "autoslic.hpp"    //  the calculation engine

#include "autosliccmd.h" //declaration of the function autosliccmd()

#define MANY_ABERR      //  define to include many aberrations

#ifdef USE_OPENMP
//#include <omp.h>
/*  get wall time for benchmarking openMP */
//#define walltim() ( omp_get_wtime() )
//double walltimer;
#endif

int autosliccmd(vector<string> params, vector<xar::Pair> defocus, xar::Pair astigm, xar::Pair xyshifts, vector<string> fileout)
{
	Counter_Obj thread_counter; // increments the thread counter on construction and decrements it on destruction
	try
	{
		printf("\nNumber of active threads (inside worker thread) = %d", thread_counter.GetCount());
		thread_counter.SetUpdated(true); // lets the main thread know that the thread counter has been updated

		string filein, filestart, filebeam, filecross, cline, description;

		const char version[] = "2-jun-2014 (ejk)";

		int lstart = 0, lpartl = 0, lbeams = 0, lwobble = 0, lcross = 0, nwobble = 1;
		int ix, iy, iz, nx, ny, nzout, i, nslic0, islice,
			ndf, nbout(0), ib, ncellx, ncelly, ncellz, NPARAM;
		int nillum, nzbeams;
		int* hbeam(0), * kbeam(0);
		int natom, * Znum, done, status, multiMode(0);
		long  ltime;

		unsigned long iseed;
		//@@@@@ start TEG code
		long iseed2;
		//@@@@@ end TEG code

		float v0, mm0, wavlen, rx, ry, ax, by, cz,
			rmin, rmax, aimin, aimax, ctiltx(0), ctilty(0),
			acmin, acmax, Cs3, Cs5, df0, sigmaf, dfdelt(0), aobj(0),
			temperature, ycross(0), dx, dy;

		float wmin, wmax, xmin, xmax, ymin, ymax, zmin, zmax;

		float* x, * y, * z, * occ, * wobble;
		float* param, * sparam;
		vector<float> vxcoor, vycoor, vzcoor, vocc, vwobble;
		vector<int> vZnum;

		double timer, deltaz, vz;
		double pi, pi180;

		//@@@@@ start TEG code
		int nmode; // the switch between multislice(0), projection(1) or 1st Born(2) approximations
		int noutput; // the switch between intensity(0), phase(1), complex amplitude(2) or 3D potential(3) output form of the result
		int nfftwinit; // the switch between copying(0) or initializing from new (1) the FFTW plan in autoslic
		int nIceSave; // the save XYZ files with ice switch: not_save(0) or save(1)
		float ctblength(0); // x and y side length of the slab (containing the sample) in Angstroms
		float ctblengthz(0); // z thickness of the slab (containing the sample) in Angstroms
		float iceThick(1); // thickness of the optional layer of amorphous ice in Angstroms
		double R1(0); // distance from the point source (0 = plane wave)
		double angleZ(0); // sample rotation angle in radians (around z axis)
		double angleY(0); // sample rotation angle in radians (around y' axis)
		double angleZ2(0); // sample rotation angle in radians (around z" axis, which is supposed to be the illumination axis)
		double dxc, dyc, dzc; // x, y and z shifts of the centre of rotation from the default position at the centre of the CT cube
		//@@@@@ end TEG code

		ofstream fp1;

		cfpix pix;		// to get results of calculation
		cfpix wave0;	// initial wavefunction (if used)
		cfpix depthpix;	// to get xz cross section results
		cfpix beams;	// to get valuse of requested beams during propagation

		floatTIFF myFile;	//  file input/output
		autoslic aslice;	// has the calculation engine

	/*  echo version date and get input file name */
#ifdef _DEBUG
		cout << "autoslic(e) version dated " << version << endl;
		cout << "Copyright (C) 1998-2013 Earl J. Kirkland" << endl;
		cout << "This program is provided AS-IS with ABSOLUTELY NO WARRANTY" << endl;
		cout << " under the GNU general public license" << endl << endl;

		cout << "perform CTEM multislice with automatic slicing and FFTW" << endl;
#ifdef USE_OPENMP
		cout << "and multithreaded using openMP" << endl;
#endif
		cout << " " << endl;
#endif


		pi = 4.0 * atan(1.0);
		pi180 = pi / 180.0;
		NPARAM = myFile.maxParam();
		param = (float*)malloc1D(NPARAM, sizeof(float), "param");
		sparam = (float*)malloc1D(NPARAM, sizeof(float), "sparam");
		for (ix = 0; ix < NPARAM; ix++) param[ix] = 0.0F;

		//@@@@@ start TEG code
		//FILE* ff0 = fopen("autoslic.txt", "rt");
		char chaa[1024], cinarg[1024];
		//if (!ff0)
		//{
		//	cout << "!!!Cannot open input parameter file!!! Input any character to exit...";
		//	cin >> chaa;
		//	exit(-1);
		//}
		//cout << "Name of file with input atomic "
		//       << "coord. in x,y,z format:"  << endl;
		//cin >> filein;
		//filein = "1grl.xyz";
		if (sscanf(params[0].data(), "%s %s", chaa, cinarg) != 2)
			throw std::runtime_error("Error reading line 1 of input parameter array.");
		filein = cinarg;

		/*  get simulation options */

			//cout << "Replicate unit cell by NCELLX,NCELLY,NCELLZ :" << endl;
			//cin >> ncellx >> ncelly >> ncellz;
			//if( ncellx < 1 ) ncellx = 1;
			//if( ncelly < 1 ) ncelly = 1;
			//if( ncellz < 1 ) ncellz = 1;
			//ncellx = 1; ncelly = 1; ncellz = 1;
		if (sscanf(params[1].data(), "%s %d %d %d", chaa, &ncellx, &ncelly, &ncellz) != 4)
			throw std::runtime_error("Error reading line 2 of input parameter array.");

		//cout << "Name of file to get binary output of multislice result:" << endl;
		//cin.getline( fileout, NCMAX );   //  ????
		//cin >> fileout ;
		//fileout = "1grl.tif";
		//if (sscanf(params[2].data(), "%s %s", chaa, cinarg) != 2)
		//	throw std::runtime_error("Error reading line 3 of input parameter array.");
		//fileout = cinarg;

		//lpartl = askYN("Do you want to include partial coherence");
		//lpartl = 0;
		if (sscanf(params[3].data(), "%s %d", chaa, &lpartl) != 2)
			throw std::runtime_error("Error reading line 4 of input parameter array.");
		if (sscanf(params[4].data(), "%s %g %g", chaa, &acmin, &acmax) != 3)
			throw std::runtime_error("Error reading line 5 of input parameter array.");
		if (sscanf(params[5].data(), "%s %g %g", chaa, &Cs3, &Cs5) != 3)
			throw std::runtime_error("Error reading line 6 of input parameter array.");
		if (sscanf(params[6].data(), "%s %g %g %g", chaa, &df0, &sigmaf, &dfdelt) != 4)
			throw std::runtime_error("Error reading line 7 of input parameter array.");
		if (sscanf(params[7].data(), "%s %g", chaa, &aobj) != 2)
			throw std::runtime_error("Error reading line 8 of input parameter array.");
		aobj = (float)fabs(aobj * 0.001F);
		if (sscanf(params[8].data(), "%s %d", chaa, &lstart) != 2)
			throw std::runtime_error("Error reading line 9 of input parameter array.");
		if (sscanf(params[9].data(), "%s %s", chaa, cinarg) != 2)
			throw std::runtime_error("Error reading line 10 of input parameter array.");
		filestart = cinarg;
		if (sscanf(params[10].data(), "%s %g", chaa, &v0) != 2)
			throw std::runtime_error("Error reading line 11 of input parameter array.");
		if (sscanf(params[11].data(), "%s %d %d", chaa, &nx, &ny) != 3)
			throw std::runtime_error("Error reading line 12 of input parameter array.");
		//if (sscanf(params[12].data(), "%s %g %g", chaa, &ctiltx, &ctilty) != 3)
		//	throw std::runtime_error("Error reading line 13 of input parameter array.");
		//@@@@@@@ TEG has reused params[12] to pass different data into this program. 
		// ctiltx and ctilty have been set to zero
		if (sscanf(params[12].data(), "%s %lg %lg %lg", chaa, &dxc, &dyc, &dzc) != 4)
			throw std::runtime_error("Error reading centre of rotation x, y and z shifts from input parameter array.");
		if (sscanf(params[13].data(), "%s %lg", chaa, &deltaz) != 2)
			throw std::runtime_error("Error reading line 14 of input parameter array.");
		if (sscanf(params[14].data(), "%s %d", chaa, &lbeams) != 2)
			throw std::runtime_error("Error reading line 15 of input parameter array.");
		if (sscanf(params[15].data(), "%s %s", chaa, cinarg) != 2)
			throw std::runtime_error("Error reading line 16 of input parameter array.");
		filebeam = cinarg;
		if (sscanf(params[16].data(), "%s %d", chaa, &nbout) != 2)
			throw std::runtime_error("Error reading line 17 of input parameter array.");
		if (sscanf(params[17].data(), "%s %d", chaa, &lwobble) != 2)
			throw std::runtime_error("Error reading line 18 of input parameter array.");
		if (lwobble != 0) lwobble = 1;
		if (sscanf(params[18].data(), "%s %g", chaa, &temperature) != 2)
			throw std::runtime_error("Error reading line 19 of input parameter array.");
		if (sscanf(params[19].data(), "%s %d", chaa, &nwobble) != 2)
			throw std::runtime_error("Error reading line 20 of input parameter array.");
		//!!! TEG NOTE: be careful about this parameter: even though this module can be called to calculate a single configuration,
		//it needs to know whether multiple configurations are being calculated by the controlling module (MultisliceK.cpp).
		//if lwobble == 1 and nwobble == 1, the code below will used the exponential model of the DW factor, instead of
		//calculating multiple configurations
		//if (sscanf(params[20].data(), "%s %d", chaa, &iseed1) != 2) // this is not used any more
		//	throw std::runtime_error("Error reading line 21 of input parameter array.");
		if (sscanf(params[21].data(), "%s %d", chaa, &lcross) != 2)
			throw std::runtime_error("Error reading line 22 of input parameter array.");
		if (sscanf(params[22].data(), "%s %s", chaa, cinarg) != 2)
			throw std::runtime_error("Error reading line 23 of input parameter array.");
		filecross = cinarg;
		if (sscanf(params[23].data(), "%s %g", chaa, &ycross) != 2)
			throw std::runtime_error("Error reading line 24 of input parameter array.");

		if (sscanf(params[24].data(), "%s %lf %lf", chaa, &angleZ, &angleY) != 3)
			throw std::runtime_error("Error reading line 25 of input parameter array.");

		if (sscanf(params[25].data(), "%s %d", chaa, &nmode) != 2)
			throw std::runtime_error("Error reading line 26 of input parameter array.");
		//if (sscanf(params[26].data(), "%s %g %g %g", chaa, &defocus_min, &defocus_max, &defocus_step) != 4)
		//	throw std::runtime_error("Error reading line 27 of input parameter array.");
		if (sscanf(params[27].data(), "%s %d", chaa, &noutput) != 2)
			throw std::runtime_error("Error reading line 28 of input parameter array.");
		bool boutTIFF;
		if (noutput == 0 || noutput == 1 || noutput == 3)
		{
			if (xar::GetFileExtension(fileout[0]) == string(".TIFF") || xar::GetFileExtension(fileout[0]) == string(".TIF")) 
				boutTIFF = true;
			else 
				if (xar::GetFileExtension(fileout[0]) == string(".GRD")) boutTIFF = false; else throw std::runtime_error("Error: incorrect filename extension.");
		}
		else
		{
			if (noutput == 2)
			{
				if (xar::GetFileExtension(fileout[0]) != string(".GRC")) throw std::runtime_error("Error: incorrect filename extension.");
			}
			else throw std::runtime_error("Error: incorrect output format.");
		}

		if (sscanf(params[28].data(), "%s %d", chaa, &nfftwinit) != 2)
			throw std::runtime_error("Error reading line 29 of input parameter array.");

		if (sscanf(params[29].data(), "%s %g", chaa, &iceThick) != 2)
			throw std::runtime_error("Error reading line 30 of input parameter array.");

		if (sscanf(params[30].data(), "%s %d", chaa, &nIceSave) != 2)
			throw std::runtime_error("Error reading line 31 of input parameter array.");

		if (sscanf(params[31].data(), "%s %lg", chaa, &R1) != 2)
			throw std::runtime_error("Error reading line 32 of input parameter array.");

		//fclose(ff0);
		//cout << "Input parameter file has been read successfully!\n";

		//acmin = acmax = 0;
		if (lpartl == 1) {
			//cout << "Illumination angle min, max in mrad:" << endl;
			//cin >> acmin >> acmax;
			acmin = acmin * 0.001F;
			acmax = acmax * 0.001F;
			//cout << "Spherical aberration Cs3, Cs5(in mm.):" << endl;
			//cin >>  Cs3 >> Cs5;
			param[pCS] = (float)(Cs3 * 1.0e7);
			param[pCS5] = (float)(Cs5 * 1.0e7);
			//cout << "Defocus, mean, standard deviation, and"
			//       " sampling size (in Angstroms) =" << endl;
			//cin >> df0 >> sigmaf >> dfdelt;
			param[pDEFOCUS] = (float)df0;
			param[pDDF] = (float)sigmaf;
			//cout << "Objective aperture (in mrad) =" << endl;
			//cin >> aobj;
			//@@@@@ end TEG code
			//aobj = (float)fabs(aobj * 0.001F);
#ifdef MANY_ABERR 
			/*   get higher order aberrations if necessary */
			cout << "type higher order aber. name (as C32a, etc.) followed\n"
				<< " by a value in mm. (END to end)" << endl;
			done = multiMode = 0;
			do {
				cin >> cline;
				if (cline.compare("END") == 0) {
					done = 1;
				}
				else {
					cin >> vz;
					status = readCnm(cline, param, vz);
					if (status < 0) {
						cout << "unrecognized aberration, exit..." << endl;
						exit(EXIT_SUCCESS);
					}
					else multiMode = 1;
				}
			} while (!done);

#endif
			lstart = 0;
		}
		else {
#ifdef _DEBUG
			cout << "NOTE, the program image must also be run." << endl;
#endif
			//@@@@@ start TEG code
			//lstart = askYN("Do you want to start from previous result");
			//lstart = 0;
			//@@@@@ end TEG code
		}

		//@@@@@ start TEG code
		if (lwobble == 1 && nwobble > 1 && !(noutput == 0 || noutput == 3))
			throw std::runtime_error("Only intensity or 3D potential output is allowed in the presence of thermal vibrations with multiple configurations.");
		if (lstart == 1) {
			//cout << "Name of file to start from:" << endl;
			//cin >> filestart;
		}
		else {
			//cout << "Incident beam energy in kev:" << endl;
			//cin >> v0;
			//v0 = 200;
			//cout << "Wavefunction size in pixels, Nx,Ny:" << endl;
			//cin >> nx >> ny;
			//nx = ny = 1024;
			//@@@@@ end TEG code
		}

		//@@@@@ start TEG code
		//cout << "Crystal tilt x,y in mrad.:" << endl;
		//cin >> ctiltx >> ctilty;
		//ctiltx = ctilty = 0;
		//@@@@@ end TEG code
		ctiltx = ctiltx / 1000;
		ctilty = ctilty / 1000;

		/*  remember that the slice thickness must be > atom size
			to use projected atomic potential */
		//@@@@@ start TEG code
		//cout << "Slice thickness (in Angstroms):" << endl;
		//cin >> deltaz;
		//deltaz = 1.3575;

		//if (deltaz < 1.0) {
		//	cout << "\nWARNING: this slice thickness is probably too thin"
		//		<< " for autoslice to work properly." << endl;
		//}
		//@@@@@ end TEG code
		if (lpartl == 0) {
			//@@@@@ start TEG code
			//lbeams = askYN("Do you want to record the (real,imag) value\n"
			//" of selected beams vs. thickness");
			//lbeams = 0;
			//@@@@@ end TEG code
			if (lbeams == 1) {
				//cout << "Name of file for beams info:" << endl;
				//cin >> filebeam;
				//cout << "Number of beams:" << endl;
				//cin >> nbout;
				if (nbout < 1) nbout = 1;
				hbeam = (int*)malloc1D(nbout, sizeof(int), "hbeam");
				kbeam = (int*)malloc1D(nbout, sizeof(int), "kbeam");
				for (ib = 0; ib < nbout; ib++) {
					cout << "Beam " << ib + 1 << ", h,k=" << endl;
					cin >> hbeam[ib] >> kbeam[ib];
				}
			}
		}

		//@@@@@ start TEG code
		//lwobble = askYN("Do you want to include thermal vibrations");
		//lwobble = 0;
		//@@@@@ end TEG code
		if (lwobble == 1 || iceThick > 0) {
			//cout << "Type the temperature in degrees K:" << endl;
			//cin >> temperature ;
			//cout << "Type number of configurations to average over:" << endl;
			//cin >>  nwobble;
			if (nwobble < 1) nwobble = 1;
			/* get random number seed from time if available
				otherwise ask for a seed */
			//@@@ TEG: note that the accuracy of the simple time() function is not enough here:
			//nanoseconds need to be used in order to create different seeds for different threads
			std::chrono::system_clock::time_point current_time = std::chrono::system_clock::now();
			std::chrono::system_clock::duration dtn = current_time.time_since_epoch();
			auto dtn_nano = std::chrono::duration_cast<std::chrono::nanoseconds>(dtn);
			ltime = (long)dtn_nano.count();
			//ltime = (long)time(NULL);
			iseed = std::abs(ltime);
			//printf("\n@@@@@ iseed = %ld", iseed);
			iseed2 = -std::abs(ltime); // gasdev requires a negative integer seed for initialization, which should not be changed between successive calls
		}
		else temperature = 0.0F;

		if (lpartl == 0) {
			//@@@@@ start TEG code
			//lcross = askYN("Do you want to output intensity vs. depth cross section");
			//lcross = 0;
			//@@@@@ end TEG code
			if (lcross == 1) {
				//cout << "Type name of file to get depth profile image:" << endl;
				//cin >> filecross;
				//cout << "Type y position of depth cross section (in Ang.):" << endl;
				//cin >> ycross;
			}
		}

		/* start timing the actual computation just for fun */

		timer = cputim();
#ifdef USE_OPENMP
		//walltimer = walltim();  /* wall time for openMP */
#endif

/*  get starting value of transmitted wavefunction if required
   (this can only be used in coherent mode)
	remember to save params for final output pix  */

		if (lstart == 1) {
			if (myFile.read(filestart.c_str()) != 1) {
				cout << "Cannot open input file: " << filestart << endl;
				exit(0);
			}

			if (myFile.getnpix() != 2) {
				cout << "Input file " << filestart << " must be complex, can't continue." << endl;
				exit(0);
			}

			nx = myFile.nx();
			ny = myFile.ny();

			nx = nx / 2;
			wave0.resize(nx, ny);

			//  save starting pix for later
			for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++) {
				wave0.re(ix, iy) = myFile(ix, iy);
				wave0.im(ix, iy) = myFile(ix + nx, iy);
			}

			//  save parameters to verify successive images are same size etc.
			for (i = 0; i < NPARAM; i++) sparam[i] = myFile.getParam(i);

			ax = sparam[pDX] * nx;
			by = sparam[pDY] * ny;
			v0 = sparam[pENERGY];
			nslic0 = (int)sparam[pNSLICES];
#ifdef _DEBUG
			cout << "Starting pix range " << sparam[pRMIN] << " to " << sparam[pRMAX]
				<< " real\n" << "           " << sparam[pIMIN] << " to "
				<< sparam[pIMAX] << " imag" << endl;
			cout << "Beam voltage = " << v0 << " kV" << endl;
			cout << "Old crystal tilt x,y = " << 1000. * sparam[pXCTILT] << ", "
				<< 1000. * sparam[pYCTILT] << " mrad" << endl;
#endif
		}
		else
		{
			nslic0 = 0;     /* end if( lstart...) */
		}

 /*  calculate relativistic factor and electron wavelength */

		mm0 = 1.0F + v0 / 511.0F;
		wavlen = (float)wavelength(v0);
#ifdef _DEBUG
		cout << "electron wavelength = " << wavlen << " Angstroms" << endl;
#endif

		/*  read in specimen coordinates and scattering factors */

		natom = ReadXYZcoord(filein.c_str(), ncellx, ncelly, ncellz,
			&ax, &by, &cz, &Znum, &x, &y, &z, &occ, &wobble,
			description);

		if (!natom)
		{
			cout << "!!!\nZero atomic coordinates read. Input any key to exit...!!!";
			cin >> chaa;
			exit(-1);
		}
#ifdef _DEBUG	
		cout << natom << " atomic coordinates read in" << endl;
		cout << description << endl;

		cout << "Size in pixels Nx, Ny= " << nx << " x " << ny << " = " << nx * ny
			<< " beams" << endl;
		cout << "Lattice constant a,b = " << ax << ", " << by << endl;
#endif

		//@@@@@ start TEG code
		//!!! Note that the CT projection simulation cube side length will be set equal to the largest dimension of the unit cell in the input XYZ file.
		//Consequently, the largest "unit cell" dimension in the input XYZ file should be large enough, so that CT rotation won't take a part of the molecule outside the simulation cube.
		if (ax > by)
			if (ax > cz) ctblength = ax; else ctblength = cz;
		else
			if (by > cz) ctblength = by; else ctblength = cz;
		ctblengthz = ctblength; // ctblengthz may be increased later if an ice layer is added

		// define incident plane or spherical wavefront
		wave0.resize(nx, ny);
		if (R1 == 0) // plane incident wave case
			for (int ix = 0; ix < nx; ix++)
				for (int iy = 0; iy < ny; iy++)
				{
					wave0.re(ix, iy) = 1.0f;
					wave0.im(ix, iy) = 0.0f;
				}
		else // spherical incident wave case
		{
			double R1in = R1 - 0.5 * ctblength; // R1 is defined relative to the centre of the object, but the incident wave will be defined at the entrance surface
			double xx, yy, yyR12;
			double R12 = R1in * R1in;
			double aR12 = 1.0 / R12;
			double tPIdlambda = 2.0 * pi / wavelength(v0);
			double xst = ctblength / nx, yst = ctblength / ny;
			int nx2 = nx / 2, ny2 = ny / 2;
			double phamax = tPIdlambda * sqrt(xst * xst * nx2 * nx2 + yst * yst * ny2 * ny2 + R12);
			double ampmin = 1.0 / sqrt(1.0 + (xst * xst * nx2 * nx2 + yst * yst * ny2 * ny2) * aR12);
			xar::XArray2D<float> amp0(ny, nx), pha0(ny, nx);
			//!!! note that Kirkland's code seems to be mostly(!) using the array(ix,iy) order (which is transposed with respect to XArray convention)
			for (int j = 0; j < ny; j++)
			{
				yy = double(j - ny2) * yst;
				yy *= yy;
				yyR12 = yy + R12;
				for (int i = 0; i < nx; i++)
				{
					xx = double(i - nx2) * xst;
					xx *= xx;
					amp0[j][i] = float(1.0 / sqrt(1.0 + (xx + yy) * aR12));
					xx += yyR12;
					pha0[j][i] = float(tPIdlambda * sqrt(xx) - phamax);
				}
			}
			//!!!NOTE: edge effects of the spherical wave create very strong artefacts, but all attempts to smooth out the edges have been unsuccessful so far
			//index_t nx8 = index_t(0.125 * nx + 0.5); // width of the edge layer for smooth transition
			//xar::XArray2DMove<float> moveamp0(amp0);
			//moveamp0.MaskSmooth(nx8, 50.0, float(ampmin)); // smoothly transition the amplitude to ampmin at the edges
			//xar::XArray2DMove<float> movepha0(pha0);
			//movepha0.MaskSmooth(nx8, 50.0, 0.0f); // smoothly transition the phase to 0 at the edges (note that phamax has been subtracted globally)
			for (ix = 0; ix < nx; ix++)
				for (iy = 0; iy < ny; iy++)
				{
					wave0.re(ix, iy) = float(amp0[iy][ix] * cos(pha0[iy][ix]));
					wave0.im(ix, iy) = float(amp0[iy][ix] * sin(pha0[iy][ix]));
				}
		}

		// save the atomic coordinates for the future
		int natom0(-1);
		if (defocus.size() > 1)
		{
			natom0 = natom;
			vxcoor.resize(natom); vycoor.resize(natom); vzcoor.resize(natom); vocc.resize(natom); vZnum.resize(natom); vwobble.resize(natom);
			for (size_t k = 0; k < natom; k++)
			{
				vxcoor[k] = x[k];
				vycoor[k] = y[k];
				vzcoor[k] = z[k];
				vocc[k] = occ[k];
				vZnum[k] = Znum[k];
				vwobble[k] = wobble[k];
			}
		}

		// rotate the sample as necessary (using extrinsic(!) Euler angles in the opposite order to the intrinsic Euler angles)
		// and also optionally apply transverse (XY) shifts of the molecule. 
		// Note that unlike the rotations, XY shifts can move the moleculue outside the containing slab.
		float xc(float(ctblength / 2.0 + dxc)), yc(float(ctblength / 2.0 + dyc)), zc(float(ctblength / 2.0 + dzc)), xxx, yyy, zzz;
		float sinZ = float(sin(angleZ)), cosZ = float(cos(angleZ)), sinY = float(sin(angleY)), cosY = float(cos(angleY));

		// start of the cycle over defocus distances
		for (size_t jjj = 0; jjj < defocus.size(); jjj++)
		{
			printf("\n z'' rotation angle = %g (deg), defocus = %g (A)", defocus[jjj].a, defocus[jjj].b);

			// restore the original atomic coordinates
			if (jjj > 0)
			{
				natom = natom0;
				for (size_t k = 0; k < natom; k++)
				{
					x[k] = vxcoor[k];
					y[k] = vycoor[k];
					z[k] = vzcoor[k];
					occ[k] = vocc[k];
					Znum[k] = vZnum[k];
					wobble[k] = vwobble[k];
				}
			}

			// rotation around Z axis (by the Z rotation angle)
			for (size_t k = 0; k < natom; k++)
			{
				xxx = xc + (x[k] - xc) * cosZ + (-y[k] + yc) * sinZ;
				yyy = yc + (x[k] - xc) * sinZ + (y[k] - yc) * cosZ;
				x[k] = xxx; y[k] = yyy;
			}

			// rotation around Y axis (by the Y' rotation angle)
			for (size_t k = 0; k < natom; k++)
			{
				xxx = xc + (x[k] - xc) * cosY + (z[k] - zc) * sinY;
				zzz = zc + (-x[k] + xc) * sinY + (z[k] - zc) * cosY;
				x[k] = xxx; z[k] = zzz;
			}
			#if(0)
			// rotation around Z axis (by the Z" rotation angle)
			if (defocus[jjj].a != 0)
			{
				angleZ2 = defocus[jjj].a * pi180;
				float sinZ2 = float(sin(angleZ2)), cosZ2 = float(cos(angleZ2));

				for (size_t k = 0; k < natom; k++)
				{
					xxx = xc + (x[k] - xc) * cosZ2 + (-y[k] + yc) * sinZ2;
					yyy = yc + (x[k] - xc) * sinZ2 + (y[k] - yc) * cosZ2;
					x[k] = xxx; y[k] = yyy;
				}
			}

			// shift along X and Y
			if (xyshifts.a != 0 || xyshifts.b != 0)
			{
				float xshift = float(xyshifts.a);
				float yshift = float(xyshifts.b);
				for (size_t k = 0; k < natom; k++)
				{
					x[k] += xshift;
					y[k] += yshift;
				}
			}
			
			// find all atoms located outside ctbox
			bool bAtomsOutsideCtbox(false);
			for (size_t k = 0; k < natom; k++) 
			{
				if (x[k] < 0 || x[k] > ctblength || y[k] < 0 || y[k] > ctblength || z[k] < 0 || z[k] > ctblength)
				{
					occ[k] = 0;
					if (!bAtomsOutsideCtbox) // do this only once when the first atom outside ctbox has been found
					{
						bAtomsOutsideCtbox = true;
						printf("\n!!!WARNING: some atoms are outside ctbox and will be removed [anlge_Z = %g, angle_Y' = %g, angle_Z'' = %g (rad)]", angleZ, angleY, defocus[jjj].a * pi180);
					}
				}
			}

			// remove all atoms that are located outside ctbox
			if (bAtomsOutsideCtbox) 
			{
				int natom1(0);
				vector<float> vxtemp(natom), vytemp(natom), vztemp(natom), vocctemp(natom), vwobbletemp(natom);
				vector<int> vZnumtemp(natom);
				for (size_t k = 0; k < natom; k++)
				{
					if (occ[k] != 0)
					{
						vxtemp[natom1] = x[k];
						vytemp[natom1] = y[k];
						vztemp[natom1] = z[k];
						vocctemp[natom1] = occ[k];
						vZnumtemp[natom1] = Znum[k];
						vwobbletemp[natom1] = wobble[k];
						natom1++;
					}
				}
				natom = natom1;
				for (size_t k = 0; k < natom1; k++)
				{
					x[k] = vxtemp[k];
					y[k] = vytemp[k];
					z[k] = vztemp[k];
					occ[k] = vocctemp[k];
					Znum[k] = vZnumtemp[k];
					wobble[k] = vwobbletemp[k];
				}
			}
			#endif
			//@@@@@ end TEG code

			/*  calculate the total specimen volume and echo */
			xmin = xmax = x[0];
			ymin = ymax = y[0];
			zmin = zmax = z[0];
			wmin = wmax = wobble[0];
			float wobbleaver = 0.0f;

			for (i = 0; i < natom; i++) {
				if (x[i] < xmin) xmin = x[i];
				if (x[i] > xmax) xmax = x[i];
				if (y[i] < ymin) ymin = y[i];
				if (y[i] > ymax) ymax = y[i];
				if (z[i] < zmin) zmin = z[i];
				if (z[i] > zmax) zmax = z[i];
				if (wobble[i] < wmin) wmin = wobble[i];
				if (wobble[i] > wmax) wmax = wobble[i];
				wobbleaver += wobble[i];
			}
			wobbleaver /= natom; // we need this value in the case of added ice, where this value is used to wobble each atom
	#ifdef _DEBUG
			cout << "Total specimen range is\n"
				<< xmin << " to " << xmax << " in x\n"
				<< ymin << " to " << ymax << " in y\n"
				<< zmin << " to " << zmax << " in z" << endl;
			if (lwobble == 1)
				cout << "Range of thermal rms displacements (300K) = "
				<< wmin << " to " << wmax << endl;
	#endif	
			//@@@@@ start TEG code mods
			if(lwobble == 1 && wmin == 0 && wmax == 0)
				throw std::runtime_error("Input XYZ file does not contain thermal vibrations for atoms.");
			// force max dimensions along xzy axes to be equal to the defined CT sample qube side length
			if (xmin < 0 || ymin < 0 || zmin < 0)
				throw std::runtime_error("Error: xmin, ymin or zmin < 0 in the XYZ file.");
			if (xmax > ctblength || ymax > ctblength || zmax > ctblength)
				throw std::runtime_error("Error: xmax, ymax or zmax in the XYZ file is larger than the defined CT sample qube side length.");
			xmin = 0; xmax = ctblength;
			ymin = 0; ymax = ctblength;
			zmin = 0; zmax = ctblength;
	
			// Optionally add amorphous ice
			// NOTE that here the ice is added only once per illumination direction, so all images at different defocus distances at this direction see the same ice
			if (iceThick > 0)
			{
				if (iceThick < ctblength) 
					throw std::runtime_error("Error: ice thickness is less than the defined CT sample qube side length.");
				else 
				{
					natom = AddIce(iceThick, ctblength, natom, &Znum, &x, &y, &z, &occ, &wobble, wobbleaver, &iseed);
					cz = ctblengthz = iceThick;
				}
				printf("\nTotal atoms in added ice plus PDB structure = %d", natom);
				// Write out the modified XYZ file with the added O and H atoms (H2O molecules of the ice)
				// NOTE that the name of the output XYZ file is created inside this function by changing the extension of the filename passed to it
				if (nIceSave == 1) xar::SaveXYZfile(fileout[0], ctblength, ctblengthz, natom, Znum, x, y, z, occ, wobble);
			}
			//@@@@@ end TEG code mods


			// ---------  setup calculation -----
			//   set calculation flags
			aslice.lbeams = lbeams;
			aslice.lcross = lcross;
			aslice.lpartl = lpartl;
			aslice.lstart = lstart;
			aslice.lwobble = lwobble;

			//   set calculation parameters (some already set above)
			//@@@@@ start TEG code
			//param[ pAX ] = ax;			// supercell size
			//param[ pBY ] = by;
			param[pAX] = ctblength;			// supercell size
			param[pBY] = ctblength;
			//@@@@@ end TEG code
			param[pNX] = (float)nx;
			param[pNY] = (float)ny;
			param[pENERGY] = v0;
			param[pDELTAZ] = (float)deltaz;	// slice thickness
			param[pOAPERT] = aobj;
			param[pXCTILT] = ctiltx;		// crystal tilt
			param[pYCTILT] = ctilty;
			param[pCAPERT] = acmax;		// condencer angles
			param[pCAPERTMIN] = acmin;
			param[pTEMPER] = (float)fabs(temperature);
			param[pNWOBBLE] = (float)nwobble;	//  number config. to average
			param[pWAVEL] = wavlen;			//  probably recal. autoslice::calculate()

			param[pMODE] = 6;  // save mode = autoslic

			if (lpartl == 1) {
				param[pDEFOCUS] = df0;
				param[pOAPERT] = aobj;
				param[pDDF] = sigmaf;
				param[pCAPERT] = acmax;

				rx = 1.0F / ax;
				ry = 1.0F / by;
				cout << "Illumination angle sampling (in mrad) = "
					<< 1000. * rx * wavlen << ", " << 1000. * ry * wavlen << "\n" << endl;
			}

			// ------- iterate the multislice algorithm proper -----------
			//@@@@@ start TEG code
			if (lwobble != 1) nwobble = 1;
			// TEG introduced a cycle over thermal vibration configurations, with each thermal configuration step potentially including multiple defocus distances

			for (index_t iwobble = 0; iwobble < nwobble; iwobble++)
			{
				if (nwobble > 1)
					printf("\nThermal configuration no. %zd, iseed2 = %ld", iwobble, iseed2);

				if (noutput != 3) // calculate multislice propagation through the molecule
					aslice.calculate(pix, wave0, depthpix, param, multiMode, natom, &iseed2,
									Znum, x, y, z, occ, wobble, beams, hbeam, kbeam, nbout, ycross, dfdelt, ctblength, ctblengthz, nfftwinit, nmode);
			//@@@@@ end TEG code

				if (lpartl == 1) {         //    with partial coherence
					nillum = aslice.nillum;
					cout << "Total number of illumination angle = "
						<< nillum << endl;
					ndf = (int)((2.5F * sigmaf) / dfdelt);  // recal same value in class
					cout << "Total number of defocus values = " << 2 * ndf + 1 << endl;
				}

				else if (lbeams == 1) {
					fp1.open(filebeam.c_str());
					if (fp1.bad()) {
						cout << "can't open file " << filebeam << endl;
						exit(0);
					}
					fp1 << " (h,k) = ";
					for (ib = 0; ib < nbout; ib++)
						fp1 << " (" << hbeam[ib] << "," << kbeam[ib] << ")";
					fp1 << endl;
					nzbeams = beams.ny();
					fp1 << "nslice, (real,imag) (real,imag) ...\n" << endl;
					for (islice = 0; islice < nzbeams; islice++) {
						fp1 << setw(5) << islice;
						for (ib = 0; ib < nbout; ib++) //  setprecision(4)
							fp1 << "  " << setw(10) << beams.re(ib, islice)		//????? "%10.6f %10.6f",
							<< "  " << setw(10) << beams.im(ib, islice);
						fp1 << endl;
					}
					fp1.close();

				} // end else 

			/*  ------------------------------------------------------
				output results and find min and max to echo
				remember that complex pix are stored in the file in FORTRAN
					order for compatibility */

				pix.findRange(rmin, rmax, aimin, aimax);

				param[pRMAX] = rmax;
				param[pIMAX] = aimax;
				param[pRMIN] = rmin;
				param[pIMIN] = aimin;

				param[pDX] = dx = (float)(ax / ((float)nx));
				param[pDY] = dy = (float)(by / ((float)ny));

				//@@@@@ start TEG code
				IXAHWave2D* ph2new = CreateWavehead2D();
				ph2new->SetData(wavlen, ymin, ymax, xmin, xmax);
				xar::XArray2D<xar::fcomplex> camp(ny, nx);
				camp.SetHeadPtr(ph2new);
				xar::XArray2DFFT<float> xafft(camp);
				float xst = (float)GetXStep(camp);
				float yst = (float)GetYStep(camp);
				xar::XArray2D<double> Xr2, Yr2, Zr2; // squared coordinate differences
				xar::XArray1D<float> x2, y2, z2; // "wobbled" coordinates of atoms

				// precalculate arrays of squares of coordinate differences for subsequent evaluation of the electrostatic potential
				if (noutput == 3)
				{
					y2.Resize(natom); x2.Resize(natom); z2.Resize(natom);
					if (lwobble == 1 && nwobble > 1) 
					{
						float scale = (float)sqrt(temperature / 300.0);
						for (size_t k = 0; k < natom; k++) 
						{
							x2[k] = x[k] + (float)(wobble[k] * gasdev(&iseed2) * scale);
							y2[k] = y[k] + (float)(wobble[k] * gasdev(&iseed2) * scale);
							z2[k] = z[k] + (float)(wobble[k] * gasdev(&iseed2) * scale);
						}
					}
					else
					{
						for (size_t k = 0; k < natom; k++)
						{
							x2[k] = x[k];
							y2[k] = y[k];
							z2[k] = z[k];
						}
					}

					Yr2.Resize(ny, natom, 0.0); Xr2.Resize(nx, natom, 0.0); Zr2.Resize(defocus.size(), natom, 0.0);
					for (iy = 0; iy < ny; iy++)
					{
						yyy = ymin + yst * iy;
						for (size_t k = 0; k < natom; k++)
							Yr2[iy][k] = (y2[k] - yyy) * (y2[k] - yyy);
					}
					for (ix = 0; ix < nx; ix++)
					{
						xxx = xmin + xst * ix;
						for (size_t k = 0; k < natom; k++)
							Xr2[ix][k] = (x2[k] - xxx) * (x2[k] - xxx);
					}
					for (size_t j = 0; j < defocus.size(); j++)
					{
						zzz = zmax + (float)defocus[j].b;
						for (size_t k = 0; k < natom; k++)
							Zr2[j][k] = (z2[k] - zzz) * (z2[k] - zzz);
					} 
				}

				float k2maxo = aobj / wavlen;
				k2maxo = k2maxo * k2maxo;
				double C3 = double(Cs3 * 1.e+7); // mm --> Angstroms
				double C5 = double(Cs5 * 1.e+7); // mm --> Angstroms
				double dblYCentre = 0.5 * (ymax + ymin) + dyc, dblXCentre = 0.5 * (xmax + xmin) + dxc; // centre of rotation
				xar::XArray2D<float> ampRe0(ny, nx), ampIm0(ny, nx), ampRe(ny, nx), ampIm(ny, nx);
				ampRe.SetHeadPtr(ph2new->Clone()); ampIm.SetHeadPtr(ph2new->Clone());
				
				// former starting place of the cycle over different defocus distances	
				if (noutput != 3) // propagate in free space
				{
					// (re)define the transmitted complex amplitude (note the transposition!)
					for (ix = 0; ix < nx; ix++)
						for (iy = 0; iy < ny; iy++)
							camp[iy][ix] = xar::fcomplex(pix.re(ix, iy), pix.im(ix, iy));

					// take magnification into account if necessary
					float MM = 1.0f; // magnification factor
					double Rprime = defocus[jjj].b - 0.5 * ctblengthz; // effective propagation distance
					#if(0)
					{
						if (R1 != 0)
						{
							MM = (float)((R1 + defocus[jjj].b) / R1);
							if (MM == 0) throw std::runtime_error("Error in autosliccmd - magnification is equal to zero.");
							Rprime = defocus[jjj].b / MM - 0.5 * ctblengthz;
							// now remove the spherical wave component from the transmitted wave
							xar::fcomplex fctemp;
							for (int ix = 0; ix < nx; ix++)
								for (int iy = 0; iy < ny; iy++)
								{
									fctemp = xar::fcomplex(wave0.re(ix, iy), wave0.im(ix, iy));
									camp[iy][ix] /= fctemp;
									//camp[iy][ix] = fctemp;
								}
						}
						//xar::XArray2D<float> inten;
						//xar::CArg(camp, inten);
						//printf("\n  Writing output file %s ...", fileout[jjj].c_str());
						//xar::XArData::WriteFileGRD(inten, fileout[jjj].data(), xar::eGRDBIN);
						//exit(-22);
					}
					#endif
					
					// propagate
					xafft.FresnelA(Rprime, astigm.a, astigm.b * pi180, true, double(k2maxo), C3, C5); // propagate to the current defocus distance
					
					#if(0)
					{
						if (R1 != 0)
						{
							IXAHWave2D* ph2new = CreateWavehead2D();
							ph2new->SetData(wavlen, ymin * MM, ymax * MM, xmin * MM, xmax * MM);
							camp.SetHeadPtr(ph2new);
							camp /= float(MM);
						}
					}
					#endif
					
					#if(0) //!!!! temporary
					{
						xafft.FresnelA(-0.5 * ctblengthz, astigm.a, astigm.b * pi180, true, double(k2maxo), C3, C5); // propagate to the current defocus distance
						xar::XA_IWFR<float> abcd;
						xar::XArray2D<float> pha0;
						xar::CArg(camp, pha0);
						abcd.ForwardPhaseCTF(pha0, defocus[jjj].b, astigm.a, astigm.b* pi180, double(k2maxo), C3, C5);
						pha0 ^= 0.5;
						xar::MakeComplex<float>(pha0, 0, camp, true);
					}
					#endif
					

					#if(1)
					// old code, implementing Z" rotation and XY shifts for the defocused images
					// the new code instead applies Z" rotation and XY shifts to the molecule (XYZ structure) prior to the illumination
					// the new code is more accurate in principle, but is slower (in the case of multiple Z" rotations / defocus distances per illumination direction)
					// the new code also can create FFT reflection (periodicity) artefacts for atoms located close to the boundaries of the ctbox
					if (defocus[jjj].a != 0 || xyshifts.a != 0 || xyshifts.b != 0) // rotate around Z, then shift along X and Y
					{
						xar::Re(camp, ampRe0); xar::Im(camp, ampIm0);
						float averRe0 = (float)ampRe0.NormAverEdge(5);
						float averIm0 = (float)ampIm0.NormAverEdge(5);
						if (defocus[jjj].a != 0)
						{
							xar::XArray2DSpln<float> xaSplnRe(ampRe0), xaSplnIm(ampIm0);
							xaSplnRe.Rotate(ampRe, defocus[jjj].a, dblYCentre, dblXCentre, averRe0); // expecting uniform background
							xaSplnIm.Rotate(ampIm, defocus[jjj].a, dblYCentre, dblXCentre, averIm0); // expecting uniform background
						}
						else
						{
							ampRe = ampRe0;
							ampIm = ampIm0;
						}
						if (xyshifts.a != 0 || xyshifts.b != 0)
						{
							xar::XArray2DMove<float> xamoveRe(ampRe);
							xamoveRe.Move((long)floor(xyshifts.b / yst + 0.5), (long)floor(xyshifts.a / xst + 0.5), averRe0);
							xar::XArray2DMove<float> xamoveIm(ampIm);
							xamoveIm.Move((long)floor(xyshifts.b / yst + 0.5), (long)floor(xyshifts.a / xst + 0.5), averIm0);
						}
						xar::MakeComplex(ampRe, ampIm, camp, false);
					}
					#endif					
				}
				else // evaluate electrostatic potential
				{
					for (iy = 0; iy < ny; iy++)
						for (ix = 0; ix < nx; ix++)
							for (size_t k = 0; k < natom; k++)
								ampRe[iy][ix] += (float)vatom(Znum[k], sqrt(Xr2[ix][k] + Yr2[iy][k] + Zr2[jjj][k])); // add up potential contributions from all atoms in the molecule

					if (lwobble == 1 && nwobble == 1)
					// TEG: apply DW factor to the calculated projected potential
					{
						//printf("\n Calculating DW factor ...");
						// Debye-Waller factor implemented according to eqs.(4.4)-(4.6) from C.J. Humphreys, Rep.Prog.Phys., vol.42, 1979, pp.1827-1887.
						double scaleTemp2 = temperature / 300.0;
						double Bdw(0); // average Debye-Waller factor for all atoms in the structure / square(freq)
						for (i = 0; i < natom; i++) Bdw += wobble[i] * wobble[i] * scaleTemp2; // wobble[i] contains RMS displacement for an atom at 300K temperature
						Bdw /= natom;
						Bdw *= 2.0 * xar::tPI * xar::tPI;

						ampIm.Fill(0);
						xar::MakeComplex(ampRe, ampIm, camp, false);

						float emM;
						double mBdw4 = -0.25 * Bdw;

						xar::MakeComplex(ampRe, ampIm, camp, false);
						xafft.FFT(true);

						double kx, ky, dkx = 1.0 / ax, dky = 1.0 / by, kx2, ky2, k2;
						int ny2 = ny / 2, nx2 = nx / 2;
						for (iy = 0; iy < ny; iy++) 
						{
							ky = (iy - ny2) * dky;
							ky2 = ky * ky;
							for (ix = 0; ix < nx; ix++) 
							{
								kx = (ix - nx2) * dkx;
								kx2 = kx * kx;
								k2 = ky2 + kx2;
								emM = (float)exp(mBdw4 * k2); // exp(-M)
								camp[iy][ix] *= emM;
							}
						
						}
						xafft.FFT(false);
						xar::Re(camp, ampRe);
					}

					if (defocus[jjj].a != 0) // rotate around Z
					{
						ampRe0 = ampRe;
						float averRe0 = (float)ampRe0.NormAverEdge(5);
						xar::XArray2DSpln<float> xaSplnRe(ampRe0);
						xaSplnRe.Rotate(ampRe, defocus[jjj].a, dblYCentre, dblXCentre, averRe0); // expecting uniform background
					}
				}

				//TIFF/GRD/GRC file output
				switch (noutput)
				{
				case 0: // intensity out
				{
					xar::XArray2D<float> inten;
					xar::Abs2(camp, inten);
					if (iwobble > 0)
					{
						xar::XArray2D<float> inten1;
						if (boutTIFF)
						{
							xar::TIFFReadFile(inten1, fileout[jjj].data());
							inten1.SetHeadPtr(ph2new->Clone());
						}
						else xar::XArData::ReadFileGRD(inten1, fileout[jjj].data(), wavlen);
						inten += inten1;
						if (iwobble == (nwobble - 1)) inten /= (float)nwobble;
					}
					printf("\n  Writing output file %s ...", fileout[jjj].c_str());
					if (boutTIFF) xar::TIFFWriteFile(inten, fileout[jjj].data(), xar::eTIFF32, false);
					else xar::XArData::WriteFileGRD(inten, fileout[jjj].data(), xar::eGRDBIN);
					break;
				}
				case 1: // phase out
				{
					xar::XArray2D<float> phase(camp.GetDim1(), camp.GetDim2(), 0.0f);
					xar::CArg(camp, phase);
					printf("\n  Writing output file %s ...", fileout[jjj].c_str());
					if (boutTIFF) xar::TIFFWriteFile(phase, fileout[jjj].data(), xar::eTIFF32, false);
					else xar::XArData::WriteFileGRD(phase, fileout[jjj].data(), xar::eGRDBIN);
					break;
				}
				case 2: // complex amplitude out
				{
					printf("\n  Writing output file %s ...", fileout[jjj].c_str());
					xar::XArData::WriteFileGRC(camp, fileout[jjj].data(), xar::eGRCBIN);
					break;
				}
				case 3: // 3D potential out
				{
					if (iwobble > 0)
					{
						if (boutTIFF)
						{
							xar::TIFFReadFile(ampRe0, fileout[jjj].data());
							ampRe0.SetHeadPtr(ph2new->Clone());
						}
						else xar::XArData::ReadFileGRD(ampRe0, fileout[jjj].data(), wavlen);
						ampRe += ampRe0;
						if (iwobble == (nwobble - 1)) ampRe /= (float)nwobble;
					}
					printf("\n  Writing output file %s ...", fileout[jjj].c_str());
					if (boutTIFF) xar::TIFFWriteFile(ampRe, fileout[jjj].data(), xar::eTIFF32, false);
					else xar::XArData::WriteFileGRD(ampRe, fileout[jjj].data(), xar::eGRDBIN);
					break;
				}
				default:
					throw std::runtime_error("Error: unknown value of the output mode parameter.");
				}

			} // end of cycle over defocus distance


		} // end of cycle over iwobble (thermal configurations)


		/*   for( ix=0; ix<NPARAM; ix++ ) myFile.setParam( ix, param[ix] );

		   if ( lpartl == 1 ) {
			   myFile.resize( nx, ny );
			   myFile.setnpix( 1 );
			   for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
				   myFile(ix,iy) = pix.re(ix,iy);
		   } else {
			   myFile.resize( 2*nx, ny );
			   myFile.setnpix( 2 );
			   for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
				   myFile(ix,iy)    = pix.re(ix,iy);
				   myFile(ix+nx,iy) = pix.im(ix,iy);
			   }
		   }

		   i = myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );
		   if( i != 1 ) cout << "autoslice cannot write TIF file " << fileout << endl;
		   cout << "pix range " << rmin << " to " << rmax << " real,\n" <<
				   "          " << aimin << " to " << aimax << " imag" << endl;
	   */
	   //@@@@@ end TEG code

	   /* ----- output depth cross section if requested ------- */
		if (lcross == 1) {
			depthpix.findRange(rmin, rmax, aimin, aimax);
			myFile.setParam(pRMAX, rmax);
			myFile.setParam(pIMAX, 0.0F);
			myFile.setParam(pRMIN, rmin);
			myFile.setParam(pIMIN, 0.0F);
			myFile.setParam(pDY, dy = (float)(deltaz));

			nzout = depthpix.ny();
			myFile.resize(nx, nzout);
			myFile.setnpix(1);
			for (ix = 0; ix < nx; ix++) for (iz = 0; iz < nzout; iz++) {
				myFile(ix, iz) = depthpix.re(ix, iz);
			}
			i = myFile.write(filecross.c_str(), rmin, rmax, aimin, aimax, dx, dy);

			if (i != 1) cout << "autoslice cannot write TIF file "
				<< filecross << endl;
			cout << "depth pix range " << rmin << " to " << rmax << " real" << endl;
		}
#ifdef _DEBUG
		cout << "Total CPU time = " << cputim() - timer << " sec." << endl;
#ifdef USE_OPENMP
		//cout << "wall time = " << walltim() - walltimer << " sec." << endl;
#endif
#endif
		//@@@@@ start TEG code
		//char a;
		//cout << "\nPress any key to exit ...";
		//cin >> a;

		if (thread_counter.GetCount() < 0)
			throw std::runtime_error("another thread has requested program termination.");
		//@@@@@ end TEG code
	}
	catch (std::runtime_error& E)
	{
#ifdef TEG_MULTITHREADED
		printf("\n!!!Exception: %s\n", E.what());
		thread_counter.SetTerminate();
#else
		std::rethrow_exception(std::current_exception());
#endif // TEG_MULTITHREADED
	}

    return 0;

} /* end main() */


