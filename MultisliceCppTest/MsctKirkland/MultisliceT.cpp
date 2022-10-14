// MultisliceCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <thread>
#include <chrono>
#include <omp.h>

#include "autosliccmd.h"
#include "pdb.h"
#include "XAHWave.h"
#include "XArray2D.h"
#include "XArray3D.h"
#include "XA_spln2.h"
#include "XA_move2.h"
#include "XA_spln3.h"
#include "XA_data.h"
#include "XA_tiff.h"

using namespace xar;

int Counter_Obj::counter = 0; // active worker thread counter (this main thread is not counted)
bool Counter_Obj::isUpdated = false; // update status of the counter object


int main(int argc, char* argv[])
{
#ifdef TEG_MULTITHREADED
	Counter_Obj thread_counter; // increments the thread counter on construction and decrements it on destruction
#endif // TEG_MULTITHREADED
	vector<string> autoslictxt(32); // 31 is the current number of input parameters; if it is changed, the corresponding changes need to be applied in autosliccmd.cpp too.

	printf("\nStarting MsctKirkland program ...");
	try
	{
		// read input parameter file
		string sInputParamFile("MsctKirkland.txt");
		if (argc > 1) sInputParamFile = argv[1];
		FILE* ff0 = fopen(sInputParamFile.c_str(), "rt");
		if (!ff0) throw std::runtime_error(string("Error: cannot open parameter file " + sInputParamFile + ".").c_str());
		else printf("\nReading input parameter file %s ...", sInputParamFile.c_str());

		char cline[1024], ctitle[1024], cparam[1024], cparam1[1024], cparam2[1024];
	
		// read and skip an arbitrary number of initial comment lines (i.e. the lines that start with // symbols)
		while (true)
		{
			fgets(cline, 1024, ff0); 
			if (!(cline[0] == '/' && cline[1] == '/')) break;
		}

		strtok(cline, "\n"); // 1st parameter: Input_file_with_atomic_numbers_and_coordinates_in_XYZ_format
		autoslictxt[0] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading input file name from the input parameter file.");
		int iModality;
		string sInputFile(cparam);
		if (GetFileExtension(sInputFile) == string(".XYZ"))
		{
			iModality = 0;
			printf("\nInput file %s shall contain contain atomic numbers and coordinates. TEM imaging is assumed.", sInputFile.c_str());
		}
		else if (GetFileExtension(sInputFile) == string(".GRD"))
		{
			iModality = 1;
			printf("\nInput files defined by the template %s shall contain 3D distribution of beta. X-ray imaging is assumed.", sInputFile.c_str());
		}
		else throw std::runtime_error("Error: unknown input file extension in the parameter file");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 2nd parameter: delta / beta
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading delta/beta from input parameter file.");
		double delta2beta = atof(cparam);
		printf("\nDelta/beta = %g.", delta2beta);
		// we ignore this parameter in the case of TEM imaging
		if (delta2beta <= 0 && iModality == 1) throw std::runtime_error("Error: delta/beta must be positive in hard X-ray imaging.");
		
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 3rd parameter: Output_files_shall_contain_intensity(0),_phase(1),_complex_amplitude(2)_or_3D_potential(3)
		autoslictxt[27] = cline; // The numbering of these parameters is 'historic', it can be changed, but then the corresponding changes need to be applied in autosliccmd.cpp too.
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading output type from the input parameter file.");
		int nOutputType = atoi(cparam);
		switch (nOutputType)
		{
		case 0:
			printf("\nOutput data type will be image intensity.");
			break;
		case 1:
			printf("\nOutput data type will be phase distribution.");
			break;
		case 2:
			printf("\nOutput data type will be complex amplitude.");
			break;
		case 3:
			printf("\nOutput data type will be 3D potential.");
			if (iModality == 1) throw std::runtime_error("Error: 3D potential cannot be output in the X-ray imaging mode.");
			break;
		default:
			throw std::runtime_error("Error: unrecognised output data type in the input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 4th parameter: Output_TIFF/GRD/GRC_filename_template
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading output file template from the input parameter file.");
		string outfilename(cparam); // output filename stub
		printf("\nOutput filename template is %s.", cparam);
		if ((GetFileExtension(outfilename) != string(".TIFF")) && (GetFileExtension(outfilename) != string(".TIF")) && (GetFileExtension(outfilename) != string(".GRD")) && (GetFileExtension(outfilename) != string(".GRC")))
			throw std::runtime_error("Error: output filename extension must be TIF, GRD or GRC.");
		if ((nOutputType == 0 || nOutputType == 1 || nOutputType == 3) && !(GetFileExtension(outfilename) == string(".GRD") || GetFileExtension(outfilename) == string(".TIF") || GetFileExtension(outfilename) == string(".TIFF")))
			throw std::runtime_error("Output type (parameter 2) is inconsistent with output filename extension in the input parameter file.");
		if (nOutputType == 2 && (GetFileExtension(outfilename) != string(".GRC")))
			throw std::runtime_error("Output type (parameter 2) is inconsistent with output filename extension in the input parameter file.");
		bool boutTIFF = (GetFileExtension(outfilename) == string(".TIFF"));

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 5th paramter: Use_multislice(0),_projection(1),_or_1st_Born(2)_approximation
		autoslictxt[25] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading calculation mode from the input parameter file.");
		int nCalculationMode = atoi(cparam);
		switch (nCalculationMode)
		{
		case 0:
			printf("\nCalculations will use the multislice model.");
			break;
		case 1:
			printf("\nCalculations will use the projection approximation.");
			break;
		case 2:
			printf("\nCalculations will use the first Born approximation.");
			break;
		default:
			throw std::runtime_error("Error: unknown calculation mode in the input parameter file.");
		}
	
		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 6th parameter: Incident__electron_or_X-ray_beam_energy_in_keV
		autoslictxt[10] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading electron / X-ray beam energy from the input parameter file.");
		double kev = atof(cparam);
		printf("\nIncident electron or X-ray beam energy is %lg.", kev);
		double wl; // X-ray wavelength
		if (iModality == 0)
		{
			double ev = kev * 1000; // accelerating voltage in eV
			constexpr double hp = 6.62607004e-34; // Planck's constant (m2 kg / s)
			constexpr double cc = 299792458; // speed of light (m / s)
			constexpr double ee = 1.602176634e-19; // electron charge (coulomb)
			constexpr double m0 = 9.1093837015e-31; // electron rest mass (kg)
			double ewl = hp * cc / sqrt(ee * ev * (2.0 * m0 * cc * cc + ee * ev));
			printf("\nElectron wavelength = %f (pm)", ewl * 1.0e+12);
		}
		else
		{
			wl = 12.4 / kev * 1.e-4; // X-ray wavelength in microns
			printf("\nIncident X-ray wavelength is %lg (A).", wl * 1.e+4);
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 7th parameter: Distance_from_point_source_in_Angstroms_(0=plane_wave)
		autoslictxt[31] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading distance from the source from the input parameter file.");
		printf("\nDistance from the source to the centre of the object is %lg Angstroms.", atof(cparam));

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 8th parameter: Wavefunction_size_in_pixels,_Nx,Ny
		autoslictxt[11] = cline;
		index_t nx, ny, nz;
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading wavefunction size from the input parameter file.");
		nx = atoi(cparam);
		ny = atoi(cparam1);
		nz = nx;
		printf("\nWavefunction size is nx = %zd, ny = %zd pixels.", nx, ny);
		if (iModality == 1) printf("\n3D distribution of beta has dimensions nx = %zd, ny = %zd, nz = %zd pixels.", nx, ny, nz);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 9th parameter: Slice_thickness_in_Angstroms
		autoslictxt[13] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading slice thickness from the input parameter file.");
		double sliceTh = atof(cparam);
		printf("\nSlice thickness is %lg Angstroms.", sliceTh);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 10th parameter:Objective_aperture_in_mrad
		autoslictxt[7] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading objective aperture from the input parameter file.");
		double aobj = atof(cparam);
		if (atof(cparam) != 0) printf("\nObjective aperture is %lg mrad.", aobj);
		else printf("\nObjective aperture is unrestricted.");
		double k2maxo = pow(aobj * 0.001f / wl, 2.0); // Fourier space bandwidth

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 11th parameter: Spherical_aberration_Cs3_and_Cs5_in_mm
		autoslictxt[5] = cline;
		if (sscanf(cline, "%s %s %s", ctitle, cparam, cparam1) != 3) throw std::runtime_error("Error reading spherical aberrations from the input parameter file.");
		double Cs3 = atof(cparam);
		double Cs5 = atof(cparam1);
		printf("\nSpherical aberrations: Cs3 = %g, Cs5 = %g (mm)", Cs3, Cs5);
		Cs3 *= 1.e+7; // mm --> Angstroms
		Cs5 *= 1.e+7; // mm --> Angstroms

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 12th paramter: Include_thermal_vibrations(1)_or_not(0)
		autoslictxt[17] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading thermal vibration option from the input parameter file.");
		int iThermalVibr = atoi(cparam);
		if (iThermalVibr && iModality) throw std::runtime_error("Error: thermal vibrations cannot be used in the X-ray imaging mode.");
		switch (iThermalVibr)
		{
		case 0:
			printf("\nThermal vibrations will not be included in the calculatioins.");
			break;
		case 1:
			printf("\nThermal vibrations will be included in the calculatioins.");
			break;
		default:
			throw std::runtime_error("Error: unknown value of the thermal vibrations switch in the input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 13th parameter: ____Temperature_in_degrees_K
		autoslictxt[18] = cline;
		double dTemperature(-1.0);
		if (iThermalVibr != 0)
		{
			if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading temperature from the input parameter file.");
			dTemperature = atof(cparam);
			printf("\nTemperature is %lg in degrees K.", dTemperature);
			if (dTemperature < 0) throw std::runtime_error("Temperature in degrees K cannot be negative.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 14th parameter: ____Number_of_configurations_to_average_over
		autoslictxt[19] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading the number of configurations parameter from the input parameter file.");
		int iNumConfig = atoi(cparam);
		printf("\nNumber of configurations to use for thermal vibration calculations is %d.", iNumConfig);
		if (iNumConfig < 0) throw std::runtime_error("Number of configurations cannot be negative.");
		if (iThermalVibr != 0 && iNumConfig > 1 && !(nOutputType == 0 || nOutputType == 3))
			throw std::runtime_error("Output type can only be 0 (intensity) or 3 (3D potential) when the thermal vibrations option with multiple configurations is enabled in the input parameter file.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 15th parameter: Ice_layer_thickness_in_Angstroms
		autoslictxt[29] = cline;
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading ice layer thickness from the input parameter file.");
		double dIceThickness = atof(cparam);
		printf("\nIce layer thickness is %lg Angstroms.", dIceThickness);
		if (dIceThickness > 0 && iModality) throw std::runtime_error("Error: amorphous ice cannot be added in the X-ray imaging mode.");
		if (dIceThickness > 0 && dTemperature > 273.15) throw std::runtime_error("Temperature cannot be higher than 273.15K in the presence of an ice layer.");

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //16th parameter: Save_XYZ_files_with_ice:_no(0), _yes(1), _1st_one_only(2)
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading save XYZ file with ice switch parameter from the input parameter file.");
		int nIceSave = atoi(cparam);
		switch (nIceSave)
		{
		case 0:
			if (dIceThickness > 0) printf("\nXYZ files with ice will not be saved.");
			autoslictxt[30] = "31.Save_XYZ_files_with_ice:_no(0),_yes(1): " + string("0");
			break;
		case 1:
			if (dIceThickness > 0) printf("\nXYZ files with ice will be saved.");
			autoslictxt[30] = "31.Save_XYZ_files_with_ice:_no(0),_yes(1): " + string("1");
			break;
		case 2:
			if (dIceThickness > 0) printf("\nOnly the first XYZ file with ice will be saved.");
			// if the value of this parameter == 2, then this array element will be written later - see below
			break;
		default:
			throw std::runtime_error("Unknown value of the save XYZ file with ice switch parameter from the input parameter file.");
		}

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 17th parameter: Text_file_with_output_rotation_angles_in_degrees_and_defocus_distances_in_Angstroms
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading the name of the file with defocus distances and rotational positions from the input parameter file.");
		printf("\nText file with output rotation angles and defocus distances is %s.", cparam);
		string DefocFile = string(cparam);
		vector<Pair> v2angles;
		vector<vector <Pair> > vvdefocus;
		vector<Pair> vastigm;
		vector<Pair> v2shifts;
		bool bRelion;
		if (GetFileExtension(DefocFile) == string(".TXT"))
		{
			bRelion = false;
			ReadDefocusParamsFile(DefocFile, v2angles, vvdefocus, true);
			vastigm.resize(v2angles.size(), xar::Pair{ 0.0, 0.0 });
			v2shifts.resize(v2angles.size(), xar::Pair{ 0.0, 0.0 });
		}
		else
			if (GetFileExtension(DefocFile) == string(".RELIONNEW"))
			{
				bRelion = true;
				//!!! note that the ReadRelionDefocusParameters currently can have only one defocus distance per illumination angle
				ReadRelionDefocusParamsFile(DefocFile, v2angles, vvdefocus, vastigm, v2shifts, true);
			}
			else throw std::runtime_error("Error: unrecognised filename extension in parameter 15 of the input parameter file.");
		//!!!NOTE that we assume that the defocus distances in the orientaion parameter file are given from the centre of the molecule,
		// and we adjust these defocus distances accordingly in the free-space propagation part of autosliccmd.cpp, where we subtract half of the z-thickness of the structure 
		//(which is already taken into account in the multislice algorithm) from each defocus distance given in the orientation file

		index_t nangles = v2angles.size(); // number of different illumination directions 
		vector<index_t> vndefocus(nangles); // vector of numbers of rotations around the optic axis = number of defocus planes at different rotation angles
		for (index_t i = 0; i < nangles; i++) vndefocus[i] = vvdefocus[i].size();
		if (GetFileExtension(outfilename) == string(".GRC") || nOutputType == 3) // check that there is only one input defocused complex amplitude file per each illumination angle
		{
			for (index_t i = 0; i < nangles; i++)
				if (vndefocus[i] != 1)
					throw std::runtime_error("Error: only one defocus distance per row is allowed in the case of output complex amplitudes or 3D potential calculations.");
		}
		vector<string> voutfilenamesTot;
		FileNames2(vndefocus, outfilename, voutfilenamesTot); // create "2D array" of output filenames

		fgets(cline, 1024, ff0); strtok(cline, "\n"); //18th parameter: Centre of rotation x, y and z shifts in Angstroms
		autoslictxt[12] = cline;
		if (sscanf(cline, "%s %s %s %s", ctitle, cparam, cparam1, cparam2) != 4)
			throw std::runtime_error("Error reading ccentre of rotation x, y and z shifts from input parameter file.");
		double dxc = atof(cparam), xc;
		double dyc = atof(cparam1), yc;
		double dzc = atof(cparam2), zc;
		printf("\nCentre of rotation shifts in Angstroms: dXc = %g, dYc = %g, dZc = %g (Angstroms)", dxc, dyc, dzc);

		fgets(cline, 1024, ff0); strtok(cline, "\n"); // 19th parameter: Number_of_worker_threads_to_launch_in_CT_simulation_mode
		if (sscanf(cline, "%s %s", ctitle, cparam) != 2) throw std::runtime_error("Error reading the number of worker threads from the input parameter file.");
		int nThreads = atoi(cparam); // number of threads to use (expected to be equal to the number of cores) 
		printf("\nMaximum number of threads to be used = %d", nThreads);
		if (nThreads < 1)
			throw std::runtime_error("Error: the number of parallel threads in input parameter file should be >= 1.");
		omp_set_num_threads(nThreads);

		fclose(ff0); // close input parameter file

		double xlo, xhi, xst, ylo, yhi, yst, zlo, zhi, zst;
		XArray3D<float> Beta3D; // 3D distribution of beta
		if (iModality == 1)
		{
			printf("\n\nReading 3D distribution of beta from files ...");

			vector<string> vinfilenamesTot; // input filenames for beta
			FileNames(1, nz, sInputFile, vinfilenamesTot); // create 1D array of input filenames to save 2D slices of beta

			Beta3D.Resize(nz, ny, nx);

			// read input 2D slices
			bool bAbort(false);
			#pragma omp parallel for shared (Beta3D)
			for (int n = 0; n < nz; n++)
			{
				if (bAbort) continue;
				try
				{
					XArray2D<float> ipIn;
					XArData::ReadFileGRD(ipIn, vinfilenamesTot[n].c_str(), wl); //	read input GRD files
					if (n == 0)
					{
						const IXAHWave2D* ph2 = dynamic_cast<const IXAHWave2D*>(ipIn.GetHeadPtr());
						ylo = ph2->GetYlo();
						yhi = ph2->GetYhi();
						yst = GetYStep(ipIn);
						xlo = ph2->GetXlo();
						xhi = ph2->GetXhi();
						xst = GetXStep(ipIn);
						zlo = xlo;
						zhi = xhi;
						zst = xst;
						yc = (ylo + yhi) / 2.0 + dyc;
						xc = (xlo + xhi) / 2.0 + dxc;
						zc = (zlo + zhi) / 2.0 + dzc;
						Beta3D.SetHeadPtr(new Wavehead3D(wl, zlo, zhi, ylo, yhi, xlo, xhi));
					}

					index_t ny0 = ipIn.GetDim1();
					index_t nx0 = ipIn.GetDim2();
					if (nx0 != nx || ny0 != ny)
							throw std::runtime_error("Error: the dimensions of the 2D array in the image file are inconsistent with the array dimensions in input parameter file.");

					for (index_t j = 0; j < ny; j++)
						for (index_t i = 0; i < nx; i++)
							Beta3D[n][j][i] = ipIn[j][i];
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
			printf("\n@@@@@ Beta3D-aver = %g", Beta3D.Norm(eNormAver));
		}

		//!!The following additional parameters are required in autosliccmd.cpp, but are not currently used by this program
		autoslictxt[1] = "2.Replicate_unit_cell_by_NCELLX,NCELLY,NCELLZ: 1 1 1";
		autoslictxt[2] = ""; // this parameter is not used any more, output filenames are passed as a separate argument vstrfileout
		autoslictxt[3] = "4.Do_you_want_to_include_partial_coherence: 0";
		autoslictxt[4] = "5.____Illumination_angle_min,_max_in_mrad: 0.0 0.0";
		autoslictxt[6] = "7.____Defocus_mean,_standard_deviation,_and_sampling_size_in_Angstroms: 0.0 0.0 0.0";
		autoslictxt[8] = "9.Do_you_want_to_start_from_previous_result: 0";
		autoslictxt[9] = "10.____Name_of_file_to_start_from: 0";
		//autoslictxt[12] = "13.Crystal_tilt_x,y_in_mrad: 0.0 0.0"; // !!! NOTE: we have reused this parameter for the centre of rotation shifts (see above) 
		autoslictxt[14] = "15.Do_you_want_to_record_the_(real,imag)_value_of_selected_beams_vs._thickness: 0";
		autoslictxt[15] = "16.____Name_of_file_for_beams_info: 0";
		autoslictxt[16] = "17.____Number_of_beams: 0";
		autoslictxt[20] = "21.____Initial_seed_for_random_number_generator: 1";
		autoslictxt[21] = "22.Do_you_want_to_output_intensity_vs._depth_cross_section: 0";
		autoslictxt[22] = "23.____Name_of_file_to_get_depth_profile_image: 0";
		autoslictxt[23] = "24.____y_position_of_depth_cross_section_in_Angstroms: 0.0";
		autoslictxt[26] = ""; // this parameter is not used any more, defocus values are passed as a separate argument vdefocus		

		// find out the number of CPU cores available in the computer
		//unsigned int nThreads = std::thread::hardware_concurrency();
		//printf("\nNumber of CPU cores detected: %d", nThreads);

		// start the execution timer
		std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
	
		// start the cycle over projection angles
		if (iModality == 0) //TEM imaging case
		{
			char bufangle[1024];
			index_t ndefcurrent(0);

			for (size_t i = 0; i < nangles; i++)
			{
				// start the execution timer
				std::chrono::system_clock::time_point start_time_angle = std::chrono::system_clock::now();

				Pair angle = v2angles[i];
				printf("\nIllumination angle: z = %g, y' = %g (degrees)", angle.a, angle.b);
				sprintf(bufangle, "%f %f", angle.a * PI180, angle.b * PI180);
				autoslictxt[24] = "25.Sample_Z_and_Y'_rotation_angles_in_radians: " + string(bufangle);

				if (nIceSave == 2)
					if (i == 0) autoslictxt[30] = "31.Save_XYZ_files_with_ice:_no(0),_yes(1): " + string("1");
					else autoslictxt[30] = "31.Save_XYZ_files_with_ice:_no(0),_yes(1): " + string("0");

				// start the cycle over defocus distances (we only create output file names in this inner cycle)
				index_t ndefocus = vndefocus[i]; // number of defocus planes at the current illumination angle
				vector<Pair> vdefocus = vvdefocus[i]; // vector of Z'' angles and defocus distances at the current illumination angle
				// Note that vdefocus[n].a should contain Z'' angle in degrees, because xar::XArray2DSpln<T>::Rotate function expects the rotation angle in degrees
				vector<string> vstrfileout(ndefocus); // vector of output filenames at the current illumination angle
				for (index_t n = 0; n < ndefocus; n++) vstrfileout[n] = voutfilenamesTot[ndefcurrent++];


				//Here we call Kirkland's autoslic at each angle
				autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 1"; // the first thread must initialize the FFTW plan, subsequent ones can copy it
#ifdef TEG_MULTITHREADED
				// NOTE that in this multithreaded model, exceptions thrown by worker threads will NOT be caught in the master thread
				if (i > 0) autoslictxt[28] = "29.Copy(0)_or_initialize(1)_FFTW_plan: 0";
				thread_counter.SetUpdated(false);
				std::thread threadObj(autosliccmd, autoslictxt, vdefocus, vastigm[i], v2shifts[i], vstrfileout);
				if (i == 0) threadObj.join(); // we need to let the first worker thread finish execution, so that it can create the FFTW "plan" to be shared  with other threads
				if (threadObj.joinable()) threadObj.detach(); // if we don't do this, threadObj will call Terminate() on the attached thread when the threadObj goes out of scope
				while (thread_counter.GetCount() >= 0 && (!thread_counter.GetUpdated() || thread_counter.GetCount() >= (int)nThreads))
					std::this_thread::sleep_for(std::chrono::milliseconds(10)); // we allow nThreads of threads to be launched
				if (thread_counter.GetCount() < 0) // meaning that a thread has requested the whole program to terminate
					throw std::runtime_error("A thread has requested the whole program to terminate.");
#else
				autosliccmd(autoslictxt, vdefocus, vstrfileout); // single-threaded execution mode
#endif // TEG_MULTITHREADED
				std::chrono::system_clock::time_point end_time_angle = std::chrono::system_clock::now();
				printf("\nExecution time for this illumination angle = %ld ms.", (long int)std::chrono::duration_cast<std::chrono::milliseconds>(end_time_angle - start_time_angle).count());
			}

#ifdef TEG_MULTITHREADED
			while (thread_counter.GetCount() > 1)
				std::this_thread::sleep_for(std::chrono::milliseconds(100)); // wait for all worker threads to finish or fail
#endif // TEG_MULTITHREADED

		}
		else // X-ray imaging case
		{
			XArray3DSpln<float> spln3d(Beta3D);
			int nThreads1 = int(getAvailableMemory() / (Beta3D.size() * sizeof(float)) * 0.8);
			nThreads1 = min(nThreads, nThreads1);
			omp_set_num_threads(nThreads1);

			if (!fftw_init_threads())
				std::runtime_error("Error: multi-threaded initialization of FFTW library failed.");
			// single Fftwd2c object is used for repeated 2D FFTs using FFTW later (the plan is created at this point and reused later)
			Fftwd2fc fftw2Df((int)ny, (int)nx, 1, true); // we don't use mutlithreaded 2D FFT here, as the code calling these functions below is already mutlithreaded

			bool bAbort(false);
			#pragma omp parallel for
			for (int na = 0; na < nangles; na++)
			{
				if (bAbort) continue;
				try
				{
					XArray3D<float> Beta3Dcopy(nz, ny, nx); // this needs to be different for each illumination angle - big problem for parallelization!!!
					Pair angle = v2angles[na];
					printf("\nIllumination angle: z = %g, y' = %g (degrees)", angle.a, angle.b);

					vector<Pair> vdefocus = vvdefocus[na];
					if (vdefocus.size() != 1) throw std::runtime_error("\nError: only the case of a single image per illumination direction has been implemented so far!");

					XArray2D<fcomplex> campOut(ny, nx, fcomplex(1.0f, 0.0f)); // "incident" complex amplitude
					campOut.SetHeadPtr(new Wavehead2D(wl, ylo, yhi, xlo, xhi));
					spln3d.Rotate3(Beta3Dcopy, angle.a * PI180, angle.b * PI180, 0.0, nThreads, zc, yc, xc, 0.0f);
					XArray3DSpln<float> spln3dA(Beta3Dcopy);
				
					spln3dA.Multislice_Hom(campOut, fftw2Df, delta2beta, sliceTh, k2maxo);

					// propagate to defocuse planes, rotate around the illumination axis and save the defocused intensities
					XArray2D<float> ampRe, ampIm;
					XArray2DFFT<float> xafft(campOut, false, true);
					//int iSign = vdefocus[n].b < 0 ? iSign = -1 : iSign = 1;
					int iSign = +1; // now we want to reproduce(!) the known aberrations, rather than compensate them
					if (bRelion)
						xafft.FresnelB(fftw2Df, vdefocus[0].b - 0.5 * nz * xst, iSign * vastigm[na].a, iSign * vastigm[na].b * PI180, true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to the current defocus distance
					else
						xafft.FresnelB(fftw2Df, vdefocus[0].b - 0.5 * nz * xst, 0, 0, true, k2maxo, iSign * Cs3, iSign * Cs5); // propagate to the current defocus distance
					if (vdefocus[0].a != 0 || bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0)) // rotate around Z'', then shift along XY, if needed
					{
						Re(campOut, ampRe); Im(campOut, ampIm);
						float averRe = float(ampRe.NormAverEdge(5));
						float averIm = float(ampIm.NormAverEdge(5));
						// rotate input defocused complex amplitude around Z''
						if (vdefocus[0].a != 0)
						{
							XArray2DSpln<float> xaSplnRe(ampRe), xaSplnIm(ampIm);
							ampRe = xaSplnRe.Rotate(vdefocus[0].a, yc, xc, averRe); // expecting uniform background
							ampIm = xaSplnIm.Rotate(vdefocus[0].a, yc, xc, averIm); // expecting uniform background
						}
						// shift along X and/or Y
						if (bRelion && (v2shifts[na].a != 0 || v2shifts[na].b != 0))
						{
							XArray2DMove<float> xamoveRe(ampRe), xamoveIm(ampIm);
							xamoveRe.Move((long)floor(v2shifts[na].b / yst + 0.5), (long)floor(v2shifts[na].a / xst + 0.5), averRe);
							xamoveIm.Move((long)floor(v2shifts[na].b / yst + 0.5), (long)floor(v2shifts[na].a / xst + 0.5), averIm);
						}
						MakeComplex(ampRe, ampIm, campOut, false);
					}

					//TIFF/GRD/GRC file output
					switch (nOutputType)
					{
					case 0: // intensity out
					{
						XArray2D<float> inten;
						Abs2(campOut, inten);
						printf("\n  Writing output file %s ...", voutfilenamesTot[na].c_str());
						if (boutTIFF) TIFFWriteFile(inten, voutfilenamesTot[na].data(), xar::eTIFF32, false);
						else XArData::WriteFileGRD(inten, voutfilenamesTot[na].data(), xar::eGRDBIN);
						break;
					}
					case 1: // phase out
					{
						XArray2D<float> phase(campOut.GetDim1(), campOut.GetDim2(), 0.0f);
						CArg(campOut, phase);
						printf("\n  Writing output file %s ...", voutfilenamesTot[na].c_str());
						if (boutTIFF) TIFFWriteFile(phase, voutfilenamesTot[na].data(), xar::eTIFF32, false);
						else XArData::WriteFileGRD(phase, voutfilenamesTot[na].data(), xar::eGRDBIN);
						break;
					}
					case 2: // complex amplitude out
					{
						printf("\n  Writing output file %s ...", voutfilenamesTot[na].c_str());
						XArData::WriteFileGRC(campOut, voutfilenamesTot[na].data(), xar::eGRCBIN);
						break;
					}
					default:
						throw std::runtime_error("Error: unknown value of the output mode parameter.");
					}
				}
				catch (std::exception& E)
				{
					printf("\n\n!!!Exception: %s\n", E.what());
					bAbort = true;
				}
			}
			if (bAbort) throw std::runtime_error("at least one thread has thrown an exception.");
		}

		std::chrono::system_clock::time_point end_time = std::chrono::system_clock::now();
		printf("\nMain program finished. Execution time = %ld s.", (long int)std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count());

	}
	catch (std::runtime_error& E)
	{
		printf("\n!!!Exception: %s\n", E.what());
	}
	catch (std::exception& E)
	{
		printf("\n\n!!!Exception: %s\n", E.what());
	}


	printf("\nPress any key to exit..."); getchar();
	return 0;
}
