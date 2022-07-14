Brief description of the PDB-XYZ-XTL file translation module "Pdb.exe"
of the Differential Holographic Tomography (DHT) software package

This program is based on A.V. Martin's PDB code. It also uses the "miniball" code written by Martin Kutz <kutz@math.fu-berlin.de>
and Kaspar Fischer <kf@iaeth.ch> for finding the minimal ball enclosing the structure
// See https://github.com/hbf/miniball for authors, theory and copyright information.
It also uses E.J. Kirkland's temsim code from https://sourceforge.net/projects/computem/. 

Pdb.exe module performs PDB and XYZ atomic structure file translation tasks, producing the output data 
in muSTEM XTL, Vesta XYZ or Kirkland XYZ file format (see details below). In particular,
the output structure files in Kirkland XYZ format can be subsequently used as input
for the MsctKirkland.exe program for forward simulation of defocused EM images. 

This program has been written in C++ and can in principle be 
compiled under any OS supporting standard C++ compilers and execution environment. 
However, at present the executable module is only available for 64-bit Windows OS. Any 
Windows x64 PC can in principle be used for running this code.

The control of execution of Pdb.exe program is managed via an editable text file 
Pdb.txt which contains all modifiable input parameters for Pdb.exe. In 
order for the program to run, the input parameter file Pdb.txt must be present in the 
same folder where Pdb.exe is started from. Alternatively, the name of a suitable input
parameter file can be presented as the first command-line argument (following the pdb.exe
name in the command line) when starting the pdb.exe program. In this case, the parameter
file can be given an arbitrary name and can be located in any place accessible from the started
instance of pdb.exe. The format of these input parameter files is fixed and must
comply exactly with the structure described below.

The format of Pdb.txt file allows any number of "comment" lines to be present at the beginning and/or
at the end of the file, each such line starting with the double forward slash symbol, "//". The comment lines 
are ignored by Pdb.exe program. Other (non-comment) lines must all have the fixed structure 
(described below) and the fixed sequential order with respect to each other. Each non-comment line consists of:
(a) An arbitrary "expression" with no white spaces, which typically represents the "name" (brief description) of a 
parameter. The contents of this parameter name expression are ignored by Pdb.exe. 
An example of the parameter name is "3.Carbon_support_thickness_and_width_in_Angstroms:".
(b) One or more alphanumerical entries separated by a single white space from each other 
and from the parameter name, each such entry representing one parameter value, for example 
" 5 10", which corresponds to two parameter values equal to 5 and 10.
Note that the number of white spaces separating the parameters must be equal to one, 
otherwise a "zero" parameter may be wrongly read in. There should be also no further white 
spaces after the last parameter in the line. Each line should be terminated by the usual "end-
of-line / new line / caret-return", etc., symbols, as common for a given OS. The number and 
types of parameters in each parameter line of the file Pdb.txt is pre-determined and 
cannot be changed, only the values of these parameters can be edited. The parameter names 
and comment lines are supposed to give explicit hints about the type and number of parameters expected
in each line. 

Here is a typical example of a valid Pdb.txt parameter file.

=======================================================================================
//***PDB, Vesta and Kirkland XYZ file conversion parameters
1.Input_file_name: ASP_NEW_Kirck_DW.xyz
2.Input_file_type_(see_below): 2
3.Carbon_support_thickness_and_width_in_Angstroms: 0 0
4.CT_sample_cube_side_length_in_Angstroms: 10.0
5.Desired_rotation_angle_around_z_axis_in_degrees: 30.0
6.Desired_rotation_angle_around_y'_axis_in_degrees: 60.0
7.Desired_rotation_angle_around_z"_axis_in_degrees: -45.0
8.Centre_structure_in_the_miniball_(1)_or_do_not_shift_(0): 1
9.Maximum_distance_in_Angstroms_to_remove_duplicates: 0.5
10.Sort_0=no_sort_1=increasing_z_2=decreasing_occupancy: 0
11.Output_file_name: ASP_NEW_Vesta.xyz
12.Out_file_type_(see_below): 1
13.Free-form_1st_line_in_output_file: ASP_test

//Input parameter file for pdb.exe program
//1st parameter contains input file name (standard PDB file or XYZ file produced by Vesta).
//2nd parameter: 0 - PDB input file, 1 - for Vesta XYZ input file, 2 - for Kirkland XYZ file, 3 - text file with orientations/defocuses.
//3rd parameter contains optional carbon support layer thickness and transversal extent (x=y) in Angstroms, separated by a single white space
//4th parameter contains the desired CT sample cube side length in Angstroms (should be large enough to contain the whole sample inside for arbitrary 3D rotations)
//5th parameter contains the desired rotation angle around z axis in degrees
//6th parameter contains the desired rotation angle around y' axis in degrees (y' is the y axis after the initial rotation around the z axis)
//7th parameter contains the desired rotation angle around z" axis in degrees (z" is the z axis after the initial rotations around the z and y' axes)
//8th parameter: 1 = shift the structure so that to centre it in the determined minimal containing ball, 0 = do not shift the structure
//9th parameter: maximum distance (in Angstroms) between two atoms to be considered "duplicates" and remove all but one; if this distance is negative, no action will be taken
//10th parameter: sort output: 0=no_sort, 1=sort_by_z_coordinates_ascending 2=sort_by_occupancy_column_values_descending (sorting according to atomic numbers is always done)
//11th parameter contains output file name
//12th parameter: 0 - muSTEM-format XTL file, 1 - Vesta XYZ file, 2 -  Kirkland-format XYZ file
//13th parameter contains the free-form info line that will be recorded into the first line of the output file

// This program reads a PDB, Vesta XYZ, Kirkland XYZ file or a text file with orientations/defocuses, 
// centers the position of the "molecule" within the "enclosing cube" [0, ctblength] x [0, ctblength] x [0, ctblength],
// optionally adds a layer of amorphous carbon immediately upstream of the "molecule" along the z coordinate,
// optionally rotates the "molecule" (together with the carbon substrate) by the given "intrinsic" Euler angles around the z, y' and z" axes, and
// optionally removes "duplicate" atoms that are located closer than the given distance to other atoms in the "molecule",
// optionally sorts all atoms in ascending order with respect to z coordinate or with respect to decreasing occupancy parameter, and
// outputs the data in the form of a muSTEM, Vesta XYZ or Kirkland XYZ file.

========================================================================================

Note that the "numbering" of parameters in the parameter names, such as "3" in 
"3.Carbon_support_thickness_and_width_in_Angstroms: 10 10", is inconsequential and is ignored by 
Pdb.exe. On the other hand, the order of the parameter lines in this file is very 
important and must be the same as in the above example. In other words, Pdb.exe 
interprets the input parameters according to the order of the non-comment lines in the 
parameter file, while ignoring any numbering that may or may not appear in the parameter 
names.

In the above example, Parameter 1 contains the name of an "atomic structure" file in 
PDB (https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format), Vesta XYZ (https://jp-minerals.org/vesta/en/),
Kirkland XYZ (E. J. Kirkland, Advanced Computing in Electron Microscopy, second edition, Springer, New York, 2010) format 
(see also https://github.com/jhgorse/kirkland/blob/master/temsim/slicelib.cpp),
or a text file with orientations/defocuses. In the latter case (text) file, the filename must have either ".txt" or
".RELION" filename extension, and the corresponding file must be in the format used in MsctKirkland and PhaseRetrieval programs
for storing image orientations and defocuses. This file may have one
of two distinct formats which is determined by the file name extenstion. 
(a) If the file name extension is ".txt", the file may contain an
arbitrary number of leading comment lines, each starting with a pair of forward slashes "//". 
Non-comment lines all have the same structure, but may contain different number of entries, 
according to the following scheme. Columns one and two always contain the rotation angles (in degrees)
around the Z and Y' axes, respectively (where Y' axis is the Y axis after the initial rotation
of the 3D space around the Z axis). These two angles are followed by an arbitrary number of pairs
of values, with each odd column (starting from the third one) containing the angle (in degrees)
of rotation around the Z'' axis (Z" axis is the Z axis after the initial rotations
of the 3D space around the Z and Y' axes, which corresponds to the optic axis = illumination direction). 
Each even column (starting from the fourth one) contains a defocus distance in 
angstroms along the Z'' axis. The number of lines in this file is not limited.
An example of a valid "rotation angles and defocus distances" file is given below. 
(b) If the file extension is ".relionnew", then this file contains an arbitrary number of lines,
each line containing exactly ten entries, separated by white spaces, with the following contents:
1st entry contains the sequential number of the line - this parameter is not used
2nd entry contains the name of the file where the data was sourced from - this parameter is not used
3rd entry is the image shift along x in angstroms 
4th entry is the image shift along y in angstroms 
5th entry is the defocus distance dx corresponding to the x coordinate
6th entry is the defocus distance dy corresponding to the y coordinate
7th entry contains the astigmatism angle (alpha) in degrees
8th entry is the rotation angle "rot" around the Z coordinate in degrees
9th entry is the rotation angle "tilt" around the Y' coordinate in degrees
10th entry is the rotation angle "psi" around the Z" coordinate in degrees
Examples of defocus/orientation files can be found in ReadmePhaseRetrieval.txt and ReadmeMsctKirkland.txt files.
The input file must be present in the same folder where Pdb.exe is started from, or, 
alternatively, the filename can include a fully specified pathname (OS specific). 
The ".txt" or ".relionnew" orientation files are usually used in this program in order to create a graphical
representation of the distribution of illumination directions on the unit sphere in 3D. This is achieved in this program
by creating a "fake" atomic structure containing only H atoms located on the surface of a 3D sphere at the points, 
corresponding to illumination directions present in the input orientation file. If the output is produced
in the form of a Vesta XYZ file (Parameter 2 set to 1), the output file with the fake H atom structure can be loaded
and viewed in Vesta. This can be useful, in particular, in order to check if there are any particular missing
regions of illumination directions in the input orientation file.

Parameter 2 can be equal to 0, 1, 2 or 3. It defines the input file format.
"0" corresponds to PDB input file format. 
"1" corresponds to Vesta XYZ input file format.
"2" corresponds to Kirkland XYZ input file format.
"3" corresponds to a text file with orientations/defocuses.

Parameter 3 contains the thickness (z dimension) and width (x and y dimensions) of the optional
amorphous carbon support layer in Angstroms. This layer is created if the thickness parameter is larger than zero.
The layer is created immediately below (along z) the atomic structure read from the input PDB or XYZ file and
is centred transversally with respect to this structure. The positions of the carbon atoms in this layer
are random and will be different for different runs of this program, even if all the input parameters are 
kept the same. The number of carbon atoms can also be slightly different between the runs. The implementation
here is borrowed from EJ Kirkland's pdb2xyz.cpp program (https://sourceforge.net/projects/computem/files/). It
uses the following parameters, in particular: CDENSITY = 2.0 (density in gm/cm^3 approx. for amorphous C),
CFILL_FRACTION = 0.9; (filling fraction for this density) and RMIN = 1.4 (min separation distance of C (in Angstroms)).
Note that amorphous ice can be also added to molecules, but this happens in subsequent simulation of defocused images
in the MsctKirkland program. Unlike the carbon support for nanoparticles, the ice layer in cryo-EM is different for
each instance of the molecule within a single dataset. Therefore, a different pseudo-random ice layer is added at
each illumination angle in MsctKirkland, rather than a single ice layer associated with the input XYZ structure.

Parameter 4 contains the side length (in Angstroms) of the virtual cube in 3D that is supposed to contain the atomic
structure inside (see details below).

Parameter 5 defines the rotation angle around the z axis in degrees. The rotation of the whole structure (including the
optionally added carbon support layer) is carried out after the centering of the structure inside the enclosing cube.
The Euler angle convention corresponds to that in Heymann et al, Journal of Structural Biology 151 (2005) 196–207.

Parameter 6 defines the rotation angle around the y' axis in degrees (y' axis corresponds to the y axis after 
the initial rotation around the z axis). The rotation of the whole structure (including the
optionally added carbon support layer) is carried out after the centering of the structure and after the initial
rotation around the z axis.

Parameter 7 defines the rotation angle around the z" axis in degrees (z" axis corresponds to the z axis after 
the initial rotations around the z and y' axes). The rotation of the whole structure (including the
optionally added carbon support layer) is carried out after the centering of the structure and after the initial
rotation around the z and y' axes.

Parameter 8 determines if the structure is shifted so that it would be centred in the determined minimal containing
ball (which contains the structure for any 3D rotation with the rotation centre at the centre of the miniball). If 
parameter 8 is equal to 1, the structure may be shifted so that it is centred in the minimal containing ball, and
the size of the containing cube may be increased as needed to contain the miniball. If parameter 8 is equal to 0, the 
structure is not shifted and, it the containing cube is not large enough to enclose the minimall ball, an error is 
displayed and the program is interrupted.

Parameter 9 defines the maximum distance (in Angstroms) between two atoms to be considered "duplicates". Such duplicates
are removed from the structure before it is written to the output file (only one atom is retained from each group of
atoms located closer than this distance to a given atom). If this distance is negative, no action will be taken.

Parameter 10 can be equal to 0, 1 or 2.
"0" corresponds to sorting of atoms in the output file according to the decreasing atomic numbers.
"1" corresponds to sorting of atoms in the output file according to the increasing z coordinates.
"2" corresponds to sorting of atoms in the output file according to the decreasing values in the occupancy column.

Parameter 11 contains the name of the output file.

Parameter 12 can be equal to 0, 1 or 2.
"0" corresponds to muSTEM XTL output file format (http://tcmp.ph.unimelb.edu.au/mustem/dist/muSTEM_v5.0_manual.pdf). 
"1" corresponds to Vesta XYZ file output format (https://jp-minerals.org/vesta/en/).
"2" corresponds to Kirkland XYZ output file format (E. J. Kirkland, Advanced Computing in Electron Microscopy,
second edition, Springer, New York, 2010) (see also https://github.com/jhgorse/kirkland/blob/master/temsim/slicelib.cpp).
Note that the optional column with RMS amplitudes of atomic vibrations is not written out. It can be added later
e.g. by editing the output file in Excel.

Parameter 13 contains a free-form information line (with no white spaces) that is written in the muSTEM, Vesta and Kirkland XYZ file formats.

================================================================================================== 
Example of a valid XYZ file in Kirkland's XYZ format (see detailed description in the 
references given above):

Aspartate_molecule
10.000000 10.000000 10.000000
7 4.386000 6.358500 5.031500 1.000000 0.1
6 4.233000 4.956500 4.621500 1.000000 0.1
6 2.835000 4.490500 4.936500 1.000000 0.1
8 2.169000 5.085500 5.751500 1.000000 0.1
6 5.242000 4.090500 5.378500 1.000000 0.1
6 6.641000 4.475500 4.969500 1.000000 0.1
8 6.812000 5.351500 4.155500 1.000000 0.1
8 7.695000 3.844500 5.508500 1.000000 0.1
8 2.329000 3.414500 4.313500 1.000000 0.1
1 3.775000 6.959500 4.498500 1.000000 0.1
1 4.225000 6.465500 6.021500 1.000000 0.1
1 4.411000 4.869500 3.549500 1.000000 0.1
1 5.122000 4.245500 6.450500 1.000000 0.1
1 5.070000 3.040500 5.141500 1.000000 0.1
1 8.572000 4.125500 5.215500 1.000000 0.1
1 1.428000 3.153500 4.549500 1.000000 0.1
-1

=========================================================================================================
Example of a valid XYZ file in Vesta's XYZ format (see detailed description in the 
references given above):

16
Aspartate_molecule
O 2.169000 6.446000 5.869000 0.000000
O 6.812000 6.712000 4.273000 0.000000
O 7.695000 5.205000 5.626000 0.000000
O 2.329000 4.775000 4.431000 0.000000
N 4.386000 7.719000 5.149000 0.000000
C 4.233000 6.317000 4.739000 0.000000
C 2.835000 5.851000 5.054000 0.000000
C 5.242000 5.451000 5.496000 0.000000
C 6.641000 5.836000 5.087000 0.000000
H 3.775000 8.320000 4.616000 0.000000
H 4.225000 7.826000 6.139000 0.000000
H 4.411000 6.230000 3.667000 0.000000
H 5.122000 5.606000 6.568000 0.000000
H 5.070000 4.401000 5.259000 0.000000
H 8.572000 5.486000 5.333000 0.000000
H 1.428000 4.514000 4.667000 0.000000

===========================================================================================================
Example of a valid file in muSTEM's XTL format (see detailed description in the 
references given above):

Aspartate_molecule
10.000000 10.000000 10.000000 90.0 90.0 90.0
200.000000 
4 
O
4 8.000000 1.0 0.0
0.216900 0.644600 0.586900
0.681200 0.671200 0.427300
0.769500 0.520500 0.562600
0.232900 0.477500 0.443100
N
1 7.000000 1.0 0.0
0.438600 0.771900 0.514900
C
4 6.000000 1.0 0.0
0.423300 0.631700 0.473900
0.283500 0.585100 0.505400
0.524200 0.545100 0.549600
0.664100 0.583600 0.508700
H
7 1.000000 1.0 0.0
0.377500 0.832000 0.461600
0.422500 0.782600 0.613900
0.441100 0.623000 0.366700
0.512200 0.560600 0.656800
0.507000 0.440100 0.525900
0.857200 0.548600 0.533300
0.142800 0.451400 0.466700

============================================================================================================