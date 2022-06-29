# Unified Tomographic Reconstruction solution
A collection of C++ source code modules for forward and inverse simulations of Unified Tomographic Reconstruction.

29 June 2022. PLEASE, NOTE THAT A NEW VERSION OF THE CODE IN THIS REPOSITORY WILL BE UPLOADED IN THE NEXT FEW DAYS. IT WILL CONTAIN A NUMBER OF IMPORTANT BUG FIXES AND UPDATES. IF YOU ARE ALREADY USING THE INITIAL VERSION OF THE CODE, IT MAY BE WORTH UPDATING TO THE NEW VERSION WHEN IT ARRIVES.

The code has been compiled and tested under Windows x64 OS (see the included Microsoft Visual Studio 2022 solution configuration files) and under Ubuntu Linux 20.04.4 LTS using gcc 9.4.0 compiler (see the makefiles in Linux_ subfolders of MultisliceCppTest/PhaseRetrieval).

The theory behind the UTR algorithm can be found in the paper: TIMUR E. GUREYEV, HAMISH G. BROWN, HARRY M. QUINEY, AND LESLIE J. ALLEN, "Fast unified reconstruction algorithm for conventional, phase-contrast and diffraction tomography", arXiv.org (2022).

This repository contains four main projects: 
(1) UTR (implemented in the UTR module) and DHT/CHR/vCTF programs (the last three methods are implemented in the DHT module) - "inverse" 3D algorithms for reconstruction of the spatial distribution of either electrostatic potential (in the case of TEM imaging) or the imaginary part beta of the complex refractive index (in the case of hard X-ray imaging) from input defocused images collected at different illumination directions (i.e. at different spatial orientation of the imaged object) in the usual manner for tomographic (CT) imaging with optional free-space propagation between the sample and the 2D detector. See usage instructions in the filese PhaseRetrieval/ReadmeUTR.txt and PhaseRetrieval/ReadmePhaseRetrieval.txt files.
(2) MsctKirkland - "forward" multislice image simulations using atomic XYZ files as input; this code is largely based on E.J. Kirkland's temsim code which is described in the following references [https://github.com/jhgorse/kirkland/tree/master/temsim. See also https://sourceforge.net/projects/computem/. 
Description of the code and the underlying theory can be found in E. J. Kirkland, Advanced  Computing in Electron Microscopy, second edition, Springer, New York, 2010.] 
See MsctKirkland/ReadmeUTR.txt file for usage instructions.
(3) pdb - simple programs for manipulation of XYZ files. See pdb/ReadmePdb.txt for usage instructions.
(4) pdb-compare - simple program for comparing pairs of XYZ files. See pdb-compare/ReadmePdb-compare.txt for usage instructions.

The folder fftwlib contains FFTW (Fast Fourier Transform) modules used by the main programs listed above. Details about FFTW libraries and how to compile them under different OS can found here: https://www.fftw.org/. The user may prefer to compile their own versions of these libraries. Note that both single-precision (float) and double-precision (double = default) 2D and 3D FFT routines may be required.

The folders XArrayLibrary and TemsimLibrary contain auxilliary modules used for compilation of the 4 main programs listed above. XArrayLibrary mostly contains general-purpose C++ templates implementing 1D, 2D and 3D matrices and mathematical operations on them. TemsimLibrary contains various modules used by MsctKirlkand, pdb and pdb-compare programs. See individual modules for details.
