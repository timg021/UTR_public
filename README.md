# Unified Tomographic Reconstruction solution
A collection of C++ source code modules for forward and inverse simulations of Unified Tomographic Reconstruction.

UPDATE 14 October 2022. Multiple updates, improvements and bug fixes applied to the previous version of the code. Both the Windows x64 and the Ubuntu Linux versions have been recompiled and tested. Main addition for this update is the GDCM library (https://github.com/malaterre/GDCM), which is used here for optionally saving output files in standard DICOM format. See details of the updated usage in the Readme*.txt files that can be found in each of the four main projects, MsctKirkland, PDB, PDB-Compare and PhaseRetrieval.

UPDATE 14 July 2022. Multiple updates, improvements and bug fixes applied to the previous version of the code. The arXiv paper describing this method has been also updated - see https://arxiv.org/abs/2206.09151 ([v4] Fri, 8 Jul 2022 03:06:41 UTC). Both the Windows x64 and the Ubuntu Linux versions have been recompiled and tested.

UPDATE 3 July 2022. Multiple updates and bug fixes have been applied. The main change, compared to the previous version, is related to the use of mutlithreading in FFTW library. The FFTW libraries have been compiled locally in Windows 10, in Linux we still use the standard installed FFTW libraries (version 3.3.5). Statically linked FFTW libraries are now used instead of DLLs in most places. The internal multithreading is still based on OMP, but has been re-modelled, resulting in imroved (x2-x3) performance. Proper (quantitatively correct) normalization and noise filtering have been implemented in UTR.

The code has been compiled and tested under Windows x64 OS (see the included Microsoft Visual Studio 2022 solution configuration files) and under Ubuntu Linux 20.04.4 LTS using gcc 9.4.0 compiler (see the makefiles in Linux_ subfolders of MultisliceCppTest/PhaseRetrieval).

The theory behind the UTR algorithm can be found in the paper: TIMUR E. GUREYEV, HAMISH G. BROWN, HARRY M. QUINEY, AND LESLIE J. ALLEN, "Fast unified reconstruction algorithm for conventional, phase-contrast and diffraction tomography", arXiv.org (2022).

This repository contains four main projects: 
(1) UTR (implemented in the UTR module) and DHT/CHR/vCTF programs (the last three methods are implemented in the DHT module) - "inverse" 3D algorithms for reconstruction of the spatial distribution of either electrostatic potential (in the case of TEM imaging) or the imaginary part beta of the complex refractive index (in the case of hard X-ray imaging) from input defocused images collected at different illumination directions (i.e. at different spatial orientation of the imaged object) in the usual manner for tomographic (CT) imaging with optional free-space propagation between the sample and the 2D detector. See usage instructions in the filese PhaseRetrieval/ReadmeUTR.txt and PhaseRetrieval/ReadmePhaseRetrieval.txt files.
(2) MsctKirkland - "forward" multislice image simulations using atomic XYZ files as input; this code is largely based on E.J. Kirkland's temsim code which is described in the following references [https://github.com/jhgorse/kirkland/tree/master/temsim. See also https://sourceforge.net/projects/computem/. 
Description of the code and the underlying theory can be found in E. J. Kirkland, Advanced  Computing in Electron Microscopy, second edition, Springer, New York, 2010.] 
See MsctKirkland/ReadmeUTR.txt file for usage instructions.
(3) pdb - simple programs for manipulation of XYZ files. See pdb/ReadmePdb.txt for usage instructions.
(4) pdb-compare - simple program for comparing pairs of XYZ files. See pdb-compare/ReadmePdb-compare.txt for usage instructions.

Details about FFTW libraries and how to compile them under different OS can found here: https://www.fftw.org/. The user may prefer to compile their own versions of these libraries. Note that both single-precision (float) and double-precision (double = default) 2D and 3D FFT routines in the multi-threaded form may be required.

The folders XArrayLibrary and TemsimLibrary contain auxilliary modules used for compilation of the 4 main programs listed above. XArrayLibrary mostly contains general-purpose C++ templates implementing 1D, 2D and 3D matrices and mathematical operations on them. TemsimLibrary contains various modules used by MsctKirlkand, pdb and pdb-compare programs. See individual modules for details.
