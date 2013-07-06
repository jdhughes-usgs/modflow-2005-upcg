modflow-2005-upcg
=================

MODFLOW-2005 version 1.8 with the UPCG solver

The UPCG solver includes support for parallel solution of MODFLOW on 
(1) multi-core CPUs using OpenMP and (2) general purpose graphical
processing units (GPGPUs) using the NVIDIA CUBLAS Library.


Documentation:

Hughes, J.D., and White, J.T., 2012, Use of general purpose graphics processing units with MODFLOW: Ground Water, doi: 10.1111/gwat.12004.


GENERAL USE:

Source code for MODFLOW--2005 with the UPCG solver is available in the code\src subdirectory. A Microsoft Visual Studio solution and associated fortran, C++, and cu projects are available in the code subdirectory.

Development and compilation of CUDA and CUBLAS functions requires installation of the CUDA toolkit which includes the necessary CUDA drivers. NVIDIA also provides GPGPU code samples in the CUDA software development kit (SDK). Current versions of each of these are available free of charge from NVIDIA at http://www.nvidia.com/content/cuda/cuda-downloads.html. A 64-bit version of the MODFLOW-2005 executable with the UPCG solver and the linked OpenMP dynamic-link library (DLL) for the Windows 7 operating system and the Tesla C2050 GPGPU is available in the bin subdirectory and only requires installation of the CUDA drivers, which are also available as a stand-alone installation from NVIDIA, and a GPGPU with NVIDIA Compute Capability 2.0 or greater. This executable was linked against version 4.1.28 of the 64-bit CUDA toolkit.

Contact Joseph D. Hughes < jdhughes (at) usgs [dot] gov >


Version History:

o Version 1.0 12/12/2012:

This version is the initial release. 

o Version 1.1 07/06/2013:

Updated MODFLOW source code to version 1.10 and modified the Visual Studio solution and upcgc C++ project to link against version 5.0 of the 64-bit CUDA toolkit. 

Modified matrix assembly in UPCG7AP subroutine (upcg7.f) to work with more general MODFLOW problems with wetting and drying, interior inactive cells, and constant head cells.


DISCLAIMER and NOTICE

Please refer to the USGS Software User Rights Notice (http://water.usgs.gov/software/help/notice/) for complete use, copyright, and distribution information. Although this software has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the software program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.  Although the software has been tested, there could be undetected errors. Users are encouraged to report any errors to these authors.  Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

