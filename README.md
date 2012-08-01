modflow-2005-upcg
=================

MODFLOW-2005 version 1.8 with the UPCG solver

The UPCG solver includes support for parallel solution of MODFLOW on 
(1) multi-core CPUs using OpenMP and (2) general purpose graphical
processing units (GPGPUs) using the NVIDIA CUBLAS Library.


Documentation:

Hughes, J.D., and White, J.T., in review, Use of general purpose graphical processing units with MODFLOW-2005: submitted to Ground Water


GENERAL USE:

Development and compilation of CUDA and CUBLAS functions requires installation of the CUDA toolkit which includes the necessary CUDA drivers. NVIDIA also provides GPGPU code samples in the CUDA software development kit (SDK). Current versions of each of these are available free of charge from NVIDIA at http://www.nvidia.com/content/cuda/cuda-downloads.html. A 64-bit version of the MODFLOW-2005 executable with the UPCG solver and the linked OpenMP dynamic-link library (DLL) for the Windows 7 operating system and the Tesla C2050 GPGPU is available in the bin subdirectory and only requires installation of the CUDA drivers, which are also available as a stand-alone installation from NVIDIA, and a GPGPU with NVIDIA Compute Capability 2.0 or greater. This executable was linked against version 4.1.28 of the 64-bit CUDA toolkit.

Contact Joseph D. Hughes < jdhughes (at) usgs [dot] gov >


DISCLAIMER and NOTICE

Please refer to the USGS Software User Rights Notice (http://water.usgs.gov/software/help/notice/) for complete use, copyright, and distribution information. The USGS provides no warranty, expressed or implied, as to the correctness of the furnished software or the suitability for any purpose. The software has been tested, but as with any complex software, there could be undetected errors. Users who find errors are requested to report them to the USGS.

References to non-USGS products, trade names, and (or) services are provided for information purposes only and do not constitute endorsement or warranty, express or implied, by the USGS, U.S. Department of Interior, or U.S. Government, as to their suitability, content, usefulness, functioning, completeness, or accuracy.

