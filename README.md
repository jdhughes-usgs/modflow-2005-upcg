modflow-2005-upcg
=================

MODFLOW-2005 version 1.10 with the UPCG solver

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

o Version 1.1 07/09/2013:

Updated MODFLOW source code to version 1.10 and modified the Visual Studio solution and upcgc C++ project to link against version 5.0 of the 64-bit CUDA toolkit. 

Modified matrix assembly in UPCG7AP subroutine (upcg7.f) to work with more general MODFLOW problems with wetting and drying, interior inactive cells, and constant head cells.

Added relaxation for MILU0 preconditioner and steady-state and transient dampening to all preconditioners.


DISCLAIMER and NOTICE

Please refer to the USGS Software User Rights Notice (http://water.usgs.gov/software/help/notice/) for complete use, copyright, and distribution information. Although this software has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the software program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.  Although the software has been tested, there could be undetected errors. Users are encouraged to report any errors to these authors.  Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

###Build Order

A makefile is not available for MODFLOW-2005-UPCG source files. However, the build order for the source files is listed below.

####Compile and build gmg.lib
  ```
  1.  solvers.c                                               
  2.  r_vector.c                                              
  3.  mf2kgmg.c                                               
  4.  ccfd.c
  ```
                                                          
####Compile and build upcgc.lib                               
  ```
  1.  cuda_kernels.cu                                         
  2.  upcgc.cpp
  ```
                                                          
####Compile MODFLOW-2005                                      
  ```
  1.  gwf2bas7.f                                              
  2.  sip7.f                                                  
  3.  gwf2swt7.f                                              
  4.  gwf2rch7.f                                              
  5.  gwf2chd7.f                                              
  6.  gwf2ibs7.f                                              
  7.  gwf2sub7.f                                              
  8.  gwf2fhb7.f                                              
  9.  de47.f                                                  
  10. gwf2bcf7.f                                              
  11. pcgn2.f90                                               
  12. gwf2str7.f                                              
  13. gwf2evt7.f                                              
  14. gwf2huf7.f                                              
  15. pcg7.f                                                  
  16. gwf2ets7.f                                              
  17. gmg7.f                                                  
  18. gwf2wel7.f                                              
  19. gwf2hfb7.f                                              
  20. gwf2riv7.f                                              
  21. gwf2drt7.f                                              
  22. obs2bas7.f                                              
  23. gwf2lpf7.f                                              
  24. upcg7.f                                                 
  25. gwf2ghb7.f                                              
  26. obs2str7.f                                              
  27. gwf2res7.f                                              
  28. gwf2drn7.f                                              
  29. gwf2lak7.f                                              
  30. gwf2gag7.f                                              
  31. obs2riv7.f                                              
  32. obs2ghb7.f                                              
  33. gwf2hydmod7.f                                           
  34. obs2drn7.f                                              
  35. gwf2mnw27.f                                             
  36. obs2chd7.f                                              
  37. gwf2mnw17.f                                             
  38. lmt7.f                                                  
  39. gwf2mnw2i7.f                                            
  40. hufutl7.f                                               
  41. upcg7polyu.f                                            
  42. upcg7polya.f                                            
  43. gwf2sfr7.f                                              
  44. upcg7lanczos.f                                          
  45. parutl7.f                                               
  46. mf2005-GPU.f                                            
  47. gwf2uzf1.f   
  ```
                                                          
####Link MODFLOW-2005 objects including  gmg.lib and upcgc.lib

