11/25/2012

The python script RunGPUTestProblem.py in the python subdirectory can be used
to create and run the model simulations evaluated in Hughes and White (2012).
The script is standard python with NUMPY except for the functionality included
in the SimpleMF.py script, which is based on standard python. 

The script was tested using python 2.6 but should work without modification if you
are using python 2.7. The script is likely to need modification for python 3.x.

To run the test problems:

  1) open a command line in the python subdirectory
  2) type "python RunGPUTestProblem.py -help" to get a list of the required and
     optional command line options 
  3) type "RunGPUTestProblem.py" with the command line options you would like to 
     evaluate.

An example of the command line options to run a 500 x 500 cell model with 1 to 3 
model layers using the MODFLOW-2005 PCG2 solver and the UPCG solver with the 
MILU0 and GLSPOLY preconditioners in serial (-cpu) and parallel (-openmp) on 
the CPU is:


    python RunGPUTestProblem.py -nx 500 -minLay 1 -maxlay 3 -pcg -upcg -cpu 
           -openmp -milu0 -glspoly


The list file for each model run is renamed and saved in the LISTArchive 
subdirectory within the subdirectory created for the model simulation. The binary
heads file for each model run can also be renamed and saved in the HEADArchive 
subdirectory within the subdirectory created for the model simulation by using
the optional -saveheads command line argument. In the example listed above the 
input and output files are created in a 00500 subdirectory that is at the same 
level as the python subdirectory containing the RunGPUTestProblem.py script. 

The directory structure that will exist after running the 200 x 200, 500 x 500,
1000 x 1000, 2000 x 2000, and 4000 x 4000 cell models for combinations of solvers,
UPCG preconditioners, and UPCG hardware solution options and saving the binary 
heads files (-saveheads) is:

    TestProblem\
        python\
        ref\
        00200\
            HEADArchive\
            LISTArchive\
        00500\
            HEADArchive\
            LISTArchive\
        01000\
            HEADArchive\
            LISTArchive\
        02000\
            HEADArchive\
            LISTArchive\
        03000\
            HEADArchive\
            LISTArchive\
            

Reference:
  Hughes, J.D. and White, J.T., 2012. Use of general purpose graphics processing
  units with MODFLOW: Ground Water, doi: 10.1111/gwat.12004.