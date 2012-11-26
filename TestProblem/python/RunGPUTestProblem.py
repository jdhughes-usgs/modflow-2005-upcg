import sys
import os
import subprocess
import shutil
import numpy as np
import SimpleMF as smf

def TestDirExist(ctest):
    #--Make sure output directories are created
    for f in ctest:
        print 'evaluating... "{0}"'.format( f )
        fa = os.path.abspath(f)
        d = os.path.dirname(fa)
        if not os.path.exists(d):
            print 'creating directory path...\n  "{0}"'.format( d )
            os.makedirs(d)


command_args = [['-nx value','required','number of rows and columns'] , \
                ['-minlay value','required','minimum number of layers to evaluate'], \
                ['-maxlay value','required','maximum number of layers to evaluate'], \
                ['-saveheads','optional','save the binary head file for each run'], \
                ['-pcg','optional','use the MODFLOW-2005 PCG2 solver'], \
                ['-gmg','optional','use the MODFLOW-2005 GMG solver'], \
                ['-upcg','optional','use the UPCG solver'], \
                ['-cpu','optional','serial UPCG solution on the CPU'], \
                ['-openmp','optional','parallel UPCG solution on the CPU'], \
                ['-gpu','optional','parallel UPCG solution on the GPU'], \
                ['-jacobi','optional','use UPCG Jacobi preconditioner'], \
                ['-milu0','optional','use UPCG MILU0 preconditioner'], \
                ['-glspoly','optional','use GLSPOLY Jacobi preconditioner'], \
                ['-help','optional','get list of available command line arguments']]
preconditioner_dict = { 'jacobi': 1, 'milu0': 3, 'glspoly': 4 }
hardware_dict = {'cpu':1, 'openmp':2, 'gpu':3 }

exe_name = os.path.join( '..','..','bin','modflow-2005-UPCG_x64.exe' )

#--base options
saveHeads = False
nx = ny = None #500
minLay = None  #1
maxLay = None  #3
#--solver options
isPCG = False
isGMG = False
isUPCG = False
#--preconditioner flags
isJacobi = False
isMILU0 = False
isGLSPOLY = False
#--hardware flags
isCPU = False
isOpenMP = False
isGPU = False
#--read command line options
narg = len(sys.argv)
iarg = 0
if narg > 1:
    while iarg < narg-1:
        iarg += 1
        basearg = sys.argv[iarg].lower()
        if basearg == '-nx' or basearg == '-ny':
            try:
                iarg += 1
                nx = int( sys.argv[iarg] )
                ny = nx
                print 'command line arg: nx = ', nx
                print 'command line arg: ny = ', ny
            except:
                print 'cannot parse command line arg: nx or ny'
        elif basearg == '-minlay':
            try:
                iarg += 1
                minLay = int( sys.argv[iarg] )
                print 'command line arg: minLay = ', minLay
            except:
                print 'cannot parse command line arg: minLay'
        elif basearg == '-maxlay':
            try:
                iarg += 1
                maxLay = int( sys.argv[iarg] )
                print 'command line arg: maxLay = ', maxLay
            except:
                print 'cannot parse command line arg: maxLay'
        elif basearg == '-saveheads':
            saveHeads = True
            print 'command line arg: -saveheads'
        #--solver flags
        elif basearg == '-upcg':
            isUPCG = True
        elif basearg == '-pcg':
            isPCG = True
        elif basearg == '-gmg':
            isGMG = True
        #--hardware flags
        elif basearg == '-cpu':
            isCPU = True
        elif basearg == '-openmp':
            isOpenMP = True
        elif basearg == '-gpu':
            isGPU = True
        #--preconditioner flags
        elif basearg == '-jacobi':
            isJacobi = True
        elif basearg == '-milu0':
            isMILU0 = True
        elif basearg == '-glspoly':
            isGLSPOLY = True
        elif 'help' in basearg:
            print 'Available command line arguments'
            for t in command_args:
                print '  {0:14s} : [{1:8s}] {2}'.format( t[0],t[1],t[2] )
            print '\n    value : specified value for a required argument' 
            sys.exit( 1 )
#--error checking
#--make sure the problem dimensions are specified
if nx == None:
    print '\nERROR\n  -nx value : command line argument not specified'
    print '     value is the number of rows and columns in the simulation'
    sys.exit( 0 )
if nx != 200 and nx != 500 and nx != 1000 and nx != 2000 and nx != 4000:
    print '\nERROR\n  Script only evaluated for -nx values of 200, 500, 1000, 2000, and 4000'
    sys.exit( 0 )
if minLay == None:
    print '\nERROR\n  -minLay value : command line argument not specified'
    print '     value is the minimum number of layers in the simulation'
    sys.exit( 0 )
if minLay < 1:
    print '\nERROR\n  the minimum number of layers (-minLay value) is 1'
    print '  Specified minLay value {0}'.format( minLay )
    sys.exit( 0 )
if maxLay == None:
    print '\nERROR\n  -maxLay value : command line argument not specified'
    print '     value is the maximum number of layers in the simulation'
    sys.exit( 0 )
if maxLay < 1 or maxLay < minLay:
    print '\nERROR\n  the maximum number of layers (-maxLay value) must be >= 1\n  and >= the minimum number of layers (-minLay)'
    print '  Specified minLay value {0}'.format( minLay )
    print '  Specified maxLay value {0}'.format( maxLay )
    sys.exit( 0 )
#--write message if -saveheads command line argument is specified
if saveHeads == True:
    print '\nBinary heads file will be saved for each simulation'
    print '  in the {0} subdirectory\n'.format( os.path.join( '..','{0:05d}'.format( nx ),'HEADArchive' ) )
#--make sure at least one solver is specified
if isUPCG == False and isPCG == False and isGMG == False:
    print '\nERROR -- nothing to do'
    print '  one or more solver keywords should be specified'
    print '    -upcg -- use the UPCG solver'
    print '    -pcg  -- use the MODFLOW-2005 PCG2 solver'
    print '    -gmg  -- use the MODFLOW-2005 GMG solver'
    sys.exit( 0 )
else:
    if isPCG == True:
        print 'using the MODFLOW-2005 PCG2 solver'
    if isGMG == True:
        print 'using the MODFLOW-2005 GMG solver'
    if isUPCG == True:
        print 'using the UPCG solver'
#--build preconditioner list
hardware_list = []
if isUPCG == True:
    if isCPU == True:
        hardware_list.append( hardware_dict['cpu'] )
    if isOpenMP == True:
        hardware_list.append( hardware_dict['openmp'] )
    if isGPU == True:
        hardware_list.append( hardware_dict['gpu'] )
    if len( hardware_list ) > 0:
        print 'UPCG hardware options to be evaluated'
        for idx in hardware_list:
            for key, value in hardware_dict.iteritems():
                if value == idx:
                    print '  {0} preconditioner'.format( key )
                    break
    else:
        print 'No UPGC hardware option specified'
        print '  possible UPGC hardware options'
        print '    -cpu    -- serial UPCG solution on the CPU'
        print '    -openmp -- parallel UPCG solution on the CPU'
        print '    -gpu    -- parallel UPCG solution on the GPU'
        sys.exit( 0 )
    #--build preconditioner list
    pc_list = []
    if isJacobi == True:
        pc_list.append( preconditioner_dict['jacobi'] )
    if isMILU0 == True:
        pc_list.append( preconditioner_dict['milu0'] )
    if isGLSPOLY == True:
        pc_list.append( preconditioner_dict['glspoly'] )
    if len( pc_list ) > 0:
        print 'UPCG preconditioners to be evaluated'
        for idx in pc_list:
            for key, value in preconditioner_dict.iteritems():
                if value == idx:
                    print '  {0} preconditioner'.format( key )
                    break
    else:
        print 'No UPGC preconditioner specified'
        print '  possible UPGC preconditioner options'
        print '    -jacobi  -- use Jacobi preconditioner'
        print '    -milu0   -- use MILU0 preconditioner'
        print '    -glspoly -- use GLSPOLY preconditioner'
        sys.exit( 0 )
#--make sure directory exists for this simulation
testname1 = os.path.join( '..','{0:05d}'.format(nx),'test.dat')
testname2 = os.path.join('..','{0:05d}'.format(nx),'LISTArchive','test.dat')
testname3 = os.path.join('..','{0:05d}'.format(nx),'HEADArchive','test.dat')
TestDirExist([testname1,testname2])
if saveHeads == True:
    TestDirExist([testname3])
#--get hydraulic conductivity for base 1000 x 1000 and current runs
k_ref = os.path.join( '..','ref','Kh_1000.ref' )
k0 = smf.loadArrayFromFile(1000,1000,k_ref)
kavg = k0.mean()
k0 = []
k_file = os.path.join( '..','ref','Kh_{0:04d}.ref'.format( nx ) )
print 'reading k data from...', k_file
kb = smf.loadArrayFromFile(ny,nx,k_file)
#--calculate grid dimensions
xlen = ylen = 1000.0
dx = xlen / float( nx )
dy = ylen / float( ny )
#constant data for all simulations
nper = 1
perlen = [100.0]
nstp = [1]
tsmult = [1.0]
top = 10
steady = [True]
ibound = 1
#--determine the number of layers to evaluate
AllLay = np.arange(minLay,maxLay+1,1)
#--run simulations for each layer in AllLay
for klay in AllLay:
    print '\nRunning...{0:02d} layer model'.format( klay )
    laytype = np.ones((klay),np.int)
    layavg  = np.zeros((klay),np.int)
    aqb = dz = 30. 
    dz = aqb / float( klay )
    zall = -np.arange(dz,30.+dz,dz)
    tb = np.append(top,zall)
    #--create well data
    nwell = 4
    welllist = np.zeros((nper, nwell, 4),np.float)
    wellQ   = -1000.
    q       = wellQ / float( nwell )
    iw  = 0
    iwr = ny / 2
    iwc = nx / 2
    for i in range(iwr,iwr+2):
        for j in range(iwc,iwc+2):
            welllist[0,iw,0] = klay
            welllist[0,iw,1] = i
            welllist[0,iw,2] = j
            welllist[0,iw,3] = q
            iw += 1
    #--create ghb data on left and right model boundaries
    Nghb = ny * 2 * klay
    lrchc = np.zeros((nper,Nghb,5),np.float)
    hb  = [10., 0.]
    jb  = [1, nx]
    sgn = [1.0, -1.0]
#    dx = 1000.0 / float( nx )
    dhdx = ( hb[0] - hb[1] ) / xlen 
    ipos = 0
    for k in range(0,klay):
        for i in range(0,ny):
            for j in range(0,2):
                h = hb[j] + sgn[j] * dhdx * dx / 2.0
                b    = min( hb[j], tb[k] ) - tb[k+1]
                c    = kavg * dx * b / dx
                lrchc[0,ipos,0] = k + 1
                lrchc[0,ipos,1] = i + 1
                lrchc[0,ipos,2] = jb[j]
                lrchc[0,ipos,3] = h
                lrchc[0,ipos,4] = c
                ipos += 1
    #--base solver data
    npcond = 1
    relax  = 1.0
    nbpol  = 0
    iprpcg = 0
    mutpcg = 1
    damp   = 1.0
    #--solver iteration data
    mxiter = 50
    iter1  = 1000
    #--solver convergence data
    hclose = 1e-3
    rscale = 1000. * 1000. / float( nx * ny )
    rclose = hclose * rscale
    #--base simulation file name
    modelname = os.path.join( '..','{0:05d}'.format(nx),'GPU_{0:05d}'.format(nx) )
    #--output file names that are common to all runs
    listname = os.path.join( '..','{0:05d}'.format(nx),'GPU_{0:05d}.list'.format(nx) )
    headname = os.path.join( '..','{0:05d}'.format(nx),'GPU_{0:05d}.hds'.format(nx) )
    #--write base MODFLOW packages
    print 'Writing',klay,'layer base MODFLOW packages - PCG with MIC'
    i = smf.write_disfile(modelname,klay,ny,nx,nper,dx,dy,top,zall,perlen,nstp,tsmult,steady)
    i = smf.write_basfile(modelname,klay,ibound,0.0)
    i = smf.write_lpffile(modelname,klay,laytype,layavg,k_file,False)
    i = smf.write_welfile(modelname,nper,nwell,welllist)
    i = smf.write_ghbfile(modelname,nper,Nghb,lrchc)
    i = smf.write_ocfile(modelname,nper,nstp)
    #--run MODFLOW model - PCG with MIC
    if isPCG == True:
        i = smf.write_pcgfile(modelname,mxiter,iter1,npcond,hclose,rclose,relax,nbpol,iprpcg,mutpcg,damp)
        i = smf.write_namfile(modelname,['LIST','DIS','BAS6','LPF','WEL','GHB','PCG','OC'])
        print 'Running',klay,'layer base model with PCG and MIC preconditioner'
        a = subprocess.Popen([exe_name,modelname+'.nam'],stdout=subprocess.PIPE)
        b = a.communicate()
        c = b[0].split('\n')
        for cc in c:
            print cc
        c = []
        #--save list file
        savelistname = os.path.join('..','{0:05d}'.format(nx),'LISTArchive','GPU_{0:05d}_{1:02d}Layers_PCG_MIC_CPU.list'.format(nx,klay))
        print 'Copying...{0} to...\n...{1}\n'.format(listname, savelistname)
        shutil.copyfile(listname,savelistname)
        #--save hds file
        if saveHeads == True:
            saveheadname = os.path.join('..','{0:05d}'.format(nx),'HEADArchive','GPU_{0:05d}_{1:02d}Layers_PCG_MIC_CPU.hds'.format(nx,klay))
            print 'Copying...{0} to...\n...{1}\n'.format(headname, saveheadname)
            shutil.copy2(headname,saveheadname)
    #GMG
    if isGMG == True:
        ioutgmg = 1 #0
        iadamp  = 0
        ism     = 1 #0 ILU(0)
        isc     = 1
        nsize   = nx*ny #*klay
        l2norm  = rclose #CalcL2NormFromRclose(nsize,rclose)
        i = smf.write_gmgfile(modelname,mxiter,iter1,hclose,l2norm,damp,iadamp,ioutgmg,ism,isc,relax)
        i = smf.write_ocfile(modelname,nper,nstp)
        i = smf.write_namfile(modelname,['LIST','DIS','BAS6','LPF','WEL','GHB','GMG','OC'])
        #--run MODFLOW model - GMG
        print 'Running',klay,'layer base model with GMG solver'
        a = subprocess.Popen([exe_name,modelname+'.nam'],stdout=subprocess.PIPE)
        b = a.communicate()
        c = b[0].split('\n')
        for cc in c:
            print cc
        c = []
        #--save list file
        savelistname = os.path.join('..','{0:05d}'.format(nx),'LISTArchive','GPU_{0:05d}_{1:02d}Layers_GMG_CPU.list'.format(nx,klay))
        print 'Copying...{0} to...\n...{1}\n'.format(listname, savelistname)
        shutil.copyfile(listname,savelistname)
        #--save hds file
        if saveHeads == True:
            saveheadname = os.path.join('..','{0:05d}'.format(nx),'HEADArchive','GPU_{0:05d}_{1:02d}Layers_GMG_CPU.hds'.format(nx,klay))
            print 'Copying...{0} to...\n...{1}\n'.format(headname, saveheadname)
            shutil.copy2(headname,saveheadname)

    #--cycle through each pre-conditioner
    if isUPCG == True:
        for cpc in pc_list:
            ipc = int( cpc )
            #--set preconditioner text string
            pc = ''
            if ipc == 1:
                pc = 'Jacobi'
            elif ipc == 3:
                pc = 'MILU0'
            elif ipc == 4:
                pc = 'GLSPOLY'
            #--cycle through each UPCG hardware option
            for chardware in hardware_list:
                iopt = int( chardware )
                #--set hardware text string
                copt = ''
                if iopt == 1:
                    copt = 'CPU'
                elif iopt == 2:
                    copt = 'CPU_OpenMP'
                elif iopt == 3:
                    copt = 'GPU'
                print '\n\nstarting...\n  ',klay,'layer model using UPCG solver and\n  ',pc,'preconditioner run on',copt
                #default options
                idegree=10
                ilanstep=-2
                itrd=-4
                itrdv=1
                #--write UPCG files
                i = smf.write_upcgfile(modelname,mxiter,iter1,ipc,iopt,idegree,ilanstep,itrd,itrdv,hclose,rclose)
                i = smf.write_ocfile(modelname,nper,nstp)
                i = smf.write_namfile(modelname,['LIST','DIS','BAS6','LPF','WEL','GHB','UPCG','OC'])
                #--run MODFLOW model - UPCG soolver
                print 'Running',klay,'layer MODFLOW model using UPCG solver and\n  ',pc,'preconditioner run on',copt
                a = subprocess.Popen([exe_name,modelname+'.nam'],stdout=subprocess.PIPE)
                b = a.communicate()
                c = b[0].split('\n')
                for cc in c:
                    print cc
                c = []
                #--save list file
                savelistname = os.path.join('..','{0:05d}'.format(nx),'LISTArchive','GPU_{0:05d}_{1:02d}Layers_UPCG_{2}_{3}.list'.format(nx,klay,pc,copt))
                print 'Copying...{0} to...\n...{1}\n'.format(listname, savelistname)
                shutil.copyfile(listname,savelistname)
                #--save hds file
                if saveHeads == True:
                    saveheadname = os.path.join('..','{0:05d}'.format(nx),'HEADArchive','GPU_{0:05d}_{1:02d}Layers_UPCG_{2}_{3}.hds'.format(nx,klay,pc,copt))
                    print 'Copying...{0} to...\n...{1}\n'.format(headname, saveheadname)
                    shutil.copy2(headname,saveheadname)
                #--summary information
                print 'finished...\n  ',klay,'layer model using UPCG solver and\n  ',pc,'preconditioner run on',copt
    #--clean up memory
    zall = []
    tb = []
    #--wells
    welllist = []
    #--ghbs
    lrchc = []

print '\n\n...Successful termination'



