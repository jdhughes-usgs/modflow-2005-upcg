#--libraries
import numpy as np
import os
import math
import sys

def loadArrayFromFile(nrow,ncol,file):
    '''
    read 2darray from file
    file(str) = path and filename
    '''
    try:
        file_in = open(file,'r')
        openFlag = True
    except:
#       assert os.path.exists(file)
        file_in = file
        openFlag = False
    
    data = np.zeros((nrow*ncol),dtype=np.float) #-1.0E+10
    data.fill( -1.0E+10 )
    d = 0
    while True:
        line = file_in.readline()
        if line is None or d == nrow*ncol:break
        raw = line.strip('\n').split()
        for a in raw:
            try:
                data[d] = float(a)
            except:
                print 'error casting to float on line: ',line
                sys.exit()
            if d == (nrow*ncol)-1:
                assert len(data) == (nrow*ncol)
                data.resize(nrow,ncol)
                return(data) 
            d += 1  
    file_in.close()
    data.resize(nrow,ncol)
    return(data)


def write_disfile(fn_path,nlay,nrow,ncol,nper,dx,dy,top,botm,perlen,nstp,tsmult,steady):
    # Open file for writing
    f_dis = open(fn_path+'.dis', 'w')
    # Item 0: heading
    f_dis.write('%s\n' % '#')
    # Item 1: NLAY, NROW, NCOL, NPER, ITMUNI, LENUNI
    f_dis.write('%10d%10d%10d%10d%10d%10d\n' % (nlay, nrow, ncol, nper, 4, 2))
    # Item 2: LAYCBD
    # Check if all items must/can be on 1 single line
    for l in range(0, nlay):
        f_dis.write('%3d' % (0))
    f_dis.write('\n')
    # Item 3: DELR
    f_dis.write('CONSTANT %10f           (5G13.0)         -1 DELR(NCOL)\n' % (dx))
    # Item 4: DELC
    f_dis.write('CONSTANT %10f           (5G13.0)         -1 DELC(NROW)\n' % (dy))
    # Item 5: Top(NCOL, NROW)
    f_dis.write('CONSTANT %10f           (5G13.0)         -1 TOP OF SYSTEM\n' % (top))
    # Item 5: BOTM(NCOL, NROW)
    for l in range(0, nlay):
        f_dis.write('CONSTANT %10f           (5G13.0)         -1 BOTTOM OF LAYER %d\n' % (botm[l],l))
    # Item 6: NPER, NSTP, TSMULT, Ss/tr
    for t in range(nper):
        f_dis.write('%14f%14d%10f' % (perlen[t],nstp[t],tsmult[t]))
        if steady[t]:
            f_dis.write('%3s\n' % 'SS')
        else:
            f_dis.write('%3s\n' % 'TR')
    f_dis.close()
    return 1


def write_basfile(fn_path,nlay,ibound=1,head=0.0):
    # Open file for writing
    f_bas = open(fn_path+'.bas6', 'w')
    # First line: heading
    f_bas.write('%s\n' % '#')
    # Second line: format specifier
    options = '  FREE'
    f_bas.write('%s\n' % options)
    # IBOUND array
    for l in range(0, nlay):
        f_bas.write('CONSTANT %10d           (10I4)           -1 IBOUND Array for Layer %d\n' % (ibound,l+1))
    # Head in inactive cells
    f_bas.write('%f\n' % -999.990000)
    # Starting heads array
    for l in range(0, nlay):
        f_bas.write('CONSTANT %10f           (5G13.0)         -1 Starting Heads in Layer %d\n' % (head,l+1))
    # Close file
    f_bas.close()
    return 1

def write_lpffile(fn_path,nlay,laytype,layavg,ck,transient):
    # Open file for writing
    f_lpf = open(fn_path+'.lpf', 'w')
    # Item 0: text
    f_lpf.write('%s\n' % '#')
    # Item 1: IBCFCB, HDRY, NPLPF
    f_lpf.write('%10d%10.1e%10d\n' % (0, -1.0e+30, 0))
    # LAYTYP array
    for l in range(0, nlay):
        f_lpf.write('%10d' % (laytype[l]))
    f_lpf.write('\n')
    # LAYAVG array
    for l in range(0, nlay):
        f_lpf.write('%10d' % (layavg[l]))
    f_lpf.write('\n')
    # CHANI array
    for l in range(0, nlay):
        f_lpf.write('%10.1e' % (1.0))
    f_lpf.write('\n')
    # LAYVKA array
    for l in range(0, nlay):
        f_lpf.write('%10d' % (0))
    f_lpf.write('\n')
    # LAYWET array
    for l in range(0, nlay):
        f_lpf.write('%10d' % (0))
    f_lpf.write('\n')
    # Item 7: WETFCT, IWETIT, IHDWET
#    iwetdry = self.laywet.sum()
#    if iwetdry > 0:
#    	f_lpf.write('%10f%10d%10d\n' % (self.wetfct, self.iwetit, self.ihdwet))
#    transient = not self.parent.get_package('DIS').steady.all()
    for l in range(nlay):
        f_lpf.write('OPEN/CLOSE %s %10f           (FREE)         -1 HK Layer %d\n' % (ck,1.0,l+1))
#    	if self.chani[i] < 1:
#        	comment = 'HANI() = Ratio of horizontal hydraulic of columns to rows of layer ' + str(i + 1)
#        	self.parent.write_array( f_lpf, self.hani[:,:,i], self.unit_number[0], True, 13, ncol, comment,ext_base='hani_'+str(i+1) )
        f_lpf.write('OPEN/CLOSE %s %10f           (FREE)         -1 VKA Layer %d\n' % (ck,1.0,l+1))
#    	if transient == True:
#    		comment = 'Ss() = Specific storage coefficient of layer ' + str(i + 1)
#    		self.parent.write_array( f_lpf, self.ss[:,:,i], self.unit_number[0], True, 13, ncol, comment,ext_base='ss_'+str(i+1) )
#    		if self.laytyp[i] !=0:
#    			comment = 'Sy() = Specific yield of layer ' + str(i + 1)
#    			self.parent.write_array( f_lpf, self.sy[:,:,i], self.unit_number[0], True, 13, ncol, comment,ext_base='sy_'+str(i+1) )
#    	if self.parent.get_package('DIS').laycbd[i] > 0:
#    		comment = 'VKCB() = Vertical hydraulic conductivity of quasi-three-dimensional confining bed of layer ' + str(i + 1)
#    		self.parent.write_array( f_lpf, self.vkcb[:,:,i], self.unit_number[0], True, 13, ncol, comment,ext_base='vkcb_'+str(i+1))
#    	if (self.laywet[i] != 0 and self.laytyp[i] != 0):
#    		comment = 'WETDRY() = Wetting threshold of layer ' + str(i + 1)
#    		self.parent.write_array( f_lpf, self.wetdry[:,:,i], self.unit_number[0], True, 13, ncol, comment,ext_base='wetdry_'+str(i+1) )
    f_lpf.close()
    return 1


def write_welfile(fn_path,nper,nwell,welllist):
    f_wel = open(fn_path+'.wel', 'w')
    f_wel.write('%s\n' % '#')
    f_wel.write('%10i%10i\n' % (nwell, 0))
    for np in range(0,nper):
        f_wel.write('%10i\n' % (nwell))
        for nw in range(0,nwell):
            f_wel.write(' %9i %9i %9i %13.6e\n' % ( welllist[np,nw,0], welllist[np,nw,1], welllist[np,nw,2], welllist[np,nw,3] ))
    f_wel.close()
    return 1
 
def write_ghbfile(fn_path,nper,nghb,ghblist):
    f_ghb = open(fn_path+'.ghb', 'w')
    f_ghb.write('%s\n' % '#')
    f_ghb.write('%10i%10i NOPRINT\n' % (nghb, 0))
    for np in range(0,nper):
        f_ghb.write('%10i\n' % (nghb))
        for ng in range(0,nghb):
            f_ghb.write(' %9i %9i %9i %13.6e %13.6e\n' % ( ghblist[np,ng,0], ghblist[np,ng,1], ghblist[np,ng,2], ghblist[np,ng,3], ghblist[np,ng,4] ))
    f_ghb.close()
    return 1

def write_pcgfile(fn_path,mxiter,iter1,npcond,hclose,rclose,relax,nbpol,iprpcg,mutpcg,damp):
    # Open file for writing
    f_pcg = open(fn_path+'.pcg', 'w')
    f_pcg.write('%s\n' % '#')
    f_pcg.write('%10i%10i%10i\n' % (mxiter,iter1,npcond))
    f_pcg.write('%10f%10f%10f%10i%10i%10i%10f\n' % (hclose,rclose,relax,nbpol,iprpcg,mutpcg,damp))
    f_pcg.close()
    return 1

def write_gmgfile(fn_path,mxiter,iter1,hclose,rclose,damp,iadamp,ioutgmg,ism,isc,relax):
    # Open file for writing
    f_gmg = open(fn_path+'.gmg', 'w')
    f_gmg.write('%s\n' % '#')
    f_gmg.write('{0:10.3e}{1:10d}{2:10.3e}{3:10d}\n'.format(rclose,iter1,hclose,mxiter))
    f_gmg.write('{0:10.3e}{1:10d}{2:10d}{3:10d}\n'.format(damp,iadamp,ioutgmg,0))
    f_gmg.write('{0:10d}{1:10d}\n'.format(ism,isc))
    if isc==4:
    	f_gmg.write('{0:10.3e}\n'.format(relax))
    f_gmg.close()
    return 1

    
def write_upcgfile(fn_path,mxiter,iter1c,npc,nopt,ndegree,nlanstep,ntrd,ntrdv,hclose,rclose):
    # Open file for writing
    f_upcg = open(fn_path+'.upcg', 'w')
    f_upcg.write('%s\n' % '#')
    f_upcg.write('#\n#--Data set 1\n')
    f_upcg.write('%10i%10i%10i%10i' % (mxiter,iter1c,npc,nopt))
    if npc == 4:
        f_upcg.write('%10i%10i' % (ndegree,nlanstep))
    if nopt == 2:
        f_upcg.write('%10i' % (ntrd))
        if ntrd < 0:
            f_upcg.write('%10i' % (ntrdv))
    f_upcg.write('\n')
        
    f_upcg.write('#\n#--Data set 2\n')
    f_upcg.write('%10f%10f\n' % (hclose,rclose))
    f_upcg.close()
    return 1

def write_ocfile(fn_path,nper,nstp):
    f_oc = open(fn_path+'.oc', 'w')
    f_oc.write('%s\n' % '#')
    f_oc.write('HEAD PRINT FORMAT   0\nHEAD SAVE UNIT    201\n')
    for np in range(0,nper):
        f_oc.write('\n')
        f_oc.write('PERIOD %10i STEP %10i\n' % (np+1,nstp[np]))
        f_oc.write(' SAVE HEAD\n')
        f_oc.write(' PRINT BUDGET\n')
    f_oc.close()
    return 1

def write_namfile(fn_path,allpkg):
    f_nam = open(fn_path+'.nam', 'w')
    f_nam.write('%s\n' % '#Name file for mf2005')
    ifile = 100
    for pkg in allpkg:
        ifile += 1
        f_nam.write('%s %i %s\n' % (pkg.upper(), ifile, fn_path+'.'+pkg.lower()) )
    #head file
    f_nam.write('DATA(BINARY) 201 %s REPLACE\n' % (fn_path+'.hds') )

    f_nam.close()
    return 1
