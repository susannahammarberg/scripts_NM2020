"""
Offline data preparation and reconstruction for standard NanoMAX ptycho.

This script is adapted for and requires ptypy 0.3. It is tested with
MPI on the compute cluster.


"""

import sys
import numpy as np
import ptypy
from ptypy.core import Ptycho
from ptypy import utils as u
from distutils.version import LooseVersion
import matplotlib.pyplot as plt

from ptypy.core import classes

import time
date_str = time.strftime("%Y%m%d_%H%M")


## simplest possible input #############################################
detector = 'merlin' # or 'merlin' or 'pilatus'
folder = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'
#scannr = 470 #Bragg InGaP center.    445#         # 397#    


scannr = 472#470
distance = 1.0
energy = 10
########################################################################

# General parameters
p = u.Param()
p.verbose_level = 3

p.run = str(scannr) + '_' + date_str 
# the io.home part tells ptypy where to save recons and dumps
p.io = u.Param()
p.io.home = './'
p.io.home = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/2drecons/'

p.io.autosave = u.Param()
p.io.autosave.interval = 1

# Scan parameters
p.scans = u.Param()
p.scans.scan00 = u.Param()
p.scans.scan00.name = 'Full'
#p.scans.scan00.coherence = u.Param()
#p.scans.scan00.coherence.num_probe_modes = 2		# Number of probe modes



p.scans.scan00.data = u.Param()
p.scans.scan00.data.name = 'NanomaxContrast' 
p.scans.scan00.data.path = folder
p.scans.scan00.data.detector = detector
p.scans.scan00.data.maskfile = {
    'merlin': 'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/merlin_mask_200430_8keV.h5', #'/data/visitors/nanomax/common/masks/merlin/latest.h5',
    'pilatus': None,
    'eiger': None,
     }[detector]
p.scans.scan00.data.scanNumber = scannr
p.scans.scan00.data.xMotor = 'npoint_buff/x' #'sx' for stepscan  
p.scans.scan00.data.yMotor = 'npoint_buff/y' #'sy'


#p.scans.scan00.data.zDetectorAngle = +0.78

# "Angle of the motor x axis relative to the lab x axis"
#TODO Check!
print('testing to use 0 here, cos the motor is moving in the lab system, maynbe. usually used 8.8')
p.scans.scan00.data.xMotorAngle = 8.8 #0 ####8.8 
#p.scans.scan00.data.shape = 110
#p.scans.scan00.data.shape = 224
##############_**************************************TGTG¤GGGGGGGGGGGGGGGGGG""""""""""""""""""""""""
#print('*****ops shape************')
p.scans.scan00.data.shape = 256# 340# 170 #512
p.scans.scan00.data.save = None
#--------- Set 445 etc
# test GainP
p.scans.scan00.data.center = (148,120+50-14)# GaInP  (med hela bilden) 
# InP close to bragg
##p.scans.scan00.data.center = (148, 345)# InPx364

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""Experiment """
#p.scans.scan00.data.center = (148, 345-256   -8)# InPx364    Use this when ptypy hacked
 
#---------
#--------- Set with bragg inp at 
#InP
#p.scans.scan00.data.center = (379,348)

#p.scans.scan00.data.center = (162,334)

# (do_transpose, do_flipud, do_fliplr)
p.scans.scan00.data.orientation = {
    'merlin': (False, False, True),
    'pilatus': None, 'eiger': None
     }[detector]
p.scans.scan00.data.distance = distance
p.scans.scan00.data.psize = {'pilatus': 172e-6,
                             'merlin': 55e-6,
                             'eiger': 75e-6 }[detector]
p.scans.scan00.data.energy = energy
p.scans.scan00.data.I0 = 'alba2/1'# None # can be like 'alba2/1'
p.scans.scan00.data.min_frames = 10
p.scans.scan00.data.load_parallel = 'all'

###does not work in master bracnh. does not rescale according to pixel size. need load_probe branch
p.scans.scan00.illumination = u.Param()

loaded_profile = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\siemensstar\scan14\np_save\focus\probe14_focus.npy')    
# save a psize, shape and the array data in the contaioner
Cprobe = ptypy.core.Container(data_dims=2, data_type='complex128')
Sprobe = Cprobe.new_storage(psize=[5.96673929e-08, 5.96673929e-08], shape=(1,128,128))
# fill storage
Sprobe.fill(loaded_profile)

p.scans.scan00.illumination.model = Sprobe

#p.scans.scan00.illumination.model = 'recon' #'stxm' 
#p.scans.scan00.illumination.recon = u.Param()
#p.scans.scan00.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan14\scan14_ML_0200.ptyr'
#p.scans.scan00.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan6\scan6_ML_0200.ptyr'
p.scans.scan00.illumination.aperture = None #need this to rescale the probe in branch 'load_probes'

# use this to out some initial intensity in the other probe modes, because as of now it is failing to do probe modes
#p.scans.scan00.illumination.diversity = u.Param()

#p.scans.scan00.illumination.diversity.noise = (10,30E-9)   #NM(.5, 1.0)
#p.scans.scan00.illumination.diversity.power = (0.8) #,0.6)#0.1#.1

#p.scans.scan00.illumination.diversity.power = (0.8,0.2) ##aaron
#p.scans.scan00.illumination.diversity.noise = (2.0,3.0)

# Reconstruction parameters
p.engines = u.Param()
p.engines.engine00 = u.Param()
p.engines.engine00.name = 'ePIE'
p.engines.engine00.numiter = 2000
p.engines.engine00.numiter_contiguous = 200
#""""OBSOBS test med  rörlig probe                                       """
p.engines.engine00.probe_update_start = 10000 ###(default 2)

#p.engines.engine00.position_refinement = u.Param()
#p.engines.engine00.position_refinement.start = 0
#p.engines.engine00.position_refinement.stop = 20
#p.engines.engine00.position_refinement.interval = 1
#p.engines.engine00.position_refinement.nshifts = 4
#p.engines.engine00.position_refinement.amplitude = 10E-9   # ??
#p.engines.engine00.position_refinement.max_shift = 40E-9  #??
#p.engines.engine00.position_refinement.record = False
#average probe in each node "Averaging seems to work the best."
#TODO evaluate this 20210915 added this
#p.engines.engine00.average_probe = True

#p.engines.engine01 = u.Param()
#p.engines.engine01.name = 'DM'
#p.engines.engine01.numiter = 500
#p.engines.engine01.numiter_contiguous = 100
#""""OBSOBS test med  rörlig probe                                       """
#p.engines.engine01.probe_update_start = 200 ###(default 2)


#TODO try this
# Clip object amplitude into this interval (tuple)
#p.engines.engine00.clip_object = 
#p.engines.engine00.beta = 0.7
#p.engines.engine00.alpha = 1E-4 #TOTRY 20220317

#does not work with ePIE (support is equal to none)
p.engines.engine00.probe_support = 0.7

P = Ptycho(p,level=5)

#%%


# run this with 1 probe mode only! Otherwise there are Nxy*n number of views

import matplotlib
matplotlib.use( 'Qt5agg' )

diff_data = P.diff.storages['S0000'].data*P.mask.storages['S0000'].data #P.diff.storages.values()[0].data*P.mask.storages.values()[0].data[0]
##
##position = [int(len(diff_data)/2)]
##
plt.figure()
plt.imshow((np.log10(sum(diff_data))),cmap='magma', interpolation='none')
plt.title('Summed intensity ()')
plt.savefig('sum_intensity')
plt.colorbar()
plt.show()


probe = np.squeeze(P.probe.storages['Sscan00G00'].data )

plt.figure()
plt.imshow(abs(loaded_profile),cmap='jet',interpolation=None)
plt.figure()
plt.imshow((abs(probe)),cmap='jet',interpolation=None)

#np.save('probe14_resapled_to_shape_170_reso_13nm',probe)
##
##obj_storage = P.obj.storages['Sscan00G00']
##
##plt.figure()
##plt.imshow(abs(np.squeeze(obj_storage.data)),cmap='magma')
##plt.show()
##
g = P.pods['P0000'].geometry
g.resolution

##
### Keep a copy of the object storage, and fill the actual one with an
### initial guess (like zeros everywhere).
##S_true = obj_storage.copy(owner=P.obj)
##
##
##z, y = S_true.grids()
##
##zcut = int(S_true.data.shape[1]/2)
##ycut = int(S_true.data.shape[2]/2)
##
##fact = 1E6
####%%
#### Control directions of motorpositions (control the order of diffraction patterns)
###obj_storage.fill(0.0)
###plt.ion()
###diff_frames=[]
###plt.ioff()
###for j, obj_view in enumerate(P.obj.views):
###    
###    plt.figure()
###    plt.imshow( abs(P.obj.views[obj_view].data[20]), cmap='RdBu_r' )
###    plt.title('pos ' + "{0:03}".format(j) )
###    plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\scanning_positions\\' + "{0:03}".format(j)) 
###    print(j)
###    
##
###%%
##
### if you dont want to start over and just continue to iteratie, just comment this part out
### zero everything
##
##print('Initialize/reset ''manual'' reconstruction script ')
### reset obj. storage
##storage_save = []
##
##
##
##obj_storage.fill(0.0)
##
##
##
###%%
##
##print('Start ''manual'' reconstruction script ')
##import time
##date_str = time.strftime("%Y%m%d_%H%M")
##algorithm = 'PIE'
##iterations = 300
##
### Here's a PIE/cPIE implementation
##if algorithm == 'PIE':
##    beta = 1.0
##    eps = 1e-3
##    
##    ferrors = []
##
##    for i in range(iterations):
##        
##        ferrors_ = []
##        
##        j = 0       
##        
##        # iterate over all views
##        for obj_view in P.obj.views:   
##
##            exit_ = P.obj.views[obj_view].data * probe
##            prop = g.propagator.fw(exit_)
##            ferrors_.append(abs(np.abs(prop)**2 - diff_data[j]))
##            #import pdb; pdb.set_trace()
##            prop[:] = np.sqrt(diff_data[j]) * np.exp(1j * np.angle(prop))
##            exit = g.propagator.bw(prop)
##            # ePIE scaling (Maiden2009)
##            #SH: det här är väl  J.M. Rodenburg and H.M.L Faulkner 2004, jag skrev exakt såhär baserat på det peket.
##            P.obj.views[obj_view].data += beta * np.conj(probe) / (np.abs(probe).max())**2 * (exit - exit_)
##            # PIE and cPIE scaling (Rodenburg2004 and Godard2011b)
##            #P.obj.views[obj_view].data += beta * np.abs(probe) / np.abs(probe).max() * np.conj(probe) / (np.abs(probe)**2 + eps) * (exit - exit_)
##            
##            #probe function is not updating
##                        
##            #update probe function
##            beta = 0.9
##            #b#eta = 0.1                       
##            probe = probe + beta *(exit -exit_) * np.conj(P.obj.views[obj_view].data)/ (np.max(abs(P.obj.views[obj_view].data))**2)
##            
##            
##            # apply constraint to 10 pixlar = 240 nm (á 24nm per pix)
##            #obj_storage.data[0][0:20-6] = 1E-9 * np.exp(1j*np.angle(obj_storage.data[0][0:20-6]))
##            #obj_storage.data[0][20+6:] = 1E-9 * np.exp(1j*np.angle( obj_storage.data[0][20+6:]))
##
##
##            j += 1
##            
##
##        ferrors.append(np.mean(ferrors_))
##
##        
##        #storage_save.append(obj_storage.copy(owner=P.obj))
##        if not (i % 10):
####            #save arrays instread of plottting to improve speed
##             storage_save.append(obj_storage.copy(owner=P.obj))
##        print('iteration ', i)
##
##
##
###%%
##date_str = time.strftime("%Y%m%d_%H%M")
##savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\2drecons\recons_with_ptypy_views' + '\\ePIE_%d_'%scannr + date_str
##import os
##if not os.path.exists(savepath):
##    os.makedirs(savepath)
##    print('new folder in this savepath was created')
##    
##    
##for ii, store in enumerate(storage_save):
###if not (ii % 10):
##    
##    #ii=
##    
##    store = storage_save[29]
##            
##    obj = np.squeeze(store.data)
##    #    obj = np.squeeze(obj_storage.data)
##    
##    fig = plt.figure()    # row col  row col
##    ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=3)
##    ax1 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=2)
##    ax2 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=2)
##    ax0.plot(ferrors,'blue',marker='.', label='Fourier error')
##    extent_zy = [ fact*y.min(), fact*y.max(), fact*z.min(), fact*z.max()]
##    
##    im1 = ax1.imshow(np.log10( np.abs(obj)) , extent=extent_zy, interpolation='none', cmap='gray')#, vmin=0.02, vmax = 0.2)
##    plt.setp(ax1, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='log amplitude')
##    plt.colorbar(im1,ax=ax1)
##    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=70)
##    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
##    #plt.savefig(savepath + '\\amp_iter%d_'%(ii))      
##    #thight layout should minimize the overlab of labels
##    plt.tight_layout() 
##
##    amplitude_range = np.abs(obj).max()-np.abs(obj).min()
##    mask_at = 0.00005 #np.abs(obj).min() + amplitude_range/2
##    mask = np.zeros((obj.shape))
##    mask[np.abs(obj) > mask_at] = 1
##    mask = np.squeeze(mask)
##    
##    im2 = ax2.imshow( mask * np.angle(obj), extent=extent_zy, interpolation='none', cmap='jet')
##    plt.setp(ax2, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='phase')
##    fig.colorbar(im2,ax=ax2 )
##    plt.tight_layout()
##    #plt.savefig(savepath + '\\obj_iter%d_'%(ii*10))
##    plt.show()
##
##    plt.figure()
##    plt.imshow(abs(probe), cmap='gray', interpolation='none')#, extent=[0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
##    plt.xlabel(' [µm]')
##    plt.ylabel(' [µm]')
##    #plt.title('Scan %d: Probe amplitude'%scan)
##    plt.colorbar()
##    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_Pamp_k%d'%(scan, k), bbox_inches='tight')
##    
##    plt.figure()                                                            # horisontalt vertikalt
##    plt.imshow(np.angle(probe), cmap='jet', interpolation='none')#, extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
##    plt.xlabel(' [µm]')
##    plt.ylabel(' [µm]')
##    #plt.title('Scan %d: Probe phase'%scan)
##    plt.colorbar()
##    ##plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_Pphase_k%d'%(scan, k), bbox_inches='tight')  
##    
##    
##    # this is not correct.. its == Q*u()
##    # this doesn make sense... but it makes the most sense if I want to compare to the XRD maps
##    
###    S_cart.data.shape
###    z_middle = int(S_cart.data.shape[2]/2)
###    y_middle = int(S_cart.data.shape[3]/2)
###    slicez = slice(z_middle-150,z_middle+150)
###    slicey = slice(y_middle-80,y_middle+80)       
##
##    #Average over the reigion of interest (that is a bit arbitrary)
##    #mean_phase_x = np.mean(np.angle(np.squeeze(S_cart.data)[xcut][slicex, slicey]))
##    #TODO not plotting strain now
##    #strain_xcut = np.angle(np.squeeze(S_cart.data)[xcut][slicex, slicey])
##    #strain_xcut = 100*(np.angle(np.squeeze(S_cart.data)[xcut][slicex, slicey]) - mean_phase_x) / mean_phase_x
##    
##    ###extent_zy_sliced = extent_zy = [fact*z.min(), fact*z.max(), fact*y.min(), fact*y.max()][slicez, slicey]
##    
###    fig, ax = plt.subplots(ncols=1)  
###    im3 = ax.imshow((np.abs(np.squeeze(S_cart.data)[xcut][slicez, slicey]) ).T,  interpolation='none', origin='lower', cmap='RdBu_r')
###    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
###    fig.colorbar(im3,ax=ax )
####    plt.savefig(savepath + '\\amp_zoomed_iter%d_'%(ii))
###    plt.show()
##    
##    
###    # mean value to compare to XRD
###    strain_ = np.mean(np.angle(np.squeeze(S_cart.data)[19-3:19+3][slicex, slicey]), axis = 0)
###    fig, ax = plt.subplots(ncols=1)  
###    im3 = ax.imshow((np.mean(np.abs(np.squeeze(S_cart.data)[19-3:19+3][slicex, slicey]), axis = 0) * strain_).T, extent=[fact* z.min(), fact*z.max(), fact*y.min(), fact*y.max()], interpolation='none', origin='lower', cmap='jet')
###    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='mean value phase \n front view')
###    fig.colorbar(im3,ax=ax )
###    
###    #mean value in amplitude
###    mean_amp = np.sum(np.abs(np.squeeze(S_cart.data)[19-3:19+3]), axis = 0)
###    fig, ax = plt.subplots(ncols=1)  
###    im3 = ax.imshow(mean_amp.T, extent=[fact* z.min(), fact*z.max(), fact*y.min(), fact*y.max()], interpolation='none', origin='lower', cmap='jet')
###    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='mean value amplitude \n front view')
###    fig.colorbar(im3,ax=ax )
##
###save the last object result
###np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\opject_ePIE_result', S_cart.data)
##
##
###plot_strain()
##
##
##        
##
