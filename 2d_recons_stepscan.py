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
#folder = '/data/visitors/nanomax/20200567/2020052308/raw/sample'
folder = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'
scannr = 395
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
p.io.home = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/2drecons/'

p.io.autosave = u.Param()
p.io.autosave.interval = 1

# Scan parameters
p.scans = u.Param()
p.scans.scan00 = u.Param()
p.scans.scan00.name = 'Full'
#p.scans.scan00.coherence = u.Param()
#p.scans.scan00.coherence.num_probe_modes = 5		# Number of probe modes

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
p.scans.scan00.data.xMotor = 'sy' #'sx' for stepscan  
p.scans.scan00.data.yMotor = 'sx' #'sy'
# "Angle of the motor x axis relative to the lab x axis"
#TODO Check!
p.scans.scan00.data.xMotorAngle = 8.8 
#p.scans.scan00.data.shape = 110
p.scans.scan00.data.shape = 198
p.scans.scan00.data.save = None
p.scans.scan00.data.center = (359,359)    

# (do_transpose, do_flipud, do_fliplr)
p.scans.scan00.data.orientation = {
    'merlin': (False, False, True),
    
    ##'merlin': (False, False, False),
    ##'merlin': (True, False, False),
    ##'merlin': (False, True, False),
    ##'merlin': (False, False, True),
    ##'merlin': (True, True, False),
    ##'merlin': (False, True, True),
    ##'merlin': (True, True, True),
    
    'pilatus': None, 'eiger': None
     }[detector]
p.scans.scan00.data.distance = distance
p.scans.scan00.data.psize = {'pilatus': 172e-6,
                             'merlin': 55e-6,
                             'eiger': 75e-6 }[detector]
p.scans.scan00.data.energy = energy
#p.scans.scan00.data.I0 = 'alba2/1'# None # can be like 'alba2/1'
p.scans.scan00.data.min_frames = 10
p.scans.scan00.data.load_parallel = 'all'

###does not work in master bracnh. does not rescale according to pixel size. need load_probe branch
##p.scans.scan00.illumination = u.Param()
##p.scans.scan00.illumination.model = 'recon' #'stxm' #recon'
##p.scans.scan00.illumination.recon = u.Param()
##p.scans.scan00.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan14\scan14_ML_0200.ptyr'
###p.scans.scan00.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan6\scan6_ML_0200.ptyr'
##p.scans.scan00.illumination.aperture = None

#######funkar nog inte med riktig probe
p.scans.scan00.illumination = u.Param()
p.scans.scan00.illumination.model = None
p.scans.scan00.illumination.aperture = u.Param() 
p.scans.scan00.illumination.aperture.form = 'circ'
p.scans.scan00.illumination.aperture.size = 450e-9 
#p.scans.scan00.illumination.propagation = u.Param()
#p.scans.scan00.illumination.propagation.focussed =
#p.scans.scan00.illumination.propagation.spot_size = 50e-9


# Reconstruction parameters
p.engines = u.Param()
p.engines.engine00 = u.Param()
p.engines.engine00.name = 'DM'
p.engines.engine00.numiter = 100
p.engines.engine00.numiter_contiguous = 1
p.engines.engine00.probe_update_start = 1

p.engines.engine00.probe_support = 0.7

##p.engines.engine01 = u.Param()
##p.engines.engine01.name = 'ML'
##p.engines.engine01.numiter = 100#100
##p.engines.engine01.numiter_contiguous = 10
##p.engines.engine01.probe_update_start = 1


#p.engines.engine00.sample_support = None # by default i think its 1e-3
# sample support option if in 3d_Bragg engine
#p.engines.engine00.sample_support = u.Param()
#p.engines.engine00.sample_support.coefficient = 0.0 # "Sample amplitude is multiplied by this value outside the support region"



if LooseVersion(ptypy.version) < LooseVersion('0.3.0'):
    raise Exception('Use ptypy 0.3.0 or better!')

P = Ptycho(p,level=5)

#%%

import matplotlib
matplotlib.use( 'Qt5agg' )

diff_data = P.diff.storages['S0000'].data*P.mask.storages['S0000'].data #P.diff.storages.values()[0].data*P.mask.storages.values()[0].data[0]

position = [int(len(diff_data)/2)]

plt.figure()
plt.imshow((np.log10(sum(diff_data))),cmap='magma', interpolation='none')
plt.title('Summed intensity ()')
plt.savefig('sum_intensity')
plt.colorbar()
plt.show()
##
##probe = np.squeeze(P.probe.storages['Sscan00G00'].data )
##
##obj_storage = P.obj.storages['Sscan00G00']
##
##plt.figure()
##plt.imshow(abs(np.squeeze(obj_storage.data)),cmap='magma')
##
##g = P.pods['P0000'].geometry
##
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
### Here's a PIE/cPIE implementation
##if algorithm == 'PIE':
##    beta = 1.0
##    eps = 1e-3
##    
##    ferrors = []
##
##    for i in range(10):
##        
##        ferrors_ = []
##        
##        j = 0       
##        
##        # iterate over all views
##        for obj_view in P.obj.views:   
##        #for j in range(len(P.obj.views)):
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
##            # probe function is not updating
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
##        storage_save.append(obj_storage.copy(owner=P.obj))
###        if not (i % 5):
####            #save arrays instread of plottting to improve speed
###             storage_save.append(obj_storage.copy(owner=P.obj))
##        print(i)
##
##
##
###%%
##
##savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons_ePIE' + '\\' + date_str
##import os
##if not os.path.exists(savepath):
##    os.makedirs(savepath)
##    print('new folder in this savepath was created')
##    
##
##for ii, store in enumerate(storage_save):
##
##    date_str = time.strftime("%Y%m%d_%H%M")
##
##
##    
##    
##    fig = plt.figure()    # row col  row col
##    ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=3)
##    ax1 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=2)
##    ax2 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=2)
##    
##    ax0.plot(ferrors,'blue',marker='.', label='Fourier error')
##    #ax0.legend(bbox_to_anchor=(0.65, 0.5), loc='center left')
##    obj = np.squeeze(obj_storage.data)
##    
##
##    extent_zy = [ fact*y.min(), fact*y.max(), fact*z.min(), fact*z.max()]
##    
##    
##    im1 = ax1.imshow( np.abs(obj) , extent=extent_zy, interpolation='none', cmap='gray')#RdBu_r')#,vmin=min_scatt, vmax=max_scatt)
##    plt.setp(ax1, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='2d recon')
##    plt.colorbar(im1,ax=ax1)
##    
##    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=70)
##    #plt.setp(ax3.xaxis.get_majorticklabels(), rotation=70)
##    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
##    #plt.savefig(savepath + '\\amp_iter%d_'%(ii))
##    
##    
##    #thight layout should minimize the overlab of labels
##    plt.tight_layout() 
###    plt.draw()
###    plt.show()
##
##    mask = np.zeros((obj_storage.shape))
##    mask[abs(obj_storage.data)>9E-8] = 1 
##    mask = np.squeeze(mask)
##    
##    #fig, ax = plt.subplots(ncols=3)
##    #plt.suptitle('Phase central cuts.')
##    im2 = ax2.imshow( 1 * np.angle(obj), extent=extent_zy, interpolation='none', cmap='jet')
##    plt.setp(ax2, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='phase')
##    fig.colorbar(im2,ax=ax2 )
##    
##
##    
##    
##    plt.tight_layout()
###    plt.draw()
###    plt.savefig(savepath + '\\phase_iter%d_'%(ii))
##    plt.show()
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
##    fig, ax = plt.subplots(ncols=1)  
##    im3 = ax.imshow((np.abs(np.squeeze(S_cart.data)[xcut][slicez, slicey]) ).T,  interpolation='none', origin='lower', cmap='RdBu_r')
##    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
##    fig.colorbar(im3,ax=ax )
###    plt.savefig(savepath + '\\amp_zoomed_iter%d_'%(ii))
##    plt.show()
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
