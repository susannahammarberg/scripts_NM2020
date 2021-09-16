"""
 Reconstruction script for NanoMAX MAy 2020
 flyscans
 
 shifting implemented
 IO normalization implemented
 
 
"""

import ptypy
from ptypy.core import Ptycho
from ptypy import utils as u
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from ptypy.experiment.nanomax3d import NanomaxBraggMay2020_shifting
import time
date_str = time.strftime("%Y%m%d_%H%M")

p = u.Param()
p.run = 'recon' + date_str





#original scans full range USE THIS WITH CURRENT SHIFTING LIST
#scans = np.arange(429,490).tolist()


# InP range
scans = np.arange(429,470).tolist()



# smaller InP range
#scans = np.arange(429,457).tolist()

#scans = np.arange(429,445+1).tolist()


#test range
#scans = np.arange(443,443+4).tolist()


# the most diffracting scans of the 170nm segment
#scans = np.arange(438,441+1).tolist()
#replacement scans
#scans = np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ]).tolist()


p.data_type = "single"   #or "double"
p.verbose_level = 3

# the io.home part tells ptypy where to save recons and dumps
p.io = u.Param()
p.io.home = './'
p.io.home = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/3drecons/'

# need this to save the dump I think
p.io.autosave = u.Param()
p.io.autosave.interval = 1 
# this requires pyzmq which i couldnt install
##p.io.autoplot = u.Param()
##p.io.autoplot.layout = 'bragg3d'
##p.io.autoplot.dump = True
##p.io.autoplot.interval = 1

p.scans = u.Param()
p.scans.scan01 = u.Param()
p.scans.scan01.name = 'Bragg3dModel'
p.scans.scan01.data = u.Param()
p.scans.scan01.data.name = 'NanomaxBraggMay2020_shifting'
###p.scans.scan01.data.name = 'NanomaxBraggMay2020'
p.scans.scan01.data.path = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'
p.scans.scan01.data.maskfile ='C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/merlin_mask_200430_8keV.h5'

p.scans.scan01.data.scans = scans
p.scans.scan01.data.xMotor = 'npoint_buff/x' 
p.scans.scan01.data.yMotor = 'npoint_buff/y' #'sy' (använde sy fram till nu 20210610)
#p.scans.scan01.data.xMotorFlipped = True

# Diffraction angle (theta, not two theta) in degrees
p.scans.scan01.data.theta_bragg = 8.8 #(#10.54 # for InP according from table value. If that is correct, I dont need 2theta. 

p.scans.scan01.data.shape = 170 
#center of RAW image
p.scans.scan01.data.cen = (148, 345)    # before: InP (148,364) GaInP:155
p.scans.scan01.data.center = None
p.scans.scan01.data.auto_center = False

##Alternatively, a 3-tuple of booleans may be provided (do_transpose, do_flipud, do_fliplr)
p.scans.scan01.data.orientation = (True, False, False)


##" without shifting"
#p.scans.scan01.data.vertical_shift = [0]*61
#p.scans.scan01.data.horizontal_shift = [0]*61
# shifting vectors for the full range original scans
p.scans.scan01.data.vertical_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\vertical_shift_vector.npy').tolist()
p.scans.scan01.data.horizontal_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\horizontal_shift_vector.npy').tolist()                               

#p.scans.scan01.data.load_parallel = None # 'all'
p.scans.scan01.data.psize = 55e-6
p.scans.scan01.data.energy = 10.0
p.scans.scan01.data.distance = 1.0
p.scans.scan01.data.I0 = 'alba2/1'
#p.scans.scan01.data.min_frames = 10
#p.scans.scan01.data.load_parallel = 'all'

p.scans.scan01.illumination = u.Param()
p.scans.scan01.illumination.model = 'recon' 
p.scans.scan01.illumination.recon = u.Param()
p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan14\scan14_ML_0200.ptyr'
#p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan13\scan13_ML_0200.ptyr'
p.scans.scan01.illumination.aperture = None #need this to rescale the probe in branch 'load_probes'

##p.scans.scan01.illumination.aperture = u.Param() 
##p.scans.scan01.illumination.aperture.form = 'circ'
##p.scans.scan01.illumination.aperture.size = 300e-9 

#p.scans.scan01.sample = u.Param()
#p.scans.scan01.sample.fill = 1E-20

p.engines = u.Param()
p.engines.engine00 = u.Param()
p.engines.engine00.name = 'DM_3dBragg'   
p.engines.engine00.numiter = 50
p.engines.engine00.probe_update_start = 100000
p.engines.engine00.numiter_contiguous = 1

p.engines.engine00.probe_support = None   #get an error without these
# sample support option if in 3d_Bragg engine
p.engines.engine00.sample_support = None #use with DM_3dBragg only

p.engines.engine00.sample_support = u.Param()
#p.engines.engine00.sample_support = True # True dont need this?
p.engines.engine00.sample_support.type = 'thinlayer' #thinlayer
p.engines.engine00.sample_support.size = 300E-9   #diameter. (long axis need fix in code). 
p.engines.engine00.sample_support.coefficient = 0.0
p.engines.engine00.sample_support.shrinkwrap = None

# need to fix which coordinates to wrap around in order to use shrinkwrap
# p.engines.engine00.sample_support.shrinkwrap = u.Param()
# p.engines.engine00.sample_support.shrinkwrap.cutoff = .3
# p.engines.engine00.sample_support.shrinkwrap.smooth = None
# p.engines.engine00.sample_support.shrinkwrap.start = 15
# p.engines.engine00.sample_support.shrinkwrap.plot = True

# prepare and run
P = Ptycho(p,level=2)



# ---------------------------------------------------------
# Load data, some metadata, geometry object, and do some initial plotting
#-----------------------------------------------------------



#%%

import matplotlib
matplotlib.use( 'Qt5agg' )

# Prepare for reconstruction


algorithm = 'PIE'


# load the data from P to run the code for 3d PIE
#import pdb; pdb.set_trace()
probe = np.squeeze(P.probe.storages['Sscan01'].data )
diff3 = P.diff.storages['S0000'].data*P.mask.storages['S0000'].data

#turn all patterns up side down
##diff3 = diff3[:,:,::-1] #correct?


plt.figure()
plt.title(str(p.scans.scan01.data.orientation))
plt.imshow(np.log10(sum(sum(diff3))),cmap='magma')
plt.show()
obj_storage = P.obj.storages['Sscan01']

#plt.figure()
#plt.title('(True, False, False)')
#plt.imshow(abs(np.squeeze(obj_storage.data)[20]),cmap='magma')

g = P.pods['P0000'].geometry


# Keep a copy of the object storage, and fill the actual one with an
# initial guess (like zeros everywhere).
S_true = obj_storage.copy(owner=P.obj)

S_true_cart = g.coordinate_shift(S_true, input_space='real', input_system='natural', keep_dims=False)
x, z, y = S_true_cart.grids()
xcut = int(S_true_cart.data.shape[1]/2)
zcut = int(S_true_cart.data.shape[2]/2)
ycut = int(S_true_cart.data.shape[3]/2)

fact = 1E6
##%%
## Control directions of motorpositions (control the order of diffraction patterns)
#obj_storage.fill(0.0)
#plt.ion()
#diff_frames=[]
#plt.ioff()
#for j, obj_view in enumerate(P.obj.views):
#    
#    plt.figure()
#    plt.imshow( abs(P.obj.views[obj_view].data[20]), cmap='RdBu_r' )
#    plt.title('pos ' + "{0:03}".format(j) )
#    plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\scanning_positions\\' + "{0:03}".format(j)) 
#    print(j)
#    

#%%

# if you dont want to start over and just continue to iteratie, just comment this part out



iterations = 10

save_recon = True
print('Initialize/reset ''manual'' reconstruction script ')
# reset obj. storage
obj_storage.fill(0.0)
storage_save = []
probes = []

#_----------------------------
# for DM!
    # we also need a constant normalization storage, which contains the
    # denominator of the DM object update equation.
    #"S=obj_storage"
    # cut htis out to be able to copy the object (cous i could not fill the copied storage
    # with the probe using views (cos it did not have views)))

#for obj_view in obj_storage.views:
#    obj_storage[obj_view] += np.abs(probe)**2
#
#Snorm = obj_storage.copy(owner=P.obj)
#
#obj_storage.fill(0.0)
#_---------------------------
#%%
# rotate diff patterns transpose images in 4D array

###diff3 = np.transpose(diff3[:,:])
###plt.imshow(sum(sum(diff3))) # why does it look like that?? its a row tocking curve, like we did in the beamtime!

#diff3 = np.transpose(diff3, [0, 1, 3, 2])
#
#plt.figure()
#plt.imshow(np.log10(sum(sum(diff3))) )

# # for flipping ud
#diff3 = diff3[:,:,::-1]
#plt.imshow(sum(sum(diff3))) 
## flipping lr (choose the axis to flip)
#diff3 = diff3[:,:,:,::-1]
#plt.imshow(sum(sum(diff3))) 

#%%
print('Start ''manual'' reconstruction script ')
date_str = time.strftime("%Y%m%d_%H%M")


# Here's a PIE/cPIE implementation
if algorithm == 'PIE':
    beta = 0.1
    eps = 1e-3
    
    ferrors = []

    for i in range(iterations):
        
        ferrors_ = []
        
        j = 0       
        
        # iterate over all views
        for obj_view in P.obj.views:   
        #for j in range(len(P.obj.views)):

            exit_ = P.obj.views[obj_view].data * probe
            prop = g.propagator.fw(exit_)
            ferrors_.append(abs(np.abs(prop)**2 - diff3[j]))
            #import pdb; pdb.set_trace()
            prop[:] = np.sqrt(diff3[j]) * np.exp(1j * np.angle(prop))
            exit = g.propagator.bw(prop)
            # ePIE scaling (Maiden2009)
            #SH: det här är väl  J.M. Rodenburg and H.M.L Faulkner 2004, jag skrev exakt såhär baserat på det peket.
            P.obj.views[obj_view].data += beta * np.conj(probe) / (np.abs(probe).max())**2 * (exit - exit_)
            # PIE and cPIE scaling (Rodenburg2004 and Godard2011b)
            #P.obj.views[obj_view].data += beta * np.abs(probe) / np.abs(probe).max() * np.conj(probe) / (np.abs(probe)**2 + eps) * (exit - exit_)
            
            #update probe function
            #probe = probe + beta *(exit -exit_) * np.conj(P.obj.views[obj_view].data)/ (np.max(abs(P.obj.views[obj_view].data))**2)
            
            #slice till någonting
            #slice1=slice(0,5)
            #slice2=slice(11,17)
            #slice till 41 rotationer test
#            slice1=slice(0,15)
#            slice2=slice(25,41)
#            obj_storage.data[0][slice1] = 1E-9 * np.exp(1j*np.angle(obj_storage.data[0][slice1]))
#            obj_storage.data[0][slice2] = 1E-9 * np.exp(1j*np.angle( obj_storage.data[0][slice2]))

            #slice till np.arange(429,490)
            # apply constraint to 10 pixlar = 240 nm (á 24nm per pix)
            #obj_storage.data[0][0:20-6] = 1E-9 * np.exp(1j*np.angle(obj_storage.data[0][0:20-6]))
            #obj_storage.data[0][20+6:] = 1E-9 * np.exp(1j*np.angle( obj_storage.data[0][20+6:]))


            j += 1
            

        ferrors.append(np.mean(ferrors_))
        
        #storage_save.append(obj_storage.copy(owner=P.obj))
        if not (i % 1):
            #save arrays instread of plottting to improve speed
             storage_save.append(obj_storage.copy(owner=P.obj))
             probes.append(probe)
             
             #TODO
             ##save the last object result
             #if save_recon == True:
             #    np.save(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons_ePIE\dumps\storage_save%d'%i, probe)

        print(i)

          
#%%
        
# PLot!

save = True

date_str = time.strftime("%Y%m%d_%H%M")

savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons_ePIE' + '\\' + date_str
import os
if not os.path.exists(savepath) and save == True:
    os.makedirs(savepath)
    print('new folder in this savepath was created')
    
extent_xy = [fact*x.min(), fact*x.max(), fact*y.min(), fact*y.max()]
extent_xz=  [fact*x.min(), fact*x.max(), fact*z.min(), fact*z.max()]
extent_zy = [fact*z.min(), fact*z.max(), fact*y.min(), fact*y.max()]

    
    
    
for ii, store in enumerate(storage_save):

 

    ###store = storage_save[10]
    S_cart = g.coordinate_shift(store, input_space='real', input_system='natural', keep_dims=False)
    x, z, y = S_cart.grids()
    
    
    fig = plt.figure()    # row col  row col
    ax0 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
    ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=1, rowspan=2)
    ax3 = plt.subplot2grid((3, 3), (1, 1), colspan=1, rowspan=2)
    ax4 = plt.subplot2grid((3, 3), (1, 2), colspan=1, rowspan=2)

    ax0.plot(ferrors,'blue',marker='.', label='Fourier error')
    ax0.legend(bbox_to_anchor=(0.65, 0.5), loc='center left')
    obj = np.squeeze(S_cart.data)
        
    im1 = ax2.imshow( np.abs(obj[:,zcut]).T, extent=extent_xy, interpolation='none', cmap='RdBu_r')
    plt.setp(ax2, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='top view')
    plt.colorbar(im1,ax=ax2)   
    im2 = ax3.imshow( np.abs(obj[:,:,ycut]).T, extent=extent_xz, interpolation='none', cmap='RdBu_r')#,vmin=min_scatt, vmax=max_scatt)
    plt.setp(ax3, ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='side view')
    plt.colorbar(im2,ax=ax3)   
    im3 = ax4.imshow( np.abs(obj[xcut]).T , extent=extent_zy, interpolation='none', cmap='RdBu_r')#,vmin=min_scatt, vmax=max_scatt)
    plt.setp(ax4, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
    plt.colorbar(im3,ax=ax4)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=70)
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=70)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    if save == True:
        plt.savefig(savepath + '\\amp_iter%d_'%(ii))
    
    
    #thight layout should minimize the overlab of labels
    plt.tight_layout() 
#    plt.draw()
#    plt.show()

    mask = np.zeros((S_cart.shape))
    mask[abs(S_cart.data)>0.005] = 1 
    mask = np.squeeze(mask)
    
    fig, ax = plt.subplots(ncols=3)
    plt.suptitle('Phase central cuts.')
    im1 = ax[0].imshow(( mask[:,zcut] * np.angle(np.squeeze(S_cart.data)[:,zcut]) ).T, extent=[fact*x.min(), fact*x.max(), fact*y.min(), fact*y.max()], interpolation='none',  cmap='jet')
    plt.setp(ax[0], ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='top view')
    fig.colorbar(im1,ax=ax[0] )    
    im2 = ax[1].imshow(( mask[:,:,ycut] * np.angle(np.squeeze(S_cart.data[0])[:,:,ycut])).T, extent=[fact*x.min(), fact*x.max(), fact*z.min(), fact*z.max()], interpolation='none', cmap='jet')
    plt.setp(ax[1], ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='side view')
    fig.colorbar(im2,ax=ax[1] )    
    im3 = ax[2].imshow(( mask[xcut] *  np.angle(np.squeeze(S_cart.data[0])[xcut])).T, extent=[fact* z.min(), fact*z.max(), fact*y.min(), fact*y.max()], interpolation='none', cmap='jet')
    plt.setp(ax[2], ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
    fig.colorbar(im3,ax=ax[2] )
    
    
    plt.tight_layout()
#    plt.draw()
    if save == True:
        plt.savefig(savepath + '\\phase_iter%d_'%(ii))
#    plt.show()
 
    
    # this is not correct.. its == Q*u()
    # this doesn make sense... but it makes the most sense if I want to compare to the XRD maps
    
    S_cart.data.shape
    z_middle = int(S_cart.data.shape[2]/2)
    y_middle = int(S_cart.data.shape[3]/2)
    slicez = slice(z_middle-150,z_middle+150)
    slicey = slice(y_middle-80,y_middle+80)       

    #Average over the reigion of interest (that is a bit arbitrary)
    #mean_phase_x = np.mean(np.angle(np.squeeze(S_cart.data)[xcut][slicex, slicey]))
    #TODO not plotting strain now
    #strain_xcut = np.angle(np.squeeze(S_cart.data)[xcut][slicex, slicey])
    #strain_xcut = 100*(np.angle(np.squeeze(S_cart.data)[xcut][slicex, slicey]) - mean_phase_x) / mean_phase_x
    
    ###extent_zy_sliced = extent_zy = [fact*z.min(), fact*z.max(), fact*y.min(), fact*y.max()][slicez, slicey]
    
    fig, ax = plt.subplots(ncols=1)  
    im3 = ax.imshow((np.abs(np.squeeze(S_cart.data)[xcut][slicez, slicey]) ).T,  interpolation='none', cmap='RdBu_r')
    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
    fig.colorbar(im3,ax=ax ); plt.show()
    if save == True:
        plt.savefig(savepath + '\\amp_zoomed_iter%d_'%(ii))
    
    fig, ax = plt.subplots(ncols=1)  
    im3 = ax.imshow((mask[xcut][slicez, slicey] * np.angle(obj[xcut][slicez, slicey]) ).T,  interpolation='none', cmap='jet')
    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
    fig.colorbar(im3,ax=ax ); plt.show()
    if save == True:
        plt.savefig(savepath + '\\phase_zoomed_iter%d_'%(ii))
    
    InP_Qvect = 18543793660.27452
    #TODO correcto? z is the second axis. and resolution is the pixel size in measuremnet grid
    dz = g.resolution[1]   
    strain = np.gradient((((np.angle(obj[xcut][slicez, slicey])))/InP_Qvect) , dz)[0]  
    masked_strain = mask[xcut][slicez, slicey]* strain

    #plot strain
    vs_max = 0.75
    vs_min =-0.75
    fig, ax = plt.subplots(ncols=1)  
    im3 = ax.imshow(np.flipud(masked_strain *100).T, cmap='RdBu_r',  interpolation='none')
    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='Strain [%] Central cut of 3D recon')
    fig.colorbar(im3,ax=ax ); plt.show()
    if save == True:
        plt.savefig(savepath + '\\strain_zoomed_iter%d_'%(ii))
        
    
#    # mean value to compare to XRD
#    strain_ = np.mean(np.angle(np.squeeze(S_cart.data)[19-3:19+3][slicex, slicey]), axis = 0)
#    fig, ax = plt.subplots(ncols=1)  
#    im3 = ax.imshow((np.mean(np.abs(np.squeeze(S_cart.data)[19-3:19+3][slicex, slicey]), axis = 0) * strain_).T, extent=[fact* z.min(), fact*z.max(), fact*y.min(), fact*y.max()], interpolation='none', origin='lower', cmap='jet')
#    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='mean value phase \n front view')
#    fig.colorbar(im3,ax=ax )
#    
#    #mean value in amplitude
#    mean_amp = np.sum(np.abs(np.squeeze(S_cart.data)[19-3:19+3]), axis = 0)
#    fig, ax = plt.subplots(ncols=1)  
#    im3 = ax.imshow(mean_amp.T, extent=[fact* z.min(), fact*z.max(), fact*y.min(), fact*y.max()], interpolation='none', origin='lower', cmap='jet')
#    plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='mean value amplitude \n front view')
#    fig.colorbar(im3,ax=ax )

    fig, ax = plt.subplots(ncols=2)
    proj= 0#int(len(probe)/2)
    
    im1 = ax[0].imshow( np.abs(probes[ii][proj]), interpolation='none', cmap='jet')#, extent=[fact*x.min(), fact*x.max(), fact*y.min(), fact*y.max()])
    plt.setp(ax[0], ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='probe amplitude')
    fig.colorbar(im1,ax=ax[0] )
    
    im2 = ax[1].imshow( np.angle(probes[ii][proj]), interpolation='none', cmap='jet')#, extent=[fact*x.min(), fact*x.max(), fact*z.min(), fact*z.max()])
    plt.setp(ax[1], ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='probe phase')
    fig.colorbar(im2,ax=ax[1] )
    plt.show()
    
    
    if save == True:
        plt.savefig(savepath + '\\probe_iter%d_'%(ii))

if save == True:
    print('Plots saved at: ', savepath)

    







