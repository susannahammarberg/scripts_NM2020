"""
 Reconstruction script for NanoMAX MAy 2020
 flyscans
 
 shifting implemented
 IO normalization not implemented
 
 
"""

import ptypy
from ptypy.core import Ptycho
from ptypy import utils as u
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from ptypy.experiment.nanomax3d import NanomaxBraggMay2020_shifting # after update need to update spec ptyScan class

p = u.Param()
p.run = 'new_roi_real_probe_and_apperture'    #shifted_support'

#scans = [444,445,446]#[np.arange(435,439)] or .tolist()#453)#[444,445,446]
#original scans subset
#scans = np.arange(434,454).tolist()    #456 #inP range
#scans = np.arange(434,437).tolist()


#original scans full range USE THIS WITH CURRENT SHIFTING LIST
scans = np.arange(429,490).tolist()

# the most diffracting scans of the 170nm segment
#scans = np.arange(438,441+1).tolist()


#replacement scans
#scans = np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ]).tolist()
p.data_type = "single"   #or "double"
p.verbose_level = 3

# use special plot layout for 3d data  (but the io.home part tells ptypy where to save recons and dumps)
p.io = u.Param()
p.io.home = './'
p.io.home = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/3drecons/'
#p.io.autosave = u.Param()
#p.io.autosave.interval = 1 # does not work
#p.io.autoplot = u.Param()
#p.io.autoplot.layout = 'bragg3d'
#p.io.autoplot.dump = True
#p.io.autoplot.interval = 1
# TODO make it create the plots


p.scans = u.Param()
p.scans.scan01 = u.Param()
p.scans.scan01.name = 'Bragg3dModel'
p.scans.scan01.data = u.Param()
p.scans.scan01.data.name = 'NanomaxBraggMay2020_shifting'
p.scans.scan01.data.path = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'
p.scans.scan01.data.maskfile ='C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/merlin_mask_200430_8keV.h5'

p.scans.scan01.data.scans = scans
p.scans.scan01.data.xMotor = 'npoint_buff/x' #'sx' for stepscan  
p.scans.scan01.data.yMotor = 'npoint_buff/y' #'sy'


#TODO Diffraction angle (theta, not two theta) in degrees
p.scans.scan01.data.theta_bragg = 11.5 #inP? 


p.scans.scan01.data.shape = 170# 512# 150
#center of RAW image
p.scans.scan01.data.cen = (148, 345)    # before: InP (148,364)
p.scans.scan01.data.center = None #(4,4)# InPx364 GaInPx155 #None # auto, you can also set (i, j) center here.

# was commented bforeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
p.scans.scan01.data.auto_center = False  
p.scans.scan01.data.orientation = 4+1+2

# shifting vectors for the full range original scans                                       
                                         #att öka flyttar ned tråden. 
p.scans.scan01.data.vertical_shift = [ 0,  1,  1,  1, 0, -1,  0,  0,  1,  2,  2,  2,  1,  1,  2,  2,  2,  #445
        1,  1,  1,  1,  2,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  3,  3,
        3,  2,  0,  1,  1,  1,  1,  2,  1,  1,  1,  0,  1,  0,  0, -1, -1,
        0,  0,  0,  0,  0, -1, -1,  0,  0,  0]


#"np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\xrf\original_scans\In_L\vertical_shift_vector.npy').tolist()
#p.scans.scan01.data.horizontal_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\xrf\original_scans\In_L\horizontal_shift_vector.npy').tolist()

                               #  minska för att flytta rl till vänster. kolla om                     #450
                              # =öka för att flytta tråd ttill höger?
                               #från 2 till 6 (tråd flyttades till vänster)
                              #öka så flyttas tråd till vänster
p.scans.scan01.data.horizontal_shift =[0, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 3, 2, 3, 3, 4, 4, 5, 5, 6,
                                                                     #472
       6, 6, 7, 6, 6, 6, 6, 5, 4, 4, 4, 5, 4, 5, 4, 5, 4, 5, 4, 4, 4, 4,
                                                      #489
       5, 5, 4, 4, 4, 4, 3, 4, 3, 3, 3, 4, 4, 3, 3, 3, 2]



#p.scans.scan01.data.load_parallel = None # 'all'
p.scans.scan01.data.psize = 55e-6
p.scans.scan01.data.energy = 10.0
p.scans.scan01.data.distance = 1.0

#p.scans.scan01.illumination = u.Param()
##p.scans.scan01.illumination.diversity = u.Param() #not working but default is none#None #to avoid rescaling issues with probe
#p.scans.scan01.illumination.model = 'recon'
#p.scans.scan01.illumination.recon = u.Param()
##p.scans.scan01.illumination.recon.label = 'scan10probe'
#p.scans.scan01.illumination.recon.rfile = 'C:/Users/Sanna/Documents/beamtime/NanoMAX062017/Analysis_ptypy/nice_probe_ptyrfiles/scan10/scan10_pilatus_ML.ptyr'
#p.scans.scan01.illumination.photons = None

#ID = "S00G00"
#p.scans.scan01.illumination.recon.layer = None



#p.scans.scan01.illumination = u.Param()
#loaded_profile = np.load(r'C:\Users\Susanna\Documents\GitHub\simulated_nanodiffraction\probe10_focus.npy')
##probe1 = np.rot90(np.copy(loaded_profile)[1:121,0:120],3)
#probe1 = np.rot90(np.copy(loaded_profile),3)
#psize=[ 1.89051824e-08,   1.85578409e-08]
##probe1 = np.load(r'F:\Susanna\ptypy_scripts\probe_191412_ePIE200_psize73nm.npy')
##psize = 7.262429299437677e-08




#Cprobe = classes.Container(data_dims=2, data_type='complex')
#Sprobe = Cprobe.new_storage(psize=psize, shape=(probe1.shape[-1]))
#Sprobe.fill(probe1)


#decomment this to use real probe
p.scans.scan01.illumination = u.Param()
p.scans.scan01.illumination.model = 'recon'#Sprobe
p.scans.scan01.illumination.recon = u.Param()
p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan14\scan14_ML_0200.ptyr'
p.scans.scan01.illumination.aperture = None 

p.scans.scan01.illumination.aperture = u.Param() 
p.scans.scan01.illumination.aperture.form = 'circ'
p.scans.scan01.illumination.aperture.size = 500e-9 

p.scans.scan01.sample = u.Param()
p.scans.scan01.sample.fill = 1e-3



p.engines = u.Param()
p.engines.engine00 = u.Param()
p.engines.engine00.name = 'DM_3dBragg'    
p.engines.engine00.numiter = 1# 20
p.engines.engine00.probe_update_start = 100000
p.engines.engine00.numiter_contiguous = 1 #5 # think this is how often you want to save
##TODO look at these!
p.engines.engine00.probe_support = None   #get an error without these 
#p.engines.engine00.sample_support = None #use with DM_3dBragg


p.engines.engine00.sample_support = u.Param()

##p.engines.engine00.sample_support = True dont need this?
p.engines.engine00.sample_support.type = 'rod_3d'
p.engines.engine00.sample_support.size = 230e-9
p.engines.engine00.sample_support.coefficient = 0.0
p.engines.engine00.sample_support.shrinkwrap = None

# need to fix which coordinates to wrap around in order to use shrinkwrap
# p.engines.engine00.sample_support.shrinkwrap = u.Param()
# p.engines.engine00.sample_support.shrinkwrap.cutoff = .3
# p.engines.engine00.sample_support.shrinkwrap.smooth = None
# p.engines.engine00.sample_support.shrinkwrap.start = 15
# p.engines.engine00.sample_support.shrinkwrap.plot = True

# prepare and run
P = Ptycho(p,level=5)

import matplotlib
matplotlib.use( 'Qt5agg' )

import logging
logging.getLogger('matplotlib.font_manager').disabled = True
diff_data = P.diff.storages['S0000'].data*P.mask.storages['S0000'].data 

# for np arrays in np format. Save outside document folder, otherwise backup will screw up
#np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\four_rot_around_170nm_InP_peak', diff_data)
np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_data', diff_data)

print('Save diff data to C:_Users_Sanna_NanoMAX_May2020_nobackup__diff_data.np')
#open this file
#loaded_array = np.load('C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\scripts\diff_data.npy')



plt.figure()
plt.imshow(np.log10(sum(sum(diff_data))),cmap='magma', interpolation='none')
plt.title('Summed intensity')
plt.colorbar()
plt.savefig('ggg2')
plt.show()
