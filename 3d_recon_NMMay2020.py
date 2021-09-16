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
#scans = np.arange(429,460).tolist()


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

#need this to save the dump I think
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
p.scans.scan01.data.path = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'
p.scans.scan01.data.maskfile ='C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/merlin_mask_200430_8keV.h5'

p.scans.scan01.data.scans = scans
p.scans.scan01.data.xMotor = 'npoint_buff/x' 
p.scans.scan01.data.yMotor = 'npoint_buff/y' #'sy' (använde sy fram till nu 20210610)
#p.scans.scan01.data.xMotorFlipped = True

# Diffraction angle (theta, not two theta) in degrees
p.scans.scan01.data.theta_bragg = 8.8 #10.54 # for InP according from table value. If that is correct, I dont need 2theta. 

p.scans.scan01.data.shape = 170 
#center of RAW image
p.scans.scan01.data.cen = (148, 345)    # before: InP (148,364) GaInP:155
p.scans.scan01.data.center = None
p.scans.scan01.data.auto_center = False

##Alternatively, a 3-tuple of booleans may be provided (do_transpose, do_flipud, do_fliplr)
p.scans.scan01.data.orientation = (True, False, False)

# shifting vectors for the full range original scans
p.scans.scan01.data.vertical_shift =   np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\vertical_shift_vector.npy').tolist()
p.scans.scan01.data.horizontal_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\horizontal_shift_vector.npy').tolist()                               

#p.scans.scan01.data.load_parallel = None # 'all'
p.scans.scan01.data.psize = 55e-6
p.scans.scan01.data.energy = 10.0
p.scans.scan01.data.distance = 1.0
p.scans.scan01.data.I0 = 'alba2/1'

p.scans.scan01.illumination = u.Param()

p.scans.scan01.illumination.model = 'recon' 
p.scans.scan01.illumination.recon = u.Param()
p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan14\scan14_ML_0200.ptyr'
#p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan13\scan13_ML_0200.ptyr'
p.scans.scan01.illumination.aperture = None 

##p.scans.scan01.illumination.aperture = u.Param() 
##p.scans.scan01.illumination.aperture.form = 'circ'
##p.scans.scan01.illumination.aperture.size = 300e-9 

p.scans.scan01.sample = u.Param()
p.scans.scan01.sample.fill = 1E-20

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
P = Ptycho(p,level=5)





import matplotlib
matplotlib.use( 'Qt5agg' )

import logging
logging.getLogger('matplotlib.font_manager').disabled = True
diff_data = P.diff.storages['S0000'].data * P.mask.storages['S0000'].data 

# for np arrays in np format. Save outside document folder, otherwise backup will screw up
#np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\four_rot_around_170nm_InP_peak', diff_data)
np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_data', diff_data)

print(r'Save diff data to C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_data.np')
#open this file
#loaded_array = np.load('C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\scripts\diff_data.npy')


# save this script in the recon folder
from shutil import copy
src = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\scripts\3d_recon_NMMay2020.py'
dst = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons' + "\\" + p.run
#dst = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons' + "\\" + p.run + "\\" + 'recon_script_' + p.run + '.py'
copy(src, dst)



plt.figure()
plt.imshow(np.log10(sum(sum(diff_data))),cmap='magma', interpolation='none')

plt.title('Summed intensity for 3 rotations (log)_orientation (False, False, False)')
plt.colorbar()
plt.savefig('ggg2')
plt.show()



