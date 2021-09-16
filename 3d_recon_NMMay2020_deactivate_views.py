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
from ptypy.experiment.nanomax3d import NanomaxBraggMay2020_shifting
import time
date_str = time.strftime("%Y%m%d_%H%M")

p = u.Param()
p.run = 'recon' + date_str    #test_deactivate_views' # 'recon' + date_str



#original scans full range USE THIS WITH CURRENT SHIFTING LIST
#scans = np.arange(429,490).tolist()

########################OBS you have to start on 429 for the shifting to work


# InP range
#scans = np.arange(429,470).tolist()

# InP tighter range
scans = np.arange(429,445+1).tolist()

# InP even tighter
#scans = np.arange(434,457+1).tolist()


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
p.scans.scan01.data.xMotor = 'npoint_buff/x' #'sx' for stepscan  
p.scans.scan01.data.yMotor = 'npoint_buff/y' #'sy'


# Diffraction angle (theta, not two theta) in degrees
p.scans.scan01.data.theta_bragg = 8.8 #10.54 # for InP according from table value. If that is correct, I dont need 2theta.       #11.5 #inP? 


# pick an even number!
p.scans.scan01.data.shape = 170# 181  # 184 includes gap.  170# 512# 150
#center of RAW image
p.scans.scan01.data.cen = (148, 345)    # before: InP (148,364) GaInPx155
p.scans.scan01.data.center = None
p.scans.scan01.data.auto_center = False


## Alternatively, a 3-tuple of booleans may be provided (do_transpose, do_flipud, do_fliplr)
p.scans.scan01.data.orientation = (True, False, False)

# shifting vectors for the full range original scans
p.scans.scan01.data.vertical_shift =   np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\vertical_shift_vector.npy').tolist()
p.scans.scan01.data.horizontal_shift = np.load(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\horizontal_shift_vector.npy').tolist()                               

#p.scans.scan01.data.load_parallel = None # 'all'
p.scans.scan01.data.psize = 55e-6
p.scans.scan01.data.energy = 10.0
p.scans.scan01.data.distance = 1.0
p.scans.scan01.data.I0 = 'alba2/1'
p.scans.scan01.data.min_frames = 10
p.scans.scan01.data.load_parallel = 'all'

p.scans.scan01.illumination = u.Param()
p.scans.scan01.illumination.model = 'recon' 
p.scans.scan01.illumination.recon = u.Param()
p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan14\scan14_ML_0200.ptyr'
#p.scans.scan01.illumination.recon.rfile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\beamtime_folder\process\ptycho\recons\scan13\scan13_ML_0200.ptyr'
p.scans.scan01.illumination.aperture = None 

##p.scans.scan01.illumination.aperture = u.Param() 
##p.scans.scan01.illumination.aperture.form = 'circ'
##p.scans.scan01.illumination.aperture.size = 50e-9 

##p.scans.scan01.sample = u.Param()
##p.scans.scan01.sample.fill = 1E-20# 1e-3 ########################OBS this fields the whole of z, not just were the views are active. maybe try to add this after level 2

p.engines = u.Param()
p.engines.engine00 = u.Param()
p.engines.engine00.name = 'ePIE'    
p.engines.engine00.numiter = 100
p.engines.engine00.probe_update_start = 100000
p.engines.engine00.numiter_contiguous = 1

#p.engines.engine00.probe_support = None   #get an error without these
# sample support option if in DM 3d_Bragg engine
#p.engines.engine00.sample_support = None #use with DM_3dBragg only

##p.engines.engine00.sample_support = u.Param()
###p.engines.engine00.sample_support = True dont need this?
##p.engines.engine00.sample_support.type = 'thinlayer'
##p.engines.engine00.sample_support.size = 200E-9   #diameter. No support for long axis. 
##p.engines.engine00.sample_support.coefficient = 0.0
##p.engines.engine00.sample_support.shrinkwrap = None

# need to fix which coordinates to wrap around in order to use shrinkwrap
# p.engines.engine00.sample_support.shrinkwrap = u.Param()
# p.engines.engine00.sample_support.shrinkwrap.cutoff = .3
# p.engines.engine00.sample_support.shrinkwrap.smooth = None
# p.engines.engine00.sample_support.shrinkwrap.start = 15
# p.engines.engine00.sample_support.shrinkwrap.plot = True

# prepare and run
P = Ptycho(p,level=2)



#____________________________________________________________________________________________________


# turn of some of the views and make a  reconstruction of a smaller region
# choose view to turn off:
turned_off_views_ID = []
############views_off = np.concatenate( [ np.arange(0,15) ,   np.arange(26,75) ,
############                              np.arange(86,135), np.arange(146,195) ,
############                              np.arange(206,255), np.arange(266,315),
############                              np.arange(326,375), np.arange(386,435), np.arange(446,480) ]).tolist()
############

# only views on InP segment on
views_off = np.concatenate( [ np.arange(0,35) ,   np.arange(46,95) ,
                              np.arange(106,155), np.arange(166,215) ,
                              np.arange(226,275), np.arange(286,335),
                              np.arange(346,395), np.arange(406,455), np.arange(466,480) ]).tolist()




for inte in views_off: 
    turned_off_views_ID.append('V'+ str('{0:04}'.format(inte)))
    
# turn of object views    
for sID, s in P.obj.S.items():
    for v in s.views:
        if v.ID in turned_off_views_ID:   #edit not in
            v.active = False

# NEed this?

## turn of diff views
for sID, s in P.diff.S.items():
    for v in s.views:
        if v.ID in turned_off_views_ID:   #edit not in    
            v.active = False



#P.run()

#____________________________________________________________________________________________________
import matplotlib
matplotlib.use( 'Qt5agg' )
#import ipdb; pdb.set_trace()


import logging
logging.getLogger('matplotlib.font_manager').disabled = True
diff_data = P.diff.storages['S0000'].data*P.mask.storages['S0000'].data


# for np arrays in np format. Save outside document folder, otherwise backup will screw up
#np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\four_rot_around_170nm_InP_peak', diff_data)
#np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_data_170nm', diff_data)
#print(r'Save diff data to C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_data_170nm.np')

#open this file
#loaded_array = np.load('...\diff_data_170nm.npy')


# save this script in the recon folder
#from shutil import copy
#src = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\scripts\3d_recon_NMMay2020_deactivate_views.py'
#dst = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons' + "\\" + p.run
#copy(src, dst)





plt.figure()
plt.imshow(np.log10(sum(sum(diff_data))),cmap='magma', interpolation='none')
plt.title('Summed intensity')
plt.colorbar()
plt.savefig('ggg2')
plt.show()



