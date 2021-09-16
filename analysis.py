"""
Author: Sanna

Plot diffraction data and make minor aaysis
after using ptypy loading class to load and prepare for ptychography reconstruction

"""

import numpy as np
import matplotlib.pyplot as plt

#scans = np.arange(429,490).tolist()


scans = np.arange(429,470).tolist()


#open this file
data = np.load(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\diff_data.npy')

#
#plot all diffraction
#################################################


plt.figure()
plt.imshow((np.log(sum(sum(data)))),cmap='magma', interpolation='none')
plt.title('Summed intensity for 3 rotations (log)')
plt.colorbar()
plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\log_61rotations')
plt.show()

#######################################################

# plot along NW long axis
plt.ioff()
for ii in range(130,165):
    
    plt.figure()
    plt.imshow(np.log10(sum(data[ii])))
    plt.title('log10 pos%d'%ii)
    plt.colorbar()
    plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\diffraction_along_long_axis\row5\pos%d'%ii)
    plt.show()


#gather motorpositions and plot them with BF map.

# function for calculating bright field of data
# input is the data as a 3D array with dimension data[scans_y * scans_x][yPixels][xPixels]
def bright_field(data,x,y):
    index = 0
    photons = np.zeros((y,x)) 
    for row in range(0,y):
        for col in range(0,x):
            photons[row,col] = sum(sum(data[index])) #/ max_intensity
            index = index+1
            
    return photons


# make BF maps
###############################################################################################
###import pdb; pdb.set_trace()
##ii=0
##for sc in scans:
##    
##    #bf_map = bright_field(data[:,ii],67,13) (all data (unshifted))
##    bf_map = bright_field(data[:,ii],60,8)
##    if sc == 429:
##        first_map = bf_map
##    fig, ax = plt.subplots(nrows=2)
##    ax[0].imshow(np.log10(first_map),cmap='jet')
##    ax[0].axvline(x=51.5,color='red')
##    ax[0].axhline(y=2,color='red')
##    ax[1].imshow(np.log10(bf_map),cmap='jet')
##    ax[1].axvline(x=51.5,color='red')
##    ax[1].axhline(y=2,color='red')
##    plt.setp(ax[0], ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]')
##    
##    #plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\BF_maps\after_shifting\S#%d'%sc)
##    plt.show()
##    
##    ii += 1
###print('saved BF images')
##################################################################################################
# make single BF image
ii = 15
#bf_map = bright_field(data[:,ii],60,8)

#________________________________________________________________________________________________
### Make single MF image from selected positions
views_off = np.concatenate( [ np.arange(0,15) ,   np.arange(26,75) ,
                              np.arange(86,135), np.arange(146,195) ,
                              np.arange(206,255), np.arange(266,315),
                              np.arange(326,375), np.arange(386,435), np.arange(446,480) ]).tolist()

##views_off = np.concatenate( [ np.arange(0,35) ,   np.arange(46,95) ,
##                              np.arange(106,155), np.arange(166,215) ,
##                              np.arange(226,275), np.arange(286,335),
##                              np.arange(346,395), np.arange(406,455), np.arange(466,480) ]).tolist()

for view in views_off:
    data[view,ii] = 0

bf_map = bright_field(data[:,ii],60,8)

##__________________________________________________________________________________________________


#runCommand('%npointflyscan sx 13.5 16.8 66 sy 2 2.6 12 .1 .003')

# plot single BF
fig, ax = plt.subplots(nrows=1)
ax.imshow((np.fliplr(bf_map)),cmap='jet', interpolation='none', origin='lower', extent= [ -50E-3*30, 50E-3*30,  -50E-3*4, 50E-3*4 ])
#ax.axvline(x=51.5,color='red')
#ax.axhline(y=2,color='red')
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='BF #S %d'%scans[ii])
plt.show()
    
#################################################################################################
#import pdb; pdb.set_trace()

# save a bragg CDI peak
###############################
##print(data.shape)
##print(data[151].shape)
### print the CDI peak
###pos=[150, 151, 152, 153, 154] # around 80 nm segment
##pos = [154, 155, 156, 157, 158, 159, 160, 161, 162, 163] # around 170 nm segment (160 seeems like the best)
##
##for d2 in pos:
##    plt.figure()
##    plt.imshow(sum(data[d2]),cmap='jet'); plt.title(d2); plt.colorbar()
 ##   #plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\diffraction_peaks_170nm_segment\total_diffraction_jet_position%d'%d2)
   ## #plt.show()

# for np arrays in np format. Save outside document folder, otherwise backup will screw up
#np.save(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\bragg_peak_80segment', data[151])
#print('saved CDI data to C:_Users_Sanna_NanoMAX_May2020_nobackup_bragg_peak_80segment')

########################################################################################################
###open this file
#loaded_array = np.load(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\bragg_peak_80segment.npy')
##loaded_array = np.load(r'C:\Users\Sanna\NanoMAX_May2020_nobackup\four_rot_around_170nm_InP_peak.npy')
##
##plt.figure()
##plt.imshow(np.log(sum(sum(loaded_array)))); plt.title('4 cetral rotations around InP of 170nm segment log')
##plt.show()
#####################################################################################################

#plot all frames in a rocking curve 80 nm segment

##scans = np.arange(429,490).tolist()
### image of all diffraction in the data
##i=0
##for dat in data[160]:
##    i += 1     
##    plt.figure()
##    plt.imshow(dat, cmap='jet'); plt.colorbar();  plt.title('Rocking curve InP Braggg peak rotation %d'%i)
##    
##    #plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\diffraction_peaks_170nm_segment\position160_rockingcurve\rot%d'%i)
##    #plt.show()
