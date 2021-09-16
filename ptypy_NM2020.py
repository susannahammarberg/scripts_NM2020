"""
 loads data via ptypy and perform xrd analysis from NanoMAX experiment
 
 
 Sanna : not at all updated to 2020! just saved the script copied from nanomax 2017 analysis
"""
import ptypy
from ptypy.core import Ptycho
from ptypy.core.classes import Container, Storage, View
from ptypy import utils as u
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from ptypy.experiment.nanomax3d import NanomaxBraggJune2017 # after update need to update spec ptyScan class
from matplotlib import ticker
import time
date_str = time.strftime("%Y%m%d") # -%H%M%S")

import sys   #to collect system path ( to collect function from another directory)
sys.path.insert(0, 'C:/Users/Sanna/Documents/python_utilities') #can I collect all functions in this folder?
sys.path.insert(0, 'F:/Susanna/Diamond_Sept2018/')
#import sannastools as st
#from movie_maker import movie_maker
import h5py
# change font in plot
from matplotlib import rc

### change the plotting config   
import matplotlib
matplotlib.use( 'Qt5agg' )

p = u.Param()
p.run = 'runrunrun'   

scans = [441,443,443] #range(192, 200+1)+range(205, 222+1) #range(192, 195)# 

p.data_type = "single"   #or "double"
# for verbose output
p.verbose_level = 5

p.scans = u.Param()
p.scans.scan01 = u.Param()
p.scans.scan01.name = 'Bragg3dModel'
p.scans.scan01.data = u.Param()
p.scans.scan01.data.name = 'NanomaxBraggJune2020'

p.scans.scan01.data.datapath ='C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/raw/'

#TODO funkar inte
p.scans.scan01.data.datafile = '%s.h5' % scans
p.scans.scan01.data.detfilepattern = 'scan_%04d_merlin.hdf5'
# not sure if this loads properly
p.scans.scan01.data.maskfile = 'C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/merlin_mask.h5'

p.scans.scan01.data.scans = scans
p.scans.scan01.data.theta_bragg = 10.91#10.91 #  10.91    #InP:10.91  GaInP: 11.09  # homogeneous wires to 11.0


#raw_center = (256,256)
# does not work with asymmetric
######################################
p.scans.scan01.data.shape = (150,150) #(512,512) #   #256 for homogeneous# needs to be an EVEN number for using shifting
p.scans.scan01.data.auto_center = False # 
p.scans.scan01.data.orientation = 4+1+2

# ptypy says: Setting center for ROI from None to [ 75.60081158  86.26238307].   but that must be in the images that iI cut out from the detector
#raw_center = (342,245) #for homogeneous wires
#raw_center = (190,330) #  for InGaP in segmented NWs
raw_center = (190,190) # for InP in segmented NWs

detind0 = int(raw_center[0] - p.scans.scan01.data.shape[0]/2)
detind1 = int(raw_center[0] + p.scans.scan01.data.shape[0]/2)
detind2 = int(raw_center[1] - p.scans.scan01.data.shape[1]/2)
detind3 = int(raw_center[1] + p.scans.scan01.data.shape[1]/2)

p.scans.scan01.data.detector_roi_indices = [detind0,detind1,detind2,detind3]  # this one should not be needed since u have shape and center...
p.scans.scan01.data.center = (raw_center[0] - detind0,raw_center[1] - detind2) # (200,270) #(512-170,245)     #(512-170,245) for 192_   #Seems like its y than x
# tprev used:  [275,425,150,300]
#calculates the center based on the first pic said Alex.  186.75310731  265.64597192] thats not wr

#p.scans.scan01.data.load_parallel = 'all'
p.scans.scan01.data.psize = 55e-6
p.scans.scan01.data.energy = 9.49
# TODO change d
p.scans.scan01.data.distance = 1.149 # =sqrt(1.065**2 + 0.43**2) is the distance along optical axis

# This shifts the entire scan (projection) in real space, in units of steps 
##############################################################################
#S segmented NWs
##############################################################################
# the smallers values here (minus values) will include the 0:th scanning position
##TODO test with no shifting and compare result: [0]*51 #
#TODO more is to do the automatic shifting for these scans!!

# shifting list THIS
#p.scans.scan01.data.vertical_shift =  [-1,-1,0,0,0,  0,0,2,1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3]
#p.scans.scan01.data.horizontal_shift =  [3,2,0,1,2,  3,4,3,4,5,  5,6,6,5,6,  5,4,7,8,8,  8,8,10,11,12,  11,12,12,11,12,  12,11,12,13,13,  14,15,14,14,14,  13,15,16,15,14,  17,19,18,18,17,   17]
#
#



# test smaller roi:
#p.scans.scan01.data.vertical_shift =  [-1,       0,0   ,  1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3]
#p.scans.scan01.data.horizontal_shift =  [0,       0,1,    1,2,  2,3,3,2,3,  2,1,4,5,5,  5,5,7,8,9,  8,9,9,8,9,  9,8,9,10,10,  11,14,11,11,11,  10,12,13,12,11,  14,16,15,15,14,   14]

#the smaller roi which fits the 45 setgment
#p.scans.scan01.data.vertical_shift =  [              0 ,   1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3]
#p.scans.scan01.data.horizontal_shift =  [           1,    1,2,  2,3,3,2,3,  2,1,4,5,5,  5,5,7,8,9,  8,9,9,8,9,  9,8,9,10,10,  11,12,11,11,11,  10,12,13,12,11,  14,15,14,14,13,   13]#,15,14,14,13,   13

# should be 0 or 1?
#p.scans.scan01.data.vertical_shift = [0] * len(scans) 
#p.scans.scan01.data.horizontal_shift = [0] * len(scans) 
            # InGaP = [116:266,234:384]
            # InP = [116:266,80:230]
#p.scans.scan01.data.detector_roi_indices = [116,266,80,230]  #[116,266,80,230]  #    # this one should not be needed since u have shape and center...

##############################################################################
# homogenius NWs
##############################################################################
## for not using shifts:            
#p.scans.scan01.data.vertical_shift = [0] * len(scans) 
#p.scans.scan01.data.horizontal_shift = [0] * len(scans) 

# referens ska vara 0, alla andra -
#open this file
## when it says horizontal i have used vertical because i am an i
## I will just have fo change the names everyw here             
#p.scans.scan01.data.horizontal_shift = list(np.load('C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\ptycho_192_222\\horizontal_shift_vector.npy'))
##vertical_shift_vector
#p.scans.scan01.data.vertical_shift = list(np.load('C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\ptycho_192_222\\vertical_shift_vector.npy'))
#
# TODO is this what I should use, or the list in the prev. row?            
#p.scans.scan01.data.vertical_shift = [ 0, 0, 1, 2, 2, 2, 3, 3, 3, 1, 0, 0, 0, 0, 0, -2, -2, -2, -3, -3, -2, -2,  -3, -2, -3, -3,  -4]
#p.scans.scan01.data.horizontal_shift = [ -8, -8, -8, -8, -8, -3, -2, -4, 0, -5, -3, -6, -3, -6, -3, -5, -5, -7, -4,  -7, -3, -6, -2, -7, -2, -6, -3 ] 

# prepare and run
P = Ptycho(p,level=2)

# ---------------------------------------------------------
# Load data, some metadata, gemetry object, and do some initial plotting
#-----------------------------------------------------------

###   use this in console for plotting      
#%matplotlib inline


#%%
# gather motorpositions from first rotation in scan list for plotting
scan_name_int = scans[2]
scan_name_string = '%d' %scan_name_int 
metadata_directory = p.scans.scan01.data.datapath + sample + '.h5'
metadata = h5py.File( metadata_directory ,'r')
motorpositions_directory = '/entry%s' %scan_name_string  

# gonphi is the rotation of the rotatation stage on which all motors samx/y/z and samsx/y/z and sample rests on.
# when gonphi is orthoganoal to the incoming beam, it is equal to 180. Read in gonphi and 
# get the list of scans numbers and theta (=180-gonphi) and gonphi
def sort_scans_after_theta():
    gonphi_list = []
    # read in all gonphi postins
    for i in range(0,len(scans)):
        scan_name_int = scans[i]
        scan_name_string = '%d' %scan_name_int 
        metadata_directory = p.scans.scan01.data.datapath + sample + '.h5'
        metadata = h5py.File( metadata_directory ,'r')
        motorpositions_directory = '/entry%s' %scan_name_string     
        motorposition_gonphi = np.array(metadata.get(motorpositions_directory + '/measurement/gonphi'))
        gonphi_list.append(motorposition_gonphi[0])

    # order the scan list after gonphi
    # first put them together
    theta_array = 180 - np.array(gonphi_list)
    zipped = list(zip(gonphi_list,theta_array,scans)).sort()
    # then sort after the first col (gonphi)
    #zipped.sort()
    return zipped
scans_sorted_theta = sort_scans_after_theta()  
      
# JW: Define your coordinate system! What are x,y,z, gonphi, and what are their directions
# calculate mean value of dy
# positive samy movement moves sample positive in y (which means the beam moves Down on the sample)
motorpositiony = np.array(metadata.get(motorpositions_directory + '/measurement/samy'))
dy = (motorpositiony[-1] - motorpositiony[0])*1./ (len(motorpositiony) -1)

# calculate mean value of dx
# instead of samx, you find the motorposition in flysca ns from 'adlink_buff' # obs a row of zeros after values in adlinkAI_buff
motorpositionx_AdLink = np.mean( np.array( metadata.get(motorpositions_directory + '/measurement/AdLinkAI_buff')), axis=0)
motorpositionx_AdLink = np.trim_zeros(motorpositionx_AdLink)
dx = (motorpositionx_AdLink[-1] - motorpositionx_AdLink[0])*1./ (len(motorpositionx_AdLink) -1 )

# calc number of rows and cols that is used from measuremtn after realigning postions
# for semgneted Nws you ned the -min,  ithink it is if you have both neg and pos values in your shift vector) Homogeneous NW dont need it
# TODO borde ju va ok att behålla även om min är 0. men måste inehålla 0, men det gör det ju alltid.
nbr_rows = len(motorpositiony) - (np.max(p.scans.scan01.data.vertical_shift) - np.min(p.scans.scan01.data.vertical_shift) )
nbr_cols = len(motorpositionx_AdLink) - (np.max(p.scans.scan01.data.horizontal_shift) - np.min(p.scans.scan01.data.horizontal_shift))

#        Nx = x_mean.shape[1] - (np.max(self.p.horizontal_shift) - np.min(self.p.horizontal_shift))  
#        Ny = x_mean.shape[0] - (np.max(self.p.vertical_shift) - np.min(self.p.vertical_shift)) 
extent_motorpos = [ 0, dx*(nbr_cols-1),0, dy*(nbr_rows-1)]
# load and look at the probe and object
#probe = P.probe.storages.values()[0].data[0]#remember, last index [0] is just for probe  
#obj = P.obj.storages.values()[0].data
# save masked diffraction patterns
# in diff data the rotations are sorted according to gonphi, starting with the LOWEST GONPHI which is the reversed of theta. Hence, diff_data[0] is is theta=12.16 and scan=222 
"                342433243222"

#TODO: remove this diffdata
diff_data = P.diff.storages['S0000'].data*P.mask.storages['S0000'].data #P.diff.storages.values()[0].data*P.mask.storages.values()[0].data[0]
#diff_data = P.diff.storages.values()[0].data[:,:,40:100]*P.mask.storages.values()[0].data[0][:,40:100]
# or load a single pod P.pods['P0896'].diff
#access the grid of the diff storage
#TODO
#a,b,c = P.diff.storages.values()[0].grids()   # but this is in units of meters in the r1r2r3 grid

# load the geometry instance to acces the geometry parameters
# one geometry connected to each POD but it this case it is the same for each pod (since its all from the same measurement, same setup. Everything is the same, psixe energy, distance etc)
g = P.pods['P0901'].geometry
# g. psize real space psize: ( Rocking curve step (in degrees) and pixel sizes (in meters))
# g.resolution: "3D sample pixel size (in meters)." "doc: Refers to the conjugate (natural) coordinate system as (r3, r1, r2)"

# plot the sum of all used diffraction images
plt.figure()
plt.imshow(np.log10(sum(sum(diff_data))),cmap='jet', interpolation='none')
plt.title('Summed intensity ()')
#plt.savefig("summed_intensity") 
#movie_maker(np.log10(abs(probe)))
#to check which axis is which I cur one of the detector images
plt.figure()
plt.title(' data matrix has indices \n diff_data[position,detector width,rotation,detector height]')
plt.imshow(np.log10(sum(sum(diff_data))[60:400][:]),cmap='jet')
plt.ylabel('should be q1 accoring to Berenguer')
plt.xlabel('should be q2')


# plot the sum of all position for one rotations
plot_rotation = round(len(scans)/2)


#TODO 
#plt.imshow(np.log10(sum(((diff_data[0:][plot_rotation])))), interpolation='none')

plt.figure()
plt.imshow(np.log10(sum(((diff_data[0:][plot_rotation])))), interpolation='none')
#plt.title('Summed intensity for all rotations for scan %d (log)'%scans_sorted_theta[plot_rotation][2])
#
##%%
## ---------------------------------------------------------
## Do bright field analysis
##-----------------------------------------------------------
#
#def bright_field(data,x,y):
#    index = 0
#    photons = np.zeros((y,x)) 
#    for row in range(0,y):
#        for col in range(0,x):
#            photons[row,col] = sum(sum(data[index])) #/ max_intensity
#            index += 1            
#    return photons
#
## do BF for all rotations
#brightfield = np.zeros((len(scans), nbr_rows, nbr_cols))
#for jj in range(0,len(scans)):
#    brightfield[jj] = bright_field(diff_data[:,jj,:,:],nbr_cols,nbr_rows)
#    #Normalize each image ( so that i can plot a better 3D image)
#    #brightfield[jj] = brightfield[jj] / brightfield[jj].max()
#
#
#def plot_bright_field():
#    interval=100 #plotting interval
#    #plot every something 2d bright fields
#    for ii in range(0,len(scans),interval):
#        plt.figure()
#        plt.imshow(np.log10(brightfield[ii]), cmap='jet', interpolation='none', extent=extent_motorpos, vmin=0, vmax=5 ) 
#        plt.title('Bright field sorted in gonphi %d'%scans_sorted_theta[ii][2])  
#        plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
#        plt.ylabel('$y$ [$\mathrm{\mu m}$]') ; plt.colorbar()
#        ##plt.savefig('C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\ptycho_192_222\BF_Merlin\\scan%d'%((scans[ii])), bbox_inches='tight')
#        #plt.savefig('C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\scan461_\BF\BF_merlin_aligned3\\log10_scan%d'%((scans_sorted_theta[ii][2])), bbox_inches='tight')
#        #plt.savefig("BF/scan%d"%scans_gonphi[ii])   
#    # plot average bright field image (average over rotation)
#    plt.figure()
#    plt.imshow(np.mean(brightfield, axis=0), cmap='jet', interpolation='none')#,extent=extent_motorpos)
#    plt.title('Average image from bright field') 
#    plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
#    plt.ylabel('$y$ [$\mathrm{\mu m}$]') 
#    #plt.savefig("C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\scan461_\BF\BF_merlin_aligned3\\Average_brightfield", bbox_inches='tight')     
#    
#    plt.figure()
#    plt.imshow(sum(brightfield), cmap='jet', interpolation='none',extent=extent_motorpos)
#    plt.title('Bright field summed over all positions') 
#    plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
#    plt.ylabel('$y$ [$\mathrm{\mu m}$]') 
##plot_bright_field()
#
#
##%%
## COM analysis
##-----------------------
#
#def COM2d(data,nbr_cols,nbr_rows):
#    
#    index = 0
#    # define a vector with length of the length of roi on the detector
#    roix = np.linspace(1, data.shape[1], data.shape[1])
#    ## define a vector with length of the height of roi on the detector
#    roiy = np.linspace(1,data.shape[2],data.shape[2])     #shape 1 or 2 ?
#    # meshgrids for center of mass calculations
#    X, Y = np.meshgrid(roix,roiy)
#    
#    COM_hor = np.zeros((nbr_rows,nbr_cols))
#    COM_ver = np.zeros((nbr_rows,nbr_cols))
#    COM_mag = np.zeros((nbr_rows,nbr_cols))
#    COM_ang = np.zeros((nbr_rows,nbr_cols))
#    
#    for row in range(0,nbr_rows):
#        for col in range(0,nbr_cols):
#            COM_hor[row,col] = sum(sum(data[index]*X))/sum(sum(data[index]))
#            COM_ver[row,col] = sum(sum(data[index]*Y))/sum(sum(data[index]))
#            if row == 0 and col == 0:
#                bkg_hor = 0#152.4#COM_hor[row,col] 
#                bkg_ver = 0#101.8#COM_ver[row,col] 
#                        # DPC in polar coordinates. r then phi:
#            COM_mag[row, col] = np.sqrt( (COM_hor[row,col]-bkg_hor)**2 + (COM_ver[row,col]-bkg_ver)**2) 
#     #       COM_ver(row_ROI(1):row_idx, col_ROI(1):col_ROI(end),scan_idx)-mean(mean(COM_ver(row_ROI(1):row_idx,col_ROI(1):col_ROI(end),scan_idx))))
#    #[COM_angle(row_ROI(1):row_idx, col_ROI(1):col_ROI(end),scan_idx), COM_magnitude(row_ROI(1):row_idx,col_ROI(1):col_ROI(end),scan_idx)] = cart2pol(COM_hor(row_ROI(1):row_idx,nbr:cols,scan_idx)-mean(mean(COM_hor(row_ROI(1):row_idx,col_ROI(1):col_ROI(end),scan_idx))), grej        
#            COM_ang[row, col] = np.arctan( COM_hor[row,col] / COM_ver[row,col])
#       #     COM_hor_tmp(col) = sum(sum(image_dat(roiy,roix).*X))/sum(sum(image_dat(roiy,roix)));
#       #     COM_ver_tmp(col) = sum(sum(image_dat(roiy,roix).*Y))/sum(sum(image_dat(roiy,roix)));
#            index += 1 
#    return COM_mag
#
## do BF for all rotations
#COM_2d = np.zeros((len(scans), nbr_rows, nbr_cols))
#for jj in range(0,len(scans)):
#    COM_2d[jj] = COM2d(diff_data[:,jj,:,:],nbr_cols,nbr_rows)
#    #Normalize each image ( so that i can plot a better 3D image)
#    #brightfield[jj] = brightfield[jj] / brightfield[jj].max()
#
#def plot_COM2d():
#    interval=1 #plotting interval
#    #plot every something 2d bright fields
#    for ii in range(0,len(scans),interval):
#        plt.figure()
#        plt.imshow(COM_2d[ii], cmap='jet', interpolation='none', extent=extent_motorpos) 
#        plt.title('COM 2d mag %d'%scans_sorted_theta[ii][2])  
#        plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
#        plt.ylabel('$y$ [$\mathrm{\mu m}$]') 
#        #plt.savefig('C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\scan461_\COM\COM_merlin_aligned\InP\\scan%d'%((scans_sorted_theta[ii][2])), bbox_inches='tight')
#        #plt.savefig("BF/scan%d"%scans_gonphi[ii])   
#    # plot average bright field image (average over rotation)
##plot_COM2d()
#
##%%
#
#"""
## Try to deconvolve the probe size from the BF signal (make a vertical line 
## scan and deconvolve with square function of size 170 nm which should be the real diamter )
## Have to make an interpolation of the signal for this to work
## ----------------------------------------------------------------------------
#col_line = sum(brightfield)[:,33]
#xx_line = np.linspace(0, extent_motorpos[3],11)
#
#def box_function(x,limit1,limit2,low,high): 
#    func = high*np.ones((len(x)))
#    func[np.where(x <= limit1)] = low
#    func[np.where(x > limit2)] = low 
#    return func
#
## wire diameter corresponds to this many pixels:
#NW_diam_pix = 0.170 /dy    
#
#limit1 = xx_line[int(len(xx_line)/2 - (np.round(NW_diam_pix) /2) )]
#limit2 = xx_line[int(len(xx_line)/2 + (np.round(NW_diam_pix)/2) )]
#boxfun = box_function(xx_line,limit1,limit2,0,col_line.max() )
#
#plt.figure()
#plt.plot(xx_line,boxfun,'b-')
#plt.plot(xx_line,col_line,'m.')
#plt.xlabel('um')
#
#
#
##from scipy.signal import deconvolve
##probe_signal = deconvolve(col_line,step)
#
#
#"""
#
#
##%%
## help function. Makes fast-check-plotting easier. 
#def imshow(data):
#    plt.figure()
#    plt.imshow(data, cmap='jet')
#    plt.colorbar()
##%%    
#    
##  corrected 20190822  
## define q1 q2 q3 + q_abs from the geometry function 
## (See "Bending and tilts in NW..." pp)
#def def_q_vectors():
#    global q3, q1, q2, q_abs    
#    #  units of reciprocal meters [m-1]
#    q_abs = 4 * np.pi / g.lam * g.sintheta
#    # (x, z, y),  (r3, r1, r2), (qx, qz, qy), or (q3, q1, q2),
#    q1 = np.linspace(-g.dq1*diff_data.shape[2]/2.+q_abs/g.costheta, g.dq1*diff_data.shape[2]/2.+q_abs/g.costheta, diff_data.shape[2]) #        ~z
#    # q3 defined as centered around 0, that means adding the component from q1
#    q3 = np.linspace(-g.dq3*diff_data.shape[1]/2. + g.sintheta*q1.min() , g.dq3*diff_data.shape[1]/2.+ g.sintheta*q1.max(), diff_data.shape[1]) #~x
#    q2 = np.linspace(-g.dq2*diff_data.shape[3]/2., g.dq2*diff_data.shape[3]/2., diff_data.shape[3]) #         ~y
#def_q_vectors()
##%%
#
## --------------------------------------------------------------
## Make a meshgrid of q3 q1 q2 and transform it to qx qz qy.
## Also define the vectors qx qz qy
##----------------------------------------------------------------
#
##TODO 
#"""
## try to load the grid for the q-vectors
## detta är real space grid hos detektorn? 0.02 osv. Kan jag transformera den till reciprok?
#aa=P.diff.storages['S0000'].grids()
#q3_te, q1_te, q2_te = g.transformed_grid(aa, input_space='real', input_system='natural')
#vilket storage har grid i reciproca rummet?
#
#Är detta rätt?    att göra ett nytt grid och allt är bra 2pi/alla värden?
#2*np.pi/aa[0][9]
#
#
#"""
#
## in the transformation is should be input and output: (qx, qz, qy), or (q3, q1, q2).
## make q-vectors into a tuple to transform to the orthogonal system; Large Q means meshgrid, small means vector
##Q3,Q2,Q1 = np.meshgrid(q3, q2, q1, indexing='ij') 
#
#Q3,Q1,Q2 = np.meshgrid(q3, q1, q2, indexing='ij') 
#
#plt.figure()
#plt.subplot(311)
#plt.imshow(Q1[1])
#plt.title('Q1', loc='left', pad =-15)
#plt.colorbar()
#
#plt.subplot(312)
#plt.title('Q2', loc='left', pad =-15)
#plt.imshow(Q2[1])
#plt.colorbar()
#
#plt.subplot(313)
#plt.title('Q3', loc='left', pad =-15)
#plt.imshow(Q3[1])
#plt.colorbar()
#
#plt.figure()
#plt.title('Q3[:,:,0]', loc='left', pad =-15)
#plt.imshow(Q3[:,:,0])
#plt.colorbar()
#
#
#tup = Q3, Q1, -Q2   
#Qx, Qz, Qy = g.transformed_grid(tup, input_space='reciprocal', input_system='natural')
#
#
#plt.figure()
#plt.subplot(311)
#plt.imshow(Qz[1])
#plt.title('Qz', loc='left', pad =-15)
#plt.colorbar()
#
#plt.subplot(312)
#plt.title('Qy', loc='left', pad =-15)
#plt.imshow(Qy[1])
#plt.colorbar()
#
#plt.subplot(313)
#plt.title('Qx', loc='left', pad =-15)
#plt.imshow(Qx[1])
#plt.colorbar()
#
##NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOTTTTTTTTTTTTTTTTTTTTTTTTTCORECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTTTTT
#qx = np.linspace(Qx.min(),Qx.max(),g.shape[0])
#qz = np.linspace(Qz.max(),Qz.min(),g.shape[1])
#qy = np.linspace(Qy.max(),Qy.min(),g.shape[2])
##qy = np.linspace(Qy[0,0,0],Qy[0,-1,0],g.shape[1])
#
#
##%%
##-----------------------------------------------------------------------------------
##TODO make a copy of the data instead of putting into diff storage. but it must be a ptypy storage as input in coordinate_shift.
## here I put the masked data in the data. for all frames. So this I need to keep before xrd analysis. (or make a new storage with masked data) 
## TODO I already did this before, maksed the data, riight? no not in the storage
#P.diff.storages['S0000'].data = P.diff.storages['S0000'].data * P.mask.storages['S0000'].data
##%%
#
##-----------------------------------------------------------------------------    
## Find max scattering position and plot that positions 3d bragg peak in 2d cuts 
## using Berenguer terminology (It is not corrected and also I change the qx axis later!)
## Plot the 'naive' Bragg peak (not skewed coordinates) in a single position in 
## 3dim.
## test trying to skew the system (the diff data) from the measurement coordinate system q3q1q2 (in reciprocal space)
## to the orthogonal reciprocal space qxqzqy with the help of the ptypy class coordinate_shift in geometry_class.py 
## (qxqzqy is what ptypy calls it, not Berenguer which has no name for it)
#
#
#"# COrrected 20190822"
#
##------------------------------------------------------------------------------
#
## check which positions has the most intensity, for a nice 3d Bragg peak plot
#pos_vect_naive = np.sum(np.sum(np.sum(diff_data, axis =1), axis =1), axis =1)
#max_pos_naive = np.argmax(pos_vect_naive)
##max_pos_naive = 524# 512for 80 segment#524 for 170 segment
#
## find where the peak is in the detector plane
#q3max = np.argmax(np.sum(np.sum(diff_data[max_pos_naive],axis=1),axis=1))
#q2max = np.argmax(np.sum(sum(diff_data[max_pos_naive]),axis=1))
#q1max = np.argmax(np.sum(sum(diff_data[max_pos_naive]),axis=0))
#    
#
## save one bragg peak to plot. set 0 values to inf   
#plot_3d_naive = np.copy(diff_data[max_pos_naive])
#plot_3d_naive[plot_3d_naive==0] = np.inf
#
## plot the 'naive' 3d peak of the most diffracting position, and centered on the peak for q1 and q2 
#factor = 1E-10 
#plt.figure()
#plt.suptitle('Naive plot of single position Bragg peak in natural coord system')
#plt.subplot(221)
## extent (left, right, bottom, top) in data coordinates
#plt.imshow(np.rot90(plot_3d_naive[q3max],k=3 ), cmap='jet', interpolation='none', extent=[q1[-1]*factor, q1[0]*factor, q2[0]*factor, q2[-1]*factor])
#plt.ylabel('$q_2$ $ (\AA ^{-1}$)') ; plt.colorbar()                                       
#plt.xlabel('$q_1$ $ (\AA ^{-1}$)')   
#
#plt.subplot(222)
#plt.imshow(plot_3d_naive[:,q1max,:], cmap='jet', interpolation='none', extent=[q2[0]*factor, q2[-1]*factor, q3[-1]*factor, q3[0]*factor])
#plt.ylabel('$q_3$ $ (\AA ^{-1}$)'); plt.colorbar() 
#plt.xlabel('$q_2$ $ (\AA ^{-1}$)') 
#
#plt.subplot(223)
#plt.imshow(plot_3d_naive[:,:,q2max], cmap='jet', interpolation='none', extent=[q1[0]*factor, q1[-1]*factor, q3[-1]*factor, q3[0]*factor])
#plt.ylabel('$q_3$ $ (\AA ^{-1}$)'); plt.colorbar()
#plt.xlabel('$q_1$ $ (\AA ^{-1}$)') 
#
#
#
#
## test the shifting and plot
#test_shift_coord = ptypy.core.geometry_bragg.Geo_Bragg.coordinate_shift(g, P.diff.storages['S0000'], input_space='reciprocal',input_system='natural', keep_dims=True,layer=max_pos_naive)
#
## make 0 values white instead of colored
## TODO : make a scatter plot instead tp check that all axis are correct, make it more visible. 
##but remember that the data is fluipped in qx(rot)axis here
#
## change font in plot
##from matplotlib import rc
##font = {'size'   : 16}
##plt.rc('font', **font)
##
#
#
#
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
##rc('text', usetex=True)
#
## change font
##plt.rcParams['font.sans-serif'] = "sans-serif"
##plt.rcParams['font.family'] = "sans-serif"  
#
#
## functions to plot and save a 3d bragg peak(cuts)
## send in the bragg peak
## the + and - on the axes are not correct anymore
#def plot_bragg_peak(container,frame):
#    
#    #note that this input container will be changed ()
#
#    #test_shift_coord.data[0] = 30*test_shift_coord.data[0]/test_shift_coord.data[0].max()
#    # make 0 values white instead of colored
#    # TODO : make a scatter plot instead tp check that all axis are correct, make it more visible. 
#    #but remember that the data is flipped in qx(rot)axis here
#    
##    q2max = np.argmax(np.sum(sum(container.data[0]),axis=0))
##    q1max = np.argmax(np.sum(sum(container.data[0]),axis=1))
#    q3max = np.argmax(np.sum(np.sum(container.data[0],axis=1),axis=1))
#    container.data[0][container.data[0]<0.05] = np.inf    
#    #test_shift_coord.data[0][test_shift_coord.data[0]<6E7] = np.inf 
#    factor = 1E-10  #if you want to plot in reciprocal m or Angstroms, user 1 or 1E-10
#    plt.figure()
#    #plt.title('Single position Bragg peak in orthogonal system \n (Berenguer terminology) qxqzqy frame:%d'%frame)    
#    #plt.subplot(221)
#    #plt.imshow((container.data[0][q3max]), cmap='jet', interpolation='none', extent=[ qy[-1]*factor, qy[0]*factor, qz[-1]*factor, qz[0]*factor])
#    # extent left, right, bottom, top
#    plt.imshow((np.rot90(container.data[0][q3max],k=3 )), cmap='jet', interpolation='none', extent=[ qz[-1]*factor, qz[0]*factor, qy[-1]*factor, qy[0]*factor])
#    plt.ylabel('$q_y$ $ (\AA ^{-1}$)')
#    plt.xlabel('$q_z$ $ (\AA ^{-1}$)')    
#    plt.tight_layout()
#    plt.colorbar(fraction=0.046, pad=0.04)
#    plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/InP_diffraction/nice_ones/170/%s_bragg_slices_pos_%d_1'%(date_str,frame)) 
#    plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/InP_diffraction/nice_ones/170/%s_bragg_slices_pos_%d_1.pdf'%(date_str,frame))
#    #plt.subplot(222)
#    plt.figure()
#    # should be x angainst y in these labels
#    plt.imshow((container.data[0][:,q1max,:]), cmap='jet', interpolation='none', extent=[qy[-1]*factor, qy[0]*factor, qx[0]*factor, qx[-1]*factor])
#    plt.ylabel('$q_x$ $ (\AA ^{-1}$)'); plt.colorbar(fraction=0.046, pad=0.04)
#    plt.xlabel('$q_y$ $ (\AA ^{-1}$)')     
#    plt.tight_layout()
#    #plt.subplot(223)
#    plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/InP_diffraction/nice_ones/170/%s_bragg_slices_pos_%d_2'%(date_str,frame)) 
#    plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/InP_diffraction/nice_ones/170/%s_bragg_slices_pos_%d_2.pdf'%(date_str,frame))
#    plt.figure()
#    plt.imshow((container.data[0][:,:,q2max]), cmap='jet', interpolation='none', extent=[qz[-1]*factor, qz[0]*factor, qx[0]*factor, qx[-1]*factor])
#    #plt.imshow((container.data[0][:,:,q2max]), cmap='jet', interpolation='none', extent=[qz[0]*factor, qz[-1]*factor, qx[0]*factor, qx[-1]*factor]) (qz (shold be like this I think))
#    plt.ylabel('$q_x$ $ (\AA ^{-1}$)'); plt.colorbar(fraction=0.046, pad=0.04)
#    plt.xlabel('$q_z$ $ (\AA ^{-1}$)') 
#    plt.tight_layout()
#    
#    container.data[0][container.data[0]==np.inf] = 0
#    
#    plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/InP_diffraction/nice_ones/170/%s_bragg_slices_pos_%d_3'%(date_str,frame)) 
#    plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/InP_diffraction/nice_ones/170/%s_bragg_slices_pos_%d_3.pdf'%(date_str,frame)) 
#
##plot_bragg_peak(test_shift_coord,max_pos_naive)
#
## x,y,z    # så jag har det
##3,2,1
##%%
#
################################################################################
## XRD analysis
################################################################################
#
## TODO: remove single photon count, if the COM is calculated for very small values like pixels with 1 photon counts, 
##then the result will be missleading. Set a threshold that keeps the resulting pixel on a mean value, like if sum(sum(sum(diffPattern)))< threshold. sum(sum(sum()))==bkg_value
## input is 4d matrix with [nbr_diffpatterns][nbr_rotations][nbr_pixels_x][nbr_pixels_y]
#def COM_voxels_reciproc(data, vect_Qx, vect_Qz, vect_Qy ):
#
#    # meshgrids for center of mass calculations in reciprocal space
#    COM_qx = np.sum(data* vect_Qx)/np.sum(data)
#    COM_qz = np.sum(data* vect_Qz)/np.sum(data)
#    COM_qy = np.sum(data* vect_Qy)/np.sum(data)
#
#   # print 'coordinates in reciprocal space:'
#   # print COM_qx, COM_qz, COM_qy
#    return COM_qx, COM_qz, COM_qy
#
## loop through all scanning postitions and move the 3D Bragg peak from the 
## natural to the orthogonal coordinate system (to be able to calculate COM)
## Calculate COM for every peak - this gives the XRD matrices
#def XRD_analysis():
#    position_idx = 0
#    XRD_qx = np.zeros((nbr_rows,nbr_cols))
#    XRD_qz = np.zeros((nbr_rows,nbr_cols))
#    XRD_qy = np.zeros((nbr_rows,nbr_cols))
#
#    for row in range(0,nbr_rows):
#        for col in range(0,nbr_cols):
#            
#            # if keep_dims is False, shouldnt the axis qz change? (q1 -->qz)
#            
#            
#            
#            data_orth_coord = ptypy.core.geometry_bragg.Geo_Bragg.coordinate_shift(g, P.diff.storages['S0000'], input_space='reciprocal',
#                         input_system='natural', keep_dims=True,
#                         layer=position_idx)         # layer is the first col in P.diff.storages.values()[0]
#      
#            
#            # TEST do analysis on the unshifted data. unshifted 
#            #COM_qx, COM_qz, COM_qy = COM_voxels_reciproc(P.diff.storages['S0000'].data[position_idx,::-1,:,:], Q3, Q1, Q2)
#            #data_orth_coord = 1
#            
#            # do the 3d COM analysis to find the orthogonal reciprocal space coordinates of each Bragg peak
#            # IMPORTATNE rotations are sorted with gonphi but for this to be correct with the coordinate system the 
#            # data should sorted with higher index = higher theta that why I fli the order here. The data in diff data are not sorted correct in qx
#
#            COM_qx, COM_qz, COM_qy = COM_voxels_reciproc(data_orth_coord.data[0][::-1,:,:], Qx, Qz, Qy)
#            
#            # insert coordinate in reciprocal space maps 
#            XRD_qx[row,col] = COM_qx
#            XRD_qz[row,col] = COM_qz
#            XRD_qy[row,col] = COM_qy
#            
#            position_idx +=1
#        
#            # plot bragg peak        
#            #plot_bragg_peak(data_orth_coord,position_idx)
##       
#
#     #plot every other 3d peak and print out the postion of the COM analysis
##            if (position_idx%500==0):
##                # TODO very har to say anything about this looking in 2d, need 3d plots!
##               
##                x_p = np.argwhere(qx>COM_qx)[0][0]
##                y_p = np.argwhere(qy>COM_qy)[0][0] #take the first value in qy where
##                z_p = np.argwhere(qz>COM_qz)[0][0]  
###                import pdb; pdb.set_trace()
##                plt.figure()
##                plt.title('COM_qz = %d, COM_qy = %d \n positions indx:%d '%((COM_qz,COM_qy,position_idx-1)))
##                plt.imshow(sum(data_orth_coord.data[0]), cmap='jet')#, extent=[ qy[0], qy[-1], qz[0], qz[-1] ])
##                # Find the coordinates of that cell closest to this value:              
##                plt.scatter(y_p, z_p, s=500, c='red', marker='x')#, extent=[ qy[0], qy[-1], qz[0], qz[-1] ])
##                #plt.title('Single Bragg peak summed in x. COM z and y found approx at red X')
###                import pdb
###                pdb.set_trace()    
#
#
#    return XRD_qx, XRD_qz, XRD_qy, data_orth_coord
#
#XRD_qx, XRD_qz, XRD_qy, data_orth_coord = XRD_analysis() # units of 1/m
#
##test plot for the coordinate system: (only works for the last position, the other peaks are not saved)
#def test_coordShift():            
#    plt.figure()
#    plt.imshow(np.log10(np.sum(diff_data[-1],axis=2)), cmap='jet', interpolation='none', extent=[ q1[0], q1[-1], q3[0], q3[-1] ])
#    plt.title('natural')
#    plt.xlabel('$q_1$ $ (\AA ^{-1}$)')   #l(' [$\mu m$]')#
#    plt.ylabel('$q_3$ $ (\AA ^{-1}$)')
#    plt.colorbar()
#    plt.figure()
#    plt.imshow(np.log10(np.sum(data_orth_coord.data[0],axis=2)), cmap='jet', interpolation='none')
#    plt.title('cartesian')
#    plt.xlabel('$q_x$ $ (\AA ^{-1}$)')
#    plt.ylabel('$q_z$ $ (\AA ^{-1}$)')
#    plt.colorbar()        
##test_coordShift()
#factor = 1E-10  #if you want to plot in m or Angstroms, user 1 or 1E-10
#
#def plot_XRD_xyz():
#    # plot reciprocal space map x y z 
#    plt.figure()
#    plt.subplot(411)
#    plt.imshow(factor*XRD_qx, cmap='jet',interpolation='none',extent=extent_motorpos)
#    plt.title('Reciprocal space map, $q_x$ $ (\AA ^{-1}$) ')
#    plt.ylabel('y [$\mu m$]')
#    plt.colorbar()
#    plt.subplot(412)
#    plt.imshow(factor*XRD_qy, cmap='jet',interpolation='none',extent=extent_motorpos) 
#    plt.title('Reciprocal space map, $q_y$ $ (\AA ^{-1}$) ')
#    plt.ylabel('y [$\mu m$]')
#    plt.colorbar()
#    plt.subplot(413)
#    plt.imshow(factor*XRD_qz, cmap='jet',interpolation='none',extent=extent_motorpos)
#    plt.title('Reciprocal space map, $q_z$ $(\AA ^{-1}$) ')
#    plt.ylabel('y [$\mu m$]')
#    plt.colorbar()
#    plt.subplot(414)
#    plt.title('Bright field (sum of all rotations)')
#    plt.imshow(sum(brightfield)/sum(brightfield).max(), cmap='jet', interpolation='none',extent=extent_motorpos)
#    plt.xlabel('x [$\mu m$]') 
#    plt.ylabel('y [$\mu m$]')
#    plt.colorbar()
##plot_XRD_xyz()
#
##%% 
##----------------------------------------------------------
## Convert q-vector from  cartesian coordinates to spherical
## (See "Bending and tilts in NW..." pp)
##----------------------------------------------------------
#
#XRD_absq =  np.sqrt(XRD_qx**2 + XRD_qy**2 + XRD_qz**2)
#XRD_alpha = np.arcsin( XRD_qy/ XRD_absq)
#XRD_beta  = np.arctan( XRD_qx / XRD_qz)
#
#
## for np arrays in np format
##np.save('XRD_absq', XRD_absq)
##np.save('filename', XRD_absq)
##np.save('filename', XRD_absq)
#
##open this file
##loaded_array = np.load('XRD_absq') #C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\ptycho_192_222\COM2d\scan222_COM_ang.npy')
#
#
##%%
##---------------------------------
## plot the XRD maps
##----------------------------------
#class Formatter(object):
#    def __init__(self, im):
#        self.im = im
#    def __call__(self, x, y):
#        z = self.im.get_array()[int(y), int(x)]
#        return 'x=%i, y=%i, z=%1.4f' % (x, y, z)
#
#
##def plot_XRD_polar():    
## cut the images in x-range:start from the first pixel: 
## remove the 1-pixel. (cheating) Did not remove it the extent, because then the 0 will not be shown in the X-scale and that will look weird
#start_cutXat = 0 
## whant to cut to the right so that the scale ends with an even number
##x-pixel nbr 67 is at 2.0194197798363955
#cutXat = nbr_cols+1 # 67
## replace the x-scales end-postion in extent_motorposition. 
#extent_motorpos_cut = np.copy(extent_motorpos)
####extent_motorpos_cut[1] = 2.0194197798363955 segmentedNW
#cutXat = 49    #    68 (for ending at 2)
#
#extent_motorpos_cut = [start_cutXat, dx*(cutXat-1),0,dy*(nbr_rows-1)]
#print('mean value BF: ')
#print( np.mean(sum(brightfield)))
## create a mask from the BF matrix, for the RSM analysis
#XRD_mask = np.copy(sum(brightfield))
#XRD_mask[XRD_mask < 24000 ] = np.nan #inP 25000  #ingap 16000   # for homo InP, use 280000
#XRD_mask[XRD_mask > 0] = 1       #make binary, all values not none to 1
#
#colormap = 'RdBu_r' #Spectral' #RdYlGn'#coolwarm' # 'bwr' #'PiYG'# #'RdYlBu' # 'PiYG'   # 'PRGn' 'BrBG' 'PuOr'
#
## plot abs q to select pixels that are 'background', not on wire, and set these pixels to NaN (make them white)
#fig=plt.figure()
##plt.suptitle(colormap)
##plt.suptitle(
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.02)
#plt.subplot(411)
#
#plt.title('Summed up intensity', loc='left', pad =-12, color ='white')
#plt.imshow(sum(brightfield[:,:,start_cutXat:cutXat])/sum(brightfield[:,:,start_cutXat:cutXat]).max(), cmap=colormap, interpolation='none',extent=extent_motorpos_cut)
#plt.ylabel('$y$ [$\mathrm{\mu}$m]')
#plt.xticks([])
##po = plt.colorbar(ticks=(10,20,30,40))#,fraction=0.046, pad=0.04) 
#po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
#
#
## if you want no mask use:
##XRD_mask = np.ones((XRD_mask.shape))
#
#
#
#plt.subplot(412)   
##calculate lattice constant a from |q|:       
## TODO: this is wrong for the homogenous wires, its not (111), for segmented InP i dont know    
#d =  2*np.pi/  XRD_absq
#a_lattice_exp = d * np.sqrt(3)                
##a_lattice_exp = np.pi*2./ (XRD_absq *np.sqrt(3))
##print 'mean lattice constant is %d' %np.mean(a_lattice_exp)
##imshow(a_lattice_exp)
##plt.title('Lattice conastant a [$\AA$]')
#mean_a = np.nanmean(XRD_mask[:,start_cutXat:cutXat]*a_lattice_exp[:,start_cutXat:cutXat])
##TODO try with reference strain equal to the center of the largest segment (for InP) # tody try with reference from the other NWs
##mean_strain = a_lattice_exp[:,start_cutXat:cutXat].max() 
#
#strain = 100*XRD_mask[:,start_cutXat:cutXat]*(a_lattice_exp[:,start_cutXat:cutXat]-mean_a)/mean_a
#
#
#plt.imshow(strain,interpolation='none',extent=extent_motorpos_cut, cmap=colormap)
##plt.imshow(XRD_mask[:,start_cutXat:cutXat]*a_lattice_exp[:,start_cutXat:cutXat], cmap='jet',interpolation='none')#,extent=extent_motorpos_cut) 
#
##plt.imshow(XRD_mask[:,start_cutXat:cutXat]*XRD_absq[:,start_cutXat:cutXat]*1E-8, cmap=colormap,interpolation='none',extent=extent_motorpos_cut) 
##plt.title('Length of Q-vector |Q|)', loc='left', pad =-12)
#plt.title('Axial Strain $\epsilon$ (%)', loc='left', pad =-12)   #
##plt.title('Lattice constant a')
#plt.ylabel('$y$ [$\mathrm{\mu}$m]')  
#plt.xticks([])
#po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
#
#plt.subplot(413)
#plt.imshow(XRD_mask[:,start_cutXat:cutXat]*1E3*XRD_alpha[:,start_cutXat:cutXat], cmap=colormap, interpolation='none',extent=extent_motorpos_cut) # not correct!
## cut in extent_motorposition. x-pixel nbr 67 is at 2.0194197798363955
#plt.title('$\\alpha$ (mrad)', loc='left', pad =-12)
#plt.ylabel('$y$ [$\mathrm{\mu}$m]')
#plt.xticks([])
#po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
#
##po = plt.colorbar(ticks=(0,1,2,3,4))
##po.set_label('Bending around $q_x$ $\degree$')
#   
#plt.subplot(414)
#plt.imshow(XRD_mask[:,start_cutXat:cutXat]*1E3*XRD_beta[:,start_cutXat:cutXat], cmap=colormap,interpolation='none',extent=extent_motorpos_cut) # not correct!
#plt.title('$\\beta$ (mrad)', loc='left', pad =-12)
#plt.ylabel('$y$ [$\mathrm{\mu}$m]')
#plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
#po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
#
##plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InGaP/%s_xrd'%(date_str))
##plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InGaP/%s_xrd.pdf'%(date_str)) 
##plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InP/%s_xrd'%(date_str))
##plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InP/%s_xrd.pdf'%(date_str)) 
#
#
## calculate mean value in each segment (InP)
#mean45 = np.nanmean(strain[:,0:10])
#mean80 = np.nanmean(strain[:,10:25])
#mean170 = np.nanmean(strain[:,26:-1])
#print( mean45,mean80,mean170)
#
#def q_abs_map():
#    q_map_temp = XRD_mask[:,start_cutXat:cutXat]*XRD_absq[:,start_cutXat:cutXat]
#    return q_map_temp
#q_map_temp = q_abs_map()    
#
##mean value of the Q_abs map
#mean_absq = np.nanmean(q_map_temp)
#print( 'Mean |q| value in the masked map' )
#print( mean_absq)
#
#plt.figure()
#plt.title('Mean |q| is: %f '%mean_absq )
#plt.imshow(q_map_temp); plt.colorbar()
##plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InP/%s_xrd_with_mean_qabs_value.pdf'%(date_str)) 
#plt.figure()
#plt.title('lattice constant a')
#plt.imshow(a_lattice_exp[:,start_cutXat:cutXat]); plt.colorbar()
##%%
#
#from scipy.optimize import curve_fit
#def line_fun(x,k,m):
#    return k*x + m
#
#
##a_InGaP_ref = 5.8060E-10
## plot 4 values compared to a table reference 
#def Ga_gradient_plot_points(a_map,indxs1,indxs2,indxs3,indxs4):
#    x = np.linspace(0,extent_motorpos_cut[1],4)
# 
#    a_inGaP_values = [] # acctually strain
#    #indxs1 = (1,1)
#    a_inGaP_values.append(1e10*a_map[indxs1]*(1.0 - 0.019))
#    a_inGaP_values.append(1e10*a_map[indxs2]*(1.0 - 0.019))
#    a_inGaP_values.append(1e10*a_map[indxs3]*(1.0 - 0.019))
#    a_inGaP_values.append(1e10*a_map[indxs4]*(1.0 - 0.019))
#    # make linear fit
#    popt, pcov = curve_fit( line_fun, x, a_inGaP_values)
#    # fitted function is
#    fit_gainp = line_fun(x,*popt)
#    
#    # calculate %Ga backwards from Vigards law but with 1.9% error
#    ga_proc = []    
#    ga_proc.append(100.0*st.Vigards_law_back(a_map[indxs1]*(1.0 - 0.019)*1e10 ))
#    ga_proc.append(100.0*st.Vigards_law_back(a_map[indxs2]*(1.0 - 0.019)*1e10 ))
#    ga_proc.append(100.0*st.Vigards_law_back(a_map[indxs3]*(1.0 - 0.019)*1e10 ))
#    ga_proc.append(100.0*st.Vigards_law_back(a_map[indxs4]*(1.0 - 0.019)*1e10 ))
#    # make linear fit
#    popt, pcov = curve_fit( line_fun, x, ga_proc)
#    # fitted function is
#    fit_proc = line_fun(x,*popt)
#            
#    fig,ax1 = plt.subplots()
#    color = 'tab:blue'
#    ax1.plot(x,a_inGaP_values,marker='o',color = color,linestyle = 'None')
#    ax1.plot(x,fit_gainp,color=color)
#    #plt.title('$\\beta$ (mrad)', loc='left', pad =-12)
#    ax1.tick_params(axis='y', labelcolor=color)
#    ax1.set_ylabel('Lattice constant $a_{GaInP}$ [$\AA$]',color=color)
#    ax1.set_xlabel('$x$ [$\mathrm{\mu m}$]') 
#    
#    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#    color2 ='g'
#    ax2.set_ylabel('Ga content  [$\%$]', color=color2)
#    ax2.tick_params(axis='y', labelcolor=color2)
#    ax2.plot(x,ga_proc,color = color2, marker = '^', linestyle = 'None')
#    ax2.plot(x,fit_proc, color = color2)
#    fig.tight_layout()  # otherwise the right y-label is slightly clipped
#Ga_gradient_plot_points(a_lattice_exp[:,start_cutXat:cutXat],indxs1=(6,2),indxs2=(6,12),indxs3=(6,24),indxs4=(6,42))
#
##%%
#
## make a plot of Ga content and lattice constant as a function of x
#def Ga_gradient_plot(a_map,row):
#    
#    # calculate %Ga backwards from Vigards law but with 1.9% error
#    ga_proc = []
#    for x in a_map[row]:
#         ga_proc.append(100.0*st.Vigards_law_back(x*(1.0 - 0.019)*1e10 ))
#        
#    fig,ax1 = plt.subplots()
#    color = 'tab:blue'
#    # multiply lattice constant this value by 1.9% since its 1.9% off:
#    #*(1.0 - 0.019)
#    ax1.plot(np.linspace(0,extent_motorpos_cut[1],a_map.shape[1]),a_map[row]*1e10*(1.0 - 0.019))
#    #plt.title('$\\beta$ (mrad)', loc='left', pad =-12)
#    ax1.tick_params(axis='y', labelcolor=color)
#    ax1.set_ylabel('Lattice constant $a_{GaInP}$ [$\AA$]',color=color)
#    ax1.set_xlabel('$x$ [$\mathrm{\mu m}$]') 
#    
#    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#    ax2.set_ylabel('Ga content  [$\%$]', color='g')
#    ax2.tick_params(axis='y', labelcolor='g')
#    ax2.plot(np.linspace(0,extent_motorpos_cut[1],a_map.shape[1]),ga_proc,'g')
#    fig.tight_layout()  # otherwise the right y-label is slightly clipped
#Ga_gradient_plot(a_lattice_exp[:,start_cutXat:cutXat],row=6)
#
##%%
#def XRD_lineplot():
#    plot = 10
#    start = 0
#    
#    x = np.linspace(0,extent_motorpos_cut[1],XRD_absq.shape[1])
#    plt.figure()
#    plt.suptitle('XRD lineplots')
#    plt.subplot(411)
#    for i in range(start,XRD_absq.shape[0]):
#        plt.plot(x,XRD_absq[i,:] ) # choose a row in the middle of the wire    
#    plt.legend()
#    plt.subplot(412)
#    for i in range(start,XRD_absq.shape[0]):
#        plt.plot(x,XRD_alpha[i,:] ) # choose a row in the middle of the wire  
#    plt.subplot(413)
#    for i in range(start,XRD_absq.shape[0]):
#        plt.plot(x,XRD_beta[i,:] ) # choose a row in the middle of the wire  
#    plt.subplot(414)
#    
#    plt.figure()
#    plt.title('BF intensity summed in y')
#    plt.plot(sum(brightfield)[plot,:] )
#    plt.xlabel('x [um]')
#
#XRD_lineplot()    
#    
##%%
#def save_np_array(nparray):
#    np.save('filename', nparray)
#    
#
#def rocking_curve_plot():
#    #458: theta= 13.1    #515 theta = 12.1
#    # Is the first point 458? they are sorted according to gonphi, right? In that case it is right.
#    # find highest intensity point
#    alla = np.sum(np.sum(np.sum(diff_data,axis=1),axis=1),axis=1)
#    index_max = np.argmax(alla)
#    
#    index_80 = 87*7+24
#    index_45 = 87*7+12
#    index_19 = 87*8+4
#    # intensity np.sum(diff_data[87*7+12])
#    
#    #which rocking curve do u want to plot:
#    index = index_19
#    
#    #theta = np.linspace(12.1,13.1,51    # dthetea =0.02   
#    plt.figure(); plt.plot((np.sum(np.sum(diff_data[index],axis=1),axis=1)))
#    plt.yscale('log'); plt.title('Rocking curve at 19nm segment')    
#    plt.ylabel('Photon counts');plt.xlabel('Rotation $\Theta$ ($\degree$)')
#    plt.grid(True)
#rocking_curve_plot()
# 
##    
#################################################################################
### Testing a few Mayavi plotting functions
#################################################################################    
###def plot_3d_isosurface():
##from mayavi import mlab   #if you can do this instde function it is ood because it changes to QT fram    
##
### check which positions has the most intensity, for a nice 3d Bragg peak plot
###TODO this is for the 'naive' system (diff_data is)
##pos_vect = np.sum(np.sum(np.sum(diff_data, axis =1), axis =1), axis =1)
##max_pos = np.argmax(pos_vect)
##
##plt.figure()
##plt.text(5, np.max(pos_vect), 'Max at: ' + str(max_pos), fontdict=None, withdash=False)
##plt.plot(pos_vect)
##plt.title('Summed intensity as a function of position')
###plt.savefig('SumI_vs_position')
##
##
### plot this position:
##pos = max_pos +0
### change the coordinate system of this data
##data_orth_coord = ptypy.core.geometry_bragg.Geo_Bragg.coordinate_shift(g, P.diff.storages.values()[0], input_space='reciprocal',
##                         input_system='natural', keep_dims=True,
##                         layer=pos) 
### check to see it is fine
##
##plt.figure()
##plt.imshow(data_orth_coord.data[0,0])
##plt.title('InP Bragg peak projection with fringes, orth coord')
##
##plot_data = diff_data[pos] #data_orth_coord.data[0]     #data[0]0
##
##
##def slice_plot():
##    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(plot_data),
##                                plane_orientation='x_axes',
##                                slice_index=10,
##                            )
##    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(plot_data),
##                                plane_orientation='y_axes',
##                                slice_index=10,
##                            )
##    mlab.outline()
##    
##slice_plot()
##
##def plot3dvolume(): #  this looks very good, but almost never works 
##    x, y, z = np.ogrid[-10:10:20j, -10:10:20j, -10:10:20j]
##    s = np.sin(x*y*z)/(x*y*z)
##    mlab.pipeline.volume(mlab.pipeline.scalar_field(s))
##    # pipeline.scalar_filed makes data on a regular grid
##    #mlab.pipeline.volume(data[max_pos], vmin=0, vmax=0.8)
###plot3dvolume()
##
### TODO axes must be wrong here. check
###def contour3d():   #iso surface
##mlab.figure()
##xmin=qx[0]*1E0; xmax = qx[-1]*1E0; ymin=qy[0]*1E0; ymax=qy[-1]*1E0; zmin=qz[0]*1E0; zmax=qz[-1]*1E0
##obj = mlab.contour3d( plot_data, contours=10, opacity=0.5, transparent=False, extent=[ qz[0], qz[-1],qy[0], qy[-1] , qx[0], qx[-1] ])  #  , vmin=0, vmax=0.8)
##mlab.axes(ranges=[xmin, xmax, ymin, ymax, zmin, zmax])
##mlab.xlabel('$Q_z$ [$\AA^{-1}$]'); mlab.ylabel('$Q_y$ [$\AA^{-1}$]'); mlab.zlabel('$Q_z$ [$\AA^{-1}$]')
###C:\Users\Sanna\Documents\Beamtime\NanoMAX062017\Analysis_ptypy\scan461_\bragg_peak_stacking\InP\
###mlab.savefig('pos_'+ str(pos) +'.jpg')
##
##
### Another way to maka an iso surface. Can also be combined with cut planes
##position = 0 
##def iso_pipeline_plot():
##    src = mlab.pipeline.scalar_field(np.log10(plot_data))  # this creates a regular space data
##    mlab.pipeline.iso_surface(src, contours=[diff_data[position].min()+0.1*diff_data[position].ptp(), ], opacity=0.5)
##    mlab.show()    
##iso_pipeline_plot()
##
## 
#################################################################################
### MOVIE makers
#################################################################################
##
### calls a movie maker function in another script
##movie_maker(abs((diff_data[:,13])),'all_poitions_S192') 
##
### alternative movie maker. 
##def movie_maker2(data, name, rotation_nbr, nbr_plots):
##    if nbr_plots == 1:
##        #figure for animation
##        fig = plt.figure()
##        # Initialize vector for animation data
##        ims = []  
##        index = 0
##        for ii in range(0,len(data)):
##
##                
##                im = plt.imshow(np.log10(data[index][rotation_nbr]), animated=True, cmap = 'jet', interpolation = 'none')#, origin='lower')
##    #        #plt.clim(0,4) to change range of colorbar
##    #        im = plt.title('Col %d'%i)    # does not work for movies
##    #        txt = plt.text(0.2,0.8,i)   # (x,y,string)
##                ims.append([im])    #ims.append([im, txt])
##        ani = animation.ArtistAnimation(fig, ims, interval=2000, blit=True,repeat_delay=0)  
##        #txt = plt.text(0.1,0.8,'row: ' + str(y) + ' col: ' + str(x) )  # (x,y,string)
##        ims.append([im])
##        #ims.append([[im],txt])
##        plt.axis('off')
##        plt.show()
##        # save animation:
##        ani.save(name +'.mp4', writer="mencoder") 
##    elif nbr_plots == 2:           
##        fig, (ax1, ax2)   = plt.subplots(2,1)  
##        ims = [] 
##        # for x and y index of 1 col
##        #index = 0
##        for index in range(5,5+21):
##            #for x in range(20,50):
##    #            im = plt.imshow(np.log10(data[index]), animated=True, cmap = 'jet', interpolation = 'none')#, origin='lower')
##                #plt..subplot(21)
##                im = ax1.imshow(np.log10(data[index][rotation_nbr]), animated=True, cmap = 'jet', interpolation = 'none')#, origin='lower')
##                #plt.subplot(22)
##                im2 = ax2.imshow(np.log10(data[index][rotation_nbr]), animated=True, cmap = 'jet', interpolation = 'none')#, origin='lower')
##                #plt.subplot(313)
##                #im3 = plt.imshow(np.log10(data[index])[2], animated=True, cmap = 'jet', interpolation = 'none')#, origin='lower')
##                index += 1
##                #plt.clim(0,4) to change range of colorbar
##                #im = plt.title('Angle %d'%i)    # does not work for movies
##     #           txt = plt.text(0.1,0.8,'row: ' + str(y) + ' col: ' + str(x) )  # (x,y,string)
##                #ims.append([[im, im2],txt])
##                ims.append([im, im2])           
##                
##        ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True,repeat_delay=0)  
##        plt.axis('on')
##        plt.show()
##        # save animation:
##        #ani.save(name +'.mp4', writer="mencoder")    
##movie_maker2(diff_data,'all_diff_data_S192',rotation_nbr=13, nbr_plots=1)
##
##plt.figure()
##plt.imshow(np.log10(diff_data[147][13]), cmap = 'jet', interpolation = 'none')#, origin='lower')
#
##%%
#
#
##def full_range_xrd_plot():
#
#
#
## replace the x-scales end-postion in extent_motorposition. 
#extent_motorpos_cut = [-extent_motorpos[1]/2, extent_motorpos[1]/2, -extent_motorpos[3]/2, extent_motorpos[3]/2]# np.copy(extent_motorpos)
## create a mask from the BF matrix, for the RSM analysis
#XRD_mask = np.copy(sum(brightfield))
#XRD_mask[XRD_mask < 24000 ] = np.nan #inP 25000  #ingap 16000   # for homo InP, use 280000
#XRD_mask[XRD_mask > 0] = 1       #make binary, all values not none to 1
#
#colormap = 'RdBu_r' #Spectral' #RdYlGn'#coolwarm' # 'bwr' #'PiYG'# #'RdYlBu' # 'PiYG'   # 'PRGn' 'BrBG' 'PuOr'
#
## plot abs q to select pixels that are 'background', not on wire, and set these pixels to NaN (make them white)
#fig=plt.figure()
##plt.suptitle(colormap)
##plt.suptitle(
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.02)
#plt.subplot(211)
#
#plt.title('Summed up intensity', loc='left', pad =-12, color ='black')
#plt.imshow(np.fliplr((sum(brightfield)/sum(brightfield).max())), cmap=colormap, interpolation='none',extent=extent_motorpos_cut)
#plt.ylabel('$y$ [$\mathrm{\mu}$m]')
##plt.xticks([])
#plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
##po = plt.colorbar(ticks=(10,20,30,40))#,fraction=0.046, pad=0.04) 
#po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
#
#
## if you want no mask use:
##XRD_mask = np.ones((XRD_mask.shape))
#
#
#
#plt.subplot(212)   
##calculate lattice constant a from |q|:       
## TODO: this is wrong for the homogenous wires, its not (111), for segmented InP i dont know    
#d =  2*np.pi/  XRD_absq
#a_lattice_exp = d * np.sqrt(3)                
##a_lattice_exp = np.pi*2./ (XRD_absq *np.sqrt(3))
##print 'mean lattice constant is %d' %np.mean(a_lattice_exp)
##imshow(a_lattice_exp)
##plt.title('Lattice conastant a [$\AA$]')
#mean_a = np.nanmean(XRD_mask*a_lattice_exp)
##TODO try with reference strain equal to the center of the largest segment (for InP) # tody try with reference from the other NWs
##mean_strain = a_lattice_exp[:,start_cutXat:cutXat].max() 
#
#strain = 100*XRD_mask*(a_lattice_exp-mean_a)/mean_a
#
#
#plt.imshow(np.fliplr(strain),interpolation='none',extent=extent_motorpos_cut, cmap=colormap)
##plt.imshow(XRD_mask[:,start_cutXat:cutXat]*a_lattice_exp[:,start_cutXat:cutXat], cmap='jet',interpolation='none')#,extent=extent_motorpos_cut) 
#
##plt.imshow(XRD_mask[:,start_cutXat:cutXat]*XRD_absq[:,start_cutXat:cutXat]*1E-8, cmap=colormap,interpolation='none',extent=extent_motorpos_cut) 
##plt.title('Length of Q-vector |Q|)', loc='left', pad =-12)
#plt.title('Axial Strain $\epsilon$ (%)', loc='left', pad =-12)   #
##plt.title('Lattice constant a')
#plt.ylabel('$y$ [$\mathrm{\mu}$m]')  
#plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
#po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
#
##plt.subplot(413)
##plt.imshow(XRD_mask*1E3*XRD_alpha, cmap=colormap, interpolation='none',extent=extent_motorpos_cut) # not correct!
### cut in extent_motorposition. x-pixel nbr 67 is at 2.0194197798363955
##plt.title('$\\alpha$ (mrad)', loc='left', pad =-12)
##plt.ylabel('$y$ [$\mathrm{\mu}$m]')
##plt.xticks([])
##po = plt.colorbar()
##tick_locator = ticker.MaxNLocator(nbins=4); po.locator = tick_locator;po.update_ticks()
##
###po = plt.colorbar(ticks=(0,1,2,3,4))
###po.set_label('Bending around $q_x$ $\degree$')
##   
##plt.subplot(414)
##plt.imshow(XRD_mask*1E3*XRD_beta, cmap=colormap,interpolation='none',extent=extent_motorpos_cut) # not correct!
##plt.title('$\\beta$ (mrad)', loc='left', pad =-12)
##plt.ylabel('$y$ [$\mathrm{\mu}$m]')
#plt.xlabel('$x$ [$\mathrm{\mu m}$]') 
##po = plt.colorbar()
#tick_locator = ticker.MaxNLocator(nbins=4); #po.locator = tick_locator;po.update_ticks()
#
##plt.savefig('C:/Users/Sanna/Documents/Beamtime/NanoMAX062017/Analysis_ptypy/scan461_/scanningXRD_InGaP/%s_xrd'%(date_str))
#
## for np arrays in np format
#np.save('BF_normalized', np.fliplr((sum(brightfield)/sum(brightfield).max())))
#np.save('extent_motorpos_cut', extent_motorpos_cut)
#np.save('strain', np.fliplr(strain))
