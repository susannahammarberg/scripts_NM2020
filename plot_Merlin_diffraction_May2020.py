# -*- coding: utf-8 -*-
"""
Created on May 2 2020
On the experiment on NanoMAX May 2020

Run this in the py36 environment where hdf5plugin is installed

Created from a copy from plot_Merlin_diffraction (from 2017)

Use section per section

Read in the Merlin  data.
Plot rocking curve 
save Merlin (so far) data as sparce matrices
Masking and choosing roi of data. 
Bright field, dark field, DPC analysis. 


@author: Susanna Hammarberg

"""


import sys   #to collect system path ( to collect function from another directory)

from numpy import fft
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import hdf5plugin #for 2020 data
import h5py
from math import floor
from math import ceil
#import scipy
from scipy import misc # to imoport image
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse 
from sys import getsizeof   #se hur mkt minne variabler tar upp  
#from mayavi import mlab
import matplotlib.animation as animation
import itertools as it # used to iterate in 2 staes in a loop
import time
date_str = time.strftime("%Y%m%d") # -%H%M%S")

def create_mask_Merlin():
    # Alex mask:
    mask_dir = 'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/'
    file_name = 'merlin_mask_200430_8keV.h5'
    with h5py.File(mask_dir  + file_name,'r' ) as fp:
        mask = np.array(fp['mask'])
#
#for i in range(0,2):
#    a[ a== a.max()] = 0
#    
    return mask

# Choose mask: gather mask or make mask
mask_Merlin = create_mask_Merlin()


#%% plot 2020
#Plot single diffraction pattern (in multiple scans if you want)
# original scans

#scans = np.arange(429,490)

#scans = np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ])

#original scans full range USE THIS WITH CURRENT SHIFTING LIST
scans = np.arange(429,490).tolist()

##mesh stepscan a bragg ptycho set that was aborted
#scans = np.arange(381,418).tolist()

scans = [340]
directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/' 
#frame = 435 # 31
for scan in scans:
    with h5py.File(directory + '000' + str(scan) +'.h5','r' ) as fp:
        data = np.array(fp['entry/measurement/merlin/frames'][frame][0:300])
        gonphi = np.array(fp['entry/snapshot/gonphi'])

    plt.figure()
    plt.imshow((mask_Merlin[0:300] * data),cmap='plasma')
    plt.title('Scan %d position %d. Masked'%(scan,frame))       
    plt.colorbar()
    #plt.savefig(r"C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\diffraction_position_%d\gonphi"%(frame) + \
#                ('%.2f'%(gonphi)).replace(".","_") + '_'  + 'S' + str(scan) + '_' + date_str)
    plt.close()

#%% Plot all data in one scan
#scans = np.arange(284,359).tolist() + np.arange(364,373).tolist()
# start at scan nbr:
scan = 537#330#320#300#33# 395#446#471 # 429


raw_slicey = slice(0,300)
#raw_slicex = slice(0,-1)
raw_slicex = slice(261,-1)

directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/000%d'%scan
with h5py.File(directory  + '.h5','r' ) as fp:
    data = np.array(fp['entry/measurement/merlin/frames'][:,raw_slicey,raw_slicex])
    #data = np.array(fp['entry/measurement/merlin/frames'])

#InP360? GaInP180?
#GaInP
xcen = 155
ycen = 148
#InP
xcen = 364
ycen = 148
size = int(256/2)
slicey = slice(ycen-size, ycen+size)
slicex = slice(xcen-size, xcen + size)

plt.figure()
plt.imshow(np.log10(mask_Merlin[raw_slicey,raw_slicex]*sum(data)), cmap='magma', interpolation='none')#magma
plt.title('Scan %d all diffraction. Masked'%(scan))       
plt.colorbar()
plt.show()

#plt.figure()
#plt.imshow((sum(mask_Merlin[slicey,slicex]*data[:,slicey,slicex])), cmap='magma', interpolation='none')#magma
#plt.title('Scan %d all diffraction. Masked'%(scan))       
#plt.colorbar()

#for ii in range(31):  
#    plt.figure()
#    plt.imshow(mask_Merlin[0:300]*(data[ii]), cmap='plasma')
#    plt.title('Scan %d all diffraction. Masked'%(scan))       
#    plt.colorbar()
##plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction_central_Bragg_slice\%s' + date_str) %(scan)
#
plt.ioff()

#%%
#do BF with this data
# function for calculating bright field of data
# input is the data as a 3D array with dimension data[scans_y * scans_x][yPixels][xPixels]
def bright_field(data,y,x):
    index = 0
    photons = np.zeros((y,x)) 
    for row in range(0,y):
        print( row)
        for col in range(0,x):
            photons[row,col] = np.sum(data[index]) #/ max_intensity
            
            #fig, ax = plt.subplots(nrows=1)
            #ax.imshow(np.log10(photons),cmap = 'jet', interpolation='none', extent= [ 0, (x-1)*50E-3,(y-1)*50E-3,0 ])
            #plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]',title=r'Bright field #S440 Position %d'%index)
            ##plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\BF_maps\before_shifting\updated_20210609\S#440_pos%s'%("{0:03}".format(index))) 
            index = index+1
            print(index)
    return photons
#y = 13#13
#x = 21#67

y= 13
x= 67#80
bf_map = bright_field((mask_Merlin[raw_slicey,raw_slicex]*data),y,x)

fig, ax = plt.subplots(nrows=1)
ax.imshow((bf_map),cmap = 'jet', interpolation='none', extent= [ 0, (x-1)*50E-3,(y-1)*50E-3,0 ])
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]')#,title=r'Bright field #S440')
#plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\BF_maps\before_shifting\updated_20210609\BF_S#440') 
plt.show()

#%%
#Make a rocking curve in a single position or all postions

rocking_curve = []
gonphis = []

# original scans
scans = np.arange(429,490)


# this is ABOUT right
#scans = np.arange(284,359).tolist() + np.arange(364,373).tolist()
#scans = np.arange(249,264).tolist()

#scans = np.arange(507,568+1).tolist()

#scans = np.arange(381,418).tolist()

#scans = np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ])
directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/' 
position = 435#160 #210# 435 #max 871
for scan in scans:
    with h5py.File(directory + '000' + str(scan) +'.h5','r' ) as fp:
        print(scan)
        #sx = np.array(fp['entry/measurement/sx'])
        #I0= np.array(fp['entry/measurement/alba/channel1'])
        #data = np.array(fp['entry/measurement/merlin/frames'][position,raw_slicey,raw_slicex])
        data = np.array(fp['entry/measurement/merlin/frames'][:,raw_slicey,raw_slicex])
        #data = np.array(fp['entry/measurement/merlin/frames'])
        rocking_curve.append(np.sum(data*mask_Merlin[raw_slicey,raw_slicex]))
        gonphi = np.array(fp['entry/snapshot/gonphi'])
        gonphis.append(gonphi)
    
plt.figure()
#plt.plot(gonphis, (rocking_curve))
#plt.xlabel('Gonphi')
plt.plot(np.arange(len(scans)), (rocking_curve),'.-')
#plt.yscale('log')
plt.xlabel('Scan #')
#locs, labels = plt.xticks()  # Get the current locations and labels.
plt.xticks(np.arange(len(scans)),scans)  # Set label locations.
plt.xticks(rotation=70)

plt.title('Rocking curve InP %d scan %d -'%(position,scans[0]))       
#plt.title('Rocking curve all positions scan %d -'%(scans[0]))      
plt.show() 
#plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\rocking_curve_pos_%d_scans_'%position + date_str)
#plt.title('Rocking curve all positions scan %d - '%(scans[0]))       
##plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\rocking_curve_allpos_stepscan_set' + date_str)


#%%

#Plot the central position in height. sum over gonphi. plot as a function of x (NW lng axis)

#original scans full range USE THIS WITH CURRENT SHIFTING LIST
scans = np.arange(429,490).tolist()
directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/' 

data_list = []
for scan in scans:
    with h5py.File(directory + '000' + str(scan) +'.h5','r' ) as fp:
        data_list.append(np.array(fp['entry/measurement/merlin/frames'][402:402+67][0:300]))

data_np = np.array(data_list)

row=0
for ii in range(0,872):
    if ii%67==0:
        print(ii)
        print('row %s (last frame)'%row)
        row+=1

for x_frame in range(0,67):
    plt.figure()
    plt.imshow( np.log10(mask_Merlin * sum(data_np)[x_frame]), cmap='plasma')
    plt.title('long axis position %d'%(x_frame))       
    plt.colorbar()
    plt.savefig(r"C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction\diffraction_along_long_axis\row5\xpos" + "{:03d}".format(x_frame))
    plt.close()
    print(x_frame)





#%%
















#Load data


# exp. parameters for conversion between detector and object plane
energy = 9.49   #keV   
wavelength = (1.23984E-9)/energy
pixel_det_Pil100K = 172E-6   # Pixel ctr to ctr distance (w) [m] #Pilatus
#z_det = 4.2   # Pil100K
#For Merlin
pixel_det = 55E-6   # Pixel ctr to ctr distance (w) [m] #Merlin
z_pilK100 = 4.2 #correct?
# object-detector distance? JWs skript säger 0.7
z_Merlin = 1.065  # inte rätt för söndagsmacrot
#z_pil1M = 1.235
epsilon = 2.220446049250313e-16



# ev matrix to rad in Pil1M data
#diffSet=np.zeros((nbr_positions, 1043, 981))  # Pil1M

# Allocate memory Pil 100K
# TODO_ cant have this it is to large
diffSet_Pil100K = np.zeros((2,nbr_rows,nbr_cols, 195, 487),dtype=np.int32)  


# Allocate memory Merlin
# u dont allocate in python! (not even np?)

# load metadata to gather motorpositions in for loop
metadata = h5py.File( metadata_directory)
# create motorpostion array to store for each scan
motorpositions_gonphi=np.zeros((nbr_rotations))
# y is sampled once per line (see JW mail called 'SV: data')
#motorpositiony = np.zeros((nbr_rotations,nbr_rows))
## there is no samsx for flyscnas (there is one but i dont know what it means. Use 'adlink_buff' for x )
#motorpositionx = np.zeros((nbr_rotations))

row_Merlin = []
list_Merlin = []# [ [], [], [], [], [] ]            #make lists inside a list li. alt kan man skriva list_Merlin.append([]) i for loopen
motorpositiony = []
motorpositionx = []
#tuple_Merlin = ((nbr_rows), )    
# maby should use tuple istead of list..
# index number called rotation
rotation = 0
# vetor defining how many pixels the diffraction pattern is shifted with (for Merlin #S458-461,496-515)
#vertical_shift = [-1,-1,0,0,0,  0,0,2,1,0,  1,1,1,0,-1,  -1,-1,-1,-1,0,  -1,-1,0,0,1,  1,-1,0,1,0,   2,0,0,1,1,  1,0,0,1,1,  1,2,2,2,4,  3,3,3,3,3,   3];

# read in and save data ROIs + motorpositions. mask data   #scan_name_int+1
for scan_number in it.chain(range(scan_name_int, scan_name_int + nbr_rotations)):#scan_name_int+1#488), range(496, 515)):
         
    # define list to save all data from 1 rotation(all rows, all postitions):
    temp_list = []

    # loop over all 16 flyscans that constitute the set
    for row in range(0, nbr_rows):  
        
        # load hdf5-file
        scan = h5py.File( directory  + str('{0:04}'.format(row)) + '.hdf5','r') # read-only

        # load and store in np array
        data_Merlin =  np.array(scan.get('/entry_0000/measurement/Merlin/data' ), dtype=np.int32)

        # Flip up and down
        # flipping the diffpatterns en och en, maybe there is an easier way
        for col in range(0, nbr_cols):
            data_Merlin[col] = np.flipud(data_Merlin[col])

        # mask the Merlin data
 #       data_Merlin = data_Merlin * mask_Merlin
        # Remove 3 hot pixels (for nice plots):
        data_Merlin[:,312,321] = data_Merlin[:,312,322] # + data_Merlin[:,312,322] + + / 4
        data_Merlin[:,294,301] = data_Merlin[:,294,302]
        data_Merlin[:,191,231] = data_Merlin[:,191,232]

        # remove the kors from the middle of the merlin detector images
#        data_Merlin[:,255,:] = data_Merlin[:,254,:]
##        #data_Merlin[:,255,:] = (data_Merlin[:,254,:] + data_Merlin[:,257,:]    )/1
##        #np.disp((data_Merlin[:,254,:] + data_Merlin[:,257,:]    ))
##        #data_Merlin[:,256,:] = (data_Merlin[:,254,:]  +data_Merlin[:,257,:]      )/2
#        data_Merlin[:,256,:] = data_Merlin[:,257,:]
#        data_Merlin[:,:,255] = data_Merlin[:,:,254]    
#        data_Merlin[:,:,256] = data_Merlin[:,:,257]    
        
        # (gör det någon skillnad om jag gör maskningen och väljer roi på denna rad istället och skippar att 
        # spara det i en 'vanlig' matris?)
        # select roi on the detector   (for the last macro)
#        roi_y = vertical_shift[rotation] + 100, vertical_shift[rotation] + 300
#        data_Merlin = data_Merlin[:, roi_y[0]:roi_y[1], 100:380] #       data_Merlin = data_Merlin[:,130:250,250:360]
#        select roi of Merlin(other scans)
        roix_start = 100# 100# 100 #140
        roix_end = 350#380#380 #340
        roiy_start = 0# 70 #0 # testar att klippa helt fel
        roiy_end = 300#130# 200
#        roix_start = 160#100# 100# 100 #140
#        roix_end = 340#350#380#380 #340
#        roiy_start = 70# 70 #0 # testar att klippa helt fel
#        roiy_end = 250
        data_Merlin = data_Merlin[:, roiy_start:roiy_end, roix_start:roix_end]
        #data_Merlin = data_Merlin[:, 70: 250, 150:330]
        #data_Merlin = data_Merlin[:, 130: 200, 200:280]
        # hade denna innan lunch
        #data_Merlin = data_Merlin[:, 140: 190, 210:270]
        plt.figure()
        plt.imshow(data_Merlin[10])
        # save all images as sparse matrices in a list M
        one_row_sparse_Merlin = [sparse.lil_matrix(data_Merlin[i]) for i in xrange(nbr_cols)]
        
        # lägg till M i en list för varje rad
        temp_list.append(one_row_sparse_Merlin)

#        #load and save transmission data from pil100K:
#        scan = h5py.File( directory_pil100K  + str('{0:04}'.format(row)) + '.hdf5','r') # read-only
#        data_pil = scan.get('/entry_0000/measurement/Pilatus/data' ) #pilatus data
#        diffSet_Pil100K[rotation][row] = np.array(data_pil)
##
##        
        #good or bad to delete inside loop? gets overwritten if not.
        #del scan, data_Merlin, data_pil
    


    # gather motor postions from metadata (h5-file). one file for each scan #S,not one file for each flyscan
    motorpositions_directory = '/entry%s' %scan_name_string  
    
    dataset_motorposition_gonphi = metadata.get(motorpositions_directory + '/measurement/gonphi')      
    dataset_motorpositiony = metadata.get(motorpositions_directory + '/measurement/samy')
    # instead of samx, you find the motorposition in flyscans from 'adlink_buff'
    dataset_motorpositionx = metadata.get(motorpositions_directory + '/measurement/AdLinkAI_buff') 
    
    # TODO: add gonphi
    motorpositions_gonphi[rotation] = np.array(dataset_motorposition_gonphi)
    
    motorpositiony.append(np.array(dataset_motorpositiony))
    ee=np.array(dataset_motorpositiony)
    # Load in x motorpositions from AdLinkAI_buff (contains additional zeros)
    motorpositionx_with_zeros = np.array(dataset_motorpositionx)
    # store the positions without the zeros
    motorpositionx.append(motorpositionx_with_zeros[:,0:nbr_cols])

    
    
    # save the whole diffraction-shebank (append one whole rotation )
    list_Merlin.append(temp_list)

    rotation = rotation + 1
    # update directories to gather next scan
    scan_name_string = '%d' %scan_number
    
    #d#irectory = 'F:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX29A_NW1/scan_0%d_merlin_'%scan_number # sundaymacro
  #  directory = 'F:/exp20170628_Wallentin_nanomax/exp20170628_Wallentin/JWX33_NW2/scan_0%d_merlin_'%scan_number 
        # inte så snyggt att deklarera om directory?
        
    np.disp('rotation:')
    np.disp(rotation)

del scan, data_Merlin, one_row_sparse_Merlin, temp_list, motorpositionx_with_zeros
#del data_pil
del dataset_motorpositionx, dataset_motorpositiony, dataset_motorposition_gonphi

dataset_title = metadata.get(motorpositions_directory + '/title')
command = np.array(dataset_title)
# this is for the SECOND rotation
print( command)
fig = plt.figure()
plt.imshow(np.log10(abs(diffSet_Pil100K[0][0][0])), cmap = 'jet', interpolation = 'none')
plt.colorbar()
# need to select roi for mask aswell as i am going to send it to ePIE 
mask_Merlin = mask_Merlin[roiy_start:roiy_end, roix_start:roix_end] 
       
fig = plt.figure()
plt.imshow(np.log10(list_Merlin[0][9][15].toarray()), cmap = 'jet', interpolation = 'none')
plt.colorbar()

#%%

## test 'unsparse' diffraction matrix
#J=[m.toarray() for m in M]
# for a single diffraction pattern:
#list_Merlin[0][0][0].toarray()
#plt.figure()
#plt.imshow(np.log10(J[9]))



def create_mask_Pil100K():
##    probe_mask = np.ones((diffSet.shape[1],diffSet.shape[2]))
    # Find all cold pixels (they show minus values)
    # my looking at the sum of all diffraction patterns from one random scan (the first)
    sumTotalDiffSet= sum(sum(diffSet_Pil100K[0]))
    probe_mask = sumTotalDiffSet > 0
    
    # remove too high intensity pixels. 
    probe_mask[84:116, 225:258] = 0 
    probe_mask[152,229] = 0 
    
    probe_mask = 1 - probe_mask
#    j=239
#    probe_mask[111,j:245] = 0 
#    probe_mask[112,j:245] = 0 
#    probe_mask[113,j:245] = 0  
#    probe_mask[114,j:245] = 0 
#    probe_mask[115,j:245] = 0 
#    probe_mask[113,242] = 0 #pixel with highest intensity
    return probe_mask

mask_Pil100K = create_mask_Pil100K()


#probe_mask_Pil100K = probe_mask_Pil100K[roiy_Pil_start:roiy_Pil_end, roix_Pil_start:roix_Pil_end] 
## apply PilK100 mask


#objectFunc, probe, ani, sse, psi, PRTF = ePIE(k, diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans)
def bright_field_analysis(data):
    
    photons = np.zeros((nbr_rows,nbr_cols)) 
    #max_intensity = np.sum(  np.sum(data,axis=1) , axis=1).max()   # sum over rows and columns not sum over different diffPatterns
    for row in range(0,nbr_rows):
        for col in range(0,nbr_cols):
            photons[row][col] = sum(sum(data[row,col,:,:])) #/ max_intensitycol
                
    return photons

def diff_phase_contrast(data):
    tempy = 0
    tempx = 0
    
    diff_phasey = np.zeros((nbr_rows,nbr_cols))
    diff_phasex = np.zeros((nbr_rows,nbr_cols))
    pol_DPC_r = np.zeros((nbr_rows,nbr_cols))
    pol_DPC_phi = np.zeros((nbr_rows,nbr_cols))
    rem_bkg_x = np.zeros((nbr_rows,nbr_cols))
    rem_bkg_y = np.zeros((nbr_rows,nbr_cols))
        
    for row in range(0, nbr_rows):
        for col in range(0, nbr_cols):
            
            # m and n pixels in the diffraction data
            for m in range(0, data.shape[2]):
                for n in range(0, data.shape[3]):
                    tempy = tempy + (m-nbr_rows/2) * data[row, col, m, n] #/ (diffSet[index, m, n]+ 2.220446049250313e-16)
                    tempx = tempx + (n-nbr_cols/2) * data[row, col, m, n]
            # spara värdet på den första pixeln:
            # detta känns onödigt krävande för då måste if satsen kollas varje gång fast jag vet vilket k jag vill ha
#            if row == 0 and col == 0:
#                bkg_x = tempx
#                bkg_y = tempy
            sum_diffSet = sum(sum(data[row,col]))
            diff_phasey[row, col] = tempy / sum_diffSet
            diff_phasex[row, col] = tempx / sum_diffSet
            rem_bkg_x[row,col] = diff_phasex[row,col] -241#- bkg_x # 68.25
            rem_bkg_y[row,col] = diff_phasey[row,col] - 92#- bkg_y # 62.2
            # DPC in polar coordinates. r then phi:
            pol_DPC_r[row, col] = np.sqrt( (rem_bkg_x[row,col])**2 + (rem_bkg_y[row,col])**2)    
            pol_DPC_phi[row, col] = np.arctan( rem_bkg_y[row,col] / rem_bkg_x[row,col])
            tempy = 0
            tempx = 0
        print( row)     
    return diff_phasex, diff_phasey, pol_DPC_r, pol_DPC_phi

#################################################################################################################################################################OBS!
diffSet_Pil100K = diffSet_Pil100K * ( 1 - mask_Pil100K)

# Choose which detector images to analyze:
# if sats so that it plots the right detector name       
detector = 2
rotation_analysis = 0
if (detector == 0):                                               

    #bright_field = bright_field_analysis(list_Merlin-...........[rotation_analysis])
    #dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(list_Merlin.........)
    analyse_detector_name_string ='Merlin'    
elif (detector == 1): 
    bright_field = bright_field_analysis(diffSet_Pil100K[rotation_analysis])
    dpc_x, dpc_y, pol_DPC_r, pol_DPC_phi = diff_phase_contrast(diffSet_Pil100K[rotation_analysis])
    analyse_detector_name_string ='Pil100K'    
else:
    print 'no analysis'
#    
def plot_analysis():
    
    plt.figure()
    plt.imshow(bright_field, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Bright field on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    #plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_transm'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')

     
    plt.figure()
    plt.imshow(dpc_x, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Horizontal DPC on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')  
    plt.colorbar()
    #plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCx'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')
    
    plt.figure()
    plt.imshow(dpc_y, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: Vertical DPC on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()
    #plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCy'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')
    
    plt.figure()
    plt.imshow(pol_DPC_r, cmap='gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC r on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    plt.xlabel('Nominal motorpositions [um]')
    plt.colorbar()
    #plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCpol_r'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')

    plt.figure()    
    plt.imshow(pol_DPC_phi, cmap = 'gray', interpolation='none')#, extent=[motorpositionx[0], motorpositionx[-1], motorpositiony[0], motorpositiony[-1] ])
    plt.title('Scan %d: DPC phi on %s'%((rotation_analysis+first_scan_nbr),analyse_detector_name_string))
    #plt.xlabel('Nominal motorpositions [um]')
    #plt.ylabel('Nominal motorpositions [um]')
    plt.colorbar()  
    #plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d_DPCpol_phi'%(rotation_analysis+first_scan_nbr), bbox_inches='tight')

# if u ran the analysis, plot
if (detector == 0) or (detector == 1): 
    plot_analysis()


roix_Pil_start = 100#150# 100# 100 #140
roix_Pil_end = 380#380#380 #340
roiy_Pil_start = 0# 70 #0 # testar att klippa helt fel
roiy_Pil_end = 200#130# 200
# choose roi on Pil100K detector 
diffSet_Pil100K = diffSet_Pil100K[:,:,:, roiy_Pil_start:roiy_Pil_end, roix_Pil_start:roix_Pil_end]
mask_Pil100K = mask_Pil100K[roiy_Pil_start:roiy_Pil_end, roix_Pil_start:roix_Pil_end]

plt.figure()
plt.imshow(np.log10(sum(sum(diffSet_Pil100K[0]))))

# När jag skickar in bilderna här är de redan klippta i roi
#def COM(rot):
roiNx= diffSet_Pil100K.shape[4]
roiNy= diffSet_Pil100K.shape[3]
    # define a vector with length of the length of roi on the detector
roix = np.linspace(1, roiNx, roiNx)
## define a vector with length of the height of roi on the detector
roiy = np.linspace(1,roiNy, roiNy)
# meshgrids for center of mass calculations
X, Y = np.meshgrid(roix,roiy)

COM_hor = np.zeros((nbr_rows,nbr_cols))
COM_ver = np.zeros((nbr_rows,nbr_cols))
COM_mag = np.zeros((nbr_rows,nbr_cols))
COM_ang = np.zeros((nbr_rows,nbr_cols))
rot = 0
for row in range(0,nbr_rows):
    for col in range(0,nbr_cols):
        COM_hor[row,col] = sum(sum(diffSet_Pil100K[rot][row][col]*X))/sum(sum(diffSet_Pil100K[rot][row][col]))
        COM_ver[row,col] = sum(sum(diffSet_Pil100K[rot][row][col]*Y))/sum(sum(diffSet_Pil100K[rot][row][col]))
        if row == 0 and col == 0:
            bkg_hor = 152.4#COM_hor[row,col] 
            bkg_ver = 101.8#COM_ver[row,col] 
                    # DPC in polar coordinates. r then phi:
        COM_mag[row, col] = np.sqrt( (COM_hor[row,col]-bkg_hor)**2 + (COM_ver[row,col]-bkg_ver)**2) 
 #       COM_ver(row_ROI(1):row_idx, col_ROI(1):col_ROI(end),scan_idx)-mean(mean(COM_ver(row_ROI(1):row_idx,col_ROI(1):col_ROI(end),scan_idx))))
#[COM_angle(row_ROI(1):row_idx, col_ROI(1):col_ROI(end),scan_idx), COM_magnitude(row_ROI(1):row_idx,col_ROI(1):col_ROI(end),scan_idx)] = cart2pol(COM_hor(row_ROI(1):row_idx,nbr:cols,scan_idx)-mean(mean(COM_hor(row_ROI(1):row_idx,col_ROI(1):col_ROI(end),scan_idx))), grej        
        COM_ang[row, col] = np.arctan( COM_hor[row,col] / COM_ver[row,col])
   #     COM_hor_tmp(col) = sum(sum(image_dat(roiy,roix).*X))/sum(sum(image_dat(roiy,roix)));
   #     COM_ver_tmp(col) = sum(sum(image_dat(roiy,roix).*Y))/sum(sum(image_dat(roiy,roix)));
        
del rot
plt.figure()
plt.imshow(COM_hor)
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

plt.rc('font', **font)
#plt.rc('ytick', labelsize=20)
#plt.rc('xtick', labelsize=20) 
plt.colorbar()

plt.figure()
plt.imshow(COM_ver)
plt.colorbar()
plt.figure()
#plt.imshow(COM_mag, cmap='jet', interpolation='none', extent=[motorpositionx[0][0][0], motorpositionx[0][0][nbr_cols-1], motorpositiony[0][0], motorpositiony[0][nbr_rows-1] ])
# till jepser utan motorpositioner:
plt.imshow(COM_mag, cmap='jet', interpolation='none', extent=[motorpositionx[0][0][0]-motorpositionx[0][0][0], motorpositionx[0][0][nbr_cols-1]-motorpositionx[0][0][0], motorpositiony[0][0]-motorpositiony[0][0], motorpositiony[0][nbr_rows-1] -motorpositiony[0][0]])
plt.colorbar()
plt.xlabel('[$\mu m$]')
plt.ylabel('[$\mu m$]')
plt.title('Center of mass')
plt.figure()
plt.imshow(COM_ang)
plt.colorbar()

def variables_scatterplot():
    #this is defined for scatter plot:
    x = np.linspace( first_scan_nbr , first_scan_nbr + nbr_rotations, nbr_rotations)
    y = np.linspace( 130, 250, 120)
    z = np.linspace( 250, 360, 110)
    X, Y, Z = np.meshgrid(y,x,z) # rörigt!
    X1, Y1 = np.meshgrid(z,y)
#variables_scatterplot()

# construct a factor I0 to compensate for intensity variation in the beam 
# where I0 = 1 is the largest value
def intens_norm_PilK100(rotation):
    
    I0 = np.zeros((nbr_rows, nbr_cols))
    for row in range(0,nbr_rows):
        for col in range (0, nbr_cols):
            I0[row][col] = sum(sum(diffSet_Pil100K[rotation][row][col]))
    I0 = I0/I0.max()
    return I0
I0 = intens_norm_PilK100(0)

# TODO: this can be for both Bragg and transmission. factor np.cos(2theta) blir ju ett för vinkel 0    
#def preps_phaseretival_2DBragg(rotation):

# each scan is on a grid of 16 rows (16 flyscans) and 101 cols 
# in total 16*101=1616 diffraction patterns. Insert these in to a 
# numpy matrix like diffSet, that I used before.
def run_ePIE(detector):
    index = 0
    if (detector == 0):
        # Run ePIE for Flyscans    
        #def run_ePIE_2D(rotation):
        rotation = 0
        # theta should be read out from gonphi for each particular scan=?
        #theta_rad = 11.1*np.pi/180 #arb
        #print theta_rad
        theta_deg = 180 - motorpositions_gonphi[rotation]
        theta_rad = (180 - motorpositions_gonphi[rotation])*np.pi/180
        print theta_deg
        #theta_rad = 11.1*np.pi/180 #arb
        # Sizes of roi of diffraction patterns (pixels used on the Merlin detector)
        Ny = list_Merlin[0][0][0].shape[0]      
        Nx = list_Merlin[0][0][0].shape[1]     
        
        # size of one pixel in objectplane
        xpixel = z_Merlin*wavelength/ (Nx * pixel_det)
        ypixel = z_Merlin*wavelength/ (Ny * pixel_det)
        
        # what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
        sizeDiffObjectx = Nx * xpixel
        sizeDiffObjecty = Ny * ypixel
           
        # hur långt motorn rör sig i x och yled: 
        # TODO: reda ut om det ska va 2 theta elelr vad. kolla ptypy kod.
        motorpositiony[rotation] = np.cos(1*theta_rad)* motorpositiony[rotation]
        motorpositionx[rotation] = np.cos(1*theta_rad)* motorpositionx[rotation]
        motorWidthx = ( motorpositionx[rotation].max() - motorpositionx[rotation].min() ) * 1E-6
        motorWidthy = ( motorpositiony[rotation].max() - motorpositiony[rotation].min() ) * 1E-6
        
        # so the size of the object function should be enough to contain: (but actually it becomes a little bit larger because i have to round of to a hole pixel)
        objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
        objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx
        
        # so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
        # should i use ceil!? or floor?
        objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
        objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)
        
        # allocate memory for object function
        objectFunc = np.zeros((int(objectFuncNy), int(objectFuncNx)))
        
        ####
        # handle motorpositions so that they are in a format that work in the ePIE algorithm
        # They should be in a vector, not a matrix, with length equal to the amount of grid positions
        #####
        # the y motorposition only has one datapoint per row so i copy it (in a confusing but good way)
        positiony = np.ones((nbr_cols,nbr_rows))
        positiony = (positiony * motorpositiony[rotation]).transpose()
        positiony = positiony.ravel()
        # convert the xpositions, which are in a matrix, to a long vector
        positionx = motorpositionx[rotation].ravel()
        
        # 'normalized' motorpostions converted to meters
        positiony = (positiony - positiony.min() ) *1E-6
        positionx = (positionx - positionx.min() ) *1E-6
        # calculate the total number of 'grid' points
        nbr_scans = nbr_cols*nbr_rows

        diffSet = np.zeros((nbr_scans,Ny, Nx))        
        for row in range(0,nbr_rows):
            for col in range(0,nbr_cols):
                # compensate for intensity variation with the factor I0
                diffSet[index] = list_Merlin[rotation][row][col].toarray() #/ I0[row][col] 
                index = index + 1
     


    elif (detector ==1):
        # Run ePIE for Flyscans    
        #def run_ePIE_2D(rotation):
        rotation = 0

        # Sizes of roi of diffraction patterns (pixels used on the Merlin detector)
        Ny = diffSet_Pil100K[0][0][0].shape[0]      
        Nx = diffSet_Pil100K[0][0][0].shape[1]     
        
        # size of one pixel in objectplane
        xpixel = z_pilK100*wavelength/ (Nx * pixel_det_Pil100K)
        ypixel = z_pilK100*wavelength/ (Ny * pixel_det_Pil100K)
        
        # what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
        sizeDiffObjectx = Nx * xpixel
        sizeDiffObjecty = Ny * ypixel
           
        # hur långt motorn rör sig i x och yled: 
        motorpositiony[rotation] = motorpositiony[rotation]
        motorpositionx[rotation] = motorpositionx[rotation]
        motorWidthx = ( motorpositionx[rotation].max() - motorpositionx[rotation].min() ) * 1E-6
        motorWidthy = ( motorpositiony[rotation].max() - motorpositiony[rotation].min() ) * 1E-6
        
        # so the size of the object function should be enough to contain: (but actually it becomes a little bit larger because i have to round of to a hole pixel)
        objectFuncSizeMaxy = motorWidthy + sizeDiffObjecty
        objectFuncSizeMaxx = motorWidthx + sizeDiffObjectx
        
        # so with a pixel-size of xpixel * ypixel, the obect function should be this many pixels:
        # should i use ceil!? or floor?
        objectFuncNy = ceil(objectFuncSizeMaxy / ypixel)
        objectFuncNx = ceil(objectFuncSizeMaxx / xpixel)
        
        # allocate memory for object function
        objectFunc = np.zeros((int(objectFuncNy), int(objectFuncNx)))
        
        ####
        # handle motorpositions so that they are in a format that work in the ePIE algorithm
        # They should be in a vector, not a matrix, with length equal to the amount of grid positions
        #####
        # the y motorposition only has one datapoint per row so i copy it (in a confusing but good way)
        positiony = np.ones((nbr_cols,nbr_rows))
        positiony = (positiony * motorpositiony[rotation]).transpose()
        positiony = positiony.ravel()
        # convert the xpositions, which are in a matrix, to a long vector
        positionx = motorpositionx[rotation].ravel()
        
        # 'normalized' motorpostions converted to meters
        positiony = (positiony - positiony.min() ) *1E-6
        positionx = (positionx - positionx.min() ) *1E-6
        
        # calculate the total number of 'grid' points
        nbr_scans = nbr_cols*nbr_rows


        diffSet = np.zeros((nbr_scans,Ny, Nx))
        
        for row in range(0,nbr_rows):
            for col in range(0,nbr_cols):
                # compensate for intensity variation with the factor I0 
                diffSet[index] = diffSet_Pil100K[rotation][row][col]
                index = index + 1
       
    return diffSet, objectFuncNx, objectFuncNy, ypixel, xpixel, positionx, positiony, objectFunc, nbr_scans, sizeDiffObjecty, sizeDiffObjectx
# choose detector from which to reconstruct
detector = 0   
diffSet, objectFuncNx, objectFuncNy, ypixel, xpixel, positionx, positiony, objectFunc, nbr_scans, sizeDiffObjecty, sizeDiffObjectx = run_ePIE(detector)

if (detector == 0):
    detector_name_string ='Merlin'  
    mask = mask_Merlin
elif (detector == 1):
    detector_name_string ='Pil100K'  
    mask = mask_Pil100K
    

#nbr of iterations
k = 1
# probe construction: create gaussian or load saved probe

probe = np.load('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan195\saved_probe\scan195_probe_k100.npy')
#probe=np.zeros((diffSet.shape[1], diffSet.shape[2]))
#sigmay = 14.1#135# 14.1# 14.1               # initial value of gaussian height     #Scan51 2 x 2 
#sigmax = 10#12#5# 10                    # initial value of gaussian width
#probe = create2Dgaussian( sigmay, sigmax, diffSet.shape[1], diffSet.shape[2])
plt.figure()
plt.imshow(abs(probe), cmap='jet')
plt.colorbar()


##############################################RÄTT MASK OCKSÅ
        
# run 2d ePIE
objectFunc, probe, ani, sse, psi, PRTF, save_G = ePIE(k, diffSet, probe, objectFuncNy, objectFuncNx, ypixel, xpixel, positiony, positionx, nbr_scans, mask)

variabel=12
#variabel = '%d' %variabel 
# save the probe
#np.save('scan195_probe_k500.npy', probe)
#
plt.figure()
plt.subplot(121)
plt.imshow(np.log10((abs(save_G))**2), cmap = 'jet', interpolation = 'none')
#plt.title('reconstructed')
plt.colorbar()
plt.subplot(122)

plt.imshow(np.log10(abs(diffSet[178])), cmap = 'jet', interpolation = 'none')
plt.colorbar()
plt.title('variable = %d'%variabel)
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d\scan%d_compareDiffPat_k%d_%s_%d'%(scan_name_int,scan_name_int, k, detector_name_string,variabel), bbox_inches='tight')
#plt.title('original=%s hu')%variabel
np.save('scan195_obj_k100.npy', objectFunc)
np.save('scan195_probe_k100.npy', probe)

    #, origin="lower"                         # sets the scale on axes.
plotting_rangex = np.linspace(0,objectFuncNx*xpixel*1E6, 349) #349=size of object in pixels
plotting_rangey = np.linspace(0,objectFuncNy*ypixel*1E6, 395) #349=size of object in pixels
#y 395
plotting_rangex = plotting_rangex[0:85]#S195[21:125]#S211[61:180]
plotting_rangey = plotting_rangey[175:250]#S195[145:220]#S211[145:240],

plotting_rangex_start = plotting_rangex[0]
plotting_rangex_end = plotting_rangex[ len(plotting_rangex)-1]
plotting_rangey_start = plotting_rangey[0]
plotting_rangey_end = plotting_rangey[ len(plotting_rangey)-1]

 
plt.figure() 
#plt.imshow( np.angle(objectFunc[105:270,40:220]), cmap='jet', interpolation='none', extent=[plotting_range, 0,objectFuncNy*ypixel*1E6])
#S195:plt.imshow( np.angle(objectFunc[145:220,21:125]), cmap='jet', interpolation='none', extent=[plotting_rangex_start, plotting_rangex_end,plotting_rangey_start, plotting_rangey_end])
plt.imshow( np.angle(objectFunc[175:250,0:85]), cmap='jet', interpolation='none', extent=[plotting_rangex_start, plotting_rangex_end,plotting_rangey_start, plotting_rangey_end])
#forscan211: plt.imshow( np.angle(objectFunc[145:240,61:180]), cmap='jet', interpolation='none', extent=[plotting_rangex_start, plotting_rangex_end,plotting_rangey_start, plotting_rangey_end])
#plt.gca().invert_yaxis() 
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Reconstructed phase')
#plt.title('Scan %d: Object phase %s'%(scan_name_int, detector_name_string))
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d\scan%d_Ophase_k%d_%s_%d'%(scan_name_int,scan_name_int, k, detector_name_string, variabel), bbox_inches='tight')
   
plt.figure()             #[120:270][40:240]                                               # horisontalt vertikalt. xpixel * size(objectfunc[xled])
#plt.imshow(abs(objectFunc[105:270,40:220]), cmap='jet', interpolation='none')#, extent=[0,objectFuncNx*xpixel*1E6, 0, objectFuncNy*ypixel*1E6])
#THIS#plt.imshow((abs(objectFunc[175:250,0:85])), cmap='jet', interpolation='none', extent=[plotting_rangex_start, plotting_rangex_end,plotting_rangey_start, plotting_rangey_end])

# test plot from 0 
plt.imshow((abs(objectFunc[175:250,0:85])), cmap='jet', interpolation='none', extent=[0, plotting_rangex_end-plotting_rangex_start,0, plotting_rangey_end-plotting_rangey_start])
##S195:plt.imshow(np.sqrt(abs(objectFunc[145:220,21:125])), cmap='jet', interpolation='none', extent=[plotting_rangex_start, plotting_rangex_end,plotting_rangey_start, plotting_rangey_end])
#forscan211: #plt.imshow((abs(objectFunc[145:240,61:180])), cmap='jet', interpolation='none', extent=[plotting_rangex_start, plotting_rangex_end,plotting_rangey_start, plotting_rangey_end])
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Reconstructed amplitude')
#plt.title('variable = %d'%variabel)
#plt.title('Scan %d: Object amplitude'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d\scan%d_Oamp_k%d_%s_%d'%(scan_name_int,scan_name_int, k, detector_name_string,variabel), bbox_inches='tight')

plt.figure()
#plt.imshow(((abs(probe[100:270,40:220]))), cmap='jet', interpolation='none', extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
plt.imshow(((abs(probe))), cmap='jet', interpolation='none', extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
#plt.title('variable = %d'%variabel)
plt.title('Probe amplitude')
#plt.title('Scan %d: Probe amplitude'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d\scan%d_Pamp_k%d_%s_%d'%(scan_name_int,scan_name_int, k, detector_name_string,variabel), bbox_inches='tight')

plt.figure()                                                            # horisontalt vertikalt
plt.imshow(np.angle(probe), cmap='jet', interpolation='none', extent=[ 0,sizeDiffObjectx*1E6, 0,sizeDiffObjecty*1E6])
plt.xlabel(' [$\mu m$]')
plt.ylabel(' [$\mu m$]')
plt.title('Probe phase')
#plt.title('variable = %d'%variabel)
#plt.title('Scan %d: Probe phase'%scan_name_int)
plt.colorbar()
plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d\scan%d_Pphase_k%d_%s_%d'%(scan_name_int,scan_name_int, k, detector_name_string, variabel), bbox_inches='tight')  
    
#    
plot_x = np.linspace(0,diffSet.shape[2]-1,diffSet.shape[2])*xpixel*1E6
plt.figure()
plt.plot(plot_x ,abs(probe.sum(axis=0)), 'b+:', label='data')                                    
#plt.plot(plot_x, gauss(xCol, *poptCol), 'r-', label='fit')
#plt.plot(plot_x, yFit, 'g', label='manual fit')
plt.xlabel(' [$\mu m$]')
plt.ylabel('Intensity')
#plt.title('Scan %d: Probe summed over all rows. FWHM: %f $\mu m$'%(scan_name_int,FWHM_col))     #horizontal line
#plt.legend()
   # plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_probe_row_lineplot_k%d'%(scan_name_int, k), bbox_inches='tight')

plot_y = np.linspace(0,diffSet.shape[1]-1,diffSet.shape[1])*ypixel*1E6
plt.figure()
plt.plot(plot_y, abs(probe.sum(axis=1)), 'b+:', label='data')  # vertical line
#plt.plot(plot_y, gauss(xRow, *poptRow),'r-', label='fit')
plt.legend() 
plt.xlabel(' [$\mu m$]')
#plt.title('Scan %d: Probe summed over all columns. FWHM: %f $\mu m$'%(scan_name_int,FWHM_row))
#    plt.savefig('dokumentering\Jespers_scans\savefig\scan%d_probe_col_lineplot_k%d'%(scan_name_int, k), bbox_inches='tight')

#plt.savefig('C:\Users\Sanna\Desktop\NanoMAX062017\Analysis\savefig\scan%d\scan%d_Original_reconstructed_diffPatt_nbr60_k%d'%(scan_name_int, scan_name_int, k), bbox_inches='tight')
#
#plt.suptitle('log10 reconstructed diffraction pattern (nbr 800)')  

#    return (objectFunc, probe)
# run ePIE 2D. use rotation as input
#objectFunc, probe = run_ePIE_2D(0)
##    #plt.gca().invert_yaxis() 

#plt.show()

# save the Bragg peak in file to open it in matlab
#scipy.io.savemat('C:/Users/Sanna/Desktop/NanoMAX062017/Bragg_peak_S458_.mat', mdict={'Braggpeak': one_position_roi})


# Look how the COM directly on the diffraction pattern varies with angle. not a good COM definition.
# remove?
def COM_variation(j, nbr_iter):

    for i in range (j,nbr_iter):
        xindex = np.argmax(np.sum(one_position[i],axis=0))
        yindex = np.argmax(np.sum(one_position[i],axis=1))
        reddot=np.zeros((512,512))
            
        # Make a centred line in x and y intersection at COM
        reddot[:,xindex] = 500000 
        reddot[yindex,:] = 500000 
        np.disp( xindex)
        plt.figure()
        noes  = ['spring', 'autumn']
        plt.imshow(np.log10(one_position[i]), cmap=noes[1] , interpolation = 'none')
        plt.imshow(np.log10(reddot))
        #plt.imshow(np.log10(one_position[1]), cmap = 'hot', interpolation = 'none')
        #plt.colorbar() funkar ej med flera imshows
        plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
        
#COM_variation(0,3)    

def plot_summed_patterns_one_point(row,col):
    # sum the diffrtaction patterns in one point
    summed_point = 0
    for j in range(0,nbr_rotations):
        summed_point  = summed_point + list_Merlin[j][row][col].toarray()  
    plt.figure()
    plt.imshow(np.log10(summed_point), cmap = 'hot', interpolation = 'none')
    plt.axis('off')
    plt.colorbar()
    plt.title('masked diffraction patterns sum of one position in %d scans'%nbr_rotations)
#plot_summed_patterns_one_point(8,49)   #inser col and row you want to plot


#            # anim TABORT
#            im = plt.imshow(abs(objectFunc), animated=True, interpolation='none', extent=[0,6.837770297837617,0,6.825238081022181])
#            
#                    ims.append([im])
#                    ani = animation.ArtistAnimation(fig, ims, interval=1500, blit=True,repeat_delay=2000)
#                    plt.show()
# plot diffraction patterns merlin + pilK100
def plotalot():
    
    
    #figure for animation
    #fig = plt.figure()
    # Initialize vector for animation data
    #ims = []    
    
    
    for col in range(0,82,1): #film går upp till 82 Detta är tidsaxeln på filmen
#        
#        # plot a single image (single rotations, single row and column)
##        plt.figure()
##        plt.imshow(np.log10(list_Merlin[i][8][49].toarray()), cmap = 'hot', interpolation = 'none')
##        plt.colorbar()
##        plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
##        
#        # plot the sum of the Bragg peak, summed in the 3 (?) direction        
#    summed_list_Merlin_col = 0
         summed_list_Merlin_1 = 0        
#        summed_list_Merlin_2 = 0        
#        summed_list_Merlin_3 = 0        
#        #TODO FIXA dessa tre plottar: (som JWs film)
         for rot in range(0,nbr_rotations-1):
            summed_list_Merlin_1 = summed_list_Merlin_1 + (list_Merlin[rot][6][col].toarray())
#        #plt.figure()    
#        #im = plt.imshow(np.log10(summed_list_Merlin_1), animated=True, cmap = 'jet', interpolation = 'none')
        #ims.append([im])
    # plotta enskilda bilder o ej animation
    plt.figure()
    plt.imshow(np.log10(summed_list_Merlin_1), cmap = 'jet', interpolation = 'none')
    plt.colorbar()
    
    #ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True,repeat_delay=200)  
    #plt.axis('off')
    #plt.show()
    # save animation:
    #ani.save('dynamic_images.mp4', writer="mencoder")
    # plot the first imjage and compare before and after the shifting

    # plot the sum of all diffraction pattterns for one row (middle of NW) for one rotation 
    #for row in range (0,nbr_rotations):
#    for col in range(0,nbr_cols):
#        summed_list_Merlin_col = summed_list_Merlin_col + (list_Merlin[0][6][col].toarray())
#        
#    plt.figure()
#    plt.imshow(np.log10(summed_list_Merlin_col), cmap = 'jet', interpolation = 'none')      
#    plt.colorbar()
            
# Dessa 2 ej rätt_
#    for j in range(0,101):
#        summed_list_Merlin_2 = summed_list_Merlin_2 + (list_Merlin[i][6][col].toarray())
#    plt.figure()    
#    plt.imshow(np.log10(summed_list_Merlin_2), cmap = 'jet', interpolation = 'none')
#    plt.colorbar()  
#
#    for j in range(0,101):
#        summed_list_Merlin_3 = summed_list_Merlin_3 + (list_Merlin[i][6][col].toarray())
#    plt.figure()    
#    plt.imshow(np.log10(summed_list_Merlin_3), cmap = 'jet', interpolation = 'none')
#    plt.colorbar()  

          
        # plot a single pil100K image
#        plt.figure()
#        plt.imshow(np.log10(diffSet_Pil100K[i][8][49]), cmap = 'hot', interpolation = 'none')
#        plt.colorbar()
#        plt.title('Scan_nbr_%d'%(first_scan_nbr+i))
        
#plotalot()

#pp=list_Merlin[0][0][0].toarray()     # OBS OBS jättestor!
#print getsizeof(pp)
#print getsizeof(list_Merlin[0][0][0])
#plt.figure()
#rowsum= np.sum(one_position,axis=1)    
#sumsum= np.sum(rowsum,axis=1)    
#xlinE = np.linspace(458,487,nbr_rotations)
#plt.plot(xlinE,sumsum,'+-')
#plt.title('Summed intensity as function of scan for one position')
#   

# TODO: gör en 3D scatter plot med färkodning för intensiteten:
    # gör meshgrids? lnispaces för xyz andvänd one_position för färg
def plot_Bragg():
#    #scatter 3 ritar ut en ring för varje punkt specificerad av vektorerna (x,y,z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #cmap = one
    ax.scatter(X, Y, Z, cmap ='jet' , c = np.log10(one_position_roi), marker ='o', alpha =0.05)
    #rätt!???:
    plt.ylabel(' rotation')
    plt.xlabel(' yDet')
    ax.set_zlabel('zDet ')
    #plt.zlabel(' zDet')
    plt.title('Bragg peak')
#plot_Bragg() 



# TODO: redo
def stage_phaseretrival_3dBragg(rotation):    
#def calc_pixelsize_Merlin():
#def calc_pixelsize_Pil100K():
    # theta should be read out from gonphi for each particular scan=?
    theta = 15*np.pi/180 #arb
    
    # Sizes of roi of diffraction patterns (pixels used on the Merlin detector)
    Ny = list_Merlin[0][0][0].shape[0]      
    Nx = list_Merlin[0][0][0].shape[1]     
    
    # define pixel sizes in reciprocal space
    dq1 = 2*np.pi*pixel_det /(wavelength * z_Merlin)
    dq2 = 2*np.pi*pixel_det /(wavelength * z_Merlin)
    
    # one pixel on the object plane correspongs to ... m in the object plane at distance z from the detector (x and y)
    # är det rätt med Nx och Ny och inte tvärtom?
    dr1 = 2*np.pi / ( Ny * dq1 * np.cos(theta))
    dr2 = 2*np.pi / ( Nx * dq2 * np.cos(theta))
    #dr3 = 2*np.pi / ( N3 * dq3 * np.cos(theta))  

    dr1 = 2*np.pi / ( Ny * dq1 * np.cos(theta))
    dr2 = 2*np.pi / ( Nx * dq2 * np.cos(theta))
    #dr3 = 2*np.pi / ( N3 * dq3 * np.cos(theta))
    
    # create realspace pixel sizes
    # hur blir det här? dr1 beror av dr3
    dx = dr3*np.cos(theta)
    dy = dr2
    #dz = dr1 + dr3*np.sin(theta)
    
    #this is the size of the object(=FOV(ja, ptycho-FOV)), right?
    # what the width of the diffraction pattern equals to in object plan (pixlar * pixelstorlekx/y)
    sizeDiffObjectx =  Nx * dx
    sizeDiffObjecty =  Ny * dy
    
 
    # calculate how long each step is in x and y OBS kan också vara minus
    stepSizex = np.zeros((nbr_rows, nbr_cols))
    stepSizey = np.zeros((nbr_rows,1))
    
    # inte tillräcklig nogrannhet i motorpos för att det ska vara nödvändigt att göra detta?
    # bara använd den stepsize du valde till scanet? testa båda sätten!
    for i in range(0,nbr_rows):   #gör 2 loops for diffrent nbr of scans in y and x . convert from microns to meters
        
        stepSizey[i] = (motorpositiony[rotation][i+1] - motorpositiony[rotation][i]) * 1E-6
        for cols in range(0,nbr_cols):
            # TODO: lots and less time to do it
            stepSizex[i][:] = (motorpositionx[rotation][i+1] - motorpositionx[rotation][i]) * 1E-6
###############################♠♠♠♠♠♠OLDOLDOLD♠♠♠     

def funkar_ej():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #X, Y, Z = axes3d.get_test_data(0.05)
    #cset = ax.contour(X, Y, one_position[1])
    #ax.contourf(X1,Y1, .1*np.sin(3*X)*np.sin(5*Y))
    
    levels = np.linspace(-1, 1, 5)
    
    ax.contourf(Y1,X1, one_position_roi[0], cmap ='hot', levels=.1*levels, zdir='z')
    ax.contourf(Y1,X1,20+ one_position_roi[8], cmap ='hot', levels=3+.1*levels, zdir='z')
    ax.contourf(Y1,X1, one_position_roi[16], cmap ='hot', levels=5+.1*levels, zdir='z')
    #ax.clabel(cset, fontsize=9, inline=1)
    
    plt.show()
