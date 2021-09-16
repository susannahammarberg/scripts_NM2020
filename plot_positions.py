import numpy as np
import matplotlib.pyplot as plt
#import nmutils
import h5py
import matplotlib.gridspec as gridspec
#import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

"""
This script visualizes the output of a ptypy run, by loading a ptyr file. (the final
recon file not the dump files

<SH> this works for python 3


print out the error! ssa?

import matplotlib as mpl
mpl.rcParams['toolbar'] = 'None'

"""
iterations = 2

folder = r'\recons'
name =  r'\orie1\orie1_DM_%04u' % iterations
# dont use this


#####

#inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons\test\test_DM_0020.ptyr'
#inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons\test\test_DM_3dBragg_0002.ptyr'
#inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons\oProbe\oProbe_DM_0003.ptyr'
#savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\plots\oProbe\oProbe_DM_0003'
#inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons\testar_save\testar_save_DM_3dBragg_0003.ptyr'

inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons' + folder + name + '.ptyr'
#savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\plots' + name
inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\recons\temp_remove\temp_remove_DM_0001.ptyr'
    

print('Opening file: \n ', inputFile)

Nx = 67
Ny = 13


### load reconstruction data
with h5py.File(inputFile, 'r') as hf:
    scanid = str(hf['content/probe'].keys()).split("['", 1)[1].split("']")[0]
    print( 'loading entry %s' % scanid)
    positions = np.array(hf.get('content/positions/%s' % scanid))


#import pdb; pdb.set_trace()
 

positions = positions[:,1:3]*1e6

# reformat to np array
#errors = np.array(errors)
##try:
##    probe = probe[0]
##    obj = obj[0]
##    psize = psize
##except IndexError:
##    raise IOError('That doesn''t look like a valid reconstruction file!')


print(positions)
print(positions.shape)
print(Nx*Ny)

x = positions[:,0].reshape(Ny,Nx)
y = positions[:,1].reshape(Ny,Nx)

fig, ax = plt.subplots(ncols=2)
ax[0].imshow(x)
plt.setp(ax[0], ylabel='y [um]', xlabel='x [um]', title='x positions in recon')
ax[1].imshow(y)
plt.setp(ax[1], ylabel='y [um]', xlabel='x [um]', title='y positions in recon')
plt.show()

fig, ax = plt.subplots()
ax.scatter(x[9,:],y[9,:], alpha=0.5)
plt.setp(ax, ylabel='y [um]', xlabel='x [um]', title='scanning positions in recon \n single row')
plt.show()

fig, ax = plt.subplots()
ax.scatter(x,y, alpha=0.5) #s = 60 (size of marker/data)
plt.setp(ax, ylabel='y [um]', xlabel='x [um]', title='scanning positions in recon \n single row')
plt.show()


##fig, ax = plt.subplots()
##ax.plot(x,y,'o')
##plt.setp(ax, ylabel='y [um]', xlabel='x [um]', title='scanning positions in recon')
##plt.show()


#import pdb; pdb.set_trace()

## means
#fig, ax = plt.subplots(ncols=3)
#ax[0].imshow(np.mean(np.angle(S_cart.data[0]), axis=1).T, extent=[x.min(), x.max(), y.min(), y.max()], interpolation='none', origin='lower', cmap='jet')
#plt.setp(ax[0], ylabel='y', xlabel='x', title='top view')
#ax[1].imshow(np.mean(np.angle(S_cart.data[0]), axis=2).T, extent=[x.min(), x.max(), z.min(), z.max()], interpolation='none', origin='lower', cmap='jet')
#plt.setp(ax[1], ylabel='z', xlabel='x', title='side view')
#ax[2].imshow(np.mean(np.angle(S_cart.data[0]), axis=0).T, extent=[ z.min(), z.max(), y.min(), y.max()], interpolation='none', origin='lower', cmap='jet')
#plt.setp(ax[2], ylabel='y', xlabel='z', title='front view')
#plt.tight_layout() 

