import numpy as np
import matplotlib.pyplot as plt
#import nmutils
import h5py
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
    
#import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
date_str = time.strftime("%Y%m%d_%H%M")
import os
from skimage.restoration import unwrap_phase  
"""
This script visualizes the output of a ptypy run, by loading a ptyr file.
"""

import matplotlib
matplotlib.use( 'Qt5agg' )

#C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/2drecons/dumps/445_20210712_1053/445_20210712_1053_EPIE_0290.ptyr
#recons/445_20210712_1053/445_20210712_1053_EPIE_0300.ptyr
iterations = 300
folder = r'\dumps'  #will only have 1 error value (the final one)
folder = r'\recons'   

#nam = '429_20220506_0734'
#nam = '430_20220505_2349'
#nam = '431_20220505_2146'
#nam = '432_20220430_1600'
#nam = '433_20220430_1346'
#nam = '434_20220202_1516' #500 iter
#nam = '434_20220430_0349'
#nam = '435_20220429_2323'
#nam = '436_20220428_2140'
#nam = '437_20220201_2126' # 500 iter
#nam = '437_20220428_1536' # test. xmotorangle = 0
#nam = '438_20210920_1139' #500 iter
#nam = '439_20220315_1120' #800 iter
#nam = '439_20220316_1452' #2000 iter
#nam = '439_20220317_1301'#_DM_0600 varying probe mkt dåligt
#nam = '439_20220330_0941'#_DM_0200
#nam = '439_20220330_1434' #DM position refinement result not so different from without pos ref
#nam = '439_20220330_1446' #DM position refinement 
#nam = '439_20220317_1532' #" DM fast then varying probe

nam = '439_20220201_1556' #500 iter THIS
#nam = '439_20220204_1423' # varying probe 

#nam = '440_20210920_1538' #200 iter
#nam = '440_20220501_0007'    # 800 iterationer (dumps)
#nam = '441_20210920_1635'

#nam = '442_20210921_1005'
#nam = '443_20210921_1124'


#nam =  '444_20210921_1339'
#nam = '445_20210712_1053' shape 170
#nam = '445_20210921_1525'
#nam = '446_20210921_2019' 
#nam = '447_20210921_2228' #300
nam = '448_20210922_0948' 
#nam = '449_20210922_1233'
#nam = '450_20210922_1548'
#nam = '451_20210922_2008'
#nam = '452_20210922_2220'
#nam = '453_20210923_1026'
#nam = '453_20220401_2135'
#GaInP                   
#nam =  '464_20220524_0946'     #3000recon
#nam = '465_20220516_2246' #GainP   5eller800    #intressant. kanske rekonstruerar segment 4
#nam = '468_20220525_2108' #5000
#nam = '470_20220516_1843' #600iter ofullständig
#nam = '470_20220525_1101' #cluster 3500
#nam = '472_20220517_0942'#_EPIE_2000
#nam = '474_20220523_1329' #2000 (cluster)



# gainp
#nam = '453_20220401_0926'


scan = int(nam[0:3])
#OBSOBS
name =  "\\" + nam + "\\" + nam + '_EPIE_%04u' % iterations
#name =  "\\" + nam + "\\" + nam + '_DM_%04u' % iterations

outputSuffix = 'png'
save = False

plt.close('all')


#inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\2drecons' + folder + name + '.ptyr'

inputFile = r'C:\Users\sanna\Phd\NanoMAX_May2020\Analysis\scans429_503\2drecons' + folder + name + '.ptyr'
savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\2drecons\plots' + name

if save == True:
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        print('new folder in this savepath was created')
    

print('Opening file: \n ', inputFile)

### load reconstruction data
with h5py.File(inputFile, 'r') as hf:
    scanid = str(list(hf['content/probe'].keys())[0])
    print('loading entry %s' % scanid)
    probe = np.array(hf.get('content/probe/%s/data' % scanid))
    obj = np.array(hf.get('content/obj/%s/data' % scanid))
    psize = np.array(hf.get('content/probe/%s/_psize' % scanid))
    energy = np.array(hf.get('content/probe/%s/_energy' % scanid))
    origin = np.array(hf.get('content/probe/%s/_origin' % scanid))

    err1 = []
    err2 = []
    for ite in range(0,iterations):
        hh = np.array(hf.get('content/runtime/iter_info/%s/error'%('{0:05}'.format(ite))))        
        if hh.size == 3:
            err1.append(hh[0])
            err2.append(hh[2])
            #print('saved error %d'%ite)
try:
    probe = probe[0]
    obj = obj[0]
    psize = psize[0]
except IndexError:
    raise IOError('That doesn''t look like a valid reconstruction file!')
print("Loaded probe %d x %d and object %d x %d, pixel size %.1f nm, energy %.2f keV"%(probe.shape + obj.shape + (psize*1e9, energy)))



def prop_probe():
    backProp = -5000.0
    forwProp = 5000.0
    steps = 4 #400
    
    probe_basez_4000 = nmutils.utils.propagateNearfield(probe, psize, -4000*1e-6, energy)
    plt.figure()
    plt.imshow(abs(probe_basez_4000[0]),cmap='RdBu_r')
    plt.axis('off')
    
    ### define distances and propagate
    dist = np.linspace(backProp, forwProp, steps) * 1e-6
    dx = dist[1] - dist[0]
    print("propagating to %d positions separated by %.1f um..."\
        % (len(dist), dx*1e6))
    ### not sure why, but the propagation goes in the other direction here!
    ### it could be a misunderstanding about motors at nanomax...
    field3d = nmutils.utils.propagateNearfield(probe, psize, -dist, energy)
    
    ### get intensities and focii
    power3d = np.abs(field3d)**2
    power_vertical = np.sum(power3d, axis=2).T
    power_horizontal = np.sum(power3d, axis=1).T
    focus_vertical_ind = np.argmax(np.sum(power_vertical**2, axis=0))
    focus_vertical_x = dist[focus_vertical_ind]
    focus_horizontal_ind = np.argmax(np.sum(power_horizontal**2, axis=0))
    focus_horizontal_x = dist[focus_horizontal_ind]
    focus_ind = np.argmax(np.sum(power3d**2, axis=(1,2)))
    focus_x = dist[focus_ind]
    focus = field3d[focus_ind]
    
    ### plot
    fig = plt.figure(figsize=(8, 10))
    outer_grid = gridspec.GridSpec(2, 2, wspace=.2, hspace=.2, height_ratios=[2,3])
    
    # probe and focus spots
    def spot_subplot(gridcell, data, shareax=None, title=''):
        subgrid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gridcell, width_ratios=[3,1], height_ratios=[1,3], hspace=.05, wspace=.05)
        lims = [-1e6*data.shape[0] * psize / 2, 1e6*data.shape[0] * psize / 2] # um
        posrange = np.linspace(lims[0], lims[1], data.shape[0])
        ax = plt.subplot(subgrid[1,0], sharex=shareax, sharey=shareax)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
        ax.imshow(nmutils.utils.complex2image(data), extent=lims+[lims[1], lims[0]], interpolation='none')
        ax_histh = plt.subplot(subgrid[0,0], sharex=ax, yticklabels=[], ylabel='Int.')
        ax_histv = plt.subplot(subgrid[1,1], sharey=ax, xticklabels=[], xlabel='Int.')
        ax_histv.plot(np.sum(np.abs(data)**2, axis=1), posrange)
        ax_histh.plot(posrange, np.sum(np.abs(data)**2, axis=0))
        ax_histh.set_title(title, x=.67)
        for tk in ax_histh.get_xticklabels(): tk.set_visible(False)
        for tk in ax_histv.get_yticklabels(): tk.set_visible(False)
    
        # FWHM:
        import scipy.interpolate
        for i in (0,1):
            y = np.sum(np.abs(data)**2, axis=i)
            edges = scipy.interpolate.UnivariateSpline(posrange, y-y.max()/2).roots()
            r1, r2 = edges[0], edges[-1]
            if i == 0:
                ax_histh.axvspan(r1, r2, fc='r', alpha=.3)
                ax_histh.text(r2, np.mean(ax_histh.get_ylim()), ' %.0f nm'%((r2-r1)*1e3), fontsize=10, rotation=-90*i)
                ax.set_xlim(np.mean([r1,r2]) + np.array([-1,1]))
            elif i == 1:
                ax_histv.axhspan(r1, r2, fc='r', alpha=.3)
                ax_histv.text(np.mean(ax_histv.get_xlim()), r1, ' %.0f nm'%((r2-r1)*1e3), fontsize=10, va='top', ha='center', rotation=-90)
                ax.set_ylim(np.mean([r1,r2]) + np.array([-1,1]))
        return ax
    a = spot_subplot(outer_grid[0,0], probe, title='Sample plane')
    a.set_ylabel('$\mu$m')
    
    b = spot_subplot(outer_grid[0,1], focus, shareax=a, title='Focal plane')
    
    # beam profiles
    subgrid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_grid[1,:], hspace=.05)
    ax_vertical = plt.subplot(subgrid[0])
    ax_horizontal = plt.subplot(subgrid[1], sharex=ax_vertical, sharey=ax_vertical)
    ax_vertical.imshow(power_vertical, cmap='gray', extent=[1e6*dist[0], 1e6*dist[-1], -1e6*psize*probe.shape[0]/2, 1e6*psize*probe.shape[1]/2], interpolation='none')
    ax_vertical.axvline(x=focus_vertical_x*1e6, lw=1, ls='--', color='r')
    ax_vertical.text(1e6*focus_vertical_x, ax_vertical.get_ylim()[0] + .2*np.diff(ax_vertical.get_ylim())[0], 
        ' %.0f um '%(1e6*focus_vertical_x), color='red', 
        ha=('right' if focus_vertical_x<focus_x else 'left'))
    ax_vertical.axvline(x=focus_x*1e6, lw=2, ls='-', color='r')
    ax_vertical.text(1e6*focus_x, ax_vertical.get_ylim()[0] + .8*np.diff(ax_vertical.get_ylim())[0],
        ' %.0f um '%(1e6*focus_x), color='red', va='top',
        ha=('right' if focus_vertical_x>focus_x else 'left'))
    ax_vertical.set_aspect('auto')
    ax_horizontal.imshow(power_horizontal, cmap='gray', extent=[1e6*dist[0], 1e6*dist[-1], -1e6*psize*probe.shape[0]/2, 1e6*psize*probe.shape[1]/2], interpolation='none')
    ax_horizontal.axvline(x=focus_horizontal_x*1e6, lw=1, ls='--', color='r')
    ax_horizontal.text(1e6*focus_horizontal_x, ax_horizontal.get_ylim()[0] + .2*np.diff(ax_horizontal.get_ylim())[0],
        ' %.0f um '%(1e6*focus_horizontal_x), color='red',
        ha=('right' if focus_horizontal_x<focus_x else 'left'))
    ax_horizontal.axvline(x=focus_x*1e6, lw=2, ls='-', color='r')
    ax_horizontal.set_aspect('auto')
    ax_horizontal.set_ylabel('$\mu$m', y=1.05)
    for tk in ax_vertical.get_xticklabels(): tk.set_visible(False)
    ax_horizontal.set_xlabel('beamline z axis ($\mu$m)', fontsize=16)
    
    lblue = (.3,.3,1.0)
    ax_vertical.axvline(x=0, lw=1, ls='--', color=lblue)
    ax_horizontal.axvline(x=0, lw=1, ls='--', color=lblue)
    ax_horizontal.text(0, ax_horizontal.get_ylim()[0] + .01*np.diff(ax_horizontal.get_ylim())[0], 
        ' sample ', color=lblue, va='bottom',
        ha=('left' if focus_x < 0 else 'right'))
    
    ax_horizontal.text(1.05, 0.5, 'horizontal focus (M2)',
         horizontalalignment='center',
         verticalalignment='center',
         rotation=90,
         transform = ax_horizontal.transAxes)
    ax_vertical.text(1.05, 0.5, 'vertical focus (M1)',
         horizontalalignment='center',
         verticalalignment='center',
         rotation=90,
         transform = ax_vertical.transAxes)
    
    if save is True:
        fn = savepath + '\\' + 'probe' + '.' + outputSuffix
        plt.savefig(fn)
        print('Saved to %s'%fn)
    
    plt.show()
#prop_probe()
#extent = 1e6 * np.array([origin[0], origin[0]+(obj.shape[1]-1)*psize, origin[1], origin[1]+(obj.shape[0]-1)*psize])'
extent = 1e6 * np.array([0, (obj.shape[1])*psize, 0, (obj.shape[0])*psize])
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""obsbosbos""""""""""""2"""
obj=np.fliplr(obj)
signflip = False

### Object
#%%

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def plot2drecons(obj, probe, extent, savepath, save, signflip):
    
    plt.close('all')
    
    fig, ax = plt.subplots(ncols=2, figsize=(10,3), sharex=True, sharey=True)
    #fig, ax = plt.subplots(ncols=2, figsize=(10,6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=.3)
    #fig.suptitle(title, fontsize=15)
       
    # amplitude
    mag = np.abs(obj)
    #slicey = slice(90,130)  #for real data recons with shape 170
    #slicex = slice(100,290) #for real data recons with shape 170
    global slicey
    global slicex
    

    slicey = slice(149-7,183-7)  #for real data recons with shape 256
    slicex = slice(195,379) #for real data recons with shape 256   (#S434) xrange 184 y range 34 # NEW ROI
    
    slicey = slice(149-4,183-4)  #for real data recons with shape 256
    slicex = slice(195,379) #for real data recons with shape 256   (#S435) xrange 184 y range 34 # NEW ROI
    
    slicey = slice(149,183)  #for real data recons with shape 256
    slicex = slice(195,379) #for real data recons with shape 256   (#S437) xrange 184 y range 34 # NEW ROI
    
    #slicey = slice(149+5,183+5)  #for real data recons with shape 256
    #slicex = slice(195,379) #for real data recons with shape 256   (#S438) xrange 184 y range 34 # NEW ROI
    
    slicey = slice(158,192)  #for real data recons with shape 256
    slicex = slice(195,379) #for real data recons with shape 256   (#S439 iaf) NY ROIT NEW ROII

    #slicey = slice(152+5,186+5)  #for real data recons with shape 256
    #slicex = slice(198,382)    #440 
    
    #Slice for plotting lineout of 170nm segmnet in 440
    #slicey = slice(152+5,186+5)  #for real data recons with shape 256
    #slicex = slice(198+30,382+30)    #440 
    
    
    #slicey = slice(152,186)  #for real data recons with shape 256
    #slicex = slice(198-4,382-4)    #441 
    
    #slicey = slice(164-12,198-12)  #for real data recons with shape 256
    #slicex = slice(210,394)    #448       184 i x 34 i y 
    
    #slicey = slice(164-8,198-8)  #for real data recons with shape 256
    #slicex = slice(210+5,394+5)    #450       184 i x 34 i y #

    #slicey = slice(164-3,198-3)  #for real data recons with shape 256
    #slicex = slice(210+8,394+8)    #451       184 i x 34 i y #
    
    
    #slicey = slice(164,198)  #for real data recons with shape 256
    #slicex = slice(210+8,394+8)    #452       184 i x 34 i y     
    
    ###slicey = slice(158,192)  #for real data recons with shape 256 v2 to match xrd
    ###slicex = slice(198,386) #for real data recons with shape 256  v2 to match xrd

    slicey = slice(132,166)     # for simulated data recons 256 (should be the same, but for now)
    slicex = slice(201+2,385+2)     # NEW ROI    #sim 33 # and 31
   
    #slicey = slice(132,166)     # for simulated data recons 256 (should be the same, but for now)
    #slicex = slice(201+2,385+2)     # NEW ROI    #sim 35
    
    #slicey = slice(132,166)     # for simulated data recons 256 (should be the same, but for now)
    #slicex = slice(201,385)     # NEW ROI    #sim 39
    
    #slicey = slice(132,166)     # for simulated data recons 256 (should be the same, but for now)
    #slicex = slice(201+3,385+3)     # NEW ROI    #sim 41

    #slicey = slice(132,166)     # for simulated data recons 256 (should be the same, but for now)
    #slicex = slice(201+7,385+7)     # NEW ROI    #sim 44
    
    #slicey = slice(164-5,198-5)  #for real data recons with shape 256
    #slicex = slice(210-10,394-10)    #452       184 i x 34 i y 

    #slicey = slice(164-10,198-10)  #for real data recons with shape 256
    #slicex = slice(210-40,394-10)    #GaINP 470

    #slicey = slice(164-10-5,198-10-5)  #for real data recons with shape 256
    #slicex = slice(210-40,394-10)    #GaINP 470

    #mag_cut = mag[mag.shape[0]//3:2*mag.shape[0]//3, mag.shape[1]//3:2*mag.shape[1]//3] # to find relevant dynamic range
    
    
    # Normalize magnitude
    mag = mag
    mag_cut = mag[slicey,slicex]/mag[slicey,slicex].max() # to find relevant dynamic range
    
    vmin = mag_cut.min()
    vmax = mag_cut.max()
    print('max ', vmax)
    print('vmin', vmin)
    
    #abs
    img = ax[0].imshow(mag, cmap='gray', extent=extent, interpolation='none')
    
    plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=70)
    ax[0].set_ylabel('$\mu$m')
    ax[0].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
    ax[0].set_title('Amplitude')
    
    
    # phase
    img = ax[1].imshow(np.angle(obj), cmap='jet', extent=extent, interpolation='none')
    #img = ax[1].imshow((np.angle(obj)), vmin=-2.5, vmax=-2, cmap='hsv', extent=extent, interpolation='none')
    plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=70 )
    ax[1].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax[1].set_title('Phase')
    
    plt.show()
    #import pdb; pdb.set_trace()
    if save is True:
        fn = savepath + '\\' + 'object' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    #zoomed in version
    global extent_zoomed
    #extent_zoomed = 1e6 * np.array([origin[0], origin[0]+(mag_cut.shape[1]-1)*psize, origin[1], origin[1]+(mag_cut.shape[0]-1)*psize])
    extent_zoomed = 1e6 * np.array([0, (mag_cut.shape[1])*psize, 0, (mag_cut.shape[0])*psize])
    
    print('*extent_zoomed*******************************',extent_zoomed)
    
    #mask with amplitude  
    amplitude_range = vmax-vmin
    #mask_at = 0.053 * (vmin + amplitude_range)
    mask_at = 0.4 * (vmin + amplitude_range)
    #mask_at = 8E-6 #simulated'
    #mask_at = 1.0 #simulated
    mask_at = 0.08 #434
    mask_at = 0.07 #435
    #mask_at = 0.025 #436
    mask_at = 0.04 #437
    #mask_at = 0.05 #439 2000
    mask_at = 0.043 #439
    #mask_at = 0.1
    #mask_at = 0.07 #441 
    #mask_at = 0.07 #442 
    #mask_at = 0.1#448
    #mask_at = 0.1#451
    #mask_at = 0.15    #sim 29
    mask_at = 0.08   #sim 31 THIS
    
    #mask_at = 0.1  #sim 33
    mask_at = 0.1  #sim 35
    
    
    #mask_at = 0.15    #sim 44
    
    #mask_at = 0.15    #470 gainp
    
    mask = np.zeros((mag_cut.shape))
    mask[(np.abs(mag_cut)) > mask_at] = 1
    #mask[np.log10(np.abs(mag_cut)) > mask_at] = 1
    mask = np.squeeze(mask)
    
    
    fig, ax = plt.subplots(ncols=2, figsize=(10,3), sharex=True, sharey=True)
    #fig, ax = plt.subplots(ncols=2, figsize=(10,6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=.3)
    #vmin = mag_cut.min()
    #vmax = mag_cut.max()
    #vmin = 0
    #vmax = 1 
    #abs
    img = ax[0].imshow(mag_cut,cmap='gray', extent=extent_zoomed, interpolation='none')#, vmin=vmin, vmax=vmax)#
    plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=70)
    ax[0].set_ylabel('$\mu$m')
    ax[0].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
    ax[0].set_title('Amplitude')
    
    # phase
    img = ax[1].imshow(mask*  np.angle(obj[slicey,slicex]), cmap='jet', extent=extent_zoomed, interpolation='none') #vmin=-np.pi, vmax=np.pi)#,
    #img = ax[1].imshow((np.angle(obj)), vmin=-2.5, vmax=-2, cmap='hsv', extent=extent, interpolation='none')
    plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=70 )
    ax[1].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax[1].set_title('Phase')
    #plt.savefig('hyyhy6')
   # plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'object_zoomed' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)

    #slicex_170 = slice(50,190)
    fig, ax = plt.subplots()
    ##amp_sqrt = (np.sqrt(mag_cut))
    #img = ax.imshow(amp_sqrt/amp_sqrt.max(), extent=extent_zoomed, interpolation='none', vmin=0, vmax=1)#
    print(' not taking sqrt*******************************************')
    img = ax.imshow(mag_cut, extent=extent_zoomed, interpolation='none', vmin=0, vmax=1)#
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(0,1))
    #plt.colorbar(img, orientation='horizontal')
    ax.set_title('Amplitude')  
    
    if save is True:
        fn = savepath + '\\' + 'temp_amp_bject_zoomed' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    
    masked_phase = mask*  np.angle(obj[slicey,slicex])
    masked_phase[masked_phase==0] = np.nan
    
    fig, ax = plt.subplots() 
    #70 segment zoomed
    img = ax.imshow( masked_phase, cmap='jet', interpolation='none')#, extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
#    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ##cb = plt.colorbar(img, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi), orientation='horizontal')
    ##cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))    
    cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude')
    plt.show()
          
    if save is True:
        fn = savepath + '\\' + 'phase_object_zoomed' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        
        
    fig, ax = plt.subplots() 
    #70 segment zoomed
    img = ax.imshow( -1* masked_phase, cmap='jet', interpolation='none')#, extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
#    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ##cb = plt.colorbar(img, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi), orientation='horizontal')
    ##cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))    
    cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude *-1')
    plt.show()
    
    #phase ramp
    phase_ramp = np.zeros((obj[slicey,slicex].shape))
    #phase_ramp[:,151:179] = np.concatenate((np.linspace(0,np.pi,num=14), np.linspace(-np.pi,0,num=14)))
    #phase_ramp[:,151:179] = np.concatenate((np.linspace(1.38,np.pi,num=14), np.linspace(-np.pi,-1.48,num=14)))
    #jättebraaaa. sen var bara ännu mer noggrann så tror jag det blir rätt
    #phase_ramp[:,153:178+1] = np.concatenate((np.linspace(1.3776,np.pi,num=12), np.linspace(-np.pi,-1.477,num=14)))
    #phase_ramp[:,153:178+1] = np.concatenate((np.linspace(0.5,np.pi,num=12), np.linspace(-np.pi,-0.5,num=14))) #bäst
    #phase_ramp[:,153:178+1] = np.concatenate((np.linspace(0.0,np.pi,num=12), np.linspace(-np.pi,-0.0,num=14)))
    # mkt bra. unwrapped phase ser ut som simulerat displacemnt field utan fel!
    phase_ramp[:,151:180+1] = np.concatenate((np.linspace(0.0,np.pi,num=14), np.linspace(-np.pi,-0.0,num=16)))
    # kan också testa modifuer pi lite i linspace
    
    ### math: np.angle(np.exp(1j*x)) = x
    #phase_ramp = np.exp(-1j*number*phase_ramp)
    #phase_ramp_conj = np.exp(1j*number*phase_ramp)
    
    #attempt to remove the phase ramp! (centered in the recostructed phase)
    fig, ax = plt.subplots() 
    img = ax.imshow(phase_ramp, cmap='jet')
    #    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ##cb = plt.colorbar(img, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi), orientation='horizontal')
    ##cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))    
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('RAMP')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'phase_RAMP' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    # multiply with complex conjugate of ramp to remove
    fig, ax = plt.subplots() 
    #img = ax.imshow( (1.0/np.exp(1j*InP_Qvect*55E-6)))*(mask* ( np.angle(obj[slicey,slicex]))), cmap='jet', interpolation='none', extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
    img = ax.imshow( mask*  np.angle(np.exp(-1j*phase_ramp)*obj[slicey,slicex]), cmap='jet')#, vmin=-4,vmax=4)  
    #    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    ##cb = plt.colorbar(img, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi), orientation='horizontal')
    ##cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))    
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude. REMOVED RAMP')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'phase_object_unramped' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    unramped_unwrapped_phase = unwrap_phase(mask*  np.angle(np.exp(-1j*phase_ramp)*obj[slicey,slicex]))
    unramped_unwrapped_phase[unramped_unwrapped_phase==0] = np.nan
    
    fig, ax = plt.subplots()
    img = ax.imshow( unramped_unwrapped_phase, cmap='jet',vmin =-4.644684461087706, vmax =9.534697854977901, interpolation='none', extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude UNRAMPED  unwrapped Vmin vmax')
    #plt.savefig('hyyhy6')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'phase_object_unramped_unwrapped' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        
    
    
    #simulated
    unwrapped_phase_sim = unwrap_phase(-1* mask* ( np.angle(obj[slicey,slicex])))
    print('*********************************************************************')
    print(' Min of this: ', unwrapped_phase_sim.min())
    print(' Max of this: ', unwrapped_phase_sim.max())
    unwrapped_phase_sim[unwrapped_phase_sim==0] = np.nan
    
     #Min of this:  -4.644684461087706
     #max of this:  9.534697854977901

    fig, ax = plt.subplots()
    #70 segment zoomed
    img = ax.imshow(unwrapped_phase_sim, cmap='jet', interpolation='none', extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Simulated phase masked with log10 amp unwrapped (without unwrapmed)')
    #plt.savefig('hyyhy6')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'phase_object_zoomed_unwrapped' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        print('shape', unwrap_phase(mask* ( np.angle(obj[slicey,slicex]))).shape)
        

    
    #calc strain and relative strain
    #--------------------------------------------------------------------    
    
    # for the simulated data I should deffenitely use the theoretical q-vect that corresponds to theta 10.54  which is the theoretical value that I have used all the way..
    # for the experimentall i am not sure. I have calibrated theta so I could i use that value. does not matter much for the strain maps acc to jesper. check
    InP_Qvect = 18543793660.27452 #this is used in simulation
    #InP_Qvect = 29543793660 #this is used in simulation
    # ska man använda den Q-vektor som sätts at theta, alltså det theta om mäts upp ~~8 deg istället för detta som motsvarar 10.54
    dz = psize
    
    #masked_phase = mask *  np.angle(np.squeeze(obj))
    # disaplacement in z, du/dz är det jag vill plotta
    # phase phi = u*Q
    #calculate the gradient before you mask other wise you
    # will get hard gradients at the mask edges naturally
    
    print('hhjjhheeeeeeeeeeeeeeeeeeeeloooooooooooooooooooooooooooooooooooooooooooooooooooooooo')
    fig, ax = plt.subplots()
    #70 segment zoomed
    #img = ax.imshow( (np.cos(np.deg2rad(11))*1.0/InP_Qvect) * unwrap_phase((mask* np.angle(obj[slicey,slicex]))), cmap='jet', interpolation='none', extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
    # not sure what the cos factor does (np.cos(np.deg2rad(11))*1.0/InP_Qvect)
    img = ax.imshow( 1.0/ InP_Qvect *  unwrap_phase(mask*  np.angle(np.exp(-1j*phase_ramp)*obj[slicey,slicex])),  cmap='jet', interpolation='none', extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Displacement field from unramped and unwrapped phase masked with log10 amplitude. Scaled to show seg 5 colors')
    #plt.savefig('hyyhy6')
    plt.show()

   
    #strain = np.diff(np.angle(obj[slicey, slicex]) /InP_Qvect , axis = 1, append=0) /dz
    plt.figure()
    plt.title('what i calculate gradient from: unwrapped phase (wthout using mask)')
    plt.imshow(unwrap_phase(np.angle(obj[slicey, slicex])), cmap='jet')  
    
    #experimetal
    strain = np.gradient(((unwrap_phase(np.angle(obj[slicey, slicex])))/InP_Qvect) , dz)[1]
    #simulation
    #strain = np.gradient(((unwrap_phase(-1*np.angle(obj[slicey, slicex])))/InP_Qvect) , dz)[1]
    
    
    #experimetal
    strain_unramped = np.gradient(((unwrap_phase(np.angle(np.exp(-1j*phase_ramp)*obj[slicey, slicex])))/InP_Qvect) , dz)[1]
    #simulation
    #strain_unramped = np.gradient(((unwrap_phase(-1*np.angle(np.exp(-1j*phase_ramp)*obj[slicey, slicex])))/InP_Qvect) , dz)[1]
    
    fig, ax = plt.subplots(nrows =3)
    normT = TwoSlopeNorm(vmin=-0.45, vcenter=0, vmax=1.46)
    ax[0].imshow(np.angle(obj[slicey, slicex]),cmap='jet')
    ax[1].imshow(1.0/ InP_Qvect *  unwrap_phase(mask*  np.angle(np.exp(-1j*phase_ramp)*obj[slicey,slicex])),vmin=-3E-10,vmax = 1E-10,cmap='jet')
    ax[2].imshow(mask* strain*100,cmap='RdBu_r', norm=normT)
    
    
    if signflip == True:
        strain *= -1
        
    global masked_strain
    masked_strain = mask* strain
    
    ##Average over the reigion of interest (that is a bit arbitrary)
    #mean_strain = np.mean(masked_strain[np.nonzero(masked_strain)])
    #print('mean value in the masked strain',mean_strain)
    
#    print('""""""""""""""""""""TEMOP TEMP remove some maxmimum pixels""""""""""""""""""""')
#    masked_strain[abs(masked_strain)==abs(masked_strain).max()] = 0
#    masked_strain[abs(masked_strain)==abs(masked_strain).max()] = 0
#    masked_strain[abs(masked_strain)==abs(masked_strain).max()] = 0
#    masked_strain[abs(masked_strain)==abs(masked_strain).max()] = 0
#    masked_strain[abs(masked_strain)==abs(masked_strain).max()] = 0
#    

    
#    remove hot or cold pixels
    #maybe it is pretty ok to use this . but i am not saying which pixels are white and which are undefined
    #for iii in range(30):
    #    masked_strain[abs(masked_strain)==abs(masked_strain).max()] = 0

    
    lineout_row = int(masked_strain.shape[0]/2)
    
    # plot strain
    fig, ax = plt.subplots(ncols=1)
    # to get to colormap at 0 at white color 
    #print('OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOBS NOT THE SAME NORM')
    #print('strain without norm')
    norm = TwoSlopeNorm(vmin=-0.45, vcenter=0, vmax=1.46) #the one used in paper
    #norm = TwoSlopeNorm(vmin=-0.1, vcenter=0, vmax=0.3) #470
    #norm = TwoSlopeNorm(vmin=-0.1, vcenter=0, vmax=0.4)
    img = ax.imshow(masked_strain*100, cmap='RdBu_r', extent=extent_zoomed, interpolation='none',  norm=norm)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    #ax.axhline(y=0.1444,color='red')
    ax.set_xlabel('$\mu$m')
    ax.set_ylabel('$\mu$m')
    divider = make_axes_locatable(ax)   
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax, ticks=(-0.4, 0, 0.4, 0.8, 1.2))
    ##cb = plt.colorbar(img, ticks=(-0.4, 0, 0.4, 0.8, 1.2), orientation='horizontal')
    #cb = plt.colorbar(img, cax=cax, ticks=((100*masked_strain).min(), 0, (100*masked_strain).max()))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Strain [%]')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'strain' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        
    #plot unramped strain:
    fig, ax = plt.subplots(ncols=1)
    # to get to colormap at 0 at white color 
    #print('OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOBS NOT THE SAME NORM')
    #print('strain without norm')
    norm = TwoSlopeNorm(vmin=-0.45, vcenter=0, vmax=1.46) #the one used in paper
    #norm = TwoSlopeNorm(vmin=-0.1, vcenter=0, vmax=0.3) #470
    #norm = TwoSlopeNorm(vmin=-0.1, vcenter=0, vmax=0.4)
    img = ax.imshow(mask*strain_unramped*100, cmap='RdBu_r', extent=extent_zoomed, interpolation='none',  norm=norm)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    #ax.axhline(y=0.1444,color='red')
    ax.set_xlabel('$\mu$m')
    ax.set_ylabel('$\mu$m')
    divider = make_axes_locatable(ax)   
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax, ticks=(-0.4, 0, 0.4, 0.8, 1.2))
    ##cb = plt.colorbar(img, ticks=(-0.4, 0, 0.4, 0.8, 1.2), orientation='horizontal')
    #cb = plt.colorbar(img, cax=cax, ticks=((100*masked_strain).min(), 0, (100*masked_strain).max()))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('UnRAMPed Strain [%]')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'unramped_strain' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    
    fig, ax = plt.subplots()
    ax.imshow(mask* np.gradient((((np.angle(obj[slicey, slicex])))/InP_Qvect) , dz)[1])
    ax.set_title('strain calc from phase without unwrapping')
    
        #import pdb; pdb.set_trace()
    print('binned strain (se what happens if lower res)')
#
    strain_bin = rebin(masked_strain[:32],(8,46))
    norm = TwoSlopeNorm(vcenter=0)
#    
    fig,ax = plt.subplots()   
    img_xx = ax.imshow(strain_bin*100, cmap='RdBu_r' , norm = norm)
    divider = make_axes_locatable(ax)  
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img_xx, cax=cax)

    lineout_xx= np.linspace(0,extent_zoomed[1],masked_strain.shape[1])
    
    lineout = (masked_strain*100)[lineout_row]
    
    plt.figure()
    plt.plot(lineout_xx, mag_cut[lineout_row])
    
    # plot lineoput
    plt.figure()
    plt.plot(lineout_xx,lineout,'-')

    if save is True:
        fn = savepath + '\\' + 'strain_lineout'
        np.save(fn,lineout)
        plt.savefig(fn)
        fn2 = savepath + '\\' + 'strain_lineout_xx'
        np.save(fn2,lineout_xx)
        
        print("Saved to %s"%fn) 
        
    fig, ax = plt.subplots(ncols=1)
    # to get to colormap at 0 at white color 
    norm = TwoSlopeNorm(vmin=-0.45, vcenter=0, vmax=1.46)
    img = ax.imshow(masked_strain*100, cmap='RdBu_r',interpolation='none',norm=norm)### extent=extent_zoomed, interpolation='none',  norm=norm)
    ax.axhline(y=lineout_row,color='red')
    
    
    if save is True:
        fn = savepath + '\\' + 'strain_lineout_map'
        plt.savefig(fn)
        print("Saved to %s"%fn) 
        
        
    #---------------------------------------------
    # plot zoomed version (zoomed on 2 largest segments)     
    #slicey2 = slice(8,50)
    #slicex2 = slice(0,100)
    
    slicey2 = slice(8,-9) #439
    slicex2 = slice(154,179)

#    slicey2 = slice(5,-5) #sim 33
#    slicex2 = slice(156,177)
    #slicey2 = slice(8,30)    #sim data
    #slicex2 = slice(12,66)    #sim data

    #slicey2 = slice(0,-1)    #sim data  245
    #slicex2 = slice(153,178+12)    #sim data 256  
    
    #for 170nm segment simulation
    #slicey2 = slice(50,150)    #sim data 245
    #slicex2 = slice(50,190)    #sim data 256
    
    #slicey2 = slice(0,-1)    #for 170nm segment simulation
    #slicex2 = slice(130,-20)    #for 170nm segment simulation

    
    masked_strain_zoomed = masked_strain[slicey2,slicex2]
    extent_zoomed2 = 1e6 * np.array([origin[0], origin[0]+(masked_strain_zoomed.shape[1]-1)*psize, origin[1], origin[1]+(masked_strain_zoomed.shape[0]-1)*psize])
   
    print('average value of strain in segment 5: ')

    pix = masked_strain_zoomed[np.nonzero(masked_strain_zoomed)] 
    print(np.average(pix*100))
    print(np.sum(masked_strain_zoomed*100)/(masked_strain_zoomed.shape[0]*masked_strain_zoomed.shape[1]))
    
    print('same average strain again')
    print(masked_strain[np.nonzero(masked_strain)].shape)
    #chose the non zero values of the cut out of the strain map. take the average
    print(np.average(100* masked_strain[slicey2,slicex2][np.nonzero(masked_strain[slicey2,slicex2])] ))
    
    print('segment 4')
    print(np.average(100* masked_strain[slicey2,slice(111,128)][np.nonzero(masked_strain[slicey2,slice(111,128)])] ))

    print('segment 3')
    print(np.average(100* masked_strain[slicey2,slice(65,100)][np.nonzero(masked_strain[slicey2,slice(65,100)])] ))

    print('segment 2')
    print(np.average(100* masked_strain[slicey2,slice(25,65)][np.nonzero(masked_strain[slicey2,slice(25,65)])] ))   

# remove hot or cold pixels
#    masked_strain_zoomed[abs(masked_strain_zoomed)==abs(masked_strain_zoomed).max()] = np.nan
#    masked_strain_zoomed[abs(masked_strain_zoomed)==abs(masked_strain_zoomed).max()] = 0
#    masked_strain_zoomed[abs(masked_strain_zoomed)==abs(masked_strain_zoomed).max()] = 0
#    masked_strain_zoomed[abs(masked_strain_zoomed)==abs(masked_strain_zoomed).max()] = 0
    
    fig, ax = plt.subplots(ncols=1)
    # to get to colormap at 0 at white color 
    norm = TwoSlopeNorm( vcenter=0)
    #img = ax.imshow(masked_strain_zoomed*100, cmap='RdBu_r', extent=extent_zoomed2, interpolation='none', vmin= -1.5, vmax=2.2)#exp#,  norm=norm)#, vmin=vs_min, vmax=vs_max)
    img = ax.imshow(masked_strain_zoomed*100, cmap='RdBu_r', extent=extent_zoomed2, interpolation='none', vmin= -1.5, vmax=2.2)#exp
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_xlabel('$\mu$m')
    ax.set_ylabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Strain [%] \n automatic colorbar range')
    #plt.savefig('strain')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'zoomed_strain' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        
    #temp displacement based on unwrapped phase from 170 semgent only
    fig, ax = plt.subplots(ncols=1)
    # to get to colormap at 0 at white color 
    #norm = DivergingNorm( vcenter=0)
    #img = ax.imshow(masked_strain_zoomed*100, cmap='RdBu_r', extent=extent_zoomed2, interpolation='none', vmin= -1.5, vmax=2.2)#exp#,  norm=norm)#, vmin=vs_min, vmax=vs_max)
     #img = ax.imshow( (np.cos(np.deg2rad(11))*1.0/InP_Qvect) * unwrap_phase( ( np.angle(obj[slicey,slicex])))
    img = ax.imshow((np.cos(np.deg2rad(11))*1.0/InP_Qvect) * unwrap_phase(np.angle(obj[slicey, slicex][slicey2,slicex2])), cmap='jet', extent=extent_zoomed2, interpolation='none')#exp
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_xlabel('$\mu$m')
    ax.set_ylabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('displacement field from cropped unwrapped strain')
    #plt.savefig('strain')
    plt.show()


"""""""""""""""""""""""""""""""""""""""""OBS run plotcode""""""""""" 


plot2drecons(obj, probe, extent,savepath, save, signflip)

#%%
def plot_recon_error(error, save):
    
    plt.figure()
    plt.title('Reconstruction error')
    plt.plot(error,'blue')
    plt.xlabel('iteration /10')
    #plt.plot(abs(ferrors),'red')
    if save is True:
        fn = savepath + '\\' + 'error' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
#plot_recon_error(err1, save)
#plot_recon_error(err2, save)
#%%
xx = np.linspace(extent_zoomed[0],extent_zoomed[1],obj[slicey,slicex].shape[1])



#plt.close('all')

def plot_line_plot(obj,xx,row,save):
    
    plt.figure()
    #Line in 2d,  mask the phase with amplitude
    #plt.title('Lineplot of phase wrapped and unwrapped')
    #plt.title('Line plot of strain [%]')


    plt.plot(xx, np.angle(obj[row]),'.-')
    plt.plot(xx, np.unwrap(np.angle(obj[row])),'r.-')
    plt.plot(xx, unwrap_phase(np.angle(obj[row])),'g.-')

    yyy = np.angle(obj[row])
    yyy[33:] -= 2*np.pi 
    yyy[54:] -= 2*np.pi 
    #plt.plot(xx, yyy,'.-')
    #plt.plot( np.angle(obj[34])-2*np.pi,'.-')
    #plt.plot( np.angle(obj[34,:]),'.-')
    plt.legend(['Recons phase wrapped','Recons phase unwrapped','Recons phase unwrapped II'],loc ='upper left')


    plt.figure()
    plt.imshow(np.angle(obj), cmap='jet', extent=extent_zoomed, interpolation='none')
    plt.colorbar()


    plt.figure()
    plt.imshow(np.abs(obj),  interpolation='none')
    
    plt.axhline(y=row,color='red')
    
    def box_function(x,limit1,limit2,low,high): 
        func = high*np.ones((len(x)))
        func[np.where(x <= limit1)] = low
        func[np.where(x > limit2)] = low 
        return func
    
    from scipy.optimize import differential_evolution
    
    y = (np.abs(obj[row,:]))/np.max(np.abs(obj[row,:]))
    #normalization has no effect
    a=abs(obj)
    plt.figure()
    plt.imshow(a)
    #a[a<0.40]=0
    plt.figure()
    plt.imshow(a)
    
    y = a[row,:]#/np.max(a[row,:])
    
 #   x,limit1,limit2,low,high
    res_170 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
                                 [[0.2, 0.3], [0.4, 0.5], [0.1, 0.4], [0.4, 1.5]])  # parameter bounds (ranges))
    res_80 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
                                 [[-0.3, -0.1], [-0.4, 0.05], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
    res_45 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
                                 [[-0.5, -0.4], [-0.4, -0.3], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
    res_19 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
                                 [[-0.80, -0.75], [-0.75, -0.63], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
    res_8 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
                                 [[-1.06, -1.015], [-1.015, -0.958], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
    
#    #temp sim
#    res_170 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
#                                 [[3.2, 3.4], [3.4, 4], [0.0, 0.4], [0, 1]])  # parameter bounds (ranges))
#    res_80 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
#                                 [[-0.3, -0.1], [-0.4, 0.05], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
#    res_45 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
#                                 [[-0.5, -0.4], [-0.4, -0.3], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
#    res_19 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
#                                 [[-0.80, -0.75], [-0.75, -0.63], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))
#    res_8 = differential_evolution(lambda p: np.sum((box_function(xx, *p) - y)**2),  # quadratic cost function
#                                 [[-1.06, -1.015], [-1.015, -0.958], [0.0, 0.1], [0.1, 1.5]])  # parameter bounds (ranges))


    #import pdb; pdb.set_trace()
    
    #plt.twinx()
    plt.figure()
    #plt.title('fitted box functions')
    plt.plot(xx, y,'b.-')
    plt.plot(xx,box_function(xx,*res_170.x),'r')
    length_170 = round(1000*(res_170.x[1]-res_170.x[0]),2) #calc length in nm
    plt.text(3.6,2.3E-5,str(length_170) + r' nm')
    #plt.text(0.46,0.7,str(length_170) + r' nm')
#    
#    plt.plot(xx,box_function(xx,*res_80.x),'r')
#    length_80 = round(1000*(res_80.x[1]-res_80.x[0]),2) #calc length in nm
#    plt.text(0.0,0.7,str(length_80) + r' nm')
#    
#    plt.plot(xx,box_function(xx,*res_45.x),'r')
#    length_45 = round(1000*(res_45.x[1]-res_45.x[0]),2) #calc length in nm
#    plt.text(-0.4,0.65,str(length_45) + r' nm')
#
#    plt.plot(xx,box_function(xx,*res_19.x),'r')
#    length_19 = round(1000*(res_19.x[1]-res_19.x[0]),2) #calc length in nm
#    plt.text(-0.7,0.6,str(length_19) + r" nm")
#    
#    plt.plot(xx,box_function(xx,*res_8.x),'r')
#    length_8 = round(1000*(res_8.x[1]-res_8.x[0]),2) #calc length in nm
#    print(res_8)
#    plt.text(-1,0.5,str(length_8) + r" nm")
#    
    plt.grid(True)
    plt.xlabel('$\mu$m')
    
    if save is True:
        fn = savepath + '\\' + 'amp_lineplot' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        
#plot_line_plot(obj[slicey,slice(198,400)], xx, row=16, save=False) #row=34
#plot_line_plot(masked_strain*100, row=34, save=False)


        
#%%
xx = np.linspace(0,obj[slicey,slicex].shape[1]*psize,obj[slicey,slicex].shape[1])
#plt.close('all')
def estimate_resolution_mtf(obj,xx,row):
    
    line = abs(obj[row])

    plt.figure()
    plt.plot(line,'r.')
    plt.plot( line,'b-')
    
    #import pdb; pdb.set_trace()
    plt.figure()
    #plt.title('fitted box functions')
    plt.plot(xx, line,'r.')
    plt.plot(xx, line,'b-')
    
    #Zoom in on right edge of 80nm segment
    s80_right_edge = line[123:148]
    xx80 = xx[123:148]
    #import pdb; pdb.set_trace()
    plt.figure()
    plt.title('edge')
    plt.plot(xx80, s80_right_edge,'r.')
    plt.plot(xx80, s80_right_edge,'b-')
    
    
    #ignore variations in amplitude inside segment
    #s80_right_edge[0:9+6] = s80_right_edge[9+6]
    
    plt.figure()
    plt.title('modified edge')
    plt.plot(xx80, s80_right_edge,'r.')
    plt.plot(xx80, s80_right_edge,'b-')
    
    #iterlude: fit edge function to modifide edge
    def edge_function(x,limit1,low,high): 
        func = high*np.ones((len(x)))
        func[np.where(x > limit1)] = low 
        return func
        
    
    from scipy.optimize import differential_evolution

    
    #x,limit1,limit2,low,high
    res_170 = differential_evolution(lambda p: np.sum((edge_function(xx80, *p) - s80_right_edge)**2),  # quadratic cost function
                                 [[xx80[0], xx80[-1]], [0.001, 0.4], [0.4, 1.5]])  # parameter bounds (ranges))
    plt.figure()
    plt.plot(xx80, s80_right_edge,'r.')
    plt.plot(xx80,edge_function(xx80,*res_170.x))
    
    grad_line = np.gradient(s80_right_edge, psize)
    
    plt.figure()
    plt.title('gradient of edge')
    plt.plot(xx80, grad_line,'r.')
    plt.plot(xx80, grad_line,'b-')
    
    #calc fft
    fft_grad_line = np.fft.fftshift(np.fft.fft(grad_line)) 
    #import pdb; pdb.set_trace()
    
    plt.figure()
    plt.title('fft (grad (right edge)) ')
    plt.plot(2*np.pi/xx80,  fft_grad_line,'r.-')
    
    plt.figure()
    plt.title('fft ( (right edge)) ')
    plt.plot(  fft_grad_line,'r.')
    plt.plot(  fft_grad_line,'b-')
    
    plt.figure()
    plt.title('Norm pos fft (grad (right edge)) ')
    #use np.fft.rfft which is a transform for real data and will return half the full FFT
    plt.plot(2*np.pi/xx80[11:],  np.fft.rfft(grad_line)/np.fft.rfft(grad_line).max(),'r.-')
    #plt.plot(2*np.pi/xx80[16:], fft_grad_line[16:]/np.fft.fft(grad_line).max(),'b-')
    # temp test skalan
    #import pdb; pdb.set_trace()
    plt.figure()
    plt.title('Norm pos fft (grad (right edge)) ')
    plt.plot(  fft_grad_line[16:]/np.fft.fft(grad_line).max(),'r.')
    plt.plot(fft_grad_line[16:]/np.fft.fft(grad_line).max(),'b-')


#plt.close('all')
lineout_row = int(obj[slicey,slicex].shape[0]/2)
#estimate_resolution_mtf(abs(obj[slicey,slicex]),xx,row=lineout_row)

#%%
#
from scipy.optimize import curve_fit
from scipy import special


#ska vara normaliserat, alltså i , ska vara minst 1


#def erf(z,a,b):
#   return a*special.erf(z)+b

def erf(x, a, b):
    return  (special.erf((x - a)/b) +1) /2

#def erf(x, mFL, a, b):
#    return mFL*special.erf((x-a)/(b*np.sqrt(2)))

#resolution estimate 2
# fit error runction to edge
line = abs(obj[slicey,slicex][lineout_row])    
xx = np.linspace(0,obj[slicey,slicex].shape[1]*psize,obj[slicey,slicex].shape[1])

#fig, ax = plt.subplots()
#ax.plot(abs(obj[slicey][lineout_row]) )

fig, ax = plt.subplots()
ax.plot(line,'+')


s170_ledge = line[154-8:154+8]
xx170 = xx[154-8:154+8]
#find 170 segment edge
fig, ax = plt.subplots()
ax.plot(s170_ledge,'+')     #S 439

#s170_ledge[15:] = s170_ledge[15]
#normalize amplitude
s170_ledge = s170_ledge/s170_ledge.max()

fig, ax = plt.subplots()
ax.set_title('normalised amplitude 170 ignoring amlitude var in segment')
ax.plot(s170_ledge)




#normalized amplitude  and flip to a left edge
s80_right_edge = (line[118:135][::-1] / line[118:135][::-1].max())
xx80 = xx[118:135]

xx_cen = ((xx80 - xx80[int(xx80.shape[0]/2)]) *1E9)
xx_cen_170 = ((xx170 - xx170[int(xx170.shape[0]/2)]) *1E9)


popt, pcov = curve_fit(erf, xx_cen, s80_right_edge,p0=[0, 50])
popt_170, pcov = curve_fit(erf, xx_cen_170, s170_ledge,p0=[0, 50])

#call function with fitted parameters
erf_fit = erf(xx_cen,*popt)

erf_fit_170 = erf(xx_cen_170,*popt_170)

#plt.close('all')

fig, ax = plt.subplots()
ax.plot(xx_cen,s80_right_edge,'r+')
ax.plot(xx_cen,erf_fit)

fig, ax = plt.subplots()
ax.plot((xx170 - xx170[0])*1E6,s170_ledge,'+')
ax.plot((xx170 - xx170[0])*1E6,erf_fit_170)
ax.legend(['Amplitude at edge','Error function fit'])

#fig, ax = plt.subplots()
#ax.plot(xx_cen,np.gradient(s80_right_edge),'r+')


#fig, ax = plt.subplots()
#ax.plot(xx_cen,erf_fit[::-1])


fig, ax = plt.subplots()
ax.set_title('80nm semgnet fit gradient')
ax.plot((xx80 - xx80[0])*1E6,np.gradient(erf_fit)/np.gradient(erf_fit).max())

fig, ax = plt.subplots()
#ax.set_title('170nm semgnet fit gradient')
ax.plot((xx170 - xx170[0])*1E6,np.gradient(erf_fit_170)/np.gradient(erf_fit_170).max(),'orange')


    