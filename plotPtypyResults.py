import numpy as np
import matplotlib.pyplot as plt
import nmutils
import h5py
import matplotlib.gridspec as gridspec
#import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import time
date_str = time.strftime("%Y%m%d_%H%M")
import os 
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
scan = 445

nam = str(scan) + '_20210921_1525' #437
#nam = str(scan) + '_20210709_1455' #442
#nam =  str(scan) + '_20210709_1605' #444
#nam =  str(scan) + '_20210712_1053' #445
#nam = str(scan) + '_20210814_1836' #446
#nam = str(scan) + '_20210814_2154' #447
#nam = str(scan) + '_20210819_1432' #444 probe not updating

name =  "\\" + nam + "\\" + nam + '_EPIE_%04u' % iterations

outputSuffix = 'png'
save = False

inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\2drecons' + folder + name + '.ptyr'
savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\2drecons\plots' + name

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
            print('saved error %d'%ite)
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
    


extent = 1e6 * np.array([origin[0], origin[0]+(obj.shape[1]-1)*psize, origin[1], origin[1]+(obj.shape[0]-1)*psize])

### Object
#%%
def plot2drecons(obj, probe, extent,savepath, save):
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
    slicey = slice(140,200)  #for real data recons with shape 170
    slicex = slice(214,438) #for real data recons with shape 170


    #slicey = slice(90,130) #for simulated data recons (should be the same, but for now)
    #slicex = slice(150,290)
    
    #mag_cut = mag[mag.shape[0]//3:2*mag.shape[0]//3, mag.shape[1]//3:2*mag.shape[1]//3] # to find relevant dynamic range
    mag_cut = mag[slicey,slicex] # to find relevant dynamic range
    #plt.figure()
    #plt.imshow(mag_cut)
    #plt.savefig('h6')
    #import pdb; pdb.set_trace()
    vmin = mag_cut.min()
    vmax = mag_cut.max()
    print('max ', vmax)
    print('vmin', vmin)
    
    #abs
    img = ax[0].imshow((mag), cmap='gray', extent=extent, vmin=vmin, vmax=vmax, interpolation='none')
    plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=70)
    ax[0].set_ylabel('$\mu$m')
    ax[0].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
    ax[0].set_title('Amplitude')
    
    # phase
    img = ax[1].imshow(np.angle(obj), vmin=-np.pi, vmax=np.pi, cmap='hsv', extent=extent, interpolation='none')
    #img = ax[1].imshow((np.angle(obj)), vmin=-2.5, vmax=-2, cmap='hsv', extent=extent, interpolation='none')
    plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=70 )
    ax[1].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax[1].set_title('Phase')
    
    #plt.show()
    #import pdb; pdb.set_trace()
    
    #zoomed in version
    global extent_zoomed
    extent_zoomed = 1e6 * np.array([origin[0], origin[0]+(mag_cut.shape[1]-1)*psize, origin[1], origin[1]+(mag_cut.shape[0]-1)*psize])
    
    #mask with amplitude  
    amplitude_range = vmax-vmin
    #mask_at = -1.5#0.05 * (vmin + amplitude_range)
    mask_at = 0.05 * (vmin + amplitude_range)
    mask = np.zeros((mag_cut.shape))
    mask[(np.abs(mag_cut)) > mask_at] = 1
    #mask[np.log10(np.abs(mag_cut)) > mask_at] = 1
    mask = np.squeeze(mask)
    
    
    fig, ax = plt.subplots(ncols=2, figsize=(10,3), sharex=True, sharey=True)
    #fig, ax = plt.subplots(ncols=2, figsize=(10,6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=.3)
    vmin = mag_cut.min()
    vmax = mag_cut.max()
    
    #abs
    img = ax[0].imshow(np.log10(mag_cut), extent=extent_zoomed, interpolation='none') #, vmin=vmin, vmax=vmax)#
    plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=70)
    ax[0].set_ylabel('$\mu$m')
    ax[0].set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)
    ax[0].set_title('Amplitude')
    
    # phase
    img = ax[1].imshow(mask* np.angle(obj[slicey,slicex]), cmap='jet', extent=extent_zoomed, interpolation='none') #vmin=-np.pi, vmax=np.pi)#,
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
    #70 segment zoomed
    img = ax.imshow( (mask* ( np.angle(obj[slicey,slicex]))), cmap='jet', interpolation='none', extent=extent_zoomed)#, vmin=-np.pi, vmax=np.pi)#)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude')
    #plt.savefig('hyyhy6')
    plt.show()
    
    
    if save is True:
        fn = savepath + '\\' + 'phase_object_zoomed' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    
    
    #calc strain and relative strain
    #--------------------------------------------------------------------
    
    
    
    InP_Qvect = 18543793660.27452
    
    # ska man använda den Q-vektor som sätts at theta, alltså det theta om mäts upp ~~8 deg istället för detta som motsvarar 10.54
    
    dz = psize
    
    #masked_phase = mask *  np.angle(np.squeeze(obj))
    # disaplacement in z, du/dz är det jag vill plotta
    # phase phi = u*Q
    #calculate the gradient before you mask other wise you
    # will get hard gradients at the mask edges naturally
    
    
    #strain = np.diff(np.angle(obj[slicey, slicex]) /InP_Qvect , axis = 1, append=0) /dz
    
    strain = np.gradient((((np.angle(obj[slicey, slicex])))/InP_Qvect) , dz)[1] 
    global masked_strain
    masked_strain = mask* strain
    
    ##Average over the reigion of interest (that is a bit arbitrary)
    #mean_strain = np.mean(masked_strain[np.nonzero(masked_strain)])
    #print('mean value in the masked strain',mean_strain)
    
    vs_max = 0.75
    vs_min =-0.75
    
    # plot strain
    fig, ax = plt.subplots(ncols=1)
    img = ax.imshow(masked_strain*100, cmap='RdBu_r')#, extent=extent_zoomed, interpolation='none')#, vmin=vs_min, vmax=vs_max)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    #plt.plot(obj.shape[1]*[1])
    plt.axhline(y=35,color='red')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Strain [%]')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'strain' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
    
    #---------------------------------------------
    # plot zoomed version (zoomed on 2 largest segments)     
    slicey2 = slice(8,50)
    slicex2 = slice(0,100)

    #slicey2 = slice(8,30) #for shape 170
    #slicex2 = slice(0,100) #for shape 170
    print('shaaaaape' , masked_strain.shape)
    
    #slicey2 = slice(8,30)    #sim data
    #slicex2 = slice(12,66)    #sim data
    
    masked_strain_zoomed = masked_strain[slicey2,slicex2]
    extent_zoomed2 = 1e6 * np.array([origin[0], origin[0]+(masked_strain_zoomed.shape[1]-1)*psize, origin[1], origin[1]+(masked_strain_zoomed.shape[0]-1)*psize])
    
    fig, ax = plt.subplots(ncols=1)
    img = ax.imshow(masked_strain_zoomed*100, cmap='RdBu_r', extent=extent_zoomed2, interpolation='none', vmin=vs_min, vmax=vs_max)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Strain [%]')
    #plt.savefig('strain')
    plt.show()
    
    if save is True:
        fn = savepath + '\\' + 'zoomed_strain' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
        
    #extent_zoomed3 =     
  
    #phase of 170 nm segment
    slicey3 = slice(15,-5)
    slicex3 = slice(28,62)
    fig, ax = plt.subplots()
    #70 segment zoomed
    img = ax.imshow( ( np.angle(obj[slicey,slicex][slicey2,slicex2][slicey3,slicex3])), cmap='jet', interpolation='none')#, vmin=-np.pi, vmax=np.pi)#)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude')
    #plt.savefig('hyyhy6')
    plt.show()
      
    #phase of 170 nm segment
    fig, ax = plt.subplots()
    #70 segment zoomed
    img = ax.imshow( unwrap_phase( np.angle(obj[slicey,slicex][slicey2,slicex2][slicey3,slicex3])), cmap='jet', interpolation='none')#, vmin=-np.pi, vmax=np.pi)#)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_ylabel('$\mu$m')
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Phase masked with log10 amplitude')
    #plt.savefig('hyyhy6')
    plt.show()
    
    #####OBSOBSOBS
    """ OBSOBSOBS TEST"""
    strain_170 = np.gradient(((    unwrap_phase(np.angle(obj[slicey,slicex][slicey2,slicex2][slicey3,slicex3])))/InP_Qvect) , dz)[1] 
    
    fig, ax = plt.subplots(ncols=1)
    img = ax.imshow(strain_170*100, cmap='RdBu_r', interpolation='none')#, vmin=vs_min)  #extent=extent_zoomed2
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
    ax.set_xlabel('$\mu$m')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(img, cax=cax)#, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
    #cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
    ax.set_title('Strain [%]')
    #plt.savefig('strain')
    plt.show()

from skimage.restoration import unwrap_phase  
#plt.close('all')      
plot2drecons(obj, probe, extent,savepath, save)

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
    
plot_recon_error(err1, save)
#plot_recon_error(err2, save)
#%%
xx = np.linspace(extent_zoomed[0],extent_zoomed[1],obj[slicey,slicex].shape[1])

def plot_line_plot(obj,save):
    
    plt.figure()
    #Line in 2d,  mask the phase with amplitude
    #plt.title('Lineplot of phase wrapped and unwrapped')
    #plt.title('Line plot of strain [%]')
    plt.plot(xx, (np.angle(obj[34,:])),'.-')
    
    #plt.plot(xx, np.unwrap(np.angle(obj[34,:])),'.-')
    yyy = np.angle(obj[34])
    yyy[33:] -= 2*np.pi 
    yyy[54:] -= 2*np.pi 
    #plt.plot(xx, yyy,'.-')
    #plt.plot( np.angle(obj[34])-2*np.pi,'.-')
    #plt.plot( np.angle(obj[34,:]),'.-')
    #plt.legend(['Recons phase wrapped','Recons phase unwrapped'],loc ='upper left')
    
    plt.twinx()
    plt.plot(xx, (np.abs(obj[34,:])),'r.-')
    #plt.legend(['Recon amplitude'],loc ='upper right')
    plt.xlabel('$\mu$m')
    
    if save is True:
        fn = savepath + '\\' + 'amp_lineplot' + '.' + outputSuffix
        plt.savefig(fn)
        print("Saved to %s"%fn)
plot_line_plot(obj[slicey,slicex],save=False)
#plot_line_plot(masked_strain*100,save=False)

obj[slicey,slicex].shape
    


