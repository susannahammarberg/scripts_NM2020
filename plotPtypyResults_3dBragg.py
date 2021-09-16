import numpy as np
import matplotlib.pyplot as plt
#import nmutils
import h5py
import matplotlib.gridspec as gridspec
#import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

"""
This script visualizes the output of a ptypy run, by loading a ptyr file.

<SH> this works for python 3


print out the error! ssa?

import matplotlib as mpl
mpl.rcParams['toolbar'] = 'None'

"""
iterations = 50
folder = r'\dumps'  #will only have 1 error value (the final one)
folder = r'\recons'   
nam =  'recon20210618_1541'
name =  "\\" + nam + "\\" + nam + '_DM_3dBragg_%04u' % iterations
###name =  "\\" + nam + "\\" + nam + '_DM_%04u' % iterations

###name = "\\" + nam + "\\" + nam +' _None_0000'


inputFile = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons' + folder + name + '.ptyr'
savepath = r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\3drecons\plots' + name


if not os.path.exists(savepath):
    os.makedirs(savepath)
    print('new folder in this savepath was created')
    

print('Opening file: \n ', inputFile)

err_types = ['00000', '00001', '00002']

### load reconstruction data
with h5py.File(inputFile, 'r') as hf:
    scanid = str(hf['content/probe'].keys()).split("['", 1)[1].split("']")[0]
    print( 'loading entry %s' % scanid)
    probe = np.array(hf.get('content/probe/%s/data' % scanid)) # 'Sscan01G00'
    obj = np.array(hf.get('content/obj/%s/data' % scanid))
    test = np.array(hf.get('content/dif/%s/data' % scanid))
    psize = np.array(hf.get('content/probe/%s/_psize' % scanid))
    energy = np.array(hf.get('content/probe/%s/_energy' % scanid))
    origin = np.array(hf.get('content/probe/%s/_origin' % scanid))

    errors = []
    for ite in range(0,iterations):
        hh = np.array(hf.get('content/runtime/iter_info/%s/error'%('{0:05}'.format(ite))))        
        if hh.size == 3:
            errors.append(hh)
            print('saved error %d'%ite)
           

# reformat to np array
errors = np.array(errors)
try:
    probe = probe[0]
    obj = obj[0]
    psize = psize
except IndexError:
    raise IOError('That doesn''t look like a valid reconstruction file!')

#TODO energy is not working
print( "Loaded probe %d x %d x %d and object %d x %d x %d"%(probe.shape + obj.shape ))
print( "pixel size %.f x %.f x %.f nm"%(psize[0]*1e9, psize[1]*1e9, psize[2]*1e9))

#print( "Loaded probe %d x %d and object %d x %d, pixel size %.1f nm"%(probe.shape[0] + obj.shape[0] + psize*1e9))

print( "object stored as (x, z, y) where x is long axis along NW and x rotation axis")
       

print('origin',origin)

#extent_zy = 1e6 * np.array([origin[1], origin[1]+(obj.shape[1]-1)*psize[1],       origin[2], origin[2]+(obj.shape[2]-1)*psize[2]])
#extent_xz = 1e6 * np.array([  origin[0], origin[0]+(obj.shape[0]-1)*psize[0],     origin[1], origin[1]+(obj.shape[1]-1)*psize[1]]) 
#extent_xy = 1e6 * np.array([origin[0], origin[0]+(obj.shape[0]-1)*psize[0],       origin[2], origin[2]+(obj.shape[2]-1)*psize[2]])

extent_zy = 1e6 * np.array([-(obj.shape[1]/2)*psize[1],(obj.shape[1]/2)*psize[1],  -(obj.shape[2])/2*psize[2],(obj.shape[2])/2*psize[2]])
extent_xz = 1e6 * np.array([ -(obj.shape[0]/2)*psize[0],(obj.shape[0]/2)*psize[0], -(obj.shape[1]/2)*psize[1],(obj.shape[1]/2)*psize[1]]) 
extent_xy = 1e6 * np.array([-(obj.shape[0]/2)*psize[0],(obj.shape[0]/2)*psize[0],  -(obj.shape[2]/2)*psize[2], (obj.shape[2]/2)*psize[2]])
#import pdb; pdb.set_trace()

extent_yz = 1e6 * np.array([      (obj.shape[2])/2*psize[2],-(obj.shape[2])/2*psize[2], -(obj.shape[1]/2)*psize[1],(obj.shape[1]/2)*psize[1]])
extent_zx = 1e6 * np.array([ -(obj.shape[1]/2)*psize[1],(obj.shape[1]/2)*psize[1], -(obj.shape[0]/2)*psize[0],(obj.shape[0]/2)*psize[0]]) 
#extent_xz = 0 1
#extent_xy = 02
#extent_zy = 12

#which one?
extent_probe = 1e6 * np.array([origin[1], origin[1]+(probe.shape[1]-1)*psize[1], origin[2], origin[2]+(probe.shape[2]-1)*psize[2]])
extent_probe = 1e6 * np.array([origin[1], origin[1]+(probe.shape[2]-1)*psize[2], origin[2], origin[2]+(probe.shape[1]-1)*psize[1]])
#Extent defines the left and right limits, and the bottom and top limits


xcut = int(obj.shape[0]/2)
zcut = int(obj.shape[1]/2)
ycut = int(obj.shape[2]/2)

#x, z, y = S_cart.grids()
#%%

#plot in a more narrow region around the object
#################################################
lenx = int(obj.shape[0])
lenz = int(obj.shape[1])
leny = int(obj.shape[2])

rangey = int(leny*0.04)
rangez = int(lenz*0.12)

#Tänk på att det fortfarande ska vara centerat kring 0
slicey = slice(ycut-rangey, ycut+rangey)
slicez = slice(zcut-10-rangez, zcut-10+rangez)

#NM
#slicey = slice(ycut-rangey-15, ycut+rangey+15)    #170x170 pixel roi
slicey = slice(ycut-rangey-50, ycut+rangey+50)
slicez = slice(zcut-rangez-85, zcut+rangez+85)   #170x170 pixel roi

slicez = slice(zcut-rangez-75, zcut+rangez+75)

###################################################################################################################OBSOBS
# for 170nm segment
#zcut = 160 #for 170x170 pixel roi on detector
#zcut = 90 # for 96x96


##_____________________________________________________________________________
print('xcut',xcut)
print('slicez',slicez)

max_scatt = np.max(abs(obj[xcut-4:xcut+4,slicez,100:115]))  #for 170x170 pixel roi
#max_scatt = np.max(abs(obj[:,slicez,:]))  #for 96x96 pixel roi

# cut out region of  object
a_vv=obj[xcut-4:xcut+4,slicez,100:115] #for 170x170 pixel roi
#a_vv=obj[xcut-4:xcut+4,slicez,15:23] #for 96x96 pixel roi
# min scatt intensity in cut out region that is not 0
#min_scatt  = np.min(abs(a_vv[np.nonzero(a_vv)]))
min_scatt = 0
print('max and min scatt intensity in NW region:')
print(max_scatt,min_scatt)

print('central  x ', xcut)
##_____________________________________________________________________________

mask = np.zeros((obj.shape))
# make mask based on a % of the intensity of the max
mask[abs(obj)>(max_scatt/50)] = 1 #0.0002
#mask[abs(obj)>(0.50E-6)] = 1 



#calc strain and relative strain
#--------------------------------------------------------------------

#Average over the reigion of interest (that is a bit arbitrary)
#TODO not plotting strain no
#strain_xcut =  np.angle(np.squeeze(obj)[xcut][slicey, slicez])
InP_Qvect = 18543793660.27452

# ska man använda den Q-vektor som sätts at theta, alltså det theta om mäts upp ~~8 deg istället för detta som motsvarar 10.54

dz = psize[1]

#masked_phase = mask *  np.angle(np.squeeze(obj))
#disaplacement in z, du/dz är det jag vill plotta
# phase phi = u*Q
strain_xcut = np.diff(np.angle(obj[xcut][slicez, slicey]) /InP_Qvect ,n=1, axis = 0, append=0) /dz

masked_strain = mask[xcut][slicez, slicey] * strain_xcut


mean_strain = np.mean(masked_strain[np.nonzero(masked_strain)])
print('mean value in the masked strain',mean_strain)

rel_strain_xcut = 100*(masked_strain - mean_strain) / mean_strain

#xmin xmax ymin  ymax 
extent_zy_cut = 1e6 * np.array([-(strain_xcut.shape[0]/2)*psize[1],(strain_xcut.shape[0]/2)*psize[1],      (strain_xcut.shape[1]-1)/2*psize[2],-(strain_xcut.shape[1]-1)/2*psize[2]])




#plot the amplitude and strain (strain, rel strain or just phase) of the zoomed in region

## Amplitude
fig, ax = plt.subplots(ncols=1)  
im3 = ax.imshow(np.abs(obj[xcut][slicez, slicey]).T, interpolation='none', cmap='RdBu_r', extent = extent_zy_cut)#, vmin=min_scatt, vmax=max_scatt)
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='Amplitude \n front view')
fig.colorbar(im3,ax=ax )
plt.draw()
plt.savefig(savepath + '\\6')
plt.show()

## Phase
fig, ax = plt.subplots()
#im2 = ax.imshow((mask[xcut][slicez, slicey] *np.angle(obj[xcut][slicez, slicey])).T, cmap='jet', interpolation='none', origin='lower', extent = extent_zy_cut)
im2 = ax.imshow((np.angle(obj[xcut][slicez, slicey])).T, cmap='jet', interpolation='none', origin='lower', extent = extent_zy_cut)#, vmax = 0.25, vmin = -1)
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='Phase\n front view. fixed scale')
#cbar = plt.colorbar(im1, ax=ax1)
#cbar.ax.set_yticklabels(['-$\pi/8$','0','-$\pi$/8'])
cbar = plt.colorbar(im2, ax=ax) #ticks=(-np.pi/8, 0, np.pi/8)
fig.tight_layout()
plt.savefig(savepath + '\\6_5_5')
plt.show()


#Interlidue---------------------------------------------------------------------
## find region of the zoomed in region
## in the other dimention (the 3rd is set by )
#fig, ax = plt.subplots(ncols=1)  
#im3 = ax.imshow(np.flipud(np.abs(obj[16:24,slicez,ycut].T)), interpolation='none', origin='lower', cmap='jet')#, extent = extent_zy_cut)
#plt.setp(ax, ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='Amplitude \n  view')
#fig.colorbar(im3,ax=ax )
#plt.draw()
#plt.show()
## in the other dimention (the 3rd is set by )
#fig, ax = plt.subplots(ncols=1)  
#im3 = ax.imshow(np.flipud(np.abs(obj[xcut,slicez,100:115].T)), interpolation='none', origin='lower', cmap='jet',vmax = max_scatt,vmin=min_scatt)#, extent = extent_zy_cut)
#plt.setp(ax, ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='Amplitude \n  view')
#fig.colorbar(im3,ax=ax )
#plt.draw()
#plt.show()
#End Interlidue---------------------------------------------------------------------



##this is good
###### strain
####strain_xcut
fig, ax = plt.subplots(ncols=1)  
im3 = ax.imshow(((masked_strain.T)), interpolation='none', origin='lower', cmap='RdBu_r', extent = extent_zy_cut, vmin=-0.002, vmax=0.0014)
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='Axial strain du/dz \n front view (masked) ')
fig.colorbar(im3,ax=ax )
#ax.axvline(x=0.5,ymin=0,ymax=1,color='white', ls = '--')
#ax.axvline(x=-0.5,ymin=0,ymax=1,color='white', ls = '--')
plt.draw()
plt.savefig(savepath + '\\6_6')
plt.show()


# medelvärdet räknas ut med de höger värderna som är i den här kartan
## Relative strain
fig, ax = plt.subplots(ncols=1)  
im3 = ax.imshow(((rel_strain_xcut.T)), interpolation='none', origin='lower', cmap='RdBu_r', extent = extent_zy_cut, vmin=-1100, vmax=700)
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='Relative axial strain [%] \n front view')
fig.colorbar(im3,ax=ax )
#plt.draw()
plt.savefig(savepath + '\\6_7')
#plt.show()
print('Files saved to %s  '%(savepath))






#fig, (ax1, ax2) = plt.subplots(1,2) ; plt.suptitle('Side view')
#ax1.imshow(np.transpose(np.log10(abs(obj[:,:,ycut]))), cmap='RdBu_r',  interpolation='none', extent = extent_xz) ; ax1.set_title('Obj amplitude')
#ax1.set_ylabel('$z$ [$\mathrm{\mu}$m]') ; ax1.set_xlabel(('$x$ [$\mathrm{\mu}$m]')) 
#ax2.imshow(np.transpose((np.angle(obj[:,:,ycut]))), cmap='jet', interpolation='none', extent = extent_xz); ax2.set_title('Obj phase')
#ax2.set_ylabel(('$z$ [$\mathrm{\mu}$m]')); ax2.set_xlabel(('$x$ [$\mathrm{\mu}$m]')); #ax2.set_xticklabels(['',''])
##cbar = plt.colorbar(ticks=(-np.pi/8, 0, np.pi/8))
##ax2.set_yticklabels(['-0.25','0.25'])
#fig.tight_layout()

#fig, (ax1, ax2) = plt.subplots(2,1); plt.suptitle('log amplitude. View from above')
#ax1.imshow(np.log10(abs(obj[:,zcut,:])), cmap='RdBu_r',  interpolation='none', extent = extent_xy) ; ax1.set_title('Obj amplitude')
#ax1.set_ylabel('$x$ [$\mathrm{\mu}$m]') ; ax1.set_xlabel(('$y$ [$\mathrm{\mu}$m]')) 
#ax2.imshow((np.angle(obj[:,zcut,:])), cmap='jet', interpolation='none', extent = extent_xy); ax2.set_title('Obj phase')
#ax2.set_ylabel('$x$ [$\mathrm{\mu}$m]') ; ax2.set_xlabel(('$y$ [$\mathrm{\mu}$m]')); #ax2.set_xticklabels(['',''])
##cbar = plt.colorbar(ticks=(-np.pi/8, 0, np.pi/8))
##ax2.set_yticklabels(['-0.25','0.25'])
#fig.tight_layout()
#

#%%
#import pdb;pdb.set_trace()

#_________________________________________________________________________________________________________
##plot the probe projections
##pcut=0
##for pcut in range(0,42):
##    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(5,5));
##    plt.suptitle('Cut %d'%pcut)
##    ax[0].imshow((abs(probe[pcut])), cmap='RdBu_r',extent = extent_probe, interpolation='none'); ax[0].set_title('Amp probe'); ax[0].set_xlabel('$\mu$m') 
##    ax[1].imshow(np.angle(probe[pcut]), extent = extent_probe, cmap='jet', interpolation='none'); ax[1].set_title('Phase probe'); ax[1].set_xlabel('$\mu$m')
##    plt.tight_layout()
##    plt.savefig(r'C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\probe\cut%d'%pcut)
###    plt.show()
#___________________________________________________________________________________________________

#Plot the probe . object is wrong
##fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(5.5,7));
##ax[0,0].imshow((abs(obj[xcut])), cmap='RdBu_r',extent = extent_zy, interpolation='none'); ax[0,0].set_title('amp obj'); ax[0,0].set_ylabel('$\mu$m'); ax[0,0].set_xlabel('$\mu$m'); 
##ax[0,1].imshow(np.angle(obj[xcut]),cmap='jet',extent = extent_zy, interpolation='none'); ax[0,1].set_title('Phase obj'); ax[0,1].set_xlabel('$\mu$m')
##ax[1, 0].imshow((abs(probe[xcut])), cmap='RdBu_r',extent = extent_probe, interpolation='none'); ax[1, 0].set_title('Amp probe'); ax[1, 0].set_xlabel('$\mu$m') 
##ax[1, 1].imshow(np.angle(probe[xcut]), extent = extent_probe, cmap='jet', interpolation='none'); ax[1, 1].set_title('Phase probe'); ax[1, 1].set_xlabel('$\mu$m')
##ax[2, 0].plot(np.log10(errors[:,0]/errors[:,0].max()),'o-')  ;  ax[2, 0].set_title('Error last iteration %f'%errors[-1][0]) ;  ax[2, 0].plot(np.log10(errors[:,2]/errors[:,2].max()),'yd'); ax[2, 0].legend(('log err1','log err3'), shadow=True)
##ax[2,1].remove()  # don't display empty ax
##plt.tight_layout()
##plt.savefig(savepath + '\\1_1')
###plt.show()
###plt.figure() ; plt.plot(np.log10(errors[:,0]/errors[:,0].max()),'o-'); plt.plot(errors[:,1]/errors[:,1].max(),'r+') ; plt.plot(errors[:,2]/errors[:,2].max(),'yd')    ; plt.legend(('err1','err2','err2'))
##
##




fig = plt.figure()    # row col  row col
ax0 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=1, rowspan=2)
ax3 = plt.subplot2grid((3, 3), (1, 1), colspan=1, rowspan=2)
ax4 = plt.subplot2grid((3, 3), (1, 2), colspan=1, rowspan=2)

ax0.plot(errors[1:,1],'blue',marker='.', label='Fourier error')
ax0.legend(bbox_to_anchor=(0.65, 0.5), loc='center left',)

#x, z, y = S_cart.grids()
#plt.suptitle('log amplitude means.', x=0.7)
im1 = ax2.imshow(np.log10(np.mean(np.abs(obj), axis=1)).T, extent=extent_xy, interpolation='none', origin='lower')
plt.setp(ax2, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='side view')
plt.colorbar(im1,ax=ax2)

im2 = ax3.imshow(np.log10(np.mean(np.abs(obj), axis=2)).T, extent=extent_xz, interpolation='none', origin='lower')
plt.setp(ax3, ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='top view')
plt.colorbar(im2,ax=ax3)

im3 = ax4.imshow(np.log10(np.mean(np.abs(obj) , axis=0)).T, extent=extent_zy, interpolation='none', origin='lower')
plt.setp(ax4, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
plt.colorbar(im3,ax=ax4)
plt.setp(ax2.xaxis.get_majorticklabels(), rotation=70)
plt.setp(ax3.xaxis.get_majorticklabels(), rotation=70)
#thight layout should minimize the overlab of labels
#plt.tight_layout() 
#plt.draw()
plt.savefig(savepath + '\\3')
#plt.show()

#______________________________________________________________________________________________________________
fig = plt.figure()    # row col  row col
ax0 = plt.subplot2grid((3, 3), (0, 0), colspan=3)
ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=1, rowspan=2)
ax3 = plt.subplot2grid((3, 3), (1, 1), colspan=1, rowspan=2)
ax4 = plt.subplot2grid((3, 3), (1, 2), colspan=1, rowspan=2)

ax0.plot(errors[1:,1],'blue',marker='.', label='Fourier error')
ax0.legend(bbox_to_anchor=(0.65, 0.5), loc='center left',)
plt.suptitle('Amplitude. Central cuts')#, x=0.4)
im1 = ax2.imshow( np.abs(obj[:,zcut]).T, extent=extent_xy, interpolation='none', origin='lower', cmap='RdBu_r')
plt.setp(ax2, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='side view')
plt.colorbar(im1,ax=ax2)

im2 = ax3.imshow( np.abs(obj[:,:,ycut]).T, extent=extent_xz, interpolation='none', origin='lower', cmap='RdBu_r')#,vmin=min_scatt, vmax=max_scatt)
plt.setp(ax3, ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='top view')
plt.colorbar(im2,ax=ax3)

im3 = ax4.imshow( np.abs(obj[xcut]).T , extent=extent_zy, interpolation='none', origin='lower', cmap='RdBu_r')#,vmin=min_scatt, vmax=max_scatt)
plt.setp(ax4, ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
plt.colorbar(im3,ax=ax4)
plt.setp(ax2.xaxis.get_majorticklabels(), rotation=70)
plt.setp(ax3.xaxis.get_majorticklabels(), rotation=70)


plt.tight_layout(rect=[0, 0.03, 1, 0.95])
#plt.tight_layout() 
plt.draw()
plt.savefig(savepath + '\\4')
plt.show()


#__________________________________________________________________________________________________
fig, ax = plt.subplots(ncols=3)
plt.suptitle('Phase. Central cuts.', x=0.8)
im1 = ax[0].imshow(( np.angle(obj[:,zcut])  ).T, extent=extent_xy, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax[0], ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='side view')
fig.colorbar(im1,ax=ax[0] )

im2 = ax[1].imshow(( np.angle(obj[:,:,ycut]) ).T, extent=extent_xz, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax[1], ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='top view')
fig.colorbar(im2,ax=ax[1] )

im3 = ax[2].imshow(( np.angle(obj[xcut]) ).T, extent=extent_zy, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax[2], ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
fig.colorbar(im3,ax=ax[2] )
plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=70)
plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=70)
plt.tight_layout() 
plt.draw()
plt.savefig(savepath + '\\5')
#plt.tight_layout() 
# dont show it just save it?
# plt.savefig("\savefig\iter%"%i)
# save as array and plot afterwards
plt.show()

fig, ax = plt.subplots(ncols=3)
plt.suptitle('Phase masked with amplitude. \n Central cuts.' , x=0.8)
im1 = ax[0].imshow(( mask[:,zcut] *np.angle(np.squeeze(obj)[:,zcut])  ).T, extent=extent_xy, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax[0], ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]', title='side view')
fig.colorbar(im1,ax=ax[0] )

im2 = ax[1].imshow((mask[:,:,ycut] *np.angle(np.squeeze(obj)[:,:,ycut]) ).T, extent=extent_xz, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax[1], ylabel=r'z [$\mu$m]', xlabel=r'x [$\mu$m]', title='top view')
fig.colorbar(im2,ax=ax[1] )

im3 = ax[2].imshow(( mask[xcut] *np.angle(np.squeeze(obj)[xcut]) ).T, extent=extent_zy, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax[2], ylabel=r'y [$\mu$m]', xlabel=r'z [$\mu$m]', title='front view')
fig.colorbar(im3,ax=ax[2] )
plt.tight_layout() 
plt.draw()
plt.savefig(savepath + '\\8')
plt.show()



#
##fig, ax = plt.subplots()
##im2 = ax.imshow((mask[:,:,ycut] *np.angle(np.squeeze(obj)[:,:,ycut])).T, extent=extent_zx, interpolation='none', origin='lower', cmap='jet')
##plt.setp(ax, ylabel=r'x [$\mu$m]', xlabel=r'z [$\mu$m]')
##fig.colorbar(im2,ax=ax)
##plt.draw()
##plt.savefig(savepath + '\\8_b')
##plt.show()

print(obj.shape)
fig, ax = plt.subplots()
im3 = ax.imshow(( mask[xcut] *np.angle(np.squeeze(obj)[xcut])).T, extent=extent_zy, interpolation='none', origin='lower', cmap='jet')
plt.setp(ax, ylabel=r'y [$\mu$m]', xlabel=r'x [$\mu$m]')
fig.colorbar(im3,ax=ax )
plt.tight_layout() 
plt.draw()
plt.savefig(savepath + '\\8_c')
plt.show()





#%%
#### define distances and propagate
#dist = np.linspace(backProp, forwProp, steps) * 1e-6
#dx = dist[1] - dist[0]
#print( "propagating to %d positions separated by %.1f um..."\
#    % (len(dist), dx*1e6))
#### not sure why, but the propagation goes in the other direction here!
#### it could be a misunderstanding about motors at nanomax...
#field3d = nmutils.utils.propagateNearfield(probe, psize, -dist, energy)
#
#### get intensities and focii
#power3d = np.abs(field3d)**2
#power_vertical = np.sum(power3d, axis=2).T
#power_horizontal = np.sum(power3d, axis=1).T
#focus_vertical_ind = np.argmax(np.sum(power_vertical**2, axis=0))
#focus_vertical_x = dist[focus_vertical_ind]
#focus_horizontal_ind = np.argmax(np.sum(power_horizontal**2, axis=0))
#focus_horizontal_x = dist[focus_horizontal_ind]
#focus_ind = np.argmax(np.sum(power3d**2, axis=(1,2)))
#focus_x = dist[focus_ind]
#focus = field3d[focus_ind]
#
#### plot
#fig = plt.figure(figsize=(8, 10))
#outer_grid = gridspec.GridSpec(2, 2, wspace=.2, hspace=.2, height_ratios=[2,3])
#
## probe and focus spots
#def spot_subplot(gridcell, data, shareax=None, title=''):
#    subgrid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gridcell, width_ratios=[3,1], height_ratios=[1,3], hspace=.05, wspace=.05)
#    lims = [-1e6*data.shape[0] * psize / 2, 1e6*data.shape[0] * psize / 2] # um
#    posrange = np.linspace(lims[0], lims[1], data.shape[0])
#    ax = plt.subplot(subgrid[1,0], sharex=shareax, sharey=shareax)
#    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
#    ax.imshow(nmutils.utils.complex2image(data), extent=lims+[lims[1], lims[0]], interpolation='none')
#    ax_histh = plt.subplot(subgrid[0,0], sharex=ax, yticklabels=[], ylabel='Int.')
#    ax_histv = plt.subplot(subgrid[1,1], sharey=ax, xticklabels=[], xlabel='Int.')
#    ax_histv.plot(np.sum(np.abs(data)**2, axis=1), posrange)
#    ax_histh.plot(posrange, np.sum(np.abs(data)**2, axis=0))
#    ax_histh.set_title(title, x=.67)
#    for tk in ax_histh.get_xticklabels(): tk.set_visible(False)
#    for tk in ax_histv.get_yticklabels(): tk.set_visible(False)
#
#    # FWHM:
#    import scipy.interpolate
#    for i in (0,1):
#        y = np.sum(np.abs(data)**2, axis=i)
#        edges = scipy.interpolate.UnivariateSpline(posrange, y-y.max()/2).roots()
#        r1, r2 = edges[0], edges[-1]
#        if i == 0:
#            ax_histh.axvspan(r1, r2, fc='r', alpha=.3)
#            ax_histh.text(r2, np.mean(ax_histh.get_ylim()), ' %.0f nm'%((r2-r1)*1e3), fontsize=10, rotation=-90*i)
#            ax.set_xlim(np.mean([r1,r2]) + np.array([-1,1]))
#        elif i == 1:
#            ax_histv.axhspan(r1, r2, fc='r', alpha=.3)
#            ax_histv.text(np.mean(ax_histv.get_xlim()), r1, ' %.0f nm'%((r2-r1)*1e3), fontsize=10, va='top', ha='center', rotation=-90)
#            ax.set_ylim(np.mean([r1,r2]) + np.array([-1,1]))
#    return ax
#a = spot_subplot(outer_grid[0,0], probe, title='Sample plane')
#a.set_ylabel('$\mu$m')
#b = spot_subplot(outer_grid[0,1], focus, shareax=a, title='Focal plane')
#
## beam profiles
#subgrid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_grid[1,:], hspace=.05)
#ax_vertical = plt.subplot(subgrid[0])
#ax_horizontal = plt.subplot(subgrid[1], sharex=ax_vertical, sharey=ax_vertical)
#ax_vertical.imshow(power_vertical, cmap='gray', extent=[1e6*dist[0], 1e6*dist[-1], -1e6*psize*probe.shape[0]/2, 1e6*psize*probe.shape[1]/2], interpolation='none')
#ax_vertical.axvline(x=focus_vertical_x*1e6, lw=1, ls='--', color='r')
#ax_vertical.text(1e6*focus_vertical_x, ax_vertical.get_ylim()[0] + .2*np.diff(ax_vertical.get_ylim())[0], 
#    ' %.0f um '%(1e6*focus_vertical_x), color='red', 
#    ha=('right' if focus_vertical_x<focus_x else 'left'))
#ax_vertical.axvline(x=focus_x*1e6, lw=2, ls='-', color='r')
#ax_vertical.text(1e6*focus_x, ax_vertical.get_ylim()[0] + .8*np.diff(ax_vertical.get_ylim())[0],
#    ' %.0f um '%(1e6*focus_x), color='red', va='top',
#    ha=('right' if focus_vertical_x>focus_x else 'left'))
#ax_vertical.set_aspect('auto')
#ax_horizontal.imshow(power_horizontal, cmap='gray', extent=[1e6*dist[0], 1e6*dist[-1], -1e6*psize*probe.shape[0]/2, 1e6*psize*probe.shape[1]/2], interpolation='none')
#ax_horizontal.axvline(x=focus_horizontal_x*1e6, lw=1, ls='--', color='r')
#ax_horizontal.text(1e6*focus_horizontal_x, ax_horizontal.get_ylim()[0] + .2*np.diff(ax_horizontal.get_ylim())[0],
#    ' %.0f um '%(1e6*focus_horizontal_x), color='red',
#    ha=('right' if focus_horizontal_x<focus_x else 'left'))
#ax_horizontal.axvline(x=focus_x*1e6, lw=2, ls='-', color='r')
#ax_horizontal.set_aspect('auto')
#ax_horizontal.set_ylabel('$\mu$m', y=1.05)
#for tk in ax_vertical.get_xticklabels(): tk.set_visible(False)
#ax_horizontal.set_xlabel('beamline z axis ($\mu$m)', fontsize=16)
#
#lblue = (.3,.3,1.0)
#ax_vertical.axvline(x=0, lw=1, ls='--', color=lblue)
#ax_horizontal.axvline(x=0, lw=1, ls='--', color=lblue)
#ax_horizontal.text(0, ax_horizontal.get_ylim()[0] + .01*np.diff(ax_horizontal.get_ylim())[0], 
#    ' sample ', color=lblue, va='bottom',
#    ha=('left' if focus_x < 0 else 'right'))
#
#ax_horizontal.text(1.05, 0.5, 'horizontal focus (M2)',
#     horizontalalignment='center',
#     verticalalignment='center',
#     rotation=90,
#     transform = ax_horizontal.transAxes)
#ax_vertical.text(1.05, 0.5, 'vertical focus (M1)',
#     horizontalalignment='center',
#     verticalalignment='center',
#     rotation=90,
#     transform = ax_vertical.transAxes)
#
#plt.suptitle(title, fontsize=20)
#if outputFile is not None:
#    fn = outputPrefix + '_probe.' + outputSuffix
#    plt.savefig(fn)
#    print( 'Saved to %s'%fn)
#
#### Object
#fig, ax = plt.subplots(ncols=2, figsize=(10,6), sharex=True, sharey=True)
#plt.subplots_adjust(wspace=.3)
#fig.suptitle('title', fontsize=20)
#extent = 1e6 * np.array([origin[0], origin[0]+(obj.shape[1]-1)*psize, origin[1], origin[1]+(obj.shape[0]-1)*psize])
##
### amplitude
#mag = np.abs(obj)
#mag_cut = mag[mag.shape[0]/3:2*mag.shape[0]/3, mag.shape[1]/3:2*mag.shape[1]/3] # to find relevant dynamic range
#vmin = mag_cut.min()
#vmax = mag_cut.max()
#img = ax[0].imshow(mag, cmap='gray', extent=extent, vmin=vmin, vmax=vmax, interpolation='none')
#plt.setp(ax[0].xaxis.get_majorticklabels(), rotation=70)
#ax[0].set_ylabel('$\mu$m')
#ax[0].set_xlabel('$\mu$m')
#divider = make_axes_locatable(ax[0])
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.colorbar(img, cax=cax)
#ax[0].set_title('Amplitude')

## phase
#img = ax[1].imshow(np.angle(obj), vmin=-np.pi, vmax=np.pi, cmap='hsv', extent=extent, interpolation='none')
##img = ax[1].imshow(np.angle(obj), vmin=-.1, vmax=.2, cmap='hsv', extent=extent, interpolation='none')
#plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=70 )
#ax[1].set_xlabel('$\mu$m')
#divider = make_axes_locatable(ax[1])
#cax = divider.append_axes("right", size="5%", pad=0.05)
#cb = plt.colorbar(img, cax=cax, ticks=(-np.pi, -np.pi/2, 0, np.pi/2, np.pi))
#cb.ax.set_yticklabels(['-$\pi$', '-$\pi/2$', '0', '$\pi/2$', '$\pi$'])
#ax[1].set_title('Phase')
#
#if outputFile is not None:
#    fn = outputPrefix + '_object.' + outputSuffix
#    plt.savefig(fn)
#    print( "Saved to %s"%fn)
#
#if interactive:
#    plt.show()
