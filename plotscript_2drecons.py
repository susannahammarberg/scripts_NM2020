# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 12:22:27 2021

use plotscrpt to plot simulated data with exp plotscript
@author: Sanna
"""

import numpy as np
import matplotlib.pyplot as plt
#from ptypy.core import View, Container, Storage, Base
import os
import sys
sys.path.append(r"C:\Users\sanna\Phd\NanoMAX_May2020\scripts_NM2020")
from plotPtypyResults import plot2drecons
from plotPtypyResults import plot_line_plot
from plotPtypyResults import estimate_resolution_mtf


import matplotlib
matplotlib.use( 'Qt5agg' )

#%%
#20220210   31

date_saved = 20220809  ###35: 20220809 #31: 20220210     #20220809#20220210#20220311#  33: 20220809#  # 20220809# 20220510 # 20220310#20220210 #20220902_projection29
##20220510_projection39
projection = 35

##nam = '2022_0310_33_
itstr = 'iter200' 


save = False
signflip = True


##openpath = r'C:\Users\Sanna\Documents\Simulations\save_simulation\NM2020_v1\recons\%s_projection%i'%(date_saved,projection)
#openpath = r'C:\Users\Sanna\Documents\Simulations\save_simulation\NM2020_updated_comsol\recons\%s_projection%i'%(date_saved,projection)

openpath = r'C:\Users\sanna\Phd\NanoMAX_May2020\Simulations\recons\%s_projection%i'%(date_saved,projection)
savepath = r'C:\Users\sanna\Phd\NanoMAX_May2020\Simulations\plots\%s_projection%i_%s'%(date_saved,projection,itstr)


#SF simulation:
#openpath = r'C:\Users\Sanna\Documents\Simulations\save_simulation\202201_NM2020_SF_simulation\recons\%s_projection%i'%(date_saved,projection)
#savepath = r'C:\Users\Sanna\Documents\Simulations\save_simulation\202201_NM2020_SF_simulation\plots\%s_projection%i_%s'%(date_saved,projection,itstr)

probe = np.load(openpath + '\\' + 'probe_%s.npy'%itstr) 
obj = np.load(openpath + '\\'  + 'object_%s.npy'%itstr, allow_pickle=True) 

x = np.squeeze(np.load(openpath + '\\' + 'x_%s.npy'%itstr)) 
y = np.squeeze(np.load(openpath + '\\' + 'y_%s.npy'%itstr)) 
z = np.squeeze(np.load(openpath + '\\' + 'z_%s.npy'%itstr)) 
#errors = np.load(openpath+ '\\errors.npy')
ferrors = np.load(openpath+ '\\ferrors.npy')


#TODO which coordinate, x,y,z?    # y och z är ungefär samma-= detector plane pix
#psize = x[1,0,0]-x[0,0,0]
psize = y[1,1,1]-y[0,0,0]
psizez = z[1,1,0]-z[0,0,0]

extent = 1e6 * np.array([0, (obj.shape[0])*psize, 0,(obj.shape[1])*psize]) # long axis, chort axis


if save == True:
    if not os.path.exists(savepath):
        os.makedirs(savepath)
        print('new folder in this savepath was created')


plt.close('all')
plot2drecons((np.rot90(obj,3)), probe, extent, savepath, save, signflip)

xx = np.linspace(extent[2],extent[3],obj.shape[0])

#plot_line_plot(np.rot90(obj,3),xx, row=162,save=False)
plt.figure()
#plt.plot(errors,'blue')
plt.plot(abs(ferrors),'red')
    
slicey = slice(137,190) #for simulated data recons 256 (should be the same, but for now)
slicex = slice(227,437)
#estimate_resolution_mtf(abs(np.rot90(obj,3)[slicey,slicex]),xx[slicex],row=24)
