"""
This module provides simulated 3D Bragg data.

SH: 2020-06-18 Written new class for loading data from NanoMAX in May2020
- normalization not working
- shifting not implemented

"""

import ptypy
from ptypy.core.data import PtyScan
import ptypy.utils as u
from ptypy import defaults_tree

from ptypy.experiment import register
from ptypy.utils.descriptor import EvalDescriptor
try:
	import hdf5plugin
except ImportError:
	print('Couldnt find hdf5plugin - better hope your h5py has bitshuffle!')
              
import h5py
import numpy as np
import os.path
import time

logger = u.verbose.logger

@register()
class NanomaxBraggJune2017(PtyScan):
    """
    Reads an early Bragg 3d ptycho format from Nanomax, multiple fly
    scans are run with rocking angle steps in between.

    Defaults:

    [name]
    default = NanomaxBraggJune2017
    type = str
    help = PtyScan subclass identifier

    [scans]
    default = []
    type = list
    help = List of scan numbers

    [theta_bragg]
    default = 0.0
    type = float
    help = Bragg angle

    [datapath]
    default = None
    type = str
    help = Path to folder containing the Sardana master file
    doc =

    [datafile]
    default = None
    type = str
    help = Sardana master file
    doc =

    [maskfile]
    default = None
    type = str
    help = Arbitrary mask file
    doc = Hdf5 file containing an array called 'mask' at the root level.

    [detfilepattern]
    default = None
    type = str
    help = Format string for detector image files
    doc = A format string with two integer fields, the first holds the scan number while the second holds the image number.

    """

    def __init__(self, pars=None, **kwargs):
        self.p = self.DEFAULT.copy(99)
        self.p.update(pars)
        self.p.update(kwargs)
        super(self.__class__, self).__init__(self.p)

    def load_common(self):
        """
        We have to communicate the number of rocking positions that the
        model should expect, otherwise it never knows when there is data
        for a complete POD. We also have to specify a single number for
        the rocking step size.
        """
        print('*** load_common')
        angles = []
        for scannr in self.p.scans:
            with h5py.File(self.p.datapath + self.p.datafile) as fp:
                angles.append(float(fp.get('entry%d'%scannr + '/measurement/gonphi').value))
        print(angles)
        step = np.mean(np.diff(sorted(angles)))
        print(step)
        return {
            'rocking_step': step,
            'n_rocking_positions': len(angles),
            'theta_bragg': self.p.theta_bragg,
            }

    def load_positions(self):
        """
        For the 3d Bragg model, load_positions returns N-by-4 positions,
        (angle, x, z, y). The angle can be relative or absolute, the
        model doesn't care, but it does have to be uniformly spaced for
        the analysis to make any sense.

        Let's load the positions (and images) in the order they were
        acquired: x fastest, then y, then scan number in the order
        provided.
        """
        print('*** load_positions')

        # first, calculate mean x and y positions for all scans, they
        # have to match anyway so may as well average them.
        x, y = [], []
        for scan in self.p.scans:
            with h5py.File(self.p.datapath + self.p.datafile) as fp:
                if scan == self.p.scans[0]:
                    # first pass: find out how many zeros to remove from the samx buffer
                    entry = 'entry%d' % scan
                    tmp = np.array(fp[entry + '/measurement/AdLinkAI_buff'])
                    for i in range(tmp.shape[1]):
                        if np.allclose(tmp[:, i:], 0.0):
                            cutoff = i
                            print('using %i samx values' % cutoff)
                            break
                x.append(np.array(fp[entry + '/measurement/AdLinkAI_buff'][:, :cutoff]))
                y.append(np.array(fp[entry + '/measurement/samy']))
        x_mean = -np.mean(x, axis=0) * 1e-6
        y_mean = np.mean(y, axis=0) * 1e-6
        Nx = x_mean.shape[1]
        Ny = x_mean.shape[0]
        Nxy = Nx * Ny
        assert Ny == y_mean.shape[0]
        print('Scan positions are Nx=%d, Ny=%d, Nxy=%d' % (Nx, Ny, Nxy))

        # save these numbers for the diff image loader
        self.Nx = Nx
        self.Ny = Ny
        self.Nxy = Nxy

        # then, go through the scans and fill in a N-by-4 array of positions per diff pattern
        x = np.empty(len(self.p.scans) * Nx * Ny, dtype=float)
        y = np.copy(x)
        theta = np.copy(x)
        z = np.zeros_like(x)
        for i in range(len(self.p.scans)):
            with h5py.File(self.p.datapath + self.p.datafile) as fp:
                entry = 'entry%d' % self.p.scans[i]
                th = fp[entry + '/measurement/gonphi'].value[0]
            x[i*Nxy:(i+1)*Nxy] = x_mean.flatten()
            y[i*Nxy:(i+1)*Nxy] = np.repeat(y_mean, Nx)
            theta[i*Nxy:(i+1)*Nxy] = th

        # adapt to our geometry
        tmp = z.copy()
        z = x   # our x motor goes toward positive z on the sample
        y = y   # our y motor goes toward more positive y on the sample
        x = tmp # this is the unimportant direction

        return np.vstack([theta, x, z, y]).T

    def load(self, indices):
        """
        This function returns diffraction image indexed from top left as
        viewed along the beam, i e they have (-q1, q2) indexing. PtyScan
        can always flip/rotate images.
        """
        print('*** load')
        raw, positions, weights = {}, {}, {}

        for ind in indices:
            scan = self.p.scans[ind // self.Nxy] # which scan to look in
            file = (ind % self.Nxy) // self.Nx   # which line of this scan
            frame = (ind % self.Nxy) % self.Nx   # which frame of this line
            with h5py.File(self.p.datapath + self.p.detfilepattern % (scan, file)) as fp:
                data = np.array(fp['entry_0000/measurement/Merlin/data'][frame])
            raw[ind] = data

        return raw, positions, weights

    def load_weight(self):
        print('*** load_weight')
        with h5py.File(self.p.maskfile) as fp:
            mask = np.array(fp['mask'])
        return mask

@register()
class NanomaxBraggMay2020(PtyScan):
    """
    This class loads data written with the nanomax pirate system,
    in a slightly matured state. Step and fly scan have the same
    format

    OR mybe I just stick to flyscans
    
    Defaults:

    [name]
    default = NanomaxBraggMay2020
    type = str
    help = PtyScan subclass identifier

    [scans]
    default = []
    type = list
    help = List of scan numbers

    [theta_bragg]
    default = 0.0
    type = float
    help = Bragg angle

    [path]
    default = None
    type = str
    help = Path to where the data is at
    doc =

    [xMotor]
    default = sx
    type = str
    help = Which x motor to use
    doc =

    [yMotor]
    default = sy
    type = str
    help = Which y motor to use
    doc =

    [xMotorFlipped]
    default = False
    type = bool
    help = Flip detector x positions
    doc =

    [yMotorFlipped]
    default = False
    type = bool
    help = Flip detector y positions
    doc =

    [xMotorAngle]
    default = 0.0
    type = float
    help = Angle of the motor x axis relative to the lab x axis
    doc =

    [yMotorAngle]
    default = 0.0
    type = float
    help = Angle of the motor y axis relative to the lab y axis
    doc =

    [maskfile]
    default = None
    type = str
    help = Arbitrary mask file
    doc = Hdf5 file containing an array called 'mask' at the root level.

    [I0]
    default = None
    type = str
    help = Normalization channel, like ni/counter1 for example
    doc =



    """
    def __init__(self, pars=None, **kwargs):
        self.p = self.DEFAULT.copy(99)
        self.p.update(pars)
        self.p.update(kwargs)
        super(self.__class__, self).__init__(self.p)

    def load_common(self):
        """
        We have to communicate the number of rocking positions that the
        model should expect, otherwise it never knows when there is data
        for a complete POD. We also have to specify a single number for
        the rocking step size.
        """
        print('*** load_common')

        angles = []
        for scan in self.p.scans:
            filename = '%06u.h5' % scan
            fullfilename = os.path.join(self.info.path, filename)
            with h5py.File(fullfilename, 'r') as fp:
                angles.append(float(fp.get('entry/snapshot/gonphi')[0]))
        print(angles)
        step = np.mean(np.diff(sorted(angles)))
        print(step)
        return {
            'rocking_step': step,
            'n_rocking_positions': len(angles),
            'theta_bragg': self.p.theta_bragg,
            }
    
    def load_positions(self):
        """
        For the 3d Bragg model, load_positions returns N-by-4 positions,
        (angle, x, z, y). The angle can be relative or absolute, the
        model doesn't care, but it does have to be uniformly spaced for
        the analysis to make any sense.

        Let's load the positions (and images) in the order they were
        acquired: x fastest, then y, then scan number in the order
        provided.
        """
        print('*** load_positions')

        #TODO: I wait a bit with the norm data
        #TODO write it so that it opens both step and fly scans. later.
        xFlipper, yFlipper = 1, 1
        if self.info.xMotorFlipped:
            xFlipper = -1
            logger.warning("note: x motor is specified as flipped")
        if self.info.yMotorFlipped:
            yFlipper = -1
            logger.warning("note: y motor is specified as flipped")

        # if the x axis is tilted, take that into account.
        xCosFactor = np.cos(self.info.xMotorAngle / 180.0 * np.pi)
        yCosFactor = np.cos(self.info.yMotorAngle / 180.0 * np.pi)
        logger.info(
            "x and y motor angles result in multiplication by %.2f, %.2f" % (xCosFactor, yCosFactor))

        # first, calculate mean x and y positions for all scans, they
        # have to match anyway so may as well average them.
        normdata, x, y = [], [], []
        for scan in self.p.scans:
            filename = '%06u.h5' % scan
            fullfilename = os.path.join(self.info.path, filename)

            #interlude to find Nx and Ny (from the first scan)
            if scan == self.p.scans[0]:
                with h5py.File(fullfilename, 'r') as fp:
                    #TODO: this maybe doesnt always exist with flyscans
                    Ny = len(np.array(fp['entry/measurement/sy']))
                    Nxy = len(np.array(fp['entry/measurement/%s' % (self.info.xMotor)]))
                Nx = int(Nxy/Ny)
        
##                # may as well get normalization data here too
##                if self.info.I0 is not None:
##                    with h5py.File(fullfilename, 'r') as hf:
##                        normdata.append(np.array(hf['entry/measurement/%s' % (self.info.I0)], dtype=float))
##                    print('*** going to normalize by channel %s' % self.info.I0)

            #Saves postions x and y as a Nscans x Ny x Nx matrix
            with h5py.File(fullfilename, 'r') as hf:
                # Just a long list of Nxy values. Save as Ny x Nx array instead
                x.append( xFlipper * xCosFactor
                     * np.array(hf['entry/measurement/%s' % (self.info.xMotor)]).reshape(Ny,Nx)  )
                y.append(yFlipper * yCosFactor
                     * np.array(hf['entry/measurement/%s' % (self.info.yMotor)]).reshape(Ny,Nx)  )

        #yes you should use the 1e-6 factor
        # Average over the 1st axis which is scans
        x_mean = np.mean(x, axis=0) * 1e-6
        y_mean = np.mean(y, axis=0) * 1e-6

        assert Ny == y_mean.shape[0]
        print('Scan positions are Nx=%d, Ny=%d, Nxy=%d' % (Nx, Ny, Nxy))

        # Save to use in load class
        self.Nxy = Nxy

        # then, go through the scans and fill in a N-by-4 array of positions per diff pattern
        x = np.empty(len(self.p.scans) * Nx * Ny, dtype=float)
        y = np.copy(x)
        theta = np.copy(x)
        z = np.zeros_like(x)
        for i in range(len(self.p.scans)):
            filename = '%06u.h5' % self.p.scans[i]
            fullfilename = os.path.join(self.info.path, filename)
            with h5py.File(fullfilename, 'r') as fp:
                th = fp['entry/snapshot/gonphi'][0]
            x[i*Nxy:(i+1)*Nxy] = x_mean.flatten()
            y[i*Nxy:(i+1)*Nxy] = y_mean.flatten()
            theta[i*Nxy:(i+1)*Nxy] = th

        # adapt to our geometry
        tmp = z.copy()
        #TODO check this
        z = x   # our x motor goes toward positive z on the sample
        y = y   # our y motor goes toward more positive y on the sample
        x = tmp # this is the unimportant direction

        #TODO in 2d script they do 'positions = -np.vstack((y, x)).T * 1e-6'

        return np.vstack([theta, x, z, y]).T


    def load(self, indices):
        """
        This function returns diffraction image indexed from top left as
        viewed along the beam, i e they have (-q1, q2) indexing. PtyScan
        can always flip/rotate images.
        """
        print('*** load')
        raw, positions, weights = {}, {}, {}

        for ind in indices:
            scan = self.p.scans[ind // self.Nxy]    # which scan to look in
            filename = '%06u.h5' % scan
            fullfilename = os.path.join(self.info.path, filename)
            frame = (ind % self.Nxy)                # which frame to save
            with h5py.File(fullfilename, 'r') as fp:
                raw[ind] = np.array(fp['entry/measurement/merlin/frames'][frame])
                #TODO add normlisation.
                # think add also normalisation between rotations? OR is that what im doing?
                #if self.info.I0:
                    #    raw[ind] = raw[ind] / self.normdata[ind]

        return raw, positions, weights

    def load_weight(self):
        """
        Loads a mask from a h5 file
        """
        print('*** load_weight')
        with h5py.File(self.p.maskfile, 'r') as fp:
            mask = np.array(fp.get('mask'))
            logger.info("loaded mask, %u x %u, sum %u, so %u masked pixels" %
                        (mask.shape + (np.sum(mask), np.prod(mask.shape)-np.sum(mask))))
        return mask

    

