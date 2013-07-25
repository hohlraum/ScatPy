# -*- coding: utf-8 -*-
"""
A set of tools for setting up and working with the program DDScat.


"""

from __future__ import division
import subprocess
import warnings
import numpy as np
import os
import os.path
import posixpath
import copy

import targets
import results
import fileio
import ranges

import config

#: DEFINE POLARIZATION STATES USING SPECTROSCOPIST'S CONVENTION
pol_cR=np.array([0, 0+1j, 1+0j])    
pol_cL=np.array([0, 1+0j, 0+1j])            
pol_lH=np.array([0, 0+0j, 1+0j])        
pol_lV=np.array([0, 1+0j, 0+0j])    

class Settings():
    '''
    DDSCAT execution parameters
    
    Most of the field names correspond to their definitions in the ddscat.par file.
    Details of the target are stored in the ```Target``` definition, not here.
   
    '''

    def __init__(self, **kwargs):        
        #: Either do or skip torque calculations (True:'DOTORQ', False:'NOTORQ')
        self.CMDTRQ= False 

        #: Solution method (PBCGS2, PBCGST, GPBICG, PETRKP, QMRCCG) 
        self.CMDSOL='PBCGS2' #: = CMDSOL*6 (PBCGS2, PBCGST, GPBICG, PETRKP, QMRCCG) -- 

        #: FFT method (GPFAFT, FFTMKL)
        self.CMDFFT='GPFAFT'

        #: Prescription for polarizabilities (GKDLDR, LATTDR)
        self.CALPHA='GKDLDR'

        #: Specify binary output
        self.CBINFLAG='NOTBIN'

        #: Initial Memory Allocation
        self.InitialMalloc=np.array([100,100,100])

        #: Either do or skip nearfield calculations (True, False)
        self.NRFLD=False 
        
        #: Fractional extension of calculated volume in (-x,+x,-y,+y,-z,+z)
        self.NRFLD_EXT=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) 

        #: Error Tolerence. MAX ALLOWED (NORM OF |G>=AC|E>-ACA|X>)/(NORM OF AC|E>)    
        self.TOL=1.00e-5
        
        #: Maximum number of iterations
        self.MXITER=600  
        
        #: Interaction cutoff parameter for PBC calculations (1e-2 is normal, 3e-3 for greater accuracy)
        self.GAMMA=1.00e-2 

        #: Angular resolution for calculation of <cos>, etc. (number of angles is proportional to [(3+x)/ETASCA]^2 )
        self.ETASCA=0.5
        
        #: Specify which output files to write (True: write ".sca" file for each target orient. False: suppress)
        self.IWRKSC=True 
        
        #: Vacuum wavelengths (first,last,how many,how=LIN,INV,LOG)
        self.wavelengths=ranges.How_Range(0.3500, 0.8000, 10, 'LIN') 

        #: Refractive index of ambient medium
        self.NAMBIENT=1.000
        
        #: Define a range of scales for the particle geometry, None indicates a single size calc
        self.scale_range=None 

        #: Define Incident Polarizations (Polarization state e01 (k along x axis)
        self.Epol=pol_lV
        
        #: Specify whether to calculate orthogonal polarization state (True, False)
        self.IORTH=2

        #: Prescribe Target Rotation beta (rotation around a1). (betamin, betamx, nbeta)
        self.beta=ranges.Lin_Range(0.,0.,1)

        #: Prescribe Target Rotation theta (angle between a1 and k). (thetamin, tetamx, ntheta)
        self.theta=ranges.Lin_Range(0.,0.,1)

        #: Prescribe Target Rotation phi (rotation of a1 around k). (phimin, phimax, nphi)
        self.phi=ranges.Lin_Range(0.,0.,1)

        #: Specify first IWAV, IRAD, IORI (0 0 0 to begin fresh)
        self.first_I=[0,   0,   0]

        #: Select Elements of S_ij Matrix to Print
        self.S_INDICES=[11, 12, 13, 14, 21, 22, 31, 41, 44]

        #: Specify reference frame for scattering ('LFRAME', 'TFRAME')
        self.CMDFRM='LFRAME'

        #: Specify Scattered Directions
        self.scat_planes=[ranges.Scat_Range(0,0,180,5), ranges.Scat_Range(90,0,180,5)]        

    def copy(self):
        """
        Make a copy of these settings.
        """
        return copy.deepcopy(self)
        

class DDscat(object):
    """
    A class for managing a DDSCAT run.
    
    Loosely a DDscat object includes two parts: a settings file and a target.
    
    All parameters are optional.
    :param folder: The subfolder in which to store files.
    :param settings: The :class:`Settings` object for this run.
    :param target: The :class:`targets.Target` for this run. Defaults to a Au sphere.
    
    
    Example::

        # Build a target
        t=CYLNDRCAP(0.100, 0.030):        
        
        # Initialize a DDscat run
        d=DDscat(target=t)
    
        # Modify run parameters as desired
        d.settings.NAMBIEND=1.5
        d.settings.Epol=np.array([0, 1j, 1])

        # Write the output and run the simulation
        d.calculate()

        # Load the results
        r=results.FolderCollection()
    
        # Plot them
        r.plot()
    """

    def __init__(self, folder=None, settings=None, target=None):
        
        if folder is None:
            self._folder='.'
        else:
            if not os.path.exists(folder):
                os.makedirs(folder)
            self._folder=folder
        
        if settings is None:
            self.settings=Settings(folder=self._folder)
        else:
            self.settings=settings.copy()
            self.settings.folder=self._folder
        
        if target is None:
            self._target=targets.Sphere(0.2, folder=self._folder)
        else:
            self._target=target
            self.target.folder=self._folder
            
        #self.results=results.Results(folder=self.folder)
        
    def __str__(self):
        """
        A string representation.
        This should be expanded
        """
        return 'DDScat Definition'

    def copy(self):
        """
        Make a copy of this run.
        """
        return copy.deepcopy(self)

    def write_submit_script(self):
        """
        Write a script for submitting the job to SGE via qsub.
        """        
            
                        
    def write(self, *args, **kwargs):
        """Write the .par file and target definitions to file.
        
        :param args: Optional arguments and kwargs are passed to the profile_script

        :param submit_script: Use True to write a .sge file for submitting the job 
                            to an SGE cluster via ```qsub```.
        
        If the current profile has a ```write_script``` then it is run once
        the ddscat.par has been written, and is called with any addition
        arguments to write.        
        """

        s=fileio.build_ddscat_par(self.settings, self.target)
        
        with open(os.path.join(self.folder, 'ddscat.par'), 'wb') as f:
            f.write(s)

        self.target.folder=self.folder
        self.target.write()

        try:
            config.config['write_script'](self, *args, **kwargs)
        except KeyError:
            pass            

    @property
    def folder(self):
        """This run's home folder"""
        return self._folder
    
    @folder.setter    
    def folder(self, newfolder):
        """Redefine the run's home folder."""
        if not os.path.exists(newfolder):
            os.makedirs(newfolder)
      
        self._folder=newfolder
        self.target.folder = newfolder

    @property
    def target(self):
        """The run's target"""
        return self._target

    @target.setter        
    def target(self, newtarget):
        """Redefine the run's target."""
        self._target=newtarget
        self._target.folder = self.folder

    def info(self):
        """Print some basic run info"""
        wave=self.settings.wavelengths.table
        N=self.target.N
        
        alpha=self.alpha()
        beta=self.beta()
        x=self.x()
        table=np.vstack([wave, x, alpha, beta])
        print "Target: "+self.target.directive
        print 'Shape: ' + str(self.target.shape)
        print 'N=%d'%N
        print 'aeff=%f'%self.target.aeff()
        print 'd=%f'%self.target.d

        print 'wave\tx\talpha\tbeta'
        print '-----------------------------------'
        for r in table.transpose():
            print r

    @property    
    def x(self):
        """Calculate the x-parameter (Userguide p8)."""
        a=self.target.aeff
        out=[]        
        for l in self.settings.wavelengths:
            out.append(2*np.pi*a/l)

        return np.asarray(out)
    
    @property    
    def mkd(self):
        """Calculate m*k*d (Userguide p8)."""
        m_dat=results.MInTable(self.target.mat_file)
        out=[]
        for l in self.settings.wavelengths:
            m=m_dat(l)
            out.append(m*2*np.pi/l*self.target.d)
            
        return np.asarray(out)

    @property        
    def alpha(self):
        """Calculate the alpha parameter (Userguide Eqn 8, p8)."""
        N=self.target.N
        m_dat=results.MInTable(self.target.mat_file)
        out=[]
        for l in self.settings.wavelengths:
            m=m_dat(l)
            out.append(9.88*l/np.abs(m['Rem']+1j*m['Imm'])*(N/10**6)**(1/3))

        return np.asarray(out)

    @property    
    def beta(self):
        """Calculateth beta parameter (Userguide Eqn 8, p8)."""
        N=self.target.N
        m_dat=results.MInTable(self.target.mat_file)
        out=[]
        for l in self.settings.wavelengths:
            m=m_dat(l)
            out.append(62/np.abs(m['Rem']+1j*m['Imm'])*(N/10**6)**(1/3))

        return np.asarray(out)
            

    def calculate(self, silent=False):
        """Start local calculation of this run.
        
        This assumes that ```ddscat``` is in the path.
        
        :param silent: If true suppresses output to the screen
        """

        self.write()

        command=os.path.join(config.config['ddscat_path'], 'ddscat')

        try:
            __IPYTHON__
        except NameError:
            pass
        else:
            silent=True
            warnings.warn('DDscat output will not be displayed in IPython')

        if silent:
            print 'Starting calculation...'
            subprocess.call(command+' 2> output.log', shell=True, cwd=self.folder)
            print 'Done!'
        else:
            subprocess.call(command + ' 2>&1 | tee output.log', shell=True, cwd=self.folder)

    def VTRconvert(self, outfile=None):
        """
        Convert the target shape into a form viewable in Paraview

        """
        self.target.VTRconvert()

    def calltarget(self):
        '''
        Executes ```calltarget``` to generate a target.ut file for builtin target geometries.
        '''
            
        self.write()

        subprocess.call(os.path.join(config.config['ddscat_path'], 'calltarget'), cwd=self.folder)



    
    
    
