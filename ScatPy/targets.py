# -*- coding: utf-8 -*-
"""
Target definitions.

"""

from __future__ import division
import subprocess
import numpy as np
import os
import os.path
import copy
import results

import utils
from config import exec_settings

#default spacing of btw dipoles (in um)
default_d=0.015 


class Target(object):
    """
    The base class for defining DDscat Target geometries.
    
    """


    def __init__(self, phys_shape, d=None, material=None, folder=None):
        '''
        Initialize the Target

        Arguments:
            phys_shape: The physical shape of the object (in microns).
            d: the dipole density
            material: the material file to use for the target (currently limited to one)
                      default is 'Au_Palik.txt'
            folder: the target working directory
        '''


        self.directive=''
        self.NCOMP=1
        if folder is None:
            self._folder='.'
        else:
            self._folder=folder
            
        if material is None:
            self.mat_file='Au_Palik.txt'
        else:
            self.mat_file=material

        if d is None:    
            self.d=default_d #dipole grid spacing in um 
        else:
            self.d=d

        self.shape=np.int32(np.ceil(np.asarray(phys_shape)/self.d))
            
        #self.aeff=0.246186 #How_Range(0.246186, 0.246186, 1, 'LIN')# = aeff (first,last,how many,how=LIN,INV,LOG)
        self.N=0

    @property
    def aeff(self):
        """Calculate the effective diameter of the target"""
                    
        return (self.N*3/4/np.pi)**(1/3)*self.d

    @property    
    def phys_shape(self):
        '''
        Returns the size of the bounding box that contains the target (in microns)
    
        '''
        return self.shape*self.d
        
    def save_str(self):
        """Return the four line string of target definition for inclusion in the ddscat.par file"""
        out='**** Target Geometry and Composition ****\n'
        out+=self.directive+'\n'
        out+=str(self.shape)[1:-1]+'\n'
        out+=str(self.NCOMP)+'\n'
        out+='\''+utils.resolve_mat_file(self.mat_file)+'\'\n'
        
        return out

    
    def write(self):
        pass

    def VTRconvert(self, outfile=None):
        """Execute VTRconvert to generate a file for viewing in Paraview"""


        if outfile is None:
            outfile='output'
            
        self.write()
            
        if os.path.exists(os.path.join(self.folder, self.fname)):
            #This assumes that vtrconvert is in your path
            #check os.environ['PATH'] to be sure
            subprocess.call(['vtrconvert', self.fname, outfile], cwd=self.folder)
        else:
            print 'No target.out file to convert'

    @property
    def folder(self):
        """The target working directory"""
        return self._folder
        
    @folder.setter
    def folder(self, newfolder):
        """Redefine the target working directory"""
        self._folder=newfolder
    
    def copy(self):
        return copy.deepcopy(self)

class Target_Builtin(Target):
    """Base class for target geometries that are built into DDSCAT"""
    def __init__(self, *args, **kwargs):
        Target.__init__(self, *args, **kwargs)
    

class RCTGLPRISM(Target_Builtin):        
    """A rectangular prism target"""

    def __init__(self, phy_shape, **kwargs):
        """
        Initialize the rectangular prism.
        
        Arguments:
        phys_shape: 3-tuple defining the the physical shape of the prism in microns
        
        **kwargs are passed to Target
        """
        Target_Builtin.__init__(self, phy_shape, **kwargs)
        self.directive='RCTGLPRSM'
        self.N=phy_shape[0]*phy_shape[1]*phy_shape[2]
        #self.d=0.015624983


class CYLNDRCAP(Target):
    """A target cylinder with hemispherical endcaps"""
    
    def __init__(self, length, radius, **kwargs):
        """
        Initialize the capped cylinder.
        
        Arguments:
            length: the length of the cylinder in microns (not including endcaps)
            rad: the radius of the cylinder
            
            **kwargs are passed to Target

        Total height of the structureis length+2*rad
        """
        Target_Builtin.__init__(self, np.asarray([length, radius*2, radius*2]), **kwargs)
        
        self.directive='CYLNDRCAP'
        
        Vcyl=self.shape[0]*(np.pi*(self.shape[1]/2)**2)
        Vsph=4/3*np.pi*(self.shape[1]/2)**3
        self.N=int(Vcyl+Vsph)
        self.length=length
        self.radius=radius
        #self.d=0.015624983
        

class ELLIPSOID(Target):        
    """
    An Ellipsoid target
    """
    
    def __init__(self, phys_shape, **kwargs):
        """
        Initialize the ellipsoid target.
        
        Arguments:
            phys_shape: 3-tuple giving the lengths of the three semi-axes

            **kwargs are passed to Target
        """        
        
        Target_Builtin.__init__(self, np.asarray(phys_shape)*2, **kwargs)
        self.directive='ELLIPSOID'
        self.N=int(4/3*np.pi*(self.shape.prod()/8))

class CYLINDER(Target):
    """A target cylinder with hemispherical endcaps"""
    
    def __init__(self, length, radius, ori,  **kwargs):
        """
        Initialize the cylinder.
        
        Arguments:
            length: the length of the cylinder in microns (not including endcaps)
            rad: the radius of the cylinder
            
            **kwargs are passed to Target

        """
        
        Target_Builtin.__init__(self, np.asarray([length, radius*2, ori]), **kwargs)
        
        self.directive='CYLINDER1'
        self.ori=ori
        
        Vcyl=self.shape[0]*(np.pi*(self.shape[1]/2)**2)
        self.N=int(Vcyl)
        self.length=length
        self.radius=radius

    def save_str(self):
        """Return the four line string of target definition for inclusion in the ddscat.par file"""
        out='**** Target Geometry and Composition ****\n'
        out+=self.directive+'\n'
        out+=str(self.shape[0:2])[1:-1]+' '+str(self.ori)+'\n'
        out+=str(self.NCOMP)+'\n'
        out+='\''+utils.resolve_mat_file(self.mat_file)+'\'\n'
        
        return out
        
class Sphere(ELLIPSOID):  
    """
    A Sphere target.

    """      
    def __init__(self, radius, **kwargs):
        """
        Initialize the sphere target.
        
        Arguments:
            radius: the radius of the sphere

            **kwargs are passed to Target
            
        """
        phys_shape=[radius]*3
        ELLIPSOID.__init__(self, phys_shape, **kwargs)
            

class IsoHomo_FROM_FILE(Target):
    '''
    Base class for an arbitrary, isotropic and homogeneous target.
    
    '''
    def __init__(self, shape, **kwargs):
        """
        Initialize a new target
        
        shape is in units of number of dipoles
        """
        Target.__init__(self, shape, **kwargs)
        self.descrip=''
        self.directive='FROM_FILE'
        self.fname='shape.dat'
        self.descriptor='FROM_FILE'
    
        self.shape=np.asarray(shape)
        self.grid=np.zeros(shape, dtype=np.bool)    
        self.update_N()
        self.a1=np.array([1,0,0])
        self.a2=np.array([0,1,0])
        self.rel_d=np.array([1,1,1])
        self.origin=self.shape/2           
    
    def update_N(self):
        """Update the number of dipoles"""
        self.N=self.grid.sum()
    
    def write(self):
        """Write the shape file."""
        self.update_N()
        with open(os.path.join(self.folder, self.fname), 'wb') as f:
            f.write(self.descriptor+'\n') 
            f.write(str(self.N)+'\n')
            f.write(str(self.a1)[1:-1]+'\n')
            f.write(str(self.a2)[1:-1]+'\n')
            f.write(str(self.rel_d)[1:-1]+'\n')
            f.write(str(self.origin)[1:-1]+'\n')
            f.write('J JX JY JZ ICOMPX,ICOMPY,ICOMPZ'+'\n')
            
            n=1
            for (ni, i) in enumerate(self.grid):
                for (nj, j) in enumerate(i):
                    for (nk, k) in enumerate(j):
                        if k:
                            f.write('%d %d  %d  %d  %d  %d  %d\n' % (n, ni, nj, nk, 1, 1, 1))
                            n+=1
    
    def VTRconvert(self, outfile=None):
        """Execute VTRConvert to generate a model file viewable in Paraview"""
        Target.VTRconvert(self, outfile)

    def show(self, *args, **kwargs):
        """
        Display the shape using mayavi
        
        """
        self.write()
        fname=os.path.join(self.folder, self.fname)
        s=results.ShapeTable(fname)
        s.show(*args, **kwargs)
        
class FROM_FILE(Target):
    pass        

class Iso_FROM_FILE(FROM_FILE):
    pass        


class Ellipsoid_FF(IsoHomo_FROM_FILE):
    """
    Build an ellipsoidal target to be loaded from file
    """

    def __init__(self, phys_shape, d=None, **kwargs):
        """
        Create a new Ellipsoid Target
        
        phys_shape is the length of the three semi-major axes in physical units
        """
        if d is None:
            d=default_d
        shape=np.int16(np.array(phys_shape) *2/d)
        IsoHomo_FROM_FILE.__init__(self, shape, **kwargs)

        #self.phys_shape=phys_shape

        self.descriptor='Ellipsoid_FF (%f, %f, %f, %f)'%(tuple(self.phys_shape)+(self.d,))            
        
        (a,b,c) = tuple(shape)
        xx,yy,zz=np.mgrid[-a:a, -b:b, -c:c] 
        dist=(xx/a)**2 + (yy/b)**2 + (zz/c)**2

        self.grid = dist<1
        
        self.N=self.grid.sum()

class Helix(IsoHomo_FROM_FILE):
    """
    A helix target
    
    """
    def __init__(self, height, pitch, major_r, minor_r, d=None, build=True, **kwargs):
        '''
        Build a Helix Target
        
        dimensions are physical, in um

        Arguments:
            height: the height of the helix (not counting it's thickness)
            pitch: the helix pitch, um/turn
            major_r: the radius of the helix sweep
            minor_r: the radius of the wire that is swept to form the helix
        
        Optional Arguments:
            d: the dipole dipole spacing, if not specified default value is used
            build: if False, delays building the helix until requested. default is True
            
            **kwargs are passed to Target
        '''

        if d is None:
            d=default_d
        shape=np.int16(np.asarray([height+2*minor_r, 2*(major_r+minor_r), 2*(major_r+minor_r)])/d)
        shape+=1
        IsoHomo_FROM_FILE.__init__(self, shape, **kwargs)

        self.height=height
        self.pitch=pitch
        self.major_r=major_r
        self.minor_r=minor_r
        self.d=d

        if build:
            self.build_helix()
        else:
            self.descriptor='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            

    def build_helix(self):

        self.descriptor='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            
        
        print 'Generating Helix...'
        #the helix dimensions in pixels
        p_height=self.height/self.d
        p_pitch=-self.pitch/self.d
        p_major_r=self.major_r/self.d
        p_minor_r=self.minor_r/self.d        

        #the sweep path
        t=np.linspace(0,1, self.height/abs(self.pitch)*360)
#        x=np.int16(p_height*t + p_minor_r)
#        y=np.int16(p_major_r * np.cos(2*np.pi* x/p_pitch) + self.origin[1])
#        z=np.int16(p_major_r * np.sin(2*np.pi* x/p_pitch) + self.origin[2])
#        
#        p=np.vstack([x,y,z]).transpose()
#        p=np.vstack([np.array(u) for u in set([tuple(l) for l in p])]) #remove duplicates

        x=p_height*t + p_minor_r
        y=p_major_r * np.cos(2*np.pi* x/p_pitch) + self.origin[1]
        z=p_major_r * np.sin(2*np.pi* x/p_pitch) + self.origin[2]
        
        p=np.vstack([x,y,z]).transpose()
#        p=np.vstack([np.array(u) for u in set([tuple(l) for l in p])]) #remove duplicates
        print 'Done constructing sweep path...'
        def dist(v):
            return np.sqrt(v[:,0]**2 +v[:,1]**2 + v[:,2]**2)
            
        def sphere_check(i, j, k):
            
            d=dist(np.array([i,j,k])-p)
            if np.any(d<=p_minor_r):
                return True
            else:
                return False
                
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for k in range(self.shape[2]):
                    self.grid[i,j,k]=sphere_check(i,j,k)

        self.update_N()

class SpheresHelix(IsoHomo_FROM_FILE):
    """
    A helix target composed of isolated spheres
    
    """
    def __init__(self, height, pitch, major_r, minor_r, n_sphere, d=None, build=True, **kwargs):
        '''
        Build a Sphere Helix Target
        
        dimensions are physical, in um

        Arguments:
            height: the height of the helix (not counting it's thickness)
            pitch: the helix pitch, um/turn
            major_r: the radius of the helix sweep
            minor_r: the radius of the spheres that compose the helix
            n_sphere: the number of spheres that compose the helix
        
        Optional Arguments:
            d: the dipole dipole spacing, if not specified default value is used
            build: if False, delays building the helix until requested. default is True
            
            **kwargs are passed to Target
        '''

        if d is None:
            d=default_d
        shape=np.int16(np.asarray([height+2*minor_r, 2*(major_r+minor_r), 2*(major_r+minor_r)])/d)
        shape+=1
        IsoHomo_FROM_FILE.__init__(self, shape, **kwargs)

        self.height=height
        self.pitch=pitch
        self.major_r=major_r
        self.minor_r=minor_r
        self.n_sphere=n_sphere
        self.d=d

        if build:
            self.build_helix()
        else:
            self.descriptor='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            

    def build_helix(self):

        self.descriptor='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            
        
        print 'Generating Helix...'
        #the helix dimensions in pixels
        p_height=self.height/self.d
        p_pitch=-self.pitch/self.d
        p_major_r=self.major_r/self.d
        p_minor_r=self.minor_r/self.d      

        #the sweep path
        t=np.linspace(0,1, self.n_sphere)
        x=p_height*t + p_minor_r
        y=p_major_r * np.cos(2*np.pi* p_height/p_pitch*t) + self.origin[1]
        z=p_major_r * np.sin(2*np.pi* p_height/p_pitch*t) + self.origin[2]
        
        p=np.vstack([x,y,z]).transpose()
#        p=np.vstack([np.array(u) for u in set([tuple(l) for l in p])]) #remove duplicates
        print 'Done constructing sweep path...'
        def dist(v):
            return np.sqrt(v[:,0]**2 +v[:,1]**2 + v[:,2]**2)
            
        def sphere_check(i, j, k):
            
            d=dist(np.array([i,j,k])-p)
            if np.any(d<=p_minor_r):
                return True
            else:
                return False
                
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for k in range(self.shape[2]):
                    self.grid[i,j,k]=sphere_check(i,j,k)

        self.update_N()


class Conical_Helix(IsoHomo_FROM_FILE):
    """
    A helix target
    
    """
    def __init__(self, height, pitch1, pitch2, major_r, minor_r, d=None, build=True, **kwargs):
        '''
        Build a Conical profile Helix Target
        
        dimensions are physical, in um

        Arguments:
            height: the height of the helix (not counting it's thickness)
            pitch1,2: the helix pitch, um/turn (if pitch1 is different from pitch2 will swipe between during the growth)
                        major_r: the radius of the helix sweep
            minor_r: the radius of the wire that is swept to form the helix
        
        Optional Arguments:
            d: the dipole dipole spacing, if not specified default value is used
            build: if False, delays building the helix until requested. default is True
            
            **kwargs are passed to Target
        '''

        if d is None:
            d=default_d
        shape=np.int16(np.asarray([height+2*minor_r, 2*(major_r+minor_r), 2*(major_r+minor_r)])/d)
        shape+=1
        IsoHomo_FROM_FILE.__init__(self, shape, **kwargs)

        self.height=height
        self.pitch1=pitch1
        self.pitch2=pitch2
        self.major_r=major_r
        self.minor_r=minor_r
        self.d=d

        if build:
            self.build_helix()
        else:
            self.descriptor='FROM_FILE_Helix (%f, %f, %f, %f, %f, %f)'%(self.height, self.pitch1, self.pitch2, self.major_r, self.minor_r, self.d)            

    def build_helix(self):

        self.descriptor='FROM_FILE_Helix (%f, %f,%f, %f, %f, %f)'%(self.height, self.pitch1, self.pitch2, self.major_r, self.minor_r, self.d)            
        
        print 'Generating Helix...'
        #the helix dimensions in pixels
        p_height=(self.height*(1-(self.pitch2-self.pitch1)/self.pitch2))/self.d
        p_pitch1=-self.pitch1/self.d
        p_pitch2=-self.pitch2/self.d
        p_major_r=self.major_r/self.d
        p_minor_r=self.minor_r/self.d        

        #the sweep path
        
        t=np.linspace(0,1, self.height/abs(self.pitch1)*360)
        k=p_height*(1+t*(p_pitch2-p_pitch1)/p_pitch2)
       
        P=t*k
        R= t*p_major_r
        
        
        
        x=P + p_minor_r
        phi=2*np.pi* t*p_height/p_pitch1
        y=R* np.cos(phi) + self.origin[1]
        z=R* np.sin(phi) + self.origin[2]


       
        p=np.vstack([x,y,z]).transpose()
        
        print 'Done constructing sweep path...'
        def dist(v):
            return np.sqrt(v[:,0]**2 +v[:,1]**2 + v[:,2]**2)
            
        def sphere_check(i, j, k):
            
            d=dist(np.array([i,j,k])-p)
            if np.any(d<=p_minor_r):
                return True
            else:
                return False
                
        for i in range(self.shape[0]):
            for j in range(self.shape[1]):
                for k in range(self.shape[2]):
                    self.grid[i,j,k]=sphere_check(i,j,k)
        
        self.update_N()

def Holify(target, radius, posns=None, num=None, seed=None):
    """
    Punches holes into a target
    
    radius: the radius of the holes to be punched (in units of d)

    Where to punch the holes can be specified by:    
    
    posns: a list of the x,y,z positions to punch the holes
    
    or    
    
    num: the number of holes. their positions are selected randomly
    seed: a seed value to use for the random number generator
    
    """
    shape=target.grid.shape
    
    if posns is None:
        if num is None:
            raise ValueError('Either posns or num must be specificed')
        
        posns=np.random.rand(num, 3)
        for (p, s) in zip(posns.T, shape):
            p *= s

    posns=np.asarray(posns, dtype=np.int16)

    r=np.ceil(radius)
    xx,yy,zz=np.mgrid[-(r+1):r+1, -(r+1):r+1, -(r+1):r+1] 
    dist=np.sqrt(xx**2 + yy**2 + zz**2)
    mask=np.array(np.where(dist<r))
    
    new_target=target.copy()
    grid=np.ravel(new_target.grid)
    print grid.sum()
    
    for p in posns:
        p=p[:, np.newaxis]
#        print p            
        mask_r=np.ravel_multi_index(mask+p, shape, mode='clip')
        print mask_r[0]
#        g_len=grid.sum()
        print grid[mask_r[0]]
        grid[mask_r] &= False
        print grid[mask_r[0]]

#        print g_len, '->', grid.sum()

    print grid.sum()
    new_target.N = new_target.grid.sum()
    return new_target