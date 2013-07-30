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
import warnings
import pdb

import utils
import fileio
import ranges

#: Default spacing between dipoles (in um)
default_d=0.015 

class Target(object):
    """
    The base class for defining DDscat Target geometries.
    '''
    Initialize the Target

    :param directive: The type of target (e.g. CYLINDER1, FROM_FILE)
    :param sh_param: The three shape parameters used by DDSCAT to describe the target
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.

    As it is, this class only creates generic targets, which contain only
    a directive, sh_param, material and aeff: the bare minimum to be useful
    necessary to write a ddscat.par file.

    Typically, this class will be subclassed to create more useful and feature-
    rich targets. Derived classes *must* provide the following attributes and
    methods:
    * sh_param(): A property that returns the three values of the SHPAR definition used by DDSCAT 
    * _calc_N(): A static method to calculate the number of dipoles based on the shape parameters
    * fromfile(): A class method to generate a working object from a ddscat.par file

    """

    def __init__(self, directive=None, sh_param=(0,0,0), material=None, aeff=None, folder=None):

        if (directive is None) or (material is None):
            # If available, settings come from default.par file       
            default = utils.resolve_profile('default.par')
            if default is not None: 
                vals = self._read_values(default)
                if directive is None:
                    directive= vals['directive']
                    sh_param = vals['sh_param']
                    material = vals['material']
                    self.aeff = aeff
                elif material is None:
                    material = vals['material']
            
        self.directive = directive
        self.sh_param = sh_param
        self.material=list(material)

        if aeff is not None:
            if isinstance(aeff, ranges.How_Range):
                self.aeff = aeff.first
            else:
                self.aeff = aeff
                
        if folder is None:
            self._folder='.'
        else:
            self._folder=folder
            
    def save_str(self):
        """Return the four line string of target definition for the ddscat.par file"""
        out='**** Target Geometry and Composition ****\n'
        out+=self.directive+'\n'
        out+=str(self.sh_param)[1:-1]+'\n'
        out+=str(len(self.material))+'\n'
        for mat in self.material:
            out+='\''+utils.resolve_mat_file(mat)+'\'\n'
        
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

    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified file.
        
        If the target directive matches the name of an existing target class 
        then object creation is delegated to that class.
        """
        t_vals = cls._read_values(fname)

        directive = t_vals['directive']
        if directive in globals():
            subclass = globals()[directive]
            try:
                return subclass.fromfile(fname)
            except NotImplementedError:
                warnings.warn('Class %s has no fromfile method. Returning a generic Target' % subclass)

        return cls(**t_vals)

    @classmethod
    def _read_values(cls, fname):
        """
        Load target definition from the specified file.        
        """

        f = open(fname, 'Ur')
        lines = [fileio._parseline(l) for l in f.readlines()]
        f.close()
    
        values = {}
        values['directive'] = lines[10] 
        values['sh_param'] = tuple(map(int, lines[11].split()))
        n_mat = int(lines[12])
        values['material'] = lines[13: 13+n_mat]

        values['aeff'] = ranges.How_Range.fromstring(lines[29+n_mat])

        return values

class Target_Builtin(Target):
    """
    Base class for target geometries that are built into DDSCAT
    
    :param d: The dipole density. Default is taken from targets.default_d.        
    """
    def __init__(self, directive, d=None, *args, **kwargs):
        Target.__init__(self, directive, *args, **kwargs)
        
        self.d = d if d else default_d

    @classmethod
    def fromfile(cls, fname):
        """
        Dummy method for loading target from file.
        
        Subclasses must implement this themselves. The form is to load the 
        target information from the par file using _read_values(). Then parse
        those values into the form used by the class's initializer..
        """
        raise NotImplementedError('Subclasses should implement this method')
        
    @property
    def N(self):
        """Calculate the numner of dipoles in the target"""
    
        return self._calc_N(self.sh_param)
        
    @property
    def aeff(self):
        """Calculate the effective diameter of the target"""
    
        return self._calc_size(N=self.N, d=self.d)           
        #return (self.N*3/4/np.pi)**(1/3)*self.d

    @staticmethod
    def _calc_size(aeff=None, N=None, d=None):
        """
        Takes two of (aeff, d, N) and returns the third.
        """

        if (aeff and d) and not N:
            return (aeff/d)**3 * (4/3*np.pi)
        elif (aeff and N) and not d:
            return aeff / (N*3/4/np.pi)**(1/3)
        elif (d and N) and not aeff:
            return (N*3/4/np.pi)**(1/3) * d            
        else:
            raise ValueError('Requires two out of three parameters')


class Periodic1D(Target):
    """Base class for 1D periodic targets"""
    pass

class Periodic2D(Target):
    """Base class for 2D periodic targets"""
    pass

class RCTGLPRSM(Target_Builtin):        
    """A rectangular prism target
    
    :param phys_shape: (length, width, height) of the prism in microns
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.
    """

    def __init__(self, phys_shape, d=None, material=None, folder=None):

        Target_Builtin.__init__(self, 'RCTGLPRSM', d=d, material=material, folder=folder)
        self.phys_shape = phys_shape

    @property
    def sh_param(self):
        """Calculate the shape parameters based on the physical shape"""
        
        return np.around(self.phys_shape/self.d).astype(int)

    @staticmethod
    def _calc_N(sh_param):
        """Calculate the number of dipoles
        
        :param sh_param: size in dipoles
        """
        return np.array(sh_param).prod()        

    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified .par file.
        """
        vals = cls._read_values(fname)
        N = cls._calc_N(vals['sh_param'])
        aeff = vals['aeff'].first
        d = cls._calc_size(aeff=aeff, N=N)
        phys_shape = np.array(vals['sh_param']) * d
        return cls(phys_shape, d, vals['material'])

class CYLNDRCAP(Target_Builtin):
    """
    Homogeneous, isotropic finite cylinder with hemispherical endcaps.

    :param length: the length of the cylinder in microns (not including endcaps)
    :param radius: the radius of the cylinder    
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.

    Total height of the structureis length+2*rad
    
    """

    def __init__(self, length, radius, d=None, material=None, folder=None):

        Target_Builtin.__init__(self, 'CYLNDRCAP', d=d, material=material, folder=folder)
        self.length=length
        self.radius=radius

    @property
    def sh_param(self):
        """Calculate the shape parameters"""
        self.sh_param = (int(round(self.length/self.d)),
                         int(round(2 * self.radius/self.d)),
                         0)
                                 
    @staticmethod
    def _calc_N(sh_param):
        """Calculate the number of dipoles
        
        :param sh_param: size in dipoles
        """
        length = sh_param[0]
        diam = sh_param[1]
        Vcyl=length*(np.pi*(diam/2)**2)
        Vsph=4/3*np.pi*(diam/2)**3
        return int(Vcyl+Vsph)        
    
    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified .par file.
        """
        vals = cls._read_values(fname)
        sh_param = vals['sh_param']
        N = cls._calc_N(sh_param)
        aeff = vals['aeff'].first
        d = cls._calc_size(aeff=aeff, N=N)
        phys_shape = np.array(sh_param) * d
        length , radius = phys_shape[0], phys_shape[1]/2
        return cls(length, radius, d, vals['material'])
        

class ELLIPSOID(Target_Builtin):        
    """
    An Ellipsoid target

    :param semiaxes: 3-tuple giving the lengths of the three semiaxes
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.
    """
    
    def __init__(self, semiaxes, d=None, material=None, folder=None):
                
        Target_Builtin.__init__(self, 'ELLIPSOID', d=d, material=material, folder=folder)
        self.semiaxes = np.array(semiaxes)

    @property
    def sh_param(self):
        """Calculate the shape parameters"""

        return 2 * np.around(self.semiaxes/self.d).astype(int)
        
    @staticmethod
    def _calc_N(sh_param):
        """Calculate the number of dipoles
        
        :param sh_param: size in dipoles
        """

        return int(4/3*np.pi*(sh_param.prod()/8))
    
    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified .par file.
        """
        vals = cls._read_values(fname)
        sh_param = vals['sh_param']
        N = cls._calc_N(sh_param)
        aeff = vals['aeff'].first
        d = cls._calc_size(aeff=aeff, N=N)
        semiaxes = np.array(sh_param) * d / 2
        return cls(semiaxes, d, vals['material'])


class Sphere(ELLIPSOID):  
    """
    A Sphere target.

    :param radius: the radius of the sphere
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.
    """      
    def __init__(self, radius, d=None, material=None, folder=None):
        ELLIPSOID.__init__(self, (radius,)*3, d=d, material=material, folder=folder)


class CYLINDER(Target_Builtin):
    """
    Homogeneous, isotropic finite cylinder    

    :param length: the length of the cylinder in microns (not including endcaps)
    :param radius: the radius of the cylinder
    :param orient: the orientation of the cylinder
                SHPAR3 = 1 for cylinder axis aˆ1 ∥ xˆTF: aˆ1 = (1, 0, 0)TF and aˆ2 = (0, 1, 0)TF;
                SHPAR3 = 2 for cylinder axis aˆ1 ∥ yˆTF: aˆ1 = (0, 1, 0)TF and aˆ2 = (0, 0, 1)TF;
                SHPAR3 = 3 for cylinder axis aˆ1 ∥ zˆTF: aˆ1 = (0, 0, 1)TF and aˆ2 = (1, 0, 0)TF in the TF.
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.

    """
    
    def __init__(self, length, radius, orient, d=None, material=None, folder=None):    

        self.length=length
        self.radius=radius
        self.orient=orient
        
        Target_Builtin.__init__(self, 'CYLINDER1', d=d, material=material, folder=folder)
        
    @property
    def sh_param(self):
        """Calculate the shape parameters"""

        return (int(round(self.length/self.d)),
                int(round(2*self.radius/self.d)),
                          self.orient)

    @staticmethod
    def _calc_N(sh_param):
        """Calculate the number of dipoles
        
        :param sh_param: size in dipoles
        """
        return int(sh_param[0]*(np.pi*(sh_param[1]/2)**2))
    
    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified .par file.
        """
        vals = cls._read_values(fname)
        sh_param = vals['sh_param']
        N = cls._calc_N(sh_param)
        aeff = vals['aeff'].first
        d = cls._calc_size(aeff=aeff, N=N)
        phys_shape = np.array(sh_param) * d
        length , radius = phys_shape[0], phys_shape[1]*2
        return cls(length, radius, d, vals['material'])
        
            

### Arbitrarily Shaped Targets

        
class FROM_FILE(Target):
    '''
    Base class for targets of arbitrary geometry.

    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.

    If anisotropic, the “microcrystals” in the target are assumed to be aligned
    with the principal axes of the dielectric tensor parallel to the TF axes.

    FROM_FILE does not inherit from Target_Builtin so for the purposes of inheritence
    it is not considered a builtin target.
    
    '''
    def __init__(self, d=None, material=None, folder=None):    
        Target.__init__(self, 'FROM_FILE', d=d, material=material, folder=folder)

        self.descrip=''
        self.fname='shape.dat'
        self.descriptor='FROM_FILE'
    
        self.grid=np.zeros((0,0,0,3), dtype=int)

        self.a1=np.array([1,0,0])
        self.a2=np.array([0,1,0])
        self.rel_d=np.array([1,1,1])
        self.origin=np.array((0,0,0))           

    @property
    def aeff(self):
        """Calculate the effective diameter of the target"""
                    
        return (self.N*3/4/np.pi)**(1/3)*self.d

    @property
    def d_shape(self):
        """
        The dimensions of the target in dipole units
        """
        return np.array(self.grid.shape[:3])
    
    @property
    def N(self):
        """The number of dipoles"""
        flat = self.grid.sum(3).astype(np.bool)        
        return flat.sum()
        
#    def refresh_N(self):
#        """Update the number of dipoles"""
#        flat = self.grid.sum(3).astype(np.bool)        
#        self.N = flat.sum()
    
    def write(self):
        """Write the shape file."""
        self.refresh_N()
        with open(os.path.join(self.folder, self.fname), 'wb') as f:
            f.write(self.descriptor+'\n') 
            f.write(str(self.N)+'\n')
            f.write(str(self.a1)[1:-1]+'\n')
            f.write(str(self.a2)[1:-1]+'\n')
            f.write(str(self.rel_d)[1:-1]+'\n')
            f.write(str(self.origin)[1:-1]+'\n')
            f.write('J JX JY JZ ICOMPX,ICOMPY,ICOMPZ'+'\n')
            
            if len(self.grid.shape) == 3: # Isotropic case
                grid = np.dstack((self.grid, self.grid, self.grid))
            else:
                grid = self.grid
            
            n=1
            for (ni, i) in enumerate(grid):
                for (nj, j) in enumerate(i):
                    for (nk, k) in enumerate(j):
                        if any(k):
                            f.write('%d %d  %d  %d  %d  %d  %d\n' % ((n, ni, nj, nk)+tuple(k)))
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

class Iso_FROM_FILE(FROM_FILE):
    '''
    Base class for targets of arbitrary geometry with isotropic materials.

    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.
    
    The major difference is that the grid in this case has only one value
    at each dipole.    
    
    '''

    def __init__(self, d=None, material=None, folder=None):    
        FROM_FILE.__init__(self, d=d, material=material, folder=folder)

        self.grid=np.zeros((0,0,0), dtype=int)


def triple(func):
    """
    Decorator to turn a function returning only one value into one that returns
    a 3-tuple of that value.
    """
    
    def triple_func(*args):
        val =func(*args)
        return (val,)*3
        
    return triple_func

def index_space(func, d, offset):
    """
    Take a function that accepts coordinates in physical units and make it
    accept units in dipole space.

    :param func: a function that accepts three arguments (x,y,z) in um.
    :param d: the dipole spacing.
    :param offset: an offset, in dipole units.
    """

    def index_func(i,j,k):
        pt = np.array((i,j,k)) * d + offset
        return func(pt[0], pt[1], pt[2])

    return index_func
    

def Target_fromfunction(func, pt1, pt2, origin=None, d=None, material=None, folder=None):
    """
    Generate a target from a function.
    
    :param func: The generating function. Returns integers corresponding to material
                 index when evaluated at each point in the target volume
    :param pt1: One corner of the target volume, (xmin, ymin, zmin), in um.
    :param pt2: The opposite corner of the target volume, (xmax, ymax, zmax), in um.
    :param origin: The origin of the target (in um)
    The shape, in units of dipoles, of the 
    
    """

    pt1=np.array(pt1)
    pt2=np.array(pt2)
    
    if origin:
        origin=np.array(origin)
    else:
        origin = 0.5 * (pt1+pt2)

    val = func(origin)
    
    target = FROM_FILE(d=d, material=material, folder=folder)        
    
    target.phys_shape = pt2 - pt1

    d_shape = np.ceil((pt2-pt1)/target.d).astype(int)
    d_pt1 = np.around(pt1/target.d).astype(int)

    grid = np.fromfunction(index_func(func, d, ), d_shape)
    
    
class Ellipsoid_FF(Iso_FROM_FILE):
    """
    Build an ellipsoidal target to be loaded from file

    :param semiaxes: is the length of the three semi-major axes in physical units
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.
    """

    def __init__(self, semiaxes, d=None, material=None, folder=None):
        """
        Create a new Ellipsoid Target
        
        """
        Iso_FROM_FILE.__init__(self, d=d, material=material, folder=folder)

        self.phys_shape = 2 *np.array(semiaxes)
        self.descriptor='Ellipsoid_FF (%f, %f, %f, %f)'%(tuple(self.phys_shape)+(self.d,))            
        
        (a,b,c) = tuple(self.d_shape)
        xx,yy,zz=np.mgrid[-a:a, -b:b, -c:c] 
        dist=(xx/a)**2 + (yy/b)**2 + (zz/c)**2

        self.grid = (dist<1).astype(int)

class Helix(Iso_FROM_FILE):
    """
    A helix target.

    dimensions are physical, in um

    :param height: the height of the helix (not counting it's thickness)
    :param pitch: the helix pitch, um/turn
    :param major_r: the radius of the helix sweep
    :param minor_r: the radius of the wire that is swept to form the helix    
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is 'Au_Palik.txt'
    :param folder: The target working directory. The default is the CWD.
    
    """
    def __init__(self, height, pitch, major_r, minor_r, d=None, material=None, folder=None):

        Iso_FROM_FILE.__init__(self, d=d, material=material, folder=folder)

        self.height=height
        self.pitch=pitch
        self.major_r=major_r
        self.minor_r=minor_r

        self._build_helix()

    def _build_helix(self):

        d_shape=(np.asarray([height+2*minor_r,
                            2*(major_r+minor_r),
                            2*(major_r+minor_r)])/self.d).astype(np.int16)
        d_shape+=1
        self.grid=np.zeros(d_shape, dtype=int)

        self.descriptor='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            
        
        print 'Generating Helix...'
        #the helix dimensions in pixels
        p_height=self.height/self.d
        p_pitch=-self.pitch/self.d
        p_major_r=self.major_r/self.d
        p_minor_r=self.minor_r/self.d        

        #the sweep path
        t=np.linspace(0,1, self.height/abs(self.pitch)*360)

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
                return 1
            else:
                return 0
                
        for i in range(d_shape[0]):
            for j in range(d_shape[1]):
                for k in range(d_shape[2]):
                    self.grid[i,j,k]=sphere_check(i,j,k)

        self.refresh_N()

class SpheresHelix(Iso_FROM_FILE):
    """
    A helix target composed of isolated spheres
    
        
    Dimensions are physical, in um
    :param height: the height of the helix (not counting it's thickness)
    :param pitch: the helix pitch, um/turn
    :param major_r: the radius of the helix sweep
    :param minor_r: the radius of the spheres that compose the helix
    :param n_sphere: the number of spheres that compose the helix
    :param d: the dipole dipole spacing, if not specified default value is used
    :param build: if False, delays building the helix until requested. default is True
        
    """
    def __init__(self, height, pitch, major_r, minor_r, n_sphere, d=None, build=True, **kwargs):

        if d is None:
            d=default_d
        d_shape=np.int16(np.asarray([height+2*minor_r, 2*(major_r+minor_r), 2*(major_r+minor_r)])/d)
        d_shape+=1
        Iso_FROM_FILE.__init__(self, d_shape, **kwargs)

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
                return 1
            else:
                return 0
                
        for i in range(self.d_shape[0]):
            for j in range(self.d_shape[1]):
                for k in range(self.d_shape[2]):
                    self.grid[i,j,k]=sphere_check(i,j,k)

        self.refresh_N()


class Conical_Helix(Iso_FROM_FILE):
    """
    A helix target

    dimensions are physical, in um

    :param height: the height of the helix (not counting it's thickness)
    :param pitch1: the starting helix pitch, um/turn (if pitch1 is different from pitch2 will swipe between during the growth)
    :param pitch2: the final helix pitch, um/turn (if pitch1 is different from pitch2 will swipe between during the growth)
    :param major_r: the radius of the helix sweep
    :param minor_r: the radius of the wire that is swept to form the helix
    :param d: the dipole dipole spacing, if not specified default value is used
    :param build: if False, delays building the helix until requested. default is True
    :param **kwargs: are passed to Target
    
    """
    def __init__(self, height, pitch1, pitch2, major_r, minor_r, d=None, build=True, **kwargs):

        if d is None:
            d=default_d
        d_shape=np.int16(np.asarray([height+2*minor_r, 2*(major_r+minor_r), 2*(major_r+minor_r)])/d)
        d_shape+=1
        Iso_FROM_FILE.__init__(self, d_shape, **kwargs)

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
                return 1
            else:
                return 0
                
        for i in range(self.d_shape[0]):
            for j in range(self.d_shape[1]):
                for k in range(self.d_shape[2]):
                    self.grid[i,j,k]=sphere_check(i,j,k)
        
        self.refresh_N()

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
    d_shape=target.grid.shape
    
    if posns is None:
        if num is None:
            raise ValueError('Either posns or num must be specificed')
        
        posns=np.random.rand(num, 3)
        for (p, s) in zip(posns.T, d_shape):
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
        mask_r=np.ravel_multi_index(mask+p, d_shape, mode='clip')
        print mask_r[0]
#        g_len=grid.sum()
        print grid[mask_r[0]]
        grid[mask_r] &= False
        print grid[mask_r[0]]

#        print g_len, '->', grid.sum()

    print grid.sum()
    new_target.N = new_target.grid.sum()
    return new_target
    