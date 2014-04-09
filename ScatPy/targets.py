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
import matplotlib as mpl

import utils
import fileio
import ranges

#: Default spacing between dipoles (in um)
default_d=0.015 

class Target(object):
    """
    The base class for defining DDscat Target geometries.

    :param directive: The type of target (e.g. 'CYLINDER1', 'FROM_FILE')
    :param sh_param: The three shape parameters used by DDSCAT to describe the target
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.

    As it is, this class creates generic targets, which contain only
    a ``directive``, ``sh_param``, ``material`` and ``aeff``: the bare minimum necessary to
    write a ddscat.par file. It is intended to provide an interface as close as
    possible to bare DDSCAT target definitions. It uses dipole units rather than
    subclasses which should use physical units.

    Typically, this class will be subclassed to create more useful and feature-
    rich targets. Derived classes *must* provide the following attributes and
    methods:
    * sh_param: A property that returns the three values of the SHPAR definition used by DDSCAT 
    * _calc_N(): A static method to calculate the number of dipoles based on the shape parameters
    * fromfile(): A class method to generate a working object from a ddscat.par file

    """

    def __init__(self, directive=None, sh_param=None, material=None, aeff=None, folder=None):

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
        
        if isinstance(material, (list, tuple)):
            self.material = list(material)
        elif isinstance(material, str):
            self.material = [material]
        else:
            raise TypeError('Material must be a string or list of strings')
            
        if sh_param:
            self.sh_param = sh_param

        if aeff is not None:
            self.aeff = aeff
                
        if folder is None:
            self._folder='.'
        else:
            self._folder=folder
            
    def save_str(self):
        """Return the multi-line target definition string for the ddscat.par file"""
        out='**** Target Geometry and Composition ****\n'
        out+=self.directive+'\n'
        out+=str(tuple(self.sh_param)).translate(None, '(),')+'\n'
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
        
        sh_param = []
        for s in lines[11].split():
            try:
                sh_param.append(int(s))
            except ValueError:
                sh_param.append(s)
        values['sh_param'] = tuple(sh_param)   

        n_mat = int(lines[12])
        values['material'] = lines[13: 13+n_mat]

        values['aeff'] = ranges.How_Range.fromstring(lines[29+n_mat])

        return values

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

    def make_periodic(self, period):
        """
        Return a periodic version of this target, if one is available.
        
        :period: The unit vector in um. A component = 0: no repetition
        """

        if isinstance(self, Periodic):
            raise TypeError('Target is already periodic')

        new_cls = None
        for old_cls in cls_conversion:
            if isinstance(self, old_cls):
                new_cls = cls_conversion[old_cls]
                break

        if new_cls is None:
            raise TypeError('No corresponding periodic class for %s.' % str(self.__class__))

        new = copy.deepcopy(self)                
        new.__class__ = new_cls
        new.period = period
        new.directive = new_cls.directive
        return new


class Target_Builtin(Target):
    """
    Base class for the standard target geometries that are built into DDSCAT.
    
    :param directive: The type of target (e.g. 'CYLINDER1', 'ELLIPSOID')
    :param d: The dipole density. Default is taken from targets.default_d.        
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.

    Typically, this class will be subclassed to create more useful and feature-
    rich targets. Classes should take their names from the name used by DDSCAT
    (e.g. RCTGLPRSM, ELLIPSOID). 
    
    Target subclasses constructed with :class:``Target_Builtin``
    are intended to work with physical units rather than dipole units used in 
    :class:``Target``. Therefore, derived classes *must* provide the following
    attributes and methods for converting between physical dimensions and the
    internal representation understood by DDSCAT:
    * sh_param: A property that returns the three values of the SHPAR definition used by DDSCAT 
    * _calc_N(): A static method to calculate the number of dipoles based on the shape parameters
    * fromfile(): A class method to generate a working object from a ddscat.par file
    """
    def __init__(self, directive, d=None, material=None, folder=None):

        Target.__init__(self, directive, material=material, folder=folder)
        
        self.d = d if d else default_d
    
    @property
    def sh_param(self):
        """
        Unimplemented method for calculating the shape parameters.

        Subclasses must implement this themselves. The method should take
        the internally stored physical dimensions (e.g. self.radius, self.length)
        and translate them into the shape parameters DDSCAT expects for
        this target type.
        
        This would be a good application for ABCs.
        """
        raise NotImplementedError('Subclasses should implement this method')

    @staticmethod
    def _calc_N(sh_param):
        """
        Unimplemented method for calculating the number of dipoles.

        Subclasses must implement this themselves. The method should accept
        DDSCAT shape parameters and return the number of dipoles for the 
        target they describe. It is a staticmethod because it must be
        accessible by the classmethod fromfile() before the new object is
        instantiated. For convenience it is wrapped by the property obj.N.
        
        This would be a good application for ABCs.
        """
        raise NotImplementedError('Subclasses should implement this method')

    @classmethod
    def fromfile(cls, fname):
        """
        Unimplemented method for loading target from file.
        
        Subclasses must implement this themselves. The form is to load the 
        target information from the par file using _read_values(). Then parse
        those values and pass them to the class's initializer.
        
        This would be a good application for ABCs.
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


class RCTGLPRSM(Target_Builtin):        
    """A rectangular prism target
    
    :param phys_shape: (length, width, height) of the prism in microns
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.
    """

    def __init__(self, phys_shape, d=None, material=None, folder=None):

        Target_Builtin.__init__(self, 'RCTGLPRSM', d=d, material=material, folder=folder)
        self.phys_shape = np.array(phys_shape)

    @property
    def sh_param(self):
        """The shape parameters based on the physical shape"""
        
        return tuple(np.around(self.phys_shape/self.d).astype(int))

    @staticmethod
    def _calc_N(sh_param):
        """Calculate the number of dipoles
        
        :param sh_param: size in dipoles
        """
        return np.array(sh_param[0:3]).prod()        

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


class Cube(RCTGLPRSM):  
    """
    A Cube target.

    :param length: the edge length of the cube
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.
    """      
    def __init__(self, length, d=None, material=None, folder=None):
        RCTGLPRSM.__init__(self, (length,)*3, d=d, material=material, folder=folder)

class CYLNDRCAP(Target_Builtin):
    """
    Homogeneous, isotropic finite cylinder with hemispherical endcaps.

    :param length: the length of the cylinder in microns (not including endcaps)
    :param radius: the radius of the cylinder    
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
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
        return (int(round(self.length/self.d)),
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
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.
    """
    
    def __init__(self, semiaxes, d=None, material=None, folder=None):
                
        Target_Builtin.__init__(self, 'ELLIPSOID', d=d, material=material, folder=folder)
        self.semiaxes = np.array(semiaxes)

    @property
    def sh_param(self):
        """Calculate the shape parameters"""

        return tuple(2 * np.around(self.semiaxes/self.d).astype(int))
        
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
                     file(s) to use for the target. Default is taken from default.par.
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
                SHPAR3 = 1 for cylinder axis a^1 | x^TF: a^1 = (1, 0, 0)TF and a^2 = (0, 1, 0)TF;
                SHPAR3 = 2 for cylinder axis a^1 | y^TF: a^1 = (0, 1, 0)TF and a^2 = (0, 0, 1)TF;
                SHPAR3 = 3 for cylinder axis a^1 | z^TF: a^1 = (0, 0, 1)TF and a^2 = (1, 0, 0)TF in the TF.
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
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
        length, radius, orient = phys_shape[0], phys_shape[1]*2, phys_shape[2]
        return cls(length, radius, orient, vals['material'])
        
            

### Arbitrarily Shaped Targets

        
class FROM_FILE(Target):
    '''
    Base class for targets of arbitrary geometry.

    :param grid: The dipole grid
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.

    If anisotropic, the “microcrystals” in the target are assumed to be aligned
    with the principal axes of the dielectric tensor parallel to the TF axes.

    FROM_FILE does not inherit from Target_Builtin so for the purposes of inheritence
    it is not considered a builtin target.
    
    '''
    def __init__(self, grid=None, d=None, material=None, folder=None):    
        Target.__init__(self, 'FROM_FILE', material=material, folder=folder)

        self.description=''
        self.fname='shape.dat'

        self.d = d if d else default_d    
    
        if grid is not None:
            self.grid=grid
        else:
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
    def sh_param(self):
        """
        The dimensions of the target in dipole units
        """
        return self.grid.shape[:3]
    
    @property
    def N(self):
        """The number of dipoles"""
        if len(self.grid.shape)==4:
            flat = self.grid.sum(3).astype(np.bool)        
        else:
            flat = self.grid.astype(np.bool)
            
        return flat.sum()
            
    def write(self):
        """Write the shape file."""
        with open(os.path.join(self.folder, self.fname), 'wb') as f:
            f.write(self.description+'\n') 
            f.write(str(self.N)+'\n')
            f.write(str(self.a1)[1:-1]+'\n')
            f.write(str(self.a2)[1:-1]+'\n')
            f.write(str(self.rel_d)[1:-1]+'\n')
            f.write(str(self.origin/self.d)[1:-1]+'\n')
            f.write('J JX JY JZ ICOMPX,ICOMPY,ICOMPZ'+'\n')
            
            table = self._grid2table()
            
            for (n, val) in enumerate(table):
                f.write('%d %d  %d  %d  %d  %d  %d\n' % ((n+1, )+tuple(val)))


    def _grid2table(self):
        """
        Convert grid into a table similar to that found in shape.dat
        """

        if len(self.grid.shape) == 4: # Anisotropic case
            grid = self.grid.sum(3)
        else:
            grid = self.grid

        table = np.zeros((self.N, 6))

        entries = np.transpose(grid.nonzero())
        
        for (n, pt) in enumerate(entries):

            val = self.grid[tuple(pt)]
            try:
                len(val)
            except TypeError:
                table[n] = np.hstack((pt,) + (val,)*3)
            else: # Anisotropic case
                table[n] = np.hstack((pt, val))
                
        return table

    @staticmethod
    def _table2grid(table):
        """
        Convert a table representation of the shape into a grid representation.

        Always returns a grid with shape (x,y,z,3) suitable for aniosotropic
        targets.
        Tables are assumed to be 1-based (i.e. the lowest index is 1)
        """
        
        x,y,z = table[:,0]-1, table[:,1]-1, table[:,2]-1
        nx, ny, nz = table[:,3], table[:,4], table[:,5]
#        ncomp = len(set(nx), set(ny), set(nz))

        shape = (max(x)+1, max(y)+1, max(z)+1, 3)
        grid = np.zeros(shape, dtype=int)
        for (comp, n) in enumerate((nx, ny, nz)):
            comp = np.array((comp,)*len(x))
            ind = np.ravel_multi_index((x,y,z, comp), shape)
            grid.put(ind, n)

        return grid

    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified .par file.
        
        Assumes that the accompanying shape.dat file is in the same folder.

        This function currently assumes that the target basis vectors are
        orthonormal.        
        """
        vals = cls._read_values(fname)

        aeff = vals['aeff'].first

        h,_ = os.path.split(fname)
        shape_file = os.path.join(h, 'shape.dat') 
        shape = results.ShapeTable(fname=shape_file)

        grid = cls._table2grid(shape.data[:,1:])        

        d = cls._calc_size(aeff=aeff, N=len(shape.data))

        ncomp = len(vals['material'])

        if ncomp == 1:
            return Iso_FROM_FILE(grid=grid, d=d, material=vals['material'])
        else:
            return cls(grid=grid, d=d, material=vals['material'])

    @classmethod
    def fromfunction(cls, func, pt1, pt2, origin=None, d=None, material=None, folder=None, **kwargs):
        """
        Generate a target from a function.
        
        :param func: The generating function. Returns integers corresponding to material
                     index when evaluated at each point in the target volume
        :param pt1: One corner of the target volume, (xmin, ymin, zmin), in um.
        :param pt2: The opposite corner of the target volume, (xmax, ymax, zmax), in um.
        :param origin: The origin of the target (in um)
        :param kwargs: Further kwargs are passed to func
        
        See numpy.fronfunction for details on the function func.        
        
        """
    
        pt1=np.array(pt1)
        pt2=np.array(pt2)
        
        if origin:
            origin=np.array(origin)
        else:
            origin = 0.5 * (pt1+pt2)
    
        val = func(*origin)
        
        try:
            len(val)
        except TypeError:
            target = Iso_FROM_FILE(d=d, material=material, folder=folder)        
        else:      
            target = FROM_FILE(d=d, material=material, folder=folder)        
        
        target.phys_shape = pt2 - pt1
    
        d_shape = np.ceil((pt2-pt1)/target.d).astype(int)
        d_pt1 = np.around(pt1/target.d).astype(int)

        def index_space(func, d, offset, **kwargs):
            """
            Take a function that accepts coordinates in physical units and make it
            accept units in dipole space.
        
            :param func: a function that accepts three arguments (x,y,z) in um.
            :param d: the dipole spacing.
            :param offset: an offset, in dipole units.
            """

            def index_func(i,j,k, **kwargs):
                x,y,z = (i + offset[0])*d , (j + offset[1])*d, (k + offset[2])*d
                return func(x, y, z, **kwargs)
        
            return index_func
    
        i_func = index_space(func, target.d, d_pt1, **kwargs)
        target.grid = np.fromfunction(i_func, d_shape)
        
        return target

    
    def VTRconvert(self, outfile=None):
        """Execute VTRConvert to generate a model file viewable in Paraview"""
        Target.VTRconvert(self, outfile)


    def show(self, *args, **kwargs):
        """
        Display the dipoles using Mayavi. 
        
        Currently assumes that all dipoles are isotropic.
        """
        from mayavi import mlab

        max_points=20000
        #This mask_point business shouldn't be necessary, but the builtin VTK
        #implementation causes my computer to segfault
        if 'mask_points' in kwargs:
            mask_points=kwargs.pop('mask_points')
        elif self.N>max_points:
            print 'Warning! Large number of datapoints in target.'
            print 'Plotting only a subset. Specify mask_points=None or an integer to force skipping value'
            mask_points=int(self.N/max_points)
        else:
            mask_points=1
        
        table = self._grid2table()        
        
        X=table[:, 0][::mask_points]
        Y=table[:, 1][::mask_points]
        Z=table[:, 2][::mask_points]
        c=table[:, 3][::mask_points]
            
        mlab.points3d(X, Y, Z, c, scale_mode='none', *args, **kwargs)
        mlab.show()
    
    

class Iso_FROM_FILE(FROM_FILE):
    '''
    Base class for targets of arbitrary geometry with isotropic materials.

    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.
    
    The major difference is that the grid in this case has only one value
    at each dipole.    
    
    '''

    def __init__(self, grid=None, d=None, material=None, folder=None):    

        if grid is None:
            grid=np.zeros((0,0,0), dtype=int)

        FROM_FILE.__init__(self, grid, d=d, material=material, folder=folder)

    
class Ellipsoid_FF(Iso_FROM_FILE):
    """
    Build an ellipsoidal target to be loaded from file

    :param semiaxes: is the length of the three semi-major axes in physical units
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.
    """

    def __init__(self, semiaxes, d=None, material=None, folder=None):
        """
        Create a new Ellipsoid Target
        
        """
        Iso_FROM_FILE.__init__(self, d=d, material=material, folder=folder)

        self.phys_shape = 2 *np.array(semiaxes)
        self.description='Ellipsoid_FF (%f, %f, %f, %f)'%(tuple(self.phys_shape)+(self.d,))            
        
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
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.
    
    """
    def __init__(self, height, pitch, major_r, minor_r, d=None, material=None, folder=None):

        Iso_FROM_FILE.__init__(self, d=d, material=material, folder=folder)

        self.height=height
        self.pitch=pitch
        self.major_r=major_r
        self.minor_r=minor_r
        self.origin = np.array((height + minor_r, minor_r+major_r, minor_r+major_r))

        self._build_helix()

    def _build_helix(self):

        d_shape=(np.asarray([self.height+2*self.minor_r,
                            2*(self.major_r+self.minor_r),
                            2*(self.major_r+self.minor_r)])/self.d).astype(np.int16)
        d_shape+=1
        self.grid=np.zeros(d_shape, dtype=int)

        self.description='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            
        
        print 'Generating Helix...'
        #the helix dimensions in pixels
        p_height=self.height/self.d
        p_pitch=-self.pitch/self.d
        p_major_r=self.major_r/self.d
        p_minor_r=self.minor_r/self.d        
        p_origin = self.origin/self.d

        #the sweep path
        t=np.linspace(0,1, self.height/abs(self.pitch)*360)

        x=p_height*t + p_minor_r
        y=p_major_r * np.cos(2*np.pi* x/p_pitch) + p_origin[1]
        z=p_major_r * np.sin(2*np.pi* x/p_pitch) + p_origin[2]
        
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


class SpheresHelix(Iso_FROM_FILE):
    """
    A helix target composed of isolated spheres
            
    Dimensions are physical, in um.
    
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
            self.description='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            

    def build_helix(self):

        self.description='FROM_FILE_Helix (%f, %f, %f, %f, %f)'%(self.height, self.pitch, self.major_r, self.minor_r, self.d)            
        
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
            self.description='FROM_FILE_Helix (%f, %f, %f, %f, %f, %f)'%(self.height, self.pitch1, self.pitch2, self.major_r, self.minor_r, self.d)            

    def build_helix(self):

        self.description='FROM_FILE_Helix (%f, %f,%f, %f, %f, %f)'%(self.height, self.pitch1, self.pitch2, self.major_r, self.minor_r, self.d)            
        
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

class Polygon(Iso_FROM_FILE):
    """
    A planar target in the shape defined by the polygon

    dimensions are physical, in um

    :param polygon: a list of vertices for the polygon
    :param thickness: the thickness of the polygon
    :param bbox: the bounding box. By default uses the smallest box that contains
                 the entire polygon. If bbox is a 2-tuple then the values specify
                 the width and height of the box and the polygon is centred within it.
                 If the bbox is a 4-tuple, then the values specify the bottom
                 left x,y coordinates and the width and height.
    :param d: the dipole dipole spacing, if not specified default value is used
    :param material: A string specifying the material
                     file to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.        
    
    Dipoles within the polygon correspond to material, outside to ambient.    
    
    """
    def __init__(self, polygon, thickness, bbox = None, d=None, material=None, folder=None):

        Iso_FROM_FILE.__init__(self, d=d, material=material, folder=folder)
        
        polygon= np.array(polygon)
        p_min = polygon.min(0)
        p_max = polygon.max(0)
        poly_centre = (p_min+p_max)/2
        poly_size = p_max-p_min
        
        if bbox is None:
            bbox_corner = p_min
            bbox_size = poly_size
        
        elif len(bbox) == 2:
            bbox_corner = poly_centre - np.array(bbox)/2
            bbox_size = np.array(bbox)

        elif len(bbox) == 4:
            bbox_corner = np.array(bbox[:2])
            bbox_size = np.array(bbox[-2:])
        else:
            raise ValueError('Argument bbox must have length 2 or 4')
    
        px_size = np.ceil((bbox_size / self.d)).astype(int)
        scale = min(px_size / bbox_size)
    
        # map the vertices into dipole space
        self.verts = (polygon - bbox_corner)*scale
        path = mpl.path.Path(self.verts)

        x, y = np.meshgrid(np.arange(px_size[0]),np.arange(px_size[1]))
        x, y = x.flatten(), y.flatten()        
        xy = np.vstack((x,y)).T
        
        grid = path.contains_points(xy)
        grid = grid.reshape(px_size).astype(int)
        
        self.grid = np.tile(grid, (round(thickness/self.d), 1, 1))
        
        self.origin = np.array(tuple(poly_centre)+(thickness/2.,))
        

### Periodic Targets

class Periodic(object):
    """Base class for periodic targets

    :param dim: The number of periodic dimensions (1 or 2)    
    :param period: A 2-tuple which defines the period of the array
                        in the x,y TF directions. (in um)

    """
    @property
    def dimension(self):
        """
        Returns the number of periodic dimensions
        """
        return sum(i != 0 for i in self.period)

class FRMFILPBC(FROM_FILE, Periodic):
    """
    Base class for periodic targets of arbitrary geometry.
    
    :param grid: The dipole grid
    :param period: A 2-tuple which defines the period of the array
                        in the x,y TF directions. (in um)
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.

    If anisotropic, the “microcrystals” in the target are assumed to be aligned
    with the principal axes of the dielectric tensor parallel to the TF axes.

    FROMFILPBC does not inherit from Target_Builtin so for the purposes of inheritence
    it is not considered a builtin target.
        
    """
    directive = 'FRMFILPBC'
    def __init__(self, grid=None, period=(0,0), d=None, material=None, folder=None, fname=None):    

        FROM_FILE.__init__(self, grid=grid, d=d, material=material, folder=folder)

#        self.directive = 'FRMFILPBC'
        self.fname = 'shape.dat' if fname is None else fname
        self.period = period

    @property
    def sh_param(self):
        """Return the shape parameter"""
        return (self.period[0]/self.d, self.period[1]/self.d, self.fname)

    @classmethod
    def fromfile(cls, fname):
        """
        Load target definition from the specified .par file.
        
        Assumes that the accompanying shape.dat file is in the same folder.

        This function currently assumes that the target basis vectors are
        orthonormal.        
        """
        vals = cls._read_values(fname)

        aeff = vals['aeff'].first

        h,_ = os.path.split(fname)
        shape_file = vals['sh_param'][2]
        shape_file = os.path.join(h, vals['sh_param'][2]) 
        shape = results.ShapeTable(fname=shape_file)

        grid = cls._table2grid(shape.data[:,1:])        

        d = cls._calc_size(aeff=aeff, N=len(shape.data))
        period = np.array(vals['sh_param'][:2] * d)

        return cls(grid=grid, period=period, d=d, material=vals['material'])


class RCTGL_PBC(RCTGLPRSM, Periodic):
    """
    A target of periodic rectangular prisms.
    
    :param phys_shape: (length, width, height) of the prism in microns
    :param period: A 2-tuple which defines the period of the array
                        in the x,y TF directions. (in um)
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.

    """
    directive = 'RCTGL_PBC'
    
    def __init__(self, phys_shape, period, d=None, material=None, folder=None, fname=None):    

        RCTGLPRSM.__init__(self, phys_shape, d=d, material=material, folder=folder)

        self.directive = 'RCTGL_PBC'
        self.period = period

    @property
    def sh_param(self):
        """Return the shape parameter"""

        return (RCTGLPRSM.sh_param.fget(self)) + (self.period[0]/self.d, self.period[1]/self.d)

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
        period = phys_shape[3:4]
        return cls(phys_shape[:3], period, d, vals['material'])


class CYLNDRPBC(CYLINDER, Periodic):
    """
    A target of periodic cylinders.
    
    :param length: the length of the cylinder in microns (not including endcaps)
    :param radius: the radius of the cylinder
    :param orient: the orientation of the cylinder
                SHPAR3 = 1 for cylinder axis aˆ1 ∥ xˆTF: aˆ1 = (1, 0, 0)TF and aˆ2 = (0, 1, 0)TF;
                SHPAR3 = 2 for cylinder axis aˆ1 ∥ yˆTF: aˆ1 = (0, 1, 0)TF and aˆ2 = (0, 0, 1)TF;
                SHPAR3 = 3 for cylinder axis aˆ1 ∥ zˆTF: aˆ1 = (0, 0, 1)TF and aˆ2 = (1, 0, 0)TF in the TF.
    :param period: The periodicity 
    :param period: A 2-tuple which defines the period of the array
                        in the x,y TF directions. (in um)
    :param d: The dipole density. Default is taken from targets.default_d.
    :param material: A string, or list of strings specifying the material
                     file(s) to use for the target. Default is taken from default.par.
    :param folder: The target working directory. The default is the CWD.

    """
    
    def __init__(self, length, radius, orient, period, d=None, material=None, folder=None):    


        CYLINDER.__init__(self, length, radius, orient, d=d, material=material, folder=folder)

        self.directive = 'CYLNDRPBC'
        self.period = period

    @property
    def sh_param(self):
        """Return the shape parameter"""

        return (CYLINDER.sh_param.fget(self)) + (self.period[0]/self.d, self.period[1]/self.d)

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
        length, radius, orient = phys_shape[0], phys_shape[1]*2, phys_shape[2]
        period = phys_shape[3:4]
        return cls(length, radius, orient, period, d, vals['material'])

#: A dict which translates between a finite isolated target class and its 
#: corresponding (semi)infite periodic partner. The keys are the classes of
#: the isolated target, and the values are the classes of the periodic ones.
cls_conversion={RCTGLPRSM:RCTGL_PBC,
                FROM_FILE:FRMFILPBC,
                CYLINDER :CYLNDRPBC}

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
    