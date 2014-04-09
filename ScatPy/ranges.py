"""
A collection of ranges found in ddscat.par


"""


from __future__ import division
import numpy as np


class How_Range():
    """
    A general range used for wavelength and scale_range (aeff) definitions.


    :param first: The first value of the range.
    :param last: The last value of the range.
    :param num: The number of steps in the range.
    :param how: An optional string defining the spacing of steps 'LIN', 'INV', 'LOG', 'TAB'.
                Default is 'LIN'.
    :param table: an optional list of table values that specify an arbitrary sequence.
    
    The range can be used as an iterator.        
    """
    def __init__(self, first, last, num, how=None, table=None):

        self.first=first
        self.last=last
        self.num=num
        self.how=how if how is not None else 'LIN'
        
        if self.how=='TAB':
            if self.tab is None:
                raise ValueError('TAB range requires table of values')
        self.build_table()

    def __str__(self):
        """ A string describing the range """
        return '%f  %f  %d  %s'%(self.first, self.last, self.num, self.how)
   
    def build_table(self):
        """ Build the internal representation of the points in the range """
        if self.how=='LIN':
            self.table=np.linspace(self.first, self.last, self.num)
        if self.how=='INV':
            l=self.table=np.linspace(1/self.first, 1/self.last, self.num)
            self.table=1/l
        if self.how=='LOG':
            self.table=np.logspace(np.log10(self.first), np.log10(self.last), self.num)

    def __iter__(self):
        self.build_table()
        self.current=0
        return self
    
    def next(self):
        if self.current==len(self.table):
            raise StopIteration
        
        self.current+=1
        return self.table[self.current-1]

    @classmethod
    def fromstring(cls, string):
        """
        Create a new How_Range based on a string of the form found in DDSCAT.par
        """

        split = string.split()
        return cls(float(split[0]), float(split[1]), int(split[2]), split[3])

class Lin_Range(How_Range):
    """
    A specialized linear range used for specifying target rotations.
    
    The range can be used as an iterator.    

    :param first: The first value of the range.
    :param last: The last value of the range.
    :param num: The number of steps in the range.
    """
    def __init__(self, first, last, num):
        How_Range.__init__(self, first, last, num, 'LIN')

    def __str__(self):
        """ A string describing the range """
        return '%f  %f  %d'%(self.first, self.last, self.num)

    @classmethod
    def fromstring(cls, string):
        """
        Create a new Lin_Range based on a string of the form found in DDSCAT.par
        """

        split = string.split()
        return cls(float(split[0]), float(split[1]), int(split[2]))

     
class Scat_Range():
    '''
    A specialist range used for specifying scattering planes.

    :param phi: The phi scattering angle.
    :param theta_min: The smallest value of theta.
    :param theta_max: The largest value of theta.
    :param d_theta: The theta stepsize.

    Cannot yet be used as an iterator    
    '''
    def __init__(self, phi, theta_min, theta_max, dtheta):

        self.phi=phi
        self.theta_min=theta_min
        self.theta_max=theta_max
        self.dtheta=dtheta

    def __str__(self):
        """ A string describing the range """
        return '%f  %f  %f  %s'%(self.phi, self.theta_min, self.theta_max, self.dtheta)


    @classmethod
    def fromstring(cls, string):
        """
        Create a new Scat_Range based on a string of the form found in DDSCAT.par
        """

        split = string.split()
        return cls(float(split[0]), float(split[1]), float(split[2]), int(split[3]))

class Scat_Range_1dPBC():
    '''
    A specialist range used for specifying scattering cones for 1D periodic targets.

    :param order: The scattering cone.
    :param zeta_min: The smallest value of zeta.
    :param zeta_max: The largest value of zeta.
    :param d_zeta: The zeta stepsize.

    Cannot yet be used as an iterator    
    '''
    def __init__(self, order, zeta_min, zeta_max, d_zeta):

        self.order=order
        self.zeta_min=zeta_min
        self.zeta_max=zeta_max
        self.d_zeta=d_zeta

    def __str__(self):
        """ A string describing the range """
        return '%f  %f  %f  %s'%(self.order, self.zeta_min, self.zeta_max, self.d_zeta)


    @classmethod
    def fromstring(cls, string):
        """
        Create a new Scat_Range_1dPBC based on a string of the form found in DDSCAT.par
        """

        split = string.split()
        return cls(float(split[0]), float(split[1]), float(split[2]), int(split[3]))


class Scat_Range_2dPBC():
    '''
    A specialist range used for specifying scattering diffraction order for 2D periodic targets.

    :param orderM: The scattering cone.
    :param orderN: The scattering cone.

    Cannot yet be used as an iterator    
    '''
    def __init__(self, orderM, orderN):

        self.orderM=orderM
        self.orderN=orderN

    def __str__(self):
        """ A string describing the range """
        return '%d  %d'%(self.orderM, self.orderN)

    @classmethod
    def fromstring(cls, string):
        """
        Create a new Scat_Range_2dPBC based on a string of the form found in DDSCAT.par
        """

        split = string.split()
        return cls(float(split[0]), float(split[1]))
