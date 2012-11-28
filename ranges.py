# -*- coding: utf-8 -*-
"""
A collection of ranges found in ddscat.par


"""


from __future__ import division
import numpy as np


class How_Range():
    """
    A general range used for wavelength and aeff definitions.
    
    The range can be used as an iterator.    
    
    """
    def __init__(self, first, last, num, how=None, table=None):
        """
        Initialize the range
        
        Arguments:
            first: the first value of the range
            last: the last value of the range
            num: the number of steps in the range
            how: an optional string defining the spacing of steps 'LIN', 'INV', 'LOG', 'TAB'
                 default is 'LIN'
            table: an optional list of table values that specify an arbitrary sequence
        
        """
        self.first=first
        self.last=last
        self.num=num
        self.how=how if how is not None else 'LIN'
        
        if self.how=='TAB':
            if self.tab is None:
                raise ValueError('TAB range requires table of values')
        self.build_table()

    def __str__(self):
        return '%f  %f  %d  %s'%(self.first, self.last, self.num, self.how)
   
    def build_table(self):
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

    

class Lin_Range(How_Range):
    """
    A specialized linear range used for specifying target rotations.
    
    The range can be used as an iterator.    
    
    """
    def __init__(self, first, last, num):
        """
        Initialize the range
        
        Arguments:
            first: the first value of the range
            last: the last value of the range
            num: the number of steps in the range
        """
        How_Range.__init__(self, first, last, num, 'LIN')

    def __str__(self):
        return '%f  %f  %d'%(self.first, self.last, self.num)

     
class Scat_Range():
    '''
    A specialist range used for specifying scattering planes.

    Does not yet include an iterator    
    '''
    def __init__(self, phi, theta_min, theta_max, dtheta):
        """
        Initialize the range
        
        Arguments:
            phi: the phi scattering angle
            theta_min: the smallest value of theta
            theta_max: the largest value of theta
            d_theta: the theta stepsize
        """

        self.phi=phi
        self.theta_min=theta_min
        self.theta_max=theta_max
        self.dtheta=dtheta

    def __str__(self):
        return '%f  %f  %f  %s'%(self.phi, self.theta_min, self.theta_max, self.dtheta)


