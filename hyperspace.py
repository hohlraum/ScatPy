# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 09:04:14 2013

@author: andrewmark
"""
import numpy as np
import glob
import re

class ORI_Space(object):
    """
    A class to represent the results for a complete set of sca
    
    
    """
    def __init__(self):
        pass
        
    def __call__(self):
        pass        
        
    def __getitem__(self, val):
    
        print val
        

def ProcessSCASpace():
    W,R,K=set(), set(), set()
    
    for f in glob.glob('*.sca'):
        [_, w, r, k, _] = re.split('w|r|k|\.',f)
        W.add(int(w))
        R.add(int(r))
        K.add(int(k))
    
    return (W,R,K)