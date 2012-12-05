# -*- coding: utf-8 -*-
"""@namespace ScatPy.utils
A few utility functions.

"""

from __future__ import division
import os
import os.path
import glob
import numpy as np
import shutil
import zipfile
from numpy.linalg import norm

from config import exec_settings

def str2pol(s):
    """
    Convert a string into a complex trhee vector indicating the polarization
    
    This uses the convention for helicity described in DDscat user guide
    THIS IS THE OPPOSITE OF THAT USED BY CD SPECTROSCOPISTS
    
    """

    s=s.lower()
    if s=='cl':
        return np.array([0, 0+1j, 1+0j])    
    elif s=='cr':
        return np.array([0, 1+0j, 0+1j])            
    elif s=='lh':
        return np.array([0, 1+0j, 0+0j])    
    elif s=='lv':
        return np.array([0, 0+0j, 1+0j])        
    else:
        raise(ValueError, 'Unknown polarization string %s'%s)        


def normalize(v):
    """
    Normalizes a 3D complex vector
    
    """
    
    return v/np.sqrt(norm(v*v.conj()))

def n_dist(v1, v2):
    """
    The distance between two nromalized complex vectors
    
    """
    v1=normalize(v1)
    v2=normalize(v2)
    
    return np.sqrt(norm(v1-v2))
    

def pol2str(v):
    """
    Convert a vector into a string
    
    This uses the convention for helicity described in DDscat user guide
    THIS IS THE OPPOSITE OF THAT USED BY CD SPECTROSCOPISTS
    
    """

    threshold=1e-6
    
    if n_dist(v, str2pol('cL'))<threshold:
        return 'cL'
    elif n_dist(v, str2pol('cR'))<threshold:
        return 'cR'
    elif n_dist(v, str2pol('lH'))<threshold:
        return 'lH'
    elif n_dist(v, str2pol('lV'))<threshold:
        return 'lV'

    else:
        raise(ValueError, 'Unknown polarization state %s'%str(v))        
    



def str_complex_v(v):
    '''
    Convert a complex three-vector into a string of three tuples
    '''
    return '(%f, %f)  (%f, %f)  (%f, %f)' % (v[0].real, v[0].imag,
                                            v[1].real, v[1].imag,
                                            v[2].real, v[2].imag)

def str2complexV(s):
    """
    Convert a string to a complex vector
    
    """

    v=np.zeros(3, dtype=complex)    
    x=s.replace('(', '').split(')') 
    for i in range(3):
        c=x[i].split(',')
        v[i]=float(c[0])+float(c[1])*1j
    
    return v                


def resolve_mat_file(material):
    '''
    Return an absolute file name pointing to a material file
    
    If the path is already absolute then return that.
    If it's relative return that with ~ expanded
    If it's only a filename, assume that file is found in the materials library
    
    '''
    path=exec_settings['path']

    if path.isabs(material):
        return material
    if path.dirname(material)<>'':
        return path.expanduser(material)
    else:
        return path.normpath(path.join(path.expanduser(exec_settings['mat_library']), material))



def compress_files(folder=None, recurse=False):
    """
    Zips all the support file output by ddscat into one archive
    
    Collects all .avg, .sca, .fml files into their own zip file.
    These can be accessed conveniently with a ZipCollection result object

    The option recurse=True does the same on all subdirectories.    
    
    """
    
    if folder is None:
        folder='.'

    if recurse:
        for (dirpath, dirnames, filenames) in os.walk(folder):
            compress_files(dirpath)
        
        return

    supdir=os.path.join(folder, 'support')

    for ext in ['avg', 'sca', 'fml']:

        if len(glob.glob1(folder, '*.'+ext)):

            zname=os.path.join(folder, 'all_'+ext+'.zip')
#            if os.path.exists(zname):
#                ans=raw_input('%s alread exists. Do you wish to overwrite? [n]\y: '%zname)
#                if ans.lower()<>'y':
#                    return
                    
            z=zipfile.ZipFile(zname, 'w', zipfile.ZIP_DEFLATED)
            
            try:
                os.mkdir(supdir)
            except OSError:
                pass
            
            for f in glob.glob(os.path.join(folder, '*.'+ext)):
                print f
                name=os.path.basename(f)
                z.write(f, name)
                shutil.move(f, supdir)

            z.close()

def MixMaterials(m1, m2, p1):

    """
    Approximate the index for alloys by a weighted average of the two components
    


    """    
    
    n=m1.copy()    
    
    for w in n.data:
        w[1:]*=p1

        w2=m2(w[0])
        w2=(1-p1)*np.array([w2['Rem'], w2['Imm']])        
        w[1:]+=w2
    
    n.fname=''
    hdr=n.header.splitlines()
    hdr[0]='Mixture. %0.1f%% %s : %s' % (p1*100, m1.fname, m2.fname)
    n.header='\n'.join(hdr)
    return n