# -*- coding: utf-8 -*-
"""
Utility functions.

"""

from __future__ import division
import os
import os.path
import glob
import numpy as np
import shutil
import zipfile
from pkg_resources import resource_exists, resource_filename, resource_listdir, cleanup_resources

from numpy.linalg import norm

import core

def gauss(x, sigma):
    
    return 1/(np.sqrt(2*np.pi) *sigma) * np.exp(-x**2/2/sigma**2)  


def r_parser(k):
    return float(k.split('r')[1].split('/')[0])

def L_parser(k):
    return float(k.split('L')[1].split('/')[0])


def weighted_gauss(c, parser, sigma):
    """
    Make a Normalized-Gaussian weighted sum of a range of spectras obtained
    with a series of simulations sweeping on a folder-collection names obtained parameters (e.g. r or L)
    from the parser functions
    
    It is written on the data structure obtained from the results.FolderCollection
    
    """

    ans=c[c.keys()[0]].copy()
    
    #Initialize new
    x_field=ans.x_field
    for k in ans.keys():
        if k != x_field:
            ans[k] *= 0
    
    gsum=0
    for spec_k in c.keys():
        x=parser(spec_k)
        for field_k in ans.keys():
            if field_k != x_field:
                g=gauss(x, sigma)
                ans[field_k] += g * c[spec_k][field_k]
                gsum+=g

    #Normalize
    for field_k in ans.keys():
        if field_k != x_field:
            ans[field_k] /= gsum
        
    return ans



def str2pol(s):
    """
    Convert a string into a complex three vector indicating the polarization
    
    Attempts to identify polarizations by name: 'cl', 'cr', 'lh', 'lv'
    
    The helicity convention used for named polarizations is determined by the
    values of the polarization constants in core.
    
    """

    s=s.lower()
    if s=='cl':
        return core.pol_cL
    elif s=='cr':
        return core.pol_cR
    elif s=='lh':
        return core.pol_lH    
    elif s=='lv':
        return core.pol_lV
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
    
    The helicity convention used for named polarizations is determined by the
    values of the polarization constants in core.
    """

    threshold=1e-6
    
    if n_dist(v, core.pol_cL)<threshold or n_dist(-v, core.pol_cL)<threshold:
        return 'cL'
    elif n_dist(v, core.pol_cR)<threshold or n_dist(-v, core.pol_cR)<threshold:
        return 'cR'
    elif n_dist(v, core.pol_lH)<threshold or n_dist(-v, core.pol_lH)<threshold:
        return 'lH'
    elif n_dist(v, core.pol_lV)<threshold or n_dist(-v, core.pol_lV)<threshold:
        return 'lV'

    else:
        raise(ValueError, 'Unknown polarization state %s'%str(v))        
    



def complexV2str(v):
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
    from core import config

    path=config['path_style']

    if path.isabs(material):
        return material
    if path.dirname(material)<>'':
        return path.expanduser(material)
    else:
        return path.normpath(path.join(path.expanduser(config['mat_library']), material))

def resolve_profile(fname):
    """
    Resolve the absolute file name of the requested profile file.

    :param fname: The name of the file to find
    :returns: Absolute path to the file. None if not found.    
    
    The search scheme is to first look for the file in the CWD, followed by
    the folder ~/.ScatPy/ and finally the subdiretory profiles/ relative to
    where the utils.py module resides.

    """

    pkg_path = os.path.dirname(__file__)

    for path in ['./', os.path.expanduser('~/.ScatPy/'), os.path.join(pkg_path, 'profiles')]:
        full_name=os.path.join(path, fname)
        if os.path.exists(full_name):            
            break
        else:
            full_name = None

    if full_name is None and resource_exists('ScatPy', os.path.join('profiles', fname)):
        return resource_filename('ScatPy', os.path.join('profiles', fname))

    return full_name

def make_profile():
    """
    Create a local profile in user's home folder.    
    """

    home_path = os.path.expanduser('~/.ScatPy/')
    try:
        os.makedirs(home_path)
    except OSError:
        raise(IOError('A folder ~./ScatPy already exists. Delete it manually before proceeding.'))

    files = resource_listdir('ScatPy', 'profiles')

    for f in files:
        src = resource_filename('ScatPy', os.path.join('profiles', f))
        shutil.copy(src, home_path)

    cleanup_resources()
    
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

    for ext in ['avg', 'sca', 'fml', 'pol*', 'E*']:

        if len(glob.glob1(folder, '*.'+ext)):

            zname=os.path.join(folder, 'all_'+ext.replace('*', 'n')+'.zip')
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