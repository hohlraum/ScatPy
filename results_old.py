# -*- coding: utf-8 -*-
"""
A set of classes for managing the output files from ddscat

"""

from __future__ import division
import numpy as np
import os
import os.path
import glob
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from collections import OrderedDict

import utils


#class Results():
#    def __init__(self, **kwargs):
##        self.log=Log()
#        self.refresh()
#        self.mtable=MTable(**kwargs)
#        self.qtable=QTable(**kwargs)
#        self.qtable2=QTable2(**kwargs)
#
#    def refresh(self):        
##        self.log.refresh()
#        self.mtable.refresh()
#        self.qtable.refresh()
#        self.qtable2.refresh()
#
#    def set_folder(self, newfolder):
#        self.mtable.set_folder(newfolder)
#        self.qtable.set_folder(newfolder)        
#        self.qtable2.set_folder(newfolder)
#    

class ResultTable():
    """
    Base class for results tables.        
    """
    def __init__(self, fname, hdr_len, c_width, folder=None):
        """
        Initialize the result table.
        
        Arguments:
            fname: the file name, resolved relative to the working folder
            hdr_len: the number of header lines to ignore
            c_width: a list of column widths
            
        Optional arguments:
            folder: the working folder from which to read the table. Default is
                    the cwd.


        
        
        Most DDscat output files have fixed column widths and cannot be relied upon 
        to have whitespace between adjacent columns. Therefore the columns' width must
        be provided in the c_width list.   

        Field names are derived from the column headings in the result table. The
        fields can be accessed like a Python dict, e.g. T['wave']                            
        """
        if folder is None:
            self.folder=''
        else:
            self.folder=folder
            
        self.fname=fname
        self.hdr_len=hdr_len
        self.c_width=c_width
        self.header=None
        self.col_lbl=None
        try:
            self.refresh()
        except (IOError):
            print 'Result Table %s not found' % utils.resolve_mat_file(os.path.join(self.folder, self.fname))
            self.data=None        
        except (ValueError):
            print 'Result Table %s format not recognized (wrong colum format?)' % utils.resolve_mat_file(os.path.join(self.folder, self.fname))

            self.data=None

    def refresh(self):
        """
        Refresh the data from the file
                
        """
        
        with open(os.path.join(self.folder, self.fname), 'Ur') as f:
            hdr=''
            for i in range(self.hdr_len):
                hdr+=f.readline()
            self.header=hdr
    
            l=f.readline()
            labels=split_string(l, self.c_width)
            for i,l in enumerate(labels):
                labels[i]=l.strip()
            self.col_lbl=labels
    
            l=f.readline()
            dat=np.asarray(map(float, split_string(l, self.c_width)))
            l=f.readline()
            while  l <> '':
                dat=np.vstack([dat, np.asarray(map(float, split_string(l, self.c_width)))])
                l=f.readline()                
            
            self.data=dat
            for (l,d) in zip(self.col_lbl, self.data.transpose()):
                setattr(self, clean_string(l), d)

    def set_folder(self, new_folder):
        """
        Change the working folder.        
        """
        self.folder=new_folder

    def plot(self, fields, **kwargs):

        ylbl=''
        for f in fields[1:]:
            ylbl+=f+', '
        ylbl=ylbl[:-2]
        
        plt.ylabel(ylbl)
        for (n,v) in dat.iteritems():
            for (i, y) in enumerate(v[1:]):
                label=n #+'_'+fields[i+1]
                if normalize==True:
                    maxy=y.max()
                else:
                    maxy=1
                plt.plot(v[0], y/maxy, label=label, lw=lw, **kwargs)
    
class AVGTable(ResultTable):
    """
    A class for reading .avg files output by DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='w000r000.avg'
        ResultTable.__init__(self, fname, 37, [6, 7,9,12,12,11,11,11,11,11,11,11], **kwargs)
        self.summary=AVGSummary(self.header)


class AVGSummary():
    """
    A class for reading the summary section of AVGTables.

    """      
    def __init__(self, hdr):
                
        hdr=hdr.splitlines()
        
        l=hdr[9].split()
        self.wave=float(l[1])
        
        l=hdr[28].split()
        self.J1_Qext=float(l[1])
        self.J1_Qabs=float(l[2])
        self.J1_Qsca=float(l[3])
        
        l=hdr[29].split()
        self.J2_Qext=float(l[1])
        self.J2_Qabs=float(l[2])
        self.J2_Qsca=float(l[3])
        



class QTable(ResultTable):
    """
    A class for reading qtable output files from DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='qtable'
        ResultTable.__init__(self, fname, 13, [10,11,11,11,11,12,11,11,6], **kwargs)

class QTable2(ResultTable):
    """
    A class for reading qtable2 output files from DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='qtable2'
        ResultTable.__init__(self, fname, 12, [10,11,12,12,12], **kwargs)


class MTable(ResultTable):
    """
    A class for reading mtables generated by DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='mtable'
        ResultTable.__init__(self, fname, 7, [10,11,10,10,10,10], **kwargs)

        
class MInTable(ResultTable):
    """
    A class for reading material input files used by DDscat.
    
    Simple file names are resolved relative to the default material_library path    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='mintable'
        ResultTable.__init__(self, fname, 2, None, **kwargs)
    
    def refresh(self):
        fname=utils.resolve_mat_file(os.path.join(self.folder, self.fname))

        with open(fname, 'Ur') as f:
            hdr=''
            for i in range(self.hdr_len):
                hdr+=f.readline()
            self.header=hdr
    
            l=f.readline()
            labels=l.split()
            for i,l in enumerate(labels):
                labels[i]=l.strip()
            self.col_lbl=labels
    
            l=f.readline()
            dat=np.asarray(map(float, l.split()))
            l=f.readline()
            while  l <> '':
                dat=np.vstack([dat, np.asarray(map(float, l.split()))])
                l=f.readline()                
            
            self.data=dat
            for (l,d) in zip(self.col_lbl, self.data.transpose()):
                setattr(self, clean_string(l), d)

            self.interps={}
            for l in self.col_lbl[1:]:
                l=clean_string(l)
                self.interps[l]=interp1d(self.wave[::-1], getattr(self, l)[::-1])                

    def __call__(self, wave):
        out={}
        for l in self.col_lbl[1:]:
            l=clean_string(l)
            out[l]=float(self.interps[l](wave))

        return out            

def collect_results_from_folders(fields, fname=None, path=None):
    '''
    plots the contents of all subfolders

    '''
    if fname is None:
        fname='qtable'

    if fname.lower()=='mtable':
        rtable=MTable
    elif fname.lower()=='qtable':
        rtable=QTable
    elif fname.lower()=='qtable2':
        rtable=QTable2
        
    if path is None:
        path='.'
        
    folders=[i for i in os.listdir(path) if os.path.isdir(i)]

    dat=OrderedDict()
    for f in folders:
        d=rtable(folder=f)
        try:
            s=np.array([getattr(d,a) for a in fields])
            print 'Reading from %s' % f
            dat[f]=s    
        except AttributeError:
            pass
        

    return dat


def dichroism_from_avg(fields, folder=None):
    
    table=[]   
    
    for files in glob.glob("w???r???.avg"):

        dat=AVGTable(files)
        
        s=dat.summary        
        
        t=[s.wave, s.J2_Qext-s.J1_Qext, s.J2_Qabs-s.J1_Qabs, s.J2_Qsca-s.J1_Qsca]
        table.append(t)
    
    
  #  return np.asarray(table)
  #  a=np.zeros(len(table), dtype=[('wave',np.float32),('CDext',np.float32),('CDabs',np.float32), ('CDsca', np.float32)])
#    return np.asarray(table, dtype='f4 f4 f4 f4')
    dat=OrderedDict()
    if 'wave' in fields:
        dat['wave']=np.asarray([t[0] for t in table])
    if 'CDext' in fields:
        dat['CDext']=np.asarray([t[1] for t in table])
    if 'CDabs' in fields:    
        dat['CDabs']=np.asarray([t[2] for t in table])
    if 'CDsca' in fields:
        dat['CDsca']=np.asarray([t[3] for t in table])

    return dat

def plot_dichroism_from_avg(fields=None, lw=None, normalize=False, **kwargs):
    

    if fields==None:
        fields=['wave', 'CDext']

    if lw==None:
        lw=2

    dat=dichroism_from_avg(fields)
    
    plt.xlabel(fields[0])

    ylbl=''
    for f in fields[1:]:
        ylbl+=f+', '
    ylbl=ylbl[:-2]
    
    plt.ylabel(ylbl)

    for f in fields[1:]:
        print f
        y=dat[f]
        label=f
        if normalize:
            maxy=y.max()
        else:
            maxy=1.
            
        plt.plot(dat[fields[0]], y/maxy, label=label, lw=lw, **kwargs)
        

def plot_folders(fields=None, fname=None, path=None, normalize=False, lw=None, **kwargs):
    
    if fields==None:
        fields=['wave', 'Q_ext']
        
    if lw==None:
        lw=2
        
    dat=collect_results_from_folders(fields, fname, path)    
    
    plt.xlabel(fields[0])
    
    ylbl=''
    for f in fields[1:]:
        ylbl+=f+', '
    ylbl=ylbl[:-2]
    
    plt.ylabel(ylbl)
    for (n,v) in dat.iteritems():
        for (i, y) in enumerate(v[1:]):
            label=n #+'_'+fields[i+1]
            if normalize==True:
                maxy=y.max()
            else:
                maxy=1
            plt.plot(v[0], y/maxy, label=label, lw=lw, **kwargs)


def plot_dichroism(fields=None, fname=None, path=None, normalize=False, lw=None, **kwargs):
    
    if fields==None:
        fields=['wave', 'Q_ext']
        
    if lw==None:
        lw=2
        
    dat=collect_results_from_folders(fields, fname, path)    
    
    plt.xlabel(fields[0])
    
    ylbl=''
    for f in fields[1:]:
        ylbl+=f+', '
    ylbl=ylbl[:-2]
    
    plt.ylabel(ylbl)

    CD=OrderedDict()
    for l in dat.keys():
        if l[-3:]=='LCP' or l[-3:]=='_cL':
            print 'Dichroism from: '+l
            if l[-3:]=='LCP':
                match='RCP'
            else:
                match='_cR'
            L=l
            print L
            rname=L[:-3]+match
            for R in dat.keys():
                if R==rname:
                    break

            if R==rname:
                cd=dat[L]-dat[R]
                cd[0]=dat[l][0]
                
                CD[L[:-3]]=cd
            else:
                print 'No complement found for %s' % L
    
    
    for (n,v) in CD.iteritems():
        for (i, y) in enumerate(v[1:]):
            label=n #+'_'+fields[i+1]
            if normalize==True:
                maxy=y.max()
            else:
                maxy=1
            plt.plot(v[0], y/maxy, label=label, lw=lw, **kwargs)

def split_string(s, widths=None):
    '''
    Splits a string into list of strings at the partition points provided by
    the sequence widths.    
    '''

    if widths is None:
        widths=[10,11,10,10,10,10]

    a=[]
    lastc=0
    for c in widths:
        a.append(s[lastc:lastc+c])
        lastc+=c

    return a

def clean_string(s):
    '''
    Remove illegal tokens from a string to return something appropriate for a python name
    '''
    for c in '<>=^+-(){}@#$%&*!?,.~':
        s=s.replace(c, '')

    return s        
