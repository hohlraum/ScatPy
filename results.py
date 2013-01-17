# -*- coding: utf-8 -*-
"""@package ScatPy.results
For reading, manipulating and plotting the output files from DDSCAT

"""

from __future__ import division
import numpy as np
import os
import os.path
import posixpath
import glob
import zipfile
import copy

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, UnivariateSpline
from collections import OrderedDict

import utils


#class Results():
#    def __init__(self, **kwargs):
#        self.log=Log()
#        self.refresh()
#        self.mtable=MTable(**kwargs)
#        self.qtable=QTable(**kwargs)
#        self.qtable2=QTable2(**kwargs)
#
#    def refresh(self):        
#        self.log.refresh()
#        self.mtable.refresh()
#        self.qtable.refresh()
#        self.qtable2.refresh()
#
#    def set_folder(self, newfolder):
#        self.mtable.set_folder(newfolder)
#        self.qtable.set_folder(newfolder)        
#        self.qtable2.set_folder(newfolder)
#    

class Table(dict):
    """
    Base class for tables
    
    """

    def __init__(self):
        dict.__init__(self)
        
        self.x_field=None
        self.default_plot_fields=None


    def plot(self, fields=None, normalize=False, smooth=False, x_field=None,
             label=None, lw=2, **kwargs):
        """
        Plot the table contents.
        
        Optional Arguments:
            fields: list of field names to plot. If absent, plots the fields found
                    in self.default_plot_fields
            normalize: if True, set the maximum of all plots to 1
            smooth: if True, applies a spline smooth to the plot
            x_field: specify the name of the x_field. Default is self.x_field
            label: a string giving a specific label for the y-axis
        
        """
        if x_field is None:
            x_field=self.x_field

        if fields is None:
            fields=self.default_plot_fields
    
        ylbl=''
        for f in fields[0:]:
            ylbl+=f+', '
        ylbl=ylbl[:-2] #trim last ', '
        
        plt.ylabel(ylbl)
      
        for f in fields:
            y=self[f]
            if label is None:
                l_text=f
            else:
                l_text=label
            
            try:
                if l_text[0]=='_':
                    l_text=l_text[1:]
            except IndexError:
                pass
            
            if normalize:
                maxy=y.max()
            else:
                maxy=1.
            
            if smooth is False:
                plt.plot(self[x_field], y/maxy, label=l_text, lw=lw, **kwargs)
            else:
                if smooth is True:
                    smooth=0.03
                    
                xnew = np.linspace(self[x_field].min(), self[x_field].max(),300)
                spline=UnivariateSpline(self[x_field], y/maxy, s=smooth)
                ynew=spline(xnew)                
                plt.plot(xnew, ynew, label=l_text, lw=lw, **kwargs)

    def scale(self, c):
        """
        Scale all fields of table by value c
        """
        
        for (k,v) in self.iteritems():
            if k==self.x_field: continue
        
            v*=c

    def export(self, fname):
        """
        Export the table as an ascii file
        
        File is written to the current working directory.
        """

        if self.x_field is not None:
            table=self[self.x_field]
            hdr=self.x_field+'\t'
        else:
            table=np.array([])
            hdr=''

        for (k,v) in self.iteritems():
            if k==self.x_field: continue
        
            hdr+=k+'\t'
            table=np.vstack([table, v])

        with open(fname, 'wt') as f:
            f.write(hdr+'\n')        
            np.savetxt(f, table.transpose(), delimiter='\t')

    def copy(self):
        return copy.deepcopy(self)


    def smooth(self, width=None, pad_sample=5):
        '''
        Broaden a spectrum by convolving with  Gaussian
        
        Spectrum is smoothed in place.
        
        sig, signal
        w, FWHM
        
        '''
        if width is None:
            width=0.030
        
        x_field=self.x_field
        xstep=self[x_field][1]-self[x_field][0] 
        
        width_steps=width/xstep #the width in number of steps
        x=np.arange(-np.ceil(width_steps), +np.ceil(width_steps)+1)
        gauss=np.sqrt(4*np.log(2)/np.pi)/float(width_steps)*np.exp(-4*np.log(2)/float(width_steps)**2 * x**2)
        gauss/=gauss.sum()
    
        data=self[x_field]
        for (k,v) in self.iteritems():
            if k==x_field: continue
            
            data=v
            
            left_pad=np.ones(width_steps+1)*(data[:pad_sample].mean())
            right_pad=np.ones(width_steps+1)*(data[-pad_sample:].mean())
            data=np.concatenate((left_pad, data, right_pad))
                        
            cv=np.convolve(data, gauss, 'same')
        
            cv=cv[width_steps+1: width_steps+1+len(v)]
                
            self[k]=cv


    def interpolate(self, n_pts=200):
        '''
        Interpolates the data onto a new (denser) xrange.

        Data is modified in place.
        
        '''

        old_x=self[self.x_field]
        new_x=np.linspace(old_x[0], old_x[-1], n_pts)        

        for (k,v) in self.iteritems():
            
            data=v
            interp=interp1d(old_x, data, kind='cubic')
                                        
            self[k]=interp(new_x)


    def __repr__(self):
        return '<%s.%s object at %s>' % (self.__class__.__module__,
                                        self.__class__.__name__,
                                        hex(id(self)))    


class ResultTable(Table):
    """
    Base class for results tables read from file.        
    """
    def __init__(self, fname, hdr_len, c_width, folder=None, zfile=None):
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
        Table.__init__(self)        
        
        if folder is None:
            self.folder=''
        else:
            self.folder=folder
         
        self.zfile=zfile
        self.fname=fname
        self.hdr_len=hdr_len
        self.c_width=c_width
        self.header=None
        self.col_lbl=None
        self.refresh()
#
#        try:
#        except (IOError):
#            print 'Result Table %s not found' % utils.resolve_mat_file(os.path.join(self.folder, self.fname))
#            self.data=None        
#        except (ValueError):
#            print 'Result Table %s format not recognized (wrong colum format?)' % utils.resolve_mat_file(os.path.join(self.folder, self.fname))
#
#            self.data=None

    def refresh(self):
        """
        Refresh the data from the file
                
        """
        
        if self.zfile:
            with zipfile.ZipFile(os.path.join(self.folder, self.zfile)) as z:
                f=z.open(self.fname, 'rU')            
                self._load(f)
                
        else:            
            with open(os.path.join(self.folder, self.fname), 'Ur') as f:
                self._load(f)

    def _load(self, f):
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
            self[l]=d

    def set_folder(self, new_folder):
        """
        Change the working folder.        
        """
        self.folder=new_folder


    
class AVGTable(ResultTable):
    """
    A class for reading .avg files output by DDscat
    
    """
    def __init__(self, fname=None, folder=None, **kwargs):
        if fname==None:
            fname='w000r000.avg'
        
        if folder==None:
            folder=''
            
        print fname
        with open(os.path.join(folder, fname), 'Ur') as f:
            
            for i in range(37):
                line=f.readline()
                if line.find('1 incident polarizations') <> -1:
                    hdr=32
                    c_widths=[3, 6, 6, 10, 10, 11, 11]
                    plot_fields=['THETA', '<|f11|^2>']
                    break
                if line.find('2 incident polarizations') <> -1:
                    hdr=37
                    c_widths=[6, 7,9,12,12,11,11,11,11,11,11,11]
                    plot_fields=['S_11']
                    break
            
        ResultTable.__init__(self, fname, hdr, c_widths, **kwargs)
        self.x_field='THETA'        
        self.default_plot_fields=plot_fields
        if hdr==32:
            self.summary=AVGSummaryTable(self.header)
        else:
            self.summary=AVGSummaryTable(self.header)

class AVGSummaryTable(Table):
    """
    A class for reading the summary section of AVGTables.

    """      
    def __init__(self, fname=None, folder=None, npol=None, zfile=None, **kwargs):
        if fname is None:
            fname='w000r000.avg'
        self.fname=fname
        
        if folder is None:
            folder=''
        self.folder=folder
        self.zfile=zfile

        Table.__init__(self)

        if npol is not None:
            self.npol=npol
        else:
            if self.zfile:
                with zipfile.ZipFile(os.path.join(self.folder, self.zfile)) as z:
                    f=z.open(self.fname, 'rU')            
                    self._find_pol(f)
            else:
                with open(os.path.join(self.folder, self.fname), 'Ur') as f:
                    self._find_pol(f)
                
        self.refresh()

    def _find_pol(self, f):

        for i in range(37):
            line=f.readline()
            if line.find('1 incident polarizations') <> -1:
                self.npol=1
                break
            if line.find('2 incident polarizations') <> -1:
                self.npol=2
                break

    def refresh(self):
 
        if self.zfile:
            with zipfile.ZipFile(os.path.join(self.folder, self.zfile)) as z:
                f=z.open(self.fname, 'rU')
                if self.npol==1:
                    self._refresh_1pol(f)
                else:
                    self._refresh_2pol(f)
        else:
            with open(os.path.join(self.folder, self.fname), 'Ur') as f:
                if self.npol==1:
                    self._refresh_1pol(f)
                else:
                    self._refresh_2pol(f)
   

    def _refresh_1pol(self, f):
                    
        hdr=''

        for i in range(32):
            hdr+=f.readline()
            
        self.header=hdr
        hdr=hdr.splitlines()

        l=hdr[9].split()
        self.wave=float(l[1])

        l=hdr[18]
        self['Epol']=np.array(utils.str2complexV(l))
                
        l=hdr[28].split()
        self['Q_ext']=np.array(float(l[1]))
        self['Q_abs']=np.array(float(l[2]))
        self['Q_sca']=np.array(float(l[3]))
    
    def _refresh_2pol(self, f):
                    
        hdr=''
        for i in range(37):
            hdr+=f.readline()
            
        self.header=hdr
        hdr=hdr.splitlines()

        l=hdr[9].split()
        self.wave=float(l[1])

        l1=hdr[18]
        l2=hdr[19]
        self['Epol']=np.array([utils.str2complexV(l1),
                              utils.str2complexV(l2)])
        
        l1=hdr[28].split()
        l2=hdr[29].split()
        self['Q_ext']=np.array([float(l1[1]), float(l2[1])])
        self['Q_abs']=np.array([float(l1[2]), float(l2[2])])
        self['Q_sca']=np.array([float(l1[3]), float(l2[3])])
            
    def dichroism(self):
        
        if self.npol==1:
            raise (ValueError, 'Table has only one polarization state')

        Epol0=utils.pol2str(self['Epol'][0])
        
        if Epol0=='cL' or Epol0=='lH':
            a=0
            b=1
        elif Epol0=='cR' or Epol0=='lV':
            a=1
            b=0
        else:
            raise(ValueError, 'Can only handle dichroism for cL, cR, lH, or lV polarizations')
        
        Qext=self['Q_ext']
        Qabs=self['Q_abs']
        Qsca=self['Q_sca']
        return [self.wave, Qext[a]-Qext[b], Qabs[a]-Qext[b], Qsca[a]-Qsca[b]]


class QTable(ResultTable):
    """
    A class for reading qtable output files from DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='qtable'
        ResultTable.__init__(self, fname, 13, [10,11,11,11,11,12,11,11,6], **kwargs)
        self.x_field='wave'
        self.default_plot_fields=['Q_ext']


class QTable2(ResultTable):
    """
    A class for reading qtable2 output files from DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='qtable2'
        ResultTable.__init__(self, fname, 12, [10,11,12,12,12], **kwargs)
        self.x_field='wave'
        self.default_plot_fields=['Q_pha']


class MTable(ResultTable):
    """
    A class for reading mtables generated by DDscat
    
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='mtable'
        ResultTable.__init__(self, fname, 7, [10,11,10,10,10,10], **kwargs)
        self.x_field='wave(um)'
        self.default_plot_fields=['wave(um)', 'Re(m)', 'Im(m)']
        
class MInTable(ResultTable):
    """
    A class for reading material input files used by DDscat.
    
    Simple file names are resolved relative to the default material_library path.
    
    This class can be called to return an interpolated value at a requested
    wavelength. e.g.:
        M=MInTable('Gold.txt')
        refr_index=M(0.750)
    """
    def __init__(self, fname=None, **kwargs):
        if fname==None:
            fname='mintable'
        ResultTable.__init__(self, fname, 2, None, **kwargs)
        self.x_field='wave'
        self.default_plot_fields=['Re(m)', 'Im(m)']
    
    def refresh(self):
        fname=utils.resolve_mat_file(posixpath.join(self.folder, self.fname))

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

    def save(self, fname):
        
        with open(fname, 'wt') as f:
            f.write(self.header+'\n')
            f.write('\t'.join(self.col_lbl)+'\n')
            np.savetxt(f, self.data, fmt='%6.4f\t%4.2f\t%5.3f')

class ResultCollection(OrderedDict):

    def __init__(self):
        """
        A collection of several results files together in one object.
    
        Results are returned as a dictionary of dictionaries where key names correspond
        to the folder names and the values are dictionaries of the requested fields
        """
        super(ResultCollection, self).__init__()    

    def __repr__(self):
        return '<%s.%s object at %s>' % (self.__class__.__module__,
                                        self.__class__.__name__,
                                        hex(id(self)))    

    
    def plot(self, fields=None, normalize=False, lw=2, **kwargs):                
        """
        Plot all of the tables in the collectin in one plot.


        """        
        
        for (key, table) in self.iteritems():
#            if fields is None:
#                flds=table.default_plot_fields
#            else:
#                flds=fields

#            label=key #+'_'+fields[i+1]
            table.plot(fields, normalize=normalize, label=key, lw=lw, **kwargs)
#                plt.plot(v[0], y/maxy, label=label, lw=lw, **kwargs)

        plt.xlabel(table.x_field)
        
        ylbl=''
        if fields is None:
            flds=table.default_plot_fields
        else:
            flds=fields
            
        for f in flds:
            ylbl+=f+', '
        ylbl=ylbl[:-2]
        
        plt.ylabel(ylbl)

    def smooth(self, width=None, pad_sample=5):
        """


        """        
        
        for v in self.itervalues():
            v.smooth(width=width, pad_sample=pad_sample)


    def dichroism(self, fields=None):
        """
        Calculate the dichroism between pairs of spectra.
        
        Looks for spectra in folders that have names suffixed with _cL and _cR
        and calculates the difference.
        """        
        
        
        if fields==None:
            fields=['Q_ext']            
            
        CD=ResultCollection()
        for L in self.keys():
            if L.endswith(('LCP', '_cL')):
                print 'Dichroism from: '+L
                if L.endswith('LCP'):
                    match='RCP'
                else:
                    match='_cR'

                R=L[:-3]+match
    
                if R in self.keys():
                    cd=Table()
                    for f in fields:
                        cd[f]=self[L][f]-self[R][f]
                    x_field=self[L].x_field
                    cd[x_field]=self[L][x_field]
                    
                    cd.x_field=x_field
                    cd.default_plot_fields=fields
                else:
                    print 'No complement found for %s' % L

                CD[L[:-3]]=cd

        return CD
        

class FolderCollection(ResultCollection):
    def __init__(self, r_type=None, path=None, recurse=True):
        """
        A collection of several results files together in one object.
    
        By default all subfolders from the current directory are scanned for
        results files. If subfolders=False then only the current directory (or
        optionally the directory specified by path will be searched for files)
    
        Results are returned as a dictionary of dictionaries where key names correspond
        to the folder names and the values are dictionaries of the requested fields

        Arguments:
            r_type: a string denoting the type of result file to load from each folder
            path: the root directory whose subfolders will be read
            recurse: True recurses all subdirectories, False loads only from the specified (or working folder)
        """
        if r_type is None:
            r_type='qtable'
    
        if r_type.lower()=='mtable':
            rtable=MTable
        elif r_type.lower()=='qtable':
            rtable=QTable
        elif r_type.lower()=='qtable2':
            rtable=QTable2
            
        if path is None:
            path='.'
            
#        folders=[i for i in os.listdir(path) if os.path.isdir(i)]

        super(FolderCollection, self).__init__()    

        folders=[]    
        if recurse:
            for root, dirs, files in os.walk(path):
                folders.append(root)                
        else:
            folders=[i for i in os.listdir(path) if os.path.isdir(i)]

        for f in folders:
            try:
                f_key=os.path.normpath(f)
                self[f_key]=rtable(folder=f)
            except (IOError):
                pass

    
class FileCollection(ResultCollection):
    def __init__(self, r_type=None, path=None):
        """
        A collection of several results files.
    
        Arguments:
            r_type: a string denoting the type of result file to load from each folder
            path: the root directory whose subfolders will be read        
        """
        if r_type is None:
            r_type='avgsummary'
    
        if r_type.lower()=='avgtable':
            rtable=AVGTable
            flt="*.avg"
        elif r_type.lower()=='avgsummary':
            rtable=AVGSummaryTable
            flt="*.avg"
            self.dichroism=self._avgsum_dichroism
        elif r_type.lower()=='fmltable':
            rtable=FMLTable
            flt='*.fml'
        elif r_type.lower()=='scatable2':
            rtable=SCATable2
            flt='*.sca'
            
        if path is None:
            path='.'

        super(FileCollection, self).__init__()    

        for f in glob.glob(os.path.join(path, flt)):
            try:
                f_key=os.path.split(f)[1]
                self[f_key]=rtable(f)
            except (IOError):
                pass

    def _avgsum_dichroism(self, fields=None):
        """
        An alternate method of calculating dichroism used for AVGSummary collections.
        
        
        """
        if fields==None:
            fields=['Q_ext', 'Q_abs', 'Q_sca']


        table=[]                
        for (k,v) in self.iteritems():
            table.append(v.dichroism())

        CD=Table()
        CD['wave']=np.asarray([t[0] for t in table])
        if 'Q_ext' in fields:
            CD['CDext']=np.asarray([t[1] for t in table])
        if 'Q_abs' in fields:    
            CD['CDabs']=np.asarray([t[2] for t in table])
        if 'Q_sca' in fields:
            CD['CDsca']=np.asarray([t[3] for t in table])
        
        CD.x_field='wave'
        CD.default_plot_fields=['CDext']
        
        return CD

    def split(self):
        """
        Split the collection into two collections corresponding to the two polarizations
        
        """

        item=self.itervalues().next()
        pol1=item['Epol'][0]
        pol2=item['Epol'][1]

        table1=[]
        table2=[]                
        for (k,v) in self.iteritems():
            m=np.array([[v.wave]*2, v['Q_ext'], v['Q_abs'], v['Q_sca']])

            table1.append(m[:,0])
            table2.append(m[:,1])

        table1=np.asarray(table1)
        table2=np.asarray(table2)

        P1=Table()
        P2=Table()

        indices={'wave':0, 'Q_ext':1, 'Q_abs':2, 'Q_sca':3}        
        for k in item.keys()+['wave']:
            if k=='Epol': continue
            idx=indices[k]
            P1[k]=table1[:,idx]
            P2[k]=table2[:,idx]
        
        P1.Epol=pol1
        P2.Epol=pol2
        
        P1.x_field='wave'
        P2.x_field='wave'
        
        P1.default_plot_fields=['Q_ext']
        P2.default_plot_fields=['Q_ext']

        C=ResultCollection()
        C['_'+utils.pol2str(P1.Epol)]=P1
        C['_'+utils.pol2str(P2.Epol)]=P2
        

        return C      
        

    def dichroism(self, fields=None):
        
        if fields==None:
            fields=['CDext']
                       
        CD=ResultCollection()
        for l in self.keys():
            print 'Dichroism from: '+l
            cd=self[l].summary.dichroism()
            CD[L[:-3]]=cd

        cd.default_plot_fields=fields
        return CD

class ZipCollection(FileCollection):
    def __init__(self, r_type=None, folder=None):
        """
        A collection of several results files archived into a single zip file.
    
        Arguments:
            r_type: a string denoting the type of result file to load from each folder
            folder: the root directory whose subfolders will be read        
        """
        if r_type is None:
            r_type='avgsummary'
    
        if r_type.lower()=='avgtable':
            rtable=AVGTable
            zname="all_avg.zip"
        elif r_type.lower()=='avgsummary':
            rtable=AVGSummaryTable
            zname="all_avg.zip"
            self.dichroism=self._avgsum_dichroism
        elif r_type.lower()=='fmltable':
            rtable=FMLTable
            zname='all_fml.zip'
        elif r_type.lower()=='scatable2':
            rtable=SCATable2
            zname='all_sca.zip'
        
        self.folder=folder
        if folder is None:
            folder=''
        
        ResultCollection.__init__(self)

#        super(ZipCollection, self).__init__()    

#        if not os.path.exists(os.path.join(folder, zname)):
#            raise (IOError, 'Archive %s not found'% os.path.join(folder, zname))
 #       else:
        with zipfile.ZipFile(os.path.join(folder, zname), 'r') as z:
            names=z.namelist()

        for n in names:
            self[n]=rtable(n, folder=self.folder, zfile=zname)



def split_string(s, widths=None):
    '''
    Splits a string into list of strings at the partition points provided by
    the sequence widths.
    
    Required to parse ddscat results files which are fixed width without whitespace
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
    Remove illegal tokens from a string to return something appropriate for a Python name
    '''
    for c in '<>=^+-(){}@#$%&*!?,.~':
        s=s.replace(c, '')

    return s        



### ===============================================================================
### DEPRECATED

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
        


