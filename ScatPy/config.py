# -*- coding: utf-8 -*-
"""
The package configuration settings.

All paths should be specified using unix style forward slashes. These will
be written to file in the correct style.

"""
import os.path
import ntpath
import posixpath

config={}

def set_config(fname=None):
    """
    Select which configuration profile to use for ScatPy

    :param fname: The name of the file that contains the configuration
                  If None then it try to load the default profile
                  
    Profiles are stored in files ending in .conf
    The search scheme is to first look for the file in the CWD, followed by
    the folder ~/.ScatPy/ and finally the subdiretory profiles/ relative to
    where the config.py module resides.
    """
    global config    
    
    if fname is None:
        fname = 'default.conf'
        
    if not fname.endswith('.conf'):
        fname += '.conf'
        
    pkg_path = os.path.dirname(__file__)

    for path in ['./', os.path.expanduser('~/.ScatPy/'), os.path.join(pkg_path, 'profiles')]:
        full_name=os.path.join(path, fname)
        if os.path.exists(full_name):            
            break
        else:
            full_name = None

    if full_name is None:
        raise(IOError('Could not find configuration profile'))

    execfile(full_name, config) 

    del config['__builtins__']
    config['file']=full_name

    # Associate the correct path style based on OS
    if config['os'].lower() == 'unix' or config['os'].lower() == 'mac':
        config['path']=posixpath
    elif config['os'].lower() == 'windows':
        config['path']=ntpath
    else:
        raise ValueError('Unknown OS: %s' % config['os'])
    

# When imported set the configuration to default
set_config(None)
