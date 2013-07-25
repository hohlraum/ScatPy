# -*- coding: utf-8 -*-

#from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name='ScatPy',
    version='0.1.0',
    author='Andrew G. Mark',
    author_email='mark@is.mpg.de',
    packages=['ScatPy'],
    url='https://github.com/hohlraum/ScatPy',
    platforms = 'All',
    license='GNU GPLv3',
    description='A Python package for setting up DDSCAT jobs and analysing the results.',
    long_description=open('README.txt').read(),
    classifiers = ['Development Status :: 3 - Alpha',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        ]        
)