#!/usr/bin/env python

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='mresvgene',
    version='1.1.0',
    description='Multi-resolution V-gene discovery library',
    long_description=open('README.rst').read(),
    url='https://github.com/pypa/X',

    author='ImmunoInformatics Group, UVigo'
           'SERGAS, Hospital Meixoeiro, Vigo',
    author_email='dnolivieri@gmail.com',
    license='MIT',
    classifiers = ['Development Status :: 3 - Alpha',
               'Programming Language :: Python :: 2.7',
               'Environment :: Console',
               'License :: OSI Approved :: BSD License',
               'License :: Free for non-commercial use',
               'Topic :: Scientific/Engineering :: Bio-Informatics',
               'Topic :: Scientific/Engineering :: Artificial Intelligence'],

    keywords='genes-finding immunogenetics mhc',
    packages=[ 'mhcfinder' ],
    install_requires=[
        'Banyan>=0.1.5',
        'pypro>=0.1.0',
        'numpy>=1.9.3',
        'matplotlib>=1.4.3',
        'biopython>=1.65',
        'scipy>=0.15.1',
        'numpy>=1.9.3'
        ], 
    

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'contig_data': ['*.json', 'contigs', 'fastas'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    data_files=[('contig_data', ['contig_data/*'])],


    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': ['bootstrap=mresvgene.mresvgene_bootstrap:main'],
    },
)
