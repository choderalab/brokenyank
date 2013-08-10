#!/usr/bin/env python

"""
setup.py - installer for YANK
"""
__author__ = "John D. Chodera"
__version__ = "1.0"

from distutils.core import setup
import sys
import os

MAJOR_VERSION_NUM = '1'
MINOR_VERSION_NUM = '0'
IS_RELEASED = False

def buildKeywordDictionary():
    keywords = dict()
    
    keywords['name'] = 'yank'
    keywords['version'] = '%s.%s' % (MAJOR_VERSION_NUM, MINOR_VERSION_NUM)
    keywords['author'] = 'John D. Chodera'
    keywords['author_email'] = 'jchodera@gmail.com'
    keywords['license'] = 'GNU General Public License v3 (GPLv3)'
    keywords['url'] = 'http://github.com/choderalab/yank'
    keywords['download_url'] = 'http://github.com/choderalab/yank'

    keywords['install_requires'] = ['simtk.openmm', 'simtk.unit']
    keywords['dependency_links'] = 

    keywords['classifiers'] = ['Development Status :: 4 - Beta',
                               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                               'Environment :: Console',
                               'Intended Audience :: Science/Research',
                               'Operating System :: POSIX',
                               'Programming Language :: Python']

    keywords['packages'] = ['yank']
    keywords['data_files'] = []
    keywords['package_data'] = {'yank' : ['data/*']}
    keywords['platforms'] = ['Linux', 'Mac OS X', 'Windows']
    keywords['description'] = "A framework for alchemical free energy calculations"
    keywords['long_description'] = """\
YANK is a framework for experimenting with GPU-accelerated alchemical free energy calculations.
"""
    
    return keywords
    
def main():
    if sys.version.info < (2, 6):
        raise Exception("YANK requires Python 2.6 or later.")

    keywords = buildKeywordDictionary()
    setup(keywords)

if __name__ == '__main__':
    main()

