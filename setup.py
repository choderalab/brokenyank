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

    keywords['requires'] = ['numpy', 'scipy', 'pymbar'] # automatically install pymbar from PyPI if not present
                                      # TODO: How can we require OpenMM?

    #keywords['dependency_links'] = ['http://github.com/SimTK/openmm/tarball/master#egg=openmm-5.2.egg'] # Specify OpenMM egg

    keywords['classifiers'] = ['Development Status :: 4 - Beta',
                               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                               'Environment :: Console',
                               'Intended Audience :: Science/Research',
                               'Operating System :: POSIX',
                               'Programming Language :: Python']

    keywords['packages'] = ['yank']
    keywords['package_dir'] = { 'yank' : 'src/yank' }
    keywords['data_files'] = []
    keywords['package_data'] = {'yank' : ['data/*/*']}
    keywords['platforms'] = ['Linux', 'Mac OS X', 'Windows']
    keywords['description'] = "A framework for alchemical free energy calculations"
    keywords['long_description'] = """\
YANK is a framework for experimenting with GPU-accelerated alchemical free energy calculations.
"""
    keywords['scripts'] = ['tools/analyze.py', 'tools/analyzeall.py', 'yank/yank.py']

    return keywords
    
def main():
    if sys.version.info < (2, 7):
        raise Exception("YANK requires Python 2.7 or later.")

    keywords = buildKeywordDictionary()
    setup(keywords)

if __name__ == '__main__':
    main()

