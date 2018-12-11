#!/home/sw229/miniconda2/bin/env python

#version: 1.0
#author: Su Wang
#Contact: wangsu0623@gmail.com

import os
import sys
from distutils.core import setup, Extension
from pkg_resources import resource_filename
from subprocess import call as subpcall
from setuptools import find_packages

if sys.version < "2.6.0" or sys.version > "2.8.0":
    print "Please use a Python with higher version than 2.6.0"
    sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
    sys.exit(1)
    exit(1)

def run_cmd(command):
    subpcall (command, shell = True)

def main():
    #compilemis()
    if not float(sys.version[:3])>=2.5:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 2.5! python 2.6.1 or newer is recommended!\n")
        sys.exit(1)
    setup(name="HiNT-Package",
          version="1.0",
          description="HiNT -- HiC for copy number vairations and translocations detection ",
          author='Su Wang',
          author_email='wangsu0623@gmail.com',
          package_dir={'HiNT' : 'HiNT'},
          install_requires=['argparse','numpy'],
          packages=['HiNT'],
          scripts=['bin/hint'],
          package_data={'HiNT':['references/*','backgroundMatrices/100kb/*','backgroundMatrices/1Mb/*','externalScripts/*']},

          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: Software Development',
            ],
          )


if __name__ == '__main__':
    main()
