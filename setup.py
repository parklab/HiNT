#!/home/sw229/miniconda2/bin/env python

#version: 2.0
#author: Su Wang
#Contact: wangsu0623@gmail.com

import os,sys,io
from distutils.core import setup, Extension
from setuptools import find_packages
from os import path

if sys.version < "3.5.0":
    print("Please use a Python with higher version than 3.5.0")
    sys.stderr.write("CRITICAL: Python version must be 3.5 or higher!\n")
    sys.exit(1)
    exit(1)

this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

#with io.open(path.join(this_directory, 'requirements.txt')) as f:
#    requires = f.read().splitlines()
#requires = [req.strip() for req in requires]

requires = [
'h5py>=2.8.0',
'numpy>=1.16.1',
'cooler>=0.8.2',
'pypairix>=0.3.0',
'pairtools>=0.2.2',
'scipy>=1.0.1',
'pandas>=0.23.0',
'scikit-learn>=0.19.1',
'multiprocess>=0.70.5',
'argparse>=1.1',
'pytabix'
]

def main():
    #compilemis()
    if not float(sys.version[:3])>=3.5:
        sys.stderr.write("CRITICAL: Python version must be greater than or equal to 3.5!\n")
        sys.exit(1)
    setup(name="HiNT-Package",
          version="2.0.7",
          description="HiNT -- HiC for copy number vairations and translocations detection ",
          long_description=long_description,
          long_description_content_type='text/markdown',
          author='Su Wang',
          author_email='wangsu0623@gmail.com',
          package_dir={'HiNT' : 'HiNT'},
          install_requires = requires,
          setup_requires = requires,
          #packages=['HiNT'],
          packages=find_packages(),
          package_data={'HiNT':['externalScripts/*']},
          include_package_data=True,
          scripts=['bin/hint'],
	  url="https://github.com/suwangbio/HiNT_py3",
          classifiers=[
            'Programming Language :: Python :: 3',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: MIT License',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX',
            'Topic :: Software Development',
            ],
          )


if __name__ == '__main__':
    main()
