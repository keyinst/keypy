#!/usr/bin/env python

from distutils.core import setup

setup(name='keypy',
      version='1.0',
      description='The KEY Institute for Brain-Mind Research - EEG Analysis Toolbox',
      author='Patricia Milz',
      author_email='patricia.milz@key.uzh.ch',
      url='https://github.com/keyinst/keypy',
      license='GPLv3',
      packages=['keypy', 'keypy.microstates', 'keypy.preprocessing', 'keypy.signal'],
      install_requires = ["numpy", "scipy", "h5py"]
     )