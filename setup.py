#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(name='seekCRIT',
      version='1.0.0',
      description='seek for circular RNA in transcriptome (identifies deferentially expressed circRNAs between two samples)',
      author='Bioinformatics Lab, University of Louisville, Kentucky Biomedical Research Infrastructure Network (KBRIN)',
      author_email='uoflbioinformatics@gmail.com',
      maintainer='Bioinformatics Lab, University of Louisville, Kentucky Biomedical Research Infrastructure Network (KBRIN)',
      maintainer_email='uoflbioinformatics@gmail.com',
      url='https://github.com/UofLBioinformatics/seekCRIT',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'pysam>=0.9.1.4',
          'numpy>=1.11.2',
          'scipy',
          'fisher',
          'mne',
          'docopt',
      ],
      entry_points={
          'console_scripts': [
              'seekCRIT.py=seekCRIT.seekCRIT:main',

          ],
      },
      )
