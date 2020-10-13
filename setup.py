#!/usr/bin/env python
from setuptools import setup, find_packages

version = '0.0.1'

setup(
    name='multiqc_alfa',
    version=version,
    author='Krish Agarwal',
    author_email='akrish136@gmail.com',
    description="MultiQC plugin for the Zavolab \
        @ University of Basel, Switzerland",
    long_description=__doc__,
    keywords='bioinformatics',
    url='https://github.com/zavolanlab/multiqc-plugins',
    download_url='',  # after releasing a version
    license='',  # to add
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'multiqc'
    ],
    entry_points={
        'multiqc.modules.v1': [
            'alfa = multiqc_alfa.modules.alfa:MultiqcModule',
        ],
    },
    classifiers=[
        'Development Status :: Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: ',  # to add
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)