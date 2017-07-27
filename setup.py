#!/usr/bin/env python

from __future__ import print_function, absolute_import
from setuptools import setup, find_packages

# Read the version from a file to be sure it is consistent with the version
# in the package.
with open('martinize/VERSION.txt') as infile:
    version = infile.readline().strip()

setup(
    name='martinize',
    version=version,

    description="A versatile tool for coarse graining molecular models.",

    url='https://github.com/Tsjerk/Martinize',

    # Author details
    author='Tsjerk A. Wassenaar',

    license='GPLv2',

    classifiers=[
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',

        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],

    install_requires=['numpy', 'simopt'],

    tests_requires=['nose'],

    packages=find_packages(),
    package_data={'martinize': ['VERSION.txt', 'quotes.txt']},

    entry_points={
        'console_scripts': [
            'martinize = martinize.cli:cli',
        ],
    },

)
