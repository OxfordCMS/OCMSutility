import sysconfig
import sys
import os
import subprocess
import re
from setuptools import setup, find_packages

setup(
    # package information
    name='ocms',
    version="0.0",
    description='ocms : Oxford Centre for Microbiome Studies apps',
    author='Nicholas Ilott, Sandi Yen, Jethro Johnson',
    license="MIT",
    platforms=["any"],
    keywords="microbiome, metagenomics, genomics",
    packages=find_packages() + find_packages("./Py_utility"),
    entry_points={
        'console_scripts': ['ocms = Py_utility.ocms:main']
    }
)
