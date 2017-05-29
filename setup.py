#!/usr/bin/env python3

from glob import glob
from setuptools import find_packages, setup

setup(
    name='binf_utils',
    version='1.0',
    author='Robert Calef',
    author_email='robert.calef@gmail.com',
    packages=find_packages(),
    scripts=[script for script in glob("bin/*")],
    description='Various bioinformatics tools and utilities',
    install_requires=[
        'argh'
    ],
    zip_safe=True,
    include_package_data = True
)
