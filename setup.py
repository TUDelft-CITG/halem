#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for Route_Optimization_in_Dynamic_Flow_Fields.
    Use setup.cfg to configure your project.

    This file was generated with PyScaffold 3.1.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
import sys

from pkg_resources import require, VersionConflict
from setuptools import setup, find_packages

try:
    require('setuptools>=38.3')
except VersionConflict:
    print("Error: version of setuptools is too old (<38.3)!")
    sys.exit(1)

requires = [
    "numpy",
    "simpy",
    "networkx",
    "shapely",
    "scipy",
    "netCDF4",
    "click",
    "matplotlib",
    "IPython",
    "geopy",
    "pint",
    "pyproj",
    "plotly",
    "simplekml",
    "nose",
    "Flask",
    "Flask-cors",
    "Dill>=0.2.8",
]

setup_requirements = [
    "pytest-runner",
]

tests_require = [
    "pytest",
    "pytest-cov",
    "pytest-timeout",
    "pytest-datadir"
]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="halem-pkg-PietervanHalem",
    version="0.1.0",
    author="Pieter van Halem",
    author_email="pieter.vanhalem@vanoord.com",
    classifiers=[
        'Development Status :: 1 - Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    description="The Route optimization for dynamic flow field searches the optimal route for a given flow field.",
    entry_points={
        'console_scripts': [
            'halem=halem.cli:cli',
        ],
    },
    install_requires=requires,
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="Route Optimization in Dynamic Flow Fields",
    packages=find_packages(include=["halem"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=tests_require,
    url="https://github.com/TUDelft-CITG/Route_optimization_in_dynamic_currents",
    zip_safe=False,
)
