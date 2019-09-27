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
    "matplotlib",
    "IPython",
    "geopy",
    "pyproj",
    "sphinx_rtd_theme",
    "simplekml",
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
    name="halem",
    version="0.3.0",
    author="Pieter van Halem",
    author_email="pieter.vanhalem@vanoord.com",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    description="The Route optimization for dynamic flow field searches the optimal route for a given flow field.",
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
