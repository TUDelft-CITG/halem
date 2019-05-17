#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for hydraulic_infrastructure_realisation.
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


setup(
    author="Pieter van Halem",
    author_email="pieter.vanhalem@vanoord.com",
    classifiers=[
        'Development Status :: 1 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="The Route optimization for dynamic flow field searches the optimal route for a given flow field.",
    entry_points={
        'console_scripts': [
            'halem=halem.cli:cli',
        ],
    },
    install_requires=requires,
    long_description="",  # README + '\n\n' + CHANGES,
    include_package_data=True,
    keywords="Route Optimization in Dynamic Flow Fields",
    name="halem",
    packages=find_packages(include=["halem"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=tests_require,
    url="https://github.com/TUDelft-CITG/Route_optimization_in_dynamic_currents",
    version="0.0.1",
    zip_safe=False,
)
