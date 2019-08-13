"""Top-level package for Halem."""

__author__ = """Pieter van Halem"""
__email__ = "pieter.vanhalem@vanoord.com"
__version__ = "0.3.0"


# from halem.Base_functions import *
from pkg_resources import get_distribution, DistributionNotFound

try:
    dist_name = "halem"
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound
