import sys

from halem.functions import (
    HALEM_co2,
    HALEM_cost,
    HALEM_func,
    HALEM_space,
    HALEM_time,
    plot_timeseries,
)
from halem.roadmap import BaseRoadmap

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "python-package-example"
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError


__all__ = [
    "BaseRoadmap",
    "plot_timeseries",
    "HALEM_func",
    "HALEM_time",
    "HALEM_space",
    "HALEM_cost",
    "HALEM_co2",
]
