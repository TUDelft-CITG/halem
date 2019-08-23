# -*- coding: utf-8 -*-

""" To integrate Halem into the discrete event simulation package OpenCLSim. """

# package(s) related to the simulation
import simpy
import networkx as nx

# spatial libraries
import pyproj
import shapely.geometry

# additional packages
import math
import datetime, time
import copy
import numpy as np
import pandas as pd

class Routeable(core.Movable):
    """
    Movement travels trough a Graph. When optimize_route == True the halem package optimizes the route for different cost functions
    if the optimization_type == 'time' then the fastest path is calculated by halem
    if the optimization_type == 'cost' then the cheapest path is calculated by halem
    if the optimization_type == 'space' then the shortest path is calculated by halem
    Optimize_Route == False the route is determind by v = s/t 
    """

    def __init__(
        self,
        route,
        optimize_route=False,
        optimization_type="time",
        loadfactors=None,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        """Initialization"""
        self.route = route
        self.optimize_route = optimize_route
        self.optimization_type = optimization_type
        self.loadfactors = loadfactors

        if self.optimize_route == True:

            if self.optimization_type == "time":
                self.optimization_func = halem.HALEM_time
            elif self.optimization_type == "space":
                self.optimization_func = halem.HALEM_space
            elif self.optimization_type == "cost":
                self.optimization_func = halem.HALEM_cost
            elif self.optimization_type == "co2":
                self.optimization_func = halem.HALEM_co2
            else:
                print("No known optimization method selected")

class HasDepthRestriction:
    """HasDepthRestriction class

    Used to add depth limits to vessels - improvement of the OpenCLSim class
    draught: should be a lambda function with input variable container.volume
    ukc: list with ukc, corresponding to wave_heights

    filling: filling degree [%]
    min_filling: minimal filling degree [%]
    max_filling: max filling degree [%]
    """

    def check_optimal_filling_Roadmap(self, loader, unloader, origin, destination):
        orig = shapely.geometry.asShape(origin.geometry)
        dest = shapely.geometry.asShape(destination.geometry)
        start = (orig.x, orig.y)
        stop = (dest.x, dest.y)
        prods = []
        tQQ = []

        for i in range(len(self.loadfactors)):
            # Departure time sailing empty
            t0 = datetime.datetime.fromtimestamp(self.env.now).strftime(
                "%d/%m/%Y %H:%M:%S"
            )

            # Duration of sailing empty
            _, TT, _ = self.optimization_func(
                stop, start, t0, self.env.Roadmap.vship[0, -1], self.env.Roadmap
            )

            # Duration of sailing empty + loading
            duration_dredging = self.loading(
                None, None, self.container.capacity * self.loadfactors[i]
            )

            TTT = self.env.now + (TT[-1] - TT[0]) + duration_dredging
            duration_sailing_empty = TT[-1] - TT[0]

            # Departure time sailing full
            t0 = datetime.datetime.fromtimestamp(TTT).strftime("%d/%m/%Y %H:%M:%S")

            # Duration of sailing full
            _, TT, _ = self.optimization_func(
                start, stop, t0, self.env.Roadmap.vship[i, -1], self.env.Roadmap
            )

            # Duration of sailing empty + loading + sailing full + unloading
            duration_unloading = self.unloading(
                None, None, self.container.capacity * self.loadfactors[i]
            )
            TTT += (TT[-1] - TT[0]) + duration_unloading
            duration_sailing_full = TT[-1] - TT[0]

            # Determine production
            prod = (self.loadfactors[i] * self.container.capacity) / (
                duration_sailing_empty
                + duration_sailing_full
                + duration_dredging
                + duration_unloading
            )
            prods.append(prod)
            tQQ.append(t0)

        optimal_loadfactor = self.loadfactors[np.argwhere(prods == max(prods))[0, 0]]
        return optimal_loadfactor