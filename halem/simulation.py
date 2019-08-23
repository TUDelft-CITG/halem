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

# HALEM
import halem.Base_functions as halem

# OpenCLSim
import openclsim.core as core


class Routeable(core.Routeable):
    """
    Movement travels trough a Graph. When optimize_route == True the halem package optimizes the route for different cost functions
    if the optimization_type == 'time' then the fastest path is calculated by halem
    if the optimization_type == 'cost' then the cheapest path is calculated by halem
    if the optimization_type == 'space' then the shortest path is calculated by halem
    Optimize_Route == False the route is determind by v = s/t 
    """

    def __init__(
        self,
        optimize_route=False,
        optimization_type="time",
        loadfactors=None,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        """Initialization"""
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

    def sailing_duration(self, origin, destination, engine_order, verbose=True):
        """ Determine the sailing duration based on the properties of the sailing route """

        # If travelling on route is required, assert environment has a graph
        assert hasattr(self.env, "FG")
        assert hasattr(self.env, "Roadmap")
        assert self.optimize_route == True

        geom = nx.get_node_attributes(self.env.FG, "geometry")
        route = self.determine_route(origin, destination)

        # Determine the duration of the following the route
        duration = 0

        for i, _ in enumerate(route):
            if i + 1 != len(route):
                orig = shapely.geometry.asShape(geom[route[i]])
                dest = shapely.geometry.asShape(geom[route[i + 1]])
                time = datetime.datetime.fromtimestamp(self.env.now)

                path, time, _ = self.optimization_func(
                    start=(orig.x, orig.y),
                    stop=(dest.x, dest.y),
                    t0=time.strftime("%d/%m/%Y %H:%M:%S"),
                    vmax=self.determine_speed(route[i], route[i + 1]),
                    Roadmap=self.env.Roadmap,
                )

                for i in range(path.shape[0]):
                    if i + 1 != path.shape[0] and verbose:
                        self.log_entry(
                            "Sailing",
                            time[i + 1],
                            0,
                            shapely.geometry.Point(path[i]),
                            self.ActivityID,
                        )

                duration += time[-1] - time[0]

        return duration


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
