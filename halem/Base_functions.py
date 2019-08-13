import halem.Mesh_maker as Mesh_maker
import halem.Functions as Functions
import halem.Calc_path as Calc_path

import numpy as np
import pickle
import matplotlib.pyplot as plt
import datetime, time
from datetime import datetime


def save_object(obj, filename):
    """This function can save the roadmap using pickle"""
    with open(filename, "wb") as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def plot_timeseries2(path, time, Roadmap, Color="r"):
    """ This function can plot the time series for the route 
    and shows a contourplot of the unsaiable areas of that route 
    
    path:       lon, lat coordinates of the route. This is in the format 
                of the output from halem.Base_functions.HALEM_func[0]
    time:       time series of the path. This is in the format 
                of the output from halem.Base_functions.HALEM_func[1]
    Roadmap:    Roadmap that is used to calculate the route. 
    Color:      Color of the plot of the time series,
                type sting, with matplotlib color"""

    dist = []
    TT = []
    D = 0
    for i in range(len(path) - 1):
        nx = (
            (Roadmap.nodes[:, 0] - path[i, 1]) ** 2
            + (Roadmap.nodes[:, 1] - path[i, 0]) ** 2
        ) ** 0.5
        idx = np.argwhere(nx == nx.min())[0][0]
        D = D + Functions.haversine(
            (path[i, 1], path[i, 0]), (path[i + 1, 1], path[i + 1, 0])
        )
        dist.append(D)
        T = Roadmap.mask[idx]
        TT.append(T)
    TT = np.array(TT)
    dist = np.array(dist)
    if Roadmap.repeat == True:
        k = Calc_path.find_k_repeat(time[0], Roadmap.t)
        plt.plot(dist, (time[:-1] - time[0]) / 3600, color=Color, label="s/t route")
        cval = np.arange(0, 1.1, 0.5)
        plt.contourf(
            dist,
            (Roadmap.t - Roadmap.t[k]) / 3600,
            np.transpose(TT),
            cval,
            colors=("cornflowerblue", "sandybrown"),
        )
        plt.contourf(
            dist,
            (Roadmap.t - Roadmap.t[k] + Roadmap.t[-1]) / 3600,
            np.transpose(TT),
            cval,
            colors=("cornflowerblue", "sandybrown"),
        )

        plt.colorbar(label="maks file, 0 = False, 1 = True")
        plt.xlabel("traveled distance [m]")
        plt.ylabel("time [h]")
        plt.ylim(0, (time[-1] - time[0]) / 3600 * 1.2)
        plt.legend(loc="best")

    else:
        plt.plot(dist, (time[:-1] - time[0]) / 3600, color=Color, label="s/t route")
        cval = np.arange(0, 1.1, 0.5)
        plt.contourf(
            dist,
            (Roadmap.t - time[0]) / 3600,
            np.transpose(TT),
            cval,
            colors=("cornflowerblue", "sandybrown"),
        )
        plt.colorbar(label="maks file, 0 = False, 1 = True")
        plt.ylim(0, (time[-1] - time[0]) / 3600 * 1.2)
        plt.xlabel("traveled distance [m]")
        plt.ylabel("time [h]")
        plt.legend(loc="best")


def HALEM_func(start, stop, t0, vmax, Roadmap, costfunction):
    """ Base of the oe lne functions halem.Base_functions.HALEM_time,
    halem.Base_functions.HALEM_cost, halem.Base_functions.HALEM_space, 
    halem.Base_functions.HALEM_co2. This function takes the pre-processing 
    file, start location, stop location, departure time, and sailing velocity 
    and returns the optimized route. 

    start:          (lon, lat) coordinates of the start location
    stop:           (lon, lat) coordinates of the destination location
    t0:             ('day'/'month'/'year' 'hour':'minute':'seconds') 
                    string that indcates the departure time
    vmax:           (N (rows) * M (columns)) numpy array that indicates the sailing 
                    velocity in deep water. For which N is the number of discretisations
                    in the load factor, and M is the number of discretisations in the 
                    dynamic sailing velocity

                    For the optimization type cost and co2 N must be larger or equal to 2.

    Roadmap:        Preprocessing file that contains the hydrodynamic properties, 
                    and vesssel parameters. Output of the function 
                    halem.Mesh_maker.Graph_flow_model
    costfunction    Costfunction of the route optimization.
                    Roadmap.weight_time returns fastest route
                    Roadmap.weight_space returns shortest route
                    Roadmap.weight_cost returns cheapest route
                    Roadmap.weight_co2 retruns least pollutant route
    """

    start = start[::-1]
    stop = stop[::-1]

    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class graph_functions_time:
        weights = costfunction[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path = Roadmap.nodes[np.array(route.route[:, 0], dtype=int)]
    time = route.route[:, 1]

    dist = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + Functions.haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
    dist = np.array(dist)
    return path[:, ::-1], time, dist


def HALEM_time(start, stop, t0, vmax, Roadmap):
    """Implementation of the function halem.Base_functions.HALEM_func() for the fastest route."""
    costfunction = Roadmap.weight_time
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def HALEM_space(start, stop, t0, vmax, Roadmap):
    """Implementation of the function halem.Base_functions.HALEM_func() for the shortest route."""
    costfunction = Roadmap.weight_space
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def HALEM_cost(start, stop, t0, vmax, Roadmap):
    """Implementation of the function halem.Base_functions.HALEM_func() for the cheapest route."""
    costfunction = Roadmap.weight_cost
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def HALEM_co2(start, stop, t0, vmax, Roadmap):
    """Implementation of the function halem.Base_functions.HALEM_func() for the least pollutant route."""
    costfunction = Roadmap.weight_co2
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)
