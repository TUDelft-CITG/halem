import halem.Mesh_maker as Mesh_maker
import halem.Functions as Functions
import halem.Calc_path as Calc_path
import halem.Flow_class as Flow_class

import math
import numpy as np
from numpy import ma
import pickle

import os
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib import animation
import datetime, time
from datetime import datetime


def save_object(obj, filename):
    with open(filename, "wb") as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def plot_timeseries(route, Roadmap):
    dist = []
    TT = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + Mesh_maker.haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
        T = Roadmap.mask[int(route.route[i, 0])]
        TT.append(T)
    TT = np.array(TT)
    dist = np.array(dist)

    plt.plot(dist, route.route[:-1, 1], "ro")
    cval = np.arange(0, 1.1, 0.5)
    plt.contourf(dist, Roadmap.t, np.transpose(TT), cval)
    plt.colorbar()
    plt.xlabel("traveled distance [m]")
    plt.ylabel("time [s]")

    plt.ylim(route.route[0, 1] - 2000, route.route[-1, 1] + 2000)


def plot_timeseries2(path, time, Roadmap, Color="r"):
    dist = []
    TT = []
    D = 0
    for i in range(len(path) - 1):
        nx = (
            (Roadmap.nodes[:, 0] - path[i, 1]) ** 2
            + (Roadmap.nodes[:, 1] - path[i, 0]) ** 2
        ) ** 0.5
        idx = np.argwhere(nx == nx.min())[0][0]
        D = D + Mesh_maker.haversine(
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


def HALEM_time(start, stop, t0, vmax, Roadmap):
    start = start[::-1]
    stop = stop[::-1]

    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])
    # arg_vship = arg_vship if arg_vship == 0 or vmax < vvmax[-1]  else arg_vship - 1

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_time[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path = Roadmap.nodes[np.array(route.route[:, 0], dtype=int)]
    time = route.route[:, 1]

    dist = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + Mesh_maker.haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
    dist = np.array(dist)
    return path[:, ::-1], time, dist


def HALEM_space(start, stop, t0, vmax, Roadmap):
    start = start[::-1]
    stop = stop[::-1]

    vvmax = Roadmap.vship[:, -1]

    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])
    # arg_vship = arg_vship if arg_vship == 0 or vmax < vvmax[-1]  else arg_vship - 1

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_space[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path = Roadmap.nodes[np.array(route.route[:, 0], dtype=int)]
    time = route.route[:, 1]

    dist = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + Mesh_maker.haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
    dist = np.array(dist)
    return path[:, ::-1], time, dist


def HALEM_cost(start, stop, t0, vmax, Roadmap):
    start = start[::-1]
    stop = stop[::-1]

    vvmax = Roadmap.vship[:, -1]

    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])
    # arg_vship = arg_vship if arg_vship == 0 or vmax < vvmax[-1]  else arg_vship - 1

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_cost[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path = Roadmap.nodes[np.array(route.route[:, 0], dtype=int)]
    time = route.route[:, 1]

    dist = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + Mesh_maker.haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
    dist = np.array(dist)
    return path[:, ::-1], time, dist


def HALEM_co2(start, stop, t0, vmax, Roadmap):
    start = start[::-1]
    stop = stop[::-1]

    vvmax = Roadmap.vship[:, -1]

    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])
    # arg_vship = arg_vship if arg_vship == 0 or vmax < vvmax[-1]  else arg_vship - 1

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_co2[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path = Roadmap.nodes[np.array(route.route[:, 0], dtype=int)]
    time = route.route[:, 1]

    dist = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + Mesh_maker.haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
    dist = np.array(dist)
    return path[:, ::-1], time, dist
