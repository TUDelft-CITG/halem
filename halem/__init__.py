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

from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = 'halem'
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound

def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def plot_timeseries(route, Roadmap):
    dist = []
    TT = []
    D = 0
    for i in range(route.route[:,0].shape[0] -1 ):
        D =D + Mesh_maker.haversine((route.y_route[i], route.x_route[i]), (route.y_route[i +1], route.x_route[i+1]))
        dist.append(D)
        T = Roadmap.mask[int(route.route[i,0])]
        TT.append(T)
    TT = np.array(TT)
    dist = np.array(dist)

    plt.plot(dist,route.route[:-1,1], 'ro')
    cval = np.arange(0,1.1, 0.5)
    plt.contourf(dist, Roadmap.t, np.transpose(TT), cval)
    plt.colorbar()
    plt.xlabel('traveled distance [m]')
    plt.ylabel('time [s]')

    plt.ylim(route.route[0,1]-2000,route.route[-1,1]+2000)

def HALEM_time(start, stop, t0, vmax, Roadmap):
        start = start[::-1]
        stop = stop[::-1]

        vvmax = Roadmap.vship[:,-1]

        vv= np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])
        class graph_functions_time:
            function_type = "time optimalisation"
            weights = Roadmap.weight_time[arg_vship].weights
            time = Roadmap.weight_time[arg_vship].weights
            vship = Roadmap.vship[arg_vship]

        route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
        path = Roadmap.nodes[np.array(route.route[:,0], dtype=int)]
        time = route.route[:,1]

        dist = []
        D = 0
        for i in range(route.route[:,0].shape[0] -1 ):
                D =D + Mesh_maker.haversine((route.y_route[i], route.x_route[i]), (route.y_route[i +1], route.x_route[i+1]))
                dist.append(D)
        dist = np.array(dist)
        return path[:,::-1], time, dist