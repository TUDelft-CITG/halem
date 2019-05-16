import RODFF.Mesh_maker as Mesh_maker
import RODFF.Functions as Functions
import RODFF.Calc_path as Calc_path
import RODFF.Flow_class as Flow_class

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
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        
def printtime(t):
    print("Route completed in",int(t/3600), "hour", 
      int((t-int(t/3600)*3600)/60), "minutes and", 
      np.round(t -int(t/3600)*3600- int((t-int(t/3600)*3600)/60)*60, 2), "seconds")
    return

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

def RODFF_time(start, stop, t0, vmax, Roadmap):
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

        path2 = np.array(route.route[:,0], dtype=int)
        path = np.zeros((len(route.route[:,0]),2))
        path[:,0] = Roadmap.nodes[path2][:,1]
        path[:,1] = Roadmap.nodes[path2][:,0]

        time = route.route[:,1]

        dist = []
        D = 0
        for i in range(route.route[:,0].shape[0] -1 ):
                D =D + Mesh_maker.haversine((route.y_route[i], route.x_route[i]), (route.y_route[i +1], route.x_route[i+1]))
                dist.append(D)
        dist = np.array(dist)
        return path, time, dist

