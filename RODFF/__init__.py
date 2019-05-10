import RODFF.Mesh_maker as Mesh_maker
import RODFF.Functions as Functions
import RODFF.Calc_path as Calc_path

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

def RODFF_time(start, stop, t0, vship, Roadmap):
        start = start[::-1]
        stop = stop[::-1]

        vv= np.abs(Roadmap.vship - vship)
        arg_vship = int(np.argwhere(vv == vv.min())[0])
        class graph_functions_time:
                function_type = "time optimalisation"
                weights = Roadmap.weight_time[arg_vship].weights
                time = Roadmap.weight_time[arg_vship].weights
               
        route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)

        path2 = Roadmap.nodes[list(map(int, route.route[:,0]))]
        path = np.zeros((path2.shape))
        path[:,0] = path2[:,1]
        path[:,1] = path2[:,0]

        time = route.route[:,1]

        dist = []
        D = 0
        for i in range(route.route[:,0].shape[0] -1 ):
                D =D + Mesh_maker.haversine((route.y_route[i], route.x_route[i]), (route.y_route[i +1], route.x_route[i+1]))
                dist.append(D)
        dist = np.array(dist)


        return path, time, dist

def RODFF_space(start, stop, t0, vship, Roadmap):
        start = start[::-1]
        stop = stop[::-1]

        vv= np.abs(Roadmap.vship - vship)
        arg_vship = int(np.argwhere(vv == vv.min())[0])

        class graph_functions_time:
                function_type = "time optimalisation"
                weights = Roadmap.weight_space[arg_vship].weights
                time = Roadmap.weight_time[arg_vship].weights
               
        route = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)

        path2 = Roadmap.nodes[list(map(int, route.route[:,0]))]
        path = np.zeros((path2.shape))
        path[:,0] = path2[:,1]
        path[:,1] = path2[:,0]

        time = route.route[:,1]

        dist = []
        D = 0
        for i in range(route.route[:,0].shape[0] -1 ):
                D =D + Mesh_maker.haversine((route.y_route[i], route.x_route[i]), (route.y_route[i +1], route.x_route[i+1]))
                dist.append(D)
        dist = np.array(dist)


        return path, time, dist


def get_sailtime(start, stop, t0, vship, Roadmap):
    start = start[::-1]
    stop = stop[::-1]

    t0 = str(datetime.fromtimestamp(t0).strftime("%d/%m/%Y %H:%M:%S"))
    #print(t0, start, stop)
    
    Path, time, dist = RODFF_time(start, stop, t0, vship, Roadmap)
    #print(Path)
    t = time[-1] - time[0]

    return t