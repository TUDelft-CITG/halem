import TomTom.Mesh_maker as Mesh_maker
import TomTom.Functions as Functions
import TomTom.Calc_path as Calc_path

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