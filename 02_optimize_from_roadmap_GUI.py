import math
import numpy as np
from numpy import ma
import netCDF4
from netCDF4 import Dataset, num2date

import TomTom
import TomTom.Mesh_maker as Mesh_maker
import TomTom.Functions as Functions
import TomTom.Calc_path as Calc_path

from matplotlib import pyplot as plt
from matplotlib import gridspec as gridspec
plt.style.use('ggplot')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import pickle
from scipy.interpolate import griddata
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ipywidgets as widgets
from tkinter import *

name_textfile_load = 'D:/DCSM-FM/A06_pieter/Roadmap_Neightbourlayers=1,Vship=0,WD_min=2'
with open(name_textfile_load, 'rb') as input:
    Roadmap = pickle.load(input)

x_r = np.arange(3,7, 0.01)
y_r = np.arange(52,55, 0.01)
y_r, x_r = np.meshgrid(y_r,x_r)
LS_r = griddata((Roadmap.flow.nodes[:,1], Roadmap.flow.nodes[:,0]), Roadmap.flow.WD[73], (x_r, y_r), method= 'linear')


class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = []
        self.ys = []
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.nodes = []

    def __call__(self, event):
        print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        self.nodes.append((event.xdata, event.ydata))

class mclass:
    def __init__(self,  window):
        self.text1 = Text(height=1)
        self.text1.insert(END, 'vship = ')
        self.text1.pack()
        self.vship = Entry()
        self.vship.pack()
        self.text2 = Text(height=1)
        self.text2.insert(END, 't0 = ')
        self.text2.pack()
        self.t0 = Entry()
        self.t0.insert(END, '01/07/2013 04:00:00')
        self.t0.pack()

        self.button1 = Button (window, text="Select Start & End", command=self.plot)
        self.button1.pack()
        self.button2 = Button (window, text="Calculate Optimal Route", command=self.calc)
        self.button2.pack()
    
    def plot(self):
        fig = plt.figure(figsize=(5, 5))
        ax = plt.subplot(projection=ccrs.PlateCarree())
        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.gridlines(color = 'grey', zorder = 3)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='palegoldenrod'))

        line, = plt.plot([], [], 'b')
        points, = plt.plot([],[], 'ro')
        self.linebuilder = LineBuilder(line)
        self.pointbuilder = LineBuilder(points)

        ax.set_extent([4.5, 6, 52.8, 53.8])
        plt.show()
    def calc(self):
        t0 = self.t0.get()
        vship = float(self.vship.get())

        vv= np.abs(Roadmap.vship - vship)
        arg_vship = int(np.argwhere(vv == vv.min()))
        print("V_ship = ", Roadmap.vship[arg_vship])

        class graph_functions_time:
            function_type = "time optimalisation"
            weights = Roadmap.weight_time[arg_vship].weights
            time = Roadmap.weight_time[arg_vship].weights
            
        class graph_functions_space:
            function_type = "length optimalisation"
            weights = Roadmap.weight_space[arg_vship].weights
            time = Roadmap.weight_time[arg_vship].weights
            
        start = (np.array(self.linebuilder.nodes)[0,:])[::-1]
        stop = (np.array(self.linebuilder.nodes)[1,:])[::-1]

        route_time = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
        route_space = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_space)

        fig = plt.figure(figsize=(8, 12))
        gs = gridspec.GridSpec(2, 2)
        ax = plt.subplot(gs[0, 0:2], projection=ccrs.Mercator())
        cval = np.arange(2,30, 0.5)
        plt.contourf(x_r,y_r, LS_r, cval, zorder = 1, transform=ccrs.PlateCarree())
        cbar = plt.colorbar()
        cbar.set_label('Depth [m]')
        plt.plot(route_space.x_route, route_space.y_route, 'y.', transform=ccrs.PlateCarree(), markersize = 5, label = 'Shortest Route')
        plt.plot(route_time.x_route, route_time.y_route, 'm.', transform=ccrs.PlateCarree(), markersize = 5, label = 'Fastest Route')
        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.gridlines(color = 'grey', zorder = 3)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor='palegoldenrod'))
        t = route_time.route[:,2][-1]
        perc = ((route_space.route[:,1][-1] - route_space.route[:,1][0]) - (route_time.route[:,1][-1] - route_time.route[:,1][0]))/ (route_space.route[:,1][-1] - route_space.route[:,1][0])*100
        plt.title("Fastest route is {} hours and {} minutes, \n Shorest route is {} km long,\n Fasest route is {} % faster than the shortest route, \n Ship speed is {} m/s \n"
                    .format(int(t/3600), int((t-int(t/3600)*3600)/60),np.round(route_space.route[:,2][-1]/1000, 2), np.round(perc,2),vship))
        ax.set_extent([4.5, 6, 52.8, 53.8])
        plt.legend(loc = 'best')
        plt.subplot(gs[1, 0])
        TomTom.plot_timeseries(route_space, Roadmap)
        plt.title('s/t Diagram \n Shortest Route')
        plt.subplot(gs[1, 1])
        TomTom.plot_timeseries(route_time, Roadmap)
        plt.title('s/t Diagram \n Fasest Route')
        plt.show()
        


window= Tk()
start= mclass (window)
window.mainloop()