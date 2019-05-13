from collections import defaultdict
import math
import numpy as np
from numpy import ma
import netCDF4
from netCDF4 import Dataset, num2date
from scipy.spatial import Delaunay
import RODFF.Functions as Functions
from scipy.signal import argrelextrema
from scipy.interpolate import griddata
import datetime, time
from datetime import datetime
import pickle
from IPython.display import clear_output

class flow_2D_FM_05nm():
    def __init__(self, name):
        nc = Dataset(name)

        t = nc.variables['time'][:]
        t0 = "22/12/2012 00:00:00"
        d = datetime.strptime(t0, "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()
        self.t = t+t0
        
        x = nc.variables['mesh2d_face_x'][:]
        y = nc.variables['mesh2d_face_y'][:]
        self.WD = nc.variables['mesh2d_waterdepth'][:,:]
        self.u = nc.variables['mesh2d_ucx'][:,:]
        self.v = nc.variables['mesh2d_ucy'][:,:]

        self.nodes = np.zeros((len(x),2))
        self.nodes[:,0] = y
        self.nodes[:,1] = x
        self.tria = Delaunay(self.nodes)

class flow_2D_FM2():
    def __init__(self, name):
        nc = Dataset(name)

        t = nc.variables['time'][:]
        t0 = "22/12/2012 00:00:00"
        d = datetime.strptime(t0, "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()
        self.t = t+t0
        
        x = nc.variables['mesh2d_face_x'][:8000]
        y = nc.variables['mesh2d_face_y'][:8000]
        self.WD = nc.variables['mesh2d_waterdepth'][:,:8000]
        self.u = nc.variables['mesh2d_ucx'][:,:8000]
        self.v = nc.variables['mesh2d_ucy'][:,:8000]

        self.nodes = np.zeros((len(x),2))
        self.nodes[:,0] = y
        self.nodes[:,1] = x
        self.tria = Delaunay(self.nodes)

class flow_2D_FM_100m_tot():
    def __init__(self, name):      
        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0008_map.nc'
        b = -35
        
        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        WD1 = nc.variables['mesh2d_waterdepth'][:,:b]
        u1 = nc.variables['mesh2d_ucx'][:,:b]
        v1 = nc.variables['mesh2d_ucy'][:,:b]
        nodes1 = np.zeros((len(x),2))
        nodes1[:,0] = y
        nodes1[:,1] = x

        print('1/7')

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0009_map.nc'
        a = 13
        b = -101
        
        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][a:b]
        y = nc.variables['mesh2d_face_y'][a:b]
        WD2 = nc.variables['mesh2d_waterdepth'][:,a:b]
        u2 = nc.variables['mesh2d_ucx'][:,a:b]
        v2 = nc.variables['mesh2d_ucy'][:,a:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x

        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        WD1 = np.concatenate((WD1, WD2, ), axis = 1)
        u1 = np.concatenate((u1, u2,), axis = 1)
        v1 = np.concatenate((v1, v2, ), axis = 1)

        print('2/7')

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0018_map.nc'
        b = -173
        
        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        WD2 = nc.variables['mesh2d_waterdepth'][:,:b]
        u2 = nc.variables['mesh2d_ucx'][:,:b]
        v2 = nc.variables['mesh2d_ucy'][:,:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x

        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        WD1 = np.concatenate((WD1, WD2, ), axis = 1)
        u1 = np.concatenate((u1, u2,), axis = 1)
        v1 = np.concatenate((v1, v2, ), axis = 1)

        print('3/7')

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0013_map.nc'
        a = -20000
        b = -200
        
        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][a:b]
        y = nc.variables['mesh2d_face_y'][a:b]
        WD2 = nc.variables['mesh2d_waterdepth'][:,a:b]
        u2 = nc.variables['mesh2d_ucx'][:,a:b]
        v2 = nc.variables['mesh2d_ucy'][:,a:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x

        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        WD1 = np.concatenate((WD1, WD2, ), axis = 1)
        u1 = np.concatenate((u1, u2,), axis = 1)
        v1 = np.concatenate((v1, v2, ), axis = 1)

        print('4/7')

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0007_map.nc'
        b = -300
        
        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        WD2 = nc.variables['mesh2d_waterdepth'][:,:b]
        u2 = nc.variables['mesh2d_ucx'][:,:b]
        v2 = nc.variables['mesh2d_ucy'][:,:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x

        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        WD1 = np.concatenate((WD1, WD2, ), axis = 1)
        u1 = np.concatenate((u1, u2,), axis = 1)
        v1 = np.concatenate((v1, v2, ), axis = 1)

        print('5/7')

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0019_map.nc'
        b = -12533
        
        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        WD2 = nc.variables['mesh2d_waterdepth'][:,:b]
        u2 = nc.variables['mesh2d_ucx'][:,:b]
        v2 = nc.variables['mesh2d_ucy'][:,:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x

        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        WD1 = np.concatenate((WD1, WD2, ), axis = 1)
        u1 = np.concatenate((u1, u2,), axis = 1)
        v1 = np.concatenate((v1, v2, ), axis = 1)

        print('6/7')
        
        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0006_map.nc'
        a = -19000 + 149 
        b = -5649

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][a:b]
        y = nc.variables['mesh2d_face_y'][a:b]
        WD2 = nc.variables['mesh2d_waterdepth'][:,a:b]
        u2 = nc.variables['mesh2d_ucx'][:,a:b]
        v2 = nc.variables['mesh2d_ucy'][:,0:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        
        t = nc.variables['time'][:]
        t0 = "22/12/2012 00:00:00"
        d = datetime.strptime(t0, "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()
        self.t = t+t0

        nodes = np.concatenate((nodes1, nodes2, ), axis = 0)
        WD = np.concatenate((WD1, WD2, ), axis = 1)
        u = np.concatenate((u1, u2,), axis = 1)
        v = np.concatenate((v1, v2, ), axis = 1)

        print('7/7')
        
        idx = np.unique(nodes, axis = 0, return_index=True)[1]       
        nodes = nodes[idx]
        WD = WD[:,idx]
        u = u[:,idx]
        v = v[:,idx]

        self.nodes = nodes
        self.WD = WD
        self.u = u 
        self.v = v
        self.tria = Delaunay(self.nodes)

class flow_3D_FM_05nm():
    def __init__(self, name):
        nc = Dataset(name)

        t = nc.variables['time'][:]
        t0 = "22/12/2011 00:00:00"
        d = datetime.strptime(t0, "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()
        self.t = t+t0
        
        x = nc.variables['mesh2d_face_x'][:8000]
        y = nc.variables['mesh2d_face_y'][:8000]

        self.nodes = np.zeros((len(x),2))
        self.nodes[:,0] = y
        self.nodes[:,1] = x

        self.WD = nc.variables['mesh2d_waterdepth'][:,:8000]
        self.u = nc.variables['mesh2d_ucx'][:,:8000, 0]
        self.v = nc.variables['mesh2d_ucy'][:,:8000, 0]

        self.tria = Delaunay(self.nodes)

class flow_NOOS():
    def __init__(self, name):
        nc = Dataset(name)
        x_domain = (250,380)                      # general-waddden sea
        y_domain = (530,760)
        # x_domain = (300,390)                      # Texel-case
        # y_domain = (650,760)

        v = nc.variables['VELV'][:,:,:]
        u = nc.variables['VELU'][:,:,:]
        d = nc.variables['SEP'][:,:,:]
        x = nc.variables['x'][:,:]
        y = nc.variables['y'][:,:]
        t = nc.variables['time'][:]
        t = t *60
        x = x[x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
        y = y[x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
        u = u[:,x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
        v = v[:,x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]
        d = d[:,x_domain[0]:x_domain[1], y_domain[0]:y_domain[1]]

        x_temp = ma.array(x.reshape(x.size))
        y_temp = ma.array(y.reshape(x.size))

        nodes = np.zeros((y_temp[y_temp.mask == False].size,2))
        nodes[:,0] = y_temp[y_temp.mask == False]
        nodes[:,1] = x_temp[y_temp.mask == False]
        print('1/3')

        bat, nodesb = self.bat()
        Db_new = griddata((nodesb[:,1],nodesb[:,0]), bat, (x,y), method='linear')

        WD = d * 0
        for i in range(d.shape[0]):
            WD[i,:,:] = d[i,:,:] - Db_new

        print('2/3')

        u_n = []
        v_n = []
        d_n = []

        for node in nodes:
            xloc = np.argwhere(x == node[1])[0,1]
            yloc = np.argwhere(y == node[0])[0,0]
            u_n.append(u[:,yloc,xloc])
            v_n.append(v[:,yloc,xloc])
            d_n.append(WD[:,yloc,xloc])

        d_n = np.array(d_n)
        d_n[d_n < -600] = 0
        v_n = np.array(v_n)
        v_n[v_n < -600] = 0
        u_n = np.array(u_n)
        u_n[u_n < -600] = 0
        
        self.nodes = nodes
        self.u = np.transpose(u_n)
        self.v = np.transpose(v_n)
        self.WD = np.transpose(d_n)
        self.tria = Delaunay(nodes)
        self.t = t
        
    def bat(self):
        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0008_map.nc'
        b = -35

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        nodes1 = np.zeros((len(x),2))
        nodes1[:,0] = y
        nodes1[:,1] = x
        bat1 = nc.variables['mesh2d_flowelem_bl'][:b]

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0009_map.nc'
        a = 13
        b = -101

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][a:b]
        y = nc.variables['mesh2d_face_y'][a:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        bat2 = nc.variables['mesh2d_flowelem_bl'][a:b]
        bat1 = np.concatenate((bat1, bat2, ), axis = 0)

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0018_map.nc'
        b = -173

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        bat2 = nc.variables['mesh2d_flowelem_bl'][:b]
        bat1 = np.concatenate((bat1, bat2, ), axis = 0)

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0013_map.nc'
        a = -20000
        b = -200

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][a:b]
        y = nc.variables['mesh2d_face_y'][a:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        bat2 = nc.variables['mesh2d_flowelem_bl'][a:b]
        bat1 = np.concatenate((bat1, bat2, ), axis = 0)

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0007_map.nc'
        b = -300

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        bat2 = nc.variables['mesh2d_flowelem_bl'][:b]
        bat1 = np.concatenate((bat1, bat2, ), axis = 0)

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0019_map.nc'
        b = -12533

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][:b]
        y = nc.variables['mesh2d_face_y'][:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        nodes1 = np.concatenate((nodes1, nodes2, ), axis = 0)
        bat2 = nc.variables['mesh2d_flowelem_bl'][:b]
        bat1 = np.concatenate((bat1, bat2, ), axis = 0)

        name = 'D:/DCSM-FM_100m/A06_pieter/DCSM-FM_100m_0006_map.nc'
        a = -19000 + 149 
        b = -5649

        nc = Dataset(name)
        x = nc.variables['mesh2d_face_x'][a:b]
        y = nc.variables['mesh2d_face_y'][a:b]
        nodes2 = np.zeros((len(x),2))
        nodes2[:,0] = y
        nodes2[:,1] = x
        nodes = np.concatenate((nodes1, nodes2, ), axis = 0)
        bat2 = nc.variables['mesh2d_flowelem_bl'][a:b]
        bat1 = np.concatenate((bat1, bat2, ), axis = 0)

        bat1 = np.array(bat1)

        idx = np.unique(nodes, axis = 0, return_index=True)[1]       
        nodes = nodes[idx]
        bat = bat1[idx]

        return bat, nodes