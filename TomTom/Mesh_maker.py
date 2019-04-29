from collections import defaultdict
import math
import numpy as np
from numpy import ma
import netCDF4
from netCDF4 import Dataset, num2date
from scipy.spatial import Delaunay
import TomTom.Functions as Functions
from scipy.signal import argrelextrema
from scipy.interpolate import griddata
import datetime, time
from datetime import datetime

class Graph():
    def __init__(self):
        """
        self.edges is a dict of all possible next nodes
        e.g. {'X': ['A', 'B', 'C', 'E'], ...}
        self.weights has all the weights between two nodes,
        with the two nodes as a tuple as the key
        e.g. {('X', 'A'): 7, ('X', 'B'): 2, ...}
        """
        self.edges = defaultdict(list)
        self.weights = {}
    
    def add_edge(self, from_node, to_node, weight):
        # Note: assumes edges are directional
        self.edges[from_node].append(to_node)
        self.weights[(from_node, to_node)] = weight

def haversine(coord1, coord2):
    R = 6372800                                      # https://janakiev.com/blog/gps-points-distance-python/
    lat1, lon1 = coord1
    lat2, lon2 = coord2                              # use the Haversine function to determine the distance between two points in the WGS84 coordinate system
    
    phi1, phi2 = math.radians(lat1), math.radians(lat2) 
    dphi       = math.radians(lat2 - lat1)
    dlambda    = math.radians(lon2 - lon1)
    
    a = math.sin(dphi/2)**2 + \
        math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2
    
    return 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))

def find_neighbors(pindex, triang):
    return triang.vertex_neighbor_vertices[1]\
           [triang.vertex_neighbor_vertices[0][pindex]:triang.vertex_neighbor_vertices[0][pindex+1]]

def find_neighbors2(index, triang, depth):
    buren = np.array([index])
    for _ in range(depth):
        for buur in buren:
            buren_temp = np.array([])
            temp = find_neighbors(int(buur), triang)
            for j in temp:
                if j in buren:
                    None
                else:
                    buren_temp = np.append(buren_temp, int(j))
            buren = np.append(buren, buren_temp)
    buren = np.delete(buren, 0)
    return buren

def FIFO_maker(y):
    y[y == np.inf ] = 10000000000
    arg = np.squeeze(argrelextrema(y, np.less))
    y_FIFO = y
    if arg.shape == ():
        loc = np.argwhere(y <= y[arg])[-2:]
        if loc.shape == (2,1):
            y_FIFO[int(loc[0]):int(loc[1])] = y[arg]
        else:
            None
    else:
        for a in arg:
            loc = np.argwhere(y <= y[a])[-2:] 
            if loc.shape == (2,1):
                y_FIFO[int(loc[0]):int(loc[1])] = y[a]
            else:
                None
    return(y_FIFO)

def closest_node(node, nodes, node_list):
    node_x = node_list[node][1]
    node_y = node_list[node][0]
    
    nodes_x = node_list[nodes][:,1]
    nodes_y = node_list[nodes][:,0]
    
    nx = ((nodes_x - node_x)**2 + (nodes_y - node_y)**2)**0.5
    pt = np.argwhere(nx == nx.min())[0][0]
    pt = nodes[pt]
    return pt

def Length_scale(node, flow, blend, nl):
    nodes = flow.nodes
    nb = find_neighbors(node, flow.tria)
    mag  = (flow.u[:,node]**2 +  flow.v[:,node]**2 )**0.5
    mag = mag.max()
    
    dvdx = []
    dudy = []

    for i in range(len(flow.u[:,node])):
        u_v = flow.u[i][nb] - flow.u[i][node]
        v_v = flow.v[i][nb] - flow.v[i][node]
        Delta = nodes[nb] - nodes[node]
        Delta_inv = np.linalg.inv(Delta[:2,:])
        _ , Dudy = Delta_inv.dot(u_v[:2])
        Dvdx, _ = Delta_inv.dot(v_v[:2])
        dudy.append(Dudy)
        dvdx.append(Dvdx)

    curl = abs((np.array(dudy) - np.array(dvdx))/len(nb)).max()
        
    LS_c = ma.array(1 / (1+curl) ** nl[0])
    LS_m = ma.array(1 / (1+mag) ** nl[1])
    LS = ma.array(blend * LS_c + (1-blend)*LS_m)
    return LS

def Get_nodes(flow, nl, dx_min, blend):
    nodes = flow.nodes
    new_nodes =[0]
    LS = []
    for i in range(len(nodes)):
        LS_node = Length_scale(i, flow, blend, nl)  
        LS.append(LS_node)
        closest_nod = closest_node(i, new_nodes, nodes)

        y_dist = nodes[closest_nod][0] - nodes[i][0]
        x_dist = nodes[closest_nod][1] - nodes[i][1]
        distu = (y_dist ** 2 + x_dist ** 2) ** 0.5

        if LS_node.mask == True:
            None
        elif distu > dx_min * LS_node:
            new_nodes.append(i)
        else:
            None

    LS = ma.array(LS, fill_value= np.nan)
    
    return new_nodes, LS

class flow_2D_FM():
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

class Graph_flow_model():
    def __init__(self, name_textfile_flow, dx_min, blend, nl, number_of_neighbor_layers, vship, Load_flow, WD_min):
        self.flow = Load_flow(name_textfile_flow)

        self.nodes_index, self.LS = Get_nodes(self.flow, nl, dx_min, blend)
        self.nodes = self.flow.nodes[self.nodes_index]

        u = np.asarray(np.transpose(self.flow.u))
        self.u = u[self.nodes_index]
        v = np.asarray(np.transpose(self.flow.v))
        self.v = v[self.nodes_index]
        WD = np.asarray(np.transpose(self.flow.WD))
        self.WD = WD[self.nodes_index]

        self.t = self.flow.t
        self.mask = np.full(self.u.shape, False)

        self.mask[self.WD < WD_min] = True

        self.tria = Delaunay(self.nodes)
        self.graph = Graph()


        'Calculate edges'
        for from_node in range(len(self.nodes)):       
            to_nodes = find_neighbors2(from_node, self.tria, number_of_neighbor_layers)
            for to_node in to_nodes:
                L = haversine(self.nodes[from_node], self.nodes[int(to_node)])
                self.graph.add_edge(from_node, int(to_node), L)


        print('1/3')
        'Calculate Weights'
        self.weight_space = []
        self.weight_time = []
        self.vship = vship
        QQ = 0
        for vs in vship:
            graph_time = Graph()
            graph_space = Graph()
            for edge in self.graph.weights:
                from_node = edge[0]
                to_node = edge[1]
                
                W = Functions.costfunction_timeseries(edge, vs, self.nodes, self.u, self.v, self.mask) + self.t
                W = FIFO_maker(W) - self.t
                graph_time.add_edge(from_node, to_node, W)
                
                L = Functions.costfunction_spaceseries(edge, vs, self.nodes, self.u, self.v, self.mask)
                L = L + np.arange(len(L))* (1/len(L))
                L = FIFO_maker(L) - np.arange(len(L))* (1/len(L))
                graph_space.add_edge(from_node, to_node, L)
            
            self.weight_space.append(graph_space)
            self.weight_time.append(graph_time)
            print(QQ, len(vship))
            QQ = QQ +1

        print("3/3")       
