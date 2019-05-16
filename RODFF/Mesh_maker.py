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
    dist = Functions.haversine(coord1, coord2)
    return dist

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

def FIFO_maker3(y,N1):
    arg = np.squeeze(argrelextrema(y, np.less))
    y_FIFO = 1*y
    if arg.shape == ():
        loc = np.argwhere(y <= y[arg])[-2:]
        if loc.shape == (2,1):
            if True in N1[int(loc[0]):int(loc[1])]:
                None
            else:
                y_FIFO[int(loc[0]):int(loc[1])] = y[arg]
        else:
            None
    else:
        for a in arg:
            loc = np.argwhere(y <= y[a])[-2:] 
            if loc.shape == (2,1):
                if True in N1[int(loc[0]):int(loc[1])]:
                    None
                else:
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

def slope(xs,ys,zs):
    tmp_A = []
    tmp_b = []
    for i in range(len(xs)):
        tmp_A.append([xs[i], ys[i], 1])
        tmp_b.append(zs[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    fit = (A.T * A).I * A.T * b
    errors = b - A * fit
    residual = np.linalg.norm(errors)
    
    return fit[0], fit[1]

def curl_func(node, flow):
    nb = find_neighbors(node, flow.tria)
    nb = np.append(nb,node)
    DUDY = []
    DVDX = []
    xs = flow.nodes[nb][:,1]
    ys = flow.nodes[nb][:,0]
    for i in range(len(flow.t)):
        u = flow.u[i,nb]
        v = flow.v[i,nb]
        dudy = float(slope(xs,ys,u)[1])
        dvdx = float(slope(xs,ys,v)[0])
        DUDY.append(dudy)
        DVDX.append(dvdx)
    return np.mean(DUDY) - np.mean(DVDX)

def Length_scale(node, flow, blend, nl):
    nb = find_neighbors(node, flow.tria)
    mag  = (flow.u[:,node]**2 +  flow.v[:,node]**2 )**0.5
    mag = mag.max()

    if len(nb) < 2 :
        return 1

    curl = abs(curl_func(node, flow))   

    LS_c = ma.array(1 / (1+curl) ** nl[0])
    LS_m = ma.array(1 / (1+mag) ** nl[1])
    LS = ma.array(blend * LS_c + (1-blend)*LS_m)
    return LS

def Get_nodes(flow, nl, dx_min, blend):
    nodes = flow.nodes
    new_nodes =[0]
    LS = []
    q = int(0)
    qq = 1
    for i in range(len(nodes)):
        q = q + int(1)
        if q == 1000:
            clear_output(wait= True)
            print(np.round(qq/ len(nodes)*100000,3), '%')
            q = int(0)
            qq +=1
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

class Graph_flow_model():
    def __init__(self, name_textfile_flow, dx_min, blend, nl, number_of_neighbor_layers, vship, Load_flow, WD_min, compute_cost, nodes_on_land):
        'Load Flow'
        flow = Load_flow(name_textfile_flow)
        print('1/4')

        'Calculate nodes and flow conditions in nodes'
        self.nodes_index, self.LS = Get_nodes(flow, nl, dx_min, blend)

        nodes = flow.nodes[self.nodes_index]
        u = np.asarray(np.transpose(flow.u))[self.nodes_index]
        v = np.asarray(np.transpose(flow.v))[self.nodes_index]
        WD = np.asarray(np.transpose(flow.WD))[self.nodes_index]

        self.nodes, self.u, self.v, self.WD = nodes_on_land(nodes,u,v,WD)

        self.tria = Delaunay(self.nodes)
        self.t = flow.t
        self.mask = np.full(self.u.shape, False)
        self.mask[self.WD < WD_min] = True
        clear_output(wait= True)
        print('2/4')

        'Calculate edges'
        self.vship = vship
        graph0 = Graph()
        for from_node in range(len(self.nodes)):       
            to_nodes = find_neighbors2(from_node, self.tria, number_of_neighbor_layers)
            for to_node in to_nodes:
                L = haversine(self.nodes[from_node], self.nodes[int(to_node)])
                graph0.add_edge(from_node, int(to_node), L)
        clear_output(wait= True)

        self.graph = Graph()
        vship1 = vship[0]
        for edge in graph0.weights:
            for i in range(len(vship1)):
                    for j in range(len(vship1)):
                        from_node = edge[0]
                        to_node = edge[1]
                        self.graph.add_edge((from_node, i), (to_node, j), 1)
        
        print('3/4')

        'Calculate Weights'
        self.weight_space = []
        self.weight_time = []
        self.weight_cost = []
        
        for vv in range(len(self.vship)):
            graph_time = Graph()
            graph_space = Graph()
            graph_cost = Graph()
            vship = self.vship[vv]
            for edge in graph0.weights:
                for i in range(len(vship)):
                        for j in range(len(vship)):
                            from_node = edge[0]
                            to_node = edge[1]
                            W = Functions.costfunction_timeseries(edge, vship[j], self.nodes, self.u, self.v, self.mask) + self.t
                            W = FIFO_maker(W) - self.t
                                            
                            L = Functions.costfunction_spaceseries(edge, vship[j], self.nodes, self.u, self.v, self.mask)
                            L = L + np.arange(len(L))* (1/len(L))
                            L = FIFO_maker(L) - np.arange(len(L))* (1/len(L))
                            euros = compute_cost(W,  vship[j] )

                            graph_time.add_edge((from_node, i), (to_node, j), W)
                            graph_space.add_edge((from_node, i), (to_node,j), L)
                            graph_cost.add_edge((from_node, i), (to_node, j), euros)
            
            self.weight_space.append(graph_space)
            self.weight_time.append(graph_time)
            self.weight_cost.append(graph_cost)
            
        clear_output(wait= True)
        print("4/4")             

class Graph_flow_model_with_indices():
    def __init__(self, name_textfile_flow, nodes_index, number_of_neighbor_layers, vship, Load_flow, WD_min, compute_cost, nodes_on_land):
        'Load Flow'
        flow = Load_flow(name_textfile_flow)
        print('1/4')

        'Calculate nodes and flow conditions in nodes'
        self.nodes_index = nodes_index

        nodes = flow.nodes[self.nodes_index]
        u = np.asarray(np.transpose(flow.u))[self.nodes_index]
        v = np.asarray(np.transpose(flow.v))[self.nodes_index]
        WD = np.asarray(np.transpose(flow.WD))[self.nodes_index]

        self.nodes, self.u, self.v, self.WD = nodes_on_land(nodes,u,v,WD)

        self.tria = Delaunay(self.nodes)
        self.t = flow.t
        self.mask = np.full(self.u.shape, False)
        self.mask[self.WD < WD_min] = True
        clear_output(wait= True)
        print('2/4')

        'Calculate edges'
        self.vship = vship
        graph0 = Graph()
        for from_node in range(len(self.nodes)):       
            to_nodes = find_neighbors2(from_node, self.tria, number_of_neighbor_layers)
            for to_node in to_nodes:
                L = haversine(self.nodes[from_node], self.nodes[int(to_node)])
                graph0.add_edge(from_node, int(to_node), L)
        clear_output(wait= True)

        self.graph = Graph()
        vship1 = vship[0]
        for edge in graph0.weights:
            for i in range(len(vship1)):
                    for j in range(len(vship1)):
                        from_node = edge[0]
                        to_node = edge[1]
                        self.graph.add_edge((from_node, i), (to_node, j), 1)
        
        print('3/4')

        'Calculate Weights'
        self.weight_space = []
        self.weight_time = []
        self.weight_cost = []
        
        for vv in range(len(self.vship)):
            graph_time = Graph()
            graph_space = Graph()
            graph_cost = Graph()
            vship = self.vship[vv]
            for edge in graph0.weights:
                for i in range(len(vship)):
                        for j in range(len(vship)):
                            from_node = edge[0]
                            to_node = edge[1]
                            W = Functions.costfunction_timeseries(edge, vship[j], self.nodes, self.u, self.v, self.mask) + self.t
                            W = FIFO_maker(W) - self.t
                                            
                            L = Functions.costfunction_spaceseries(edge, vship[j], self.nodes, self.u, self.v, self.mask)
                            L = L + np.arange(len(L))* (1/len(L))
                            L = FIFO_maker(L) - np.arange(len(L))* (1/len(L))
                            euros = compute_cost(W,  vship[j] )

                            graph_time.add_edge((from_node, i), (to_node, j), W)
                            graph_space.add_edge((from_node, i), (to_node,j), L)
                            graph_cost.add_edge((from_node, i), (to_node, j), euros)
            
            self.weight_space.append(graph_space)
            self.weight_time.append(graph_time)
            self.weight_cost.append(graph_cost)
            
        clear_output(wait= True)
        print("4/4")       