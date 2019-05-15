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

def Length_scale(node, flow, blend, nl):
    nodes = flow.nodes
    nb = find_neighbors(node, flow.tria)
    mag  = (flow.u[:,node]**2 +  flow.v[:,node]**2 )**0.5
    mag = mag.max()

    if len(nb) == 0 :
        return 1
    
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
        graph0 = Graph()
        for from_node in range(len(self.nodes)):       
            to_nodes = find_neighbors2(from_node, self.tria, number_of_neighbor_layers)
            for to_node in to_nodes:
                L = haversine(self.nodes[from_node], self.nodes[int(to_node)])
                graph0.add_edge(from_node, int(to_node), L)
        clear_output(wait= True)

        self.graph = Graph()
        vship = vship[0]
        for edge in graph0.weights:
            for i in range(len(vship)):
                    for j in range(len(vship)):
                        from_node = edge[0]
                        to_node = edge[1]
                        self.graph.add_edge((from_node, i), (to_node, j), 1)
        
        print('3/4')

        'Calculate Weights'
        self.weight_space = []
        self.weight_time = []
        self.weight_cost = []
        self.vship = vship

        for vv in range(len(self.vship[:,-1])):
            graph_time = Graph()
            graph_space = Graph()
            graph_cost = Graph()
            vship = self.vship[vv]
            print(vship)
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