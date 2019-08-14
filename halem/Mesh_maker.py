from collections import defaultdict
import math
import numpy as np
from numpy import ma
import netCDF4
from netCDF4 import Dataset, num2date
from scipy.spatial import Delaunay
import halem.Functions as Functions
from scipy.signal import argrelextrema
from scipy.interpolate import griddata
import datetime, time
from datetime import datetime
import pickle
from IPython.display import clear_output


class Graph_flow_model:
    """Pre-processing function fir the HALEM optimizations. In this fucntion the hydrodynamic
    model and the vessel properties are transformed into weights for the Time dependend Dijkstra
    function.
    
    name_textfile_flow:             string that gives the location of the hydrodynamic model 
                                    in the directory.
    dx_min:                         float, minimal spatial resolution. Parameter of the lengt scale
                                    function concerning the node reduction
    blend:                          blend factor between the verticity and magnitude of the flow.
                                    Parameter of the lengt scale function concerning the node reduction
    nl:                             float (nl_c, nl_m) Non linearity factor consisting out of two numbers
                                    nl_c non-linearity factor for the corticity, nl_m non-linearity factor
                                    for the magnitude of the flow. Parameter of the lengt scale function 
                                    concerning the node reduction
    number_of_neighbor_layers:      number of neigbouring layers for which edges are created. increasing this
                                    number results in a higher directional resolution. 
    vship:                           (N (rows) * M (columns)) numpy array that indicates the sailing 
                                    velocity in deep water. For which N is the number of discretisations
                                    in the load factor, and M is the number of discretisations in the 
                                    dynamic sailing velocity
                                    For the optimization type cost and co2 N must be larger or equal to 2.
    WD_min:                         numpy array with the draft of the vessel. Numpy array has the shape of 
                                    the number of discretisations in the dynamic sailing velocity
    WVPI:                            Numpy array with the total weight of the 
    Load_flow:                      Class that contains the output of the hydrodynamic model. An example is
                                    is provided on https://halem.readthedocs.io/en/latest/examples.html
                                    class must have the following instances. 
                                    u: numpy array with shape (N, M)
                                    v: numpy array with shape (N, M)
                                    WD: numpy array with shape (N, M)
                                    nodes: numpy array with shape (N, 2) (lat, lon)
                                    t: numpy array with shape M (seconds since 01-01-1970 00:00:00)
                                    tria: triangulation of the nodes (output of scipy.spatial.Delaunay(nodes))
                                    in which N is the number of nodes of the hydrodynamic model, and 
                                    M is the number of time steps of the hydrodynamic model
    compute_cost:                   Lambda function that returns the cost for sailing based on the travel
                                    time and the travel velocity.
    compute_co2:                    Lambda function that returns the emmision for sailing based on the travel
                                    time and the travel velocity.      
    WWL:                            Width over Water Line of the vessel in meters 
    LWL:                            Length over Water Line of the vessel in meters   
    ukc:                            Minimal needed under keel clearance in  meters. 
    nodes_on_land:                  Function that adds hydrodynamic conditions on land to if nodes on land are 
                                    not included in the hydrodynamic model
    repeat:                         Indicator if the roadmap can be repeated (True / False) True for 
                                    hydrodynamic models based on a tidal analysis
    optimization_type:              list of optimization types. Excluding one or more not needed optimization 
                                    types can significantly decrease the size of the preprocessing file   
    nodes_index:                    Numpy array that contains the indices of the nodes of the reduced hydrodynamic model.
                                    nodes_index is the output of Roadmap.nodes_index. This option allows you to skip the 
                                    node reduction step if this is already done.                                                           
    """

    def __init__(
        self,
        name_textfile_flow,
        dx_min,
        blend,
        nl,
        number_of_neighbor_layers,
        vship,
        Load_flow,
        WD_min,
        WVPI,
        compute_cost=None,
        compute_co2=None,
        WWL=20,
        LWL=80,
        ukc=1.5,
        nodes_on_land=Functions.nodes_on_land_None,
        repeat=False,
        optimization_type=["time", "space", "cost", "co2"],
        nodes_index=np.array([None]),
    ):
        def compute_cost_f(week_rate, fuel_rate):
            second_rate = week_rate / 7 / 24 / 60 / 60
            return lambda travel_time, speed: (
                travel_time * second_rate + fuel_rate * travel_time * speed ** 3
            )

        def compute_co2_f(fuel_rate):
            return lambda travel_time, speed: (fuel_rate * travel_time * speed ** 3)

        if compute_cost == None:
            compute_cost = compute_cost_f(700_000, 0.0008)
        if compute_co2 == None:
            compute_co2 = compute_co2_f(1)

        self.WWL = WWL
        self.LWL = LWL
        self.ukc = ukc
        self.WVPI = WVPI
        self.repeat = repeat
        self.vship = vship

        # 'Load Flow'
        flow = Load_flow(name_textfile_flow)  # ABC van maken
        print("1/4")

        # 'Calculate nodes and flow conditions in nodes'
        if nodes_index.all() == None:
            reduces_nodes = node_reduction(flow, nl, dx_min, blend)
            self.nodes_index = reduces_nodes.new_nodes
            self.LS = reduces_nodes.LS
        else:
            self.nodes_index = nodes_index

        nodes = flow.nodes[self.nodes_index]
        u = np.asarray(np.transpose(flow.u))[self.nodes_index]
        v = np.asarray(np.transpose(flow.v))[self.nodes_index]
        WD = np.asarray(np.transpose(flow.WD))[self.nodes_index]

        self.nodes, self.u, self.v, self.WD = nodes_on_land(nodes, u, v, WD)

        self.tria = Delaunay(self.nodes)
        self.t = flow.t
        self.mask = np.full(self.u.shape, False)
        self.mask[self.WD < WD_min.max() + ukc] = True
        self.WD_min = WD_min
        clear_output(wait=True)
        print("2/4")

        # 'Calculate edges'
        graph0 = Graph()
        for from_node in range(len(self.nodes)):
            to_nodes = Functions.find_neighbors2(
                from_node, self.tria, number_of_neighbor_layers
            )
            for to_node in to_nodes:
                L = Functions.haversine(self.nodes[from_node], self.nodes[int(to_node)])
                graph0.add_edge(from_node, int(to_node), L)
        clear_output(wait=True)

        self.graph = Graph()
        vship1 = vship[0]
        for edge in graph0.weights:
            for i in range(len(vship1)):
                for j in range(len(vship1)):
                    from_node = edge[0]
                    to_node = edge[1]
                    self.graph.add_edge((from_node, i), (to_node, j), 1)

        print("3/4")

        # 'Calculate Weights'

        if self.repeat == True:
            calc_weights = self.calc_weights_time
        else:
            calc_weights = self.calc_weights_time

        self.weight_space = []  # Moet een Dict worden
        self.weight_time = []
        self.weight_cost = []
        self.weight_co2 = []

        for vv in range(len(self.vship)):
            graph_time = Graph()
            graph_space = Graph()
            graph_cost = Graph()
            graph_co2 = Graph()
            vship = self.vship[vv]
            WD_min = self.WD_min[vv]
            WVPI = self.WVPI[vv]
            for edge in graph0.weights:
                for i in range(len(vship)):
                    for j in range(len(vship)):
                        from_node = edge[0]
                        to_node = edge[1]

                        L, W, euros, co2 = calc_weights(
                            edge,
                            i,
                            j,
                            vship,
                            WD_min,
                            WVPI,
                            self,
                            compute_cost,
                            compute_co2,
                            number_of_neighbor_layers,
                        )

                        graph_time.add_edge((from_node, i), (to_node, j), W)
                        graph_space.add_edge((from_node, i), (to_node, j), L)
                        graph_cost.add_edge((from_node, i), (to_node, j), euros)
                        graph_co2.add_edge((from_node, i), (to_node, j), co2)

            if "space" in optimization_type:
                self.weight_space.append(graph_space)
            if "time" in optimization_type:
                self.weight_time.append(graph_time)
            if "cost" in optimization_type:
                self.weight_cost.append(graph_cost)
            if "co2" in optimization_type:
                self.weight_co2.append(graph_co2)

            clear_output(wait=True)
            print(np.round((vv + 1) / len(self.vship) * 100, 2), "%")

        clear_output(wait=True)
        print("4/4")

    def calc_weights_time(
        self,
        edge,
        i,
        j,
        vship,
        WD_min,
        WVPI,
        self_f,
        compute_cost,
        compute_co2,
        number_of_neighbor_layers,
    ):
        """Function that retruns the weight of an arc"""
        from_node = edge[0]
        W = (
            Functions.costfunction_timeseries(
                edge,
                vship[j],
                WD_min,
                self_f,
                WVPI,
                number_of_neighbor_layers,
                self_f.tria,
            )
            + self_f.t
        )
        W = self.FIFO_maker2(W, self_f.mask[from_node]) - self_f.t

        L = Functions.costfunction_spaceseries(
            edge, vship[j], WD_min, self_f, WVPI, number_of_neighbor_layers, self_f.tria
        )
        L = L + np.arange(len(L)) * (1 / len(L))
        L = self.FIFO_maker2(L, self_f.mask[from_node]) - np.arange(len(L)) * (
            1 / len(L)
        )
        euros = compute_cost(W, vship[j])
        co2 = compute_co2(W, vship[j])

        return L, W, euros, co2

    def FIFO_maker2(self, y, N1):
        """Makes a FIFO time series from a Non-FIFO time series
        y:             Time series
        N1:            Mask file of the time series
        """
        arg = np.squeeze(argrelextrema(y, np.less))
        if arg.shape == ():
            arg = np.array([arg])
        else:
            None
        y_FIFO = 1 * y
        for a in arg:
            loc = np.argwhere(y[: a + 1] <= y[a])[-2:]
            if loc.shape == (2, 1):
                if True in N1[int(loc[0]) : int(loc[1])]:
                    None
                else:
                    y_FIFO[int(loc[0]) : int(loc[1])] = y[a]
        return y_FIFO


class Graph:
    """class that contains the nodes, arcs, and weights for the time-dependent, 
    directional, weighted, and Non-FIFO graph of the route optimization problem.
    This class is used multiple times in the halem.Mesh_maker.Graph_flow_model() 
    function"""

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


class node_reduction:
    def __init__(self, flow, nl, dx_min, blend):
        self.new_nodes, self.LS = self.Get_nodes(flow, nl, dx_min, blend)

    def Get_nodes(self, flow, nl, dx_min, blend):
        nodes = flow.nodes
        new_nodes = [0]
        LS = []
        q = int(0)
        qq = 1
        for i in range(len(nodes)):
            q = q + int(1)
            if q == 1000:
                clear_output(wait=True)
                print(np.round(qq / len(nodes) * 100000, 3), "%")
                q = int(0)
                qq += 1
            LS_node = self.Length_scale(i, flow, blend, nl)
            LS.append(LS_node)
            closest_nod = self.closest_node(i, new_nodes, nodes)

            y_dist = nodes[closest_nod][0] - nodes[i][0]
            x_dist = nodes[closest_nod][1] - nodes[i][1]
            distu = (y_dist ** 2 + x_dist ** 2) ** 0.5

            if distu > dx_min * LS_node:
                new_nodes.append(i)

        LS = ma.array(LS, fill_value=np.nan)

        return new_nodes, LS

    def Length_scale(self, node, flow, blend, nl):
        nb = Functions.find_neighbors(node, flow.tria)
        mag = (flow.u[:, node] ** 2 + flow.v[:, node] ** 2) ** 0.5
        mag = mag.max()
        curl = abs(self.curl_func(node, flow))
        LS_c = ma.array(1 / (1 + curl) ** nl[0])
        LS_m = ma.array(1 / (1 + mag) ** nl[1])
        LS = ma.array(blend * LS_c + (1 - blend) * LS_m)
        return LS

    def curl_func(self, node, flow):
        """"""

        nb = Functions.find_neighbors(node, flow.tria)
        nb = np.append(nb, node)
        DUDY = []
        DVDX = []
        xs = flow.nodes[nb][:, 1]
        ys = flow.nodes[nb][:, 0]
        for i in range(len(flow.t)):
            u = flow.u[i, nb]
            v = flow.v[i, nb]
            dudy = float(self.slope(xs, ys, u)[1])
            dvdx = float(self.slope(xs, ys, v)[0])
            DUDY.append(dudy)
            DVDX.append(dvdx)
        DUDY = np.array(DUDY)
        DVDX = np.array(DVDX)
        curl = (np.abs(DUDY - DVDX)).max()

        return curl

    def closest_node(self, node, nodes, node_list):
        """Finds the closest node for a subset of nodes in a set of node, based on WGS84 coordinates.

        node:           considered node
        nodes:          indices of the subset
        node_list:      total list of the nodes
        """

        node_x = node_list[node][1]
        node_y = node_list[node][0]

        nodes_x = node_list[nodes][:, 1]
        nodes_y = node_list[nodes][:, 0]

        nx = ((nodes_x - node_x) ** 2 + (nodes_y - node_y) ** 2) ** 0.5
        pt = np.argwhere(nx == nx.min())[0][0]
        pt = nodes[pt]
        return pt

    def slope(self, xs, ys, zs):
        """Function for the slope of a plane in x and y direction. 
        Used to calculate the  curl of the flow for the node reduction step"""

        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        fit = (A.T * A).I * A.T * b

        return fit[0], fit[1]
