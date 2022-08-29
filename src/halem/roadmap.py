from abc import ABC, abstractmethod
from collections import defaultdict

import numpy as np
import scipy.spatial
from numpy import ma
from scipy.signal import argrelextrema

import halem.functions as functions


class Graph:
    def __init__(self):
        self.edges = defaultdict(list)
        self.weights = {}

    def add_edge(self, from_node, to_node, weight):
        # Note: assumes edges are directional
        self.edges[from_node].append(to_node)
        self.weights[(from_node, to_node)] = weight


class NodeReduction:
    def __init__(
        self,
        dx_min,
        blend,
        nl,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.dx_min = dx_min
        self.blend = blend
        self.nl = nl

    def get_nodes(self):
        nodes = self.nodes
        new_nodes = [0]
        LS = []
        for i in range(len(nodes)):
            LS_node = self.length_scale(i)
            LS.append(LS_node)
            closest_nod = self.closest_node(i, new_nodes, nodes)

            y_dist = nodes[closest_nod][0] - nodes[i][0]
            x_dist = nodes[closest_nod][1] - nodes[i][1]
            distu = (y_dist**2 + x_dist**2) ** 0.5

            if distu > self.dx_min * LS_node:
                new_nodes.append(i)
        return new_nodes, ma.array(LS, fill_value=np.nan)

    def length_scale(self, node):
        mag = (self.u[:, node] ** 2 + self.v[:, node] ** 2) ** 0.5
        mag = mag.max()
        curl = abs(self.curl_func(node))
        LS_c = ma.array(1 / (1 + curl) ** self.nl[0])
        LS_m = ma.array(1 / (1 + mag) ** self.nl[1])
        LS = ma.array(self.blend * LS_c + (1 - self.blend) * LS_m)
        return LS

    def curl_func(self, node):
        nb = functions.find_neighbors(node, self.tria)
        nb = np.append(nb, node)
        DUDY = []
        DVDX = []
        xs = self.nodes[nb][:, 1]
        ys = self.nodes[nb][:, 0]
        for i in range(len(self.t)):
            u = self.u[i, nb]
            v = self.v[i, nb]
            dudy = float(self.slope(xs, ys, u)[1])
            dvdx = float(self.slope(xs, ys, v)[0])
            DUDY.append(dudy)
            DVDX.append(dvdx)
        DUDY = np.array(DUDY)
        DVDX = np.array(DVDX)
        curl = (np.abs(DUDY - DVDX)).max()

        return curl

    def closest_node(self, node, nodes, node_list):
        """Finds the closest node for a subset of nodes in a set of node.

        based on WGS84 coordinates.

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


class BaseRoadmap(ABC, NodeReduction):
    def __init__(
        self,
        number_of_neighbor_layers,
        vship,
        WD_min,
        WVPI,
        repeat=False,
        WWL=20,
        LWL=80,
        ukc=1.5,
        optimization_type=None,
        nodes_index=None,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.optimization_type = optimization_type or ["time", "space", "cost", "co2"]
        self.nodes_index = nodes_index
        self.number_of_neighbor_layers = number_of_neighbor_layers
        self.WD_min = WD_min
        self.WWL = WWL
        self.LWL = LWL
        self.ukc = ukc
        self.WVPI = WVPI
        self.vship = vship
        self.repeat = repeat

    @staticmethod
    def nodes_on_land(nodes, u, v, WD):
        """Standard function that returns itself"""
        return nodes, u, v, WD

    @staticmethod
    def compute_cost(travel_time, speed):
        week_rate = 700_000
        fuel_rate = 0.008
        second_rate = week_rate / 7 / 24 / 60 / 60
        return travel_time * second_rate + fuel_rate * travel_time * speed**3

    @staticmethod
    def compute_co2(travel_time, speed):
        fuel_rate = 1
        return fuel_rate * travel_time * speed**3

    @abstractmethod
    def load():
        return {
            "time": np.array([]),
            "nodes": np.array([]),
            "u": np.array([]),
            "v": np.array([]),
            "water_depth": np.array([]),
        }

    def load_hydrodynamic(self):
        data = self.load()

        self.t = data["time"]
        self.nodes = data["nodes"]
        self.u = data["u"]
        self.v = data["v"]
        self.WD = data["water_depth"]

        assert isinstance(self.t, np.ndarray)
        assert isinstance(self.nodes, np.ndarray)
        assert isinstance(self.u, np.ndarray)
        assert isinstance(self.v, np.ndarray)
        assert isinstance(self.WD, np.ndarray)

        assert len(self.nodes.shape) == 2
        len_nodes = len(self.nodes)

        assert len(self.t.shape) == 1
        len_t = len(self.t)

        assert self.u.shape == (len_t, len_nodes)
        assert self.v.shape == (len_t, len_nodes)
        assert self.WD.shape == (len_t, len_nodes)

        assert np.issubdtype(self.nodes.dtype, np.number)
        assert np.issubdtype(self.t.dtype, np.number)
        assert np.issubdtype(self.u.dtype, np.number)
        assert np.issubdtype(self.v.dtype, np.number)
        assert np.issubdtype(self.WD.dtype, np.number)

        self.tria = scipy.spatial.Delaunay(self.nodes)

    def parse(self):
        self.load_hydrodynamic()

        # 'Calculate nodes and flow conditions in nodes'
        if self.nodes_index is None:
            self.nodes_index, self.LS = self.get_nodes()
        else:
            self.nodes_index = self.nodes_index
            self.LS = None

        nodes = self.nodes[self.nodes_index]
        u = np.asarray(np.transpose(self.u))[self.nodes_index]
        v = np.asarray(np.transpose(self.v))[self.nodes_index]
        WD = np.asarray(np.transpose(self.WD))[self.nodes_index]

        self.nodes, self.u, self.v, self.WD = self.nodes_on_land(nodes, u, v, WD)
        self.tria = scipy.spatial.Delaunay(self.nodes)
        self.mask = np.full(self.u.shape, False)
        self.mask[self.WD < self.WD_min.max() + self.ukc] = True

        # 'Calculate edges'
        graph0 = Graph()
        for from_node in range(len(self.nodes)):
            to_nodes = functions.find_neighbors2(
                from_node, self.tria, self.number_of_neighbor_layers
            )
            for to_node in to_nodes:
                L = functions.haversine(self.nodes[from_node], self.nodes[int(to_node)])
                graph0.add_edge(from_node, int(to_node), L)

        self.graph = Graph()
        vship1 = self.vship[0]
        for edge in graph0.weights:
            for i in range(len(vship1)):
                for j in range(len(vship1)):
                    from_node = edge[0]
                    to_node = edge[1]
                    self.graph.add_edge((from_node, i), (to_node, j), 1)

        # 'Calculate Weights'
        calc_weights = self.calc_weights_time
        self.weight_space = []
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
                            self.compute_cost,
                            self.compute_co2,
                            self.number_of_neighbor_layers,
                        )

                        graph_time.add_edge((from_node, i), (to_node, j), W)
                        graph_space.add_edge((from_node, i), (to_node, j), L)
                        graph_cost.add_edge((from_node, i), (to_node, j), euros)
                        graph_co2.add_edge((from_node, i), (to_node, j), co2)

            if "space" in self.optimization_type:
                self.weight_space.append(graph_space)
            if "time" in self.optimization_type:
                self.weight_time.append(graph_time)
            if "cost" in self.optimization_type:
                self.weight_cost.append(graph_cost)
            if "co2" in self.optimization_type:
                self.weight_co2.append(graph_co2)

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
            functions.costfunction_timeseries(
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
        W = self.fifo_maker(W, self_f.mask[from_node]) - self_f.t

        L = functions.costfunction_spaceseries(
            edge, vship[j], WD_min, self_f, WVPI, number_of_neighbor_layers, self_f.tria
        )
        L = L + np.arange(len(L)) * (1 / len(L))
        L = self.fifo_maker(L, self_f.mask[from_node]) - np.arange(len(L)) * (
            1 / len(L)
        )
        euros = compute_cost(W, vship[j])
        co2 = compute_co2(W, vship[j])

        return L, W, euros, co2

    def fifo_maker(self, y, N1):
        """Makes a FIFO time series from a Non-FIFO time series
        y:  Time series
        N1: Mask file of the time series
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
