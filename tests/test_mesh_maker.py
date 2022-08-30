import numpy as np
from scipy.signal import argrelextrema
from scipy.spatial import Delaunay

from halem.roadmap import BaseRoadmap, Graph, NodeReduction


def test_Graph():
    node1 = 1
    node2 = 2
    node3 = 3
    weight = np.pi
    G = Graph()

    G.add_edge(node1, node2, weight)
    assert G.weights[1, 2] == weight
    assert G.edges[node1] == [node2]
    assert G.edges[node2] == []
    assert G.edges[node3] == []

    G.add_edge(node1, node3, weight)
    assert G.weights[1, 3] == weight
    assert G.edges[node1] == [node2, node3]
    assert G.edges[node2] == []
    assert G.edges[node3] == []

    G.add_edge(node2, node1, weight)
    assert G.weights[2, 1] == weight
    assert G.edges[node1] == [node2, node3]
    assert G.edges[node2] == [node1]
    assert G.edges[node3] == []

    G.add_edge(node2, node3, weight)
    assert G.weights[2, 3] == weight
    assert G.edges[node1] == [node2, node3]
    assert G.edges[node2] == [node1, node3]
    assert G.edges[node3] == []

    G.add_edge(node3, node1, weight)
    assert G.weights[3, 1] == weight
    assert G.edges[node1] == [node2, node3]
    assert G.edges[node2] == [node1, node3]
    assert G.edges[node3] == [node1]

    G.add_edge(node3, node2, weight)
    assert G.weights[3, 2] == weight
    assert G.edges[node1] == [node2, node3]
    assert G.edges[node2] == [node1, node3]
    assert G.edges[node3] == [node1, node2]


class flow_class:
    def load(self):
        t = np.arange(0, 10)

        x = np.arange(0, 10, 0.5)
        y = np.arange(10, 20, 0.5)
        yy, xx = np.meshgrid(y, x)
        xx = xx.reshape(xx.size)
        yy = yy.reshape(yy.size)
        nodes = np.zeros((len(xx), 2))
        nodes[:, 1] = xx
        nodes[:, 0] = yy

        WD = np.ones((len(t), len(nodes))) * 100
        u = [np.sin(np.pi * nodes[:, 1] / 5)]
        v = [np.cos(np.pi * nodes[:, 1] / 5)]
        u = u
        v = v

        for _ in range(len(t) - 1):
            u = np.concatenate((u, u))
            v = np.concatenate((v, v))

        return {
            "v": v,
            "u": u,
            "water_depth": WD,
            "time": t,
            "nodes": nodes,
        }


class flow_class2:
    def load(self):
        t = np.arange(0, 10)

        x = np.arange(0, 10, 0.5)
        y = np.arange(10, 20, 0.5)
        yy, xx = np.meshgrid(y, x)
        xx = xx.reshape(xx.size)
        yy = yy.reshape(yy.size)
        nodes = np.zeros((len(xx), 2))
        nodes[:, 1] = xx
        nodes[:, 0] = yy

        WD = np.ones((len(t), len(nodes))) * 100
        u = [np.sin(np.pi * nodes[:, 0] / 5)]
        v = [np.cos(np.pi * nodes[:, 0] / 5)]
        u = u
        v = v

        for _ in range(len(t) - 1):
            u = np.concatenate((u, u))
            v = np.concatenate((v, v))

        return {
            "v": v,
            "u": u,
            "water_depth": WD,
            "time": t,
            "nodes": nodes,
        }


class flow_class3:
    def load(self):
        t = np.arange(0, 10)

        x = np.arange(0, 10, 0.5)
        y = np.arange(10, 20, 0.5)
        yy, xx = np.meshgrid(y, x)
        xx = xx.reshape(xx.size)
        yy = yy.reshape(yy.size)
        nodes = np.zeros((len(xx), 2))
        nodes[:, 1] = xx
        nodes[:, 0] = yy

        WD = np.ones((len(t), len(nodes))) * 100
        u = [np.sin(np.pi * nodes[:, 0] / 5)]
        v = [-np.cos(np.pi * nodes[:, 1] / 5)]
        u = u
        v = v

        for _ in range(len(t) - 1):
            u = np.concatenate((u, u))
            v = np.concatenate((v, v))

        return {
            "v": v,
            "u": u,
            "water_depth": WD,
            "time": t,
            "nodes": nodes,
        }


class flow_class4:
    def load(self):
        t = np.arange(0, 5)

        x = np.linspace(0, 100, 1100)
        y = np.linspace(0, 100, 1100)
        yy, xx = np.meshgrid(y, x)
        xx = xx.reshape(xx.size)
        yy = yy.reshape(yy.size)
        nodes = np.zeros((len(xx), 2))
        nodes[:, 1] = xx
        nodes[:, 0] = yy

        WD = np.ones((len(t), len(nodes))) * 100
        u = [np.sin(np.pi * nodes[:, 1] / 5)]
        v = [np.cos(np.pi * nodes[:, 1] / 5)]
        u = u
        v = v

        for _ in range(len(t) - 1):
            u = np.concatenate((u, u))
            v = np.concatenate((v, v))

        return {
            "v": v,
            "u": u,
            "water_depth": WD,
            "time": t,
            "nodes": nodes,
        }


class SetFlow:
    def set_flow(self):
        data = self.load()

        self.t = data["time"]
        self.nodes = data["nodes"]
        self.u = data["u"]
        self.v = data["v"]
        self.WD = data["water_depth"]
        self.tria = Delaunay(self.nodes)


def test_fifo_maker():
    x = np.arange(0, 2 * np.pi, 0.01)
    y = 2 * np.sin(x) + x
    N1 = np.full(len(y), False)
    y = BaseRoadmap.fifo_maker(y, N1)
    loc_min = argrelextrema(y, np.less)
    assert len(loc_min[0]) == 0

    x = np.arange(0, 4 * np.pi, 0.01)
    y = 2 * np.sin(x) + x
    y = BaseRoadmap.fifo_maker(y, N1)
    loc_min = argrelextrema(y, np.less)
    assert len(loc_min[0]) == 0

    x = np.arange(0, 2 * np.pi, 0.01)
    y = 2 * np.sin(x) + x
    N1 = np.full(len(y), True)
    y = BaseRoadmap.fifo_maker(y, N1)
    loc_min = argrelextrema(y, np.less)
    assert len(loc_min[0]) == 1


def test_closest_node():
    nodes = np.array([(0, 0), (-1, -1), (-2, 2), (-2, -2), (2, 2), (2, -2), (0, 1)])
    node = 0
    node_list = np.arange(1, 5, dtype=int)

    cn = NodeReduction.closest_node(node, node_list, nodes)
    assert cn == 1


def test_length_scale():
    blend = 0
    nl = (1, 1)
    flow = type("flow", (flow_class, NodeReduction, SetFlow), {})(
        dx_min=0.1, blend=blend, nl=nl
    )
    flow.set_flow()

    for i in range(len(flow.nodes)):
        ls = flow.length_scale(i)
        assert ls == 0.5


def test_length_scale2():
    blend = 1
    nl = (1, 1)
    error = 0
    flow = type("flow", (flow_class, NodeReduction, SetFlow), {})(
        dx_min=0.1, blend=blend, nl=nl
    )
    flow.set_flow()

    for i in range(len(flow.nodes)):
        ls = flow.length_scale(i)
        C = np.pi / 5 * np.sin(2 * np.pi * flow.nodes[i, 1] / 10)
        LS = 1 / (1 + abs(C))
        assert abs(LS - ls) < 0.2 * LS
        e = abs(LS - ls) / LS
        if e > error:
            error = e


def test_length_scale3():
    blend = 1
    nl = (1, 1)
    error = 0
    flow = type("flow", (flow_class2, NodeReduction, SetFlow), {})(
        dx_min=0.1, blend=blend, nl=nl
    )
    flow.set_flow()

    for i in range(len(flow.nodes)):
        ls = flow.length_scale(i)
        C = np.pi / 5 * np.cos(2 * np.pi * flow.nodes[i, 0] / 10)
        LS = 1 / (1 + abs(C))
        assert abs(LS - ls) < 0.2 * LS
        e = abs(LS - ls) / LS
        if e > error:
            error = e


def test_length_scale4():
    blend = 1
    nl = (1, 1)
    error = 0
    flow = type("flow", (flow_class3, NodeReduction, SetFlow), {})(
        dx_min=0.1, blend=blend, nl=nl
    )
    flow.set_flow()
    for i in range(len(flow.nodes)):
        ls = flow.length_scale(i)
        C = (
            np.pi
            / 5
            * (
                np.cos(2 * np.pi * flow.nodes[i, 0] / 10)
                - np.sin(2 * np.pi * flow.nodes[i, 1] / 10)
            )
        )
        LS = 1 / (1 + abs(C))
        assert abs(LS - ls) < 0.2
        e = abs(LS - ls) / LS
        if e > error:
            error = e


def test_node_reduction():
    blend = 0
    nl = (1, 1)
    dx_min = 0.1

    flow = type("flow", (flow_class, NodeReduction, SetFlow), {})(
        dx_min=dx_min, blend=blend, nl=nl
    )
    flow.set_flow()
    new_nodes, LS = flow.get_nodes()

    assert len(new_nodes) == 400


def test_node_reduction2():
    blend = 0
    nl = (1, 1)
    dx_min = 1

    flow = type("flow", (flow_class, NodeReduction, SetFlow), {})(
        dx_min=dx_min, blend=blend, nl=nl
    )
    flow.set_flow()
    new_nodes, LS = flow.get_nodes()

    assert len(new_nodes) == 200
    assert LS.shape == (400,)
