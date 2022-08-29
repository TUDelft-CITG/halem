# import numpy as np
# from IPython.display import clear_output
# from scipy.signal import argrelextrema
# from scipy.spatial import Delaunay

# import halem.mesh_maker as mesh_maker
# from halem.mesh_maker import NodeReduction


# class flow_class:
#     def __init__(self, name="maaktnietuit"):
#         self.t = np.arange(0, 10)

#         x = np.arange(0, 10, 0.5)
#         y = np.arange(10, 20, 0.5)
#         yy, xx = np.meshgrid(y, x)
#         xx = xx.reshape(xx.size)
#         yy = yy.reshape(yy.size)
#         self.nodes = np.zeros((len(xx), 2))
#         self.nodes[:, 1] = xx
#         self.nodes[:, 0] = yy
#         self.tria = Delaunay(self.nodes)

#         self.WD = np.ones((len(self.t), len(self.nodes))) * 100
#         self.u = [np.sin(np.pi * self.nodes[:, 1] / 5)]
#         self.v = [np.cos(np.pi * self.nodes[:, 1] / 5)]
#         u = self.u
#         v = self.v

#         for _ in range(len(self.t) - 1):
#             self.u = np.concatenate((self.u, u))
#             self.v = np.concatenate((self.v, v))


# class flow_class2:
#     def __init__(self):
#         self.t = np.arange(0, 10)

#         x = np.arange(0, 10, 0.5)
#         y = np.arange(10, 20, 0.5)
#         yy, xx = np.meshgrid(y, x)
#         xx = xx.reshape(xx.size)
#         yy = yy.reshape(yy.size)
#         self.nodes = np.zeros((len(xx), 2))
#         self.nodes[:, 1] = xx
#         self.nodes[:, 0] = yy
#         self.tria = Delaunay(self.nodes)

#         self.WD = np.ones((len(self.t), len(self.nodes))) * 100
#         self.u = [np.sin(np.pi * self.nodes[:, 0] / 5)]
#         self.v = [np.cos(np.pi * self.nodes[:, 0] / 5)]
#         u = self.u
#         v = self.v

#         for _ in range(len(self.t) - 1):
#             self.u = np.concatenate((self.u, u))
#             self.v = np.concatenate((self.v, v))


# class flow_class3:
#     def __init__(self):
#         self.t = np.arange(0, 10)

#         x = np.arange(0, 10, 0.5)
#         y = np.arange(10, 20, 0.5)
#         yy, xx = np.meshgrid(y, x)
#         xx = xx.reshape(xx.size)
#         yy = yy.reshape(yy.size)
#         self.nodes = np.zeros((len(xx), 2))
#         self.nodes[:, 1] = xx
#         self.nodes[:, 0] = yy
#         self.tria = Delaunay(self.nodes)

#         self.WD = np.ones((len(self.t), len(self.nodes))) * 100
#         self.u = [np.sin(np.pi * self.nodes[:, 0] / 5)]
#         self.v = [-np.cos(np.pi * self.nodes[:, 1] / 5)]
#         u = self.u
#         v = self.v

#         for _ in range(len(self.t) - 1):
#             self.u = np.concatenate((self.u, u))
#             self.v = np.concatenate((self.v, v))


# class flow_class4:
#     def __init__(self, name="maaktnietuit"):
#         self.t = np.arange(0, 5)

#         x = np.linspace(0, 100, 1100)
#         y = np.linspace(0, 100, 1100)
#         yy, xx = np.meshgrid(y, x)
#         xx = xx.reshape(xx.size)
#         yy = yy.reshape(yy.size)
#         self.nodes = np.zeros((len(xx), 2))
#         self.nodes[:, 1] = xx
#         self.nodes[:, 0] = yy
#         self.tria = Delaunay(self.nodes)

#         self.WD = np.ones((len(self.t), len(self.nodes))) * 100
#         self.u = [np.sin(np.pi * self.nodes[:, 1] / 5)]
#         self.v = [np.cos(np.pi * self.nodes[:, 1] / 5)]
#         u = self.u
#         v = self.v

#         for _ in range(len(self.t) - 1):
#             self.u = np.concatenate((self.u, u))
#             self.v = np.concatenate((self.v, v))


# def test_Graph():
#     node1 = 1
#     node2 = 2
#     node3 = 3
#     weight = np.pi
#     G = mesh_maker.Graph()

#     G.add_edge(node1, node2, weight)
#     assert G.weights[1, 2] == weight
#     assert G.edges[node1] == [node2]
#     assert G.edges[node2] == []
#     assert G.edges[node3] == []

#     G.add_edge(node1, node3, weight)
#     assert G.weights[1, 3] == weight
#     assert G.edges[node1] == [node2, node3]
#     assert G.edges[node2] == []
#     assert G.edges[node3] == []

#     G.add_edge(node2, node1, weight)
#     assert G.weights[2, 1] == weight
#     assert G.edges[node1] == [node2, node3]
#     assert G.edges[node2] == [node1]
#     assert G.edges[node3] == []

#     G.add_edge(node2, node3, weight)
#     assert G.weights[2, 3] == weight
#     assert G.edges[node1] == [node2, node3]
#     assert G.edges[node2] == [node1, node3]
#     assert G.edges[node3] == []

#     G.add_edge(node3, node1, weight)
#     assert G.weights[3, 1] == weight
#     assert G.edges[node1] == [node2, node3]
#     assert G.edges[node2] == [node1, node3]
#     assert G.edges[node3] == [node1]

#     G.add_edge(node3, node2, weight)
#     assert G.weights[3, 2] == weight
#     assert G.edges[node1] == [node2, node3]
#     assert G.edges[node2] == [node1, node3]
#     assert G.edges[node3] == [node1, node2]


# def test_fifo_maker():
#     x = np.arange(0, 2 * np.pi, 0.01)
#     y = 2 * np.sin(x) + x
#     N1 = np.full(len(y), False)
#     y = mesh_maker.GraphFlowModel.fifo_maker(NodeReduction, y, N1)
#     loc_min = argrelextrema(y, np.less)
#     assert len(loc_min[0]) == 0

#     x = np.arange(0, 4 * np.pi, 0.01)
#     y = 2 * np.sin(x) + x
#     y = mesh_maker.GraphFlowModel.fifo_maker(NodeReduction, y, N1)
#     loc_min = argrelextrema(y, np.less)
#     assert len(loc_min[0]) == 0

#     x = np.arange(0, 2 * np.pi, 0.01)
#     y = 2 * np.sin(x) + x
#     N1 = np.full(len(y), True)
#     y = mesh_maker.GraphFlowModel.fifo_maker(NodeReduction, y, N1)
#     loc_min = argrelextrema(y, np.less)
#     assert len(loc_min[0]) == 1


# def test_closest_node():
#     nodes = np.array([(0, 0), (-1, -1), (-2, 2), (-2, -2), (2, 2), (2, -2), (0, 1)])
#     node = 0
#     node_list = np.arange(1, 5, dtype=int)

#     cn = NodeReduction.closest_node(NodeReduction, node, node_list, nodes)
#     assert cn == 1


# def test_length_scale():
#     flow = flow_class()
#     blend = 0
#     nl = (1, 1)
#     NR = NodeReduction(flow, nl, 0.1, blend)

#     for i in range(len(flow.nodes)):
#         ls = NR.length_scale(i, flow, blend, nl)
#         assert ls == 0.5

#     blend = 1
#     nl = (1, 1)
#     error = 0
#     NR = NodeReduction(flow, nl, 0.1, blend)

#     for i in range(len(flow.nodes)):
#         ls = NR.length_scale(i, flow, blend, nl)
#         C = np.pi / 5 * np.sin(2 * np.pi * flow.nodes[i, 1] / 10)
#         LS = 1 / (1 + abs(C))
#         assert abs(LS - ls) < 0.2 * LS
#         e = abs(LS - ls) / LS
#         if e > error:
#             error = e
#     # print(error)

#     flow = flow_class2()
#     blend = 1
#     nl = (1, 1)
#     error = 0
#     NR = NodeReduction(flow, nl, 0.1, blend)

#     for i in range(len(flow.nodes)):
#         ls = NR.length_scale(i, flow, blend, nl)
#         C = np.pi / 5 * np.cos(2 * np.pi * flow.nodes[i, 0] / 10)
#         LS = 1 / (1 + abs(C))
#         assert abs(LS - ls) < 0.2 * LS
#         e = abs(LS - ls) / LS
#         if e > error:
#             error = e
#     # print(error)

#     flow = flow_class3()
#     blend = 1
#     nl = (1, 1)
#     error = 0
#     NR = NodeReduction(flow, nl, 0.1, blend)

#     for i in range(len(flow.nodes)):
#         ls = NR.length_scale(i, flow, blend, nl)
#         C = (
#             np.pi
#             / 5
#             * (
#                 np.cos(2 * np.pi * flow.nodes[i, 0] / 10)
#                 - np.sin(2 * np.pi * flow.nodes[i, 1] / 10)
#             )
#         )
#         LS = 1 / (1 + abs(C))
#         assert abs(LS - ls) < 0.2
#         e = abs(LS - ls) / LS
#         if e > error:
#             error = e


# def test_node_reduction():
#     flow = flow_class()
#     blend = 0
#     nl = (1, 1)
#     dx_min = 0.1

#     reduced_nodes = NodeReduction(flow, nl, dx_min, blend)

#     assert len(reduced_nodes.new_nodes) == 400

#     flow = flow_class()
#     blend = 0
#     nl = (1, 1)
#     dx_min = 1

#     reduced_nodes = NodeReduction(flow, nl, dx_min, blend)

#     assert len(reduced_nodes.new_nodes) == 200
#     assert reduced_nodes.LS.shape == (400,)


# def test_GraphFlowModel():
#     name_textfile_flow = "maaktnietuit"
#     load_flow = flow_class
#     blend = 0
#     nl = (1, 1)
#     dx_min = 0.5
#     vship = np.array([[4], [5]])
#     WD_min = np.array([1, 1])
#     WVPI = np.array([5000, 6000])
#     number_of_neighbor_layers = 1

#     Roadmap = mesh_maker.GraphFlowModel(
#         name_textfile_flow,
#         dx_min,
#         blend,
#         nl,
#         number_of_neighbor_layers,
#         vship,
#         load_flow,
#         WD_min,
#         WVPI,
#     )

#     clear_output()

#     assert Roadmap.v.shape == (400, 10)
#     assert Roadmap.t.shape[0] == 10


# def test_GraphFlowModel_with_indices():
#     nodes_index = np.loadtxt("tests/Data/idx.csv", dtype=int)
#     name_textfile_flow = "maaktnietuit"
#     load_flow = flow_class
#     blend = 0
#     nl = (1, 1)
#     dx_min = 0.5
#     vship = np.array([[4], [5]])
#     WD_min = np.array([1, 1])
#     WVPI = np.array([5000, 6000])
#     number_of_neighbor_layers = 1

#     Roadmap = mesh_maker.GraphFlowModel(
#         name_textfile_flow,
#         dx_min,
#         blend,
#         nl,
#         number_of_neighbor_layers,
#         vship,
#         load_flow,
#         WD_min,
#         WVPI,
#         nodes_index=nodes_index,
#     )

#     clear_output()

#     assert Roadmap.v.shape == (400, 10)
#     assert Roadmap.t.shape[0] == 10


# def test_GraphFlowModel_repeat():
#     name_textfile_flow = "maaktnietuit"
#     load_flow = flow_class
#     blend = 0
#     nl = (1, 1)
#     dx_min = 0.5
#     vship = np.array([[4], [5]])
#     WD_min = np.array([1, 1])
#     WVPI = np.array([5000, 6000])
#     number_of_neighbor_layers = 1

#     Roadmap = mesh_maker.GraphFlowModel(
#         name_textfile_flow,
#         dx_min,
#         blend,
#         nl,
#         number_of_neighbor_layers,
#         vship,
#         load_flow,
#         WD_min,
#         WVPI,
#         repeat=True,
#     )

#     clear_output()

#     assert Roadmap.v.shape == (400, 10)
#     assert Roadmap.t.shape[0] == 10


# def test_percentageprinter():
#     class flow_class:
#         def __init__(self):
#             x = np.arange(0, 1100)
#             y = 0 * np.arange(0, 1100)
#             y[::2] = 5

#             nodes = np.zeros((x.size, 2))
#             nodes[:, 1] = x.reshape(x.size)
#             nodes[:, 0] = y.reshape(x.size)
#             tria = Delaunay(nodes)
#             self.t = np.arange(3)
#             self.nodes = nodes
#             blank = np.zeros((len(self.t), len(nodes)))
#             self.tria = tria
#             self.u = blank
#             self.v = blank
#             self.WD = blank

#     f = flow_class()
#     NodeReduction(f, (0, 0), 1, 0)
