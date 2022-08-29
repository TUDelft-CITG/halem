import halem.Base_functions as halem
import halem.Mesh_maker as Mesh_maker
import halem.Functions as Functions
import halem.Calc_path as Calc_path
import datetime

import pytest
import numpy as np
from scipy.spatial import Delaunay
from IPython.display import clear_output
import os
import pickle


class flow_class:
    def __init__(self, name="maaktnietuit"):
        self.t = np.arange(0, 100) + 1558077464
        self.nodes = np.array([(0, 0), (0, 0.001), (0.001, 0.001), (0, 0.003)])

        self.tria = Delaunay(self.nodes)

        self.WD = np.ones((len(self.t), len(self.nodes))) * 100
        self.u = np.ones((len(self.t), len(self.nodes))) * 0
        self.v = np.ones((len(self.t), len(self.nodes))) * 0
        self.u[:, 2] = 10


name_textfile_flow = "maaktnietuit"
Load_flow = flow_class
blend = 0
nl = (1, 1)
dx_min = 0.0000001
vship = np.array([[3, 4], [4, 5]])
WD_min = np.array([1, 1])
WVPI = np.array([5000, 7000])
ukc = 0


nodes_on_land = Functions.nodes_on_land_None
number_of_neighbor_layers = 1

Roadmap = Mesh_maker.Graph_flow_model(
    name_textfile_flow,
    dx_min,
    blend,
    nl,
    number_of_neighbor_layers,
    vship,
    Load_flow,
    WD_min,
    WVPI,
)

Roadmap2 = Mesh_maker.Graph_flow_model(
    name_textfile_flow,
    dx_min,
    blend,
    nl,
    number_of_neighbor_layers,
    vship,
    Load_flow,
    WD_min,
    WVPI,
    repeat=True,
)
clear_output()


def test_find_startstop():
    nodes = np.array([(0, 0), (0, 1), (1, 0), (1, 1)])
    start = (0.5, 0.5)
    start = Calc_path.Has_route.find_startstop(Calc_path.Has_route, start, nodes)
    assert start == 0
    start = (0.1, 0)
    start = Calc_path.Has_route.find_startstop(Calc_path.Has_route, start, nodes)
    assert start == 0
    start = (0.1, 1)
    start = Calc_path.Has_route.find_startstop(Calc_path.Has_route, start, nodes)
    assert start == 1
    start = (1.1, 1)
    start = Calc_path.Has_route.find_startstop(Calc_path.Has_route, start, nodes)
    assert start == 3


def test_find_k():
    ts = 100
    t = np.arange(0, 200, 0.33)
    k = Calc_path.Has_route.find_k_time(Calc_path.Has_route, t, ts)
    assert k == 303


def test_dijstra():
    vmax = 5
    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])
    t0 = 1558077464
    t0 = datetime.datetime.fromtimestamp(t0)
    t0 = datetime.datetime.strftime(t0, "%d/%m/%Y %H:%M:%S")

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_time[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    class graph_functions_space:
        function_type = "time optimalisation"
        weights = Roadmap.weight_space[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    TT = Calc_path.Has_route((0, 0), (0, 3), Roadmap, t0, graph_functions_time)
    SS = Calc_path.Has_route((0, 0), (0, 3), Roadmap, t0, graph_functions_space)

    time_path = TT.route
    space_path = SS.route

    clear_output()

    assert time_path[1][0] == 2
    assert space_path[1][0] == 1
    assert np.array(time_path)[1, -1] == 1


def test_Has_route():
    start = (0.0001, 0.0001)
    stop = (0.0001, 0.003001)
    t0 = "17/05/2019 9:18:15"
    vmax = 5
    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_time[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    class graph_functions_space:
        function_type = "time optimalisation"
        weights = Roadmap.weight_space[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route_time = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    route_space = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_space)

    assert route_space.route[1, 0] == 1
    assert route_time.route[1, 0] == 2
    assert route_time.y_route[1] == 0.001
    assert route_space.y_route[1] == 0


def test_save_obj():
    halem.save_object(Roadmap, "tests/Data/Roadmap")
    assert os.path.exists("tests/Data/Roadmap")
    with open("tests/Data/Roadmap", "rb") as input:
        Roadmap_load = pickle.load(input)
    os.remove("tests/Data/Roadmap")
    np.testing.assert_array_equal(Roadmap_load.nodes, Roadmap.nodes)
    np.testing.assert_array_equal(Roadmap_load.u, Roadmap.u)
    np.testing.assert_array_equal(Roadmap_load.v, Roadmap.v)
    np.testing.assert_array_equal(Roadmap_load.WD, Roadmap.WD)


def test_find_k_repeat():
    ts = 100
    t = np.arange(0, 200, 0.33)
    k = Calc_path.Has_route.find_k_repeat(Calc_path.Has_route, ts, t)
    assert k == 303

    ts = np.inf
    t = np.arange(0, 200, 0.33)
    k = Calc_path.Has_route.find_k_repeat(Calc_path.Has_route, ts, t)
    assert k == len(t) - 1


def test_HALEM_time():
    start = (0.0001, 0.0001)
    stop = (0.0001, 0.003001)
    t0 = "17/05/2019 9:18:15"
    vmax = 5
    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_time[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route_time = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path, _, _ = halem.HALEM_time(start[::-1], stop[::-1], t0, vmax, Roadmap)

    np.testing.assert_array_equal(
        Roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
    )


def test_HALEM_space():
    start = (0.0001, 0.0001)
    stop = (0.0001, 0.003001)
    t0 = "17/05/2019 9:18:15"
    vmax = 5
    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_space[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route_time = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path, _, _ = halem.HALEM_space(start[::-1], stop[::-1], t0, vmax, Roadmap)

    np.testing.assert_array_equal(
        Roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
    )


def test_HALEM_cost():
    start = (0.0001, 0.0001)
    stop = (0.0001, 0.003001)
    t0 = "17/05/2019 9:18:15"
    vmax = 5
    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_cost[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route_time = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path, _, _ = halem.HALEM_cost(start[::-1], stop[::-1], t0, vmax, Roadmap)

    np.testing.assert_array_equal(
        Roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
    )


def test_HALEM_co2():
    start = (0.0001, 0.0001)
    stop = (0.0001, 0.003001)
    t0 = "17/05/2019 9:18:15"
    vmax = 5
    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class graph_functions_time:
        function_type = "time optimalisation"
        weights = Roadmap.weight_co2[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route_time = Calc_path.Has_route(start, stop, Roadmap, t0, graph_functions_time)
    path, _, _ = halem.HALEM_co2(start[::-1], stop[::-1], t0, vmax, Roadmap)

    np.testing.assert_array_equal(
        Roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
    )


def test_plot():
    start = (0.0001, 0.0001)
    stop = (0.0001, 0.003001)
    t0 = "17/05/2019 9:18:15"
    vmax = 5

    path, time, _ = halem.HALEM_co2(start[::-1], stop[::-1], t0, vmax, Roadmap)
    halem.plot_timeseries2(path, time, Roadmap, Color="r")

    path, time, _ = halem.HALEM_co2(start[::-1], stop[::-1], t0, vmax, Roadmap2)
    halem.plot_timeseries2(path, time, Roadmap2, Color="r")
