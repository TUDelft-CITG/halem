import datetime

import numpy as np
import pytest

import halem
import halem.path_finder as path_finder


class Roadmap(halem.BaseRoadmap):
    def load(self):
        t = np.arange(0, 100) + 1558077464
        nodes = np.array([(0, 0), (0, 0.001), (0.001, 0.001), (0, 0.003)])

        wd = np.ones((len(t), len(nodes))) * 100
        u = np.ones((len(t), len(nodes))) * 0
        v = np.ones((len(t), len(nodes))) * 0
        u[:, 2] = 10

        return {
            "v": v,
            "u": u,
            "water_depth": wd,
            "time": t,
            "nodes": nodes,
        }


@pytest.fixture(scope="module")
def roadmap():
    blend = 0
    nl = (1, 1)
    dx_min = 0.0000001
    vship = np.array([[3, 4], [4, 5]])
    WD_min = np.array([1, 1])
    WVPI = np.array([5000, 7000])
    number_of_neighbor_layers = 1

    roadmap = Roadmap(
        dx_min=dx_min,
        blend=blend,
        nl=nl,
        number_of_neighbor_layers=number_of_neighbor_layers,
        vship=vship,
        WD_min=WD_min,
        WVPI=WVPI,
    )
    roadmap.parse()
    yield roadmap


@pytest.fixture(scope="module")
def roadmap2():
    blend = 0
    nl = (1, 1)
    dx_min = 0.0000001
    vship = np.array([[3, 4], [4, 5]])
    WD_min = np.array([1, 1])
    WVPI = np.array([5000, 7000])
    number_of_neighbor_layers = 1

    roadmap2 = Roadmap(
        dx_min=dx_min,
        blend=blend,
        nl=nl,
        number_of_neighbor_layers=number_of_neighbor_layers,
        vship=vship,
        WD_min=WD_min,
        WVPI=WVPI,
        repeat=True,
    )
    roadmap2.parse()
    yield roadmap2


class TestCalcPath:
    def test_find_startstop(self):
        nodes = np.array([(0, 0), (0, 1), (1, 0), (1, 1)])
        start = (0.5, 0.5)
        start = path_finder.PathFinder.find_startstop(
            path_finder.PathFinder, start, nodes
        )
        assert start == 0
        start = (0.1, 0)
        start = path_finder.PathFinder.find_startstop(
            path_finder.PathFinder, start, nodes
        )
        assert start == 0
        start = (0.1, 1)
        start = path_finder.PathFinder.find_startstop(
            path_finder.PathFinder, start, nodes
        )
        assert start == 1
        start = (1.1, 1)
        start = path_finder.PathFinder.find_startstop(
            path_finder.PathFinder, start, nodes
        )
        assert start == 3

    def test_find_k(self):
        ts = 100
        t = np.arange(0, 200, 0.33)
        k = path_finder.PathFinder.find_k_time(path_finder.PathFinder, t, ts)
        assert k == 303

    def test_dijstra(self, roadmap):
        vmax = 5
        vvmax = roadmap.vship[:, -1]
        vv = np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])
        t0 = 1558077464
        t0 = datetime.datetime.fromtimestamp(t0)
        t0 = datetime.datetime.strftime(t0, "%d/%m/%Y %H:%M:%S")

        class GraphFunctionsTime:
            function_type = "time optimalisation"
            weights = roadmap.weight_time[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        class GraphFunctionsSpace:
            function_type = "time optimalisation"
            weights = roadmap.weight_space[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        TT = path_finder.PathFinder((0, 0), (0, 3), roadmap, t0, GraphFunctionsTime)
        SS = path_finder.PathFinder((0, 0), (0, 3), roadmap, t0, GraphFunctionsSpace)

        time_path = TT.route
        space_path = SS.route

        assert time_path[1][0] == 2
        assert space_path[1][0] == 1
        assert np.array(time_path)[1, -1] == 1

    def test_PathFinder(self, roadmap):
        start = (0.0001, 0.0001)
        stop = (0.0001, 0.003001)
        t0 = "17/05/2019 9:18:15"
        vmax = 5
        vvmax = roadmap.vship[:, -1]
        vv = np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])

        class GraphFunctionsTime:
            function_type = "time optimalisation"
            weights = roadmap.weight_time[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        class GraphFunctionsSpace:
            function_type = "time optimalisation"
            weights = roadmap.weight_space[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        route_time = path_finder.PathFinder(
            start, stop, roadmap, t0, GraphFunctionsTime
        )
        route_space = path_finder.PathFinder(
            start, stop, roadmap, t0, GraphFunctionsSpace
        )

        assert route_space.route[1, 0] == 1
        assert route_time.route[1, 0] == 2
        assert route_time.y_route[1] == 0.001
        assert route_space.y_route[1] == 0

    def test_find_k_repeat(self):
        ts = 100
        t = np.arange(0, 200, 0.33)
        k = path_finder.PathFinder.find_k_repeat(path_finder.PathFinder, ts, t)
        assert k == 303

        ts = np.inf
        t = np.arange(0, 200, 0.33)
        k = path_finder.PathFinder.find_k_repeat(path_finder.PathFinder, ts, t)
        assert k == len(t) - 1

    def test_HALEM_time(self, roadmap):
        start = (0.0001, 0.0001)
        stop = (0.0001, 0.003001)
        t0 = "17/05/2019 9:18:15"
        vmax = 5
        vvmax = roadmap.vship[:, -1]
        vv = np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])

        class GraphFunctionsTime:
            function_type = "time optimalisation"
            weights = roadmap.weight_time[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        route_time = path_finder.PathFinder(
            start, stop, roadmap, t0, GraphFunctionsTime
        )
        path, _, _ = halem.HALEM_time(start[::-1], stop[::-1], t0, vmax, roadmap)

        np.testing.assert_array_equal(
            roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
        )

    def test_HALEM_space(self, roadmap):
        start = (0.0001, 0.0001)
        stop = (0.0001, 0.003001)
        t0 = "17/05/2019 9:18:15"
        vmax = 5
        vvmax = roadmap.vship[:, -1]
        vv = np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])

        class GraphFunctionsTime:
            function_type = "time optimalisation"
            weights = roadmap.weight_space[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        route_time = path_finder.PathFinder(
            start, stop, roadmap, t0, GraphFunctionsTime
        )
        path, _, _ = halem.HALEM_space(start[::-1], stop[::-1], t0, vmax, roadmap)

        np.testing.assert_array_equal(
            roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
        )

    def test_HALEM_cost(self, roadmap):
        start = (0.0001, 0.0001)
        stop = (0.0001, 0.003001)
        t0 = "17/05/2019 9:18:15"
        vmax = 5
        vvmax = roadmap.vship[:, -1]
        vv = np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])

        class GraphFunctionsTime:
            function_type = "time optimalisation"
            weights = roadmap.weight_cost[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        route_time = path_finder.PathFinder(
            start, stop, roadmap, t0, GraphFunctionsTime
        )
        path, _, _ = halem.HALEM_cost(start[::-1], stop[::-1], t0, vmax, roadmap)

        np.testing.assert_array_equal(
            roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
        )

    def test_HALEM_co2(self, roadmap):
        start = (0.0001, 0.0001)
        stop = (0.0001, 0.003001)
        t0 = "17/05/2019 9:18:15"
        vmax = 5
        vvmax = roadmap.vship[:, -1]
        vv = np.abs(vvmax - vmax)
        arg_vship = int(np.argwhere(vv == vv.min())[0])

        class GraphFunctionsTime:
            function_type = "time optimalisation"
            weights = roadmap.weight_co2[arg_vship].weights
            time = roadmap.weight_time[arg_vship].weights
            vship = roadmap.vship[arg_vship]

        route_time = path_finder.PathFinder(
            start, stop, roadmap, t0, GraphFunctionsTime
        )
        path, _, _ = halem.HALEM_co2(start[::-1], stop[::-1], t0, vmax, roadmap)

        np.testing.assert_array_equal(
            roadmap.nodes[np.array(route_time.route[:, 0], dtype=int)], path[:, ::-1]
        )

    def test_plot(self, roadmap, roadmap2):
        start = (0.0001, 0.0001)
        stop = (0.0001, 0.003001)
        t0 = "17/05/2019 9:18:15"
        vmax = 5

        path, time, _ = halem.HALEM_co2(start[::-1], stop[::-1], t0, vmax, roadmap)
        halem.plot_timeseries(path, time, roadmap, Color="r")

        path, time, _ = halem.HALEM_co2(start[::-1], stop[::-1], t0, vmax, roadmap2)
        halem.plot_timeseries(path, time, roadmap2, Color="r")
