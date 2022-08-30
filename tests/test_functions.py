import geopy.distance
import numpy as np
from scipy.spatial import Delaunay

import halem.functions as functions


def coord_a():
    return (0, 0)


def coord_b():
    return (0, 1)


def vship():
    return 5


def v(mag):
    v = np.array([[0, 0, 1, -1]])
    v = mag * np.concatenate((v, v), axis=0)
    return v


def u(mag):
    u = np.array([[1, -1, 0, 0]])
    u = mag * np.concatenate((u, u), axis=0)
    return u


class flow:
    def __init__(self, mag, name="maaktnietuit"):
        self.t = np.arange(0, 100) + 1558077464
        self.nodes = np.array([(0, 0), (0, 1), (1, 1), (0, 3)])

        self.tria = Delaunay(self.nodes)
        self.WD = np.ones((len(self.t), len(self.nodes))) * 1000
        self.u = np.ones((len(self.t), len(self.nodes))) * mag
        self.v = np.ones((len(self.t), len(self.nodes))) * 0
        self.WWL = 1
        self.LWL = 1
        self.ukc = 1


def test_haversine():
    dist = functions.haversine(coord_a(), coord_a())
    dist1 = functions.haversine(coord_a(), coord_b())
    dist2 = geopy.distance.geodesic(coord_a(), coord_b()).m

    assert dist == 0
    assert abs(dist1 - dist2) < 0.01 * dist1


def test_costfunction_time():
    mag = 3
    WD_min = 1
    edge = (0, 1)
    WVPI = 1
    L = functions.costfunction_timeseries(
        edge, vship(), WD_min, flow(3), WVPI, 1, flow(3).tria
    )

    VSHIP = functions.squat(
        flow(3).WD[0], WD_min, vship(), flow(3).LWL, flow(3).WWL, flow(3).ukc, WVPI
    )

    VV = np.array([VSHIP[0] + mag, VSHIP[0] + mag, VSHIP[0] + mag, VSHIP[0] + mag])
    dist1 = functions.haversine(coord_a(), coord_b())
    dist = dist1 / VV

    np.testing.assert_array_equal(L, dist)


def test_costfunction_space():
    WD_min = 1
    edge = (0, 1)
    WVPI = 1
    L = functions.costfunction_spaceseries(
        edge, vship(), WD_min, flow(3), WVPI, 1, flow(3).tria
    )
    dist1 = functions.haversine(coord_a(), coord_b())
    dist = dist1

    np.testing.assert_array_equal(L, dist)


def test_find_neighbor():
    nodes = [
        (3, 3),
        (2, 2),
        (2, 4),
        (4, 2),
        (4, 4),
        (1, 1),
        (1, 3),
        (1, 5),
        (3, 1),
        (3, 5),
        (5, 1),
        (5, 3),
        (5, 5),
        (0, 0),
        (0, 2),
        (0, 4),
        (0, 6),
        (2, 0),
        (4, 0),
        (2, 6),
        (4, 6),
        (6, 0),
        (6, 2),
        (6, 4),
        (6, 6),
    ]

    tria = Delaunay(nodes)
    nb = functions.find_neighbors(0, tria)

    assert len(nb) == 4
    for i in range(1, 5):
        assert i in nb


def test_find_neighbor2():
    nodes = [
        (3, 3),
        (2, 2),
        (2, 4),
        (4, 2),
        (4, 4),
        (1, 1),
        (0.9, 3),
        (1, 5),
        (3, 1),
        (3, 5.1),
        (5, 1),
        (5, 3),
        (5, 5),
        (0, 0),
        (-0.1, 2),
        (-0.1, 4),
        (0, 6),
        (2, 0),
        (4, 0),
        (2, 6.1),
        (4, 6.1),
        (6, 0),
        (6, 2),
        (6, 4.1),
        (6, 6),
    ]
    tria = Delaunay(nodes)

    nb = functions.find_neighbors2(0, tria, 0)
    assert len(nb) == 0

    nb = functions.find_neighbors2(0, tria, 1)
    assert len(nb) == 4
    for i in range(1, 5):
        assert i in nb

    nb = functions.find_neighbors2(0, tria, 2)
    assert len(nb) == 12
    for i in range(1, 13):
        assert i in nb

    nb = functions.find_neighbors2(0, tria, 3)
    assert len(nb) == 24
    for i in range(1, 25):
        assert i in nb


def test_inbetweenpoints():
    x = range(0, 5)
    y = range(0, 5)
    x, y = np.meshgrid(x, y)

    nodes = np.zeros((x.size, 2))
    nodes[:, 1] = x.reshape(x.size)
    nodes[:, 0] = y.reshape(x.size)
    tria = Delaunay(nodes)
    IB = functions.inbetweenpoints(5, 18, 3, tria)

    np.testing.assert_array_equal(IB, np.array([5, 18, 7, 11, 12, 16, 17]))
