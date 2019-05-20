import halem.Mesh_maker as Mesh_maker
import halem.Functions as Functions
import halem.Calc_path as Calc_path
import halem.Flow_class as Flow_class

import pytest
import numpy as np
import geopy.distance

def coord_a():
    return (0,0)

def coord_b():
    return (0,1)

def vship():
    return 5

def v(mag):
    v = np.array([[0,0,1,-1]])
    v = mag * np.concatenate((v,v), axis = 0)
    return v

def u(mag):
    u = np.array([[1,-1,0,0]])
    u = mag * np.concatenate((u,u), axis = 0)
    return u


def test_haversine():
    dist = Functions.haversine(coord_a(), coord_a())
    dist1 = Functions.haversine(coord_a(), coord_b())
    dist2 = geopy.distance.geodesic(coord_a(), coord_b()).m

    assert dist == 0 
    assert abs(dist1 - dist2) < 0.01*dist1
    
def test_costfunction_space():
    edge = (0,1)
    nodes = [coord_a(), coord_b()]
    mask = np.full((u(1).shape), False)
    L = Functions.costfunction_spaceseries(edge, vship(), nodes, u(1), v(1), mask)
    dist = Functions.haversine(coord_a(), coord_b()) * np.ones(u(1).shape[1])
    np.testing.assert_array_equal(L,dist)
    
def test_costfunction_time():
    mag = 3
    edge = (0,1)
    nodes = [coord_a(), coord_b()]
    mask = np.full((u(mag).shape), False)
    L = Functions.costfunction_timeseries(edge, vship(), nodes, u(mag), v(mag), mask)
    VV = np.array([vship() + mag, vship() - mag, (vship()**2 - mag**2)**0.5, (vship()**2 - mag**2)**0.5 ])
    dist1 = Functions.haversine(coord_a(), coord_b())
    dist = dist1/VV
    np.testing.assert_array_equal(L,dist)