import numpy as np
from halem.particle_tracking import *


def test_particletracking():
    class flow:
        def __init__(self, N=100):
            x = np.linspace(0, 0.1, N)
            y = np.linspace(0, 0.1, N)
            t = np.linspace(0, 10, 1000) * 362
            yy, _, xx = np.meshgrid(y, t, x)
            y, x = np.meshgrid(y, x)

            u = np.cos(np.pi * xx * 100)
            v = -np.cos(np.pi * yy * 100)

            self.x = x
            self.y = y
            self.t = t
            self.u = u
            self.v = v

    f = flow()
    x0 = 0.003
    y0 = 0.003
    T = f.t
    pt = particle_tracking(f)
    RK4 = pt.RK4(x0, y0, T)
    T = np.linspace(0, 10, 1000) * 242
    EF = pt.EF(x0, y0, T)
    assert (EF[0][0] - EF[0][-1]) < EF[0][0] * 0.01
    assert (RK4[0][0] - RK4[0][-1]) < RK4[0][0] * 0.01


def test_time_out_of_range():
    class flow:
        def __init__(self, N=100):
            x = np.linspace(0, 0.1, N)
            y = np.linspace(0, 0.1, N)
            t = np.linspace(0, 10, 1000) * 362
            yy, tt, xx = np.meshgrid(y, t, x)
            y, x = np.meshgrid(y, x)

            u = np.cos(np.pi * xx * 100)
            v = -np.cos(np.pi * yy * 100)

            self.x = x
            self.y = y
            self.t = t
            self.u = u
            self.v = v

    f = flow()
    x0 = 0.003
    y0 = 0.003
    T = f.t + 10000
    pt = particle_tracking(f)
    RK4 = pt.RK4(x0, y0, T)
    EF = pt.EF(x0, y0, T)
    assert len(EF[0]) == 1
    assert len(RK4[0]) == 1
    assert len(EF[1]) == 1
    assert len(RK4[1]) == 1
