from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import Delaunay

import halem


class RoadMap(halem.BaseRoadmap):
    def load(self):
        d = datetime.strptime("23/03/2019 00:00:00", "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()
        x = np.arange(0, 1, 0.04)
        y = np.arange(0, 1, 0.04)
        t = np.arange(t0, (t0 + 2 * 30 * 60 * 30), 30 * 60)

        y, x = np.meshgrid(y, x)
        y = y.reshape(y.size)
        x = x.reshape(y.size)

        nodes = np.zeros((len(x), 2))
        nodes[:, 0] = y
        nodes[:, 1] = x

        u = []
        v = []
        for node in nodes:
            ut = 0 * t + 2 * np.cos(np.pi * (node[0]))
            vt = 0 * t - 2 * np.cos(np.pi * (node[1]))
            u.append(ut)
            v.append(vt)

        v = np.transpose(np.array(v))
        u = np.transpose(np.array(u))
        wd = u * 0 + 20

        return {
            "v": v,
            "u": u,
            "water_depth": wd,
            "time": t,
            "nodes": nodes,
        }


def compute_cost(week_rate, fuel_rate):
    second_rate = week_rate / 7 / 24 / 60 / 60
    return lambda travel_time, speed: (
        travel_time * second_rate + fuel_rate * travel_time * speed**3
    )


def compute_co2(fuel_rate):
    return lambda travel_time, speed: (fuel_rate * travel_time * speed**3)


def main():
    roadmap = RoadMap(
        number_of_neighbor_layers=2,
        vship=np.array([[1, 4, 8], [1, 5, 10]]),
        WD_min=np.array([8, 7]),
        WVPI=np.array([4000, 10000]),
        compute_cost=compute_cost(700_000, 0.008),
        compute_co2=compute_co2(1),
        WWL=40,
        LWL=80,
        ukc=1.5,
        nl=(3, 2.5),
        dx_min=0.1,
        blend=1,
    )

    roadmap.parse()


if __name__ == "__main__":
    main()
