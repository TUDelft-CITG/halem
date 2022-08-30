import math

import matplotlib.pyplot as plt
import numpy as np

import halem.path_finder as path_finder


def plot_timeseries(path, time, Roadmap, Color="r", range_CP=5):
    """This function can plot the time series for the route
    and shows a contourplot of the unsaiable areas of that route.

    path:       lon, lat coordinates of the route.
                This is in the format of the output
                from halem.HALEM_func[0]
    time:       time series of the path.
                This is in the format of the output
                from halem.HALEM_func[1]
    Roadmap:    Roadmap that is used to calculate the route.
    Color:      Color of the plot of the time series.
                Type sting, with matplotlib color"""

    dist = []
    TT = []
    D = 0
    for i in range(len(path) - 1):
        nx = (
            (Roadmap.nodes[:, 0] - path[i, 1]) ** 2
            + (Roadmap.nodes[:, 1] - path[i, 0]) ** 2
        ) ** 0.5
        idx = np.argwhere(nx == nx.min())[0][0]
        D = D + haversine((path[i, 1], path[i, 0]), (path[i + 1, 1], path[i + 1, 0]))
        dist.append(D)
        T = Roadmap.mask[idx]
        TT.append(T)
    TT = np.array(TT)
    dist = np.array(dist)
    if Roadmap.repeat is True:
        k = path_finder.PathFinder.find_k_repeat(
            path_finder.PathFinder,
            time[0],
            Roadmap.t,
        )
        plt.plot(dist, (time[:-1] - time[0]) / 3600, color=Color, label="s/t route")
        cval = np.arange(0, 1.1, 0.5)

        for J in range(range_CP):
            plt.contourf(
                dist,
                (Roadmap.t - Roadmap.t[k] + (Roadmap.t[-1] - Roadmap.t[0]) * J) / 3600,
                np.transpose(TT),
                cval,
                colors=("cornflowerblue", "sandybrown"),
            )

        plt.colorbar(label="maks file, 0 = False, 1 = True")
        plt.xlabel("traveled distance [m]")
        plt.ylabel("time [h]")
        plt.ylim(0, (time[-1] - time[0]) / 3600 * 1.2)
        plt.legend(loc="best")

    else:
        plt.plot(dist, (time[:-1] - time[0]) / 3600, color=Color, label="s/t route")
        cval = np.arange(0, 1.1, 0.5)
        plt.contourf(
            dist,
            (Roadmap.t - time[0]) / 3600,
            np.transpose(TT),
            cval,
            colors=("cornflowerblue", "sandybrown"),
        )
        plt.colorbar(label="maks file, 0 = False, 1 = True")
        plt.ylim(0, (time[-1] - time[0]) / 3600 * 1.2)
        plt.xlabel("traveled distance [m]")
        plt.ylabel("time [h]")
        plt.legend(loc="best")


def HALEM_func(start, stop, t0, vmax, Roadmap, costfunction):
    """Base of the functions HALEM_time,
    HALEM_cost, HALEM_space,
    HALEM_co2. This function takes the pre-processing
    file, start location, stop location, departure time, and sailing velocity
    and returns the optimized route.

    start:          (lon, lat) coordinates of the start location
    stop:           (lon, lat) coordinates of the destination location
    t0:             string that indcates the departure time
                    ('day'/'month'/'year' 'hour':'minute':'seconds')
    vmax:           (N (rows) * M (columns)) numpy array that indicates the sailing
    velocity in deep water.
                    For which N is the number of discretisations
                    in the load factor, and M is the number of discretisations in the
                    dynamic sailing velocity

                    For the optimization type cost and co2 N must be larger or
                    equal to 2.

    Roadmap:        Preprocessing file that contains the hydrodynamic properties.
    costfunction    Costfunction of the route optimization.
                    Roadmap.weight_time returns fastest route
                    Roadmap.weight_space returns shortest route
                    Roadmap.weight_cost returns cheapest route
                    Roadmap.weight_co2 retruns least pollutant route
    """

    start = start[::-1]
    stop = stop[::-1]

    vvmax = Roadmap.vship[:, -1]
    vv = np.abs(vvmax - vmax)
    arg_vship = int(np.argwhere(vv == vv.min())[0])

    class GraphFunctionsTime:
        weights = costfunction[arg_vship].weights
        time = Roadmap.weight_time[arg_vship].weights
        vship = Roadmap.vship[arg_vship]

    route = path_finder.PathFinder(start, stop, Roadmap, t0, GraphFunctionsTime)
    path = Roadmap.nodes[np.array(route.route[:, 0], dtype=int)]
    time = route.route[:, 1]

    dist = []
    D = 0
    for i in range(route.route[:, 0].shape[0] - 1):
        D = D + haversine(
            (route.y_route[i], route.x_route[i]),
            (route.y_route[i + 1], route.x_route[i + 1]),
        )
        dist.append(D)
    dist = np.array(dist)
    return path[:, ::-1], time, dist


def HALEM_time(start, stop, t0, vmax, Roadmap):
    """Implementation of the function HALEM_func() for the fastest route.

    start:          (lon, lat) coordinates of the start location
    stop:           (lon, lat) coordinates of the destination location
    t0:             string that indcates the departure time
                    ('day'/'month'/'year' 'hour':'minute':'seconds')
    vmax:           (N (rows) * M (columns)) numpy array that indicates the sailing
    velocity in deep water.
                    For which N is the number of discretisations
                    in the load factor, and M is the number of discretisations in the
                    dynamic sailing velocity

                    For the optimization type cost and co2 N must be larger or
                    equal to 2.

    Roadmap:        Preprocessing file that contains the hydrodynamic properties.
    """
    costfunction = Roadmap.weight_time
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def HALEM_space(start, stop, t0, vmax, Roadmap):
    """Implementation of the function HALEM_func() for the shortest route.

    start:          (lon, lat) coordinates of the start location
    stop:           (lon, lat) coordinates of the destination location
    t0:             string that indcates the departure time
                    ('day'/'month'/'year' 'hour':'minute':'seconds')
    vmax:           (N (rows) * M (columns)) numpy array that indicates the sailing
    velocity in deep water.
                    For which N is the number of discretisations
                    in the load factor, and M is the number of discretisations in the
                    dynamic sailing velocity

                    For the optimization type cost and co2 N must be larger or
                    equal to 2.

    Roadmap:        Preprocessing file that contains the hydrodynamic properties."""
    costfunction = Roadmap.weight_space
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def HALEM_cost(start, stop, t0, vmax, Roadmap):
    """Implementation of the function HALEM_func() for the cheapest route.

    start:          (lon, lat) coordinates of the start location
    stop:           (lon, lat) coordinates of the destination location
    t0:             string that indcates the departure time
                    ('day'/'month'/'year' 'hour':'minute':'seconds')
    vmax:           (N (rows) * M (columns)) numpy array that indicates the sailing
    velocity in deep water.
                    For which N is the number of discretisations
                    in the load factor, and M is the number of discretisations in the
                    dynamic sailing velocity

                    For the optimization type cost and co2 N must be larger or
                    equal to 2.

    Roadmap:        Preprocessing file that contains the hydrodynamic properties."""
    costfunction = Roadmap.weight_cost
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def HALEM_co2(start, stop, t0, vmax, Roadmap):
    """Implementation of the function HALEM_func() for the least pollutant route.

    start:          (lon, lat) coordinates of the start location
    stop:           (lon, lat) coordinates of the destination location
    t0:             string that indcates the departure time
                    ('day'/'month'/'year' 'hour':'minute':'seconds')
    vmax:           (N (rows) * M (columns)) numpy array that indicates the sailing
    velocity in deep water.
                    For which N is the number of discretisations
                    in the load factor, and M is the number of discretisations in the
                    dynamic sailing velocity

                    For the optimization type cost and co2 N must be larger or
                    equal to 2.

    Roadmap:        Preprocessing file that contains the hydrodynamic properties."""
    costfunction = Roadmap.weight_co2
    return HALEM_func(start, stop, t0, vmax, Roadmap, costfunction)


def find_neighbors(pindex, triang):
    """Function that can find the neighbours of a Delauney mesh.

    pindex:         Index of the considered node.
    triang:         Triangulation generated with scipy.spatial.Delaunay()
    """
    return triang.vertex_neighbor_vertices[1][
        triang.vertex_neighbor_vertices[0][pindex] : triang.vertex_neighbor_vertices[0][
            pindex + 1
        ]
    ]


def find_neighbors2(index, triang, depth):
    """Function that can find the neighbours of a Delauney mesh, for
    multiple layers of neighbours.

    pindex:         Index of the considered node.
    triang:         Triangulation generated with scipy.spatial.Delaunay()
    Depth:          Number of neigbouring layers (nb)
    """
    buren = np.array([index])
    for _ in range(depth):
        for buur in buren:
            buren_temp = np.array([])
            temp = find_neighbors(int(buur), triang)
            for j in temp:
                if j in buren:
                    None
                else:
                    buren_temp = np.append(buren_temp, int(j))
            buren = np.append(buren, buren_temp)
    buren = np.delete(buren, 0)
    return buren


def squat(h, T, V_max, LWL, WWL, ukc, WVPI):
    """Function for reducing the sailing velocity in deep water to the sailing velocity
    in shallow unconfined waters.

    h:                              Array of the water depth in meters
    V_max:                          Sailing velocity in deep water in meters per second
    WWL:                            Width over Water Line of the vessel in meters
    LWL:                            Length over Water Line of the vessel in meters
    ukc:                            Minimal needed under keel clearance in  meters.
    T:                              numpy array with the draft of the vessel. Numpy
                                    array has the shape of
                                    the number of discretisations in the dynamic
                                    sailing velocity in meters
    WVPI:                           total weight of the the vessel in tf

    V:                              Array of sailing velocities reduced for squat,
                                    corresponding to the input arrat h.
    """
    Gamma_water = 1025
    b = 9 * WWL
    g = 9.81
    ghv2 = g * h / (V_max**2)
    squat_max = h - T - ukc
    CB = WVPI / (LWL * WWL * T * Gamma_water)
    AsAc = (WWL * T) / (b * h - WWL * T)
    RH = (b * h - WWL * T) / (b + 2 * h + 2 * T + WWL)

    VhV1 = (
        1.002
        + 0.005 * (np.sqrt(WWL * T) / RH)
        - 0.1159 * (np.sqrt(WWL * T) / RH) ** 2
        + 0.0191 * (np.sqrt(WWL * T) / RH) ** 3
    )
    V1Vinf = ((np.exp(ghv2) - np.exp(-ghv2)) / (np.exp(ghv2) + np.exp(-ghv2))) ** 0.5

    V_grens = V_max * V1Vinf * VhV1
    V_squat_max = np.zeros(len(h))
    V = np.zeros(len(h))
    for i in range(len(h)):
        V_squat_max[i] = (
            0
            if squat_max[i] < 0
            else (squat_max[i] * 30 / CB / (AsAc[i] ** (2 / 3))) ** (1 / 2.08)
        )
        V[i] = V_grens[i] if V_squat_max[i] > V_grens[i] else V_squat_max[i]

    return V


def inbetweenpoints(start, stop, LL, tria):
    """This node returns the nodes of influence for a specific arc. This function
    retruns the start and stop node plus the nodes in between the start and stop
    node. This function makes sure the route does not jump over hydrodynamic features
    when the neightbouring layers are higher than one.

    start:      (int) index of the start node
    stop:       (int) index of the destination node
    LL:         (int) number of neighbouring layers.
    tria:       triangulation of the nodes (output of scipy.spatial.Delaunay(nodes)
    """

    nodes = [start, stop]
    for L in range(1, LL):
        L = L + 1
        NB = np.array(find_neighbors2(start, tria, L - 1), dtype=int)
        if stop not in NB:
            add_nodes = set(find_neighbors2(stop, tria, L - 1)) & set(
                find_neighbors2(start, tria, L - 1)
            )
            add_nodes = np.array(list(add_nodes), dtype=int)
            nodes = np.concatenate((nodes, add_nodes))
    return nodes


def haversine(coord1, coord2):
    """use the Haversine function to determine the distance between two points
    in the WGS84 coordinate system. Returns the distance between the two points
    in meters.
    Source: https://janakiev.com/blog/gps-points-distance-python/

    coord1:     (lat, lon) coordinates of first point
    coord2:     (lat, lon) coordinates of second point
    """
    R = 6372800
    lat1, lon1 = coord1
    lat2, lon2 = coord2

    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)

    a = (
        math.sin(dphi / 2) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    )

    return 2 * R * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def costfunction_timeseries(edge, V_max, WD_min, flow, WVPI, L, tria):
    """Function that returns the time series of the weights of a specific edge.

    edge:       (int) cosidered edge. edge: index of the location node
                in Roadmap.nodes
    V_max:      Shipping velocity in deep water  in meters per second
    WD_min:     minimal needed draft in meters
    flow:       Class that contains the hydrodynamic conditions
    WVPI:       Weight of the vessel in tf
    L:         (int) number of neighbouring layers.
    tria:       triangulation of the nodes (output of scipy.spatial.Delaunay(nodes)
    """

    xfrom = flow.nodes[edge[0]][1]
    yfrom = flow.nodes[edge[0]][0]
    xto = flow.nodes[edge[1]][1]
    yto = flow.nodes[edge[1]][0]

    IB = inbetweenpoints(edge[0], edge[1], L, tria)

    v_w = flow.v[IB[0]]
    u_w = flow.u[IB[0]]
    WD_W = flow.WD[IB[0]]
    for i in range(1, len(IB)):
        v_w = v_w + flow.v[IB[i]]
        u_w = u_w + flow.u[IB[i]]

        # WD_W = WD_W + flow.WD[IB[i]]
        WD_W = np.minimum(WD_W, flow.WD[IB[i]])

    # WD_W= WD_W / len(IB)
    v_w = v_w / len(IB)
    u_w = u_w / len(IB)
    U_w = (u_w**2 + v_w**2) ** 0.5

    vship = squat(WD_W, WD_min, V_max, flow.LWL, flow.WWL, flow.ukc, WVPI)
    # vship = V_max + 0 * WD_W

    alpha1 = np.arctan2((yto - yfrom), (xto - xfrom))
    alpha2 = np.arctan2(v_w, u_w) - alpha1

    s_t1 = U_w * np.cos(alpha2)
    s_t2 = vship**2 - (U_w * np.sin(alpha2)) ** 2
    s_t = np.array(
        [s_t1[i] + s_t2[i] ** 0.5 if s_t2[i] > 0 else 0 for i in range(len(s_t1))]
    )

    u_t = np.cos(alpha1) * (s_t)
    v_t = np.sin(alpha1) * (s_t)

    L = haversine((yfrom, xfrom), (yto, xto))
    U_t = (u_t**2 + v_t**2) ** 0.5
    t = np.array([L / U_t1 if U_t1 > 0 else np.inf for U_t1 in U_t])

    t[U_t == np.inf] = np.inf
    t[np.isnan(t)] = np.inf
    t[WD_W < WD_min + flow.ukc] = np.inf
    t[WD_W < WD_min + flow.ukc] = np.inf
    t[(U_w * np.sin(alpha2)) ** 2 > vship**2] = np.inf
    t[np.isnan(s_t)] = np.inf
    t[s_t < 0] = np.inf
    return np.array(t)


def costfunction_spaceseries(edge, V_max, WD_min, flow, WVPI, L, tria):
    """Function that returns the time series of the weights of a specifiv edge.

    edge:       (int) cosidered edge. edge: index of the location node
                in Roadmap.nodes
    V_max:      Shipping velocity in deep water  in meters per second
    WD_min:     minimal needed draft in meters
    flow:       Class that contains the hydrodynamic conditions
    WVPI:       Weight of the vessel in tf
    L:         (int) number of neighbouring layers.
    tria:       triangulation of the nodes (output of scipy.spatial.Delaunay(nodes)
    """

    xfrom = flow.nodes[edge[0]][1]
    yfrom = flow.nodes[edge[0]][0]
    xto = flow.nodes[edge[1]][1]
    yto = flow.nodes[edge[1]][0]
    IB = inbetweenpoints(edge[0], edge[1], L, tria)
    v_w = flow.v[IB[0]]
    u_w = flow.u[IB[0]]
    WD_W = flow.WD[IB[0]]
    for i in range(1, len(IB)):
        v_w = v_w + flow.v[IB[i]]
        u_w = u_w + flow.u[IB[i]]

        # WD_W = WD_W + flow.WD[IB[i]]
        WD_W = np.minimum(WD_W, flow.WD[IB[i]])

    # WD_W= WD_W / len(IB)
    v_w = v_w / len(IB)
    u_w = u_w / len(IB)
    U_w = (u_w**2 + v_w**2) ** 0.5

    vship = squat(WD_W, WD_min, V_max, flow.LWL, flow.WWL, flow.ukc, WVPI)

    alpha1 = np.arctan2((yto - yfrom), (xto - xfrom))
    alpha2 = np.arctan2(v_w, u_w) - alpha1

    s_t1 = U_w * np.cos(alpha2)
    s_t2 = vship**2 - (U_w * np.sin(alpha2)) ** 2
    s_t = np.array(
        [s_t1[i] + s_t2[i] ** 0.5 if s_t2[i] > 0 else 0 for i in range(len(s_t1))]
    )

    u_t = np.cos(alpha1) * (s_t)
    v_t = np.sin(alpha1) * (s_t)

    L = haversine((yfrom, xfrom), (yto, xto))
    U_t = (u_t**2 + v_t**2) ** 0.5
    t = np.array([L / U_t1 if U_t1 > 0 else np.inf for U_t1 in U_t])

    t[U_t == np.inf] = np.inf
    t[np.isnan(t)] = np.inf
    t[WD_W < WD_min + flow.ukc] = np.inf
    t[WD_W < WD_min + flow.ukc] = np.inf
    t[(U_w * np.sin(alpha2)) ** 2 > vship**2] = np.inf
    t[np.isnan(s_t)] = np.inf
    t[s_t < 0] = np.inf
    t[t != np.inf] = L

    return np.array(t)
