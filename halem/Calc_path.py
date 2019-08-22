import halem.Functions as Functions
from datetime import datetime
import datetime, time
import numpy as np

# Added for openclsim integration
import openclsim.core as core

class Has_route:
    """ This class contains the code for calculating the optimal route from the pre-proccessed Roadmap

    start:              start location (lat, lon)
    stop:               destination location (lat, lon)
    Roadmap:            Preprocessing file
    graph_functions:    class that selects the correct weights from the Roadmap.
    """

    def __init__(self, start, stop, Roadmap, t0, graph_functions):
        d = datetime.datetime.strptime(t0, "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()

        start = self.find_startstop(start, Roadmap.nodes)
        stop = self.find_startstop(stop, Roadmap.nodes)

        self.start = (start, 0)
        self.stop = (stop, 0)

        self.route = np.array(
            self.dijsktra(Roadmap, self.start, self.stop, t0, graph_functions)
        )

        self.x_route = np.zeros(len(self.route[:, 0]))
        self.y_route = np.zeros(len(self.x_route))
        self.t_route = self.route[:, 1]
        for i in range(len(self.x_route)):
            self.x_route[i] = Roadmap.nodes[int(self.route[i, 0])][1]
            self.y_route[i] = Roadmap.nodes[int(self.route[i, 0])][0]

        self.sailing_time = self.t_route[-1]

    def dijsktra(self, Roadmap, initial, end, t0, graph_functions):  # Typefout

        shortest_paths = {initial: (None, 0)}
        time_paths = {initial: (None, t0)}
        current_node = initial
        visited = set()
        find_k = self.find_k_time if Roadmap.repeat == False else self.find_k_repeat

        while current_node != end:
            visited.add(current_node)
            destinations = Roadmap.graph.edges[current_node]
            weight_to_current_node = shortest_paths[current_node][1]
            time_to_current_node = time_paths[current_node][1]

            for next_node in destinations:
                k = find_k(time_to_current_node, Roadmap.t)
                weight = (
                    weight_to_current_node
                    + graph_functions.weights[(current_node, next_node)][k]
                )
                time = (
                    time_to_current_node
                    + graph_functions.time[(current_node, next_node)][k]
                )

                if next_node not in shortest_paths:
                    shortest_paths[next_node] = (current_node, weight)
                    time_paths[next_node] = (current_node, time)
                else:
                    current_shortest_weight = shortest_paths[next_node][1]
                    if current_shortest_weight > weight:
                        shortest_paths[next_node] = (current_node, weight)
                        time_paths[next_node] = (current_node, time)

            next_destinations = {
                node: shortest_paths[node]
                for node in shortest_paths
                if node not in visited
            }
            current_node = min(next_destinations, key=lambda k: next_destinations[k][1])

        path = []
        while current_node is not None:
            path.append(
                (
                    current_node[0],
                    time_paths[current_node][1],
                    shortest_paths[current_node][1],
                    current_node[1],
                )
            )
            next_node = shortest_paths[current_node][0]
            current_node = next_node
        path = path[::-1]
        return path

    def find_startstop(self, start, nodes):
        node_x = start[1]
        node_y = start[0]

        nodes_x = nodes[:, 1]
        nodes_y = nodes[:, 0]

        nx = ((nodes_x - node_x) ** 2 + (nodes_y - node_y) ** 2) ** 0.5
        pt = np.argwhere(nx == nx.min())[0][0]

        return pt

    def find_k_time(self, t, ts):
        QQ = abs(ts - t)
        k = np.argwhere(QQ == QQ.min())[0][0]
        return k

    def find_k_repeat(self, t, ts):
        if t == np.inf:
            return len(ts) - 1
        else:
            ts = ts - ts[0]
            t = t - ts[0]
            N = int(t / ts[-1])
            t = t - N * ts[-1]

            QQ = abs(ts - t)
            k = np.argwhere(QQ == QQ.min())[0][0]
            return k

class Routeable(core.Movable):
    """
    Movement travels trough a Graph. When optimize_route == True the halem package optimizes the route for different cost functions
    if the optimization_type == 'time' then the fastest path is calculated by halem
    if the optimization_type == 'cost' then the cheapest path is calculated by halem
    if the optimization_type == 'space' then the shortest path is calculated by halem
    Optimize_Route == False the route is determind by v = s/t 
    """

    def __init__(
        self,
        route,
        optimize_route=False,
        optimization_type="time",
        loadfactors=None,
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        """Initialization"""
        self.route = route
        self.optimize_route = optimize_route
        self.optimization_type = optimization_type
        self.loadfactors = loadfactors

        if self.optimize_route == True:

            if self.optimization_type == "time":
                self.optimization_func = halem.HALEM_time
            elif self.optimization_type == "space":
                self.optimization_func = halem.HALEM_space
            elif self.optimization_type == "cost":
                self.optimization_func = halem.HALEM_cost
            elif self.optimization_type == "co2":
                self.optimization_func = halem.HALEM_co2
            else:
                print("No known optimization method selected")

class HasDepthRestriction:
    """HasDepthRestriction class

    Used to add depth limits to vessels - improvement of the OpenCLSim class
    draught: should be a lambda function with input variable container.volume
    ukc: list with ukc, corresponding to wave_heights

    filling: filling degree [%]
    min_filling: minimal filling degree [%]
    max_filling: max filling degree [%]
    """

    def check_optimal_filling_Roadmap(self, loader, unloader, origin, destination):
        orig = shapely.geometry.asShape(origin.geometry)
        dest = shapely.geometry.asShape(destination.geometry)
        start = (orig.x, orig.y)
        stop = (dest.x, dest.y)
        prods = []
        tQQ = []

        for i in range(len(self.loadfactors)):
            # Departure time sailing empty
            t0 = datetime.datetime.fromtimestamp(self.env.now).strftime(
                "%d/%m/%Y %H:%M:%S"
            )

            # Duration of sailing empty
            _, TT, _ = self.optimization_func(
                stop, start, t0, self.env.Roadmap.vship[0, -1], self.env.Roadmap
            )

            # Duration of sailing empty + loading
            duration_dredging = self.loading(
                None, None, self.container.capacity * self.loadfactors[i]
            )

            TTT = self.env.now + (TT[-1] - TT[0]) + duration_dredging
            duration_sailing_empty = TT[-1] - TT[0]

            # Departure time sailing full
            t0 = datetime.datetime.fromtimestamp(TTT).strftime("%d/%m/%Y %H:%M:%S")

            # Duration of sailing full
            _, TT, _ = self.optimization_func(
                start, stop, t0, self.env.Roadmap.vship[i, -1], self.env.Roadmap
            )

            # Duration of sailing empty + loading + sailing full + unloading
            duration_unloading = self.unloading(
                None, None, self.container.capacity * self.loadfactors[i]
            )
            TTT += (TT[-1] - TT[0]) + duration_unloading
            duration_sailing_full = TT[-1] - TT[0]

            # Determine production
            prod = (self.loadfactors[i] * self.container.capacity) / (
                duration_sailing_empty
                + duration_sailing_full
                + duration_dredging
                + duration_unloading
            )
            prods.append(prod)
            tQQ.append(t0)

        optimal_loadfactor = self.loadfactors[np.argwhere(prods == max(prods))[0, 0]]
        return optimal_loadfactor