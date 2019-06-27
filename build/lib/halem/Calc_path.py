from collections import defaultdict
import math
import numpy as np
from numpy import ma
import datetime, time
from datetime import datetime
import halem.Functions as Functions
import halem.Mesh_maker as Mesh_maker

def find_startstop(start, nodes):
    node_x = start[1]
    node_y = start[0]
    
    nodes_x = nodes[:, 1]
    nodes_y = nodes[:, 0]
    
    nx = ((nodes_x - node_x)**2 + (nodes_y - node_y)**2)**0.5
    pt = np.argwhere(nx == nx.min())[0][0]

    return pt

def find_k_time(t, ts):
    QQ = abs(ts - t)
    k = np.argwhere(QQ == QQ.min())[0][0]
    return k

def find_k_repeat(t, ts):
    if t == np.inf:
        return len(ts) -1
    else:
        ts = ts - ts[0]
        t = t - ts[0]
        N = int(t/ts[-1])
        t = t - N * ts[-1]
        
        QQ = abs(ts - t)
        k = np.argwhere(QQ == QQ.min())[0][0]
        return k

def dijsktra(graph, initial, end, t0, graph_functions):           # Typefout

    shortest_paths = {initial: (None, 0)}
    time_paths = {initial: (None, t0)}
    current_node = initial
    visited = set()
    Graph_data = graph
    graph = graph.graph

    find_k = find_k_time if Graph_data.repeat == False else find_k_repeat

    while current_node != end:
        visited.add(current_node)
        destinations = graph.edges[current_node]
        weight_to_current_node = shortest_paths[current_node][1]
        time_to_current_node = time_paths[current_node][1]

        for next_node in destinations:
            k = find_k(time_to_current_node, Graph_data.t)
            weight = weight_to_current_node + graph_functions.weights[(current_node, next_node)][k]
            time = time_to_current_node + graph_functions.time[(current_node, next_node)][k]

            if next_node not in shortest_paths:
                shortest_paths[next_node] = (current_node, weight)
                time_paths[next_node] = (current_node, time)
            else:
                current_shortest_weight = shortest_paths[next_node][1]
                if current_shortest_weight > weight:
                    shortest_paths[next_node] = (current_node, weight)
                    time_paths[next_node] = (current_node, time)
        
        next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
        if not next_destinations:
            return "Route Not Possible"
        # next node is the destination with the lowest weight
        current_node = min(next_destinations, key=lambda k: next_destinations[k][1])
    
    # Work back through destinations in shortest path
    path = []
    while current_node is not None:
        path.append((current_node[0], time_paths[current_node][1], shortest_paths[current_node][1], current_node[1]))
        next_node = shortest_paths[current_node][0]
        current_node = next_node
    # Reverse path
    path = path[::-1]   
    return path

class Has_route:
    def __init__(self,start, stop, graph, t0, graph_functions):
        d = datetime.strptime(t0, "%d/%m/%Y %H:%M:%S")
        t0 = d.timestamp()

        start = find_startstop(start, graph.nodes)
        stop =  find_startstop(stop, graph.nodes)

        self.start = (start, 0)
        self.stop = (stop, 0)

        self.route= np.array(dijsktra(graph, self.start, self.stop, t0, graph_functions))

        self.x_route = np.zeros(len(self.route[:,0]))
        self.y_route = np.zeros(len(self.x_route))
        self.t_route = self.route[:,1]
        for i in range(len(self.x_route)):
            self.x_route[i] = graph.nodes[int(self.route[i,0])][1]
            self.y_route[i] = graph.nodes[int(self.route[i,0])][0]
        
        self.sailing_time = self.t_route[-1]

        for i in range(100):    # Moet verwijderd worden in de clean-up
           self.x_route = np.append(self.x_route, self.x_route[-1])
           self.y_route = np.append(self.y_route, self.y_route[-1])
           self.t_route = np.append(self.t_route, (self.t_route[-1] + 10*60))
            
