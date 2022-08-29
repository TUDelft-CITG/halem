import numpy as np
from IPython.display import clear_output
from numpy import ma

import halem.functions as functions


class NodeReduction:
    """This class can reduce the number of gridpoints of the hydrodynamic model. This
    is done based on the vorticity and the magnitude of the flow. The nodes are pruned
    based on a length scale. The formula for this length scale is:
    LS / ∆min = α(1+|∇×u|)^−βc+(1−α)(1+|u|)^−βm.
    With: LS = resulting length scale, α = blend factor between the curl and
    the magnitude method, ∆min = minimal length scale, βc = non linearity
    parameter for the method with  the curl of the flow, βm = non linearity parameter
    for the method with  the magnitude of the flow, and u = the velocity vector
    of the flow.

    flow:       class that contains the hydrodynamic properties.
                class must have the following instances.
                u: numpy array with shape (N, M)
                v: numpy array with shape (N, M)
                WD: numpy array with shape (N, M)
                nodes: numpy array with shape (N, 2) (lat, lon)
                t: numpy array with shape M (seconds since 01-01-1970 00:00:00)
                tria: triangulation of the nodes (output of scipy.spatial.Delaunay)
                in which N is the number of nodes of the hydrodynamic model, and
                M is the number of time steps of the hydrodynamic model
    dx_min:     float, minimal spatial resolution.
                Parameter of the lengt scale function concerning the node reduction
    blend:      blend factor between the verticity and magnitude of the flow.
                Parameter of the lengt scale function concerning the node reduction
    nl:         float (nl_c, nl_m)
                Non linearity factor consisting out of two numbers
                nl_c non-linearity factor for the corticity, nl_m non-linearity factor
                for the magnitude of the flow. Parameter of the lengt scale function
                concerning the node reduction
    number_of_neighbor_layers:   number of neigbouring layers for which edges are
                created. increasing this number results in a higher directional
                resolution.
    """

    def __init__(self, flow, nl, dx_min, blend):
        self.new_nodes, self.LS = self.get_nodes(flow, nl, dx_min, blend)

    def get_nodes(self, flow, nl, dx_min, blend):
        nodes = flow.nodes
        new_nodes = [0]
        LS = []
        q = int(0)
        qq = 1
        for i in range(len(nodes)):
            q = q + int(1)
            if q == 1000:
                clear_output(wait=True)
                print(np.round(qq / len(nodes) * 100000, 3), "%")
                q = int(0)
                qq += 1
            LS_node = self.length_scale(i, flow, blend, nl)
            LS.append(LS_node)
            closest_nod = self.closest_node(i, new_nodes, nodes)

            y_dist = nodes[closest_nod][0] - nodes[i][0]
            x_dist = nodes[closest_nod][1] - nodes[i][1]
            distu = (y_dist**2 + x_dist**2) ** 0.5

            if distu > dx_min * LS_node:
                new_nodes.append(i)

        LS = ma.array(LS, fill_value=np.nan)

        return new_nodes, LS

    def length_scale(self, node, flow, blend, nl):
        mag = (flow.u[:, node] ** 2 + flow.v[:, node] ** 2) ** 0.5
        mag = mag.max()
        curl = abs(self.curl_func(node, flow))
        LS_c = ma.array(1 / (1 + curl) ** nl[0])
        LS_m = ma.array(1 / (1 + mag) ** nl[1])
        LS = ma.array(blend * LS_c + (1 - blend) * LS_m)
        return LS

    def curl_func(self, node, flow):
        nb = functions.find_neighbors(node, flow.tria)
        nb = np.append(nb, node)
        DUDY = []
        DVDX = []
        xs = flow.nodes[nb][:, 1]
        ys = flow.nodes[nb][:, 0]
        for i in range(len(flow.t)):
            u = flow.u[i, nb]
            v = flow.v[i, nb]
            dudy = float(self.slope(xs, ys, u)[1])
            dvdx = float(self.slope(xs, ys, v)[0])
            DUDY.append(dudy)
            DVDX.append(dvdx)
        DUDY = np.array(DUDY)
        DVDX = np.array(DVDX)
        curl = (np.abs(DUDY - DVDX)).max()

        return curl

    def closest_node(self, node, nodes, node_list):
        """Finds the closest node for a subset of nodes in a set of node.

        based on WGS84 coordinates.

        node:           considered node
        nodes:          indices of the subset
        node_list:      total list of the nodes
        """

        node_x = node_list[node][1]
        node_y = node_list[node][0]

        nodes_x = node_list[nodes][:, 1]
        nodes_y = node_list[nodes][:, 0]

        nx = ((nodes_x - node_x) ** 2 + (nodes_y - node_y) ** 2) ** 0.5
        pt = np.argwhere(nx == nx.min())[0][0]
        pt = nodes[pt]
        return pt

    def slope(self, xs, ys, zs):
        """Function for the slope of a plane in x and y direction.
        Used to calculate the  curl of the flow for the node reduction step"""

        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.matrix(tmp_b).T
        A = np.matrix(tmp_A)
        fit = (A.T * A).I * A.T * b

        return fit[0], fit[1]
