from functools import cache
from math import sqrt, floor, ceil, log2

import numpy as np
import networkx as nx
from sklearn.cluster import KMeans
from networkx.algorithms.approximation.steinertree import metric_closure

DELTA = 1E-2  # controls the precision of distance calculations


def radius(center, n, dist):
    ret = 0
    for v in range(n):
        if v == center: continue
        ret = max(ret, dist(center, v))
    return ret


def trivial_ftp_ub(n, dist):
    minimum_depth = ceil(log2(n))  # minimum depth for a balanced wake-up tree

    diameter = 0
    for u in range(n):
        for v in range(u + 1, n):
            diameter = max(diameter, dist(u, v))
    return diameter * minimum_depth


def min_dist(n, dist):
    ret = float('inf')
    for u in range(n):
        for v in range(u + 1, n):
            ret = min(ret, dist(u, v))
    return ret


class Dist:
    def __init__(self, n):
        self.n = n
        self.func = None

    def __call__(self, u, v):
        return self.func(u, v)

    def __len__(self):
        return self.n

    def __getitem__(self, points):  # restrict to a sub-set of points
        ret = Dist(len(points))

        def aux(u, v):
            u, v = points[u], points[v]
            return self.func(u, v)

        ret.func = aux
        return ret


class L2Norm(Dist):
    def __init__(self, coords, delta=DELTA):
        super().__init__(len(coords))
        self.coords = coords
        self.delta = delta

        @cache
        def aux(u, v):
            d = sqrt((self.coords[u][0] - self.coords[v][0]) ** 2 + (self.coords[u][1] - self.coords[v][1]) ** 2)
            return floor(d / self.delta)  # scale up so we can treat rational values as integers

        self.func = aux


class GraphDist(Dist):
    def __init__(self, edges, delta=DELTA):
        self.graph = metric_closure(nx.Graph(edges))
        super().__init__(self.graph.number_of_nodes())
        self.delta = delta

        @cache
        def aux(u, v):
            if u == v: return 0
            d = self.graph.get_edge_data(u, v)['distance']
            return floor(d / self.delta)  # scale up so we can treat rational values as integers

        self.func = aux


def greedy_solution(source, n, dist):
    # todo: optimize this function
    if n <= 1: return [], 0

    # Connect the source to the closest node:
    min_target = None
    min_arc_weight = float('inf')
    for v in range(n):
        if v == source: continue
        aux = dist(source, v)
        if aux < min_arc_weight:
            min_target = v
            min_arc_weight = aux
    node_activation = {source: 0, min_target: min_arc_weight}
    sol_edges = [(source, min_target)]
    makespan = min_arc_weight

    # Init the degree map (-1 are not yet border nodes and -2 saturated nodes):
    degree = {v: -1 for v in range(n)}
    degree[source] = -2
    degree[min_target] = 0

    border = [min_target]
    while len(sol_edges) < n - 1:
        min_arc = (None, None)
        min_makespan = float('inf')
        for u in border:
            for v in range(n):
                if v == source or v == u: continue
                if degree[v] == -1:  # not yet added
                    aux = node_activation[u] + dist(u, v)
                    if aux < min_makespan:
                        min_arc = (u, v)
                        min_makespan = aux
        u, v = min_arc
        sol_edges.append((u, v))
        node_activation[v] = min_makespan
        makespan = max(makespan, min_makespan)
        degree[u] += 1
        if degree[u] == 3 - (u != source):
            border.remove(u)
        degree[v] = 0
        border.append(v)
    return sol_edges, makespan


def normalize(coords, eps):
    grid_dim = ceil(1.0 / eps)
    xs, ys = zip(*coords)  # gather all points x and y coordinates

    # Stick the bounding box corner to the origin:
    min_x_coord = min(xs)
    min_y_coord = min(ys)
    xs = [-min_x_coord + x for x in xs]
    ys = [-min_y_coord + y for y in ys]
    bounding_box_size = max(*xs, *ys) + 1  # add one to avoid boundary points

    # Resize the points to fit in a ceil(1 / eps) x ceil(1 / eps) square:
    factor = grid_dim / bounding_box_size
    return [(factor * x, factor * y) for x, y in zip(xs, ys)], factor


def discretize(names, coords, eps):  # the source is assumed to be the last node
    grid_dim = ceil(1.0 / eps)
    grid = [[[] for _ in range(grid_dim)] for _ in range(grid_dim)]

    source = 0
    for v, (x, y) in zip(names, coords):
        if v == source:
            grid[floor(x)][floor(y)].insert(0, source)  # make sure the source is a representative
        else:
            grid[floor(x)][floor(y)].append(v)

    grid_map = {}
    rep_names, rep_coords, rep_degrees = [], [], []
    source_cell = None
    for line in grid:
        for cell in line:
            if len(cell) == 0: continue

            for v in cell: grid_map[v] = cell

            if cell[0] == source:
                source_cell = cell
                continue

            rep_names.append(cell[0])
            rep_coords.append(coords[cell[0]])
            rep_degrees.append(min(len(cell) + 1, grid_dim ** 2 - 1))

    # Make the whole instance source the representatives source also:
    rep_names.insert(0, source)
    rep_coords.insert(0, coords[source])
    rep_degrees.insert(0, min(len(source_cell), grid_dim ** 2 - 1))

    return rep_names, rep_coords, rep_degrees, grid_map


def calc_height(root, sol_dg, dist):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(dist(root, v) + calc_height(v, sol_dg, dist) for v in children)


def calc_depth(root, sol_dg):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(1 + calc_depth(v, sol_dg) for v in children)


def spectral_clustering(graph, k, dim=2):
    coords_dict = nx.spectral_layout(graph, dim=dim)
    coors_list = np.array([coords_dict[node] for node in graph.nodes])
    kmeans = KMeans(n_clusters=k, random_state=0).fit(coors_list)
    return [kmeans.labels_[node] for node in graph.nodes]
