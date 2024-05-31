import random
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


def random_furthest_point(centers, space):
    dist_to_center = {v: float('inf') for v in space if v not in centers}

    max_dist = 0
    for v in dist_to_center.keys():
        for u in centers:
            dist_to_center[v] = min(dist_to_center[v], space(u, v))
        max_dist = max(max_dist, dist_to_center[v])

    candidates = [v for v in dist_to_center.keys() if dist_to_center[v] == max_dist]
    return random.choices(candidates)[0]


def random_closest_center(point, centers, space):
    dist_to_centers = {v: space(v, point) for v in centers}
    min_center_dist = min(dist_to_centers.values())
    candidates = [v for v in dist_to_centers.keys() if dist_to_centers[v] == min_center_dist]
    return random.choices(candidates)[0]


def clusters(space, centers_names):
    cluster_map = {v: [v] for v in centers_names}
    for v in space:
        if v not in centers_names:
            close_center = random_closest_center(v, centers_names, space)
            cluster_map[close_center].append(v)
            cluster_map[v] = cluster_map[close_center]
    rep_degrees = [len(cluster_map[v]) for v in centers_names]

    return rep_degrees, cluster_map


def partition(space, labels):
    parts = {l: [] for l in labels}
    for v, l in enumerate(labels):
        parts[l].append(v)

    centers_names, rep_degrees, cluster_map = [], [], {}
    for l, part in parts.items():
        part.sort()
        centers_names.append(part[0])
        rep_degrees.append(len(part))
        for v in part: cluster_map[v] = part
    return centers_names, rep_degrees, cluster_map


def max_cluster_radius(space, centers_names, cluster_map):
    ret = 0
    for c in centers_names:
        for u in cluster_map[c]:
            for v in cluster_map[c]:
                ret = max(ret, space(u, v))
    return ret


class MetricSpace:
    def __init__(self, n):
        self.n = n
        self.func = None

    def __call__(self, u, v):
        return self.func(u, v)

    def __len__(self):
        return self.n

    def __getitem__(self, points):  # restrict to a sub-set of points
        ret = MetricSpace(len(points))

        def aux(u, v):
            u, v = points[u], points[v]
            return self.func(u, v)

        ret.func = aux
        return ret

    def __iter__(self):
        yield from range(self.n)


class L2Norm(MetricSpace):
    def __init__(self, coords, delta=DELTA):
        super().__init__(len(coords))
        self.coords = coords
        self.delta = delta

        @cache
        def aux(u, v):
            d = sqrt((self.coords[u][0] - self.coords[v][0]) ** 2 + (self.coords[u][1] - self.coords[v][1]) ** 2)
            return floor(d / self.delta)  # scale up so we can treat rational values as integers

        self.func = aux


class GraphDist(MetricSpace):
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


def k_center_clustering(space: MetricSpace, eps):
    source = 0
    n_clusters = min(space.n, ceil(1.0 / eps))

    centers_names = [source]
    for _ in range(n_clusters - 1):
        center = random_furthest_point(centers_names, space)
        centers_names.append(center)

    rep_degrees, cluster_map = clusters(space, centers_names)
    return centers_names, rep_degrees, cluster_map


def spectral_clustering(space: GraphDist, eps):
    n_clusters = min(space.n, ceil(1.0 / eps))
    coords_dict = nx.spectral_layout(space.graph, dim=3)
    coors_list = [coords_dict[node] for node in space.graph.nodes]
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(coors_list)

    centers_names, rep_degrees, cluster_map = partition(space, kmeans.labels_)
    return centers_names, rep_degrees, cluster_map


def clusterize(space: MetricSpace, eps, method='k_center'):
    if method == 'spectral':
        return spectral_clustering(space, eps)
    elif method == 'k_center':
        return k_center_clustering(space, eps)


def calc_height(root, sol_dg, dist):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(dist(root, v) + calc_height(v, sol_dg, dist) for v in children)


def calc_depth(root, sol_dg):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(1 + calc_depth(v, sol_dg) for v in children)
