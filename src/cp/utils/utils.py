import random
from math import floor, ceil, log2

import numpy as np
import networkx as nx
from sklearn import cluster

from src.cp.utils.constants import *
from src.cp.utils.metric_space import *

random.seed(SEED)


def radius(center, space):
    ret = 0
    for v in range(space.n):
        if v == center: continue
        ret = max(ret, space(center, v))
    return ret


def trivial_ftp_ub(n, dist):
    minimum_depth = ceil(log2(n))  # minimum depth for a balanced wake-up tree

    diameter = 0
    for u in range(n):
        for v in range(u + 1, n):
            diameter = max(diameter, dist(u, v))
    return diameter * minimum_depth


def min_dist(space):
    ret = float('inf')
    for u in range(space.n):
        for v in range(u + 1, space.n):
            ret = min(ret, space(u, v))
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

    # Make point zero the center of its own cluster:
    centers_names.remove(cluster_map[0][0])  # remove the previous center
    centers_names.append(0)
    centers_names.sort()
    cluster_map[0].sort()

    rep_degrees = [min(len(cluster_map[0]), space.n - 1)]
    for v in centers_names[1:]:
        rep_degrees.append(min(len(cluster_map[v]) + 1, space.n - 1))

    return centers_names, rep_degrees, cluster_map


def partition(space, labels):
    parts = {l: [] for l in labels}
    for v, l in enumerate(labels):
        parts[l].append(v)

    centers_names, cluster_map = [], {}
    for l, part in parts.items():
        part.sort()
        centers_names.append(part[0])
        for v in part: cluster_map[v] = part

    # Make point zero the center of its own cluster:
    centers_names.remove(cluster_map[0][0])  # remove the previous center
    centers_names.append(0)
    centers_names.sort()
    cluster_map[0].sort()

    rep_degrees = [min(len(cluster_map[0]), space.n - 1)]
    for v in centers_names[1:]:
        rep_degrees.append(min(len(cluster_map[v]) + 1, space.n - 1))

    return centers_names, rep_degrees, cluster_map


def max_cluster_radius(space, centers_names, cluster_map):
    ret = 0
    for c in centers_names:
        for u in cluster_map[c]:
            for v in cluster_map[c]:
                ret = max(ret, space(u, v))
    return ret


def greedy_solution(source, space):
    # todo: optimize this function
    if space.n <= 1: return [], 0

    # Connect the source to the closest node:
    min_target = None
    min_arc_weight = float('inf')
    for v in range(space.n):
        if v == source: continue
        aux = space(source, v)
        if aux < min_arc_weight:
            min_target = v
            min_arc_weight = aux
    node_activation = {source: 0, min_target: min_arc_weight}
    sol_edges = [(source, min_target)]
    makespan = min_arc_weight

    # Init the degree map (-1 are not yet border nodes and -2 saturated nodes):
    degree = {v: -1 for v in range(space.n)}
    degree[source] = -2
    degree[min_target] = 0

    border = [min_target]
    while len(sol_edges) < space.n - 1:
        min_arc = (None, None)
        min_makespan = float('inf')
        for u in border:
            for v in range(space.n):
                if v == source or v == u: continue
                if degree[v] == -1:  # not yet added
                    aux = node_activation[u] + space(u, v)
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
    bounding_box_size = 1.01 * max(*xs, *ys)  # add 1% to avoid boundary points

    # Resize the points to fit in a ceil(1 / eps) x ceil(1 / eps) square:
    factor = grid_dim / bounding_box_size
    return [(factor * x, factor * y) for x, y in zip(xs, ys)], factor


def discretize(coords, eps):  # the source is assumed to be the last node
    grid_dim = ceil(1.0 / eps)
    grid = [[[] for _ in range(grid_dim)] for _ in range(grid_dim)]

    source = 0
    for v, (x, y) in enumerate(coords):
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

    _, rep_degrees, cluster_map = clusters(space, centers_names)
    return centers_names, rep_degrees, cluster_map


def spectral_clustering(space: GraphDist, eps):
    n_clusters = min(space.n, ceil(1.0 / eps))
    coords_dict = nx.spectral_layout(space.graph, dim=3)
    coors_list = [coords_dict[node] for node in space.graph.nodes]
    kmeans = cluster.KMeans(n_clusters=n_clusters).fit(coors_list)

    centers_names, rep_degrees, cluster_map = partition(space, kmeans.labels_)
    return centers_names, rep_degrees, cluster_map


def graph_to_edge_matrix(graph):
    edge_mat = np.zeros((len(graph), len(graph)), dtype=int)
    for u in graph:
        for v in graph.neighbors(u):
            edge_mat[u][v] = 1
        edge_mat[u][u] = 1
    return edge_mat


def k_means_clustering(space, eps):
    n_clusters = min(space.n, ceil(1.0 / eps))
    if isinstance(space, GraphDist):
        edge_mat = graph_to_edge_matrix(space.original_graph)
        kmeans = cluster.KMeans(n_clusters=n_clusters).fit(edge_mat)

        centers_names = []
        for label, c in enumerate(kmeans.cluster_centers_):
            cluster_nodes = np.argwhere(kmeans.labels_ == label)
            best_center_pos = np.argmax(c[cluster_nodes])
            best_center = cluster_nodes[best_center_pos]
            centers_names.append(int(best_center))

        centers_names, rep_degrees, cluster_map = clusters(space, centers_names)
    else:
        coords_array = np.array(space.coords)
        kmeans = cluster.KMeans(n_clusters=n_clusters).fit(coords_array)
        centers_names, rep_degrees, cluster_map = partition(space, kmeans.labels_)

    return centers_names, rep_degrees, cluster_map


def clusterize(space: MetricSpace, eps, method='k_means'):
    if method == 'spectral':
        return spectral_clustering(space, eps)
    elif method == 'k_center':
        return k_center_clustering(space, eps)
    elif method == 'k_means':
        return k_means_clustering(space, eps)


def calc_height(root, sol_dg, dist):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(dist(root, v) + calc_height(v, sol_dg, dist) for v in children)


def calc_depth(root, sol_dg):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(1 + calc_depth(v, sol_dg) for v in children)
