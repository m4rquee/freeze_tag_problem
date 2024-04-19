from functools import cache
from math import sqrt, floor, ceil, log2

def radius(n, center, dist):
    ret = 0
    for v in range(n):
        ret = max(ret, dist(center, v))
    return ret


def trivial_ub(n, dist):
    minimum_depth = ceil(log2(n))

    diameter = 0
    for u in range(n):
        for v in range(n):
            if u != v:
                diameter = max(diameter, dist(u, v))
    return diameter * minimum_depth


def l2_norm(coords, delta):
    @cache
    def aux(u, v):
        d = sqrt((coords[u][0] - coords[v][0]) ** 2 + (coords[u][1] - coords[v][1]) ** 2)
        return floor(d / delta)  # scale up so we can treat rational values as integers

    return aux


def greedy_solution(source, n, dist, names):
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
    sol_edges = [(names[source], names[min_target])]
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
        sol_edges.append((names[u], names[v]))
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


def discretize(names, name_to_i, coords, eps):  # the source is assumed to be the last node
    grid_dim = ceil(1.0 / eps)
    grid = [[[] for _ in range(grid_dim)] for _ in range(grid_dim)]

    source = names[-1]
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
            rep_coords.append(coords[name_to_i[cell[0]]])
            rep_degrees.append(min(len(cell) + 1, grid_dim ** 2 - 1))

    # Make the whole instance source the representatives source also:
    rep_names.append(source)
    rep_coords.append(coords[-1])
    rep_degrees.append(min(len(source_cell), grid_dim ** 2 - 1))

    return rep_names, rep_coords, rep_degrees, grid_map


def calc_height(root, names_to_i, sol_dg, dist):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(dist(names_to_i[root], names_to_i[v]) + calc_height(v, names_to_i, sol_dg, dist) for v in children)


def calc_depth(root, names_to_i, sol_dg):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(1 + calc_depth(v, names_to_i, sol_dg) for v in children)
