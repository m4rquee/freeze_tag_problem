from functools import cache
from math import sqrt, floor, ceil, log2


def trivial_ub(n, dist):
    minimum_depth = ceil(log2(n))

    diameter = float('-inf')
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


def normalize(coords, eps):
    grid_dim = ceil(1.0 / eps)
    xs, ys = zip(*coords)  # gather all points x and y coordinates

    # Stick the bounding box corner to the origin:
    min_x_coord = min(xs)
    min_y_coord = min(ys)
    xs = [-min_x_coord + x for x in xs]
    ys = [-min_y_coord + y for y in ys]

    # Resize the points to fit in a ceil(1 / eps) x ceil(1 / eps) square:
    factor = grid_dim / (max(xs + ys) + 1)  # add one to avoid boundary points
    return [(factor * x, factor * y) for x, y in zip(xs, ys)], factor


def discretize(names, coords, source, eps):
    num_cells = ceil(1 / eps)
    grid = [[[] for _ in range(num_cells)] for _ in range(num_cells)]
    for i, (x, y) in enumerate(coords):
        if i == source:
            grid[floor(x)][floor(y)].insert(0, i)  # make sure the source is a representative
        else:
            grid[floor(x)][floor(y)].append(i)

    disc_names, disc_coords, degrees = [], [], []
    source_name, source_coord, source_degree = None, None, None
    for line in grid:
        for cell in line:
            if len(cell) == 0:
                continue

            if cell[0] == source:
                source_name = names[cell[0]]
                source_coord = coords[cell[0]]
                source_degree = min(len(cell), num_cells ** 2 - 1)
            else:
                disc_names.append(names[cell[0]])
                disc_coords.append(coords[cell[0]])
                degrees.append(min(len(cell) + 1, num_cells ** 2 - 1))

    disc_names.append(source_name)
    disc_coords.append(source_coord)
    degrees.append(source_degree)
    source_i = len(disc_names) - 1
    return disc_names, disc_coords, degrees, source_i, grid


def calc_height(root, names_to_i, sol_dg, dist):
    children = list(sol_dg.neighbors(root))
    if len(children) == 0: return 0
    return max(dist(names_to_i[root], names_to_i[v]) + calc_height(v, names_to_i, sol_dg, dist) for v in children)
