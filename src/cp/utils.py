from functools import cache
from math import sqrt, floor, ceil


def l2_norm(coords, delta):
    def aux(u, v):
        d = sqrt((coords[u][0] - coords[v][0]) ** 2 + (coords[u][1] - coords[v][1]) ** 2)
        return floor(d / delta)  # scale up so we can treat rational values as integers

    return aux


def normalize(coords, num_cells):
    xs, ys = zip(*coords)
    multiplier = num_cells / (max(xs + ys) + 1)
    return [(multiplier * x, multiplier * y) for x, y in coords]


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
    return disc_names, disc_coords, degrees, source_i
