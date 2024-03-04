from sys import argv
from math import ceil, log2

import networkx as nx
from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.solvers import solve_bdhst
from src.cp.reading import read_tsplib_graph
from src.cp.plotting import plot_solution, plot_grid
from src.cp.utils import trivial_ub, l2_norm, normalize, discretize

EPS = float(argv[1])
MAX_TIME = int(argv[2])

delta = 1E-4

# Setup:
names, coords = read_tsplib_graph()
coords, factor = normalize(coords, EPS)
n = len(names)
source = names[n - 1]
names_to_i = {name: i for i, name in enumerate(names)}
DG = nx.complete_graph(names, nx.DiGraph)

# Discretized BDHST solving:
names, coords, degrees, source_i = discretize(names, coords, names_to_i[source], EPS)
n = len(names)
dist = l2_norm(coords, delta)
UB = trivial_ub(n, dist)
hop_depth = min(n - 1, ceil(log2(n)) ** 2)
status, model, solver, depth, d_v, x_e = solve_bdhst(names, dist, degrees, MAX_TIME, UB, hop_depth)

status = status == cp_model.FEASIBLE or status == cp_model.OPTIMAL
if status:
    to_orig = delta / factor

    # Discretized BDHST solution:
    print('\n\nDiscretized BDHST solution:')
    depth = solver.Value(depth)
    print(f'  number of nodes        : {n}')
    print(f'  solution solution depth: {to_orig * depth:.2f}')
    print(f'  d_v = (', end='')
    node_colors = []
    for v in range(n - 1):
        depth_v = solver.Value(d_v[v])
        node_colors.append('cyan' if depth_v == depth else 'black')
        print(f'{to_orig * depth_v:.2f}', end=', ')
    print(f'{to_orig * solver.Value(d_v[-1])})')
    node_colors.append('red')

    print('  edges: ', end='')
    sol_edges = []
    for u in range(n):
        for v in range(n - 1):
            if u == v: continue
            if solver.Value(x_e[u][v]):
                print(f'{names[u]}-{names[v]}; ', end='')
                sol_edges.append((names[u], names[v]))

    # Solution plotting:
    plt.figure(figsize=(10, 6))

    coords_dict = {names[i]: c for i, c in enumerate(coords)}
    plot_solution(DG, sol_edges, coords_dict, names, node_colors, 'green', style='dotted', node_size=10)

    plot_grid(EPS)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
else:
    print('No solution found.')
