from sys import argv
from math import ceil

import networkx as nx
from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.reading import read_tsplib_graph
from src.cp.solvers import solve_ftp, solve_bdhst
from src.cp.plotting import plot_solution, plot_grid
from src.cp.utils import l2_norm, normalize, discretize

EPS = float(argv[1])
MAX_TIME = int(argv[2])

delta = 1E-4

# Setup:
names, coords = read_tsplib_graph()
coords = normalize(coords, ceil(1.0 / EPS))
n = len(names)
source = n - 1

DG = nx.complete_graph(n, nx.DiGraph)
dist = l2_norm(coords, delta)

# FTP solving:
status, solver, depth, d_v, x_e = solve_ftp(names, source, dist, MAX_TIME)

names_h, coords_h, degrees, source_i = discretize(names, coords, source, EPS)
n_h = len(names_h)
dist_h = l2_norm(coords_h, delta)
status_h, solver_h, depth_h, d_v_h, x_e_h = solve_bdhst(names_h, source_i, dist_h, degrees, MAX_TIME)

if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
    print('Found a feasible solution.')
    depth = solver.Value(depth)
    print(f'  solution depth: {delta * depth:.2f}')
    print(f'  d_v = (', end='')
    node_colors = []
    for v in range(n - 1):
        depth_v = solver.Value(d_v[v])
        node_colors.append('cyan' if depth_v == depth else 'gray')
        print(f'{delta * depth_v:.2f}', end=', ')
    print(f'{delta * solver.Value(d_v[source])})')
    node_colors.append('magenta')

    print('  edges: ', end='')
    sol_edges = []
    for u in range(n):
        for v in range(n - 1):
            if u == v:
                continue
            if solver.Value(x_e[u][v]):
                print(f'{names[u]}-{names[v]}; ', end='')
                sol_edges.append((u, v))

    depth_h = solver_h.Value(depth_h)
    print(f'\n\n  Higher solution depth: {delta * depth_h:.2f}')
    print(f'  d_v = (', end='')
    node_colors_h = []
    for v in range(n_h - 1):
        depth_v = solver_h.Value(d_v_h[v])
        node_colors_h.append('blue' if depth_v == depth_h else 'black')
        print(f'{delta * depth_v:.2f}', end=', ')
    print(f'{delta * solver_h.Value(d_v[source_i])})')
    node_colors_h.append('red')

    print('  edges: ', end='')
    sol_edges_h = []
    for u in range(n_h):
        for v in range(n_h - 1):
            if u == v:
                continue
            if solver_h.Value(x_e_h[u][v]):
                print(f'{names_h[u]}-{names_h[v]}; ', end='')
                sol_edges_h.append((int(names_h[u]) - 1, int(names_h[v]) - 1))

    plt.figure(figsize=(10, 6))
    plot_solution(DG, sol_edges, coords, names, node_colors, 'red')
    coords_h = {int(names_h[i]) - 1: v for i, v in enumerate(coords_h)}
    # names_h = {int(names_h[i]) - 1: v for i, v in enumerate(names_h)}
    plot_solution(DG, sol_edges_h, coords_h, names_h, node_colors_h, 'green', style='dotted')

    plot_grid(EPS)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
else:
    print('No solution found.')

# Statistics:
print('\nStatistics')
print(f'  status   : {solver.StatusName(status)}')
print(f'  conflicts: {solver.NumConflicts()}')
print(f'  branches : {solver.NumBranches()}')
print(f'  wall time: {solver.WallTime()} s')
