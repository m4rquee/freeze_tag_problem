from sys import argv

import networkx as nx
from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.reading import read_tsplib_graph
from src.cp.solvers import solve_ftp, solve_bdhst
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
dist = l2_norm(coords, delta)
UB = trivial_ub(n, dist)

# Full FTP solving:
status, model, solver, depth, d_v, x_e = solve_ftp(names, dist, MAX_TIME, UB)

# Discretized BDHST solving:
names_d, coords_d, degrees_d, source_i = discretize(names, coords, names_to_i[source], EPS)
n_d = len(names_d)
dist_d = l2_norm(coords_d, delta)
UB_d = trivial_ub(n_d, dist_d)
status_d, model_d, solver_d, depth_d, d_v_d, x_e_d = solve_bdhst(names_d, dist_d, degrees_d, MAX_TIME, UB_d)

status = status == cp_model.FEASIBLE or status == cp_model.OPTIMAL
status_d = status_d == cp_model.FEASIBLE or status_d == cp_model.OPTIMAL
if status and status_d:
    to_orig = delta / factor

    # Full FTP solution:
    print('Full FTP solution:')
    depth = solver.Value(depth)
    print(f'  solution depth: {to_orig * depth:.2f}')
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

    # Discretized BDHST solution:
    print('\n\nDiscretized BDHST solution:')
    depth_d = solver_d.Value(depth_d)
    print(f'  solution solution depth: {to_orig * depth_d:.2f}')
    print(f'  d_v = (', end='')
    node_colors_d = []
    for v in range(n_d - 1):
        depth_v = solver_d.Value(d_v_d[v])
        node_colors_d.append('cyan' if depth_v == depth_d else 'white')
        print(f'{to_orig * depth_v:.2f}', end=', ')
    print(f'{to_orig * solver_d.Value(d_v_d[-1])})')
    node_colors_d.append('red')

    print('  edges: ', end='')
    sol_edges_d = []
    for u in range(n_d):
        for v in range(n_d - 1):
            if u == v: continue
            if solver_d.Value(x_e_d[u][v]):
                print(f'{names_d[u]}-{names_d[v]}; ', end='')
                sol_edges_d.append((names_d[u], names_d[v]))

    # Solution plotting:
    plt.figure(figsize=(10, 6))

    coords_dict = {names[i]: c for i, c in enumerate(coords)}
    plot_solution(DG, sol_edges, coords_dict, names, node_colors, 'red', style='dashed', node_size=40)

    coords_dict_d = {names_d[i]: c for i, c in enumerate(coords_d)}
    plot_solution(DG, sol_edges_d, coords_dict_d, names_d, node_colors_d, 'green', style='dotted', node_size=10)

    plot_grid(EPS)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
else:
    print('No solution found.')
