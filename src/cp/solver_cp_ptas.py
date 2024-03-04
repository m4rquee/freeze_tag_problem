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

# Upper level discretized BDHST solving:
d_names, d_coords, d_degrees, d_source_i, d_grid = discretize(names, coords, names_to_i[source], EPS)
d_n = len(d_names)
d_dist = l2_norm(d_coords, delta)
d_UB = trivial_ub(d_n, d_dist)
hop_depth = min(d_n - 1, ceil(log2(d_n)) ** 2)
d_status, d_model, d_solver, d_depth, d_d_v, d_x_e = solve_bdhst(d_names, d_dist, d_degrees, MAX_TIME, d_UB, hop_depth)

d_status = d_status == cp_model.FEASIBLE or d_status == cp_model.OPTIMAL
if not d_status:
    exit('Could not find any solution to the upper level discretization!')

to_orig = delta / factor

# Upper level discretized BDHST solution:
print('\n\nUpper level discretized BDHST solution:')
d_depth = d_solver.Value(d_depth)
print(f'  number of nodes        : {d_n}')
print(f'  solution solution depth: {to_orig * d_depth:.2f}')
print(f'  d_v = (', end='')
node_colors = []
for v in range(d_n - 1):
    depth_v = d_solver.Value(d_d_v[v])
    node_colors.append('cyan' if depth_v == d_depth else 'black')
    print(f'{to_orig * depth_v:.2f}', end=', ')
print(f'{to_orig * d_solver.Value(d_d_v[-1])})')
node_colors.append('red')

print('  edges: ', end='')
sol_edges = []
for u in range(d_n):
    for v in range(d_n - 1):
        if u == v: continue
        if d_solver.Value(d_x_e[u][v]):
            print(f'{d_names[u]}-{d_names[v]}; ', end='')
            sol_edges.append((d_names[u], d_names[v]))

# Solution plotting:
plt.figure(figsize=(10, 6))

coords_dict = {d_names[i]: c for i, c in enumerate(d_coords)}
plot_solution(DG, sol_edges, coords_dict, d_names, node_colors, 'green', style='dotted', node_size=10)

plot_grid(EPS)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
