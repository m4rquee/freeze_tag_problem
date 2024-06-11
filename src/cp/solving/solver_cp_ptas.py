from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.utils.utils import *
from src.cp.utils.cli import read_arguments
from src.cp.utils.instance import read_tsplib_2d_graph
from src.cp.utils.plotting import plot_graph, plot_grid
from src.cp.solving.solvers import solve_bdhst, solve_ftp_inner


MAX_TIME, EPS = read_arguments(int, float)
TOTAL_TIME = MAX_TIME
delta = 1E-4

# Setup:
names, coords = read_tsplib_2d_graph()
coords, factor = normalize(coords, EPS)
n = len(names)
source = 0

# Print instance info:
print('Freeze-Tag instance information:')
print(f'\tTime limit = {MAX_TIME}s')
print('\tNumber of nodes =', n)
print('\tSource =', source)

# Upper level discretized BDHST solving:
d_names, d_coords, d_degrees, grid_map = discretize(names, coords, EPS)
d_n = len(d_names)
d_dist = L2Norm(d_coords, delta)
d_sol_edges, d_UB = greedy_solution(source, d_n, d_dist)
hop_depth = 0.0  # min(d_n - 1, ceil(log2(d_n) ** 2))

# Print upper level instance info:
print('Upper level discretized BDHST information:')
print('\tNumber of nodes =', d_n)
print('\tHop depth =', hop_depth if hop_depth > 0 else "no limit")

d_status, _, d_solver, d_depth, d_d_v, d_x_e = \
    solve_bdhst(d_names, d_dist, d_degrees, MAX_TIME, 0, d_UB, hop_depth, True, 'Freeze-Tag Problem', 0, d_sol_edges)

d_status = d_status == cp_model.FEASIBLE or d_status == cp_model.OPTIMAL
if not d_status:
    exit('Could not find any solution to the upper level discretization!')

to_orig = delta / factor

# Upper level discretized BDHST solution:
print('Upper level discretized BDHST solution:')
d_depth = d_solver.Value(d_depth)
print(f'  number of nodes: {d_n}')
print(f'  solution height: {to_orig * d_depth:.2f}')

d_sol_edges = []
for u in range(d_n):
    for v in range(1, d_n):
        if u == v: continue
        if d_solver.Value(d_x_e[u][v]):
            d_sol_edges.append((d_names[u], d_names[v]))
d_tree = nx.DiGraph(d_sol_edges)  # upper level tree
d_hop_depth = calc_depth(source, d_tree)
print(f'  hop depth      : {d_hop_depth}\n')

# Inner FTPs solving:
MAX_TIME -= d_solver.WallTime()
sol_edges = []
dist = L2Norm(coords, delta)
_, MAX_TIME = solve_ftp_inner(sol_edges, d_tree, source, dist, grid_map, delta, MAX_TIME)
print(f'Finished solving all inner cell problems with {MAX_TIME:.1f}s out of {TOTAL_TIME:.1f}s remaining...')
tree = nx.DiGraph(sol_edges)  # full solution tree

makespan = calc_height(source, tree, dist)
hop_depth = calc_depth(source, tree)
source_radius = radius(source, n, dist)
biggest_cell = max_cluster_radius(dist, d_names, grid_map)
d_lb = d_depth - d_hop_depth * biggest_cell  # todo: only count the cells that were actually expanded
lb = max(source_radius, d_lb)

# Final solution:
print('\nFreeze-Tag solution:')
print(f'  number of nodes  : {n}')
print(f'  solution makespan: {to_orig * makespan:.2f}')
print(f'  source radius    : {to_orig * source_radius:.2f}')
print(f'  d_lower bound    : {to_orig * d_lb:.2f}')
print(f'  gap              : {100 * (makespan - lb) / lb:.2f}%')
print(f'  hop depth        : {hop_depth}')
print(f'  time to solve    : {TOTAL_TIME - MAX_TIME:.2f}s')

# Solution plotting:
plt.figure(figsize=(8, 6))

coords_dict = {names[i]: c for i, c in enumerate(coords)}
# sol_edges = [(i, j) for i, j in sol_edges if coords[i] != coords[j]] # hide self loops
node_colors = ['black' if source != node else 'red' for node in tree.nodes]
plot_graph(tree, sol_edges, coords_dict, node_colors, 'green', style='solid', node_size=40)
plot_graph(d_tree, d_sol_edges, coords_dict, 'white', 'gray', style='dotted', node_size=10,
              connectionstyle='arc3,rad=0.1')

plot_grid(EPS)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
