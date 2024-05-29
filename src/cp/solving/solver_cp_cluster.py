from sys import argv

from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.utils import *
from src.cp.reading import *
from src.cp.plotting import plot_solution
from src.cp.solving.solvers import solve_bdhst, solve_ftp_inner

MAX_TIME = int(argv[1])
EPS = float(argv[2])
TYPE = argv[3]
TOTAL_TIME = MAX_TIME
delta = 1E-4

# Setup:
if TYPE == 'tsplib_hcp':
    names, edges = read_tsplib_hcp_graph()
elif 'gnp' in TYPE:
    names, edges = gnp_graph(*TYPE.split('-')[1::])
elif 'ws' in TYPE:
    names, edges = ws_graph(*TYPE.split('-')[1::])
elif 'tree' in TYPE:
    names, edges = tree_graph(*TYPE.split('-')[1::])
elif 'regular' in TYPE:
    names, edges = regular_graph(*TYPE.split('-')[1::])
else:  # if TYPE == 'dig':
    names, edges = read_dig_graph()

n = len(names)
source = 0
space = GraphDist(edges, delta)  # the underlying graph metric space

# Print instance info:
print('Freeze-Tag instance information:')
print(f'\tTime limit = {MAX_TIME}s')
print('\tNumber of nodes =', n)
print('\tSource =', source)

# Upper level discretized BDHST solving:
d_names, d_degrees, cluster_map = clusterize(space, EPS)
d_n = len(d_names)
d_space = space[d_names]
d_sol_edges, d_UB = greedy_solution(source, d_n, d_space)
hop_depth = 0.0  # min(d_n - 1, ceil(log2(d_n) ** 2))

# Print upper level instance info:
print('Upper level discretized BDHST information:')
print('\tNumber of nodes =', d_n)
print('\tHop depth =', hop_depth if hop_depth > 0 else "no limit")

d_status, _, d_solver, d_depth, d_d_v, d_x_e = \
    solve_bdhst(d_names, d_space, d_degrees, MAX_TIME, 0, d_UB, hop_depth, True, 'Freeze-Tag Problem', 0, d_sol_edges)

d_status = d_status == cp_model.FEASIBLE or d_status == cp_model.OPTIMAL
if not d_status:
    exit('Could not find any solution to the upper level discretization!')

# Upper level discretized BDHST solution:
print('Upper level discretized BDHST solution:')
d_depth = d_solver.Value(d_depth)
print(f'  number of nodes: {d_n}')
print(f'  solution height: {delta * d_depth:.2f}')

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
_, MAX_TIME = solve_ftp_inner(sol_edges, d_tree, source, space, cluster_map, delta, MAX_TIME)
print(f'Finished solving all inner cell problems with {MAX_TIME:.1f}s out of {TOTAL_TIME:.1f}s remaining...')
tree = nx.DiGraph(sol_edges)  # full solution tree

makespan = calc_height(source, tree, space)
hop_depth = calc_depth(source, tree)
source_radius = radius(source, n, space)
biggest_cell = max_cluster_radius(space, d_names, cluster_map)
d_lb = d_depth - d_hop_depth * biggest_cell  # todo: only count the cells that were actually expanded
lb = max(source_radius, d_lb)

# Final solution:
print('\nFreeze-Tag solution:')
print(f'  number of nodes  : {n}')
print(f'  solution makespan: {delta * makespan:.2f}')
print(f'  source radius    : {delta * source_radius:.2f}')
print(f'  max cluster size : {delta * biggest_cell:.2f}')
print(f'  d_lower bound    : {delta * d_lb:.2f}')
print(f'  gap              : {100 * (makespan - lb) / lb:.2f}%')
print(f'  hop depth        : {hop_depth}')
print(f'  time to solve    : {TOTAL_TIME - MAX_TIME:.2f}s')

# Solution plotting:
plt.figure(figsize=(10, 6))

whole_graph = nx.empty_graph(n)
whole_graph.add_edges_from(edges)
coords_dict = nx.spring_layout(whole_graph)

# Full solution:
plot_solution(tree, coords_dict, 'White', 'green', style='solid', node_size=100,
              connectionstyle='arc3,rad=0.15', with_labels=True, verticalalignment='bottom')

# Domain graph:
node_colors = [cluster_map[v][0] for v in range(n)]
plot_solution(whole_graph, coords_dict, node_colors, 'black', style='dotted', node_size=100)

# Upper level solution:
node_colors = ['white' if source != node else 'red' for node in d_names]
plot_solution(d_tree, coords_dict, node_colors, 'red', style='dotted', node_size=25,
              connectionstyle='arc3,rad=0.1')

plt.gca().set_aspect('equal', adjustable='box')
plt.show()
