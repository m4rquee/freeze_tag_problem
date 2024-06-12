from ortools.sat.python import cp_model

from src.cp.utils.utils import *
from src.cp.utils.cli import read_arguments
from src.cp.utils.instance import get_instance
from src.cp.utils.plotting import plot_solution
from src.cp.solving.solvers import Solution, solve_bdmhst, solve_ftp_inner


# Setup:
MAX_TIME, EPS, TYPE = read_arguments(int, float, str)
TOTAL_TIME = MAX_TIME
space = get_instance(TYPE)
source = 0

# Print instance info:
print('Freeze-Tag instance information:')
print(f'\tTime limit = {MAX_TIME}s')
print('\tNumber of nodes =', space.n)
if isinstance(space, GraphDist):
    print('\tNumber of edges =', space.original_graph.number_of_edges())
print('\tSource =', source)

# Upper level discretized BDMHST solving:
d_names, d_degrees, cluster_map = clusterize(space, EPS)
d_n = len(d_names)
d_space = space[d_names]
d_sol_edges, d_UB = greedy_solution(source, d_space)
hop_depth = 0.0  # min(d_n - 1, ceil(log2(d_n) ** 2))

# Print upper level instance info:
print('Upper level discretized BDMHST information:')
print('\tNumber of nodes =', d_n)
print('\tHop depth =', hop_depth if hop_depth > 0 else "no limit")

d_bdmhst = solve_bdmhst(d_names, d_space, d_degrees, MAX_TIME, 0, d_UB, hop_depth, True, 'Freeze-Tag Problem', 0, d_sol_edges)

if d_bdmhst.status != cp_model.FEASIBLE and d_bdmhst.status != cp_model.OPTIMAL:
    exit('Could not find any solution to the upper level discretization!')

# Upper level discretized BDMHST solution:
print('Upper level discretized BDMHST solution:')
d_depth = d_bdmhst.solver.Value(d_bdmhst.depth)
print(f'  number of nodes: {d_n}')
print(f'  solution height: {DELTA * d_depth:.2f}')

d_sol_edges = []
for u in range(d_n):
    for v in range(1, d_n):
        if u == v: continue
        if d_bdmhst.solver.Value(d_bdmhst.x_e[u][v]):
            d_sol_edges.append((d_names[u], d_names[v]))
d_tree = nx.DiGraph(d_sol_edges)  # upper level tree
d_hop_depth = calc_depth(source, d_tree)
print(f'  hop depth      : {d_hop_depth}\n')

# Inner FTPs solving:
MAX_TIME -= d_bdmhst.solver.WallTime()
sol_edges = []
_, MAX_TIME = solve_ftp_inner(sol_edges, d_tree, source, space, cluster_map, MAX_TIME)
print(f'Finished solving all inner cell problems with {MAX_TIME:.1f}s out of {TOTAL_TIME:.1f}s remaining...')
tree = nx.DiGraph(sol_edges)  # full solution tree

makespan = calc_height(source, tree, space)
hop_depth = calc_depth(source, tree)
source_radius = radius(source, space)
biggest_cell = max_cluster_radius(space, d_names, cluster_map)
d_lb = d_depth - d_hop_depth * biggest_cell  # todo: only count the cells that were actually expanded
min_edge = min_dist(space)
lb = max(source_radius, d_lb, min_edge * ceil(log2(space.n)))

# Final solution:
_, UB = greedy_solution(source, space)
print('\nFreeze-Tag solution:')
print(f'  number of nodes  : {space.n}')
if isinstance(space, GraphDist):
    print('\tnumber of edges =', space.original_graph.number_of_edges())
print(f'  solution makespan: {DELTA * makespan:.2f}')
print(f'  source radius    : {DELTA * source_radius:.2f}')
print(f'  max cluster size : {DELTA * biggest_cell:.2f}')
print(f'  d_lower bound    : {DELTA * d_lb:.2f}')
print(f'  final lb bound   : {DELTA * lb:.2f}')
print(f'  gap              : {100 * (makespan - lb) / lb:.2f}%')
print(f'  greedy gap       : {100 * (UB - lb) / lb:.2f}%')
print(f'  hop depth        : {hop_depth}')
print(f'  time to solve    : {TOTAL_TIME - MAX_TIME:.2f}s')

# Solution plotting:
node_colors = [cluster_map[v][0] for v in space]
d_node_colors = ['white' if source != node else 'red' for node in d_names]

upper_solution = Solution(d_tree, d_node_colors)
solution = Solution(tree, node_colors)
plot_solution('dot', space, solution, upper_solution)
