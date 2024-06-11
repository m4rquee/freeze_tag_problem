from sys import argv

from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.utils.utils import *
from src.cp.utils.instance import *
from src.cp.utils.cli import read_arguments
from src.cp.solving.solvers import solve_ftp
from src.cp.utils.plotting import plot_graph


MAX_TIME, TYPE = read_arguments(int, str)
TYPE = argv[2]

if TYPE == 'tsplib_2d':
    names, coords = read_tsplib_2d_graph()
else:
    names, edges = get_instance(TYPE)
n = len(names)
source = 0
if TYPE == 'tsplib_2d':
    dist = L2Norm(coords)
else:
    dist = GraphDist(n, edges)
sol_edges, UB = greedy_solution(source, n, dist)
source_radius = radius(source, n, dist)

# Print instance info:
print('Freeze-Tag instance information:')
print(f'\tTime limit = {MAX_TIME}s')
print('\tNumber of nodes =', n)
print('\tSource =', source)
print(f'\tSource radius = {DELTA * source_radius:.2f}')
print(f'\tGreedy bound = {DELTA * UB:.2f}')
min_edge = min_dist(n, dist)
LB = max(source_radius, min_edge * ceil(log2(n)))
print(f'\tStart lower bound = {DELTA * LB:.2f}')
print(f'\tInitial gap = {100 * (UB - LB) / LB:.2f}%')

# FTP solving:
status, model, solver, depth, d_v, x_e = solve_ftp(names, dist, MAX_TIME, LB, UB, True, sol_edges)

if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
    # FTP solution:
    print('Full FTP solution:')
    depth = solver.Value(depth)
    print(f'  solution makespan: {DELTA * depth:.2f}')
    lb = solver.BestObjectiveBound()
    print(f'  gap: {100 * (depth - lb) / lb:.2f}%')
    print(f'  greedy gap: {100 * (UB - lb) / LB:.2f}%')
    print(f'  d_v: (', end='')
    for v in range(n - 1):
        depth_v = solver.Value(d_v[v])
        print(f'{DELTA * depth_v:.2f}', end=', ')
    print(f'{DELTA * solver.Value(d_v[-1])})')

    print('  edges: ', end='')
    sol_edges = []
    for u in range(n):
        for v in range(1, n):
            if u == v: continue
            if solver.Value(x_e[u][v]):
                print(f'{names[u]}-{names[v]}; ', end='')
                sol_edges.append((names[u], names[v]))
    tree = nx.DiGraph(sol_edges)
    hop_depth = calc_depth(source, tree)
    print(f'\n  hop depth: {hop_depth}')

    node_colors = []
    for v in tree.nodes:
        if v == source:
            node_colors.append('red')
        else:
            node_colors.append('cyan' if solver.Value(d_v[v]) == depth else 'black')

    # Solution plotting:
    plt.figure(figsize=(10, 6))
    if TYPE != 'tsplib_2d':
        coords_dict = nx.nx_agraph.graphviz_layout(dist.original_graph, prog='dot')  # nx.spectral_layout(dist.original_graph, dim=2)
        plot_graph(dist.original_graph, coords_dict, 'white', 'gray', style='dotted', node_size=40)
        plot_graph(tree, coords_dict, node_colors, 'green', style='solid', node_size=40,
                      connectionstyle='arc3,rad=0.1')
    else:
        coords_dict = {names[i]: c for i, c in enumerate(coords)}
        plot_graph(tree, coords_dict, node_colors, 'green', style='solid', node_size=40)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
else:
    print('No solution found.')
