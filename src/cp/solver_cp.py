from sys import argv

from matplotlib import pyplot as plt
from ortools.sat.python import cp_model

from src.cp.utils import *
from src.cp.reading import *
from src.cp.solvers import solve_ftp
from src.cp.plotting import plot_solution

MAX_TIME = int(argv[1])
TYPE = argv[2]

delta = 1E-2

# Setup:
if TYPE == 'tsplib_2d':
    names, coords = read_tsplib_2d_graph()
elif TYPE == 'tsplib_hcp':
    names, edges = read_tsplib_hcp_graph()
elif 'gnp' in TYPE:
    names, edges = gnp_graph(*TYPE.split('-')[1::])
elif 'ws' in TYPE:
    names, edges = ws_graph(*TYPE.split('-')[1::])
else:  # if TYPE == 'dig':
    names, edges = read_dig_graph()

n = len(names)
source = names[n - 1]
names_to_i = {name: i for i, name in enumerate(names)}
if TYPE == 'tsplib_2d':
    dist = l2_norm(coords, delta)
else:
    dist = graph_dist(edges, names, delta)
sol_edges, UB = greedy_solution(n - 1, n, dist, names)
source_radius = radius(n, n - 1, dist)

# Print instance info:
print('Freeze-Tag instance information:')
print(f'\tTime limit = {MAX_TIME}s')
print('\tNumber of nodes =', n)
print('\tSource =', source)
print(f'\tSource radius = {delta * source_radius:.2f}')
print(f'\tGreedy bound = {delta * UB:.2f}')
min_edge = min_dist(n, dist)
LB = max(source_radius, min_edge * ceil(log2(n)))
print(f'\tInitial gap = {100 * (UB - LB) / LB:.2f}%')

# FTP solving:
status, model, solver, depth, d_v, x_e = solve_ftp(names, dist, MAX_TIME, LB, UB, True, sol_edges)

if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
    # FTP solution:
    print('Full FTP solution:')
    depth = solver.Value(depth)
    print(f'  solution makespan: {delta * depth:.2f}')
    lb = solver.BestObjectiveBound()
    print(f'  gap: {100 * (depth - lb) / lb:.2f}%')
    print(f'  d_v: (', end='')
    for v in range(n - 1):
        depth_v = solver.Value(d_v[v])
        print(f'{delta * depth_v:.2f}', end=', ')
    print(f'{delta * solver.Value(d_v[-1])})')

    print('  edges: ', end='')
    sol_edges = []
    for u in range(n):
        for v in range(n - 1):
            if u == v: continue
            if solver.Value(x_e[u][v]):
                print(f'{names[u]}-{names[v]}; ', end='')
                sol_edges.append((names[u], names[v]))
    tree = nx.DiGraph(sol_edges)
    hop_depth = calc_depth(source, names_to_i, tree)
    print(f'\n  hop depth: {hop_depth}')

    node_colors = []
    for v in tree.nodes:
        if v == source:
            node_colors.append('red')
        else:
            node_colors.append('cyan' if solver.Value(d_v[names_to_i[v]]) == depth else 'black')

    # Solution plotting:
    plt.figure(figsize=(10, 6))
    if TYPE != 'tsplib_2d':
        whole_graph = nx.Graph(edges)
        coords_dict = nx.nx_agraph.graphviz_layout(whole_graph, prog='dot')
        plot_solution(whole_graph, edges, coords_dict, names, 'white', 'gray', style='dotted', node_size=40)
        plot_solution(tree, sol_edges, coords_dict, names, node_colors, 'green', style='solid', node_size=40,
                      connectionstyle='arc3,rad=0.1')
    else:
        coords_dict = {names[i]: c for i, c in enumerate(coords)}
        plot_solution(tree, sol_edges, coords_dict, names, node_colors, 'green', style='solid', node_size=40)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
else:
    print('No solution found.')
