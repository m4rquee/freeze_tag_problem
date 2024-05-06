from matplotlib import pyplot as plt

from src.cp.utils import *
from src.cp.plotting import plot_solution
from src.cp.reading import read_tsplib_2d_graph

delta = 1E-2

# Setup:
names, coords = read_tsplib_2d_graph()
n = len(names)
source = names[n - 1]
names_to_i = {name: i for i, name in enumerate(names)}
dist = l2_norm(coords, delta)

# Print instance info:
print('Freeze-Tag instance information:')
print('\tNumber of nodes =', n)
print('\tSource =', source)

# FTP solving:
sol_edges, makespan = greedy_solution(n - 1, n, dist, names)

print('\nGreedy FTP solution:')
print(f'\tsolution makespan: {delta * makespan:.2f}')
tree = nx.DiGraph(sol_edges)
hop_depth = calc_depth(source, names_to_i, tree)
print(f'\thop depth: {hop_depth}')
source_radius = radius(n, n - 1, dist)
print(f'\tsource radius: {delta * source_radius:.2f}')
print(f'\tgap: {100 * (makespan - source_radius) / source_radius:.2f}%')

# Solution plotting:
plt.figure(figsize=(10, 6))

node_colors = ['black' if source != node else 'red' for node in tree.nodes]
coords_dict = {names[i]: c for i, c in enumerate(coords)}
plot_solution(tree, sol_edges, coords_dict, names, node_colors, 'green', style='solid', node_size=40)

plt.gca().set_aspect('equal', adjustable='box')
plt.show()
