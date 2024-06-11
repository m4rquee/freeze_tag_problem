from matplotlib import pyplot as plt

from src.cp.utils.utils import *
from src.cp.utils.plotting import plot_graph
from src.cp.utils.reading import read_tsplib_2d_graph

delta = 1E-2

# Setup:
names, coords = read_tsplib_2d_graph()
n = len(names)
source = 0
dist = L2Norm(coords, delta)

# Print instance info:
print('Freeze-Tag instance information:')
print('\tNumber of nodes =', n)
print('\tSource =', source)

# FTP solving:
sol_edges, makespan = greedy_solution(source, n, dist)

print('\nGreedy FTP solution:')
print(f'\tsolution makespan: {delta * makespan:.2f}')
tree = nx.DiGraph(sol_edges)
hop_depth = calc_depth(source, tree)
print(f'\thop depth: {hop_depth}')
source_radius = radius(source, n, dist)
print(f'\tsource radius: {delta * source_radius:.2f}')
print(f'\tgap: {100 * (makespan - source_radius) / source_radius:.2f}%')

# Solution plotting:
plt.figure(figsize=(10, 6))

node_colors = ['black' if source != node else 'red' for node in tree.nodes]
coords_dict = {names[i]: c for i, c in enumerate(coords)}
plot_graph(tree, coords_dict, node_colors, 'green', style='solid', node_size=40)

plt.gca().set_aspect('equal', adjustable='box')
plt.show()
