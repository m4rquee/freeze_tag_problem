from src.cp.utils.utils import *
from src.cp.utils.cli import read_arguments
from src.cp.solving.solvers import Solution
from src.cp.utils.instance import get_instance
from src.cp.utils.plotting import plot_solution

# Setup:
TYPE, = read_arguments(str)
space = get_instance(TYPE)
source = 0

# Print instance info:
print('Freeze-Tag instance information:')
print('\tNumber of nodes =', space.n)
print('\tSource =', source)

# FTP solving:
sol_edges, makespan = greedy_solution(source, space)

print('\nGreedy FTP solution:')
print(f'\tsolution makespan: {DELTA * makespan:.2f}')
tree = nx.DiGraph(sol_edges)
hop_depth = calc_depth(source, tree)
print(f'\thop depth: {hop_depth}')
source_radius = radius(source, space)
print(f'\tsource radius: {DELTA * source_radius:.2f}')
print(f'\tgap: {100 * (makespan - source_radius) / source_radius:.2f}%')

# Solution plotting:
node_colors = ['black' if source != node else 'red' for node in tree.nodes]
plot_solution('dot', space, Solution(tree, node_colors))
