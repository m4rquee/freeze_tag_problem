from ortools.sat.python import cp_model

from src.cp.utils.utils import *
from src.cp.utils.instance import *
from src.cp.utils.cli import read_arguments
from src.cp.utils.plotting import plot_solution
from src.cp.solving.solvers import get_solution, solve_ftp

# Setup:
MAX_TIME, TYPE = read_arguments(int, str)
space = get_instance(TYPE)
source = 0
sol_edges, UB = greedy_solution(source, space)
source_radius = radius(source, space)

# Print instance info:
print('Freeze-Tag instance information:')
print(f'\tTime limit = {MAX_TIME}s')
print('\tNumber of nodes =', space.n)
print('\tSource =', source)
print(f'\tSource radius = {DELTA * source_radius:.2f}')
print(f'\tGreedy bound = {DELTA * UB:.2f}')
min_edge = min_dist(space)
LB = max(source_radius, min_edge * ceil(log2(space.n)))
print(f'\tLower bound = {DELTA * LB:.2f}')
print(f'\tInitial gap = {100 * (UB - LB) / LB:.2f}%')

# FTP solving:
names = list(space)
bdmhst = solve_ftp(names, space, MAX_TIME, LB, UB, True, sol_edges)

if bdmhst.status == cp_model.FEASIBLE or bdmhst.status == cp_model.OPTIMAL:
    solution = get_solution(space, bdmhst)
    plot_solution('dot', space, solution)
else:
    print('No solution found.')
