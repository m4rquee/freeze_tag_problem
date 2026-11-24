import argparse
from math import ceil, log2

from ortools.sat.python import cp_model

from src.cp.solving.solvers import gap_percentage, get_solution, solve_ftp
from src.cp.utils.constants import DELTA
from src.cp.utils.instance import get_instance
from src.cp.utils.plotting import plot_solution
from src.cp.utils.utils import greedy_solution, min_dist, radius


def main():
    parser = argparse.ArgumentParser(description='Solve freeze-tag with the base CP-SAT model.')
    parser.add_argument('seconds', type=int, help='CP-SAT search time limit')
    parser.add_argument('type', help='Instance type, e.g. dig, tsplib_2d, or tree-2-2')
    parser.add_argument('--extended-cuts', action='store_true', help='Enable two-cycle and leaf cuts')
    args = parser.parse_args()
    space = get_instance(args.type)
    if space is None:
        parser.error(f'unknown instance type: {args.type}')
    if space.n < 1:
        parser.error('the instance must contain a root')

    source = 0
    sol_edges, upper = greedy_solution(source, space)
    source_radius = radius(source, space)
    lower = max(source_radius, min_dist(space) * ceil(log2(space.n))) if space.n > 1 else 0

    print('Freeze-Tag instance information:')
    print(f'\tTime limit = {args.seconds}s')
    print('\tNumber of nodes =', space.n)
    print('\tSource =', source)
    print(f'\tSource radius = {DELTA * source_radius:.2f}')
    print(f'\tGreedy bound = {DELTA * upper:.2f}')
    print(f'\tLower bound = {DELTA * lower:.2f}')
    print(f'\tInitial gap = {gap_percentage(upper, lower):.2f}%')
    if args.extended_cuts:
        print('\tExtended cuts = enabled')

    bdmhst = solve_ftp(
        list(space), space, args.seconds, lower, upper, True, sol_edges,
        extended_cuts=args.extended_cuts,
    )
    if bdmhst.status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        solution = get_solution(space, bdmhst)
        plot_solution('dot', space, solution)
    else:
        print('No solution found.')


if __name__ == '__main__':
    main()
