from ortools.sat.python import cp_model

from src.cp.solving.ftp_cuts import add_ftp_extended_cuts
from src.cp.utils.utils import *


def solve_bdmhst(names, dist, degrees, max_time, lb=0, ub=float('inf'), hop_depth=0, log=False, name='BDMHST Problem',
                gap=0.0, init_sol=None, *, extended_cuts=False, workers=8):
    n = len(names)
    source = 0
    if extended_cuts and list(degrees) != [1] + [2] * (n - 1):
        raise ValueError('FTP extended cuts require root capacity one and nonroot capacity two')

    # Model creation:
    model = cp_model.CpModel()

    # Creates the variables:
    d_v = [model.NewIntVar(0, 0, f'd_{source}')] + \
          [model.NewIntVar(0, ub, f'd_{names[i]}') for i in range(1, n)]

    x_e = [[0] + [model.NewBoolVar(f'x_({names[u]},{names[v]})') if v != u else 0 for v in range(1, n)]
           for u in range(n)]
    depth = model.NewIntVar(lb, ub, 'depth')

    # Creates the constraints:
    out_degree_sum = sum(x_e[source])
    if n > 1:
        model.Add(out_degree_sum <= degrees[source])
        model.Add(out_degree_sum >= 1)  # the out-degree of the root is at least one
    # the in-degree of the root is zero (there is no in-arc variable for the source)

    model.AddMaxEquality(depth, d_v)  # the depth of the tree is the maximum of the depth of each of its nodes

    # The in-degree is one and the out-degree is at most degree[v] for each internal node v:
    for v in range(1, n):
        out_degree_sum = sum(x_e[v])
        model.Add(out_degree_sum <= degrees[v])  # the out-degree is at most degree[v] for each internal node v
        model.Add(sum(x_e[u][v] for u in range(n)) == 1)  # the in-degree of internal nodes is one

    # A node depth is its parents depth plus the edge to it:
    for v in range(1, n):
        for u in range(n):
            if u == v:  continue
            model.Add(d_v[v] == d_v[u] + dist(u, v)).OnlyEnforceIf(x_e[u][v])
            # Additional lb restriction to help the solver:
            model.Add(d_v[v] >= dist(source, u) + dist(u, v)).OnlyEnforceIf(x_e[u][v])

    if hop_depth > 0:  # should control the hop depth
        hop_d_v = [model.NewIntVar(0, 0, f'hop_d_{source}')] + \
                  [model.NewIntVar(0, hop_depth, f'hop_d_{names[i]}') for i in range(1, n)]

        # A node hop depth is its parents hop depth plus one:
        for v in range(1, n):
            for u in range(n):
                if u == v:  continue
                model.Add(hop_d_v[v] == hop_d_v[u] + 1).OnlyEnforceIf(x_e[u][v])

    if extended_cuts:
        arcs = {(u, v): x_e[u][v] for u in range(n) for v in range(1, n) if u != v}
        add_ftp_extended_cuts(model, n, arcs)

    # Provides a warm start to the solver:
    if init_sol is not None:
        for u, v in init_sol:
            model.AddHint(x_e[u][v], 1)

    # Creates a solver and solves the model:
    model.Minimize(depth)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = max_time
    solver.parameters.log_search_progress = log
    solver.parameters.name = name
    solver.parameters.relative_gap_limit = gap
    status = solver.Solve(model)

    return BDMHST(source, ub, status, model, solver, depth, d_v, x_e)


def solve_ftp(names, dist, max_time, lb, ub, log=False, init_sol=None, *, extended_cuts=False):
    degrees = [1] + (len(names) - 1) * [2]
    return solve_bdmhst(names, dist, degrees, max_time, lb, ub, log=log, name='Freeze-Tag Problem',
                       init_sol=init_sol, extended_cuts=extended_cuts, workers=16)


def solve_ftp_inner(sol_edges, d_tree, source, dist, cluster_map, max_time):
    cluster = cluster_map[source]  # current cluster
    source_father = next(d_tree.predecessors(source), None)  # the father of this cluster's representative

    # Gather the representatives for the children cells:
    leaves = []
    for child_cell_representative in d_tree.neighbors(source):
        child_cell = cluster_map[child_cell_representative]
        rep, max_time = solve_ftp_inner(sol_edges, d_tree, child_cell[0], dist, cluster_map, max_time)
        leaves.append(rep)  # save the new dynamically chosen representative

    cell_names = cluster + leaves
    degrees = (len(cluster) - 1) * [2] + len(leaves) * [0]

    # Solve the sub-problem with source_father as source to let the solver figure out the best cluster representative:
    if source_father is not None:  # this is not the root cluster
        cell_names.insert(0, source_father)
        degrees.insert(0, 2)
    degrees.insert(0, 1)
    n = len(cell_names)

    p = 100.0 * len(sol_edges) / (len(dist) - 1)
    print(f'Solving inner cluster with {n} points - {p:.2f}% done with {max_time:.1f}s remaining...', end='\r')

    # Solve the sub-problem:
    local_dist = dist[cell_names]

    source_radius = radius(0, local_dist)
    min_edge = min_dist(local_dist)
    LB = min_edge * ceil(log2(n))
    cell_sol_edges, UB = greedy_solution(0, local_dist)  # it is not yet valid because there are fixed leaves
    UB += 2 * source_radius  # account for the leaves by adding a relocation cost

    bdmhst = solve_bdmhst(cell_names, local_dist, degrees, max_time, LB, UB, 0, False, 'Freeze-Tag Problem', 0,
        cell_sol_edges)
    status = bdmhst.status == cp_model.FEASIBLE or bdmhst.status == cp_model.OPTIMAL
    max_time -= bdmhst.solver.WallTime()
    if not status:
        exit(f'\nCould not find any solution to the inner cluster!')

    # Gather the sub-problem solution edges:
    new_source = None
    for u in range(n):
        for v in range(1, n):
            if u == v: continue
            if bdmhst.solver.Value(bdmhst.x_e[u][v]):
                if cell_names[u] == source_father:
                    # Do not connect the new source to source_father, as the connection between source_father's cluster
                    # and the current one will be defined later:
                    new_source = cell_names[v]  # this cluster's new source selected by the solver
                else:
                    sol_edges.append((cell_names[u], cell_names[v]))

    return new_source, max_time


class BDMHST():
    def __init__(self, source, UB, status, model, solver, depth, d_v, x_e):
        self.source = source
        self.UB = UB
        self.status = status
        self.model = model
        self.solver = solver
        self.depth = depth
        self.d_v = d_v
        self.x_e = x_e


class Solution():
    def __init__(self, tree, node_colors):
        self.tree = tree
        self.node_colors = node_colors


def gap_percentage(upper, lower):
    if lower == 0:
        return 0.0 if upper == 0 else float('inf')
    return 100 * (upper - lower) / lower


def get_solution(space, bdmhst):
    print('Full FTP solution:')
    depth = bdmhst.solver.Value(bdmhst.depth)
    print(f'  solution makespan: {DELTA * depth:.2f}')
    lb = bdmhst.solver.BestObjectiveBound()
    print(f'  gap: {gap_percentage(depth, lb):.2f}%')
    print(f'  greedy gap: {gap_percentage(bdmhst.UB, lb):.2f}%')
    print(f'  d_v: (', end='')
    for v in range(space.n - 1):
        depth_v = bdmhst.solver.Value(bdmhst.d_v[v])
        print(f'{DELTA * depth_v:.2f}', end=', ')
    print(f'{DELTA * bdmhst.solver.Value(bdmhst.d_v[-1])})')

    print('  edges: ', end='')
    sol_edges = []
    for u in range(space.n):
        for v in range(1, space.n):
            if u == v: continue
            if bdmhst.solver.Value(bdmhst.x_e[u][v]):
                print(f'{u}-{v}; ', end='')
                sol_edges.append((u, v))
    tree = nx.DiGraph(sol_edges)
    tree.add_node(bdmhst.source)
    hop_depth = calc_depth(bdmhst.source, tree)
    print(f'\n  hop depth: {hop_depth}')

    node_colors = []
    for v in tree.nodes:
        if v == bdmhst.source:
            node_colors.append('red')
        else:
            node_colors.append('cyan' if bdmhst.solver.Value(bdmhst.d_v[v]) == depth else 'black')

    return Solution(tree, node_colors)
