from ortools.sat.python import cp_model

from src.cp.utils import l2_norm, trivial_ub


def solve_bdhst(names, dist, degrees, max_time, ub, hop_depth=0, log=False, name='BDHST Problem', gap=0.0, init_sol=None):
    n = len(names)
    source = names[-1]
    names_to_i = {name: i for i, name in enumerate(names)}

    # Model creation:
    model = cp_model.CpModel()

    # Creates the variables:
    d_v = [model.NewIntVar(0, ub, f'd_{names[i]}') for i in range(n - 1)] + \
          [model.NewIntVar(0, 0, f'd_{source}')]
    x_e = [[model.NewBoolVar(f'x_({names[u]},{names[v]})') if v != u else None for v in range(n - 1)]
           for u in range(n)]
    depth = model.NewIntVar(0, ub, 'depth')

    # Creates the constraints:
    out_degree_sum = sum(x_e[-1])
    model.Add(out_degree_sum <= degrees[-1])
    model.Add(out_degree_sum >= 1)  # the out-degree of the root is at least one
    # the in-degree of the root is zero (there is no in-arc variable for the source)

    model.AddMaxEquality(depth, d_v)  # the depth of the tree is the maximum of the depth of each of its nodes

    # The in-degree is one and the out-degree is at most degree[v] for each internal node v:
    for v in range(n - 1):
        out_degree_sum = sum(uv for uv in x_e[v] if uv is not None)
        model.Add(out_degree_sum <= degrees[v])  # the out-degree is at most degree[v] for each internal node v
        model.Add(sum(x_e[j][v] for j in range(n) if j != v) == 1)  # the in-degree of internal nodes is one

    # A node depth is its parents depth plus the edge to it:
    for v in range(n - 1):
        for u in range(n):
            if u == v:  continue
            model.Add(d_v[v] == d_v[u] + dist(u, v)).OnlyEnforceIf(x_e[u][v])
            # Additional lb restriction to help the solver:
            model.Add(d_v[v] >= dist(-1, u) + dist(u, v)).OnlyEnforceIf(x_e[u][v])

    if hop_depth > 0:  # should control the hop depth
        hop_d_v = [model.NewIntVar(0, hop_depth, f'hop_d_{names[i]}') for i in range(n - 1)] + \
                  [model.NewIntVar(0, 0, f'hop_d_{source}')]

        # A node hop depth is its parents hop depth plus one:
        for v in range(n - 1):
            for u in range(n):
                if u == v:  continue
                model.Add(hop_d_v[v] == hop_d_v[u] + 1).OnlyEnforceIf(x_e[u][v])

    # Provides a warm start to the solver:
    if init_sol is not None:
        for u, v in init_sol:
            model.AddHint(x_e[names_to_i[u]][names_to_i[v]], 1)

    # Creates a solver and solves the model:
    model.Minimize(depth)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 8
    solver.parameters.max_time_in_seconds = max_time
    solver.parameters.log_search_progress = log
    solver.parameters.name = name
    solver.parameters.relative_gap_limit = gap
    status = solver.Solve(model)

    return status, model, solver, depth, d_v, x_e


def solve_ftp(names, dist, max_time, ub, log=False, init_sol=None):
    degrees = (len(names) - 1) * [2] + [1]
    return solve_bdhst(names, dist, degrees, max_time, ub, 0, log, 'Freeze-Tag Problem', 0, init_sol)


def solve_ftp_inner(sol_edges, d_tree, names_to_i, source, coords, grid_map, delta, max_time):
    cell = grid_map[source]  # current cell
    source_father = next(d_tree.predecessors(source), None)  # the father of this cell's representative

    # Gather the representatives for the children cells:
    leaves = []
    for child_cell_representative in d_tree.neighbors(source):
        child_cell = grid_map[child_cell_representative]
        rep, max_time = solve_ftp_inner(sol_edges, d_tree, names_to_i, child_cell[0], coords, grid_map, delta, max_time)
        leaves.append(rep)  # save the new dynamically chosen representative

    cell_names = leaves + cell[::-1]  # the source will be the last node
    degrees = len(leaves) * [0] + (len(cell) - 1) * [2]

    # Solve the sub-problem with source_father as source to let the solver figure out the best cell representative:
    if source_father is not None:  # this is not the root cell
        cell_names.append(source_father)
        degrees.append(2)
    degrees.append(1)
    n = len(cell_names)

    p = 100.0 * len(sol_edges) / (len(coords) - 1)
    print(f'Solving inner cell with {n} points - {p:.2f}% done with {max_time:.1f}s remaining...', end='\r')

    # Solve the sub-problem:
    cell_coords = [coords[names_to_i[v]] for v in cell_names]
    dist = l2_norm(cell_coords, delta)
    UB = trivial_ub(n, dist)
    status, _, solver, _, _, x_e = solve_bdhst(cell_names, dist, degrees, max_time, UB, 0, False, 'Freeze-Tag Problem')
    status = status == cp_model.FEASIBLE or status == cp_model.OPTIMAL
    max_time -= solver.WallTime()
    if not status:
        exit(f'\nCould not find any solution to the inner cell!')

    # Gather the sub-problem solution edges:
    new_source = None
    for u in range(n):
        for v in range(n - 1):
            if u == v: continue
            if solver.Value(x_e[u][v]):
                if cell_names[u] == source_father:
                    # Do not connect the new source to source_father, as the connection between source_father's cell and
                    # the current one will be defined later:
                    new_source = cell_names[v]  # this cell's new source selected by the solver
                else:
                    sol_edges.append((cell_names[u], cell_names[v]))

    return new_source, max_time
