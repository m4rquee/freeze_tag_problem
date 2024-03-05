from ortools.sat.python import cp_model

from src.cp.utils import l2_norm, trivial_ub


def solve_bdhst(names, dist, degrees, max_time, ub, hop_depth=0, log=False, name='BDHST Problem', gap=0.0):
    n = len(names)
    source = names[-1]

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
        model.Add(d_v[v] >= dist(-1, v))  # additional lb restriction to help the solver

    if hop_depth > 0:  # should control the hop depth
        hop_d_v = [model.NewIntVar(0, hop_depth, f'hop_d_{names[i]}') for i in range(n - 1)] + \
                  [model.NewIntVar(0, 0, f'hop_d_{source}')]

        # A node hop depth is its parents hop depth plus one:
        for v in range(n - 1):
            for u in range(n):
                if u == v:  continue
                model.Add(hop_d_v[v] == hop_d_v[u] + 1).OnlyEnforceIf(x_e[u][v])

    # Creates a solver and solves the model:
    model.Minimize(depth)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 4
    solver.parameters.max_time_in_seconds = max_time
    solver.parameters.log_search_progress = log
    solver.parameters.name = name
    solver.parameters.relative_gap_limit = gap
    status = solver.Solve(model)

    return status, model, solver, depth, d_v, x_e


def solve_ftp(names, dist, max_time, ub, log=False):
    degrees = (len(names) - 1) * [2] + [1]
    return solve_bdhst(names, dist, degrees, max_time, ub, 0, log, 'Freeze-Tag Problem')


def solve_ftp_inner(dg, names, names_to_i, source, coords, grid, delta, max_time):
    print()
    sol_edges = []
    for i, line in enumerate(grid):
        for j, cell in enumerate(line):
            if len(cell) == 0: continue

            root_name = names[cell[0]]
            if len(cell) == 1:
                original_edges = [(root_name, v) for v in dg.neighbors(root_name)]
                sol_edges.extend(original_edges)
                continue

            leaves = list(dg.neighbors(root_name))
            cell_names = leaves + [names[k] for k in cell[::-1]]
            degrees = len(leaves) * [0] + (len(cell) - 1) * [2]
            degrees.append(1 if root_name == source else 2)

            n = len(cell_names)
            cell_nodes_coords = [coords[names_to_i[v]] for v in cell_names]
            dist = l2_norm(cell_nodes_coords, delta)
            UB = trivial_ub(n, dist)

            print(f'Solving inner level cell ({i}; {j}) with {n} points...', end='\r', flush=True)
            status, _, solver, _, _, x_e = \
                solve_bdhst(cell_names, dist, degrees, max_time, UB, 0, False,
                            'Freeze-Tag Problem', 0.0)
            status = status == cp_model.FEASIBLE or status == cp_model.OPTIMAL
            if not status:
                exit(f'Could not find any solution to the inner level cell ({i}; {j})!')

            max_time -= solver.WallTime()
            for u in range(n):
                for v in range(n - 1):
                    if u == v: continue
                    if solver.Value(x_e[u][v]):
                        sol_edges.append((cell_names[u], cell_names[v]))

    print('\n')
    return sol_edges
