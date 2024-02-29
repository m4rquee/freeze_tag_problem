from ortools.sat.python import cp_model


def solve_bdhst(names, dist, degrees, max_time, ub, name='BDHST Problem'):
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

    # Creates a solver and solves the model:
    model.Minimize(depth)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 4
    solver.parameters.max_time_in_seconds = max_time
    solver.parameters.log_search_progress = True
    solver.parameters.name = name
    status = solver.Solve(model)

    return status, model, solver, depth, d_v, x_e


def solve_ftp(names, dist, max_time, ub):
    degrees = (len(names) - 1) * [2] + [1]
    return solve_bdhst(names, dist, degrees, max_time, ub, 'Freeze-Tag Problem')
