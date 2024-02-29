from ortools.sat.python import cp_model

inf = int(1E10)


def solve_ftp(names, source, dist, max_time):
    n = len(names)

    # Model creation:
    model = cp_model.CpModel()

    # Creates the variables:
    d_v = [model.NewIntVar(0, inf, f'd_{names[i]}') for i in range(n - 1)] + \
          [model.NewIntVar(0, 0, f'd_{names[source]}')]
    x_e = [[model.NewBoolVar(f'x_({names[u]},{names[v]})') if v != u else None for v in range(n - 1)] for u in range(n)]
    depth = model.NewIntVar(0, inf, 'depth')

    # Creates the constraints:
    model.AddMaxEquality(depth, d_v)  # the depth of the tree is the maximum of the depth of each of its nodes

    model.Add(sum(x_e[source]) == 1)  # the out-degree of the root is one
    # the in-degree of the root is zero (there is no in-arc variable for the source)

    # The in-degree is one and the out-degree is at most two for each internal node:
    for v in range(n - 1):
        model.Add(sum(uv for uv in x_e[v] if uv is not None) <= 2)  # the out-degree of internal nodes is at most two
        model.Add(sum(x_e[j][v] for j in range(n) if j != v) == 1)  # the in-degree of internal nodes is one

    # A node depth is its parents depth plus the edge to it:
    for v in range(n - 1):
        for u in range(n):
            if u == v:  continue
            model.Add(d_v[v] == d_v[u] + dist(u, v)).OnlyEnforceIf(x_e[u][v])
        model.Add(d_v[v] >= dist(source, v))

    # Creates a solver and solves the model:
    model.Minimize(depth)
    print('Solving the model...')
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 4
    solver.parameters.max_time_in_seconds = max_time

    return solver.Solve(model), solver, depth, d_v, x_e


def solve_bdhst(names, source_i, dist, degrees, max_time):
    n = len(names)

    # Model creation:
    model = cp_model.CpModel()

    # Creates the variables:
    d_v = [model.NewIntVar(0, inf, f'd_{names[i]}') for i in range(n - 1)] + \
          [model.NewIntVar(0, 0, f'd_{names[source_i]}')]
    x_e = [[model.NewBoolVar(f'x_({names[u]},{names[v]})') if v != u else None for v in range(n - 1)] for u in range(n)]
    depth = model.NewIntVar(0, inf, 'depth')

    # Creates the constraints:
    model.AddMaxEquality(depth, d_v)  # the depth of the tree is the maximum of the depth of each of its nodes

    out_degree_sum = sum(x_e[source_i])
    model.Add(out_degree_sum <= degrees[source_i])
    model.Add(out_degree_sum >= 1)  # the out-degree of the root is at least one
    # the in-degree of the root is zero (there is no in-arc variable for the source)

    # The in-degree is one and the out-degree is at most two for each internal node:
    for v in range(n - 1):
        out_degree_sum = sum(uv for uv in x_e[v] if uv is not None)
        model.Add(out_degree_sum <= degrees[v])
        model.Add(sum(x_e[j][v] for j in range(n) if j != v) == 1)  # the in-degree of internal nodes is one

    # A node depth is its parents depth plus the edge to it:
    for v in range(n - 1):
        for u in range(n):
            if u == v:  continue
            model.Add(d_v[v] == d_v[u] + dist(u, v)).OnlyEnforceIf(x_e[u][v])
        model.Add(d_v[v] >= dist(source_i, v))

    # The total number of edges should be n-1:
    edge_sum = 0
    for n_v in x_e:
        for e in n_v:
            if e is not None:
                edge_sum += e
    model.Add(edge_sum == n - 1)

    # Creates a solver and solves the model:
    model.Minimize(depth)
    print('Solving the model...')
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 4
    solver.parameters.max_time_in_seconds = max_time

    return solver.Solve(model), solver, depth, d_v, x_e
