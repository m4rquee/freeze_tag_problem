"""Optional cuts for FTP with one root child and binary branching elsewhere."""


def add_ftp_extended_cuts(model, n, arcs):
    """Add the two-cycle and leaf cuts, and return the leaf variables."""
    outgoing = [[] for _ in range(n)]
    for (parent, child), arc in arcs.items():
        outgoing[parent].append(arc)
        if parent < child and (child, parent) in arcs:
            model.add(arc + arcs[child, parent] <= 1).with_name(f'two_cycle_{parent}_{child}')

    leaves = [model.new_bool_var(f'leaf_{i}') for i in range(n)]
    for i, leaf in enumerate(leaves):
        children = sum(outgoing[i])
        capacity = 1 if i == 0 else 2
        model.add(children <= capacity * (1 - leaf)).with_name(f'leaf_degree_upper_{i}')
        model.add(children >= 1 - leaf).with_name(f'leaf_degree_lower_{i}')
    # For n > 1, k internal vertices can supply at most 1 + 2(k - 1) of the n - 1 arcs.
    model.add(sum(leaves) <= max(1, n // 2)).with_name('leaf_count')
    return leaves
