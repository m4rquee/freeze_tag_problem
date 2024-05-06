import networkx as nx


def read_tsplib_2d_graph():
    while input() != 'NODE_COORD_SECTION': pass

    names = []
    coords = []
    while (line := input()) != 'EOF':
        name, x, y = line.split()
        coord = (float(x), float(y))
        # if coord in coords: continue  # ignore duplicates
        names.append(int(name))
        coords.append(coord)
    # todo: change the other functions instead of inverting the output
    return names[::-1], coords[::-1]


def read_dig_graph():
    input()  # skip header
    _, m = map(int, input().split()[:-1])
    input()  # skip header

    names = []
    while (line := input()) != 'tail head weight':
        name, x, y = line.split()
        names.append(int(name))

    edges = []
    for _ in range(m):
        line = input()
        tail, head, weight = line.split()
        edges.append((int(tail), int(head), {'weight': float(weight)}))
    # todo: change the other functions instead of inverting the output
    return names[::-1], edges


def read_tsplib_hcp_graph():
    while 'DIMENSION' not in (line := input()): pass

    n = int(line.split()[-1])
    names = list(range(1, n + 1))

    while input() != 'EDGE_DATA_SECTION': pass

    edges = []
    while (line := input()) != 'EOF':
        tail, head = line.split()
        edges.append((int(tail), int(head), {'weight': 1}))
    # todo: change the other functions instead of inverting the output
    return names[::-1], edges


def gnp_graph(n, p):
    n = int(n)
    p = float(p)
    aux = nx.erdos_renyi_graph(n, p)
    # todo: change the other functions instead of inverting the output
    return list(aux.nodes)[::-1], list(aux.edges)


def ws_graph(n, k, p):
    n = int(n)
    p = float(p)
    k = int(k)
    aux = nx.connected_watts_strogatz_graph(n, k, p)
    # todo: change the other functions instead of inverting the output
    return list(aux.nodes)[::-1], list(aux.edges)
