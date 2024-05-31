from math import log
import networkx as nx


def read_tsplib_2d_graph():
    while input() != 'NODE_COORD_SECTION': pass

    coords = []
    while (line := input()) != 'EOF':
        _, x, y = line.split()
        coord = (float(x), float(y))
        # if coord in coords: continue  # ignore duplicates
        coords.append(coord)
    n = len(coords)
    return list(range(n)), coords


def read_dig_graph():
    input()  # skip header
    n, m = map(int, input().split()[:-1])
    input()  # skip header

    while input() != 'tail head weight': pass

    edges = []
    for _ in range(m):
        line = input()
        tail, head, weight = line.split()
        edges.append((int(tail), int(head), {'weight': float(weight)}))
    return list(range(n)), edges


def read_tsplib_hcp_graph():
    while 'DIMENSION' not in (line := input()): pass

    n = int(line.split()[-1])

    while input() != 'EDGE_DATA_SECTION': pass

    edges = []
    while (line := input()) != 'EOF':
        tail, head = line.split()
        edges.append((int(tail) - 1, int(head) - 1, {'weight': 1}))
    return list(range(n)), edges


def gnp_graph(n, p, max_iter=1000):
    n = int(n)
    p = float(p)  # 2 * log(n) / n
    graph = None
    for _ in range(max_iter):
        aux = nx.fast_gnp_random_graph(n, p)
        if nx.is_connected(aux):
            graph = aux
            break
    if graph is None: exit(f'Could not generate a connected graph after {max_iter} trials!')
    return list(graph.nodes), list(graph.edges)


def ws_graph(n, k, p):
    n = int(n)
    p = float(p)
    k = int(k)
    aux = nx.connected_watts_strogatz_graph(n, k, p)
    return list(aux.nodes), list(aux.edges)


def tree_graph(r, h):
    r = int(r)
    h = int(h)
    aux = nx.balanced_tree(r, h)
    return list(aux.nodes), list(aux.edges)


def regular_graph(d, n):
    d = int(d)
    n = int(n)
    aux = nx.random_regular_graph(d, n)
    return list(aux.nodes), list(aux.edges)
