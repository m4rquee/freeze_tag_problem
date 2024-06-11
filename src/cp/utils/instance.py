from math import log
import networkx as nx
from random import sample


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
    return list(graph.nodes), list(graph.edges(data=True))


def ws_graph(n, k, p):
    n = int(n)
    p = float(p)
    k = int(k)
    aux = nx.connected_watts_strogatz_graph(n, k, p)
    return list(aux.nodes), list(aux.edges(data=True))


def tree_graph(r, h):
    r = int(r)
    h = int(h)
    aux = nx.balanced_tree(r, h)
    return list(aux.nodes), list(aux.edges(data=True))


def regular_graph(d, n):
    d = int(d)
    n = int(n)
    aux = nx.random_regular_graph(d, n)
    return list(aux.nodes), list(aux.edges(data=True))


def lobster_graph(n, p1, p2):
    n = int(n)
    p1 = float(p1)
    p2 = float(p2)
    aux = nx.random_lobster(n, p1, p2)
    return list(aux.nodes), list(aux.edges(data=True))


def combine_graphs(g, h, k=1):
    prod = nx.empty_graph(len(g) * len(h))
    for e in g.edges:
        for i in sample(range(len(h)), k):
            prod.add_edge(e[0] * len(h) + i, e[1] * len(h) + i)
    for e in h.edges:
        for v in g:
            prod.add_edge(v * len(h) + e[0], v * len(h) + e[1])
    return list(prod.nodes), list(prod.edges(data=True))


def get_instance(instance_type):
    instance, *args = instance_type.split('-')
    match instance:
        case 'tsplib_hcp': return read_tsplib_hcp_graph()
        case 'gnp': return gnp_graph(*args)
        case 'ws': return ws_graph(*args)
        case 'tree': return tree_graph(*args)
        case 'regular': return regular_graph(*args)
        case 'lobster': return lobster_graph(*args)
        case _: return read_dig_graph()  # if it is 'dig'
