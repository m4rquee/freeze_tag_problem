import sys
from math import cos, sin, pi
from os.path import join
from random import shuffle, sample


def print_nodes(n, out_file):
    step = 2 * pi / n
    for i in range(n): print(i, cos(i * step), sin(i * step), file=out_file)


def print_edges(edges, out_file):
    for u, v in edges: print(u, v, 1, file=out_file)


def print_graph(n, m, edges, base_path):
    with open(join(base_path, f'unit-graph-n{n}-m{m}.dig'), 'w+') as out_file:
        print('nnodes narcs type', file=out_file)
        print(n, m, 'digraph', file=out_file)
        print('nodename posx posy', file=out_file)
        print_nodes(n, out_file)
        print('tail head weight', file=out_file)
        print_edges(edges, out_file)


def gen_graph(n, m, base_path):
    nodes = list(range(n))
    shuffle(nodes)
    edges = set()
    for i in range(n - 1):  # make sure the graph is connected
        edges.add((nodes[i], nodes[i + 1]))
    while len(edges) < m:  # add the remaining edges
        e = sample(nodes, 2)
        e.sort()
        edges.add(tuple(e))
    print_graph(n, m, edges, base_path)


N = int(sys.argv[1])  # number of nodes
M = int(sys.argv[2])  # number of edges
Base_path = sys.argv[3]
gen_graph(N, M, Base_path)
