import sys
from os.path import join
from random import shuffle, random


def print_nodes(n, out_file):
    for i in range(n): print(i, n, n * n, file=out_file)


def print_edges(edges, out_file):
    for u, v in edges: print(u, v, 1, file=out_file)


def print_graph(n, p, m, edges, base_path):
    with open(join(base_path, f'gnp-graph-n{n}-p{p}.dig'), 'w+') as out_file:
        print('nnodes narcs type', file=out_file)
        print(n, m, 'digraph', file=out_file)
        print('nodename posx posy', file=out_file)
        print_nodes(n, out_file)
        print('tail head weight', file=out_file)
        print_edges(edges, out_file)


def gen_graph(n, p, base_path):
    nodes = list(range(n))
    shuffle(nodes)
    edges = set()
    for i in range(n - 1):  # make sure the graph is connected
        edges.add(tuple(sorted((nodes[i], nodes[i + 1]))))
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if random() >= 1 - p:
                edges.add(tuple(sorted((i, j))))
    print_graph(n, p, len(edges), edges, base_path)


N = int(sys.argv[1])  # number of nodes
P = float(sys.argv[2])  # probability that a given edge is present
Base_path = sys.argv[3]
gen_graph(N, P, Base_path)
