from math import ceil

import networkx as nx
from matplotlib import pyplot as plt


def plot_solution(dg, sol_edges, node_coords, node_names, node_colors, edge_color="black", **kwds):
    aux = dg.edge_subgraph(sol_edges).subgraph(node_names)
    # todo: fix networkx random node shuffling
    nx.draw(aux, pos=node_coords, node_color=node_colors, edge_color=edge_color, **kwds)


def plot_grid(eps):
    num_cells = ceil(1 / eps)
    for i in range(num_cells + 1):
        plt.axvline(x=i, linestyle='--', color='gray', alpha=0.5)
        plt.axhline(y=i, linestyle='--', color='gray', alpha=0.5)
