from math import ceil

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt


def plot_graph(graph, node_coords, node_color='black', edge_color='black', **kwds):
    nx.draw(graph, pos=node_coords, node_color=node_color, edge_color=edge_color, **kwds)


def plot_grid(eps):
    num_cells = ceil(1 / eps)
    for i in range(num_cells + 1):
        plt.axvline(x=i, linestyle='--', color='gray', alpha=0.5)
        plt.axhline(y=i, linestyle='--', color='gray', alpha=0.5)


def plot_graph3d(graph, node_coords, node_color='black', edge_color='black', **kwds):
    cf = plt.gcf()
    cf.set_facecolor('w')
    if cf.axes:
        ax = cf.gca()
    else:
        ax = cf.add_subplot(111, projection='3d')

    # Plot the nodes:
    node_xyz = np.array([node_coords[v] for v in graph])
    ax.scatter(*node_xyz.T, c=node_color, **kwds)

    # Plot the edges:
    edge_xyz = np.array([(node_coords[u], node_coords[v]) for u, v in graph.edges()])
    for edge_coords in edge_xyz:
        ax.plot(*edge_coords.T, color=edge_color)
    
    ax.set_axis_off()
    plt.draw_if_interactive()
