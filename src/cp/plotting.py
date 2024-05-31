from math import ceil

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt


def plot_solution(dg, node_coords, node_colors, edge_color='black', **kwds):
    nx.draw(dg, pos=node_coords, node_color=node_colors, edge_color=edge_color, **kwds)


def plot_grid(eps):
    num_cells = ceil(1 / eps)
    for i in range(num_cells + 1):
        plt.axvline(x=i, linestyle='--', color='gray', alpha=0.5)
        plt.axhline(y=i, linestyle='--', color='gray', alpha=0.5)


def plot_solution3d(ax, dg, node_coords, node_colors, edge_color='black', **kwds):
    node_xyz = np.array([node_coords[v] for v in sorted(dg)])

    # Plot the nodes
    ax.scatter(*node_xyz.T, c=node_colors, **kwds)

    edge_xyz = np.array([(node_coords[u], node_coords[v]) for u, v in dg.edges()])
    # Plot the edges
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color=edge_color)
