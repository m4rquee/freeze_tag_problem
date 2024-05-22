from math import ceil

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt


def plot_solution(dg, sol_edges, node_coords, node_names, node_colors, edge_color="black", **kwds):
    nx.draw(dg, pos=node_coords, node_color=node_colors, edge_color=edge_color, **kwds)


def plot_grid(eps):
    num_cells = ceil(1 / eps)
    for i in range(num_cells + 1):
        plt.axvline(x=i, linestyle='--', color='gray', alpha=0.5)
        plt.axhline(y=i, linestyle='--', color='gray', alpha=0.5)


def plot_solution3d(dg, node_coords, node_names, node_colors, edge_color='black', **kwds):
    node_xyz = np.array([node_coords[v] for v in sorted(dg)])

    # Create the 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the nodes - alpha is scaled by "depth" automatically
    ax.scatter(*node_xyz.T, s=100, ec="w", c=node_colors)

    edge_xyz = np.array([(node_coords[u], node_coords[v]) for u, v in dg.edges()])
    # Plot the edges
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color="tab:gray")
