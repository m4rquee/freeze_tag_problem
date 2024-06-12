from math import ceil

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

from src.cp.utils.metric_space import GraphDist

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


def plot_solution(layout, space, solution, upper_solution=None, eps=None):
    plt.figure(figsize=(10, 6))
    if isinstance(space, GraphDist):
        match layout:
            case 'spring': coords_dict = nx.spring_layout(space.original_graph)
            case 'spectral': coords_dict = nx.spectral_layout(space.original_graph)
            case _: coords_dict = nx.nx_agraph.graphviz_layout(space.original_graph, prog=layout)

        plot_graph(space.original_graph, coords_dict, 'white', 'gray', style='dashed', node_size=40)
    else:
        coords_dict = dict(zip(space, space.coords))

    plot_graph(solution.tree, coords_dict, solution.node_colors, 'green', style='solid', node_size=40,
                      connectionstyle='arc3,rad=0.15')
    if upper_solution is not None:
        plot_graph(upper_solution.tree, coords_dict, upper_solution.node_colors, 'red', style='dotted',
                   node_size=10, connectionstyle='arc3,rad=0.1')

    if eps is not None: plot_grid(eps)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
