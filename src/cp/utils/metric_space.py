from functools import cache
from math import sqrt, floor

import networkx as nx
from networkx.algorithms.approximation.steinertree import metric_closure

from src.cp.utils.constants import DELTA


class MetricSpace:
    def __init__(self, n, func=None, delta=DELTA):
        self.n = n
        self.func = func
        self.delta = delta

    def __call__(self, u, v):
        return floor(self.func(u, v) / self.delta)  # scale up so we can treat rational values as integers

    def __len__(self):
        return self.n

    def __getitem__(self, points):  # restrict to a sub-set of points
        def aux(u, v):
            u, v = points[u], points[v]
            return self.func(u, v)

        return MetricSpace(len(points), aux, self.delta)

    def __iter__(self):
        yield from range(self.n)


class L2Norm(MetricSpace):
    def __init__(self, coords, delta=DELTA):
        @cache
        def aux(u, v):
            d = sqrt((self.coords[u][0] - self.coords[v][0]) ** 2 + (self.coords[u][1] - self.coords[v][1]) ** 2)
            return d

        super().__init__(len(coords), aux, delta)
        self.coords = coords


class GraphDist(MetricSpace):
    def __init__(self, delta=DELTA):
        super().__init__(None, None, delta)

    def __init_base__(self, n):
        self.n = n
        @cache
        def aux(u, v):
            if u == v: return 0
            d = self.graph.get_edge_data(u, v)['distance']
            return d
        self.func = aux

    def init_from_edges(self, n, edges):
        self.__init_base__(n)
        self.original_graph = nx.empty_graph(n)  # create an empty graph first to preserve node order
        self.original_graph.add_edges_from(edges)
        self.graph = metric_closure(self.original_graph)

    def init_from_graph(self, graph):
        self.init_from_edges(len(graph), graph.edges(data=True))
