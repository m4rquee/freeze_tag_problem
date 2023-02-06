import sys
from os.path import join
from heapq import heappop, heappush
from dataclasses import dataclass, field


def partitions(n, k, min_part=1):
    if k < 1: return
    if k == 1:
        if n >= min_part: yield n,
        return
    for i in range(min_part, n + 1):
        for result in partitions(n - i, k - 1, i):
            yield result + (i,)


def all_partitions(n, min_k, max_k, min_part=1):
    for k in range(min_k, max_k + 1):
        yield from partitions(n, k, min_part)


def print_tree(root, base_path, deg, n, count):
    with open(join(base_path, f'tree-deg{deg}-n{n}-{count}.g'), 'w+') as out_file:
        print('nnodes nedges type', file=out_file)
        print(N, N - 1, 'graph', file=out_file)
        print('nodename posx posy', file=out_file)
        root.print_nodes(out_file)
        print('endpoint1 endpoint2 weight', file=out_file)
        root.print_edges(out_file)


@dataclass(order=True)
class Node:
    sort_index: int = field(init=False, repr=False)

    code: int
    x: float
    y: float
    max_degree: int
    size: int
    depth: int
    min_degree: int = 1
    children_gen = iter(())
    children = []

    def __post_init__(self):
        self.sort_index = -self.code

    def __str__(self):
        if len(self.children) > 0:
            return f"{self.code}->[{', '.join(map(str, self.children))}]"
        return str(self.code)

    def print_nodes(self, out_file):
        print(self.code, self.x, self.y, file=out_file)
        for child in self.children: child.print_nodes(out_file)

    def print_edges(self, out_file):
        for child in self.children:
            print(self.code, child.code, 1, file=out_file)
            child.print_edges(out_file)

    def init(self):
        if self.size > 1:
            # Iterate over all ordered partitions of the descendents into up to d subtrees:
            self.children_gen = all_partitions(self.size - 1, self.min_degree, self.max_degree)
        return self

    def gen_children(self):  # return the next subtree rooted at this node
        if self.size == 1: return []  # nothing to do

        bags = next(self.children_gen, None)
        if bags is None: return []  # generated all subtrees

        base_code = self.code + 1
        self.children = []
        for bag in bags:
            x, y = base_code, self.depth + 1
            self.children.append(Node(base_code, x, y, self.max_degree, bag, self.depth + 1).init())
            base_code += bag
        return self.children


def gen_trees(n, deg, min_root_degree, base_path):
    root = Node(0, 0, 0, deg, n, 0, min_root_degree).init()
    node_heap = [root]
    count = 0  # count the number of generated trees

    while True:
        # Backtrack until an expandable node is found:
        removed_cache = []
        curr_node = heappop(node_heap)
        while not curr_node.gen_children():
            removed_cache.append(curr_node)
            if len(node_heap) == 0: return  # generated all trees
            curr_node = heappop(node_heap)

        # Add all nodes back:
        # Keep only the non descendents of curr_node and reset then:
        removed_cache = [r.init() for r in removed_cache if r.code >= curr_node.code + curr_node.size]
        heappush(node_heap, curr_node)
        leafs = curr_node.children + removed_cache
        for leaf in leafs: heappush(node_heap, leaf)

        # Grow a new tree from curr_node's id onwards:
        while len(node_heap) < n:
            next_leafs = []
            for leaf in leafs:
                for child in leaf.gen_children():  # expand the nodes children list into the heap
                    heappush(node_heap, child)
                    next_leafs.append(child)
            leafs = next_leafs
        print(root)  # has a full-grown n node tree
        count += 1
        print_tree(root, base_path, deg, n, count)


N = max(3, int(sys.argv[1]))  # must be at least 3
Degree = int(sys.argv[2])
Base_path = sys.argv[3]
gen_trees(N, Degree, 2, Base_path)
