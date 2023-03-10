import sys
from typing import List
from os.path import join
from random import randrange
from dataclasses import dataclass
from collections import defaultdict
from more_itertools import distinct_permutations as permutations


def ordered_partitions(n, k, max_part):
    partition = [0] * (k + 2)
    curr_k = 1
    y = n - 1
    while curr_k != 0:
        x = partition[curr_k - 1] + 1
        curr_k -= 1
        while 2 * x <= y and curr_k < k:
            partition[curr_k] = x
            y -= x
            curr_k += 1
        last = curr_k + 1
        while x <= y:
            partition[curr_k] = x
            partition[last] = y
            if curr_k == k - 2:
                ret = partition[:k]
                if ret[-1] <= max_part: yield from permutations(ret)
            x += 1
            y -= 1
        partition[curr_k] = x + y
        y = x + y - 1
        if curr_k == k - 1:
            ret = partition[:k]
            if ret[-1] <= max_part: yield from permutations(ret)


def print_tree(root, base_path, deg, n, count):
    with open(join(base_path, f'tree-deg{deg}-n{n}-{count}.dig'), 'w+') as out_file:
        print('nnodes narcs type', file=out_file)
        print(N, N - 1, 'digraph', file=out_file)
        print('nodename posx posy', file=out_file)
        root.print_nodes(out_file)
        print('tail head weight', file=out_file)
        root.print_edges(out_file)


@dataclass
class Node:
    code: int
    x: int
    y: int
    num_children: int
    children: List = None

    def __post_init__(self):
        self.children = []

    def __str__(self):
        if self.num_children > 0:
            return f"{self.code}->[{', '.join(map(str, self.children))}]"
        return str(self.code)

    def print_nodes(self, out_file):
        print(self.code, self.x, self.y, file=out_file)
        for child in self.children: child.print_nodes(out_file)

    def print_edges(self, out_file):
        for child in self.children:
            print(self.code, child.code, 1, file=out_file)
            child.print_edges(out_file)


def gen_tree_from_seq(d_seq, root_degree):
    code = 0
    queue = [root := Node(code, 0, 0, root_degree)]
    nodes_per_depth = defaultdict(int)
    while queue:
        first = queue.pop(0)
        if first.num_children > len(d_seq): return None  # this sequence is invalid
        curr_child_deg = float('inf')
        for _ in range(first.num_children):
            code += 1
            y = first.y + 1
            x = nodes_per_depth[y]
            nodes_per_depth[y] += 1
            new_child_deg = d_seq.pop(0) - 1
            if new_child_deg > curr_child_deg: return None  # consider only non-increasing sequences (skip isomorphisms)
            curr_child_deg = new_child_deg
            first.children.append(new_child := Node(code, x, y, new_child_deg))
            if new_child.num_children > 0: queue.append(new_child)  # it's not a leaf and so there is nothing to expand
    if len(d_seq) != 0: return None  # this sequence is invalid
    return root


def gen_trees(n, max_degree, root_degree, base_path):
    m = N - 1  # number of edges
    degree_sum = 2 * m - root_degree  # the sum of the degree sequence equals 2M
    count = 0
    printed = 0
    random_num = randrange(n)
    for d_seq in ordered_partitions(degree_sum, n - 1, max_degree):
        count += 1
        if count % n ** 2 == random_num:  # skip most of the generated sequences
            d_seq = list(reversed(d_seq))  # look first for tree that are more dense closer to the root
            print(f'\rdegree sequence count = {count}; generated trees = {printed}...', end='')
            if (tree := gen_tree_from_seq(d_seq, root_degree)) is not None:
                printed += 1
                print_tree(tree, base_path, max_degree, n, printed)
                yield tree
                if printed == 1: break
                random_num = randrange(n)
    print('\r', end='')
    yield count


N = int(sys.argv[1])  # number of tree nodes
Max_degree = int(sys.argv[2])
Root_degree = int(sys.argv[3])
Base_path = sys.argv[4]
print(f'Generating trees with {N} nodes, maximum degree {Max_degree} and root degree {Root_degree}:')
print(*gen_trees(N, Max_degree, Root_degree, Base_path), sep='\n')
print('trees were generated')
