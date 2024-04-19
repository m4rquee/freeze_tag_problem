def read_tsplib_graph():
    while input() != "NODE_COORD_SECTION":
        pass

    names = []
    coords = []
    while (line := input()) != "EOF":
        name, x, y = line.split()
        coord = (float(x), float(y))
        # if coord in coords: continue  # ignore duplicates
        names.append(int(name))
        coords.append(coord)
    # todo: change the other functions instead of inverting the output
    return names[::-1], coords[::-1]


def read_dig_graph():
    input()  # skip header
    _, m = map(int, input().split()[:-1])
    input()  # skip header

    names = []
    coords = []
    while (line := input()) != "tail head weight":
        name, x, y = line.split()
        coord = (float(x), float(y))
        names.append(int(name))
        coords.append(coord)

    edges = []
    for _ in range(m):
        line = input()
        tail, head, weight = line.split()
        edges.append((int(tail), int(head), {'weight': float(weight)}))
    # todo: change the other functions instead of inverting the output
    return names[::-1], coords[::-1], edges