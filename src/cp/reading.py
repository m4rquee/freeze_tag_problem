def read_tsplib_graph():
    while input() != "NODE_COORD_SECTION":
        pass

    names = []
    coords = []
    while (line := input()) != "EOF":
        name, x, y = line.split()
        names.append(int(name))
        coords.append((float(x), float(y)))
    return names, coords
