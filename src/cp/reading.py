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
    return names, coords
