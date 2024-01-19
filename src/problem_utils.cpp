#include "problem_utils.hpp"

Problem_Instance::Problem_Instance(const string &filename, int time_limit, bool calc_clojure, bool tsplib)
    : time_limit(time_limit), vname(g), weight(g), px(g), py(g), original(g), solution(g), node_activation(g) {
    read_instance(filename, tsplib);
    source_radius = graph_radius = -1;
    if (calc_clojure) clojure();

    // Initialize the attributes:
    nnodes = countNodes(g);
    narcs = countArcs(g);
    source = INVALID;
#ifndef BDHST
    if (tsplib) source = Digraph::nodeFromId(g.maxNodeId());
    else
        source = Digraph::nodeFromId(0);
#endif
    for (ArcIt e(g); e != INVALID; ++e) {
        arc_map[g.source(e)][g.target(e)] = e;
        solution[e] = false;
    }
    for (DNodeIt v(g); v != INVALID; ++v) node_activation[v] = -1;
}

void Problem_Instance::start_counter() { start = chrono::system_clock::now(); }

void Problem_Instance::stop_counter() { stop = chrono::system_clock::now(); }

void Problem_Instance::print_instance() {
#ifdef BDHST
    cout << "Bounded Degree Minimum Height Spanning Tree instance information:" << endl;
#else
    cout << "Freeze-Tag instance information:" << endl;
#endif
    cout << "\tTime limit = " << time_limit << "s" << endl;
    cout << "\tNumber of nodes = " << nnodes << endl;
    cout << "\tNumber of edges = " << narcs << endl;
    if (source != INVALID) {
        cout << "\tSource = " << vname[source] << endl;
        cout << "\tSource radius = " << source_radius * MY_EPS << endl;
    }
    cout << "\tGraph radius = " << graph_radius * MY_EPS << endl;
    cout << endl;
}

inline int euclidean_distance(DNode &u, DNode &v, DNodePosMap &posx, DNodePosMap &posy) {
    auto dist = sqrt(pow(posx[u] - posx[v], 2) + pow(posy[u] - posy[v], 2));
    return floor(dist / MY_EPS);// scale up so we can treat rational values as integers
}

// todo: move this function to mygraphlib
void Problem_Instance::read_tsplib_instance(const string &filename) {
    ifstream file;
    file.open(filename.c_str());
    if (!file) {
        stringstream buffer;
        buffer << "Error: Could not open file " << filename << "." << endl;
        throw runtime_error(buffer.str());
    }

    string line;
    do {// skip the header lines
        getline(file, line);
    } while (line != "NODE_COORD_SECTION");

    int id;
    double x, y;
    for (getline(file, line); line != "EOF"; getline(file, line)) {
        sscanf(line.c_str(), "%d %lf %lf", &id, &x, &y);
        auto u = g.addNode();
        vname[u] = to_string(id - 1);
        px[u] = x;
        py[u] = y;
        DNodeIt v(g);
        for (++v; v != INVALID; ++v) weight[g.addArc(u, v)] = euclidean_distance(u, v, px, py);
    }
}

void Problem_Instance::read_instance(const string &filename, bool tsplib) {
    if (tsplib) {
        read_tsplib_instance(filename);
    } else {
        ReadDigraph(filename, g, vname, px, py, weight);
        for (ArcIt e(g); e != INVALID; ++e)
            weight[e] = floor(weight[e] / MY_EPS);// scale up so we can treat rational values as integers
    }

    // Ensure the digraph is strongly connected by adding each arc reverse:
    int i = 0;
    vector<Arc> original_edges(countArcs(g));
    for (ArcIt e(g); e != INVALID; ++e, i++) original[original_edges[i] = e] = true;// save the original arcs
    for (i--; i >= 0; i--) {
        auto e = original_edges[i];
        auto rev = findArc(g, g.target(e), g.source(e));
        if (rev == INVALID) {// if the reverse of arc e is not already present
            rev = g.addArc(g.target(e), g.source(e));
            original[rev] = false;
            weight[rev] = weight[e];
        }
    }
}

void Problem_Instance::clojure() {
    int curr_radius;
    graph_radius = INT_MAX;
    // Run a Dijkstra using each node as a source and then add an arc to each node with the distance:
    DijkstraSolver dijkstra_solver(g, weight);
    for (DNodeIt u(g); u != INVALID; ++u) {
        curr_radius = 0;
        dijkstra_solver.run(u);
        for (DNodeIt v(g); v != INVALID; ++v) {
            if (u == v) continue;

            int dist = dijkstra_solver.dist(v);
            auto aux = findArc(g, u, v);
            curr_radius = max(curr_radius, dist);// update the current radius
#ifndef BDHST
            if (v == source) source_radius = max(source_radius, dist);// update the source radius
#endif
            if (aux == INVALID) {// if this arc is not present, then add
                aux = g.addArc(u, v);
                original[aux] = false;// the tslib instances represent complete graphs
            }
            weight[aux] = dist;// update the weight to reflect the minimum cost path
        }
        graph_radius = min(graph_radius, curr_radius);// update the graph's radius
    }
}

bool Problem_Instance::view_solution(double LB, double UB, const string &msg, bool only_active_edges) {
    DigraphAttributes GA(g, vname, px, py);
    GA.SetDigraphAttrib("splines=true");
    GA.SetDefaultDNodeAttrib("color=gray style=filled shape=circle fixedsize=true");
    GA.SetDefaultArcAttrib("color=black arrowhead=none fontcolor=red");
    if (only_active_edges) GA.SetDefaultArcAttrib("style=invis");

    int i = 0;
    ArcVector used_arcs(nnodes - 1);
    for (ArcIt e(g); e != INVALID; ++e) {
        if (solution[e]) used_arcs[i++] = e;
    }
    for (i = 0; i < nnodes - 1; i++) {
        auto e = used_arcs[i];
        if (nnodes < 100) GA.SetLabel(e, weight[e] * MY_EPS);
        if (original[e]) e = g.addArc(g.source(e), g.target(e));// duplicate if already exists
        GA.SetColor(e, "red");
        GA.SetAttrib(e, "style=dashed arrowhead=normal");
    }
    auto max_height = 0;
    for (DNodeIt v(g); v != INVALID; ++v) max_height = max(max_height, node_activation[v]);
    for (DNodeIt v(g); v != INVALID; ++v)
        if (node_activation[v] == max_height) GA.SetColor(v, "cyan");// highlight the deepest nodes
    GA.SetColor(source, "pink");
#ifdef BDHST
    GA.SetLabel("Tree rooted at node " + vname[source] + " of height " + to_string(UB) + ". LB = " + to_string(LB) +
                ". " + msg);
#else
    GA.SetLabel("Scheduling starting from node " + vname[source] + " of makespan " + to_string(UB) +
                ". LB = " + to_string(LB) + ". " + msg);
#endif
    GA.View();
    return true;
}
