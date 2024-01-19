#include "problem_utils.hpp"

Problem_Instance::Problem_Instance(const string &filename, int time_limit, bool calc_clojure, bool tsplib)
    : time_limit(time_limit), vname(g), weight(g), px(g), py(g), original(g), solution(g), node_activation(g) {
    if (!read_instance(filename, calc_clojure, tsplib)) {
        cout << "Error while reding the input graph." << endl;
        exit(EXIT_FAILURE);
    }

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
    cout << "Bounded Degree Minimum Height Spanning Tree graph information" << endl;
#else
    cout << "Freeze-Tag graph information" << endl;
#endif
    cout << "\tTime limit = " << time_limit << "s" << endl;
    cout << "\tNumber of nodes = " << nnodes << endl;
    if (source != INVALID) cout << "\tSource = " << vname[source] << endl;
    if (source != INVALID) cout << "\tSource radius = " << radius * MY_EPS << endl;
    cout << endl;
}

inline int node_distance(DNode &u, DNode &v, DNodePosMap &posx, DNodePosMap &posy) {
    auto dist = sqrt(pow(posx[u] - posx[v], 2) + pow(posy[u] - posy[v], 2));
    return floor(dist / MY_EPS);// scale up to treat rational values as integers
}

bool ReadTSPLIBDigraph(const string &filename, Digraph &g, DNodeStringMap &vname, DNodePosMap &posx, DNodePosMap &posy,
                       ArcIntMap &weight) {
    ifstream file;
    file.open(filename.c_str());
    if (!file) {
        cout << "Error: Could not open file " << filename << "." << endl;
        exit(EXIT_FAILURE);
    }

    try {
        string line;
        do {// skip the header lines
            std::getline(file, line);
        } while (line != "NODE_COORD_SECTION");

        int id;
        double x, y;
        for (std::getline(file, line); line != "EOF"; std::getline(file, line)) {
            sscanf(line.c_str(), "%d %lf %lf", &id, &x, &y);
            auto u = g.addNode();
            vname[u] = to_string(id - 1);
            posx[u] = x;
            posy[u] = y;
            DNodeIt v(g);
            for (++v; v != INVALID; ++v) weight[g.addArc(u, v)] = node_distance(u, v, posx, posy);
        }
    } catch (...) { return false; }
    return true;
}

bool Problem_Instance::read_instance(const string &filename, bool calc_clojure, bool tsplib) {
    if (tsplib) {
        if (!ReadTSPLIBDigraph(filename, g, vname, px, py, weight)) return false;
    } else {
        if (!ReadDigraph(filename, g, vname, px, py, weight)) return false;
        // Scale up to treat rational values as integers:
        for (ArcIt e(g); e != INVALID; ++e) weight[e] = (int) (weight[e] / MY_EPS);
    }
    nnodes = countNodes(g);
#ifndef BDHST
    if (tsplib) source = Digraph::nodeFromId(g.maxNodeId());
    else
        source = Digraph::nodeFromId(0);
#else
    source = INVALID;
#endif

    // Ensure the digraph is strongly connected by adding each arc reverse:
    int i = 0;
    vector<Arc> original_edges(countArcs(g));
    for (ArcIt e(g); e != INVALID; ++e, i++) original[original_edges[i] = e] = true;// mark the originals
    for (i--; i >= 0; i--) {
        auto e = original_edges[i];
        auto rev = findArc(g, g.target(e), g.source(e));
        if (rev == INVALID) {
            rev = g.addArc(g.target(e), g.source(e));
            original[rev] = false;
            weight[rev] = weight[e];
        }
    }

    if (!calc_clojure) return true;

    radius = 0;
#ifdef BDHST
    int curr_radius;
    radius = INT_LEAST32_MAX;
#endif
    // Run a Dijkstra using each node as source and them add an arc to each node with the distance:
    DijkstraSolver dijkstra_solver(g, weight);
    for (DNodeIt u(g); u != INVALID; ++u) {
#ifdef BDHST
        curr_radius = 0.0;
#endif
        if (!tsplib) dijkstra_solver.run(u);
        for (DNodeIt v(g); v != INVALID; ++v) {
            if (u == v) continue;

            int dist;
            auto aux = findArc(g, u, v);
            if (tsplib) dist = weight[aux];
            else
                dist = dijkstra_solver.dist(v);
#ifdef BDHST
            curr_radius = max(curr_radius, dist);
#else
            if (v == source) radius = max(radius, dist);
#endif
            if (tsplib) continue;// no need to execute the following
            if (aux == INVALID) {// if this arc is not present then add it with the distance as weight
                aux = g.addArc(u, v);
                original[aux] = tsplib;// the tslib instances represents complete graphs
            }
            weight[aux] = dist;
        }
#ifdef BDHST
        radius = min(radius, curr_radius);
#endif
    }
    return true;
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
