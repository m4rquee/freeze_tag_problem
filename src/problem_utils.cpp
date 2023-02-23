#include "problem_utils.hpp"

Problem_Instance::Problem_Instance(Digraph &graph, DNodeStringMap &vvname, DNodePosMap &posx, DNodePosMap &posy,
                                   DNode &sourcenode, int &nnodes, int &time_limit, ArcIntMap &pweight,
                                   ArcBoolMap &poriginal, int &psource_radius)
    : g(graph), vname(vvname), px(posx), py(posy), nnodes(nnodes), source(sourcenode), time_limit(time_limit),
      weight(pweight), original(poriginal), solution(g), source_radius(psource_radius), node_height(g) {
    for (ArcIt e(g); e != INVALID; ++e) {
        arc_map[g.source(e)][g.target(e)] = e;
        solution[e] = false;
    }
}

void Problem_Instance::start_counter() { start = chrono::system_clock::now(); }

void Problem_Instance::stop_counter() { stop = chrono::system_clock::now(); }

void PrintInstanceInfo(Problem_Instance &P) {
#ifdef BDHST
    cout << "Bounded Degree Minimum Height Spanning Tree graph information" << endl;
#else
    cout << "Freeze-Tag graph information" << endl;
#endif
    cout << "\tTime limit = " << P.time_limit << "s" << endl;
    cout << "\tNumber of nodes = " << P.nnodes << endl;
    if (P.source != INVALID) cout << "\tSource = " << P.vname[P.source] << endl;
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

bool ReadProblemGraph(const string &filename, Digraph &g, DNodeStringMap &vname, DNodePosMap &posx, DNodePosMap &posy,
                      DNode &source, int &nnodes, ArcIntMap &weight, ArcBoolMap &original, int &source_radius,
                      bool calc_clojure, bool tsplib) {
    if (tsplib) {
        if (!ReadTSPLIBDigraph(filename, g, vname, posx, posy, weight)) return false;
    } else {
        if (!ReadDigraph(filename, g, vname, posx, posy, weight)) return false;
        // Scale up to treat rational values as integers:
        for (ArcIt e(g); e != INVALID; ++e) weight[e] = (int) (weight[e] / MY_EPS);
    }
    nnodes = countNodes(g);
#ifndef BDHST
    source = g.nodeFromId(g.maxNodeId());
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

    source_radius = 0;
#ifdef BDHST
    int curr_radius;
    source_radius = INT_LEAST32_MAX;
#endif
    // Run a Dijkstra using each node as source and them add an arc to each node with the distance:
    DijkstraSolver dijkstra_solver(g, weight);
    for (DNodeIt u(g); u != INVALID; ++u) {
#ifdef BDHST
        curr_radius = 0.0;
#endif
        dijkstra_solver.run(u);
        for (DNodeIt v(g); v != INVALID; ++v) {
            if (u == v) continue;
#ifdef BDHST
            curr_radius = max(curr_radius, dijkstra_solver.dist(v));
#else
            if (v == source) source_radius = max(source_radius, dijkstra_solver.dist(v));
#endif
            auto aux = findArc(g, u, v);
            if (aux == INVALID) {// if this arc is not present then add it with the distance as weight
                aux = g.addArc(u, v);
                original[aux] = tsplib;// the tslib instances represents complete graphs
            }
            weight[aux] = dijkstra_solver.dist(v);
        }
#ifdef BDHST
        source_radius = min(source_radius, curr_radius);
#endif
    }
    return true;
}

bool ViewProblemSolution(Problem_Instance &P, double &LB, double &UB, const string &msg, bool &only_active_edges) {
    DigraphAttributes GA(P.g, P.vname, P.px, P.py);
    GA.SetDefaultDNodeAttrib("color=LightGray style=filled width=0.1 height=0.1 fixedsize=false");

    int i = 0;
    ArcVector used_arcs(P.nnodes - 1);
    for (ArcIt e(P.g); e != INVALID; ++e) {
        if (P.solution[e]) used_arcs[i++] = e;
        if (!P.original[e]) continue;
        GA.SetColor(e, only_active_edges ? "#00000000" : "#00000070");
        GA.SetAttrib(e, "style=dotted");
    }
    for (i = 0; i < P.nnodes - 1; i++) {
        auto e = used_arcs[i];
        auto weight = P.weight[e];
        if (P.original[e]) e = P.g.addArc(P.g.source(e), P.g.target(e));// duplicate if already exists
        GA.SetColor(e, "#ff000070");
        GA.SetAttrib(e, "style=dashed splines=true");
        if (P.nnodes < 100) GA.SetLabel(e, weight * MY_EPS);
    }
    for (DNodeIt v(P.g); v != INVALID; ++v)
        if (P.node_height[v] * MY_EPS == UB) GA.SetColor(v, "Cyan");// highlight the deepest nodes
    GA.SetColor(P.source, "Red");
#ifdef BDHST
    GA.SetLabel("Tree rooted at node " + P.vname[P.source] + " of height " + DoubleToString(UB) +
                ". LB = " + DoubleToString(LB) + ". " + msg);
#else
    GA.SetLabel("Scheduling starting from node " + P.vname[P.source] + " of makespan " + DoubleToString(UB) +
                ". LB = " + DoubleToString(LB) + ". " + msg);
#endif
    GA.View();
    return true;
}
