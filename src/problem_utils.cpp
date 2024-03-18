#include "problem_utils.hpp"

Problem_Instance::Problem_Instance(const string &filename, int time_limit, bool calc_clojure, bool tsplib)
    : time_limit(time_limit), vname(g), weight(g, MY_INF), px(g), py(g), original(g), solution(g, false),
      node_activation(g), tsplib(tsplib), complete(calc_clojure) {
    read_instance(filename, tsplib);

    // Initialize the attributes:
    source_radius = graph_radius = graph_diameter = -1;
    source = INVALID;
#ifndef BDHST
    source = Digraph::nodeFromId(0);
#endif
    if (calc_clojure) clojure();
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
    cout << "\tNumber of original edges = " << narcs << endl;
    if (source != INVALID) {
        cout << "\tSource = " << vname[source] << endl;
        cout << "\tSource radius = " << source_radius * MY_EPS << endl;
    }
    cout << "\tGraph radius = " << graph_radius * MY_EPS << endl;
    cout << "\tGraph diameter = " << graph_diameter * MY_EPS << endl;
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
        vname[u] = to_string(id);
        px[u] = x;
        py[u] = y;
        DNodeIt v(g);
        for (++v; v != INVALID; ++v) weight[g.addArc(u, v)] = euclidean_distance(u, v, px, py);
    }
}

void Problem_Instance::read_instance(const string &filename, bool tsplib) {
    if (tsplib) {
        read_tsplib_instance(filename);
        narcs = 0;
    } else {
        ReadDigraph(filename, g, vname, px, py, weight);
        for (ArcIt e(g); e != INVALID; ++e)
            weight[e] = floor(weight[e] / MY_EPS);// scale up so we can treat rational values as integers
        narcs = countArcs(g);
    }
    nnodes = countNodes(g);

    // Ensure the digraph is strongly connected by adding each arc reverse:
    int i = 0;
    vector<Arc> original_edges(countArcs(g));
    for (ArcIt e(g); e != INVALID; ++e, i++) original[original_edges[i] = e] = !tsplib;// save the original arcs
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
            curr_radius = max(curr_radius, dist);// update the current radius
            graph_diameter = max(graph_diameter, dist);
#ifndef BDHST
            if (u == source) source_radius = max(source_radius, dist);// update the source radius
#endif

            // Add a directed both way connection between u and v:
            auto aux = findArc(g, u, v);
            if (aux == INVALID) {// if this arc is not present, then add
                aux = g.addArc(u, v);
                original[aux] = false;
            }
            weight[aux] = dist;// update the weight to reflect the cost of a minimum path
            auto rev = findArc(g, v, u);
            if (rev == INVALID) {// if the reverse of this arc is not present, then add
                rev = g.addArc(v, u);
                original[rev] = false;
            }
            weight[rev] = dist;
        }
        graph_radius = min(graph_radius, curr_radius);// update the graph's radius
    }
}

void Problem_Instance::view_solution(double LB, double UB, const string &msg, bool only_active_edges) {
    DigraphAttributes GA(g, vname, px, py);
    GA.SetDefaultDNodeAttrib("color=gray style=filled shape=circle fixedsize=true");
    GA.SetDefaultArcAttrib("color=black arrowhead=none fontcolor=black style=invis fontsize=10");

    ArcVector used_arcs;// save used arcs to avoid arc duplication inside interation
    for (ArcIt e(g); e != INVALID; ++e) {
        if (solution[e]) used_arcs.push_back(e);
        if (!only_active_edges && original[e]) GA.SetAttrib(e, "style=solid");
    }
    for (auto &e: used_arcs) {
        int e_weight = weight[e];
        auto rev = findArc(g, g.target(e), g.source(e));
        auto rev_is_orig = rev != INVALID && original[rev];
        if (original[e] || rev_is_orig) e = g.addArc(g.source(e), g.target(e));// duplicate if already exists
        GA.SetColor(e, "red");
        GA.SetAttrib(e, "style=dashed arrowhead=normal");
        if (nnodes < 100) {
            std::ostringstream oss;
            oss << std::setprecision(ceil(log10(UB / MY_EPS))) << e_weight * MY_EPS;
            GA.SetLabel(e, oss.str());
        }
    }
    auto max_height = 0;
    for (DNodeIt v(g); v != INVALID; ++v) max_height = max(max_height, node_activation[v]);
    for (DNodeIt v(g); v != INVALID; ++v)
        if (node_activation[v] == max_height) GA.SetColor(v, "cyan");// highlight the deepest nodes
    GA.SetColor(source, "pink");
#ifdef BDHST
    GA.SetLabel("Tree rooted at node " + vname[source] + " of height " + to_string(UB) + ". LB = " + to_string(LB) +
                "." + msg);
#else
    GA.SetLabel("Scheduling starting from node " + vname[source] + " of makespan " + to_string(UB) +
                ". LB = " + to_string(LB) + "." + msg);
#endif
    MY_GRAPHLIB_PARAMETERS.graphviz_drawing_program = tsplib ? "neato" : "dot";
    GA.View();
}

double greedy_solution(Problem_Instance &P, int max_degree) {
    // Connect the source to some closest node:
    Arc min_arc = INVALID;
    int min_arc_weight = GRB_MAXINT;
    auto source = P.source != INVALID ? P.source : Digraph::nodeFromId(0);
    for (OutArcIt e(P.g, source); e != INVALID; ++e)
        if (P.weight[e] < min_arc_weight) {
            min_arc = e;
            min_arc_weight = P.weight[e];
        }
    P.solution[min_arc] = true;
    P.node_activation[source] = 0;
    auto target = P.g.target(min_arc);
    P.node_activation[target] = min_arc_weight;
    P.solution_makespan = min_arc_weight;

    // Init the degree map (-1 are not yet border nodes and -2 saturated nodes):
    DNodeIntMap degree(P.g);
    for (DNodeIt v(P.g); v != INVALID; ++v) degree[v] = -1;
#ifdef BDHST
    degree[source] = max_degree == 1 ? -2 : 1;
#else
    degree[source] = -2;
#endif
    degree[target] = 0;

    DNodeVector border;
    border.push_back(target);
    for (int i = 1; i < P.nnodes - 1; i++) {
        min_arc_weight = GRB_MAXINT;
        for (DNode u: border)
            for (OutArcIt e(P.g, u); e != INVALID; ++e) {
                auto v = P.g.target(e);
                if (degree[v] == -1)// not yet added
                    if (P.weight[e] < min_arc_weight) {
                        min_arc = e;
                        min_arc_weight = P.weight[e];
                    }
            }
        auto u = P.g.source(min_arc), v = P.g.target(min_arc);
        P.solution[min_arc] = true;
        P.node_activation[v] = P.node_activation[u] + min_arc_weight;
        P.solution_makespan = max(P.solution_makespan, P.node_activation[v]);
        degree[u]++;
        if (degree[u] == max_degree - (u != source)) {// becomes saturated with max_degree-1 children
            auto aux = remove(border.begin(), border.end(), u);
            border.erase(aux, border.end());
        }
        degree[v] = 0;
        border.push_back(v);
    }
    return P.solution_makespan;
}

LocalSearchCB::LocalSearchCB(Problem_Instance &_P, Digraph::ArcMap<GRBVar> &_x_e)
    : P(_P), x_e(_x_e), arc_value(P.g), node_depth(P.g), node_degree(P.g), node_depth_heap(P.nnodes - 1),
      available_fathers(P.nnodes - 1) {
    for (ArcIt e(P.g); e != INVALID; ++e)// starts all arcs values
        arc_value[e] = P.solution[e];
    for (DNodeIt v(P.g); v != INVALID; ++v)// starts all node's depth
        node_depth[v] = P.node_activation[v];
}

double LocalSearchCB::init() {
    where = GRB_CB_START;
    callback();
    return P.solution_makespan;
}

// Get the child's shallowest ancestor that still can be children of the new_father while improving its own cost:
DNode LocalSearchCB::get_shallowest_ancestor(const DNode &father, const DNode &child, const DNode &new_father) {
    // Cannot move a node to itself or its current father:
    if (child == new_father || father == new_father) return INVALID;

    DNode grandfather = get_father(father);
    double new_father_h = node_depth[new_father];
    if (grandfather != P.source) {// reached the root and so moving the father will disconnect the tree
        auto it = P.arc_map[new_father].find(father);
        if (it != P.arc_map[new_father].end())// if this arc exists
            if (new_father_h + P.weight[P.arc_map[new_father][father]] < node_depth[father])
                return get_shallowest_ancestor(grandfather, father, new_father);// can keep going up
    }
    auto it = P.arc_map[new_father].find(child);
    if (it != P.arc_map[new_father].end())// if this arc exists
        if (new_father_h + P.weight[P.arc_map[new_father][child]] < node_depth[child])
            return child;// found the shallowest ancestor that makes an improving swap
    return INVALID;      // there is no improving swap
}

void LocalSearchCB::calc_depth(DNode &v, double v_depth) {
    for (OutArcIt e(P.g, v); e != INVALID; ++e) {// todo: bug with the online edge update
        (this->*set_solution)(x_e[e], arc_value[e]);
        auto rev = findArc(P.g, P.g.target(e), P.g.source(e));
        (this->*set_solution)(x_e[rev], false);
        if (arc_value[e]) {
            auto child = P.g.target(e);
            calc_depth(child, node_depth[child] = v_depth + P.weight[e]);
        }
    }
}

void LocalSearchCB::heap_init() {
    node_depth_heap.clear();
    available_fathers.clear();
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        if (v == P.source) {
            node_degree[v] = 1;
            continue;// the source will not be swapped, nor be receiving new children
        }
        node_degree[v] = 0;
        for (OutArcIt e(P.g, v); e != INVALID; ++e) node_degree[v] += arc_value[e];
        node_depth_heap.emplace_back(node_depth[v], v);
        if (node_degree[v] >= 2) continue;// cannot receive new children
        available_fathers.emplace_back(ii_pair(-node_depth[v], -node_degree[v]), v);
    }
    make_heap(node_depth_heap.begin(), node_depth_heap.end());
    make_heap(available_fathers.begin(), available_fathers.end());
}

void LocalSearchCB::callback() {
    // Get the correct function to obtain the values of the lp variables:
    if (where == GRB_CB_MIPSOL) {// all variables are integer
        solution_value = &LocalSearchCB::getSolution;
        set_solution = &LocalSearchCB::setSolution;
        for (ArcIt e(P.g); e != INVALID; ++e)// saves all arcs values
            arc_value[e] = (this->*solution_value)(x_e[e]) >= 1 - MY_EPS;
    } else if (where == GRB_CB_START)
        set_solution = &LocalSearchCB::setSolutionStart;
    else
        return;// this callback does not take advantage of the other options

    heap_init();
    bool found_new_sol = false;
    while (!available_fathers.empty()) {
        iin_triple popped_triple = available_fathers.front();
        DNode &new_father = popped_triple.second;
        pop_heap(available_fathers.begin(), available_fathers.end());// moves the triple to the end
        available_fathers.pop_back();                                // deletes the triple from the heap

        // Check if it can make a deep node child of the candidate father: -----------------------------------------
        in_pair popped_pair;
        DNode shallowest_ancestor = INVALID, child;
        do {
            if (node_depth_heap.empty()) break;// cannot make a node child of this new father while improving the cost
            popped_pair = node_depth_heap.front();
            child = shallowest_ancestor = popped_pair.second;
            pop_heap(node_depth_heap.begin(), node_depth_heap.end());
            node_depth_heap.pop_back();
            shallowest_ancestor = get_shallowest_ancestor(get_father(child), child, new_father);
        } while (shallowest_ancestor == INVALID);

        if (shallowest_ancestor == INVALID) continue;// could not find an improving swap with this candidate father
        found_new_sol = true;

        // Make the shallowest_ancestor child of the new_father: ---------------------------------------------------
        Arc rev;
        for (InArcIt e(P.g, shallowest_ancestor); e != INVALID; ++e)
            if (arc_value[e]) {// remove the arc to the shallowest_ancestor from its old father
                (this->*set_solution)(x_e[e], arc_value[e] = false);
                rev = findArc(P.g, P.g.target(e), P.g.source(e));
                (this->*set_solution)(x_e[rev], false);
                break;
            }
        // Connect shallowest_ancestor to the new_father:
        auto new_arc = P.arc_map[new_father][shallowest_ancestor];
        (this->*set_solution)(x_e[new_arc], arc_value[new_arc] = true);
        rev = findArc(P.g, P.g.target(new_arc), P.g.source(new_arc));
        (this->*set_solution)(x_e[rev], false);

        node_depth[shallowest_ancestor] = node_depth[new_father] + P.weight[new_arc];
        calc_depth(shallowest_ancestor, node_depth[shallowest_ancestor]);// update the whole subtree
        heap_init();
    }

    // Has found an improving solution and so informs it:
    if (found_new_sol) {
        auto new_sol_h = 0.0;
        for (DNodeIt v(P.g); v != INVALID; ++v) new_sol_h = max(new_sol_h, node_depth[v]);
        if (new_sol_h < P.solution_makespan) {// copy the solution if its better:
            for (ArcIt e(P.g); e != INVALID; ++e) P.solution[e] = arc_value[e];
            for (DNodeIt v(P.g); v != INVALID; ++v) P.node_activation[v] = (int) node_depth[v];
            cout << "\n→ Found a solution with local search of height " << MY_EPS * new_sol_h << " over "
                 << MY_EPS * P.solution_makespan;
            P.solution_makespan = new_sol_h;
        } else
            cout << "\n→ Found a solution with local search of same height, but better total cost";
        cout << "\n\n";
        useSolution();// informs gurobi of this solution
    }
}
