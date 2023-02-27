#include "gurobi_c++.h"
#include "problem_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

const long unsigned seed = 42;// seed to the random number generator

class LocalSearchCB : public GRBCallback {
    Problem_Instance &P;

    Digraph::ArcMap<GRBVar> &x_e;
    Digraph::NodeMap<GRBVar> &h_v;
    GRBVar &height;

    ArcBoolMap arc_value;
    DNodeValueMap node_height;
    DNodeIntMap node_degree;

    double (GRBCallback::*solution_value)(GRBVar) = nullptr;

    typedef pair<int, int> key;
    typedef pair<key, DNode> entry;
    typedef pair<int, DNode> rank;

    vector<rank> node_rank;
    vector<entry> available_fathers;

public:
    LocalSearchCB(Problem_Instance &_P, Digraph::ArcMap<GRBVar> &_x_e, Digraph::NodeMap<GRBVar> &_h_v, GRBVar &_height)
        : P(_P), x_e(_x_e), h_v(_h_v), height(_height), arc_value(P.g), node_height(P.g), node_degree(P.g),
          node_rank(P.nnodes - 1), available_fathers(P.nnodes - 1) {}

protected:
    inline DNode get_father(const DNode &node) {
        for (InArcIt e(P.g, node); e != INVALID; ++e)
            if (arc_value[e]) return P.g.source(e);
        return INVALID;
    }

    DNode get_shallowest_ancestor(const DNode &father, const DNode &child, const DNode &new_father) {
        // Cannot move a node to itself, or to its current father:
        if (child == new_father || father == new_father) return INVALID;

        DNode next_father = get_father(father);
        double &new_father_h = node_height[new_father];
        if (next_father != P.source &&// reached the root and so moving the father will disconnect the tree
            new_father_h + P.weight[P.arc_map[new_father][father]] < node_height[father])
            return get_shallowest_ancestor(next_father, father, new_father);// can keep going up
        if (new_father_h + P.weight[P.arc_map[new_father][child]] < node_height[child])
            return child;// found the shallowest ancestor that makes an improving swap
        return INVALID;  // there is no improving swap
    }

    void calc_height(DNode &v, double &v_height) {
        if (node_degree[v] > 0)
            for (OutArcIt e(P.g, v); e != INVALID; ++e)
                if (arc_value[e]) {
                    auto child = P.g.target(e);
                    calc_height(child, node_height[child] = v_height + P.weight[e]);
                }
    }

    void heap_init() {
        node_rank.clear();
        available_fathers.clear();
        for (DNodeIt v(P.g); v != INVALID; ++v) {
            node_degree[v] = 0;
            for (OutArcIt e(P.g, v); e != INVALID; ++e) node_degree[v] += arc_value[e];
            if (v == P.source) continue;// the source will not be swapped
            node_rank.push_back(rank(node_height[v], v));
            available_fathers.push_back(entry(key(-node_degree[v], node_height[v]), v));
        }
        make_heap(node_rank.begin(), node_rank.end());
        make_heap(available_fathers.begin(), available_fathers.end());
    }

    void callback() override {
        // Get the correct function to obtain the values of the lp variables:
        if (where == GRB_CB_MIPSOL)// all variables are integer
            solution_value = &LocalSearchCB::getSolution;
        else
            return;// this code do not take advantage of the other options

        for (ArcIt e(P.g); e != INVALID; ++e)// saves all arcs values
            arc_value[e] = (this->*solution_value)(x_e[e]) >= 1 - MY_EPS;
        for (DNodeIt v(P.g); v != INVALID; ++v)// saves all nodes height
            node_height[v] = (this->*solution_value)(h_v[v]);

        heap_init();
        auto newSolH = MY_INF;
        while (!available_fathers.empty()) {
            entry &popped_entry = available_fathers.front();
            DNode &new_father = popped_entry.second;
            pop_heap(available_fathers.begin(), available_fathers.end());// moves the entry to the end
            available_fathers.pop_back();                                // deletes the entry from the heap

            if (node_degree[new_father] == 2) break;// no father with available slot

            // Can make a deep node child of this one: -----------------------------------------------------------------
            rank popped_rank;
            DNode shallowest_ancestor, child;
            do {
                if (node_rank.empty()) break;// cannot make a node child of this new father while improving the cost
                popped_rank = node_rank.front();
                child = shallowest_ancestor = popped_rank.second;
                pop_heap(node_rank.begin(), node_rank.end());
                node_rank.pop_back();
            } while ((shallowest_ancestor = get_shallowest_ancestor(get_father(child), child, new_father)) == INVALID);

            if (shallowest_ancestor == INVALID) continue;// could not find an improving swap with this father

            // Make shallowest_ancestor child of the new_father:
            auto new_arc = P.arc_map[new_father][shallowest_ancestor];
            setSolution(x_e[new_arc], arc_value[new_arc] = true);// connect shallowest_ancestor to the new_father
            for (InArcIt e(P.g, shallowest_ancestor); e != INVALID; ++e)
                if (arc_value[e]) {// remove the arc from the old father to the shallowest_ancestor
                    setSolution(x_e[e], arc_value[e] = false);
                    break;
                }
            node_height[shallowest_ancestor] = node_height[new_father] + P.weight[new_arc];
            calc_height(shallowest_ancestor, node_height[shallowest_ancestor]);// update the whole subtree

            // Reset the heaps and set the solution:
            heap_init();
            newSolH = 0.0;
            for (DNodeIt v(P.g); v != INVALID; ++v)
                if (v != P.source) {
                    newSolH = max(newSolH, node_height[v]);
                    setSolution(h_v[v], node_height[v]);
                }
        }

        if (newSolH < MY_INF) {
            setSolution(height, newSolH);
            auto oldSolH = MY_EPS * (this->*solution_value)(height);
            newSolH *= MY_EPS;
            cout << "\n→ Solution found with local search of height " << newSolH << " over " << oldSolH << "\n\n";
            useSolution();// informs gurobi of this solution
        }
    }
};

int calc_height(Problem_Instance &P, DNode &root) {
    DNode left = INVALID, right = INVALID;
    for (OutArcIt e(P.g, root); e != INVALID; ++e)
        if (P.solution[e]) {
            if (left == INVALID) left = P.g.target(e);
            else {
                right = P.g.target(e);
                break;
            }
        }
    if (left == INVALID) return 0;
    int left_h = P.weight[P.arc_map[root][left]] + calc_height(P, left);
    if (right == INVALID) return left_h;
    int right_h = P.weight[P.arc_map[root][right]] + calc_height(P, right);
    return max(left_h, right_h);
}

double greedy_solution(Problem_Instance &P) {
    // Connect the source to the closest node:
    Arc min_arc = INVALID;
    double min_arc_weight = MY_INF;
    for (DNodeIt v(P.g); v != INVALID; ++v)
        if (v != P.source && P.weight[P.arc_map[P.source][v]] < min_arc_weight) {
            min_arc = P.arc_map[P.source][v];
            min_arc_weight = P.weight[min_arc];
        }
    P.solution[min_arc] = true;

    // Init the degree map (-1 are not yet added nodes and -2 saturated nodes):
    DNodeIntMap degree(P.g);
    for (DNodeIt v(P.g); v != INVALID; ++v) degree[v] = -1;
    degree[P.source] = -2;
    degree[P.g.target(min_arc)] = 0;

    DNodeVector added;
    added.push_back(P.g.target(min_arc));
    for (int i = 1; i < P.nnodes - 1; i++) {
        min_arc_weight = MY_INF;
        for (DNode &u: added)
            if (degree[u] >= 0)// already added and not saturated
                for (DNodeIt v(P.g); v != INVALID; ++v)
                    if (degree[v] == -1)// not yet added
                        if (P.weight[P.arc_map[u][v]] < min_arc_weight) {
                            min_arc = P.arc_map[u][v];
                            min_arc_weight = P.weight[min_arc];
                        }
        P.solution[min_arc] = true;
        auto u = P.g.source(min_arc), v = P.g.target(min_arc);
        if (++degree[u] == 2) degree[u] *= -1;// becomes saturated with two children
        degree[v]++;
        added.push_back(v);
    }
    return calc_height(P, P.source);
}

bool solve(Problem_Instance &P, double &LB, double &UB, int max_degree = 3) {
    P.start_counter();

    // Calculates the best know objective bounds:
    int MAX_EDGE = 0, DIAMETER = 0;
    for (ArcIt e(P.g); e != INVALID; ++e) {
        if (P.original[e]) MAX_EDGE = max(MAX_EDGE, P.weight[e]);
        DIAMETER = max(DIAMETER, P.weight[e]);
    }
#ifdef BDHST
    double auxUB = P.source_radius * log(P.nnodes) / log(max_degree - 1);
#else
    double auxUB = 2 * P.source_radius * log2(P.nnodes);
#endif
    LB = max(LB / MY_EPS, (double) P.source_radius);// the source radius if there is one or the graph`s radius
    UB = max(LB, min(UB / MY_EPS, ceil(auxUB)));
    // Construct an initial greedy solution:
#ifndef BDHST
    UB = min(UB, greedy_solution(P));
#endif

    cout << "Set parameter LB to value " << LB << endl;
    cout << "Set parameter UB to value " << UB << endl;
    auto COST_MULTIPLIER = pow(10, ceil(log10((double) P.nnodes * UB)));// make a tens power
    cout << "Set parameter COST_MULTIPLIER to value " << COST_MULTIPLIER << endl;

    // Gurobi ILP problem setup:
    auto *env = new GRBEnv();
    env->set(GRB_IntParam_Seed, seed);
    env->set(GRB_DoubleParam_TimeLimit, P.time_limit);
    env->set(GRB_DoubleParam_Cutoff, (UB + 1) * COST_MULTIPLIER);// set the best know UB
    GRBModel model = GRBModel(*env);
#ifdef BDHST
    model.set(GRB_StringAttr_ModelName, "BDHST Problem");
#else
    model.set(GRB_StringAttr_ModelName, "Freeze-Tag Problem");
#endif
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_OptimalityTol, MY_EPS);

    // ILP solver parameters: ----------------------------------------------------
    if (P.nnodes >= 300) model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    else
        model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_OPTIMALITY);
    model.set(GRB_IntParam_Cuts, GRB_CUTS_AGGRESSIVE);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    model.set(GRB_DoubleParam_Heuristics, 0.25);

    // ILP problem variables: ----------------------------------------------------
    Digraph::ArcMap<GRBVar> x_e(P.g);                                          // if arc e is present in the tree
    Digraph::NodeMap<GRBVar> h_v(P.g);                                         // height of node v
    auto height = model.addVar(LB, UB, COST_MULTIPLIER, GRB_INTEGER, "height");// tree's height
#ifdef BDHST
    Digraph::NodeMap<GRBVar> r_v(P.g);// if node v is the root
#endif

    for (ArcIt e(P.g); e != INVALID; ++e) {
        char name[100];
        sprintf(name, "x_(%s,%s)", P.vname[P.g.source(e)].c_str(), P.vname[P.g.target(e)].c_str());
        x_e[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
#ifndef BDHST
        x_e[e].set(GRB_DoubleAttr_Start, P.solution[e]);
#endif
    }
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        char name[100];
        sprintf(name, "h_%s", P.vname[v].c_str());
        h_v[v] = model.addVar(0.0, UB, 1.0, GRB_INTEGER, name);
#ifdef BDHST
        sprintf(name, "s_%s", P.vname[v].c_str());
        r_v[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
#endif
    }
    model.update();// run update to use model inserted variables

    // ILP problem restrictions: -------------------------------------------------
    cout << "Adding the model restrictions:" << endl;

    // A node height is at least the source distance to it:
    int constrCount = 0;
    if (P.source != INVALID) {
        for (DNodeIt v(P.g); v != INVALID; ++v, constrCount++)
            if (v != P.source) model.addConstr(h_v[v] >= P.weight[P.arc_map[P.source][v]]);
        cout << "-> a node height is at least the source distance to it - " << constrCount - 1 << " constrs" << endl;
    }

    // There is only one root:
#ifdef BDHST
    GRBLinExpr one_source_expr;
    for (DNodeIt v(P.g); v != INVALID; ++v) one_source_expr += r_v[v];
    model.addConstr(one_source_expr == 1);
    cout << "-> there is only one root - " << 1 << " constrs" << endl;
#endif

    // The source is the root and has out-degree one:
#ifndef BDHST
    GRBLinExpr s_out_degree_expr;
    for (OutArcIt e(P.g, P.source); e != INVALID; ++e) s_out_degree_expr += x_e[e];
    model.addConstr(s_out_degree_expr == 1);
    GRBLinExpr s_in_degree_expr;
    for (InArcIt e(P.g, P.source); e != INVALID; ++e) s_in_degree_expr += x_e[e];
    model.addConstr(s_in_degree_expr == 0);
    cout << "-> the source is the root and has out-degree one - " << 2 << " constrs" << endl;
#endif

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v) {
#ifndef BDHST
        if (v == P.source) continue;
#endif
        GRBLinExpr out_degree_expr, in_degree_expr;
        for (OutArcIt e(P.g, v); e != INVALID; ++e) out_degree_expr += x_e[e];
        for (InArcIt e(P.g, v); e != INVALID; ++e) in_degree_expr += x_e[e];
#ifdef BDHST
        // The in-degree is one and the out-degree is at most max_degree - 1 for each internal node, and
        // zero and max_degree respectively for the root:
        model.addConstr(out_degree_expr <= max_degree - 1 + r_v[v]);
        model.addConstr(in_degree_expr == 1 - r_v[v]);
#else
        // The in-degree is one and the out-degree is at most two for each internal node:
        model.addConstr(out_degree_expr <= 2);
        model.addConstr(in_degree_expr == 1);
#endif
        constrCount += 2;
    }
#ifdef BDHST
    cout << "-> the in-degree is one and the out-degree is at most max_degree - 1 for each internal node, and"
         << "\n   zero and max_degree respectively for the root - " << constrCount << " constrs" << endl;
#else
    cout << "-> the in-degree is one and the out-degree is at most two for each internal node - " << constrCount
         << " constrs" << endl;
#endif

    GRBLinExpr arc_sum_expr;
    for (ArcIt e(P.g); e != INVALID; ++e) arc_sum_expr += x_e[e];
    model.addConstr(arc_sum_expr == P.nnodes - 1);
    cout << "-> the number of edges is n-1 for any tree - " << 1 << " constrs" << endl;

    if (P.source != INVALID) {
        model.addConstr(h_v[P.source] == 0);
        constrCount = 1;
    }
#ifdef BDHST
    else {
        constrCount = 0;
        for (DNodeIt v(P.g); v != INVALID; ++v, constrCount++) model.addConstr(h_v[v] <= UB * (1 - r_v[v]));
    }
#endif
    cout << "-> the root is at height zero - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v, constrCount++) model.addConstr(height >= h_v[v]);
    cout << "-> the tree height is the maximum of each of its node's height - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v)
        if (v != P.source)
            for (InArcIt e(P.g, v); e != INVALID; ++e) {
                constrCount++;
                model.addConstr(h_v[v] >= h_v[P.g.source(e)] + P.weight[e] + (P.weight[e] + UB) * (x_e[e] - 1));
                if (P.nnodes > 500) continue;
                constrCount++;
                model.addConstr(h_v[v] <= h_v[P.g.source(e)] + P.weight[e] + (P.weight[e] + UB) * (1 - x_e[e]));
            }
    cout << "-> a node height is its parents height plus the edge to it - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt u(P.g); u != INVALID; ++u)
        for (OutArcIt e(P.g, u); e != INVALID; ++e) {
            DNode v = P.g.target(e);
            if (P.g.id(u) < P.g.id(v)) {
                constrCount++;
                model.addConstr(x_e[e] + x_e[findArc(P.g, v, u)] <= 1);
            }
        }
    cout << "-> only one of a parallel arc pair is allowed - " << constrCount << " constrs" << endl;

    // Callback setup: -----------------------------------------------------------
#ifndef BDHST
    LocalSearchCB cb(P, x_e, h_v, height);
    model.setCallback(&cb);
#endif

    // ILP solving: --------------------------------------------------------------
    model.optimize();// trys to solve optimally within the time limit
    P.stop_counter();

    LB = max(LB, ceil(model.get(GRB_DoubleAttr_ObjBound)) / COST_MULTIPLIER) * MY_EPS;
    bool improved = model.get(GRB_IntAttr_SolCount) > 0;
    if (improved) {// a better solution was found
        UB = ceil(height.get(GRB_DoubleAttr_X)) * MY_EPS;
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)// solved optimally
            LB = UB;
    } else
        return false;

    // Display the best know solution:
    int i = 0;
    cout << endl << "Tree edges: ";
    for (ArcIt e(P.g); e != INVALID; ++e) {
        bool active = x_e[e].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
        P.solution[e] = active;
        if (active) cout << P.vname[P.g.source(e)].c_str() << '-' << P.vname[P.g.target(e)].c_str() << ";";
        if (!P.original[e] && !active) P.g.erase(e);
    }
    cout << endl << "Nodes height: ";
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        int node_height = ceil(h_v[v].get(GRB_DoubleAttr_X));
        P.node_height[v] = node_height;
        cout << P.vname[v].c_str() << '-' << node_height * MY_EPS << ";";
#ifdef BDHST
        if (r_v[v].get(GRB_DoubleAttr_X) >= 1 - MY_EPS) P.source = v;
#endif
    }
    cout << endl;
    cout << "New LB: " << LB << endl;

    return improved;
}

int main(int argc, char *argv[]) {
    int maxtime;
    bool only_active_edges = true, tsplib = false;
    Digraph g;// graph declaration
    string graph_filename, source_node_name;
    DNodeStringMap vname(g); // name of graph nodes
    DNodePosMap px(g), py(g);// xy-coordinates for each node
    DNodeColorMap vcolor(g); // color of nodes
    ArcStringMap ename(g);   // name for graph arcs
    ArcColorMap ecolor(g);   // color of arcs
    ArcValueMap lpvar(g);    // used to obtain the contents of the LP variables
    ArcIntMap weight(g);     // arc weights
    ArcBoolMap original(g);  // if an arc is original
    vector<DNode> V;

    set_pdfreader("xdg-open");// the Linux will choose the default one
    if (argc < 3) {
        cout << endl
             << "Integer Linear Program for the Freeze-Tag Problem using the Gurobi solver;" << endl
             << "Usage: " << argv[0]
             << "  <ftp_graph_filename> <maximum_time_sec> <ftp_graph_filename> [-tsplib=true|false]"
                " [-only_active_edges=true|false] [LB] [UB]"
             << endl
             << endl;
        cout << "Example:" << endl
             << "\t" << argv[0] << " "
             << "../instances/small/tree-003.g 10" << endl
             << endl;
        exit(0);
    }

    graph_filename = argv[1];
    maxtime = atoi(argv[2]);
    if (argc >= 4) tsplib = strcmp(argv[3], "-tsplib=true") == 0;
    if (argc >= 5) only_active_edges = strcmp(argv[3], "-only_active_edges=true") == 0;
    MY_EPS = 1E-2;
    double LB = 0, UB = MY_INF;// consider MY_INF as infinity.
    if (argc >= 6) LB = atof(argv[5]);
    if (argc >= 7) UB = atof(argv[6]);

    DNode source;
    int nnodes, source_radius;
    if (!ReadProblemGraph(graph_filename, g, vname, px, py, source, nnodes, weight, original, source_radius, true,
                          tsplib)) {
        cout << "Error while reding the input graph." << endl;
        exit(EXIT_FAILURE);
    }

    Problem_Instance P(g, vname, px, py, source, nnodes, maxtime, weight, original, source_radius);
    PrintInstanceInfo(P);

    try {
        if (solve(P, LB, UB)) {
            char msg[100];
            sprintf(msg, " Gap of %.2f%%.", 100 * (UB - LB) / UB);
            ViewProblemSolution(P, LB, UB, msg, only_active_edges);
            cout << "UB cost: " << UB << endl;
        }
    } catch (std::exception &e) {
        cout << "UB cost: " << UB << endl;
        cerr << "\nException: " << e.what() << endl;
        return 1;
    } catch (GRBException &e) {
        cout << "UB cost: " << UB << endl;
        cerr << "\nGRBException: " << e.getMessage() << endl;
        return 1;
    }
    return 0;
}
