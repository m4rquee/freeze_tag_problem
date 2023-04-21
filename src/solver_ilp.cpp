#include "gurobi_c++.h"
#include "problem_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

int seed = 42;// seed to the random number generator

class LocalSearchCB : public GRBCallback {
    Problem_Instance &P;

    Digraph::ArcMap<GRBVar> &x_e;
    Digraph::NodeMap<GRBVar> &h_v;
    GRBVar &height;

    ArcBoolMap arc_value;
    DNodeValueMap node_height;
    DNodeIntMap node_degree;

    double (GRBCallback::*solution_value)(GRBVar) = nullptr;

    typedef pair<int, int> h_d_pair;
    typedef pair<h_d_pair, DNode> h_d_triple;
    typedef pair<int, DNode> h_pair;

    vector<h_pair> node_height_heap;     // max heap ordered by nodes height
    vector<h_d_triple> available_fathers;// min heap ordered by nodes height and then their degree

public:
    LocalSearchCB(Problem_Instance &_P, Digraph::ArcMap<GRBVar> &_x_e, Digraph::NodeMap<GRBVar> &_h_v, GRBVar &_height)
        : P(_P), x_e(_x_e), h_v(_h_v), height(_height), arc_value(P.g), node_height(P.g), node_degree(P.g),
          node_height_heap(P.nnodes - 1), available_fathers(P.nnodes - 1) {}

protected:
    inline DNode get_father(const DNode &v) {
        for (InArcIt e(P.g, v); e != INVALID; ++e)
            if (arc_value[e]) return P.g.source(e);
        return INVALID;
    }

    // Get the child's shallowest ancestor that still can be children of the new_father while improving its own cost:
    DNode get_shallowest_ancestor(const DNode &father, const DNode &child, const DNode &new_father) {
        // Cannot move a node to itself or its current father:
        if (child == new_father || father == new_father) return INVALID;

        DNode grandfather = get_father(father);
        double new_father_h = node_height[new_father];
        if (grandfather != P.source &&// reached the root and so moving the father will disconnect the tree
            new_father_h + P.weight[P.arc_map[new_father][father]] < node_height[father])
            return get_shallowest_ancestor(grandfather, father, new_father);// can keep going up
        if (new_father_h + P.weight[P.arc_map[new_father][child]] < node_height[child])
            return child;// found the shallowest ancestor that makes an improving swap
        return INVALID;  // there is no improving swap
    }

    void calc_height(DNode &v, double v_height) {
        if (node_degree[v] > 0)
            for (OutArcIt e(P.g, v); e != INVALID; ++e)
                if (arc_value[e]) {
                    auto child = P.g.target(e);
                    calc_height(child, node_height[child] = v_height + P.weight[e]);
                }
    }

    void heap_init() {
        node_height_heap.clear();
        available_fathers.clear();
        for (DNodeIt v(P.g); v != INVALID; ++v) {
            if (v == P.source) {
                node_degree[v] = 1;
                continue;// the source will not be swapped, nor be receiving new children
            }
            node_degree[v] = 0;
            for (OutArcIt e(P.g, v); e != INVALID; ++e) node_degree[v] += arc_value[e];
            node_height_heap.push_back(h_pair(node_height[v], v));
            if (node_degree[v] == 2) continue;// cannot receive new children
            available_fathers.push_back(h_d_triple(h_d_pair(-node_height[v], -node_degree[v]), v));
        }
        make_heap(node_height_heap.begin(), node_height_heap.end());
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
        auto found_better_sol = false;
        while (!available_fathers.empty()) {
            h_d_triple popped_triple = available_fathers.front();
            DNode &new_father = popped_triple.second;
            pop_heap(available_fathers.begin(), available_fathers.end());// moves the triple to the end
            available_fathers.pop_back();                                // deletes the triple from the heap

            // Check if it can make a deep node child of the candidate father: -----------------------------------------
            h_pair popped_pair;
            DNode shallowest_ancestor = INVALID, child;
            do {
                if (node_height_heap.empty())
                    break;// cannot make a node child of this new father while improving the cost
                popped_pair = node_height_heap.front();
                child = shallowest_ancestor = popped_pair.second;
                pop_heap(node_height_heap.begin(), node_height_heap.end());
                node_height_heap.pop_back();
                shallowest_ancestor = get_shallowest_ancestor(get_father(child), child, new_father);
            } while (shallowest_ancestor == INVALID);

            if (shallowest_ancestor == INVALID) continue;// could not find an improving swap with this candidate father
            found_better_sol = true;

            // Make shallowest_ancestor child of the new_father:
            auto new_arc = P.arc_map[new_father][shallowest_ancestor];
            setSolution(x_e[new_arc], arc_value[new_arc] = true);// connect shallowest_ancestor to the new_father
            for (InArcIt e(P.g, shallowest_ancestor); e != INVALID; ++e)
                if (arc_value[e]) {// remove the arc to the shallowest_ancestor from its old father
                    setSolution(x_e[e], arc_value[e] = false);
                    break;
                }
            node_height[shallowest_ancestor] = node_height[new_father] + P.weight[new_arc];
            calc_height(shallowest_ancestor, node_height[shallowest_ancestor]);// update the whole subtree
            heap_init();
        }


        // Has found an improving solution and so informs it:
        if (found_better_sol) {
            auto new_sol_height = 0.0;
            for (DNodeIt v(P.g); v != INVALID; ++v) {
                new_sol_height = max(new_sol_height, node_height[v]);
                setSolution(h_v[v], node_height[v]);
            }

            auto old_sol_h = MY_EPS * (this->*solution_value)(height);
            setSolution(height, new_sol_height);
            new_sol_height *= MY_EPS;
            if (old_sol_h != new_sol_height)
                cout << "\n→ Found solution with local search of height " << new_sol_height << " over " << old_sol_h;
            else
                cout << "\n→ Found solution with local search of same height, but better total cost";
            cout << "\n\n";
            useSolution();// informs gurobi of this solution
        }
    }
};

int calc_height(Problem_Instance &P, DNode &root) {
    P.node_height[root] = 0;
    for (OutArcIt e(P.g, root); e != INVALID; ++e)
        if (P.solution[e]) {
            auto child = P.g.target(e);
            int child_height = P.weight[P.arc_map[root][child]] + calc_height(P, child);
            P.node_height[root] = max(P.node_height[root], child_height);
        }
    return P.node_height[root];
}

double greedy_solution(Problem_Instance &P, int max_degree) {
    // Connect the source to some closest node:
    Arc min_arc = INVALID;
    double min_arc_weight = MY_INF;
    auto source = P.source != INVALID ? P.source : Digraph::nodeFromId(rand() % countNodes(P.g));
    for (OutArcIt e(P.g, source); e != INVALID; ++e)
        if (P.weight[e] < min_arc_weight) {
            min_arc = e;
            min_arc_weight = P.weight[e];
        }
    P.solution[min_arc] = true;

    // Init the degree map (-1 are not yet border nodes and -2 saturated nodes):
    DNodeIntMap degree(P.g);
    for (DNodeIt v(P.g); v != INVALID; ++v) degree[v] = -1;
#ifdef BDHST
    degree[source] = max_degree == 1 ? -2 : 1;
#else
    degree[source] = -2;
#endif
    degree[P.g.target(min_arc)] = 0;

    DNodeVector border;
    border.push_back(P.g.target(min_arc));
    for (int i = 1; i < P.nnodes - 1; i++) {
        min_arc_weight = MY_INF;
        for (DNode u: border)
            for (OutArcIt e(P.g, u); e != INVALID; ++e) {
                auto v = P.g.target(e);
                if (degree[v] == -1)// not yet added
                    if (P.weight[e] < min_arc_weight) {
                        min_arc = e;
                        min_arc_weight = P.weight[e];
                    }
            }
        P.solution[min_arc] = true;
        auto u = P.g.source(min_arc), v = P.g.target(min_arc);
        degree[u]++;
        if (degree[u] == max_degree - (u != source)) {// becomes saturated with max_degree-1 children
            auto aux = remove(border.begin(), border.end(), u);
            border.erase(aux, border.end());
        }
        degree[v] = 0;
        border.push_back(v);
    }
    return calc_height(P, source);
}

bool solve(Problem_Instance &P, double &LB, double &UB, int max_degree = 3) {
    P.start_counter();

    // Calculates the best know objective bounds:
    int MAX_EDGE = 0, DIAMETER = 0;
    for (ArcIt e(P.g); e != INVALID; ++e) {
        if (P.original[e]) MAX_EDGE = max(MAX_EDGE, P.weight[e]);
        DIAMETER = max(DIAMETER, P.weight[e]);
    }
    double minimum_depth = ceil(log(P.nnodes) / log(max_degree - 1));
    auto auxUB = MAX_EDGE * minimum_depth;
    LB = max(LB / MY_EPS, (double) P.source_radius);
    UB = max(LB, min(UB / MY_EPS, auxUB));
    // Construct an initial greedy solution:
    auto greedy_sol = greedy_solution(P, max_degree);
    UB = min(UB, greedy_sol);

    cout << "Set parameter LB to value " << LB * MY_EPS << endl;
    cout << "Set parameter UB to value " << UB * MY_EPS << endl;
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
    height.set(GRB_DoubleAttr_Start, greedy_sol);

    for (ArcIt e(P.g); e != INVALID; ++e) {
        char name[100];
        sprintf(name, "x_(%s,%s)", P.vname[P.g.source(e)].c_str(), P.vname[P.g.target(e)].c_str());
        x_e[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        x_e[e].set(GRB_DoubleAttr_Start, P.solution[e]);
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
                if (P.nnodes > 500) continue;// reduce the model size for big instances
                constrCount++;
                model.addConstr(h_v[v] <= h_v[P.g.source(e)] + P.weight[e] + (P.weight[e] + UB) * (1 - x_e[e]));
            }
    cout << "-> a node height is its parents height plus the edge to it - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt u(P.g); u != INVALID; ++u)
        for (OutArcIt e(P.g, u); e != INVALID; ++e) {
            DNode v = P.g.target(e);
            if (Digraph::id(u) < Digraph::id(v)) {
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
    if (!improved) return false;

    // A better solution was found:
    cout << endl << "Nodes height: ";
    UB = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        int node_height = ceil(h_v[v].get(GRB_DoubleAttr_X));
        P.node_height[v] = node_height;
        cout << P.vname[v].c_str() << '-' << node_height * MY_EPS << ";";
#ifdef BDHST
        if (r_v[v].get(GRB_DoubleAttr_X) >= 1 - MY_EPS) P.source = v;
#endif
        UB = max(UB, (double) node_height);
    }

    UB *= MY_EPS;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)// solved optimally
        LB = UB;

    // Display the best know solution:
    int i = 0;
    cout << endl << "Tree edges: ";
    for (ArcIt e(P.g); e != INVALID; ++e) {
        bool active = x_e[e].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
        P.solution[e] = active;
        if (active) cout << P.vname[P.g.source(e)].c_str() << '-' << P.vname[P.g.target(e)].c_str() << ";";
        if (!P.original[e] && !active) P.g.erase(e);
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

    using namespace std::chrono;
    seed = (int) duration_cast<minutes>(system_clock::now().time_since_epoch()).count();
    srand(seed);
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
