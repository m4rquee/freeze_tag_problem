#include "problem_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

int seed = 42;// seed to the random number generator

bool solve(Problem_Instance &P, double &LB, double &UB, int max_degree = 3) {
    P.start_counter();

    // Calculates the best known objective bounds: -------------------------------
    double minimum_depth = ceil(log(P.nnodes) / log(2));
    LB = max(LB / MY_EPS, (double) P.source_radius);

    auto auxUB = P.graph_diameter * minimum_depth;
    bool improved = auxUB < UB / MY_EPS;
    UB = max(LB, min(UB / MY_EPS, auxUB));
    // Construct an initial greedy solution:
    auto greedy_UB = greedy_solution(P, max_degree);
    improved |= greedy_UB < UB;
    UB = min(UB, greedy_UB);

    // Gurobi ILP problem setup: -------------------------------------------------
    auto *env = new GRBEnv();
    env->set(GRB_IntParam_Seed, seed);
    env->set(GRB_DoubleParam_TimeLimit, P.time_limit);
    GRBModel model = GRBModel(*env);
    model.set(GRB_StringAttr_ModelName, "Freeze-Tag Problem");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_DoubleParam_OptimalityTol, MY_EPS);

    // ILP solver parameters: ----------------------------------------------------
    /*model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_OPTIMALITY);
    model.set(GRB_IntParam_Cuts, GRB_CUTS_AGGRESSIVE);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    model.set(GRB_DoubleParam_Heuristics, 0.25);*/
    model.set(GRB_IntParam_Threads, 4);

    // ILP problem variables: ----------------------------------------------------
    Digraph::ArcMap<GRBVar> x_e(P.g); // if arc e is present in the solution tree
    Digraph::NodeMap<GRBVar> d_v(P.g);// depth of node v
    Digraph::NodeMap<GRBVar> l_v(P.g);// if node v is a leaf
    GRBVar depth;                     // the solution tree depth

    // Callback setup: -----------------------------------------------------------
    LocalSearchCB cb(P, x_e);
    model.setCallback(&cb);
    UB = min(UB, cb.init());// try to improve the initial solution with local search

    // Cutoff and bounds setup: --------------------------------------------------
    cout << "Set parameter LB to value " << LB * MY_EPS << endl;
    cout << "Set parameter UB to value " << UB * MY_EPS << endl;
    auto cutoff = UB + 1;
    model.set(GRB_DoubleParam_Cutoff, cutoff);// set the best know UB

    // ILP problem variables startup: --------------------------------------------
    depth = model.addVar(LB, UB, 1.0, GRB_INTEGER, "depth");
    depth.set(GRB_DoubleAttr_Start, P.solution_makespan);

    for (ArcIt e(P.g); e != INVALID; ++e) {
        char name[100];
        sprintf(name, "x_(%s,%s)", P.vname[P.g.source(e)].c_str(), P.vname[P.g.target(e)].c_str());
        x_e[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        x_e[e].set(GRB_DoubleAttr_Start, P.solution[e]);
    }
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        char name[100];
        sprintf(name, "d_%s", P.vname[v].c_str());
        d_v[v] = model.addVar(0.0, UB, 0.0, GRB_INTEGER, name);
        d_v[v].set(GRB_DoubleAttr_Start, P.node_activation[v]);

        sprintf(name, "l_%s", P.vname[v].c_str());
        l_v[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
    }
    model.update();// needed before using the model inserted variables

    // ILP problem restrictions: -------------------------------------------------
    cout << "Adding the model restrictions:" << endl;

    // There is at most floor(n/2) leaves in a binary tree:
    GRBLinExpr leaves_expr;
    for (DNodeIt v(P.g); v != INVALID; ++v) leaves_expr += l_v[v];
    model.addConstr(leaves_expr <= floor(P.nnodes / 2));
    cout << "-> there is at most floor(n/2) leaves in a binary tree - " << 1 << " constrs" << endl;

    //  A node depth is at least the source distance to its father plus the distance from its father to it:
    int constrCount = 0;
    if (P.source != INVALID)// additional cutting planes
        for (DNodeIt u(P.g); u != INVALID; ++u)
            if (u != P.source)
                for (DNodeIt v(P.g); v != INVALID; ++v)
                    if (v != u && v != P.source) {
                        auto su = P.arc_map[P.source][u], sv = P.arc_map[P.source][v], uv = P.arc_map[u][v];
                        model.addConstr(d_v[v] >=
                                        x_e[uv] * (P.weight[su] + P.weight[uv]) + (1 - x_e[uv]) * P.weight[sv]);
                        constrCount++;
                    }
    cout << "-> a node depth is at least the source distance to its father plus the distance from its father to it - "
         << constrCount << " constrs" << endl;

    // The source is the root and has out-degree one:
    GRBLinExpr s_out_degree_expr;
    for (OutArcIt e(P.g, P.source); e != INVALID; ++e) s_out_degree_expr += x_e[e];
    model.addConstr(s_out_degree_expr == 1);
    GRBLinExpr s_in_degree_expr;
    for (InArcIt e(P.g, P.source); e != INVALID; ++e) s_in_degree_expr += x_e[e];
    model.addConstr(s_in_degree_expr == 0);
    cout << "-> the source is the root and has out-degree one - " << 2 << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        if (v == P.source) continue;

        GRBLinExpr out_degree_expr, in_degree_expr;
        for (OutArcIt e(P.g, v); e != INVALID; ++e) out_degree_expr += x_e[e];
        for (InArcIt e(P.g, v); e != INVALID; ++e) in_degree_expr += x_e[e];

        // The in-degree is one and the out-degree is at most two for each internal node:
        model.addConstr(out_degree_expr <= 2 * (1 - l_v[v]));
        model.addConstr(out_degree_expr >= (1 - l_v[v]));
        model.addConstr(in_degree_expr == 1);

        constrCount += 3;
    }
    cout << "-> the in-degree is one and the out-degree is at most two for each internal node - " << constrCount
         << " constrs" << endl;

    if (P.source != INVALID) {
        model.addConstr(d_v[P.source] == 0);
        constrCount = 1;
    }
    cout << "-> the root is at depth zero - " << constrCount << " constrs" << endl;

    // The depth of the tree is the maximum of the depth of each of its nodes:
    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v, constrCount++) model.addConstr(depth >= d_v[v]);
    cout << "-> the depth of the tree is the maximum of the depth of each of its nodes: - " << constrCount << " constrs"
         << endl;

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v)
        if (v != P.source)
            for (InArcIt e(P.g, v); e != INVALID; ++e) {
                constrCount++;
                if (UB <= P.weight[e]) {// cannot use this edge
                    model.addConstr(x_e[e] == 0);
                    continue;
                }
                model.addConstr(d_v[v] >= d_v[P.g.source(e)] + P.weight[e] + (P.weight[e] + UB) * (x_e[e] - 1) -
                                                  (P.weight[e] + UB) * l_v[P.g.source(e)]);
            }
    cout << "-> a node depth is its parents depth plus the edge to it - " << constrCount << " constrs" << endl;

    // ILP solving: --------------------------------------------------------------
    model.write("gurobi_model.lp");
    model.optimize();// trys to solve optimally within the time limit
    P.stop_counter();

    LB = max(LB, model.get(GRB_DoubleAttr_ObjBound)) * MY_EPS;
    improved |= model.get(GRB_IntAttr_SolCount) > 0;
    if (!improved) return improved;

    // A better solution was found:
    cout << endl << "Nodes depth: ";
    UB = ceil(depth.get(GRB_DoubleAttr_X)) * MY_EPS;
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        int node_depth = ceil(d_v[v].get(GRB_DoubleAttr_X));
        P.node_activation[v] = node_depth;
        cout << P.vname[v].c_str() << '-' << node_depth * MY_EPS << ";";
    }

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)// solved optimally
        LB = UB;

    // Display the best know solution:
    ArcVector to_remove;// save non-original and non-used arcs to remove latter
    cout << endl << "Tree edges: ";
    for (ArcIt e(P.g); e != INVALID; ++e) {
        bool active = x_e[e].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
        P.solution[e] = active;
        if (active) cout << P.vname[P.g.source(e)].c_str() << '-' << P.vname[P.g.target(e)].c_str() << ";";
        else if (!P.original[e])
            to_remove.push_back(e);
    }
    // todo: use lemon's graph windows
    for (auto &e: to_remove) P.g.erase(e);// remove non-original and non-used arcs
    cout << endl;
    cout << "New LB: " << LB << endl;

    return improved;
}

using namespace std::chrono;

int main(int argc, char *argv[]) {
    int time_limit;
    bool only_active_edges = true, tsplib = false;
    string filename, source_node_name;

    if (argc < 3) {
        cout << "Integer Linear Program for the Freeze-Tag Problem using the Gurobi solver;" << endl
             << "Usage:" << endl
             << "\t" << argv[0]
             << " <ftp_graph_filename> <maximum_time_sec> [-tsplib=true|false]"
                " [-only_active_edges=true|false] [LB] [UB]"
             << endl;
        cout << "Example:" << endl << "\t" << argv[0] << " ../instances/tsp_lib/test10.tsp 10 -tsplib=true" << endl;
        exit(0);
    }

    seed = (int) duration_cast<minutes>(system_clock::now().time_since_epoch()).count();
    srand(seed);
    filename = argv[1];
    time_limit = atoi(argv[2]);
    if (argc >= 4) tsplib = strcmp(argv[3], "-tsplib=true") == 0;
    if (argc >= 5) only_active_edges = strcmp(argv[4], "-only_active_edges=true") == 0;
    MY_EPS = 1E-2;
    double LB = 0, UB = MY_INF;// consider MY_INF as infinity.
    if (argc >= 6) LB = atof(argv[5]);
    if (argc >= 7) UB = atof(argv[6]);

    try {
        Problem_Instance P(filename, time_limit, true, tsplib);
        P.print_instance();
        if (solve(P, LB, UB)) {
            cout << "UB cost: " << UB << endl;
            char msg[100];
            sprintf(msg, " Gap of %.2f%%.", 100 * (UB - LB) / UB);
            P.view_solution(LB, UB, msg, only_active_edges);
        } else
            cout << "No solution found." << endl;
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
