#include "gurobi_c++.h"
#include "problem_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

int seed = 42;// seed to the random number generator

bool solve(Problem_Instance &P, double &LB, double &UB) {
    P.start_counter();

    // Calculates the best known objective bounds: -------------------------------
    double minimum_depth = ceil(log(P.nnodes) / log(2));
    LB = max(LB / MY_EPS, (double) P.source_radius);
    auto auxUB = P.graph_diameter * minimum_depth;
    bool improved = auxUB < UB / MY_EPS;
    UB = max(LB, min(UB / MY_EPS, auxUB));

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
    Digraph::ArcMap<GRBVar> y_e(P.g);              // if arc e is present in the solution tree
    Digraph::ArcMap<GRBVar> x_e(P.g);              // number of nodes in the sutree rooted at j (e=ij)
    Digraph::ArcMap<map<DNode, GRBVar>> v_e_k(P.g);// if e is on the path from root to k
    GRBVar depth;                                  // the solution tree depth

    // Cutoff and bounds setup: --------------------------------------------------
    cout << "Set parameter LB to value " << LB * MY_EPS << endl;
    cout << "Set parameter UB to value " << UB * MY_EPS << endl;
    auto cutoff = UB + 1;
    model.set(GRB_DoubleParam_Cutoff, cutoff);// set the best know UB

    // ILP problem variables startup: --------------------------------------------
    depth = model.addVar(LB, UB, 1.0, GRB_INTEGER, "depth");

    for (ArcIt e(P.g); e != INVALID; ++e) {
        char name[100];
        auto u_name = P.vname[P.g.source(e)].c_str(), v_name = P.vname[P.g.target(e)].c_str();
        sprintf(name, "y_(%s,%s)", u_name, v_name);
        y_e[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);

        sprintf(name, "x_(%s,%s)", u_name, v_name);
        x_e[e] = model.addVar(0.0, P.nnodes, 0.0, GRB_INTEGER, name);

        for (DNodeIt k(P.g); k != INVALID; ++k) {
            sprintf(name, "v_(%s,%s)^%s", u_name, v_name, P.vname[k].c_str());
            v_e_k[e][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        }
    }
    model.update();// needed before using the model inserted variables

    // ILP problem restrictions: -------------------------------------------------
    cout << "Adding the model restrictions:" << endl;

    // The number of nodes in a subtree rooted at some node equals the number of nodes in each subtree rooted at its
    // children plus one:
    int constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        if (v == P.source) continue;
        GRBLinExpr out_expr, in_expr;
        for (OutArcIt e(P.g, v); e != INVALID; ++e) out_expr += x_e[e];
        for (InArcIt e(P.g, v); e != INVALID; ++e) in_expr += x_e[e];
        // The in-degree is one and the out-degree is at most two for each internal node:
        model.addConstr(out_expr + 1 == in_expr);
        constrCount++;
    }
    cout << "-> the number of nodes in a subtree rooted at some node equals the number of nodes in each subtree rooted "
            "at its children plus one - "
         << constrCount << " constrs" << endl;

    // Every node is in the subtree rooted at the source:
    GRBLinExpr root_out_expr;
    for (OutArcIt e(P.g, P.source); e != INVALID; ++e) root_out_expr += x_e[e];
    model.addConstr(root_out_expr == P.nnodes - 1);
    cout << "-> every node is in the subtree rooted at the source: - " << 1 << " constrs" << endl;

    // Only active arcs can have subtree rooted at its head:
    constrCount = 0;
    for (ArcIt e(P.g); e != INVALID; ++e, constrCount += 2) {
        model.addConstr(x_e[e] >= y_e[e]);
        model.addConstr(x_e[e] <= (P.nnodes - 1) * y_e[e]);
    }
    cout << "-> only active arcs can have subtree rooted at its head - " << constrCount << " constrs" << endl;

    // The out-degree of each node is at most two if it is internal, otherwise, it's at most one:
    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v, constrCount++) {
        GRBLinExpr out_degree_expr;
        for (OutArcIt e(P.g, v); e != INVALID; ++e) out_degree_expr += y_e[e];
        model.addConstr(out_degree_expr <= 1 + (v != P.source));
    }
    cout << "-> The out-degree of each node is at most two if it is internal, otherwise, it's at most one "
         << constrCount << " constrs" << endl;

    // The number of edges is n-1 for any tree:
    GRBLinExpr arc_sum_expr;
    for (ArcIt e(P.g); e != INVALID; ++e) arc_sum_expr += x_e[e];
    model.addConstr(arc_sum_expr == P.nnodes - 1);
    cout << "-> the number of edges is n-1 for any tree - " << 1 << " constrs" << endl;

    // If edge ij is present, then it is on the path from the root to j:
    constrCount = 0;
    for (ArcIt e(P.g); e != INVALID; ++e, constrCount++) {
        auto v = P.g.target(e);
        model.addConstr(v_e_k[e][v] == y_e[e]);
    }
    cout << "-> if edge ij is present, then it is on the path from the root to j - " << constrCount << " constrs"
         << endl;

    // The number of root-to-node paths traversing an edge equals the number of nodes in the subtree rooted in said
    // arc's head:
    constrCount = 0;
    for (ArcIt e(P.g); e != INVALID; ++e, constrCount++) {
        GRBLinExpr arc_usage_expr;
        for (DNodeIt k(P.g); k != INVALID; ++k) {
            if (k == P.source) continue;
            arc_usage_expr += v_e_k[e][k];
        }
        model.addConstr(arc_usage_expr == x_e[e]);
    }
    cout << "-> the number of root-to-node paths traversing an edge equals the number of nodes in the subtree rooted "
            "in said arc's head - "
         << constrCount << " constrs" << endl;

    // The number of root-to-node paths going through a node does not change if that node is not the target:
    constrCount = 0;
    for (DNodeIt k(P.g); k != INVALID; ++k) {
        if (k == P.source) continue;
        for (DNodeIt i(P.g); i != INVALID; ++i) {
            if (i == P.source || i == k) continue;
            GRBLinExpr out_paths_expr, in_paths_expr;

            for (OutArcIt e(P.g, i); e != INVALID; ++e) out_paths_expr += v_e_k[e][k];
            for (InArcIt e(P.g, i); e != INVALID; ++e) in_paths_expr += v_e_k[e][k];

            model.addConstr(out_paths_expr == in_paths_expr);
            constrCount++;
        }
    }
    cout << "-> the number of root-to-node paths going through a node does not change if that node is not the target - "
         << constrCount << " constrs" << endl;

    // The depth of the tree is the maximum of the depth of each of its nodes:
    constrCount = 0;
    for (DNodeIt k(P.g); k != INVALID; ++k) {
        if (k == P.source) continue;

        GRBLinExpr path_cost_expr;
        for (ArcIt e(P.g); e != INVALID; ++e)
            // The weight of the path from root-to-k is the sum of the arc weights of the arcs used to reach k:
            path_cost_expr += P.weight[e] * v_e_k[e][k];

        model.addConstr(depth >= path_cost_expr);
        constrCount++;
    }
    cout << "-> the depth of the tree is the maximum of the depth of each of its nodes: - " << constrCount << " constrs"
         << endl;

    // Only one of a parallel arc pair is allowed:
    constrCount = 0;
    for (DNodeIt u(P.g); u != INVALID; ++u)// additional cutting planes
        for (OutArcIt e(P.g, u); e != INVALID; ++e) {
            DNode v = P.g.target(e);
            if (Digraph::id(u) < Digraph::id(v)) {// avoid adding duplicates
                constrCount++;
                model.addConstr(y_e[e] + y_e[findArc(P.g, v, u)] <= 1);
            }
        }
    cout << "-> only one of a parallel arc pair is allowed - " << constrCount << " constrs" << endl;

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

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)// solved optimally
        LB = UB;

    // Display the best know solution:
    ArcVector to_remove;// save non-original and non-used arcs to remove latter
    cout << endl << "Tree edges: ";
    for (ArcIt e(P.g); e != INVALID; ++e) {
        bool active = y_e[e].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
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
