#include "ftp_utils.hpp"
#include "gurobi_c++.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

const long unsigned seed = 42;// seed to the random number generator

bool solve(FTP_Instance &P, double &LB, double &UB) {
    P.start_counter();

    // Calculates the best know objective bounds:
    LB = max(LB, P.source_radius);
    auto MAX_EDGE = 0.0;
    for (ArcIt e(P.g); e != INVALID; ++e)
        if (P.original[e]) MAX_EDGE = max(MAX_EDGE, P.weight[e]);
    MY_INF = P.nnodes * MAX_EDGE;
    UB = MAX_EDGE * min(UB, ceil(2 * MAX_EDGE * log2(P.nnodes)));
    cout << "Set parameter MAX_EDGE to value " << MAX_EDGE << endl;

    // Gurobi ILP problem setup:
    auto *env = new GRBEnv();
    env->set(GRB_IntParam_Seed, seed);
    env->set(GRB_DoubleParam_TimeLimit, P.time_limit);
    env->set(GRB_DoubleParam_Cutoff, UB);// set the best know UB
    GRBModel model = GRBModel(*env);
    model.set(GRB_StringAttr_ModelName, "Freeze-Tag Problem");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // ILP solver parameters: ----------------------------------------------------
    if (P.nnodes >= 20) {// focus only on new UBs
        model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
        model.set(GRB_IntParam_Cuts, GRB_CUTS_AGGRESSIVE);
        model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
        model.set(GRB_IntParam_MinRelNodes, 500);
        model.set(GRB_IntParam_PumpPasses, 10);
        model.set(GRB_IntParam_ZeroObjNodes, 500);
        model.set(GRB_DoubleParam_Heuristics, 0.01);
    }

    // ILP problem variables: ----------------------------------------------------
    Digraph::ArcMap<GRBVar> x_e(P.g);                                       // if arc e is present in the wake-up tree
    Digraph::NodeMap<GRBVar> t_v(P.g);                                      // activation time of node v
    auto makespan = model.addVar(LB, UB, MAX_EDGE, GRB_INTEGER, "makespan");// wake-up tree's makespan

    for (ArcIt e(P.g); e != INVALID; ++e) {
        char name[100];
        sprintf(name, "x_(%s,%s)", P.vname[P.g.source(e)].c_str(), P.vname[P.g.target(e)].c_str());
        x_e[e] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY, name);
    }
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        char name[100];
        sprintf(name, "t_%s", P.vname[v].c_str());
        t_v[v] = model.addVar(0.0, UB, 0.0, GRB_INTEGER, name);
    }
    model.update();// run update to use model inserted variables

    // ILP problem restrictions: -------------------------------------------------
    cout << "Adding the model restrictions:" << endl;

    // The source is the root and has out-degree one:
    GRBLinExpr s_out_degree_expr;
    for (OutArcIt e(P.g, P.source); e != INVALID; ++e) s_out_degree_expr += x_e[e];
    model.addConstr(s_out_degree_expr == 1);
    GRBLinExpr s_in_degree_expr;
    for (InArcIt e(P.g, P.source); e != INVALID; ++e) s_in_degree_expr += x_e[e];
    model.addConstr(s_in_degree_expr == 0);
    cout << "-> the source is the root and has out-degree one - " << 2 << " constrs" << endl;

    int constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        if (v == P.source) continue;
        GRBLinExpr out_degree_expr, in_degree_expr;
        for (OutArcIt e(P.g, v); e != INVALID; ++e) out_degree_expr += x_e[e];
        for (InArcIt e(P.g, v); e != INVALID; ++e) in_degree_expr += x_e[e];
        // The in-degree is one and the out-degree is at most two for each internal node:
        model.addConstr(out_degree_expr <= 2);
        model.addConstr(in_degree_expr == 1);
        constrCount += 2;
    }
    cout << "-> the in-degree is one and the out-degree is at most two for each internal node - " << constrCount
         << " constrs" << endl;

    GRBLinExpr arc_sum_expr;
    for (ArcIt e(P.g); e != INVALID; ++e) arc_sum_expr += x_e[e];
    model.addConstr(arc_sum_expr == P.nnodes - 1);
    cout << "-> the number of arcs is n-1 for any tree - " << 1 << " constrs" << endl;

    model.addConstr(t_v[P.source] == 0);
    cout << "-> the root begins activated - " << 1 << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v, constrCount++) model.addConstr(makespan >= t_v[v]);
    cout << "-> makespan is the last node is activation time - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (DNodeIt v(P.g); v != INVALID; ++v)
        for (InArcIt e(P.g, v); e != INVALID; ++e, constrCount++)
            model.addConstr(t_v[v] >= t_v[P.g.source(e)] + P.weight[e] + MY_INF * (x_e[e] - 1));
    cout << "-> a node is activated at its parent time plus the time to get to it - " << constrCount << " constrs"
         << endl;

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


    // ILP solving: --------------------------------------------------------------
    model.optimize();// trys to solve optimally within the time limit

    LB = max(LB, floor(model.get(GRB_DoubleAttr_ObjBound) / MAX_EDGE));
    bool improved = model.get(GRB_IntAttr_SolCount) > 0;
    if (improved) {// a better solution was found
        UB = floor(model.get(GRB_DoubleAttr_ObjVal) / MAX_EDGE);
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)// solved optimally
            LB = UB;
    }
    cout << "New LB - " << LB << endl;

    // Display the best know solution:
    int i = 0;
    cout << endl << "Tree edges: ";
    for (ArcIt e(P.g); e != INVALID; ++e) {
        bool active = x_e[e].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
        if (active) {
            cout << P.vname[P.g.source(e)].c_str() << '-' << P.vname[P.g.target(e)].c_str() << ";";
            P.solution[i++] = e;
        }
        if (!P.original[e] && !active) P.g.erase(e);
    }
    cout << endl << "Nodes wakeup time: ";
    for (DNodeIt v(P.g); v != INVALID; ++v) {
        int activation_time = ceil(t_v[v].get(GRB_DoubleAttr_X));
        cout << P.vname[v].c_str() << '-' << activation_time << ";";
    }
    cout << endl;

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
    ArcValueMap weight(g);   // arc weights
    ArcBoolMap original(g);  // if an arc is original
    vector<DNode> V;

    set_pdfreader("xdg-open");// the Linux will choose the default one
    if (argc < 3) {
        cout << endl
             << "Integer Linear Program for the Freeze-Tag Problem using the Gurobi solver;" << endl
             << "Usage: " << argv[0] << "  <ftp_graph_filename> <maximum_time_sec>" << endl
             << endl;
        cout << "Example:" << endl
             << "\t" << argv[0] << " "
             << "../instances/small/tree-003.g 10" << endl
             << endl;
        exit(0);
    }

    graph_filename = argv[1];
    maxtime = atoi(argv[2]);
    if (argc >= 4) tsplib = atoi(argv[3]);
    if (argc >= 5) only_active_edges = atoi(argv[4]);
    MY_EPS = 1E-1;
    double LB = 0, UB = MY_INF;// consider MY_INF as infinity.
    if (argc >= 6) LB = atof(argv[5]);
    if (argc >= 7) UB = atof(argv[6]);
    DNode source;

    int nnodes;
    double source_radius;
    if (!ReadFTPGraph(graph_filename, g, vname, px, py, source, nnodes, weight, original, source_radius, true, tsplib)) {
        cout << "Error while reding the input graph." << endl;
        exit(EXIT_FAILURE);
    }

    FTP_Instance P(g, vname, px, py, source, nnodes, maxtime, weight, original, source_radius);
    PrintInstanceInfo(P);

    try {
        if (solve(P, LB, UB)) {
            ViewFTPSolution(P, LB, UB, " Best solution found.", only_active_edges);
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
