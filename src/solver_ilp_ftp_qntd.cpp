#include "gurobi_c++.h"
#include "problem_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

const long unsigned seed = 42;// seed to the random number generator

bool solve(Problem_Instance &P, double &LB, double &UB) {
    P.start_counter();

    // Compute the makespan upper bound:
    const int max_D = P.nnodes - 1;                     // maximum possible diameter/degree
    int T_MAX = (int) (2.0 * max_D * ceil(log2(max_D)));// makespan UB (greedy tree based scheduling)
    // T_MAX += max_D;                                     // make room for source reunion
    UB = min(UB, (double) T_MAX);
    T_MAX = UB;

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
        model.set(GRB_IntParam_PumpPasses, 5000);
        model.set(GRB_IntParam_ZeroObjNodes, 500);
        model.set(GRB_DoubleParam_Heuristics, 0.3);
    }

    // ILP problem variables: ----------------------------------------------------
    auto q_t_v = new Graph::NodeMap<GRBVar> *[T_MAX];// quantity of active robots in node v at time t
    auto a_t_v = new Graph::NodeMap<GRBVar> *[T_MAX];// if node v has already been visited at time t
    auto f_t = new GRBVar[T_MAX];                    // if the scheduling is finished by time t

    for (int t = 0; t < T_MAX; t++) {
        char name[100];
        sprintf(name, "f_%d", t);
        f_t[t] = model.addVar(0.0, t != 0, -1.0, GRB_BINARY, name);

        q_t_v[t] = new Graph::NodeMap<GRBVar>(P.g);
        a_t_v[t] = new Graph::NodeMap<GRBVar>(P.g);
        for (NodeIt v(P.g); v != INVALID; ++v) {
            bool is_source = v == P.source;// the source node is always active
            sprintf(name, "q_%d_%s", t, P.vname[v].c_str());
            (*q_t_v[t])[v] = model.addVar(0.0, P.nnodes, 0.0, GRB_INTEGER, name);
            sprintf(name, "a_%d_%s", t, P.vname[v].c_str());
            (*a_t_v[t])[v] = model.addVar(is_source, 1.0, 0.0, GRB_BINARY, name);
        }
    }
    model.update();// run update to use model inserted variables

    // ILP problem restrictions: -------------------------------------------------
    cout << "Adding the model restrictions:" << endl;

    /*model.addConstr((*q_t_v[T_MAX - 1])[P.source] == P.nnodes);
    cout << "->  in the end all robots must get together in source - " << 1 << " constrs" << endl;*/

    int constrCount = 0;
    for (NodeIt v(P.g); v != INVALID; ++v, constrCount += 2) {
        bool is_source = v == P.source;
        model.addConstr((*q_t_v[0])[v] == is_source);
        if (is_source) {
            constrCount--;
            continue;
        }
        model.addConstr((*a_t_v[0])[v] == 0);
    }
    cout << "-> in the beginning there is only one active robot at the source - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (int t = 0; t < T_MAX - 1; t++, constrCount++) {
        model.addConstr(f_t[t + 1] >= f_t[t]);
        for (NodeIt v(P.g); v != INVALID; ++v, constrCount++) {
            if (v == P.source) {
                constrCount--;
                continue;
            }
            model.addConstr((*a_t_v[t + 1])[v] >= (*a_t_v[t])[v]);
        }
    }
    cout << "-> after activated/finished remains activated/finished - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (int t = 0; t < T_MAX - 1; t++)
        for (NodeIt v(P.g); v != INVALID; ++v, constrCount++) {
            if (v == P.source) {
                constrCount--;
                continue;
            }
            GRBLinExpr neighborhood_sum = 0;
            for (IncEdgeIt e(P.g, v); e != INVALID; ++e) {
                auto u = P.g.oppositeNode(v, e);
                neighborhood_sum += (*q_t_v[t])[u];
            }
            model.addConstr((*a_t_v[t + 1])[v] <= (*a_t_v[t])[v] + neighborhood_sum);
        }
    cout << "-> a node becomes active only after receiving an active robot - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (int t = 0; t < T_MAX - 1; t++)
        for (NodeIt v(P.g); v != INVALID; ++v, constrCount++)
            model.addConstr((*q_t_v[t + 1])[v] >= 2 * ((*a_t_v[t + 1])[v] - (*a_t_v[t])[v]));
    cout << "-> at activation there are at least two active robots in a node - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (NodeIt v(P.g); v != INVALID; ++v) {
        if (v == P.source) continue;
        for (int t = 1; t < T_MAX; t++, constrCount++) model.addConstr((*q_t_v[t])[v] <= P.nnodes * (*a_t_v[t])[v]);
    }
    cout << "-> a node can only have active robots after being activated - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (int t = 1; t < T_MAX; t++, constrCount++) {
        GRBLinExpr expr = 0;
        for (NodeIt v(P.g); v != INVALID; ++v) expr += (*a_t_v[t])[v];
        model.addConstr(expr >= P.nnodes * f_t[t]);
    }
    cout << "-> if finished then all robots are active - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (int t = 1; t < T_MAX; t++, constrCount++) {
        GRBLinExpr expr = 0;
        // active robots + one if the node has an inactive one
        for (NodeIt v(P.g); v != INVALID; ++v) expr += (*q_t_v[t])[v] + (1 - (*a_t_v[t])[v]);
        model.addConstr(expr == P.nnodes);
    }
    cout << "-> the total number of robots remains fixed - " << constrCount << " constrs" << endl;

    constrCount = 0;
    for (int t = 0; t < T_MAX - 1; t++)
        for (NodeIt v(P.g); v != INVALID; ++v, constrCount++) {
            GRBLinExpr neighborhood_sum = 0;
            for (IncEdgeIt e(P.g, v); e != INVALID; ++e) {
                auto u = P.g.oppositeNode(v, e);
                neighborhood_sum += (*q_t_v[t])[u];
            }
            GRBLinExpr new_active = (*a_t_v[t + 1])[v] - (*a_t_v[t])[v];// a node gain an active robot after activation
            model.addConstr((*q_t_v[t + 1])[v] <= (*q_t_v[t])[v] + neighborhood_sum + new_active);
        }
    cout << "-> the robots moves according to the nodes neighborhoods - " << constrCount << " constrs" << endl;

    // ILP solving: --------------------------------------------------------------
    model.optimize();// trys to solve optimally within the time limit

    LB = max(LB, model.get(GRB_DoubleAttr_ObjBound) + T_MAX);// the cost constant is added
    bool improved = model.get(GRB_IntAttr_SolCount) > 0;
    if (improved) {                                      // a better solution was found
        UB = model.getObjective().getValue() + T_MAX;    // the cost constant is added
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)// solved optimally
            LB = UB;
    }
    cout << "New LB - " << LB << endl;

    // Display the best know solution:
    cout << endl;
    for (int t = 0; t < T_MAX; t++) {
        bool active = f_t[t].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
        cout << "- t = " << t << " - finished = " << active << endl;
        for (NodeIt v(P.g); v != INVALID; ++v, constrCount++) {
            active = (*a_t_v[t])[v].get(GRB_DoubleAttr_X) >= 1 - MY_EPS;
            cout << "a_" << P.vname[v].c_str() << " = " << active << " ";
        }
        cout << endl;
        for (NodeIt v(P.g); v != INVALID; ++v, constrCount++) {
            int q = (int) (*q_t_v[t])[v].get(GRB_DoubleAttr_X);
            cout << "q_" << P.vname[v].c_str() << " = " << q << " ";
        }
        cout << endl << endl;
    }

    return improved;
}

int main(int argc, char *argv[]) {
    int maxtime;
    Graph g;// graph declaration
    string graph_filename, source_node_name;
    NodeStringMap vname(g); // name of graph nodes
    NodePosMap px(g), py(g);// xy-coordinates for each node
    NodeColorMap vcolor(g); // color of nodes
    EdgeStringMap ename(g); // name for graph edges
    EdgeColorMap ecolor(g); // color of edges
    EdgeValueMap lpvar(g);  // used to obtain the contents of the LP variables
    vector<Node> V;

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
    MY_EPS = 1E-1;
    double LB = 0, UB = MY_INF;// consider MY_INF as infinity.
    if (argc >= 4) LB = atof(argv[3]);
    if (argc >= 5) UB = atof(argv[4]);
    Node source;

    int nnodes;
    if (!ReadProblemGraph(graph_filename, g, vname, px, py, source, nnodes)) {
        cout << "Error while reding the input graph." << endl;
        exit(EXIT_FAILURE);
    }

    Problem_Instance P(g, vname, px, py, source, nnodes, maxtime);
    PrintInstanceInfo(P);

    try {
        if (solve(P, LB, UB)) {
            ViewProblemSolution(P, LB, UB, " Best solution found.", true);
            cout << "cost: " << UB << endl;
        }
    } catch (std::exception &e) {
        cout << "cost: " << UB << endl;
        cerr << "\nException: " << e.what() << endl;
        return 1;
    } catch (GRBException &e) {
        cout << "cost: " << UB << endl;
        cerr << "\nGRBException: " << e.getMessage() << endl;
        return 1;
    }
    return 0;
}