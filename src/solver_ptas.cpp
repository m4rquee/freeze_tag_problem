#include "problem_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

#define EPS 1E-1

bool tree_path(Problem_Instance &P, DNode u, DNode &v, vector<DNode> &path) {
    if (u == v) {
        path.push_back(v);
        return true;
    }
    for (OutArcIt e(P.g, u); e != INVALID; ++e) {
        if (P.weight[e] == 1) {// original edge
            auto aux = P.g.target(e);
            if (tree_path(P, aux, v, path)) {
                path.push_back(u);
                return true;
            }
        }
    }
    return false;
}

bool solve(Problem_Instance &P, double &LB, double &UB) {
    P.start_counter();

    const auto m = 1 / EPS;
    DNode u, v;
    int DIAMETER = 0;
    for (ArcIt e(P.g); e != INVALID; ++e) {
        if (P.weight[e] > DIAMETER) {
            DIAMETER = max(DIAMETER, P.weight[e]);
            u = P.g.source(e);
            v = P.g.target(e);
        }
    }

    vector<DNode> path(P.nnodes);
    tree_path(P, u, v, path);

    P.stop_counter();
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
             << "2+e Polynomial Time Approximation Scheme for the Freeze-Tag Problem;" << endl
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
    MY_EPS = 1;
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
            sprintf(msg, "");
            ViewProblemSolution(P, LB, UB, msg, only_active_edges);
            cout << "UB cost: " << UB << endl;
        }
    } catch (std::exception &e) {
        cout << "UB cost: " << UB << endl;
        cerr << "\nException: " << e.what() << endl;
        return 1;
    }
    return 0;
}
