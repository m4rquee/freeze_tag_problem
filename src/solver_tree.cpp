#include "gurobi_c++.h"
#include "problem_utils.hpp"
#include <cstdlib>
#include <iostream>
#include <lemon/list_graph.h>
#include <string>

int seed = 42;// seed to the random number generator

int solve(Problem_Instance &P, int max_degree = 3) {
    int makespan = 0;


    return makespan;
}

int main(int argc, char *argv[]) {
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
    if (argc < 2) {
        cout << endl
             << "Solver for the Freeze-Tag Problem on unweighted trees;" << endl
             << "Usage: " << argv[0] << "  <ftp_graph_filename>" << endl
             << endl;
        cout << "Example:" << endl << "\t" << argv[0] << " ../instances/small/tree-deg3-n7-1.g" << endl << endl;
        exit(0);
    }

    using namespace std::chrono;
    seed = (int) duration_cast<minutes>(system_clock::now().time_since_epoch()).count();
    srand(seed);
    graph_filename = argv[1];

    DNode source;
    int nnodes, source_radius;
    if (!ReadProblemGraph(graph_filename, g, vname, px, py, source, nnodes, weight, original, source_radius, true)) {
        cout << "Error while reding the input graph." << endl;
        exit(EXIT_FAILURE);
    }

    Problem_Instance P(g, vname, px, py, source, nnodes, -1, weight, original, source_radius);
    PrintInstanceInfo(P);

    try {
        int makespan = solve(P);
        ViewProblemSolution(P, makespan, makespan, "", false);
        cout << "makespan cost: " << makespan << endl;
    } catch (std::exception &e) {
        cerr << "\nException: " << e.what() << endl;
        return 1;
    } catch (GRBException &e) {
        cerr << "\nGRBException: " << e.getMessage() << endl;
        return 1;
    }
    return 0;
}
