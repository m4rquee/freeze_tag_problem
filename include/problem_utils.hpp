#ifndef PROBLEM_UTILS_DEFINE
#define PROBLEM_UTILS_DEFINE

#include "mygraphlib.hpp"
#include <chrono>
#include <lemon/dijkstra.h>

using namespace lemon;

typedef chrono::time_point<chrono::system_clock> time_point;
typedef Dijkstra<Digraph, ArcIntMap> DijkstraSolver;

class Problem_Instance {// Problem_Instance has all relevant information in one class
private:
    bool read_instance(const string &filename, bool calc_clojure, bool tsplib);

public:
    Problem_Instance(const string &filename, int time_limit, bool calc_clojure = false, bool tsplib = false);

    // Auxiliary function:
    void print_instance();
    bool view_solution(double LB, double UB, const string &msg, bool only_active_edges);
    void start_counter();
    void stop_counter();

    // Time keeping:
    time_point start;
    time_point stop;
    unsigned int time_limit;

    // Graph instance:
    Digraph g;           // underlying digraph
    DNodeStringMap vname;// node name
    ArcIntMap weight;    // arc weight
    DNodePosMap px;      // node x position
    DNodePosMap py;      // node y position
    unsigned int nnodes; // number of nodes
    unsigned int narcs;  // number of arcs
    DNode source;        // Freeze-Tag's source node

    // Auxiliary structures:
    ArcBoolMap original;// weather some arc was added in the clojure calculation
    ArcBoolMap solution;// if the arc is present in the found solution
    int radius;         // the source radius if there is one, or the graph's radius
    // todo: use the AdjacencyMatrix class
    map<DNode, map<DNode, Arc>> arc_map;// graph's adjacency matrix
    DNodeIntMap node_activation;        // node activation time
    int solution_makespan;              // makespan of the found solution
};

#endif// PROBLEM_UTILS_DEFINE
