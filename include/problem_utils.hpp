#ifndef PROBLEM_UTILS_DEFINE
#define PROBLEM_UTILS_DEFINE

#include "mygraphlib.hpp"

#include "gurobi_c++.h"
#include <chrono>
#include <lemon/dijkstra.h>

#define GRB_CB_START -1

using namespace lemon;

typedef chrono::time_point<chrono::system_clock> time_point;
typedef Dijkstra<Digraph, ArcIntMap> DijkstraSolver;

class Problem_Instance {// Problem_Instance has all relevant information in one class
private:
    void read_instance(const string &filename, bool tsplib);
    void read_tsplib_instance(const string &filename);

public:
    Problem_Instance(const string &filename, int time_limit, bool calc_clojure = false, bool tsplib = false);

    // Auxiliary function:
    void print_instance();
    void view_solution(double LB, double UB, const string &msg, bool only_active_edges);
    void start_counter();
    void stop_counter();
    void clojure();

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
    int source_radius;  // the source radius
    int graph_radius;   // the graph's radius
    int graph_diameter; // the graph's diameter
    // todo: use the AdjacencyMatrix class
    map<DNode, map<DNode, Arc>> arc_map;// graph's adjacency matrix
    DNodeIntMap node_activation;        // node activation time
    int solution_makespan = -1;         // makespan of the found solution
    bool tsplib;                        // the input follows the tsplib pattern
    bool complete;                      // the graph will be complete
};

double greedy_solution(Problem_Instance &P, int max_degree = 3);

class LocalSearchCB : public GRBCallback {
    Problem_Instance &P;

    Digraph::ArcMap<GRBVar> &x_e;

    ArcBoolMap arc_value;
    DNodeValueMap node_depth;
    DNodeIntMap node_degree;

    double (GRBCallback::*solution_value)(GRBVar) = nullptr;
    void (LocalSearchCB::*set_solution)(GRBVar, double) = nullptr;

    typedef pair<int, int> ii_pair;         // int, int pair
    typedef pair<ii_pair, DNode> iin_triple;// int, int, node triple
    typedef pair<int, DNode> in_pair;       // int, node pair

    vector<in_pair> node_depth_heap;     // max heap ordered by node's depth
    vector<iin_triple> available_fathers;// min heap ordered by node's depth and then their degree

public:
    LocalSearchCB(Problem_Instance &_P, Digraph::ArcMap<GRBVar> &_x_e);
    double init();

protected:
    inline void setSolutionStart(GRBVar v, double val) {}

    inline DNode get_father(const DNode &v) {
        for (InArcIt e(P.g, v); e != INVALID; ++e)
            if (arc_value[e]) return P.g.source(e);
        return INVALID;
    }

    bool is_ancestor(const DNode &u, const DNode &v);
    bool improving_swap(const DNode &child, const DNode &new_father);
    void calc_depth(DNode &v, double v_depth);
    void heap_init();

    void callback() override;
};

#endif// PROBLEM_UTILS_DEFINE
