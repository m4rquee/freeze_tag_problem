#ifndef Problem_UTILS_DEFINE
#define Problem_UTILS_DEFINE

#include "mygraphlib.hpp"
#include <chrono>
#include <lemon/dijkstra.h>

using namespace lemon;
using namespace std;

typedef chrono::time_point<chrono::system_clock> time_point;
typedef vector<Arc> ArcVector;
typedef vector<DNode> DNodeVector;
typedef Dijkstra<Digraph, ArcIntMap> DijkstraSolver;

// Problem_Instance put all relevant information in one class.
class Problem_Instance {
public:
    Problem_Instance(Digraph &graph, DNodeStringMap &vvname, DNodePosMap &posx, DNodePosMap &posy, DNode &sourcenode,
                     int nnodes, int time_limit, ArcIntMap &weight, ArcBoolMap &original, int source_radius);
    void start_counter();
    void stop_counter();

    Digraph &g;
    DNodeStringMap &vname;
    DNodePosMap &px;
    DNodePosMap &py;
    const int nnodes;
    DNode &source;
    time_point start;
    time_point stop;
    const int time_limit;
    ArcIntMap &weight;
    ArcBoolMap &original;
    ArcBoolMap solution;
    int source_radius;                  // the source radius if there is one, or the graph`s radius
    map<DNode, map<DNode, Arc>> arc_map;// used for fast arc lookup
    DNodeIntMap node_height;
    double solution_height;
};

void PrintInstanceInfo(Problem_Instance &P);

bool ReadProblemGraph(const string &filename, Digraph &g, DNodeStringMap &vname, DNodePosMap &posx, DNodePosMap &posy,
                      DNode &source, int &nnodes, ArcIntMap &weight, ArcBoolMap &original, int &source_radius,
                      bool calc_clojure = false, bool tsplib = false);

bool ViewProblemSolution(Problem_Instance &P, double LB, double UB, const string &msg, bool only_active_edges);

#endif// Problem_UTILS_DEFINE
