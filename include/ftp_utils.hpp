#ifndef FTP_UTILS_DEFINE
#define FTP_UTILS_DEFINE

#include "mygraphlib.hpp"
#include <chrono>
#include <lemon/min_cost_arborescence.h>
#include <queue>

#define ELAPSED ((chrono::system_clock::now() - P.start).count() / 1E9)
#define _NEW_UB_MESSAGE(SOL, MSG)                                                                                      \
    {                                                                                                                  \
        PrintSolution(P, SOL, MSG);                                                                                    \
        printf("custo: %05.5f - %02.2f%% ótimo\n", UB, 100 * LB / UB);                                                 \
    }
#define NEW_UB_MESSAGE(SOL) _NEW_UB_MESSAGE(SOL, "\nNovo UB.")

using namespace lemon;
using namespace std;

typedef chrono::time_point<chrono::system_clock> time_point;
typedef vector<Node> NodeVector;

// FTP_Instance put all relevant information in one class.
class FTP_Instance {
public:
    FTP_Instance(Graph &graph, NodeStringMap &vvname, NodePosMap &posx, NodePosMap &posy, Node &sourcenode, int &nnodes,
                 int &time_limit);
    void start_counter();

    Graph &g;
    NodeStringMap &vname;
    NodePosMap &px;
    NodePosMap &py;
    const int nnodes;
    Node &source;
    time_point start;
    const int time_limit;
};

void PrintInstanceInfo(FTP_Instance &P);
void PrintSolution(FTP_Instance &P, NodeVector &Sol, const string &msg);

bool ReadFTPGraph(const string &filename, Graph &g, NodeStringMap &vname, NodePosMap &posx, NodePosMap &posy,
                  Node &source, int &nnodes);

bool ViewFTPSolution(FTP_Instance &P, double &LB, double &UB, const string &msg);

#endif// FTP_UTILS_DEFINE
