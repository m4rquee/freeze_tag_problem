#include "ftp_utils.hpp"

FTP_Instance::FTP_Instance(Graph &graph, NodeStringMap &vvname, NodePosMap &posx, NodePosMap &posy, Node &sourcenode,
                           int &nnodes, int &time_limit)
    : g(graph), vname(vvname), px(posx), py(posy), nnodes(nnodes), source(sourcenode), time_limit(time_limit) {}

void FTP_Instance::start_counter() { start = chrono::system_clock::now(); }

void PrintInstanceInfo(FTP_Instance &P) {
    cout << endl << endl;
    cout << "Pickup Delivery Graph Informations" << endl;
    cout << "\tTime limit = " << P.time_limit << "s" << endl;
    cout << "\tSource = " << P.vname[P.source] << endl;
    cout << endl;
}

void PrintSolution(FTP_Instance &P, NodeVector &Sol, const string &msg) {
    // Print the solution to the terminal:
    cout << msg << endl << "\t";
    cout << P.vname[Sol[0]];
    for (int i = 1; i < P.nnodes; i++) cout << "-->" << P.vname[Sol[i]];
    cout << endl;
}

bool ReadFTPGraph(const string &filename, Graph &g, NodeStringMap &vname, NodePosMap &posx, NodePosMap &posy,
                  Node &source, int &nnodes) {
    ReadGraph(filename, g, vname, posx, posy, nullptr);
    nnodes = countNodes(g);
    Node DN[nnodes];
    int i = 0;
    for (NodeIt v(g); v != INVALID; ++v, i++) DN[nnodes - i - 1] = v;
    source = DN[0];
    return true;
}

bool ViewFTPSolution(FTP_Instance &P, double &LB, double &UB, const string &msg) {
    GraphAttributes GA(P.g, P.vname, P.px, P.py);
    GA.SetDefaultNodeAttrib("color=LightGray style=filled width=0.2 height=0.2 fixedsize=true");
    for (EdgeIt e(P.g); e != INVALID; ++e) GA.SetColor(e, "Invis");
    GA.SetColor(P.source, "Red");// source and target are painted in White
    GA.SetShape(P.source, "star");
    GA.SetLabel("Scheduling starting from node " + P.vname[P.source] + " of value " + DoubleToString(UB) +
                ". LB = " + DoubleToString(LB) + ". " + msg);
    GA.View();
    return true;
}
