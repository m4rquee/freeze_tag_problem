#include "ftp_utils.hpp"

FTP_Instance::FTP_Instance(Digraph &graph, DNodeStringMap &vvname, DNodePosMap &posx, DNodePosMap &posy,
                           DNode &sourcenode, int &nnodes, int &time_limit, ArcValueMap &pweight, ArcBoolMap &poriginal)
    : g(graph), vname(vvname), px(posx), py(posy), nnodes(nnodes), source(sourcenode), time_limit(time_limit),
      weight(pweight), original(poriginal) {
    solution = new Arc[nnodes];
}

FTP_Instance::~FTP_Instance() { delete[] solution; }

void FTP_Instance::start_counter() { start = chrono::system_clock::now(); }

void PrintInstanceInfo(FTP_Instance &P) {
    cout << "Freeze-Tag graph information" << endl;
    cout << "\tTime limit = " << P.time_limit << "s" << endl;
    cout << "\tSource = " << P.vname[P.source] << endl;
    cout << endl;
}

bool ReadFTPGraph(const string &filename, Digraph &g, DNodeStringMap &vname, DNodePosMap &posx, DNodePosMap &posy,
                  DNode &source, int &nnodes, ArcValueMap &weight, ArcBoolMap &original, bool calc_clojure) {
    if (!ReadDigraph(filename, g, vname, posx, posy, weight)) return false;
    nnodes = countNodes(g);
    source = GetDNodeByName(g, vname, "0");
    if (calc_clojure) {
        int i = 0;
        vector<Arc> original_edges;
        for (ArcIt e(g); e != INVALID; ++e, i++) {
            original_edges.resize(i + 1);
            original_edges[i] = e;
            original[e] = true;
        }
        for (i--; i >= 0; i--) {// make each edge a reverse arc pair
            auto e = original_edges[i];
            auto rev = g.addArc(g.target(e), g.source(e));
            original[rev] = false;
            weight[rev] = weight[e];
        }

        // Run a Dijkstra using each node as source and them add an arc to each node with the distance:
        DijkstraSolver dijkstra_test(g, weight);
        for (DNodeIt u(g); u != INVALID; ++u) {
            dijkstra_test.run(u);
            for (DNodeIt v(g); v != INVALID; ++v) {
                if (u == v) continue;
                auto aux = findArc(g, u, v);
                if (aux == INVALID) {// if this arc is not present then add it with the distance as weight
                    aux = g.addArc(u, v);
                    weight[aux] = dijkstra_test.dist(v);
                    original[aux] = false;
                }
            }
        }
    }
    return true;
}

bool ViewFTPSolution(FTP_Instance &P, double &LB, double &UB, const string &msg) {
    DigraphAttributes GA(P.g, P.vname, P.px, P.py);
    GA.SetDefaultDNodeAttrib("color=LightGray style=filled width=0.2 height=0.2 fixedsize=true");
    for (ArcIt e(P.g); e != INVALID; ++e) {
        if (!P.original[e]) continue;
        GA.SetColor(e, "#00000070");
        GA.SetAttrib(e, "style=dotted");
    }
    for (int i = 0; i < P.nnodes - 1; i++) {
        auto e = P.solution[i];
        if (P.original[e]) e = P.g.addArc(P.g.source(e), P.g.target(e));// duplicate if already exists
        GA.SetColor(e, "#ff000070");
        GA.SetAttrib(e, "style=dashed splines=true");
    }
    GA.SetColor(P.source, "Red");
    GA.SetLabel("Scheduling starting from node " + P.vname[P.source] + " of value " + DoubleToString(UB) +
                ". LB = " + DoubleToString(LB) + ". " + msg);
    GA.View();
    return true;
}
