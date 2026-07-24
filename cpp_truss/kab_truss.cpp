// kab_truss.cpp — C++ port of code/kab_truss.py (a>0 projection path) + func projection graph.
//
// Computes the (k,a,b)-truss and outputs the surviving vertex set and
// hyperedge-id set, IDENTICAL to Python kab_truss.run(H, k, a, b, network)
// for a>0 (order-independent unique fixpoint).
//
// Projection graph (func.build_projection_graph):
//   node  = an unordered vertex-pair {x,y} that co-occurs in some hyperedge (pid).
//   omega[pid] = #hyperedges containing that pair.
//   edge {p,q} exists iff the two pairs share exactly one vertex (their union is
//     3 vertices x,y,z); sigma(p,q) = #hyperedges containing all of x,y,z.
//   A "triangle" of the projection graph = 3 pairs {xy,xz,yz} of 3 vertices.
//   TN(u,v) = the third pair completing the triangle (or none).
// Algorithm (kab_truss.run, a>0):
//   G = 2k-core of the projection graph.
//   valid triangle {u,v,w}:  sigma>=a  AND  omega[u]+omega[v]+omega[w]-2*sigma >= b
//     (== #he containing all 3 >= a  AND  #he containing >=2 of the 3 >= b).
//   sup[pid] = #valid triangles through pid.  Peel pairs with sup<k (node-truss).
//   Output: nodes = vertices of surviving pairs; hids = hyperedges containing
//   some surviving pair.  All rules monotone => unique fixpoint (matches Python).
//
// Usage: ./kab_truss <network.hyp> <k> <a> <b> <out.txt>   (requires a>=1)
//   out.txt line 1: surviving original vertex ids, sorted
//           line 2: surviving hyperedge ids (1-based, file order), sorted
//
// Memory note: stores omega + per-edge sigma + projection adjacency, but NOT the
// per-pair hyperedge-id SETS that the Python/networkx version keeps (the final
// hid union is recovered by scanning hyperedges against surviving pairs). This
// removes the dominant memory sink that OOMs the Python version.

#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

using namespace std;

static inline uint64_t pk(int a, int b) {           // pack two dense ids, a<b
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: %s <network.hyp> <k> <a> <b> <out.txt>\n", argv[0]);
        return 1;
    }
    const string in_path = argv[1];
    const long K = atol(argv[2]), A = atol(argv[3]), B = atol(argv[4]);
    const string out_path = argv[5];
    if (A < 1) { fprintf(stderr, "this port implements the a>=1 path only\n"); return 2; }

    // ---- load hyperedges (mirror func.Hypergraph.load_from_file) ----
    vector<vector<int>> H;
    {
        ifstream fin(in_path);
        if (!fin) { fprintf(stderr, "cannot open %s\n", in_path.c_str()); return 1; }
        string line;
        while (getline(fin, line)) {
            vector<int> vs; int val = 0; bool have = false, neg = false;
            for (size_t i = 0; i <= line.size(); ++i) {
                char c = (i < line.size()) ? line[i] : ',';
                if (c == '-' && !have) { neg = true; have = true; }
                else if (c >= '0' && c <= '9') { val = val * 10 + (c - '0'); have = true; }
                else { if (have) { vs.push_back(neg ? -val : val); val = 0; have = false; neg = false; } }
            }
            sort(vs.begin(), vs.end());
            vs.erase(unique(vs.begin(), vs.end()), vs.end());
            if (vs.size() > 1) H.push_back(move(vs));
        }
    }
    const int M = (int)H.size();

    // ---- compress vertex ids ----
    unordered_map<int,int> comp; comp.reserve(1 << 16);
    vector<int> id2orig;
    auto cid = [&](int v) { auto it = comp.find(v); if (it != comp.end()) return it->second;
        int n = (int)id2orig.size(); comp.emplace(v, n); id2orig.push_back(v); return n; };
    for (auto& e : H) for (int& v : e) v = cid(v);

    // ---- build pair-nodes (pid), omega ----
    unordered_map<uint64_t,int> pair2id; pair2id.reserve((size_t)M * 4);
    vector<pair<int,int>> pairVerts;      // pid -> (x,y) dense, x<y
    vector<long long> omega;
    auto getpid = [&](int x, int y) -> int {
        uint64_t key = pk(x, y);
        auto it = pair2id.find(key);
        if (it != pair2id.end()) return it->second;
        int id = (int)pairVerts.size();
        pair2id.emplace(key, id);
        pairVerts.push_back({min(x,y), max(x,y)});
        omega.push_back(0);
        return id;
    };
    for (auto& e : H) {
        int s = (int)e.size();
        for (int i = 0; i < s; ++i)
            for (int j = i + 1; j < s; ++j)
                omega[getpid(e[i], e[j])]++;
    }
    const int P = (int)pairVerts.size();

    // ---- projection edges + sigma (via vertex triples) ----
    unordered_map<uint64_t,long long> sigma; sigma.reserve((size_t)P * 2);
    for (auto& e : H) {
        int s = (int)e.size();                       // e sorted ascending (dense ids)
        for (int i = 0; i < s; ++i)
            for (int j = i + 1; j < s; ++j)
                for (int l = j + 1; l < s; ++l) {
                    int pij = pair2id[pk(e[i], e[j])];
                    int pil = pair2id[pk(e[i], e[l])];
                    int pjl = pair2id[pk(e[j], e[l])];
                    sigma[pk(pij, pil)]++;
                    sigma[pk(pij, pjl)]++;
                    sigma[pk(pil, pjl)]++;
                }
    }

    // ---- adjacency from sigma keys ----
    vector<vector<int>> adj(P);
    {
        vector<int> deg(P, 0);
        for (auto& kv : sigma) { deg[(int)(kv.first >> 32)]++; deg[(int)(kv.first & 0xffffffffu)]++; }
        for (int u = 0; u < P; ++u) adj[u].reserve(deg[u]);
        for (auto& kv : sigma) {
            int a = (int)(kv.first >> 32), b = (int)(kv.first & 0xffffffffu);
            adj[a].push_back(b); adj[b].push_back(a);
        }
        for (int u = 0; u < P; ++u) sort(adj[u].begin(), adj[u].end());
    }
    auto hasEdge = [&](int u, int v) {
        return binary_search(adj[u].begin(), adj[u].end(), v);
    };
    // TN: third pair completing the triangle with u,v (share exactly one vertex)
    auto TN = [&](int u, int v) -> int {
        int ux = pairVerts[u].first, uy = pairVerts[u].second;
        int vx = pairVerts[v].first, vy = pairVerts[v].second;
        int ou, ov;                                   // the non-shared endpoints
        // exactly-one-common check
        int common = -1, ncom = 0;
        if (ux == vx || ux == vy) { common = ux; ncom++; }
        if (uy == vx || uy == vy) { common = uy; ncom++; }
        if (ncom != 1) return -1;
        ou = (ux == common) ? uy : ux;
        ov = (vx == common) ? vy : vx;
        auto it = pair2id.find(pk(ou, ov));
        return (it != pair2id.end()) ? it->second : -1;
    };

    vector<char> removed(P, 0);

    // ---- 2k-core of the projection graph ----
    {
        const long thr = 2 * K;
        vector<int> deg(P);
        vector<int> Q;
        for (int u = 0; u < P; ++u) { deg[u] = (int)adj[u].size(); if (deg[u] < thr) { removed[u] = 1; Q.push_back(u); } }
        for (size_t qi = 0; qi < Q.size(); ++qi) {
            int u = Q[qi];
            for (int v : adj[u]) if (!removed[v]) { if (--deg[v] < thr) { removed[v] = 1; Q.push_back(v); } }
        }
    }

    // ---- support: #valid triangles through each pair (on the 2k-core) ----
    vector<long long> sup(P, 0);
    {
        unordered_set<uint64_t> visited; visited.reserve(sigma.size() * 2);
        for (auto& kv : sigma) {
            int u = (int)(kv.first >> 32), v = (int)(kv.first & 0xffffffffu);
            if (removed[u] || removed[v]) continue;
            if (visited.count(kv.first)) continue;
            int w = TN(u, v);
            if (w < 0 || removed[w]) continue;        // (do not mark visited — mirrors Python)
            visited.insert(pk(u, v)); visited.insert(pk(u, w)); visited.insert(pk(v, w));
            long long sig = kv.second;                // sigma(u,v)
            if (sig >= A && (omega[u] + omega[v] + omega[w] - 2 * sig) >= B) {
                sup[u]++; sup[v]++; sup[w]++;
            }
        }
    }

    // ---- peel pairs with sup < k (node-truss on valid triangles) ----
    {
        vector<int> Q;
        for (int u = 0; u < P; ++u) if (!removed[u] && sup[u] < K) Q.push_back(u);
        for (size_t qi = 0; qi < Q.size(); ++qi) {
            int u = Q[qi];
            if (removed[u]) continue;
            for (int v : adj[u]) {
                if (removed[v]) continue;
                int w = TN(u, v);
                if (w < 0 || removed[w]) continue;
                if (!hasEdge(u, w) || !hasEdge(v, w)) continue;
                long long sig = sigma[pk(u, v)];
                if (sig >= A && (omega[u] + omega[v] + omega[w] - 2 * sig) >= B) {
                    sup[v]--;
                    if (sup[v] < K) Q.push_back(v);
                }
            }
            removed[u] = 1;
        }
    }

    // ---- collect result ----
    unordered_set<uint64_t> survPairs; survPairs.reserve(P);
    vector<char> node_seen(id2orig.size(), 0);
    for (int u = 0; u < P; ++u) {
        if (removed[u]) continue;
        int x = pairVerts[u].first, y = pairVerts[u].second;
        survPairs.insert(pk(x, y));
        node_seen[x] = node_seen[y] = 1;
    }
    // hids: a hyperedge is kept iff it contains some surviving pair
    vector<int> hids;
    for (int i = 0; i < M; ++i) {
        auto& e = H[i]; int s = (int)e.size(); bool keep = false;
        for (int p = 0; p < s && !keep; ++p)
            for (int q = p + 1; q < s; ++q)
                if (survPairs.count(pk(e[p], e[q]))) { keep = true; break; }
        if (keep) hids.push_back(i + 1);
    }
    vector<int> nodes;
    for (size_t u = 0; u < node_seen.size(); ++u) if (node_seen[u]) nodes.push_back(id2orig[u]);
    sort(nodes.begin(), nodes.end());
    // hids already ascending by construction

    ofstream fout(out_path);
    for (size_t i = 0; i < nodes.size(); ++i) fout << (i ? " " : "") << nodes[i];
    fout << "\n";
    for (size_t i = 0; i < hids.size(); ++i) fout << (i ? " " : "") << hids[i];
    fout << "\n";
    return 0;
}
