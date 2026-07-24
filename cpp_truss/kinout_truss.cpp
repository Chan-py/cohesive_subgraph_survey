// kinout_truss.cpp — C++ port of code/kinout_truss.py (+ func clique-graph peeling).
//
// Computes the (k_in, k_out)-in/out-truss of a hypergraph and outputs the
// surviving vertex set and hyperedge-id set, IDENTICAL to the Python
// kinout_truss.run() result (order-independent unique fixpoint).
//
// Definitions (derived from the Python code, verified against it):
//   Build the clique graph: an edge {u,v} exists iff some (alive) hyperedge
//   contains both; mult(u,v) = #alive hyperedges containing both u and v.
//   For an incidence (e,u,v) with u,v in hyperedge e:
//     sup_in[(e,u,v)]      = |e| - 2                       (triangles inside e)
//     S(u,v)               = sum_{w in N(u)&N(v)} mult(u,w)*mult(v,w)
//     sup_out[(e,u,v)]     = S(u,v) - (|e|-2)
//     sup_node_out[(e,u)]  = max_{v in e\{u}} S(u,v) - (|e|-2)
//   Peel to the maximal sub-hypergraph where every alive hyperedge e has
//   |e| >= k_in+2 AND every incidence (e,u) has sup_node_out[(e,u)] >= k_out.
//   Both rules are monotone => unique fixpoint (matches Python regardless of
//   the intra-round peeling order used there).
//
// Usage: ./kinout_truss <network.hyp> <k_in> <k_out> <out.txt>
//   out.txt line 1: surviving original vertex ids, sorted, space-separated
//           line 2: surviving hyperedge ids (1-based, file order), sorted
//
// Memory: only the clique graph (mult map + adjacency + per-edge S) per round,
// rebuilt each round from the current hyperedges — comparable to the Python
// clique graph but far more compact (no networkx object overhead).

#include <cstdio>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s <network.hyp> <k_in> <k_out> <out.txt>\n", argv[0]);
        return 1;
    }
    const string in_path  = argv[1];
    const long   k_in     = atol(argv[2]);
    const long   k_out    = atol(argv[3]);
    const string out_path = argv[4];

    // ---- load hyperedges (mirror func.Hypergraph.load_from_file) ----
    // each line -> set of ints; skip if <2 distinct; hid = 1-based over kept lines.
    vector<vector<int>> H;                 // H[i] = sorted-unique vertices; hid = i+1
    vector<int> orig_hid;                  // orig_hid[i] = i+1 (kept for clarity)
    {
        ifstream fin(in_path);
        if (!fin) { fprintf(stderr, "cannot open %s\n", in_path.c_str()); return 1; }
        string line;
        while (getline(fin, line)) {
            vector<int> vs;
            // split on commas / whitespace
            int val = 0; bool have = false; bool neg = false;
            for (size_t i = 0; i <= line.size(); ++i) {
                char c = (i < line.size()) ? line[i] : ',';
                if (c == '-' && !have) { neg = true; have = true; }
                else if (c >= '0' && c <= '9') { val = val * 10 + (c - '0'); have = true; }
                else { if (have) { vs.push_back(neg ? -val : val); val = 0; have = false; neg = false; } }
            }
            sort(vs.begin(), vs.end());
            vs.erase(unique(vs.begin(), vs.end()), vs.end());
            if (vs.size() > 1) { H.push_back(move(vs)); orig_hid.push_back((int)H.size()); }
        }
    }
    const int M = (int)H.size();

    // ---- compress vertex ids to dense 0..n-1 ----
    unordered_map<int,int> comp; comp.reserve(1 << 16);
    vector<int> id2orig;
    auto cid = [&](int v) -> int {
        auto it = comp.find(v);
        if (it != comp.end()) return it->second;
        int nid = (int)id2orig.size(); comp.emplace(v, nid); id2orig.push_back(v); return nid;
    };
    for (auto& e : H) for (int& v : e) v = cid(v);   // rewrite H to dense ids
    const int N = (int)id2orig.size();

    vector<char> alive(M, 1);

    auto key = [](int a, int b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
    };

    // ---- peeling rounds ----
    bool changed = true;
    while (changed) {
        changed = false;

        // (1) build clique graph from alive hyperedges: mult(u,v)
        unordered_map<uint64_t,int> mult;
        mult.reserve((size_t)M * 4);
        for (int i = 0; i < M; ++i) {
            if (!alive[i]) continue;
            auto& e = H[i];
            const int s = (int)e.size();
            for (int a = 0; a < s; ++a)
                for (int b = a + 1; b < s; ++b)
                    ++mult[key(e[a], e[b])];
        }

        // adjacency with per-neighbor multiplicity, sorted by neighbor id
        vector<vector<pair<int,int>>> adj(N);   // adj[u] = list of (w, mult(u,w))
        {
            vector<int> deg(N, 0);
            for (auto& kv : mult) {
                int a = (int)(kv.first >> 32), b = (int)(kv.first & 0xffffffffu);
                ++deg[a]; ++deg[b];
            }
            for (int u = 0; u < N; ++u) adj[u].reserve(deg[u]);
            for (auto& kv : mult) {
                int a = (int)(kv.first >> 32), b = (int)(kv.first & 0xffffffffu);
                adj[a].push_back({b, kv.second});
                adj[b].push_back({a, kv.second});
            }
            for (int u = 0; u < N; ++u)
                sort(adj[u].begin(), adj[u].end());
        }

        // (2) S(u,v) for every clique edge = sum over common nbrs w of mult(u,w)*mult(v,w)
        unordered_map<uint64_t,long long> Sval;
        Sval.reserve(mult.size() * 2);
        for (auto& kv : mult) {
            int u = (int)(kv.first >> 32), v = (int)(kv.first & 0xffffffffu);
            // merge-intersect adj[u], adj[v] by neighbor id
            const auto& au = adj[u]; const auto& av = adj[v];
            size_t i = 0, j = 0; long long S = 0;
            while (i < au.size() && j < av.size()) {
                if (au[i].first < av[j].first) ++i;
                else if (au[i].first > av[j].first) ++j;
                else { S += (long long)au[i].second * (long long)av[j].second; ++i; ++j; }
            }
            Sval[kv.first] = S;
        }

        // (3a) in-peel: mark hyperedges with |e| < k_in+2 (snapshot sizes)
        // (3b) out-peel: remove vertex u from e if sup_node_out[(e,u)] < k_out
        //      using snapshot S/sizes; hyperedges killed by in-peel are skipped.
        vector<vector<int>> newH = H;      // will hold vertex removals
        for (int i = 0; i < M; ++i) {
            if (!alive[i]) continue;
            auto& e = H[i];
            const int s = (int)e.size();
            if (s - 2 < k_in) { alive[i] = 0; changed = true; continue; }   // in-peel
            // out-peel: keep vertex u iff max_{v in e\{u}} S(u,v) - (s-2) >= k_out
            vector<int> kept;
            kept.reserve(s);
            for (int a = 0; a < s; ++a) {
                int u = e[a];
                long long bestS = -1;
                for (int b = 0; b < s; ++b) {
                    if (b == a) continue;
                    auto it = Sval.find(key(u, e[b]));
                    long long S = (it != Sval.end()) ? it->second : 0;
                    if (S > bestS) bestS = S;
                }
                long long node_out = bestS - (long long)(s - 2);
                if (node_out >= k_out) kept.push_back(u);
                else changed = true;
            }
            if ((int)kept.size() < s) {
                if (kept.size() < 2) { alive[i] = 0; }   // hyperedge dies (<2 vertices)
                else newH[i] = move(kept);
            }
        }
        H.swap(newH);
    }

    // ---- collect result ----
    vector<char> node_seen(N, 0);
    vector<int> hids;
    for (int i = 0; i < M; ++i) {
        if (!alive[i]) continue;
        if (H[i].size() < 2) continue;
        hids.push_back(orig_hid[i]);
        for (int v : H[i]) node_seen[v] = 1;
    }
    vector<int> nodes;
    for (int u = 0; u < N; ++u) if (node_seen[u]) nodes.push_back(id2orig[u]);
    sort(nodes.begin(), nodes.end());
    sort(hids.begin(), hids.end());

    ofstream fout(out_path);
    for (size_t i = 0; i < nodes.size(); ++i) fout << (i ? " " : "") << nodes[i];
    fout << "\n";
    for (size_t i = 0; i < hids.size(); ++i) fout << (i ? " " : "") << hids[i];
    fout << "\n";
    return 0;
}
