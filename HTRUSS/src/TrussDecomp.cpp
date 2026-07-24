// TrussDecomp.cpp
//
// Computes, for every ORIGINAL hyperedge, its hyper-truss "trussness" value and
// writes a CSV with columns: edge_id,trussness
//
// This is a faithful C++ port of the Python reference implementation in
//   code/hyper_k_trussness.py  (functions get_support() + run())
// The triangle enumeration is the same hyper-triangle definition used by the
// published HTRUSS code (code/HTRUSS/src/Baseline.cpp): for edge a, gather
// worklist of b<a sharing a node; (a,b,c) is a triangle when source[b],source[c]
// are disjoint AND b,c intersect. We enumerate over the SAME reduced graph that
// the Python builds (uniqueEdges + core reduction), so the multiset of triangles
// and hence htn/adj are identical.
//
// At each found triangle (a,b,c) with per-edge weights w_a,w_b,w_c we accumulate
// (keyed by representative ORIGINAL edge ids), exactly as the Python get_support:
//   htn[a] += w_b*w_c;  htn[b] += w_a*w_c;  htn[c] += w_a*w_b;
//   adj[a][b]+=w_c; adj[b][a]+=w_c; adj[a][c]+=w_b; adj[c][a]+=w_b;
//   adj[b][c]+=w_a; adj[c][b]+=w_a;
// Then we peel exactly as Python run().
//
// Edge-id semantics: original edge ids are 1-based and skip input lines with < 2
// distinct nodes, matching func.Hypergraph.add_hyperedge (id = count+1).
// htn/adj are keyed by ONE representative original id per group of merged
// duplicate edges (matching Python's new_old_map). Duplicates that were merged
// away and edges removed by the reduction are "not in htn" -> trussness 0,
// exactly as the Python does.
//
// CLI:  ./TrussDecomp.exe <path-to-network.hyp> <output.csv>

#include <bits/stdc++.h>
#include "hypergraph.hpp"
using namespace std;

// ---------------------------------------------------------------------------
// Replicate CPython set iteration order for a set of nonnegative ints.
// hash(i)==i for nonneg ints. CPython set uses open addressing with a linear
// probe run of LINEAR_PROBES then a perturbed jump; iteration yields entries in
// table-slot order. The Python reference's uniqueEdges picks, among the edges
// sharing a node-set key, the LAST one in this set-iteration order to be the
// group's representative (new2old_eid[ne] is overwritten per element). We
// reproduce the table layout so the representative matches EXACTLY.
// ---------------------------------------------------------------------------
namespace pyset {
  constexpr size_t LINEAR_PROBES = 9;
  constexpr size_t PERTURB_SHIFT = 5;
  struct Entry { long long key; bool used; };
  struct PySet {
    vector<Entry> table; size_t mask, fill, used;
    PySet() { init(8); }
    void init(size_t sz){ table.assign(sz, Entry{0,false}); mask=sz-1; fill=0; used=0; }
    void add(long long key){
      if (insert(key)) {
        if (fill*5 >= (mask+1)*3) {
          size_t target = (used>50000? used*2 : used*4);
          size_t ns=8; while (ns<=target) ns<<=1;
          rehash(ns);
        }
      }
    }
    bool insert(long long key){
      size_t perturb=(size_t)key;
      size_t i=(size_t)key & mask;
      while (true){
        for (size_t j=0;j<=LINEAR_PROBES;j++){
          size_t idx=(i+j)&mask;
          Entry &e=table[idx];
          if (!e.used){ e.used=true; e.key=key; used++; fill++; return true; }
          if (e.key==key) return false;
        }
        perturb>>=PERTURB_SHIFT;
        i=(i*5+1+perturb)&mask;
      }
    }
    void rehash(size_t ns){
      vector<Entry> old=table; init(ns);
      for (auto &e:old) if (e.used) insert(e.key);
    }
    vector<long long> iter() const {
      vector<long long> out;
      for (auto &e:table) if (e.used) out.push_back(e.key);
      return out;
    }
  };
  // Return the LAST element in CPython set-iteration order of the given keys,
  // inserted in the order provided (matches key2hids[key].add(hid) scan order).
  inline long long lastInIter(const vector<long long>& keysInOrder){
    PySet s; for (long long k:keysInOrder) s.add(k);
    auto it=s.iter(); return it.empty()? (keysInOrder.empty()?-1:keysInOrder.back()) : it.back();
  }
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << "Usage: " << argv[0] << " <path-to-network.hyp> <output.csv>" << endl;
    return -1;
  }
  const string inPath = argv[1];
  const string outPath = argv[2];

  auto time_start = chrono::steady_clock::now();

  // -------------------------------------------------------------------------
  // Read the ORIGINAL file EXACTLY like func.Hypergraph.load_from_file:
  //   split each line on non-digits, take the SET of ints; keep only lines with
  //   > 1 distinct node; assign 1-based ids in file order.
  // -------------------------------------------------------------------------
  vector<vector<vid_t>> origKey; // origKey[i] = sorted node set of original edge (i+1)
  {
    ifstream fin(inPath);
    if (!fin.is_open()) { cerr << "Cannot open " << inPath << endl; return -1; }
    string line;
    while (getline(fin, line)) {
      vector<vid_t> nodes;
      long long v = 0; bool inNum = false;
      for (char c : line) {
        if (c >= '0' && c <= '9') { v = v * 10 + (c - '0'); inNum = true; }
        else { if (inNum) { nodes.push_back((vid_t)v); v = 0; inNum = false; } }
      }
      if (inNum) nodes.push_back((vid_t)v);
      sort(nodes.begin(), nodes.end());
      nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
      if (nodes.size() > 1) origKey.push_back(std::move(nodes));
    }
  }
  const long long E_orig = (long long)origKey.size(); // ids 1..E_orig

  // number of ORIGINAL distinct nodes (== len(G.nodes) in Python run())
  long long nNodes;
  {
    unordered_set<vid_t> vs;
    for (auto &nk : origKey) for (vid_t v : nk) vs.insert(v);
    nNodes = (long long)vs.size();
  }

  // -------------------------------------------------------------------------
  // Literal port of the Python reduction (get_support prefix):
  //   while True:
  //     uniqueEdges (merge duplicate node-sets, sum weights, drop <2-node edges)
  //     core        (remove nodes with weighted-deg<2; cascade on unweighted
  //                  incident-edge-count<2; edges are NEVER deleted here, only
  //                  shrunk by node removal)
  //     stop when node-count and edge-count both stop changing.
  // We carry a representative ORIGINAL id per current edge through the whole
  // process (mirroring old2new_eid/new2old_eid -> new_old_map). The chain keeps
  // ONE representative original id per merged group.
  //
  // Data structures mirror func.Hypergraph:
  //   hyperedges: eid -> set<node>       (we use sorted vector as the node-set)
  //   nodes:      nid -> set<eid>        (incident edges)
  //   weight:     eid -> w
  //   rep:        eid -> representative original id
  // eids are dense indices into vectors; "removed" edges/nodes tracked by flag.
  // -------------------------------------------------------------------------

  // Current graph edges (after each uniqueEdges rebuild these are re-indexed 0..).
  vector<vector<vid_t>> H;    // node-set per edge (sorted)
  vector<long long>     W;    // weight per edge
  vector<long long>     REP;  // representative original id per edge

  // uniqueEdges: LITERAL port of Python uniqueEdges.
  //   For each CURRENT edge with >=2 nodes, group by node-set key.
  //   new weight[key] = NUMBER OF CURRENT EDGES with that key  (norm2weight is
  //   incremented by +1 per edge -- it does NOT sum prior weights, so weights
  //   are RECOMPUTED every round, not accumulated).
  //   representative: new2old_eid chains to ONE original id per key; on a merge
  //   the LAST edge processed wins. (Only matters for identical duplicate edges;
  //   see report caveat -- exact winner among true duplicates is Python set-hash
  //   dependent, but the trussness VALUE it carries is well-defined.)
  // firstCall: on the very first uniqueEdges, key2hids stores ORIGINAL ids and
  // new2old_eid is empty, so the group members for the pyset are original ids
  // (== REP[e], since REP is initialised to the original id). On later calls,
  // key2hids stores CURRENT edge ids (== H index e), and new2old_eid[old]=REP[e].
  // Either way, the pyset must be built over the SAME integer values Python puts
  // in key2hids[key], in ascending scan order. On the first call those values are
  // the original ids; on later calls they are the current H indices (0-based).
  bool firstCall = true;
  auto uniqueEdges = [&]() {
    // Group current edges (>=2 nodes) by node-set key. Keep, per key:
    //   - new weight = count of current edges
    //   - the list of member "python ids" in ascending scan order (for pyset)
    //   - a map from that python id -> REP (current representative original id)
    map<vector<vid_t>, size_t> idx;
    vector<vector<vid_t>> H2;
    vector<long long> W2;
    vector<vector<long long>> memberPyIds;   // per new edge: python ids in scan order
    vector<vector<long long>> memberReps;     // parallel: REP of each member
    for (size_t e = 0; e < H.size(); e++) {
      if (H[e].size() < 2) continue;
      long long pyid = firstCall ? REP[e] : (long long)e; // key2hids stores this
      auto it = idx.find(H[e]);
      if (it == idx.end()) {
        size_t j = H2.size();
        idx[H[e]] = j;
        H2.push_back(H[e]);
        W2.push_back(1);
        memberPyIds.push_back({pyid});
        memberReps.push_back({REP[e]});
      } else {
        size_t j = it->second;
        W2[j] += 1;
        memberPyIds[j].push_back(pyid);
        memberReps[j].push_back(REP[e]);
      }
    }
    // Determine representative per new edge = REP of the member that is LAST in
    // CPython set-iteration order over its python ids.
    vector<long long> R2(H2.size());
    for (size_t j = 0; j < H2.size(); j++) {
      if (memberPyIds[j].size() == 1) { R2[j] = memberReps[j][0]; continue; }
      long long winner = pyset::lastInIter(memberPyIds[j]);
      // find winner's REP
      long long rep = memberReps[j].back();
      for (size_t t = 0; t < memberPyIds[j].size(); t++)
        if (memberPyIds[j][t] == winner) { rep = memberReps[j][t]; break; }
      R2[j] = rep;
    }
    H.swap(H2); W.swap(W2); REP.swap(R2);
    firstCall = false;
  };

  // core: LITERAL port of the Python core().
  //   Q = nodes with weighted incident degree < 2  (computed ONCE at entry)
  //   while Q: pop n; for incident e with <2 nodes, add e to EQ; del_node(n);
  //            while EQ: pop e (never deleted); if EQ fired, rescan ALL live
  //                      nodes and add those with unweighted incident-edge
  //                      count < 2 to Q.
  //   del_node(n): remove n from every incident edge (shrink); n no longer a
  //                node. Other nodes' incidence lists are NOT modified (Python
  //                only removes n from hyperedges[e], not from other nodes'
  //                node.Edge), so a live node's incident-edge COUNT is fixed
  //                for the duration of this core() call.
  // Empirically EQ never fires on our data; we still port it faithfully.
  auto core = [&]() {
    size_t E = H.size();
    unordered_map<vid_t, vector<size_t>> inc;     // node -> incident edge indices (fixed)
    for (size_t e = 0; e < E; e++)
      for (vid_t v : H[e]) inc[v].push_back(e);
    unordered_map<vid_t, long long> wdeg;         // weighted incident degree
    for (auto &kv : inc) {
      long long sw = 0;
      for (size_t e : kv.second) sw += W[e];
      wdeg[kv.first] = sw;
    }
    unordered_set<vid_t> deadV;
    unordered_set<vid_t> inQ;
    vector<vid_t> Q;
    for (auto &kv : wdeg) if (kv.second < 2) { Q.push_back(kv.first); inQ.insert(kv.first); }

    // helper: current node-count of an edge (H already reflects shrinks)
    while (!Q.empty()) {
      vid_t n = Q.back(); Q.pop_back(); inQ.erase(n);
      if (deadV.count(n)) continue;
      bool eqFired = false;
      for (size_t e : inc[n]) {
        if (H[e].size() < 2) eqFired = true;  // edge already shrunk below 2 nodes
      }
      // del_node(n): shrink incident edges, mark n dead
      deadV.insert(n);
      for (size_t e : inc[n]) {
        auto &he = H[e];
        auto it = lower_bound(he.begin(), he.end(), n);
        if (it != he.end() && *it == n) he.erase(it);
      }
      // EQ loop: only rescans (adds to Q) when EQ was non-empty
      if (eqFired) {
        for (auto &kv : inc) {
          vid_t v = kv.first;
          if (deadV.count(v) || inQ.count(v)) continue;
          if ((long long)kv.second.size() < 2) { Q.push_back(v); inQ.insert(v); }
        }
      }
    }
    (void)E;
  };

  // -------------------------------------------------------------------------
  // Initialize H/W/REP from original edges, then run reduction to fixpoint.
  // Python compares len(G.nodes) and len(G.hyperedges) each round.
  // -------------------------------------------------------------------------
  H.reserve(origKey.size());
  for (long long i = 0; i < E_orig; i++) {
    H.push_back(origKey[i]);
    W.push_back(1);
    REP.push_back(i + 1); // 1-based original id
  }

  auto countNodes = [&]() {
    unordered_set<vid_t> vs;
    for (auto &he : H) for (vid_t v : he) vs.insert(v);
    return vs.size();
  };
  auto countEdges = [&]() {
    // Python len(G.hyperedges): after uniqueEdges, edges are the merged ones.
    // We measure at the top of the loop after previous round; to mirror exactly
    // we count current H entries that are still >=1 (present). Python keeps all
    // edges in the dict (even shrunk); but the fixpoint check is right after
    // uniqueEdges+core. We instead count via the same pattern used below.
    return H.size();
  };

  while (true) {
    size_t cntV = countNodes();
    size_t cntE = countEdges();
    uniqueEdges();
    core();
    size_t nV = countNodes();
    size_t nE = countEdges();
    if (nV == cntV && nE == cntE) break;
  }
  // Python does NOT run an extra uniqueEdges after the loop; the loop breaks only
  // when the last core removed nothing, so H is already clean (all edges >=2
  // nodes). Enumeration below also skips any <2-node edge defensively.

  const long long Ered = (long long)H.size();

  // -------------------------------------------------------------------------
  // Build a compact HyperGraph for triangle enumeration; remap node ids to a
  // dense range and sort by degree (sortV) as the enumeration expects b<a to
  // range over shared-node edges. Enumeration is order-independent for the
  // resulting htn/adj multiset, so we may use the Baseline order directly.
  // -------------------------------------------------------------------------
  HyperGraph Gr;
  {
    unordered_map<vid_t, vid_t> remap;
    vid_t nid = 0;
    // deterministic: assign ids in ascending original node order
    {
      set<vid_t> vs;
      for (long long e = 0; e < Ered; e++) for (vid_t v : H[e]) vs.insert(v);
      for (vid_t v : vs) remap[v] = nid++;
    }
    Gr.E.reserve(Ered);
    for (long long e = 0; e < Ered; e++) {
      HyperEdge he; he.weight = (eid_t)W[e];
      he.resize(H[e].size());
      for (size_t j = 0; j < H[e].size(); j++) he[j] = remap[H[e][j]];
      sort(he.begin(), he.end());
      Gr.E.push_back(std::move(he));
    }
    Gr.V.assign(nid, {});
    for (eid_t e = 0; e < (eid_t)Ered; e++) for (auto &v : Gr.E[e]) Gr.V[v].push_back(e);
    for (auto &vv : Gr.V) sort(vv.begin(), vv.end());
  }

  // rep per reduced edge index (Gr.E index == H index; sortV reorders NODES only)
  vector<long long> repOf = REP;
  Gr.sortV();

  // -------------------------------------------------------------------------
  // Triangle enumeration (same definition as Baseline.cpp) + htn/adj accumulate.
  // -------------------------------------------------------------------------
  const eid_t Ecur = Gr.sizeE();
  unordered_map<long long, long long> htn;
  unordered_map<long long, unordered_map<long long, long long>> adj;
  htn.reserve(Ecur * 2);
  adj.reserve(Ecur * 2);

  vector<vector<vid_t>> source(Ecur);
  vector<eid_t> worklist;

  for (eid_t a = 0; a < Ecur; a++) {
    for (const auto &v : Gr.E[a]) {
      for (const auto &b : Gr.V[v]) {
        if (b >= a) break;
        if (source[b].empty()) worklist.emplace_back(b);
        source[b].emplace_back(v);
      }
    }
    for (size_t wi = 0; wi < worklist.size(); wi++) {
      const auto b = worklist[wi];
      for (size_t wj = 0; wj < wi; wj++) {
        const auto c = worklist[wj];
        if (!vectorIntersect(source[b], source[c]) && hyperedge_intersection(Gr.E[b], Gr.E[c])) {
          long long wa = (long long)Gr.E[a].weight;
          long long wb = (long long)Gr.E[b].weight;
          long long wc = (long long)Gr.E[c].weight;
          long long ra = repOf[a], rb = repOf[b], rc = repOf[c];
          htn[ra] += wb * wc;
          htn[rb] += wa * wc;
          htn[rc] += wa * wb;
          adj[ra][rb] += wc; adj[rb][ra] += wc;
          adj[ra][rc] += wb; adj[rc][ra] += wb;
          adj[rb][rc] += wa; adj[rc][rb] += wa;
        }
      }
    }
    for (const auto &w : worklist) source[w].clear();
    worklist.clear();
  }

  // -------------------------------------------------------------------------
  // Peeling, exactly as Python run().
  // -------------------------------------------------------------------------
  vector<long long> trussness(E_orig + 1, -1); // 1-based; -1 = unassigned
  vector<char> alive(E_orig + 1, 0);
  for (long long hid = 1; hid <= E_orig; hid++) alive[hid] = 1;
  long long aliveCount = E_orig;

  for (long long hid = 1; hid <= E_orig; hid++)
    if (htn.find(hid) == htn.end()) trussness[hid] = 0;

  long long kmax = (nNodes * (nNodes - 1)) / 2;
  long long k = 3;
  for (k = 3; k <= kmax; k++) {
    if (aliveCount == 0) break;
    deque<long long> q;
    for (long long e = 1; e <= E_orig; e++) {
      if (!alive[e]) continue;
      auto it = htn.find(e);
      long long hv = (it == htn.end()) ? 0 : it->second;
      if (hv < k - 2) q.push_back(e);
    }
    while (!q.empty()) {
      long long e = q.front(); q.pop_front();
      if (!alive[e]) continue;
      alive[e] = 0; aliveCount--;
      if (trussness[e] == -1) trussness[e] = k - 1;
      auto ait = adj.find(e);
      if (ait != adj.end()) {
        for (auto &pr : ait->second) {
          long long nb = pr.first, tcnt = pr.second;
          if (nb >= 1 && nb <= E_orig && alive[nb]) {
            htn[nb] -= tcnt;
            if (htn[nb] < k - 2) q.push_back(nb);
            auto nit = adj.find(nb);
            if (nit != adj.end()) nit->second.erase(e);
          }
        }
        adj.erase(ait);
      }
      htn.erase(e);
    }
  }
  for (long long e = 1; e <= E_orig; e++) if (alive[e]) trussness[e] = k;
  for (long long e = 1; e <= E_orig; e++) if (trussness[e] == -1) trussness[e] = 0;

  // -------------------------------------------------------------------------
  // Write CSV: edge_id,trussness for ALL original ids (1..E_orig), sorted.
  // -------------------------------------------------------------------------
  {
    ofstream fout(outPath);
    fout << "edge_id,trussness\n";
    for (long long e = 1; e <= E_orig; e++)
      fout << e << "," << trussness[e] << "\n";
  }

  cout << "Original edges: " << E_orig
       << " | reduced edges: " << Ered
       << " | original nodes: " << nNodes << endl;
  cout << "Output written to " << outPath << endl;
  cout << "Total runtime: "
       << chrono::duration<double>{chrono::steady_clock::now() - time_start}.count()
       << " s" << endl;
  return 0;
}
