import math
import re
import os
import csv
import time
import networkx as nx
from itertools import combinations
from collections import defaultdict
from statistics import mean, stdev

import func


# ─────────────────────────────────────────────
# Existing helpers (kept for backward compat)
# ─────────────────────────────────────────────

def num_nodes(G):
    return len(G.nodes)

def num_hyperedges(G):
    return len(G.hyperedges)

def density(G):
    if len(G.nodes) == 0:
        return 0
    return len(G.hyperedges) / len(G.nodes)

def stat_deg(G):
    if not G.nodes:
        return 0, 0, 0.0, 0.0
    degrees = [len(node.Edge) for node in G.nodes.values()]
    return min(degrees), max(degrees), mean(degrees), (stdev(degrees) if len(degrees) > 1 else 0.0)

def stat_card(G):
    if not G.hyperedges:
        return 0, 0, 0.0, 0.0
    sizes = [len(e) for e in G.hyperedges.values()]
    return min(sizes), max(sizes), mean(sizes), (stdev(sizes) if len(sizes) > 1 else 0.0)


# ─────────────────────────────────────────────
# Section 4.2  Basic size measures
# ─────────────────────────────────────────────

def node_retention_ratio(H_orig, H_result):
    n = len(H_orig.nodes)
    return len(H_result.nodes) / n if n > 0 else 0.0

def edge_retention_ratio(H_orig, H_result):
    m = len(H_orig.hyperedges)
    return len(H_result.hyperedges) / m if m > 0 else 0.0


# ─────────────────────────────────────────────
# Section 4.3  Degree / cardinality measures
# ─────────────────────────────────────────────

def avg_degree(H):
    if not H.nodes:
        return 0.0
    return sum(len(v.Edge) for v in H.nodes.values()) / len(H.nodes)

def avg_cardinality(H):
    if not H.hyperedges:
        return 0.0
    return sum(len(e) for e in H.hyperedges.values()) / len(H.hyperedges)

def total_incidence(H):
    return sum(len(e) for e in H.hyperedges.values())


# ─────────────────────────────────────────────
# Section 4.4  Density measures
# ─────────────────────────────────────────────

def incidence_density(H):
    n = len(H.nodes)
    m = len(H.hyperedges)
    if n == 0 or m == 0:
        return 0.0
    return total_incidence(H) / (n * m)

def density_ratio(H_orig, H_result):
    orig_d = density(H_orig)
    if orig_d == 0:
        return 0.0
    return density(H_result) / orig_d


# ─────────────────────────────────────────────
# Section 4.5  Neighbor measures
# ─────────────────────────────────────────────

def avg_num_neighbors(H):
    if not H.nodes:
        return 0.0
    total = 0
    for nid, node in H.nodes.items():
        nbrs = set()
        for eid in node.Edge:
            if eid in H.hyperedges:
                for v in H.hyperedges[eid]:
                    if v != nid:
                        nbrs.add(v)
        total += len(nbrs)
    return total / len(H.nodes)


# ─────────────────────────────────────────────
# Section 4.6  Support measures
# ─────────────────────────────────────────────

def _build_co(H):
    co = defaultdict(int)
    for edge_nodes in H.hyperedges.values():
        nodes = sorted(edge_nodes)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                co[(nodes[i], nodes[j])] += 1
    return co

def avg_support(H):
    co = _build_co(H)
    if not co:
        return 0.0
    return sum(co.values()) / len(co)

def max_support(H):
    co = _build_co(H)
    return max(co.values()) if co else 0


# ─────────────────────────────────────────────
# Section 4.7  DCS measures
# ─────────────────────────────────────────────

def avg_dcs(H, c=1.0):
    dcs = defaultdict(float)
    for edge_nodes in H.hyperedges.values():
        d = len(edge_nodes)
        if d < 2:
            continue
        w = 1.0 / (d ** c)
        nodes = sorted(edge_nodes)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                dcs[(nodes[i], nodes[j])] += w
    if not dcs:
        return 0.0
    return sum(dcs.values()) / len(dcs)


# ─────────────────────────────────────────────
# Section 4.8  Connected component measures
# ─────────────────────────────────────────────

def _clique_graph(H):
    G = nx.Graph()
    G.add_nodes_from(H.nodes.keys())
    for edge_nodes in H.hyperedges.values():
        nodes = list(edge_nodes)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                G.add_edge(nodes[i], nodes[j])
    return G

def num_connected_components(H):
    if not H.nodes:
        return 0
    return nx.number_connected_components(_clique_graph(H))

def avg_connected_component_size(H):
    if not H.nodes:
        return 0.0
    comps = list(nx.connected_components(_clique_graph(H)))
    return sum(len(c) for c in comps) / len(comps) if comps else 0.0

def largest_connected_component_ratio(H):
    if not H.nodes:
        return 0.0
    comps = list(nx.connected_components(_clique_graph(H)))
    return max(len(c) for c in comps) / len(H.nodes) if comps else 0.0


# ─────────────────────────────────────────────
# Section 4.9  Cut-based measures
# ─────────────────────────────────────────────

def _split_cut(H_orig, S):
    """Split-based cut weight: sum_e |e∩S|*|e\\S| / |e|."""
    cut = 0.0
    for edge_nodes in H_orig.hyperedges.values():
        inter = sum(1 for v in edge_nodes if v in S)
        diff = len(edge_nodes) - inter
        if inter > 0 and diff > 0:
            cut += inter * diff / len(edge_nodes)
    return cut

def _vol(H_orig, S):
    return sum(len(H_orig.nodes[v].Edge) for v in S if v in H_orig.nodes)

def conductance(H_orig, H_result):
    if not H_result.nodes or not H_orig.nodes:
        return 0.0
    S = set(H_result.nodes.keys())
    if len(S) == len(H_orig.nodes):
        return 0.0
    cut = _split_cut(H_orig, S)
    vol_S = _vol(H_orig, S)
    vol_Sbar = _vol(H_orig, set(H_orig.nodes.keys()) - S)
    denom = min(vol_S, vol_Sbar)
    return cut / denom if denom > 0 else 0.0

def normalized_cut(H_orig, H_result):
    if not H_result.nodes or not H_orig.nodes:
        return 0.0
    S = set(H_result.nodes.keys())
    if len(S) == len(H_orig.nodes):
        return 0.0
    cut = _split_cut(H_orig, S)
    vol_S = _vol(H_orig, S)
    vol_Sbar = _vol(H_orig, set(H_orig.nodes.keys()) - S)
    r1 = cut / vol_S if vol_S > 0 else 0.0
    r2 = cut / vol_Sbar if vol_Sbar > 0 else 0.0
    return r1 + r2


# ─────────────────────────────────────────────
# Section 4.10  Modularity (Kaminski et al. 2019)
# ─────────────────────────────────────────────

def modularity(H_orig, H_result):
    """
    Two-partition hypergraph modularity.
    Partition: S = H_result nodes, S_bar = remaining nodes in H_orig.
    Q = (1/|E|) * sum over {S, S_bar} of [e_H(A) - sum_{d>=2} |E_d|*(vol(A)/vol(V))^d]
    """
    if not H_result.nodes or not H_orig.nodes:
        return 0.0
    m_total = len(H_orig.hyperedges)
    if m_total == 0:
        return 0.0

    S = set(H_result.nodes.keys())
    Sbar = set(H_orig.nodes.keys()) - S

    vol_S    = sum(len(H_orig.nodes[v].Edge) for v in S    if v in H_orig.nodes)
    vol_Sbar = sum(len(H_orig.nodes[v].Edge) for v in Sbar if v in H_orig.nodes)
    vol_V = vol_S + vol_Sbar
    if vol_V == 0:
        return 0.0

    # Count hyperedges by size in H_orig
    edge_size_count = defaultdict(int)
    for edge_nodes in H_orig.hyperedges.values():
        edge_size_count[len(edge_nodes)] += 1

    # e_H(A): number of hyperedges with ALL nodes in A
    e_S    = sum(1 for en in H_orig.hyperedges.values() if all(v in S    for v in en))
    e_Sbar = sum(1 for en in H_orig.hyperedges.values() if all(v in Sbar for v in en))

    # Expected contribution
    expected = 0.0
    for d, cnt in edge_size_count.items():
        if d >= 2:
            expected += cnt * ((vol_S / vol_V) ** d + (vol_Sbar / vol_V) ** d)

    return ((e_S + e_Sbar) - expected) / m_total


# ─────────────────────────────────────────────
# Aggregated compute functions
# ─────────────────────────────────────────────

def compute_size_measures(H_orig, H_result):
    """EQ1/EQ2: basic size measures only."""
    return {
        "num_nodes":      num_nodes(H_result),
        "num_hyperedges": num_hyperedges(H_result),
    }

def compute_all_measures(H_orig, H_result, c_dcs=1.0):
    """EQ3/EQ4: all effectiveness measures, in canonical output order."""
    ZERO_KEYS = [
        "density", "density_ratio",
        "avg_degree", "avg_num_neighbors",
        "avg_support", "modularity",
        "conductance", "normalized_cut",
        "avg_connected_component_size",
        "avg_dcs",
        "num_nodes", "num_hyperedges", "avg_cardinality",
    ]

    if num_nodes(H_result) == 0:
        return {k: 0.0 for k in ZERO_KEYS}

    return {
        "density":                     density(H_result),
        "density_ratio":               density_ratio(H_orig, H_result),
        "avg_degree":                  avg_degree(H_result),
        "avg_num_neighbors":           avg_num_neighbors(H_result),
        "avg_support":                 avg_support(H_result),
        "modularity":                  modularity(H_orig, H_result),
        "conductance":                 conductance(H_orig, H_result),
        "normalized_cut":              normalized_cut(H_orig, H_result),
        "avg_connected_component_size": avg_connected_component_size(H_result),
        "avg_dcs":                     avg_dcs(H_result, c=c_dcs),
        "num_nodes":                   num_nodes(H_result),
        "num_hyperedges":              num_hyperedges(H_result),
        "avg_cardinality":             avg_cardinality(H_result),
    }
