import time
from collections import defaultdict
from itertools import combinations

import func


def initial_support(G):
    sup_in = defaultdict(int)
    sup_out = defaultdict(int)

    for u, v in list(G.edges()):
        u, v = func.get_key(u, v)
        Nu = set(G.neighbors(u))
        Nv = set(G.neighbors(v))
        com_nb = Nu & Nv
        for e in G[u][v]['hids']:
            for w in com_nb:
                huw = G[u][w]['hids']
                hvw = G[v][w]['hids']
                inn = 0
                if (e in huw) and (e in hvw):
                    sup_in[(e, u, v)] += 1
                    inn = 1
                sup_out[(e, u, v)] += (len(huw) * len(hvw) - inn)
    sup_node_out = defaultdict(int)
    for (e, u, v), sup in sup_out.items():
        if sup_node_out[(e, u)] < sup:
            sup_node_out[(e, u)] = sup
        if sup_node_out[(e, v)] < sup:
            sup_node_out[(e, v)] = sup
    return sup_in, sup_out, sup_node_out


def _drop_isolated_nodes(G):
    for n in [n for n in G.nodes() if G.degree(n) == 0]:
        G.remove_node(n)


def _remove_hyperedge(G, hyperedges, e):
    """Remove hyperedge e from the clique graph and the hyperedges dict, keeping the two in sync."""
    for x, y in combinations(hyperedges[e], 2):
        a, b = func.get_key(x, y)
        if not G.has_edge(a, b):
            continue
        if e in G[a][b]['hids']:
            G[a][b]['hids'].discard(e)
            if not G[a][b]['hids']:
                G.remove_edge(a, b)
    del hyperedges[e]


def _remove_vertex_from_hyperedge(G, hyperedges, e, u):
    """Remove vertex u from hyperedge e, dropping the now-stale incidences from the clique graph."""
    comembers = [w for w in hyperedges[e] if w != u]
    hyperedges[e].discard(u)
    for w in comembers:
        a, b = func.get_key(u, w)
        if G.has_edge(a, b) and e in G[a][b]['hids']:
            G[a][b]['hids'].discard(e)
            if not G[a][b]['hids']:
                G.remove_edge(a, b)
    if len(hyperedges[e]) < 2:
        _remove_hyperedge(G, hyperedges, e)


def _violating_triples(G, hyperedges, sup, threshold):
    """All (e, u, v) incidences currently in the graph whose support is below threshold.

    Iterating the graph (rather than the support defaultdict) is essential: incidences with
    support 0 are absent from the defaultdict but must still be enqueued for peeling.
    """
    Q = set()
    for u, v in list(G.edges()):
        a, b = func.get_key(u, v)
        for e in G[a][b]['hids']:
            if e in hyperedges and sup[(e, a, b)] < threshold:
                Q.add((e, a, b))
    return Q


def peel(G, hyperedges, sup_in, sup_out, sup_node_out, k_in, k_out):
    change = True
    while change:
        change = False

        # ---- in-peeling: drop whole hyperedges whose in-support is too small ----
        Q_in = _violating_triples(G, hyperedges, sup_in, k_in)
        while Q_in:
            e, u, v = Q_in.pop()
            if e not in hyperedges:
                continue
            _remove_hyperedge(G, hyperedges, e)
            change = True

        # ---- out-peeling: drop vertices whose node out-support is too small ----
        Q_out = _violating_triples(G, hyperedges, sup_out, k_out)
        while Q_out:
            e, u, v = Q_out.pop()
            if e not in hyperedges:
                continue
            # peel the endpoint with the smaller node out-support
            if sup_node_out[(e, u)] > sup_node_out[(e, v)]:
                u, v = v, u
            if sup_node_out[(e, u)] >= k_out:
                continue
            if u in hyperedges[e]:
                _remove_vertex_from_hyperedge(G, hyperedges, e, u)
                change = True

        _drop_isolated_nodes(G)

        # recompute support from scratch for the next round
        if change:
            sup_in, sup_out, sup_node_out = initial_support(G)

    for u, v in list(G.edges()):
        if not G[u][v]['hids']:
            G.remove_edge(u, v)
    _drop_isolated_nodes(G)
    nodes = {n for n in G.nodes() if G.degree[n] != 0}
    hids = set().union(*(G[u][v]["hids"] for (u, v) in G.edges())) if G.number_of_edges() else set()
    return nodes, hids


def run(H, network, k_in, k_out):
    G, proj_time = func.load_clique_graph(H, network)
    hyperedges = H.hyperedges
    sup_in, sup_out, sup_node_out, sup_time = func.load_support_inout(G, network)
    s_time = time.time()

    nodes, hids = peel(G, hyperedges, sup_in, sup_out, sup_node_out, k_in, k_out)

    e_time = time.time() - s_time + proj_time + sup_time
    return nodes, hids, e_time
