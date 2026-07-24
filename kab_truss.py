import time
import networkx as nx
from collections import deque

import func

def initial_support(G, a, b):
    sup = {}
    deg = G.degree
    for u, v in list(G.edges()):
        key = (min(u, v), max(u, v))
        if deg[u] > deg[v]:
            u, v = v, u
        Nu = set(G.neighbors(u))
        Nv = set(G.neighbors(v))
        for w in (Nu & Nv):
            if G.has_edge(u, w) and G.has_edge(v, w):
                Huv = G[u][v]["hids"]
                Hvw = G[w][v]["hids"]
                Huw = G[u][w]["hids"]
                if len(Huv & Hvw & Huw) < a: continue
                if len(Huv | Hvw | Huw) < b: continue
                sup[key] = sup.get(key, 0) + 1
    return sup

def ce_truss(H, k, a, b, network):
    G, proj_time = func.load_clique_graph(H, network)
    s_time = time.time()
    sup = initial_support(G, a, b)

    Q = deque()
    for (u, v) in G.edges():
        key = (min(u, v), max(u, v))
        if sup.get(key, 0) < k:
            Q.append(key)

    removed_edges = 0
    while Q:
        u, v = Q.popleft()
        if not G.has_edge(u, v):
            continue
        Huv = G[u][v]["hids"]
        G.remove_edge(u, v)
        removed_edges += 1
        if G.degree(u) > G.degree(v):
            u, v = v, u

        Nu = set(G.neighbors(u))
        Nv = set(G.neighbors(v))
        for w in (Nu & Nv):
            if G.has_edge(u, w) and G.has_edge(v, w):
                Hvw = G[w][v]["hids"]
                Huw = G[u][w]["hids"]
                if len(Huv & Hvw & Huw) < a: continue
                if len(Huv | Hvw | Huw) < b: continue
                k_uw = (min(u, w), max(u, w))
                k_vw = (min(v, w), max(v, w))
                sup[k_uw] = sup.get(k_uw, 0) - 1
                if sup[k_uw] < k:
                    Q.append(k_uw)

                sup[k_vw] = sup.get(k_vw, 0) - 1
                if sup[k_vw] < k:
                    Q.append(k_vw)

    e_time = time.time() - s_time + proj_time
    G.remove_nodes_from([n for n, d in G.degree() if d == 0])
    hids = set().union(*(G[u][v]["hids"] for (u, v) in G.edges()))
    return G.nodes, hids, e_time

def TN(u, v, id2pair, pair2id):
    pu = set(id2pair[u])
    pv = set(id2pair[v]) 
    inter = pu & pv
    if len(inter) != 1:
        return None
    a = (pu - inter).pop()
    b = (pv - inter).pop()
    target = tuple(sorted((a, b)))
    return pair2id.get(target, None)

def get_support_ab(G, a, b, id2pair, pair2id):
    visited = {}
    for u, v in G.edges():
        visited[(min(u, v), max(u, v))] = False
    sup = {}
    for u in G.nodes():
        for v in G.neighbors(u):
            if not visited[(min(u, v), max(u, v))]:
                w = TN(u, v, id2pair, pair2id)
                if w is None or not G.has_node(w): continue
                if G.has_node(u) and G.has_node(v):
                    if G.has_edge(u, v):
                        if G[u][v]['sigma'] >= a and (G.nodes[u]['omega'] + G.nodes[v]['omega'] + G.nodes[w]['omega'] - 2 * G[u][v]['sigma']) >= b:
                            sup[u] = sup.get(u, 0) + 1
                            sup[v] = sup.get(v, 0) + 1
                            sup[w] = sup.get(w, 0) + 1
                key = sorted((u, v, w))
                visited[(key[0], key[1])] = True
                visited[(key[0], key[2])] = True
                visited[(key[1], key[2])] = True
    return sup



def run(H, k, a, b, network):
    if a == 0:
        return ce_truss(H, k, a, b, network)
    G1, proj_time, id2pair, pair2id = func.load_projection_graph(H, network)

    s_time = time.time()
    Q = set()
    G = nx.k_core(G1, 2*k)
    
    sup = get_support_ab(G, a, b, id2pair, pair2id)

    Q = {u for u in G.nodes() if sup.get(u, 0) < k}
    while Q:
        # print(len(Q))
        u = Q.pop()
        if u not in G: continue
        for v in list(G.neighbors(u)):
            if v not in sup: continue
            w = TN(u, v, id2pair, pair2id)
            if w is None: continue
            if G.has_node(u) and G.has_node(v) and G.has_node(w):
                if G.has_edge(u, v) and G.has_edge(u, w) and G.has_edge(v, w):
                    if G[u][v]['sigma'] >= a and (G.nodes[u]['omega'] + G.nodes[v]['omega'] + G.nodes[w]['omega'] - 2 * G[u][v]['sigma']) >= b:
                        sup[v] -= 1
                    if sup[v] < k: Q.add(v)
        G.remove_node(u)
        # print('remove', id2pair[u])
    e_time = time.time() - s_time + proj_time
    nodes = set()
    hids = set()
    iddd = set()
    for u in G.nodes():
        nodes.update(id2pair[u])
        iddd.add(id2pair[u])
        hids.update(G.nodes[u]['hids'])
    # print(e_time)
    # print(nodes)
    # print(hyperedges)
    # print(hids)
    # print(iddd)
    return nodes, hids, e_time