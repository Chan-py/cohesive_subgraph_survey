import time
import copy

import func
import kg_core

def updatesuptable(ng):
    M = {}
    for node in ng:
        occurence = ng[node]
        if occurence not in M:
            M[occurence] = 0
        M[occurence] += 1
    return M

def getK(M):
    if len(M) == 0:
        return 0
    sum = 0
    for occurence in M:
        sum += M[occurence]
    return sum

def getEdgeLB(v, M, g):
    c = v.EdgeCnt
    T = {}
    for s in M:
        if s - c >= g:
            T[s - c] = M[s]
    v.EdgeCnt = 0
    return getK(T), T

def getNodeLB(v, M):
    c = v.NodeCnt
    keys = list(M.keys())
    keys.sort()
    while c != 0:
        if len(M) == 0:
            return 0
        mk = keys.pop()
        if c > M[mk]:
            c -= M[mk]
            del M[mk]
        else:
            M[mk] -= c
            if M[mk] == 0:
                del M[mk]
            break
    v.NodeCnt = 0
    return getK(M)

def run(G, k, g, p):
    G2 = copy.deepcopy(G)
    G1, S, kg_time = kg_core.run(G2, k, g)
    s_time = time.time()

    PQ = func.PriorityQueue()
    VC = set()
    EC = set(G1.hyperedges.keys())
    M = {}
    for v in S:
        M[v] = {g: S[v]}
    while EC:
        for edge in EC:
            len_G1 = len(G1.hyperedges[edge])
            if len_G1 / len(G.hyperedges[edge]) < p:
                for v in G1.hyperedges[edge]:
                    if len_G1 -1 == 0:
                        break
                    G1.nodes[v].EdgeCnt = G1.nodes[v].EdgeCnt + 1
                    VC.add(v)
                G1.del_edge(edge)
        EC.clear()
        for v in VC:
            a, M[v] = getEdgeLB(G1.nodes[v], M[v], g)
            if a < k:
                ng = func.getgNbrMap(G1, v, g)
                if len(ng) < k:
                    for e in G1.nodes[v].Edge:
                        EC.add(e)
                    G1.del_node(v)
                    if PQ.contains(v):
                        PQ.remove(v)
                    del M[v]
                    for w in ng:
                        G1.nodes[w].NodeCnt = G1.nodes[w].NodeCnt + 1
                        if PQ.contains(w):
                            PQ.remove(w)
                        PQ.push(w, G1.nodes[w].NodeCnt)
                else:
                    M[v] = updatesuptable(ng)
        VC.clear()
        while not PQ.empty():
            v = PQ.pop()
            a = getNodeLB(G1.nodes[v], M[v])
            if a < k:
                ng = func.getgNbrMap(G1, v, g)
                if len(ng) < k:
                    for e in G1.nodes[v].Edge:
                        EC.add(e)
                    G1.del_node(v)
                    if PQ.contains(v):
                        PQ.remove(v)
                    del M[v]
                    for w in ng:
                        G1.nodes[w].NodeCnt = G1.nodes[w].NodeCnt + 1
                        if PQ.contains(w):
                            PQ.remove(w)
                        PQ.push(w, G1.nodes[w].NodeCnt)
                else:
                    M[v] = updatesuptable(ng)

    e_time = time.time() - s_time + kg_time
    return G1, e_time




