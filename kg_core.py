import copy
import queue
import time

import func

def run(G, k, g):
    s_time = time.time()

    S = {}
    VQ = queue.Queue()
    VQ1 = set()
    for node in G.nodes:
        ng = func.getgNbrMap(G, node, g)
        S[node] = len(ng)
        if S[node] < k:
            VQ.put(node)
            VQ1.add(node)
    while not VQ.empty():
        v = VQ.get()
        VQ1.remove(v)
        ng = func.getgNbrMap(G, v, g)
        G.del_node(v)
        del S[v]
        for w in ng:
            if w not in VQ1:
                S[w] -= 1
                if S[w] < k:
                    VQ.put(w)
                    VQ1.add(w)

    e_time = time.time() - s_time
    return G, S, e_time