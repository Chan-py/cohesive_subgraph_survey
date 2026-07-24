import queue
import time

def getsNbrMap(G1, map, node, s, c):
    cnt = {}
    for hyperedge in G1.nodes[node].Edge:
        for neighbor in G1.hyperedges[hyperedge]:
            if neighbor != node:
                if neighbor not in cnt:
                    cnt[neighbor] = 0
                cnt[neighbor] += map[hyperedge]
    ng = {node: count for node, count in cnt.items() if count >= s}
    return ng

def get_map(G, c):
    map = {}
    for id, edge_set in G.hyperedges.items():
        map[id] = 1/(len(edge_set))**c
    return map

def run(G1, k, s, c):
    s_time = time.time()
    map = get_map(G1, c)
    # ss = {}
    # for node in G1.nodes:
    #     cnt = getsNbrMap(G1, map, node, 0, 1)
    #     for n, v in cnt.items():
    #         ss[(min(node, n),max(node, n))] = v
    # print(ss)
    # exit()
    CS = {}
    VQ = queue.Queue()
    VQ1 = set()
    for node in G1.nodes:
        nb = getsNbrMap(G1, map, node, s, c)
        CS[node] = len(nb)
        if CS[node] < k:
            VQ.put(node)
            VQ1.add(node)
    while not VQ.empty():
        v = VQ.get()
        VQ1.remove(v)
        nb = getsNbrMap(G1, map, v, s, c)
        G1.del_node(v)
        del CS[v]
        for w in nb:
            if w not in VQ1:
                CS[w] -= 1
                if CS[w] < k:
                    VQ.put(w)
                    VQ1.add(w)

    e_time = time.time() - s_time
    return G1.nodes, e_time