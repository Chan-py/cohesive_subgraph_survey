import time

import func

def run(G, dataset):
    s_time = time.time()
    maxKey = 0
    deg = {}
    flag = {}
    core = {}

    for node_id, node in G.nodes.items():
        s = len(node.Edge)
        deg[node_id] = s
        flag[node_id] = 0
        if s > maxKey:
            maxKey = s
    bin_count = [0] * (maxKey + 1)
    for node_id in G.nodes:
        bin_count[deg[node_id]] += 1
    start = 0
    for k in range(maxKey + 1):
        num = bin_count[k]
        bin_count[k] = start
        start += num
    vert = [None] * len(G.nodes) 
    pos = {}
    for node_id in G.nodes:
        k = deg[node_id] 
        pos[node_id] = bin_count[k]
        vert[pos[node_id]] = node_id
        bin_count[k] += 1  
    
    for k in range(maxKey, 0, -1):
        bin_count[k] = bin_count[k - 1]
    bin_count[0] = 0

    edge_removed = {h: False for h in G.hyperedges}
    for i in range(len(G.nodes)):
        v = vert[i]
        flag[v] = 1
        core[v] = deg[v]

        for h in G.nodes[v].Edge:
            if edge_removed[h]:
                continue
            edge_removed[h] = True

            for u in G.hyperedges[h]:
                if flag[u] == 0 and deg[u] > deg[v]:
                    ku = deg[u]
                    pu = pos[u]
                    pw = bin_count[ku]
                    w = vert[pw]
                    if u != w:
                        vert[pu], vert[pw] = vert[pw], vert[pu]
                        pos[u], pos[w] = pw, pu
                    deg[u] -= 1
                    bin_count[ku] += 1
    e_time = time.time() - s_time
    func.save_time("k_hypercoreness", dataset, e_time)
    return core, e_time
