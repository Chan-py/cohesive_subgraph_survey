import time
import networkx as nx

def run(G, k, q):
    s_time = time.time()

    VC = {u for u in G.nodes if len(G.nodes[u].Edge) < k}
    EC = {eid for eid in G.hyperedges if len(G.hyperedges[eid]) < q}

    while VC or EC:
        while VC:
            v = VC.pop()
            if v not in G.nodes:
                continue

            incident_edges = list(G.nodes[v].Edge)
            G.del_node(v)

            for eid in incident_edges:
                if eid in G.hyperedges and len(G.hyperedges[eid]) < q:
                    EC.add(eid)

        while EC:
            e = EC.pop()
            if e not in G.hyperedges:
                continue

            incident_nodes = list(G.hyperedges[e])
            G.del_edge(e)

            for u in incident_nodes:
                if u in G.nodes and len(G.nodes[u].Edge) < k:
                    VC.add(u)

    e_time = time.time() - s_time
    return set(G.nodes), e_time
