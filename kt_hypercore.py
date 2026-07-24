import time
import pickle
from pathlib import Path
import math

from kt_hypercoreness import preprocess_from_hypergraph

def run(G, k, t, dataset):
    data_dir = '../cache/kt_pre'
    ds_dir = Path(data_dir) / dataset
    i2edges_path = ds_dir / "i2edges.pkl"
    v2edges_path = ds_dir / "v2edges.pkl"
    if not (i2edges_path.exists() and v2edges_path.exists()):
        preprocess_from_hypergraph(G, data_dir, dataset)

    with (i2edges_path).open('rb') as f:
        i2edges = pickle.load(f)   # {edge_id: set(nodes)}
    with (v2edges_path).open('rb') as f:
        v2edges = pickle.load(f)   # {node_id: set(edge_ids)}

    s_time = time.time()
    i2th = dict()
    for i_e, e in i2edges.items():
        i2th[i_e] = max(math.ceil(t * len(e)), 2)
    nodes_to_remove = set()
    for v, E_v in v2edges.items():
        if len(E_v) < k:
            nodes_to_remove.add(v)
    while nodes_to_remove:
        edges_to_check = set()
        for v in nodes_to_remove:
            edges_to_check.update(v2edges[v])
            del v2edges[v]
        assert edges_to_check.issubset(i2edges)
        nodes_to_remove_next = set()
        for i_e in edges_to_check:
            i2edges[i_e] -= nodes_to_remove
            if len(i2edges[i_e]) < i2th[i_e]:
                for v in i2edges[i_e]:
                    if v not in v2edges:
                        continue
                    if len(v2edges[v])-1 < k:
                        nodes_to_remove_next.add(v)
                    v2edges[v].remove(i_e)
                del i2edges[i_e]
        nodes_to_remove = nodes_to_remove_next.copy()
    e_time = time.time() - s_time
    print(len(v2edges.keys()))
    print(len(i2edges.keys()))
    return v2edges.keys(), i2edges.keys(), e_time

# def run(G, k, t, dataset):
#     EQ = set()
#     VQ = set()
#     D = {eid: len(edge) for eid, edge in G.hyperedges.items()}
#     # i2th = {}
#     # for i_e, e in G.hyperedges.items():
#     #     i2th[i_e] = max(math.ceil(t * len(e)), 2)

#     for node in G.nodes:
#         if len(G.nodes[node].Edge) < k:
#             EQ.update(G.nodes[node].Edge)
#             VQ.add(node)
#     while EQ:
#         while VQ:
#             node = VQ.pop()
#             G.del_node(node)

#         eid = EQ.pop()
#         if  len(G.hyperedges[eid])/D[eid] < t:
#             node_set = G.hyperedges[eid]
#             G.del_edge(eid)
#             for n in node_set:
#                 if len(G.nodes[n].Edge)-1 < k:
#                     VQ.add(n)
#                     EQ.update(G.nodes[n].Edge)
#     return G.nodes, G.hyperedges, 0
