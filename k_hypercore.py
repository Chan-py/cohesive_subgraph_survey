import csv
import time
from pathlib import Path

import func
import k_hypercoreness
import etc

def run(network, dataset, k):
    s_time = time.time()
    decomp_time = func.read_decomp_from_csv(dataset, "k_hypercoreness")
    e_dec = 0
    if not decomp_time:
        dec = time.time()
        G = func.Hypergraph()
        G.load_from_file(network)

        coreness, decomp_time = k_hypercoreness.run(G, dataset)
        etc.dict_to_csv(coreness, 'node_id', 'coreness', f'../output/{dataset}_k_hypercoreness.csv')
        e_dec = time.time() - dec
    file_path = Path(f"../output/{dataset}_k_hypercoreness.csv")
    nodes = set()
    with open(file_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                kk = int(row["coreness"])
                if kk >= k:
                    node = int(row["node_id"])
                    nodes.add(node)
    e_time = (time.time() - s_time) + decomp_time - e_dec
    return nodes, e_time

# def run(network, dataset, k):
#     # G = copy.deepcopy(G_)
#     G = func.Hypergraph()
#     G.load_from_file(network)
#     R = {v for v in G.nodes if len(G.nodes[v].Edge) < k}
#     remove_nodes = set()

#     while R:
#         remove_edge = set()
#         for v in R:
#             for eid in G.nodes[v].Edge:
#                 remove_edge.add(eid)
#             G.del_node(v)
#             remove_nodes.add(v)
#         R_ = set()
#         for eid in remove_edge:
#             for u in G.hyperedges[eid]:
#                 if u in remove_nodes: continue
#                 if len(G.nodes[u].Edge)-1 < k:
#                     R_.add(u)
#             G.del_edge(eid)
#         R = R_

#     return G, 0