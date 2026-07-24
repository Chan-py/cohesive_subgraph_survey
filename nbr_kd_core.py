import csv
from collections import defaultdict
import time
from pathlib import Path
import re
import os
import copy
import queue    

import func
import nbr_kd_coreness

def run(network, dataset, k, d, algorithm):
    ns = network.split('/')
    s_time = time.time()
    decomp_time = func.read_decomp_from_csv(dataset, f"{algorithm}ness")
    e_dec = 0
    if not decomp_time:
        dec = time.time()
        # dataset path relative to datasets/ (any depth: real/<name> OR synthetic/eq3/<name>)
        di = ns.index("datasets")
        network = "/".join(ns[di + 1:-1])
        decomp_time = nbr_kd_coreness.run(network, f"{algorithm}ness")
        e_dec = time.time() - dec
    file_path = Path(f"../output/{dataset}_{algorithm}ness.csv")

    if algorithm == "nbr_k_core":
        # nbr[k] = {node, ...}
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

    elif algorithm == "kd_core":
        # kd[k][d] = {node, ...}
        kd = defaultdict(lambda: defaultdict(set))
        PAIR_RE = re.compile(r'\(\s*(\d+)\s*,\s*(\d+)\s*\)')
        nodes = set()
        with open(file_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                node = int(row["node_id"])
                coreness = row['coreness']
                for k_str, d_str in PAIR_RE.findall(coreness):
                    if int(k_str) >= k and int(d_str) >= d:
                        nodes.add(node)
                        continue

        e_time = (time.time() - s_time) + decomp_time
        return nodes, e_time

def kdCore(G, k, d):
    G1 = copy.deepcopy(G)
    VQ = queue.Queue()

    for node in G1.nodes:
        # degree < d인 노드 제거
        if len(G1.nodes[node].Edge) < d:
            VQ.put(node)
        # num_of_neighbors < k인 노드 제거
        else:
            neighbors = set()
            for edge in G1.nodes[node].Edge:
                for neighbor in G1.hyperedges[edge]:
                    if neighbor == node:
                        continue
                    neighbors.add(neighbor)
            if len(neighbors) < k:
                VQ.put(node)

    NS = set()  # set of neighbors
    while not VQ.empty():
        v = VQ.get()
        ss = set()
        for edge in G1.nodes[v].Edge:
            for neighbor in G1.hyperedges[edge]:
                if neighbor == v:
                    continue
                NS.add(neighbor)
                # G1.nodes[neighbor].Edge.remove(edge)
                ss.add(edge)
        for edge in ss:
            G1.del_edge(edge)

        G1.del_node(v)


        # k,d constraint 재확인
        if VQ.empty():
            for u in NS:
                if u in G1.nodes:
                    if len(G1.nodes[u].Edge) < d:
                        VQ.put(u)
                    else:
                        neighbors = set()
                        for edge in G1.nodes[u].Edge:
                            for neighbor in G1.hyperedges[edge]:
                                if neighbor == u:
                                    continue
                                neighbors.add(neighbor)
                        if len(neighbors) < k:
                            VQ.put(u)
            NS = set()

    return G1