import pickle
from pathlib import Path
import numpy as np
import time
import math

import func

def preprocess_from_hypergraph(H, out_dir, ds_name):
    out_path = Path(out_dir) / ds_name
    out_path.mkdir(parents=True, exist_ok=True)

    i2edges = {eid: set(nodes) for eid, nodes in H.hyperedges.items()}
    v2edges = {vid: set(node.Edge) for vid, node in H.nodes.items()}

    with (out_path / "i2edges.pkl").open("wb") as f:
        pickle.dump(i2edges, f)
    with (out_path / "v2edges.pkl").open("wb") as f:
        pickle.dump(v2edges, f)

def coreness(ds_dir):
    s_time = time.time()
    results = {}
    with (ds_dir / 'i2edges.pkl').open('rb') as f:
        i2edges = pickle.load(f)
    with (ds_dir / 'v2edges.pkl').open('rb') as f:
        v2edges = pickle.load(f)
    n, m = len(v2edges), len(i2edges)
    possible_t = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for i, t in enumerate(possible_t):
        v2coreness_t = dict()
        if i:
            with (ds_dir / 'i2edges.pkl').open('rb') as f:
                i2edges = pickle.load(f)
            with (ds_dir / 'v2edges.pkl').open('rb') as f:
                v2edges = pickle.load(f)
        i2th = dict()
        for i_e, e in i2edges.items():
            i2th[i_e] = max(math.ceil(t * len(e)), 2)
        k = None
        nodes_to_remove = set()
        while i2edges:
            if not nodes_to_remove:
                k = min(len(E_v) for E_v in v2edges.values()) + 1
                nodes_to_remove = set()
                for v, E_v in v2edges.items():
                    if len(E_v) == k - 1:
                        nodes_to_remove.add(v)
                        v2coreness_t[v] = k - 1
                continue
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
                        if len(v2edges[v]) == k:
                            nodes_to_remove_next.add(v)
                            v2coreness_t[v] = k - 1
                        v2edges[v].remove(i_e)
                    del i2edges[i_e]
            nodes_to_remove = nodes_to_remove_next.copy()
        assert len(v2coreness_t) == n
        results[t] = v2coreness_t
    e_time = time.time() - s_time
    return results, e_time

def invert_core_dict(results):
    node_to_pairs = {}
    for t, v2k in results.items():
        for v, k in v2k.items():
            if v not in node_to_pairs:
                node_to_pairs[v] = set()
            node_to_pairs[v].add((k, round(float(t),2)))

    D = {}
    for v, pairs in node_to_pairs.items():
        keep = set(pairs)
        for (k1,t1) in pairs:
            for (k2,t2) in pairs:
                if (k1, t1) == (k2, t2): continue
                if k2 <= k1 and t2 <= t1 and (k1, t1) != (k2, t2):
                    if (k2, t2) in keep:
                        keep.remove((k2, t2))
        D[v] = keep
    return D

def run(H, dataset):
    data_dir = '../cache/kt_pre'
    ds_dir = Path(data_dir) / dataset
    i2edges_path = ds_dir / "i2edges.pkl"
    v2edges_path = ds_dir / "v2edges.pkl"
    if not (i2edges_path.exists() and v2edges_path.exists()):
        preprocess_from_hypergraph(H, data_dir, dataset)

    results, e_time = coreness(ds_dir)
    D = invert_core_dict(results)
    func.save_time("kt_hypercoreness", dataset, e_time)
    return D, e_time
