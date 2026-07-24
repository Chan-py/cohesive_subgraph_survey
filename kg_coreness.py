import time
import time
from collections import defaultdict, deque

import func

def getNbrMap(hypergraph, v, g):
    neighbor_counts = {}
    for hyperedge in hypergraph.nodes[v]['hyperedges']:
        # Increment the count for each neighbor in the hyperedge
        for neighbor in hyperedge:
            if neighbor != v:
                neighbor_counts[neighbor] = neighbor_counts.get(neighbor, 0) + 1

    filtered_neighbors = {neighbor: count for neighbor, count in neighbor_counts.items() if count >= g}

    return filtered_neighbors

def reorganise(cores):
    D = {}

    for (k, g) in cores.keys():
        V = cores[(k, g)].copy()
        if (k + 1, g) in cores:
            V -= cores[(k + 1, g)]
        if (k, g + 1) in cores:
            V -= cores[(k, g + 1)]
        if V:
            D[(k, g)] = V
    return D

def decompose(hypergraph):
    D = {}
    g = 1

    while True:
        k = 0
        VQ = deque()
        queued = set()
        prev = set(hypergraph.nodes())
        C = {}
        T = defaultdict(set)
        H = set(hypergraph.nodes())

        for v in hypergraph.nodes():
            nbr_map = getNbrMap(hypergraph, v, g)
            C[v] = len(nbr_map)
            T[C[v]].add(v)

        total_nodes = {node for count, nodes in T.items() if count != 0 for node in nodes}
        if not total_nodes:
            break

        while H:
            k += 1
            for count_val in sorted(T.keys()):
                if count_val >= k:
                    break
                for node in T[count_val]:
                    if node not in queued:
                        VQ.append(node)
                        queued.add(node)
                T[count_val].clear()

            while VQ:
                current = VQ.popleft()
                queued.discard(current)
                if current not in H:
                    continue
                H.remove(current)
                del C[current]

                nbr_current = [w for w in getNbrMap(hypergraph, current, g) if w in C]
                for w in nbr_current:
                    if w in queued:
                        continue
                    old_count = C[w]
                    if w in T[old_count]:
                        T[old_count].remove(w)
                    C[w] = old_count - 1
                    new_count = C[w]
                    T[new_count].add(w)
                    if new_count < k and w not in queued:
                        VQ.append(w)
                        queued.add(w)
            if k > 1:
                D[(k - 1, g)] = prev - H
            prev = H.copy()
        g += 1

    D = reorganise(D)
    return D

def run(H, dataset):
    G = func.hypergraph_to_networkx(H)
    s_time = time.time()
    # D: {(k,g): set(nodes)}
    D = decompose(G)
    e_time = time.time() - s_time
    D = invert_core_dict(D)
    func.save_time("kg_coreness", dataset, e_time)
    return D, e_time

def invert_core_dict(D):
    node2core = {}
    for (k, g), nodes in D.items():
        for v in nodes:
            if v not in node2core:
                node2core[v] = set()
            node2core[v].add((k, g))
    return node2core