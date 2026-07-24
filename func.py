from sortedcontainers import SortedDict
import re
import networkx as nx
import time
from itertools import combinations
import os
import ast
from collections import defaultdict
import pickle
from pathlib import Path
import csv, sys

def _bump_csv_limit():
    max_int = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_int)
            break
        except OverflowError:
            max_int //= 10

# 전처리 캐시 루트 (datasets/ 밖에 분리 보관)
CACHE_DIR = "../cache"

def cache_path(network, cache_name):
    """원본 network.hyp 경로 -> 해당 데이터셋의 캐시 network.hyp 경로.

    예) ../datasets/synthetic/synthetic_p2_n10k_r0/network.hyp
        -> ../cache/kab_proj/synthetic_p2_n10k_r0/network.hyp
    """
    name = Path(network).parent.name        # 데이터셋 폴더명 (type 폴더는 제외)
    return str(Path(CACHE_DIR) / cache_name / name / "network.hyp")

class Node:
    def __init__(self, node_id):
        self.id = node_id
        self.NodeCnt = 0 
        self.EdgeCnt = 0
        self.Edge = set()

class Hypergraph:
    def __init__(self):
        self.nodes = {}
        self.hyperedges = {}
        self.weight = {}

    def add_hyperedge(self, edge_nodes):
        if len(edge_nodes) < 2: return
        hyperedge_id = len(self.hyperedges) + 1
        self.hyperedges[hyperedge_id] = edge_nodes

        for node in edge_nodes:
            if node not in self.nodes:
                self.nodes[node] = Node(node)
            self.nodes[node].Edge.add(hyperedge_id)

    def load_from_file(self, file_path):
        with open(file_path, 'r') as file:
            for line in file:
                current_nodes = {int(node.strip()) for node in re.split(r'[,\s]+', line.strip())}
                if len(current_nodes) > 1:
                    self.add_hyperedge(current_nodes)

    def del_node(self, node):
        if node in self.nodes:
            for hyperedge in self.nodes[node].Edge:
                self.hyperedges[hyperedge].remove(node)
            del self.nodes[node]
    
    def del_node_in_edge(self, node, edge):
        if node in self.nodes:
            if edge in self.nodes[node].Edge:
                self.hyperedges[edge].remove(node)
                self.nodes[node].Edge.remove(edge)
                return True
            else:
                return False
        return False

    def del_edge(self, edge):
        if edge in self.hyperedges:
            for node in self.hyperedges[edge]:
                self.nodes[node].Edge.remove(edge)
            del self.hyperedges[edge]
            return True
        else: return False

    def trans_bipartite(self):
        B = nx.Graph()
        node_ids = {f'n{key}': key for key in self.nodes.keys()}
        hyperedge_ids = {f'h{key}': key for key in self.hyperedges.keys()}
        B.add_nodes_from(hyperedge_ids, bipartite=0)
        B.add_nodes_from(node_ids, bipartite=1)
        edge_list = [(f'h{id}', f'n{node}') for id, edge_nodes in self.hyperedges.items() for node in edge_nodes]
        B.add_edges_from(edge_list)
        return B

    def trans_clique(self):
        G = nx.Graph()
        for edge_nodes in self.hyperedges.values():
            for u, v in combinations(edge_nodes, 2):
                G.add_edge(u, v)
        return G
    

class PriorityQueue:
    def __init__(self):
        self.data = SortedDict()
        self.node_to_priority = {}

    def push(self, node, priority):
        if priority not in self.data:
            self.data[priority] = set()
        self.data[priority].add(node)
        self.node_to_priority[node] = priority

    def pop(self):
        if not self.data:
            raise IndexError("pop from empty priority queue")
        highest_priority = self.data.peekitem(-1)[0]
        node = self.data[highest_priority].pop()
        del self.node_to_priority[node]
        if not self.data[highest_priority]:
            del self.data[highest_priority]
        return node

    def remove(self, node):
        priority = self.node_to_priority[node]
        nodes = self.data[priority]
        nodes.remove(node)
        if not nodes:
            del self.data[priority]
        del self.node_to_priority[node]

    def empty(self):
        return len(self.data) == 0

    def contains(self, node):
        return node in self.node_to_priority

def getgNbrMap(G, node, g):
    cnt = {}
    for hyperedge in G.nodes[node].Edge:
        for neighbor in G.hyperedges[hyperedge]:
            if neighbor != node:
                if neighbor not in cnt:
                    cnt[neighbor] = 0
                cnt[neighbor] += 1
    ng = {node: count for node, count in cnt.items() if count >= g}
    return ng

def get_partially_subgraph(G, node_set, min_size=2):
    H = Hypergraph()
    S = set(node_set)
    for edge_nodes in G.hyperedges.values():
        sub = edge_nodes & S
        if len(sub) >= min_size:
            H.add_hyperedge(sub)
    return H

def prune_small_edges(G, min_size=2):
    # Def 8 floor: a partially-induced hyperedge must retain >= min_size members.
    H = Hypergraph()
    for edge_nodes in G.hyperedges.values():
        if len(edge_nodes) >= min_size:
            H.add_hyperedge(edge_nodes)
    return H

def get_strongly_subgraph(G, node_set):
    H = Hypergraph()
    S = set(node_set)
    for edge_nodes in G.hyperedges.values():
        if edge_nodes.issubset(S):
            H.add_hyperedge(edge_nodes)
    return H

def get_edge_subgraph(G, HE):
    H = Hypergraph()
    for hid, edge_nodes in G.hyperedges.items():
        if hid in HE:
            H.add_hyperedge(edge_nodes)
    return H

def get_node_edge_supgraph(G, node_set, HE):
    H = Hypergraph()
    S = set(node_set)
    for hid, edge_nodes in G.hyperedges.items():
        if hid in HE:
            sub = edge_nodes & S
            if len(sub) > 1:
                H.add_hyperedge(sub)
    return H

def hypergraph_to_networkx(H):
    G = nx.Graph()
    for _, nodes in H.hyperedges.items():
        hedge = set(nodes)
        for v in hedge:
            if v not in G:
                G.add_node(v, hyperedges = [])
            G.nodes[v]['hyperedges'].append(hedge)
    return G

def get_common_hyperedges(u, v, H):
    comm = H.nodes[u].Edge & H.nodes[v].Edge
    return comm

def parse_hids(s):
    # 문자열 "{1,2,3}" / "{  }" / "" / None 모두 안전 처리
    if not isinstance(s, str):
        return set()
    s = s.strip()
    if s.startswith("{") and s.endswith("}"):
        s = s[1:-1].strip()
    if not s:
        return set()
    return set(map(int, filter(None, (x.strip() for x in s.split(",")))))

def load_projection_graph(H, network):
    _bump_csv_limit()
    rp = cache_path(network, 'kab_proj')
    if not os.path.exists(rp):
        os.makedirs(os.path.dirname(rp), exist_ok=True)
        G, proj_time, id2pair, pair2id = build_projection_graph(H.hyperedges, rp)
        with open(rp.replace("network.hyp", "time.txt"), "w") as f:
            f.write(str(proj_time))
    else:
        G = nx.Graph()
        with open(rp.replace("network.hyp", "omega.csv"), "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                node = int(row["node"])        
                weight = int(row["omega"])
                hids = parse_hids(row['hids'])
                if node not in G:
                    G.add_node(node)
                G.nodes[node]["omega"] = weight
                G.nodes[node]["hids"] = hids
            G.add_edges_from(nx.read_edgelist(rp, nodetype=int, data=(("sigma", int),)).edges(data=True))


        id2pair = {}
        pair2id = {}
        with open(rp.replace("network.hyp", "id2pair.csv"), "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                id2pair[int(row['id'])] = tuple(ast.literal_eval(row['pair']))
                pair2id[tuple(ast.literal_eval(row['pair']))] = int(row['id'])
        with open(rp.replace("network.hyp", "time.txt"), "r") as f:
            proj_time = float(f.readline().strip())
    return G, proj_time, id2pair, pair2id

def build_projection_graph(hyperedges, path):
    s_time = time.time()
    G = nx.Graph()
    omega = defaultdict(int)
    hids = defaultdict(set)
    id = 0
    id2pair = {}
    pair2id = {}
    for hid, hedge in hyperedges.items():
        if len(hedge) < 2: continue

        # Add nodes
        arr = sorted(hedge)
        nodes = []
        for node in combinations(arr, 2):
            node = tuple(sorted(node))
            if node in pair2id:
                pid = pair2id[node]
            else: 
                id += 1
                pid = id
                pair2id[node] = pid
                id2pair[pid] = node
            # Update omega
            omega[pid] += 1
            hids[pid].add(hid)
            # if pid in node:
            #     if pid == 2:
            #         print("pp", hids[pid])
            # if pid ==2:
            #     print('tmp', hids[pid])
            if pid not in G:                
                G.add_node(pid, omega=omega[pid], hids=hids[pid])
            else:
                G.nodes[pid]['omega'] = omega[pid]
                G.nodes[pid]['hids'] = hids[pid]
            nodes.append(pid)
        
        # Update edges
        for node1, node2 in combinations(nodes, 2):
            node1, node2 = min(node1, node2), max(node1, node2)
            if len(set(id2pair[node1] + id2pair[node2])) != 3: continue

            if G.has_edge(node1, node2):
                G[node1][node2]['sigma'] += 1
            else:
                G.add_edge(node1, node2, sigma=1)
    # for node in G.nodes:
    #     print(node, G.nodes[node]['hids'])
    # for hi in hids:
    #     print(hi, hids[hi])
    e_time = time.time() - s_time
    # exit()
    save_dict2csv(path.replace('network.hyp', 'id2pair.csv'), ["id", "pair"], id2pair)
    save_dict2csv2(path.replace('network.hyp', 'omega.csv'), ["node", "omega", "hids"], omega, hids)
    nx.write_edgelist(G, path, data=['sigma'])
    return G, e_time, id2pair, pair2id

def save_dict2csv(path, header, dict):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header) 
        for a, b in dict.items():
            writer.writerow([a, b])

def save_dict2csv2(path, header, dict1, dict2):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header) 
        for a, b in dict1.items():
            writer.writerow([a, b, dict2[a]])

def build_clique_graph(hyperedges):
    s_time = time.time()
    G = nx.Graph()

    for hid, nodes in hyperedges.items():
        arr = sorted(nodes)
        if len(arr) < 2:
            continue
        for u, v in combinations(arr, 2):
            if G.has_edge(u, v):
                G[u][v]["hids"].add(hid)
            else:
                G.add_edge(u, v, hids={hid})

    e_time = time.time() - s_time
    return G, e_time

def load_clique_graph(H, network):
    rp = cache_path(network, 'kab_clique')
    if not os.path.exists(rp):
        os.makedirs(os.path.dirname(rp), exist_ok=True)
        G, proj_time = build_clique_graph(H.hyperedges)
        with open(rp.replace('network.hyp', 'network.pkl'), "wb") as f:
            pickle.dump(G, f) 
        with open(rp.replace("network.hyp", "time.txt"), "w") as f:
            f.write(str(proj_time))
    else:
        with open(rp.replace('network.hyp', 'network.pkl'), "rb") as f:
            G = pickle.load(f) 
        with open(rp.replace("network.hyp", "time.txt"), "r") as f:
            proj_time = float(f.readline().strip())
    return G, proj_time

def read_decomp_from_csv(dataset_name, algorithm):
    csv_path = "../output/time.csv"
    latest_time = None

    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        # Expect header: algorithm,dataset,time
        for row in reader:
            if not row:
                continue
            algo = (row.get("algorithm") or "").strip()
            ds   = (row.get("dataset") or "").strip()
            tval = (row.get("time") or "").strip()
            if algo == algorithm and ds == dataset_name:
                # keep last matching (file is append-only; last is most recent)
                latest_time = tval

    if latest_time is None:
        return False

    try:
        return float(latest_time)
    except ValueError:
        raise ValueError(f"Invalid numeric value for time in CSV (algorithm='{algorithm}', dataset='{dataset_name}'): '{latest_time}'")

def initial_inout_support(G):
    s_time = time.time()
    sup_in = defaultdict(int)
    sup_out = defaultdict(int)

    for u, v in list(G.edges()):
        u, v = get_key(u,v)
        Nu = set(G.neighbors(u))
        Nv = set(G.neighbors(v))
        com_nb = Nu & Nv
        for e in G[u][v]['hids']:
            for w in com_nb:
                huw = G[u][w]['hids']
                hvw = G[v][w]['hids']
                inn = 0
                if (e in huw) and (e in hvw):
                    sup_in[(e, u, v)] += 1
                    inn = 1
                sup_out[(e, u, v)] += (len(huw)*len(hvw)-inn)
    sup_node_out = defaultdict(int)
    for (e, u, v), sup in sup_out.items():
        if sup_node_out[(e, u)] < sup:
            sup_node_out[(e, u)] = sup
        if sup_node_out[(e, v)] < sup:
            sup_node_out[(e, v)] = sup
    e_time = time.time() - s_time
    return sup_in, sup_out, sup_node_out, e_time

def get_key(u, v):
    if u > v:
        return v, u
    return u, v


def load_support_inout(G, network):
    rp = cache_path(network, 'kinout_support')
    if not os.path.exists(rp):
        os.makedirs(os.path.dirname(rp), exist_ok=True)
        sup_in, sup_out, sup_node_out, sup_time = initial_inout_support(G)
        # --- sup_in.csv ---
        with open(rp.replace('network.hyp', 'sup_in.csv'), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["e", "u", "v", "sup_in"])
            for (e, u, v), val in sup_in.items():
                writer.writerow([e, u, v, val])

        # --- sup_out.csv ---
        with open(rp.replace('network.hyp', 'sup_out.csv'), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["e", "u", "v", "sup_out"])
            for (e, u, v), val in sup_out.items():
                writer.writerow([e, u, v, val])

        # --- sup_node_out.csv ---
        with open(rp.replace('network.hyp', 'sup_node_out.csv'), "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["e", "u", "sup_node_out"])
            for (e, u), val in sup_node_out.items():
                writer.writerow([e, u, val])
            with open(rp.replace("network.hyp", "time.txt"), "w") as f:
                f.write(str(sup_time))
    else:
        sup_in = defaultdict(int)
        sup_out = defaultdict(int)
        sup_node_out = defaultdict(int)

        # --- sup_in.csv ---
        with open(rp.replace('network.hyp', 'sup_in.csv'), "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sup_in[(int(row["e"]), int(row["u"]), int(row["v"]))] = int(row["sup_in"])

        # --- sup_out.csv ---
        with open(rp.replace('network.hyp', 'sup_out.csv'), "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                e = int(row["e"])
                u = int(row["u"])
                v = int(row["v"])
                sup_out[(int(row["e"]), int(row["u"]), int(row["v"]))] = int(row["sup_out"])

        # --- sup_node_out.csv ---
        with open(rp.replace('network.hyp', 'sup_node_out.csv'), "r", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sup_node_out[(int(row["e"]), int(row["u"]))] = int(row["sup_node_out"])

        # --- time.txt ---
        time_path = rp.replace("network.hyp", "time.txt")
        sup_time = 0.0
        if os.path.exists(time_path):
            with open(time_path, "r") as f:
                sup_time = float(f.read().strip() or 0)
    return sup_in, sup_out, sup_node_out, sup_time

def update_support_inout(G, sup_in, sup_out, sup_node_out, affected_nodes):
    # print(affected_nodes)/
    if not affected_nodes:
        return sup_in, sup_out, sup_node_out

    affected_nodes = set(affected_nodes)

    # 1) 영향 쌍(엣지) 수집
    affected_pairs = set()

    # 1-1) affected node에 '붙은' 엣지들
    for u in list(affected_nodes):
        if not G.has_node(u):
            continue
        for v in G.neighbors(u):
            if G.has_edge(u, v):
                a, b = (u, v) if u < v else (v, u)
                affected_pairs.add((a, b))

    for w in list(affected_nodes):
        if not G.has_node(w):
            continue
        nbrs = [x for x in G.neighbors(w) if G.has_node(x)]
        for x, y in combinations(nbrs, 2):
            a, b = (x, y) if x < y else (y, x)
            if G.has_edge(a, b):
                affected_pairs.add((a, b))

    def _keys_for_pair_in(d, u, v):
        return [(e, a, b) for (e, a, b) in d.keys() if (a, b) == (u, v)]

    for (u, v) in list(affected_pairs):
        if not (G.has_node(u) and G.has_node(v) and G.has_edge(u, v)):
            for key in _keys_for_pair_in(sup_in, u, v):
                sup_in.pop(key, None)
                sup_out.pop(key, None)
            for key in _keys_for_pair_in(sup_out, u, v):
                sup_in.pop(key, None)
                sup_out.pop(key, None)
            continue

        Huv = G[u][v].get('hids', set())

        for (e, a, b) in _keys_for_pair_in(sup_in, u, v):
            if e not in Huv:
                sup_in.pop((e, a, b), None)
                sup_out.pop((e, a, b), None)
        for (e, a, b) in _keys_for_pair_in(sup_out, u, v):
            if e not in Huv:
                sup_in.pop((e, a, b), None)
                sup_out.pop((e, a, b), None)
            
        Nu, Nv = G[u], G[v]
        if len(Nu) > len(Nv):
            uu, vv = v, u
            Nu, Nv = Nv, Nu
        else:
            uu, vv = u, v
        com_nb = [w for w in Nu if w in Nv]

        # 2-3) 현재 살아있는 e만 재계산
        for e in list(Huv):
            cnt_in = 0
            sum_out = 0
            for w in com_nb:
                if not (G.has_edge(uu, w) and G.has_edge(vv, w)):
                    continue
                Huw = G[uu][w].get('hids', set())
                Hvw = G[vv][w].get('hids', set())
                inn = (e in Huw) and (e in Hvw)
                if inn:
                    cnt_in += 1
                # 내부 1개 제외
                sum_out += (len(Huw) * len(Hvw) - (1 if inn else 0))
            if cnt_in == 0 and sum_out == 0:    
                sup_in.pop((e, a, b), None)
                sup_out.pop((e, a, b), None)
            else: 
                key = get_key(uu, vv)
                sup_in[(e, key[0], key[1])] = cnt_in
                sup_out[(e, key[0], key[1])] = sum_out
    

    sup_node_out = defaultdict(int)
    for (e, u, v), sup in sup_out.items():
        if sup_node_out[(e, u)] < sup:
            sup_node_out[(e, u)] = sup
        if sup_node_out[(e, v)] < sup:
            sup_node_out[(e, v)] = sup

    return sup_in, sup_out, sup_node_out
    
def save_time(algorithm, dataset, e_time):
    if dataset[:2] == "ex": return
    # 경로 설정
    output_dir = Path("../output")
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "time.csv"

    # 파일 존재 여부 확인
    file_exists = csv_path.exists()

    # CSV 작성
    with open(csv_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        # 파일이 없으면 헤더 추가
        if not file_exists:
            writer.writerow(["algorithm", "dataset", "time"])

        # 새 행 추가
        writer.writerow([algorithm, dataset, e_time])