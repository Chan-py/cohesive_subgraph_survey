import time
from collections import defaultdict, Counter, deque
from math import comb
import csv
import sys
from pathlib import Path
import copy

sys.setrecursionlimit(10**7)

import func

class TNode:
    __slots__ = ("id", "children", "label", "parent", "depth", "mask", "eid", "weight", "min_hid", "max_hid", "wl_begin", "wl_end", "jump")
    def __init__(self, id, label=None, parent=None, depth=0, mask=1, eid=None, weight=0):
        self.id = id               
        self.children = {}          
        self.label = label        
        self.parent = parent      
        self.depth = depth       
        self.mask = mask             
        self.eid = eid               
        self.weight = weight      
        self.min_hid = None      
        self.max_hid = None         
        self.wl_begin = 0
        self.wl_end = 0
        self.jump = None
        
class Tree:
    def __init__(self):
        self.root = TNode(id=-1, depth=0, mask=1)
        self.root.parent = self.root

    def add_child(self, parent, key, label=None):
        if key not in parent.children:
            child = TNode(
                id=key,
                label=label,
                parent=parent,
                depth=parent.depth + 1,
                mask=(1 if parent is None else parent.mask << 1)
            )
            parent.children[key] = child
        return parent.children[key]
    
    def dfs_leaf(self, node, ancestor, leaves):
        node.jump = ancestor
        cur_min = node.eid if node.eid is not None else None
        cur_max = node.eid if node.eid is not None else None

        if node.eid is not None:
            leaves.append(node.eid)
            
        children = sorted(node.children.values(), key=lambda n: n.id)
        for i, child in enumerate(children):
            nxt_ancestor = children[i+1] if i+1 < len(children) else ancestor
            self.dfs_leaf(child, nxt_ancestor, leaves)
            
            if child.min_hid is not None:
                cur_min = child.min_hid if cur_min is None else min(cur_min, child.min_hid)
            if child.max_hid is not None:
                cur_max = child.max_hid if cur_max is None else max(cur_max, child.max_hid)

        node.min_hid = cur_min
        node.max_hid = cur_max
            
    def find_leaf(self, root=None):     
        if root is None:
            root = self.root
        leaves = []
        self.dfs_leaf(root, root, leaves)
        return leaves

    def visualize(self, root=None):
        if root is None:
            root = self.root
        lines = []
        items = list(root.children.items())
        for i, (k, child) in enumerate(items):
            is_last = (i == len(items) - 1)
            branch = "└─ " if is_last else "├─ "
            tag = f" [eid={child.eid}, w={child.weight}]" if child.eid is not None else ""
            lines.append(branch + str(k) + tag)
            self._viz(child, "   " if is_last else "│  ", lines)
        return "\n".join(lines) if lines else "(empty)"

    def _viz(self, node, prefix, lines):
        items = list(node.children.items())
        for i, (k, child) in enumerate(items):
            is_last = (i == len(items) - 1)
            branch = "└─ " if is_last else "├─ "
            tag = f" [eid={child.eid}, w={child.weight}]" if child.eid is not None else ""
            lines.append(prefix + branch + str(k) + tag)
            self._viz(child, prefix + ("   " if is_last else "│  "), lines)

def sortV(G):
    degree = {nid: len(G.nodes[nid].Edge) for nid in G.nodes}
    # print('d', degree)
    sorted_nodes = sorted(G.nodes.keys(), key=lambda x: -degree[x])
    # print(sorted_nodes)

    old_to_new = {old_id: new_id for new_id, old_id in enumerate(sorted_nodes)}

    new_nodes = {}
    for old_id, new_id in old_to_new.items():
        node_obj = G.nodes[old_id]
        node_obj.id = new_id
        new_nodes[new_id] = node_obj
    G.nodes = new_nodes

    new_hyperedges = {}
    for eid, nodeset in G.hyperedges.items():
        remapped = sorted(old_to_new[v] for v in nodeset)
        new_hyperedges[eid] = remapped
    G.hyperedges = new_hyperedges

def build_prefixForest(G):
    F = Tree()
    root = F.root
    
    for eid in sorted(G.hyperedges.keys()):        
        path = sorted(G.hyperedges[eid])
        if len(path) < 2:
            continue
        node = root
        for u in path:
            node = F.add_child(node, u)

        if node.eid is None:
            node.eid = eid
        else:
            G.weight[node.eid] += G.weight[eid]
        node.label = node.eid
    return F

def collect(n, edge2Node):
        if n.eid is not None:
            edge2Node[n.eid] = n
        for ch in n.children.values():
            collect(ch, edge2Node)

def edgeReorder(G, F):
    E_old = len(G.hyperedges)
    edgeNid = [-1] * E_old  
    edge2Node = []          

    eNid = 0
    def buildTreeDFS(node, ancestor):
        nonlocal eNid
        node.children = dict(sorted(node.children.items(), key=lambda kv: kv[0]))
        if node.eid is not None:
            old = node.eid
            node.eid = eNid
            edgeNid[old] = eNid
            edge2Node.append(node)   
            node.min_hid = node.eid
            node.max_hid = node.eid
            eNid += 1
        node.jump = ancestor
        if not node.children:
            return
        kids = list(node.children.values())
        for i in range(len(kids)-1):
            buildTreeDFS(kids[i], kids[i+1])
        buildTreeDFS(kids[-1], ancestor)
        node.min_hid = min([x for x in [node.min_hid, kids[0].min_hid] if x is not None])
        node.max_hid = max([x for x in [node.max_hid, kids[-1].max_hid] if x is not None])

    for _, ch in sorted(F.root.children.items(), key=lambda kv: kv[0]):
        buildTreeDFS(ch, F.root)

    if eNid > 0:
        F.root.min_hid = 0
        F.root.max_hid = eNid - 1
    else:
        F.root.min_hid = F.root.max_hid = None

    new_hyperedges = {}
    new_weight = {}
    for old, nodes in G.hyperedges.items():
        new = edgeNid[old]
        if new != -1:
            new_hyperedges[new] = sorted(nodes)
            new_weight[new] = G.weight.get(old, 1)
    G.hyperedges = new_hyperedges
    G.weight = new_weight

    for node_obj in G.nodes.values():
        node_obj.Edge = [edgeNid[e] for e in node_obj.Edge if edgeNid[e] != -1]
        node_obj.Edge.sort()

    E = len(G.hyperedges)
    if len(edge2Node) != E:
        for old, nodes in new_hyperedges.items():
            pass
        edge2Node = edge2Node[:E]

    old2new_map = {old: edgeNid[old] for old in range(E_old) if edgeNid[old] != -1}

    return edge2Node, old2new_map


def nodeIntersect(a, b):
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i] < b[j]:
            i += 1
        elif b[j] < a[i]:
            j += 1
        else:
            return True 
    return False

def uniqueEdges(G, old2new_eid, new2old_eid):
    norm2weight = {}
    key2hids = {}
    for hid, nodes in G.hyperedges.items():
        if len(nodes) < 2: continue
        key = tuple(sorted(nodes))
        norm2weight[key] = norm2weight.get(key, 0) + 1
        key2hids.setdefault(key, set()).add(hid)

    key_list = sorted(norm2weight.keys())
    new_hyperedges = {}
    new_weight = {}
    key2neweid = {}
    for eid, key in enumerate(key_list):
        new_hyperedges[eid] = [*key]
        new_weight[eid] = norm2weight[key]
        key2neweid[key] = eid

    old2new_eid2 = {}
    new2old_eid2 = {}
    for key, hids in key2hids.items():
        ne = key2neweid[key]
        for old in hids:
            if old in new2old_eid:
                old2new_eid2[new2old_eid[old]] = ne
                new2old_eid2[ne] = new2old_eid[old]
            else:
                old2new_eid2[old] = ne
                new2old_eid2[ne] = old

    G.hyperedges = new_hyperedges
    G.weight = new_weight

    for node in G.nodes.values():
        node.Edge = []
    for eid, nodes in G.hyperedges.items():
        for v in nodes:
            G.nodes[v].Edge.append(eid)
    for node in G.nodes.values():
        node.Edge.sort()
    return G, old2new_eid2, new2old_eid2

def core(G):
    Q = {nid for nid, node in G.nodes.items() if sum(G.weight[eid] for eid in node.Edge) < 2}
    EQ = set()
    while Q:
        n = Q.pop()
        for e in list(G.nodes[n].Edge):
            if len(G.hyperedges.get(e, set())) < 2:
                EQ.add(e)
        G.del_node(n) 

        while EQ:
            e = EQ.pop()
            if e not in G.hyperedges:
                continue
            for vn, node in G.nodes.items():
                if len(node.Edge) < 2:
                    Q.add(vn)
    return G
def get_support(G):
    triangle_count = 0
    triangle = set()
    htn = Counter()
    adj = defaultdict(Counter) 
    old2new_eid = {}
    new2old_eid = {}
    while True:
        cntV = len(G.nodes)
        cntE = len(G.hyperedges)
        G, old2new_eid, new2old_eid = uniqueEdges(G, old2new_eid, new2old_eid)
        G = core(G)
        if len(G.nodes) == cntV and len(G.hyperedges) == cntE: break

    sortV(G)
    F = build_prefixForest(G)
    edge2Node, old2new_map = edgeReorder(G, F)
    new_old_map = {old2new_map[v]:idx for idx, v in old2new_eid.items()}
    E = len(G.hyperedges)
    source = [[] for _ in range(E)]
    work_list = []
    time_stamp = [0 for _ in range(E)]
    stamp = 0

    for node in sorted(F.root.children.values(), key=lambda n: n.id):
        treeId = node.min_hid
        work_list.clear()

        while node is not F.root:
            node.wl_begin = len(work_list)
            for e in G.nodes[node.id].Edge:
                if e > node.max_hid: break
                if (not source[e]) and (e < treeId):
                    work_list.append(e)
                source[e].append(node.id)
            node.wl_end = len(work_list)

            if node.children:
                nxt = node.children[min(node.children.keys())]
            else:
                nxt = node.jump

            if node.eid is not None:
                nb = node
                while nb != F.root:
                    for wb in range(nb.wl_begin, nb.wl_end):
                        b = work_list[wb]
                        ttc = 0
                        stamp += 1
                        apt = node
                        bpt = edge2Node[b]
                        while bpt != F.root:
                            while (apt != F.root) and (apt.id > bpt.id):
                                apt = apt.parent
                            if (apt != F.root) and (apt.id == bpt.id): 
                                bpt = bpt.parent
                                continue
                            for c in G.nodes[bpt.id].Edge:
                                if c>= b: break
                                if source[c] and (time_stamp[c] != stamp):
                                    time_stamp[c] = stamp
                                    if not nodeIntersect(source[b], source[c]):
                                        ttc += G.weight[node.eid]*G.weight[b]*G.weight[c]
                                        key = tuple(sorted(list([new_old_map[node.eid], new_old_map[b], new_old_map[c]])))
                                        triangle.add(key)

                                        aa = new_old_map[node.eid]
                                        bb = new_old_map[b]
                                        cc = new_old_map[c]
                                        adj[aa][bb] += G.weight[c]
                                        adj[bb][aa] += G.weight[c]
                                        adj[aa][cc] += G.weight[b]
                                        adj[cc][aa] += G.weight[b]
                                        adj[bb][cc] += G.weight[node.eid]
                                        adj[cc][bb] += G.weight[node.eid]
                                        htn[aa] += G.weight[b]*G.weight[c]
                                        htn[bb] += G.weight[node.eid]*G.weight[c]
                                        htn[cc] += G.weight[b]*G.weight[node.eid]
                            bpt = bpt.parent
                        triangle_count += ttc
                    nb = nb.parent
                for b in range(treeId, node.eid):
                    ttc = 0
                    apt = node
                    bpt = edge2Node[b]
                    while bpt != F.root:
                        while apt.id > bpt.id:
                            apt = apt.parent
                        if apt is bpt: break
                        if apt.id == bpt.id: 
                            bpt = bpt.parent
                            continue
                        for wi in range(bpt.wl_begin, bpt.wl_end):
                            c = work_list[wi]
                            if source[c] and (not nodeIntersect(source[b], source[c])):
                                ttc += G.weight[node.eid]*G.weight[b]*G.weight[c]
                                key = tuple(sorted(list([new_old_map[node.eid], new_old_map[b], new_old_map[c]])))
                                triangle.add(key)

                                aa = new_old_map[node.eid]
                                bb = new_old_map[b]
                                cc = new_old_map[c]
                                adj[aa][bb] += G.weight[c]
                                adj[bb][aa] += G.weight[c]
                                adj[aa][cc] += G.weight[b]
                                adj[cc][aa] += G.weight[b]
                                adj[bb][cc] += G.weight[node.eid]
                                adj[cc][bb] += G.weight[node.eid]
                                htn[aa] += G.weight[b]*G.weight[c]
                                htn[bb] += G.weight[node.eid]*G.weight[c]
                                htn[cc] += G.weight[b]*G.weight[node.eid]
                                # print(a)
                                # print(b)
                                # print(c)
                                
                        bpt = bpt.parent
                    triangle_count += ttc    
            if not node.children:
                n = node
                while n is not nxt.parent:
                    for wi in range(n.wl_begin, n.wl_end):
                        source[work_list[wi]].clear()
                    n = n.parent

                n = nxt.parent
                while n is not F.root:
                    for wi in range(n.wl_begin, n.wl_end):
                        while source[work_list[wi]] and (source[work_list[wi]][-1] > nxt.parent.id):
                            source[work_list[wi]].pop()
                    n = n.parent

                if nxt is F.root:
                    for e in range(treeId, node.max_hid+1):
                        source[e].clear()
                else:
                    for e in range(treeId, nxt.min_hid):
                        while source[e] and (source[e][-1] > nxt.parent.id):
                            source[e].pop()
            node = nxt
    # print(triangle_count)
    # print(triangle)
    # print(htn)
    return htn, adj

def run(G, dataset):
    s_time = time.time()
    G1 = copy.deepcopy(G)
    htn, adj = get_support(G1)
    n = len(G.nodes)
    trussness = {}
    alive = set(G.hyperedges.keys()) 
    for hid in alive:
        if hid not in htn:
            trussness[hid] = 0
    for k in range(3, comb(n, 2)+1):
        if len(alive) == 0: break
        q = deque([e for e in list(alive) if htn.get(e, 0) < k - 2])
        while q:
            e = q.popleft()
            if e not in alive: continue

            alive.remove(e)
            if e not in trussness:
                trussness[e] = k-1

            for nb, tcnt in list(adj.get(e, {}).items()):
                if nb in alive:
                    htn[nb] -= tcnt
                    if htn[nb] < k-2:
                        q.append(nb)
                    adj[nb].pop(e, None)
            adj.pop(e, None)
            htn.pop(e, None)
    for e in alive:
        trussness[e] = k
    e_time = time.time() - s_time
    func.save_time("hyper_k_trussness", dataset, e_time)
    return trussness, e_time