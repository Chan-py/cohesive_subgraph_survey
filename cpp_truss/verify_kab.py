#!/usr/bin/env python3
"""kab_truss C++ vs Python 검증 (a>0 경로).
소규모 하이퍼그래프(정점 적게 → triple 반복 → valid triangle 발생) × 여러 (k,a,b).
Python kab_truss.run() 와 C++ ./kab_truss 의 (nodes, hids) 집합 일치 확인.
code/ 에서 실행. 전부 소규모 → 메모리 안전.
"""
import os, sys, random, subprocess, shutil
from pathlib import Path

CODE = Path("/home/seungchan/cohesive_survey/code")
sys.path.insert(0, str(CODE))
import func
import kab_truss

HERE = Path(__file__).resolve().parent
TESTD = HERE / "tests"
EXE  = HERE / "kab_truss"
CACHE = Path("/home/seungchan/cohesive_survey/cache")

# (k, a, b) — real groups are k=3 with (2,3)(4,5)(6,7)(8,9); plus extras. a>=1.
PARAMS = [(3,2,3),(3,4,5),(3,6,7),(3,8,9),
          (1,1,1),(2,1,2),(2,2,2),(1,2,4),(2,3,5),(3,1,1),(1,5,10),(2,2,8)]


def gen_graph(name, n, m, cmin, cmax, seed):
    random.seed(seed)
    d = TESTD / name; d.mkdir(parents=True, exist_ok=True)
    lines = []
    for _ in range(m):
        c = random.randint(cmin, cmax)
        verts = random.sample(range(1, n + 1), min(c, n))
        lines.append(",".join(map(str, verts)))
    (d / "network.hyp").write_text("\n".join(lines) + "\n")
    return d / "network.hyp"


def clear_cache(name):
    for sub in ("kab_clique", "kinout_support", "kab_proj"):
        p = CACHE / sub / name
        if p.exists():
            shutil.rmtree(p, ignore_errors=True)


def py_run(net, k, a, b):
    name = Path(net).parent.name
    clear_cache(name)
    H = func.Hypergraph(); H.load_from_file(str(net))
    nodes, hids, _ = kab_truss.run(H, k, a, b, str(net))
    return set(int(x) for x in nodes), set(int(x) for x in hids)


def cpp_run(net, k, a, b):
    out = HERE / "tests" / "_cpp_kab_out.txt"
    subprocess.run([str(EXE), str(net), str(k), str(a), str(b), str(out)], check=True)
    t = out.read_text().split("\n")
    nodes = set(int(x) for x in t[0].split()) if len(t) > 0 and t[0].strip() else set()
    hids  = set(int(x) for x in t[1].split()) if len(t) > 1 and t[1].strip() else set()
    return nodes, hids


def main():
    graphs = [
        gen_graph("vtest_kab_01", 10,  60, 3, 5, 11),   # small n, dense triples
        gen_graph("vtest_kab_02",  8,  80, 3, 6, 12),
        gen_graph("vtest_kab_03", 14, 100, 3, 5, 13),
        gen_graph("vtest_kab_04", 12,  70, 4, 6, 14),   # high cardinality
        gen_graph("vtest_kab_05", 16, 120, 3, 4, 15),
        gen_graph("vtest_kab_06",  9,  50, 3, 7, 16),
        gen_graph("vtest_kab_07", 20, 150, 3, 6, 17),
        gen_graph("vtest_kab_08",  6,  40, 3, 5, 18),   # tiny, very dense
    ]
    real_contact = Path("/home/seungchan/cohesive_survey/datasets/real/contact/network.hyp")
    if real_contact.exists():
        graphs.append(real_contact)

    total = fail = nonempty = 0
    for net in graphs:
        gname = Path(net).parent.name
        for (k, a, b) in PARAMS:
            total += 1
            pn, ph = py_run(net, k, a, b)
            cn, ch = cpp_run(net, k, a, b)
            ok = (pn == cn and ph == ch)
            if pn or ph:
                nonempty += 1
            if not ok:
                fail += 1
                print(f"[FAIL] {gname:14s} k={k} a={a} b={b}")
                print(f"       nodes py={len(pn)} cpp={len(cn)} only_py={sorted(pn-cn)[:8]} only_cpp={sorted(cn-pn)[:8]}")
                print(f"       hids  py={len(ph)} cpp={len(ch)} only_py={sorted(ph-ch)[:8]} only_cpp={sorted(ch-pn)[:8]}")
            else:
                print(f"[ ok ] {gname:14s} k={k} a={a} b={b}  nodes={len(pn):4d} hids={len(ph):4d}")

    print("=" * 60)
    print(f"total={total}  pass={total-fail}  fail={fail}  (non-empty results={nonempty})")
    for net in graphs:
        clear_cache(Path(net).parent.name)
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
