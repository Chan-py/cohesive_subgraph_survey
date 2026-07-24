#!/usr/bin/env python3
"""kinout_truss C++ vs Python 검증.
- 소규모 랜덤 하이퍼그래프 여러 개 생성(고정 시드) + 여러 (k_in,k_out) 조합.
- Python kinout_truss.run() 와 C++ ./kinout_truss 의 (nodes, hids) 집합이 정확히 일치하는지 비교.
반드시 code/ 디렉토리에서 실행 (func 의 상대경로 캐시 사용). 메모리 안전: 전부 소규모.
"""
import os, sys, random, subprocess, shutil
from pathlib import Path

CODE = Path("/home/seungchan/cohesive_survey/code")
sys.path.insert(0, str(CODE))
import func
import kinout_truss

HERE = Path(__file__).resolve().parent
TESTD = HERE / "tests"
EXE  = HERE / "kinout_truss"
CACHE = Path("/home/seungchan/cohesive_survey/cache")

PARAMS = [(0,0),(1,1),(2,2),(2,3),(2,6),(3,3),(3,5),(4,3),(4,6),(5,4),(6,10),(1,10)]


def gen_graph(name, n, m, cmin, cmax, seed):
    random.seed(seed)
    d = TESTD / name
    d.mkdir(parents=True, exist_ok=True)
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


def py_run(net, k_in, k_out):
    name = Path(net).parent.name
    clear_cache(name)                     # fresh state (peel mutates H/G)
    H = func.Hypergraph()
    H.load_from_file(str(net))
    nodes, hids, _ = kinout_truss.run(H, str(net), k_in, k_out)
    return set(int(x) for x in nodes), set(int(x) for x in hids)


def cpp_run(net, k_in, k_out):
    out = HERE / "tests" / "_cpp_out.txt"
    subprocess.run([str(EXE), str(net), str(k_in), str(k_out), str(out)], check=True)
    txt = out.read_text().split("\n")
    nodes = set(int(x) for x in txt[0].split()) if len(txt) > 0 and txt[0].strip() else set()
    hids  = set(int(x) for x in txt[1].split()) if len(txt) > 1 and txt[1].strip() else set()
    return nodes, hids


def main():
    graphs = [
        gen_graph("vtest_ki_01", 20,  40, 2, 4, 1),
        gen_graph("vtest_ki_02", 30,  80, 2, 5, 2),
        gen_graph("vtest_ki_03", 15,  60, 3, 6, 3),   # dense small
        gen_graph("vtest_ki_04", 50, 120, 2, 3, 4),   # sparse-ish
        gen_graph("vtest_ki_05", 25, 100, 2, 6, 5),
        gen_graph("vtest_ki_06", 40,  50, 4, 6, 6),   # high cardinality
        gen_graph("vtest_ki_07", 12,  30, 2, 5, 7),   # tiny
        gen_graph("vtest_ki_08", 60, 150, 2, 4, 8),
    ]
    # real small dataset (kinout cache here is unused by running experiments)
    real_contact = Path("/home/seungchan/cohesive_survey/datasets/real/contact/network.hyp")
    if real_contact.exists():
        graphs.append(real_contact)

    total = fail = 0
    nonempty = 0
    for net in graphs:
        gname = Path(net).parent.name
        for (ki, ko) in PARAMS:
            total += 1
            pn, ph = py_run(net, ki, ko)
            cn, ch = cpp_run(net, ki, ko)
            ok = (pn == cn and ph == ch)
            if pn or ph:
                nonempty += 1
            if not ok:
                fail += 1
                print(f"[FAIL] {gname:14s} k_in={ki} k_out={ko}")
                print(f"       nodes  py={len(pn)} cpp={len(cn)}  only_py={sorted(pn-cn)[:8]} only_cpp={sorted(cn-pn)[:8]}")
                print(f"       hids   py={len(ph)} cpp={len(ch)}  only_py={sorted(ph-ch)[:8]} only_cpp={sorted(ch-ph)[:8]}")
            else:
                print(f"[ ok ] {gname:14s} k_in={ki} k_out={ko}  nodes={len(pn):4d} hids={len(ph):4d}")

    print("=" * 60)
    print(f"total={total}  pass={total-fail}  fail={fail}  (non-empty results={nonempty})")
    # cleanup test caches
    for net in graphs:
        clear_cache(Path(net).parent.name)
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
