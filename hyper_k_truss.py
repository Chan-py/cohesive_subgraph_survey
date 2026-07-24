import csv
import os
import time
import subprocess
from pathlib import Path

import func

# Hyper-k-truss uses the C++ HTRUSS reference implementation (paper:
# "Truss Decomposition in Hypergraphs"). We now use TrussDecomp_pf.exe: a
# drop-in with BIT-IDENTICAL per-edge trussness (verified) but a much faster
# pipeline — lazy-min-heap peeling (O(E+adj) instead of O(max_htn*E)) and
# OpenMP-parallel PushForward triangle enumeration. Same CLI/output as the
# original TrussDecomp.exe. (The old Python decomposition in
# hyper_k_trussness.py is no longer used.)
_HTRUSS = Path(__file__).resolve().parent / "HTRUSS"
_EXE    = _HTRUSS / "TrussDecomp_pf.exe"
# Cap enumeration threads so a hyper job doesn't monopolise the box when other
# datasets/jobs run concurrently under the experiment scheduler.
_OMP_THREADS = os.environ.get("HYPER_OMP_THREADS", "16")


def _ensure_exe():
    if not _EXE.exists():
        subprocess.run(
            ["g++", "-std=c++17", "-O3", "-fopenmp", "-Iinclude",
             "src/TrussDecomp_pf.cpp", "-o", "TrussDecomp_pf.exe"],
            cwd=str(_HTRUSS), check=True)


def run(network, dataset, k):
    """Return the set of hyperedge ids with hyper-truss trussness >= k.

    network : path to the dataset's network.hyp
    dataset : dataset basename (used for the trussness cache file)
    k       : truss threshold
    """
    s_time = time.time()
    csv_path = Path(f"../output/{dataset}_hyper_k_trussness.csv")
    decomp_time = func.read_decomp_from_csv(dataset, "hyper_k_trussness")
    e_dec = 0
    if not decomp_time or not csv_path.exists():
        _ensure_exe()
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        dec = time.time()
        # C++ writes edge_id,trussness for every original hyperedge id.
        _env = dict(os.environ, OMP_NUM_THREADS=_OMP_THREADS)
        subprocess.run([str(_EXE), str(network), str(csv_path)],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, env=_env)
        decomp_time = time.time() - dec
        func.save_time("hyper_k_trussness", dataset, decomp_time)
        e_dec = time.time() - dec

    HE = set()
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if int(row["trussness"]) >= k:
                HE.add(int(row["edge_id"]))

    e_time = (time.time() - s_time) + decomp_time - e_dec
    return HE, e_time
