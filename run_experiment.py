#!/usr/bin/env python3
"""
Hypergraph cohesive subgraph benchmark — single-job runner.

Runs exactly one (dataset, model, param_group) experiment and appends
one row to the output CSV.  All looping and parallelism is handled by
the calling shell script.

Usage (always run inside myenv):
    conda run -n myenv python run_experiment.py \\
        --eq EQ4 \\
        --dataset real/walmart \\
        --model kinout_truss_kout3 \\
        --param-group q2 \\
        --output ../results/eq3_real_tmp/walmart_kinout_truss_kout3_q2.csv
"""

import argparse
import copy
import time
import csv
import os
import traceback
import threading
import multiprocessing as mp
import psutil
import resource
from pathlib import Path

import func
import eval as ev

import k_hypercore
import nbr_kd_core
import kt_hypercore
import kq_core
import kg_core
import kgp_core
import ks_core
import kab_truss
import kinout_truss
import hyper_k_truss
import nbr_kd_coreness
import cpp_truss_wrap   # C++ ports of kab_truss / kinout_truss (drop-in, verified)

# ─────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────

DATASETS_DIR = Path("../datasets")
RESULTS_DIR  = Path("../results")

# ─────────────────────────────────────────────
# Model parameter table
# ─────────────────────────────────────────────

# Unified single-row design (see hypergraph_experiment_plan.md §3.1):
#   models with an extra parameter fix k=4 and sweep that parameter over q1..q4;
#   pure-k models (k_hypercore, nbr_k_core, hyper_k_truss) keep their k sweep;
#   truss models sweep both new params (kab diagonal b=a+1; kinout 2x2 grid).
MODEL_PARAMS = {
    "k_hypercore":   {"q1":{"k":2},    "q2":{"k":4},    "q3":{"k":6},    "q4":{"k":8}},
    "nbr_k_core":    {"q1":{"k":2},    "q2":{"k":4},    "q3":{"k":6},    "q4":{"k":8}},
    "hyper_k_truss": {"q1":{"k":4},    "q2":{"k":6},    "q3":{"k":8},    "q4":{"k":10}},

    "kd_core":       {"q1":{"k":4,"d":2}, "q2":{"k":4,"d":4}, "q3":{"k":4,"d":6}, "q4":{"k":4,"d":8}},
    # second (k,d)-core variant: d fixed at 4, sweep k (mirrors the pure-k grid)
    "kd_core_ksweep":{"q1":{"k":2,"d":4}, "q2":{"k":4,"d":4}, "q3":{"k":6,"d":4}, "q4":{"k":8,"d":4}},
    "kt_hypercore":  {"q1":{"k":4,"t":0.2},"q2":{"k":4,"t":0.4},"q3":{"k":4,"t":0.6},"q4":{"k":4,"t":0.8}},
    "kq_core":       {"q1":{"k":4,"q":2}, "q2":{"k":4,"q":4}, "q3":{"k":4,"q":6}, "q4":{"k":4,"q":8}},
    "kg_core":       {"q1":{"k":4,"g":2}, "q2":{"k":4,"g":4}, "q3":{"k":4,"g":6}, "q4":{"k":4,"g":8}},
    "kgp_core":      {"q1":{"k":4,"g":4,"p":0.2},"q2":{"k":4,"g":4,"p":0.4},"q3":{"k":4,"g":4,"p":0.6},"q4":{"k":4,"g":4,"p":0.8}},
    "ks_core":       {"q1":{"k":4,"s":0.2,"c":1.0},"q2":{"k":4,"s":0.4,"c":1.0},"q3":{"k":4,"s":0.6,"c":1.0},"q4":{"k":4,"s":0.8,"c":1.0}},

    "kab_truss":     {"q1":{"k":3,"a":2,"b":3},"q2":{"k":3,"a":4,"b":5},"q3":{"k":3,"a":6,"b":7},"q4":{"k":3,"a":8,"b":9}},
    "kinout_truss":  {"q1":{"k_in":2,"k_out":3},"q2":{"k_in":2,"k_out":6},"q3":{"k_in":4,"k_out":3},"q4":{"k_in":4,"k_out":6}},
}

_ALGO_BASE = {
    "kd_core_ksweep":      "kd_core",
    "kd_core_d4":          "kd_core",
    "kd_core_d6":          "kd_core",
    "kt_hypercore_k4":     "kt_hypercore",
    "kt_hypercore_k6":     "kt_hypercore",
    "kq_core_q4":          "kq_core",
    "kq_core_q6":          "kq_core",
    "kg_core_g4":          "kg_core",
    "kg_core_g6":          "kg_core",
    "kgp_core_g4":         "kgp_core",
    "kgp_core_g6":         "kgp_core",
    "ks_core_k4":          "ks_core",
    "ks_core_k6":          "ks_core",
    "kab_truss_b3":        "kab_truss",
    "kab_truss_b6":        "kab_truss",
    "kinout_truss_kout3":  "kinout_truss",
    "kinout_truss_kout6":  "kinout_truss",
}

# CSV column order
BASE_COLS = [
    "experiment_id", "eq", "dataset_type", "dataset_name",
    "generator_setting", "n", "model", "parameter_group", "parameters",
    "repeat_id", "seed",
    "runtime_sec_total", "runtime_sec_model_reported", "peak_memory_mb",
]
EFFECTIVENESS_COLS = [
    "density", "density_ratio",
    "avg_degree", "avg_num_neighbors",
    "avg_support", "modularity",
    "conductance", "normalized_cut",
    "avg_connected_component_size",
    "avg_dcs",
    "num_nodes", "num_hyperedges", "avg_cardinality",
]
TAIL_COLS = ["status", "error_message"]
ALL_COLS  = BASE_COLS + EFFECTIVENESS_COLS + TAIL_COLS

JOB_TIMEOUT_SEC     = int(os.environ.get("JOB_TIMEOUT_SEC", 24 * 3600))
JOB_MEMORY_LIMIT_MB = int(os.environ.get("JOB_MEMORY_LIMIT_MB", 50_000))
# Kill the job if *system-wide* memory usage exceeds this percentage of total RAM.
JOB_MEMORY_LIMIT_PERCENT = float(os.environ.get("JOB_MEMORY_LIMIT_PERCENT", 90))


# ─────────────────────────────────────────────
# Peak memory tracker
# ─────────────────────────────────────────────

class PeakMemoryTracker:
    def __init__(self):
        self._proc   = psutil.Process(os.getpid())
        self._peak   = 0.0
        self._stop   = False
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self._stop:
            try:
                rss = self._proc.memory_info().rss
                # include child processes (e.g. the C++ TrussDecomp subprocess)
                for c in self._proc.children(recursive=True):
                    try:
                        rss += c.memory_info().rss
                    except Exception:
                        pass
                rss = rss / 1024 / 1024
                if rss > self._peak:
                    self._peak = rss
            except Exception:
                pass
            time.sleep(0.05)

    def __enter__(self):
        self._peak = self._proc.memory_info().rss / 1024 / 1024
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop = True
        self._thread.join()

    @property
    def peak_mb(self):
        return self._peak


# ─────────────────────────────────────────────
# Model runner
# ─────────────────────────────────────────────

def _run_model(algorithm, G_orig, params, network_path, dataset_name):
    algorithm = _ALGO_BASE.get(algorithm, algorithm)
    G = copy.deepcopy(G_orig)

    if algorithm == "k_hypercore":
        C, rt = k_hypercore.run(network_path, dataset_name, params["k"])
        GC = func.get_strongly_subgraph(G_orig, C)

    elif algorithm == "nbr_k_core":
        C, rt = nbr_kd_core.run(network_path, dataset_name, params["k"], 0, "nbr_k_core")
        GC = func.get_strongly_subgraph(G_orig, C)

    elif algorithm == "kd_core":
        t0 = time.time()
        GC = nbr_kd_core.kdCore(G, params["k"], params["d"])
        rt = time.time() - t0

    elif algorithm == "kt_hypercore":
        C, HE, rt = kt_hypercore.run(G, params["k"], params["t"], dataset_name)
        GC = func.get_node_edge_supgraph(G_orig, C, HE)

    elif algorithm == "kq_core":
        C, rt = kq_core.run(G, params["k"], params["q"])
        # Def 13: every retained hyperedge must keep >= q members. The peeled
        # edge set equals {e & V' : |e & V'| >= q}, so rebuild with floor q.
        GC = func.get_partially_subgraph(G_orig, C, min_size=max(2, params["q"]))

    elif algorithm == "kg_core":
        result, _, rt = kg_core.run(G, params["k"], params["g"])
        GC = func.get_partially_subgraph(G_orig, result.nodes)

    elif algorithm == "kgp_core":
        GC, rt = kgp_core.run(G, params["k"], params["g"], params["p"])
        # Def 8 floor: drop size-<2 remnant edges the peeler leaves behind
        # (p-qualified single-survivor edges are not valid partial hyperedges).
        GC = func.prune_small_edges(GC, min_size=2)

    elif algorithm == "ks_core":
        C, rt = ks_core.run(G, params["k"], params["s"], params.get("c", 1.0))
        GC = func.get_partially_subgraph(G_orig, C)

    elif algorithm == "kab_truss":
        if params["a"] == 0:                      # a==0 (ce_truss) path: C++ port not implemented
            C, HE, rt = kab_truss.run(G, params["k"], params["a"], params["b"], network_path)
        else:
            C, HE, rt = cpp_truss_wrap.run_kab(network_path, params["k"], params["a"], params["b"])
        GC = func.get_node_edge_supgraph(G_orig, C, HE)

    elif algorithm == "kinout_truss":
        C, HE, rt = cpp_truss_wrap.run_kinout(network_path, params["k_in"], params["k_out"])
        GC = func.get_node_edge_supgraph(G_orig, C, HE)

    elif algorithm == "hyper_k_truss":
        HE, rt = hyper_k_truss.run(network_path, dataset_name, params["k"])
        GC = func.get_edge_subgraph(G_orig, HE)

    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")

    return GC, rt


# ─────────────────────────────────────────────
# Child-process worker
# ─────────────────────────────────────────────

def _job_worker(queue, model, G_orig, params, network_path, dataset_name, compute_quality):
    try:
        with PeakMemoryTracker() as mem:
            t0 = time.time()
            GC, model_rt = _run_model(model, G_orig, params, network_path, dataset_name)
            t_total = time.time() - t0
        # Peak memory = max of three complementary sources:
        #   self_hwm  : kernel high-water mark of THIS worker (ru_maxrss).
        #               Deterministic — never misses, unlike the 50ms sampler,
        #               which gets GIL-starved during CPU-bound pure-python
        #               phases and then records only the post-load baseline.
        #   child_hwm : kernel peak of any reaped child (the C++ truss exes),
        #               which the parent's own counters do not include.
        #   mem.peak_mb: the sampling tracker — kept as the safety net for the
        #               CONCURRENT parent+child sum (e.g. kab_truss ~7.6GB =
        #               2GB python + 5.7GB C++ at the same moment); it works
        #               reliably there because the parent blocks on the child
        #               (GIL released).
        self_hwm  = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        child_hwm = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0
        measures = (ev.compute_all_measures(G_orig, GC) if compute_quality
                    else ev.compute_size_measures(G_orig, GC))
        queue.put({
            "ok": True,
            "runtime_sec_total":          round(t_total, 4),
            "runtime_sec_model_reported": round(model_rt, 4) if model_rt else "",
            "peak_memory_mb":             round(max(mem.peak_mb, self_hwm, child_hwm), 2),
            "measures":                   measures,
        })
    except Exception:
        queue.put({"ok": False, "error": traceback.format_exc(limit=3).replace("\n", " | ")})


# ─────────────────────────────────────────────
# Single experiment
# ─────────────────────────────────────────────

def _parse_n(dataset_name):
    """Target node count from names like "..._n10k_r0" / "..._n40000_r0".
    Supports k/m suffixes (10k -> 10000, 2m -> 2000000); 0 if no match (e.g. real)."""
    import re as _re
    try:
        m_size = _re.search(r"_n(\d+)([kKmM]?)", dataset_name)
        if m_size:
            base = int(m_size.group(1))
            suf  = m_size.group(2).lower()
            return base * (1000 if suf == "k" else 1_000_000 if suf == "m" else 1)
    except Exception:
        pass
    return 0


def _base_row(eq, dataset_type, dataset_name, generator_setting,
              model, param_group, params, repeat_id, seed):
    """Build a CSV row dict with metadata filled and all measures blank."""
    row = {c: "" for c in ALL_COLS}
    row.update({
        "experiment_id":     f"{eq}_{dataset_name}_{model}_{param_group}_r{repeat_id}",
        "eq":                eq,
        "dataset_type":      dataset_type,
        "dataset_name":      dataset_name,
        "generator_setting": generator_setting,
        "n":                 _parse_n(dataset_name),
        "model":             model,
        "parameter_group":   param_group,
        "parameters":        str(params),
        "repeat_id":         repeat_id,
        "seed":              seed,
    })
    return row


def run_one(eq, dataset_type, dataset_name, generator_setting,
            model, param_group, params, repeat_id, seed,
            compute_quality, writer, csvfile):
    exp_id = f"{eq}_{dataset_name}_{model}_{param_group}_r{repeat_id}"

    network_path = str(DATASETS_DIR / dataset_name / "network.hyp")
    if not os.path.exists(network_path):
        print(f"  [SKIP] {network_path} not found")
        return

    row = _base_row(eq, dataset_type, dataset_name, generator_setting,
                    model, param_group, params, repeat_id, seed)

    G_orig = func.Hypergraph()
    G_orig.load_from_file(network_path)

    queue = mp.Queue()
    proc  = mp.Process(target=_job_worker,
                       args=(queue, model, G_orig, params, network_path,
                             dataset_name.split("/")[-1], compute_quality))
    proc.start()

    kill_reason = [None]

    def _mem_watch():
        try:
            child_proc = psutil.Process(proc.pid)
            while proc.is_alive():
                # kill if *system-wide* memory usage exceeds the threshold.
                # .percent is based on `available`, so reclaimable buff/cache
                # is NOT counted as pressure.
                if psutil.virtual_memory().percent > JOB_MEMORY_LIMIT_PERCENT:
                    kill_reason[0] = "oom"
                    for gc in child_proc.children(recursive=True):
                        try:
                            gc.kill()
                        except Exception:
                            pass
                    proc.kill()
                    return
                time.sleep(0.5)
        except Exception:
            pass

    watcher = threading.Thread(target=_mem_watch, daemon=True)
    watcher.start()

    proc.join(JOB_TIMEOUT_SEC)

    if proc.is_alive():
        kill_reason[0] = "timeout"
        proc.kill()
        proc.join()

    if kill_reason[0]:
        row["status"]            = kill_reason[0]
        row["runtime_sec_total"] = "inf"
        print(f"  [{kill_reason[0].upper()}] {exp_id}")
    elif not queue.empty():
        res = queue.get_nowait()
        if res["ok"]:
            row["runtime_sec_total"]          = res["runtime_sec_total"]
            row["runtime_sec_model_reported"] = res["runtime_sec_model_reported"]
            row["peak_memory_mb"]             = res["peak_memory_mb"]
            row.update(res["measures"])
            row["status"] = "ok"
            print(f"  [OK] {exp_id}  nodes={row['num_nodes']}  rt={res['runtime_sec_total']}s")
        else:
            row["status"]        = "failed"
            row["error_message"] = res["error"]
            print(f"  [FAIL] {exp_id}: {res['error'][:80]}")
    else:
        row["status"]        = "failed"
        row["error_message"] = f"worker died (exitcode={proc.exitcode})"
        print(f"  [FAIL] {exp_id}: worker died (exitcode={proc.exitcode})")

    writer.writerow([row[c] for c in ALL_COLS])
    csvfile.flush()


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def _derive_metadata(eq, dataset):
    """Derive dataset_type, generator_setting, repeat_id, seed, compute_quality."""
    dataset_type    = "real" if dataset.startswith("real/") else "synthetic"
    compute_quality = eq in ("EQ3", "EQ4")

    if dataset_type == "real":
        return dataset_type, "real", 0, 0, compute_quality

    basename = dataset.split("/")[-1]
    r_idx    = basename.rfind("_r")
    if r_idx >= 0:
        repeat_id         = int(basename[r_idx + 2:])
        generator_setting = basename[:r_idx]
    else:
        repeat_id         = 0
        generator_setting = basename
    seed = 42 + repeat_id
    return dataset_type, generator_setting, repeat_id, seed, compute_quality


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eq",          required=True, choices=["EQ1", "EQ2", "EQ3", "EQ4"])
    parser.add_argument("--dataset",     required=True,
                        help="Dataset path, e.g. real/walmart or synthetic/gamma_2_r0")
    parser.add_argument("--model",       required=True,
                        help="Model name, e.g. kinout_truss_kout3")
    parser.add_argument("--param-group", required=True, choices=["q1", "q2", "q3", "q4"],
                        help="Parameter group to run")
    parser.add_argument("--output",      required=True,
                        help="Output CSV path (row is appended if file exists)")
    parser.add_argument("--mark-skipped", action="store_true",
                        help="Do not run; just append a row with status=skipped "
                             "(used by the runner to record early-skipped jobs).")
    args = parser.parse_args()

    if args.model not in MODEL_PARAMS:
        print(f"Unknown model: {args.model}")
        print(f"Available: {list(MODEL_PARAMS)}")
        raise SystemExit(1)

    dataset_type, generator_setting, repeat_id, seed, compute_quality = \
        _derive_metadata(args.eq, args.dataset)
    params = MODEL_PARAMS[args.model][args.param_group]

    out_path = args.output
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    file_exists = os.path.exists(out_path)
    with open(out_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(ALL_COLS)

        if args.mark_skipped:
            row = _base_row(args.eq, dataset_type, args.dataset, generator_setting,
                            args.model, args.param_group, params, repeat_id, seed)
            row["status"] = "skipped"
            writer.writerow([row[c] for c in ALL_COLS])
            f.flush()
            print(f"  [SKIPPED] {row['experiment_id']}")
            print(f"Done → {out_path}")
            return

        print(f"=== {args.eq} | {args.dataset} | {args.model} | {args.param_group} "
              f"| params={params} ===")
        run_one(
            eq=args.eq,
            dataset_type=dataset_type,
            dataset_name=args.dataset,
            generator_setting=generator_setting,
            model=args.model,
            param_group=args.param_group,
            params=params,
            repeat_id=repeat_id,
            seed=seed,
            compute_quality=compute_quality,
            writer=writer,
            csvfile=f,
        )

    print(f"Done → {out_path}")


if __name__ == "__main__":
    main()
