# Cohesive Subgraph Models on Hypergraphs

> **Note for authors:** fill in the paper title, authors, venue, and BibTeX below before publishing.

Reference implementations and a unified benchmark harness for a family of
**cohesive subgraph models on hypergraphs**, covering both the *core* family
(degeneracy-style peeling) and the *truss* family (edge/pair-support peeling).

If you use this code, please cite:

```bibtex
@inproceedings{TODO,
  title     = {TODO},
  author    = {TODO},
  booktitle = {TODO},
  year      = {2026}
}
```

---

## Models

All models expose a single `run(...)` entry point and return the surviving
subgraph together with the model-reported runtime. `k` is the common degree
threshold; the second parameter (when present) is the model-specific one.

### Core family

| Model | Module | `run(...)` signature | Extra parameter |
|-------|--------|----------------------|-----------------|
| k-hypercore            | `k_hypercore.py`   | `run(network_path, dataset_name, k)`        | — |
| neighbour k-core       | `nbr_kd_core.py`   | `run(network_path, dataset_name, k, d, "nbr_k_core")` | — |
| (k,d)-core             | `nbr_kd_core.py`   | `run(network_path, dataset_name, k, d, "kd_core")`    | `d` (neighbour degree) |
| (k,t)-hypercore        | `kt_hypercore.py`  | `run(G, k, t, dataset_name)`                | `t` (fraction) |
| (k,q)-core             | `kq_core.py`       | `run(G, k, q)`                              | `q` (cardinality) |
| (k,g)-core             | `kg_core.py`       | `run(G, k, g)`                              | `g` (co-occurrence) |
| (k,g,p)-core           | `kgp_core.py`      | `run(G, k, g, p)`                           | `g`, `p` (p-rule filter) |
| (k,s)-core             | `ks_core.py`       | `run(G, k, s, c=1.0)`                       | `s` (DCS threshold), `c` |

### Truss family

| Model | Module | `run(...)` signature | Extra parameter |
|-------|--------|----------------------|-----------------|
| hyper k-truss   | `hyper_k_truss.py` | `run(network_path, dataset_name, k)`       | — |
| (k,a,b)-truss   | `kab_truss.py`     | `run(G, k, a, b, network_path)`            | `a`, `b` (pair support) |
| (k_in,k_out)-truss | `kinout_truss.py` | `run(G, k_in, k_out, ...)`             | `k_in`, `k_out` |

### Decomposition / coreness variants

Standalone scripts that compute the full decomposition (per-node coreness /
per-edge trussness) rather than a single-threshold subgraph:
`k_hypercoreness.py`, `kt_hypercoreness.py`, `kg_coreness.py`,
`hyper_k_trussness.py`, and the C++ `nbr_kd_coreness/`.

---

## Repository layout

```
.
├── run_experiment.py        # unified single-job runner (one dataset × model × params → one CSV row)
├── func.py                  # Hypergraph class, I/O, projection/clique graphs, caching
├── eval.py                  # subgraph quality measures (density, retention, conductance, …)
├── etc.py                   # small helpers
├── cpp_truss_wrap.py        # thin wrapper around the C++ (k,a,b)- / (k_in,k_out)-truss ports
├── <model>.py               # the model implementations listed above
├── cpp_truss/               # C++ ports of kab_truss / kinout_truss (+ verification scripts)
├── HTRUSS/                  # C++ hyper k-truss decomposition (third-party, see Attribution)
└── nbr_kd_coreness/         # C++ (k,d)-core / neighbour-coreness decomposition
```

---

## Installation

Python 3.10+ is recommended.

```bash
conda create -n hcs python=3.10 -y
conda activate hcs
pip install -r requirements.txt
```

### Building the C++ components

The C++ implementations are optional accelerators; the pure-Python models run
without them. Build them if you want the fast paths used in the paper.

```bash
# (k,a,b)-truss and (k_in,k_out)-truss ports
g++ -std=c++17 -O2 -o cpp_truss/kab_truss    cpp_truss/kab_truss.cpp
g++ -std=c++17 -O2 -o cpp_truss/kinout_truss cpp_truss/kinout_truss.cpp

# hyper k-truss (HTRUSS) — parallel PushForward variant used in the paper
cd HTRUSS && cmake -B build && cmake --build build     # see HTRUSS/README.md

# (k,d)-core / neighbour-coreness
cd nbr_kd_coreness && bash runkd.sh
```

---

## Data

Datasets are **not** bundled in this repository (size and licensing). The
real-world hypergraphs used in the paper follow the standard node-labeled
hypergraph format and can be obtained from the public benchmark collections,
e.g. Austin R. Benson's hypergraph datasets: <https://www.cs.cornell.edu/~arb/data/>.
Synthetic instances are generated with the H-ABCD generator
(<https://github.com/bkamins/ABCDHypergraphGenerator.jl>).

**Expected on-disk layout** (relative to the runner):

```
../datasets/real/<name>/network.hyp     # one hyperedge per line
../datasets/synthetic/<...>/network.hyp
../cache/                               # auto-generated preprocessing caches
../results/                             # experiment output CSVs
```

`network.hyp` format — one hyperedge per line, node ids separated by commas or
whitespace (hyperedges with fewer than 2 nodes are dropped on load):

```
2,5
1,2
11,8
```

---

## Usage

### Single experiment via the runner

```bash
python run_experiment.py \
    --eq EQ3 \
    --dataset real/walmart \
    --model kq_core \
    --param-group q2 \
    --output ../results/walmart_kq_core_q2.csv
```

Available models and their `q1…q4` parameter groups are defined in
`MODEL_PARAMS` at the top of `run_experiment.py`.

### Calling a model directly

```python
import func, eval as ev
import kq_core

H = func.Hypergraph()
H.load_from_file("../datasets/real/walmart/network.hyp")
G = func.hypergraph_to_networkx(H)

C, runtime = kq_core.run(G, k=4, q=4)   # surviving subgraph + runtime
print(ev.compute_all_measures(G, C))
```

---

## Attribution

`HTRUSS/` contains a hyper k-truss decomposition framework derived from
[VeryLargeGraph/HTRUSS](https://github.com/VeryLargeGraph/HTRUSS)
(MIT License, © 2021 Haizs Chen); its `LICENSE` is retained in that directory.
`HTRUSS/src/TrussDecomp_pf.cpp` (parallel PushForward + heap peeling) is our
addition on top of that framework.

## License

This project is released under the MIT License (see `LICENSE`). The bundled
`HTRUSS/` directory keeps its own upstream MIT license.
