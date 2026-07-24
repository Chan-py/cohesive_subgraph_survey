# cpp_truss — C++ ports of kab_truss / kinout_truss

Fast, low-memory C++ reimplementations of the Python models in
`code/kab_truss.py` and `code/kinout_truss.py`, mirroring the pattern already
used for `hyper_k_truss` (HTRUSS `TrussDecomp.exe`). Each binary reads a
`network.hyp` directly and writes the surviving vertex set + hyperedge-id set,
**identical** to the corresponding Python `run(...)` output.

Kept fully isolated from the running experiment pipeline: nothing here imports or
edits `run_experiment.py` / `kab_truss.py` / `kinout_truss.py`.

## Build
```
g++ -std=c++17 -O2 -o kinout_truss kinout_truss.cpp
g++ -std=c++17 -O2 -o kab_truss    kab_truss.cpp
```
(Compilation peak RSS ~150 MB each — negligible.)

## Usage
```
./kinout_truss <network.hyp> <k_in> <k_out> <out.txt>
./kab_truss    <network.hyp> <k> <a> <b>   <out.txt>      # a>=1 path only
```
`out.txt`: line 1 = surviving original vertex ids (sorted), line 2 = surviving
hyperedge ids (1-based, file order, sorted).

## Correctness
Both models peel to a **unique fixpoint**, so the result set is independent of
peeling order — the C++ result must equal the Python result exactly.

- `verify_kinout.py` — 8 random hypergraphs + real `contact`, 12 `(k_in,k_out)`
  combos each → **108/108 match** (nodes & hids sets).
- `verify_kab.py`   — 8 random hypergraphs + real `contact`, 12 `(k,a,b)` combos
  each → **108/108 match**.

Run from the `code/` dir:
```
python cpp_truss/verify_kinout.py
python cpp_truss/verify_kab.py
```

## Key equivalences used (derived from the Python code, verified)
- **kinout**: `sup_in[(e,u,v)] = |e|-2` (in-peel ⇔ hyperedge size `< k_in+2`);
  `sup_node_out[(e,u)] = max_{v∈e\{u}} S(u,v) - (|e|-2)` with
  `S(u,v)=Σ_{w∈N(u)∩N(v)} mult(u,w)·mult(v,w)`.
- **kab (a>0)**: projection graph on vertex-pairs; a triangle {u,v,w} is valid iff
  `sigma>=a` and `omega[u]+omega[v]+omega[w]-2·sigma>=b`
  (i.e. `#he⊇{x,y,z} >= a` and `#he containing ≥2 of x,y,z >= b`);
  node-truss peel of pairs with `< k` valid triangles, after a `2k`-core prefilter.

## Memory
The C++ kab port does **not** materialise the per-pair hyperedge-id sets that the
networkx version keeps (the dominant OOM sink); the final hid union is recovered
by scanning hyperedges against the surviving pair set. Only `omega` (ints),
per-edge `sigma`, and projection adjacency are held.
