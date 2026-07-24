import subprocess
import os
from pathlib import Path
import csv

import func

# === Base paths ===
BASE_DIR = Path("./nbr_kd_coreness").resolve()   # project root
SRC_DIR = BASE_DIR                               # contains main.cpp, etc.
OUT_DIR = Path("../output").resolve()                             # txt files saved here

def build(binary_path):
    bin_path = SRC_DIR / binary_path
    # Skip recompile if the binary is already up-to-date w.r.t. all sources/headers.
    # This avoids a "Text file busy" (Errno 26) race when parallel jobs would each
    # rebuild the same binary while another process is executing it.
    srcs = [SRC_DIR / s for s in ("main.cpp", "algorithms.cpp", "hypergraph.cpp")]
    if bin_path.exists():
        newest_src = max(p.stat().st_mtime for p in srcs + list(SRC_DIR.glob("*.h")))
        if bin_path.stat().st_mtime >= newest_src:
            return str(bin_path)
    cmd = ["g++", "-std=c++17", "-O3", "main.cpp", "algorithms.cpp", "hypergraph.cpp", "-I.", "-o", str(bin_path)]
    subprocess.run(cmd, check=True, cwd=SRC_DIR)
    return str(bin_path)


def run(network, algorithm, binary_path="main"):
    exe = build(binary_path)

    # Run the C++ executable
    subprocess.run(
        [str(Path(exe).resolve()), network, algorithm],
        cwd=SRC_DIR,
        check=True
    )

    dataset_name = network.split('/')[-1]   # last component = dataset id (matches C++ getDatasetName)
    return func.read_decomp_from_csv(dataset_name, algorithm)
