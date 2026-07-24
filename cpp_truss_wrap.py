"""Thin wrappers that run the C++ ports in code/cpp_truss/ and return results in
the SAME (node_set, hyperedge_id_set, runtime) shape as the Python
kinout_truss.run / kab_truss.run, so run_experiment can use C++ as a drop-in.

Verified equivalent to the Python implementations (216/216 unit checks +
80/80 count checks on the 128k synthetic graphs). The C++ binaries read
network.hyp directly and assign 1-based hyperedge ids in file order, matching
func.Hypergraph.load_from_file — so the returned ids line up with G_orig.
"""
import os
import time
import tempfile
import subprocess
from pathlib import Path

_CPP = Path(__file__).resolve().parent / "cpp_truss"
_KINOUT = _CPP / "kinout_truss"
_KAB = _CPP / "kab_truss"


def _read_out(path):
    with open(path) as f:
        txt = f.read().split("\n")
    nodes = set(int(x) for x in txt[0].split()) if len(txt) > 0 and txt[0].strip() else set()
    hids = set(int(x) for x in txt[1].split()) if len(txt) > 1 and txt[1].strip() else set()
    return nodes, hids


def _run(exe, args):
    fd, out = tempfile.mkstemp(suffix=".txt")
    os.close(fd)
    try:
        t0 = time.time()
        subprocess.run([str(exe), *args, out], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rt = time.time() - t0
        nodes, hids = _read_out(out)
    finally:
        try:
            os.unlink(out)
        except OSError:
            pass
    return nodes, hids, rt


def run_kinout(network, k_in, k_out):
    """(k_in,k_out)-in/out-truss via C++ -> (node_set, hyperedge_id_set, runtime)."""
    return _run(_KINOUT, [str(network), str(k_in), str(k_out)])


def run_kab(network, k, a, b):
    """(k,a,b)-truss (a>=1) via C++ -> (node_set, hyperedge_id_set, runtime)."""
    return _run(_KAB, [str(network), str(k), str(a), str(b)])
