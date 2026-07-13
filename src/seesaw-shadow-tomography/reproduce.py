#!/usr/bin/env python3
"""
Reproduce benchmark results from Section IV.B of the supplemental material.

For each graph in the selected benchmark dataset, runs the three-stage pipeline:
  1. seesaw_oracle  — see-saw for fixed weight vector
  2. mwu_ensemble   — multiplicative weight update outer loop
  3. solve_lp        — linear program for best convex combination

Usage:
    python reproduce.py                       # 8-vertex graphs (default)
    python reproduce.py --dataset 9           # 9-vertex graphs only
    python reproduce.py --dataset all         # both 8- and 9-vertex
    python reproduce.py --eta auto            # compute eta from formula
    python reproduce.py --trials 3            # best of 3 runs per graph
    python reproduce.py --dry-run             # print plan without executing

NOTE: 9-vertex graphs use 9 Pauli strings → 512x512 density matrices and
      T=10000 iterations. Expect significantly longer runtime. The default
      (--dataset 8) processes only the faster 8-vertex graphs.
"""

import argparse
import ast
import os
import sys
import time
import numpy as np
from math import sqrt, log

# Import from the local module
from seesaw import seesaw_oracle, mwu_ensemble, compute_lp_coefficients, solve_lp

# ----------------------------------------------------------------------
# Defaults (matching the paper: T=10000, tol≈5e-6, max_iter≈150)
# ----------------------------------------------------------------------
DEFAULTS = {
    "T": 10000,
    "eta": 0.02,
    "seesaw_max_iter": 150,
    "seesaw_tol": 5e-6,
    "trials": 1,
    "dataset": "8",
    "data_dir": "data",
    "output": None,          # set from dataset + data_dir in main()
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Reproduce seesaw shadow tomography benchmark results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reproduce.py                       # 8-vertex graphs (default)
  python reproduce.py --dataset 9           # 9-vertex graphs only
  python reproduce.py --dataset all         # both 8- and 9-vertex
  python reproduce.py --dry-run             # preview without running
  python reproduce.py --eta auto --trials 3 # auto eta, 3 runs per graph
""")

    p.add_argument("--T", type=int, default=DEFAULTS["T"],
                   help=f"MWU iterations (default: {DEFAULTS['T']})")
    p.add_argument("--eta", type=str, default=str(DEFAULTS["eta"]),
                   help="MWU learning rate; 'auto' for eta = sqrt(8 ln m / T) "
                        f"(default: {DEFAULTS['eta']})")
    p.add_argument("--seesaw-max-iter", type=int,
                   default=DEFAULTS["seesaw_max_iter"],
                   help="Oracle inner iterations "
                        f"(default: {DEFAULTS['seesaw_max_iter']})")
    p.add_argument("--seesaw-tol", type=float,
                   default=DEFAULTS["seesaw_tol"],
                   help="Oracle convergence tolerance "
                        f"(default: {DEFAULTS['seesaw_tol']})")
    p.add_argument("--trials", type=int, default=DEFAULTS["trials"],
                   help="Repetitions per graph with different seeds; "
                        f"best delta is kept (default: {DEFAULTS['trials']})")
    p.add_argument("--dataset", type=str, default=DEFAULTS["dataset"],
                   choices=["8", "9", "all"],
                   help="Which dataset to process: 8 (default, 9 graphs, "
                        "fast), 9 (295 graphs, slow), all (304 graphs)")
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would run without executing")
    p.add_argument("--output", type=str, default=DEFAULTS["output"],
                   help="Output file path (default: auto from --dataset)")
    p.add_argument("--data-dir", type=str, default=DEFAULTS["data_dir"],
                   help=f"Directory for input files "
                        f"(default: {DEFAULTS['data_dir']})")

    args = p.parse_args()

    # Resolve data_dir relative to script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    args.data_dir_abs = os.path.join(script_dir, args.data_dir)

    # Auto output naming from dataset
    if args.output is None:
        if args.dataset == "8":
            args.output = os.path.join(args.data_dir_abs, "newres8.txt")
        elif args.dataset == "9":
            args.output = os.path.join(args.data_dir_abs, "newres9.txt")
        else:  # all
            args.output = os.path.join(args.data_dir_abs, "newres_all.txt")

    # Parse eta
    if args.eta.lower() == "auto":
        args.eta = None  # computed per-graph
    else:
        args.eta = float(args.eta)

    return args


def load_test(data_dir, dataset):
    """Load test file → list of (graph6_code, paulis)."""
    path = os.path.join(data_dir, f"test{dataset}.txt")
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = ast.literal_eval(line)
            entries.append((entry[0], entry[2]))  # (graph6, paulis)
    return entries


def load_ndelta(data_dir, dataset):
    """Load ndelta file → list of (index, old_value)."""
    path = os.path.join(data_dir, f"ndelta{dataset}_3_hard.txt")
    with open(path) as f:
        content = f.read()
    data = ast.literal_eval(content)
    return [(item[0], item[2]) for item in data]  # (index, old_value)


def compute_eta(m, T):
    """Learning rate from the paper: eta = sqrt(8 ln m / T)."""
    return sqrt(8 * log(m) / T)


def process_one_graph(idx, old_val, paulis, graph6, args):
    """Run the full pipeline for one graph, return (graph6, old_val, best_delta)."""
    m = len(paulis)
    eta = args.eta if args.eta is not None else compute_eta(m, args.T)
    best_delta = -np.inf

    for run in range(args.trials):
        seed = idx * 1000 + run
        np.random.seed(seed)

        try:
            ensemble = mwu_ensemble(paulis, args.T, eta,
                                   seesaw_max_iter=args.seesaw_max_iter,
                                   seesaw_tol=args.seesaw_tol)
            rhos = [rho for rho, _ in ensemble]
            sigmas = [sigma for _, sigma in ensemble]
            e_rho, e_sigma, _ = compute_lp_coefficients(rhos, sigmas, paulis)
            delta_opt, _, _ = solve_lp(e_rho, e_sigma)
            if delta_opt > best_delta:
                best_delta = delta_opt
        except Exception as e:
            print(f"  {graph6} run {run+1} error: {e}", file=sys.stderr)
            continue

    if best_delta == -np.inf:
        best_delta = float('nan')
    return (graph6, old_val, best_delta)


def process_dataset(dataset, args):
    """Process one dataset (8 or 9), return list of result tuples."""
    test8_entries = load_test(args.data_dir_abs, dataset)
    ndelta_entries = load_ndelta(args.data_dir_abs, dataset)

    print(f"Dataset {dataset}: {len(test8_entries)} graphs, "
          f"{len(ndelta_entries)} benchmark entries")
    print(f"Parameters: T={args.T}, "
          f"eta={'auto' if args.eta is None else args.eta}, "
          f"seesaw_max_iter={args.seesaw_max_iter}, "
          f"seesaw_tol={args.seesaw_tol}, trials={args.trials}")
    print()

    if args.dry_run:
        print(f"DRY RUN — would process these graphs:")
        for idx, old_val in ndelta_entries:
            if idx < 1 or idx > len(test8_entries):
                print(f"  [{idx}] (out of range, skipping)")
                continue
            graph6, paulis = test8_entries[idx - 1]
            m = len(paulis)
            eta = args.eta if args.eta is not None else compute_eta(m, args.T)
            n_qubits = len(paulis[0]) if paulis else "?"
            print(f"  [{idx}] {graph6}  {n_qubits} qubits  m={m}  "
                  f"eta={eta:.6f}")
        print()
        return []

    results = []
    t_start = time.time()
    n_total = len(ndelta_entries)

    for i, (idx, old_val) in enumerate(ndelta_entries):
        if idx < 1 or idx > len(test8_entries):
            print(f"  [{idx}] (out of range, skipping)")
            continue

        graph6, paulis = test8_entries[idx - 1]
        n_qubits = len(paulis[0])
        dim = 2 ** n_qubits

        print(f"[{i+1}/{n_total}] {graph6} "
              f"({n_qubits} qubits, dim={dim}x{dim})...",
              end=" ", flush=True)

        graph6_out, old_val_out, new_val = process_one_graph(
            idx, old_val, paulis, graph6, args)
        results.append((graph6_out, old_val_out, new_val))
        elapsed = time.time() - t_start
        print(f"delta = {new_val:.8f}  [{elapsed:.0f}s]")

    return results


def main():
    args = parse_args()

    # Timing warning for 9-vertex and all
    if args.dataset in ("9", "all") and not args.dry_run:
        print("⚠  WARNING: Processing 9-vertex graphs (512x512 density matrices, "
              "T=10000 iterations) may take several hours.\n"
              "   Use --dataset 8 for the faster 8-vertex benchmark (default).\n"
              "   Use --dry-run to preview without running.\n"
              "   Press Ctrl+C to abort, or wait 5 seconds to continue...\n",
              flush=True)
        try:
            time.sleep(5)
        except KeyboardInterrupt:
            print("\nAborted.")
            sys.exit(0)

    all_results = []

    if args.dataset in ("8", "all"):
        all_results.extend(process_dataset("8", args))

    if args.dataset in ("9", "all"):
        all_results.extend(process_dataset("9", args))

    if args.dry_run:
        if args.dataset == "8":
            print(f"Would write results to {args.output}")
        elif args.dataset == "9":
            print(f"Would write results to {args.output}")
        else:
            print(f"Would write results to {args.output}")
        return

    # Write output
    if all_results:
        with open(args.output, "w") as f:
            for graph6, old_val, new_val in all_results:
                f.write(f"[{graph6}, {old_val}, {new_val}]\n")
        print(f"\nDone. {len(all_results)} results written to {args.output}")
    else:
        print("\nNo results to write.")


if __name__ == "__main__":
    main()
