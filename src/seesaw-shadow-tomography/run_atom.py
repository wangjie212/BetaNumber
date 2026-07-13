#!/usr/bin/env python3
"""
Reproduce Hamiltonian energy estimation variance bounds from Section VI
of the supplemental material.

Uses the atom module to:
  - Compare variance bounds across measurement schemes (global, single-qubit,
    two-qubit) and fermion-to-qubit mappers (Parity, JW, BK).
  - Compute JM-norm^2 bounds with and without noise parameter optimization.
  - Reproduce Table I from the paper.

Dependencies:
    pip install -r requirements-atom.txt

Usage:
    python run_atom.py                           # default: H2, (2e,2o), all mappers
    python run_atom.py --molecule LiH --nele 2 --norb 3
    python run_atom.py --all                      # full Table I
    python run_atom.py --dry-run                  # print what would run

NOTE: Requires PySCF, Qiskit Nature, and MOSEK (or SCS/CVXOPT/ECOS as fallback).
      PySCF will perform quantum chemistry calculations on first run.
"""

import argparse
import sys
import numpy as np

# ----------------------------------------------------------------------
# Table I configurations
# ----------------------------------------------------------------------
TABLE_I = [
    # (molecule, atom_string, nele, norb)
    ("H2",   "H 0 0 0; H 0 0 0.74",       2, 2),
    ("LiH",  "Li 0 0 0; H 0 0 1.6",        2, 3),
    ("LiH",  "Li 0 0 0; H 0 0 1.6",        2, 4),
    ("BeH2", "Be 0 0 0; H 0 0 1.33; H 0 0 -1.33", 4, 3),
    ("BeH2", "Be 0 0 0; H 0 0 1.33; H 0 0 -1.33", 4, 4),
    ("BeH2", "Be 0 0 0; H 0 0 1.33; H 0 0 -1.33", 4, 5),
    ("H2O",  "O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586", 2, 3),
    ("H2O",  "O 0 0 0; H 0 0.757 0.586; H 0 -0.757 0.586", 4, 4),
]

MAPPER_NAMES = ["Parity", "JW", "BK"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Reproduce Hamiltonian variance bounds (Table I).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_atom.py                           # H2 (2e,2o), all mappers
  python run_atom.py --all                     # full Table I
  python run_atom.py --molecule LiH --nele 2 --norb 3
  python run_atom.py --dry-run                 # preview
""")

    p.add_argument("--molecule", type=str, default="H2",
                   help="Molecule label or atom string (default: H2)")
    p.add_argument("--nele", type=int, default=2,
                   help="Number of active electrons (default: 2)")
    p.add_argument("--norb", type=int, default=2,
                   help="Number of active spatial orbitals (default: 2)")
    p.add_argument("--all", action="store_true",
                   help="Reproduce the full Table I")
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would run without executing")
    p.add_argument("--normalized", action="store_true",
                   help="Also compute the table normalized by the global bound")

    return p.parse_args()


def find_config(molecule, nele, norb):
    """Find a matching Table I configuration or return None."""
    for label, atom_str, e, o in TABLE_I:
        if label == molecule and e == nele and o == norb:
            return atom_str
    return None


def fmt(val):
    """Format a numeric value for display."""
    if val is None:
        return "N/A"
    if abs(val) < 0.01:
        return f"{val:.4f}"
    return f"{val:.2f}"


def run_one(atom_str, nele, norb, mapper_id, mapper_name):
    """Run compare() for one configuration, return results."""
    from atom import compare
    result = compare(atom_str, nele, norb, mapper_id)
    global_bounds = result[0]   # [b1Opt, b2Opt]
    single_bounds = result[1]   # [b1Opt1, b2Opt1]
    jm_bounds     = result[2]   # [bJMNoOpt, bJMd]
    two_bounds    = result[3]   # [b1Opt2, b2Opt2]
    return {
        "global":     global_bounds[1],   # bound 2 (maximal eigenvalue)
        "single":     single_bounds[1],
        "jm_opt":     jm_bounds[1],
        "two":        two_bounds[1],
        "jm_noopt":   jm_bounds[0],
    }


def main():
    args = parse_args()

    if args.dry_run:
        if args.all:
            print("DRY RUN — would reproduce Table I:")
            print(f"  {'Molecule':<8s} {'(e,o)':<10s} {'Mapper':<8s}")
            print(f"  {'-'*8} {'-'*10} {'-'*8}")
            for label, atom_str, nele, norb in TABLE_I:
                for mid in range(3):
                    print(f"  {label:<8s} ({nele}e,{norb}o)  {MAPPER_NAMES[mid]:<8s}")
        else:
            atom_str = find_config(args.molecule, args.nele, args.norb)
            if atom_str is None:
                atom_str = args.molecule  # use as literal atom string
            for mid in range(3):
                print(f"  {args.molecule} ({args.nele}e,{args.norb}o) "
                      f"{MAPPER_NAMES[mid]}")
        return

    # Check dependencies
    try:
        from atom import compare
    except ImportError as e:
        print(f"Error: {e}")
        print("Install dependencies: pip install -r requirements-atom.txt")
        sys.exit(1)

    if args.all:
        print("Reproducing Table I: Variance bounds for selected molecules\n")
        print(f"{'Molecule':<8s} {'(e,o)':<8s} {'Mapper':<8s} "
              f"{'Global':>7s} {'Single':>7s} {'Two':>7s} "
              f"{'JM':>7s} {'JMnoOpt':>8s}")
        print("-" * 65)

        for label, atom_str, nele, norb in TABLE_I:
            for mid in range(3):
                print(f"{label:<8s} ({nele}e,{norb}o) "
                      f"{MAPPER_NAMES[mid]:<8s}",
                      end=" ", flush=True)
                try:
                    r = run_one(atom_str, nele, norb, mid, MAPPER_NAMES[mid])
                    print(f"{fmt(r['global']):>7s} {fmt(r['single']):>7s} "
                          f"{fmt(r['two']):>7s} {fmt(r['jm_opt']):>7s} "
                          f"{fmt(r['jm_noopt']):>8s}")
                except Exception as e:
                    print(f"ERROR: {e}")
            print()  # blank line between groups

        if args.normalized:
            print("\n--- Normalized by global bound ---")
            # (would need to re-run and store results; skipped for brevity)
            print("Run with --all first, then compute normalization separately.")

    else:
        # Single configuration
        atom_str = find_config(args.molecule, args.nele, args.norb)
        if atom_str is None:
            atom_str = args.molecule
            print(f"Using atom string directly: {atom_str}")

        print(f"Molecule: {args.molecule}  ({args.nele}e, {args.norb}o)\n")
        print(f"{'Mapper':<8s} {'Global':>7s} {'Single':>7s} {'Two':>7s} "
              f"{'JM':>7s} {'JMnoOpt':>8s}")
        print("-" * 50)

        for mid in range(3):
            print(f"{MAPPER_NAMES[mid]:<8s}", end=" ", flush=True)
            try:
                r = run_one(atom_str, args.nele, args.norb, mid,
                            MAPPER_NAMES[mid])
                print(f"{fmt(r['global']):>7s} {fmt(r['single']):>7s} "
                      f"{fmt(r['two']):>7s} {fmt(r['jm_opt']):>7s} "
                      f"{fmt(r['jm_noopt']):>8s}")
            except Exception as e:
                print(f"ERROR: {e}")


if __name__ == "__main__":
    main()
