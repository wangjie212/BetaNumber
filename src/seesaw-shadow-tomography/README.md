# Seesaw Shadow Tomography

Reproducible implementation of the see-saw + multiplicative weight update
algorithm for constructing optimal measurement strategies in Pauli shadow
tomography, from Section IV.B of the supplemental material for
*"Optimal strategies for shadow tomography with limited resources"*
(Xu, Ye, Wang, Koßmann, Winter, Schwonnek, 2026).

## What this does

Given a set of Pauli observables {S₁, …, S_m}, the algorithm computes a
lower bound on the sample complexity parameter **δ** and constructs an
explicit measurement ensemble that achieves it — via a three-stage pipeline:

1. **Oracle** (see-saw) — for a fixed weight vector *w*, find the pair of
   states (ρ, σ) that maximizes Σᵢ wᵢ Tr(Sᵢ ρ) Tr(Sᵢ σ).
2. **Multiplicative weight update (MWU)** — iteratively update *w* to build
   an ensemble E = {(ρ_t, σ_t)} of state pairs.
3. **Linear program** — solve for the best convex combination of the
   ensemble to obtain the tightest possible δ bound.

## Requirements

### See-saw / MWU (Section IV)

- Python 3.9+
- numpy ≥ 1.21
- scipy ≥ 1.7

```bash
pip install -r requirements.txt
```

### Hamiltonian energy estimation (Section VI)

Additional dependencies for the `atom` module and `run_atom.py`:

- cvxpy, igraph, PySCF, Qiskit, Qiskit Nature

```bash
pip install -r requirements-atom.txt
```

**Note:** MOSEK requires a license. The cvxpy solver falls back to SCS, CVXOPT,
or ECOS (bundled with cvxpy) if MOSEK is unavailable.

## Quick start

Reproduce the benchmark results on 9 ℏ-imperfect graphs with 8 vertices:

```bash
python reproduce.py
```

This runs the full pipeline with default parameters (T = 10000, trials = 1)
on the 8-vertex dataset and writes results to `data/newres8.txt`.

### Options

```bash
python reproduce.py --dataset 9            # 9-vertex graphs (295 graphs, slower)
python reproduce.py --dataset all          # both 8- and 9-vertex graphs
python reproduce.py --eta auto             # compute η from √(8 ln m / T)
python reproduce.py --trials 3             # best of 3 runs per graph
python reproduce.py --dry-run              # print plan without executing
python reproduce.py --T 1000               # fewer MWU iterations (faster, less precise)
```

Full usage:

```
usage: reproduce.py [-h] [--T T] [--eta ETA] [--seesaw-max-iter N]
                    [--seesaw-tol FLOAT] [--trials N]
                    [--dataset {8,9,all}] [--dry-run]
                    [--output PATH] [--data-dir PATH]
```

### Runtime note

The 9-vertex dataset involves 512×512 density matrices (9 qubits) and
T=10000 MWU iterations per graph across 295 graphs. **Expect several hours
of runtime.** The default `--dataset 8` processes only the 9 faster
8-vertex graphs (256×256 matrices). Use `--dry-run` to preview before
committing to a long run.

## Files

| File | Description |
|---|---|
| `seesaw.py` | Core module: `seesaw_oracle`, `mwu_ensemble`, `compute_lp_coefficients`, `solve_lp` |
| `reproduce.py` | CLI driver to reproduce benchmark results |
| `data/test8.txt` | 18 frustration graphs with 8 vertices (graph6 codes + Pauli encodings) |
| `data/ndelta8_3_hard.txt` | Benchmark indices for 8-vertex graphs (9 entries, 1-based into test8.txt) |
| `data/newres8.txt` | Pre-computed reference output for 8-vertex graphs |
| `data/test9.txt` | 1419 frustration graphs with 9 vertices |
| `data/ndelta9_3_hard.txt` | Benchmark indices for 9-vertex graphs (295 entries) |
| `data/newres9.txt` | Pre-computed reference output for 9-vertex graphs |
| `requirements.txt` | Python dependencies (see-saw/MWU) |
| `atom.py` | Hamiltonian variance bounds module (Section VI) |
| `run_atom.py` | CLI driver for variance bound benchmarks (Table I) |
| `requirements-atom.txt` | Additional dependencies for Hamiltonian estimation |

## Data sources and the "hard" label

The files with the suffix `_hard` contain the benchmark index sets selected
for this study. The "hard" label refers to graphs whose sample complexity
parameter δ is not trivially determined by the fractional chromatic number
— i.e., ℏ-imperfect frustration graphs for which the see-saw/MWU method is
necessary.

The old bounds (second column in the ndelta and newres files) come from the
**SDP hierarchy based on state-polynomial optimization (SPO)** described in
Section III.A of the supplemental material. These provide converging upper
bounds on δ, which the see-saw/MWU lower bounds are compared against.

## Algorithm overview

The sample complexity of Pauli shadow tomography is governed by:

```
δ({S_i}) = min_{w ∈ Δ_m} max_{ρ,σ} Σ_i w_i Tr(S_i ρ) Tr(S_i σ)
```

### Stage 1 — See-saw oracle (Eq. 49)

For a fixed weight vector *w*, solve the inner maximization:

```
β(G, w) = max_{ρ,σ} Σ_i w_i Tr(S_i ρ) Tr(S_i σ)
```

1. Initialise ρ as a random pure state.
2. Form M₁ = Σ_i w_i Tr(S_i ρ) S_i; set σ to the maximal eigenstate of M₁.
3. Form M₂ = Σ_i w_i Tr(S_i σ) S_i; set ρ to the maximal eigenstate of M₂.
4. Repeat 2–3 until convergence (tolerance ~10⁻⁵).

### Stage 2 — Multiplicative weight update (§IV.B.2)

- Start with uniform weights wᵢ = 1/m.
- For T iterations, set η = √(8 ln m / T):
  1. Call the oracle for current weights → (ρ_t, σ_t).
  2. Compute cᵢ = Tr(Sᵢ ρ_t) Tr(Sᵢ σ_t).
  3. Update wᵢ ← wᵢ exp(−η cᵢ), renormalise.
- Collect ensemble E = {(ρ_t, σ_t)}.

### Stage 3 — Linear program (Eq. 51)

Find the best convex combination of the ensemble:

```
max  μ
s.t.  M = Σ_t (a_t ρ_t⊗ρ_t + b_t σ_t⊗σ_t),
      Tr(M (S_i⊗S_i)) ≥ μ,
      a_t, b_t ≥ 0,  Σ_t (a_t + b_t) = 1.
```

The optimal μ is the δ bound. The solution (a_t, b_t) defines the measurement
strategy: for each ensemble pair, measure with probability proportional to
a_t (via ρ_t) or b_t (via σ_t).

## Data format

### Pauli encoding

Pauli strings are encoded as integer vectors:

| Integer | Pauli |
|---|---|
| 0 | I |
| 1 | X |
| 2 | Y |
| 3 | Z |

Example: `[1, 0, 3]` → X ⊗ I ⊗ Z

### Input: test{8,9}.txt

One line per graph, each a Python literal:

```
[graph6_code, old_bound, [[1,0,0], [0,1,0], ...]]
```

- `graph6_code`: graph6 string identifier for the frustration graph
- `old_bound`: upper bound on δ from the SPO-based SDP hierarchy
- Pauli list: integer-encoded Pauli strings defining the graph realization

### Input: ndelta{8,9}_3_hard.txt

A Python list of lists:

```
[[index, v1, old_bound], ...]
```

- `index`: 1-based index into the corresponding test file
- `v1`: internal parameter
- `old_bound`: SPO-based SDP upper bound (same as in test file)

### Output: newres{8,9}.txt

One line per graph:

```
[graph6_code, old_bound, new_delta]
```

- `new_delta`: lower bound on δ computed by the see-saw/MWU/LP pipeline

## Reproducing results

To reproduce the 8-vertex benchmark:

```bash
python reproduce.py                     # default: --dataset 8
```

To reproduce the 9-vertex benchmark (warning: several hours):

```bash
python reproduce.py --dataset 9
```

To verify against pre-computed output:

```bash
python reproduce.py --dataset 8
diff data/newres8.txt data/newres8_ref.txt   # compare (if backed up)
```

For faster but approximate results (e.g., for testing):

```bash
python reproduce.py --T 1000 --output data/newres8_quick.txt
```

## Hamiltonian energy estimation (Section VI)

The `atom` module and `run_atom.py` implement the variance bounds for
Hamiltonian energy estimation described in Section VI of the supplemental
material.

### What it computes

For a molecular Hamiltonian $H = \sum_i c_i S_i$ expressed in Pauli strings,
the module computes variance bounds across multiple measurement schemes:

- **Bound 1** (weighted): $\left(\sum_i |c_i| / \sqrt{T_i}\right)^2$
- **Bound 2** (maximal eigenvalue): largest eigenvalue of the effective
  Hamiltonian $\sum_{v,w} c_v c_w \frac{\sum_{e \supseteq \{v,w\}} n_e}{N_v N_w} S_v S_w$
- **JM norm²**: joint-measurement norm squared, with and without noise
  parameter optimization
- Measurement schemes: **global**, **single-qubit**, and **two-qubit**
- Fermion-to-qubit mappers: **Parity**, **Jordan–Wigner**, **Bravyi–Kitaev**

### Usage

```bash
python run_atom.py                           # H2, (2e,2o), all mappers
python run_atom.py --all                     # full Table I
python run_atom.py --molecule LiH --nele 2 --norb 3
python run_atom.py --dry-run                 # preview
```

### Module functions

`atom.py` exports general-purpose functions usable independently of
the driver:

| Function | Description |
|---|---|
| `bounds(coeffs, eds, strings)` | All variance bounds for a given Hamiltonian |
| `variance_bound_weighted(coeffs, weights, subsets)` | Weighted variance bound |
| `maximal_eigenvalue(c_list, subsets, weights, pauli_strings)` | Max-eigenvalue bound |
| `solve_minimize_cv_over_sqrt_z(c_list, I_list)` | Optimal context weights via convex optimization |
| `compare(atomstr, nele, norb, mapperId)` | End-to-end: molecule → Hamiltonian → all bounds |
| `compute_JM_norm_original(coeffs, paulis, noise_params)` | JM norm² |
| `minimize_noise_params(coeffs, paulis)` | Optimize noise parameters for JM norm |

## Citation

If you use this code, please cite:

> Z.-P. Xu, Q. Ye, J. Wang, G. Koßmann, A. Winter, and R. Schwonnek,
> "Optimal strategies for shadow tomography with limited resources" (2026).
