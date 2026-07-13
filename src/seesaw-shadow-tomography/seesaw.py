"""
Seesaw + multiplicative weight update algorithm for shadow tomography.

Implements the three-stage pipeline from Section IV.B of the supplemental
material for "Optimal strategies for shadow tomography with limited resources":

  1. seesaw_oracle    — See-saw optimisation for a fixed weight vector (Eq. 49)
  2. mwu_ensemble     — Multiplicative weight update outer loop (§IV.B.2)
  3. solve_lp         — Linear program for the best convex combination (Eq. 51)

Notation:
  {S_1, ..., S_m}  — Pauli strings (n-qubit)
  w ∈ Δ_m          — probability weights
  β(G, w)          — weighted beta number
  δ({S_i})         — sample complexity parameter
  E = {(ρ_t, σ_t)} — ensemble of state pairs
"""

import numpy as np
from functools import reduce
from typing import List, Tuple
from scipy.optimize import linprog

# ----------------------------------------------------------------------
# Pauli and state helpers
# ----------------------------------------------------------------------

def pauli_from_ints(vec: List[int]) -> np.ndarray:
    """Convert integer-encoded Pauli string to tensor product matrix.

    0 -> I, 1 -> X, 2 -> Y, 3 -> Z

    Example: [1, 0, 3] -> X ⊗ I ⊗ Z
    """
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    pauli_map = {0: I, 1: X, 2: Y, 3: Z}
    matrices = [pauli_map[v] for v in vec]
    return reduce(np.kron, matrices)


def random_pure_state(dim: int) -> np.ndarray:
    """Random pure state density matrix |ψ⟩⟨ψ|."""
    psi = np.random.randn(dim) + 1j * np.random.randn(dim)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def max_eigenstate(mat: np.ndarray) -> np.ndarray:
    """Density matrix of the eigenvector with largest eigenvalue."""
    eigvals, eigvecs = np.linalg.eigh(mat)
    return np.outer(eigvecs[:, np.argmax(eigvals)],
                    eigvecs[:, np.argmax(eigvals)].conj())


# ----------------------------------------------------------------------
# Stage 1 — Oracle: see-saw for a fixed weight vector  (Eq. 49)
# ----------------------------------------------------------------------

def seesaw_oracle(pauli_strings: List[List[int]], w: List[float],
                  max_iter: int = 100, tol: float = 1e-5):
    """See-saw optimisation for a fixed probability vector w.

    Solves the weighted beta number:
        β(G, w) = max_{ρ,σ} Σ_i w_i Tr(S_i ρ) Tr(S_i σ)

    Oracle steps:
      1. Initialise ρ as a random pure state.
      2. M₁ = Σ_i w_i Tr(S_i ρ) S_i  →  σ = max eigenstate of M₁
      3. M₂ = Σ_i w_i Tr(S_i σ) S_i  →  ρ = max eigenstate of M₂
      4. Repeat 2–3 until convergence.

    Args:
        pauli_strings: list of integer-encoded Pauli vectors
            e.g. [[1,0,0], [0,1,0], [3,0,0], ...]
        w: probability weight vector (length m = len(pauli_strings))
        max_iter: maximum see-saw iterations
        tol: convergence tolerance on objective change

    Returns:
        (rho, sigma, objective) where objective = Σ_i w_i Tr(S_i ρ) Tr(S_i σ)
    """
    S = [pauli_from_ints(s) for s in pauli_strings]
    m = len(S)
    dim = S[0].shape[0]

    rho = random_pure_state(dim)
    sigma = None
    prev_obj = -np.inf

    for _ in range(max_iter):
        # Step 2: update sigma
        rho_exps = np.array([np.trace(S[i] @ rho).real for i in range(m)])
        M = sum(w[i] * rho_exps[i] * S[i] for i in range(m))
        sigma = max_eigenstate(M)

        # Step 3: update rho
        sigma_exps = np.array([np.trace(S[i] @ sigma).real for i in range(m)])
        N = sum(w[i] * sigma_exps[i] * S[i] for i in range(m))
        rho_new = max_eigenstate(N)

        # Objective value
        objective = sum(w[i] * rho_exps[i] * sigma_exps[i] for i in range(m))

        if abs(objective - prev_obj) < tol:
            rho = rho_new
            break
        prev_obj = objective
        rho = rho_new

    # Final objective with converged states
    final_rho_exps = np.array([np.trace(S[i] @ rho).real for i in range(m)])
    final_sigma_exps = np.array([np.trace(S[i] @ sigma).real for i in range(m)])
    final_objective = sum(w[i] * final_rho_exps[i] * final_sigma_exps[i]
                          for i in range(m))

    return rho, sigma, final_objective


# ----------------------------------------------------------------------
# Stage 2 — Multiplicative weight update  (§IV.B.2)
# ----------------------------------------------------------------------

def mwu_ensemble(pauli_strings: List[List[int]], T: int, eta: float,
                 seesaw_max_iter: int = 100, seesaw_tol: float = 1e-5):
    """Multiplicative weight update outer loop.

    For t = 1..T:
      1. Call seesaw_oracle with current weights w to get (ρ_t, σ_t).
      2. Compute c_i = Tr(S_i ρ_t) Tr(S_i σ_t).
      3. Update w_i ← w_i exp(-η c_i), renormalise.

    Args:
        pauli_strings: list of integer-encoded Pauli vectors
        T: number of MWU iterations
        eta: learning rate (paper recommends η = √(8 ln m / T))
        seesaw_max_iter: oracle inner iterations
        seesaw_tol: oracle convergence tolerance

    Returns:
        ensemble: list of (rho, sigma) pairs, length T
    """
    m = len(pauli_strings)

    # Start with uniform weights
    w = np.ones(m) / m
    ensemble = []

    for _ in range(T):
        rho, sigma, _ = seesaw_oracle(pauli_strings, w,
                                      max_iter=seesaw_max_iter,
                                      tol=seesaw_tol)
        ensemble.append((rho, sigma))

        # Correlation coefficients
        S = [pauli_from_ints(s) for s in pauli_strings]
        rho_exps = np.array([np.trace(S[i] @ rho).real for i in range(m)])
        sigma_exps = np.array([np.trace(S[i] @ sigma).real for i in range(m)])
        c = rho_exps * sigma_exps

        # Multiplicative weight update
        w *= np.exp(-eta * c)
        w /= np.sum(w)

    return ensemble


# ----------------------------------------------------------------------
# Stage 3 — Linear program for the best convex combination  (Eq. 51)
# ----------------------------------------------------------------------

def compute_lp_coefficients(ensemble_rho, ensemble_sigma, pauli_strings):
    """Extract LP coefficients from the ensemble.

    Computes e_ρ[t,i] = (Tr(S_i ρ_t))²  and  e_σ[t,i] = (Tr(S_i σ_t))².

    Returns:
        (e_rho, e_sigma, m)  where each is shape (T, m)
    """
    S = [pauli_from_ints(s) for s in pauli_strings]
    m = len(S)
    T = len(ensemble_rho)
    e_rho = np.zeros((T, m))
    e_sigma = np.zeros((T, m))
    for t, (rho, sigma) in enumerate(zip(ensemble_rho, ensemble_sigma)):
        for i, Si in enumerate(S):
            e_rho[t, i] = (np.trace(Si @ rho).real) ** 2
            e_sigma[t, i] = (np.trace(Si @ sigma).real) ** 2
    return e_rho, e_sigma, m


def solve_lp(e_rho: np.ndarray, e_sigma: np.ndarray):
    """Solve LP for the best convex combination (Eq. 51).

        max  μ
        s.t. M = Σ_t (a_t ρ_t⊗ρ_t + b_t σ_t⊗σ_t),
             Tr(M (S_i⊗S_i)) ≥ μ,
             a_t, b_t ≥ 0,  Σ_t (a_t + b_t) = 1.

    Returns:
        (δ_opt, a_opt, b_opt) — optimal δ, a-weights, b-weights
    """
    T, m = e_rho.shape
    n_vars = 2 * T + 1          # a_1..a_T, b_1..b_T, δ
    c = np.zeros(n_vars)
    c[-1] = -1.0                # minimise -δ  →  maximise δ

    # Inequality constraints: Σ_t (a_t e_ρ[t,i] + b_t e_σ[t,i]) ≥ δ
    A_ub = np.zeros((m, n_vars))
    for i in range(m):
        for t in range(T):
            A_ub[i, t] = -e_rho[t, i]
            A_ub[i, T + t] = -e_sigma[t, i]
        A_ub[i, -1] = 1.0
    b_ub = np.zeros(m)

    # Equality: Σ_t (a_t + b_t) = 1
    A_eq = np.zeros((1, n_vars))
    for t in range(T):
        A_eq[0, t] = 1.0
        A_eq[0, T + t] = 1.0
    b_eq = np.array([1.0])

    bounds = [(0, None)] * (2 * T) + [(None, None)]
    result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                     bounds=bounds, method='highs')
    if not result.success:
        raise RuntimeError(f"LP failed: {result.message}")

    x = result.x
    return x[-1], x[:T], x[T:2 * T]
