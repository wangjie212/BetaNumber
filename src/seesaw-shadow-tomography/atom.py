import numpy as np
import itertools
from functools import reduce
from collections import defaultdict
from typing import List, Iterable

import igraph as ig
import cvxpy as cp
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix, kron
from scipy.optimize import minimize

from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.mappers import ParityMapper, JordanWignerMapper, BravyiKitaevMapper
from qiskit.quantum_info import SparsePauliOp
# ============================================================
# Shared Pauli primitives
# ============================================================

PAULI_DENSE = [
    np.array([[1, 0], [0, 1]], dtype=complex),   # I
    np.array([[0, 1], [1, 0]], dtype=complex),   # X
    np.array([[0, -1j], [1j, 0]], dtype=complex),# Y
    np.array([[1, 0], [0, -1]], dtype=complex),  # Z
]

PAULI_SPARSE = [csr_matrix(p) for p in PAULI_DENSE]

PRODUCT_TABLE = {
    (0,0):(0, 1),  (0,1):(1, 1),   (0,2):(2, 1),   (0,3):(3, 1),
    (1,0):(1, 1),  (1,1):(0, 1),   (1,2):(3, 1j),  (1,3):(2,-1j),
    (2,0):(2, 1),  (2,1):(3,-1j),  (2,2):(0, 1),   (2,3):(1, 1j),
    (3,0):(3, 1),  (3,1):(2, 1j),  (3,2):(1,-1j),  (3,3):(0, 1)
}

# ============================================================
# Pauli utilities
# ============================================================

def pauli_string_to_sparse(pauli):
    mat = csr_matrix(np.array([[1.0]], dtype=complex))
    for p in pauli:
        mat = kron(mat, PAULI_SPARSE[p], format='csr')
    return mat

def pauli_matrix(pauli_str):
    result = PAULI_DENSE[pauli_str[0]]
    for p in pauli_str[1:]:
        result = np.kron(result, PAULI_DENSE[p])
    return result

def mul_pauli(p1, p2):
    """Multiply two Pauli strings using PRODUCT_TABLE."""
    res = [0] * len(p1)
    phase = 1.0 + 0.0j
    for i, (a, b) in enumerate(zip(p1, p2)):
        c, ph = PRODUCT_TABLE[(a, b)]
        res[i] = c
        phase *= ph
    return res, phase

# ============================================================
# Graph / string utilities
# ============================================================

def pairwise_map_product(list1, list2):
    if len(list1) != len(list2):
        raise ValueError("两个列表长度必须相同")
    result = 1
    for i, j in zip(list1, list2):
        if i != 0 and j != 0 and i != j:
            result *= -1
    return result

def gheds(strings):
    return [(i, j) for i in range(len(strings))
            for j in range(i+1, len(strings))
            if pairwise_map_product(strings[i], strings[j]) == -1]

def find_maximal_independent_sets(edges: list[tuple], vertex_range: range | list):
    vertices = list(vertex_range)
    unknown = {v for e in edges for v in e} - set(vertices)
    if unknown:
        raise ValueError(f"Edge endpoints not in vertex_range: {unknown}")
    label_to_idx = {v: i for i, v in enumerate(vertices)}
    idx_edges = [(label_to_idx[u], label_to_idx[v]) for u, v in edges]
    G = ig.Graph(n=len(vertices), edges=idx_edges)
    mis_list = [set(vertices[i] for i in mis) for mis in G.maximal_independent_vertex_sets()]
    return sorted(mis_list, key=lambda s: (-len(s), sorted(s)))

def labels1(eds, vs):
    (i, j) = eds[0]
    label = []
    for k in vs:
        if k == i:
            label.append(1)
        elif k == j:
            label.append(3)
        else:
            ti = ([k, i] in eds or [i, k] in eds)
            tj = ([k, j] in eds or [j, k] in eds)
            label.append([0, 1, 3, 2][ti*2 + tj])
    return label

def labels(preeds, vs, lbs):
    eds = [list(it) for it in preeds]
    if len(eds) == 0:
        return [[lb[i] for lb in lbs] for i in range(len(vs))]
    lb = labels1(eds, vs)
    (i, j) = eds[0]
    newlbs = lbs + [lb]
    neweds = []
    for i2 in vs:
        for j2 in vs[i2+1::]:
            t1 = [i2, j2] in eds or [j2, i2] in eds
            t2 = lb[i2]*lb[j2] != 0 and lb[i2] != lb[j2]
            if t1 + t2 == 1:
                neweds.append([i2, j2])
    return labels(neweds, vs, newlbs)

def minRep(strings):
    eds = gheds(strings)
    return labels(eds, range(len(strings)), [])

def toOps(lb):
    return [reduce(np.kron, [PAULI_DENSE[i] for i in li]) for li in lb]

# ============================================================
# Noise / Hamiltonian utilities
# ============================================================

def eta_P(pauli, noise_params=None):
    if noise_params is None:
        return (1 / np.sqrt(3)) ** sum(1 for p in pauli if p != 0)
    prod = 1.0
    for i, p in enumerate(pauli):
        if p == 1:   prod *= noise_params[i]['x']
        elif p == 2: prod *= noise_params[i]['y']
        elif p == 3: prod *= noise_params[i]['z']
    return prod

def f_PQ(P, Q):
    for pi, qi in zip(P, Q):
        if pi != 0 and qi != 0 and pi != qi:
            return 0
    return 1

def eta_PQ_original(P, Q, noise_params=None):
    prod_pauli = [PRODUCT_TABLE[(pi, qi)][0] for pi, qi in zip(P, Q)]
    return eta_P(prod_pauli, noise_params)

def compute_JM_norm_original(coeffs, paulis, noise_params=None):
    n = len(paulis[0])
    dim = 1 << n
    H_eff_coeffs = {}

    for coeff_P, P in zip(coeffs, paulis):
        etaP = eta_P(P, noise_params)
        for coeff_Q, Q in zip(coeffs, paulis):
            if not f_PQ(P, Q):
                continue
            etaQ = eta_P(Q, noise_params)
            etaPQ_val = eta_PQ_original(P, Q, noise_params)
            kappa = (etaPQ_val / (etaP * etaQ)) * coeff_P * coeff_Q

            prod_pauli, phase = [], 1.0 + 0j
            for pi, qi in zip(P, Q):
                res_type, ph = PRODUCT_TABLE[(pi, qi)]
                prod_pauli.append(res_type)
                phase *= ph

            key = tuple(prod_pauli)
            H_eff_coeffs[key] = H_eff_coeffs.get(key, 0) + kappa * phase

    H_sparse = csr_matrix((dim, dim), dtype=complex)
    for pauli, coeff in H_eff_coeffs.items():
        if abs(coeff) > 1e-12:
            H_sparse += coeff * pauli_string_to_sparse(pauli)

    H_sparse = (H_sparse + H_sparse.conj().T) * 0.5

    if dim <= 64:
        eigvals = np.linalg.eigvalsh(H_sparse.toarray())
        eigval_maxabs = np.max(np.abs(eigvals))
    else:
        eigvals = eigsh(H_sparse, k=1, which='LM', return_eigenvectors=False)
        eigval_maxabs = abs(eigvals[0])

    return np.sqrt(eigval_maxabs)

# ============================================================
# Shared helper: compute vertex degrees N[v] and pairwise sums
# ============================================================

def _compute_N_and_pair_sum(n_vertices, subsets, weights=None):
    N = [0.0] * n_vertices
    pair_sum = defaultdict(float)
    for e, w_e in zip(subsets, (weights if weights is not None else itertools.repeat(1))):
        e_list = list(e)
        for v in e_list:
            N[v] += w_e
        k = len(e_list)
        for i in range(k):
            for j in range(i + 1, k):
                v, w = e_list[i], e_list[j]
                if v > w:
                    v, w = w, v
                pair_sum[(v, w)] += w_e
    return N, pair_sum

# ============================================================
# Variance bound functions
# ============================================================

def variance_bound_weighted(coeffs: List[float], weights: List[float], subsets: Iterable[Iterable[int]]) -> float:
    subsets_list = list(subsets)
    if len(weights) != len(subsets_list):
        raise ValueError("Number of weights must equal number of subsets")
    n_vertices = len(coeffs)
    N, pair_sum = _compute_N_and_pair_sum(n_vertices, subsets_list, weights)
    for v in range(n_vertices):
        if N[v] == 0.0 and coeffs[v] != 0.0:
            raise ValueError(f"Vertex {v} has N_v=0 but c_v={coeffs[v]} != 0 → expression is infinite")
    total = sum(coeffs[v]**2 / N[v] for v in range(n_vertices) if N[v] != 0.0)
    for (v, w), s in pair_sum.items():
        if N[v] != 0.0 and N[w] != 0.0:
            total += 2.0 * coeffs[v] * coeffs[w] * s / (N[v] * N[w])
    return total

# ============================================================
# Maximal eigenvalue (kept as requested)
# ============================================================

def maximal_eigenvalue(c_list, subsets, weights, pauli_strings):
    """
    Compute the maximal eigenvalue of
        H = sum_{v,w} c_v c_w * (sum_{e supseteq {v,w}} n_e)/(N_v N_w) * S_v S_w

    Parameters
    ----------
    c_list        : list of float
    subsets       : list of list of int
    weights       : list of float
    pauli_strings : list of list of int (entries in {0,1,2,3})

    Returns
    -------
    float : largest eigenvalue of H
    """
    n_vertices = len(c_list)
    if len(pauli_strings) != n_vertices:
        raise ValueError("len(pauli_strings) must equal len(c_list)")
    if len(subsets) != len(weights):
        raise ValueError("Number of subsets must equal number of weights")

    N, pair_sum = _compute_N_and_pair_sum(n_vertices, subsets, weights)

    for v in range(n_vertices):
        if N[v] == 0.0 and c_list[v] != 0.0:
            raise ValueError(f"Vertex {v} has N_v=0 but c_v={c_list[v]} != 0")

    L = len(pauli_strings[0])
    op_dict = defaultdict(complex)

    # Identity (diagonal) terms v == w
    identity_term = sum(
        c_list[v]**2 / N[v] for v in range(n_vertices) if N[v] != 0.0
    )
    if identity_term != 0.0:
        op_dict[tuple([0]*L)] += identity_term

    # Off-diagonal pairs v < w
    for v in range(n_vertices):
        if N[v] == 0.0:
            continue
        for w in range(v+1, n_vertices):
            if N[w] == 0.0:
                continue
            s = pair_sum.get((v, w), 0.0)
            if s == 0.0:
                continue
            M = c_list[v] * c_list[w] * s / (N[v] * N[w])
            prod_str, phase = mul_pauli(pauli_strings[v], pauli_strings[w])
            real_phase = phase.real
            if real_phase != 0.0:
                op_dict[tuple(prod_str)] += 2.0 * M * real_phase

    dim = 1 << L
    H = np.zeros((dim, dim), dtype=complex)
    for ps, coeff in op_dict.items():
        if coeff != 0.0:
            H += coeff * pauli_matrix(ps)

    eigvals = np.linalg.eigvalsh(H)
    return float(np.max(eigvals.real))

# ============================================================
# Optimisation
# ============================================================

def solve_minimize_cv_over_sqrt_z(c_list, I_list, eps=1e-8):
    c_list = np.real_if_close(c_list).astype(float)
    if np.any(np.iscomplex(c_list)):
        raise ValueError("c_list contains non-negligible imaginary parts. Please provide real coefficients.")
    if np.any(c_list < 0):
        raise ValueError("This formulation requires all coefficients c_v >= 0. Please adjust your data or take absolute values.")
    n_elem = len(c_list)
    n_sets = len(I_list)
    for i, I in enumerate(I_list):
        for v in I:
            if not (0 <= v < n_elem):
                raise ValueError(f"Element {v} in I_list[{i}] is out of range [0, {n_elem-1}]")
    x = cp.Variable(n_sets, nonneg=True)
    y = cp.Variable(n_elem, nonneg=True)
    constraints = [cp.sum(x) == 1]
    for v in range(n_elem):
        indices = [i for i, I in enumerate(I_list) if v in I]
        if indices:
            constraints.append(cp.square(y[v]) <= cp.sum([x[i] for i in indices]))
        else:
            constraints.append(y[v] == 0)
    objective = cp.Minimize(cp.sum([c_list[v] * cp.inv_pos(y[v] + eps) for v in range(n_elem)]))
    problem = cp.Problem(objective, constraints)
    solved = False
    for solver_name in ['MOSEK', 'SCS', 'CVXOPT', 'ECOS']:
        try:
            problem.solve(solver=solver_name, verbose=False)
            if problem.status in ["optimal", "optimal_inaccurate"]:
                solved = True
                break
        except cp.SolverError:
            continue
        except Exception as e:
            print(f"Solver {solver_name} failed: {e}")
            continue
    if not solved:
        raise RuntimeError("No suitable solver found. Please install SCS, CVXOPT, or ECOS.")
    if problem.status not in ["optimal", "optimal_inaccurate"]:
        print(f"Warning: Solver status = {problem.status}")
    return x.value, y.value**2, problem

def strategy(eds, coeffs):
    n = len(coeffs)
    miss = find_maximal_independent_sets(eds, range(n))
    c_list_abs = np.abs(coeffs)
    _, _, prob = solve_minimize_cv_over_sqrt_z(c_list_abs, miss)
    return prob.value**2

# ============================================================
# Noise parameter optimisation for JM norm
# ============================================================

def _noise_params_to_vector(noise_params):
    vec = []
    for d in noise_params:
        x, y, z = d['x'], d['y'], d['z']
        r = np.sqrt(x*x + y*y + z*z)
        if r < 1e-12:
            theta = 0.0
            phi = 0.0
        else:
            theta = np.arccos(z / r)
            phi = np.arctan2(y, x)
        vec.extend([r, theta, phi])
    return np.array(vec)

def _vector_to_noise_params(vec, n_qubits):
    noise_params = []
    for i in range(n_qubits):
        r = vec[3*i]
        theta = vec[3*i+1]
        phi = vec[3*i+2]
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        noise_params.append({'x': x, 'y': y, 'z': z})
    return noise_params

def _objective(flat_vec, coeffs, paulis, n_qubits):
    noise_params = _vector_to_noise_params(flat_vec, n_qubits)
    return compute_JM_norm_original(coeffs, paulis, noise_params)

def minimize_noise_params(coeffs, paulis, initial_noise_params=None, method='L-BFGS-B'):
    n_qubits = len(paulis[0])
    if initial_noise_params is None:
        init_val = 1.0 / np.sqrt(3)
        initial_noise_params = [{'x': init_val, 'y': init_val, 'z': init_val} for _ in range(n_qubits)]
    else:
        if len(initial_noise_params) != n_qubits:
            raise ValueError("len(initial_noise_params) must equal number of qubits")
    x0 = _noise_params_to_vector(initial_noise_params)
    bounds = []
    for _ in range(n_qubits):
        bounds.extend([(0.0, 1.0), (0.0, np.pi), (None, None)])
    res = minimize(_objective, x0, args=(coeffs, paulis, n_qubits),
                   method=method, bounds=bounds, options={'maxiter': 500, 'ftol': 1e-8})
    if not res.success:
        print(f"Warning: optimisation did not converge: {res.message}")
    return _vector_to_noise_params(res.x, n_qubits), res.fun

def compute_diagonal_sum(coeffs, paulis, noise_params=None):
    total = 0.0
    for c, P in zip(coeffs, paulis):
        eta = eta_P(P, noise_params)
        if eta == 0.0:
            return float('inf')
        total += (c * c) / (eta * eta)
    return total

def minimize_diagonal_sum(coeffs, paulis, initial_noise_params=None, method='L-BFGS-B'):
    n_qubits = len(paulis[0])
    if initial_noise_params is None:
        init_val = 1.0 / np.sqrt(3)
        initial_noise_params = [{'x': init_val, 'y': init_val, 'z': init_val} for _ in range(n_qubits)]
    x0 = _noise_params_to_vector(initial_noise_params)
    bounds = [(0.0, 1.0), (0.0, np.pi), (None, None)] * n_qubits
    def obj(flat_vec):
        return compute_diagonal_sum(coeffs, paulis, _vector_to_noise_params(flat_vec, n_qubits))
    res = minimize(obj, x0, method=method, bounds=bounds, options={'maxiter': 500})
    return _vector_to_noise_params(res.x, n_qubits), res.fun


def bounds(coeffs,eds,strings):
    n = len(coeffs)
    miss = find_maximal_independent_sets(eds,range(n))
    c_list_abs = np.abs(coeffs)
    x_opt, z_opt, prob = solve_minimize_cv_over_sqrt_z(c_list_abs, miss)
    b1NoOpt = variance_bound_weighted(c_list_abs, [1/len(miss)]*len(miss), miss) # no optimize of the coefficients of each context
    b1Opt = variance_bound_weighted(c_list_abs, x_opt, miss) # it works finally!!!
    b2NoOpt = maximal_eigenvalue(coeffs, miss, [1/len(miss)]*len(miss), strings)
    b2Opt = maximal_eigenvalue(coeffs, miss, x_opt, strings)
    print("bound 1 No Opt, bound 1 Opt, bound 2 No Opt, bound 2 Opt")
    return [b1NoOpt,b1Opt,b2NoOpt,b2Opt]


def compare(atomstr,nele,norb,mapperId=0):
    driver = PySCFDriver(
    atom=atomstr,
    unit=DistanceUnit.ANGSTROM,
    basis="sto3g",)
    problem = driver.run()
    # Active space: 4 electrons in 4 spatial orbitals (valence space for BeH₂)
    transformer = ActiveSpaceTransformer(num_electrons=nele, num_spatial_orbitals=norb)
    reduced_problem = transformer.transform(problem)
    second_q_op = reduced_problem.hamiltonian.second_q_op()
    # Select mapper based on mapperId
    if mapperId == 0:
        mapper = ParityMapper()
    elif mapperId == 1:
        mapper = JordanWignerMapper()
    elif mapperId == 2:
        mapper = BravyiKitaevMapper()         
    hamiltonian = mapper.map(second_q_op)
    hamiltonian_simp = SparsePauliOp(hamiltonian).simplify()
    #print(f"Number of Pauli terms: {len(hamiltonian_simp)}")
    terms = hamiltonian_simp.paulis[1:]
    coeffs = hamiltonian_simp.coeffs[1:]
    char_to_num = {'I': 0, 'X': 1, 'Y': 2, 'Z': 3}
    strings = [[char_to_num.get(i) for i in str(it)] for it in terms]
    eds = gheds(strings)
    #print(f"Length of Pauli strings: {len(strings[0])}")
    n = len(coeffs)
    miss = find_maximal_independent_sets(eds,range(n))
    c_list_abs = np.abs(coeffs)
    x_opt, z_opt, prob = solve_minimize_cv_over_sqrt_z(c_list_abs, miss)
    b1Opt = variance_bound_weighted(c_list_abs, x_opt, miss) 
    b2Opt = maximal_eigenvalue(coeffs, miss, x_opt, strings)
    # measurements on single qubits
    stringss = [[[it[k]] for it in strings] for k in range(len(strings[0]))]
    alleds = sum((gheds(s) for s in stringss), [])
    miss2 = find_maximal_independent_sets(alleds,range(n))
    x_opt, z_opt, prob = solve_minimize_cv_over_sqrt_z(c_list_abs, miss2)
    b1Opt1 = variance_bound_weighted(c_list_abs, x_opt, miss2)
    b2Opt1 = maximal_eigenvalue(coeffs, miss2, x_opt, strings)
    bJMNoOpt = compute_JM_norm_original(coeffs, strings)**2
    nparas = minimize_diagonal_sum(np.real(coeffs), strings, initial_noise_params=None, method='L-BFGS-B')[0]
    bJMd = compute_JM_norm_original(coeffs, strings, nparas)**2
    l = len(strings[0])
    parts = [(2*i,2*i+1) for i in range(int(l/2))]
    stringss = [[[it[i] for i in part] for it in strings] for part in parts]
    alleds = sum((gheds(s) for s in stringss), [])
    miss2 = find_maximal_independent_sets(alleds,range(n))
    x_opt, z_opt, prob = solve_minimize_cv_over_sqrt_z(c_list_abs, miss2)
    b1Opt2 = variance_bound_weighted(c_list_abs, x_opt, miss2)
    b2Opt2 = maximal_eigenvalue(coeffs, miss2, x_opt, strings)
    #print("bound 1 Opt Global, bound 2 Opt Global, bound 1 Opt Single, bound 2 Opt Single, bound JM No Opt, bound JM Opt, bound 1 Opt two, bound 2 Opt Two")
    return [[b1Opt,b2Opt],[b1Opt1,b2Opt1],[bJMNoOpt,bJMd],[b1Opt2,b2Opt2]]