import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import root
from typing import Dict, List, Tuple
from collections import defaultdict

from reaction import Reaction, ReactionNetwork
from kinetics import KineticsCalculator


class GolgiModel:
    """
    Solves the steady-state mass balance equations for 4 cascaded CSTRs.
    """
    def __init__(
        self,
        network: ReactionNetwork,
        kinetics: KineticsCalculator,
        tau: np.ndarray,
        enzyme_dist: List[Dict[str, float]],
        donor_conc: Dict[str, float],
        initial_feed: np.ndarray,
    ):
        self.network = network
        self.kinetics = kinetics
        self.tau = tau
        self.enzyme_dist = enzyme_dist          # [{enzyme_name: concentration}, ...] one dict per compartment
        self.donor_conc = donor_conc            # {cosubstrate name: concentration}
        self.initial_feed = initial_feed
        self.n_structs = len(network.structures)
        self.n_compartments = len(tau)

        # ---- helpers for the sparse extended-system solver ----
        enzyme_names_set = sorted(
            {rxn.enzyme.enzyme_name for rxn in self.network.reactions.values()}
        )
        self.enzyme_names = enzyme_names_set
        self.n_enzymes = len(enzyme_names_set)
        self.enzyme_idx = {name: i for i, name in enumerate(enzyme_names_set)}

        # Per-enzyme list of reaction ids (deterministic ordering)
        self.enzyme_rxn_ids: Dict[str, List[int]] = defaultdict(list)
        for rxn_id, rxn in self.network.reactions.items():
            self.enzyme_rxn_ids[rxn.enzyme.enzyme_name].append(rxn_id)
        for name in self.enzyme_rxn_ids:
            self.enzyme_rxn_ids[name].sort()

        self.block_size = self.n_structs + self.n_enzymes
        self.n_total = self.block_size * self.n_compartments

    def _calc_single_rate(
        self,
        rxn: Reaction,
        c_j: np.ndarray,
        competitor_sums: Dict[str, float],
        enzyme_conc_j: Dict[str, float],
    ) -> float:
        """
        Computes the rate of a single reaction:
        r = kf * [Et] * [donor] * [Pi] / (Kmi * (Kmd + [donor]) * (1 + sum(Pj)/Kmj))
        """
        enzyme = rxn.enzyme
        sub_id = rxn.substrate_id
        kmi = rxn.km_value
        kf = enzyme.kf
        kmd = enzyme.Kmd

        et = enzyme_conc_j.get(enzyme.enzyme_name, 0.0)
        donor = self.donor_conc.get(enzyme.cosubstrate, 1.0)
        pi = c_j[sub_id] if 0 <= sub_id < len(c_j) else 0.0
        cs = competitor_sums.get(enzyme.enzyme_name, 1.0)

        denom = kmi * (kmd + donor) * cs
        if denom <= 0 or et <= 0:
            return 0.0
        return kf * et * donor * pi / denom

    def _calculate_net_rates(
        self,
        c_j: np.ndarray,
        competitor_sums: Dict[str, float],
        enzyme_conc_j: Dict[str, float],
        compartment_idx: int,
    ) -> np.ndarray:
        """
        Computes the net production rate for each glycoform in the specified compartment.
        Caches all reaction rates first, then aggregates by structure using
        Glycoform._consuming_reactions / _producing_reactions.
        """
        rxn_rates = {}
        for rxn_id, rxn in self.network.reactions.items():
            rxn_rates[rxn_id] = self._calc_single_rate(
                rxn, c_j, competitor_sums, enzyme_conc_j
            )

        r = np.zeros(self.n_structs)
        for i, glycan in self.network.structures.items():
            for rxn_id in glycan._consuming_reactions:
                r[i] -= rxn_rates.get(rxn_id, 0.0)
            for rxn_id in glycan._producing_reactions:
                r[i] += rxn_rates.get(rxn_id, 0.0)

        return r

    def residual(self, x: np.ndarray) -> np.ndarray:
        """
        Legacy residual F(x) on the concentration-only variable layout.
        x : shape (n_structs * n_compartments,)
        """
        c = x.reshape(self.n_compartments, self.n_structs).T
        F = np.zeros_like(c)

        for j in range(self.n_compartments):
            c_j = c[:, j]
            c_prev = c[:, j - 1] if j > 0 else self.initial_feed

            competitor_sums = self.kinetics.compute_all_competitor_sums(c_j)
            enzyme_conc_j = self.enzyme_dist[j]
            r = self._calculate_net_rates(c_j, competitor_sums, enzyme_conc_j, j)

            F[:, j] = c_j - (c_prev + self.tau[j] * r)

        return F.flatten()

    def solve(
        self,
        method: str = 'krylov',
        tol: float = 1e-6,
        max_iter: int = 100,
        verbose: bool = True,
        **kwargs,
    ) -> np.ndarray:
        """
        Legacy wrapper around scipy.optimize.root.
        """
        x0 = np.zeros(self.n_structs * self.n_compartments)
        for j in range(self.n_compartments):
            x0[j * self.n_structs : (j + 1) * self.n_structs] = self.initial_feed

        options = {'maxiter': max_iter, **kwargs}
        if verbose:
            options['disp'] = True

        result = root(self.residual, x0, method=method, tol=tol, options=options)

        if result.success:
            c = result.x.reshape(self.n_structs, self.n_compartments)
            return np.maximum(c, 0)
        raise RuntimeError(f"Steady-state solver did not converge: {result.message}")

    def solve_damped(
        self,
        tol: float = 1e-6,
        max_iter: int = 200,
        damping: float = 0.3,
        verbose: bool = True,
    ) -> np.ndarray:
        """
        Damped Picard iteration (no Jacobian, memory-friendly).
        """
        x = np.zeros(self.n_structs * self.n_compartments)
        for j in range(self.n_compartments):
            x[j * self.n_structs : (j + 1) * self.n_structs] = self.initial_feed

        for it in range(max_iter):
            F = self.residual(x)
            norm = np.linalg.norm(F)
            if verbose and it % 8 == 0:
                print(f"Iter {it}: ||F|| = {norm:.6e}")
            if norm < tol:
                if verbose:
                    print(f"Converged at iteration {it}, ||F|| = {norm:.6e}")
                break
            x = x - damping * F
            x = np.maximum(x, 0)

        return x.reshape(self.n_structs, self.n_compartments).T

    # =====================================================================
    # Constrained Newton-Raphson with auxiliary Michaelis-Menten denominator
    # variables and a sparse direct linear solver.
    #
    # Reformulation (Krambeck & Betenbaugh, 2005):
    #
    # Original reaction rate for reaction r with enzyme e and substrate s:
    #     r = alpha_r * c[s] / (1 + sum_{r' in enzyme_e} c[s_r'] / Km_r')
    #
    # where alpha_r = kf * Et_e[j] * Donor / (Km_r * (Kmd_e + Donor)) is
    # parameter-only. We promote the competitor sum for each enzyme e to an
    # auxiliary variable d_e and add one equation per enzyme per compartment
    # to define it. The reaction rate becomes
    #     r = alpha_r * c[s] / d_e
    # so each r depends on only one c and one d, keeping the Jacobian sparse.
    #
    # Variables (one block of size S+E per compartment, S=n_structs, E=n_enzymes):
    #     [c_1, ..., c_S, d_1, ..., d_E]      per compartment
    #
    # Equations (same block layout):
    #     F_c(i, j) = c[i, j] - c[i, j-1] - tau_j * r_net_i(c_j, d_j) = 0
    #     F_d(e, j) = d[e, j] - compute_all_competitor_sums(c_j)[e]   = 0
    # =====================================================================

    def _init_sparse_structures(self):
        """Precompute indexing structures and the rate coefficient alpha
        used by the sparse Newton solver. Called lazily by solve_sparse.
        """
        # Reactions that consume / produce each structure
        self._consuming_rxns: Dict[int, List[int]] = defaultdict(list)
        self._producing_rxns: Dict[int, List[int]] = defaultdict(list)
        for rxn_id, rxn in self.network.reactions.items():
            self._consuming_rxns[rxn.substrate_id].append(rxn_id)
            self._producing_rxns[rxn.product_id].append(rxn_id)

        # Per-enzyme list of reactions and unique substrate ids, plus the
        # sum of 1/Km for each (enzyme, substrate) pair. The latter is the
        # Jacobian entry for the denominator equation F_d/dc[s].
        self._enzyme_rxns_list: List[List[int]] = [[] for _ in range(self.n_enzymes)]
        self._enzyme_substrates: List[List[int]] = [[] for _ in range(self.n_enzymes)]
        self._enzyme_substrate_invkm: List[List[float]] = [[] for _ in range(self.n_enzymes)]
        for rxn_id, rxn in self.network.reactions.items():
            e_idx = self.enzyme_idx[rxn.enzyme.enzyme_name]
            self._enzyme_rxns_list[e_idx].append(rxn_id)
            if rxn.km_value <= 0:
                continue
            sub_id = rxn.substrate_id
            inv_km = 1.0 / rxn.km_value
            if sub_id in self._enzyme_substrates[e_idx]:
                # Same substrate appears multiple times (rare); accumulate
                k = self._enzyme_substrates[e_idx].index(sub_id)
                self._enzyme_substrate_invkm[e_idx][k] += inv_km
            else:
                self._enzyme_substrates[e_idx].append(sub_id)
                self._enzyme_substrate_invkm[e_idx].append(inv_km)

        # Per-enzyme set of substrates that are consumed / produced by struct i,
        # cached so the inner Jacobian loop is cheap.
        self._struct_cons_eidx: Dict[int, List[int]] = {}
        self._struct_prod_eidx: Dict[int, List[int]] = {}
        self._struct_prod_subid: Dict[int, List[int]] = {}
        for i in range(self.n_structs):
            cons = self._consuming_rxns.get(i, [])
            self._struct_cons_eidx[i] = [
                self.enzyme_idx[self.network.reactions[r].enzyme.enzyme_name]
                for r in cons
            ]
            prod = self._producing_rxns.get(i, [])
            self._struct_prod_eidx[i] = [
                self.enzyme_idx[self.network.reactions[r].enzyme.enzyme_name]
                for r in prod
            ]
            self._struct_prod_subid[i] = [
                self.network.reactions[r].substrate_id for r in prod
            ]

        # Precompute alpha[r, j] = kf * Et_e[j] * Donor / (Km_r * (Kmd_e + Donor))
        # This depends only on (reaction, compartment) parameters.
        n_rxns = len(self.network.reactions)
        self._alpha = np.zeros((n_rxns, self.n_compartments))
        for rxn_id, rxn in self.network.reactions.items():
            enzyme = rxn.enzyme
            donor = self.donor_conc.get(enzyme.cosubstrate, 1.0)
            const_denom = rxn.km_value * (enzyme.Kmd + donor)
            if const_denom <= 0:
                continue
            for j in range(self.n_compartments):
                Et = self.enzyme_dist[j].get(enzyme.enzyme_name, 0.0)
                if Et > 0:
                    self._alpha[rxn_id, j] = enzyme.kf * Et * donor / const_denom

        # Jacobian sparsity pattern (COO triplet arrays).
        self._build_sparsity_pattern()

    def _build_sparsity_pattern(self):
        """Build COO (row, col) arrays for the Jacobian non-zero pattern.

        For each compartment j the block structure is::

            F_c   =  J_cc   J_cd     c[j]    F_c
                    J_dc   J_dd     d[j]    F_d

        with J_dc coupling every substrate of an enzyme to its denominator
        equation, J_cd coupling every enzyme touching a structure to the
        structure's mass balance, and off-diagonal blocks across
        compartments only on the c -> F_c linkage to the previous compartment.
        """
        S, E, N = self.n_structs, self.n_enzymes, self.n_compartments
        block = S + E

        rows: List[int] = []
        cols: List[int] = []

        for j in range(N):
            # F_c[i, j] rows
            for i in range(S):
                row_c = j * block + i
                # dF_c/dc[i, j]
                rows.append(row_c); cols.append(j * block + i)
                # dF_c/dc[i, j-1]
                if j > 0:
                    rows.append(row_c); cols.append((j - 1) * block + i)
                # dF_c/ddenom[e, j] for every enzyme touching struct i
                involved = set(self._struct_cons_eidx.get(i, []))
                involved.update(self._struct_prod_eidx.get(i, []))
                for e in involved:
                    rows.append(row_c); cols.append(j * block + S + e)

            # F_d[e, j] rows
            for e in range(E):
                row_d = j * block + S + e
                # dF_d/ddenom[e, j]
                rows.append(row_d); cols.append(j * block + S + e)
                # dF_d/dc[s, j] for every substrate of enzyme e
                for s in self._enzyme_substrates[e]:
                    rows.append(row_d); cols.append(j * block + s)

        self._jac_rows = np.asarray(rows, dtype=np.int32)
        self._jac_cols = np.asarray(cols, dtype=np.int32)
        self._jac_nnz = len(rows)

    def _residual_extended(self, x: np.ndarray) -> np.ndarray:
        """Compute F(x) for the extended [c, d] variable layout.

        Mass-balance equations use the auxiliary denominator d[j] as the
        competitor sum when evaluating reaction rates (this keeps the
        analytical Jacobian consistent with the residual). The denominator
        equations F_d then enforce that d[e, j] matches
        ``compute_all_competitor_sums(c_j)[e]``.
        """
        S, E, N = self.n_structs, self.n_enzymes, self.n_compartments
        block = S + E
        F = np.zeros(N * block)

        for j in range(N):
            c_j = x[j * block : j * block + S]
            d_j = x[j * block + S : (j + 1) * block]

            # Denominator equations: enforce d[e, j] = competitive sum(c_j)[e]
            competitor_sums = self.kinetics.compute_all_competitor_sums(c_j)
            for e_idx, e_name in enumerate(self.enzyme_names):
                F[j * block + S + e_idx] = d_j[e_idx] - competitor_sums[e_name]

            # Mass-balance equations: rates use d_j (not competitor_sums)
            # so that dF/dc, dF/ddenom match the closed-form Jacobian.
            d_as_dict = {name: d_j[idx] for idx, name in enumerate(self.enzyme_names)}
            enzyme_conc_j = self.enzyme_dist[j]
            net_rates = self._calculate_net_rates(c_j, d_as_dict, enzyme_conc_j, j)
            c_prev = x[(j - 1) * block : (j - 1) * block + S] if j > 0 else self.initial_feed
            F[j * block : j * block + S] = c_j - c_prev - self.tau[j] * net_rates

        return F

    def _jacobian_extended(self, x: np.ndarray) -> sp.csr_matrix:
        """Analytical sparse Jacobian of the extended residual.

        With r_r = alpha_r * c[s_r] / d[e_r] and r_net_i = sum(cons) - sum(prod),
        the derivatives reduce to a few closed-form terms per row:

          dF_c/dc[i, j]      = 1 - tau * sum_{r in cons(i)} alpha_r / d[e_r]
          dF_c/dc[i, j-1]    = -1
          dF_c/ddenom[e, j]  = tau * ( sum_{r in cons(i), e_r=e} alpha_r * c[s_r]
                                       - sum_{r in prod(i), e_r=e} alpha_r * c[s_r] ) / d[e]^2

          dF_d/ddenom[e, j]  = 1
          dF_d/dc[s, j]      = -1/Km_r for each r in enzyme_e with substrate s
        """
        S, E, N = self.n_structs, self.n_enzymes, self.n_compartments
        block = S + E
        data = np.zeros(self._jac_nnz)
        cursor = 0

        for j in range(N):
            c_j = x[j * block : j * block + S]
            d_j = x[j * block + S : (j + 1) * block]
            tau_j = self.tau[j]

            # F_c[i, j] rows
            for i in range(S):
                # dF_c/dc[i, j]
                val = 1.0
                for k, r in enumerate(self._consuming_rxns.get(i, [])):
                    e_idx = self._struct_cons_eidx[i][k]
                    d_e = d_j[e_idx]
                    if d_e > 1e-30:
                        val += tau_j * self._alpha[r, j] / d_e
                data[cursor] = val
                cursor += 1

                # dF_c/dc[i, j-1]
                if j > 0:
                    data[cursor] = -1.0
                    cursor += 1

                # dF_c/ddenom[e, j]
                cons_eidx = self._struct_cons_eidx.get(i, [])
                prod_eidx = self._struct_prod_eidx.get(i, [])
                prod_subid = self._struct_prod_subid.get(i, [])
                involved = set(cons_eidx)
                involved.update(prod_eidx)
                for e_idx in involved:
                    d_e = d_j[e_idx]
                    if d_e <= 1e-30:
                        data[cursor] = 0.0
                        cursor += 1
                        continue
                    inv_d2 = 1.0 / (d_e * d_e)
                    val = 0.0
                    for k, r in enumerate(self._consuming_rxns.get(i, [])):
                        if cons_eidx[k] == e_idx:
                            val -= tau_j * self._alpha[r, j] * c_j[i] * inv_d2
                    for k, r in enumerate(self._producing_rxns.get(i, [])):
                        if prod_eidx[k] == e_idx:
                            val += tau_j * self._alpha[r, j] * c_j[prod_subid[k]] * inv_d2
                    data[cursor] = val
                    cursor += 1

            # F_d[e, j] rows
            for e_idx in range(E):
                # dF_d/ddenom[e, j]
                data[cursor] = 1.0
                cursor += 1

                # dF_d/dc[s, j]
                inv_km_row = self._enzyme_substrate_invkm[e_idx]
                for k in range(len(self._enzyme_substrates[e_idx])):
                    data[cursor] = -inv_km_row[k]
                    cursor += 1

        J = sp.coo_matrix(
            (data, (self._jac_rows, self._jac_cols)),
            shape=(self.n_total, self.n_total),
        ).tocsr()
        return J

    def solve_sparse(
        self,
        tol: float = 1e-6,
        max_iter: int = 50,
        damping: float = 1.0,
        line_search: bool = True,
        min_step: float = 0.01,
        verbose: bool = True,
    ) -> np.ndarray:
        """Constrained Newton-Raphson with sparse direct linear solve.

        Reformulates the steady-state system by promoting each Michaelis-Menten
        competitor sum to an auxiliary variable (one per enzyme per compartment),
        so each reaction rate depends only on its substrate concentration and
        that one denominator. This keeps the Jacobian sparse even for very
        large networks.

        Parameters
        ----------
        tol : float
            Convergence tolerance on the L2 norm of the residual.
        max_iter : int
            Maximum Newton iterations.
        damping : float
            Initial step size multiplier; halved during line search until the
            residual norm stops growing.
        line_search : bool
            If True, perform a backtracking line search that rejects any step
            that does not strictly decrease ||F||.
        min_step : float
            Smallest step size the line search will try before giving up.
        verbose : bool
            Print progress every few iterations.

        Returns
        -------
        np.ndarray
            Concentration matrix of shape (n_structs, n_compartments).
        """
        self._init_sparse_structures()

        S, E, N = self.n_structs, self.n_enzymes, self.n_compartments
        block = S + E

        # ---- initial guess ----
        # c: broadcast the feed through every compartment;
        # d: set to the baseline competitor sum so F_d starts at ~0.
        x = np.zeros(block * N)
        baseline_comp = self.kinetics.compute_all_competitor_sums(self.initial_feed)
        baseline_d = np.array([baseline_comp[name] for name in self.enzyme_names])
        for j in range(N):
            x[j * block : j * block + S] = self.initial_feed
            x[j * block + S : (j + 1) * block] = baseline_d

        # ---- Newton iteration ----
        converged = False
        norm = np.inf
        for it in range(max_iter):
            F = self._residual_extended(x)
            norm = np.linalg.norm(F)
            if verbose and (it < 5 or it % 5 == 0):
                print(f"Newton iter {it}: ||F|| = {norm:.6e}")
            if norm < tol:
                converged = True
                if verbose:
                    print(f"Converged at iter {it}, ||F|| = {norm:.6e}")
                break

            J = self._jacobian_extended(x)
            try:
                dx = spla.spsolve(J, -F)
            except RuntimeError as e:
                if verbose:
                    print(f"spsolve failed at iter {it}: {e}; using lsqr fallback")
                dx = spla.lsqr(J, -F, atol=tol * 0.1, btol=tol * 0.1)[0]

            # Damped Newton step with backtracking line search.
            # Clamp only the concentration portion to >= 0; the denominator
            # variables d_j stay as the solver chooses.
            step = damping
            best_step = 0.0
            best_norm = norm
            best_x = x.copy()
            while step >= min_step:
                x_trial = x + step * dx
                x_trial[: N * S] = np.maximum(x_trial[: N * S], 0.0)
                norm_trial = np.linalg.norm(self._residual_extended(x_trial))
                if norm_trial < best_norm:
                    best_step = step
                    best_norm = norm_trial
                    best_x = x_trial
                    break  # first improvement is good enough; finer steps rarely help
                step *= 0.5

            if best_step == 0.0:
                if verbose:
                    print(f"Newton iter {it}: line search failed at ||F|| = {norm:.6e}")
                raise RuntimeError(
                    f"Line search could not find a descent step at iter {it} "
                    f"(||F|| = {norm:.6e})"
                )

            x = best_x

        if not converged:
            if verbose:
                print(f"Did not converge in {max_iter} iterations, ||F|| = {norm:.6e}")
            raise RuntimeError(
                f"Sparse Newton solver did not converge in {max_iter} iterations "
                f"(final ||F|| = {norm:.6e})"
            )

        concs = np.zeros((S, N))
        for j in range(N):
            concs[:, j] = x[j * block : j * block + S]
        return concs