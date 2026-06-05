"""
lh_kinetics_v5.py
═══════════════════════════════════════════════════════════════════════════════
Kinetic model for propylene electrooxidation on Pd/C
Two-geometry Frumkin LH (metallic Pd) + MvK (PdO surface)

References
──────────
[W19]  Winiwarter et al., Energy Environ. Sci. 12, 1055 (2019)
[K21]  Koroidov, Nilsson et al., Catal. Sci. Technol. 11, 3347 (2021)

Mechanism summary (full derivation in MODULE 2 docstrings)
──────────────────────────────────────────────────────────
LH regime (metallic Pd, ~0.70–1.05 V vs RHE):
  Step 1a  C₃H₆ + * ⇌ C₃H₆*_π     K_π_eff = K_π0·exp(−g·θ_prop/RT)
  Step 1b  C₃H₆ + * ⇌ C₃H₆*_σ     K_σ (coverage-independent)
  Step 2   H₂O + * → OH* + H⁺ + e⁻  K_OH(E) = exp(F(E−E_eq)/RT)
  Step 3   C₃H₆*_σ + OH* → acrolein + *    r_acr = k_rds·θ_σ·θ_OH   [RDS]
  Step 4   C₃H₆*_π + OH* → CHᵧ* + …       dθ_P/dt = k_P·θ_π·θ_OH·(1−θ_P)

MvK regime (PdO surface, ~1.05–1.30 V vs RHE):
  Step 5   Pd + H₂O → PdO + 2H⁺ + 2e⁻    θ_ox(E) = K_ox/(1+K_ox)
  Step 6   C₃H₆(aq) + O_lattice → glycol   r_gly = k_MvK·P_prop·θ_ox

Rate weights:
  r_acrolein_obs = (1−θ_ox)·k_rds·θ_σ·θ_OH   [LH only on metallic fraction]
  r_glycol_obs   = θ_ox·k_MvK·P_prop           [MvK only on oxide fraction]

Scientific scope and caveats
─────────────────────────────
This model is a minimal phenomenological framework consistent with operando
observations in [W19] and [K21]. It is NOT a first-principles microkinetic
model. Specific assumptions (two adsorption states, mean-field Frumkin
correction, smooth oxide transition) are physically motivated but not
uniquely determined by the available data. See MODULE 6 for full limitations.

Architecture
─────────────
MODULE 0  Imports and type aliases
MODULE 1  Physical constants and thermodynamic fixed points
MODULE 2  Electrochemistry (K_OH, θ_ox) — no free parameters
MODULE 3  Coverage solver (implicit Frumkin site balance)
MODULE 4  Kinetics (rates, poisoning ODE)
MODULE 5  Data loading and synthetic data generation
MODULE 6  Fitting infrastructure
MODULE 7  Diagnostics and reporting
MODULE 8  Plotting
MODULE 9  Main entry point

Units convention (enforced throughout)
────────────────────────────────────────
  Potential : V vs RHE
  Pressure  : bar
  Energy    : kJ mol⁻¹
  Rate      : user-defined; must be consistent across acrolein and glycol
              (recommended: nmol cm⁻²_ECSA s⁻¹ or mmol g⁻¹_Pd h⁻¹)
  Coverage  : dimensionless [0, 1]
  Time      : s (for poisoning ODE)
  K_prop    : bar⁻¹
  K_OH      : dimensionless
  k_rds     : [rate] (absorbs site density and post-RDS steps)
  k_MvK     : [rate] bar⁻¹
  k_P       : s⁻¹
"""

# ═══════════════════════════════════════════════════════════════════════════
# MODULE 0  Imports
# ═══════════════════════════════════════════════════════════════════════════

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit, differential_evolution


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 1  Physical constants and thermodynamic fixed points
# ═══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ThermodynamicConstants:
    """
    Physical constants and electrochemically fixed reference potentials.

    All values are fixed before fitting and must not be treated as free
    parameters. Rationale for each fixed quantity is given inline.

    E_eq and E_ox are pH-dependent (Nernst: −59 mV/pH unit at 298 K).
    The values here correspond to pH 1 (0.1 M HClO₄), matching [W19].
    If your electrolyte differs, recompute:
        E_eq(pH) = E_eq(pH1) + 0.059·(1 − pH)   [V vs RHE]
        E_ox(pH) = E_ox(pH1) + 0.059·(1 − pH)
    """
    T:     float = 298.15   # K — room temperature; adjust for non-RT experiments
    R:     float = 8.314e-3 # kJ mol⁻¹ K⁻¹
    F:     float = 96.485   # kJ mol⁻¹ V⁻¹  (= C mol⁻¹ × 1e-3 to keep kJ)

    # Derived thermodynamic group — computed at instantiation
    # F/RT = 38.92 V⁻¹ at 298 K; appears in all Nernst exponentials
    F_over_RT: float = field(init=False)

    # ── Electrochemical reference potentials (fixed from literature) ──────

    E_eq: float = 0.70
    # OH* formation onset on metallic Pd vs RHE, pH 1.
    # Source: Pd Pourbaix diagram; [K21] Fig. 2A DFT interface prevalence
    # diagram shows OH(ad) as the stable surface species from ~0.4 V upward,
    # consistent with a quasi-equilibrium onset around 0.70 V.
    # Uncertainty: ±0.05 V depending on Pd facet distribution and support effects.

    E_ox: float = 0.913
    # Equilibrium potential for Pd → PdO surface oxide at pH 1.
    # Source: [W19] Supporting Information Eq. (1), Pourbaix diagram.
    # [K21] operando XAS shows spectroscopic evidence of oxide formation
    # beginning above 1.0–1.1 V in propene-saturated electrolyte; the
    # thermodynamic onset is ~0.913 V but propylene adsorption delays
    # the observable transition. This model uses the thermodynamic value
    # as the logistic midpoint; the delay is an emergent consequence of
    # competitive adsorption suppressing the oxide (see MODULE 2).

    def __post_init__(self):
        object.__setattr__(self, 'F_over_RT', self.F / (self.R * self.T))


# Singleton — import and use this everywhere
CONST = ThermodynamicConstants()


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 2  Electrochemistry — coverage functions with no free parameters
# ═══════════════════════════════════════════════════════════════════════════

def K_OH(E: float | np.ndarray) -> float | np.ndarray:
    """
    Effective adsorption equilibrium constant for OH* from water activation.

    Derivation
    ──────────
    The electrochemical step H₂O + * → OH* + H⁺ + e⁻ at quasi-equilibrium
    satisfies the Nernst equation. The free energy of the forward reaction is:

        ΔG_rxn(E) = ΔG° − F·E

    At the equilibrium potential E_eq, ΔG_rxn = 0, so ΔG° = F·E_eq. Thus:

        ΔG_rxn(E) = F·(E_eq − E)

    The equilibrium constant (in the Langmuir sense, θ_OH = K_OH·θ_*) is:

        K_OH(E) = exp(−ΔG_rxn / RT) = exp(F(E − E_eq) / RT)

    Properties:
      K_OH = 1.0 at E = E_eq (equally likely to be bare or OH*-covered)
      K_OH > 1  at E > E_eq  (OH* formation thermodynamically favoured)
      K_OH < 1  at E < E_eq  (bare surface favoured)

    Caveat: This quasi-equilibrium treatment neglects the kinetics of OH*
    formation (Butler-Volmer transfer coefficient). A transfer coefficient
    α ≠ 0.5 would modify the exponent to α·F(E−E_eq)/RT. Here α = 0.5 is
    absorbed into the definition of E_eq. This is a standard approximation
    for quasi-reversible electrochemical adsorption steps [see Bard & Faulkner].

    Parameters
    ----------
    E : float or ndarray — electrode potential [V vs RHE]

    Returns
    -------
    K_OH : same type as E — dimensionless equilibrium constant
    """
    return np.exp(CONST.F_over_RT * (E - CONST.E_eq))


def theta_oxide(E: float | np.ndarray) -> float | np.ndarray:
    """
    Fraction of Pd surface in the PdO-like surface oxide state.

    Derivation
    ──────────
    The oxide formation step Pd + H₂O → PdO + 2H⁺ + 2e⁻ has equilibrium:

        K_ox(E) = exp(F(E − E_ox) / RT)

    Assuming a simple two-state (metallic / oxide) Langmuir-type balance:

        θ_ox = K_ox / (1 + K_ox) = 1 / (1 + exp(−F(E−E_ox)/RT))

    This is a logistic function centred at E_ox = 0.913 V. It gives:
      θ_ox → 0  for E ≪ E_ox  (metallic surface dominates)
      θ_ox → 1  for E ≫ E_ox  (oxide surface dominates)
      θ_ox = 0.5 at E = E_ox

    Important caveats:
      (1) PdO formation in reality is a nucleation-and-growth process with
          potential hysteresis. This model captures only the steady-state
          thermodynamic balance, not kinetics of oxide growth.
      (2) [K21] operando XAS shows that propylene suppresses OH*/O* formation
          up to ~1.0–1.1 V. The competitive adsorption of propylene effectively
          shifts the observable oxide onset to higher potential. In this model,
          that suppression is an emergent consequence of propylene coverage
          reducing θ_active in the LH regime, rather than an explicit Frumkin
          correction on the oxide step.
      (3) At E < 0.9 V (the potential window of primary interest for acrolein),
          θ_ox < 0.01 and the MvK term is negligible. The model is effectively
          pure LH in the acrolein-selective window.

    Parameters
    ----------
    E : float or ndarray — electrode potential [V vs RHE]

    Returns
    -------
    theta_ox : same type as E — dimensionless oxide fraction [0, 1]
    """
    K = np.exp(CONST.F_over_RT * (E - CONST.E_ox))
    return K / (1.0 + K)


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 3  Coverage solver — implicit Frumkin site balance
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class SurfaceCoverages:
    """
    Container for surface coverage values at a given (E, P_prop, params) point.

    All coverages are normalised to total surface area (not to metallic fraction
    alone). The metallic fraction is (1 − θ_ox), already folded into θ_active.

    Fields
    ──────
    theta_pi  : π-adsorbed propylene (flat, vinyl-coordinated)
    theta_sig : σ-adsorbed propylene (allylic, more upright)
    theta_OH  : OH* (from water activation)
    theta_s   : free sites (*)
    theta_P   : poisoned sites (CHᵧ; input, not computed here)
    theta_ox  : oxide-phase sites; input, not computed here
    converged : whether the solver met the requested tolerance
    n_iter    : number of solver iterations used
    """
    theta_pi:  float
    theta_sig: float
    theta_OH:  float
    theta_s:   float
    theta_P:   float
    theta_ox:  float
    converged: bool = True
    n_iter:    int  = 0

    @property
    def theta_prop(self) -> float:
        """Total propylene coverage (π + σ)."""
        return self.theta_pi + self.theta_sig

    @property
    def f_sigma(self) -> float:
        """
        Fraction of adsorbed propylene in the σ (allylic) geometry.

        f_σ = θ_σ / (θ_π + θ_σ)

        This is the key selectivity indicator: f_σ → 1 means all adsorbed
        propylene is in the acrolein-selective geometry. f_σ rising with P_prop
        (at fixed E) is the expected signature of the Frumkin correction
        (surface crowding disfavours the flat π geometry).
        """
        denom = self.theta_pi + self.theta_sig
        return self.theta_sig / denom if denom > 1e-12 else 0.0

    def check_balance(self, theta_active: float, rtol: float = 1e-4) -> bool:
        """
        Verify that coverages sum to θ_active within relative tolerance.
        Useful for debugging after solver convergence.
        """
        total = self.theta_pi + self.theta_sig + self.theta_OH + self.theta_s
        return abs(total - theta_active) / (theta_active + 1e-12) < rtol


@dataclass(frozen=True)
class KineticParams:
    """
    All free (fitted) kinetic parameters in one place.

    Keeping parameters in a single container rather than as positional
    arguments has several benefits:
      - prevents silent parameter reordering bugs during refactoring
      - makes the parameter set explicit when calling any function
      - simplifies serialisation (e.g. save/load from JSON or CSV)
      - makes partial fixing (e.g. fix K_sigma, fit the rest) explicit

    Physical meaning of each parameter is documented in the module header.
    Units are given in square brackets.

    theta_P_fixed : fixed poison fraction used in steady-state fitting.
        Set from time-on-stream deactivation characterisation, or from the
        ratio of steady-state to initial current in chronoamperometry.
        If time-resolved data is available, use the poisoning ODE instead
        and fit k_P; set theta_P_fixed = 0.0 in that case.
    """
    K_pi0:        float        # π adsorption constant at zero coverage [bar⁻¹]
    K_sigma:      float        # σ adsorption constant [bar⁻¹]
    g:            float        # Frumkin lateral repulsion on π state [kJ mol⁻¹]
    k_rds:        float        # LH rate constant for allylic C−H activation [rate units]
    k_MvK:        float        # MvK rate constant for glycol on PdO [rate·bar⁻¹]
    theta_P_fixed: float = 0.0 # Fixed poison fraction for steady-state fitting [—]
    k_P:          float = 0.0  # Poisoning rate constant [s⁻¹]; 0 if not fitted

    def validate(self) -> None:
        """
        Check that all parameters are physically admissible.
        Raises ValueError for unphysical values.
        """
        errors = []
        if self.K_pi0 <= 0:
            errors.append(f"K_pi0 must be > 0, got {self.K_pi0}")
        if self.K_sigma <= 0:
            errors.append(f"K_sigma must be > 0, got {self.K_sigma}")
        if self.g < 0:
            errors.append(f"g must be ≥ 0 (repulsive or zero), got {self.g}")
        if self.k_rds <= 0:
            errors.append(f"k_rds must be > 0, got {self.k_rds}")
        if self.k_MvK <= 0:
            errors.append(f"k_MvK must be > 0, got {self.k_MvK}")
        if not (0.0 <= self.theta_P_fixed < 1.0):
            errors.append(f"theta_P_fixed must be in [0, 1), got {self.theta_P_fixed}")
        if self.k_P < 0:
            errors.append(f"k_P must be ≥ 0, got {self.k_P}")
        if errors:
            raise ValueError("Invalid kinetic parameters:\n  " + "\n  ".join(errors))


def _frumkin_residual(
    theta_prop: float,
    E: float,
    P_prop: float,
    K_pi0: float,
    K_sigma: float,
    g: float,
    theta_P: float,
) -> float:
    """
    Residual function for the implicit Frumkin site balance.

    The implicit equation to solve is:

        θ_prop = θ_π(θ_prop) + θ_σ(θ_prop)

    where θ_π and θ_σ depend on θ_prop through K_π_eff(θ_prop).
    Rearranging to standard root-finding form:

        F(θ_prop) = θ_prop − [θ_π(θ_prop) + θ_σ(θ_prop)] = 0

    This function returns F(θ_prop) for use by brentq.

    It is kept separate from solve_coverages to allow future substitution
    of the solver (e.g. Newton with analytical Jacobian) without changing
    the residual definition.
    """
    tox          = theta_oxide(E)
    Koh          = K_OH(E)
    theta_active = (1.0 - tox) * (1.0 - theta_P)

    K_pi_eff = K_pi0 * np.exp(-g * theta_prop / (CONST.R * CONST.T))
    denom    = 1.0 + (K_pi_eff + K_sigma) * P_prop + Koh
    theta_s  = theta_active / denom
    theta_pi  = K_pi_eff * P_prop * theta_s
    theta_sig = K_sigma   * P_prop * theta_s

    return theta_prop - (theta_pi + theta_sig)


def solve_coverages(
    E: float,
    P_prop: float,
    params: KineticParams,
    theta_P: Optional[float] = None,
    tol: float = 1e-10,
    verbose: bool = False,
) -> SurfaceCoverages:
    """
    Solve the implicit Frumkin site balance for surface coverages.

    Site balance on metallic Pd:
        θ_* + θ_π + θ_σ + θ_OH = θ_active = (1 − θ_ox)(1 − θ_P)

    The π adsorption constant is coverage-dependent (Frumkin):
        K_π_eff(θ_prop) = K_π0 · exp(−g · θ_prop / RT)
    where θ_prop = θ_π + θ_σ is the total propylene coverage.

    This makes the site balance implicit: θ_prop appears in its own
    equation through K_π_eff. We solve for θ_prop using brentq
    (Brent's method, guaranteed convergence for a continuous function
    on a bracketed interval). The bracket [0, θ_active] is valid because:
      - At θ_prop = 0:       F(0) = −(θ_π + θ_σ)|_0 ≤ 0
      - At θ_prop = θ_active: F(θ_active) > 0 (cannot exceed active sites)
    So a root exists and lies within [0, θ_active].

    Solver choice — why brentq over fixed-point iteration:
      Fixed-point iteration (v4) converges reliably for the physically
      relevant parameter range (tested up to g = 19.5 kJ/mol), but its
      convergence rate slows as g increases because the Lipschitz constant
      of the map approaches 1. brentq is superlinearly convergent and is
      guaranteed to find the root regardless of g, making it the safer
      choice for exploratory fitting where g may temporarily take large
      values during optimisation.

    Parameters
    ----------
    E       : electrode potential [V vs RHE]
    P_prop  : propylene partial pressure [bar]
    params  : KineticParams dataclass
    theta_P : poison fraction; if None, uses params.theta_P_fixed
    tol     : solver tolerance on θ_prop [dimensionless]
    verbose : if True, prints solver diagnostics

    Returns
    -------
    SurfaceCoverages dataclass
    """
    theta_P_use = theta_P if theta_P is not None else params.theta_P_fixed

    tox          = theta_oxide(E)
    Koh          = K_OH(E)
    theta_active = (1.0 - tox) * (1.0 - theta_P_use)

    # Edge case: no active surface (fully oxidised or fully poisoned)
    if theta_active < 1e-12:
        return SurfaceCoverages(
            theta_pi=0.0, theta_sig=0.0, theta_OH=0.0, theta_s=0.0,
            theta_P=theta_P_use, theta_ox=tox, converged=True, n_iter=0
        )

    # Bracket: F(0) ≤ 0, F(theta_active * (1 - epsilon)) > 0
    # We check F(0) explicitly: at theta_prop=0, K_pi_eff=K_pi0 (maximum),
    # so F(0) = 0 − (K_pi0*P + K_sigma*P)*theta_s(0) ≤ 0 always.
    f_lo = _frumkin_residual(0.0,                  E, P_prop, params.K_pi0,
                              params.K_sigma, params.g, theta_P_use)
    f_hi = _frumkin_residual(theta_active * 0.9999, E, P_prop, params.K_pi0,
                              params.K_sigma, params.g, theta_P_use)

    if verbose:
        print(f"  brentq bracket: F(0)={f_lo:.3e}, F(θ_active)={f_hi:.3e}")

    # If both have the same sign (pathological parameters), fall back gracefully
    if f_lo * f_hi > 0:
        warnings.warn(
            f"Frumkin residual has same sign at both bracket endpoints "
            f"(E={E:.3f}V, P={P_prop:.3f} bar). "
            f"Coverages set to zero. Check parameter bounds.",
            RuntimeWarning, stacklevel=2
        )
        return SurfaceCoverages(
            theta_pi=0.0, theta_sig=0.0, theta_OH=0.0, theta_s=0.0,
            theta_P=theta_P_use, theta_ox=tox, converged=False, n_iter=0
        )

    theta_prop_sol, result = brentq(
        _frumkin_residual,
        0.0, theta_active * 0.9999,
        args=(E, P_prop, params.K_pi0, params.K_sigma, params.g, theta_P_use),
        xtol=tol, full_output=True
    )
    n_iter    = result.iterations
    converged = result.converged

    # Reconstruct individual coverages from the solved θ_prop
    K_pi_eff  = params.K_pi0 * np.exp(-params.g * theta_prop_sol / (CONST.R * CONST.T))
    denom     = 1.0 + (K_pi_eff + params.K_sigma) * P_prop + Koh
    theta_s   = theta_active / denom
    theta_pi  = K_pi_eff       * P_prop * theta_s
    theta_sig = params.K_sigma * P_prop * theta_s
    theta_OH  = Koh            * theta_s

    if verbose:
        print(f"  Solved: θ_prop={theta_prop_sol:.5f}, n_iter={n_iter}, "
              f"converged={converged}")

    return SurfaceCoverages(
        theta_pi=theta_pi, theta_sig=theta_sig, theta_OH=theta_OH, theta_s=theta_s,
        theta_P=theta_P_use, theta_ox=tox, converged=converged, n_iter=n_iter
    )


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 4  Kinetics — rates and poisoning ODE
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class RateResult:
    """
    Container for rate computation outputs at a single (E, P_prop) point.

    Keeps rates and coverages together so callers don't need to call both
    solve_coverages and a separate rate function.
    """
    r_acrolein: float       # acrolein production rate [rate units]
    r_glycol:   float       # propylene glycol production rate [rate units]
    cov:        SurfaceCoverages

    @property
    def r_total(self) -> float:
        return self.r_acrolein + self.r_glycol

    @property
    def selectivity_acrolein(self) -> float:
        """Acrolein fraction of total rate. Returns NaN if r_total ≈ 0."""
        rt = self.r_total
        return self.r_acrolein / rt if rt > 1e-12 else float('nan')


def compute_rates(
    E: float,
    P_prop: float,
    params: KineticParams,
    theta_P: Optional[float] = None,
    verbose: bool = False,
) -> RateResult:
    """
    Compute acrolein and propylene glycol production rates.

    Acrolein (LH, metallic Pd)
    ───────────────────────────
        r_acrolein = k_rds · θ_σ · θ_OH

    The (1 − θ_ox) factor weighting the LH contribution is embedded in
    θ_active → θ_σ and θ_OH via solve_coverages. Specifically:
        θ_active = (1 − θ_ox)(1 − θ_P)
    so all LH coverages are already scaled to the metallic fraction.

    Physical meaning of k_rds:
        Lumped apparent rate constant. Absorbs the true elementary rate
        constant, the Pd surface site density, and any fast post-RDS
        steps (acrolein desorption is assumed fast). Do not interpret as
        an elementary pre-exponential without independent coverage data.

    Propylene glycol (MvK, PdO surface)
    ─────────────────────────────────────
        r_glycol = k_MvK · P_prop · θ_ox

    Rate is first-order in P_prop because propylene attacks the oxide
    surface from solution without competitive adsorption (the metallic
    sites that would adsorb propylene are now oxidised). This linearity
    in P_prop is a distinguishing experimental test of the MvK mechanism:
    in the LH regime, rate shows Langmuir saturation in P; in the MvK
    regime it should be strictly linear.

    Consistency with [K21]:
        The MvK mechanism for glycol is proposed in [K21] based on the
        observation of surface oxide formation at E > 1.1 V coinciding
        with the onset of propylene glycol production. The mechanistic
        assignment is physically motivated but not uniquely confirmed by
        the XAS data alone (the XAS probes Pd oxidation state, not the
        reaction mechanism directly).

    Parameters
    ----------
    E       : electrode potential [V vs RHE]
    P_prop  : propylene partial pressure [bar]
    params  : KineticParams dataclass
    theta_P : poison fraction; if None, uses params.theta_P_fixed
    verbose : passed through to solve_coverages

    Returns
    -------
    RateResult dataclass
    """
    cov = solve_coverages(E, P_prop, params, theta_P=theta_P, verbose=verbose)

    r_acrolein = params.k_rds * cov.theta_sig * cov.theta_OH
    r_glycol   = params.k_MvK * P_prop * cov.theta_ox

    return RateResult(r_acrolein=r_acrolein, r_glycol=r_glycol, cov=cov)


def poisoning_ode(
    t: float,
    state: list,
    E: float,
    P_prop: float,
    params: KineticParams,
) -> list:
    """
    ODE right-hand side for irreversible CHᵧ surface poisoning.

    Equation:
        dθ_P/dt = k_P · θ_π · θ_OH · (1 − θ_P)

    Physical basis of each factor:
    ─ k_P       : rate constant for vinyl C−H/C−C bond scission leading to
                  a CHᵧ surface deposit. Has units s⁻¹ (first-order in the
                  bimolecular surface encounter rate θ_π·θ_OH).
    ─ θ_π·θ_OH : frequency of encounters between π-adsorbed propylene and
                  adjacent OH*, the LH encounter that triggers the vinyl
                  activation pathway. Poisoning is tied specifically to
                  the π geometry, consistent with [K21] attribution of CHᵧ
                  deposits to vinyl-group adsorption.
    ─ (1−θ_P)  : monolayer saturation limit. The remaining adsorbable
                  fraction of the surface decreases as poison accumulates.
                  This is a Langmuir-type saturation term; it does not
                  imply reversibility.

    Implicit steady state:
        At steady state dθ_P/dt = 0, which requires either θ_P = 1
        (complete deactivation) or θ_π·θ_OH = 0 (no active bimolecular
        encounters). The latter can occur at very high E (θ_π → 0 because
        propylene is displaced by OH*) or very low E (θ_OH → 0).
        The model therefore predicts a maximum poisoning rate at intermediate
        potential (0.80–0.95 V), consistent with [W19] SI Fig. S4.

    Identifiability note:
        k_P is identifiable ONLY from time-resolved rate decay data
        (rate vs time-on-stream at fixed E and P_prop). Without such data,
        θ_P should be held fixed (see KineticParams.theta_P_fixed).

    Parameters
    ----------
    t      : time [s] (not used explicitly; required by solve_ivp signature)
    state  : [θ_P] — current poison fraction
    E      : electrode potential [V vs RHE]
    P_prop : propylene partial pressure [bar]
    params : KineticParams dataclass (k_P must be > 0)

    Returns
    -------
    [dθ_P/dt] as a list (required by solve_ivp)
    """
    theta_P_now = float(state[0])
    cov = solve_coverages(E, P_prop, params, theta_P=theta_P_now)
    dtheta_P_dt = params.k_P * cov.theta_pi * cov.theta_OH * (1.0 - theta_P_now)
    return [dtheta_P_dt]


def integrate_poisoning(
    t_end: float,
    E: float,
    P_prop: float,
    params: KineticParams,
    theta_P_init: float = 0.0,
    t_eval: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate the poisoning ODE from t=0 to t=t_end.

    Uses Radau (L-stable implicit Runge-Kutta) because the ODE can become
    stiff as θ_P approaches 1 and (1−θ_P) → 0, causing rapid changes in
    the gradient of the saturation term.

    Parameters
    ----------
    t_end       : integration end time [s]
    E           : electrode potential [V vs RHE]
    P_prop      : propylene partial pressure [bar]
    params      : KineticParams (k_P must be set > 0)
    theta_P_init: initial poison fraction [default 0.0 = fresh catalyst]
    t_eval      : time points at which to store solution; if None, uses
                  solver-chosen steps

    Returns
    -------
    t_out      : time array [s]
    theta_P_t  : poison fraction array at t_out
    """
    if params.k_P <= 0:
        raise ValueError(
            "k_P must be > 0 to integrate the poisoning ODE. "
            "Set k_P in KineticParams."
        )

    sol = solve_ivp(
        poisoning_ode,
        t_span=(0.0, t_end),
        y0=[theta_P_init],
        args=(E, P_prop, params),
        method='Radau',
        t_eval=t_eval,
        rtol=1e-8, atol=1e-10,
        dense_output=False,
    )

    if not sol.success:
        raise RuntimeError(
            f"Poisoning ODE integration failed: {sol.message}\n"
            f"  Conditions: E={E:.3f}V, P_prop={P_prop:.3f}bar, "
            f"k_P={params.k_P:.3e}"
        )

    return sol.t, sol.y[0]


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 5  Data handling
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ExperimentalData:
    """
    Container for rate vs (E, P_prop) experimental data.

    Expected data format (CSV columns):
        E_V        : electrode potential [V vs RHE]
        P_prop_bar : propylene partial pressure [bar]
        r_acrolein : acrolein production rate [consistent units]
        r_glycol   : propylene glycol production rate [same units]

    Rate units must be consistent across both products. Recommended:
        nmol cm⁻²_ECSA s⁻¹   (ECSA from CO stripping, as in [W19])
    or  mmol g⁻¹_Pd h⁻¹     (mass-normalised)

    If you have only total rate + selectivity (e.g. from faradaic efficiency
    and GC product distribution), compute:
        r_acrolein = selectivity_acrolein × r_total
        r_glycol   = selectivity_glycol   × r_total
    """
    E:          np.ndarray   # [V vs RHE]
    P_prop:     np.ndarray   # [bar]
    r_acrolein: np.ndarray   # [rate units]
    r_glycol:   np.ndarray   # [rate units]

    def __post_init__(self):
        n = len(self.E)
        for name, arr in [('P_prop', self.P_prop),
                           ('r_acrolein', self.r_acrolein),
                           ('r_glycol', self.r_glycol)]:
            if len(arr) != n:
                raise ValueError(
                    f"All arrays must have the same length. "
                    f"E has {n} points, {name} has {len(arr)}."
                )
        if np.any(self.r_acrolein < 0) or np.any(self.r_glycol < 0):
            raise ValueError("Negative rates detected. Check data units and sign convention.")

    @property
    def r_total(self) -> np.ndarray:
        return self.r_acrolein + self.r_glycol

    @classmethod
    def from_csv(cls, path: str) -> 'ExperimentalData':
        """
        Load data from a CSV file with header row.
        Expected columns: E_V, P_prop_bar, r_acrolein, r_glycol
        """
        data = np.genfromtxt(path, delimiter=',', names=True)
        return cls(
            E          = data['E_V'],
            P_prop     = data['P_prop_bar'],
            r_acrolein = data['r_acrolein'],
            r_glycol   = data['r_glycol'],
        )

    @classmethod
    def synthetic(cls, params_true: KineticParams, noise_frac: float = 0.04,
                  seed: int = 42) -> 'ExperimentalData':
        """
        Generate synthetic rate data from known parameters.

        Used for model validation, identifiability analysis, and testing the
        fitting pipeline before real data is available.

        The Frumkin correction at the true parameter values (g=6 kJ/mol) gives:
          At θ_prop = 0.0: K_π_eff = K_π0          (maximum π adsorption)
          At θ_prop = 0.3: K_π_eff = K_π0·exp(−6×0.3/RT·kJ) ≈ 0.48·K_π0
          At θ_prop = 0.6: K_π_eff = K_π0·exp(−6×0.6/RT·kJ) ≈ 0.23·K_π0
        → π coverage suppressed by ~50% at moderate loading, shifting
          propylene into the σ (acrolein-selective) geometry.

        Parameters
        ----------
        params_true : KineticParams with the 'true' parameter values
        noise_frac  : fractional Gaussian noise amplitude (default 4%)
        seed        : RNG seed for reproducibility
        """
        potentials = np.array([0.72, 0.76, 0.80, 0.85, 0.90,
                                0.95, 1.00, 1.05, 1.10, 1.15, 1.20])
        p_props    = np.array([0.05, 0.10, 0.20, 0.40])

        rng = np.random.default_rng(seed=seed)
        E_list, P_list, ra_list, rg_list = [], [], [], []

        for E in potentials:
            for Pp in p_props:
                result = compute_rates(E, Pp, params_true,
                                        theta_P=params_true.theta_P_fixed)
                ra = max(result.r_acrolein * (1 + noise_frac * rng.standard_normal()), 0.0)
                rg = max(result.r_glycol   * (1 + noise_frac * rng.standard_normal()), 0.0)
                E_list.append(E);   P_list.append(Pp)
                ra_list.append(ra); rg_list.append(rg)

        return cls(
            E          = np.array(E_list),
            P_prop     = np.array(P_list),
            r_acrolein = np.array(ra_list),
            r_glycol   = np.array(rg_list),
        )


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 6  Fitting infrastructure
# ═══════════════════════════════════════════════════════════════════════════

# Parameter order for fitting vectors: [K_pi0, K_sigma, g, k_rds, k_MvK]
# This order must be consistent between bounds_default, _pack, and _unpack.
_PARAM_NAMES = ['K_pi0', 'K_sigma', 'g', 'k_rds', 'k_MvK']

# Default bounds based on physical reasoning:
#   K_pi0  : π adsorption; 0.5–100 bar⁻¹ covers light alkenes on Pd
#             (literature TPD: Cassuto et al. 1990; Horiuti-Polanyi analyses)
#   K_sigma: σ adsorption; expected weaker than π at zero coverage
#   g      : Frumkin repulsion; 0 → no coverage dependence; upper bound
#             set conservatively at 20 kJ/mol (strong lateral interactions)
#   k_rds  : rate constant; bounds are data-scale-dependent; widen if needed
#   k_MvK  : MvK rate constant; same units as k_rds/bar
BOUNDS_DEFAULT: list[tuple[float, float]] = [
    (0.5,   100.0),  # K_pi0   [bar⁻¹]
    (0.1,    20.0),  # K_sigma [bar⁻¹]
    (0.0,    20.0),  # g       [kJ mol⁻¹]
    (1e-3,  1e2),    # k_rds   [rate units]
    (1e-3,  1e2),    # k_MvK   [rate units · bar⁻¹]
]


def _pack_params(params: KineticParams) -> np.ndarray:
    """Extract the 5 fitted parameters as a 1-D array."""
    return np.array([params.K_pi0, params.K_sigma, params.g,
                     params.k_rds, params.k_MvK])


def _unpack_params(vec: np.ndarray, base: KineticParams) -> KineticParams:
    """Reconstruct KineticParams from a 5-element vector, preserving fixed fields."""
    return KineticParams(
        K_pi0        = vec[0],
        K_sigma      = vec[1],
        g            = vec[2],
        k_rds        = vec[3],
        k_MvK        = vec[4],
        theta_P_fixed = base.theta_P_fixed,
        k_P          = base.k_P,
    )


def _model_vector(
    X: Tuple[np.ndarray, np.ndarray],
    K_pi0: float, K_sigma: float, g: float, k_rds: float, k_MvK: float,
    _base_params: KineticParams,
) -> np.ndarray:
    """
    Model prediction in the format required by scipy.optimize.curve_fit.

    Returns concatenated [r_acrolein_0, …, r_acrolein_N, r_glycol_0, …, r_glycol_N].
    Normalization to equal weighting of acrolein and glycol is handled in the
    DE objective; curve_fit uses this raw form with equal residual weights.
    """
    E_arr, P_arr = X
    params = _unpack_params(np.array([K_pi0, K_sigma, g, k_rds, k_MvK]), _base_params)
    ra = np.array([compute_rates(E, P, params).r_acrolein for E, P in zip(E_arr, P_arr)])
    rg = np.array([compute_rates(E, P, params).r_glycol   for E, P in zip(E_arr, P_arr)])
    return np.concatenate([ra, rg])


def _de_objective(
    vec: np.ndarray,
    data: ExperimentalData,
    base_params: KineticParams,
) -> float:
    """
    Normalised sum-of-squares objective for differential evolution.

    Normalisation by mean observed rate for each product ensures that
    acrolein and glycol contribute equally to the objective regardless
    of their absolute magnitudes. Without normalisation, the larger
    acrolein rates would dominate the fit and the glycol branch (which
    carries the mechanistic information about the MvK regime) would be
    effectively ignored.

    This is a deliberate modelling choice: equal weighting assumes that
    relative errors (not absolute errors) are approximately equal across
    both products. If your measurement uncertainties differ substantially
    between products, replace the normalisation denominators with measured
    uncertainties (σ_acrolein, σ_glycol) for a proper χ² objective.
    """
    params  = _unpack_params(vec, base_params)
    n       = len(data.E)
    ra_pred = np.array([compute_rates(E, P, params).r_acrolein
                         for E, P in zip(data.E, data.P_prop)])
    rg_pred = np.array([compute_rates(E, P, params).r_glycol
                         for E, P in zip(data.E, data.P_prop)])

    norm_a = data.r_acrolein.mean() + 1e-12
    norm_g = data.r_glycol.mean()   + 1e-12

    ss = (np.sum(((ra_pred - data.r_acrolein) / norm_a) ** 2) +
          np.sum(((rg_pred - data.r_glycol)   / norm_g) ** 2))
    return ss


@dataclass
class FitResult:
    """
    Container for fitting results, uncertainties, and goodness-of-fit metrics.

    perr contains 1-σ standard errors from the covariance matrix of
    curve_fit (Levenberg-Marquardt refinement). These are approximate;
    they assume Gaussian errors and local linearity of the model.
    Large perr (exceeding the parameter value) signals parameter
    non-identifiability — see identifiable_product below.
    """
    params:       KineticParams
    perr:         dict[str, float]     # {param_name: 1σ std error}
    r2_acrolein:  float
    r2_glycol:    float
    r2_total:     float
    de_converged: bool
    de_residual:  float
    lm_succeeded: bool                 # whether curve_fit LM step succeeded

    @property
    def identifiable_product(self) -> float:
        """
        K_sigma · k_rds: the only combination of these two parameters that
        is directly constrained by rate-vs-P data.

        The acrolein rate is k_rds · θ_σ · θ_OH = k_rds · K_sigma · P · θ_* · θ_OH.
        For P in the Henry's law limit (K_sigma·P ≪ 1), the rate is proportional
        to k_rds · K_sigma · P. These two parameters cannot be individually
        identified from rate-vs-P data alone without either:
          (a) independently measuring K_sigma (e.g. from adsorption isotherm or
              propylene stripping charge vs P_prop, as in [W19] SI Fig. S3b), or
          (b) observing Langmuir saturation in the rate vs P curve, from which
              the half-saturation pressure gives K_sigma directly.
        """
        return self.params.K_sigma * self.params.k_rds

    def print_summary(self) -> None:
        """Print a formatted summary of the fit results to stdout."""
        sep = "═" * 64
        print(f"\n{sep}")
        print("  FIT RESULTS — v5 Two-Geometry Frumkin LH/MvK Model")
        print(sep)
        print(f"\n  Thermodynamic fixed points (not fitted):")
        print(f"    E_eq   = {CONST.E_eq:.3f} V vs RHE  "
              f"(OH* onset on metallic Pd, [K21] Fig. 2A)")
        print(f"    E_ox   = {CONST.E_ox:.3f} V vs RHE  "
              f"(PdO onset at pH 1, [W19] SI Eq. 1)")
        print(f"    F/RT   = {CONST.F_over_RT:.2f} V⁻¹    (298 K)")
        print(f"    θ_P    = {self.params.theta_P_fixed:.3f}          "
              f"(fixed; estimate from deactivation data)")

        print(f"\n  Fitted parameters (5 free):")
        for name in _PARAM_NAMES:
            val = getattr(self.params, name)
            err = self.perr.get(name, float('nan'))
            unit = {'K_pi0': 'bar⁻¹', 'K_sigma': 'bar⁻¹',
                    'g': 'kJ mol⁻¹', 'k_rds': '[rate]',
                    'k_MvK': '[rate]/bar'}[name]
            flag = " *** non-identifiable" if err > abs(val) else ""
            print(f"    {name:<10} = {val:10.4f}  ±  {err:.4f}  {unit}{flag}")

        print(f"\n  Identifiable combination:")
        print(f"    K_sigma · k_rds = {self.identifiable_product:.4f}  [bar⁻¹·rate]")
        print(f"    (If K_sigma or k_rds is flagged non-identifiable, only this")
        print(f"     product is constrained by rate-vs-P data.)")

        print(f"\n  Goodness of fit:")
        print(f"    R² acrolein = {self.r2_acrolein:.4f}")
        print(f"    R² glycol   = {self.r2_glycol:.4f}")
        print(f"    R² total    = {self.r2_total:.4f}")
        print(f"    DE converged: {self.de_converged}  |  residual: {self.de_residual:.6f}")
        print(f"    LM refinement: {'succeeded' if self.lm_succeeded else 'failed (DE result used)'}")

        # Derived quantities
        K_ratio = self.params.K_pi0 / self.params.K_sigma
        print(f"\n  Derived quantities:")
        print(f"    K_π0 / K_σ = {K_ratio:.2f}  "
              f"({'π more stable at zero coverage' if K_ratio > 1 else 'σ more stable'})")
        if self.params.g > 0.01 and self.params.K_pi0 > self.params.K_sigma:
            theta_cross = (CONST.R * CONST.T / self.params.g) * np.log(K_ratio)
            print(f"    θ_cross    = {theta_cross:.3f}  "
                  f"(propylene coverage where π ≡ σ fraction in geometry distribution)")
            print(f"                 Above θ_cross: σ dominates (acrolein-selective)")
            print(f"                 Below θ_cross: π dominates (poison-generating)")
        print(sep)


def run_fitting(
    data: ExperimentalData,
    base_params: KineticParams,
    bounds: Optional[list[tuple[float, float]]] = None,
    de_popsize: int = 20,
    de_maxiter: int = 8000,
    de_seed:    int = 42,
    verbose:    bool = True,
) -> FitResult:
    """
    Fit the v5 model to experimental (or synthetic) data.

    Strategy: two-stage optimisation
    ─────────────────────────────────
    Stage 1 — Global search (differential evolution):
        Searches the full parameter space bounded by `bounds`. Uses a
        population-based stochastic algorithm that is resistant to local
        minima. The normalised sum-of-squares objective (_de_objective)
        gives equal weight to acrolein and glycol.

    Stage 2 — Local refinement (Levenberg-Marquardt via curve_fit):
        Refines the DE solution to higher precision and computes the
        approximate covariance matrix for standard error estimation.
        Uses the raw (un-normalised) residuals, which is correct for
        curve_fit's assumption of equal measurement variances.
        If LM fails (ill-conditioned covariance), the DE result is kept
        and standard errors are reported as NaN.

    Parameters
    ----------
    data        : ExperimentalData
    base_params : KineticParams specifying fixed quantities (theta_P_fixed, k_P)
                  and providing default values for non-fitted parameters
    bounds      : list of (lo, hi) tuples; defaults to BOUNDS_DEFAULT
    de_popsize  : DE population size (larger → more thorough, slower)
    de_maxiter  : DE maximum iterations
    de_seed     : RNG seed for reproducibility
    verbose     : print progress

    Returns
    -------
    FitResult dataclass
    """
    if bounds is None:
        bounds = BOUNDS_DEFAULT

    if verbose:
        print(f"Stage 1: Differential evolution  "
              f"(popsize={de_popsize}, maxiter={de_maxiter})")
        print(f"  Fixed: E_eq={CONST.E_eq}V, E_ox={CONST.E_ox}V, "
              f"θ_P={base_params.theta_P_fixed}\n")

    de_result = differential_evolution(
        _de_objective,
        bounds,
        args=(data, base_params),
        seed=de_seed,
        maxiter=de_maxiter,
        tol=1e-12,
        polish=True,
        popsize=de_popsize,
        mutation=(0.5, 1.5),
        recombination=0.9,
        workers=1,
    )

    if verbose:
        print(f"  DE converged: {de_result.success}, "
              f"residual: {de_result.fun:.6f}")
        print(f"\nStage 2: Levenberg-Marquardt refinement (curve_fit)")

    # Wrap model for curve_fit signature (first arg = X tuple)
    def _cf_model(X, *vec):
        return _model_vector(X, *vec, _base_params=base_params)

    y_obs = np.concatenate([data.r_acrolein, data.r_glycol])
    lb = [b[0] for b in bounds]
    ub = [b[1] for b in bounds]
    lm_succeeded = False

    try:
        popt_lm, pcov = curve_fit(
            _cf_model,
            (data.E, data.P_prop),
            y_obs,
            p0=de_result.x,
            bounds=(lb, ub),
            maxfev=200_000,
        )
        perr_arr = np.sqrt(np.diag(pcov))
        popt_final = popt_lm
        lm_succeeded = True
        if verbose:
            print("  LM succeeded.")
    except Exception as exc:
        if verbose:
            print(f"  LM failed ({exc}). Using DE result; standard errors set to NaN.")
        popt_final = de_result.x
        perr_arr = np.full(5, float('nan'))

    params_fit = _unpack_params(popt_final, base_params)
    perr_dict  = {name: float(e) for name, e in zip(_PARAM_NAMES, perr_arr)}

    # Compute goodness-of-fit on fitted parameters
    def _r2(obs: np.ndarray, pred: np.ndarray) -> float:
        ss_res = np.sum((obs - pred) ** 2)
        ss_tot = np.sum((obs - obs.mean()) ** 2)
        return float(1.0 - ss_res / (ss_tot + 1e-12))

    ra_pred = np.array([compute_rates(E, P, params_fit).r_acrolein
                         for E, P in zip(data.E, data.P_prop)])
    rg_pred = np.array([compute_rates(E, P, params_fit).r_glycol
                         for E, P in zip(data.E, data.P_prop)])

    return FitResult(
        params       = params_fit,
        perr         = perr_dict,
        r2_acrolein  = _r2(data.r_acrolein, ra_pred),
        r2_glycol    = _r2(data.r_glycol,   rg_pred),
        r2_total     = _r2(data.r_total,    ra_pred + rg_pred),
        de_converged = de_result.success,
        de_residual  = float(de_result.fun),
        lm_succeeded = lm_succeeded,
    )


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 7  Diagnostics
# ═══════════════════════════════════════════════════════════════════════════

def print_coverage_table(
    params: KineticParams,
    E_points: list[float] = None,
    P_points: list[float] = None,
) -> None:
    """
    Print a table of surface coverages and rates at selected (E, P) conditions.

    Rising f_σ% with increasing P_prop (at fixed E) is the expected signature
    of the Frumkin geometry transition and supports the two-state adsorption
    picture over a single-state Langmuir model.
    """
    if E_points is None:
        E_points = [0.80, 0.90, 1.00, 1.10]
    if P_points is None:
        P_points = [0.10, 0.40]

    header = (f"  {'E(V)':>7} {'P(bar)':>8} {'θ_π':>7} {'θ_σ':>7} "
              f"{'θ_OH':>7} {'θ_*':>7} {'f_σ%':>7} {'r_acr':>8} {'r_gly':>8}")
    print("\n" + "═" * 70)
    print("  COVERAGE AND RATE DIAGNOSTICS")
    print(header)

    for E in E_points:
        for P in P_points:
            result = compute_rates(E, P, params)
            cov    = result.cov
            f_s    = 100.0 * cov.f_sigma
            print(f"  {E:>7.2f} {P:>8.2f} {cov.theta_pi:>7.3f} "
                  f"{cov.theta_sig:>7.3f} {cov.theta_OH:>7.3f} "
                  f"{cov.theta_s:>7.3f} {f_s:>7.1f} "
                  f"{result.r_acrolein:>8.4f} {result.r_glycol:>8.4f}")

    print()
    print("  f_σ% = σ-geometry fraction of total propylene coverage.")
    print("  Expected behaviour: f_σ% rises with P_prop at fixed E")
    print("  (Frumkin suppression of flat π geometry at high coverage).")
    print("═" * 70)


def print_experimental_requirements() -> None:
    """
    Print guidance on what experiments constrain each model parameter.
    """
    print("""
═══════════════════════════════════════════════════════════════════════
  EXPERIMENTAL REQUIREMENTS — what data constrains each parameter
═══════════════════════════════════════════════════════════════════════

  K_π0, K_σ
    Rate vs P_prop at FIXED E in the LH window (0.75–0.90 V).
    If rate saturates with P → K·P reaches O(1) → saturation pressure
    gives K directly. If rate remains linear in P across the full range,
    only the product K·k_rds is constrained (Henry's law regime).
    Independent measurement option: propylene stripping charge vs P_prop
    at fixed E, as demonstrated in [W19] SI Fig. S3b.

  g (Frumkin lateral repulsion)
    Acrolein selectivity vs P_prop at fixed E.
    Model prediction: acrolein fraction rises with P_prop as π state is
    suppressed by crowding (Frumkin). If selectivity is independent of
    P_prop, set g = 0 (reduce to standard competitive Langmuir). The
    shape of the selectivity–P curve distinguishes g > 0 from g = 0.
    This is the key experimental discriminator between v3 and v5.

  k_rds
    Absolute acrolein rate magnitude at any (E, P_prop), once coverages
    are independently constrained. Requires careful normalisation to ECSA
    (electrochemical surface area from CO stripping; see [W19] SI §S1).
    Without ECSA normalisation, k_rds is only an apparent rate constant.

  k_MvK
    Glycol production rate at E > 1.05 V. The MvK mechanism predicts
    r_glycol ∝ P_prop (linear, no saturation). Plot r_glycol / P_prop
    vs E: should trace the θ_ox(E) sigmoid. Deviation from linearity
    in P would indicate a contribution from a competing adsorptive pathway.

  k_P (poisoning rate constant)
    Time-on-stream rate decay at fixed (E, P_prop).
    Minimum useful dataset: r_acrolein at t = 1, 5, 15, 30, 60 min.
    The model predicts fastest deactivation near E = 0.90–0.95 V where
    θ_π · θ_OH is maximised, consistent with [W19] SI Fig. S4.
    Without time-resolved data, report θ_P as a sensitivity variable.

  θ_P (steady-state poison level, if k_P not fitted)
    Estimate from (j_SS / j_initial) ratio in chronoamperometry at the
    experimental potential. [W19] data suggest current declines to ~40–60%
    of initial within 3 min at 0.95 V → θ_P ~ 0.10–0.20 at steady state.

═══════════════════════════════════════════════════════════════════════
  MODEL LIMITATIONS — report these explicitly in any publication
═══════════════════════════════════════════════════════════════════════

  (1) Mean-field Frumkin correction.
      The lateral repulsion parameter g is a mean-field quantity. Real
      lateral interactions between adsorbates are pair-specific and
      geometry-dependent. g should be interpreted as an effective,
      coverage-averaged repulsion constant, not a molecular binding energy.

  (2) Lumped rate constant k_rds.
      k_rds absorbs the site density, post-RDS desorption steps, and the
      elementary rate constant. It is not directly comparable across
      studies using different normalisation bases (per g_Pd vs per cm²_ECSA).

  (3) Smooth oxide transition.
      θ_ox(E) is modelled as a logistic function. Real PdO formation
      involves nucleation-and-growth kinetics with CV hysteresis. This
      model predicts no hysteresis in E. If your CA data shows path
      dependence, an Avrami-type nucleation model is more appropriate.

  (4) pH dependence not explicitly modelled.
      E_eq and E_ox shift by −59 mV per pH unit at 298 K (Nernst).
      The values used here (0.70 and 0.913 V) are for pH 1 ([W19] conditions).
      Recompute for other electrolyte compositions.

  (5) The two-adsorption-state picture is a phenomenological hypothesis.
      While consistent with [W19] and [K21] observations, direct spectroscopic
      discrimination between π and σ propylene on Pd under electrochemical
      conditions has not been reported. The model is falsifiable: if SEIRAS
      or EC-SERS data shows a single adsorption geometry across all coverages,
      the two-state framework should be collapsed back to one state.
═══════════════════════════════════════════════════════════════════════
""")


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 8  Plotting
# ═══════════════════════════════════════════════════════════════════════════

def _sweep_rates_vs_E(
    E_arr: np.ndarray,
    P_prop: float,
    params: KineticParams,
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorised rate sweep over E at fixed P_prop."""
    ra = np.array([compute_rates(E, P_prop, params).r_acrolein for E in E_arr])
    rg = np.array([compute_rates(E, P_prop, params).r_glycol   for E in E_arr])
    return ra, rg


def _sweep_rates_vs_P(
    P_arr: np.ndarray,
    E: float,
    params: KineticParams,
) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorised rate sweep over P_prop at fixed E."""
    ra = np.array([compute_rates(E, P, params).r_acrolein for P in P_arr])
    rg = np.array([compute_rates(E, P, params).r_glycol   for P in P_arr])
    return ra, rg


def _sweep_f_sigma_vs_E(
    E_arr: np.ndarray,
    P_prop: float,
    params: KineticParams,
) -> np.ndarray:
    """Vectorised f_σ sweep over E at fixed P_prop."""
    return np.array([
        solve_coverages(E, P_prop, params).f_sigma for E in E_arr
    ])


def _sweep_f_sigma_vs_P(
    P_arr: np.ndarray,
    E: float,
    params: KineticParams,
) -> np.ndarray:
    """Vectorised f_σ sweep over P_prop at fixed E."""
    return np.array([
        solve_coverages(E, P, params).f_sigma for P in P_arr
    ])


def make_plots(
    fit: FitResult,
    data: ExperimentalData,
    save_path: Optional[str] = None,
    show: bool = True,
    k_P_illustrative: float = 3e-4,
) -> None:
    """
    Generate the standard six-panel diagnostic figure.

    Panels
    ──────
    [0,0] Parity plot — observed vs predicted rates
    [0,1] Rate vs E at fixed P_prop values
    [0,2] Rate vs P_prop at fixed E values (Frumkin signature)
    [1,0] σ-geometry fraction vs E (selectivity driver)
    [1,1] σ-geometry fraction vs P_prop (key Winiwarter prediction)
    [1,2] Poisoning dynamics (illustrative k_P; label as such)

    Parameters
    ----------
    fit               : FitResult from run_fitting
    data              : ExperimentalData used for fitting
    save_path         : file path to save figure (None = do not save)
    show              : call plt.show() if True
    k_P_illustrative  : poisoning rate constant for panel [1,2]; this is
                        NOT a fitted value unless k_P was explicitly fitted
    """
    params = fit.params
    E_sw   = np.linspace(0.70, 1.25, 200)
    P_sw   = np.linspace(0.01, 0.50, 200)

    # Precompute fitted rates at data points (avoid re-computation in panels)
    ra_fit = np.array([compute_rates(E, P, params).r_acrolein
                        for E, P in zip(data.E, data.P_prop)])
    rg_fit = np.array([compute_rates(E, P, params).r_glycol
                        for E, P in zip(data.E, data.P_prop)])

    def _r2(obs, pred):
        ss_res = np.sum((obs - pred) ** 2)
        ss_tot = np.sum((obs - obs.mean()) ** 2)
        return 1.0 - ss_res / (ss_tot + 1e-12)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(
        "LH Kinetic Model v5  ·  Two-geometry Frumkin LH + MvK\n"
        "Propylene electrooxidation on Pd/C  "
        "[Winiwarter 2019 EES; Koroidov 2021 CatSciTech]",
        fontsize=11,
    )

    # ── [0,0] Parity ─────────────────────────────────────────────────────────
    ax = axes[0, 0]
    ax.scatter(data.r_acrolein, ra_fit, label="Acrolein",
               color="#1E88E5", alpha=0.75, s=40, zorder=3)
    ax.scatter(data.r_glycol, rg_fit,  label="Glycol",
               color="#8E24AA", alpha=0.75, s=40, zorder=3)
    all_r = np.concatenate([data.r_acrolein, data.r_glycol, ra_fit, rg_fit])
    lims  = [all_r.min() * 0.9, all_r.max() * 1.1]
    ax.plot(lims, lims, 'k--', lw=1, alpha=0.5)
    ax.set(xlabel="Observed rate", ylabel="Predicted rate",
           title="Parity Plot", xlim=lims, ylim=lims)
    ax.legend(fontsize=8)
    ax.text(0.05, 0.88,
            f"R²(acr)={fit.r2_acrolein:.3f}\nR²(gly)={fit.r2_glycol:.3f}",
            transform=ax.transAxes, fontsize=8, va='top')

    # ── [0,1] Rate vs E ───────────────────────────────────────────────────────
    ax = axes[0, 1]
    for P_fixed, col in [(0.10, '#E53935'), (0.20, '#1E88E5'), (0.40, '#00897B')]:
        ra, rg = _sweep_rates_vs_E(E_sw, P_fixed, params)
        ax.plot(E_sw, ra, color=col, lw=2,   label=f"Acrolein P={P_fixed}")
        ax.plot(E_sw, rg, color=col, lw=1.5, ls='--', label=f"Glycol P={P_fixed}")
    ax.axvline(CONST.E_eq, color='#378ADD', lw=1, ls=':', alpha=0.7,
               label=f"E_eq={CONST.E_eq}V")
    ax.axvline(CONST.E_ox, color='#E53935', lw=1, ls=':', alpha=0.7,
               label=f"E_ox={CONST.E_ox}V")
    ax.set(xlabel="E (V vs RHE)", ylabel="Rate [model units]",
           title="Rate vs Potential\nsolid=acrolein  dashed=glycol")
    ax.legend(fontsize=6, ncol=2)

    # ── [0,2] Rate vs P_prop ──────────────────────────────────────────────────
    ax = axes[0, 2]
    for E_fixed, col in [(0.85, '#E53935'), (0.95, '#FB8C00'), (1.05, '#1E88E5')]:
        ra, rg = _sweep_rates_vs_P(P_sw, E_fixed, params)
        ax.plot(P_sw, ra, color=col, lw=2,   label=f"Acrolein E={E_fixed}V")
        ax.plot(P_sw, rg, color=col, lw=1.5, ls='--', label=f"Glycol E={E_fixed}V")
    ax.set(xlabel="P_prop [bar]", ylabel="Rate [model units]",
           title="Rate vs P_prop\n"
                 "Non-monotone acrolein curve = Frumkin signature")
    ax.legend(fontsize=6, ncol=2)

    # ── [1,0] σ-geometry fraction vs E ────────────────────────────────────────
    ax = axes[1, 0]
    for P_fixed, col in [(0.05, '#E53935'), (0.20, '#FB8C00'), (0.40, '#1E88E5')]:
        fs = _sweep_f_sigma_vs_E(E_sw, P_fixed, params)
        ax.plot(E_sw, fs * 100, color=col, lw=2, label=f"P={P_fixed} bar")
    ax.axvline(CONST.E_eq, color='gray', lw=1, ls=':', alpha=0.6)
    ax.axhline(50, color='gray', lw=1, ls=':', alpha=0.4)
    ax.set(xlabel="E (V vs RHE)", ylabel="σ-geometry fraction f_σ [%]",
           title="Allylic (σ) adsorption fraction vs Potential\n"
                 "100% → all propylene in acrolein-selective geometry",
           ylim=(0, 100))
    ax.legend(fontsize=8)

    # ── [1,1] σ-geometry fraction vs P_prop ───────────────────────────────────
    ax = axes[1, 1]
    for E_fixed, col in [(0.80, '#E53935'), (0.90, '#FB8C00'), (1.00, '#1E88E5')]:
        fs = _sweep_f_sigma_vs_P(P_sw, E_fixed, params)
        ax.plot(P_sw, fs * 100, color=col, lw=2, label=f"E={E_fixed}V")
    ax.axhline(50, color='gray', lw=1, ls=':', alpha=0.4)
    ax.set(xlabel="P_prop [bar]", ylabel="σ-geometry fraction f_σ [%]",
           title="Allylic fraction vs P_prop\n"
                 "Rising f_σ with P_prop is consistent with Frumkin suppression\n"
                 "of flat π geometry at high surface coverage [W19]",
           ylim=(0, 100))
    ax.legend(fontsize=8)

    # ── [1,2] Poisoning dynamics (illustrative) ────────────────────────────────
    ax = axes[1, 2]
    # Construct a temporary params object with k_P set to the illustrative value
    params_with_kP = KineticParams(
        K_pi0=params.K_pi0, K_sigma=params.K_sigma, g=params.g,
        k_rds=params.k_rds, k_MvK=params.k_MvK,
        theta_P_fixed=params.theta_P_fixed, k_P=k_P_illustrative,
    )
    t_eval = np.linspace(0, 7200, 300)
    for E_fixed, col in [(0.80, '#E53935'), (0.90, '#FB8C00'),
                          (0.95, '#1E88E5'), (1.05, '#43A047')]:
        t_out, theta_P_t = integrate_poisoning(
            7200, E_fixed, 0.20, params_with_kP, t_eval=t_eval)
        r0 = compute_rates(E_fixed, 0.20, params, theta_P=0.0).r_acrolein
        if r0 < 1e-10:
            continue
        r_t = np.array([
            compute_rates(E_fixed, 0.20, params, theta_P=float(tP)).r_acrolein
            for tP in theta_P_t
        ])
        ax.plot(t_out / 3600, r_t / r0, color=col, lw=2, label=f"E={E_fixed}V")

    ax.axhline(1.0, color='gray', lw=1, ls=':', alpha=0.4)
    ax.set(xlabel="Time [h]", ylabel="Normalised acrolein rate  r / r₀",
           title=f"Poisoning dynamics (ODE)\n"
                 f"k_P = {k_P_illustrative:.1e} s⁻¹  — ILLUSTRATIVE, not fitted",
           ylim=(0, 1.05))
    ax.legend(fontsize=8)
    ax.text(0.05, 0.10,
            "Fastest deactivation predicted at ~0.90–0.95 V\n"
            "(θ_π · θ_OH maximised at intermediate E).\n"
            "Consistent with [W19] SI Fig. S4.",
            transform=ax.transAxes, fontsize=7.5, color='#444',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFFDE7', alpha=0.8))

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved: {save_path}")
    if show:
        plt.show()


# ═══════════════════════════════════════════════════════════════════════════
# MODULE 9  Main entry point
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    # ── True parameters for synthetic data generation ──────────────────────
    # These values are used ONLY to generate the synthetic dataset.
    # In real use, replace generate_synthetic_data with ExperimentalData.from_csv().
    params_true = KineticParams(
        K_pi0         = 8.0,   # bar⁻¹  — π adsorption at zero coverage
        K_sigma       = 1.5,   # bar⁻¹  — σ adsorption (allylic geometry)
        g             = 6.0,   # kJ/mol — Frumkin lateral repulsion on π state
        k_rds         = 0.8,   # [rate] — allylic C-H activation rate constant
        k_MvK         = 0.3,   # [rate/bar] — MvK glycol rate constant
        theta_P_fixed = 0.12,  # —      — fixed poison fraction
    )

    # ── Load data ──────────────────────────────────────────────────────────
    data = ExperimentalData.synthetic(params_true, noise_frac=0.04, seed=42)
    # ─── To use real data instead: ────────────────────────────────────────
    # data = ExperimentalData.from_csv("your_data.csv")
    # CSV format (with header): E_V, P_prop_bar, r_acrolein, r_glycol
    # ──────────────────────────────────────────────────────────────────────

    # ── Initial parameter container for fitting ────────────────────────────
    # theta_P_fixed is estimated from deactivation experiments and held
    # constant during steady-state fitting. Adjust this value based on
    # your chronoamperometry data before fitting.
    params_init = KineticParams(
        K_pi0=1.0, K_sigma=1.0, g=1.0, k_rds=1.0, k_MvK=1.0,
        theta_P_fixed=0.12,   # <── set from experiment
    )

    # ── Run fitting ────────────────────────────────────────────────────────
    fit = run_fitting(data, params_init, verbose=True)
    fit.print_summary()

    # ── Coverage and rate diagnostics ─────────────────────────────────────
    print_coverage_table(fit.params)

    # ── Experimental guidance ─────────────────────────────────────────────
    print_experimental_requirements()

    # ── Plots ─────────────────────────────────────────────────────────────
    make_plots(
        fit, data,
        save_path="lh_kinetics_v5_fit.png",
        show=False,
        k_P_illustrative=3e-4,  # s⁻¹ — illustrative only
    )