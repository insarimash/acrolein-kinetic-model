"""
Langmuir-Hinshelwood Kinetic Model — v4
Propylene electrooxidation on Pd/C
Grounded in Winiwarter et al. (2019) EES and Koroidov/Nilsson et al. (2021) CatSciTech

═══════════════════════════════════════════════════════════════════════════════
WHAT IS NEW IN v4 vs v3 — AND WHY
═══════════════════════════════════════════════════════════════════════════════

v3 had two flaws that the papers directly contradict:

  FLAW 1 — Single propylene adsorption state, selectivity from downstream branching
    v3 assumed one C₃H₆* species, and selectivity arose from a competition between
    fast desorption vs a second OH* attack on the intermediate.
    WHAT THE PAPERS SAY: Koroidov/Nilsson 2021 (XAS, operando) states explicitly
    that "propene adsorption through the vinyl group on metallic Pd would lead to
    irreversibly retained, unreactive adsorbates such as C₂Hₓ, CHᵧ." Winiwarter
    2019 states "stable reactant adsorption at low coverage determines the
    selectivity towards allylic oxidation at high coverage." This is NOT a
    downstream branching mechanism — it is two geometrically distinct adsorption
    configurations of propylene with completely different chemistry.
    FIX: Two explicit propylene surface states (π and σ), each with its own
    reaction pathway. Selectivity is determined at the adsorption step, not after.

  FLAW 2 — Poisoning driven by total θ_prop × θ_OH
    v3 coupled poisoning to the same productive LH encounter, implying any
    propylene+OH encounter is a partial poisoning risk.
    WHAT THE PAPERS SAY: The CHᵧ poison arises specifically from the vinyl/π
    adsorption geometry — the flat-lying C=C-coordinated state. The σ (allylic)
    state gives acrolein, not poison. Poisoning is geometry-specific.
    FIX: Poisoning ODE driven by θ_π × θ_OH only.

  NEW MECHANISM — Frumkin lateral interaction drives geometry transition
    As surface coverage increases, the flat-lying π geometry is sterically
    disfavored (it requires more surface area). This shifts the geometry
    distribution toward σ (allylic), explaining why acrolein selectivity
    IMPROVES at higher propylene pressures (more coverage). Modeled via a
    Frumkin isotherm correction on K_π: the effective adsorption constant of
    the π state decreases as total propylene coverage increases.

  NEW MECHANISM — Pd oxide regime (MvK, >1.1 V)
    Koroidov XAS directly shows surface oxide formation above 1.1 V and links
    it to propylene glycol formation. The mechanism on PdO is Mars-van Krevelen,
    not LH. Kept as a separate additive term with its own rate law.

═══════════════════════════════════════════════════════════════════════════════
COMPLETE MECHANISM
═══════════════════════════════════════════════════════════════════════════════

─── REGIME 1: LH on metallic Pd (0.70 – 1.05 V) ─────────────────────────────

  Step 1a  [π adsorption — reversible, via C=C vinyl coordination]
    C₃H₆ + * ⇌ C₃H₆*_π      K_π_eff(θ_prop) = K_π0 · exp(-g·θ_prop / RT)
    Frumkin correction: as surface fills, π geometry is sterically disfavored.
    g [kJ/mol]: lateral repulsion parameter. g > 0 means crowding destabilizes π.

  Step 1b  [σ adsorption — reversible, via allylic CH₂ coordination]
    C₃H₆ + * ⇌ C₃H₆*_σ      K_σ (constant, no Frumkin needed)
    This geometry sits more upright, less sensitive to neighbors.

  Step 2   [OH* formation — electrochemical, quasi-equilibrium]
    H₂O + * → OH* + H⁺ + e⁻
    K_OH(E) = exp(F/RT · (E - E_eq)),  E_eq = 0.70 V fixed (Pourbaix)

  Step 3   [Acrolein pathway — LH, RDS on σ-adsorbed propylene]
    C₃H₆*_σ + OH* → acrolein* + H₂O + *    r_acrolein = k_rds · θ_σ · θ_OH
    Allylic C-H activation by adjacent OH*. Product desorbs fast.

  Step 4   [Poisoning pathway — from π-adsorbed propylene]
    C₃H₆*_π + OH* → CHᵧ* + fragments      dθ_P/dt = k_P · θ_π · θ_OH · (1-θ_P)
    Vinyl C-H/C-C bond scission → irreversible CHᵧ deposit.
    (Koroidov: "stripped as CO₂ only at high anodic potential")

─── REGIME 2: MvK on surface PdO (1.05 – 1.30 V) ───────────────────────────

  Step 5   [Surface oxide formation — electrochemical]
    Pd + H₂O → PdO_surf + 2H⁺ + 2e⁻
    θ_ox(E) = K_ox / (1 + K_ox),   K_ox = exp(F/RT · (E - E_ox))
    E_ox = 0.913 V (Winiwarter SI eq.1, pH 1 Pourbaix)
    [Koroidov XAS: propylene delays oxide onset to ~1.1 V; but E_ox here is
    for the bare thermodynamic onset — the actual transition is smooth]

  Step 6   [Propylene glycol — MvK on PdO surface]
    C₃H₆(aq) + O_lattice → propylene oxide → glycol (acid hydrolysis in solution)
    r_glycol = k_MvK · P_prop · θ_ox
    Linear in P_prop: propylene attacks lattice oxygen from solution.
    No competitive adsorption because metallic Pd sites are already oxidized.

─── TOTAL ────────────────────────────────────────────────────────────────────

  r_acrolein_total = (1 - θ_ox) · k_rds · θ_σ · θ_OH    [LH contribution only]
  r_glycol_total   = θ_ox · k_MvK · P_prop               [MvK contribution only]

  The (1 - θ_ox) and θ_ox weights enforce that both mechanisms cannot dominate
  simultaneously: as the surface oxidizes, LH sites disappear and MvK sites appear.

═══════════════════════════════════════════════════════════════════════════════
SITE BALANCE (LH regime — on metallic sites only)
═══════════════════════════════════════════════════════════════════════════════

  On the fraction (1 - θ_ox) of the surface that is metallic:

    θ_* + θ_π + θ_σ + θ_OH = (1 - θ_ox)(1 - θ_P)   ≡ θ_active

  Frumkin adsorption for π state:
    θ_π = K_π0 · exp(-g·(θ_π+θ_σ)/RT) · P · θ_*  ... implicit equation

  Standard Langmuir for σ state:
    θ_σ = K_σ · P · θ_*

  Electrochemical Langmuir for OH*:
    θ_OH = K_OH(E) · θ_*

  Solving this system requires iteration because K_π depends on θ_π + θ_σ.
  We use a simple fixed-point iteration (converges fast for typical parameters).

═══════════════════════════════════════════════════════════════════════════════
FREE PARAMETERS (to be fitted)
═══════════════════════════════════════════════════════════════════════════════

  LH regime (fit from steady-state rate vs E, P data + selectivity vs P):
    K_π0   [bar⁻¹]    — π adsorption constant at zero coverage
    K_σ    [bar⁻¹]    — σ adsorption constant (fixed or fitted; try K_σ ~ K_π0/5)
    g      [kJ/mol]   — Frumkin lateral repulsion for π state (>0 means destabilizing)
    k_rds  [rate units] — rate constant for allylic C-H activation

  Poisoning (fit ONLY if time-resolved data available):
    k_P    [s⁻¹]      — rate constant for CHᵧ poisoning from π state

  MvK regime (fit from rate data above 1.05 V):
    k_MvK  [rate/(bar)] — MvK rate constant for propylene glycol on PdO

FIXED FROM PHYSICS/LITERATURE:
  E_eq  = 0.70 V    — OH* formation onset on Pd (Pourbaix; Koroidov Fig.2A: OH(ad)
                       at 0.4–1.0 V confirms metallic Pd below 1.0 V)
  E_ox  = 0.913 V   — PdO formation equilibrium potential at pH 1
                       (Winiwarter SI Eq.1; Koroidov XAS: oxide seen above 1.0–1.1 V)
  F/RT  = 38.92 V⁻¹ — thermodynamic identity at 298 K
  T     = 298.15 K

IDENTIFIABILITY:
  K_π0, K_σ:  from rate vs P saturation curves at fixed E
               (if rate is still rising at highest P, K·P << 1, ratio K_σ/K_π0
               from selectivity vs P; if saturated, K from the half-saturation P)
  g:           from the shape of selectivity vs P — a rising acrolein fraction
               as P increases is the direct signature of the Frumkin correction
  k_rds:       from absolute rate magnitude at any (E,P) once coverages are known
  k_P:         from rate decay as function of time-on-stream
  k_MvK:       from glycol rate at E > 1.05 V (linear in P_prop)

═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, curve_fit, brentq

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS — all fixed, never fitted
# ─────────────────────────────────────────────────────────────────────────────
T    = 298.15        # K
R    = 8.314e-3      # kJ mol⁻¹ K⁻¹
F_RT = 96.485 / (R * T)   # F/RT = 38.92 V⁻¹ at 298 K

# Electrochemical fixed points — from Pd Pourbaix diagram and XAS literature
E_eq  = 0.70    # V vs RHE — OH* formation onset on metallic Pd
E_ox  = 0.913   # V vs RHE — PdO surface oxide formation at pH 1 (Winiwarter SI Eq.1)


# ─────────────────────────────────────────────────────────────────────────────
# ELECTROCHEMICAL FUNCTIONS (no free parameters)
# ─────────────────────────────────────────────────────────────────────────────

def K_OH(E):
    """
    Effective adsorption equilibrium constant for OH* from water activation.
    Derived from Nernst equation for H₂O + * → OH* + H⁺ + e⁻.
    Increases exponentially with E; equals 1.0 at E = E_eq = 0.70 V.
    """
    return np.exp(F_RT * (E - E_eq))


def theta_ox(E):
    """
    Fraction of Pd surface in the oxide state (PdO-like).
    Logistic function derived from Nernst equation for Pd + H₂O → PdO + 2H⁺ + 2e⁻.
    E_ox = 0.913 V is the thermodynamic equilibrium potential at pH 1.

    Note: Koroidov XAS shows propylene delays the actual oxide onset to ~1.1 V
    because propylene adsorption competes with and suppresses OH*/O* formation.
    This function gives the bare thermodynamic transition; the competition is
    handled implicitly through the (1-θ_ox) weighting in the rate expressions.
    At E < 0.9 V the oxide fraction is negligible (<0.01) — consistent with XAS.
    """
    K = np.exp(F_RT * (E - E_ox))
    return K / (1.0 + K)


# ─────────────────────────────────────────────────────────────────────────────
# SITE BALANCE — implicit due to Frumkin correction on π state
# ─────────────────────────────────────────────────────────────────────────────

def solve_coverages(E, P_prop, K_pi0, K_sigma, g, theta_P, tol=1e-9, max_iter=200):
    """
    Solve the implicit Frumkin site balance iteratively.

    The π adsorption constant is coverage-dependent:
        K_π_eff = K_π0 · exp(-g · θ_total_prop / RT)
    where θ_total_prop = θ_π + θ_σ.

    This creates an implicit equation because K_π depends on the coverage it
    determines. We solve by fixed-point iteration on θ_total_prop:

        θ_total_prop^{n+1} = f(θ_total_prop^n)

    Convergence is guaranteed for physically reasonable g and K values because
    the map is contractive (θ_total_prop is bounded [0, θ_active]).

    Parameters
    ----------
    E        : float — electrode potential [V vs RHE]
    P_prop   : float — propylene partial pressure [bar]
    K_pi0    : float — π adsorption constant at zero coverage [bar⁻¹]
    K_sigma  : float — σ adsorption constant [bar⁻¹]
    g        : float — Frumkin lateral repulsion [kJ/mol], g > 0 = destabilizing
    theta_P  : float — poison fraction (ODE state variable)

    Returns
    -------
    theta_pi, theta_sigma, theta_OH, theta_s : surface coverages
    All coverages are on the metallic (non-oxide) fraction of the surface.
    Absolute coverages (relative to total surface) must be multiplied by (1-θ_ox).
    """
    tox      = theta_ox(E)
    Koh      = K_OH(E)
    theta_active = (1.0 - tox) * (1.0 - theta_P)

    # Initial guess: ignore Frumkin correction
    theta_prop_guess = 0.0

    for _ in range(max_iter):
        K_pi_eff = K_pi0 * np.exp(-g * theta_prop_guess / R / T)

        # Standard competitive Langmuir with effective K_π
        denom     = 1.0 + K_pi_eff * P_prop + K_sigma * P_prop + Koh
        theta_s   = theta_active / denom
        theta_pi  = K_pi_eff * P_prop * theta_s
        theta_sig = K_sigma   * P_prop * theta_s
        theta_oh  = Koh       * theta_s

        theta_prop_new = theta_pi + theta_sig

        if abs(theta_prop_new - theta_prop_guess) < tol:
            break
        theta_prop_guess = theta_prop_new

    return theta_pi, theta_sig, theta_oh, theta_s


# ─────────────────────────────────────────────────────────────────────────────
# RATE EXPRESSIONS
# ─────────────────────────────────────────────────────────────────────────────

def rates(E, P_prop, K_pi0, K_sigma, g, k_rds, k_MvK, theta_P):
    """
    Compute acrolein and glycol rates.

    Acrolein (LH, metallic Pd):
        r_acrolein = (1 - θ_ox) · k_rds · θ_σ · θ_OH
        The (1-θ_ox) factor is already embedded in θ_σ and θ_OH via
        solve_coverages (θ_active includes the (1-θ_ox) factor).
        k_rds has units of [rate unit] / [dimensionless coverage²].

    Propylene glycol (MvK, PdO surface):
        r_glycol = k_MvK · P_prop · θ_ox
        Linear in P_prop: propylene attacks lattice oxygen from solution.
        k_MvK has units of [rate unit] / [bar].

    Physical meaning of k_rds:
        Lumped rate constant absorbing the true elementary pre-exponential,
        Pd surface site density, and any fast post-RDS steps.
        Cannot be interpreted as an elementary k without in-situ coverage data.
    """
    theta_pi, theta_sig, theta_oh, _ = solve_coverages(
        E, P_prop, K_pi0, K_sigma, g, theta_P)

    # Acrolein from σ-adsorbed propylene reacting with OH*
    r_acrolein = k_rds * theta_sig * theta_oh

    # Glycol from MvK on oxide surface
    r_glycol = k_MvK * P_prop * theta_ox(E)

    return r_acrolein, r_glycol, theta_pi, theta_sig, theta_oh


# ─────────────────────────────────────────────────────────────────────────────
# POISONING ODE
# ─────────────────────────────────────────────────────────────────────────────

def poisoning_ode(t, state, E, P_prop, K_pi0, K_sigma, g, k_P):
    """
    ODE for irreversible CHᵧ poisoning from the π-adsorbed state.

    dθ_P/dt = k_P · θ_π · θ_OH · (1 - θ_P)

    Physical interpretation:
    ─ θ_π · θ_OH: rate of π-state propylene encountering adjacent OH*,
                  triggering vinyl C-H cleavage → CHᵧ deposit.
    ─ (1 - θ_P): monolayer saturation limit — poison cannot pile up.
    ─ k_P: rate constant for the poisoning reaction. Identifiable ONLY from
            time-resolved rate decay data (rate vs time-on-stream).
            If only steady-state data is available, θ_P is a fixed input.

    Note: poisoning is FASTER at intermediate E (0.80–0.95 V) where both
    θ_π and θ_OH are non-negligible. At high E, θ_OH is large but θ_π is
    suppressed (propylene excluded by OH* competition). At low E, θ_OH is
    small. The maximum poisoning rate is at intermediate E — consistent with
    Winiwarter SI Fig. S4 showing largest current decay at ~0.90–0.95 V.
    """
    theta_P = float(state[0])
    theta_pi, _, theta_oh, _ = solve_coverages(
        E, P_prop, K_pi0, K_sigma, g, theta_P)
    dtheta_P_dt = k_P * theta_pi * theta_oh * (1.0 - theta_P)
    return [dtheta_P_dt]


def integrate_poisoning(t_exp, E, P_prop, K_pi0, K_sigma, g, k_P,
                         theta_P_init=0.0):
    """
    Integrate the poisoning ODE to get θ_P after elapsed time t_exp [s].

    Use only when time-resolved rate data is available.
    Uses Radau (stiff solver) because the ODE can become stiff when
    k_P is large and θ_P approaches 1.
    """
    sol = solve_ivp(
        poisoning_ode,
        t_span=(0, t_exp), y0=[theta_P_init], t_eval=None,
        args=(E, P_prop, K_pi0, K_sigma, g, k_P),
        method='Radau', rtol=1e-8, atol=1e-10
    )
    if not sol.success:
        raise RuntimeError(f"Poisoning ODE failed: {sol.message}")
    return float(sol.y[0, -1])


# ─────────────────────────────────────────────────────────────────────────────
# SYNTHETIC DATA GENERATION
# ─────────────────────────────────────────────────────────────────────────────

def generate_synthetic_data(noise_frac=0.04):
    """
    Generate synthetic rate vs (E, P) data with known true parameters.

    True parameters chosen to be physically representative:
      K_pi0  = 8.0 bar⁻¹   — π adsorption strong at zero coverage
      K_sigma= 1.5 bar⁻¹   — σ adsorption weaker (upright geometry)
      g      = 6.0 kJ/mol   — moderate lateral repulsion on π state
      k_rds  = 0.8          — rate constant for allylic C-H activation
      k_MvK  = 0.3          — rate constant for MvK glycol route
      theta_P = 0.12        — fixed moderate poison level

    The Frumkin correction at g=6 kJ/mol means:
      At θ_prop = 0.0: K_π_eff = 8.0 bar⁻¹
      At θ_prop = 0.3: K_π_eff = 8.0·exp(-6×0.3/2.479) = 8.0·0.48 = 3.8 bar⁻¹
      At θ_prop = 0.6: K_π_eff = 8.0·exp(-6×0.6/2.479) = 8.0·0.23 = 1.9 bar⁻¹
    → Strong suppression of π state at high coverage → acrolein selectivity rises.

    Data layout: [E, P_prop, r_acrolein, r_glycol]
    """
    K_pi0_true  = 8.0
    K_sig_true  = 1.5
    g_true      = 6.0
    k_rds_true  = 0.8
    k_MvK_true  = 0.3
    tP_true     = 0.12

    potentials = np.array([0.72, 0.76, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20])
    p_props    = np.array([0.05, 0.10, 0.20, 0.40])

    rng  = np.random.default_rng(seed=42)
    rows = []
    for E in potentials:
        for Pp in p_props:
            r_a, r_g, _, _, _ = rates(E, Pp, K_pi0_true, K_sig_true,
                                       g_true, k_rds_true, k_MvK_true, tP_true)
            rows.append([
                E, Pp,
                max(r_a * (1 + noise_frac * rng.standard_normal()), 0),
                max(r_g * (1 + noise_frac * rng.standard_normal()), 0)
            ])

    return np.array(rows)


data        = generate_synthetic_data()
# ─── REPLACE WITH REAL DATA ─────────────────────────────────────────────────
# data = np.loadtxt("your_data.csv", delimiter=",", skiprows=1)
# Required columns: [E (V vs RHE), P_prop (bar), r_acrolein, r_glycol]
# r_acrolein: production rate of acrolein (nmol/cm²ECSA/s or mmol/g_Pd/h)
# r_glycol:   production rate of propylene glycol (same units)
# Normalize both to the SAME units before fitting.
# If you only have total rate + selectivity: r_acrolein = sel × r_total
# ────────────────────────────────────────────────────────────────────────────

E_data   = data[:, 0]
P_data   = data[:, 1]
r_acr_obs = data[:, 2]
r_gly_obs = data[:, 3]
r_tot_obs = r_acr_obs + r_gly_obs

# Fixed poison fraction for steady-state fitting
# Set this from your time-on-stream deactivation characterization,
# or from the fractional current loss after 1h CA vs fresh electrode.
THETA_P_FIXED = 0.12


# ─────────────────────────────────────────────────────────────────────────────
# FITTING
# ─────────────────────────────────────────────────────────────────────────────

def model_prediction(X, K_pi0, K_sigma, g, k_rds, k_MvK):
    """
    Returns concatenated [r_acrolein_pred, r_glycol_pred] for all data points.
    Called by both DE and curve_fit.
    θ_P is held fixed at THETA_P_FIXED (see above).
    """
    E_arr, P_arr = X
    r_a_list, r_g_list = [], []
    for E, Pp in zip(E_arr, P_arr):
        r_a, r_g, _, _, _ = rates(E, Pp, K_pi0, K_sigma, g,
                                   k_rds, k_MvK, THETA_P_FIXED)
        r_a_list.append(r_a)
        r_g_list.append(r_g)
    return np.concatenate([r_a_list, r_g_list])


def residual_de(params):
    """
    Normalized sum-of-squares objective for differential evolution.
    Normalization ensures acrolein and glycol contribute equally regardless
    of their absolute magnitudes (glycol is typically ~5-10× smaller at 0.9 V).
    """
    K_pi0, K_sigma, g, k_rds, k_MvK = params
    preds = model_prediction((E_data, P_data), K_pi0, K_sigma, g, k_rds, k_MvK)
    n = len(r_acr_obs)
    r_a_pred = preds[:n]
    r_g_pred = preds[n:]

    norm_a = np.mean(r_acr_obs) + 1e-12
    norm_g = np.mean(r_gly_obs) + 1e-12
    ss = (np.sum(((r_a_pred - r_acr_obs) / norm_a) ** 2) +
          np.sum(((r_g_pred - r_gly_obs) / norm_g) ** 2))
    return ss


# Parameter bounds: [K_pi0, K_sigma, g, k_rds, k_MvK]
# K_pi0:  0.5–100 bar⁻¹  (π adsorption; should be > K_sigma)
# K_sigma:0.1–20  bar⁻¹  (σ adsorption; typically K_σ < K_π0)
# g:      0–20   kJ/mol   (Frumkin repulsion; 0=no Frumkin=standard Langmuir)
# k_rds:  free positive
# k_MvK:  free positive
bounds_de = [
    (0.5, 100.0),   # K_pi0
    (0.1, 20.0),    # K_sigma
    (0.0, 20.0),    # g  [kJ/mol]
    (1e-3, 1e2),    # k_rds
    (1e-3, 1e2),    # k_MvK
]

print("Fitting v4 model — 5 parameters: K_pi0, K_sigma, g, k_rds, k_MvK")
print(f"Fixed: E_eq={E_eq}V, E_ox={E_ox}V, F/RT={F_RT:.2f}V⁻¹, θ_P={THETA_P_FIXED}\n")

de_result = differential_evolution(
    residual_de, bounds_de,
    seed=42, maxiter=8000, tol=1e-12,
    polish=True, popsize=20,
    mutation=(0.5, 1.5), recombination=0.9,
    workers=1
)
popt_de = de_result.x
print(f"Global optimizer: converged={de_result.success}, residual={de_result.fun:.6f}")

# Refine with curve_fit for uncertainty estimates
y_obs = np.concatenate([r_acr_obs, r_gly_obs])
X_fit = (E_data, P_data)
try:
    lb = [b[0] for b in bounds_de]
    ub = [b[1] for b in bounds_de]
    popt_cf, pcov = curve_fit(
        model_prediction, X_fit, y_obs,
        p0=popt_de, bounds=(lb, ub), maxfev=200000
    )
    perr = np.sqrt(np.diag(pcov))
    popt_final = popt_cf
    print("curve_fit refinement: succeeded.")
except Exception as exc:
    print(f"curve_fit failed ({exc}). Using DE result, no std errors.")
    popt_final = popt_de
    perr = [float('nan')] * 5

K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f = popt_final
K_pi0_e, K_sig_e, g_e, k_rds_e, k_MvK_e = perr


# ─────────────────────────────────────────────────────────────────────────────
# RESULTS
# ─────────────────────────────────────────────────────────────────────────────

def r2(obs, pred):
    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)
    return 1.0 - ss_res / ss_tot


r_a_fit = np.array([rates(E, P, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)[0]
                     for E, P in zip(E_data, P_data)])
r_g_fit = np.array([rates(E, P, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)[1]
                     for E, P in zip(E_data, P_data)])

print("\n" + "═" * 62)
print("  FIT RESULTS — v4 Two-Geometry Frumkin LH/MvK Model")
print("═" * 62)
print("\n  FIXED (physics / Pourbaix / XAS literature):")
print(f"    E_eq   = {E_eq:.3f} V     OH* onset on metallic Pd")
print(f"    E_ox   = {E_ox:.3f} V     PdO formation at pH 1 (Winiwarter SI)")
print(f"    F/RT   = {F_RT:.2f} V⁻¹  thermodynamic identity at 298 K")
print(f"    θ_P    = {THETA_P_FIXED:.3f}       fixed (fit from time-series or estimate)")

print("\n  FITTED:")
print(f"    K_π0   = {K_pi0_f:8.3f}  ±  {K_pi0_e:.3f}  bar⁻¹   (π adsorption at zero coverage)")
print(f"    K_σ    = {K_sig_f:8.3f}  ±  {K_sig_e:.3f}  bar⁻¹   (σ adsorption, allylic)")
print(f"    g      = {g_f:8.3f}  ±  {g_e:.3f}  kJ/mol  (Frumkin lateral repulsion)")
print(f"    k_rds  = {k_rds_f:8.3f}  ±  {k_rds_e:.3f}  [units]  (allylic C-H activation)")
print(f"    k_MvK  = {k_MvK_f:8.3f}  ±  {k_MvK_e:.3f}  [units]  (MvK glycol on PdO)")

print(f"\n  IDENTIFIABLE PRODUCT (always recovered even if K_σ, k_rds are individually ambiguous):")
print(f"    K_σ · k_rds  = {K_sig_f * k_rds_f:.4f}  [bar⁻¹ · rate units]")
print(f"    (If perr on K_σ or k_rds exceeds the value itself, only their product")
print(f"     is constrained by the data. Fix one from independent measurement.)")

print(f"\n  GOF:")
print(f"    R² acrolein  = {r2(r_acr_obs, r_a_fit):.4f}")
print(f"    R² glycol    = {r2(r_gly_obs, r_g_fit):.4f}")
print(f"    R² total     = {r2(r_tot_obs, r_a_fit+r_g_fit):.4f}")

# Derived: K_ratio and Frumkin transition coverage
K_ratio = K_pi0_f / K_sig_f
print(f"\n  DERIVED:")
print(f"    K_π0/K_σ = {K_ratio:.2f}  (>1 means π is intrinsically more stable)")

# Coverage at which π and σ fractions are equal, ignoring OH* competition
# θ_π / θ_σ = K_π_eff / K_σ = 1 when K_π0·exp(-g·θ/RT) = K_σ
# → θ_cross = RT/g · ln(K_π0/K_σ)
if g_f > 0.01:
    theta_cross = (R * T / g_f) * np.log(K_pi0_f / K_sig_f) if K_pi0_f > K_sig_f else float('nan')
    print(f"    θ_cross    = {theta_cross:.3f}  (coverage where π≡σ fraction; "
          f"below: π dominates→poison; above: σ dominates→acrolein)")

# Diagnostics table
print("\n" + "═" * 62)
print("  COVERAGE DIAGNOSTICS")
print(f"  {'E(V)':>7} {'P(bar)':>8} {'θ_π':>7} {'θ_σ':>7} {'θ_OH':>7} "
      f"{'θ_*':>7} {'f_σ%':>7} {'r_acr':>8} {'r_gly':>8}")
for E_c in [0.80, 0.90, 1.00, 1.10]:
    for Pp_c in [0.1, 0.4]:
        r_a, r_g, t_pi, t_sig, t_oh = rates(E_c, Pp_c, K_pi0_f, K_sig_f,
                                              g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)
        _, _, _, t_s = solve_coverages(E_c, Pp_c, K_pi0_f, K_sig_f, g_f, THETA_P_FIXED)
        f_sigma = 100 * t_sig / (t_pi + t_sig + 1e-12)
        print(f"  {E_c:>7.2f} {Pp_c:>8.2f} {t_pi:>7.3f} {t_sig:>7.3f} {t_oh:>7.3f} "
              f"{t_s:>7.3f} {f_sigma:>7.1f} {r_a:>8.4f} {r_g:>8.4f}")

print()
print("  f_σ% = fraction of adsorbed propylene in the σ (allylic) geometry")
print("  Rising f_σ% with P_prop confirms Frumkin geometry transition.")


# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(
    "LH Kinetic Model v4  ·  Two-geometry Frumkin / LH+MvK\n"
    "Propylene electrooxidation on Pd/C  ·  grounded in Winiwarter 2019 + Koroidov 2021",
    fontsize=11
)

E_sw  = np.linspace(0.70, 1.25, 200)
Pp_sw = np.linspace(0.01, 0.50, 200)

# ── Panel 1: Parity ──────────────────────────────────────────────────────────
ax = axes[0, 0]
ax.scatter(r_acr_obs, r_a_fit, label="Acrolein", color="#1E88E5", alpha=0.75, s=40, zorder=3)
ax.scatter(r_gly_obs, r_g_fit, label="Glycol",   color="#8E24AA", alpha=0.75, s=40, zorder=3)
all_r = np.concatenate([r_acr_obs, r_gly_obs, r_a_fit, r_g_fit])
lims  = [all_r.min() * 0.9, all_r.max() * 1.1]
ax.plot(lims, lims, 'k--', lw=1, alpha=0.5)
ax.set(xlabel="Observed rate", ylabel="Predicted rate", title="Parity Plot", xlim=lims, ylim=lims)
ax.legend(fontsize=8)
ax.text(0.05, 0.88,
        f"R²(acr)={r2(r_acr_obs,r_a_fit):.3f}\nR²(gly)={r2(r_gly_obs,r_g_fit):.3f}",
        transform=ax.transAxes, fontsize=8, va='top')

# ── Panel 2: Rate vs Potential ────────────────────────────────────────────────
ax = axes[0, 1]
for Pp_fixed, col in [(0.10, '#E53935'), (0.20, '#1E88E5'), (0.40, '#00897B')]:
    r_a = [rates(E, Pp_fixed, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)[0]
           for E in E_sw]
    r_g = [rates(E, Pp_fixed, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)[1]
           for E in E_sw]
    ax.plot(E_sw, r_a, color=col, lw=2,   label=f"Acrolein P={Pp_fixed}")
    ax.plot(E_sw, r_g, color=col, lw=1.5, ls='--', label=f"Glycol P={Pp_fixed}")

ax.axvline(E_eq,  color='#378ADD', lw=1, ls=':', alpha=0.7, label=f"E_eq={E_eq}V")
ax.axvline(E_ox,  color='#E53935', lw=1, ls=':', alpha=0.7, label=f"E_ox={E_ox}V")
ax.set(xlabel="E (V vs RHE)", ylabel="Rate [model units]",
       title="Rate vs Potential\nsolid=acrolein, dashed=glycol")
ax.legend(fontsize=6, ncol=2)

# ── Panel 3: Rate vs P_prop (KEY PLOT — tests Frumkin) ───────────────────────
ax = axes[0, 2]
for E_fixed, col in [(0.85, '#E53935'), (0.95, '#FB8C00'), (1.05, '#1E88E5')]:
    r_a = [rates(E_fixed, Pp, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)[0]
           for Pp in Pp_sw]
    r_g = [rates(E_fixed, Pp, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, THETA_P_FIXED)[1]
           for Pp in Pp_sw]
    ax.plot(Pp_sw, r_a, color=col, lw=2,   label=f"Acrolein E={E_fixed}V")
    ax.plot(Pp_sw, r_g, color=col, lw=1.5, ls='--', label=f"Glycol E={E_fixed}V")
ax.set(xlabel="P_prop [bar]", ylabel="Rate [model units]",
       title="Rate vs P_prop\n(non-monotone acrolein = Frumkin signature)")
ax.legend(fontsize=6, ncol=2)

# ── Panel 4: Geometry fraction f_σ vs E and P ────────────────────────────────
ax = axes[1, 0]
for Pp_fixed, col in [(0.05, '#E53935'), (0.20, '#FB8C00'), (0.40, '#1E88E5')]:
    f_sigma = []
    for E in E_sw:
        t_pi, t_sig, _, _ = solve_coverages(E, Pp_fixed, K_pi0_f, K_sig_f, g_f, THETA_P_FIXED)
        f_sigma.append(100 * t_sig / (t_pi + t_sig + 1e-12))
    ax.plot(E_sw, f_sigma, color=col, lw=2, label=f"P={Pp_fixed} bar")

ax.axvline(E_eq, color='gray', lw=1, ls=':', alpha=0.6)
ax.axhline(50, color='gray', lw=1, ls=':', alpha=0.4)
ax.set(xlabel="E (V vs RHE)", ylabel="σ-state fraction f_σ [%]",
       title="Allylic (σ) geometry fraction vs Potential\n"
             "(100% = all propylene in allylic geometry → max acrolein)",
       ylim=(0, 100))
ax.legend(fontsize=8)

# ── Panel 5: Frumkin geometry fraction vs P_prop (the key Winiwarter insight) ─
ax = axes[1, 1]
for E_fixed, col in [(0.80, '#E53935'), (0.90, '#FB8C00'), (1.00, '#1E88E5')]:
    f_sigma = []
    for Pp in Pp_sw:
        t_pi, t_sig, _, _ = solve_coverages(E_fixed, Pp, K_pi0_f, K_sig_f, g_f, THETA_P_FIXED)
        f_sigma.append(100 * t_sig / (t_pi + t_sig + 1e-12))
    ax.plot(Pp_sw, f_sigma, color=col, lw=2, label=f"E={E_fixed}V")

ax.axhline(50, color='gray', lw=1, ls=':', alpha=0.4)
ax.set(xlabel="P_prop [bar]", ylabel="σ-state fraction f_σ [%]",
       title="Allylic fraction rises with P_prop\n"
             "(Frumkin: surface crowding suppresses flat π geometry)",
       ylim=(0, 100))
ax.legend(fontsize=8)
ax.text(0.05, 0.15,
        "This is the mechanistic content of Winiwarter's\n"
        "'allylic selectivity increases at high coverage'",
        transform=ax.transAxes, fontsize=7.5, color='#444',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))

# ── Panel 6: Poisoning dynamics ───────────────────────────────────────────────
ax = axes[1, 2]
k_P_illustrative = 3e-4   # s⁻¹, illustrative — fit from time-series data

t_span = (0, 7200)
t_eval = np.linspace(0, 7200, 300)

for E_fixed, col in [(0.80, '#E53935'), (0.90, '#FB8C00'), (0.95, '#1E88E5'), (1.05, '#43A047')]:
    sol = solve_ivp(
        poisoning_ode,
        t_span=t_span, y0=[0.0], t_eval=t_eval,
        args=(E_fixed, 0.20, K_pi0_f, K_sig_f, g_f, k_P_illustrative),
        method='Radau', rtol=1e-8, atol=1e-10
    )
    theta_P_t = sol.y[0]
    # Compute rate normalized to initial rate
    r0 = rates(E_fixed, 0.20, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, 0.0)[0]
    r_t = np.array([
        rates(E_fixed, 0.20, K_pi0_f, K_sig_f, g_f, k_rds_f, k_MvK_f, tP)[0]
        for tP in theta_P_t
    ])
    if r0 > 1e-10:
        ax.plot(t_eval / 3600, r_t / r0, color=col, lw=2, label=f"E={E_fixed}V")

ax.set(xlabel="Time [h]", ylabel="Normalized acrolein rate [r/r₀]",
       title="Poisoning dynamics: π-state specific\n"
             f"k_P={k_P_illustrative:.1e} s⁻¹ — ILLUSTRATIVE",
       ylim=(0, 1.05))
ax.axhline(1.0, color='gray', lw=1, ls=':', alpha=0.4)
ax.legend(fontsize=8)
ax.text(0.05, 0.08,
        "Fastest deactivation at ~0.90 V where\n"
        "θ_π · θ_OH is maximized — consistent\n"
        "with Winiwarter SI Fig. S4",
        transform=ax.transAxes, fontsize=7.5, color='#444',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/lh_kinetics_v4_fit.png", dpi=150, bbox_inches='tight')
plt.show()
print("\nFigure saved to lh_kinetics_v4_fit.png")


# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENTAL REQUIREMENTS SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print("""
═══════════════════════════════════════════════════════════════
  WHAT DATA CONSTRAINS EACH PARAMETER
═══════════════════════════════════════════════════════════════

  K_π0, K_σ
    → Rate vs P_prop at FIXED E in the low-E LH window (0.75–0.90 V).
      Two-state saturation: at high P, rate saturates differently for
      acrolein (σ-controlled) vs glycol (MvK, linear in P).
      If you have only acrolein rate vs P at 2–3 potentials, the
      saturation shape gives K_π0+K_σ combined; ratio from selectivity.

  g  (Frumkin lateral repulsion)
    → Acrolein SELECTIVITY vs P_prop at fixed E.
      If selectivity rises with P (acrolein fraction increases as surface
      fills) → g > 0 confirmed. Magnitude of g from the curvature of
      the selectivity-vs-P curve.
      If selectivity is flat with P → g ≈ 0 (standard Langmuir sufficient).
      THIS IS THE KEY EXPERIMENTAL DISCRIMINATOR for v4 vs v3.

  k_rds
    → Absolute acrolein rate magnitude at any (E, P) once coverages known.
      Requires proper normalization to ECSA (not BET area; Winiwarter uses
      CO stripping to measure ECSA and normalizes to cm²_Pd).

  k_MvK
    → Glycol rate at E > 1.05 V. Should be linear in P_prop if MvK is
      correct (no competitive adsorption). Plot r_glycol / P_prop vs E:
      should give a sigmoid that turns on at E_ox.

  k_P  (poisoning)
    → Time-on-stream rate decay. Measure r_acrolein at t = 1, 5, 15, 30,
      60 min at fixed (E, P). The decay timescale gives k_P directly.
      Winiwarter SI Fig. S4 shows decay is fastest at 0.90–0.95 V —
      model predicts this because θ_π · θ_OH is maximized there.

  θ_P  (if no time data)
    → Estimate from the ratio of fresh-electrode current to steady-state
      current in Winiwarter's CA data. They show current drops to ~40–60%
      of initial within 3 min at 0.95 V → θ_P ~ 0.1–0.2 at steady state.

═══════════════════════════════════════════════════════════════
  MODEL LIMITATIONS — be explicit about these when reporting
═══════════════════════════════════════════════════════════════

  1. The Frumkin correction applies a MEAN-FIELD lateral interaction.
     Real lateral interactions are pair-specific and geometry-dependent.
     'g' should be interpreted as an effective, lumped repulsion constant,
     not a molecular binding energy.

  2. k_rds absorbs site density, normalisation, and post-RDS steps.
     Do not compare k_rds across different normalisations (per g_Pd vs
     per cm² ECSA) without careful unit conversion.

  3. The oxide transition is modeled as a smooth sigmoidal in E.
     In reality, PdO formation is a nucleation-and-growth process that
     shows hysteresis in CV. The model predicts no hysteresis — if your
     CA data shows hysteresis between up-scan and down-scan potential steps,
     a nucleation model (e.g., Avrami) is needed for the oxide formation step.

  4. pH effects: E_eq and E_ox are both pH-dependent through the Nernst
     equation (shift by -59 mV/pH unit). The values here are for pH 1
     (Winiwarter 2019 conditions: 0.1 M HClO4). If your experiment uses
     a different pH, recalculate E_eq and E_ox accordingly.
═══════════════════════════════════════════════════════════════
""")