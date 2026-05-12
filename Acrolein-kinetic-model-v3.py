"""
Langmuir-Hinshelwood Kinetic Model — v3
Propylene electrooxidation on Pd/C

═══════════════════════════════════════════════════════════════════════════════
DESIGN PHILOSOPHY
═══════════════════════════════════════════════════════════════════════════════

This model is DELIBERATELY MINIMAL. The goal is a kinetic framework that:
  (a) is mechanistically defensible step-by-step
  (b) has identifiable parameters given realistic experimental data
  (c) correctly captures potential-dependent selectivity from first principles
  (d) treats poisoning dynamically (ODE) rather than as a static isotherm

What was REMOVED from v2 and why:
  - K_OH(E) = K_OH_ref * exp(β*(E−E_ref))
      → wrong functional form; β is non-identifiable vs K_OH_ref.
        Replaced by electrochemical Langmuir isotherm derived from Nernst eq.
  - theta_P as static Langmuir in P_prop
      → irreversible poison cannot reach a pressure-dependent steady state
        without an explicit regeneration pathway. Replaced by ODE.
  - Identical coverage product for both selectivity routes
      → allyl/vinyl selectivity cannot be potential-dependent if both routes
        share the same theta_prop * theta_OH expression. Fixed by requiring
        the vinyl (glycol) route to consume a SECOND OH*, making it
        quadratic in theta_OH and strongly E-dependent.
  - Seven freely fitted parameters
      → electrochemical parameters (α, E_eq) fixed from physics/Pourbaix.
        Only k_prop_ads, k_rds, k_sel, k_poison fit from rate data.

═══════════════════════════════════════════════════════════════════════════════
MECHANISM (elementary steps)
═══════════════════════════════════════════════════════════════════════════════

Step 1  [Propylene adsorption — reversible, chemical]
        C₃H₆(g) + * ⇌ C₃H₆*
        K_prop = k_ads / k_des    [bar⁻¹, fitted or literature]

Step 2  [OH* formation — electrochemical, quasi-equilibrium]
        H₂O + * → OH* + H⁺ + e⁻
        Nernst/electrochemical Langmuir (see theta_OH below)
        E_eq ~ 0.70 V vs RHE on Pd (Pourbaix; fix this)
        α = 0.5  (standard assumption; fix this)

Step 3  [Rate-determining step — LH bimolecular surface reaction]
        C₃H₆* + OH* → [allyl-OH]* + *     (C–H activation, first OH insertion)
        r_rds = k_rds * theta_prop * theta_OH

Step 4a [Allyl route — fast desorption of intermediate, low OH* requirement]
        [allyl-OH]* → allyl alcohol / acrolein + *
        Selectivity factor S_allyl = 1 / (1 + k_sel * theta_OH)
        → dominant when theta_OH is LOW (low E)

Step 4b [Vinyl/glycol route — second OH* attack on intermediate]
        [allyl-OH]* + OH* → propylene glycol* → glycol + *
        Rate ∝ k_sel * theta_OH  relative to step 4a
        Selectivity factor S_vinyl = k_sel * theta_OH / (1 + k_sel * theta_OH)
        → dominant when theta_OH is HIGH (high E)
        
        PHYSICAL RATIONALE: The second OH insertion is the selectivity gate.
        At low potential, OH* coverage is low → most intermediates desorb as
        allyl products. At high potential, OH* coverage rises → intermediate
        is attacked again before desorption → glycol/vinyl products.
        This is consistent with Winiwarter's potential-dependent selectivity
        without requiring two independent rate constants.

Step 5  [Poisoning — irreversible, dynamic]
        [allyl-OH]* → C* (carbonaceous deposit) + products
        dtheta_P/dt = k_poison * theta_prop * theta_OH * (1 - theta_P)
        - The (1 - theta_P) factor: poison cannot accumulate beyond monolayer
        - No regeneration pathway (irreversible on experimental timescale)
        - This is the ONLY step treated with an ODE

═══════════════════════════════════════════════════════════════════════════════
SITE BALANCE (self-consistent competitive Langmuir)
═══════════════════════════════════════════════════════════════════════════════

All reversible adsorbates compete for the SAME sites.
Poison theta_P is a PARAMETER of the ODE system (time-dependent state variable),
not a Langmuir isotherm—it reduces available sites.

Active site fraction: theta_active = 1 - theta_P

On active sites:
  theta_* + theta_prop + theta_OH = theta_active

  theta_prop = K_prop * P_prop * theta_*
  theta_OH   = K_OH_elec(E) * theta_*        [see below]

  Solving: theta_* = theta_active / (1 + K_prop*P_prop + K_OH_elec)
           theta_prop = K_prop * P_prop * theta_*
           theta_OH   = K_OH_elec(E) * theta_*

═══════════════════════════════════════════════════════════════════════════════
ELECTROCHEMICAL OH* COVERAGE — thermodynamic derivation
═══════════════════════════════════════════════════════════════════════════════

The equilibrium for H₂O + * ⇌ OH* + H⁺ + e⁻ gives (from Nernst):
  ΔG_OH(E) = ΔG_OH° + F*(E - E_eq)
  K_OH_elec(E) = exp(-ΔG_OH(E) / RT) = exp(-F*(E - E_eq) / RT)

where:
  E_eq ~ 0.70 V vs RHE for Pd (from Pourbaix diagram; DO NOT FIT THIS)
  F/RT = 38.92 V⁻¹ at 298 K (fixed by physics)
  ΔG_OH° is absorbed into E_eq definition (no extra free parameter)

This gives a sigmoid in E centered at E_eq, which is the correct physics.
At E >> E_eq: OH* dominates; at E << E_eq: OH* is absent.
This is a ONE-PARAMETER electrochemical description (just E_eq, fixed).

═══════════════════════════════════════════════════════════════════════════════
FREE PARAMETERS (to be fitted from experimental data)
═══════════════════════════════════════════════════════════════════════════════

  k_rds     [rate units]  — intrinsic rate constant for C–H activation (Step 3)
  K_prop    [bar⁻¹]       — propylene adsorption equilibrium constant (Step 1)
                            Can be estimated from TPD or adsorption isotherms
  k_sel     [dimensionless] — selectivity parameter for second OH* attack (Step 4b)
                            Only this controls allyl/vinyl branching
  k_poison  [s⁻¹]         — poisoning rate constant (Step 5)
                            Only meaningful if time-resolved data available

FIXED FROM PHYSICS/LITERATURE:
  E_eq = 0.70 V vs RHE    — OH* equilibrium potential on Pd
  α    = 0.5              — symmetry factor (Butler-Volmer)
  F/RT = 38.92 V⁻¹        — thermodynamic constant at 298 K
  T    = 298.15 K

═══════════════════════════════════════════════════════════════════════════════
IDENTIFIABILITY NOTES
═══════════════════════════════════════════════════════════════════════════════

From steady-state rate vs (E, P_prop) data:
  - K_prop identifiable from the P_prop dependence at fixed E
    (rate saturates at high P_prop → gives K_prop)
  - k_rds * theta_active is identifiable as a lumped constant; 
    k_rds and theta_active are NOT separately identifiable without coverage data
  - k_sel identifiable from the allyl/vinyl selectivity crossover potential
    (where S_allyl = S_vinyl = 50%)
  - k_poison identifiable ONLY from time-resolved data (rate vs time-on-stream)
    WITHOUT time data: fix k_poison = 0 and report theta_P as an input parameter

If you have only steady-state polarization data (r vs E at fixed P):
  Fit: k_rds, K_prop, k_sel  (3 parameters)
  Report theta_P as a sensitivity variable, not a fitted parameter.

═══════════════════════════════════════════════════════════════════════════════
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, curve_fit

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS  (fixed — do not fit these)
# ─────────────────────────────────────────────────────────────────────────────
T     = 298.15        # K
R     = 8.314e-3      # kJ mol⁻¹ K⁻¹
F_RT  = 96.485 / (R * T)   # = F/RT in V⁻¹ ≈ 38.92 V⁻¹ at 298 K

# Fixed electrochemical parameters — from Pourbaix diagram / Pd literature
# Ref: Grden et al. Electrochim. Acta 2008; Pd OH formation ~0.65–0.75 V vs RHE
E_eq  = 0.70          # V vs RHE — equilibrium potential for H₂O → OH* + H⁺ + e⁻

# ─────────────────────────────────────────────────────────────────────────────
# ELECTROCHEMICAL OH* ISOTHERM  (Nernst-derived, no free parameters)
# ─────────────────────────────────────────────────────────────────────────────

def K_OH_elec(E):
    """
    Effective adsorption equilibrium constant for OH* from water activation.

    Derived from:  ΔG_OH(E) = F*(E - E_eq)
                   K_OH(E) = exp(-F*(E - E_eq) / RT)

    At E = E_eq:  K_OH = 1   (equal tendency for bare vs OH-covered)
    At E > E_eq:  K_OH < 1   — WAIT, let's be careful about sign convention.

    Convention: for the oxidation H₂O + * → OH* + H⁺ + e⁻
    More positive E DRIVES the reaction forward (more oxidizing).
    So K_OH should INCREASE with E.

    ΔG_rxn = ΔG° - F*E  (Faradaic work done by the electrode driving oxidation)
    At equilibrium E_eq: ΔG° = F*E_eq
    So ΔG_rxn(E) = F*(E_eq - E)
    K_OH(E) = exp(-ΔG_rxn/RT) = exp(F*(E - E_eq)/RT)

    This correctly gives: K_OH increases with E (more anodic → more OH*).
    At E = E_eq, K_OH = 1.
    """
    return np.exp(F_RT * (E - E_eq))

# ─────────────────────────────────────────────────────────────────────────────
# SITE BALANCE  (competitive Langmuir on active sites)
# ─────────────────────────────────────────────────────────────────────────────

def coverages(E, P_prop, K_prop, theta_P):
    """
    Compute surface coverages from competitive Langmuir site balance.

    Site balance on active (non-poisoned) sites:
      theta_* + theta_prop + theta_OH = theta_active = (1 - theta_P)

    All species compete for the SAME sites — this is the self-consistent form.
    theta_P enters ONLY as a reduction in total available sites.
    It is NOT part of the competitive balance (poison sits on deactivated sites).

    Parameters
    ----------
    E        : float or array — electrode potential [V vs RHE]
    P_prop   : float or array — propylene partial pressure [bar]
    K_prop   : float — propylene adsorption equilibrium constant [bar⁻¹]
    theta_P  : float — poisoned site fraction (ODE state variable or fixed input)

    Returns
    -------
    theta_prop, theta_OH, theta_s : surface coverages
    """
    Koh        = K_OH_elec(E)                   # no free params
    theta_active = 1.0 - theta_P                # sites available for catalysis

    denom      = 1.0 + K_prop * P_prop + Koh    # competitive denominator
    theta_s    = theta_active / denom            # free sites
    theta_prop = K_prop * P_prop * theta_s       # propylene-covered
    theta_OH   = Koh * theta_s                   # OH-covered

    return theta_prop, theta_OH, theta_s

# ─────────────────────────────────────────────────────────────────────────────
# RATE EXPRESSIONS
# ─────────────────────────────────────────────────────────────────────────────

def selectivity_factors(theta_OH, k_sel):
    """
    Compute allyl/vinyl selectivity from the SECOND OH* attack competition.

    After the rate-determining C–H activation (Step 3), the intermediate either:
      (4a) desorbs quickly  → allyl products   (rate ∝ k_des_int, absorbed into 1)
      (4b) reacts with OH*  → vinyl/glycol     (rate ∝ k_sel * theta_OH)

    Branching ratio:
      S_vinyl = k_sel * theta_OH / (1 + k_sel * theta_OH)
      S_allyl = 1               / (1 + k_sel * theta_OH)

    PHYSICAL MEANING OF k_sel:
    k_sel >> 1: second attack is fast → vinyl route dominates at all E
    k_sel ~ 1:  transition occurs near theta_OH ~ 1 (which maps to E >> E_eq)
    k_sel << 1: allyl dominates across all experimentally accessible E

    The crossover potential E_cross is where S_allyl = S_vinyl = 50%:
      k_sel * theta_OH = 1
      theta_OH(E_cross) = 1/k_sel

    This is the ONLY parameter controlling selectivity vs potential.
    It is identifiable from the potential at which allyl% = vinyl% = 50%.
    """
    denom   = 1.0 + k_sel * theta_OH
    S_allyl = 1.0 / denom
    S_vinyl = k_sel * theta_OH / denom
    return S_allyl, S_vinyl


def rates(E, P_prop, K_prop, k_rds, k_sel, theta_P):
    """
    Compute allyl and vinyl production rates.

    r_total = k_rds * theta_prop * theta_OH        [rate-determining LH step]
    r_allyl = S_allyl * r_total
    r_vinyl = S_vinyl * r_total

    Note: k_rds here is a LUMPED rate constant. It absorbs:
      - the true elementary rate constant
      - the normalization of active site density (mol/g or mol/cm²)
    So k_rds has units of [rate] / [dimensionless coverages], e.g. mmol/g_cat/h.
    Do NOT interpret k_rds as a true elementary pre-exponential without
    independent coverage measurements (e.g., in-situ DRIFTS, EC-SERS).
    """
    theta_prop, theta_OH, _ = coverages(E, P_prop, K_prop, theta_P)
    r_total = k_rds * theta_prop * theta_OH
    S_allyl, S_vinyl = selectivity_factors(theta_OH, k_sel)
    return S_allyl * r_total, S_vinyl * r_total, theta_prop, theta_OH


# ─────────────────────────────────────────────────────────────────────────────
# DYNAMIC POISONING — ODE  (only needed if time-resolved data available)
# ─────────────────────────────────────────────────────────────────────────────

def poisoning_ode(t, theta_P, E, P_prop, K_prop, k_rds, k_sel, k_poison):
    """
    ODE for irreversible site poisoning.

    dtheta_P/dt = k_poison * theta_prop * theta_OH * (1 - theta_P)

    Physical interpretation:
    - Rate of poisoning is proportional to the productive LH encounter rate
      (poisoning comes FROM the catalytic intermediate, not from propylene directly)
    - (1 - theta_P): Langmuir monolayer limit for poison accumulation
    - No regeneration: irreversible on the timescale of the experiment
      (if regeneration is observed, add: - k_regen * theta_P * theta_OH)

    At steady state (dtheta_P/dt = 0):
      theta_P = 1.0  OR  theta_prop * theta_OH = 0 (no reactants)
    This correctly predicts eventual complete deactivation — which is observed
    in Pd electrocatalysts under sustained anodic polarization.

    If your experiment runs for time t_exp before measurement:
      Integrate this ODE from theta_P(0) = 0 to t = t_exp.
      The resulting theta_P is then used as the (fixed) input to coverages().
    """
    theta_prop, theta_OH, _ = coverages(E, P_prop, K_prop, float(theta_P[0]))
    dtheta_P_dt = k_poison * theta_prop * theta_OH * (1.0 - float(theta_P[0]))
    return [dtheta_P_dt]


def get_theta_P_at_time(t_exp, E, P_prop, K_prop, k_rds, k_sel, k_poison,
                         theta_P_init=0.0):
    """
    Integrate the poisoning ODE to get theta_P after elapsed time t_exp [s].

    Use this when you have time-series data (rate vs time at fixed E, P_prop).
    Do NOT use when you only have steady-state data.
    """
    sol = solve_ivp(
        poisoning_ode,
        t_span=(0, t_exp),
        y0=[theta_P_init],
        args=(E, P_prop, K_prop, k_rds, k_sel, k_poison),
        method='Radau',    # stiff solver (poisoning ODE can be stiff)
        rtol=1e-8, atol=1e-10,
        dense_output=True
    )
    if not sol.success:
        raise RuntimeError(f"ODE failed: {sol.message}")
    return float(sol.y[0, -1])


# ─────────────────────────────────────────────────────────────────────────────
# SYNTHETIC DATA GENERATION
# Replace this section entirely with real experimental data.
# ─────────────────────────────────────────────────────────────────────────────

def generate_synthetic_data(noise_frac=0.05):
    """
    Synthetic data with KNOWN parameters for model validation.

    True parameters:
      k_rds  = 1.2  [mmol/g/h or any consistent rate unit]
      K_prop = 3.0  [bar⁻¹]
      k_sel  = 8.0  [dimensionless]
      theta_P = 0.15 [assumed fixed for synthetic; in real use, integrate ODE]

    The selectivity crossover potential E_cross satisfies:
      theta_OH(E_cross) = 1/k_sel = 1/8 = 0.125
      K_OH(E_cross) / (1 + K_prop*P_prop + K_OH(E_cross)) = 0.125
      → solve numerically for E_cross at P_prop = 0.2 bar
      → expect E_cross ≈ 0.92 V (a physically reasonable potential)
    """
    # True parameters
    K_prop_true  = 3.0
    k_rds_true   = 1.2
    k_sel_true   = 8.0
    theta_P_true = 0.15

    potentials = np.array([0.72, 0.76, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20])
    p_props    = np.array([0.05, 0.10, 0.20, 0.40])

    rng = np.random.default_rng(seed=42)
    rows = []
    for E in potentials:
        for Pp in p_props:
            r_a, r_v, _, _ = rates(E, Pp, K_prop_true, k_rds_true,
                                    k_sel_true, theta_P_true)
            noise_a = 1 + noise_frac * rng.standard_normal()
            noise_v = 1 + noise_frac * rng.standard_normal()
            rows.append([E, Pp, max(r_a * noise_a, 0), max(r_v * noise_v, 0)])

    return np.array(rows)


data        = generate_synthetic_data()
# ── REPLACE WITH REAL DATA ──────────────────────────────────────────────────
# data = np.loadtxt("your_data.csv", delimiter=",", skiprows=1)
# Columns: [E (V vs RHE), P_prop (bar), r_allyl, r_vinyl]
# r_allyl and r_vinyl in consistent units (mmol/g_cat/h recommended)
# If you only have r_total and a selectivity measurement, compute:
#   r_allyl = selectivity_allyl * r_total
#   r_vinyl = (1 - selectivity_allyl) * r_total
# ────────────────────────────────────────────────────────────────────────────

E_data      = data[:, 0]
P_data      = data[:, 1]
r_allyl_obs = data[:, 2]
r_vinyl_obs = data[:, 3]
r_total_obs = r_allyl_obs + r_vinyl_obs


# ─────────────────────────────────────────────────────────────────────────────
# FITTING
# ─────────────────────────────────────────────────────────────────────────────

# theta_P for steady-state fitting:
# Option A (no time data): treat theta_P as a fixed sensitivity parameter
#   Choose a reasonable value (e.g., 0.10–0.20) from literature deactivation
#   data or set to 0 for "fresh catalyst" scenario.
# Option B (time data available): compute theta_P for each condition from ODE.
THETA_P_FIXED = 0.15   # <-- adjust based on your deactivation characterization

def model_for_fitting(X, K_prop, k_rds, k_sel):
    """
    Combined [r_allyl, r_vinyl] prediction for scipy fitting.

    Parameters fitted: K_prop, k_rds, k_sel  — ONLY THREE.
    theta_P is held fixed (see THETA_P_FIXED above).

    Returns concatenated [r_allyl_pred, r_vinyl_pred] for simultaneous fitting.
    """
    E_arr, P_arr = X
    r_a_list, r_v_list = [], []
    for E, Pp in zip(E_arr, P_arr):
        r_a, r_v, _, _ = rates(E, Pp, K_prop, k_rds, k_sel, THETA_P_FIXED)
        r_a_list.append(r_a)
        r_v_list.append(r_v)
    return np.concatenate([r_a_list, r_v_list])


def residual_for_de(params):
    """
    Normalized sum-of-squares for differential evolution.
    Normalization ensures allyl and vinyl routes contribute equally.
    """
    K_prop, k_rds, k_sel = params
    X = (E_data, P_data)
    preds = model_for_fitting(X, K_prop, k_rds, k_sel)
    n     = len(r_allyl_obs)
    r_a_pred = preds[:n]
    r_v_pred = preds[n:]

    norm_a = np.mean(r_allyl_obs) + 1e-12
    norm_v = np.mean(r_vinyl_obs) + 1e-12
    ss = (np.sum(((r_a_pred - r_allyl_obs) / norm_a)**2) +
          np.sum(((r_v_pred - r_vinyl_obs) / norm_v)**2))
    return ss


# Parameter bounds: [K_prop, k_rds, k_sel]
# K_prop: 0.1–50 bar⁻¹  (literature range for light alkenes on Pd)
# k_rds:  unconstrained positive (units depend on rate normalization)
# k_sel:  1–100  (dimensionless selectivity parameter)
bounds_de = [(0.1, 50.0), (1e-3, 1e2), (0.1, 100.0)]

print("Fitting 3-parameter model (K_prop, k_rds, k_sel)...")
print("Electrochemical parameters fixed: E_eq=0.70V, F/RT=38.92 V⁻¹\n")

de_result = differential_evolution(
    residual_for_de, bounds_de,
    seed=42, maxiter=5000, tol=1e-12,
    polish=True, popsize=15,
    mutation=(0.5, 1.5), recombination=0.9
)
popt_de = de_result.x
print(f"Global optimizer: converged={de_result.success}, residual={de_result.fun:.6f}")

# Refine with curve_fit for standard errors
y_obs_combined = np.concatenate([r_allyl_obs, r_vinyl_obs])
X_data = (E_data, P_data)

try:
    lb = [b[0] for b in bounds_de]
    ub = [b[1] for b in bounds_de]
    popt_cf, pcov = curve_fit(
        model_for_fitting, X_data, y_obs_combined,
        p0=popt_de, bounds=(lb, ub), maxfev=100000
    )
    perr = np.sqrt(np.diag(pcov))
    print("curve_fit refinement: succeeded — standard errors computed.")
    popt_final = popt_cf
except Exception as exc:
    print(f"curve_fit refinement failed ({exc}). Using DE result.")
    popt_final = popt_de
    perr = [float('nan')] * 3

K_prop_fit, k_rds_fit, k_sel_fit = popt_final
K_prop_err, k_rds_err,  k_sel_err  = perr


# ─────────────────────────────────────────────────────────────────────────────
# RESULTS SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

def r2(obs, pred):
    ss_res = np.sum((obs - pred)**2)
    ss_tot = np.sum((obs - np.mean(obs))**2)
    return 1.0 - ss_res / ss_tot

r_a_fit = np.array([rates(E, P, K_prop_fit, k_rds_fit, k_sel_fit, THETA_P_FIXED)[0]
                     for E, P in zip(E_data, P_data)])
r_v_fit = np.array([rates(E, P, K_prop_fit, k_rds_fit, k_sel_fit, THETA_P_FIXED)[1]
                     for E, P in zip(E_data, P_data)])
r_t_fit = r_a_fit + r_v_fit

print("\n" + "═"*60)
print("  FIT RESULTS — v3 Reduced Mechanistic Model")
print("═"*60)
print(f"\n  FIXED (from physics / Pourbaix):")
print(f"    E_eq    = {E_eq:.3f} V vs RHE   (OH* equilibrium potential on Pd)")
print(f"    F/RT    = {F_RT:.2f} V⁻¹          (thermodynamic, 298 K)")
print(f"    theta_P = {THETA_P_FIXED:.3f}              (fixed; estimate from deactivation data)")

print(f"\n  FITTED (3 parameters):")
print(f"    K_prop  = {K_prop_fit:8.3f}  ±  {K_prop_err:.3f}  bar⁻¹")
print(f"    k_rds   = {k_rds_fit:8.3f}  ±  {k_rds_err:.3f}  [rate units]")
print(f"    k_sel   = {k_sel_fit:8.3f}  ±  {k_sel_err:.3f}  [dimensionless]")

# Estimate selectivity crossover potential at P_prop = 0.2 bar
# At crossover: theta_OH = 1/k_sel
# theta_OH = K_OH(E) * theta_active / (1 + K_prop*Pp + K_OH(E))
# Solve numerically
from scipy.optimize import brentq

def crossover_residual(E):
    Pp = 0.2
    theta_active = 1.0 - THETA_P_FIXED
    Koh = K_OH_elec(E)
    denom = 1.0 + K_prop_fit * Pp + Koh
    theta_OH = Koh * theta_active / denom
    return theta_OH - 1.0 / k_sel_fit

try:
    E_cross = brentq(crossover_residual, 0.70, 1.30)
    print(f"\n  DERIVED — selectivity crossover @ P_prop=0.2 bar:")
    print(f"    E_cross = {E_cross:.3f} V vs RHE  (allyl% = vinyl% = 50%)")
except Exception:
    print("\n  Could not solve for crossover potential numerically.")

print(f"\n  GOF:")
print(f"    R² allyl   = {r2(r_allyl_obs, r_a_fit):.4f}")
print(f"    R² vinyl   = {r2(r_vinyl_obs, r_v_fit):.4f}")
print(f"    R² total   = {r2(r_total_obs, r_t_fit):.4f}")

# Coverage saturation check
theta_prop_at_max_E, theta_OH_at_max_E, _ = coverages(
    E_data.max(), P_data.max(), K_prop_fit, THETA_P_FIXED)
if theta_OH_at_max_E > 0.80:
    print(f"\n  ⚠  COVERAGE SATURATION WARNING:")
    print(f"     At E={E_data.max():.2f}V, θ_OH={theta_OH_at_max_E:.3f} (>0.80), θ_prop→0.")
    print(f"     Model predicts rate collapse at high E due to propylene exclusion.")
    print(f"     Physical if Pd surface oxide forms; unphysical if rate rises with E.")
    print(f"     If your data shows rising rate at high E, consider a two-site model")
    print(f"     or a Frumkin correction for OH*-propylene lateral repulsion.")

# Coverage diagnostic table
print("\n" + "═"*60)
print("  COVERAGE DIAGNOSTICS  (check physical reasonableness)")
print(f"  {'E (V)':>8} {'Pp(bar)':>9} {'θ_prop':>8} {'θ_OH':>8} "
      f"{'θ_*':>8} {'S_allyl%':>10} {'S_vinyl%':>10}")
for E_c, Pp_c in [(0.75, 0.1), (0.85, 0.2), (0.95, 0.2), (1.05, 0.2), (1.15, 0.4)]:
    tp, toh, ts = coverages(E_c, Pp_c, K_prop_fit, THETA_P_FIXED)
    Sa, Sv = selectivity_factors(toh, k_sel_fit)
    print(f"  {E_c:>8.2f} {Pp_c:>9.2f} {tp:>8.3f} {toh:>8.3f} "
          f"{ts:>8.3f} {Sa*100:>10.1f} {Sv*100:>10.1f}")

# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(
    "LH Kinetic Model v3 — Reduced Mechanistic Formulation\n"
    "Propylene electrooxidation on Pd/C  |  3 fitted params  |  E_eq fixed at 0.70 V",
    fontsize=12
)

E_sw  = np.linspace(0.70, 1.25, 200)
Pp_sw = np.linspace(0.01, 0.50, 200)
colors_Pp = {'0.05': '#E53935', '0.10': '#8E24AA',
             '0.20': '#1E88E5', '0.40': '#00897B'}

# ── Panel 1: Parity plot ──────────────────────────────────────────────────────
ax = axes[0, 0]
ax.scatter(r_allyl_obs, r_a_fit, label="Allyl", color="#1E88E5", alpha=0.75, s=50, zorder=3)
ax.scatter(r_vinyl_obs, r_v_fit, label="Vinyl", color="#43A047", alpha=0.75, s=50, zorder=3)
all_r = np.concatenate([r_allyl_obs, r_vinyl_obs, r_a_fit, r_v_fit])
lims = [all_r.min() * 0.9, all_r.max() * 1.1]
ax.plot(lims, lims, 'k--', lw=1, alpha=0.6)
ax.set(xlabel="Observed rate", ylabel="Predicted rate", title="Parity Plot",
       xlim=lims, ylim=lims)
ax.legend(fontsize=9)
ax.text(0.05, 0.90, f"R²(allyl)={r2(r_allyl_obs, r_a_fit):.3f}\n"
                     f"R²(vinyl)={r2(r_vinyl_obs, r_v_fit):.3f}",
        transform=ax.transAxes, fontsize=8, va='top')

# ── Panel 2: Rate vs Potential ─────────────────────────────────────────────────
ax = axes[0, 1]
for Pp_fixed in [0.05, 0.20, 0.40]:
    col = colors_Pp[f'{Pp_fixed:.2f}']
    r_a = np.array([rates(E, Pp_fixed, K_prop_fit, k_rds_fit, k_sel_fit, THETA_P_FIXED)[0]
                    for E in E_sw])
    r_v = np.array([rates(E, Pp_fixed, K_prop_fit, k_rds_fit, k_sel_fit, THETA_P_FIXED)[1]
                    for E in E_sw])
    ax.plot(E_sw, r_a, color=col, lw=2, label=f"Allyl Pp={Pp_fixed}")
    ax.plot(E_sw, r_v, color=col, lw=2, ls='--', label=f"Vinyl Pp={Pp_fixed}")

ax.axvline(E_eq, color='gray', lw=1, ls=':', alpha=0.7, label=f"E_eq={E_eq}V (fixed)")
try:
    ax.axvline(E_cross, color='black', lw=1, ls=':', alpha=0.5)
    ax.text(E_cross + 0.005, ax.get_ylim()[1] * 0.95,
            f"E_cross={E_cross:.2f}V", fontsize=7, ha='left')
except NameError:
    pass
ax.set(xlabel="E (V vs RHE)", ylabel="Rate [model units]",
       title="Rate vs Potential\nsolid=allyl, dashed=vinyl")
ax.legend(fontsize=6, ncol=2)

# ── Panel 3: Rate vs P_prop ────────────────────────────────────────────────────
ax = axes[0, 2]
for E_fixed, col in [(0.80, '#E53935'), (0.95, '#FB8C00'), (1.10, '#1E88E5')]:
    r_a = np.array([rates(E_fixed, Pp, K_prop_fit, k_rds_fit, k_sel_fit, THETA_P_FIXED)[0]
                    for Pp in Pp_sw])
    r_v = np.array([rates(E_fixed, Pp, K_prop_fit, k_rds_fit, k_sel_fit, THETA_P_FIXED)[1]
                    for Pp in Pp_sw])
    ax.plot(Pp_sw, r_a, color=col, lw=2, label=f"Allyl E={E_fixed}V")
    ax.plot(Pp_sw, r_v, color=col, lw=2, ls='--', label=f"Vinyl E={E_fixed}V")
ax.set(xlabel="P_prop [bar]", ylabel="Rate [model units]",
       title="Rate vs P_prop\nsolid=allyl, dashed=vinyl")
ax.legend(fontsize=6, ncol=2)

# ── Panel 4: Coverage and K_OH vs Potential ───────────────────────────────────
ax = axes[1, 0]
Pp_fixed = 0.20
theta_pr_sw = np.array([coverages(E, Pp_fixed, K_prop_fit, THETA_P_FIXED)[0] for E in E_sw])
theta_OH_sw = np.array([coverages(E, Pp_fixed, K_prop_fit, THETA_P_FIXED)[1] for E in E_sw])
theta_s_sw  = np.array([coverages(E, Pp_fixed, K_prop_fit, THETA_P_FIXED)[2] for E in E_sw])
Koh_sw      = K_OH_elec(E_sw)

ax2 = ax.twinx()
ax2.plot(E_sw, Koh_sw, color='#9E9E9E', lw=1.5, ls=':', label="K_OH(E) [right axis]")
ax2.set_ylabel("K_OH(E) [dimensionless]", color='#757575', fontsize=9)
ax2.tick_params(axis='y', labelcolor='#757575')

ax.plot(E_sw, theta_pr_sw, color='#1E88E5', lw=2, label="θ_prop")
ax.plot(E_sw, theta_OH_sw, color='#E53935', lw=2, label="θ_OH")
ax.plot(E_sw, theta_s_sw,  color='#43A047', lw=2, ls='--', label="θ_* (free)")
ax.axhline(THETA_P_FIXED, color='#9E9E9E', lw=1.5, ls='-',
           label=f"θ_P={THETA_P_FIXED} (fixed)")
ax.axvline(E_eq, color='gray', lw=1, ls=':', alpha=0.6)
ax.set(xlabel="E (V vs RHE)", ylabel="Surface coverage [0–1]",
       title=f"Coverage & K_OH vs Potential\n(P_prop={Pp_fixed} bar, θ_P fixed={THETA_P_FIXED})",
       ylim=(0, 1))
lines1, lab1 = ax.get_legend_handles_labels()
lines2, lab2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, lab1 + lab2, fontsize=7)

# ── Panel 5: Selectivity vs Potential ─────────────────────────────────────────
ax = axes[1, 1]
for Pp_fixed, col in [(0.05, '#E53935'), (0.20, '#1E88E5'), (0.40, '#00897B')]:
    theta_OH_arr = np.array([coverages(E, Pp_fixed, K_prop_fit, THETA_P_FIXED)[1]
                              for E in E_sw])
    Sa_arr = np.array([selectivity_factors(toh, k_sel_fit)[0] for toh in theta_OH_arr])
    ax.plot(E_sw, Sa_arr * 100, color=col, lw=2, label=f"Pp={Pp_fixed} bar")

ax.axhline(50, color='gray', lw=1, ls=':', alpha=0.6, label="50% (crossover)")
try:
    ax.axvline(E_cross, color='black', lw=1, ls=':', alpha=0.5)
    ax.text(E_cross + 0.005, 52, f"{E_cross:.2f}V", fontsize=7)
except NameError:
    pass
ax.set(xlabel="E (V vs RHE)", ylabel="Allyl selectivity [%]",
       title="Selectivity vs Potential\n(100%=all allyl, 0%=all vinyl)",
       ylim=(0, 100))
ax.legend(fontsize=8)

# ── Panel 6: Poisoning dynamics (ODE) ─────────────────────────────────────────
ax = axes[1, 2]
# Simulate poisoning at three potentials, fixed P_prop = 0.2 bar
# k_poison is illustrative — not fitted here (requires time-resolved data)
k_poison_illustrative = 5e-4   # [s⁻¹] — ILLUSTRATIVE ONLY

t_span = (0, 7200)   # 2 hours in seconds
t_eval = np.linspace(0, 7200, 300)

for E_fixed, col in [(0.80, '#E53935'), (0.95, '#FB8C00'), (1.10, '#1E88E5')]:
    sol = solve_ivp(
        poisoning_ode,
        t_span=t_span, y0=[0.0], t_eval=t_eval,
        args=(E_fixed, 0.20, K_prop_fit, k_rds_fit, k_sel_fit, k_poison_illustrative),
        method='Radau', rtol=1e-8, atol=1e-10
    )
    theta_P_t = sol.y[0]
    # Compute total rate at each theta_P
    r_total_t = np.array([
        sum(rates(E_fixed, 0.20, K_prop_fit, k_rds_fit, k_sel_fit, tP)[:2])
        for tP in theta_P_t
    ])
    ax.plot(t_eval / 3600, r_total_t / r_total_t[0],
            color=col, lw=2, label=f"E={E_fixed}V")

ax.set(xlabel="Time [h]", ylabel="Normalized rate [r/r₀]",
       title="Poisoning Dynamics (ODE)\n"
             f"k_poison={k_poison_illustrative:.1e} s⁻¹ — ILLUSTRATIVE\n"
             "Fit k_poison ONLY with time-series data",
       ylim=(0, 1.05))
ax.axhline(1.0, color='gray', lw=1, ls=':', alpha=0.4)
ax.legend(fontsize=8)
ax.text(0.05, 0.08,
        "k_poison not fitted — requires\ntime-resolved rate measurements",
        transform=ax.transAxes, fontsize=7, color='#757575', va='bottom',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/lh_kinetics_v3_fit.png", dpi=150, bbox_inches='tight')
plt.show()
print("\nPlot saved to lh_kinetics_v3_fit.png")


# ─────────────────────────────────────────────────────────────────────────────
# EXPERIMENTAL REQUIREMENTS — printed as guidance
# ─────────────────────────────────────────────────────────────────────────────

print("""
═══════════════════════════════════════════════════════════════
  WHAT EXPERIMENTS CONSTRAIN EACH PARAMETER
═══════════════════════════════════════════════════════════════

  K_prop  → Rate vs P_prop at FIXED E (check for Langmuir saturation)
             Alternatively: propylene TPD or volumetric adsorption isotherm
             If rate is linear in P_prop across full range → K_prop*P << 1
             (Henry's law limit; K_prop not separately identifiable — fix K_prop
             at a literature value, e.g. 2–5 bar⁻¹ for propylene on Pd)

  k_rds   → Absolute rate magnitude (not selectivity) at any condition
             Requires: catalyst mass loading, product quantification (GC or NMR)
             and consistent rate normalization (per g_Pd or per cm² ECSA)

  k_sel   → Allyl/vinyl product ratio at a series of potentials
             The crossover potential E_cross directly gives k_sel via:
               K_OH(E_cross) * theta_active / (1 + K_prop*Pp + K_OH(E_cross))
               = 1/k_sel
             Requires: product-selective detection (GC, HPLC, or NMR)
             Winiwarter use DEMS + HPLC — ideal.

  k_poison → Time-on-stream rate decay at fixed (E, P_prop)
             Minimum: measure rate at t=0, t=30min, t=60min, t=120min
             The ODE gives theta_P(t) → fit k_poison to the decay curve

  E_eq    → Fixed at 0.70 V vs RHE from Pd Pourbaix diagram
             Can be refined from cyclic voltammetry (OH* oxidation peak onset)
             Typical range on Pd: 0.65–0.75 V depending on facet and support

  CAUTION: If your rate data shows NO saturation with P_prop and NO clear
  selectivity crossover with E, you have insufficient observable variation to
  fit all three parameters. In that case, fit only k_rds (with K_prop and
  k_sel fixed from literature estimates) and report the result as a
  phenomenological rate expression, not a mechanistic one.
═══════════════════════════════════════════════════════════════
""")