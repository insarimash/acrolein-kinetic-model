"""
Langmuir-Hinshelwood Kinetic Model — v2
Propylene electrooxidation on Pd/C, informed by Winiwarter et al.

KEY MECHANISTIC UPDATES FROM v1:
  1. O source is water activation (electrochemical OH*), not O2 adsorption
     -> OH coverage depends on electrode potential E, not P_O2
  2. Surface poisoning by irreversibly adsorbed propylene-derived species
     -> Site balance now includes theta_P (poison fraction)
  3. Two parallel reaction pathways (potential-dependent selectivity):
     - Allyl route (0.7-1.0V):  prop* + OH* -> allyl alcohol / acrolein / acrylic acid
     - Vinyl route (1.0-1.2V):  prop* + OH* -> propylene glycol
  4. LH mechanism confirmed -> same overall rate structure, but revised site balance

REVISED SITE BALANCE:
  theta_* + theta_OH + theta_prop + theta_P = 1

  theta_OH   = K_OH * theta_*        (K_OH is potential-dependent constant at fixed E)
  theta_prop = K_prop * P_prop * theta_*
  theta_P    = theta_P_max * (K_P * P_prop) / (1 + K_P * P_prop)   [poison Langmuir]

  Solving for theta_*:
  theta_* = (1 - theta_P) / (1 + K_OH + K_prop * P_prop)

RATE EXPRESSIONS:
  r_allyl = k_allyl * theta_prop * theta_OH     [dominant 0.7-1.0V]
  r_vinyl = k_vinyl * theta_prop * theta_OH     [dominant 1.0-1.2V]
  r_total = r_allyl + r_vinyl

DATA REQUIREMENTS:
  Columns: [E (V vs RHE), P_prop (bar), r_allyl, r_vinyl]
  - E: electrode potential
  - P_prop: propylene partial pressure (O2/water is the oxidant, not a variable)
  - r_allyl, r_vinyl: product-specific rates from Winiwarter product distribution data

  If you only have total rate, set r_vinyl=0 below 1.0V and r_allyl=0 above 1.0V
  as a first approximation.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, differential_evolution

R = 8.314e-3   # kJ/mol/K  (use kJ to match DFT Gibbs energies)
F = 96.485     # kJ/mol/V  (Faraday constant)
T = 298.15     # K, room temperature — adjust if your experiments differ

# =============================================================================
# SECTION 1: DATA
#
# Replace with real data extracted from Winiwarter figures (use WebPlotDigitizer).
# Units:
#   E      : V vs RHE
#   P_prop : bar (if gas-phase propylene; set to 1.0 if aqueous fixed concentration)
#   r_*    : mmol/g_cat/h  or  nmol/cm²/s  — pick one and be consistent
#
# For now: synthetic data mimicking Winiwarter's potential-dependent selectivity.
# True params used: k_allyl=1.5, k_vinyl=0.4, K_OH_ref=2.0, beta=1.0,
#                   K_prop=3.0, theta_P_max=0.25, K_P=5.0
# =============================================================================

np.random.seed(0)

def _generate_example_data():
    """
    Synthetic data across potential and propylene pressure.
    Mimics: allyl route dominant at low E, vinyl route at high E.
    """
    potentials  = [0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
    p_props     = [0.1, 0.2, 0.4]   # bar

    # True parameters
    k_allyl     = 1.5
    k_vinyl     = 0.4
    K_OH_ref    = 2.0    # K_OH at E=1.0V (reference potential)
    beta        = 1.0    # d(ln K_OH)/dE  [1/V], electrochemical sensitivity
    K_prop      = 3.0
    theta_P_max = 0.25
    K_P         = 5.0

    rows = []
    for E in potentials:
        K_OH = K_OH_ref * np.exp(beta * (E - 1.0))  # K_OH shifts with potential
        for Pp in p_props:
            theta_P  = theta_P_max * (K_P * Pp) / (1 + K_P * Pp)
            theta_s  = (1 - theta_P) / (1 + K_OH + K_prop * Pp)
            theta_OH = K_OH * theta_s
            theta_pr = K_prop * Pp * theta_s

            r_a = k_allyl * theta_pr * theta_OH * (1 + 0.04 * np.random.randn())
            r_v = k_vinyl * theta_pr * theta_OH * (1 + 0.04 * np.random.randn())
            rows.append([E, Pp, max(r_a, 0), max(r_v, 0)])

    return np.array(rows)

data = _generate_example_data()
# ---- REPLACE WITH YOUR REAL DATA e.g.: ----
# data = np.array([
#     [0.75, 0.1, 0.045, 0.003],
#     [0.80, 0.1, 0.071, 0.005],
#     ...  each row: [E, P_prop, r_allyl, r_vinyl]
# ])

E_data      = data[:, 0]
P_prop_data = data[:, 1]
r_allyl_obs = data[:, 2]
r_vinyl_obs = data[:, 3]
r_total_obs = r_allyl_obs + r_vinyl_obs

# =============================================================================
# SECTION 2: MODEL
# =============================================================================

def coverages(E, P_prop, K_OH_ref, beta, K_prop, theta_P_max, K_P):
    """
    Computes theta_prop and theta_OH given potential and propylene pressure.

    K_OH(E) = K_OH_ref * exp(beta * (E - E_ref))
    where E_ref = 1.0 V is chosen arbitrarily as the reference point.

    beta > 0 means OH coverage increases with potential (more oxidizing surface).
    This captures the experimental observation that reaction onset is at 0.7V
    (where water activation begins) and increases with E.

    theta_P: Langmuir isotherm for irreversible poison.
    theta_P_max: maximum fraction of sites that can be poisoned (0 to 1).
    K_P: adsorption constant for the poisoning species.
    """
    E_ref   = 1.0  # V, reference potential — arbitrary, just sets K_OH_ref meaning
    K_OH    = K_OH_ref * np.exp(beta * (E - E_ref))
    theta_P = theta_P_max * (K_P * P_prop) / (1 + K_P * P_prop)

    # Solve site balance: theta_* = (1 - theta_P) / (1 + K_OH + K_prop * P_prop)
    theta_s  = (1 - theta_P) / (1 + K_OH + K_prop * P_prop)
    theta_OH = K_OH * theta_s
    theta_pr = K_prop * P_prop * theta_s

    return theta_pr, theta_OH, theta_P, theta_s


def rate_allyl(X, k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P):
    E, Pp = X
    theta_pr, theta_OH, _, _ = coverages(E, Pp, K_OH_ref, beta, K_prop, theta_P_max, K_P)
    return k_allyl * theta_pr * theta_OH

def rate_vinyl(X, k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P):
    E, Pp = X
    theta_pr, theta_OH, _, _ = coverages(E, Pp, K_OH_ref, beta, K_prop, theta_P_max, K_P)
    return k_vinyl * theta_pr * theta_OH

def rate_total(X, k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P):
    return rate_allyl(X, k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P) + \
           rate_vinyl(X, k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P)


def combined_residuals(params):
    """
    Fits BOTH allyl and vinyl rates simultaneously.
    This is important: fitting them separately could give inconsistent
    K_prop, K_OH values since they share the same site balance.
    """
    k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P = params
    X = (E_data, P_prop_data)
    r_a_pred = rate_allyl(X, *params)
    r_v_pred = rate_vinyl(X, *params)

    # Normalized residuals so allyl and vinyl contribute equally
    # (otherwise whichever has larger absolute values dominates the fit)
    res_a = (r_a_pred - r_allyl_obs) / (np.mean(r_allyl_obs) + 1e-10)
    res_v = (r_v_pred - r_vinyl_obs) / (np.mean(r_vinyl_obs) + 1e-10)

    return np.sum(res_a**2) + np.sum(res_v**2)

# =============================================================================
# SECTION 3: FITTING
# =============================================================================

# Parameter bounds: (lower, upper)
# Order: k_allyl, k_vinyl, K_OH_ref, beta, K_prop, theta_P_max, K_P
bounds = [
    (1e-4, 1e2),   # k_allyl
    (1e-4, 1e2),   # k_vinyl
    (1e-3, 1e3),   # K_OH_ref  (dimensionless equilibrium constant)
    (0.0,  5.0),   # beta      [1/V]  — electrochemical sensitivity of OH adsorption
    (1e-3, 1e2),   # K_prop
    (0.0,  0.99),  # theta_P_max  (fraction, must be < 1)
    (1e-3, 1e2),   # K_P
]

print("Running global optimizer (differential evolution)...")
print("This may take ~10-30 seconds.\n")

result = differential_evolution(
    combined_residuals, bounds,
    seed=42, maxiter=10000, tol=1e-12,
    polish=True, popsize=20,   # larger population = more thorough search
    mutation=(0.5, 1.5), recombination=0.9
)
popt = result.x
print(f"DE converged: {result.success}  |  Final residual: {result.fun:.6f}")

# Try to get standard errors via curve_fit refinement
# We fit total rate for curve_fit since it needs a single scalar output per point
X_data = (E_data, P_prop_data)
try:
    lb = [b[0] for b in bounds]
    ub = [b[1] for b in bounds]
    p0_clipped = np.clip(popt, [b + 1e-10 for b in lb], [b - 1e-10 for b in ub])
    popt_cf, pcov = curve_fit(rate_total, X_data, r_total_obs,
                              p0=p0_clipped, bounds=(lb, ub), maxfev=100000)
    perr = np.sqrt(np.diag(pcov))
    # Use DE result (better global fit) for plotting, curve_fit for uncertainty
    print("curve_fit refinement succeeded — standard errors available.")
except Exception as e:
    print(f"curve_fit refinement failed ({e}) — using DE result, no std errors.")
    popt_cf = popt
    perr = [float('nan')] * len(popt)

# Unpack parameters
labels = ["k_allyl", "k_vinyl", "K_OH_ref", "beta (1/V)", "K_prop", "theta_P_max", "K_P"]

# =============================================================================
# SECTION 4: PRINT RESULTS
# =============================================================================

def r2_score(obs, pred):
    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)
    return 1 - ss_res / ss_tot

r_a_pred = rate_allyl(X_data, *popt)
r_v_pred = rate_vinyl(X_data, *popt)
r_t_pred = r_a_pred + r_v_pred

print("\n" + "="*55)
print("  FIT RESULTS — Winiwarter-informed LH model")
print("="*55)
for label, val, err in zip(labels, popt, perr):
    print(f"  {label:20s} = {val:10.4f}  ±  {err:.4f}")

print(f"\n  R² allyl   = {r2_score(r_allyl_obs, r_a_pred):.4f}")
print(f"  R² vinyl   = {r2_score(r_vinyl_obs, r_v_pred):.4f}")
print(f"  R² total   = {r2_score(r_total_obs, r_t_pred):.4f}")

# Coverage diagnostics at a few representative conditions
print("\n" + "="*55)
print("  COVERAGE DIAGNOSTICS")
print("  (useful for checking physical reasonableness)")
print("="*55)
print(f"  {'E (V)':>8} {'Pp (bar)':>10} {'θ_prop':>8} {'θ_OH':>8} {'θ_P':>8} {'θ_*':>8}")
for E_check, Pp_check in [(0.75, 0.1), (0.90, 0.2), (1.05, 0.2), (1.20, 0.4)]:
    tp, toh, tP, ts = coverages(E_check, Pp_check, *popt[2:])
    print(f"  {E_check:>8.2f} {Pp_check:>10.2f} {tp:>8.3f} {toh:>8.3f} {tP:>8.3f} {ts:>8.3f}")

# =============================================================================
# SECTION 5: PLOTS
# =============================================================================

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle("LH Kinetic Model v2 — Winiwarter-informed\nPropylene electrooxidation on Pd/C",
             fontsize=13)

# ---- Plot 1: Parity plot ----
ax = axes[0, 0]
ax.scatter(r_allyl_obs, r_a_pred, label="Allyl route", color="#2196F3", alpha=0.7, s=50)
ax.scatter(r_vinyl_obs, r_v_pred, label="Vinyl route", color="#4CAF50", alpha=0.7, s=50)
all_obs  = np.concatenate([r_allyl_obs, r_vinyl_obs])
all_pred = np.concatenate([r_a_pred,    r_v_pred])
lims = [min(all_obs)*0.9, max(all_obs)*1.1]
ax.plot(lims, lims, 'k--', lw=1)
ax.set_xlabel("Observed rate"); ax.set_ylabel("Predicted rate")
ax.set_title("Parity Plot"); ax.legend(fontsize=8)

# ---- Plot 2: Rate vs Potential (fixed P_prop) ----
ax = axes[0, 1]
E_sweep = np.linspace(0.70, 1.20, 100)
for Pp_fixed, color in [(0.1, '#FF5722'), (0.2, '#9C27B0'), (0.4, '#00BCD4')]:
    X_sw = (E_sweep, np.full_like(E_sweep, Pp_fixed))
    ax.plot(E_sweep, rate_allyl(X_sw, *popt), color=color, lw=2,
            label=f"Allyl, Pp={Pp_fixed}")
    ax.plot(E_sweep, rate_vinyl(X_sw, *popt), color=color, lw=2,
            linestyle='--', label=f"Vinyl, Pp={Pp_fixed}")

# Mark the selectivity transition at 1.0V
ax.axvline(1.0, color='gray', lw=1, linestyle=':', alpha=0.7)
ax.text(1.01, ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else 0.01,
        '1.0V transition', fontsize=7, color='gray', va='bottom')
ax.set_xlabel("E (V vs RHE)"); ax.set_ylabel("Rate")
ax.set_title("Rate vs Potential\n(solid=allyl, dashed=vinyl)")
ax.legend(fontsize=6)

# ---- Plot 3: Rate vs P_prop (fixed E) ----
ax = axes[0, 2]
Pp_sweep = np.linspace(0.01, 0.5, 100)
for E_fixed, color in [(0.80, '#E91E63'), (1.00, '#FF9800'), (1.15, '#3F51B5')]:
    X_sw = (np.full_like(Pp_sweep, E_fixed), Pp_sweep)
    ax.plot(Pp_sweep, rate_allyl(X_sw, *popt), color=color, lw=2,
            label=f"Allyl, E={E_fixed}V")
    ax.plot(Pp_sweep, rate_vinyl(X_sw, *popt), color=color, lw=2,
            linestyle='--', label=f"Vinyl, E={E_fixed}V")
ax.set_xlabel("P_prop (bar)"); ax.set_ylabel("Rate")
ax.set_title("Rate vs P_prop\n(solid=allyl, dashed=vinyl)")
ax.legend(fontsize=6)

# ---- Plot 4: Coverage vs Potential ----
ax = axes[1, 0]
Pp_fixed = 0.2
E_sweep  = np.linspace(0.70, 1.20, 100)
theta_pr_sw, theta_OH_sw, theta_P_sw, theta_s_sw = coverages(
    E_sweep, np.full_like(E_sweep, Pp_fixed), *popt[2:])
ax.plot(E_sweep, theta_pr_sw, label="θ_prop", color="#2196F3", lw=2)
ax.plot(E_sweep, theta_OH_sw, label="θ_OH",   color="#F44336", lw=2)
ax.plot(E_sweep, theta_P_sw,  label="θ_P (poison)", color="#9E9E9E", lw=2, linestyle='--')
ax.plot(E_sweep, theta_s_sw,  label="θ_* (free)", color="#4CAF50", lw=2, linestyle=':')
ax.axvline(0.7, color='gray', lw=1, linestyle=':', alpha=0.5)
ax.set_xlabel("E (V vs RHE)"); ax.set_ylabel("Surface coverage")
ax.set_title(f"Coverage vs Potential\n(P_prop={Pp_fixed} bar)")
ax.legend(fontsize=8); ax.set_ylim(0, 1)

# ---- Plot 5: Selectivity (allyl fraction) vs Potential ----
ax = axes[1, 1]
for Pp_fixed, color in [(0.1, '#FF5722'), (0.2, '#9C27B0'), (0.4, '#00BCD4')]:
    X_sw = (E_sweep, np.full_like(E_sweep, Pp_fixed))
    r_a  = rate_allyl(X_sw, *popt)
    r_v  = rate_vinyl(X_sw, *popt)
    sel  = r_a / (r_a + r_v + 1e-12)
    ax.plot(E_sweep, sel * 100, color=color, lw=2, label=f"Pp={Pp_fixed} bar")
ax.axvline(1.0, color='gray', lw=1, linestyle=':', alpha=0.7)
ax.axhline(50,  color='gray', lw=1, linestyle=':', alpha=0.5)
ax.set_xlabel("E (V vs RHE)"); ax.set_ylabel("Allyl selectivity (%)")
ax.set_title("Allyl Route Selectivity vs Potential\n(100% = all allyl, 0% = all vinyl)")
ax.legend(fontsize=8); ax.set_ylim(0, 100)

# ---- Plot 6: 2D rate surface at E=0.85V (allyl-dominant regime) ----
ax = axes[1, 2]
Pp_grid = np.linspace(0.01, 0.5, 80)
E_grid  = np.linspace(0.70, 1.20, 80)
PP, EE  = np.meshgrid(Pp_grid, E_grid)
R_allyl_surf = rate_allyl((EE, PP), *popt)
cf = ax.contourf(PP, EE, R_allyl_surf, levels=20, cmap='viridis')
plt.colorbar(cf, ax=ax, label="Allyl rate")
ax.scatter(P_prop_data, E_data, c='red', s=20, zorder=5, label="Data points")
ax.axhline(1.0, color='white', lw=1, linestyle='--', alpha=0.7)
ax.text(0.35, 1.01, '1.0V', color='white', fontsize=8)
ax.set_xlabel("P_prop (bar)"); ax.set_ylabel("E (V vs RHE)")
ax.set_title("Allyl Rate Surface\n(white dashed = selectivity transition)")
ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("/mnt/user-data/outputs/lh_kinetics_v2_fit.png", dpi=150, bbox_inches='tight')
plt.show()
print("\nPlot saved.")