"""
Langmuir-Hinshelwood Kinetic Model Fitting
Propylene oxidation to acrolein over Pd/C

Model:
    r = (k * K_prop * P_prop * K_O * P_O2) / (1 + K_prop*P_prop + K_O*P_O2)^2

    where k, K_prop, K_O are the fit parameters at a given temperature.
    (Assumes: molecular O2 adsorption — simple first model. See note on dissociative.)

Usage:
    1. Replace the EXAMPLE DATA section with your real data (or literature data).
    2. Run: python lh_kinetics.py
    3. Check the parity plot and residuals — if bad, try switching to the dissociative model.

Requirements: numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, differential_evolution
from scipy.stats import pearsonr

# =============================================================================
# SECTION 1: DATA
# Replace this with your actual data from literature or experiments.
# Units: pressures in bar (or atm, just be consistent), rate in mmol/g_cat/h
# =============================================================================

# Each row: [P_prop (bar), P_O2 (bar), rate (mmol/g/h)]
# This is SYNTHETIC example data — replace with real literature values
# Generated from the LH model with k=2.0, K_prop=3.0, K_O=1.5, plus 5% noise

np.random.seed(42)

def _generate_example_data():
    k_true, Kp_true, Ko_true = 2.0, 3.0, 1.5
    conditions = [
        (0.05, 0.10), (0.05, 0.20), (0.05, 0.40),
        (0.10, 0.10), (0.10, 0.20), (0.10, 0.40),
        (0.20, 0.10), (0.20, 0.20), (0.20, 0.40),
        (0.30, 0.10), (0.30, 0.20), (0.30, 0.40),
        (0.40, 0.10), (0.40, 0.20), (0.40, 0.40),
    ]
    data = []
    for Pp, Po in conditions:
        r_true = (k_true * Kp_true * Pp * Ko_true * Po) / (1 + Kp_true*Pp + Ko_true*Po)**2
        r_noisy = r_true * (1 + 0.05 * np.random.randn())
        data.append([Pp, Po, r_noisy])
    return np.array(data)

data = _generate_example_data()
# ---- REPLACE ABOVE WITH e.g.: ----
# data = np.array([
#     [0.05, 0.10, 0.082],
#     [0.10, 0.10, 0.134],
#     ...
# ])

P_prop = data[:, 0]
P_O2   = data[:, 1]
r_obs  = data[:, 2]

# =============================================================================
# SECTION 2: MODEL DEFINITIONS
# =============================================================================

def rate_LH_molecular(P, k, K_prop, K_O):
    """
    Standard dual-site LH, molecular O2 adsorption.
    P = (P_prop, P_O2) tuple
    """
    Pp, Po = P
    num = k * K_prop * Pp * K_O * Po
    den = (1 + K_prop * Pp + K_O * Po) ** 2
    return num / den


def rate_LH_dissociative(P, k, K_prop, K_O):
    """
    LH with dissociative O2 adsorption: O2 + 2* -> 2O*
    theta_O = sqrt(K_O * P_O2) * theta_*
    More physically realistic for Pd.
    """
    Pp, Po = P
    sqrt_Ko_Po = np.sqrt(K_O * Po)
    num = k * K_prop * Pp * sqrt_Ko_Po
    den = (1 + K_prop * Pp + sqrt_Ko_Po) ** 2
    return num / den


def rate_power_law(P, k, a, b):
    """
    Simple power law: r = k * P_prop^a * P_O2^b
    Useful as a sanity check and for extracting apparent reaction orders.
    """
    Pp, Po = P
    return k * (Pp ** a) * (Po ** b)

# =============================================================================
# SECTION 3: FITTING
# Uses differential evolution (global optimizer) first, then curve_fit to
# refine and get parameter covariance / standard errors.
# This avoids getting stuck in local minima — important for LH models.
# =============================================================================

def fit_model(rate_func, P_data, r_data, bounds, model_name="Model"):
    """
    Two-stage fit:
      1. Differential evolution (global, no initial guess needed)
      2. curve_fit (local, gives covariance -> standard errors)
    """
    print(f"\n{'='*55}")
    print(f"  Fitting: {model_name}")
    print(f"{'='*55}")

    P_tuple = (P_data[:, 0], P_data[:, 1])

    # Stage 1: global search
    def residuals_sq(params):
        r_pred = rate_func(P_tuple, *params)
        return np.sum((r_pred - r_data) ** 2)

    result = differential_evolution(residuals_sq, bounds, seed=42,
                                    maxiter=5000, tol=1e-10, polish=True)
    p0 = result.x

    # Stage 2: refine with curve_fit to get covariance
    try:
        lb = [b[0] for b in bounds]
        ub = [b[1] for b in bounds]
        # Clip DE result to bounds before passing as p0 to curve_fit
        p0_clipped = np.clip(p0, [b*(1+1e-8) for b in lb], [b*(1-1e-8) for b in ub])
        popt, pcov = curve_fit(rate_func, P_tuple, r_data,
                               p0=p0_clipped, bounds=(lb, ub),
                               maxfev=50000)
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError:
        print("  curve_fit did not converge; using DE result (no std errors)")
        popt = p0
        perr = [float('nan')] * len(p0)

    r_pred = rate_func(P_tuple, *popt)
    ss_res = np.sum((r_pred - r_data) ** 2)
    ss_tot = np.sum((r_data - np.mean(r_data)) ** 2)
    r2 = 1 - ss_res / ss_tot
    rmse = np.sqrt(np.mean((r_pred - r_data) ** 2))

    return popt, perr, r_pred, r2, rmse


# --- Fit all three models ---

# Parameter bounds: (lower, upper) for each param — adjust if fit is struggling
bounds_LH   = [(1e-4, 1e3), (1e-4, 1e3), (1e-4, 1e3)]   # k, K_prop, K_O
bounds_pl   = [(1e-4, 1e2), (0.0, 2.0), (0.0, 2.0)]       # k, a, b

results = {}

popt, perr, r_pred, r2, rmse = fit_model(
    rate_LH_molecular, data, r_obs, bounds_LH, "LH Molecular O2")
results["LH Molecular"] = dict(
    popt=popt, perr=perr, r_pred=r_pred, r2=r2, rmse=rmse,
    labels=["k", "K_prop", "K_O"], func=rate_LH_molecular)

popt, perr, r_pred, r2, rmse = fit_model(
    rate_LH_dissociative, data, r_obs, bounds_LH, "LH Dissociative O2")
results["LH Dissociative"] = dict(
    popt=popt, perr=perr, r_pred=r_pred, r2=r2, rmse=rmse,
    labels=["k", "K_prop", "K_O"], func=rate_LH_dissociative)

popt, perr, r_pred, r2, rmse = fit_model(
    rate_power_law, data, r_obs, bounds_pl, "Power Law")
results["Power Law"] = dict(
    popt=popt, perr=perr, r_pred=r_pred, r2=r2, rmse=rmse,
    labels=["k", "alpha (prop order)", "beta (O2 order)"], func=rate_power_law)

# =============================================================================
# SECTION 4: PRINT RESULTS TABLE
# =============================================================================

print("\n\n" + "="*55)
print("  SUMMARY OF FIT RESULTS")
print("="*55)

for name, res in results.items():
    print(f"\n  {name}")
    print(f"  R² = {res['r2']:.4f}    RMSE = {res['rmse']:.4f}")
    for label, val, err in zip(res['labels'], res['popt'], res['perr']):
        print(f"    {label:30s} = {val:10.4f}  ±  {err:.4f}")

print("\n  Model discrimination rule of thumb:")
print("  R² > 0.99 and physically sensible K values (>0) = acceptable fit")
print("  If multiple models fit equally well, prefer mechanistically motivated one")

# =============================================================================
# SECTION 5: PLOTS
# =============================================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle("LH Kinetic Model Fitting — Propylene Oxidation / Pd/C", fontsize=13)

colors = {"LH Molecular": "#2196F3", "LH Dissociative": "#4CAF50", "Power Law": "#FF9800"}

# Plot 1: Parity plot (observed vs predicted) for all models
ax = axes[0]
for name, res in results.items():
    ax.scatter(r_obs, res['r_pred'], label=f"{name} (R²={res['r2']:.3f})",
               alpha=0.7, color=colors[name], s=50)
lims = [min(r_obs)*0.9, max(r_obs)*1.1]
ax.plot(lims, lims, 'k--', lw=1, label="Perfect fit")
ax.set_xlabel("Observed rate")
ax.set_ylabel("Predicted rate")
ax.set_title("Parity Plot")
ax.legend(fontsize=7)
ax.set_xlim(lims); ax.set_ylim(lims)

# Plot 2: Rate vs P_prop at fixed P_O2 (median value)
ax = axes[1]
Po_fixed = np.median(P_O2)
Pp_range = np.linspace(P_prop.min(), P_prop.max(), 100)
P_sweep = (Pp_range, np.full_like(Pp_range, Po_fixed))

# Actual data points at ~this P_O2
mask = np.abs(P_O2 - Po_fixed) < 0.05
ax.scatter(P_prop[mask], r_obs[mask], color='black', zorder=5,
           label="Data", s=60, marker='o')

for name, res in results.items():
    r_sweep = res['func'](P_sweep, *res['popt'])
    ax.plot(Pp_range, r_sweep, color=colors[name], label=name, lw=2)

ax.set_xlabel("P_prop (bar)")
ax.set_ylabel("Rate")
ax.set_title(f"Rate vs P_prop\n(P_O2 ≈ {Po_fixed:.2f} bar fixed)")
ax.legend(fontsize=7)

# Plot 3: Rate vs P_O2 at fixed P_prop (median value)
ax = axes[2]
Pp_fixed = np.median(P_prop)
Po_range = np.linspace(P_O2.min(), P_O2.max(), 100)
P_sweep2 = (np.full_like(Po_range, Pp_fixed), Po_range)

mask2 = np.abs(P_prop - Pp_fixed) < 0.05
ax.scatter(P_O2[mask2], r_obs[mask2], color='black', zorder=5,
           label="Data", s=60, marker='o')

for name, res in results.items():
    r_sweep2 = res['func'](P_sweep2, *res['popt'])
    ax.plot(Po_range, r_sweep2, color=colors[name], label=name, lw=2)

ax.set_xlabel("P_O2 (bar)")
ax.set_ylabel("Rate")
ax.set_title(f"Rate vs P_O2\n(P_prop ≈ {Pp_fixed:.2f} bar fixed)")
ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("lh_kinetics_fit.png", dpi=150, bbox_inches='tight')
plt.show()
print("\nPlot saved to lh_kinetics_fit.png")

# =============================================================================
# SECTION 6: SENSITIVITY / RATE SURFACE (useful for catalyst design intuition)
# =============================================================================

fig2, ax2 = plt.subplots(figsize=(7, 5))

# Use best-fitting LH model
best_model_name = max(results, key=lambda n: results[n]['r2'])
best = results[best_model_name]
print(f"\nBest model by R²: {best_model_name}")

Pp_grid = np.linspace(0.01, 0.5, 80)
Po_grid = np.linspace(0.01, 0.5, 80)
PP, PO = np.meshgrid(Pp_grid, Po_grid)
R_surf = best['func']((PP, PO), *best['popt'])

cf = ax2.contourf(PP, PO, R_surf, levels=20, cmap='viridis')
plt.colorbar(cf, ax=ax2, label="Rate")
ax2.scatter(P_prop, P_O2, c='red', s=40, zorder=5, label="Data points")
ax2.set_xlabel("P_prop (bar)")
ax2.set_ylabel("P_O2 (bar)")
ax2.set_title(f"Rate Surface — {best_model_name}\n(useful for finding optimal operating P)")
ax2.legend()
plt.tight_layout()
plt.savefig("lh_rate_surface.png", dpi=150, bbox_inches='tight')
plt.show()
print("Rate surface saved to lh_rate_surface.png")