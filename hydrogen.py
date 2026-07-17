"""
Analytical (Gaussian-expansion) model of the neutral-hydrogen X-ray atomic
scattering factor, following the method of:

    Thorkildsen, G. (2023). "New benchmarks in the modelling of X-ray atomic
    form factors." Acta Cryst. A79, 318-330. https://doi.org/10.1107/S2053273323003996

The model (paper's Eq. 5), based on the inverse Mott-Bethe formula:

    f_X(s; Z) = Z - 8*pi^2*a0*s^2 * [ alpha + sum_i c_i * exp(-d_i * s^2) ]

is fitted by nonlinear least squares to reference f0(s) data over
s = sin(theta)/lambda in [0, 6] Angstrom^-1, step 0.01 (601 points), following
the paper's Search -> Expand -> Verify strategy: start with a small number of
Gaussians n, refine from many random starts, and grow n until the fit stops
improving.

For neutral hydrogen, the reference data is the EXACT non-relativistic
quantum-mechanical result (unlike heavier atoms, which require numerical
Hartree-Fock densities and thus tabulated f0 data):

    f0(s) = 1 / [1 + (2*pi*a0*s)^2]^2

derived from the Fourier transform of the 1s ground-state electron density,
rho(r) = (1/(pi*a0^3)) * exp(-2r/a0).
"""

import numpy as np
from scipy.optimize import curve_fit

# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------
A0 = 0.5291772109  # Bohr radius, Angstrom (CODATA)
Z = 1               # neutral hydrogen: 1 electron

# ---------------------------------------------------------------------
# s grid: 0.00 to 6.00 Angstrom^-1, step 0.01  (601 points)
# ---------------------------------------------------------------------
s = np.round(np.arange(0, 601) * 0.01, 2)

# ---------------------------------------------------------------------
# Reference ("exact") neutral-hydrogen form factor
# ---------------------------------------------------------------------
def f0_exact(s, a0=A0):
    return 1.0 / (1.0 + (2 * np.pi * a0 * s) ** 2) ** 2

f0_data = f0_exact(s)


# ---------------------------------------------------------------------
# Model: inverse Mott-Bethe, MB[nG + alpha]  (paper's Eq. 5)
# ---------------------------------------------------------------------
def mb_model(s, params, n, Z=Z, a0=A0):
    alpha = params[0]
    c = np.asarray(params[1:1 + n])
    d = np.asarray(params[1 + n:1 + 2 * n])
    gaussians = np.exp(-np.outer(s ** 2, d))       # shape (len(s), n)
    bracket = alpha + gaussians @ c
    return Z - 8 * np.pi ** 2 * a0 * s ** 2 * bracket


def make_curvefit_func(n):
    """Wrap mb_model for scipy.optimize.curve_fit's flat-parameter signature."""
    def func(s, *params):
        return mb_model(s, params, n)
    return func


# ---------------------------------------------------------------------
# Fitting procedure: Search -> Expand -> Verify
# ---------------------------------------------------------------------
def random_start(n, rng):
    """Random initial guesses for {alpha, c_1..c_n, d_1..d_n}, in the spirit
    of the paper's RandomReal[...] seeding, spread over increasing d ranges."""
    alpha0 = rng.uniform(0.001, 0.010)
    c0 = rng.uniform(0.05, 1.0, size=n)
    # spread d over widening decades so Search can find well separated terms
    edges = np.geomspace(0.02, 60.0, n + 1)
    d0 = np.array([rng.uniform(edges[i], edges[i + 1]) for i in range(n)])
    return np.concatenate(([alpha0], c0, d0))


def fit_n_gaussians(n, n_tries=40, seed=0):
    """Search: many random-start refinements; keep the best converged fit
    satisfying the paper's acceptance conditions (c_i > 0, d_i in [0.01,1000],
    d's not too close together, alpha >= 0)."""
    rng = np.random.default_rng(seed)
    func = make_curvefit_func(n)

    lower = [0.0] + [0.0] * n + [0.01] * n
    upper = [np.inf] + [np.inf] * n + [1000.0] * n

    best = None
    for _ in range(n_tries):
        p0 = random_start(n, rng)
        try:
            popt, _ = curve_fit(func, s, f0_data, p0=p0,
                                 bounds=(lower, upper), maxfev=20000)
        except RuntimeError:
            continue

        alpha, c, d = popt[0], popt[1:1 + n], popt[1 + n:1 + 2 * n]
        if alpha < 0 or np.any(c <= 0):
            continue
        order = np.argsort(d)
        d_sorted = d[order]
        if n > 1 and np.min(d_sorted[1:] / d_sorted[:-1]) <= 1.5:
            continue  # d's too close -> reject (paper's condition)

        resid = f0_data - mb_model(s, popt, n)
        mae = np.mean(np.abs(resid))
        if best is None or mae < best[0]:
            c_sorted, d_sorted = c[order], d[order]
            best = (mae, alpha, c_sorted, d_sorted, popt)

    return best  # (mae, alpha, c, d, raw_popt) or None


def run_search_and_expand(n_max=14, mae_tol=1e-7):
    """Expand: grow n = 1, 2, 3, ... until MAE stops improving meaningfully
    or reaches the target tolerance."""
    results = {}
    prev_mae = np.inf
    for n in range(1, n_max + 1):
        result = fit_n_gaussians(n)
        if result is None:
            print(f"  n={n}: no converged/accepted fit found")
            continue
        mae, alpha, c, d, popt = result
        results[n] = result
        print(f"  n={n} Gaussians: MAE = {mae:.3e}  alpha = {alpha:.6e}")
        if mae < mae_tol or mae > 0.95 * prev_mae:
            prev_mae = mae
            break
        prev_mae = mae
    return results


# ---------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------
if __name__ == "__main__":
    print("Fitting MB[nG+alpha] model to exact neutral-hydrogen f0(s), "
          "s in [0,6] A^-1, step 0.01\n")
    print("Search & Expand:")
    results = run_search_and_expand(n_max=20, mae_tol=1e-9)

    n_best = max(results.keys())
    mae, alpha, c, d, popt = results[n_best]

    print(f"\nFinal accepted model: MB[{n_best}G + alpha]")
    print(f"  alpha = {alpha:.8e}")
    print(f"  {'i':>2}  {'c_i':>16}  {'d_i (A^2)':>16}")
    for i, (ci, di) in enumerate(zip(c, d), start=1):
        print(f"  {i:>2}  {ci:>16.8e}  {di:>16.8e}")

    f_model = mb_model(s, popt, n_best)
    resid = f0_data - f_model
    print(f"\nFit quality over s = [0.00, 6.00] A^-1, step 0.01 (601 pts):")
    print(f"  Mean absolute error  <|df0|>_s = {np.mean(np.abs(resid)):.4e}")
    print(f"  Max absolute error   |df0|_max = {np.max(np.abs(resid)):.4e}")
    print(f"  RMS error                       = {np.sqrt(np.mean(resid**2)):.4e}")

    # Save full table: s, exact f0, model f0, difference
    out = np.column_stack([s, f0_data, f_model, resid])
    header = "s(1/A)      f0_exact        f0_model        f0_exact-f0_model"
    np.savetxt("hydrogen_form_factor_fit.csv", out,
               header=header, fmt="%.6f  %.10e  %.10e  %.3e", comments="")
    print("\nSaved full s-grid comparison table to hydrogen_form_factor_fit.csv")

    # Save coefficients to their own small file too
    with open("hydrogen_MB_coefficients.txt", "w") as fh:
        fh.write("Neutral hydrogen (Z=1) inverse Mott-Bethe MB[nG+alpha] model\n")
        fh.write("f_X(s) = Z - 8*pi^2*a0*s^2*[alpha + sum_i c_i*exp(-d_i*s^2)]\n")
        fh.write(f"a0 = {A0} Angstrom,  s range [0,6] A^-1, step 0.01\n\n")
        fh.write(f"n = {n_best}\n")
        fh.write(f"alpha = {alpha:.12e}\n")
        for i, (ci, di) in enumerate(zip(c, d), start=1):
            fh.write(f"c_{i} = {ci:.12e}   d_{i} = {di:.12e} A^2\n")
        fh.write(f"\nMean absolute error = {np.mean(np.abs(resid)):.4e}\n")
        fh.write(f"Max absolute error  = {np.max(np.abs(resid)):.4e}\n")
    print("Saved coefficients to hydrogen_MB_coefficients.txt")
