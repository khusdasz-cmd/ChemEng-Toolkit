"""Pyrolysis kinetics analysis using the Coats-Redfern method.

Supports single-heating-rate analysis with 15 reaction mechanism models,
segmented (alpha-dependent) activation energy calculation, and Ea vs. alpha
visualization.

Reaction models available:
    - F1, F2, F3: Reaction order models
    - R2, R3: Phase boundary-controlled
    - D1-D4: Diffusion models
    - A2-A4: Avrami-Erofeev nucleation
    - P2-P4: Power law models
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

from chemeng_toolkit.utils.plot_helpers import COLORS, SAMPLE_LABELS, setup_plot_style

R = 8.314  # Universal gas constant, J/(mol.K)

# Default heating rate (°C/min) — adjust if different
BETA = 10.0


# ---------------------------------------------------------------------------
# Reaction mechanism models g(alpha)
# ---------------------------------------------------------------------------

REACTION_MODELS = [
    # Reaction order models
    ("F1", "First-order", lambda a: -np.log(1 - a)),
    ("F2", "Second-order", lambda a: 1 / (1 - a) - 1),
    ("F3", "Third-order", lambda a: ((1 / (1 - a)) ** 2 - 1) / 2),
    # Phase boundary
    ("R2", "Contracting cylinder", lambda a: 1 - (1 - a) ** (1 / 2)),
    ("R3", "Contracting sphere", lambda a: 1 - (1 - a) ** (1 / 3)),
    # Diffusion
    ("D1", "1D diffusion", lambda a: a ** 2),
    ("D2", "2D diffusion", lambda a: (1 - a) * np.log(1 - a) + a),
    ("D3", "3D diffusion (Jander)", lambda a: (1 - (1 - a) ** (1 / 3)) ** 2),
    ("D4", "3D diffusion (Ginstling)", lambda a: 1 - 2 * a / 3 - (1 - a) ** (2 / 3)),
    # Avrami-Erofeev
    ("A2", "Avrami-Erofeev (n=2)", lambda a: (-np.log(1 - a)) ** (1 / 2)),
    ("A3", "Avrami-Erofeev (n=3)", lambda a: (-np.log(1 - a)) ** (1 / 3)),
    ("A4", "Avrami-Erofeev (n=4)", lambda a: (-np.log(1 - a)) ** (1 / 4)),
    # Power law
    ("P2", "Power law (n=2)", lambda a: a ** (1 / 2)),
    ("P3", "Power law (n=3)", lambda a: a ** (1 / 3)),
    ("P4", "Power law (n=4)", lambda a: a ** (1 / 4)),
]


# ---------------------------------------------------------------------------
# Coats-Redfern (global)
# ---------------------------------------------------------------------------


def coats_redfern(temp, weight, t_range=(150, 550), alpha_range=(0.05, 0.90),
                  heating_rate=BETA):
    """Perform Coats-Redfern kinetic analysis on a single sample.

    Parameters
    ----------
    temp : np.ndarray
        Temperature array (°C).
    weight : np.ndarray
        Weight array (%).
    t_range : tuple of (float, float)
        Temperature range for analysis (default: 150-550°C).
    alpha_range : tuple of (float, float)
        Conversion range (default: 0.05-0.90).
    heating_rate : float
        Heating rate in °C/min (default: 10.0).

    Returns
    -------
    list of dict or None
        Each entry: code, name, E_kJ, lnA, R2, n_points.
        Sorted by R² descending. Returns None if analysis fails.
    """
    w_norm = weight / weight[0] * 100

    mask = (temp >= t_range[0]) & (temp <= t_range[1])
    if not mask.any():
        return None

    T_K = temp[mask] + 273.15
    w = w_norm[mask]

    w0 = w[0]
    w_inf = w[-1]
    alpha = (w0 - w) / (w0 - w_inf)

    a_mask = (alpha >= alpha_range[0]) & (alpha <= alpha_range[1])
    if a_mask.sum() < 5:
        return None

    T_fit = T_K[a_mask]
    alpha_fit = alpha[a_mask]

    results = []
    for code, name, g_func in REACTION_MODELS:
        try:
            g_val = g_func(alpha_fit)
            valid = np.isfinite(g_val) & (g_val > 0)
            if valid.sum() < 5:
                continue

            x = 1 / T_fit[valid]
            y = np.log(g_val[valid] / T_fit[valid] ** 2)

            slope, intercept, r_val, _, _ = linregress(x, y)
            r2 = r_val ** 2

            E = -slope * R / 1000  # kJ/mol
            A = np.exp(intercept) * heating_rate * E * 1000 / R

            results.append({
                "code": code,
                "name": name,
                "E_kJ": E,
                "lnA": np.log(A) if A > 0 else np.nan,
                "R2": r2,
                "n_points": int(valid.sum()),
            })
        except Exception:
            continue

    if not results:
        return None

    results.sort(key=lambda x: x["R2"], reverse=True)
    return results


def run_coats_redfern_batch(data_dict, t_range=(150, 550), alpha_range=(0.05, 0.90)):
    """Run Coats-Redfern analysis on multiple samples.

    Parameters
    ----------
    data_dict : dict of {str: dict}
        Sample label -> {"temp": ndarray, "weight": ndarray, ...}
    t_range, alpha_range : tuples
        As in coats_redfern().

    Returns
    -------
    best_models : dict of {str: dict}
        Sample label -> best model result.
    f1_results : dict of {str: dict}
        Sample label -> F1 model result.
    """
    best_models = {}
    f1_results = {}

    print("=" * 80)
    print("Coats-Redfern Kinetic Analysis")
    print("=" * 80)
    print(f"Heating rate: {BETA}°C/min, Temperature: {t_range[0]}-{t_range[1]}°C")
    print(f"Conversion: alpha = {alpha_range[0]}-{alpha_range[1]}")
    print("=" * 80)

    for label in data_dict:
        d = data_dict[label]
        results = coats_redfern(d["temp"], d["weight"], t_range, alpha_range)
        if results is None:
            print(f"\n  {label}: analysis failed")
            continue

        print(f"\n  ----- {label} -----")
        print(f"  {'Model':<8} {'E (kJ/mol)':<14} {'ln(A/s-1)':<14} {'R2':<8}")
        print(f"  {'-' * 48}")
        for r in results[:5]:  # top 5
            marker = " *" if r is results[0] else ""
            print(f"  {r['code']:<6s} {r['E_kJ']:<12.2f} {r['lnA']:<12.2f} "
                  f"{r['R2']:<8.4f}{marker}")

        best_models[label] = results[0]
        for r in results:
            if r["code"] == "F1":
                f1_results[label] = r
                break

    # Summary table
    print("\n" + "=" * 80)
    print("Best-fit Model Summary")
    print("=" * 80)
    print(f"  {'Sample':<10} {'Model':<8} {'E (kJ/mol)':<14} {'ln(A/s-1)':<14} {'R2':<8}")
    print(f"  {'-' * 56}")
    for label, best in best_models.items():
        print(f"  {label:<10s} {best['code']:<6s} {best['E_kJ']:<12.2f} "
              f"{best['lnA']:<12.2f} {best['R2']:<8.4f}")

    if f1_results:
        print("\n" + "=" * 80)
        print("Activation Energy via First-Order Model (F1)")
        print("=" * 80)
        print(f"  {'Sample':<10} {'E (kJ/mol)':<14} {'ln(A/s-1)':<14} {'R2':<8}")
        for label, r in f1_results.items():
            print(f"  {label:<10s} {r['E_kJ']:<12.2f} {r['lnA']:<12.2f} {r['R2']:<8.4f}")

    return best_models, f1_results


# ---------------------------------------------------------------------------
# Segmented Coats-Redfern (alpha-dependent Ea)
# ---------------------------------------------------------------------------


def segmented_coats_redfern(temp, weight, alpha_values=None, half_window=0.05,
                            t_range=(200, 500)):
    """Calculate alpha-dependent activation energy using segmented F1 model.

    Useful for constructing Ea vs. alpha plots to study reaction mechanisms.

    Parameters
    ----------
    temp : np.ndarray
        Temperature array (°C).
    weight : np.ndarray
        Weight array (%).
    alpha_values : list of float, optional
        Alpha values to evaluate (default: 0.1 to 0.8 step 0.1).
    half_window : float
        Half-width of the alpha window (default: 0.05).
    t_range : tuple of (float, float)
        Temperature range for analysis.

    Returns
    -------
    list of dict
        Each entry: {"alpha": float, "E_kJ": float, "lnA": float, "R2": float}
    """
    if alpha_values is None:
        alpha_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    w_norm = weight / weight[0] * 100

    mask = (temp >= t_range[0]) & (temp <= t_range[1])
    T_K = temp[mask] + 273.15
    w = w_norm[mask]

    w0 = w[0]
    w_inf = w[-1]
    alpha = (w0 - w) / (w0 - w_inf)

    results = []
    for av in alpha_values:
        lo = max(0.02, av - half_window)
        hi = min(0.95, av + half_window)
        am = (alpha >= lo) & (alpha <= hi)
        if am.sum() < 5:
            results.append({"alpha": av, "E_kJ": None, "lnA": None, "R2": None})
            continue

        T_fit = T_K[am]
        a_fit = alpha[am]
        g_val = -np.log(1 - a_fit)

        valid = (g_val > 0) & np.isfinite(g_val)
        if valid.sum() < 5:
            results.append({"alpha": av, "E_kJ": None, "lnA": None, "R2": None})
            continue

        x = 1 / T_fit[valid]
        y = np.log(g_val[valid] / T_fit[valid] ** 2)

        slope, intercept, r_val, _, _ = linregress(x, y)
        E = -slope * R / 1000
        A = np.exp(intercept) * BETA * E * 1000 / R

        results.append({
            "alpha": av,
            "E_kJ": E,
            "lnA": np.log(A) if A > 0 else None,
            "R2": r_val ** 2,
        })

    return results


def run_segmented_batch(data_dict, alpha_values=None, half_window=0.05,
                        t_range=(200, 500)):
    """Run segmented Coats-Redfern on multiple samples and print results table.

    Parameters
    ----------
    data_dict : dict of {str: dict}
        Sample label -> data.
    alpha_values, half_window, t_range : as in segmented_coats_redfern().

    Returns
    -------
    dict of {str: list of dict}
        Sample label -> segmented results.
    """
    if alpha_values is None:
        alpha_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

    all_results = {}

    print("\nSegmented Coats-Redfern (F1 model) — alpha-dependent Ea")
    print("=" * 70)

    for label in data_dict:
        d = data_dict[label]
        results = segmented_coats_redfern(
            d["temp"], d["weight"], alpha_values, half_window, t_range
        )
        all_results[label] = results

    # Print results table
    header = f"{'alpha':<8}"
    for label in data_dict:
        header += f"  {label:<14}"
    print(header)
    print("-" * 70)

    for av in alpha_values:
        row = f"{av:<8.1f}"
        for label in data_dict:
            for r in all_results[label]:
                if r["alpha"] == av and r["E_kJ"] is not None:
                    row += f"  {r['E_kJ']:<14.2f}"
                    break
            else:
                row += f"  {'N/A':<14}"
        print(row)

    print("-" * 70)
    avg_row = f"{'Mean':<8}"
    for label in data_dict:
        vals = [r["E_kJ"] for r in all_results[label]
                if r["E_kJ"] is not None and r["E_kJ"] > 0]
        if vals:
            avg_row += f"  {np.mean(vals):<14.2f}"
        else:
            avg_row += f"  {'N/A':<14}"
    print(avg_row)

    return all_results


# ---------------------------------------------------------------------------
# Ea vs. alpha plotting
# ---------------------------------------------------------------------------


def plot_ea_vs_alpha(ea_data, labels=None, output_path=None, show=True):
    """Plot activation energy vs. conversion rate (Ea vs. alpha).

    Parameters
    ----------
    ea_data : dict of {str: list}
        Sample label -> list of Ea values (at standard alpha points).
    labels : list of str, optional
        Sample labels for the legend (default: sorted keys of ea_data).
    output_path : str or Path, optional
        If provided, save figure to this path.
    show : bool
        Whether to display the plot.
    """
    if labels is None:
        labels = sorted(ea_data.keys())

    alpha = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    markers = ["o", "s", "D", "^"]

    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8, 5.5))

    for idx, label in enumerate(labels):
        vals = ea_data[label]
        color = COLORS[idx % len(COLORS)]
        marker = markers[idx % len(markers)]
        valid = [(a, v) for a, v in zip(alpha, vals)
                 if v is not None and not np.isnan(v)]
        if not valid:
            continue
        valid_alpha, valid_vals = zip(*valid)
        ax.plot(valid_alpha, valid_vals, color=color, linewidth=1.8,
                marker=marker, markersize=7, markerfacecolor="white",
                markeredgewidth=1.5, label=label, zorder=5)
        for a, v in zip(valid_alpha, valid_vals):
            ax.annotate(f"{v:.1f}", xy=(a, v), xytext=(5, 8),
                        textcoords="offset points", fontsize=7.5,
                        color=color, ha="center", va="bottom")

    ax.set_xlabel("Conversion (alpha)")
    ax.set_ylabel("Activation Energy (kJ/mol)")
    ax.set_xlim(0.05, 0.85)
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, linestyle=":", alpha=0.3, linewidth=0.5)

    fig.tight_layout()

    if output_path:
        from pathlib import Path
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    if show:
        plt.show()
    plt.close(fig)
