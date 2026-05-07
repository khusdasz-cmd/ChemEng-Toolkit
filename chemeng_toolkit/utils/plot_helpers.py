"""Shared plotting utilities and color schemes for ChemEng-Toolkit."""

import numpy as np
import matplotlib.pyplot as plt


# Color scheme used throughout the toolkit
COLORS = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"]

# Standard sample labels for fractionated lignin
SAMPLE_LABELS = [">100k", "100k-20k", "20k-10k", "10k-5k"]

# FT-IR peak assignments: (wavenumber, label)
FTIR_PEAKS = [
    (3400, "O-H"),
    (2947, "C-H"),
    (1682, "C=O"),
    (1600, "Aromatic ring"),
    (1512, "Aromatic ring"),
    (1404, "C-H/OCH3"),
    (1327, "S ring"),
    (1265, "G ring"),
    (1126, "S-H"),
    (1034, "C-O"),
]


def setup_plot_style(usetex: bool = False):
    """Configure matplotlib for publication-quality plots.

    Parameters
    ----------
    usetex : bool
        Whether to use LaTeX for text rendering (requires LaTeX installed).
    """
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 14,
        "legend.fontsize": 10,
        "lines.linewidth": 1.8,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "axes.unicode_minus": False,
    })
    if usetex:
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Times New Roman"],
        })


def downsample_data(temp, data, target_pts=2000):
    """Downsample data to a uniform temperature grid via interpolation.

    Parameters
    ----------
    temp : np.ndarray
        Temperature array.
    data : np.ndarray
        Data array to downsample.
    target_pts : int
        Number of points in the output grid.

    Returns
    -------
    t_grid : np.ndarray
        Uniformly spaced temperature grid.
    data_interp : np.ndarray
        Interpolated data on the grid.
    """
    t_min, t_max = temp.min(), temp.max()
    t_grid = np.linspace(t_min, t_max, target_pts)
    data_interp = np.interp(t_grid, temp, data)
    return t_grid, data_interp
