"""FT-IR spectroscopy data loading, processing, and visualization.

Supports loading CSV-format FT-IR data, batch processing of fractionated
samples, and stacked spectral plotting with peak annotation.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from chemeng_toolkit.utils.plot_helpers import (
    setup_plot_style,
)


# Standard FT-IR peak assignments: (wavenumber, label)
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


class FTIRProcessor:
    """Load, process, and visualize FT-IR spectral data.

    Parameters
    ----------
    files : list of (str, str, str), optional
        List of (filename, label, hex_color) tuples.
        Example: [("sample.CSV", "L5", "#1f77b4"), ...]
    data_dir : str or Path, optional
        Base directory for data files.
    """

    def __init__(self, files=None, data_dir=None):
        self.files = files or []
        self.data_dir = Path(data_dir) if data_dir else None
        self.data = []  # list of (label, (wn, abs), color)

    # ----------------------------------------------------------------
    # Data loading
    # ----------------------------------------------------------------

    def load_csv(self, filepath):
        """Load a single FT-IR CSV file.

        Expects two columns: wavenumber, absorbance (no header).

        Parameters
        ----------
        filepath : str or Path
            Path to the CSV file.

        Returns
        -------
        (wn, abs) : tuple of np.ndarray
            Wavenumber and absorbance arrays.
        """
        filepath = Path(filepath)
        df = pd.read_csv(filepath, header=None, names=["wn", "abs"])
        return df["wn"].values.astype(float), df["abs"].values.astype(float)

    def load_batch(self, file_list=None):
        """Load all FT-IR files configured for this processor.

        Parameters
        ----------
        file_list : list of (str, str, str), optional
            Override self.files: (filename, label, hex_color).

        Returns
        -------
        list of (label, (wn, abs), color)
        """
        if file_list is not None:
            self.files = file_list

        self.data = []
        for filename, label, color in self.files:
            if self.data_dir:
                filepath = self.data_dir / filename
            else:
                filepath = Path(filename)

            if not filepath.exists():
                print(f"  Warning: {filepath} not found, skipping")
                continue

            wn, abs_data = self.load_csv(filepath)
            self.data.append((label, (wn, abs_data), color))
            print(f"  Loaded: {filename} ({len(wn)} points)")

        if len(self.data) >= 2:
            # Validate wavenumber grid consistency
            wn_ref = self.data[0][1][0]
            for label, (wn, _), _ in self.data[1:]:
                if not np.allclose(wn, wn_ref):
                    print(f"  Warning: {label} wavenumber grid mismatch!")

        return self.data

    # ----------------------------------------------------------------
    # Plotting
    # ----------------------------------------------------------------

    def plot_stacked(self, normalize=False, offset_step=None,
                     title=None, output_path=None, show=True):
        """Plot FT-IR spectra with vertical offset stacking.

        Each spectrum is vertically offset for clear visual comparison.
        Major absorption peaks are annotated with assignments.

        Parameters
        ----------
        normalize : bool
            If True, normalize each spectrum to its maximum absorbance.
        offset_step : float, optional
            Vertical offset between spectra. Auto-calculated if None.
        title : str, optional
            Plot title.
        output_path : str or Path, optional
            If provided, save figure to this path.
        show : bool
            Whether to display the plot.
        """
        if not self.data:
            print("  No data loaded. Call load_batch() first.")
            return

        if offset_step is None:
            offset_step = 1.3 if normalize else 0.8

        setup_plot_style()
        fig, ax = plt.subplots(figsize=(12, 6.5))

        data_list = self.data
        n_curves = len(data_list)

        # Plot spectra with vertical offset
        for i, (label, (wn, abs_data), color) in enumerate(data_list):
            y = abs_data / abs_data.max() if normalize else abs_data
            offset = i * offset_step
            ax.plot(wn, y + offset, color=color, label=label, linewidth=1.2)

        # Zero-point lines
        for i, (label, _, color) in enumerate(data_list):
            offset = i * offset_step
            ax.axhline(y=offset, color=color, linestyle=":", alpha=0.3, linewidth=0.6)
            ax.text(wn[0] - 30, offset, f"  {label}", fontsize=9, color=color,
                    ha="right", va="center", fontweight="bold")

        # Peak annotation lines
        y_span = (n_curves - 1) * offset_step + (
            1.0 if normalize else max(d[1][1].max() for d in data_list)
        )
        for pos, _ in FTIR_PEAKS:
            ax.axvline(x=pos, color="gray", linestyle="--", alpha=0.35, linewidth=0.7)

        # IR convention: decreasing wavenumber
        ax.invert_xaxis()

        # Peak labels at top
        y_min, y_max = 0, y_span + offset_step * 0.15
        ax.set_ylim(0, y_max)
        label_y = y_max - (y_max - 0) * 0.04
        for pos, label_text in FTIR_PEAKS:
            ax.text(pos, label_y, f"{label_text}\n{pos}", fontsize=7, color="gray",
                    ha="center", va="top", linespacing=1.2)

        # Hide Y axis (no quantitative meaning after offset)
        ax.set_yticks([])
        ax.set_yticklabels([])

        ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
        ax.set_ylabel("Absorbance (offset)")
        if title:
            ax.set_title(title)

        ax.legend(loc="upper right", framealpha=0.9, edgecolor="gray",
                  title="Fraction")
        ax.grid(True, axis="x", linestyle=":", alpha=0.5)
        ax.grid(False, axis="y")

        fig.tight_layout()
        self._save_or_show(fig, output_path, show)

    def plot_individual(self, normalize=False, output_path=None, show=True):
        """Plot individual spectra in separate subplots."""
        n = len(self.data)
        if n == 0:
            print("  No data loaded.")
            return

        setup_plot_style()
        fig, axes = plt.subplots(n, 1, figsize=(10, 2.5 * n), sharex=True)

        if n == 1:
            axes = [axes]

        for ax, (label, (wn, abs_data), color) in zip(axes, self.data):
            y = abs_data / abs_data.max() if normalize else abs_data
            ax.plot(wn, y, color=color, linewidth=1.2)
            ax.set_ylabel("Absorbance")
            ax.set_title(label)
            ax.invert_xaxis()
            for pos, peak_label in FTIR_PEAKS:
                ax.axvline(x=pos, color="gray", linestyle=":", alpha=0.3)
            ax.grid(True, axis="x", linestyle=":", alpha=0.5)

        axes[-1].set_xlabel(r"Wavenumber (cm$^{-1}$)")
        fig.tight_layout()
        self._save_or_show(fig, output_path, show)

    # ----------------------------------------------------------------
    # Internal
    # ----------------------------------------------------------------

    @staticmethod
    def _save_or_show(fig, output_path, show):
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"  Saved: {output_path}")
        if show:
            plt.show()
        plt.close(fig)
