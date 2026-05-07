"""TGA/DSC data loading, processing, and visualization.

Supports loading .xls files from TA Instruments, batch processing of
fractionated samples, TG/DTG plotting, decomposition stage analysis,
and feature parameter extraction (T_5%, T_10%, T_max, residue).
"""

from pathlib import Path
from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks

from chemeng_toolkit.utils.plot_helpers import (
    COLORS,
    SAMPLE_LABELS,
    setup_plot_style,
    downsample_data,
)


class TGAProcessor:
    """Load, process, and visualize TGA/DTG data from multiple samples.

    Parameters
    ----------
    samples : dict of {str: str or Path}
        Mapping of sample label -> folder path containing .xls data files.
        Example: {">100k": "data/tga_sample/gt_100k", ...}
    sheet_name : str
        Excel sheet name containing TGA data.
    data_start_row : int
        Row index (0-based) where numeric data begins.
    column_map : dict
        Mapping of {"time": int, "temp": int, "weight": int, "dtg": int}
        specifying column indices in the sheet.
    """

    def __init__(
        self,
        samples=None,
        sheet_name="Merged 1",
        data_start_row=3,
        column_map=None,
    ):
        self.samples = samples or {}
        self.sheet_name = sheet_name
        self.data_start_row = data_start_row
        self.column_map = column_map or {
            "time": 0,
            "temp": 1,
            "weight": 3,
            "dtg": 4,
        }
        self.data = {}  # {label: {"temp": ndarray, "weight": ndarray, ...}}

    # ----------------------------------------------------------------
    # Data loading
    # ----------------------------------------------------------------

    def load_from_xls(self, filepath):
        """Load a single TGA .xls file and return a data dict.

        Reads data row by row (rows 3+), parsing temperature, weight, and
        DTG columns while skipping malformed rows — same logic as the
        original TA Instruments .xls reader.

        Parameters
        ----------
        filepath : str or Path
            Path to the .xls file.

        Returns
        -------
        dict or None
            Dict with keys "time", "temp", "weight", "dtg" as numpy arrays,
            or None if loading failed.
        """
        filepath = Path(filepath)
        if not filepath.exists():
            print(f"  Warning: {filepath} not found")
            return None

        try:
            df = pd.read_excel(
                filepath,
                sheet_name=self.sheet_name,
                header=None,
                skiprows=self.data_start_row,
            )
        except Exception as e:
            print(f"  Error reading {filepath}: {e}")
            return None

        col = self.column_map
        time_arr, temp_arr, weight_pct, dtg_arr = [], [], [], []
        for r in range(len(df)):
            try:
                ti = float(df.iloc[r, col["time"]])
                t = float(df.iloc[r, col["temp"]])
                w = float(df.iloc[r, col["weight"]])
                d = float(df.iloc[r, col["dtg"]])
                time_arr.append(ti)
                temp_arr.append(t)
                weight_pct.append(w)
                dtg_arr.append(d)
            except (ValueError, TypeError):
                pass

        data = {
            "time": np.array(time_arr),
            "temp": np.array(temp_arr),
            "weight": np.array(weight_pct),
            "dtg": np.array(dtg_arr),
        }

        if len(data["temp"]) == 0:
            print(f"  Warning: no valid temperature data in {filepath}")
            return None

        return data

    def load_batch(self, folder_map):
        """Load TGA data from multiple subfolders.

        Parameters
        ----------
        folder_map : list of (str, str)
            List of (folder_name, sample_label) pairs.
            Example: [("gt_100k", ">100k"), ("100k-20k", "100k-20k"), ...]
        """
        base_dir = None
        if self.samples:
            # Use first sample's parent as base for relative folder_map paths
            first_path = Path(list(self.samples.values())[0])
            if not first_path.is_absolute():
                base_dir = Path.cwd()

        for folder, label in folder_map:
            if self.samples and label in self.samples:
                path = Path(self.samples[label])
            else:
                path = Path(folder)

            # Support glob pattern to find .xls files
            files = glob(str(path / "*.xls"))
            if not files:
                print(f"  No .xls files in {folder}")
                continue

            data = self.load_from_xls(files[0])
            if data is not None:
                self.data[label] = data
                print(f"  {label}: {len(data['temp'])} data points loaded")

        return self.data

    # ----------------------------------------------------------------
    # Analysis
    # ----------------------------------------------------------------

    def smooth(self, x, window=51):
        """Apply moving average smoothing."""
        if len(x) < window:
            window = len(x) // 2 * 2 + 1
        return uniform_filter1d(x, size=window)

    def find_peaks(self, temp, dtg):
        """Find all significant DTG peaks.

        Parameters
        ----------
        temp : np.ndarray
            Temperature array.
        dtg : np.ndarray
            DTG array (negative for mass loss).

        Returns
        -------
        list of dict
            Each entry: {"t_peak": float, "dtg_peak": float}
        """
        dtg_pos = -dtg  # make positive for peak detection
        step = max(1, len(dtg_pos) // 10000)
        t_short = temp[::step]
        d_short = dtg_pos[::step]

        d_smooth = self.smooth(d_short, window=21)
        peaks_idx, _ = find_peaks(
            d_smooth, height=0.03, prominence=0.02, width=3, distance=30
        )

        result = []
        for pi in peaks_idx:
            search_t = t_short[pi]
            search_range = (temp >= search_t - 20) & (temp <= search_t + 5)
            if search_range.any():
                raw_idx = np.where(search_range)[0]
                exact_idx = raw_idx[np.argmax(dtg_pos[search_range])]
            else:
                exact_idx = np.argmin(np.abs(temp - search_t))

            result.append({
                "t_peak": float(temp[exact_idx]),
                "dtg_peak": float(dtg_pos[exact_idx]),
            })

        result.sort(key=lambda x: x["t_peak"])
        return result

    def find_decomposition_stages(self, temp, weight, dtg):
        """Identify decomposition stages and compute characteristic parameters.

        Parameters
        ----------
        temp : np.ndarray
            Temperature array (°C).
        weight : np.ndarray
            Weight array (%).
        dtg : np.ndarray
            DTG array (%/°C).

        Returns
        -------
        stages : list of dict
            Each stage: name, t_start, t_end, t_peak, mass_loss, etc.
        summary : dict
            Overall parameters: residue, total_loss, T_5%, T_10%.
        """
        dtg_smooth = self.smooth(dtg, window=101)
        w_norm = weight / weight[0] * 100
        stages = []

        # Stage I: Dehydration (30-150°C)
        mask1 = (temp >= 30) & (temp <= 150)
        if mask1.any():
            idx1 = np.where(mask1)[0]
            peak_idx1 = idx1[np.argmin(dtg_smooth[mask1])]
            stages.append({
                "name": "I Dehydration",
                "t_start": float(temp[idx1[0]]),
                "t_end": float(temp[idx1[-1]]),
                "t_peak": float(temp[peak_idx1]),
                "dtg_peak": -dtg_smooth[peak_idx1],
                "w_start": w_norm[idx1[0]],
                "w_end": w_norm[idx1[-1]],
                "mass_loss": w_norm[idx1[0]] - w_norm[idx1[-1]],
            })

        # Stage II: Main pyrolysis (150-600°C)
        mask_main = (temp >= 150) & (temp <= 600)
        if mask_main.any():
            idx_main = np.where(mask_main)[0]
            peak_main_idx = idx_main[np.argmin(dtg_smooth[mask_main])]
            t_peak_main = temp[peak_main_idx]

            pre_mask = (temp >= 150) & (temp <= t_peak_main)
            idx_pre = np.where(pre_mask)[0]
            oc = idx_pre[dtg_smooth[pre_mask] > -0.02]
            t_start2 = temp[oc[-1]] if len(oc) > 0 else temp[idx_pre[0]]

            post_mask = (temp >= t_peak_main) & (temp <= 600)
            idx_post = np.where(post_mask)[0]
            ec = idx_post[dtg_smooth[post_mask] > -0.02]
            t_end2 = temp[ec[0]] if len(ec) > 0 else temp[idx_post[-1]]

            stages.append({
                "name": "II Main pyrolysis",
                "t_start": float(t_start2),
                "t_end": float(t_end2),
                "t_peak": float(t_peak_main),
                "dtg_peak": -dtg_smooth[peak_main_idx],
                "w_start": float(np.interp(t_start2, temp, w_norm)),
                "w_end": float(np.interp(t_end2, temp, w_norm)),
                "mass_loss": float(
                    np.interp(t_start2, temp, w_norm)
                    - np.interp(t_end2, temp, w_norm)
                ),
            })

            # Stage III: Secondary pyrolysis
            mask3 = (temp >= t_end2) & (temp <= 600)
            if mask3.any() and np.sum(mask3) > 10:
                idx3 = np.where(mask3)[0]
                min_dtg = dtg_smooth[mask3].min()
                if min_dtg < -0.03:
                    peak3_idx = idx3[np.argmin(dtg_smooth[mask3])]
                    t_peak3 = temp[peak3_idx]
                else:
                    t_peak3 = None
                stages.append({
                    "name": "III Secondary pyrolysis",
                    "t_start": float(temp[idx3[0]]),
                    "t_end": float(temp[idx3[-1]]),
                    "t_peak": float(t_peak3) if t_peak3 else None,
                    "dtg_peak": -min_dtg if min_dtg < -0.03 else 0,
                    "w_start": w_norm[idx3[0]],
                    "w_end": w_norm[idx3[-1]],
                    "mass_loss": w_norm[idx3[0]] - w_norm[idx3[-1]],
                })

        # Stage IV: Carbonization (600-800°C)
        mask4 = (temp >= 600) & (temp <= 800)
        if mask4.any():
            idx4 = np.where(mask4)[0]
            stages.append({
                "name": "IV Carbonization",
                "t_start": float(temp[idx4[0]]),
                "t_end": float(temp[idx4[-1]]),
                "t_peak": None,
                "dtg_peak": 0,
                "w_start": w_norm[idx4[0]],
                "w_end": w_norm[idx4[-1]],
                "mass_loss": w_norm[idx4[0]] - w_norm[idx4[-1]],
            })

        residue = float(w_norm[-1])
        summary = {
            "residue": residue,
            "total_loss": 100 - residue,
            "T_5%": float(np.interp(95, w_norm[::-1], temp[::-1])),
            "T_10%": float(np.interp(90, w_norm[::-1], temp[::-1])),
        }
        return stages, summary

    def summary_table(self, print_table=True):
        """Print or return TGA feature parameter summary for all samples.

        Parameters
        ----------
        print_table : bool
            If True, print the table to stdout.

        Returns
        -------
        dict of {str: dict}
            Sample label -> {T_5%, T_10%, T_max, residue, ...}
        """
        results = {}
        header = (
            f"{'Sample':<12} {'T_5% (°C)':<12} {'T_10% (°C)':<13} "
            f"{'T_max (°C)':<12} {'Residue (%)':<12}"
        )
        sep = "-" * 65

        if print_table:
            print("\n" + "=" * 65)
            print("TGA Feature Parameter Summary")
            print("=" * 65)
            print(header)
            print(sep)

        for label in self.data:
            d = self.data[label]
            temp, weight, dtg = d["temp"], d["weight"], d["dtg"]
            stages, summary = self.find_decomposition_stages(temp, weight, dtg)

            t_max = None
            for s in stages:
                if s["t_peak"] and "Main pyrolysis" in s["name"]:
                    t_max = s["t_peak"]
                    break
            if t_max is None and len(stages) > 1:
                t_max = stages[1].get("t_peak")

            results[label] = {
                "T_5%": summary["T_5%"],
                "T_10%": summary["T_10%"],
                "T_max": t_max,
                "residue": summary["residue"],
                "total_loss": summary["total_loss"],
            }

            if print_table:
                tmax_str = f"{t_max:<12.1f}" if t_max else f"{'N/A':<12}"
                print(
                    f"{label:<12} {summary['T_5%']:<12.1f} {summary['T_10%']:<13.1f} "
                    f"{tmax_str} {summary['residue']:<12.2f}"
                )

        if print_table:
            print(sep)
        return results

    # ----------------------------------------------------------------
    # Plotting
    # ----------------------------------------------------------------

    def plot_tg(self, output_path=None, show=True):
        """Plot TG curves for all loaded samples.

        Parameters
        ----------
        output_path : str or Path, optional
            If provided, save the figure to this path.
        show : bool
            Whether to display the plot.
        """
        setup_plot_style()
        fig, ax = plt.subplots(figsize=(8, 5.5))

        for idx, label in enumerate(self.data):
            d = self.data[label]
            color = COLORS[idx % len(COLORS)]
            t_grid, w_grid = downsample_data(d["temp"], d["weight"])
            ax.plot(t_grid, w_grid, color=color, linewidth=1.8, label=label)

        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Weight residual ratio (%)")
        ax.set_xlim(30, 800)
        ax.set_ylim(30, 105)
        ax.legend(fontsize=11, loc="lower left")
        ax.grid(True, linestyle=":", alpha=0.3, linewidth=0.5)

        fig.tight_layout()
        self._save_or_show(fig, output_path, show)

    def plot_dtg(self, output_path=None, show=True):
        """Plot DTG curves for all loaded samples."""
        setup_plot_style()
        fig, ax = plt.subplots(figsize=(8, 5.5))

        for idx, label in enumerate(self.data):
            d = self.data[label]
            color = COLORS[idx % len(COLORS)]
            t_grid, dtg_grid = downsample_data(d["temp"], d["dtg"])
            ax.plot(t_grid, dtg_grid, color=color, linewidth=1.8, label=label)

        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("DTG (%/°C)")
        ax.set_xlim(30, 800)
        ax.legend(fontsize=11, loc="upper right")
        ax.axhline(y=0, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.grid(True, linestyle=":", alpha=0.3, linewidth=0.5)

        fig.tight_layout()
        self._save_or_show(fig, output_path, show)

    def plot_dual_axis(self, label, output_path=None, show=True):
        """Plot dual-axis TGA+DTG for a single sample with peak annotations.

        Parameters
        ----------
        label : str
            Sample label to plot.
        output_path : str or Path, optional
        show : bool
        """
        if label not in self.data:
            print(f"  Sample '{label}' not loaded")
            return

        d = self.data[label]
        temp, weight, dtg = d["temp"], d["weight"], d["dtg"]
        w_norm = weight / weight[0] * 100
        stages, summary = self.find_decomposition_stages(temp, weight, dtg)
        peaks = self.find_peaks(temp, dtg)

        tga_color = "#2166AC"
        dtg_color = "#D6604D"

        setup_plot_style()
        fig, ax1 = plt.subplots(figsize=(9, 6))
        ax2 = ax1.twinx()

        t_plot, w_plot = downsample_data(temp, w_norm, 3000)
        _, dtg_plot = downsample_data(temp, dtg, 3000)

        # TGA (left axis)
        ax1.plot(t_plot, w_plot, color=tga_color, linewidth=2.2, label="TG", zorder=5)
        ax1.set_xlim(0, 800)
        ax1.set_ylim(0, 105)
        ax1.set_xlabel("Temperature (°C)")
        ax1.set_ylabel("Weight residual ratio (%)", color=tga_color)
        ax1.tick_params(axis="y", colors=tga_color)
        ax1.spines["left"].set_color(tga_color)

        # DTG (right axis, positive direction)
        ax2.plot(t_plot, -dtg_plot, color=dtg_color, linewidth=1.8, label="DTG")
        ax2.set_ylabel("DTG (%/°C)", color=dtg_color)
        dtg_max = max(-np.nanmin(dtg_plot) * 1.2, 0.6)
        ax2.set_ylim(0, dtg_max)
        ax2.tick_params(axis="y", colors=dtg_color)
        ax2.spines["right"].set_color(dtg_color)
        ax2.axhline(y=0, color="gray", linewidth=0.8, alpha=0.4)

        # Stage annotations
        stage_colors = ["#4DAF4A", "#377EB8", "#984EA3", "#FF7F00"]
        for i, stage in enumerate(stages):
            sc = stage_colors[i % len(stage_colors)]
            if i == 0 or stage["t_start"] != stages[i - 1]["t_end"]:
                ax1.axvline(x=stage["t_start"], color=sc, linestyle="--",
                            linewidth=0.8, alpha=0.5)
            ax1.axvline(x=stage["t_end"], color=sc, linestyle="--",
                        linewidth=0.8, alpha=0.5)
            t_mid = (stage["t_start"] + stage["t_end"]) / 2
            y_text = [100, 94, 88, 82][i % 4]
            ax1.text(t_mid, y_text, f'{stage["name"]}\n{stage["mass_loss"]:.2f}%',
                     fontsize=8.5, ha="center", va="top", color=sc,
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                               edgecolor=sc, alpha=0.85, linewidth=0.8))

        # DTG peak annotations
        for p in peaks:
            t_p = p["t_peak"]
            d_p = p["dtg_peak"]
            if t_p < 35:
                continue
            t_str = f"{t_p:g}"
            d_str = f"{d_p:g}"
            label_text = f"{t_str}°C\n{d_str}%/°C"

            if t_p < 100:
                offset_x, offset_y = 0, 30
                ha = "center"
            elif t_p < 350:
                offset_x, offset_y = -40, 20
                ha = "right"
            else:
                offset_x, offset_y = 40, 20
                ha = "left"

            ax2.annotate(label_text, xy=(t_p, d_p),
                         xytext=(offset_x, offset_y),
                         textcoords="offset points",
                         fontsize=7.5, fontweight="bold",
                         color=dtg_color, ha=ha, va="bottom",
                         arrowprops=dict(arrowstyle="->", color=dtg_color,
                                        linewidth=1.0, alpha=0.7),
                         bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                   edgecolor=dtg_color, alpha=0.8, linewidth=0.5))

        # T_5%, T_10% annotations
        for kt, kl in [(summary["T_5%"], "T5%"), (summary["T_10%"], "T10%")]:
            if np.isfinite(kt):
                kw = np.interp(kt, temp, w_norm)
                ax1.plot(kt, kw, "o", color=tga_color, markersize=4, zorder=6)
                ax1.annotate(f"{kl}\n{kt:g}°C, {kw:g}%", xy=(kt, kw),
                             xytext=(15, -15), textcoords="offset points",
                             fontsize=7, color=tga_color, ha="left", va="top",
                             arrowprops=dict(arrowstyle="->", color=tga_color,
                                            linewidth=0.8, alpha=0.6))
                ax1.axvline(x=kt, color="gray", linestyle=":", linewidth=0.7, alpha=0.4)

        # Residue annotation
        residue_val = summary["residue"]
        ax1.plot(800, residue_val, "s", color=tga_color, markersize=6, zorder=6)
        ax1.annotate(f"Residue {residue_val:g}%", xy=(800, residue_val),
                     xytext=(-5, -20), textcoords="offset points",
                     fontsize=10, fontweight="bold", color=tga_color,
                     ha="right", va="top",
                     bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFFFCC",
                               edgecolor=tga_color, alpha=0.9, linewidth=1.0))

        # Legend
        legend_elements = [
            plt.Line2D([0], [0], color=tga_color, linewidth=2.2, label="TG"),
            plt.Line2D([0], [0], color=dtg_color, linewidth=1.8, label="DTG"),
        ]
        ax1.legend(handles=legend_elements, fontsize=11, loc="upper right",
                   frameon=True, fancybox=False, edgecolor="gray")

        ax1.set_title(f"TGA-DTG: {label}")
        ax1.spines["top"].set_visible(False)
        ax2.spines["top"].set_visible(False)

        fig.tight_layout()
        self._save_or_show(fig, output_path, show)

    def plot_overlay(self, output_path=None, show=True):
        """Overlay TGA and DTG curves for all samples on a dual-axis plot."""
        setup_plot_style()
        fig, ax1 = plt.subplots(figsize=(9, 6))
        ax2 = ax1.twinx()

        dtg_colors = COLORS

        for idx, label in enumerate(self.data):
            d = self.data[label]
            color = dtg_colors[idx % len(dtg_colors)]
            t_smooth, w_smooth = downsample_data(d["temp"], d["weight"], 3000)
            _, dtg_smooth = downsample_data(d["temp"], d["dtg"], 3000)

            ax1.plot(t_smooth, w_smooth, color=color, linewidth=1.6,
                     label=f"{label} (TG)")
            ax2.plot(t_smooth, -dtg_smooth, color=color, linewidth=1.2,
                     linestyle="--", alpha=0.7, label=f"{label} (DTG)")

        ax1.set_xlim(0, 800)
        ax1.set_ylim(0, 105)
        ax1.set_xlabel("Temperature (°C)")
        ax1.set_ylabel("Weight residual ratio (%)")
        ax1.tick_params(axis="y")

        ax2.set_ylabel("DTG (%/°C)", color="#666666")
        ax2.set_ylim(0, 0.6)
        ax2.tick_params(axis="y")

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=10,
                   loc="lower left", ncol=2, frameon=True, edgecolor="gray")

        ax1.spines["top"].set_visible(False)
        ax2.spines["top"].set_visible(False)

        fig.tight_layout()
        self._save_or_show(fig, output_path, show)

    # ----------------------------------------------------------------
    # Internal
    # ----------------------------------------------------------------

    @staticmethod
    def _save_or_show(fig, output_path, show):
        """Save figure and/or display it."""
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"  Saved: {output_path}")
        if show:
            plt.show()
        plt.close(fig)
