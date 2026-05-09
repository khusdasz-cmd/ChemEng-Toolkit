"""Generate example figures for ChemEng-Toolkit GitHub homepage."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chemeng_toolkit.thermal_analysis import (
    TGAProcessor,
    FTIRProcessor,
    segmented_coats_redfern,
    plot_ea_vs_alpha,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
OUT = SCRIPT_DIR

# ── TGA Figures ────────────────────────────────────────────────────────────

folder_map = [
    (os.path.join(REPO_DIR, "data/tga_sample/gt_100k.xls"), ">100k"),
    (os.path.join(REPO_DIR, "data/tga_sample/100k-20k.xls"), "100k-20k"),
    (os.path.join(REPO_DIR, "data/tga_sample/20k-10k.xls"), "20k-10k"),
    (os.path.join(REPO_DIR, "data/tga_sample/10k-5k.xls"), "10k-5k"),
]
processor = TGAProcessor()
data = {}
for path, label in folder_map:
    d = processor.load_from_xls(path)
    if d:
        data[label] = d
processor.data = data

print("Generating TG plot...")
processor.plot_tg(output_path=os.path.join(OUT, "tga_tg_curves.png"), show=False)

print("Generating DTG plot...")
processor.plot_dtg(output_path=os.path.join(OUT, "tga_dtg_curves.png"), show=False)

print("Generating dual-axis plot...")
processor.plot_dual_axis(">100k", output_path=os.path.join(OUT, "tga_dual_axis.png"), show=False)

print("Generating Ea vs alpha plot...")
seg_results = {}
for label in data:
    d = data[label]
    seg_results[label] = segmented_coats_redfern(d["temp"], d["weight"])
ea_data = {label: [r["E_kJ"] for r in results] for label, results in seg_results.items()}
plot_ea_vs_alpha(ea_data, output_path=os.path.join(OUT, "kinetics_ea_alpha.png"), show=False)

# ── FT-IR Figures ──────────────────────────────────────────────────────────

print("Generating FT-IR stacked plot...")
ftir = FTIRProcessor(data_dir=os.path.join(REPO_DIR, "data/ftir_sample"))
ftir.load_batch([
    ("gt_100k.CSV", "gt-100k", "#E41A1C"),
    ("MBL-20-100K.CSV", "20-100k", "#377EB8"),
    ("MBL-10-20K.CSV", "10-20k", "#4DAF4A"),
    ("MBL-5-10K.CSV", "5-10k", "#984EA3"),
])
ftir.plot_stacked(title="FT-IR Spectra of Lignin Samples",
                  output_path=os.path.join(OUT, "ftir_stacked.png"), show=False)

print("Done — all figures saved to examples/")
