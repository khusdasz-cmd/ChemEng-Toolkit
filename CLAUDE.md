# ChemEng-Toolkit — CLAUDE.md

## Project Overview

Python toolkit for chemical engineering data analysis — TGA/DSC processing, FT-IR spectroscopy, pyrolysis kinetics (Coats-Redfern), and distillation column design. Built as an open-source portfolio project for undergraduate research applications.

- **GitHub**: https://github.com/khusdasz-cmd/ChemEng-Toolkit
- **Author**: Lin Haokang (khusdasz@gmail.com)
- **Python**: >= 3.9
- **Status**: Active development

## Build & Install

```bash
# Editable install (development)
pip install -e . --no-build-isolation

# Or with conda
conda install -c conda-forge numpy pandas matplotlib scipy xlrd openpyxl
pip install -e . --no-build-isolation
```

## Project Structure

```
ChemEng-Toolkit/
├── chemeng_toolkit/
│   ├── thermal_analysis/
│   │   ├── tga_processor.py      # TGAProcessor class: load .xls, TG/DTG plots, peak finding, stage analysis, summary table
│   │   ├── kinetics.py           # coats_redfern(), run_coats_redfern_batch(), segmented_coats_redfern(), plot_ea_vs_alpha()
│   │   └── ftir_processor.py     # FTIRProcessor class: load CSV, stacked offset plots, peak annotations
│   └── utils/
│       └── plot_helpers.py       # COLORS, FTIR_PEAKS, setup_plot_style(), downsample_data()
├── examples/
│   ├── tga_example.ipynb         # TGA pipeline walkthrough
│   └── ftir_example.ipynb        # FT-IR pipeline walkthrough
└── data/
    ├── tga_sample/               # Sample .xls files (TA Instruments format)
    └── ftir_sample/              # Sample .CSV files (wavenumber, absorbance)
```

## Architecture Decisions

| Decision | Rationale |
|---|---|
| **Class-based API** for data loading/plotting (`TGAProcessor`, `FTIRProcessor`) | Encapsulates shared state (samples, data, color mapping) |
| **Functional API** for analysis (`coats_redfern()`, `segmented_coats_redfern()`) | Pure functions with no side effects — easy to chain and test |
| **English-only** output (labels, docstrings, variables) | International portfolio project |
| **Default matplotlib** (no Chinese font config) | Removes Windows-specific font dependency |
| **pyproject.toml** over setup.py | Modern Python packaging standard |
| **Scipy for smoothing** (not custom) | Leverages validated numerical methods |

## Coding Conventions

- **Language**: English everywhere — code, comments, docstrings, plot labels
- **Docstrings**: Google-style, one-line summary + Args/Returns sections
- **Naming**: `snake_case` for functions/variables, `PascalCase` for classes, `UPPER_CASE` for constants
- **Types**: Use type hints for all public function signatures
- **Imports**: Standard lib → third-party → local; absolute imports preferred
- **Matplotlib**: Always use `plt.ioff()` + `plt.close()` in library code (never leave figures open)
- **Error handling**: Validate at public API boundaries; trust internal calls

## Key Modules Reference

### tga_processor.py
- `TGAProcessor(samples=None, sheet_name="Merged 1", ...)` — constructor
- `.load_from_xls(filepath)` → dict with keys: time, temp, weight, dtg
- `.load_batch(folder_map)` — load multiple samples from a dict of `{label: folder_path}`
- `.plot_tg()`, `.plot_dtg()`, `.plot_dual_axis(label)`, `.plot_overlay()`
- `.find_decomposition_stages(temp, weight, dtg)` — identify stage boundaries from DTG peaks
- `.summary_table()` — extract T_5%, T_10%, T_max, residue

### kinetics.py
- `coats_redfern(temp, weight, t_range, alpha_range, heating_rate)` — fits 15 models, returns list sorted by R²
- `run_coats_redfern_batch(data_dict)` — batch analysis returning best models + F1 results
- `segmented_coats_redfern(temp, weight, alpha_values, half_window, t_range)` — alpha-dependent Ea using F1 model
- `plot_ea_vs_alpha(ea_data, labels)` — publication-quality Ea vs α plot
- 15 models: F1-F3, R2-R3, D1-D4, A2-A4, P2-P4

### ftir_processor.py
- `FTIRProcessor(files=None, data_dir=None)` — constructor
- `.load_csv(filepath)` → (wavenumber, absorbance) tuple
- `.load_batch(file_list)` — load list of `(path, label, color)` tuples
- `.plot_stacked(normalize, offset_step, title)` — vertical offset stacked spectra
- `.plot_individual(normalize)` — separate subplots per sample

### plot_helpers.py
- `COLORS = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', ...]`
- `FTIR_PEAKS = [(3400,"O-H"), (2947,"C-H"), (1682,"C=O"), (1600,"Aromatic"), ...]`
- `setup_plot_style()` — sets matplotlib rcParams for publication quality
- `downsample_data(x, y, max_points=5000)` — evenly spaced downsampling
- `SAMPLE_LABELS = {">100k", "100k-20k", "20k-10k", "10k-5k"}`

## Git Workflow

- **Remote**: `origin` → https://github.com/khusdasz-cmd/ChemEng-Toolkit.git
- **Branch**: work on feature branches, merge to main
- **Commit style**: Concise English, present tense ("Add X", "Fix Y", "Refactor Z")
- **Avoid**: Force push to main, --no-verify hooks

## Environment Notes

- **OS**: Windows 10, paths use backslashes externally but forward slashes in Python
- **Shell**: Git Bash (available in PATH)
- **Package manager**: pip (with `--no-build-isolation` flag needed for editable installs due to proxy)
- **Proxy**: System may have proxy restrictions — pip install may need `--proxy` or `--no-build-isolation`
- **Jupyter**: Notebooks in `examples/` run with kernel from the same environment

## Common Tasks

```bash
# Install locally
pip install -e . --no-build-isolation

# Run all examples (from repo root)
jupyter nbconvert --execute examples/tga_example.ipynb --to notebook --output tga_example_output.ipynb

# Check imports
python -c "from chemeng_toolkit.thermal_analysis import TGAProcessor, FTIRProcessor; print('OK')"

# Git push
git push origin main
```
