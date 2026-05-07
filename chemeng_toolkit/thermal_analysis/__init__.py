"""Thermal analysis modules: TGA/DSC processing, kinetics, FT-IR spectroscopy."""

from .tga_processor import TGAProcessor
from .kinetics import coats_redfern, run_coats_redfern_batch, segmented_coats_redfern, plot_ea_vs_alpha
from .ftir_processor import FTIRProcessor

__all__ = [
    "TGAProcessor",
    "coats_redfern",
    "run_coats_redfern_batch",
    "segmented_coats_redfern",
    "plot_ea_vs_alpha",
    "FTIRProcessor",
]
