"""Temperature-dependent PDF fitting for SrFe2As2 phase transition
study.

Fits each temperature's PDF data with both tetragonal and orthorhombic
structures to track the structural phase transition.
"""

# --------------------------------------------------------------
# First, let's import everything we'll need
from pdf_fit import make_pdf_recipe, refine_recipe_sequentially
from diffpy.srfit.fitbase import FitResults
from diffpy.utils.parsers.loaddata import loadData
import numpy as np
from pathlib import Path
import glob
import matplotlib.pyplot as plt
from bg_mpl_stylesheets.styles import all_styles

plt.style.use(all_styles["bg-style"])

# --------------------------------------------------------------
# Configuration: set up file paths and fitting parameters
DATA_DIR = (
    Path(__file__).parent / "pdf" / "ch07StructuralPhaseTransitions" / "data"
)
TETRAGONAL_CIF = str(DATA_DIR / "SrFe2As2_tetragonal.cif")
ORTHORHOMBIC_CIF = str(DATA_DIR / "SrFe2As2_orthorhombic.cif")
REFINEMENT_TAGS = [
    "scale",
    "lattice_params",
    "adps",
    "atomic_positions",
    "delta2",
]
R_MIN, R_MAX = 0.5, 50

# --------------------------------------------------------------
# Before we start, lets write a function that plots all the gr files in a waterfall plot
# to see what the data look like


def get_sorted_gr_files(data_dir):
    """Get all .gr files sorted by filename."""
    return sorted(glob.glob(str(data_dir / "*.gr")))


def get_temperature_from_filename(filepath):
    """Extract temperature value from filename like
    'SrFe2As2_150K.gr'."""
    stem = Path(filepath).stem
    temp_str = stem.split("_")[-1].replace("K", "")
    return int(temp_str)


def plot_waterfall(gr_files, offset_step=5):
    """Plot all .gr files in a waterfall plot to visualize data."""
    plt.figure(figsize=(8, 6))
    offset = 0
    for gr_file in gr_files:
        temp = get_temperature_from_filename(gr_file)
        data = loadData(gr_file)
        r = data[:, 0]
        g = data[:, 1]
        plt.plot(r, g + offset, label=f"{temp} K")
        offset += offset_step  # Adjust offset for visibility
    plt.xlabel("r")
    plt.ylabel("G")
    plt.legend()
    plt.show()


gr_files = get_sorted_gr_files(DATA_DIR)
# plot_waterfall(gr_files)
# --------------------------------------------------------------
# Our process to plot U11 vs temperature requires:
# 1. Extracting temperature from filename
# 2. Fitting each temperature's data
# 3. Collecting results (Rw and U11 values)
# 4. Plotting the results


def fit_single_file(gr_file, cif_file, rmin, rmax):
    """Fit a single temperature's data and return results."""
    recipe = make_pdf_recipe(gr_file, cif_file, rmin=rmin, rmax=rmax)
    refine_recipe_sequentially(recipe, REFINEMENT_TAGS)
    results = FitResults(recipe)
    return recipe, results


def save_fit_results(fitresults, filename, data_dir):
    results_dir = data_dir / "fit-results"
    results_dir.mkdir(exist_ok=True)
    fitresults.saveResults(results_dir / f"{filename}.res")


def main():
    gr_files = get_sorted_gr_files(DATA_DIR)
    plot_waterfall(gr_files)
    cif_files = [ORTHORHOMBIC_CIF, TETRAGONAL_CIF]
    for gr in gr_files:
        for cif in cif_files:
            recipe, results = fit_single_file(gr, cif, R_MIN, R_MAX)
            temp = get_temperature_from_filename(gr)
            phase = "ortho" if cif == ORTHORHOMBIC_CIF else "tet"
            save_fit_results(results, f"{phase}_{temp}K", DATA_DIR)


if __name__ == "__main__":
    main()
