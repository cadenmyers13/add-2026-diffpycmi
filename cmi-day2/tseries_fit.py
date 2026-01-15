"""Please see cmi-day1 for a more complete version."""

from pdf_fit import make_pdf_recipe, run_staged_refinement
from scatter_data_fit import plot_recipe
from diffpy.srfit.fitbase import FitResults
from diffpy.utils.parsers.loaddata import loadData
from pathlib import Path
import matplotlib.pyplot as plt

data_dir = Path(__file__).parent / "pdf" / "ch07StructuralPhaseTransitions" / "data"
ortho_path = str(data_dir / "SrFe2As2_orthorhombic.cif")
tet_cif = str(data_dir / "SrFe2As2_tetragonal.cif")


def get_sorted_gr_files(data_dir):
    gr_files = sorted(data_dir.glob("*.gr"))
    return gr_files

def main():
    gr_files = get_sorted_gr_files(data_dir)
    cif_files = [ortho_path, tet_cif]
    for file in gr_files:
        for cif in cif_files:
            recipe = make_pdf_recipe(file, cif)
            run_staged_refinement(recipe)

if __name__ == "__main__":
    main()