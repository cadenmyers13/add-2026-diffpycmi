from pdf_fit import (
    make_pdf_recipe,
    run_staged_refinement,
    path_to_data_dir,
    cif_path,
)
from diffpy.srfit.fitbase import FitResults
from diffpy.utils.parsers.loaddata import loadData
import matplotlib.pyplot as plt

ortho_cif = cif_path
tet_cif = path_to_data_dir / "SrFe2As2_tetragonal.cif"


def get_sorted_gr_files(data_dir):
    gr_files = sorted(data_dir.glob("*.gr"))
    return gr_files


def plot_waterfall(gr_files):
    offset = 0
    for file in gr_files:
        data = loadData(str(file))
        r = data[:, 0]
        g = data[:, 1]
        plt.plot(r, g + offset)
        offset += 5
    plt.show()


def save_fit_results(fitresults, filename, data_dir=path_to_data_dir):
    results_dir = data_dir / "fit-results"
    results_dir.mkdir(exist_ok=True)
    fitresults.saveResults(str(results_dir / f"{filename}.res"))


def main():
    gr_files = get_sorted_gr_files(path_to_data_dir)
    plot_waterfall(gr_files)
    cif_files = [ortho_cif, tet_cif]
    results_dict = {}
    for file in gr_files:
        for cif in cif_files:
            recipe = make_pdf_recipe(str(file), str(cif))
            recipe.fithooks[0].verbose = 0
            run_staged_refinement(recipe)
            fit_results = FitResults(recipe)
            phase = "ortho" if cif == ortho_cif else "tet"
            save_fit_results(fit_results, f"{phase}_{str(file.stem)}")
            print(file.stem)
            print(phase)
            print(fit_results)


if __name__ == "__main__":
    main()
