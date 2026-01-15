import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import (
    Profile,
    FitContribution,
    FitRecipe,
    FitResults
)

from scipy.optimize import leastsq
from bg_mpl_stylesheets.styles import all_styles
plt.style.use(all_styles["bg-style"])

from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.structure import constrainAsSpaceGroup
from diffpy.structure.parsers import getParser
from diffpy.utils.parsers.loaddata import loadData
from pathlib import Path

def plot_gr(gr_file):
    """Plots the given .gr file.
    
    Parameters
    ----------
    gr_file: str
        The path to .gr file.
        
    Return
    ------
    None
    """
    r, g = loadData(gr_file, unpack=True)
    plt.plot(r, g, label=Path(gr_file).stem)

def make_pdf_recipe(gr_file, cif_file):
    # 1
    profile = Profile()
    pdf_parser = PDFParser()
    pdf_parser.parseFile(gr_file)
    profile.loadParsedData(pdf_parser)
    profile.setCalculationRange(0.5, 20, dx=.05)

    # 2
    contribution_pdf = FitContribution("crystal_fit")
    contribution_pdf.setProfile(profile)
    cif_parser = getParser("cif")
    structure = cif_parser.parseFile(cif_file)
    spacegroup = cif_parser.spacegroup.short_name

    model_pdf = PDFGenerator("G")
    model_pdf.setStructure(structure)

    metadata = loadData(gr_file, headers=True)
    qmax = metadata.get("qmax")
    qmin = metadata.get("qmin")
    model_pdf.setQmax(qmax)
    model_pdf.setQmin(qmin)

    qdamp = 0.0349
    qbroad = 0.0176
    model_pdf.qdamp.value = qdamp
    model_pdf.qbroad.value = qbroad

    contribution_pdf.addProfileGenerator(model_pdf)
    contribution_pdf.setEquation("s*G")

    # contribution_pdf.show()
    # 3
    recipe_pdf = FitRecipe()
    recipe_pdf.addContribution(contribution_pdf)

    sg_parameters = constrainAsSpaceGroup(model_pdf.phase, spacegroup)

    for lat_param in sg_parameters.latpars:
        recipe_pdf.addVar(lat_param, tag="lattice_params")
    for adp_param in sg_parameters.adppars:
        recipe_pdf.addVar(adp_param, tag="adps")
    for pos_param in sg_parameters.xyzpars:
        recipe_pdf.addVar(pos_param, tag="atomic_coordinates")

    recipe_pdf.addVar(recipe_pdf.crystal_fit.s, 1)
    recipe_pdf.addVar(model_pdf.delta2, 1.5)
    return recipe_pdf

from scatter_data_fit import refine_linear_recipe

def plot_recipe_and_print_results(recipe, figsize=(8, 6), offset_scale=1.0):
    results = FitResults(recipe)
    print(results)
    rw = results.rw
    for name, contrib in recipe._contributions.items():
        profile = contrib.profile
        x = profile.x
        yobs = profile.y
        ycalc = profile.ycalc
        diff = yobs - ycalc
        base_offset = min(yobs.min(), ycalc.min()) - 0.1 * (yobs.max() - yobs.min())
        offset = base_offset * offset_scale
        plt.figure(figsize=figsize)
        plt.plot(x, yobs, "o")
        plt.plot(x, ycalc, label=f"Rw={rw:.4f}")
        plt.plot(x, diff + offset)
        plt.axhline(offset, color="black")
        plt.legend()
        plt.show()
    return results

def run_staged_refinement(recipe):
    recipe.fix("all")
    tags = ["lattice_params", "s", "adps", "atomic_coordinates", "delta2", "all"]
    for tag in tags:
        recipe.free(tag)
        refine_linear_recipe(recipe)

def main():
    data_dir = Path(__file__).parent / "pdf" / "ch07StructuralPhaseTransitions" / "data"
    gr150_path = str(data_dir / "SrFe2As2_150K.gr")
    ortho_path = str(data_dir / "SrFe2As2_orthorhombic.cif")
    recipe = make_pdf_recipe(gr150_path, ortho_path)
    run_staged_refinement(recipe)
    plot_recipe_and_print_results(recipe)

    return

if __name__ == "__main__":
    main()