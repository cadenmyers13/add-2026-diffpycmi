import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import (
    Profile, FitContribution, FitRecipe, FitResults
    )
from diffpy.cmi.fit_tools import optimize_recipe
from bg_mpl_stylesheets.styles import all_styles

plt.style.use(all_styles["bg-style"])

from pathlib import Path

from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.structure import constrainAsSpaceGroup
from diffpy.structure.parsers import getParser
from diffpy.utils.parsers.loaddata import loadData

path_to_data_dir = Path(__file__).parent / "pdf" / "ch07StructuralPhaseTransitions" / "data"
cif_path = str(Path(path_to_data_dir) / "SrFe2As2_orthorhombic.cif")
gr150_path = str(Path(path_to_data_dir) / "SrFe2As2_150K.gr")

def make_pdf_recipe(gr_path, cif_path, rmin=.5, rmax=50):
    # Profile
    profile=Profile()
    pdf_parser = PDFParser()
    pdf_parser.parseFile(gr_path)
    profile.loadParsedData(pdf_parser)
    profile.setCalculationRange(xmin=rmin, xmax=rmax, dx=.05)

    # fitcontribution
    contribution = FitContribution("crystal_fit")
    contribution.setProfile(profile)
    cif_parser = getParser("cif")
    crystal_structure = cif_parser.parseFile(cif_path)
    spacegroup = cif_parser.spacegroup.short_name

    model_pdf = PDFGenerator("G")
    model_pdf.setStructure(crystal_structure)

    metadata = loadData(gr_path, headers=True)
    qmax = metadata['qmax']
    qmin = metadata['qmin']
    model_pdf.setQmax(qmax)
    model_pdf.setQmin(qmin)

    qdamp = 0.0349
    qbroad = 0.0176
    model_pdf.qdamp.value = qdamp
    model_pdf.qbroad.value = qbroad
    contribution.addProfileGenerator(model_pdf)
    contribution.setEquation("scale*G")

    # FitRecipe
    recipe = FitRecipe()
    recipe.addContribution(contribution)

    recipe.addVar(contribution.scale, 1, tag="scale")

    sg_parameters = constrainAsSpaceGroup(model_pdf.phase, spacegroup)
    for lat_param in sg_parameters.latpars:
        recipe.addVar(lat_param, tag="lattice_params")
    for adp_param in sg_parameters.adppars:
        recipe.addVar(adp_param, tag="adps")
    for pos_param in sg_parameters.xyzpars:
        recipe.addVar(pos_param, tag="atomic_positions")

    recipe.addVar(model_pdf.delta2, 1.5, tag="delta2")

    return recipe

def refine_recipe_sequentially(recipe, tags):
    recipe.fix("all")
    for tag in tags:
        recipe.free(tag)
        optimize_recipe(recipe)

from scatter_data_fit import plot_recipe_and_print_results

def main():
    pdf_recipe = make_pdf_recipe(gr150_path, cif_path, rmax=10)
    pdf_recipe.show()
    tags = ["lattice_params", "scale", "adps", "atomic_positions", "delta2", "all"]
    refine_recipe_sequentially(pdf_recipe, tags)
    plot_recipe_and_print_results(pdf_recipe)
    plt.show()

if __name__ == "__main__":
    main()

