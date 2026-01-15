import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import (
    Profile,
    FitContribution,
    FitRecipe,
    FitResults,
)
from bg_mpl_stylesheets.styles import all_styles

plt.style.use(all_styles["bg-style"])

from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.structure import constrainAsSpaceGroup
from diffpy.structure.parsers import getParser
from diffpy.utils.parsers.loaddata import loadData

from pathlib import Path

path_to_data_dir = Path(__file__).parent / "pdf" / "ch07StructuralPhaseTransitions" / "data"
gr150_path = str(path_to_data_dir / "SrFe2As2_150K.gr")
cif_path = str(path_to_data_dir / "SrFe2As2_orthorhombic.cif")

def make_pdf_recipe(gr_path, cif_path):
    # 1
    profile = Profile()
    pdf_parser = PDFParser()
    pdf_parser.parseFile(gr150_path)
    profile.loadParsedData(pdf_parser)

    g = profile.y
    r = profile.x

    profile.setCalculationRange(xmin=.5, xmax=20, dx=.05)
    
    # 2
    pdf_contribution = FitContribution("crystal_fit")
    pdf_contribution.setProfile(profile)
    cif_parser = getParser("cif")
    ortho_structure = cif_parser.parseFile(cif_path)
    spacegroup = cif_parser.spacegroup.short_name
    print(ortho_structure)
    print(spacegroup)

    model_pdf = PDFGenerator("G")
    model_pdf.setStructure(ortho_structure)

    metadata = loadData(gr_path, headers=True)
    qmax = metadata.get('qmax')
    qmin = metadata.get('qmin')
    model_pdf.setQmax(qmax)
    model_pdf.setQmin(qmin)

    qdamp = 0.0349
    qbroad = 0.0176
    model_pdf.qdamp.value = qdamp
    model_pdf.qbroad.value = qbroad
    pdf_contribution.addProfileGenerator(model_pdf)
    pdf_contribution.setEquation("s*G")

    # 3
    pdf_recipe = FitRecipe()
    pdf_recipe.addContribution(pdf_contribution)

    sg_parameters = constrainAsSpaceGroup(model_pdf.phase, spacegroup)

    for lat_param in sg_parameters.latpars:
        pdf_recipe.addVar(lat_param, tag="lattice_params")
    for adp_param in sg_parameters.adppars:
        pdf_recipe.addVar(adp_param, tag="adps")
    for pos_param in sg_parameters.xyzpars:
        pdf_recipe.addVar(pos_param, tag="atomic_coordinates")

    pdf_recipe.addVar(pdf_contribution.s, 1, tag="scale")
    pdf_recipe.addVar(model_pdf.delta2, 1.5, tag="delta2")

    return pdf_recipe

from scatter_data_fit import refine_recipe, plot_recipe

def run_staged_refinement(recipe):
    recipe.fix("all")
    tags = ["lattice_params", "scale", "adps", "atomic_coordinates", "delta2", "all"]
    for tag in tags:
        recipe.free(tag)
        refine_recipe(recipe)

def main():
    recipe = make_pdf_recipe(gr150_path, cif_path)
    run_staged_refinement(recipe)
    plot_recipe(recipe)
    return

if __name__ == "__main__":
    main()