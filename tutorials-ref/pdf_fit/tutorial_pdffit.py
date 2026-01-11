"""Now that we've gotten familiar with diffpy.cmi, lets start fitting
PDFs. We will start with the PDF of SrFe2As2.

To obtain the data for this example, run `cmi copy
ch07StructuralPhaseTransitions`
"""

# -------------------------------------------------------
# First thing we must do is import all necessary packages.
# This includes everything from the previous example plus
# additional imports necessary for PDF fitting.


import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import (
    Profile,
    FitContribution,
    FitRecipe,
    FitResults,
)
from diffpy.cmi.fit_tools import optimize_recipe
from bg_mpl_stylesheets.styles import all_styles

plt.style.use(all_styles["bg-style"])

from pathlib import Path

# these are the new tools for this PDF fit
from diffpy.srfit.pdf import PDFGenerator, PDFParser
from diffpy.srfit.structure import constrainAsSpaceGroup
from diffpy.structure.parsers import getParser
from diffpy.utils.parsers.loaddata import loadData

# -------------------------------------------------------
# First thing is to load in and plot the PDF as a sanity check.
# We can get the path to the data using pathlib.Path
gr_path = str(Path(__file__).parent / "data" / "SrFe2As2_150K.gr")

# -------------------------------------------------------
# Like we did with the cubic fit example, we will first instatiate
# a Profile object. Then, we parse and load in the PDF
# and link it to the Profile object.

profile = Profile()
# parse the pdf
pdf_parser = PDFParser()
pdf_parser.parseFile(gr_path)
# load pdf to profile object
profile.loadParsedData(pdf_parser)

# Now for the sanity check, plot it from the profile object
g = profile.yobs
r = profile.xobs

plt.plot(r, g, "o")
plt.title("PDF from Profile Object")
plt.show()

# As you can see, the signal goes all the way up to r=100A.
# We can limit the range of the data we want to fit by setting
# the calculation range on the profile object.

rmin = 0.5
rmax = 50
rstep = 0.05
profile.setCalculationRange(xmin=rmin, xmax=rmax, dx=rstep)

# -------------------------------------------------------
# Now that we've loaded our data and added it to the profile,
# we can start working on our model.
# Like last time, we will first instatiate the FitContribution and give it
# a name.
contribution = FitContribution("crystal_fit")

# Last time, the models we were fitting were nth order equations.
# This time, we are fitting a PDF based on a crystal structure.
# Therefore, we need to generate a model PDF from a structural model.

# This process includes several steps, so we will take things slowly.
# 1. Get the path to the cif file
# 2. Load the cif parser
# 3. Parse the cif file to get the structure and spacegroup
# 4. Generate the model PDF from the structure using PDFGenerator
# 5. Add the model PDF to the FitContribution object

# First off, our structural model is in the form of a cif file, which is
# a standardized format for reading in crystal structures.
# To do anything with the cif file, we first must obtain the path
cif_path = str(Path(__file__).parent / "data" / "SrFe2As2_tetragonal.cif")

# To get all the important information from the cif file in a usable format,
# we will use the getParser method from diffpy.structure
cif_parser = getParser("cif")
# to get the structure run this
tetragonal_structure = cif_parser.parseFile(cif_path)
# to get the spacegroup run this, we will save this for later
spacegroup = cif_parser.spacegroup.short_name

# Now, we want to generate our model PDF from the cif file.
# To do this, we will use the PDFGenerator() object.
model_pdf = PDFGenerator("G1")
model_pdf.setStructure(tetragonal_structure)

# When obtaining your PDF experimentally, there are 
# two parameter that are determined from your experimental setup
# These are Qmin and Qmax.

# Using loadData from diffpy.utils.parsers.loaddata, we can
# obtain these parameters from the header of the .gr file.
metadata = loadData(gr_path, headers=True)
qmax = metadata["qmax"]
qmin = metadata["qmin"]
# Now we can set these parameters on the model PDF object
model_pdf.setQmax(qmax)
model_pdf.setQmin(qmin)

# In addition to Qmin and Qmax, there are two more parameters
# that are determined through a calibrant sample (Ni or CeO2 usually).
# These are Qdamp and Qbroad. These values are given to us in this example,
# but usually you get this from refining on these parameters in a Ni or CeO2
# fit.
qdamp = 0.0349
qbroad = 0.0176
model_pdf.qdamp.value = qdamp
model_pdf.qbroad.value = qbroad

# With the PDFGenerator object created, we can now add it to
# the FitContribution object and link the profile to the contribution.
contribution.addProfileGenerator(model_pdf)
contribution.setProfile(profile)

# To scale the generated PDF to the observed PDF, we need to add
# a scale factor to the contribution in the form of an equation.
contribution.setEquation("s1*G1")

contribution.show()
# -------------------------------------------------------
# Great, we have now loaded our data, created a model and linked them
# together. Now, we need to set up our FitRecipe to perform the fit.
recipe = FitRecipe()
recipe.addContribution(contribution)

# now, we can start adding these parameters as variables to the recipe
recipe.addVar(contribution.s1, 1.0, tag="scale")

sg_parameters = constrainAsSpaceGroup(model_pdf.phase, spacegroup)

for lat_param in sg_parameters.latpars:
    recipe.addVar(lat_param, tag="lattice_params")
for adp_param in sg_parameters.adppars:
    recipe.addVar(adp_param, tag="adps")
for pos_param in sg_parameters.xyzpars:
    recipe.addVar(pos_param, tag="atomic_positions")

recipe.addVar(model_pdf.delta2, 1.5, tag="delta2")

recipe.show()

# -------------------------------------------------------
# Now that we have our recipe set up, we can perform the fit.
# We will do this sequentially, first fitting only the scale factor,
# then fitting the lattice parameters, then the atomic positions, etc

recipe.fix("all")
tags = ["lattice_params", "scale", "adps", "atomic_positions", "delta2", "all"]
for tag in tags:
    recipe.free(tag)
    optimize_recipe(recipe)

# -------------------------------------------------------
results = FitResults(recipe)
rw = results.rw
robs = recipe.crystal_fit.profile.x
gobs = recipe.crystal_fit.profile.y
gcalc = recipe.crystal_fit.profile.ycalc
plt.plot(robs, gobs, "o", )
plt.plot(robs, gcalc, "-", label=f"Rw={rw:.4f}")
# plot the difference below the fit
diff = gobs - gcalc
offset = np.min(gobs) - 0.5
plt.plot(robs, diff + offset, "-")
plt.hlines(offset, xmin=rmin, xmax=rmax, colors="k")
plt.title("PDF Fit of tetragonal SrFe2As2 at 150K")
plt.xlabel("r (A)")
plt.ylabel("G(r)")
plt.legend()
plt.show()

# -------------------------------------------------------
# Just like we did with the line fit, lets wrap this all in a function

def make_recipe(gr_path, cif_path, rmin, rmax):
    profile = Profile()
    pdf_parser = PDFParser()
    pdf_parser.parseFile(gr_path)
    profile.loadParsedData(pdf_parser)
    profile.setCalculationRange(xmin=rmin, xmax=rmax)

    contribution = FitContribution("crystal_fit")
    cif_parser = getParser("cif")
    structure = cif_parser.parseFile(cif_path)
    spacegroup = cif_parser.spacegroup.short_name
    model_pdf = PDFGenerator("G1")
    model_pdf.setStructure(structure)
    metadata = loadData(gr_path, headers=True)
    qmax = metadata["qmax"]
    qmin = metadata["qmin"]
    model_pdf.setQmax(qmax)
    model_pdf.setQmin(qmin)
    qdamp = 0.0349
    qbroad = 0.0176
    model_pdf.qdamp.value = qdamp
    model_pdf.qbroad.value = qbroad
    contribution.addProfileGenerator(model_pdf)
    contribution.setProfile(profile)
    contribution.setEquation("s1*G1")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(contribution.s1, 1.0, tag="scale")
    sg_parameters = constrainAsSpaceGroup(model_pdf.phase, spacegroup)
    for lat_param in sg_parameters.latpars:
        recipe.addVar(lat_param, tag="lattice_params")
    for adp_param in sg_parameters.adppars:
        recipe.addVar(adp_param, tag="adps")
    for pos_param in sg_parameters.xyzpars:
        recipe.addVar(pos_param, tag="atomic_positions")
    recipe.addVar(model_pdf.delta2, 1.5, tag="delta2")

    return recipe

def refine_recipe(recipe):
    recipe.fix("all")
    tags = ["lattice_params", "scale", "adps", "atomic_positions", "delta2", "all"]
    for tag in tags:
        recipe.free(tag)
        optimize_recipe(recipe)

def plot_recipe_and_print_results(recipe, figsize=(8, 6)):
    results = FitResults(recipe)
    print(results)
    rw = results.rw
    for name, contrib in recipe._contributions.items():
        profile = contrib.profile
        xobs = profile.x
        yobs = profile.y
        ycalc = profile.ycalc
        plt.figure(figsize=figsize)
        plt.plot(xobs, yobs, "o")
        plt.plot(xobs, ycalc, label=f"Rw={rw:.4f}")
        diff = yobs - ycalc
        offset = np.min(yobs) - 0.5
        plt.plot(xobs, diff + offset, "-")
        plt.hlines(offset, xmin=np.min(xobs), xmax=np.max(xobs), colors="k")
        plt.title("PDF Fit Result")
        plt.xlabel("r (A)")
        plt.ylabel("G(r)")
        plt.legend()
        plt.show()

# Now we can use these functions to do the fit again

def main():
    recipe = make_recipe(gr_path, cif_path, 0, 30)
    refine_recipe(recipe)
    plot_recipe_and_print_results(recipe)

if __name__ == "__main__":
    main()