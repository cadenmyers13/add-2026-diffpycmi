"""Now that we've gotten familiar with diffpy.cmi, lets start fitting
PDFs. We will start with the PDF of SrFe2As2.

To obtain the data for this example, run `cmi copy
ch07StructuralPhaseTransitions`
"""

# -------------------------------------------------------
# First thing we must do is import all necessary packages.
# This includes everything from the previous example plus
# additional imports necessary for PDF fitting.
# Currently, we do not have all necessary packages installed,
# install these packages with `cmi install pdf`.

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

# -------------------------------------------------------
# Now for the sanity check, plot it from the profile object
g = profile.yobs
r = profile.xobs

plt.plot(r, g, "o")
plt.show()

# -------------------------------------------------------
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

# -------------------------------------------------------
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

# -------------------------------------------------------
# Now, we want to generate our model PDF from the cif file.
# To do this, we will use the PDFGenerator() object.
model_pdf = PDFGenerator("G1")
model_pdf.setStructure(tetragonal_structure)
# -------------------------------------------------------
