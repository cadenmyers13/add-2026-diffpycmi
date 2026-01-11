"""
To teach you the basics of diffpy.cmi, we are going to
fit the data that you saw during Simon's presentation.
"""

# 1. First we will import everything we need to do the fit.
#    To ensure we have everything we need, we will install the
#    "plotting" pack of diffpy.cmi.
import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import Profile, FitContribution, FitRecipe, FitResults
from bg_mpl_stylesheets.styles import all_styles
from pathlib import Path
plt.style.use(all_styles["bg-style"])

# 2. Next we will create some synthetic data to fit.
#    We will do this by creating a cubic function and adding
#    some noise to it.
def cubic_function(x, a=-.0015, b=.04, c=.14, d=.5):
    y = a*x**3 + b*x**2 + c*x + d
    return y

savedir = Path(__file__).parent
xmin = 0
xmax = 21.5
x = np.linspace(xmin, xmax, 100)
yobs = cubic_function(x)
noise = np.random.normal(0, 1, size=yobs.shape) * .8
yobs += noise


figsize=(8,6)
plt.figure(figsize=figsize)
plt.legend()
plt.tick_params(labelbottom=False, labelleft=False)
plt.savefig(savedir / "scatterplot-blank.pdf", dpi=300)
plt.show()

plt.figure(figsize=figsize)
plt.plot(x, yobs, 'o')
plt.savefig(savedir / "scatterplot-no-labels.pdf", dpi=300)
plt.show()

plt.figure(figsize=figsize)
plt.plot(x, yobs, 'o')
plt.xlabel("X")
plt.ylabel("Y")
plt.savefig(savedir / "scatterplot-with-labels.pdf", dpi=300)
plt.show()

def plot_recipe(recipe):
    results = FitResults(recipe)
    rw = results.rw
    for name, contrib in recipe._contributions.items():
        profile = contrib.profile
        x = profile.xobs
        yobs = profile.yobs
        ycalc = profile.ycalc
        plt.figure(figsize=figsize)
        plt.plot(x, yobs, 'o')
        plt.plot(x, ycalc, label=f"Rw={rw:.4f}")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.legend()

def make_recipe_linear(x, yobs):
    profile = Profile()
    profile.setObservedProfile(x, yobs)

    contribution = FitContribution("linear_fit")
    contribution.setProfile(profile)
    contribution.setEquation("m * x + b ")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.linear_fit.m, 0.0)
    recipe.addVar(recipe.linear_fit.b, 0.0)

    return recipe

def make_recipe_parabolic(x, yobs):
    profile = Profile()
    profile.setObservedProfile(x, yobs)

    contribution = FitContribution("parabolic_fit")
    contribution.setProfile(profile)
    contribution.setEquation("a * x**2 + b * x + c ")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.parabolic_fit.a, 0.0)
    recipe.addVar(recipe.parabolic_fit.b, 0.0)
    recipe.addVar(recipe.parabolic_fit.c, 0.0)

    return recipe

def make_recipe_cubic(x, yobs):
    profile = Profile()
    profile.setObservedProfile(x, yobs)

    contribution = FitContribution("cubic_fit")
    contribution.setProfile(profile)
    contribution.setEquation("a * x**3 + b * x**2 + c * x + d ")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.cubic_fit.a, 0.0)
    recipe.addVar(recipe.cubic_fit.b, 0.0)
    recipe.addVar(recipe.cubic_fit.c, 0.0)
    recipe.addVar(recipe.cubic_fit.d, 0.0)

    return recipe

def refine_and_plot(recipe):
    from diffpy.cmi.fit_tools import optimize_recipe
    optimize_recipe(recipe)
    plot_recipe(recipe)

if __name__ == "__main__":
    linear_recipe = make_recipe_linear(x, yobs)
    refine_and_plot(linear_recipe)
    plt.savefig(savedir / "scatterplot-linear-fit.pdf", dpi=300)
    plt.show()
    parabolic_recipe = make_recipe_parabolic(x, yobs)
    refine_and_plot(parabolic_recipe)
    plt.savefig(savedir / "scatterplot-parabolic-fit.pdf", dpi=300)
    plt.show()

    cubic_recipe = make_recipe_cubic(x, yobs)
    refine_and_plot(cubic_recipe)
    plt.savefig(savedir / "scatterplot-cubic-fit.pdf", dpi=300)
    plt.show()
