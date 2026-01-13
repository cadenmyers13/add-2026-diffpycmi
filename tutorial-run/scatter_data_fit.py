import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import (
    Profile, FitContribution, FitRecipe, FitResults
    )
from diffpy.cmi.fit_tools import optimize_recipe
from scipy.optimize import leastsq
from bg_mpl_stylesheets.styles import all_styles

plt.style.use(all_styles["bg-style"])

def generate_synthetic_data(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    noise = np.random.normal(0, 1, size=y.shape) * .8
    y_observed = y + noise
    return y_observed

x_observed = np.linspace(0, 21.5, 100)
y_observed = generate_synthetic_data(x_observed, -.0015, .04, .14, .5)

def make_linear_recipe(x, yobs):
    # 1. Profile()
    profile = Profile()
    profile.setObservedProfile(x_observed, y_observed)

    #2. FitContribution
    contribution = FitContribution("linear_fit")
    contribution.setProfile(profile)
    contribution.setEquation("m * x + b")

    #3. FitRecipe()
    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.linear_fit.m, 1)
    recipe.addVar(recipe.linear_fit.b, 1)

    return recipe


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
    return results



def main():
    recipe = make_linear_recipe(x_observed, y_observed)
    optimize_recipe(recipe)
    plot_recipe_and_print_results(recipe)
    plt.show()

if __name__ == "__main__":
    main()