import numpy as np
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

def generate_synthetic_data(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    noise = np.random.normal(0, 1, size=y.shape) * .8
    y_observed = y + noise
    return y_observed

def make_linear_recipe(x, yobs):
    # 1
    profile = Profile()
    profile.setObservedProfile(x, yobs)
    
    # 2
    line_contribution = FitContribution("linear_fit")
    line_contribution.setProfile(profile)

    line_contribution.setEquation("m*x + b")

    # 3
    linear_recipe = FitRecipe()
    linear_recipe.addContribution(line_contribution)

    linear_recipe.addVar(linear_recipe.linear_fit.m, 1)
    linear_recipe.addVar(linear_recipe.linear_fit.b, .5)

    return linear_recipe

def refine_linear_recipe(recipe):
    residual = recipe.residual
    variable_values = recipe.values
    recipe.fithooks[0].verbose = 3
    leastsq(residual, variable_values)
    return

def main():
    x = np.linspace(0, 22, 101)
    yobs = generate_synthetic_data(x, -.0015, .04, .14, .5)
    linear_recipe = make_linear_recipe(x, yobs)
    refine_linear_recipe(linear_recipe)
    fit_results = FitResults(linear_recipe)
    print(fit_results)
    rw = fit_results.rw

    ycalc = linear_recipe.linear_fit.profile.ycalc
    plt.plot(x, yobs, 'o')
    plt.plot(x, ycalc, label=f"Rw={rw}")
    plt.legend()
    plt.show()
    return

if __name__ == "__main__":
    main()