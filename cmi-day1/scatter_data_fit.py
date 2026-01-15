import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import (
    Profile,
    FitContribution,
    FitRecipe,
    FitResults,
)
from scipy.optimize import leastsq
from bg_mpl_stylesheets.styles import all_styles

plt.style.use(all_styles["bg-style"])


def generate_synthetic_data(x, a, b, c, d):
    y = a * x**3 + b * x**2 + c * x + d
    noise = np.random.normal(0, 1, size=y.shape) * 0.8
    y_observed = y + noise
    return y_observed


x = np.linspace(0, 21.5, 100)
yobs = generate_synthetic_data(x, -0.0015, 0.04, 0.14, 0.5)

# 1
scatter = Profile()
scatter.setObservedProfile(x, yobs)

# 2
line = FitContribution("linear_fit")
line.setProfile(scatter)

line.setEquation("m*x + b")

# 3
linear_recipe = FitRecipe()
linear_recipe.addContribution(line)

linear_recipe.addVar(linear_recipe.linear_fit.m, 1)
linear_recipe.addVar(linear_recipe.linear_fit.b, 1)

# 4
residual = linear_recipe.residual
values = linear_recipe.values
leastsq(residual, values)
fit_results = FitResults(linear_recipe)


def make_linear_recipe(x, yobs):

    profile = Profile()
    profile.setObservedProfile(x, yobs)

    contribution = FitContribution("linear_fit")
    contribution.setProfile(profile)
    contribution.setEquation("m*x + b")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.linear_fit.m, 1)
    recipe.addVar(recipe.linear_fit.b, 1)
    return recipe


def plot_fit(linear_recipe, fit_results):
    rw = fit_results.rw
    x = linear_recipe.linear_fit.profile.x
    ycalc = linear_recipe.linear_fit.profile.ycalc
    yobs = linear_recipe.linear_fit.profile.y
    plt.plot(x, yobs, "o")
    plt.plot(x, ycalc, label=f"Rw={round(rw, 4)}")
    plt.legend()
    plt.show()
    return


def refine_recipe(recipe):
    residual = recipe.residual
    values = recipe.values
    leastsq(residual, values)
    fit_results = FitResults(recipe)
    return fit_results


def plot_recipe(recipe, figsize=(8, 6), offset_scale=1.0):
    results = FitResults(recipe)
    rw = results.rw
    for name, contrib in recipe._contributions.items():
        profile = contrib.profile
        x = profile.x
        yobs = profile.y
        ycalc = profile.ycalc
        diff = yobs - ycalc
        base_offset = min(yobs.min(), ycalc.min()) - 0.1 * (
            yobs.max() - yobs.min()
        )
        offset = base_offset * offset_scale
        plt.figure(figsize=figsize)
        plt.plot(x, yobs, "o")
        plt.plot(x, ycalc, label=f"Rw={rw:.4f}")
        plt.plot(x, diff + offset)
        plt.axhline(offset, color="black")
        plt.legend()
        plt.show()
    return results


def main():
    recipe = make_linear_recipe(x, yobs)
    fit_results = refine_recipe(recipe)
    plot_recipe(recipe)
    return


if __name__ == "__main__":
    main()
