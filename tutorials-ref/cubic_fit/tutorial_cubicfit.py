"""
To teach you the basics of diffpy.cmi, we are going to
fit the data that you saw during Simon's presentation.
"""
# ------------------------------------------------------
# 1. First we will import everything we need to do the fit.
#    To ensure we have everything we need, we will install the
#    "plotting" pack of diffpy.cmi by running `cmi install plotting`.
import numpy as np
import matplotlib.pyplot as plt
from diffpy.srfit.fitbase import Profile, FitContribution, FitRecipe, FitResults
from scipy.optimize import leastsq
from bg_mpl_stylesheets.styles import all_styles
plt.style.use(all_styles["bg-style"])
 
# ------------------------------------------------------
# 2. Next we will create some synthetic data to fit.
#    We will do this by creating a cubic function and adding
#    some noise to it.
def generate_synthetic_data(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    noise = np.random.normal(0,1, size=y.shape) * .8
    y_observed = y + noise
    return y_observed

xobs = np.linspace(0, 21.5, 100)
yobs = generate_synthetic_data(xobs, -.0015, .04, .14, .5)

# ------------------------------------------------------
# lets plot the data to see what this looks like
plt.plot(xobs, yobs, "o")
plt.show()
# ------------------------------------------------------
# 3. As humans, we have inherent bias in how we might chose to go about
#    fitting this data. Some of us may be biased to fit a line, others
#    a parabola, and others a cubic. Because we created the data using
#    a cubic function, we already are biased to see the underlying cubic
#    curvature of the plot.
#    Lets ignore our bias and use our scientific reason to find the function
#    that best fits the model.

# First, lets try to fit the data with a linear function. 
# This is where we will first use objects from the diffpy.cmi world

# Instatiate a Profile object and set the observed data.
profile = Profile()
profile.setObservedProfile(xobs, yobs)

# we can access our data directly from the profile now.
# this can be done like so
x_from_profile = profile.xobs
y_from_profile = profile.yobs

# lets plot them to see if we get what we expect
plt.plot(x_from_profile, y_from_profile, 'o')
plt.title("Data from Profile object")
plt.show()

# ------------------------------------------------------
# 4. Next, we create a FitContribution object to hold the profile.
#    We will give it the name "linear_fit".
contribution = FitContribution("linear_fit")
# we must link the data to the fit contribution
contribution.setProfile(profile)
contribution.show()

# ------------------------------------------------------
# 5. To see what the contribution looks like, we can call show().
# and now we define our linear function in the form of a string
contribution.setEquation("m * x + b ")
contribution.show()
# ------------------------------------------------------
# We can access the data and equation parameters from contribution
# as well
y_from_contrib = contribution.profile.yobs
x_from_contrib = contribution.profile.xobs
m_from_contrib = contribution.m.value
b_from_contrib = contribution.b.value


print(f"m_from_contrib={m_from_contrib}")
print(f"b_from_contrib={b_from_contrib}")
plt.plot(x_from_contrib, y_from_contrib, 'o')
plt.title("Data from FitContribution object")
plt.show()
# ------------------------------------------------------
# 6. Now we create a FitRecipe to hold the contribution.
recipe = FitRecipe()
recipe.addContribution(contribution)
recipe.show()
# ------------------------------------------------------
# Once again, we can access the data and equation parameters
# from recipe as well
y_from_recipe = recipe.linear_fit.profile.yobs
x_from_recipe = recipe.linear_fit.profile.xobs
m_from_recipe = recipe.linear_fit.m.value
b_from_recipe = recipe.linear_fit.b.value
# additionally we can obtain all values from the recipe
values_from_recipe = recipe.values

print(f"values_from_recipe={values_from_recipe}")
print(f"m_from_recipe={m_from_recipe}")
print(f"b_from_recipe={b_from_recipe}")
plt.plot(x_from_recipe, y_from_recipe, 'o')
plt.title("Data from FitRecipe object")
plt.show()

# ------------------------------------------------------
# To tell the engine what values we want to refine,
# We need to add the variables with the addVar() method.
recipe.addVar(recipe.linear_fit.m, 1)
recipe.addVar(recipe.linear_fit.b, 1)

recipe.show()

# ------------------------------------------------------
# 8. Now that the recipe has been created, we can fit our model
# to the data. This is done through the least_squares function from
# scipy.optimize. We need to provide a value for the function to minimize
# and the initial values for the variables to refine.
residual = recipe.residual
values = recipe.values
leastsq(residual, values)

# ------------------------------------------------------
# 9. Finally, we can obtain the results and plot them.
fit_results = FitResults(recipe)
print(fit_results)


rw = fit_results.rw
xobs = recipe.linear_fit.profile.xobs
ycalc = recipe.linear_fit.profile.ycalc
yobs = recipe.linear_fit.profile.yobs

plt.plot(xobs, yobs, 'o', label='Observed Data')
plt.plot(xobs, ycalc, '-', label='Fitted Line')
plt.title(f'Fitted Line: rw={rw:.4f}')
plt.legend()
plt.show()
# ------------------------------------------------------
# 10. To keep these steps organized, we will now wrap them
#     in functions. In general, when you create a function,
#     you want to make it do one thing only. So we will
#     create one function to make the recipe and one to plot
#     the results. diffpy.cmi.fit_tools has a function to
#     optimize the recipe, so we will use that to optimize 
#    the recipe

#     First, let's create the function to make
#     the linear recipe.
def make_linear_recipe(xobs, yobs):
    profile = Profile()
    profile.setObservedProfile(xobs, yobs)

    contribution = FitContribution("linear_fit")
    contribution.setProfile(profile)
    contribution.setEquation("m * x + b ")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.linear_fit.m, 0.0)
    recipe.addVar(recipe.linear_fit.b, 0.0)

    return recipe

# ------------------------------------------------------
# 11. The function to optimize the recipe is already
#     available in diffpy.cmi.fit_tools, so we will use that.
from diffpy.cmi.fit_tools import optimize_recipe

# ------------------------------------------------------
# 12. Now, let's create the function to plot the results.

def plot_recipe_and_print_results(recipe, figsize=(8,6)):
    results = FitResults(recipe)
    print(results)
    rw = results.rw
    for name, contrib in recipe._contributions.items():
        profile = contrib.profile
        xobs = profile.xobs
        yobs = profile.yobs
        ycalc = profile.ycalc
        plt.figure(figsize=figsize)
        plt.plot(xobs, yobs, 'o')
        plt.plot(xobs, ycalc, label=f"Rw={rw:.4f}")
        plt.legend()
        plt.show()
# ------------------------------------------------------
# 13. Finally we can wrap everything together in a final
#     function that takes a recipe, optimizes it, and plots
#     and prints the results. Typically, it is common practice
#     to refer to a function that runs all your functions as main().

def main(recipe):
    optimize_recipe(recipe)
    plot_recipe_and_print_results(recipe)

# ------------------------------------------------------
# 14. Now we can create the linear recipe, optimize it,
#     and plot the results.
linear_recipe = make_linear_recipe(xobs, yobs)
main(linear_recipe)



# ------------------------------------------------------
# 14. Great! Now that we have fit the data with a linear function,
#     let's try to fit it with a parabolic function.
#     We will create a new function to make the parabolic recipe.
def make_parabolic_recipe(xobs, yobs):
    profile = Profile()
    profile.setObservedProfile(xobs, yobs)

    contribution = FitContribution("parabolic_fit")
    contribution.setProfile(profile)
    contribution.setEquation("a * x**2 + b * x + c ")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    recipe.addVar(recipe.parabolic_fit.a, 1)
    recipe.addVar(recipe.parabolic_fit.b, 1)
    recipe.addVar(recipe.parabolic_fit.c, 1)

    return recipe

# ------------------------------------------------------
# 15. Now we can create the parabolic recipe, optimize it,
#     and plot the results.

parabolic_recipe = make_parabolic_recipe(xobs, yobs)
main(parabolic_recipe)

# ------------------------------------------------------
# 16. Finally, let's try to fit the data with a cubic function.
#     The cubic function is our ground truth, so we expect
#     this fit to be very good.

def make_cubic_recipe(xobs, yobs):
    profile = Profile()
    profile.setObservedProfile(xobs, yobs)

    contribution = FitContribution("cubic_fit")
    contribution.setProfile(profile)
    contribution.setEquation("a * x**3 + b * x**2 + c * x + d ")

    recipe = FitRecipe()
    recipe.addContribution(contribution)
    cubic_variables = ["a", "b", "c", "d"]
    for var in cubic_variables:
        recipe.addVar(getattr(recipe.cubic_fit, var), 1)
    return recipe

# ------------------------------------------------------
# 17. Now we can create the cubic recipe, optimize it,
#     and plot the results.

cubic_recipe = make_cubic_recipe(xobs, yobs)
main(cubic_recipe)
