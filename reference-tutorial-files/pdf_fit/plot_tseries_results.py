from tseries_fit import (
    get_sorted_gr_files,
    get_temperature_from_filename,
    DATA_DIR,
)
import matplotlib.pyplot as plt

# Consistent colors for each element across both phases
ELEMENT_COLORS = {
    "Sr": "C0",  # Blue
    "Fe": "C1",  # Orange
    "As": "C2",  # Green
}

# --------------------------------------------------------------

from pathlib import Path


def parse_res_file(results_path, u11_mapping):
    """Parse a .res file to extract Rw and U11 values.

    Parameters
    ----------
    results_path : str or Path
        Path to the .res file.
    u11_mapping : dict[str, str]
        Mapping from element name to atom index, e.g. {"As": "0", "Fe": "4", "Sr": "12"}

    Returns
    -------
    rw : float
    u11_vals : dict[str, float]
        Keys are element names, values are U11 values.
    """
    results_path = Path(results_path)
    rw = None
    u11_vals = {}
    with open(results_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Rw"):
                parts = line.split()
                if len(parts) >= 2:
                    rw = float(parts[1])
            else:
                for element, idx in u11_mapping.items():
                    if line.startswith(f"U11_{idx}"):
                        u11_vals[element] = float(line.split()[1])
    if rw is None:
        raise ValueError(f"Rw not found in {results_path}")
    for element in u11_mapping:
        if element not in u11_vals:
            raise ValueError(
                f"U11_{u11_mapping[element]} not found in {results_path}"
            )
    return rw, u11_vals


def get_fit_results_orthorhombic(gr_file):
    """Returns temp, rw, u11_vals for orthorhombic phase by parsing the
    saved .res file."""
    temp = get_temperature_from_filename(gr_file)
    results_path = DATA_DIR / "fit-results" / f"ortho_{temp}K.res"
    u11_mapping = {"As": "0", "Fe": "4", "Sr": "12"}
    rw, u11_vals = parse_res_file(results_path, u11_mapping)
    return temp, rw, u11_vals


def get_fit_results_tetragonal(gr_file):
    """Returns temp, rw, u11_vals for tetragonal phase by parsing the
    saved .res file."""
    temp = get_temperature_from_filename(gr_file)
    results_path = DATA_DIR / "fit-results" / f"tet_{temp}K.res"
    u11_mapping = {"Sr": "0", "Fe": "2", "As": "6"}
    rw, u11_vals = parse_res_file(results_path, u11_mapping)
    return temp, rw, u11_vals


# --------------------------------------------------------------
# Now lets write our plotting function


def plot_all_fit_results(gr_files):
    """Create a 2x2 plot of Rw and U11 vs Temperature for both phases.

    Top row: Rw (Orthorhombic left, Tetragonal right)
    Bottom row: U11 (Orthorhombic left, Tetragonal right)
    """
    # Collect results for orthorhombic
    ortho_temps, ortho_rw, ortho_u11 = [], [], {"As": [], "Fe": [], "Sr": []}
    for gr_file in gr_files:
        temp, rw, u11 = get_fit_results_orthorhombic(gr_file)
        ortho_temps.append(temp)
        ortho_rw.append(rw)
        for element in u11:
            ortho_u11[element].append(u11[element])

    # Collect results for tetragonal
    tet_temps, tet_rw, tet_u11 = [], [], {"Sr": [], "Fe": [], "As": []}
    for gr_file in gr_files:
        temp, rw, u11 = get_fit_results_tetragonal(gr_file)
        tet_temps.append(temp)
        tet_rw.append(rw)
        for element in u11:
            tet_u11[element].append(u11[element])

    # Sort temperatures (in case files are out of order)
    ortho_sorted = sorted(
        zip(
            ortho_temps,
            ortho_rw,
            ortho_u11["As"],
            ortho_u11["Fe"],
            ortho_u11["Sr"],
        )
    )
    tet_sorted = sorted(
        zip(tet_temps, tet_rw, tet_u11["Sr"], tet_u11["Fe"], tet_u11["As"])
    )

    ortho_temps, ortho_rw, ortho_u11_As, ortho_u11_Fe, ortho_u11_Sr = zip(
        *ortho_sorted
    )
    tet_temps, tet_rw, tet_u11_Sr, tet_u11_Fe, tet_u11_As = zip(*tet_sorted)

    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)

    # ---- Rw plots ----
    axes[0, 0].plot(ortho_temps, ortho_rw, "o-")
    axes[0, 0].set_title("Orthorhombic Phase")
    axes[0, 0].set_ylabel("Rw")
    axes[0, 0].grid(alpha=0.3)

    axes[0, 1].plot(tet_temps, tet_rw, "o-")
    axes[0, 1].set_title("Tetragonal Phase")
    axes[0, 1].set_ylabel("Rw")
    axes[0, 1].grid(alpha=0.3)

    # ---- U11 plots ----
    # Orthorhombic (element order: Sr, Fe, As)
    for element, vals in zip(
        ["Sr", "Fe", "As"], [ortho_u11_Sr, ortho_u11_Fe, ortho_u11_As]
    ):
        axes[1, 0].plot(
            ortho_temps,
            vals,
            "o-",
            color=ELEMENT_COLORS[element],
            label=element,
        )
    axes[1, 0].set_xlabel("Temperature (K)")
    axes[1, 0].set_ylabel("U11")
    axes[1, 0].legend()
    axes[1, 0].grid(alpha=0.3)

    # Tetragonal (element order: Sr, Fe, As)
    for element, vals in zip(
        ["Sr", "Fe", "As"], [tet_u11_Sr, tet_u11_Fe, tet_u11_As]
    ):
        axes[1, 1].plot(
            tet_temps, vals, "o-", color=ELEMENT_COLORS[element], label=element
        )
    axes[1, 1].set_xlabel("Temperature (K)")
    axes[1, 1].set_ylabel("U11")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


def main():
    gr_files = get_sorted_gr_files(DATA_DIR)
    plot_all_fit_results(gr_files)


if __name__ == "__main__":
    main()


# def get_fit_results_orthorhombic(gr_file):
#     """
#     Gets the fit results from of a given gr_file for the orthorhombic phase fit.

#     Returns
#     -------
#     temp : int
#     rw : float
#     u11_vals : dict[str, float]
#         Keys are atom indices as strings, e.g. {"0": val, "4": val, "12": val}
#     """
#     # Create a recipe for this temperature
#     dummy_recipe = make_pdf_recipe(gr_file, ORTHORHOMBIC_CIF)
#     # Extract temperature from filename
#     temp = get_temperature_from_filename(gr_file)
#     # Path to saved results
#     results_dir = DATA_DIR / "fit-results"
#     phase = "ortho"
#     results_path = results_dir / f"{phase}_{temp}K.res"
#     # Initialize recipe variables from results file
#     initializeRecipe(dummy_recipe, str(results_path))
#     # Wrap in FitResults (not strictly needed for U11, but keeps API consistent)
#     results = FitResults(dummy_recipe)
#     # Extract Rw
#     rw = results.rw
#     # Extract U11 values with element labels as keys
#     u11_indices = {"As": "0", "Fe": "4", "Sr": "12"}
#     u11_vals = {}
#     for element, idx in u11_indices.items():
#         u11_vals[element] = getattr(dummy_recipe, f"U11_{idx}").value
#     return temp, rw, u11_vals


# def get_fit_results_tetragonal(gr_file):
#     """
#     Gets the fit results from of a given gr_file for the orthorhombic phase fit.

#     Returns
#     -------
#     temp : int
#     rw : float
#     u11_vals : dict[str, float]
#         Keys are element names ("Sr", "Fe", "As")
#     """
#     # Create a recipe for this temperature
#     dummy_recipe = make_pdf_recipe(gr_file, str(TETRAGONAL_CIF))
#     # Extract temperature from filename
#     temp = get_temperature_from_filename(gr_file)
#     # Path to saved results
#     results_dir = DATA_DIR / "fit-results"
#     phase = "tet"
#     results_path = results_dir / f"{phase}_{temp}K.res"
#     # Initialize recipe variables from results file
#     initializeRecipe(dummy_recipe, str(results_path))
#     # Wrap in FitResults
#     results = FitResults(dummy_recipe)
#     # Extract Rw
#     rw = results.rw
#     # Extract U11 values with element labels as keys
#     u11_indices = {"Sr": "0", "Fe": "2", "As": "6"}
#     u11_vals = {}
#     for element, idx in u11_indices.items():
#         u11_vals[element] = getattr(dummy_recipe, f"U11_{idx}").value
#     return temp, rw, u11_vals
