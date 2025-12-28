"""Plotting semivariograms"""

# MIT License

# Copyright (c) 2025 Guillaume Rongier

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import numpy as np
import matplotlib.pyplot as plt


################################################################################
# Variogram


def plot_semivariogram(
    direction,
    exp_vario,
    vario_model_values=None,
    distance_unit=None,
    ncols=4,
    basesize=(9, 2.5),
):
    """
    Plot an experimental semivariogram and a semivariogram model if available.

    Parameters
    ----------
    direction : Direction or array-like of shape (n_dir,)
        Direction(s) along which the experimental semi variogram was computed.
    exp_vario : pd.DataFrame
        Experimental semivariogram.
    vario_model_values : pd.DataFrame, default=None
        Values of the semivariogram model for each direction.
    distance_unit : str, default=None
        Unit for the distance.
    ncols : int, default=4
        Number of columns in the plot.
    basesize : array-like of shape (2,), default=(9, 2.5)
        Size of a row of the figure in inches.
    """
    if isinstance(direction, (tuple, list)) == False:
        direction = [direction]
    if distance_unit is not None:
        distance_unit = " (" + distance_unit + ")"
    else:
        distance_unit = ""

    # Deduce the number of rows based on the number of variogram directions
    nrows = int(np.ceil(len(direction) / ncols))
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharey=True,
        figsize=(basesize[0], basesize[1] * nrows),
        layout="constrained",
    )
    if axs.ndim == 1:
        axs = axs[np.newaxis]

    for j in range(nrows):
        for i in range(ncols):
            # Check to plot only if there is a direction left to plot
            k = j * ncols + i
            if k < len(direction):
                # Add the direction as title
                axs[j, i].set_title(
                    "Direction "
                    + "{:g}".format(direction[k].azimuth)
                    + " | Dip "
                    + "{:g}".format(direction[k].dip)
                )
                # Plot the experimental variogram for the direction
                axs[j, i].scatter(
                    exp_vario.loc[exp_vario["Direction"] == k + 1, "Distance"],
                    exp_vario.loc[exp_vario["Direction"] == k + 1, "Value"],
                    s=10,
                )
                if vario_model_values is not None:
                    # Plot the corresponding values for the variogram model
                    is_valid = (
                        vario_model_values["Azimuth"] == direction[k].azimuth
                    ) & (vario_model_values["Dip"] == direction[k].dip)
                    axs[j, i].plot(
                        vario_model_values.loc[is_valid, "Distance"],
                        vario_model_values.loc[is_valid, "Semivariance"],
                        c="tab:orange",
                    )
                # Add some legend elements
                axs[j, i].set_xlabel("Distance" + distance_unit)
                if i == 0:
                    axs[j, i].set_ylabel("Semivariance")
            else:
                axs[j, i].axis("off")

    plt.show()
