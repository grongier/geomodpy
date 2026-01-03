"""Plotting grids"""

# MIT License

# Copyright (c) 2025 Guillaume Rongier

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import ipywidgets as widgets


################################################################################
# Slices through 3D structured grids


def plot_3d_slices(
    ds,
    var,
    i=None,
    j=None,
    k=None,
    x="X",
    y="Y",
    z="Z",
    min_aspect=1.0,
    plot_slices=True,
    cmap=None,
    use_log=False,
    norm=None,
    clabel=None,
    spatial_unit=None,
    figsize=None,
):
    """
    Plots slices from a variable of a 3D regular, structured grid in a Xarray
    DataSet or a DataArray along each axis.

    Parameters
    ----------
    ds : xr.DataSet or xr.DataArray of shape (nz, ny, nx)
        3D structured grid containing the variable to plot.
    var : str or int
        Name of the variable to plot.
    i : int, default=None
        Index of the slice perpendicular to the x axis to plot. If None, uses
        half the number of cells along the x axis.
    j : int, default=None
        Index of the slice perpendicular to the y axis to plot. If None, uses
        half the number of cells along the y axis.
    k : int, default=None
        Index of the slice perpendicular to the z axis to plot. If None, uses
        half the number of cells along the z axis.
    x : str, default='X'
        Name of the x-axis coordinates in `ds`.
    y : str, default='Y'
        Name of the x-axis coordinates in `ds`.
    z : str, default='Z'
        Name of the x-axis coordinates in `ds`.
    min_aspect : float, default=1.
        Minimum aspect factor for the vertical slices.
    plot_slices : bool, default=True
        If True, displays lines that shows where the slices are located.
    cmap : str or Colormap, default=None
        Colormap instance or registered colormap name used to map scalar data to
        colors.
    use_log : bool, default=False
        If True, use a log normalizer for the colorbar.
    norm : Normalize, default=None
        Normalization method used to scale scalar data to the [0, 1] range
        before mapping to colors using cmap. By default, a linear scaling is
        used, mapping the lowest value to 0 and the highest to 1.
    clabel : str, default=None
        Label of the colorbar.
    spatial_unit : str, default=None
        Unit for the spatial dimensions (x, y, z).
    figsize : array-like of shape (2,), default=None
        Size of the figure in inches.
    """
    fig = plt.figure(figsize=figsize)
    fig.canvas.header_visible = False

    if isinstance(ds, xr.DataArray) == False or isinstance(var, int):
        ds = ds[var]

    i = ds.shape[0] // 2 if i is None else i
    j = ds.shape[1] // 2 if j is None else j
    k = ds.shape[2] // 2 if k is None else k
    if norm is None:
        vmin = min(ds[i].min(), ds[:, j].min(), ds[..., k].min())
        vmax = max(ds[i].max(), ds[:, j].max(), ds[..., k].max())
        if vmax > 0.0 and vmin < 0.0:
            vmax = max(abs(vmin), vmax)
            vmin = -vmax
        if use_log == True and vmin > 0.0 and vmax > 0.0:
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if spatial_unit is not None:
        spatial_unit = " (" + spatial_unit + ")"
    else:
        spatial_unit = ""

    ax = fig.add_gridspec(top=0.75, right=0.75).subplots()
    ds[..., k].plot(x=x, y=y, ax=ax, norm=norm, cmap=cmap, add_colorbar=False)
    if plot_slices == True:
        if ds[x].ndim == 1 and ds[y].ndim == 1:
            ax.axvline(ds[x][i], c="k", ls="--", lw=0.75)
            ax.axhline(ds[y][j], c="k", ls="--", lw=0.75)
        else:
            ax.plot(ds[x][i, :, k], ds[y][i, :, k], c="k", ls="--", lw=0.75)
            ax.plot(ds[x][:, j, k], ds[y][:, j, k], c="k", ls="--", lw=0.75)
    ax.set(aspect="equal", title=None, xlabel=x + spatial_unit, ylabel=y + spatial_unit)

    aspect = (ds[z].max() - ds[z].min()) / (ds[y].max() - ds[y].min())
    aspect = aspect if aspect > min_aspect else min_aspect
    ratio = ds.shape[2] / ds.shape[1] if ds.shape[2] < ds.shape[1] else 1.0
    ax_slicey = ax.inset_axes([0, 1.0 + ratio * 0.1, 1, ratio * aspect], sharex=ax)
    ax_slicey.tick_params(axis="x", labelbottom=False)
    ds[:, j].plot(x=x, y=z, ax=ax_slicey, norm=norm, cmap=cmap, add_colorbar=False)
    if plot_slices == True:
        if ds[x].ndim == 1 and ds[z].ndim == 1:
            ax_slicey.axvline(ds[x][i], c="k", ls="--", lw=0.75)
            ax_slicey.axhline(ds[z][k], c="k", ls="--", lw=0.75)
        else:
            ax_slicey.plot(ds[x][i, j, :], ds[z][i, j, :], c="k", ls="--", lw=0.75)
            ax_slicey.plot(ds[x][:, j, k], ds[z][:, j, k], c="k", ls="--", lw=0.75)
    ax_slicey.set(title=None, xlabel=None, ylabel=z + spatial_unit)

    aspect = (ds[z].max() - ds[z].min()) / (ds[x].min() - ds[x].max())
    aspect = aspect if aspect > min_aspect else min_aspect
    ratio = ds.shape[2] / ds.shape[0] if ds.shape[2] < ds.shape[0] else 1.0
    ax_slicex = ax.inset_axes([1.0 + ratio * 0.1, 0, ratio * aspect, 1], sharey=ax)
    ax_slicex.tick_params(axis="y", labelleft=False)
    im = ds[i].T.plot(x=z, y=y, ax=ax_slicex, norm=norm, cmap=cmap, add_colorbar=False)
    if plot_slices == True:
        if ds[y].ndim == 1 and ds[z].ndim == 1:
            ax_slicex.axhline(ds[y][j], c="k", ls="--", lw=0.75)
            ax_slicex.axvline(ds[z][k], c="k", ls="--", lw=0.75)
        else:
            ax_slicex.plot(ds[z][i, j, :], ds[y][i, j, :], c="k", ls="--", lw=0.75)
            ax_slicex.plot(ds[z][i, :, k], ds[y][i, :, k], c="k", ls="--", lw=0.75)
    ax_slicex.set(title=None, xlabel=z + spatial_unit, ylabel=None)

    ax_cbar = ax.inset_axes([1.0 + ratio * (0.15 + aspect + 0.05), 0.25, 0.03, 0.5])
    fig.colorbar(im, cax=ax_cbar, label=clabel)

    plt.show()


def plot_3d_slices_interactive(
    ds,
    x="U",
    y="V",
    z="W",
    min_aspect=1.0,
    plot_slices=True,
    cmap=None,
    clabel=None,
    spatial_unit=None,
    figsize=None,
):
    """
    Interactively plots slices from a variable of a 3D regular, structured grid in
    a Xarray DataSet or a DataArray along each axis.

    Parameters
    ----------
    ds : xr.DataSet or xr.DataArray of shape (nz, ny, nx)
        3D regular, structured grid containing the variable to plot.
    x : str, default='U'
        Name of the x-axis coordinates in `ds`.
    y : str, default='V'
        Name of the x-axis coordinates in `ds`.
    z : str, default='W'
        Name of the x-axis coordinates in `ds`.
    min_aspect : float, default=1.
        Minimum aspect factor for the vertical slices.
    plot_slices : bool, default=True
        If True, displays lines that shows where the slices are located.
    cmap : str or Colormap, default=None
        Colormap instance or registered colormap name used to map scalar data to
        colors.
    clabel : str, default=None
        Label of the colorbar.
    spatial_unit : str, default=None
        Unit for the spatial dimensions (x, y, z).
    figsize : array-like of shape (2,), default=None
        Size of the figure in inches.

    Returns
    -------
    interactive_plot : Output
        Interactive plot.
    """
    fig = plt.figure(figsize=figsize)
    fig.canvas.header_visible = False

    var = list(ds.keys())
    i = ds[var[0]].shape[0] // 2
    j = ds[var[0]].shape[1] // 2
    k = ds[var[0]].shape[2] // 2
    vmin = min(ds[var[0]][i].min(), ds[var[0]][:, j].min(), ds[var[0]][..., k].min())
    vmax = max(ds[var[0]][i].max(), ds[var[0]][:, j].max(), ds[var[0]][..., k].max())
    if vmax > 0.0 and vmin < 0.0:
        vmax = max(abs(vmin), vmax)
        vmin = -vmax
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if spatial_unit is not None:
        spatial_unit = " (" + spatial_unit + ")"
    else:
        spatial_unit = ""

    args = {"continuous_update": False}
    wdg_var = widgets.Dropdown(
        options=var, description="Variable", value=var[0], **args
    )
    wdg_log = widgets.Checkbox(value=False, description="Log", **args)
    wdg_i = widgets.IntSlider(
        value=i, min=0, max=ds[var[0]].shape[0] - 1, step=1, description="i", **args
    )
    wdg_j = widgets.IntSlider(
        value=j, min=0, max=ds[var[0]].shape[1] - 1, step=1, description="j", **args
    )
    wdg_k = widgets.IntSlider(
        value=k, min=0, max=ds[var[0]].shape[2] - 1, step=1, description="k", **args
    )

    ax = fig.add_gridspec(top=0.75, right=0.75).subplots()
    im1 = ds[var[0]][..., k].plot.imshow(
        x=x, y=y, ax=ax, norm=norm, cmap=cmap, add_colorbar=False
    )
    if plot_slices == True:
        ln1 = ax.axvline(ds[x][i], c="k", ls="--", lw=0.75)
        ln2 = ax.axhline(ds[y][j], c="k", ls="--", lw=0.75)
    ax.set(aspect="equal", title=None, xlabel=x + spatial_unit, ylabel=y + spatial_unit)

    aspect = (ds[z].max() - ds[z].min()) / (ds[y].max() - ds[y].min())
    aspect = aspect if aspect > min_aspect else min_aspect
    ratio = (
        ds[var[0]].shape[2] / ds[var[0]].shape[1]
        if ds[var[0]].shape[2] < ds[var[0]].shape[1]
        else 1.0
    )
    ax_slicey = ax.inset_axes([0, 1.0 + ratio * 0.1, 1, ratio * aspect], sharex=ax)
    ax_slicey.tick_params(axis="x", labelbottom=False)
    im2 = ds[var[0]][:, j].plot.imshow(
        x=x, y=z, ax=ax_slicey, norm=norm, cmap=cmap, add_colorbar=False
    )
    if plot_slices == True:
        ln3 = ax_slicey.axvline(ds[x][i], c="k", ls="--", lw=0.75)
        ln4 = ax_slicey.axhline(ds[z][k], c="k", ls="--", lw=0.75)
    ax_slicey.set(title=None, xlabel=None, ylabel="z" + spatial_unit)

    aspect = (ds[z].max() - ds[z].min()) / (ds[x].min() - ds[x].max())
    aspect = aspect if aspect > min_aspect else min_aspect
    ratio = (
        ds[var[0]].shape[2] / ds[var[0]].shape[0]
        if ds[var[0]].shape[2] < ds[var[0]].shape[0]
        else 1.0
    )
    ax_slicex = ax.inset_axes([1.0 + ratio * 0.1, 0, ratio * aspect, 1], sharey=ax)
    ax_slicex.tick_params(axis="y", labelleft=False)
    im3 = ds[var[0]][i].T.plot.imshow(
        x=z, y=y, ax=ax_slicex, norm=norm, cmap=cmap, add_colorbar=False
    )
    if plot_slices == True:
        ln5 = ax_slicex.axhline(ds[y][j], c="k", ls="--", lw=0.75)
        ln6 = ax_slicex.axvline(ds[z][k], c="k", ls="--", lw=0.75)
    ax_slicex.set(title=None, xlabel=z + spatial_unit, ylabel=None)

    ax_cbar = ax.inset_axes([1.0 + ratio * (0.15 + aspect + 0.05), 0.25, 0.03, 0.5])
    fig.colorbar(im3, cax=ax_cbar, label=clabel)

    def update(var, i, j, k, log):
        vmin = min(ds[var][i].min(), ds[var][:, j].min(), ds[var][..., k].min())
        vmax = max(ds[var][i].max(), ds[var][:, j].max(), ds[var][..., k].max())
        if vmax > 0.0 and vmin < 0.0:
            vmax = max(abs(vmin), vmax)
            vmin = -vmax
        if log == True and vmin > 0.0 and vmax > 0.0:
            norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        im1.set_data(ds[var][..., k].T)
        im1.set_norm(norm)
        im2.set_data(ds[var][:, j].T)
        im2.set_norm(norm)
        im3.set_data(ds[var][i])
        im3.set_norm(norm)
        if plot_slices == True:
            ln1.set_xdata([ds[x][i]])
            ln2.set_ydata([ds[y][j]])
            ln3.set_xdata([ds[x][i]])
            ln4.set_ydata([ds[z][k]])
            ln5.set_ydata([ds[y][j]])
            ln6.set_xdata([ds[z][k]])

    out = widgets.interactive_output(
        update, {"var": wdg_var, "i": wdg_i, "j": wdg_j, "k": wdg_k, "log": wdg_log}
    )
    return widgets.VBox(
        [widgets.HBox([wdg_var, wdg_log]), widgets.HBox([wdg_i, wdg_j, wdg_k]), out]
    )
