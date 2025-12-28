"""Gaussian simulation"""

# MIT License

# Copyright (c) 2022-2025 Guillaume Rongier

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
from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.utils import check_random_state
from sklearn.utils.validation import validate_data

from ..analysis.variography import VarioModel
from ..utils import create_dataset


################################################################################
# Gaussian simulation


def build_coords_of_cartesian_grid(shape, spacing, origin):
    """
    Build coordinates of a regular structured grid from its shape, spacing, and
    origin.

    Parameters
    ----------
    shape : array-like of shape (n_dims,)
        Shape of the grid along each axis.
    spacing : array-like of shape (n_dims,)
        Cell size along the x, y, and z axes.
    origin : array-like of shape (n_dims,)
        Coordinates at the center of the cell at the bottom left of the grid.

    Returns
    -------
    coords : ndarray of shape (size_1, ..., size_n, n_dims)
        Coordinates of the grid cells.
    """
    return np.stack(
        np.meshgrid(
            *[
                np.linspace(o, o + l * (s - 1), s)
                for s, l, o in zip(shape, spacing, origin)
            ],
            indexing="ij",
        ),
        axis=-1,
    )


class GaussianSimulator(RegressorMixin, BaseEstimator):
    """
    Gaussian Simulation.

    Parameters
    ----------
    kernel : kernel instance from sklearn or VarioModel, default=VarioModel()
        The kernel specifying the covariance function or the semivariogram model.
    sampler : str, default='cholesky_decomposition'
        Sampling method to use, which can be:
            - 'cholesky_decomposition' for the Cholesky decomposition (default).
            - 'eigen_decomposition' for the eigen-decomposition, which is slower.
            - 'circulant_embedding' for the circulant embedding approach, which
              is faster but requires to sample or predict on a 3D regular
              structured grid.
    y_transformer : object, default=None
        Estimator object such as derived from `TransformerMixin`. The transformer is
        restricting y to be a numpy array.
    return_dataset : bool, default=True
        If true, returns a Xarray DataSet, otherwise returns a NumPy array.
    variable_name : string, default='Variable'
        Name of the simulated variable. Only used if `return_dataset` is true.

    Attributes
    ----------
    X_train_ : array-like of shape (n_samples, n_features)
        Feature vectors or other representations of training data (also
        required for prediction).
    y_train_ : array-like of shape (n_samples,)
        Target values in training data (also required for prediction).
    kernel_ : kernel instance from sklearn or VarioModel
        The kernel used for prediction.
    y_transformer_ : object
        Transformer used in :meth:`fit`, :meth:`sample_y`, and :meth:`predict`.

    Source
    ------
    https://github.com/driftingtides/hyvr/blob/master/hyvr/utils.py

    References
    ----------
    Alabert , F. (1987)
        The practice of fast conditional simulations through the LU decomposition of the covariance matrix
        https://doi.org/10.1007/BF00897191
    Dietrich, C.R., Newsam, G.N. (1993)
        A Fast and Exact Method for Multidimensional Gaussian Stochastic Simulations
        https://doi.org/10.1029/93WR01070
    """

    def __init__(
        self,
        kernel=VarioModel(),
        sampler="cholesky_decomposition",
        y_transformer=None,
        return_dataset=False,
        variable_name="Variable",
    ):

        self.kernel = kernel
        self.sampler = sampler
        self.y_transformer = y_transformer
        self.return_dataset = return_dataset
        self.variable_name = variable_name

    def fit(self, X, y):
        """
        Fit the Gaussian simulator.

        Contrary to scikit-learn's GaussianProcessRegressor, this function does
        not fit the kernel's parameters. It simply stores the conditioning data
        for sampling and prediction.

        Conditioning with the circulant embedding approach isn't implemented yet.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Feature vectors or other representations of training data.
        y : array-like of shape (n_samples,)
            Target values.

        Returns
        -------
        self : object
            GaussianSimulator class instance.
        """
        X, y = validate_data(self, X, y)

        self.X_train_ = X
        if self.y_transformer is not None:
            self.y_transformer_ = clone(self.y_transformer)
            self.y_train_ = self.y_transformer_.fit_transform(y.reshape(-1, 1)).squeeze(
                axis=1
            )
        else:
            self.y_train_ = y

        if (
            self.sampler == "cholesky_decomposition"
            or self.sampler == "eigen_decomposition"
        ):
            if isinstance(self.kernel, VarioModel) == True:
                self.kernel_ = self.kernel.covariance_matrix
            else:
                self.kernel_ = self.kernel
            self._C_data = self.kernel_(X)

        return self

    def _decompose(self, C):
        """
        Decompose a covariance matrix.
        """
        if self.sampler == "cholesky_decomposition":
            L = np.linalg.cholesky(C)
        else:
            # Compute the eigenvalues and eigenvectors
            evals, evecs = np.linalg.eigh(C)
            # Construct C, so C*C^t = cov
            L = np.dot(evecs, np.diag(np.sqrt(evals)))

        return L

    def _sample_y_lusim(self, X, n_samples, rng):
        """
        Draw samples using the decomposition of the covariance matrix.
        """
        if isinstance(X, dict) == True:
            _X = build_coords_of_cartesian_grid(X["shape"], X["spacing"], X["origin"])
            _X = _X.reshape(-1, len(X["shape"]))
            shape = X["shape"]
        else:
            _X = X
            shape = (len(_X),)
        _X = validate_data(self, _X, reset=False)

        if hasattr(self, "X_train_"):
            n_data = len(self.X_train_)
            C = np.empty((n_data + len(_X), n_data + len(_X)))
            C[:n_data, :n_data] = self._C_data
            C[n_data:, n_data:] = self.kernel_(_X)
            C[:n_data, n_data:] = self.kernel_(self.X_train_, _X)
            C[n_data:, :n_data] = C[:n_data, n_data:].T
        else:
            C = self.kernel_(_X)
        L = self._decompose(C)

        w = rng.normal(size=(len(_X), n_samples))
        if hasattr(self, "X_train_"):
            _w = np.linalg.solve(L[:n_data, :n_data], self.y_train_)
            _w = np.repeat(_w[:, None], n_samples, axis=-1)
            y_samples = L[n_data:] @ np.concatenate((_w, w), axis=0)
        else:
            y_samples = L @ w
        y_samples = y_samples.reshape(shape + (n_samples,))

        return y_samples

    def _build_embedding(self, shape, spacing):
        """
        Build the embedding for a 3D regular structured grid from its shape and spacing.
        """
        embedding = np.stack(
            np.meshgrid(
                *[
                    np.concatenate(
                        (np.linspace(0.0, l * (s - 1), s), np.linspace(-l * s, -l, s))
                    )
                    for s, l in zip(shape, spacing)
                ],
                indexing="ij",
            ),
            axis=-1,
        )
        embedding = self.kernel_(embedding.reshape(-1, 3), np.zeros((1, 3))).reshape(
            embedding.shape[:-1]
        )

        return np.fft.fftn(embedding) / embedding.size

    def _sample_y_circemb(self, shape, spacing, n_samples, rng):
        """
        Draw samples using the circulant embedding approach.
        """
        embedding = self._build_embedding(shape, spacing)

        noise = rng.randn(2, n_samples // 2, *embedding.shape)
        y_samples = (noise[0] + 1j * noise[1]) * np.sqrt((n_samples // 2) * embedding)
        y_samples = np.concatenate(
            (
                np.real(np.fft.ifftn(y_samples * embedding.size)),
                np.imag(np.fft.ifftn(y_samples * embedding.size)),
            )
        )
        y_samples = np.moveaxis(y_samples[:, : shape[0], : shape[1], : shape[2]], 0, -1)

        return y_samples

    def sample_y(self, X, n_samples=1, random_state=42):
        """
        Draw samples from the Gaussian simulator.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or dictionary
            Query points where the Gaussian simulation is run. If a dictionary,
            then the query points are located on a 3D regular structured grid
            defined by the following keys:
                - shape : int or array-like (x, y, z)
                    Shape of the grid along the x, y, and z axes.
                - spacing : float or array-like (x, y, z)
                    Cell size along the x, y, and z axes.
                - origin : float or array-like (x, y, z)
                    Coordinates at the center of the cell at the bottom left of
                    the grid along the x, y, and z axes.
        n_samples : int, default=100
            Number of samples drawn from the Gaussian simulator per query point.
        random_state : int, RandomState instance or None, default=42
            Determines random number generation to randomly draw samples.
            Pass an int for reproducible results across multiple function
            calls.

        Returns
        -------
        y_samples : ndarray of shape (n_samples_X, n_samples)
            Values of `n_samples` realizations evaluated at the query points.
        """
        rng = check_random_state(random_state)

        if hasattr(self, "kernel_") == False:
            if isinstance(self.kernel, VarioModel) == True:
                self.kernel_ = self.kernel.covariance_matrix
            else:
                self.kernel_ = self.kernel

        if self.sampler == "circulant_embedding" and isinstance(X, dict) == True:
            y_samples = self._sample_y_circemb(X["shape"], X["spacing"], n_samples, rng)
        else:
            y_samples = self._sample_y_lusim(X, n_samples, rng)

        if hasattr(self, "y_transformer_") == True:
            y_samples = self.y_transformer_.inverse_transform(
                y_samples.reshape(-1, 1)
            ).reshape(y_samples.shape)

        if self.return_dataset == True and isinstance(X, dict) == True:
            return create_dataset(
                np.moveaxis(y_samples, -1, 0),
                X["spacing"],
                X["origin"],
                n_samples,
                var_name=self.variable_name,
            )
        return y_samples

    def predict(self, X, return_std=False, n_samples=100, random_state=42):
        """
        Predict using the Gaussian simulator.

        We can also predict based on an unfitted model by using the Gaussian
        process prior. In addition to the mean of the predictive distribution,
        optionally also returns its standard deviation (`return_std=True`).

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or dictionary
            Query points where the Gaussian simulation is run. If a dictionary,
            then the query points are located on a 3D regular structured grid
            defined by the following keys:
                - shape : int or array-like (x, y, z)
                    Shape of the grid along the x, y, and z axes.
                - spacing : float or array-like (x, y, z)
                    Cell size along the x, y, and z axes.
                - origin : float or array-like (x, y, z)
                    Coordinates at the center of the cell at the bottom left of
                    the grid along the x, y, and z axes.
        return_std : bool, default=False
            If True, the standard-deviation of the predictive distribution at
            the query points is returned along with the mean.
        n_samples : int, default=100
            Number of samples drawn from the Gaussian simulator per query point.
        random_state : int, RandomState instance or None, default=42
            Determines random number generation to randomly draw samples.
            Pass an int for reproducible results across multiple function
            calls.

        Returns
        -------
        y_mean : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Mean of predictive distribution at query points.
        y_std : ndarray of shape (n_samples,) or (n_samples, n_targets), optional
            Standard deviation of predictive distribution at query points.
            Only returned when `return_std` is True.
        """
        y_samples = self.sample_y(X, n_samples=n_samples, random_state=random_state)

        if self.return_dataset == True and isinstance(X, dict) == True:
            y_pred = y_samples.mean("Realization")
            if return_std == True:
                y_pred["Standard_deviation"] = y_samples[self.variable_name].std(
                    "Realization"
                )
            return y_pred
        else:
            if return_std == False:
                return np.mean(y_samples, axis=-1)
            else:
                return np.mean(y_samples, axis=-1), np.std(y_samples, axis=-1)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.requires_fit = False
        return tags
