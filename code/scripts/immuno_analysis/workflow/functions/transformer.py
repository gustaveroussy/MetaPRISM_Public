# -*- coding: utf-8 -*-
"""
@modified: Jan 06 2021
@created: Jan 05 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France


Transformer sklearn classes (feature selectors).
"""

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler

# Feature selection/extraction =========================================================================================

def mean_abs_dev(X):
    return np.mean(np.abs(X-np.mean(X)))

class ScalerPandas(BaseEstimator, TransformerMixin):
    def __init__(self, center_func=None, scale_func=None, clip=None):
        self.center_func = center_func
        self.scale_func = scale_func
        self.clip = clip

    def fit(self, X, y=None):
        if isinstance(X, pd.DataFrame):
            self._X_columns = X.columns
            X_t = X
        else:
            X_t = pd.DataFrame(X)

        if self.center_func is None:
            self._centers = pd.Series(0, index=X_t.columns)
        else:
            self._centers = X_t.apply(self.center_func, axis=0)

        if self.scale_func is None:
            self._scales = pd.Series(1, index=X_t.columns)
        else:
            self._scales = X_t.apply(self.scale_func, axis=0)

        return self

    def transform(self, X):
        if isinstance(X, pd.DataFrame):
            X_t = X
        else:
            X_t = pd.DataFrame(X, index=self._X_columns)

        X_t = (X_t - self._centers) / self._scales

        if self.clip is not None:
            X_t = X_t.clip(-self.clip, self.clip)

        if isinstance(X, pd.DataFrame):
            return pd.DataFrame(X_t, index=X.index)
        else:
            return X_t

    def fit_transform(self, X, y=None):
        return self.fit(X,y).transform(X)


class SubpartsScalerPandas(BaseEstimator, TransformerMixin):
    """
    Performs a scaling using independent ScalerPandas on independent subparts of the intput data.
    Each subpart consist of a set of rows sharing an identical label in the 2nd and last level of the multi-index of
    the input.
    """
    def __init__(self, center_func=None, scale_func=None, clip=None):
        self.center_func = center_func
        self.scale_func = scale_func
        self.clip = clip

    def fit(self, X, y=None):
        self._labels_subparts = X.index.get_level_values(-1).unique().tolist()
        self._scaler_subparts = {}
        idx = pd.IndexSlice
        for label_subpart in self._labels_subparts:
            X_subpart = X.loc[idx[:,label_subpart],:]
            scaler_subpart = ScalerPandas(center_func=self.center_func, scale_func=self.scale_func, clip=self.clip)
            scaler_subpart = scaler_subpart.fit(X_subpart)
            self._scaler_subparts[label_subpart] = scaler_subpart.fit(X_subpart)

        return self

    def transform(self, X):
        X_subparts = []
        idx = pd.IndexSlice
        for label_subpart in X.index.get_level_values(-1).unique().tolist():
            X_subpart = X.loc[idx[:,label_subpart],:].copy()
            scaler_subpart = self._scaler_subparts[label_subpart]
            X_subpart = scaler_subpart.transform(X_subpart)
            X_subparts.append(X_subpart)

        X_t = pd.concat(X_subparts, axis=0)
        X_t = X_t.loc[X.index]
        X_t.index = X_t.index.droplevel(-1)

        return X_t

    def fit_transform(self, X, y=None):
        return self.fit(X,y).transform(X)
