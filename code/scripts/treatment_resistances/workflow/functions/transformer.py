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

from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.linear_model import LogisticRegression

import warnings

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


class CustomFold():
    """
    Build an object that can return a list of train and test indices from the `split` method
    similarly to instances of StratifiedKFold or KFold from sklearn. This class allows to specify
    a custom sequence of indices stored in a pandas dataframe.
    """
    def __init__(self, df_indices):
        self.df_indices = df_indices
        self.n_splits = len(df_indices["Fold"].unique())

    def split(self, X=None, y=None, groups=None):
        folds = sorted(self.df_indices["Fold"].unique())
        indices = []
        for fold in folds:
            mask_fold = self.df_indices["Fold"]==fold
            train_index = self.df_indices.loc[~mask_fold, "Index"].unique()
            test_index = self.df_indices.loc[mask_fold, "Index"].unique()
            indices.append((train_index, test_index))
        return indices

    def get_n_splits(self, X=None, y=None, groups=None):
        return self.n_splits



class LassoFeatureSelector(BaseEstimator, TransformerMixin):
    def __init__(self, param_grid, kfold_cv=None, scoring="roc_auc", n_jobs=1, verbose=1):
        self.param_grid = param_grid
        self.kfold_cv = kfold_cv
        self.scoring = scoring
        self.n_jobs = n_jobs
        self.verbose = verbose

        if self.kfold_cv is None:
            self.kfold_cv = StratifiedKFold(n_splits=5, shuffle=False)


    def fit(self, X, y):
        covariates = list(X.columns)
        self.df_selection_ = pd.DataFrame({"Covariate": covariates})
        self.df_selection_["Selection"] = "In"

        self.selector_ = GridSearchCV(
            estimator=LogisticRegression(penalty="l1"),
            param_grid=self.param_grid,
            cv=self.kfold_cv,
            scoring=self.scoring,
            return_train_score=True,
            n_jobs=self.n_jobs,
            verbose=self.verbose)

        self.selector_ = self.selector_.fit(X,y)

        mask_in = self.selector_.best_estimator_.coef_.flatten() != 0
        if sum(mask_in)==0:
            warnings.warn("0 feature passed the selection by the Lasso. All features are retained.", UserWarning)
        else:
            self.df_selection_.loc[~mask_in, "Selection"] = "Out"

        self.covariates_out_ = self.df_selection_.loc[self.df_selection_["Selection"]=="Out", "Covariate"].tolist()
        self.covariates_in_ = self.df_selection_.loc[self.df_selection_["Selection"]=="In", "Covariate"].tolist()

        return self


    def transform(self, X):
        return X[self.covariates_in_].copy()

    def fit_transform(self, X, y):
        return self.fit(X,y).transform(X)
