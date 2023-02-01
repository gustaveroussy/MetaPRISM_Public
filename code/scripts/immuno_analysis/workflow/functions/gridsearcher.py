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

Extended gridsearch with options to save best estimators and plot learning curves.
"""

import joblib
import os
import numpy as np
import pandas as pd
pd.options.display.float_format = '{:,.3f}'.format
import yaml

from sklearn.base import BaseEstimator, ClassifierMixin, is_classifier, is_regressor, clone
from sklearn.metrics import classification_report, confusion_matrix

class GridSearchCVExtended(BaseEstimator, ClassifierMixin):
    """
    GridSearchCVExtends extends the GridSearchCV class so that estimators for the best parameters fit on each split are
    saved as well as the cross-validation scores.
    """
    def __init__(self, grid, save_best_estimators_cv=False, save_best_estimator_train=True, save_best_params=True,
                 save_best_cv_scores=True, save_folder="models", verbose=True):
        """
        Parameters
        ----------
        grid: a GridSearchCV object
        save_best_estimators_cv: should the best estimators on each split be saved?
        save_best_estimator_train: should the best estimator refit on the train be saved?
        save_best_params: should the best parameters of the grid be saved?
        save_best_cv_scores: should the best cross-validation scores be saved?
        save_folder: estimators will be saved in subfolder `estimators` while scores will be saved in subfolder `scores`
        verbose: should intermediate messages be printed?
        """
        self.grid = grid
        self.save_best_estimators_cv = save_best_estimators_cv
        self.save_best_estimator_train = save_best_estimator_train
        self.save_best_cv_scores = save_best_cv_scores
        self.save_best_params = save_best_params
        self.save_folder = save_folder
        self.verbose = verbose

    def _print_clf_report(self, y_true, y_pred, labels):
        np_cm = confusion_matrix(y_true=y_true, y_pred=y_pred, labels=labels)
        df_cm = pd.DataFrame(np_cm)
        df_cm.index = pd.MultiIndex.from_product([["True"], labels])
        df_cm.columns = pd.MultiIndex.from_product([["Predicted"], labels])

        print("Confusion matrix", flush=True)
        print(df_cm, flush=True)
        print("\nClassification report", flush=True)
        print(classification_report(y_true=y_true, y_pred=y_pred, digits=3, zero_division=0), flush=True)

    def _save_estimator(self, estimator, filename):
        folder = os.path.join(self.save_folder, "estimators")
        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)
            if self.verbose:
                print("created folder %s" % folder, flush=True)

        if not filename.endswith(".joblib"):
            filename += ".joblib"

        filepath = os.path.join(folder, filename)
        with open(filepath, "wb") as file:
            joblib.dump(estimator, file)
        if self.verbose:
            print("saved estimator at %s" % filepath, flush=True)


    def _prepare_save_folder(self):
        folder = os.path.join(self.save_folder, "best")
        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)
            if self.verbose:
                print("\ncreated folder %s" % folder, flush=True)

        return folder

    def _save_scores(self, scores, filename):
        folder = self._prepare_save_folder()

        if not filename.endswith(".tsv"):
            filename += ".tsv"

        filepath = os.path.join(folder, filename)
        scores.to_csv(filepath, sep="\t", index=False)
        if self.verbose:
            print("saved scores at %s" % filepath, flush=True)


    def _save_params(self, params, filename):
        folder = self._prepare_save_folder()

        if not filename.endswith(".yaml"):
            filename += ".yaml"

        filepath = os.path.join(folder, filename)
        with open(filepath, "w") as file:
            yaml.dump(params, file, default_flow_style=False)

        if self.verbose:
            print("saved params at %s" % filepath, flush=True)

    def _convert_to_python_types(self, dt):
        dt_transformed = dt.copy()
        for key, val in dt.items():
            if np.issubdtype(type(val), np.inexact):
                dt_transformed[key] = float(val)
            elif np.issubdtype(type(val), np.integer):
                dt_transformed[key] = int(val)
            elif np.issubdtype(type(val), np.character):
                dt_transformed[key] = str(val)

        return dt_transformed

    def fit(self, X, y):
        self.best_estimators_ = {}
        self.cv_scores_ = pd.DataFrame()

        if isinstance(X, pd.DataFrame):
            X = X.values
        if isinstance(y, pd.Series):
            y = y.values

        self.grid = self.grid.fit(X,y)
        self.best_estimators_["train"] = self.grid.best_estimator_

        if self.verbose:
            print("\nParam grid provided to GridSearchCV:\n", flush=True)
            for key, val in self.grid.param_grid.items():
                print("\t - %s: %s" % (key, val), flush=True)
            print("\n")

        if self.save_best_estimator_train:
            self._save_estimator(self.best_estimators_["train"], "estimator_train")

        for i, (idx_train, idx_test) in enumerate(self.grid.cv.split(X,y)):
            split_name = "split%d" % i

            if self.verbose:
                print("\nSplit [%d/%d]" % (i+1, self.grid.cv.n_splits), flush=True)
                print("*"*60, flush=True)

            # Fit the estimator with best parameters on the split
            # NB: "Clone does a deep copy of the model in an estimator
            # without actually copying attached data. It yields a new estimator
            # with the same parameters that has not been fit on any data."
            best_estimator = clone(self.grid.best_estimator_)
            best_estimator = best_estimator.fit(X[idx_train], y[idx_train])
            self.best_estimators_[split_name] = best_estimator

            if self.save_best_estimators_cv:
                self._save_estimator(self.best_estimators_[split_name], "estimator_%s" % split_name)

            for idx, name in zip([idx_train, idx_test], ["train", "test"]):
                row_cv_scores = {"split": split_name, "type": name}

                y_true = y[idx]
                y_pred = best_estimator.predict(X[idx])

                if self.verbose:
                    print("\n\nCross-validation report %s" % name, flush=True)
                    print("-" * 40, flush=True)
                    if is_classifier(best_estimator):
                        self._print_clf_report(y_true=y_true, y_pred=y_pred, labels=np.unique(y))
                    elif is_regressor(best_estimator):
                        print("regression report NOT IMPLEMENTED yet")
                    else:
                        print("UNKNOWN estimator type")

                score_summary = []
                for score_name, scorer in self.grid.scoring.items():
                    row_cv_scores[score_name] = scorer(estimator=best_estimator,
                                                       X=X[idx],
                                                       y_true=y_true)
                    score_summary.append("%s: %.3g" % (score_name, row_cv_scores[score_name]))

                if self.verbose:
                    print("Summary %s" % name, flush=True)
                    print(" | ".join(sorted(score_summary)), flush=True)

                self.cv_scores_ = self.cv_scores_.append(row_cv_scores, ignore_index=True)

        if self.save_best_cv_scores:
            self._save_scores(self.cv_scores_, "best_cv_scores")

        if self.save_best_params:
            best_params = self._convert_to_python_types(self.grid.best_params_)
            self._save_params(best_params, "best_params")

        if self.verbose:
            print("\n")
            print("=" * 80)
            print("Global results\n")

            print("Best params:", flush=True)
            [print("\t%s: %s" % (key,value), flush=True) for key,value in self.grid.best_params_.items()]
            print("Best score: %.5g\n" % self.grid.best_score_, flush=True)

            print("CV scores:")

            print(self.cv_scores_.to_string(na_rep=""))

        return self

    def predict(self, X):
        return self.grid.predict(X)

    def predict_proba(self, X):
        return self.grid.predict_proba(X)

    def fit_predict(self, X,y):
        return self.fit(X,y).predict(X)

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self
