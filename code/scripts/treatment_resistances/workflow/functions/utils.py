import numpy as np
import pandas as pd
from pathlib import Path
import re
import sys

import statsmodels.api as sm
import statsmodels.sandbox.stats.multicomp  as  mp


from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import RFE
from sklearn.model_selection import GridSearchCV, StratifiedKFold

from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, roc_auc_score
from sklearn.metrics import make_scorer

sys.path.append("workflow/functions")
from getter import get_estimator, get_param_grid
from gridsearcher import GridSearchCVExtended
from transformer import CustomFold, LassoFeatureSelector


def load_XY(outcome_var, outcome_rgx, dat_table, cov_table):
    df_dat = pd.read_table(dat_table)
    df_cov = pd.read_table(cov_table)

    covs_codes = df_cov["Code"]

    X = df_dat[covs_codes]
    Y = df_dat[outcome_var].apply(lambda x: int(bool(re.search(outcome_rgx, x))))

    return X,Y


def get_kfold_cv(df_indices, run):
    df_indices_run = df_indices.loc[df_indices["Run"]==run]
    df_indices_run = df_indices_run.rename(columns={"Split": "Fold"})
    kfold_cv = CustomFold(df_indices=df_indices_run)
    return kfold_cv


def split_train_test(X, Y, df_indices, b):
    df_indices_b = df_indices.loc[df_indices["Bootstrap"]==b]
    indices_train = df_indices_b.loc[df_indices_b["Split"]=="Train"].Index.values
    indices_test = df_indices_b.loc[df_indices_b["Split"]=="Test"].Index.values
    X_train, Y_train = X.loc[indices_train], Y.loc[indices_train]
    X_test, Y_test = X.loc[indices_test], Y.loc[indices_test]

    df_indices_b_cv = df_indices_b.loc[~df_indices_b["Split"].isin(["Train", "Test"])]
    df_indices_b_cv = df_indices_b_cv.rename(columns={"Split": "Fold"})
    kfold_cv_train = CustomFold(df_indices=df_indices_b_cv)
    return X_train, Y_train, X_test, Y_test, kfold_cv_train


def select_features(X, Y, args, kfold_cv=None):
    covariates = list(X.columns)
    df_selection = pd.DataFrame({"Covariate": covariates})
    df_selection["Selection"] = "In"

    if args.selection.endswith("_fdr"):
        if args.model not in ["LogisticRegressionStandard"]:
            raise ValueError("-ERROR: features selection method %s not applicable for model %s" % (args.selection,
                                                                                                   args.model))
        else:
            alpha = float(args.selection.split("%_")[0])/100
            get_pval = lambda cov: sm.Logit(Y, sm.add_constant(X[[cov]])).fit(method="bfgs").pvalues[cov]
            pvalues = Parallel(n_jobs=args.threads)(delayed(get_pval)(cov) for cov in covariates)
            mask_in, qvalues, _, _ = mp.multipletests(pvalues, method="fdr_bh", alpha=alpha)
            df_selection["p-value"] = pvalues
            df_selection["q-value"] = qvalues
    elif "_rfe" in args.selection:
        if "%_rfe" in args.selection:
            n_features_to_select = float(args.selection.split("%_")[0])/100
        elif "n_rfe" in args.selection:
            n_features_to_select = int(args.selection.split("n_")[0])

        if args.selection.endswith("_rfe"):
            estimator = get_estimator(args.model)
        else:
            model_rfe = args.selection.split("_rfe_")[1]
            estimator = get_estimator(model_rfe)

        selector = RFE(estimator=get_estimator(args.model),
                       n_features_to_select=n_features_to_select)
        selector = selector.fit(X,Y)
        mask_in = selector.support_
    elif args.selection == "lasso":
        selector = GridSearchCV(
            estimator=get_estimator("LogisticRegression"),
            param_grid=get_param_grid("LogisticRegressionLasso", args.models_grids),
            cv=kfold_cv,
            scoring="roc_auc",
            return_train_score=True,
            verbose=1,
            n_jobs=args.threads)

        selector = selector.fit(X,Y)
        mask_in = selector.best_estimator_.coef_.flatten() != 0
    elif args.selection=="None":
        mask_in = pd.Series(True, index=df_selection.index)
    else:
        raise ValueError("-ERROR: unsupported feature selection method. Choose x%/n_ref or x%_fdr (linear model only)")

    if sum(mask_in)==0:
        print("-WARNING: 0 feature pass the selection criteria %s. All features are retained." % args.selection)
    else:
        df_selection.loc[~mask_in, "Selection"] = "Out"

    return df_selection


def build_grid(args, estimator, param_grid, kfold_cv, verbose=1, n_jobs=1):
    scoring = {"accuracy": make_scorer(accuracy_score),
               "recall": make_scorer(recall_score, average="macro", zero_division=0),
               "precision": make_scorer(precision_score, average="macro", zero_division=0),
               "f1_weighted": make_scorer(f1_score, average="weighted", zero_division=0),
               "f1_macro": make_scorer(f1_score, average="macro", zero_division=0),
               "roc_auc": make_scorer(roc_auc_score)}

    grid = GridSearchCVExtended(
        grid = GridSearchCV(
            estimator=estimator,
            param_grid=param_grid,
            cv=kfold_cv,
            scoring=scoring,
            refit="roc_auc",
            return_train_score=True,
            verbose=verbose,
            n_jobs=n_jobs
        ),
        cv_tests=[],
        save_best_estimators_cv=False,
        save_best_estimator_train=False,
        save_best_cv_scores=False,
        save_best_params=False,
        save_folder=str(Path(args.output_pipe).parent.absolute()),
        verbose=verbose)

    return grid


def build_pipe(args, grid, kfold_cv=None, verbose=1, n_jobs=1):
    if args.selection == "lasso":
        param_grid = get_param_grid("LogisticRegressionLasso", args.selections_grids)
        pipe = Pipeline([
            ('selector', LassoFeatureSelector(kfold_cv=kfold_cv, param_grid=param_grid,
                                              n_jobs=n_jobs, scoring="roc_auc", verbose=verbose)),
            ('gridsearch', grid),
        ])
    else:
        raise NotImplementedError("-ERROR: selection %s not implemented yet" % (args.selection))

    return pipe


def extract_metrics_train(pipe, args):
    grid = pipe.steps[-1][1]
    df_met = grid.cv_scores_
    df_met["Selection"] = args.selection
    df_met["Model"] = args.model
    df_met["Features"] = args.features

    return df_met


def extract_metrics_test(pipe, args, X_test, Y_test, grid_run="Test"):
    grid = pipe.steps[-1][1]
    scoring = grid.grid.scoring
    scores = {}
    for name, scorer in scoring.items():
        scores[name] = [scorer(pipe, X_test, Y_test)]

    df_met = pd.DataFrame(scores)
    df_met["Grid_Run"] = grid_run
    df_met["Selection"] = args.selection
    df_met["Model"] = args.model
    df_met["Features"] = args.features

    return df_met



def extract_coefficients_train(pipe, args, include_internal_splits=False, include_global=True):
    df_cov = pd.DataFrame()
    grid = pipe.steps[-1][1]
    model = grid.best_estimators_["train"]

    selector = pipe.steps[0][1]
    covariates_in = selector.covariates_in_
    covariates_out = selector.covariates_out_

    if hasattr(model, "coef_"):
        attr = "coef_"
        attr_name = "Coefficient"
    elif hasattr(model, "feature_importances_"):
        attr = "feature_importances_"
        attr_name = "Feature_Importance"
    else:
        attr = None
        attr_name = None

    if attr is not None:
        for grid_run, model in grid.best_estimators_.items():
            if not include_internal_splits and grid_run.startswith("split"):
                pass
            elif not include_global and grid_run=="train":
                pass
            else:
                df_cov_run = pd.DataFrame({"Covariate": covariates_in, attr_name: list(getattr(model, attr).flatten())})
                df_cov_run["Selected"] = 1
                df_cov_run["Selection"] = args.selection
                df_cov_run["Grid_Run"] = grid_run
                df_cov_run["Model"] = args.model
                df_cov_run["Features"] = args.features
                for param_name, param_value in grid.grid.best_params_.items():
                    df_cov_run[param_name] = param_value

                # if any covariate was excluded in the feature selection step, record it here as 0
                if len(covariates_out) > 0:
                    df_cov_out = pd.DataFrame({"Covariate": covariates_out})
                    df_cov_out[attr_name] = 0
                    df_cov_out["Selected"] = 0
                    df_cov_out["Selection"] = args.selection
                    df_cov_out["Grid_Run"] = grid_run
                    df_cov_out["Model"] = args.model
                    df_cov_out["Features"] = args.features
                    df_cov_run = pd.concat((df_cov_run, df_cov_out))

                df_cov = pd.concat((df_cov, df_cov_run), axis=0)

    return df_cov
