# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 20:29 2019

@author: Yoann Pradat

The present file implements usual statistical tests used in survival analysis.
These include the logrank test and all of its variants.
"""

import numpy as np
import pandas as pd

from scipy.stats import chi2, norm
from sksurv.nonparametric import kaplan_meier_estimator
from typing import Union, List, Tuple

# type aliases
Table = Union[List[List[float]],np.ndarray,pd.core.frame.DataFrame]
Vector = Union[List[float],np.ndarray,pd.core.frame.DataFrame]


def logrank_test(dfs: List[pd.core.frame.DataFrame], evnt_cols:List[str], time_cols:List[str],
                 type: str='CMH_logrank') -> Tuple[float, int, float]:
    """
    Cochran-Mantel-Haenszel
        # Suppose we have k (2x2) tables, all independent and we want to test for
        # a common group effect. This is the Cochran-Mantel-Haenszel that, in case
        # of survival analysis, tests for odds ratio being 1 at all times
        #
        # Null hypothesis H0: odds ratio is 1 for all tables

    Linear Rank
        # The linear rank version of the logrank test is based on adding up "scores" for one of the 2 groups
        # The "score" is based on the Nelson-Aalen estimator
        #
        # Null hypothesis H0: odds ratio is 1 for all tables

    Paramters
    --------
    dfs: List
        List of pandas dataframes with evnt_col and time_col column names.
    evnt_cols: str
        Name of the columns with events.
    time_cols: str
        Name of the columns with times.
    type: str
        The name of the statistical test desired.

    Returns
    -------
    statistic, df, p_value: float, int, float
        The value of the statistic, the number of degrees of freedom and the associated p-value.
    """
    list_data_y = []
    for df, evnt_col, time_col in zip(dfs, evnt_cols, time_cols):
        df_slct = df[[evnt_col, time_col]].dropna()
        data_y = np.array([(event, time) for event, time in zip(df_slct[evnt_col], df_slct[time_col])],
                          dtype=[('event', '?'), ('time', float)])
        list_data_y.append(data_y)
    n_groups = len(list_data_y)
    n = sum([l.shape[0] for l in list_data_y])
    df = n_groups - 1

    # Row vector d sums number of events (deaths) for first n_groups groups per each distinct time
    d = np.zeros((1, n_groups - 1))
    # Row vector d sums expected number of events (deaths) for first n_groups-1 groups over each distinct timei
    E = np.zeros((1, n_groups - 1))
    # Matrix V will sum variances of number of events (deaths) for first n_groups-1 groups over distinct times
    V = np.zeros((n_groups - 1, n_groups - 1))
    # List of list W will record scores for each group at each distinct time
    W = []
    # r_t stores at each time t the number at risk in each group
    r_t = [l.shape[0] for l in list_data_y]

    # Lbda is the Nelson-Aalen estimate of the cumulative hazard for all groups combined
    # It is used as a score for the linear rank test statistic
    lbda = 0

    for t in sorted(np.unique(np.concatenate([l['time'] for l in list_data_y]).flatten())):
        d_t = []
        E_t = []
        W_t = []

        for k in range(n_groups):
            d_t.append(list_data_y[k][(list_data_y[k]['time'] == t) & (list_data_y[k]['event'] == True)].shape[0])
        for k in range(n_groups):
            E_t.append(sum(d_t) * r_t[k] / sum(r_t))

        d += np.array([d_t[:-1]])
        E += np.array([E_t[:-1]])
        lbda += sum(d_t) / sum(r_t)

        for k in range(n_groups):
            w_k = []
            for s in list_data_y[k][list_data_y[k]['time'] == t]:
                if s['event'] == True:
                    w_k.append(1 - lbda)
                else:
                    w_k.append(-lbda)
            W_t.append(w_k)
        W.append(W_t)

        if sum(r_t) > 1:
            V_t = []
            for i in range(n_groups):
                r_i = []
                for j in range(n_groups):
                    if i != j:
                          # Pay attention to the - sign for covariance term
                        r_i.append(-(sum(r_t) - sum(d_t)) * sum(d_t) * r_t[i] * r_t[j] / ((sum(r_t) - 1) * sum(r_t) ** 2))
                    else:
                        r_i.append((sum(r_t) - sum(d_t)) * sum(d_t) * r_t[i] * (sum(r_t) - r_t[i]) / ((sum(r_t) - 1) * sum(r_t) ** 2))
                V_t.append(r_i)
            V_t = np.array(V_t)
            V += (V_t[:-1, :-1])
        r_t = [r - e[e['time'] == t].shape[0] for r, e in zip(r_t, list_data_y)]  # Update number of patients at risk in each group

    if type == 'CMH_logrank':
        statistic = np.float((d - E).dot(np.linalg.inv(V).dot((d - E).T)))
        p_value = 1 - chi2.cdf(statistic, df)

    elif type == 'Linear_rank_logrank':
        if n_groups > 2:
            raise ValueError("Linear_rank_logrank statistic is not yet implemented for comparison between more than 2 groups.")
        else:
            S = np.concatenate(pd.DataFrame(W).iloc[:, 0]).sum()
            V = list_data_y[0].shape[0] * list_data_y[1].shape[0] * sum([sum(np.concatenate(pd.DataFrame(W).iloc[:, i]) ** 2) for i in range(n_groups)]) / (n * (n - 1))
            statistic = S ** 2 / V
            p_value = 1 - chi2.cdf(statistic, df)

    else:
        raise ValueError("Unsupported type %s. Choose one of: 'CHM_logrank', 'linear_rank_logrank'." % type)

    return statistic, df, p_value


def wilcoxon_test(X: Vector, Y: Vector, type: str='normal_2_sided') -> Tuple[float, float]:
    """
    The present implements the Mann-Whitney for of the Wilcoxon test with
    null hypothesis that samples X and Y are from the same distribution.

    Parameters
    ----------
    X: array-like
        Vector of X samples
    Y: array-like
        Vector of Y samples
    type: str
        The name of the statistical test desired.

    Returns
    -------
    statistic, p_value: float, float
        The value of the statistic and the associated p-value.
    """
    df_X = pd.DataFrame(X, columns=['value'])
    df_X.loc[:, 'group'] = 'X'
    df_Y = pd.DataFrame(Y, columns=['value'])
    df_Y.loc[:, 'group'] = 'Y'

    df_XY = pd.concat((df_X, df_Y))
    df_XY.sort_values(by='value', inplace=True)
    df_XY.reset_index(inplace=True, drop=True)

    df_XY.loc[:, 'rank'] = df_XY.index + 1

    series_sum_ranks = df_XY.groupby('value').sum()
    series_counts = df_XY.groupby('value').count()['rank']

    df_ranks = pd.DataFrame(series_sum_ranks.divide(series_counts, axis=0)).reset_index()
    del df_XY['rank']

    df_XY = df_XY.merge(df_ranks, on='value')

    m = len(X)
    n = len(Y)

    S = df_XY.loc[df_XY.group == 'X', 'rank'].sum()
    E = m * (m + n + 1) / 2.
    V = m * n * (m + n + 1) / 12.

    if type == 'normal_2_sided':
        statistic = (S - E) / np.sqrt(V)
        p_value = 2 * norm.cdf(-np.abs(statistic), loc=0, scale=1.0)

    elif type == 'normal_1_sided':
        statistic = (S - E) / np.sqrt(V)
        p_value = norm.cdf(-np.abs(statistic), loc=0, scale=1.0)

    return statistic, p_value


def surv_wilcoxon_test(dfs: List[pd.core.frame.DataFrame], type: str='Gehan') -> Tuple[float, float]:
    """
    The present implements generalized Wilcoxon tests extended to censored data.
    Gehan's and Peto-Prentice's versions are implemented for now.

    Paramters
    --------
    dfs: List
        List of pandas dataframes with 'os_status' and 'os' column names.
    type: str
        The name of the statistical test desired.

    Returns
    -------
    statistic, p_value: float,, float
        The value of the statistic and the associated p-value.
    """
    list_df_groups = []

    for i, df in enumerate(dfs):
        df_slct = df[['os_status', 'os']].dropna()
        data_y = np.array([(status, time) for status, time in zip(df_slct.os_status, df_slct.os)],
                          dtype=[('os_status', '?'), ('os', np.float)])
        df_group = pd.DataFrame(data_y)
        df_group.loc[:, 'group'] = i
        list_df_groups.append(df_group)

    df_groups = pd.concat(list_df_groups)
    df_groups.sort_values(by='os', inplace=True)
    df_groups.reset_index(drop=True, inplace=True)

    if type == 'Gehan':
        df_groups.loc[:, 'U'] = np.nan
        for i in range(df_groups.shape[0]):
            U_i = 0
            t_i = df_groups.iloc[i].os
            status_i = df_groups.iloc[i].os_status

            for j in range(df_groups.shape[0]):
                t_j = df_groups.iloc[j].os
                status_j = df_groups.iloc[j].os_status
                if (t_i > t_j and status_i == True and status_j == True) or (t_i >= t_j and status_i == False and status_j == True):
                    U_i += 1
                elif (t_i < t_j and status_i == True and status_j == True) or (t_i <= t_j and status_j == False and status_i == True):
                    U_i -= 1
                else:
                    U_i += 0

            df_groups.iloc[i, -1] = U_i

        m = list_df_groups[0].shape[0]
        n = list_df_groups[1].shape[0]

        var = m * n * df_groups.U.apply(lambda x: x ** 2).sum() / ((m + n) * (m + n - 1))
        statistic = (df_groups.loc[df_groups.group == 0].U.sum()) ** 2 / var
        p_value = 1 - chi2.cdf(statistic, 1)
    elif type == 'Peto-Peto-Prentice':
        df_groups.loc[:, 'W'] = np.nan
        time, survival_prob = kaplan_meier_estimator(df_groups.os_status.values, df_groups.os.values)
        df_km = pd.DataFrame({'os': time, 'prob': survival_prob})
        df_groups = df_groups.merge(df_km, on='os', how='left')
        survival_prob_dict = {k: v for k, v in zip(time, survival_prob)}

        t_ = df_groups.os.min()
        df_groups.loc[(df_groups.os == t_) & (df_groups.os_status == True), 'W'] = survival_prob_dict[t_]
        df_groups.loc[(df_groups.os == t_) & (df_groups.os_status == False), 'W'] = survival_prob_dict[t_] - 1

        for t in df_groups.os.unique()[1:]:
            df_groups.loc[(df_groups.os == t) & (df_groups.os_status == True), 'W'] = survival_prob_dict[t] + survival_prob_dict[t_] - 1
            df_groups.loc[(df_groups.os == t) & (df_groups.os_status == False), 'W'] = survival_prob_dict[t] - 1
            t_ = t

        m = list_df_groups[0].shape[0]
        n = list_df_groups[1].shape[0]

        var = m * n * df_groups.W.apply(lambda x: x ** 2).sum() / ((m + n) * (m + n - 1))
        statistic = (df_groups.loc[df_groups.group == 0].W.sum()) ** 2 / var
        p_value = 1 - chi2.cdf(statistic, 1)

    return statistic, p_value
