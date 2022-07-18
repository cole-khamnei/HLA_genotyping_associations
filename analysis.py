import os
import sys

import numpy as np
import pandas as pd
import scipy
import scipy.stats as stats
import KDEpy

import matplotlib.pyplot as plt

from typing import Any, List, Optional, Tuple, Union

import constants, utilities

if utilities.is_jupyter_notebook():
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

sys.path.append(constants.GRANTHAM_DISTANCE_PATH)

import grantham_distance as gd

########################################################################################################################
### Constants ###
########################################################################################################################


########################################################################################################################
### Simple Statistics ###
########################################################################################################################


def calculate_OR(disease: np.ndarray, exposure: np.ndarray,
                 disease_value: Any = True, exposure_value: Any = True,
                 CI: float = 95) -> Tuple[float, float, float, float, int]:
    """ Calculates the odds, p value, lower and upper bounds."""

    exposure, no_exposure = exposure == exposure_value, exposure != exposure_value
    disease, no_disease = disease == disease_value, disease != disease_value

    disease_exposure = np.count_nonzero(disease & exposure)
    disease_no_exposure = np.count_nonzero(disease & no_exposure)
    no_disease_exposure = np.count_nonzero(no_disease & exposure)
    no_disease_no_exposure = np.count_nonzero(no_disease & no_exposure)

    contingency_table = [[disease_exposure, no_disease_exposure], [disease_no_exposure, no_disease_no_exposure]]

    odds_ratio, p_value = stats.fisher_exact(contingency_table)

    if CI:
        standard_error = np.sqrt(np.sum(1 / np.array(contingency_table).ravel()))
        assert 100 > CI > 1, f"Given CI percent ({CI}%) is in invalid. Must be in (1, 100)"
        z_score = stats.norm.ppf(CI / 100)

        upper_bound = np.exp(np.log(odds_ratio) + z_score * standard_error)
        lower_bound = np.exp(np.log(odds_ratio) - z_score * standard_error)
    else:
        upper_bound, lower_bound = None, None

    return odds_ratio, p_value, disease_exposure, lower_bound, upper_bound


def variable_OR_test(illness: Union[str, np.ndarray], variable: Union[str, np.ndarray],
                     data: Optional[pd.DataFrame] = None, variable_baseline: Any = None,
                     CI: Optional[float] = 95, as_text: bool = False) -> str:
    """ """

    if isinstance(variable, str):
        assert data is not None
        variable_name = variable
        variable = data[variable]
    else:
        variable_name = ""

    if isinstance(illness, str):
        assert data is not None
        illness = data[illness]

    nan_indices = pd.isnull(illness) | pd.isnull(variable)
    illness, variable = illness[~nan_indices], variable[~nan_indices]

    text = ""
    variable_values, counts = np.unique(variable[~pd.isnull(variable)], return_counts=True)

    if variable_baseline is None:
        variable_baseline = variable_values[np.argmax(counts)]

    odds_ratio, p_value, N, lower_bound, upper_bound = calculate_OR(illness, variable != variable_baseline, CI=CI)

    if len(variable_values) > 2:
        p_value_str = f"{p_value:.2}" if "e" in str(p_value) else f"{p_value:.3f}"
        line = f"{variable_name} != {variable_baseline}: OR: {odds_ratio:.3f} p-value: {p_value_str}"
        if CI:
            line += f" 95% CI: {lower_bound:.2f} - {upper_bound:.2f}"
        line += f" N: {N:,}"
        text += line + "\n" if p_value > .05 else line + " ***\n"

    for variable_value in sorted(variable_values):
        if variable_value == variable_baseline:
            continue

        odds_ratio, p_value, N, lower_bound, upper_bound = calculate_OR(illness, variable == variable_value, CI=CI)

        p_value_str = f"{p_value:.2}" if "e" in str(p_value) else f"{p_value:.3f}"
        line = f"{variable_name} == {variable_value}: OR: {odds_ratio:.3f} p-value: {p_value_str}"
        if CI:
            line += f" 95% CI: {lower_bound:.2f} - {upper_bound:.2f}"
        line += f" N: {N:,}"
        text += line + "\n" if p_value > .05 else line + " ***\n"

    return text


def variable_OR_plot(data: pd.DataFrame, illness: str, illness_feature: str, ax=None,
                     x: str = "grantham_divergence", OR_variable: str = "zygosity", bw: float = .2,
                     no_illness_label: Optional[str] = None, title: str = ""):
    """"""

    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 6))
    else:
        fig = None

    if isinstance(illness_feature, str):
        illness_values = data[illness_feature] == illness
    else:
        illness_values = illness_feature

    if not no_illness_label:
        no_illness_label = "no " + illness

    illness_labels = np.array([no_illness_label, illness])[1 * illness_values]
    utilities.kde_plot(data=data, x=x, hue=illness_labels, ax=ax, bw=bw)

    if OR_variable:
        OR_variables = [OR_variable] if isinstance(OR_variable, str) else OR_variable
        text = ""
        for OR_variable in OR_variables:
            text += variable_OR_test(illness=illness_values, variable=OR_variable, data=data, CI=None)

        ax.text(.01, .70, text, transform=ax.transAxes, va="top", ha="left", fontdict={'family': 'monospace'})

    disease_distances = data[x].loc[illness_values]
    no_disease_distances = data[x].loc[~illness_values]
    ks_stat, ks_p_value = stats.kstest(disease_distances, no_disease_distances)

    title = f"{utilities.titleize(illness)}{title}"
    title += f"\nKS Test: Stat={ks_stat:.03}, p={ks_p_value:.03}"
    if ks_p_value < 0.05:
        title += " **"

    ax.set_title(title)

    ax.legend(title="Status", loc="upper left")
    utilities.add_plt_labels(ax, x=x, y="Density")

    return fig, ax


########################################################################################################################
### HLA grantham distance analysis ###
########################################################################################################################

grantham_distance = gd.GranthamDistance(gd.GRANTHAM_DISTANCE_MATRIX_PATH)
memoized_grantham_distance = utilities.MemoizedFunction(grantham_distance.sequence_pair_distance)


def grantham_linkage(p1, p2):
    """"""
    inter_person_distances = np.array([[memoized_grantham_distance._call(seq1, seq2) for seq2 in p2] for seq1 in p1])
    row_ind, col_ind = scipy.optimize.linear_sum_assignment(inter_person_distances)
    return inter_person_distances[row_ind, col_ind].sum()


memoized_grantham_linkage = utilities.MemoizedFunction(grantham_linkage)


########################################################################################################################
### Helpers ###
########################################################################################################################


def get_base_feature(feature: str) -> str:
    """ returns the features base name"""
    return "_".join(feature.split("_")[:-1]) if "." in feature else feature


def get_multiple_features_from_base_feature(biobank_data: pd.DataFrame, base_feature: str) -> List[str]:
    """ Finds all the features that match with the base feature"""
    return [feature for feature in biobank_data.columns if get_base_feature(feature) == base_feature]


def get_illness_value(data: pd.DataFrame, illness: str, base_feature: str, fuzzy: bool = False) -> np.ndarray:
    """ Finds all the individuals with a specific illness"""
    assert base_feature in data, f"Invalid base_feature: {base_feature}"

    features = get_multiple_features_from_base_feature(data, base_feature)

    if not isinstance(illness, str) or fuzzy:
        fuzzy = True
        illness_tokens = illness

        def fuzzy_test(s, illness_token):
            return illness_token in s if not pd.isnull(s) else False

    illness_value = np.zeros(len(data))
    for feature in features:
        if not fuzzy:
            illness_value = (data[feature] == illness) | illness_value
        else:
            for illness_token in illness_tokens:
                illness_value = data[feature].apply(fuzzy_test, args=(illness_token)) | illness_value

    return illness_value


def get_illness_value_dx_age(data: pd.DataFrame, illness: str, base_feature: str, fuzzy: bool = False
                             ) -> np.ndarray:
    """ Finds all the individuals with a specific illness"""

    assert base_feature in data, f"Invalid base_feature: {base_feature}"
    features = get_multiple_features_from_base_feature(data, base_feature)

    base_feature_dx_age = base_feature.strip("_code") + "_dx_age_interpolated"
    assert base_feature_dx_age in data, f"Dx age feature not found in data '{base_feature_dx_age}'"
    dx_age_features = get_multiple_features_from_base_feature(data, base_feature_dx_age)

    if not isinstance(illness, str) or fuzzy:
        fuzzy = True
        illness_tokens = illness

        def fuzzy_test(s, illness_token):
            return illness_token in s if not pd.isnull(s) else False

    illness_value = np.zeros(len(data))
    illness_dx_ages = -99 * np.ones(len(data))

    for feature, feature_dx_age in zip(features, dx_age_features):
        if not fuzzy:
            feature_values = data[feature] == illness
        else:
            feature_values = np.zeros(len(data))
            for illness_token in illness_tokens:
                feature_values = data[feature].apply(fuzzy_test, args=(illness_token)) | feature_values

        ages = data[feature_dx_age].values[feature_values]
        illness_dx_ages[feature_values] = ages
        illness_value = feature_values | illness_value

    return illness_value, illness_dx_ages


def create_sampled_distribution(sampling: np.ndarray, max_value: Optional[float] = None,
                                min_value: Optional[float] = None, bw: float = 0.05,
                                n_sample: int = 1000, force_kde: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """"""
    max_value = max_value if max_value is not None else np.max(sampling) + .05 * np.abs(np.max(sampling))
    min_value = min_value if min_value is not None else np.min(sampling) - .05 * np.abs(np.min(sampling))

    assert max_value >= np.max(sampling)
    assert min_value <= np.min(sampling)

    if not force_kde and all(value % 1 == 0 for value in sampling):
        max_value = int(np.ceil(max_value)) + 1
        sampling = sampling.astype(int)
        return np.arange(max_value), np.bincount(sampling, minlength=max_value) / len(sampling)

    evaluate_grid = np.linspace(min_value, max_value, n_sample)

    bw = bw * np.std(sampling)
    bw = 0.001 if bw == 0 else bw
    kde = KDEpy.FFTKDE(bw=bw, kernel="gaussian")
    y = kde.fit(sampling)(evaluate_grid)
    return evaluate_grid, y / np.sum(y)


def probability_y_greater_than_x(x_sampling: np.ndarray, y_sampling: np.ndarray) -> float:
    """"""
    max_value = np.ceil(max(np.max(x_sampling), np.max(y_sampling))) + 1
    min_value = np.floor(min(np.min(x_sampling), np.min(y_sampling))) - 1

    y_values, y_probs = create_sampled_distribution(y_sampling, min_value=min_value, max_value=max_value)
    x_values, x_probs = create_sampled_distribution(x_sampling, min_value=min_value, max_value=max_value)

    assert all(x_values == y_values)

    return np.sum(np.append([0], np.cumsum(x_probs))[:-1] * y_probs)  # honk for vector ops


def preceding_event_test(precursor_dx_ages, illness_dx_ages):
    """"""
    precursor_dx_ages = precursor_dx_ages.values if isinstance(precursor_dx_ages, pd.Series) else precursor_dx_ages
    illness_dx_ages = illness_dx_ages.values if isinstance(illness_dx_ages, pd.Series) else illness_dx_ages

    illness_values = illness_dx_ages >= 0
    precursor_values = precursor_dx_ages >= 0

    illness_and_infection = precursor_values & illness_values
    n_with_illness_and_infection = np.count_nonzero(illness_and_infection)
    n_with_illness_and_prior_infection = np.count_nonzero(illness_and_infection & (illness_dx_ages > precursor_dx_ages))

    illness_ages = illness_dx_ages[illness_values]
    precursor_ages = precursor_dx_ages[precursor_values]

    p = probability_y_greater_than_x(precursor_ages, illness_ages)

    if n_with_illness_and_infection <= 2:
        return p, None

    r = stats.binomtest(n_with_illness_and_prior_infection, n_with_illness_and_infection, p=p)
    return p, r.pvalue


########################################################################################################################
### End ###
########################################################################################################################
