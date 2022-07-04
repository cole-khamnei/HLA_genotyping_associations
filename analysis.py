import os

import numpy as np
import pandas as pd
import scipy.stats as stats

import matplotlib.pyplot as plt

from typing import Any, List, Optional, Tuple, Union

import constants, utilities

if utilities.is_jupyter_notebook():
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm


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

    disease_exposure = np.sum(disease & exposure)
    disease_no_exposure = np.sum(disease & no_exposure)
    no_disease_exposure = np.sum(no_disease & exposure)
    no_disease_no_exposure = np.sum(no_disease & no_exposure)

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
                    data: Optional[pd.DataFrame] = None, variable_baseline:
                    Any = None, CI: Optional[float] = 95, as_text: bool = False) -> str:
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


def variable_OR_plot(data: pd.DataFrame, illness: str, illness_feature: str, ax = None,
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

        ax.text(.01, .70, text, transform=ax.transAxes, va="top", ha="left", fontdict={'family' : 'monospace'})

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
### HLA tools ###
########################################################################################################################


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

    if not isinstance(illness, str):
        fuzzy = True
        illness_tokens = illness
    
    illness_value = np.zeros(len(data))
    for feature in features:
        if not fuzzy:
            illness_value = (data[feature] == illness) | illness_value
        else:
            for illness_token in illness_tokens:
                illness_value = (data[feature].apply(lambda s: illness_token in s if not pd.isnull(s) else False)) | illness_value

    return illness_value


########################################################################################################################
### End ###
########################################################################################################################
