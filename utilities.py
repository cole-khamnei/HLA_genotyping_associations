import time
import sys
from gc import get_referents
import six

from typing import Any, Iterable, List, Optional, Tuple, Union
from types import ModuleType, FunctionType

import matplotlib.pyplot as plt
import multiprocess as mp
import numpy as np
import pandas as pd
import seaborn as sns

import KDEpy

from thefuzz import fuzz

########################################################################################################################
### Multiprocessing ###
########################################################################################################################


########################################################################################################################
### Random ###
########################################################################################################################


class Timer:
    def __init__(self, print_on_exit: bool = True):
        self.start_time = time.time()
        self.final_time = False
        self.print_on_exit = print_on_exit

    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.get_elapsed_time()
        self.final_time = True
        if self.print_on_exit:
            self.print_time()

    def get_elapsed_time(self):
        """"""
        if not self.final_time:
            self.elapsed_time = time.time() - self.start_time

    def print_time(self):
        """"""
        self.get_elapsed_time()
        print(f"Elapsed time: {self.elapsed_time:,.04f} seconds")


def fuzzy_index_search(term: str, descriptions: Iterable[str], fuzzy_threshold: int = 95,
                       and_search: bool = False) -> List[int]:
    """ Searches a list of descriptions and returns any with a fuzzy index above a threshold."""

    if isinstance(term, str):
        term = term.lower()
        return [i for (i, desc) in enumerate(descriptions) if term in desc.lower()]
        # return [i for (i, desc) in enumerate(descriptions) if fuzz.partial_ratio(term, desc.lower()) > fuzzy_threshold]

    index_sets = [set(fuzzy_index_search(term_i, descriptions, fuzzy_threshold=fuzzy_threshold)) for term_i in term]

    final_indices = index_sets[0]
    for indices in index_sets[1:]:
        if and_search:
            final_indices = final_indices.intersection(indices)
        else:
            final_indices = final_indices.union(indices)

    return sorted(final_indices)


def exclude(main_list: List[Any], exclude_list: List[Any]) -> List[Any]:
    """ Excludes item found in another list"""
    return [item for item in main_list if item not in exclude_list]


########################################################################################################################
### Get size ###
########################################################################################################################


"""
https://stackoverflow.com/questions/449560/how-do-i-determine-the-size-of-an-object-in-python
"""

# Custom objects know their class.
# Function objects seem to know way too much, including modules.
# Exclude modules as well.
BLACKLIST = type, ModuleType, FunctionType


def get_size(obj):
    """sum size of object & members."""

    if isinstance(obj, BLACKLIST):
        raise TypeError('get_size() does not take argument of type: ' + str(type(obj)))

    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


########################################################################################################################
### Random ###
########################################################################################################################


def multiprocess_pool(function, param_sets: list, n_processes: int = 4):
    """ Runs a function in parallel using the pool multiprocess"""

    with mp.Pool(n_processes) as pool:
        results = [pool.apply_async(function, param_set) for param_set in param_sets]
        return [result.get() for result in results]


########################################################################################################################
### jupyternotebook utilities ###
########################################################################################################################


def is_jupyter_notebook():
    """
    https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    """
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter


def tqdm_import():
    global tqdm
    if is_jupyter_notebook():
        from tqdm.notebook import tqdm
    else:
        from tqdm import tqdm


def display_df(df: pd.DataFrame, max_rows: Optional[int] = None, max_columns: Optional[int] = None) -> None:
    """ Simple function to easily display dataframes"""

    with pd.option_context('display.max_rows', max_rows, 'display.max_columns', max_columns):
        display(df)


########################################################################################################################
### Plotting Utilities ###
########################################################################################################################


def titleize(label: str) -> str:
    """ Makes label into title format"""
    return label.replace("_", " ").title()


def add_plt_labels(ax, x: Optional[str] = None, y: Optional[str] = None, title: Optional[str] = None, **kwargs) -> None:
    """ Adds plot labels"""
    if x: 
        ax.set_xlabel(titleize(x))

    if y:
        ax.set_ylabel(titleize(y))

    if title:
        ax.set_title(titleize(title))


def create_subplot(n_plots: int, ncols: int = 2, width: float = 16, height_per: float = 3, hspace: float = 0.3, wspace: float = .2):
    """ this function is such a homie"""

    assert n_plots > 1, f"n_plots must be greater than 1, given: {n_plots}"

    n_rows = int(np.ceil(n_plots / ncols))
    offset = 2 * max(n_rows // 3, 1)

    fig = plt.figure(figsize=(width, offset + height_per * n_rows))
    gs = fig.add_gridspec(n_rows, ncols)

    axes = []
    i = 1
    for i in range(n_plots // ncols * ncols):
        axes.append(fig.add_subplot(gs[i // ncols, i % ncols]))

    if i != n_plots - 1:
        axes.append(fig.add_subplot(gs[n_rows - 1, :]))
    fig.tight_layout()
    plt.subplots_adjust(wspace=wspace, hspace=hspace)

    return fig, axes


def kde_smooth(values: np.ndarray, bw: Union[str, float] = "silverman", kernel: str = "gaussian",
               n_grid: int = 400, clip: Optional[Tuple[float, float]] = None, mean: bool = True):
    """"""

    if isinstance(bw, float):
        bw = bw * np.std(values)

    values = np.array(values)
    values = values[~pd.isnull(values)]

    assert len(values) >= 2, f"Values array must contained at least 2 non-null values."

    values_min, range_ = np.min(values), np.ptp(values)
    value_lower_lim, value_upper_lim = values_min - 0.05 * range_, values_min + 1.05 * range_
    
    kde = KDEpy.FFTKDE(bw=bw, kernel=kernel)
    x, y = kde.fit(values)(n_grid)
    
    if clip:
        clip_min, clip_max = clip
        start_index = 0
        multiplier = 1
        append_values = [values]
        if clip_min is not None:
            append_values.append(2 * clip_min - values)
            multiplier += 1
            start_index = n_grid
        else:
            clip_min = - np.inf
        
        if clip_max is not None:
            append_values.append(2 * clip_max - values)
            multiplier += 1
        else:
            clip_max = np.inf
        
        values = np.concatenate(append_values)

        x_fft, y_fft = KDEpy.FFTKDE(bw=bw, kernel=kernel).fit(values)(n_grid * multiplier)
        y_fft = y_fft * multiplier
        
        if mean:
            x_fft, y_fft = x, np.sqrt(y * y_fft[start_index: start_index + n_grid])
        else:
            x_fft, y_fft = x, y_fft[start_index: start_index + n_grid]


        clipped_indices = (x_fft <= clip_min) | (x_fft >= clip_max)
        return x_fft[~clipped_indices], y_fft[~clipped_indices]

    else:
        return x, y


def single_kde_plot(x: Union[np.ndarray, str], data: Optional[dict] = None,
                    shade: bool = True, bw: str = "silverman", kernel: str = "gaussian", mean: bool = True,
                    clip: Optional[Tuple[float, float]] = None, ax = None, **params):
    """"""
    if isinstance(x, str):
        assert data is not None, "Must provide 'data'"
        x = data[x]
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 6))

    if (~pd.isnull(x)).sum() < 2:
        ax.plot([], [])
        return ax

    label = params.pop("label", None)
    fill_label, line_label = (label, None) if shade else (None, label)

    x_, y_ = kde_smooth(x, bw=bw, mean=mean, clip=clip, kernel=kernel)
    
    baseline, = ax.plot(x_, y_, **params, linewidth=1, label=line_label)

    if shade:
        shade_alpha = params.get("alpha", .2)
        ax.fill_between(x_, y1=y_, alpha=shade_alpha, facecolor=baseline.get_color(), label=fill_label)
    
    return ax


def kde_plot(data, x, hue: str = None, ax=None, label: str = None, labels: list = None,
             shade: bool = True, alpha: float = .6, bw: float = .35, **params):
    """ plots a kdeplot"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    if hue:
        group_variable = hue
        group_values, counts = np.unique(data[group_variable], return_counts=True)
        sort_index = np.argsort(group_values)
        group_values, counts = group_values[sort_index], counts[sort_index]

        if not labels:
            labels = [f"{str(group_value).title()} ( N = {count:,})" for group_value, count in zip(group_values, counts)]

        i = 0
        for group_value, label in zip(group_values, labels):
            group_data = data.loc[data[group_variable] == group_value]
            single_kde_plot(data=group_data, x=x, ax=ax, shade=shade, alpha=alpha, bw=bw, label=label,
                            zorder=i, **params)
            i -= 1

    else:
        single_kde_plot(data=data, x=x, ax=ax, shade=shade, alpha=alpha, bw=bw, label=label)
    if label:
        ax.legend()

    return ax


def sns_wrapper(function, data, x, hue: str = None, ax=None, label: str = None, labels: list = None,
                alpha: float = .6, **params):
    """ plots a kdeplot"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    if hue:
        group_variable = hue
        group_values, counts = np.unique(data[group_variable], return_counts=True)
        sort_index = np.argsort(group_values)
        group_values, counts = group_values[sort_index], counts[sort_index]

        if not labels:
            labels = [f"{str(group_value).title()} ( N = {count:,})"
                      for group_value, count in zip(group_values, counts)]

        i = 0
        for group_value, label in zip(group_values, labels):
            group_data = data.loc[data[group_variable] == group_value]
            function(data=group_data, x=x, ax=ax, alpha=alpha, label=label, zorder=i, **params)
            i -= 1

    else:
        function(data=data, x=x, ax=ax, alpha=alpha, label=label, **params)

    if label:
        ax.legend(title=hue.title())

    return ax


def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, fig=None, **kwargs):
    """ stolen from stack overflow """
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in  six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])

    return fig, ax


########################################################################################################################
### End ###
########################################################################################################################
