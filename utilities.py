import time

import multiprocess as mp
import pandas as pd

from thefuzz import fuzz

from typing import Any, List, Optional, Iterable

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
        return [i for (i, desc) in enumerate(descriptions) if fuzz.partial_ratio(term, desc.lower()) > fuzzy_threshold]

    index_sets = [set(fuzzy_index_search(term_i, descriptions, fuzzy_threshold=fuzzy_threshold)) for term_i in term]

    final_indices = index_sets[0]
    for indices in index_sets[1:]:
        if and_search:
            final_indices = final_indices.intersection(indices)
        else:
            final_indices = final_indices.union(indices)

    return sorted(final_indices)


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


def add_plt_labels(ax, x: str, y: str, title: Optional[str] = None, **kwargs) -> None:
    """ Adds plot labels"""
    ax.set_xlabel(titleize(x))
    ax.set_ylabel(titleize(y))

    if title:
        ax.set_title(titleize(title))


########################################################################################################################
### End ###
########################################################################################################################
