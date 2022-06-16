import pandas as pd

from typing import Any, Optional

########################################################################################################################
### Constants ###
########################################################################################################################


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


def notebook_tqdm_import():
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
