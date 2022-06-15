import pandas as pd

from typing import Any, Optional

####################################################################################################
### Constants ###
####################################################################################################


def display_df(df: pd.DataFrame, max_rows: Optional[int] = None, max_columns: Optional[int] = None) -> None:
    """ Simple function to easily display dataframes"""
    
    with pd.option_context('display.max_rows', max_rows, 'display.max_columns', max_columns):
        display(df)


####################################################################################################
### Plotting Utilities ###
####################################################################################################


def titleize(label: str) -> str:
    """ Makes label into title format"""
    return label.replace("_", " ").title()


def add_plt_labels(ax, x: str, y: str, title: Optional[str] = None, **kwargs) -> None:
    """ Adds plot labels"""
    ax.set_xlabel(titleize(x))
    ax.set_ylabel(titleize(y))
    
    if title:
        ax.set_title(titleize(title))


####################################################################################################
### End ###
####################################################################################################
