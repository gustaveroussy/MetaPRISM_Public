"""
The :mod:`pyprism.util` module includes functions for performing technical operations like setting the current working
directory to a specific folder or saving results in images and tables.
"""

from ._setwd import setwd_to_data, setwd_to_results, setwd_to_scripts

from ._utils import explode_df, rm_duplicates_list

__all__ = [
    'setwd_to_data',
    'setwd_to_results',
    'setwd_to_scripts',
    'explode_df',
    'rm_duplicates_list'
]
