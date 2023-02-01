"""
The :mod:`pyprism.survival` module includes standard functions for survival analyses.
"""

from .kaplan_meier import PlotKMConfig, plot_km, PlotCompareKMConfig, plot_compare_km

__all__ = [
    'PlotKMConfig',
    'plot_km',
    'PlotCompareKMConfig',
    'plot_compare_km'
]
