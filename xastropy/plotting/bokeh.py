"""
#;+ 
#; NAME:
#; utils
#;    Version 1.0
#;
#; PURPOSE:
#;    Plotting utilities
#;   24-Nov-2014 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

try:
    import bokeh
except ImportError:
    print('-----------------------------------------------------------')
    print('-----------------------------------------------------------')
    print('WARNING: Not loading bokeh module in xastropy.plotting   \n Install bokeh if you want it')
    print('-----------------------------------------------------------')


from bokeh.io import output_notebook, show
from bokeh.plotting import figure
#from bokeh.models import Range1d

output_notebook()

def plot_spec_notebook(spec, title=None):
    """ Simple spectrum plot in a Notebook

    spec : XSpectrum1D
    title : str, optional

    """
    # Lya
    p = figure(plot_width=900, plot_height=500, title=title)
    # Data
    p.line(spec.dispersion.value, spec.flux.value, color='black', line_width=2)
    p.line(spec.dispersion.value, spec.sig, color='red', line_width=0.5)
    # Labels
    p.xaxis.axis_label = "Wavelength"
    p.yaxis.axis_label = "Flux"
    show(p)
