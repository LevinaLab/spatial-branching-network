'''
This script contain the basic plotting functions and plotting parameters (e.g., color, font, etc.) for the figures.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import farrow_and_ball as fab #for fancy colors

from mpmath import *

#------------------- fonts
plt.rcParams["axes.edgecolor"] = "k"
plt.rcParams["axes.facecolor"] = "w"
plt.rcParams["axes.linewidth"] = "1"
plt.rcParams.update({'font.size': 10})
plt.rcParams["font.family"] = "Helvetica"

plt.rcParams['pdf.fonttype'] = 42 # prepare as vector graphic
plt.rcParams['ps.fonttype'] = 42

alphFont = 12 # font size for figure alphabets
titleFont = 10 # font size of figure titles

#------------------- colors

#cK[k][L] gives the color that corresponds to k, L case

reds = fab.get_palette(fab.BaseColorPalette.REDS)
greens = fab.get_palette(fab.BaseColorPalette.GREENS)
blues = fab.get_palette(fab.BaseColorPalette.BLUES)


"""
cK = {1: {8: '#806934', 32: '#BF9D4E', 128: '#FFD269'},
      2: {8: '#462352', 32: '#A251BD', 128: '#D86DFC'},
      3: {8: '#265C36', 32: '#409C5B', 128: '#5ADB81'}}
"""
cK = {1: {8: reds[4], 16:reds[4], 32: reds[2], 128: reds[0], 64:reds[2], 512: reds[0], -1:reds[3]},
      2: {8: greens[3], 16:greens[3], 32: greens[1], 128: greens[0]},
      3: {8: blues[4], 16:blues[4], 32: blues[2], -1:blues[3], 128: blues[1]}}

cBP = 'gray' # mean-field BP color
cExponent = 'gray' # color of exponent line for avalanches

#cR = {8: '#344B80', 32: '#4E70BF', 128: '#6997FF'}
#cR = {8: '#313F66', 32: '#566EB3', 128: '#7A9EFF'}

def get_color_shades(color, ncolors=3, sat_inc=0.1, lum_inc=0.0):
    """
    Creates a set of color shades based on a given color
    color: the color, in any format that matplotlib can recognize
    ncolors: number of colors to generate
    sat_inc: change of saturation in each color. Can be negative
    lum_inc: change of luminance in each color. Can be negative.
    """
    rgb = mcolors.to_rgb(color)
    hsv = mcolors.rgb_to_hsv(rgb)
    
    new_colors = np.empty((ncolors,3))
    
    for j in range(ncolors):
        saturation = np.clip(hsv[1] + sat_inc*j, 0.0, 1.0)
        luminance =  np.clip(hsv[2] + sat_inc*j, 0.0, 1.0)
        new_colors[j] = mcolors.hsv_to_rgb([hsv[0], saturation, luminance])
        
    
    return new_colors
    
cR = cK[1]
# colors for rewiring
# cRWL1 = '#344B80'
# cRWL2 = '#4E70BF'
# cRWL3 = '#6997FF'

#------------------- legend properties
hL = 1.7 # legend handle length
hPad = 0.5 # legend handle textpad
ls = '-'

#------------------- axis, figure size properties (in cm)
fig_width_2col = 21
fig_width_1col = 11

avsize_lable = r'$S$'
prob_lable = r'$P(S)$'


#------------------ others
lw = 2 # line width

#------------------ utils functions

def to_inches(cm):
    """
    Convert cm to inches
    """
    return cm/2.54

def one_col_fig(height=fig_width_1col/1.618):
    """
    Returns a tuple (w,h), where w is the width if a single-column graph.
    Height by default is the golden ratio of the width, but can be chosen (in cm)
    """
    width = to_inches(fig_width_1col)
    height = to_inches(height)
    return (width, height)

def two_col_fig(height=fig_width_2col/1.618):
    """
    Returns a tuple (w,h), where w is the width if a double-column graph.
    Height by default is the golden ratio of the width, but can be chosen (in cm)
    """
    width = to_inches(fig_width_2col)
    height = to_inches(height)
    return (width, height)

def despine(axs):
    """
    Despine both a single axe or an array of axes
    """
    
    if type(axs) is np.ndarray:
        for ax in np.ravel(axs):
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
    else:
        axs.spines['right'].set_visible(False)
        axs.spines['top'].set_visible(False)

def label_axes(axs, textpos, uppercase=False, bracket=True):
    """
    Fast way to label all the diagrams
    """
    
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    if uppercase:
        alphabet = alphabet.upper()

    if bracket:
        axlabel = "({0})"
    else:
        axlabel = "{0}"
    
    for i,ax in enumerate(axs.flatten()):
        ax.text(textpos[0], textpos[1], axlabel.format(alphabet[i]), color='k', transform=ax.transAxes, weight="bold")


#------------------ plotting functions
def MLE_powerlaw_1Dgridsearch(data, alpha_range = np.arange(0.01,2.2,0.01), xmin = 10, xmax = 10**2,\
                              plot_it = 0, cFit= None, label_fit= None):
    """Estimate power-law exponent with MLE.

    Parameters
    -----------
    data : 1d array
        avalanche size or time distribution.
    alpha_range : 1d array
        range for the grid search of exponent.
    xmin: int
        minimum avalanche size for fitting the power law.
    xmax : int
        maximum avalanche size for fitting the power law.
    plot_it, cFit, label_fit : 
        parameters for plotting the log-likelihood versis alpha
       
    Returns
    -------
    best_alpha : float
        Estimated power-law exponent.  
    MLLE: float
        Maximum log-likelihood
    """
    
    
    def L(alpha,xmin,xmax,data):
        n = len(data)
        return -n * ln(zeta(alpha,xmin,0,method = 'borwein')- zeta(alpha,xmax,0,method = 'borwein'))- alpha*np.sum(np.log(data))

    data = data[data <= xmax]
    data = data[(data>=xmin)]
    log_like = np.array([L(k,xmin,xmax,data) for k in alpha_range])
    best_alpha = alpha_range[np.where(log_like == max(log_like))[0]]
    
    if plot_it:
        plt.plot(alpha_range,log_like,linewidth =2, color = cFit, label = label_fit)
        plt.xlabel(r'$\alpha$')
        plt.ylabel('log-likelihood')
    return best_alpha[0], max(log_like)

def xmaxRatio(av, ratio):
    """ 
    compute maximum avalanche size based on the given ratio of distribution.
    av: avalanche sizes
    ratio: defines percentile of distribution (number between [0,1])
    """
    data_list = np.ndarray.tolist(av)
    uniq_list = np.sort(data_list)

    N = int(len(uniq_list)*(ratio))-1
    xmax =int(uniq_list[N])
    return xmax


def plot_av(ax, av, minlogbin = 9, maxlogbin = 10**8, nlogbings = 70, col = 'k', lw = 2, ls = '-', label ='av',\
            plot_fit = 0, ls_powerlaw = '-', \
            alpha_range = np.arange(0.01,2.2,0.01), xmin = 10, xmax_percentile = 0.96, fit_zeroValue = 1, fitlinepos = 1.0):
    """Plot avalanche size/time duration, with plotting the fit on top of the fitted region.

    Parameters
    -----------
    ax : object
        plotting axis.
    av : 1d array
        avalanche size or time distribution.
    minlogbin : int
        minimum avalanche size to start log binning.
    maxlogbin : int
        maximum avalanche size for log binning.
    nlogbings : int
        number of log bins for histogram. 
    col, lw, ls, label : 
        plotting parameters
    plot_fit: boealian
        If plot the power-law fit line.
    fit_zeroValue: float
        The value of fitted power-law line at avalanche-size = 1.
    ls_powerlaw: string
        Line style for power-law line
    alpha_range, xmin, xmax:
        Parameters for fitting avalanches with power law. You can set the xmax using the percentiles of distribution.
       
    Returns
    -------
    alpha : float
        Estimated power-law exponent.     
    """
    # log binning for large avalanches
    logmin = np.log10(minlogbin)
    logmax = np.log10(maxlogbin)
    bins_log = np.logspace(logmin, logmax, nlogbings)
    bins_log = bins_log[bins_log>minlogbin]

    # linear binning for small avalanches
    minx, maxx, nx = 1,minlogbin,1
    bins_lin = np.arange(minx, maxx, nx)

    bins = np.concatenate((bins_lin,bins_log))
    dist = np.histogram(av, bins=bins, density=True)[0]
    ax.loglog(bins[:-1][dist>0], dist[dist>0], color=col, lw=lw, ls = ls, label=label)
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    if plot_fit:
        xmax = xmaxRatio(av, xmax_percentile)
        alpha, MLLE = MLE_powerlaw_1Dgridsearch(av, alpha_range, xmin, xmax)
        normalization = fit_zeroValue
        bins = np.arange(xmin, xmax+1)
        ax.plot(bins, fitlinepos * normalization*bins**(-alpha), color=cExponent, lw=1.5, ls = ls_powerlaw)
        return alpha
    else: return None
    
    

def plot_av_truncatedPowerLaw(ax, av, minlogbin = 9, maxlogbin = 10**8, nlogbings = 70, col = 'k', lw = 2, ls = '-', label ='av',\
            plot_fit = 0, ls_powerlaw = '-', \
            alpha_range = np.arange(0.01,2.2,0.01), xmin = 10, xmax_percentile = 0.95):
    """Plot avalanche size/time duration, with plotting the fit as the correct truncated power-law on top of the fitted region.

    Parameters
    -----------
    ax : object
        plotting axis.
    av : 1d array
        avalanche size or time distribution.
    minlogbin : int
        minimum avalanche size to start log binning.
    maxlogbin : int
        maximum avalanche size for log binning.
    nlogbings : int
        number of log bins for histogram. 
    col, lw, ls, label : 
        plotting parameters
    plot_fit: boealian
        If plot the power-law fit line.
    ls_powerlaw: string
        Line style for power-law line
    alpha_range, xmin, xmax:
        Parameters for fitting avalanches with power law. You can set the xmax using the percentiles of distribution.
       
    Returns
    -------
    alpha : float
        Estimated power-law exponent.     
    """
    # log binning for large avalanches
    logmin = np.log10(minlogbin)
    logmax = np.log10(maxlogbin)
    bins_log = np.logspace(logmin, logmax, nlogbings)
    bins_log = bins_log[bins_log>minlogbin]

    # linear binning for small avalanches
    minx, maxx, nx = 1,minlogbin,1
    bins_lin = np.arange(minx, maxx, nx)

    bins = np.concatenate((bins_lin,bins_log))
    dist = np.histogram(av, bins=bins, density=True)[0]
    ax.loglog(bins[:-1][dist>0], dist[dist>0], color=col, lw=lw, ls = ls, label=label)
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    if plot_fit:
        xmax = xmaxRatio(av, xmax_percentile)
        alpha, MLLE = MLE_powerlaw_1Dgridsearch(av, alpha_range, xmin, xmax)
        normalization = 1/(zeta(alpha,xmin,0,method = 'borwein')- zeta(alpha,xmax,0,method = 'borwein'))
        bins = np.arange(xmin, xmax+1)
        ax.plot(bins, normalization*bins**(-alpha), color=cExponent, lw=1, ls = ls_powerlaw)
        return alpha
    else: return None


def plot_av_general(ax, av, minlogbin = 9, maxlogbin = 10**8, nlogbings = 70, col = 'k', lw = 2, ls = '-', label ='av',\
            plot_fit = 0, fit_zeroValue = 1, ls_powerlaw = '-', \
            alpha_range = np.arange(0.01,2.2,0.01), xmin = 10, xmax = 10**2):
    """Plot avalanche size/time duration.

    Parameters
    -----------
    ax : object
        plotting axis.
    av : 1d array
        avalanche size or time distribution.
    minlogbin : int
        minimum avalanche size to start log binning.
    maxlogbin : int
        maximum avalanche size for log binning.
    nlogbings : int
        number of log bins for histogram. 
    col, lw, ls, label : 
        plotting parameters
    plot_fit: boealian
        If plot the power-law fit line.
    fit_zeroValue: float
        The value of fitted power-law line at avalanche-szie = 1.
    ls_powerlaw: string
        Line style for power-law line
    alpha_range, xmin, xmax:
        Parameters for fitting avalanches with power law. You can set the xmax using the percentiles of distribution.
       
    Returns
    -------
    alpha : float
        Estimated power-law exponent.     
    """
    # log binning for large avalanches
    logmin = np.log10(minlogbin)
    logmax = np.log10(maxlogbin)
    bins_log = np.logspace(logmin, logmax, nlogbings)
    bins_log = bins_log[bins_log>minlogbin]

    # linear binning for small avalanches
    minx, maxx, nx = 1,minlogbin,1
    bins_lin = np.arange(minx, maxx, nx)

    bins = np.concatenate((bins_lin,bins_log))
    dist = np.histogram(av, bins=bins, density=True)[0]
    ax.loglog(bins[:-1][dist>0], dist[dist>0], color=col, lw=lw, ls = ls, label=label)
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    if plot_fit:
        alpha, MLLE = MLE_powerlaw_1Dgridsearch(av, alpha_range, xmin, xmax)
        ax.plot(bins[:-1], fit_zeroValue*bins[:-1]**(-alpha), color=cExponent, lw=1, ls = ls_powerlaw)
        return alpha[0]
    else: return None




def plot_av_gamma(ax, av, minlogbin = 9, maxlogbin = 10**8, nlogbings = 70, col = 'k', lw = 2, ls = '-', label ='av',\
            plot_fit = 0, ls_powerlaw = '-', \
            alpha_range = np.arange(-2.5, -0.01, 0.01), xmin = 10, xmax = 1e3, fit_zeroValue = 1, fitlinepos = 1.0):
    """Plot avalanche size/time duration, with plotting the fit on top of the fitted region.

    Parameters
    -----------
    ax : object
        plotting axis.
    av : 1d array
        avalanche size or time distribution.
    minlogbin : int
        minimum avalanche size to start log binning.
    maxlogbin : int
        maximum avalanche size for log binning.
    nlogbings : int
        number of log bins for histogram. 
    col, lw, ls, label : 
        plotting parameters
    plot_fit: boealian
        If plot the power-law fit line.
    fit_zeroValue: float
        The value of fitted power-law line at avalanche-size = 1.
    ls_powerlaw: string
        Line style for power-law line
    alpha_range, xmin, xmax:
        Parameters for fitting avalanches with power law. You can set the xmax using the percentiles of distribution.
       
    Returns
    -------
    alpha : float
        Estimated power-law exponent.     
    """
    
    xaxis = np.arange(av.size)
    yaxis = av.copy()
    yaxis[yaxis <= 1e-18] = None
    ax.loglog(xaxis, yaxis, color=col, lw=lw, ls = ls, label=label)
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    if plot_fit:
        bins = np.arange(xmin, xmax+1)
        ax.plot(bins, fitlinepos * bins**plot_fit, color=cExponent, lw=1.5, ls = ls_powerlaw)
    return 