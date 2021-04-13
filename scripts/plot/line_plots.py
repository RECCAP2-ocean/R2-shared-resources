from matplotlib import pyplot as plt


def style_line_subplot(ax, add_zero_line=True, xlim=None, y_range=None):
    import numpy as np
    import pandas as pd
    
    if ax is None:
        ax = plt.gca()
    plt.sca(ax)
    
    plt.xticks(rotation=0, ha='center')
    plt.xlabel('')
    plt.title('')
    
    if ax.get_ylim()[0] < 0 < ax.get_ylim()[1]:
        ax.axhline(0, color='k', ls='--', lw=0.5, zorder=0)
    
    if xlim is None:
        xlim = ax.get_xlim()
    ax.set_xticks(np.arange('1980', '2020', 5, dtype='datetime64[Y]'))
    ax.set_xticklabels(np.arange(1980, 2020, 5))
    ax.set_xlim(*xlim)
    
    if y_range is not None:
        center = np.mean(ax.get_ylim())
        upp = center + y_range/2
        low = center - y_range/2
        ax.set_ylim(low, upp)
        
    return ax


def plot_ensemble_lines(
    da, 
    ax=None, 
    dim='variable', 
    draw_mean=True, 
    draw_members=True, 
    draw_stdev=True, 
    y_range=None, 
    color='C0', 
    name='', 
    **kwargs
):
    if ax is None:
        ax = plt.gca()
    
    x_dim = list(set(da.dims) - set([dim]))[0]
    x = da[x_dim].values
        
    ens = da
    avg = da.mean(dim)
    
    if draw_stdev:   
        std = da.std(dim)
        upr = avg + std
        lwr = avg - std
        # drawn first to maintain zorder
        ax.fill_between(x, upr, lwr, color=color, alpha=0.3, lw=0, **kwargs)
    
    if draw_members:
        ens.plot(c=color, lw=0.5, alpha=0.6, ax=ax, add_legend=False, hue=dim, **kwargs)
    
    if draw_mean:
        avg.plot(c=color, lw=2.5, alpha=1.0, ax=ax, label=name, **kwargs)
    
    xlim = x.min(), x.max()
    ax = style_line_subplot(ax, add_zero_line=False, xlim=xlim, y_range=y_range)
    
    return ax