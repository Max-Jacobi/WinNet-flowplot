def plot_magn(ax, x=True, y=False, alpha=.5, color='k'):
    '''
    routine for ploting lines at magic numbers
    '''
    ixymagic = [2, 8, 20, 28, 50, 82, 126]
    for mn in ixymagic:
        ax.axhline(y=mn-.5, linestyle='--', color=color, linewidth=.5)
        ax.axhline(y=mn+.5, linestyle='--', color=color, linewidth=.5)
        ax.axvline(x=mn-.5, linestyle='--', color=color, linewidth=.5)
        ax.axvline(x=mn+.5, linestyle='--', color=color, linewidth=.5)
        ax.fill_between(ax.get_xlim(), [mn+.5, mn+.5], [mn-.5, mn-.5], color='k', alpha=.1)
        ax.fill_betweenx(ax.get_ylim(), [mn+.5, mn+.5], [mn-.5, mn-.5], color='k', alpha=.1)
