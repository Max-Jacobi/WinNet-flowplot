from ... import os, np, St_N, St_Z

def TimeDensTempBox(ax, file, **kwargs):
    '''
    plot a box with metainfo from flowfile or snapshot 
    '''
    text =  r'time: {:5.2e}$\,$ms'+'\n'
    text += r'dens: {:5.2e}$\,$g/cm$^3$'+'\n'
    text += r'temp: {:5.2f}$\,$GK'
    text = text.format(file.time*1e-3, file.dens, file.temp)
    text = ax.text(0.02, 0.98, text,
                   horizontalalignment='left',
                   verticalalignment='top',
                   transform=ax.transAxes,
                   bbox=dict(facecolor='w', edgecolor='k', boxstyle='round'),
                   **kwargs)
    return text

def plotMagicNumbers(ax, x=True, y=True, alpha=.2, color='k', lw=.1):
    '''
    routine for ploting lines at magic numbers
    '''
    ixymagic = [2, 8, 20, 28, 50, 82, 126]
    for mn in ixymagic:
        if y:
            ax.axhline(y=mn-.5, color=color, lw=lw)
            ax.axhline(y=mn+.5, color=color, lw=lw)
            ax.fill_between(ax.get_xlim(), [mn+.5, mn+.5], [mn-.5, mn-.5], color=color, alpha=.3*alpha)
        if x:
            ax.axvline(x=mn-.5, color=color, lw=lw)
            ax.axvline(x=mn+.5, color=color, lw=lw)
            ax.fill_betweenx(ax.get_ylim(), [mn+.5, mn+.5], [mn-.5, mn-.5], color=color, alpha=.3*alpha)

def plotStableIsotopes(ax, **kwargs):
    '''
    routine for ploting the positions of stable isotopes
    '''
    standart = {'color': 'k',
                'marker': 'o'}
    for k, s in standart.items():
        if k not in kwargs.keys():
            kwargs[k] = s
        
    ax.scatter(St_N, St_Z, **kwargs)
    return

def ScanDir(path):
    '''
    Scans a directory for snapshots or flowfiles
    Arguments:
    path    - path to directory
    returns (filepaths, filetimes, filepaths)
    '''
    filepaths = []
    filetimes = []
    filetemps = []
    for file in os.listdir(path):
        with open('{}/{}'.format(path, file), 'r') as ff:
            header = ff.readline()
            if 'dt' in header:
                header = ff.readline()
                time, dt, temp, dens = np.array(header.split()).astype(float)
            else:
                header = ff.readline()
                time, temp, dens = np.array(header.split()).astype(float)
        filetimes.append(time)
        filetemps.append(temp)
        filepaths.append('{}/{}'.format(path, file))
        sort = np.argsort(filetimes)[::-1]
        filepaths=np.array(filepaths)[sort]
        filetimes=np.array(filetimes)[sort]
        filetemps=np.array(filetemps)[sort]
    return filepaths, filetimes, filepaths
