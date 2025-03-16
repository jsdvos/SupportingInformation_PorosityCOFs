import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from yaff import System, log, angstrom, nanometer
log.set_level(0)

def plot_geo(dfs, tp, keys, ylabel = None, figname = None):
    fig, ax = plt.subplots()
    data = [dfs[key][0] for key in keys]
    positions = [tp['_'.join(key.split('_')[1:])] for key in keys]
    data = [x for _, x in sorted(zip(positions, data))]
    positions = np.array([x for x, _ in sorted(zip(positions, data))])
    parts = ax.violinplot(data, positions = 100*positions, showextrema = False, widths = 5)
    ax.plot(100*positions, np.mean(data, axis = 1), 'o-', c = 'k')
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor('C0')
        pc.set_edgecolor('black')
        pc.set_alpha(0.5)
    ax.set_xlabel('TP content [%]')
    if ylabel:
        ax.set_ylabel(ylabel)
    fig.savefig(os.path.join('../figs/GeometryDistributions', figname))
    plt.close()

if __name__ == '__main__':
    # Get TFB content of all structures
    tp_content = {}
    for struct in os.listdir('../calculations/0K'):
        sys = System.from_file('../calculations/0K/{}/frame0/frame0.chk'.format(struct))
        n_tp = 0.0
        n_tfb = 0.0
        for ffatype_id, ffatype in enumerate(sys.ffatypes):
            if ffatype == 'O_TP':
                n_tp = (sys.ffatype_ids == ffatype_id).sum()
            elif ffatype == 'H_C_C2_TPB':
                n_tfb = (sys.ffatype_ids == ffatype_id).sum()
        tp_content[struct] = float(n_tp)/float(n_tp + n_tfb)
    
    # Iterate over properties to plot
    for key in ['asa_m', 'av_m', 'volume', 'rho']:
        # Load data
        dfs = {}
        for temp in os.listdir('../data'):
            if temp == 'exp': continue
            for struct in os.listdir(os.path.join('../data', temp)):
                df = pd.read_csv(os.path.join('../data', temp, struct, key + '.dat'), header = None, delim_whitespace = True)
                if key == 'volume':
                    df[0] = df[0]*angstrom**3/(nanometer)**3
                dfs['{}_{}'.format(temp, struct)] = df
        
        # Make figures for each temperature
        for temp in [77, 293, 400, 500]:
            keys = [k for k in dfs.keys() if k.startswith(str(temp))]
            ylabel = {
                    'asa_m': r'Gravimetric accessible surface area [m$^2$/g]',
                    'av_m': r'Gravimetric pore volume [cm$^3$/g]',
                    'volume': r'Unit cell volume [nm$^3$]',
                    'rho': 'Mass density [g/cm$^3$]',
                    }[key]
            plot_geo(dfs, tp_content, keys, ylabel = ylabel, figname = '{}_{}K.pdf'.format(key, temp))

