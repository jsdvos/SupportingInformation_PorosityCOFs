import os

import pandas as pd
import matplotlib.pyplot as plt

def plot_pxrds(dfs, keys):
    plt.figure()
    xmax = min(max(dfs[key]['2theta']) for key in keys)

    for i, key in enumerate(keys):
        # Get correct PXRD format
        df = dfs[key]
        x = [x_ for x_ in df['2theta'] if x_ < xmax]
        y = df['Intensity'][:len(x)]
        y -= min(y)
        y /= max(y)
        y *= 0.95
        y += (len(keys)-i)*0.5

        # Get label
        if 'exp' in key:
            exp = True
            tp = float(key.split('_')[1][:-2])/100
        elif 'none' in key:
            exp = False
            tp = 1.0
        elif 'all' in key:
            exp = False
            tp = 0.0
        elif 'full' in key:
            exp = False
            tp = 1.0 - float(key[-2:])/96
        if 'exp' in key:
            label = 'Experimental (TP content = {:.2f}%)'.format(100*tp)
        else:
            label = 'Theoretical (TP content = {:.2f}%)'.format(100*tp)
        
        # Plot
        plt.plot(x, y, label = label)


    plt.xlabel(r'2$\theta$ [$^\circ$]')
    plt.ylabel('Intensity')
    plt.yticks([])
    plt.legend()
    figname = '_'.join(keys) + '.pdf'
    plt.savefig(os.path.join('../figs/PXRDs', figname))
    plt.close()

if __name__ == '__main__':
    # Load data
    dfs = {}
    for temp in ['exp', '293K']:
        if temp == 'exp':
            for struct in os.listdir(os.path.join('../data', temp)):
                if not struct.endswith('.dat'): continue
                df = pd.read_csv(os.path.join('../data', temp, struct))
                dfs['{}_{}'.format(temp, struct.split('.')[0])] = df
        else:
            for struct in os.listdir(os.path.join('../data', temp)):
                df = pd.read_csv(os.path.join('../data', temp, struct, 'output.dat'))
                if len(df['2theta']) == 0: continue
                dfs['{}_{}'.format(temp, struct)] = df
    
    # Plot PXRDs
    plot_pxrds(dfs, ['exp_0TP', '293K_ABCDEF_all'])
    plot_pxrds(dfs, ['exp_37TP', '293K_ABCDEF_full58'])
    plot_pxrds(dfs, ['exp_93TP', '293K_ABCDEF_full10'])
    plot_pxrds(dfs, ['exp_100TP', '293K_ABCDEF_none'])
