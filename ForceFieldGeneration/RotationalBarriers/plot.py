import os

import numpy as np
import pandas as pd
import h5py as h5
import matplotlib.pyplot as plt

from molmod.units import kjmol

def get_ai_scan(folder, label):
    table_label = 'table_{}'.format('_'.join(label.split('_')[:2]))
    fn = os.path.join(folder, label, '{}_hdf5/{}/pyirontable.csv'.format(table_label, table_label))
    df = pd.read_csv(fn)
    return df

def get_ff_scan(df, folder, label):
    for key in ['nodih', 'wdih', 'polysix']:
        scan_label = 'scan_{}_{}'.format('_'.join(label.split('_')[:2]), key)
        fn = os.path.join(folder, '{}_hdf5/{}/output.h5'.format(scan_label, scan_label))
        f = h5.File(fn, mode = 'r')
        df[key] = f['trajectory']['epot'][:]
    return df

def process_data(df):
    df['old_ff'] = df['wdih']
    df['new_ff'] = df['nodih'] + df['polysix']
    for key in ['E', 'nodih', 'wdih', 'polysix', 'old_ff', 'new_ff']:
        df[key] -= min(df[key])

def plot(df, figname):
    plt.figure()
    result = {i: [df[key].values[i] for key in ['old_ff', 'new_ff', 'E']] for i in range(len(df))}
    angles = np.arange(0, 181)
    old_ff = []
    new_ff = []
    ai = []
    for i in angles:
        if i <= 90:
            old_ff.append(result[i][0])
            new_ff.append(result[i][1])
            ai.append(result[i][2])
        else:
            old_ff.append(result[180-i][0])
            new_ff.append(result[180-i][1])
            ai.append(result[180-i][2])
    old_ff = np.array(old_ff)
    new_ff = np.array(new_ff)
    ai = np.array(ai)
    plt.plot(angles, old_ff/kjmol, c = 'C0', label = 'Initial FF')
    plt.plot(angles, new_ff/kjmol, c = 'C1', label = 'Fitted FF')
    plt.plot(angles, ai/kjmol, c = 'C3', linestyle = '--', label = 'AI reference')
    plt.ylabel('Energy [kJ/mol]')
    plt.xlabel(r'Dihedral angle [$^\circ$]')
    plt.legend()
    plt.xlim(0, 180)
    ymin, ymax = plt.ylim()
    plt.ylim(0, ymax)
    plt.savefig(figname, bbox_inches = 'tight')
    plt.clf()

if __name__ == '__main__':
    folder = '/arcanine/scratch/gent/vo/000/gvo00003/vsc42354/BACKUPs/Rasha/barriers/rasha'
    for label in ['PA11_amine_2_0_10_11', 'PA12_amine_3_1_18_19', 'TP_amine_9_12_18_29', 'PA12_imine_2_0_10_11', 'PA22_imine_2_0_10_11', 'TPB_imine_15_16_13_12']:
        df = get_ai_scan(folder, label)
        df = get_ff_scan(df, folder, label)
        process_data(df)
        print(label)
        print(df)
        plot(df, 'figs/{}.pdf'.format(label))

