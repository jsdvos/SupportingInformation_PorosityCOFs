import os
import numpy as np
import h5py

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.stats import gaussian_kde

from yaff import System, CVCOMProjection, log
log.set_level(0)
from molmod import MolecularGraph
from molmod.periodic import periodic as pt
from molmod.units import angstrom, kelvin
from molmod.constants import boltzmann

def init_system(fn_h5):
    system = System.from_file(fn_h5).supercell(1,1,1)
    graph = MolecularGraph(system.bonds, system.numbers)
    indices = graph.independent_vertices
    z_orders = {}
    for index in indices:
        zs = [system.pos[i][2] for i in index]
        z_orders[tuple(index)] = np.mean(zs)
    sorted_indices = [np.array(list(i)) for i in sorted(z_orders.keys(), key = lambda x: z_orders[x])]
    return system, sorted_indices

def get_tfb_content(system, indices):
    n_tp = 0.0
    n_tfb = 0.0
    ffatypes = []
    count = 0
    for i in indices:
        ffatype = system.get_ffatype(i)
        if ffatype == 'O_TP':
            n_tp += 1
            if count % 3 == 0:
                ffatypes.append(ffatype)
            count += 1
        elif ffatype == 'H_C_C2_TPB':
            n_tfb += 1
            if count % 3 == 0:
                ffatypes.append(ffatype)
            count += 1
    if len(ffatypes) == 16:
        tp_ = 0
        tfb_ = 0
        mix_ = 0
        for i in range(8):
            if ffatypes[i] == ffatypes[i+8] == 'O_TP':
                tp_ += 1
            elif ffatypes[i] == ffatypes[i+8] == 'H_C_C2_TPB':
                tfb_ += 1
            else:
                mix_ += 1
        return float(tfb_)/8
    else:
        return float(n_tfb)/float(n_tfb + n_tp)



if __name__ == '__main__':
    for temp in [77]:
        src = '../../MDrun/{}K'.format(temp)
        for struct in os.listdir(src):
            # Load trajectory and atom indices of every layer
            fn_h5 = os.path.join(src, struct, 'traj.h5')
            system, sorted_indices = init_system(fn_h5)
            
            # Calculate interlayer distance (z) and layer offset (d) for every sampled frame
            z = [] # Interlayer distance
            d = [] # Layer offset
            tfb = [] # nuber of TFB-pairs in neighboring layers
            f = h5py.File(fn_h5, 'r')
            for count in range(7500, 15001, 15):
                # Get system geometry
                pos = f['trajectory']['pos'][count]
                rvecs = f['trajectory']['cell'][count]
                system.pos[:] = pos
                system.cell.update_rvecs(rvecs)

                # Iterate over neighboring layers
                for i in range(12):
                    j = i + 1
                    if j == 12:
                        # Correct for the last layer
                        sys = system.supercell(1,1,2)
                        j -= 12
                        indices_i = sorted_indices[i]
                        indices_j = sorted_indices[j] + system.natom
                    else:
                        sys = system
                        indices_i = sorted_indices[i]
                        indices_j = sorted_indices[j]

                    tfb_content = get_tfb_content(system, np.concatenate([sorted_indices[i], sorted_indices[j]]))
                    x = CVCOMProjection(sys, groups = [indices_i, indices_j], index = 0).compute()
                    y = CVCOMProjection(sys, groups = [indices_i, indices_j], index = 1).compute()
                    z.append(CVCOMProjection(sys, groups = [indices_i, indices_j], index = 2).compute())
                    d.append(np.sqrt(x**2 + y**2))
                    tfb.append(tfb_content) # Number of TFB-pairs in neighboring layers (i/8, i: 0->8)
            
            z = np.array(z)
            d = np.array(d)
            tfb = np.array(tfb)


            # Plot distributions of interlayer distance and layer offset
            cmap = plt.cm.get_cmap('viridis')
            eps = 0.01
            for key, values in [('z', z), ('d', d)]:
                bw = {'z': 0.025, 'd': 0.13}[key]
                x = np.linspace(min(values), max(values), 1000)
                ys = []
                count = 0
                # Create a distribution for each possible number of tfb-pairs
                for i in range(9):
                    inds = np.where((tfb > float(i)/8-eps) & (tfb < float(i)/8+eps))
                    count += len(inds[0])
                    data = values[inds]
                    if len(data) == 0:
                        ys.append(np.full(x.shape, -6012*12)) # Plot out of figure box but include for cmap normalization
                    else:
                        density = gaussian_kde(data, bw_method = bw/data.std())
                        to_plot = density(x)
                        to_plot *= len(inds[0]) # Normalize according to number of occurences
                        ys.append(to_plot)
                assert count == len(tfb) == len(z) == len(d)
                max_y = np.max(ys)
                segs = [np.column_stack([x/angstrom, y/count]) for y in ys]

                # Plot the distributions
                fig, ax = plt.subplots()
                xmin = np.min(x/angstrom)
                xmax = np.max(x/angstrom)
                x_range = xmax - xmin
                if key == 'd':
                    ax.set_xlim(1, 5)
                elif key == 'z':
                    ax.set_xlim(2.7, 3.5)
                line_segments = LineCollection(segs, array = np.array([100*float(i)/8 for i in range(9)]))
                ax.add_collection(line_segments)

                # Overall distribution
                density = gaussian_kde(values, bw_method = bw/values.std())
                to_plot = density(x)
                plt.plot(x/angstrom, to_plot, 'k--')


                ymax = max(to_plot)
                ax.set_ylim(0, ymax*1.05)
                if key == 'z':
                    ax.set_xlabel(r'Interlayer distance [$\AA$]')
                elif key == 'd':
                    ax.set_xlabel(r'Layer offset [$\AA$]')
                ax.set_ylabel('Probability')
                ax.set_yticks([])
                ax.set_title('TFB content = {:.2f}%'.format(100*get_tfb_content(system, np.arange(system.natom))))
                plt.savefig('../figs/LayerDistributions/{}_{}_{}K.pdf'.format(struct, key, temp), bbox_inches = 'tight')
                axcb = fig.colorbar(line_segments)
                axcb.set_label('Percentage of TFB pairs per 2 layers [%]')
                plt.savefig('../figs/LayerDistributions/{}_{}_{}K_cmap.pdf'.format(struct, key, temp), bbox_inches = 'tight')
                plt.close()
