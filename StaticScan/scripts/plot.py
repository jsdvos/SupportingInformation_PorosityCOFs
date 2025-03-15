import os

import numpy as np

import matplotlib.pyplot as plt

from yaff import System, CVCOMProjection, log
from yaff.pes.ff import ForceField
from yaff.sampling.dof import StrainCellDOF
from yaff.sampling.opt import CGOptimizer
from yaff.pes.ext import Cell
log.set_level(0)
from molmod import MolecularGraph
#from molmod.periodic import periodic as pt
from molmod.units import angstrom, kjmol
#from molmod.constants import boltzmann

# Step 2: scan 2layer system

for cof in ['ABCDEF_none', 'ABCDEF_all']:
    folder = '../data/{}'.format(cof)
    for fn in os.listdir(folder):
        # Load optimized system
        struct = fn.split('.')[0]
        if not (fn.startswith('opt_2layers') and fn.endswith('.chk')): continue
        fn = os.path.join(folder, fn)
        system = System.from_file(fn)
        graph = MolecularGraph(system.bonds, system.numbers)
        indices = graph.independent_vertices
        z_orders = {}
        for index in indices:
            zs = [system.pos[i][2] for i in index]
            z_orders[tuple(index)] = np.mean(zs)
        sorted_indices = [np.array(list(i)) for i in sorted(z_orders.keys(), key = lambda x: z_orders[x])]

        # Load FF
        pars = os.path.join(folder, 'pars.txt')
        ff_kwargs = {
            'rcut': 15*angstrom,
            'alpha_scale': 3.2,
            'gcut_scale': 1.5,
            'smooth_ei': True,
            'tailcorrections': True
        }
        ff = ForceField.generate(system, pars, **ff_kwargs)

        # Store equilibrium geometry
        pos0 = system.pos.copy()
        rvecs0 = system.cell.rvecs.copy()

        es = []
        e_covs = []
        e_eis = []
        e_vdws = []
        zs = []
        ds = []
        # Iterate over different interlayer distances
        for dz in np.arange(-1.0, 10.0, 0.1):
            # Change system geometry
            pos = pos0.copy()
            for i in sorted_indices[1]:
                pos[i][2] += dz
            rvecs = rvecs0.copy()
            rvecs[2, 2] += 2*dz

            # Calculate FF energy
            ff.update_pos(pos)
            ff.update_rvecs(rvecs)
            e = ff.compute()
            e_cov = 0.0
            e_ei = 0.0
            e_vdw = 0.0
            for part in ff.parts:
                if 'valence' in part.name:
                    e_cov += part.energy
                elif 'ei' in part.name or 'ewald' in part.name:
                    e_ei += part.energy
                elif 'mm3' in part.name:
                    e_vdw += part.energy
                else:
                    raise NotImplementedError('{} not recognized'.format(part.name))
            es.append(e)
            e_covs.append(e_cov)
            e_eis.append(e_ei)
            e_vdws.append(e_vdw)
            
            # Calculate layer characteristics (layer offset and interlayer distance)
            z = []
            d = []
            for i in range(2):
                j = i + 1
                if j == 2:
                    sys = ff.system.supercell(1,1,2)
                    j -= 2
                    indices_i = sorted_indices[i]
                    indices_j = sorted_indices[j] + system.natom
                else:
                    sys = ff.system
                    indices_i = sorted_indices[i]
                    indices_j = sorted_indices[j]
            
                x = CVCOMProjection(sys, groups = [indices_i, indices_j], index = 0).compute()
                y = CVCOMProjection(sys, groups = [indices_i, indices_j], index = 1).compute()
                z.append(CVCOMProjection(sys, groups = [indices_i, indices_j], index = 2).compute())
                d.append(np.sqrt(x**2 + y**2))
                ds.append(np.sqrt(x**2 + y**2))
            
            zs.append(np.mean(z))

        # Convert to numpy arrays
        zs = np.array(zs)
        e_covs = np.array(e_covs)
        e_eis = np.array(e_eis)
        e_vdws = np.array(e_vdws)
        es = np.array(es)

        # Shift energies to become zero at "infinity"
        e_covs -= e_covs[-1]
        e_eis -= e_eis[-1]
        e_vdws -= e_vdws[-1]
        es -= es[-1]

        # Plot
        fig, ax = plt.subplots()
        ax.plot(zs/angstrom, e_eis/kjmol, label = 'Electrostatic energy')
        ax.plot(zs/angstrom, e_vdws/kjmol, label = 'Van der Waals energy')
        ax.plot(zs/angstrom, es/kjmol, 'k--', label = 'Total energy')
        ax.set_xlabel(r'Interlayer distance [$\AA$]')
        ax.set_ylabel(r'Energy [kJ/mol]')
        ax.legend()
        ax.set_title('Layer offset = {:.3f}A'.format(np.mean(ds)/angstrom))
        fig.savefig('../figs/{}_z_{}.pdf'.format(cof, struct))

