import numpy as np

from yaff import System, CVCOMProjection, log
from yaff.pes.ff import ForceField
from yaff.sampling.dof import StrainCellDOF
from yaff.sampling.opt import CGOptimizer
from yaff.pes.ext import Cell
#log.set_level(0)
from molmod import MolecularGraph
#from molmod.periodic import periodic as pt
from molmod.units import angstrom, kjmol
#from molmod.constants import boltzmann

# Step 1: create 2layer system

for struct in ['ABCDEF_all', 'ABCDEF_none']:
    system = System.from_file('../../MolecularDynamics/MDanalysis/calculations/0K/{}/frame0/frame0.chk'.format(struct))
    graph = MolecularGraph(system.bonds, system.numbers)
    indices = graph.independent_vertices
    z_orders = {}
    for index in indices:
        zs = [system.pos[i][2] for i in index]
        z_orders[tuple(index)] = np.mean(zs)
    sorted_indices = [np.array(list(i)) for i in sorted(z_orders.keys(), key = lambda x: z_orders[x])]
    
    if struct == 'ABCDEF_all':
        subsys = system.subsystem(np.concatenate([sorted_indices[1], sorted_indices[2]]))
    elif struct == 'ABCDEF_none':
        subsys = system.subsystem(np.concatenate([sorted_indices[0], sorted_indices[1]]))
    rvecs = subsys.cell.rvecs.copy()
    rvecs[2, 2] /= 6
    subsys.cell = Cell(rvecs)
    subsys.to_file('../data/{}/nonopt_2layers.chk'.format(struct))

    pars = '../data/{}/pars.txt'.format(struct)
    ff_kwargs = {
            'rcut': 15*angstrom,
            'alpha_scale': 3.2,
            'gcut_scale': 1.5,
            'smooth_ei': True,
            'tailcorrections': True
    }
    ff = ForceField.generate(subsys, pars, **ff_kwargs)
    dof = StrainCellDOF(ff)
    opt = CGOptimizer(dof)
    while True:
        opt.run(100000)
        if opt.dof.converged:
            opt.dof.ff.system.to_file('../data/{}/opt_2layers.chk'.format(struct))
            break

