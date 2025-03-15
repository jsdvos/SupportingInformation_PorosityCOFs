from yaff import System
from yaff.pes.ff import ForceField
from yaff.sampling.dof import CartesianDOF, StrainCellDOF
from yaff.sampling.opt import CGOptimizer

from molmod.units import angstrom

sys = System.from_file('ABCDEF_full29.chk')
ff = ForceField.generate(sys, 'pars.txt', rcut = 15.0*angstrom,
        alpha_scale = 2.86, gcut_scale = 1.0, smooth_ei = True, tailcorrections = True)
dof = StrainCellDOF(ff)
opt = CGOptimizer(dof)
while True:
    opt.run(10000)
    opt.dof.ff.system.to_file('ABCDEF_full29_save{}.chk'.format(opt.counter))
    if opt.dof.converged:
        opt.dof.ff.system.to_file('ABCDEF_full29_opt.chk')
        break
    if not opt.counter % 10000 == 0:
        break

