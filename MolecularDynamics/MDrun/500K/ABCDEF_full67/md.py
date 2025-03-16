import numpy as np
import h5py as h5

from yaff import System, log
from yaff.pes.ff import ForceField
from yaff.external.liblammps import swap_noncovalent_lammps

from yaff.sampling.io import HDF5Writer
from yaff.sampling.verlet import VerletIntegrator, VerletScreenLog
from yaff.sampling.nvt import NHCThermostat
from yaff.sampling.npt import MTKBarostat, TBCombination

from molmod.units import angstrom, femtosecond, kelvin, bar, atm, kjmol

mpi = False
if mpi:
    from mpi4py import MPI

    # Setup MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    # Turn off logging for all processes, except one
    log.set_level(log.silent)
    if rank==0: log.set_level(log.medium)

def get_ff_kwargs():
    ff_kwargs = {
            'rcut': 15*angstrom,
            'alpha_scale': 3.2,
            'gcut_scale': 1.5,
            'smooth_ei': True,
            'tailcorrections': True
    }
    return ff_kwargs

def load_ff(sys, fn_pars, use_lammps = True):
    ff_kwargs = get_ff_kwargs()
    ff = ForceField.generate(sys, fn_pars, **ff_kwargs)
    if use_lammps:
        fn_sys = 'system.dat' # LAMMPS System file
        fn_table = 'table.dat' # LAMMPS force field tabulation file
        # Tabulate the non-bonded interactions
        # Bonded interactions remain calculated by Yaff
        try:
            if mpi:
                ff_lammps = swap_noncovalent_lammps(ff, fn_system = fn_sys,
                        fn_table = fn_table, comm = comm)
            else:
                ff_lammps = swap_noncovalent_lammps(ff, fn_system = fn_sys,
                        fn_table = fn_table)
        except ValueError:
            sys_large = sys.supercell(2,2,2)
            ff_large = ForceField.generate(sys_large, fn_pars, **ff_kwargs)
            if mpi:
                ff_large_lammps = swap_noncovalent_lammps(ff_large,
                        fn_system = 'system_large.dat', fn_table = fn_table, comm = comm)
            else:
                ff_large_lammps = swap_noncovalent_lammps(ff_large,
                        fn_system = 'system_large.dat', fn_table = fn_table)
            ff_lammps = swap_noncovalent_lammps(ff, fn_system = fn_sys, fn_table = fn_table)
        gpos, vtens = np.zeros((sys.natom, 3)), np.zeros((3, 3))
        gpos_lammps, vtens_lammps = np.zeros((sys.natom, 3)), np.zeros((3, 3))
        e = ff.compute(gpos, vtens)
        e_lammps = ff_lammps.compute(gpos_lammps, vtens_lammps)
        p = np.trace(vtens)/3.0/ff.system.cell.volume
        p_lammps = np.trace(vtens_lammps)/3.0/ff.system.cell.volume
        print("E(Yaff) = %12.3f E(LAMMPS) = %12.3f deltaE = %12.3e kJ/mol"%(e/kjmol,e_lammps/kjmol,(e_lammps-e)/kjmol))
        print("P(Yaff) = %12.3f P(LAMMPS) = %12.3f deltaP = %12.3e bar"%(p/bar,p_lammps/bar,(p_lammps-p)/bar))
        return ff_lammps
    else:
        return ff

def load_integrator(ff, temp = 300.0*kelvin, press = 1*atm):
    vsl = VerletScreenLog(step = 20) # Print information to screen every 10 steps
    hdf5_writer = HDF5Writer(h5.File('traj.h5', mode='w'), step = 20) # HDF5 file with trajectory data
    thermo = NHCThermostat(temp, timecon=100*femtosecond) # Thermostat with timeconstant of 100fs
    baro = MTKBarostat(ff, temp = temp, press = press, timecon=1000*femtosecond) # Barostat
    tbc = TBCombination(thermo, baro) # Combine thermostat and barostat
    verlet = VerletIntegrator(ff, 0.5*femtosecond, temp0 = temp, hooks=[vsl, hdf5_writer, tbc]) # Verlet integrator
    return verlet

if __name__ == '__main__':
    sys = System.from_file('init.chk')
    ff = load_ff(sys, 'pars.txt')
    verlet = load_integrator(ff, temp = 500.0*kelvin)
    verlet.run(300000)


