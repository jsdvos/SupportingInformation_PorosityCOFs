#!/bin/sh
#
#PBS -N _opt_ABCDEF_full77
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -m n

date

# Set up input
ORIGDIR=$PBS_O_WORKDIR
cd $ORIGDIR

# Load modules
module load LAMMPS/3Mar2020-foss-2019b-Python-3.7.4-kokkos # Load LAMMPS, also loads yaff

# Run
python opt.py > opt.log

date

