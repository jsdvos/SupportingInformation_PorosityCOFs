#!/bin/sh
#
#PBS -N _md_ABCDEF_full19_293K
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -m n

date

# Set up input
ORIGDIR=$PBS_O_WORKDIR
WORKDIR=/local/$PBS_JOBID

if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR; fi
cd $WORKDIR

cp ${ORIGDIR}/md.py $WORKDIR
cp ${ORIGDIR}/pars.txt $WORKDIR
cp ${ORIGDIR}/init.chk $WORKDIR

# Copy back results every half hour
( while true; do
        sleep 1800
        cp ${WORKDIR}/traj.h5 ${ORIGDIR}
        cp ${WORKDIR}/md.log ${ORIGDIR}
  done ) &


# Load modules
module load LAMMPS/3Mar2020-foss-2019b-Python-3.7.4-kokkos # Load LAMMPS, also loads yaff

# Run
python md.py > md.log

# Copy back results
cp ${WORKDIR}/traj.h5 ${ORIGDIR}
cp ${WORKDIR}/md.log ${ORIGDIR}

# Finalize
rm -rf $WORKDIR

date

