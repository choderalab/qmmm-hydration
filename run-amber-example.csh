#!/bin/tcsh
#  Batch script for mpirun job on cbio cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=08:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q gpu
#
# nodes: number of nodes
#   ppn: how many cores per node to use (must use one per GPU)
#   gpus: number of GPUs per node
#   shared keyword to allow GPUs to be switched to shared mode
#PBS -l nodes=1:ppn=2:gpus=2:shared
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N ambertest
#
# specify email
#PBS -M jchodera@gmail.com
#
# mail settings
#PBS -m n
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
# DOES NOT WORK WITH PBS -k oe FOR IMMEDIATE OUTPUT
#PBS -o /cbio/jclab/home/chodera/vvvr/openmm/timescale-correction/vvvr.out

cd /cbio/jclab/projects/collab/cisborn
source amber/amber.csh

cd ${PBS_O_WORKDIR}/molecules/benzene/solvent


date

setenv name system
#mpirun -np 1 $TCBin -UseMPI > mpi1.tc_job.dat & 
mpirun -np 1 sander.MPI -O -i $name.in -o $name.out -p $name.prmtop -c $name.inpcrd -x $name.mdcrd -r $name.rst -inf $name.info

date

