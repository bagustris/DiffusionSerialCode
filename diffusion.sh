#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#date
#/usr/bin/openssl speed
mpirun -np 2 diffusion.x 
#date
