#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#date
#/usr/bin/openssl speed
mpirun -np 2 Diffusion_serial_code/diffusion.x 
#date
