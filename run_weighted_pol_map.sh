#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00:00
#PBS -N weightmap
#PBS -m abe
#PBS -M michael.peel@manchester.ac.uk

source ~/.bashrc

cd /home/mpeel/git/astrocode/smoothmaps/

echo "Running job in"
echo `pwd`

python run_weighted_pol_map.py
