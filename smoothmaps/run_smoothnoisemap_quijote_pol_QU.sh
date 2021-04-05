#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l walltime=21:00:00:00
#PBS -N qt_smoothmaps
#PBS -m abe
#PBS -M michael.peel@manchester.ac.uk

source ~/.bashrc

cd /home/mpeel/git/astrocode/smoothmaps/

echo "Running job in"
echo `pwd`

python run_smoothnoisemap_quijote_pol_QU.py
