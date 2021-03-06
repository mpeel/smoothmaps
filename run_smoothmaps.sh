#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l walltime=7:00:00:00
#PBS -N smoothmaps
#PBS -m abe
#PBS -M michael.peel@manchester.ac.uk

source ~/.bashrc

cd /home/mpeel/git/astrocode/smoothmaps/

echo "Running job in"
echo `pwd`

# python run_smoothmaps_wmap_pol.py
python run_smoothmaps_quijote.py
# python run_smoothmaps_planck2015_pol.py
# python run_smoothmaps_planck2018_pol.py
# python run_smoothmaps_planck2020_pol.py
