#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00:00
#PBS -N combinemaps
#PBS -m abe
#PBS -M michael.peel@manchester.ac.uk

source ~/.bashrc

cd /home/mpeel/git/astrocode/smoothmaps/

echo "Running job in"
echo `pwd`

# python combine_smoothed_maps_20.py
python combine_smoothed_maps.py
