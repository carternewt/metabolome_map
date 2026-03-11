#!/bin/bash
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=carter.newton@uga.edu
#SBATCH --mail-type=START,END,FAIL
#SBATCH --error=/work/lylab/cjn40747/metabolome/logs/%j.err
#SBATCH --output=/work/lylab/cjn40747/metabolome/logs/%j.out

OUT='/work/lylab/cjn40747/metabolome'
HOME='/home/cjn40747/metabolome_map'

ml purge
ml CarveMe/1.6.6

carve -o $OUT/F7_5.xml -g LB -v --fbc2 $OUT/F7_5.faa

ml purge 

source activate cobra
python cobra_analysis.py