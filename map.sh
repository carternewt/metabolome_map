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

find $OUT/type_strains/ -name *.faa -type f | while read -r file; do
    base=$(basename "$file" .faa)
    carve -o $OUT/"$base"_GG.xml -v -g GG --mediadb $HOME/GG_medium.tsv --fbc2 $file
    carve -o $OUT/"$base"_CDB.xml -v -g CDB --mediadb $HOME/GG_medium.tsv --fbc2 $file
done

ml purge 
source activate cobra

find $OUT -name *_CDB.xml -type f | while read -r file; do
    python cobra_analysis_CDB.py $file
done

find $OUT -name *_GG.xml -type f | while read -r file; do
    python cobra_analysis_GG.py $file
done