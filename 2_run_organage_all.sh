#!/bin/bash
#SBATCH -c 2
#SBATCH --output=organage.out
#SBATCH --error=organage.err
#SBATCH -p high_p
#SBATCH --mem=10G

#setwd("/data/Epic/subprojects/Somalogic/work/Oliver")

echo 'Hello'
python3.9 --version
python3.9 py_script_organage_EPIC_all.py



# system("sbatch run_organage_all.sh")
# system("squeue -u oliver.robinson")

