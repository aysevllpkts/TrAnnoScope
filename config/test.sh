#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time 12:00:00
#SBATCH --partition short
#SBATCH --mem 1G
#SBATCH --account MolGen 
#SBATCH -J subsetting
#SBATCH -e subsetting_%A.e
#SBATCH -o subsetting_%A.o

echo "Hello world!"
