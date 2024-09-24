#!/bin/bash

# Prepare diamond index for database 
INPUT=$1
OUTPUT=$2 

diamond makedb -p 8 --in $IN -d $OUT

echo "Diamond index created for $INPUT."
