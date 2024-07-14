#!/bin/bash

input=$1
output=$2
sample=$3
log=$4


#sed  '/^>/ s/\//_/g; /^>/ s/$/_${sample}/' $input > $output 2> ${log}
#echo "gene_trans_map were created successfully!"

#!/bin/bash

sed  "/^>/ s/\//_/g; /^>/ s/$/_${sample}/" $input > $output 2> $log
echo "Fasta header were modified successfully!"
