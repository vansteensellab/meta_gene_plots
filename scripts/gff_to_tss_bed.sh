#!/usr/bin/env sh
## select all transcript start sites
is_zipped=$(file $1 | grep 'gzip' | wc -l)

if [ $is_zipped -eq 0 ]
then
    cat="cat"
else
    cat="zcat"
fi

$cat $1 | \
    sed -rn "s/ID=([^;]+)/\1\t/p" | sed -rn "s/gene_id=([^;]+)/\t\1\t/p" | \
    awk -vOFS='\t' '{
            if ($3=="transcript"){
                if ($7 == "+"){
                    tss=$4 - 1
                } else {
                    tss=$5 - 1
                }
                print $1 , tss , tss + 1, $9 , "." , $7 , $11
            }
        }' | bedtools sort -i - | gzip -c
