#!/usr/bin/awk -f
BEGIN {
    OFS="\t"
}
{
    name=NF>3?$4:"."
    score=NF>4?$5:"."
    if ($2 > 0){
        print $1, $2-1, $2, name, score, "+"
    }
    print $1, $3-1, $3, name, score, "-"
}
