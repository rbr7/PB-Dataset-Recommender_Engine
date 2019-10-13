#!/bin/bash
cat ./modules/*.*gmt > temp1
n_modules=$(cat temp1|wc -l)
cat $1 temp1 > temp2
mv temp2 $1
rm temp1
echo $n_modules "new modules have been added"
