#!/bin/bash
cat $1 module-info.txt  > temp
mv temp $1
n_datasets=$(cat module-info.txt|wc -l)
echo $n_datasets "datasets have been added"
