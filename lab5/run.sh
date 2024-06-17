#!/bin/bash

for i in 1 2 4 8 16
do
    echo -e "#################"
    echo -e "$i processes"
    echo -e "#################"
    
    mpirun -n $i ./"a.out" 270 271 900
    echo -e "\n--------------------"
done
