#!/bin/bash

prog="omp-lab"
cpp=".cpp"
out=".out"

static="static"
dynamic="dynamic"
guided="guided"

export OMP_SCHEDULE
export OMP_NUM_THREADS

g++ -fopenmp -o $prog$out $prog$cpp

for mode in $static $dynamic $guided
do
    OMP_SCHEDULE=$mode
    echo -e "#################"
    echo -e "$mode"
    echo -e "#################"
    for i in 1 2 4 8 12 16
    do
        echo -e "program: $i threads"
        OMP_NUM_THREADS=$i
        ./$prog$out 
        echo -e "\n--------------------"
    done
done

export OMP_SCHEDULE
export OMP_NUM_THREADS=12

for mode in $static $dynamic $guided
do
    echo -e "#################"
    echo -e "$mode"
    echo -e "#################"
    for i in 1 2 4 8 16 100
    do
        echo -e "program: $i chunks"
        OMP_SCHEDULE=$mode,$i
        ./$prog$out 
        echo -e "\n--------------------"
    done
done

