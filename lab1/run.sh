#!/bin/bash

prog1="opp-lab1-seq"
prog2="opp-lab1-point-to-point"
prog3="opp-lab1-coll-comm"

mpicxx $prog1".cpp"
mpicxx $prog1".cpp" -o $prog1".out" 

mpicxx $prog2".cpp" 
mpicxx $prog2".cpp" -o $prog2".out"

mpicxx $prog3".cpp" 
mpicxx $prog3".cpp" -o $prog3".out"

echo -e "seq program"

./$prog1".out" 
./$prog1".out" 
./$prog1".out" 

for i in 2 4 8 16 24:
do
    echo -e "point to point program: $i processes"
    mpiexec -n $i ./$prog2".out" 
    mpiexec -n $i ./$prog2".out" 
    mpiexec -n $i ./$prog2".out" 

    echo -e "comm program: $i processes"
    mpiexec -n $i ./$prog3".out" 
    mpiexec -n $i ./$prog3".out" 
    mpiexec -n $i ./$prog3".out" 
done
