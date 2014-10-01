#!/bin/bash
# Script to compile, build and run code using cvode package

g++ -I"/usr/local/include" -c $1 -o $2.o
g++ -o $2 $2.o -lsundials_cvode -lsundials_nvecserial -lm -lgsl -lgslcblas
./$2
