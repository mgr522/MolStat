#!/bin/bash
# Script to compile, build and run code using cvode package

g++ -I"/usr/local/include" -c $1 -o $2.o
g++ -o $2 $2.o "/usr/local/lib/libsundials_cvode.a" "/usr/local/lib/libsundials_nvecserial.a" -lm -lgsl -lgslcblas
./$2
