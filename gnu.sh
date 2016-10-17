#!/bin/bash

gfortran tp1.f90
./a.out
gnuplot -e "set style data boxes; p 'histogram.dat'; pause -1"

exit 0
