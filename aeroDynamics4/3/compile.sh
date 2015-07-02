#! /bin/sh

gfortran -o $1.exe $1.f90 && ./$1.exe
