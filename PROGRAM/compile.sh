#!/bin/bash

#Compliler for C and Fortran 
comp1='gcc'
comp2='gfortran'

# Optimization Flags



# Warning and Error flags gfortran
gflags_hard='-Wno-tabs -Wall -Wextra -Warray-temporaries -fbounds-check -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan'
gflags_soft='-Wall -fbounds-check -Wno-tabs'
# Warning and Error flags intel compiler
iflags_hard='-check all -fpe0 -warn -traceback -debug extended'
iflags_soft=''

flags=
#$gflags_hard
#opt=
opt='-O'

$comp1 $opt -c secs.c
$comp2 $opt -c $flags ran2.f
$comp2 $opt -c $flags r1279.f90
$comp2 $opt -c $flags read_data.f90
$comp2 $opt -c $flags init_vars.f90
$comp2 $opt -c $flags routines.f90
$comp2 $opt -c $flags main.f90
#$comp2 mainSQA.o secs.o ran2.o r1279.o initialization.o init_vars.o initial_config.o thermo.o  $opt $flags -o r_mainSQA
$comp2 main.o secs.o ran2.o r1279.o read_data.o init_vars.o routines.o  $opt $flags -o r_main
##gfortran mainSQA.o dummy.x -o r_mainSQA
rm *.o
rm *.mod