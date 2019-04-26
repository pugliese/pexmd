#!/bin/bash

start=$SECONDS

/home/pugliese/LJ/fermi_gas_red.e $1 $2

duration=$(( SECONDS - start ))

echo -e "$duration" >> wall_time.file

notify facupugliese95@gmail.com Pauli D*=$1 rho*=$2 terminado
