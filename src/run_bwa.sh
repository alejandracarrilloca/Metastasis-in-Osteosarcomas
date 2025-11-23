#!/bin/bash
bwa mem -t 10 -o $4 $1 $2 $3
# $1: Index generado
# $2: Archivo 1
# $3: Archivo 2
# $4: Salida
