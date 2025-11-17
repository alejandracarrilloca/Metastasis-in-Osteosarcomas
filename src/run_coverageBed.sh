  GNU nano 5.6.1                                      run_coverageBed.sh
#!/bin/bash
#  Este script calcula las coberturas de cada archivo bam
#  Recibe el  archivo bam  y su archivo gff con las anotaciones
#  ejemplo de ejecucion:
# qsub  -l h_vmem=8G /home/vjimenez/bin/run_coverageBed.sh dmel-all-r6.65.OnlyGenes.gff BWA/sampleGSM461177.sort.bam

echo "inicia trabajando con: $2 ";date

echo "coverageBed -a $1 –b  $2 >  ${2/bam/count.txt}"
coverageBed -a $1 -b  $2 >  ${2/bam/count.txt}

echo "termina con salida en ${2/bam/count.txt} ";date
