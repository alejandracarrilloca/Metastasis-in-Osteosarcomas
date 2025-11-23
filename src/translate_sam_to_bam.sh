#!/bin/bash
# Este script convierte de SAM a BAM, ordena y calcula las estadísticas de ensamblado
# Recibe el archivo SAM a convertir
# ejemplo de ejecución:
# qsub  -l h_vmem=8G /home/vjimenez/bin/translate_sam_to_bam.sh

echo "inicia sam to bam: $1 "; date

# Convertir SAM a BAM y ordenar
samtools view -bS "$1" | samtools sort - -o "${1/.sam/.sort.bam}"

# Si quieres conservar el SAM, NO lo borramos
# if [ -e "${1/.sam/.sort.bam}" ]; then
#     rm "$1"
# fi

echo "termina sam to bam $1 "; date
