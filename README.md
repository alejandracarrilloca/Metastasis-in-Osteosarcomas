# Biological Activity of Biomarkers Associated With Metastasis in Osteosarcoma Cell Lines

26/08/2025

* Dulce Alejandra Carrillo Carlos [carrillo@lcg.unam.mx](mailto:carrillo@lcg.unam.mx)
* Addiel Antonio Platas Renteral [addielpr@lcg.unam.mx](mailto:addielpr@lcg.unam.mx)
* Andrea Villarruel García [andreavg@lcg.unam.mx](mailto:andreavg@lcg.unam.mx)

This analysis replicates the following article: [Biological Activity of Biomarkers Associated With Metastasis in Osteosarcoma Cell Lines](https://doi.org/10.1002/cam4.70391).

The data that support the findings of this study are openly available in **ENA EMBL-­EBI** ([https://www.ebi.ac.uk/ena/browser/home](https://www.ebi.ac.uk/ena/browser/home)) with Project IDs: **PRJEB65286**, **PRJNA698672**, and **PRJNA51801**.

---

## Planteamiento del problema

El osteosarcoma es el  cáncer más común entre adolescentes y adultos jóvenes. Aunque el tratamiento convencional incluye cirujía y quimioterapia, la **metástasis pulmonar** continua siendo La principal causa de muerte en pacientes con este tipo de cáncer. 

A nivel molecular, la metástasis implica una serie de cambios en la **expresión génica**, incluyendo vías relacionadas relacionadas con migración, adhesión y remodelación del microambiente tumoral. Sin embargo, los **biomarcadores moleculares** asociados a estos procesos no están completamente caracterizados, lo que limita el desarrollo de nuevas estategias terapéuticas.  

Estudiar los perfiles de expresión génica asociados con la invasión y migración celular es crucial para:

* Identificar **biomarcadores** que permitan pronósticos más precisos.
* Comprender los mecanismos moleculares de la metástasis.

Este proyecto busca **replicar y validar** los hallazgos del artículo mediante un pipeline bioinformático de análisis de RNA-seq y análisis de expresión diferencial.

---

## Calendario de trabajo

| Actividad                                                                | Fecha          | Responsable    | Entregable                     |
| ------------------------------------------------------------------------ | -------------- | -------------- | ------------------------------ |
| Descripción del proyecto                                                 | 28 ago – 2 sep | Todo el equipo | Documento markdown             |
| Descarga y control de calidad de datos                                   | 5 – 20 sep     | Andrea         | Reporte QC (FastQC)            |
| Alineamiento y cuantificación de lecturas                                | 21 – 30 sep    | Antonio        | Archivos, tabla de conteos     |
| Análisis de expresión diferencial (DESeq2)                               | 1 – 15 oct     | Alejandra      | Script, tablas de DEGs       |
| Anotación funcional                                                      | 16 – 25 oct    | Alejandra      | Reporte en markdown            |
| Integración de resultados y gráficas                                     | 26 oct – 5 nov | Antonio        | Figuras en repositorio         |
| Redacción del reporte de resultados                                      | 6 – 15 nov     | Andrea         | Documento markdown             |
| Presentación del proyecto                                                | Diciembre      | Todo el equipo | Repositorio GitHub (release)   |

---

## Metodología

**Preguntas de investigación**

1. ¿Qué genes están diferencialmente expresados entre células osteoblásticas normales y líneas celulares de osteosarcoma metastásicas?

**Pipeline propuesto**

1. **Localización y descarga de datos** (ENA, proyectos PRJEB65286, PRJNA698672, PRJNA51801).
2. **Control de calidad de lecturas FASTQ** (FastQC, Trim Galore).
3. **Alineamiento**  (STAR).
4. **Generación de tabla de conteos**
5. **Análisis de expresión diferencial** (DESeq2)
6. **Anotación funcional** 
7. **Visualización de resultados**
8. **Integración y discusión de resultados** 

## Resultados esperados

* Identificación de genes diferencialmente expresados entre células normales y osteosarcoma metastásico.
* Detección de biomarcadores relevantes en la metástasis. 
* Reproducción de resultados clave reportados en el artículo, con figuras. 
* Generación de un repositorio reproducible en GitHub. 

