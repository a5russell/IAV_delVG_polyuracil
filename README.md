# IAV delVG polyuracil
The code contained within this repository describe the analysis of data, and generation of figures, in our attempts to understand positions in delVGs that alter the innate immune response, with a focus on the polyuracil tract.

Analysis are split across several different jupyter notebooks, with descriptions below, and figures inline.

All intermediate files less than 10mb were pushed to this repository to aid in understanding our analysis.
Larger files, such as sequencing, will need to be downloaded and reconstructed in the apppropriate folders in order to re-run this analysis.
Specifically, place NGS files within appropriate directories within Sequencing to rerun all sequencing analyses.

Data may be found in GEO (including some processed files) at accession GSE276697, GSE350335, and GSE350336.

## Directories

The following directories, and their purpose, exist within this repository.

- <b>Database</b>       Repository for influenza genomic sequences and the same with adapters for appropriate mapping as well as sequence with plasmid backbone. A/WSN/1933 BLAST database and STAR indicies provided. Additionally contains PB1 sequences used for analysis of natural diversity.
-  <b>Figures</b>       Figures, main and some supplemental.
- <b>SuppFigures</b>    Additional supplemental figures
- <b>Results</b>        Final datafiles for analyses after processing. Most critical ones additionally on GEO.
- <b>Scripts</b>        Short scripts written for this analysis. Seperated from jupyter notebooks for readability and portability.
- <b>Sequencing</b>     Folder that would contain NGS samples. Folder achitecture essential to run off pipeline.
  

## Dependencies

Code within this repository was run with the following tools, and versions, installed and available from PATH. Specific websites and documentation are provided where available. 

- <b>Python</b>      run with version 3.12.2 . Available from https://www.python.org/downloads/release/python-370/
- <b>Trimmomatic</b> run with version 0.39. Available from http://www.usadellab.org/cms/?page=trimmomatic
- <b>STAR</b>        run with version 2.7.11b. Available from https://github.com/alexdobin/STAR
- <b>Samtools</b>    run with version 1.19. Available from http://www.htslib.org/
- <b>FastQC</b>      run with version 0.12.1. Available from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- <b>CellRanger</b> run with version 1.3.1. Available from 10x genomics.
  

The following python packages and versions were used. All were installed using conda. (https://docs.conda.io/en/latest/)
- <b>numpy</b>       run with version 1.26.4. (https://numpy.org/)
- <b>matplotlib</b>  run with version 3.8.3. (https://matplotlib.org/)
- <b>seaborn</b>     run with version 0.13.2. (https://seaborn.pydata.org/)
- <b>pandas</b>      run with version 2.2.1.(https://pandas.pydata.org/)
- <b>scipy</b>       run with version 1.12.0. (https://www.scipy.org/)
- <b>statsmodels</b> run with version 0.14.1. (https://www.statsmodels.org/stable/index.html)
- <b>PyDESeq2</b>    run with version 0.5.4. (https://pydeseq2.readthedocs.io/en/stable/index.html)

## Jupyter notebooks

General descriptions of pipelines within each notebook described below.

### Assigment_alignment_DMS_libraries.ipynb

Aligning DMS libraries to appropriate templates and base calling.

### DMS_analysis.ipynb

Visualization and statistical analysis of DMS libraries.

### Natural_analysis.ipynb

Alignment, visualization, and analysis of resequencing from Mendes and Russell 2021 at 5' ends of vRNA of PB2, PB1, and PA.

### total_RNAseq.ipynb

Alignment, analysis, and visualization of total RNAseq of PB1_445_350 U16A, U20A, and U20 with baloxavir conditions.

### scRNAseq.ipynb

Analysis and visualization (downstream of cellranger count) of scRNAseq data from WT and PB1_177_385 infections in A549 cells.

### qPCR_titers.ipynb

Analysis and visualization of qPCR and titer data presented in the manuscript.

