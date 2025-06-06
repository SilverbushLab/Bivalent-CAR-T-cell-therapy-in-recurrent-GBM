CSF and infusion product analysis:

Fastq files can be obtained from the GEO repository with the accession code GSE296419. This repository contains raw fastq files, Cellbender-processed count matrices, and cell type annotations. 

A Docker image is available which generates an environment suitable for running all code for generating CSF and infusion-product related figures. This can be accessed at the Silverbush Lab DockerHub account (https://hub.docker.com/repository/docker/silverbushdana/jupyter_seurat_infercnv_v2/general) or use: 

docker pull silverbushdana/jupyter_seurat_infercnv_v2

There are separate preprocessing (pp) scripts for CSF and infusion product samples. Raw data (fastq) and processed data (Cellbender-processed .tsv files) can be obtained from Gene Expression Omnibus at accession code GSE296419 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296419). These data are private until publication of the associated manuscript, but can be accessed using the token provided alongside manuscript submission.

Preprocessing files for CSF and infusion product samples should be run in the following order: Load, Scrublet, Clustering, TCR integration. These files output a Seurat object (.rds file) that can be used as input to figure-generating scripts.




Tumor analysis: 

Code for tumor cell analysis is run separately from the CSF and Infusion Product analysis.

Use the included dockerfile in "Dockerfile_tumor_analysis" create an environment containing required packages. This includes both R and python packages. Alternatively, the corresponding image can be obtained from the Silverbush Lab DockerHub page (https://hub.docker.com/repository/docker/silverbushdana/e20_with_rstudio_jupyter_updated/general) using:

docker pull silverbushdana/e20_with_rstudio_jupyter_updated

Fastq files can be obtained from the GEO repository with the accession code GSE296419. This repository contains raw fastq files, Cellbender-processed count matrices, and cell type annotations. Either fastq files or cellbender-preprocessed files can be used as input into preprocessing script 1. Cellbender-preprocessed files may be preferred, as Cellbender is highly computationally intensive. 

Raw or preprocessed files will need to be moved into an directory for import into R or python.

Inside the container, use the 4 tumor_preprocessing files to import and preprocess fastq files and output a useable Seurat object as a .rds file. This .rds file is used as input to the figure-generating .rmd files. These files are numbered according to the correct running order.


Many scripts pull from files in the resources folder. There should be no need to edit files in resources.