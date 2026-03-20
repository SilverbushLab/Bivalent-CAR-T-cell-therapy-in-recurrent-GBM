Welcome to the Silverbush Lab Github page containing code for data exploration and figure reproduction from our manuscript "The critical role of the host endogenous immune compartment after intracerebroventricular CAR T cell therapy in recurrent GBM".

Data availability

	Raw data are protected to preserve patient privacy. Processed data can be downloaded from GEO with accession code GSE296419. Processed data are available in the form of counts  (processed for ambient RNA correction and empty droplet removal with Cellbender 0.3.0, .h5 file), contigs containing TCR clonotype data (csv.gz), per-cell annotation files, and fully annotated analysis-ready Seurat objects (.csv.gz). Counts and TCR clonotype data are available individually for each sample. Seurat objects are available as Supplementary .rds files and are labelled according to sample type. Note that the Infusion Product sample for Patient 14 is in Seurat object (.rds). Bulk RNA sequencing data are also in supplementary files: as raw counts (GSE296419_FFPE_bulkRNA_raw_counts.tsv.gz) and TPM (GSE296419_FFPE_bulkRNA_tpm_gene_symbol_id_dedup.tsv.gz). 

	GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9543458

	Seurat objects:
	CSF Seurat object: GSE296419_CSF_seurat_final.rds
	Tumor Seurat object: GSE296419_srt_tumor_final.rds
	IP Seurat object (P1-P8): GSE296419_srt_IP_disc.rds
	IP Seurat object(P14): GSE296419_srt_P14_IP.rds


Code availability:
	Code is provided to reproduce figures from our manuscript. Code is separated by relevant tissue type into CSF, Infusion product, and tumor folders. Preprocessing code is also available to generate Seurat objects from counts. Please find the associated preprocessing code for each tissue type in their respective folders. 


Docker:
Docker images are available which generate an environment suitable for running all code for generating figures. There are two separate Docker images: One for the CSF/Infusion Product and one for tumor samples. 

	CSF and Infusion Product:
	Can be accessed at the Silverbush Lab DockerHub account (https://hub.docker.com/repository/docker/silverbushdana/jupyter_seurat_infercnv_v2/general) or use: 

	docker pull silverbushdana/jupyter_seurat_infercnv_v2

	Tumor:
	Can be accessed at the Silverbush Lab DockerHub account (https://hub.docker.com/repository/docker/silverbushdana/jupyter_seurat_infercnv_v2/general) or use: 

	docker pull silverbushdana/jupyter_seurat_infercnv_v2

	This docker image can also be generated from scratch using the dockerfile contained in the "tumor_analysis" folder


Resources:
	Many scripts pull from files in the resources folder. There should be no need to edit files in resources. When loading resources, ensure your path is correct.


