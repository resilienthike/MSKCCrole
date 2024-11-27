# OVERVIEW:

This folder contains a project using an edited version of the scGEAToolbox (Cai (2020) Bioinformatics) to analyze single-cell RNA-sequencing data.

# Dataset Description:
The dataset used in this homework comes from the single-cell RNA-sequencing of glioblastoma tumors described in Neftel et. al. 
An integrative Model of Cellular States, Plasticity and Genetics for Glioblastoma. (2019) Cell, DOI: 10.1016/j.cell.2019.06.024 and was downloaded from the NCBI Gene Expression Omnibus (GEO) repository under the accession number GSE131928. 
Each sample in the dataset refers to the single-cell transcriptomes of cells from the tumor of one
patient. Tumors were dissociated, cells sorted into 2 fractions using the pan-immune marker
CD45 and subjected to single-cell RNA-sequencing using the SMART-seq2 protocol. For a given
sample, the majority of cells are malignant, and a subset are normal, tumor-associated cells.
Transcripts per million were converted to inferred mRNA counts per cell using Census (Qiu et al.
(2017) Nature Methods) as implemented in the Monocle2 function relative2abs().
