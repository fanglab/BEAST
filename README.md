# BEAST
Bacterial Epigenomics Analysis SuiTe

## Description

With only a few thousand bacterial methylomes published to date, it is becoming increasingly evident that epigenetic regulation of gene expression is highly prevalent across bacterial species. Despite the exciting prospects for studying epigenetic regulation, our ability to comprehensively analyze bacterial epigenomes is limited by a bottleneck in integratively characterizing methylation events, methylation motifs, transcriptomic data, and functional genomics data. In this regard, we provided the first comprehensive comparative analysis of a large collection of epigenomes in a single bacterial species, as well as a detailed roadmap that can be used by the scientific community to leverage the current status quo of epigenetic analyses [(Oliveira et al., 2019, Nat. Microbiology)](https://www.nature.com/articles/s41564-019-0613-4). 

PLEASE NOTE: BEAST includes a set of R and wrapper shell scripts that allow to reconstitute the major analyses steps of this publication. However, this pipeline requires multiple dependencies, which have to be installed prior the use of BEAST. The user has two possibilities:

* For those having a basic programming knowledge, they are ore more than welcome to install these dependencies and run BEAST (a detailed step-by-step tutorial is nevertheless provided in the [Documentation](#documentation)). 
* For newcomers or those less familiar with installing some of the dependencies, we have provided links to Docker containers with pre-installed dependencies [here](https://beast-docs.readthedocs.io/en/latest/usage.html#running-beast-via-docker).

## Sections
Details on the requirements and usage of all shell wrapper scripts can be accessed by **scriptname -h (or -help)**.


### 1) Motif refining tool for the SMRT-seq pipeline
* **Motif_Refine.R** performs the refining of methylation motifs estimated by the SMRT portal pipeline. Takes a Parameter file as input (provided as example). Each genome to be analysed should be placed in a separate folder containing the genome fasta file, and the output files of the SMRT pipeline (Modifications and Motif files).

### 2) Analysis of methylation motif enrichment / depletion using Markov models and a multiscale representation framework
* **GO_Abundance.sh** maps methylation motifs in a FASTA file and computes their [scores of exceptionality](https://www.worldscientific.com/doi/abs/10.1142/9789814327732_0002) using Markov models.
* **GO_MSR.sh** highlights chromosomal regions with enrichment / depletion of a given signal (methylation motifs in this case) using the [Multiscale Signal Representation (MSR)](https://www.ncbi.nlm.nih.gov/pubmed/24727652) method. An example parameters.txt is also given. 

 * **Plot_MSR.R** is used to plot the pruned MSR output.

### 3) Analysis of orthologous conserved / variable methylation motifs
* **GO_ConsVar.sh** performs multiple whole-genome alignment and looks at the conservation of methylation motifs across genomes.

### 4) Analysis of Transcription Factor Binding Sites (TFBS) and Transcription Start Sites (TSS) in bacterial genomes
* **GO_TFBS.sh** takes a TFBS multifasta file, and computes the corresponding PSSM and TFBS hits in a given genome. An example of a TFBS multifasta file (XylR.fasta) is provided.
 * **GO_PSSM.R** is used to build a PSSM.

* **GO_TSS.sh** estimates transcription start sites (TSS) through the [reconstruction of a transcriptional landscape](https://www.ncbi.nlm.nih.gov/pubmed/24470570) from RNA-Seq data. An example *.chrom.sizes file in IGV format is provided.


### 5) Pipeline to compute differential expression genes from RNA-seq data
* **GO_GetCounts.sh** performs RNA-seq paired-end read cleaning, mapping, and counting for differential expression (DE) analysis.
 * **GO_DE.R** performs the DE analysis. A CDIF.count.txt and a colData.csv files are provided as examples.


### 6) Pipeline to compute homologous recombination and horizontal gene transfer
* **GO_HR.sh** computes [homologous recombination](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041) events given an ordered core genome alignment and corresponding phylogenetic tree in Newick format.
 * **GO_TTR-MBL.R** is used to compute transition / transversion ratios (TTR) and mean branch lengths (MBL).

* **GO_HGT.sh** uses a pan-genome matrix and phylogenetic tree to perform [ancestral reconstruction](https://www.ncbi.nlm.nih.gov/pubmed/20551134) and infer family and lineage specific characteristic along the species' tree. An example pan-genome matrix and phylogenetic tree are provided.
 
 
### 7) Detection of mobile genetic elements and defense systems
* Restriction modification systems were detected using dedicated HMM profiles and scripts previously published in [Oliveira et al, 2014](https://www.ncbi.nlm.nih.gov/pubmed/25120263) and available [here](https://github.com/pedrocas81).
* **GO_CRISPRs.sh** computes [CRISPRs](https://www.ncbi.nlm.nih.gov/pubmed/17577412) from a FASTA file and parses the output into a TAB format.
* **GO_Prophages_Integrons.sh** computes prophages using [Phage Finder](https://academic.oup.com/nar/article/34/20/5839/3100473) and integrons using [Integron Finder](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4889954/).


## Documentation
For a comprehensive guide on each section and corresponding scripts, please see the documentation available [here](https://beast-docs.readthedocs.io/en/latest/).
