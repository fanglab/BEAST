=====
Usage
=====

Section 1
=============

------------------------------------------------
Motif refining
------------------------------------------------

* A separate folder should be created for each genome to be analyzed. Each folder should contain 3 files: the genome FASTA file (genome.fasta), and the output files of the SMRT pipeline (modifications.csv and motif_summary.csv files). An output folder should also be created.
* Manually add the paths to each of the above folders in Parameter_File.txt. Mean methylation motif fraction is already given for a sample of 50 common motifs. This motif list can be tailored as desired.

Fundamental dependencies
------------------------
* Recent version of `R <http://cran.fhcrc.org>`_

Example usage
--------------
.. code-block:: console

	$ Rscript Motif_Refine.R Parameter_File.txt

This program will take a genome fasta file and SMRT-seq kinetic data (modifications.csv and motif_summary.csv) 
and output a list of refined (trimmed) methylated motifs.

One data directory will be created for each genome folder. This folder will have output files for each motif indicated in the Parameter_File.txt. 
The latter include one exp.csv and one con.csv for motif statistics respectively for the leading and lagging strands, graphical distribution 
of scores and IPD ratios for each base of the motif, and expdat2.csv with kinetic information for all possible motif combinations. 


Section 2
=============

------------------------------------------------
Motif count exceptionalities using Markov models
------------------------------------------------

Fundamental dependencies
------------------------
* `R'MES <http://migale.jouy.inra.fr/?q=fr/rmes>`_
* `SeqKit <https://github.com/shenwei356/seqkit>`_

Example usage
-------------
.. code-block:: console

	$ ./GO_Abundance.sh <*.fasta> <motif> <word_length> <markov_model_order> <approximation method>
	
This wrapper script takes a genome FASTA file and computes the statistical over/under-representation of a given motif under Markov models of order N.

<word_length>: length l of the word to search

<markov_model_order>: From 1 to a maximum of l-2

<approximation_method>: gauss (gaussian), compoundpoisson (Poisson)
	
For example:
./GO_Abundance.sh genome.fasta CAAAAA 6 4 compoundpoisson

------------------------------------------------
Multi-Scale Representation (MSR) of methylation motifs
------------------------------------------------

Fundamental dependencies
------------------------
* `SeqKit <https://github.com/shenwei356/seqkit>`_
* MSR standalone executable (msr_runtime_SIGNAL and wrapper run_msr_runtime_SIGNAL.sh). Can be found `here <https://github.com/tknijnen/MSR>`_.
* MATLAB Compiler Runtime (MCR) version 8.1 (R2013a). Can be found `here <https://mathworks.com/products/compiler/matlab-runtime.html>`_.

Example usage
--------------
.. code-block:: console

	$ ./GO_MSR.sh <*.fasta> <motif> <mcr_directory> <parameter_file> <MSR_output_filename (same as in parameter file)>

This wrapper script takes (1) a genome FASTA file and (2) a methylation motif, builds a signal file, runs the MSR pipeline, parses, and plots its output.


Section 3
=============

------------------------------------------------
Conservation of methylation motifs across genomes
------------------------------------------------

Fundamental dependencies
------------------------
* `ProgressiveMauve <http://darlinglab.org/mauve/mauve.html>`_
* `SeqKit <https://github.com/shenwei356/seqkit>`_
* `convertAlignment.pl <https://github.com/lskatz/lskScripts>`_
* `Bedtools <https://bedtools.readthedocs.io/en/latest/>`_
* `jvarkit <https://github.com/lindenb/jvarkit>`_
* `VCFtools <https://vcftools.github.io/index.html>`_


Example usage
--------------
.. code-block:: console

	$ ./GO_ConsVar.sh <minimal length of LCB> <number of genomes to align> <species_prefix> <MAUVE_DIR> <motif>

This wrapper performs multiple whole-genome alignment and computes orthologous (conserved/variable) and non-orthologous methylation motifs across genomes. The list of genome FASTA files should be placed in the same folder as GO_ConsVar.sh.

To align 10 genomes of C. difficile with a minimum length of local collinear blocks (LCB) of 50 bp and compute conserved CAAAAA motifs: 

./GO_ConsVar.sh 50 10 CDIF /path/to/mauve/ CAAAAA

It outputs 4 main files: 

PREFIX.Indels.txt: All orthologous variable motif positions harbouring indels.

PREFIX.SNPs.txt: All orthologous variable motif positions harbouring SNPs.

PREFIX.Conserved.txt: All orthologous conserved motif positions. 

PREFIX.NonOrthologous.txt: All non-orthologous positions. Usually are found within MGEs.

Section 4
=============

------------------------------------------------
TFBS mapping
------------------------------------------------

Fundamental dependencies
------------------------
* Recent version of `R <http://cran.fhcrc.org>`_
* Bioconductor
* MAST `(MEME Suite) <http://meme-suite.org/doc/download.html>`_


Example usage
--------------

Install the latest release of R, then get the latest version of Bioconductor by starting R and entering the commands:

.. code-block:: console

	$ if (!requireNamespace("BiocManager", quietly = TRUE))
	     install.packages("BiocManager")
	  BiocManager::install()
	
Install the Biostrings package:

.. code-block:: console

	$ BiocManager::install("Biostrings")

Then run:

.. code-block:: console

	$ ./GO_TFBS.sh <TFBS_multifasta> <*.fasta> <TF_name>

This wrapper takes a TFBS multifasta, computes a PSSM, and corresponding TFBS hits in a given FASTA file. To compute a hit list of XylR TFBSs in genome.fasta using the multifasta XylR.fasta (provided):

./GO_TFBS XylR.fasta genome.fasta XylR

The output file XylR_TFBS.txt will contain all raw TFBS hits (no thresholds introduced).

------------------------------------------------
TSS mapping
------------------------------------------------

Fundamental dependencies
------------------------
* `Samtools <http://samtools.sourceforge.net>`_
* `IGVtools <https://software.broadinstitute.org/software/igv/igvtools>`_
* `Parseq <http://www.lgm.upmc.fr/parseq/>`_
* `GSL <https://www.gnu.org/software/gsl/>`_

Example usage
--------------

.. code-block:: console

	$ ./GO_TSS.sh <bam_file> <chrom_size_file> <*.fas> <counts_folder_path> <results_folder_path> <parseq_folder_path>


This wrapper runs the Parseq program for reconstruction of the transcriptional landscape from RNA-Seq data, and infers TSSs from abrupt shifts in transcription levels. It takes as input a bam file, a chromosome fasta file, and a chromosome size file (provided).

Section 5
=============

------------------------------------------------
Differential gene expression analysis
------------------------------------------------

Fundamental dependencies
------------------------
* `Java <https://www.java.com/download/>`_
* `BWA <http://bio-bwa.sourceforge.net/>`_
* `AdapterRemoval <https://github.com/MikkelSchubert/adapterremoval>`_
* `Trimmomatic <https://github.com/timflutre/trimmomatic>`_
* `SortMeRNA <https://github.com/biocore/sortmerna>`_
* `Samtools <http://samtools.sourceforge.net>`_
* `featureCounts <http://subread.sourceforge.net>`_


Example usage
--------------

.. code-block:: console

	$ ./GO_GetCounts.sh <*.fastq1 file> <*.fastq2 file> <*.fasta reference file> <*.SAF annotation file> <PATH to SILVA rRNA_databases folder> <PATH to SILVA index folder> <PATH to Trimmomatic.jar> <PATH to adapters *.fa> <file prefix>

This wrapper performs RNA-seq paired-end read cleaning and mapping for differential expression analysis.
It takes as input the FASTQ files, a reference FASTA file for read mapping, and a SAF annotation file.

The SAF annotation format (example below) has five required columns, including GeneID, Chr, Start, End and Strand. These columns can be in any order. More columns can be included in the annotation. Columns are tab-delimited. Column names are case insensitive. GeneID column may contain integers or character strings. Chromosomal names included in the Chr column must match those used included in the mapping results, otherwise reads will fail to be assigned. 

.. code-block:: console

	GeneID	Chr	Start	End	Strand
	Gene_A	chr1	134	1376	+
	Gene_B	chr1	4031	4528	+
	Gene_C	chr1	4909	5313	-
	Gene_D	chr1	9034	9267	-
	...


The adapters_*.fa file is a multifasta file with adapter sequences identified in the FASTQ files, for example, via the AdapterRemoval tool:

.. code-block:: console

	$ AdapterRemoval --identify-adapters --file1 <*.fastq1 file> --file2 <*.fastq2 file>

DE analysis is performed by DESeq2 through the dependent GO_DE.R script. This produces the following output files:

outinfo.txt: outputs for each gene, information on non-normalised and normalised read counts, log2FC, P value, and FDR.

MyData.csv: Subset of outinfo.txt containing log2FC, standard error, P value, and FDR.

plotMA.pdf and plotMAShrunk.pdf: Regular and shrunken MA plots (log2FC vs mean of normalized read counts).

PCA.pdf and dendrogram.pdf: Builds a dendrogram to evaluate sample clustering and a 2D PCA plot to check for potential outliers.

Section 6
=============

------------------------------------------------
Gene flux analysis - Horizontal Gene Transfer (HGT)
------------------------------------------------

Fundamental dependencies
------------------------
* `Java <https://www.java.com/download/>`_
* `Count <http://www.iro.umontreal.ca/~csuros/gene_content/count.html>`_

Example usage
--------------

.. code-block:: console

	$ ./GO_HGT.sh <Pan_genome_matrix.csv> <Newick_tree> <Species_prefix> <Posterior_gain_probability>

This wrapper runs Count to perform ancestral reconstruction and infer family and lineage specific characteristics along the evolutionary tree. It takes as input a pan-genome matrix file (example provided for 45 C. difficile genomes), which can be obtained, for example, with `Roary <https://sanger-pathogens.github.io/Roary/>`_. It also requires the species tree in Newick format (example provided). The user is also required to specify the posterior gain probability for the family sizes at inner nodes.

The final output file *.Gains.out contains the sum of gene families acquired at each tip of the phylogenetic tree.

------------------------------------------------
Gene flux analysis - Homologous Recombination (HR)
------------------------------------------------

Fundamental dependencies
------------------------
* Recent version of `R <http://cran.fhcrc.org>`_
* `ClonalFrameML <https://github.com/xavierdidelot/ClonalFrameML>`_

Example usage
--------------

.. code-block:: console

	$ ./GO_HR.sh <ordered_core_alignment> <Newick_tree> <Species_prefix>

This wrapper runs ClonalFrameML given an ordered core genome alignment and corresponding phylogenetic tree in Newick format. The ordered core genome alignment can be extracted from a progressiveMauve alignment using stripSubsetLCBs as described `here <https://github.com/xavierdidelot/ClonalOrigin/wiki/Usage>`_.

The output files include .smout.log.txt for the standard model run, .pbmout.log.txt for the per-branch model run (recombination parameters are estimated per branch), a list of core sites (Core_Sites.txt), and a pdf with a graphical representation of the latter.

Section 7
=============

------------------------------------------------
CRISPR detection
------------------------------------------------

Fundamental dependencies
------------------------
* `Java <https://www.java.com/download/>`_
* `CRT <https://github.com/xavierdidelot/ClonalFrameML>`_

Example usage
--------------

.. code-block:: console

	$ ./GO_CRISPRs.sh <*.fasta> <CRT_filename.jar> <minNR> <minRL> <maxRL> <minSL> <maxSL> <searchWL>
	
	minNR: minimum number of repeats a CRISPR must contain; default 3
	minRL: minimum length of a CRISPR's repeated region; default 19
	maxRL: maximum length of a CRISPR's repeated region; default 38
	minSL: minimum length of a CRISPR's non-repeated region (or spacer region); default 19
	maxSL: maximum length of a CRISPR's non-repeated region (or spacer region); default 48
	searchWL: length of search window used to discover CRISPRs; (range: 6-9)

This wrapper runs CRT on a FASTA file and parses the output file (.crispr_raw) into a tab delimited output (.crispr_parsed).

------------------------------------------------
Prophage detection
------------------------------------------------

Fundamental dependencies
------------------------
* `Phage Finder <http://phage-finder.sourceforge.net>`_
* `EDirect <https://dataguide.nlm.nih.gov/edirect/install.html>`_

Example usage
--------------

.. code-block:: console

	$ ./GO_Prophages.sh <genome_accession_number> <prefix>


This wrapper searches prophages in full genome sequences using Phage Finder.


------------------------------------------------
Integron detection
------------------------------------------------

Fundamental dependencies
------------------------
* `Integron Finder <https://github.com/gem-pasteur/Integron_Finder>`_

Example usage
--------------

.. code-block:: console

	$ ./GO_Integrons.sh <*.fasta>

This wrapper computes integrons from a FASTA file using Integron Finder under default conditions. A list of optional arguments are detailed here, and can be added to the script as 

By default, integron_finder will output 3 files under Results_Integron_Finder_mysequences:

mysequences.integrons : A file with all integrons and their elements detected in all sequences in the input file.
mysequences.summary : A summary file with the number and type of integrons per sequence.
integron_finder.out : A copy of standard output.

