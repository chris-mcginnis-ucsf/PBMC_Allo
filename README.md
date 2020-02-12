# Analysis and visualization code for McGinnis et al., "No detectable alloreactive transcriptional responses during donor-multiplexed single-cell RNA sequencing of peripheral blood mononuclear cells." 

![alternativetext](schematic.png)

PreProcessing.R : (i) MULTI-seq and Cell Hashing FASTQ pre-processing, (ii) Gene expression data pre-processing using Seurat, (iii) Gene expression data QC, (iv) DoubletFinder, (v) Cell type annotation, (vi) Figure S1 visualizations.

SampleClassification.R : (i) deMULTIplex classification for MULTI-seq and Cell Hashing count matrices, (ii) demuxEM classification for MULTI-seq and Cell Hashing count matrices, (iii) souporcell in silico genotyping, (iv) Figure 2 visualizations.

TrimaAnalysis.R : (i) Sex-corrected marker analysis between Trima and Ficoll PBMCs, (ii) Figure S3 visualizations.

Alloreactivity.R : (i) CD4+ T-cell and classical monocyte sub-clustering, (ii) Jensen-Shannon Divergence computation, (iii) hierarchical clustering of JSD results, (iv) Label permutation test, (v) Figure S2 and Figure 3 visualizations.

 
