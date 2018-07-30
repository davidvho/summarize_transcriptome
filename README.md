# summarize_transcriptome
Python script to get useful statistics from a transcriptome and its annotation

Inputs: TRINITY transcripome assembly (.fasta) and TRINOTATE annotation file

Note that the Trinotate annotation file is given to the user as a .xls file from the software, please convert to a tab-delimited file before proceeding.

The summary statistics you get from the Trinity assembly are as follows:

No. of genes
No. of transcripts
No. of total assembled bases
Average transcript length (bp)
Min gene length (bp)
Max gene length (bp)
Number of genes > 1 Kb
Number of genes > 5 Kb
Number of genes > 10 Kb
Transcript N50 (bp)
GC content (%)

From the Trinotate file input, you will get output files with values you can use for plotting and GO terms for each gene.
