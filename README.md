# summarize_transcriptome
Python script to get useful statistics from a transcriptome and its annotation

Inputs: TRINITY transcripome assembly (.fasta) and TRINOTATE annotation file

Links to Trinity and Trinotate:

https://github.com/trinityrnaseq/trinityrnaseq/wiki
https://github.com/Trinotate/Trinotate.github.io/wiki

Note that the Trinotate annotation file is given to the user as a .xls file from the software, please convert to a tab-delimited file before proceeding.

The summary statistics you get from the Trinity assembly are as follows:

No. of genes<br>
No. of transcripts<br>
No. of total assembled bases<br>
Average transcript length (bp)<br>
Min gene length (bp)<br>
Max gene length (bp)<br>
Number of genes > 1 Kb<br>
Number of genes > 5 Kb<br>
Number of genes > 10 Kb<br>
Transcript N50 (bp)<br>
GC content (%)<br>

From the Trinotate file input, you will get output files with values you can use for plotting and GO terms for each gene.

Last updated: 31.07.2018
