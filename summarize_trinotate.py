##################################
## ARGPARSE ######################
##################################

import argparse
import numpy as np
import collections

def mm():
	parser = argparse.ArgumentParser(description='Given a transcriptome assembled from TRINITY (.fasta) and an annotation file from TRINOTATE, this script will output summary statistics for the transcriptome assembly and produce 4 files from the Trinotate annotation that can be used for plotting (like in R): (1) Isoform count distribution, (2) Organism BLAST hits, (3) First GO term for each gene, (4) All GO terms. Important note: Trinotate output needs to be converted into a tab-delimited file before use.')
	parser.add_argument("-t1", "--transcriptome", help="Required. FASTA format.", required=True, type=str)
	parser.add_argument("-t2", "--trinotate", help="Required. Be sure to convert the .xls from Trinotate into a tab-delimited file.", required=True, type=str)
	return parser.parse_args()

args = mm()

# File
trinity_file = args.transcriptome
trinotate_file = args.trinotate

##################################
### FUNCTIONS ##################################
##################################

def GCcontent(sequence):
    GCcount = 0
    for letter in sequence:
        if letter == "G" or letter == "C":
            GCcount += 1
    return GCcount

def trinity_stats(file):
    '''Output statistics about genome-guided Trinity transcriptome assembly'''
    with open(file) as trinity:
        final = {}

        for line in trinity:
            line = line.strip()
            if line.startswith(">"):
                ID = line
                final[ID] = " "
            else:
                final[ID] += line
                
        trans_no = 0
        uniq_gene = {}
        base_no = 0
        length_trans = {}
        GC = 0
        for key, val in final.items():
            length_trans[key] = len(val)
            trans_no += 1
            header = key.split("_")
            gene = header[2]+"_"+header[3]
            if gene in uniq_gene:
                    uniq_gene[gene] += 1
            else:
                uniq_gene[gene] = 1
            
            base_no += len(val)
            
            GC += GCcontent(val)
        
        lengths = list(length_trans.values())
        avg_trans = np.mean(lengths)
        array_lengths = np.array(lengths)
        
        print("No. of genes |", len(uniq_gene))
        print("No. of transcripts |", trans_no)
        print("No. of total assembled bases |", base_no)
        print("Average transcript length (bp) |", round(avg_trans,2))
        print("Min gene length (bp) |", round(np.min(lengths), 2))        
        print("Max gene length (bp) |", round(np.max(lengths), 2))
        print("Number of genes > 1 Kb |", sum(array_lengths > 1000))
        print("Number of genes > 5 Kb |", sum(array_lengths > 5000))
        print("Number of genes > 10 Kb |", sum(array_lengths > 10000))
        print("Transcript N50 (bp) |", np.median(lengths))
        print("GC content (%) |", round(GC/base_no*100, 2))

def isoform_no(file):
    '''From a genome-guided transcriptome assembly, outputs the frequency of isoforms per gene'''
    
    with open(file) as trinity:
        final = {}
        for line in trinity:
            line = line.strip()
            if line.startswith(">"):
                ID = line
                final[ID] = " "
            else:
                final[ID] += line
                
        uniq_gene = {}
        for key, val in final.items():
            header = key.split("_")
            gene = header[2]+"_"+header[3]
            if gene in uniq_gene:
                    uniq_gene[gene] += 1
            else:
                uniq_gene[gene] = 1
    
    return dict(Counter(uniq_gene.values()))

def BLAST_organisms(file):
	with open(file) as fh:
		first_line = fh.readline()
		species = {}
		genes = {}
		for line in fh:
			items = line.strip().split("\t")
			gene = items[0]
			blast_hit = items[2]
			if gene not in genes:                         # just do it onces for a gene
				genes[gene] = "done"
				if not blast_hit.startswith("."):
					taxonomy = blast_hit.split(";")
					genus = taxonomy[-1]                  # last item in a list
					if genus not in species:
						species[genus] = 1
					else:
						species[genus] += 1
	with open("blast_species.txt", "w") as blast:
		blast.write("genus"+"\t"+"count"+"\n")
		for k,v in species.items():
			blast.write(k+"\t"+str(v)+"\n")


##################################
### MAIN ##################################
##################################

### Trinity-transcriptome stats

trinity_stats(trinity_file)

#### ISOFORM file

iso_no = isoform_no(trinity_file)

with open("isoform_count.txt", "w") as fh:
    fh.write("Isoform_count"+"\t"+"Gene_no"+"\n")
    for k,v in iso_no.items():
        fh.write(str(k)+"\t"+str(v)+"\n")

#### BLAST file
BLAST_organisms(trinotate_file)

#### TXT file with the first GO term of each gene

with open(trinotate_file) as fh:
    g_g = {}
    for line in fh:
        gene_go = line.strip().split("\t")
        GO = gene_go[12].split("^")[0]
            
        if gene_go[0] not in g_g:
            if "." not in GO:
                g_g[gene_go[0]] = GO

with open("first_GO_term_of_gene.txt", "w") as file:
    for v,k in g_g.items():
        file.write(k+"\n")

#### TXT file with all GO terms for genes

with open(trinotate_file) as fh, open("all_GO_terms_all_genes.txt", "w") as write:
    first_line = fh.readline()
    go_num = 0
    for line in fh:
        gene_go = line.strip().split("\t")
        GO_term = gene_go[12].split("`")
        if "." not in GO_term:
            go_num += len(GO_term)
            for item in GO_term:
                what = item.split("^")
                write.write(item[:10]+"\t"+what[1]+"\t"+what[2]+"\n")
                
print("Files written: (1) Isoform count distribution, (2) Organism BLAST hits, (3) First GO term for each gene, (4) All GO terms")