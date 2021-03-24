import os
import pandas as pd
import sys


"""USAGE:  python Generate_single_cell_matrix.py [base_folder] [fastq_1] [fastq_2]

#Files you need in base folder:

#1) t2g.py: downloaded from https://github.com/pachterlab/kallisto-transcriptome-indices/releases
#2) Barcodes (10xv3_whitelist.txt): downloaded from https://github.com/BUStools/getting_started/releases
#3) transcripts_to_genes.txt - generated in preprocessing
#4) Transcriptome index (.index) -- generated in preprocessing

Results will create an output directory called 'bus_output'

"""


"""Assign variables"""

base_folder = sys.argv[1]
fastq_1 = sys.argv[2]
fastq_2 = sys.argv[3]


"""Navigate to base folder and create output folders"""

os.chdir(base_folder)

sample = fastq_1.split("/")[-1].split(".fastq")[0]
    
if not os.path.isdir("{0}/bus_output".format(base_folder)):
    os.system("mkdir bus_output")

if not os.path.isdir("{0}/bus_output/{1}".format(base_folder,sample)):
    os.system("mkdir bus_output/{0}".format(sample))
    

"""Do preprocessing, if necessary"""
preprocessing = True

if preprocessing:
    
    os.system("kallisto index -i hg38.gencode_v24.basic_ccds_nopar.transcripts.fa.index /lab/Page_lab-users/Alex/gtex/index/hg38.gencode_v24.basic_ccds_nopar.transcripts.fa")
    os.system("python t2g.py < /lab/Page_lab-users/Alex/gtex/index/gencode.v24.annotation.basic_ccds_nopar.gtf > transcripts_to_genes.txt")


"""Define functions"""
        

def main(base_folder, fastq_1, fastq_2):
    
    os.system("kallisto bus -i hg38.gencode_v24.basic_ccds_nopar.transcripts.fa.index -o bus_output/{0}/ -x 10xv3 -t 10 {1} {2}".format(sample, fastq_1, fastq_2))


    #Process results
    os.chdir("{0}/bus_output/{1}".format(base_folder,sample))

    os.system("bustools correct -w ../../10xv3_whitelist.txt -o output.correct.bus output.bus")
    
    os.system("bustools sort -t 4 -o output.correct.sort.bus output.correct.bus")

    os.system("mkdir eqcount")
    os.system("mkdir genecount")

    os.system("bustools count -o eqcount/tcc -g ../../transcripts_to_genes.txt -e matrix.ec -m -t transcripts.txt output.correct.sort.bus")

    os.system("bustools count -o genecount/gene -g ../../transcripts_to_genes.txt -e matrix.ec -m -t transcripts.txt --genecounts output.correct.sort.bus")

    #Revise gene.genes.txt to have gene symbols instead of Ensembl IDs
    substitute_gene_symbols(base_folder, sample)
    

def substitute_gene_symbols(base_folder, sample):
    
    os.chdir("{0}/bus_output/{1}/genecount/".format(base_folder,sample))
    t2g = pd.read_csv("{0}/transcripts_to_genes.txt".format(base_folder), sep = "\t", header=None)
    t2g = t2g.drop_duplicates([1])
    
    ensg = []
    symbols = []
    
    for i in range(0,len(t2g.index)):
        symbol = t2g.iloc[i,2]
        ensg.append(t2g.iloc[i,1])
        if symbol not in symbols:
            symbols.append(symbol)
        else:
            symbols.append("{0}.1".format(symbol))
            
    t2g_dict = dict(zip(ensg,symbols)) 
    with open("gene.genes.txt","r") as old:
        with open("gene.genes.symbols.txt","w") as new:
            for line in old:
                new.write("{0}\n".format(t2g_dict[line.strip()]))
                
    os.system("mv gene.genes.txt gene.genes.original.txt")
    os.system("mv gene.genes.symbols.txt gene.genes.txt")



"""Run"""
if __name__ == "__main__":
    
    main(base_folder,fastq_1, fastq_2)
        
        
