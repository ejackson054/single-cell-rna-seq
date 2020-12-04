import os
import pandas as pd

#Files you need in base folder:

#1) Transcriptome (.fasta)
#2) Barcodes (10xv3_whitelist.txt: downloaded from https://github.com/BUStools/getting_started/releases
#3) t2g.py: downloaded from https://github.com/pachterlab/kallisto-transcriptome-indices/releases


"""Assign variables"""

fastq_folder = "/lab/solexa_public/Page/201118_WIGTC-NOVASEQ1A_AHFWHGDSXY/FASTQ"
base_folder = "/lab/page_human_data/emily/single_cell"


"""Merge FASTQ files from different lanes"""

#Identify unique samples
samples = []

for i in os.listdir(fastq_folder):
    
    samples.append("{0}/{1}".format(fastq_folder,i.split("_L")[0]))
    
samples = list(set(samples))


#Merge samples into single FASTQ files
#os.system("mkdir {0}/fastq".format(base_folder))

#for sample in samples:
    
#    reads_1 = "{0}_L001_R1_001.fastq.gz {0}_L002_R1_001.fastq.gz {0}_L003_R1_001.fastq.gz".format(sample)
#    reads_2 = "{0}_L001_R2_001.fastq.gz {0}_L002_R2_001.fastq.gz {0}_L003_R2_001.fastq.gz".format(sample)
    
#    name = sample.split("/")[-1]
    
#    os.system("LSB_JOB_REPORT_MAIL=N bsub 'cat {0} > {1}/fastq/{2}_R1.fastq.gz'".format(reads_1, base_folder, name))
#    os.system("LSB_JOB_REPORT_MAIL=N bsub 'cat {0} > {1}/fastq/{2}_R2.fastq.gz'".format(reads_2, base_folder, name))


"""Run bustools"""

#Generate raw output
os.chdir(base_folder)


#os.system("kallisto index -i Homo_sapiens.GRCh38.cdna.v100.all.fa.index Homo_sapiens.GRCh38.cdna.v100.all.fa")

#os.system("python t2g.py --use_version < /lab/page/emily/Homo_sapiens.GRCh38.100.gtf > transcripts_to_genes.txt")


fastq_list = sorted(os.listdir("{0}/fastq".format(base_folder)))

to_drop = ["Undetermined_S0_R1.fastq.gz", "Undetermined_S0_R2.fastq.gz","L16_2577_S2_R1.fastq.gz","L16_2577_S2_R2.fastq.gz"]
fastq_list = [i for i in fastq_list if i not in to_drop]

fastq_list = " ".join(["{0}/fastq/{1}".format(base_folder,i) for i in fastq_list])

os.system("kallisto bus -i Homo_sapiens.GRCh38.cdna.v100.all.fa.index -o bus_output/ -x 10xv3 -t 10 {0}".format(fastq_list))


#Process results
os.chdir("bus_output")

os.system("bustools correct -w ../10xv3_whitelist.txt -o output.correct.bus output.bus")

os.system("bustools sort -t 4 -o output.correct.sort.bus output.correct.bus")

os.system("mkdir eqcount")
os.system("mkdir genecount")

os.system("bustools count -o eqcount/tcc -g ../transcripts_to_genes.txt -e matrix.ec -m -t transcripts.txt output.correct.sort.bus")

os.system("bustools count -o genecount/gene -g ../transcripts_to_genes.txt -e matrix.ec -m -t transcripts.txt --genecounts output.correct.sort.bus")

#Revise gene.genes.txt to have gene names
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
        
        