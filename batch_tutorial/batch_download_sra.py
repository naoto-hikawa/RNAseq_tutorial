# import some necessary tools
import subprocess
import os

# obtain environmental variables set on bash/zsh
SEQ_HOME = os.environ['SEQ_HOME']

# CHANGE HERE AND ENTER THE WHERE YOU WANT TO DOWNLOAD, DEFAULT IS $SEQ_HOME/ncbi/public/files/
dl_dir = SEQ_HOME + "/ncbi/public/files/"

# MAKE LIST WITH ACCESSION NUMBERS
sra_list = []

# WHAT IS THE FIRST ACCESSION NUMBER?
num = 9678461

for i in range(4):  # HOW MANY SAMPLES?
    sra = "SRR" + str(num)
    num = num + 1
    sra_list.append(sra)

for sra in sra_list:
    print("Downloading: " + sra)
    prefetch = "prefetch -v " + sra
    subprocess.call(prefetch, shell=True)

    print("Generating fastq for: " + sra)
    fasterq_dump = "fasterq-dump --split-3 -p -e 16 -O " + dl_dir + " -t " + SEQ_HOME + "/ncbi/ " + sra
    subprocess.call(fasterq_dump, shell=True)
