# import some necessary tools
import subprocess
import os

# obtain environmental variables set on bash/zsh
SEQ_HOME = os.environ['SEQ_HOME']

# CHANGE HERE AND ENTER THE DIRECTORY WITH FASTQ FILES, OR USE THE DEFAULT $SEQ_HOME/fastq/
work_dir = SEQ_HOME + "/downstream_practice/fastq/"

# list all files in working directory
file_names = os.listdir(work_dir)

sample_names = []

# filter unnecessary files such as .DS_Store
for name in file_names:
    if name.endswith(".gz"):
        sample_names.append(name)

sample_names = sorted(sample_names)

# create list of aligned bam so we count all alignments into 1 file
bam_list = ""

# loop through the list of fastq.gz
for file_name in sample_names:

    # define the actual file name without extension, remove ".fastq.gz"
    name = file_name[0:-9]

    # where did we put our fastq files?
    file_dir = work_dir + file_name

    # define result directory
    result_dir = SEQ_HOME + "/downstream_practice/results/batch/" + name

    # define Quality Check directory
    fastqc_dir = result_dir + "/qc"

    # make QC directory if it does not exist
    if os.path.exists(fastqc_dir):
        pass
    else:
        os.makedirs(fastqc_dir)

    # where to put the trimmed fastq?
    trimmed_file = result_dir + "/" + "trimmed.fastq.gz"

    # define command as a string so we can send to shell (Terminal)
    cutadapt = "cutadapt -m 20 -O 18 -a 'polyA=A{18}' -a 'QUALITY=G{18}' -n 2 " + file_dir +  \
               " -j 4 | " \
               "cutadapt -m 20 -O 3 --nextseq-trim=20 " \
               "-a 'r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000' " \
               "-a 'r2adapter=A{18}AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=3;max_error_rate=0.100000' " \
               "-n 2 - -j 4 | " \
               "cutadapt -m 20 -O 20 -g 'r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20' " \
               "--discard-trimmed -o " + trimmed_file + " - -j 4"

    # run cutadapt
    subprocess.call(cutadapt, shell=True)

    # define fastqc command, both original and trimmed files.
    fastqc = SEQ_HOME + "/tools/FastQC.app/Contents/MacOS/fastqc -o " + fastqc_dir + \
             " -t 8 --nogroup " + file_dir + " " + trimmed_file

    # run fastqc
    subprocess.call(fastqc, shell=True)

    # define where to put STAR results
    STAR_dir = result_dir + "/STAR_align"

    # make directory if it does not exist
    if os.path.exists(STAR_dir):
        pass
    else:
        os.makedirs(STAR_dir)

    # define STAR command
    STAR = "STAR --genomeDir " + SEQ_HOME + "/genome_index " \
           "--runThreadN 4 " \
           "--readFilesIn " + trimmed_file + " --readFilesCommand gunzip -c " \
           "--outFilterType BySJout " \
           "--outFilterMultimapNmax 20 " \
           "--alignIntronMin 20 " \
           "--alignIntronMax 1000000 " \
           "--alignMatesGapMax 1000000 " \
           "--alignSJoverhangMin 8 " \
           "--alignSJDBoverhangMin 1 " \
           "--outSAMattributes All " \
           "--outSAMtype BAM SortedByCoordinate " \
           "--outFileNamePrefix " + STAR_dir + "/"

    # run STAR
    subprocess.call(STAR, shell=True)

    # where is the sorted bam file after STAR alignment?
    sorted_bam = STAR_dir + "/Aligned.sortedByCoord.out.bam"

    # define samtools command to index the bam file
    samtools = "samtools index " + sorted_bam

    # call samtools
    subprocess.call(samtools, shell=True)

    # add the bam file absolute path to the list to save which files to count later
    bam_list = bam_list + sorted_bam + " "

# when all STAR and indexing are done, define where to put read count files
read_dir = SEQ_HOME + "/downstream_practice/results/batch/read_counts"

# make directory if it does not exist
if os.path.exists(read_dir):
    pass
else:
    os.makedirs(read_dir)

# define featureCounts command, notice that bam_list is used which means all files will be counted to 1 file
featureCounts = SEQ_HOME + "/tools/subread-2.0.3-macOS-x86_64/bin/featureCounts " \
                "-a " + SEQ_HOME + "/annotations/dmel-all-r6.44.gtf " \
                "-o " + read_dir + "/readcount.txt " + bam_list + "-T 4"

# run featureCounts
subprocess.call(featureCounts, shell=True)

# define cut command. columns 7~(7+sample number-1) are the sample counts, so we will cut these out
sample_num = len(sample_names)
range_end = 7 + sample_num - 1
cut = "cut -f 1,7-" + str(range_end) + " " + \
      read_dir + "/readcount.txt > " + \
      read_dir + "/readcount_clean.txt"

# run cut
subprocess.call(cut, shell=True)
