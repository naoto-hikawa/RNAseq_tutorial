if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

# Install stringr package
install.packages("stringr")                        
library("stringr") 

# Choose the clean.txt file from popup window
get_file=file.choose(new=F)

# Retrieve file name
name_file=basename(get_file)

# Retrieve parent directory name
pathname=substr(get_file,1,nchar(get_file)-nchar(name_file))

# Set retrieved directory name as working directory
setwd(pathname)　

# Read the file as table
c <- read.table(name_file, header=F, sep="\t", quote="")

# How many columns including the gene names?
numcol=ncol(c)

# We want to change the sample names so they do not include the whole path name
cnames=c[1,]

# Loop through each sample
for (i in 2:numcol) {
  # Start of actual sample name is after '/results/batch/'
  start_pos = str_locate(pattern = "downstream_practice/results/batch/", cnames[i]) [2] +1
  
  # End of actual sample name is before '/STAR_align/'
  end_pos = str_locate(pattern = "/STAR_align/Aligned.sortedByCoord.out.bam", cnames[i]) [1] -1
  
  # Fix the sample name column so that only sample name remains
  cnames[i]=substr(cnames[i], start_pos, end_pos)
}

# Make the col names as the truncated version, include Geneid
colnames(c)<-cnames

# When we read the table, first row was the header section and we don't need it anymore
c <- c[2:nrow(c),]

# It's not absolutely crucial but we can order the geneid 
c_ordered<-c[order(c$Geneid),]

# Write the cleaner version to the working directory
write.table(c_ordered, "cleaner.txt", sep="\t", append=F, quote=F, row.names=F)
