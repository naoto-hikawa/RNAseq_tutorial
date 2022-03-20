# install DESeq2
BiocManager::install("DESeq2", force=TRUE)
library(DESeq2)

# Get file and set working directory
get_file=file.choose(new=F)
name_file=basename(get_file)
pathname=substr(get_file,1,nchar(get_file)-nchar(name_file))
setwd(pathname)　
name=substr(name_file,1,nchar(name_file)-4)

# Read the modified version of read counts
c <- read.table(name_file, header=T, sep="\t", quote="")
c <- data.frame(c, row.names = 1)
count <- as.matrix(c)

# Grouping: how to compare?
group <- data.frame(con = factor(c(rep("Group1", 4), rep("Group2", 4))))  #Change here

# DESeq analysis
deseq <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~con) 
deseq <- DESeq(deseq,fitType='local')

# Conversion table of Flybase ID, gene symbol, and CG numbers
data_geneid <- rownames(c)

# Place your downloaded gene_conversion_r6.44.csv in $SEQ_HOME
conversion_table <- read.table("gene_conversion_r6.44.csv", 
                               header=T, sep=",", quote="", stringsAsFactors=F)

df_conversion = as.data.frame(matrix(ncol = 2, nrow = 0))
for (i in 1:nrow(c)) {
  temp = as.data.frame(matrix(ncol = 2, nrow = 0))
  colnum = which(conversion_table$GENE_ID==data_geneid[i])
  temp <- conversion_table[colnum,2:3]
  df_conversion <- rbind(df_conversion, temp)
}

data_geneid <- cbind(data_geneid, df_conversion)
colnames(data_geneid) <- c("GENE_ID", "SYMBOL", "ANNOTATION_ID")

# Extract result information
res <- results(deseq, cooksCutoff=F, independentFiltering = F)

# Bind the geneid table with results
all_stats <- cbind(data_geneid, res)

# order by FDR and save all stats information including non-DEGs
resOrdered <- all_stats[order(all_stats$padj),]
stats=paste(name, "_all_stats.txt", sep="")
write.table(resOrdered, file=paste("results/", stats, sep=""), sep="\t", append=F, quote=F, row.names=F)

# Keep all genes with FDR < 0.05 AND baseMean >10
resSig <- subset(resOrdered, padj < 0.05) 
resSIG_FC_up <- subset(resSig, (log2FoldChange) > 0.58) 
resSIG_FC_down <- subset(resSig, (log2FoldChange) < -0.58) 
resSIG_FC <- subset(resSig, abs(log2FoldChange) > 0.58) 

# Write all DEGs to text and R data file
DEG_res=paste(name, "_all_DEGs.txt", sep="")
DEG_res_rda=paste(name, "_all_DEGs.rda", sep="")
write.table(resSig, file=paste("results/", DEG_res, sep=""), sep="\t", append=F, quote=F, row.names=F)
save(resSig, file=paste("results/", DEG_res_rda, sep=""))

# Write DEGs with absolute Fold Change > 1.5
sig_FC=paste(name, "_DEGs_FC.txt", sep="")
write.table(resSIG_FC, file=paste("results/", sig_FC, sep=""), sep="\t", append=F, quote=F, row.names=F)

# Write DEGs with Fold Change > 1.5
sig_FC_up=paste(name, "_DEGs_FC_up.txt", sep="")
write.table(resSIG_FC_up, file=paste("results/", sig_FC_up, sep=""), sep="\t", append=F, quote=F, row.names=F)

# Write DEGs with inverse Fold Change > 1.5
sig_FC_down=paste(name, "_DEGs_FC_down.txt", sep="")
write.table(resSIG_FC_down, file=paste("results/", sig_FC_down, sep=""), sep="\t", append=F, quote=F, row.names=F)

# Create normalized count data file
DESeq2_normalize <- counts(deseq, normalized=TRUE) 
DESeq2_normalize <- cbind(data_geneid, DESeq2_normalize)
normalized_file=paste(name, "_normalized.txt", sep="")
write.table(DESeq2_normalize, file=paste("results/", normalized_file, sep=""), sep="\t", append=F, quote=F, row.names=F) 

# Create file suited for GSEA, filter expressed genes only
GSEA = data.frame(matrix(nrow=nrow(all_stats), ncol = 0))
GSEA$Geneid <- all_stats$ANNOTATION_ID
GSEA$description <- all_stats$SYMBOL
GSEA <- cbind(GSEA, DESeq2_normalize[4:11])
GSEA$baseMean <- all_stats$baseMean

GSEASig <- subset(GSEA, baseMean >= 10)

GSEA_DEG =paste(name, "_GSEA_DEG.txt", sep="")
write.table(GSEASig[1:10], file=paste("results/", GSEA_DEG, sep=""), sep="\t", append=F, quote=F, row.names=F)

### Code written/edited by Naoto Hikawa, 2022.03.20
