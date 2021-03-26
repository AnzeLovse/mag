log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

counts <- read.csv(snakemake@input[["counts"]], sep="\t", row.names="Geneid")
counts <- subset(counts, select = -c(Chr, Start, End, Strand, Length))

samples <- read.csv(snakemake@params[["samples"]], sep="\t", stringsAsFactors = TRUE)
units <- read.csv(snakemake@params[["units"]], sep="\t", stringsAsFactors = TRUE)

coldata <- merge(samples, units, by='sample')
coldata <- subset(coldata, select=-c(fq1, fq2))
coldata[] <- lapply(coldata, as.factor)
rownames(coldata) <- paste(coldata$sample, coldata$replicate, sep=".")


ddsFull <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# remove uninformative columns
ddsFull <- ddsFull[rowSums(counts(ddsFull)) > 1, ]
dds <- DESeq(ddsFull)
saveRDS(dds, file = snakemake@output[[1]])