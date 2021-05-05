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
time <- paste("t", coldata$time, sep="")
rownames(coldata) <- paste(coldata$sample, time, coldata$replicate, sep=".")

# Set the last factor of contrast as a reference level.
coldata$condition <- relevel(coldata$condition, ref = snakemake@params[["contrast"]][2])

# Reorder columns in counts to match 'coldata'.
counts <- counts[, rownames(coldata)]
ddsFull <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)

# remove uninformative columns
ddsFull <- ddsFull[rowSums(counts(ddsFull)) > 1, ]
dds <- DESeq(ddsFull)
saveRDS(dds, file = snakemake@output[["rds"]])

# Contrast has three levels factor name and two levels.
contrast <- c("condition", snakemake@params[["contrast"]])
result <- results(dds, contrast=contrast)

# Coef has to be separated with underscores e.g. "condition_p7_vs_empty"
condition <- paste(snakemake@params[["contrast"]], collapse="_vs_")
coef <- paste("condition", condition, sep="_")

result <- lfcShrink(dds, coef=coef, res=result)
result <- result[order(result$padj),]

svg(snakemake@output[["ma_plot"]])
plotMA(result, ylim=c(-2,2))
dev.off()

out.table <- as.data.frame(result)
out.table <- cbind(geneid = rownames(out.table), out.table)

write.table(out.table, file=snakemake@output[["table"]], sep="\t", quote=FALSE, row.names=FALSE)
