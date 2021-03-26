log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")


dds <- readRDS(snakemake@input[[1]])

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
