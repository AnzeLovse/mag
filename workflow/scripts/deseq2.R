log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

counts <- read.csv(snakemake@input[["counts"]], sep="\t", row.names="Geneid")
counts <- subset(counts, select = -c(Chr, Start, End, Strand, Length))

samples <- read.csv(snakemake@params[["samples"]], sep="\t", stringsAsFactors = TRUE)
units <- read.csv(snakemake@params[["units"]], sep="\t", stringsAsFactors = TRUE)

coldata <- merge(unique(samples), units, by='sample')
coldata <- subset(coldata, select=-c(fq1, fq2))
coldata[] <- lapply(coldata, as.factor)
rownames(coldata) <- paste(coldata$sample, paste("t",coldata$time,sep=""), coldata$replicate, sep=".")

# Speciffic steps needed for our analysis
#--------------------------------------------------
counts$G_M.t0.1 <- counts$G.t0.1
counts$G_M.t0.2 <- counts$G.t0.2

zero_rows <- data.frame(rbind(
  c("G_M", "mitomycin", 1, 0),
  c("G_M", "mitomycin", 2, 0)
))
colnames(zero_rows) <- colnames(coldata)
rownames(zero_rows) <- c("G_M.t0.1", "G_M.t0.2")
coldata <- rbind(coldata, zero_rows)
#--------------------------------------------------

# Filter rows with less than 5 counts in at least 2 samples
filter <- apply(counts, 1, function(x) length(x[x>5])>=2)
filtered <- counts[filter,]

# Create a DESeq data set with a full desing.
ddsFull <- DESeqDataSetFromMatrix(
    countData = filtered,
    colData = coldata,
    design = as.formula(snakemake@params[["model"]]))

# Use Likelihood-ratio test on a reduced model.
dds <- DESeq(ddsFull, test="LRT", reduced = as.formula(snakemake@params[["reduced_model"]]))
saveRDS(dds, file = snakemake@output[["rds"]])

alpha <- as.numeric(snakemake@params[["alpha"]])
result <- results(dds, alpha=alpha)

# Construct a list of shrunken LFC for the coefficients of the model.
logfc.list <- list()
logfc.list[[1]] <- data.frame(result)
for (i in 2:length(resultsNames(dds))) {
    coeff = resultsNames(dds)[i]
    res.t <- results(dds, name=coeff, test="Wald", alpha = alpha)
    res_shrunken <- lfcShrink(dds, coef=coeff, type="apeglm", res = res.t)

    out.table <- as.data.frame(res_shrunken)
    out.table <- cbind(geneid = rownames(out.table), out.table)
    write.table(
        out.table,
        file=paste("deseq2/", coeff, "_shrunken.tsv", sep=""),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    svg(paste("deseq2/", coeff, "ma_plot.svg", sep=""))
    plotMA(result, ylim=c(-2,2))
    dev.off()

    # Add the shrunken log2FC to the list.
    lfc <- res_shrunken[, "log2FoldChange", drop=FALSE]
    colnames(lfc) <- c(paste(coeff, "sLFC", sep="."))
    logfc.list[[i]] <- lfc
}

# Write combined results in the main output.
out.result <- do.call(cbind, logfc.list)
out.result <- out.result[order(out.result$padj),]
out.result <- cbind(geneid = rownames(out.result), out.result)
  write.table(
    out.result,
    file=snakemake@output[["table"]],
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )

# Use variance stabilising transformation and then plot PCA.
vsd <- vst(ddsFull, blind=FALSE)
svg(snakemake@output[["pca_plot"]])
plotPCA(vsd, intgroup=c("condition", "time"), ntop = 500)
dev.off()
