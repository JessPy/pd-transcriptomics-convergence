#!/usr/bin/env Rscript

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(openxlsx)

# Usage: Rscript de_analysis.R <counts.tsv> <metadata.csv> <outdir> <gene_name> <id_to_symbol.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript de_analysis.R <counts.tsv> <metadata.csv> <outdir> <gene_name> <id_to_symbol.tsv>")
}

counts_file <- args[1]
meta_file <- args[2]
outdir <- args[3]
gene_name <- args[4]
mapping_file <- args[5]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
counts <- round(counts)

meta <- read.csv(meta_file)
if ("Run" %in% colnames(meta)) {
    rownames(meta) <- meta$Run
} else {
    stop("Metadata must have a 'Run' column.")
}

# Sync samples
common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

# Load Mapping
mapping <- read.table(mapping_file, header = FALSE, sep = "\t", col.names = c("gene_id", "symbol"))
rownames(mapping) <- mapping$gene_id

# 2. Experimental Design
if (gene_name == "SNCA") {
    meta$genotype <- factor(gsub(" ", "_", meta$genotype))
    meta$comparison <- factor(meta$comparison)
    meta$Group <- factor(paste0(meta$genotype, "_", meta$comparison))
    design_formula <- ~ Group
} else if (gene_name == "GBA1" || gene_name == "GAB1") {
    # Specific design for GBA1 with three levels
    meta$comparison <- factor(meta$comparison, levels = c(
        "GBA1 +/+ (Controles Saudáveis / Wild Type)",
        "GBA1 IVS/+ (Mutante Heterozigoto)",
        "GBA1 IVS/IVS (Mutante Homozigoto Nulo)"
    ))
    design_formula <- ~ comparison
} else {
    meta$condition <- factor(meta$condition, levels = c("Control", "Mutant"))
    design_formula <- ~ condition
}

# 3. Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

# 4. Extract Results
results_list <- list()
if (gene_name == "SNCA") {
    # All groups: Triplication_SNCA_NSC, Triplication_SNCA_WTPFF, Triplication_SNCA_PDamp, Triplication_SNCA_MSAamp
    #             Control_NSC, Control_WTPFF, Control_PDamp, Control_MSAamp
    
    # Baseline Genotype Effect
    results_list[["Genotype_Triplication_vs_Control_NSC"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_NSC", "Control_NSC"))
    
    # Treatment effects in Control
    results_list[["Treatment_Control_WTPFF_vs_NSC"]] <- results(dds, contrast=c("Group", "Control_WTPFF", "Control_NSC"))
    results_list[["Treatment_Control_PDamp_vs_NSC"]] <- results(dds, contrast=c("Group", "Control_PDamp", "Control_NSC"))
    results_list[["Treatment_Control_MSAamp_vs_NSC"]] <- results(dds, contrast=c("Group", "Control_MSAamp", "Control_NSC"))
    
    # Treatment effects in Triplication
    results_list[["Treatment_Triplication_WTPFF_vs_NSC"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_WTPFF", "Triplication_SNCA_NSC"))
    results_list[["Treatment_Triplication_PDamp_vs_NSC"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_PDamp", "Triplication_SNCA_NSC"))
    results_list[["Treatment_Triplication_MSAamp_vs_NSC"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_MSAamp", "Triplication_SNCA_NSC"))

    # Genotype effect under each treatment
    results_list[["Genotype_Triplication_vs_Control_WTPFF"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_WTPFF", "Control_WTPFF"))
    results_list[["Genotype_Triplication_vs_Control_PDamp"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_PDamp", "Control_PDamp"))
    results_list[["Genotype_Triplication_vs_Control_MSAamp"]] <- results(dds, contrast=c("Group", "Triplication_SNCA_MSAamp", "Control_MSAamp"))

} else if (gene_name == "GBA1" || gene_name == "GAB1") {
    results_list[["Heterozygote_vs_Control"]] <- results(dds, contrast=c("comparison", "GBA1 IVS/+ (Mutante Heterozigoto)", "GBA1 +/+ (Controles Saudáveis / Wild Type)"))
    results_list[["Homozygote_vs_Control"]] <- results(dds, contrast=c("comparison", "GBA1 IVS/IVS (Mutante Homozigoto Nulo)", "GBA1 +/+ (Controles Saudáveis / Wild Type)"))
    results_list[["Homozygote_vs_Heterozygote"]] <- results(dds, contrast=c("comparison", "GBA1 IVS/IVS (Mutante Homozigoto Nulo)", "GBA1 IVS/+ (Mutante Heterozigoto)"))
} else {
    results_list[["Main_Effect_Mutant_vs_Control"]] <- results(dds)
}

# 5. Export Plots
vsd <- vst(dds, blind=FALSE)

# PCA Plot
p_pca <- plotPCA(vsd, intgroup=if(gene_name=="SNCA") c("genotype", "comparison") else if (gene_name == "GBA1" || gene_name == "GAB1") "comparison" else "condition") +
  theme_minimal() + labs(title = paste("PCA -", gene_name))
ggsave(file.path(outdir, "pca_plot.png"), p_pca, width=8, height=6, dpi=300)

for (res_name in names(results_list)) {
    res <- results_list[[res_name]]
    res_df <- as.data.frame(res)
    res_df$symbol <- mapping[rownames(res_df), "symbol"]
    res_df <- res_df[, c("symbol", setdiff(colnames(res_df), "symbol"))]
    
    # ── Excel Export with Tabs ──────────────────────────────────────────────
    res_df$status <- "Not Significant"
    res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up-regulated"
    res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down-regulated"
    
    up_genes <- res_df[res_df$status == "Up-regulated", ]
    down_genes <- res_df[res_df$status == "Down-regulated", ]
    
    # Order by significance
    res_df <- res_df[order(res_df$padj), ]
    up_genes <- up_genes[order(up_genes$padj), ]
    down_genes <- down_genes[order(down_genes$padj), ]

    wb <- createWorkbook()
    addWorksheet(wb, "All_Genes")
    addWorksheet(wb, "Up_regulated")
    addWorksheet(wb, "Down_regulated")
    
    writeData(wb, "All_Genes", res_df, rowNames = TRUE)
    writeData(wb, "Up_regulated", up_genes, rowNames = TRUE)
    writeData(wb, "Down_regulated", down_genes, rowNames = TRUE)
    
    saveWorkbook(wb, file.path(outdir, paste0("DE_results_", res_name, ".xlsx")), overwrite = TRUE)
    # Also save CSV for compatibility
    write.csv(res_df, file.path(outdir, paste0("DE_results_", res_name, ".csv")))

    # ── Volcano Plot ─────────────────────────────────────────────────────────
    res_df$label <- ifelse(res_df$padj < 1e-10 & abs(res_df$log2FoldChange) > 2, res_df$symbol, "")
    
    p_volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=status, label=label)) +
        geom_point(alpha=0.4) +
        geom_text_repel(max.overlaps = 15) +
        theme_minimal() +
        scale_color_manual(values=c("Down-regulated" = "blue", "Not Significant" = "grey", "Up-regulated" = "red")) +
        labs(title = paste("Volcano -", res_name),
             subtitle = "Thresholds: padj < 0.05 & |log2FC| > 1",
             color = "Status") +
        geom_vline(xintercept = c(-1, 1), linetype="dashed", color="black", alpha=0.3) +
        geom_hline(yintercept = -log10(0.05), linetype="dashed", color="black", alpha=0.3)
    
    ggsave(file.path(outdir, paste0("volcano_", res_name, ".png")), p_volcano, width=8, height=6, dpi=300)
}

# Heatmap with Symbols
top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
heatmap_mat <- assay(vsd)[top_genes,]
rownames(heatmap_mat) <- mapping[rownames(heatmap_mat), "symbol"]

ann_cols <- colnames(colData(dds))
ann_cols <- ann_cols[!ann_cols %in% c("sizeFactor", "Group", "replaceable")]
annotation_df <- as.data.frame(colData(dds)[, ann_cols, drop=FALSE])

png(file.path(outdir, "heatmap_top50_var.png"), width=10, height=12, units="in", res=300)
pheatmap(heatmap_mat, 
         annotation_col = annotation_df,
         show_colnames = FALSE, main = "Top 50 Most Variable Genes (Symbols)")
dev.off()

print(paste("Analysis complete for", gene_name))
