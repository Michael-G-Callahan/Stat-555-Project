library(DESeq2)
library(knitr)
library(ggplot2)
library(edgeR)
library(VennDiagram)
library(gridExtra)
library(EDASeq)
library(gprofiler2)
library(here)


get_DE_results <- function(selected_lines, colData, countData){

  selected_samples <- colData %>%
    filter(group %in% selected_lines) %>%
    mutate(group = factor(group, levels = selected_lines)) %>%
    arrange(group) %>%
    pull(source_name)

  designFormula <- "~ group"

  #create a DESeq dataset object from the count matrix and the colData
  dds <- DESeqDataSetFromMatrix(
    countData = countData[,selected_samples],
    colData = colData[match(selected_samples, colData$source_name),],
    design = as.formula(designFormula)
  )

  #For each gene, we count the total number of reads for that gene in all samples
  #and remove those that don't have at least 1 read.
  dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]

  dds <- DESeq(dds)

  #compute the contrast for the 'group' variable where 'CTRL'
  #samples are used as the control group.
  DEresults = results(dds, contrast = c("group", selected_lines))
  DEresults = DEresults[order(DEresults$pvalue),]

  ### -- Diagnostic Plots --
  #MA plot
  DESeq2::plotMA(object = dds, ylim = c(-5, 5))


  #Pvalue plot
  pval_plot <- ggplot(data = as.data.frame(DEresults), aes(x = padj)) +
    geom_histogram(bins = 100)
  plot(pval_plot)

  #Normalized counts plot
  new_x_labels = sub("_RNA", "", selected_samples)

  plotRLE(DESeq2::counts(dds, normalized = TRUE),
          outline=FALSE, ylim=c(-4, 4),
          col = as.numeric(colData$group),
          main = 'Normalized Counts',
          xaxt = "n"
          )
  axis(side = 1, at = 1:4, labels = new_x_labels)

  current_project_path <- here()
  DE_dir_path <- file.path(current_project_path, "Data",
                           paste0("DE_results_",selected_lines[1], "_", selected_lines[2],".Rdata"))
  save(DEresults, file = DE_dir_path)

  #Volcano plot
  DE_data <- DEresults@listData %>%
    as.data.frame %>%
    mutate(diffexpressed = "NO") %>%
    mutate(diffexpressed = ifelse(log2FoldChange > 0.6 & padj < 0.01, "UP", diffexpressed)) %>%
    mutate(diffexpressed = ifelse(log2FoldChange < -0.6 & padj < 0.01, "DOWN", diffexpressed)) %>%
    mutate(gene_symbol = DEresults@rownames)


  volcano_plot <- ggplot(data = DE_data, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 2) +
    coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) +
    scale_color_manual(values = c("#00AFBB", "grey", "#FFDB6D"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    ggtitle(paste("Volcano Plot:", selected_lines[1], "vs.", selected_lines[2]))
  plot(volcano_plot)

  return(DEresults)
}

get_lv_results <- function(selected_lines, colData, countData){
  #Create DGEList object
  d0 <- DGEList(countData[,selected_samples])

  #Add normalizing factors
  d0 <- calcNormFactors(d0)

  #Drop low-expressed genes
  cutoff <- 5
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,]
  dim(d) # number of genes left

  group = as.character(colData[match(selected_samples, colData$source_name),]$group)

  mm <- model.matrix(~group)
  colnames(mm) <- gsub("group", "", colnames(mm))

  y <- voom(d, mm, plot = T)
  fit <- lmFit(y, mm)
  tmp <- eBayes(fit)

  top.table <- topTable(tmp, sort.by = "P", n = Inf)


  pval_plot <- ggplot(data = top.table, aes(x = adj.P.Val)) +
    geom_histogram(bins = 100)
  plot(pval_plot)

  return(top.table)
}

compare_de_lv <- function(DEresults, lv_results){
  current_project_path <- here()
  venn_dir_path <- file.path(current_project_path, "Data",
                             paste0("venn_diagramm",selected_lines[1], "_", selected_lines[2],".png"))

  #Get DE significant genes
  DE <- DEresults[!is.na(DEresults$padj),]
  lowest_pvalue_indexes <- order(DE@listData$pvalue)[1:1000]
  DE2_dif_ex_genes = DE@rownames[lowest_pvalue_indexes]

  #Get LV significant genes
  lv_dif_ex_genes = rownames(head(lv_results, 1000))

  # Create a Venn diagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.plot <- venn.diagram(
    x = list(lv_dif_ex_genes, DE2_dif_ex_genes),
    category.names = c("Limma-voom" , "DESEQ2"),
    filename = venn_dir_path,
    output=TRUE,
    disable.logging = TRUE
  )
}

get_GO_results <- function(DE_results){
  #remove genes with NA values
  DE <- DEresults[!is.na(DEresults$padj),]
  #select genes with adjusted p-values below 0.1
  DE <- DE[DE$padj < 0.1,]
  #select genes with log2 fold change above 1 (two-fold change)
  up_reg_DE <- DE[DE$log2FoldChange > 1,]
  dn_reg_DE <- DE[DE$log2FoldChange < -1,]

  #get the list of upregulated genes of interest
  up_reg_gene_names <- rownames(up_reg_DE)
  up_reg_gene_names <- sapply(strsplit(up_reg_gene_names, "\\."), function(x) x[[1]])
  up_reg_gene_names <- unique(up_reg_gene_names)

  #get the list of downregulated genes of interest
  dn_reg_gene_names <- rownames(dn_reg_DE)
  dn_reg_gene_names <- sapply(strsplit(dn_reg_gene_names, "\\."), function(x) x[[1]])
  dn_reg_gene_names <- unique(dn_reg_gene_names)

  #calculate enriched GO terms
  up_go_response <- gost(query = up_reg_gene_names,
                         organism = 'mmusculus',
                         sources = c("GO"))
  dn_go_response <- gost(query = dn_reg_gene_names,
                         organism = 'mmusculus',
                         sources = c("GO"))

  # gostplot(up_go_response, capped=FALSE)
  up_go_results = up_go_response$result
  dn_go_results = dn_go_response$result

  # https://mikulskibartosz.name/f1-score-explained
  up_go_results = up_go_results %>%
    mutate(f1 = 2*(precision * recall) / (precision + recall)) %>%
    mutate(f2 = 2*(precision * recall) / (400*precision + recall)) %>%
    arrange(desc(f1))

  dn_go_results = dn_go_results %>%
    mutate(f1 = 2*(precision * recall) / (precision + recall)) %>%
    mutate(f2 = 2*(precision * recall) / (400*precision + recall)) %>%
    arrange(desc(f1))

  current_project_path <- here()
  up_go_path <- file.path(current_project_path, "Data",
                             paste0("GO_term_up_",selected_lines[1], "_", selected_lines[2],".csv"))

  dn_go_path <- file.path(current_project_path, "Data",
                          paste0("GO_term_dn_",selected_lines[1], "_", selected_lines[2],".csv"))

  up_go_results <- subset(up_go_results, select = -parents)
  dn_go_results <- subset(dn_go_results, select = -parents)

  write.csv(up_go_results, file = up_go_path)
  write.csv(dn_go_results, file = dn_go_path)

  return(list(up_reg = up_go_results, down_reg = dn_go_results))
}
