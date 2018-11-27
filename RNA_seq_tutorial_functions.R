## collections of functions needed for RNA-seq analysis

aggregate_func = function(x) {
  
  ## to be used inside sapply(df, aggregate_func) in order to:
    # concatenate characters into strings
    # calculate mean for numerics
    # NA will stay NA
  
  if (x %in% c(Inf, -Inf)) toString(mean(x))
  else if (is.numeric(x)) 
    round(mean(x),2)
  else if (x %in% c(TRUE, FALSE)) x[[1]][1]
  else if (!is.na(x))
    toString(x)
  else x
}


data_diagnosis <- function(counts) {
  
  require(e1071)
  require(reshape2)
  require(ggplot2)
  require(cowplot)
  
  # prints skewness & kurtosis of a RNA-seq raw counts
  # display a density and a boxplot on the same canvas
  
  counts = data.frame(counts)
  cat("skewness is:\n", sapply(counts, skewness), "\n")
  cat("kurtosis is:\n", sapply(counts, kurtosis), "\n")
  df = melt(counts, id.vars = NULL)
  density_plot <- ggplot(df, aes(x = value, color = variable, fill = variable)) + theme(legend.position="none") +
    geom_density(alpha = 0.2, size = 1.25) + ggtitle("density plot")
  
  boxplot <- ggplot(df, aes(x = variable, y = value, fill = variable)) + geom_boxplot() + 
    ylab("counts") + ggtitle("boxplot")
  
  plot_grid(density_plot, boxplot, labels = "AUTO", align = "v")
}


draw_PCA <- function(counts, ntop=20, title="PCA analysis", palatte="Spectral") {
  
  require(graphics)
  require(RColorBrewer)
  
  # draw PCA plot for a maximum of 10 samples
  
  counts = as.matrix(counts)
  colors = colorRampPalette(brewer.pal(11, palatte))(10)
  rv = rowVars(counts)
  select = order(rv, decreasing = TRUE)[1:ntop]
  pca = prcomp(t(counts[select,]))
  variance = ((pca$sdev^2) / (sum(pca$sdev^2)))*100
  min = min(pca$x[,1:2]) * 1.1
  max = max(pca$x[,1:2]) * 1.1
  plot(pca$x, type="n", main = title,
       xlab = paste("PC1: ", round(variance[1], 2), "%"),
       ylab = paste("PC2: ", round(variance[2], 2), "%"),
       xlim = c(min, max),
       ylim = c(min, max),
       cex.lab = 1.5,
       cex.main = 2)
  points(pca$x, col = colors, pch = 16, cex = 2)
  text(pca$x[,1], pca$x[,2], labels = colnames(counts), cex = 1.5, pos = 4)
}


draw_corr_heatmap <- function(counts, show_cellnote = TRUE, palatte="GnBu",
                              title = "clustering sample-to-sample\n distance"){
  require(gplots)
  require(graphics)
  require(RColorBrewer)
  
  # this function is modified from Dr. Istvan Albert's scripts
  # draw heatmap showing correlation of samples, cell notes are displayed by default
  
  mat.dist = counts
  mat.dist = as.matrix(dist(t(mat.dist)))
  mat.dist = mat.dist/max(mat.dist)
  hmcol = colorRampPalette(brewer.pal(9, palatte))(16)
  if (show_cellnote == TRUE)
    heatmap.2(mat.dist, col = rev(hmcol),
              cellnote = round(mat.dist,2),
              notecex = 1.0,
              notecol = "black",
              main = title,
              trace = "none", density.info = "none",
              keysize = 1.0,
              key.par = list(mar=c(3,1,3,1))
    )
  else 
    heatmap.2(mat.dist,
              col = rev(hmcol),
              main = title,
              trace = "none", density.info = "none",
              keysize = 1.0,
              key.par = list(mar=c(3,1,3,1))
    )
}

join_id <- function(dflist){
  
  require(plyr)
  
  # join DE analysis data by id to find common
  # Please make sure all data frames in the list contains an "id" column!
  ## Usage e.g. dflist = list(DESeq=df1, DESeq2=df2)
  # the double square bracket "[[]]" helps extract an element within a list directly
  
  n = length(dflist)
  merged_table = join_all(dflist, by = "id", type = "full")
  merged_table = subset(merged_table, select = c("id"))
  for (i in 1:n){
    query_df = dflist[c(i)]
    merged_table[,names(dflist[c(i)])] = merged_table$id %in% dflist[[i]]$id
  }
  rownames(merged_table) = merged_table$id
  
  common_id <- merged_table[rowSums(merged_table[2:ncol(merged_table)]) == n,]$id
  
  summary = rowSums(merged_table[,2:ncol(merged_table)])
  
  message <- paste(c("There are ", toString(length(summary[summary==4])), " genes DE in all 4 methods ",  
                    toString(length(summary[summary==3])), " genes in 3, ", 
                    toString(length(summary[summary==2])), " genes in 2, ",
                    toString(length(summary[summary==1])), " genes in 1."), collapse = "")
  
  list(merged_table=merged_table, common_id=common_id, message=message)
}


draw_venndiagram <- function(dflist, merged_table, alpha = 0.3, title = "Venn diagram"){
  
  require(VennDiagram)
  
  # draws a venndiagram 
  n = length(dflist)
  vdlist <- dflist
  
  for (i in 1:n){
    vdlist[[i]] <- rownames(merged_table[merged_table[,names(dflist[c(i)])] == TRUE,])
  }
  
  grid.newpage()
  vp <- venn.diagram(vdlist, fill = 2:5, alpha = alpha, filename = NULL, main = title)
  grid.draw(vp)
  
}

assign_color <- function(x, threshold = 0.05) {
  
  # This function is called upon when drawing MA and volcano plot in order to differentiate color of DE (red) & non-DE (black) dots
  
    if (is.na(x)) "black"
    else if (x >= threshold) "black"
    else "red"
}

DE_heatmap <- function(nc, common_id=NULL, use_jitters=TRUE, km=2,
                       cluster_columns=FALSE, title="heatmap"){
  
  require(ComplexHeatmap)
  require(RColorBrewer)
  require(graphics)
  
  # Draws a heatmap showing normalized read counts in each library for all DEGs
  # Provides a good way of visualizing similarities among biological replicates
  # Most importantly, it gives an idea if the normalization is good
  
  if (!is.null(common_id)) nc = subset(nc, id %in% common_id)
  gene = nc[,1]
  vals = as.matrix(nc[,2:11])
  if (use_jitters == TRUE) {vals = jitter(vals, factor = 1, amount = 0.00001)}
  score = NULL
  for (i in 1:nrow(vals)) {
    row=vals[i,]
    zscore = (row-mean(row))/sd(row)
    score = rbind(score,zscore)
  }
  row.names(score) = gene
  zscore=score
  mat = as.matrix(zscore)
  colors = colorRampPalette(c("green","black","red"),space="rgb")(256)
  Heatmap(mat,
          km = km,
          col = colors,
          column_title = title,
          name = "",
          cluster_columns = cluster_columns,
          show_row_names = FALSE,
          show_heatmap_legend = FALSE,
          row_dend_side = "left"
  )
}

draw_MA <- function(res, type = "DESeq", title = paste(type, "MA plot"),
                    pch = 16, cex = 0.5, xlab = "log normalized counts",
                    ylab = "logFC", xlim = c(-2, 18), ylim = c(-5,5)) {
  
  
  res = data.frame(res)
  if (type == "DESeq" | type == "DESeq2") {
    x = log2(res$baseMean)
    y = res$log2FoldChange
    z = res$padj
  } else if (type == "edgeR") {
    x = res$logCPM
    y = res$logFC
    z = res$FDR
  } else if (type == "limma") {
    x = res$AveExpr
    y = res$logFC
    z = res$adj.P.Val
  } else return ("wrong type")
  plot(x , y,
       main = title,
       pch = pch, cex = cex,
       xlab = xlab, 
       ylab = ylab,
       xlim = xlim,
       ylim = ylim,
       col = sapply(z, assign_color),
       cex.lab = 1.5,
       cex.main = 2
       )
}
  

  

draw_volcano <- function(res, type = "DESeq", 
                         pch = 16, cex = 0.5,
                         xlim = c(-4,4), ylim = c(0,20),
                         plotly = FALSE) {
    res = data.frame(res)
    if (type == "DESeq") {
      x = res$log2FoldChange
      y = res$pval
      z = res$padj
    } else if (type == "DESeq2") {
      x = res$log2FoldChange
      y = res$pvalue
      z = res$padj
    } else if (type == "edgeR") {
      x = res$logFC
      y = res$PValue
      z = res$FDR
    } else if (type == "limma") {
      x = res$logFC
      y = res$P.Value
      z = res$adj.P.Val
    } else return ("wrong type")
      plot(x, -log10(y),
           main = paste("Volcano plot for", type),
           pch = pch, cex = cex,
           xlab = expression(log[2]~fold~change),
           ylab = expression(-log[10]~pvalue),
           xlim = xlim,
           ylim = ylim,
           col = sapply(z, assign_color),
           cex.lab = 1.5,
           cex.main = 2
    )
}

