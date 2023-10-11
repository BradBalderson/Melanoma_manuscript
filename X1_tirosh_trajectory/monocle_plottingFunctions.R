
utils::globalVariables(c("Pseudotime", "value", "ids", "prin_graph_dim_1", "prin_graph_dim_2", "State", 
                         "value", "feature_label", "expectation", "colInd", "rowInd", "value", 
                         "source_prin_graph_dim_1", "source_prin_graph_dim_2"))

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' Plots clusters of cells .
#'
#' @param cds CellDataSet for the experiment
#' @param x the column of reducedDimS(cds) to plot on the horizontal axis
#' @param y the column of reducedDimS(cds) to plot on the vertical axis
#' @param color_by the cell attribute (e.g. the column of pData(cds)) to map to each cell's color
#' @param markers a gene name or gene id to use for setting the size of each cell in the plot
#' @param show_cell_names draw the name of each cell in the plot
#' @param cell_size The size of the point for each cell
#' @param cell_name_size the size of cell name labels
#' @param ... additional arguments passed into the scale_color_viridis function
#' @return a ggplot2 plot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom viridis scale_color_viridis
#' @export
#' @examples
#' \dontrun{
#' library(HSMMSingleCell)
#' HSMM <- load_HSMM()
#' HSMM <- reduceD
#' plot_cell_clusters(HSMM)
#' plot_cell_clusters(HSMM, color_by="Pseudotime")
#' plot_cell_clusters(HSMM, markers="MYH3")
#' }
plot_cell_clustersAlt <- function(cds, 
                               x=1, 
                               y=2, 
                               color_by="Cluster", 
                               markers=NULL, 
                               show_cell_names=FALSE, 
                               cell_size=1.5,
                               cell_name_size=2, 
                               ...){
  if (is.null(cds@reducedDimA) | length(pData(cds)$Cluster) == 0){
    stop("Error: Clustering is not performed yet. Please call clusterCells() before calling this function.")
  }
  
  gene_short_name <- NULL
  sample_name <- NULL
  data_dim_1 <- NULL
  data_dim_2 <- NULL
  
  #TODO: need to validate cds as ready for this plot (need mst, pseudotime, etc)
  lib_info <- pData(cds)
  
  tSNE_dim_coords <- reducedDimA(cds)
  data_df <- data.frame(t(tSNE_dim_coords[c(x,y),]))
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- colnames(cds)
  data_df <- merge(data_df, lib_info, by.x="sample_name", by.y="row.names")
  
  markers_exprs <- NULL
  if (is.null(markers) == FALSE){
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1){
      cds_subset <- cds[row.names(markers_fData),]
      if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        integer_expression <- TRUE
      }
      else {
        integer_expression <- FALSE
        
      }
      if (integer_expression) {
        cds_exprs <- exprs(cds_subset)
        
        if (is.null(sizeFactors(cds_subset))) {
          stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
        }
        cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
        
        cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
      }
      else {
        cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
      }
      markers_exprs <- cds_exprs
      #markers_exprs <- reshape2::melt(as.matrix(cds_exprs))
      colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y="row.names")
      #print (head( markers_exprs[is.na(markers_exprs$gene_short_name) == FALSE,]))
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name", by.y="cell_id")
    
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + facet_wrap(~feature_label) 
  }else{
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) 
  }
  
  # FIXME: setting size here overrides the marker expression funtionality. 
  # Don't do it!
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0){
    g <- g + geom_point(aes(color=log10(value + 0.1)), size=I(cell_size), na.rm = TRUE) + 
      scale_color_viridis(name = paste0("log10(value + 0.1)"), ...)
  }else {
    g <- g + geom_point(aes_string(color = color_by), size=I(cell_size), na.rm = TRUE)
  }
  
  g <- g + 
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() + 
    xlab(paste("Component", x)) + 
    ylab(paste("Component", y)) +
    #theme(legend.position="top", legend.key.height=grid::unit(0.35, "in")) +
    #guides(color = guide_legend(label.position = "top")) +
    #theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(text = element_text(size = 5))
  g
}

#'  Create a heatmap to demonstrate the bifurcation of gene expression along two branchs
#'  
#'  @description returns a heatmap that shows changes in both lineages at the same time. 
#'  It also requires that you choose a branch point to inspect. 
#'  Columns are points in pseudotime, rows are genes, and the beginning of pseudotime is in the middle of the heatmap. 
#'  As you read from the middle of the heatmap to the right, you are following one lineage through pseudotime. As you read left, the other. 
#'  The genes are clustered hierarchically, so you can visualize modules of genes that have similar lineage-dependent expression patterns.
#'
#' @param cds_subset CellDataSet for the experiment (normally only the branching genes detected with branchTest)
#' @param branch_point The ID of the branch point to visualize. Can only be used when reduceDimension is called with method = "DDRTree".
#' @param branch_states The two states to compare in the heatmap. Mutually exclusive with branch_point. 
#' @param branch_labels The labels for the branchs. 
#' @param cluster_rows Whether to cluster the rows of the heatmap.
#' @param hclust_method The method used by pheatmap to perform hirearchical clustering of the rows. 
#' @param num_clusters Number of clusters for the heatmap of branch genes
#' @param hmcols The color scheme for drawing the heatmap.
#' @param branch_colors The colors used in the annotation strip indicating the pre- and post-branch cells.
#' @param add_annotation_row Additional annotations to show for each row in the heatmap. Must be a dataframe with one row for each row in the fData table of cds_subset, with matching IDs.
#' @param add_annotation_col Additional annotations to show for each column in the heatmap. Must be a dataframe with one row for each cell in the pData table of cds_subset, with matching IDs.
#' @param show_rownames Whether to show the names for each row in the table.
#' @param use_gene_short_name Whether to use the short names for each row. If FALSE, uses row IDs from the fData table.
#' @param scale_max The maximum value (in standard deviations) to show in the heatmap. Values larger than this are set to the max.
#' @param scale_min The minimum value (in standard deviations) to show in the heatmap. Values smaller than this are set to the min.
#' @param norm_method Determines how to transform expression values prior to rendering
#' @param trend_formula A formula string specifying the model used in fitting the spline curve for each gene/feature.
#' @param return_heatmap Whether to return the pheatmap object to the user. 
#' @param cores Number of cores to use when smoothing the expression curves shown in the heatmap.
#' @param ... Additional arguments passed to buildBranchCellDataSet
#' @return A list of heatmap_matrix (expression matrix for the branch committment), ph (pheatmap heatmap object),
#' annotation_row (annotation data.frame for the row), annotation_col (annotation data.frame for the column). 
#' @import pheatmap
#' @importFrom stats sd as.dist cor cutree
#' @export
#'
plot_genes_branched_heatmapAlt <- function(cds_subset, 
                                        
                                        branch_point=1,
                                        branch_states=NULL,
                                        branch_labels = c("Cell fate 1", "Cell fate 2"), 
                                        cluster_rows = TRUE,
                                        hclust_method = "ward.D2", 
                                        num_clusters = 6,
                                        hmcols = NULL, 
                                        branch_colors = c('#979797', '#F05662', '#7990C8'), 
                                        add_annotation_row = NULL,
                                        add_annotation_col = NULL,
                                        show_rownames = FALSE, 
                                        use_gene_short_name = TRUE,
                                        scale_max=3, 
                                        scale_min=-3, 
                                        norm_method = c("log", "vstExprs"), 
                                        
                                        trend_formula = '~sm.ns(Pseudotime, df=3) * Branch',
                                        
                                        return_heatmap=FALSE,
                                        cores = 1, ...) {
  
  cds <- NA
  new_cds <- buildBranchCellDataSet(cds_subset, 
                                    branch_states=branch_states, 
                                    branch_point=branch_point, 
                                    progenitor_method = 'duplicate',
                                    ...)
  
  new_cds@dispFitInfo <- cds_subset@dispFitInfo
  
  if(is.null(branch_states)) {
    progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, 'State']
    branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
  }
  
  col_gap_ind <- 101
  # newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  # newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100))
  
  newdataA <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))   
  newdataB <- data.frame(Pseudotime = seq(0, 100,
                                          length.out = 100), 
                         Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
  
  BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores=cores, trend_formula = trend_formula,  
                                    relative_expr = T, new_data = rbind(newdataA, newdataB))
  
  BranchA_exprs <- BranchAB_exprs[, 1:100]
  BranchB_exprs <- BranchAB_exprs[, 101:200]
  
  #common_ancestor_cells <- row.names(pData(new_cds)[duplicated(pData(new_cds)$original_cell_id),])
  common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State == setdiff(pData(new_cds)$State, branch_states),])
  BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime'])))
  BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells, 'Pseudotime']))
  BranchB_num <- BranchA_num
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs') {
    BranchA_exprs <- vstExprs(new_cds, expr_matrix=BranchA_exprs)
    BranchB_exprs <- vstExprs(new_cds, expr_matrix=BranchB_exprs)
  }     
  else if(norm_method == 'log') {
    BranchA_exprs <- log10(BranchA_exprs + 1)
    BranchB_exprs <- log10(BranchB_exprs + 1)
  }
  
  heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1], BranchB_exprs)
  
  heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
  heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
  heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE,]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix>scale_max] = scale_max
  heatmap_matrix[heatmap_matrix<scale_min] = scale_min
  
  heatmap_matrix_ori <- heatmap_matrix
  heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch 
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  exp_rng <- range(heatmap_matrix) #bks is based on the expression range
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
  if(is.null(hmcols)) {
    hmcols <- blue2green2red(length(bks) - 1)
  }
  
  # prin  t(hmcols)
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=TRUE, 
                 show_rownames=F, 
                 show_colnames=F, 
                 #scale="row",
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 color=hmcols
                 #color=hmcols#,
                 # filename="expression_pseudotime_pheatmap.pdf",
  )
  #save(heatmap_matrix, row_dist, num_clusters, hmcols, ph, branchTest_df, qval_lowest_thrsd, branch_labels, BranchA_num, BranchP_num, BranchB_num, file = 'heatmap_matrix')
  
  #annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  
  if(!is.null(add_annotation_row)) {
    annotation_row <- add_annotation_row
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)), "Cell Type" = c(rep(branch_labels[1], BranchA_num),
                                                                                      rep("Pre-branch",  2 * BranchP_num),
                                                                                      rep(branch_labels[2], BranchB_num)))
  
  colnames(annotation_col) <- "Cell Type"  
  
  if(!is.null(add_annotation_col)) {
    annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])  
  }
  
  names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  
  annotation_colors=list("Cell Type"=branch_colors)
  
  names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
  
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  print(head(annotation_row))
  row.names(annotation_row) <- row_ann_labels
  
  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols=FALSE, 
                     cluster_rows=FALSE, 
                     show_rownames=show_rownames, 
                     show_colnames=F, 
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     annotation_colors=annotation_colors,
                     gaps_col = col_gap_ind,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=hmcols, 
                     border_color = NA,
                     silent=TRUE)
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(list(BranchA_exprs = BranchA_exprs, BranchB_exprs = BranchB_exprs, heatmap_matrix = heatmap_matrix, 
                heatmap_matrix_ori = heatmap_matrix_ori, ph = ph, col_gap_ind = col_gap_ind, row_dist = row_dist, hmcols = hmcols, 
                annotation_colors = annotation_colors, annotation_row = annotation_row, annotation_col = annotation_col, 
                ph_res = ph_res))
  }
}
