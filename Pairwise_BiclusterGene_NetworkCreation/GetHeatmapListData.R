#####################################################################################################################
#
#  GetHeatmapListData.R
#
#  Extract data used to generate the heatmap at the end of the Seurat biclustering process for analysis and 
#  plotting.  Most code was take from Seurat library GetHeatmapData function.
#
#  Date         Name          Comment
#  7/10/2018  Kari Palmier   File Created
#  8/29/2018  Kari Palmier   File commented
#
#
#####################################################################################################################
library(reshape)
library(dplyr)

SetIfIsNull <- function(x, default) {
  if(is.null(x = x)){
    return(default)
  } else {
    return(x)
  }
}

GetHeatmapListData <- function(
  object,
  data.use = NULL,
  use.scaled = TRUE,
  cells.use = NULL,
  genes.use = NULL,
  disp.min = -2.5,
  disp.max = 2.5,
  group.by = "ident",
  group.order = NULL,
  assay.type = "RNA"
) {
  
  if (is.null(x = data.use)) {
    if (use.scaled) {
      data.use <- GetAssayData(object,assay.type = assay.type,slot = "scale.data")
    } else {
      data.use <- GetAssayData(object,assay.type = assay.type,slot = "data")
    }
  }
  
  # note: data.use should have cells as column names, genes as row names
  cells.use <- SetIfIsNull(x = cells.use, default = object@cell.names)
  cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  
  genes.use <- SetIfIsNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  } else {
    cells.ident <- factor(x = FetchData(
      object = object,
      cells.use = cells.use,
      vars.all = group.by
    )[, 1])
    names(x = cells.ident) <- cells.use
  }
  
  cells.ident <- factor(
    x = cells.ident,
    labels = intersect(x = levels(x = cells.ident), y = cells.ident)
  )
  
  data.use <- data.use[genes.use, cells.use, drop = FALSE]
  if ((!use.scaled)) {
    data.use = as.matrix(x = data.use)
    if (disp.max==2.5) disp.max = 10;
  }
  
  data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  
  data.use <- as.data.frame(x = t(x = data.use))
  
  data.use$cell <- rownames(x = data.use)
  
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  
  data.use %>% melt(id.vars = "cell") -> data.use
  
  names(x = data.use)[names(x = data.use) == 'variable'] <- 'gene'
  
  names(x = data.use)[names(x = data.use) == 'value'] <- 'expression'
  
  data.use$ident <- cells.ident[data.use$cell]
  
  if(!is.null(group.order)) {
    if(length(group.order) == length(levels(data.use$ident)) && all(group.order %in% levels(data.use$ident))) {
      data.use$ident <- factor(data.use$ident, levels = group.order)
    }
    else {
      stop("Invalid group.order")
    }
  }
  
  breaks <- seq(
    from = min(data.use$expression),
    to = max(data.use$expression),
    length = length(x = PurpleAndYellow()) + 1
  )
  
  data.use$gene <- with(
    data = data.use,
    expr = factor(x = gene, levels = rev(x = unique(x = data.use$gene)))
  )
  
  data.use$cell <- with(
    data = data.use,
    expr = factor(x = cell, levels = cells.use)
  )

  temp_data = list(heatmap_cell = data.use$cell,
                   heatmap_gene = data.use$gene,
                   heatmap_expression = data.use$expression
                   
   )
  
  return(temp_data)
  
}
