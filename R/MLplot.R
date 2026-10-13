#' Archimedean Spiral Heatmap for Machine Learning Performance
#'
#' Visualize machine learning model performance across training and validation
#' datasets using an Archimedean spiral heatmap. The top-performing models are
#' highlighted with radial bars.
#'
#' @param Object A machine learning result object containing performance matrices.
#' @param features Character. Performance metric to visualize, such as `"AUC"`,
#'   `"F1"`, `"ACC"`, `"Youden"`. `"RECALL"`, `"NPV"`, `"PPV"`, `"Sensitive"`, `"Specific"`, `"Kappa"`
#' @param top_n Integer. Number of top-performing models to highlight. Default is 10.
#' @param group_color Character vector specifying colors for datasets. If `NULL`,
#'   colors are automatically generated using the AAAS palette.
#' @param bar_col Character vector of three colors used for the top-model bar
#'   gradient. Default is `c("#A681C4", "#F8B538", "#E04646")`.
#' @param pheat_col Character vector of three colors used for the heatmap gradient.
#'   Default is `c("#5B9DD5", "white", "#D85E87")`.
#' @param circle Numeric. Number of spiral rotations. Default is 2.
#' @param plot_size Numeric. Size of the spiral plotting region in millimeters.
#'   Default is 140.
#' @param plot_x Numeric. Horizontal position of the spiral plotting region in
#'   normalized parent coordinates. Default is 0.4.
#' @param plot_y Numeric. Vertical position of the spiral plotting region in
#'   normalized parent coordinates. Default is 0.6.
#' @param bar_h Numeric. Baseline height of the bars highlighting the top models.
#'   Default is 1.5.
#' @param order Character vector specifying the display order of datasets.
#'   By default, the order is obtained from the selected performance matrix.
#'
#' @returns Invisibly returns the processed model-performance results while
#'   drawing an Archimedean spiral heatmap.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' MLplot(Object, features = "AUC")
#'
#' MLplot(
#'   Object,
#'   features = "AUC",
#'   top_n = 20,
#'   circle = 2,
#'   group_color = c("#5A6E9B", "#62C3C5", "#F8B538")
#' )
#' }
MLplot <- function(Object,features,top_n=10,group_color = NULL,bar_col=NULL,pheat_col=NULL,
                   circle=2, plot_size = 140,plot_x = 0.4,
                   plot_y = 0.6,bar_h=1.5,
                   order=colnames(do.call(rbind,Object[[features]]))){
  data <- Object[[features]]
  data <- do.call(rbind,data)
  data <- data[,c("Train",setdiff(order,"Train")), drop = FALSE]
  data <- data.frame(data)
  data <-na.omit(data)
  method_name  <- data.frame(method=rownames(data))
  get_gene_set <- function(x) {
    ifelse(grepl("_", x),trimws(gsub("_", " ", sub("_NN-MLP.*", "", x))),
           trimws(sub("\\+.*", "", x)))}
  method_name$Select  <- get_gene_set(method_name$method)
  get_method <- function(x) {
    ifelse(grepl("_NN-MLP", x),
           trimws(gsub("_", " ", sub("^.*?_(NN-MLP.*)$", "\\1", x))),
           trimws(sub("^.*?\\+", "", x)))}
  method_name$model <- get_method(method_name$method)
  method_name$Select_method <- sub("\\s+\\d+\\s+genes$", "", method_name$Select)
  method_name$Gene_number <- as.numeric(sub(".*\\s+(\\d+)\\s+genes$", "\\1", method_name$Select))
  clean_string <- function(text) {
    text <- trimws(text)
    text <- sub("\\(.*", "", text)
    text <- sub(" .*", "", text)
    text}
  method_name$model2 <-  clean_string(method_name$model)
  method_name$simple <- paste0(method_name$Select,"+",method_name$model2)
  if (ncol(data) < 2) {
    stop("MLplot requires at least one validation dataset.")
  }
  data_cleaned <- data
  data_cleaned$Mean <- rowMeans(data_cleaned[, -1, drop = FALSE],na.rm = TRUE)
  data_cleaned <- tibble::rownames_to_column(data_cleaned, var = "method")
  final_result <- dplyr::inner_join(data_cleaned, method_name, by = "method")
  final_result <- dplyr::group_by(final_result, rlang::.data$simple)
  final_result <- dplyr::arrange(final_result,
                                 dplyr::desc(rlang::.data$Mean),
                                 rlang::.data$method,
                                 .by_group = TRUE)
  final_result <- dplyr::slice(final_result, 1)
  final_result <- dplyr::ungroup(final_result)
  final_result <- tibble::column_to_rownames(final_result, var = "method")

  df_rank <- final_result
  rownames(df_rank) <- df_rank$simple
  df_rank <- df_rank[,c("Train",setdiff(order,"Train"),"Mean"),
                     drop = FALSE]
  df_rank <- df_rank[base::order(df_rank$Mean, decreasing = F),,drop = FALSE]
  mat <- as.matrix(df_rank )
  mat1 <-mat[,-ncol(mat), drop = FALSE]
  #######
  top_n <- min(top_n, nrow(df_rank))
  top_index <- seq.int(from = nrow(df_rank) - top_n + 1,to = nrow(df_rank))
  top_data <- df_rank[top_index, -ncol(mat), drop = FALSE]
  top_value <- df_rank $Mean[top_index]
  top_name <- rownames(df_rank )[top_index]
  #####
  bar_baseline <- ncol(mat1)
  bar_min <- min(top_value, na.rm = TRUE)
  bar_max <- max(top_value, na.rm = TRUE)
  bar_height <- if (bar_min == bar_max) {
    rep(bar_h, length(top_value))
  } else {
    (top_value - bar_min) /
      (bar_max - bar_min) * 2 + bar_h
  }

  bar_end <- bar_baseline + bar_height
  group_x <- nrow(mat1) + 1
  #####
  if (is.null(group_color)) {group_color <- ggsci::pal_aaas()(ncol(mat1))
  }
  if (is.null(bar_col)) {bar_col <- c("#A681C4", "#F8B538", "#E04646")
  }
  if (is.null(pheat_col)) {pheat_col <- c("#5B9DD5", "white", "#D85E87")
  }
  #####
  if (bar_min == bar_max) {
    eps <- max(abs(bar_min) * 0.001, 1e-6)
    bar_col_fun <- circlize::colorRamp2(c(bar_min - eps, bar_min, bar_min + eps),bar_col)
  } else {bar_col_fun <- circlize::colorRamp2(c(bar_min, mean(c(bar_min,
                                                                bar_max)),
                                                bar_max),
                                              bar_col)
  }
  #####
  mat_min <- min(mat1, na.rm = TRUE)
  mat_max <- max(mat1, na.rm = TRUE)
  mat_mid <- mean(mat1, na.rm = TRUE)
  if (mat_min == mat_max) {eps <- max(abs(mat_min) * 0.001, 1e-6)
    mat_breaks <- c(mat_min - eps, mat_min, mat_min + eps)
  } else {
    mat_breaks <- c(mat_min, mat_mid, mat_max)
  }
  col_fun <- circlize::colorRamp2(mat_breaks, pheat_col)
  #####
  spiralize::spiral_initialize(
    xlim = c(0.5, group_x + 0.5),
    start = 360,
    end = 360 * circle + 45,
    scale_by = "curve_length",
    padding = grid::unit(5, "mm"),
    vp_param = list(
      x = grid::unit(plot_x, "npc"),
      y = grid::unit(plot_y, "npc"),
      width = grid::unit(plot_size, "mm"),
      height = grid::unit(plot_size, "mm")))
  spiralize::spiral_track(ylim = c(0, ncol(mat1)),background = FALSE)
  #####
  for (j in seq_len(ncol(mat1))) {
    spiralize::spiral_rect(xleft = seq_len(nrow(mat1)) - 0.5,
                           ybottom = j - 1,xright = seq_len(nrow(mat1)) + 0.5,
                           ytop = j,
                           gp = grid::gpar(fill = col_fun(mat1[, j]),
                                           col = "black",lwd = 0.5))
  }

  #####
  for (j in seq_len(ncol(mat1))) {
    spiralize::spiral_rect(xleft = group_x - 0.5,ybottom = j - 1,
                           xright = group_x + 0.5,ytop = j,
                           gp = grid::gpar(fill = group_color[j],
                                           col = "black",lwd = 0.5))
  }

  #####
  spiralize::spiral_bars(pos = top_index,value = bar_end,
                         baseline = ncol(mat1),bar_width = 0.8,
                         gp = grid::gpar(fill = bar_col_fun(top_value),
                                         col = "white",lwd = 0.5)
  )

  #####
  spiralize::spiral_text(x = top_index,y = bar_end,
                         text = sprintf("%.2f", top_value),
                         facing = "clockwise",nice_facing = TRUE,
                         vjust = 0,offset = grid::unit(2, "mm"),
                         hjust = 0,
                         gp = grid::gpar(fontsize = 7,col = "#B2182B",fontface = "bold"))

  #####
  spiralize::spiral_text(x = top_index,y = bar_end,text = top_name,
                         facing = "clockwise",nice_facing = TRUE,
                         vjust = 0,offset = grid::unit(10, "mm"),hjust = 0,
                         gp = grid::gpar(fontsize = 7,col = "black"))
  #####
  lgd1 <- ComplexHeatmap::Legend(title = features,col_fun = col_fun)
  lgd2 <- ComplexHeatmap::Legend(title = "Top Models",col_fun = bar_col_fun)
  lgd_group <- ComplexHeatmap::Legend(title = "Dataset",labels = colnames(mat1),
                                      legend_gp = grid::gpar(fill = group_color[seq_along(colnames(mat1))],
                                                             col = "black",lwd = 0.5),
                                      type = "grid",
                                      grid_width = grid::unit(5, "mm"),
                                      grid_height = grid::unit(5, "mm"),
                                      labels_gp = grid::gpar(fontsize = 9),
                                      title_gp = grid::gpar(fontsize = 10,
                                                            fontface = "bold"),
                                      gap = grid::unit(2, "mm"))
  ComplexHeatmap::draw(ComplexHeatmap::packLegend(lgd_group, lgd1, lgd2),
                       x = grid::unit(0, "npc"),y = grid::unit(0.5, "npc"))
  invisible(df_rank)
}
