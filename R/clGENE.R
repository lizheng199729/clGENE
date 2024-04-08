#'The potential molecular mechanisms behind the similar phenotypes of different diseases were explored
#' This function accepts a variable number of gene expression matrices and performs
#' a series of analyses including normalization, PCA, and plotting of PCA results.
#' Optionally, it can identify and export closely related genes based on Euclidean distance.
#'
#' @param ... Variable argument list, each representing a gene expression matrix.
#'
#' @return A list of PCA result objects for each gene expression matrix.
#'
#' @examples
#' # Assuming 'matrix1' and 'matrix2' are available gene expression matrices:
#' result <- clGENE(matrix1, matrix2)
#'
#' @export
#########################################
select_top_n<-function(scores,n_top){
  d <- data.frame(
    x   = data.table::copy(scores),
    indice=seq(1,length(scores)))

  data.table::setDT(d)
  data.table::setorder(d,-x)
  n_top_indice<-d$indice[1:n_top]
  return(n_top_indice)
}

#########################
# L1 penalty function
penalty_function_L1 <- function(scores, penalty_factor) {
  penalized_scores <- sapply(scores, function(x) {
    sign(x) * max(0, abs(x) - penalty_factor)
  })
  return(penalized_scores)
}
######################################################
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
######################################################
ssGED <- function(data, group_info, penalty_factor = 0.1) {
  normalized_data <- t(apply(data, 1, normalize))
  genexcell <- normalized_data
  n_sample <- ncol(normalized_data)
  n_gene <- nrow(normalized_data)
  gene_name <- rownames(genexcell)

  groups_order <- sort(unique(group_info$group))
  n_cluster <- length(groups_order)

  cluster_mat <- matrix(0, nrow = n_cluster, ncol = n_sample)
  order_i <- 1
  for (group_i in groups_order) {
    idx_i <- group_info == group_i
    cluster_mat[order_i, idx_i] <- 1
    order_i <- order_i + 1
  }

  genexcell <- Matrix(as.matrix(genexcell), sparse = TRUE)
  cosine_sim <- proxyC::simil(genexcell, cluster_mat, method = "cosine", drop0 = TRUE)

  # Apply L1 penalty to cosine similarity
  penalized_cosine_sim <- apply(cosine_sim, 2, function(column) penalty_function_L1(column, penalty_factor))

  # Initializes the resulting data frame
  rank_stats_names <- data.frame(matrix(NA, n_gene, length(groups_order),
                                        dimnames = list(seq(1, n_gene), groups_order)),
                                 stringsAsFactors = FALSE)
  rank_stats_scores <- data.frame(matrix(NA, n_gene, length(groups_order),
                                         dimnames = list(seq(1, n_gene), groups_order)),
                                  stringsAsFactors = FALSE)

  order_i <- 1
  for (group_i in groups_order) {
    idx_i <- group_info == group_i
    scores <- penalized_cosine_sim[, order_i]

    # The select_top_n function was used for scoring
    global_indices <- select_top_n(scores, n_gene)
    rank_stats_names[, order_i] <- gene_name[global_indices]
    rank_stats_scores[, order_i] <- scores[global_indices]

    order_i <- order_i + 1
  }

  colnames(rank_stats_names) <- groups_order
  colnames(rank_stats_scores) <- groups_order

  ranks_stats <- list(names = rank_stats_names,
                      scores = rank_stats_scores)


  first_5 <- apply(ranks_stats$name, 2, head, n = 20)
  merged_column <- c(first_5)
  TOPexp<- ncol(normalized_data)
  TOPgene <-normalized_data[merged_column, (1:TOPexp)]
  return(ranks_stats)
  }
######################################################
# Define a function that transforms the input matrix into a new S4 object
clGENE <- function(..., threshold = 0.5, export_close_genes = TRUE) {
  matrices <- list(...)

  # 初始化结果列表
  results <- list()
  close_genes_list <- list() # 新增一个列表来保存靠近基因对

  # 对每个矩阵执行操作
  for (i in seq_along(matrices)) {
    matrix <- matrices[[i]]
    normalized_data <- normalize(matrix)
    # 提取第一行
    first_row <- matrix[1, , drop = FALSE]
    column_df <- as.data.frame(t(first_row))
    # 从原矩阵中去掉第一行
    matrix <- matrix[-1, , drop = FALSE]
    ranks_stats <- ssGED(matrix, column_df)
    # 将第一行和剩余矩阵保存为一个列表
    first_20 <- apply(ranks_stats$name, 2, head, n = 200)
    merged_column <- c(first_20)
    TOPexp <- ncol(matrix)
    TOPgene <- normalized_data[merged_column, (1:TOPexp)]
    n_rows <- nrow(TOPgene)

    # 创建一个序列，表示每个子矩阵的起始行号
    row_indices <- seq(1, n_rows, 200)

    # 使用 lapply 和 split 函数分割矩阵
    list_of_matrices <- lapply(row_indices, function(idx) {
      # 计算子矩阵的结束行号
      end_row <- min(idx + 199, n_rows)

      # 提取子矩阵并返回
      TOPgene[idx:end_row, ]
    })

    list_of_pca_results <- lapply(list_of_matrices, function(sub_matrix) {
      # 转置矩阵，使行变为列
      transposed_matrix <- sub_matrix


      # 检查转置后的矩阵是否包含足够的数据进行PCA
      if(nrow(transposed_matrix) > 1) {
        #删除常数行列
        std_devs_cols <- apply(transposed_matrix, 2, sd)
        transposed_matrix <- transposed_matrix[, std_devs_cols != 0, drop = FALSE]
        std_devs_rows <- apply(transposed_matrix, 1, sd)
        transposed_matrix <- transposed_matrix[std_devs_rows != 0, , drop = FALSE]
        # 应用PCA
        pca_result <- prcomp(transposed_matrix, scale. = TRUE) # scale. = TRUE 表示数据标准化
        return(pca_result)
      } else {
        # 如果转置后只有一行，PCA无法应用
        return(NULL)
      }
    })
    results[[i]] <- list_of_pca_results
  }

  # 计算基因之间的相近性并导出结果
  for (i in seq_along(results)) {
    pca_results <- results[[i]]
    for (j in seq_along(pca_results)) {
      pca_result <- pca_results[[j]]
      if (!is.null(pca_result)) {
        # 提取坐标和基因名
        coords <- pca_result$x[, 1:2] # 只考虑前两个主成分
        genes <- rownames(pca_result$x)

        # 计算所有点对之间的欧氏距离
        dist_matrix <- as.matrix(dist(coords))

        # 寻找靠近的基因对
        close_genes <- which(dist_matrix < threshold, arr.ind = TRUE)

        # 创建一个数据框来存储靠近的基因对
        close_genes_df <- data.frame(
          gene1 = genes[close_genes[, 1]],
          gene2 = genes[close_genes[, 2]],
          distance = dist_matrix[close_genes]
        )

        # 过滤掉重复的组合和自身对比
        close_genes_df <- close_genes_df[close_genes_df$gene1 != close_genes_df$gene2, ]
        close_genes_df <- close_genes_df[!duplicated(t(apply(close_genes_df, 1, sort))), ]

        # 将结果添加到列表中
        close_genes_list[[paste("matrix", i, "comparison", j)]] <- close_genes_df

        # 如果需要导出靠近基因对
        if (export_close_genes) {
          write.csv(close_genes_df, sprintf("close_genes_matrix_%d_comparison_%d.csv", i, j), row.names = FALSE)
        }
      }
    }
  }

  # 返回结果
  return(list(results = results, close_genes_list = close_genes_list))
}

################################################
