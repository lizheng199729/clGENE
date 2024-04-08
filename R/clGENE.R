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
clGENE <- function(...){
  matrices <- list(...)

  # Initializes the result list
  results <- list()
  grouped_pca_results <- list()
  close_genes_list <- list()
  # Perform an operation on each matrix
  for (i in seq_along(matrices)) {
    matrix <- matrices[[i]]
    normalized_data <- normalize(matrix)
    # Extract the first row
    first_row <- matrix[1, , drop = FALSE]
    column_df <- as.data.frame(t(first_row))
    # Remove the first row from the original matrix
    matrix <- matrix[-1, , drop = FALSE]
    ranks_stats<-ssGED(matrix, column_df)
    # Save the first row and the remaining matrix as a list
    first_20 <- apply(ranks_stats$name, 2, head, n = 20)
    merged_column <- c(first_20)
    TOPexp<- ncol(normalized_data)
    TOPgene <-normalized_data[merged_column, (1:TOPexp)]
    n_rows <- nrow(TOPgene)

    # Create a sequence representing the starting row number of each submatrix
    row_indices <- seq(1, n_rows, 20)

    # The matrix is split using the lapply and split functions
    list_of_matrices <- lapply(row_indices, function(i) {
      # Calculates the end row number of the submatrix
      end_row <- min(i + 19, n_rows)

      # Extract the submatrix and return it
      TOPgene[i:end_row, ]
    })
    list_of_pca_results <- lapply(list_of_matrices, function(sub_matrix) {
      # Transpose the matrix so that the rows become columns
      transposed_matrix <- t(sub_matrix)

      # The transposed matrix was checked to see if it contained sufficient data for PCA
      # At least two rows, two columns of the original matrix, are required to perform PCA
      if(nrow(transposed_matrix) > 1) {
        # Applying PCA
        pca_result <- prcomp(transposed_matrix, scale. = TRUE) # scale. = TRUE表示数据标准化
        return(pca_result)
      } else {
        # PCA cannot be applied if there is only one line after transposition
        return(NULL)
      }
    })
    results[[i]] <- list_of_pca_results
    grouped_pca_results <- list()

    for(i in seq_along(results)) {
      pca_results <- results[[i]]

      for(pca_result in pca_results) {
        if(!is.null(pca_result)) {
          # Extract the name of the current PCA result
          name <- rownames(pca_result$x)[1]

          # If the matrix of that name is not already grouped, initialize the grouping
          if(!name %in% names(grouped_pca_results)) {
            grouped_pca_results[[name]] <- list()
          }

          # The current PCA result is added to the corresponding grouping
          grouped_pca_results[[name]] <- c(grouped_pca_results[[name]], list(pca_result))
        }
      }
    }

    # Mapping analysis was performed for each group
    for(name in names(grouped_pca_results)) {
      pca_group <- grouped_pca_results[[name]]

      # Create a blank artboard to draw multiple PCA results
      plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "PC1", ylab = "PC2", main = paste("PCA for", name))

      for(pca_result in pca_group) {
        # Scatter plots were drawn for each PCA result
        points(pca_result$x[,1], pca_result$x[,2], col = rainbow(length(pca_result$x[,1])), pch = 16, cex = 0.6)

        # Add labels so that you know which list each point comes from
        text(pca_result$x[,1], pca_result$x[,2], labels = seq_along(pca_result$x[,1]), pos = 3, cex = 0.5)
      }
    }
    for(name in names(grouped_pca_results)) {
      pca_group <- grouped_pca_results[[name]]

      # Create a blank artboard to draw multiple PCA results
      plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "PC1", ylab = "PC2", main = paste("PCA for", name))

      for(pca_result in pca_group) {
        # Scatter plots were drawn for each PCA result
        points(pca_result$x[,1], pca_result$x[,2], col = rainbow(length(pca_result$x[,1])), pch = 16, cex = 0.6)
        # Add labels so that you know which list each point comes from
        text(pca_result$x[,1], pca_result$x[,2], labels = seq_along(pca_result$x[,1]), pos = 3, cex = 0.5)
      }

      # If needed to export close gene pairs
      if(export_close_genes) {
        # Coordinates and gene names were extracted
        coords <- pca_result$x[, 1:2] # Suppose we consider only the first two principal components
        genes <- rownames(coords)

        # Compute the Euclidean distance between all pairs of points
        dist_matrix <- as.matrix(dist(coords))

        # Initializes a data frame to store pairs of genes that are close together
        close_genes <- data.frame(gene1 = character(), gene2 = character(), distance = numeric())

        # The distance matrix was traversed to find close gene pairs
        for(j in 1:(nrow(dist_matrix)-1)) {
          for(k in (j+1):nrow(dist_matrix)) {
            if(dist_matrix[j, k] < threshold) {
              close_genes <- rbind(close_genes, data.frame(gene1 = genes[j], gene2 = genes[k], distance = dist_matrix[j, k]))
            }
          }
        }

        # Close gene pairs are added to the list
        close_genes_list[[paste(name, i)]] <- close_genes
      }
    }
  }

  # If needed to export close gene pairs
  if(export_close_genes) {
    # Close gene pairs for each grouping were derived
    for(name in names(close_genes_list)) {
      write.csv(close_genes_list[[name]], sprintf("%s_close_genes.csv", name), row.names = FALSE)
    }
  }

  return(results)
}


################################################
