#' @title Compute Cellular Activity Integration Score (CAIScore)
#'
#' @description
#' This function calculates a composite score named **CAIScore** based on five widely-used gene set scoring algorithms
#' — AUCell, UCell, singscore, ssGSEA, and AddModuleScore — to assess biological signal intensity at the single-cell or sample level.
#' It supports both Seurat objects and gene expression matrices, and allows for flexible scaling and summarization strategies.
#'
#' This function supports both Seurat objects and standard gene expression matrices, providing a standardized, scalable, and reproducible framework for quantifying gene set activity.
#' By integrating multiple orthogonal scoring algorithms, it enables robust gene signature evaluation across diverse transcriptomic contexts
#' —particularly in single-cell or bulk RNA-seq analyses where consistent and interpretable quantification of biological pathways or functional modules is essential.
#'
#' @param expr A Seurat object (with `SCT` or `RNA` assay) or a numeric gene expression matrix with gene symbols as row names and samples/cells as column names.
#'
#' @param geneset A gene set or a list of gene sets used for scoring. Acceptable formats include:
#'   - A character vector (representing a single gene set),
#'   - A named list of character vectors (multiple gene sets).
#'
#' @param scaling Character string indicating the scaling method to apply to the individual algorithm scores before integration. Options:
#'   - `"zscore"`: Z-score normalization across cells/samples,
#'   - `"minmax"`: Rescale scores to the range [0, 1],
#'   - `"none"`: No scaling (raw scores used).
#'
#' @param summary_method Character string specifying how to aggregate scores from the five methods. Options:
#'   - `"mean"` (default),
#'   - `"sum"`,
#'   - `"median"`.
#'
#' @param auc_nCores Integer. Number of CPU cores to use in AUCell calculation. Default is 1.
#'
#' @param kcdf_type Character string specifying the kernel type for ssGSEA via GSVA. Options:
#'   - `"Gaussian"`: Suitable for normalized expression data,
#'   - `"Poisson"`: Recommended for raw count data.
#' @param AMS_nbin Integer, number of bins used in Seurat::AddModuleScore (default: 24). 
#' Useful to reduce this (e.g. to 10) when sample size is small to avoid cut_number errors.
#'
#' @return If `expr` is a Seurat object, returns the object with added metadata columns for each individual score and the final integrated `CAIScore`.
#' If `expr` is a matrix, returns a `data.frame` containing one row per sample/cell with scores from all five methods and the final integrated score.
#'
#' @details
#' The **CAIScore** is derived by integrating the outputs from five orthogonal gene signature scoring methods:
#'
#' - **AUCell**: Area Under the Curve-based method evaluating gene set enrichment in ranked expression profiles.
#' - **UCell**: A rank-based, normalization-agnostic signature scoring approach optimized for single-cell data.
#' - **singscore**: Uses relative ranks of gene expression to compute signature scores in an unbiased, reproducible manner.
#' - **ssGSEA** (via GSVA): Computes single-sample GSEA scores based on expression-level integration.
#' - **AddModuleScore**: Seurat-native module scoring function for evaluating predefined gene sets.
#'
#' The aggregated CAIScore enhances signal robustness by leveraging the strengths of each algorithm. It offers a unified score that better reflects biological activity across diverse transcriptomic contexts.
#'
#' @importFrom Seurat GetAssayData AddModuleScore AddMetaData
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @importFrom UCell ScoreSignatures_UCell
#' @importFrom GSVA gsva
#' @importFrom singscore simpleScore
#' @importFrom matrixStats rowMedians
#' @import Matrix
#'
#' @examples
#' # Example 1: Using a Seurat object
#' data("pbmc_small")
#' geneset <- list(Tcell_Signature = c("CD3D", "CD3E", "CD3G"))
#' pbmc_small <- CAIScore(
#'   expr = pbmc_small,
#'   geneset = geneset,
#'   scaling = "zscore",
#'   summary_method = "mean",
#'   auc_nCores = 2,
#'   kcdf_type = "Gaussian"
#' )
#'
#' # Example 2: Using an expression matrix
#' mat <- matrix(rnorm(1000), nrow = 100)
#' rownames(mat) <- paste0("Gene", 1:100)
#' colnames(mat) <- paste0("Sample", 1:10)
#' geneset <- list(CustomSet = paste0("Gene", 1:10))
#' score_result <- CAIScore(
#'   expr = mat,
#'   geneset = geneset,
#'   scaling = "minmax",
#'   summary_method = "sum",
#'   auc_nCores = 4,
#'   kcdf_type = "Gaussian"
#' )
#'
#' @author
#' Algorithm design: Ruiqiu Chen
#' Package development and integration: Ruiqiu Chen
#'
#' @seealso
#' \code{\link{AUCell}},
#' \code{\link{UCell}},
#' \code{\link{gsva}},
#' \code{\link{singscore}},
#' \code{\link[Seurat]{AddModuleScore}}
#'
#' @export
CAIScore <- function(expr,
                     geneset,
                     scaling = c("zscore", "minmax", "none"),
                     summary_method = c("mean", "sum", "median"),
                     auc_nCores = 1,
                     kcdf_type = c("Gaussian", "Poisson"),
                     AMS_nbin = 10 ) {

  ## ----------------------------
  ## Section 1: Argument validation and preprocessing
  ## ----------------------------
  if (missing(expr) || missing(geneset)) {
    stop("Both 'expr' and 'geneset' must be provided.")
  }

  scaling <- match.arg(scaling)
  summary_method <- match.arg(summary_method)
  kcdf_type <- match.arg(kcdf_type)

  # Convert vector-based gene set to list format
  if (!is.list(geneset)) {
    if (is.atomic(geneset)) {
      geneset <- list(Signature = unname(geneset))
    } else {
      stop("'geneset' must be a list or an atomic vector.")
    }
  }
  if (is.null(names(geneset))) {
    names(geneset) <- paste0("Signature_", seq_along(geneset))
  }

  ## ----------------------------
  ## Section 2: Handle input (Seurat or matrix)
  ## ----------------------------
  is_seurat <- inherits(expr, "Seurat")
  expr_mat <- if (is_seurat) {
    if ("SCT" %in% names(expr@assays)) {
      expr_mat <- Seurat::GetAssayData(expr, assay = "SCT", slot = "data")
    } else if ("RNA" %in% names(expr@assays)) {
      expr_mat <- Seurat::GetAssayData(expr, assay = "RNA", slot = "data")
    } else {
      stop("No valid assay (SCT or RNA) found in Seurat object.")
    }
  } else {
    expr <- as.matrix(expr)
    if (is.null(rownames(expr))) stop("Expression matrix must have rownames.")
    expr
  }

  ## ----------------------------
  ## Section 3: Gene set filtering
  ## ----------------------------
  geneset <- lapply(geneset, function(gs) {
    valid <- intersect(rownames(expr_mat), gs)
    if (length(valid) < 2) {
      warning("A gene set contains <2 valid genes and will be excluded.")
    }
    valid
  })
  geneset <- geneset[lengths(geneset) >= 2]
  if (length(geneset) == 0) stop("No valid gene sets after filtering.")

  sample_ids <- colnames(expr_mat)
  score_list <- list()

  ## ----------------------------
  ## Section 4: Scoring using multiple algorithms
  ## ----------------------------
  message("[CAIScore] Computing AUCell scores...")
  rankings <- AUCell::AUCell_buildRankings(expr_mat, nCores = auc_nCores, plotStats = FALSE)
  auc_res <- AUCell::AUCell_calcAUC(geneset, rankings, aucMaxRank = ceiling(0.2 * nrow(expr_mat)))
  auc_scores <- t(AUCell::getAUC(auc_res))
  colnames(auc_scores) <- paste0("AUCell_", names(geneset))
  score_list$AUCell <- rowMeans(auc_scores, na.rm = TRUE)

  message("[CAIScore] Computing UCell scores...")
  ucell_scores <- UCell::ScoreSignatures_UCell(matrix = expr_mat, features = geneset)
  colnames(ucell_scores) <- paste0("UCell_", names(geneset))
  score_list$UCell <- rowMeans(ucell_scores, na.rm = TRUE)

  message("[CAIScore] Computing singscore...")
  singscore_wrapper <- function(expr, gs) {
    ranked <- apply(expr, 2, rank)
    colMeans(ranked[gs, , drop = FALSE]) / nrow(expr)
  }
  singscore_mat <- sapply(geneset, function(gs) singscore_wrapper(expr_mat, gs))
  if (is.vector(singscore_mat)) singscore_mat <- matrix(singscore_mat, ncol = 1)
  colnames(singscore_mat) <- paste0("singscore_", names(geneset))
  score_list$singscore <- rowMeans(singscore_mat, na.rm = TRUE)

  message("[CAIScore] Computing ssGSEA scores...")
  ssgsea_mat <- GSVA::gsva(expr_mat, gset.idx.list = geneset, method = "ssgsea", kcdf = kcdf_type, verbose = FALSE)
  ssgsea_mat <- t(ssgsea_mat)
  colnames(ssgsea_mat) <- paste0("ssGSEA_", names(geneset))
  score_list$ssGSEA <- rowMeans(ssgsea_mat, na.rm = TRUE)

  message("[CAIScore] Computing AddModuleScore via Seurat...")
  if (is_seurat) {
    tmp_obj <- Seurat::AddModuleScore(expr, features = geneset, name = "AMS", nbin = AMS_nbin)
    ams_scores <- tmp_obj@meta.data[, grep("^AMS", colnames(tmp_obj@meta.data)), drop = FALSE]
    score_list$AddModuleScore <- rowMeans(ams_scores)
  } else {
    ams_mat <- sapply(geneset, function(gs) Matrix::colMeans(expr_mat[gs, , drop = FALSE]))
    if (is.vector(ams_mat)) ams_mat <- matrix(ams_mat, ncol = 1)
    colnames(ams_mat) <- paste0("AddModuleScore_", names(geneset))
    score_list$AddModuleScore <- rowMeans(ams_mat, na.rm = TRUE)
  }

  ## ----------------------------
  ## Section 5: Merge and scale scores
  ## ----------------------------
  merged_scores <- data.frame(
    AUCell = score_list$AUCell,
    UCell = score_list$UCell,
    singscore = score_list$singscore,
    ssGSEA = score_list$ssGSEA,
    AddModuleScore = score_list$AddModuleScore,
    row.names = sample_ids
  )

  # 标准化处理 merged_scores
  merged_scores_scaled <- switch(
    scaling,
    "zscore" = scale(merged_scores),
    "minmax" = apply(merged_scores, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (diff(rng) == 0) return(rep(0, length(x)))
      (x - rng[1]) / diff(rng)
    }),
    "none" = merged_scores
  )
  merged_scores_scaled <- as.data.frame(merged_scores_scaled)

  # 汇总计算最终 CAIScore
  merged_scores_scaled$CAIScore <- switch(
    summary_method,
    "mean" = rowMeans(merged_scores_scaled, na.rm = TRUE),
    "sum" = rowSums(merged_scores_scaled, na.rm = TRUE),
    "median" = matrixStats::rowMedians(as.matrix(merged_scores_scaled), na.rm = TRUE)
  )

  ## ----------------------------
  ## Section 6: Output
  ## ----------------------------
  if (is_seurat) {
    expr <- Seurat::AddMetaData(expr, merged_scores_scaled)
    return(expr)
  } else {
    expr_df <- as.data.frame(t(expr_mat))

    if (!identical(rownames(expr_df), rownames(merged_scores_scaled))) {
      if (!all(rownames(merged_scores_scaled) %in% rownames(expr_df))) {
        stop("Merged score samples not found in expression matrix.")
      }
      expr_df <- expr_df[rownames(merged_scores_scaled), , drop = FALSE]
    }

    score_result <- cbind(merged_scores_scaled, expr_df)
    score_cols <- colnames(merged_scores_scaled)
    score_result <- score_result[, c(score_cols, setdiff(colnames(score_result), score_cols))]
    return(score_result)
  }
}
