# ---- 结点（knots）分位数布点 ----
make_knots <- function(x, k) {
  uq <- sort(unique(x))
  if (length(uq) <= k) return(uq)
  as.numeric(quantile(uq, probs = seq(0, 1, length.out = k)))
}

# ---- 安全 GAM 拟合 ----
safe_gam <- function(form, dat, fam, knots = NULL, method = "REML", select = FALSE) {
  try(
    mgcv::gam(
      formula = form, family = fam, data = dat,
      knots = knots, method = method, select = select
    ),
    silent = TRUE
  )
}

# ---- 稀疏/全零基因预过滤（一次性），并输出 message ----
# counts: G x N (genes x spots)
prefilter_genes_sum_nnz <- function(counts, min_nnz = 100, min_sum = 0) {
  nnz_per_gene <- Matrix::rowSums(counts > 0)
  sum_per_gene <- Matrix::rowSums(counts)

  keep <- which(nnz_per_gene >= min_nnz & sum_per_gene > min_sum)
  drop <- setdiff(seq_len(nrow(counts)), keep)

  # 统计：全零（sum=0）与过稀疏（sum>0 但 nnz<阈值）
  all_zero  <- which(sum_per_gene == 0)
  too_sparse_only <- setdiff(drop, all_zero)

  message(
    "Prefilter summary: drop=", length(drop),
    " (all-zero=", length(all_zero), ", too-sparse=", length(too_sparse_only), "). ",
    "Keep=", length(keep), " genes."
  )
  keep
}

# ---- 只按“总表达量”选择 Top N 基因（HEG / Top expressed）----
select_top_expressed_genes <- function(counts, genes, top_n = 1000) {
  gene_sum <- Matrix::rowSums(counts)              # 每个基因 across 所有 spots 的总 counts
  ord <- order(gene_sum, decreasing = TRUE)
  keep <- head(ord, min(top_n, length(ord)))
  message("Selected top ", length(keep), " expressed genes by total counts.")
  list(counts = counts[keep, , drop = FALSE], genes = genes[keep])
}

# ---- 单基因拟合（支持 knots + 稀疏跳过）----
fit_one_gene <- function(j, counts, genes, covar, distance_vars,
                         fam_nb, use_interaction = FALSE, k = 6L,
                         min_nnz_gene = 1L, knots_list = NULL) {
  # 只取一行避免稠密化整块
  y <- as.numeric(counts[j, , drop = FALSE])

  # 稀疏基因过滤（行内保险）
  if (sum(y > 0) < min_nnz_gene || sum(y) == 0) {
    return(c(gene = genes[j], status = "too_sparse",
             dev_expl = NA, p_int = NA))
  }

  dat <- data.frame(
    y   = y,
    off = covar$log_off,
    covar[distance_vars],
    check.names = FALSE
  )
  # 标准化距离 + 守护 offset
  dat[distance_vars] <- lapply(dat[distance_vars], function(x) as.numeric(scale(x)))
  dat$off[!is.finite(dat$off)] <- min(dat$off[is.finite(dat$off)], na.rm = TRUE)
  ok <- is.finite(dat$off)
  for (v in distance_vars) ok <- ok & is.finite(dat[[v]])
  dat <- dat[ok, , drop = FALSE]

  # baseline
  m0 <- safe_gam(y ~ offset(off), dat, fam = fam_nb, method = "REML")

  # additive
  smooth_terms <- paste0("s(", distance_vars, ", k=", k, ", bs='cr')", collapse = " + ")
  form_add <- as.formula(paste("y ~ offset(off) +", smooth_terms))
  m_add <- safe_gam(form_add, dat, fam = fam_nb, method = "REML", knots = knots_list)

  if (inherits(m_add, "try-error") || inherits(m0, "try-error")) {
    return(c(gene = genes[j], status = "fit_failed", dev_expl = NA, p_int = NA))
  }

  st <- summary(m_add)$s.table
  p_vals <- rep(NA_real_, length(distance_vars))
  for (i in seq_along(distance_vars)) {
    ix <- grep(distance_vars[i], rownames(st))
    if (length(ix)) p_vals[i] <- st[ix[1], "p-value"]
  }
  names(p_vals) <- paste0("p_", distance_vars)

  dev_expl <- 1 - deviance(m_add) / deviance(m0)

  # 可选交互
  p_int <- NA_real_
  if (use_interaction && length(distance_vars) >= 2) {
    pairs <- combn(distance_vars, 2)
    ti_terms <- apply(pairs, 2, function(v) sprintf("ti(%s,%s,k=%d,bs=c('cr','cr'))", v[1], v[2], k))
    form_full <- as.formula(paste("y ~ offset(off) +", smooth_terms, "+", paste(ti_terms, collapse = " + ")))
    m_full <- safe_gam(form_full, dat, fam = fam_nb, method = "ML", knots = knots_list)
    m_add_ml <- try(update(m_add, method = "ML"), silent = TRUE)
    if (!inherits(m_full, "try-error") && !inherits(m_add_ml, "try-error")) {
      p_int <- tryCatch(anova(m_add_ml, m_full, test = "Chisq")[2, "Pr(>Chi)"], error = function(e) NA_real_)
    }
  }

  c(gene = genes[j], status = "ok", dev_expl = dev_expl, p_int = p_int, p_vals)
}

# ---- 批量分析：支持 HEG 子集 + 预过滤 + 结点 + 信息打印 ----
analyze_spatial_genes <- function(counts, genes, covar, distance_vars,
                                  fam_nb = mgcv::nb(link = "log"),
                                  use_interaction = FALSE, k = 6L,
                                  min_nnz = 100, min_sum = 0,
                                  min_nnz_gene = 1L,
                                  top_n_expressed = NULL) {
  G_all <- nrow(counts)

  # 仅跑 Top expressed（可选）
  if (!is.null(top_n_expressed)) {
    sel <- select_top_expressed_genes(counts, genes, top_n = top_n_expressed)
    counts <- sel$counts
    genes  <- sel$genes
  }

  # 预过滤（一次性）
  keep <- prefilter_genes_sum_nnz(counts, min_nnz = min_nnz, min_sum = min_sum)
  if (length(keep) == 0) {
    message("No genes left after prefilter. Stopping.")
    return(data.frame(gene = character(), status = character()))
  }

  # 生成 knots（每个 dist_* 一个结点向量）
  knots_list <- lapply(setNames(distance_vars, distance_vars), function(v) {
    make_knots(covar[[v]], k)
  })

  # 统计跳过/失败
  skipped_sparse <- 0L
  failed_fit     <- 0L

  res_list <- pbapply::pblapply(keep, function(j) {
    out <- fit_one_gene(
      j, counts, genes, covar, distance_vars, fam_nb,
      use_interaction = use_interaction, k = k,
      min_nnz_gene = min_nnz_gene, knots_list = knots_list
    )
    status <- as.character(out["status"])
    if (identical(status, "too_sparse")) skipped_sparse <<- skipped_sparse + 1L
    if (identical(status, "fit_failed")) failed_fit   <<- failed_fit + 1L
    out
  })

  message("Skip (all-zero / too-sparse in per-gene check): ", skipped_sparse, " genes")
  message("Fit failed (try-error): ", failed_fit, " genes")
  message("Analyze done. Input genes: ", G_all,
          ", kept for modeling: ", length(keep),
          ", returned rows: ", length(res_list), ".")

  res <- as.data.frame(do.call(rbind, res_list), stringsAsFactors = FALSE)
  num_cols <- c("dev_expl", "p_int", paste0("p_", distance_vars))
  for (nm in intersect(num_cols, names(res))) res[[nm]] <- as.numeric(res[[nm]])
  res
}

# ---- 计算 log offset（counts 与 covar$barcode 对齐）----
compute_log_offset <- function(counts, barcodes, pseudocount = 1e-8, normalize = TRUE) {
  if (!is.matrix(counts) && !"dgCMatrix" %in% class(counts)) {
    counts <- as.matrix(counts)
  }
  ord <- match(barcodes, colnames(counts))
  if (anyNA(ord)) stop("Some barcodes in covar not found in counts colnames.")
  counts_aln <- counts[, ord, drop = FALSE]

  libsize <- Matrix::colSums(counts_aln)
  if (normalize) {
    sf <- as.numeric(libsize / stats::median(libsize[libsize > 0]))
  } else {
    sf <- as.numeric(libsize)
  }
  log_off <- log(pmax(sf, pseudocount))
  list(size_factor = sf, log_off = log_off, counts_aligned = counts_aln)
}


# exportToolsUI_global <- function(id) {
#   ns <- NS(id)
#   tagList(
#     hr(), h4("Export & tools"),
#     fluidRow(
#       column(3, actionButton(ns("clear_sel"), "Clear selection")),
#       column(3, downloadButton(ns("dl_sel"), "Download selected CSV")),
#       column(3, downloadButton(ns("dl_multi_distance"), "Download all structure distances")),
#       column(3, downloadButton(ns("dl_gam_results"), "Download GAM results"))  # ← 合并进来
#     )
#   )
# }


# -- Server 模块 --
gammmServer <- function(id, rv, multi_distance_df_r, log1p_r = reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {

    `%||%` <- function(a,b) if (is.null(a) || length(a)==0 || (is.atomic(a) && !isTRUE(is.finite(a)))) b else a

    # 兜底：如果项目里没有 compute_log_offset，就临时内置一个
    .compute_log_offset <- function(counts, barcodes, normalize = FALSE) {
      ord <- match(barcodes, colnames(counts))
      ok  <- !is.na(ord)
      counts_aligned <- counts[, ord[ok], drop = FALSE]
      sf <- rep(1, ncol(counts_aligned))
      if (isTRUE(normalize)) {
        lib <- Matrix::colSums(counts_aligned)
        med <- stats::median(lib[lib > 0])
        if (is.finite(med) && med > 0) sf <- as.numeric(lib / med)
      }
      list(
        counts_aligned = counts_aligned,
        size_factor    = sf,
        log_off        = log(pmax(sf, 1e-8)),
        ok_mask        = ok
      )
    }

    # 把 analyze_spatial_genes 作为依赖函数检查一下
    if (!is.function(get0("analyze_spatial_genes", mode = "function"))) {
      stop("需要函数 analyze_spatial_genes()，请确保已 source 进来。")
    }

    output$dl_gam_results <- downloadHandler(
      filename = function() paste0("gam_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) {
        withProgress(message = "Running GAM fitting ...", value = 0, {
          req(rv$counts, rv$gene_names, rv$pos)
          covar <- as.data.frame(multi_distance_df_r())   # <- 从外部 reactive 拿 covariates
          stopifnot("barcode" %in% names(covar))

          # 兼容你项目里的 compute_log_offset；没有就用兜底版
          off_fun <- get0("compute_log_offset", mode = "function") %||% .compute_log_offset
          off <- off_fun(rv$counts, covar$barcode, normalize = isTRUE(log1p_r()))
          covar$size_factor   <- off$size_factor
          covar$log_off       <- off$log_off
          counts_aligned      <- off$counts_aligned

          # 取 dist_* 作为平滑变量
          distance_vars <- grep("^dist_", names(covar), value = TRUE)
          if (length(distance_vars) == 0) stop("No dist_* columns found, 请先生成距离。")

          # 参数：若外层 UI 没提供，就用默认值
          k_val       <- (get0("input", ifnotfound = NULL)$svg_k %||% 6L)      # 你也可改成固定 6L
          min_nnz_val <- (get0("input", ifnotfound = NULL)$svg_min_nnz %||% 100L)
          topN_val    <- (get0("input", ifnotfound = NULL)$svg_topN %||% 1000L)

          res <- analyze_spatial_genes(
            counts            = counts_aligned,
            genes             = rv$gene_names,
            covar             = covar,
            distance_vars     = distance_vars,
            fam_nb            = mgcv::nb(link = "log"),
            use_interaction   = FALSE,
            k                 = as.integer(k_val),
            min_nnz           = as.integer(min_nnz_val),
            min_sum           = 0,
            min_nnz_gene      = 1L,
            top_n_expressed   = as.integer(topN_val)
          )

          utils::write.csv(res, file, row.names = FALSE)
        })
      }
    )

      gam_df <- reactive({
        req(input$gam_csv)  # 注意：只有当 fileInput(ns("gam_csv")) 在 gamUI 里时，这里才看得到
        df <- tryCatch(
          utils::read.csv(input$gam_csv$datapath, stringsAsFactors = FALSE, check.names = FALSE),
          error = function(e) NULL
        )
        req(!is.null(df))
        nm <- trimws(names(df)); nm <- sub("^ï\\.+", "", nm); names(df) <- nm
        if (!"gene" %in% names(df)) stop("CSV missing 'gene' column")
        df
      })
    return(list(gam_df = gam_df))
  })
}
