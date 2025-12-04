#' 日志信息输出函数
#'
#' @param type 日志级别，可选"INFO", "WARNING", "ERROR"
#' @param msg 日志内容（必填）
#' @return 无返回值，直接输出日志到控制台
#' @export
log_message <- function(type = c("INFO", "WARNING", "ERROR"), msg = NULL) {
  if (is.null(msg)) {
    stop("参数 'msg' 不能为空，请提供日志内容！")
  }
  type <- match.arg(type)
  message(sprintf("[%s] %s", type, msg))
}

#' dmr分析统计主函数
#'
#' @param output 输出文件夹路径
#' @param windows 检测窗口长度向量
#' @param compare 分组比较矩阵，每列为一组比较
#' @param BSobj DSS::BSseq对象
#' @param delta DML/DMR差异阈值（默认0.1）
#' @param p.threshold 显著性阈值（默认0.05）
#' @return 命名list，每个元素为对应比对条件和窗口的dmr数据框
#' @export
analyze_data <- function(output, windows, compare, BSobj, delta = 0.1, p.threshold = 0.05) {
  # 检查输出目录
  if (!dir.exists(output)) dir.create(output, recursive = TRUE, showWarnings = FALSE)

  results <- data.frame(
    name = character(),
    delta = numeric(),
    p.threshold = numeric(),
    window = numeric(),
    dml_rows = integer(),
    dmr_rows = integer(),
    stringsAsFactors = FALSE
  )

  dmr_all <- list()  # 多window/比对统一收集

  for (window in windows) {
    for (i in seq_len(ncol(compare))) {
      comparison <- compare[, i]
      comp_name <- paste0(comparison[1], "_vs_", comparison[2])
      log_message("INFO", sprintf("processing: %s; window: %s", comp_name, window))

      # 统计检验，对WGBS数据进行smooth
      BS_test <- DSS::DMLtest(BSobj, group1 = comparison[1], group2 = comparison[2], smoothing = TRUE)
      dmls <- DSS::callDML(BS_test, delta = delta, p.threshold = p.threshold)
      dmrs <- DSS::callDMR(BS_test, delta = delta, p.threshold = p.threshold, minlen = window)

      # 输出dml文件
      dml_file <- file.path(output, sprintf("%s_dmls_%s.bed", comp_name, window))
      utils::write.table(dmls, dml_file, quote = FALSE, row.names = FALSE)
      log_message("INFO", sprintf("保存dml文件为：%s", dml_file))

      # 检查dmr结果并处理
      if (is.null(dmrs) || (is.data.frame(dmrs) && nrow(dmrs) == 0)) {
        dmrs <- data.frame()
        log_message("WARNING", sprintf("找不到dmr [比对: %s; window: %s]", comp_name, window))
      }

      # 导出dmr文件，空结果也输出
      dmr_file <- file.path(output, sprintf("%s_dmrs_%s.bed", comp_name, window))
      utils::write.table(dmrs, dmr_file, quote = FALSE, row.names = FALSE)
      log_message("INFO", sprintf("保存dmr文件为：%s", dmr_file))

      # 保存结果至list，key为比对+窗口
      list_key <- paste(comp_name, window, sep = "_win")
      dmr_all[[list_key]] <- dmrs

      # 统计行数并添加到结果表
      results <- dplyr::bind_rows(results, data.frame(
        name = comp_name,
        delta = delta,
        p.threshold = p.threshold,
        window = window,
        dml_rows = if (!is.null(dmls) && is.data.frame(dmls)) nrow(dmls) else 0,
        dmr_rows = if (!is.null(dmrs) && is.data.frame(dmrs)) nrow(dmrs) else 0,
        stringsAsFactors = FALSE
      ))
    }
  }

  # 保存统计汇总
  summary_file <- file.path(output, "summary.csv")
  utils::write.csv(results, summary_file, row.names = FALSE, quote = FALSE)
  log_message("INFO", sprintf("保存统计汇总为：%s", summary_file))

  return(dmr_all)
}