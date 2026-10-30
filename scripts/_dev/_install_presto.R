options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
lib <- "D:/R/library/4.6"
res <- tryCatch({
  install.packages("presto", lib = lib)
  "install-call-done"
}, error = function(e) paste("ERR", conditionMessage(e)),
   warning = function(w) paste("WARN", conditionMessage(w)))
cat(res, "\n")
ok <- requireNamespace("presto", quietly = TRUE)
cat("presto available:", ok, "\n")
cat("libpaths:", paste(.libPaths(), collapse = ";"), "\n")
