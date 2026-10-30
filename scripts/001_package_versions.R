## 001_package_versions.R -- inventory of R packages actually used by the pipeline
## writes docs/package_versions.tsv and docs/package_versions.md
source("E:/another data file/GEObulk analysis project/scripts/00_setup.R")
use("data.table")

sc <- list.files(DIR$scripts, pattern = "\\.(R|r)$", full.names = TRUE, recursive = FALSE)
sc <- sc[!grepl("_dev", sc, fixed = TRUE)]
txt <- vapply(sc, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1))

pat_lib <- "(?:library|require)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)"
pat_ns  <- "([A-Za-z][A-Za-z0-9._]*):::?"
hits_lib <- unlist(regmatches(txt, gregexpr(pat_lib, txt, perl = TRUE)))
hits_lib <- unique(sub(pat_lib, "\\1", hits_lib, perl = TRUE))
hits_ns  <- unlist(regmatches(txt, gregexpr(pat_ns, txt, perl = TRUE)))
hits_ns  <- unique(sub(pat_ns, "\\1", hits_ns, perl = TRUE))

base_pkgs <- rownames(installed.packages(priority = "base"))
cand <- unique(c(hits_lib, hits_ns))
cand <- cand[nzchar(cand)]
ip <- installed.packages()[, c("Package", "Version", "LibPath")]
cand <- intersect(cand, ip[, "Package"])

ver <- data.table(package = cand,
                  version = ip[cand, "Version"],
                  lib = ip[cand, "LibPath"],
                  base = cand %in% base_pkgs)
setorder(ver, -base, package)

d <- new_dir(DIR$root, "docs")
fwrite(ver, file.path(d, "package_versions.tsv"), sep = "\t")

md <- c("# R package versions used by the GEObulk heart-failure pipeline",
        "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               " | R ", getRversion(), " | platform ", R.version$platform),
        "",
        "This table is a lightweight substitute for `renv.lock`. `renv` is not",
        "installed on the analysis machine and CRAN access was unavailable, so",
        "version pinning is recorded here instead. To convert it into a real",
        "lockfile on a networked machine: install `renv`, then run",
        "`renv::init(bare = TRUE)` followed by `renv::snapshot()`.",
        "",
        paste0("Total packages referenced by the scripts: ", nrow(ver)),
        "",
        "| package | version | base |",
        "| --- | --- | --- |",
        paste0("| ", ver$package, " | ", ver$version, " | ", ifelse(ver$base, "yes", ""), " |"))
writeLines(md, file.path(d, "package_versions.md"))

cat("scanned scripts:", length(sc), "\n")
cat("packages written:", nrow(ver), "\n")
print(head(ver, 20))
write_session(file.path(DIR$results, "_setup_test"), "001_package_versions")
