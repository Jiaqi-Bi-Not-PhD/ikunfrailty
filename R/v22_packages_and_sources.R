## ============================================================
## Package and legacy-source loading for V2.2.
## Legacy files are loaded only for family skeleton/FamEvent support and
## kinship mechanics. V2.2 files override the missing-data and pooling methods.
## ============================================================

v22_set_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    BLAS_NUM_THREADS = "1",
    LAPACK_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  root <- if (exists("v22_code_root", inherits = TRUE)) {
    get("v22_code_root", inherits = TRUE)
  } else {
    getwd()
  }
  lib_candidates <- unique(c(
    file.path(root, "r_libs_4.5"),
    file.path(getwd(), "r_libs_4.5"),
    file.path(getwd(), "Code V2.2 - frailtypack PDMI congenial May 27; exact-slice", "r_libs_4.5"),
    file.path(dirname(root), "r_libs_4.5"),
    file.path(dirname(root), "Code V2.2 - frailtypack PDMI congenial May 27; exact-slice", "r_libs_4.5")
  ))
  lib_candidates <- lib_candidates[dir.exists(lib_candidates)]
  if (length(lib_candidates)) {
    old <- Sys.getenv("R_LIBS_USER", "")
    Sys.setenv(R_LIBS_USER = paste(c(lib_candidates, old[nzchar(old)]), collapse = .Platform$path.sep))
    .libPaths(unique(c(lib_candidates, .libPaths())))
  }
  invisible(TRUE)
}

v22_require_packages <- function(include_frailtypack = TRUE, include_smcfcs = FALSE) {
  if (isTRUE(include_smcfcs)) {
    stop("SMCFCS comparators are not included in ikunfrailty; only proposed PDMI is supported.")
  }
  required <- c("survival")
  if (isTRUE(include_frailtypack)) required <- c(required, "frailtypack")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required package(s) not installed: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

v22_source_if_exists <- function(path) {
  if (file.exists(path)) {
    source(path)
    return(TRUE)
  }
  FALSE
}

v22_load_family_sources <- function(code_root = v22_find_code_root(getwd())) {
  dep_dir <- file.path(code_root, "dependencies")
  rels <- c(
    "Delete Males.R",
    "FamEvent/R/cumhaz.R",
    "FamEvent/R/hazards.R",
    "FamEvent/R/gh.R",
    "FamEvent/R/penmodel.R",
    "FamEvent/R/loglik_frailty.R",
    "FamEvent/R/dlaplace.R",
    "FamEvent/R/laplace.R",
    "FamEvent/R/familyDesign.R",
    "FamEvent/R/fgeneZX.R",
    "FamEvent/R/Pgene.R",
    "FamEvent/R/surv.dist.R",
    "FamEvent/R/survp.dist.R",
    "FamEvent/R/inv.surv.R",
    "FamEvent/R/inv2.surv.R",
    "FamEvent/R/parents.g.R",
    "FamEvent/R/kids.g.R",
    "FamEvent/R/simfam.R",
    "famevent_namespace_fallbacks.R",
    "familyStructure_1to20.R"
  )
  loaded <- vapply(file.path(dep_dir, rels), v22_source_if_exists, logical(1))
  invisible(loaded)
}
