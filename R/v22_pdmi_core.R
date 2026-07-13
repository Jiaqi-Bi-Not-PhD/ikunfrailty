## ============================================================
## SMCFCS-style public PDMI interface for frailtypack-congenial MI.
## May 27 exact-slice technical note: posterior parameter draws plus Rubin pooling.
## ============================================================

v22_pdmi_frailtypack <- function(dat, K,
                                 missing_type = c("continuous", "binary", "joint"),
                                 prior_version = c("C-O", "C-R", "B-O", "B-R", "J-O", "J-R"),
                                 M = 20L,
                                 numit = NULL,
                                 config = v22_default_config(),
                                 seed = NULL,
                                 pedigree_dat = NULL) {
  missing_type <- match.arg(missing_type)
  prior_version <- match.arg(prior_version)
  cfg <- config
  cfg$M_imp_pdmi <- as.integer(M)
  cfg$M_imp_cong <- as.integer(M)
  cfg$pdmi_numit <- as.integer(numit %||% cfg$pdmi_numit %||% 10L)

  if (missing_type == "continuous") {
    if (!prior_version %in% c("C-O", "C-R")) {
      stop("Continuous PDMI prior_version must be C-O or C-R.")
    }
    return(v22_run_continuous_pdmi(dat, K, prior_version, cfg, seed))
  }

  if (missing_type == "joint") {
    if (!prior_version %in% c("J-O", "J-R")) {
      stop("Joint PDMI prior_version must be J-O or J-R.")
    }
    return(v22_run_joint_pdmi(dat, K, prior_version, pedigree_dat, cfg, seed))
  }

  if (!prior_version %in% c("B-O", "B-R")) {
    stop("Binary PDMI prior_version must be B-O or B-R.")
  }
  v22_run_binary_pdmi(dat, K, prior_version, pedigree_dat, cfg, seed)
}
