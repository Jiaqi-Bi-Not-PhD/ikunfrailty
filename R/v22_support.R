## Minimal support retained from V2.2 simulation helpers.
## The package excludes data generation, CCA, SMCFCS, and result-writer code.

v22_check_popplus_support <- function(dat) {
  blocks <- v22_family_blocks(dat)
  checks <- lapply(blocks, function(idx) {
    p <- which(as.numeric(dat$proband[idx]) == 1)
    if (length(p) != 1L) return(FALSE)
    j <- idx[p]
    isTRUE(as.numeric(dat$mgene[j]) == 1) && isTRUE(as.numeric(dat$status[j]) == 1)
  })
  all(unlist(checks))
}
