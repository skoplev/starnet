nextflow.enable.dsl=2

params.out = 'OUT'
params.gt_gex_stats = 'IN/gt_gex_stats.tsv'

workflow {
  gt_gex_stats = Channel.fromPath(params.gt_gex_stats)

  db_effects_outcomes(gt_gex_stats)
  effects_exposures(gt_gex_stats, db_effects_outcomes.out)
  mr_twoSamples(db_effects_outcomes.out, effects_exposures.out)
}

process db_effects_outcomes {
  executor 'local'

  input: path('gt_gex_stats.tsv')
  output: path('effects_outcomes.RData')

  """
#!/usr/bin/env Rscript

library(tidyverse)
library(RPostgres)
library(dbplyr)

DB <- DBI::dbConnect(RPostgres::Postgres(),
   dbname = Sys.getenv('STATGEN_DB'),
   host = Sys.getenv('STATGEN_HOST'),
   port = as.numeric(Sys.getenv('STATGEN_PORT')),
   user = Sys.getenv('STATGEN_USER')
  )

message(" fetching variants")
variants <-
    read_tsv("gt_gex_stats.tsv",
             col_types = cols_only(gene_id = 'c', variant_id = 'c', p.value = 'd')) %>%
    rename(exposure = gene_id, pvalue = p.value) %>%
    filter(pvalue <= 0.001) %>%
    distinct()
str(variants)

message(" splitting in chunks")
## do this a few exposures at a time, otherwise the db might complain...
CHUNK_SIZE <- 100L
EXPOSURES <- sort(unique(variants\$exposure))
NUM_CHUNKS <- ceiling(length(EXPOSURES) / CHUNK_SIZE)
EXPOSURES_CHUNKS <- split(EXPOSURES, floor(seq(1, length(EXPOSURES)) %% NUM_CHUNKS))
str(EXPOSURES_CHUNKS)
VARIANTS_CHUNKS <- map(EXPOSURES_CHUNKS, ~ filter(variants, exposure %in% .))
str(EXPOSURES_CHUNKS)

message(" fetching from db")
get_effects_outcomes_chunk <- function(V) {
  DB %>%
    tbl(in_schema("stats", "stats_gwas")) %>%
    inner_join(V, by = c(variant_id = "variant_id"), copy = TRUE, suffix = c("", ".exposure")) %>%
    select(variant_id, exposure, study_id, trait_id, beta, beta_se, pvalue, pvalue.exposure) %>%
    collect() %>%
    group_by(exposure, study_id, trait_id) %>%
    arrange(pvalue.exposure) %>%
    filter(row_number() == 1) %>%
    select(-pvalue.exposure) %>%
    ungroup()
}

effects_outcomes <-
  map_dfr(VARIANTS_CHUNKS, get_effects_outcomes_chunk)

message(" saving results to disk")
save(effects_outcomes, file = "effects_outcomes.RData")

message("analysis completed.")
  """
}

process effects_exposures {
  executor 'local'

  input:
    path('gt_gex_stats.tsv')
    path('effects_outcomes.RData')
  output:
    path('effects_exposures.RData')

  """
  #!/usr/bin/env Rscript
library(tidyverse)

load("effects_outcomes.RData")

effects_exposures <-
    read_tsv("gt_gex_stats.tsv", col_types = cols_only(gene_id = 'c', variant_id = 'c', estimate = 'd', std.error = 'd', p.value = 'd')) %>%
    rename(exposure = gene_id, instrument = variant_id, beta = estimate, beta_se = std.error, pvalue = p.value) %>%
    semi_join(effects_outcomes, by = c(exposure = 'exposure', instrument = 'variant_id'))

save(effects_exposures, file = "effects_exposures.RData")
  """
}

process mr_twoSamples {
  
  executor 'local'

  publishDir params.out, mode: 'copy'

  input:
    path('effects_outcomes.RData')
    path('effects_exposures.RData')

  output:
    path('mr_twoSamples.tsv')

"""
#!/usr/bin/env Rscript

library(tidyverse)
library(MendelianRandomization)

load("effects_outcomes.RData")
load("effects_exposures.RData")

data <-
    effects_outcomes %>%
    inner_join(effects_exposures,
        by = c(variant_id = "instrument", exposure = "exposure"),
        suffix = c(".outcome", ".exposure")) %>%
    group_by(exposure, study_id, trait_id) %>%
    arrange(pvalue.exposure) %>%
    filter(row_number() == 1) %>%
    ungroup()

get_mr <- function(d) {
  tryCatch({
    with(d,
      mr_input(
        by = beta.outcome, byse = beta_se.outcome,
        bx = beta.exposure, bxse = beta_se.exposure
      ) %>%
      mr_allmethods(method = "ivw")
    )@Values %>%
    rename(pvalue = `P-value`)
  }, error = function(e) NULL)
}

mr_twoSamples <-
    data %>%
    group_by(exposure, study_id, trait_id) %>%
    nest() %>%
    mutate(fit = map(data, get_mr)) %>%
    select(-data) %>%
    unnest(fit) %>%
    filter(Method == "IVW") %>%
    select(-Method) %>%
    ungroup() %>%
    arrange(pvalue) %>%
    write_tsv("mr_twoSamples.tsv")
"""
}
