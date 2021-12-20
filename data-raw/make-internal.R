### Internal datasets
library(devtools)
library(usethis)

### Comorbidity maps
source("data-raw/make-mapping.R")

### Comorbidity weights
source("data-raw/make-weights.R")

### Export
usethis::use_data(
  .maps,
  .weights,
  Elixhauser2020Formats,
  Elixhauser2021Formats,
  Elixhauser2022Formats,
  internal = TRUE,
  overwrite = TRUE
)
