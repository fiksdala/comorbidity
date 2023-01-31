### Downloads and parses the 2023 formatting sas program

destfile_path = "AHRQ-Elixhauser/sas-parse/icd10cm_2023_1/CMR_v2023-1.zip"

# Download file
download.file(
  url = "https://www.hcup-us.ahrq.gov/toolssoftware/comorbidityicd10/CMR_v2023-1.zip",
  destfile = destfile_path
)
# Unzip
unzip(destfile_path,
      exdir = 'AHRQ-Elixhauser/sas-parse/icd10cm_2023_1/ElixhauserComorbidity_v2023-1')

# Get raw SAS code line-by-line
raw_format = readLines(
  "AHRQ-Elixhauser/sas-parse/icd10cm_2023_1/ElixhauserComorbidity_v2023-1/CMR_Format_Program_v2023-1.sas"
)

# Remove quotes, commas, and whitespace
trim_format = trimws(gsub(',', '', gsub('"', "", raw_format)))
# Remove 'proc' lines
trim_format = trim_format[!grepl('Proc', trim_format)]
# Remove 'run' lines
trim_format = trim_format[!grepl(';', trim_format)]
# Remove 'other' lines
trim_format = trim_format[!grepl('other', trim_format)]

# Separate vector by blank line
format_list = split(
  trim_format[trim_format!=''],
  cumsum(trim_format=="")[trim_format!='']
)

# Remove extraneous
format_list = format_list[3:length(format_list)] # Header stuff

# Get comfmt vs. poaxmpt entries
poaxmpt_filter <- unlist(
  lapply(format_list, function(x){grepl(' = 1', tail(x, 1))})
)

# Get names for comfmt
comfmt_names <- as.vector(
  unlist(
    lapply(
      format_list[!poaxmpt_filter],
      function(x){
        strsplit(tail(x, 1), ' = ')[[1]][2]
      }
    )
  )
)
# Get values for each group in comfmt
comfmt_values <- list()
# First element is unique as it contains comfmt header
comfmt_values <- append(
  comfmt_values,
  list(
    sapply(
      strsplit(format_list[!poaxmpt_filter][[1]][-1], ' = '),
      function(x){
        x[[1]]
      }
    )
  )
)

# Subsequent elements follow the same pattern
comfmt_values <- append(
  comfmt_values,
  lapply(
    format_list[!poaxmpt_filter][-1],
    function(x){
      # Splits the last element that contains " = NAME"
      sapply(
        # Takes the first element of the previous split
        strsplit(x, ' = '),
        function(y){
          y[[1]]
        }
      )
    }
  )
)

# Add names to values
names(comfmt_values) <- comfmt_names

# Get poaxmpt names
poaxmpt_names <- as.vector(
  sapply(
    format_list[poaxmpt_filter],
    function(x){
      strsplit(x[1], '\\$')[[1]][2]
    }
  )
)

# Get poaxmpt values
poaxmpt_values <- lapply(
  format_list[poaxmpt_filter],
  function(x){
    gsub(' = 1', '', x[-1])
  }
)

# Add names to values
names(poaxmpt_values) <- tolower(poaxmpt_names)

# Create complete ElixhauserAHRQ2023Map
ElixhauserAHRQ2023Map <- list(
  comfmt = comfmt_values
)
ElixhauserAHRQ2023Map <- append(
  ElixhauserAHRQ2023Map,
  poaxmpt_values
)

ElixhauserAHRQ2023PreExclusion = c(
  'AIDS',
  'ALCOHOL',
  'ANEMDEF',
  'AUTOIMMUNE',
  'BLDLOSS',
  'CANCER_LYMPH',
  'CANCER_LEUK',
  'CANCER_METS',
  'CANCER_NSITU',
  'CANCER_SOLID',
  'CBVD_SQLA',
  'CBVD_POA',
  'CBVD_NPOA',
  'CBVD',
  'HF',
  'COAG',
  'DEMENTIA',
  'DEPRESS',
  'DIAB_UNCX',
  'DIAB_CX',
  'DRUG_ABUSE',
  'HTN_CX',
  'HTN_UNCX',
  'LIVER_MLD',
  'LIVER_SEV',
  'LUNG_CHRONIC',
  'NEURO_MOVT',
  'NEURO_OTH',
  'NEURO_SEIZ',
  'OBESE',
  'PARALYSIS',
  'PERIVASC',
  'PSYCHOSES',
  'PULMCIRC',
  'RENLFL_MOD',
  'RENLFL_SEV',
  'THYROID_HYPO',
  'THYROID_OTH',
  'ULCER_PEPTIC',
  'VALVE',
  'WGHTLOSS'
)

# Define and save final comorbidities in AHRQ format:
ElixhauserAHRQ2023Abbr <- c(
  'AIDS',
  'ALCOHOL',
  'ANEMDEF',
  'AUTOIMMUNE',
  'BLDLOSS',
  'CANCER_LYMPH',
  'CANCER_LEUK',
  'CANCER_METS',
  'CANCER_NSITU',
  'CANCER_SOLID',
  'CBVD',
  'HF',
  'COAG',
  'DEMENTIA',
  'DEPRESS',
  'DIAB_UNCX',
  'DIAB_CX',
  'DRUG_ABUSE',
  'HTN_CX',
  'HTN_UNCX',
  'LIVER_MLD',
  'LIVER_SEV',
  'LUNG_CHRONIC',
  'NEURO_MOVT',
  'NEURO_OTH',
  'NEURO_SEIZ',
  'OBESE',
  'PARALYSIS',
  'PERIVASC',
  'PSYCHOSES',
  'PULMCIRC',
  'RENLFL_MOD',
  'RENLFL_SEV',
  'THYROID_HYPO',
  'THYROID_OTH',
  'ULCER_PEPTIC',
  'VALVE',
  'WGHTLOSS'
)

# Save list of format objects
Elixhauser2023Formats = list(
  ElixhauserAHRQ2023Map = ElixhauserAHRQ2023Map,
  ElixhauserAHRQ2023Abbr = ElixhauserAHRQ2023Abbr,
  ElixhauserAHRQ2023PreExclusion = ElixhauserAHRQ2023PreExclusion
)

# Remove .zip file
file.remove(
  destfile_path
)

# Remove unzipped folder
unlink(
  'AHRQ-Elixhauser/sas-parse/icd10cm_2023_1/ElixhauserComorbidity_v2023-1',
  recursive = T
)
