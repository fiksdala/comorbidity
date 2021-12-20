library(stringr)
library(data.table)

download.file(
  url = "https://www.hcup-us.AHRQ.gov/toolssoftware/comorbidityicd10/comformat_icd10cm_2020_1.txt",
  destfile = "data-raw/icd10cm_2020_1/comformat_icd10cm_2020_1.txt"
)


# Example pattern to extract:
# "D500",
# "O9081",
# "O99011",
# "O99012",
# "O99013",
# "O99019",
# "O9902",
# "O9903"="BLDLOSS"        /*Blood loss anemia*/

# Read in data
sas_path <- "data-raw/icd10cm_2020_1/comformat_icd10cm_2020_1.txt"
sas_AHRQ_raw <- readLines(sas_path)

# Make list of lists for AHRQ codes to compare
make_sas_list = function(sas_AHRQ_raw){
  # Assigns ICD-10 codes to comorbidty labels from sas file located here:
  # https://www.hcup-us.AHRQ.gov/toolssoftware/comorbidityicd10/comformat_icd10cm_2020_1.txt
  # Omits /**** ICD-10 MS-DRG V37 Formats ****/

  # Clean up readlines
  sas_AHRQ_prep <- sas_AHRQ_raw[sas_AHRQ_raw!=""] # Remove empty lines
  sas_AHRQ_prep <- sas_AHRQ_prep[-(1:18)] # First 18 elements are extraneous
  sas_AHRQ_prep <- lapply(
    sas_AHRQ_prep,
    function(x) str_split(x,"\\/\\*")[[1]][1]
  ) # Drop sas comments
  sas_AHRQ_prep <- unlist(sas_AHRQ_prep)
  sas_AHRQ_prep <- str_trim(sas_AHRQ_prep) # Trim white space
  sas_AHRQ_prep <- str_replace_all(sas_AHRQ_prep, '\\"', "") # Remove extraneous characters
  sas_AHRQ_prep <- str_replace_all(sas_AHRQ_prep, ',', "") # Remove extraneous characters

  AHRQ_list = list() # create empty list
  temp_list = c() # placeholder for codes
  for(l in sas_AHRQ_prep){
    if(grepl("=", l, fixed=T)){
      split_l = str_split(l,'=')[[1]]
      temp_list = append(temp_list, split_l[1])
      # AHRQ_list[[split_l[2]]] = str_c(temp_list, collapse="|")
      # Must have ^ so that regex doesn't search for within-code substrings
      AHRQ_list[[split_l[2]]] = paste0("^", str_c(temp_list, collapse="|^"))
      temp_list = c()
    } else {
      temp_list = append(temp_list, l)
    }
    # Omit everything after wghtloss
    if(l == "R636=WGHTLOSS"){
      break
    }
  }
  AHRQ_list # return the list
}
icd10cm_2020_1_lofregex = make_sas_list(sas_AHRQ_raw)

# Remove comformat_icd10cm_2020_1.txt
file.remove(sas_path)
