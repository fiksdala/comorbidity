#' @title Comorbidity mapping.
#'
#' @description Maps comorbidity conditions using algorithms from the Charlson and the Elixhauser comorbidity scores.
#'
#' @param x A tidy `data.frame` (or a `data.table`; `tibble`s are supported too) with one column containing an individual ID and a column containing all diagnostic codes.
#' Extra columns other than ID and codes are discarded.
#' Column names must be syntactically valid names, otherwise they are forced to be so by calling the [make.names()] function.
#' @param id Column of `x` containing the individual ID.
#' @param code Column of `x` containing diagnostic codes.
#' Codes must be in upper case with no punctuation in order to be properly recognised.
#' @param map The mapping algorithm to be used (values are case-insensitive).
#' Possible values are the Charlson score with either ICD-10 or ICD-9-CM codes (`charlson_icd10_quan`, `charlson_icd9_quan`) and the Elixhauser score, again using either ICD-10 or ICD-9-CM (`elixhauser_icd10_quan`, `elixhauser_icd9_quan`).
#' These mapping are based on the paper by Quan et al. (2011).
#' It is also possible to obtain a Swedish (`charlson_icd10_se`) or Australian (`charlson_icd10_am`) modification of the Charlson score using ICD-10 codes.
#' @param assign0 Apply a hierarchy of comorbidities: should a comorbidity be present in a patient with different degrees of severity, then the milder form will be assigned a value of 0.
#' By doing this, a type of comorbidity is not counted more than once in each patient.
#' The comorbidities that are affected by this argument are:
#' * "Mild liver disease" (`mld`) and "Moderate/severe liver disease" (`msld`) for the Charlson score;
#' * "Diabetes" (`diab`) and "Diabetes with complications" (`diabwc`) for the Charlson score;
#' * "Cancer" (`canc`) and "Metastatic solid tumour" (`metacanc`) for the Charlson score;
#' * "Hypertension, uncomplicated" (`hypunc`) and "Hypertension, complicated" (`hypc`) for the Elixhauser score;
#' * "Diabetes, uncomplicated" (`diabunc`) and "Diabetes, complicated" (`diabc`) for the Elixhauser score;
#' * "Solid tumour" (`solidtum`) and "Metastatic cancer" (`metacanc`) for the Elixhauser score.
#'
#' @param labelled Attach labels to each comorbidity, compatible with the RStudio viewer via the [utils::View()] function.
#' Defaults to `TRUE`.
#' @param tidy.codes Tidy diagnostic codes?
#' If `TRUE`, all codes are converted to upper case and all non-alphanumeric characters are removed using the regular expression \code{[^[:alnum:]]}.
#' Defaults to `TRUE`.
#'
#' @return A data frame with `id`, columns relative to each comorbidity domain, comorbidity score, weighted comorbidity score, and categorisations of such scores, with one row per individual.
#'
#' For the Charlson score, the following variables are included in the dataset:
#' * The `id` variable as defined by the user;
#' * `ami`, for acute myocardial infarction;
#' * `chf`, for congestive heart failure;
#' * `pvd`, for peripheral vascular disease;
#' * `cevd`, for cerebrovascular disease;
#' * `dementia`, for dementia;
#' * `copd`, chronic obstructive pulmonary disease;
#' * `rheumd`, for rheumatoid disease;
#' * `pud`, for peptic ulcer disease;
#' * `mld`, for mild liver disease;
#' * `diab`, for diabetes without complications;
#' * `diabwc`, for diabetes with complications;
#' * `hp`, for hemiplegia or paraplegia;
#' * `rend`, for renal disease;
#' * `canc`, for cancer (any malignancy);
#' * `msld`, for moderate or severe liver disease;
#' * `metacanc`, for metastatic solid tumour;
#' * `aids`, for AIDS/HIV;
#'
#' Conversely, for the Elixhauser score the dataset contains the following variables:
#' * The `id` variable as defined by the user;
#' * `chf`, for congestive heart failure;
#' * `carit`, for cardiac arrhythmias;
#' * `valv`, for valvular disease;
#' * `pcd`, for pulmonary circulation disorders;
#' * `pvd`, for peripheral vascular disorders;
#' * `hypunc`, for hypertension, uncomplicated;
#' * `hypc`, for hypertension, complicated;
#' * `para`, for paralysis;
#' * `ond`, for other neurological disorders;
#' * `cpd`, for chronic pulmonary disease;
#' * `diabunc`, for diabetes, uncomplicated;
#' * `diabc`, for diabetes, complicated;
#' * `hypothy`, for hypothyroidism;
#' * `rf`, for renal failure;
#' * `ld`, for liver disease;
#' * `pud`, for peptic ulcer disease, excluding bleeding;
#' * `aids`, for AIDS/HIV;
#' * `lymph`, for lymphoma;
#' * `metacanc`, for metastatic cancer;
#' * `solidtum`, for solid tumour, without metastasis;
#' * `rheumd`, for rheumatoid arthritis/collaged vascular disease;
#' * `coag`, for coagulopathy;
#' * `obes`, for obesity;
#' * `wloss`, for weight loss;
#' * `fed`, for fluid and electrolyte disorders;
#' * `blane`, for blood loss anaemia;
#' * `dane`, for deficiency anaemia;
#' * `alcohol`, for alcohol abuse;
#' * `drug`, for drug abuse;
#' * `psycho`, for psychoses;
#' * `depre`, for depression;
#'
#' Labels are presented to the user when using the RStudio viewer (e.g. via the [utils::View()] function) for convenience.
#'
#' @details
#' The ICD-10 and ICD-9-CM coding for the Charlson and Elixhauser scores is based on work by Quan _et al_. (2005).
#' ICD-10 and ICD-9 codes must be in upper case and with alphanumeric characters only in order to be properly recognised; set `tidy.codes = TRUE` to properly tidy the codes automatically.
#' A message is printed to the R console when non-alphanumeric characters are found.
#'
#' @references Quan H, Sundararajan V, Halfon P, Fong A, Burnand B, Luthi JC, et al. _Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data_. Medical Care 2005; 43(11):1130-1139.
#' @references Charlson ME, Pompei P, Ales KL, et al. _A new method of classifying prognostic comorbidity in longitudinal studies: development and validation_. Journal of Chronic Diseases 1987; 40:373-383.
#' @references Ludvigsson JF, Appelros P, Askling J et al. _Adaptation of the Charlson Comorbidity Index for register-based research in Sweden_. Clinical Epidemiology 2021; 13:21-41.
#' @references Sundararajan V, Henderson T, Perry C, Muggivan A, Quan H, Ghali WA. _New ICD-10 version of the Charlson comorbidity index predicted in-hospital mortality_. Journal of Clinical Epidemiology 2004; 57(12):1288-1294.
#' @examples
#' set.seed(1)
#' x <- data.frame(
#'   id = sample(1:15, size = 200, replace = TRUE),
#'   code = sample_diag(200),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Charlson score based on ICD-10 diagnostic codes:
#' comorbidity(x = x, id = "id", code = "code", map = "charlson_icd10_quan", assign0 = FALSE)
#'
#' # Elixhauser score based on ICD-10 diagnostic codes:
#' comorbidity(x = x, id = "id", code = "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
#' @export

# comorbidity <- function(x, id, code, map, assign0, labelled = TRUE, tidy.codes = TRUE) {
comorbidity <-  function(x, id, code, map, assign0, factorise = FALSE, labelled = TRUE, tidy.codes = TRUE, drg = NULL, icd_rank = NULL, poa = NULL, year = NULL, quarter = NULL, icd10cm_vers = NULL){
  # set.seed(1)
    # x <- data.frame(
    #   id = 1,
    #   code = sample_diag(10),
    #   stringsAsFactors = FALSE
    # )
    # xa <- data.frame(id = 1:2, code = NA_character_)
    # xm <- rbind(x, xa)
    #
    # x = xm
    # id = "id"
    # code = "code"
    # map = "charlson_icd10_quan"
    # assign0 = FALSE
    # labelled = T
    # tidy.codes = T
  ### Check arguments
  arg_checks <- checkmate::makeAssertCollection()
  # x must be a data.frame (or a data.table)
  checkmate::assert_multi_class(x, classes = c("data.frame", "data.table", "tbl", "tbl_df"), add = arg_checks)
  # id, code, map must be a single string value
  checkmate::assert_string(id, add = arg_checks)
  checkmate::assert_string(code, add = arg_checks)
  checkmate::assert_string(map, add = arg_checks)
  # map must be one of the supported; case insensitive
  map <- tolower(map)
  #checkmate::assert_choice(map, choices = names(.maps), add = arg_checks)

  checkmate::assert_choice(map, choices = c(names(.maps),
                                            'elixhauser_ahrq_2020',
                                            'elixhauser_ahrq_2021',
                                            'elixhauser_ahrq_2022'),
											  add = arg_checks)

  # assign0, labelled, tidy.codes must be a single boolean value
  checkmate::assert_logical(assign0, add = arg_checks)
  checkmate::assert_logical(labelled, len = 1, add = arg_checks)
  checkmate::assert_logical(tidy.codes, len = 1, add = arg_checks)
  # force names to be syntactically valid:
  if (any(names(x) != make.names(names(x)))) {
    names(x) <- make.names(names(x))
    warning("Names of the input dataset 'x' have been modified by make.names(). See ?make.names() for more details.", call. = FALSE)
  }
  if (id != make.names(id)) {
    id <- make.names(id)
    warning("The input 'id' string has been modified by make.names(). See ?make.names() for more details.", call. = FALSE)
  }
  if (code != make.names(code)) {
    code <- make.names(code)
    warning("The input 'id' string has been modified by make.names(). See ?make.names() for more details.", call. = FALSE)
  }
  # id, code must be in x
  checkmate::assert_subset(id, choices = names(x), add = arg_checks)
  checkmate::assert_subset(code, choices = names(x), add = arg_checks)
  # Report if there are any errors
  if (!arg_checks$isEmpty()) checkmate::reportAssertions(arg_checks)

  ##############################################################################
  ## Isolate scores computed within comorbidity() vs. within other functions
  if (map %in% names(.maps)) {

  ### Tidy codes if required
  if (tidy.codes) x <- .tidy(x = x, code = code)

  ### Create regex from a list of codes
  regex <- lapply(X = .maps[[map]], FUN = .codes_to_regex)

  ### Subset only 'id' and 'code' columns
  if (data.table::is.data.table(x)) {
    mv <- c(id, code)
    x <- x[, ..mv]
  } else {
    x <- x[, c(id, code)]
  }

  ### Turn x into a DT
  data.table::setDT(x)

  ### Deal with missing codes
  mvb <- id
  backup <- x[, ..mvb]
  backup <- unique(backup)
  x <- stats::na.omit(x)
  # If there are no rows left (= user passed missing data only), then error:
  if (nrow(x) == 0) stop("No non-missing data, please check your input data", call. = FALSE)

  ### Get list of unique codes used in dataset that match comorbidities
  ..cd <- unique(x[[code]])
  loc <- sapply(X = regex, FUN = function(p) stringi::stri_subset_regex(str = ..cd, pattern = p))
  loc <- utils::stack(loc)
  data.table::setDT(loc)
  data.table::setnames(x = loc, new = c(code, "ind"))

  ### Merge list with original data.table (data.frame)
  x <- merge(x, loc, all.x = TRUE, allow.cartesian = TRUE, by = code)
  x[, (code) := NULL]
  x <- unique(x)

  ### Spread wide
  mv <- c(id, "ind")
  xin <- x[, ..mv]
  xin[, value := 1L]
  x <- data.table::dcast.data.table(xin, stats::as.formula(paste(id, "~ ind")), fill = 0)
  if (!is.null(x[["NA"]])) x[, `NA` := NULL]

  ### Restore IDs
  x <- merge(x, backup, by = id, all.y = TRUE, )
  data.table::setnafill(x = x, type = "const", fill = 0L)

  ### Add missing columns
  for (col in names(regex)) {
    if (is.null(x[[col]])) x[, (col) := 0L]
  }
  data.table::setcolorder(x, c(id, names(regex)))

  ### Assign zero-values to avoid double-counting comorbidities, if requested
  if (assign0) {
    x <- .assign0(x = x, map = map)
  }

  ### Turn internal DT into a DF
  data.table::setDF(x)
} else {

    if (map == 'elixhauser_ahrq_2020') {
      x <- get_ahrq_2020(x, id, code, assign0, drg, icd_rank)

    } else if (map == 'elixhauser_ahrq_2021') {
      x <- get_ahrq_2021(
        df = x,
        patient_id = id,
        icd_code = code,
        icd_seq = icd_rank,
        poa_code = poa,
        year = year,
        quarter = quarter,
        icd10cm_vers = icd10cm_vers # If NULL, vers derived from year/quarter columns
      )

    } else {
      x <- get_ahrq_2022(
        df = x,
        patient_id = id,
        icd_code = code,
        icd_seq = icd_rank,
        poa_code = poa,
        year = year,
        quarter = quarter,
        icd10cm_vers = icd10cm_vers # If NULL, vers derived from year/quarter columns
      )
    }

}
  ### Check output for possible unknown-state errors
  .check_output(x = x, id = id)

  ### Label variables for RStudio viewer if requested
  if (labelled) x <- .labelled(x = x, map = map)

  ### Return it, after adding class 'comorbidity' and some attributes
  class(x) <- c("comorbidity", class(x))
  attr(x = x, which = "map") <- map
  return(x)
}
