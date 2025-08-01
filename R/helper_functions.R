.onAttach <- function(libname, pkgname) {
  startup <- function() {
    cli::cli_inform(c(
      "v" = "{.strong {.pkg invivoPKfit} v{utils::packageVersion('invivoPKfit')} has been successfully loaded!}",
      "i" = paste("For help, please run:", cli::col_cyan("help('invivoPKfit')"))
    ), class = "packageStartupMessage"
    )
  }

  startup()
}


#' Converting common quotation mark symbols into basic ASCII forms.
#'
#' `force_ascii_quotes()` Takes curly quotes commonly found in chemical names and
#' converts them to simple ASCII equivalent. This can be important for
#' cross-compatibility with other databases and systems that use simple
#' apostrophes `'` and quotes `"` instead of curly quotes or prime symbols.
#'
#'
#' @param x A character vector that may have some special characters.
#'
#' @returns Character vector with converted characters.
#'
#' @examples
#' try(cat(invivoPKfit:::force_ascii_quotes("3,3’,4,4′-tetrachloro-5-biphenylol")))
#'
force_ascii_quotes <- function(x) {
  single_quote_chars <- c("\x91", "\x92", "\xB4")
  double_quote_chars <- c("\x93", "\x94")
  prime_char <- c("\u2032")
  double_prime_char <- c("\u2033")

  # Converting the 'escaped' characters to printed characters
  Encoding(single_quote_chars) <- "latin1"
  Encoding(double_quote_chars) <- "latin1"
  Encoding(prime_char) <- "UTF-8"
  Encoding(double_prime_char) <- "UTF-8"

  # Adding the prime character as single quote
  single_quote_chars <- append(single_quote_chars, prime_char)
  double_quote_chars <- append(double_quote_chars, double_prime_char)

  # Make the regular expressions
  sq_regex <- paste0("(",
                     paste(single_quote_chars, collapse = "|"),
                     ")")
  dq_regex <- paste0("(",
                     paste(double_quote_chars, collapse = "|"),
                     ")")

  temp <- gsub(sq_regex, "'", x)
  temp <- gsub(dq_regex, '"', temp)
  temp <- gsub("(, | ,)", ",", temp)
  # Encoding(temp) <- "latin1"
  return(temp)
}

#' Converting comma-delimited string representations of numeric vectors to numerics
#'
#' @param x A character vector that represents numeric value(s)
#'
#' @return A numeric vector of values or `NA_real_`
#'
#' @examples
#' invivoPKfit:::string2num(" ")
#' invivoPKfit:::string2num("22.2")
#' invivoPKfit:::string2num("12.3,4.1,10.1")
#'
string2num <- function(x) {

  if (is.na(x)) {
    x <- ""
  }

  stopifnot(is.character(x))

  num <- trimws(x)

  if (!nzchar(num)) {
    return(NA_real_)
  }

  num <- strsplit(num, split = ",", fixed = TRUE) |> unlist() |> as.numeric()

  return(num)
}


#' Converts numeric vector to comma-delimited string representation
#'
#' @param x A numeric vector.
#' @param return.na A logical values that determines whether to return " " or NA when
#' NA values are encountered. Default: `FALSE`.
#' @return A string with comma-delimited values of `x` or an 'empty' string.
#' @examples
#' invivoPKfit:::num2string(c(1,2.3))
#' invivoPKfit:::num2string(NA_real_)
#' invivoPKfit:::num2string(NA_real_, return.na = TRUE)
#'
#'
num2string <- function(x, return.na = FALSE) {
  stopifnot(is.numeric(x) | is.na(x))

  if (length(x) == 1 && is.na(x)) {
    if (return.na) {
      return(NA)
    } else {
    return(" ")
    }
  }

  return(paste(x, collapse = ","))
}

#' Creates a string representation of a fraction
#'
#' Parentheses will enclose the numerator or denominator if there is already
#' a division symbol "/" in the expression (if it is an expression or character vector).
#'
#' @param x Numerator.
#' @param y Denominator.
#'
#' @return A character vector.
paste_fraction <- function(x, y) {
  paste_enclose <- function(string0) paste0("(", string0, ")")
  px <- x
  py <- y
  if (grepl("/", x, fixed = TRUE)) px <- paste_enclose(x)
  if (grepl("/", y, fixed = TRUE)) py <- paste_enclose(y)
  paste0(px, "/", py)
}


#' Wrapper for [dplyr::case_match()] that accepts named lists as arguments
#' @param .x A vector to match against. In a data context, a column.
#' @param .lst A named list in the form `old_value = new_value`.
#' @param .default Value to be used when there are no values. Suggest setting it to `.x`
#' @return A vector of replaced values
case_match_lst <- function(.x, .lst, .default = NULL) {

  stopifnot(rlang::is_named(.lst))

  dplyr::case_match(
    .x,
    !!!lapply(paste0("\"", names(.lst),"\"~\"", .lst, "\""), stats::as.formula),
    .default = .default
  )
}


#' Wrapper for [stringr::str_replace_all()] which replaces whole words.
#' @param .x A vector to match against.
#' @param .lst A named vector in the form `old_value = new_value`.
#' @return A character vector with values substituted in.
str_replace_lst <- function(.x, .lst) {

  stopifnot(rlang::is_named(.lst))

  stringr::str_replace_all(.x, setNames(.lst, paste0("^", names(.lst), "$")))
}


#' Creating a simple test CvT dataset
#'
#' @param params A named numeric vector of model parameter values.
#' @param time A numeric vector of times, reflecting the time point when
#'  concentration is measured after the corresponding single bolus dose. Must be
#'  same length as `dose` and `iv.dose`, or length 1.
#' @param dose A numeric vector of doses, reflecting single bolus doses
#'  administered at time 0. Must be same length as `time` and `iv.dose`, or
#'  length 1. In this model, it is expected that this value represents a measurement
#'  of radioactive particles from a radiolabeling experiment.
#' @param route A character vector, reflecting the route of administration of
#'  each single bolus dose: `'oral'` or `'iv'`.  Must be same length as `time`
#'  and `dose`, or length 1.
#' @param medium A character vector reflecting the medium in which each resulting
#'  concentration is to be calculated: "blood" or "plasma". Default is "plasma".
#'  Must be same length as `time` and `dose`, or length 1.
#' @param N Numeric, positive and non-zero integer. Number of individual subjects.
#' @param var Numeric between 0 and 1. Describes variation in the measurements.
#'
#' @return A data frame with concentration over time data.

pseudo_cvt <- function(
    params = c(
      Clint = 10,
      Q_gfr = 0.31,
      Q_totli = 0.743,
      Fup = 0.2,
      Vdist = 1.2,
      Fgutabs = 0.75,
      kgutabs = 0.3,
      Rblood2plasma = 0.8,
      Frec = 0.95
    ),
    time = seq(0, 30, by = 0.5),
    dose = 100,
    route = c("oral", "iv"),
    medium = c("blood", "plasma", "excreta"),
    N = 4,
    var = 0.5) {

  n_routes <- length(route)
  n_media <- length(medium)
  n_subjects <- N
  timepoints <- time
  n_timepoints <- length(timepoints)
  total_combinations <- n_routes * n_media * n_subjects * n_timepoints

  if (total_combinations > 10000) {
    if (n_subjects > 5) {
      n_subjects <- 3
      message("Setting subjects to 3.")
    }

    ideal_timepoints <- 10000 %/% (n_media * n_subjects * n_routes)
    min_tp <- min(time)
    max_tp <- max(time)
    timepoints <- seq(min_tp, max_tp, length.out = ideal_timepoints)
    message("Reduced number of timepoints to ", ideal_timepoints)
  }


  tmp_subject_time <- data.frame(
    Time = rep(time, each = n_subjects),
    Subject_ID = rep(seq(1, n_subjects), n_timepoints)
  )

  tmp_route_media <- data.frame(
    Route = rep(route, each = n_media),
    Medium = rep(medium, n_routes)
  )

  tmp_cvt <- tidyr::crossing(tmp_subject_time, tmp_route_media)

  tmp_cvt$Conc <- with(
    tmp_cvt,
    cp_1comp_rad(params, time = Time, dose,
                 route = Route, medium = Medium)
  )

  tmp_cvt <- tmp_cvt |>
    dplyr::rowwise() |>
    dplyr::mutate(Conc = pmin(dose,
                              pmax(sqrt(.Machine$double.eps),
                                   rnorm(1, mean = Conc, sd = var*log2(dose))
                              )
    )) |>
    ungroup()

  return(tmp_cvt)
}


#' Creates a [cli::cli_fmt()] output for the pattern "name(x), length = "
#'
#' @param x A named vector with a single value for each element.
#' @param extra A character vector that should be printed between each name and value.
#' @return A cli::cli_fmt output list.
#' @author Gilberto Padilla Mercado
#'

post_name_value <- function(x, extra = "") {
  cli::cli_fmt({
    cli::cli_ul()
    for (k in seq_along(x)) {
      cli::cli_li(paste("{names(x[k])},", extra, "{x[k]}"))
    }
    cli::cli_end()
  })
}


