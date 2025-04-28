
#' Converting common quotation mark symbols into basic ASCII forms.
#'
#' `force_ascii()` Takes curly quotes commonly found in chemical names and
#' converts them to simple ASCII equivalent. This can be important for
#' cross-compatibility with other databases and systems that use simple
#' apostrophes `'` and quotes `"` instead of curly quotes or prime symbols.
#'
#'
#' @param x A character vector that may have some special characters.
#'
#' @returns Character vector with converted characters.
#' @export
#'
#' @examples
#' try(cat(force_ascii("3,3’,4,4′-tetrachloro-5-biphenylol")))
#'
force_ascii <- function(x) {
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
                     paste0(single_quote_chars, collapse = "|"),
                     ")")
  dq_regex <- paste0("(",
                     paste0(double_quote_chars, collapse = "|"),
                     ")")

  temp <- gsub(sq_regex, "'", x)
  temp <- gsub(dq_regex, '"', temp)
  temp <- gsub("(, | ,)", ",", temp)
  # Encoding(temp) <- "latin1"
  return(temp)
}
