# This is a background job script to speed up analysis time
# This is meant to be run after `setup` chunk in RestrictiveClearance_comparison.Rmd
devtools::load_all()
library(httk)
httk::load_dawson2021()

# Convenience function that takes named list and outputs data.frame with names as column
nlist2df <- function(.list = list(),
                     .name = "Chemical") {
  stopifnot(!is.null(names(.list))) # must be named list

  this_df <-  do.call(rbind,
                      lapply(.list,
                             \(x) {
                               as.data.frame(x, row.names = NULL)
                             }
                      )
  )
  this_df[.name] <- rownames(this_df)
  rownames(this_df) <- NULL

  return(this_df)
}

# From distribution of values (median, lower 95% & upper 95% CI)
# get the coefficient of variation or return a default value
calc_coefvar <- function(x, N = 3, default = 0.1) {
  if (is.na(x)) return(default)

  stopifnot(is.character(x), length(x) == 1)

  num_vec <- as.numeric(unlist(strsplit(x, split = ",")))
  if (length(num_vec) == 1 || num_vec[1] == 0) return(default)

  # Length must be greater than one
  sem = abs(diff(num_vec[2:3]))/(2 * 1.96)
  sdev = sem * sqrt(N)

  return(sdev/num_vec[1])

}

# Need to use a model that results in different values for restrictive & non-restrictive clearance
get_common_chems <- function(species = "human") {

  species_chems <- httk::get_cheminfo(info = "DTXSID",
                                species = species,
                                model = "gas_pbtk",
                                median.only = TRUE,
                                physchem.exclude = TRUE,
                                suppress.messages = TRUE)
  intersect(species_chems, cvt[['analyzed_chem_dtxsid']])
}

# Get the 95% Confidence interval from the simulations
get_95CI_mc <- function(mat, suffix = "_R") {
  # Generate Column Names and Extract Plasma Concentrations
  CPLASMA_DF <- as.data.frame(
    mat[,
        c(as.list(quantile(Cplasma, c(0.025, 0.5, 0.975)))),
        by = time
    ]
  )
  colnames(CPLASMA_DF) <- c("Time", paste0(c("l95", "MED", "u95"), suffix))
  CPLASMA_DF
}
# Similar to above but extract the means and standard deviation

get_meansd_mc <- function(lst, suffix = "_R") {
  # Generate Column Names and Extract Plasma Concentrations
  CPLASMA_MEANS <- as.data.frame(lst[['means']])[c("time", "Cplasma")]
  colnames(CPLASMA_MEANS) <- c("Time", paste0("MEAN", suffix))
  CPLASMA_SDS <- as.data.frame(lst[['sds']])[c("time", "Cplasma")]
  colnames(CPLASMA_SDS) <- c("Time", paste0("SD", suffix))
  merge(CPLASMA_MEANS, CPLASMA_SDS)
}


# Calculate the overlap between two lines defined by a matrix 'mat'
# common axis is mat[,1] aka first column
# overlap defined by places where mat[,3] > mat[,2]
calc_overlap <- function(mat) {
  ids <- colnames(mat)
  stopifnot(
    "Time" %in% ids,
    NCOL(mat) == 3,
    "Time" == ids[1]
  )

  auc_out <- numeric(1L)

  message("Calculating overlap between ",
          ids[2], " and ", ids[3], " when ",
          ids[2], " < ", ids[3])

  LTH = mat[, 2] < mat[, 3]
  EQT = mat[, 2] == mat[, 3]

  if (!any(LTH)) {
    auc_out <- 0
    return(auc_out)
  }

  for (i in seq(1, (NROW(mat) - 1))) {

    itv <- c(i, i + 1) # interval
    matflor <- mat[itv, c(1,2)] # line of "lower" interval
    matceil <- mat[itv, c(1,3)] # line of "upper" interval

    if (!any(LTH[itv])) { # no overlap or same line along interval
      next
    }

    if (all(LTH[itv])) { # overlap, interval forms quadrilateral
      # QUAD
      this_area <- sum(
        abs(det(
          cbind(rbind(matflor, matceil[1, ]), 1)
        )),
        abs(det(
          cbind(rbind(matflor[2, ], matceil), 1)
        ))
      ) / 2
      if (this_area <= sqrt(.Machine$double.neg.eps)) this_area <- 0

      auc_out <- auc_out + this_area
      next

    }

    if (any(EQT[itv])) { # forms triangle
      # TRIANGLE
      Nmidpt <- !EQT[itv]
      this_area <- abs(det(cbind(rbind(matflor, matceil[Nmidpt,]), 1)) / 2)
      if (this_area <= sqrt(.Machine$double.neg.eps)) this_area <- 0

      auc_out <- auc_out + this_area
      message("Triangle Area:\t", this_area)
      next

    } else { # forms triangle, but midpoint is unmeasured
      # MIDPOINT
      del1 <- diff(matflor) # row-vector (dx, dy1)
      del2 <- diff(matceil) # row-vector (dx, dy2)
      #midpoint x: (b_2 - b_1) / (m_1 - m_2)
      midpt <- (matceil[1,2] - matflor[1,2]) / abs(det(rbind(del1,del2)))
      midpt <- c( midpt, del2[1] * midpt + matceil[2, 1])
      # TRIANGLE
      Nmidpt <- LTH[itv] # Not the midpoint but where less than
      this_area <- abs(det(
        cbind(rbind(midpt, matflor[Nmidpt, ], matceil[Nmidpt, ]), 1)
      ) / 2)
      if (this_area <= sqrt(.Machine$double.neg.eps)) this_area <- 0

      message("Triangle Area using midpoint:\t", this_area)
      auc_out <- auc_out + abs(this_area)

    }
  }
  return(auc_out)
}


# Parameterize both restrictive and non-restrictive gas_pbtk
standard_parameterization <- function(dtxsid, species = "Human",
                                      default_to_human = FALSE,
                                      force_human = FALSE,
                                      par_fun = "parameterize_gas_pbtk") {
  # Check the dtxsid
  chem_check <- dtxsid %in% get_cheminfo(info = "dtxsid", species = species,
                                         model = "gas_pbtk",
                                         suppress.messages = TRUE)
  if (!chem_check) stop(paste(dtxsid, "is not available for gas_pbtk model"))

  message("Parameterizing gas_pbtk model")
  par_rest <- do.call(
    par_fun,
    args = list(
      dtxsid = dtxsid,
      species = species,
      default.to.human = default_to_human,
      force.human.clint.fup = force_human,
      Caco2.options = list(overwrite.invivo = TRUE),
      restrictive.clearance = TRUE,
      suppress.messages = TRUE
    )
  )
  par_nonrest <- do.call(
    par_fun,
    args = list(
      dtxsid = dtxsid,
      species = species,
      default.to.human = default_to_human,
      force.human.clint.fup = force_human,
      Caco2.options = list(overwrite.invivo = TRUE),
      restrictive.clearance = FALSE,
      suppress.messages = TRUE
    )
  )

  list(par_rest = par_rest, par_nonrest = par_nonrest)
}

# Get the httk predictions for restrictive and nonrestrictive clearance
standard_httk_preds <- function(dtxsid,
                                parameters = NULL,
                                species = "Human",
                                default_to_human = FALSE,
                                force_human = FALSE,
                                dose = numeric(),
                                route = "iv",
                                times = numeric()) {
  IV.DOSE <- ifelse(route == "iv", TRUE, FALSE)
  # Decided to do both Clearances in one go
  # parameterization to be able to explicitly force human clint and fup

  if (is.null(parameters) || !all(c("par_rest", "par_nonrest") %in% names(parameters))) {
    parameters <- standard_parameterization(
      dtxsid = dtxsid,
      species = species,
      default_to_human = default_to_human,
      force_human = force_human
    )
  }

  message("Solving numerically...")
  numerical_solution_rest <- solve_gas_pbtk(
    parameters = parameters[["par_rest"]],
    species = species,
    times = times,
    restrictive.clearance = TRUE,
    dose = dose,
    iv.dose = IV.DOSE,
    exp.conc = 0, # This is needed for gas_pbtk
    default.to.human = default_to_human,
    input.units = "mg/kg", # Needed for gas_pbtk
    output.units = "mg/L",
    monitor.vars = "Cplasma",
    Caco2.options = list(overwrite.invivo = TRUE),
    suppress.messages = TRUE
  ) |>
    suppressWarnings() |>
    as.data.frame()

  numerical_solution_nonrest <- solve_gas_pbtk(
    parameters = parameters[["par_nonrest"]],
    species = species,
    times = times,
    restrictive.clearance = FALSE,
    dose = dose,
    iv.dose = IV.DOSE,
    exp.conc = 0,
    default.to.human = default_to_human,
    input.units = "mg/kg",
    output.units = "mg/L",
    monitor.vars = "Cplasma",
    Caco2.options = list(overwrite.invivo = TRUE),
    suppress.messages = TRUE
  ) |>
    suppressWarnings() |>
    as.data.frame()

  numerical_solution_rest <- numerical_solution_rest[c("time", "Cplasma")]
  colnames(numerical_solution_rest) <- c("Time", "Sol_R")
  numerical_solution_nonrest <- numerical_solution_nonrest[c("time", "Cplasma")]
  colnames(numerical_solution_nonrest) <- c("Time", "Sol_N")

  return(merge(numerical_solution_rest, numerical_solution_nonrest))
}

# Get httk predictions with MC sampling of select parameters
standard_httk_mc <- function(
    # Parameterizing args
  dtxsid,
  parameters = NULL,
  species = "Human",
  default_to_human = FALSE,
  force_human = FALSE,
  # Solving args
  dose = numeric(),
  route = "oral",
  times = numeric(),
  # MC sampler args
  vary.params = list(Clint = 0.1),
  censored.params = list(Funbound.plasma = list(cv = 0.1, lod = 1E-4)),
  samples = 1000
) {
  IV.DOSE <- ifelse(route == "iv", TRUE, FALSE)

  if (is.null(parameters) || !all(c("par_rest", "par_nonrest") %in% names(parameters))) {
    parameters <- standard_parameterization(
      dtxsid = dtxsid,
      species = species,
      default_to_human = default_to_human,
      force_human = force_human
    )
  }

  message("Solving Monte-Carlo Samples...")
  message("Restrictive...")
  # Make the MC sampled simulations
  mc_solution_rest <- calc_mc_tk(
    parameters = parameters[["par_rest"]],
    samples = samples,
    species = species,
    model = "gas_pbtk",
    httkpop = FALSE,
    output.units = "mg/L",
    invitrouv = FALSE,
    parameterize.args.list = list(
      restrictive.clearance = TRUE,
      default.to.human = default_to_human,
      force.human.clint.fup = force_human,
      # adjusted.Funbound.plasma = FALSE,
      # adjusted.Clint = FALSE,
      Caco2.options = list(overwrite.invivo = TRUE)
    ),
    solvemodel.arg.list = list(
      times = times,
      restrictive.clearance = TRUE,
      dose = dose,
      iv.dose = IV.DOSE,
      default.to.human = default_to_human,
      exp.conc = 0,
      input.units = "mg/kg",
      output.units = "mg/L",
      monitor.vars = "Cplasma"
    ),
    propagate.invitrouv.arg.list = list(
      restrictive.clearance = TRUE,
      default.to.human = default_to_human,
      force.human.clint.fup = force_human,
      Caco2.options = list(overwrite.invivo = TRUE)
    ),
    Caco2.options = list(overwrite.invivo = TRUE),
    vary.params = vary.params,
    censored.params = censored.params,
    return.all.sims = TRUE
  ) |>
    suppressWarnings()

  mc_solution_rest <- get_95CI_mc(mc_solution_rest[['sims']], suffix = "_R")

  message("Nonrestrictive...")
  mc_solution_nonrest <- calc_mc_tk(
    parameters = parameters[["par_nonrest"]],
    samples = samples,
    species = species,
    model = "gas_pbtk",
    httkpop = FALSE,
    output.units = "mg/L",
    invitrouv = FALSE,
    parameterize.args.list = list(
      restrictive.clearance = FALSE,
      default.to.human = default_to_human,
      force.human.clint.fup = force_human,
      # adjusted.Funbound.plasma = FALSE,
      # adjusted.Clint = FALSE,
      Caco2.options = list(overwrite.invivo = TRUE)
    ),
    solvemodel.arg.list = list(
      times = times,
      restrictive.clearance = FALSE,
      dose = dose,
      iv.dose = IV.DOSE,
      default.to.human = default_to_human,
      exp.conc = 0,
      input.units = "mg/kg",
      output.units = "mg/L",
      monitor.vars = "Cplasma"
    ),
    propagate.invitrouv.arg.list = list(
      restrictive.clearance = FALSE,
      default.to.human = default_to_human,
      force.human.clint.fup = force_human,
      Caco2.options = list(overwrite.invivo = TRUE)
    ),
    Caco2.options = list(overwrite.invivo = TRUE),
    vary.params = vary.params,
    censored.params = censored.params,
    return.all.sims = TRUE
  ) |>
    suppressWarnings()

  mc_solution_nonrest <- get_95CI_mc(mc_solution_nonrest[['sims']], suffix = "_N")

  return(merge(mc_solution_rest, mc_solution_nonrest))
}

# Get all simulations (restrictive and nonrestrictive)
conc_fun_wsamples <- function(
    # Parameterizing args
  dtxsid,
  species = "Human",
  default_to_human = FALSE,
  force_human = FALSE,
  # Solving args
  dose = 20,
  route = "iv",
  times = seq(0, 48, by = 1)/24,
  # MC sampler args
  vary.params = NULL,
  censored.params = NULL,
  samples = 1000
) {
  # Decided to do both Clearances in one go
  # parameterization to be able to explicitly force human clint and fup

  # Parameterize
  gas_pars <- standard_parameterization(
    dtxsid = dtxsid,
    species = species,
    default_to_human = default_to_human,
    force_human = force_human
  )

  # SOLVE NUMERICALLY
  numerical_solution <- standard_httk_preds(parameters = gas_pars,
                                            species = species,
                                            default_to_human = default_to_human,
                                            force_human = force_human,
                                            dose = dose,
                                            route = route,
                                            times = times)

  mc_solutions <- standard_httk_mc(parameters = gas_pars,
                                   species = species,
                                   default_to_human = default_to_human,
                                   force_human = force_human,
                                   dose = dose,
                                   route = route,
                                   times = times,
                                   vary.params = vary.params,
                                   censored.params = censored.params,
                                   samples = samples)

  # Merge the solutions
  tmp <- merge(numerical_solution, mc_solutions)
  tmp <- subset(tmp, Time %in% times)

  # Time from days to hours
  tmp$Time <- janitor::round_to_fraction(tmp$Time * 24, denominator = 60)

  return(tmp)
}


# Meant to take in a data.frame with Chemical, Species, Dose, Route, Time columns
# Chemical, Species, Dose, and Route must be unique values
# Times is expected to be in hour, will be converted to days in function
cvt_httk_predictions <- function(
  x,
  vary.params = VARPARS,
  censored.params = CENSORPARS,
  default_to_human = FALSE,
  force_human = FALSE
) {

  stopifnot(c(
    "DTXSID", "Species", "Dose",
    "Media", "Route", "Time"
    ) %in% names(x))

  dtxsid <- unique(x[["DTXSID"]])
  species <- tools::toTitleCase(unique(x[['Species']]))
  media <- unique(x[['Media']])
  route <- unique(x[['Route']])
  dose <- unique(x[['Dose']])
  times <- unique(x[['Time']])/24 # Multiple values (in hours)

  is_unique <- function(x) {
    length(unique(x)) == 1
  }

  stopifnot(is_unique(dtxsid), is_unique(species),
            is_unique(dose), is_unique(route))

  invivo_preds <- conc_fun_wsamples(
    dtxsid = dtxsid,
    species = species,
    dose = dose,
    route = route,
    times = times,
    samples = 1000,
    force_human = force_human,
    vary.params = vary.params,
    censored.params = censored.params
  ) |> try()

  if (inherits(invivo_preds, "try-error")) return(NULL)

  # If measured media is blood need to divide CPLASMA results by Rblood2plasma
  if (media == "blood") {
    Rb2p <- standard_parameterization(
      dtxsid = dtxsid,
      species = species,
      force_human = force_human
    )[['par_rest']][['Rblood2plasma']]

    cols2convert <- apply(
      expand.grid(c("Sol", "l95", "MED", "u95"), c("_R", "_N")),
      MARGIN = 1,
      paste0, collapse = ""
    )

    invivo_preds <- invivo_preds |>
      dplyr::mutate(dplyr::across(dplyr::all_of(cols2convert), function(x) x / Rb2p))
  }

  cbind(
    data.frame(
      DTXSID = dtxsid,
      Species = tolower(species),
      Dose = dose,
      Route = route,
      Media = media
    ),
    as.data.frame(censored.params),
    as.data.frame(vary.params),
    invivo_preds
  )

}

