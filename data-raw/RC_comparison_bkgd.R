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


# Parameterize the parameters
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

# Meant to take in a data.frame with Chemical, Species, Dose, Route, Time columns
# Chemical, Species, Dose, and Route must be unique values
# Times is expected to be in hour, will be converted to days in function
cvt_httk_predictions <- function(
    # Parameterizing args
  cvt.data,
  default_to_human = FALSE,
  force_human = FALSE
) {

  stopifnot(c("Chemical", "Species", "Dose", "Route", "Time") %in% names(cvt.data))

  dtxsid <- unique(cvt.data[["Chemical"]])
  species <- unique(cvt.data[["Species"]])
  dose <- unique(cvt.data[["Dose"]])
  route <- unique(cvt.data[["Route"]])
  times <- unique(cvt.data[["Time"]])/24

  is_unique <- function(x) {
    length(unique(x)) == 1
  }

  stopifnot(is_unique(dtxsid), is_unique(species),
            is_unique(dose), is_unique(route))

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

  return(numerical_solution)

}


#
# parameterize_all <- function(species = "human", all = FALSE) {
#   if (isTRUE(all)) {
#     sp_chems <- httk::get_cheminfo(info = "DTXSID",
#                                    species = species,
#                                    model = "gas_pbtk",
#                                    median.only = TRUE,
#                                    default.to.human = FALSE,
#                                    physchem.exclude = TRUE,
#                                    suppress.messages = TRUE)
#   } else {
#     sp_chems <- get_common_chems(species = species)
#   }
#   function(human_clint_fup = FALSE) {
#     if (isTRUE(human_clint_fup)) {
#       if (isTRUE(all)) {
#         sp_chems <- httk::get_cheminfo(info = "DTXSID",
#                                        species = "human",
#                                        model = "gas_pbtk",
#                                        median.only = TRUE,
#                                        physchem.exclude = TRUE,
#                                        suppress.messages = TRUE)
#       } else {
#         sp_chems <- get_common_chems(species = "human")
#       }
#     } else {
#       human_clint_fup <- FALSE
#     }
#     message(glue::glue("Total number of chemicals loaded: {length(sp_chems)}"))
#     message(glue::glue("Fup & Clint set to Human values? {human_clint_fup}"))
#
#     # Get restrictive and non-restrictive parameters
#     params_r <- tryCatch(
#       expr = {
#         sapply(sp_chems,
#                \(x) {
#                  httk::parameterize_gas_pbtk(
#                    dtxsid = x,
#                    species = species,
#                    restrictive.clearance = TRUE,
#                    default.to.human = human_clint_fup,
#                    force.human.clint.fup = human_clint_fup,
#                    suppress.messages = TRUE)
#                },
#                simplify = FALSE,
#                USE.NAMES = TRUE
#         )
#       }, error = function(msg) {
#         message("Not possible to parameterize_3comp2 for all dtxsid ",
#                 "try to use default.to.human = TRUE.")
#       }
#     )
#
#     params_nr <- tryCatch(
#       expr = {
#         sapply(sp_chems,
#                \(x) {
#                  httk::parameterize_gas_pbtk(
#                    dtxsid = x,
#                    species = species,
#                    restrictive.clearance = FALSE,
#                    default.to.human = human_clint_fup,
#                    force.human.clint.fup = human_clint_fup,
#                    suppress.messages = TRUE)
#                },
#                simplify = FALSE,
#                USE.NAMES = TRUE
#         )
#       }, error = function(msg) {
#         message("Not possible to parameterize_3comp2 for all dtxsid ",
#                 "try to use default.to.human = TRUE.")
#       }
#     )
#
#     if (any(is.na(params_r)) || any(is.na(params_nr))) {
#       stop("Params not calculated for all DTXSIDs!")
#     }
#
#     # Collate the parameters in a data.frame
#     extr_params_r <- nlist2df(params_r)
#     extr_params_nr <- nlist2df(params_nr)
#
#     # Specify changed parameters between restrictive and non-restrictive
#     names(extr_params_r)[which(
#       names(extr_params_r) %in% c("Clmetabolismc","Fabsgut")
#     )] <- paste(
#       names(extr_params_r)[which(
#         names(extr_params_r) %in% c("Clmetabolismc","Fabsgut")
#       )],
#       "restrictive",
#       sep = ".")
#
#     names(extr_params_nr)[which(
#       names(extr_params_nr) %in% c("Clmetabolismc","Fabsgut")
#     )] <- paste(
#       names(extr_params_nr)[which(
#         names(extr_params_nr) %in% c("Clmetabolismc","Fabsgut")
#       )],
#       "nonrestrictive",
#       sep = ".")
#
#     # Get the additional Vdist (in case we need it) for each
#     # This is the same for either clearance type (in this model)
#     vdist_params <- data.frame(
#       Chemical = names(params_r),
#       Vdist_r = vapply(params_r,
#                        \(x) {
#                          httk::calc_vdist(
#                            parameters = x,
#                            species = species,
#                            default.to.human = FALSE,
#                            suppress.messages = TRUE)
#                        }, FUN.VALUE = numeric(1),
#                        USE.NAMES = TRUE
#       )
#     )
#
#     extr_params <- cbind(
#       data.frame(forced_human_values = human_clint_fup),
#       merge(
#         merge(extr_params_r, extr_params_nr),
#         vdist_params
#       )
#     )
#
#     # Calculate restrictive or nonrestrictive clearance using the
#     # appropriate parameterization
#     restrictive_clearance <- data.frame(
#       Chemical = names(params_r),
#       CLtot_restrictive = vapply(params_r,
#                                  \(x) {
#                                    httk::calc_total_clearance(
#                                      parameters = x,
#                                      species = species,
#                                      model = "gas_pbtk",
#                                      suppress.messages = TRUE,
#                                      restrictive.clearance = TRUE
#                                    )
#                                  },
#                                  FUN.VALUE = double(1)
#       )
#     )
#     nonrestrictive_clearance <- data.frame(
#       Chemical = names(params_nr),
#       CLtot_nonrestrictive = vapply(params_nr,
#                                  \(x) {
#                                    httk::calc_total_clearance(
#                                      parameters = x,
#                                      species = species,
#                                      model = "gas_pbtk",
#                                      suppress.messages = TRUE,
#                                      restrictive.clearance = FALSE
#                                    )
#                                  },
#                                  FUN.VALUE = double(1)
#       )
#     )
#
#     # Get half-lives of each chemical
#     restrictive_halflife <- data.frame(
#       Chemical = names(params_nr),
#       halflife_restrictive = vapply(params_r,
#                                     \(x) {
#                                       httk::calc_half_life(
#                                         parameters = x,
#                                         species = species,
#                                         model = "gas_pbtk",
#                                         suppress.messages = TRUE,
#                                         restrictive.clearance = TRUE
#                                       )
#                                     },
#                                     FUN.VALUE = double(1)
#       )
#     )
#
#     nonrestrictive_halflife <- data.frame(
#       Chemical = names(params_nr),
#       halflife_nonrestrictive = vapply(params_nr,
#                                        \(x) {
#                                          httk::calc_half_life(
#                                            parameters = x,
#                                            species = species,
#                                            model = "gas_pbtk",
#                                            suppress.messages = TRUE,
#                                            restrictive.clearance = FALSE
#                                          )
#                                        },
#                                        FUN.VALUE = double(1)
#       )
#     )
#
#     # Assemble the clearance and halflife data.frame
#     clearance_df <- merge(
#       data.frame(
#         Chemical = sp_chems,
#         Species = species
#       ),
#       merge(
#         merge(restrictive_clearance, nonrestrictive_clearance),
#         merge(restrictive_halflife,nonrestrictive_halflife)
#       )
#     )
#
#     # Return a list
#     return(list(param_df = extr_params,
#                 clearance_df = clearance_df,
#                 named_params_restrictive = params_r,
#                 named_params_nonrestrictive = params_nr,
#                 species = unique(species)))
#   }
#   # Implicitly returns the function above
# }
#
# parameterize_rat <- parameterize_all(species = "rat")
# parameterize_human <- parameterize_all(species = "human")
#
# get_httk_preds <- function(parameters, pk_obj, species = "human") {
#   pk_df <- unique.data.frame(
#     subset(
#       get_data(pk_obj),
#       subset = (
#         Media %in% "plasma" &
#           Species %in% species &
#           Chemical %in% names(parameters) &
#           exclude %in% FALSE &
#           Detect %in% TRUE
#       ),
#       select = c(
#         Chemical, Species, Dose, Route, Media, Reference,
#         N_Subjects, pLOQ, Time, Conc, Conc_trans, Conc_SD,
#         data_sigma_group, Detect, exclude
#       )
#     )
#   )
#   pk_data <- unique.data.frame(
#     pk_df[
#       c("Chemical", "Species", "Dose", "Route", "Reference",
#         "pLOQ", "Time")
#     ]
#   )
#   group_factor <- factor(
#     with(pk_data, paste(Chemical, Species, Dose, Reference, pLOQ, Route))
#   )
#   pk_data <- split(pk_data, group_factor)
#
#   pk_dlist <- lapply(pk_data,
#                      \(x) {
#                        this_params <- parameters[[unique(x$Chemical)]]
#
#                        tmp_solution <- suppressWarnings(
#                          httk::solve_3comp2(
#                            parameters = this_params,
#                            dose = unique(x$Dose),
#                            exp.conc = 0,
#                            iv.dose = unique("iv" %in% x$Route),
#                            times = unique(c(0, x$Time/24)),
#                            input.units = "mg/kg",
#                            output.units = "mg/L",
#                            species = unique(stringr::str_to_title(x$Species)),
#                            suppress.messages = TRUE
#                          )
#                        )
#
#                        retdf <- data.frame(
#                          Time = round(tmp_solution[, "time"]*24, digits = 5),
#                          httk_Preds = pmax(
#                            tmp_solution[, "Cplasma"],
#                            unique(x$pLOQ)
#                          )
#                        )
#
#                        # This resolves the need to match the times
#                        merge(x, retdf)
#                      }
#   )
#   return(merge(pk_df, do.call(rbind, pk_dlist)))
# }
#
# css_httk_batch <- function(param_list,
#                            restrictive = TRUE) {
#   stopifnot(is.logical(restrictive))
#
#   css_list_r <- vapply(
#     param_list,
#     \(x) {
#       calc_analytic_css_3comp2(parameters = x,
#                                restrictive.clearance = restrictive)
#     },
#     FUN.VALUE = numeric(1L)
#   )
#   css_list_names <- names(param_list)
#
#   css_df <- data.frame(Chemical = css_list_names,
#                        Css_httk = css_list_r)
#   if (restrictive) {
#     css_df <- rename(css_df, Css_httk_r = Css_httk)
#   } else {
#     css_df <- rename(css_df, Css_httk_nr = Css_httk)
#   }
#   return(css_df)
# }
#
#
# filter_targets <- function(.data, .targets = target_chems_human) {
#   stopifnot("Chemical" %in% names(.data))
#   dplyr::filter(.data = .data, Chemical %in% .targets)
#
# }
