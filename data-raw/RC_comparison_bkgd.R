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

# Need to use a model that results in different values for restrictive & non-restrictive clearance

get_common_chems <- function(species = "human") {

  species_chems <- httk::get_cheminfo(info = "DTXSID",
                                species = species,
                                model = "3compartment2",
                                median.only = TRUE,
                                physchem.exclude = TRUE,
                                suppress.messages = TRUE)
  chems <- subset(
    cvt,
    analyte_dtxsid %in% species_chems
  )
  remove(list = "species_chems")
  chems <- unique(chems[["analyte_dtxsid"]])
  return(chems)
}


# Parameterize the parameters
parameterize_all <- function(species = "human", all = FALSE) {
  if (isTRUE(all)) {
    sp_chems <- httk::get_cheminfo(info = "DTXSID",
                                   species = species,
                                   model = "3compartment2",
                                   median.only = TRUE,
                                   default.to.human = FALSE,
                                   physchem.exclude = TRUE,
                                   suppress.messages = TRUE)
  } else {
    sp_chems <- get_common_chems(species = species)
  }
  function(human_clint_fup = FALSE) {
    if (isTRUE(human_clint_fup)) {
      if (isTRUE(all)) {
        sp_chems <- httk::get_cheminfo(info = "DTXSID",
                                       species = "human",
                                       model = "3compartment2",
                                       median.only = TRUE,
                                       physchem.exclude = TRUE,
                                       suppress.messages = TRUE)
      } else {
        sp_chems <- get_common_chems(species = "human")
      }
    } else {
      human_clint_fup <- FALSE
    }
    message(glue::glue("Total number of chemicals loaded: {length(sp_chems)}"))
    message(glue::glue("Fup & Clint set to Human values? {human_clint_fup}"))

    # Get restrictive and non-restrictive parameters
    params_r <- tryCatch(
      expr = {
        sapply(sp_chems,
               \(x) {
                 httk::parameterize_3comp2(
                   dtxsid = x,
                   species = species,
                   restrictive.clearance = TRUE,
                   default.to.human = human_clint_fup,
                   force.human.clint.fup = human_clint_fup,
                   suppress.messages = TRUE)
               },
               simplify = FALSE,
               USE.NAMES = TRUE
        )
      }, error = function(msg) {
        message("Not possible to parameterize_3comp2 for all dtxsid ",
                "try to use default.to.human = TRUE.")
      }
    )
    browser()

    params_nr <- tryCatch(
      expr = {
        sapply(sp_chems,
               \(x) {
                 httk::parameterize_3comp2(
                   dtxsid = x,
                   species = species,
                   restrictive.clearance = FALSE,
                   default.to.human = human_clint_fup,
                   force.human.clint.fup = human_clint_fup,
                   suppress.messages = TRUE)
               },
               simplify = FALSE,
               USE.NAMES = TRUE
        )
      }, error = function(msg) {
        message("Not possible to parameterize_3comp2 for all dtxsid ",
                "try to use default.to.human = TRUE.")
      }
    )

    if (any(is.na(params_r)) || any(is.na(params_nr))) {
      stop("Params not calculated for all DTXSIDs!")
    }

    # Collate the parameters in a data.frame
    extr_params_r <- nlist2df(params_r)
    extr_params_nr <- nlist2df(params_nr)

    # Specify changed parameters between restrictive and non-restrictive
    names(extr_params_r)[which(
      names(extr_params_r) %in% c("Clmetabolismc","Fabsgut")
    )] <- paste(
      names(extr_params_r)[which(
        names(extr_params_r) %in% c("Clmetabolismc","Fabsgut")
      )],
      "restrictive",
      sep = ".")

    names(extr_params_nr)[which(
      names(extr_params_nr) %in% c("Clmetabolismc","Fabsgut")
    )] <- paste(
      names(extr_params_nr)[which(
        names(extr_params_nr) %in% c("Clmetabolismc","Fabsgut")
      )],
      "nonrestrictive",
      sep = ".")

    # Get the additional Vdist (in case we need it) for each
    # This is the same for either clearance type (in this model)
    vdist_params <- data.frame(
      Chemical = names(params_r),
      Vdist_r = vapply(params_r,
                       \(x) {
                         httk::calc_vdist(
                           parameters = x,
                           species = species,
                           default.to.human = FALSE,
                           suppress.messages = TRUE)
                       }, FUN.VALUE = numeric(1),
                       USE.NAMES = TRUE
      )
    )

    extr_params <- cbind(
      data.frame(forced_human_values = human_clint_fup),
      merge(
        merge(extr_params_r, extr_params_nr),
        vdist_params
      )
    )

    # Calculate restrictive or nonrestrictive clearance using the
    # appropriate parameterization
    restrictive_clearance <- data.frame(
      Chemical = names(params_r),
      CLtot_restrictive = vapply(params_r,
                                 \(x) {
                                   httk::calc_total_clearance(
                                     parameters = x,
                                     species = species,
                                     model = "3compartment2",
                                     suppress.messages = TRUE,
                                     restrictive.clearance = TRUE
                                   )
                                 },
                                 FUN.VALUE = double(1)
      )
    )
    nonrestrictive_clearance <- data.frame(
      Chemical = names(params_nr),
      CLtot_nonrestrictive = vapply(params_nr,
                                 \(x) {
                                   httk::calc_total_clearance(
                                     parameters = x,
                                     species = species,
                                     model = "3compartment2",
                                     suppress.messages = TRUE,
                                     restrictive.clearance = FALSE
                                   )
                                 },
                                 FUN.VALUE = double(1)
      )
    )

    # Get half-lives of each chemical
    restrictive_halflife <- data.frame(
      Chemical = names(params_nr),
      halflife_restrictive = vapply(params_r,
                                    \(x) {
                                      httk::calc_half_life(
                                        parameters = x,
                                        species = species,
                                        model = "3compartment2",
                                        suppress.messages = TRUE,
                                        restrictive.clearance = TRUE
                                      )
                                    },
                                    FUN.VALUE = double(1)
      )
    )

    nonrestrictive_halflife <- data.frame(
      Chemical = names(params_nr),
      halflife_nonrestrictive = vapply(params_nr,
                                       \(x) {
                                         httk::calc_half_life(
                                           parameters = x,
                                           species = species,
                                           model = "3compartment2",
                                           suppress.messages = TRUE,
                                           restrictive.clearance = FALSE
                                         )
                                       },
                                       FUN.VALUE = double(1)
      )
    )

    # Assemble the clearance and halflife data.frame
    clearance_df <- merge(
      data.frame(
        Chemical = sp_chems,
        Species = species
      ),
      merge(
        merge(restrictive_clearance, nonrestrictive_clearance),
        merge(restrictive_halflife,nonrestrictive_halflife)
      )
    )

    # Return a list
    return(list(param_df = extr_params,
                clearance_df = clearance_df,
                named_params_restrictive = params_r,
                named_params_nonrestrictive = params_nr,
                species = unique(species)))
  }
  # Implicitly returns the function above
}

parameterize_rat <- parameterize_all(species = "rat")
parameterize_human <- parameterize_all(species = "human")

get_httk_preds <- function(parameters, pk_obj, species = "human") {
  pk_df <- unique.data.frame(
    subset(
      get_data(pk_obj),
      subset = (
        Media %in% "plasma" &
          Species %in% species &
          Chemical %in% names(parameters) &
          exclude %in% FALSE &
          Detect %in% TRUE
      ),
      select = c(
        Chemical, Species, Dose, Route, Media, Reference,
        N_Subjects, pLOQ, Time, Conc, Conc_trans, Conc_SD,
        data_sigma_group, Detect, exclude
      )
    )
  )
  pk_data <- unique.data.frame(
    pk_df[
      c("Chemical", "Species", "Dose", "Route", "Reference",
        "pLOQ", "Time")
    ]
  )
  group_factor <- factor(
    with(pk_data, paste(Chemical, Species, Dose, Reference, pLOQ, Route))
  )
  pk_data <- split(pk_data, group_factor)

  pk_dlist <- lapply(pk_data,
                     \(x) {
                       this_params <- parameters[[unique(x$Chemical)]]

                       tmp_solution <- suppressWarnings(
                         httk::solve_3comp2(
                           parameters = this_params,
                           dose = unique(x$Dose),
                           exp.conc = 0,
                           iv.dose = unique("iv" %in% x$Route),
                           times = unique(c(0, x$Time/24)),
                           input.units = "mg/kg",
                           output.units = "mg/L",
                           species = unique(stringr::str_to_title(x$Species)),
                           suppress.messages = TRUE
                         )
                       )

                       retdf <- data.frame(
                         Time = round(tmp_solution[, "time"]*24, digits = 5),
                         httk_Preds = pmax(
                           tmp_solution[, "Cplasma"],
                           unique(x$pLOQ)
                         )
                       )

                       # This resolves the need to match the times
                       merge(x, retdf)
                     }
  )
  return(merge(pk_df, do.call(rbind, pk_dlist)))
}

css_httk_batch <- function(param_list,
                           restrictive = TRUE) {
  stopifnot(is.logical(restrictive))

  css_list_r <- vapply(
    param_list,
    \(x) {
      calc_analytic_css_3comp2(parameters = x,
                               restrictive.clearance = restrictive)
    },
    FUN.VALUE = numeric(1L)
  )
  css_list_names <- names(param_list)

  css_df <- data.frame(Chemical = css_list_names,
                       Css_httk = css_list_r)
  if (restrictive) {
    css_df <- rename(css_df, Css_httk_r = Css_httk)
  } else {
    css_df <- rename(css_df, Css_httk_nr = Css_httk)
  }
  return(css_df)
}


filter_targets <- function(.data) {
  stopifnot("Chemical" %in% names(.data))
  dplyr::filter(.data = .data, Chemical %in% target_chems)
}
