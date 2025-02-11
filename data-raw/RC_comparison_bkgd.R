# This is a background job script to speed up analysis time
# This is meant to be run after `setup` chunk in RestrictiveClearance_comparison.Rmd
devtools::load_all()
library(httk)
httk::load_dawson2021()
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
parameterize_all <- function(species = "human") {
  sp_chems <- get_common_chems(species = species)
  function(human_clint_fup = FALSE) {
    if (isTRUE(human_clint_fup)) {
      sp_chems <- get_common_chems(species = "human")
    } else {
      human_clint_fup <- FALSE
    }
    message(glue::glue("Total number of chemicals loaded: {length(sp_chems)}"))
    message(glue::glue("Fup & Clint set to Human values? {human_clint_fup}"))
    # Get restrictive and non-restrictive parameters
    params_r <- tryCatch(
      expr = {
        lapply(sp_chems,
               \(x) {
                 httk::parameterize_3comp2(
                   dtxsid = x,
                   species = species,
                   restrictive.clearance = TRUE,
                   default.to.human = human_clint_fup,
                   force.human.clint.fup = human_clint_fup,
                   suppress.messages = TRUE)
               }
        )
      }, error = function(msg) {
        message("Not possible to parameterize_gas_pbtk for all dtxsid ",
                "try to use default.to.human = TRUE.")
      }
    )

    params_nr <- tryCatch(
      expr = {
        lapply(sp_chems,
               \(x) {
                 httk::parameterize_3comp2(
                   dtxsid = x,
                   species = species,
                   restrictive.clearance = FALSE,
                   default.to.human = human_clint_fup,
                   force.human.clint.fup = human_clint_fup,
                   suppress.messages = TRUE)
               }
        )
      }, error = function(msg) {
        message("Not possible to parameterize_gas_pbtk for all dtxsid ",
                "try to use default.to.human = TRUE.")
      }
    )

    if (any(is.na(params_r)) || any(is.na(params_nr))) {
      stop("Params not calculated for all DTXSIDs!")
    }

    # Get the additional Vdist (in case we need it) for each
    # This is the same for either clearance type (in this model)
    vdist_params_r <- data.frame(
      Vdist_r = vapply(params_r,
                     \(x) {
                       httk::calc_vdist(
                         parameters = x,
                         species = species,
                         default.to.human = human_clint_fup,
                         suppress.messages = TRUE)
                     }, FUN.VALUE = numeric(1)
      )
    )

    # Collate the parameters in a data.frame
    extr_params_r <- do.call(rbind,
                           lapply(params_r, as.data.frame,
                                  row.names = NULL))

    names(extr_params_r)[which(
      names(extr_params_r) %in% c("Clmetabolismc","Fabsgut")
    )] <- paste(
      names(extr_params_r)[which(
        names(extr_params_r) %in% c("Clmetabolismc","Fabsgut")
      )],
      "restrictive",
      sep = ".")

    extr_params_nr <- do.call(rbind,
                           lapply(params_nr, as.data.frame,
                                  row.names = NULL))

    names(extr_params_nr)[which(
      names(extr_params_nr) %in% c("Clmetabolismc","Fabsgut")
    )] <- paste(
      names(extr_params_nr)[which(
        names(extr_params_nr) %in% c("Clmetabolismc","Fabsgut")
      )],
                                     "nonrestrictive",
                                     sep = ".")
    extr_params <- cbind(
      data.frame(
        Chemical = sp_chems,
        Species = species,
        forced_human_values = human_clint_fup),
      merge(extr_params_r, extr_params_nr),
    vdist_params_r
    )

    # Calculate restrictive or nonrestrictive clearance using the
    # appropriate parameterization
    restrictive_clearance <- vapply(params_r,
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

    nonrestrictive_clearance <- vapply(params_nr,
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


    # Get half-lives of each chemical
    restrictive_halflife <- vapply(params_r,
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

    nonrestrictive_halflife <- vapply(params_nr,
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

    # Assemble the clearance and halflife data.frame
    clearance_df <- data.frame(
      Chemical = sp_chems,
      Species = species,
      restrictive_clearance,
      nonrestrictive_clearance,
      restrictive_halflife,
      nonrestrictive_halflife
    )

    # Give the parameters the chemical names
    names(params_nr) <- sp_chems
    names(params_r) <- sp_chems

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

rat_pars_rat <- parameterize_rat()
rat_pars_human <- parameterize_rat(human_clint_fup = TRUE)

#human_pars_human <- parameterize_human()


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
        N_Subjects, pLOQ, Time, Conc
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

