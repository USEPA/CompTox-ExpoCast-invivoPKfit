# This is a background job script to speed up analysis time
# This is meant to be run after `setup` chunk in RestrictiveClearance_comparison.Rmd
devtools::load_all()
library(httk)
httk::load_dawson2021()

get_common_chems <- function(species = "human") {

  species_chems <- get_cheminfo(info = "DTXSID",
                                species = species,
                                model = "gas_pbtk",
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

parameterize_all_gas <- function(species = "human") {
  sp_chems <- get_common_chems(species = species)
  function(human_clint_fup = FALSE) {
    if (isTRUE(human_clint_fup)) {
      sp_chems <- get_common_chems(species = "human")
    } else {
      human_clint_fup <- FALSE
    }
  message(glue::glue("Total number of chemicals loaded: {length(sp_chems)}"))
  message(glue::glue("Fup & Clint set to Human values? {human_clint_fup}"))
    params <- tryCatch(
      expr = {
        lapply(sp_chems,
               \(x) {
                 httk::parameterize_gas_pbtk(
                   dtxsid = x,
                   species = species,
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

    if (any(is.na(params))) stop("Params not calculated for all DTXSIDs!")

    extr_params <- do.call(rbind,
                           lapply(params, as.data.frame,
                                  row.names = NULL))
    extr_params <- cbind(
      data.frame(
        Chemical = sp_chems,
        Species = species,
        forced_human_values = human_clint_fup),
      extr_params)

    restrictive_clearance <- vapply(params,
                                    \(x) {
                                      httk::calc_total_clearance(
                                        parameters = x,
                                        species = species,
                                        suppress.messages = TRUE,
                                        restrictive.clearance = TRUE
                                      )
                                    },
                                    FUN.VALUE = double(1)
    )

    nonrestrictive_clearance <- vapply(params,
                                       \(x) {
                                         httk::calc_total_clearance(
                                           parameters = x,
                                           species = species,
                                           suppress.messages = TRUE,
                                           restrictive.clearance = FALSE
                                         )
                                       },
                                       FUN.VALUE = double(1)
    )
    clearance_df <- data.frame(
      Chemical = sp_chems,
      Species = species,
      restrictive_clearance,
      nonrestrictive_clearance
    )

    names(params) <- sp_chems

    return(list(param_df = extr_params,
                clearance_df = clearance_df,
                named_params = params,
                species = unique(species)))
  }
  # Implicitly returns the function above
}

parameterize_rat <- parameterize_all_gas(species = "rat")
parameterize_human <- parameterize_all_gas(species = "human")

rat_pars_rat <- parameterize_rat()
rat_pars_human <- parameterize_rat(human_clint_fup = TRUE)

human_pars_human <- parameterize_human()


get_httk_preds <- function(parameters, pk_obj, species = "human") {
  this_species = parameters$species
  all_params = parameters$named_params
  if (this_species != species) {
    message("Species do not match! ", "Parameters from species:",
            this_species, " and ", species, "data will be used.")
  }
  pk_df <- unique.data.frame(
    subset(
      get_data(pk_obj),
      subset = (
        Media %in% "plasma" &
          Species %in% species &
          Chemical %in% names(all_params) &
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
                       this_params <- all_params[[unique(x$Chemical)]]

                       tmp_solution <- httk::solve_gas_pbtk(
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

