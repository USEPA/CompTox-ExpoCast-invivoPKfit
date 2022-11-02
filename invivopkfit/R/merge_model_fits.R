#' Merge tables from fit_all from different models and select best fitting model
#'
#' This takes multiple tables of pharmcokinetic parameters (one table per model)
#' genetated by \code{\link{fit_all}}
#' and merges them into a table containing one row per chemical-species
#' combination. Each row contains all the model
#' parameters estimated for that model, and a set of optimal parameters as
#' identified from the AIC. This function assigns each chemical a column
#' "Model" containing the name of he winning (lowest AIC) model and, from that
#' winning model, reports "Vdist", "kelim", "thalf", "kgutabs", "Fbioavail".
#'
#' @param fit_list a named list of PK parameter tables from \code{\link{fit_all}}
#' -- the names should be the names of the models
#'
#' @param compound.col Column indicating chemicalname
#'
#' @param dtxsid.col Column indicating DSSTox Structure ID
#'
#' @param cas.col Column indicating Chemical Abstracts Service Registry Number
#'
#' @param species.col Column indicating species
#'
#' @param param.value.type.col Column indicating type of estimate
#'
#' @param param.value.type Type of estimate to be included
#'
#' @param data.analyzed.col Column indicating which datasets were used
#'
#' @return A data.frame with one row per chemical-species combination
#'
#' @author John Wambaugh
#'
#' @export merge_model_fits
merge_model_fits <- function(fit.list,
        compound.col="Compound",
        dtxsid.col="DTXSID",
        cas.col="CAS",
        species.col="Species",
        media.col="Media",
        param.value.type.col="param.value.type",
        data.analyzed.col="Data.Analyzed",
        param.value.type="Fitted geometric mean")
{
  main.table <- data.frame(
        DTXSID=NA,
        Compound=NA,
        CAS=NA,
        Species=NA,
        Media=NA,
        Reference=NA,
        AIC.best=NA,
        Model=NA,
        Vdist=NA,
        kelim=NA,
        kgutabs=NA,
        Fgutabs=NA)
  for (this.table.name in sort(unique(names(fit.list))))
  {
    this.table <- fit.list[[this.table.name]]
    this.table$Css <- 1/this.table$CLtot

    # select for the correct parameter type
    if (param.value.type.col %in% colnames(this.table))
    {
      this.table <- this.table[this.table[,param.value.type.col] ==
                           param.value.type,]
    } else stop(paste(param.value.type,
                      "not found in column",
                      param.value.type.col,
                      "for table",
                      this.table.name))

    # Check for pathological fits, set AIC to Inf:
    if ("kelim" %in% colnames(this.table))
    {
     # Negative kelim is bad news:
      badfits <- this.table$kelim
      badfits[is.na(badfits)] <- -0.1
      badfits <- badfits<0
      this.table[badfits,"AIC"] <- Inf
    }
    if ("Vdist" %in% colnames(this.table))
    {
     # Really large Vd is bad news:
      badfits <- this.table$Vdist
      badfits[is.na(badfits)] <- 100000
      badfits <- badfits>10000
      this.table[badfits,"AIC"] <- Inf
    }
    if ("k12" %in% colnames(this.table))
    {
      # Negative k12 indicates bad fit:
      badfits <- this.table$k12
      badfits[is.na(badfits)] <- -0.1
      badfits <- badfits<0
      this.table[badfits,"AIC"] <- Inf
    }

    for (this.dtxsid in sort(unique(this.table[,dtxsid.col])))
    {
      this.subset <- subset(this.table, this.table[,dtxsid.col] == this.dtxsid)
      this.compound <- this.subset[1,compound.col]
      this.cas <- this.subset[1,cas.col]

      for (this.medium in sort(unique(this.subset[,media.col])))
      {
        this.medium.subset <- subset(this.subset, this.subset[, media.col] ==
                                                  this.medium)

        for (this.species in sort(unique(this.medium.subset[,species.col])))
        {
          this.species.subset <- subset(this.medium.subset,
                                        this.medium.subset[,species.col] ==
                                        this.species)

          # If more than one source, require that the joint analysis worked:
          if (dim(this.species.subset)[1]>1)
          {
            browser()
            if ("Joint Analysis" %in% this.species.subset[,data.analyzed.col])
            {
              this.species.subset <- subset(this.species.subset,
                                      this.species.subset[,data.analyzed.col] ==
                                      "Joint Analysis")
            } else {
              this.species.subset <- apply(this.species.subset, 2, mean)
              this.species.subset[,data.analyzed.col] <-
                paste("Mean: ",
                  this.species.subset[,data.analyzed.col],
                  collapse=", ")
            }
          }

          # See if this chemical is already in the main fittable:
          if (this.dtxsid %in% main.table$DTXSID)
          {
            this.index <- which(main.table$DTXSID == this.dtxsid)
          } else this.index <- dim(main.table)[1]+1
          # See if this chemical-species combo is already in the main fittable:
          if (this.species %in% main.table[this.index,species.col])
          {
            this.index <- this.index & main.table$Species==this.species
          } else this.index <- dim(main.table)[1]+1

          model.postfix <- substr(this.table.name,1,5)

          model.params <- colnames(this.table)
          model.params <- model.params[!(model.params%in% c(
            "DTXSID",
            "Compound",
            "CAS",
            "Media"))]
          # Copy parameters into main table:
          for (this.param in model.params)
            main.table[this.index, paste(this.param, model.postfix)] <-
            this.species.subset[,this.param]

        } # Close species loop
      } # Close medium looop
    } # Close chemical loop
  } # Close model loop

  # Now pick the best AIC
  AIC.cols <- regexpr("AIC",colnames(main.table))
  for (this.col in AIC.cols)
    main.table[is.na(this.col), this.col] <- Inf

  for (this.chemical in sort(unique(main.table$DTXSID)))
  {
    this.subset <- subset(main.table, DTXSID==this.chemical)
    for (this.species in sort(unique(this.subset$Species)))
    {
      this.species.subset <- subset(this.subset, Species==this.species)
      # Which row of main.table do we want?
      this.index <- which(main.table$DTXSID == this.chemical &
                          main.table$Species == this.species)

      AICs <- this.species.subset[,AIC.cols]
      if (all(AICs==Inf))
      {
        # No models fit:
        main.table[this.index, "Model"] <- "None"
        main.table[this.index, "AIC.Best"] <- NA
      } else {
        # Or load all the parameters:
        min.AIC <- min(AICs)
        model.postfix <- colnames(AICS)[AICs==min.AIC]
        model.postfix <- stsplit(".", model.postfix)[[1]][2]

        main.table[this.index, "Model"] <- paste("AIC",model.postfix)
        main.table[this.index, "AIC.Best"] <- min.AIC
        for (this.param in c("Vdist", "kelim", "kgutabs", "Fgutabs"))
          if (paste(this.param, model.postfix) %in% colnames(main.table))
            main.table[this.index, this.param] <- paste(this.param, model.postfix)
      }
    }
  }

  # Calculate half-life:
  main.table$halflife <- log(2)/main.table$kelim

  # Set reasonable sig figs:
  for (this.col in c(
    "AIC.flat", "AIC.1comp", "AIC.2comp", "AIC.best", "Vdist", "kelim", "kgutabs",
    "Fgutabs", "kgutabs.1comp","kgutabs.2comp",
    "Vdist.1comp", "kelim.1comp", "Fgutabs.1comp", "Vdist.2comp", "V1.2comp",
    "k12.2comp", "k21.2comp", "kelim.2comp", "Fgutabs.2comp",
    "Ralphatokelim.2comp","Fbetaofalpha.2comp","alpha.2comp","beta.2comp",
    "halflife"))
    main.table[,this.col] <- signif(main.table[,this.col], 3)

  return(main.table)
}
