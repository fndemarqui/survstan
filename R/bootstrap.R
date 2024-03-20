
rename_mf <- function(mf){
  vars <- names(mf)
  off <- vars[grep("offset", vars)]
  substr(off, start = 8, stop=nchar(off)-1)
  vars[grep("offset", vars)] <- substr(off, start = 8, stop=nchar(off)-1)
  names(mf) <- vars
  return(mf)
}

bootstrap <- function(object, nboot, cores, ...){

  # cl <- parallel::makeCluster(cores)
  # doParallel::registerDoParallel(cl)
  if(cores > 1){
    future::plan(future::multisession, workers = cores)
  }

  survreg <- object$survreg
  baseline <- object$baseline
  formula <- object$formula
  formula <- stats::update(formula, survival::Surv(time, status) ~ .)
  mf <- object$mf
  resp <- stats::model.response(mf)

  time <- resp[,1]
  status <- resp[,2]

  offset <- stats::model.offset(mf)
  if(!is.null(offset)){
    mf <- rename_mf(mf)
  }

  mf <- mf %>%
    dplyr::select(-dplyr::starts_with("Surv("))

  data <- data.frame(
    time = time,
    status = status
  ) %>%
    dplyr::bind_cols(mf)

  n <- object$n
  p <- object$p
  tau <- object$tau


  index <- 1:n
  index1 <- which(status==1)
  index2 <- which(status==0)
  n1 <- length(index1)
  n2 <- length(index2)
  pars_hat <- estimates(object)
  k <- length(estimates(object))

  if(cores>1){
    pars <- foreach::foreach(
      b = 1:nboot, .combine = rbind,
      .options.future = list(seed = TRUE)
    ) %dofuture% {
      samp1 <- sample(index1, size=n1, replace=TRUE)
      samp2 <- sample(index2, size=n2, replace=TRUE)
      samp <- c(samp1, samp2)
      mydata <- dplyr::slice(data, samp)
      switch(survreg,
        "yp" = suppressWarnings({invisible(utils::capture.output(fit <- survstan::ypreg(formula, data=mydata, baseline = baseline)))}),
        "eh" = suppressWarnings({invisible(utils::capture.output(fit <- survstan::ehreg(formula, data=mydata, baseline = baseline)))}),
        "ah" = suppressWarnings({invisible(utils::capture.output(fit <- survstan::ahreg(formula, data=mydata, baseline = baseline)))})
      )
      if(!is(object, "try-error")){
        survstan::estimates(fit)
      }
    }
  }else{
    pars <- foreach::foreach(
      b = 1:nboot, .combine = rbind
    ) %do% {
      samp1 <- sample(index1, size=n1, replace=TRUE)
      samp2 <- sample(index2, size=n2, replace=TRUE)
      samp <- c(samp1, samp2)
      mydata <- dplyr::slice(data, samp)
      switch(survreg,
        "yp" = suppressWarnings({invisible(utils::capture.output(fit <- survstan::ypreg(formula, data=mydata, baseline = baseline)))}),
        "eh" = suppressWarnings({invisible(utils::capture.output(fit <- survstan::ehreg(formula, data=mydata, baseline = baseline)))}),
        "ah" = suppressWarnings({invisible(utils::capture.output(fit <- survstan::ahreg(formula, data=mydata, baseline = baseline)))})
      )
      if(!is(object, "try-error")){
        survstan::estimates(fit)
      }
    }
  }

  return(pars)
}
