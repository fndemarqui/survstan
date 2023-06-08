


#---------------------------------------------
#' Generic S3 method ggresiduals
#' @aliases ggresiduals
#' @export
#' @param object a fitted model object.
#' @details Generic method to plot residuals of survival models.
#' @param ... further arguments passed to or from other methods.
#' @return the desired residual plot.
#'
ggresiduals <- function(object, ...) UseMethod("ggresiduals")


#' ggresiduals method for survstan models
#' @aliases ggresiduals.survstan
#' @export
#' @param object a fitted model object of the class survstan.
#' @details This function produces residuals plots of Cox-Snell residuals, martingale residuals and deviance residuals.
#' @param type type of residuals used in the plot: coxsnell (default), martingale and deviance.
#' @param ... further arguments passed to or from other methods.
#' @return the desired residual plot.
#' @examples
#' \donttest{
#' library(survstan)
#' ovarian$rx <- as.factor(ovarian$rx)
#' fit <- aftreg(Surv(futime, fustat) ~ age + rx, data = ovarian, baseline = "weibull", init = 0)
#' ggresiduals(fit, type = "coxsnell")
#' ggresiduals(fit, type = "martingale")
#' ggresiduals(fit, type = "deviance")
#' }
#'
ggresiduals.survstan <- function(object, type = c("coxsnell", "martingale", "deviance"), ...){

  type <- tolower(type)
  type <- match.arg(type)
  rcs <- object$residuals
  event <- object$event
  baseline <- object$baseline
  p <- object$p


  if(type == "coxsnell"){
    ekm <- survfit(Surv(rcs, event)~1)
    tb <- data.frame(
      time = ekm$time,
      ekm = ekm$surv,
      surv = exp(-ekm$time)
    )

    ggplot(data = tb, aes(x=.data$ekm, y=.data$surv)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, color = "blue") +
      labs(x = "Kaplan-Meier", y = baseline) +
      xlim(0, 1) +
      ylim(0, 1)

  }else if(type == "martingale"){
    mf <- stats::model.frame(object)
    time <- stats::model.response(mf)[, 1]
    martingale <- event - rcs
    if(p==0){
      message("Martingale residuals plots available only for models with at least one covariate!")
    }else{
      labels <- names(mf)[-1]
      tb <- mf %>%
        select(all_of(labels))
      labels <- names(tb)
      q <- ncol(tb)

      plots <- list()
      for(j in 1:q){
        df <- data.frame(
          x = tb[,j],
          martingale = martingale
        )
        if(is.factor(df$x)){
          plots[[j]] <- ggplot(data = df, aes(x = .data$x, y = .data$martingale)) +
            geom_jitter(position=position_jitter(0.1)) +
            xlab(labels[j])
        }else{
          plots[[j]] <- ggplot(data = df, aes(x = .data$x, y = .data$martingale)) +
            geom_point() +
            geom_smooth(se = FALSE) +
            xlab(labels[j])
        }
      }
      if(q==1){
        gridExtra::grid.arrange(grobs = plots)
      }else{
        gridExtra::grid.arrange(grobs = plots, ncol = 2)
      }
    }
  }else{
    mf <- stats::model.frame(object)
    time <- stats::model.response(mf)[, 1]
    martingale <- event - rcs
    df <- data.frame(
      obs = 1:length(time),
      deviance = sign(martingale)*sqrt((-2*(martingale + event*log(event - martingale))))
    )

    ggplot(data = df, aes(x = .data$obs, y = .data$deviance)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 0 , color = "blue")

  }
}
