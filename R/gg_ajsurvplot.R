#' Plot the Aalen-Johansen curve with ggplot for competing risks
#' @description
#' Extend \link[survminer:ggsurvplot]{ggsurvplot} functionality to plot competing risk curve
#'
#' @param formula A standard model formula, with survival on the left and covariates on the right
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model
#' @param weights,subset,na.action,etype,count,id,timefix parameter passed to finegray model. See \code{\link[survival:finegray]{finegray}},
#' @param ... optional parameters passed to \code{ggsurvplot}.
#' @export
gg_ajsurvplot <- function(formula, data, weights, subset, na.action, etype, count, id, timefix, ...){
  dot <- list(...)
  fgargs <- match.call()[-1]
  fgargs <- as.list(fgargs[setdiff(names(fgargs), names(dot))])
  fg <- do.call(survival::finegray, fgargs)
  # browser()
  fml <- force(update(formula, Surv(fgstart,fgstop,fgstatus)~.))
  sf <- eval(substitute(survival::survfit(fml, data = fg,  weights=fgwt), list(fml=fml)))
  survminer::ggsurvplot(
    sf,
    data=fg,
    ...
  )
}
