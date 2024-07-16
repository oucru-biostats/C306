# create survival comparison table ----------------------------------------

#' Summarize results for a Cox survival model or restricted mean survival time with the treatment arm as the main covariate
#'
#' @description A function to summarize results for a survival model with the treatment arm (variable "arm") as the main covariate
#'
#' @param model a formula which can be used to fit the survival model. This formula can include other covariates than arm BUT arm must be the first covariate in the model.
#' @param data a data frame to fit the survival model.
#' @param add.risk [\code{TRUE}] a logical value specifies whether the event probability ("absolute risk") at time "infinity" should be displayed.
#' @param time [\code{Inf}] the truncation time, affecting the descriptive and the RMST model, set to \code{Inf} to perform analyses at maximum time available
#' (minimax of the observed time across two arms in RMST model)
#' @param reference.arm [\code{B}] reference arm, default to the second arm ("B"), change to "A" for base on the first arm
#' @param compare.method [\code{cox}] a string, either "cox" for CoxPH model, "cuminc" for cumulative incidence, or "rmst" for restricted mean survival time.
#' Note that if "cox" is specified and model is a mstate model, a Fine-Gray model is used.
#' If CoxPH is preferred, used Surv(t, ev == 'event-of-interest') on the LHS.
#'
#' @param compare.args: a list of additional args for compare.methods,
#'
#'  For compare.method = 'cox', it is
#'  `add.prop.haz.test` [\code{TRUE}]: a logical value specifies whether a test for proportional hazards should be added,, additional args are fed directly to `survival::coxph`.
#'  only when model is a mstate model, `cause`, default to whatever the first cause.
#'
#'  For compare.method = 'cuminc', args are fed to \code{\link[eventglm:cumincglm]{cumincglm}}
#'  `type`: [\code{diff}] a string, "diff" for difference in cumulative incidence, "ratio' for ratio of cumulative incidence,
#'  other optional args include: model.censoring, formula.censoring, ipcw.method. See \code{\link[eventglm:cumincglm]{cumincglm}} for more details.
#'
#'  For compare.method = 'rmst', args are fed to \code{\link[eventglm:rmeanglm]{rmeanglm}} ,
#'  `type`: [\code{diff}] a string, "diff" for difference in RMST, "ratio' for ratio of RMST, "lost.diff" for difference in restricted mean time lost (RMTL, = -diff), and "lost.ratio" for ratio of RMTL
#'  other optional args include: model.censoring, formula.censoring, ipcw.method. See \code{\link[eventglm:rmeanglm]{rmeanglm}} for more details.
#' @param add.prop.haz.test [\code{TRUE}] (legacy, depricated), please move this to compare.args
#' @param medsum [\code{TRUE}] a logical value, specifying whether median (IQR) of time to event should be described.
#' @param p.compare [\code{TRUE}] a logical value, specifying whether we should report p-value for the main comparison
#' @param digits [\code{2}] a number specifies number of significant digits for numeric statistics.
#' @param pdigits [\code{3}] a number specifies number of significant digits for p value.
#' @param pcutoff [\code{0.001}] a number specifies threshold value of p value to be displayed as "< pcutoff".
#' @param footer a [\code{NULL}] vector of strings to be used as footnote of table.
#' @param flextable [\code{TRUE}] a logical value specifies whether output will be a flextable-type table.
#' @param bg [\code{#F2EFEE}] a character specifies color of the odd rows in the body of flextable-type table.
#' @return a flextable-type table or a list with values/headers/footers
#'
#' @author This function was originally written by Marcel Wolbers. Lam Phung Khanh did some modification.
#' @import survival
#' @export
sstable.survcomp <- function(
    model, data, add.risk = TRUE,
    time = Inf,
    reference.arm = c('B', 'A'),
    compare.method = c('cox', 'rmst', 'cuminc'),
    compare.args = list(),
    add.prop.haz.test = TRUE, medsum = TRUE,
    p.compare = TRUE,
    digits = 2, pdigits = 3, pcutoff = 0.001, footer = NULL, flextable = TRUE, bg = "#F2EFEE"){
  requireNamespace("survival")
  
  # Initialise the table --------------------------
  ## strip the tibble class which causes issue - trinhdhk
  data <- as.data.frame(data)
  
  mf <- model.frame(update(model, .~1), data=data)
  NAs <- attr(mf, 'na.action') |> unname()
  
  if (length(NAs)){
    warning(sprintf('Missing values on observation(s) %s',
                    paste(NAs,collapse=', ')))
    # data <- data[seq_len(nrow(data))[-NAs],]
    data <- dplyr::slice(data, seq_len(nrow(data))[-NAs])
  }
  
  compare.method <- match.arg(compare.method)
  reference.arm <- match.arg(reference.arm)
  
  # arm.var <- if (length(model[[3]]) > 1) {deparse(model[[3]][[2]])} else {deparse(model[[3]])}
  arm.var <- formula.tools::rhs.vars(model)[[1]]
  if (!inherits(data[, arm.var], "factor")) data[, arm.var] <- factor(data[, arm.var])
  ### [trinhdhk] 2024-04: reverse level of y for better summary
  data[, arm.var] <- ._lv_rev_(data[, arm.var])
  arm.names <- levels(data[, arm.var])
  if (length(arm.names) > 2) stop('At the moment, only two-arm studies are supported. A lot of code has covered more than 2 arms but please think about what should be compared in those cases (chunk test or pairwise test for instance) before implementing / trinhdhk.')
  
  
  # Prepare survfit
  time2 <- if (time == Inf) .Machine$integer.max else time
  fit.surv0 <- survival::survfit(update(model, new = as.formula(paste0(". ~ ", arm.var))), data = data)
  # [Trinhdhk] Use integer.max instead of Inf b/c summary.survfit does not want Inf anymore. 05/24
  
  # Check for competing risk aka mstate -trinhdhk
  ms <- inherits(fit.surv0, 'survfitms')
  if (ms) {
    all_causes <- attr(mf[[1]], 'states')
    if(length(unique(mf[[1]][,2]))<2)
      stop('You wanted to use a mstate model. But only one event happened.')
    
    if (is.null(compare.args$cause)) compare.args$cause <- all_causes[[1]]
    if (!compare.args$cause %in% all_causes)
      stop('Specified cause is not available in the list of causes.
If you are running this in survcomp.subgroup, perhaps in one subgroup an event did not happen. I don\'t know how to fix this.')
    this_cause <- which(all_causes == compare.args$cause)+1
  }
  
  
  # Table header
  header1 <- c(paste(arm.names, " (n=", table(data[, arm.var]), ")", sep = ""), "Comparison")
  compare.stat <- switch(compare.method,
                         cox = 'HR',
                         cuminc = if (is.null(compare.args$type)) 'Cumul.inc difference'
                         else switch(compare.args$type,
                                     diff = 'Cumul.inc difference',
                                     ratio = 'Cumul.inc ratio',
                                     stop('Illegal type for cumulative incidence comparison model')),
                         rmst = if (is.null(compare.args$type)) 'RMST difference'
                         else switch(compare.args$type,
                                     diff = 'RMST difference',
                                     lost.diff = 'RMTL difference',
                                     ratio = 'RMST ratio',
                                     lost.ratio = 'RMTL ratio',
                                     stop('Illegal type for RMST comparison model')))
  
  header2 <- c(rep(ifelse(add.risk, "events/n (risk [%])", "events/n"), length(arm.names)), paste(compare.stat, if (p.compare) "(95%CI); p-value" else "(95%CI)"))
  header <- rbind(header1, header2)
  result <- rbind(header, "")
  
  ## footer
  compare.note <- switch(
    compare.stat,
    'HR' = 'HR = hazard ratio',
    'Cumul.inc difference' = 'Cumul.inc : Cumulative incidence',
    'Cumul.inc ratio' = 'Cumul.inc : Cumulative incidence',
    'RMST difference' = 'RMST: Restricted mean survival time',
    'RMTL difference' = 'RMTL: Restricted mean time loss',
    'RMST ratio' = 'RMST: Restricted mean survival time',
    'RMTL ratio' = 'RMTL: Restricted mean time lost'
  )
  compare.name <- switch(compare.method,
                         cuminc='Generalized linear models for cumulative incidence',
                         cox=if (ms) 'Fine-Gray model' else 'Cox proportional hazards model',
                         rmst='Restricted mean survival time model')
  footer <- c(
    paste0(compare.note, "; IQR = interquartile range."),
    paste(compare.stat, if (p.compare) "and p value", "were based on",
          paste0(compare.name, '.')),
    footer)
  
  
  # Descriptive analysis ---------------------------
  # add number of events and risks
  
  fit.surv <- summary(fit.surv0, time = time2, extend = TRUE)
  
  # [Trinhdhk] This is crap as always returns at inf
  # if (length(unique(data[, arm.var])) < length(arm.names)) {
  #   tmp <- fit.surv$table
  #   dim(tmp) <- c(length(unique(data[, arm.var])), length(tmp))
  #   colnames(tmp) <- names(fit.surv$table)
  # } else {
  #   tmp <- fit.surv$table
  # }
  # browser()
  n.event <- if (ms) fit.surv$n.event[,this_cause] else fit.surv$n.event
  events.n <- paste(n.event, fit.surv$n, sep = "/")
  if (add.risk) {
    cumhaz <- if (ms) fit.surv$cumhaz[, this_cause-1] else fit.surv$cumhaz
    events.n <- paste(events.n,
                      " (", formatC(100*(cumhaz), digits, format = "f"), ")", sep="")
  }
  idx <- which(arm.names %in% unique(data[, arm.var]))
  result[3, 1:length(arm.names)] <- rep("-", length(arm.names))
  result[3, idx] <- events.n
  
  
  
  # Comparison --------------------------------------
  # Re-base the arm factor
  if (reference.arm == 'B') data[, arm.var] <- ._lv_rev_(data[, arm.var])
  if (compare.method == "cox"){
    if (!is.null(compare.args$add.prop.haz.test))
      add.prop.haz.test <- compare.args$add.prop.haz.test
    if (length(events.n) < length(arm.names)) {
      result[3, length(arm.names) + 1] <- "-"
      if (add.prop.haz.test){result <- cbind(result, c("Test for proportional hazards", "p-value", "-"))}
    } else {
      compare.args$add.prop.haz.test <- NULL
      compare.args$formula <- model
      compare.args$data <- data
      
      # add HR, CI, p-value
      fit.coxph <- if (!ms)
        do.call(survival::coxph, compare.args) #survival::coxph(model, data)
      else{
        compare.args$etype <- compare.args$cause
        compare.args$cause <- NULL
        fml <- force(update(model, Surv(fgstart,fgstop,fgstatus)~.))
        fg <- do.call(survival::finegray, compare.args)
        eval(substitute(survival::coxph(fml, data = fg,  weights=fgwt), list(fml=fml)))
      }
      est <- coef(fit.coxph)[1:(length(arm.names)-1)] #trinhdhk: if there were more covariables than just arm, subsetting this to only get the arm
      hr <- formatC(exp(est), digits, format = "f")
      
      result[3, length(arm.names) + 1] <-
        ifelse(is.na(est),"-",
               
               {
                 # browser()
                 se <- sqrt(fit.coxph$var[[1]]) # trinhdhk: if there were more covariables than just arm, $var returns a var-cov mat, this first element is the var of arm (given arm is the first covariable)
                 z <- abs(est/se)
                 pval <- format.pval((1 - pnorm(z)) * 2, eps = pcutoff, digits = pdigits)
                 ci <- sapply(est, \(.est) {
                   paste(formatC(exp(c(.est - qnorm(0.975) * se, .est + qnorm(0.975) * se)), digits, format = "f"), collapse = ", ")
                 })
                 if (p.compare){
                   hr.ci.p <- paste(hr, " (", ci, "); p=", pval, sep = "")
                 } else {
                   hr.ci.p <- paste(hr, " (", ci, ")", sep = "")
                 }
                 
                 hr.ci.p
               })
      
      
      # add test for proportional hazards
      if (add.prop.haz.test){
        if (is.na(est)) {
          result <- cbind(result, c("Test for proportional hazards", "p-value", "-"))
        } else {
          attr(model, ".Environment") <- environment() # needed for cox.zph to work
          p.prop.haz <-
            if(!ms){
              cox <- do.call(survival::coxph, compare.args)
              survival::cox.zph(cox)$table[1, "p"]
            } else{
              fml <- force(update(model, Surv(fgstart,fgstop,fgstatus)~.))
              fg <- do.call(survival::finegray, compare.args)
              cox <- eval(substitute(survival::coxph(fml, data = fg,  weights=fgwt), list(fml=fml)))
              survival::cox.zph(cox)$table[1, "p"]
            }
          if (is.na(p.prop.haz)) {
            p.prop.haz <- "-"
          } else {
            p.prop.haz <- format.pval(p.prop.haz, eps = pcutoff, digits = pdigits)
          }
          result <- cbind(result, c("Test for proportional hazards", "p-value", p.prop.haz))
        }
      }
    }
  } else {
    # requireNamespace('eventglm')
    if (length(events.n) < length(arm.names)) {
      result[3, length(arm.names) + 1] <- "-"
    } else {
      compare.args$time <- time
      compare.args$add.prop.haz.test <- NULL
      compare.args$formula <- model
      compare.args$data <- data
      type <- if (is.null(compare.args$type)) 'diff' else compare.args$type
      compare.args$type <- NULL
      compare.args$link <- switch(type, "diff" = 'identity', "ratio" = 'log', 'lost.ratio' = 'log')
      if (is.null(compare.args$model.censoring)) compare.args$model.censoring <- 'stratified'
      if (is.null(compare.args$formula.censoring)) compare.args$formula.censoring <- as.formula(paste0("~`", arm.var, '`'))
      # Extract the time from Surv object with the arm included
      mf <- model.frame(update(model, new = as.formula(paste0(". ~`", arm.var, '`'))), data = data)
      mf <- cbind(unclass(mf[,1]), mf[,2])
      # Get max of stop time, is the minimum of last observed time between two arms.
      # The stop time for right cens data is the first column ($time), for interval-cens data it is the second column ($stop)
      # With laziness it is coded as mf[, ncol(mf)-2] (the last two columns are the event status from Surv object, and the arm)
      # browser()
      minimax.time <- min(by(mf[, ncol(mf)-2], mf[, ncol(mf)], max))
      max.time <- max(mf[, ncol(mf)-2])
      # tau for rmst is capped as minamax.time
      if (is.null(compare.args$time)) compare.args$time <- minimax.time
      else if (compare.args$time > minimax.time) {
        if (time != Inf) warning(
          sprintf('Specified time (tau) is later than maximum observed time between two arms (%.2f). Set too %.2f', minimax.time, minimax.time)
        )
        compare.args$time <- minimax.time
      }
      # If type == lost.ratio then reverse the time
      # i.e. Surv(tau - ev_time, ev) ~ .
      if (type == 'lost.ratio') {
        if (compare.method == 'cuminc') stop('type lost.ratio is meaningless in cuminc comparison.')
        model.lhs <- formula.tools::lhs(model)
        for (i in 2:(length(model.lhs)-1)) {
          old <- model.lhs[[i]];
          model.lhs[[i]] <- substitute(max.time - old, list(old=old, max.time=max.time))
        }
        formula.tools::lhs(model) <- model.lhs
        # browser()
        compare.args$formula <- model
      }
      
      # Perform rmeanglm
      # rmeanglm <- eventglm::rmeanglm
      fitter <- if (compare.method == 'rmst') eventglm::rmeanglm else eventglm::cumincglm
      fit.rmst <- do.call(fitter, compare.args)
      est <- coef(fit.rmst)[2:(length(arm.names))] # get the coef for arm
      invlink <- fit.rmst$family$linkinv
      if (type=='lost.diff') invlink <- function(x) -x # loss.diff=-diff
      diff <- formatC(invlink(est), digits, format = "f")
      
      result[3, length(arm.names) + 1] <-
        ifelse(is.na(est),"-",
               {
                 summary.rmst <- summary(fit.rmst)$coefficients
                 p <- summary.rmst[2:length(arm.names), 'Pr(>|z|)']
                 pval <- format.pval(p, eps = pcutoff, digits = pdigits)
                 cf <- confint(fit.rmst)[2:length(arm.names),, drop=FALSE]
                 ci <- apply(cf, 1,
                             \(.cf) paste(formatC(sort(invlink(.cf)), digits, format = "f"), collapse = ", ")
                 )
                 if (is.na(p)) diff.ci.p <- paste(diff, '(-)')
                 else if (p.compare){
                   diff.ci.p <- paste(diff, " (", ci, "); p=", pval, sep = "")
                 } else {
                   diff.ci.p <- paste(diff, " (", ci, ")", sep = "")
                 }
                 diff.ci.p
               })
    }
    
  }
  
  rownames(result) <- NULL
  
  # add label
  #browser()
  varname <- all.vars(model)[2:1]
  varlabel <- sapply(varname, function(x){
    ifelse(is.null(attr(data[, x], "label")), x, attr(data[, x], "label"))
  })
  # output
  output <- if (medsum & !ms) {
    result <- rbind(result, "", "", "", "")
    # add median (IQR) of time-to-event
    qfit <- as.matrix(do.call(cbind, quantile(fit.surv0, probs = c(0.5, 0.25, 0.75))))
    result[5:7, 1:length(arm.names)] <- apply(formatC(qfit, digits, format = "f"), 1, function(x){
      sapply(1:3, function(z) paste(x[z], " (", x[z + 3], ", ", x[z + 6], ")", sep = ""))
    })
    cbind(c("Endpoint", "", varlabel, "- Median (95%CI)", "- Lower IQR (95%CI)", "- Upper IQR (95%CI)"), result)
  } else {
    output <- cbind(c("Endpoint", "", varlabel[1]), result)
  }
  rownames(output) <- NULL
  value <- output[-c(1:2), , drop = FALSE]
  #browser()
  ## flextable
  if (flextable) {
    requireNamespace("flextable")
    requireNamespace("officer")
    
    ## main table
    tab <- flextable::flextable(as.data.frame(value))
    
    ## header
    header1 <- output[1, ]; header2 <- output[2, ]
    header2[1] <- header1[1]
    assign("tab",
           eval(parse(text = paste0("flextable::set_header_labels(tab,",
                                    paste(paste0("V", 1:length(header1)), paste0("'", header1, "'"), sep = "=", collapse = ","),
                                    ")"))))
    assign("tab",
           eval(parse(text = paste0("flextable::add_header(tab,",
                                    paste(paste0("V", 1:length(header1)), paste0("'", header2, "'"), sep = "=", collapse = ","),
                                    ", top = FALSE)"))))
    
    tab <- flextable::merge_v(tab, part = "header")
    
    if (any(value == "-")) footer <- c(footer, "- : value cannot be estimated.")
    
    for (k in (1:length(footer))) {
      tab <- flextable::add_footer(tab, V1 = footer[k], top = FALSE)
      tab <- flextable::merge_at(tab, i = k, j = 1:length(header1), part = "footer")
    }
    
    ## format
    ### width
    tab <- flextable::autofit(tab)
    ### alignment
    tab <- flextable::align(tab, j = 1, align = "left", part = "all")
    ### faces of header
    tab <- flextable::bold(tab, part = "header")
    ### background
    tab <- flextable::bg(tab, i = seq(from = 1, to = nrow(value), by = 2), j = 1:length(header1), bg = bg, part = "body")
    ### border
    tabbd <- officer::fp_border(color="black", width = 1.5)
    tab <- flextable::border_remove(tab)
    tab <- flextable::hline(tab, border = tabbd, part = "header")
    tab <- flextable::hline_top(tab, border = tabbd, part = "all")
    tab <- flextable::hline_bottom(tab, border = tabbd, part = "body")
    
  } else {
    tab <- list(table = output,
                footer = footer)
    class(tab$table) <- c('survcomp_tbl', 'ss_tbl', class(tab$table))
  }
  if (!flextable) class(tab) <- c('ss_survcomp','ss_obj')
  return(tab)
}