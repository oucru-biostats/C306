#' Summarize results for a survival model by treatment arm and subgroup
#'
#' @description A function to summarize results for a survival model by treatment arm (variable "arm") and subgroup.
#'
#' @param base.model a formula from which sub-group specific estimates are extracted (!! arm must be the first covariate in the model).
#' @param overall.model an optional one-sided formula for additional terms for overall population only
#' @param subgroup.model a formula of the form "~subgrouping.variable1+subgrouping.variable2" (!! subgrouping.variable must be factors and there should be nothing on the left-hand side of the formula).
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
#' @param p.compare [\code{TRUE}] a logical value, specifying whether we should report p-value for the main comparison
#' @param digits [\code{2}] a number specifies number of significant digits for numeric statistics.
#' @param pdigits [\code{3}] a number specifies number of significant digits for p value.
#' @param pcutoff [\code{0.001}] a number specifies threshold value of p value to be displayed as "< pcutoff".
#' @param footer a [\code{NULL}] vector of strings to be used as footnote of table.
#' @param flextable [\code{TRUE}] a logical value specifies whether output will be a flextable-type table.
#' @param bg [\code{#F2EFEE}] a character specifies color of the odd rows in the body of flextable-type table.
#' @param overall [\code{TRUE}] where to print overall model
#' @param ... arguments that are passed to sstable.survcomp
#'
#' @return a flextable-type table or a list with values/headers/footers
#'
#' @author This function was originally written by Marcel Wolbers. Trinh Dong and Lam Phung Khanh did some modification.
#' @import survival
#' @export
sstable.survcomp.subgroup <- function(base.model, subgroup.model, overall.model, data,
                                      time = Inf,
                                      reference.arm = c('B', 'A'),
                                      compare.method = c('cox', 'rmst', 'cuminc'),
                                      compare.args = list(),
                                      p.compare = TRUE,
                                      digits = 2, pdigits = 3, pcutoff = 0.001, footer = NULL, flextable = TRUE, bg = "#F2EFEE", overall = TRUE,...){
  
  
  requireNamespace("survival")
  
  ## strip the tibble class which causes issue - trinhdhk
  data <- as.data.frame(data)
  NAs <- model.frame(update(base.model, .~1), data=data) |> attr('na.action') |> unname()
  if (length(NAs)){
    warning(sprintf('Missing values on observation(s) %s',
                    paste(NAs,collapse=', ')))
    # data <- data[seq_len(nrow(data))[-NAs],]
    data <- dplyr::slice(data, seq_len(nrow(data))[-NAs])
  }
  compare.method <- match.arg(compare.method)
  fit.surv0 <- survfit(update(base.model, .~1), data=data)
  ms <- inherits(fit.surv0, 'survfitms')
  
  # arm.var <- if (length(base.model[[3]]) > 1) {
  #   deparse(base.model[[3]][[2]])
  # } else {
  #     deparse(base.model[[3]])
  #   }
  arm.var <- formula.tools::rhs.vars(base.model)[[1]]
  if (!inherits(data[, arm.var], "factor")) data[, arm.var] <- factor(data[, arm.var])
  
  # result in entire population
  if (!missing(overall.model)){
    added.terms <- terms(overall.model) |> attr('term.labels')
    overall.model <- update(base.model,
                            as.formula(paste('. ~ . +',
                                             paste(added.terms,
                                                   collapse='+'))))
  } else overall.model <- base.model
  result <- sstable.survcomp(model = overall.model, data = data, time=time, reference.arm=reference.arm,
                             medsum = FALSE, digits = digits,
                             compare.method = compare.method, compare.args = compare.args,
                             p.compare = p.compare,
                             pdigits = 3, pcutoff = pcutoff, flextable = FALSE, ...)$table[,-1]
  result <- cbind(c("Subgroup", "", "All patients"), result, c("Test for heterogeneity", "p-value", ""))
  
  # Preparation of models and data
  subgroup.char <- all.vars(subgroup.model)
  
  #browser()
  for (k in 1:length(subgroup.char)){
    main.model <- update(base.model, as.formula(paste(". ~ . +", subgroup.char[k], sep = "")))
    ia.model <- update(base.model, as.formula(paste(". ~ . +`", arm.var,  "`*", subgroup.char[k], sep = "")))
    data$.subgroup.var <- data[, subgroup.char[k]]
    if (!inherits(data[, subgroup.char[k]], 'factor'))  data[, subgroup.char[k]] <- factor(data[, subgroup.char[k]])
    factor.levels <- levels(data[, subgroup.char[k]])
    compare.args.subgroup <- compare.args
    # Add interaction test for heterogeneity
    result <- rbind(result, "")
    result[nrow(result), 1] <- ifelse(is.null(attr(data[, subgroup.char[k]], "label")),
                                      subgroup.char[k], attr(data[, subgroup.char[k]], "label"))
    ia.pval <-
      if (compare.method == 'cox'){
        if (ms){
          etype <- compare.args$cause
          fg.ia <- survival::finegray(ia.model, data=data, etype=etype)
          fg.main <- survival::finegray(main.model, data=data, etype=etype)
          
          anova(
            survival::coxph(update(ia.model, Surv(fgstart,fgstop,fgstatus)~.),
                            data=fg.ia, weights=fgwt),
            survival::coxph(update(main.model, Surv(fgstart,fgstop,fgstatus)~.),
                            data=fg.main, weights=fgwt),
            test = "Chisq")[2, "Pr(>|Chi|)"]
        } else {
          anova(survival::coxph(ia.model, data = data),
                survival::coxph(main.model, data = data),
                test = "Chisq")[2, "Pr(>|Chi|)"]
        }
      } else {
        ia.args <- compare.args
        ia.args$add.prop.haz.test <- NULL
        # ia.args$formula <- ia.model
        ia.args$data <- data
        
        type <- if (is.null(ia.args$type)) 'diff' else ia.args$type
        ia.args$type <- NULL
        ia.args$link <- switch(type, "diff" = 'identity', "ratio" = 'log', 'lost.ratio' = 'log')
        ia.args$model.censoring <- if (is.null(ia.args$model.censoring)) 'stratified'
        ia.args$formula.censoring <-
          if (is.null(ia.args$formula.censoring)) update(main.model, NULL ~.) else
            update(ia.args$formula.censoring, as.formula(paste("NULL ~ . +", subgroup.char[k], sep = "")))
        mf <- model.frame(update(base.model, new = as.formula(paste0(". ~`", arm.var, '`'))), data = data)
        mf <- cbind(unclass(mf[,1]), mf[,2])
        # Get max of stop time, is the minimum of last observed time between two arms.
        # The stop time for right cens data is the first column ($time), for interval-cens data it is the second column ($stop)
        # With laziness it is coded as mf[, ncol(mf)-2] (the last two columns are the event status from Surv object, and the arm)
        # browser()
        minimax.time <- min(by(mf[, ncol(mf)-2], mf[, ncol(mf)], max))
        
        # tau for rmst is capped as minamax.time
        if (is.null(ia.args$time)) ia.args$time <- minimax.time
        else if (ia.args$time > minimax.time) {
          ia.args$time <- minimax.time
        }
        
        # Remove missing value from predictor bc this causes problem in censoring model
        NAs <- model.frame(main.model, data=data) |> attr('na.action') |> unname()
        if (length(NAs)){
          warning(sprintf('Missing %s on observation(s) %s',
                          subgroup.char[k],
                          paste(NAs,collapse=', ')))
          # data <- data[seq_len(nrow(data))[-NAs],]
          ia.args$data <- dplyr::slice(data, seq_len(nrow(data))[-NAs])
        }
        
        # browsere)
        # Perform multivariate wald test for interaction term
        tryCatch({
          fitter <- switch(compare.method,
                           'rmst' = eventglm::rmeanglm,
                           'cuminc' = eventglm::cumincglm)
          ia.fit <-
            do.call(fitter, append(ia.args, c(formula=ia.model)))
          main.terms <- colnames(model.matrix(main.model, data=ia.args$data))
          ia.terms <- coef(ia.fit)
          test.terms <- which(!names(ia.terms) %in% main.terms)
          test <- aod::wald.test(vcov(ia.fit), b=ia.terms, Terms=test.terms)
          test$result$chi2['P']
        },
        error=\(e) NA
        )
        
      }
    # browser()
    
    result[nrow(result), ncol(result)] <-
      if (is.na(ia.pval)) '-' else
        format.pval(ia.pval, digits = pdigits, eps = pcutoff)
    
    # Add results for each subgroup level
    for (j in 1:length(factor.levels)){
      result <- rbind(result, "")
      result[nrow(result), 1] <- paste("-", factor.levels[j])
      d.subgroup <- subset(data, .subgroup.var == factor.levels[j])
      # if (compare.method=='rmst'){
      #   compare.args$model.censoring <- if (is.null(compare.args$model.censoring)) 'stratified'
      #   compare.args$formula.censoring <-
      #     if (is.null(compare.args$formula.censoring)) update(main.model, NULL ~.) else
      #       update(compare.args$formula.censoring, as.formula(paste("NULL ~ . +", subgroup.char[k], sep = "")))
      # }
      # browser()
      result[nrow(result), 2:(ncol(result) - 1)] <- sstable.survcomp(model = base.model, data = d.subgroup, time=time,
                                                                     compare.method = compare.method, compare.args = compare.args,
                                                                     p.compare = p.compare,
                                                                     medsum = FALSE, digits = digits, pdigits = pdigits, pcutoff = pcutoff,
                                                                     flextable = FALSE, ...)$table[3,-1]
    }
  }
  ## footer
  compare.stat <- switch(compare.method,
                         cox = 'HR',
                         cuminc = if (is.null(compare.args$type)) 'Cumul.inc difference'
                         else switch(compare.args$type,
                                     diff = 'Cumul.inc difference', ratio = 'Cumul.inc ratio',
                                     stop('Illegal type for cumulative incidence comparison model')),
                         rmst = if (is.null(compare.args$type)) 'RMST difference'
                         else switch(compare.args$type,
                                     diff = 'RMST difference',
                                     lost.diff = 'RMTL difference',
                                     ratio = 'RMST ratio',
                                     lost.ratio = 'RMTL ratio',
                                     stop('Illegal type for RMST comparison model')))
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
    compare.note,
    paste(compare.stat, "and p value were based on", paste0(compare.name, '.')),
    sprintf(
      "Test for heterogeneity is an %s test for interaction between treatment effect and each subgroup in the survival model not including other variables.",
      if (compare.method == 'cox') 'Likelihood-ratio' else 'multivariate Wald'
    ),
    footer)
  
  # flextable
  if (!overall) result <- result[-3,]
  if (flextable) {
    requireNamespace("flextable")
    requireNamespace("officer")
    
    ## main table
    value <- result[-c(1, 2), ]
    tab <- flextable::flextable(as.data.frame(value))
    
    ## header
    header1 <- result[1, ]; header2 <- result[2, ]
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
    
    ## footer
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
    tab <- flextable::bg(tab, i = seq(from = 1, to = nrow(result[-c(1:2), ]), by = 2), j = 1:length(header1), bg = bg, part = "body")
    ### border
    tabbd <- officer::fp_border(color="black", width = 1.5)
    tab <- flextable::border_remove(tab)
    tab <- flextable::hline(tab, border = tabbd, part = "header")
    tab <- flextable::hline_top(tab, border = tabbd, part = "all")
    tab <- flextable::hline_bottom(tab, border = tabbd, part = "body")
    
  } else {
    tab <- list(table = result,
                footer = footer)
    class(tab$table) <- c('survcomp_tbl', 'ss_tbl', class(tab$table))
  }
  if (!flextable) class(tab) <- c('ss_survcomp','ss_obj')
  return(tab)
}

#' @export
print.ss_tbl <- function(sstable, pretty=getOption('ss_pretty.print', TRUE)){
  if (!pretty) return(print.default(sstable))
  tryCatch(
    {
      print(huxtable::as_hux(sstable) |> huxtable::theme_article())
    },
    error = function(e) {print.default(sstable)}
  )
  invisible(sstable)
}

