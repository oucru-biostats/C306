## Quick sstable drawing helpers
## Author: trinhdhk
## Day first written: Sep 20 2019
## Latest build: Sep 27 2019
## Ver 0.1.0.092019

## Check legitibility
ss_legit <- function(sstable){
  # requireNamespace('tibble')
  if (!is.matrix(sstable) & !is.data.frame(sstable)) stop('sstable must be of type matrix or data frame.')
  return(as.matrix(sstable))
}
#' Set template for sstable
#'
#' @description This function set the template for sstable
#' @param sstable a data frame following sstable's grammar
#' @param template
#' An accepted template for sstable: must be either 'baseline', 'survcomp', or 'ae'.
#'
#' If NA, the existing template in the sstable is kept as-is. It no template built-in, will return a no-template sstable instead.
#' @return a matrix of class ss_tbl
#' @export
ss_template <- function(sstable, template =  c('baseline', 'survcomp', 'ae')){
  sstable <- ss_legit(sstable)
  ss.class <- class(sstable)
  if (missing(template)) template <- NA
  if (!is.na(template)) template <- match.arg(template)
  if (is.na(template)){
    template <- c('baseline', 'survcomp', 'ae')
    template.match <- match(template, ss.class)
    template.match <- template.match[!is.na(template.match)]
    if (length(template.match) > 1) template.match <- template.match[[1]]
    template <- template[template.match]
  }

  template <- paste0(template, '_tbl')
  class(sstable) <- unique(c(template, 'ss_tbl', 'matrix', class(sstable)))
  return(sstable)
}

#' Designation of header rows for custom sstable
#'
#' @description This function set the designated rows as header rows of a sstable
#' @param sstable a data frame following sstable's grammar
#' @param rows a numeric vector
#' @return a matrix of class ss_tbl
#' @export
ss_header <- function(sstable, rows) {
  sstable <- ss_legit(sstable)
  if (!length(rownames(sstable))) rownames(sstable) <- seq_len(nrow(sstable))
  rownames(sstable)[rows] <- paste("header", seq_along(rows), sep=".")
  ss_table <- ss_template(sstable, template = NA)
  sstable
}

#' Designation of body part for custom sstable
#'
#' @description This function set the designated rows as body part of a sstable
#' @param sstable a data frame following sstable's grammar
#' @param rows a numeric vector
#' @return a matrix of class ss_tbl
#' @export
ss_body <- function(sstable, rows){
  sstable <- ss_legit(sstable)
  if (!length(rownames(sstable))) rownames(sstable) <- seq_len(nrow(sstable))
  rownames(sstable)[rows] <- paste("body", seq_along(rows), sep=".")
  ss_table <- ss_template(sstable, template = NA)
  sstable
}

#' Designation of section title rows for custom sstable
#'
#' @description This function set the designated rows as section title rows of an sstable
#' @param sstable a data frame following sstable's grammar
#' @param rows a numeric vector
#' @return a matrix of class ss_tbl
#' @export
ss_section <- function(sstable, rows){
  sstable <- ss_legit(sstable)
  if (!length(rownames(sstable))) rownames(sstable) <- seq_len(nrow(sstable))
  rownames(sstable)[rows] <- paste("section", seq_along(rows), sep=".")
  ss_table <- ss_template(sstable, template = NA)
  sstable
}

#' Quick format a sstable in preparation for create flextable/huxtable
#'
#' @description This function is the combination of ss_header, ss_body, and ss_section.
#' @param sstable a data frame following sstable's grammar
#' @param header a numeric vector that will be passed to ss_header
#' @param section a numeric vector that will be passed to ss_section
#' @param body a numeric vector that will be passed to ss_body
#' @param template
#' An accepted template for sstable: must be either 'baseline', 'survcomp', or 'ae'.
#'
#' If NA, the existing template in the sstable is kept as-is. It no template built-in, will return a no-template sstable instead.
#' @param .guess a logical value. If TRUE, the function will try to guess to which part each row belongs.
#' @return a matrix of class ss_tbl
#' @export
ss_format <- function(sstable, header = c(), section = c(), body = c(), template =  c('baseline', 'survcomp', 'ae'), .guess= TRUE){
  if (missing(template)) template <- NA
  if (!is.na(template)) template <- match.arg(template)
  sstable <- ss_template(sstable, template = template)
  if (.guess) {
    no.guess.in.table <- which(grepl(
      '(header)|(section)|(body)',
      rownames(sstable),
      perl = TRUE
    ))
    no.guess <- union(no.guess.in.table, na.omit(as.numeric(c(header, section, body))))
    # filter out from guess what have been defined
    guess <- ss_guess_format(sstable)
    h <- which(guess == 'header')
    s <- which(guess == 'section')
    b <- which(guess == 'body')
    # filling the undecided row with auto guess
    header <- unique(c(h[!h %in% no.guess], header))
    # if (length(colnames(sstable))) header <- c('colnames', header)
    section <- unique(c(s[!s %in% no.guess], section))
    body <- unique(c(b[!b %in% no.guess], body))
  }
  if (length(header)) sstable <- ss_header(sstable, rows = header)
  if (length(body)) sstable <- ss_body(sstable, rows = body)
  if (length(section)) sstable <- ss_section(sstable, rows = section)
  if (sum(sapply(c('header', 'section', 'body'), grepl, x=rownames(sstable), fixed=TRUE)) != nrow(sstable))
    stop('Some row are not classified as header, section, or body. Perhaps you want to set .guess to TRUE?')
  class(sstable) <- c('formatted_sstable', class(sstable))
  sstable
}

# An incomplete algorithm to guess the format of each sstable row.
ss_guess_format <- function(sstable){
  UseMethod('ss_guess_format')
}

ss_guess_format.default <- function(sstable){
  guess <- sapply(seq_len(nrow(sstable)),
                  function(i){
                    r <- sstable[i,]
                    if (isTRUE(grepl("header", rownames(sstable)[i]))) return("header")
                    if (all(is.na(suppressWarnings(as.numeric(r)))) && i <= 2) return("header")
                    if (sum(r == '', na.rm=TRUE) == ncol(sstable)-1) return("section")
                    return("body")
                  })
  return(guess)
}

ss_guess_format.ae_tbl <- function(sstable){
  guess <- sapply(seq_len(nrow(sstable)),
                  function(i){
                    if (i <= 2) return('header')
                    r <- sstable[i, ]
                    if (sum(r == '', na.rm=TRUE) >= 2) return('section')
                    return('body')
                  })
  return(guess)
}

ss_guess_format.baseline_tbl <- function(sstable){
  guess <- c(rep('header',2), rep('body', nrow(sstable)-2))
  return(guess)
}

ss_guess_format.survcomp_tbl <- function(sstable){
  guess <- c(rep('header',2), rep('body', nrow(sstable)-2))
  return(guess)
}

ss_guess_format.summary_tbl <- function(sstable){
  guess <- c('header', rep('body', nrow(sstable)-1))
  return(guess)
}

#' Create summary table using flextable package.
#'
#' @description This function generate a flextable from a sstable.
#' @param sstable a data frame following sstable's grammar
#' @param footer a character vector each of which is the footnote of the flextable
#' @param bg a character string that defines stripped background color of the flextable
#' @param ... additional parameters that will be passed to ss_format if the sstable has yet to be formatted.
#' @seealso \link[flextable]{flextable}
#' @return an object of class flextable
#' @export
ss_flextable <- function(sstable, ...){
  UseMethod('ss_flextable')
}

#' @rdname ss_flextable
#' @method ss_flextable list
#' @param add_footer additional footer lines to be appended to object footers
#' @export
ss_flextable.list <- function(sstable, add_footer = NULL,...){
  ss_flextable(sstable$table, footer = c(sstable$footer, add_footer), ...)
}

#' @rdname ss_flextable
#' @method ss_flextable ss_obj
#' @export
ss_flextable.ss_obj <- ss_flextable.list

#' @rdname ss_flextable
#' @export
ss_flextable.default <- function(sstable, footer = NULL, bg = "#F2EFEE", ...){
  requireNamespace("flextable")
  requireNamespace("officer")
  if (!inherits(sstable, 'formatted_sstable')) sstable <- ss_format(sstable, ..., .guess = TRUE)
  header <- which(grepl("header", rownames(sstable)))
  body <- which(grepl("body", rownames(sstable)))
  section <-  which(grepl("section", rownames(sstable)))

  colnames(sstable) <- NULL
  sstable <- as.data.frame(sstable, stringsAsFactors = FALSE)
  col.even <- ncol(sstable) %% 2 == 0

  ## headers
  ss.header <- sstable[header,]
  names(ss.header) <- colnames(sstable)
  # browser()

  ss.header2 <-
      lapply(seq_along(header),
             function(i){
               h <- ss.header[i, ]
               if (col.even) col <- seq_len((ncol(sstable)-1))[-1]
               else col <- seq_len(ncol(sstable))[-1]
               right_fill <- h[[max(col)]] != "" & h[[1]] == ""

               h <- sapply(seq_along(h),
                      function(j){
                        if (h[[j]] == "") {
                          if (j == ncol(sstable) && i>1) return(ss.header[i-1,j])
                          if (j == 1){
                            if (i > 1) return(ss.header[i-1,j]) #lag 1 row
                            return(h[[j]]) #return ""
                          }
                          if (right_fill) return(h[[j+1]])
                          return(h[[j-1]])
                        }
                        return(h[[j]])
                      })
               names(h) <- paste0('V',seq_len(ncol(sstable)))
               return(h)
             }
    )
  # browser()

  ## Create the flextable
  ft <- flextable::flextable(sstable[-header,, drop = FALSE])
  # browser()

  ## header format
  ft <- flextable::set_header_labels(ft, values = ss.header2[[1]])
  if (length(ss.header2) > 1)
    for (i in seq_along(ss.header2)[-1])
      ft <- flextable::add_header_row(ft, values = ss.header2[[i]], top = F)

  ft <- flextable::merge_h(ft, i = seq_len(length(ss.header2)-1), part = "header")
  ft <- flextable::merge_v(ft, part = "header")

  ## footer format
  for (k in (seq_along(footer))) {
    ft <- flextable::add_footer(ft, V1 = footer[k], top = FALSE)
    ft <- flextable::merge_at(ft, i = k, j = seq_len(ncol(sstable)), part = "footer")
  }

  ## section format
  for (k in section){

    ### merging cells that from the left if the whole row is empty
    if (all(sstable[k, -1] == '') %in% TRUE)
      ft <- flextable::merge_at(ft, i = k, j = seq_len(ncol(sstable)), part = 'body')
    else ft <- flextable::merge_h(ft, i = k , part = 'body')

    ft <- flextable::bold(ft, i = k-length(header), j = 1, part = 'body')
  }

  ## format flextable
  ## width
  ft <- flextable::autofit(ft)
  ### alignment
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")

  ft_sstheme(ft, bg = "#F2EFEE" )
}


#' Flextable theming for sstable
#'
#' @description Theming flextable for sstable
#' @param ft flextable
#' @export
ft_sstheme <- function(ft, bg = "#F2EFEE"){
  ### faces of header
  ft <- flextable::bold(ft, part = "header")
   ### background
  ft <- flextable::bg(ft, i = seq(from = 1, to = nrow(ft$body$dataset), by = 2), j = seq_len(ncol(ft$body$dataset)),
                      bg = bg, part = "body")
  ### border
  tabbd <- officer::fp_border(color="black", width = 1.5)
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline(ft, border = tabbd, part = "header")
  ft <- flextable::hline_top(ft, border = tabbd, part = "all")
  ft <- flextable::hline_bottom(ft, border = tabbd, part = "body")
  return(ft)
}

#' Create summary table using huxtable package.
#'
#' @description This function generate a huxtable from a sstable.
#' @param sstable a data frame following sstable's grammar
#' @param footer a character vector each of which is the footnote of the flextable
#' @param caption a string containing table caption. Default is NULL
#' @param caption_pos
#' A length-one character vector,
#' one of "top", "bottom", "topleft", "topcenter", "topright", "bottomleft", "bottomcenter", "bottomright".
#' Default is "bottomcenter".
#'
#' See also \link[huxtable]{caption_pos}
#' @param bg a character vector that defines background color of the flextable. If length(bg) >= 2, the table will have stripe background, otherwise plain.
#' @param border_width a number that defines huxtable border width
#' @param border_color a character string that defines huxtable border color
#' @param wrap a logical value. Default is FALSE. If TRUE, long texts would be wrapped, long lines would be broken.
#' This applies to the whole table. Please use huxtable::set_wrap instead if you want a cell-wise approach.
#' See also \link[huxtable]{wrap}
#' @param ... additional parameters that will be passed to ss_format if the sstable has yet to be formatted.
#' @return an object of class huxtable
#' @seealso \link[huxtable]{huxtable}
#' @export
ss_huxtable <- function(sstable,...){
  UseMethod('ss_huxtable')
}

#' @rdname ss_huxtable
#' @param add_footer additional footer lines to be appended to object footers
#' @export
ss_huxtable.list <- function(sstable, add_footer = NULL,...){
  ss_huxtable(sstable$table, footer = c(sstable$footer, add_footer), ...)
}

#' @rdname ss_huxtable
#' @export
ss_huxtable.default <- function(sstable, footer = NULL,
                        caption = NULL, caption_pos = c("top", "bottom", "topleft", "topcenter", "topright",
                                                        "bottomleft", "bottomcenter", "bottomright"),
                        bg = c(grey(.95), 'white'), border_width=0.8, border_color = grey(.75), wrap = FALSE,...){
  requireNamespace('huxtable')
  if (missing(caption_pos)) caption_pos <- 'bottomcenter' else caption_pos <- match.arg(caption_pos)
  if (!inherits(sstable, 'formatted_sstable')) sstable <- ss_format(sstable, ..., .guess = TRUE)
  header <- which(grepl("header", rownames(sstable)))
  body <- which(grepl("body", rownames(sstable)))
  section <-  which(grepl("section", rownames(sstable)))

  colnames(sstable) <- NULL
  sstable <- as.data.frame(sstable, stringsAsFactors = FALSE)
  col.even <- ncol(sstable) %% 2 == 0

  ## headers pre-processing
  ss.header <- sstable[header,]
  names(ss.header) <- colnames(sstable)
  # browser()

  ss.header2 <-
    lapply(seq_along(header),
           function(i){
             h <- ss.header[i, ]
             if (col.even) col <- 2:(ncol(sstable)-1)
             else col <- 2:ncol(sstable)
             right_fill <- h[[max(col)]] != "" & h[[2]] == ""

             h <- sapply(seq_along(h),
                         function(j){
                           if (h[[j]] == "") {
                             if (j == ncol(sstable) && i>1) return(ss.header[i-1,j])
                             if (j == 1){
                               if (i > 1) return(ss.header[i-1,j]) #lag 1 row
                               return(h[[j]]) #return ""
                             }
                             if (right_fill) return(h[[j+1]])
                             return(h[[j-1]])
                           }
                           return(h[[j]])
                         })
             return(h)
           }
    )

  sstable.headless <- sstable[-header,]
  colnames(sstable.headless) <- ss.header2[[length(ss.header2)]]

  ## create table
  ht <- huxtable::as_hux(sstable.headless, add_colnames = TRUE, add_rownames = FALSE)

  ## insert header
  if (nrow(ss.header)>1){
    for (i in seq_along(ss.header2)[-length(ss.header2)]){
      ht <- huxtable::insert_row(ht, ss.header2[[i]])
    }
  }

  ## format header
  # browser()

  # merge cols that have same values
  for (i in seq_along(ss.header2)){
    which.start <- 1
    for (j in seq_along(ss.header2[[i]])){
      if (j > 1){
        if (ss.header2[[i]][[j]] != ss.header2[[i]][[j-1]])  which.start <- c(which.start, j)
      }
      ht <- huxtable::set_align(ht, i, j, value = 'center')
    }

    which.start <- sort(unique(c(which.start, ncol(sstable)+1)))

    for (k in seq_along(which.start)){
      if (k < length(which.start) & which.start[k+1] - which.start[k] > 1){
        ht <- huxtable::merge_cells(ht, i, which.start[k]:(which.start[k+1]-1))
      }
    }
  }
  # merge rows that have same values
  for (i in seq_len(ncol(sstable))){
    which.start <- 1
    for (j in seq_along(ss.header2)){
      if (j > 1){
        if (ss.header2[[j]][[i]] != ss.header2[[j-1]][[i]]) which.start <- c(which.start, j)
      }
    }

    which.start <- sort(unique(c(which.start, length(ss.header2)+1)))

    for (k in seq_along(which.start)){
      if (k < length(which.start) & which.start[k+1] - which.start[k] > 1){
        ht <- huxtable::merge_cells(ht, which.start[k]:(which.start[k+1]-1), i)
      }
    }
  }

  ## format sections
  for (i in section){
    ht <- huxtable::set_bold(ht, i, huxtable::everywhere, value = TRUE)
    ### merging cells that from the left if the whole row is empty
    if (all(sstable[i, -1] == ''))
      ht <- huxtable::merge_cells(ht, i, huxtable::everywhere)
  }

  ## footer
  for (i in seq_along(footer)){
    ht <- huxtable::add_footnote(ht, footer[i], border = if (i>1) 0 else border_width)
  }

  ## caption
  if (length(caption)){
    huxtable::caption(ht) <- caption
    huxtable::caption_pos(ht) <- caption_pos
  }

  ## format huxtable
  ht <- huxtable::set_width(ht, 1)

  ht <- ht_theme_markdown(ht, header_rows = header, header_cols = NULL,
                          border_width = border_width,
                          border_color = border_color,
                          bg = bg)

  ## format wrapping
  huxtable::wrap(ht) <- wrap

  return(ht)
}


#' A stripe theme for huxtable object
#'
#' @description This function provides a markdown stripe theme for huxtable object
#' @param ht an object of class huxtable
#' @param header_rows a numeric vector that delimits the header zone.
#' @param bg a character vector that defines background color of the flextable. If length(bg) >= 2, the table will have stripe background, otherwise plain.
#' @param border_width a number that defines huxtable border width
#' @param border_color a character string that defines huxtable border color
#' @return an object of class huxtable
#' @export
ht_theme_markdown <- function(ht, header_rows = 1:2, header_cols = NULL,
                              border_width=0.8, border_color = grey(0.75), bg = c(grey(.95),'white')){
  huxtable::top_border(ht)[1, ] <- border_width
  huxtable::bottom_border(ht)[nrow(ht),] <- border_width

  for (header_row in header_rows){
    huxtable::bold(ht)[header_row, ] <- TRUE
  }

  huxtable::bottom_border(ht)[max(header_rows), ] <- border_width
  huxtable::top_border(ht)[max(header_rows)+1, ] <- border_width

  for (header_col in header_cols){
    huxtable::bold(ht)[, header_col] <- TRUE
  }
  ht <- huxtable::set_all_border_colors(ht, border_color)
  if (length(bg > 1)) ht <- huxtable::map_background_color(ht, do.call(huxtable::by_rows, as.list(bg)))
  else huxtable::background_color(ht) <- bg
  ht
}

#' A kable-esque theme for huxtable object
#'
#' @description This function provides a kable-esque theme for huxtable object
#' @param ht an object of class huxtable
#' @param header_rows a numeric vector that delimits the header zone.
#' @param bg a character vector that defines background color of the flextable. If length(bg) >= 2, the table will have stripe background, otherwise plain.
#' @param border_width a number that defines huxtable border width
#' @param border_color a character string that defines huxtable border color
#' @return an object of class huxtable
#' @export
ht_theme_kable <- function(ht, header_rows = 1:2, header_cols = NULL,
                           border_width = 1, border_color = '#dddddd', bg = c('white', '#f9f9f9')){
  huxtable::top_border(ht)[1, ] <- 0
  huxtable::bottom_border(ht)[nrow(ht),] <- border_width

  for (header_row in header_rows){
    huxtable::bold(ht)[header_row, ] <- TRUE
    for (col in seq_len(ncol(ht))){
      if (!is.na(ht[header_row, col]) & stringr::str_trim(ht[header_row, col]) != ''){
        huxtable::bottom_border(ht)[header_row, col] <- border_width
      }
    }
  }
  # browser()

  huxtable::bottom_border(ht)[max(header_rows):nrow(ht), ] <- border_width
  huxtable::top_border(ht)[max(header_rows)+1, ] <- border_width*2

  for (header_col in header_cols){
    huxtable::bold(ht)[, header_col] <- TRUE
  }

  ht <- huxtable::set_all_border_colors(ht, border_color)
  if (length(bg > 1)){
    ht <- huxtable::map_background_color(ht, do.call(huxtable::by_rows, append(as.list(bg), c(from = max(header_rows)+1))))
    huxtable::background_color(ht)[header_rows,] <- bg[[1]]
  } else huxtable::background_color(ht) <- bg
  ht
}

#' Coerce an object to sstable
#' @description A function to coerce objects to a sstable
#' @param x An object, usually a named list of length 2 whose names are 'table' and 'footer', or a data.frame/matrix (optionally with attribute "footer")
#' @param flextable logical value specifying whether to return a flextable. Default is FALSE
#' @param include_footnote logical value specifying whether to include footnote in the output. Default is FALSE
#' @param ... additional parameters passed to \link{ss_flextable}
#' @return A matrix of class ss_tbl if flextable == FALSE, otherwise a flextable
#' @export
#'
as_sstable <- function(x,...){
  UseMethod('as_sstable')
}

#' @rdname as_sstable
#' @export
as_sstable.default <- function(x, flextable = FALSE, ...){
  out <- list()
  out$table <- as.matrix(x$table)
  out$table <- rbind(colnames(out$table), out$table)
  out$table <- cbind(rownames(out$table), out$table)
  colnames(out$table) <- rownames(out$table) <- NULL
  if (length(attr(x, 'footer'))) out$footer <- attr(x, 'footer')

  class(out$table) <- c('ss_tbl', 'matrix')
  if (flextable) return(ss_flextable(out, ...))
  return(out)
}

#' @rdname as_sstable
#' @export
as_sstable.list <- function(x, flextable = FALSE, ...){
  out <- list()

  if (is.null(x$header))
    return(as_sstable.default(x, flextable = flextable, ...))
  header <- if (inherits(x$header, 'matrix')) x$header else do.call(rbind, x$header)
  rownames(header) = paste0('header', seq_len(nrow(header)))

  if (inherits(x$body, 'list')){
    body <- lapply(seq_along(x$body), function(i){
      cont <- x$body[[i]]
      title <- names(x$body)[[i]]
      if (is.null(title)) title <- ''
      title <- rep(title, ncol(cont))
      sec <- rbind(title, cont)
      rownames(sec) = c(paste0('section', i), paste0('body',i,  seq_len(nrow(cont))))
      if (title == '' & i == 1) sec <- sec[-1,]
      sec
    })
    body <- do.call(rbind, body)
  } else {
      body <- as.matrix(x$body)
      rownames(body) <- paste0('body', seq_len(nrow(body)))
  }

  out$table <- rbind(header, body)
  out$footer <- x$footer
  class(out$table) <- c('ss_tbl', 'matrix')

  if (flextable) return(ss_flextable(out, ...))
  return(out)
}

#' @rdname as_sstable
#' @export
as_sstable.logist_summary <- function(x, include_footnote = TRUE, flextable = FALSE, ...){
  out <- list()
  out$table <- as.matrix(x)
  out$table <- rbind(colnames(out$table), out$table)
  out$table <- cbind(rownames(out$table), out$table)
  colnames(out$table) <- rownames(out$table) <- NULL
  out$footer <- attr(x, 'footer')
  class(out$table) <- c('summary_tbl', 'ss_tbl', 'matrix')

  if (flextable){
    logist_summary.sstable <- ss_flextable(out$table, footer = out$footer, ...)
    return(logist_summary.sstable)
  }

  if (include_footnote) return(out)
  out$table
}

#' @rdname as_sstable
#' @param include_footsnote logical value specifying whether to include footnote in the output. Default is FALSE
#' @export
as_sstable.overlap_summary <- function(x, include_footnote = TRUE, flextable = FALSE, ...){
  out <- list()
  out$table <- as.matrix(x)
  out$table <- rbind(colnames(out$table), out$table)
  out$table <- cbind(rownames(out$table), out$table)
  colnames(out$table) <- rownames(out$table) <- NULL
  out$footer <- attr(x, 'footer')
  class(out$table) <- c('summary_tbl', 'ss_tbl', 'matrix')

  if (flextable){
    overlap_summary.sstable <- ss_flextable(out$table, footer = out$footer, ...)
    return(overlap_summary.sstable)
  }

  if (include_footnote) return(out)
  out$table
}

#' @rdname as_sstable
#' @export
as_sstable.subgroup_logist_summary <- function(x, include_footnote = TRUE, flextable = FALSE, ...)
{
  tables <- sapply(x, as_sstable, include_footnote = FALSE, ..., simplify = FALSE)
  footer <- attr(x[[1]], 'footer')
  first_header <- c('', tables[[1]][1,-1])
  tables <- lapply(tables, `[`, -1,)
  headers <- matrix(rep('', length(tables) * ncol(tables[[1]])), ncol = ncol(tables[[1]]))
  headers[,1] <- paste(attr(x, 'base_var'), '=', names(tables))
  # browser()
  table <- do.call(rbind,
                   lapply(seq_along(tables),
                          function(i) rbind(headers[i,], tables[[i]])))
  table <- rbind(first_header, table)
  sstable <- ss_format(table, header = 1, section = (seq_len(length(tables)) - 1)*(nrow(tables[[1]])+1)+2)
  if (flextable) return(ss_flextable(sstable, footer = footer))
  if (include_footnote) return(list(table = sstable, footer = footer))
  sstable
}


#' Row-binding two sstables
#' @description
#' Row-binding two sstables
#' @param tbl1,tbl2, two object of class sstables of the same template. Only export for technical use.
#' @param footer [NULL] default is the joint footer of two ss_tbl
#' @param ... args passed to the downstream method
#' @export
rbind.ss_obj <- function(tbl1, tbl2, footer=NULL, ...){
  if (!identical(class(tbl1$table), class(tbl2$table)))
    stop('tbl1 and tbl2 must be of the same class.')
  if (ncol(tbl1$table) != ncol(tbl2$table))
    stop('Unmatched number of columns between tbl1 and tbl2.')
  NextMethod('rbind')
}

._do_rbind <- function(tbl1, tbl2, header=1){
  new_tbl <- list()
  new_tbl$table <- rbind(tbl1$table, tbl2$table[-header,])
  new_tbl$footer <- tbl1$footer
  class(new_tbl$table) <- class(tbl1$table)
  new_tbl
}

#' @rdname rbind.ss_obj
#' @method rbind ss_ae
#' @export
rbind.ss_ae <- function(tbl1, tbl2, footer=NULL){
  if (!all.equal(tbl1$footer, tbl2$footer) & is.null(footer))
    stop('Two footers mismatched. Perhaps two tbl are using different config? \n
         To ignore this, set a specific footer.')
  out <- ._do_rbind(tbl1, tbl2, header=1:3)
  if (!is.null(footer)) out$footer <- footer
  out
}

#' @rdname rbind.ss_obj
#' @method rbind ss_survcomp
#' @export
rbind.ss_survcomp <- function(tbl1, tbl2, footer=NULL){
  if (!all.equal(tbl1$footer, tbl2$footer) & is.null(footer))
    stop('Two footers mismatched. Perhaps two tbl are using different config? \n
         To ignore this, set a specific footer.')
  out <- ._do_rbind(tbl1, tbl2, header=1:3)
  if (!is.null(footer)) out$footer <- footer
  out
}

#' @rdname rbind.ss_obj
#' @method rbind ss_baseline
#' @export
rbind.ss_baseline <- function(tbl1, tbl2, footer=NULL){
  if (!all.equal(tbl1$footer, tbl2$footer) & is.null(footer))
    stop('Two footers mismatched. Perhaps two tbl are using different config? \n
         To ignore this, set a specific footer.')
  out <- ._do_rbind(tbl1, tbl2, header=1)
  if (!is.null(footer)) out$footer <- footer
  out
}
