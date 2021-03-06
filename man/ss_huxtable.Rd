% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssformat.R
\name{ss_huxtable}
\alias{ss_huxtable}
\title{Create summary table using huxtable package.}
\usage{
ss_huxtable(
  sstable,
  footer = NULL,
  caption = NULL,
  caption_pos = c("top", "bottom", "topleft", "topcenter", "topright", "bottomleft",
    "bottomcenter", "bottomright"),
  bg = c(grey(0.95), "white"),
  border_width = 0.8,
  border_color = grey(0.75),
  wrap = FALSE,
  ...
)
}
\arguments{
\item{sstable}{a data frame following sstable's grammar}

\item{footer}{a character vector each of which is the footnote of the flextable}

\item{caption}{a string containing table caption. Default is NULL}

\item{caption_pos}{A length-one character vector,
one of "top", "bottom", "topleft", "topcenter", "topright", "bottomleft", "bottomcenter", "bottomright".
Default is "bottomcenter".

See also \link[huxtable]{caption_pos}}

\item{bg}{a character vector that defines background color of the flextable. If length(bg) >= 2, the table will have stripe background, otherwise plain.}

\item{border_width}{a number that defines huxtable border width}

\item{border_color}{a character string that defines huxtable border color}

\item{wrap}{a logical value. Default is FALSE. If TRUE, long texts would be wrapped, long lines would be broken.
This applies to the whole table. Please use huxtable::set_wrap instead if you want a cell-wise approach.
See also \link[huxtable]{wrap}}

\item{...}{additional parameters that will be passed to ss_format if the sstable has yet to be formatted.}
}
\value{
an object of class huxtable
}
\description{
This function generate a huxtable from a sstable.
}
\seealso{
\link[huxtable]{huxtable}
}
