% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_mult.R
\name{summarise_f}
\alias{summarise_f}
\alias{summarise_f.data.frame}
\alias{summarise_f.data.table}
\alias{summarise_f.tbl}
\alias{summarise_f.grouped_df}
\title{Reduce multiple values by summarisation functions.}
\usage{
summarise_f(.data, ...)

\method{summarise_f}{data.frame}(.data, ...)

\method{summarise_f}{data.table}(.data, ...)

\method{summarise_f}{tbl}(.data, ...)

\method{summarise_f}{grouped_df}(.data, ...)
}
\arguments{
\item{.data}{A data.frame, data.table, or tbl}

\item{...}{Summarisation functions.
Each of which must return a named list, whose names stand for new variables' names.
Anonymous functions must be wrapped in parentheses.}
}
\value{
a data.frame, data.table, or a tbl
}
\description{
Create variables summarising exisiting variables, by using summarisation function.
}
\examples{
df <- data.frame(baz = 1:10, foo = c(rep(1, 5), rep(2, 5)), bar =4:13)
df \%>\%
   dplyr::group_by(foo) \%>\%
   summarise_f((function(baz) list(a = cumsum(baz), b = baz*2))(baz), rangemisc::overlap_collapse(bar, baz))
}
