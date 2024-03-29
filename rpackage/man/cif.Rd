% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cif.R
\name{cif}
\alias{cif}
\title{Cumulative incidence function (CIF) with inverse probability weighting (IPW)}
\usage{
cif(
  times,
  event,
  eoi = 1,
  group = NULL,
  weight = NULL,
  alpha = 0.05,
  entry = NULL
)
}
\arguments{
\item{times}{time to event (any-cause) or censoring}

\item{event}{integer denoting the event type, where 0=censored. If there is only 1 event type (values take 0,1), results are equivalent to a non-competing risk setting}

\item{eoi}{the event of interest, which should be an event value (defaults to 1)}

\item{group}{exposure/treatment group. If left blank, function calculates the CIF in the overall sample.}

\item{weight}{user can specify weights to adjust for confounding, i.e. inverse probability weights, overlap weights, or survey weights. If left blank, function calculates unadjusted results.}

\item{alpha}{level of significance (defaults to 0.05).}

\item{entry}{optional argument to allow delayed entry times (i.e., age as the time scale)}
}
\value{
IPW (or unadjusted) CIF in each group, confidence intervals, and timepoints for plotting cumulative incidence curves.
}
\description{
Cumulative incidence function (CIF) with inverse probability weighting (IPW)
}
