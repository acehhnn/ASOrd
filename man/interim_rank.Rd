% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interim_testing.R
\name{interim_rank}
\alias{interim_rank}
\title{Use rank sum statistics in the interim analysis}
\usage{
interim_rank(K, k, y1, y0, A, alpha = 0.05, t = NULL, func = "obf", gamma = 0)
}
\arguments{
\item{K}{number of tests including the final analysis}

\item{k}{order of the test among all tests}

\item{y1}{outcomes of currently involved members in the new treatment group}

\item{y0}{outcomes of currently involved members in the control group}

\item{A}{number of levels of the outcomes}

\item{alpha}{overall Type I error}

\item{t}{a sequence of fractions in the interim analysis}

\item{func}{choosing an alpha spending function}

\item{gamma}{parameter used in gamma family alpha spending function}
}
\value{
Results of hypothesis testing
}
\description{
Use rank sum statistics in the interim analysis
}
