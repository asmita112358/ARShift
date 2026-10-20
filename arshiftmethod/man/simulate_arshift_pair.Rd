\name{simulate_arshift_pair}
\alias{simulate_arshift_pair}
\title{Generate a two-taxa marked point pattern for ARShift tests}
\usage{
simulate_arshift_pair(
  win = spatstat.geom::square(1),
  model = c("thomas", "cauchy"),
  dependant = FALSE,
  shift_radius = 0.2,
  kappa = c(10, 10),
  mu = c(20, 20),
  scale = c(0.05, 0.05),
  seed = NULL
)
}
\arguments{
\item{win}{Simulation window.}
\item{model}{Point process family.}
\item{dependant}{If \code{TRUE}, create dependency by shifting one pattern.}
\item{shift_radius}{Shift radius used when \code{dependant = TRUE}.}
\item{kappa}{Cluster intensity parameters for each taxon.}
\item{mu}{Offspring parameters for each taxon.}
\item{scale}{Cluster scale parameters for each taxon.}
\item{seed}{Optional random seed.}
}
\value{A marked \code{ppp} object with mark levels \code{1} and \code{2}.}
\description{Utility simulation generator for quick reproducible ARShift experiments.}
