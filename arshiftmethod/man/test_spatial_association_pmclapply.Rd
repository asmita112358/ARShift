\name{test_spatial_association_pmclapply}
\alias{test_spatial_association_pmclapply}
\title{Run ARShift spatial association test with permutation resampling}
\usage{
test_spatial_association_pmclapply(
  data,
  base_taxa = 1,
  shift_taxa = 2,
  r = NULL,
  n_perm = 499,
  k = 1.2,
  bw = "silverman",
  type = "inhom",
  lite = FALSE,
  thin_prob = 0.2,
  include_RL = FALSE,
  jump_radius = NULL,
  correction = "isotropic",
  workers = max(1L, parallel::detectCores() - 1L)
)
}
\arguments{
\item{data}{A marked \code{ppp} point pattern.}
\item{base_taxa}{Base mark label.}
\item{shift_taxa}{Shifted mark label.}
\item{r}{Distance vector. If \code{NULL}, chosen from window geometry.}
\item{n_perm}{Number of retained valid permutations.}
\item{k}{Multiplier controlling attempted permutations.}
\item{bw}{Bandwidth rule name or numeric value for variance correction.}
\item{type}{Either \code{"inhom"} or \code{"hom"}.}
\item{lite}{Whether to thin points before testing.}
\item{thin_prob}{Thinning probability when \code{lite = TRUE}.}
\item{include_RL}{Whether to include random-labelling baseline p-value.}
\item{jump_radius}{Shift radius for non-toroidal shifts.}
\item{correction}{K-function correction mode.}
\item{workers}{Number of workers for permutation attempts.}
}
\value{A named list with p-values and permutation diagnostics.}
\description{Extracted and hardened implementation of the ARShift fast wrapper method.}
