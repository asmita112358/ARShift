\name{ARShift}
\alias{ARShift}
\title{Compute ARShift area-corrected p-value}
\usage{
ARShift(
  data,
  base_taxa = 1,
  shift_taxa = 2,
  r = NULL,
  n_perm = 499,
  k = 1.2,
  type = "inhom",
  lite = FALSE,
  thin_prob = 0.2,
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
\item{type}{Either \code{"inhom"} or \code{"hom"}.}
\item{lite}{Whether to thin points before testing.}
\item{thin_prob}{Thinning probability when \code{lite = TRUE}.}
\item{jump_radius}{Shift radius for non-toroidal shifts.}
\item{correction}{K-function correction mode.}
\item{workers}{Number of workers for permutation attempts.}
}
\value{A list containing \code{pval_Kcross_area} and permutation diagnostics.}
\description{Extracted implementation of the ARShift area-corrected p-value calculation from \code{test_spatial_association_pmclapply}.}
