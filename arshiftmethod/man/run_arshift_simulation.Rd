\name{run_arshift_simulation}
\alias{run_arshift_simulation}
\title{Run repeated ARShift simulations and summarize rejection rates}
\usage{
run_arshift_simulation(
  n_sim = 20,
  generator = NULL,
  alpha = 0.05,
  seed = 123,
  test_args = list(base_taxa = 1, shift_taxa = 2, n_perm = 39, type = "inhom", workers = 1)
)
}
\arguments{
\item{n_sim}{Number of simulation replicates.}
\item{generator}{Function that returns one marked \code{ppp} dataset per replicate.}
\item{alpha}{Significance threshold for rejection summaries.}
\item{seed}{Optional random seed.}
\item{test_args}{Named arguments passed to \code{test_spatial_association_pmclapply}.}
}
\value{A data frame of replicate p-values and rejection flags with summary attributes.}
\description{Convenience runner for reproducible small simulation studies.}
