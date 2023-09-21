
# BuscarR

BuscarR, short for Biased Urn based Stopping Criteria for technology
Assisted Review in R, is a package implementing the stopping criteria
described in Callaghan, Max, and Finn Müller-Hansen. 2020. “Statistical
Stopping Criteria for Automated Screening in Systematic Reviews.”
Systematic Reviews, September. <https://doi.org/10.21203/rs.2.18218/v2>.

The stopping criteria offers a *reliable* method to decide when to stop
screening during machine learning-prioritised screening for a systematic
review, or for document review in ediscovery.

Such a method is *vital* in order to use machine learning safely (by
estimating the risk of missing documents). Our method describes
statistical confidence levels that any given level of recall has *not*
been achieved. As such, it offers a robust data-driven method to decide
when to stop, and a clear and transparent way of communicating the risks
of stopping at any given point. Read `vignette("stopping-criteria")` for
more details on how this works in the package, or check out the original
[paper](https://doi.org/10.21203/rs.2.18218/v2) for a full description
of the criteria, including theory and evaluations.

## Installation

<div class=".pkgdown-devel">

``` r
# Install development version from GitHub
devtools::install_github("mcallaghan/buscarR")
```

</div>

## Usage

BuscarR requires a dataframe with a row for each document in a query,
and two columns: `relevant`, and `seen`. `relevant` should contain a 1,
where a documents has been included by a human, or a 0 where it has been
excluded by a human. Any documents not yet screened should have an NA.
All documents screened by a human should contain a 1 in the `seen`
column, and documents not yet seen by a human should contain a 0. To get
a p-score for a null hypothesis that a given recall target has been
missed, use the `calculate_h0` function on your dataset.

``` r
library(buscarR)
calculate_h0(df, recall_target=0.95)
```

If the p score is below 1 minus your confidence level, you can reject
the null hypothesis that the target has been missed, and stop screening.
`vignette("stopping-criteria")` contains more examples, explanations,
and additional functions.
