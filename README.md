# KODAMAextra

`KODAMAextra` contains focused visualization and trajectory helpers that
complement KODAMA:

- `plot_slide()` draws labels across one or more tissue sections.
- `volume_rendering()` renders labelled three-dimensional regions.
- `new_trajectory()` interactively defines and summarizes a trajectory.
- `read_annotations()` reads the package's supported annotation CSV format.

Numerical functionality now lives in the package that owns it:

- use `KODAMA::passing.message()` for spatial message passing;
- use `fastEmbedR::graph_cluster(method = "louvain")` for Louvain;
- use `fastEmbedR::graph_cluster(method = "leiden")` for Leiden;
- use `fastEmbedR::graph_cluster(method = "walktrap")` for Walktrap.

Install the current development version with:

```r
remotes::install_github("tkcaccia/KODAMAextra")
```
