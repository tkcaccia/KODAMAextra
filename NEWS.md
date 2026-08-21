# KODAMAextra 1.2.3

- Removed `passing.message()`; use `KODAMA::passing.message()`.
- Removed the Louvain and Leiden resolution wrappers; use
  `fastEmbedR::graph_cluster()` with `method = "louvain"` or `"leiden"`.
- Removed incomplete, unexported development helpers.
- Reduced the package to its maintained plotting, rendering, trajectory, and
  annotation utilities and removed unused compiled/parallel dependencies.
