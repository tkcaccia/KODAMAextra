test_that("the public API contains only maintained utilities", {
  expect_setequal(
    getNamespaceExports("KODAMAextra"),
    c("new_trajectory", "plot_slide", "read_annotations", "volume_rendering")
  )
  expect_false(exists("passing.message", asNamespace("KODAMAextra"), inherits = FALSE))
  expect_false(exists("louvain", asNamespace("KODAMAextra"), inherits = FALSE))
  expect_false(exists("leiden", asNamespace("KODAMAextra"), inherits = FALSE))
})

test_that("plot_slide returns a ggplot without attaching ggplot2", {
  xy <- cbind(x = rep(1:4, 2), y = rep(1:2, each = 4))
  plot <- plot_slide(
    xy, rep(c("A", "B"), each = 4), factor(rep(c("one", "two"), 4))
  )
  expect_s3_class(plot, "ggplot")
})
