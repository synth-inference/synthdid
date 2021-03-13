test_that("panel.matrices works as expected", {
  data("california_prop99")
  panel = california_prop99
  panel.shuffled = panel[sample(1:nrow(panel)), ]
  out = panel.matrices(panel)
  out.shuffled = panel.matrices(panel.shuffled)
  expect_equal(out, out.shuffled)

  # Removing one (unit, year) entry causes an unbalanced panel error
  expect_error(panel.matrices(panel[-10, ]), "Input `panel` must be a balanced panel data set.")

  # Duplicating some units causes an unbalanced panel error
  expect_error(panel.matrices(rbind(panel, panel[5:10, ])), "Input `panel` must be a balanced panel data set.")
})
