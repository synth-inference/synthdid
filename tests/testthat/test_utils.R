test_that("make.panel works as expected", {
  data("california_prop99")
  panel = california_prop99
  panel.shuffled = panel[sample(1:nrow(panel)), ]
  out = make.panel(panel)
  out.shuffled = make.panel(panel.shuffled)
  expect_equal(out, out.shuffled)

  # Removing one (unit, year) entry causes an unbalanced panel error
  expect_error(make.panel(panel[-10, ]), "Input `panel` must be a balaned panel data set.")

  # Duplicating some units causes an unbalanced panel error
  expect_error(make.panel(rbind(panel, panel[5:10, ])), "Input `panel` must be a balaned panel data set.")
})
