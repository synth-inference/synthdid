test_that("panel.matrices works as expected", {
  data("california_prop99")
  panel = california_prop99
  panel.shuffled = panel[sample(1:nrow(panel)), ]
  covariates = as.data.frame(matrix(runif(nrow(panel) * 5), nrow(panel), 5))
  panel.aug = cbind(covariates, panel)
  panel.aug = panel.aug[sample(1:ncol(panel.aug))]

  # Expected output
  out = panel.matrices(panel)

  # Is the same if input data rows are shuffled
  out.shuffled = panel.matrices(panel.shuffled)
  expect_equal(out, out.shuffled)

  # In the California data set the last unit in the reshape data is California
  expect_equal(row.names(out$Y)[nrow(out$Y)], "California")

  # Mixed column input works
  out.int.string <- panel.matrices(panel, unit = "State", time = "Year", outcome = 3, treatment = "treated")
  expect_equal(out, out.int.string)

  # Passing in a data.frame with more columns
  out.aug <- panel.matrices(panel.aug, unit = "State", time = "Year", outcome = "PacksPerCapita", treatment = "treated")
  expect_equal(out, out.aug)
  out.aug2 <- panel.matrices(cbind(covariates, panel), unit = 6, time = 7, outcome = "PacksPerCapita", treatment = 9)
  expect_equal(out, out.aug2)

  # A year column in Date format is fine
  panel.date = panel
  panel.date$Year = as.Date(panel.date$Year, origin = "1970-12-30")
  out.date = panel.matrices(panel.date)
  expect_equal(unname(out.date$Y), unname(out$Y))
  expect_equal(out.date$N0, out$N0)
  expect_equal(out.date$T0, out$T0)
  expect_equal(unname(out.date$W), unname(out$W))

  # Removing one (unit, year) entry causes an unbalanced panel error
  expect_error(panel.matrices(panel[-10, ]), "Input `panel` must be a balanced panel: it must have an observation for every unit at every time.")

  # Duplicating some units causes an unbalanced panel error
  expect_error(panel.matrices(rbind(panel, panel[5:10, ])), "Input `panel` must be a balanced panel: it must have an observation for every unit at every time.")

  # If we make Kansas a treated state, it is the last row in the reshaped matrix
  panel[panel$State =="Kansas" & panel$Year >=1989, "treated"] = 1
  out = panel.matrices(panel)
  expect_equal(row.names(out$Y)[(nrow(out$Y)-1):nrow(out$Y)], c("California", "Kansas"))

  # Treating it one year earlier causes an error
  panel[panel$State =="Kansas" & panel$Year >=1988, "treated"] = 1
  expect_error(panel.matrices(panel), "The package cannot use this data. Treatment adoption is not simultaneous.")
})
