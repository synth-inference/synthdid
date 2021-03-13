test_that("panel.matrices works as expected", {
  data("california_prop99")
  panel = california_prop99
  panel.shuffled = panel[sample(1:nrow(panel)), ]
  out = panel.matrices(panel)
  out.shuffled = panel.matrices(panel.shuffled)
  expect_equal(out, out.shuffled)
  # In the California data set the last unit in the reshape data is California
  expect_equal(row.names(out$Y)[nrow(out$Y)], "California")

  # Removing one (unit, year) entry causes an unbalanced panel error
  expect_error(panel.matrices(panel[-10, ]), "Input `panel` must be a balanced panel data set.")

  # Duplicating some units causes an unbalanced panel error
  expect_error(panel.matrices(rbind(panel, panel[5:10, ])), "Input `panel` must be a balanced panel data set.")

  # If we make Kansas a treated state, it is the last row in the reshaped matrix
  panel[panel$State =="Kansas" & panel$Year >=1989, "treated"] = 1
  out = panel.matrices(panel)
  expect_equal(row.names(out$Y)[(nrow(out$Y)-1):nrow(out$Y)], c("California", "Kansas"))

  # Treating it one year earlier causes an error
  panel[panel$State =="Kansas" & panel$Year >=1988, "treated"] = 1
  expect_error(panel.matrices(panel), "The package cannot use this data. Treatment adoption is not simultaneous.")
})
