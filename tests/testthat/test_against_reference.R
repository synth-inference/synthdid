tol=.03
min.decrease = 1e-6
max.iter = 1e5
if(TRUE) {
test_that("synthdid point estimate agrees with the reference implementation", {
    if (!requireNamespace("CVXR", quietly = TRUE)) { stop("Package CVXR must be installed to run this test.") }
    data(california_prop99)
    setup = panel.matrices(california_prop99)
    expect_equal(c(synthdid_estimate(setup$Y, setup$N0, setup$T0, min.decrease=min.decrease, max.iter=max.iter)),
		 c(synthdid:::synthdid.reference(setup$Y, setup$N0, setup$T0)), tol=tol)
})

test_that("sc point estimate agrees with the reference implementation", {
    if (!requireNamespace("CVXR", quietly = TRUE)) { stop("Package CVXR must be installed to run this test.") }
    data(california_prop99)
    setup = panel.matrices(california_prop99)
    expect_equal(c(sc_estimate(setup$Y, setup$N0, setup$T0, min.decrease=min.decrease, max.iter=max.iter)),
		 c(synthdid:::sc.reference(setup$Y, setup$N0, setup$T0)), tol=tol)
})

test_that("did point estimate agrees with the reference implementation", {
    if (!requireNamespace("CVXR", quietly = TRUE)) { stop("Package CVXR must be installed to run this test.") }
    data(california_prop99)
    setup = panel.matrices(california_prop99)
    expect_equal(c(did_estimate(setup$Y, setup$N0, setup$T0)),
		 c(synthdid:::did.reference(setup$Y, setup$N0, setup$T0)), tol=tol)
})
}
