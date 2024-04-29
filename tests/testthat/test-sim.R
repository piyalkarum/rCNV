test_that("simulation works", {
  expect_type(sim.als(n=50,nrun=100), "list")
})

