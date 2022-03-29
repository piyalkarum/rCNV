test_that("simulation works", {
  expect_is(sim.als(n=50,nrun=100), class="list")
})

