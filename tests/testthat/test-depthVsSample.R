test_that("Sims works", {
  expect_is(depthVsSample(cov.len=50,sam.len=100), class="matrix")
})
