test_that("Sims works", {
  expect_type(depthVsSample(cov.len=50,sam.len=100), "double")
})
