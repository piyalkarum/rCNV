test_that("Infor works", {
  expect_type(allele.info(ADtable,x.norm=ADnorm), "list")
})
