test_that("Infor works", {
  expect_is(allele.info(ADtable,x.norm=ADnorm), class="data.frame")
})
