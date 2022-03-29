test_that("stats works", {
  expect_is(vcf.stat(readVCF(paste0(path.package("rCNV"), "/example.raw.vcf.gz"))), class="data.frame")
})
