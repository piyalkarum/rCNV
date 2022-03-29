test_that("readVCF works", {
  expect_is(readVCF(paste0(path.package("rCNV"), "/example.raw.vcf.gz")), class="list")
})
