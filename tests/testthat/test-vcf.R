test_that("stats works", {
  expect_type(vcf.stat(readVCF(vcf.file.path=paste0(path.package("rCNV"), "/example.raw.vcf.gz"))), "data.frame")
})
