test_that("readVCF works", {
  expect_type(readVCF(vcf.file.path=paste0(path.package("rCNV"), "/example.raw.vcf.gz")), "list")
})
