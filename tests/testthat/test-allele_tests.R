test_that("AlleleFreq works", {
  expect_type(allele.freq(hetTgen(readVCF(vcf.file.path=paste0(path.package("rCNV"), "/example.raw.vcf.gz")),"GT"),"ind"), "data.frame")
})
