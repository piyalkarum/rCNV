test_that("AlleleFreq works", {
  expect_is(allele.freq(hetTgen(readVCF(paste0(path.package("rCNV"), "/example.raw.vcf.gz")),"GT"),"ind"), class="data.frame")
})
