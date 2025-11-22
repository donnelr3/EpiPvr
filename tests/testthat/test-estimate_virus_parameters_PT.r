#########################################
## estimate_virus_parameters_PT ########
#########################################

test_that("PT virus parameter estimation output dimensions are appropriate", {
  # Use precomputed deterministic data to avoid long MCMC runs during CRAN checks

  test_data_path <- system.file("extdata/estimate_data_test_PT.rda", package = "EpiPvr")
  load(test_data_path)

  expect_equal(dim(estimate_data_test_PT$array1), c(1000,1,4))  # Expected shape
  expect_equal(length(estimate_data_test_PT$array2), 10)  # Expected shape
  expect_equal(length(estimate_data_test_PT$array3), 4)  # Expected shape
  expect_equal(length(estimate_data_test_PT$array4), 4)  # Expected shape

  expect_true(all(identical(dimnames(estimate_data_test_PT$array1)$variable,c("al[1]", "be[1]", "mu[1]", "lat[1]"))))
})


test_that("PT virus parameter estimation output content is appropriate", {
  # Use precomputed deterministic data to avoid long MCMC runs during CRAN checks
  test_data_path <- system.file("extdata/estimate_data_test_PT.rda", package = "EpiPvr")
  load(test_data_path)

  valmat1=rbind(estimate_data_test_PT$array2[[1]],estimate_data_test_PT$array2[[2]],estimate_data_test_PT$array2[[3]],estimate_data_test_PT$array2[[4]])
  valmat2=rbind(estimate_data_test_PT$array3[[1]],estimate_data_test_PT$array3[[2]],estimate_data_test_PT$array3[[3]],estimate_data_test_PT$array3[[4]])
  valmat3=rbind(estimate_data_test_PT$array4[[1]],estimate_data_test_PT$array4[[2]],estimate_data_test_PT$array4[[3]],estimate_data_test_PT$array4[[4]])

  exclude <- c("par1[1]", "par3[1]", "par6[1]", "lp__")

  expect_true(all(estimate_data_test_PT$array>=0))  # Expected output (estimates)
  expect_true(all(estimate_data_test_PT$array1>=0))  # Expected output (estimates)
  expect_true(all(valmat1[, !(valmat1[1, ] %in% exclude)]>=0))  # Expected output (assay data with estimated percentiles)
  expect_true(all(valmat2[, !(valmat2[1, ] %in% exclude)]>=0))  # Expected output (assay data with estimated percentiles)
  expect_true(all(valmat3[, !(valmat3[1, ] %in% exclude)]>=0))  # Expected output (assay data with estimated percentiles)
})

