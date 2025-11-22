#########################################
## estimate_virus_parameters_SPT ########
#########################################

test_that("SPT virus parameter estimation output dimensions are appropriate", {
  # Use precomputed deterministic data to avoid long MCMC runs during CRAN checks
  test_data_path <- system.file("extdata/estimate_data_test_SPT.rda", package = "EpiPvr")
  load(test_data_path)

  expect_equal(dim(estimate_data_test_SPT$array1), c(1000,1,3))  # Expected shape
  expect_equal(length(estimate_data_test_SPT$array2), 10)  # Expected shape
  expect_equal(length(estimate_data_test_SPT$array3), 4)  # Expected shape

  expect_true(all(identical(dimnames(estimate_data_test_SPT$array1)$variable,c("al[1]", "be[1]", "mu[1]"))))
})

test_that("SPT virus parameter estimation output content is appropriate", {
  # Use precomputed deterministic data to avoid long MCMC runs during CRAN checks
  test_data_path <- system.file("extdata/estimate_data_test_SPT.rda", package = "EpiPvr")
  load(test_data_path)

  valmat1=rbind(estimate_data_test_SPT$array2[[1]],estimate_data_test_SPT$array2[[2]],estimate_data_test_SPT$array2[[3]],estimate_data_test_SPT$array2[[4]])
  valmat2=rbind(estimate_data_test_SPT$array3[[1]],estimate_data_test_SPT$array3[[2]],estimate_data_test_SPT$array3[[3]],estimate_data_test_SPT$array3[[4]])

  exclude <- c("XX", "lp__")

  expect_true(all(estimate_data_test_SPT$array>=0))  # Expected output (estimates)
  expect_true(all(estimate_data_test_SPT$array1>=0))  # Expected output (estimates)
  expect_true(all(valmat1[, !(valmat1[1, ] %in% exclude)]>=0))  # Expected output (assay data with estimated percentiles)
  expect_true(all(valmat2>=0))  # Expected output (assay data with estimated percentiles)
})

