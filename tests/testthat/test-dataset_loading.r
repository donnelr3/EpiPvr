#############################################
#### ap_data_sim_SPT ########################
#############################################

test_that("SPT data has the correct basic dimensions", {
  data("ap_data_sim_SPT", package = "EpiPvr")
  expect_equal(length(attributes(ap_data_sim_SPT)), 4)
  expect_equal(length(ap_data_sim_SPT), 5)
})

test_that("SPT data has the correct internal categories", {
  data("ap_data_sim_SPT", package = "EpiPvr")
  expect_true(grepl(names(ap_data_sim_SPT)[1],"d_AAP"))
  expect_true(grepl(names(ap_data_sim_SPT)[2],"d_IAP"))
  expect_true(grepl(names(ap_data_sim_SPT)[3],"d_durations"))
  expect_true(grepl(names(ap_data_sim_SPT)[4],"d_vectorspp"))
  expect_true(grepl(names(ap_data_sim_SPT)[5],"d_virusType"))
})


test_that("SPT data assay structures are correct size", {
  data("ap_data_sim_SPT", package = "EpiPvr")
  tm1=dim(ap_data_sim_SPT$d_AAP)
  expect_equal(tm1[1],3)
  expect_true(tm1[2]>0)
  tm2=dim(ap_data_sim_SPT$d_IAP)
  expect_equal(tm2[1],3)
  expect_true(tm2[2]>0)
  expect_equal(dim(ap_data_sim_SPT$d_durations),c(2,2))
  expect_equal(length(ap_data_sim_SPT$d_vectorspp),1)
  expect_equal(length(ap_data_sim_SPT$d_virusType),1)
})

test_that("SPT data non-attributes has appropriate content", {
  data("ap_data_sim_SPT", package = "EpiPvr")
  expect_true(grepl(rownames(ap_data_sim_SPT$d_AAP)[1],"T_vec"))
  expect_true(grepl(rownames(ap_data_sim_SPT$d_AAP)[2],"R_vec"))
  expect_true(grepl(rownames(ap_data_sim_SPT$d_AAP)[3],"I_vec"))
  expect_true(grepl(rownames(ap_data_sim_SPT$d_IAP)[1],"T_vec"))
  expect_true(grepl(rownames(ap_data_sim_SPT$d_IAP)[2],"R_vec"))
  expect_true(grepl(rownames(ap_data_sim_SPT$d_IAP)[3],"I_vec"))
  expect_true(all(ap_data_sim_SPT$d_AAP>=0))
  expect_true(all(ap_data_sim_SPT$d_IAP>=0))
  expect_equal(diag(ap_data_sim_SPT$d_duration),c(-1,-1))
  expect_true(sum(ap_data_sim_SPT$d_duration)>-2)
  expect_true(ap_data_sim_SPT$d_vectorspp>0)
})

test_that("SPT data attributes has appropriate names", {
  data("ap_data_sim_SPT", package = "EpiPvr")
  expect_true(any(grepl("alpha",names(attributes(ap_data_sim_SPT)))))
  expect_true(any(grepl("mu",names(attributes(ap_data_sim_SPT)))))
  expect_true(any(grepl("beta",names(attributes(ap_data_sim_SPT)))))
})

test_that("SPT data attributes has appropriate content", {
  data("ap_data_sim_SPT", package = "EpiPvr")
  expect_true(attributes(ap_data_sim_SPT)$"alpha">0)
  expect_true(attributes(ap_data_sim_SPT)$"mu">0)
  expect_true(attributes(ap_data_sim_SPT)$"beta">0)
})

#############################################
#### ap_data_sim_PT #########################
#############################################

test_that("PT data has the correct basic dimensions", {
  data("ap_data_sim_PT", package = "EpiPvr")
  expect_equal(length(attributes(ap_data_sim_PT)), 5)
  expect_equal(length(ap_data_sim_PT), 6)
})

test_that("PT data has the correct internal dimensions", {
  data("ap_data_sim_PT", package = "EpiPvr")
  expect_true(grepl(names(ap_data_sim_PT)[1],"d_AAP"))
  expect_true(grepl(names(ap_data_sim_PT)[2],"d_LAP"))
  expect_true(grepl(names(ap_data_sim_PT)[3],"d_IAP"))
  expect_true(grepl(names(ap_data_sim_PT)[4],"d_durations"))
  expect_true(grepl(names(ap_data_sim_PT)[5],"d_vectorspp"))
  expect_true(grepl(names(ap_data_sim_PT)[6],"d_virusType"))
})


test_that("PT data assay structures are correct size", {
  data("ap_data_sim_PT", package = "EpiPvr")
  expect_equal(dim(ap_data_sim_PT$d_AAP)[1],3)
  expect_true(dim(ap_data_sim_PT$d_AAP)[2]>0)
  expect_equal(dim(ap_data_sim_PT$d_LAP)[1],3)
  expect_true(dim(ap_data_sim_PT$d_LAP)[2]>0)
  expect_equal(dim(ap_data_sim_PT$d_IAP)[1],3)
  expect_true(dim(ap_data_sim_PT$d_IAP)[2]>0)
  expect_equal(dim(ap_data_sim_PT$d_durations),c(3,3))
  expect_equal(length(ap_data_sim_PT$d_vectorspp),1)
  expect_equal(length(ap_data_sim_PT$d_virusType),1)
})

test_that("PT data non-attributes has appropriate content", {
  data("ap_data_sim_PT", package = "EpiPvr")
  expect_true(grepl(rownames(ap_data_sim_PT$d_AAP)[1],"T_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_AAP)[2],"R_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_AAP)[3],"I_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_LAP)[1],"T_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_LAP)[2],"R_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_LAP)[3],"I_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_IAP)[1],"T_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_IAP)[2],"R_vec"))
  expect_true(grepl(rownames(ap_data_sim_PT$d_IAP)[3],"I_vec"))
  expect_true(all(ap_data_sim_PT$d_AAP>=0))
  expect_true(all(ap_data_sim_PT$d_LAP>=0))
  expect_true(all(ap_data_sim_PT$d_IAP>=0))
  expect_equal(diag(ap_data_sim_PT$d_duration),c(-1,-1,-1))
  expect_true(sum(ap_data_sim_PT$d_duration)>-3)
  expect_true(ap_data_sim_PT$d_vectorspp>0)
})

test_that("PT data attributes has appropriate names", {
  data("ap_data_sim_PT", package = "EpiPvr")
  expect_true(any(grepl("alpha",names(attributes(ap_data_sim_PT)))))
  expect_true(any(grepl("mu",names(attributes(ap_data_sim_PT)))))
  expect_true(any(grepl("gamma",names(attributes(ap_data_sim_PT)))))
  expect_true(any(grepl("beta",names(attributes(ap_data_sim_PT)))))
})

test_that("PT data attributes has appropriate content", {
  data("ap_data_sim_PT", package = "EpiPvr")
  expect_true(attributes(ap_data_sim_PT)$"alpha">0)
  expect_true(attributes(ap_data_sim_PT)$"mu">0)
  expect_true(attributes(ap_data_sim_PT)$"beta">0)
  expect_true(attributes(ap_data_sim_PT)$"gamma">0)
})


