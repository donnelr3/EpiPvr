##########################
### AP_assay_simulator ###
##########################

# dependencies: AP_insect_simulator inoc_durtn_calculator.R

# does it perform its function? test return vals
test_that("AP_assay_simulator produces valid insect SPT inoculation data", {
  T_vec_in=seq(0,100,by=10)
  R_vec_in=rep(10,length(T_vec_in))
  I_vec_in=rep(0,length(T_vec_in))
  assayinput=rbind(T_vec_in,R_vec_in,I_vec_in)
  assayoutput=AP_assay_simulator(assayinput,
                                 c(100,100,100),
                                 10,rep(0.1,4),0,'PT')
  expect_equal(dim(assayinput),dim(assayoutput)) # check structure
  tmpdiff=assayoutput-assayinput
  expect_equal(sum(tmpdiff[1:2,]==0)*1,dim(assayinput)[2]*2) # check only TPIs (final row) changed
  expect_true(all(tmpdiff-floor(tmpdiff)==0)) # check integer
  expect_true(all(tmpdiff>=0)) # check positive
  expect_true(all((tmpdiff[3,]<=assayinput[2,]))) # check TPIs less than numReps
})


# test edge cases=0 t_A>T_A, t_L>T_A
test_that("AP_insect_simulator produces appropriate response to unusual rate input", {
  T_vec_in=seq(0,100,by=10)
  R_vec_in=rep(10,length(T_vec_in))
  I_vec_in=rep(0,length(T_vec_in))
  assayinput=rbind(T_vec_in,R_vec_in,I_vec_in)
  ###
  expect_warning(AP_assay_simulator(assayinput,
                                    c(100,100,100),
                                    10, c(0,rep(0.1,3)),0,'PT'), regexp="infinite time sample")
  expect_true(is.numeric(suppressWarnings(AP_assay_simulator(assayinput,
                                                             c(100,100,100),
                                                             10, c(0,rep(0.1,3)),0,'PT')))) # check numeric anyway suppressing warning as tested explicitly in line above
  ###
  expect_error(AP_assay_simulator(assayinput,
                                    c(100,100,100),
                                    10, c(0,rep(0.1,3)),0,'SPT'),'input error! SPT virus should have NA gamma')
  expect_true(is.numeric(suppressWarnings(AP_assay_simulator(assayinput,
                                                             c(100,100,100),
                                                             10, c(0,0.1,NA,0.1),0,'SPT')))) # check numeric anyway suppressing warning as tested explicitly in line above
  ###
  expect_error(AP_assay_simulator(assayinput,
                                    c(100,100,100),
                                    10, c(rep(0.1,3),0),0,'SPT'), regexp="input error! SPT virus should have NA gamma")
  expect_warning(AP_assay_simulator(assayinput,
                                    c(100,100,100),
                                    10, c(rep(0.1,3),0),0,'PT'), regexp="infinite time sample")
  expect_true(is.numeric(suppressWarnings(AP_assay_simulator(assayinput,
                                                             c(100,100,100),
                                                             10, c(0.1,0.1,NA,0),0,'SPT')))) # check numeric anyway suppressing warning as tested explicitly in line above
  ###
  expect_warning(AP_assay_simulator(assayinput,
                                    c(100,100,100),
                                    10, c(0.1,0,NA,0.1),0,'SPT'), regexp="infinite time sample")
  expect_true(is.numeric(suppressWarnings(AP_assay_simulator(assayinput,
                                                             c(100,100,100),
                                                             10, c(0.1,0,NA,0.1),0,'SPT')))) # check numeric anyway suppressing warning as tested explicitly in line above
})



# test bad and large input
test_that("AP_insect_simulator produces appropriate response to unusual durations input", {
  T_vec_in=seq(0,100,by=10)
  R_vec_in=rep(10,length(T_vec_in))
  I_vec_in=rep(0,length(T_vec_in))
  assayinput=rbind(T_vec_in,R_vec_in,I_vec_in)
  assayinputBad=assayinput
  assayinputBad[1,1]=0/0
  expect_error(AP_assay_simulator(assayinputBad,
                                  c(100,100,100),
                                  10, rep(0.1,4),0,'SPT'), regexp="missing or NaN values")
  expect_error(AP_assay_simulator(assayinput,
                                  c(0/0,100,100),
                                  10, rep(0.1,4),0,'SPT'), regexp="missing or NaN values")
  expect_error(AP_assay_simulator(assayinput,
                                  c(100,100,100),
                                  10, c(0/0,rep(0.1,3)),0,'SPT'), regexp="input error! SPT virus should have NA gamma")
  expect_error(AP_assay_simulator(assayinput,
                                  c(100,100,100),
                                  10, c(0/0,rep(0.1,3)),0,'PT'), regexp="undefined exponential rate")
  ###
  T_vec_in=seq(0,1000,by=1000)
  R_vec_in=rep(10,length(T_vec_in))
  I_vec_in=rep(0,length(T_vec_in))
  assayinput=rbind(T_vec_in,R_vec_in,I_vec_in)
  assayoutput=AP_assay_simulator(assayinput,
                                 c(100,100,100),
                                 10, c(0.1,0.1,NA,0.1),0,'SPT')
  expect_equal(dim(assayinput),dim(assayoutput)) # check structure
  tmpdiff=assayoutput-assayinput
  expect_equal(sum((tmpdiff[1:2,]==0)*1),dim(assayinput)[2]*2) # check only TPIs (final row) changed
  expect_true((all(tmpdiff[3,]%%1==0))&&(all(tmpdiff[3,]>=0))) # check positive integer
  expect_true(all(tmpdiff[3,]<=assayinput[2,])) # check TPIs less than numReps
})

