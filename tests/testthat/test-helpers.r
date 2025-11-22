
  ##########################
  ## AP_insect_simulator ###
  ##########################

  # dependencies: inoc_durtn_calculator

  # does it perform its function? test return vals
  test_that("AP_insect_simulator produces valid SPT insect inoculation data", {
    fixedComp_tmp=c(100,100,100)
    smarkpams_tmp=c(0.1/60,0.1/60,0.1/60,0.1/60)
    expect_error(AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'SPT'),
                 'input error! SPT virus should have NA gamma')
  })

  # does it perform its function? test return vals
  test_that("AP_insect_simulator produces valid SPT insect inoculation data", {
    fixedComp_tmp=c(100,100,100)
    smarkpams_tmp=c(0.1/60,0.1/60,NA,0.1/60)
    ins_out=AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'SPT')
    expect_true(ins_out[1] %in% c(0,1)) # check binary
    expect_true(ins_out[2]%%1==0) # check integer
  })


  test_that("AP_insect_simulator produces valid PT insect inoculation data", {
    fixedComp_tmp=c(100,100,100)
    smarkpams_tmp=c(0.1/60,0.1/60,0.1/60,0.1/60)
    ins_out=AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'PT')
    expect_true(ins_out[1] %in% c(0,1)) # check binary
    expect_true(ins_out[2]==0) # check zero as do not expect re-accquisition for PT virus
  })

  # test edge cases=0 t_A>T_A, t_L>T_A
  test_that("AP_insect_simulator produces valid insect inoculation data", {
    fixedComp_tmp=c(0,100,100)                      #test zero acquisition period
    smarkpams_tmp=c(0.1/60,0.1/60,0.1/60,0.1/60)
    ins_out=AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'PT')
    expect_equal(ins_out[1],0) # confirm zero
    expect_equal(ins_out[2],0) # confirm zero
    fixedComp_tmp=c(100,100,0)                      #test zero inoculation period
    smarkpams_tmp=c(0.1/60,0.1/60,0.1/60,0.1/60)
    ins_out=AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'PT')
    expect_equal(ins_out[1],0) # confirm zero
    expect_true(is.numeric(ins_out[2])) # confirm numeric
    fixedComp_tmp=c(100,100,100)
    smarkpams_tmp=c(0,0.1/60,0.1/60,0.1/60)          #test zero acquisition
    ins_out=suppressWarnings(AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'PT')) # suppress known zero rate warnings as tested separately in test_that("AP_insect_simulator problem exponential rates inform user in gen_exp()", {
    expect_equal(ins_out[1],0) # confirm zero
    expect_equal(ins_out[2],0) # confirm zero
    smarkpams_tmp=c(0.1/60,0.1/60,0.1/60,0)          #test zero recovery
    ins_out=suppressWarnings(AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'PT'))  # suppress known zero rate warnings as tested separately in test_that("AP_insect_simulator problem exponential rates inform user in gen_exp()", {
    expect_true(ins_out[1] %in% c(0,1)) # check binary
    expect_true(is.numeric(ins_out[2])) # confirm numeric
    smarkpams_tmp=c(0.1/60,0,0.1/60,0.1/60)          #test zero inoculation
    ins_out=suppressWarnings(AP_insect_simulator(fixedComp_tmp, smarkpams_tmp,0,'PT'))  # suppress known zero rate warnings as tested separately in test_that("AP_insect_simulator problem exponential rates inform user in gen_exp()", {
    expect_equal(ins_out[1],0) # confirm zero
    expect_true(is.numeric(ins_out[2])) # confirm numeric
  })

  # test difficult random sampling cases
  test_that("AP_insect_simulator problem exponential rates inform user in gen_exp()", {
    expect_warning(gen_exp(1,0,'PT'), regexp="infinite time sample")
    expect_warning(gen_exp(1,1/0,'PT'), regexp="zero time sample")
    expect_error(gen_exp(1,NA,'PT'), regexp="missing exponential rate")
    expect_error(gen_exp(1,NaN,'PT'), regexp="undefined exponential rate")
  })


  #####################################
  ## inoc_durtn_calculator ############
  #####################################

  # does inoc_durtn_calculator perform its function?
  test_that("inoc_durtn_calculator produces valid duration data", {
    # Now call the nested function directly
    fixedComp_tmp = c(100, 100, 100)
    smarkpams_tmp = c(90, 90, 90)
    dur_out = inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp, 0)
    # Check that the output is numeric
    expect_true(is.numeric(dur_out), info = "Output should be numeric")
    # Check that the output is non-negative
    expect_true(all(dur_out >= 0), info = "Output should be non-negative")
  })

  # test edge cases=0 t_A>T_A, t_L>T_A
  test_that("inoc_durtn_calculator produces appropriate durations edge cases", {
    fixedComp_tmp=c(100,100,100)
    smarkpams_tmp1=c(110,90,90)
    smarkpams_tmp2=c(90,120,90)
    smarkpams_tmp3=c(90,90,10)
    smarkpams_tmp2b=c(90,220,90)
    expect_message(inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp1,1), regexp="no acquisition in acquisition period")
    expect_message(inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp2,1), regexp="no latency pass in latent period")
    expect_message(inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp3,1), regexp="loss of pathogen prior to inoculation period")
    # several cases should lead to zero duration
    expect_equal(inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp1,0),0)
    expect_equal(inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp2b,0),0)
    expect_equal(inoc_durtn_calculator(fixedComp_tmp, smarkpams_tmp3,0),0)
  })


  #####################################
  ## solveInoculumStatesBP.R ##########
  #####################################



