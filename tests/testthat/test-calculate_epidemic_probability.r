#####################################
## calculate_epidemic_probability ###
#####################################

# dependencies: solveInoculumStatesBP.R

calculate_epidemic_probability(numberInsects=3,start_intervl=(1/10),localParameters=rep(0.1,5),virusParameters=rep(0.1,3),thresh=10^(-7))

# does it perform its function? Done
test_that("calculate_epidemic_probability produces valid epidemic probabilities for SPT edge values", {
  qm_out <- calculate_epidemic_probability(6, (1/10), rep(0.1,5), c(0.0001, 0.0001, 60),10^(-7)) # unrealistically short retention; unrealistically slow acquisition and inocuation
  expect_true(is.numeric(qm_out))
  expect_true(all(qm_out >= 0 & qm_out <= 1)) # Probabilities should be in [0,1]
  qm_out <- calculate_epidemic_probability(6, (1/10), rep(0.1,5), c(1, 1, 60),10^(-7)) # unrealistically short retention; more typical acquisition and inocuation
  expect_true(is.numeric(qm_out))
  expect_true(all(qm_out >= 0 & qm_out <= 1)) # Probabilities should be in [0,1]
})

test_that("calculate_epidemic_probability produces valid epidemic probabilities for PT edge values", {
  qm_out <- calculate_epidemic_probability(1, (1/10), rep(0.1,5), c(0.0001, 0.0001, 0),10^(-7))  # zero rate retention loss; unrealistically slow acquisition and inocuation
  expect_true(is.numeric(qm_out))
  expect_true(all(qm_out >= 0 & qm_out <= 1)) # Probabilities should be in [0,1]
  qm_out <- calculate_epidemic_probability(1, (1/10), rep(0.1,5), c(1, 1, 0),10^(-7))  # zero rate retention loss; more typical acquisition and inocuation
  expect_true(is.numeric(qm_out))
  expect_true(all(qm_out >= 0 & qm_out <= 1)) # Probabilities should be in [0,1]
})

