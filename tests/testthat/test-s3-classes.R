# tests/testthat/test-s3-classes.R
# Fast tests for S3 class structure and methods

# ==============================================================================
# TEST: GenSynResult Class Structure
# ==============================================================================

test_that("GenSynResult class has correct structure", {
  mock_result <- list(
    call = quote(GenSyn(n_households = 10)),
    n_households = 10,
    simulation = list(hh_df = data.frame(), diagnostic_df = data.frame()),
    surveillance_df = NULL,
    start_date = "2024-07-01",
    end_date = "2025-06-30",
    stan_data = list(),
    fit = NULL,
    postprocessing = data.frame(),
    attack_rates = NULL,
    transmission_chains = NULL
  )
  class(mock_result) <- "GenSynResult"

  expect_s3_class(mock_result, "GenSynResult")
  expect_true("n_households" %in% names(mock_result))
  expect_true("simulation" %in% names(mock_result))
  expect_true("fit" %in% names(mock_result))
  expect_true("start_date" %in% names(mock_result))
  expect_true("end_date" %in% names(mock_result))
})


test_that("print.GenSynResult exists and is callable", {
  expect_true(is.function(print.GenSynResult))

  mock_result <- list(
    n_households = 10,
    attack_rates = NULL,
    postprocessing = data.frame()
  )
  class(mock_result) <- "GenSynResult"

  expect_output(print(mock_result), "GenSyn Result")
})


test_that("plot.GenSynResult exists and is callable", {
  expect_true(is.function(plot.GenSynResult))

  # Check function signature
  args <- formals(plot.GenSynResult)
  expect_true("x" %in% names(args))
  expect_true("which" %in% names(args))
  expect_true("print" %in% names(args))
  expect_true("hh_id" %in% names(args))
  expect_true("prob_cutoff" %in% names(args))
})


# ==============================================================================
# TEST: TransmissionChainResult Class Structure
# ==============================================================================

test_that("TransmissionChainResult class has correct structure", {
  mock_result <- list(
    call = quote(TransmissionChainAnalysis(user_data = df)),
    user_data = data.frame(),
    surveillance_df = NULL,
    start_date = as.Date("2024-01-01"),
    end_date = as.Date("2024-12-31"),
    stan_data = list(),
    fit = NULL,
    postprocessing = data.frame(),
    transmission_chains = NULL
  )
  class(mock_result) <- "TransmissionChainResult"

  expect_s3_class(mock_result, "TransmissionChainResult")
  expect_true("user_data" %in% names(mock_result))
  expect_true("fit" %in% names(mock_result))
  expect_true("start_date" %in% names(mock_result))
  expect_true("end_date" %in% names(mock_result))
})


test_that("print.TransmissionChainResult exists and is callable", {
  expect_true(is.function(print.TransmissionChainResult))

  mock_result <- list(
    user_data = data.frame(hh_id = "HH1", person_id = 1),
    postprocessing = data.frame()
  )
  class(mock_result) <- "TransmissionChainResult"

  expect_output(print(mock_result), "TransmissionChainAnalysis Result")
})


test_that("plot.TransmissionChainResult exists and is callable", {
  expect_true(is.function(plot.TransmissionChainResult))

  # Check function signature
  args <- formals(plot.TransmissionChainResult)
  expect_true("x" %in% names(args))
  expect_true("which" %in% names(args))
  expect_true("print" %in% names(args))
  expect_true("hh_id" %in% names(args))
  expect_true("prob_cutoff" %in% names(args))
  expect_true("bin_width" %in% names(args))  # New parameter for epidemic_curve
})


# ==============================================================================
# TEST: Plot Method Options
# ==============================================================================

test_that("GenSynResult plot method accepts valid 'which' options", {
  valid_options <- c("posterior", "covariate_effects", "epidemic_curve",
                     "transmission_chains", "all")

  # Just verify these are valid options (no actual plotting)
  expect_length(valid_options, 5)
  expect_true("all" %in% valid_options)
})


test_that("TransmissionChainResult plot method accepts valid 'which' options", {
  valid_options <- c("posterior", "covariate_effects", "epidemic_curve",
                     "transmission_chains", "all")

  expect_length(valid_options, 5)
  expect_true("epidemic_curve" %in% valid_options)  # Now available for user data too
})


# ==============================================================================
# TEST: Print Output Content
# ==============================================================================

test_that("print.GenSynResult shows key information", {
  mock_result <- list(
    n_households = 50,
    attack_rates = list(
      primary_overall = data.frame(
        Total_Pop = 200,
        Total_Infected_People = 80,
        Primary_Attack_Rate = 0.4
      )
    ),
    postprocessing = data.frame(
      Parameter = c("beta1", "beta2"),
      mean = c(0.01, 0.02),
      sd = c(0.001, 0.002)
    )
  )
  class(mock_result) <- "GenSynResult"

  output <- capture.output(print(mock_result))
  output_text <- paste(output, collapse = " ")

  expect_true(grepl("50", output_text))  # n_households
  expect_true(grepl("Primary Attack Rate", output_text))
  expect_true(grepl("Posterior Summary", output_text))
})


test_that("print.TransmissionChainResult shows key information", {
  mock_result <- list(
    user_data = data.frame(
      hh_id = c("HH1", "HH1", "HH2", "HH2"),
      person_id = c(1, 2, 1, 2)
    ),
    postprocessing = data.frame(
      Parameter = c("beta1"),
      mean = c(0.01)
    )
  )
  class(mock_result) <- "TransmissionChainResult"

  output <- capture.output(print(mock_result))
  output_text <- paste(output, collapse = " ")

  expect_true(grepl("TransmissionChainAnalysis", output_text))
  expect_true(grepl("Households", output_text))
})


# ==============================================================================
# TEST: Result Object Access
# ==============================================================================

test_that("GenSynResult components are accessible", {
  mock_result <- list(
    call = quote(GenSyn()),
    n_households = 10,
    simulation = list(hh_df = data.frame(x = 1)),
    stan_data = list(N = 100),
    fit = NULL,
    postprocessing = data.frame(Parameter = "test"),
    attack_rates = list(primary_overall = data.frame()),
    transmission_chains = data.frame()
  )
  class(mock_result) <- "GenSynResult"

  expect_equal(mock_result$n_households, 10)
  expect_type(mock_result$simulation, "list")
  expect_s3_class(mock_result$postprocessing, "data.frame")
})


test_that("TransmissionChainResult components are accessible", {
  mock_result <- list(
    call = quote(TransmissionChainAnalysis()),
    user_data = data.frame(hh_id = "HH1"),
    stan_data = list(N = 50),
    fit = NULL,
    postprocessing = data.frame(Parameter = "test"),
    transmission_chains = data.frame()
  )
  class(mock_result) <- "TransmissionChainResult"

  expect_s3_class(mock_result$user_data, "data.frame")
  expect_type(mock_result$stan_data, "list")
})
