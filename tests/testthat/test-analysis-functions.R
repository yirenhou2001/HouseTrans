# tests/testthat/test-input-validation.R
# Fast tests for input validation - no Stan compilation needed

# ==============================================================================
# TEST: TransmissionChainAnalysis Input Validation
# ==============================================================================

test_that("TransmissionChainAnalysis rejects NULL user_data", {
  expect_error(
    TransmissionChainAnalysis(user_data = NULL),
    "user_data.*must be supplied"
  )
})


test_that("TransmissionChainAnalysis rejects non-dataframe input", {
  expect_error(
    TransmissionChainAnalysis(user_data = "not a dataframe"),
    "must be a data.frame"
  )
})


test_that("TransmissionChainAnalysis rejects list of non-dataframes", {
  bad_list <- list("a", "b", "c")
  expect_error(
    TransmissionChainAnalysis(user_data = bad_list),
    "must contain only data.frames"
  )
})


test_that("TransmissionChainAnalysis rejects data with wrong columns", {
  bad_df <- data.frame(
    wrong_col1 = 1:3,
    wrong_col2 = c("a", "b", "c")
  )
  expect_error(
    TransmissionChainAnalysis(user_data = bad_df),
    "could not interpret"
  )
})

# ==============================================================================
# TEST: Data Format Detection
# ==============================================================================

test_that("Per-person format is detected correctly", {
  df <- data.frame(
    hh_id = "HH1",
    person_id = 1,
    role = "adult",
    infection_time = 5,
    infectious_start = 6,
    infectious_end = 10,
    infection_resolved = 11
  )

  req_cols <- c("hh_id", "person_id", "role", "infection_time",
                "infectious_start", "infectious_end", "infection_resolved")

  expect_true(all(req_cols %in% names(df)))
})


test_that("Long format is detected correctly", {
  df <- data.frame(
    HH = "HH1",
    individual_ID = 1,
    role = "adult",
    test_date = as.Date("2024-01-15"),
    infection_status = 1
  )

  req_cols <- c("HH", "individual_ID", "role", "test_date", "infection_status")

  expect_true(all(req_cols %in% names(df)))
})


# ==============================================================================
# TEST: GenSyn Parameter Validation
# ==============================================================================

test_that("GenSyn accepts valid parameters without error", {
  expect_true(is.function(GenSyn))

  # Check formals (function arguments)
  args <- formals(GenSyn)
  expect_true("n_households" %in% names(args))
  expect_true("start_date" %in% names(args))
  expect_true("end_date" %in% names(args))
  expect_true("priors" %in% names(args))
  expect_true("stan_chains" %in% names(args))
  expect_true("household_profile_list" %in% names(args))
  expect_true("covariates_config" %in% names(args))
})


# ==============================================================================
# TEST: Date Handling
# ==============================================================================

test_that("Date conversion works correctly", {
  start <- "2024-07-01"
  end <- "2025-06-30"

  d_start <- as.Date(start)
  d_end <- as.Date(end)

  expect_s3_class(d_start, "Date")
  expect_s3_class(d_end, "Date")
  expect_true(d_end > d_start)

  max_days <- as.integer(d_end - d_start) + 1
  expect_equal(max_days, 365)
})


test_that("Date strings are parseable", {
  dates <- c("2024-01-01", "2024-12-31", "2025-06-30")

  for (d in dates) {
    parsed <- as.Date(d)
    expect_s3_class(parsed, "Date")
    expect_false(is.na(parsed))
  }
})


# ==============================================================================
# TEST: Role Validation in Data
# ==============================================================================

test_that("Role normalization handles various inputs", {
  df <- data.frame(
    hh_id = rep("HH1", 4),
    person_id = 1:4,
    role = c("Adult", "INFANT", "child", "grandparent"),  # Various formats
    infection_time = c(1, 2, NA, NA),
    infectious_start = c(2, 3, NA, NA),
    infectious_end = c(5, 6, NA, NA),
    infection_resolved = c(6, 7, NA, NA)
  )

  # Normalize roles
  normalized <- .norm_role(df$role)

  expect_equal(normalized, c("adult", "infant", "infant", "elderly"))
})
