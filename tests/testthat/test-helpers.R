# tests/testthat/test-helpers.R
# Fast unit tests for helper functions - no Stan compilation needed

# ==============================================================================
# TEST: Role Normalization
# ==============================================================================

test_that(".norm_role normalizes role labels correctly", {
  # Test basic normalization
  expect_equal(.norm_role("Adult"), "adult")
  expect_equal(.norm_role("INFANT"), "infant")
  expect_equal(.norm_role("  toddler  "), "toddler")

  # Test synonyms
  expect_equal(.norm_role("child"), "infant")
  expect_equal(.norm_role("baby"), "infant")
  expect_equal(.norm_role("sibling"), "toddler")
  expect_equal(.norm_role("kid"), "toddler")
  expect_equal(.norm_role("parent"), "adult")
  expect_equal(.norm_role("mother"), "adult")
  expect_equal(.norm_role("father"), "adult")
  expect_equal(.norm_role("elder"), "elderly")
  expect_equal(.norm_role("grandparent"), "elderly")

  # Test vector input
  roles <- c("Adult", "child", "ELDERLY", "sibling")
  expected <- c("adult", "infant", "elderly", "toddler")
  expect_equal(.norm_role(roles), expected)
})


# ==============================================================================
# TEST: Null Coalescing Operator
# ==============================================================================

test_that("%||% returns correct values", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
  expect_equal("a" %||% "b", "a")
  expect_equal(NULL %||% NULL, NULL)
})


# ==============================================================================
# TEST: Viral Load Trajectory
# ==============================================================================

test_that("simulate_viral_load_trajectory returns valid output", {
  t <- 0:10
  result <- simulate_viral_load_trajectory(t, v_p = 5, t_p = 4, lambda_g = 2, lambda_d = 1)

  expect_length(result, length(t))
  expect_true(all(is.finite(result)))
  expect_true(all(result >= -9))
  # Peak should be around t_p
  expect_true(which.max(result) >= 3 && which.max(result) <= 6)
})


test_that("simulate_Ct_trajectory returns valid output", {
  t <- 0:10
  result <- simulate_Ct_trajectory(t, Cpeak = 25, r = 2, d = 1.5, t_peak = 4)

  expect_length(result, length(t))
  expect_true(all(is.finite(result)))
  # Ct values should be positive
  expect_true(all(result > 0))
  # Minimum (peak) should be around t_peak
  expect_true(which.min(result) >= 4 && which.min(result) <= 6)
})


# ==============================================================================
# TEST: Parameter Drawing
# ==============================================================================

test_that("draw_random_VL_params returns correct structure", {
  params <- draw_random_VL_params("adult")

  expect_type(params, "list")
  expect_named(params, c("v_p", "t_p", "lambda_g", "lambda_d"))
  expect_true(all(sapply(params, is.numeric)))
})


test_that("draw_random_VL_params works for all roles", {
  for (role in c("adult", "infant", "toddler", "elderly")) {
    params <- draw_random_VL_params(role)
    expect_type(params, "list")
    expect_true(all(c("v_p", "t_p", "lambda_g", "lambda_d") %in% names(params)))
  }
})


test_that("draw_random_Ct_params returns correct structure", {
  params <- draw_random_Ct_params("infant")

  expect_type(params, "list")
  expect_named(params, c("Cpeak", "r", "d", "t_peak"))
  expect_true(all(sapply(params, is.numeric)))
})


test_that("draw_random_Ct_params works for all roles", {
  for (role in c("adult", "infant", "toddler", "elderly")) {
    params <- draw_random_Ct_params(role)
    expect_type(params, "list")
    expect_true(all(c("Cpeak", "r", "d", "t_peak") %in% names(params)))
  }
})


test_that("draw_random_VL_params fails for unknown role", {
  expect_error(draw_random_VL_params("unknown_role"))
})


test_that("draw_random_Ct_params fails for unknown role", {
  expect_error(draw_random_Ct_params("unknown_role"))
})


# ==============================================================================
# TEST: Household Generation (Updated for new API)
# ==============================================================================

test_that("generate_household_roles returns valid roles", {
  set.seed(123)
  roles <- generate_household_roles()

  expect_type(roles, "character")
  expect_true(length(roles) >= 1)  # At least someone in household
  expect_true(all(roles %in% c("adult", "infant", "toddler", "elderly")))
})


test_that("generate_household_roles respects custom profile", {
  set.seed(456)
  # Force: 1 adult, 1 infant, no siblings, no elderly
  profile <- list(
    prob_adults = c(0, 1, 0),     # Always 1 adult
    prob_infant = 1.0,            # Always 1 infant
    prob_siblings = c(1, 0, 0),   # Always 0 toddlers
    prob_elderly = c(1, 0, 0)     # Always 0 elderly
  )
  roles <- generate_household_roles(profile)

  expect_equal(sum(roles == "adult"), 1)
  expect_equal(sum(roles == "infant"), 1)
  expect_equal(sum(roles == "toddler"), 0)
  expect_equal(sum(roles == "elderly"), 0)
})


test_that("generate_household_roles can have no infant", {
  set.seed(789)
  profile <- list(
    prob_adults = c(0, 0, 1),    # Always 2 adults
    prob_infant = 0.0,           # Never an infant
    prob_siblings = c(1, 0, 0),  # No toddlers
    prob_elderly = c(1, 0, 0)    # No elderly
  )
  roles <- generate_household_roles(profile)

  expect_equal(sum(roles == "adult"), 2)
  expect_equal(sum(roles == "infant"), 0)
})


test_that("generate_household_roles can have 0 adults", {
  set.seed(101)
  profile <- list(
    prob_adults = c(1, 0, 0),    # Always 0 adults
    prob_infant = 1.0,           # Always 1 infant
    prob_siblings = c(1, 0, 0),  # No toddlers
    prob_elderly = c(0, 0, 1)    # Always 2 elderly (grandparents raising child)
  )
  roles <- generate_household_roles(profile)

  expect_equal(sum(roles == "adult"), 0)
  expect_equal(sum(roles == "infant"), 1)
  expect_equal(sum(roles == "elderly"), 2)
})


test_that("generate_household_roles can have 2 adults and 2 toddlers", {
  set.seed(202)
  profile <- list(
    prob_adults = c(0, 0, 1),    # Always 2 adults
    prob_infant = 1.0,           # Always 1 infant
    prob_siblings = c(0, 0, 1),  # Always 2 toddlers
    prob_elderly = c(1, 0, 0)    # No elderly
  )
  roles <- generate_household_roles(profile)

  expect_equal(sum(roles == "adult"), 2)
  expect_equal(sum(roles == "infant"), 1)
  expect_equal(sum(roles == "toddler"), 2)
  expect_equal(sum(roles == "elderly"), 0)
  expect_equal(length(roles), 5)
})


test_that("generate_household_roles uses defaults when NULL passed", {
  set.seed(303)
  roles <- generate_household_roles(NULL)

  expect_type(roles, "character")
  expect_true(length(roles) >= 1)
})


# ==============================================================================
# TEST: Prior Parsing
# ==============================================================================

test_that(".parse_prior handles NULL input", {
  result <- .parse_prior(NULL, 1, c(0, 1))

  expect_equal(result$type, 1)
  expect_equal(result$params, c(0, 1))
})


test_that(".parse_prior parses normal distribution", {
  p <- list(dist = "normal", params = c(-5, 2))
  result <- .parse_prior(p, 2, c(0, 1))

  expect_equal(result$type, 1L)  # Normal = 1
  expect_equal(result$params, c(-5, 2))
})


test_that(".parse_prior parses uniform distribution", {
  p <- list(dist = "uniform", params = c(0, 10))
  result <- .parse_prior(p, 1, c(0, 1))

  expect_equal(result$type, 2L)  # Uniform = 2
  expect_equal(result$params, c(0, 10))
})


test_that(".parse_prior parses lognormal distribution", {
  p <- list(dist = "lognormal", params = c(1.5, 0.5))
  result <- .parse_prior(p, 1, c(0, 1))

  expect_equal(result$type, 3L)  # LogNormal = 3
  expect_equal(result$params, c(1.5, 0.5))
})


test_that(".parse_prior is case insensitive", {
  p1 <- list(dist = "NORMAL", params = c(0, 1))
  p2 <- list(dist = "Normal", params = c(0, 1))
  p3 <- list(dist = "normal", params = c(0, 1))

  expect_equal(.parse_prior(p1, 1, c(0, 1))$type, 1L)
  expect_equal(.parse_prior(p2, 1, c(0, 1))$type, 1L)
  expect_equal(.parse_prior(p3, 1, c(0, 1))$type, 1L)
})


# ==============================================================================
# TEST: Default Parameters Exist
# ==============================================================================

test_that("default parameters are defined", {
  expect_true(exists("default_VL_params"))
  expect_true(exists("default_Ct_params"))
  expect_true(exists("default_priors"))
  expect_true(exists("default_household_profile"))

  # Check structure of VL params
  expect_true(all(c("adult", "infant", "toddler", "elderly") %in% names(default_VL_params)))

  # Check structure of Ct params
  expect_true(all(c("adult", "infant", "toddler", "elderly") %in% names(default_Ct_params)))

  # Check structure of household profile (new API)
  expect_true(all(c("prob_adults", "prob_infant", "prob_siblings", "prob_elderly") %in%
                    names(default_household_profile)))
})


test_that("default_household_profile has valid probabilities", {
  # prob_adults should sum to 1
  expect_equal(sum(default_household_profile$prob_adults), 1)

  # prob_siblings should sum to 1
  expect_equal(sum(default_household_profile$prob_siblings), 1)

  # prob_elderly should sum to 1
  expect_equal(sum(default_household_profile$prob_elderly), 1)

  # prob_infant should be between 0 and 1
  expect_true(default_household_profile$prob_infant >= 0 &&
                default_household_profile$prob_infant <= 1)
})


test_that("default_priors has expected keys", {
  expected_keys <- c("beta1", "beta2", "alpha", "covariates",
                     "gen_shape", "gen_rate", "ct50", "slope")
  expect_true(all(expected_keys %in% names(default_priors)))
})


# ==============================================================================
# TEST: Curve Calculation Helper
# ==============================================================================

test_that(".calc_curve_val works for Ct data", {
  params <- list(Cpeak = 25, r = 2, d = 1.5, t_peak = 5)

  # Before peak
  val_before <- .calc_curve_val(t = 2, p = params, is_ct = TRUE)
  expect_true(val_before > params$Cpeak)  # Ct higher before peak

  # At peak
  val_at_peak <- .calc_curve_val(t = 5, p = params, is_ct = TRUE)
  expect_equal(val_at_peak, params$Cpeak)

  # After peak
  val_after <- .calc_curve_val(t = 8, p = params, is_ct = TRUE)
  expect_true(val_after > params$Cpeak)  # Ct higher after peak
})


test_that(".calc_curve_val works for VL data", {
  params <- list(v_p = 5, t_p = 4, lambda_g = 2, lambda_d = 1)

  # Calculate values at different times
  val_early <- .calc_curve_val(t = 1, p = params, is_ct = FALSE)
  val_peak <- .calc_curve_val(t = 4, p = params, is_ct = FALSE)
  val_late <- .calc_curve_val(t = 10, p = params, is_ct = FALSE)

  # Peak should be highest
  expect_true(val_peak > val_early)
  expect_true(val_peak > val_late)
})
