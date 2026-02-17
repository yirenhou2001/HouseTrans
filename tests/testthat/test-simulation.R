# tests/testthat/test-simulation.R
# Fast tests for simulation-related functions - no Stan compilation needed

# ==============================================================================
# TEST: Household Generation Integration
# ==============================================================================

test_that("generate_household_roles produces reproducible results", {
  set.seed(42)
  roles1 <- generate_household_roles()

  set.seed(42)
  roles2 <- generate_household_roles()

  expect_equal(roles1, roles2)
})


test_that("generate_household_roles default produces expected composition", {
  # Run multiple times to check distribution
  set.seed(999)
  role_counts <- list(adult = 0, infant = 0, toddler = 0, elderly = 0)

  for (i in 1:100) {
    roles <- generate_household_roles()
    for (r in roles) {
      role_counts[[r]] <- role_counts[[r]] + 1
    }
  }

  # With defaults, should always have 2 adults (prob_adults = c(0,0,1))
  # and always 1 infant (prob_infant = 1.0)
  # So expect ~200 adults and ~100 infants over 100 households
  expect_true(role_counts$adult >= 150)  # Should be close to 200
  expect_true(role_counts$infant >= 80)  # Should be close to 100
})


# ==============================================================================
# TEST: Viral Load Trajectory Shapes
# ==============================================================================

test_that("VL trajectory has correct shape (rise then fall)", {
  t <- 0:20
  vl <- simulate_viral_load_trajectory(t, v_p = 5, t_p = 5, lambda_g = 2, lambda_d = 1)

  # Find peak
  peak_idx <- which.max(vl)

  # Values before peak should be increasing
  if (peak_idx > 2) {
    expect_true(vl[peak_idx] > vl[1])
  }

  # Values after peak should be decreasing
  if (peak_idx < length(vl) - 1) {
    expect_true(vl[peak_idx] > vl[length(vl)])
  }
})


test_that("Ct trajectory has correct shape (fall then rise)", {
  t <- 0:20
  ct <- simulate_Ct_trajectory(t, Cpeak = 20, r = 3, d = 2, t_peak = 5)

  # Find minimum (peak infection = lowest Ct)
  min_idx <- which.min(ct)

  # Values before minimum should be decreasing (higher Ct)
  if (min_idx > 2) {
    expect_true(ct[1] > ct[min_idx])
  }

  # Values after minimum should be increasing (lower viral load = higher Ct)
  if (min_idx < length(ct) - 1) {
    expect_true(ct[length(ct)] > ct[min_idx])
  }
})


# ==============================================================================
# TEST: Role-Specific Parameter Differences
# ==============================================================================

test_that("Different roles have different VL parameters", {
  adult_params <- draw_random_VL_params("adult")
  infant_params <- draw_random_VL_params("infant")

  # Infants typically have higher peak viral load
  expect_true(infant_params$v_p > adult_params$v_p)
})


test_that("Different roles have different Ct parameters", {
  adult_params <- draw_random_Ct_params("adult")
  infant_params <- draw_random_Ct_params("infant")

  # Both should have valid parameters
  expect_true(adult_params$Cpeak > 0)
  expect_true(infant_params$Cpeak > 0)
})


# ==============================================================================
# TEST: Edge Cases in Household Generation
# ==============================================================================

test_that("Empty household is possible with extreme profile", {
  profile <- list(
    prob_adults = c(1, 0, 0),    # 0 adults
    prob_infant = 0.0,           # 0 infants
    prob_siblings = c(1, 0, 0),  # 0 toddlers
    prob_elderly = c(1, 0, 0)    # 0 elderly
  )

  set.seed(123)
  roles <- generate_household_roles(profile)

  expect_equal(length(roles), 0)
})


test_that("Large household is possible", {
  profile <- list(
    prob_adults = c(0, 0, 1),    # 2 adults
    prob_infant = 1.0,           # 1 infant
    prob_siblings = c(0, 0, 1),  # 2 toddlers
    prob_elderly = c(0, 0, 1)    # 2 elderly
  )

  set.seed(123)
  roles <- generate_household_roles(profile)

  expect_equal(length(roles), 7)  # 2 + 1 + 2 + 2
})


# ==============================================================================
# TEST: Priors Configuration
# ==============================================================================

test_that("default_priors can be modified", {
  custom_priors <- default_priors
  custom_priors$beta1 <- list(dist = "uniform", params = c(-10, 0))

  expect_equal(custom_priors$beta1$dist, "uniform")
  expect_equal(custom_priors$beta1$params, c(-10, 0))

  # Original should be unchanged
  expect_equal(default_priors$beta1$dist, "normal")
})


test_that("all default priors have valid structure", {
  for (prior_name in names(default_priors)) {
    prior <- default_priors[[prior_name]]

    expect_true("dist" %in% names(prior),
                info = paste("Prior", prior_name, "missing 'dist'"))
    expect_true("params" %in% names(prior),
                info = paste("Prior", prior_name, "missing 'params'"))
    expect_true(prior$dist %in% c("normal", "uniform", "lognormal"),
                info = paste("Prior", prior_name, "has invalid dist"))
  }
})


# ==============================================================================
# TEST: Trajectory Parameter Bounds
# ==============================================================================

test_that("VL parameters are within reasonable bounds", {
  for (role in c("adult", "infant", "toddler", "elderly")) {
    params <- draw_random_VL_params(role)

    expect_true(params$v_p > 0 && params$v_p < 10,
                info = paste(role, "v_p out of bounds"))
    expect_true(params$t_p > 0 && params$t_p < 20,
                info = paste(role, "t_p out of bounds"))
    expect_true(params$lambda_g > 0,
                info = paste(role, "lambda_g should be positive"))
    expect_true(params$lambda_d > 0,
                info = paste(role, "lambda_d should be positive"))
  }
})


test_that("Ct parameters are within reasonable bounds", {
  for (role in c("adult", "infant", "toddler", "elderly")) {
    params <- draw_random_Ct_params(role)

    expect_true(params$Cpeak > 15 && params$Cpeak < 45,
                info = paste(role, "Cpeak out of bounds"))
    expect_true(params$t_peak > 0 && params$t_peak < 15,
                info = paste(role, "t_peak out of bounds"))
    expect_true(params$r > 0,
                info = paste(role, "r should be positive"))
    expect_true(params$d > 0,
                info = paste(role, "d should be positive"))
  }
})
