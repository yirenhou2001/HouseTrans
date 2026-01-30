test_that("simulate_multiple_households_comm returns expected structure", {
  skip_on_cran()
  set.seed(1)

  sim <- simulate_multiple_households_comm(
    n_households = 3, alpha_comm_by_role = 5e-3,
    max_days = 15, verbose = FALSE
  )

  expect_type(sim, "list")
  expect_true(all(c("hh_df", "households") %in% names(sim)))
  expect_s3_class(sim$hh_df, "data.frame")
  expect_true(length(sim$households) == 3)

  req_cols <- c("hh_id","person_id","role","infection_time",
                "infectious_start","infectious_end","detection_time",
                "infection_resolved")
  expect_true(all(req_cols %in% names(sim$households[[1]])))

  # Attributes
  expect_true(!is.null(attr(sim$households[[1]], "test_days")) ||
                !is.null(attr(sim$households[[1]], "params")))
})
