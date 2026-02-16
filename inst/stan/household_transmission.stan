data {
  int<lower=1> N;                          // Total number of person-episodes

int<lower=1> T;                          // Maximum time horizon (days)
  int<lower=1> H;                          // Number of households
  int<lower=1> R;                          // Number of roles (4: adult, infant, toddler, elderly)
  real delta;                              // Household size scaling exponent

  // --- FLAGS ---
  int<lower=0, upper=1> use_vl_data;       // 1 = use viral load data, 0 = use generation curve only
  int<lower=0, upper=1> vl_type;           // 1 = log10 viral load, 0 = Ct values
  int<lower=0, upper=1> use_curve_logic;   // 1 = use estimated generation curve

  // --- INDIVIDUAL DATA ---
  array[N] int<lower=1, upper=H> hh_id;    // Household ID for each person-episode
  array[N] int<lower=1, upper=R> role_id;  // Role ID for each person-episode

  array[N, T] int<lower=0, upper=1> I;     // Infection indicator (1 on day of infection)
  array[N, T] int<lower=0, upper=1> Y;     // Infectious indicator (1 during infectious period)
  array[N, T] real V;                      // Viral load values (log10 or Ct)
  array[N] int<lower=1, upper=T> start_risk; // Day person becomes at risk
  array[N] int<lower=1> p_id;              // Person ID (same across episodes for reinfections)

  // --- HOUSEHOLD STRUCTURE ---
  array[H] int<lower=1> hh_size_people;    // Number of unique people per household
  int<lower=1> hh_max_size;                // Maximum household size
  array[H, hh_max_size] int<lower=0, upper=N> hh_members; // Member indices per household

  // --- SEASONALITY ---
  matrix[T, R] seasonal_forcing_mat;       // Seasonal forcing by day and role

  // --- REFERENCE VALUES ---
  real<lower=0> reference_phi;             // Reference susceptibility (usually 1.0)
  real<lower=0> reference_kappa;           // Reference infectivity (usually 1.0)

  // --- COVARIATES ---
  int<lower=0> K_susc;                     // Number of susceptibility covariates
  matrix[N, K_susc] X_susc;                // Susceptibility covariate matrix

  int<lower=0> K_inf;                      // Number of infectivity covariates
  matrix[N, K_inf] X_inf;                  // Infectivity covariate matrix

  // --- FLEXIBLE PRIORS ---
  // Prior types: 1=Normal, 2=Uniform, 3=LogNormal

  // Beta 1 (Time/baseline weight)
  int<lower=0> prior_beta1_type;
  vector[2] prior_beta1_params;

  // Beta 2 (Viral load weight)
  int<lower=0> prior_beta2_type;
  vector[2] prior_beta2_params;

  // Alpha (Community transmission rate)
  int<lower=0> prior_alpha_type;
  vector[2] prior_alpha_params;

  // Covariate priors
  int<lower=0> prior_cov_type;
  vector[2] prior_cov_params;

  // Viral dynamics priors
  int<lower=0> prior_shape_type;
  vector[2] prior_shape_params;

  int<lower=0> prior_rate_type;
  vector[2] prior_rate_params;

  int<lower=0> prior_ct50_type;
  vector[2] prior_ct50_params;

  int<lower=0> prior_slope_type;
  vector[2] prior_slope_params;
}

transformed data {
  array[N] int infection_day;

  for (n in 1:N) {
    infection_day[n] = 0;
    for (t in start_risk[n]:T) {
      if (I[n, t] == 1) {
        infection_day[n] = t;
        break;
      }
    }
  }
}

parameters {
  vector[R-1] log_phi_by_role_raw;         // Log susceptibility multipliers (relative to role 1)
  vector[R-1] log_kappa_by_role_raw;       // Log infectivity multipliers (relative to role 1)

  real log_beta1;                          // Log baseline transmission coefficient
  real log_beta2;                          // Log viral load transmission coefficient
  real log_alpha_comm;                     // Log community transmission rate

  real<lower=1.0, upper=20.0> gen_shape;   // Generation interval shape parameter
  real<lower=0.1, upper=5.0> gen_rate;     // Generation interval rate parameter
  real<lower=0> Ct50;                      // Ct value at 50% infectivity
  real<lower=0> slope_ct;                  // Slope of Ct-infectivity relationship

  vector[K_susc] beta_susc;                // Susceptibility covariate coefficients
  vector[K_inf] beta_inf;                  // Infectivity covariate coefficients
}

transformed parameters {
  vector<lower=0>[R] phi_by_role;          // Susceptibility by role
  vector<lower=0>[R] kappa_by_role;        // Infectivity by role
  vector[T] g_curve_est;                   // Estimated generation curve
  real<lower=0> alpha_comm = exp(log_alpha_comm);

  real beta1 = exp(log_beta1);
  real beta2 = exp(log_beta2);

  // Build role-specific parameters (role 1 is reference)
  phi_by_role[1] = reference_phi;
  for (r in 2:R) {
    phi_by_role[r] = reference_phi * exp(log_phi_by_role_raw[r-1]);
  }

  kappa_by_role[1] = reference_kappa;
  for (r in 2:R) {
    kappa_by_role[r] = reference_kappa * exp(log_kappa_by_role_raw[r-1]);
  }

  // Build generation curve from Gamma distribution
  {
    vector[T] raw_curve;
    for (d in 1:T) {
      raw_curve[d] = exp(gamma_lpdf(d | gen_shape, gen_rate));
    }
    g_curve_est = raw_curve / max(raw_curve);
  }

  // Calculate viral load contribution term
  matrix[N, T] V_term_calc = rep_matrix(0.0, N, T);
  if (use_vl_data == 1) {
    for (n in 1:N) {
      for (t in 1:T) {
        if (Y[n, t] == 1) {
          real val = V[n, t];
          if (vl_type == 1) {
            // Log10 viral load: power function
            V_term_calc[n, t] = pow(fmax(0.0, val) / Ct50, slope_ct);
          } else {
            // Ct values: logistic function (lower Ct = higher VL = more infectious)
            V_term_calc[n, t] = inv_logit((Ct50 - val) / slope_ct);
          }
        }
      }
    }
  }
}

model {
  // --- ROLE PRIORS (Fixed Normal) ---
  log_phi_by_role_raw ~ normal(0, 1);
  log_kappa_by_role_raw ~ normal(0, 1);

  // --- TRANSMISSION PRIORS (Flexible) ---
  if (prior_beta1_type == 1) {
    log_beta1 ~ normal(prior_beta1_params[1], prior_beta1_params[2]);
  } else if (prior_beta1_type == 2) {
    log_beta1 ~ uniform(prior_beta1_params[1], prior_beta1_params[2]);
  }

  if (prior_beta2_type == 1) {
    log_beta2 ~ normal(prior_beta2_params[1], prior_beta2_params[2]);
  } else if (prior_beta2_type == 2) {
    log_beta2 ~ uniform(prior_beta2_params[1], prior_beta2_params[2]);
  }

  if (prior_alpha_type == 1) {
    log_alpha_comm ~ normal(prior_alpha_params[1], prior_alpha_params[2]);
  } else if (prior_alpha_type == 2) {
    log_alpha_comm ~ uniform(prior_alpha_params[1], prior_alpha_params[2]);
  }

  // --- COVARIATE PRIORS ---
  if (K_susc > 0) {
    if (prior_cov_type == 1) {
      beta_susc ~ normal(prior_cov_params[1], prior_cov_params[2]);
    } else if (prior_cov_type == 2) {
      beta_susc ~ uniform(prior_cov_params[1], prior_cov_params[2]);
    }
  }

  if (K_inf > 0) {
    if (prior_cov_type == 1) {
      beta_inf ~ normal(prior_cov_params[1], prior_cov_params[2]);
    } else if (prior_cov_type == 2) {
      beta_inf ~ uniform(prior_cov_params[1], prior_cov_params[2]);
    }
  }

  // --- VIRAL DYNAMICS PRIORS (Flexible) ---
  if (prior_shape_type == 1) {
    gen_shape ~ normal(prior_shape_params[1], prior_shape_params[2]);
  } else if (prior_shape_type == 2) {
    gen_shape ~ uniform(prior_shape_params[1], prior_shape_params[2]);
  } else if (prior_shape_type == 3) {
    gen_shape ~ lognormal(prior_shape_params[1], prior_shape_params[2]);
  }

  if (prior_rate_type == 1) {
    gen_rate ~ normal(prior_rate_params[1], prior_rate_params[2]);
  } else if (prior_rate_type == 2) {
    gen_rate ~ uniform(prior_rate_params[1], prior_rate_params[2]);
  } else if (prior_rate_type == 3) {
    gen_rate ~ lognormal(prior_rate_params[1], prior_rate_params[2]);
  }

  if (prior_ct50_type == 1) {
    Ct50 ~ normal(prior_ct50_params[1], prior_ct50_params[2]);
  } else if (prior_ct50_type == 2) {
    Ct50 ~ uniform(prior_ct50_params[1], prior_ct50_params[2]);
  } else if (prior_ct50_type == 3) {
    Ct50 ~ lognormal(prior_ct50_params[1], prior_ct50_params[2]);
  }

  if (prior_slope_type == 1) {
    slope_ct ~ normal(prior_slope_params[1], prior_slope_params[2]);
  } else if (prior_slope_type == 2) {
    slope_ct ~ uniform(prior_slope_params[1], prior_slope_params[2]);
  } else if (prior_slope_type == 3) {
    slope_ct ~ lognormal(prior_slope_params[1], prior_slope_params[2]);
  }

  // --- LIKELIHOOD ---
  for (n in 1:N) {
    int t_stop = (infection_day[n] == 0) ? T : infection_day[n];

    // Calculate effective susceptibility with covariates
    real log_susc_mod = 0.0;
    if (K_susc > 0) {
      log_susc_mod = dot_product(X_susc[n], beta_susc);
    }
    real phi_eff = phi_by_role[role_id[n]] * exp(log_susc_mod);

    for (t in start_risk[n]:t_stop) {
      real lambda = 0.0;

      // Community transmission
      lambda += phi_eff * alpha_comm * seasonal_forcing_mat[t, role_id[n]];

      // Household transmission
      int h = hh_id[n];
      real scaling_h = pow(1.0 / max(hh_size_people[h], 1), delta);

      for (m_idx in 1:hh_max_size) {
        int m = hh_members[h, m_idx];
        if (m == 0) break;
        if (p_id[m] == p_id[n]) continue;  // Same person (different episode)
        if (m == n) continue;               // Same person-episode
        if (Y[m, t] == 0) continue;         // Not infectious

        // Calculate effective infectivity with covariates
        real log_inf_mod = 0.0;
        if (K_inf > 0) {
          log_inf_mod = dot_product(X_inf[m], beta_inf);
        }
        real kappa_eff = kappa_by_role[role_id[m]] * exp(log_inf_mod);

        // Generation curve component
        real g_t = 0.0;
        int m_infection_day = infection_day[m];
        if (m_infection_day != 0 && t >= m_infection_day) {
          int dt = t - m_infection_day + 1;
          if (use_curve_logic == 1) {
            if (dt <= T) g_t = g_curve_est[dt];
          } else {
            g_t = 1.0;
          }
        }

        // Viral load component
        real v_comp = 0.0;
        if (use_vl_data == 1) {
          v_comp = V_term_calc[m, t];
        }

        // Combined transmission term
        real term_combined;
        if (use_vl_data == 0) {
          term_combined = beta1 + (beta2 * g_t);
        } else {
          term_combined = beta1 + (beta2 * v_comp);
        }

        real h_mt = scaling_h * kappa_eff * term_combined;
        lambda += phi_eff * h_mt;
      }

      // Clamp lambda to avoid numerical issues
      lambda = fmin(fmax(lambda, 1e-12), 1e6);

      // Bernoulli likelihood
      int outcome = (t == infection_day[n]) ? 1 : 0;
      target += bernoulli_lpmf(outcome | 1 - exp(-lambda));
    }
  }
}
