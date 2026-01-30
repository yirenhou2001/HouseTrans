 // scenario3_combined.stan
data {
  int<lower=1> N;                   
  int<lower=1> T;                   
  int<lower=1> H;                   
  int<lower=1> R;        
  real delta;

  int<lower=1, upper=H> hh_id[N];      
  int<lower=1, upper=R> role_id[N];    

  int<lower=0, upper=1> I[N, T];  
  int<lower=-1> time_since_infection[N, T];
  real V_term[N, T];               

  real<lower=0> alpha_comm_by_role;        

  int<lower=1> hh_size[H];              
  int<lower=1> hh_max_size;              
  int<lower=1> hh_members[H, hh_max_size];  

  real<lower=0> reference_phi;
  real<lower=0> reference_kappa;
  
  matrix[T, R] seasonal_forcing_mat;
  vector[T] g_profile; 
}

parameters {
  vector[R-1] log_phi_by_role_raw;        
  vector[R-1] log_kappa_by_role_raw;      
  
  real<lower=0> beta1; 
  real<lower=0> beta2;
  real<lower=0.5, upper=10> latent_mean;
}

transformed parameters {
  vector<lower=0>[R] phi_by_role;
  vector<lower=0>[R] kappa_by_role;
  
  real latent_shape = 3.0; 
  real latent_rate = latent_shape / latent_mean;
  vector[T] latent_gate;
  
  for (t in 1:T) {
    latent_gate[t] = gamma_cdf(t | latent_shape, latent_rate);
  }

  phi_by_role[1] = reference_phi;
  for (r in 2:R) phi_by_role[r] = reference_phi * exp(log_phi_by_role_raw[r-1]);

  kappa_by_role[1] = reference_kappa;
  for (r in 2:R) kappa_by_role[r] = reference_kappa * exp(log_kappa_by_role_raw[r-1]);
}

model {
  // Priors
  log_phi_by_role_raw ~ normal(0, 1); 
  log_kappa_by_role_raw ~ normal(0, .5);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 0.5);
  latent_mean ~ normal(2.0, 1.0); 
  
  // Likelihood
  for (n in 1:N) {
    for (t in 1:T) {
      if (t == 1 || sum(I[n, 1:(t-1)]) == 0) {  
        
        real lambda = 0.0;
        
        lambda += phi_by_role[role_id[n]] * alpha_comm_by_role * seasonal_forcing_mat[t, role_id[n]];

        int h = hh_id[n]; 
        real scaling_h = pow(1.0 / max(hh_size[h], 1), delta);

        for (m_idx in 1:hh_size[h]) { 
          int m = hh_members[h, m_idx];
          if (m == n) continue;
          
          int tau = time_since_infection[m, t];
          if (tau < 0) continue; 
          
          // SCENARIO 3 LOGIC:
          // We apply the Soft Gate to BOTH terms, because even if you have VL,
          // if you are technically in the "latent" phase, you shouldn't transmit 
          // (though VL usually handles this itself, the gate enforces consistency for the g(t) term).
          
          real prob_infectious = latent_gate[tau + 1];
          
          real term1 = beta1;// * g_profile[tau + 1];
          real term2 = beta2 * V_term[m, t];
          
          real h_mt = scaling_h * kappa_by_role[role_id[m]] * prob_infectious * (term1 + term2);
          lambda += phi_by_role[role_id[n]] * h_mt; 
        }

        lambda = fmin(fmax(lambda, 1e-12), 1e6);
        target += bernoulli_lpmf(I[n, t] | 1 - exp(-lambda));
      }
    }
  }
}

