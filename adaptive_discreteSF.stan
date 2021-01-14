functions {

  real calculate_beta(real t, real[] beta_params) {
    // standardize time so that values are between -1 and 1 (approximately)
    // TODO: standardization is hardcoded!!
    real standardized_t = (t - 68.0) / 40.0;
    //--------------> simplificamos beta y usamos fabs para asegurar
    // que sea positivo
    return fabs(
      beta_params[1] +
      beta_params[2] * standardized_t
      );
  }

  real[, ] sim_sird(real[] y0_, real[] theta, int N, int T) {
    real sird[T, 4];
    real gamma = theta[1];
    real delta = theta[2];
    real beta_params[2];
    real beta;
    real S;
    real I;
    real R;
    real D;

    beta_params[1] = theta[3];
    beta_params[2] = theta[4];
    //beta_params[3] = theta[5]; 
    //beta_params[4] = theta[6];


    sird[1, 1] = y0_[1];
    sird[1, 2] = y0_[2];
    sird[1, 3] = y0_[3];
    sird[1, 4] = y0_[4];

    for (t in 2:T) {
      S = sird[t - 1, 1];
      I = sird[t - 1, 2];
      R = sird[t - 1, 3];
      D = sird[t - 1, 4];

      beta = calculate_beta(t, beta_params);
      sird[t, 1] = S - beta * I * S / N;
      sird[t, 2] = I + beta * I * S / N - gamma * I - delta * I;
      sird[t, 3] = R + gamma * I;
      sird[t, 4] = D + delta * I;
    }

    return sird;
  }
}

data {
  int<lower=1> n_days;
  real t0;
  real ts[n_days];
  int N;
  int cum_deaths[n_days];
}
parameters {
  real<lower=0> gamma;
  real<lower=0> i0;
  real<lower=0> r0;
  real<lower=0, upper=1> p;
  real<lower=0, upper=1> p2;
  real beta_params[2];
  real<lower=0> delta;
  real<lower=0> phi_inv;
  real<lower=0> phi_death_inv;
}
transformed parameters{
  // -----------> Añadimos esta restriccion para que y sea siempre positivo
  real<lower=0> y[n_days, 4];
  real phi = 1. / phi_inv;
  real phi_death = 1. / phi_death_inv;
  {
    real theta[4];
    real y0_[4];
    theta[1] = gamma;
    theta[2] = delta;
    theta[3] = beta_params[1];
    theta[4] = beta_params[2];
    //theta[5] = beta_params[3];
    //theta[6] = beta_params[4];

    y0_[1] = N - i0 - r0 -1;
    y0_[2] = i0;
    y0_[3] = r0;
    y0_[4] = 1;

    y = sim_sird(y0_, theta, N, n_days);
  }
}
model {
  //priors
  for (i in 1:2) {
    beta_params[i] ~ normal(0, 1);
  }
  // ----> No ponemos énfasis en gamma, así que se deja un prior con
  // una sd alta
  gamma ~ normal(0.1, 2);
  // --> Valores para delta obtenidos de https://www.medrxiv.org/content/10.1101/2020.08.06.20169722v1.full.pdf
  //(0.78, 0.89)
  delta ~ normal(0.00835, 0.00055);  
  phi_inv ~ exponential(5);
  phi_death_inv ~ exponential(5);
  // ------------------> Permitimos a stan inferir condiciones iniciales.
  // Usamos priors con gran sd para no dar mucha info a prior.
  i0 ~ normal(100, 500);
  r0 ~ normal(10, 100);
  
  cum_deaths ~ neg_binomial_2(col(to_matrix(y), 4), phi_death);
}
generated quantities {
  real R0[n_days];
  real beta[n_days];
  real pred_deaths[n_days];
  real real_infected[n_days];

  real recovery_time = 1 / gamma;
  for (t in 1:n_days) {
    beta[t] = calculate_beta(t, beta_params);
    R0[t] = beta[t] / (delta + gamma);
    real_infected[t] = y[t, 2];
  }
  pred_deaths = neg_binomial_2_rng(col(to_matrix(y), 4), phi_death);
}
