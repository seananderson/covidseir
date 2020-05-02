functions{
  real[] sir(real t,        // time (actual time; not an increment starting at 1)
             real[] state,  // state
             real[] theta,  // parameters
             real[] x_r,    // data (real)
             int[]  x_i) {  // data (integer)
    real S     = state[1];
    real E1    = state[2];
    real E2    = state[3];
    real I     = state[4];
    real Q     = state[5];
    real R     = state[6];
    real Sd    = state[7];
    real E1d   = state[8];
    real E2d   = state[9];
    real Id    = state[10];
    real Qd    = state[11];
    real Rd    = state[12];

    real N     = x_r[1];
    real D     = x_r[2];
    real k1    = x_r[3];
    real k2    = x_r[4];
    real q     = x_r[5];
    real r     = x_r[6];
    real ur    = x_r[7];
    real f1    = x_r[8];

    real start_decline = x_r[9];
    real end_decline = x_r[10];
    real fixed_f_forecast = x_r[11];
    real last_day_obs = x_r[12];
    real day_start_fixed_f_forecast = x_r[13];

    real R0 = theta[1];

    int n_f = x_i[1]; // the number of f parameters in time (after f1)
    int f_seg_id[n_f]; // a lookup vector to grab the appropriate f parameter in time

    real f; // will store the f value for this time point

    real dydt[12];

    // integer version of the day for this time point:
    // (must use this workaround since floor() in Stan produces a real number
    // that can't be used to index an array)
    int day;
    day = 1;
    while ((day + 1) < floor(t)) day = day + 1;

    for (i in 1:n_f) {
      // `i + 1` because of number of x_i before `f_seg_id`
      // `+ 1` at end because of number of thetas (here just R0) before f thetas
      f_seg_id[i] = x_i[i + 1] + 1;
    }

    f = f1; // business as usual before physical distancing
    if (t < start_decline) {
      f = f1;
    }
    if (t >= start_decline && t < end_decline) { // allow a ramp-in of physical distancing
      f = theta[f_seg_id[day]] + (end_decline - t) *
          (f1 - theta[f_seg_id[day]]) / (end_decline - start_decline);
    }
    if (t >= end_decline) {
      f = theta[f_seg_id[day]];
    }
    if (t >= day_start_fixed_f_forecast && fixed_f_forecast != 0) {
      f = fixed_f_forecast;
    }

    dydt[1]  = -(R0/(D+1/k2)) * (I + E2 + f*(Id+E2d)) * S/N - r*S + ur*Sd;
    dydt[2]  = (R0/(D+1/k2)) * (I + E2 + f*(Id+E2d)) * S/N - k1*E1 -r*E1 + ur*E1d;
    dydt[3]  = k1*E1 - k2*E2 - r*E2 + ur*E2d;
    dydt[4]  = k2*E2 - q*I - I/D - r*I + ur*Id;
    dydt[5]  = q*I - Q/D - r*Q + ur*Qd;
    dydt[6]  = I/D + Q/D - r*R + ur*Rd;

    dydt[7]  = -(f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N + r*S - ur*Sd;
    dydt[8]  = (f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N - k1*E1d +r*E1 - ur*E1d;
    dydt[9]  = k1*E1d - k2*E2d + r*E2 - ur*E2d;
    dydt[10] = k2*E2d - q*Id - Id/D + r*I - ur*Id;
    dydt[11] = q*Id - Qd/D + r*Q - ur*Qd;
    dydt[12] = Id/D + Qd/D + r*R - ur*Rd;

    return dydt;
  }
}
data {
  int<lower=1> T;     // number of time steps
  int<lower=1> N;     // number of days
  int<lower=1> J;     // number of response data timeseries
  int<lower=1> S;     // number of physical distancing segments
  real y0[12];        // initial state
  real t0;            // first time step
  real time[T];       // time increments
  int days[N];        // day increments
  int last_day_obs;   // last day of observed data; days after this are projections
  int daily_cases[last_day_obs,J]; // daily new case counts
  int sampFrac_seg[N]; // optional index of estimated sample fractions for 1st timeseries
  int<lower=1, upper=4> sampFrac_type; // Type of sample fraction: fixed, estimated, rw, segmented
  real x_r[13];       // data for ODEs (real numbers)
  int n_x_i;          // the number of x_i values
  int x_i[n_x_i];    // data for ODEs (integer numbers)
  real sampFrac[N,J];   // fraction of cases sampled per time step
  real delayScale[J];    // Weibull parameter for delay in becoming a case count
  real delayShape[J];    // Weibull parameter for delay in becoming a case count
  int time_day_id[N]; // last time increment associated with each day
  int time_day_id0[N];// first time increment for Weibull integration of case counts
  real R0_prior[2];   // lognormal log mean and SD for R0 prior
  real phi_prior;     // SD of normal prior on 1/sqrt(phi) [NB2(mu, phi)]
  real f2_prior[2];   // beta prior for f2
  real sampFrac2_prior[2];   // beta prior for sampFrac2
  int<lower=0, upper=1> priors_only; // logical: include likelihood or just priors?
  int<lower=0, upper=J> est_phi; // estimate NB phi?
  int<lower=0, upper=N> n_sampFrac2; // number of sampFrac2
  int<lower=0, upper=1> obs_model; // observation model: 0 = Poisson, 1 = NB2
  real<lower=0> rw_sigma; // specified random walk standard deviation
  real ode_control[3]; // vector of ODE control numbers
}
parameters {
 real R0; // Stan ODE solver seems to be more efficient without this bounded at > 0
 real<lower=0, upper=1> f_s[S]; // strength of social distancing for segment `s`
 real<lower=0> phi[est_phi]; // NB2 (inverse) dispersion; `est_phi` turns on/off
 real<lower=0, upper=1> sampFrac2[n_sampFrac2];
}
transformed parameters {
  real dx = time[2] - time[1]; // time increment
  real ft[T]; // container for the lambda function at time t
  real lambda_d[N,J]; // estimated daily cases for each day
  real sum_ft_inner; // holds a temporary calculation
  real eta[N,J]; // expected value on link scale (log)
  real k2; // from ODE
  real E2; // from ODE exposed and symptomatic
  real E2d; // from ODE exposed and symptomatic and distancing
  real theta[2]; // gathers up the parameters (which come in with various limits)
  real y_hat[T,12]; // predicted states for each time t from ODE
  real this_samp; // holds the sample fraction for a given day

  theta[1] = R0;
  for (s in 1:S) {
    // `s + 1` because of number of thetas before f_s (just R0)
    theta[s + 1] = f_s[s];
  }

  y_hat = integrate_ode_rk45(sir, y0, t0, time, theta, x_r, x_i,
                             ode_control[1], ode_control[2], ode_control[3]);

  // Calculating the expected case counts given the delays in reporting:
  for (j in 1:J) { // data_type increment
    for (n in 1:N) { // day increment
      this_samp = sampFrac[n,j];
      if (n_sampFrac2 > 1 && j == 1) {
        if (sampFrac_type != 4) { // anything but segmented sampFrac2
          if (n <= last_day_obs) {
            this_samp = sampFrac2[n];
          }
          if (n > last_day_obs && j == 1) {
            this_samp = sampFrac2[n_sampFrac2]; // forecast with last value
          }
        } else { // segmented
          this_samp = sampFrac2[sampFrac_seg[n]];
        }
      }
      if (n_sampFrac2 == 1 && j == 1) {
        this_samp = sampFrac2[1];
      }
      for (t in 1:T) {
        ft[t] = 0; // initialize across the full 1:T
      }
      // a fancy way of moving across a window of time:
      for (t in time_day_id0[n]:time_day_id[n]) { // t is an increment here
        k2 = x_r[4];
        E2 = y_hat[t,3];
        E2d = y_hat[t,9];

        ft[t] = this_samp * k2 * (E2 + E2d) *
        exp(weibull_lpdf(time[time_day_id[n]] - time[t] | delayShape[j], delayScale[j]));
      }
      sum_ft_inner = 0; // initialize
      for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
        sum_ft_inner += ft[t];
      }
      // trapezoid integration:
      lambda_d[n,j] = 0.5 * dx *
      (ft[time_day_id0[n]] + 2 * sum_ft_inner + ft[time_day_id[n]]);
      eta[n,j] = log(lambda_d[n,j]);
    }
  }
}
model {
  // priors:
  if (est_phi > 0 && obs_model == 1) { // NB2
  // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  // D(expression(1/sqrt(x)), "x"); log(0.5 * x^-0.5/sqrt(x)^2
  for (j in 1:J) {
    1/sqrt(phi[j]) ~ normal(0, phi_prior);
    target += log(0.5) - 1.5 * log(phi[j]); // Jacobian adjustment
  }
  }
  R0 ~ lognormal(R0_prior[1], R0_prior[2]);
  for (s in 1:S) {
    f_s[s] ~ beta(f2_prior[1], f2_prior[2]); // FIXME: allow separate priors
  }
  if (n_sampFrac2 > 0 && sampFrac_type != 4) { // sampFrac estimated but not segmented
    sampFrac2[1] ~ beta(sampFrac2_prior[1], sampFrac2_prior[2]);
    if (n_sampFrac2 > 1) {
      for (n in 2:n_sampFrac2) {
        sampFrac2[n] ~ normal(sampFrac2[n - 1], rw_sigma); // RW
      }
    }
  }
  if (n_sampFrac2 > 0 && sampFrac_type != 4) { // sampFrac segmented
    for (n in 1:n_sampFrac2) {
      sampFrac2[n] ~ beta(sampFrac2_prior[1], sampFrac2_prior[2]);
    }
  }

  // data likelihood:
  if (!priors_only) { // useful to turn off for prior predictive checks
    for (n in 1:last_day_obs) {
      for (j in 1:J) {
        if (daily_cases[n,j] != 9999999) { // NA magic number
        if (obs_model == 0) {
          daily_cases[n,j] ~ poisson_log(eta[n,j]);
        } else if (obs_model == 1) {
          daily_cases[n,j] ~ neg_binomial_2_log(eta[n,j], phi[j]);
        }
        }
      }
    }
  }
}
generated quantities{
  int y_rep[N,J]; // posterior predictive replicates
  for (j in 1:J) {
    for (n in 1:N) {
      if (obs_model == 0) {
        y_rep[n,j] = poisson_log_rng(eta[n,j]);
      } else if (obs_model == 1) {
        y_rep[n,j] = neg_binomial_2_log_rng(eta[n,j], phi[1]);
      }
    }
  }
}

