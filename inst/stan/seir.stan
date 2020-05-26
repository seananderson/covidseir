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
    real ud     = x_r[6];
    // real ur    = x_r[7];
    real f0    = x_r[8];
    real f_ramp_rate = x_r[9];
    real imported_cases = x_r[10];
    real imported_window = x_r[11];

    // real start_decline = x_r[9];
    // real end_decline = x_r[10];

    real last_day_obs = x_i[1];

    real R0 = theta[1];
    real start_decline = theta[2];
    real end_decline = theta[3];
    real ur = theta[4];
    real f1 = theta[5];

    int n_f = x_i[2]; // the number of f parameters in time (after f0)
    int f_seg_id[n_f]; // a lookup vector to grab the appropriate f parameter in time

    real f; // will store the f value for this time point
    real introduced; // will store the introduced cases

    real dydt[12];

    real X;
    // integer version of the day for this time point:
    // (must use this workaround since floor() in Stan produces a real number
    // that can't be used to index an array)
    int day;
    day = 1;
    while ((day + 1) < floor(t)) day = day + 1;
    // int time_int;
    // time_int = 1;
    // while ((time_int + 1) < floor(t)) time_int = time_int + 1;
    // this_time_id = f_seg_id[time_int];

    for (i in 1:n_f) {
      // `i + 2` because of number of x_i before `f_seg_id`
      // `+ 4` at end because of number of thetas before f_s thetas
      f_seg_id[i] = x_i[i + 2] + 4;
    }

    f = f0; // business as usual before physical distancing
    if (t < start_decline) {
      f = f0;
    }
     // allow a ramp-in of physical distancing:
    if (f_ramp_rate == 0.0) {
      if (t >= start_decline && t < end_decline) {
      f = f1 + (end_decline - t) *
      (f0 - f1) / (end_decline - start_decline);
      }
    } else {
      if (t >= start_decline && t < end_decline) {
      X = (f0 - f1) / (exp(f_ramp_rate * (end_decline - start_decline)) - 1);
      f = f0 - X * (exp(f_ramp_rate * (t - start_decline)) - 1);
      }
    }
    if (t >= end_decline) {
      f = theta[f_seg_id[day]]; // the respective f segment
    }
    if (t > last_day_obs && t <= (last_day_obs + imported_window)) {
      introduced = imported_cases / imported_window;
    } else {
      introduced = 0.0;
    }

    dydt[1]  = -(R0/(D+1/k2)) * (I + E2 + f*(Id+E2d)) * S/N - ud*S + ur*Sd;
    dydt[2]  = (R0/(D+1/k2)) * (I + E2 + f*(Id+E2d)) * S/N - k1*E1 -ud*E1 + ur*E1d;
    dydt[3]  = k1*E1 - k2*E2 - ud*E2 + ur*E2d + introduced;
    dydt[4]  = k2*E2 - q*I - I/D - ud*I + ur*Id;
    dydt[5]  = q*I - Q/D - ud*Q + ur*Qd;
    dydt[6]  = I/D + Q/D - ud*R + ur*Rd;

    dydt[7]  = -(f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N + ud*S - ur*Sd;
    dydt[8]  = (f*R0/(D+1/k2)) * (I+E2 + f*(Id+E2d)) * Sd/N - k1*E1d +ud*E1 - ur*E1d;
    dydt[9]  = k1*E1d - k2*E2d + ud*E2 - ur*E2d;
    dydt[10] = k2*E2d - q*Id - Id/D + ud*I - ur*Id;
    dydt[11] = q*Id - Qd/D + ud*Q - ur*Qd;
    dydt[12] = Id/D + Qd/D + ud*R - ur*Rd;

    return dydt;
  }
}
data {
  int<lower=1> T;     // number of time steps
  int<lower=1> N;     // number of days
  int<lower=1> J;     // number of response data timeseries
  int<lower=1> S;     // number of physical distancing segments
  real y0_vars[10];   // initial state
  real t0;            // first time step
  real time[T];       // time increments
  int days[N];        // day increments
  int last_day_obs;   // last day of observed data; days after this are projections
  int daily_cases[last_day_obs,J]; // daily new case counts
  int samp_frac_seg[N]; // optional index of estimated sample fractions for 1st timeseries
  int<lower=1, upper=4> samp_frac_type; // Type of sample fraction: fixed, estimated, rw, segmented
  int n_x_r;         // the number of x_r values
  real x_r[n_x_r];   // data for ODEs (real numbers)
  int n_x_i;         // the number of x_i values
  int x_i[n_x_i];    // data for ODEs (integer numbers)
  real samp_frac_fixed[N,J];   // fraction of cases sampled per time step
  real delay_scale[J];    // Weibull parameter for delay in becoming a case
  real delay_shape[J];    // Weibull parameter for delay in becoming a case
  int time_day_id[N]; // last time increment associated with each day
  int time_day_id0[N];// first time increment for Weibull integration of cases
  real R0_prior[2];   // lognormal log mean and SD for R0 prior
  real i0_prior[2];   // lognormal log mean and SD for i0 prior
  real e_prior[2];   // beta prior on fraction social distancing (e)
  real phi_prior;     // SD of normal prior on 1/sqrt(phi) [NB2(mu, phi)]
  real f_prior[2];   // beta prior for f2
  real samp_frac_prior[2];   // beta prior for samp_frac
  real start_decline_prior[2];   // prior for start_decline day
  real end_decline_prior[2];   // prior for end_decline day
  int<lower=0, upper=1> priors_only; // logical: include likelihood or just priors?
  int<lower=0, upper=J> est_phi; // estimate NB phi?
  int<lower=0, upper=N> n_samp_frac; // number of samp_frac
  int<lower=0, upper=1> obs_model; // observation model: 0 = Poisson, 1 = NB2
  real<lower=0> rw_sigma; // specified random walk standard deviation
  int<lower=0, upper=1> contains_NAs; // Logical: contains NA values?
  real ode_control[3]; // vector of ODE control numbers
}
parameters {
 real<lower=0, upper=x_r[1]> i0; // incidence at initial time point (default -30 days); upper = N
 real R0; // Stan ODE solver seems to be more efficient without this bounded at > 0
 real<lower=0> start_decline;
 real<lower=0> end_decline;
 real<lower=0, upper=1> f_s[S]; // strength of social distancing for segment `s`
 real<lower=0> phi[est_phi]; // NB2 (inverse) dispersion; `est_phi` turns on/off
 real<lower=0, upper=1> samp_frac[n_samp_frac];
 real<lower=0, upper=1> ur; // reasonable for ud ~= 0.1; about 10% FIXME
}
transformed parameters {
  real dx = time[2] - time[1]; // time increment
  real ft[T]; // container for the lambda function at time t
  real mu[N,J]; // estimated daily cases for each day
  real sum_ft_inner; // holds a temporary calculation
  real eta[N,J]; // expected value on link scale (log)
  real k2; // from ODE
  real E2; // from ODE exposed and symptomatic
  real E2d; // from ODE exposed and symptomatic and distancing
  real theta[S + 4]; // gathers parameters (which come with various limits); + 4 is number of thetas before f_s
  real y_hat[T,12]; // predicted states for each time t from ODE
  real this_samp; // holds the sample fraction for a given day
  real y0[12]; // initial states
  real fsi; // fraction social distancing
  real nsi; // fraction not social distancing
  real ud = x_r[6]; // grab this rate input
  //real ur = x_r[7]; // grab this rate input
  real N_pop = x_r[1]; // grab population size

  fsi = ud / (ud + ur); // fraction social distancing
  nsi = 1 - fsi; // fraction not social distancing

  // set up the initial state:
  y0[1] = nsi * (N_pop - i0); // S = nsi * (pars[["N"]] - i0)
  y0[2] = y0_vars[1] * nsi * i0; // E1 = E1_frac * nsi * i0
  y0[3] = y0_vars[2] * nsi * i0; // E2 = E2_frac * nsi * i0
  y0[4] = y0_vars[3] * nsi * i0; // I = I_frac * nsi * i0
  y0[5] = y0_vars[4]; // Q = Q_num
  y0[6] = y0_vars[5]; // R = R_num
  y0[7] = fsi * (N_pop - i0); // Sd = fsi * (pars[["N"]] - i0)
  y0[8] = y0_vars[6] * fsi * i0; // E1d = E1d_frac * fsi * i0
  y0[9] = y0_vars[7] * fsi * i0; // E2d = E2d_frac * fsi * i0
  y0[10] = y0_vars[8] * fsi * i0; // Id = Id_frac * fsi * i0
  y0[11] = y0_vars[9]; // Qd = Qd_num
  y0[12] = y0_vars[10]; // Rd = Rd_num

  theta[1] = R0;
  theta[2] = start_decline;
  theta[3] = end_decline;
  theta[4] = ur;
  for (s in 1:S) {
    // `s + 4` because of number of thetas before f_s
    theta[s + 4] = f_s[s];
  }

  y_hat = integrate_ode_rk45(sir, y0, t0, time, theta, x_r, x_i,
  // y_hat = integrate_ode_bdf(sir, y0, t0, time, theta, x_r, x_i,
                             ode_control[1], ode_control[2], ode_control[3]);

  // Calculating the expected case counts given the delays in reporting:
  for (j in 1:J) { // data_type increment
    for (n in 1:N) { // day increment
      this_samp = samp_frac_fixed[n,j];
      if (n_samp_frac > 1 && j == 1) {
        if (samp_frac_type != 4) { // anything but segmented samp_frac
          if (n <= last_day_obs) {
            this_samp = samp_frac[n];
          }
          if (n > last_day_obs && j == 1) {
            this_samp = samp_frac[n_samp_frac]; // forecast with last value
          }
        } else { // segmented
          this_samp = samp_frac[samp_frac_seg[n]];
        }
      }
      if (n_samp_frac == 1 && j == 1) {
        this_samp = samp_frac[1];
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
        exp(weibull_lpdf(time[time_day_id[n]] - time[t] | delay_shape[j], delay_scale[j]));
      }
      sum_ft_inner = 0; // initialize
      for (t in (time_day_id0[n] + 1):(time_day_id[n] - 1)) {
        sum_ft_inner += ft[t];
      }
      // trapezoid integration:
      mu[n,j] = 0.5 * dx *
      (ft[time_day_id0[n]] + 2 * sum_ft_inner + ft[time_day_id[n]]);
      eta[n,j] = log(mu[n,j]);
    }
  }
}
model {
  // priors:
  if (est_phi > 0 && obs_model == 1) { // NB2
  // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  // log(abs(D(expression(1/sqrt(x)), "x"))); log(0.5 * x^-0.5/sqrt(x)^2
  // log absolute derivative of the transform
  for (j in 1:J) {
    1/sqrt(phi[j]) ~ normal(0, phi_prior);
    target += log(0.5) - 1.5 * log(phi[j]); // Jacobian adjustment
  }
  }
  R0 ~ lognormal(R0_prior[1], R0_prior[2]);
  i0 ~ lognormal(i0_prior[1], i0_prior[2]);
  fsi ~ beta(e_prior[1], e_prior[2]); // two names
  // log(abs(D(expression(ud / (ud + ur)), "ur"))); log(abs(-(ud/(ud + ur)^2))
  // log absolute derivative of the transform
  target += log(ud) - 2 * log(ud + ur); // Jacobian adjustment

  start_decline ~ lognormal(start_decline_prior[1], start_decline_prior[2]);
  end_decline ~ lognormal(end_decline_prior[1], end_decline_prior[2]);
  for (s in 1:S) {
    f_s[s] ~ beta(f_prior[1], f_prior[2]); // FIXME: allow separate priors
  }
  if (n_samp_frac > 0 && samp_frac_type != 4) { // samp_frac estimated but not segmented
    samp_frac[1] ~ beta(samp_frac_prior[1], samp_frac_prior[2]);
    if (n_samp_frac > 1) {
      for (n in 2:n_samp_frac) {
        samp_frac[n] ~ normal(samp_frac[n - 1], rw_sigma); // RW
      }
    }
  }
  if (n_samp_frac > 0 && samp_frac_type != 4) { // samp_frac segmented
    for (n in 1:n_samp_frac) {
      samp_frac[n] ~ beta(samp_frac_prior[1], samp_frac_prior[2]);
    }
  }

  // data likelihood:
  if (!priors_only) { // useful to turn off for prior predictive checks
    if (contains_NAs) { // Not vectorized to 'easily' deal with NAs:
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
    } else { // No NAs; vectorized for increased efficiency:
      for (j in 1:J) {
        if (obs_model == 0) {
          daily_cases[1:last_day_obs,j] ~ poisson_log(eta[1:last_day_obs,j]);
        } else if (obs_model == 1) {
          daily_cases[1:last_day_obs,j] ~ neg_binomial_2_log(eta[1:last_day_obs,j], phi[j]);
        }
      }
    }
  }
}
generated quantities{
  real e; // renamed fsi
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
  e = fsi; // renamed
}
