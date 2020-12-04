# covidseir

# covidseir 0.0.1.9001

* Add optional projection plot coloured by Rt.

# covidseir 0.0.1.9000

* Added `post2prior()`, which allows passing the posterior from one model fit
  into priors of a second model fit. This allows the model to be fit in blocks.
  There is a new vignette showing a full example.

* Removed the `f_ramp_rate` argument and option from `fit_seir()`.

* Allow not fitting the initial social distancing ramp and allow specifying the
  NB2 phi prior as a log normal distribution to allow passing the posterior as a
  prior.

# covidseir 0.0.0.9010

* Decrease default `time_increment` to 0.25 and make ODE control
  more precise again after issues in leeway paper.

* Add initialization option (MAP or random).

* Stop making projections parallel for now (possible issue in leeway paper).

# covidseir 0.0.0.9009

* Add estimation of fraction distancing (`e`) back into the model.

# covidseir 0.0.0.9007

* Add a vignette

* Add VB and MAP algorithm options

# covidseir 0.0.0.9006

* Change default `time_increment` to 0.5 now that the pandemic
  has been going on for many months. You could even increase this
  to 1 for more speed.

* Remove estimation of fraction distancing (e and e_prior). There isn't
  generally information to inform this parameter and it slows everything down.

# covidseir 0.0.0.9005

* Fix a bug in the imported cases during the projection period.

# covidseir 0.0.0.9004

* Add the ability to incorporate imported cases during the projection period.

# covidseir 0.0.0.9003

* Some of the default arguments have changed to accommodate the new
  estimated parameters and to make the model easier to work with. In particular,
  population size is now specified via the `N_pop` argument and `i0` is now
  `i0_prior` and requires a vector of length 2 to specify the prior.

* Add estimation of the fraction distancing (parameter e). In most cases there
  will not be information in the data to inform this parameter; however,
  the prior allows for the inclusion of uncertainty on this fraction.

* Add estimation of i0, the total number of infected people at the initial
  point in time (default -30 days), rather than requiring it to be a fixed value.

* Add estimation of the starting and ending dates of the initial ramp in of
  physical_distancing.

# covidseir 0.0.0.9002

* `forecast_seir()` has been renamed to `project_seir()`.

* `plot_projection()` and `tidy_sier()` have been added. See examples
  in `project_seir()`.

# covidseir 0.0.0.9001

* `fit_seir()` now allows for multiple estimated blocks (segments) of fractions
  of normal contacts (f) through time. The f values are indexed by `s` for
  segment. I.e., `f_s`. They currently all share the same prior.

* `fit_seir()` now allows for estimated blocks (segments) of sampling fractions
  for the first data type. The intent is that this is used for the reported
  cases in which the sampling fraction is unknown and can be estimated if
  another type of case data (e.g., hospitalizations) is included with a much
  better-known fraction of positive cases represented.

* Forecasting is now done with the same Stan model as the fitting via
  `forecast_seir()`. This eliminates approximately half of the code base (the
  R reimplementation), means that all of the model fitting options also work
  with the forecasting, makes it much easier to add functionality and harder to
  make mistakes, and generates much faster forecasts. The f values are
  represented as a vector and so can represent any desired pattern.

* `fit_seir()` now allows for NA values in the case data. These are omitted
  from the likelihood. If no NA values are present then a vectorized likelihood
  is used, which is slightly faster.

* The output from `forecast_seir()` includes helpful metadata such as which
  time steps are forecasted and which time steps use a fixed (vs. estimated)
  f value.

* `fit_seir()` now uses a time step of 0.2 by default, which seems to give
  similar results with faster fits.

* Fixed f values can be set to start on an arbitrary date in the future.

* `project_fit()` has been renamed to `forecast_seir()`.

* Arguments to the functions and within the Stan model have been standardized
  and cleaned up.

* Added continuous integration testing via Travis
  <https://travis-ci.org/github/seananderson/covidseir/builds>.

* Added many examples to `fit_seir()` and `forecast_seir()`. Removed the
  vignette for now as it was out of date and is currently replaced by those
  examples.

* The previous plotting function has been removed for now. It was bloated and
  included many dependencies and wasn't working with all options. A simpler
  version will likely be added back later.

* Documentation has been cleaned up throughout.


