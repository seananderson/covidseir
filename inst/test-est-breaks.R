dat <- structure(list(value = c(1, 1, 4, 6, 4, 6, 1, 3, 5, 11, 20, 12,
  25, 12, 27, 46, 41, 60, 47, 62, 46, 43, 88, 59, 91, 62, 62, 37,
  30, 82, 56, 47, 55, 31, 21, 13, 51, 40, 33, 34, 41, 28, 15, 30,
  44, 14, 51, 25, 30, 15, 34, 63, 33, 37, 72, 61, 15, 58, 41, 26,
  29, 23, 33, 21, 16, 24, 29, 21, 14, 8, 17, 10, 13, 16, 14, 21,
  8, 11, 5, 19, 11, 22, 12, 4, 9, 8, 7, 10, 5, 11, 9, 16, 4, 23,
  2, 4, 12, 7, 9, 10, 14, 12, 17, 14, 14, 9, 11, 19, 9, 5, 9, 6,
  17, 11, 12, 21, 12, 9, 13, 3, 12, 17, 10, 9, 12, 16, 7, 10, 19,
  18, 22, 25, 23, 20, 13, 19, 24, 33, 44, 19, 29, 33, 34, 32, 30,
  30, 21, 21, 24, 47, 26, 39, 52, 29, 48, 28, 43, 46, 47, 51, 38,
  42, 46, 86, 73, 90, 100, 86, 51, 77, 65, 84, 90, 115, 79, 81),
  day = 1:176), row.names = c(NA, -176L),
  class = "data.frame")
dat <- dplyr::as_tibble(dat)
dat$date <- lubridate::ymd("2020-03-01") + dat$day - 1
samp_frac <- c(rep(0.14, 13), rep(0.21, 40 - 13), rep(0.21, 11))
samp_frac <- c(samp_frac, rep(0.37, nrow(dat) - length(samp_frac)))

plot(dat$value)

f_seg <- c(0, rep(1, nrow(dat) - 1))
day_new_f <- which(dat$date == lubridate::ymd("2020-05-01"))
f_seg[seq(day_new_f, length(f_seg))] <- 1
day_ch <- which(dat$date == lubridate::ymd("2020-06-01"))
f_seg[seq(day_ch, length(f_seg))] <- 2
f_seg

fit <- covidseir::fit_seir(
  daily_cases = dat$value,
  samp_frac_fixed = samp_frac,
  f_seg = f_seg,
  time_increment = 1,
  i0_prior = c(log(8), 1),
  start_decline_prior = c(log(15), 0.1),
  end_decline_prior = c(log(22), 0.1),
  f_break_prior = rbind(c(log(80), 0.1)),
  f_prior = cbind(c(0.4, 0.6), c(0.2, 0.2)),
  R0_prior = c(log(2.6), 0.2),
  N_pop = 5.1e6,
  iter = 140,
  chains = 1,
  fit_type = "optimizing"
)
fit


# ON:

make_f_seg <- function(.dat, .date = "2020-05-01") {
  f_seg <- c(0L, rep(1L, nrow(.dat) - 1))
  day_new_f <- which(.dat$date == lubridate::ymd(.date))
  f_seg[seq(day_new_f, length(f_seg))] <- 2L
  f_seg
}

dat <- structure(list(date = structure(c(18322, 18323, 18324, 18325,
  18326, 18327, 18328, 18329, 18330, 18331, 18332, 18333, 18334,
  18335, 18336, 18337, 18338, 18339, 18340, 18341, 18342, 18343,
  18344, 18345, 18346, 18347, 18348, 18349, 18350, 18351, 18352,
  18353, 18354, 18355, 18356, 18357, 18358, 18359, 18360, 18361,
  18362, 18363, 18364, 18365, 18366, 18367, 18368, 18369, 18370,
  18371, 18372, 18373, 18374, 18375, 18376, 18377, 18378, 18379,
  18380, 18381, 18382, 18383, 18384, 18385, 18386, 18387, 18388,
  18389, 18390, 18391, 18392, 18393, 18394, 18395, 18396, 18397,
  18398, 18399, 18400, 18401, 18402, 18403, 18404, 18405, 18406,
  18407, 18408, 18409, 18410, 18411, 18412, 18413, 18414, 18415,
  18416, 18417, 18418, 18419, 18420, 18421, 18422, 18423, 18424,
  18425, 18426, 18427, 18428, 18429, 18430, 18431, 18432, 18433,
  18434, 18435, 18436, 18437, 18438, 18439, 18440, 18441, 18442,
  18443, 18444, 18445, 18446, 18447, 18448, 18449, 18450, 18451,
  18452, 18453, 18454, 18455, 18456, 18457, 18458, 18459, 18460,
  18461, 18462, 18463, 18464, 18465, 18466, 18467, 18468, 18469,
  18470, 18471, 18472, 18473, 18474, 18475, 18476, 18477, 18478,
  18479, 18480, 18481, 18482, 18483, 18484, 18485, 18486, 18487,
  18488, 18489, 18490, 18491, 18492, 18493, 18494, 18495, 18496,
  18497, 18498, 18499, 18500, 18501, 18502, 18503, 18504, 18505,
  18506, 18507, 18508, 18509, 18510, 18511, 18512), class = "Date"),
  day = 1:191, value = c(4, 0, 5, 0, 2, 6, 0, 4, 3, 2, 5, 17,
    20, 24, 42, 32, 12, 25, 44, 60, 59, 48, 78, 85, 100, 170,
    225, 241, 301, 441, 350, 437, 437, 398, 335, 420, 404, 546,
    447, 559, 407, 489, 433, 494, 560, 522, 519, 626, 561, 531,
    675, 492, 712, 577, 490, 428, 487, 456, 458, 463, 391, 552,
    438, 243, 549, 380, 449, 481, 430, 320, 315, 362, 344, 357,
    359, 393, 373, 330, 308, 502, 439, 404, 497, 470, 338, 419,
    326, 318, 383, 371, 327, 372, 464, 413, 352, 343, 416, 251,
    473, 312, 279, 259, 237, 212, 246, 175, 207, 190, 196, 188,
    209, 181, 157, 224, 238, 201, 178, 115, 164, 252, 238, 151,
    124, 159, 165, 116, 124, 184, 84, 164, 121, 119, 81, 66,
    256, 83, 104, 121, 143, 199, 147, 154, 187, 136, 143, 208,
    128, 104, 150, 125, 93, 102, 118, 109, 74, 70, 126, 75, 86,
    91, 77, 59, 135, 91, 73, 72, 103, 82, 48, 131, 117, 91, 97,
    142, 78, 92, 157, 86, 125, 138, 116, 114, 88, 156, 135, 140,
    150, 128, 157, 126, 185)), row.names = c(NA, -191L), class = c("tbl_df",
      "tbl", "data.frame"))

dat <- dat[1:170,]
plot(dat$value, type = "l")
f_on <- make_f_seg(dat, .date = "2020-05-01")
day_ch <- which(dat$date == lubridate::ymd("2020-06-01"))
day_ch2 <- which(dat$date == lubridate::ymd("2020-07-06"))
day_ch3 <- which(dat$date == lubridate::ymd("2020-07-15"))
day_ch4 <- which(dat$date == lubridate::ymd("2020-08-01"))
f_on[seq(day_ch, day_ch2 - 1)] <- 3L
f_on[seq(day_ch2, length(dat$value))] <- 4L
f_on[seq(day_ch3, length(dat$value))] <- 5L
# f_on[seq(day_ch4, length(dat$value))] <- 6L
f_on[-c(1:10)] <- 2
f_on[length(f_on)] <- 3
f_on
plot(dat$value, type = "l")
abline(v = 8)
abline(v = 23)
abline(v = 23 + 47)
abline(v = 47 + 50)

fit_on <- covidseir::fit_seir(
  daily_cases = dat$value,
  samp_frac_fixed = rep(0.2, nrow(dat)),
  i0_prior = c(log(1), 1),
  start_decline_prior = c(log(8), 0.1),
  end_decline_prior = c(log(23), 0.1),
  f_break_prior = rbind(c(log(30), 0.1), c(log(35), 0.1)),
  N_pop = 14.5e6,
  f_seg = f_on,
  f_prior = rbind(c(0.4, 0.2), c(0.4, 0.2), c(0.4, 0.2)),
  iter = 400,
  # ode_control = c(1e-8, 1e-7, 1e7),
  time_increment = 1,
  fit_type = "optimizing"
)
fit_on

# hist(fit_on$post$f_breaks[,1])
# hist(fit_on$post$f_breaks[,2])
# hist(fit_on$post$end_decline + fit_on$post$f_breaks[,1])
# hist(fit_on$post$end_decline + fit_on$post$f_breaks[,1] + fit_on$post$f_breaks[,2])

library(ggplot2)
p <- project_seir(fit_on, iter = 1:100)
tidy_seir(p) %>%
  plot_projection(obs_dat = dat) +
  geom_vline(data =
      data.frame(day = fit_on$post$end_decline + fit_on$post$f_breaks[,1]),
    aes(xintercept = day), alpha = 0.05, col = "red") +
  geom_vline(data =
      data.frame(day = fit_on$post$end_decline + fit_on$post$f_breaks[,1] + fit_on$post$f_breaks[,2]),
    aes(xintercept = day), alpha = 0.05, col = "purple")

