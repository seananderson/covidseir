
#' function create vaccine sequence based on schedule
#' @description Function converts a vaccine schedule into a list for each age group
#'              with the vaccination group scheduled in that age. Vaccination
#'              schedule can either be defined for each age group as a vector representing
#'              the start and end points of the vaccine schedule or as a list representing a
#'              series of vaccine roll-outs for that age group with their start end date and proportion
#'              of that population vaccinated.
#' @param vaccination_schedule list of age_groups. Each item in list is either:
#' \describe{
#'    \item{`NULL`}{Group will not be vaccinated}
#'    \item{`character`}{vector of two characters denoting the start and end
#'                       times in ymd format that age-group is vaccinated}
#'    \item{`list`}{numeric list where each item has a `start_date`, an `end_date`,
#'                  and a `proportion` denoting the proportion of that age-group
#'                  to vaccinate at that time. }
#' }
#' @param population_sizes list of age_groups with numeric values denoting
#'   population size in that age group
#' @param min_date string for min date
#' @param max_date string in ymd for max date
#' @param hesitancy numeric denoting proportion who will not take vaccine
#' @param immunity_delay Delay between vaccination and immunological impacts
#'   (in days)
#' @examples
#' # example vaccination schedule that can be used.
#' vaccination_schedule <- list(
#'   "< 2" = NULL,
#'   "2 - 5" = NULL,
#'   "6 - 17" = NULL,
#'   "18 - 24" = c("2021-09-01", "2021-10-01"),
#'   "25 - 34" = c("2021-08-01", "2021-09-01"),
#'   "35 - 44" = c("2021-07-01", "2021-08-01"),
#'   "45 - 54" = c("2021-07-01", "2021-08-01"),
#'   "55 - 64" = c("2021-06-01", "2021-08-01"),
#'   "65 - 74" = c("2021-05-01", "2021-06-01"),
#'   "> 75" = c("2021-03-01", "2021-05-01")
#' )
#' @return  list of all age groups with a vector of length min_date to max_date
#' where each day denotes the vaccination rate
#' @export
create_vaccine_seqs <- function(vaccination_schedule, population_sizes,
                                min_date = "2021-02-01", max_date = "2021-10-01",
                                hesitancy = 0.15, immunity_delay = 14) {
  dates <- seq(lubridate::ymd(min_date), lubridate::ymd(max_date), by = "day")
  vaccination_seqs <- list()
  for (age_group in names(vaccination_schedule)) {
    if (typeof(vaccination_schedule[[age_group]]) == "character") {

      # age group vaccine start date (adjusting for dealy in immunity build up)
      start_date <- lubridate::ymd(vaccination_schedule[[age_group]][1]) +
        lubridate::days(immunity_delay)
      end_date <- lubridate::ymd(vaccination_schedule[[age_group]][2]) +
        lubridate::days(immunity_delay)

      # get population size
      pop_size <- population_sizes %>%
        dplyr::filter(age_group == !!age_group) %>%
        dplyr::pull(size)

      # create when vaccination is occurring
      vaccine_rollout <- (dates >= start_date) & (dates <= end_date)
      # total scheduled vaccine days
      total_days <- sum(vaccine_rollout)
      # vaccine roll-out is total population size with hesitants removed divided among all vaccine days
      vaccine_rollout <- vaccine_rollout * (pop_size / total_days) * (1 - hesitancy)

      # age group vaccine end date

      vaccination_seqs[[age_group]] <- vaccine_rollout
    } else if (typeof(vaccination_schedule[[age_group]]) == "list") {

      # get list of vaccine schedule for specific age_group
      age_group_schedule <- vaccination_schedule[[age_group]]

      # get population size
      pop_size <- population_sizes %>%
        dplyr::filter(age_group == !!age_group) %>%
        dplyr::pull(size)

      # create empty vaccine roll-out to populate with age-specific schedule
      vaccine_rollout <- rep(0, length(dates))

      for (i in 1:length(age_group_schedule)) {
        start_date <- age_group_schedule[[i]]$start_date
        end_date <- age_group_schedule[[i]]$end_date
        proportion <- age_group_schedule[[i]]$proportion

        # age group vaccine start date (adjusting for dealy in immunity build up)
        start_date <- lubridate::ymd(start_date) + lubridate::days(immunity_delay)
        end_date <- lubridate::ymd(end_date) + lubridate::days(immunity_delay)

        # create when vaccination is occurring
        vaccine_days <- (dates >= start_date) & (dates <= end_date)
        # total scheduled vaccine days
        total_days <- sum(vaccine_days)
        # vaccine roll-out is total population size with hesitants removed divided among all vaccine days
        vaccine_rollout <- vaccine_rollout + proportion * vaccine_days * (pop_size / total_days) * (1 - hesitancy)
      }

      vaccination_seqs[[age_group]] <- vaccine_rollout
    } else if (is.null(vaccination_schedule[[age_group]])) {
      vaccination_seqs[[age_group]] <- rep(0, length(dates))
    }
  }

  return(vaccination_seqs)
}

#' function create vaccine sequence
#' @param vaccination_schedule list of age_groups either NULL if not vaccination
#'   in group or vector of length 2 containing string of dates in ymd format
#'   expressing the start and stop of vaccination
#' @param population_sizes tibble of age_groups with numeric values denoting
#'   population size in that age group
#' @param age_susceptibility tibble with columns age_group and susceptibility.
#'   Denotes each age group susceptible to infection
#'   thus effectively weights the contacts
#' @param contact_rates tibble of age_groups with numeric values denoting
#'   contact rates in that age group
#' @param min_date string for min date
#' @param max_date string in ymd for max date
#' @param hesitancy numeric denoting proportion who will not take vaccine
#' @param immunity_delay Delay between vaccination and immunological impacts
#'   (in days)
#' @return  list of all age groups with a vector of length min_date to max_date
#'   where each day denotes the vaccination rate
#' @export
create_adjusted_vaccination_rollout <- function(vaccination_schedule, population_sizes,
                                                age_susceptibility, contact_rates,
                                                min_date = "2021-02-01",
                                                max_date = "2021-10-01",
                                                hesitancy = 0.15,
                                                immunity_delay = 14,
                                                tb_eff = 0.8) {

  # get vaccination sequence
  vaccination_seqs <- create_vaccine_seqs(vaccination_schedule, population_sizes,
    min_date = min_date,
    max_date = max_date,
    hesitancy = hesitancy,
    immunity_delay = immunity_delay
  )

  # create empty array of correct size
  total_vac_rate <- rep(0, length(vaccination_seqs[[1]]))

  # get total activity (of total population)
  total_contacts <- sum(population_sizes$size * contact_rates$contacts * age_susceptibility$susceptibility)

  # get total population size
  total_population_size <- sum(population_sizes$size)

  for (age_group in names(vaccination_seqs)) {

    # get population size
    pop_size <- population_sizes %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(size)

    # get susceptibility
    susceptibility <- age_susceptibility %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(susceptibility)

    # get contacts
    contacts <- contact_rates %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(contacts)

    # get adjusted vaccine rate for age group
    adj_vac <- total_population_size * vaccination_seqs[[age_group]] * contacts * susceptibility / total_contacts

    # adjust to reduce due to imperfect vaccine
    adj_vac <- tb_eff * adj_vac

    # add to total vaccination rate
    total_vac_rate <- total_vac_rate + adj_vac
  }

  return(total_vac_rate)
}

#' get proportion not vaccinated in each age group
#' @inheritParams create_adjusted_vaccination_rollout
#' @return  list of all age groups with a vector of length min_date to max_date
#' where each day denotes the proportion not vaccinated
#' @noRd
get_prop_not_vaccinated <- function(vaccination_seqs, population_sizes) {
  prop_not_vaccinated <- list()
  for (age_group in names(vaccination_seqs)) {
    # get population size
    pop_size <- population_sizes %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(size)

    # get total not vac to given date
    total_not_vac <- pop_size - cumsum(vaccination_seqs[[age_group]])

    prop_not_vaccinated[[age_group]] <- total_not_vac / pop_size
  }
  return(prop_not_vaccinated)
}


#' function create outcome rate adjusted for who is vaccinated at a given time
#' Outcomes can either be hospitalizations or deaths, update the outcome rate
#' parameter to change this
#' @param outcome_rates tibble of age_groups and rate denoting either hosp or
#'   deaths by age group
#' @inheritParams create_adjusted_vaccination_rollout
#' @param ve_outcome Vaccine efficacy for outcome
#' @return  list of all age groups with a vector of length min_date to max_date
#' where each day denotes the vaccination rate
#' @export
create_adjusted_outcome_rate <- function(outcome_rates, vaccination_schedule,
                                         population_sizes,
                                         age_susceptibility, contact_rates,
                                         min_date = "2021-02-01",
                                         max_date = "2021-10-01",
                                         hesitancy = 0.15, immunity_delay = 14,
                                         ve_outcome = 0.95) {

  # get vaccination sequence
  vaccination_seqs <- create_vaccine_seqs(vaccination_schedule, population_sizes,
    min_date = min_date, max_date = max_date,
    hesitancy = hesitancy,
    immunity_delay = immunity_delay
  )

  # create empty array of correct size
  total_rate <- rep(0, length(vaccination_seqs[[1]]))

  # get proportion not vaccinated for each age group
  prop_not_vaccinated <- get_prop_not_vaccinated(vaccination_seqs, population_sizes)

  # get total activity (of total population)
  total_contacts <- sum(population_sizes$size * contact_rates$contacts * age_susceptibility$susceptibility)

  # get total population size
  total_population_size <- sum(population_sizes$size)

  for (age_group in names(vaccination_seqs)) {
    # get population size
    pop_size <- population_sizes %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(size)

    # get susceptibility
    susceptibility <- age_susceptibility %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(susceptibility)

    # get contacts
    contacts <- contact_rates %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(contacts)

    # get total not vac to given date
    prop_not_vac <- prop_not_vaccinated[[age_group]]

    # get rate
    outcome_rate <- outcome_rates %>%
      dplyr::filter(age_group == !!age_group) %>%
      dplyr::pull(rate)

    # adjusted rate composition of prop vaccinated and vac efficacy
    adjusted_outcome_rate <- outcome_rate * prop_not_vac +
      (1 - ve_outcome) * (1 - prop_not_vac) * outcome_rate

    # get adjusted vaccine rate for age group
    adj_rate <- adjusted_outcome_rate * pop_size * contacts * susceptibility / total_contacts

    # add to total vaccination rate
    total_rate <- total_rate + adj_rate
  }

  return(total_rate)
}
