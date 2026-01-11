#' The `codyna` Package.
#'
#' @description Performs analysis of complex dynamic systems with a focus on the
#' temporal unfolding of patterns, changes, and state transitions in
#' behavioral data. The package supports both time series and sequence data
#' and provides tools for the analysis and visualization of complexity,
#' pattern identification, trends, regimes, sequence typology as well as
#' early warning signals.
#'
#' @author Santtu Tikka and Mohammed Saqr
#'
"_PACKAGE"

#' Example Data on Student Engagement
#'
#' Students' engagement states (Active / Average / Disengaged)
#' throughout a whole study program. The data was generated synthetically
#' based on the article "The longitudinal association between engagement and
#' achievement varies by time, students' profiles, and achievement state:
#' A full program study". Used also in the `tna` package.
#'
#' @encoding UTF-8
#' @source \doi{10.1016/j.compedu.2023.104787}
#' @format An `stslist` object (sequence data).
#' @references
#' Tikka S, López-Pernas S, Saqr M (2025).
#' "tna: An R Package for Transition Network Analysis."
#' _Applied Psychological Measurement_. \doi{10.1177/01466216251348840}
#'
"engagement"

#' Example Data on Group Regulation
#'
#' Students' regulation during collaborative learning. Students' interactions
#' were coded as: "adapt", "cohesion", "consensus", "coregulate", "discuss",
#' "emotion", "monitor", "plan", "synthesis". Used also in the `tna` package.
#'
#' @encoding UTF-8
#' @source The data was generated synthetically.
#' @format A `data.frame` object.
#' @references
#' Tikka S, López-Pernas S, Saqr M (2025).
#' "tna: An R Package for Transition Network Analysis."
#' _Applied Psychological Measurement_. \doi{10.1177/01466216251348840}
#'
"group_regulation"


#' Ecological Momentary Assessment (EMA) Data
#'
#' Example data for complex adaptive systems perspective to behavior change
#' research. The dataset consists of 20 individuals with 9 self-report
#' variables (and time of response) each. For more information on the data,
#' please see
#' \url{https://heinonmatti.github.io/complexity-behchange/dataset-info.html}
#'
#' @source \url{https://github.com/heinonmatti/complexity-behchange}
#' @format A `data.frame` object.
#'
"ema"
