#' California proposition 99
#'
#' A dataset containing per-capita cigarette consumption (in packs).
#' In year 1989 California imposed a Tobacco tax. The column `treated` is 1 from then on for California.
#'
#' @docType data
#' @name california_prop99
#'
#' @format A data frame with 1209 rows and 4 variables:
#' \describe{
#'   \item{State}{US state name, character string}
#'   \item{Year}{Year, integer}
#'   \item{PacksPerCapita}{per-capita cigarette consumption, numeric}
#'   \item{treated}{the treatmed indicator 0: control, 1: treated, numeric}
#' }
#' @source Abadie, Alberto, Alexis Diamond, and Jens Hainmueller.
#'  "Synthetic control methods for comparative case studies: Estimating the effect of California’s tobacco control program."
#'   Journal of the American statistical Association 105, no. 490 (2010): 493-505.
#'
#' @usage data(california_prop99)
#'
#' @examples
#' \donttest{
#' # Load tobacco sales in long panel format.
#' data("california_prop99")
#' # Transform to N*T matrix format required for synthdid,
#' # where N is the number of units and T the time periods.
#' setup <- panel.matrices(california_prop99)
#' }
#'
NULL

#' PENN
#'
#' @docType data
#' @name PENN
#'
#' @format A data frame with 3219 rows and 5 variables.
#' \describe{
#'   \item{country}{country}
#'   \item{year}{year}
#'   \item{log_gdp}{log_gdp}
#'   \item{dem}{dem}
#'   \item{educ}{educ}
#' }
#'
#' @usage data(PENN)
#'
NULL

#' CPS
#'
#' @docType data
#' @name CPS
#'
#' @format A data frame with 2000 rows and 8 variables.
#' \describe{
#'   \item{state}{state}
#'   \item{year}{year}
#'   \item{log_wage}{log_wage}
#'   \item{hours}{hours}
#'   \item{urate}{urate}
#'   \item{min_wage}{min_wage}
#'   \item{open_carry}{open_carry}
#'   \item{abort_ban}{abort_ban}
#' }
#'
#' @usage data(CPS)
#'
NULL
