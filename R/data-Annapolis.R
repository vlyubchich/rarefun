#' Annapolis daily Daymet time series (1995-2024)
#'
#' Daily maximum temperature and precipitation for a single Daymet pixel near
#' Annapolis, Maryland (USA). These data were retrieved via the package \code{daymetr}
#' and pared down for compact examples in this package.
#' \itemize{
#'   \item Source: Daymet (ORNL DAAC / NASA ESDIS) via the package \code{daymetr}.
#'   \item Units: \code{tmax} in degrees Celsius; \code{pcp} in millimeters per day.
#'   \item Coverage: \code{1995-01-01} through \code{2024-12-31} (inclusive; 365-day years, excluding 29 February), one location (Annapolis, MD; approx. 38.9784° N, 76.4922° W).
#' }
#'
#' @format A tibble (data frame) with the following columns:
#' \itemize{
#'   \item \code{date} Date. Calendar date.
#'   \item \code{tmax} numeric. Daily maximum 2 m air temperature (°C).
#'   \item \code{pcp} numeric. Daily total precipitation (mm/day).
#' }
#'
#' @source Thornton, M. M., Shrestha, R., Wei, Y., Thornton, P. E., & Kao, S.-C. (2022). 
#' Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 4 R1 (Version 4.1). 
#' ORNL Distributed Active Archive Center. \url{https://doi.org/10.3334/ORNLDAAC/2129} Date Accessed: 2025-10-31.
#'
#' @docType data
#' @keywords datasets
#' @name Annapolis
#' @examples
#' data(Annapolis)
#' head(Annapolis)
#'
#' # Quick check of data
#' summary(Annapolis)
#'
#' # Minimal plot (if graphics device available)
#' plot(Annapolis$date, Annapolis$tmax, type = "l",
#'      xlab = "Date", ylab = "Tmax (°C)",
#'      las = 1)
"Annapolis"
