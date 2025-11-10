#' Detect Heatwaves (or Cold Spells) Using heatwaveR
#'
#' This is a thin wrapper around heatwaveR that detects discrete, prolonged,
#' anomalous high-temperature events (Hobday et al., 2016). It constructs a
#' seasonal climatology and a percentile-based threshold, and then identifies
#' events that exceed the threshold for a minimum duration.
#'
#' The function requires daily data with a Date column. For general anomaly
#' detection in arbitrary time series or at other sampling resolutions, see
#' [rare_residuals()].
#'
#' @param data A data.frame (or tibble) containing a Date column and a numeric
#'   temperature column.
#' @param date_col Name of the Date column (default: "date").
#' @param value_col Name of the numeric temperature column (default: "temp").
#' @param climatology_period Length-2 character vector with start and end dates
#'   (inclusive) for the fixed long-term climatology period, e.g.,
#'   c("1981-01-01","2010-12-31"). Defaults to Hobday-style 30-year period if
#'   present in the data; otherwise the range of available dates will be used.
#' @param pctile Percentile threshold for events (default: 90 for 90th percentile).
#' @param min_duration Minimum consecutive days above threshold to define an
#'   event (default: 5).
#' @param cold_spells Logical; if TRUE, detect cold-spells (events below the
#'   threshold) instead of heatwaves (default: FALSE).
#' @param category Logical; if TRUE, also compute event categories per Hobday
#'   et al. (2018) on the detected events (default: TRUE).
#' @param ... Additional arguments passed to heatwaveR::ts2clm() or
#'   heatwaveR::detect_event().
#'
#' @return A list with class "rare_heatwaves" containing:
#'   - data: a data.frame with the full time series, seasonal climatology,
#'     threshold, and boolean indicators for threshold and event criteria
#'     (column names follow heatwaveR, but include a convenience column
#'     `is_event`).
#'   - events: a data.frame summarising detected events (start, peak, end,
#'     duration, intensities, rates, etc.).
#'   - params: a list of key parameters used (pctile, min_duration,
#'     climatology_period, cold_spells, category).
#'
#' @details
#' This function uses heatwaveR::ts2clm() to compute a seasonal climatology and a
#' percentile-based threshold over a fixed climatology window, then
#' heatwaveR::detect_event() to detect events. It assumes daily data. If your
#' data are not daily, aggregate or resample to daily before calling this
#' function to ensure meaningful climatologies and durations in days.
#'
#' @references
#' Hobday, A. J., Oliver, E. C. J., Sen Gupta, A., Benthuysen, J. A., Burrows,
#' M. T., Donat, M. G., Holbrook, N. J., Moore, P. J., Thomsen, M. S., Wernberg,
#' T., & Smale, D. A. (2016). A hierarchical approach to defining marine
#' heatwaves. Progress in Oceanography, 141, 227–238.
#'
#' Hobday, A. J., Oliver, E. C. J., Sen Gupta, A., Benthuysen, J. A., Burrows,
#' M. T., Donat, M. G., Holbrook, N. J., Moore, P. J., Thomsen, M. S., Wernberg,
#' T., & Smale, D. A. (2018). Categorizing and naming marine heatwaves.
#' Oceanography, 31(2), 162–173.
#'
#' @examples
#' \dontrun{
#'   # Example with Annapolis data, use the first 10 years for climatology
#'   data(Annapolis)
#'   summary(Annapolis)
#'
#'   # Detect heatwaves
#'   hw <- rare_heatwaves(Annapolis, date_col = "date", value_col = "tmax",
#'                        climatology_period = c("1995-01-01", "2004-12-31"),
#'                        cold_spells = FALSE,
#'                        pctile = 90, min_duration = 5)
#'
#'   # Detected events
#'   head(hw$events)
#'   tail(hw$events)
#'
#'   # Flagged days in the full series
#'   table(hw$data$event)
#'
#'   # Detect cold spells
#'   cs <- rare_heatwaves(Annapolis, date_col = "date", value_col = "tmax",
#'                     climatology_period = c("1995-01-01", "2004-12-31"),
#'                     cold_spells = TRUE,
#'                     pctile = 10, min_duration = 5)
#'
#'   # Detected events
#'   head(cs$events)
#'   tail(cs$events)
#'
#'   # Flagged days in the full series
#'   table(cs$data$event)
#' }
#'
#' @export
rare_heatwaves <- function(data,
                           date_col = "date",
                           value_col = "value",
                           climatology_period = NULL,
                           pctile = 90,
                           min_duration = 5,
                           cold_spells = FALSE,
                           # category = TRUE,
                           ...) {
  if (!requireNamespace("heatwaveR", quietly = TRUE)) {
    stop("The 'heatwaveR' package is required for rare_heatwaves(). Please install it: install.packages('heatwaveR')",
         call. = FALSE)
  }
  stopifnot(is.data.frame(data))
  # Allow column identifiers as strings or bare names
  date_col <- if (is.character(date_col)) date_col else as.character(substitute(date_col))
  value_col <- if (is.character(value_col)) value_col else as.character(substitute(value_col))
  if (!all(c(date_col, value_col) %in% names(data))) {
    stop(sprintf("Columns '%s' and/or '%s' not found in 'data'", date_col, value_col), call. = FALSE)
  }
  # Try to coerce Date column to Date class
  data[[date_col]] <- as.Date(data[[date_col]])
  # Ensure Date type for date column
  if (!inherits(data[[date_col]], "Date")) {
    stop(sprintf("Column '%s' must be of class Date or coercible via as.Date()", date_col), call. = FALSE)
  }
  # Default climatology period: prefer 30-year window if available, else full range
  d_min <- min(data[[date_col]], na.rm = TRUE)
  d_max <- max(data[[date_col]], na.rm = TRUE)
  if (is.null(climatology_period)) {
    start_year <- as.integer(format(d_min, "%Y"))
    end_year <- as.integer(format(d_max, "%Y"))
    if (end_year - start_year + 1 >= 30) {
      climatology_period <- c(sprintf("%d-01-01", start_year), sprintf("%d-12-31", start_year + 29))
    } else {
      climatology_period <- c(as.character(d_min), as.character(d_max))
    }
  }
  # Build climatology
  # heatwaveR::ts2clm() uses eval(substitute(x), data). To program with column names,
  # we must pass the literal symbols for the columns (not variables holding symbols).
  # Construct the call with do.call so that x/y are the actual column symbols.
  extra_args <- list(...)
  ts2_args <- c(list(
    data = data,
    x = as.name(date_col),
    y = as.name(value_col),
    climatologyPeriod = climatology_period,
    pctile = pctile
  ), extra_args)
  clim_df <- do.call(heatwaveR::ts2clm, ts2_args)
  # Detect events (heatwaves by default; cold_spells flips sense)
  # Normalize column names in the climatology to heatwaveR defaults to avoid NSE issues.
  # clim_df <- clim#$climatology
  nm <- names(clim_df)
  # Ensure time column 't'
  if (!("t" %in% nm)) {
    if (date_col %in% nm) {
      clim_df$`t` <- clim_df[[date_col]]
    } else {
      stop("rare_heatwaves(): Could not find a time column. Expected either 't' or the provided date_col in ts2clm output. Available columns: ",
           paste(nm, collapse = ", "))
    }
  }
  # Ensure value column 'temp'
  if (!("temp" %in% nm)) {
    if (value_col %in% nm) {
      clim_df$`temp` <- clim_df[[value_col]]
    } else {
      stop("rare_heatwaves(): Could not find a value column. Expected either 'temp' or the provided value_col in ts2clm output. Available columns: ",
           paste(nm, collapse = ", "))
    }
  }
  # seas and thresh are produced by ts2clm; if missing, error clearly
  if (!("seas" %in% names(clim_df)) || !("thresh" %in% names(clim_df))) {
    stop("rare_heatwaves(): ts2clm output is missing 'seas' and/or 'thresh' columns. Check inputs and heatwaveR version.")
  }
  # Now rely on detect_event defaults (x=t, y=temp, seasClim=seas, threshClim=thresh)
  detect_args <- c(list(
    data = clim_df,
    minDuration = min_duration,
    coldSpells = cold_spells
  ), extra_args)
  ev <- do.call(heatwaveR::detect_event, detect_args)
  # # Optionally add categories
  # if (isTRUE(category)) {
  #   ev <- heatwaveR::category(ev)
  # }
  # Convenience: boolean event flag on full series
  dat <- ev$climatology
  if (!"event" %in% names(dat) && "event" %in% names(ev$climatology)) {
    # ensure event column exists; detect_event should add it
    dat$event <- ev$climatology$event
  }
  res <- list(
    data = dat,
    events = ev$event,
    params = list(pctile = pctile,
                  min_duration = min_duration,
                  climatology_period = climatology_period,
                  cold_spells = cold_spells,
                  category = category)
  )
  class(res) <- c("rare_heatwaves", class(res))
  res
}
