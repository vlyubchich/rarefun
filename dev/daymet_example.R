# Script to create a small Daymet-derived daily time series for examples
#
# Do not run this script during R CMD check.
# It is intended to be run manually to regenerate the dataset in data/.
#
# Requirements:
# - Packages: daymetr, dplyr
# - Internet access to ORNL DAAC Daymet web services
#
# Dataset goal: 30+ years of daily data for a single poin: daily maximum air temperature
#   (tmax) and precipitation (pcp) to keep package size small.
#
# Provenance and citation:
# - Data source:
#   Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 4 R1 https://doi.org/10.3334/ORNLDAAC/2129
#
# Re-run steps:
# - Set a point of interest (lat/lon) and a date range available in Daymet.
# - Download single-pixel data, simplify to tibble, keep date + tmax + pcp.
# - (Optionally) compute tmean from tmax/tmin.
# - Save compressed dataset to data/.

library(daymetr)
library(dplyr)
library(tidyr)

# Choose Annapolis, Maryland
site <- "Annapolis_MD"
lat <- 38.9784
lon <- -76.4922
end_year <- 2024
start_year <- end_year - 45 + 1  # 45 years total

# Download single-pixel Daymet data as a tibble
# This returns a nested list; with simplify=TRUE we get a tibble in $data
dm0 <- daymetr::download_daymet(site = site, lat = lat, lon = lon,
                                start = start_year, end = end_year,
                                internal = TRUE, simplify = TRUE, silent = FALSE, force = FALSE)

# The tibble columns are standardized; compute Date and keep only needed fields
# Daymet returns year and yday (day-of-year). Create a proper Date column.
Annapolis <- dm0 %>%
    tidyr::pivot_wider(names_from = measurement, values_from = value) %>%
    dplyr::mutate(date = as.Date(paste(year, yday, sep = "-"), format = "%Y-%j")) %>%
    dplyr::rename(pcp = prcp..mm.day.,
                  tmax = tmax..deg.c.) %>%
    dplyr::select(date, pcp, tmax)

# Save to data/ as an internal example dataset for the package
usethis::use_data(Annapolis, overwrite = TRUE, compress = "xz")
