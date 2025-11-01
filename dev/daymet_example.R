# Script to create a small Daymet-derived daily time series for examples
#
# Do not run this script during R CMD check.
# It is intended to be run manually to regenerate the dataset in data/.
#
# Requirements:
# - Packages: daymetr, dplyr
# - Internet access to ORNL DAAC Daymet web services
#
# Dataset goal:
# - 30-year daily series for a single point, keeping only daily maximum air temperature
#   (tmax) and precipitation (pcp) to keep package size small.
# - Columns: date (Date), temp (numeric), lat, lon, site (chr)
# - Saved as data/daymet_tmax_1990_2019.rda (compress = "xz")
#
# Provenance and citation:
# - Data source: Daymet (ORNL DAAC / NASA ESDIS). Please cite per
#   https://daymet.ornl.gov/citations and include dataset DOI if applicable.
# - NASA ESDIS open data guidance: https://earthdata.nasa.gov/earth-observation-data/data-use-policy
#
# Re-run steps:
# - Set a point of interest (lat/lon) and a 30-year range available in Daymet.
# - Download single-pixel data, simplify to tibble, keep date + tmax.
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
start_year <- end_year - 29  # 30 years total

# Download single-pixel Daymet data as a tibble
# This returns a nested list; with simplify=TRUE we get a tibble in $data
dm0 <- daymetr::download_daymet(site = site, lat = lat, lon = lon,
                                start = start_year, end = end_year,
                                internal = TRUE, simplify = TRUE, silent = FALSE, force = FALSE)

# The tibble columns are standardized; compute Date and keep only needed fields
# Daymet returns year and yday (day-of-year). Create a proper Date column.
Annapolis <- dm0 %>%
    tidyr::pivot_wider(names_from = measurement, values_from = value) %>%
    # To create Date without Feb 29, first use a non-leap year, then replace with actual year
    dplyr::mutate(date = as.Date(paste(2001, yday, sep = "-"), format = "%Y-%j")) %>%
    dplyr::mutate(date = as.Date(paste0(year, "-", format(date, "%m-%d")))) %>%
    dplyr::rename(pcp = prcp..mm.day.,
                  tmax = tmax..deg.c.) %>%
    dplyr::select(date, pcp, tmax)

# Save to data/ as an internal example dataset for the package
usethis::use_data(Annapolis, overwrite = TRUE, compress = "xz")
