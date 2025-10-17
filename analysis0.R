rm(list = ls())

# Packages and functions ----
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)

library(ggplot2)
theme_set(theme_light())
library(patchwork)
library(plot.matrix)

library(geosphere)
library(ncdf4)
library(readr)
library(isotree)
library(dbscan)
library(lubridate)
library(quantmod)

# List files in /code_fun starting with "fun_" and ending with ".R"
r_files <- list.files(
    path = "code_fun",
    pattern = "^fun_.*\\.R$",
    full.names = TRUE
)
r_files

# Source each file
# for (file in r_files) {
#     source(file)
# }
source("code_fun/fun_match_spatial_points.R")
source("code_fun/fun_extract_ts.R")
source("code_fun/fun_kneedle.R")
source("code_fun/fun_rare_dbscan.R")
source("code_fun/fun_rare_iforest.R")
source("code_fun/fun_rare_residuals.R")
source("code_fun/fun_contingency_chisq.R")


# Data ----

## SeineSurvey_sbass ----

# Load and bind datasets
# this is an efficient replacement of multiple data.table::fread()
SeineSurvey <- rbindlist(lapply(c("dataderived/MDNRseineSurveyRE_sbassAM_revJun20.csv",
                                  "dataderived/VIMSseineCPUERE_sbassAM_revJul8.csv")
                                , fread), fill = TRUE)
SeineSurvey <- setnames(SeineSurvey, c("STATION", "LAT", "LON", "YEAR", "MONTH", "DAY", "CPUE"),
                        c("id", "lat", "lon", "year", "month", "day", "value"))


# Filter by species
SeineSurvey_sbass <- SeineSurvey[SPP_AGE == "STRIPED BASS YOY" & STATIONTYPE == "F"][, variable := "sbass_CPUE"]

# Calculate data availability
presence_absence <- SeineSurvey_sbass %>%
    group_by(year, id) %>%
    summarise(presence = any(!is.na(value))) %>%
    pivot_wider(names_from = id, values_from = presence, values_fill = 0)
presence_absence_matrix <- as.matrix(presence_absence[, -1]) # Remove the year column for the matrix
rownames(presence_absence_matrix) <- presence_absence$year
availability <- apply(presence_absence_matrix, 2, mean)
summary(availability)

# Plot presence/absence matrix
# plot(t(presence_absence_matrix),
#      las = 1,
#      xlab = "Year", ylab = "")
# Reshape to long format
presence_absence_long <- presence_absence %>%
    pivot_longer(cols = -year, names_to = "id", values_to = "Presence") %>%
    mutate(Presence = factor(Presence,
                             levels = c(FALSE, TRUE),
                             labels = c("Absent", "Present")),
           id = factor(id, levels = unique(id)))

# Create the heatmap using ggplot2
ggplot(presence_absence_long, aes(x = factor(year), y = id, fill = Presence)) +
    geom_tile(color = "white", linewidth = 0.5) + # Add white borders to tiles
    scale_fill_manual(values = c("Absent" = "lightgray", "Present" = "steelblue")) +
    labs(
        x = "Year",
        y = "Seine survey station ID",
        title = "Presence/absence of catches across years"
    ) +
    theme_minimal() + # Use a minimal theme
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_blank()
    )
# Save
ggsave("images/SeineSurvey_sbass_presence_absence.png",
       bg = "white",
       width = 10, height = 12)

# Keep only stations with at least 80% of non-missing data
id2keep <- names(availability)[availability >= 0.8]
SeineSurvey_sbass <- SeineSurvey_sbass[id %in% id2keep]

# Remove time dimension by selecting unique locations, rename to c("id", "lat", "lon")
SeineSurvey_sbass_stations <- unique(SeineSurvey_sbass[, .(id, lat, lon)])


## Daymet ----

# daymet_tmp <- data.table::fread("./dataderived/daymet/daymet_prcp_2021.csv",
#                                 na.strings = c("NA", "NaN", "N/A", ""))
# daymet_tmp <- daymet_tmp[DoY == 1, .(lat = Lat, lon = Lon)][, id := 1:.N]


## CBEFS ----
# Load the output of "1_CBEFS.R"
CBEFS_tmp <- data.table::fread("./dataderived/CBEFS_SampleCBEFS.csv",
                               na.strings = c("NA", "NaN", "N/A", ""))
# level 1 = surface
# level 20 = bottom



# Spatial matches ----
match_SeineSurvey_sbass_daymet <- match_spatial_points(SeineSurvey_sbass_stations,
                                                       daymet_tmp,
                                                       fast = FALSE)
match_SeineSurvey_sbass_CBEFS <- match_spatial_points(SeineSurvey_sbass_stations,
                                                      CBEFS_tmp,
                                                      fast = FALSE)
readr::write_csv(match_SeineSurvey_sbass_CBEFS,
                 "dataderived/match_SeineSurvey_sbass_CBEFS.csv")

# Plots to check
# hist(match_SeineSurvey_sbass_CBEFS$distance_m, 30)
ggplot(match_SeineSurvey_sbass_CBEFS, aes(x = distance_m)) +
    geom_histogram(binwidth = 500, fill = "steelblue", color = "black", boundary = 0) +
    labs(
        x = "Distance, m",
        y = "Frequency"
    )

ggplot(CBEFS_tmp, aes(x = lon, y = lat)) +
    geom_point() +
    geom_point(data = SeineSurvey_sbass_stations, aes(lon, lat), col = "cornflowerblue") +
    geom_point(data = SeineSurvey_sbass_stations %>% dplyr::filter(id == 106),
               aes(lon, lat), col = "red", pch = 19, cex = 3)
ggsave("images/map_SeineSurvey_sbass_CBEFS.png", width = 3, height = 4)



# Extract ts ----
# Output a data frame with columns: "id", "lat", "lon", "time", "variable", "value"

# SeineSurvey_sbass is fast, do it here
ts_SeineSurvey_sbass <- SeineSurvey_sbass %>%
    dplyr::filter(id %in% match_SeineSurvey_sbass_CBEFS$id_SeineSurvey_sbass_stations) %>%
    dplyr::mutate(time = as.Date(paste(year, month, day, sep = "-"))) %>%
    dplyr::select(id, lat, lon, time, variable, value)


# CBEFS is slower, run on cluster "1_extract_CBEFS_SeineSurvey_sbass.R":
# cd RareEvents/
# sbatch --nodes=1 --mem=0 --time=1-23 R CMD BATCH "--vanilla --no-save --no-restore" code_fun/1_extract_CBEFS_SeineSurvey_sbass.R
# squeue -u lyubchich
ts_CBEFS <- data.table::fread("./dataderived/CBEFS_tsfor_SeineSurvey_sbass_CBEFS.csv",
                              na.strings = c("NA", "NaN", "N/A", ""))
# table(ts_CBEFS$variable)


# Rare events ----

# Apply rare_residuals() to each time series grouped by id
set.seed(123)
res_SeineSurvey_sbass <- ts_SeineSurvey_sbass[, .(data = list(rare_residuals(.SD, iforest_args = list(threshold = 0.7)))),
                                              by = .(id, variable)]
# set.seed(123)
# res_CBEFS <- ts_CBEFS[, .(data = list(rare_residuals(.SD,
#                                                      seasonal = TRUE, period = 365,
#                                                      iforest_args = list(threshold = 0.7)))),
#                       by = .(id, variable)]
# saveRDS(res_CBEFS, file = "dataderived/tmp_res_CBEFS.rds")
res_CBEFS <- readRDS("dataderived/tmp_res_CBEFS.rds")

# Plot some examples fish
ex = 1
ex = 2
ex = which(res_SeineSurvey_sbass$id == 106)
p1 <- res_SeineSurvey_sbass$data[[ex]]$data %>%
    dplyr::mutate(Anomaly = factor(direction,
                                   levels = c(1, 0, -1),
                                   labels = c("Positive", "0", "Negative"))) %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    labs(x = "Year", y = "CPUE") +
    geom_point(data = . %>% filter(is_anomaly_iforest),
               aes(x = time, y = value, color = Anomaly)) +
    ggtitle(paste(res_SeineSurvey_sbass$variable[[ex]], "at", res_SeineSurvey_sbass$id[[ex]], "iForest"))
p2 <- res_SeineSurvey_sbass$data[[ex]]$data %>%
    dplyr::mutate(Anomaly = factor(direction,
                                   levels = c(1, 0, -1),
                                   labels = c("Positive", "0", "Negative"))) %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    labs(x = "Year", y = "CPUE") +
    geom_point(data = . %>% filter(is_anomaly_dbscan),
               aes(x = time, y = value, color = Anomaly)) +
    ggtitle(paste(res_SeineSurvey_sbass$variable[[ex]], "at", res_SeineSurvey_sbass$id[[ex]], "DBSCAN"))
p1 + p2 +
    plot_layout(ncol = 2, byrow = TRUE, axis_titles = "collect", guides = "collect")

ggsave("images/ts_SeineSurvey_sbass_106.png", width = 9, height = 3)


# Plot some examples CBEFS
ex = 39
ex = 40
ex_sbass_match = match_SeineSurvey_sbass_CBEFS %>%
    dplyr::filter(id_SeineSurvey_sbass_stations == 106) %>%
    dplyr::pull(id_CBEFS_tmp)
ex_list <- which(res_CBEFS$id == ex_sbass_match)

for (ex in ex_list) {
    p1 <- res_CBEFS$data[[ex]]$data %>%
        dplyr::mutate(Anomaly = factor(sign(res_CBEFS$data[[ex]]$data$residual),
                                       levels = c(1, 0, -1),
                                       labels = c("Positive", "0", "Negative"))) %>%
        ggplot(aes(x = time, y = value)) +
        geom_line() +
        labs(x = "Time", y = res_CBEFS$variable[[ex]]) +
        geom_point(data = . %>% filter(is_anomaly_iforest),
                   aes(x = time, y = value, color = Anomaly)) +
        ggtitle(paste(#res_CBEFS$variable[[ex]], "at",
            res_CBEFS$id[[ex]], "iForest"))
    p2 <- res_CBEFS$data[[ex]]$data %>%
        dplyr::mutate(Anomaly = factor(sign(res_CBEFS$data[[ex]]$data$residual),
                                       levels = c(1, 0, -1),
                                       labels = c("Positive", "0", "Negative"))) %>%
        ggplot(aes(x = time, y = value)) +
        geom_line() +
        labs(x = "Time", y = res_CBEFS$variable[[ex]]) +
        geom_point(data = . %>% filter(is_anomaly_dbscan),
                   aes(x = time, y = value, color = Anomaly)) +
        ggtitle(paste(#res_CBEFS$variable[[ex]], "at",
            res_CBEFS$id[[ex]], "DBSCAN")) +
        theme(plot.title = element_text(size = 12))
    p1 + p2 +
        plot_layout(ncol = 2, byrow = TRUE, axis_titles = "collect", guides = "collect")

    ggsave(paste0("images/ts_CBEFS_", res_CBEFS$variable[[ex]], "_for_SeineSurvey_sbass_106.png"),
           width = 9, height = 3)
}

# Cross rare events ----

# Count rare events in fisheries
rare_SeineSurvey_sbass <- tibble()
for (i in 1:nrow(res_SeineSurvey_sbass)) {
    # Data for the given combination of location id and variable
    d <- res_SeineSurvey_sbass$data[[i]]$data %>%
        dplyr::mutate(month = lubridate::month(time)) %>%
        dplyr::filter(month > 6) %>%
        dplyr::group_by(year) %>%
        dplyr::summarise(is_anomaly_iforest = sum(is_anomaly_iforest),
                         is_anomaly_dbscan = sum(is_anomaly_dbscan),
                         direction = mean(sign(residual)) %>% sign(),
                         .groups = "drop") %>%
        dplyr::mutate(id = res_SeineSurvey_sbass$id[i],
                      variable = res_SeineSurvey_sbass$variable[i])

    # Append data
    rare_SeineSurvey_sbass <- rbind(rare_SeineSurvey_sbass, d)
}


# Count rare events in env. data in April-June
rare_CBEFS <- tibble()
for (i in 1:nrow(res_CBEFS)) {
    # Data for the given combination of location id and variable
    d <- res_CBEFS$data[[i]]$data %>%
        dplyr::mutate(month = lubridate::month(time)) %>%
        dplyr::filter(month %in% c(4, 5, 6)) %>%
        dplyr::group_by(year, month) %>%
        dplyr::summarise(is_anomaly_iforest = sum(is_anomaly_iforest),
                         is_anomaly_dbscan = sum(is_anomaly_dbscan),
                         direction = mean(sign(residual)) %>% sign(),
                         .groups = "drop") %>%
        dplyr::mutate(id = res_CBEFS$id[i],
                      variable = res_CBEFS$variable[i])

    # Append data
    rare_CBEFS <- rbind(rare_CBEFS, d)
}

# Aggregate months by year
rare_CBEFS_agg <- rare_CBEFS %>%
    group_by(id, variable, year) %>%
    dplyr::summarise(is_anomaly_iforest = sum(is_anomaly_iforest),
                     is_anomaly_dbscan = sum(is_anomaly_dbscan),
                     direction = mean(sign(direction)) %>% sign(),
                     .groups = "drop")

# Percent of years labeled as rare by each method
(rare_CBEFS_agg$is_anomaly_iforest > 0) %>% mean() * 100
# 2.546296
(rare_CBEFS_agg$is_anomaly_dbscan > 0) %>% mean() * 100
# 34.06339

# Wide format for variables in columns
rare_CBEFS_agg_wide_positive <- rare_CBEFS_agg %>%
    # Analyze positive env. anomalies, assign negative to 0 to disregard them
    dplyr::mutate(is_anomaly_iforest = ifelse(direction > 0, is_anomaly_iforest, 0),
                  is_anomaly_dbscan = ifelse(direction > 0, is_anomaly_dbscan, 0)) %>%
    tidyr::pivot_wider(id_cols = c(id, year),
                       names_from = variable,
                       values_from = c(is_anomaly_iforest, is_anomaly_dbscan),
                       names_glue = "{.value}_{variable}")

rare_CBEFS_agg_wide_negative <- rare_CBEFS_agg %>%
    # Analyze negative env. anomalies, assign positive to 0 to disregard them
    dplyr::mutate(is_anomaly_iforest = ifelse(direction > 0, 0, is_anomaly_iforest),
                  is_anomaly_dbscan = ifelse(direction > 0, 0, is_anomaly_dbscan)) %>%
    tidyr::pivot_wider(id_cols = c(id, year),
                       names_from = variable,
                       values_from = c(is_anomaly_iforest, is_anomaly_dbscan),
                       names_glue = "{.value}_{variable}")


# Combine fish and env data
D_positive <- dplyr::right_join(match_SeineSurvey_sbass_CBEFS,
                      rare_SeineSurvey_sbass %>% dplyr::select(-variable),
                      by = c("id_SeineSurvey_sbass_stations" = "id")) %>%
    dplyr::left_join(rare_CBEFS_agg_wide_positive %>% dplyr::rename(id_CBEFS_tmp = id),
                     by = c("id_CBEFS_tmp", "year")) %>%
    tibble()
D_split_positive <- D_positive %>%
    group_by(id_SeineSurvey_sbass_stations) %>%
    group_split()

D_negative <- dplyr::right_join(match_SeineSurvey_sbass_CBEFS,
                                rare_SeineSurvey_sbass %>% dplyr::select(-variable),
                                by = c("id_SeineSurvey_sbass_stations" = "id")) %>%
    dplyr::left_join(rare_CBEFS_agg_wide_negative %>% dplyr::rename(id_CBEFS_tmp = id),
                     by = c("id_CBEFS_tmp", "year")) %>%
    tibble()
D_split_negative <- D_negative %>%
    group_by(id_SeineSurvey_sbass_stations) %>%
    group_split()

# Identify columns for iforest and dbscan
iforest_cols <- colnames(D_positive)[grepl("iforest", colnames(D_positive)) & colnames(D_positive) != "is_anomaly_iforest"]
dbscan_cols <- colnames(D_positive)[grepl("dbscan", colnames(D_positive)) & colnames(D_positive) != "is_anomaly_dbscan"]

# Compute contingency tables and chi-square tests for iforest
iforest_results_positive <- map(D_split_positive, function(group) {
    id <- unique(group$id_SeineSurvey_sbass_stations)
    map(iforest_cols, function(col) {
        as_tibble(contingency_chisq(
            x = group$is_anomaly_iforest > 0,
            y = group[[col]] > 0,
            id = id
        )[c("id", "chisq_statistic", "p_value")]) %>%
            mutate(comparison = paste("is_anomaly_iforest vs", col))
    }) %>%
        bind_rows()
}) %>%
    bind_rows()

# Compute contingency tables and chi-square tests for dbscan
dbscan_results <- map(D_split, function(group) {
    id <- unique(group$id_SeineSurvey_sbass_stations)
    map(dbscan_cols, function(col) {
        as_tibble(contingency_chisq(
            x = group$is_anomaly_dbscan > 0,
            y = group[[col]] > 0,
            id = id
        )[c("id", "chisq_statistic", "p_value")]) %>%
            mutate(comparison = paste("is_anomaly_dbscan vs", col))
    }) %>%
        bind_rows()
}) %>%
    bind_rows()

# Combine results
final_results <- bind_rows(iforest_results, dbscan_results) %>%
    select(id, comparison, chisq_statistic, p_value)

# Print or inspect results
print(final_results)

# Show results for the station
final_results %>%
    dplyr::filter(id == 106) %>%
    readr::write_csv("dataderived/chisq_SeineSurvey_sbass_CBEFS_106.csv")

# Summarize proportions of significant variables across methods
final_results %>%
    group_by(comparison) %>%
    summarise(signif_prop = mean(p_value < 0.05, na.rm = TRUE)) %>%
    readr::write_csv("dataderived/chisq_SeineSurvey_sbass_CBEFS.csv")


# Regression modeling ----

## Compile data for regression modeling ----

# Monthly maxima of CBEFS variables
ts_CBEFS_monthly_max <- ts_CBEFS[,
                                 .(value = max(value, na.rm = TRUE)),
                                 by = .(variable, id,
                                        year = as.numeric(format(time, "%Y")),
                                        month = as.numeric(format(time, "%m")))
]
ts_CBEFS_monthly_min <- ts_CBEFS[,
                                 .(value = min(value, na.rm = TRUE)),
                                 by = .(variable, id,
                                        year = as.numeric(format(time, "%Y")),
                                        month = as.numeric(format(time, "%m")))
]
ts_CBEFS_monthly_avg <- ts_CBEFS[,
                                 .(value = mean(value, na.rm = TRUE)),
                                 by = .(variable, id,
                                        year = as.numeric(format(time, "%Y")),
                                        month = as.numeric(format(time, "%m")))
]


# Put the data in wider format
ts_CBEFS_monthly_max_wide <- ts_CBEFS_monthly_max %>%
    dplyr::filter(month %in% c(4, 5, 6)) %>%
    tidyr::pivot_wider(id_cols = c(id, year),
                       names_from = c(variable, month),
                       values_from = value,
                       names_prefix = "CBEFS_max_")
ts_CBEFS_monthly_min_wide <- ts_CBEFS_monthly_min %>%
    dplyr::filter(month %in% c(4, 5, 6)) %>%
    tidyr::pivot_wider(id_cols = c(id, year),
                       names_from = c(variable, month),
                       values_from = value,
                       names_prefix = "CBEFS_min_")
ts_CBEFS_monthly_avg_wide <- ts_CBEFS_monthly_avg %>%
    dplyr::filter(month %in% c(4, 5, 6)) %>%
    tidyr::pivot_wider(id_cols = c(id, year),
                       names_from = c(variable, month),
                       values_from = value,
                       names_prefix = "CBEFS_avg_")

# Merge fisheries and environmental data

D_regression <- SeineSurvey_sbass %>%
    dplyr::select(-variable, -SPP_AGE, -CATCH, -STATIONTYPE) %>%
    # Average fish data within certain months
    dplyr::filter((6 < month) & (month < 11)) %>%
    group_by(year, month, id) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::select(-month) %>%
    dplyr::left_join(match_SeineSurvey_sbass_CBEFS %>%
                        dplyr::select(id_SeineSurvey_sbass_stations, id_CBEFS_tmp),
                      by = c("id" = "id_SeineSurvey_sbass_stations")) %>%
    dplyr::left_join(ts_CBEFS_monthly_max_wide,
                     by = c("id_CBEFS_tmp" = "id", "year")) %>%
    dplyr::left_join(ts_CBEFS_monthly_min_wide,
                     by = c("id_CBEFS_tmp" = "id", "year")) %>%
    dplyr::left_join(ts_CBEFS_monthly_avg_wide,
                     by = c("id_CBEFS_tmp" = "id", "year"))  %>%
    dplyr::select(-id_CBEFS_tmp)


## Cross-validation ----

library(ranger)
source("code_fun/fun_gof_qr.R")


# Create a cross-validation fold column
K = 10
set.seed(123)
D_regression <- D_regression %>%
    dplyr::mutate(fold = sample(1:K, n(), replace = TRUE)) %>%
    # Add lagged values for each id
    dplyr::group_by(id) %>%
    dplyr::mutate(value_lag1 = dplyr::lag(value, n = 1, order_by = year)) %>%
    dplyr::ungroup()

# Save data for regression modeling
readr::write_csv(D_regression, "dataderived/D_regression.csv")

# Check the distribution of folds
table(D_regression$fold)

# Define quantiles for quantile regression
qs <- c(0.05, 0.5, 0.95)

results <- tibble()
for (k in 1:K) { # k = 1
    # Split data into training and testing sets
    train_data <- D_regression %>%
        # dplyr::select(-id, year) %>%
        dplyr::filter(fold != k) %>%
        dplyr::select(-fold) %>%
        dplyr::arrange(id, year)
    test_data <- D_regression %>%
        # dplyr::select(-id, year) %>%
        dplyr::filter(fold == k) %>%
        dplyr::select(-fold) %>%
        dplyr::arrange(id, year)
    test_y <- test_data %>% pull(value)

    # PCA
    train_pca <- prcomp(train_data %>% dplyr::select(-any_of(c("value", "value_lag1", "id", "year"))),
                        center = TRUE, scale. = TRUE)
    test_pca_scores <- predict(train_pca,
                        newdata = test_data %>% dplyr::select(-any_of(c("value", "value_lag1", "id", "year"))))
    train_pca_scores <- as_tibble(predict(train_pca))
    # Proportion of variance explained
    cum_var <- summary(train_pca)$importance[2, ] %>% cumsum
    # Select PCs explaining at least 80% variance
    n_pcs <- min(which(cum_var >= 0.8))
    # # Subset the first n_pcs components
    train_pca_scores_subset <- tibble(train_pca_scores[, 1:n_pcs, drop = FALSE],
                                      value = train_data$value,
                                      value_lag1 = train_data$value_lag1,
                                      id = train_data$id,
                                      year = train_data$year)
    test_pca_scores_subset <- test_pca_scores[, 1:n_pcs, drop = FALSE] %>%
        as_tibble() %>%
        mutate(id = test_data$id,
               year = test_data$year)

    # Fit the model on the training set
    model_qrf <- ranger(value ~ .,
                        data = train_data,
                        num.trees = 1000,
                        min.node.size = 5,
                        keep.inbag = TRUE,
                        quantreg = TRUE,
                        # importance = "impurity",
                        num.threads = 4
    )
    model_qr <- rq(value ~ .,
                   data = train_data,
                   method = "br",
                   tau = qs)
    model_qrf_pca <- ranger(value ~ .,
                        data = train_pca_scores_subset,
                        num.trees = 1000,
                        min.node.size = 5,
                        keep.inbag = TRUE,
                        quantreg = TRUE,
                        # importance = "impurity",
                        num.threads = 4
    )
    model_qr_pca <- rq(value ~ .,
                   data = train_pca_scores_subset,
                   method = "br",
                   tau = qs)

    # Make predictions on the test set
    predict_qrf <- predict(model_qrf, test_data, type = "quantiles", quantiles = qs)$predictions
    predict_qr <- predict(model_qr, newdata = test_data, type = "quantiles")
    predict_qrf_pca <- predict(model_qrf_pca, test_pca_scores_subset, type = "quantiles", quantiles = qs)$predictions
    predict_qr_pca <- predict(model_qr_pca, newdata = test_pca_scores_subset, type = "quantiles")

    # Calculate performance metrics
    perf_qrf <- gof_qr(test_y, predict_qrf)
    perf_qr <- gof_qr(test_y, predict_qr)
    perf_qrf_pca <- gof_qr(test_y, predict_qrf_pca)
    perf_qr_pca <- gof_qr(test_y, predict_qr_pca)

    results <- results %>%
        bind_rows(tibble(
            fold = k,
            quantile = qs,
            method = "QRF",
            R1 = perf_qrf["R1", ],
            ATWE = perf_qrf["ATWE", ]
        )) %>%
        bind_rows(tibble(
            fold = k,
            quantile = qs,
            method = "QR",
            R1 = perf_qr["R1", ],
            ATWE = perf_qr["ATWE", ]
        )) %>%
        bind_rows(tibble(
            fold = k,
            quantile = qs,
            method = "QRF_PCA",
            R1 = perf_qrf_pca["R1", ],
            ATWE = perf_qrf_pca["ATWE", ]
        )) %>%
        bind_rows(tibble(
            fold = k,
            quantile = qs,
            method = "QR_PCA",
            R1 = perf_qr_pca["R1", ],
            ATWE = perf_qr_pca["ATWE", ]
        ))
}
# results

# Save results
# readr::write_csv(results, "dataderived/CVresults_qrf_qr.csv")
results <- readr::read_csv("dataderived/CVresults_qrf_qr.csv")

# Boxplot results (facets by R1 or ATWE; plot boxplots for different quantiles from different
# methods together, grouped by quantile, so that we can see how the quantiles compare across methods)

results %>%
    dplyr::mutate(quantile = factor(quantile, levels = qs)) %>%
    pivot_longer(cols = c(R1, ATWE), names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = quantile, y = value, fill = method)) +
    geom_boxplot() +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = "Quantile", y = "Value", fill = "Method") +
    ggtitle("Cross-validation results") +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()
    )
ggsave("images/CVresults_qrf_qr.png", width = 10, height = 6, bg = "white")


library(gamlss)
?gamlss.family
set.seed(123)
model_gamlss1 <- gamlss::gamlss(
    value ~ .,
    data = train_pca_scores_subset %>% mutate(id = as.factor(id)),
    control = gamlss.control(trace = TRUE,
                             n.cyc = 100,
                             c.crit = 0.01),
    family = ZAGA #ZAGA #SEP2 #SEP1 # WEI GG
)
plot(model_gamlss1)


set.seed(123)
model_gamlss1 <- gamlss.inf::gamlssZadj(
    y = value,
    mu.formula =  ~ pb(PC1) + pb(PC2) + pb(PC3) + pb(PC4) + pb(PC5) + pb(PC6) + pb(PC7) + id,
    data = train_pca_scores_subset,
    control = gamlss.control(trace = TRUE,
                             n.cyc = 50,
                             c.crit = 0.01),
    family = GG #WEI2 GA GG LOGNO
)
plot(model_gamlss1)
AIC(model_gamlss1)
term.plot(model_gamlss1)


library(gamlss.inf)
set.seed(123)
model_gamlss1 <- gamlss.inf::gamlssZadj(
    y = value,
    mu.formula = ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + random(id) + pb(value_lag1),
    sigma.formula = ~id,
    data = train_pca_scores_subset %>% mutate(id = as.factor(id)) %>% na.omit(),
    control = gamlss.control(trace = TRUE,
                             n.cyc = 50,
                             c.crit = 0.01),
    family = GG #WEI2 GA GG LOGNO
)

# Incorporate random effects
model_gamlss1 <- gamlss.inf::gamlssZadj(
    y = value,
    mu.formula = ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + random(id),
    sigma.formula = ~id + PC1,
    correlation = corARMA(p = 2, q = 2, form = ~ year | id),
    data = train_pca_scores_subset %>% mutate(id = as.factor(id)),
    control = gamlss.control(trace = TRUE,
                             n.cyc = 50,
                             c.crit = 0.01),
    family = GG #WEI2 GA GG LOGNO
)

model_gamlss1 <- gamlss.inf::gamlssZadj(
    y = value,
    mu.formula = ~ re(random=~1|id),
    sigma.formula = ~id,
    correlation = corARMA(p = 2, q = 2, form = ~ year | id),
    data = train_pca_scores_subset %>% mutate(id = as.factor(id)),
    control = gamlss.control(trace = TRUE,
                             n.cyc = 50,
                             c.crit = 0.01),
    family = GG #WEI2 GA GG LOGNO
)



plot(model_gamlss1)
plot(model_gamlss1, ts = TRUE)
summary(getSmo(model_gamlss1$dist))

ranef(model_gamlss1)


summary(model_gamlss1)


str(model_gamlss1)
summary(model_gamlss1$dist)
model_gamlss1$dist$mu.coefSmo

model_gamlss1$dist$mu.coefficients


# library(gamlssx)
# model_gamlss1 <- gamlssx::fitGEV(
#     (value + 0.01)*1000 ~ .,
#     data = train_pca_scores_subset,
#     # mu.link = "log",
#     # scoring = "quasi",
#     control = gamlss.control(trace = TRUE,
#                              n.cyc = 100,
#                              c.crit = 0.01)
# )

# Summary of the selected model
summary(model_gamlss1)

# ANOVA table for forward selection steps
model_gamlss1$anova

# Plot diagnostics for the selected model
plot(model_gamlss1)



# ## Quantile random forest ----
#
# # library(grf)
# library(ranger)
# source("code_fun/fun_gof_qr.R")
#
# # Fit a quantile random forest model
# model_qrf <- ranger(value ~ .,
#     data = D_regression %>% dplyr::select(-id, -id_CBEFS_tmp, -month),
#     num.trees = 1000,
#     min.node.size = 5,
#     keep.inbag = TRUE,
#     quantreg = TRUE,
#     importance = "impurity",
#     num.threads = 4
# )
#
# # Print model summary
# print(model_qrf)
#
# Plot variable importance sorted
model_qrf <- ranger(value ~ .,
                    data = train_data,
                    num.trees = 1000,
                    min.node.size = 5,
                    keep.inbag = TRUE,
                    quantreg = TRUE,
                    importance = "impurity",
                    num.threads = 4
)
qrf_importance <- tibble(variable = names(model_qrf$variable.importance),
                          importance = model_qrf$variable.importance) %>%
    arrange(desc(importance))

ggplot(qrf_importance, aes(x = reorder(variable, importance), y = importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "Variable", y = "Importance") +
    ggtitle("Variable Importance in QRF Model")
#
#
# ## Quantile regression ----
# library(quantreg)
# # Fit a quantile regression model
# model_qr <- rq(value ~ .,
#     data = D_regression %>% dplyr::select(-id, -id_CBEFS_tmp, -month),
#     tau = c(0.05, 0.5, 0.95),
#     method = "fn",
#     na.action = na.omit
# )
# # Print model summary
# summary(model_qr)
#
# # Plot the coefficients for different quantiles
# coef_df <- as.data.frame(coef(model_qr)) %>%
#     tibble::rownames_to_column(var = "variable") %>%
#     pivot_longer(-variable, names_to = "quantile", values_to = "coefficient") %>%
#     mutate(quantile = factor(gsub("tau= ", "", quantile))) %>%
#     dplyr::filter(variable != "(Intercept)")
# ggplot(coef_df, aes(x = variable, y = coefficient, fill = quantile)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     coord_flip() +
#     labs(x = "Variable", y = "Coefficient") +
#     ggtitle("Quantile Regression Coefficients")
#
#
# ## GAM / GAMLSS / GEV ----
#
# # Historgram of the response variable
# ggplot(D_regression, aes(x = value)) +
#     geom_histogram(bins = 30, fill = "steelblue", color = "black", boundary = 0) +
#     labs(x = "Value", y = "Frequency") +
#     # log scale for better visibility
#     scale_x_log10() +
#     ggtitle("Histogram of the Response Variable (CPUE)")
#
# # Select predictors with monotonic relationships (significant Spearman correlation)
# predictors <- D_regression %>%
#     dplyr::select(-id, -id_CBEFS_tmp, -month, -value) %>%
#     dplyr::select_if(~ is.numeric(.)) %>%
#     colnames()
#
# # Calculate Spearman correlations
# SP <- lapply(predictors, function(p) {
#     tmp <- cor.test(D_regression$value + 0.01, D_regression[[p]], method = "spearman")
#     tibble(correlation = tmp$estimate, p_value = tmp$p.value, variable = p)
# }) %>%
#     bind_rows()
#
# # Filter significant correlations
# predictors_signif <- SP %>%
#     dplyr::filter(p_value < 0.05 & abs(correlation) > 0.05) %>%
#     dplyr::pull(variable)
#
# # PCA of predictors, for PCA regression
# library(ggfortify)
# pca <- prcomp(D_regression %>% dplyr::select(all_of(predictors)), scale. = TRUE)
# autoplot(pca, data = D_regression, colour = 'value') +
#     labs(title = "PCA of Predictors") +
#     theme(legend.position = "bottom")
#
# # Extract PC scores
# pca_scores <- as_tibble(predict(pca))
#
# # Summary of PCA
# summary_pca <- summary(pca)
# cum_var <- cumsum(summary_pca$importance[2, ])  # Proportion of variance explained
# n_pcs <- min(which(cum_var >= 0.8))  # Select PCs explaining at least 80% variance
# cat("Number of PCs explaining at least 80% variance:", n_pcs, "\n")
#
# # Subset the first n_pcs components
# pca_scores_subset <- pca_scores[, 1:n_pcs, drop = FALSE]
#
# # Combine PC scores with response variable
# D_regression_pca <- tibble(value = D_regression$value, pca_scores_subset)
#
#
#
# # Select a GAMLSS using AIC
#
library(gamlss)
#
# # Create a formula for the model, wrap predictors with pb()
# formula_gamlss <- as.formula(paste("value + 0.01 ~",
#                                    paste(paste0("pb(", names(pca_scores_subset), ")"), collapse = " + ")))
# scope <- as.formula(paste("~", paste(paste0("pb(", names(pca_scores_subset), ")"), collapse = " + ")))
#
# save.image("dataderived/tmp_image_gamlss.RData")
#
# # Fit an initial minimal model (intercept-only or no predictors)
# model_gamlss0 <- gamlss(
#     formula = value + 0.01 ~ 1,  # Minimal model with intercept only
#     data = D_regression_pca,
#     family = LOGNO2,
#     control = gamlss.control(trace = FALSE)
# )
#
# # Perform forward selection with stepGAIC
# model_gamlss1 <- gamlss::stepGAIC(
#     model_gamlss0,
#     scope = scope,
#     direction = "forward",
#     parallel = "multicore",
#     ncpus = 8
# )
# # Selected model:  value + 0.01 ~ pb(PC3) + pb(PC6) + pb(PC5) + pb(PC7) + pb(PC1) + pb(PC2)
#
model_gamlss1 <- gamlss::gamlss(
    value + 0.01 ~  pb(PC3) + pb(PC6) + pb(PC5) + pb(PC7) + pb(PC1) + id,
    data = train_pca_scores_subset,
    family = GG #SEP2 #SEP1 # WEI
)


model_gamlss1 <- gamlss::gamlss(
    value + 0.01 ~  pb(PC1) + pb(PC2)  + id,
    data = train_pca_scores_subset,
    family = GG #SEP2 #SEP1 # WEI
)


# Summary of the selected model
summary(model_gamlss1)

# ANOVA table for forward selection steps
model_gamlss1$anova

# Plot diagnostics for the selected model
plot(model_gamlss1)


#
# library(gamlssx)
# model_gamlssGEV1 <- gamlssx::fitGEV(
#     value + 0.01 ~ pb(PC3, max.df = 5) + pb(PC6, max.df = 5) + pb(PC5, max.df = 5) +
#         pb(PC7, max.df = 5) + pb(PC1, max.df = 5) + pb(PC2, max.df = 5),
#     data = train_pca_scores_subset,
#     control = gamlss.control(trace = TRUE,
#                              n.cyc = 100,
#                              c.crit = 0.1)
# )
#
# # Summary of the selected model
# summary(model_gamlssGEV1)
#
# # ANOVA table for forward selection steps
# model_gamlssGEV1$anova
#
# # Plot diagnostics for the selected model
# plot(model_gamlssGEV1)









