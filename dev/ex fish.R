library(rarefun)
?rarefun::rare_residuals
vignette(package = "rarefun")

devtools::install_github("atsa-es/atsalibrary")
library(atsalibrary)
atsalibrary::sockeye

library(dplyr)

x <- sockeye |>
  dplyr::ungroup() |>
  dplyr::filter(region == "NAKNEK") |>
  dplyr::transmute(value = spawners, 
                   time = brood_year) |> 
  na.omit()
result <- rare_residuals(x)

plot(result$data$time, result$data$value, type = "l", main = "Time Series")
points(result$data$time[result$data$is_anomaly_iforest],
       result$data$value[result$data$is_anomaly_iforest], col = "red", pch = 19)


