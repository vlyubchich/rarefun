#' Match Spatial Points between Two Datasets
#'
#' Matches each point in the first dataset (`A`) to the nearest point in the second
#' dataset (`B`) based on spatial coordinates (latitude and longitude). The function
#' supports both fast Euclidean distance-based matching (using `RANN::nn2`) and
#' slower, more accurate geodesic distance-based matching (using `geosphere::distGeo`).
#' Geodesic distances are always reported in the output, regardless of the matching
#' method. An optional maximum distance threshold can be applied to filter matches.
#'
#' @param A A data frame containing the first set of points with columns `id`
#'   (character or numeric identifier), `lat` (latitude in degrees, -90 to 90), and
#'   `lon` (longitude in degrees, -180 to 180).
#' @param B A data frame containing the second set of points with columns `id`
#'   (character or numeric identifier), `lat` (latitude in degrees, -90 to 90), and
#'   `lon` (longitude in degrees, -180 to 180).
#' @param max_dist Numeric, the maximum geodesic distance (in meters) for a valid
#'   match. Points with distances exceeding this threshold are assigned `NA` in the
#'   output. Default is `Inf` (no threshold).
#' @param fast Logical, indicating whether to use fast Euclidean distance-based
#'   matching (`TRUE`, default) via `RANN::nn2` or slower geodesic distance-based
#'   matching (`FALSE`) via `geosphere::distGeo`. Regardless of the method, final
#'   distances are computed as geodesic distances.
#' @param ... Additional arguments passed to `RANN::nn2` (if `fast = TRUE`).
#'
#' @return A data frame with one row per point in `A`, containing the following columns:
#'   \itemize{
#'     \item \code{id_A}: Identifier of the point from `A`.
#'     \item \code{lat_A}: Latitude of the point from `A`.
#'     \item \code{lon_A}: Longitude of the point from `A`.
#'     \item \code{id_B}: Identifier of the nearest point from `B` (or `NA` if no match within `max_dist`).
#'     \item \code{lat_B}: Latitude of the nearest point from `B` (or `NA`).
#'     \item \code{lon_B}: Longitude of the nearest point from `B` (or `NA`).
#'     \item \code{distance_m}: Geodesic distance (in meters) to the nearest point in `B` (or `NA` if no match).
#'   }
#'   Column names are dynamically prefixed with the names of the input data frames
#'   (e.g., `id_Aa` if `A` is named `Aa`).
#'
#' @details
#' The function matches each point in `A` to the nearest point in `B` based on spatial
#' coordinates. If `fast = TRUE`, it uses Euclidean distances for matching (faster but
#' less accurate for large distances) via `RANN::nn2`. If `fast = FALSE`, it computes
#' geodesic distances for matching using `geosphere::distGeo`, which is more accurate
#' but slower. In both cases, the reported `distance_m` is the geodesic distance.
#' If `max_dist` is specified, matches exceeding this distance are replaced with `NA`.
#' The function handles empty inputs by returning an empty data frame with appropriate
#' column names.
#'
#' @examples
#' \dontrun{
#' # Example datasets
#' Aa <- data.frame(
#'   id = c("fish1", "fish2", "fish3"),
#'   lat = c(40.7128, 40.7228, 40.7328),
#'   lon = c(-74.0060, -74.0160, -74.0260)
#' )
#' Bb <- data.frame(
#'   id = c("grid1", "grid2"),
#'   lat = c(40.7000, 40.7500),
#'   lon = c(-74.0000, -74.0200)
#' )
#'
#' # Fast matching (Euclidean-based)
#' match_spatial_points(Aa, Bb, fast = TRUE)
#'
#' # Accurate matching (geodesic-based)
#' match_spatial_points(Aa, Bb, fast = FALSE)
#'
#' # Apply maximum distance threshold
#' match_spatial_points(Aa, Bb, max_dist = 2000)
#' }
#'
#' @importFrom geosphere distGeo
#' @importFrom RANN nn2
#' @importFrom stats setNames
#' @export
#'
match_spatial_points <- function(A, B, max_dist = Inf, fast = TRUE, ...) {

  # Capture input data frame names
  name_A <- deparse(substitute(A))
  name_B <- deparse(substitute(B))

  # Input validation
  stopifnot(
    is.data.frame(A) && is.data.frame(B),
    all(c("id", "lat", "lon") %in% names(A)),
    all(c("id", "lat", "lon") %in% names(B)),
    is.numeric(A$lat) && is.numeric(A$lon),
    is.numeric(B$lat) && is.numeric(B$lon),
    all(is.finite(A$lat)) && all(is.finite(A$lon)),
    all(is.finite(B$lat)) && all(is.finite(B$lon)),
    all(A$lat >= -90 & A$lat <= 90),
    all(B$lat >= -90 & B$lat <= 90),
    all(A$lon >= -180 & A$lon <= 180),
    all(B$lon >= -180 & B$lon <= 180),
    is.numeric(max_dist) && length(max_dist) == 1 && max_dist >= 0
  )

  # Handle empty inputs
  if (nrow(A) == 0 || nrow(B) == 0) {
    return(data.frame(
      setNames(
        data.frame(
          character(0), numeric(0), numeric(0),
          character(0), numeric(0), numeric(0), numeric(0),
          stringsAsFactors = FALSE
        ),
        c(
          paste0("id_", name_A), paste0("lat_", name_A), paste0("lon_", name_A),
          paste0("id_", name_B), paste0("lat_", name_B), paste0("lon_", name_B),
          "distance_m"
        )
      )
    ))
  }

  if (fast) {
    # Find the nearest neighbor using RANN (based on Euclidean distance, not geodesic)
    # k = 1 for the single nearest neighbor
    nn_results <- RANN::nn2(query = A[, c("lon", "lat")],
                            data =  B[, c("lon", "lat")],
                            k = 1,
                            ...)

    # Get the indices of the nearest neighbors in B
    nearest_b_idx <- nn_results$nn.idx[, 1]

    # Populate B's information into the match_df
    match_df <- data.frame(
      id_A = A$id,
      lat_A = A$lat,
      lon_A = A$lon,
      id_B = B$id[nearest_b_idx],
      lat_B = B$lat[nearest_b_idx],
      lon_B = B$lon[nearest_b_idx],
      distance_m = geosphere::distGeo(
        p1 = A[, c("lon", "lat")],
        p2 = B[nearest_b_idx, c("lon", "lat")]
      ),
      stringsAsFactors = FALSE
    )

  } else {
    # Use geodesic distances (more accurate but slow)

    # Initialize output data frame with A's information
    n_A <- nrow(A)
    match_df <- data.frame(
      id_A = A$id,
      lat_A = A$lat,
      lon_A = A$lon,
      id_B = NA_character_,
      lat_B = NA_real_,
      lon_B = NA_real_,
      distance_m = NA_real_,
      stringsAsFactors = FALSE
    )

    # For each element in A, calculate distances and find nearest points
    for (i in 1:n_A) {

      # Coordinates of the current point from A
      point_A <- c(A$lon[i], A$lat[i])

      # Calculate geodesic distances to all points in B
      distances <- geosphere::distGeo(
        point_A,
        B[, c("lon", "lat")]
      )

      # Find index of minimum distance
      min_idx <- which.min(distances)

      # Fill in matching information
      match_df$id_B[i] <- B$id[min_idx]
      match_df$lat_B[i] <- B$lat[min_idx]
      match_df$lon_B[i] <- B$lon[min_idx]
      match_df$distance_m[i] <- distances[min_idx]
    }
  }

  # Rename columns to include input data frame names
  colnames(match_df) <- c(
    paste0("id_", name_A),
    paste0("lat_", name_A),
    paste0("lon_", name_A),
    paste0("id_", name_B),
    paste0("lat_", name_B),
    paste0("lon_", name_B),
    "distance_m"
  )

  # Apply maximum distance filter
  match_df[[paste0("id_", name_B)]][match_df$distance_m > max_dist] <- NA_character_
  match_df[[paste0("lat_", name_B)]][match_df$distance_m > max_dist] <- NA_real_
  match_df[[paste0("lon_", name_B)]][match_df$distance_m > max_dist] <- NA_real_
  match_df$distance_m[match_df$distance_m > max_dist] <- NA_real_

  return(match_df)
}
