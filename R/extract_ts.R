extract_ts <- function(x,
                       dataname = c("daymet", "CBEFS"),
                       vars2use = NULL) {
    # Extract name
    xname <- deparse1(substitute(x))

    # Match dataname to one of the described datasets
    dataname <- match.arg(dataname)

    # Check that the locations are not duplicated, to speedup data extractions
    x_nodup <- x %>%
        filter(!duplicated(id))

    # Initialize an empty data.table outside the loop to accumulate results
    D <- data.table()

    if (dataname == "daymet") {

    } else if (dataname == "CBEFS") {
        if (is.null(vars2use)) {
            vars2use <- c("atten_depth", "diss_o2", "ph", "salinity", "temperature")
        }

        files_CBEFS <- list.files("dataraw/CBEFS",
                                  pattern = "genny_nesslage\\d{4}\\.nc",
                                  full.names = TRUE) %>%
            sort()

        # Loop through each NetCDF file
        for (fl in files_CBEFS) {

            # fl = files_CBEFS[1]; v = vars2use[6]; i = 3
            print(Sys.time())
            print(fl)

            # Open the NetCDF file
            nc <- nc_open(fl)

            # Loop through each cell to extract
            for (i in 1:nrow(x_nodup)) {
                d <- extract_data_per_cell(nc,
                                           cell_row = x_nodup$coord_row[i],
                                           cell_col = x_nodup$coord_col[i],
                                           variable = vars2use)
                D <- rbind(D, d)
            }

            # Close the NetCDF file
            nc_close(nc)
        }

        # Fix dates in the combined file, see nc$dim$time
        D <- D %>%
            mutate(time = as.Date(time, origin = "1975-01-01"))

    }
    return(D)
}

# Help function for CBEFS data extraction when need to search non-empty data from
# neighboring grid cells
matrix_neighbors_excluding_center <- function(mat_nrow, mat_ncol, x, y, radius = 1) {
    # Get indices of neighboring cells within the specified radius
    x_range <- pmax(1, x - radius):pmin(mat_nrow, x + radius)
    y_range <- pmax(1, y - radius):pmin(mat_ncol, y + radius)

    # Exclude the current cell (x, y) and its immediate neighbors
    x_range0 <- pmax(1, x - (radius - 1)):pmin(mat_nrow, x + radius - 1)
    y_range0 <- pmax(1, y - (radius - 1)):pmin(mat_ncol, y + radius - 1)

    # Report matrix indices of the resulting neighborhood
    indices <- expand.grid(x = x_range, y = y_range) %>%
        as_tibble() %>%
        mutate(ii = paste0("x", x, "_y", y))
    indices0 <- expand.grid(x = x_range0, y = y_range0) %>%
        as_tibble() %>%
        mutate(ii = paste0("x", x, "_y", y))

    indices %>%
        filter(!is.element(ii, indices0$ii)) %>%
        select(x, y)
}

# Help function for CBEFS data extraction
extract_data_per_cell <- function(nc_file,
                                  cell_row = 1L,
                                  cell_col = 1L,
                                  variable = c("atten_depth", "diss_o2", "ph",
                                               "salinity", "temperature", "wave_height"))
{
    Dcell <- data.table()
    for (v in variable) {

        # Extracting wave_height (3-dimensional variable)
        if (v == "wave_height") { # 3-dimensional variable
            x <- ncvar_get(nc_file, v,
                           start = c(cell_row, cell_col, 1),
                           count = c(1, 1, -1))

            # Based on the email from Pierre St-Laurent <pst-laurent@vims.edu> on 2024-07-10,
            # "interpolate" by using non-zero information from nearby cells
            if (all(x == 0)) {

                # Radius to search around
                rr = 1

                # while all wave heights are 0
                while (all(x == 0) | all(is.na(x))) {

                    # select neighbors around the current cell
                    nei <- matrix_neighbors_excluding_center(mat_nrow = nc_file$dim$x$len,
                                                             mat_ncol = nc_file$dim$y$len,
                                                             cell_row, cell_col,
                                                             radius = rr)

                    # start checking the neighboring cells until find one with non-zeros
                    i_nei = 1
                    while ((all(x == 0) | all(is.na(x))) & (i_nei <= nrow(nei))) {
                        x <- ncvar_get(nc_file, v,
                                       start = c(nei$x[i_nei], nei$y[i_nei], 1),
                                       count = c(1, 1, -1))
                        i_nei <- i_nei + 1
                    }

                    # if still not found, increase radius
                    rr <- rr + 1
                }
                print(paste0("Found a replacement wave height for cell [", cell_row, ", ", cell_col,
                             "] within radius ", rr))
            }

            # Create data.table for this subset
            X <- data.table(variable = v,
                            id = paste0("CBEFS_", cell_row, "_", cell_col), # Stations_nodup$CBEFS_id[i],
                            time = nc_file$dim$time$vals,
                            # layer = 0L,
                            value = x)
        }

        # Extracting 3-dimensional variables (2D + time) with surface level only (1L = Surface)
        else if (v == "atten_depth") {

            # Extract for the surface level, invert to calculate "diffuse attenuation coefficient at the surface"
            x <- 1 / ncvar_get(nc_file, v,
                           start = c(cell_row, cell_col, 2, 1),
                           count = c(1, 1, 1, -1))

            # Create data.table for this subset
            X <- data.table(variable = "DiffuseAttenuation",
                            id = paste0("CBEFS_", cell_row, "_", cell_col), #Stations_nodup$CBEFS_id[i],
                            time = nc_file$dim$time$vals,
                            # layer = 1L,
                            value = x)
        }

        # Extracting 4-dimensional variables with surface and bottom levels (1L = Surface, 20L = Bottom)
        else {

            # Extract for both levels at once
            x <- ncvar_get(nc_file, v,
                           start = c(cell_row, cell_col, 1, 1),
                           count = c(1, 1, -1, -1))

            # Create data.table for this subset by combining the levels
            X <- data.table(variable = v,
                            id = paste0("CBEFS_", cell_row, "_", cell_col), #Stations_nodup$CBEFS_id[i],
                            time = c(nc_file$dim$time$vals, nc_file$dim$time$vals),
                            layer = rep(c(20L, 1L), each = nc_file$dim$time$len),
                            value = c(x[1,], x[2,])) %>%
                dplyr::mutate(variable = paste(variable, layer, sep = "_")) %>%
                dplyr::select(-layer)
        }

        # Bind rows to Dcell
        Dcell <- rbind(Dcell, X)

    }
    return(Dcell)
}
