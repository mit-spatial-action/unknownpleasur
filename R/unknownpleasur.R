get_dims <- function(sdf, n, type = "horizontal") {
  #' Generate list of dimension-related quantities.
  #'
  #' @param bbox Object of type `bbox`.
  #' @param n Number of regularly-spaced lines to generate.
  #' @param type String indicating direction of lines---options are "horizontal" and "vertical".
  #' @returns A list of interval-relevant values.
  if (type == "horizontal") {
    axes <- list("int_dim" = "y", "edge_dim" = "x")
  } else if (type == "vertical") {
    axes <- list("int_dim" = "x", "edge_dim" = "y")
  } else {
    stop("Invalid line type - expects horizontal or vertical.")
  }
  bbox <- sf::st_bbox(sdf)
  int_min <- unname(bbox[paste0(axes$int_dim, "min")])
  int_max <- unname(bbox[paste0(axes$int_dim, "max")])
  edge_min <- unname(bbox[paste0(axes$edge_dim, "min")])
  edge_max <- unname(bbox[paste0(axes$edge_dim, "max")])
  list(
    "axes" = axes,
    "int_min" = int_min, 
    "int_max" = int_max, 
    "edge_min" = edge_min, 
    "edge_max" = edge_max, 
    "interval" = (int_max - int_min) / (n - 1),
    "n" = n,
    "type" = type
    )
}

st_regular_lines <- function(df, dims, mask = TRUE) {
  #' Generate n regularly-spaced lines over `df` extent.
  #'
  #' @param df An `sf` dataframe.
  #' @param dims Object returned by `get_dims()`.
  #' @param mask Whether resulting lines should be clipped to df extent.
  #' @returns A sf dataframe of regularly spaced POLYLINEs.
  #' @export
  message("Generating regularly spaced lines over input bbox...")
  bbox <- df %>%
    sf::st_bbox()
  line_positions <- seq(
      from = dims$int_min, 
      to = dims$int_max, 
      by = dims$interval
    )
  lines <- sf::st_sfc(crs = sf::st_crs(df))
  for (line in line_positions) {
    if (dims$type == "vertical") {
      coords <- c(line, dims$edge_max, line, dims$edge_min)
    } else if (dims$type == "horizontal") {
      coords <- c(dims$edge_min, line, dims$edge_max, line)
    } else {
      stop("Invalid line type - expects horizontal or vertical.")
    }
      
    line_const <- sf::st_sfc(
      sf::st_linestring(
        matrix(
          coords, 
          ncol = 2, 
          byrow = TRUE
          )
      ),
      crs = sf::st_crs(df)
    )
    lines <- sf::st_sfc(
        c(lines, line_const),
        crs = sf::st_crs(df)
    )
  }
  sf <- sf::st_as_sf(data.frame(id = 1:dims$n, geometry = lines)) 
  if (mask) {
    message("Clipping lines to input layer extent...")
    sf %>%
      sf::st_intersection(
        sf::st_geometry(df)
      ) %>%
      dplyr::rowwise() %>%
      dplyr::filter(st_geometry_type(geometry) != "POINT") %>%
      dplyr::ungroup() %>% 
      sf::st_cast("MULTILINESTRING") %>%
      sf::st_cast("LINESTRING", warn = FALSE) %>%
      dplyr::mutate(
        id = dplyr::row_number()
      )
  } else {
    sf
  }
}

st_unknown_pleasures <- function(
    lines, 
    raster, 
    dims,
    max = FALSE,
    sample_size = 250,
    bleed_factor = 1.5,
    mode = "planar",
    polygon = TRUE) {
  #' Transforms regularly-spaced lines into an "Unknown Pleasures"-esque set of regular section cuts based on raster value.
  #'
  #' @param lines A spatial dataframe containing regularly spaced lines.
  #' @param dims An object returned by the `get_dims()` function.
  #' @param sample_size Interval, in meters, to sample along lines.
  #' @param bleed_factor How much should maxima bleed into the area of the line above or below.
  #' @param mode If `planar`, results will be planar offset lines. If `xyz`, lines will be offset on `LINESTRING` z axis.
  #' @param polygon If `TRUE`, outputs polygons (or closed linestrings if paired with `mode = "xyz"`. If `FALSE`, outputs lines with no baseline.
  #' @returns A spatial dataframe containing "Unknown Pleasures" features.
  #' @export
  message("Sampling along lines...")
  elevated_lines <- lines %>%
    dplyr::mutate(
      geometry = sf::st_line_sample(
        geometry, 
        density = (1 / units::as_units(sample_size, "m")),
        type = "regular"
      )
    ) %>%
    sf::st_cast("POINT", warn = FALSE)
    
  message("Extracting elevations at sample points...")
  elevated_lines$elev <- raster::extract(raster, elevated_lines)
  elevated_lines
  if (!max) {
    max <- max(abs(elevated_lines$elev), na.rm = TRUE)
  }
  scale <- dims$interval / max
  elevated_lines <- elevated_lines %>%
    tidyr::drop_na(elev) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      elev_scaled = (elev + (0.01 * elev)) * (scale * bleed_factor)
    ) %>%
    tidyr::drop_na(elev_scaled)
  if (mode == "planar") {
    message("Performing Affine Transform on points...")
    elevated_lines <- elevated_lines %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        geometry = ifelse(
          (dims$type == "horizontal"),
          geometry + c(0, elev_scaled),
          geometry + c(elev_scaled, 0)
          ),
        coords = ifelse(
          (dims$type == "horizontal"),
          sf::st_coordinates(geometry)[,1],
          sf::st_coordinates(geometry)[,2]
        )
      ) %>%
      dplyr::ungroup() %>%
      sf::st_set_crs(sf::st_crs(lines)) %>%
      dplyr::group_by(id)
      if (dims$type == "vertical") {
        elevated_lines <- elevated_lines %>%
          dplyr::arrange(coords, .by_group = TRUE) %>%
          dplyr::ungroup()
      } else if (dims$type == "horizontal") {
        elevated_lines <- elevated_lines %>%
          dplyr::arrange(dplyr::desc(coords), .by_group = TRUE) %>%
          dplyr::ungroup()
      }
  } else if (mode == "xyz") {
    message("Attaching Z values to component points...")
    elevated_lines <- elevated_lines %>%
      dplyr::mutate(
        x = sf::st_coordinates(.)[,1],
        y = sf::st_coordinates(.)[,2]
      ) %>%
      sf::st_drop_geometry() %>%
      sf::st_as_sf(
        coords = c("x", "y", "elev_scaled"),
        dim = "XYZ",
        crs = sf::st_crs(lines)
      )
  }
  elevated_lines <- elevated_lines %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(do_union = FALSE) %>%
    sf::st_cast("LINESTRING")
  
  if (polygon) {
    message("Building closed loops...")
    lines <- lines %>%
      dplyr::filter(
        id %in% dplyr::pull(elevated_lines, id)
      )
    if (mode == "xyz") {
      lines <- lines %>%
        sf::st_zm(
          drop = FALSE, 
          what = "Z"
        ) %>%
        sf::st_reverse()
    }
    elevated_lines <- elevated_lines %>%
      dplyr::bind_rows(lines) %>%
      sf::st_cast("POINT", warn = FALSE) %>%
      dplyr::group_by(id) %>%
      dplyr::summarize(do_union = FALSE) %>%
      sf::st_cast("POLYGON") %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(id))
    if (mode == "xyz") {
      warning("Can't build polygons in XYZ mode---returning POLYLINES instead.")
      elevated_lines <- elevated_lines %>%
        sf::st_cast("POINT", warn = FALSE) %>%
        dplyr::mutate(
          z = sf::st_coordinates(.)[,3]
        ) %>%
        sf::st_zm(drop = TRUE) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          geometry = ifelse(
            (dims$type == "horizontal"),
            geometry + c(0, z),
            geometry + c(z, 0)
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(do_union = FALSE) %>%
        sf::st_cast("POLYGON") %>%
        dplyr::ungroup() %>%
        sf::st_buffer(0.0) %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON", warn = FALSE) %>%
        st_buffer(0.0) %>%
        sf::st_make_valid() %>%
        dplyr::mutate(
          id = dplyr::row_number()
        ) %>%
        sf::st_cast("POINT", warn = FALSE) %>%
        dplyr::mutate(
          x = sf::st_coordinates(.)[,1],
          y = sf::st_coordinates(.)[,2]
        ) %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(id) 
      if (dims$type == "horizontal") {
        elevated_lines <- elevated_lines %>%
          dplyr::mutate(
            z = y - min(y),
            y = y - z
          )
      } else {
        elevated_lines <- elevated_lines %>%
          dplyr::mutate(
            z = x - min(x),
            x = x - z
          )
      }
      elevated_lines %>%
        sf::st_as_sf(
          coords = c("x", "y", "z"),
          dim = "XYZ",
          crs = sf::st_crs(lines)
        ) %>%
        dplyr::group_by(id) %>%
        dplyr::filter(n() >= 4) %>%
        dplyr::summarize(do_union = FALSE) %>%
        sf::st_cast("POLYGON") %>%
        sf::st_cast("LINESTRING", warn = FALSE)
    } else {
      message("Returning polygons.")
      elevated_lines %>%
        sf::st_buffer(0.0)
    }
    
  } else {
    elevated_lines
  }
}