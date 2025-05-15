#' Generate list of dimension-related quantities.
#'
#' @param df `sf` object.
#' @param n Number of regularly-spaced lines to generate.
#' @param type String indicating direction of lines---options are "horizontal" and "vertical".
#' @returns A list of interval-relevant values.
#' 
#' @export
up_get_dims <- function(df, n, type = "horizontal") {
  if (type == "horizontal") {
    axes <- list("int_dim" = "y", "edge_dim" = "x")
  } else if (type == "vertical") {
    axes <- list("int_dim" = "x", "edge_dim" = "y")
  } else {
    stop("Invalid line type - expects horizontal or vertical.")
  }
  bbox <- sf::st_bbox(df)
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

#' Interpolate polygons using inverse distance weighting
#'
#' @param x An `sf` object containing POLYGONs.
#' @param field Character. Name of field whose values should be interpolated.
#' @param nmax Maximum number of adjacent observations.
#' @param maxdist Maximum distance to be counted.
#' @param idp Inverse distance power. Defaults to 2.
#'
#' @returns A `SpatRaster`.
#' @export
up_interpolate_vector <- function(x, field, nmax = 15, maxdist = Inf, idp = 2) {
  pts <- x |>
    sf::st_point_on_surface() |>
    tidyr::drop_na(dplyr::all_of(field))
  
  raster <- tracts |>
    terra::ext() |>
    terra::rast(resolution = 250, crs = terra::crs(x))
  
  grid <- raster |>
    terra::as.points() |> 
    sf::st_as_sf()
  
  grid[[field]] <- gstat::gstat(
    formula = as.formula(glue::glue("{field} ~ 1")),
    data = pts,
    nmax = nmax,
    maxdist = maxdist,
    set = list(idp = idp)
  ) |>
    predict(newdata = grid) |>
    dplyr::pull(var1.pred)
  
  grid |>
    terra::vect() |>
    terra::rasterize(y = raster, field = field) |>
    terra::mask(terra::vect(x))
}

#' Mask lines by shapes.
#'
#' @param x An `sf` object to mask.
#' @param y An `sf` object to serve as mask/
#'
#' @returns An `sf` object.
#' @export
up_mask_lines <- function(x, y) {
  x |>
    sf::st_intersection(
      y |> 
        sf::st_union()
    ) |>
    dplyr::filter(sf::st_geometry_type(.data$geometry) != "POINT") |>
    sf::st_cast("MULTILINESTRING") |>
    sf::st_cast("LINESTRING", warn = FALSE) |>
    dplyr::mutate(
      id = dplyr::row_number()
    )
}

#' Generate regularly-spaced lines an an angle over an extent.
#'
#' @param x An `sf` dataframe.
#' @param cellsize A spacing distance in linear units.
#' @param angle Angle in degrees. 0 is horizontal.
#' @returns `sf` object of regularly spaced `LINESTRING`s.
#' 
#' @export
up_regular_lines <- function(x, 
                             cellsize,
                             angle = 0) {
  
  x |>
    up_rotate_extent(-angle) |>
    sf::st_make_grid(cellsize = cellsize, what="corners") |>
    sf::st_coordinates() |>
    tibble::as_tibble() |>
    sf::st_as_sf(
      coords=c("X", "Y"),
      crs = sf::st_crs(x),
      remove=FALSE
    ) |>
    dplyr::summarize(
      geometry = sf::st_cast(sf::st_union(geometry), "LINESTRING"),
      .by = Y
    ) |>
    dplyr::rename(
      id = Y
    ) |>
    up_rotate_extent(angle)
}

#' Drape Geometries Using Column Value
#'
#' @param x An `sf` object.
#' @param col Name of column containing offset values.
#'
#' @returns An `sf` object.
#' @export
up_drape_lines <- function(x, col) {
  x |>
    dplyr::rowwise() |>
    dplyr::mutate(
      geometry = list(sf::st_point(c(sf::st_coordinates(geometry)[1:2], .data[[col]])))
    ) |>
    dplyr::ungroup() |>
    sf::st_as_sf() |>
    dplyr::summarize(
      geometry = sf::st_cast(sf::st_union(geometry), "LINESTRING"),
      do_union = FALSE,
      .by = id
    ) |>
    sf::st_set_crs(sf::st_crs(x))
}

#' Offset Points based on a Column Value and an Angle
#'
#' @param x An `sf` object containing points.
#' @param angle Angle in degrees.
#' @param col Name of column containing offset values.
#'
#' @returns `sf` object containing offset points.
#' @export
up_offset_lines <- function(x, angle, col) {
  x |>
    up_rotate_extent(-angle) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      geometry = .data$geometry + as.numeric(c(0, .data[[col]]))
    ) |>
    dplyr::ungroup() |>
    dplyr::summarize(
      geometry = sf::st_cast(sf::st_union(geometry), "LINESTRING"),
      do_union = FALSE,
      .by = id
    ) |>
    up_rotate_extent(angle) |>
    sf::st_set_crs(sf::st_crs(x))
}

#' Construct a Rotation Matrix.
#'
#' @param a Angle in radians.
#'
#' @returns A rotation matrix.
up_rotate = function(a) {
  matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
}

#' Rotate a Given `sf` Object by an Arbitrary Angle.
#'
#' @param x An `sf` object.
#' @param angle Angle in degrees.
#'
#' @returns A rotated `sf` object.
#' @export
up_rotate_extent <- function(x, angle) {
  angle <- angle * pi / 180
  
  bbox <- sf::st_bbox(x)
  center <- c(mean(bbox[c("xmin", "xmax")]), mean(bbox[c("ymin", "ymax")]))
  
  x |>
    dplyr::mutate(
      geometry = (sf::st_geometry(x) - center) * up_rotate(angle) + center
    ) |>
    sf::st_set_crs(sf::st_crs(x))
}


#' Sample Along Lines
#'
#' @param x An `sf` object containing lines.
#' @param interval 
#'
#' @returns An `sf` object containing sampled points.
#' @export
up_sample_lines <- function(x, interval) {
  x |>
    dplyr::mutate(
      geometry = sf::st_line_sample(geometry, density = 1 / interval)
    ) |>
    dplyr::filter(lengths(geometry) > 0) |>
    sf::st_cast("POINT", warn = FALSE)
}

#' Elevate Lines Based on Raster Values
#'
#' @param x `sf` object containing `LINESTRING`s.
#' @param interval Line spacing. Can be a `units` object.
#' @param raster `terra` `Spatraster` from which to extract elevations.
#' @param angle Angle of lines in degrees.
#' @param scale 
#' @param factor Z exaggeration.
#' @param mode One of `"xyz"`, or `"planar"`. If `"xyz"` returns linestrings
#' with Z component. If `"planar"` returns lines offset by value.
#'
#' @returns `sf` object containing `LINESTRING`s.
#' @export
up_elevate <- function(
    x, 
    interval, 
    raster,
    angle = 90, 
    scale = "spacing",
    factor = 1,
    mode = "xyz") {
  
  x <- x |>
    up_sample_lines(interval) |>
    dplyr::mutate(
      z = terra::extract(raster, sf::st_as_sf(.data$geometry))[,2]
    )
  
  if (scale == "interval") {
    scaler <- (interval / max(abs(x$z), na.rm = TRUE))
  } else if (scale == "actual") {
    scaler <- 1
  }
  
  x <- x |>
    tidyr::drop_na(z) |>
    dplyr::mutate(
      z = .data$z * scaler * factor
    )
  
  if (mode == "xyz") {
    x <- x |>
      up_drape_lines(col = "z")
  }
  else if (mode == "planar") {
    x <- x |>
      up_offset_lines(angle = angle, col = "z")
  }
  x
}

#' Transforms regularly-spaced lines into an "Unknown Pleasures"-esque set of regular section cuts based on raster value.
#'
#' @param lines `sf` object containing regularly spaced lines.
#' @param raster `raster` object.
#' @param dims An object returned by the `up_get_dims()` function.
#' @param max Numeric.
#' @param interval Interval. Can be `units` object. Otherwise will be in unit of CRS.
#' @param bleed_factor How much should maxima bleed into the area of the line above or below.
#' @param mode If `planar`, results will be planar offset lines. If `xyz`, lines will be offset on `LINESTRING` z axis.
#' @param polygon If `TRUE`, outputs polygons (or closed linestrings if paired with `mode = "xyz"`. If `FALSE`, outputs lines with no baseline.
#' 
#' @returns `sf` object  containing "Unknown Pleasures" features.
#' 
#' @export
up_unknown_pleasures <- function(
    lines, 
    raster, 
    dims,
    interval,
    max = NULL,
    bleed_factor = 1.5,
    mode = "planar",
    polygon = TRUE) {
  
  line_points <- up_sample_lines(lines, interval)
    
  message("Extracting elevations at sample points...")
  line_points$elev <- terra::extract(raster, line_points)
  
  if (is.null(max)) {
    max <- max(abs(line_points$elev), na.rm = TRUE)
  }
  
  line_points <- line_points  |>
    tidyr::drop_na(.data$elev) |>
    # Filter out lines with only one point.
    dplyr::group_by(.data$id) |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::ungroup() |>
    # Scale elvation.
    dplyr::mutate(
      elev_scaled = (.data$elev + (0.01 * .data$elev)) * (dims$interval / max * bleed_factor)
    ) |>
    tidyr::drop_na(.data$elev_scaled)
  if (mode == "planar") {
    message("Performing Affine Transform on points...")
    line_points <- line_points |>
      dplyr::rowwise() |>
      dplyr::mutate(
        geometry = dplyr::case_when(
          dims$type == "horizontal" ~ .data$geometry + c(0, .data$elev_scaled),
          dims$type == "vertical" ~ .data$geometry + c(.data$elev_scaled, 0),
          .default = stop("Invalid type.")
        ),
        coords = dplyr::case_when(
          dims$type == "horizontal" ~ sf::st_coordinates(.data$geometry)[,1],
          dims$type == "vertical" ~ sf::st_coordinates(.data$geometry)[,2],
          .default = stop("Invalid type.")
        )
      ) |>
      dplyr::ungroup() |>
      sf::st_set_crs(sf::st_crs(lines)) |>
      dplyr::group_by(.data$id)
      if (dims$type == "vertical") {
        line_points <- line_points|>
          dplyr::arrange(.data$coords, .by_group = TRUE) |>
          dplyr::ungroup()
      } else if (dims$type == "horizontal") {
        line_points <- line_points |>
          dplyr::arrange(dplyr::desc(.data$coords), .by_group = TRUE) |>
          dplyr::ungroup()
      }
  } else if (mode == "xyz") {
    message("Attaching Z values to component points...")
    line_points <- line_points |>
      dplyr::mutate(
        x = sf::st_coordinates(.data$geometry)[,1],
        y = sf::st_coordinates(.data$geometry)[,2]
      ) |>
      sf::st_drop_geometry() |>
      sf::st_as_sf(
        coords = c("x", "y", "elev_scaled"),
        dim = "XYZ",
        crs = sf::st_crs(lines)
      )
  }
  line_points <- line_points |>
    dplyr::group_by(.data$id) |>
    dplyr::summarize(do_union = FALSE) |>
    sf::st_cast("LINESTRING")
  
  if (polygon) {
    message("Building closed loops...")
    lines <- lines |>
      dplyr::filter(
        .data$id %in% dplyr::pull(line_points, .data$id)
      )
    if (mode == "xyz") {
      lines <- lines |>
        sf::st_zm(
          drop = FALSE, 
          what = "Z"
        ) |>
        sf::st_reverse()
    }
    line_points <- line_points |>
      dplyr::bind_rows(lines) |>
      sf::st_cast("POINT", warn = FALSE) |>
      dplyr::group_by(.data$id) |>
      dplyr::summarize(do_union = FALSE) |>
      sf::st_cast("POLYGON") |>
      dplyr::ungroup() |>
      dplyr::arrange(dplyr::desc(.data$id))
    if (mode == "xyz") {
      warning("Can't build polygons in XYZ mode---returning POLYLINES instead.")
      line_points <- line_points |>
        sf::st_cast("POINT", warn = FALSE) |>
        dplyr::mutate(
          z = sf::st_coordinates(.data$geometry)[,3]
        ) |>
        sf::st_zm(drop = TRUE) |>
        dplyr::rowwise() |>
        dplyr::mutate(
          geometry = ifelse(
            (dims$type == "horizontal"),
            .data$geometry + c(0, .data$z),
            .data$geometry + c(.data$z, 0)
          )
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(.data$id) |>
        dplyr::summarize(do_union = FALSE) |>
        sf::st_cast("POLYGON") |>
        dplyr::ungroup() |>
        sf::st_buffer(0.0) |>
        sf::st_cast("MULTIPOLYGON") |>
        sf::st_cast("POLYGON", warn = FALSE) |>
        sf::st_buffer(0.0) |>
        sf::st_make_valid() |>
        dplyr::mutate(
          id = dplyr::row_number()
        ) |>
        sf::st_cast("POINT", warn = FALSE) |>
        dplyr::mutate(
          x = sf::st_coordinates(.data$geometry)[,1],
          y = sf::st_coordinates(.data$geometry)[,2]
        ) |>
        sf::st_drop_geometry() |>
        dplyr::group_by(.data$id) 
      if (dims$type == "horizontal") {
        line_points <- line_points |>
          dplyr::mutate(
            z = .data$y - min(.data$y),
            y = .data$y - .data$z
          )
      } else {
        line_points <- line_points |>
          dplyr::mutate(
            z = .data$x - min(.data$x),
            x = .data$x - .data$z
          )
      }
      line_points <- line_points |>
        sf::st_as_sf(
          coords = c("x", "y", "z"),
          dim = "XYZ",
          crs = sf::st_crs(lines)
        ) |>
        dplyr::group_by(.data$id) |>
        dplyr::filter(dplyr::n() >= 4) |>
        dplyr::summarize(do_union = FALSE) |>
        sf::st_cast("POLYGON") |>
        sf::st_cast("LINESTRING", warn = FALSE)
    } else {
      message("Returning polygons.")
      line_points <- line_points |>
        sf::st_buffer(0.0)
    }
    
  } else {
    line_points
  }
}