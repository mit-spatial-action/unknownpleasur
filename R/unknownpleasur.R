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

#' Section Lines to Polygons
#'
#' @param x `sf` object containing sectional lines.
#' @param baselines `sf` object containing baselines.
#' @param id_col Name of id column.
#' @param mode If `planar`, results will be `"POLYGON"`s. If `xyz`, results will be `"LINESTRINGS"``.
#'
#' @returns An `sf` object.
#' @export
up_polygonize <- function(x, baselines, id_col, mode = "planar") {
  baselines <- baselines |>
    sf::st_reverse() |>
    sf::st_cast("POINT", warn = FALSE)
  if (mode == "xyz") {
    baselines <- baselines |>
      dplyr::rowwise() |>
      dplyr::mutate(
        geometry = list(sf::st_point(c(sf::st_coordinates(geometry)[1:2], 0)))
      ) |>
      dplyr::ungroup() |>
      sf::st_set_crs(sf::st_crs(baselines))
  }
  x <- x |>
    dplyr::filter(
      .data[[id_col]] %in% dplyr::pull(baselines, dplyr::all_of(id_col)),
    ) |>
    sf::st_cast("POINT", warn = FALSE) |>
    dplyr::bind_rows(baselines) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::filter(
      dplyr::n() >= 4
    ) |>
    dplyr::summarize(do_union = FALSE)
  
  if (mode == "planar") {
    x <- x |>
      sf::st_cast("POLYGON")
  } else if (mode == "xyz") {
    x <- x |>
      sf::st_cast("LINESTRING")
  }
  
  x |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data[[id_col]]))
}

#' Transforms regularly-spaced lines into an "Unknown Pleasures"-esque set of regular section cuts based on raster value.
#'
#' @param x `sf` object that will be treated as an extent..
#' @param raster `raster` object.
#' @param interval Line and sample spacing. Can be `units` object. Otherwise will be in unit of CRS.
#' @param line_angle Angle (in degrees) at which to create lines.
#' @param elev_angle Angle at which to offset points.
#' @param scale One of "actual" or "interval." If "actual", use actual raster values. If "interval," base scale on interval between lines.
#' @param factor Z exaggeration, essentially.
#' @param mode If `planar`, results will be planar offset lines. If `xyz`, lines will be offset on `LINESTRING` z axis.
#' @param polygon If `TRUE`, outputs polygons (or closed linestrings if paired with `mode = "xyz"`. If `FALSE`, outputs lines with no baseline.
#' 
#' @returns `sf` object  containing "Unknown Pleasures" features.
#' 
#' @export
up_unknown_pleasures <- function(
    x, 
    raster,
    interval,
    line_angle,
    scale,
    factor,
    elev_angle = line_angle,
    mode = "planar",
    mask = TRUE,
    polygon = TRUE) {
  
  lines <- x |>
    up_regular_lines(interval, angle = line_angle) 
  
  if (mask) {
    lines <- lines |>
      up_mask_lines(x)
  }
  
  elev <- lines |>
    up_elevate(
      interval, 
      raster = dem, 
      angle = elev_angle, 
      scale = scale, 
      factor = factor, 
      mode = mode
      )
  
  if (polygon) {
    elev <- elev |>
      up_polygonize(elev, "id", mode="planar")
  }
  
  elev
}