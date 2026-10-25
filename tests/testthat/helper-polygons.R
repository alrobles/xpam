# helper-polygons.R -- shared polygon fixtures for test-pam_polygons.R
# and test-pam_poly_multi.R.
#
# Geometry convention: the env grid is nc x nr over extent [0, nc] x [0, nr]
# (1 unit per pixel, north-up, no CRS).  Test rectangles are INSET by 0.1
# inside the intended cell block, so the ALL_TOUCHED and pixel-centre rules
# agree exactly on the covered cell set -- which is then computable in
# plain R without depending on terra::rasterize's edge conventions.

# Cells (0-based raster order, row*nc + col) covered by the rectangle whose
# pixel-space bounds are (cx0, cy0_top)-(cx1, cy1_bot) in 0-based
# col/row-index terms: covers cols [cx0, cx1) and rows [ry0, ry1).
rect_cells <- function(cx0, cx1, ry0, ry1, nc) {
  as.vector(outer(ry0:(ry1 - 1), cx0:(cx1 - 1),
                  function(r, c) r * nc + c))
}
# WKT rectangle in raster coordinates: x in [cx0, cx1], y from the TOP of
# the grid (raster row ry <-> y in [nr - ry - 1, nr - ry]).
rect_wkt <- function(cx0, cx1, ry0, ry1, nr) {
  x0 <- cx0 + 0.1; x1 <- cx1 - 0.1
  y0 <- (nr - ry1) + 0.1; y1 <- (nr - ry0) - 0.1
  sprintf("POLYGON ((%g %g, %g %g, %g %g, %g %g, %g %g))",
          x0, y0, x1, y0, x1, y1, x0, y1, x0, y0)
}

read_xpoc <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  expect_equal(readChar(con, 4, useBytes = TRUE), "XPOC")
  expect_equal(readBin(con, "integer", 1, size = 4), 1L)
  np  <- readBin(con, "integer", 1, size = 8)   # header stores uint64
  nsp <- readBin(con, "integer", 1, size = 8)
  cells <- readBin(con, "double", np, size = 8)
  sps   <- readBin(con, "integer", np, size = 4)
  expect_equal(length(cells), np)
  expect_false(is.unsorted(cells))          # the OccSource contract
  list(cell = cells, sp = sps, n_species = nsp)
}

make_world <- function() {
  nr <- 20; nc <- 30
  r <- xpam_sim_env(nr, nc, 2L, seed = 1L)
  fenv <- tempfile(fileext = ".tif")
  terra::writeRaster(r, fenv, overwrite = TRUE)
  # species cell blocks (0-based, half-open [lo, hi))
  blocks <- list(
    sp1 = c(2, 6, 3, 8),      # cols 2:6, rows 3:8
    sp2 = c(10, 14, 12, 18),
    sp3 = c(4, 9, 6, 10),     # overlaps sp1 partially
    sp5 = c(20, 25, 1, 4))
  wkts <- vapply(blocks, function(b) rect_wkt(b[1], b[2], b[3], b[4], nr),
                 character(1))
  # sp4 is a multipolygon: two disjoint blocks
  ring <- function(cx0, cx1, ry0, ry1)
    sub("^POLYGON ", "", rect_wkt(cx0, cx1, ry0, ry1, nr))
  mp <- paste0("MULTIPOLYGON (", ring(15, 17, 2, 5), ",",
               ring(24, 27, 15, 18), ")")
  df <- data.frame(
    sciname = c(names(blocks), "sp4", "unmatched", NA),
    wkt = c(unname(wkts), mp,
            rect_wkt(0, 2, 0, 2, nr),          # unmatched species name
            rect_wkt(0, 2, 0, 2, nr)),          # NA field
    stringsAsFactors = FALSE)
  v <- terra::vect(df, geom = "wkt")
  fv <- tempfile(fileext = ".gpkg")
  terra::writeVector(v, fv, overwrite = TRUE)

  occ <- do.call(rbind, lapply(names(blocks), function(nm) {
    b <- blocks[[nm]]
    data.frame(cell = rect_cells(b[1], b[2], b[3], b[4], nc) + 1L,
               species = nm)
  }))
  occ <- rbind(occ, data.frame(
    cell = c(rect_cells(15, 17, 2, 5, nc), rect_cells(24, 27, 15, 18, nc)) + 1L,
    species = "sp4"))
  list(r = r, fenv = fenv, fv = fv, occ = occ, nr = nr, nc = nc)
}
