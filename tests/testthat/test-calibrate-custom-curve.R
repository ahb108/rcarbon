test_that("calibrate.UncalGrid accepts a custom curve given as a matrix", {
  curve <- rcarbon:::read_cal_curve_from_file("intcal20")
  expect_true(is.matrix(curve))

  grid <- as.UncalGrid(data.frame(CRA = 4000:6000, PrDens = 1))

  res_chr <- calibrate(grid, calCurves = "intcal20",
                       timeRange = c(8000, 3000), verbose = FALSE)
  res_mat <- calibrate(grid, calCurves = curve,
                       timeRange = c(8000, 3000), verbose = FALSE)
  res_df  <- calibrate(grid, calCurves = as.data.frame(curve),
                       timeRange = c(8000, 3000), verbose = FALSE)

  expect_s3_class(res_mat, "CalGrid")
  expect_equal(res_mat$calBP, res_chr$calBP)
  expect_equal(res_mat$PrDens, res_chr$PrDens)
  expect_equal(res_df$PrDens, res_chr$PrDens)
})

test_that("calibrate.UncalGrid(type = 'full') respects timeRange", {
  grid <- as.UncalGrid(data.frame(CRA = 4400:4600, PrDens = 1))
  res <- calibrate(grid, timeRange = c(8000, 3000), type = "full",
                   verbose = FALSE)
  expect_s3_class(res, "CalGrid")
  expect_true(all(res$calBP <= 8000 & res$calBP >= 3000))
})

test_that("calibrate.UncalGrid(type = 'full') uses the requested curve", {
  grid <- as.UncalGrid(data.frame(CRA = 4400:4600, PrDens = 1))
  res_int <- calibrate(grid, calCurves = "intcal20", type = "full",
                       timeRange = c(8000, 3000), verbose = FALSE)
  res_sh  <- calibrate(grid, calCurves = "shcal20", type = "full",
                       timeRange = c(8000, 3000), verbose = FALSE)
  expect_false(isTRUE(all.equal(res_int$PrDens, res_sh$PrDens)))
})

test_that("a custom curve covering timeRange leaves timeRange unchanged", {
  curve <- rcarbon:::read_cal_curve_from_file("intcal20")
  expect_warning(
    x <- calibrate(4500, 30, calCurves = curve, timeRange = c(8000, 3000),
                   verbose = FALSE),
    NA
  )
  expect_equal(x$metadata$StartBP, 8000)
  expect_equal(x$metadata$EndBP, 3000)
})

test_that("timeRange is only shrunk where the custom curve ends", {
  curve <- rcarbon:::read_cal_curve_from_file("intcal20")
  sub <- curve[curve[, 1] >= 2000 & curve[, 1] <= 10000, ]

  expect_warning(
    x <- calibrate(4500, 30, calCurves = sub, timeRange = c(12000, 0),
                   verbose = FALSE),
    "timeRange value not supported"
  )
  expect_equal(x$metadata$StartBP, max(sub[, 1]))
  expect_equal(x$metadata$EndBP, min(sub[, 1]))
})

