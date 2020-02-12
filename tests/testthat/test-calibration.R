test_that("multiplication works", {
  expect_equal_to_reference(rcarbon::calibrate(5000,50), "cal_5000_50.RDS")
  expect_equal_to_reference(rcarbon::calibrate(5000,50, ncores=2), "cal_5000_50.RDS")

})
