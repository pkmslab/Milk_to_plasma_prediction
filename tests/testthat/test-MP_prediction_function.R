aspirin <- digest::digest(MP_prediction_function("aspirin")) # digesting the result generated from using aspirin
custom <- digest::digest(MP_prediction_function("custom")) # digesting the result generated from using the custom feature

test_that("expect_equal test comparing the input of aspirin", {
  expect_equal(digest::digest(MP_prediction_function("aspirin")), aspirin)
})

test_that("expect_equal test comparing the input of custom", {
  expect_equal(digest::digest(MP_prediction_function("custom")), custom)
})
