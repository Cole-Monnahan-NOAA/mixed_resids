## Unit tests for functions from resid_fns.R
library(testthat)

test_that('simple, all',{
  runExample('simple')
  osa.simple.fg <- oneStepPredict(obj, observation.name = "x", method="fullGaussian")
  osa.simple.osg <- oneStepPredict(obj, observation.name = "x", 
                                   data.term.indicator = 'keep', method="oneStepGaussian")
  osa.simple.gen <- oneStepPredict(obj, observation.name = "x", 
                                   data.term.indicator = 'keep', method="oneStepGeneric")
  osa.simple.cdf <- oneStepPredict(obj, observation.name = "x", 
                                   data.term.indicator = 'keep', method="cdf")

  #test against each other
  expect_equal(osa.simple.fg$residual, osa.simple.osg$residual)
  expect_equal(osa.simple.fg$residual, osa.simple.gen$residual)
  expect_equal(osa.simple.fg$residual, osa.simple.cdf$residual)
  expect_equal(osa.simple.osg$residual, osa.simple.gen$residual)
  expect_equal(osa.simple.osg$residual, osa.simple.cdf$residual)
  expect_equal(osa.simple.gen$residual, osa.simple.cdf$residual)
  
})