## Unit tests for functions from resid_fns.R
source('../../R/resid_fns.R')

context('osa function tests')
#use oneStepPredict example from TMB documentation

test_that('ar1xar1 cdf',{
  runExample("ar1xar1")
  osa.ar1xar1 <- oneStepPredict(obj, "N", "keep", method="cdf", discrete=TRUE, subset = 1:100)
  osa.myfun <- calculate.osa(obj, "cdf", "N", 'keep', Discrete = TRUE, Subset = 1:100)
  expect_equal(osa.ar1xar1$residual, osa.myfun$cdf)
})

test_that('simple, all',{
  runExample('simple')
  osa.simple.fg <- oneStepPredict(obj, observation.name = "x", method="fullGaussian")
  osa.simple.osg <- oneStepPredict(obj, observation.name = "x", 
                                   data.term.indicator = 'keep', method="oneStepGaussian")
  osa.simple.gen <- oneStepPredict(obj, observation.name = "x", 
                                   data.term.indicator = 'keep', method="oneStepGeneric")
  osa.simple.cdf <- oneStepPredict(obj, observation.name = "x", 
                                   data.term.indicator = 'keep', method="cdf")
  #gen not working?
  osa.myfun <- calculate.osa(obj,  methods=c("fg", "osg", "cdf"), "x")
  expect_equal(osa.simple.fg$residual, osa.myfun$fg)
  expect_equal(osa.simple.osg$residual, osa.myfun$osg)
  expect_equal(osa.simple.cdf$residual, osa.myfun$cdf)
  #test against each other
  # expect_equal(osa.simple.fg$residual, osa.simple.osg$residual)
  # expect_equal(osa.simple.fg$residual, osa.simple.gen$residual)
  # expect_equal(osa.simple.fg$residual, osa.simple.cdf$residual)
  # expect_equal(osa.simple.osg$residual, osa.simple.gen$residual)
  # expect_equal(osa.simple.osg$residual, osa.simple.cdf$residual)
  # expect_equal(osa.simple.gen$residual, osa.simple.cdf$residual)
  
})