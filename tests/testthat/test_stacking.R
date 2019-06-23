test_that('stacking function', {
  data(Env)
  data(Occurrences)
  Occ1 <- subset(Occurrences, Occurrences$SPECIES == 'elliptica')
  Occ2 <- subset(Occurrences, Occurrences$SPECIES == 'gracilis')
  ESDM1 <- ensemble_modelling(c('CTA', 'KSVM'), Occ1, Env, rep = 1,
                             Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                             name = 'elliptica', ensemble.thresh = 0,
                             verbose = FALSE)
  ESDM2 <- ensemble_modelling(c('CTA', 'KSVM'), Occ2, Env, rep = 1,
                             Xcol = 'LONGITUDE', Ycol = 'LATITUDE',
                             name = 'gracilis', ensemble.thresh = 0,
                             verbose = FALSE)
  SSDM <- stacking(ESDM1, ESDM2, verbose = FALSE)
  expect_is(SSDM, 'Stacked.SDM')
})
