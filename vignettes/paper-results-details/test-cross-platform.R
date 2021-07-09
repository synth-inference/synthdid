dir.1 = 'simulations-local'
dir.2 = 'simulations'
differences = do.call(rbind, lapply(intersect(list.files(dir.1), list.files(dir.2)), function(filename) {
    x=readRDS(sprintf('%s/%s', dir.1, filename))
    y=readRDS(sprintf('%s/%s', dir.2, filename))
    data.frame(estimator=x$estimator, estimate=x$estimate - y$estimate)
}))
max(abs(differences[differences$estimator != 'mc', 'estimate']))
quantile(abs(differences[differences$estimator == 'mc', 'estimate']), seq(0,10)/10)


